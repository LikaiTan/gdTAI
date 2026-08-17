#!/usr/bin/env python3
"""Refine the V4.2 development pool to T/NK cells and reintegrate it.

This workflow is a non-destructive analysis sidecar. It uses the completed
V4.2 development integration only as a first-pass map, applies a conservative
cell-level T/NK gate to raw UMI counts, recomputes HVGs within the retained
pool, and then refits scVI and Leiden clustering. It does not define NK truth,
fit a gdTAI classifier, or mutate any source H5AD.
"""

from __future__ import annotations

import argparse
import gc
import html
import itertools
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Sequence

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import anndata as ad  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import seaborn as sns  # noqa: E402
from scipy import sparse  # noqa: E402
from sklearn.metrics import adjusted_rand_score  # noqa: E402

from workflows.gdtai.gdtai_v4_2_integration_core import (  # noqa: E402
    atomic_write_json,
    h5ad_obs_frame,
    h5ad_var_names,
    input_file_state,
    normalize_text,
    read_json,
    read_sparse_rows_genes,
    resolve,
    sha256_file,
    source_balanced_sample_indices,
)
from workflows.gdtai.run_gdtai_v4_2_nk_reference_integration import (  # noqa: E402
    atomic_write_h5ad,
    build_development_obs,
    development_roles,
    fit_stage,
    require_gpu,
)


DEFAULT_CONFIG = PROJECT_ROOT / "configs/models/gdtai/v4_2_tnk_reintegration.json"
SCRIPT_PATH = Path(__file__).resolve()
OTHER_HVG_EXCLUDED = re.compile(r"^(?:MT-|RPS|RPL|IG[HKL])", flags=re.IGNORECASE)
TCR_VJD = re.compile(r"^TR[ABGD](?:V|J|D)", flags=re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stage",
        choices=[
            "validate",
            "prepare",
            "fit",
            "cluster",
            "boundary",
            "report",
            "all",
        ],
        default="validate",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def load_runtime_config(path: Path) -> dict[str, Any]:
    config = read_json(path)
    base = read_json(config["base_execution_config"])
    for key in ("preflight", "inputs", "current_atlas_recovery"):
        config[key] = base[key]
    config["staging"]["batch_key_hierarchy"] = base["staging"]["batch_key_hierarchy"]
    config["_config_path"] = str(resolve(path))
    config["_config_sha256"] = sha256_file(path)
    config["_base_config_sha256"] = sha256_file(config["base_execution_config"])
    return config


def stage_paths(config: dict[str, Any]) -> dict[str, Path]:
    ssd = Path(config["resources"]["ssd_root"])
    outputs = config["outputs"]
    return {
        "ssd": ssd,
        "staged_h5ad": ssd / "tnk_refined_hvg_counts.h5ad",
        "marker_counts": ssd / "tnk_refined_marker_counts.npz",
        "model": ssd / "scvi_model",
        "latent": ssd / "X_scVI.npy",
        "partitions": ssd / "cluster_partitions.npz",
        "umap_sample": ssd / "diagnostic_umap_sample.npz",
        "boundary_partitions": ssd / "nk_boundary_partitions.npz",
        "boundary_umap_sample": ssd / "nk_boundary_umap_sample.npz",
        "boundary_tcr_flags": ssd / "nk_boundary_productive_tcr_flags.npz",
        "table": resolve(outputs["table_dir"]),
        "figure": resolve(outputs["figure_dir"]),
        "log": resolve(outputs["log_dir"]),
        "static": resolve(outputs["static_dir"]),
    }


def ensure_paths_and_resources(
    config: dict[str, Any], paths: dict[str, Path], stage: str
) -> dict[str, float]:
    paths["ssd"].mkdir(parents=True, exist_ok=True)
    for key in ("table", "figure", "log", "static"):
        paths[key].mkdir(parents=True, exist_ok=True)
    floor_key = {
        "prepare": "minimum_ssd_free_gib_prepare",
        "fit": "minimum_ssd_free_gib_fit",
        "cluster": "minimum_ssd_free_gib_cluster",
        "boundary": "minimum_ssd_free_gib_cluster",
        "report": "minimum_ssd_free_gib_cluster",
    }.get(stage, "minimum_ssd_free_gib_cluster")
    free = shutil.disk_usage(paths["ssd"]).free / 2**30
    required = float(config["resources"][floor_key])
    if free < required:
        raise RuntimeError(
            f"SSD free space {free:.1f} GiB is below the {required:.1f}-GiB floor for {stage}"
        )
    return {"observed_free_gib": free, "required_free_gib": required}


def marker_genes(config: dict[str, Any]) -> list[str]:
    gate = config["tnk_gate"]
    genes: list[str] = ["CD247", "LCK", "TRAT1", "BCL11B", "THEMIS"]
    for key in (
        "t_lineage_genes",
        "gd_lineage_genes",
        "nk_receptor_genes",
        "nk_adaptor_genes",
        "shared_cytotoxic_genes",
        "myeloid_context_genes",
    ):
        genes.extend(gate[key])
    return list(dict.fromkeys(genes))


def matrix_keys(config: dict[str, Any]) -> dict[str, str]:
    return read_json(config["preflight"]["config"])["matrix_keys"]


def load_first_pass_contract(
    config: dict[str, Any], obs: pd.DataFrame
) -> tuple[np.ndarray, dict[str, str]]:
    first = config["first_pass"]
    staged = Path(first["staged_h5ad"])
    latent = Path(first["latent"])
    partitions = Path(first["partitions"])
    for path in (staged, latent, partitions):
        if not path.is_file():
            raise FileNotFoundError(f"Missing first-pass artifact: {path}")
    backed = ad.read_h5ad(staged, backed="r")
    try:
        if backed.n_obs != obs.shape[0]:
            raise RuntimeError(
                "First-pass sidecar and reconstructed development obs differ in size"
            )
        if not backed.obs_names.equals(obs.index):
            raise RuntimeError(
                "First-pass sidecar and reconstructed development obs differ in order"
            )
    finally:
        backed.file.close()
    saved = np.load(partitions, allow_pickle=False)
    names = saved["global_names"].astype(str).tolist()
    run = str(first["cluster_run"])
    if run not in names:
        raise KeyError(f"First-pass cluster run is absent: {run}")
    labels = saved["global_labels"][:, names.index(run)].astype(np.int32)
    if labels.shape[0] != obs.shape[0]:
        raise RuntimeError("First-pass labels differ from development obs in size")
    provenance = {
        "first_pass_staged_h5ad_sha256": sha256_file(staged),
        "first_pass_latent_sha256": sha256_file(latent),
        "first_pass_partitions_sha256": sha256_file(partitions),
        "first_pass_cluster_run": run,
    }
    return labels, provenance


def extract_marker_counts(
    config: dict[str, Any], roles: pd.DataFrame, obs: pd.DataFrame, genes: Sequence[str]
) -> sparse.csr_matrix:
    keys = matrix_keys(config)
    blocks: list[sparse.csr_matrix] = []
    ordered: list[pd.DataFrame] = []
    for role in roles.itertuples(index=False):
        subset = obs.loc[obs["input_cohort_id"].eq(role.cohort_id)].sort_values(
            "input_row"
        )
        blocks.append(
            read_sparse_rows_genes(
                role.path,
                keys.get(role.cohort_id, keys["default"]),
                genes,
                rows=subset["input_row"].to_numpy(np.int64),
                row_chunk_size=int(config["staging"]["row_chunk_size"]),
            )
        )
        ordered.append(subset)
    ordered_obs = pd.concat(ordered, axis=0, copy=False)
    if not ordered_obs.index.equals(obs.index):
        raise RuntimeError("Marker extraction order differs from development obs")
    result = sparse.vstack(blocks, format="csr", dtype=np.float32)
    if result.shape != (obs.shape[0], len(genes)):
        raise RuntimeError("Marker extraction shape is inconsistent")
    return result


def hit_count(
    matrix: sparse.csr_matrix, genes: Sequence[str], selected: Sequence[str]
) -> np.ndarray:
    lookup = {gene: index for index, gene in enumerate(genes)}
    columns = [lookup[gene] for gene in selected]
    return np.asarray((matrix[:, columns] > 0).sum(axis=1)).ravel().astype(np.int16)


def detected(matrix: sparse.csr_matrix, genes: Sequence[str], gene: str) -> np.ndarray:
    index = list(genes).index(gene)
    return np.asarray(matrix[:, index].toarray()).ravel() > 0


def nonempty_cdr3(values: pd.Series) -> np.ndarray:
    return (
        ~values.astype("string")
        .fillna("")
        .str.strip()
        .str.lower()
        .isin(["", "nan", "none", "na"])
        .to_numpy()
    )


def restore_productive_tcr(obs: pd.DataFrame, roles: pd.DataFrame) -> np.ndarray:
    """Restore TRA/TRB evidence omitted from the original current-atlas sidecar."""
    productive = obs["productive_tcr_anchor"].astype(bool).to_numpy().copy()
    current_role = roles.loc[roles["cohort_id"].eq("current_atlas")]
    if current_role.shape[0] != 1:
        raise RuntimeError("Exactly one current-atlas role is required")
    current = h5ad_obs_frame(current_role.iloc[0]["path"], ["TRA_cdr3", "TRB_cdr3"])
    current_productive = nonempty_cdr3(current["TRA_cdr3"]) | nonempty_cdr3(
        current["TRB_cdr3"]
    )
    current_rows = np.flatnonzero(
        obs["input_cohort_id"].astype(str).eq("current_atlas").to_numpy()
    )
    input_rows = obs.iloc[current_rows]["input_row"].to_numpy(np.int64)
    productive[current_rows] |= current_productive[input_rows]
    primary = obs["primary_nk_anchor"].astype(bool).to_numpy()
    doublet = obs["doublet_flag_effective"].astype(bool).to_numpy()
    productive &= ~primary & ~doublet
    return productive


def apply_tnk_gate(
    config: dict[str, Any],
    obs: pd.DataFrame,
    genes: Sequence[str],
    counts: sparse.csr_matrix,
    first_labels: np.ndarray,
) -> tuple[np.ndarray, pd.DataFrame]:
    gate = config["tnk_gate"]
    t_hits = hit_count(counts, genes, gate["t_lineage_genes"])
    gd_hits = hit_count(counts, genes, gate["gd_lineage_genes"])
    nk_hits = hit_count(counts, genes, gate["nk_receptor_genes"])
    adaptor_hits = hit_count(counts, genes, gate["nk_adaptor_genes"])
    cytotoxic_hits = hit_count(counts, genes, gate["shared_cytotoxic_genes"])
    myeloid_hits = hit_count(counts, genes, gate["myeloid_context_genes"])
    cd3_hits = hit_count(counts, genes, ["CD3D", "CD3E", "CD3G"])
    delta_v = [
        gene
        for gene in gate["gd_lineage_genes"]
        if gene not in {"TRDC", "TRGC1", "TRGC2"}
    ]
    delta_v_hits = hit_count(counts, genes, delta_v)
    trdc = detected(counts, genes, "TRDC")

    primary_nk = obs["primary_nk_anchor"].astype(bool).to_numpy()
    if "productive_tcr_repaired" not in obs:
        raise KeyError(
            "productive_tcr_repaired is required before applying the T/NK gate"
        )
    productive = obs["productive_tcr_repaired"].astype(bool).to_numpy()
    protected = np.zeros(obs.shape[0], dtype=bool)
    if bool(gate["protect_primary_nk_anchors"]):
        protected |= primary_nk
    if bool(gate["protect_productive_tcr_cells"]):
        protected |= productive

    t_evidence = t_hits >= 1
    gd_evidence = gd_hits >= 1
    nk_receptor_evidence = nk_hits >= int(gate["minimum_nk_receptor_hits"])
    nk_cytotoxic_evidence = (
        cytotoxic_hits >= int(gate["minimum_shared_cytotoxic_hits_with_adaptor"])
    ) & (adaptor_hits >= int(gate["minimum_nk_adaptor_hits_with_cytotoxic"]))
    myeloid_conflict = (
        (myeloid_hits >= int(gate["minimum_myeloid_context_hits_for_conflict"]))
        & ~t_evidence
        & ~gd_evidence
        & ~nk_receptor_evidence
    )
    keep = (
        protected
        | t_evidence
        | gd_evidence
        | nk_receptor_evidence
        | nk_cytotoxic_evidence
    )
    keep &= ~(myeloid_conflict & ~protected)
    doublet = obs["doublet_flag_effective"].astype(bool).to_numpy()
    if bool(gate["remove_flagged_doublets"]):
        keep &= ~doublet

    reason = np.full(obs.shape[0], "removed_no_direct_tnk_evidence", dtype=object)
    reason[myeloid_conflict & ~protected] = "removed_myeloid_conflict"
    reason[keep & nk_cytotoxic_evidence] = "shared_cytotoxic_plus_nk_adaptor"
    reason[keep & nk_receptor_evidence] = "nk_receptor"
    reason[keep & gd_evidence] = "gd_lineage"
    reason[keep & t_evidence] = "t_lineage"
    reason[keep & productive] = "productive_tcr_protected"
    reason[keep & primary_nk] = "primary_nk_protected"
    reason[doublet & bool(gate["remove_flagged_doublets"])] = "removed_flagged_doublet"

    phenotype = np.full(obs.shape[0], "TRDC-/delta-V-", dtype=object)
    phenotype[~trdc & (delta_v_hits > 0)] = "TRDC-/delta-V+"
    phenotype[trdc & (delta_v_hits > 0)] = "TRDC+/delta-V+"
    phenotype[trdc & (delta_v_hits == 0) & (cd3_hits >= 2)] = (
        "TRDC+/delta-V-/CD3-strong"
    )
    phenotype[trdc & (delta_v_hits == 0) & (cd3_hits <= 1)] = "TRDC+/delta-V-/CD3-weak"

    audit = pd.DataFrame(
        {
            "tnk_gate_keep": keep,
            "tnk_gate_reason": reason,
            "first_pass_cluster": first_labels,
            "t_lineage_hits": t_hits,
            "gd_lineage_hits": gd_hits,
            "nk_receptor_hits": nk_hits,
            "nk_adaptor_hits": adaptor_hits,
            "shared_cytotoxic_hits": cytotoxic_hits,
            "myeloid_context_hits": myeloid_hits,
            "cd3_d_e_g_hits": cd3_hits,
            "delta_v_hits": delta_v_hits,
            "trdc_detected": trdc,
            "trdc_delta_v_cd3_phenotype": phenotype,
            "productive_tcr_protected": productive,
            "primary_nk_protected": primary_nk,
        },
        index=obs.index,
    )
    if not keep.any():
        raise RuntimeError("T/NK gate retained zero cells")
    if bool(gate["remove_flagged_doublets"]) and np.any(keep & doublet):
        raise RuntimeError("Flagged doublets passed a gate configured to remove them")
    if np.any(primary_nk & ~doublet & ~keep) or np.any(productive & ~doublet & ~keep):
        raise RuntimeError("A non-doublet protected anchor was removed")
    return keep, audit


def hvg_eligible_gene(gene: str) -> bool:
    return not OTHER_HVG_EXCLUDED.match(gene) and not TCR_VJD.match(gene)


def select_refined_hvgs(
    config: dict[str, Any],
    roles: pd.DataFrame,
    obs: pd.DataFrame,
) -> tuple[list[str], pd.DataFrame, pd.DataFrame]:
    common = set.intersection(
        *(
            set(h5ad_var_names(row.path).astype(str))
            for row in roles.itertuples(index=False)
        )
    )
    eligible = sorted(gene for gene in common if hvg_eligible_gene(gene))
    hvg = config["hvg"]
    n_top = int(hvg["n_top_genes"])
    if len(eligible) < n_top:
        raise RuntimeError(
            "Eligible common-gene universe is smaller than requested HVGs"
        )
    forced = list(dict.fromkeys(hvg["forced_lineage_genes"]))
    missing_forced = sorted(set(forced) - set(eligible))
    if missing_forced:
        raise RuntimeError(f"Forced lineage genes are unavailable: {missing_forced}")

    sampled = source_balanced_sample_indices(
        obs["source_gse_id"].astype(str).to_numpy(),
        int(hvg["sample_cap_per_source"]),
        int(config["random_seed"]),
    )
    sample_obs = obs.iloc[sampled].copy()
    keys = matrix_keys(config)
    blocks: list[sparse.csr_matrix] = []
    ordered: list[pd.DataFrame] = []
    for role in roles.itertuples(index=False):
        subset = sample_obs.loc[
            sample_obs["input_cohort_id"].eq(role.cohort_id)
        ].sort_values("input_row")
        blocks.append(
            read_sparse_rows_genes(
                role.path,
                keys.get(role.cohort_id, keys["default"]),
                eligible,
                rows=subset["input_row"].to_numpy(np.int64),
                row_chunk_size=int(config["staging"]["row_chunk_size"]),
            )
        )
        ordered.append(subset)
    sample_obs = pd.concat(ordered, axis=0, copy=False)
    matrix = sparse.vstack(blocks, format="csr", dtype=np.float32)
    if matrix.shape != (sample_obs.shape[0], len(eligible)):
        raise RuntimeError("HVG sample matrix is misaligned")
    sample_adata = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(
            {hvg["batch_key"]: sample_obs[hvg["batch_key"]].astype(str).to_numpy()},
            index=sample_obs.index,
        ),
        var=pd.DataFrame(index=pd.Index(eligible, name="gene")),
    )
    sc.pp.highly_variable_genes(
        sample_adata,
        flavor=str(hvg["flavor"]),
        n_top_genes=n_top,
        batch_key=str(hvg["batch_key"]),
        subset=False,
        inplace=True,
    )
    table = sample_adata.var.reset_index()
    rank = pd.to_numeric(table["highly_variable_rank"], errors="coerce")
    ranked = (
        table.assign(_rank=rank)
        .sort_values(["_rank", "gene"], na_position="last", kind="mergesort")["gene"]
        .astype(str)
        .tolist()
    )
    selected = ranked[:n_top]
    selected_set = set(selected)
    forced_set = set(forced)
    for gene in forced:
        if gene in selected_set:
            continue
        removable = next(
            candidate for candidate in reversed(selected) if candidate not in forced_set
        )
        selected.remove(removable)
        selected_set.remove(removable)
        selected.append(gene)
        selected_set.add(gene)
    selected = sorted(selected_set)
    if len(selected) != n_top or not forced_set.issubset(selected_set):
        raise RuntimeError(
            "Forced-HVG replacement did not produce the requested feature set"
        )
    table["forced_lineage_gene"] = table["gene"].isin(forced_set)
    table["selected_for_refined_scvi"] = table["gene"].isin(selected_set)
    source_sample = (
        sample_obs.groupby("source_gse_id", observed=True)
        .size()
        .rename("n_sampled")
        .reset_index()
    )
    del sample_adata, matrix, blocks
    gc.collect()
    return selected, table, source_sample


def atomic_sparse_npz(matrix: sparse.csr_matrix, path: Path) -> None:
    temporary = path.with_name(path.name + ".partial.npz")
    temporary.unlink(missing_ok=True)
    sparse.save_npz(temporary, matrix, compressed=True)
    temporary.replace(path)


def source_state_equal(before: pd.DataFrame, after: pd.DataFrame) -> pd.DataFrame:
    result = before[["cohort_id", "path", "size_bytes", "mtime_ns"]].merge(
        after[["cohort_id", "size_bytes", "mtime_ns"]],
        on="cohort_id",
        suffixes=("_before", "_after"),
    )
    result["unchanged"] = result["size_bytes_before"].eq(
        result["size_bytes_after"]
    ) & result["mtime_ns_before"].eq(result["mtime_ns_after"])
    return result


def prepare_stage(
    config: dict[str, Any], paths: dict[str, Path], overwrite: bool
) -> dict[str, Any]:
    for path in (paths["staged_h5ad"], paths["marker_counts"]):
        if path.exists() and not overwrite:
            raise FileExistsError(f"Refined staging artifact exists: {path}")
    roles = development_roles(config)
    before = input_file_state(roles)
    obs = build_development_obs(config, roles)
    obs["productive_tcr_repaired"] = restore_productive_tcr(obs, roles)
    first_labels, provenance = load_first_pass_contract(config, obs)
    genes = marker_genes(config)
    marker = extract_marker_counts(config, roles, obs, genes)
    keep, gate_audit = apply_tnk_gate(config, obs, genes, marker, first_labels)
    refined_obs = obs.loc[keep].copy()
    refined_audit = gate_audit.loc[keep].copy()
    refined_obs[refined_audit.columns] = refined_audit

    selected, hvg_table, source_sample = select_refined_hvgs(config, roles, refined_obs)
    keys = matrix_keys(config)
    blocks: list[sparse.csr_matrix] = []
    ordered: list[pd.DataFrame] = []
    for role in roles.itertuples(index=False):
        subset = refined_obs.loc[
            refined_obs["input_cohort_id"].eq(role.cohort_id)
        ].sort_values("input_row")
        blocks.append(
            read_sparse_rows_genes(
                role.path,
                keys.get(role.cohort_id, keys["default"]),
                selected,
                rows=subset["input_row"].to_numpy(np.int64),
                row_chunk_size=int(config["staging"]["row_chunk_size"]),
            )
        )
        ordered.append(subset)
    refined_obs = pd.concat(ordered, axis=0, copy=False)
    counts = sparse.vstack(blocks, format="csr", dtype=np.float32)
    if counts.shape != (refined_obs.shape[0], len(selected)):
        raise RuntimeError("Refined staged count matrix is misaligned")

    keep_columns = [
        "original_cell_id",
        "input_cohort_id",
        "input_row",
        "source_gse_id",
        "donor_id",
        "sample_id",
        "library_id",
        "technology_simple",
        "integration_batch",
        "integration_batch_level",
        "primary_nk_anchor",
        "productive_tcr_anchor",
        "productive_tcr_repaired",
        "doublet_flag_effective",
        "allow_pseudolabel_training",
        "candidate_eligible",
        *gate_audit.columns.tolist(),
    ]
    staged_obs = refined_obs[keep_columns].copy()
    categorical = [
        "input_cohort_id",
        "source_gse_id",
        "donor_id",
        "sample_id",
        "library_id",
        "technology_simple",
        "integration_batch",
        "integration_batch_level",
        "tnk_gate_reason",
        "trdc_delta_v_cd3_phenotype",
    ]
    for column in categorical:
        staged_obs[column] = pd.Categorical(
            staged_obs[column].astype("string").fillna("")
        )
    var = hvg_table.set_index("gene").loc[selected].copy()
    adata = ad.AnnData(X=counts, obs=staged_obs, var=var)
    adata.uns["gdtai_v4_2_tnk_reintegration"] = {
        "protocol_version": config["protocol_version"],
        "config_sha256": config["_config_sha256"],
        "script_sha256": sha256_file(SCRIPT_PATH),
        "raw_counts": True,
        "source_h5ad_mutation": False,
        **provenance,
    }
    ad.settings.allow_write_nullable_strings = True
    atomic_write_h5ad(adata, paths["staged_h5ad"], config["staging"]["compression"])
    marker_refined = marker[keep].tocsr()
    atomic_sparse_npz(marker_refined, paths["marker_counts"])

    gate_manifest = obs[
        ["original_cell_id", "input_cohort_id", "source_gse_id", "input_row"]
    ].join(gate_audit)
    removed = gate_manifest.loc[~keep].copy()
    removed.to_csv(paths["table"] / "tnk_gate_removed_cells.csv.gz", compression="gzip")
    gate_reason = (
        gate_manifest.groupby(["tnk_gate_keep", "tnk_gate_reason"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    gate_reason.to_csv(paths["table"] / "tnk_gate_reason_counts.csv", index=False)
    gate_source = (
        gate_manifest.groupby(["source_gse_id", "tnk_gate_keep"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    totals = gate_source.groupby("source_gse_id", observed=True)["n_cells"].transform(
        "sum"
    )
    gate_source["fraction_of_source"] = gate_source["n_cells"] / totals
    gate_source.to_csv(paths["table"] / "tnk_gate_counts_by_source.csv", index=False)
    phenotype = (
        gate_manifest.groupby(
            ["trdc_delta_v_cd3_phenotype", "tnk_gate_keep"], observed=True
        )
        .size()
        .rename("n_cells")
        .reset_index()
    )
    phenotype.to_csv(paths["table"] / "tnk_gate_phenotype_counts.csv", index=False)
    hvg_table.to_csv(paths["table"] / "refined_hvg_selection.csv", index=False)
    source_sample.to_csv(paths["table"] / "refined_hvg_source_sample.csv", index=False)
    pd.DataFrame({"gene": genes, "marker_index": np.arange(len(genes))}).to_csv(
        paths["table"] / "marker_panel.csv", index=False
    )
    before.to_csv(paths["table"] / "source_h5ad_state_before.csv", index=False)
    after = input_file_state(roles)
    unchanged = source_state_equal(before, after)
    unchanged.to_csv(paths["table"] / "source_h5ad_state_after.csv", index=False)
    if not unchanged["unchanged"].all():
        raise RuntimeError("A source H5AD changed during T/NK refinement staging")

    summary = {
        "stage": "prepare",
        "status": "PASS",
        "n_input_cells": int(obs.shape[0]),
        "n_retained_cells": int(adata.n_obs),
        "n_removed_cells": int(obs.shape[0] - adata.n_obs),
        "retained_fraction": float(adata.n_obs / obs.shape[0]),
        "n_hvgs": int(adata.n_vars),
        "n_forced_lineage_genes": int(
            adata.var["forced_lineage_gene"].astype(bool).sum()
        ),
        "matrix_nnz": int(adata.X.nnz),
        "marker_matrix_nnz": int(marker_refined.nnz),
        "staged_h5ad": str(paths["staged_h5ad"]),
        "staged_h5ad_sha256": sha256_file(paths["staged_h5ad"]),
        "marker_counts_sha256": sha256_file(paths["marker_counts"]),
        "primary_nk_anchors_retained": int(adata.obs["primary_nk_anchor"].sum()),
        "productive_tcr_cells_retained": int(
            adata.obs["productive_tcr_protected"].sum()
        ),
        "flagged_doublets_retained": int(adata.obs["doublet_flag_effective"].sum()),
        "source_h5ads_unchanged": True,
        **provenance,
    }
    atomic_write_json(paths["log"] / "prepare_summary.json", summary)
    return summary


def label_codes(values: pd.Series) -> np.ndarray:
    return pd.Categorical(values).codes.astype(np.int32)


def visualization_sample(obs: pd.DataFrame, config: dict[str, Any]) -> np.ndarray:
    clustering = config["clustering"]
    maximum = int(clustering["umap_sample_max_cells"])
    sources = obs["source_gse_id"].astype(str).to_numpy()
    n_sources = max(1, np.unique(sources).size)
    cap = min(
        int(clustering["umap_sample_cap_per_source"]),
        max(1, maximum // n_sources),
    )
    selected = source_balanced_sample_indices(sources, cap, int(config["random_seed"]))
    mandatory = np.flatnonzero(obs["primary_nk_anchor"].astype(bool).to_numpy())
    selected = np.union1d(selected, mandatory)
    if selected.size > maximum:
        optional = np.setdiff1d(selected, mandatory, assume_unique=True)
        n_optional = max(0, maximum - mandatory.size)
        rng = np.random.default_rng(int(config["random_seed"]) + 991)
        optional = np.sort(rng.choice(optional, size=n_optional, replace=False))
        selected = np.union1d(mandatory, optional)
    return selected.astype(np.int64)


def cluster_stage(
    config: dict[str, Any], paths: dict[str, Path], overwrite: bool
) -> dict[str, Any]:
    for path in (paths["partitions"], paths["umap_sample"]):
        if path.exists() and not overwrite:
            raise FileExistsError(f"Refined clustering artifact exists: {path}")
    gpu = require_gpu(config)
    import rapids_singlecell as rsc

    backed = ad.read_h5ad(paths["staged_h5ad"], backed="r")
    obs = backed.obs.copy()
    backed.file.close()
    latent = np.load(paths["latent"], mmap_mode="r")
    if latent.shape[0] != obs.shape[0]:
        raise RuntimeError("Refined latent and staged obs differ in size")
    graph = ad.AnnData(
        X=sparse.csr_matrix((obs.shape[0], 1), dtype=np.float32), obs=obs
    )
    graph.obsm["X_scVI"] = np.asarray(latent)
    cl = config["clustering"]
    rsc.pp.neighbors(
        graph,
        n_neighbors=int(cl["n_neighbors"]),
        use_rep="X_scVI",
        random_state=int(config["random_seed"]),
        algorithm=cl["algorithm"],
        algorithm_kwds=cl["algorithm_kwds"],
    )
    labels: list[np.ndarray] = []
    names: list[str] = []
    for resolution in cl["resolutions"]:
        for seed in cl["seeds"]:
            key = f"leiden_r{float(resolution):g}_s{int(seed)}"
            rsc.tl.leiden(
                graph,
                resolution=float(resolution),
                random_state=int(seed),
                key_added=key,
            )
            labels.append(label_codes(graph.obs[key]))
            names.append(key)
    label_matrix = np.column_stack(labels).astype(np.int32)
    temporary = paths["partitions"].with_name(paths["partitions"].name + ".partial.npz")
    np.savez_compressed(
        temporary,
        labels=label_matrix,
        names=np.asarray(names),
        staged_h5ad_sha256=np.asarray(sha256_file(paths["staged_h5ad"])),
        latent_sha256=np.asarray(sha256_file(paths["latent"])),
    )
    temporary.replace(paths["partitions"])

    sample_indices = visualization_sample(obs, config)
    sample = ad.AnnData(
        X=sparse.csr_matrix((sample_indices.size, 1), dtype=np.float32),
        obs=obs.iloc[sample_indices].copy(),
    )
    sample.obsm["X_scVI"] = np.asarray(latent[sample_indices])
    rsc.pp.neighbors(
        sample,
        n_neighbors=int(cl["n_neighbors"]),
        use_rep="X_scVI",
        random_state=int(config["random_seed"]),
        algorithm=cl["algorithm"],
        algorithm_kwds=cl["algorithm_kwds"],
    )
    rsc.tl.umap(
        sample,
        min_dist=float(cl["umap_min_dist"]),
        spread=float(cl["umap_spread"]),
        random_state=int(config["random_seed"]),
    )
    coordinates = np.asarray(sample.obsm["X_umap"], dtype=np.float32)
    if (
        coordinates.shape != (sample_indices.size, 2)
        or not np.isfinite(coordinates).all()
    ):
        raise RuntimeError("Diagnostic UMAP is invalid")
    review = str(cl["review_run"])
    if review not in names:
        raise KeyError(f"Review run is absent: {review}")
    review_labels = label_matrix[:, names.index(review)]
    np.savez_compressed(
        paths["umap_sample"],
        sample_indices=sample_indices,
        X_umap=coordinates,
        refined_cluster=review_labels[sample_indices],
        source_gse_id=np.asarray(
            obs.iloc[sample_indices]["source_gse_id"].astype(str).to_numpy(),
            dtype=str,
        ),
        first_pass_cluster=obs.iloc[sample_indices]["first_pass_cluster"].to_numpy(
            np.int32
        ),
        phenotype=np.asarray(
            obs.iloc[sample_indices]["trdc_delta_v_cd3_phenotype"]
            .astype(str)
            .to_numpy(),
            dtype=str,
        ),
        primary_nk_anchor=obs.iloc[sample_indices]["primary_nk_anchor"]
        .astype(bool)
        .to_numpy(),
        productive_tcr=obs.iloc[sample_indices]["productive_tcr_protected"]
        .astype(bool)
        .to_numpy(),
    )

    cluster_counts = (
        pd.DataFrame(
            {
                "cluster": review_labels,
                "source_gse_id": obs["source_gse_id"].astype(str).to_numpy(),
                "phenotype": obs["trdc_delta_v_cd3_phenotype"].astype(str).to_numpy(),
                "primary_nk_anchor": obs["primary_nk_anchor"].astype(bool).to_numpy(),
                "productive_tcr": obs["productive_tcr_protected"]
                .astype(bool)
                .to_numpy(),
            }
        )
        .groupby("cluster", observed=True)
        .agg(
            n_cells=("cluster", "size"),
            n_primary_nk=("primary_nk_anchor", "sum"),
            n_productive_tcr=("productive_tcr", "sum"),
        )
        .reset_index()
    )
    cluster_counts.to_csv(paths["table"] / "refined_cluster_counts.csv", index=False)
    pd.crosstab(obs["first_pass_cluster"], review_labels).to_csv(
        paths["table"] / "first_pass_to_refined_cluster_counts.csv"
    )
    summary = {
        "stage": "cluster",
        "status": "PASS",
        "n_cells": int(obs.shape[0]),
        "n_runs": len(names),
        "review_run": review,
        "n_review_clusters": int(np.unique(review_labels).size),
        "n_umap_sample": int(sample_indices.size),
        "partitions_sha256": sha256_file(paths["partitions"]),
        "umap_sample_sha256": sha256_file(paths["umap_sample"]),
        "gpu": gpu,
        "cpu_fallback": False,
    }
    atomic_write_json(paths["log"] / "cluster_summary.json", summary)
    return summary


def boundary_visualization_sample(
    obs: pd.DataFrame, config: dict[str, Any]
) -> np.ndarray:
    boundary = config["nk_boundary"]
    maximum = int(boundary["umap_sample_max_cells"])
    sources = obs["source_gse_id"].astype(str).to_numpy()
    selected = source_balanced_sample_indices(
        sources,
        int(boundary["umap_sample_cap_per_source"]),
        int(config["random_seed"]) + 1701,
    )
    mandatory = np.flatnonzero(obs["primary_nk_anchor"].astype(bool).to_numpy())
    selected = np.union1d(selected, mandatory)
    if selected.size > maximum:
        optional = np.setdiff1d(selected, mandatory, assume_unique=True)
        n_optional = max(0, maximum - mandatory.size)
        rng = np.random.default_rng(int(config["random_seed"]) + 1702)
        optional = np.sort(rng.choice(optional, size=n_optional, replace=False))
        selected = np.union1d(mandatory, optional)
    return selected.astype(np.int64)


def boundary_stage(
    config: dict[str, Any], paths: dict[str, Path], overwrite: bool
) -> dict[str, Any]:
    for path in (paths["boundary_partitions"], paths["boundary_umap_sample"]):
        if path.exists() and not overwrite:
            raise FileExistsError(f"NK-boundary artifact exists: {path}")
    gpu = require_gpu(config)
    import rapids_singlecell as rsc

    backed = ad.read_h5ad(paths["staged_h5ad"], backed="r")
    obs = backed.obs.copy()
    backed.file.close()
    latent = np.load(paths["latent"], mmap_mode="r")
    saved = np.load(paths["partitions"], allow_pickle=False)
    global_names = saved["names"].astype(str).tolist()
    boundary = config["nk_boundary"]
    parent_run = str(boundary["parent_review_run"])
    if parent_run not in global_names:
        raise KeyError(f"Parent review run is absent: {parent_run}")
    parent_labels = saved["labels"][:, global_names.index(parent_run)].astype(np.int32)
    requested_parents = np.asarray(boundary["parent_clusters"], dtype=np.int32)
    boundary_indices = np.flatnonzero(np.isin(parent_labels, requested_parents))
    if boundary_indices.size == 0:
        raise RuntimeError("The configured NK boundary contains zero cells")
    boundary_obs = obs.iloc[boundary_indices].copy()
    boundary_obs["parent_refined_cluster"] = parent_labels[boundary_indices]

    graph = ad.AnnData(
        X=sparse.csr_matrix((boundary_indices.size, 1), dtype=np.float32),
        obs=boundary_obs,
    )
    graph.obsm["X_scVI"] = np.asarray(latent[boundary_indices])
    clustering = config["clustering"]
    rsc.pp.neighbors(
        graph,
        n_neighbors=int(boundary["n_neighbors"]),
        use_rep="X_scVI",
        random_state=int(config["random_seed"]),
        algorithm=clustering["algorithm"],
        algorithm_kwds=clustering["algorithm_kwds"],
    )
    names: list[str] = []
    labels: list[np.ndarray] = []
    for resolution in boundary["resolutions"]:
        for seed in boundary["seeds"]:
            key = f"boundary_r{float(resolution):g}_s{int(seed)}"
            rsc.tl.leiden(
                graph,
                resolution=float(resolution),
                random_state=int(seed),
                key_added=key,
            )
            names.append(key)
            labels.append(label_codes(graph.obs[key]))
    label_matrix = np.column_stack(labels).astype(np.int32)
    temporary = paths["boundary_partitions"].with_name(
        paths["boundary_partitions"].name + ".partial.npz"
    )
    np.savez_compressed(
        temporary,
        boundary_indices=boundary_indices.astype(np.int64),
        parent_labels=parent_labels[boundary_indices],
        labels=label_matrix,
        names=np.asarray(names),
        staged_h5ad_sha256=np.asarray(sha256_file(paths["staged_h5ad"])),
        latent_sha256=np.asarray(sha256_file(paths["latent"])),
        global_partitions_sha256=np.asarray(sha256_file(paths["partitions"])),
    )
    temporary.replace(paths["boundary_partitions"])

    review = str(boundary["review_run"])
    if review not in names:
        raise KeyError(f"Boundary review run is absent: {review}")
    review_labels = label_matrix[:, names.index(review)]
    sample_local = boundary_visualization_sample(boundary_obs, config)
    sample = ad.AnnData(
        X=sparse.csr_matrix((sample_local.size, 1), dtype=np.float32),
        obs=boundary_obs.iloc[sample_local].copy(),
    )
    sample.obsm["X_scVI"] = np.asarray(latent[boundary_indices[sample_local]])
    rsc.pp.neighbors(
        sample,
        n_neighbors=int(boundary["n_neighbors"]),
        use_rep="X_scVI",
        random_state=int(config["random_seed"]) + 1703,
        algorithm=clustering["algorithm"],
        algorithm_kwds=clustering["algorithm_kwds"],
    )
    rsc.tl.umap(
        sample,
        min_dist=float(boundary["umap_min_dist"]),
        spread=float(boundary["umap_spread"]),
        random_state=int(config["random_seed"]) + 1704,
    )
    coordinates = np.asarray(sample.obsm["X_umap"], dtype=np.float32)
    if (
        coordinates.shape != (sample_local.size, 2)
        or not np.isfinite(coordinates).all()
    ):
        raise RuntimeError("NK-boundary diagnostic UMAP is invalid")
    np.savez_compressed(
        paths["boundary_umap_sample"],
        sample_local_indices=sample_local,
        sample_global_indices=boundary_indices[sample_local],
        X_umap=coordinates,
        boundary_cluster=review_labels[sample_local],
        parent_refined_cluster=parent_labels[boundary_indices[sample_local]],
        source_gse_id=np.asarray(
            boundary_obs.iloc[sample_local]["source_gse_id"].astype(str).to_numpy(),
            dtype=str,
        ),
        phenotype=np.asarray(
            boundary_obs.iloc[sample_local]["trdc_delta_v_cd3_phenotype"]
            .astype(str)
            .to_numpy(),
            dtype=str,
        ),
        primary_nk_anchor=boundary_obs.iloc[sample_local]["primary_nk_anchor"]
        .astype(bool)
        .to_numpy(),
        productive_tcr=boundary_obs.iloc[sample_local]["productive_tcr_protected"]
        .astype(bool)
        .to_numpy(),
    )

    audit = pd.DataFrame(
        {
            "boundary_cluster": review_labels,
            "parent_refined_cluster": parent_labels[boundary_indices],
            "source_gse_id": boundary_obs["source_gse_id"].astype(str).to_numpy(),
            "phenotype": boundary_obs["trdc_delta_v_cd3_phenotype"]
            .astype(str)
            .to_numpy(),
            "primary_nk_anchor": boundary_obs["primary_nk_anchor"]
            .astype(bool)
            .to_numpy(),
            "productive_tcr": boundary_obs["productive_tcr_protected"]
            .astype(bool)
            .to_numpy(),
        }
    )
    counts = (
        audit.groupby("boundary_cluster", observed=True)
        .agg(
            n_cells=("boundary_cluster", "size"),
            n_primary_nk=("primary_nk_anchor", "sum"),
            n_productive_tcr=("productive_tcr", "sum"),
            n_sources=("source_gse_id", "nunique"),
        )
        .reset_index()
    )
    source_counts = (
        audit.groupby(["boundary_cluster", "source_gse_id"], observed=True)
        .agg(
            n_cells=("boundary_cluster", "size"),
            n_primary_nk=("primary_nk_anchor", "sum"),
            n_productive_tcr=("productive_tcr", "sum"),
        )
        .reset_index()
    )
    dominant = (
        source_counts.sort_values(
            ["boundary_cluster", "n_cells", "source_gse_id"],
            ascending=[True, False, True],
        )
        .groupby("boundary_cluster", observed=True)
        .first()
        .reset_index()[["boundary_cluster", "source_gse_id", "n_cells"]]
        .rename(
            columns={
                "source_gse_id": "dominant_source",
                "n_cells": "n_dominant_source",
            }
        )
    )
    counts = counts.merge(dominant, on="boundary_cluster", how="left", validate="1:1")
    counts["primary_nk_fraction"] = counts["n_primary_nk"] / counts["n_cells"]
    counts["productive_tcr_fraction"] = counts["n_productive_tcr"] / counts["n_cells"]
    counts["dominant_source_fraction"] = counts["n_dominant_source"] / counts["n_cells"]
    counts.to_csv(paths["table"] / "nk_boundary_cluster_counts.csv", index=False)
    source_counts.to_csv(
        paths["table"] / "nk_boundary_cluster_counts_by_source.csv", index=False
    )
    pd.crosstab(audit["parent_refined_cluster"], audit["boundary_cluster"]).to_csv(
        paths["table"] / "nk_boundary_parent_transition_counts.csv"
    )
    (
        audit.groupby(["boundary_cluster", "phenotype"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .to_csv(paths["table"] / "nk_boundary_phenotype_counts.csv", index=False)
    )
    summary = {
        "stage": "boundary",
        "status": "PASS_REVIEW_REQUIRED",
        "parent_review_run": parent_run,
        "parent_clusters": requested_parents.astype(int).tolist(),
        "n_boundary_cells": int(boundary_indices.size),
        "n_primary_nk": int(audit["primary_nk_anchor"].sum()),
        "primary_nk_recall": float(
            audit["primary_nk_anchor"].sum()
            / obs["primary_nk_anchor"].astype(bool).sum()
        ),
        "n_productive_tcr": int(audit["productive_tcr"].sum()),
        "n_runs": len(names),
        "review_run": review,
        "n_review_clusters": int(np.unique(review_labels).size),
        "n_umap_sample": int(sample_local.size),
        "boundary_partitions_sha256": sha256_file(paths["boundary_partitions"]),
        "boundary_umap_sample_sha256": sha256_file(paths["boundary_umap_sample"]),
        "gpu": gpu,
        "cpu_fallback": False,
    }
    atomic_write_json(paths["log"] / "boundary_summary.json", summary)
    return summary


def plotting_defaults() -> None:
    sns.set_theme(style="white", context="paper")
    plt.rcParams.update(
        {
            "figure.dpi": 180,
            "savefig.dpi": 300,
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "axes.linewidth": 0.7,
        }
    )


def clean_umap(ax: plt.Axes, title: str) -> None:
    ax.set_title(title)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(length=2)


def save_figure(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_clusters(sample: dict[str, np.ndarray], path: Path) -> None:
    xy = sample["X_umap"]
    labels = sample["refined_cluster"].astype(int)
    fig, ax = plt.subplots(figsize=(8.2, 6.7))
    ax.scatter(
        xy[:, 0],
        xy[:, 1],
        c=labels,
        cmap="tab20",
        s=0.35,
        linewidths=0,
        rasterized=True,
    )
    for cluster in np.unique(labels):
        center = np.median(xy[labels == cluster], axis=0)
        ax.text(
            center[0],
            center[1],
            str(cluster),
            ha="center",
            va="center",
            fontsize=7,
            bbox={
                "boxstyle": "round,pad=0.14",
                "fc": "white",
                "ec": "none",
                "alpha": 0.8,
            },
        )
    clean_umap(ax, "T/NK-restricted scVI and Leiden clusters")
    save_figure(fig, path)


def plot_categories(
    sample: dict[str, np.ndarray], values: np.ndarray, title: str, path: Path
) -> None:
    xy = sample["X_umap"]
    categories = pd.Series(values.astype(str)).fillna("").unique().tolist()
    phenotype_order = [
        "TRDC-/delta-V-",
        "TRDC-/delta-V+",
        "TRDC+/delta-V+",
        "TRDC+/delta-V-/CD3-strong",
        "TRDC+/delta-V-/CD3-weak",
    ]
    if set(categories).issubset(phenotype_order):
        ordered = [category for category in phenotype_order if category in categories]
        color_lookup = {
            "TRDC-/delta-V-": "#c9ced3",
            "TRDC-/delta-V+": "#5f9e6e",
            "TRDC+/delta-V+": "#2a6fbb",
            "TRDC+/delta-V-/CD3-strong": "#8d5aa6",
            "TRDC+/delta-V-/CD3-weak": "#d4552d",
        }
        colors = [color_lookup[category] for category in ordered]
    else:
        ordered = sorted(categories)
        colors = sns.color_palette("tab20", n_colors=max(3, len(categories)))
    fig, ax = plt.subplots(figsize=(8.4, 6.7))
    ax.scatter(
        xy[:, 0], xy[:, 1], s=0.25, color="#d5d9dd", linewidths=0, rasterized=True
    )
    for color, category in zip(colors, ordered, strict=False):
        mask = values.astype(str) == category
        ax.scatter(
            xy[mask, 0],
            xy[mask, 1],
            s=0.45,
            color=color,
            linewidths=0,
            rasterized=True,
            label=f"{category} ({mask.sum():,})",
        )
    clean_umap(ax, title)
    ax.legend(
        loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, markerscale=5
    )
    save_figure(fig, path)


def plot_features(
    sample: dict[str, np.ndarray],
    expression: np.ndarray,
    genes: Sequence[str],
    selected: Sequence[str],
    title: str,
    path: Path,
) -> None:
    lookup = {gene: index for index, gene in enumerate(genes)}
    chosen = [gene for gene in selected if gene in lookup]
    ncols = 4
    nrows = math.ceil(len(chosen) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(12.0, 2.65 * nrows), squeeze=False)
    xy = sample["X_umap"]
    for ax, gene in zip(axes.ravel(), chosen, strict=False):
        values = expression[:, lookup[gene]]
        order = np.argsort(values, kind="mergesort")
        ax.scatter(
            xy[order, 0],
            xy[order, 1],
            c=values[order],
            cmap="viridis",
            s=0.32,
            linewidths=0,
            rasterized=True,
        )
        clean_umap(ax, f"{gene}  n+={(values > 0).sum():,}")
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in axes.ravel()[len(chosen) :]:
        ax.axis("off")
    fig.suptitle(title, fontsize=12, y=1.002)
    save_figure(fig, path)


def cluster_marker_table(
    labels: np.ndarray,
    counts: sparse.csr_matrix,
    genes: Sequence[str],
) -> pd.DataFrame:
    clusters = np.unique(labels)
    code = {cluster: index for index, cluster in enumerate(clusters)}
    rows = np.asarray([code[value] for value in labels], dtype=np.int32)
    membership = sparse.csr_matrix(
        (np.ones(labels.size, dtype=np.float32), (rows, np.arange(labels.size))),
        shape=(clusters.size, labels.size),
    )
    detected_counts = (membership @ (counts > 0)).toarray().astype(np.float64)
    sizes = np.asarray(membership.sum(axis=1)).ravel()
    fractions = detected_counts / sizes[:, None]
    records = []
    for row, cluster in enumerate(clusters):
        for column, gene in enumerate(genes):
            records.append(
                {
                    "cluster": int(cluster),
                    "n_cluster": int(sizes[row]),
                    "gene": gene,
                    "n_detected": int(detected_counts[row, column]),
                    "detected_fraction": float(fractions[row, column]),
                }
            )
    return pd.DataFrame(records)


def plot_cluster_heatmap(
    table: pd.DataFrame,
    genes: Sequence[str],
    path: Path,
    title: str = "Raw-UMI marker detection by refined cluster",
    ylabel: str = "Refined cluster",
) -> None:
    matrix = table.pivot(index="cluster", columns="gene", values="detected_fraction")
    matrix = matrix.reindex(columns=[gene for gene in genes if gene in matrix.columns])
    fig, ax = plt.subplots(figsize=(13.2, max(5.0, 0.27 * matrix.shape[0] + 1.5)))
    sns.heatmap(
        matrix,
        cmap="mako",
        vmin=0,
        vmax=1,
        linewidths=0.15,
        linecolor="white",
        cbar_kws={"label": "Fraction detected"},
        ax=ax,
    )
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    save_figure(fig, path)


def plot_boundary_clusters(sample: dict[str, np.ndarray], path: Path) -> None:
    xy = sample["X_umap"]
    labels = sample["boundary_cluster"].astype(int)
    fig, ax = plt.subplots(figsize=(8.2, 6.7))
    colors = sns.color_palette("colorblind", n_colors=np.unique(labels).size)
    for color, cluster in zip(colors, np.unique(labels), strict=False):
        mask = labels == cluster
        ax.scatter(
            xy[mask, 0],
            xy[mask, 1],
            s=0.38,
            color=color,
            linewidths=0,
            rasterized=True,
            label=f"{cluster} ({mask.sum():,})",
        )
        center = np.median(xy[mask], axis=0)
        ax.text(
            center[0],
            center[1],
            str(cluster),
            ha="center",
            va="center",
            fontsize=8,
            bbox={"boxstyle": "round,pad=0.14", "fc": "white", "ec": "none"},
        )
    clean_umap(ax, "Second-pass clustering inside the NK-like boundary")
    ax.legend(
        loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, markerscale=5
    )
    save_figure(fig, path)


def plot_boundary_evidence(cluster_counts: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 6.2))
    sizes = 70 + 900 * np.sqrt(
        cluster_counts["n_cells"] / cluster_counts["n_cells"].max()
    )
    plotted = ax.scatter(
        100 * cluster_counts["fraction_productive_ab"],
        100 * cluster_counts["primary_nk_fraction"],
        s=sizes,
        c=100 * cluster_counts["dominant_source_fraction"],
        cmap="viridis_r",
        vmin=0,
        vmax=100,
        edgecolor="white",
        linewidth=0.7,
    )
    for row in cluster_counts.itertuples(index=False):
        ax.annotate(
            str(row.boundary_cluster),
            (100 * row.fraction_productive_ab, 100 * row.primary_nk_fraction),
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold",
        )
    ax.set_xlabel("Harmonized productive TRA or TRB in cluster (%)")
    ax.set_ylabel("Primary NK annotation anchors in cluster (%)")
    ax.set_title("Boundary annotation conflict; bubble area reflects cluster size")
    ax.spines[["top", "right"]].set_visible(False)
    fig.colorbar(plotted, ax=ax, label="Dominant dataset contribution (%)")
    save_figure(fig, path)


def plot_boundary_source_heatmap(source_counts: pd.DataFrame, path: Path) -> None:
    top_sources = (
        source_counts.groupby("source_gse_id", observed=True)["n_cells"]
        .sum()
        .nlargest(15)
        .index
    )
    selected = source_counts.loc[source_counts["source_gse_id"].isin(top_sources)]
    matrix = selected.pivot_table(
        index="boundary_cluster",
        columns="source_gse_id",
        values="n_cells",
        aggfunc="sum",
        fill_value=0,
    ).reindex(columns=top_sources)
    matrix = matrix.div(matrix.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    fig, ax = plt.subplots(figsize=(12.5, 5.4))
    sns.heatmap(
        matrix,
        cmap="crest",
        vmin=0,
        vmax=max(0.25, float(matrix.to_numpy().max())),
        cbar_kws={"label": "Fraction of boundary cluster"},
        ax=ax,
    )
    ax.set_title("Top dataset contributions to each boundary cluster")
    ax.set_xlabel("")
    ax.set_ylabel("Boundary cluster")
    ax.tick_params(axis="x", rotation=45)
    save_figure(fig, path)


def boundary_stability_table(
    names: Sequence[str], labels: np.ndarray
) -> tuple[pd.DataFrame, pd.DataFrame]:
    records: list[dict[str, Any]] = []
    for resolution in sorted(
        {float(re.search(r"boundary_r([0-9.]+)_", name).group(1)) for name in names}
    ):
        indices = [
            index
            for index, name in enumerate(names)
            if name.startswith(f"boundary_r{resolution:g}_")
        ]
        for left, right in itertools.combinations(indices, 2):
            records.append(
                {
                    "resolution": resolution,
                    "run_a": names[left],
                    "run_b": names[right],
                    "n_clusters_a": int(np.unique(labels[:, left]).size),
                    "n_clusters_b": int(np.unique(labels[:, right]).size),
                    "adjusted_rand_index": float(
                        adjusted_rand_score(labels[:, left], labels[:, right])
                    ),
                }
            )
    pairs = pd.DataFrame(records)
    summary = (
        pairs.groupby("resolution", observed=True)
        .agg(
            n_seed_pairs=("adjusted_rand_index", "size"),
            mean_adjusted_rand_index=("adjusted_rand_index", "mean"),
            min_adjusted_rand_index=("adjusted_rand_index", "min"),
            max_adjusted_rand_index=("adjusted_rand_index", "max"),
        )
        .reset_index()
    )
    return pairs, summary


def boundary_productive_chain_flags(
    config: dict[str, Any],
    paths: dict[str, Path],
    obs: pd.DataFrame,
    boundary_indices: np.ndarray,
) -> dict[str, np.ndarray]:
    cache = paths["boundary_tcr_flags"]
    boundary_hash = sha256_file(paths["boundary_partitions"])
    if cache.is_file():
        saved = np.load(cache, allow_pickle=False)
        if str(saved["boundary_partitions_sha256"]) == boundary_hash and np.array_equal(
            saved["boundary_indices"], boundary_indices
        ):
            return {
                chain: saved[f"productive_{chain.lower()}"].astype(bool)
                for chain in ("TRA", "TRB", "TRD")
            }

    chain_columns = ["TRA_cdr3", "TRB_cdr3", "TRD_cdr3"]
    boundary_obs = obs.iloc[boundary_indices].copy()
    flags = {
        chain: np.zeros(boundary_indices.size, dtype=bool)
        for chain in ("TRA", "TRB", "TRD")
    }
    covered = np.zeros(boundary_indices.size, dtype=bool)
    roles = development_roles(config)
    for role in roles.itertuples(index=False):
        local = np.flatnonzero(
            boundary_obs["input_cohort_id"].astype(str).eq(role.cohort_id).to_numpy()
        )
        if local.size == 0:
            continue
        if role.cohort_id == "current_atlas":
            target = boundary_obs.iloc[local][["source_gse_id"]].copy()
            source_rows = boundary_obs.iloc[local]["input_row"].to_numpy(np.int64)
            source_ids = h5ad_obs_frame(role.path, ["original_cell_id"])
            target["original_cell_id"] = source_ids.iloc[source_rows][
                "original_cell_id"
            ].to_numpy()
            target["source_gse_id"] = normalize_text(target["source_gse_id"])
            target["original_cell_id"] = normalize_text(target["original_cell_id"])
            target_keys = (
                target["source_gse_id"] + "||" + target["original_cell_id"]
            ).to_numpy()
            if pd.Index(target_keys).duplicated().any():
                raise RuntimeError("Current-atlas boundary TCR keys are not unique")
            target_set = set(target_keys)
            blocks: list[pd.DataFrame] = []
            usecols = ["source_gse_id", "original_cell_id", *chain_columns]
            for item in config["current_atlas_recovery"]["metadata_sources"]:
                metadata_path = resolve(item["path"])
                header = pd.read_csv(metadata_path, nrows=0)
                missing = sorted(set(usecols) - set(header.columns))
                if missing:
                    raise KeyError(
                        f"Recovery metadata {metadata_path} lacks TCR fields: {missing}"
                    )
                for chunk in pd.read_csv(
                    metadata_path,
                    usecols=usecols,
                    dtype="string",
                    chunksize=500_000,
                    low_memory=False,
                ):
                    source = normalize_text(chunk["source_gse_id"])
                    cell = normalize_text(chunk["original_cell_id"])
                    keys = source + "||" + cell
                    selected = keys.isin(target_set)
                    if selected.any():
                        block = chunk.loc[selected, chain_columns].copy()
                        block.index = pd.Index(keys.loc[selected].to_numpy())
                        blocks.append(block)
            if not blocks:
                raise RuntimeError(
                    "No current-atlas boundary cells matched TCR metadata"
                )
            metadata = pd.concat(blocks, axis=0, copy=False)
            if metadata.index.duplicated().any():
                duplicated = int(metadata.index.duplicated(keep=False).sum())
                raise RuntimeError(
                    f"Current-atlas boundary TCR metadata has {duplicated} duplicate keys"
                )
            missing_keys = pd.Index(target_keys).difference(metadata.index)
            if not missing_keys.empty:
                raise RuntimeError(
                    f"Current-atlas boundary TCR metadata missed {missing_keys.size} cells"
                )
            aligned = metadata.reindex(target_keys)
        else:
            source = h5ad_obs_frame(role.path, chain_columns)
            rows = boundary_obs.iloc[local]["input_row"].to_numpy(np.int64)
            aligned = source.iloc[rows][chain_columns].copy()
        for chain, column in zip(("TRA", "TRB", "TRD"), chain_columns, strict=True):
            flags[chain][local] = nonempty_cdr3(aligned[column])
        covered[local] = True
    if not covered.all():
        raise RuntimeError(
            f"Productive-chain extraction missed {(~covered).sum():,} boundary cells"
        )
    temporary = cache.with_name(cache.name + ".partial.npz")
    np.savez_compressed(
        temporary,
        boundary_indices=boundary_indices.astype(np.int64),
        productive_tra=flags["TRA"],
        productive_trb=flags["TRB"],
        productive_trd=flags["TRD"],
        boundary_partitions_sha256=np.asarray(boundary_hash),
        definition=np.asarray(
            "nonempty harmonized productive-filtered chain-specific CDR3"
        ),
    )
    temporary.replace(cache)
    return flags


def plot_boundary_productive_chains(
    sample: dict[str, np.ndarray],
    sample_flags: dict[str, np.ndarray],
    path: Path,
) -> None:
    xy = sample["X_umap"]
    colors = {"TRA": "#2774ae", "TRB": "#16827c", "TRD": "#cf4b32"}
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.7), sharex=True, sharey=True)
    for ax, chain in zip(axes, ("TRA", "TRB", "TRD"), strict=True):
        positive = sample_flags[chain]
        ax.scatter(
            xy[:, 0],
            xy[:, 1],
            s=0.3,
            color="#d4d8dc",
            linewidths=0,
            rasterized=True,
        )
        ax.scatter(
            xy[positive, 0],
            xy[positive, 1],
            s=0.8,
            color=colors[chain],
            linewidths=0,
            rasterized=True,
        )
        clean_umap(ax, f"Productive {chain}  n={positive.sum():,}")
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(
        "Productive TCR-chain evidence inside the NK-like boundary",
        fontsize=12,
        y=1.01,
    )
    save_figure(fig, path)


def plot_transition(path_csv: Path, path_png: Path) -> None:
    table = pd.read_csv(path_csv, index_col=0)
    normalized = table.div(table.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    fig, ax = plt.subplots(figsize=(11.5, 8.0))
    sns.heatmap(
        normalized,
        cmap="rocket_r",
        vmin=0,
        vmax=max(0.25, float(normalized.to_numpy().max())),
        cbar_kws={"label": "Fraction of first-pass cluster"},
        ax=ax,
    )
    ax.set_title(
        "First-pass cluster redistribution after T/NK-restricted HVG selection"
    )
    ax.set_xlabel("Refined cluster")
    ax.set_ylabel("First-pass cluster")
    save_figure(fig, path_png)


def format_table(frame: pd.DataFrame, max_rows: int = 40) -> str:
    display = frame.head(max_rows).copy()
    for column in display.select_dtypes(include=["float"]).columns:
        display[column] = display[column].map(lambda value: f"{value:.4f}")
    return display.to_html(index=False, border=0, classes="data-table")


def write_report_pdf(index: Path, output: Path) -> bool:
    candidates = [
        Path("/usr/bin/google-chrome"),
        Path("/usr/bin/google-chrome-stable"),
        Path("/usr/bin/chromium"),
        Path("/usr/bin/chromium-browser"),
    ]
    chrome = next((path for path in candidates if path.exists()), None)
    if chrome is None:
        return False
    command = [
        str(chrome),
        "--headless",
        "--disable-gpu",
        "--no-sandbox",
        "--print-to-pdf-no-header",
        f"--print-to-pdf={output}",
        index.resolve().as_uri(),
    ]
    completed = subprocess.run(
        command, check=False, capture_output=True, text=True, timeout=300
    )
    return (
        completed.returncode == 0
        and output.is_file()
        and output.stat().st_size > 10_000
    )


def report_stage(config: dict[str, Any], paths: dict[str, Path]) -> dict[str, Any]:
    plotting_defaults()
    prepare = read_json(paths["log"] / "prepare_summary.json")
    fit = read_json(paths["log"] / "fit_summary.json")
    cluster = read_json(paths["log"] / "cluster_summary.json")
    boundary_summary = read_json(paths["log"] / "boundary_summary.json")
    backed = ad.read_h5ad(paths["staged_h5ad"], backed="r")
    obs = backed.obs.copy()
    backed.file.close()
    saved = np.load(paths["partitions"], allow_pickle=False)
    names = saved["names"].astype(str).tolist()
    review = str(config["clustering"]["review_run"])
    labels = saved["labels"][:, names.index(review)].astype(np.int32)
    sample_npz = np.load(paths["umap_sample"], allow_pickle=False)
    sample = {key: np.asarray(sample_npz[key]) for key in sample_npz.files}
    marker = sparse.load_npz(paths["marker_counts"]).tocsr()
    marker_table = pd.read_csv(paths["table"] / "marker_panel.csv")
    genes = marker_table["gene"].astype(str).tolist()
    if marker.shape != (obs.shape[0], len(genes)):
        raise RuntimeError("Marker matrix is inconsistent at report stage")
    sample_expression = np.log1p(marker[sample["sample_indices"]].toarray()).astype(
        np.float32
    )

    boundary_saved = np.load(paths["boundary_partitions"], allow_pickle=False)
    boundary_indices = boundary_saved["boundary_indices"].astype(np.int64)
    boundary_names = boundary_saved["names"].astype(str).tolist()
    boundary_review = str(config["nk_boundary"]["review_run"])
    if boundary_review not in boundary_names:
        raise KeyError(f"Boundary review run is absent: {boundary_review}")
    boundary_labels = boundary_saved["labels"][
        :, boundary_names.index(boundary_review)
    ].astype(np.int32)
    boundary_stability_pairs, boundary_stability = boundary_stability_table(
        boundary_names, boundary_saved["labels"].astype(np.int32)
    )
    boundary_stability_pairs.to_csv(
        paths["table"] / "nk_boundary_seed_stability_pairs.csv", index=False
    )
    boundary_stability.to_csv(
        paths["table"] / "nk_boundary_seed_stability_summary.csv", index=False
    )
    review_resolution = float(
        re.search(r"boundary_r([0-9.]+)_", boundary_review).group(1)
    )
    review_stability = boundary_stability.loc[
        np.isclose(boundary_stability["resolution"], review_resolution)
    ].iloc[0]
    boundary_sample_npz = np.load(paths["boundary_umap_sample"], allow_pickle=False)
    boundary_sample = {
        key: np.asarray(boundary_sample_npz[key]) for key in boundary_sample_npz.files
    }
    productive_chain_flags = boundary_productive_chain_flags(
        config, paths, obs, boundary_indices
    )
    sample_chain_flags = {
        chain: values[boundary_sample["sample_local_indices"]]
        for chain, values in productive_chain_flags.items()
    }
    productive_ab = productive_chain_flags["TRA"] | productive_chain_flags["TRB"]
    sample_productive_ab = sample_chain_flags["TRA"] | sample_chain_flags["TRB"]
    chain_summary_records = [
        {
            "chain": chain,
            "n_boundary_productive": int(values.sum()),
            "fraction_boundary_productive": float(values.mean()),
            "n_umap_sample_productive": int(sample_chain_flags[chain].sum()),
            "fraction_umap_sample_productive": float(sample_chain_flags[chain].mean()),
        }
        for chain, values in productive_chain_flags.items()
    ]
    chain_summary_records.append(
        {
            "chain": "TRA_or_TRB",
            "n_boundary_productive": int(productive_ab.sum()),
            "fraction_boundary_productive": float(productive_ab.mean()),
            "n_umap_sample_productive": int(sample_productive_ab.sum()),
            "fraction_umap_sample_productive": float(sample_productive_ab.mean()),
        }
    )
    chain_summary = pd.DataFrame(chain_summary_records)
    chain_summary.to_csv(
        paths["table"] / "nk_boundary_productive_chain_counts.csv", index=False
    )
    chain_audit = pd.DataFrame(
        {
            "boundary_cluster": boundary_labels,
            "source_gse_id": obs.iloc[boundary_indices]["source_gse_id"]
            .astype(str)
            .to_numpy(),
            **{
                f"productive_{chain.lower()}": values
                for chain, values in productive_chain_flags.items()
            },
            "productive_ab": productive_ab,
        }
    )
    chain_by_cluster = (
        chain_audit.groupby("boundary_cluster", observed=True)
        .agg(
            n_cells=("boundary_cluster", "size"),
            n_productive_tra=("productive_tra", "sum"),
            n_productive_trb=("productive_trb", "sum"),
            n_productive_trd=("productive_trd", "sum"),
            n_productive_ab=("productive_ab", "sum"),
        )
        .reset_index()
    )
    for chain in ("tra", "trb", "trd", "ab"):
        chain_by_cluster[f"fraction_productive_{chain}"] = (
            chain_by_cluster[f"n_productive_{chain}"] / chain_by_cluster["n_cells"]
        )
    chain_by_cluster.to_csv(
        paths["table"] / "nk_boundary_productive_chain_counts_by_cluster.csv",
        index=False,
    )
    chain_by_source = (
        chain_audit.groupby("source_gse_id", observed=True)
        .agg(
            n_cells=("source_gse_id", "size"),
            n_productive_tra=("productive_tra", "sum"),
            n_productive_trb=("productive_trb", "sum"),
            n_productive_trd=("productive_trd", "sum"),
            n_productive_ab=("productive_ab", "sum"),
        )
        .reset_index()
    )
    chain_by_source.to_csv(
        paths["table"] / "nk_boundary_productive_chain_counts_by_source.csv",
        index=False,
    )
    boundary_marker = marker[boundary_indices].tocsr()
    boundary_sample_expression = np.log1p(
        marker[boundary_sample["sample_global_indices"]].toarray()
    ).astype(np.float32)
    boundary_marker_table = cluster_marker_table(
        boundary_labels, boundary_marker, genes
    )
    boundary_marker_table.to_csv(
        paths["table"] / "nk_boundary_cluster_marker_detection.csv", index=False
    )
    boundary_counts = pd.read_csv(paths["table"] / "nk_boundary_cluster_counts.csv")
    boundary_counts["n_sidecar_productive_tcr"] = boundary_counts["n_productive_tcr"]
    boundary_counts["sidecar_productive_tcr_fraction"] = (
        boundary_counts["n_productive_tcr"] / boundary_counts["n_cells"]
    )
    boundary_counts = boundary_counts.merge(
        chain_by_cluster.drop(columns="n_cells"),
        on="boundary_cluster",
        how="left",
        validate="1:1",
    )
    review_annotations = {
        int(key): str(value)
        for key, value in config["nk_boundary"]["review_annotations"].items()
    }
    observed_boundary_clusters = set(boundary_counts["boundary_cluster"].astype(int))
    if set(review_annotations) != observed_boundary_clusters:
        raise RuntimeError(
            "Boundary review annotations do not exactly cover observed clusters"
        )
    boundary_counts["review_annotation"] = boundary_counts["boundary_cluster"].map(
        review_annotations
    )
    boundary_counts.to_csv(
        paths["table"] / "nk_boundary_cluster_review.csv", index=False
    )
    boundary_source = pd.read_csv(
        paths["table"] / "nk_boundary_cluster_counts_by_source.csv"
    )
    core_clusters = sorted(
        int(cluster)
        for cluster in config["nk_boundary"]["previous_candidate_core_clusters"]
    )
    core_rows = boundary_counts.loc[
        boundary_counts["boundary_cluster"].isin(core_clusters)
    ]
    n_primary_nk_total = int(obs["primary_nk_anchor"].astype(bool).sum())
    n_core_cells = int(core_rows["n_cells"].sum())
    n_core_primary_nk = int(core_rows["n_primary_nk"].sum())
    n_core_productive_tcr = int(core_rows["n_productive_tcr"].sum())
    core_primary_nk_recall = n_core_primary_nk / n_primary_nk_total
    core_mask = np.isin(boundary_labels, core_clusters)
    n_core_productive_ab = int(productive_ab[core_mask].sum())
    n_core_productive_trd = int(productive_chain_flags["TRD"][core_mask].sum())
    boundary_primary_nk = (
        obs.iloc[boundary_indices]["primary_nk_anchor"].astype(bool).to_numpy()
    )
    n_primary_nk_productive_ab = int((boundary_primary_nk & productive_ab).sum())

    cluster_marker = cluster_marker_table(labels, marker, genes)
    cluster_marker.to_csv(
        paths["table"] / "refined_cluster_marker_detection.csv", index=False
    )
    cluster_summary = pd.read_csv(paths["table"] / "refined_cluster_counts.csv")
    cluster_summary["primary_nk_fraction"] = (
        cluster_summary["n_primary_nk"] / cluster_summary["n_cells"]
    )
    cluster_summary["productive_tcr_fraction"] = (
        cluster_summary["n_productive_tcr"] / cluster_summary["n_cells"]
    )
    phenotype_cluster = (
        pd.DataFrame(
            {
                "cluster": labels,
                "phenotype": obs["trdc_delta_v_cd3_phenotype"].astype(str).to_numpy(),
            }
        )
        .groupby(["cluster", "phenotype"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    phenotype_cluster.to_csv(
        paths["table"] / "refined_phenotype_counts_by_cluster.csv", index=False
    )

    plot_clusters(sample, paths["figure"] / "umap_refined_clusters.png")
    plot_categories(
        sample,
        sample["phenotype"],
        "TRDC, delta-V, and CD3 phenotype on the refined T/NK UMAP",
        paths["figure"] / "umap_trdc_delta_v_cd3_phenotype.png",
    )
    plot_categories(
        sample,
        sample["source_gse_id"],
        "Dataset source on the refined T/NK UMAP",
        paths["figure"] / "umap_source_gse_id.png",
    )
    plot_features(
        sample,
        sample_expression,
        genes,
        [
            "CD3D",
            "CD3E",
            "CD3G",
            "CD247",
            "TRAC",
            "TRBC1",
            "TRBC2",
            "TRDC",
            "TRGC1",
            "TRGC2",
            "LCK",
            "TRAT1",
        ],
        "CD3 and TCR constant-chain expression",
        paths["figure"] / "feature_umap_cd3_tcr_constant.png",
    )
    plot_features(
        sample,
        sample_expression,
        genes,
        config["tnk_gate"]["gd_lineage_genes"],
        "Gamma-delta TCR expression",
        paths["figure"] / "feature_umap_gdt_genes.png",
    )
    plot_features(
        sample,
        sample_expression,
        genes,
        [
            *config["tnk_gate"]["nk_receptor_genes"],
            *config["tnk_gate"]["nk_adaptor_genes"],
        ],
        "NK-lineage receptor and adaptor context",
        paths["figure"] / "feature_umap_nk_lineage.png",
    )
    plot_features(
        sample,
        sample_expression,
        genes,
        config["tnk_gate"]["myeloid_context_genes"],
        "Myeloid contamination context",
        paths["figure"] / "feature_umap_myeloid_context.png",
    )
    heatmap_genes = [
        "CD3D",
        "CD3E",
        "CD3G",
        "CD247",
        "TRAC",
        "TRBC1",
        "TRBC2",
        "TRDC",
        "TRGC1",
        "TRGC2",
        "TRDV1",
        "TRDV2",
        "TRDV3",
        "KLRD1",
        "NCR1",
        "FCGR3A",
        "KLRC1",
        "KLRF1",
        "S1PR5",
        "FCER1G",
        "TYROBP",
        "LST1",
        "AIF1",
        "CTSS",
        "LILRB1",
        "LILRB2",
        "NKG7",
        "GNLY",
        "PRF1",
    ]
    plot_cluster_heatmap(
        cluster_marker,
        heatmap_genes,
        paths["figure"] / "refined_cluster_marker_detection_heatmap.png",
    )
    plot_transition(
        paths["table"] / "first_pass_to_refined_cluster_counts.csv",
        paths["figure"] / "first_pass_to_refined_cluster_heatmap.png",
    )
    plot_boundary_clusters(
        boundary_sample, paths["figure"] / "nk_boundary_umap_clusters.png"
    )
    plot_categories(
        boundary_sample,
        boundary_sample["phenotype"],
        "TRDC, delta-V, and CD3 phenotype inside the NK-like boundary",
        paths["figure"] / "nk_boundary_umap_phenotype.png",
    )
    plot_categories(
        boundary_sample,
        boundary_sample["parent_refined_cluster"].astype(str),
        "Parent refined cluster inside the NK-like boundary",
        paths["figure"] / "nk_boundary_umap_parent_cluster.png",
    )
    plot_features(
        boundary_sample,
        boundary_sample_expression,
        genes,
        [
            "CD3D",
            "CD3E",
            "CD3G",
            "TRAC",
            "TRBC1",
            "TRBC2",
            "TRDC",
            "TRDV1",
            "TRDV2",
            "KLRD1",
            "FCER1G",
            "TYROBP",
        ],
        "T, gamma-delta, and NK evidence inside the NK-like boundary",
        paths["figure"] / "nk_boundary_feature_umap.png",
    )
    plot_boundary_productive_chains(
        boundary_sample,
        sample_chain_flags,
        paths["figure"] / "nk_boundary_productive_tcr_umap.png",
    )
    plot_cluster_heatmap(
        boundary_marker_table,
        heatmap_genes,
        paths["figure"] / "nk_boundary_marker_detection_heatmap.png",
        title="Raw-UMI marker detection by NK-boundary subcluster",
        ylabel="Boundary cluster",
    )
    plot_boundary_evidence(
        boundary_counts, paths["figure"] / "nk_boundary_evidence.png"
    )
    plot_boundary_source_heatmap(
        boundary_source, paths["figure"] / "nk_boundary_source_heatmap.png"
    )

    gate_reason = pd.read_csv(paths["table"] / "tnk_gate_reason_counts.csv")
    gate_source = pd.read_csv(paths["table"] / "tnk_gate_counts_by_source.csv")
    phenotype = pd.read_csv(paths["table"] / "tnk_gate_phenotype_counts.csv")
    source_kept = gate_source.loc[
        gate_source["tnk_gate_keep"].astype(bool)
    ].sort_values("n_cells", ascending=False)
    boundary_display = boundary_counts[
        [
            "boundary_cluster",
            "review_annotation",
            "n_cells",
            "n_primary_nk",
            "n_sidecar_productive_tcr",
            "n_productive_ab",
            "n_productive_trd",
            "dominant_source",
            "primary_nk_fraction",
            "fraction_productive_ab",
            "dominant_source_fraction",
        ]
    ].sort_values("boundary_cluster")
    prefix = (
        "../../Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_tnk_reintegration"
    )
    css = """
@page { size: A4 landscape; margin: 10mm; }
body { font-family: Arial, sans-serif; color:#17212b; max-width:1280px; margin:0 auto; padding:24px; line-height:1.45; }
h1 { font-size:28px; margin-bottom:4px; } h2 { font-size:19px; border-bottom:2px solid #d9e1e8; padding-bottom:5px; margin-top:30px; }
.status { border-left:5px solid #147d64; background:#eef8f5; padding:12px 16px; margin:16px 0; }
.warning { border-left:5px solid #c47a00; background:#fff7e8; padding:12px 16px; margin:16px 0; }
.grid { display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:16px; }
figure { margin:8px 0; break-inside:avoid; } img { width:100%; height:auto; border:1px solid #d8dee4; }
.data-table { border-collapse:collapse; width:100%; font-size:8.5px; } .data-table th,.data-table td { border:1px solid #cfd7df; padding:4px 6px; text-align:right; }
.data-table th:first-child,.data-table td:first-child { text-align:left; } .table-wrap { overflow-x:auto; break-inside:avoid; }
code { background:#eef1f4; padding:1px 4px; } .small { color:#53606d; font-size:9px; }
"""
    report_html = f"""<!doctype html><html><head><meta charset="utf-8"><title>gdTAI V4.2 T/NK-restricted reintegration</title><style>{css}</style></head><body>
<h1>gdTAI V4.2 T/NK-restricted reintegration</h1>
<p class="small">Generated 2026-08-17. This is a development-sidecar analysis, not a classifier result or model promotion.</p>
<div class="status"><b>Refined T/NK pool:</b> {prepare["n_retained_cells"]:,} of {prepare["n_input_cells"]:,} cells retained ({prepare["retained_fraction"]:.2%}); {prepare["n_removed_cells"]:,} cells removed. scVI used {prepare["n_hvgs"]:,} features and GPU only.</div>
<div class="warning"><b>Updated conclusion: no NK training core is accepted.</b> Harmonized productive-chain metadata identify TRA or TRB in {productive_ab.sum():,}/{boundary_summary["n_boundary_cells"]:,} boundary cells ({productive_ab.mean():.2%}) and in {n_primary_nk_productive_ab:,}/{boundary_summary["n_primary_nk"]:,} primary NK annotation anchors ({n_primary_nk_productive_ab / boundary_summary["n_primary_nk"]:.2%}). The previous clusters 0/3/5 candidate-core interpretation is withdrawn.</div>
<h2>1. What changed</h2>
<p>The inputs were already high-recall T/NK candidate objects, but they retained ambiguous off-target cells. This second pass preserves all non-doublet productive-TCR cells and primary NK anchors, then requires direct T-lineage, gamma-delta, NK-receptor, or multigene NK-adaptor-plus-cytotoxic evidence. Cells with broad myeloid evidence and no direct T/gdT/NK-receptor support are removed.</p>
<div class="warning"><b>TRDC+ / delta-V- / weak-CD3 is not used as an NK label.</b> It is retained as a phenotype for unsupervised review. Dropout, ambient RNA, and true NK biology remain competing explanations.</div>
<div class="table-wrap">{format_table(gate_reason)}</div>
<h2>2. Recomputed HVGs and integration</h2>
<p>HVGs were recomputed only after T/NK subsetting with source-balanced sampling. TCR V/J/D genes were excluded from the HVG ranking to limit clonotype-driven structure. A fixed 27-gene CD3/TCR-constant and T/NK lineage panel was forced into the 4,000-feature set. The fit completed in {fit["elapsed_seconds"] / 60:.1f} minutes on {html.escape(fit["gpu"]["device_name"])}.</p>
<div class="grid"><figure><img src="{prefix}/umap_refined_clusters.png"><figcaption>Unsupervised refined clusters.</figcaption></figure><figure><img src="{prefix}/umap_source_gse_id.png"><figcaption>Dataset mixing remains visible and must be considered in annotation.</figcaption></figure></div>
<h2>3. TRDC, delta-V, and CD3 phenotype</h2>
<figure><img src="{prefix}/umap_trdc_delta_v_cd3_phenotype.png"><figcaption>The proposed NK-like phenotype is displayed but not treated as truth.</figcaption></figure>
<div class="table-wrap">{format_table(phenotype)}</div>
<h2>4. CD3 and TCR genes</h2>
<figure><img src="{prefix}/feature_umap_cd3_tcr_constant.png"></figure>
<figure><img src="{prefix}/feature_umap_gdt_genes.png"></figure>
<h2>5. NK and contamination context</h2>
<div class="grid"><figure><img src="{prefix}/feature_umap_nk_lineage.png"></figure><figure><img src="{prefix}/feature_umap_myeloid_context.png"></figure></div>
<h2>6. Cluster-level evidence</h2>
<figure><img src="{prefix}/refined_cluster_marker_detection_heatmap.png"></figure>
<div class="table-wrap">{format_table(cluster_summary.sort_values("primary_nk_fraction", ascending=False))}</div>
<h2>7. Second-pass NK-boundary clustering</h2>
<p>Refined clusters 9 and 18 contain {boundary_summary["n_boundary_cells"]:,} cells and recover {boundary_summary["n_primary_nk"]:,}/{n_primary_nk_total:,} primary NK annotation anchors ({boundary_summary["primary_nk_recall"]:.2%}). The original sidecar mask counted only {boundary_summary["n_productive_tcr"]:,} productive-T controls because many chain fields had not been propagated into <code>TNK_cleaned.h5ad</code>. The harmonized metadata audit now identifies {productive_ab.sum():,} cells with productive TRA or TRB ({productive_ab.mean():.2%}) and {productive_chain_flags["TRD"].sum():,} with productive TRD. Nine GPU Leiden runs were performed inside this boundary; the frozen review partition <code>{html.escape(boundary_review)}</code> contains {boundary_summary["n_review_clusters"]} subclusters. At resolution {review_resolution:g}, the mean adjusted Rand index across three seed pairs is {review_stability["mean_adjusted_rand_index"]:.3f} (minimum {review_stability["min_adjusted_rand_index"]:.3f}).</p>
<div class="grid"><figure><img src="{prefix}/nk_boundary_umap_clusters.png"><figcaption>Second-pass unsupervised subclusters.</figcaption></figure><figure><img src="{prefix}/nk_boundary_umap_phenotype.png"><figcaption>The weak-CD3 phenotype is enriched but is not an NK label.</figcaption></figure></div>
<figure><img src="{prefix}/nk_boundary_productive_tcr_umap.png"><figcaption>Chain-specific productive TCR evidence is defined by a nonempty harmonized productive-filtered CDR3, not by RNA expression. Absence of a plotted chain is not informative in libraries without the corresponding V(D)J assay.</figcaption></figure>
<div class="table-wrap">{format_table(chain_summary)}</div>
<div class="warning"><b>The productive-chain overlay overturns the earlier boundary interpretation.</b> Harmonized TRA/TRB calls are a strict superset of the source-H5AD calls: all 6,307 source-H5AD-positive boundary cells agree, while 364,091 additional current-atlas cells have productive TRA/TRB metadata that had not been propagated into the H5AD. In addition, {n_primary_nk_productive_ab:,}/{boundary_summary["n_primary_nk"]:,} ({n_primary_nk_productive_ab / boundary_summary["n_primary_nk"]:.2%}) primary NK annotation anchors carry productive TRA or TRB, so those annotations cannot independently validate NK identity.</div>
<div class="grid"><figure><img src="{prefix}/nk_boundary_evidence.png"></figure><figure><img src="{prefix}/nk_boundary_source_heatmap.png"></figure></div>
<div class="table-wrap">{format_table(boundary_stability)}</div>
<h2>8. Revised NK-boundary interpretation</h2>
<p>Boundary clusters {", ".join(map(str, core_clusters))} were initially considered a transcriptomic NK-like core because they combine weak <code>CD3D/CD3G</code>, high NK receptors, <code>FCER1G/TYROBP</code>, and cytotoxic genes with low <code>LST1/AIF1</code>. The productive-chain audit rejects that interpretation: {n_core_productive_ab:,}/{n_core_cells:,} cells ({n_core_productive_ab / n_core_cells:.2%}) have productive TRA or TRB, while only {n_core_productive_trd:,} have productive TRD. The previous candidate-core designation is therefore withdrawn; these clusters are NK-like transcriptomically but TCR-alpha-beta dominated.</p>
<div class="warning"><b>No boundary subcluster is promoted as NK training truth.</b> The transcriptomic similarity of cytotoxic alpha-beta T, gamma-delta T, and NK cells cannot be resolved by <code>TRDC</code>, cytotoxic genes, <code>FCER1G/TYROBP</code>, or any single marker. The TCR metadata conflict and the compromised primary NK anchors must be repaired before a new classifier iteration.</div>
<div class="table-wrap">{format_table(boundary_display)}</div>
<figure><img src="{prefix}/nk_boundary_marker_detection_heatmap.png"></figure>
<figure><img src="{prefix}/nk_boundary_feature_umap.png"></figure>
<h2>9. Difference from the first pass</h2>
<figure><img src="{prefix}/first_pass_to_refined_cluster_heatmap.png"><figcaption>Rows are first-pass clusters and columns are refined clusters; values are row fractions.</figcaption></figure>
<h2>10. Retained cells by source</h2><div class="table-wrap">{format_table(source_kept, max_rows=80)}</div>
<h2>11. Interpretation boundary</h2>
<p>The integration identifies coherent NK-like and T-like transcriptomic neighborhoods, but the harmonized chain audit shows that the selected boundary is predominantly productive TCR-alpha-beta. The existing primary NK annotation anchors are also compromised by productive TRA/TRB evidence. No boundary cluster is suitable as NK training truth until TCR metadata are propagated consistently, anchor conflicts are resolved, and expression-independent NK labels are reconstructed.</p>
<h2>Reproducibility</h2><ul>
<li>Config SHA-256: <code>{html.escape(config["_config_sha256"])}</code></li><li>Script SHA-256: <code>{sha256_file(SCRIPT_PATH)}</code></li>
<li>Staged H5AD SHA-256: <code>{prepare["staged_h5ad_sha256"]}</code></li><li>Latent SHA-256: <code>{fit["latent_sha256"]}</code></li>
<li>Partitions SHA-256: <code>{cluster["partitions_sha256"]}</code></li><li>Boundary partitions SHA-256: <code>{boundary_summary["boundary_partitions_sha256"]}</code></li><li>Source H5AD mutation: <b>none</b>; classifier fitting: <b>none</b>; GitHub push: <b>none</b>.</li></ul>
</body></html>"""
    index = paths["static"] / "index.html"
    index.write_text(report_html, encoding="utf-8")
    pdf = paths["static"] / "gdtai_v4_2_tnk_reintegration_report.pdf"
    pdf_ok = write_report_pdf(index, pdf)
    summary = {
        "stage": "report",
        "status": "PASS_AUDIT_NK_CORE_REJECTED",
        "n_input_cells": prepare["n_input_cells"],
        "n_retained_cells": prepare["n_retained_cells"],
        "n_removed_cells": prepare["n_removed_cells"],
        "n_refined_clusters": cluster["n_review_clusters"],
        "n_boundary_clusters": boundary_summary["n_review_clusters"],
        "n_nk_boundary_cells": boundary_summary["n_boundary_cells"],
        "n_boundary_productive_tra": int(productive_chain_flags["TRA"].sum()),
        "n_boundary_productive_trb": int(productive_chain_flags["TRB"].sum()),
        "n_boundary_productive_trd": int(productive_chain_flags["TRD"].sum()),
        "previous_candidate_core_rejected": True,
        "n_previous_candidate_core_cells": n_core_cells,
        "n_previous_candidate_core_primary_nk": n_core_primary_nk,
        "previous_candidate_core_primary_nk_recall": core_primary_nk_recall,
        "n_previous_candidate_core_sidecar_productive_tcr": n_core_productive_tcr,
        "n_previous_candidate_core_harmonized_productive_ab": n_core_productive_ab,
        "n_primary_nk_anchor_harmonized_productive_ab": n_primary_nk_productive_ab,
        "boundary_review_mean_seed_ari": float(
            review_stability["mean_adjusted_rand_index"]
        ),
        "n_umap_sample": cluster["n_umap_sample"],
        "html_report": str(index),
        "pdf_report": str(pdf) if pdf_ok else "not_rendered",
        "pdf_rendered": pdf_ok,
        "source_h5ad_mutation_performed": False,
        "classifier_fitting_performed": False,
        "github_push_performed": False,
        "config_sha256": config["_config_sha256"],
        "script_sha256": sha256_file(SCRIPT_PATH),
    }
    atomic_write_json(paths["log"] / "summary.json", summary)
    return summary


def validate_stage(config: dict[str, Any], paths: dict[str, Path]) -> dict[str, Any]:
    roles = development_roles(config)
    first = config["first_pass"]
    missing = [
        path
        for path in first.values()
        if isinstance(path, str)
        and path.endswith((".h5ad", ".npy", ".npz"))
        and not Path(path).is_file()
    ]
    genes = marker_genes(config)
    common = set.intersection(
        *(
            set(h5ad_var_names(row.path).astype(str))
            for row in roles.itertuples(index=False)
        )
    )
    absent = sorted(set(genes) - common)
    if missing or absent:
        raise RuntimeError(
            f"Validation failed; missing artifacts={missing}; absent markers={absent}"
        )
    return {
        "stage": "validate",
        "status": "PASS",
        "protocol_version": config["protocol_version"],
        "development_cohorts": int(roles.shape[0]),
        "development_cells_expected": int(roles["expected_cells"].sum()),
        "marker_genes": len(genes),
        "common_genes": len(common),
        "forced_lineage_genes": len(config["hvg"]["forced_lineage_genes"]),
        "config_sha256": config["_config_sha256"],
        "script_sha256": sha256_file(SCRIPT_PATH),
        "ssd_root": str(paths["ssd"]),
    }


def main() -> None:
    args = parse_args()
    config = load_runtime_config(args.config)
    paths = stage_paths(config)
    validation = validate_stage(config, paths)
    if args.stage == "validate":
        print(json.dumps(validation, indent=2, sort_keys=True))
        return
    stages = (
        [args.stage]
        if args.stage != "all"
        else ["prepare", "fit", "cluster", "boundary", "report"]
    )
    results: dict[str, Any] = {"validation": validation}
    stage_summary_files = {
        "prepare": paths["log"] / "prepare_summary.json",
        "fit": paths["log"] / "fit_summary.json",
        "cluster": paths["log"] / "cluster_summary.json",
        "boundary": paths["log"] / "boundary_summary.json",
        "report": paths["log"] / "summary.json",
    }
    for stage_name, summary_path in stage_summary_files.items():
        if summary_path.is_file():
            results[stage_name] = {
                "recovered_from_stage_summary": True,
                **read_json(summary_path),
            }
    for stage in stages:
        resources = ensure_paths_and_resources(config, paths, stage)
        if stage == "prepare":
            result = prepare_stage(config, paths, args.overwrite)
        elif stage == "fit":
            result = fit_stage(config, paths, args.overwrite)
        elif stage == "cluster":
            result = cluster_stage(config, paths, args.overwrite)
        elif stage == "boundary":
            result = boundary_stage(config, paths, args.overwrite)
        elif stage == "report":
            result = report_stage(config, paths)
        else:
            raise AssertionError(stage)
        results[stage] = {"resources": resources, **result}
        print(json.dumps(result, indent=2, sort_keys=True))
    atomic_write_json(paths["log"] / "run_summary.json", results)


if __name__ == "__main__":
    main()
