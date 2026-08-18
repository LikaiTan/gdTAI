#!/usr/bin/env python3
"""Rebuild the full T/NK atlas from frozen historical and extension inputs.

The workflow is deliberately stage-addressable. Source H5AD files are opened
read-only; all large intermediates and the final H5AD are written to the
approved SSD mirror. The validated TCR-repair sidecar is not consumed here.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import html
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Iterable

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

os.environ.setdefault("ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS", "1")

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from workflows.gdtai.gdtai_v4_2_integration_core import (
    attach_recovery_metadata,
    axis_names,
    obs_values,
    read_sparse_rows_genes,
)
from workflows.integration.phase4_gdt_module_scoring import (
    CHUNK_SIZE,
    CTRL_SIZE,
    N_BINS,
    PACKAGE_ZIP,
    PHASE4_SCALED_SCORE_COLUMNS,
    PHASE4_SCORE_COLUMNS,
    RANDOM_STATE,
    TARGET_SUM,
    add_scaled_scores,
    append_obs_columns_in_place,
    compute_gene_means,
    compute_scores,
    find_module_genes,
    pick_control_genes,
)


DEFAULT_CONFIG = PROJECT_ROOT / "configs/datasets/full_atlas_rebuild.json"
SCRIPT_PATH = Path(__file__).resolve()
HIGH_CARDINALITY_TEXT = {
    "source_obs_name",
    "original_cell_id",
    "barcode",
    "barcode_core",
    *{f"{chain}_{field}" for chain in ("TRA", "TRB", "TRG", "TRD") for field in ("cdr3", "cdr3_nt", "clone_id")},
}
TCR_HVG_PATTERN = re.compile(r"^TR[ABGD][VDJ][0-9A-Z-]*$", re.IGNORECASE)
OTHER_HVG_PATTERN = re.compile(
    r"^(MT-|RPS|RPL|MRPS|MRPL|LINC|MIR|MIRLET|SNHG|SNORA|SNORD|SCARNA|RNU|RN7|RNA5|RNA18|RNA28|YRNA|Y_RNA|AC[0-9]+|AL[0-9]+|AP[0-9]+|BX[0-9]+|RP[0-9]+-|IG[HKL][VDJ])",
    re.IGNORECASE,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stage",
        choices=["preflight", "prepare", "fit", "embed", "assemble", "score", "report", "all"],
        default="preflight",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--send-email", action="store_true")
    return parser.parse_args()


def resolve(path: str | Path) -> Path:
    value = Path(path)
    return value if value.is_absolute() else PROJECT_ROOT / value


def read_config(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        config = json.load(handle)
    if len(config["inputs"]) != int(config["expected_input_cohorts"]):
        raise ValueError("Input-cohort count differs from the frozen contract")
    if sum(int(row["expected_cells_effective"]) for row in config["inputs"]) != int(config["expected_total_cells"]):
        raise ValueError("Input cell counts do not sum to expected_total_cells")
    return config


def sha256_file(path: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".partial")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def atomic_h5ad(adata: ad.AnnData, path: Path, compression: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    temporary.unlink(missing_ok=True)
    adata.write_h5ad(temporary, compression=compression)
    temporary.replace(path)


def output_paths(config: dict[str, Any]) -> dict[str, Path]:
    out = config["outputs"]
    result = {key: resolve(value) for key, value in out.items()}
    for key in ("ssd_root", "prepared_dir", "table_dir", "figure_dir", "log_dir", "report_dir"):
        result[key].mkdir(parents=True, exist_ok=True)
    return result


def input_path(row: dict[str, Any]) -> Path:
    return resolve(row["path"]).resolve()


def input_matrix_key(row: dict[str, Any]) -> str:
    mode = row["matrix_mode"]
    if mode == "layer_counts":
        return "layers/counts"
    if mode in {"X_counts", "seurat_lognormalize_inverse"}:
        return "X"
    raise ValueError(f"Unknown matrix_mode: {mode}")


def matrix_group(handle: h5py.File, key: str) -> h5py.Group:
    value: Any = handle
    for component in key.split("/"):
        value = value[component]
    if not isinstance(value, h5py.Group):
        raise TypeError(f"Expected sparse H5AD matrix group at {key}")
    return value


def matrix_shape(group: h5py.Group) -> tuple[int, int]:
    if "shape" in group.attrs:
        return tuple(int(value) for value in group.attrs["shape"])
    return int(group["indptr"].shape[0] - 1), int(group.attrs.get("n_cols", 0))


def source_values(path: Path, fallback: str) -> np.ndarray:
    with h5py.File(path, "r") as handle:
        if "source_gse_id" not in handle["obs"]:
            return np.repeat(fallback, len(axis_names(handle, "obs"))).astype(object)
        values = np.asarray(obs_values(handle, "source_gse_id"), dtype=object)
    values = np.asarray([str(value).strip() if value is not None else "" for value in values], dtype=object)
    values[values == ""] = fallback
    return values


def sampled_matrix_audit(path: Path, key: str, normalized_expected: bool) -> dict[str, Any]:
    with h5py.File(path, "r") as handle:
        group = matrix_group(handle, key)
        data = group["data"]
        n = min(100_000, int(data.shape[0]))
        values = np.asarray(data[:n])
    finite = bool(np.isfinite(values).all())
    integer_fraction = float(np.mean(np.isclose(values, np.rint(values)))) if values.size else 1.0
    if normalized_expected:
        state_ok = finite and integer_fraction < 0.99
    else:
        state_ok = finite and integer_fraction >= 0.999
    return {
        "sampled_values": int(values.size),
        "sampled_finite": finite,
        "sampled_integer_fraction": integer_fraction,
        "matrix_state_ok": state_ok,
    }


def preflight_stage(config_path: Path, config: dict[str, Any], paths: dict[str, Path]) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []
    target_genes: pd.Index | None = None
    total = 0
    all_sources: set[str] = set()
    input_states: list[dict[str, Any]] = []

    recovery = config["core_recovery"]
    checks.append({"check": "config_expected_total", "status": "PASS", "detail": config["expected_total_cells"]})
    free_gib = shutil.disk_usage(paths["ssd_root"]).free / 1024**3
    checks.append(
        {
            "check": "ssd_free_space",
            "status": "PASS" if free_gib >= float(config["resources"]["minimum_ssd_free_gib"]) else "FAIL",
            "detail": f"{free_gib:.1f} GiB free",
        }
    )
    for dependency in [recovery["row_exclusion_manifest"], *[item["path"] for item in recovery["metadata_sources"]]]:
        path = resolve(dependency)
        checks.append({"check": f"dependency::{dependency}", "status": "PASS" if path.exists() else "FAIL", "detail": str(path)})
    if sha256_file(resolve(recovery["row_exclusion_manifest"])) != recovery["row_exclusion_manifest_sha256"]:
        checks.append({"check": "core_exclusion_manifest_sha256", "status": "FAIL", "detail": "checksum mismatch"})
    else:
        checks.append({"check": "core_exclusion_manifest_sha256", "status": "PASS", "detail": recovery["row_exclusion_manifest_sha256"]})
    for item in recovery["metadata_sources"]:
        path = resolve(item["path"])
        observed = sha256_file(path) if path.exists() else "missing"
        checks.append(
            {
                "check": f"core_metadata_sha256::{item['path']}",
                "status": "PASS" if observed == item["sha256"] else "FAIL",
                "detail": observed,
            }
        )

    for row in config["inputs"]:
        cohort = row["cohort_id"]
        path = input_path(row)
        if not path.exists():
            checks.append({"check": f"exists::{cohort}", "status": "FAIL", "detail": str(path)})
            continue
        stat = path.stat()
        input_states.append({"cohort_id": cohort, "path": str(path), "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns})
        observed_sha = sha256_file(path)
        expected_sha = row["expected_sha256"]
        sha_ok = expected_sha != "PENDING_FULL_SHA256" and observed_sha == expected_sha
        checks.append({"check": f"sha256::{cohort}", "status": "PASS" if sha_ok else "FAIL", "detail": observed_sha})
        with h5py.File(path, "r") as handle:
            obs_names = axis_names(handle, "obs")
            var_names = pd.Index(axis_names(handle, "var").astype(str), dtype="string")
            key = input_matrix_key(row)
            shape = matrix_shape(matrix_group(handle, key))
        expected_raw = int(row["expected_cells_raw"])
        checks.append({"check": f"shape::{cohort}", "status": "PASS" if shape[0] == expected_raw else "FAIL", "detail": str(shape)})
        checks.append({"check": f"obs_unique::{cohort}", "status": "PASS" if pd.Index(obs_names).is_unique else "FAIL", "detail": len(obs_names)})
        if target_genes is None:
            target_genes = var_names
            checks.append({"check": "target_gene_count", "status": "PASS" if len(target_genes) == int(config["expected_genes"]) else "FAIL", "detail": len(target_genes)})
            checks.append({"check": "target_gene_unique", "status": "PASS" if target_genes.is_unique else "FAIL", "detail": int(target_genes.duplicated().sum())})
        overlap = int(target_genes.isin(var_names).sum())
        matrix_audit = sampled_matrix_audit(path, key, row["matrix_mode"] == "seurat_lognormalize_inverse")
        checks.append({"check": f"matrix_state::{cohort}", "status": "PASS" if matrix_audit["matrix_state_ok"] else "FAIL", "detail": matrix_audit["sampled_integer_fraction"]})
        if row.get("apply_core_recovery"):
            excluded = set(pd.read_csv(resolve(recovery["row_exclusion_manifest"]), usecols=["obs_name"])["obs_name"].astype(str))
            intersect = int(pd.Index(obs_names.astype(str)).isin(excluded).sum())
            checks.append({"check": "core_exclusions", "status": "PASS" if intersect == int(recovery["expected_intersecting_exclusions"]) else "FAIL", "detail": intersect})
        effective = int(row["expected_cells_effective"])
        total += effective
        fallback = row.get("fallback_source_gse_id", cohort)
        sources = source_values(path, fallback)
        source_unique = sorted(set(sources.astype(str)))
        all_sources.update(source_unique)
        rows.append(
            {
                "cohort_id": cohort,
                "path": str(path),
                "cells_raw": shape[0],
                "cells_effective": effective,
                "genes": shape[1],
                "genes_overlapping_target": overlap,
                "matrix_key": key,
                "matrix_mode": row["matrix_mode"],
                "source_accessions": ";".join(source_unique),
                "n_source_accessions": len(source_unique),
                "atlas_input_role": row["atlas_input_role"],
                "model_evaluation_role_frozen": row["model_evaluation_role_frozen"],
                "sha256": observed_sha,
                **matrix_audit,
            }
        )

    checks.append({"check": "effective_cell_total", "status": "PASS" if total == int(config["expected_total_cells"]) else "FAIL", "detail": total})
    checks.append({"check": "input_cohort_total", "status": "PASS" if len(rows) == int(config["expected_input_cohorts"]) else "FAIL", "detail": len(rows)})
    checks.append({"check": "tcr_sidecar_not_propagated", "status": "PASS" if config["source_policy"]["propagate_repaired_tcr_sidecar"] is False else "FAIL", "detail": config["source_policy"]["tcr_sidecar_status"]})
    checks_df = pd.DataFrame(checks)
    inputs_df = pd.DataFrame(rows)
    checks_df.to_csv(paths["table_dir"] / "preflight_checks.csv", index=False)
    inputs_df.to_csv(paths["table_dir"] / "input_contract.csv", index=False)
    pd.DataFrame(input_states).to_csv(paths["table_dir"] / "source_file_state_before.csv", index=False)
    status = "PASS" if checks_df["status"].eq("PASS").all() else "FAIL"
    summary = {
        "stage": "preflight",
        "status": status,
        "protocol_version": config["protocol_version"],
        "config_sha256": sha256_file(config_path),
        "script_sha256": sha256_file(SCRIPT_PATH),
        "n_checks": int(len(checks_df)),
        "n_passed": int(checks_df["status"].eq("PASS").sum()),
        "n_input_cohorts": int(len(rows)),
        "n_source_accessions": int(len(all_sources)),
        "expected_total_cells": int(config["expected_total_cells"]),
        "expected_genes": int(config["expected_genes"]),
        "tcr_sidecar_propagated": False,
    }
    atomic_json(paths["log_dir"] / "preflight_summary.json", summary)
    if status != "PASS":
        failures = checks_df.loc[checks_df["status"].ne("PASS")].to_dict("records")
        raise RuntimeError(f"Full-atlas preflight failed: {failures[:10]}")
    return summary


def normalized_text(series: pd.Series | Iterable[Any], index: pd.Index) -> pd.Series:
    result = pd.Series(series, index=index, dtype="string").fillna("").str.strip()
    return result.mask(result.str.lower().isin({"nan", "none", "na", "n/a", "null"}), "")


def first_text(obs: pd.DataFrame, names: list[str], default: str = "") -> pd.Series:
    result = pd.Series("", index=obs.index, dtype="string")
    for name in names:
        if name not in obs:
            continue
        values = normalized_text(obs[name], obs.index)
        take = result.eq("") & values.ne("")
        result.loc[take] = values.loc[take]
    return result.mask(result.eq(""), default)


def bool_values(obs: pd.DataFrame, names: list[str]) -> pd.Series:
    result = pd.Series(False, index=obs.index, dtype=bool)
    for name in names:
        if name not in obs:
            continue
        values = obs[name]
        if pd.api.types.is_bool_dtype(values):
            result |= values.fillna(False).astype(bool)
        else:
            text = normalized_text(values, obs.index).str.lower()
            result |= text.isin({"true", "1", "yes", "y", "productive", "paired"})
    return result


def standardize_obs(obs: pd.DataFrame, row: dict[str, Any]) -> pd.DataFrame:
    cohort = row["cohort_id"]
    fallback = row.get("fallback_source_gse_id", cohort)
    source_obs_name = pd.Series(obs.index.astype(str), index=obs.index, dtype="string")
    original = first_text(obs, ["original_cell_id", "cell_id", "raw_obs_name"], default="")
    original = original.mask(original.eq(""), source_obs_name)
    source = first_text(obs, ["source_gse_id", "source_accession", "gse_id", "GSE"], default=fallback)
    out = pd.DataFrame(index=pd.Index([f"{cohort}::{value}" for value in source_obs_name], name="cell_id"))
    if not out.index.is_unique:
        raise RuntimeError(f"Deterministic atlas cell IDs are duplicated in {cohort}")
    out["source_obs_name"] = source_obs_name.to_numpy()
    out["original_cell_id"] = original.to_numpy()
    out["source_gse_id"] = source.to_numpy()
    out["input_cohort_id"] = cohort
    out["atlas_input_role"] = row["atlas_input_role"]
    out["model_evaluation_role_frozen"] = row["model_evaluation_role_frozen"]

    text_specs = {
        "source_accession": ["source_accession", "source_gse_id"],
        "project_name": ["project_name", "project name", "source_name", "reference"],
        "sample_id": ["sample_id", "biological_sample", "sampleid", "sample", "Sample Name", "orig.ident", "cell_to_sample_ID"],
        "library_id": ["library_id", "library", "Sample Name", "sample", "sample_id", "batch"],
        "donor_id": ["donor_id", "donor", "patient_id", "patient_ID", "patient"],
        "barcode": ["barcode", "barcodes", "cell_barcode", "barcode_original", "raw_barcode"],
        "barcode_core": ["barcode_core", "tcr_barcode"],
        "tissue_original": ["tissue_original", "tissue", "Tissue", "tissue_region", "sample_tissue"],
        "tissue_harmonized": ["tissue_harmonized", "tissue_corrected", "tissue", "sample_tissue"],
        "specimen_context": ["specimen_context", "sample_tissue_status", "diagnosis_tissue", "type"],
        "diagnosis": ["diagnosis", "disease", "condition", "colon_status", "group"],
        "condition": ["condition", "diagnosis", "disease", "colon_status", "group"],
        "treatment": ["treatment", "timepoint_group", "efficacy"],
        "technology_simple": ["technology_simple", "library_chemistry", "assay_type"],
        "assay_type": ["assay_type", "technology_simple"],
        "input_population": ["input_population", "subset"],
        "source_cell_type": ["cell_type_level2", "cell_type", "major_cluster", "phase1_annotation_label", "cell_marker", "cluster_ID"],
        "source_cell_type_level1": ["cell_type_level1", "major_cluster", "phase1_annotation_label"],
        "sorted_gdt": ["Sorted_gdT"],
        "tcr_chain_mode": ["tcr_chain_mode", "tcr_status"],
        "tcr_schema_provenance": ["tcr_schema_provenance", "TCR_TRA_source"],
    }
    for target, aliases in text_specs.items():
        out[target] = first_text(obs, aliases).to_numpy()
    if cohort.startswith("GDT") or cohort == "MalteGDT":
        out["sorted_gdt"] = "True"
        out["input_population"] = "purified_gdt"

    chain_aliases = {
        "cdr3": ["{c}_cdr3", "{c}_cdr3_aa"],
        "cdr3_nt": ["{c}_cdr3_nt"],
        "v": ["{c}_v", "{c}_v_gene"],
        "d": ["{c}_d", "{c}_d_gene"],
        "j": ["{c}_j", "{c}_j_gene"],
        "clone_id": ["{c}_clone_id", "tcr_clone_id"],
        "umis": ["{c}_umis", "TCR_{c}_umi_count"],
        "reads": ["{c}_reads", "TCR_{c}_read_count"],
        "c_gene": ["{c}_c_gene"],
    }
    detected: dict[str, pd.Series] = {}
    for chain in ("TRA", "TRB", "TRG", "TRD"):
        for field, aliases in chain_aliases.items():
            out[f"{chain}_{field}"] = first_text(obs, [name.format(c=chain) for name in aliases]).to_numpy()
        detected[chain] = (
            normalized_text(out[f"{chain}_cdr3"], out.index).ne("")
            | bool_values(obs, [f"has_{chain}", f"productive_{chain}", f"TCR_{chain}_productive"]).to_numpy()
        )
        out[f"has_{chain}"] = detected[chain].to_numpy(dtype=bool)
    out["has_TRA_TRB_paired"] = (detected["TRA"] & detected["TRB"]).to_numpy(dtype=bool)
    out["has_TRG_TRD_paired"] = (detected["TRG"] & detected["TRD"]).to_numpy(dtype=bool)
    out["has_any_ab_tcr"] = (detected["TRA"] | detected["TRB"]).to_numpy(dtype=bool)
    out["has_any_gd_tcr"] = (detected["TRG"] | detected["TRD"]).to_numpy(dtype=bool)
    out["TCRseq"] = np.where(out["has_any_ab_tcr"] | out["has_any_gd_tcr"], "yes", "no")

    numeric_specs = {
        "total_counts": ["total_counts", "n_counts", "published_n_counts", "nCount_RNA"],
        "n_genes_by_counts": ["n_genes_by_counts", "n_genes", "published_n_genes", "nFeature_RNA"],
        "pct_counts_mt": ["pct_counts_mt", "pct_mito"],
        "pct_counts_ribo": ["pct_counts_ribo"],
    }
    for target, aliases in numeric_specs.items():
        value = pd.Series(np.nan, index=obs.index, dtype=np.float32)
        for alias in aliases:
            if alias in obs:
                candidate = pd.to_numeric(obs[alias], errors="coerce").astype(np.float32)
                value = value.fillna(candidate)
        out[target] = value.to_numpy()

    source_values_out = normalized_text(out["source_gse_id"], out.index)
    sample_values = normalized_text(out["sample_id"], out.index)
    library_values = normalized_text(out["library_id"], out.index)
    donor_values = normalized_text(out["donor_id"], out.index)
    batch = source_values_out.copy()
    level = pd.Series("source_gse_id", index=out.index, dtype="string")
    for values, name in [(donor_values, "donor_id"), (sample_values, "sample_id"), (library_values, "library_id")]:
        take = values.ne("")
        batch.loc[take] = source_values_out.loc[take] + "::" + name + "::" + values.loc[take]
        level.loc[take] = name
    out["integration_batch"] = batch.to_numpy()
    out["integration_batch_level"] = level.to_numpy()
    out["input_row"] = np.arange(len(out), dtype=np.int64)

    for column in out.columns:
        if column in HIGH_CARDINALITY_TEXT or pd.api.types.is_bool_dtype(out[column]) or pd.api.types.is_numeric_dtype(out[column]):
            continue
        text_values = normalized_text(out[column], out.index).fillna("").astype(str)
        out[column] = pd.Categorical(text_values)
    return out


def choose_gene_names(adata: ad.AnnData, target: pd.Index) -> None:
    best = pd.Index(adata.var_names.astype(str))
    best_overlap = int(best.isin(target).sum())
    for column in ("gene_symbol", "feature_name"):
        if column not in adata.var:
            continue
        candidate = pd.Index(adata.var[column].astype(str))
        overlap = int(candidate.isin(target).sum())
        if overlap > best_overlap:
            best, best_overlap = candidate, overlap
    adata.var_names = best.astype(str)
    if adata.var_names.has_duplicates:
        adata.var_names_make_unique()


def select_count_matrix(adata: ad.AnnData, row: dict[str, Any]) -> tuple[sparse.csr_matrix, dict[str, Any]]:
    mode = row["matrix_mode"]
    if mode == "layer_counts":
        matrix = adata.layers["counts"]
    else:
        matrix = adata.X
    matrix = matrix.tocsr().astype(np.float32) if sparse.issparse(matrix) else sparse.csr_matrix(np.asarray(matrix), dtype=np.float32)
    detail: dict[str, Any] = {"matrix_mode": mode}
    if mode == "seurat_lognormalize_inverse":
        count_column = next((name for name in ("nCount_RNA", "total_counts") if name in adata.obs), None)
        if count_column is None:
            raise KeyError("HRA005041 inverse transformation requires nCount_RNA or total_counts")
        probe = matrix[: min(1000, matrix.shape[0])].copy()
        probe.data = np.expm1(probe.data)
        probe_sum = np.asarray(probe.sum(axis=1)).ravel()
        median = float(np.median(probe_sum))
        if not 9500 <= median <= 10500:
            raise RuntimeError(f"Seurat LogNormalize probe median is {median:.2f}, expected approximately 10000")
        counts = pd.to_numeric(adata.obs[count_column], errors="coerce").to_numpy(dtype=np.float32)
        if not np.isfinite(counts).all():
            raise RuntimeError(f"{count_column} contains non-finite values")
        matrix = matrix.copy()
        matrix.data = np.expm1(matrix.data)
        matrix = matrix.multiply(counts[:, None] / 10_000.0).tocsr()
        matrix.data = np.rint(matrix.data).astype(np.float32, copy=False)
        matrix.eliminate_zeros()
        detail.update({"inverse_count_column": count_column, "probe_expm1_sum_median": median})
    sample = matrix.data[: min(100_000, matrix.nnz)]
    if not np.isfinite(sample).all() or (sample.size and np.mean(np.isclose(sample, np.rint(sample))) < 0.999):
        raise RuntimeError(f"Selected count matrix for {row['cohort_id']} is not finite integer-like count data")
    return matrix, detail


def align_matrix(matrix: sparse.csr_matrix, source_genes: pd.Index, target_genes: pd.Index) -> tuple[sparse.csr_matrix, int]:
    lookup = {gene: position for position, gene in enumerate(source_genes.astype(str))}
    shared_target_positions = np.asarray([i for i, gene in enumerate(target_genes.astype(str)) if gene in lookup], dtype=np.int64)
    source_positions = np.asarray([lookup[str(target_genes[i])] for i in shared_target_positions], dtype=np.int64)
    aligned = matrix[:, source_positions].tocsr()
    aligned.indices = shared_target_positions[aligned.indices]
    aligned._shape = (matrix.shape[0], len(target_genes))
    aligned.sort_indices()
    return aligned.astype(np.float32, copy=False), int(len(target_genes) - len(shared_target_positions))


def verify_preflight(paths: dict[str, Path]) -> dict[str, Any]:
    path = paths["log_dir"] / "preflight_summary.json"
    if not path.exists():
        raise RuntimeError("Run the preflight stage first")
    summary = json.loads(path.read_text(encoding="utf-8"))
    if summary.get("status") != "PASS":
        raise RuntimeError("Preflight did not pass")
    return summary


def prepared_path(paths: dict[str, Path], cohort: str) -> Path:
    return paths["prepared_dir"] / f"{cohort}.h5ad"


def prepare_one(config: dict[str, Any], row: dict[str, Any], target_genes: pd.Index, paths: dict[str, Path], overwrite: bool) -> dict[str, Any]:
    cohort = row["cohort_id"]
    destination = prepared_path(paths, cohort)
    print(f"[{time.strftime('%F %T')}] prepare input {cohort}", flush=True)
    if destination.exists() and not overwrite:
        with h5py.File(destination, "r") as handle:
            shape = matrix_shape(matrix_group(handle, "X"))
        if shape == (int(row["expected_cells_effective"]), int(config["expected_genes"])):
            return {"cohort_id": cohort, "status": "REUSED", "cells": shape[0], "genes": shape[1], "prepared_h5ad": str(destination), "prepared_sha256": sha256_file(destination)}
        raise RuntimeError(f"Existing prepared file has wrong shape: {destination}: {shape}")

    source = input_path(row)
    adata = ad.read_h5ad(source)
    source_obs = adata.obs.copy()
    if row.get("apply_core_recovery"):
        exclusion = set(pd.read_csv(resolve(config["core_recovery"]["row_exclusion_manifest"]), usecols=["obs_name"])["obs_name"].astype(str))
        keep = ~pd.Index(adata.obs_names.astype(str)).isin(exclusion)
        adata = adata[keep].copy()
        source_obs = adata.obs.copy()
        recovery_frame = source_obs.copy()
        recovery_frame["source_original_cell_id"] = first_text(recovery_frame, ["original_cell_id"], default="").mask(lambda value: value.eq(""), pd.Series(recovery_frame.index.astype(str), index=recovery_frame.index)).to_numpy()
        recovery_config = {"current_atlas_recovery": config["core_recovery"]}
        source_obs = attach_recovery_metadata(recovery_frame, recovery_config)
        adata.obs = source_obs
    if adata.n_obs != int(row["expected_cells_effective"]):
        raise RuntimeError(f"Effective cell count mismatch for {cohort}: {adata.n_obs}")
    choose_gene_names(adata, target_genes)
    matrix, detail = select_count_matrix(adata, row)
    obs = standardize_obs(adata.obs.copy(), row)
    aligned, missing = align_matrix(matrix, pd.Index(adata.var_names.astype(str)), target_genes)
    prepared = ad.AnnData(X=aligned, obs=obs, var=pd.DataFrame(index=target_genes.astype(str)))
    prepared.uns["full_atlas_input_provenance"] = {
        "protocol_version": config["protocol_version"],
        "source_h5ad": str(source),
        "source_sha256": row["expected_sha256"],
        "matrix_mode": row["matrix_mode"],
        "repaired_tcr_sidecar_propagated": False,
        **detail,
    }
    atomic_h5ad(prepared, destination, config["staging"]["prepared_compression"])
    result = {
        "cohort_id": cohort,
        "status": "BUILT",
        "cells": int(prepared.n_obs),
        "genes": int(prepared.n_vars),
        "matrix_nnz": int(prepared.X.nnz),
        "missing_target_genes": missing,
        "prepared_h5ad": str(destination),
        "prepared_sha256": sha256_file(destination),
    }
    del adata, source_obs, matrix, aligned, prepared
    gc.collect()
    print(f"[{time.strftime('%F %T')}] prepared {cohort}: {result['cells']:,} cells", flush=True)
    return result


def hvg_excluded(gene: str) -> bool:
    value = str(gene).upper()
    return bool(TCR_HVG_PATTERN.match(value) or OTHER_HVG_PATTERN.match(value))


def select_hvgs(config: dict[str, Any], prepared: list[tuple[dict[str, Any], Path]], target_genes: pd.Index, paths: dict[str, Path]) -> list[str]:
    cap = int(config["hvg"]["sample_cap_per_source"])
    seed = int(config["random_seed"])
    blocks: list[sparse.csr_matrix] = []
    sample_batches: list[np.ndarray] = []
    sample_rows: list[dict[str, Any]] = []
    for offset, (row, path) in enumerate(prepared):
        print(f"[{time.strftime('%F %T')}] HVG sample {row['cohort_id']}", flush=True)
        backed = ad.read_h5ad(path, backed="r")
        sources = backed.obs["source_gse_id"].astype(str).to_numpy()
        selected_parts: list[np.ndarray] = []
        for source_offset, source in enumerate(sorted(set(sources))):
            candidates = np.flatnonzero(sources == source)
            if candidates.size > cap:
                rng = np.random.default_rng(seed + offset * 104729 + source_offset * 1009)
                candidates = np.sort(rng.choice(candidates, size=cap, replace=False))
            selected_parts.append(candidates)
            sample_rows.append({"cohort_id": row["cohort_id"], "source_gse_id": source, "n_available": int(np.sum(sources == source)), "n_sampled": int(candidates.size)})
        selected = np.sort(np.concatenate(selected_parts)) if selected_parts else np.asarray([], dtype=np.int64)
        blocks.append(read_sparse_rows_genes(path, "X", target_genes.tolist(), rows=selected, row_chunk_size=25_000))
        sample_batches.append(sources[selected])
        backed.file.close()
    matrix = sparse.vstack(blocks, format="csr", dtype=np.float32)
    batches = np.concatenate(sample_batches)
    sample = ad.AnnData(X=matrix, obs=pd.DataFrame({"source_gse_id": pd.Categorical(batches)}), var=pd.DataFrame(index=target_genes))
    print(f"[{time.strftime('%F %T')}] selecting HVGs on {sample.n_obs:,} source-balanced cells", flush=True)
    candidate_n = min(sample.n_vars, max(int(config["hvg"]["n_top_genes"]) * 3, 12_000))
    sc.pp.highly_variable_genes(sample, flavor=config["hvg"]["flavor"], n_top_genes=candidate_n, batch_key=config["hvg"]["batch_key"], subset=False)
    rank = pd.to_numeric(sample.var.get("highly_variable_rank"), errors="coerce")
    eligible = pd.Series([not hvg_excluded(gene) for gene in sample.var_names], index=sample.var_names)
    ranked = rank[eligible & rank.notna()].sort_values().index.astype(str).tolist()
    n_top = int(config["hvg"]["n_top_genes"])
    if len(ranked) < n_top:
        raise RuntimeError(f"Only {len(ranked)} eligible ranked HVGs are available")
    forced = [gene for gene in config["hvg"]["forced_context_genes"] if gene in target_genes and not hvg_excluded(gene)]
    selected = ranked[:n_top]
    forced_set = set(forced)
    for gene in forced:
        if gene not in selected:
            replacement = next(
                (position for position in range(len(selected) - 1, -1, -1) if selected[position] not in forced_set),
                None,
            )
            if replacement is None:
                raise RuntimeError("No replaceable HVG remains while adding forced context genes")
            selected[replacement] = gene
    selected = list(dict.fromkeys(selected))
    for gene in ranked:
        if len(selected) >= n_top:
            break
        if gene not in selected:
            selected.append(gene)
    hvg_table = pd.DataFrame({"gene": sample.var_names.astype(str), "highly_variable_rank": rank.to_numpy(), "excluded": ~eligible.to_numpy(), "selected": sample.var_names.astype(str).isin(selected), "forced_context": sample.var_names.astype(str).isin(forced)})
    hvg_table.to_csv(paths["table_dir"] / "hvg_selection.csv", index=False)
    pd.DataFrame(sample_rows).to_csv(paths["table_dir"] / "hvg_source_balanced_sample.csv", index=False)
    (paths["table_dir"] / "hvg_genes.txt").write_text("\n".join(selected) + "\n", encoding="utf-8")
    del sample, matrix, blocks
    gc.collect()
    return selected


def build_hvg_h5ad(config: dict[str, Any], prepared: list[tuple[dict[str, Any], Path]], genes: list[str], paths: dict[str, Path], overwrite: bool) -> None:
    destination = paths["hvg_h5ad"]
    if destination.exists() and not overwrite:
        with h5py.File(destination, "r") as handle:
            shape = matrix_shape(matrix_group(handle, "X"))
        if shape == (int(config["expected_total_cells"]), len(genes)):
            return
        raise RuntimeError(f"Existing HVG H5AD has wrong shape: {shape}")
    matrices: list[sparse.csr_matrix] = []
    metadata: list[pd.DataFrame] = []
    keep_obs = ["source_obs_name", "original_cell_id", "source_gse_id", "input_cohort_id", "atlas_input_role", "model_evaluation_role_frozen", "sample_id", "library_id", "donor_id", "integration_batch", "integration_batch_level", "input_row"]
    for row, path in prepared:
        print(f"[{time.strftime('%F %T')}] stage 4,000 HVGs from {row['cohort_id']}", flush=True)
        matrices.append(read_sparse_rows_genes(path, "X", genes, rows=None, row_chunk_size=25_000))
        backed = ad.read_h5ad(path, backed="r")
        frame = backed.obs[keep_obs].copy()
        backed.file.close()
        metadata.append(frame)
    matrix = sparse.vstack(matrices, format="csr", dtype=np.float32)
    obs = pd.concat(metadata, axis=0, copy=False)
    if matrix.shape != (int(config["expected_total_cells"]), len(genes)) or len(obs) != matrix.shape[0] or not obs.index.is_unique:
        raise RuntimeError("HVG matrix/metadata dimensions or cell IDs violate the frozen contract")
    adata = ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    adata.uns["full_atlas_rebuild"] = {"protocol_version": config["protocol_version"], "raw_counts": True, "tcr_sidecar_propagated": False}
    atomic_h5ad(adata, destination, config["staging"]["prepared_compression"])
    print(f"[{time.strftime('%F %T')}] wrote HVG H5AD: {destination}", flush=True)
    del adata, matrix, matrices, obs, metadata
    gc.collect()


def prepare_stage(config_path: Path, config: dict[str, Any], paths: dict[str, Path], overwrite: bool) -> dict[str, Any]:
    verify_preflight(paths)
    target_path = input_path(config["inputs"][0])
    with h5py.File(target_path, "r") as handle:
        target_genes = pd.Index(axis_names(handle, "var").astype(str), name="gene")
    before = pd.read_csv(paths["table_dir"] / "source_file_state_before.csv")
    results = []
    prepared = []
    for row in config["inputs"]:
        result = prepare_one(config, row, target_genes, paths, overwrite)
        results.append(result)
        prepared.append((row, prepared_path(paths, row["cohort_id"])))
    result_df = pd.DataFrame(results)
    result_df.to_csv(paths["table_dir"] / "prepared_inputs.csv", index=False)
    genes = select_hvgs(config, prepared, target_genes, paths)
    build_hvg_h5ad(config, prepared, genes, paths, overwrite)
    after_rows = []
    for row in config["inputs"]:
        path = input_path(row)
        stat = path.stat()
        after_rows.append({"cohort_id": row["cohort_id"], "path": str(path), "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns})
    after = pd.DataFrame(after_rows)
    unchanged = before.merge(after, on=["cohort_id", "path"], suffixes=("_before", "_after"), validate="one_to_one")
    unchanged["unchanged"] = unchanged["size_bytes_before"].eq(unchanged["size_bytes_after"]) & unchanged["mtime_ns_before"].eq(unchanged["mtime_ns_after"])
    unchanged.to_csv(paths["table_dir"] / "source_file_state_after_prepare.csv", index=False)
    if not unchanged["unchanged"].all():
        raise RuntimeError("At least one source H5AD changed during prepare")
    with h5py.File(paths["hvg_h5ad"], "r") as handle:
        hvg_shape = matrix_shape(matrix_group(handle, "X"))
    summary = {
        "stage": "prepare",
        "status": "PASS",
        "n_cells": hvg_shape[0],
        "n_hvgs": hvg_shape[1],
        "n_prepared_inputs": len(prepared),
        "hvg_h5ad": str(paths["hvg_h5ad"]),
        "hvg_h5ad_sha256": sha256_file(paths["hvg_h5ad"]),
        "source_h5ads_unchanged": True,
        "tcr_sidecar_propagated": False,
        "config_sha256": sha256_file(config_path),
    }
    atomic_json(paths["log_dir"] / "prepare_summary.json", summary)
    return summary


def require_gpu() -> dict[str, Any]:
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", "0")
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    import torch

    if not torch.cuda.is_available():
        raise RuntimeError("A CUDA GPU is required; CPU fallback is forbidden")
    return {"name": torch.cuda.get_device_name(0), "memory_bytes": int(torch.cuda.get_device_properties(0).total_memory)}


def fit_stage(config: dict[str, Any], paths: dict[str, Path], overwrite: bool) -> dict[str, Any]:
    if not (paths["log_dir"] / "prepare_summary.json").exists():
        raise RuntimeError("Prepare stage has not passed")
    if paths["latent_npy"].exists() and not overwrite:
        raise FileExistsError(paths["latent_npy"])
    gpu = require_gpu()
    import scvi

    adata = ad.read_h5ad(paths["hvg_h5ad"])
    if adata.shape != (int(config["expected_total_cells"]), int(config["hvg"]["n_top_genes"])):
        raise RuntimeError(f"Unexpected HVG H5AD shape: {adata.shape}")
    scvi.settings.seed = int(config["random_seed"])
    scvi.model.SCVI.setup_anndata(adata, batch_key="integration_batch")
    model = scvi.model.SCVI(adata, n_hidden=int(config["scvi"]["n_hidden"]), n_layers=int(config["scvi"]["n_layers"]), n_latent=int(config["scvi"]["n_latent"]), gene_likelihood=config["scvi"]["gene_likelihood"])
    started = time.time()
    model.train(max_epochs=int(config["scvi"]["max_epochs"]), batch_size=int(config["scvi"]["batch_size"]), accelerator="gpu", devices=1, early_stopping=bool(config["scvi"]["early_stopping"]), early_stopping_patience=int(config["scvi"]["early_stopping_patience"]))
    latent = model.get_latent_representation(batch_size=int(config["scvi"]["batch_size"])).astype(np.float32)
    if latent.shape != (adata.n_obs, int(config["scvi"]["n_latent"])) or not np.isfinite(latent).all():
        raise RuntimeError("scVI latent matrix failed shape/finite validation")
    temporary = paths["latent_npy"].with_name(paths["latent_npy"].name + ".partial")
    with temporary.open("wb") as handle:
        np.save(handle, latent)
    temporary.replace(paths["latent_npy"])
    if paths["scvi_model"].exists() and overwrite:
        shutil.rmtree(paths["scvi_model"])
    model.save(paths["scvi_model"], overwrite=overwrite, save_anndata=False)
    history = pd.DataFrame({key: np.ravel(value) for key, value in model.history.items() if np.asarray(value).ndim <= 2})
    history.to_csv(paths["table_dir"] / "scvi_training_history.csv", index=False)
    summary = {"stage": "fit", "status": "PASS", "elapsed_seconds": time.time() - started, "n_cells": int(adata.n_obs), "n_latent": int(latent.shape[1]), "latent_sha256": sha256_file(paths["latent_npy"]), "gpu": gpu, "cpu_fallback": False}
    atomic_json(paths["log_dir"] / "fit_summary.json", summary)
    return summary


def embed_stage(config: dict[str, Any], paths: dict[str, Path], overwrite: bool) -> dict[str, Any]:
    if paths["umap_npy"].exists() and not overwrite:
        raise FileExistsError(paths["umap_npy"])
    gpu = require_gpu()
    import rapids_singlecell as rsc

    latent = np.load(paths["latent_npy"], mmap_mode="r")
    if latent.shape != (int(config["expected_total_cells"]), int(config["scvi"]["n_latent"])):
        raise RuntimeError(f"Unexpected latent shape: {latent.shape}")
    adata = ad.AnnData(X=np.asarray(latent))
    started = time.time()
    rsc.get.anndata_to_GPU(adata)
    rsc.pp.neighbors(adata, n_neighbors=int(config["embedding"]["n_neighbors"]), metric=config["embedding"]["metric"], use_rep="X")
    rsc.tl.leiden(adata, resolution=float(config["embedding"]["leiden_resolution"]), random_state=int(config["random_seed"]))
    rsc.tl.umap(adata, random_state=int(config["random_seed"]))
    rsc.get.anndata_to_CPU(adata)
    umap = np.asarray(adata.obsm["X_umap"], dtype=np.float32)
    leiden = adata.obs["leiden"].astype(str).to_numpy(dtype="U32")
    for target, values in [(paths["umap_npy"], umap), (paths["leiden_npy"], leiden)]:
        temporary = target.with_name(target.name + ".partial")
        with temporary.open("wb") as handle:
            np.save(handle, values)
        temporary.replace(target)
    pd.Series(leiden).value_counts().rename_axis("leiden").reset_index(name="n_cells").to_csv(paths["table_dir"] / "leiden_counts.csv", index=False)
    summary = {"stage": "embed", "status": "PASS", "elapsed_seconds": time.time() - started, "n_cells": int(len(leiden)), "n_clusters": int(pd.Series(leiden).nunique()), "umap_sha256": sha256_file(paths["umap_npy"]), "leiden_sha256": sha256_file(paths["leiden_npy"]), "gpu": gpu}
    atomic_json(paths["log_dir"] / "embed_summary.json", summary)
    return summary


def assemble_stage(config_path: Path, config: dict[str, Any], paths: dict[str, Path], overwrite: bool) -> dict[str, Any]:
    destination = paths["final_h5ad"]
    if destination.exists() and not overwrite:
        raise FileExistsError(destination)
    latent = np.load(paths["latent_npy"], mmap_mode="r")
    umap = np.load(paths["umap_npy"], mmap_mode="r")
    leiden = np.load(paths["leiden_npy"], mmap_mode="r")
    inputs = []
    for row in config["inputs"]:
        print(f"[{time.strftime('%F %T')}] load prepared full-gene input {row['cohort_id']}", flush=True)
        inputs.append(ad.read_h5ad(prepared_path(paths, row["cohort_id"])))
    started = time.time()
    merged = ad.concat(inputs, axis=0, join="inner", merge="same", index_unique=None)
    print(f"[{time.strftime('%F %T')}] assembled sparse matrix {merged.n_obs:,} x {merged.n_vars:,}", flush=True)
    del inputs
    gc.collect()
    merged.X = merged.X.tocsr().astype(np.float32, copy=False)
    if merged.shape != (int(config["expected_total_cells"]), int(config["expected_genes"])) or not merged.obs_names.is_unique:
        raise RuntimeError(f"Assembled atlas violates shape/ID contract: {merged.shape}, unique={merged.obs_names.is_unique}")
    if latent.shape[0] != merged.n_obs or umap.shape != (merged.n_obs, 2) or leiden.shape[0] != merged.n_obs:
        raise RuntimeError("Embedding arrays do not align to assembled cells")
    merged.obsm["X_scVI"] = np.asarray(latent, dtype=np.float32)
    merged.obsm["X_umap"] = np.asarray(umap, dtype=np.float32)
    merged.obs["leiden"] = pd.Categorical(np.asarray(leiden).astype(str))
    merged.uns["full_atlas_rebuild"] = {
        "protocol_version": config["protocol_version"],
        "config_sha256": sha256_file(config_path),
        "script_sha256": sha256_file(SCRIPT_PATH),
        "raw_counts_in_X": True,
        "n_input_cohorts": len(config["inputs"]),
        "historical_cells": int(config["expected_historical_cells"]),
        "extension_cells": int(config["expected_extension_cells"]),
        "copd_cells": int(config["expected_copd_cells"]),
        "repaired_tcr_sidecar_propagated": False,
        "source_policy": config["source_policy"],
    }
    temporary = destination.with_name(destination.name + ".partial")
    temporary.unlink(missing_ok=True)
    print(f"[{time.strftime('%F %T')}] writing final compressed H5AD", flush=True)
    merged.write_h5ad(temporary, compression=config["staging"]["final_compression"])
    with h5py.File(temporary, "r") as handle:
        observed = matrix_shape(matrix_group(handle, "X"))
        has_latent = "X_scVI" in handle["obsm"]
        has_umap = "X_umap" in handle["obsm"]
    if observed != merged.shape or not has_latent or not has_umap:
        raise RuntimeError("Written atlas failed post-write shape/embedding checks")
    temporary.replace(destination)
    link = paths["nfs_link"]
    link.parent.mkdir(parents=True, exist_ok=True)
    if link.exists() and not link.is_symlink():
        raise RuntimeError(f"Refusing to replace a non-symlink canonical path: {link}")
    link_tmp = link.with_name(link.name + ".partial")
    link_tmp.unlink(missing_ok=True)
    link_tmp.symlink_to(destination)
    link_tmp.replace(link)
    counts = merged.obs.groupby(["input_cohort_id", "source_gse_id"], observed=True).size().rename("n_cells").reset_index()
    counts.to_csv(paths["table_dir"] / "atlas_cells_by_input_and_source.csv", index=False)
    summary = {"stage": "assemble", "status": "PASS", "elapsed_seconds": time.time() - started, "n_cells": int(merged.n_obs), "n_genes": int(merged.n_vars), "n_source_accessions": int(merged.obs["source_gse_id"].nunique()), "n_input_cohorts": int(merged.obs["input_cohort_id"].nunique()), "matrix_nnz": int(merged.X.nnz), "final_h5ad": str(destination), "nfs_link": str(link), "final_size_bytes": destination.stat().st_size, "tcr_sidecar_propagated": False}
    atomic_json(paths["log_dir"] / "assemble_summary.json", summary)
    return summary


def score_stage(config: dict[str, Any], paths: dict[str, Path]) -> dict[str, Any]:
    target = paths["final_h5ad"]
    with h5py.File(target, "r") as handle:
        var_names = pd.Index(axis_names(handle, "var").astype(str), dtype="string")
        n_obs = len(axis_names(handle, "obs"))
        n_vars = len(var_names)
        if all(column in handle["obs"] for column in PHASE4_SCORE_COLUMNS.values()):
            return {"stage": "score", "status": "REUSED", "n_cells": n_obs}
        leiden = np.asarray(obs_values(handle, "leiden"), dtype=object).astype(str)
    modules = find_module_genes(var_names)
    gene_means = compute_gene_means(target, n_obs, n_vars, CHUNK_SIZE)
    controls = {name: pick_control_genes(genes, var_names, gene_means, random_state=RANDOM_STATE) for name, genes in modules.items()}
    module_idx = {name: var_names.get_indexer(genes).astype(np.int32) for name, genes in modules.items()}
    control_idx = {name: var_names.get_indexer(genes).astype(np.int32) for name, genes in controls.items()}
    codes, _ = pd.factorize(leiden, sort=True)
    started = time.time()
    scores, _, _ = compute_scores(target, n_obs=n_obs, n_vars=n_vars, chunk_size=CHUNK_SIZE, module_gene_idx=module_idx, module_ctrl_idx=control_idx, leiden_codes=codes, marker_idx=np.asarray([], dtype=np.int32))
    scores, scaling = add_scaled_scores(scores)
    columns = {column: scores[name] for name, column in {**PHASE4_SCORE_COLUMNS, **PHASE4_SCALED_SCORE_COLUMNS}.items()}
    payload = {
        "package_source": str(PACKAGE_ZIP),
        "integrated_h5ad": str(target),
        "scoring_mode": "temporary_normalize_total_log1p_on_raw_count_X",
        "continuous_only": True,
        "random_state": RANDOM_STATE,
        "target_sum": TARGET_SUM,
        "ctrl_size": CTRL_SIZE,
        "n_bins": N_BINS,
        "module_genes": {name: [str(gene) for gene in genes] for name, genes in modules.items()},
        "control_genes": {name: [str(gene) for gene in genes] for name, genes in controls.items()},
        "score_columns": {**PHASE4_SCORE_COLUMNS, **PHASE4_SCALED_SCORE_COLUMNS},
        "scaled_score_ranges": scaling,
        "full_atlas_rebuild": True,
    }
    append_obs_columns_in_place(target, columns, payload)
    summary_rows = []
    for name, column in {**PHASE4_SCORE_COLUMNS, **PHASE4_SCALED_SCORE_COLUMNS}.items():
        value = scores[name]
        summary_rows.append({"score": column, "mean": float(np.mean(value)), "sd": float(np.std(value)), "median": float(np.median(value)), "p01": float(np.quantile(value, 0.01)), "p99": float(np.quantile(value, 0.99))})
    pd.DataFrame(summary_rows).to_csv(paths["table_dir"] / "phase4_score_summary.csv", index=False)
    summary = {"stage": "score", "status": "PASS", "elapsed_seconds": time.time() - started, "n_cells": n_obs, "score_columns": list(columns), "phase4_module": {name: list(map(str, genes)) for name, genes in modules.items()}}
    atomic_json(paths["log_dir"] / "score_summary.json", summary)
    return summary


def stratified_plot_sample(obs: pd.DataFrame, cap: int, seed: int) -> np.ndarray:
    if len(obs) <= cap:
        return np.arange(len(obs), dtype=np.int64)
    rng = np.random.default_rng(seed)
    selected = []
    for _, group in obs.groupby("source_gse_id", observed=True, sort=True):
        target = max(500, int(round(cap * len(group) / len(obs))))
        positions = obs.index.get_indexer(group.index)
        selected.append(rng.choice(positions, size=min(target, len(positions)), replace=False))
    result = np.unique(np.concatenate(selected))
    if len(result) > cap:
        result = np.sort(rng.choice(result, size=cap, replace=False))
    return result


def plot_embedding(umap: np.ndarray, labels: np.ndarray, title: str, path: Path, seed: int, cap: int = 300_000) -> None:
    frame = pd.DataFrame({"label": pd.Categorical(labels)})
    selected = stratified_plot_sample(frame.rename(columns={"label": "source_gse_id"}), cap, seed)
    values = frame["label"].astype(str).to_numpy()[selected]
    categories = sorted(set(values))
    cmap = plt.get_cmap("tab20")
    fig, ax = plt.subplots(figsize=(13, 9), constrained_layout=True)
    for i, category in enumerate(categories):
        take = values == category
        ax.scatter(umap[selected[take], 0], umap[selected[take], 1], s=0.5, alpha=0.55, linewidths=0, rasterized=True, color=cmap(i % 20), label=category)
    ax.set_title(title)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(markerscale=8, fontsize=7, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left", ncol=2 if len(categories) > 24 else 1)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def report_stage(config: dict[str, Any], paths: dict[str, Path], send_email: bool) -> dict[str, Any]:
    target = paths["final_h5ad"]
    adata = ad.read_h5ad(target, backed="r")
    obs = adata.obs.copy()
    umap = np.asarray(adata.obsm["X_umap"])
    plot_embedding(umap, obs["source_gse_id"].astype(str).to_numpy(), "Full T/NK atlas by source", paths["figure_dir"] / "umap_by_source.png", int(config["random_seed"]))
    plot_embedding(umap, obs["leiden"].astype(str).to_numpy(), "Full T/NK atlas by unsupervised cluster", paths["figure_dir"] / "umap_by_leiden.png", int(config["random_seed"]) + 1)
    adata.file.close()
    by_input = obs.groupby(["input_cohort_id", "atlas_input_role", "model_evaluation_role_frozen"], observed=True).size().rename("n_cells").reset_index().sort_values("n_cells", ascending=False)
    by_source = obs.groupby("source_gse_id", observed=True).size().rename("n_cells").reset_index().sort_values("n_cells", ascending=False)
    by_tissue = obs.groupby("tissue_harmonized", observed=True).size().rename("n_cells").reset_index().sort_values("n_cells", ascending=False)
    by_input.to_csv(paths["table_dir"] / "report_cells_by_input.csv", index=False)
    by_source.to_csv(paths["table_dir"] / "report_cells_by_source.csv", index=False)
    by_tissue.to_csv(paths["table_dir"] / "report_cells_by_tissue.csv", index=False)
    assemble = json.loads((paths["log_dir"] / "assemble_summary.json").read_text())
    fit = json.loads((paths["log_dir"] / "fit_summary.json").read_text())
    embed = json.loads((paths["log_dir"] / "embed_summary.json").read_text())
    image_source = os.path.relpath(paths["figure_dir"] / "umap_by_source.png", paths["report_dir"])
    image_cluster = os.path.relpath(paths["figure_dir"] / "umap_by_leiden.png", paths["report_dir"])
    css = """
    body{font-family:Arial,sans-serif;max-width:1180px;margin:32px auto;color:#17202a;line-height:1.45;padding:0 20px}
    h1,h2{color:#123b4a} .metrics{display:grid;grid-template-columns:repeat(4,1fr);gap:10px}.metric{border-top:4px solid #2a788e;background:#f4f7f8;padding:14px}.metric b{font-size:24px;display:block}
    figure{margin:24px 0}img{max-width:100%;height:auto}table{border-collapse:collapse;width:100%;font-size:12px}th,td{border-bottom:1px solid #d8dee2;padding:6px;text-align:left}th{background:#eaf1f3;position:sticky;top:0}.note{background:#fff4d6;border-left:5px solid #e2a400;padding:12px}.page-break{break-before:page}@media print{body{max-width:none}.metrics{grid-template-columns:repeat(4,1fr)}table{font-size:9px}thead{display:table-header-group}tr{break-inside:avoid}img{max-height:170mm;object-fit:contain}}
    """
    report = f"""<!doctype html><html><head><meta charset='utf-8'><title>Full T/NK Atlas Rebuild</title><style>{css}</style></head><body>
    <h1>Full T/NK Atlas Rebuild</h1><p>Rebuilt from the frozen historical atlas inputs, eight independently filtered extension cohorts, and the COPD BALF/BLOOD cohort.</p>
    <div class='metrics'><div class='metric'><b>{assemble['n_cells']:,}</b>cells</div><div class='metric'><b>{assemble['n_genes']:,}</b>genes</div><div class='metric'><b>{assemble['n_input_cohorts']}</b>input cohorts</div><div class='metric'><b>{assemble['n_source_accessions']}</b>source accessions</div></div>
    <h2>Scope and provenance</h2><p>The historical component contains {config['expected_historical_cells']:,} cells. The eight new extension cohorts add {config['expected_extension_cells']:,} cells, and COPD BALF/BLOOD adds {config['expected_copd_cells']:,} cells. Frozen model-validation roles are retained as metadata and were not used as biological labels during integration.</p>
    <div class='note'><b>TCR metadata boundary:</b> the separately validated 14-source repaired TCR sidecar was not propagated. This rebuild preserves the source H5AD TCR fields and records that limitation explicitly.</div>
    <h2>Method</h2><p>All inputs were converted to sparse raw-count space, aligned to the historical 27,413-gene universe, and assigned deterministic cohort-prefixed cell IDs. The HRA005041 Seurat LogNormalize matrix was inverted using its per-cell count totals, matching the historical workflow. Four thousand source-balanced HVGs were selected after excluding TCR V/J/D, mitochondrial, ribosomal, immunoglobulin, and common noncoding features; a fixed T/NK context panel was retained. scVI used a 30-dimensional latent space on an A100 GPU, followed by RAPIDS 15-neighbor graph construction, Leiden clustering at resolution 1.0, and UMAP.</p>
    <p>scVI runtime: {fit['elapsed_seconds']/60:.1f} minutes. Embedding runtime: {embed['elapsed_seconds']/60:.1f} minutes. CPU fallback: {fit['cpu_fallback']}.</p>
    <h2>Integrated structure</h2><figure><img src='{html.escape(image_source)}'><figcaption>UMAP colored by source accession; a stratified sample is rendered for legibility.</figcaption></figure>
    <figure><img src='{html.escape(image_cluster)}'><figcaption>UMAP colored by unsupervised Leiden cluster.</figcaption></figure>
    <h2 class='page-break'>Input composition</h2>{by_input.to_html(index=False, border=0)}
    <h2>Source composition</h2>{by_source.to_html(index=False, border=0)}
    <h2>Tissue composition</h2>{by_tissue.head(80).to_html(index=False, border=0)}
    <h2>Canonical outputs</h2><p>H5AD: <code>{html.escape(str(target))}</code><br>NFS link: <code>{html.escape(str(paths['nfs_link']))}</code><br>Model: <code>{html.escape(str(paths['scvi_model']))}</code></p>
    </body></html>"""
    index = paths["report_dir"] / "index.html"
    index.write_text(report, encoding="utf-8")
    pdf = paths["report_dir"] / "full_atlas_rebuild_report.pdf"
    subprocess.run(["google-chrome", "--headless", "--no-sandbox", "--disable-gpu", f"--print-to-pdf={pdf}", str(index)], check=True)
    summary = {"stage": "report", "status": "PASS", "html": str(index), "pdf": str(pdf), "pdf_size_bytes": pdf.stat().st_size, "email_sent": False}
    if send_email:
        sender = Path.home() / ".codex/skills/email-to-likai/scripts/send_email_to_likai.py"
        subprocess.run([sys.executable, str(sender), "--subject", "Full T/NK atlas rebuild completed", "--body", f"The full atlas rebuild passed QC with {assemble['n_cells']:,} cells. The historical atlas, eight extension cohorts, and COPD BALF/BLOOD are included. The repaired TCR sidecar remains unpropagated. Canonical H5AD: {target}", "--attachment", str(pdf)], check=True)
        summary["email_sent"] = True
    atomic_json(paths["log_dir"] / "report_summary.json", summary)
    return summary


def run_stage(args: argparse.Namespace) -> None:
    config_path = resolve(args.config)
    config = read_config(config_path)
    paths = output_paths(config)
    stages = ["preflight", "prepare", "fit", "embed", "assemble", "score", "report"] if args.stage == "all" else [args.stage]
    for stage in stages:
        started = time.time()
        print(f"[{time.strftime('%F %T')}] START {stage}", flush=True)
        if stage == "preflight":
            result = preflight_stage(config_path, config, paths)
        elif stage == "prepare":
            result = prepare_stage(config_path, config, paths, args.overwrite)
        elif stage == "fit":
            result = fit_stage(config, paths, args.overwrite)
        elif stage == "embed":
            result = embed_stage(config, paths, args.overwrite)
        elif stage == "assemble":
            result = assemble_stage(config_path, config, paths, args.overwrite)
        elif stage == "score":
            result = score_stage(config, paths)
        else:
            result = report_stage(config, paths, args.send_email)
        print(json.dumps(result, indent=2, sort_keys=True), flush=True)
        print(f"[{time.strftime('%F %T')}] END {stage} elapsed={time.time()-started:.1f}s", flush=True)


if __name__ == "__main__":
    run_stage(parse_args())
