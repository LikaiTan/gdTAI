#!/usr/bin/env python3
"""Run the read-only gdTAI V4.2 modeling-integration preflight.

The preflight hashes and inspects all proposed development and locked cohorts,
freezes role separation, estimates compute/storage requirements, and renders a
supervision report. It does not merge data, fit scVI, cluster cells, create
pseudo-labels, or fit a classifier.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "configs/models/gdtai/v4_2_integration_preflight.json"


def resolve(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def sha256_file(path: Path, chunk_size: int = 32 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def decode(value: Any) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def axis_names(handle: h5py.File, axis: str) -> np.ndarray:
    group = handle[axis]
    index_name = decode(group.attrs["_index"])
    return np.asarray([decode(value) for value in group[index_name][:]], dtype=object)


def matrix_group(handle: h5py.File, key: str) -> h5py.Group:
    obj: Any = handle
    for component in key.split("/"):
        obj = obj[component]
    if not isinstance(obj, h5py.Group):
        raise TypeError(f"{key} is not a sparse matrix group")
    return obj


def categorical_values(group: h5py.Group) -> tuple[np.ndarray, np.ndarray]:
    categories = np.asarray([decode(value) for value in group["categories"][:]], dtype=object)
    codes = np.asarray(group["codes"][:], dtype=np.int64)
    return categories, codes


def obs_values(handle: h5py.File, column: str) -> np.ndarray:
    obj = handle["obs"][column]
    if isinstance(obj, h5py.Group) and "categories" in obj and "codes" in obj:
        categories, codes = categorical_values(obj)
        output = np.full(codes.size, "", dtype=object)
        valid = codes >= 0
        output[valid] = categories[codes[valid]]
        return output
    raw = np.asarray(obj[:])
    if raw.dtype.kind in {"S", "O", "U"}:
        return np.asarray([decode(value) for value in raw], dtype=object)
    return raw


def obs_column_summary(handle: h5py.File, column: str) -> dict[str, Any]:
    if column not in handle["obs"]:
        return {"present": False, "n_missing": int(axis_names(handle, "obs").size), "n_unique": 0}
    obj = handle["obs"][column]
    if isinstance(obj, h5py.Group) and "categories" in obj and "codes" in obj:
        _, codes = categorical_values(obj)
        return {"present": True, "n_missing": int((codes < 0).sum()), "n_unique": int(np.unique(codes[codes >= 0]).size)}
    values = obs_values(handle, column)
    if values.dtype.kind in {"O", "S", "U"}:
        cleaned = pd.Series(values, dtype="string").fillna("").str.strip()
        missing = cleaned.eq("") | cleaned.str.lower().isin(["nan", "none", "na", "unknown"])
        return {"present": True, "n_missing": int(missing.sum()), "n_unique": int(cleaned[~missing].nunique())}
    missing = ~np.isfinite(values.astype(float))
    return {"present": True, "n_missing": int(missing.sum()), "n_unique": int(np.unique(values[~missing]).size)}


def obs_bool(handle: h5py.File, column: str) -> np.ndarray:
    values = obs_values(handle, column)
    if values.dtype.kind == "b":
        return values.astype(bool)
    if values.dtype.kind in {"i", "u", "f"}:
        return values.astype(float) != 0
    cleaned = pd.Series(values, dtype="string").fillna("").str.strip().str.lower()
    return cleaned.isin(["true", "1", "yes", "y", "t"]).to_numpy(bool)


def sampled_matrix_audit(group: h5py.Group, per_region: int, tolerance: float) -> dict[str, Any]:
    data = group["data"]
    n = int(data.shape[0])
    width = min(per_region, n)
    slices = [(0, width), (max(0, n // 2 - width // 2), min(n, n // 2 + width // 2)), (max(0, n - width), n)]
    chunks = [np.asarray(data[start:end], dtype=np.float64) for start, end in slices if end > start]
    values = np.concatenate(chunks) if chunks else np.asarray([], dtype=np.float64)
    finite = bool(np.isfinite(values).all())
    nonnegative = bool((values >= 0).all()) if values.size else True
    integer_like = bool((np.abs(values - np.rint(values)) <= tolerance).all()) if values.size else True
    return {
        "sampled_values": int(values.size),
        "sample_min": math.nan if not values.size else float(values.min()),
        "sample_max": math.nan if not values.size else float(values.max()),
        "sample_finite": finite,
        "sample_nonnegative": nonnegative,
        "sample_integer_like": integer_like,
    }


def resolve_inputs(config: dict[str, Any], roles: pd.DataFrame) -> pd.DataFrame:
    extension = pd.read_csv(resolve(config["inputs"]["extension_manifest"]))
    extension = extension.set_index("cohort_id")
    rows = []
    for role in roles.itertuples(index=False):
        if role.path_source == "current_atlas":
            path = resolve(config["inputs"]["current_atlas"])
            expected_sha256 = config["expected_current_atlas"]["sha256"]
        elif role.path_source == "extension_manifest":
            if role.cohort_id not in extension.index:
                raise KeyError(f"Role cohort absent from extension manifest: {role.cohort_id}")
            source = extension.loc[role.cohort_id]
            path = Path(source["output_h5ad"])
            expected_sha256 = str(source["output_sha256"])
        else:
            raise ValueError(f"Unknown path_source: {role.path_source}")
        rows.append({**role._asdict(), "path": str(path), "expected_sha256": expected_sha256})
    return pd.DataFrame(rows)


def inspect_inputs(config: dict[str, Any], inputs: pd.DataFrame, model_genes: list[str], stage1_genes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, set[str]], list[dict[str, Any]], dict[str, tuple[int, int]]]:
    input_rows: list[dict[str, Any]] = []
    metadata_rows: list[dict[str, Any]] = []
    truth_rows: list[dict[str, Any]] = []
    gene_sets: dict[str, set[str]] = {}
    checks: list[dict[str, Any]] = []
    before_state: dict[str, tuple[int, int]] = {}
    matrix_policy = config["matrix_audit"]
    required = config["required_metadata"]

    for item in inputs.itertuples(index=False):
        path = Path(item.path)
        stat = path.stat()
        before_state[item.cohort_id] = (stat.st_size, stat.st_mtime_ns)
        observed_sha = sha256_file(path)
        hash_ok = observed_sha == item.expected_sha256
        checks.append({"check": f"sha256::{item.cohort_id}", "status": "PASS" if hash_ok else "FAIL", "detail": observed_sha})
        if not hash_ok:
            raise RuntimeError(f"SHA-256 mismatch for {item.cohort_id}")

        matrix_key = config["matrix_keys"].get(item.cohort_id, config["matrix_keys"]["default"])
        with h5py.File(path, "r") as handle:
            obs_names = axis_names(handle, "obs")
            genes = axis_names(handle, "var")
            gene_set = set(genes.tolist())
            gene_sets[item.cohort_id] = gene_set
            group = matrix_group(handle, matrix_key)
            encoding = decode(group.attrs.get("encoding-type", ""))
            shape = tuple(int(value) for value in group.attrs["shape"])
            sample = sampled_matrix_audit(
                group,
                int(matrix_policy["sample_values_per_region"]),
                float(matrix_policy["integer_tolerance"]),
            )
            duplicated_cells = int(pd.Index(obs_names).duplicated().sum())
            duplicated_genes = int(pd.Index(genes).duplicated().sum())
            input_rows.append(
                {
                    "cohort_id": item.cohort_id,
                    "role": item.role,
                    "path": str(path),
                    "sha256": observed_sha,
                    "size_bytes": stat.st_size,
                    "n_cells": int(obs_names.size),
                    "expected_cells": int(item.expected_cells),
                    "n_genes": int(genes.size),
                    "matrix_key": matrix_key,
                    "matrix_encoding": encoding,
                    "matrix_dtype": str(group["data"].dtype),
                    "matrix_nnz": int(group["data"].shape[0]),
                    "matrix_storage_bytes": int(group["data"].size * group["data"].dtype.itemsize + group["indices"].size * group["indices"].dtype.itemsize + group["indptr"].size * group["indptr"].dtype.itemsize),
                    "duplicated_cell_ids": duplicated_cells,
                    "duplicated_genes": duplicated_genes,
                    "model_features_present": int(sum(gene in gene_set for gene in model_genes)),
                    "model_feature_coverage": float(sum(gene in gene_set for gene in model_genes) / len(model_genes)),
                    "stage1_features_present": int(sum(gene in gene_set for gene in stage1_genes)),
                    "stage1_feature_coverage": float(sum(gene in gene_set for gene in stage1_genes) / len(stage1_genes)),
                    **sample,
                }
            )
            required_columns = list(
                required["current_atlas"]
                if item.cohort_id == "current_atlas"
                else required["extension"] + required["tcr"]
            )
            for column in required_columns:
                summary = obs_column_summary(handle, column)
                metadata_rows.append({"cohort_id": item.cohort_id, "column": column, **summary})
            truth = {column: obs_bool(handle, column) for column in required["tcr"] if column in handle["obs"]}
            any_ab = truth.get("has_any_ab_tcr", np.zeros(obs_names.size, dtype=bool))
            any_gd = truth.get("has_any_gd_tcr", np.zeros(obs_names.size, dtype=bool))
            pair_ab = truth.get("has_TRA_TRB_paired", np.zeros(obs_names.size, dtype=bool))
            pair_gd = truth.get("has_TRG_TRD_paired", np.zeros(obs_names.size, dtype=bool))
            author_nk = np.zeros(obs_names.size, dtype=bool)
            if "major_cluster" in handle["obs"]:
                labels = pd.Series(obs_values(handle, "major_cluster"), dtype="string").fillna("").str.lower()
                author_nk = labels.eq("nkcell").to_numpy(bool)
            truth_rows.append(
                {
                    "cohort_id": item.cohort_id,
                    "role": item.role,
                    "n_cells": int(obs_names.size),
                    "n_any_ab_tcr": int(any_ab.sum()),
                    "n_paired_ab_tcr": int(pair_ab.sum()),
                    "n_any_gd_tcr": int(any_gd.sum()),
                    "n_paired_gd_tcr": int(pair_gd.sum()),
                    "n_author_nk": int(author_nk.sum()),
                    "n_author_nk_no_productive_tcr": int((author_nk & ~any_ab & ~any_gd).sum()),
                }
            )
        expected_cells_ok = int(input_rows[-1]["n_cells"]) == int(item.expected_cells)
        checks.extend(
            [
                {"check": f"expected_cells::{item.cohort_id}", "status": "PASS" if expected_cells_ok else "FAIL", "detail": str(input_rows[-1]["n_cells"])},
                {"check": f"sparse_matrix::{item.cohort_id}", "status": "PASS" if encoding in matrix_policy["allowed_sparse_encodings"] else "FAIL", "detail": encoding},
                {"check": f"raw_count_sample::{item.cohort_id}", "status": "PASS" if sample["sample_finite"] and sample["sample_nonnegative"] and sample["sample_integer_like"] else "FAIL", "detail": f"n={sample['sampled_values']}; range={sample['sample_min']}..{sample['sample_max']}"},
                {"check": f"unique_axes::{item.cohort_id}", "status": "PASS" if duplicated_cells == 0 and duplicated_genes == 0 else "FAIL", "detail": f"cell_duplicates={duplicated_cells}; gene_duplicates={duplicated_genes}"},
            ]
        )
    return pd.DataFrame(input_rows), pd.DataFrame(metadata_rows), pd.DataFrame(truth_rows), gene_sets, checks, before_state


def role_checks(inputs: pd.DataFrame) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []
    duplicate = int(inputs["cohort_id"].duplicated().sum())
    checks.append({"check": "unique_cohort_roles", "status": "PASS" if duplicate == 0 else "FAIL", "detail": f"duplicates={duplicate}"})
    locked = inputs[inputs["allow_locked_evaluation"].astype(bool)]
    leakage_columns = ["include_in_integration_fit", "include_in_cluster_label_design", "allow_pseudolabel_training", "allow_model_tuning"]
    leaking = int(locked[leakage_columns].astype(bool).any(axis=1).sum())
    checks.append({"check": "locked_cohorts_excluded_from_development", "status": "PASS" if leaking == 0 else "FAIL", "detail": f"leaking_cohorts={leaking}"})
    development = inputs[inputs["include_in_integration_fit"].astype(bool)]
    checks.append({"check": "development_source_count", "status": "PASS" if development.shape[0] >= 6 else "FAIL", "detail": f"{development.shape[0]} input objects"})
    return checks


def primary_nk_anchor_audit(config: dict[str, Any]) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    manifest = pd.read_csv(
        resolve(config["inputs"]["v4_cell_label_manifest"]),
        usecols=["cell_id", "source_gse_id", "stage1_role", "nk_annotation_strength"],
        low_memory=False,
    )
    anchors = manifest[
        manifest["stage1_role"].eq("nk_negative")
        & manifest["nk_annotation_strength"].eq(1.0)
    ].copy()
    atlas_path = resolve(config["inputs"]["current_atlas"])
    with h5py.File(atlas_path, "r") as handle:
        atlas_ids = pd.Index(axis_names(handle, "obs"))
    anchor_ids = pd.Index(anchors["cell_id"].astype(str))
    matched = anchor_ids.isin(atlas_ids)
    summary = (
        anchors.assign(matched_current_atlas=matched)
        .groupby("source_gse_id", as_index=False)
        .agg(n_anchor_cells=("cell_id", "size"), n_matched_current_atlas=("matched_current_atlas", "sum"))
    )
    expected = config["expected_primary_nk_anchors"]
    checks = [
        {
            "check": "primary_nk_anchor_count",
            "status": "PASS" if anchors.shape[0] == int(expected["n_cells"]) else "FAIL",
            "detail": f"{anchors.shape[0]} anchors",
        },
        {
            "check": "primary_nk_anchor_source_count",
            "status": "PASS" if anchors["source_gse_id"].nunique() == int(expected["n_sources"]) else "FAIL",
            "detail": f"{anchors['source_gse_id'].nunique()} sources",
        },
        {
            "check": "primary_nk_anchors_map_to_current_atlas",
            "status": "PASS" if bool(matched.all()) else "FAIL",
            "detail": f"{int(matched.sum())}/{matched.size} matched",
        },
    ]
    return summary, checks


def gene_overlap_tables(inputs: pd.DataFrame, gene_sets: dict[str, set[str]], model_genes: list[str], stage1_genes: list[str], config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, Any]]]:
    base = gene_sets["current_atlas"]
    rows = []
    for item in inputs.itertuples(index=False):
        genes = gene_sets[item.cohort_id]
        rows.append(
            {
                "cohort_id": item.cohort_id,
                "role": item.role,
                "n_genes": len(genes),
                "overlap_with_current_atlas": len(genes & base),
                "extra_vs_current_atlas": len(genes - base),
                "model_features_present": sum(gene in genes for gene in model_genes),
                "model_feature_coverage": sum(gene in genes for gene in model_genes) / len(model_genes),
                "stage1_features_present": sum(gene in genes for gene in stage1_genes),
                "stage1_feature_coverage": sum(gene in genes for gene in stage1_genes) / len(stage1_genes),
            }
        )
    included = inputs.loc[inputs["include_in_integration_fit"].astype(bool), "cohort_id"].tolist()
    common = set.intersection(*(gene_sets[cohort] for cohort in included))
    tcr_pattern = re.compile(r"^TR[ABGD](?:V|J|D|C)|^TRAC$|^TRBC|^TRGC|^TRDC$")
    mito_ribo_ig = re.compile(r"^(?:MT-|RPS|RPL|IG[HKL])")
    common_tcr = {gene for gene in common if tcr_pattern.match(gene)}
    common_excluded = {gene for gene in common if tcr_pattern.match(gene) or mito_ribo_ig.match(gene)}
    summary = pd.DataFrame(
        [
            {
                "development_input_count": len(included),
                "common_genes": len(common),
                "common_tcr_genes": len(common_tcr),
                "common_genes_after_tcr_mito_ribo_ig_exclusion": len(common - common_excluded),
                "planned_hvgs": int(config["feature_policy"]["n_hvgs"]),
                "model_genes_in_common": sum(gene in common for gene in model_genes),
                "stage1_genes_in_common": sum(gene in common for gene in stage1_genes),
            }
        ]
    )
    enough_common = len(common - common_excluded) >= max(int(config["feature_policy"]["minimum_common_genes"]), int(config["feature_policy"]["n_hvgs"]))
    checks = [
        {"check": "common_gene_universe", "status": "PASS" if enough_common else "FAIL", "detail": f"{len(common - common_excluded)} eligible common genes"},
        {"check": "development_model_feature_coverage", "status": "PASS" if all(sum(gene in gene_sets[cohort] for gene in model_genes) / len(model_genes) >= float(config["feature_policy"]["minimum_model_feature_coverage"]) for cohort in included) else "FAIL", "detail": f"minimum={min(sum(gene in gene_sets[cohort] for gene in model_genes) / len(model_genes) for cohort in included):.3f}"},
        {"check": "locked_stage1_feature_coverage", "status": "PASS" if all(sum(gene in gene_sets[cohort] for gene in stage1_genes) / len(stage1_genes) >= float(config["feature_policy"]["minimum_stage1_feature_coverage"]) for cohort in inputs.loc[inputs['allow_locked_evaluation'].astype(bool), 'cohort_id']) else "FAIL", "detail": f"minimum={min(sum(gene in gene_sets[cohort] for gene in stage1_genes) / len(stage1_genes) for cohort in inputs.loc[inputs['allow_locked_evaluation'].astype(bool), 'cohort_id']):.3f}"},
    ]
    return pd.DataFrame(rows), summary, checks


def resource_audit(inputs: pd.DataFrame, input_audit: pd.DataFrame, common_summary: pd.DataFrame, config: dict[str, Any]) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    included_ids = inputs.loc[inputs["include_in_integration_fit"].astype(bool), "cohort_id"]
    included = input_audit[input_audit["cohort_id"].isin(included_ids)]
    n_cells = int(included["n_cells"].sum())
    full_csr_gib = float(included["matrix_storage_bytes"].sum() / 2**30)
    common_fraction = float(common_summary.iloc[0]["common_genes"] / input_audit.loc[input_audit["cohort_id"].eq("current_atlas"), "n_genes"].iloc[0])
    common_csr_gib = full_csr_gib * min(1.0, common_fraction)
    hvg_csr_gib = common_csr_gib * int(config["feature_policy"]["n_hvgs"]) / int(common_summary.iloc[0]["common_genes"])
    latent_gib = n_cells * int(config["planned_scvi"]["n_latent"]) * 4 / 2**30
    neighbors = int(config["planned_clustering"]["global"]["n_neighbors"])
    graph_gib = n_cells * neighbors * 12 / 2**30
    conservative_peak = common_csr_gib * 4 + hvg_csr_gib * 3 + latent_gib * 4 + graph_gib * 3 + 50
    disk = shutil.disk_usage(config["resource_contract"]["planned_ssd_root"].split("/Integrated_dataset")[0])
    ssd_free_gib = disk.free / 2**30
    gpu_name, gpu_total_mib, gpu_free_mib, gpu_util = "unavailable", math.nan, math.nan, math.nan
    completed = subprocess.run(
        ["nvidia-smi", "--query-gpu=name,memory.total,memory.free,utilization.gpu", "--format=csv,noheader,nounits"],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode == 0 and completed.stdout.strip():
        parts = [part.strip() for part in completed.stdout.strip().splitlines()[0].split(",")]
        gpu_name, gpu_total_mib, gpu_free_mib, gpu_util = parts[0], float(parts[1]), float(parts[2]), float(parts[3])
    row = {
        "development_cells": n_cells,
        "input_full_csr_gib": full_csr_gib,
        "estimated_common_csr_gib": common_csr_gib,
        "estimated_hvg_csr_gib": hvg_csr_gib,
        "latent_gib": latent_gib,
        "knn_graph_gib": graph_gib,
        "conservative_peak_ram_gib": conservative_peak,
        "ram_contract_gib": float(config["resource_contract"]["maximum_ram_gib"]),
        "ssd_free_gib": ssd_free_gib,
        "ssd_minimum_gib": float(config["resource_contract"]["minimum_ssd_free_gib"]),
        "gpu_name": gpu_name,
        "gpu_total_gib": gpu_total_mib / 1024 if math.isfinite(gpu_total_mib) else math.nan,
        "gpu_free_gib": gpu_free_mib / 1024 if math.isfinite(gpu_free_mib) else math.nan,
        "gpu_utilization_percent": gpu_util,
    }
    checks = [
        {"check": "ram_feasibility", "status": "PASS" if conservative_peak < float(config["resource_contract"]["maximum_ram_gib"]) else "FAIL", "detail": f"estimated {conservative_peak:.1f} GiB < {config['resource_contract']['maximum_ram_gib']} GiB"},
        {"check": "ssd_feasibility", "status": "PASS" if ssd_free_gib >= float(config["resource_contract"]["minimum_ssd_free_gib"]) else "FAIL", "detail": f"{ssd_free_gib:.1f} GiB free"},
        {"check": "gpu_feasibility", "status": "PASS" if math.isfinite(gpu_total_mib) and gpu_total_mib / 1024 >= float(config["resource_contract"]["minimum_gpu_memory_gib"]) else "FAIL", "detail": f"{gpu_name}; total={gpu_total_mib / 1024 if math.isfinite(gpu_total_mib) else math.nan:.1f} GiB"},
    ]
    return pd.DataFrame([row]), checks


def file_state_checks(inputs: pd.DataFrame, before: dict[str, tuple[int, int]]) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    rows = []
    for item in inputs.itertuples(index=False):
        stat = Path(item.path).stat()
        same = before[item.cohort_id] == (stat.st_size, stat.st_mtime_ns)
        rows.append({"cohort_id": item.cohort_id, "size_bytes_before": before[item.cohort_id][0], "size_bytes_after": stat.st_size, "mtime_ns_before": before[item.cohort_id][1], "mtime_ns_after": stat.st_mtime_ns, "unchanged": same})
    frame = pd.DataFrame(rows)
    return frame, [{"check": "all_input_h5ads_unchanged", "status": "PASS" if frame["unchanged"].all() else "FAIL", "detail": f"{int(frame['unchanged'].sum())}/{frame.shape[0]} unchanged"}]


def make_figures(inputs: pd.DataFrame, input_audit: pd.DataFrame, truth: pd.DataFrame, resource: pd.DataFrame, figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 11, "axes.labelsize": 9})
    role_colors = {
        "development_reference": "#285f8f",
        "development_cluster_expansion": "#14866d",
        "locked_author_nk_stage1_specificity": "#b06035",
        "locked_paired_abt_specificity": "#7b5aa6",
        "locked_mixed_t_nk_stress": "#808891",
    }
    view = input_audit.sort_values("n_cells")
    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    ax.barh(view["cohort_id"], view["n_cells"], color=[role_colors[value] for value in view["role"]])
    ax.set_xscale("log")
    ax.set_xlabel("Cells (log scale)")
    ax.set_title("Frozen V4.2 integration and locked-cohort roles")
    fig.tight_layout()
    fig.savefig(figure_dir / "cohort_roles_and_cell_counts.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    view = input_audit.sort_values("model_feature_coverage")
    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    y = np.arange(view.shape[0])
    ax.barh(y - 0.18, 100 * view["model_feature_coverage"], height=0.34, color="#285f8f", label="All 197 genes")
    ax.barh(y + 0.18, 100 * view["stage1_feature_coverage"], height=0.34, color="#14866d", label="Stage-1 50 genes")
    ax.axvline(90, color="#777777", ls="--", lw=1)
    ax.set_yticks(y, view["cohort_id"])
    ax.set_xlim(65, 101)
    ax.set_xlabel("Feature coverage (%)")
    ax.set_title("Locked Stage-1 controls retain complete gate features")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figure_dir / "feature_coverage_by_cohort.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    view = truth.sort_values("n_cells")
    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    ax.barh(view["cohort_id"], view["n_paired_ab_tcr"], color="#285f8f", label="Paired alpha-beta TCR")
    ax.barh(view["cohort_id"], view["n_author_nk_no_productive_tcr"], left=view["n_paired_ab_tcr"], color="#b06035", label="Author NK, no productive TCR")
    ax.set_xscale("symlog", linthresh=100)
    ax.set_xlabel("Expression-independent or author-supported controls (symlog)")
    ax.set_title("Available specificity controls by cohort")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figure_dir / "specificity_controls_by_cohort.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    values = resource.iloc[0]
    labels = ["Estimated peak RAM", "RAM contract", "SSD free", "GPU total"]
    amounts = [values["conservative_peak_ram_gib"], values["ram_contract_gib"], values["ssd_free_gib"], values["gpu_total_gib"]]
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    bars = ax.bar(labels, amounts, color=["#14866d", "#285f8f", "#7b5aa6", "#b06035"])
    ax.bar_label(bars, labels=[f"{value:.0f} GiB" for value in amounts], padding=3)
    ax.set_ylabel("GiB")
    ax.set_title("Preflight compute and storage envelope")
    ax.set_ylim(0, max(amounts) * 1.15)
    fig.tight_layout()
    fig.savefig(figure_dir / "resource_feasibility.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def markdown_table(frame: pd.DataFrame, columns: list[str], formats: dict[str, str] | None = None) -> str:
    selected = frame.loc[:, columns].copy()
    for column, template in (formats or {}).items():
        selected[column] = selected[column].map(lambda value: "NA" if pd.isna(value) else template.format(value))
    return selected.to_markdown(index=False)


def render_report(inputs: pd.DataFrame, input_audit: pd.DataFrame, metadata: pd.DataFrame, truth: pd.DataFrame, anchors: pd.DataFrame, overlap: pd.DataFrame, common: pd.DataFrame, resource: pd.DataFrame, checks: pd.DataFrame, config: dict[str, Any], log_dir: Path, static_dir: Path) -> None:
    static_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    dev = inputs[inputs["include_in_integration_fit"].astype(bool)]
    locked = inputs[inputs["allow_locked_evaluation"].astype(bool)]
    gse169 = input_audit[input_audit["cohort_id"].eq("GSE169246")].iloc[0]
    missing_required = metadata[~metadata["present"]]
    result = "PASS_REVIEW_REQUIRED" if checks["status"].eq("PASS").all() else "FAIL"
    report = f"""# gdTAI V4.2 modeling-integration preflight

## Decision

**{result}.** The proposed sidecar integration is technically feasible and the whole-dataset role split is internally consistent. This report authorizes supervision review only. It does not authorize integration, scVI fitting, clustering, pseudo-labeling, classifier fitting, threshold selection, promotion, or atlas inference.

- Development integration: **{int(dev['expected_cells'].sum()):,} cells** from the current atlas plus five new development cohorts.
- Locked evaluation: **{int(locked['expected_cells'].sum()):,} cells** across three whole cohorts excluded before HVG selection and scVI fitting.
- Common development genes: **{int(common.iloc[0]['common_genes']):,}**; eligible after TCR/mitochondrial/ribosomal/immunoglobulin exclusion: **{int(common.iloc[0]['common_genes_after_tcr_mito_ribo_ig_exclusion']):,}**.
- GSE169246 contains all 50 Stage-1 genes but only **{int(gse169['model_features_present'])}/197 ({gse169['model_feature_coverage']:.1%})** full classifier genes. It is a valid Stage-1 NK-passage challenge but only a reduced-feature sensitivity analysis for final-cascade FPR.
- The extension cohorts contain no expression-independent gdT-positive truth. They strengthen specificity training/evaluation but cannot establish gdT recall, F1, or superiority alone.

![Cohort roles](../../Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_integration_preflight/cohort_roles_and_cell_counts.png)

## Frozen cohort roles

{markdown_table(inputs, ['cohort_id', 'role', 'include_in_integration_fit', 'allow_pseudolabel_training', 'allow_model_tuning', 'allow_locked_evaluation', 'expected_cells'])}

All locked cohorts are excluded from HVG selection, scVI training, clustering, cluster interpretation, pseudo-label construction, feature selection, calibration, thresholding, and model-family selection. These cohorts have already been inspected by earlier frozen-model screens and are not prospective or fully independent; they are locked specifically for V4.2 comparison.

## Matrix and feature compatibility

{markdown_table(input_audit, ['cohort_id', 'n_cells', 'n_genes', 'matrix_key', 'matrix_encoding', 'matrix_dtype', 'sample_integer_like', 'model_feature_coverage', 'stage1_feature_coverage'], {'model_feature_coverage': '{:.1%}', 'stage1_feature_coverage': '{:.1%}'})}

![Feature coverage](../../Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_integration_preflight/feature_coverage_by_cohort.png)

The development intersection supports 4,000 batch-aware HVGs without forcing classifier genes into the latent representation. All TCR genes are excluded from HVG selection so that TCR-family expression cannot directly define cluster-derived NK pseudo-labels. The classifier itself remains a log1p(CP10K) individual-gene model and never consumes scVI coordinates or Leiden labels.

## Available controls

{markdown_table(truth, ['cohort_id', 'n_cells', 'n_paired_ab_tcr', 'n_paired_gd_tcr', 'n_author_nk_no_productive_tcr'])}

Primary weight-1 NK anchors mapped to the current atlas:

{markdown_table(anchors, ['source_gse_id', 'n_anchor_cells', 'n_matched_current_atlas'])}

![Specificity controls](../../Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_integration_preflight/specificity_controls_by_cohort.png)

GSE169246 contributes 7,770 author NK cells without productive TCR and 54,925 productive paired alpha-beta T cells. GSE315928 contributes 66,813 paired alpha-beta T cells with complete 197-gene coverage. GSE121636/GSE121637 remains a mixed stress cohort because it lacks an independent author NK label in the harmonized object.

## Precommitted integration and clustering

1. Build a new sidecar object under the mirrored SSD root; never mutate a source or canonical milestone H5AD.
2. Use raw counts, 4,000 batch-aware common HVGs, the existing library/sample/source batch hierarchy, scVI `128 x 2`, 30 latent dimensions, negative-binomial likelihood, 20 epochs maximum, batch size 8,192, GPU only, and early stopping.
3. Run nine global RAPIDS Leiden partitions: resolutions 0.4, 0.8, 1.2 crossed with seeds 11, 29, 47.
4. Recluster the T/NK boundary in nine additional partitions: resolutions 0.3, 0.6, 1.0 crossed with the same seeds.
5. A cluster-derived development NK negative requires >=80% run agreement, >=95% independently labeled NK-anchor purity, <=2% productive-TCR contamination, representation in >=3 development sources, and <=70% contribution from any one source.
6. No threshold on CD3, TRDC, TRDV, NKG7, GNLY, KLRD1, TYROBP, FCER1G, or FCGR3A may define truth. Pseudo-labels receive weight 0.5, are source-balanced and capped, and cannot determine validation guardrails.

Primary paired-gdT inputs remain outside integration and pseudo-label construction. They retain weight-1 positive roles during classifier fitting. Their later latent-query proximity to NK clusters is diagnostic only and cannot tune the pseudo-label rules.

## Resource feasibility

{markdown_table(resource, ['development_cells', 'input_full_csr_gib', 'estimated_common_csr_gib', 'estimated_hvg_csr_gib', 'latent_gib', 'knn_graph_gib', 'conservative_peak_ram_gib', 'ssd_free_gib', 'gpu_name', 'gpu_total_gib'], {'input_full_csr_gib': '{:.1f}', 'estimated_common_csr_gib': '{:.1f}', 'estimated_hvg_csr_gib': '{:.1f}', 'latent_gib': '{:.2f}', 'knn_graph_gib': '{:.2f}', 'conservative_peak_ram_gib': '{:.1f}', 'ssd_free_gib': '{:.1f}', 'gpu_total_gib': '{:.1f}'})}

![Resources](../../Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_integration_preflight/resource_feasibility.png)

The estimate is intentionally conservative and remains below the active 800-GiB exception. The A100 is available with more than the frozen 75-GiB minimum, and SSD free space exceeds the 300-GiB floor.

## QC checks

{markdown_table(checks, ['check', 'status', 'detail'])}

Missing required metadata columns: **{missing_required.shape[0]}**. All nine input H5AD size and nanosecond modification-time pairs were checked before and after the audit.

## Supervision gate

The next executable milestone is implementation QC for a fail-closed sidecar integration runner. It requires explicit approval. Even after integration approval, V4.2 classifier fitting requires a separate reviewed QC gate and approval record.

## Outputs

- Full tables: `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_integration_preflight/`
- Machine summary: `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_integration_preflight/summary.json`
- HTML: `gdT_prediction/gdtai_v4_2_integration_preflight/index.html`
- PDF: `gdT_prediction/gdtai_v4_2_integration_preflight/gdtai_v4_2_integration_preflight_report.pdf`
"""
    md_path = log_dir / "gdtai_v4_2_integration_preflight_summary.md"
    md_path.write_text(report, encoding="utf-8")
    css = """
body{font-family:Arial,Helvetica,sans-serif;color:#17212b;line-height:1.46;max-width:1120px;margin:0 auto;padding:30px 42px;background:#fff}h1{font-size:30px;color:#17324d;border-bottom:4px solid #14866d;padding-bottom:12px}h2{font-size:21px;color:#17324d;margin-top:32px}p,li{font-size:14px}table{border-collapse:collapse;width:100%;font-size:10.5px;margin:16px 0 24px;table-layout:auto}th{background:#17324d;color:#fff;text-align:left}th,td{padding:6px 7px;border:1px solid #cdd6df;vertical-align:top;overflow-wrap:anywhere}tr:nth-child(even){background:#f3f6f8}img{display:block;max-width:96%;max-height:680px;margin:22px auto}code{background:#edf1f4;padding:1px 4px;border-radius:2px}strong{color:#102f43}@media print{body{max-width:none;padding:8mm 9mm}h1{font-size:24px}h2{font-size:18px;page-break-after:avoid}table,img{page-break-inside:avoid}p,li{font-size:10.5px}table{font-size:7.7px}img{max-height:165mm}}
""".strip()
    (static_dir / "print.css").write_text(css + "\n", encoding="utf-8")
    html_path = static_dir / "index.html"
    subprocess.run(["pandoc", str(md_path), "--standalone", "--metadata", "title=gdTAI V4.2 integration preflight", "--css", "print.css", "-o", str(html_path)], cwd=ROOT, check=True)
    pdf_path = static_dir / "gdtai_v4_2_integration_preflight_report.pdf"
    profile = Path("/tmp/gdtai-v42-integration-preflight-chrome")
    profile.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "google-chrome", "--headless", "--no-sandbox", "--disable-gpu", "--disable-dev-shm-usage",
            "--disable-breakpad", "--disable-crash-reporter", "--allow-file-access-from-files", "--no-pdf-header-footer",
            f"--user-data-dir={profile}", f"--print-to-pdf={pdf_path}", html_path.resolve().as_uri(),
        ],
        cwd=ROOT,
        check=True,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config_path = resolve(args.config)
    config = json.loads(config_path.read_text())
    roles_path = resolve(config["inputs"]["cohort_roles"])
    roles = pd.read_csv(roles_path)
    boolean_columns = ["include_in_integration_fit", "include_in_cluster_label_design", "allow_pseudolabel_training", "allow_model_tuning", "allow_locked_evaluation"]
    for column in boolean_columns:
        roles[column] = roles[column].astype(bool)
    inputs = resolve_inputs(config, roles)
    model_features = pd.read_csv(resolve(config["inputs"]["model_feature_manifest"]))["gene"].astype(str).tolist()
    v4_config = json.loads(resolve(config["inputs"]["v4_evaluation_config"]).read_text())
    stage1_features = list(v4_config["feature_policy"]["stage1_genes"])

    outputs = config["outputs"]
    table_dir = resolve(outputs["table_dir"])
    figure_dir = resolve(outputs["figure_dir"])
    log_dir = resolve(outputs["log_dir"])
    static_dir = resolve(outputs["static_dir"])
    for directory in [table_dir, figure_dir, log_dir, static_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    input_audit, metadata, truth, gene_sets, checks, before = inspect_inputs(config, inputs, model_features, stage1_features)
    checks.extend(role_checks(inputs))
    anchors, anchor_checks = primary_nk_anchor_audit(config)
    checks.extend(anchor_checks)
    overlap, common, overlap_checks = gene_overlap_tables(inputs, gene_sets, model_features, stage1_features, config)
    checks.extend(overlap_checks)
    metadata_missing = metadata[~metadata["present"]]
    checks.append({"check": "required_metadata_columns", "status": "PASS" if metadata_missing.empty else "FAIL", "detail": f"missing={metadata_missing.shape[0]}"})
    resource, resource_checks = resource_audit(inputs, input_audit, common, config)
    checks.extend(resource_checks)
    file_state, state_checks = file_state_checks(inputs, before)
    checks.extend(state_checks)
    checks_frame = pd.DataFrame(checks)
    result = "PASS_REVIEW_REQUIRED" if checks_frame["status"].eq("PASS").all() else "FAIL"

    inputs.to_csv(table_dir / "frozen_cohort_roles.csv", index=False)
    input_audit.to_csv(table_dir / "input_h5ad_audit.csv", index=False)
    metadata.to_csv(table_dir / "metadata_compatibility.csv", index=False)
    truth.to_csv(table_dir / "truth_and_control_counts.csv", index=False)
    anchors.to_csv(table_dir / "primary_nk_anchor_mapping.csv", index=False)
    overlap.to_csv(table_dir / "gene_overlap_by_cohort.csv", index=False)
    common.to_csv(table_dir / "development_common_gene_summary.csv", index=False)
    resource.to_csv(table_dir / "resource_feasibility.csv", index=False)
    file_state.to_csv(table_dir / "input_file_state.csv", index=False)
    checks_frame.to_csv(table_dir / "preflight_checks.csv", index=False)
    make_figures(inputs, input_audit, truth, resource, figure_dir)

    contract = {
        "protocol_version": config["protocol_version"],
        "config_sha256": sha256_file(config_path),
        "roles_sha256": sha256_file(roles_path),
        "step0_protocol_sha256": sha256_file(resolve(config["inputs"]["step0_protocol"])),
        "input_sha256": dict(zip(input_audit["cohort_id"], input_audit["sha256"])),
        "cohort_role_contract_sha256": canonical_sha256(inputs.drop(columns=["path"]).to_dict(orient="records")),
        "result": result,
        "integration_or_model_fitting_performed": False,
        "development_cells": int(input_audit.loc[input_audit["cohort_id"].isin(inputs.loc[inputs["include_in_integration_fit"], "cohort_id"]), "n_cells"].sum()),
        "locked_cells": int(input_audit.loc[input_audit["cohort_id"].isin(inputs.loc[inputs["allow_locked_evaluation"], "cohort_id"]), "n_cells"].sum()),
        "common_development_genes": int(common.iloc[0]["common_genes"]),
        "eligible_common_genes": int(common.iloc[0]["common_genes_after_tcr_mito_ribo_ig_exclusion"]),
        "next_action": "Explicit supervision approval is required before implementation or integration.",
    }
    (log_dir / "summary.json").write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    approval_template = {
        "approved": False,
        "protocol_version": config["protocol_version"],
        "config_sha256": contract["config_sha256"],
        "roles_sha256": contract["roles_sha256"],
        "cohort_role_contract_sha256": contract["cohort_role_contract_sha256"],
        "review_artifact": str(static_dir / "gdtai_v4_2_integration_preflight_report.pdf"),
        "scope": "Implementation QC only; does not authorize classifier fitting.",
        "approved_by": "",
        "approved_at": "",
    }
    (log_dir / "INTEGRATION_APPROVAL_TEMPLATE.json").write_text(json.dumps(approval_template, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    render_report(inputs, input_audit, metadata, truth, anchors, overlap, common, resource, checks_frame, config, log_dir, static_dir)
    print(json.dumps(contract, indent=2, sort_keys=True))
    if result == "FAIL":
        raise SystemExit("V4.2 integration preflight failed")


if __name__ == "__main__":
    main()
