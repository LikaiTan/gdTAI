#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Build a curated gdT biology atlas from gdTAI predictions.

The workflow is read-only with respect to the source H5AD. It starts from the
gdTAI v2 NK-optimized prediction mask, removes likely false positives, adds
evidence-backed false negatives, scores phenotype marker programs and V-gene
RNA expression, and renders a static HTML/PDF report.
"""

from __future__ import annotations

# TNK_WORKFLOW_BOOTSTRAP
import sys as _tnk_sys
from pathlib import Path as _TnkPath

_TNK_PROJECT_ROOT = _TnkPath(__file__).resolve().parents[2]
for _tnk_path in (
    _TNK_PROJECT_ROOT,
    _TNK_PROJECT_ROOT / "src",
    _TNK_PROJECT_ROOT / "workflows" / "integration",
    _TNK_PROJECT_ROOT / "workflows" / "intake",
    _TNK_PROJECT_ROOT / "workflows" / "metadata",
    _TNK_PROJECT_ROOT / "workflows" / "analysis",
    _TNK_PROJECT_ROOT / "workflows" / "gdtai",
    _TNK_PROJECT_ROOT / "workflows" / "gdt_atlas",
    _TNK_PROJECT_ROOT / "workflows" / "reporting",
    _TNK_PROJECT_ROOT / "workflows" / "maintenance",
):
    _tnk_value = str(_tnk_path)
    if _tnk_value not in _tnk_sys.path:
        _tnk_sys.path.insert(0, _tnk_value)


import argparse
import html
import json
import logging
import math
import subprocess
from pathlib import Path
from typing import Any

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from run_gdtai_nk_optimized_gdt_atlas import (
    DEFAULT_H5AD,
    DEFAULT_METADATA,
    DEFAULT_PREDICTIONS,
    DEFAULT_RULES,
    PROJECT_ROOT,
    add_harmonized_metadata,
    apply_metadata_lookup,
    clean_series,
    dataframe_to_html,
    dataframe_to_markdown,
    load_metadata_lookup,
    read_bool_obs,
    read_obs_column,
    read_string_dataset,
)


TABLE_DIR = PROJECT_ROOT / "Integrated_dataset/tables/gdT_atlas_curated_phenotypes"
FIGURE_DIR = PROJECT_ROOT / "Integrated_dataset/figures/gdT_atlas_curated_phenotypes"
LOG_DIR = PROJECT_ROOT / "Integrated_dataset/logs/gdT_atlas_curated_phenotypes"
STATIC_DIR = PROJECT_ROOT / "gdT_atlas/curated_phenotypes"
REPORT_HTML = STATIC_DIR / "index.html"
REPORT_PDF = STATIC_DIR / "gdT_atlas_curated_phenotypes_report.pdf"
SUMMARY_MD = LOG_DIR / "gdt_atlas_curated_phenotypes_summary.md"
MANIFEST_JSON = LOG_DIR / "gdt_atlas_curated_phenotypes_manifest.json"
RUN_LOG = LOG_DIR / "run.log"

DEFAULT_MARKERS = PROJECT_ROOT / "configs" / "gdt_atlas" / "phenotype_marker_panels.json"
DEFAULT_FN_TABLE = (
    PROJECT_ROOT
    / "Integrated_dataset/tables/gdT_prediction/gdtai_annotation_specific_tp_fn_audit/"
    / "validation_positive_tp_fn_cell_table.csv.gz"
)
PREDICTION_KEY = "annotation_specific_pred"
EXPECTED_N_OBS = 5_128_904
EXPECTED_PREDICTED = 359_857

OBS_COLUMNS = [
    "Sorted_gdT",
    "TCRseq",
    "TRA_cdr3",
    "TRB_cdr3",
    "TRG_cdr3",
    "TRD_cdr3",
    "TRA_v",
    "TRB_v",
    "TRG_v",
    "TRD_v",
    "condition",
    "donor_id",
    "has_TRA_TRB_paired",
    "has_TRG_TRD_paired",
    "has_any_ab_tcr",
    "has_any_gd_tcr",
    "library_id",
    "sample_id",
    "sampleid",
    "simple_annotation_plus6",
    "source_gse_id",
    "tissue",
    "tissue_corrected",
    "phase4_trd_score",
    "phase4_trab_score",
    "phase4_trd_minus_trab",
]

BOOL_COLUMNS = {
    "Sorted_gdT",
    "TCRseq",
    "has_TRA_TRB_paired",
    "has_TRG_TRD_paired",
    "has_any_ab_tcr",
    "has_any_gd_tcr",
}
FLOAT_COLUMNS = {"phase4_trd_score", "phase4_trab_score", "phase4_trd_minus_trab"}
INVALID_STRINGS = {"", "nan", "none", "na", "n/a", "<na>", "null", "unknown"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build curated gdT phenotype atlas report.")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--prediction-npz", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--rules", type=Path, default=DEFAULT_RULES)
    parser.add_argument("--marker-config", type=Path, default=DEFAULT_MARKERS)
    parser.add_argument("--metadata", type=Path, action="append", default=DEFAULT_METADATA)
    parser.add_argument("--fn-table", type=Path, default=DEFAULT_FN_TABLE)
    parser.add_argument("--expression-row-block", type=int, default=5000)
    parser.add_argument("--umap-background", type=int, default=300_000)
    parser.add_argument("--umap-primary", type=int, default=160_000)
    parser.add_argument("--seed", type=int, default=20260618)
    parser.add_argument("--no-pdf", action="store_true")
    return parser.parse_args()


def ensure_dirs() -> None:
    for path in [TABLE_DIR, FIGURE_DIR, LOG_DIR, STATIC_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def setup_logging() -> None:
    ensure_dirs()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(RUN_LOG, mode="w", encoding="utf-8"), logging.StreamHandler()],
        force=True,
    )


def load_obs(path: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    stat = path.stat()
    with h5py.File(path, "r") as handle:
        n_obs = int(handle["obs"]["_index"].shape[0])
        obs_index = read_string_dataset(handle["obs"]["_index"])
        data: dict[str, Any] = {"obs_name": obs_index, "obs_row": np.arange(n_obs, dtype=np.int64)}
        for column in OBS_COLUMNS:
            if column in BOOL_COLUMNS:
                data[column] = read_bool_obs(handle, column)
            else:
                data[column] = read_obs_column(handle, column)
        has_umap = "obsm" in handle and "X_umap" in handle["obsm"]
        has_scvi = "obsm" in handle and "X_scVI" in handle["obsm"]
    df = pd.DataFrame(data)
    for column in df.columns:
        if column in BOOL_COLUMNS or column == "obs_row":
            continue
        if column in FLOAT_COLUMNS:
            df[column] = pd.to_numeric(df[column], errors="coerce").fillna(0.0).astype(np.float32)
        elif df[column].dtype == object or str(df[column].dtype).startswith("string"):
            df[column] = clean_series(df[column]).astype(str)
    info = {
        "input_h5ad": str(path),
        "n_obs": n_obs,
        "has_umap": has_umap,
        "has_scvi": has_scvi,
        "size_before": stat.st_size,
        "mtime_ns_before": stat.st_mtime_ns,
    }
    return df, info


def load_prediction_arrays(path: Path, n_obs: int) -> dict[str, np.ndarray]:
    with np.load(path) as arrays:
        out = {key: np.asarray(arrays[key]) for key in arrays.files}
    if PREDICTION_KEY not in out:
        raise KeyError(f"Missing prediction key: {PREDICTION_KEY}")
    pred = out[PREDICTION_KEY].astype(bool)
    if pred.shape[0] != n_obs:
        raise RuntimeError(f"Prediction length {pred.shape[0]:,} does not match n_obs {n_obs:,}")
    if n_obs != EXPECTED_N_OBS:
        raise RuntimeError(f"Unexpected n_obs {n_obs:,}; expected {EXPECTED_N_OBS:,}")
    if int(pred.sum()) != EXPECTED_PREDICTED:
        raise RuntimeError(f"Unexpected predicted count {int(pred.sum()):,}; expected {EXPECTED_PREDICTED:,}")
    out[PREDICTION_KEY] = pred
    for key in ["current_pred", "original_pred"]:
        if key in out:
            out[key] = out[key].astype(bool)
    return out


def nonempty_text(series: pd.Series) -> pd.Series:
    return clean_series(series).astype(str).str.lower().map(lambda x: x not in INVALID_STRINGS)


def load_gold_fn_rows(path: Path, n_obs: int) -> np.ndarray:
    mask = np.zeros(n_obs, dtype=bool)
    if not path.exists():
        logging.warning("FN audit table not found: %s", path)
        return mask
    cols = pd.read_csv(path, nrows=0).columns.tolist()
    usecols = [col for col in ["obs_row", "outcome", "validation_component"] if col in cols]
    if "obs_row" not in usecols or "outcome" not in usecols:
        logging.warning("FN audit table lacks obs_row/outcome: %s", path)
        return mask
    table = pd.read_csv(path, usecols=usecols)
    fn = table["outcome"].astype(str).eq("FN")
    rows = pd.to_numeric(table.loc[fn, "obs_row"], errors="coerce").dropna().astype(np.int64)
    rows = rows[(rows >= 0) & (rows < n_obs)]
    mask[rows.to_numpy()] = True
    return mask


def build_curation(df: pd.DataFrame, pred: np.ndarray, fn_table: Path) -> pd.DataFrame:
    n_obs = df.shape[0]
    has_trg = nonempty_text(df["TRG_cdr3"])
    has_trd = nonempty_text(df["TRD_cdr3"])
    has_ab = df["has_any_ab_tcr"].astype(bool).to_numpy()
    has_paired_ab = df["has_TRA_TRB_paired"].astype(bool).to_numpy()
    has_paired_gd = df["has_TRG_TRD_paired"].astype(bool).to_numpy()
    has_any_gd = df["has_any_gd_tcr"].astype(bool).to_numpy()
    strong_gd_evidence = has_paired_gd | has_any_gd | has_trg.to_numpy() | has_trd.to_numpy()
    strong_gd_evidence |= df["phase4_trd_score"].to_numpy(dtype=float) > 0.35

    is_nk = df["simple_annotation_plus6"].astype(str).str.upper().str.contains("NK", regex=False).to_numpy()
    paired_ab_no_gd = pred & has_paired_ab & ~strong_gd_evidence
    nk_no_gd = pred & is_nk & ~strong_gd_evidence

    sorted_gold = df["Sorted_gdT"].astype(bool).to_numpy()
    paired_gd_gold = has_paired_gd & ~has_ab
    audit_fn_gold = load_gold_fn_rows(fn_table, n_obs)
    known_fn_gold = (~pred) & (sorted_gold | paired_gd_gold | audit_fn_gold)
    known_fn_silver = (~pred) & (~known_fn_gold) & ((has_any_gd | has_trg.to_numpy() | has_trd.to_numpy()) & ~has_ab)

    out = pd.DataFrame(
        {
            "gdtai_predicted_gdT": pred,
            "strong_gd_evidence": strong_gd_evidence,
            "fp_paired_TCRAB_no_gd_evidence": paired_ab_no_gd,
            "fp_NK_annotation_no_gd_evidence": nk_no_gd,
            "fn_gold_sorted_or_paired_or_validation": known_fn_gold,
            "fn_silver_unpaired_gd_evidence": known_fn_silver,
        }
    )
    out["gdt_atlas_removed_fp_pre_expression"] = out[
        ["fp_paired_TCRAB_no_gd_evidence", "fp_NK_annotation_no_gd_evidence"]
    ].any(axis=1)
    out["gdt_atlas_added_fn"] = out["fn_gold_sorted_or_paired_or_validation"] | out["fn_silver_unpaired_gd_evidence"]
    return out


def read_var_names(path: Path) -> list[str]:
    with h5py.File(path, "r") as handle:
        return [str(x) for x in read_string_dataset(handle["var"]["_index"])]


def ordered_unique(items: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for item in items:
        if item not in seen:
            out.append(item)
            seen.add(item)
    return out


def load_selected_expression(
    h5ad_path: Path,
    selected_rows: np.ndarray,
    genes: list[str],
    block_rows: int,
) -> tuple[pd.DataFrame, np.ndarray, list[str], np.ndarray]:
    var_names = read_var_names(h5ad_path)
    gene_to_idx = {gene: i for i, gene in enumerate(var_names)}
    present = [gene for gene in genes if gene in gene_to_idx]
    missing = [gene for gene in genes if gene not in gene_to_idx]
    if missing:
        logging.warning("Missing marker genes: %s", ", ".join(missing))
    if not present:
        raise RuntimeError("No requested marker genes are present in H5AD")

    selected_rows = np.asarray(selected_rows, dtype=np.int64)
    n_selected = int(selected_rows.size)
    n_genes = len(present)
    counts = np.zeros((n_selected, n_genes), dtype=np.float32)
    total_counts = np.zeros(n_selected, dtype=np.float32)

    with h5py.File(h5ad_path, "r") as handle:
        n_obs = int(handle["obs"]["_index"].shape[0])
        indptr = handle["X"]["indptr"][:]
        indices_ds = handle["X"]["indices"]
        data_ds = handle["X"]["data"]
        selected_mask = np.zeros(n_obs, dtype=bool)
        row_to_out = np.full(n_obs, -1, dtype=np.int32)
        selected_mask[selected_rows] = True
        row_to_out[selected_rows] = np.arange(n_selected, dtype=np.int32)
        gene_col = np.full(len(var_names), -1, dtype=np.int16)
        gene_col[[gene_to_idx[g] for g in present]] = np.arange(n_genes, dtype=np.int16)

        for block_start in range(0, n_obs, block_rows):
            block_end = min(n_obs, block_start + block_rows)
            if not selected_mask[block_start:block_end].any():
                continue
            start = int(indptr[block_start])
            end = int(indptr[block_end])
            if end <= start:
                continue
            idx = indices_ds[start:end]
            data = data_ds[start:end].astype(np.float32, copy=False)
            lengths = np.diff(indptr[block_start : block_end + 1]).astype(np.int64, copy=False)
            local_rows = np.arange(block_start, block_end, dtype=np.int32)
            row_ids = np.repeat(local_rows, lengths)

            in_selected = selected_mask[row_ids]
            if in_selected.any():
                out_rows = row_to_out[row_ids[in_selected]]
                np.add.at(total_counts, out_rows, data[in_selected])

            cols = gene_col[idx]
            keep = in_selected & (cols >= 0)
            if keep.any():
                out_rows = row_to_out[row_ids[keep]]
                out_cols = cols[keep].astype(np.int64, copy=False)
                counts[out_rows, out_cols] = data[keep]
            if block_start % (block_rows * 50) == 0:
                logging.info("Expression extraction rows %s/%s", f"{block_start:,}", f"{n_obs:,}")

    norm = np.log1p(counts / np.maximum(total_counts[:, None], 1.0) * 10_000.0)
    expr = pd.DataFrame(norm, columns=present)
    return expr, total_counts, present, np.asarray(missing, dtype=object)


def score_programs(expr: pd.DataFrame, marker_config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    score_cols: dict[str, pd.Series] = {}
    manifest_rows = []
    for program, genes in marker_config["programs"].items():
        present = [gene for gene in genes if gene in expr.columns]
        missing = [gene for gene in genes if gene not in expr.columns]
        manifest_rows.append(
            {"program": program, "present_genes": ",".join(present), "missing_genes": ",".join(missing), "n_present": len(present)}
        )
        if present:
            score_cols[f"score_{program}"] = expr[present].mean(axis=1)
        else:
            score_cols[f"score_{program}"] = pd.Series(0.0, index=expr.index)
    return pd.DataFrame(score_cols), pd.DataFrame(manifest_rows)


def robust_z(scores: pd.DataFrame, mask: np.ndarray, columns: list[str]) -> pd.DataFrame:
    out = pd.DataFrame(index=scores.index)
    for column in columns:
        values = scores.loc[mask, column].to_numpy(dtype=float)
        median = float(np.nanmedian(values)) if values.size else 0.0
        mad = float(np.nanmedian(np.abs(values - median))) if values.size else 0.0
        denom = 1.4826 * mad if mad > 1e-6 else float(np.nanstd(values) + 1e-6)
        out[column.replace("score_", "z_")] = (scores[column].to_numpy(dtype=float) - median) / denom
    return out


def assign_phenotypes(scores: pd.DataFrame, primary_mask: np.ndarray) -> pd.DataFrame:
    phenotype_cols = [
        "score_cytotoxic",
        "score_naive_memory",
        "score_tissue_resident",
        "score_il17_inflammatory",
        "score_ifn_response",
        "score_proliferating",
        "score_activation_exhaustion",
        "score_nk_like",
        "score_regulatory_penalty",
        "score_stress_low_quality",
    ]
    z = robust_z(scores, primary_mask, [c for c in phenotype_cols if c in scores.columns])
    rename = {
        "z_cytotoxic": "cytotoxic",
        "z_naive_memory": "naive_memory",
        "z_tissue_resident": "tissue_resident",
        "z_il17_inflammatory": "il17_inflammatory",
        "z_ifn_response": "ifn_response",
        "z_proliferating": "proliferating",
        "z_activation_exhaustion": "activation_exhaustion",
        "z_nk_like": "nk_like",
        "z_regulatory_penalty": "regulatory_like",
        "z_stress_low_quality": "stress_low_quality",
    }
    z_for_call = z.rename(columns=rename)
    call_cols = list(z_for_call.columns)
    arr = z_for_call[call_cols].to_numpy(dtype=float)
    order = np.argsort(arr, axis=1)
    best_idx = order[:, -1]
    second_idx = order[:, -2] if arr.shape[1] > 1 else best_idx
    best = arr[np.arange(arr.shape[0]), best_idx]
    second = arr[np.arange(arr.shape[0]), second_idx]
    state = np.asarray([call_cols[i] for i in best_idx], dtype=object)
    state[best < 0.25] = "gdT_other"

    if "score_myeloid_b_contamination" in scores:
        contam = (
            (scores["score_myeloid_b_contamination"].to_numpy(dtype=float) > 0.6)
            & (scores.get("score_core_t_gdt", pd.Series(0, index=scores.index)).to_numpy(dtype=float) < 0.4)
        )
        state[contam] = "contaminant_like"
    confidence = best - second
    out = pd.DataFrame({"gdt_phenotype_state": state, "gdt_phenotype_confidence": confidence})
    return pd.concat([out, z], axis=1)


def apply_expression_fp_flags(curation: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    no_gd = ~curation["strong_gd_evidence"].to_numpy(dtype=bool)
    predicted = curation["gdtai_predicted_gdT"].to_numpy(dtype=bool)
    core = scores.get("score_core_t_gdt", pd.Series(0, index=scores.index)).to_numpy(dtype=float)
    nk_like = scores.get("score_nk_like", pd.Series(0, index=scores.index)).to_numpy(dtype=float)
    reg = scores.get("score_regulatory_penalty", pd.Series(0, index=scores.index)).to_numpy(dtype=float)
    contam = scores.get("score_myeloid_b_contamination", pd.Series(0, index=scores.index)).to_numpy(dtype=float)
    trdc = scores.get("score_core_t_gdt", pd.Series(0, index=scores.index)).to_numpy(dtype=float)

    curation["fp_expression_low_cd3_no_gd"] = predicted & no_gd & (core < 0.15)
    curation["fp_expression_nk_like_no_gd"] = predicted & no_gd & (nk_like > 1.0) & (core < 0.45)
    curation["fp_expression_regulatory_no_gd"] = predicted & no_gd & (reg > 0.7) & (trdc < 0.55)
    curation["fp_expression_myeloid_b_no_gd"] = predicted & no_gd & (contam > 0.6)
    fp_cols = [
        "gdt_atlas_removed_fp_pre_expression",
        "fp_expression_low_cd3_no_gd",
        "fp_expression_nk_like_no_gd",
        "fp_expression_regulatory_no_gd",
        "fp_expression_myeloid_b_no_gd",
    ]
    curation["gdt_atlas_removed_fp"] = curation[fp_cols].any(axis=1)
    curation["gdt_atlas_primary_cell"] = (
        curation["gdtai_predicted_gdT"] & ~curation["gdt_atlas_removed_fp"]
    ) | curation["gdt_atlas_added_fn"]

    tier = np.full(curation.shape[0], "excluded_uncertain", dtype=object)
    tier[curation["gdtai_predicted_gdT"].to_numpy(dtype=bool) & ~curation["gdt_atlas_removed_fp"].to_numpy(dtype=bool)] = "predicted_core"
    tier[curation["fn_silver_unpaired_gd_evidence"].to_numpy(dtype=bool)] = "known_fn_silver"
    tier[curation["fn_gold_sorted_or_paired_or_validation"].to_numpy(dtype=bool)] = "known_fn_gold"
    tier[curation["gdt_atlas_removed_fp"].to_numpy(dtype=bool)] = "removed_fp"
    curation["gdt_atlas_evidence_tier"] = tier
    curation["gdt_atlas_fp_reason"] = fp_reason(curation)
    return curation


def fp_reason(curation: pd.DataFrame) -> np.ndarray:
    reason = np.full(curation.shape[0], "", dtype=object)
    reason[curation["fp_paired_TCRAB_no_gd_evidence"].to_numpy(dtype=bool)] = "paired_TCRAB_no_gd_evidence"
    reason[curation["fp_NK_annotation_no_gd_evidence"].to_numpy(dtype=bool)] = "NK_annotation_no_gd_evidence"
    reason[curation.get("fp_expression_low_cd3_no_gd", False).to_numpy(dtype=bool)] = "expression_low_CD3_no_gd"
    reason[curation.get("fp_expression_nk_like_no_gd", False).to_numpy(dtype=bool)] = "expression_NK_like_no_gd"
    reason[curation.get("fp_expression_regulatory_no_gd", False).to_numpy(dtype=bool)] = "expression_regulatory_no_gd"
    reason[curation.get("fp_expression_myeloid_b_no_gd", False).to_numpy(dtype=bool)] = "expression_myeloid_B_no_gd"
    return reason


def summarize_curation(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label, mask in [
        ("raw_gdtai_predicted", df["gdtai_predicted_gdT"]),
        ("removed_fp", df["gdt_atlas_removed_fp"]),
        ("added_fn_gold", df["fn_gold_sorted_or_paired_or_validation"]),
        ("added_fn_silver", df["fn_silver_unpaired_gd_evidence"]),
        ("curated_primary", df["gdt_atlas_primary_cell"]),
    ]:
        rows.append({"category": label, "n_cells": int(mask.sum()), "fraction_of_total": float(mask.mean())})
    return pd.DataFrame(rows)


def summarize_group(df: pd.DataFrame, group_cols: list[str], path: Path) -> pd.DataFrame:
    out = (
        df.groupby(group_cols, dropna=False)
        .agg(
            total_cells=("gdtai_predicted_gdT", "size"),
            raw_predicted_gdT=("gdtai_predicted_gdT", "sum"),
            removed_fp=("gdt_atlas_removed_fp", "sum"),
            added_fn=("gdt_atlas_added_fn", "sum"),
            curated_primary_gdT=("gdt_atlas_primary_cell", "sum"),
        )
        .reset_index()
    )
    out["curated_fraction_of_total"] = out["curated_primary_gdT"] / out["total_cells"].replace(0, np.nan)
    out["removed_fp_fraction_of_raw_prediction"] = out["removed_fp"] / out["raw_predicted_gdT"].replace(0, np.nan)
    out = out.sort_values(["curated_primary_gdT", "total_cells"], ascending=False)
    out.to_csv(path, index=False)
    return out


def summarize_phenotypes(df: pd.DataFrame, primary: pd.DataFrame) -> dict[str, pd.DataFrame]:
    results: dict[str, pd.DataFrame] = {}
    results["phenotype_counts"] = (
        primary.groupby("gdt_phenotype_state", dropna=False)
        .agg(
            curated_cells=("gdt_atlas_primary_cell", "size"),
            mean_confidence=("gdt_phenotype_confidence", "mean"),
            n_sources=("source_gse_id", "nunique"),
            n_tissues=("tissue_site_gdt_atlas", "nunique"),
        )
        .reset_index()
        .sort_values("curated_cells", ascending=False)
    )
    for name, cols in {
        "phenotype_by_tissue": ["gdt_phenotype_state", "tissue_site_gdt_atlas"],
        "phenotype_by_disease": ["gdt_phenotype_state", "disease_status_gdt_atlas"],
        "phenotype_by_source": ["gdt_phenotype_state", "source_gse_id"],
        "phenotype_by_evidence_tier": ["gdt_phenotype_state", "gdt_atlas_evidence_tier"],
    }.items():
        results[name] = (
            primary.groupby(cols, dropna=False)
            .size()
            .reset_index(name="n_cells")
            .sort_values("n_cells", ascending=False)
        )
    for name, table in results.items():
        table.to_csv(TABLE_DIR / f"{name}.csv", index=False)
    return results


def summarize_scores(df: pd.DataFrame, score_cols: list[str], group_col: str, path: Path) -> pd.DataFrame:
    out = df.groupby(group_col, dropna=False)[score_cols].mean().reset_index()
    counts = df.groupby(group_col, dropna=False).size().reset_index(name="n_cells")
    out = counts.merge(out, on=group_col, how="left").sort_values("n_cells", ascending=False)
    out.to_csv(path, index=False)
    return out


def read_full_umap(h5ad_path: Path) -> np.ndarray | None:
    with h5py.File(h5ad_path, "r") as handle:
        if "obsm" not in handle or "X_umap" not in handle["obsm"]:
            return None
        obj = handle["obsm"]["X_umap"]
        if isinstance(obj, h5py.Dataset):
            return obj[:, :2].astype(np.float32, copy=False)
        if isinstance(obj, h5py.Group) and "data" in obj:
            return obj["data"][:, :2].astype(np.float32, copy=False)
    return None


def sample_indices(indices: np.ndarray, n: int, seed: int) -> np.ndarray:
    indices = np.asarray(indices, dtype=np.int64)
    if indices.size <= n:
        return indices
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(indices, size=n, replace=False))


def plot_curation_flow(summary: pd.DataFrame, path: Path) -> Path:
    order = ["raw_gdtai_predicted", "removed_fp", "added_fn_gold", "added_fn_silver", "curated_primary"]
    plot = summary.set_index("category").loc[order].reset_index()
    fig, ax = plt.subplots(figsize=(8.6, 4.6), constrained_layout=True)
    colors = ["#3b82f6", "#dc2626", "#16a34a", "#84cc16", "#7c3aed"]
    ax.bar(plot["category"], plot["n_cells"], color=colors)
    ax.set_ylabel("cells")
    ax.set_title("Curated gdT atlas curation flow")
    ax.tick_params(axis="x", rotation=25)
    for i, value in enumerate(plot["n_cells"]):
        ax.text(i, value, f"{int(value):,}", ha="center", va="bottom", fontsize=9)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_umap_status(umap_full: np.ndarray | None, df: pd.DataFrame, path: Path, background_n: int, primary_n: int, seed: int) -> Path | None:
    if umap_full is None:
        return None
    bg = sample_indices(df.index.to_numpy(dtype=np.int64), background_n, seed)
    removed = sample_indices(df.index[df["gdt_atlas_removed_fp"]].to_numpy(dtype=np.int64), primary_n // 3, seed + 1)
    added = sample_indices(df.index[df["gdt_atlas_added_fn"]].to_numpy(dtype=np.int64), primary_n // 3, seed + 2)
    primary = sample_indices(df.index[df["gdt_atlas_primary_cell"]].to_numpy(dtype=np.int64), primary_n, seed + 3)
    rows = np.unique(np.concatenate([bg, removed, added, primary]))
    umap = umap_full[rows]
    pos = pd.Series(np.arange(rows.size), index=rows)
    fig, ax = plt.subplots(figsize=(7.2, 6.4), constrained_layout=True)
    ax.scatter(umap[pos.loc[bg], 0], umap[pos.loc[bg], 1], s=1, c="#d1d5db", alpha=0.15, linewidths=0, rasterized=True)
    if removed.size:
        ax.scatter(umap[pos.loc[removed], 0], umap[pos.loc[removed], 1], s=3, c="#dc2626", alpha=0.65, label="removed FP", linewidths=0, rasterized=True)
    if primary.size:
        ax.scatter(umap[pos.loc[primary], 0], umap[pos.loc[primary], 1], s=2, c="#2563eb", alpha=0.45, label="curated gdT", linewidths=0, rasterized=True)
    if added.size:
        ax.scatter(umap[pos.loc[added], 0], umap[pos.loc[added], 1], s=5, c="#16a34a", alpha=0.8, label="added FN", linewidths=0, rasterized=True)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("Curated gdT atlas on 5 million cell T/NK UMAP")
    ax.legend(markerscale=4, frameon=False)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_umap_category(umap_full: np.ndarray | None, df: pd.DataFrame, category: str, path: Path, n: int, seed: int) -> Path | None:
    if umap_full is None:
        return None
    primary_idx = df.index[df["gdt_atlas_primary_cell"]].to_numpy(dtype=np.int64)
    rows = sample_indices(primary_idx, n, seed)
    umap = umap_full[rows]
    plot = df.loc[rows, [category]].copy()
    cats = plot[category].astype(str).fillna("unknown")
    top = cats.value_counts().head(14).index.tolist()
    cats = cats.where(cats.isin(top), "other")
    palette = sns.color_palette("tab20", n_colors=cats.nunique())
    color_map = dict(zip(cats.value_counts().index, palette))
    colors = cats.map(color_map)
    fig, ax = plt.subplots(figsize=(7.4, 6.4), constrained_layout=True)
    ax.scatter(umap[:, 0], umap[:, 1], s=2, c=list(colors), alpha=0.65, linewidths=0, rasterized=True)
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color_map[c], label=c, markersize=5)
        for c in color_map
    ]
    ax.legend(handles=handles, frameon=False, fontsize=7, loc="best", ncol=1)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(category.replace("_", " "))
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_heatmap(table: pd.DataFrame, index: str, columns: str, values: str, path: Path, title: str, top_n: int = 20) -> Path:
    mat = table.pivot_table(index=index, columns=columns, values=values, aggfunc="sum", fill_value=0)
    if mat.empty:
        return path
    mat = mat.loc[mat.sum(axis=1).sort_values(ascending=False).head(top_n).index]
    fig, ax = plt.subplots(figsize=(max(7.0, 0.45 * mat.shape[1] + 3), max(4.2, 0.32 * mat.shape[0] + 1.8)), constrained_layout=True)
    sns.heatmap(np.log10(mat + 1), cmap="viridis", ax=ax, cbar_kws={"label": "log10(cells + 1)"})
    ax.set_title(title)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_score_heatmap(score_table: pd.DataFrame, group_col: str, score_cols: list[str], path: Path, title: str) -> Path:
    mat = score_table.set_index(group_col)[score_cols]
    mat.columns = [c.replace("score_", "") for c in mat.columns]
    fig, ax = plt.subplots(figsize=(max(9.0, 0.42 * mat.shape[1] + 3), max(4.5, 0.34 * mat.shape[0] + 1.8)), constrained_layout=True)
    sns.heatmap(mat, cmap="mako", ax=ax, cbar_kws={"label": "mean log-normalized score"})
    ax.set_title(title)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_bar(table: pd.DataFrame, label_col: str, value_col: str, path: Path, title: str, top_n: int = 25) -> Path:
    plot = table.sort_values(value_col, ascending=False).head(top_n)
    fig, ax = plt.subplots(figsize=(8.4, max(4, 0.3 * plot.shape[0] + 1.2)), constrained_layout=True)
    y = np.arange(plot.shape[0])
    ax.barh(y, plot[value_col], color="#0f766e")
    ax.set_yticks(y, labels=plot[label_col].astype(str))
    ax.invert_yaxis()
    ax.set_title(title)
    ax.set_xlabel(value_col.replace("_", " "))
    for i, value in enumerate(plot[value_col].astype(int)):
        ax.text(value, i, f" {value:,}", va="center", fontsize=8)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_v_gene_panels(umap_full: np.ndarray | None, df: pd.DataFrame, expr: pd.DataFrame, selected_rows: np.ndarray, genes: list[str], path: Path, seed: int) -> Path | None:
    if umap_full is None:
        return None
    row_series = pd.Series(np.arange(selected_rows.size), index=selected_rows)
    primary_rows = df.index[df["gdt_atlas_primary_cell"]].intersection(selected_rows)
    rows = sample_indices(primary_rows.to_numpy(dtype=np.int64), 80_000, seed)
    if rows.size == 0:
        return None
    umap = umap_full[rows]
    expr_pos = row_series.loc[rows].to_numpy(dtype=int)
    present = [gene for gene in genes if gene in expr.columns]
    if not present:
        return None
    ncols = min(3, len(present))
    nrows = int(math.ceil(len(present) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.4 * ncols, 3.8 * nrows), constrained_layout=True)
    axes_arr = np.asarray(axes).reshape(-1)
    for ax, gene in zip(axes_arr, present):
        values = expr.iloc[expr_pos][gene].to_numpy(dtype=float)
        sc = ax.scatter(umap[:, 0], umap[:, 1], c=values, s=2, cmap="magma", linewidths=0, rasterized=True)
        ax.set_title(gene)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    for ax in axes_arr[len(present) :]:
        ax.axis("off")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def render_pdf(skip: bool) -> bool:
    if skip:
        return False
    try:
        subprocess.run(
            [
                "google-chrome",
                "--headless",
                "--disable-gpu",
                "--no-sandbox",
                f"--print-to-pdf={REPORT_PDF}",
                str(REPORT_HTML.resolve()),
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        return REPORT_PDF.exists()
    except Exception as exc:
        logging.warning("PDF render failed: %s", exc)
        return False


def rel(path: Path) -> str:
    return html.escape(str(path.relative_to(PROJECT_ROOT)))


def write_report(
    overall: pd.DataFrame,
    curation_by_source: pd.DataFrame,
    curation_by_tissue: pd.DataFrame,
    phenotype_counts: pd.DataFrame,
    score_by_pheno: pd.DataFrame,
    v_by_pheno: pd.DataFrame,
    figures: list[Path],
    manifest: dict[str, Any],
) -> None:
    lines = [
        "# Curated gdT Phenotype Atlas",
        "",
        f"- Input H5AD: `{manifest['input_h5ad']}`",
        f"- Total cells: `{manifest['n_obs']:,}`",
        f"- Raw gdTAI NK-optimized predictions: `{manifest['raw_predicted_gdT']:,}`",
        f"- Removed likely FP cells: `{manifest['removed_fp']:,}`",
        f"- Added FN cells: `{manifest['added_fn']:,}`",
        f"- Curated primary gdT cells: `{manifest['curated_primary_gdT']:,}`",
        f"- H5AD unchanged: `{manifest['h5ad_unchanged']}`",
        "",
        "## Curation Policy",
        "",
        "- Primary biology uses gdTAI predicted cells after deterministic FP removal plus evidence-backed FN add-back.",
        "- TCR repertoire is not compared because most cells lack TCR-seq; V-gene RNA expression is analyzed instead.",
        "- Silver add-back cells are retained but labeled separately in `gdt_atlas_evidence_tier`.",
        "",
        "## Overall",
        dataframe_to_markdown(overall),
        "",
        "## Outputs",
        f"- Tables: `{TABLE_DIR}`",
        f"- Figures: `{FIGURE_DIR}`",
        f"- HTML: `{REPORT_HTML}`",
        f"- PDF: `{REPORT_PDF}`",
    ]
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    figure_html = "\n".join(
        f"<figure><img src='../../{rel(fig)}'><figcaption>{html.escape(fig.name)}</figcaption></figure>"
        for fig in figures
    )
    css = """
    @page{size:A4 landscape;margin:8mm}
    body{font-family:Arial,Helvetica,sans-serif;margin:24px auto;max-width:1360px;color:#17202a;line-height:1.45}
    h1{font-size:30px;margin:0 0 8px} h2{font-size:21px;margin-top:24px}
    .grid{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:16px 0}
    .metric{border:1px solid #d9e2ec;background:#f8fafc;padding:12px}.value{font-size:24px;font-weight:bold}
    .note{background:#fff7ed;border-left:4px solid #f97316;padding:10px;margin:12px 0}
    table{border-collapse:collapse;width:100%;font-size:10px;table-layout:fixed}
    th,td{border:1px solid #d8dee5;padding:4px 5px;text-align:left;vertical-align:top;overflow-wrap:anywhere}
    th{background:#eef2f6} img{max-width:100%;height:auto;border:1px solid #d8dee5}
    figure{break-inside:avoid;margin:18px 0}
    """
    html_text = f"""<!doctype html><html><head><meta charset='utf-8'>
<title>Curated gdT Phenotype Atlas</title><style>{css}</style></head><body>
<h1>Curated gdT Phenotype Atlas</h1>
<p>This report curates gdTAI v2 NK-optimized predictions by removing likely false positives and adding evidence-backed false negatives, then focuses on gdT phenotypes and V-gene RNA expression.</p>
<div class='grid'>
  <div class='metric'><div>Raw gdTAI predictions</div><div class='value'>{manifest['raw_predicted_gdT']:,}</div></div>
  <div class='metric'><div>Removed FP</div><div class='value'>{manifest['removed_fp']:,}</div></div>
  <div class='metric'><div>Added FN</div><div class='value'>{manifest['added_fn']:,}</div></div>
  <div class='metric'><div>Curated gdT</div><div class='value'>{manifest['curated_primary_gdT']:,}</div></div>
</div>
<section class='note'>TCR repertoire comparisons are intentionally omitted because most cells lack TCR-seq data. The report analyzes TRD/TRG V-gene RNA expression instead.</section>
<h2>Figures</h2>{figure_html}
<h2>Curation Summary</h2>{dataframe_to_html(overall, 20)}
<h2>Curated gdT By Source</h2>{dataframe_to_html(curation_by_source, 35)}
<h2>Curated gdT By Tissue</h2>{dataframe_to_html(curation_by_tissue, 35)}
<h2>Phenotype Counts</h2>{dataframe_to_html(phenotype_counts, 40)}
<h2>Marker Program Scores By Phenotype</h2>{dataframe_to_html(score_by_pheno, 40)}
<h2>V-Gene RNA Expression By Phenotype</h2>{dataframe_to_html(v_by_pheno, 40)}
</body></html>"""
    REPORT_HTML.write_text(html_text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    setup_logging()
    marker_config = json.loads(args.marker_config.read_text(encoding="utf-8"))
    metadata_rules = json.loads(args.rules.read_text(encoding="utf-8"))

    logging.info("Loading obs metadata")
    df, info = load_obs(args.input_h5ad)
    pred_arrays = load_prediction_arrays(args.prediction_npz, int(info["n_obs"]))
    pred = pred_arrays[PREDICTION_KEY]
    curation = build_curation(df, pred, args.fn_table)

    logging.info("Applying harmonized metadata")
    lookups = load_metadata_lookup(args.metadata)
    df = apply_metadata_lookup(df, lookups)
    df = add_harmonized_metadata(df, metadata_rules)

    logging.info("Extracting marker expression for predicted and add-back candidates")
    marker_genes = ordered_unique(
        [gene for genes in marker_config["programs"].values() for gene in genes]
        + marker_config["v_gene_expression"]
    )
    analysis_mask = pred | curation["gdt_atlas_added_fn"].to_numpy(dtype=bool)
    analysis_rows = np.flatnonzero(analysis_mask)
    expr, total_counts, present_genes, missing_genes = load_selected_expression(
        args.input_h5ad,
        analysis_rows,
        marker_genes,
        args.expression_row_block,
    )
    pd.DataFrame({"present_gene": present_genes}).to_csv(TABLE_DIR / "marker_present_genes.csv", index=False)
    pd.DataFrame({"missing_gene": missing_genes}).to_csv(TABLE_DIR / "marker_missing_genes.csv", index=False)

    scores, marker_manifest = score_programs(expr, marker_config)
    marker_manifest.to_csv(TABLE_DIR / "marker_program_manifest.csv", index=False)

    full_scores = pd.DataFrame(index=np.arange(df.shape[0]))
    for column in scores.columns:
        full_scores[column] = np.nan
    full_scores.loc[analysis_rows, scores.columns] = scores.to_numpy()
    full_scores = full_scores.fillna(0.0)
    curation = apply_expression_fp_flags(curation, full_scores)
    primary_mask = curation["gdt_atlas_primary_cell"].to_numpy(dtype=bool)
    analysis_primary_mask = primary_mask[analysis_rows]
    phenotype = assign_phenotypes(scores, analysis_primary_mask)
    phenotype_full = pd.DataFrame(
        {
            "gdt_phenotype_state": np.full(df.shape[0], "not_curated_gdT", dtype=object),
            "gdt_phenotype_confidence": np.zeros(df.shape[0], dtype=np.float32),
        }
    )
    phenotype_full.loc[analysis_rows, "gdt_phenotype_state"] = phenotype["gdt_phenotype_state"].to_numpy(dtype=object)
    phenotype_full.loc[analysis_rows, "gdt_phenotype_confidence"] = phenotype[
        "gdt_phenotype_confidence"
    ].to_numpy(dtype=np.float32)

    logging.info("Combining curation and phenotype metadata")
    for column in curation.columns:
        df[column] = curation[column].to_numpy()
    df["gdt_phenotype_state"] = phenotype_full["gdt_phenotype_state"].to_numpy(dtype=object)
    df["gdt_phenotype_confidence"] = phenotype_full["gdt_phenotype_confidence"].to_numpy(dtype=np.float32)
    df["total_counts_selected_rows"] = np.nan
    df.loc[analysis_rows, "total_counts_selected_rows"] = total_counts
    for column in scores.columns:
        df[column] = full_scores[column].to_numpy(dtype=np.float32)

    expr_index = pd.Series(np.arange(analysis_rows.size), index=analysis_rows)
    v_genes = [gene for gene in marker_config["v_gene_expression"] if gene in expr.columns]
    for gene in v_genes:
        df[f"expr_{gene}"] = np.nan
        df.loc[analysis_rows, f"expr_{gene}"] = expr[gene].to_numpy(dtype=np.float32)

    logging.info("Writing tables")
    curation_summary = summarize_curation(df)
    curation_summary.to_csv(TABLE_DIR / "curation_summary.csv", index=False)
    curation_by_source = summarize_group(df, ["source_gse_id"], TABLE_DIR / "curation_by_source.csv")
    curation_by_tissue = summarize_group(df, ["tissue_site_gdt_atlas"], TABLE_DIR / "curation_by_tissue.csv")
    curation_by_disease = summarize_group(df, ["disease_status_gdt_atlas"], TABLE_DIR / "curation_by_disease.csv")
    summarize_group(df, ["source_gse_id", "tissue_site_gdt_atlas"], TABLE_DIR / "curation_by_source_tissue.csv")

    removed_cols = [
        "obs_name",
        "source_gse_id",
        "sample_id",
        "library_id",
        "simple_annotation_plus6",
        "tissue_site_gdt_atlas",
        "disease_status_gdt_atlas",
        "gdt_atlas_fp_reason",
        "strong_gd_evidence",
        "has_TRA_TRB_paired",
        "has_TRG_TRD_paired",
        "has_any_ab_tcr",
        "has_any_gd_tcr",
        "phase4_trd_score",
        "phase4_trab_score",
        "phase4_trd_minus_trab",
        "score_core_t_gdt",
        "score_nk_like",
        "score_regulatory_penalty",
        "score_myeloid_b_contamination",
    ]
    df.loc[df["gdt_atlas_removed_fp"], removed_cols].to_csv(TABLE_DIR / "removed_fp_cells.csv.gz", index=False)
    df.loc[df["gdt_atlas_added_fn"], removed_cols + ["gdt_atlas_evidence_tier"]].to_csv(
        TABLE_DIR / "added_fn_cells.csv.gz", index=False
    )
    primary_cols = [
        "obs_name",
        "source_gse_id",
        "sample_id",
        "library_id",
        "donor_id",
        "simple_annotation_plus6",
        "tissue_site_gdt_atlas",
        "disease_status_gdt_atlas",
        "condition_detail_gdt_atlas",
        "age_group_gdt_atlas",
        "sex_gdt_atlas",
        "gdt_atlas_evidence_tier",
        "gdt_phenotype_state",
        "gdt_phenotype_confidence",
        "phase4_trd_score",
        "phase4_trab_score",
        "phase4_trd_minus_trab",
    ]
    score_cols = [col for col in df.columns if col.startswith("score_")]
    v_cols = [f"expr_{gene}" for gene in v_genes]
    df.loc[df["gdt_atlas_primary_cell"], primary_cols + score_cols + v_cols].to_csv(
        TABLE_DIR / "curated_gdt_cell_metadata.csv.gz", index=False
    )

    primary = df.loc[df["gdt_atlas_primary_cell"]].copy()
    phenotype_tables = summarize_phenotypes(df, primary)
    score_by_pheno = summarize_scores(
        primary,
        score_cols,
        "gdt_phenotype_state",
        TABLE_DIR / "marker_program_scores_by_phenotype.csv",
    )
    v_by_pheno = summarize_scores(
        primary,
        v_cols,
        "gdt_phenotype_state",
        TABLE_DIR / "v_gene_expression_by_phenotype.csv",
    )
    summarize_scores(primary, v_cols, "tissue_site_gdt_atlas", TABLE_DIR / "v_gene_expression_by_tissue.csv")
    summarize_scores(primary, score_cols, "tissue_site_gdt_atlas", TABLE_DIR / "marker_program_scores_by_tissue.csv")

    sample_keys = ["source_gse_id", "sample_id", "library_id", "tissue_site_gdt_atlas", "disease_status_gdt_atlas"]
    sample_abund = (
        df.groupby(sample_keys, dropna=False)
        .agg(total_cells=("gdtai_predicted_gdT", "size"), curated_gdT=("gdt_atlas_primary_cell", "sum"))
        .reset_index()
    )
    sample_abund["curated_gdT_fraction"] = sample_abund["curated_gdT"] / sample_abund["total_cells"].replace(0, np.nan)
    sample_abund.to_csv(TABLE_DIR / "sample_level_curated_gdt_abundance.csv", index=False)

    logging.info("Plotting figures")
    figures: list[Path] = []
    umap_full = read_full_umap(args.input_h5ad)
    figures.append(plot_curation_flow(curation_summary, FIGURE_DIR / "curation_flow.png"))
    umap_status = plot_umap_status(
        umap_full, df, FIGURE_DIR / "whole_atlas_umap_curated_status.png", args.umap_background, args.umap_primary, args.seed
    )
    if umap_status is not None:
        figures.append(umap_status)
    for category, filename in [
        ("gdt_phenotype_state", "curated_gdt_umap_by_phenotype.png"),
        ("gdt_atlas_evidence_tier", "curated_gdt_umap_by_evidence_tier.png"),
        ("tissue_site_gdt_atlas", "curated_gdt_umap_by_tissue.png"),
        ("disease_status_gdt_atlas", "curated_gdt_umap_by_disease.png"),
    ]:
        out = plot_umap_category(umap_full, df, category, FIGURE_DIR / filename, args.umap_primary, args.seed + len(figures))
        if out is not None:
            figures.append(out)
    figures.append(plot_bar(curation_by_source, "source_gse_id", "curated_primary_gdT", FIGURE_DIR / "curated_gdt_by_source.png", "Curated gdT cells by source"))
    figures.append(plot_bar(curation_by_tissue, "tissue_site_gdt_atlas", "curated_primary_gdT", FIGURE_DIR / "curated_gdt_by_tissue.png", "Curated gdT cells by tissue"))
    figures.append(plot_bar(curation_by_disease, "disease_status_gdt_atlas", "curated_primary_gdT", FIGURE_DIR / "curated_gdt_by_disease.png", "Curated gdT cells by disease"))
    figures.append(plot_heatmap(phenotype_tables["phenotype_by_tissue"], "gdt_phenotype_state", "tissue_site_gdt_atlas", "n_cells", FIGURE_DIR / "phenotype_by_tissue_heatmap.png", "Curated gdT phenotypes by tissue"))
    figures.append(plot_heatmap(phenotype_tables["phenotype_by_disease"], "gdt_phenotype_state", "disease_status_gdt_atlas", "n_cells", FIGURE_DIR / "phenotype_by_disease_heatmap.png", "Curated gdT phenotypes by disease"))
    figures.append(plot_score_heatmap(score_by_pheno, "gdt_phenotype_state", score_cols, FIGURE_DIR / "marker_program_scores_by_phenotype.png", "Marker program scores by phenotype"))
    figures.append(plot_score_heatmap(v_by_pheno, "gdt_phenotype_state", v_cols, FIGURE_DIR / "v_gene_expression_by_phenotype.png", "TRD/TRG V-gene RNA expression by phenotype"))
    v_panel = plot_v_gene_panels(
        umap_full,
        df,
        expr,
        analysis_rows,
        marker_config["primary_v_gene_panels"],
        FIGURE_DIR / "v_gene_expression_umap_panels.png",
        args.seed,
    )
    if v_panel is not None:
        figures.append(v_panel)

    stat_after = args.input_h5ad.stat()
    manifest = {
        **info,
        "prediction_npz": str(args.prediction_npz),
        "prediction_key": PREDICTION_KEY,
        "raw_predicted_gdT": int(pred.sum()),
        "removed_fp": int(df["gdt_atlas_removed_fp"].sum()),
        "added_fn": int(df["gdt_atlas_added_fn"].sum()),
        "added_fn_gold": int(df["fn_gold_sorted_or_paired_or_validation"].sum()),
        "added_fn_silver": int(df["fn_silver_unpaired_gd_evidence"].sum()),
        "curated_primary_gdT": int(df["gdt_atlas_primary_cell"].sum()),
        "size_after": stat_after.st_size,
        "mtime_ns_after": stat_after.st_mtime_ns,
        "h5ad_unchanged": bool(info["size_before"] == stat_after.st_size and info["mtime_ns_before"] == stat_after.st_mtime_ns),
        "html": str(REPORT_HTML),
        "pdf": str(REPORT_PDF),
        "tables": str(TABLE_DIR),
        "figures": str(FIGURE_DIR),
    }
    write_report(
        curation_summary,
        curation_by_source,
        curation_by_tissue,
        phenotype_tables["phenotype_counts"],
        score_by_pheno,
        v_by_pheno,
        figures,
        manifest,
    )
    pdf_ok = render_pdf(args.no_pdf)
    manifest["pdf_rendered"] = bool(pdf_ok)
    MANIFEST_JSON.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    logging.info("Wrote report: %s", REPORT_HTML)
    logging.info("H5AD unchanged: %s", manifest["h5ad_unchanged"])


if __name__ == "__main__":
    main()
