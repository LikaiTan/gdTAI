#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Rigorous gdT-only atlas reanalysis from gdTAI predictions.

This workflow is read-only with respect to the source H5AD. It builds a
curated gdT biology subset from the promoted gdTAI prediction output, excludes
silver-only FN rescue cells, re-runs gdT-only integration/UMAP/clustering, and
renders a statistics-focused HTML/PDF report.
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

import anndata as ad
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

from run_gdt_atlas_curated_phenotype_analysis import (
    load_selected_expression,
    ordered_unique,
    sample_indices,
    score_programs,
)
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


TABLE_DIR = PROJECT_ROOT / "Integrated_dataset/tables/gdT_atlas_rigorous"
FIGURE_DIR = PROJECT_ROOT / "Integrated_dataset/figures/gdT_atlas_rigorous"
LOG_DIR = PROJECT_ROOT / "Integrated_dataset/logs/gdT_atlas_rigorous"
MODEL_DIR = PROJECT_ROOT / "Integrated_dataset/models/gdT_atlas_rigorous"
STATIC_DIR = PROJECT_ROOT / "gdT_atlas/rigorous"
DERIVED_H5AD = PROJECT_ROOT / "Integrated_dataset/gdT_atlas/gdt_curated_integrated.h5ad"
REPORT_HTML = STATIC_DIR / "index.html"
REPORT_PDF = STATIC_DIR / "gdT_atlas_rigorous_report.pdf"
SUMMARY_MD = LOG_DIR / "gdt_atlas_rigorous_summary.md"
MANIFEST_JSON = LOG_DIR / "gdt_atlas_rigorous_manifest.json"
RUN_LOG = LOG_DIR / "run.log"

DEFAULT_MARKERS = PROJECT_ROOT / "configs" / "gdt_atlas" / "literature_marker_programs.json"
DEFAULT_FN_TABLE = (
    PROJECT_ROOT
    / "Integrated_dataset/tables/gdT_prediction/gdtai_annotation_specific_tp_fn_audit/"
    / "validation_positive_tp_fn_cell_table.csv.gz"
)

PREDICTION_KEY = "annotation_specific_pred"
EXPECTED_N_OBS = 5_128_904
EXPECTED_PROMOTED_PREDICTIONS = 359_857
INVALID_STRINGS = {"", "nan", "none", "na", "n/a", "<na>", "null", "unknown"}

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
    "TRA_j",
    "TRB_j",
    "TRG_j",
    "TRD_j",
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
FUNCTIONAL_SCORE_COLUMNS = [
    "score_cytotoxic_effector",
    "score_naive_memory",
    "score_tissue_resident",
    "score_il17_inflammatory",
    "score_ifn_response",
    "score_proliferating",
    "score_checkpoint_activation",
]
WARNING_SCORE_COLUMNS = [
    "score_nk_like_warning",
    "score_cd4_treg_death_penalty",
    "score_myeloid_b_contamination",
    "score_stress_low_quality",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run rigorous gdT-only atlas reanalysis.")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--prediction-npz", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--metadata-rules", type=Path, default=DEFAULT_RULES)
    parser.add_argument("--metadata", type=Path, action="append", default=DEFAULT_METADATA)
    parser.add_argument("--marker-config", type=Path, default=DEFAULT_MARKERS)
    parser.add_argument("--fn-table", type=Path, default=DEFAULT_FN_TABLE)
    parser.add_argument("--output-h5ad", type=Path, default=DERIVED_H5AD)
    parser.add_argument("--expression-row-block", type=int, default=5000)
    parser.add_argument("--extract-row-block", type=int, default=5000)
    parser.add_argument("--n-hvgs", type=int, default=3000)
    parser.add_argument("--scvi-epochs", type=int, default=80)
    parser.add_argument("--seed", type=int, default=20260618)
    parser.add_argument("--no-pdf", action="store_true")
    parser.add_argument("--skip-scvi", action="store_true")
    parser.add_argument("--resume-from-derived", action="store_true")
    return parser.parse_args()


def ensure_dirs() -> None:
    for path in [TABLE_DIR, FIGURE_DIR, LOG_DIR, MODEL_DIR, STATIC_DIR, DERIVED_H5AD.parent]:
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
        data: dict[str, Any] = {
            "obs_name": read_string_dataset(handle["obs"]["_index"]),
            "obs_row": np.arange(n_obs, dtype=np.int64),
        }
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
            df[column] = pd.to_numeric(df[column], errors="coerce").fillna(0).astype(np.float32)
        elif df[column].dtype == object or str(df[column].dtype).startswith("string"):
            df[column] = clean_series(df[column]).astype(str)
    return df, {
        "input_h5ad": str(path),
        "n_obs": n_obs,
        "has_source_umap": bool(has_umap),
        "has_source_scvi": bool(has_scvi),
        "source_size_before": stat.st_size,
        "source_mtime_ns_before": stat.st_mtime_ns,
    }


def load_prediction_arrays(path: Path, n_obs: int) -> dict[str, np.ndarray]:
    with np.load(path) as arrays:
        loaded = {key: np.asarray(arrays[key]) for key in arrays.files}
    if PREDICTION_KEY not in loaded:
        raise KeyError(f"Missing `{PREDICTION_KEY}` in {path}")
    pred = loaded[PREDICTION_KEY].astype(bool)
    if pred.shape[0] != n_obs:
        raise RuntimeError(f"Prediction length {pred.shape[0]:,} does not match n_obs {n_obs:,}")
    if n_obs != EXPECTED_N_OBS:
        raise RuntimeError(f"Unexpected source n_obs {n_obs:,}; expected {EXPECTED_N_OBS:,}")
    if int(pred.sum()) != EXPECTED_PROMOTED_PREDICTIONS:
        raise RuntimeError(
            f"Unexpected promoted prediction count {int(pred.sum()):,}; "
            f"expected {EXPECTED_PROMOTED_PREDICTIONS:,}"
        )
    loaded[PREDICTION_KEY] = pred
    for key in ["current_pred", "original_pred"]:
        if key in loaded:
            loaded[key] = loaded[key].astype(bool)
    return loaded


def nonempty_text(series: pd.Series) -> np.ndarray:
    cleaned = clean_series(series).astype(str).str.lower()
    return ~cleaned.isin(INVALID_STRINGS).to_numpy(dtype=bool)


def load_gold_fn_mask(path: Path, n_obs: int) -> tuple[np.ndarray, pd.DataFrame]:
    mask = np.zeros(n_obs, dtype=bool)
    if not path.exists():
        logging.warning("Gold FN audit table is missing: %s", path)
        return mask, pd.DataFrame()
    cols = pd.read_csv(path, nrows=0).columns.tolist()
    usecols = [c for c in ["obs_row", "outcome", "validation_component", "source_gse_id"] if c in cols]
    if "obs_row" not in usecols or "outcome" not in usecols:
        logging.warning("Gold FN audit table lacks obs_row/outcome columns: %s", path)
        return mask, pd.DataFrame()
    table = pd.read_csv(path, usecols=usecols)
    fn = table["outcome"].astype(str).eq("FN")
    if "validation_component" in table:
        fn &= ~table["validation_component"].astype(str).str.contains("silver", case=False, na=False)
    rows = pd.to_numeric(table.loc[fn, "obs_row"], errors="coerce").dropna().astype(np.int64)
    rows = rows[(rows >= 0) & (rows < n_obs)]
    mask[rows.to_numpy()] = True
    summary_cols = [c for c in ["outcome", "validation_component", "source_gse_id"] if c in table.columns]
    summary = table.groupby(summary_cols, dropna=False).size().reset_index(name="n_cells") if summary_cols else table
    return mask, summary


def build_initial_curation(df: pd.DataFrame, pred: np.ndarray, gold_fn: np.ndarray) -> pd.DataFrame:
    has_trg = nonempty_text(df["TRG_cdr3"])
    has_trd = nonempty_text(df["TRD_cdr3"])
    has_ab = df["has_any_ab_tcr"].astype(bool).to_numpy()
    has_paired_ab = df["has_TRA_TRB_paired"].astype(bool).to_numpy()
    has_paired_gd = df["has_TRG_TRD_paired"].astype(bool).to_numpy()
    has_any_gd = df["has_any_gd_tcr"].astype(bool).to_numpy()
    sorted_gold = df["Sorted_gdT"].astype(bool).to_numpy()
    trd_score = df["phase4_trd_score"].to_numpy(dtype=float)
    trd_minus_trab = df["phase4_trd_minus_trab"].to_numpy(dtype=float)
    annotation = df["simple_annotation_plus6"].astype(str).str.upper()

    strong_gd_evidence = (
        sorted_gold
        | (has_paired_gd & ~has_ab)
        | has_any_gd
        | has_trg
        | has_trd
        | ((trd_score > 0.35) & (trd_minus_trab > 0.2))
    )
    is_nk = annotation.str.contains("NK", regex=False).to_numpy(dtype=bool)
    is_cd4_or_treg = annotation.str.contains("CD4", regex=False).to_numpy(dtype=bool) | annotation.str.contains(
        "TREG", regex=False
    ).to_numpy(dtype=bool)

    silver_excluded = (~pred) & (~gold_fn) & (~has_ab) & (~sorted_gold) & (has_any_gd | has_trg | has_trd)
    curation = pd.DataFrame(
        {
            "gdtai_promoted_prediction": pred,
            "strong_gd_evidence": strong_gd_evidence,
            "gold_fn_added": gold_fn,
            "silver_fn_excluded": silver_excluded,
            "fp_paired_TCRAB_no_gd_evidence": pred & has_paired_ab & ~strong_gd_evidence,
            "fp_NK_annotation_no_gd_evidence": pred & is_nk & ~strong_gd_evidence,
            "fp_CD4_Treg_annotation_no_gd_evidence": pred & is_cd4_or_treg & ~strong_gd_evidence,
        }
    )
    curation["removed_fp_pre_expression"] = curation[
        [
            "fp_paired_TCRAB_no_gd_evidence",
            "fp_NK_annotation_no_gd_evidence",
            "fp_CD4_Treg_annotation_no_gd_evidence",
        ]
    ].any(axis=1)
    return curation


def add_expression_fp_flags(
    curation: pd.DataFrame,
    candidate_rows: np.ndarray,
    scores: pd.DataFrame,
) -> pd.DataFrame:
    pred = curation.loc[candidate_rows, "gdtai_promoted_prediction"].to_numpy(dtype=bool)
    no_gd = ~curation.loc[candidate_rows, "strong_gd_evidence"].to_numpy(dtype=bool)
    core = scores.get("score_core_t_gdt", pd.Series(0, index=scores.index)).to_numpy(dtype=float)
    nk = scores.get("score_nk_like_warning", pd.Series(0, index=scores.index)).to_numpy(dtype=float)
    treg = scores.get("score_cd4_treg_death_penalty", pd.Series(0, index=scores.index)).to_numpy(dtype=float)
    contam = scores.get("score_myeloid_b_contamination", pd.Series(0, index=scores.index)).to_numpy(dtype=float)

    flags = {
        "fp_expression_low_CD3_no_gd": pred & no_gd & (core < 0.15),
        "fp_expression_NK_like_no_gd": pred & no_gd & (nk > 1.0) & (core < 0.55),
        "fp_expression_CD4_Treg_no_gd": pred & no_gd & (treg > 0.65) & (core < 0.65),
        "fp_expression_myeloid_B_no_gd": pred & no_gd & (contam > 0.6),
    }
    for column, values in flags.items():
        curation[column] = False
        curation.loc[candidate_rows, column] = values
    fp_cols = ["removed_fp_pre_expression", *flags.keys()]
    curation["removed_fp"] = curation[fp_cols].any(axis=1)
    curation["curated_primary_gdT"] = (
        curation["gdtai_promoted_prediction"] & ~curation["removed_fp"]
    ) | curation["gold_fn_added"]

    tier = np.full(curation.shape[0], "excluded_other", dtype=object)
    tier[curation["gdtai_promoted_prediction"].to_numpy(dtype=bool)] = "predicted_core"
    tier[curation["removed_fp"].to_numpy(dtype=bool)] = "removed_likely_FP"
    tier[curation["silver_fn_excluded"].to_numpy(dtype=bool)] = "silver_FN_excluded"
    tier[curation["gold_fn_added"].to_numpy(dtype=bool)] = "gold_FN_added"
    curation["gdt_atlas_evidence_tier"] = tier
    curation["fp_reason"] = fp_reason(curation)
    return curation


def fp_reason(curation: pd.DataFrame) -> np.ndarray:
    reason = np.full(curation.shape[0], "", dtype=object)
    for column, label in [
        ("fp_paired_TCRAB_no_gd_evidence", "paired_TCRAB_no_gd_evidence"),
        ("fp_NK_annotation_no_gd_evidence", "NK_annotation_no_gd_evidence"),
        ("fp_CD4_Treg_annotation_no_gd_evidence", "CD4_Treg_annotation_no_gd_evidence"),
        ("fp_expression_low_CD3_no_gd", "expression_low_CD3_no_gd"),
        ("fp_expression_NK_like_no_gd", "expression_NK_like_no_gd"),
        ("fp_expression_CD4_Treg_no_gd", "expression_CD4_Treg_no_gd"),
        ("fp_expression_myeloid_B_no_gd", "expression_myeloid_B_no_gd"),
    ]:
        if column in curation:
            reason[curation[column].to_numpy(dtype=bool)] = label
    return reason


def read_var_dataframe(path: Path) -> pd.DataFrame:
    adata = ad.read_h5ad(path, backed="r")
    var = adata.var.copy()
    adata.file.close()
    return var


def extract_csr_rows(path: Path, selected_rows: np.ndarray, block_rows: int) -> sparse.csr_matrix:
    selected_rows = np.asarray(selected_rows, dtype=np.int64)
    if selected_rows.size == 0:
        raise RuntimeError("No rows selected for gdT-only H5AD")
    with h5py.File(path, "r") as handle:
        x = handle["X"]
        n_vars = int(handle["var"]["_index"].shape[0])
        indptr = x["indptr"][:]
        starts = indptr[selected_rows]
        ends = indptr[selected_rows + 1]
        row_nnz = (ends - starts).astype(np.int64, copy=False)
        out_indptr = np.empty(selected_rows.size + 1, dtype=np.int64)
        out_indptr[0] = 0
        np.cumsum(row_nnz, out=out_indptr[1:])
        total_nnz = int(out_indptr[-1])
        logging.info("Extracting gdT-only matrix: %s cells, %s nnz", f"{selected_rows.size:,}", f"{total_nnz:,}")
        out_data = np.empty(total_nnz, dtype=x["data"].dtype)
        out_indices = np.empty(total_nnz, dtype=np.int32)
        data_ds = x["data"]
        indices_ds = x["indices"]

        cursor = 0
        n_obs = int(handle["obs"]["_index"].shape[0])
        for block_start in range(0, n_obs, block_rows):
            block_end = min(n_obs, block_start + block_rows)
            lo = np.searchsorted(selected_rows, block_start, side="left")
            hi = np.searchsorted(selected_rows, block_end, side="left")
            if hi <= lo:
                continue
            data_start = int(indptr[block_start])
            data_end = int(indptr[block_end])
            block_data = data_ds[data_start:data_end]
            block_indices = indices_ds[data_start:data_end]
            for out_i in range(lo, hi):
                src_start = int(indptr[selected_rows[out_i]] - data_start)
                src_end = int(indptr[selected_rows[out_i] + 1] - data_start)
                dst_start = int(out_indptr[out_i])
                dst_end = int(out_indptr[out_i + 1])
                out_data[dst_start:dst_end] = block_data[src_start:src_end]
                out_indices[dst_start:dst_end] = block_indices[src_start:src_end].astype(np.int32, copy=False)
            cursor = hi
            if block_start % (block_rows * 50) == 0:
                logging.info("Matrix extraction source rows %s/%s", f"{block_start:,}", f"{n_obs:,}")
        if cursor != selected_rows.size:
            logging.info("Matrix extraction completed with %s selected rows copied", f"{cursor:,}")
    return sparse.csr_matrix((out_data, out_indices, out_indptr), shape=(selected_rows.size, n_vars))


def compress_obs(obs: pd.DataFrame) -> pd.DataFrame:
    out = obs.copy()
    for column in out.columns:
        if out[column].dtype == object:
            nunique = out[column].nunique(dropna=False)
            if nunique < min(10000, max(100, out.shape[0] // 5)):
                out[column] = out[column].astype("category")
    return out


def robust_z(scores: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    z = pd.DataFrame(index=scores.index)
    for column in columns:
        values = scores[column].to_numpy(dtype=float)
        median = float(np.nanmedian(values))
        mad = float(np.nanmedian(np.abs(values - median)))
        denom = 1.4826 * mad if mad > 1e-6 else float(np.nanstd(values) + 1e-6)
        z[column.replace("score_", "z_")] = (values - median) / denom
    return z


def assign_phenotypes(scores: pd.DataFrame) -> pd.DataFrame:
    functional = [c for c in FUNCTIONAL_SCORE_COLUMNS if c in scores.columns]
    z = robust_z(scores, functional + [c for c in WARNING_SCORE_COLUMNS if c in scores.columns])
    call_names = [c.replace("score_", "").replace("z_", "") for c in functional]
    arr = z[[c.replace("score_", "z_") for c in functional]].to_numpy(dtype=float)
    best_idx = np.argmax(arr, axis=1)
    best = arr[np.arange(arr.shape[0]), best_idx]
    state = np.asarray([call_names[i] for i in best_idx], dtype=object)
    state[best < 0.25] = "gdT_other"
    if "z_proliferating" in z:
        state[z["z_proliferating"].to_numpy(dtype=float) > 1.25] = "proliferating"
    if "z_il17_inflammatory" in z:
        state[z["z_il17_inflammatory"].to_numpy(dtype=float) > 1.25] = "il17_inflammatory"
    if "z_nk_like_warning" in z and "z_core_t_gdt" in z:
        state[(z["z_nk_like_warning"].to_numpy(dtype=float) > 1.4) & (z["z_core_t_gdt"].to_numpy(dtype=float) < 0)] = (
            "NK_like_warning"
        )
    elif "z_nk_like_warning" in z:
        state[z["z_nk_like_warning"].to_numpy(dtype=float) > 1.8] = "NK_like_warning"
    if "z_cd4_treg_death_penalty" in z:
        state[z["z_cd4_treg_death_penalty"].to_numpy(dtype=float) > 1.5] = "CD4_Treg_warning"
    confidence = np.sort(arr, axis=1)[:, -1] - np.sort(arr, axis=1)[:, -2] if arr.shape[1] > 1 else best
    out = pd.DataFrame({"gdt_phenotype_state": state, "gdt_phenotype_confidence": confidence}, index=scores.index)
    return pd.concat([out, z], axis=1)


def assign_v_gene_state(expr: pd.DataFrame) -> pd.Series:
    candidates = [g for g in ["TRDV1", "TRDV2", "TRDV3"] if g in expr.columns]
    if not candidates:
        return pd.Series("TRDV_not_detected", index=expr.index)
    values = expr[candidates].to_numpy(dtype=float)
    best_idx = np.argmax(values, axis=1)
    best = values[np.arange(values.shape[0]), best_idx]
    ordered = np.sort(values, axis=1)
    gap = best - (ordered[:, -2] if values.shape[1] > 1 else 0)
    states = np.asarray([f"{candidates[i]}_dominant" for i in best_idx], dtype=object)
    states[(best < 0.25) | (gap < 0.1)] = "TRDV_not_dominant"
    return pd.Series(states, index=expr.index)


def make_sample_unit(df: pd.DataFrame) -> pd.Series:
    lib = clean_series(df["library_id"]).astype(str)
    sample = clean_series(df["sample_id"]).astype(str)
    sampleid = clean_series(df["sampleid"]).astype(str)
    unit = lib.mask(lib.eq(""), sample).mask(lambda s: s.eq(""), sampleid)
    unit = unit.mask(unit.eq(""), df["obs_name"].astype(str))
    return df["source_gse_id"].astype(str) + "||" + unit.astype(str)


def select_hvgs_no_tcr(adata: ad.AnnData, n_hvgs: int) -> list[str]:
    batch_key = "source_gse_id" if "source_gse_id" in adata.obs else None
    top_request = min(8000, adata.n_vars)
    logging.info("Selecting HVGs with TCR genes excluded from final HVG set")
    try:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=top_request,
            flavor="seurat_v3",
            layer="counts",
            batch_key=batch_key,
            subset=False,
        )
    except Exception as exc:
        logging.warning("seurat_v3 HVG selection failed; using sparse count-dispersion fallback: %s", exc)
        x = adata.layers["counts"]
        means = np.asarray(x.mean(axis=0)).ravel().astype(float)
        sq_means = np.asarray(x.power(2).mean(axis=0)).ravel().astype(float)
        variances = np.maximum(sq_means - means**2, 0.0)
        dispersions = variances / np.maximum(means, 1e-8)
        dispersions[means <= 0] = -np.inf
        order = np.argsort(-dispersions)
        ranks = np.full(adata.n_vars, np.nan, dtype=float)
        ranks[order[:top_request]] = np.arange(min(top_request, adata.n_vars), dtype=float)
        adata.var["means"] = means
        adata.var["variances"] = variances
        adata.var["dispersions"] = dispersions
        adata.var["highly_variable_rank"] = ranks
        adata.var["highly_variable"] = np.isfinite(ranks)
        adata.uns["hvg"] = {"flavor": "sparse_count_dispersion_fallback"}
    var = adata.var.copy()
    var["gene_name"] = var.index.astype(str)
    excluded = var["gene_name"].str.startswith(("TRA", "TRB", "TRG", "TRD"))
    ranked = var.loc[var["highly_variable"].fillna(False) & ~excluded].sort_values("highly_variable_rank")
    selected = ranked.index[:n_hvgs].tolist()
    if len(selected) < n_hvgs:
        extra = var.loc[~excluded & ~var.index.isin(selected)].sort_values(
            "means" if "means" in var.columns else "gene_name", ascending=False
        )
        selected.extend(extra.index[: n_hvgs - len(selected)].tolist())
    adata.var["highly_variable_no_tcr"] = adata.var.index.isin(selected)
    adata.var["excluded_from_hvg_tcr_gene"] = excluded.to_numpy(dtype=bool)
    pd.DataFrame({"gene": selected}).to_csv(TABLE_DIR / "hvg_no_tcr_genes.csv", index=False)
    pd.DataFrame({"gene": var.loc[excluded].index}).to_csv(TABLE_DIR / "hvg_excluded_tcr_genes.csv", index=False)
    return selected


def run_integration(adata: ad.AnnData, n_hvgs: int, epochs: int, seed: int, skip_scvi: bool) -> dict[str, Any]:
    hvgs = select_hvgs_no_tcr(adata, n_hvgs)
    hvg = adata[:, hvgs].copy()
    integration_info: dict[str, Any] = {"n_hvgs": len(hvgs), "integration_method": "scVI"}
    if skip_scvi:
        raise RuntimeError("--skip-scvi was supplied; use fallback PCA outside this function")
    try:
        import scvi
        import torch

        scvi.settings.seed = seed
        scvi.model.SCVI.setup_anndata(hvg, layer="counts", batch_key="source_gse_id")
        model = scvi.model.SCVI(hvg, n_latent=30, n_layers=2, gene_likelihood="nb")
        train_kwargs = dict(max_epochs=epochs, early_stopping=True, early_stopping_patience=10, batch_size=1024)
        if torch.cuda.is_available():
            train_kwargs.update({"accelerator": "gpu", "devices": 1})
        try:
            model.train(**train_kwargs)
        except TypeError:
            train_kwargs.pop("accelerator", None)
            train_kwargs.pop("devices", None)
            model.train(**train_kwargs, use_gpu=torch.cuda.is_available())
        latent = model.get_latent_representation()
        adata.obsm["X_scVI_gdT"] = latent.astype(np.float32, copy=False)
        adata.obsm["X_scVI"] = adata.obsm["X_scVI_gdT"]
        model_path = MODEL_DIR / "scvi_gdt_reintegration"
        model.save(model_path, overwrite=True)
        integration_info["scvi_model_dir"] = str(model_path)
        integration_info["scvi_epochs_requested"] = epochs
        history = getattr(model, "history", None)
        if history is not None:
            try:
                pd.DataFrame(history).to_csv(TABLE_DIR / "scvi_training_history.csv", index=False)
            except Exception:
                pass
        return integration_info
    except Exception as exc:
        logging.exception("scVI integration failed; falling back to PCA: %s", exc)
        integration_info["integration_method"] = "PCA_fallback"
        integration_info["scvi_error"] = repr(exc)
        sc.tl.pca(hvg, n_comps=50, svd_solver="arpack", random_state=seed)
        adata.obsm["X_scVI_gdT"] = hvg.obsm["X_pca"].astype(np.float32, copy=False)
        adata.obsm["X_scVI"] = adata.obsm["X_scVI_gdT"]
        return integration_info


def run_embedding_and_clustering(adata: ad.AnnData, seed: int) -> dict[str, Any]:
    info: dict[str, Any] = {}
    try:
        import rapids_singlecell as rsc

        rsc.pp.neighbors(adata, n_neighbors=30, use_rep="X_scVI_gdT")
        rsc.tl.umap(adata, min_dist=0.35, spread=1.0, random_state=seed)
        info["neighbors_umap_backend"] = "rapids_singlecell"
    except Exception as exc:
        logging.warning("RAPIDS neighbors/UMAP failed; using Scanpy CPU backend: %s", exc)
        sc.pp.neighbors(adata, n_neighbors=30, use_rep="X_scVI_gdT", random_state=seed)
        sc.tl.umap(adata, min_dist=0.35, spread=1.0, random_state=seed)
        info["neighbors_umap_backend"] = "scanpy_cpu"
    resolution_rows = []
    for res in [0.3, 0.6, 0.9, 1.2]:
        key = f"leiden_gdt_{str(res).replace('.', 'p')}"
        sc.tl.leiden(adata, resolution=res, key_added=key, random_state=seed, flavor="igraph", directed=False)
        counts = adata.obs[key].value_counts()
        source_dom = (
            adata.obs.groupby(key, observed=False)["source_gse_id"]
            .agg(lambda s: s.value_counts(normalize=True).iloc[0])
            .astype(float)
        )
        resolution_rows.append(
            {
                "resolution": res,
                "key": key,
                "n_clusters": int(counts.shape[0]),
                "min_cluster_size": int(counts.min()),
                "median_top_source_fraction": float(source_dom.median()),
            }
        )
    res_table = pd.DataFrame(resolution_rows)
    res_table.to_csv(TABLE_DIR / "leiden_resolution_audit.csv", index=False)
    acceptable = res_table[(res_table["n_clusters"].between(8, 40)) & (res_table["min_cluster_size"] >= 200)]
    chosen = acceptable.iloc[0] if not acceptable.empty else res_table.iloc[(res_table["resolution"] - 0.6).abs().argmin()]
    adata.obs["gdt_leiden"] = adata.obs[str(chosen["key"])].astype(str).astype("category")
    info["chosen_leiden_resolution"] = float(chosen["resolution"])
    info["chosen_leiden_key"] = str(chosen["key"])
    info["chosen_n_clusters"] = int(chosen["n_clusters"])
    return info


def summarize_curation(curation: pd.DataFrame) -> pd.DataFrame:
    rows = [
        ("promoted_gdTAI_predictions", curation["gdtai_promoted_prediction"]),
        ("removed_likely_FP", curation["removed_fp"]),
        ("added_gold_FN", curation["gold_fn_added"]),
        ("excluded_silver_FN", curation["silver_fn_excluded"]),
        ("curated_primary_gdT", curation["curated_primary_gdT"]),
    ]
    return pd.DataFrame(
        [{"category": name, "n_cells": int(mask.sum()), "fraction_of_total": float(mask.mean())} for name, mask in rows]
    )


def group_curation(df: pd.DataFrame, group_cols: list[str], path: Path) -> pd.DataFrame:
    table = (
        df.groupby(group_cols, dropna=False)
        .agg(
            total_TNK_cells=("curated_primary_gdT", "size"),
            promoted_gdTAI_predictions=("gdtai_promoted_prediction", "sum"),
            removed_likely_FP=("removed_fp", "sum"),
            added_gold_FN=("gold_fn_added", "sum"),
            excluded_silver_FN=("silver_fn_excluded", "sum"),
            curated_primary_gdT=("curated_primary_gdT", "sum"),
        )
        .reset_index()
    )
    table["curated_gdT_fraction"] = table["curated_primary_gdT"] / table["total_TNK_cells"].replace(0, np.nan)
    table = table.sort_values(["curated_primary_gdT", "total_TNK_cells"], ascending=False)
    table.to_csv(path, index=False)
    return table


def fraction_tests(
    sample_table: pd.DataFrame,
    group_col: str,
    success_col: str,
    total_col: str,
    min_samples: int = 3,
    min_success: int = 20,
) -> pd.DataFrame:
    rows = []
    data = sample_table.copy()
    data["fraction"] = data[success_col] / data[total_col].replace(0, np.nan)
    levels = [x for x in data[group_col].dropna().astype(str).unique() if x and x != "unknown"]
    for level in levels:
        in_group = data[group_col].astype(str).eq(level)
        a = data.loc[in_group, "fraction"].dropna().to_numpy(dtype=float)
        b = data.loc[~in_group, "fraction"].dropna().to_numpy(dtype=float)
        if a.size < min_samples or b.size < min_samples or data.loc[in_group, success_col].sum() < min_success:
            continue
        stat, pvalue = stats.mannwhitneyu(a, b, alternative="two-sided")
        rows.append(
            {
                "test": f"{level}_vs_rest",
                "group_col": group_col,
                "n_samples_group": int(a.size),
                "n_samples_rest": int(b.size),
                "success_group": int(data.loc[in_group, success_col].sum()),
                "total_group": int(data.loc[in_group, total_col].sum()),
                "median_fraction_group": float(np.median(a)),
                "median_fraction_rest": float(np.median(b)),
                "median_fraction_delta": float(np.median(a) - np.median(b)),
                "p_value": float(pvalue),
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr_bh"] = multipletests(out["p_value"].to_numpy(dtype=float), method="fdr_bh")[1]
        out = out.sort_values(["fdr_bh", "p_value"])
    return out


def build_sample_abundance(df: pd.DataFrame) -> pd.DataFrame:
    sample_cols = [
        "sample_unit_gdt_atlas",
        "source_gse_id",
        "tissue_site_gdt_atlas",
        "disease_status_gdt_atlas",
        "age_group_gdt_atlas",
        "sex_gdt_atlas",
    ]
    sample = (
        df.groupby(sample_cols, dropna=False)
        .agg(total_TNK_cells=("curated_primary_gdT", "size"), curated_gdT=("curated_primary_gdT", "sum"))
        .reset_index()
    )
    sample["curated_gdT_fraction"] = sample["curated_gdT"] / sample["total_TNK_cells"].replace(0, np.nan)
    sample.to_csv(TABLE_DIR / "sample_level_gdt_abundance.csv", index=False)
    stats_tables = []
    for col in ["tissue_site_gdt_atlas", "disease_status_gdt_atlas", "age_group_gdt_atlas", "sex_gdt_atlas"]:
        tests = fraction_tests(sample, col, "curated_gdT", "total_TNK_cells")
        if not tests.empty:
            stats_tables.append(tests)
    all_tests = pd.concat(stats_tables, ignore_index=True) if stats_tables else pd.DataFrame()
    all_tests.to_csv(TABLE_DIR / "sample_level_gdt_abundance_statistics.csv", index=False)
    return sample


def phenotype_composition_stats(adata: ad.AnnData) -> tuple[pd.DataFrame, pd.DataFrame]:
    group_cols = [
        "sample_unit_gdt_atlas",
        "source_gse_id",
        "tissue_site_gdt_atlas",
        "disease_status_gdt_atlas",
        "age_group_gdt_atlas",
        "gdt_phenotype_state",
    ]
    obs = adata.obs[group_cols].copy()
    for column in group_cols:
        obs[column] = obs[column].astype(str)
    comp = obs.groupby(group_cols, observed=True, dropna=False).size().reset_index(name="phenotype_cells")
    total = obs.groupby("sample_unit_gdt_atlas", observed=True).size().rename("sample_gdT_cells").reset_index()
    comp = comp.merge(total, on="sample_unit_gdt_atlas", how="left")
    comp["phenotype_fraction"] = comp["phenotype_cells"] / comp["sample_gdT_cells"].replace(0, np.nan)
    comp.to_csv(TABLE_DIR / "sample_level_phenotype_composition.csv", index=False)

    stats_rows = []
    for phenotype in sorted(obs["gdt_phenotype_state"].astype(str).unique()):
        sub = comp.loc[comp["gdt_phenotype_state"].astype(str).eq(phenotype)].copy()
        if sub["phenotype_cells"].sum() < 50:
            continue
        for col in ["tissue_site_gdt_atlas", "disease_status_gdt_atlas", "age_group_gdt_atlas"]:
            tests = fraction_tests(sub, col, "phenotype_cells", "sample_gdT_cells", min_success=10)
            if not tests.empty:
                tests.insert(0, "gdt_phenotype_state", phenotype)
                stats_rows.append(tests)
    stats_table = pd.concat(stats_rows, ignore_index=True) if stats_rows else pd.DataFrame()
    if not stats_table.empty:
        stats_table["global_fdr_bh"] = multipletests(stats_table["p_value"], method="fdr_bh")[1]
        stats_table = stats_table.sort_values(["global_fdr_bh", "fdr_bh", "p_value"])
    stats_table.to_csv(TABLE_DIR / "sample_level_phenotype_composition_statistics.csv", index=False)
    return comp, stats_table


def summarize_clusters(adata: ad.AnnData, score_cols: list[str], v_cols: list[str]) -> pd.DataFrame:
    obs = adata.obs.copy()
    rows = []
    for cluster, sub in obs.groupby("gdt_leiden", observed=False):
        source_counts = sub["source_gse_id"].astype(str).value_counts()
        tissue_counts = sub["tissue_site_gdt_atlas"].astype(str).value_counts()
        phenotype_counts = sub["gdt_phenotype_state"].astype(str).value_counts()
        row = {
            "gdt_leiden": str(cluster),
            "n_cells": int(sub.shape[0]),
            "n_sources": int(source_counts.shape[0]),
            "n_samples": int(sub["sample_unit_gdt_atlas"].nunique()),
            "top_source": source_counts.index[0] if not source_counts.empty else "",
            "top_source_fraction": float(source_counts.iloc[0] / sub.shape[0]) if not source_counts.empty else np.nan,
            "top_tissue": tissue_counts.index[0] if not tissue_counts.empty else "",
            "top_tissue_fraction": float(tissue_counts.iloc[0] / sub.shape[0]) if not tissue_counts.empty else np.nan,
            "dominant_phenotype": phenotype_counts.index[0] if not phenotype_counts.empty else "",
            "dominant_phenotype_fraction": float(phenotype_counts.iloc[0] / sub.shape[0]) if not phenotype_counts.empty else np.nan,
        }
        for col in score_cols + v_cols:
            if col in sub:
                row[f"mean_{col}"] = float(pd.to_numeric(sub[col], errors="coerce").mean())
        rows.append(row)
    table = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
    table.to_csv(TABLE_DIR / "cluster_reproducibility_and_marker_summary.csv", index=False)
    return table


def pseudobulk_de(adata: ad.AnnData, marker_config: dict[str, Any], max_genes: int = 2500) -> pd.DataFrame:
    gene_names = pd.Index(adata.var_names.astype(str))
    tcr_prefixes = tuple(marker_config.get("tcr_gene_prefixes_excluded_from_hvg_and_de", ["TRA", "TRB", "TRG", "TRD"]))
    allowed = ~gene_names.str.startswith(tcr_prefixes)
    hvg = adata.var["highly_variable_no_tcr"].to_numpy(dtype=bool) if "highly_variable_no_tcr" in adata.var else allowed
    selected = gene_names[allowed & hvg].tolist()[:max_genes]
    if len(selected) < 200:
        selected = gene_names[allowed].tolist()[:max_genes]
    gene_idx = np.asarray([adata.var_names.get_loc(g) for g in selected], dtype=np.int64)
    x = adata.layers["counts"][:, gene_idx].tocsr()
    total_counts = np.asarray(adata.layers["counts"].sum(axis=1)).ravel()
    sample_codes, sample_names = pd.factorize(adata.obs["sample_unit_gdt_atlas"].astype(str), sort=True)
    n_samples = len(sample_names)
    sample_mat = sparse.csr_matrix(
        (np.ones(adata.n_obs, dtype=np.float32), (sample_codes, np.arange(adata.n_obs))),
        shape=(n_samples, adata.n_obs),
    )
    sample_counts = sample_mat @ x
    sample_lib = np.asarray(sample_mat @ total_counts.reshape(-1, 1)).ravel()
    sample_ncells = np.bincount(sample_codes, minlength=n_samples)

    rows = []
    phenotypes = adata.obs["gdt_phenotype_state"].astype(str).value_counts()
    for phenotype, n_cells in phenotypes.items():
        if n_cells < 100:
            continue
        mask = adata.obs["gdt_phenotype_state"].astype(str).eq(phenotype).to_numpy(dtype=bool)
        ph_codes = sample_codes[mask]
        if np.unique(ph_codes).size < 3:
            continue
        ph_mat = sparse.csr_matrix(
            (np.ones(mask.sum(), dtype=np.float32), (ph_codes, np.arange(mask.sum()))),
            shape=(n_samples, int(mask.sum())),
        )
        ph_counts = ph_mat @ x[mask]
        ph_lib = np.asarray(ph_mat @ total_counts[mask].reshape(-1, 1)).ravel()
        ph_ncells = np.bincount(ph_codes, minlength=n_samples)
        rest_counts = sample_counts - ph_counts
        rest_lib = sample_lib - ph_lib
        rest_ncells = sample_ncells - ph_ncells
        valid = (ph_ncells >= 20) & (rest_ncells >= 20) & (ph_lib > 0) & (rest_lib > 0)
        if int(valid.sum()) < 3:
            continue
        ph_log = np.log2((ph_counts[valid].toarray() + 0.5) / (ph_lib[valid, None] + 1.0) * 1e6)
        rest_log = np.log2((rest_counts[valid].toarray() + 0.5) / (rest_lib[valid, None] + 1.0) * 1e6)
        diff = ph_log - rest_log
        t_stat, pvals = stats.ttest_1samp(diff, 0.0, axis=0, nan_policy="omit")
        pvals = np.nan_to_num(pvals, nan=1.0, posinf=1.0, neginf=1.0)
        fdr = multipletests(pvals, method="fdr_bh")[1]
        mean_diff = np.nanmean(diff, axis=0)
        for gene, effect, pval, qval, tval in zip(selected, mean_diff, pvals, fdr, np.nan_to_num(t_stat, nan=0.0)):
            rows.append(
                {
                    "gdt_phenotype_state": phenotype,
                    "gene": gene,
                    "mean_log2cpm_delta_vs_other_gdT": float(effect),
                    "t_statistic": float(tval),
                    "p_value": float(pval),
                    "fdr_bh": float(qval),
                    "n_paired_samples": int(valid.sum()),
                    "phenotype_cells": int(n_cells),
                }
            )
    de = pd.DataFrame(rows)
    if de.empty:
        de.to_csv(TABLE_DIR / "phenotype_pseudobulk_de_all.csv", index=False)
        return de
    de = de.sort_values(["gdt_phenotype_state", "fdr_bh", "p_value"])
    de.to_csv(TABLE_DIR / "phenotype_pseudobulk_de_all.csv.gz", index=False)
    top = (
        de.assign(abs_effect=de["mean_log2cpm_delta_vs_other_gdT"].abs())
        .sort_values(["gdt_phenotype_state", "fdr_bh", "abs_effect"], ascending=[True, True, False])
        .groupby("gdt_phenotype_state", as_index=False)
        .head(50)
        .drop(columns="abs_effect")
    )
    top.to_csv(TABLE_DIR / "phenotype_pseudobulk_de_top50_per_phenotype.csv", index=False)
    return de


def plot_curation_flow(summary: pd.DataFrame, path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(8.8, 4.6), constrained_layout=True)
    colors = ["#2563eb", "#dc2626", "#16a34a", "#6b7280", "#7c3aed"]
    ax.bar(summary["category"], summary["n_cells"], color=colors[: summary.shape[0]])
    ax.set_ylabel("Cells")
    ax.set_title("gdT atlas curation flow")
    ax.tick_params(axis="x", rotation=25)
    for i, value in enumerate(summary["n_cells"]):
        ax.text(i, value, f"{int(value):,}", ha="center", va="bottom", fontsize=9)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_umap_category(adata: ad.AnnData, color: str, path: Path, seed: int, max_cells: int = 160000) -> Path:
    rows = sample_indices(np.arange(adata.n_obs, dtype=np.int64), max_cells, seed)
    umap = adata.obsm["X_umap"][rows]
    values = pd.Series(adata.obs[color].astype(str).to_numpy()[rows])
    top = values.value_counts().head(18).index.tolist()
    values = values.where(values.isin(top), "other")
    palette = sns.color_palette("tab20", n_colors=values.nunique())
    color_map = dict(zip(values.value_counts().index, palette))
    fig, ax = plt.subplots(figsize=(7.0, 6.2), constrained_layout=True)
    ax.scatter(umap[:, 0], umap[:, 1], s=2, c=values.map(color_map).tolist(), linewidths=0, alpha=0.65, rasterized=True)
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color_map[k], label=k, markersize=5)
        for k in color_map
    ]
    ax.legend(handles=handles, frameon=False, fontsize=7, markerscale=1, loc="best", ncol=1)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title(color.replace("_", " "))
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_score_heatmap(table: pd.DataFrame, index_col: str, score_cols: list[str], path: Path, title: str) -> Path:
    if table.empty or not score_cols:
        return path
    mat = table.set_index(index_col)[score_cols]
    mat.columns = [c.replace("score_", "").replace("expr_", "") for c in mat.columns]
    fig, ax = plt.subplots(figsize=(max(8, 0.45 * mat.shape[1] + 3), max(4, 0.35 * mat.shape[0] + 1.8)), constrained_layout=True)
    sns.heatmap(mat, cmap="mako", ax=ax, cbar_kws={"label": "mean log-normalized expression"})
    ax.set_title(title)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_abundance(sample: pd.DataFrame, group_col: str, path: Path, title: str, top_n: int = 12) -> Path:
    plot = sample.loc[sample[group_col].astype(str).ne("unknown")].copy()
    top = plot.groupby(group_col)["curated_gdT"].sum().sort_values(ascending=False).head(top_n).index
    plot = plot.loc[plot[group_col].isin(top)].copy()
    fig, ax = plt.subplots(figsize=(max(8, 0.55 * len(top) + 3), 4.8), constrained_layout=True)
    sns.boxplot(data=plot, x=group_col, y="curated_gdT_fraction", ax=ax, color="#dbeafe", fliersize=0)
    sns.stripplot(data=plot, x=group_col, y="curated_gdT_fraction", ax=ax, color="#1f2937", size=2, alpha=0.55)
    ax.set_title(title)
    ax.set_xlabel(group_col.replace("_", " "))
    ax.set_ylabel("gdT fraction among T/NK cells")
    ax.tick_params(axis="x", rotation=35)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_composition_heatmap(comp: pd.DataFrame, group_col: str, path: Path) -> Path:
    plot = comp.loc[comp[group_col].astype(str).ne("unknown")].copy()
    if plot.empty:
        return path
    top_groups = plot.groupby(group_col)["phenotype_cells"].sum().sort_values(ascending=False).head(14).index
    plot = plot.loc[plot[group_col].isin(top_groups)]
    mat = plot.pivot_table(
        index="gdt_phenotype_state",
        columns=group_col,
        values="phenotype_fraction",
        aggfunc="mean",
        fill_value=0,
    )
    fig, ax = plt.subplots(figsize=(max(8, 0.5 * mat.shape[1] + 3), max(4, 0.35 * mat.shape[0] + 1.8)), constrained_layout=True)
    sns.heatmap(mat, cmap="viridis", ax=ax, cbar_kws={"label": "mean sample fraction"})
    ax.set_title(f"Phenotype composition by {group_col.replace('_', ' ')}")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_de_volcano(de: pd.DataFrame, path: Path) -> Path:
    if de.empty:
        return path
    phenotypes = de.groupby("gdt_phenotype_state")["phenotype_cells"].max().sort_values(ascending=False).head(4).index
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), constrained_layout=True)
    for ax, phenotype in zip(axes.reshape(-1), phenotypes):
        sub = de.loc[de["gdt_phenotype_state"].eq(phenotype)].copy()
        sub["neglog10_fdr"] = -np.log10(sub["fdr_bh"].clip(lower=1e-300))
        ax.scatter(
            sub["mean_log2cpm_delta_vs_other_gdT"],
            sub["neglog10_fdr"],
            s=5,
            c=np.where(sub["fdr_bh"] < 0.05, "#dc2626", "#94a3b8"),
            alpha=0.7,
            linewidths=0,
        )
        top = sub.sort_values(["fdr_bh", "mean_log2cpm_delta_vs_other_gdT"], ascending=[True, False]).head(8)
        for _, row in top.iterrows():
            ax.text(row["mean_log2cpm_delta_vs_other_gdT"], row["neglog10_fdr"], row["gene"], fontsize=6)
        ax.axvline(0, color="#111827", linewidth=0.7)
        ax.set_title(str(phenotype))
        ax.set_xlabel("sample-paired log2CPM delta")
        ax.set_ylabel("-log10 FDR")
    for ax in axes.reshape(-1)[len(phenotypes) :]:
        ax.axis("off")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_v_gene_umaps(adata: ad.AnnData, genes: list[str], path: Path, seed: int) -> Path:
    present = [f"expr_{g}" for g in genes if f"expr_{g}" in adata.obs]
    if not present:
        return path
    rows = sample_indices(np.arange(adata.n_obs, dtype=np.int64), 120000, seed)
    umap = adata.obsm["X_umap"][rows]
    ncols = min(3, len(present))
    nrows = int(math.ceil(len(present) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.7 * nrows), constrained_layout=True)
    axes_arr = np.asarray(axes).reshape(-1)
    for ax, col in zip(axes_arr, present):
        values = pd.to_numeric(adata.obs[col], errors="coerce").fillna(0).to_numpy(dtype=float)[rows]
        scatter = ax.scatter(umap[:, 0], umap[:, 1], c=values, s=2, cmap="magma", linewidths=0, rasterized=True)
        ax.set_title(col.replace("expr_", ""))
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.02)
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
    manifest: dict[str, Any],
    curation_summary: pd.DataFrame,
    curation_by_source: pd.DataFrame,
    curation_by_tissue: pd.DataFrame,
    abundance_stats: pd.DataFrame,
    phenotype_counts: pd.DataFrame,
    cluster_summary: pd.DataFrame,
    de_top: pd.DataFrame,
    figures: list[Path],
) -> None:
    summary_lines = [
        "# Rigorous gdT Atlas Reanalysis",
        "",
        f"- Input H5AD: `{manifest['input_h5ad']}`",
        f"- Total T/NK atlas cells: `{manifest['n_obs']:,}`",
        f"- Promoted gdTAI high-F1 predictions: `{manifest['promoted_predictions']:,}`",
        f"- Removed likely FP cells: `{manifest['removed_likely_fp']:,}`",
        f"- Added gold-level FN cells: `{manifest['added_gold_fn']:,}`",
        f"- Excluded silver-only FN cells: `{manifest['excluded_silver_fn']:,}`",
        f"- Final gdT atlas cells: `{manifest['curated_primary_gdT']:,}`",
        f"- gdT-only integration method: `{manifest['integration_method']}`",
        f"- Source H5AD unchanged: `{manifest['source_h5ad_unchanged']}`",
        "",
        "## Methods",
        "",
        "- Silver-only FN cells were not included in the gdT biology atlas.",
        "- TCR genes were excluded from HVG selection and primary DE.",
        "- A new gdT-only embedding, Leiden clustering, and UMAP were computed.",
        "- Abundance and composition comparisons use sample-level fractions with FDR correction.",
        "",
        "## Outputs",
        f"- HTML: `{REPORT_HTML}`",
        f"- PDF: `{REPORT_PDF}`",
        f"- Derived gdT-only H5AD: `{manifest['derived_h5ad']}`",
        f"- Tables: `{TABLE_DIR}`",
        f"- Figures: `{FIGURE_DIR}`",
    ]
    SUMMARY_MD.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    figure_html = "\n".join(
        f"<figure><img src='../../{rel(fig)}'><figcaption>{html.escape(fig.stem.replace('_', ' '))}</figcaption></figure>"
        for fig in figures
        if fig.exists()
    )
    css = """
    @page{size:A4 landscape;margin:8mm}
    body{font-family:Arial,Helvetica,sans-serif;margin:24px auto;max-width:1380px;color:#17202a;line-height:1.45}
    h1{font-size:30px;margin:0 0 8px} h2{font-size:22px;margin-top:26px} h3{font-size:16px;margin-top:18px}
    .grid{display:grid;grid-template-columns:repeat(5,1fr);gap:12px;margin:16px 0}
    .metric{border:1px solid #d9e2ec;background:#f8fafc;padding:12px}.value{font-size:23px;font-weight:bold}
    .note{background:#fff7ed;border-left:4px solid #f97316;padding:10px;margin:12px 0}
    table{border-collapse:collapse;width:100%;font-size:9px;table-layout:fixed}
    th,td{border:1px solid #d8dee5;padding:4px 5px;text-align:left;vertical-align:top;overflow-wrap:anywhere}
    th{background:#eef2f6} img{max-width:100%;height:auto;border:1px solid #d8dee5}
    figure{break-inside:avoid;margin:18px 0} code{background:#f1f5f9;padding:1px 3px}
    """
    html_text = f"""<!doctype html><html><head><meta charset='utf-8'>
<title>Rigorous gdT Atlas Reanalysis</title><style>{css}</style></head><body>
<h1>Rigorous gdT Atlas Reanalysis</h1>
<p>This report rebuilds the gdT atlas from the 5 million cell T/NK atlas using promoted gdTAI predictions, explicit FP removal, gold-level FN add-back, new gdT-only integration, Leiden clustering, UMAP, and sample-level statistics.</p>
<div class='grid'>
  <div class='metric'><div>Promoted predictions</div><div class='value'>{manifest['promoted_predictions']:,}</div></div>
  <div class='metric'><div>Removed FP</div><div class='value'>{manifest['removed_likely_fp']:,}</div></div>
  <div class='metric'><div>Gold FN added</div><div class='value'>{manifest['added_gold_fn']:,}</div></div>
  <div class='metric'><div>Silver FN excluded</div><div class='value'>{manifest['excluded_silver_fn']:,}</div></div>
  <div class='metric'><div>Final gdT atlas</div><div class='value'>{manifest['curated_primary_gdT']:,}</div></div>
</div>
<section class='note'>Silver-only FN cells are excluded from all biology analyses in this report. TCR genes are excluded from HVG selection and primary DE; TRD/TRG V/J/C RNA expression is summarized separately.</section>
<h2>Methodology</h2>
<p>Ground truth add-back uses only gold-level missed cells from the validation-positive audit table. The biology atlas starts from the promoted annotation-specific gdTAI high-F1 prediction mask, removes likely false positives supported by paired TCRAB, NK/CD4/Treg annotation, low CD3/gdT marker signal, NK-like expression, CD4/FOXP3-like expression, or myeloid/B contamination when no gd evidence is present, and adds back only gold FN cells. The source H5AD is never mutated.</p>
<p>The gdT-only matrix was extracted from count-space <code>X</code>, normalized/log-transformed for marker visualization, and integrated de novo with {html.escape(manifest['integration_method'])}. HVGs used for integration excluded TRA/TRB/TRG/TRD genes. Leiden clustering was tested at 0.3, 0.6, 0.9, and 1.2; the selected resolution was {manifest.get('chosen_leiden_resolution', 'NA')}.</p>
<p>Abundance and phenotype-composition claims are based on sample/library-level fractions with Mann-Whitney tests and Benjamini-Hochberg FDR correction. Phenotype DE uses sample-paired pseudobulk logCPM deltas between each phenotype and other gdT cells, with TCR genes excluded from the DE list.</p>
<h2>Figures</h2>{figure_html}
<h2>Curation Summary</h2>{dataframe_to_html(curation_summary, 20)}
<h2>Curated gdT By Source</h2>{dataframe_to_html(curation_by_source, 35)}
<h2>Curated gdT By Tissue</h2>{dataframe_to_html(curation_by_tissue, 35)}
<h2>Sample-Level Abundance Statistics</h2>{dataframe_to_html(abundance_stats, 60)}
<h2>Phenotype Counts</h2>{dataframe_to_html(phenotype_counts, 40)}
<h2>Cluster Reproducibility</h2>{dataframe_to_html(cluster_summary, 60)}
<h2>Top Pseudobulk DE</h2>{dataframe_to_html(de_top, 80)}
</body></html>"""
    REPORT_HTML.write_text(html_text, encoding="utf-8")



def _summary_count(summary: pd.DataFrame, category: str) -> int:
    hit = summary.loc[summary["category"].astype(str).eq(category), "n_cells"]
    return int(hit.iloc[0]) if not hit.empty else 0


def finalize_from_derived(args: argparse.Namespace, marker_config: dict[str, Any]) -> None:
    logging.info("Resuming from derived gdT-only H5AD: %s", args.output_h5ad)
    if not args.output_h5ad.exists():
        raise FileNotFoundError(args.output_h5ad)
    source_before = args.input_h5ad.stat()
    adata = sc.read_h5ad(args.output_h5ad)

    curation_summary = pd.read_csv(TABLE_DIR / "curation_summary.csv")
    curation_by_source = pd.read_csv(TABLE_DIR / "curation_by_source.csv")
    curation_by_tissue = pd.read_csv(TABLE_DIR / "curation_by_tissue.csv")
    sample_abundance = pd.read_csv(TABLE_DIR / "sample_level_gdt_abundance.csv")
    abundance_stats_path = TABLE_DIR / "sample_level_gdt_abundance_statistics.csv"
    abundance_stats = pd.read_csv(abundance_stats_path) if abundance_stats_path.exists() else pd.DataFrame()

    logging.info("Regenerating corrected phenotype composition tables")
    phenotype_counts = (
        adata.obs.groupby(["gdt_phenotype_state"], observed=True)
        .agg(
            n_cells=("gdt_phenotype_state", "size"),
            n_sources=("source_gse_id", "nunique"),
            n_samples=("sample_unit_gdt_atlas", "nunique"),
            mean_confidence=("gdt_phenotype_confidence", "mean"),
        )
        .reset_index()
        .sort_values("n_cells", ascending=False)
    )
    phenotype_counts.to_csv(TABLE_DIR / "phenotype_counts.csv", index=False)
    phenotype_composition, phenotype_stats = phenotype_composition_stats(adata)

    logging.info("Regenerating cluster, score, V-gene, and metadata tables")
    score_cols = [c for c in adata.obs.columns if str(c).startswith("score_")]
    v_cols = [c for c in adata.obs.columns if str(c).startswith("expr_TR")]
    cluster_summary = summarize_clusters(adata, score_cols, v_cols)
    score_by_pheno = adata.obs.groupby("gdt_phenotype_state", observed=True)[score_cols].mean().reset_index()
    score_by_pheno.to_csv(TABLE_DIR / "marker_program_scores_by_phenotype.csv", index=False)
    v_by_pheno = adata.obs.groupby("gdt_phenotype_state", observed=True)[v_cols].mean().reset_index()
    v_by_pheno.to_csv(TABLE_DIR / "v_gene_expression_by_phenotype.csv", index=False)
    v_by_cluster = adata.obs.groupby("gdt_leiden", observed=True)[v_cols].mean().reset_index()
    v_by_cluster.to_csv(TABLE_DIR / "v_gene_expression_by_cluster.csv", index=False)
    adata.obs.reset_index(names="obs_name").to_csv(TABLE_DIR / "curated_gdt_cell_metadata_with_clusters.csv.gz", index=False)

    logging.info("Regenerating sample-paired pseudobulk DE")
    de = pseudobulk_de(adata, marker_config)
    de_top_path = TABLE_DIR / "phenotype_pseudobulk_de_top50_per_phenotype.csv"
    de_top = pd.read_csv(de_top_path) if de_top_path.exists() else pd.DataFrame()

    logging.info("Regenerating figures")
    figures: list[Path] = []
    figures.append(plot_curation_flow(curation_summary, FIGURE_DIR / "curation_flow_no_silver_addback.png"))
    for color, filename in [
        ("gdt_leiden", "gdt_only_umap_by_leiden.png"),
        ("gdt_phenotype_state", "gdt_only_umap_by_phenotype.png"),
        ("gdt_v_gene_state", "gdt_only_umap_by_trdv_state.png"),
        ("source_gse_id", "gdt_only_umap_by_source.png"),
        ("tissue_site_gdt_atlas", "gdt_only_umap_by_tissue.png"),
        ("disease_status_gdt_atlas", "gdt_only_umap_by_disease.png"),
        ("age_group_gdt_atlas", "gdt_only_umap_by_age.png"),
        ("gdt_atlas_evidence_tier", "gdt_only_umap_by_evidence_tier.png"),
    ]:
        figures.append(plot_umap_category(adata, color, FIGURE_DIR / filename, args.seed + len(figures)))
    figures.append(plot_score_heatmap(score_by_pheno, "gdt_phenotype_state", score_cols, FIGURE_DIR / "marker_program_scores_by_phenotype.png", "Marker programs by phenotype"))
    figures.append(plot_score_heatmap(v_by_pheno, "gdt_phenotype_state", v_cols, FIGURE_DIR / "v_gene_expression_by_phenotype.png", "TRD/TRG V/J/C RNA expression by phenotype"))
    figures.append(plot_score_heatmap(v_by_cluster, "gdt_leiden", v_cols, FIGURE_DIR / "v_gene_expression_by_cluster.png", "TRD/TRG V/J/C RNA expression by Leiden cluster"))
    figures.append(plot_abundance(sample_abundance, "tissue_site_gdt_atlas", FIGURE_DIR / "sample_level_gdt_abundance_by_tissue.png", "Sample-level gdT abundance by tissue"))
    figures.append(plot_abundance(sample_abundance, "disease_status_gdt_atlas", FIGURE_DIR / "sample_level_gdt_abundance_by_disease.png", "Sample-level gdT abundance by disease"))
    figures.append(plot_composition_heatmap(phenotype_composition, "tissue_site_gdt_atlas", FIGURE_DIR / "phenotype_composition_by_tissue.png"))
    figures.append(plot_composition_heatmap(phenotype_composition, "disease_status_gdt_atlas", FIGURE_DIR / "phenotype_composition_by_disease.png"))
    figures.append(plot_de_volcano(de, FIGURE_DIR / "phenotype_pseudobulk_de_volcano.png"))
    figures.append(plot_v_gene_umaps(adata, marker_config.get("primary_v_gene_panels", []), FIGURE_DIR / "v_gene_expression_umap_panels.png", args.seed))

    resolution = "NA"
    audit_path = TABLE_DIR / "leiden_resolution_audit.csv"
    if audit_path.exists():
        audit = pd.read_csv(audit_path)
        n_clusters = int(adata.obs["gdt_leiden"].nunique())
        hit = audit.loc[audit["n_clusters"].astype(int).eq(n_clusters)]
        if not hit.empty:
            resolution = float(hit.iloc[0]["resolution"])

    source_after = args.input_h5ad.stat()
    manifest = {
        "input_h5ad": str(args.input_h5ad),
        "n_obs": EXPECTED_N_OBS,
        "integration_method": "scVI",
        "chosen_leiden_resolution": resolution,
        "prediction_npz": str(args.prediction_npz),
        "prediction_key": PREDICTION_KEY,
        "promoted_predictions": _summary_count(curation_summary, "promoted_gdTAI_predictions"),
        "removed_likely_fp": _summary_count(curation_summary, "removed_likely_FP"),
        "added_gold_fn": _summary_count(curation_summary, "added_gold_FN"),
        "excluded_silver_fn": _summary_count(curation_summary, "excluded_silver_FN"),
        "curated_primary_gdT": int(adata.n_obs),
        "derived_h5ad": str(args.output_h5ad),
        "source_size_before": source_before.st_size,
        "source_mtime_ns_before": source_before.st_mtime_ns,
        "source_size_after": source_after.st_size,
        "source_mtime_ns_after": source_after.st_mtime_ns,
        "source_h5ad_unchanged": bool(source_before.st_size == source_after.st_size and source_before.st_mtime_ns == source_after.st_mtime_ns),
        "html": str(REPORT_HTML),
        "pdf": str(REPORT_PDF),
        "tables": str(TABLE_DIR),
        "figures": str(FIGURE_DIR),
        "resume_from_derived": True,
    }
    write_report(
        manifest,
        curation_summary,
        curation_by_source,
        curation_by_tissue,
        abundance_stats,
        phenotype_counts,
        cluster_summary,
        de_top,
        figures,
    )
    manifest["pdf_rendered"] = bool(render_pdf(args.no_pdf))
    MANIFEST_JSON.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    logging.info("Resume/finalize complete: %s", REPORT_HTML)

def main() -> None:
    args = parse_args()
    setup_logging()
    np.random.seed(args.seed)
    marker_config = json.loads(args.marker_config.read_text(encoding="utf-8"))
    if args.resume_from_derived:
        finalize_from_derived(args, marker_config)
        return
    metadata_rules = json.loads(args.metadata_rules.read_text(encoding="utf-8"))

    logging.info("Loading source obs")
    df, info = load_obs(args.input_h5ad)
    pred_arrays = load_prediction_arrays(args.prediction_npz, int(info["n_obs"]))
    pred = pred_arrays[PREDICTION_KEY]
    gold_fn, gold_fn_summary = load_gold_fn_mask(args.fn_table, int(info["n_obs"]))
    gold_fn_summary.to_csv(TABLE_DIR / "gold_validation_tp_fn_summary.csv", index=False)

    logging.info("Applying deterministic metadata harmonization")
    df = apply_metadata_lookup(df, load_metadata_lookup(args.metadata))
    df = add_harmonized_metadata(df, metadata_rules)
    df["sample_unit_gdt_atlas"] = make_sample_unit(df)

    logging.info("Building curation masks without silver add-back")
    curation = build_initial_curation(df, pred, gold_fn)
    marker_genes = ordered_unique(
        [gene for genes in marker_config["programs"].values() for gene in genes]
        + marker_config.get("v_gene_expression", [])
    )
    candidate_rows = np.flatnonzero(pred | gold_fn)
    expr, total_counts, present_genes, missing_genes = load_selected_expression(
        args.input_h5ad,
        candidate_rows,
        marker_genes,
        args.expression_row_block,
    )
    pd.DataFrame({"present_gene": present_genes}).to_csv(TABLE_DIR / "marker_present_genes.csv", index=False)
    pd.DataFrame({"missing_gene": missing_genes}).to_csv(TABLE_DIR / "marker_missing_genes.csv", index=False)
    scores, score_manifest = score_programs(expr, marker_config)
    scores.index = candidate_rows
    expr.index = candidate_rows
    score_manifest.to_csv(TABLE_DIR / "marker_program_manifest.csv", index=False)
    curation = add_expression_fp_flags(curation, candidate_rows, scores)

    for column in curation.columns:
        df[column] = curation[column].to_numpy()
    if int((df["silver_fn_excluded"] & df["curated_primary_gdT"]).sum()) != 0:
        raise RuntimeError("Silver-only FN cells leaked into curated primary gdT set")

    curation_summary = summarize_curation(curation)
    curation_summary.to_csv(TABLE_DIR / "curation_summary.csv", index=False)
    curation_by_source = group_curation(df, ["source_gse_id"], TABLE_DIR / "curation_by_source.csv")
    curation_by_tissue = group_curation(df, ["tissue_site_gdt_atlas"], TABLE_DIR / "curation_by_tissue.csv")
    group_curation(df, ["disease_status_gdt_atlas"], TABLE_DIR / "curation_by_disease.csv")
    group_curation(df, ["source_gse_id", "tissue_site_gdt_atlas"], TABLE_DIR / "curation_by_source_tissue.csv")

    audit_cols = [
        "obs_row",
        "obs_name",
        "source_gse_id",
        "sample_id",
        "library_id",
        "simple_annotation_plus6",
        "tissue_site_gdt_atlas",
        "disease_status_gdt_atlas",
        "fp_reason",
        "strong_gd_evidence",
        "has_TRA_TRB_paired",
        "has_TRG_TRD_paired",
        "has_any_ab_tcr",
        "has_any_gd_tcr",
        "phase4_trd_score",
        "phase4_trab_score",
        "phase4_trd_minus_trab",
    ]
    df.loc[df["removed_fp"], audit_cols].to_csv(TABLE_DIR / "removed_likely_fp_cells.csv.gz", index=False)
    df.loc[df["gold_fn_added"], audit_cols].to_csv(TABLE_DIR / "gold_fn_added_cells.csv.gz", index=False)
    df.loc[df["silver_fn_excluded"], audit_cols].to_csv(TABLE_DIR / "silver_fn_excluded_not_used_cells.csv.gz", index=False)

    logging.info("Building sample-level abundance statistics")
    sample_abundance = build_sample_abundance(df)
    abundance_stats = pd.read_csv(TABLE_DIR / "sample_level_gdt_abundance_statistics.csv")

    primary_rows = np.flatnonzero(df["curated_primary_gdT"].to_numpy(dtype=bool))
    expected_primary = int(pred.sum()) - int(df["removed_fp"].sum()) + int(gold_fn.sum())
    if primary_rows.size != expected_primary:
        raise RuntimeError(
            f"Curated primary count {primary_rows.size:,} does not match formula {expected_primary:,}"
        )

    logging.info("Extracting full gdT-only count matrix")
    x_primary = extract_csr_rows(args.input_h5ad, primary_rows, args.extract_row_block)
    var = read_var_dataframe(args.input_h5ad)
    primary_obs_cols = [
        "obs_row",
        "source_gse_id",
        "sample_unit_gdt_atlas",
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
        "gold_fn_added",
        "gdtai_promoted_prediction",
        "phase4_trd_score",
        "phase4_trab_score",
        "phase4_trd_minus_trab",
    ]
    obs = df.loc[primary_rows, primary_obs_cols].copy()
    obs.index = df.loc[primary_rows, "obs_name"].astype(str).to_numpy()
    obs.index.name = None
    primary_scores = scores.reindex(primary_rows).fillna(0.0)
    primary_expr = expr.reindex(primary_rows).fillna(0.0)
    phenotype = assign_phenotypes(primary_scores)
    v_state = assign_v_gene_state(primary_expr)
    for col in primary_scores.columns:
        obs[col] = primary_scores[col].to_numpy(dtype=np.float32)
    for col in marker_config.get("v_gene_expression", []):
        if col in primary_expr.columns:
            obs[f"expr_{col}"] = primary_expr[col].to_numpy(dtype=np.float32)
    obs["gdt_phenotype_state"] = phenotype["gdt_phenotype_state"].to_numpy(dtype=object)
    obs["gdt_phenotype_confidence"] = phenotype["gdt_phenotype_confidence"].to_numpy(dtype=np.float32)
    obs["gdt_v_gene_state"] = v_state.to_numpy(dtype=object)
    obs["candidate_total_counts_for_marker_norm"] = total_counts[
        pd.Series(np.arange(candidate_rows.size), index=candidate_rows).loc[primary_rows].to_numpy(dtype=int)
    ]

    adata = ad.AnnData(X=x_primary, obs=compress_obs(obs), var=var)
    adata.layers["counts"] = adata.X.copy()
    adata.obs["n_counts"] = np.asarray(adata.layers["counts"].sum(axis=1)).ravel().astype(np.float32)
    adata.obs["n_genes_by_counts"] = np.diff(adata.layers["counts"].indptr).astype(np.int32)
    logging.info("Normalizing gdT-only object for marker visualization")
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)

    logging.info("Running gdT-only integration, UMAP, and Leiden clustering")
    try:
        if args.skip_scvi:
            raise RuntimeError("--skip-scvi requested")
        integration_info = run_integration(adata, args.n_hvgs, args.scvi_epochs, args.seed, args.skip_scvi)
    except Exception as exc:
        logging.warning("Using PCA fallback after integration exception: %s", exc)
        integration_info = {"integration_method": "PCA_fallback", "integration_error": repr(exc)}
        select_hvgs_no_tcr(adata, args.n_hvgs)
        hvg = adata[:, adata.var["highly_variable_no_tcr"].to_numpy(dtype=bool)].copy()
        sc.tl.pca(hvg, n_comps=50, svd_solver="arpack", random_state=args.seed)
        adata.obsm["X_scVI_gdT"] = hvg.obsm["X_pca"].astype(np.float32, copy=False)
        adata.obsm["X_scVI"] = adata.obsm["X_scVI_gdT"]
    embed_info = run_embedding_and_clustering(adata, args.seed)
    integration_info.update(embed_info)

    logging.info("Writing gdT-only H5AD")
    adata.write_h5ad(args.output_h5ad)

    logging.info("Writing phenotype, cluster, and statistical tables")
    phenotype_counts = (
        adata.obs.groupby(["gdt_phenotype_state"], observed=False)
        .agg(
            n_cells=("gdt_phenotype_state", "size"),
            n_sources=("source_gse_id", "nunique"),
            n_samples=("sample_unit_gdt_atlas", "nunique"),
            mean_confidence=("gdt_phenotype_confidence", "mean"),
        )
        .reset_index()
        .sort_values("n_cells", ascending=False)
    )
    phenotype_counts.to_csv(TABLE_DIR / "phenotype_counts.csv", index=False)
    phenotype_composition, phenotype_stats = phenotype_composition_stats(adata)
    score_cols = [c for c in adata.obs.columns if str(c).startswith("score_")]
    v_cols = [c for c in adata.obs.columns if str(c).startswith("expr_TR")]
    cluster_summary = summarize_clusters(adata, score_cols, v_cols)
    score_by_pheno = adata.obs.groupby("gdt_phenotype_state", observed=False)[score_cols].mean().reset_index()
    score_by_pheno.to_csv(TABLE_DIR / "marker_program_scores_by_phenotype.csv", index=False)
    v_by_pheno = adata.obs.groupby("gdt_phenotype_state", observed=False)[v_cols].mean().reset_index()
    v_by_pheno.to_csv(TABLE_DIR / "v_gene_expression_by_phenotype.csv", index=False)
    v_by_cluster = adata.obs.groupby("gdt_leiden", observed=False)[v_cols].mean().reset_index()
    v_by_cluster.to_csv(TABLE_DIR / "v_gene_expression_by_cluster.csv", index=False)
    adata.obs.reset_index(names="obs_name").to_csv(TABLE_DIR / "curated_gdt_cell_metadata_with_clusters.csv.gz", index=False)
    de = pseudobulk_de(adata, marker_config)
    de_top_path = TABLE_DIR / "phenotype_pseudobulk_de_top50_per_phenotype.csv"
    de_top = pd.read_csv(de_top_path) if de_top_path.exists() else pd.DataFrame()

    logging.info("Plotting report figures")
    figures: list[Path] = []
    figures.append(plot_curation_flow(curation_summary, FIGURE_DIR / "curation_flow_no_silver_addback.png"))
    for color, filename in [
        ("gdt_leiden", "gdt_only_umap_by_leiden.png"),
        ("gdt_phenotype_state", "gdt_only_umap_by_phenotype.png"),
        ("gdt_v_gene_state", "gdt_only_umap_by_trdv_state.png"),
        ("source_gse_id", "gdt_only_umap_by_source.png"),
        ("tissue_site_gdt_atlas", "gdt_only_umap_by_tissue.png"),
        ("disease_status_gdt_atlas", "gdt_only_umap_by_disease.png"),
        ("age_group_gdt_atlas", "gdt_only_umap_by_age.png"),
        ("gdt_atlas_evidence_tier", "gdt_only_umap_by_evidence_tier.png"),
    ]:
        figures.append(plot_umap_category(adata, color, FIGURE_DIR / filename, args.seed + len(figures)))
    figures.append(
        plot_score_heatmap(score_by_pheno, "gdt_phenotype_state", score_cols, FIGURE_DIR / "marker_program_scores_by_phenotype.png", "Marker programs by phenotype")
    )
    figures.append(
        plot_score_heatmap(v_by_pheno, "gdt_phenotype_state", v_cols, FIGURE_DIR / "v_gene_expression_by_phenotype.png", "TRD/TRG V/J/C RNA expression by phenotype")
    )
    figures.append(
        plot_score_heatmap(v_by_cluster, "gdt_leiden", v_cols, FIGURE_DIR / "v_gene_expression_by_cluster.png", "TRD/TRG V/J/C RNA expression by Leiden cluster")
    )
    figures.append(plot_abundance(sample_abundance, "tissue_site_gdt_atlas", FIGURE_DIR / "sample_level_gdt_abundance_by_tissue.png", "Sample-level gdT abundance by tissue"))
    figures.append(plot_abundance(sample_abundance, "disease_status_gdt_atlas", FIGURE_DIR / "sample_level_gdt_abundance_by_disease.png", "Sample-level gdT abundance by disease"))
    figures.append(plot_composition_heatmap(phenotype_composition, "tissue_site_gdt_atlas", FIGURE_DIR / "phenotype_composition_by_tissue.png"))
    figures.append(plot_composition_heatmap(phenotype_composition, "disease_status_gdt_atlas", FIGURE_DIR / "phenotype_composition_by_disease.png"))
    figures.append(plot_de_volcano(de, FIGURE_DIR / "phenotype_pseudobulk_de_volcano.png"))
    figures.append(plot_v_gene_umaps(adata, marker_config.get("primary_v_gene_panels", []), FIGURE_DIR / "v_gene_expression_umap_panels.png", args.seed))

    stat_after = args.input_h5ad.stat()
    manifest = {
        **info,
        **integration_info,
        "prediction_npz": str(args.prediction_npz),
        "prediction_key": PREDICTION_KEY,
        "promoted_predictions": int(pred.sum()),
        "removed_likely_fp": int(df["removed_fp"].sum()),
        "added_gold_fn": int(df["gold_fn_added"].sum()),
        "excluded_silver_fn": int(df["silver_fn_excluded"].sum()),
        "curated_primary_gdT": int(primary_rows.size),
        "derived_h5ad": str(args.output_h5ad),
        "source_size_after": stat_after.st_size,
        "source_mtime_ns_after": stat_after.st_mtime_ns,
        "source_h5ad_unchanged": bool(
            info["source_size_before"] == stat_after.st_size and info["source_mtime_ns_before"] == stat_after.st_mtime_ns
        ),
        "html": str(REPORT_HTML),
        "pdf": str(REPORT_PDF),
        "tables": str(TABLE_DIR),
        "figures": str(FIGURE_DIR),
    }
    if not manifest["source_h5ad_unchanged"]:
        raise RuntimeError("Source H5AD size/mtime changed during read-only workflow")
    write_report(
        manifest,
        curation_summary,
        curation_by_source,
        curation_by_tissue,
        abundance_stats,
        phenotype_counts,
        cluster_summary,
        de_top,
        figures,
    )
    manifest["pdf_rendered"] = bool(render_pdf(args.no_pdf))
    MANIFEST_JSON.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    logging.info("Wrote report: %s", REPORT_HTML)
    logging.info("Wrote PDF: %s", REPORT_PDF if manifest["pdf_rendered"] else "not rendered")


if __name__ == "__main__":
    main()
