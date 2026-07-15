#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""DEG reproducibility for TRD-high no-abTCR versus expanded negative cells.

Definitions:
- Start from the whole plus6 object and subset to phase4_trd_score > 0.1.
- High group: phase4_trd_minus_trab > 0.65 and no productive alpha-beta TCR
  evidence, represented by has_any_ab_tcr == False.
- Low group: phase4_trd_minus_trab < 0.35 and any productive alpha-beta TCR
  evidence, represented by has_any_ab_tcr == True. This includes single TRA,
  single TRB, or paired TRA/TRB evidence.
- Expanded negative group for datasets without TCRseq metadata: phase4_trab_score
  > -0.05 and no detected expression of any TRDV* or TRGV* gene.

The input H5AD is opened read-only. Expression is extracted as temporary
log1p(CP10K) from count-space X for the selected rows only. TCR genes are
excluded from all DEG outputs.
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


import html
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Any

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import pearsonr, spearmanr, t as student_t

from run_gdt_prediction_package_evaluation import (
    DEFAULT_INPUT_H5AD,
    FIGURE_DIR,
    LOG_DIR,
    STATIC_ASSET_DIR,
    STATIC_DIR,
    TABLE_DIR,
    clean_group_values,
    dataframe_to_html,
    dataframe_to_markdown,
    normalize_strings,
    read_bool_obs,
    read_float_obs,
    read_nonempty_string_mask,
    read_obs_column,
    read_string_dataset,
    required_columns_present,
)


OUTPUT_PREFIX = "trd_gt0p1_high_no_ab_vs_expanded_negative_no_tcr_genes"
TARGET_SUM = 10_000.0
TRD_SCORE_MIN = 0.1
TRD_MINUS_HIGH = 0.65
TRD_MINUS_LOW = 0.35
NO_TCRSEQ_NEGATIVE_TRAB_MIN = -0.05
GROUP_HIGH = "TRDminus_gt0p65_no_productive_abTCR"
GROUP_LOW = "expanded_negative"
MIN_GROUP_CELLS_FOR_SOURCE_DE = 50
DE_FDR_CUTOFF = 0.05
DE_MIN_ABS_DELTA = 0.10
DE_MIN_MAX_PCT_DETECTED = 0.05
STRICT_MIN_SOURCES_TESTED = 3
STRICT_MIN_SIGN_FRACTION = 0.75
STRICT_MIN_NOMINAL_SOURCES = 3
FIGURE_DPI = 300

LOG_PATH = LOG_DIR / f"{OUTPUT_PREFIX}_deg_reproducibility.log"
SUMMARY_MD = LOG_DIR / f"{OUTPUT_PREFIX}_deg_reproducibility_report.md"
REPORT_HTML = STATIC_DIR / f"{OUTPUT_PREFIX}_deg_reproducibility_report.html"
REPORT_PDF = STATIC_DIR / f"{OUTPUT_PREFIX}_deg_reproducibility_report.pdf"

CELL_COUNTS_SOURCE_CSV = TABLE_DIR / f"{OUTPUT_PREFIX}_cell_counts_by_source.csv"
CELL_COUNTS_TISSUE_CSV = TABLE_DIR / f"{OUTPUT_PREFIX}_cell_counts_by_tissue.csv"
SELECTION_RULE_COUNTS_CSV = TABLE_DIR / f"{OUTPUT_PREFIX}_selection_rule_counts.csv"
TCR_GENE_EXCLUSION_CSV = TABLE_DIR / f"{OUTPUT_PREFIX}_excluded_tcr_genes.csv"
GLOBAL_DEG_CSV = TABLE_DIR / f"{OUTPUT_PREFIX}_global_deg.csv"
SOURCE_DEG_CSV = TABLE_DIR / f"{OUTPUT_PREFIX}_within_source_deg.csv"
REPRO_CSV = TABLE_DIR / f"{OUTPUT_PREFIX}_deg_reproducibility.csv"
TOP_REPRO_CSV = TABLE_DIR / f"{OUTPUT_PREFIX}_top_reproducible_deg.csv"
SOURCE_PEARSON_CSV = TABLE_DIR / f"{OUTPUT_PREFIX}_source_delta_pearson_correlation.csv"
SOURCE_SPEARMAN_CSV = TABLE_DIR / f"{OUTPUT_PREFIX}_source_delta_spearman_correlation.csv"

VOLCANO_PNG = FIGURE_DIR / f"{OUTPUT_PREFIX}_global_deg_volcano.png"
SOURCE_CORR_PNG = FIGURE_DIR / f"{OUTPUT_PREFIX}_source_delta_correlation_heatmap.png"
TOP_HEATMAP_PNG = FIGURE_DIR / f"{OUTPUT_PREFIX}_top_reproducible_logfc_heatmap.png"
SOURCE_COUNTS_PNG = FIGURE_DIR / f"{OUTPUT_PREFIX}_cell_counts_by_source.png"


def setup_logging() -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    handlers = [logging.FileHandler(LOG_PATH, mode="w", encoding="utf-8"), logging.StreamHandler()]
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s", handlers=handlers, force=True)


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=np.float64)
    out = np.full(p.shape, np.nan, dtype=np.float64)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    valid_p = p[valid]
    order = np.argsort(valid_p)
    ranked = valid_p[order]
    n = ranked.size
    adjusted = ranked * n / np.arange(1, n + 1, dtype=np.float64)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    valid_out = np.empty_like(valid_p)
    valid_out[order] = adjusted
    out[valid] = valid_out
    return out


TCR_GENE_PREFIXES = (
    "TRAV",
    "TRAJ",
    "TRAC",
    "TRBV",
    "TRBJ",
    "TRBC",
    "TRGV",
    "TRGJ",
    "TRGC",
    "TRDV",
    "TRDD",
    "TRDJ",
    "TRDC",
)


def is_tcr_gene(gene: str) -> bool:
    return str(gene).upper().startswith(TCR_GENE_PREFIXES)


def read_bool_obs_if_present(handle: h5py.File, key: str, n_obs: int) -> np.ndarray:
    if key not in handle["obs"]:
        return np.zeros(n_obs, dtype=bool)
    return read_bool_obs(handle, key)


def read_nonempty_string_mask_if_present(handle: h5py.File, key: str, n_obs: int) -> np.ndarray:
    if key not in handle["obs"]:
        return np.zeros(n_obs, dtype=bool)
    return read_nonempty_string_mask(handle, key)


def build_no_tcrseq_source_mask(handle: h5py.File, source: np.ndarray) -> tuple[np.ndarray, list[str]]:
    n_obs = source.shape[0]
    any_tcr_metadata = (
        read_bool_obs_if_present(handle, "has_any_ab_tcr", n_obs)
        | read_bool_obs_if_present(handle, "has_any_gd_tcr", n_obs)
        | read_bool_obs_if_present(handle, "has_TRA_TRB_paired", n_obs)
        | read_bool_obs_if_present(handle, "has_TRG_TRD_paired", n_obs)
        | read_nonempty_string_mask_if_present(handle, "TRA_cdr3", n_obs)
        | read_nonempty_string_mask_if_present(handle, "TRB_cdr3", n_obs)
        | read_nonempty_string_mask_if_present(handle, "TRG_cdr3", n_obs)
        | read_nonempty_string_mask_if_present(handle, "TRD_cdr3", n_obs)
    )
    source_has_tcrseq = pd.Series(any_tcr_metadata).groupby(source).any()
    no_tcrseq_sources = sorted(source_has_tcrseq.index[~source_has_tcrseq.to_numpy(dtype=bool)].astype(str).tolist())
    return np.isin(source, no_tcrseq_sources), no_tcrseq_sources


def detect_gene_expression_for_rows(
    handle: h5py.File,
    selected_idx: np.ndarray,
    gene_idx: np.ndarray,
    *,
    label: str,
) -> np.ndarray:
    detected = np.zeros(selected_idx.size, dtype=bool)
    if selected_idx.size == 0 or gene_idx.size == 0:
        return detected
    gene_idx = np.asarray(np.sort(gene_idx), dtype=np.int32)
    x_group = handle["X"]
    indptr_ds = x_group["indptr"]
    indices_ds = x_group["indices"]
    data_ds = x_group["data"]
    for out_row, obs_idx in enumerate(selected_idx):
        start = int(indptr_ds[obs_idx])
        end = int(indptr_ds[obs_idx + 1])
        if end <= start:
            continue
        row_indices = indices_ds[start:end].astype(np.int32, copy=False)
        row_data = data_ds[start:end]
        gene_pos = np.isin(row_indices, gene_idx, assume_unique=False)
        detected[out_row] = bool(gene_pos.any() and np.any(row_data[gene_pos] > 0))
        if out_row and out_row % 20000 == 0:
            logging.info("Checked %s expression for %s / %s no-TCRseq negative candidates", label, f"{out_row:,}", f"{selected_idx.size:,}")
    return detected


def welch_deg_from_log_matrix(matrix: sparse.csr_matrix, group: np.ndarray, var_names: pd.Index) -> pd.DataFrame:
    """Return Welch-style DE stats for high group minus low group on log1p(CP10K)."""
    mask_high = group == GROUP_HIGH
    mask_low = group == GROUP_LOW
    n_high = int(mask_high.sum())
    n_low = int(mask_low.sum())
    if n_high == 0 or n_low == 0:
        raise ValueError(f"Both groups are required for DE; got high={n_high}, low={n_low}")

    x_high = matrix[mask_high]
    x_low = matrix[mask_low]
    sum_high = np.asarray(x_high.sum(axis=0)).ravel().astype(np.float64)
    sum_low = np.asarray(x_low.sum(axis=0)).ravel().astype(np.float64)
    sumsq_high = np.asarray(x_high.multiply(x_high).sum(axis=0)).ravel().astype(np.float64)
    sumsq_low = np.asarray(x_low.multiply(x_low).sum(axis=0)).ravel().astype(np.float64)
    mean_high = sum_high / n_high
    mean_low = sum_low / n_low
    var_high = np.maximum((sumsq_high - (sum_high * sum_high / n_high)) / max(n_high - 1, 1), 0.0)
    var_low = np.maximum((sumsq_low - (sum_low * sum_low / n_low)) / max(n_low - 1, 1), 0.0)
    se2 = var_high / n_high + var_low / n_low
    delta = mean_high - mean_low
    t_stat = np.divide(delta, np.sqrt(se2), out=np.zeros_like(delta), where=se2 > 0)
    numerator = se2 * se2
    denominator = (var_high * var_high) / (n_high * n_high * max(n_high - 1, 1)) + (
        var_low * var_low
    ) / (n_low * n_low * max(n_low - 1, 1))
    dof = np.divide(numerator, denominator, out=np.ones_like(numerator), where=denominator > 0)
    p_values = 2.0 * student_t.sf(np.abs(t_stat), dof)
    p_values[~np.isfinite(p_values)] = 1.0
    padj = bh_fdr(p_values)

    pct_high = np.asarray((x_high > 0).sum(axis=0)).ravel().astype(np.float64) / n_high
    pct_low = np.asarray((x_low > 0).sum(axis=0)).ravel().astype(np.float64) / n_low

    out = pd.DataFrame(
        {
            "gene": var_names.astype(str),
            "n_high_no_ab": n_high,
            "n_expanded_negative": n_low,
            "mean_log1p_cp10k_high_no_ab": mean_high,
            "mean_log1p_cp10k_expanded_negative": mean_low,
            "delta_log1p_cp10k_high_no_ab_minus_expanded_negative": delta,
            "pct_detected_high_no_ab": pct_high,
            "pct_detected_expanded_negative": pct_low,
            "t_stat": t_stat,
            "welch_df": dof,
            "p_value": p_values,
            "padj_bh": padj,
        }
    )
    out["max_pct_detected"] = out[["pct_detected_high_no_ab", "pct_detected_expanded_negative"]].max(axis=1)
    out["deg_pass"] = (
        (out["padj_bh"] < DE_FDR_CUTOFF)
        & (out["delta_log1p_cp10k_high_no_ab_minus_expanded_negative"].abs() >= DE_MIN_ABS_DELTA)
        & (out["max_pct_detected"] >= DE_MIN_MAX_PCT_DETECTED)
    )
    out = out.sort_values(
        ["padj_bh", "delta_log1p_cp10k_high_no_ab_minus_expanded_negative"],
        ascending=[True, False],
    ).reset_index(drop=True)
    out["rank_global"] = np.arange(1, out.shape[0] + 1, dtype=int)
    return out


def extract_selected_log_matrix(handle: h5py.File, selected_idx: np.ndarray) -> tuple[sparse.csr_matrix, pd.Index]:
    var_names = pd.Index(read_string_dataset(handle["var"]["_index"]), dtype="string")
    x_group = handle["X"]
    indptr_ds = x_group["indptr"]
    indices_ds = x_group["indices"]
    data_ds = x_group["data"]
    data_parts: list[np.ndarray] = []
    indices_parts: list[np.ndarray] = []
    indptr = np.zeros(selected_idx.size + 1, dtype=np.int64)

    for out_row, obs_idx in enumerate(selected_idx):
        start = int(indptr_ds[obs_idx])
        end = int(indptr_ds[obs_idx + 1])
        if end <= start:
            indptr[out_row + 1] = indptr[out_row]
            continue
        row_indices = indices_ds[start:end].astype(np.int32, copy=False)
        row_data = data_ds[start:end].astype(np.float32, copy=False)
        row_sum = float(np.sum(row_data, dtype=np.float64))
        if row_sum > 0:
            row_data = np.log1p(row_data * (TARGET_SUM / row_sum)).astype(np.float32, copy=False)
        else:
            row_data = np.zeros_like(row_data, dtype=np.float32)
        data_parts.append(row_data)
        indices_parts.append(row_indices)
        indptr[out_row + 1] = indptr[out_row] + row_data.size
        if out_row and out_row % 10000 == 0:
            logging.info("Extracted %s / %s selected rows", f"{out_row:,}", f"{selected_idx.size:,}")

    data = np.concatenate(data_parts) if data_parts else np.array([], dtype=np.float32)
    indices = np.concatenate(indices_parts) if indices_parts else np.array([], dtype=np.int32)
    matrix = sparse.csr_matrix((data, indices, indptr), shape=(selected_idx.size, len(var_names)), dtype=np.float32)
    tcr_gene_mask = np.asarray([is_tcr_gene(gene) for gene in var_names.astype(str)], dtype=bool)
    excluded_tcr_genes = pd.DataFrame({"gene": var_names[tcr_gene_mask].astype(str)})
    excluded_tcr_genes.to_csv(TCR_GENE_EXCLUSION_CSV, index=False)
    if tcr_gene_mask.any():
        matrix = matrix[:, ~tcr_gene_mask].tocsr()
        var_names = var_names[~tcr_gene_mask]
    logging.info("Excluded %s TCR genes from DEG matrix", f"{int(tcr_gene_mask.sum()):,}")
    return matrix, var_names


def build_subset_metadata(handle: h5py.File) -> pd.DataFrame:
    required_columns_present(handle)
    source = clean_group_values(read_obs_column(handle, "source_gse_id"))
    n_obs = source.shape[0]
    library_id = normalize_strings(read_obs_column(handle, "library_id"))
    sample_id = normalize_strings(read_obs_column(handle, "sample_id"))
    tissue_key = "tissue_corrected" if "tissue_corrected" in handle["obs"] else "tissue"
    tissue = clean_group_values(read_obs_column(handle, tissue_key))
    has_any_ab_tcr = read_bool_obs(handle, "has_any_ab_tcr")
    trd_score = read_float_obs(handle, "phase4_trd_score")
    trab_score = read_float_obs(handle, "phase4_trab_score")
    trd_minus_trab = read_float_obs(handle, "phase4_trd_minus_trab")
    no_tcrseq_source, no_tcrseq_sources = build_no_tcrseq_source_mask(handle, source)
    var_names = pd.Index(read_string_dataset(handle["var"]["_index"]), dtype="string")
    trdv_gene_idx = np.flatnonzero(np.asarray(var_names.astype(str).str.upper().str.startswith("TRDV"), dtype=bool))
    trgv_gene_idx = np.flatnonzero(np.asarray(var_names.astype(str).str.upper().str.startswith("TRGV"), dtype=bool))

    base = trd_score > TRD_SCORE_MIN
    preliminary_high = base & (trd_minus_trab > TRD_MINUS_HIGH) & (~has_any_ab_tcr)
    low_any_ab = base & (trd_minus_trab < TRD_MINUS_LOW) & has_any_ab_tcr
    no_tcrseq_negative_candidate = base & no_tcrseq_source & (trab_score > NO_TCRSEQ_NEGATIVE_TRAB_MIN)
    candidate_idx = np.flatnonzero(no_tcrseq_negative_candidate).astype(np.int64)
    candidate_has_trdv = detect_gene_expression_for_rows(handle, candidate_idx, trdv_gene_idx, label="TRDV")
    candidate_has_trgv = detect_gene_expression_for_rows(handle, candidate_idx, trgv_gene_idx, label="TRGV")
    no_tcrseq_negative_no_trdv_trgv = np.zeros(n_obs, dtype=bool)
    no_tcrseq_negative_no_trdv_trgv[candidate_idx] = ~(candidate_has_trdv | candidate_has_trgv)
    low_no_tcrseq = no_tcrseq_negative_candidate & no_tcrseq_negative_no_trdv_trgv
    low = low_any_ab | low_no_tcrseq
    overlap_high_low = preliminary_high & low
    high = preliminary_high & (~low)
    selected = high | low
    if not selected.any():
        raise RuntimeError("No cells passed the requested group definitions.")
    group = np.where(high[selected], GROUP_HIGH, GROUP_LOW)
    selection_rule = np.full(n_obs, "", dtype=object)
    selection_rule[high] = "high_trdminus_gt0p65_no_abTCR"
    selection_rule[low_any_ab] = "negative_trdminus_lt0p35_any_abTCR"
    selection_rule[low_no_tcrseq] = "negative_no_tcrseq_trab_gt_neg0p05_no_TRDV_TRGV_expression"

    metadata = pd.DataFrame(
        {
            "obs_index": np.flatnonzero(selected).astype(np.int64),
            "group": group,
            "selection_rule": selection_rule[selected],
            "source_gse_id": source[selected],
            "library_id": library_id[selected],
            "sample_id": sample_id[selected],
            "tissue": tissue[selected],
            "phase4_trd_score": trd_score[selected],
            "phase4_trab_score": trab_score[selected],
            "phase4_trd_minus_trab": trd_minus_trab[selected],
            "has_any_ab_tcr": has_any_ab_tcr[selected],
            "source_has_tcrseq_metadata": (~no_tcrseq_source[selected]),
        }
    )
    metadata.attrs["n_no_tcrseq_sources"] = len(no_tcrseq_sources)
    metadata.attrs["no_tcrseq_sources"] = ", ".join(no_tcrseq_sources)
    metadata.attrs["n_no_tcrseq_negative_candidates_before_TCRV_filter"] = int(no_tcrseq_negative_candidate.sum())
    metadata.attrs["n_no_tcrseq_negative_candidates_with_TRDV_expression"] = int(candidate_has_trdv.sum())
    metadata.attrs["n_no_tcrseq_negative_candidates_with_TRGV_expression"] = int(candidate_has_trgv.sum())
    metadata.attrs["n_no_tcrseq_negative_candidates_with_TRDV_or_TRGV_expression_removed"] = int((candidate_has_trdv | candidate_has_trgv).sum())
    metadata.attrs["n_high_precedence_removed_by_expanded_negative"] = int(overlap_high_low.sum())
    metadata.attrs["n_trdv_genes_used_for_negative_filter"] = int(trdv_gene_idx.size)
    metadata.attrs["n_trgv_genes_used_for_negative_filter"] = int(trgv_gene_idx.size)
    return metadata


def write_count_tables(metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    source_counts = (
        metadata.groupby(["source_gse_id", "group"], as_index=False)
        .agg(n_cells=("group", "size"))
        .sort_values(["source_gse_id", "group"])
    )
    source_counts.to_csv(CELL_COUNTS_SOURCE_CSV, index=False)
    tissue_counts = (
        metadata.groupby(["tissue", "group"], as_index=False)
        .agg(n_cells=("group", "size"))
        .sort_values(["tissue", "group"])
    )
    tissue_counts.to_csv(CELL_COUNTS_TISSUE_CSV, index=False)
    rule_counts = (
        metadata.groupby(["selection_rule", "group"], as_index=False)
        .agg(n_cells=("group", "size"))
        .sort_values(["group", "selection_rule"])
    )
    rule_counts.to_csv(SELECTION_RULE_COUNTS_CSV, index=False)
    return source_counts, tissue_counts


def compute_within_source_deg(matrix: sparse.csr_matrix, metadata: pd.DataFrame, var_names: pd.Index) -> pd.DataFrame:
    rows = []
    for source, idx in metadata.groupby("source_gse_id", sort=True).groups.items():
        idx_array = np.asarray(list(idx), dtype=np.int64)
        groups = metadata.loc[idx_array, "group"].to_numpy(dtype=object)
        n_high = int((groups == GROUP_HIGH).sum())
        n_low = int((groups == GROUP_LOW).sum())
        if n_high < MIN_GROUP_CELLS_FOR_SOURCE_DE or n_low < MIN_GROUP_CELLS_FOR_SOURCE_DE:
            logging.info("Skipping source %s for within-source DE: high=%s low=%s", source, n_high, n_low)
            continue
        source_deg = welch_deg_from_log_matrix(matrix[idx_array], groups, var_names)
        source_deg.insert(0, "source_gse_id", source)
        rows.append(source_deg)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, axis=0, ignore_index=True)


def build_reproducibility_table(global_deg: pd.DataFrame, source_deg: pd.DataFrame) -> pd.DataFrame:
    base_cols = [
        "gene",
        "rank_global",
        "delta_log1p_cp10k_high_no_ab_minus_expanded_negative",
        "padj_bh",
        "p_value",
        "mean_log1p_cp10k_high_no_ab",
        "mean_log1p_cp10k_expanded_negative",
        "pct_detected_high_no_ab",
        "pct_detected_expanded_negative",
        "deg_pass",
    ]
    repro = global_deg[base_cols].rename(
        columns={
            "delta_log1p_cp10k_high_no_ab_minus_expanded_negative": "global_delta_log1p",
            "padj_bh": "global_padj_bh",
            "p_value": "global_p_value",
            "mean_log1p_cp10k_high_no_ab": "global_mean_log1p_high_no_ab",
            "mean_log1p_cp10k_expanded_negative": "global_mean_log1p_expanded_negative",
            "pct_detected_high_no_ab": "global_pct_detected_high_no_ab",
            "pct_detected_expanded_negative": "global_pct_detected_expanded_negative",
            "deg_pass": "global_deg_pass",
        }
    )
    if source_deg.empty:
        repro["n_sources_tested"] = 0
        return repro

    source_delta = source_deg.pivot(index="gene", columns="source_gse_id", values="delta_log1p_cp10k_high_no_ab_minus_expanded_negative")
    source_p = source_deg.pivot(index="gene", columns="source_gse_id", values="p_value")
    source_padj = source_deg.pivot(index="gene", columns="source_gse_id", values="padj_bh")
    for source in source_delta.columns:
        repro = repro.merge(
            source_delta[[source]].rename(columns={source: f"{source}_delta_log1p"}),
            left_on="gene",
            right_index=True,
            how="left",
        )
        repro = repro.merge(
            source_p[[source]].rename(columns={source: f"{source}_p_value"}),
            left_on="gene",
            right_index=True,
            how="left",
        )
        repro = repro.merge(
            source_padj[[source]].rename(columns={source: f"{source}_padj_bh"}),
            left_on="gene",
            right_index=True,
            how="left",
        )

    source_cols = [f"{source}_delta_log1p" for source in source_delta.columns]
    deltas = repro[source_cols].to_numpy(dtype=np.float64)
    global_sign = np.sign(repro["global_delta_log1p"].to_numpy(dtype=np.float64))[:, None]
    valid = np.isfinite(deltas)
    same = valid & (np.sign(deltas) == global_sign) & (global_sign != 0)
    repro["n_sources_tested"] = valid.sum(axis=1)
    repro["n_sources_same_direction"] = same.sum(axis=1)
    repro["source_sign_concordance_fraction"] = np.divide(
        repro["n_sources_same_direction"],
        repro["n_sources_tested"],
        out=np.zeros(repro.shape[0], dtype=np.float64),
        where=repro["n_sources_tested"].to_numpy(dtype=int) > 0,
    )
    nominal_replicated = np.zeros(repro.shape[0], dtype=int)
    fdr_replicated = np.zeros(repro.shape[0], dtype=int)
    for source in source_delta.columns:
        same_dir = np.sign(repro[f"{source}_delta_log1p"].to_numpy(dtype=np.float64)) == np.sign(
            repro["global_delta_log1p"].to_numpy(dtype=np.float64)
        )
        nominal_replicated += ((repro[f"{source}_p_value"].to_numpy(dtype=np.float64) < 0.05) & same_dir).astype(int)
        fdr_replicated += ((repro[f"{source}_padj_bh"].to_numpy(dtype=np.float64) < 0.05) & same_dir).astype(int)
    repro["n_sources_nominal_p_lt_0p05_same_direction"] = nominal_replicated
    repro["n_sources_fdr_lt_0p05_same_direction"] = fdr_replicated
    repro["reproducible_deg_strict"] = (
        repro["global_deg_pass"]
        & (repro["n_sources_tested"] >= STRICT_MIN_SOURCES_TESTED)
        & (repro["source_sign_concordance_fraction"] >= STRICT_MIN_SIGN_FRACTION)
        & (repro["n_sources_nominal_p_lt_0p05_same_direction"] >= STRICT_MIN_NOMINAL_SOURCES)
    )
    repro = repro.sort_values(
        [
            "reproducible_deg_strict",
            "n_sources_nominal_p_lt_0p05_same_direction",
            "source_sign_concordance_fraction",
            "global_padj_bh",
            "global_delta_log1p",
        ],
        ascending=[False, False, False, True, False],
    ).reset_index(drop=True)
    return repro


def compute_source_correlations(source_deg: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    if source_deg.empty:
        return pd.DataFrame(), pd.DataFrame(), {}
    source_delta = source_deg.pivot(index="gene", columns="source_gse_id", values="delta_log1p_cp10k_high_no_ab_minus_expanded_negative")
    pearson = source_delta.corr(method="pearson")
    spearman = source_delta.corr(method="spearman")
    pearson.to_csv(SOURCE_PEARSON_CSV)
    spearman.to_csv(SOURCE_SPEARMAN_CSV)

    def offdiag_values(frame: pd.DataFrame) -> np.ndarray:
        if frame.shape[0] < 2:
            return np.array([], dtype=np.float64)
        mask = ~np.eye(frame.shape[0], dtype=bool)
        return frame.to_numpy(dtype=np.float64)[mask]

    pearson_vals = offdiag_values(pearson)
    spearman_vals = offdiag_values(spearman)
    stats: dict[str, Any] = {
        "n_sources_with_both_groups": int(source_delta.shape[1]),
        "n_genes_compared": int(source_delta.shape[0]),
        "median_pairwise_pearson_r": float(np.nanmedian(pearson_vals)) if pearson_vals.size else None,
        "median_pairwise_spearman_r": float(np.nanmedian(spearman_vals)) if spearman_vals.size else None,
        "min_pairwise_pearson_r": float(np.nanmin(pearson_vals)) if pearson_vals.size else None,
        "max_pairwise_pearson_r": float(np.nanmax(pearson_vals)) if pearson_vals.size else None,
    }
    return pearson, spearman, stats


def plot_volcano(global_deg: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 6.0), constrained_layout=True)
    x = global_deg["delta_log1p_cp10k_high_no_ab_minus_expanded_negative"].to_numpy(dtype=np.float64)
    y = -np.log10(np.maximum(global_deg["padj_bh"].to_numpy(dtype=np.float64), 1e-300))
    sig = global_deg["deg_pass"].to_numpy(dtype=bool)
    ax.scatter(x[~sig], y[~sig], s=4, c="#b8b8b8", alpha=0.45, linewidths=0, rasterized=True)
    ax.scatter(x[sig], y[sig], s=6, c="#b51f2a", alpha=0.7, linewidths=0, rasterized=True)
    ax.axvline(DE_MIN_ABS_DELTA, color="black", linestyle="--", lw=0.9)
    ax.axvline(-DE_MIN_ABS_DELTA, color="black", linestyle="--", lw=0.9)
    ax.axhline(-np.log10(DE_FDR_CUTOFF), color="black", linestyle=":", lw=0.9)
    top = global_deg.loc[sig].head(15)
    for _, row in top.iterrows():
        ax.text(
            row["delta_log1p_cp10k_high_no_ab_minus_expanded_negative"],
            -np.log10(max(row["padj_bh"], 1e-300)),
            row["gene"],
            fontsize=7,
        )
    ax.set_xlabel("Mean log1p(CP10K) difference: high no-abTCR - expanded negative")
    ax.set_ylabel("-log10(BH FDR)")
    ax.set_title("TRDscore>0.1 high no-abTCR vs expanded negative global DE")
    fig.savefig(VOLCANO_PNG, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_source_correlation_heatmap(pearson: pd.DataFrame) -> None:
    if pearson.empty:
        return
    fig, ax = plt.subplots(figsize=(max(6, 0.55 * pearson.shape[1] + 2), max(5.5, 0.5 * pearson.shape[0] + 2)), constrained_layout=True)
    im = ax.imshow(pearson.to_numpy(dtype=np.float64), cmap="vlag" if "vlag" in plt.colormaps() else "coolwarm", vmin=-1, vmax=1)
    ax.set_xticks(np.arange(pearson.shape[1]), labels=pearson.columns.tolist(), rotation=45, ha="right")
    ax.set_yticks(np.arange(pearson.shape[0]), labels=pearson.index.tolist())
    ax.set_title("Pairwise source correlation of gene-level delta log1p")
    for i in range(pearson.shape[0]):
        for j in range(pearson.shape[1]):
            ax.text(j, i, f"{pearson.iloc[i, j]:.2f}", ha="center", va="center", fontsize=7)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Pearson r")
    fig.savefig(SOURCE_CORR_PNG, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_top_heatmap(repro: pd.DataFrame) -> None:
    source_cols = [column for column in repro.columns if column.endswith("_delta_log1p") and column != "global_delta_log1p"]
    heat_cols = ["global_delta_log1p", *source_cols]
    top = repro.loc[repro["reproducible_deg_strict"]].head(35)
    if top.empty:
        top = repro.loc[repro["global_deg_pass"]].head(35)
    if top.empty:
        return
    mat = top.set_index("gene")[heat_cols]
    fig, ax = plt.subplots(figsize=(max(7, 0.75 * len(heat_cols)), max(6.5, 0.23 * mat.shape[0] + 1.5)), constrained_layout=True)
    vmax = float(np.nanquantile(np.abs(mat.to_numpy(dtype=np.float64)), 0.95))
    vmax = max(vmax, 0.1)
    im = ax.imshow(mat.to_numpy(dtype=np.float64), aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(len(heat_cols)), labels=[col.replace("_delta_log1p", "") for col in heat_cols], rotation=45, ha="right")
    ax.set_yticks(np.arange(mat.shape[0]), labels=mat.index.tolist())
    ax.set_title("Top reproducible/global DEG logFC")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Delta log1p(CP10K)")
    fig.savefig(TOP_HEATMAP_PNG, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_source_counts(source_counts: pd.DataFrame) -> None:
    pivot = source_counts.pivot(index="source_gse_id", columns="group", values="n_cells").fillna(0)
    pivot["total"] = pivot.sum(axis=1)
    pivot = pivot.sort_values("total", ascending=False).head(30)
    plot_cols = [column for column in [GROUP_HIGH, GROUP_LOW] if column in pivot.columns]
    fig, ax = plt.subplots(figsize=(9.5, max(4.8, 0.27 * pivot.shape[0] + 1.2)), constrained_layout=True)
    y = np.arange(pivot.shape[0])
    left = np.zeros(pivot.shape[0], dtype=np.float64)
    colors = {GROUP_HIGH: "#b51f2a", GROUP_LOW: "#2d6f9f"}
    labels = {GROUP_HIGH: "TRD-TRAB > 0.65, no abTCR", GROUP_LOW: "Expanded negative"}
    for col in plot_cols:
        values = pivot[col].to_numpy(dtype=np.float64)
        ax.barh(y, values, left=left, color=colors.get(col, "#777777"), label=labels.get(col, col))
        left += values
    ax.set_yticks(y, labels=pivot.index.tolist())
    ax.invert_yaxis()
    ax.set_xlabel("Cells")
    ax.set_title("Selected TRDscore>0.1 cells by source")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(SOURCE_COUNTS_PNG, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def copy_static_assets() -> list[str]:
    STATIC_ASSET_DIR.mkdir(parents=True, exist_ok=True)
    copied: list[str] = []
    for path in [VOLCANO_PNG, SOURCE_CORR_PNG, TOP_HEATMAP_PNG, SOURCE_COUNTS_PNG]:
        if path.exists():
            shutil.copyfile(path, STATIC_ASSET_DIR / path.name)
            copied.append(path.name)
    return copied


def write_markdown_report(
    input_h5ad: Path,
    n_obs: int,
    metadata: pd.DataFrame,
    source_counts: pd.DataFrame,
    tissue_counts: pd.DataFrame,
    global_deg: pd.DataFrame,
    source_deg: pd.DataFrame,
    repro: pd.DataFrame,
    corr_stats: dict[str, Any],
) -> None:
    n_high = int((metadata["group"] == GROUP_HIGH).sum())
    n_low = int((metadata["group"] == GROUP_LOW).sum())
    global_pass_n = int(global_deg["deg_pass"].sum())
    strict_n = int(repro["reproducible_deg_strict"].sum()) if "reproducible_deg_strict" in repro else 0
    sources = sorted(source_deg["source_gse_id"].unique().tolist()) if not source_deg.empty else []
    top_strict = repro.loc[repro["reproducible_deg_strict"]].head(30)
    top_global = global_deg.head(30)
    top_high = (
        repro.loc[repro["reproducible_deg_strict"]]
        .sort_values("global_delta_log1p", ascending=False)
        .head(15)
        [["gene", "global_delta_log1p", "source_sign_concordance_fraction", "n_sources_nominal_p_lt_0p05_same_direction"]]
    )
    top_low = (
        repro.loc[repro["reproducible_deg_strict"]]
        .sort_values("global_delta_log1p", ascending=True)
        .head(15)
        [["gene", "global_delta_log1p", "source_sign_concordance_fraction", "n_sources_nominal_p_lt_0p05_same_direction"]]
    )
    rule_counts = (
        metadata.groupby(["selection_rule", "group"], as_index=False)
        .agg(n_cells=("group", "size"))
        .sort_values(["group", "selection_rule"])
    )
    lines = [
        "# TRDscore>0.1 high no-abTCR vs expanded negative DEG reproducibility",
        "",
        "## Definition",
        "",
        f"- Input H5AD: `{input_h5ad}`",
        f"- Whole plus6 cells: `{n_obs:,}`",
        f"- Initial subset: `phase4_trd_score > {TRD_SCORE_MIN}`.",
        f"- High group: `phase4_trd_minus_trab > {TRD_MINUS_HIGH}` and `has_any_ab_tcr == False`.",
        f"- Expanded negative comparator includes TCRseq cells with `phase4_trd_minus_trab < {TRD_MINUS_LOW}` and `has_any_ab_tcr == True`, including single TRA, single TRB, or paired TRA/TRB.",
        f"- Expanded negative comparator also includes no-TCRseq-source cells with `phase4_trab_score > {NO_TCRSEQ_NEGATIVE_TRAB_MIN}` and no detected expression of any `TRDV*` or `TRGV*` gene.",
        "- Expanded-negative assignment has precedence over high-score assignment if a no-TCRseq cell satisfies both masks.",
        "- TCR genes were excluded from the DEG matrix and DEG lists.",
        "- Expression was calculated as temporary `log1p(CP10K)` from count-space `X`; the H5AD was not modified.",
        "",
        "## Summary",
        "",
        f"- Selected cells: `{metadata.shape[0]:,}`",
        f"- High no-abTCR cells: `{n_high:,}`",
        f"- Expanded negative cells: `{n_low:,}`",
        f"- Global DEG-pass genes: `{global_pass_n:,}` using FDR < `{DE_FDR_CUTOFF}`, abs delta >= `{DE_MIN_ABS_DELTA}`, max detected fraction >= `{DE_MIN_MAX_PCT_DETECTED}`.",
        f"- Strict reproducible genes: `{strict_n:,}` using global DEG pass, >= `{STRICT_MIN_SOURCES_TESTED}` sources tested, >= `{STRICT_MIN_SIGN_FRACTION:.0%}` same direction, and >= `{STRICT_MIN_NOMINAL_SOURCES}` same-direction nominal source tests.",
        f"- Sources used for within-source reproducibility: `{', '.join(sources) if sources else 'none'}`.",
        f"- Median pairwise source Pearson r: `{corr_stats.get('median_pairwise_pearson_r')}`",
        f"- Median pairwise source Spearman r: `{corr_stats.get('median_pairwise_spearman_r')}`",
        f"- No-TCRseq sources detected: `{metadata.attrs.get('n_no_tcrseq_sources')}`.",
        f"- No-TCRseq negative candidates before TCRV-expression filter: `{metadata.attrs.get('n_no_tcrseq_negative_candidates_before_TCRV_filter')}`.",
        f"- No-TCRseq negative candidates with any TRDV expression: `{metadata.attrs.get('n_no_tcrseq_negative_candidates_with_TRDV_expression')}`.",
        f"- No-TCRseq negative candidates with any TRGV expression: `{metadata.attrs.get('n_no_tcrseq_negative_candidates_with_TRGV_expression')}`.",
        f"- No-TCRseq negative candidates removed because any TRDV or TRGV gene was detected: `{metadata.attrs.get('n_no_tcrseq_negative_candidates_with_TRDV_or_TRGV_expression_removed')}`.",
        f"- High-score cells reassigned/excluded by expanded-negative precedence: `{metadata.attrs.get('n_high_precedence_removed_by_expanded_negative')}`.",
        f"- TRDV genes used for no-expression filter: `{metadata.attrs.get('n_trdv_genes_used_for_negative_filter')}`.",
        f"- TRGV genes used for no-expression filter: `{metadata.attrs.get('n_trgv_genes_used_for_negative_filter')}`.",
        "",
        "## Top Strict Reproducible High-Group-Enriched Genes",
        "",
        dataframe_to_markdown(top_high),
        "",
        "## Top Strict Reproducible Expanded-Negative-Enriched Genes",
        "",
        dataframe_to_markdown(top_low),
        "",
        "## Selection Rule Counts",
        "",
        dataframe_to_markdown(rule_counts),
        "",
        "## Cell Counts By Source",
        "",
        dataframe_to_markdown(source_counts, max_rows=40),
        "",
        "## Cell Counts By Tissue",
        "",
        dataframe_to_markdown(tissue_counts, max_rows=40),
        "",
        "## Top Strict Reproducible DEGs",
        "",
        dataframe_to_markdown(top_strict, max_rows=30),
        "",
        "## Top Global DEGs",
        "",
        dataframe_to_markdown(top_global, max_rows=30),
        "",
        "## Output Files",
        "",
        f"- Global DEG table: `{GLOBAL_DEG_CSV}`",
        f"- Within-source DEG table: `{SOURCE_DEG_CSV}`",
        f"- Reproducibility table: `{REPRO_CSV}`",
        f"- Excluded TCR genes: `{TCR_GENE_EXCLUSION_CSV}`",
        f"- Static HTML report: `{REPORT_HTML}`",
        f"- Static PDF report: `{REPORT_PDF}`",
        "",
    ]
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def image_section(title: str, filename: str) -> str:
    return (
        f"<section class='figure-card'><h3>{html.escape(title)}</h3>"
        f"<img src='assets/{html.escape(filename)}' alt='{html.escape(title)}'></section>"
    )


def write_html_report(
    input_h5ad: Path,
    n_obs: int,
    metadata: pd.DataFrame,
    source_counts: pd.DataFrame,
    tissue_counts: pd.DataFrame,
    global_deg: pd.DataFrame,
    repro: pd.DataFrame,
    corr_stats: dict[str, Any],
    copied_assets: list[str],
) -> None:
    n_high = int((metadata["group"] == GROUP_HIGH).sum())
    n_low = int((metadata["group"] == GROUP_LOW).sum())
    strict_n = int(repro["reproducible_deg_strict"].sum()) if "reproducible_deg_strict" in repro else 0
    global_pass_n = int(global_deg["deg_pass"].sum())
    top_high = (
        repro.loc[repro["reproducible_deg_strict"]]
        .sort_values("global_delta_log1p", ascending=False)
        .head(20)
        [["gene", "global_delta_log1p", "source_sign_concordance_fraction", "n_sources_nominal_p_lt_0p05_same_direction"]]
    )
    top_low = (
        repro.loc[repro["reproducible_deg_strict"]]
        .sort_values("global_delta_log1p", ascending=True)
        .head(20)
        [["gene", "global_delta_log1p", "source_sign_concordance_fraction", "n_sources_nominal_p_lt_0p05_same_direction"]]
    )
    rule_counts = (
        metadata.groupby(["selection_rule", "group"], as_index=False)
        .agg(n_cells=("group", "size"))
        .sort_values(["group", "selection_rule"])
    )
    titles = {
        VOLCANO_PNG.name: "Global DEG volcano",
        SOURCE_CORR_PNG.name: "Source delta correlation heatmap",
        TOP_HEATMAP_PNG.name: "Top reproducible DEG logFC heatmap",
        SOURCE_COUNTS_PNG.name: "Selected cell counts by source",
    }
    figure_html = "\n".join(image_section(titles.get(name, name), name) for name in copied_assets)
    css = """
    :root{--bg:#f6f7f8;--paper:#ffffff;--ink:#1e252b;--muted:#5d6975;--line:#d9dee3;--accent:#b51f2a;--blue:#2d6f9f}
    *{box-sizing:border-box} body{margin:0;background:var(--bg);color:var(--ink);font-family:Arial,Helvetica,sans-serif;line-height:1.55}
    .page{width:min(1240px,calc(100vw - 32px));margin:22px auto 44px}
    .hero,.section{background:var(--paper);border:1px solid var(--line);padding:24px;margin-top:18px}.hero{border-top:6px solid var(--accent)}
    h1{margin:0 0 10px;font-size:32px;line-height:1.15;letter-spacing:0} h2{margin:0 0 12px;font-size:23px} h3{margin:16px 0 8px;font-size:16px}
    p{margin:0 0 12px;color:var(--muted)} code{background:#eef1f4;padding:1px 5px;border-radius:4px}
    .metric-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:10px;margin-top:16px}
    .metric{border:1px solid var(--line);background:#fbfcfd;padding:11px 13px}.metric .value{display:block;font-size:23px;font-weight:700}.metric .label{display:block;font-size:12px;text-transform:uppercase;color:var(--muted)}
    .figure-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(430px,1fr));gap:16px}.figure-card{border:1px solid var(--line);background:#fbfcfd;padding:12px}.figure-card img{width:100%;height:auto;border:1px solid var(--line);background:white}
    table{border-collapse:collapse;width:100%;font-size:12px;margin:8px 0 16px} th,td{border:1px solid var(--line);padding:5px 7px;text-align:left;vertical-align:top} th{background:#eef1f4}.empty-table{font-style:italic}
    """
    html_parts = [
        "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'>",
        "<title>TRDscore high no-abTCR DEG reproducibility</title>",
        f"<style>{css}</style></head><body><main class='page'>",
        "<section class='hero'>",
        "<h1>TRDscore&gt;0.1 High No-abTCR vs Expanded Negative DEG Reproducibility</h1>",
        f"<p>Read-only whole-plus6 analysis from <code>{html.escape(str(input_h5ad))}</code>.</p>",
        "<div class='metric-grid'>",
        f"<div class='metric'><span class='value'>{n_obs:,}</span><span class='label'>Plus6 cells</span></div>",
        f"<div class='metric'><span class='value'>{metadata.shape[0]:,}</span><span class='label'>Selected cells</span></div>",
        f"<div class='metric'><span class='value'>{n_high:,}</span><span class='label'>High no-abTCR</span></div>",
        f"<div class='metric'><span class='value'>{n_low:,}</span><span class='label'>Expanded negative</span></div>",
        f"<div class='metric'><span class='value'>{global_pass_n:,}</span><span class='label'>Global DEG-pass genes</span></div>",
        f"<div class='metric'><span class='value'>{strict_n:,}</span><span class='label'>Strict reproducible genes</span></div>",
        "</div></section>",
        "<section class='section'><h2>Definitions</h2>",
        f"<p>Initial subset: <code>phase4_trd_score &gt; {TRD_SCORE_MIN}</code>. High group: <code>phase4_trd_minus_trab &gt; {TRD_MINUS_HIGH}</code> and <code>has_any_ab_tcr == False</code>. Expanded negative includes TCRseq cells with <code>phase4_trd_minus_trab &lt; {TRD_MINUS_LOW}</code> and <code>has_any_ab_tcr == True</code>, plus no-TCRseq-source cells with <code>phase4_trab_score &gt; {NO_TCRSEQ_NEGATIVE_TRAB_MIN}</code> and no detected <code>TRDV*</code> or <code>TRGV*</code> expression.</p>",
        "<p>Expanded-negative assignment has precedence over high-score assignment for no-TCRseq cells satisfying both masks. TCR genes were excluded from the DEG matrix and DEG lists.</p>",
        "<p>Expression was calculated as temporary <code>log1p(CP10K)</code> from count-space <code>X</code>; no H5AD was modified.</p>",
        f"<p>Strict reproducibility requires global DEG pass, at least {STRICT_MIN_SOURCES_TESTED} tested source datasets, at least {STRICT_MIN_SIGN_FRACTION:.0%} same-direction source effects, and at least {STRICT_MIN_NOMINAL_SOURCES} same-direction nominal source tests.</p>",
        f"<p>Median pairwise source Pearson r: <code>{corr_stats.get('median_pairwise_pearson_r')}</code>; median pairwise source Spearman r: <code>{corr_stats.get('median_pairwise_spearman_r')}</code>.</p>",
        f"<p>No-TCRseq negative candidates before TCRV-expression filtering: <code>{metadata.attrs.get('n_no_tcrseq_negative_candidates_before_TCRV_filter')}</code>; removed for any TRDV or TRGV expression: <code>{metadata.attrs.get('n_no_tcrseq_negative_candidates_with_TRDV_or_TRGV_expression_removed')}</code>; high-score cells reassigned/excluded by negative precedence: <code>{metadata.attrs.get('n_high_precedence_removed_by_expanded_negative')}</code>.</p>",
        "</section>",
        f"<section class='section'><h2>Figures</h2><div class='figure-grid'>{figure_html}</div></section>",
        "<section class='section'><h2>Top Strict Reproducible High-Group-Enriched Genes</h2>",
        dataframe_to_html(top_high),
        "</section>",
        "<section class='section'><h2>Top Strict Reproducible Expanded-Negative-Enriched Genes</h2>",
        dataframe_to_html(top_low),
        "</section>",
        "<section class='section'><h2>Selection Rule Counts</h2>",
        dataframe_to_html(rule_counts),
        "</section>",
        "<section class='section'><h2>Cell Counts</h2><h3>By Source</h3>",
        dataframe_to_html(source_counts, max_rows=60),
        "<h3>By Tissue</h3>",
        dataframe_to_html(tissue_counts, max_rows=60),
        "</section>",
        "<section class='section'><h2>Top Global DEGs</h2>",
        dataframe_to_html(global_deg.head(50), max_rows=50),
        "</section>",
        "</main></body></html>",
    ]
    REPORT_HTML.write_text("\n".join(html_parts), encoding="utf-8")


def render_pdf() -> None:
    chrome_candidates = [
        Path("/usr/bin/google-chrome"),
        Path("/usr/bin/google-chrome-stable"),
        shutil.which("google-chrome"),
        shutil.which("google-chrome-stable"),
    ]
    chrome = next((Path(path) for path in chrome_candidates if path and Path(path).exists()), None)
    if chrome is None:
        logging.warning("Skipping PDF export because google-chrome was not found")
        return
    subprocess.run(
        [
            str(chrome),
            "--headless",
            "--disable-gpu",
            "--no-sandbox",
            "--print-to-pdf-no-header",
            f"--print-to-pdf={REPORT_PDF}",
            str(REPORT_HTML),
        ],
        check=True,
    )


def main() -> None:
    setup_logging()
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    STATIC_DIR.mkdir(parents=True, exist_ok=True)
    input_h5ad = DEFAULT_INPUT_H5AD.resolve()
    stat_before = input_h5ad.stat()

    with h5py.File(input_h5ad, "r") as handle:
        n_obs = int(handle["obs"]["_index"].shape[0])
        logging.info("Building requested TRDscore>0.1 high/expanded-negative metadata from %s cells", f"{n_obs:,}")
        metadata = build_subset_metadata(handle)
        source_counts, tissue_counts = write_count_tables(metadata)
        logging.info(
            "Selected %s cells: high no-abTCR=%s, expanded negative=%s",
            f"{metadata.shape[0]:,}",
            f"{int((metadata['group'] == GROUP_HIGH).sum()):,}",
            f"{int((metadata['group'] == GROUP_LOW).sum()):,}",
        )

        logging.info("Extracting log1p(CP10K) matrix for selected rows")
        matrix, var_names = extract_selected_log_matrix(handle, metadata["obs_index"].to_numpy(dtype=np.int64))

    logging.info("Computing global high no-abTCR vs expanded negative DE")
    global_deg = welch_deg_from_log_matrix(matrix, metadata["group"].to_numpy(dtype=object), var_names)
    global_deg.to_csv(GLOBAL_DEG_CSV, index=False)

    logging.info("Computing within-source DE for reproducibility")
    source_deg = compute_within_source_deg(matrix, metadata, var_names)
    source_deg.to_csv(SOURCE_DEG_CSV, index=False)

    logging.info("Building reproducibility table")
    repro = build_reproducibility_table(global_deg, source_deg)
    repro.to_csv(REPRO_CSV, index=False)
    top_repro = repro.loc[repro["reproducible_deg_strict"]].head(300)
    if top_repro.empty:
        top_repro = repro.loc[repro["global_deg_pass"]].head(300)
    top_repro.to_csv(TOP_REPRO_CSV, index=False)

    logging.info("Computing source delta correlations")
    pearson, _spearman, corr_stats = compute_source_correlations(source_deg)

    logging.info("Rendering figures")
    plot_volcano(global_deg)
    plot_source_correlation_heatmap(pearson)
    plot_top_heatmap(repro)
    plot_source_counts(source_counts)

    copied_assets = copy_static_assets()
    write_markdown_report(input_h5ad, n_obs, metadata, source_counts, tissue_counts, global_deg, source_deg, repro, corr_stats)
    write_html_report(input_h5ad, n_obs, metadata, source_counts, tissue_counts, global_deg, repro, corr_stats, copied_assets)
    render_pdf()

    stat_after = input_h5ad.stat()
    if (stat_before.st_size != stat_after.st_size) or (stat_before.st_mtime_ns != stat_after.st_mtime_ns):
        raise RuntimeError("Input H5AD changed during read-only DEG analysis.")

    logging.info("Saved markdown report: %s", SUMMARY_MD)
    logging.info("Saved HTML report: %s", REPORT_HTML)
    if REPORT_PDF.exists():
        logging.info("Saved PDF report: %s", REPORT_PDF)
    logging.info("Saved global DEG table: %s", GLOBAL_DEG_CSV)
    logging.info("Saved reproducibility table: %s", REPRO_CSV)


if __name__ == "__main__":
    main()
