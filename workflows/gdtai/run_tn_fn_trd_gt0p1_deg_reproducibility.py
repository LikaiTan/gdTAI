#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""DEG and cross-dataset reproducibility for TRD-high TN versus FN cells.

Definitions:
- TN: corrected abT_gold cells predicted negative by the current conservative
  gdT classifier.
- FN: corrected gdT_gold cells predicted negative by the same classifier.
- Analysis subset: TN/FN cells with phase4_trd_score > 0.1.

The input H5AD is opened read-only. Expression is extracted as temporary
log1p(CP10K) from count-space X for the selected rows only.
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


import json
import logging
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
    TABLE_DIR,
    build_corrected_tcr_evidence,
    build_sublibrary_labels,
    clean_group_values,
    dataframe_to_markdown,
    make_truth_labels,
    normalize_strings,
    read_bool_obs,
    read_float_obs,
    read_nonempty_string_mask,
    read_obs_column,
    read_string_dataset,
    required_columns_present,
)


OUTPUT_PREFIX = "tn_fn_trd_gt0p1"
TARGET_SUM = 10_000.0
TRD_SCORE_MIN = 0.1
MIN_GROUP_CELLS_FOR_SOURCE_DE = 50
DE_FDR_CUTOFF = 0.05
DE_MIN_ABS_DELTA = 0.10
DE_MIN_MAX_PCT_DETECTED = 0.05
FIGURE_DPI = 300

LOG_PATH = LOG_DIR / f"{OUTPUT_PREFIX}_deg_reproducibility.log"
SUMMARY_MD = LOG_DIR / f"{OUTPUT_PREFIX}_deg_reproducibility_summary.md"
COUNTS_CSV = TABLE_DIR / f"{OUTPUT_PREFIX}_cell_counts_by_source.csv"
GLOBAL_DEG_CSV = TABLE_DIR / f"{OUTPUT_PREFIX}_global_deg.csv"
SOURCE_DEG_CSV = TABLE_DIR / f"{OUTPUT_PREFIX}_within_source_deg.csv"
REPRO_CSV = TABLE_DIR / f"{OUTPUT_PREFIX}_deg_reproducibility.csv"
TOP_REPRO_CSV = TABLE_DIR / f"{OUTPUT_PREFIX}_top_reproducible_deg.csv"
VOLCANO_PNG = FIGURE_DIR / f"{OUTPUT_PREFIX}_global_deg_volcano.png"
REPRO_SCATTER_PNG = FIGURE_DIR / f"{OUTPUT_PREFIX}_source_logfc_reproducibility.png"
TOP_HEATMAP_PNG = FIGURE_DIR / f"{OUTPUT_PREFIX}_top_reproducible_logfc_heatmap.png"


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


def welch_deg_from_log_matrix(
    matrix: sparse.csr_matrix,
    group: np.ndarray,
    var_names: pd.Index,
    *,
    label_a: str = "FN",
    label_b: str = "TN",
) -> pd.DataFrame:
    """Return Welch-style DE stats for label_a minus label_b on log1p(CP10K)."""
    mask_a = group == label_a
    mask_b = group == label_b
    n_a = int(mask_a.sum())
    n_b = int(mask_b.sum())
    if n_a == 0 or n_b == 0:
        raise ValueError(f"Both groups are required for DE; got {label_a}={n_a}, {label_b}={n_b}")

    x_a = matrix[mask_a]
    x_b = matrix[mask_b]
    sum_a = np.asarray(x_a.sum(axis=0)).ravel().astype(np.float64)
    sum_b = np.asarray(x_b.sum(axis=0)).ravel().astype(np.float64)
    sumsq_a = np.asarray(x_a.multiply(x_a).sum(axis=0)).ravel().astype(np.float64)
    sumsq_b = np.asarray(x_b.multiply(x_b).sum(axis=0)).ravel().astype(np.float64)
    mean_a = sum_a / n_a
    mean_b = sum_b / n_b
    var_a = np.maximum((sumsq_a - (sum_a * sum_a / n_a)) / max(n_a - 1, 1), 0.0)
    var_b = np.maximum((sumsq_b - (sum_b * sum_b / n_b)) / max(n_b - 1, 1), 0.0)
    se2 = var_a / n_a + var_b / n_b
    delta = mean_a - mean_b
    t_stat = np.divide(delta, np.sqrt(se2), out=np.zeros_like(delta), where=se2 > 0)
    numerator = se2 * se2
    denominator = (var_a * var_a) / (n_a * n_a * max(n_a - 1, 1)) + (var_b * var_b) / (n_b * n_b * max(n_b - 1, 1))
    dof = np.divide(numerator, denominator, out=np.ones_like(numerator), where=denominator > 0)
    p_values = 2.0 * student_t.sf(np.abs(t_stat), dof)
    p_values[~np.isfinite(p_values)] = 1.0
    padj = bh_fdr(p_values)

    pct_a = np.asarray((x_a > 0).sum(axis=0)).ravel().astype(np.float64) / n_a
    pct_b = np.asarray((x_b > 0).sum(axis=0)).ravel().astype(np.float64) / n_b

    out = pd.DataFrame(
        {
            "gene": var_names.astype(str),
            "n_FN": n_a,
            "n_TN": n_b,
            "mean_log1p_cp10k_FN": mean_a,
            "mean_log1p_cp10k_TN": mean_b,
            "delta_log1p_cp10k_FN_minus_TN": delta,
            "pct_detected_FN": pct_a,
            "pct_detected_TN": pct_b,
            "t_stat": t_stat,
            "welch_df": dof,
            "p_value": p_values,
            "padj_bh": padj,
        }
    )
    out["max_pct_detected"] = out[["pct_detected_FN", "pct_detected_TN"]].max(axis=1)
    out["deg_pass"] = (
        (out["padj_bh"] < DE_FDR_CUTOFF)
        & (out["delta_log1p_cp10k_FN_minus_TN"].abs() >= DE_MIN_ABS_DELTA)
        & (out["max_pct_detected"] >= DE_MIN_MAX_PCT_DETECTED)
    )
    out = out.sort_values(["padj_bh", "delta_log1p_cp10k_FN_minus_TN"], ascending=[True, False]).reset_index(drop=True)
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
        if out_row and out_row % 5000 == 0:
            logging.info("Extracted %s / %s selected rows", f"{out_row:,}", f"{selected_idx.size:,}")

    data = np.concatenate(data_parts) if data_parts else np.array([], dtype=np.float32)
    indices = np.concatenate(indices_parts) if indices_parts else np.array([], dtype=np.int32)
    matrix = sparse.csr_matrix((data, indices, indptr), shape=(selected_idx.size, len(var_names)), dtype=np.float32)
    return matrix, var_names


def build_subset_metadata(handle: h5py.File) -> pd.DataFrame:
    required_columns_present(handle)
    source = clean_group_values(read_obs_column(handle, "source_gse_id"))
    library_id = normalize_strings(read_obs_column(handle, "library_id"))
    sample_id = normalize_strings(read_obs_column(handle, "sample_id"))
    tissue = clean_group_values(read_obs_column(handle, "tissue_corrected" if "tissue_corrected" in handle["obs"] else "tissue"))
    has_TRA_TRB_paired = read_bool_obs(handle, "has_TRA_TRB_paired")
    has_TRG_TRD_paired_raw = read_bool_obs(handle, "has_TRG_TRD_paired")
    has_any_ab_tcr = read_bool_obs(handle, "has_any_ab_tcr")
    has_any_gd_tcr_raw = read_bool_obs(handle, "has_any_gd_tcr")
    trg_nonempty = read_nonempty_string_mask(handle, "TRG_cdr3")
    trd_nonempty = read_nonempty_string_mask(handle, "TRD_cdr3")
    trd_score = read_float_obs(handle, "phase4_trd_score")
    trab_score = read_float_obs(handle, "phase4_trab_score")
    trd_minus_trab = read_float_obs(handle, "phase4_trd_minus_trab")

    corrected_has_any_gd_tcr, corrected_has_TRG_TRD_paired, _audit = build_corrected_tcr_evidence(
        source=source,
        has_TRA_TRB_paired=has_TRA_TRB_paired,
        has_TRG_TRD_paired_raw=has_TRG_TRD_paired_raw,
        has_any_ab_tcr=has_any_ab_tcr,
        has_any_gd_tcr_raw=has_any_gd_tcr_raw,
        trg_nonempty=trg_nonempty,
        trd_nonempty=trd_nonempty,
    )
    is_gdtcr_sublib, _sublib = build_sublibrary_labels(
        source=source,
        library_id=library_id,
        sample_id=sample_id,
        corrected_has_any_gd_tcr=corrected_has_any_gd_tcr,
        trg_nonempty=trg_nonempty,
        trd_nonempty=trd_nonempty,
    )
    class_code, _conflicts = make_truth_labels(
        source=source,
        is_gdtcr_sequenced_sublibrary=is_gdtcr_sublib,
        has_TRA_TRB_paired=has_TRA_TRB_paired,
        corrected_has_TRG_TRD_paired=corrected_has_TRG_TRD_paired,
        has_any_ab_tcr=has_any_ab_tcr,
        corrected_has_any_gd_tcr=corrected_has_any_gd_tcr,
        trd_nonempty=trd_nonempty,
    )

    selected_json = json.loads((LOG_DIR / "gdT_prediction_selected_model.json").read_text(encoding="utf-8"))
    model = selected_json["selected_model"]
    threshold = float(selected_json["selected_threshold_on_score"])
    score = {
        "TRD_score_high": trd_score,
        "TRAB_score_low": -trab_score,
        "TRD_minus_TRAB_high": trd_minus_trab,
    }[model]
    predicted_positive = score >= threshold
    tn = (class_code == 1) & (~predicted_positive) & (trd_score > TRD_SCORE_MIN)
    fn = (class_code == 2) & (~predicted_positive) & (trd_score > TRD_SCORE_MIN)
    selected = tn | fn
    group = np.where(fn[selected], "FN", "TN")

    metadata = pd.DataFrame(
        {
            "obs_index": np.flatnonzero(selected).astype(np.int64),
            "group": group,
            "source_gse_id": source[selected],
            "library_id": library_id[selected],
            "sample_id": sample_id[selected],
            "tissue": tissue[selected],
            "phase4_trd_score": trd_score[selected],
            "phase4_trab_score": trab_score[selected],
            "phase4_trd_minus_trab": trd_minus_trab[selected],
            "classifier_model": model,
            "classifier_threshold": threshold,
        }
    )
    return metadata


def compute_within_source_deg(matrix: sparse.csr_matrix, metadata: pd.DataFrame, var_names: pd.Index) -> pd.DataFrame:
    rows = []
    for source, idx in metadata.groupby("source_gse_id", sort=True).groups.items():
        idx_array = np.asarray(list(idx), dtype=np.int64)
        groups = metadata.loc[idx_array, "group"].to_numpy(dtype=object)
        n_fn = int((groups == "FN").sum())
        n_tn = int((groups == "TN").sum())
        if n_fn < MIN_GROUP_CELLS_FOR_SOURCE_DE or n_tn < MIN_GROUP_CELLS_FOR_SOURCE_DE:
            logging.info("Skipping source %s for within-source DE: FN=%s TN=%s", source, n_fn, n_tn)
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
        "delta_log1p_cp10k_FN_minus_TN",
        "padj_bh",
        "p_value",
        "mean_log1p_cp10k_FN",
        "mean_log1p_cp10k_TN",
        "pct_detected_FN",
        "pct_detected_TN",
        "deg_pass",
    ]
    repro = global_deg[base_cols].rename(
        columns={
            "delta_log1p_cp10k_FN_minus_TN": "global_delta_log1p",
            "padj_bh": "global_padj_bh",
            "p_value": "global_p_value",
            "mean_log1p_cp10k_FN": "global_mean_log1p_FN",
            "mean_log1p_cp10k_TN": "global_mean_log1p_TN",
            "pct_detected_FN": "global_pct_detected_FN",
            "pct_detected_TN": "global_pct_detected_TN",
            "deg_pass": "global_deg_pass",
        }
    )
    if source_deg.empty:
        repro["n_sources_tested"] = 0
        return repro

    source_delta = source_deg.pivot(index="gene", columns="source_gse_id", values="delta_log1p_cp10k_FN_minus_TN")
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
        & (repro["n_sources_tested"] >= 2)
        & (repro["source_sign_concordance_fraction"] == 1.0)
        & (repro["n_sources_nominal_p_lt_0p05_same_direction"] >= 2)
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


def plot_volcano(global_deg: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 6.0), constrained_layout=True)
    x = global_deg["delta_log1p_cp10k_FN_minus_TN"].to_numpy(dtype=np.float64)
    y = -np.log10(np.maximum(global_deg["padj_bh"].to_numpy(dtype=np.float64), 1e-300))
    sig = global_deg["deg_pass"].to_numpy(dtype=bool)
    ax.scatter(x[~sig], y[~sig], s=4, c="#b8b8b8", alpha=0.5, linewidths=0, rasterized=True)
    ax.scatter(x[sig], y[sig], s=6, c="#b51f2a", alpha=0.75, linewidths=0, rasterized=True)
    ax.axvline(DE_MIN_ABS_DELTA, color="black", linestyle="--", lw=0.9)
    ax.axvline(-DE_MIN_ABS_DELTA, color="black", linestyle="--", lw=0.9)
    ax.axhline(-np.log10(DE_FDR_CUTOFF), color="black", linestyle=":", lw=0.9)
    top = global_deg.loc[sig].head(15)
    for _, row in top.iterrows():
        ax.text(row["delta_log1p_cp10k_FN_minus_TN"], -np.log10(max(row["padj_bh"], 1e-300)), row["gene"], fontsize=7)
    ax.set_xlabel("Mean log1p(CP10K) difference: FN - TN")
    ax.set_ylabel("-log10(BH FDR)")
    ax.set_title("TRD>0.1 FN vs TN global DE")
    fig.savefig(VOLCANO_PNG, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_reproducibility(repro: pd.DataFrame) -> dict[str, float | int | None]:
    delta_cols = [column for column in repro.columns if column.endswith("_delta_log1p") and column != "global_delta_log1p"]
    stats: dict[str, float | int | None] = {"n_sources_with_both_groups": len(delta_cols)}
    if len(delta_cols) < 2:
        return stats
    x_col, y_col = delta_cols[:2]
    plot_df = repro.dropna(subset=[x_col, y_col]).copy()
    x = plot_df[x_col].to_numpy(dtype=np.float64)
    y = plot_df[y_col].to_numpy(dtype=np.float64)
    stats["n_genes_compared"] = int(plot_df.shape[0])
    stats["pearson_r"] = float(pearsonr(x, y).statistic) if plot_df.shape[0] > 2 else None
    stats["spearman_r"] = float(spearmanr(x, y).statistic) if plot_df.shape[0] > 2 else None

    fig, ax = plt.subplots(figsize=(6.3, 6.0), constrained_layout=True)
    strict = plot_df["reproducible_deg_strict"].to_numpy(dtype=bool)
    ax.scatter(x[~strict], y[~strict], s=4, c="#9ca3af", alpha=0.45, linewidths=0, rasterized=True)
    ax.scatter(x[strict], y[strict], s=8, c="#b51f2a", alpha=0.8, linewidths=0, rasterized=True)
    lim = np.nanquantile(np.abs(np.concatenate([x, y])), 0.995)
    lim = max(float(lim), 0.25)
    ax.plot([-lim, lim], [-lim, lim], color="black", linestyle="--", lw=1)
    ax.axhline(0, color="black", lw=0.7)
    ax.axvline(0, color="black", lw=0.7)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(x_col.replace("_delta_log1p", " delta log1p"))
    ax.set_ylabel(y_col.replace("_delta_log1p", " delta log1p"))
    ax.set_title(
        "Within-source logFC reproducibility\n"
        f"Pearson r={stats['pearson_r']:.3f}, Spearman r={stats['spearman_r']:.3f}"
    )
    fig.savefig(REPRO_SCATTER_PNG, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    return stats


def plot_top_heatmap(repro: pd.DataFrame) -> None:
    source_cols = [column for column in repro.columns if column.endswith("_delta_log1p") and column != "global_delta_log1p"]
    heat_cols = ["global_delta_log1p", *source_cols]
    top = repro.loc[repro["reproducible_deg_strict"]].head(30)
    if top.empty:
        top = repro.loc[repro["global_deg_pass"]].head(30)
    if top.empty:
        return
    mat = top.set_index("gene")[heat_cols]
    fig, ax = plt.subplots(figsize=(max(5.5, 1.2 * len(heat_cols)), max(6, 0.25 * mat.shape[0] + 1.5)), constrained_layout=True)
    vmax = float(np.nanquantile(np.abs(mat.to_numpy(dtype=np.float64)), 0.95))
    vmax = max(vmax, 0.1)
    im = ax.imshow(mat.to_numpy(dtype=np.float64), aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(len(heat_cols)), labels=[col.replace("_delta_log1p", "") for col in heat_cols], rotation=45, ha="right")
    ax.set_yticks(np.arange(mat.shape[0]), labels=mat.index.tolist())
    ax.set_title("Top reproducible/global DEG logFC: FN - TN")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Delta log1p(CP10K)")
    fig.savefig(TOP_HEATMAP_PNG, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def write_summary(
    metadata: pd.DataFrame,
    global_deg: pd.DataFrame,
    source_deg: pd.DataFrame,
    repro: pd.DataFrame,
    repro_stats: dict[str, Any],
) -> None:
    counts = pd.read_csv(COUNTS_CSV)
    strict_n = int(repro["reproducible_deg_strict"].sum()) if "reproducible_deg_strict" in repro else 0
    global_pass_n = int(global_deg["deg_pass"].sum())
    lines = [
        "# TRD>0.1 TN vs FN DEG reproducibility",
        "",
        "Definitions:",
        "",
        "- TN: corrected `abT_gold` cells predicted negative by the current conservative gdT classifier.",
        "- FN: corrected `gdT_gold` cells predicted negative by the same classifier.",
        f"- Additional subset: `phase4_trd_score > {TRD_SCORE_MIN}`.",
        "- Expression: temporary `log1p(CP10K)` from count-space `X`; the H5AD was not modified.",
        "",
        f"- Cells analyzed: `{metadata.shape[0]:,}`",
        f"- TN cells: `{int((metadata['group'] == 'TN').sum()):,}`",
        f"- FN cells: `{int((metadata['group'] == 'FN').sum()):,}`",
        f"- Global DEG pass genes: `{global_pass_n:,}` using FDR < `{DE_FDR_CUTOFF}`, abs delta >= `{DE_MIN_ABS_DELTA}`, max detected fraction >= `{DE_MIN_MAX_PCT_DETECTED}`.",
        f"- Strict reproducible genes across datasets with both TN and FN: `{strict_n:,}`.",
        f"- Within-source reproducibility sources: `{', '.join(sorted(source_deg['source_gse_id'].unique())) if not source_deg.empty else 'none'}`.",
        f"- Pearson r across source logFC: `{repro_stats.get('pearson_r')}`",
        f"- Spearman r across source logFC: `{repro_stats.get('spearman_r')}`",
        "",
        "## Cell counts",
        "",
        dataframe_to_markdown(counts),
        "",
        "## Top strict reproducible DEGs",
        "",
        dataframe_to_markdown(repro.loc[repro["reproducible_deg_strict"]].head(30)),
        "",
        "## Top global DEGs",
        "",
        dataframe_to_markdown(global_deg.head(30)),
        "",
    ]
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    setup_logging()
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    input_h5ad = DEFAULT_INPUT_H5AD.resolve()
    stat_before = input_h5ad.stat()

    with h5py.File(input_h5ad, "r") as handle:
        logging.info("Building TN/FN TRD>0.1 metadata")
        metadata = build_subset_metadata(handle)
        if metadata.empty:
            raise RuntimeError("No TN/FN cells passed the TRD score filter.")
        counts = (
            metadata.groupby(["source_gse_id", "group"], as_index=False)
            .agg(n_cells=("group", "size"))
            .sort_values(["source_gse_id", "group"])
        )
        counts.to_csv(COUNTS_CSV, index=False)
        logging.info("Selected %s cells: %s TN, %s FN", f"{metadata.shape[0]:,}", int((metadata["group"] == "TN").sum()), int((metadata["group"] == "FN").sum()))

        logging.info("Extracting log1p(CP10K) matrix for selected rows")
        matrix, var_names = extract_selected_log_matrix(handle, metadata["obs_index"].to_numpy(dtype=np.int64))

    logging.info("Computing global FN vs TN DE")
    global_deg = welch_deg_from_log_matrix(matrix, metadata["group"].to_numpy(dtype=object), var_names)
    global_deg.to_csv(GLOBAL_DEG_CSV, index=False)

    logging.info("Computing within-source DE for reproducibility")
    source_deg = compute_within_source_deg(matrix, metadata, var_names)
    source_deg.to_csv(SOURCE_DEG_CSV, index=False)

    repro = build_reproducibility_table(global_deg, source_deg)
    repro.to_csv(REPRO_CSV, index=False)
    top_repro = repro.loc[repro["reproducible_deg_strict"]].head(200)
    if top_repro.empty:
        top_repro = repro.loc[repro["global_deg_pass"]].head(200)
    top_repro.to_csv(TOP_REPRO_CSV, index=False)

    logging.info("Rendering DEG figures")
    plot_volcano(global_deg)
    repro_stats = plot_reproducibility(repro)
    plot_top_heatmap(repro)
    write_summary(metadata, global_deg, source_deg, repro, repro_stats)

    stat_after = input_h5ad.stat()
    if (stat_before.st_size != stat_after.st_size) or (stat_before.st_mtime_ns != stat_after.st_mtime_ns):
        raise RuntimeError("Input H5AD changed during read-only DEG analysis.")

    logging.info("Saved global DEG table: %s", GLOBAL_DEG_CSV)
    logging.info("Saved reproducibility table: %s", REPRO_CSV)
    logging.info("Saved summary: %s", SUMMARY_MD)


if __name__ == "__main__":
    main()
