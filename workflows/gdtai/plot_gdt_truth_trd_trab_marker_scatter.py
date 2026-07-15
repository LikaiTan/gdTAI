#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Plot TRD-vs-TRAB score space for gold/silver truth labels colored by markers."""

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


from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from run_gdt_prediction_package_evaluation import (
    DEFAULT_INPUT_H5AD,
    FIGURE_DIR,
    TABLE_DIR,
    build_corrected_tcr_evidence,
    build_sublibrary_labels,
    categorical_truth,
    clean_group_values,
    make_truth_labels,
    normalize_strings,
    read_bool_obs,
    read_float_obs,
    read_nonempty_string_mask,
    read_obs_column,
    read_string_dataset,
    required_columns_present,
)


MARKER_GENES = ["CD4", "CD8A", "CD8B", "FOXP3"]
OUTPUT_PNG = FIGURE_DIR / "trd_vs_trab_gold_silver_marker_expression.png"
COUNTS_CSV = TABLE_DIR / "trd_vs_trab_gold_silver_marker_scatter_counts.csv"
EXPR_SUMMARY_CSV = TABLE_DIR / "trd_vs_trab_gold_silver_marker_expression_summary.csv"
TARGET_SUM = 10_000.0
FIGURE_DPI = 300


def read_var_names(handle: h5py.File) -> pd.Index:
    return pd.Index(read_string_dataset(handle["var"]["_index"]), dtype="string")


def extract_log1p_normalized_marker_expression(
    handle: h5py.File,
    cell_idx: np.ndarray,
    marker_genes: list[str],
) -> pd.DataFrame:
    """Extract per-cell log1p(CP10K) marker expression without mutating the H5AD."""
    var_names = read_var_names(handle)
    gene_idx = var_names.get_indexer(pd.Index(marker_genes, dtype="string"))
    if np.any(gene_idx < 0):
        missing = [gene for gene, idx in zip(marker_genes, gene_idx.tolist(), strict=False) if idx < 0]
        raise ValueError(f"Missing marker genes in {DEFAULT_INPUT_H5AD}: {missing}")

    x_group = handle["X"]
    indptr = x_group["indptr"]
    indices = x_group["indices"]
    data = x_group["data"]
    marker_lookup = {int(var_idx): pos for pos, var_idx in enumerate(gene_idx.tolist())}
    out = np.zeros((cell_idx.size, len(marker_genes)), dtype=np.float32)

    for out_row, obs_idx in enumerate(cell_idx):
        start = int(indptr[obs_idx])
        end = int(indptr[obs_idx + 1])
        if end <= start:
            continue
        row_indices = indices[start:end]
        row_data = data[start:end]
        row_sum = float(np.sum(row_data, dtype=np.float64))
        if row_sum <= 0:
            continue
        positions = np.searchsorted(row_indices, gene_idx)
        valid = positions < row_indices.shape[0]
        if np.any(valid):
            valid_positions = positions[valid]
            valid_gene_idx = gene_idx[valid]
            valid[valid] = row_indices[valid_positions] == valid_gene_idx
        if not np.any(valid):
            continue
        scale = TARGET_SUM / row_sum
        for marker_var_idx, position in zip(gene_idx[valid], positions[valid], strict=False):
            out[out_row, marker_lookup[int(marker_var_idx)]] = np.log1p(float(row_data[position]) * scale)

    return pd.DataFrame(out, columns=marker_genes)


def build_truth_labels(handle: h5py.File) -> tuple[np.ndarray, pd.Categorical]:
    required_columns_present(handle)
    source = clean_group_values(read_obs_column(handle, "source_gse_id"))
    library_id = normalize_strings(read_obs_column(handle, "library_id"))
    sample_id = normalize_strings(read_obs_column(handle, "sample_id"))
    has_TRA_TRB_paired = read_bool_obs(handle, "has_TRA_TRB_paired")
    has_TRG_TRD_paired_raw = read_bool_obs(handle, "has_TRG_TRD_paired")
    has_any_ab_tcr = read_bool_obs(handle, "has_any_ab_tcr")
    has_any_gd_tcr_raw = read_bool_obs(handle, "has_any_gd_tcr")
    trg_nonempty = read_nonempty_string_mask(handle, "TRG_cdr3")
    trd_nonempty = read_nonempty_string_mask(handle, "TRD_cdr3")

    corrected_has_any_gd_tcr, corrected_has_TRG_TRD_paired, _audit_df = build_corrected_tcr_evidence(
        source=source,
        has_TRA_TRB_paired=has_TRA_TRB_paired,
        has_TRG_TRD_paired_raw=has_TRG_TRD_paired_raw,
        has_any_ab_tcr=has_any_ab_tcr,
        has_any_gd_tcr_raw=has_any_gd_tcr_raw,
        trg_nonempty=trg_nonempty,
        trd_nonempty=trd_nonempty,
    )
    is_gdtcr_sequenced_sublibrary, _sublibrary_summary = build_sublibrary_labels(
        source=source,
        library_id=library_id,
        sample_id=sample_id,
        corrected_has_any_gd_tcr=corrected_has_any_gd_tcr,
        trg_nonempty=trg_nonempty,
        trd_nonempty=trd_nonempty,
    )
    class_code, _conflict_df = make_truth_labels(
        source=source,
        is_gdtcr_sequenced_sublibrary=is_gdtcr_sequenced_sublibrary,
        has_TRA_TRB_paired=has_TRA_TRB_paired,
        corrected_has_TRG_TRD_paired=corrected_has_TRG_TRD_paired,
        has_any_ab_tcr=has_any_ab_tcr,
        corrected_has_any_gd_tcr=corrected_has_any_gd_tcr,
        trd_nonempty=trd_nonempty,
    )
    return class_code, categorical_truth(class_code)


def write_summary_tables(plot_df: pd.DataFrame) -> None:
    counts = (
        plot_df.groupby(["truth_class", "source_gse_id"], observed=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["truth_class", "n_cells", "source_gse_id"], ascending=[True, False, True])
    )
    counts.to_csv(COUNTS_CSV, index=False)

    summary_rows = []
    for truth_class, group in plot_df.groupby("truth_class", observed=False):
        for gene in MARKER_GENES:
            values = group[gene].to_numpy(dtype=np.float32, copy=False)
            summary_rows.append(
                {
                    "truth_class": truth_class,
                    "gene": gene,
                    "n_cells": int(values.size),
                    "mean_log1p_cp10k": float(np.mean(values)),
                    "median_log1p_cp10k": float(np.median(values)),
                    "pct_detected": float(np.mean(values > 0.0)),
                    "q90_log1p_cp10k": float(np.quantile(values, 0.90)),
                    "q99_log1p_cp10k": float(np.quantile(values, 0.99)),
                }
            )
    pd.DataFrame(summary_rows).to_csv(EXPR_SUMMARY_CSV, index=False)


def plot_marker_scatter(plot_df: pd.DataFrame) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    x = plot_df["phase4_trab_score"].to_numpy(dtype=np.float32, copy=False)
    y = plot_df["phase4_trd_score"].to_numpy(dtype=np.float32, copy=False)
    xlim = np.quantile(x[np.isfinite(x)], [0.001, 0.999])
    ylim = np.quantile(y[np.isfinite(y)], [0.001, 0.999])
    line_x = np.linspace(float(xlim[0]), float(xlim[1]), 200)
    cutoff = 0.37491393089294434

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 11.5), constrained_layout=True)
    axes = axes.ravel()
    truth_counts = plot_df["truth_class"].value_counts().to_dict()
    subtitle = ", ".join(f"{label}: {count:,}" for label, count in truth_counts.items())

    for ax, gene in zip(axes, MARKER_GENES, strict=False):
        values = plot_df[gene].to_numpy(dtype=np.float32, copy=False)
        order = np.argsort(values, kind="mergesort")
        vmax = float(np.quantile(values[np.isfinite(values)], 0.995))
        scatter = ax.scatter(
            x[order],
            y[order],
            c=values[order],
            s=1.0,
            cmap="viridis",
            vmin=0.0,
            vmax=max(vmax, 1e-6),
            linewidths=0,
            alpha=0.75,
            rasterized=True,
        )
        ax.plot(line_x, line_x + cutoff, color="black", linestyle="--", linewidth=1.0, label="TRD - TRAB = 0.3749")
        ax.set_xlim(float(xlim[0]), float(xlim[1]))
        ax.set_ylim(float(ylim[0]), float(ylim[1]))
        ax.set_title(f"{gene} expression")
        ax.set_xlabel("phase4_trab_score")
        ax.set_ylabel("phase4_trd_score")
        ax.grid(alpha=0.15)
        ax.legend(frameon=False, fontsize=8, loc="upper left")
        cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("log1p(CP10K)")

    fig.suptitle(
        "TRD vs TRAB score space for corrected gdT/abT truth labels\n"
        f"Cells shown: {plot_df.shape[0]:,}; {subtitle}. Axes clipped to 0.1-99.9% quantiles.",
        fontsize=14,
    )
    fig.savefig(OUTPUT_PNG, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    input_h5ad = DEFAULT_INPUT_H5AD.resolve()
    stat_before = input_h5ad.stat()
    with h5py.File(input_h5ad, "r") as handle:
        class_code, truth = build_truth_labels(handle)
        selected_idx = np.flatnonzero(class_code > 0).astype(np.int64)
        source = clean_group_values(read_obs_column(handle, "source_gse_id"))[selected_idx]
        trd_score = read_float_obs(handle, "phase4_trd_score")[selected_idx]
        trab_score = read_float_obs(handle, "phase4_trab_score")[selected_idx]
        expr_df = extract_log1p_normalized_marker_expression(handle, selected_idx, MARKER_GENES)

    plot_df = pd.DataFrame(
        {
            "truth_class": pd.Series(np.asarray(truth[selected_idx], dtype=object)).astype(str),
            "source_gse_id": source,
            "phase4_trd_score": trd_score,
            "phase4_trab_score": trab_score,
        }
    )
    plot_df = pd.concat([plot_df, expr_df], axis=1)
    plot_marker_scatter(plot_df)
    write_summary_tables(plot_df)

    stat_after = input_h5ad.stat()
    if (stat_before.st_size != stat_after.st_size) or (stat_before.st_mtime_ns != stat_after.st_mtime_ns):
        raise RuntimeError("Input H5AD changed during read-only plotting.")

    print(f"saved_png\t{OUTPUT_PNG}")
    print(f"saved_csv\t{COUNTS_CSV}")
    print(f"saved_csv\t{EXPR_SUMMARY_CSV}")
    print(f"n_cells_plotted\t{plot_df.shape[0]}")


if __name__ == "__main__":
    main()
