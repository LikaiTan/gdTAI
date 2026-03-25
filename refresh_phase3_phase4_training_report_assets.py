#!/usr/bin/env python3
"""Refresh report-specific Phase 3/4 assets for the training-report HTML.

This helper generates only the report assets that were missing or mismatched:
- Phase 3 UMAP by tissue
- Raw `TRD score > 0.1` threshold summaries by tissue and by GSE
- Raw-threshold barplots derived from those summaries
"""

from __future__ import annotations

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
FIGURE_DIR = OUTPUT_ROOT / "figures"
TABLE_DIR = OUTPUT_ROOT / "tables"
INTEGRATED_H5AD = PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "integrated.h5ad"

RANDOM_SEED = 1
UMAP_SAMPLE_SIZE = 180_000
FIGURE_DPI = 300
BASE_WIDTH = 12
BASE_HEIGHT = 4
ROW_HEIGHT = 0.28


def read_strings(dataset: h5py.Dataset) -> np.ndarray:
    """Read a string-like HDF5 dataset into a plain object array."""
    values = dataset[:]
    return np.asarray(
        [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values],
        dtype=object,
    )


def ensure_output_dirs() -> None:
    """Create NFS-side output directories used by the report assets."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)


def sample_phase3_tissue_frame() -> pd.DataFrame:
    """Sample UMAP coordinates and tissue labels for the Phase 3 tissue plot."""
    with h5py.File(INTEGRATED_H5AD, "r") as handle:
        umap = np.asarray(handle["obsm"]["X_umap"][:], dtype=np.float32)
        tissues = read_strings(handle["obs"]["tissue_corrected"])
    n_obs = umap.shape[0]
    rng = np.random.default_rng(RANDOM_SEED)
    sample_n = min(UMAP_SAMPLE_SIZE, n_obs)
    sample_idx = np.sort(rng.choice(n_obs, size=sample_n, replace=False))
    return pd.DataFrame(
        {
            "umap_1": umap[sample_idx, 0],
            "umap_2": umap[sample_idx, 1],
            "tissue_corrected": tissues[sample_idx],
        }
    )


def write_phase3_umap_by_tissue() -> None:
    """Render a Phase 3 UMAP colored by corrected tissue."""
    frame = sample_phase3_tissue_frame()
    tissues = sorted(frame["tissue_corrected"].astype(str).unique().tolist())
    palette = dict(zip(tissues, sns.color_palette("tab20", n_colors=len(tissues)), strict=False))
    fig, ax = plt.subplots(figsize=(9.2, 7.4))
    sns.scatterplot(
        data=frame,
        x="umap_1",
        y="umap_2",
        hue="tissue_corrected",
        palette=palette,
        s=4,
        linewidth=0,
        alpha=0.72,
        rasterized=True,
        ax=ax,
    )
    ax.set_title("Phase 3 UMAP by Tissue", fontsize=16)
    ax.set_xlabel("UMAP1", fontsize=11)
    ax.set_ylabel("UMAP2", fontsize=11)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    legend = ax.get_legend()
    if legend is not None:
        legend.set_title("Tissue")
        for text in legend.get_texts():
            text.set_fontsize(8)
        legend.get_title().set_fontsize(9)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "phase3_umap_by_tissue.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def read_raw_trd_threshold_frames() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarize cells with raw `TRD score > 0.1` by tissue and by GSE."""
    with h5py.File(INTEGRATED_H5AD, "r") as handle:
        obs = handle["obs"]
        trd = np.asarray(obs["phase4_trd_score"][:], dtype=np.float32)
        tissues = read_strings(obs["tissue_corrected"])
        gses = read_strings(obs["source_gse_id"])
    mask = trd > 0.1
    tissue_frame = (
        pd.Series(tissues[mask], name="tissue_corrected")
        .value_counts(dropna=False)
        .rename_axis("tissue_corrected")
        .reset_index(name="cell_n")
        .sort_values(["cell_n", "tissue_corrected"], ascending=[False, True])
        .reset_index(drop=True)
    )
    gse_frame = (
        pd.Series(gses[mask], name="source_gse_id")
        .value_counts(dropna=False)
        .rename_axis("source_gse_id")
        .reset_index(name="cell_n")
        .sort_values(["cell_n", "source_gse_id"], ascending=[False, True])
        .reset_index(drop=True)
    )
    return tissue_frame, gse_frame


def write_raw_trd_threshold_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Write raw `TRD score > 0.1` summary tables to NFS."""
    tissue_frame, gse_frame = read_raw_trd_threshold_frames()
    tissue_frame.to_csv(TABLE_DIR / "phase4_trd_gt_0p1_by_tissue.csv", index=False)
    gse_frame.to_csv(TABLE_DIR / "phase4_trd_gt_0p1_by_gse.csv", index=False)
    return tissue_frame, gse_frame


def annotate_bars(ax: plt.Axes, values: pd.Series) -> None:
    """Write value labels at the right edge of each horizontal bar."""
    max_value = float(values.max()) if len(values) else 0.0
    offset = max(max_value * 0.01, 1.0)
    for patch, value in zip(ax.patches, values, strict=False):
        ax.text(
            patch.get_width() + offset,
            patch.get_y() + patch.get_height() / 2.0,
            f"{int(value):,}",
            va="center",
            ha="left",
            fontsize=8,
        )


def write_horizontal_barplot(
    frame: pd.DataFrame,
    *,
    group_col: str,
    title: str,
    ylabel: str,
    color: str,
    output_name: str,
) -> None:
    """Render one sorted horizontal barplot."""
    n_rows = max(len(frame), 1)
    height = max(BASE_HEIGHT, n_rows * ROW_HEIGHT + 1.5)
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(BASE_WIDTH, height))
    sns.barplot(data=frame, x="cell_n", y=group_col, color=color, ax=ax)
    annotate_bars(ax, frame["cell_n"])
    ax.set_title(title, fontsize=15, pad=12)
    ax.set_xlabel("Cell Count", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.tick_params(axis="y", labelsize=9)
    ax.tick_params(axis="x", labelsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / output_name, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    """Refresh the missing report assets on NFS."""
    ensure_output_dirs()
    write_phase3_umap_by_tissue()
    tissue_frame, gse_frame = write_raw_trd_threshold_tables()
    write_horizontal_barplot(
        tissue_frame,
        group_col="tissue_corrected",
        title="Raw TRD > 0.1 Cells by Tissue",
        ylabel="Tissue",
        color="#0f766e",
        output_name="phase4_trd_gt_0p1_by_tissue_barplot.png",
    )
    write_horizontal_barplot(
        gse_frame,
        group_col="source_gse_id",
        title="Raw TRD > 0.1 Cells by GSE",
        ylabel="GSE",
        color="#8a5a00",
        output_name="phase4_trd_gt_0p1_by_gse_barplot.png",
    )


if __name__ == "__main__":
    main()
