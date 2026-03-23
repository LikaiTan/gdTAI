#!/usr/bin/env python3
"""Render Phase 4 threshold count barplots from exported summary tables.

This helper intentionally reads the compact CSV summaries rather than the
large integrated H5AD. It produces publication-readable PNG barplots on NFS
for the four threshold/grouping combinations requested by the user.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
TABLE_DIR = OUTPUT_ROOT / "tables"
FIGURE_DIR = OUTPUT_ROOT / "figures"
LOG_DIR = OUTPUT_ROOT / "logs"

FIGURE_DPI = 300
BASE_WIDTH = 12
BASE_HEIGHT = 4
ROW_HEIGHT = 0.28

PLOT_SPECS = [
    {
        "table": TABLE_DIR / "phase4_trab_minus_trd_lt_neg0p6_by_tissue.csv",
        "group_col": "tissue_corrected",
        "title": "TRAB - TRD < -0.6 Cells by Tissue",
        "xlabel": "Cell Count",
        "ylabel": "Tissue",
        "color": "#1f4e79",
        "output": FIGURE_DIR / "phase4_trab_minus_trd_lt_neg0p6_by_tissue_barplot.png",
    },
    {
        "table": TABLE_DIR / "phase4_trab_minus_trd_lt_neg0p6_by_gse.csv",
        "group_col": "source_gse_id",
        "title": "TRAB - TRD < -0.6 Cells by GSE",
        "xlabel": "Cell Count",
        "ylabel": "GSE",
        "color": "#7a2e2e",
        "output": FIGURE_DIR / "phase4_trab_minus_trd_lt_neg0p6_by_gse_barplot.png",
    },
    {
        "table": TABLE_DIR / "phase4_trd_score_scaled_gt_0p1_by_tissue.csv",
        "group_col": "tissue_corrected",
        "title": "Scaled TRD > 0.1 Cells by Tissue",
        "xlabel": "Cell Count",
        "ylabel": "Tissue",
        "color": "#0f766e",
        "output": FIGURE_DIR / "phase4_trd_score_scaled_gt_0p1_by_tissue_barplot.png",
    },
    {
        "table": TABLE_DIR / "phase4_trd_score_scaled_gt_0p1_by_gse.csv",
        "group_col": "source_gse_id",
        "title": "Scaled TRD > 0.1 Cells by GSE",
        "xlabel": "Cell Count",
        "ylabel": "GSE",
        "color": "#8a5a00",
        "output": FIGURE_DIR / "phase4_trd_score_scaled_gt_0p1_by_gse_barplot.png",
    },
]


def configure_logging() -> None:
    """Configure concise logging for the plotting step."""
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(LOG_DIR / "plot_phase4_threshold_barplots.log", mode="a", encoding="utf-8"),
            logging.StreamHandler(),
        ],
        force=True,
    )


def ensure_output_dirs() -> None:
    """Create required NFS-side output directories."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)


def load_sorted_table(path: Path, group_col: str) -> pd.DataFrame:
    """Load one threshold summary table and sort it from high to low count."""
    frame = pd.read_csv(path)
    required = {group_col, "cell_n"}
    missing = required.difference(frame.columns)
    if missing:
        raise KeyError(f"{path} is missing required columns: {sorted(missing)}")
    frame[group_col] = frame[group_col].fillna("").astype(str)
    frame["cell_n"] = pd.to_numeric(frame["cell_n"], errors="raise")
    return frame.sort_values(["cell_n", group_col], ascending=[False, True]).reset_index(drop=True)


def annotate_bars(ax: plt.Axes, values: pd.Series) -> None:
    """Add count labels to the right edge of each horizontal bar."""
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


def render_barplot(frame: pd.DataFrame, *, group_col: str, title: str, xlabel: str, ylabel: str, color: str, output: Path) -> None:
    """Render one descending horizontal barplot to PNG."""
    n_rows = max(len(frame), 1)
    height = max(BASE_HEIGHT, n_rows * ROW_HEIGHT + 1.5)
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(BASE_WIDTH, height))
    sns.barplot(data=frame, x="cell_n", y=group_col, color=color, ax=ax)
    annotate_bars(ax, frame["cell_n"])
    ax.set_title(title, fontsize=15, pad=12)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.tick_params(axis="y", labelsize=9)
    ax.tick_params(axis="x", labelsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(output, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    logging.info("Wrote %s", output)


def main() -> None:
    """Render the four requested Phase 4 threshold barplots."""
    configure_logging()
    ensure_output_dirs()
    for spec in PLOT_SPECS:
        logging.info("Rendering barplot from %s", spec["table"])
        frame = load_sorted_table(spec["table"], spec["group_col"])
        render_barplot(
            frame,
            group_col=spec["group_col"],
            title=spec["title"],
            xlabel=spec["xlabel"],
            ylabel=spec["ylabel"],
            color=spec["color"],
            output=spec["output"],
        )


if __name__ == "__main__":
    main()
