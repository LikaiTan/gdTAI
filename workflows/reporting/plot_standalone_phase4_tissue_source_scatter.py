#!/usr/bin/env python3
"""Plot raw TRAB-vs-TRD scatter colored by tissue and source-like metadata."""

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
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = PROJECT_ROOT / "downloads" / "per_gse_h5ad_with_metadata" / "HRA005041_T_cells_subset.h5ad"
DEFAULT_ALIAS = "HRA005041"
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset" / "figures"
FIGURE_DPI = 300


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--alias", type=str, default=DEFAULT_ALIAS)
    parser.add_argument("--tissue-col", type=str, default="tissue")
    parser.add_argument("--source-col", type=str, default="condition")
    parser.add_argument("--max-categories", type=int, default=24)
    return parser.parse_args()


def compress_categories(values: pd.Series, max_categories: int) -> pd.Series:
    """Keep the most common categories and collapse the rest to Other."""
    values = values.fillna("").astype(str).str.strip()
    values = values.mask(values.eq(""), "unknown")
    counts = values.value_counts(dropna=False)
    if len(counts) <= max_categories:
        return values
    keep = set(counts.head(max_categories - 1).index.astype(str))
    return values.where(values.isin(keep), "Other")


def resolve_plot_df(adata: ad.AnnData, tissue_col: str, source_col: str, max_categories: int) -> tuple[pd.DataFrame, str]:
    """Extract the required plotting frame and a note about the source panel."""
    required = ["phase4_trab_score", "phase4_trd_score"]
    for col in required:
        if col not in adata.obs.columns:
            raise KeyError(f"Missing required score column `{col}` in {adata.filename}.")

    if tissue_col not in adata.obs.columns:
        raise KeyError(f"Missing requested tissue column `{tissue_col}` in {adata.filename}.")

    chosen_source = source_col
    if chosen_source not in adata.obs.columns or adata.obs[chosen_source].astype(str).nunique(dropna=False) <= 1:
        for fallback in ["condition", "Stage", "Donor", "sample_id", "Sample_name"]:
            if fallback in adata.obs.columns and adata.obs[fallback].astype(str).nunique(dropna=False) > 1:
                chosen_source = fallback
                break

    if chosen_source not in adata.obs.columns:
        raise KeyError("Could not resolve a usable source-like metadata column for the source panel.")

    frame = pd.DataFrame(
        {
            "phase4_trab_score": adata.obs["phase4_trab_score"].to_numpy(dtype=np.float32, copy=True),
            "phase4_trd_score": adata.obs["phase4_trd_score"].to_numpy(dtype=np.float32, copy=True),
            "tissue_panel": compress_categories(adata.obs[tissue_col], max_categories),
            "source_panel": compress_categories(adata.obs[chosen_source], max_categories),
        }
    )
    return frame, chosen_source


def categorical_palette(categories: list[str]) -> dict[str, tuple[float, float, float]]:
    """Build a deterministic categorical palette."""
    if len(categories) <= 10:
        colors = sns.color_palette("tab10", n_colors=len(categories))
    elif len(categories) <= 20:
        colors = sns.color_palette("tab20", n_colors=len(categories))
    else:
        colors = sns.color_palette("husl", n_colors=len(categories))
    return dict(zip(categories, colors, strict=True))


def draw_panel(ax: plt.Axes, frame: pd.DataFrame, color_col: str, title: str) -> None:
    """Draw one categorical scatter panel."""
    categories = sorted(pd.Index(frame[color_col].astype(str)).unique().tolist())
    palette = categorical_palette(categories)
    for category in categories:
        subset = frame.loc[frame[color_col].eq(category)]
        ax.scatter(
            subset["phase4_trab_score"],
            subset["phase4_trd_score"],
            s=2,
            c=[palette[category]],
            linewidths=0,
            rasterized=True,
            alpha=0.65,
            label=category,
        )
    ax.set_title(title)
    ax.set_xlabel("Raw TRAB score")
    ax.set_ylabel("Raw TRD score")


def write_tissue_facet_plot(frame: pd.DataFrame, figure_dir: Path) -> Path:
    """Write a facet-by-tissue raw TRAB-vs-TRD scatter plot."""
    tissue_counts = frame["tissue_panel"].value_counts()
    tissues = tissue_counts.index.astype(str).tolist()
    n_panels = len(tissues)
    n_cols = 5 if n_panels > 20 else 4
    n_rows = int(np.ceil(n_panels / n_cols))
    x_min = float(frame["phase4_trab_score"].min())
    x_max = float(frame["phase4_trab_score"].max())
    y_min = float(frame["phase4_trd_score"].min())
    y_max = float(frame["phase4_trd_score"].max())

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.2 * n_cols, 3.6 * n_rows), constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()
    point_color = "#2F6690"

    for ax, tissue in zip(axes, tissues, strict=False):
        subset = frame.loc[frame["tissue_panel"].eq(tissue)]
        ax.scatter(
            subset["phase4_trab_score"],
            subset["phase4_trd_score"],
            s=2,
            c=point_color,
            linewidths=0,
            rasterized=True,
            alpha=0.65,
        )
        ax.set_title(f"{tissue}\n(n={len(subset):,})", fontsize=10)
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_xlabel("Raw TRAB score")
        ax.set_ylabel("Raw TRD score")

    for ax in axes[n_panels:]:
        ax.axis("off")

    out_path = figure_dir / "phase4_trab_vs_trd_facet_by_tissue.png"
    fig.savefig(out_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main() -> None:
    args = parse_args()
    figure_dir = OUTPUT_ROOT / f"{args.alias}_phase4"
    figure_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input_h5ad, backed="r")
    frame, chosen_source = resolve_plot_df(adata, args.tissue_col, args.source_col, args.max_categories)
    adata.file.close()

    sns.set_theme(style="whitegrid", context="talk")
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
    draw_panel(axes[0], frame, "tissue_panel", f"Tissue ({args.tissue_col})")
    draw_panel(axes[1], frame, "source_panel", f"Source ({chosen_source})")

    handles0, labels0 = axes[0].get_legend_handles_labels()
    handles1, labels1 = axes[1].get_legend_handles_labels()
    axes[0].legend(handles0, labels0, loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=True, ncol=2, fontsize=8)
    axes[1].legend(handles1, labels1, loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=True, ncol=1, fontsize=8)

    out_path = figure_dir / "phase4_trab_vs_trd_tissue_source.png"
    fig.savefig(out_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    facet_path = write_tissue_facet_plot(frame, figure_dir)
    print(out_path)
    print(facet_path)


if __name__ == "__main__":
    main()
