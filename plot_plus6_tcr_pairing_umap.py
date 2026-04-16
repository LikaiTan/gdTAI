#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Plot plus6 UMAP panels for paired abTCR, paired gdTCR, and paired-both cells."""

from __future__ import annotations

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent
INPUT_H5AD = PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "integrated_plus6.h5ad"
FIGURE_DIR = PROJECT_ROOT / "Integrated_dataset" / "figures" / "plus6"
TABLE_DIR = PROJECT_ROOT / "Integrated_dataset" / "tables" / "plus6"
OUTPUT_PNG = FIGURE_DIR / "plus6_umap_paired_tcr_doublets.png"
OUTPUT_CSV = TABLE_DIR / "plus6_umap_paired_tcr_doublets_counts.csv"

SAMPLE_N = 300_000
RANDOM_SEED = 0
FIGURE_DPI = 300


def choose_sample_indices(n_obs: int, sample_n: int, seed: int) -> np.ndarray:
    if n_obs <= sample_n:
        return np.arange(n_obs, dtype=np.int64)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(np.arange(n_obs), size=sample_n, replace=False))


def load_plot_frame(h5ad_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    with h5py.File(h5ad_path, "r") as handle:
        n_obs = int(handle["obs"]["_index"].shape[0])
        umap_all = np.asarray(handle["obsm"]["X_umap"][:], dtype=np.float32)
        paired_ab_all = np.asarray(handle["obs"]["has_TRA_TRB_paired"][:], dtype=bool)
        paired_gd_all = np.asarray(handle["obs"]["has_TRG_TRD_paired"][:], dtype=bool)
        idx = choose_sample_indices(n_obs, SAMPLE_N, RANDOM_SEED)
        umap = umap_all[idx]
        paired_ab = paired_ab_all[idx]
        paired_gd = paired_gd_all[idx]

    doublet = paired_ab & paired_gd
    plot_df = pd.DataFrame(
        {
            "umap1": umap[:, 0],
            "umap2": umap[:, 1],
            "paired_ab": paired_ab,
            "paired_gd": paired_gd,
            "doublet": doublet,
        }
    )
    count_df = pd.DataFrame(
        {
            "category": ["paired_TRA_TRB", "paired_TRG_TRD", "paired_both_doublet_proxy"],
            "sampled_cells": [int(paired_ab.sum()), int(paired_gd.sum()), int(doublet.sum())],
            "sampled_total": [int(plot_df.shape[0])] * 3,
        }
    )
    return plot_df, count_df


def draw_binary_panel(ax: plt.Axes, plot_df: pd.DataFrame, column: str, title: str) -> None:
    false_df = plot_df.loc[~plot_df[column]]
    true_df = plot_df.loc[plot_df[column]]
    ax.scatter(
        false_df["umap1"],
        false_df["umap2"],
        s=2,
        c="#4C78A8",
        linewidths=0,
        rasterized=True,
        alpha=0.35,
        label="others",
    )
    ax.scatter(
        true_df["umap1"],
        true_df["umap2"],
        s=2,
        c="#D62728",
        linewidths=0,
        rasterized=True,
        alpha=0.75,
        label=column,
    )
    ax.set_title(title)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(loc="upper right", frameon=True, fontsize=8)


def main() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    plot_df, count_df = load_plot_frame(INPUT_H5AD)
    count_df.to_csv(OUTPUT_CSV, index=False)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
    draw_binary_panel(axes[0], plot_df, "paired_ab", "paired TRA/TRB")
    draw_binary_panel(axes[1], plot_df, "paired_gd", "paired TRG/TRD")
    draw_binary_panel(axes[2], plot_df, "doublet", "paired both (doublet proxy)")
    fig.savefig(OUTPUT_PNG, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"saved_png\t{OUTPUT_PNG}")
    print(f"saved_csv\t{OUTPUT_CSV}")


if __name__ == "__main__":
    main()
