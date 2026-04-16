#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Plot a plus6 UMAP highlighting Sorted_gdT cells."""

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
OUTPUT_PNG = FIGURE_DIR / "plus6_umap_sorted_gdt_highlight.png"
OUTPUT_CSV = TABLE_DIR / "plus6_umap_sorted_gdt_counts.csv"

SAMPLE_N = 300_000
RANDOM_SEED = 0
FIGURE_DPI = 300


def choose_sample_indices(n_obs: int, sample_n: int, seed: int) -> np.ndarray:
    if n_obs <= sample_n:
        return np.arange(n_obs, dtype=np.int64)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(np.arange(n_obs), size=sample_n, replace=False))


def read_boolean_obs(handle: h5py.File, column: str) -> np.ndarray:
    obj = handle["obs"][column]
    if isinstance(obj, h5py.Group):
        categories = np.asarray(
            [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in obj["categories"][:]],
            dtype=object,
        )
        codes = obj["codes"][:]
        out = np.full(codes.shape, "", dtype=object)
        valid = codes >= 0
        out[valid] = categories[codes[valid]]
        lowered = np.char.lower(out.astype(str))
        return np.asarray(lowered == "true", dtype=bool)
    raw = obj[:]
    if raw.dtype == bool:
        return np.asarray(raw, dtype=bool)
    lowered = np.char.lower(raw.astype(str))
    return np.asarray(lowered == "true", dtype=bool)


def main() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    with h5py.File(INPUT_H5AD, "r") as handle:
        n_obs = int(handle["obs"]["_index"].shape[0])
        umap_all = np.asarray(handle["obsm"]["X_umap"][:], dtype=np.float32)
        sorted_gdt_all = read_boolean_obs(handle, "Sorted_gdT")
        idx = choose_sample_indices(n_obs, SAMPLE_N, RANDOM_SEED)
        umap = umap_all[idx]
        sorted_gdt = sorted_gdt_all[idx]

    plot_df = pd.DataFrame(
        {
            "umap1": umap[:, 0],
            "umap2": umap[:, 1],
            "sorted_gdt": sorted_gdt,
        }
    )
    pd.DataFrame(
        {
            "category": ["Sorted_gdT_true"],
            "sampled_cells": [int(sorted_gdt.sum())],
            "sampled_total": [int(plot_df.shape[0])],
        }
    ).to_csv(OUTPUT_CSV, index=False)

    false_df = plot_df.loc[~plot_df["sorted_gdt"]]
    true_df = plot_df.loc[plot_df["sorted_gdt"]]
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    ax.scatter(
        false_df["umap1"],
        false_df["umap2"],
        s=2,
        c="#4C78A8",
        linewidths=0,
        rasterized=True,
        alpha=0.30,
        label="others",
    )
    ax.scatter(
        true_df["umap1"],
        true_df["umap2"],
        s=2,
        c="#D62728",
        linewidths=0,
        rasterized=True,
        alpha=0.80,
        label="Sorted_gdT = True",
    )
    ax.set_title("plus6 UMAP highlighting Sorted_gdT cells")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(loc="upper right", frameon=True, fontsize=8)
    fig.savefig(OUTPUT_PNG, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"saved_png\t{OUTPUT_PNG}")
    print(f"saved_csv\t{OUTPUT_CSV}")


if __name__ == "__main__":
    main()
