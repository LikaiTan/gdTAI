#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Plot a 15-panel T/NK marker UMAP figure from the current integrated object.

This helper reads the Phase 3 integrated milestone from the SSD-backed mirrored
tree, samples a fixed number of cells for plotting readability, extracts a
small marker panel directly from the H5AD CSR matrix, and saves the resulting
PNG figure plus a small gene-order table back to NFS.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Config
INTEGRATED_H5AD = Path("/ssd/tnk_phase3/Integrated_dataset/integrated.h5ad")
FIGURE_DIR = Path("Integrated_dataset/figures")
TABLE_DIR = Path("Integrated_dataset/tables")
OUTPUT_PNG = FIGURE_DIR / "phase3_tnk_marker_umap_panel.png"
OUTPUT_GENE_CSV = TABLE_DIR / "phase3_tnk_marker_umap_panel_genes.csv"

# Balanced panel across pan-T, alpha-beta T, gamma-delta T, and NK programs.
MARKER_GENES = [
    "CD3D",
    "CD3E",
    "TRAC",
    "TRBC1",
    "TRBC2",
    "IL7R",
    "LTB",
    "TRDC",
    "TRGC1",
    "TRGC2",
    "NKG7",
    "KLRD1",
    "KLRC1",
    "CTSW",
    "CCL5",
]

QC_SAMPLE_MAX_CELLS = 30_000
RANDOM_SEED = 0
FIGSIZE = (18, 22)
POINT_SIZE = 1.2
CMAP = "viridis"


def choose_sample_indices(n_obs: int, max_cells: int, seed: int) -> np.ndarray:
    """Choose a reproducible cell sample for readable multi-panel UMAPs."""
    if n_obs <= max_cells:
        return np.arange(n_obs, dtype=np.int64)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(n_obs, size=max_cells, replace=False).astype(np.int64))


def decode_strings(raw: np.ndarray) -> list[str]:
    """Decode HDF5 string arrays into plain Python strings."""
    values: list[str] = []
    for item in raw:
        if isinstance(item, bytes):
            values.append(item.decode("utf-8"))
        else:
            values.append(str(item))
    return values


def extract_marker_expression(
    h5ad_path: Path,
    sample_idx: np.ndarray,
    marker_genes: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    """Extract UMAP coordinates and marker expression for the sampled cells.

    The H5AD stores X as CSR. For a small marker panel, row-wise sparse lookup
    over the sampled cells is more memory-efficient than loading the full AnnData.
    """
    with h5py.File(h5ad_path, "r") as handle:
        var_names = decode_strings(handle["var"]["_index"][...])
        gene_to_idx = {gene: idx for idx, gene in enumerate(var_names)}
        missing = [gene for gene in marker_genes if gene not in gene_to_idx]
        if missing:
            raise ValueError(f"Missing marker genes in {h5ad_path}: {missing}")

        marker_var_idx = np.array([gene_to_idx[gene] for gene in marker_genes], dtype=np.int64)
        marker_lookup = {var_idx: pos for pos, var_idx in enumerate(marker_var_idx)}

        umap = np.asarray(handle["obsm"]["X_umap"][sample_idx], dtype=np.float32)
        indptr = handle["X"]["indptr"]
        indices = handle["X"]["indices"]
        data = handle["X"]["data"]

        expression = np.zeros((sample_idx.shape[0], len(marker_genes)), dtype=np.float32)
        for out_row, cell_idx in enumerate(sample_idx):
            start = int(indptr[cell_idx])
            end = int(indptr[cell_idx + 1])
            row_indices = indices[start:end]
            row_data = data[start:end]
            keep_mask = np.isin(row_indices, marker_var_idx, assume_unique=False)
            if not np.any(keep_mask):
                continue
            for gene_idx, value in zip(row_indices[keep_mask], row_data[keep_mask], strict=False):
                expression[out_row, marker_lookup[int(gene_idx)]] = float(value)

    return umap, expression


def plot_marker_panel(umap: np.ndarray, expression: np.ndarray, marker_genes: list[str]) -> None:
    """Render one 15-panel UMAP PNG with log1p marker expression overlays."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    plot_values = np.log1p(np.maximum(expression, 0.0))
    n_rows, n_cols = 5, 3
    fig, axes = plt.subplots(n_rows, n_cols, figsize=FIGSIZE, constrained_layout=True)
    axes = axes.ravel()

    for ax, gene, values in zip(axes, marker_genes, plot_values.T, strict=False):
        order = np.argsort(values)
        sorted_xy = umap[order]
        sorted_values = values[order]
        vmax = float(np.percentile(sorted_values, 99.0)) if np.any(sorted_values > 0) else 1.0
        if vmax <= 0:
            vmax = 1.0

        scatter = ax.scatter(
            sorted_xy[:, 0],
            sorted_xy[:, 1],
            c=sorted_values,
            s=POINT_SIZE,
            cmap=CMAP,
            vmin=0.0,
            vmax=vmax,
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(f"{gene} (log1p counts)", fontsize=11)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.02)
        cbar.ax.tick_params(labelsize=8)

    fig.suptitle("Phase 3 T/NK Marker UMAP Panel", fontsize=18, y=1.01)
    fig.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(
        {
            "panel_order": np.arange(1, len(marker_genes) + 1, dtype=int),
            "gene": marker_genes,
            "source_h5ad": str(INTEGRATED_H5AD),
            "sampled_cells": int(umap.shape[0]),
        }
    ).to_csv(OUTPUT_GENE_CSV, index=False)


def main() -> None:
    """Build the sampled marker panel and write outputs to NFS."""
    if not INTEGRATED_H5AD.exists():
        raise FileNotFoundError(f"Integrated H5AD not found: {INTEGRATED_H5AD}")

    with h5py.File(INTEGRATED_H5AD, "r") as handle:
        n_obs = int(handle["obs"]["_index"].shape[0])
    sample_idx = choose_sample_indices(n_obs=n_obs, max_cells=QC_SAMPLE_MAX_CELLS, seed=RANDOM_SEED)
    umap, expression = extract_marker_expression(
        h5ad_path=INTEGRATED_H5AD,
        sample_idx=sample_idx,
        marker_genes=MARKER_GENES,
    )
    plot_marker_panel(umap=umap, expression=expression, marker_genes=MARKER_GENES)
    print(f"saved_png\t{OUTPUT_PNG}")
    print(f"saved_gene_table\t{OUTPUT_GENE_CSV}")
    print(f"sampled_cells\t{sample_idx.shape[0]}")


if __name__ == "__main__":
    main()
