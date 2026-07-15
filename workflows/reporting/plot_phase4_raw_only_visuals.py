#!/usr/bin/env python3
"""Generate raw-only Phase 4 visualizations for the training report.

These plots intentionally avoid scaled-score panels. They use only the raw
Phase 4 score columns stored in the integrated milestone.
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


from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from phase4_gdt_module_scoring import (
    CHUNK_SIZE,
    attach_tcr_presence_flags,
    build_tcr_presence_lookup,
    extract_log1p_gene_expression_for_sample,
    load_obs_join_fields_for_sample,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
FIGURE_DIR = OUTPUT_ROOT / "figures"
INTEGRATED_H5AD = Path("/ssd/tnk_phase3/Integrated_dataset/integrated.h5ad")
RANDOM_SEED = 1
SAMPLE_SIZE = 180_000
GENE_PANEL_SAMPLE_SIZE = 60_000
FIGURE_DPI = 300
RAW_GENE_PANEL = [
    "TRDC",
    "TRDV1",
    "TRDV2",
    "TRGC1",
    "TRGC2",
    "TRGV9",
    "CD4",
    "CD8A",
    "CD8B",
    "CD3D",
    "NKG7",
    "TRBC1",
]


def sample_obs_frame() -> tuple[pd.DataFrame, np.ndarray]:
    """Read the raw score columns and UMAP coordinates for a sampled cell set."""
    with h5py.File(INTEGRATED_H5AD, "r") as handle:
        obs = handle["obs"]
        umap = np.asarray(handle["obsm"]["X_umap"][:], dtype=np.float32)
        trab = np.asarray(obs["phase4_trab_score"][:], dtype=np.float32)
        trd = np.asarray(obs["phase4_trd_score"][:], dtype=np.float32)
        trd_minus_trab = np.asarray(obs["phase4_trd_minus_trab"][:], dtype=np.float32)
        n_obs = umap.shape[0]
        rng = np.random.default_rng(RANDOM_SEED)
        sample_n = min(SAMPLE_SIZE, n_obs)
        sample_idx = np.sort(rng.choice(n_obs, size=sample_n, replace=False))
        frame = pd.DataFrame(
            {
                "umap_1": umap[sample_idx, 0],
                "umap_2": umap[sample_idx, 1],
                "phase4_trab_score": trab[sample_idx],
                "phase4_trd_score": trd[sample_idx],
                "phase4_trd_minus_trab": trd_minus_trab[sample_idx],
            }
        )
    frame["phase4_trab_minus_trd"] = -frame["phase4_trd_minus_trab"]
    return frame, sample_idx


def attach_sample_tcr_flags(frame: pd.DataFrame, sample_idx: np.ndarray) -> pd.DataFrame:
    """Attach paired TRA/TRB versus no-TCR flags for the sampled cells."""
    join_df = load_obs_join_fields_for_sample(
        INTEGRATED_H5AD,
        sample_idx,
        ["project name", "sampleid", "barcodes"],
    )
    lookup_df = build_tcr_presence_lookup()
    merged = attach_tcr_presence_flags(join_df, lookup_df)
    return pd.concat([frame.reset_index(drop=True), merged.reset_index(drop=True)], axis=1)


def add_shared_scatter_style(ax: plt.Axes) -> None:
    """Apply a consistent minimalist style to one scatter axis."""
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#b7b2aa")
    ax.spines["bottom"].set_color("#b7b2aa")


def write_raw_umap_overlays(frame: pd.DataFrame) -> None:
    """Write raw-only UMAP overlays for TRAB, TRD, and TRD-TRAB."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    specs = [
        ("phase4_trab_score", "Raw TRAB score", "viridis"),
        ("phase4_trd_score", "Raw TRD score", "magma"),
        ("phase4_trd_minus_trab", "Raw TRD - TRAB", "coolwarm"),
    ]
    for ax, (column, title, cmap) in zip(axes, specs, strict=False):
        points = ax.scatter(
            frame["umap_1"],
            frame["umap_2"],
            c=frame[column],
            s=2,
            cmap=cmap,
            linewidths=0,
            alpha=0.75,
            rasterized=True,
        )
        add_shared_scatter_style(ax)
        ax.set_title(title, fontsize=12)
        plt.colorbar(points, ax=ax, fraction=0.046, pad=0.03)
    fig.suptitle("Phase 4 Raw Score UMAP Overlays", fontsize=15)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "phase4_umap_raw_score_overlays.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def write_raw_trab_trd_scatter(frame: pd.DataFrame) -> None:
    """Write a raw-only TRAB-versus-TRD scatter colored by raw TRD-minus-TRAB."""
    fig, ax = plt.subplots(figsize=(6.8, 5.8))
    points = ax.scatter(
        frame["phase4_trab_score"],
        frame["phase4_trd_score"],
        c=frame["phase4_trd_minus_trab"],
        cmap="coolwarm",
        s=4,
        linewidths=0,
        alpha=0.65,
        rasterized=True,
    )
    ax.set_title("Raw TRAB-versus-TRD Score Space", fontsize=15)
    ax.set_xlabel("Raw TRAB score", fontsize=11)
    ax.set_ylabel("Raw TRD score", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.colorbar(points, ax=ax, fraction=0.046, pad=0.03, label="Raw TRD - TRAB")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "phase4_trab_vs_trd_raw_only.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def write_raw_paired_tcr_scatter(frame: pd.DataFrame) -> None:
    """Write a raw-only TRAB-versus-TRD scatter for paired TRA/TRB versus no TCR."""
    paired_df = frame.loc[
        frame["has_paired_tra_trb_cdr3"] | frame["has_no_tra_trb_cdr3"]
    ].copy()
    paired_df["tcr_pairing_group"] = np.where(
        paired_df["has_paired_tra_trb_cdr3"],
        "Paired TRA/TRB",
        "No TCR",
    )
    palette = {"No TCR": "#D1495B", "Paired TRA/TRB": "#0077B6"}
    fig, ax = plt.subplots(figsize=(6.8, 5.8))
    for group_name in ["No TCR", "Paired TRA/TRB"]:
        group_df = paired_df.loc[paired_df["tcr_pairing_group"] == group_name]
        ax.scatter(
            group_df["phase4_trab_score"],
            group_df["phase4_trd_score"],
            s=4,
            linewidths=0,
            alpha=0.72,
            rasterized=True,
            color=palette[group_name],
            label=group_name,
        )
    ax.set_title("Raw TRAB-versus-TRD Space with Paired TRA/TRB Context", fontsize=15)
    ax.set_xlabel("Raw TRAB score", fontsize=11)
    ax.set_ylabel("Raw TRD score", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(title="TRA/TRB CDR3", frameon=False, fontsize=9, title_fontsize=10)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / "phase4_trab_vs_trd_paired_tratrb_vs_no_tcr_raw_only.png",
        dpi=FIGURE_DPI,
        bbox_inches="tight",
    )
    plt.close(fig)


def write_raw_gene_panel(frame: pd.DataFrame, sample_idx: np.ndarray) -> None:
    """Write a 4x3 raw TRAB-versus-TRD gene-expression panel."""
    rng = np.random.default_rng(RANDOM_SEED)
    panel_n = min(GENE_PANEL_SAMPLE_SIZE, sample_idx.size)
    panel_pos = np.sort(rng.choice(sample_idx.size, size=panel_n, replace=False))
    panel_idx = sample_idx[panel_pos]
    panel_frame = frame.iloc[panel_pos].reset_index(drop=True)
    expr_df = extract_log1p_gene_expression_for_sample(
        INTEGRATED_H5AD,
        panel_idx,
        RAW_GENE_PANEL,
        CHUNK_SIZE,
    )
    plot_df = pd.concat([panel_frame, expr_df.reset_index(drop=True)], axis=1)

    fig, axes = plt.subplots(4, 3, figsize=(12, 15), constrained_layout=True)
    for ax, gene_name in zip(axes.flat, RAW_GENE_PANEL, strict=False):
        scatter = ax.scatter(
            plot_df["phase4_trab_score"],
            plot_df["phase4_trd_score"],
            c=plot_df[gene_name],
            cmap="viridis",
            s=3,
            linewidths=0,
            alpha=0.7,
            rasterized=True,
        )
        ax.set_title(gene_name, fontsize=11)
        ax.set_xlabel("Raw TRAB score", fontsize=9)
        ax.set_ylabel("Raw TRD score", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.03)

    fig.suptitle("Selected Marker Genes on Raw TRAB-versus-TRD Space", fontsize=16)
    fig.savefig(FIGURE_DIR / "phase4_trab_vs_trd_selected_marker_panel_raw_only.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def crop_existing_trg_panel() -> None:
    """Crop the raw portion from the existing TRG raw-vs-scaled panel."""
    trg_path = FIGURE_DIR / "phase4_trab_vs_trd_trgc_trgv9_expression.png"

    trg_image = plt.imread(trg_path)
    trg_height = trg_image.shape[0]
    trg_raw_only = trg_image[: trg_height // 2, :, :]
    plt.imsave(
        FIGURE_DIR / "phase4_trab_vs_trd_trgc_trgv9_expression_raw_only.png",
        trg_raw_only,
    )


def main() -> None:
    """Generate the raw-only Phase 4 visuals."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    frame, sample_idx = sample_obs_frame()
    tcr_frame = attach_sample_tcr_flags(frame, sample_idx)
    write_raw_umap_overlays(frame)
    write_raw_trab_trd_scatter(frame)
    write_raw_paired_tcr_scatter(tcr_frame)
    write_raw_gene_panel(frame, sample_idx)
    crop_existing_trg_panel()


if __name__ == "__main__":
    main()
