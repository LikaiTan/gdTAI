#!/usr/bin/env python3
"""Run Phase 4-style TRD/TRAB scoring on one standalone H5AD.

This wrapper reuses the package-faithful module definitions and sparse chunk
scoring logic from `phase4_gdt_module_scoring.py`, but it does not mutate the
input H5AD. It is intended for one-off review on pre-merge or external files.
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


import argparse
import logging
import time
from pathlib import Path

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from phase4_gdt_module_scoring import (
    CHUNK_SIZE,
    CTRL_SIZE,
    FIGURE_DPI,
    MARKER_GENES,
    N_BINS,
    PACKAGE_ZIP,
    PHASE4_SCORE_COLUMNS,
    PHASE4_SCALED_SCORE_COLUMNS,
    PLOT_SAMPLE_SIZE,
    RANDOM_STATE,
    SCATTER_COLOR_GENES,
    SCATTER_PLOT_SAMPLE_SIZE,
    TARGET_SUM,
    TOP_CELL_N,
    add_scaled_scores,
    append_obs_columns_in_place,
    compute_gene_means,
    compute_scores,
    downsample_indices,
    extract_log1p_gene_expression_for_sample,
    find_module_genes,
    pick_control_genes,
)


OUTPUT_ROOT = _TNK_PROJECT_ROOT / "Integrated_dataset"
DEFAULT_INPUT = _TNK_PROJECT_ROOT / "downloads/per_gse_h5ad_with_metadata/HRA005041_T_cells_subset.h5ad"
DEFAULT_ALIAS = "HRA005041"


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description="Run Phase 4-style TRD/TRAB scoring on one standalone H5AD.")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT, help="Input H5AD path.")
    parser.add_argument("--alias", type=str, default=DEFAULT_ALIAS, help="Short dataset alias for output naming.")
    parser.add_argument("--chunk-size", type=int, default=CHUNK_SIZE, help="Sparse chunk size.")
    parser.add_argument("--plot-sample-size", type=int, default=PLOT_SAMPLE_SIZE, help="Cells to sample for plots.")
    parser.add_argument(
        "--scatter-sample-size",
        type=int,
        default=SCATTER_PLOT_SAMPLE_SIZE,
        help="Cells to sample for the TRAB-vs-TRD scatter panel.",
    )
    parser.add_argument("--top-cell-n", type=int, default=TOP_CELL_N, help="Top cells to export by TRD-TRAB.")
    parser.add_argument(
        "--write-back",
        action="store_true",
        help="Append Phase 4 score columns back into the input H5AD in place.",
    )
    return parser.parse_args()


def configure_logging(log_path: Path) -> None:
    """Configure console and file logging."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    handlers = [
        logging.FileHandler(log_path, mode="w", encoding="utf-8"),
        logging.StreamHandler(),
    ]
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
        force=True,
    )


def resolve_umap(adata: ad.AnnData) -> np.ndarray | None:
    """Resolve an existing UMAP embedding from either `obsm` or `obs` columns."""
    if "X_umap" in adata.obsm:
        umap = np.asarray(adata.obsm["X_umap"])
        if umap.ndim == 2 and umap.shape[1] >= 2:
            return umap[:, :2].astype(np.float32, copy=False)
    if {"UMAP_1", "UMAP_2"}.issubset(adata.obs.columns):
        return adata.obs[["UMAP_1", "UMAP_2"]].to_numpy(dtype=np.float32, copy=True)
    return None


def resolve_leiden_labels(adata: ad.AnnData) -> np.ndarray:
    """Resolve one clustering label column for summary figures."""
    for column in ["leiden", "Cluster", "Recluster", "tnk_seurat_clusters"]:
        if column in adata.obs.columns:
            return adata.obs[column].astype(str).to_numpy()
    return np.repeat("all", adata.n_obs)


def resolve_marker_indices(var_names: pd.Index) -> tuple[np.ndarray, list[str]]:
    """Return available marker indices without failing on absent markers."""
    available_markers = [gene for gene in MARKER_GENES if gene in var_names]
    marker_idx = var_names.get_indexer(pd.Index(available_markers, dtype="string")).astype(np.int32, copy=False)
    return marker_idx, available_markers


def select_plot_sample(labels: np.ndarray, sample_size: int, random_state: int) -> np.ndarray:
    """Sample cells for plots while preserving small groups when available."""
    rng = np.random.default_rng(random_state)
    label_series = pd.Series(labels, dtype="string")
    label_counts = label_series.value_counts(sort=False)
    total_cells = int(label_counts.sum())
    target = min(sample_size, total_cells)
    sampled_idx: list[np.ndarray] = []
    label_values = label_series.to_numpy()
    for label, count in label_counts.items():
        label_idx = np.flatnonzero(label_values == label)
        if count <= 0:
            continue
        label_target = int(round(target * (count / total_cells)))
        label_target = max(500, label_target)
        label_target = min(label_target, int(count))
        sampled_idx.append(rng.choice(label_idx, size=label_target, replace=False))
    sample = np.concatenate(sampled_idx)
    if len(sample) > target:
        sample = rng.choice(sample, size=target, replace=False)
    return np.sort(sample.astype(np.int64, copy=False))


def build_output_paths(alias: str) -> dict[str, Path]:
    """Build output directories for one standalone run."""
    table_dir = OUTPUT_ROOT / "tables" / f"{alias}_phase4"
    figure_dir = OUTPUT_ROOT / "figures" / f"{alias}_phase4"
    log_path = OUTPUT_ROOT / "logs" / f"{alias}_phase4.log"
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    return {"table_dir": table_dir, "figure_dir": figure_dir, "log_path": log_path}


def write_score_summary(scores: dict[str, np.ndarray], table_dir: Path) -> pd.DataFrame:
    """Write overall score summary statistics."""
    rows = []
    score_column_map = {**PHASE4_SCORE_COLUMNS, **PHASE4_SCALED_SCORE_COLUMNS}
    for module_name, column_name in score_column_map.items():
        values = scores[module_name]
        rows.append(
            {
                "score_name": column_name,
                "n_cells": values.size,
                "mean": float(np.mean(values)),
                "std": float(np.std(values)),
                "min": float(np.min(values)),
                "p01": float(np.quantile(values, 0.01)),
                "p05": float(np.quantile(values, 0.05)),
                "median": float(np.quantile(values, 0.50)),
                "p95": float(np.quantile(values, 0.95)),
                "p99": float(np.quantile(values, 0.99)),
                "max": float(np.max(values)),
            }
        )
    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(table_dir / "phase4_score_summary.csv", index=False)
    return summary_df


def write_module_membership(modules: dict[str, pd.Index], table_dir: Path) -> None:
    """Write one table listing the genes in each module."""
    rows = []
    for module_name, genes in modules.items():
        for gene in genes:
            rows.append({"module_name": module_name, "gene_symbol": str(gene)})
    pd.DataFrame(rows).to_csv(table_dir / "phase4_module_gene_membership.csv", index=False)


def write_top_cells(
    *,
    obs_names: np.ndarray,
    sample_names: np.ndarray | None,
    donor_names: np.ndarray | None,
    scores: dict[str, np.ndarray],
    top_cell_n: int,
    table_dir: Path,
) -> None:
    """Write the top cells ranked by TRD minus TRAB."""
    top_n = min(top_cell_n, obs_names.size)
    top_idx = np.argpartition(scores["trd_minus_trab"], -top_n)[-top_n:]
    top_idx = top_idx[np.argsort(scores["trd_minus_trab"][top_idx])[::-1]]
    top_df = pd.DataFrame(
        {
            "cell_id": obs_names[top_idx],
            "sample_name": sample_names[top_idx] if sample_names is not None else "",
            "donor": donor_names[top_idx] if donor_names is not None else "",
            "phase4_tra_score": scores["tra"][top_idx],
            "phase4_trb_score": scores["trb"][top_idx],
            "phase4_trab_score": scores["trab"][top_idx],
            "phase4_trd_score": scores["trd"][top_idx],
            "phase4_trd_minus_trab": scores["trd_minus_trab"][top_idx],
            "phase4_trab_score_scaled": scores["trab_scaled"][top_idx],
            "phase4_trd_score_scaled": scores["trd_scaled"][top_idx],
            "phase4_trd_minus_trab_scaled": scores["trd_minus_trab_scaled"][top_idx],
        }
    )
    top_df.to_csv(table_dir / "phase4_top_cells_by_trd_minus_trab.csv", index=False)


def build_cluster_summary(
    *,
    labels: np.ndarray,
    scores: dict[str, np.ndarray],
    cluster_counts: np.ndarray,
    marker_detection_counts: np.ndarray,
    available_markers: list[str],
    table_dir: Path,
) -> pd.DataFrame:
    """Build one cluster-level summary table."""
    summary_df = pd.DataFrame(
        {
            "cluster": pd.Categorical(labels),
            **{column_name: scores[module_name] for module_name, column_name in {**PHASE4_SCORE_COLUMNS, **PHASE4_SCALED_SCORE_COLUMNS}.items()},
        }
    )
    cluster_summary = summary_df.groupby("cluster", observed=True).agg(
        n_cells=("cluster", "size"),
        phase4_tra_score_mean=("phase4_tra_score", "mean"),
        phase4_trb_score_mean=("phase4_trb_score", "mean"),
        phase4_trab_score_mean=("phase4_trab_score", "mean"),
        phase4_trd_score_mean=("phase4_trd_score", "mean"),
        phase4_trd_minus_trab_mean=("phase4_trd_minus_trab", "mean"),
        phase4_trd_minus_trab_median=("phase4_trd_minus_trab", "median"),
        phase4_trd_score_scaled_mean=("phase4_trd_score_scaled", "mean"),
        phase4_trab_score_scaled_mean=("phase4_trab_score_scaled", "mean"),
        phase4_trd_minus_trab_scaled_mean=("phase4_trd_minus_trab_scaled", "mean"),
        phase4_trd_minus_trab_scaled_median=("phase4_trd_minus_trab_scaled", "median"),
    ).reset_index()
    if available_markers:
        marker_fraction_df = pd.DataFrame({"cluster": cluster_summary["cluster"].astype(str)})
        cluster_categories = cluster_summary["cluster"].astype(str).to_numpy()
        for marker_pos, marker_gene in enumerate(available_markers):
            marker_fraction_df[f"{marker_gene}_detected_fraction"] = (
                marker_detection_counts[:, marker_pos] / np.maximum(cluster_counts, 1)
            )
        marker_fraction_df["cluster"] = cluster_categories
        cluster_summary["cluster"] = cluster_categories
        cluster_summary = cluster_summary.merge(marker_fraction_df, on="cluster", how="left", validate="one_to_one")
    cluster_summary = cluster_summary.sort_values("phase4_trd_minus_trab_median", ascending=False)
    cluster_summary.to_csv(table_dir / "phase4_cluster_score_summary.csv", index=False)
    return cluster_summary


def write_scored_sample_export(sample_df: pd.DataFrame, table_dir: Path) -> None:
    """Write the sampled cells used for figures."""
    sample_df.to_csv(table_dir / "phase4_plot_sample_scores.csv.gz", index=False, compression="gzip")


def paired_flag_from_obs(obs: pd.DataFrame, pair_col: str, chain_a: str, chain_b: str) -> np.ndarray:
    """Resolve one paired-TCR boolean vector from canonical obs columns."""
    if pair_col in obs.columns:
        values = obs[pair_col]
        if pd.api.types.is_bool_dtype(values):
            return values.to_numpy(dtype=bool, copy=False)
        text = values.astype(str).str.strip().str.lower()
        return text.isin({"true", "1", "yes", "y"}).to_numpy(dtype=bool, copy=False)
    return (
        obs.get(f"{chain_a}_cdr3", pd.Series("", index=obs.index)).astype(str).str.strip().ne("")
        & obs.get(f"{chain_b}_cdr3", pd.Series("", index=obs.index)).astype(str).str.strip().ne("")
    ).to_numpy(dtype=bool, copy=False)


def write_paired_tcr_scatter_panel(sample_df: pd.DataFrame, figure_dir: Path) -> None:
    """Write raw TRAB-vs-TRD scatter panels colored by paired abTCR and gdTCR."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    panel_specs = [
        ("has_TRA_TRB_paired", "Paired abTCR"),
        ("has_TRG_TRD_paired", "Paired gdTCR"),
    ]

    for ax, (column, title) in zip(axes, panel_specs):
        paired_mask = sample_df[column].to_numpy(dtype=bool, copy=False)
        ax.scatter(
            sample_df.loc[~paired_mask, "phase4_trab_score"],
            sample_df.loc[~paired_mask, "phase4_trd_score"],
            s=3,
            c="#2F6690",
            linewidths=0,
            rasterized=True,
            label="Other",
            alpha=0.65,
        )
        ax.scatter(
            sample_df.loc[paired_mask, "phase4_trab_score"],
            sample_df.loc[paired_mask, "phase4_trd_score"],
            s=3,
            c="#C1121F",
            linewidths=0,
            rasterized=True,
            label=title,
            alpha=0.75,
        )
        ax.set_title(f"{title} (n={int(paired_mask.sum()):,})")
        ax.set_xlabel("Raw TRAB score")
        ax.set_ylabel("Raw TRD score")
        ax.legend(loc="best", frameon=True)

    fig.savefig(figure_dir / "phase4_trab_vs_trd_paired_ab_gd_tcr.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def build_uns_payload(
    *,
    input_h5ad: Path,
    modules: dict[str, pd.Index],
    module_controls: dict[str, pd.Index],
    scaling_stats: dict[str, dict[str, float]],
) -> dict:
    """Build a standalone Phase 4 provenance payload for H5AD write-back."""
    return {
        "package_source": str(PACKAGE_ZIP),
        "integrated_h5ad": str(input_h5ad),
        "scoring_mode": "temporary_normalize_total_log1p_on_count_space_X",
        "continuous_only": True,
        "random_state": RANDOM_STATE,
        "target_sum": TARGET_SUM,
        "ctrl_size": CTRL_SIZE,
        "n_bins": N_BINS,
        "module_genes": {name: [str(gene) for gene in genes] for name, genes in modules.items()},
        "control_genes": {name: [str(gene) for gene in genes] for name, genes in module_controls.items()},
        "score_columns": {**PHASE4_SCORE_COLUMNS, **PHASE4_SCALED_SCORE_COLUMNS},
        "scaled_score_columns": PHASE4_SCALED_SCORE_COLUMNS,
        "scaled_score_ranges": scaling_stats,
        "scANVI_usage": "not_applicable_standalone_dataset",
    }


def write_figures(
    *,
    sample_df: pd.DataFrame,
    cluster_summary: pd.DataFrame,
    available_markers: list[str],
    figure_dir: Path,
) -> None:
    """Write standalone Phase 4-style figures."""
    sns.set_theme(style="whitegrid", context="talk")

    dist_fig, dist_axes = plt.subplots(2, 3, figsize=(18, 10), constrained_layout=True)
    dist_specs = [
        ("phase4_trd_score", "TRD score"),
        ("phase4_trab_score", "TRAB score"),
        ("phase4_trd_minus_trab", "TRD - TRAB"),
        ("phase4_trd_score_scaled", "TRD score scaled 0-1"),
        ("phase4_trab_score_scaled", "TRAB score scaled 0-1"),
        ("phase4_trd_minus_trab_scaled", "Scaled TRD - TRAB"),
    ]
    for ax, (column, title) in zip(dist_axes.flatten(), dist_specs):
        sns.histplot(sample_df[column], bins=100, ax=ax, color="#2F6690")
        ax.set_title(title)
        ax.set_xlabel(column)
    dist_fig.savefig(figure_dir / "phase4_score_distributions.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(dist_fig)

    if {"umap1", "umap2"}.issubset(sample_df.columns):
        umap_fig, umap_axes = plt.subplots(2, 3, figsize=(18, 10), constrained_layout=True)
        for ax, (column, title) in zip(umap_axes.flatten(), dist_specs):
            scatter = ax.scatter(
                sample_df["umap1"],
                sample_df["umap2"],
                c=sample_df[column],
                cmap="viridis",
                s=3,
                linewidths=0,
                rasterized=True,
            )
            ax.set_title(title)
            ax.set_xlabel("UMAP1")
            ax.set_ylabel("UMAP2")
            umap_fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
        umap_fig.savefig(figure_dir / "phase4_umap_score_overlays.png", dpi=FIGURE_DPI, bbox_inches="tight")
        plt.close(umap_fig)

    if cluster_summary["cluster"].nunique() > 1:
        heatmap_cols = [
            "phase4_tra_score_mean",
            "phase4_trb_score_mean",
            "phase4_trab_score_mean",
            "phase4_trd_score_mean",
            "phase4_trd_minus_trab_median",
            "phase4_trd_score_scaled_mean",
            "phase4_trab_score_scaled_mean",
            "phase4_trd_minus_trab_scaled_median",
        ]
        heatmap_df = cluster_summary.set_index("cluster")[heatmap_cols]
        heatmap_fig, heatmap_ax = plt.subplots(figsize=(10, 10), constrained_layout=True)
        sns.heatmap(heatmap_df, cmap="vlag", center=0, ax=heatmap_ax)
        heatmap_ax.set_title("Phase 4 Score Summary by Cluster")
        heatmap_fig.savefig(figure_dir / "phase4_cluster_score_summary.png", dpi=FIGURE_DPI, bbox_inches="tight")
        plt.close(heatmap_fig)

        if available_markers:
            marker_cols = [
                "phase4_trd_score_mean",
                "phase4_trab_score_mean",
                "phase4_trd_minus_trab_median",
            ] + [f"{marker}_detected_fraction" for marker in available_markers]
            marker_df = cluster_summary.set_index("cluster")[marker_cols]
            marker_fig, marker_ax = plt.subplots(figsize=(12, 10), constrained_layout=True)
            sns.heatmap(marker_df, cmap="rocket", ax=marker_ax)
            marker_ax.set_title("Phase 4 Scores and TCR Marker Detection by Cluster")
            marker_fig.savefig(figure_dir / "phase4_marker_score_comparison.png", dpi=FIGURE_DPI, bbox_inches="tight")
            plt.close(marker_fig)


def write_trab_trd_scatter_panel(sample_df: pd.DataFrame, figure_dir: Path) -> None:
    """Write the TRAB-vs-TRD scatter panel for the selected sample."""
    color_genes = [gene for gene in SCATTER_COLOR_GENES if gene in sample_df.columns]
    raw_color_specs = [("phase4_trab_minus_trd", "TRAB - TRD")] + [(gene, gene) for gene in color_genes]
    scaled_color_specs = [("phase4_trab_minus_trd_scaled", "Scaled TRAB - TRD")] + [(gene, gene) for gene in color_genes]
    fig, axes = plt.subplots(2, len(raw_color_specs), figsize=(4 * len(raw_color_specs), 8), constrained_layout=True)

    def draw_row(row_axes: np.ndarray, x_col: str, y_col: str, color_specs: list[tuple[str, str]], row_title: str) -> None:
        for ax, (color_col, title) in zip(row_axes, color_specs):
            cmap = "coolwarm" if "TRAB - TRD" in title else "viridis"
            scatter = ax.scatter(
                sample_df[x_col],
                sample_df[y_col],
                c=sample_df[color_col],
                cmap=cmap,
                s=3,
                linewidths=0,
                rasterized=True,
            )
            ax.set_title(title)
            ax.set_xlabel(x_col)
            ax.set_ylabel(y_col)
            fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
        row_axes[0].annotate(
            row_title,
            xy=(-0.38, 0.5),
            xycoords="axes fraction",
            rotation=90,
            va="center",
            ha="center",
            fontsize=16,
            fontweight="bold",
        )

    draw_row(axes[0], "phase4_trab_score", "phase4_trd_score", raw_color_specs, "Raw scores")
    draw_row(
        axes[1],
        "phase4_trab_score_scaled",
        "phase4_trd_score_scaled",
        scaled_color_specs,
        "Scaled scores",
    )
    fig.savefig(figure_dir / "phase4_trab_vs_trd_scatter_panel.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def write_qc_summary(
    *,
    alias: str,
    input_h5ad: Path,
    modules: dict[str, pd.Index],
    available_markers: list[str],
    score_summary: pd.DataFrame,
    cluster_summary: pd.DataFrame,
    output_paths: dict[str, Path],
    elapsed_seconds: float,
) -> None:
    """Write one markdown QC summary for the standalone run."""
    summary_path = output_paths["table_dir"] / "phase4_qc_summary.md"
    top_cluster = cluster_summary.iloc[0]
    lines = [
        f"# Phase 4 Standalone QC Summary: {alias}",
        "",
        f"- Input H5AD: `{input_h5ad}`",
        f"- Cells: `{int(score_summary['n_cells'].iloc[0]):,}`",
        f"- Module gene counts: TRA={len(modules['tra'])}, TRB={len(modules['trb'])}, TRAB={len(modules['trab'])}, TRD={len(modules['trd'])}",
        f"- Available marker genes for cluster QC: `{', '.join(available_markers) if available_markers else 'none'}`",
        f"- Top cluster by median TRD-TRAB: `{top_cluster['cluster']}` with median `{float(top_cluster['phase4_trd_minus_trab_median']):.4f}`",
        f"- Runtime seconds: `{elapsed_seconds:.1f}`",
        "",
        "## Output files",
        "",
        f"- Tables: `{output_paths['table_dir']}`",
        f"- Figures: `{output_paths['figure_dir']}`",
        f"- Log: `{output_paths['log_path']}`",
    ]
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """Run standalone Phase 4 scoring on one H5AD without mutating it."""
    args = parse_args()
    output_paths = build_output_paths(args.alias)
    configure_logging(output_paths["log_path"])
    start_time = time.time()

    input_h5ad = args.input_h5ad.resolve()
    logging.info("Using standalone H5AD at %s", input_h5ad)
    adata = ad.read_h5ad(input_h5ad, backed="r")
    n_obs = adata.n_obs
    n_vars = adata.n_vars
    obs_names = adata.obs_names.to_numpy(dtype=str, copy=True)
    var_names = pd.Index(adata.var_names.astype(str), dtype="string")
    obs_frame = adata.obs.copy()
    leiden_labels = resolve_leiden_labels(adata)
    umap = resolve_umap(adata)
    sample_names = (
        obs_frame["sample_id"].astype(str).to_numpy()
        if "sample_id" in obs_frame.columns
        else (obs_frame["Sample_name"].astype(str).to_numpy() if "Sample_name" in obs_frame.columns else None)
    )
    donor_names = (
        obs_frame["donor_id"].astype(str).to_numpy()
        if "donor_id" in obs_frame.columns
        else (obs_frame["Donor"].astype(str).to_numpy() if "Donor" in obs_frame.columns else None)
    )
    paired_ab_flags = paired_flag_from_obs(obs_frame, "has_TRA_TRB_paired", "TRA", "TRB")
    paired_gd_flags = paired_flag_from_obs(obs_frame, "has_TRG_TRD_paired", "TRG", "TRD")
    adata.file.close()

    logging.info("Input matrix shape: n_obs=%s, n_vars=%s", n_obs, n_vars)

    modules = find_module_genes(var_names)
    for module_name, genes in modules.items():
        logging.info("Module %s contains %s genes", module_name, len(genes))

    gene_means = compute_gene_means(input_h5ad, n_obs, n_vars, args.chunk_size)
    gene_pool = pd.Index(var_names, dtype="string")
    module_controls = {
        module_name: pick_control_genes(gene_list, gene_pool, gene_means, random_state=RANDOM_STATE)
        for module_name, gene_list in modules.items()
    }
    module_gene_idx = {name: var_names.get_indexer(genes).astype(np.int32) for name, genes in modules.items()}
    module_ctrl_idx = {
        name: var_names.get_indexer(ctrl_genes).astype(np.int32) for name, ctrl_genes in module_controls.items()
    }
    marker_idx, available_markers = resolve_marker_indices(var_names)

    label_codes, _ = pd.factorize(leiden_labels, sort=True)
    scores, cluster_counts, marker_detection_counts = compute_scores(
        input_h5ad,
        n_obs=n_obs,
        n_vars=n_vars,
        chunk_size=args.chunk_size,
        module_gene_idx=module_gene_idx,
        module_ctrl_idx=module_ctrl_idx,
        leiden_codes=label_codes.astype(np.int32, copy=False),
        marker_idx=marker_idx,
    )
    scores, scaling_stats = add_scaled_scores(scores)

    score_summary = write_score_summary(scores, output_paths["table_dir"])
    write_module_membership(modules, output_paths["table_dir"])
    cluster_summary = build_cluster_summary(
        labels=leiden_labels,
        scores=scores,
        cluster_counts=cluster_counts,
        marker_detection_counts=marker_detection_counts,
        available_markers=available_markers,
        table_dir=output_paths["table_dir"],
    )
    write_top_cells(
        obs_names=obs_names,
        sample_names=sample_names,
        donor_names=donor_names,
        scores=scores,
        top_cell_n=args.top_cell_n,
        table_dir=output_paths["table_dir"],
    )

    sample_idx = select_plot_sample(leiden_labels, args.plot_sample_size, RANDOM_STATE)
    sample_df = pd.DataFrame(
        {
            "cell_id": obs_names[sample_idx],
            "cluster": leiden_labels[sample_idx],
            "phase4_trd_score": scores["trd"][sample_idx],
            "phase4_trab_score": scores["trab"][sample_idx],
            "phase4_trd_minus_trab": scores["trd_minus_trab"][sample_idx],
            "phase4_trd_score_scaled": scores["trd_scaled"][sample_idx],
            "phase4_trab_score_scaled": scores["trab_scaled"][sample_idx],
            "phase4_trd_minus_trab_scaled": scores["trd_minus_trab_scaled"][sample_idx],
            "has_TRA_TRB_paired": paired_ab_flags[sample_idx],
            "has_TRG_TRD_paired": paired_gd_flags[sample_idx],
        }
    )
    sample_df["phase4_trab_minus_trd"] = sample_df["phase4_trab_score"] - sample_df["phase4_trd_score"]
    sample_df["phase4_trab_minus_trd_scaled"] = (
        sample_df["phase4_trab_score_scaled"] - sample_df["phase4_trd_score_scaled"]
    )
    if umap is not None:
        sample_df["umap1"] = umap[sample_idx, 0]
        sample_df["umap2"] = umap[sample_idx, 1]

    scatter_idx = downsample_indices(sample_idx, args.scatter_sample_size, RANDOM_STATE)
    scatter_selector = np.searchsorted(sample_idx, scatter_idx)
    scatter_df = sample_df.iloc[scatter_selector].reset_index(drop=True)
    available_color_genes = [gene for gene in SCATTER_COLOR_GENES if gene in var_names]
    if available_color_genes:
        expr_df = extract_log1p_gene_expression_for_sample(
            input_h5ad,
            scatter_idx,
            available_color_genes,
            args.chunk_size,
        )
        scatter_df = pd.concat([scatter_df, expr_df], axis=1)

    write_scored_sample_export(sample_df, output_paths["table_dir"])
    write_figures(
        sample_df=sample_df,
        cluster_summary=cluster_summary,
        available_markers=available_markers,
        figure_dir=output_paths["figure_dir"],
    )
    write_trab_trd_scatter_panel(scatter_df, output_paths["figure_dir"])
    write_paired_tcr_scatter_panel(scatter_df, output_paths["figure_dir"])

    if args.write_back:
        append_obs_columns_in_place(
            input_h5ad,
            {
                **{column_name: scores[module_name] for module_name, column_name in PHASE4_SCORE_COLUMNS.items()},
                **{column_name: scores[module_name] for module_name, column_name in PHASE4_SCALED_SCORE_COLUMNS.items()},
            },
            build_uns_payload(
                input_h5ad=input_h5ad,
                modules=modules,
                module_controls=module_controls,
                scaling_stats=scaling_stats,
            ),
        )
        logging.info("Appended Phase 4 score columns back into %s", input_h5ad)

    elapsed_seconds = time.time() - start_time
    write_qc_summary(
        alias=args.alias,
        input_h5ad=input_h5ad,
        modules=modules,
        available_markers=available_markers,
        score_summary=score_summary,
        cluster_summary=cluster_summary,
        output_paths=output_paths,
        elapsed_seconds=elapsed_seconds,
    )
    logging.info("Standalone Phase 4 run for %s completed in %.1f seconds", args.alias, elapsed_seconds)


if __name__ == "__main__":
    main()
