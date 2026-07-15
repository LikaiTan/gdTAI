#!/usr/bin/env python3
"""Run Leiden DEG and gamma-delta T-marker enrichment analysis.

This helper keeps the current integrated milestone unchanged and writes a
post-Phase-4 cluster-analysis package to NFS. The analysis has two layers:

1. one-vs-rest DEG on a balanced per-cluster sample of the integrated object
2. all-cell mean/detection summaries for TRDC/TRDV1/TRDV2/TRGV9

The DEG layer is sampling-based for practicality at project scale. The focused
gamma-delta marker enrichment layer is computed on all cells.
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


import logging
from pathlib import Path

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse

from phase4_gdt_module_scoring import TARGET_SUM, build_csr_chunk, normalize_log1p_chunk


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
TABLE_DIR = OUTPUT_ROOT / "tables"
FIGURE_DIR = OUTPUT_ROOT / "figures"
LOG_DIR = OUTPUT_ROOT / "logs"

INTEGRATED_H5AD = PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "integrated.h5ad"
PHASE4_LEIDEN_SCORE_CSV = TABLE_DIR / "phase4_leiden_score_summary.csv"

DEG_FULL_CSV = TABLE_DIR / "phase3_leiden_deg_full.csv"
DEG_TOP_CSV = TABLE_DIR / "phase3_leiden_deg_top_markers.csv"
GDT_ENRICHMENT_CSV = TABLE_DIR / "phase3_leiden_gdt_marker_enrichment.csv"
CLUSTER_PRIORITY_CSV = TABLE_DIR / "phase3_leiden_gdt_cluster_priority.csv"
SUMMARY_MD = LOG_DIR / "phase3_leiden_deg_summary.md"
RUN_LOG = LOG_DIR / "phase3_leiden_deg.log"

TOP_MARKER_HEATMAP_PNG = FIGURE_DIR / "phase3_leiden_deg_top_markers_heatmap.png"
GDT_EXPR_PNG = FIGURE_DIR / "phase3_gdt_marker_expression_by_cluster.png"
GDT_DETECT_PNG = FIGURE_DIR / "phase3_gdt_marker_detection_by_cluster.png"
GDT_ENRICHED_UMAP_PNG = FIGURE_DIR / "phase3_gdt_enriched_clusters_on_umap.png"

FOCUS_GENES = ["TRDC", "TRDV1", "TRDV2", "TRGV9"]
CHUNK_SIZE = 50_000
MAX_DEG_CELLS_PER_CLUSTER = 5_000
TOP_MARKERS_PER_CLUSTER = 25
HEATMAP_TOP_MARKERS_PER_CLUSTER = 3
HEATMAP_MAX_GENES = 40
MIN_ADJ_P = 1e-10
MIN_Z = 1.0
FIGURE_DPI = 300
RANDOM_SEED = 1
UMAP_BACKGROUND_MAX = 120_000
UMAP_MAX_PER_ENRICHED_CLUSTER = 12_000


def configure_logging() -> None:
    """Configure file and console logging."""
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(RUN_LOG, mode="a", encoding="utf-8"),
            logging.StreamHandler(),
        ],
        force=True,
    )


def ensure_output_dirs() -> None:
    """Ensure NFS output directories exist."""
    for path in [TABLE_DIR, FIGURE_DIR, LOG_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def cluster_sort_key(value: object) -> tuple[int, str]:
    """Sort Leiden labels numerically when possible."""
    text = str(value)
    return (0, f"{int(text):06d}") if text.isdigit() else (1, text)


def read_string_column(obs_group: h5py.Group, column: str) -> np.ndarray:
    """Read a string obs column as a plain object array."""
    values = obs_group[column][:]
    return np.asarray(
        [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values],
        dtype=object,
    )


def load_obs_frame() -> pd.DataFrame:
    """Load only the obs fields needed for DEG sampling and cluster summaries."""
    with h5py.File(INTEGRATED_H5AD, "r") as handle:
        obs = handle["obs"]
        frame = pd.DataFrame(
            {
                "cell_id": read_string_column(obs, "_index"),
                "leiden": read_string_column(obs, "leiden"),
            }
        )
    frame["leiden"] = frame["leiden"].astype(str)
    return frame


def sample_cells_by_cluster(obs_frame: pd.DataFrame) -> np.ndarray:
    """Build a balanced per-cluster sample for one-vs-rest DEG."""
    rng = np.random.default_rng(RANDOM_SEED)
    sampled: list[np.ndarray] = []
    leiden_series = obs_frame["leiden"]
    for cluster in sorted(leiden_series.unique(), key=cluster_sort_key):
        idx = np.flatnonzero(leiden_series.to_numpy() == cluster)
        if idx.size <= MAX_DEG_CELLS_PER_CLUSTER:
            sampled.append(idx)
        else:
            sampled.append(np.sort(rng.choice(idx, size=MAX_DEG_CELLS_PER_CLUSTER, replace=False)))
    sample_idx = np.concatenate(sampled)
    sample_idx.sort()
    logging.info("Balanced DEG sample includes %s cells across %s clusters", sample_idx.size, len(sampled))
    return sample_idx


def load_sampled_anndata(obs_frame: pd.DataFrame, sample_idx: np.ndarray) -> ad.AnnData:
    """Read the balanced DEG sample from the integrated CSR matrix."""
    logging.info("Loading balanced DEG sample from integrated H5AD")
    with h5py.File(INTEGRATED_H5AD, "r") as handle:
        var_names = pd.Index(handle["var"]["_index"].asstr()[:], dtype="object")
        n_obs = int(handle["obs"]["_index"].shape[0])
        n_vars = int(var_names.size)
        x_group = handle["X"]
        sample_mask = np.zeros(n_obs, dtype=bool)
        sample_mask[sample_idx] = True
        blocks: list[sparse.csr_matrix] = []
        kept_rows = 0
        for start in range(0, n_obs, CHUNK_SIZE):
            end = min(start + CHUNK_SIZE, n_obs)
            chunk = build_csr_chunk(x_group, start, end, n_vars)
            chunk_mask = sample_mask[start:end]
            if np.any(chunk_mask):
                selected = chunk[chunk_mask, :]
                blocks.append(selected)
                kept_rows += selected.shape[0]
            if start == 0 or end == n_obs or (start // CHUNK_SIZE) % 10 == 0:
                logging.info("Loaded DEG chunk %s-%s / %s", start, end, n_obs)
    if kept_rows != sample_idx.size:
        raise ValueError(f"Expected {sample_idx.size} sampled cells but loaded {kept_rows}.")
    sampled_x = sparse.vstack(blocks, format="csr")
    sampled_obs = obs_frame.iloc[sample_idx].reset_index(drop=True).copy()
    sampled_obs.index = sampled_obs["cell_id"].astype(str)
    sampled_var = pd.DataFrame(index=var_names.astype(str))
    return ad.AnnData(X=sampled_x, obs=sampled_obs[["leiden"]], var=sampled_var)


def run_deg(adata_deg: ad.AnnData) -> pd.DataFrame:
    """Run one-vs-rest Wilcoxon DEG on the sampled AnnData."""
    logging.info("Normalizing balanced DEG sample")
    sc.pp.normalize_total(adata_deg, target_sum=TARGET_SUM)
    sc.pp.log1p(adata_deg)
    adata_deg.obs["leiden"] = pd.Categorical(
        adata_deg.obs["leiden"].astype(str),
        categories=sorted(adata_deg.obs["leiden"].astype(str).unique(), key=cluster_sort_key),
    )

    logging.info("Running rank_genes_groups on sampled DEG matrix")
    sc.tl.rank_genes_groups(
        adata_deg,
        groupby="leiden",
        method="wilcoxon",
        reference="rest",
        pts=True,
        tie_correct=True,
        key_added="phase3_leiden_deg",
        n_genes=adata_deg.n_vars,
    )

    groups = [str(group) for group in adata_deg.uns["phase3_leiden_deg"]["names"].dtype.names]
    frames: list[pd.DataFrame] = []
    for group in groups:
        frame = sc.get.rank_genes_groups_df(adata_deg, group=group, key="phase3_leiden_deg")
        frame.insert(0, "leiden", group)
        frame.insert(1, "rank", np.arange(1, frame.shape[0] + 1, dtype=int))
        frames.append(frame)
    deg = pd.concat(frames, ignore_index=True)
    deg = deg.rename(
        columns={
            "names": "gene",
            "scores": "score",
            "logfoldchanges": "logfoldchange",
            "pvals": "p_value",
            "pvals_adj": "p_value_adj",
            "pct_nz_group": "pct_nz_cluster",
            "pct_nz_reference": "pct_nz_rest",
        }
    )
    return deg


def compute_focus_marker_summary(obs_frame: pd.DataFrame) -> pd.DataFrame:
    """Compute all-cell cluster mean and detection for the focused gamma-delta genes."""
    logging.info("Computing all-cell cluster summaries for %s", ", ".join(FOCUS_GENES))
    with h5py.File(INTEGRATED_H5AD, "r") as handle:
        var_names = pd.Index(handle["var"]["_index"].asstr()[:], dtype="object")
        missing = [gene for gene in FOCUS_GENES if gene not in var_names]
        if missing:
            raise KeyError(f"Missing focus genes in integrated milestone: {missing}")
        gene_idx = np.asarray([int(var_names.get_loc(gene)) for gene in FOCUS_GENES], dtype=np.int32)
        n_obs = int(handle["obs"]["_index"].shape[0])
        n_vars = int(var_names.size)
        cluster_order = sorted(obs_frame["leiden"].unique(), key=cluster_sort_key)
        cluster_to_code = {cluster: code for code, cluster in enumerate(cluster_order)}
        cluster_codes = obs_frame["leiden"].map(cluster_to_code).to_numpy(dtype=np.int32, copy=False)
        cluster_counts = np.bincount(cluster_codes, minlength=len(cluster_order)).astype(np.int64)
        expr_sum = np.zeros((len(cluster_order), len(FOCUS_GENES)), dtype=np.float64)
        detect_sum = np.zeros((len(cluster_order), len(FOCUS_GENES)), dtype=np.int64)

        x_group = handle["X"]
        for start in range(0, n_obs, CHUNK_SIZE):
            end = min(start + CHUNK_SIZE, n_obs)
            chunk = build_csr_chunk(x_group, start, end, n_vars)
            normalized = normalize_log1p_chunk(chunk.copy(), TARGET_SUM)
            normalized_focus = normalized[:, gene_idx]
            detected_focus = chunk[:, gene_idx].copy()
            detected_focus.data[:] = 1.0
            codes = cluster_codes[start:end]
            for code in np.unique(codes):
                mask = codes == code
                expr_sum[code] += np.asarray(normalized_focus[mask, :].sum(axis=0)).ravel()
                detect_sum[code] += np.asarray(detected_focus[mask, :].sum(axis=0)).ravel().astype(np.int64)
            if start == 0 or end == n_obs or (start // CHUNK_SIZE) % 10 == 0:
                logging.info("Processed focus-gene chunk %s-%s / %s", start, end, n_obs)

    records: list[dict[str, object]] = []
    safe_counts = np.maximum(cluster_counts, 1)
    for cluster, code in cluster_to_code.items():
        for gene_pos, gene in enumerate(FOCUS_GENES):
            records.append(
                {
                    "leiden": cluster,
                    "gene": gene,
                    "cluster_cells": int(cluster_counts[code]),
                    "mean_log1p_expr": float(expr_sum[code, gene_pos] / safe_counts[code]),
                    "detected_fraction": float(detect_sum[code, gene_pos] / safe_counts[code]),
                }
            )
    return pd.DataFrame(records)


def build_enrichment_tables(deg: pd.DataFrame, focus_summary: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Merge DEG plus all-cell marker summaries and compute enrichment calls."""
    focus_deg = deg.loc[deg["gene"].isin(FOCUS_GENES), ["leiden", "gene", "rank", "score", "logfoldchange", "p_value_adj"]].copy()
    merged = focus_summary.merge(focus_deg, on=["leiden", "gene"], how="left", validate="one_to_one")

    merged["mean_expr_rank"] = merged.groupby("gene")["mean_log1p_expr"].rank(method="dense", ascending=False)
    merged["detected_fraction_rank"] = merged.groupby("gene")["detected_fraction"].rank(method="dense", ascending=False)
    merged["mean_expr_z"] = merged.groupby("gene")["mean_log1p_expr"].transform(
        lambda series: (series - series.mean()) / (series.std(ddof=0) if float(series.std(ddof=0)) > 0 else 1.0)
    )
    merged["detected_fraction_z"] = merged.groupby("gene")["detected_fraction"].transform(
        lambda series: (series - series.mean()) / (series.std(ddof=0) if float(series.std(ddof=0)) > 0 else 1.0)
    )
    merged["deg_positive"] = merged["logfoldchange"].fillna(0.0) > 0
    merged["deg_significant"] = merged["p_value_adj"].fillna(1.0) <= MIN_ADJ_P
    merged["marker_enriched"] = (
        merged["deg_positive"]
        & merged["deg_significant"]
        & (merged["mean_expr_z"] >= MIN_Z)
        & (merged["detected_fraction_z"] >= MIN_Z)
    )

    phase4_leiden = pd.read_csv(PHASE4_LEIDEN_SCORE_CSV, dtype={"leiden": str})
    cluster_priority = (
        merged.groupby("leiden", as_index=False)
        .agg(
            cluster_cells=("cluster_cells", "first"),
            enriched_markers_n=("marker_enriched", "sum"),
            mean_expr_z_sum=("mean_expr_z", "sum"),
            detected_fraction_z_sum=("detected_fraction_z", "sum"),
        )
        .assign(priority_score=lambda df: 10.0 * df["enriched_markers_n"] + df["mean_expr_z_sum"] + df["detected_fraction_z_sum"])
    )
    cluster_priority = cluster_priority.merge(
        phase4_leiden[
            [
                "leiden",
                "phase4_trd_score_mean",
                "phase4_trab_score_mean",
                "phase4_trd_minus_trab_mean",
                "phase4_trd_minus_trab_median",
                "TRDC_detected_fraction",
            ]
        ],
        on="leiden",
        how="left",
        validate="one_to_one",
    )
    cluster_priority = cluster_priority.sort_values(
        ["enriched_markers_n", "priority_score", "cluster_cells"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    return merged, cluster_priority


def build_top_marker_table(deg: pd.DataFrame) -> pd.DataFrame:
    """Extract the top positive markers per cluster."""
    positive = deg.loc[deg["logfoldchange"] > 0].copy()
    positive = positive.sort_values(["leiden", "rank"], ascending=[True, True])
    return positive.groupby("leiden", as_index=False, group_keys=False).head(TOP_MARKERS_PER_CLUSTER).reset_index(drop=True)


def plot_top_marker_heatmap(deg_sample: ad.AnnData, top_markers: pd.DataFrame) -> None:
    """Plot a cluster-by-marker heatmap from sampled DEG data."""
    selected_genes: list[str] = []
    for cluster in sorted(top_markers["leiden"].astype(str).unique(), key=cluster_sort_key):
        genes = top_markers.loc[top_markers["leiden"].astype(str) == cluster, "gene"].astype(str).tolist()
        for gene in genes[:HEATMAP_TOP_MARKERS_PER_CLUSTER]:
            if gene not in selected_genes:
                selected_genes.append(gene)
            if len(selected_genes) >= HEATMAP_MAX_GENES:
                break
        if len(selected_genes) >= HEATMAP_MAX_GENES:
            break
    matrix = pd.DataFrame(index=deg_sample.obs["leiden"].cat.categories, columns=selected_genes, dtype=float)
    for cluster in matrix.index:
        cluster_view = deg_sample[deg_sample.obs["leiden"] == cluster, :]
        for gene in selected_genes:
            gene_idx = int(cluster_view.var_names.get_loc(gene))
            matrix.loc[cluster, gene] = float(np.asarray(cluster_view.X[:, gene_idx].mean()).ravel()[0])
    sns.set_theme(style="white")
    fig, ax = plt.subplots(figsize=(max(12, 0.4 * len(selected_genes)), 10), constrained_layout=True)
    sns.heatmap(matrix.astype(float), cmap="rocket", ax=ax)
    ax.set_title("Phase 3 Leiden Top DEG Marker Heatmap")
    ax.set_xlabel("Top sampled DEG markers")
    ax.set_ylabel("Leiden cluster")
    fig.savefig(TOP_MARKER_HEATMAP_PNG, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_focus_heatmaps(enrichment: pd.DataFrame, cluster_priority: pd.DataFrame) -> None:
    """Plot cluster-level expression and detection heatmaps for focused markers."""
    ordered_clusters = cluster_priority["leiden"].astype(str).tolist()
    expr = (
        enrichment.pivot(index="leiden", columns="gene", values="mean_log1p_expr")
        .reindex(ordered_clusters)
        .loc[:, FOCUS_GENES]
    )
    detect = (
        enrichment.pivot(index="leiden", columns="gene", values="detected_fraction")
        .reindex(ordered_clusters)
        .loc[:, FOCUS_GENES]
    )

    sns.set_theme(style="white")
    fig_expr, ax_expr = plt.subplots(figsize=(8, 10), constrained_layout=True)
    sns.heatmap(expr.astype(float), cmap="magma", ax=ax_expr)
    ax_expr.set_title("Gamma-delta Marker Mean Expression by Leiden Cluster")
    ax_expr.set_xlabel("Marker gene")
    ax_expr.set_ylabel("Leiden cluster")
    fig_expr.savefig(GDT_EXPR_PNG, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig_expr)

    fig_detect, ax_detect = plt.subplots(figsize=(8, 10), constrained_layout=True)
    sns.heatmap(detect.astype(float), cmap="viridis", ax=ax_detect)
    ax_detect.set_title("Gamma-delta Marker Detection Fraction by Leiden Cluster")
    ax_detect.set_xlabel("Marker gene")
    ax_detect.set_ylabel("Leiden cluster")
    fig_detect.savefig(GDT_DETECT_PNG, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig_detect)


def plot_enriched_clusters_umap(obs_frame: pd.DataFrame, cluster_priority: pd.DataFrame) -> None:
    """Highlight enriched clusters on the integrated UMAP embedding."""
    enriched_clusters = set(cluster_priority.loc[cluster_priority["enriched_markers_n"] > 0, "leiden"].astype(str))
    if not enriched_clusters:
        logging.warning("No enriched clusters were called; plotting only the background UMAP sample.")
    rng = np.random.default_rng(RANDOM_SEED)
    leiden = obs_frame["leiden"].astype(str).to_numpy()
    enriched_mask = np.isin(leiden, list(enriched_clusters))
    background_idx = np.flatnonzero(~enriched_mask)
    if background_idx.size > UMAP_BACKGROUND_MAX:
        background_idx = np.sort(rng.choice(background_idx, size=UMAP_BACKGROUND_MAX, replace=False))
    highlight_parts: list[np.ndarray] = []
    for cluster in sorted(enriched_clusters, key=cluster_sort_key):
        cluster_idx = np.flatnonzero(leiden == cluster)
        if cluster_idx.size > UMAP_MAX_PER_ENRICHED_CLUSTER:
            cluster_idx = np.sort(rng.choice(cluster_idx, size=UMAP_MAX_PER_ENRICHED_CLUSTER, replace=False))
        highlight_parts.append(cluster_idx)
    highlight_idx = np.concatenate(highlight_parts) if highlight_parts else np.empty(0, dtype=np.int64)
    use_idx = np.unique(np.concatenate([background_idx, highlight_idx]))

    with h5py.File(INTEGRATED_H5AD, "r") as handle:
        umap = np.asarray(handle["obsm"]["X_umap"][use_idx, :], dtype=np.float32)
    plot_frame = pd.DataFrame(
        {
            "UMAP1": umap[:, 0],
            "UMAP2": umap[:, 1],
            "leiden": leiden[use_idx],
        }
    )
    plot_frame["is_enriched_cluster"] = plot_frame["leiden"].isin(enriched_clusters)

    fig, ax = plt.subplots(figsize=(9, 8), constrained_layout=True)
    background = plot_frame.loc[~plot_frame["is_enriched_cluster"]]
    highlighted = plot_frame.loc[plot_frame["is_enriched_cluster"]]
    ax.scatter(background["UMAP1"], background["UMAP2"], s=2, c="#c9c9c9", alpha=0.35, linewidths=0)
    palette = sns.color_palette("tab10", n_colors=max(highlighted["leiden"].nunique(), 1))
    for color, cluster in zip(palette, sorted(highlighted["leiden"].astype(str).unique(), key=cluster_sort_key), strict=False):
        cluster_frame = highlighted.loc[highlighted["leiden"].astype(str) == cluster]
        ax.scatter(cluster_frame["UMAP1"], cluster_frame["UMAP2"], s=4, alpha=0.85, linewidths=0, label=f"Leiden {cluster}", color=color)
    ax.set_title("Leiden UMAP Highlighting Gamma-delta-Enriched Clusters")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    if highlighted["leiden"].nunique() > 0:
        ax.legend(markerscale=3, fontsize=9, frameon=False, loc="best")
    fig.savefig(GDT_ENRICHED_UMAP_PNG, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def write_summary(top_markers: pd.DataFrame, enrichment: pd.DataFrame, cluster_priority: pd.DataFrame, sample_n: int) -> None:
    """Write a concise Markdown summary of the cluster DEG package."""
    per_gene_top: list[str] = []
    for gene in FOCUS_GENES:
        gene_frame = enrichment.loc[enrichment["gene"] == gene].sort_values(
            ["marker_enriched", "mean_log1p_expr", "detected_fraction"],
            ascending=[False, False, False],
        )
        top = gene_frame.iloc[0]
        per_gene_top.append(
            f"- `{gene}`: top cluster `Leiden {top['leiden']}` "
            f"(mean={top['mean_log1p_expr']:.3f}, detected_fraction={top['detected_fraction']:.3f}, "
            f"deg_rank={int(top['rank']) if pd.notna(top['rank']) else 'NA'}, enriched={bool(top['marker_enriched'])})"
        )

    priority_preview = cluster_priority.head(10)
    lines = [
        "# Phase 3 Leiden DEG and Gamma-delta Marker Enrichment",
        "",
        f"- Integrated source: `{INTEGRATED_H5AD}`",
        f"- DEG mode: balanced one-vs-rest Wilcoxon on `{sample_n:,}` sampled cells",
        f"- All-cell gamma-delta marker summaries were computed for `{', '.join(FOCUS_GENES)}`",
        f"- Leiden clusters reviewed: `{cluster_priority.shape[0]}`",
        "",
        "## Top cluster per focused marker",
        "",
        *per_gene_top,
        "",
        "## Highest-priority clusters",
        "",
    ]
    for _, row in priority_preview.iterrows():
        lines.append(
            f"- `Leiden {row['leiden']}`: enriched_markers_n={int(row['enriched_markers_n'])}, "
            f"priority_score={row['priority_score']:.3f}, "
            f"phase4_trd_minus_trab_mean={row['phase4_trd_minus_trab_mean']:.3f}, "
            f"cells={int(row['cluster_cells']):,}"
        )
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """Run the full Leiden DEG and gamma-delta marker enrichment package."""
    configure_logging()
    ensure_output_dirs()

    obs_frame = load_obs_frame()
    sample_idx = sample_cells_by_cluster(obs_frame)
    adata_deg = load_sampled_anndata(obs_frame, sample_idx)
    deg = run_deg(adata_deg)
    top_markers = build_top_marker_table(deg)
    focus_summary = compute_focus_marker_summary(obs_frame)
    enrichment, cluster_priority = build_enrichment_tables(deg, focus_summary)

    logging.info("Writing cluster-analysis tables")
    deg.to_csv(DEG_FULL_CSV, index=False)
    top_markers.to_csv(DEG_TOP_CSV, index=False)
    enrichment.to_csv(GDT_ENRICHMENT_CSV, index=False)
    cluster_priority.to_csv(CLUSTER_PRIORITY_CSV, index=False)

    logging.info("Writing cluster-analysis figures")
    plot_top_marker_heatmap(adata_deg, top_markers)
    plot_focus_heatmaps(enrichment, cluster_priority)
    plot_enriched_clusters_umap(obs_frame, cluster_priority)
    write_summary(top_markers, enrichment, cluster_priority, sample_n=int(sample_idx.size))
    logging.info("Done")


if __name__ == "__main__":
    main()
