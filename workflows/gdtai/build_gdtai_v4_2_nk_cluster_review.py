#!/usr/bin/env python3
"""Build marker, annotation, and unsupervised-cluster review of V4.2 NK cells.

This workflow is read-only with respect to every H5AD. It reuses the saved
V4.2 scVI latent-space UMAP sample and an existing unsupervised Leiden
partition, then summarizes NK- and T-lineage expression without assigning a
new training label or fitting a classifier.
"""

from __future__ import annotations

import html
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import anndata as ad  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

from workflows.gdtai.gdtai_v4_2_integration_core import (  # noqa: E402
    h5ad_obs_frame,
    read_json,
    sha256_file,
)


BASE_NAME = "gdtai_v4_2_nk_cluster_review"
TABLE_DIR = PROJECT_ROOT / f"Integrated_dataset/tables/gdT_prediction/{BASE_NAME}"
FIGURE_DIR = PROJECT_ROOT / f"Integrated_dataset/figures/gdT_prediction/{BASE_NAME}"
LOG_DIR = PROJECT_ROOT / f"Integrated_dataset/logs/gdT_prediction/{BASE_NAME}"
STATIC_DIR = PROJECT_ROOT / f"gdT_prediction/{BASE_NAME}"

SSD_DIR = Path("/ssd/tnk_phase3/Integrated_dataset/gdtai_v4_2_nk_reference")
STAGED_H5AD = SSD_DIR / "development_hvg_counts.h5ad"
PARTITIONS = SSD_DIR / "cluster_partitions.npz"
DIAGNOSTIC_TABLE = (
    PROJECT_ROOT
    / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_nk_reference_diagnostics/diagnostic_visualization_sample.csv.gz"
)
DIAGNOSTIC_COORDINATES = (
    PROJECT_ROOT
    / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_nk_reference_diagnostics/diagnostic_umap_coordinates.npz"
)
REPAIR_CANDIDATES = (
    PROJECT_ROOT
    / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_nk_definition_repair/candidate_nk_evidence.csv.gz"
)
V4_LABELS = (
    PROJECT_ROOT
    / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_preflight/cell_label_manifest.csv.gz"
)
CURRENT_ATLAS = PROJECT_ROOT / "Integrated_dataset/TNK_cleaned.h5ad"
FIT_SUMMARY = (
    PROJECT_ROOT
    / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_nk_reference/fit_summary.json"
)
CLUSTER_SUMMARY = (
    PROJECT_ROOT
    / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_nk_reference/cluster_summary.json"
)

RANDOM_SEED = 20260817
CLUSTER_RUN = "global_r0.8_s29"
EXPRESSION_CHUNK = 10_000

T_LINEAGE_GENES = ["CD3E", "TRAT1", "BCL11B", "LCK", "CD247", "THEMIS"]
NK_LINEAGE_GENES = [
    "FCER1G",
    "TYROBP",
    "KLRD1",
    "NCR1",
    "FCGR3A",
    "KLRC1",
    "KLRF1",
    "S1PR5",
]
CYTOTOXIC_GENES = ["NKG7", "GNLY", "PRF1", "GZMB", "XCL1", "XCL2"]
MYELOID_CONTEXT_GENES = ["LST1", "AIF1", "CTSS", "LILRB1", "LILRB2", "S100A8", "S100A9"]
MARKER_GENES = (
    T_LINEAGE_GENES + NK_LINEAGE_GENES + MYELOID_CONTEXT_GENES + CYTOTOXIC_GENES
)

ANNOTATION_ORDER = [
    "Unlabeled / other",
    "Productive TRA/TRB T",
    "Prior scANVI NK only",
    "Primary NK: scANVI + author",
    "New NK_CONFIDENT",
]
ANNOTATION_COLORS = {
    "Unlabeled / other": "#cbd2d6",
    "Productive TRA/TRB T": "#247ba0",
    "Prior scANVI NK only": "#d9a441",
    "Primary NK: scANVI + author": "#2a7f62",
    "New NK_CONFIDENT": "#7a3e9d",
}


def configure_plotting() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 13,
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def finish_figure(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=260, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def clean_umap_axis(ax: plt.Axes, title: str) -> None:
    ax.set_title(title)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)


def nonempty(values: pd.Series) -> np.ndarray:
    return (
        ~values.astype("string")
        .fillna("")
        .str.strip()
        .str.lower()
        .isin(["", "nan", "none", "na"])
        .to_numpy()
    )


def restore_productive_t(obs: pd.DataFrame) -> np.ndarray:
    primary = obs["primary_nk_anchor"].astype(bool).to_numpy()
    doublet = obs["doublet_flag_effective"].astype(bool).to_numpy()
    productive = obs["productive_tcr_anchor"].astype(bool).to_numpy().copy()
    current_rows = np.flatnonzero(
        obs["input_cohort_id"].astype(str).eq("current_atlas").to_numpy()
    )
    current = h5ad_obs_frame(CURRENT_ATLAS, ["TRA_cdr3", "TRB_cdr3"])
    current_productive = nonempty(current["TRA_cdr3"]) | nonempty(current["TRB_cdr3"])
    input_rows = obs.iloc[current_rows]["input_row"].to_numpy(np.int64)
    productive[current_rows] |= current_productive[input_rows]
    productive &= ~primary & ~doublet
    return productive


def load_prior_scanvi_nk(sample_obs: pd.DataFrame) -> np.ndarray:
    labels = pd.read_csv(
        V4_LABELS,
        usecols=[
            "cell_id",
            "expression_input_id",
            "stage1_role",
            "nk_annotation_strength",
        ],
    )
    labels = labels[
        labels["expression_input_id"].eq("current_atlas")
        & labels["stage1_role"].eq("nk_negative")
    ]
    if labels["cell_id"].duplicated().any():
        raise RuntimeError(
            "Prior scANVI NK manifest contains duplicate current-atlas cells"
        )
    strength = labels.set_index(labels["cell_id"].astype(str))["nk_annotation_strength"]
    mapped = strength.reindex(sample_obs["original_cell_id"].astype(str)).fillna(0.0)
    return mapped.to_numpy(np.float32)


def extract_marker_expression(
    sample_indices: np.ndarray,
    marker_genes: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    cache_path = TABLE_DIR / "diagnostic_sample_marker_expression.npz"
    if cache_path.exists():
        cached = np.load(cache_path, allow_pickle=False)
        if (
            np.array_equal(cached["sample_indices"], sample_indices)
            and cached["genes"].astype(str).tolist() == marker_genes
        ):
            return (
                np.asarray(
                    cached["log1p_cp10k_integration_features"], dtype=np.float32
                ),
                np.asarray(cached["detected"], dtype=bool),
            )

    backed = ad.read_h5ad(STAGED_H5AD, backed="r")
    missing = [gene for gene in marker_genes if gene not in backed.var_names]
    if missing:
        backed.file.close()
        raise KeyError(f"Staged scVI matrix lacks requested marker genes: {missing}")
    columns = np.asarray([backed.var_names.get_loc(gene) for gene in marker_genes])
    values = np.zeros((sample_indices.size, len(marker_genes)), dtype=np.float32)
    totals = np.zeros(sample_indices.size, dtype=np.float64)
    for start in range(0, sample_indices.size, EXPRESSION_CHUNK):
        end = min(sample_indices.size, start + EXPRESSION_CHUNK)
        rows = sample_indices[start:end]
        block = backed.X[rows, :].tocsr()
        totals[start:end] = np.asarray(block.sum(axis=1)).ravel()
        values[start:end] = block[:, columns].toarray().astype(np.float32, copy=False)
        if end % 50_000 == 0 or end == sample_indices.size:
            print(
                f"Marker extraction: {end:,}/{sample_indices.size:,} cells", flush=True
            )
    backed.file.close()
    if (totals <= 0).any():
        raise RuntimeError(
            "Diagnostic sample contains cells with zero staged-feature counts"
        )
    detected = values > 0
    normalized = np.log1p(values * (10_000.0 / totals[:, None])).astype(np.float32)
    np.savez_compressed(
        cache_path,
        sample_indices=sample_indices,
        genes=np.asarray(marker_genes, dtype="U"),
        log1p_cp10k_integration_features=normalized,
        detected=detected,
    )
    return normalized, detected


def robust_program_score(expression: np.ndarray, columns: list[int]) -> np.ndarray:
    panel = expression[:, columns].astype(np.float32, copy=True)
    for column in range(panel.shape[1]):
        positive = panel[:, column][panel[:, column] > 0]
        ceiling = float(np.quantile(positive, 0.99)) if positive.size else 1.0
        panel[:, column] = np.clip(panel[:, column] / max(ceiling, 1e-6), 0.0, 1.0)
    return panel.mean(axis=1, dtype=np.float64).astype(np.float32)


def sampling_weights(obs: pd.DataFrame, sample: pd.DataFrame) -> np.ndarray:
    primary = obs["primary_nk_anchor"].astype(bool).to_numpy()
    productive = obs["productive_tcr_anchor"].astype(bool).to_numpy() & ~primary
    candidate = (
        obs["candidate_eligible"].astype(bool).to_numpy() & ~primary & ~productive
    )
    role = np.full(obs.shape[0], "Other", dtype=object)
    role[candidate] = "Eligible candidate"
    role[productive] = "Productive TCR anchor"
    role[primary] = "Primary NK anchor"
    full = pd.DataFrame(
        {
            "source_gse_id": obs["source_gse_id"].astype(str).to_numpy(),
            "diagnostic_role": role,
        }
    )
    full_count = full.value_counts(sort=False)
    sample_count = sample[["source_gse_id", "diagnostic_role"]].value_counts(sort=False)
    weights = np.ones(sample.shape[0], dtype=np.float64)
    for key, n_sample in sample_count.items():
        mask = sample["source_gse_id"].eq(key[0]) & sample["diagnostic_role"].eq(key[1])
        weights[mask.to_numpy()] = float(full_count.get(key, n_sample) / n_sample)
    return weights


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    return float(np.average(values, weights=weights)) if values.size else float("nan")


def cluster_marker_summary(
    sample: pd.DataFrame,
    expression: np.ndarray,
    detected: np.ndarray,
    weights: np.ndarray,
    t_score: np.ndarray,
    nk_score: np.ndarray,
    myeloid_score: np.ndarray,
    cytotoxic_score: np.ndarray,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for cluster in sorted(sample["cluster"].unique(), key=int):
        mask = sample["cluster"].eq(cluster).to_numpy()
        local_weights = weights[mask]
        row: dict[str, Any] = {
            "cluster": str(cluster),
            "n_sample": int(mask.sum()),
            "weighted_t_program": weighted_mean(t_score[mask], local_weights),
            "weighted_nk_program": weighted_mean(nk_score[mask], local_weights),
            "weighted_myeloid_program": weighted_mean(
                myeloid_score[mask], local_weights
            ),
            "weighted_cytotoxic_program": weighted_mean(
                cytotoxic_score[mask], local_weights
            ),
        }
        row["nk_minus_t_program"] = (
            row["weighted_nk_program"] - row["weighted_t_program"]
        )
        row["nk_specific_margin"] = row["weighted_nk_program"] - max(
            row["weighted_t_program"], row["weighted_myeloid_program"]
        )
        for index, gene in enumerate(MARKER_GENES):
            row[f"{gene}_detected_fraction"] = weighted_mean(
                detected[mask, index].astype(np.float32), local_weights
            )
            row[f"{gene}_mean"] = weighted_mean(expression[mask, index], local_weights)
        rows.append(row)
    result = pd.DataFrame(rows)
    result["nk_likeness_rank"] = (
        result["nk_specific_margin"].rank(method="min", ascending=False).astype(int)
    )
    return result


def full_cluster_summary(
    obs: pd.DataFrame,
    labels: np.ndarray,
    productive_t: np.ndarray,
    confident_rows: np.ndarray,
) -> pd.DataFrame:
    primary = obs["primary_nk_anchor"].astype(bool).to_numpy()
    candidate = obs["candidate_eligible"].astype(bool).to_numpy()
    confident = np.zeros(obs.shape[0], dtype=bool)
    confident[confident_rows] = True
    source = obs["source_gse_id"].astype(str).to_numpy()
    rows: list[dict[str, Any]] = []
    for cluster in range(int(labels.max()) + 1):
        mask = labels == cluster
        candidate_mask = mask & candidate
        values, counts = np.unique(source[candidate_mask], return_counts=True)
        if counts.size:
            dominant_index = int(np.argmax(counts))
            dominant_source = str(values[dominant_index])
            dominant_fraction = float(counts[dominant_index] / counts.sum())
        else:
            dominant_source = ""
            dominant_fraction = float("nan")
        rows.append(
            {
                "cluster": str(cluster),
                "n_cells": int(mask.sum()),
                "n_primary_nk": int((mask & primary).sum()),
                "primary_nk_fraction": float(primary[mask].mean()),
                "n_productive_t": int((mask & productive_t).sum()),
                "productive_t_fraction": float(productive_t[mask].mean()),
                "n_eligible_candidates": int(candidate_mask.sum()),
                "n_confident_nk": int((mask & confident).sum()),
                "candidate_source_count": int(counts.size),
                "dominant_candidate_source": dominant_source,
                "dominant_candidate_source_fraction": dominant_fraction,
            }
        )
    return pd.DataFrame(rows)


def make_feature_panel(
    coordinates: np.ndarray,
    expression: np.ndarray,
    genes: list[str],
    title: str,
    path: Path,
    columns: int,
) -> None:
    lookup = {gene: index for index, gene in enumerate(MARKER_GENES)}
    rows = int(np.ceil(len(genes) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(4.0 * columns, 3.35 * rows))
    axes = np.asarray(axes).reshape(-1)
    for ax, gene in zip(axes, genes, strict=False):
        values = expression[:, lookup[gene]]
        positive = values[values > 0]
        ceiling = float(np.quantile(positive, 0.995)) if positive.size else 1.0
        order = np.argsort(values, kind="mergesort")
        plotted = ax.scatter(
            coordinates[order, 0],
            coordinates[order, 1],
            c=np.clip(values[order], 0, ceiling),
            cmap="viridis",
            vmin=0,
            vmax=ceiling,
            s=0.6,
            alpha=0.75,
            linewidths=0,
            rasterized=True,
        )
        clean_umap_axis(ax, gene)
        colorbar = fig.colorbar(plotted, ax=ax, shrink=0.65, pad=0.01)
        colorbar.ax.tick_params(labelsize=7)
    for ax in axes[len(genes) :]:
        ax.axis("off")
    fig.suptitle(title, fontsize=16, y=1.005)
    finish_figure(fig, path)


def make_program_figure(
    coordinates: np.ndarray,
    t_score: np.ndarray,
    nk_score: np.ndarray,
    myeloid_score: np.ndarray,
    cytotoxic_score: np.ndarray,
    path: Path,
) -> None:
    values = [
        t_score,
        nk_score,
        myeloid_score,
        cytotoxic_score,
        nk_score - t_score,
        nk_score - np.maximum(t_score, myeloid_score),
    ]
    titles = [
        "T-lineage program",
        "NK-lineage program",
        "Myeloid contamination context",
        "Shared cytotoxic program",
        "NK minus T lineage balance",
        "NK-specific margin",
    ]
    fig, axes = plt.subplots(2, 3, figsize=(14.2, 8.0))
    for ax, local, title in zip(axes.flat, values, titles, strict=True):
        order = np.argsort(local, kind="mergesort")
        if "balance" in title or "margin" in title:
            limit = float(np.quantile(np.abs(local), 0.995))
            norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
            plotted = ax.scatter(
                coordinates[order, 0],
                coordinates[order, 1],
                c=local[order],
                cmap="coolwarm",
                norm=norm,
                s=0.7,
                alpha=0.75,
                linewidths=0,
                rasterized=True,
            )
        else:
            ceiling = float(np.quantile(local, 0.995))
            plotted = ax.scatter(
                coordinates[order, 0],
                coordinates[order, 1],
                c=np.clip(local[order], 0, ceiling),
                cmap="viridis",
                vmin=0,
                vmax=ceiling,
                s=0.7,
                alpha=0.75,
                linewidths=0,
                rasterized=True,
            )
        clean_umap_axis(ax, title)
        fig.colorbar(plotted, ax=ax, shrink=0.68, pad=0.01)
    finish_figure(fig, path)


def make_annotation_umap(
    sample: pd.DataFrame, coordinates: np.ndarray, path: Path
) -> None:
    fig, ax = plt.subplots(figsize=(9.0, 7.0))
    values = sample["annotation_evidence"].astype(str).to_numpy()
    for label in ANNOTATION_ORDER:
        mask = values == label
        if not mask.any():
            continue
        background = label == "Unlabeled / other"
        ax.scatter(
            coordinates[mask, 0],
            coordinates[mask, 1],
            c=ANNOTATION_COLORS[label],
            s=0.45 if background else 1.2,
            alpha=0.10 if background else 0.68,
            linewidths=0,
            rasterized=True,
            label=f"{label} ({mask.sum():,})",
        )
    clean_umap_axis(ax, "Recoverable prior annotation evidence on the scVI UMAP")
    ax.legend(frameon=False, fontsize=8, markerscale=5, loc="best")
    finish_figure(fig, path)


def cluster_colors(n_clusters: int) -> list[Any]:
    maps = [plt.get_cmap("tab20"), plt.get_cmap("tab20b")]
    return [maps[index // 20](index % 20) for index in range(n_clusters)]


def make_cluster_umap(
    sample: pd.DataFrame, coordinates: np.ndarray, path: Path
) -> None:
    clusters = sample["cluster"].astype(str).to_numpy()
    levels = sorted(np.unique(clusters), key=int)
    colors = cluster_colors(len(levels))
    fig, ax = plt.subplots(figsize=(10.0, 7.6))
    for color, level in zip(colors, levels, strict=True):
        mask = clusters == level
        ax.scatter(
            coordinates[mask, 0],
            coordinates[mask, 1],
            c=[color],
            s=0.7,
            alpha=0.55,
            linewidths=0,
            rasterized=True,
        )
        center = np.median(coordinates[mask], axis=0)
        ax.text(
            center[0],
            center[1],
            level,
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold",
            bbox={
                "boxstyle": "round,pad=0.16",
                "fc": "white",
                "ec": "#6d7a83",
                "alpha": 0.8,
            },
        )
    clean_umap_axis(ax, f"Unsupervised scVI Leiden clusters: {CLUSTER_RUN}")
    finish_figure(fig, path)


def make_cluster_heatmap(summary: pd.DataFrame, path: Path) -> None:
    ordered = summary.sort_values("nk_likeness_rank").copy()
    genes = T_LINEAGE_GENES + NK_LINEAGE_GENES + MYELOID_CONTEXT_GENES + CYTOTOXIC_GENES
    matrix = np.asarray(
        [
            [row[f"{gene}_detected_fraction"] for gene in genes]
            for _, row in ordered.iterrows()
        ]
    )
    fig, ax = plt.subplots(figsize=(14.0, 8.4))
    image = ax.imshow(matrix, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_yticks(np.arange(ordered.shape[0]), [f"c{x}" for x in ordered["cluster"]])
    ax.set_xticks(np.arange(len(genes)), genes, rotation=45, ha="right")
    ax.axvline(len(T_LINEAGE_GENES) - 0.5, color="white", linewidth=2)
    ax.axvline(
        len(T_LINEAGE_GENES) + len(NK_LINEAGE_GENES) - 0.5, color="white", linewidth=2
    )
    ax.axvline(
        len(T_LINEAGE_GENES) + len(NK_LINEAGE_GENES) + len(MYELOID_CONTEXT_GENES) - 0.5,
        color="white",
        linewidth=2,
    )
    ax.set_title("Sampling-weighted marker detection by unsupervised cluster")
    colorbar = fig.colorbar(image, ax=ax, shrink=0.82, pad=0.01)
    colorbar.set_label("Estimated cells with detected expression")
    finish_figure(fig, path)


def make_cluster_score_plot(summary: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 6.6))
    size = 30 + 280 * np.sqrt(summary["primary_nk_fraction"].to_numpy())
    plotted = ax.scatter(
        summary["weighted_t_program"],
        summary["weighted_nk_program"],
        c=summary["weighted_myeloid_program"],
        s=size,
        cmap="magma",
        alpha=0.82,
        edgecolors="#24333d",
        linewidths=0.4,
    )
    for row in summary.itertuples(index=False):
        ax.text(
            row.weighted_t_program,
            row.weighted_nk_program,
            str(row.cluster),
            fontsize=7,
            ha="center",
            va="center",
            color="white"
            if row.weighted_myeloid_program
            > summary["weighted_myeloid_program"].median()
            else "#17212b",
        )
    ax.set_xlabel("T-lineage program")
    ax.set_ylabel("NK-lineage program")
    ax.set_title("Cluster lineage balance; point size reflects primary-NK fraction")
    colorbar = fig.colorbar(plotted, ax=ax, pad=0.02)
    colorbar.set_label("Myeloid contamination context")
    ax.grid(color="#e2e8eb", linewidth=0.6)
    finish_figure(fig, path)


def table_html(frame: pd.DataFrame) -> str:
    return frame.to_html(index=False, border=0, classes="data-table", escape=True)


def render_report(
    summary: dict[str, Any],
    cluster_table: pd.DataFrame,
    tier_summary: pd.DataFrame,
) -> Path:
    prefix = f"../../Integrated_dataset/figures/gdT_prediction/{BASE_NAME}"
    top = cluster_table.sort_values("nk_likeness_rank").head(10).copy()
    display = top[
        [
            "cluster",
            "nk_likeness_rank",
            "n_cells",
            "n_primary_nk",
            "n_productive_t",
            "n_eligible_candidates",
            "n_confident_nk",
            "weighted_t_program",
            "weighted_nk_program",
            "weighted_myeloid_program",
            "weighted_cytotoxic_program",
            "nk_minus_t_program",
            "nk_specific_margin",
            "dominant_candidate_source",
            "dominant_candidate_source_fraction",
        ]
    ].copy()
    for column in [
        "weighted_t_program",
        "weighted_nk_program",
        "weighted_myeloid_program",
        "weighted_cytotoxic_program",
        "nk_minus_t_program",
        "nk_specific_margin",
    ]:
        display[column] = display[column].map(lambda value: f"{value:.3f}")
    display["dominant_candidate_source_fraction"] = display[
        "dominant_candidate_source_fraction"
    ].map(lambda value: "" if pd.isna(value) else f"{value:.1%}")
    tier_display = tier_summary.copy()
    tier_display["maximum_source_fraction"] = tier_display[
        "maximum_source_fraction"
    ].map(lambda value: f"{value:.1%}")

    document = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>gdTAI V4.2 NK and T lineage cluster review</title><style>
:root{{--ink:#17212b;--navy:#17324d;--teal:#197278;--amber:#b86e16;--line:#ccd6dc;--soft:#f3f6f7}}
*{{box-sizing:border-box}}body{{font-family:Arial,Helvetica,sans-serif;color:var(--ink);line-height:1.48;max-width:1160px;margin:0 auto;padding:30px 42px;background:white}}h1{{font-size:29px;color:var(--navy);border-bottom:4px solid var(--teal);padding-bottom:12px}}h2{{font-size:20px;color:var(--navy);margin-top:31px}}p,li{{font-size:14px}}.status{{border-left:6px solid var(--teal);background:#eef8f6;padding:14px 18px;margin:18px 0}}.warning{{border-left:6px solid var(--amber);background:#fff9ef;padding:13px 17px;margin:18px 0}}.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin:18px 0}}.card{{border:1px solid var(--line);padding:12px;background:var(--soft)}}.card b{{display:block;font-size:21px;color:var(--navy)}}.card span{{font-size:11px;color:#52636e}}.grid{{display:grid;grid-template-columns:1fr 1fr;gap:18px}}figure{{margin:18px 0;break-inside:avoid}}img{{display:block;width:100%;height:auto;margin:0 auto}}figcaption{{font-size:11px;color:#52636e;text-align:center;margin-top:7px}}.table-wrap{{overflow-x:auto;margin:14px 0 24px}}table{{border-collapse:collapse;width:100%;font-size:9px}}th{{background:var(--navy);color:white;text-align:left}}th,td{{padding:5px 6px;border:1px solid var(--line);vertical-align:top;overflow-wrap:anywhere}}tr:nth-child(even){{background:var(--soft)}}code{{background:#edf1f4;padding:1px 4px}}
@media(max-width:800px){{body{{padding:18px}}.cards,.grid{{grid-template-columns:1fr}}}}
@media print{{@page{{size:A4 landscape;margin:8mm}}body{{max-width:none;padding:0}}h1{{font-size:22px}}h2{{font-size:16px;break-after:avoid}}p,li{{font-size:9.5px}}.cards{{grid-template-columns:repeat(4,1fr)}}.card b{{font-size:15px}}.card span{{font-size:7px}}figure{{break-inside:avoid}}figure img{{max-height:165mm;object-fit:contain}}table{{font-size:6.4px}}tr{{break-inside:avoid}}.table-wrap{{overflow:visible}}}}
</style></head><body>
<h1>gdTAI V4.2 NK and T lineage cluster review</h1>
<p>This read-only review uses the completed development-only scVI representation. UMAP coordinates and marker plots use the deterministic {summary["n_sample_cells"]:,}-cell diagnostic sample; cluster cell and anchor counts use all {summary["n_sidecar_cells"]:,} cells. No locked cohort was included.</p>
<div class="status"><b>Result: PASS_VISUAL_REVIEW_READY.</b> The feature plots, recoverable annotation evidence, and independent Leiden clusters are ready for biological review before changing the NK label definition.</div>
<div class="cards"><div class="card"><b>{summary["n_clusters"]}</b><span>unsupervised clusters</span></div><div class="card"><b>{summary["n_primary_nk"]:,}</b><span>primary NK anchors</span></div><div class="card"><b>{summary["n_productive_t"]:,}</b><span>restored productive-T anchors</span></div><div class="card"><b>{summary["n_confident_nk"]:,}</b><span>current NK_CONFIDENT cells</span></div></div>

<h2>1. Annotation evidence on the scVI UMAP</h2><p>The missing historical integrated H5AD prevents recovery of its complete full-cell scANVI annotation column. The plot therefore shows only traceable evidence: prior scANVI-NK cells retained in the frozen V4 manifest, dual scANVI-plus-author NK anchors, restored productive-TRA/TRB T anchors, and the 469 current two-evidence NK candidates.</p><figure><img src="{prefix}/umap_annotation_evidence.png" alt="Annotation evidence UMAP"><figcaption>These are provenance classes, not a new annotation model.</figcaption></figure>

<h2>2. Unsupervised structure</h2><div class="grid"><figure><img src="{prefix}/umap_unsupervised_clusters.png" alt="Unsupervised cluster UMAP"><figcaption>Frozen global Leiden partition <code>{CLUSTER_RUN}</code>.</figcaption></figure><figure><img src="{prefix}/cluster_lineage_score_scatter.png" alt="Cluster lineage score"><figcaption>Scores combine multiple genes; color shows myeloid contamination context and point size shows primary-NK fraction.</figcaption></figure></div>

<h2>3. T-cell lineage feature plots</h2><figure><img src="{prefix}/feature_umap_t_lineage.png" alt="T lineage markers"><figcaption>Expression is log1p CP10K within the fixed 4,000 scVI integration features.</figcaption></figure>
<h2>4. NK lineage feature plots</h2><figure><img src="{prefix}/feature_umap_nk_lineage.png" alt="NK lineage markers"><figcaption>NK evidence requires a multigene program; no single receptor is treated as definitive.</figcaption></figure>
<h2>5. Non-T/NK contamination context</h2><p><code>FCER1G</code> and <code>TYROBP</code> are not sufficient for NK identity because myeloid cells also express them. The following panel is used only to recognize obvious non-T/NK clusters; these genes are not proposed as gdTAI classifier features.</p><figure><img src="{prefix}/feature_umap_myeloid_context.png" alt="Myeloid context markers"><figcaption>Clusters with broad LST1/AIF1/CTSS/LILRB/S100 support should not be called NK solely from adaptor expression.</figcaption></figure>
<h2>6. Shared cytotoxic program</h2><div class="warning"><b>Do not use cytotoxicity as NK truth.</b> NKG7, GNLY, PRF1, GZMB, XCL1, and XCL2 can also be high in CD8 and gamma-delta T cells.</div><figure><img src="{prefix}/feature_umap_shared_cytotoxic.png" alt="Cytotoxic markers"><figcaption>These markers are shown as phenotype context only.</figcaption></figure>

<h2>7. Multigene programs and cluster heatmap</h2><figure><img src="{prefix}/umap_lineage_programs.png" alt="Lineage programs"><figcaption>Each program is the mean of gene-wise robust-scaled expression. NK-specific margin subtracts the stronger of T-lineage or myeloid-context evidence.</figcaption></figure><figure><img src="{prefix}/cluster_marker_detection_heatmap.png" alt="Cluster marker heatmap"><figcaption>Sampling-weighted marker detection estimates compensate for the diagnostic sample's source and evidence-role oversampling.</figcaption></figure>

<h2>8. Clusters ranked for review</h2><p>The table ranks clusters by the continuous NK-specific margin: NK lineage minus the stronger of T-lineage or myeloid-context evidence. It does not impose a new hard NK cutoff. Exact anchor, candidate, and source-composition counts are included so expression, annotation provenance, and source bias can be considered together.</p><div class="table-wrap">{table_html(display)}</div>

<h2>9. Recommended decision</h2><div class="status"><b>Strict core recommendation: retain {summary["n_recommended_core_nk"]:,} cells.</b> These cells pass both independent cell-level evidences and fall in cluster 19, which contains {summary["core_cluster_primary_nk_fraction_of_all"]:.1%} of all primary NK anchors. Hold the remaining {summary["n_current_confident_outside_core"]:,} current calls outside cluster 19 rather than using them as strict training NK labels.</div><p>No automatic expansion is recommended yet. Cluster 19 contains plausible dropout-rescue tiers, but its eligible candidates are {summary["core_cluster_dominant_candidate_source_fraction"]:.1%} from one source, above the 70% source cap. Cluster 9 is a biologically plausible secondary NK-like cluster with better source balance, but it has only {summary["secondary_cluster_primary_nk"]:,} primary NK anchors and weak latent-anchor passage. Both require source-resolved review before use.</p><div class="table-wrap">{table_html(tier_display)}</div>

<h2>Decision boundary</h2><p><b>No H5AD or frozen label manifest is changed by this report.</b> A defensible expansion should require coherent NK-lineage expression, low T-lineage and myeloid support, localization to reproducible unsupervised clusters, and acceptable source composition. Cytotoxicity alone, one receptor, or the historical scANVI label alone is insufficient.</p>
<h2>Reproducibility</h2><ul><li>Latent SHA-256: <code>{html.escape(summary["latent_sha256"])}</code></li><li>Partition SHA-256: <code>{html.escape(summary["partitions_sha256"])}</code></li><li>Cluster run: <code>{CLUSTER_RUN}</code></li><li>Marker-expression cache SHA-256: <code>{html.escape(summary["marker_cache_sha256"])}</code></li><li>Script SHA-256: <code>{html.escape(summary["script_sha256"])}</code></li><li>H5AD mutation: <b>none</b>; classifier fitting: <b>none</b>; GitHub push: <b>none</b>.</li></ul>
</body></html>"""
    path = STATIC_DIR / "index.html"
    path.write_text(document, encoding="utf-8")
    return path


def export_pdf(html_path: Path) -> Path:
    chrome = shutil.which("google-chrome") or shutil.which("google-chrome-stable")
    if chrome is None:
        raise RuntimeError("Google Chrome is required for PDF export")
    output = STATIC_DIR / f"{BASE_NAME}_report.pdf"
    subprocess.run(
        [
            chrome,
            "--headless",
            "--no-sandbox",
            "--disable-gpu",
            "--print-to-pdf-no-header",
            f"--print-to-pdf={output}",
            html_path.resolve().as_uri(),
        ],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if not output.exists() or output.stat().st_size < 10_000:
        raise RuntimeError("PDF export did not produce a valid report")
    return output


def main() -> None:
    configure_plotting()
    for directory in [TABLE_DIR, FIGURE_DIR, LOG_DIR, STATIC_DIR]:
        directory.mkdir(parents=True, exist_ok=True)

    fit_summary = read_json(FIT_SUMMARY)
    cluster_summary = read_json(CLUSTER_SUMMARY)
    if sha256_file(SSD_DIR / "X_scVI.npy") != fit_summary["latent_sha256"]:
        raise RuntimeError("Saved scVI latent checksum differs from fit summary")
    if sha256_file(PARTITIONS) != cluster_summary["partitions_sha256"]:
        raise RuntimeError(
            "Saved cluster partition checksum differs from cluster summary"
        )

    diagnostic = pd.read_csv(DIAGNOSTIC_TABLE)
    coordinates_saved = np.load(DIAGNOSTIC_COORDINATES, allow_pickle=False)
    sample_indices = coordinates_saved["sample_indices"].astype(np.int64)
    coordinates = coordinates_saved["X_umap"].astype(np.float32)
    if diagnostic.shape[0] != sample_indices.size or coordinates.shape != (
        sample_indices.size,
        2,
    ):
        raise RuntimeError("Diagnostic table, sample index, and UMAP sizes differ")

    backed = ad.read_h5ad(STAGED_H5AD, backed="r")
    obs = backed.obs.copy()
    sample_obs = backed.obs.iloc[sample_indices].copy()
    backed.file.close()
    if not np.array_equal(
        diagnostic["integration_cell_id"].astype(str).to_numpy(),
        sample_obs.index.astype(str).to_numpy(),
    ):
        raise RuntimeError("Diagnostic sample IDs do not match the staged sidecar rows")
    if (
        obs["input_cohort_id"]
        .astype(str)
        .isin(["GSE169246", "GSE315928", "GSE121636_GSE121637"])
        .any()
    ):
        raise RuntimeError("A locked cohort leaked into the sidecar review")

    saved = np.load(PARTITIONS, allow_pickle=False)
    names = saved["global_names"].astype(str).tolist()
    if CLUSTER_RUN not in names:
        raise KeyError(f"Requested cluster run is absent: {CLUSTER_RUN}")
    full_labels = saved["global_labels"][:, names.index(CLUSTER_RUN)].astype(np.int32)
    sample_labels = full_labels[sample_indices]
    diagnostic["cluster"] = sample_labels.astype(str)

    productive_t = restore_productive_t(obs)
    candidates = pd.read_csv(REPAIR_CANDIDATES)
    confident_rows = candidates.loc[
        candidates["nk_confident"].astype(bool), "global_row"
    ].to_numpy(np.int64)
    if np.unique(confident_rows).size != 469:
        raise RuntimeError("The repaired NK_CONFIDENT row count differs from 469")

    prior_scanvi_strength = load_prior_scanvi_nk(sample_obs)
    sample_primary = sample_obs["primary_nk_anchor"].astype(bool).to_numpy()
    sample_productive = productive_t[sample_indices]
    sample_confident = np.isin(sample_indices, confident_rows)
    annotation = np.full(sample_indices.size, "Unlabeled / other", dtype=object)
    annotation[sample_productive] = "Productive TRA/TRB T"
    annotation[prior_scanvi_strength == 0.5] = "Prior scANVI NK only"
    annotation[sample_primary] = "Primary NK: scANVI + author"
    annotation[sample_confident] = "New NK_CONFIDENT"
    diagnostic["annotation_evidence"] = annotation

    expression, detected = extract_marker_expression(sample_indices, MARKER_GENES)
    lookup = {gene: index for index, gene in enumerate(MARKER_GENES)}
    t_score = robust_program_score(
        expression, [lookup[gene] for gene in T_LINEAGE_GENES]
    )
    nk_score = robust_program_score(
        expression, [lookup[gene] for gene in NK_LINEAGE_GENES]
    )
    myeloid_score = robust_program_score(
        expression, [lookup[gene] for gene in MYELOID_CONTEXT_GENES]
    )
    cytotoxic_score = robust_program_score(
        expression, [lookup[gene] for gene in CYTOTOXIC_GENES]
    )
    weights = sampling_weights(obs, diagnostic)

    marker_summary = cluster_marker_summary(
        diagnostic,
        expression,
        detected,
        weights,
        t_score,
        nk_score,
        myeloid_score,
        cytotoxic_score,
    )
    exact_summary = full_cluster_summary(obs, full_labels, productive_t, confident_rows)
    cluster_table = exact_summary.merge(
        marker_summary, on="cluster", how="left", validate="one_to_one"
    )
    candidate_cluster = candidates.copy()
    candidate_cluster["cluster"] = full_labels[
        candidate_cluster["global_row"].to_numpy(np.int64)
    ].astype(str)
    candidate_cluster_summary = (
        candidate_cluster.groupby("cluster", observed=True)
        .agg(
            n_candidates=("integration_cell_id", "size"),
            n_strict_expression=("strict_expression_nk", "sum"),
            n_latent=("latent_pass", "sum"),
            n_confident=("nk_confident", "sum"),
            median_latent_nk_fraction=("latent_nk_fraction", "median"),
            median_nearest_nk_distance=("nearest_nk_distance", "median"),
        )
        .reset_index()
    )
    cluster_value = candidate_cluster["cluster"].astype(str)
    expression_pass = candidate_cluster["strict_expression_nk"].astype(bool)
    latent_pass = candidate_cluster["latent_pass"].astype(bool)
    confident_pass = candidate_cluster["nk_confident"].astype(bool)
    candidate_cluster["review_tier"] = np.select(
        [
            cluster_value.eq("19") & confident_pass,
            cluster_value.eq("19") & latent_pass & ~expression_pass,
            cluster_value.eq("19") & expression_pass & ~latent_pass,
            cluster_value.eq("9") & (expression_pass | latent_pass),
            confident_pass,
        ],
        [
            "A_core_both_cluster19",
            "B_cluster19_latent_only_dropout_review",
            "C_cluster19_expression_only_distance_review",
            "D_cluster9_secondary_nk_review",
            "E_current_confident_outside_core_hold",
        ],
        default="F_no_nk_assignment",
    )
    tier_rows: list[dict[str, Any]] = []
    for tier, frame in candidate_cluster.groupby("review_tier", observed=True):
        source_counts = frame["source_gse_id"].astype(str).value_counts()
        tier_rows.append(
            {
                "review_tier": tier,
                "n_cells": int(frame.shape[0]),
                "n_sources": int(source_counts.size),
                "dominant_source": str(source_counts.index[0]),
                "maximum_source_fraction": float(
                    source_counts.iloc[0] / frame.shape[0]
                ),
                "recommended_use": {
                    "A_core_both_cluster19": "strict low-weight NK training reference",
                    "B_cluster19_latent_only_dropout_review": "dropout-rescue review only",
                    "C_cluster19_expression_only_distance_review": "source/distance review only",
                    "D_cluster9_secondary_nk_review": "secondary-cluster review only",
                    "E_current_confident_outside_core_hold": "hold out of strict NK training",
                    "F_no_nk_assignment": "do not assign NK",
                }[tier],
            }
        )
    tier_summary = pd.DataFrame(tier_rows).sort_values("review_tier")

    make_feature_panel(
        coordinates,
        expression,
        T_LINEAGE_GENES,
        "T-cell lineage markers on the scVI UMAP",
        FIGURE_DIR / "feature_umap_t_lineage.png",
        columns=3,
    )
    make_feature_panel(
        coordinates,
        expression,
        NK_LINEAGE_GENES,
        "NK-lineage markers on the scVI UMAP",
        FIGURE_DIR / "feature_umap_nk_lineage.png",
        columns=4,
    )
    make_feature_panel(
        coordinates,
        expression,
        MYELOID_CONTEXT_GENES,
        "Myeloid contamination context on the scVI UMAP",
        FIGURE_DIR / "feature_umap_myeloid_context.png",
        columns=4,
    )
    make_feature_panel(
        coordinates,
        expression,
        CYTOTOXIC_GENES,
        "Shared cytotoxic markers on the scVI UMAP",
        FIGURE_DIR / "feature_umap_shared_cytotoxic.png",
        columns=3,
    )
    make_program_figure(
        coordinates,
        t_score,
        nk_score,
        myeloid_score,
        cytotoxic_score,
        FIGURE_DIR / "umap_lineage_programs.png",
    )
    make_annotation_umap(
        diagnostic, coordinates, FIGURE_DIR / "umap_annotation_evidence.png"
    )
    make_cluster_umap(
        diagnostic, coordinates, FIGURE_DIR / "umap_unsupervised_clusters.png"
    )
    make_cluster_heatmap(
        cluster_table, FIGURE_DIR / "cluster_marker_detection_heatmap.png"
    )
    make_cluster_score_plot(
        cluster_table, FIGURE_DIR / "cluster_lineage_score_scatter.png"
    )

    diagnostic_output = diagnostic[
        [
            "integration_cell_id",
            "source_gse_id",
            "diagnostic_role",
            "cluster",
            "annotation_evidence",
            "UMAP1",
            "UMAP2",
        ]
    ].copy()
    diagnostic_output["sampling_weight"] = weights
    diagnostic_output["t_lineage_program"] = t_score
    diagnostic_output["nk_lineage_program"] = nk_score
    diagnostic_output["myeloid_context_program"] = myeloid_score
    diagnostic_output["cytotoxic_program"] = cytotoxic_score
    diagnostic_output.to_csv(
        TABLE_DIR / "diagnostic_sample_cluster_annotation_scores.csv.gz",
        index=False,
        compression="gzip",
    )
    cluster_table.to_csv(TABLE_DIR / "cluster_lineage_evidence.csv", index=False)
    candidate_cluster_summary.to_csv(
        TABLE_DIR / "candidate_evidence_by_cluster.csv", index=False
    )
    candidate_cluster.to_csv(
        TABLE_DIR / "candidate_cluster_review_tiers.csv.gz",
        index=False,
        compression="gzip",
    )
    tier_summary.to_csv(TABLE_DIR / "candidate_review_tier_summary.csv", index=False)
    pd.DataFrame(
        {
            "gene": MARKER_GENES,
            "program": ["T_lineage"] * len(T_LINEAGE_GENES)
            + ["NK_lineage"] * len(NK_LINEAGE_GENES)
            + ["myeloid_context_only"] * len(MYELOID_CONTEXT_GENES)
            + ["shared_cytotoxic_audit_only"] * len(CYTOTOXIC_GENES),
        }
    ).to_csv(TABLE_DIR / "marker_panel.csv", index=False)

    core_cluster = cluster_table.loc[cluster_table["cluster"].eq("19")].iloc[0]
    secondary_cluster = cluster_table.loc[cluster_table["cluster"].eq("9")].iloc[0]
    n_recommended_core = int(
        tier_summary.loc[
            tier_summary["review_tier"].eq("A_core_both_cluster19"), "n_cells"
        ].iloc[0]
    )
    n_confident_outside = int(confident_rows.size - n_recommended_core)
    summary = {
        "result": "PASS_VISUAL_REVIEW_READY",
        "n_sidecar_cells": int(obs.shape[0]),
        "n_sample_cells": int(sample_indices.size),
        "n_clusters": int(np.unique(full_labels).size),
        "cluster_run": CLUSTER_RUN,
        "n_primary_nk": int(obs["primary_nk_anchor"].astype(bool).sum()),
        "n_productive_t": int(productive_t.sum()),
        "n_confident_nk": int(confident_rows.size),
        "recommended_core_cluster": "19",
        "n_recommended_core_nk": n_recommended_core,
        "n_current_confident_outside_core": n_confident_outside,
        "core_cluster_primary_nk_fraction_of_all": float(
            core_cluster["n_primary_nk"] / obs["primary_nk_anchor"].astype(bool).sum()
        ),
        "core_cluster_dominant_candidate_source_fraction": float(
            core_cluster["dominant_candidate_source_fraction"]
        ),
        "secondary_cluster": "9",
        "secondary_cluster_primary_nk": int(secondary_cluster["n_primary_nk"]),
        "n_sample_prior_scanvi_only": int((prior_scanvi_strength == 0.5).sum()),
        "n_sample_primary_nk": int(sample_primary.sum()),
        "n_sample_productive_t": int(sample_productive.sum()),
        "n_sample_confident_nk": int(sample_confident.sum()),
        "top_nk_like_clusters": cluster_table.sort_values("nk_likeness_rank")
        .head(6)["cluster"]
        .astype(str)
        .tolist(),
        "latent_sha256": fit_summary["latent_sha256"],
        "partitions_sha256": cluster_summary["partitions_sha256"],
        "marker_cache_sha256": sha256_file(
            TABLE_DIR / "diagnostic_sample_marker_expression.npz"
        ),
        "script_sha256": sha256_file(Path(__file__)),
        "locked_cohorts_included": False,
        "source_h5ad_mutation_performed": False,
        "classifier_fitting_performed": False,
        "nk_labels_changed": False,
        "github_push_performed": False,
        "expression_scale": "log1p_cp10k_within_4000_scvi_integration_features",
        "sampling_weighted_cluster_expression": True,
        "historical_full_scanvi_annotation_available": False,
    }
    (LOG_DIR / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    html_path = render_report(summary, cluster_table, tier_summary)
    pdf_path = export_pdf(html_path)
    print(
        json.dumps(
            {
                "result": summary["result"],
                "html": str(html_path),
                "pdf": str(pdf_path),
                "top_nk_like_clusters": summary["top_nk_like_clusters"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
