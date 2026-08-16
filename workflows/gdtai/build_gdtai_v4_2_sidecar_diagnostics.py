#!/usr/bin/env python3
"""Visualize the V4.2 development-only scVI sidecar and NK-label failure mode."""

from __future__ import annotations

import html
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

_inherited_mps = os.environ.get("CUDA_MPS_PIPE_DIRECTORY", "")
if _inherited_mps == "/tmp/nvidia-mps" or (
    _inherited_mps and not _inherited_mps.startswith("/tmp/gdtai-v42-diagnostics-mps-")
):
    os.environ.pop("CUDA_MPS_PIPE_DIRECTORY", None)
os.environ.setdefault("CUDA_MPS_PIPE_DIRECTORY", f"/tmp/gdtai-v42-diagnostics-mps-{os.getuid()}-{os.getpid()}")
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "0")
os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
os.environ.setdefault("PYTHONHASHSEED", "20260816")

import torch  # Load torch's CUDA runtime before RAPIDS libraries.
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rapids_singlecell as rsc
from scipy import sparse

from workflows.gdtai.gdtai_v4_2_integration_core import read_json, resolve, sha256_file


CONFIG_PATH = PROJECT_ROOT / "configs/models/gdtai/v4_2_integration_execution.json"
REPORT_SCRIPT = Path(__file__).resolve()
BASE_NAME = "gdtai_v4_2_nk_reference_diagnostics"
TABLE_DIR = PROJECT_ROOT / f"Integrated_dataset/tables/gdT_prediction/{BASE_NAME}"
FIGURE_DIR = PROJECT_ROOT / f"Integrated_dataset/figures/gdT_prediction/{BASE_NAME}"
LOG_DIR = PROJECT_ROOT / f"Integrated_dataset/logs/gdT_prediction/{BASE_NAME}"
STATIC_DIR = PROJECT_ROOT / f"gdT_prediction/{BASE_NAME}"
RANDOM_SEED = 20260816
SAMPLE_TARGET = 250_000

COHORT_LABELS = {
    "current_atlas": "Existing atlas reference",
    "GSE114724": "GSE114724",
    "GSE159251": "GSE159251",
    "GSE292700": "GSE292700",
    "GSE294273_GSE294274": "GSE294273/4",
    "GSE296954": "GSE296954",
}
COHORT_COLORS = {
    "Existing atlas reference": "#8c979e",
    "GSE114724": "#16697a",
    "GSE159251": "#d98f39",
    "GSE292700": "#8f5da2",
    "GSE294273/4": "#2a9d6f",
    "GSE296954": "#c64b44",
}
CANDIDATE_COLORS = {
    "GSE114724": "#16697a",
    "GSE159251": "#d98f39",
    "GSE292700": "#8f5da2",
    "GSE294273": "#2a9d6f",
    "GSE296954": "#c64b44",
}
ROLE_COLORS = {
    "Other": "#cbd2d6",
    "Eligible candidate": "#d98f39",
    "Productive TCR anchor": "#247ba0",
    "Primary NK anchor": "#2a7f62",
}


def direct_cuda_environment() -> None:
    inherited = os.environ.get("CUDA_MPS_PIPE_DIRECTORY", "")
    if inherited == "/tmp/nvidia-mps" or (
        inherited and not inherited.startswith("/tmp/gdtai-v42-diagnostics-mps-")
    ):
        os.environ.pop("CUDA_MPS_PIPE_DIRECTORY", None)
    os.environ.setdefault("CUDA_MPS_PIPE_DIRECTORY", f"/tmp/gdtai-v42-diagnostics-mps-{os.getuid()}-{os.getpid()}")
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", "0")
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    os.environ.setdefault("PYTHONHASHSEED", str(RANDOM_SEED))


def balanced_sample(
    indices: np.ndarray,
    groups: np.ndarray,
    target: int,
    rng: np.random.Generator,
) -> np.ndarray:
    indices = np.asarray(indices, dtype=np.int64)
    if indices.size <= target:
        return indices
    levels = np.sort(np.unique(groups[indices]))
    quota = max(1, target // max(1, levels.size))
    selected: list[np.ndarray] = []
    leftovers: list[np.ndarray] = []
    for level in levels:
        group_indices = indices[groups[indices] == level]
        shuffled = rng.permutation(group_indices)
        selected.append(shuffled[: min(quota, shuffled.size)])
        leftovers.append(shuffled[min(quota, shuffled.size) :])
    result = np.concatenate(selected) if selected else np.empty(0, dtype=np.int64)
    remaining = target - result.size
    if remaining > 0:
        pool = np.concatenate(leftovers) if leftovers else np.empty(0, dtype=np.int64)
        result = np.concatenate([result, rng.choice(pool, size=remaining, replace=False)])
    return np.sort(result[:target])


def make_visualization_sample(obs: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(RANDOM_SEED)
    source = obs["source_gse_id"].astype(str).to_numpy()
    primary = obs["primary_nk_anchor"].astype(bool).to_numpy()
    productive = obs["productive_tcr_anchor"].astype(bool).to_numpy() & ~primary
    candidate = obs["candidate_eligible"].astype(bool).to_numpy() & ~primary & ~productive
    other = ~(primary | productive | candidate)

    selected_primary = np.flatnonzero(primary)
    selected_productive = balanced_sample(np.flatnonzero(productive), source, 60_000, rng)
    selected_candidate = balanced_sample(np.flatnonzero(candidate), source, 60_000, rng)
    remaining = SAMPLE_TARGET - selected_primary.size - selected_productive.size - selected_candidate.size
    if remaining <= 0:
        raise RuntimeError("Declared diagnostic sample target is smaller than fixed anchor samples")
    selected_other = balanced_sample(np.flatnonzero(other), source, remaining, rng)
    sample_indices = np.sort(
        np.concatenate([selected_primary, selected_productive, selected_candidate, selected_other])
    )
    if np.unique(sample_indices).size != SAMPLE_TARGET:
        raise RuntimeError("Diagnostic sample is not unique or did not reach its declared size")

    role = np.full(obs.shape[0], "Other", dtype=object)
    role[candidate] = "Eligible candidate"
    role[productive] = "Productive TCR anchor"
    role[primary] = "Primary NK anchor"
    return sample_indices, role


def partition_labels(saved: Any, run: str, n_cells: int) -> np.ndarray:
    global_names = saved["global_names"].astype(str).tolist()
    boundary_names = saved["boundary_names"].astype(str).tolist()
    if run in global_names:
        return np.asarray(saved["global_labels"][:, global_names.index(run)], dtype=np.int32)
    if run not in boundary_names:
        raise KeyError(f"Unknown partition run: {run}")
    labels = np.full(n_cells, -1, dtype=np.int32)
    labels[saved["boundary_indices"].astype(np.int64)] = saved["boundary_labels"][:, boundary_names.index(run)]
    return labels


def near_cluster_diagnostics(
    saved: Any,
    near: pd.DataFrame,
    obs: pd.DataFrame,
) -> tuple[np.ndarray, pd.DataFrame]:
    source = obs["source_gse_id"].astype(str).to_numpy()
    candidate = obs["candidate_eligible"].astype(bool).to_numpy()
    votes = np.zeros(obs.shape[0], dtype=np.int8)
    rows: list[dict[str, Any]] = []
    for record in near.itertuples(index=False):
        labels = partition_labels(saved, str(record.run), obs.shape[0])
        member = labels == int(record.cluster)
        votes += member.astype(np.int8)
        candidate_member = member & candidate
        values, counts = np.unique(source[candidate_member], return_counts=True)
        total = int(counts.sum())
        for value, count in zip(values, counts, strict=True):
            rows.append(
                {
                    "run": str(record.run),
                    "cluster": int(record.cluster),
                    "source_gse_id": str(value),
                    "n_candidate_cells": int(count),
                    "candidate_fraction": float(count / total) if total else 0.0,
                }
            )
    return votes, pd.DataFrame(rows)


def to_scipy_csr(matrix: Any) -> sparse.csr_matrix:
    if hasattr(matrix, "get"):
        matrix = matrix.get()
    return sparse.csr_matrix(matrix)


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


def scatter_categories(
    coordinates: np.ndarray,
    values: np.ndarray,
    order: list[str],
    colors: dict[str, str],
    title: str,
    path: Path,
    point_size: float = 1.2,
) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 6.5))
    for value in order:
        mask = values == value
        if not mask.any():
            continue
        alpha = 0.18 if value == "Other" else 0.58
        size = 0.55 if value == "Other" else point_size
        ax.scatter(
            coordinates[mask, 0],
            coordinates[mask, 1],
            s=size,
            alpha=alpha,
            c=colors[value],
            edgecolors="none",
            rasterized=True,
            label=f"{value} ({mask.sum():,})",
        )
    clean_umap_axis(ax, title)
    ax.legend(frameon=False, markerscale=5, fontsize=8, loc="best")
    finish_figure(fig, path)


def make_figures(
    sample: pd.DataFrame,
    cohort_summary: pd.DataFrame,
    near_source: pd.DataFrame,
    membership: pd.DataFrame,
    mixing_summary: pd.DataFrame,
) -> None:
    configure_plotting()
    coordinates = sample[["UMAP1", "UMAP2"]].to_numpy()

    cohort_values = sample["cohort_label"].to_numpy(str)
    scatter_categories(
        coordinates,
        cohort_values,
        list(COHORT_COLORS),
        COHORT_COLORS,
        "V4.2 development sidecar by input cohort",
        FIGURE_DIR / "diagnostic_umap_by_cohort.png",
        point_size=1.0,
    )

    role_values = sample["diagnostic_role"].to_numpy(str)
    scatter_categories(
        coordinates,
        role_values,
        ["Other", "Eligible candidate", "Productive TCR anchor", "Primary NK anchor"],
        ROLE_COLORS,
        "Anchor classes and candidates (anchor class is cohort-confounded)",
        FIGURE_DIR / "diagnostic_umap_by_anchor_role.png",
        point_size=1.35,
    )

    candidate_mask = sample["diagnostic_role"].eq("Eligible candidate").to_numpy()
    fig, ax = plt.subplots(figsize=(8.0, 6.5))
    ax.scatter(coordinates[:, 0], coordinates[:, 1], s=0.45, c="#d5dadd", alpha=0.10, edgecolors="none", rasterized=True)
    for source_name, color in CANDIDATE_COLORS.items():
        mask = candidate_mask & sample["source_gse_id"].eq(source_name).to_numpy()
        ax.scatter(
            coordinates[mask, 0], coordinates[mask, 1], s=1.7, c=color,
            alpha=0.65, edgecolors="none", rasterized=True,
            label=f"{source_name} ({mask.sum():,})",
        )
    clean_umap_axis(ax, "Eligible pseudo-NK candidates by source")
    ax.legend(frameon=False, markerscale=5, fontsize=8, loc="best")
    finish_figure(fig, FIGURE_DIR / "diagnostic_umap_candidate_sources.png")

    fig, ax = plt.subplots(figsize=(8.0, 6.5))
    zero = sample["near_cluster_memberships"].eq(0).to_numpy()
    ax.scatter(coordinates[zero, 0], coordinates[zero, 1], s=0.45, c="#d5dadd", alpha=0.10, edgecolors="none", rasterized=True)
    nonzero = ~zero
    plotted = ax.scatter(
        coordinates[nonzero, 0], coordinates[nonzero, 1],
        c=sample.loc[nonzero, "near_cluster_memberships"], cmap="viridis",
        vmin=1, vmax=6, s=1.7, alpha=0.72, edgecolors="none", rasterized=True,
    )
    colorbar = fig.colorbar(plotted, ax=ax, shrink=0.72, pad=0.02)
    colorbar.set_label("Memberships among 6 near-qualifying clusters")
    clean_umap_axis(ax, "Stability of the anchor-defined NK-like region")
    finish_figure(fig, FIGURE_DIR / "diagnostic_umap_near_cluster_membership.png")

    cohort_order = list(COHORT_LABELS.values())
    ordered = cohort_summary.set_index("cohort_label").reindex(cohort_order).reset_index()
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.3), gridspec_kw={"width_ratios": [0.92, 1.55]})
    bars = axes[0].barh(ordered["cohort_label"], ordered["n_cells"], color=[COHORT_COLORS[x] for x in ordered["cohort_label"]])
    axes[0].set_xscale("log")
    axes[0].bar_label(bars, labels=[f"{value:,}" for value in ordered["n_cells"]], padding=3, fontsize=8)
    axes[0].set_title("Full sidecar cell counts")
    axes[0].set_xlabel("Cells, log scale")
    role_columns = ["primary_nk_fraction", "productive_tcr_fraction", "candidate_fraction"]
    labels = ["Primary NK anchors", "Productive TCR anchors", "Eligible candidates"]
    colors = [ROLE_COLORS["Primary NK anchor"], ROLE_COLORS["Productive TCR anchor"], ROLE_COLORS["Eligible candidate"]]
    x = np.arange(ordered.shape[0])
    width = 0.25
    for offset, column, label, color in zip([-1, 0, 1], role_columns, labels, colors, strict=True):
        axes[1].bar(x + offset * width, 100 * ordered[column], width=width, label=label, color=color)
    axes[1].set_xticks(x, ordered["cohort_label"], rotation=25, ha="right")
    axes[1].set_ylabel("Cells (%)")
    axes[1].set_title("Anchor roles are confounded with cohort lane")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].grid(axis="y", color="#e1e7ea", linewidth=0.6)
    finish_figure(fig, FIGURE_DIR / "diagnostic_cohort_composition.png")

    pivot = near_source.pivot_table(
        index=["run", "cluster"], columns="source_gse_id", values="candidate_fraction", fill_value=0.0
    )
    pivot = pivot.reindex(columns=list(CANDIDATE_COLORS), fill_value=0.0)
    fig, ax = plt.subplots(figsize=(9.2, 4.4))
    bottom = np.zeros(pivot.shape[0])
    y = np.arange(pivot.shape[0])
    for source_name in pivot.columns:
        values = 100 * pivot[source_name].to_numpy()
        ax.barh(y, values, left=bottom, color=CANDIDATE_COLORS[source_name], label=source_name)
        bottom += values
    row_labels = [f"{run} / c{cluster}" for run, cluster in pivot.index]
    ax.set_yticks(y, row_labels)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Eligible candidate composition (%)")
    ax.set_title("Every NK-like cluster is dominated by GSE292700 candidates")
    ax.axvline(70, color="#b3261e", linestyle="--", linewidth=1.2, label="70% source cap")
    ax.legend(frameon=False, fontsize=8, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.36))
    finish_figure(fig, FIGURE_DIR / "diagnostic_near_cluster_source_composition.png")

    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    for role, color in ROLE_COLORS.items():
        subset = membership[membership["diagnostic_role"].eq(role)]
        ax.plot(
            subset["near_cluster_memberships"], 100 * subset["fraction"],
            marker="o", linewidth=2, color=color, label=role,
        )
    ax.set_xticks(range(7))
    ax.set_xlabel("Memberships among 6 near-qualifying clusters")
    ax.set_ylabel("Cells within role (%)")
    ax.set_title("The NK-anchor region is stable but cohort-confounded")
    ax.grid(color="#e1e7ea", linewidth=0.6)
    ax.legend(frameon=False, fontsize=8)
    finish_figure(fig, FIGURE_DIR / "diagnostic_near_cluster_membership_by_role.png")

    rng = np.random.default_rng(RANDOM_SEED)
    fig, ax = plt.subplots(figsize=(9.0, 4.6))
    violin_data = []
    labels = []
    for cohort in cohort_order:
        values = sample.loc[sample["cohort_label"].eq(cohort), "local_same_source_fraction"].to_numpy()
        if values.size > 5_000:
            values = rng.choice(values, size=5_000, replace=False)
        violin_data.append(values)
        labels.append(cohort)
    parts = ax.violinplot(violin_data, showmedians=True, showextrema=False)
    for body, cohort in zip(parts["bodies"], cohort_order, strict=True):
        body.set_facecolor(COHORT_COLORS[cohort])
        body.set_edgecolor("none")
        body.set_alpha(0.75)
    parts["cmedians"].set_color("#17212b")
    ax.set_xticks(np.arange(1, len(labels) + 1), labels, rotation=25, ha="right")
    ax.set_ylim(-0.02, 1.02)
    ax.set_ylabel("Same-source fraction among sampled kNN edges")
    ax.set_title("Local source retention in the diagnostic latent-space sample")
    ax.grid(axis="y", color="#e1e7ea", linewidth=0.6)
    finish_figure(fig, FIGURE_DIR / "diagnostic_local_source_retention.png")


def table_html(frame: pd.DataFrame, columns: list[str] | None = None) -> str:
    shown = frame if columns is None else frame[columns]
    return shown.to_html(index=False, border=0, classes="data-table", escape=True)


def render_report(
    summary: dict[str, Any],
    cohort_summary: pd.DataFrame,
    near_display: pd.DataFrame,
    mixing_display: pd.DataFrame,
) -> Path:
    figure_prefix = f"../../Integrated_dataset/figures/gdT_prediction/{BASE_NAME}"
    cohort_table = cohort_summary.copy()
    for column in ["primary_nk_fraction", "productive_tcr_fraction", "candidate_fraction"]:
        cohort_table[column] = cohort_table[column].map(lambda value: f"{value:.2%}")
    cohort_table = cohort_table.rename(
        columns={
            "cohort_label": "Cohort",
            "n_cells": "Cells",
            "n_primary_nk": "Primary NK",
            "primary_nk_fraction": "Primary NK %",
            "n_productive_tcr": "Productive TCR",
            "productive_tcr_fraction": "Productive TCR %",
            "n_candidate": "Eligible candidates",
            "candidate_fraction": "Candidate %",
        }
    )
    document = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>gdTAI V4.2 development sidecar diagnostics</title><style>
:root{{--ink:#17212b;--navy:#17324d;--teal:#197278;--red:#b3261e;--amber:#b86e16;--line:#ccd6dc;--soft:#f3f6f7}}
*{{box-sizing:border-box}}body{{font-family:Arial,Helvetica,sans-serif;color:var(--ink);line-height:1.48;max-width:1120px;margin:0 auto;padding:30px 42px;background:white}}h1{{font-size:29px;color:var(--navy);border-bottom:4px solid var(--teal);padding-bottom:12px}}h2{{font-size:20px;color:var(--navy);margin-top:31px}}p,li{{font-size:14px}}.status{{border-left:6px solid var(--amber);background:#fff9ef;padding:14px 18px;margin:18px 0}}.danger{{border-left:6px solid var(--red);background:#fff4f2;padding:13px 17px;margin:18px 0}}.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin:18px 0}}.card{{border:1px solid var(--line);padding:12px;background:var(--soft)}}.card b{{display:block;font-size:21px;color:var(--navy)}}.card span{{font-size:11px;color:#52636e}}.grid{{display:grid;grid-template-columns:1fr 1fr;gap:18px}}figure{{margin:18px 0;break-inside:avoid}}img{{display:block;width:100%;height:auto;margin:0 auto}}figcaption{{font-size:11px;color:#52636e;text-align:center;margin-top:7px}}.table-wrap{{overflow-x:auto;margin:14px 0 24px}}table{{border-collapse:collapse;width:100%;font-size:10px}}th{{background:var(--navy);color:white;text-align:left}}th,td{{padding:6px 7px;border:1px solid var(--line);vertical-align:top;overflow-wrap:anywhere}}tr:nth-child(even){{background:var(--soft)}}code{{background:#edf1f4;padding:1px 4px}}.warning{{color:var(--red);font-weight:bold}}
@media(max-width:800px){{body{{padding:18px}}.cards,.grid{{grid-template-columns:1fr}}}}
@media print{{@page{{size:A4;margin:9mm}}body{{max-width:none;padding:0}}h1{{font-size:23px}}h2{{font-size:17px;break-after:avoid}}p,li{{font-size:10.2px}}.cards{{grid-template-columns:repeat(4,1fr)}}.card b{{font-size:16px}}.card span{{font-size:8px}}.grid{{display:block}}figure{{break-inside:avoid}}figure img{{max-height:145mm;object-fit:contain}}table{{font-size:7.2px}}tr{{break-inside:avoid}}.table-wrap{{overflow:visible}}}}
</style></head><body>
<h1>gdTAI V4.2 development sidecar diagnostics</h1>
<p>This is a diagnostic visualization of the development-only scVI sidecar. It is not a new canonical atlas milestone and contains no locked validation cohort. Full-data counts use all {summary['n_cells']:,} cells; UMAP and local-neighbor plots use a deterministic, source-balanced {summary['n_visualization_cells']:,}-cell sample.</p>
<div class="status"><b>Interpretation:</b> the integration contains an NK-anchor-rich region, but the current method cannot confidently extend that region to unlabeled NK cells. Source balance, boundary localization, and anchor/cohort confounding all prevent a defensible pseudo-label set.</div>
<div class="cards"><div class="card"><b>{summary['n_cells']:,}</b><span>development cells</span></div><div class="card"><b>{summary['n_new_cohort_cells']:,}</b><span>new-cohort cells</span></div><div class="card"><b>{summary['n_primary_nk']:,}</b><span>primary NK anchors</span></div><div class="card"><b>{summary['n_selected_pseudo_nk']:,}</b><span>accepted pseudo-NK cells</span></div></div>

<h2>What was integrated</h2><figure><img src="{figure_prefix}/diagnostic_cohort_composition.png" alt="Cohort composition"><figcaption>Full-data cohort sizes and evidence classes. Candidate eligibility is intentionally restricted to the five new development cohorts.</figcaption></figure><div class="table-wrap">{table_html(cohort_table, ['Cohort','Cells','Primary NK','Primary NK %','Productive TCR','Productive TCR %','Eligible candidates','Candidate %'])}</div><div class="danger"><b>Critical confounding:</b> all {summary['n_primary_nk']:,} primary NK anchors are from the existing atlas reference, whereas all {summary['n_productive_tcr']:,} productive-TCR anchors are from the five new cohorts. The clustering cannot cleanly distinguish lineage structure from old-versus-new cohort structure.</div>

<h2>Latent-space structure</h2><div class="grid"><figure><img src="{figure_prefix}/diagnostic_umap_by_cohort.png" alt="UMAP by cohort"><figcaption>Diagnostic scVI UMAP by input cohort.</figcaption></figure><figure><img src="{figure_prefix}/diagnostic_umap_by_anchor_role.png" alt="UMAP by evidence role"><figcaption>Primary NK and productive-TCR anchors are expression-independent; candidates are no-TCR, non-doublet cells from permitted development cohorts.</figcaption></figure></div>
<div class="grid"><figure><img src="{figure_prefix}/diagnostic_umap_candidate_sources.png" alt="Candidate sources"><figcaption>Eligible candidates show strong source-dependent localization.</figcaption></figure><figure><img src="{figure_prefix}/diagnostic_umap_near_cluster_membership.png" alt="Near cluster membership"><figcaption>Membership frequency across the six clusters that passed anchor purity and TCR contamination criteria.</figcaption></figure></div>

<h2>Why pseudo-NK selection failed</h2><p>All six NK-anchor-rich clusters were low in productive-TCR contamination, but their eligible candidate cells were {summary['near_source_min']:.2%}-{summary['near_source_max']:.2%} GSE292700. This violates the frozen 70% maximum-source rule, so no cluster qualified and no candidate cell received a valid consensus label. Because the two anchor classes come from disjoint cohort lanes, anchor purity alone cannot prove lineage identity.</p><figure><img src="{figure_prefix}/diagnostic_near_cluster_source_composition.png" alt="Near cluster source composition"><figcaption>The dashed line is the frozen 70% maximum contribution from one candidate source.</figcaption></figure><figure><img src="{figure_prefix}/diagnostic_near_cluster_membership_by_role.png" alt="Membership by role"><figcaption>The anchor-defined region is stable across runs, but cohort-confounded anchors and source-skewed candidates prevent a defensible pseudo-label set.</figcaption></figure><div class="table-wrap">{table_html(near_display)}</div>

<h2>Integration mixing diagnostic</h2><p>The same-source-neighbor fraction is a descriptive diagnostic, not a biological pass/fail metric: true tissue, disease, and phenotype differences can legitimately retain source structure. It nevertheless shows where new cohorts remain locally isolated after scVI.</p><figure><img src="{figure_prefix}/diagnostic_local_source_retention.png" alt="Local source retention"><figcaption>Same-GSE fraction among sampled 30-nearest-neighbor graph edges.</figcaption></figure><div class="table-wrap">{table_html(mixing_display)}</div>

<h2>Conclusion</h2><p><b>An NK-anchor-enriched region is visible, but it cannot yet establish NK identity for unlabeled cells.</b> The next iteration must first break the anchor/cohort confounding, then use source-balanced candidate sampling and local anchor-neighborhood propagation with explicit rejection of source-specific neighborhoods. It should not loosen the source cap based on these results.</p><p>Classifier fitting performed: <b>no</b>. Source H5AD mutation: <b>none</b>. GitHub publication: <b>disabled until V4 is finished</b>.</p>
<h2>Reproducibility</h2><ul><li>Latent SHA-256: <code>{html.escape(summary['latent_sha256'])}</code></li><li>Partitions SHA-256: <code>{html.escape(summary['partitions_sha256'])}</code></li><li>Diagnostic sample seed: <code>{RANDOM_SEED}</code></li><li>Report script SHA-256: <code>{html.escape(summary['report_script_sha256'])}</code></li></ul>
</body></html>"""
    output = STATIC_DIR / "index.html"
    output.write_text(document, encoding="utf-8")
    return output


def main() -> None:
    direct_cuda_environment()
    if not torch.cuda.is_available():
        raise RuntimeError("A CUDA GPU is required; CPU fallback is prohibited")
    config = read_json(CONFIG_PATH)
    ssd_dir = Path(config["resources"]["ssd_root"])
    staged_path = ssd_dir / "development_hvg_counts.h5ad"
    latent_path = ssd_dir / "X_scVI.npy"
    partitions_path = ssd_dir / "cluster_partitions.npz"
    base_table_dir = resolve(config["outputs"]["table_dir"])
    base_log_dir = resolve(config["outputs"]["log_dir"])
    for directory in [TABLE_DIR, FIGURE_DIR, LOG_DIR, STATIC_DIR]:
        directory.mkdir(parents=True, exist_ok=True)

    fit_summary = read_json(base_log_dir / "fit_summary.json")
    cluster_summary = read_json(base_log_dir / "cluster_summary.json")
    consensus_summary = read_json(base_log_dir / "consensus_summary.json")
    if sha256_file(latent_path) != fit_summary["latent_sha256"]:
        raise RuntimeError("Latent checksum differs from the completed fit summary")
    if sha256_file(partitions_path) != cluster_summary["partitions_sha256"]:
        raise RuntimeError("Partition checksum differs from the completed cluster summary")

    backed = ad.read_h5ad(staged_path, backed="r")
    obs = backed.obs.copy()
    backed.file.close()
    if obs["input_cohort_id"].astype(str).isin(["GSE169246", "GSE315928", "GSE121636_GSE121637"]).any():
        raise RuntimeError("A locked cohort leaked into the diagnostic sidecar")
    latent = np.load(latent_path, mmap_mode="r")
    saved = np.load(partitions_path, allow_pickle=False)
    near = pd.read_csv(base_table_dir / "near_qualifying_clusters.csv")
    votes, near_source = near_cluster_diagnostics(saved, near, obs)
    sample_indices, full_role = make_visualization_sample(obs)

    sample_obs = obs.iloc[sample_indices].copy()
    sample_obs["diagnostic_role"] = full_role[sample_indices]
    sample_obs["cohort_label"] = sample_obs["input_cohort_id"].astype(str).map(COHORT_LABELS)
    sample_obs["near_cluster_memberships"] = votes[sample_indices]
    diagnostic = ad.AnnData(
        X=sparse.csr_matrix((sample_indices.size, 1), dtype=np.float32),
        obs=sample_obs,
    )
    diagnostic.obsm["X_scVI"] = np.asarray(latent[sample_indices], dtype=np.float32)
    rsc.pp.neighbors(
        diagnostic,
        n_neighbors=30,
        n_pcs=None,
        use_rep="X_scVI",
        random_state=RANDOM_SEED,
        algorithm="ivfflat",
        algorithm_kwds={"n_lists": 4096, "n_probes": 80},
    )
    rsc.tl.umap(diagnostic, min_dist=0.3, spread=1.0, random_state=RANDOM_SEED)
    coordinates = diagnostic.obsm["X_umap"]
    if hasattr(coordinates, "get"):
        coordinates = coordinates.get()
    coordinates = np.asarray(coordinates, dtype=np.float32)
    if coordinates.shape != (SAMPLE_TARGET, 2) or not np.isfinite(coordinates).all():
        raise RuntimeError("Diagnostic UMAP is missing, malformed, or non-finite")

    distances = to_scipy_csr(diagnostic.obsp["distances"])
    row_indices = np.repeat(np.arange(distances.shape[0]), np.diff(distances.indptr))
    sample_sources = sample_obs["source_gse_id"].astype(str).to_numpy()
    same_source = sample_sources[row_indices] == sample_sources[distances.indices]
    same_counts = np.bincount(row_indices, weights=same_source.astype(np.float32), minlength=distances.shape[0])
    degrees = np.diff(distances.indptr)
    same_fraction = np.divide(same_counts, degrees, out=np.zeros(distances.shape[0]), where=degrees > 0)

    sample = pd.DataFrame(
        {
            "integration_cell_id": sample_obs.index.astype(str),
            "input_cohort_id": sample_obs["input_cohort_id"].astype(str).to_numpy(),
            "cohort_label": sample_obs["cohort_label"].astype(str).to_numpy(),
            "source_gse_id": sample_sources,
            "diagnostic_role": sample_obs["diagnostic_role"].astype(str).to_numpy(),
            "near_cluster_memberships": sample_obs["near_cluster_memberships"].to_numpy(np.int8),
            "local_same_source_fraction": same_fraction,
            "UMAP1": coordinates[:, 0],
            "UMAP2": coordinates[:, 1],
        }
    )

    grouped = obs.assign(cohort_label=obs["input_cohort_id"].astype(str).map(COHORT_LABELS)).groupby("cohort_label", observed=True)
    cohort_summary = grouped.agg(
        n_cells=("input_cohort_id", "size"),
        n_primary_nk=("primary_nk_anchor", "sum"),
        n_productive_tcr=("productive_tcr_anchor", "sum"),
        n_candidate=("candidate_eligible", "sum"),
    ).reset_index()
    for count_column, fraction_column in [
        ("n_primary_nk", "primary_nk_fraction"),
        ("n_productive_tcr", "productive_tcr_fraction"),
        ("n_candidate", "candidate_fraction"),
    ]:
        cohort_summary[fraction_column] = cohort_summary[count_column] / cohort_summary["n_cells"]

    membership_rows: list[dict[str, Any]] = []
    for role in ROLE_COLORS:
        role_votes = votes[full_role == role]
        for value in range(7):
            count = int((role_votes == value).sum())
            membership_rows.append(
                {
                    "diagnostic_role": role,
                    "near_cluster_memberships": value,
                    "n_cells": count,
                    "fraction": float(count / role_votes.size) if role_votes.size else 0.0,
                }
            )
    membership = pd.DataFrame(membership_rows)
    mixing_summary = (
        sample.groupby("cohort_label", observed=True)["local_same_source_fraction"]
        .agg(n_sample_cells="size", median="median", mean="mean", p90=lambda values: values.quantile(0.9))
        .reset_index()
    )

    sample.to_csv(TABLE_DIR / "diagnostic_visualization_sample.csv.gz", index=False, compression="gzip")
    cohort_summary.to_csv(TABLE_DIR / "diagnostic_cohort_summary.csv", index=False)
    near_source.to_csv(TABLE_DIR / "diagnostic_near_cluster_source_composition.csv", index=False)
    membership.to_csv(TABLE_DIR / "diagnostic_near_cluster_membership_by_role.csv", index=False)
    mixing_summary.to_csv(TABLE_DIR / "diagnostic_local_source_retention.csv", index=False)
    np.savez_compressed(TABLE_DIR / "diagnostic_umap_coordinates.npz", sample_indices=sample_indices, X_umap=coordinates)

    make_figures(sample, cohort_summary, near_source, membership, mixing_summary)
    gse292700_near = near_source[near_source["source_gse_id"].eq("GSE292700")]
    summary = {
        "result": "PASS_DIAGNOSTIC_ONLY",
        "n_cells": int(obs.shape[0]),
        "n_existing_atlas_cells": int(obs["input_cohort_id"].astype(str).eq("current_atlas").sum()),
        "n_new_cohort_cells": int((~obs["input_cohort_id"].astype(str).eq("current_atlas")).sum()),
        "n_visualization_cells": int(sample.shape[0]),
        "n_primary_nk": int(obs["primary_nk_anchor"].sum()),
        "n_productive_tcr": int(obs["productive_tcr_anchor"].sum()),
        "n_candidate_eligible": int(obs["candidate_eligible"].sum()),
        "n_selected_pseudo_nk": int(consensus_summary["n_selected_pseudo_nk"]),
        "n_near_clusters": int(near.shape[0]),
        "near_source_min": float(gse292700_near["candidate_fraction"].min()),
        "near_source_max": float(gse292700_near["candidate_fraction"].max()),
        "locked_cohorts_included": False,
        "primary_nk_anchors_from_existing_atlas_fraction": float(
            obs.loc[obs["primary_nk_anchor"].astype(bool), "input_cohort_id"]
            .astype(str)
            .eq("current_atlas")
            .mean()
        ),
        "productive_tcr_anchors_from_new_cohorts_fraction": float(
            (~obs.loc[obs["productive_tcr_anchor"].astype(bool), "input_cohort_id"]
            .astype(str)
            .eq("current_atlas"))
            .mean()
        ),
        "anchor_cohort_confounding_detected": True,
        "classifier_fitting_performed": False,
        "source_h5ad_mutation_performed": False,
        "github_push_performed": False,
        "gpu": torch.cuda.get_device_name(0),
        "cpu_fallback": False,
        "latent_sha256": fit_summary["latent_sha256"],
        "partitions_sha256": cluster_summary["partitions_sha256"],
        "report_script_sha256": sha256_file(REPORT_SCRIPT),
        "random_seed": RANDOM_SEED,
    }
    (LOG_DIR / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    near_display = near[[
        "run", "cluster", "n_primary_nk_anchors", "n_productive_tcr_anchors",
        "anchor_nk_purity", "productive_tcr_contamination",
        "dominant_candidate_source", "maximum_candidate_source_fraction",
    ]].copy()
    for column in ["anchor_nk_purity", "productive_tcr_contamination", "maximum_candidate_source_fraction"]:
        near_display[column] = near_display[column].map(lambda value: f"{value:.2%}")
    near_display = near_display.rename(
        columns={
            "run": "Run",
            "cluster": "Cluster",
            "n_primary_nk_anchors": "NK anchors",
            "n_productive_tcr_anchors": "TCR anchors",
            "anchor_nk_purity": "NK purity",
            "productive_tcr_contamination": "TCR/cell",
            "dominant_candidate_source": "Dominant source",
            "maximum_candidate_source_fraction": "Source share",
        }
    )
    mixing_display = mixing_summary.copy()
    for column in ["median", "mean", "p90"]:
        mixing_display[column] = mixing_display[column].map(lambda value: f"{value:.2%}")
    mixing_display = mixing_display.rename(
        columns={
            "cohort_label": "Cohort",
            "n_sample_cells": "Sample cells",
            "median": "Median same-source",
            "mean": "Mean same-source",
            "p90": "90th percentile",
        }
    )
    report_html = render_report(summary, cohort_summary, near_display, mixing_display)

    chrome = shutil.which("google-chrome") or shutil.which("google-chrome-stable")
    if chrome is None:
        raise RuntimeError("Google Chrome is required for PDF export")
    profile = Path("/tmp/gdtai-v42-sidecar-diagnostics-chrome")
    profile.mkdir(parents=True, exist_ok=True)
    pdf_path = STATIC_DIR / f"{BASE_NAME}_report.pdf"
    subprocess.run(
        [
            chrome, "--headless", "--no-sandbox", "--disable-gpu",
            "--disable-dev-shm-usage", "--disable-breakpad", "--disable-crash-reporter",
            "--allow-file-access-from-files", "--no-pdf-header-footer",
            f"--user-data-dir={profile}", f"--print-to-pdf={pdf_path}",
            report_html.resolve().as_uri(),
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
