#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Build a detailed, read-only overview of the corrected full T/NK atlas."""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
import subprocess
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_H5AD = Path(
    "/ssd/tnk_phase3/Integrated_dataset/full_atlas/tcr_corrected/"
    "integrated_full_atlas.h5ad"
)
REPORT_DIR = ROOT / "full_atlas_overview"
FIGURE_DIR = ROOT / "Integrated_dataset/figures/full_atlas_overview"
TABLE_DIR = ROOT / "Integrated_dataset/tables/full_atlas_overview"
LOG_DIR = ROOT / "Integrated_dataset/logs/full_atlas_overview"

SEED = 20260818
PLOT_SAMPLE_CAP = 350_000
EXPRESSION_SAMPLE_CAP = 120_000
TCR_BACKGROUND_CAP = 300_000
AB_POSITIVE_CAP = 160_000
FIGURE_DPI = 220

MARKER_GENES = [
    "CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2",
    "CD4", "IL7R", "CCR7", "CD8A", "CD8B", "CCL5",
    "FOXP3", "IL2RA", "CTLA4", "TRDC", "TRDV1", "TRDV2",
    "TRGC1", "TRGC2", "NKG7", "KLRD1", "GNLY", "FCER1G",
]

SIGNATURES = {
    "Pan-T": ["CD3D", "CD3E", "CD3G", "TRAC"],
    "Alpha-beta T": ["TRAC", "TRBC1", "TRBC2"],
    "CD4/naive": ["CD4", "IL7R", "CCR7"],
    "CD8/cytotoxic": ["CD8A", "CD8B", "CCL5"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4"],
    "Gamma-delta": ["TRDC", "TRDV1", "TRDV2", "TRGC1", "TRGC2"],
    "NK/cytotoxic": ["NKG7", "KLRD1", "GNLY", "FCER1G"],
}

DISPLAY_ORDER = ["CD4 T", "CD8 T", "Treg", "MAIT", "NK", "gdT", "Other T", "Other", "Unannotated"]
DISPLAY_COLORS = {
    "CD4 T": "#2f6f9f", "CD8 T": "#2d8b57", "Treg": "#8a5aa5",
    "MAIT": "#d18f2b", "NK": "#b04a4a", "gdT": "#d62728",
    "Other T": "#5f6b73", "Other": "#a8a8a8", "Unannotated": "#dedede",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--skip-pdf", action="store_true")
    return parser.parse_args()


def decode(values: np.ndarray) -> np.ndarray:
    if values.dtype.kind not in "SO":
        return values
    return np.asarray(
        [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values],
        dtype=object,
    )


def obs_values(handle: h5py.File, key: str) -> np.ndarray:
    node = handle["obs"][key]
    if isinstance(node, h5py.Group):
        categories = decode(node["categories"][:])
        codes = np.asarray(node["codes"][:], dtype=np.int64)
        if categories.dtype.kind in "biuf":
            out = np.full(len(codes), np.nan, dtype=np.float64)
        else:
            out = np.full(len(codes), "", dtype=object)
        valid = (codes >= 0) & (codes < len(categories))
        out[valid] = categories[codes[valid]]
        return out
    return decode(node[:])


def axis_names(handle: h5py.File, axis: str) -> np.ndarray:
    key = handle[axis].attrs["_index"]
    if isinstance(key, bytes):
        key = key.decode("utf-8")
    return decode(handle[axis][str(key)][:])


def clean_text(values: np.ndarray) -> pd.Series:
    return pd.Series(values, dtype="string").fillna("").str.strip()


def source_label_group(label: str) -> str:
    value = str(label).strip().lower().replace("γ", "gamma").replace("δ", "delta")
    if value in {"", "false", "unknown", "unannotated", "nsclc_scrna", "na", "nan", "none"}:
        return "Unannotated"
    if any(token in value for token in ["gdt", "gamma_delta", "gamma delta", "gamma-delta"]):
        return "gdT"
    if "treg" in value or "regulatory" in value:
        return "Treg"
    if "mait" in value:
        return "MAIT"
    if value.startswith("nk") or "| nk |" in value or "b-nk" in value or "natural killer" in value:
        return "NK"
    if "cd8" in value:
        return "CD8 T"
    if "cd4" in value or "tfh" in value:
        return "CD4 T"
    if value in {"t", "tcell", "t_cell", "t cell"} or value.startswith("t |") or " t cell" in value:
        return "Other T"
    return "Other"


def source_balanced_sample(source: np.ndarray, cap: int, seed: int) -> np.ndarray:
    n = len(source)
    if n <= cap:
        return np.arange(n, dtype=np.int64)
    rng = np.random.default_rng(seed)
    codes, uniques = pd.factorize(pd.Series(source, dtype="string"), sort=True)
    counts = np.bincount(codes[codes >= 0], minlength=len(uniques))
    floor = min(1_000, max(1, cap // max(1, 4 * len(uniques))))
    quotas = np.minimum(counts, floor)
    remaining = cap - int(quotas.sum())
    if remaining > 0:
        weights = np.sqrt(counts.astype(float))
        weights[counts <= quotas] = 0
        if weights.sum() > 0:
            extra = np.floor(remaining * weights / weights.sum()).astype(int)
            quotas = np.minimum(counts, quotas + extra)
    selected: list[np.ndarray] = []
    for code in range(len(uniques)):
        rows = np.flatnonzero(codes == code)
        take = min(len(rows), int(quotas[code]))
        if take:
            selected.append(rng.choice(rows, size=take, replace=False))
    result = np.unique(np.concatenate(selected))
    if len(result) < cap:
        pool = np.setdiff1d(np.arange(n, dtype=np.int64), result, assume_unique=True)
        result = np.concatenate([result, rng.choice(pool, size=cap - len(result), replace=False)])
    elif len(result) > cap:
        result = rng.choice(result, size=cap, replace=False)
    return np.sort(result.astype(np.int64))


def safe_fraction(numerator: pd.Series, denominator: pd.Series) -> pd.Series:
    return numerator.astype(float).div(denominator.replace(0, np.nan)).fillna(0.0)


def exact_tables(handle: h5py.File, obs: dict[str, np.ndarray]) -> dict[str, pd.DataFrame]:
    n = len(obs["source_gse_id"])
    source = clean_text(obs["source_gse_id"]).replace("", "unresolved")
    annotation = clean_text(obs["source_cell_type_level1"]).map(source_label_group)
    truth = pd.DataFrame({
        "source_gse_id": source,
        "annotation_display_group": annotation,
        "has_TRA": np.asarray(obs["has_TRA"], dtype=bool),
        "has_TRB": np.asarray(obs["has_TRB"], dtype=bool),
        "paired_TRA_TRB": np.asarray(obs["has_TRA_TRB_paired"], dtype=bool),
        "has_TRG": np.asarray(obs["has_TRG"], dtype=bool),
        "has_TRD": np.asarray(obs["has_TRD"], dtype=bool),
        "paired_TRG_TRD": np.asarray(obs["has_TRG_TRD_paired"], dtype=bool),
        "any_ab_tcr": np.asarray(obs["has_any_ab_tcr"], dtype=bool),
        "any_gd_tcr": np.asarray(obs["has_any_gd_tcr"], dtype=bool),
        "sidecar_covered": np.asarray(obs["tcr_sidecar_covered_v2"], dtype=bool),
    })
    by_source = truth.groupby("source_gse_id", observed=True, sort=False).agg(
        n_cells=("source_gse_id", "size"),
        n_sidecar_covered=("sidecar_covered", "sum"),
        n_has_TRA=("has_TRA", "sum"), n_has_TRB=("has_TRB", "sum"),
        n_any_ab_tcr=("any_ab_tcr", "sum"), n_paired_TRA_TRB=("paired_TRA_TRB", "sum"),
        n_has_TRG=("has_TRG", "sum"), n_has_TRD=("has_TRD", "sum"),
        n_any_gd_tcr=("any_gd_tcr", "sum"), n_paired_TRG_TRD=("paired_TRG_TRD", "sum"),
    ).reset_index()
    for numerator in ["n_sidecar_covered", "n_any_ab_tcr", "n_paired_TRA_TRB", "n_any_gd_tcr", "n_paired_TRG_TRD"]:
        by_source[numerator.replace("n_", "fraction_")] = safe_fraction(by_source[numerator], by_source["n_cells"])
    by_source = by_source.sort_values("n_cells", ascending=False).reset_index(drop=True)

    annotation_counts = truth.groupby("annotation_display_group", observed=True).size().rename("n_cells").reset_index()
    annotation_counts["fraction"] = annotation_counts["n_cells"] / n
    annotation_counts["annotation_display_group"] = pd.Categorical(
        annotation_counts["annotation_display_group"], DISPLAY_ORDER, ordered=True
    )
    annotation_counts = annotation_counts.sort_values("annotation_display_group").reset_index(drop=True)

    cross = truth.groupby(["source_gse_id", "annotation_display_group"], observed=True).size().rename("n_cells").reset_index()
    cross["source_total"] = cross.groupby("source_gse_id")["n_cells"].transform("sum")
    cross["fraction_within_source"] = cross["n_cells"] / cross["source_total"]

    dimensions = []
    for key in ["input_cohort_id", "source_gse_id", "tissue_harmonized_v2", "specimen_context_harmonized_v2", "diagnosis", "technology_simple", "leiden"]:
        values = clean_text(obs[key]).replace("", "unresolved")
        counts = values.value_counts(dropna=False).rename_axis("value").reset_index(name="n_cells")
        counts.insert(0, "dimension", key)
        counts["fraction"] = counts["n_cells"] / n
        dimensions.append(counts)
    dimension_counts = pd.concat(dimensions, ignore_index=True)

    missing_rows = []
    for key in ["source_gse_id", "sample_id_harmonized_v2", "library_id_harmonized_v2", "donor_id_harmonized_v2", "tissue_harmonized_v2", "specimen_context_harmonized_v2", "diagnosis"]:
        values = clean_text(obs[key]).str.lower()
        missing = values.isin({"", "na", "nan", "none", "unknown", "unresolved", "unassigned"})
        missing_rows.append({"field": key, "n_missing_or_unresolved": int(missing.sum()), "fraction": float(missing.mean())})
    metadata_missingness = pd.DataFrame(missing_rows)

    overall = pd.DataFrame([{
        "n_cells": n,
        "n_genes": int(handle["X"].attrs["shape"][1]),
        "n_source_accessions": int(source.nunique()),
        "n_input_cohorts": int(clean_text(obs["input_cohort_id"]).nunique()),
        "n_leiden_clusters": int(clean_text(obs["leiden"]).nunique()),
        "n_tissues": int(clean_text(obs["tissue_harmonized_v2"]).replace("unresolved", "").replace("", np.nan).nunique()),
        "n_any_ab_tcr": int(truth["any_ab_tcr"].sum()),
        "n_paired_TRA_TRB": int(truth["paired_TRA_TRB"].sum()),
        "n_any_gd_tcr": int(truth["any_gd_tcr"].sum()),
        "n_paired_TRG_TRD": int(truth["paired_TRG_TRD"].sum()),
        "n_sidecar_covered": int(truth["sidecar_covered"].sum()),
    }])
    return {
        "overall": overall, "by_source": by_source, "annotation_counts": annotation_counts,
        "annotation_by_source": cross, "dimension_counts": dimension_counts,
        "metadata_missingness": metadata_missingness, "truth": truth,
    }


def extract_marker_expression(
    handle: h5py.File, rows: np.ndarray, marker_genes: list[str]
) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    var_names = pd.Index(axis_names(handle, "var").astype(str))
    present = [gene for gene in marker_genes if gene in var_names]
    marker_idx = var_names.get_indexer(present).astype(np.int64)
    order = np.argsort(marker_idx)
    sorted_idx = marker_idx[order]
    sorted_pos = np.arange(len(marker_idx), dtype=np.int64)[order]
    result = np.zeros((len(rows), len(present)), dtype=np.float32)
    row_sums = np.zeros(len(rows), dtype=np.float64)
    detected_genes = np.zeros(len(rows), dtype=np.int32)
    indptr, indices, data = handle["X/indptr"], handle["X/indices"], handle["X/data"]
    for out_row, row in enumerate(rows):
        start, end = int(indptr[row]), int(indptr[row + 1])
        row_idx = indices[start:end]
        row_data = data[start:end]
        row_sums[out_row] = float(np.sum(row_data, dtype=np.float64))
        detected_genes[out_row] = int(np.count_nonzero(row_data))
        positions = np.searchsorted(row_idx, sorted_idx)
        valid = positions < len(row_idx)
        matched = np.zeros(len(sorted_idx), dtype=bool)
        matched[valid] = row_idx[positions[valid]] == sorted_idx[valid]
        if np.any(matched):
            result[out_row, sorted_pos[matched]] = row_data[positions[matched]]
    scale = 10_000.0 / np.maximum(row_sums, 1.0)
    result = np.log1p(result * scale[:, None])
    return present, result, row_sums, detected_genes


def style_axes(ax: plt.Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)


def plot_composition(by_source: pd.DataFrame, annotation: pd.DataFrame) -> list[Path]:
    outputs = []
    frame = by_source.sort_values("n_cells")
    fig, ax = plt.subplots(figsize=(10, 10), constrained_layout=True)
    ax.barh(frame["source_gse_id"], frame["n_cells"], color="#356a7c")
    ax.set_xlabel("Cells")
    ax.set_title("Atlas cells by source accession")
    ax.xaxis.set_major_formatter(lambda x, _pos: f"{x/1e6:.1f}M" if x >= 1e6 else f"{x/1e3:.0f}k")
    style_axes(ax)
    path = FIGURE_DIR / "cells_by_source.png"
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig); outputs.append(path)

    frame = annotation.copy()
    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    colors = [DISPLAY_COLORS.get(str(x), "#999999") for x in frame["annotation_display_group"]]
    ax.bar(frame["annotation_display_group"].astype(str), frame["n_cells"], color=colors)
    ax.set_ylabel("Cells")
    ax.set_title("Source-author annotation display groups")
    ax.tick_params(axis="x", rotation=35)
    style_axes(ax)
    path = FIGURE_DIR / "annotation_group_counts.png"
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig); outputs.append(path)
    return outputs


def plot_tcr_coverage(by_source: pd.DataFrame) -> Path:
    frame = by_source.sort_values("n_cells")
    fig, axes = plt.subplots(1, 2, figsize=(14, 11), constrained_layout=True, sharey=True)
    y = np.arange(len(frame))
    axes[0].barh(y - 0.18, frame["fraction_any_ab_tcr"], height=0.34, color="#356a7c", label="Any TRA/TRB")
    axes[0].barh(y + 0.18, frame["fraction_paired_TRA_TRB"], height=0.34, color="#78a9bb", label="Paired TRA/TRB")
    axes[1].barh(y - 0.18, frame["fraction_any_gd_tcr"], height=0.34, color="#9b3a47", label="Any TRG/TRD")
    axes[1].barh(y + 0.18, frame["fraction_paired_TRG_TRD"], height=0.34, color="#d7838c", label="Paired TRG/TRD")
    for ax, title in zip(axes, ["Alpha-beta TCR metadata", "Gamma-delta TCR metadata"]):
        ax.set_yticks(y, frame["source_gse_id"])
        ax.set_xlim(0, 1)
        ax.set_xlabel("Fraction of all atlas cells in source")
        ax.set_title(title)
        ax.legend(frameon=False, loc="lower right")
        style_axes(ax)
    path = FIGURE_DIR / "tcr_coverage_by_source.png"
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def categorical_umap(umap: np.ndarray, rows: np.ndarray, labels: np.ndarray, title: str, path: Path) -> Path:
    frame_labels = np.asarray(labels[rows], dtype=object)
    categories = [x for x in DISPLAY_ORDER if x in set(frame_labels)] if "annotation" in path.name else sorted(set(frame_labels))
    cmap = plt.get_cmap("tab20")
    fig, ax = plt.subplots(figsize=(11, 8), constrained_layout=True)
    for i, category in enumerate(categories):
        take = frame_labels == category
        color = DISPLAY_COLORS.get(category, cmap(i % 20))
        ax.scatter(umap[rows[take], 0], umap[rows[take], 1], s=0.4, alpha=0.6, linewidths=0, rasterized=True, color=color, label=category)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2"); ax.set_title(title)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal", adjustable="datalim")
    ax.legend(markerscale=8, fontsize=7, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left", ncol=2 if len(categories) > 20 else 1)
    for spine in ax.spines.values(): spine.set_visible(False)
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def tcr_umap(umap: np.ndarray, background: np.ndarray, masks: list[tuple[str, np.ndarray]], path: Path) -> Path:
    rng = np.random.default_rng(SEED + 11)
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    for ax, (title, mask) in zip(axes.ravel(), masks, strict=True):
        ax.scatter(umap[background, 0], umap[background, 1], s=0.3, color="#d8d8d8", alpha=0.35, linewidths=0, rasterized=True)
        positives = np.flatnonzero(mask)
        cap = AB_POSITIVE_CAP if "TRA/TRB" in title else len(positives)
        if len(positives) > cap:
            positives = rng.choice(positives, size=cap, replace=False)
        ax.scatter(umap[positives, 0], umap[positives, 1], s=0.65, color="#b42338", alpha=0.68, linewidths=0, rasterized=True)
        ax.set_title(f"{title}\n{int(mask.sum()):,} cells ({mask.mean():.2%})")
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values(): spine.set_visible(False)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def marker_umaps(umap: np.ndarray, rows: np.ndarray, genes: list[str], expression: np.ndarray) -> Path:
    selected = ["CD3D", "TRAC", "TRBC1", "CD4", "CD8A", "FOXP3", "TRDC", "TRDV1", "TRDV2", "TRGC1", "NKG7", "KLRD1", "GNLY", "FCER1G"]
    selected = [gene for gene in selected if gene in genes]
    ncol = 4; nrow = int(np.ceil(len(selected) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(14, 3.3 * nrow), constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, gene in zip(axes, selected, strict=False):
        values = expression[:, genes.index(gene)]
        order = np.argsort(values)
        vmax = max(float(np.quantile(values, 0.99)), 0.25)
        scatter = ax.scatter(umap[rows[order], 0], umap[rows[order], 1], c=values[order], s=0.45, cmap="magma", norm=Normalize(0, vmax), linewidths=0, rasterized=True)
        ax.set_title(gene)
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values(): spine.set_visible(False)
        fig.colorbar(scatter, ax=ax, fraction=0.04, pad=0.01).ax.tick_params(labelsize=6)
    for ax in axes[len(selected):]: ax.axis("off")
    fig.suptitle("Signature genes on the integrated UMAP (log1p normalized counts)", fontsize=15)
    path = FIGURE_DIR / "signature_gene_umaps.png"
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    return path


def signature_dotplot(leiden: np.ndarray, rows: np.ndarray, genes: list[str], expression: np.ndarray) -> tuple[Path, pd.DataFrame]:
    records = []
    clusters = sorted(set(leiden[rows]), key=lambda x: int(x) if str(x).isdigit() else str(x))
    for cluster in clusters:
        take = leiden[rows] == cluster
        for signature, members in SIGNATURES.items():
            positions = [genes.index(gene) for gene in members if gene in genes]
            values = expression[take][:, positions]
            records.append({
                "leiden": str(cluster), "signature": signature,
                "mean_log1p_normalized_expression": float(values.mean()) if values.size else np.nan,
                "mean_gene_detection_fraction": float((values > 0).mean()) if values.size else np.nan,
                "sampled_cells": int(take.sum()), "genes_present": ", ".join([genes[p] for p in positions]),
            })
    table = pd.DataFrame(records)
    matrix_mean = table.pivot(index="leiden", columns="signature", values="mean_log1p_normalized_expression").loc[[str(x) for x in clusters]]
    matrix_det = table.pivot(index="leiden", columns="signature", values="mean_gene_detection_fraction").loc[[str(x) for x in clusters]]
    long = matrix_mean.stack().rename("mean").reset_index().merge(matrix_det.stack().rename("det").reset_index(), on=["leiden", "signature"])
    x_order = list(SIGNATURES); y_order = [str(x) for x in clusters]
    x_map = {x: i for i, x in enumerate(x_order)}; y_map = {y: i for i, y in enumerate(y_order)}
    fig, ax = plt.subplots(figsize=(11, 10), constrained_layout=True)
    scatter = ax.scatter(long["signature"].map(x_map), long["leiden"].map(y_map), s=20 + 420 * long["det"], c=long["mean"], cmap="viridis", linewidths=0.3, edgecolors="white")
    ax.set_xticks(range(len(x_order)), x_order, rotation=35, ha="right")
    ax.set_yticks(range(len(y_order)), y_order)
    ax.invert_yaxis(); ax.set_xlabel(""); ax.set_ylabel("Leiden cluster")
    ax.set_title("Lineage-signature expression by unsupervised cluster")
    fig.colorbar(scatter, ax=ax, label="Mean log1p normalized expression")
    for size, label in [(0.05, "5%"), (0.25, "25%"), (0.50, "50%")]:
        ax.scatter([], [], s=20 + 420 * size, color="#777777", label=label)
    ax.legend(title="Mean gene detection", frameon=False, bbox_to_anchor=(1.17, 0.35))
    style_axes(ax)
    path = FIGURE_DIR / "signature_dotplot_by_leiden.png"
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    return path, table


def qc_plot(
    source: np.ndarray,
    rows: np.ndarray,
    obs: dict[str, np.ndarray],
    matrix_total_counts: np.ndarray,
    matrix_detected_genes: np.ndarray,
) -> tuple[Path, pd.DataFrame]:
    frame = pd.DataFrame({
        "source_gse_id": source[rows], "total_counts": matrix_total_counts,
        "n_genes_by_counts": matrix_detected_genes,
        "pct_counts_mt": np.asarray(obs["pct_counts_mt"])[rows], "pct_counts_ribo": np.asarray(obs["pct_counts_ribo"])[rows],
    })
    summary = frame.groupby("source_gse_id", observed=True).agg(
        sampled_cells=("source_gse_id", "size"), median_total_counts=("total_counts", "median"),
        median_n_genes=("n_genes_by_counts", "median"), median_pct_mt=("pct_counts_mt", "median"),
        median_pct_ribo=("pct_counts_ribo", "median"),
    ).reset_index().sort_values("median_n_genes")
    fig, axes = plt.subplots(1, 3, figsize=(15, 10), constrained_layout=True, sharey=True)
    y = np.arange(len(summary))
    for ax, column, title in zip(axes, ["median_total_counts", "median_n_genes", "median_pct_mt"], ["Median total counts", "Median detected genes", "Median mitochondrial %"]):
        ax.scatter(summary[column], y, s=22, color="#356a7c")
        ax.set_yticks(y, summary["source_gse_id"])
        ax.set_xlabel(title); style_axes(ax)
    fig.suptitle("Source-level QC medians (source-balanced descriptive sample)")
    path = FIGURE_DIR / "qc_medians_by_source.png"
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    return path, summary


def image_card(title: str, filename: str, caption: str) -> str:
    return f"<figure><h3>{html.escape(title)}</h3><img src='assets/{html.escape(filename)}' alt='{html.escape(title)}'><figcaption>{html.escape(caption)}</figcaption></figure>"


def table_html(frame: pd.DataFrame, rows: int = 50) -> str:
    shown = frame.head(rows).copy()
    for column in shown.select_dtypes(include="float").columns:
        shown[column] = shown[column].map(lambda value: f"{value:.4g}" if pd.notna(value) else "")
    return f"<div class='table-wrap'>{shown.to_html(index=False, border=0, escape=True)}</div>"


def render_report(h5ad: Path, tables: dict[str, pd.DataFrame], figure_paths: list[Path], input_signature: dict[str, object], skip_pdf: bool) -> tuple[Path, Path | None]:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    asset_dir = REPORT_DIR / "assets"; asset_dir.mkdir(parents=True, exist_ok=True)
    for source in figure_paths:
        target = asset_dir / source.name
        if target.exists() or target.is_symlink(): target.unlink()
        os.symlink(os.path.relpath(source, asset_dir), target)
    overall = tables["overall"].iloc[0]
    source_table = tables["by_source"].copy()
    source_table = source_table[["source_gse_id", "n_cells", "n_any_ab_tcr", "n_paired_TRA_TRB", "n_any_gd_tcr", "n_paired_TRG_TRD", "fraction_any_ab_tcr", "fraction_any_gd_tcr"]]
    dimensions = tables["dimension_counts"]
    tissue = dimensions[dimensions["dimension"] == "tissue_harmonized_v2"]
    specimen = dimensions[dimensions["dimension"] == "specimen_context_harmonized_v2"]
    diagnosis = dimensions[dimensions["dimension"] == "diagnosis"]
    body = f"""<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>
<title>Full T/NK Atlas Overview</title><style>
:root{{--ink:#20252b;--muted:#5d6670;--line:#d8dde2;--paper:#fff;--bg:#eef2f3;--red:#9b3a47;--teal:#356a7c}}*{{box-sizing:border-box}}body{{margin:0;background:var(--bg);color:var(--ink);font-family:Arial,Helvetica,sans-serif;line-height:1.5;letter-spacing:0}}main{{width:min(1240px,calc(100% - 32px));margin:18px auto 48px}}header{{background:#18333d;color:white;padding:34px 38px;border-radius:6px}}h1{{font-size:40px;margin:0 0 10px}}header p{{max-width:920px;margin:6px 0;color:#e5eef1}}.metrics{{display:grid;grid-template-columns:repeat(5,1fr);gap:10px;margin-top:22px}}.metric{{border-top:2px solid #77aabb;padding-top:9px}}.metric b{{display:block;font-size:23px}}.metric span{{font-size:12px;color:#d7e4e8}}section{{background:var(--paper);padding:26px 30px;margin-top:16px;border:1px solid var(--line);border-radius:6px}}h2{{font-size:25px;margin:0 0 8px}}h3{{font-size:16px;margin:0 0 8px}}p.note{{background:#f2f6f7;border-left:4px solid var(--teal);padding:12px 14px;color:var(--muted)}}.grid{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px}}figure{{margin:14px 0 0;border-top:1px solid var(--line);padding-top:14px;break-inside:avoid}}figure img{{display:block;width:100%;height:auto}}figcaption{{font-size:12px;color:var(--muted);margin-top:7px}}.table-wrap{{overflow:auto;max-height:560px;border:1px solid var(--line)}}table{{border-collapse:collapse;width:100%;font-size:11px}}th,td{{padding:6px 7px;border-bottom:1px solid #e4e7e9;text-align:left;white-space:nowrap}}th{{position:sticky;top:0;background:#edf2f4}}code{{font-size:12px;word-break:break-all}}@media(max-width:850px){{.metrics{{grid-template-columns:repeat(2,1fr)}}.grid{{grid-template-columns:1fr}}}}@media print{{body{{background:white}}main{{width:100%;margin:0}}header,section{{border-radius:0;box-shadow:none}}.table-wrap{{overflow:visible;max-height:none}}table{{font-size:8px}}.page-break{{break-before:page}}}}
@media print{{h2,h3{{break-after:avoid-page}}}}</style></head><body><main><header><h1>Full T/NK Atlas Overview</h1><p>A read-only overview of the rebuilt 5.9-million-cell atlas using the validated metadata- and TCR-corrected candidate. This report separates source-author annotations, unsupervised clusters, TCR metadata, and expression evidence.</p><div class='metrics'>
<div class='metric'><b>{int(overall.n_cells):,}</b><span>cells</span></div><div class='metric'><b>{int(overall.n_source_accessions):,}</b><span>source accessions</span></div><div class='metric'><b>{int(overall.n_genes):,}</b><span>genes</span></div><div class='metric'><b>{int(overall.n_paired_TRA_TRB):,}</b><span>paired TRA/TRB</span></div><div class='metric'><b>{int(overall.n_paired_TRG_TRD):,}</b><span>paired TRG/TRD</span></div></div></header>
<section><h2>Interpretation boundaries</h2><p class='note'><b>Annotation:</b> the atlas lacks a complete scANVI label column. “Annotation display groups” below are deterministic broad groupings of available source-author labels, not model-derived consensus labels. Leiden clusters and marker expression are shown independently.</p><p class='note'><b>TCR:</b> exact counts use the validated repaired sidecar. Absence of TRG/TRD is not evidence against gamma-delta identity in sources without gamma-delta V(D)J sequencing. Coverage percentages use all atlas cells in each source as denominator.</p><p class='note'><b>Expression:</b> X contains raw counts. Marker panels use log1p(10,000 × gene count / cell total counts) on a fixed, source-balanced sample of {EXPRESSION_SAMPLE_CAP:,} cells. These panels are descriptive and are not cell labels.</p></section>
<section><h2>Dataset composition</h2><div class='grid'>{image_card('Cells by source', 'cells_by_source.png', 'Exact all-cell counts for all 40 source accessions.')}{image_card('Author-label display groups', 'annotation_group_counts.png', 'Exact all-cell counts after deterministic broad grouping of available source-author labels.')}</div><h3>Source-level exact counts</h3>{table_html(source_table, 45)}</section>
<section class='page-break'><h2>Integrated structure and annotation evidence</h2><div class='grid'>{image_card('Unsupervised clusters', 'umap_by_leiden.png', 'Source-balanced display sample colored by the 33 Leiden clusters from the rebuilt integration.')}{image_card('Source-author annotation groups', 'umap_by_annotation.png', 'The same manifold colored only by broad groupings of source-author labels; unannotated cells remain explicit.')}</div></section>
<section><h2>TCR coverage</h2><div class='grid'>{image_card('TCR metadata on UMAP', 'tcr_coverage_umaps.png', 'Exact positive totals are stated in panel titles; alpha-beta positives are display-sampled when needed, while all gamma-delta positives are shown.')}{image_card('TCR coverage by source', 'tcr_coverage_by_source.png', 'Exact fractions of all retained atlas cells in each source carrying repaired productive-chain evidence.')}</div>{table_html(source_table, 45)}</section>
<section class='page-break'><h2>Signature gene expression</h2>{image_card('Marker expression UMAPs', 'signature_gene_umaps.png', 'Source-balanced cells; expression is temporary library-size normalization of raw counts followed by log1p.')}{image_card('Signature dot plot by Leiden cluster', 'signature_dotplot_by_leiden.png', 'Dot color is mean expression and dot size is mean per-gene detection across each curated signature.')}</section>
<section><h2>QC and metadata coverage</h2>{image_card('QC medians by source', 'qc_medians_by_source.png', 'Source-level medians computed on the same source-balanced descriptive sample; exact sample sizes are exported.')}
<div class='grid'><div><h3>Tissue</h3>{table_html(tissue, 35)}</div><div><h3>Specimen context</h3>{table_html(specimen, 20)}</div><div><h3>Diagnosis</h3>{table_html(diagnosis, 30)}</div><div><h3>Missing or unresolved metadata</h3>{table_html(tables['metadata_missingness'], 20)}</div></div></section>
<section><h2>Provenance and reproducibility</h2><p>Input H5AD: <code>{html.escape(str(h5ad))}</code></p><p>Input size/mtime signature: <code>{html.escape(json.dumps(input_signature, sort_keys=True))}</code></p><p>All exact count tables, sampled QC summaries, the signature matrix, and the report manifest are exported under <code>Integrated_dataset/tables/full_atlas_overview/</code> and <code>Integrated_dataset/logs/full_atlas_overview/</code>. No H5AD field or canonical link was modified.</p></section>
</main></body></html>"""
    index = REPORT_DIR / "index.html"; index.write_text(body, encoding="utf-8")
    pdf = None
    if not skip_pdf:
        pdf = REPORT_DIR / "full_atlas_overview_report.pdf"
        subprocess.run(["google-chrome", "--headless", "--no-sandbox", "--disable-gpu", "--print-to-pdf-no-header", f"--print-to-pdf={pdf}", str(index)], check=True)
    return index, pdf


def main() -> None:
    args = parse_args(); h5ad = args.h5ad.resolve()
    for directory in [FIGURE_DIR, TABLE_DIR, LOG_DIR, REPORT_DIR]: directory.mkdir(parents=True, exist_ok=True)
    before = {"size": h5ad.stat().st_size, "mtime_ns": h5ad.stat().st_mtime_ns}
    needed = ["source_gse_id", "input_cohort_id", "source_cell_type_level1", "leiden", "tissue_harmonized_v2", "specimen_context_harmonized_v2", "diagnosis", "technology_simple", "sample_id_harmonized_v2", "library_id_harmonized_v2", "donor_id_harmonized_v2", "has_TRA", "has_TRB", "has_TRA_TRB_paired", "has_TRG", "has_TRD", "has_TRG_TRD_paired", "has_any_ab_tcr", "has_any_gd_tcr", "tcr_sidecar_covered_v2", "total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo"]
    with h5py.File(h5ad, "r") as handle:
        missing = [key for key in needed if key not in handle["obs"]]
        if missing: raise KeyError(f"Missing required obs columns: {missing}")
        obs = {key: obs_values(handle, key) for key in needed}
        tables = exact_tables(handle, obs)
        source = clean_text(obs["source_gse_id"]).replace("", "unresolved").to_numpy(dtype=object)
        leiden = clean_text(obs["leiden"]).to_numpy(dtype=object)
        annotation = clean_text(obs["source_cell_type_level1"]).map(source_label_group).to_numpy(dtype=object)
        plot_rows = source_balanced_sample(source, PLOT_SAMPLE_CAP, SEED)
        expression_rows = source_balanced_sample(source, EXPRESSION_SAMPLE_CAP, SEED + 1)
        umap = np.asarray(handle["obsm/X_umap"][:], dtype=np.float32)
        genes, expression, matrix_total_counts, matrix_detected_genes = extract_marker_expression(
            handle, expression_rows, MARKER_GENES
        )

    for name, frame in tables.items():
        if name != "truth": frame.to_csv(TABLE_DIR / f"{name}.csv", index=False)
    pd.DataFrame({"row_index": plot_rows, "source_gse_id": source[plot_rows]}).to_csv(TABLE_DIR / "umap_display_sample.csv.gz", index=False, compression="gzip")
    pd.DataFrame({"row_index": expression_rows, "source_gse_id": source[expression_rows]}).to_csv(TABLE_DIR / "expression_display_sample.csv.gz", index=False, compression="gzip")

    figures = plot_composition(tables["by_source"], tables["annotation_counts"])
    figures.append(plot_tcr_coverage(tables["by_source"]))
    figures.append(categorical_umap(umap, plot_rows, leiden, "Unsupervised Leiden clusters", FIGURE_DIR / "umap_by_leiden.png"))
    figures.append(categorical_umap(umap, plot_rows, annotation, "Source-author annotation display groups", FIGURE_DIR / "umap_by_annotation.png"))
    rng = np.random.default_rng(SEED + 2)
    background = np.sort(rng.choice(np.arange(len(source)), size=min(TCR_BACKGROUND_CAP, len(source)), replace=False))
    truth = tables["truth"]
    figures.append(tcr_umap(umap, background, [
        ("Any productive TRA/TRB", truth["any_ab_tcr"].to_numpy()),
        ("Paired productive TRA/TRB", truth["paired_TRA_TRB"].to_numpy()),
        ("Any productive TRG/TRD", truth["any_gd_tcr"].to_numpy()),
        ("Paired productive TRG/TRD", truth["paired_TRG_TRD"].to_numpy()),
    ], FIGURE_DIR / "tcr_coverage_umaps.png"))
    figures.append(marker_umaps(umap, expression_rows, genes, expression))
    dot_path, signature_table = signature_dotplot(leiden, expression_rows, genes, expression); figures.append(dot_path)
    signature_table.to_csv(TABLE_DIR / "signature_summary_by_leiden.csv", index=False)
    qc_path, qc_table = qc_plot(
        source, expression_rows, obs, matrix_total_counts, matrix_detected_genes
    ); figures.append(qc_path)
    qc_table.to_csv(TABLE_DIR / "qc_summary_by_source.csv", index=False)

    after = {"size": h5ad.stat().st_size, "mtime_ns": h5ad.stat().st_mtime_ns}
    if before != after: raise RuntimeError("Input H5AD size/mtime changed during read-only report generation")
    index, pdf = render_report(h5ad, tables, figures, before, args.skip_pdf)
    manifest = {
        "status": "PASS", "input_h5ad": str(h5ad), "input_signature_before": before,
        "input_signature_after": after, "input_unchanged": before == after,
        "html": str(index), "pdf": str(pdf) if pdf else None,
        "figures": [str(path) for path in figures],
        "tables": [str(path) for path in sorted(TABLE_DIR.glob("*"))],
        "plot_sample_cells": int(len(plot_rows)), "expression_sample_cells": int(len(expression_rows)),
        "random_seed": SEED, "expression_scale": "log1p(10000 * raw_gene_count / total_counts)",
    }
    (LOG_DIR / "report_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    summary = ["# Full T/NK atlas overview", "", "- Status: `PASS`", f"- Cells: `{int(tables['overall'].iloc[0]['n_cells']):,}`", f"- HTML: `{index}`", f"- PDF: `{pdf}`", f"- Input H5AD unchanged: `{before == after}`"]
    (LOG_DIR / "summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
