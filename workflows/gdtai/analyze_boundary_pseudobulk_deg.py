#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Paired pseudobulk DEG for productive-AB-enriched versus depleted NK-boundary clusters."""

from __future__ import annotations

import argparse
import gzip
import html
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse
from scipy.io import mmwrite
from scipy.stats import fisher_exact


ROOT = Path(__file__).resolve().parents[2]
CORRECTED = Path("/ssd/tnk_phase3/Integrated_dataset/full_atlas/tcr_corrected/integrated_full_atlas.h5ad")
STAGED = Path("/ssd/tnk_phase3/Integrated_dataset/gdtai_v4_2_tnk_reintegration/tnk_refined_hvg_counts.h5ad")
BOUNDARY = Path("/ssd/tnk_phase3/Integrated_dataset/gdtai_v4_2_tnk_reintegration/nk_boundary_partitions.npz")
BOUNDARY_UMAP = Path("/ssd/tnk_phase3/Integrated_dataset/gdtai_v4_2_tnk_reintegration/nk_boundary_umap_sample.npz")
R_SCRIPT = ROOT / "workflows/gdtai/run_boundary_pseudobulk_deg.R"
OUTPUT_NAME = "gdtai_v4_2_boundary_pseudobulk_deg"
TABLE_DIR = ROOT / f"Integrated_dataset/tables/gdT_prediction/{OUTPUT_NAME}"
FIGURE_DIR = ROOT / f"Integrated_dataset/figures/gdT_prediction/{OUTPUT_NAME}"
LOG_DIR = ROOT / f"Integrated_dataset/logs/gdT_prediction/{OUTPUT_NAME}"
REPORT_DIR = ROOT / f"gdT_prediction/{OUTPUT_NAME}"
ASSET_DIR = REPORT_DIR / "assets"
REVIEW_PARTITION = "boundary_r0.4_s29"
MIN_CELLS_PER_GROUP = 20
BLOCK_ROWS = 4096
FIGURE_DPI = 260
TCR_COMPLEX_GENES = {"CD3D", "CD3E", "CD3G", "CD247", "TRAT1"}
TCR_LOCUS_PATTERN = re.compile(r"^TR[ABDG][CVDJ]")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse-pseudobulk", action="store_true")
    return parser.parse_args()


def decode(values: np.ndarray) -> np.ndarray:
    if values.dtype.kind not in {"O", "S", "U"}:
        return values
    return np.asarray([value.decode() if isinstance(value, bytes) else str(value) for value in values], dtype=object)


def obs_values(obs: h5py.Group, name: str) -> np.ndarray:
    node = obs[name]
    if isinstance(node, h5py.Group) and {"categories", "codes"}.issubset(node):
        category_node = node["categories"]
        if isinstance(category_node, h5py.Group) and {"values", "mask"}.issubset(category_node):
            categories = decode(category_node["values"][:])
            categories[np.asarray(category_node["mask"][:], dtype=bool)] = ""
        else:
            categories = decode(category_node[:])
        codes = np.asarray(node["codes"][:], dtype=np.int64)
        result = np.full(len(codes), "", dtype=object)
        valid = (codes >= 0) & (codes < len(categories))
        result[valid] = categories[codes[valid]]
        return result
    if isinstance(node, h5py.Group) and {"values", "mask"}.issubset(node):
        result = decode(node["values"][:])
        result[np.asarray(node["mask"][:], dtype=bool)] = ""
        return result
    return decode(node[:])


def axis_names(handle: h5py.File, axis: str) -> np.ndarray:
    key = handle[axis].attrs["_index"]
    if isinstance(key, bytes):
        key = key.decode()
    return decode(handle[axis][str(key)][:])


def invalid_text(values: np.ndarray) -> np.ndarray:
    lowered = np.char.lower(np.asarray(values, dtype=str))
    return np.isin(lowered, ["", "na", "nan", "none", "unknown", "unresolved", "unassigned"])


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty_like(ranked)
    result[order] = np.minimum(ranked, 1.0)
    return result


def map_staged_rows_to_atlas(staged_rows: np.ndarray) -> np.ndarray:
    with h5py.File(STAGED, "r") as staged, h5py.File(CORRECTED, "r") as atlas:
        staged_source = obs_values(staged["obs"], "source_gse_id")
        staged_original = obs_values(staged["obs"], "original_cell_id")
        staged_cohort = obs_values(staged["obs"], "input_cohort_id")
        atlas_source = obs_values(atlas["obs"], "source_gse_id")
        atlas_original = obs_values(atlas["obs"], "original_cell_id")
        atlas_source_obs = obs_values(atlas["obs"], "source_obs_name")
    mapped = np.full(len(staged_rows), -1, dtype=np.int64)
    selected_source = staged_source[staged_rows]
    selected_original = staged_original[staged_rows]
    selected_cohort = staged_cohort[staged_rows]
    recovered = selected_cohort == "current_atlas"
    for source in sorted(set(selected_source[recovered])):
        local = np.flatnonzero(recovered & (selected_source == source))
        atlas_positions = np.flatnonzero(atlas_source == source)
        lookup = dict(zip(atlas_source_obs[atlas_positions], atlas_positions, strict=False))
        mapped[local] = np.fromiter((lookup.get(cell, -1) for cell in selected_original[local]), dtype=np.int64)
    for source in sorted(set(selected_source)):
        local = np.flatnonzero((selected_source == source) & ~recovered)
        if len(local) == 0:
            continue
        atlas_positions = np.flatnonzero(atlas_source == source)
        lookup = dict(zip(atlas_original[atlas_positions], atlas_positions, strict=False))
        mapped[local] = np.fromiter((lookup.get(cell, -1) for cell in selected_original[local]), dtype=np.int64)
    if np.any(mapped < 0) or len(np.unique(mapped)) != len(mapped):
        raise RuntimeError("Boundary staged-to-atlas mapping is incomplete or non-unique")
    return mapped


def cluster_enrichment(labels: np.ndarray, productive_ab: np.ndarray) -> pd.DataFrame:
    total = len(labels)
    total_positive = int(productive_ab.sum())
    rows = []
    for cluster in sorted(set(labels)):
        mask = labels == cluster
        n_cells = int(mask.sum())
        n_positive = int(productive_ab[mask].sum())
        odds, p_value = fisher_exact(
            [[n_positive, n_cells - n_positive], [total_positive - n_positive, (total - n_cells) - (total_positive - n_positive)]]
        )
        rows.append({
            "boundary_cluster": int(cluster), "n_cells": n_cells,
            "n_productive_TRA_or_TRB": n_positive,
            "fraction_productive_TRA_or_TRB": n_positive / n_cells,
            "odds_ratio_vs_rest": float(odds), "fisher_p_value": float(p_value),
        })
    result = pd.DataFrame(rows)
    result["fisher_fdr"] = bh_fdr(result["fisher_p_value"].to_numpy())
    result["cluster_group"] = np.where(
        (result["odds_ratio_vs_rest"] > 1) & (result["fisher_fdr"] < 0.05), "enriched",
        np.where((result["odds_ratio_vs_rest"] < 1) & (result["fisher_fdr"] < 0.05), "depleted", "excluded"),
    )
    return result


def prepare_boundary_metadata() -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, dict[str, object]]:
    with np.load(BOUNDARY) as archive:
        boundary_indices = np.asarray(archive["boundary_indices"], dtype=np.int64)
        names = np.asarray(archive["names"]).astype(str)
        labels = np.asarray(archive["labels"][:, np.flatnonzero(names == REVIEW_PARTITION)[0]], dtype=np.int32)
    atlas_rows = map_staged_rows_to_atlas(boundary_indices)
    with h5py.File(CORRECTED, "r") as handle:
        obs = handle["obs"]
        any_ab_all = np.asarray(obs["has_any_ab_tcr"][:], dtype=bool)
        source_all = obs_values(obs, "source_gse_id")
        display_all = obs_values(obs, "source_accession_harmonized_v2")
        sample_all = obs_values(obs, "sample_id_harmonized_v2")
        sample_fallback_all = obs_values(obs, "sample_id")
        library_all = obs_values(obs, "library_id_harmonized_v2")
        donor_all = obs_values(obs, "donor_id_harmonized_v2")
    productive_ab = any_ab_all[atlas_rows]
    enrichment = cluster_enrichment(labels, productive_ab)
    group_map = enrichment.set_index("boundary_cluster")["cluster_group"].to_dict()
    cluster_group = np.asarray([group_map[int(cluster)] for cluster in labels], dtype=object)
    source = np.asarray(source_all[atlas_rows], dtype=object)
    display = np.asarray(display_all[atlas_rows], dtype=object)
    display[invalid_text(display)] = source[invalid_text(display)]
    sample = np.asarray(sample_all[atlas_rows], dtype=object)
    fallback = np.asarray(sample_fallback_all[atlas_rows], dtype=object)
    library = np.asarray(library_all[atlas_rows], dtype=object)
    donor = np.asarray(donor_all[atlas_rows], dtype=object)
    bad = invalid_text(sample); sample[bad] = fallback[bad]
    bad = invalid_text(sample); sample[bad] = library[bad]
    bad = invalid_text(sample); sample[bad] = donor[bad]
    if np.any(invalid_text(sample)):
        raise RuntimeError("Some boundary cells lack a usable biological sample key")
    pair_key = np.char.add(np.char.add(source.astype(str), "::"), sample.astype(str))
    frame = pd.DataFrame({
        "boundary_local_row": np.arange(len(labels), dtype=np.int64), "atlas_row": atlas_rows,
        "boundary_cluster": labels, "cluster_group": cluster_group,
        "productive_TRA_or_TRB": productive_ab, "source_gse_id": source,
        "dataset_name": display, "sample_id_for_pseudobulk": sample, "pair_key": pair_key,
    })
    counts = frame.loc[frame["cluster_group"].isin(["enriched", "depleted"])].groupby(
        ["pair_key", "cluster_group"], observed=True
    ).size().unstack(fill_value=0)
    eligible = counts.index[(counts.get("enriched", 0) >= MIN_CELLS_PER_GROUP) & (counts.get("depleted", 0) >= MIN_CELLS_PER_GROUP)]
    selected = frame.loc[frame["pair_key"].isin(eligible) & frame["cluster_group"].isin(["enriched", "depleted"])].copy()
    pair_lookup = selected.drop_duplicates("pair_key").set_index("pair_key")[["source_gse_id", "dataset_name", "sample_id_for_pseudobulk"]]
    pair_rows = []
    for pair in sorted(eligible):
        for group in ("enriched", "depleted"):
            subset = selected.loc[(selected["pair_key"] == pair) & (selected["cluster_group"] == group)]
            pair_rows.append({
                "pseudobulk_id": f"{pair}::{group}", "pair_key": pair, "cluster_group": group,
                "source_gse_id": pair_lookup.loc[pair, "source_gse_id"],
                "dataset_name": pair_lookup.loc[pair, "dataset_name"],
                "sample_id_for_pseudobulk": pair_lookup.loc[pair, "sample_id_for_pseudobulk"],
                "n_cells": len(subset),
            })
    pseudobulk_meta = pd.DataFrame(pair_rows)
    pb_lookup = pseudobulk_meta.reset_index().set_index(["pair_key", "cluster_group"])["index"].to_dict()
    selected["pseudobulk_index"] = [pb_lookup[(pair, group)] for pair, group in zip(selected["pair_key"], selected["cluster_group"], strict=True)]
    summary = {
        "n_boundary_cells": len(frame), "n_corrected_productive_ab": int(productive_ab.sum()),
        "boundary_productive_ab_fraction": float(productive_ab.mean()),
        "enriched_clusters": enrichment.loc[enrichment["cluster_group"] == "enriched", "boundary_cluster"].astype(int).tolist(),
        "depleted_clusters": enrichment.loc[enrichment["cluster_group"] == "depleted", "boundary_cluster"].astype(int).tolist(),
        "n_eligible_pairs": len(eligible), "n_eligible_datasets": int(pseudobulk_meta["dataset_name"].nunique()),
        "n_selected_cells": len(selected), "minimum_cells_per_group": MIN_CELLS_PER_GROUP,
    }
    return selected, pseudobulk_meta, atlas_rows, {"summary": summary, "enrichment": enrichment, "all_boundary": frame}


def is_tcr_related(gene: str) -> tuple[bool, str]:
    if TCR_LOCUS_PATTERN.match(gene):
        return True, "TRA/TRB/TRD/TRG receptor locus"
    if gene in TCR_COMPLEX_GENES:
        return True, "TCR-CD3 receptor complex"
    return False, ""


def aggregate_pseudobulks(selected: pd.DataFrame, meta: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    before = {"size": CORRECTED.stat().st_size, "mtime_ns": CORRECTED.stat().st_mtime_ns}
    with h5py.File(CORRECTED, "r", rdcc_nbytes=1024 * 1024 * 1024, rdcc_nslots=2_000_003) as handle:
        genes = axis_names(handle, "var").astype(str)
        shape = tuple(int(value) for value in handle["X"].attrs["shape"])
        indptr = np.asarray(handle["X/indptr"][:], dtype=np.int64)
        positions = selected["atlas_row"].to_numpy(dtype=np.int64)
        pb_ids = selected["pseudobulk_index"].to_numpy(dtype=np.int64)
        order = np.argsort(positions)
        positions = positions[order]; pb_ids = pb_ids[order]
        if len(np.unique(positions)) != len(positions):
            raise RuntimeError("Selected atlas rows are not unique")
        counts = np.zeros((len(meta), shape[1]), dtype=np.float64)
        block_ids = np.unique(positions // BLOCK_ROWS)
        started = time.time()
        for block_number, block_id in enumerate(block_ids, start=1):
            row_start = int(block_id * BLOCK_ROWS); row_end = min(shape[0], row_start + BLOCK_ROWS)
            left = np.searchsorted(positions, row_start); right = np.searchsorted(positions, row_end)
            local_rows = positions[left:right] - row_start
            local_pb = pb_ids[left:right]
            data_start, data_end = int(indptr[row_start]), int(indptr[row_end])
            block = sparse.csr_matrix(
                (handle["X/data"][data_start:data_end], handle["X/indices"][data_start:data_end], indptr[row_start:row_end + 1] - data_start),
                shape=(row_end - row_start, shape[1]),
            )
            selected_block = block[local_rows]
            assignment = sparse.csr_matrix(
                (np.ones(len(local_rows), dtype=np.float64), (local_pb, np.arange(len(local_rows)))),
                shape=(len(meta), len(local_rows)),
            )
            aggregated = (assignment @ selected_block).tocoo()
            np.add.at(counts, (aggregated.row, aggregated.col), aggregated.data)
            if block_number % 100 == 0 or block_number == len(block_ids):
                print(f"aggregate blocks {block_number}/{len(block_ids)} elapsed={time.time()-started:.1f}s", flush=True)
    after = {"size": CORRECTED.stat().st_size, "mtime_ns": CORRECTED.stat().st_mtime_ns}
    if before != after:
        raise RuntimeError("Corrected atlas size/mtime changed during aggregation")
    expected_cells = selected.groupby("pseudobulk_index").size().reindex(range(len(meta)), fill_value=0).to_numpy()
    if not np.array_equal(expected_cells, meta["n_cells"].to_numpy()):
        raise RuntimeError("Pseudobulk cell-count audit failed")
    if not np.allclose(counts, np.rint(counts), atol=1e-5):
        raise RuntimeError("Raw pseudobulk counts are not integer-like")
    excluded_rows = []
    keep = np.ones(len(genes), dtype=bool)
    for index, gene in enumerate(genes):
        exclude, reason = is_tcr_related(str(gene))
        if exclude:
            keep[index] = False
            excluded_rows.append({"gene": gene, "reason": reason})
    filtered = np.rint(counts[:, keep]).astype(np.int64)
    if np.any(filtered.sum(axis=1) <= 0):
        raise RuntimeError("A pseudobulk has zero non-TCR library size")
    return filtered, genes[keep], pd.DataFrame(excluded_rows)


def write_r_inputs(counts: np.ndarray, genes: np.ndarray, meta: pd.DataFrame) -> None:
    matrix_path = TABLE_DIR / "pseudobulk_counts.mtx.gz"
    with gzip.open(matrix_path, "wb", compresslevel=6) as handle:
        mmwrite(handle, sparse.csr_matrix(counts), field="integer", symmetry="general")
    (TABLE_DIR / "genes.txt").write_text("\n".join(genes.astype(str)) + "\n", encoding="utf-8")
    meta.to_csv(TABLE_DIR / "pseudobulk_metadata.csv", index=False)


def save_figure(fig: plt.Figure, name: str) -> Path:
    path = FIGURE_DIR / name
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def make_figures(enrichment: pd.DataFrame, meta: pd.DataFrame, deg: pd.DataFrame) -> list[Path]:
    sns.set_theme(style="whitegrid", font_scale=0.9)
    figures = []
    colors = {"enriched": "#B23A48", "depleted": "#287271", "excluded": "#9CA3AF"}
    frame = enrichment.sort_values("boundary_cluster")
    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    ax.bar(frame["boundary_cluster"].astype(str), frame["fraction_productive_TRA_or_TRB"], color=[colors[value] for value in frame["cluster_group"]])
    overall = frame["n_productive_TRA_or_TRB"].sum() / frame["n_cells"].sum()
    ax.axhline(overall, color="#111827", linestyle="--", linewidth=1.2, label=f"Boundary mean {overall:.2%}")
    ax.set(xlabel="Frozen boundary subcluster", ylabel="Fraction with corrected productive TRA or TRB", title="Productive alpha-beta TCR enrichment defines comparison groups")
    ax.legend(frameon=False); ax.spines[["top", "right"]].set_visible(False)
    figures.append(save_figure(fig, "cluster_productive_ab_enrichment.png"))

    pairs = meta.drop_duplicates("pair_key")["dataset_name"].value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(9, 6.5), constrained_layout=True)
    ax.barh(pairs.index, pairs.values, color="#3A6B7E")
    ax.set(xlabel="Paired biological samples", title="Pseudobulk pairs retained by dataset")
    ax.spines[["top", "right"]].set_visible(False)
    figures.append(save_figure(fig, "paired_samples_by_dataset.png"))

    plot = deg.copy()
    plot["minus_log10_fdr"] = -np.log10(np.maximum(plot["paired_fdr"], 1e-300))
    plot["status"] = np.where(plot["robust_cross_dataset"], "Robust cross-dataset", np.where(plot["paired_significant"], "Paired only", "Not significant"))
    fig, ax = plt.subplots(figsize=(9, 6.5), constrained_layout=True)
    for status, color, size, alpha in [("Not significant", "#B8BEC4", 8, 0.35), ("Paired only", "#D99B45", 11, 0.65), ("Robust cross-dataset", "#A82035", 15, 0.8)]:
        take = plot["status"] == status
        ax.scatter(plot.loc[take, "mean_log2FC_enriched_vs_depleted"], plot.loc[take, "minus_log10_fdr"], s=size, alpha=alpha, linewidths=0, color=color, label=status, rasterized=True)
    label_rows = plot.loc[plot["robust_cross_dataset"]].sort_values("paired_fdr").head(18)
    if len(label_rows) < 10:
        label_rows = plot.loc[plot["paired_significant"]].sort_values("paired_fdr").head(18)
    for row in label_rows.itertuples():
        ax.text(row.mean_log2FC_enriched_vs_depleted, row.minus_log10_fdr, row.gene, fontsize=7)
    ax.axvline(-0.25, color="#6B7280", linestyle="--"); ax.axvline(0.25, color="#6B7280", linestyle="--")
    ax.set(xlabel="Mean paired log2FC: enriched - depleted", ylabel="-log10 paired FDR", title="Non-TCR pseudobulk differential expression")
    ax.legend(frameon=False); ax.spines[["top", "right"]].set_visible(False)
    figures.append(save_figure(fig, "deg_volcano.png"))

    source_effect = pd.read_csv(TABLE_DIR / "source_mean_log2fc.csv.gz").set_index("gene")
    top = deg.loc[deg["robust_cross_dataset"]].sort_values("paired_fdr").head(30)["gene"].tolist()
    if len(top) < 15:
        top = deg.loc[deg["paired_significant"]].sort_values("paired_fdr").head(30)["gene"].tolist()
    heat = source_effect.loc[top]
    vmax = max(0.5, float(np.nanquantile(np.abs(heat.to_numpy()), 0.98)))
    fig, ax = plt.subplots(figsize=(13, max(5.5, 0.28 * len(top))), constrained_layout=True)
    sns.heatmap(heat, cmap="vlag", norm=TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax), ax=ax, cbar_kws={"label": "Mean paired log2FC"})
    ax.set(xlabel="Dataset (>=5 paired samples)", ylabel="Gene", title="Cross-dataset direction of top non-TCR DEGs")
    figures.append(save_figure(fig, "top_deg_dataset_effect_heatmap.png"))

    with np.load(BOUNDARY_UMAP) as archive:
        xy = np.asarray(archive["X_umap"], dtype=np.float32)
        clusters = np.asarray(archive["boundary_cluster"], dtype=int)
    group_map = frame.set_index("boundary_cluster")["cluster_group"].to_dict()
    groups = np.asarray([group_map[int(value)] for value in clusters], dtype=object)
    fig, ax = plt.subplots(figsize=(8.2, 6.5), constrained_layout=True)
    for group in ["depleted", "enriched"]:
        take = groups == group
        ax.scatter(xy[take, 0], xy[take, 1], s=0.45, linewidths=0, alpha=0.55, color=colors[group], label=f"{group} clusters", rasterized=True)
    ax.set_title("Frozen NK-boundary partition used for pseudobulk DEG"); ax.set_xticks([]); ax.set_yticks([])
    ax.legend(frameon=False, markerscale=7)
    for spine in ax.spines.values(): spine.set_visible(False)
    figures.append(save_figure(fig, "boundary_umap_comparison_groups.png"))
    return figures


def table_html(frame: pd.DataFrame, max_rows: int = 30) -> str:
    shown = frame.head(max_rows).copy()
    for column in shown.select_dtypes(include="float").columns:
        shown[column] = shown[column].map(lambda value: f"{value:.4g}" if pd.notna(value) else "")
    return f"<div class='table-wrap'>{shown.to_html(index=False, border=0, escape=True)}</div>"


def render_report(summary: dict[str, object], enrichment: pd.DataFrame, meta: pd.DataFrame, deg: pd.DataFrame, excluded: pd.DataFrame, figures: list[Path]) -> tuple[Path, Path]:
    for directory in (REPORT_DIR, ASSET_DIR): directory.mkdir(parents=True, exist_ok=True)
    for source in figures:
        target = ASSET_DIR / source.name
        if target.exists() or target.is_symlink(): target.unlink()
        os.symlink(os.path.relpath(source, ASSET_DIR), target)
    top_up = deg.loc[deg["paired_significant"] & (deg["mean_log2FC_enriched_vs_depleted"] > 0)].sort_values("paired_fdr")
    top_down = deg.loc[deg["paired_significant"] & (deg["mean_log2FC_enriched_vs_depleted"] < 0)].sort_values("paired_fdr")
    page = f"""<!doctype html><html><head><meta charset='utf-8'><title>NK-boundary pseudobulk DEG</title><style>
body{{font-family:Arial,sans-serif;margin:0;background:#eef2f3;color:#20262c;line-height:1.5}}main{{max-width:1180px;margin:18px auto}}header,section{{background:white;border:1px solid #d7dde1;padding:26px 30px;margin:14px 0;border-radius:6px}}header{{background:#173640;color:white}}h1{{font-size:36px;margin:0 0 8px}}h2{{margin:0 0 10px}}.metrics{{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin-top:20px}}.metric{{border-top:2px solid #79aebb;padding-top:8px}}.metric b{{display:block;font-size:23px}}.note{{background:#f0f5f6;border-left:4px solid #357180;padding:12px}}.warning{{background:#fff4e5;border-left:4px solid #c77b1a;padding:12px}}.grid{{display:grid;grid-template-columns:1fr 1fr;gap:18px}}figure{{margin:10px 0;break-inside:avoid}}img{{width:100%}}figcaption{{font-size:12px;color:#66717a}}.table-wrap{{overflow:auto;max-height:520px;border:1px solid #d7dde1}}table{{border-collapse:collapse;width:100%;font-size:11px}}th,td{{padding:6px;border-bottom:1px solid #e4e8eb;text-align:left;white-space:nowrap}}th{{background:#eef3f5}}code{{word-break:break-all}}@media print{{body{{background:white}}main{{max-width:none;margin:0}}header,section{{border-radius:0}}h2,h3{{break-after:avoid-page}}.table-wrap{{overflow:visible;max-height:none}}table{{font-size:8px}}.page-break{{break-before:page}}}}
</style></head><body><main><header><h1>NK-Boundary Pseudobulk Differential Expression</h1><p>Corrected productive TRA/TRB enrichment defines frozen transcriptomic cluster groups; paired biological-sample pseudobulks provide the statistical unit.</p><div class='metrics'><div class='metric'><b>{summary['n_boundary_cells']:,}</b><span>boundary cells</span></div><div class='metric'><b>{summary['n_eligible_pairs']:,}</b><span>paired samples</span></div><div class='metric'><b>{summary['n_eligible_datasets']:,}</b><span>datasets</span></div><div class='metric'><b>{int(deg['robust_cross_dataset'].sum()):,}</b><span>robust DEGs</span></div></div></header>
<section><h2>Design and interpretation</h2><p class='note'>Frozen boundary subclusters were tested against the rest of the 475,953-cell boundary using corrected productive TRA-or-TRB metadata and two-sided Fisher tests with BH correction. Enriched clusters are {summary['enriched_clusters']}; depleted clusters are {summary['depleted_clusters']}. No cluster is literally receptor-free.</p><p>Raw counts were summed separately for enriched and depleted cluster sets within each biological sample. A pair was retained only when both groups contained at least {summary['minimum_cells_per_group']} cells. TMM-normalized log2-CPM differences were tested across paired samples with limma empirical Bayes. A dataset-macro analysis first averaged paired effects within datasets having at least five pairs.</p><p class='warning'><b>TCR exclusion:</b> {len(excluded):,} genes were removed before normalization and testing: all TRA/TRB/TRD/TRG receptor-locus genes plus CD3D, CD3E, CD3G, CD247, and TRAT1. Downstream signaling genes such as LCK, ZAP70, LAT, and LCP2 were retained.</p></section>
<section><h2>Cluster definition</h2>{table_html(enrichment, 20)}<div class='grid'><figure><img src='assets/cluster_productive_ab_enrichment.png'></figure><figure><img src='assets/boundary_umap_comparison_groups.png'></figure></div></section>
<section><h2>Pseudobulk coverage</h2><figure><img src='assets/paired_samples_by_dataset.png'></figure>{table_html(pd.read_csv(TABLE_DIR/'paired_samples_by_dataset.csv').sort_values('n_paired_samples',ascending=False), 40)}</section>
<section class='page-break'><h2>Differential expression</h2><div class='grid'><figure><img src='assets/deg_volcano.png'></figure><figure><img src='assets/top_deg_dataset_effect_heatmap.png'></figure></div><h3>Top enriched-cluster genes</h3>{table_html(top_up[['gene','mean_log2FC_enriched_vs_depleted','paired_fdr','dataset_macro_mean_log2FC','dataset_macro_fdr','source_sign_consistency','robust_cross_dataset']],30)}<h3>Top depleted-cluster genes</h3>{table_html(top_down[['gene','mean_log2FC_enriched_vs_depleted','paired_fdr','dataset_macro_mean_log2FC','dataset_macro_fdr','source_sign_consistency','robust_cross_dataset']],30)}</section>
<section><h2>Output and claim boundary</h2><p>Full results: <code>{TABLE_DIR/'deg_results.csv.gz'}</code></p><p>These DEGs describe transcriptomic differences between frozen boundary cluster sets associated with different corrected productive-alpha-beta TCR frequencies. They do not prove that depleted clusters are pure NK cells, and they are not new training labels.</p><p>The corrected atlas, boundary partition, cell membership, UMAP, and receptor metadata were read only.</p></section></main></body></html>"""
    index = REPORT_DIR / "index.html"; index.write_text(page, encoding="utf-8")
    pdf = REPORT_DIR / "boundary_pseudobulk_deg_report.pdf"
    subprocess.run(["google-chrome", "--headless", "--no-sandbox", "--disable-gpu", "--print-to-pdf-no-header", f"--print-to-pdf={pdf}", str(index)], check=True)
    return index, pdf


def main() -> None:
    args = parse_args()
    for directory in (TABLE_DIR, FIGURE_DIR, LOG_DIR, REPORT_DIR, ASSET_DIR): directory.mkdir(parents=True, exist_ok=True)
    input_signatures = {str(path): {"size": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns} for path in (CORRECTED, STAGED, BOUNDARY, BOUNDARY_UMAP)}
    selected, meta, _atlas_rows, payload = prepare_boundary_metadata()
    enrichment = payload["enrichment"]
    summary = payload["summary"]
    enrichment.to_csv(TABLE_DIR / "cluster_productive_ab_enrichment.csv", index=False)
    selected.drop(columns="pseudobulk_index").to_csv(TABLE_DIR / "selected_boundary_cells.csv.gz", index=False, compression="gzip")
    meta.to_csv(TABLE_DIR / "pseudobulk_metadata.csv", index=False)
    matrix_path = TABLE_DIR / "pseudobulk_counts.mtx.gz"
    genes_path = TABLE_DIR / "genes.txt"
    excluded_path = TABLE_DIR / "excluded_tcr_genes.csv"
    if args.reuse_pseudobulk and matrix_path.exists() and genes_path.exists() and excluded_path.exists():
        excluded = pd.read_csv(excluded_path)
    else:
        counts, genes, excluded = aggregate_pseudobulks(selected, meta)
        write_r_inputs(counts, genes, meta)
        excluded.to_csv(excluded_path, index=False)
    subprocess.run(["Rscript", str(R_SCRIPT), str(TABLE_DIR), str(TABLE_DIR)], check=True)
    deg = pd.read_csv(TABLE_DIR / "deg_results.csv.gz")
    figures = make_figures(enrichment, meta, deg)
    index, pdf = render_report(summary, enrichment, meta, deg, excluded, figures)
    final_signatures = {str(path): {"size": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns} for path in (CORRECTED, STAGED, BOUNDARY, BOUNDARY_UMAP)}
    if input_signatures != final_signatures:
        raise RuntimeError("A read-only input changed during pseudobulk analysis")
    run_summary = {
        "status": "PASS_BOUNDARY_PSEUDOBULK_DEG", **summary,
        "n_tcr_genes_excluded": len(excluded), "n_genes_tested": len(deg),
        "n_paired_significant": int(deg["paired_significant"].sum()),
        "n_robust_cross_dataset": int(deg["robust_cross_dataset"].sum()),
        "input_signatures": input_signatures, "inputs_unchanged": True,
        "html": str(index), "pdf": str(pdf), "figures": [str(path) for path in figures],
    }
    (LOG_DIR / "run_summary.json").write_text(json.dumps(run_summary, indent=2) + "\n", encoding="utf-8")
    (LOG_DIR / "summary.md").write_text(
        "# NK-boundary pseudobulk DEG\n\n" + "\n".join(f"- {key}: `{value}`" for key, value in run_summary.items() if key not in {"input_signatures", "figures"}) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(run_summary, indent=2))


if __name__ == "__main__":
    main()
