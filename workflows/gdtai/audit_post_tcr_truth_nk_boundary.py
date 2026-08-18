#!/usr/bin/env python3
"""Regenerate truth/control and frozen NK-boundary audits after TCR repair."""

from __future__ import annotations

import hashlib
import html
import json
import shutil
import subprocess
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
CORRECTED = Path("/ssd/tnk_phase3/Integrated_dataset/full_atlas/tcr_corrected/integrated_full_atlas.h5ad")
LEGACY = Path("/ssd/tnk_phase3/Integrated_dataset/full_atlas/metadata_corrected/integrated_full_atlas.h5ad")
STAGED = Path("/ssd/tnk_phase3/Integrated_dataset/gdtai_v4_2_tnk_reintegration/tnk_refined_hvg_counts.h5ad")
BOUNDARY = Path("/ssd/tnk_phase3/Integrated_dataset/gdtai_v4_2_tnk_reintegration/nk_boundary_partitions.npz")
BOUNDARY_UMAP = Path("/ssd/tnk_phase3/Integrated_dataset/gdtai_v4_2_tnk_reintegration/nk_boundary_umap_sample.npz")
TABLE_DIR = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_post_tcr_audit"
FIGURE_DIR = ROOT / "Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_post_tcr_audit"
LOG_DIR = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_post_tcr_audit"
REPORT_DIR = ROOT / "gdT_prediction/gdtai_v4_2_post_tcr_audit"
ASSET_DIR = REPORT_DIR / "assets"
CORRECTED_SHA = "d32c9d2bdb955b12e1eafbed8322f8cb965cf3a225191e612b53f3d3783480d5"


def sha256_file(path: Path, chunk: int = 64 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def decode(values: np.ndarray) -> np.ndarray:
    return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in values], dtype=object)


def obs_values(obs: h5py.Group, name: str) -> np.ndarray:
    node = obs[name]
    if isinstance(node, h5py.Group) and {"categories", "codes"}.issubset(node):
        category_node = node["categories"]
        if isinstance(category_node, h5py.Group) and {"values", "mask"}.issubset(category_node):
            categories = decode(category_node["values"][:])
            categories[np.asarray(category_node["mask"][:], dtype=bool)] = ""
        else:
            categories = decode(category_node[:])
        codes = node["codes"][:]
        result = np.full(len(codes), "", dtype=object)
        present = codes >= 0
        result[present] = categories[codes[present]]
        return result
    if isinstance(node, h5py.Group) and {"values", "mask"}.issubset(node):
        result = decode(node["values"][:])
        result[np.asarray(node["mask"][:], dtype=bool)] = ""
        return result
    values = node[:]
    return decode(values) if values.dtype.kind in {"O", "S", "U"} else np.asarray(values)


def read_positions(dataset: h5py.Dataset, positions: np.ndarray) -> np.ndarray:
    order = np.argsort(positions)
    sorted_values = dataset[positions[order]]
    result = np.empty_like(sorted_values)
    result[order] = sorted_values
    return result


def map_staged_rows_to_atlas(staged_rows: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
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
        lookup = {
            cell_id: pos
            for cell_id, pos in zip(atlas_source_obs[atlas_positions], atlas_positions)
        }
        mapped[local] = np.fromiter(
            (lookup.get(cell_id, -1) for cell_id in selected_original[local]),
            dtype=np.int64,
        )
    for source in sorted(set(selected_source)):
        local = np.flatnonzero((selected_source == source) & ~recovered)
        if len(local) == 0:
            continue
        atlas_positions = np.flatnonzero(atlas_source == source)
        lookup = {cell_id: pos for cell_id, pos in zip(atlas_original[atlas_positions], atlas_positions)}
        mapped[local] = np.fromiter((lookup.get(cell_id, -1) for cell_id in selected_original[local]), dtype=np.int64)
    if np.any(mapped < 0) or len(set(mapped)) != len(mapped):
        raise RuntimeError("Staged-to-atlas source/cell mapping is incomplete or non-unique")
    return mapped, selected_source


def save(fig: plt.Figure, name: str) -> None:
    path = FIGURE_DIR / name
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    shutil.copy2(path, ASSET_DIR / name)


def truth_tables() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    source_rows = []
    tissue_rows = []
    totals = {"gdT_gold": 0, "abT_gold": 0, "strict_gdT_receptor": 0, "dual_tcr": 0, "conflict": 0}
    with h5py.File(CORRECTED, "r") as handle:
        obs = handle["obs"]
        sources = obs_values(obs, "source_gse_id")
        tissues = obs_values(obs, "tissue_harmonized_v2")
        sorted_gdt = np.char.lower(obs_values(obs, "sorted_gdt").astype(str)) == "true"
        any_ab = np.asarray(obs["has_any_ab_tcr"][:], dtype=bool)
        any_gd = np.asarray(obs["has_any_gd_tcr"][:], dtype=bool)
        paired_ab = np.asarray(obs["has_TRA_TRB_paired"][:], dtype=bool)
        paired_gd = np.asarray(obs["has_TRG_TRD_paired"][:], dtype=bool)
    strict_gd = paired_gd & ~any_ab
    gdt_gold = sorted_gdt | strict_gd
    conflict = gdt_gold & paired_ab & ~any_gd
    abt_gold = paired_ab & ~any_gd & ~gdt_gold
    dual = any_ab & any_gd
    for key, mask in (("gdT_gold", gdt_gold), ("abT_gold", abt_gold), ("strict_gdT_receptor", strict_gd), ("dual_tcr", dual), ("conflict", conflict)):
        totals[key] = int(mask.sum())
    for source in sorted(set(sources)):
        mask = sources == source
        source_rows.append({"source_gse_id": source, "n_cells": int(mask.sum()), **{key: int((mask & value).sum()) for key, value in (("gdT_gold", gdt_gold), ("abT_gold", abt_gold), ("strict_gdT_receptor", strict_gd), ("dual_tcr", dual), ("conflict", conflict))}})
    for tissue in sorted(set(tissues)):
        mask = tissues == tissue
        tissue_rows.append({"tissue_harmonized_v2": tissue, "n_cells": int(mask.sum()), "gdT_gold": int((mask & gdt_gold).sum()), "abT_gold": int((mask & abt_gold).sum())})
    return pd.DataFrame(source_rows), pd.DataFrame(tissue_rows), totals


def main() -> None:
    for path in (TABLE_DIR, FIGURE_DIR, LOG_DIR, REPORT_DIR, ASSET_DIR):
        path.mkdir(parents=True, exist_ok=True)
    if sha256_file(CORRECTED) != CORRECTED_SHA:
        raise RuntimeError("Corrected candidate checksum mismatch")
    with np.load(BOUNDARY) as archive:
        boundary_indices = np.asarray(archive["boundary_indices"], dtype=np.int64)
        boundary_labels = np.asarray(archive["labels"], dtype=np.int32)
        boundary_names = np.asarray(archive["names"]).astype(str)
    review_col = int(np.flatnonzero(boundary_names == "boundary_r0.4_s29")[0])
    review_cluster = boundary_labels[:, review_col]
    with h5py.File(STAGED, "r") as handle:
        primary_all = np.asarray(handle["obs"]["primary_nk_anchor"][:], dtype=bool)
    primary_rows = np.flatnonzero(primary_all)
    needed = np.unique(np.concatenate([boundary_indices, primary_rows]))
    mapped_needed, _ = map_staged_rows_to_atlas(needed)
    map_lookup = dict(zip(needed, mapped_needed))
    boundary_atlas = np.fromiter((map_lookup[row] for row in boundary_indices), dtype=np.int64)
    primary_atlas = np.fromiter((map_lookup[row] for row in primary_rows), dtype=np.int64)

    def receptor_at(path: Path, positions: np.ndarray) -> dict[str, np.ndarray]:
        with h5py.File(path, "r") as handle:
            obs = handle["obs"]
            return {
                name: np.asarray(obs[name][:], dtype=bool)[positions]
                for name in ("has_TRA", "has_TRB", "has_TRD", "has_any_ab_tcr", "has_any_gd_tcr", "has_TRA_TRB_paired")
            } | {"source": obs_values(obs, "source_gse_id")[positions]}

    old_primary = receptor_at(LEGACY, primary_atlas)
    new_primary = receptor_at(CORRECTED, primary_atlas)
    primary_source_rows = []
    for source in sorted(set(new_primary["source"])):
        mask = new_primary["source"] == source
        primary_source_rows.append({
            "source_gse_id": source,
            "n_primary_nk": int(mask.sum()),
            "legacy_any_ab": int((mask & old_primary["has_any_ab_tcr"]).sum()),
            "corrected_any_ab": int((mask & new_primary["has_any_ab_tcr"]).sum()),
            "corrected_any_gd": int((mask & new_primary["has_any_gd_tcr"]).sum()),
            "corrected_no_productive_tcr": int((mask & ~new_primary["has_any_ab_tcr"] & ~new_primary["has_any_gd_tcr"]).sum()),
        })
    primary_df = pd.DataFrame(primary_source_rows)
    primary_df.to_csv(TABLE_DIR / "primary_nk_anchor_tcr_conflicts.csv", index=False)

    old_boundary = receptor_at(LEGACY, boundary_atlas)
    new_boundary = receptor_at(CORRECTED, boundary_atlas)
    boundary_rows = []
    for cluster in sorted(set(review_cluster)):
        mask = review_cluster == cluster
        boundary_rows.append({
            "boundary_cluster": int(cluster),
            "n_cells": int(mask.sum()),
            "legacy_any_ab": int((mask & old_boundary["has_any_ab_tcr"]).sum()),
            "corrected_any_ab": int((mask & new_boundary["has_any_ab_tcr"]).sum()),
            "corrected_paired_ab": int((mask & new_boundary["has_TRA_TRB_paired"]).sum()),
            "corrected_any_gd": int((mask & new_boundary["has_any_gd_tcr"]).sum()),
        })
    boundary_df = pd.DataFrame(boundary_rows)
    boundary_df.to_csv(TABLE_DIR / "boundary_tcr_by_cluster.csv", index=False)
    source_df, tissue_df, truth_totals = truth_tables()
    source_df.to_csv(TABLE_DIR / "corrected_truth_by_source.csv", index=False)
    tissue_df.to_csv(TABLE_DIR / "corrected_truth_by_tissue.csv", index=False)

    plot_primary = primary_df.sort_values("corrected_any_ab")
    y = np.arange(len(plot_primary))
    fig, ax = plt.subplots(figsize=(8.8, 4.6))
    ax.barh(y + 0.18, plot_primary["legacy_any_ab"], height=0.34, color="#9CA3AF", label="Legacy")
    ax.barh(y - 0.18, plot_primary["corrected_any_ab"], height=0.34, color="#087E8B", label="Corrected")
    ax.set_yticks(y, plot_primary["source_gse_id"]); ax.set_xlabel("Primary NK anchors with productive TRA or TRB")
    ax.set_title("Validated TCR repair changes frozen NK-anchor conflicts"); ax.grid(axis="x", color="#E5E7EB"); ax.legend(frameon=False)
    ax.spines[["top", "right", "left"]].set_visible(False); save(fig, "primary_nk_tcr_conflicts.png")

    plot_boundary = boundary_df.sort_values("boundary_cluster")
    x = np.arange(len(plot_boundary))
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    ax.bar(x - 0.2, plot_boundary["legacy_any_ab"] / plot_boundary["n_cells"], width=0.4, color="#9CA3AF", label="Legacy")
    ax.bar(x + 0.2, plot_boundary["corrected_any_ab"] / plot_boundary["n_cells"], width=0.4, color="#087E8B", label="Corrected")
    ax.set_xticks(x, plot_boundary["boundary_cluster"]); ax.set_xlabel("Frozen boundary subcluster"); ax.set_ylabel("Fraction with productive TRA or TRB")
    ax.set_title("Same boundary and clustering, corrected receptor evidence"); ax.grid(axis="y", color="#E5E7EB"); ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False); save(fig, "boundary_tcr_conflicts_by_cluster.png")

    with np.load(BOUNDARY_UMAP) as archive:
        local = np.asarray(archive["sample_local_indices"], dtype=np.int64)
        xy = np.asarray(archive["X_umap"], dtype=np.float32)
        old_sample = np.asarray(archive["productive_tcr"], dtype=bool)
    corrected_sample = new_boundary["has_any_ab_tcr"][local]
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.0))
    for ax, values, title in ((axes[0], old_sample, "Legacy productive-TCR overlay"), (axes[1], corrected_sample, "Corrected productive TRA/TRB")):
        ax.scatter(xy[~values, 0], xy[~values, 1], s=0.35, c="#D5D8DC", rasterized=True)
        ax.scatter(xy[values, 0], xy[values, 1], s=0.55, c="#D1495B", rasterized=True)
        ax.set_title(title); ax.set_xticks([]); ax.set_yticks([]); ax.spines[:].set_visible(False)
    fig.suptitle("Frozen NK-boundary UMAP after receptor metadata repair", fontsize=14); fig.tight_layout(); save(fig, "boundary_umap_legacy_vs_corrected_tcr.png")

    corrected_primary_ab = int(new_primary["has_any_ab_tcr"].sum())
    corrected_primary_gd = int(new_primary["has_any_gd_tcr"].sum())
    corrected_primary_clean = int((~new_primary["has_any_ab_tcr"] & ~new_primary["has_any_gd_tcr"]).sum())
    summary = {
        "status": "PASS_POST_TCR_TRUTH_NK_AUDIT",
        "candidate_sha256": CORRECTED_SHA,
        "n_primary_nk_anchors": len(primary_rows),
        "n_primary_nk_legacy_any_ab": int(old_primary["has_any_ab_tcr"].sum()),
        "n_primary_nk_prior_unsafe_overlay_any_ab": 11_526,
        "n_primary_nk_corrected_any_ab": corrected_primary_ab,
        "n_primary_nk_corrected_any_gd": corrected_primary_gd,
        "n_primary_nk_corrected_no_productive_tcr": corrected_primary_clean,
        "n_boundary_cells": len(boundary_indices),
        "n_boundary_legacy_any_ab": int(old_boundary["has_any_ab_tcr"].sum()),
        "n_boundary_prior_unsafe_overlay_any_ab": 373_339,
        "n_boundary_corrected_any_ab": int(new_boundary["has_any_ab_tcr"].sum()),
        "n_boundary_corrected_paired_ab": int(new_boundary["has_TRA_TRB_paired"].sum()),
        "n_boundary_corrected_any_gd": int(new_boundary["has_any_gd_tcr"].sum()),
        "truth": truth_totals,
        "integration_or_clustering_rerun": False,
        "canonical_symlink_switched": False,
    }
    (LOG_DIR / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    rows = "".join(f"<tr><td>{html.escape(str(r.source_gse_id))}</td><td>{r.n_primary_nk:,}</td><td>{r.legacy_any_ab:,}</td><td>{r.corrected_any_ab:,}</td><td>{r.corrected_no_productive_tcr:,}</td></tr>" for r in primary_df.itertuples(index=False))
    report = f"""<!doctype html><html><head><meta charset='utf-8'><title>Post-TCR truth and NK audit</title><style>
    @page{{size:A4 landscape;margin:10mm}}body{{font:14px/1.45 Arial;color:#17202A;margin:0}}main{{max-width:1160px;margin:auto}}section{{padding:26px 34px;border-bottom:1px solid #D5DBDB;page-break-inside:avoid}}h1{{background:#123B45;color:white;padding:34px;margin:0}}h2{{color:#123B45}}.grid{{display:grid;grid-template-columns:repeat(4,1fr);gap:10px}}.card{{border-left:4px solid #18A999;background:#F2F7F7;padding:12px}}.big{{font-size:23px;font-weight:bold}}img{{width:100%;max-height:560px;object-fit:contain}}table{{border-collapse:collapse;width:100%;font-size:11px}}th,td{{border:1px solid #D5DBDB;padding:6px;text-align:right}}th{{background:#123B45;color:white}}th:first-child,td:first-child{{text-align:left}}.page{{page-break-before:always}}.note{{background:#EAF5F3;border-left:4px solid #18A999;padding:12px}}
    .primary-figure{max-height:350px}
    </style></head><body><main><h1>Post-repair TCR truth and NK-boundary audit</h1><section><div class='grid'>
    <div class='card'><div class='big'>{len(primary_rows):,}</div>frozen primary NK anchors</div><div class='card'><div class='big'>{corrected_primary_ab:,}</div>corrected AB conflicts</div><div class='card'><div class='big'>{corrected_primary_clean:,}</div>clean primary NK anchors</div><div class='card'><div class='big'>{truth_totals['gdT_gold']:,}</div>rebuilt gdT gold</div></div>
    <p class='note'>Cell membership, scVI, UMAP, boundary definition, and clustering are frozen. Only receptor evidence changed. The prior unsafe harmonized overlay marked 11,526 primary NK anchors and 373,339 boundary cells as alpha-beta TCR-positive; validated repair reduces these to {corrected_primary_ab:,} and {int(new_boundary['has_any_ab_tcr'].sum()):,}.</p></section>
    <section><h2>Primary NK-anchor audit</h2><img class='primary-figure' src='assets/primary_nk_tcr_conflicts.png'><table><thead><tr><th>Source</th><th>Primary NK</th><th>Legacy AB</th><th>Corrected AB</th><th>Corrected no TCR</th></tr></thead><tbody>{rows}</tbody></table></section>
    <section class='page'><h2>Frozen boundary by subcluster</h2><img src='assets/boundary_tcr_conflicts_by_cluster.png'></section>
    <section class='page'><h2>Boundary UMAP</h2><img src='assets/boundary_umap_legacy_vs_corrected_tcr.png'><p>The corrected overlay uses validated productive TRA/TRB calls. Lack of TRD remains uninformative for alpha-beta-only VDJ libraries.</p></section>
    <section class='page'><h2>Truth reconstruction</h2><p>Gold gdT includes sorted-gdT cells plus paired TRG/TRD cells without alpha-beta TCR. Gold abT requires paired TRA/TRB, no gamma-delta TCR, and no overlap with gdT gold. Silver positives are not used.</p><p><b>gdT gold:</b> {truth_totals['gdT_gold']:,}; <b>abT gold:</b> {truth_totals['abT_gold']:,}; <b>strict receptor gdT:</b> {truth_totals['strict_gdT_receptor']:,}; <b>dual TCR:</b> {truth_totals['dual_tcr']:,}; <b>rule conflicts:</b> {truth_totals['conflict']:,}.</p></section>
    <section><h2>Decision</h2><p>This audit determines whether the previous NK-reference exclusions were driven by faulty receptor joins. Only dual-annotation NK anchors with no corrected productive TCR are eligible for the next training review. No model was trained or promoted, and the canonical atlas was not switched.</p></section>
    </main></body></html>"""
    html_path = REPORT_DIR / "index.html"; html_path.write_text(report, encoding="utf-8")
    pdf_path = REPORT_DIR / "post_tcr_truth_nk_audit.pdf"
    subprocess.run(["google-chrome", "--headless", "--no-sandbox", "--disable-gpu", "--disable-crash-reporter", f"--print-to-pdf={pdf_path}", f"file://{html_path}"], check=True)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
