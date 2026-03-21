#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Phase 2 merged cleanup for the TNK integration workflow.

This script performs a stricter second-pass cleanup on the merged
`Integrated_dataset/TNK_candidates.h5ad` milestone. The implementation stays
conservative: it removes cells only when merged-context evidence supports an
obvious off-target identity or a severe low-quality pattern, then reapplies the
canonical gene filter of at least 500 expressing cells.

The filtered result is written to `Integrated_dataset/TNK_cleaned.h5ad`.
QC tables, figures, and a markdown summary are written under
`Integrated_dataset/`.
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp


# Config
INPUT_H5AD = Path("Integrated_dataset/TNK_candidates.h5ad")
OUTPUT_H5AD = Path("Integrated_dataset/TNK_cleaned.h5ad")
OUTPUT_ROOT = Path("Integrated_dataset")
TABLE_DIR = OUTPUT_ROOT / "tables"
FIGURE_DIR = OUTPUT_ROOT / "figures"
LOG_DIR = OUTPUT_ROOT / "logs"
TEMP_OUTPUT_H5AD = OUTPUT_ROOT / "TNK_cleaned.phase2_tmp.h5ad"

CHUNK_SIZE = 20_000
MIN_GENE_CELLS = 500
SAMPLE_EVERY_N = 35

PHASE2_SUMMARY_CSV = TABLE_DIR / "phase2_cleanup_summary.csv"
PHASE2_GSE_CSV = TABLE_DIR / "phase2_gse_before_after.csv"
PHASE2_REMOVED_CELLS_CSV = TABLE_DIR / "phase2_removed_cells.csv"
PHASE2_GENE_CSV = TABLE_DIR / "phase2_gene_detection_summary.csv"
PHASE2_SAMPLED_CELLS_CSV = TABLE_DIR / "phase2_sampled_cells_for_qc.csv"
PHASE2_QC_MD = LOG_DIR / "phase2_qc_summary.md"

PHASE2_GSE_PNG = FIGURE_DIR / "phase2_gse_cell_retention.png"
PHASE2_REMOVE_PNG = FIGURE_DIR / "phase2_removal_reasons.png"
PHASE2_SCORE_PNG = FIGURE_DIR / "phase2_offtarget_vs_tnk.png"
PHASE2_GENE_PNG = FIGURE_DIR / "phase2_gene_detection_distribution.png"

KEY_TNK_GENES = [
    "TRAC",
    "TRBC1",
    "TRBC2",
    "TRDC",
    "TRGC1",
    "TRGC2",
    "CD3D",
    "CD3E",
    "CD3G",
    "NKG7",
    "KLRD1",
]

REQUIRED_OBS_COLUMNS = [
    "source_gse_id",
    "phase1_selection_reason",
    "phase1_t_score",
    "phase1_gd_score",
    "phase1_nk_score",
    "phase1_contam_score",
    "phase1_t_hits",
    "phase1_gd_hits",
    "phase1_nk_hits",
    "phase1_nk_strong_hits",
]

MARKER_GROUPS = {
    "myeloid": ["LYZ", "FCER1G", "CTSS", "AIF1", "LST1", "SAT1", "CST3"],
    "bcell": ["MS4A1", "CD79A", "CD74", "HLA-DRA"],
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "erythroid": ["HBB", "HBA1", "HBA2"],
    "platelet": ["PF4", "PPBP"],
    "mitochondrial": [],
}


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        default=INPUT_H5AD,
        help="Path to the merged candidate milestone.",
    )
    parser.add_argument(
        "--output-h5ad",
        type=Path,
        default=OUTPUT_H5AD,
        help="Path to the cleaned Phase 2 milestone.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=CHUNK_SIZE,
        help="Rows per sparse chunk during evaluation.",
    )
    parser.add_argument(
        "--min-gene-cells",
        type=int,
        default=MIN_GENE_CELLS,
        help="Retain genes expressed in at least this many kept cells.",
    )
    parser.add_argument(
        "--sample-every-n",
        type=int,
        default=SAMPLE_EVERY_N,
        help="Keep every Nth evaluated cell in a light QC sample table.",
    )
    return parser.parse_args()


def setup_logging() -> None:
    """Configure concise logging."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )


def ensure_output_dirs() -> None:
    """Create required output directories."""
    for path in (TABLE_DIR, FIGURE_DIR, LOG_DIR):
        path.mkdir(parents=True, exist_ok=True)


def load_required_obs(input_h5ad: Path) -> pd.DataFrame:
    """Load only the metadata required for Phase 2 decisions."""
    logging.info("Loading required Phase 1 metadata from %s", input_h5ad)
    adata = ad.read_h5ad(input_h5ad, backed="r")
    missing = [column for column in REQUIRED_OBS_COLUMNS if column not in adata.obs.columns]
    if missing:
        if getattr(adata, "file", None) is not None:
            adata.file.close()
        raise ValueError(f"Missing required Phase 2 obs columns: {missing}")

    obs = adata.obs[REQUIRED_OBS_COLUMNS].copy()
    obs.insert(0, "obs_name", obs.index.astype(str))
    if getattr(adata, "file", None) is not None:
        adata.file.close()

    numeric_columns = [
        "phase1_t_score",
        "phase1_gd_score",
        "phase1_nk_score",
        "phase1_contam_score",
        "phase1_t_hits",
        "phase1_gd_hits",
        "phase1_nk_hits",
        "phase1_nk_strong_hits",
    ]
    for column in numeric_columns:
        obs[column] = pd.to_numeric(obs[column], errors="coerce").fillna(0)

    obs["source_gse_id"] = obs["source_gse_id"].fillna("unknown_gse").astype(str)
    obs["phase1_selection_reason"] = (
        obs["phase1_selection_reason"].fillna("").astype(str).str.strip()
    )
    obs["phase2_tnk_support"] = obs[
        ["phase1_t_score", "phase1_gd_score", "phase1_nk_score"]
    ].max(axis=1)
    obs["phase2_total_hits"] = (
        obs["phase1_t_hits"]
        + obs["phase1_gd_hits"]
        + obs["phase1_nk_hits"]
        + obs["phase1_nk_strong_hits"]
    )
    return obs


def build_marker_index(var_names: pd.Index) -> dict[str, np.ndarray]:
    """Resolve marker genes to integer column indices once."""
    marker_index: dict[str, np.ndarray] = {}
    mt_idx = [i for i, gene in enumerate(var_names) if gene.startswith("MT-")]

    for label, genes in MARKER_GROUPS.items():
        if label == "mitochondrial":
            marker_index[label] = np.asarray(mt_idx, dtype=np.int64)
            continue
        indices = [var_names.get_loc(gene) for gene in genes if gene in var_names]
        marker_index[label] = np.asarray(indices, dtype=np.int64)
    return marker_index


def safe_row_sums(matrix: sp.spmatrix | np.ndarray) -> np.ndarray:
    """Return row sums as a flat float array."""
    return np.asarray(matrix.sum(axis=1)).ravel().astype(np.float64, copy=False)


def compute_marker_score(
    matrix_chunk: sp.spmatrix | np.ndarray,
    gene_idx: np.ndarray,
    total_counts: np.ndarray,
) -> np.ndarray:
    """Compute a compact per-cell log-normalized mean marker score."""
    if gene_idx.size == 0:
        return np.zeros(matrix_chunk.shape[0], dtype=np.float32)

    marker_counts = safe_row_sums(matrix_chunk[:, gene_idx])
    mean_marker_counts = marker_counts / float(gene_idx.size)
    normalized = 1e4 * mean_marker_counts / np.clip(total_counts, 1.0, None)
    return np.log1p(normalized).astype(np.float32, copy=False)


def choose_primary_reason(reason_matrix: dict[str, np.ndarray]) -> np.ndarray:
    """Choose one primary removal reason per cell from the boolean reason matrix."""
    n_cells = len(next(iter(reason_matrix.values())))
    reason = np.full(n_cells, "kept", dtype=object)
    priority = [
        "low_quality",
        "extreme_phase1_contam",
        "myeloid_conflict",
        "bcell_conflict",
        "epithelial_conflict",
        "erythroid_conflict",
        "platelet_conflict",
    ]
    for label in priority:
        mask = reason_matrix[label] & (reason == "kept")
        reason[mask] = label
    return reason


def evaluate_cells_and_genes(
    input_h5ad: Path,
    obs: pd.DataFrame,
    chunk_size: int,
    sample_every_n: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame, pd.DataFrame, np.ndarray]:
    """Evaluate all cells, decide removals, and count kept-cell gene detection.

    Returns:
    - keep mask
    - gene detection counts across kept cells
    - gene names
    - removed cell QC table
    - sampled cell QC table
    - per-cell removal reason array
    """
    adata = ad.read_h5ad(input_h5ad, backed="r")
    marker_index = build_marker_index(adata.var_names)
    gene_detect_counts = np.zeros(adata.n_vars, dtype=np.int64)
    gene_names = adata.var_names.to_numpy(copy=True)
    keep_mask = np.ones(adata.n_obs, dtype=bool)
    reason_array = np.full(adata.n_obs, "kept", dtype=object)
    removed_parts: list[pd.DataFrame] = []
    sampled_parts: list[pd.DataFrame] = []
    n_chunks = (adata.n_obs + chunk_size - 1) // chunk_size

    for chunk_id, start in enumerate(range(0, adata.n_obs, chunk_size), start=1):
        end = min(start + chunk_size, adata.n_obs)
        matrix_chunk = adata.X[start:end]
        obs_chunk = obs.iloc[start:end].copy()

        total_counts = safe_row_sums(matrix_chunk)
        if sp.issparse(matrix_chunk):
            n_genes = np.asarray(matrix_chunk.getnnz(axis=1)).ravel().astype(np.int32, copy=False)
        else:
            n_genes = np.count_nonzero(matrix_chunk, axis=1).astype(np.int32, copy=False)

        pct_mt = np.zeros(matrix_chunk.shape[0], dtype=np.float32)
        mt_idx = marker_index["mitochondrial"]
        if mt_idx.size > 0:
            mt_counts = safe_row_sums(matrix_chunk[:, mt_idx])
            pct_mt = (mt_counts / np.clip(total_counts, 1.0, None)).astype(np.float32, copy=False)

        score_data: dict[str, np.ndarray] = {}
        for label in ["myeloid", "bcell", "epithelial", "erythroid", "platelet"]:
            score_data[label] = compute_marker_score(
                matrix_chunk=matrix_chunk,
                gene_idx=marker_index[label],
                total_counts=total_counts,
            )

        off_target_scores = np.vstack(
            [
                score_data["myeloid"],
                score_data["bcell"],
                score_data["epithelial"],
                score_data["erythroid"],
                score_data["platelet"],
            ]
        ).T
        off_target_labels = np.array(
            ["myeloid", "bcell", "epithelial", "erythroid", "platelet"],
            dtype=object,
        )
        dominant_off_target_idx = off_target_scores.argmax(axis=1)
        dominant_off_target_score = off_target_scores[
            np.arange(off_target_scores.shape[0]),
            dominant_off_target_idx,
        ]
        dominant_off_target_label = off_target_labels[dominant_off_target_idx]

        tnk_support = obs_chunk["phase2_tnk_support"].to_numpy(dtype=np.float32, copy=False)
        total_hits = obs_chunk["phase2_total_hits"].to_numpy(dtype=np.float32, copy=False)
        phase1_contam = obs_chunk["phase1_contam_score"].to_numpy(dtype=np.float32, copy=False)

        low_quality = (
            (n_genes < 100)
            | (total_counts < 150)
            | ((pct_mt >= 0.30) & (n_genes < 200))
        )
        extreme_phase1_contam = (
            (phase1_contam >= 1.75)
            & (tnk_support <= 0.50)
            & (total_hits == 0)
        )

        shared_conflict_gate = (
            (dominant_off_target_score >= 1.20)
            & (tnk_support <= 0.75)
            & (total_hits <= 1)
            & (dominant_off_target_score >= (tnk_support + 0.35))
        )
        myeloid_conflict = shared_conflict_gate & (dominant_off_target_label == "myeloid")
        bcell_conflict = shared_conflict_gate & (dominant_off_target_label == "bcell")
        epithelial_conflict = (
            (score_data["epithelial"] >= 1.00)
            & (tnk_support <= 0.75)
            & (total_hits <= 1)
            & (score_data["epithelial"] >= (tnk_support + 0.35))
        )
        erythroid_conflict = (
            (score_data["erythroid"] >= 1.00)
            & (tnk_support <= 0.75)
            & (total_hits <= 1)
            & (score_data["erythroid"] >= (tnk_support + 0.35))
        )
        platelet_conflict = (
            (score_data["platelet"] >= 1.00)
            & (tnk_support <= 0.75)
            & (total_hits <= 1)
            & (score_data["platelet"] >= (tnk_support + 0.35))
        )

        reason_matrix = {
            "low_quality": low_quality,
            "extreme_phase1_contam": extreme_phase1_contam,
            "myeloid_conflict": myeloid_conflict,
            "bcell_conflict": bcell_conflict,
            "epithelial_conflict": epithelial_conflict,
            "erythroid_conflict": erythroid_conflict,
            "platelet_conflict": platelet_conflict,
        }
        chunk_reason = choose_primary_reason(reason_matrix)
        chunk_keep = chunk_reason == "kept"

        keep_mask[start:end] = chunk_keep
        reason_array[start:end] = chunk_reason

        if np.any(~chunk_keep):
            removed = obs_chunk.loc[~chunk_keep, ["obs_name", "source_gse_id"]].copy()
            removed["phase2_remove_reason"] = chunk_reason[~chunk_keep]
            removed["phase2_n_counts"] = total_counts[~chunk_keep]
            removed["phase2_n_genes"] = n_genes[~chunk_keep]
            removed["phase2_pct_mt"] = pct_mt[~chunk_keep]
            removed["phase2_tnk_support"] = tnk_support[~chunk_keep]
            removed["phase2_phase1_contam_score"] = phase1_contam[~chunk_keep]
            removed["phase2_myeloid_score"] = score_data["myeloid"][~chunk_keep]
            removed["phase2_bcell_score"] = score_data["bcell"][~chunk_keep]
            removed["phase2_epithelial_score"] = score_data["epithelial"][~chunk_keep]
            removed["phase2_erythroid_score"] = score_data["erythroid"][~chunk_keep]
            removed["phase2_platelet_score"] = score_data["platelet"][~chunk_keep]
            removed["phase2_dominant_offtarget_label"] = dominant_off_target_label[~chunk_keep]
            removed["phase2_dominant_offtarget_score"] = dominant_off_target_score[~chunk_keep]
            removed_parts.append(removed)

        sampled_index = np.arange(start, end, sample_every_n, dtype=np.int64) - start
        if sampled_index.size > 0:
            sampled = obs_chunk.iloc[sampled_index][["obs_name", "source_gse_id"]].copy()
            sampled["phase2_remove_reason"] = chunk_reason[sampled_index]
            sampled["phase2_keep"] = chunk_keep[sampled_index]
            sampled["phase2_n_counts"] = total_counts[sampled_index]
            sampled["phase2_n_genes"] = n_genes[sampled_index]
            sampled["phase2_pct_mt"] = pct_mt[sampled_index]
            sampled["phase2_tnk_support"] = tnk_support[sampled_index]
            sampled["phase2_phase1_contam_score"] = phase1_contam[sampled_index]
            sampled["phase2_myeloid_score"] = score_data["myeloid"][sampled_index]
            sampled["phase2_bcell_score"] = score_data["bcell"][sampled_index]
            sampled["phase2_epithelial_score"] = score_data["epithelial"][sampled_index]
            sampled["phase2_erythroid_score"] = score_data["erythroid"][sampled_index]
            sampled["phase2_platelet_score"] = score_data["platelet"][sampled_index]
            sampled["phase2_dominant_offtarget_score"] = dominant_off_target_score[sampled_index]
            sampled_parts.append(sampled)

        kept_chunk = matrix_chunk[chunk_keep]
        if kept_chunk.shape[0] > 0:
            if sp.issparse(kept_chunk):
                gene_detect_counts += np.asarray(kept_chunk.getnnz(axis=0)).ravel()
            else:
                gene_detect_counts += np.count_nonzero(kept_chunk, axis=0)

        if chunk_id == 1 or chunk_id % 20 == 0 or chunk_id == n_chunks:
            logging.info(
                "Phase 2 evaluation progress: chunk %s / %s | removed so far=%s",
                chunk_id,
                n_chunks,
                int((~keep_mask[:end]).sum()),
            )

    if getattr(adata, "file", None) is not None:
        adata.file.close()

    removed_cells = (
        pd.concat(removed_parts, ignore_index=True)
        if removed_parts
        else pd.DataFrame(
            columns=[
                "obs_name",
                "source_gse_id",
                "phase2_remove_reason",
                "phase2_n_counts",
                "phase2_n_genes",
                "phase2_pct_mt",
                "phase2_tnk_support",
                "phase2_phase1_contam_score",
                "phase2_myeloid_score",
                "phase2_bcell_score",
                "phase2_epithelial_score",
                "phase2_erythroid_score",
                "phase2_platelet_score",
                "phase2_dominant_offtarget_label",
                "phase2_dominant_offtarget_score",
            ]
        )
    )
    sampled_cells = pd.concat(sampled_parts, ignore_index=True)
    return keep_mask, gene_detect_counts, gene_names, removed_cells, sampled_cells, reason_array


def write_subset_h5ad(
    input_h5ad: Path,
    output_h5ad: Path,
    keep_mask: np.ndarray,
    gene_keep_mask: np.ndarray,
) -> None:
    """Write the filtered milestone to a temporary file."""
    if output_h5ad.exists():
        output_h5ad.unlink()

    logging.info(
        "Writing Phase 2 output to %s with shape (%s, %s)",
        output_h5ad,
        int(keep_mask.sum()),
        int(gene_keep_mask.sum()),
    )
    adata = ad.read_h5ad(input_h5ad, backed="r")
    subset = adata[keep_mask, gene_keep_mask]
    subset.write_h5ad(output_h5ad, convert_strings_to_categoricals=False)
    if getattr(adata, "file", None) is not None:
        adata.file.close()


def validate_filtered_h5ad(output_h5ad: Path, expected_shape: tuple[int, int]) -> None:
    """Validate the temporary output before replacing the canonical milestone."""
    adata = ad.read_h5ad(output_h5ad, backed="r")
    actual_shape = adata.shape
    if getattr(adata, "file", None) is not None:
        adata.file.close()

    if actual_shape != expected_shape:
        raise ValueError(f"Phase 2 shape mismatch: expected {expected_shape}, got {actual_shape}")


def build_gse_summary(obs: pd.DataFrame, keep_mask: np.ndarray) -> pd.DataFrame:
    """Summarize before/after cell counts by GSE."""
    before = obs.groupby("source_gse_id").size().rename("cells_before")
    after = obs.loc[keep_mask].groupby("source_gse_id").size().rename("cells_after")
    summary = pd.concat([before, after], axis=1).fillna(0)
    summary["cells_before"] = summary["cells_before"].astype(int)
    summary["cells_after"] = summary["cells_after"].astype(int)
    summary["cells_removed"] = summary["cells_before"] - summary["cells_after"]
    summary["retained_fraction"] = np.where(
        summary["cells_before"] > 0,
        summary["cells_after"] / summary["cells_before"],
        0.0,
    )
    return summary.reset_index().sort_values(
        ["cells_before", "source_gse_id"],
        ascending=[False, True],
    ).reset_index(drop=True)


def build_gene_summary(
    gene_names: np.ndarray,
    gene_detect_counts: np.ndarray,
    gene_keep_mask: np.ndarray,
) -> pd.DataFrame:
    """Build the per-gene detection summary."""
    gene_summary = pd.DataFrame(
        {
            "gene": gene_names,
            "cells_detected": gene_detect_counts,
            "keep_gene": gene_keep_mask,
        }
    )
    gene_summary["key_tnk_gene"] = gene_summary["gene"].isin(KEY_TNK_GENES)
    return gene_summary.sort_values(
        ["keep_gene", "cells_detected", "gene"],
        ascending=[False, False, True],
    ).reset_index(drop=True)


def write_tables(
    obs: pd.DataFrame,
    keep_mask: np.ndarray,
    reason_array: np.ndarray,
    gse_summary: pd.DataFrame,
    gene_summary: pd.DataFrame,
    removed_cells: pd.DataFrame,
    sampled_cells: pd.DataFrame,
    min_gene_cells: int,
) -> pd.DataFrame:
    """Write Phase 2 QC tables and return the top-level summary."""
    reason_counts = (
        pd.Series(reason_array, name="phase2_remove_reason")
        .value_counts()
        .rename_axis("phase2_remove_reason")
        .reset_index(name="cells")
    )
    reason_lookup = dict(zip(reason_counts["phase2_remove_reason"], reason_counts["cells"]))

    key_gene_status = gene_summary.loc[gene_summary["key_tnk_gene"], ["gene", "keep_gene"]]
    cleanup_summary = pd.DataFrame(
        [
            {
                "cells_before": int(len(obs)),
                "cells_after": int(keep_mask.sum()),
                "cells_removed_total": int((~keep_mask).sum()),
                "removed_low_quality": int(reason_lookup.get("low_quality", 0)),
                "removed_extreme_phase1_contam": int(reason_lookup.get("extreme_phase1_contam", 0)),
                "removed_myeloid_conflict": int(reason_lookup.get("myeloid_conflict", 0)),
                "removed_bcell_conflict": int(reason_lookup.get("bcell_conflict", 0)),
                "removed_epithelial_conflict": int(reason_lookup.get("epithelial_conflict", 0)),
                "removed_erythroid_conflict": int(reason_lookup.get("erythroid_conflict", 0)),
                "removed_platelet_conflict": int(reason_lookup.get("platelet_conflict", 0)),
                "genes_before": int(len(gene_summary)),
                "genes_after": int(gene_summary["keep_gene"].sum()),
                "genes_removed_lt_min_cells": int((~gene_summary["keep_gene"]).sum()),
                "min_gene_cells": int(min_gene_cells),
                "key_tnk_genes_kept_n": int(key_gene_status["keep_gene"].sum()),
                "key_tnk_genes_dropped_n": int((~key_gene_status["keep_gene"]).sum()),
            }
        ]
    )

    cleanup_summary.to_csv(PHASE2_SUMMARY_CSV, index=False)
    gse_summary.to_csv(PHASE2_GSE_CSV, index=False)
    removed_cells.to_csv(PHASE2_REMOVED_CELLS_CSV, index=False)
    gene_summary.to_csv(PHASE2_GENE_CSV, index=False)
    sampled_cells.to_csv(PHASE2_SAMPLED_CELLS_CSV, index=False)
    return cleanup_summary


def plot_gse_retention(gse_summary: pd.DataFrame) -> None:
    """Plot before/after cell counts by GSE."""
    ordered = gse_summary.sort_values("cells_before", ascending=True)
    height = max(8.0, 0.32 * len(ordered) + 1.5)
    fig, ax = plt.subplots(figsize=(11, height))
    y = np.arange(len(ordered))
    ax.barh(y, ordered["cells_before"], color="#d9e0e8", label="Before Phase 2")
    ax.barh(y, ordered["cells_after"], color="#255f85", label="After Phase 2")
    ax.set_yticks(y)
    ax.set_yticklabels(ordered["source_gse_id"], fontsize=9)
    ax.set_xlabel("Cells")
    ax.set_title("Phase 2 Cell Retention By GSE")
    ax.set_xscale("log")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(PHASE2_GSE_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_removal_reasons(removed_cells: pd.DataFrame) -> None:
    """Plot counts for each Phase 2 removal reason."""
    counts = (
        removed_cells["phase2_remove_reason"].value_counts()
        .rename_axis("reason")
        .reset_index(name="cells")
    )
    if counts.empty:
        counts = pd.DataFrame({"reason": ["no_cells_removed"], "cells": [0]})

    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    ax.bar(counts["reason"], counts["cells"], color="#9f3d2f")
    ax.set_ylabel("Cells removed")
    ax.set_title("Phase 2 Cell Removal Reasons")
    ax.tick_params(axis="x", labelrotation=20)
    fig.tight_layout()
    fig.savefig(PHASE2_REMOVE_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_score_scatter(sampled_cells: pd.DataFrame) -> None:
    """Plot T/NK support against dominant off-target score on a QC sample."""
    fig, ax = plt.subplots(figsize=(7.5, 6.2))
    kept = sampled_cells["phase2_keep"].to_numpy(dtype=bool, copy=False)
    ax.scatter(
        sampled_cells.loc[kept, "phase2_tnk_support"],
        sampled_cells.loc[kept, "phase2_dominant_offtarget_score"],
        s=4,
        alpha=0.14,
        color="#4a8c63",
        label="Kept",
        rasterized=True,
    )
    ax.scatter(
        sampled_cells.loc[~kept, "phase2_tnk_support"],
        sampled_cells.loc[~kept, "phase2_dominant_offtarget_score"],
        s=6,
        alpha=0.35,
        color="#b24a3a",
        label="Removed",
        rasterized=True,
    )
    ax.axvline(0.75, color="#555555", linestyle="--", linewidth=1.2)
    ax.axhline(1.20, color="#555555", linestyle="--", linewidth=1.2)
    ax.set_xlabel("Phase 2 T/NK support")
    ax.set_ylabel("Dominant off-target score")
    ax.set_title("Phase 2 Off-Target Conflict Review")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(PHASE2_SCORE_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_gene_detection_distribution(gene_summary: pd.DataFrame, min_gene_cells: int) -> None:
    """Plot the gene detection distribution and the retention threshold."""
    detected = np.clip(gene_summary["cells_detected"].to_numpy(), 1, None)
    bins = np.logspace(0, np.log10(max(detected.max(), min_gene_cells)), 60)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.hist(detected, bins=bins, color="#4f7d8f", alpha=0.85)
    ax.axvline(min_gene_cells, color="#9f3d2f", linestyle="--", linewidth=2)
    ax.set_xscale("log")
    ax.set_xlabel("Cells expressing gene")
    ax.set_ylabel("Genes")
    ax.set_title("Phase 2 Gene Detection Distribution")
    ax.text(
        min_gene_cells * 1.05,
        ax.get_ylim()[1] * 0.9,
        f"Threshold = {min_gene_cells} cells",
        color="#9f3d2f",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(PHASE2_GENE_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_qc_summary(
    input_h5ad: Path,
    output_h5ad: Path,
    cleanup_summary: pd.DataFrame,
    gse_summary: pd.DataFrame,
    gene_summary: pd.DataFrame,
    removed_cells: pd.DataFrame,
    min_gene_cells: int,
) -> None:
    """Write the markdown QC summary for Phase 2."""
    summary_row = cleanup_summary.iloc[0]
    key_gene_status = gene_summary.loc[gene_summary["key_tnk_gene"], ["gene", "keep_gene"]]
    kept_key_genes = key_gene_status.loc[key_gene_status["keep_gene"], "gene"].tolist()
    dropped_key_genes = key_gene_status.loc[~key_gene_status["keep_gene"], "gene"].tolist()
    touched_gses = gse_summary.loc[gse_summary["cells_removed"] > 0, "source_gse_id"].tolist()
    removal_counts = (
        removed_cells["phase2_remove_reason"].value_counts().to_dict()
        if not removed_cells.empty
        else {}
    )
    removal_fraction = 0.0
    if summary_row["cells_before"] > 0:
        removal_fraction = summary_row["cells_removed_total"] / summary_row["cells_before"]

    lines = [
        "# Phase 2 QC Summary",
        "",
        "## Scope",
        f"- Input milestone: `{input_h5ad}`",
        f"- Output milestone: `{output_h5ad}`",
        "- Phase 2 goal: merged-context second-pass cleanup without aggressive loss of plausible T/NK states",
        f"- Gene filter: keep genes expressed in at least {min_gene_cells} cleaned cells",
        "",
        "## Cell cleanup",
        f"- Cells before: {int(summary_row['cells_before']):,}",
        f"- Cells after: {int(summary_row['cells_after']):,}",
        f"- Total cells removed: {int(summary_row['cells_removed_total']):,}",
        f"- Removal fraction: {removal_fraction:.4%}",
        f"- Removed for low quality: {int(removal_counts.get('low_quality', 0)):,}",
        f"- Removed for extreme carried Phase 1 contamination: {int(removal_counts.get('extreme_phase1_contam', 0)):,}",
        f"- Removed for myeloid conflict: {int(removal_counts.get('myeloid_conflict', 0)):,}",
        f"- Removed for B-cell conflict: {int(removal_counts.get('bcell_conflict', 0)):,}",
        f"- Removed for epithelial conflict: {int(removal_counts.get('epithelial_conflict', 0)):,}",
        f"- Removed for erythroid conflict: {int(removal_counts.get('erythroid_conflict', 0)):,}",
        f"- Removed for platelet conflict: {int(removal_counts.get('platelet_conflict', 0)):,}",
        f"- GSEs with any removed cells: {', '.join(touched_gses) if touched_gses else 'none'}",
        "",
        "## Gene cleanup",
        f"- Genes before: {int(summary_row['genes_before']):,}",
        f"- Genes after: {int(summary_row['genes_after']):,}",
        f"- Genes removed for detection < {min_gene_cells}: {int(summary_row['genes_removed_lt_min_cells']):,}",
        f"- Key T/NK genes kept: {', '.join(kept_key_genes) if kept_key_genes else 'none'}",
        f"- Key T/NK genes dropped: {', '.join(dropped_key_genes) if dropped_key_genes else 'none'}",
        "",
        "## Outputs",
        f"- Summary table: `{PHASE2_SUMMARY_CSV}`",
        f"- GSE before/after table: `{PHASE2_GSE_CSV}`",
        f"- Removed cells table: `{PHASE2_REMOVED_CELLS_CSV}`",
        f"- Gene detection table: `{PHASE2_GENE_CSV}`",
        f"- QC sample table: `{PHASE2_SAMPLED_CELLS_CSV}`",
        f"- Figures: `{PHASE2_GSE_PNG}`, `{PHASE2_REMOVE_PNG}`, `{PHASE2_SCORE_PNG}`, `{PHASE2_GENE_PNG}`",
        "",
        "## QC conclusion",
        "- This Phase 2 run remains marker-conflict driven and avoids broad cluster-level pruning.",
        "- If off-target removals stay small and key T/NK genes are preserved, the cleaned milestone is acceptable for Phase 3 integration.",
        "",
    ]
    PHASE2_QC_MD.write_text("\n".join(lines), encoding="utf-8")


def replace_canonical_output(temp_h5ad: Path, output_h5ad: Path) -> None:
    """Atomically replace the canonical Phase 2 output after validation."""
    logging.info("Replacing %s with validated Phase 2 output", output_h5ad)
    os.replace(temp_h5ad, output_h5ad)


def main() -> None:
    """Run the Phase 2 merged cleanup workflow."""
    args = parse_args()
    setup_logging()
    ensure_output_dirs()

    obs = load_required_obs(args.input_h5ad)
    (
        keep_mask,
        gene_detect_counts,
        gene_names,
        removed_cells,
        sampled_cells,
        reason_array,
    ) = evaluate_cells_and_genes(
        input_h5ad=args.input_h5ad,
        obs=obs,
        chunk_size=args.chunk_size,
        sample_every_n=args.sample_every_n,
    )

    gene_keep_mask = gene_detect_counts >= args.min_gene_cells
    if not np.any(gene_keep_mask):
        raise ValueError("The Phase 2 gene filter removed every gene.")

    gse_summary = build_gse_summary(obs=obs, keep_mask=keep_mask)
    gene_summary = build_gene_summary(
        gene_names=gene_names,
        gene_detect_counts=gene_detect_counts,
        gene_keep_mask=gene_keep_mask,
    )
    cleanup_summary = write_tables(
        obs=obs,
        keep_mask=keep_mask,
        reason_array=reason_array,
        gse_summary=gse_summary,
        gene_summary=gene_summary,
        removed_cells=removed_cells,
        sampled_cells=sampled_cells,
        min_gene_cells=args.min_gene_cells,
    )

    plot_gse_retention(gse_summary)
    plot_removal_reasons(removed_cells)
    plot_score_scatter(sampled_cells)
    plot_gene_detection_distribution(gene_summary, args.min_gene_cells)

    write_subset_h5ad(
        input_h5ad=args.input_h5ad,
        output_h5ad=TEMP_OUTPUT_H5AD,
        keep_mask=keep_mask,
        gene_keep_mask=gene_keep_mask,
    )
    validate_filtered_h5ad(
        output_h5ad=TEMP_OUTPUT_H5AD,
        expected_shape=(int(keep_mask.sum()), int(gene_keep_mask.sum())),
    )
    replace_canonical_output(TEMP_OUTPUT_H5AD, args.output_h5ad)

    write_qc_summary(
        input_h5ad=args.input_h5ad,
        output_h5ad=args.output_h5ad,
        cleanup_summary=cleanup_summary,
        gse_summary=gse_summary,
        gene_summary=gene_summary,
        removed_cells=removed_cells,
        min_gene_cells=args.min_gene_cells,
    )

    logging.info(
        "Phase 2 complete: cells_before=%s cells_after=%s genes_after=%s",
        int(cleanup_summary.iloc[0]["cells_before"]),
        int(cleanup_summary.iloc[0]["cells_after"]),
        int(cleanup_summary.iloc[0]["genes_after"]),
    )


if __name__ == "__main__":
    main()
