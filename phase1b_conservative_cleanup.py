#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Phase 1b conservative cleanup for the merged T/NK candidate milestone.

This step is intentionally conservative. It removes only:
1. obvious non-T/NK cells that slipped through annotation-only retention
2. an optional, very strict high-confidence T/NK doublet pattern
3. genes expressed in fewer than the configured number of cells

The milestone file remains `Integrated_dataset/TNK_candidates.h5ad`.
The script writes a temporary H5AD, validates it, and then replaces the
canonical milestone in place.
"""

from __future__ import annotations

import argparse
import gc
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
OUTPUT_ROOT = Path("Integrated_dataset")
TABLE_DIR = OUTPUT_ROOT / "tables"
FIGURE_DIR = OUTPUT_ROOT / "figures"
LOG_DIR = OUTPUT_ROOT / "logs"
TEMP_OUTPUT_H5AD = OUTPUT_ROOT / "TNK_candidates.phase1b_tmp.h5ad"

CONTAM_THRESHOLD = 1.5
DOUBLET_CONTAM_THRESHOLD = 2.0
DOUBLET_T_SCORE_THRESHOLD = 1.5
DOUBLET_NK_SCORE_THRESHOLD = 1.5
MIN_GENE_CELLS = 500
CHUNK_SIZE = 25_000

PHASE1B_SUMMARY_CSV = TABLE_DIR / "phase1b_cleanup_summary.csv"
PHASE1B_GSE_CSV = TABLE_DIR / "phase1b_gse_before_after.csv"
PHASE1B_REMOVED_CELLS_CSV = TABLE_DIR / "phase1b_removed_cells.csv"
PHASE1B_GENE_CSV = TABLE_DIR / "phase1b_gene_detection_summary.csv"
PHASE1B_QC_MD = LOG_DIR / "phase1b_qc_summary.md"

PHASE1B_GSE_PNG = FIGURE_DIR / "phase1b_gse_cell_retention.png"
PHASE1B_REMOVE_PNG = FIGURE_DIR / "phase1b_cell_removal_reasons.png"
PHASE1B_GENE_PNG = FIGURE_DIR / "phase1b_gene_detection_distribution.png"

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
    "phase1_annotation_label",
    "phase1_contam_score",
    "phase1_t_score",
    "phase1_gd_score",
    "phase1_nk_score",
    "phase1_t_hits",
    "phase1_gd_hits",
    "phase1_nk_hits",
    "phase1_nk_strong_hits",
]


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        default=INPUT_H5AD,
        help="Path to the merged Phase 1 candidate milestone.",
    )
    parser.add_argument(
        "--min-gene-cells",
        type=int,
        default=MIN_GENE_CELLS,
        help="Drop genes expressed in fewer than this many kept cells.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=CHUNK_SIZE,
        help="Row chunk size for backed gene-detection counting.",
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


def load_phase1_obs(input_h5ad: Path) -> pd.DataFrame:
    """Load only the obs columns needed for Phase 1b decisions."""
    logging.info("Loading Phase 1 metadata from %s", input_h5ad)
    adata = ad.read_h5ad(input_h5ad, backed="r")
    missing = [col for col in REQUIRED_OBS_COLUMNS if col not in adata.obs.columns]
    if missing:
        if getattr(adata, "file", None) is not None:
            adata.file.close()
        raise ValueError(f"Missing required Phase 1 obs columns: {missing}")

    obs = adata.obs[REQUIRED_OBS_COLUMNS].copy()
    if getattr(adata, "file", None) is not None:
        adata.file.close()

    for column in [
        "phase1_contam_score",
        "phase1_t_score",
        "phase1_gd_score",
        "phase1_nk_score",
        "phase1_t_hits",
        "phase1_gd_hits",
        "phase1_nk_hits",
        "phase1_nk_strong_hits",
    ]:
        obs[column] = pd.to_numeric(obs[column], errors="coerce").fillna(0)

    obs["phase1_selection_reason"] = (
        obs["phase1_selection_reason"].fillna("").astype(str).str.strip()
    )
    obs["phase1_annotation_label"] = (
        obs["phase1_annotation_label"].fillna("").astype(str).str.strip()
    )
    obs["source_gse_id"] = obs["source_gse_id"].fillna("unknown_gse").astype(str)
    obs["phase1_total_hits"] = (
        obs["phase1_t_hits"]
        + obs["phase1_gd_hits"]
        + obs["phase1_nk_hits"]
        + obs["phase1_nk_strong_hits"]
    )
    return obs


def build_phase1b_decisions(obs: pd.DataFrame) -> pd.DataFrame:
    """Create conservative removal decisions from Phase 1 metadata.

    The main removal rule only targets cells that were kept solely by
    annotation, have no T/NK marker support, and carry a strong contaminant
    score. A much stricter T/NK doublet pattern is tracked separately.
    """
    obvious_non_tnk = (
        (obs["phase1_contam_score"] >= CONTAM_THRESHOLD)
        & obs["phase1_selection_reason"].eq("annotation_only")
        & obs["phase1_total_hits"].eq(0)
        & obs["phase1_t_score"].eq(0)
        & obs["phase1_gd_score"].eq(0)
        & obs["phase1_nk_score"].eq(0)
    )

    high_conf_doublet = (
        (obs["phase1_contam_score"] >= DOUBLET_CONTAM_THRESHOLD)
        & (obs["phase1_t_score"] >= DOUBLET_T_SCORE_THRESHOLD)
        & (obs["phase1_nk_score"] >= DOUBLET_NK_SCORE_THRESHOLD)
    )

    decisions = obs.copy()
    decisions["phase1b_remove_obvious_non_tnk"] = obvious_non_tnk
    decisions["phase1b_remove_high_conf_doublet"] = high_conf_doublet
    decisions["phase1b_remove"] = obvious_non_tnk | high_conf_doublet
    decisions["phase1b_remove_reason"] = np.where(
        obvious_non_tnk,
        "obvious_non_tnk",
        np.where(high_conf_doublet, "high_conf_doublet", "kept"),
    )
    return decisions


def count_cells_per_gene(
    input_h5ad: Path,
    cell_keep_mask: np.ndarray,
    chunk_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Count how many kept cells express each gene without loading all data.

    This uses backed row slices so the large matrix stays on disk. Only a
    moderate sparse chunk is materialized at a time.
    """
    adata = ad.read_h5ad(input_h5ad, backed="r")
    gene_names = adata.var_names.to_numpy(copy=True)
    gene_detect_counts = np.zeros(adata.n_vars, dtype=np.int64)
    n_obs = adata.n_obs
    n_chunks = (n_obs + chunk_size - 1) // chunk_size

    for chunk_idx, start in enumerate(range(0, n_obs, chunk_size), start=1):
        end = min(start + chunk_size, n_obs)
        local_keep = cell_keep_mask[start:end]
        if not np.any(local_keep):
            continue

        matrix_chunk = adata.X[start:end]
        matrix_chunk = matrix_chunk[local_keep]

        if sp.issparse(matrix_chunk):
            gene_detect_counts += np.asarray(matrix_chunk.getnnz(axis=0)).ravel()
        else:
            gene_detect_counts += np.count_nonzero(matrix_chunk, axis=0)

        if chunk_idx == 1 or chunk_idx % 20 == 0 or chunk_idx == n_chunks:
            logging.info(
                "Gene detection count progress: chunk %s / %s",
                chunk_idx,
                n_chunks,
            )

    if getattr(adata, "file", None) is not None:
        adata.file.close()
    return gene_detect_counts, gene_names


def write_subset_h5ad(
    input_h5ad: Path,
    output_h5ad: Path,
    cell_keep_mask: np.ndarray,
    gene_keep_mask: np.ndarray,
) -> None:
    """Write the filtered candidate object to a temporary H5AD on disk."""
    if output_h5ad.exists():
        output_h5ad.unlink()

    logging.info(
        "Writing filtered milestone to %s with shape (%s, %s)",
        output_h5ad,
        int(cell_keep_mask.sum()),
        int(gene_keep_mask.sum()),
    )
    adata = ad.read_h5ad(input_h5ad, backed="r")
    subset = adata[cell_keep_mask, gene_keep_mask]
    subset.write_h5ad(output_h5ad, convert_strings_to_categoricals=False)
    if getattr(adata, "file", None) is not None:
        adata.file.close()


def validate_filtered_h5ad(
    output_h5ad: Path,
    expected_n_obs: int,
    expected_n_vars: int,
) -> None:
    """Validate that the temporary filtered milestone has the expected shape."""
    filtered = ad.read_h5ad(output_h5ad, backed="r")
    actual_shape = filtered.shape
    if getattr(filtered, "file", None) is not None:
        filtered.file.close()

    if actual_shape != (expected_n_obs, expected_n_vars):
        raise ValueError(
            "Filtered H5AD shape mismatch: "
            f"expected {(expected_n_obs, expected_n_vars)}, got {actual_shape}"
        )


def build_gse_summary(decisions: pd.DataFrame, kept_decisions: pd.DataFrame) -> pd.DataFrame:
    """Summarize cell counts before and after cleanup for each GSE."""
    before = decisions.groupby("source_gse_id").size().rename("cells_before")
    after = kept_decisions.groupby("source_gse_id").size().rename("cells_after")
    summary = pd.concat([before, after], axis=1).fillna(0)
    summary["cells_before"] = summary["cells_before"].astype(int)
    summary["cells_after"] = summary["cells_after"].astype(int)
    summary["cells_removed"] = summary["cells_before"] - summary["cells_after"]
    summary["retained_fraction"] = np.where(
        summary["cells_before"] > 0,
        summary["cells_after"] / summary["cells_before"],
        0.0,
    )
    summary = summary.reset_index()
    return summary.sort_values(
        ["cells_before", "source_gse_id"],
        ascending=[False, True],
    ).reset_index(drop=True)


def build_gene_summary(
    gene_names: np.ndarray,
    gene_detect_counts: np.ndarray,
    gene_keep_mask: np.ndarray,
) -> pd.DataFrame:
    """Build the per-gene detection summary table."""
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
    decisions: pd.DataFrame,
    kept_decisions: pd.DataFrame,
    gse_summary: pd.DataFrame,
    gene_summary: pd.DataFrame,
    min_gene_cells: int,
) -> pd.DataFrame:
    """Write Phase 1b QC tables and return the top-level cleanup summary."""
    removal_counts = (
        decisions["phase1b_remove_reason"]
        .value_counts()
        .rename_axis("phase1b_remove_reason")
        .reset_index(name="cells")
    )
    removal_lookup = dict(zip(removal_counts["phase1b_remove_reason"], removal_counts["cells"]))

    key_gene_status = gene_summary.loc[gene_summary["key_tnk_gene"], ["gene", "cells_detected", "keep_gene"]]
    key_gene_kept = key_gene_status.loc[key_gene_status["keep_gene"], "gene"].tolist()
    key_gene_dropped = key_gene_status.loc[~key_gene_status["keep_gene"], "gene"].tolist()

    cleanup_summary = pd.DataFrame(
        [
            {
                "cells_before": int(len(decisions)),
                "cells_after": int(len(kept_decisions)),
                "cells_removed_total": int(decisions["phase1b_remove"].sum()),
                "removed_obvious_non_tnk": int(removal_lookup.get("obvious_non_tnk", 0)),
                "removed_high_conf_doublet": int(removal_lookup.get("high_conf_doublet", 0)),
                "genes_before": int(len(gene_summary)),
                "genes_after": int(gene_summary["keep_gene"].sum()),
                "genes_removed_lt_min_cells": int((~gene_summary["keep_gene"]).sum()),
                "min_gene_cells": int(min_gene_cells),
                "key_tnk_genes_kept_n": int(len(key_gene_kept)),
                "key_tnk_genes_dropped_n": int(len(key_gene_dropped)),
            }
        ]
    )

    removed_cells = decisions.loc[decisions["phase1b_remove"]].copy()
    removed_cells.insert(0, "obs_name", removed_cells.index.astype(str))

    cleanup_summary.to_csv(PHASE1B_SUMMARY_CSV, index=False)
    gse_summary.to_csv(PHASE1B_GSE_CSV, index=False)
    removed_cells.to_csv(PHASE1B_REMOVED_CELLS_CSV, index=False)
    gene_summary.to_csv(PHASE1B_GENE_CSV, index=False)
    return cleanup_summary


def plot_gse_retention(gse_summary: pd.DataFrame) -> None:
    """Plot before/after cell counts by GSE."""
    ordered = gse_summary.sort_values("cells_before", ascending=True)
    height = max(8.0, 0.32 * len(ordered) + 1.5)
    fig, ax = plt.subplots(figsize=(11, height))
    y = np.arange(len(ordered))

    ax.barh(y, ordered["cells_before"], color="#c7d3dd", label="Before Phase 1b")
    ax.barh(y, ordered["cells_after"], color="#2f6b8a", label="After Phase 1b")
    ax.set_yticks(y)
    ax.set_yticklabels(ordered["source_gse_id"], fontsize=9)
    ax.set_xlabel("Cells")
    ax.set_title("Phase 1b Cell Retention By GSE")
    ax.set_xscale("log")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(PHASE1B_GSE_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_removal_reasons(decisions: pd.DataFrame) -> None:
    """Plot the small set of Phase 1b cell removal reasons."""
    counts = (
        decisions.loc[decisions["phase1b_remove"], "phase1b_remove_reason"]
        .value_counts()
        .rename_axis("reason")
        .reset_index(name="cells")
    )
    if counts.empty:
        counts = pd.DataFrame({"reason": ["no_cells_removed"], "cells": [0]})

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(counts["reason"], counts["cells"], color=["#a53f2b"] * len(counts))
    ax.set_ylabel("Cells removed")
    ax.set_title("Phase 1b Cell Removal Reasons")
    ax.tick_params(axis="x", labelrotation=15)
    fig.tight_layout()
    fig.savefig(PHASE1B_REMOVE_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_gene_detection_distribution(gene_summary: pd.DataFrame, min_gene_cells: int) -> None:
    """Plot the gene detection distribution and the retention threshold."""
    fig, ax = plt.subplots(figsize=(9, 5.5))
    clipped = np.clip(gene_summary["cells_detected"].to_numpy(), 1, None)
    bins = np.logspace(0, np.log10(max(clipped.max(), min_gene_cells)), 60)

    ax.hist(clipped, bins=bins, color="#4c8c5b", alpha=0.85)
    ax.axvline(min_gene_cells, color="#a53f2b", linestyle="--", linewidth=2)
    ax.set_xscale("log")
    ax.set_xlabel("Cells expressing gene")
    ax.set_ylabel("Genes")
    ax.set_title("Phase 1b Gene Detection Distribution")
    ax.text(
        min_gene_cells * 1.05,
        ax.get_ylim()[1] * 0.9,
        f"Threshold = {min_gene_cells} cells",
        color="#a53f2b",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(PHASE1B_GENE_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_qc_summary(
    input_h5ad: Path,
    cleanup_summary: pd.DataFrame,
    gse_summary: pd.DataFrame,
    gene_summary: pd.DataFrame,
    min_gene_cells: int,
) -> None:
    """Write the Phase 1b QC markdown summary."""
    summary_row = cleanup_summary.iloc[0]
    key_gene_status = gene_summary.loc[gene_summary["key_tnk_gene"], ["gene", "cells_detected", "keep_gene"]]
    kept_key_genes = key_gene_status.loc[key_gene_status["keep_gene"], "gene"].tolist()
    dropped_key_genes = key_gene_status.loc[~key_gene_status["keep_gene"], "gene"].tolist()
    touched_gses = gse_summary.loc[gse_summary["cells_removed"] > 0, "source_gse_id"].tolist()

    lines = [
        "# Phase 1b QC Summary",
        "",
        "## Scope",
        f"- Input milestone: `{input_h5ad}`",
        "- Phase 1b goal: conservative first-pass cleanup only",
        f"- Gene filter added per user instruction: remove genes expressed in fewer than {min_gene_cells} cells",
        "",
        "## Cell cleanup",
        f"- Cells before: {int(summary_row['cells_before']):,}",
        f"- Cells after: {int(summary_row['cells_after']):,}",
        f"- Total cells removed: {int(summary_row['cells_removed_total']):,}",
        f"- Removed as obvious non-T/NK: {int(summary_row['removed_obvious_non_tnk']):,}",
        f"- Removed as high-confidence doublets: {int(summary_row['removed_high_conf_doublet']):,}",
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
        f"- Summary table: `{PHASE1B_SUMMARY_CSV}`",
        f"- GSE before/after table: `{PHASE1B_GSE_CSV}`",
        f"- Removed cell table: `{PHASE1B_REMOVED_CELLS_CSV}`",
        f"- Gene detection table: `{PHASE1B_GENE_CSV}`",
        f"- Figures: `{PHASE1B_GSE_PNG}`, `{PHASE1B_REMOVE_PNG}`, `{PHASE1B_GENE_PNG}`",
        "",
        "## QC conclusion",
        "- Phase 1b remains conservative: only high-confidence obvious contaminants are removed at the cell level.",
        "- No broad T/NK-like population is removed by this rule set.",
        "- This milestone is ready for user QC review before any Phase 1c / Phase 2 work.",
        "",
    ]
    PHASE1B_QC_MD.write_text("\n".join(lines), encoding="utf-8")


def replace_canonical_milestone(temp_h5ad: Path, canonical_h5ad: Path) -> None:
    """Atomically replace the canonical milestone after validation."""
    logging.info("Replacing %s with validated Phase 1b output", canonical_h5ad)
    os.replace(temp_h5ad, canonical_h5ad)


def main() -> None:
    """Run Phase 1b cleanup and rewrite the milestone in place."""
    args = parse_args()
    setup_logging()
    ensure_output_dirs()

    decisions = build_phase1b_decisions(load_phase1_obs(args.input_h5ad))
    cell_keep_mask = (~decisions["phase1b_remove"]).to_numpy()
    kept_decisions = decisions.loc[~decisions["phase1b_remove"]].copy()

    logging.info(
        "Phase 1b cell removal counts: obvious_non_tnk=%s, high_conf_doublet=%s",
        int(decisions["phase1b_remove_obvious_non_tnk"].sum()),
        int(decisions["phase1b_remove_high_conf_doublet"].sum()),
    )

    gene_detect_counts, gene_names = count_cells_per_gene(
        input_h5ad=args.input_h5ad,
        cell_keep_mask=cell_keep_mask,
        chunk_size=args.chunk_size,
    )
    gene_keep_mask = gene_detect_counts >= args.min_gene_cells
    if not np.any(gene_keep_mask):
        raise ValueError("The configured gene filter removed every gene.")

    gse_summary = build_gse_summary(decisions, kept_decisions)
    gene_summary = build_gene_summary(gene_names, gene_detect_counts, gene_keep_mask)
    cleanup_summary = write_tables(
        decisions=decisions,
        kept_decisions=kept_decisions,
        gse_summary=gse_summary,
        gene_summary=gene_summary,
        min_gene_cells=args.min_gene_cells,
    )

    plot_gse_retention(gse_summary)
    plot_removal_reasons(decisions)
    plot_gene_detection_distribution(gene_summary, args.min_gene_cells)
    write_qc_summary(
        input_h5ad=args.input_h5ad,
        cleanup_summary=cleanup_summary,
        gse_summary=gse_summary,
        gene_summary=gene_summary,
        min_gene_cells=args.min_gene_cells,
    )

    write_subset_h5ad(
        input_h5ad=args.input_h5ad,
        output_h5ad=TEMP_OUTPUT_H5AD,
        cell_keep_mask=cell_keep_mask,
        gene_keep_mask=gene_keep_mask,
    )
    validate_filtered_h5ad(
        output_h5ad=TEMP_OUTPUT_H5AD,
        expected_n_obs=int(cell_keep_mask.sum()),
        expected_n_vars=int(gene_keep_mask.sum()),
    )
    replace_canonical_milestone(TEMP_OUTPUT_H5AD, args.input_h5ad)

    logging.info(
        "Phase 1b completed: cells %s -> %s, genes %s -> %s",
        int(cleanup_summary.iloc[0]["cells_before"]),
        int(cleanup_summary.iloc[0]["cells_after"]),
        int(cleanup_summary.iloc[0]["genes_before"]),
        int(cleanup_summary.iloc[0]["genes_after"]),
    )
    gc.collect()


if __name__ == "__main__":
    main()
