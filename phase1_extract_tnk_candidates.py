#!/usr/bin/env python3
"""Phase 1 coarse T/NK extraction for approved Category A datasets."""

from __future__ import annotations

import math
import re
import shutil
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns


# Config
PROJECT_ROOT = Path(__file__).resolve().parent
AUDIT_CSV = PROJECT_ROOT / "Integrated_dataset" / "tables" / "phase0_dataset_audit.csv"
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
FIGURE_DIR = OUTPUT_ROOT / "figures"
TABLE_DIR = OUTPUT_ROOT / "tables"
LOG_DIR = OUTPUT_ROOT / "logs"
TMP_DIR = OUTPUT_ROOT / "_tmp_phase1_candidates"
TNK_CANDIDATES_H5AD = OUTPUT_ROOT / "TNK_candidates.h5ad"

PHASE1_SUMMARY_CSV = TABLE_DIR / "phase1_categoryA_selection_summary.csv"
PHASE1_MARKER_CSV = TABLE_DIR / "phase1_categoryA_marker_availability.csv"
PHASE1_QC_MD = LOG_DIR / "phase1_qc_summary.md"

ANNOTATION_COLUMN_HINTS = (
    "celltype",
    "cell_type",
    "annotation",
    "annot",
    "label",
    "cluster",
    "lineage",
    "major",
    "sub_cell_type",
    "pred",
    "identity",
    "ident",
)
DONOR_HINTS = ("donor", "patient", "subject", "individual")
SAMPLE_HINTS = ("sample", "orig.ident", "specimen", "biosample")
LIBRARY_HINTS = ("library", "library_id", "batch", "lane", "channel")

ANNOTATION_POSITIVE_PATTERNS = [
    re.compile(pattern, flags=re.IGNORECASE)
    for pattern in [
        r"(?:^|[^A-Za-z])T(?:[^A-Za-z]|$)",
        r"T/?NK",
        r"CD3",
        r"CD4",
        r"CD8",
        r"NK",
        r"MAIT",
        r"gamma.?delta",
        r"gdt",
        r"TRDV",
        r"TRGC",
        r"TRDC",
        r"Treg",
    ]
]
ANNOTATION_NEGATIVE_PATTERNS = [
    re.compile(pattern, flags=re.IGNORECASE)
    for pattern in [
        r"(?:^|[^A-Za-z])B(?:[^A-Za-z]|$)",
        r"plasma",
        r"M/DC",
        r"MNP",
        r"myeloid",
        r"mono",
        r"macroph",
        r"DC",
        r"pDC",
        r"cDC",
        r"eryth",
        r"RBC",
        r"mast",
        r"fibro",
        r"strom",
        r"epithelial",
        r"doublet",
    ]
]

T_MARKERS = ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "IL7R", "LTB"]
GD_MARKERS = ["TRDC", "TRGC1", "TRGC2", "TRDV1", "TRDV2", "TRDV3"]
NK_MARKERS = ["NKG7", "KLRD1", "GNLY", "PRF1", "CTSW", "KLRC1", "XCL1", "XCL2"]
NK_STRONG_MARKERS = ["KLRD1", "GNLY", "PRF1", "CTSW", "XCL1", "XCL2"]
CONTAM_MARKERS = [
    "MS4A1",
    "CD79A",
    "CD79B",
    "MZB1",
    "LST1",
    "S100A8",
    "S100A9",
    "FCER1A",
    "C1QA",
    "C1QB",
    "EPCAM",
    "KRT8",
    "KRT18",
    "COL1A1",
    "DCN",
    "HBA1",
    "HBB",
    "PPBP",
]


def canonical_gene_symbol(name: object) -> str:
    """Normalize common prefixed feature names to a marker-matchable symbol.

    Public datasets sometimes store genes as `GRCh38_TRAC`, `hg38_TRAC`, or
    similar genome-prefixed names. Phase 1 marker matching should treat those as
    `TRAC` rather than missing markers.
    """
    text = str(name).strip()
    if not text:
        return ""

    upper = text.upper()
    prefix_patterns = (
        r"^GRCH\d+[_:-]+",
        r"^HG\d+[_:-]+",
        r"^MM\d+[_:-]+",
        r"^GRCM\d+[_:-]+",
    )
    for pattern in prefix_patterns:
        cleaned = re.sub(pattern, "", upper)
        if cleaned != upper:
            return cleaned
    return upper


def sanitize_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Convert extension dtypes to portable pandas types before writing H5AD."""
    clean = df.copy()
    for column in clean.columns:
        series = clean[column]
        if isinstance(series.dtype, pd.CategoricalDtype):
            clean[column] = series.astype(object)
        elif isinstance(series.dtype, pd.StringDtype):
            clean[column] = series.astype(object)
    return clean


def ensure_output_dirs() -> None:
    """Ensure the canonical output directories exist."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)


def configure_plotting() -> None:
    """Set publication-readable plotting defaults."""
    sns.set_theme(style="whitegrid")
    plt.rcParams.update(
        {
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 9,
        }
    )


def load_category_a_registry() -> pd.DataFrame:
    """Return the approved Category A datasets from the Phase 0 audit."""
    audit = pd.read_csv(AUDIT_CSV)
    category_a = audit.loc[audit["phase0_category"] == "A", ["gse_id", "h5ad_path"]].copy()
    if category_a.empty:
        raise ValueError("No Category A datasets found in phase0_dataset_audit.csv")
    return category_a.sort_values("gse_id").reset_index(drop=True)


def find_annotation_columns(obs_columns: list[str]) -> list[str]:
    """Pick likely annotation columns from obs metadata."""
    selected = []
    for column in obs_columns:
        lower = column.lower()
        if any(hint in lower for hint in ANNOTATION_COLUMN_HINTS):
            selected.append(column)
    return selected


def first_matching_column(obs_columns: list[str], hints: tuple[str, ...]) -> str:
    """Return the first obs column whose name matches one of the requested hints."""
    for hint in hints:
        for column in obs_columns:
            if hint in column.lower():
                return column
    return ""


def obs_label_series(obs: pd.DataFrame, column_name: str) -> pd.Series:
    """Return a portable string label series for a chosen metadata column."""
    if not column_name:
        return pd.Series([""] * obs.shape[0], index=obs.index, dtype=object)

    values = obs[column_name].astype(str).replace({"nan": "", "<NA>": "", "None": ""})
    return values.astype(object)


def marker_positions_from_var_names(var_names: pd.Index, markers: list[str]) -> np.ndarray:
    """Return first-occurrence positions for requested markers, case-insensitively."""
    first_pos: dict[str, int] = {}
    for idx, name in enumerate(var_names.astype(str)):
        canonical = canonical_gene_symbol(name)
        if canonical and canonical not in first_pos:
            first_pos[canonical] = idx

    positions = [first_pos[canonical_gene_symbol(gene)] for gene in markers if canonical_gene_symbol(gene) in first_pos]
    return np.asarray(positions, dtype=int)


def normalized_marker_scores(adata: ad.AnnData, markers: list[str]) -> tuple[np.ndarray, int]:
    """Compute a mean log1p-normalized marker score without densifying full X."""
    if adata.n_obs == 0:
        return np.zeros(0, dtype=np.float32), 0

    positions = marker_positions_from_var_names(pd.Index(adata.var_names), markers)
    if len(positions) == 0:
        return np.zeros(adata.n_obs, dtype=np.float32), 0

    matrix = adata.X[:, positions]
    total_counts = np.ravel(np.asarray(adata.X.sum(axis=1))).astype(np.float32)
    total_counts[total_counts == 0] = 1.0

    if sp.issparse(matrix):
        matrix = matrix.tocsr().astype(np.float32, copy=False)
        scale = sp.diags(1e4 / total_counts, format="csr")
        normalized = scale @ matrix
        normalized.data = np.log1p(normalized.data)
        scores = np.ravel(np.asarray(normalized.mean(axis=1))).astype(np.float32)
    else:
        dense = np.asarray(matrix, dtype=np.float32)
        dense = np.log1p((dense / total_counts[:, None]) * 1e4)
        scores = dense.mean(axis=1).astype(np.float32)
    return scores, len(positions)


def positive_marker_hits(adata: ad.AnnData, markers: list[str]) -> tuple[np.ndarray, int]:
    """Count how many markers are detected per cell."""
    if adata.n_obs == 0:
        return np.zeros(0, dtype=np.int16), 0

    positions = marker_positions_from_var_names(pd.Index(adata.var_names), markers)
    if len(positions) == 0:
        return np.zeros(adata.n_obs, dtype=np.int16), 0

    matrix = adata.X[:, positions]
    if sp.issparse(matrix):
        hits = np.ravel((matrix > 0).sum(axis=1)).astype(np.int16)
    else:
        hits = np.sum(np.asarray(matrix) > 0, axis=1).astype(np.int16)
    return hits, len(positions)


def annotation_masks(obs: pd.DataFrame, annotation_columns: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Return positive and explicitly negative annotation masks."""
    if not annotation_columns:
        empty = np.zeros(obs.shape[0], dtype=bool)
        return empty, empty

    keep = np.zeros(obs.shape[0], dtype=bool)
    negative = np.zeros(obs.shape[0], dtype=bool)
    for column in annotation_columns:
        values = obs[column].astype(str)
        col_keep = np.zeros(obs.shape[0], dtype=bool)
        col_negative = np.zeros(obs.shape[0], dtype=bool)
        for pattern in ANNOTATION_POSITIVE_PATTERNS:
            col_keep |= values.str.contains(pattern, na=False).to_numpy()
        for pattern in ANNOTATION_NEGATIVE_PATTERNS:
            col_negative |= values.str.contains(pattern, na=False).to_numpy()
        keep |= col_keep
        negative |= col_negative & ~col_keep
    return keep, negative


def selection_reason(annotation_keep: np.ndarray, marker_keep: np.ndarray) -> np.ndarray:
    """Encode why a cell was retained."""
    reasons = np.full(annotation_keep.shape[0], "marker_only", dtype=object)
    reasons[annotation_keep & ~marker_keep] = "annotation_only"
    reasons[annotation_keep & marker_keep] = "annotation_and_marker"
    return reasons


def extract_one_dataset(gse_id: str, h5ad_path: str) -> tuple[ad.AnnData, dict, list[dict]]:
    """Run coarse T/NK extraction for a single dataset."""
    adata = ad.read_h5ad(h5ad_path)
    adata.var_names_make_unique()

    annotation_columns = find_annotation_columns(list(adata.obs.columns))
    ann_keep, ann_negative = annotation_masks(adata.obs, annotation_columns)
    donor_column = first_matching_column(list(adata.obs.columns), DONOR_HINTS)
    sample_column = first_matching_column(list(adata.obs.columns), SAMPLE_HINTS)
    library_column = first_matching_column(list(adata.obs.columns), LIBRARY_HINTS)
    annotation_label_column = annotation_columns[0] if annotation_columns else ""

    t_hits, t_present = positive_marker_hits(adata, T_MARKERS)
    gd_hits, gd_present = positive_marker_hits(adata, GD_MARKERS)
    nk_hits, nk_present = positive_marker_hits(adata, NK_MARKERS)
    nk_strong_hits, nk_strong_present = positive_marker_hits(adata, NK_STRONG_MARKERS)

    t_score, t_score_present = normalized_marker_scores(adata, T_MARKERS)
    gd_score, gd_score_present = normalized_marker_scores(adata, GD_MARKERS)
    nk_score, nk_score_present = normalized_marker_scores(adata, NK_MARKERS)
    contam_score, contam_present = normalized_marker_scores(adata, CONTAM_MARKERS)

    t_rule = t_hits >= 1
    gd_rule = gd_hits >= 1
    nk_rule = (nk_hits >= 2) & (nk_strong_hits >= 1)
    score_rule = (t_score >= 0.60) | (gd_score >= 0.35) | (nk_score >= 0.80)
    marker_keep = t_rule | gd_rule | nk_rule | score_rule

    # A weak contamination override is only applied when there is no annotation support.
    strong_contam = contam_score >= 1.25
    weak_lymphoid = (t_hits + gd_hits + nk_hits) <= 1
    marker_keep &= ~(strong_contam & weak_lymphoid & ~ann_keep)

    keep_mask = ann_keep | (marker_keep & ~ann_negative)
    reasons = selection_reason(ann_keep, marker_keep)

    adata.obs["original_cell_id"] = adata.obs_names.astype(str)
    adata.obs["source_gse_id"] = gse_id
    adata.obs["phase1_annotation_keep"] = ann_keep
    adata.obs["phase1_marker_keep"] = marker_keep
    adata.obs["phase1_keep"] = keep_mask
    adata.obs["phase1_selection_reason"] = reasons
    adata.obs["phase1_t_hits"] = t_hits
    adata.obs["phase1_gd_hits"] = gd_hits
    adata.obs["phase1_nk_hits"] = nk_hits
    adata.obs["phase1_nk_strong_hits"] = nk_strong_hits
    adata.obs["phase1_t_score"] = t_score
    adata.obs["phase1_gd_score"] = gd_score
    adata.obs["phase1_nk_score"] = nk_score
    adata.obs["phase1_contam_score"] = contam_score
    adata.obs_names = pd.Index([f"{gse_id}__{cell_id}" for cell_id in adata.obs_names.astype(str)])

    subset = adata[keep_mask].copy()
    if not sp.issparse(subset.X):
        subset.X = sp.csr_matrix(np.asarray(subset.X))
    elif not sp.isspmatrix_csr(subset.X):
        subset.X = subset.X.tocsr()
    subset.obs = sanitize_dataframe(
        pd.DataFrame(
            {
                "original_cell_id": subset.obs["original_cell_id"].astype(object),
                "source_gse_id": subset.obs["source_gse_id"].astype(object),
                "phase1_donor_label": obs_label_series(subset.obs, donor_column),
                "phase1_sample_label": obs_label_series(subset.obs, sample_column),
                "phase1_library_label": obs_label_series(subset.obs, library_column),
                "phase1_annotation_label": obs_label_series(subset.obs, annotation_label_column),
                "phase1_annotation_keep": subset.obs["phase1_annotation_keep"].astype(bool),
                "phase1_marker_keep": subset.obs["phase1_marker_keep"].astype(bool),
                "phase1_keep": subset.obs["phase1_keep"].astype(bool),
                "phase1_selection_reason": subset.obs["phase1_selection_reason"].astype(object),
                "phase1_t_hits": subset.obs["phase1_t_hits"].astype(np.int16),
                "phase1_gd_hits": subset.obs["phase1_gd_hits"].astype(np.int16),
                "phase1_nk_hits": subset.obs["phase1_nk_hits"].astype(np.int16),
                "phase1_nk_strong_hits": subset.obs["phase1_nk_strong_hits"].astype(np.int16),
                "phase1_t_score": subset.obs["phase1_t_score"].astype(np.float32),
                "phase1_gd_score": subset.obs["phase1_gd_score"].astype(np.float32),
                "phase1_nk_score": subset.obs["phase1_nk_score"].astype(np.float32),
                "phase1_contam_score": subset.obs["phase1_contam_score"].astype(np.float32),
            },
            index=subset.obs_names,
        )
    )
    subset.var = pd.DataFrame(index=pd.Index(subset.var_names.astype(str), name=subset.var_names.name))
    subset.uns = {}
    subset.obsm = {}
    subset.varm = {}
    subset.obsp = {}
    subset.varp = {}

    summary = {
        "gse_id": gse_id,
        "input_n_obs": int(adata.n_obs),
        "input_n_vars": int(adata.n_vars),
        "candidate_n_obs": int(subset.n_obs),
        "candidate_fraction": float(subset.n_obs / adata.n_obs) if adata.n_obs else math.nan,
        "annotation_positive_n": int(ann_keep.sum()),
        "marker_positive_n": int(marker_keep.sum()),
        "annotation_and_marker_n": int((ann_keep & marker_keep).sum()),
        "annotation_only_n": int((ann_keep & ~marker_keep).sum()),
        "marker_only_n": int((~ann_keep & marker_keep).sum()),
        "annotation_columns": ";".join(annotation_columns),
        "mean_t_score_candidates": float(subset.obs["phase1_t_score"].mean()) if subset.n_obs else math.nan,
        "mean_nk_score_candidates": float(subset.obs["phase1_nk_score"].mean()) if subset.n_obs else math.nan,
        "mean_gd_score_candidates": float(subset.obs["phase1_gd_score"].mean()) if subset.n_obs else math.nan,
        "mean_contam_score_candidates": float(subset.obs["phase1_contam_score"].mean()) if subset.n_obs else math.nan,
    }

    marker_rows = [
        {"gse_id": gse_id, "marker_group": "T", "genes_present": int(t_present), "score_genes_present": int(t_score_present)},
        {"gse_id": gse_id, "marker_group": "GD", "genes_present": int(gd_present), "score_genes_present": int(gd_score_present)},
        {"gse_id": gse_id, "marker_group": "NK", "genes_present": int(nk_present), "score_genes_present": int(nk_score_present)},
        {"gse_id": gse_id, "marker_group": "NK_STRONG", "genes_present": int(nk_strong_present), "score_genes_present": int(nk_strong_present)},
        {"gse_id": gse_id, "marker_group": "CONTAM", "genes_present": int(contam_present), "score_genes_present": int(contam_present)},
    ]
    return subset, summary, marker_rows


def write_temp_subset(subset: ad.AnnData, gse_id: str) -> Path:
    """Write one dataset-specific candidate subset to a temp H5AD."""
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    tmp_path = TMP_DIR / f"{gse_id}_phase1_candidates.h5ad"
    if tmp_path.exists():
        tmp_path.unlink()
    subset.write_h5ad(tmp_path)
    return tmp_path


def concat_temp_subsets_on_disk(temp_files: dict[str, Path]) -> None:
    """Merge per-dataset candidate subsets into the canonical milestone H5AD."""
    if TNK_CANDIDATES_H5AD.exists():
        TNK_CANDIDATES_H5AD.unlink()
    ad.experimental.concat_on_disk(
        temp_files,
        TNK_CANDIDATES_H5AD,
        join="outer",
        merge=None,
        uns_merge=None,
        label="source_gse_id_concat",
        index_unique=None,
        fill_value=0,
        max_loaded_elems=50_000_000,
    )


def validate_merged_candidates(summary_df: pd.DataFrame) -> None:
    """Validate the merged candidate object against the per-dataset summary."""
    adata = ad.read_h5ad(TNK_CANDIDATES_H5AD, backed="r")
    try:
        expected_n_obs = int(summary_df["candidate_n_obs"].sum())
        if int(adata.n_obs) != expected_n_obs:
            raise ValueError(f"TNK_candidates.h5ad n_obs mismatch: {adata.n_obs} != {expected_n_obs}")

        if "source_gse_id" not in adata.obs.columns:
            raise ValueError("Merged candidate object is missing obs['source_gse_id']")

        observed = adata.obs["source_gse_id"].value_counts().sort_index()
        expected = summary_df.set_index("gse_id")["candidate_n_obs"].sort_index()
        aligned_index = expected.index.union(observed.index)
        observed = observed.reindex(aligned_index, fill_value=0).astype(int)
        expected = expected.reindex(aligned_index, fill_value=0).astype(int)
        if not observed.equals(expected):
            raise ValueError(
                "Per-dataset candidate counts in TNK_candidates.h5ad do not match the Phase 1 summary"
            )
    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()


def write_figures(summary_df: pd.DataFrame) -> None:
    """Write Phase 1 QC figures."""
    configure_plotting()

    plot_df = summary_df.sort_values("candidate_n_obs", ascending=False).copy()

    fig1, ax1 = plt.subplots(figsize=(9, 5.5))
    sns.barplot(data=plot_df, x="gse_id", y="candidate_n_obs", color="#2f6c8f", ax=ax1)
    ax1.set_title("Phase 1 Candidate Yield by Dataset")
    ax1.set_xlabel("GSE ID")
    ax1.set_ylabel("Retained candidate cells")
    ax1.tick_params(axis="x", rotation=35)
    fig1.tight_layout()
    fig1.savefig(FIGURE_DIR / "phase1_categoryA_candidate_yield.png")
    plt.close(fig1)

    fig2, ax2 = plt.subplots(figsize=(9, 5.5))
    sns.barplot(data=plot_df, x="gse_id", y="candidate_fraction", color="#4b9e5f", ax=ax2)
    ax2.set_title("Phase 1 Candidate Fraction by Dataset")
    ax2.set_xlabel("GSE ID")
    ax2.set_ylabel("Retained fraction")
    ax2.tick_params(axis="x", rotation=35)
    fig2.tight_layout()
    fig2.savefig(FIGURE_DIR / "phase1_categoryA_candidate_fraction.png")
    plt.close(fig2)

    score_plot = summary_df[["gse_id", "mean_t_score_candidates", "mean_nk_score_candidates", "mean_gd_score_candidates"]].copy()
    score_plot = score_plot.melt(id_vars="gse_id", var_name="score_type", value_name="mean_score")
    fig3, ax3 = plt.subplots(figsize=(10, 5.5))
    sns.barplot(data=score_plot, x="gse_id", y="mean_score", hue="score_type", ax=ax3)
    ax3.set_title("Phase 1 Mean Marker Support in Retained Candidates")
    ax3.set_xlabel("GSE ID")
    ax3.set_ylabel("Mean score")
    ax3.tick_params(axis="x", rotation=35)
    fig3.tight_layout()
    fig3.savefig(FIGURE_DIR / "phase1_categoryA_marker_support.png")
    plt.close(fig3)


def write_qc_summary(summary_df: pd.DataFrame) -> None:
    """Write a concise Phase 1 QC summary."""
    total_input = int(summary_df["input_n_obs"].sum())
    total_candidates = int(summary_df["candidate_n_obs"].sum())
    total_annotation = int(summary_df["annotation_positive_n"].sum())
    total_marker = int(summary_df["marker_positive_n"].sum())

    with PHASE1_QC_MD.open("w", encoding="utf-8") as handle:
        handle.write("# Phase 1 QC Summary\n\n")
        handle.write("## Scope\n\n")
        handle.write("- Input set: Category A datasets only\n")
        handle.write(f"- Datasets processed: {summary_df.shape[0]}\n")
        handle.write(f"- Total input cells: {total_input}\n")
        handle.write(f"- Total retained candidate cells: {total_candidates}\n")
        handle.write(
            f"- Overall retained fraction: {total_candidates / total_input:.4f}\n\n"
            if total_input
            else "- Overall retained fraction: NA\n\n"
        )

        handle.write("## Selection Signals\n\n")
        handle.write(f"- Cells retained by annotation support: {total_annotation}\n")
        handle.write(f"- Cells retained by marker support: {total_marker}\n")
        handle.write(
            f"- Datasets with annotation columns used: {(summary_df['annotation_columns'] != '').sum()}\n\n"
        )

        handle.write("## Per-dataset yield\n\n")
        for _, row in summary_df.sort_values("candidate_n_obs", ascending=False).iterrows():
            handle.write(
                f"- {row['gse_id']}: retained {int(row['candidate_n_obs'])} / {int(row['input_n_obs'])} "
                f"({row['candidate_fraction']:.4f}); annotation columns=`{row['annotation_columns'] or 'none'}`\n"
            )

        handle.write("\n## QC Conclusion\n\n")
        handle.write(
            "Phase 1 Category A coarse extraction is complete. Review TNK_candidates.h5ad, "
            "the per-dataset selection table, and the Phase 1 figures before deciding whether "
            "to proceed to Phase 1b or adjust thresholds.\n"
        )


def main() -> None:
    """Run Phase 1 coarse extraction for the approved Category A datasets."""
    ensure_output_dirs()
    registry = load_category_a_registry()

    summary_rows: list[dict] = []
    marker_rows: list[dict] = []
    temp_files: dict[str, Path] = {}

    for _, row in registry.iterrows():
        gse_id = row["gse_id"]
        print(f"Processing {gse_id} ...", flush=True)
        subset, summary, markers = extract_one_dataset(gse_id, row["h5ad_path"])
        temp_files[gse_id] = write_temp_subset(subset, gse_id)
        summary_rows.append(summary)
        marker_rows.extend(markers)

    summary_df = pd.DataFrame(summary_rows).sort_values("gse_id").reset_index(drop=True)
    marker_df = pd.DataFrame(marker_rows).sort_values(["gse_id", "marker_group"]).reset_index(drop=True)

    concat_temp_subsets_on_disk(temp_files)
    validate_merged_candidates(summary_df)

    summary_df.to_csv(PHASE1_SUMMARY_CSV, index=False)
    marker_df.to_csv(PHASE1_MARKER_CSV, index=False)
    write_figures(summary_df)
    write_qc_summary(summary_df)
    shutil.rmtree(TMP_DIR, ignore_errors=True)

    print("\nWrote:")
    print(TNK_CANDIDATES_H5AD)
    print(PHASE1_SUMMARY_CSV)
    print(PHASE1_MARKER_CSV)
    print(PHASE1_QC_MD)
    for figure_name in [
        "phase1_categoryA_candidate_yield.png",
        "phase1_categoryA_candidate_fraction.png",
        "phase1_categoryA_marker_support.png",
    ]:
        print(FIGURE_DIR / figure_name)


if __name__ == "__main__":
    main()
