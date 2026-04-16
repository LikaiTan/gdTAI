#!/usr/bin/env python3
"""Build the parallel plus6 integrated milestone.

This pipeline keeps the current validated 25-GSE milestone untouched and builds
an independent `plus6` lineage on the mirrored SSD tree. It handles:

1. input compatibility prep for the six incoming datasets
2. merged prephase3 construction
3. fresh scVI integration + Leiden + UMAP
4. Phase 4 rescoring on the merged object
5. simple marker/TCR-based annotation
6. HTML/PDF report generation and optional email delivery
"""

from __future__ import annotations

import argparse
import html
import logging
import math
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Iterable

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata._io.specs import write_elem
from scipy import sparse
from sklearn.utils import sparsefuncs
from pandas.api.types import CategoricalDtype

from phase4_gdt_module_scoring import (
    CHUNK_SIZE,
    CTRL_SIZE,
    FIGURE_DPI,
    N_BINS,
    PACKAGE_ZIP,
    PHASE4_SCALED_SCORE_COLUMNS,
    PHASE4_SCORE_COLUMNS,
    RANDOM_STATE,
    SCATTER_PLOT_SAMPLE_SIZE,
    TARGET_SUM,
    add_scaled_scores,
    append_obs_columns_in_place,
    compute_gene_means,
    compute_scores,
    extract_log1p_gene_expression_for_sample,
    find_module_genes,
    pick_control_genes,
)


PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
TABLE_DIR = OUTPUT_ROOT / "tables" / "plus6"
FIGURE_DIR = OUTPUT_ROOT / "figures" / "plus6"
LOG_DIR = OUTPUT_ROOT / "logs" / "plus6"
MODEL_DIR = OUTPUT_ROOT / "models" / "phase3_plus6_scvi"

BASE_INTEGRATED = PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "integrated.h5ad"
PLUS6_PREP_DIR = PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "plus6_inputs"
PLUS6_SOURCE_COPY_DIR = PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "plus6_source_inputs"
PLUS6_PREPHASE3 = PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "integrated_plus6_prephase3.h5ad"
PLUS6_INTEGRATED = PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "integrated_plus6.h5ad"
PLUS6_REPORT_MD = OUTPUT_ROOT / "plus6_profile_report.md"
PLUS6_REPORT_HTML = OUTPUT_ROOT / "plus6_profile_report.html"
PLUS6_REPORT_PDF = OUTPUT_ROOT / "plus6_profile_report.pdf"

LOG_FILE = LOG_DIR / "plus6_parallel_pipeline.log"
PREP_SUMMARY_CSV = TABLE_DIR / "plus6_input_compatibility_summary.csv"
PREP_QC_MD = LOG_DIR / "plus6_input_compatibility_summary.md"
PHASE3_SUMMARY_CSV = TABLE_DIR / "plus6_phase3_summary.csv"
PHASE3_HVG_FILTER_CSV = TABLE_DIR / "plus6_phase3_hvg_filter_summary.csv"
PHASE3_HVG_GENES_TXT = TABLE_DIR / "plus6_phase3_hvg_genes.txt"
PHASE3_LEIDEN_COUNTS_CSV = TABLE_DIR / "plus6_phase3_leiden_counts.csv"
PHASE4_SCORE_SUMMARY_CSV = TABLE_DIR / "plus6_phase4_score_summary.csv"
PHASE4_GSE_SUMMARY_CSV = TABLE_DIR / "plus6_phase4_gse_summary.csv"
PHASE4_TISSUE_SUMMARY_CSV = TABLE_DIR / "plus6_phase4_tissue_summary.csv"
PHASE4_TRD_OVER_TRAB_BY_TISSUE_CSV = TABLE_DIR / "plus6_phase4_trd_over_trab_by_tissue.csv"
PHASE4_TRD_OVER_TRAB_BY_GSE_CSV = TABLE_DIR / "plus6_phase4_trd_over_trab_by_gse.csv"
PHASE4_TRD_GT_0P1_BY_TISSUE_CSV = TABLE_DIR / "plus6_phase4_trd_gt_0p1_by_tissue.csv"
PHASE4_TRD_GT_0P1_BY_GSE_CSV = TABLE_DIR / "plus6_phase4_trd_gt_0p1_by_gse.csv"
ANNOTATION_CLUSTER_SUMMARY_CSV = TABLE_DIR / "plus6_simple_annotation_cluster_summary.csv"
ANNOTATION_CHANGES_CSV = TABLE_DIR / "plus6_simple_annotation_changes.csv"
ANNOTATION_GSE_COUNTS_CSV = TABLE_DIR / "plus6_simple_annotation_by_gse.csv"
ANNOTATION_TISSUE_COUNTS_CSV = TABLE_DIR / "plus6_simple_annotation_by_tissue.csv"
ANNOTATION_SUMMARY_MD = LOG_DIR / "plus6_simple_annotation_summary.md"

LEIDEN_RESOLUTION = 1.0
HVG_N = 4000
SCVI_MAX_EPOCHS = 200
SCVI_BATCH_SIZE = 8192
LEGEND_SAMPLE_LIMIT = 300_000

NEW_INPUTS = {
    "HRA005041": PROJECT_ROOT / "downloads" / "per_gse_h5ad_with_metadata" / "HRA005041_T_cells_subset.h5ad",
    "GDT_2020AUG_woCOV": PROJECT_ROOT / "newdata" / "Sorted_gdT" / "GDT_2020AUG_woCOV_sorted_gdt.h5ad",
    "GDTlung2023july_7p": PROJECT_ROOT / "newdata" / "Sorted_gdT" / "GDTlung2023july_7p_sorted_gdt.h5ad",
    "MalteGDT": PROJECT_ROOT / "newdata" / "Sorted_gdT" / "MalteGDT_sorted_gdt.h5ad",
    "GSE144469": PROJECT_ROOT / "newdata" / "GSE144469_tnk_subset.h5ad",
    "GSE206325": PROJECT_ROOT / "newdata" / "GSE206325_integrated.h5ad",
}
SSD_LOCAL_SOURCE_ALIASES = {"HRA005041", "GSE144469", "GSE206325"}
PREP_COUNT_MATRIX_SOURCE_HINTS = {
    "current_integrated_base": "X",
    "HRA005041": "pseudo_counts_from_seurat_lognormalize",
    "GDT_2020AUG_woCOV": "X_integer_like",
    "GDTlung2023july_7p": "X_integer_like",
    "MalteGDT": "X_integer_like",
    "GSE144469": "layers[counts]",
    "GSE206325": "layers[counts]",
}

BASE_KEEP_COLUMNS = [
    "original_cell_id",
    "source_gse_id",
    "project name",
    "sampleid",
    "sample_id",
    "library_id",
    "barcodes",
    "barcode",
    "barcode_core",
    "donor_id",
    "tissue",
    "tissue_corrected",
    "condition",
    "technology_simple",
    "assay_type",
    "input_population",
    "Sorted_gdT",
    "tcr_chain_mode",
    "TCRseq",
    "TRA_cdr3",
    "TRA_v",
    "TRA_d",
    "TRA_j",
    "TRA_cdr3_nt",
    "TRA_clone_id",
    "TRA_umis",
    "TRA_reads",
    "TRA_c_gene",
    "TRB_cdr3",
    "TRB_v",
    "TRB_d",
    "TRB_j",
    "TRB_cdr3_nt",
    "TRB_clone_id",
    "TRB_umis",
    "TRB_reads",
    "TRB_c_gene",
    "TRG_cdr3",
    "TRG_v",
    "TRG_d",
    "TRG_j",
    "TRG_cdr3_nt",
    "TRG_clone_id",
    "TRG_umis",
    "TRG_reads",
    "TRG_c_gene",
    "TRD_cdr3",
    "TRD_v",
    "TRD_d",
    "TRD_j",
    "TRD_cdr3_nt",
    "TRD_clone_id",
    "TRD_umis",
    "TRD_reads",
    "TRD_c_gene",
]

PHASE4_OLD_COLUMNS = [
    "phase4_tra_score",
    "phase4_trb_score",
    "phase4_trab_score",
    "phase4_trd_score",
    "phase4_trd_minus_trab",
    "phase4_trab_score_scaled",
    "phase4_trd_score_scaled",
    "phase4_trd_minus_trab_scaled",
]

OLD_INTEGRATION_COLUMNS = [
    "phase3_batch_key",
    "phase3_batch_level",
    "leiden",
    "simple_annotation",
    "simple_annotation_plus6",
    "scanvi_detailed_label",
    "scanvi_label_confidence",
    "scanvi_tnk_superclass",
    "scanvi_reference_other_flag",
    "scanvi_transfer_method",
]

ANNOTATION_GENES = [
    "TRDC",
    "TRDV1",
    "TRDV2",
    "TRGV9",
    "FOXP3",
    "IL2RA",
    "CTLA4",
    "TIGIT",
    "IKZF2",
    "NKG7",
    "GNLY",
    "KLRD1",
    "FCER1G",
    "PRF1",
    "CTSW",
    "CD3D",
    "CD3E",
    "TRAC",
    "TRBC1",
    "TRBC2",
    "CD4",
    "IL7R",
    "LTB",
    "MAL",
    "CD8A",
    "CD8B",
    "CCL5",
]

GD_ANNOTATION_ALLOWED_GENES = ["TRDC", "TRDV1", "TRDV2", "TRGV9"]
TREG_DIAGNOSTIC_GENES = [
    "FOXP3",
    "IL2RA",
    "CTLA4",
    "CD4",
    "IL7R",
    "LTB",
    "CD8A",
    "CD8B",
    "NKG7",
    "GNLY",
]
TREG_CLUSTER_DOTPLOT_GENES = [
    "FOXP3",
    "IL2RA",
    "CTLA4",
    "CD4",
    "IL7R",
    "LTB",
    "CD8A",
    "CD8B",
    "CCL5",
    "NKG7",
    "GNLY",
    "KLRD1",
]

MISSING_STRING = ""
WRITE_SAFE_TEXT_COLUMNS = [
    "source_gse_id",
    "original_cell_id",
    "project name",
    "sampleid",
    "barcodes",
    "sample_id",
    "library_id",
    "donor_id",
    "tissue",
    "tissue_corrected",
    "technology_simple",
    "assay_type",
    "input_population",
    "Sorted_gdT",
    "tcr_chain_mode",
    "TCRseq",
    "simple_annotation_plus6",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the plus6 parallel integrated milestone.")
    parser.add_argument(
        "--stage",
        choices=["prepare", "integrate", "score", "annotate", "report", "all"],
        default="all",
        help="Pipeline stage to run.",
    )
    parser.add_argument(
        "--send-email",
        action="store_true",
        help="Send the final PDF report to Likai during the report stage.",
    )
    return parser.parse_args()


def configure_logging() -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    PLUS6_PREP_DIR.mkdir(parents=True, exist_ok=True)
    PLUS6_SOURCE_COPY_DIR.mkdir(parents=True, exist_ok=True)
    handlers = [
        logging.FileHandler(LOG_FILE, mode="a", encoding="utf-8"),
        logging.StreamHandler(),
    ]
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
        force=True,
    )


def normalize_text_series(series: pd.Series) -> pd.Series:
    series = series.astype("string")
    series = series.fillna("")
    series = series.str.strip()
    series = series.mask(series.eq(""), pd.NA)
    return series


def make_obs_write_safe(adata: ad.AnnData) -> None:
    for column in WRITE_SAFE_TEXT_COLUMNS:
        if column in adata.obs.columns:
            normalized = normalize_text_series(pd.Series(adata.obs[column], index=adata.obs_names))
            adata.obs[column] = pd.Series(
                np.where(normalized.isna(), "", normalized.astype(str)),
                index=adata.obs_names,
                dtype=object,
            )

    for column in adata.obs.columns:
        series = adata.obs[column]
        if isinstance(series.dtype, CategoricalDtype) or pd.api.types.is_string_dtype(series.dtype):
            normalized = normalize_text_series(pd.Series(series, index=adata.obs_names))
            adata.obs[column] = pd.Series(
                np.where(normalized.isna(), "", normalized.astype(str)),
                index=adata.obs_names,
                dtype=object,
            )


def make_var_write_safe(adata: ad.AnnData) -> None:
    adata.var_names = pd.Index(adata.var_names.astype(str), dtype=object)
    adata.var.index = adata.var_names
    for column in adata.var.columns:
        series = adata.var[column]
        if isinstance(series.dtype, CategoricalDtype) or pd.api.types.is_string_dtype(series.dtype):
            normalized = normalize_text_series(pd.Series(series, index=adata.var_names))
            adata.var[column] = pd.Series(
                np.where(normalized.isna(), "", normalized.astype(str)),
                index=adata.var_names,
                dtype=object,
            )


def build_hvg_exclusion_frame(var_names: pd.Index) -> pd.DataFrame:
    symbols = pd.Series(var_names.astype(str), index=var_names.astype(str), dtype="string")
    upper = symbols.str.upper()

    mitochondrial = upper.str.startswith("MT-")
    ribosomal = upper.str.startswith(("RPS", "RPL", "MRPS", "MRPL"))
    noncoding_named = upper.str.startswith(
        (
            "LINC",
            "MIR",
            "MIRLET",
            "SNHG",
            "SNORA",
            "SNORD",
            "SCARNA",
            "RNU",
            "RN7",
            "RNA5",
            "RNA18",
            "RNA28",
            "YRNA",
            "Y_RNA",
            "RMRP",
            "RPPH1",
            "XIST",
            "TSIX",
            "NEAT1",
            "MALAT1",
        )
    )
    noncoding_clone = upper.str.match(r"^(AC[0-9]+|AL[0-9]+|AP[0-9]+|BX[0-9]+|RP[0-9]+-|CTA-|CTB-|CTC-|CTD-|CTE-|CU-|DKFZP)")
    noncoding_like = noncoding_named | noncoding_clone
    excluded = mitochondrial | ribosomal | noncoding_like

    frame = pd.DataFrame(
        {
            "gene": symbols.to_numpy(),
            "exclude_from_hvg": excluded.to_numpy(),
            "exclude_mitochondrial": mitochondrial.to_numpy(),
            "exclude_ribosomal": ribosomal.to_numpy(),
            "exclude_noncoding_like": noncoding_like.to_numpy(),
        }
    )
    frame["exclude_reason"] = np.select(
        [
            frame["exclude_mitochondrial"] & frame["exclude_ribosomal"] & frame["exclude_noncoding_like"],
            frame["exclude_mitochondrial"] & frame["exclude_ribosomal"],
            frame["exclude_mitochondrial"] & frame["exclude_noncoding_like"],
            frame["exclude_ribosomal"] & frame["exclude_noncoding_like"],
            frame["exclude_mitochondrial"],
            frame["exclude_ribosomal"],
            frame["exclude_noncoding_like"],
        ],
        [
            "mitochondrial,ribosomal,noncoding_like",
            "mitochondrial,ribosomal",
            "mitochondrial,noncoding_like",
            "ribosomal,noncoding_like",
            "mitochondrial",
            "ribosomal",
            "noncoding_like",
        ],
        default="",
    )
    return frame


def text_from_obs(obs: pd.DataFrame, candidates: Iterable[str], default: str = MISSING_STRING) -> pd.Series:
    for column in candidates:
        if column in obs.columns:
            return normalize_text_series(obs[column])
    return pd.Series(default, index=obs.index, dtype="string")


def ensure_string_column(obs: pd.DataFrame, column: str, values: pd.Series | str) -> None:
    if isinstance(values, pd.Series):
        obs[column] = normalize_text_series(values).fillna(MISSING_STRING).astype(str)
    else:
        obs[column] = np.repeat(str(values), obs.shape[0])


def derive_barcodes_from_obs_names(obs_names: pd.Index) -> pd.Series:
    values = []
    for obs_name in obs_names.astype(str):
        token = obs_name.split("__")[-1]
        token = token.split("_")[-1]
        values.append(token)
    return pd.Series(values, index=obs_names, dtype="string")


def derive_barcode_core(barcodes: pd.Series) -> pd.Series:
    stripped = normalize_text_series(barcodes).fillna("")
    stripped = stripped.str.replace(r"-\d+$", "", regex=True)
    stripped = stripped.str.replace(r"_\d+$", "", regex=True)
    stripped = stripped.str.replace(r"[^ACGTN]", "", regex=True)
    return stripped.astype("string")


def bool_from_text(values: pd.Series, truthy: set[str] | None = None) -> pd.Series:
    truthy = truthy or {"true", "1", "yes", "y", "productive", "paired"}
    normalized = values.astype(str).str.strip().str.lower()
    return normalized.isin(truthy)


def choose_var_names(adata: ad.AnnData, base_genes: set[str]) -> None:
    if adata.var_names is not None and len(adata.var_names):
        current_overlap = sum(g in base_genes for g in adata.var_names.astype(str))
    else:
        current_overlap = -1

    best_names = pd.Index(adata.var_names.astype(str), dtype="string")
    best_overlap = current_overlap
    for candidate in ["gene_symbol", "feature_name", "_index"]:
        if candidate in adata.var.columns:
            series = normalize_text_series(pd.Series(adata.var[candidate].astype(str), index=adata.var.index))
            series = series.fillna("")
            overlap = int(series.isin(base_genes).sum())
            if overlap > best_overlap:
                best_overlap = overlap
                best_names = pd.Index(series.to_numpy(), dtype="string")

    adata.var_names = best_names.astype(str)
    adata.var_names_make_unique()


def move_counts_to_x(adata: ad.AnnData, alias: str) -> str:
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
        return "layers[counts]"
    x = adata.X
    if sparse.issparse(x):
        data = x.data
        if data.size == 0:
            return "X_empty"
        integer_like = np.allclose(data[: min(100000, data.size)], np.round(data[: min(100000, data.size)]))
        return "X_integer_like" if integer_like else "X_noninteger_like"
    arr = np.asarray(x)
    integer_like = np.allclose(arr[: min(1000, arr.shape[0]), : min(1000, arr.shape[1])], np.round(arr[: min(1000, arr.shape[0]), : min(1000, arr.shape[1])]))
    if not integer_like:
        logging.warning("%s lacks an explicit counts layer and X looks normalized; keeping X as-is for merge.", alias)
    return "X"


def infer_seurat_pseudocounts(adata: ad.AnnData, alias: str) -> str:
    if alias != "HRA005041":
        return ""
    if "nCount_RNA" not in adata.obs.columns and "total_counts" not in adata.obs.columns:
        return ""
    x = adata.X.tocsr() if sparse.issparse(adata.X) else sparse.csr_matrix(np.asarray(adata.X))
    probe_rows = min(1000, x.shape[0])
    probe = x[:probe_rows].copy()
    probe.data = np.expm1(probe.data)
    row_sums = np.asarray(probe.sum(axis=1)).ravel()
    median_sum = float(np.median(row_sums)) if row_sums.size else 0.0
    if not (9500.0 <= median_sum <= 10500.0):
        logging.warning(
            "%s did not match Seurat LogNormalize expectation closely enough for pseudo-count inversion (median expm1 row sum %.2f)",
            alias,
            median_sum,
        )
        return ""

    ncount_col = "nCount_RNA" if "nCount_RNA" in adata.obs.columns else "total_counts"
    ncount = pd.to_numeric(adata.obs[ncount_col], errors="coerce").to_numpy(dtype=np.float32)
    if np.isnan(ncount).any():
        raise ValueError(f"{alias} has missing {ncount_col} values; cannot infer pseudo-counts.")
    logging.info(
        "Inferring pseudo-counts for %s from Seurat LogNormalize using %s and scale factor 10000",
        alias,
        ncount_col,
    )
    x = x.copy()
    x.data = np.expm1(x.data)
    x = x.multiply(ncount[:, None] / 1e4)
    x.data = np.rint(x.data).astype(np.float32, copy=False)
    x.eliminate_zeros()
    adata.X = x.tocsr()
    adata.uns["plus6_pseudocounts_inferred"] = {
        "method": "seurat_lognormalize_inverse",
        "scale_factor": 10000.0,
        "ncount_column": ncount_col,
        "probe_row_expm1_sum_median": median_sum,
        "warning": "pseudo-counts inferred from normalized data; acceptable for integration sensitivity, not downstream DE truth",
    }
    return "pseudo_counts_from_seurat_lognormalize"


def drop_stale_columns(obs: pd.DataFrame) -> pd.DataFrame:
    drop_cols = [column for column in PHASE4_OLD_COLUMNS + OLD_INTEGRATION_COLUMNS if column in obs.columns]
    if drop_cols:
        obs = obs.drop(columns=drop_cols)
    return obs


def blank_series(index: pd.Index) -> pd.Series:
    return pd.Series(np.repeat(MISSING_STRING, len(index)), index=index, dtype="string")


def harmonize_common_tcr_fields(obs: pd.DataFrame) -> pd.DataFrame:
    for chain in ("TRA", "TRB", "TRG", "TRD"):
        for suffix in ("cdr3", "v", "d", "j", "cdr3_nt", "clone_id", "umis", "reads", "c_gene"):
            column = f"{chain}_{suffix}"
            if column not in obs.columns:
                obs[column] = blank_series(obs.index)
            else:
                obs[column] = normalize_text_series(obs[column]).fillna(MISSING_STRING).astype(str)

    has_tra = obs["TRA_cdr3"].ne("")
    has_trb = obs["TRB_cdr3"].ne("")
    has_trg = obs["TRG_cdr3"].ne("")
    has_trd = obs["TRD_cdr3"].ne("")

    if "has_TRA" in obs.columns:
        has_tra = has_tra | bool_from_text(normalize_text_series(obs["has_TRA"]).fillna("false"))
    if "has_TRB" in obs.columns:
        has_trb = has_trb | bool_from_text(normalize_text_series(obs["has_TRB"]).fillna("false"))
    if "has_TRG" in obs.columns:
        has_trg = has_trg | bool_from_text(normalize_text_series(obs["has_TRG"]).fillna("false"))
    if "has_TRD" in obs.columns:
        has_trd = has_trd | bool_from_text(normalize_text_series(obs["has_TRD"]).fillna("false"))

    obs["has_TRA_TRB_paired"] = pd.Series((has_tra & has_trb).to_numpy(dtype=bool), index=obs.index)
    obs["has_TRG_TRD_paired"] = pd.Series((has_trg & has_trd).to_numpy(dtype=bool), index=obs.index)
    obs["has_any_ab_tcr"] = pd.Series((has_tra | has_trb).to_numpy(dtype=bool), index=obs.index)
    obs["has_any_gd_tcr"] = pd.Series((has_trg | has_trd).to_numpy(dtype=bool), index=obs.index)

    if "TCRseq" not in obs.columns:
        obs["TCRseq"] = np.where(obs["has_any_ab_tcr"] | obs["has_any_gd_tcr"], "yes", "no")
    else:
        normalized = normalize_text_series(obs["TCRseq"]).fillna("")
        yes_mask = bool_from_text(normalized) | normalized.str.lower().isin({"yes", "detected", "available"})
        obs["TCRseq"] = np.where(yes_mask | obs["has_any_ab_tcr"] | obs["has_any_gd_tcr"], "yes", "no")

    return obs


def standardize_obs_names(obs: pd.DataFrame, source_gse_id: str) -> pd.Index:
    original = normalize_text_series(text_from_obs(obs, ["original_cell_id", "cell_id", "obs_name", "raw_obs_name"]))
    if original.isna().all():
        original = pd.Series(obs.index.astype(str), index=obs.index, dtype="string")
    original = original.fillna(pd.Series(obs.index.astype(str), index=obs.index, dtype="string"))
    obs["original_cell_id"] = original.astype(str)
    prefixed = pd.Index([f"{source_gse_id}__{value}" for value in obs["original_cell_id"].astype(str)], dtype=object)
    if prefixed.has_duplicates:
        dup_counts = pd.Series(prefixed.astype(str)).groupby(prefixed.astype(str)).cumcount()
        prefixed = pd.Index(
            [
                f"{name}__dup{count}" if count > 0 else name
                for name, count in zip(prefixed.astype(str), dup_counts.astype(int), strict=False)
            ],
            dtype=object,
        )
    return prefixed


def standardize_generic_metadata(adata: ad.AnnData, alias: str) -> tuple[ad.AnnData, dict[str, object]]:
    obs = adata.obs.copy()
    obs = drop_stale_columns(obs)

    source_gse_id = text_from_obs(obs, ["source_gse_id", "gse_id", "GSE"], default=alias).fillna(alias)
    nonempty_gse = source_gse_id.replace("", pd.NA).dropna()
    if nonempty_gse.empty:
        source_gse_id = pd.Series(np.repeat(alias, obs.shape[0]), index=obs.index, dtype="string")
        multi_source_gse = False
        source_gse_id_value = alias
    else:
        multi_source_gse = nonempty_gse.nunique() > 1
        source_gse_id_value = str(nonempty_gse.iloc[0]) if not multi_source_gse else alias

    project_name = text_from_obs(obs, ["project name", "Project", "source_name", "reference"], default=source_gse_id_value)
    sample_id = text_from_obs(obs, ["sample_id", "sampleid", "sample_ID", "sample", "orig.ident", "cell_to_sample_ID", "Sample Name"], default=source_gse_id_value)
    sampleid = text_from_obs(obs, ["sampleid", "sample_id", "sample_ID", "sample", "orig.ident", "cell_to_sample_ID", "Sample Name"], default=sample_id.iloc[0] if len(sample_id) else source_gse_id_value)
    library_id = text_from_obs(obs, ["library_id", "Sample Name", "sample", "sample_id", "sampleid", "orig.ident", "batch"], default=sample_id.iloc[0] if len(sample_id) else source_gse_id_value)
    donor_id = text_from_obs(obs, ["donor_id", "donor", "Donor", "patient_ID", "patient_id", "patient_ID", "patient"], default=MISSING_STRING)
    tissue = text_from_obs(obs, ["tissue_corrected", "tissue", "Tissue", "Tissue_merged", "tissue_region", "source_name"], default=MISSING_STRING)
    tissue_corrected = text_from_obs(obs, ["tissue_corrected", "tissue", "Tissue", "Tissue_merged", "tissue_region"], default=MISSING_STRING)
    condition = text_from_obs(obs, ["condition", "colon_status", "disease", "group", "Stage", "treatment"], default=MISSING_STRING)
    assay_type = text_from_obs(obs, ["assay_type"], default="scRNA")
    technology_simple = text_from_obs(obs, ["technology_simple", "library_chemistry"], default="10x_scRNA")
    input_population = text_from_obs(obs, ["input_population"], default="tnk_subset")
    sorted_gdt = text_from_obs(obs, ["Sorted_gdT"], default="False")
    tcr_chain_mode = text_from_obs(obs, ["tcr_chain_mode"], default="none")

    barcodes = text_from_obs(obs, ["barcodes", "barcode", "cell_barcode", "cellbarcode"], default=MISSING_STRING)
    missing_barcodes = barcodes.fillna("").eq("")
    if missing_barcodes.any():
        derived = derive_barcodes_from_obs_names(obs.index)
        barcodes = barcodes.mask(missing_barcodes, derived)
    barcode = text_from_obs(obs, ["barcode", "barcodes", "cell_barcode", "cellbarcode"], default=MISSING_STRING)
    missing_barcode = barcode.fillna("").eq("")
    if missing_barcode.any():
        barcode = barcode.mask(missing_barcode, barcodes)
    barcode_core = text_from_obs(obs, ["barcode_core"], default=MISSING_STRING)
    missing_core = barcode_core.fillna("").eq("")
    if missing_core.any():
        barcode_core = barcode_core.mask(missing_core, derive_barcode_core(barcodes))

    ensure_string_column(obs, "source_gse_id", source_gse_id)
    ensure_string_column(obs, "project name", project_name)
    ensure_string_column(obs, "sample_id", sample_id)
    ensure_string_column(obs, "sampleid", sampleid)
    ensure_string_column(obs, "library_id", library_id)
    ensure_string_column(obs, "donor_id", donor_id)
    ensure_string_column(obs, "tissue", tissue)
    ensure_string_column(obs, "tissue_corrected", tissue_corrected)
    ensure_string_column(obs, "condition", condition)
    ensure_string_column(obs, "assay_type", assay_type)
    ensure_string_column(obs, "technology_simple", technology_simple)
    ensure_string_column(obs, "input_population", input_population)
    ensure_string_column(obs, "Sorted_gdT", sorted_gdt)
    ensure_string_column(obs, "tcr_chain_mode", tcr_chain_mode)
    ensure_string_column(obs, "barcodes", barcodes)
    ensure_string_column(obs, "barcode", barcode)
    ensure_string_column(obs, "barcode_core", barcode_core)

    obs = harmonize_common_tcr_fields(obs)
    if multi_source_gse:
        original = normalize_text_series(text_from_obs(obs, ["original_cell_id", "cell_id", "obs_name", "raw_obs_name"]))
        if original.isna().all():
            original = pd.Series(obs.index.astype(str), index=obs.index, dtype="string")
        original = original.fillna(pd.Series(obs.index.astype(str), index=obs.index, dtype="string"))
        obs["original_cell_id"] = original.astype(str)
        obs.index = pd.Index(obs.index.astype(str), dtype=object)
    else:
        obs.index = standardize_obs_names(obs, source_gse_id_value)

    # Dataset-specific repairs.
    if alias.startswith("GDT") or alias == "MalteGDT":
        obs["Sorted_gdT"] = "True"
        obs["input_population"] = "purified_gdt"
        if obs["tcr_chain_mode"].eq("").all():
            obs["tcr_chain_mode"] = "gd_only"

    if alias == "HRA005041":
        if obs["input_population"].eq("").all():
            obs["input_population"] = "tnk_subset"
        if obs["tcr_chain_mode"].eq("").all():
            obs["tcr_chain_mode"] = "ab_and_gd"

    if alias == "GSE144469":
        obs["source_gse_id"] = "GSE144469"
        obs["project name"] = "GSE144469"
        obs["input_population"] = "tnk_subset"
        obs["technology_simple"] = "10x_scRNA"
        obs["tcr_chain_mode"] = "partial_metadata"
        if "TRD_cdr3" not in obs.columns:
            obs["TRD_cdr3"] = blank_series(obs.index)
        if "TRD_v" in adata.obs.columns:
            obs["TRD_v"] = normalize_text_series(adata.obs["TRD_v"]).fillna(MISSING_STRING).astype(str)
        if "TRG_v" in adata.obs.columns:
            obs["TRG_v"] = normalize_text_series(adata.obs["TRG_v"]).fillna(MISSING_STRING).astype(str)

    if alias == "GSE206325":
        obs["source_gse_id"] = "GSE206325"
        if obs["project name"].eq("").all():
            obs["project name"] = "GSE206325"
        if obs["input_population"].eq("").all():
            obs["input_population"] = "tnk_subset"
        if obs["technology_simple"].eq("").all():
            obs["technology_simple"] = "10x_scRNA"

    keep_columns = BASE_KEEP_COLUMNS + [
        "has_TRA_TRB_paired",
        "has_TRG_TRD_paired",
        "has_any_ab_tcr",
        "has_any_gd_tcr",
    ]
    for column in keep_columns:
        if column not in obs.columns:
            if column.startswith("has_"):
                obs[column] = False
            else:
                obs[column] = MISSING_STRING
    obs = obs.loc[:, keep_columns].copy()

    adata.obs = obs
    adata.obsm.clear()
    adata.varm.clear()
    adata.obsp.clear()
    adata.varp.clear()
    adata.uns.clear()
    adata.raw = None

    prep_stats = {
        "alias": alias,
        "source_gse_id": source_gse_id_value if not multi_source_gse else "multiple_existing_gses",
        "n_cells": int(adata.n_obs),
        "n_genes_before_align": int(adata.n_vars),
        "n_samples": int(pd.Series(obs["sample_id"]).replace("", pd.NA).dropna().nunique()),
        "n_paired_ab": int(obs["has_TRA_TRB_paired"].sum()),
        "n_paired_gd": int(obs["has_TRG_TRD_paired"].sum()),
        "sorted_gdt_flag": bool(pd.Series(obs["Sorted_gdT"]).astype(str).str.lower().eq("true").all()),
    }
    return adata, prep_stats


def subset_to_base_genes(adata: ad.AnnData, base_var_names: pd.Index) -> tuple[ad.AnnData, int]:
    current = pd.Index(adata.var_names.astype(str), dtype="string")
    shared = base_var_names.intersection(current)
    missing_count = int(base_var_names.size - shared.size)
    adata = adata[:, shared].copy()
    adata = adata[:, base_var_names.intersection(adata.var_names)].copy()
    # Reindex to base order and add zero-filled genes that are missing from this dataset.
    if adata.var_names.equals(base_var_names):
        return adata, missing_count

    shared_mask = base_var_names.isin(adata.var_names)
    shared_target = base_var_names[shared_mask]
    adata = adata[:, shared_target].copy()
    if shared_target.size < base_var_names.size:
        zero_block = sparse.csr_matrix((adata.n_obs, base_var_names.size - shared_target.size), dtype=np.float32)
        combined = sparse.hstack([adata.X.tocsr(), zero_block], format="csr")
        filler_var = pd.DataFrame(index=base_var_names[~shared_mask])
        adata = ad.AnnData(X=combined, obs=adata.obs.copy(), var=pd.concat([adata.var.copy(), filler_var], axis=0))
        adata.var_names = pd.Index(list(shared_target.astype(str)) + list(base_var_names[~shared_mask].astype(str)), dtype="string")
    adata = adata[:, base_var_names].copy()
    return adata, missing_count


def prepare_one_input(alias: str, path: Path, base_var_names: pd.Index) -> dict[str, object]:
    out_path = PLUS6_PREP_DIR / f"{alias}_plus6_prepped.h5ad"
    if out_path.exists():
        return summarize_existing_prepped(alias, out_path)
    logging.info("Preparing input %s from %s", alias, path)
    adata = ad.read_h5ad(path)
    choose_var_names(adata, set(base_var_names.astype(str)))
    count_source = move_counts_to_x(adata, alias)
    pseudo_source = infer_seurat_pseudocounts(adata, alias)
    if pseudo_source:
        count_source = pseudo_source
    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(np.asarray(adata.X), dtype=np.float32)
    else:
        adata.X = adata.X.tocsr().astype(np.float32)
    adata, prep_stats = standardize_generic_metadata(adata, alias)
    adata.var = pd.DataFrame(index=adata.var_names.astype(str))
    adata, missing_count = subset_to_base_genes(adata, base_var_names)
    make_obs_write_safe(adata)
    make_var_write_safe(adata)
    adata.write_h5ad(out_path, compression="gzip")
    prep_stats["count_matrix_source"] = count_source
    prep_stats["n_genes_after_align"] = int(adata.n_vars)
    prep_stats["n_missing_vs_base"] = missing_count
    prep_stats["output_h5ad"] = str(out_path)
    return prep_stats


def prepare_base(base_var_names: pd.Index) -> dict[str, object]:
    out_path = PLUS6_PREP_DIR / "current_integrated_base_plus6_prepped.h5ad"
    if out_path.exists():
        return summarize_existing_prepped("current_integrated_base", out_path)
    logging.info("Preparing current integrated base %s", BASE_INTEGRATED)
    adata = ad.read_h5ad(BASE_INTEGRATED)
    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(np.asarray(adata.X), dtype=np.float32)
    else:
        adata.X = adata.X.tocsr().astype(np.float32)
    adata.obs = drop_stale_columns(adata.obs.copy())
    adata, prep_stats = standardize_generic_metadata(adata, "current_integrated_base")
    adata.var = pd.DataFrame(index=adata.var_names.astype(str))
    adata = adata[:, base_var_names].copy()
    make_obs_write_safe(adata)
    make_var_write_safe(adata)
    adata.write_h5ad(out_path, compression="gzip")
    prep_stats["count_matrix_source"] = "X"
    prep_stats["n_genes_after_align"] = int(adata.n_vars)
    prep_stats["n_missing_vs_base"] = 0
    prep_stats["output_h5ad"] = str(out_path)
    return prep_stats


def summarize_existing_prepped(alias: str, path: Path) -> dict[str, object]:
    logging.info("Reusing existing prepped input for %s from %s", alias, path)
    with h5py.File(path, "r") as handle:
        n_cells = _read_h5ad_axis_len(handle, "obs")
        n_genes = _read_h5ad_axis_len(handle, "var")
        source_values = _obs_used_unique_values(handle, "source_gse_id")
        sample_values = _obs_used_unique_values(handle, "sample_id")
        source_label = "multiple_existing_gses" if source_values.size > 1 else str(source_values[0]) if source_values.size == 1 else alias
        n_samples = int(sample_values.size)
        paired_ab = _obs_true_count(handle, "has_TRA_TRB_paired")
        paired_gd = _obs_true_count(handle, "has_TRG_TRD_paired")
        sorted_flag = _obs_all_true(handle, "Sorted_gdT")

    stats = {
        "alias": alias,
        "source_gse_id": source_label,
        "n_cells": int(n_cells),
        "n_genes_before_align": int(_source_n_vars(alias)),
        "n_samples": n_samples,
        "n_paired_ab": paired_ab,
        "n_paired_gd": paired_gd,
        "sorted_gdt_flag": bool(sorted_flag),
        "count_matrix_source": PREP_COUNT_MATRIX_SOURCE_HINTS.get(alias, "cached_prepped"),
        "n_genes_after_align": int(n_genes),
        "n_missing_vs_base": int(_source_missing_vs_base(alias)),
        "output_h5ad": str(path),
    }
    return stats


def ensure_ssd_source_copy(alias: str, path: Path) -> Path:
    if alias not in SSD_LOCAL_SOURCE_ALIASES:
        return path
    target = PLUS6_SOURCE_COPY_DIR / path.name
    if target.exists() and target.stat().st_size == path.stat().st_size:
        return target
    logging.info("Copying %s to SSD source cache %s", path, target)
    subprocess.run(["rsync", "-a", str(path), str(target)], check=True)
    return target


def stage_prepare() -> None:
    logging.info("Stage prepare: harmonizing current base plus six new inputs")
    base_backed = ad.read_h5ad(BASE_INTEGRATED, backed="r")
    base_var_names = pd.Index(base_backed.var_names.astype(str), dtype="string")
    prep_rows = [prepare_base(base_var_names)]
    for alias, path in NEW_INPUTS.items():
        prep_rows.append(prepare_one_input(alias, ensure_ssd_source_copy(alias, path), base_var_names))

    prep_df = pd.DataFrame(prep_rows).sort_values(["alias"])
    prep_df.to_csv(PREP_SUMMARY_CSV, index=False)
    total_cells = int(prep_df["n_cells"].sum())
    lines = [
        "# plus6 input compatibility summary",
        "",
        f"- Base milestone: `{BASE_INTEGRATED}`",
        f"- Prepared inputs: `{len(prep_df)}`",
        f"- Total cells after merge prep: `{total_cells:,}`",
        "",
        "## Prepared datasets",
        "",
    ]
    for row in prep_df.itertuples(index=False):
        lines.extend(
            [
                f"### {row.alias}",
                f"- source_gse_id: `{row.source_gse_id}`",
                f"- cells: `{row.n_cells:,}`",
                f"- samples: `{row.n_samples:,}`",
                f"- count matrix source: `{row.count_matrix_source}`",
                f"- missing genes vs base: `{row.n_missing_vs_base:,}`",
                f"- paired abTCR: `{row.n_paired_ab:,}`",
                f"- paired gdTCR: `{row.n_paired_gd:,}`",
                f"- prepped H5AD: `{row.output_h5ad}`",
                "",
            ]
        )
    PREP_QC_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_prepped_inputs() -> list[ad.AnnData]:
    paths = [PLUS6_PREP_DIR / "current_integrated_base_plus6_prepped.h5ad"] + [
        PLUS6_PREP_DIR / f"{alias}_plus6_prepped.h5ad" for alias in NEW_INPUTS
    ]
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing prepared plus6 inputs: {missing}")
    return [ad.read_h5ad(path) for path in paths]


def build_batch_key(adata: ad.AnnData) -> None:
    def clean(column: str) -> pd.Series:
        if column not in adata.obs.columns:
            return pd.Series(pd.NA, index=adata.obs_names, dtype="string")
        return normalize_text_series(pd.Series(adata.obs[column], index=adata.obs_names))

    source_gse = clean("source_gse_id").fillna("unknown_gse")
    library_id = clean("library_id")
    sampleid = clean("sampleid")
    sample_id = clean("sample_id")

    batch_level = np.full(adata.n_obs, "gse", dtype=object)
    batch_key = source_gse.to_numpy(dtype=object, copy=True)
    library_mask = library_id.notna()
    batch_key[library_mask] = (source_gse[library_mask] + "|library|" + library_id[library_mask]).to_numpy(dtype=object)
    batch_level[library_mask] = "library_id"
    sampleid_mask = (~library_mask) & sampleid.notna()
    batch_key[sampleid_mask] = (source_gse[sampleid_mask] + "|sampleid|" + sampleid[sampleid_mask]).to_numpy(dtype=object)
    batch_level[sampleid_mask] = "sampleid"
    sample_mask = (~library_mask) & (~sampleid_mask) & sample_id.notna()
    batch_key[sample_mask] = (source_gse[sample_mask] + "|sample|" + sample_id[sample_mask]).to_numpy(dtype=object)
    batch_level[sample_mask] = "sample_id"
    adata.obs["phase3_batch_key"] = pd.Categorical(batch_key)
    adata.obs["phase3_batch_level"] = pd.Categorical(batch_level)


def select_hvgs(adata: ad.AnnData, n_hvgs: int) -> tuple[list[str], str, pd.DataFrame]:
    import scanpy as sc

    hvg_batch_key = "source_gse_id" if "source_gse_id" in adata.obs.columns else "phase3_batch_key"
    exclusion = build_hvg_exclusion_frame(adata.var_names)
    exclude_mask = exclusion["exclude_from_hvg"].to_numpy(dtype=bool, copy=False)
    allowed_mask = ~exclude_mask
    candidate_hvgs = min(adata.n_vars, max(n_hvgs * 3, n_hvgs + 5000))

    try:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=candidate_hvgs,
            flavor="seurat_v3",
            batch_key=hvg_batch_key,
            subset=False,
        )
        hvg_method = "seurat_v3"
    except Exception as exc:
        logging.warning("seurat_v3 HVG selection failed for plus6 (%s); falling back to cell_ranger", exc)
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=candidate_hvgs,
            flavor="cell_ranger",
            subset=False,
        )
        hvg_method = "cell_ranger_fallback"

    if "highly_variable_rank" in adata.var.columns:
        rank = pd.to_numeric(adata.var["highly_variable_rank"], errors="coerce")
        score = rank.copy()
        score[~allowed_mask] = np.inf
        selected = score.nsmallest(n_hvgs).index.astype(str).tolist()
    else:
        if "variances_norm" in adata.var.columns:
            score = pd.to_numeric(adata.var["variances_norm"], errors="coerce").fillna(-np.inf)
        else:
            raise ValueError("No HVG ranking metric found for plus6 run.")
        score[~allowed_mask] = -np.inf
        selected = score.nlargest(n_hvgs).index.astype(str).tolist()

    if len(selected) < n_hvgs:
        raise ValueError("Unable to recover enough allowed HVGs for plus6.")

    adata.var["phase3_hvg_excluded"] = exclude_mask
    adata.var["phase3_hvg_exclusion_reason"] = exclusion["exclude_reason"].to_numpy()
    adata.var["highly_variable"] = adata.var_names.isin(selected)

    summary = pd.DataFrame(
        [
            {
                "genes_total": int(adata.n_vars),
                "genes_allowed_for_hvg": int(allowed_mask.sum()),
                "genes_excluded_total": int(exclude_mask.sum()),
                "genes_excluded_mitochondrial": int(exclusion["exclude_mitochondrial"].sum()),
                "genes_excluded_ribosomal": int(exclusion["exclude_ribosomal"].sum()),
                "genes_excluded_noncoding_like": int(exclusion["exclude_noncoding_like"].sum()),
                "hvgs_selected": int(len(selected)),
            }
        ]
    )
    return selected, hvg_method, summary


def run_scvi_integration(adata: ad.AnnData):
    import scvi

    build_batch_key(adata)
    hvg_genes, hvg_method, hvg_summary = select_hvgs(adata, HVG_N)
    hvg_summary.to_csv(PHASE3_HVG_FILTER_CSV, index=False)
    PHASE3_HVG_GENES_TXT.write_text("\n".join(hvg_genes) + "\n", encoding="utf-8")
    adata_train = adata[:, hvg_genes].copy()
    scvi.model.SCVI.setup_anndata(adata_train, batch_key="phase3_batch_key")
    model = scvi.model.SCVI(
        adata_train,
        n_hidden=128,
        n_layers=2,
        n_latent=30,
        gene_likelihood="nb",
    )
    model.train(
        max_epochs=SCVI_MAX_EPOCHS,
        batch_size=SCVI_BATCH_SIZE,
        accelerator="gpu",
        devices=1,
        early_stopping=True,
        early_stopping_patience=5,
    )
    model.save(MODEL_DIR, overwrite=True, save_anndata=False)
    adata.obsm["X_scVI"] = model.get_latent_representation(adata_train).astype(np.float32, copy=False)
    return model, adata_train, hvg_method


def run_rapids_embedding(adata: ad.AnnData) -> None:
    import rapids_singlecell as rsc

    latent_obs = pd.DataFrame(index=adata.obs_names.copy())
    latent_adata = ad.AnnData(X=adata.obsm["X_scVI"].copy(), obs=latent_obs)
    rsc.get.anndata_to_GPU(latent_adata)
    rsc.pp.neighbors(latent_adata, n_neighbors=15, metric="euclidean")
    rsc.tl.leiden(latent_adata, resolution=LEIDEN_RESOLUTION)
    rsc.tl.umap(latent_adata)
    rsc.get.anndata_to_CPU(latent_adata)
    adata.obs["leiden"] = latent_adata.obs["leiden"].astype(str).to_numpy()
    adata.obsm["X_umap"] = np.asarray(latent_adata.obsm["X_umap"]).astype(np.float32, copy=False)


def write_phase3_figures(adata: ad.AnnData) -> None:
    sns.set_theme(style="white", context="talk")
    plot_df = sample_plot_df(adata)
    for column, out_name in [
        ("source_gse_id", "plus6_phase3_umap_by_gse.png"),
        ("tissue", "plus6_phase3_umap_by_tissue.png"),
        ("leiden", "plus6_phase3_umap_by_leiden.png"),
    ]:
        save_categorical_umap_with_legend(plot_df, column, FIGURE_DIR / out_name, f"plus6 Phase 3 UMAP by {column}")


def stage_integrate() -> None:
    logging.info("Stage integrate: building plus6 prephase3 and running fresh scVI")
    adatas = load_prepped_inputs()
    merged = ad.concat(adatas, axis=0, join="outer", merge="same", uns_merge="same", index_unique=None)
    merged.X = merged.X.tocsr().astype(np.float32)
    merged.write_h5ad(PLUS6_PREPHASE3, compression="gzip")

    start = time.time()
    model, adata_train, hvg_method = run_scvi_integration(merged)
    run_rapids_embedding(merged)
    make_obs_write_safe(merged)
    merged.write_h5ad(PLUS6_INTEGRATED, compression="gzip")
    elapsed = time.time() - start

    leiden_counts = merged.obs["leiden"].astype(str).value_counts().rename_axis("leiden").reset_index(name="n_cells")
    leiden_counts.to_csv(PHASE3_LEIDEN_COUNTS_CSV, index=False)
    summary = pd.DataFrame(
        [
            {
                "n_cells": int(merged.n_obs),
                "n_genes": int(merged.n_vars),
                "n_batches": int(merged.obs["phase3_batch_key"].nunique()),
                "n_gses": int(merged.obs["source_gse_id"].astype(str).nunique()),
                "n_leiden": int(merged.obs["leiden"].astype(str).nunique()),
                "hvg_method": hvg_method,
                "runtime_seconds": elapsed,
            }
        ]
    )
    summary.to_csv(PHASE3_SUMMARY_CSV, index=False)
    write_phase3_figures(merged)
    del model, adata_train, merged


def read_string_dataset(dataset: h5py.Dataset) -> np.ndarray:
    values = dataset[:]
    return np.asarray([value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values], dtype=object)


def _read_h5ad_axis_len(handle: h5py.File, axis: str) -> int:
    group = handle[axis]
    index_name = group.attrs.get("_index")
    if index_name is not None and index_name in group:
        return int(group[index_name].shape[0])
    return int(group.shape[0])


def _used_categorical_values(group: h5py.Group) -> np.ndarray:
    categories = read_string_dataset(group["categories"])
    codes = group["codes"][:]
    valid_codes = codes[codes >= 0]
    if valid_codes.size == 0:
        return np.asarray([], dtype=object)
    return categories[np.unique(valid_codes)]


def _obs_used_unique_values(handle: h5py.File, column: str) -> np.ndarray:
    if column not in handle["obs"]:
        return np.asarray([], dtype=object)
    obj = handle["obs"][column]
    if isinstance(obj, h5py.Group) and obj.attrs.get("encoding-type") == "categorical":
        values = _used_categorical_values(obj)
    else:
        if isinstance(obj, h5py.Dataset):
            raw = obj[:]
        else:
            return np.asarray([], dtype=object)
        if raw.dtype.kind in {"S", "O", "U"}:
            values = np.unique(np.asarray([v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in raw], dtype=object))
        else:
            values = np.unique(raw)
    values = np.asarray([str(v) for v in values if str(v) != ""], dtype=object)
    return values


def read_obs_values(handle: h5py.File, column: str) -> np.ndarray:
    if column not in handle["obs"]:
        raise KeyError(f"obs column not found: {column}")
    obj = handle["obs"][column]
    if isinstance(obj, h5py.Group) and obj.attrs.get("encoding-type") == "categorical":
        categories = read_string_dataset(obj["categories"])
        codes = obj["codes"][:]
        out = np.full(codes.shape, "", dtype=object)
        valid = codes >= 0
        out[valid] = categories[codes[valid]]
        return out
    if isinstance(obj, h5py.Dataset):
        raw = obj[:]
        if raw.dtype.kind in {"S", "O", "U"}:
            return np.asarray([v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in raw], dtype=object)
        return np.asarray(raw)
    raise TypeError(f"Unsupported obs column storage for {column}: {type(obj)}")


def _obs_true_count(handle: h5py.File, column: str) -> int:
    if column not in handle["obs"]:
        return 0
    obj = handle["obs"][column]
    if isinstance(obj, h5py.Dataset):
        return int(np.asarray(obj[:], dtype=bool).sum())
    if isinstance(obj, h5py.Group) and obj.attrs.get("encoding-type") == "categorical":
        values = _used_categorical_values(obj)
        return int(np.sum(np.char.lower(values.astype(str)) == "true"))
    return 0


def _obs_all_true(handle: h5py.File, column: str) -> bool:
    if column not in handle["obs"]:
        return False
    obj = handle["obs"][column]
    if isinstance(obj, h5py.Group) and obj.attrs.get("encoding-type") == "categorical":
        values = _used_categorical_values(obj)
        if values.size == 0:
            return False
        lowered = np.char.lower(values.astype(str))
        return bool(np.all(lowered == "true"))
    if isinstance(obj, h5py.Dataset):
        arr = np.asarray(obj[:])
        if arr.dtype == bool:
            return bool(arr.all())
        lowered = np.char.lower(arr.astype(str))
        return bool(np.all(lowered == "true"))
    return False


def _source_path_for_alias(alias: str) -> Path:
    if alias == "current_integrated_base":
        return BASE_INTEGRATED
    return NEW_INPUTS[alias]


def _source_n_vars(alias: str) -> int:
    path = _source_path_for_alias(alias)
    with h5py.File(path, "r") as handle:
        return _read_h5ad_axis_len(handle, "var")


def _source_missing_vs_base(alias: str) -> int:
    if alias == "current_integrated_base":
        return 0
    source_path = _source_path_for_alias(alias)
    with h5py.File(BASE_INTEGRATED, "r") as base_handle, h5py.File(source_path, "r") as src_handle:
        base_var = read_string_dataset(base_handle["var"][base_handle["var"].attrs["_index"]])
        src_var = read_string_dataset(src_handle["var"][src_handle["var"].attrs["_index"]])
    return int(len(set(base_var.tolist()) - set(src_var.tolist())))


def write_plus6_phase4_tables(adata: ad.AnnData, scores: dict[str, np.ndarray]) -> None:
    score_rows = []
    for module_name, column_name in {**PHASE4_SCORE_COLUMNS, **PHASE4_SCALED_SCORE_COLUMNS}.items():
        values = scores[module_name]
        score_rows.append(
            {
                "score_name": column_name,
                "n_cells": values.size,
                "mean": float(np.mean(values)),
                "std": float(np.std(values)),
                "min": float(np.min(values)),
                "median": float(np.quantile(values, 0.50)),
                "p95": float(np.quantile(values, 0.95)),
                "p99": float(np.quantile(values, 0.99)),
                "max": float(np.max(values)),
            }
        )
    pd.DataFrame(score_rows).to_csv(PHASE4_SCORE_SUMMARY_CSV, index=False)

    gse_df = pd.DataFrame(
        {
            "source_gse_id": adata.obs["source_gse_id"].astype(str).to_numpy(),
            "tissue": adata.obs.get("tissue_corrected", adata.obs.get("tissue", "")).astype(str).to_numpy(),
            "phase4_trab_score": scores["trab"],
            "phase4_trd_score": scores["trd"],
            "phase4_trd_minus_trab": scores["trd_minus_trab"],
        }
    )
    gse_summary = (
        gse_df.groupby("source_gse_id", observed=True)
        .agg(
            n_cells=("source_gse_id", "size"),
            phase4_trab_score_mean=("phase4_trab_score", "mean"),
            phase4_trd_score_mean=("phase4_trd_score", "mean"),
            phase4_trd_minus_trab_mean=("phase4_trd_minus_trab", "mean"),
            phase4_trd_minus_trab_median=("phase4_trd_minus_trab", "median"),
        )
        .reset_index()
        .sort_values("n_cells", ascending=False)
    )
    gse_summary.to_csv(PHASE4_GSE_SUMMARY_CSV, index=False)

    tissue_summary = (
        gse_df.groupby("tissue", observed=True)
        .agg(
            n_cells=("tissue", "size"),
            phase4_trab_score_mean=("phase4_trab_score", "mean"),
            phase4_trd_score_mean=("phase4_trd_score", "mean"),
            phase4_trd_minus_trab_mean=("phase4_trd_minus_trab", "mean"),
            phase4_trd_minus_trab_median=("phase4_trd_minus_trab", "median"),
        )
        .reset_index()
        .sort_values("n_cells", ascending=False)
    )
    tissue_summary.to_csv(PHASE4_TISSUE_SUMMARY_CSV, index=False)

    high_conf = gse_df.loc[gse_df["phase4_trd_minus_trab"] > 0.6].copy()
    broad = gse_df.loc[gse_df["phase4_trd_score"] > 0.1].copy()

    (
        high_conf.groupby("tissue", observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values("n_cells", ascending=False)
        .to_csv(PHASE4_TRD_OVER_TRAB_BY_TISSUE_CSV, index=False)
    )
    (
        high_conf.groupby("source_gse_id", observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values("n_cells", ascending=False)
        .to_csv(PHASE4_TRD_OVER_TRAB_BY_GSE_CSV, index=False)
    )
    (
        broad.groupby("tissue", observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values("n_cells", ascending=False)
        .to_csv(PHASE4_TRD_GT_0P1_BY_TISSUE_CSV, index=False)
    )
    (
        broad.groupby("source_gse_id", observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values("n_cells", ascending=False)
        .to_csv(PHASE4_TRD_GT_0P1_BY_GSE_CSV, index=False)
    )


def write_plus6_phase4_figures(adata: ad.AnnData, scores: dict[str, np.ndarray]) -> None:
    sns.set_theme(style="white", context="talk")
    plot_df = sample_plot_df(adata, sample_n=SCATTER_PLOT_SAMPLE_SIZE)

    for column, title, out_name, cmap in [
        ("phase4_trd_score", "plus6 UMAP by TRD score", "plus6_phase4_umap_trd_score.png", "viridis"),
        ("phase4_trab_score", "plus6 UMAP by TRAB score", "plus6_phase4_umap_trab_score.png", "viridis"),
        ("phase4_trd_minus_trab", "plus6 UMAP by TRD - TRAB", "plus6_phase4_umap_trd_minus_trab.png", "coolwarm"),
    ]:
        fig, ax = plt.subplots(figsize=(12, 10), constrained_layout=True)
        scatter = ax.scatter(
            plot_df["umap1"],
            plot_df["umap2"],
            c=plot_df[column],
            cmap=cmap,
            s=4,
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(title)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
        fig.savefig(FIGURE_DIR / out_name, dpi=FIGURE_DPI, bbox_inches="tight")
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for ax, column, title, color in [
        (axes[0], "has_TRA_TRB_paired", "Paired abTCR", "#C1121F"),
        (axes[1], "has_TRG_TRD_paired", "Paired gdTCR", "#C1121F"),
    ]:
        mask = plot_df[column].to_numpy(dtype=bool, copy=False)
        ax.scatter(plot_df.loc[~mask, "phase4_trab_score"], plot_df.loc[~mask, "phase4_trd_score"], s=3, c="#2F6690", linewidths=0, alpha=0.55, rasterized=True, label="Other")
        ax.scatter(plot_df.loc[mask, "phase4_trab_score"], plot_df.loc[mask, "phase4_trd_score"], s=3, c=color, linewidths=0, alpha=0.75, rasterized=True, label=title)
        ax.set_title(title)
        ax.set_xlabel("Raw TRAB score")
        ax.set_ylabel("Raw TRD score")
        ax.legend(loc="best", frameon=True)
    fig.savefig(FIGURE_DIR / "plus6_phase4_trab_vs_trd_paired_tcr.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
    scatter = ax.scatter(
        plot_df["phase4_trab_score"],
        plot_df["phase4_trd_score"],
        c=plot_df["phase4_trd_minus_trab"],
        cmap="coolwarm",
        s=4,
        linewidths=0,
        rasterized=True,
    )
    ax.set_title("plus6 Raw TRAB-versus-TRD score space")
    ax.set_xlabel("Raw TRAB score")
    ax.set_ylabel("Raw TRD score")
    fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04, label="TRD - TRAB")
    fig.savefig(FIGURE_DIR / "plus6_phase4_trab_vs_trd_raw.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)
    for ax, column, title, color in [
        (axes[0], "phase4_trab_score", "TRAB score", "#2F6690"),
        (axes[1], "phase4_trd_score", "TRD score", "#C1121F"),
        (axes[2], "phase4_trd_minus_trab", "TRD - TRAB", "#6C757D"),
    ]:
        sns.histplot(plot_df[column], bins=80, ax=ax, color=color, edgecolor=None)
        ax.set_title(title)
    fig.savefig(FIGURE_DIR / "plus6_phase4_score_distributions.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)

    high_conf_tissue = pd.read_csv(PHASE4_TRD_OVER_TRAB_BY_TISSUE_CSV)
    high_conf_gse = pd.read_csv(PHASE4_TRD_OVER_TRAB_BY_GSE_CSV)
    broad_tissue = pd.read_csv(PHASE4_TRD_GT_0P1_BY_TISSUE_CSV)
    broad_gse = pd.read_csv(PHASE4_TRD_GT_0P1_BY_GSE_CSV)
    horizontal_barplot(high_conf_tissue, "tissue", "n_cells", "High-Confidence TRD-over-TRAB Candidates by Tissue", FIGURE_DIR / "plus6_phase4_trd_over_trab_by_tissue_barplot.png")
    horizontal_barplot(high_conf_gse, "source_gse_id", "n_cells", "High-Confidence TRD-over-TRAB Candidates by GSE", FIGURE_DIR / "plus6_phase4_trd_over_trab_by_gse_barplot.png")
    horizontal_barplot(broad_tissue, "tissue", "n_cells", "Broad TRD-Enriched Candidates by Tissue", FIGURE_DIR / "plus6_phase4_trd_gt_0p1_by_tissue_barplot.png")
    horizontal_barplot(broad_gse, "source_gse_id", "n_cells", "Broad TRD-Enriched Candidates by GSE", FIGURE_DIR / "plus6_phase4_trd_gt_0p1_by_gse_barplot.png")


def build_phase4_uns_payload(modules: dict[str, pd.Index], controls: dict[str, pd.Index], scaling_stats: dict[str, dict[str, float]]) -> dict:
    return {
        "package_source": str(PACKAGE_ZIP),
        "integrated_h5ad": str(PLUS6_INTEGRATED),
        "scoring_mode": "temporary_normalize_total_log1p_on_count_space_X",
        "continuous_only": True,
        "random_state": RANDOM_STATE,
        "target_sum": TARGET_SUM,
        "ctrl_size": CTRL_SIZE,
        "n_bins": N_BINS,
        "module_genes": {name: [str(gene) for gene in genes] for name, genes in modules.items()},
        "control_genes": {name: [str(gene) for gene in genes] for name, genes in controls.items()},
        "score_columns": {**PHASE4_SCORE_COLUMNS, **PHASE4_SCALED_SCORE_COLUMNS},
        "scaled_score_ranges": scaling_stats,
        "plus6_lineage": True,
    }


def stage_score() -> None:
    logging.info("Stage score: rescoring integrated_plus6 with Phase 4")
    adata = ad.read_h5ad(PLUS6_INTEGRATED)
    var_names = pd.Index(adata.var_names.astype(str), dtype="string")
    n_obs, n_vars = adata.n_obs, adata.n_vars
    module_gene_sets = find_module_genes(var_names)
    gene_means = compute_gene_means(PLUS6_INTEGRATED, n_obs, n_vars, CHUNK_SIZE)
    module_controls = {
        name: pick_control_genes(genes, var_names, gene_means, random_state=RANDOM_STATE)
        for name, genes in module_gene_sets.items()
    }
    module_gene_idx = {name: var_names.get_indexer(genes).astype(np.int32, copy=False) for name, genes in module_gene_sets.items()}
    module_ctrl_idx = {name: var_names.get_indexer(genes).astype(np.int32, copy=False) for name, genes in module_controls.items()}
    leiden_labels = adata.obs["leiden"].astype(str).to_numpy()
    leiden_codes, _ = pd.factorize(leiden_labels, sort=True)

    marker_genes = [gene for gene in ANNOTATION_GENES if gene in var_names]
    marker_idx = var_names.get_indexer(pd.Index(marker_genes, dtype="string")).astype(np.int32, copy=False)
    scores, _, _ = compute_scores(
        PLUS6_INTEGRATED,
        n_obs=n_obs,
        n_vars=n_vars,
        chunk_size=CHUNK_SIZE,
        module_gene_idx=module_gene_idx,
        module_ctrl_idx=module_ctrl_idx,
        leiden_codes=leiden_codes,
        marker_idx=marker_idx,
    )
    scores, scaling_stats = add_scaled_scores(scores)
    score_columns = {column_name: scores[module_name] for module_name, column_name in {**PHASE4_SCORE_COLUMNS, **PHASE4_SCALED_SCORE_COLUMNS}.items()}
    append_obs_columns_in_place(PLUS6_INTEGRATED, score_columns, build_phase4_uns_payload(module_gene_sets, module_controls, scaling_stats))
    adata = ad.read_h5ad(PLUS6_INTEGRATED)
    write_plus6_phase4_tables(adata, scores)
    write_plus6_phase4_figures(adata, scores)


def extract_log1p_gene_expression_all(h5ad_path: Path, gene_names: list[str], chunk_size: int = CHUNK_SIZE) -> tuple[np.ndarray, list[str]]:
    with h5py.File(h5ad_path, "r") as handle:
        var_names = pd.Index(read_string_dataset(handle["var"]["_index"]), dtype="string")
        gene_idx = var_names.get_indexer(pd.Index(gene_names, dtype="string"))
        present_mask = gene_idx >= 0
        present_idx = gene_idx[present_mask].astype(np.int32, copy=False)
        n_obs = int(handle["obs"]["_index"].shape[0])
        expr = np.zeros((n_obs, len(gene_names)), dtype=np.float32)
        x_group = handle["X"]
        for start in range(0, n_obs, chunk_size):
            end = min(start + chunk_size, n_obs)
            indptr = x_group["indptr"][start : end + 1].astype(np.int64, copy=False)
            offset = int(indptr[0])
            indptr = indptr - offset
            nnz = int(indptr[-1])
            data = x_group["data"][offset : offset + nnz].astype(np.float32, copy=False)
            indices = x_group["indices"][offset : offset + nnz].astype(np.int32, copy=False)
            chunk = sparse.csr_matrix((data, indices, indptr), shape=(end - start, var_names.size))
            chunk = chunk.tocsr()
            row_sums = np.asarray(chunk.sum(axis=1)).ravel().astype(np.float32, copy=False)
            scale = np.zeros_like(row_sums, dtype=np.float32)
            mask = row_sums > 0
            scale[mask] = TARGET_SUM / row_sums[mask]
            sparsefuncs.inplace_row_scale(chunk, scale)
            np.log1p(chunk.data, out=chunk.data)
            if present_idx.size:
                values = chunk[:, present_idx].toarray().astype(np.float32, copy=False)
                expr[start:end, np.flatnonzero(present_mask)] = values
    return expr, list(gene_names)


def append_text_obs_columns(h5ad_path: Path, columns: dict[str, np.ndarray]) -> None:
    with h5py.File(h5ad_path, "r+") as handle:
        obs_group = handle["obs"]
        column_order = list(obs_group.attrs["column-order"])
        for final_name, values in columns.items():
            tmp_name = f"__{final_name}_tmp"
            for existing in (tmp_name, final_name):
                if existing in obs_group:
                    del obs_group[existing]
            write_elem(obs_group, tmp_name, np.asarray(values, dtype=object))
            obs_group.move(tmp_name, final_name)
            if final_name not in column_order:
                column_order.append(final_name)
        obs_group.attrs["column-order"] = np.asarray(column_order, dtype=object)


def classify_clusters(adata: ad.AnnData, marker_matrix: np.ndarray, marker_genes: list[str]) -> tuple[pd.DataFrame, np.ndarray]:
    leiden = adata.obs["leiden"].astype(str).to_numpy()
    source_gse_id = adata.obs["source_gse_id"].astype(str).to_numpy()
    phase4_trd_score = adata.obs["phase4_trd_score"].to_numpy(dtype=np.float32)
    phase4_trab_score = adata.obs["phase4_trab_score"].to_numpy(dtype=np.float32)
    phase4_trd_minus_trab = adata.obs["phase4_trd_minus_trab"].to_numpy(dtype=np.float32)
    has_TRA_TRB_paired = adata.obs["has_TRA_TRB_paired"].to_numpy(dtype=bool)
    has_TRG_TRD_paired = adata.obs["has_TRG_TRD_paired"].to_numpy(dtype=bool)
    gene_avail = {gene for gene in marker_genes if gene in ANNOTATION_GENES}
    gene_to_idx = {gene: idx for idx, gene in enumerate(marker_genes)}

    def mean_available(row_values: np.ndarray, genes: list[str]) -> float:
        avail = [gene_to_idx[gene] for gene in genes if gene in gene_avail]
        if not avail:
            return 0.0
        return float(np.mean(row_values[avail]))

    summaries = []
    cluster_labels = {}
    for cluster in sorted(pd.unique(leiden)):
        mask = leiden == cluster
        marker_med = np.median(marker_matrix[mask, :], axis=0)
        row = {
            "leiden": cluster,
            "n_cells": int(mask.sum()),
            "n_gses": int(pd.unique(source_gse_id[mask]).size),
            "paired_ab_fraction": float(has_TRA_TRB_paired[mask].mean()),
            "paired_gd_fraction": float(has_TRG_TRD_paired[mask].mean()),
            "phase4_trd_score_median": float(np.median(phase4_trd_score[mask])),
            "phase4_trab_score_median": float(np.median(phase4_trab_score[mask])),
            "phase4_trd_minus_trab_median": float(np.median(phase4_trd_minus_trab[mask])),
            "gdt_marker_mean": mean_available(marker_med, GD_ANNOTATION_ALLOWED_GENES),
            "treg_marker_mean": mean_available(marker_med, ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2"]),
            "nk_marker_mean": mean_available(marker_med, ["NKG7", "GNLY", "KLRD1", "FCER1G", "PRF1", "CTSW"]),
            "tcell_marker_mean": mean_available(marker_med, ["CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2"]),
            "cd4_marker_mean": mean_available(marker_med, ["CD4", "IL7R", "LTB", "MAL"]),
            "cd8_marker_mean": mean_available(marker_med, ["CD8A", "CD8B", "CCL5", "PRF1"]),
        }
        gdt_rule = (
            row["paired_gd_fraction"] >= 0.10
            or (row["phase4_trd_score_median"] >= 0.12 and row["phase4_trd_minus_trab_median"] >= 0.08 and row["gdt_marker_mean"] >= 0.15)
        )
        nk_rule = row["nk_marker_mean"] >= 0.25 and row["tcell_marker_mean"] < 0.10
        treg_rule = (
            row["tcell_marker_mean"] >= 0.10
            and row["treg_marker_mean"] >= 0.25
            and row["cd4_marker_mean"] >= row["cd8_marker_mean"]
            and row["nk_marker_mean"] < 0.35
        )
        cd8_rule = row["tcell_marker_mean"] >= 0.10 and row["cd8_marker_mean"] >= row["cd4_marker_mean"] + 0.03
        cd4_rule = row["tcell_marker_mean"] >= 0.08 and row["cd4_marker_mean"] >= row["cd8_marker_mean"]
        row["gdt_rule"] = bool(gdt_rule)
        row["nk_rule"] = bool(nk_rule)
        row["treg_rule"] = bool(treg_rule)
        row["cd8_rule"] = bool(cd8_rule)
        row["cd4_rule"] = bool(cd4_rule)
        row["treg_cd4_ge_cd8"] = bool(row["cd4_marker_mean"] >= row["cd8_marker_mean"])
        row["treg_nk_low"] = bool(row["nk_marker_mean"] < 0.35)
        if gdt_rule:
            label = "gdT_cell"
        elif nk_rule:
            label = "NK_cell"
        elif treg_rule:
            label = "Treg"
        elif cd8_rule:
            label = "CD8_T"
        elif cd4_rule:
            label = "CD4_T"
        else:
            label = "other"
        row["simple_annotation_plus6"] = label
        summaries.append(row)
        cluster_labels[cluster] = label

    summary_df = pd.DataFrame(summaries).sort_values(
        ["simple_annotation_plus6", "phase4_trd_minus_trab_median", "n_cells"],
        ascending=[True, False, False],
    )
    labels = pd.Series(leiden).map(cluster_labels).fillna("other").to_numpy(dtype=object)
    return summary_df, labels


def choose_plot_indices(n_obs: int, sample_n: int = LEGEND_SAMPLE_LIMIT, seed: int = RANDOM_STATE) -> np.ndarray:
    rng = np.random.default_rng(seed)
    sample_n = min(sample_n, n_obs)
    return np.sort(rng.choice(np.arange(n_obs), size=sample_n, replace=False))


def get_old_cluster_labels(adata: ad.AnnData) -> pd.DataFrame:
    if "simple_annotation_plus6" not in adata.obs.columns:
        return pd.DataFrame(columns=["leiden", "old_label"])
    old = (
        adata.obs[["leiden", "simple_annotation_plus6"]]
        .copy()
        .assign(leiden=lambda x: x["leiden"].astype(str), old_label=lambda x: x["simple_annotation_plus6"].astype(str))
        .drop(columns=["simple_annotation_plus6"])
    )
    old = old.groupby("leiden", observed=True)["old_label"].agg(lambda s: sorted(pd.unique(s.astype(str)))[0]).reset_index()
    return old


def build_annotation_change_table(cluster_summary: pd.DataFrame, old_cluster_labels: pd.DataFrame) -> pd.DataFrame:
    cluster_summary = cluster_summary.copy()
    old_cluster_labels = old_cluster_labels.copy()
    cluster_summary["leiden"] = cluster_summary["leiden"].astype(str)
    old_cluster_labels["leiden"] = old_cluster_labels["leiden"].astype(str)
    change_df = cluster_summary.merge(old_cluster_labels, on="leiden", how="left")
    change_df["old_label"] = change_df["old_label"].fillna("NA")
    change_df = change_df.rename(columns={"simple_annotation_plus6": "new_label"})
    change_df["changed"] = change_df["old_label"] != change_df["new_label"]
    ordered_cols = [
        "leiden",
        "n_cells",
        "old_label",
        "new_label",
        "changed",
        "treg_marker_mean",
        "cd4_marker_mean",
        "cd8_marker_mean",
        "nk_marker_mean",
        "tcell_marker_mean",
        "paired_gd_fraction",
        "phase4_trd_minus_trab_median",
        "gdt_rule",
        "nk_rule",
        "treg_rule",
        "cd8_rule",
        "cd4_rule",
    ]
    existing = [col for col in ordered_cols if col in change_df.columns]
    remaining = [col for col in change_df.columns if col not in existing]
    return change_df[existing + remaining].sort_values(["changed", "n_cells"], ascending=[False, False])


def relabel_cluster_summary(cluster_summary: pd.DataFrame) -> pd.DataFrame:
    df = cluster_summary.copy()
    df["gdt_rule"] = (
        (df["paired_gd_fraction"] >= 0.10)
        | (
            (df["phase4_trd_score_median"] >= 0.12)
            & (df["phase4_trd_minus_trab_median"] >= 0.08)
            & (df["gdt_marker_mean"] >= 0.15)
        )
    )
    df["nk_rule"] = (df["nk_marker_mean"] >= 0.25) & (df["tcell_marker_mean"] < 0.10)
    df["treg_rule"] = (
        (df["tcell_marker_mean"] >= 0.10)
        & (df["treg_marker_mean"] >= 0.25)
        & (df["cd4_marker_mean"] >= df["cd8_marker_mean"])
        & (df["nk_marker_mean"] < 0.35)
    )
    df["cd8_rule"] = (df["tcell_marker_mean"] >= 0.10) & (df["cd8_marker_mean"] >= df["cd4_marker_mean"] + 0.03)
    df["cd4_rule"] = (df["tcell_marker_mean"] >= 0.08) & (df["cd4_marker_mean"] >= df["cd8_marker_mean"])
    df["treg_cd4_ge_cd8"] = df["cd4_marker_mean"] >= df["cd8_marker_mean"]
    df["treg_nk_low"] = df["nk_marker_mean"] < 0.35

    labels = np.full(df.shape[0], "other", dtype=object)
    labels[df["gdt_rule"].to_numpy(dtype=bool)] = "gdT_cell"
    labels[(labels == "other") & df["nk_rule"].to_numpy(dtype=bool)] = "NK_cell"
    labels[(labels == "other") & df["treg_rule"].to_numpy(dtype=bool)] = "Treg"
    labels[(labels == "other") & df["cd8_rule"].to_numpy(dtype=bool)] = "CD8_T"
    labels[(labels == "other") & df["cd4_rule"].to_numpy(dtype=bool)] = "CD4_T"
    df["simple_annotation_plus6"] = labels
    return df.sort_values(
        ["simple_annotation_plus6", "phase4_trd_minus_trab_median", "n_cells"],
        ascending=[True, False, False],
    )


def legacy_relabel_cluster_summary(cluster_summary: pd.DataFrame) -> pd.DataFrame:
    df = cluster_summary.copy()
    labels = np.full(df.shape[0], "other", dtype=object)
    gdt_rule = (
        (df["paired_gd_fraction"] >= 0.10)
        | (
            (df["phase4_trd_score_median"] >= 0.12)
            & (df["phase4_trd_minus_trab_median"] >= 0.08)
            & (df["gdt_marker_mean"] >= 0.15)
        )
    )
    treg_rule = (df["treg_marker_mean"] >= 0.18) & (df["tcell_marker_mean"] >= 0.10)
    nk_rule = (df["nk_marker_mean"] >= 0.25) & (df["tcell_marker_mean"] < 0.10)
    cd8_rule = (df["tcell_marker_mean"] >= 0.10) & (df["cd8_marker_mean"] >= df["cd4_marker_mean"] + 0.03)
    cd4_rule = (df["tcell_marker_mean"] >= 0.08) & (df["cd4_marker_mean"] >= df["cd8_marker_mean"])
    labels[gdt_rule.to_numpy(dtype=bool)] = "gdT_cell"
    labels[(labels == "other") & treg_rule.to_numpy(dtype=bool)] = "Treg"
    labels[(labels == "other") & nk_rule.to_numpy(dtype=bool)] = "NK_cell"
    labels[(labels == "other") & cd8_rule.to_numpy(dtype=bool)] = "CD8_T"
    labels[(labels == "other") & cd4_rule.to_numpy(dtype=bool)] = "CD4_T"
    return pd.DataFrame({"leiden": df["leiden"].astype(str), "old_label": labels})


def make_feature_plot_panel(plot_df: pd.DataFrame, genes: list[str], out_path: Path, title: str, n_cols: int = 2) -> None:
    n_rows = int(math.ceil(len(genes) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 5 * n_rows), constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, gene in zip(axes, genes, strict=False):
        values = plot_df[gene].to_numpy(dtype=np.float32)
        order = np.argsort(values)
        sorted_xy = plot_df.loc[order, ["umap1", "umap2"]].to_numpy(dtype=np.float32)
        sorted_values = values[order]
        vmax = float(np.percentile(sorted_values, 99.0)) if np.any(sorted_values > 0) else 1.0
        if vmax <= 0:
            vmax = 1.0
        scatter = ax.scatter(
            sorted_xy[:, 0],
            sorted_xy[:, 1],
            c=sorted_values,
            s=4,
            cmap="viridis",
            vmin=0.0,
            vmax=vmax,
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(gene)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.02)
        cbar.ax.tick_params(labelsize=8)
    for ax in axes[len(genes):]:
        ax.axis("off")
    fig.suptitle(title, fontsize=18)
    fig.savefig(out_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def write_treg_cd8_diagnostic_figures(
    adata: ad.AnnData,
    marker_matrix: np.ndarray,
    marker_genes: list[str],
    cluster_summary: pd.DataFrame,
    change_df: pd.DataFrame,
) -> None:
    sns.set_theme(style="white", context="talk")
    idx = choose_plot_indices(adata.n_obs)
    plot_df = pd.DataFrame(
        {
            "umap1": adata.obsm["X_umap"][idx, 0],
            "umap2": adata.obsm["X_umap"][idx, 1],
            "simple_annotation_plus6": adata.obs["simple_annotation_plus6"].astype(str).to_numpy()[idx],
        }
    )
    gene_to_idx = {gene: i for i, gene in enumerate(marker_genes)}
    sampled_expr = pd.DataFrame(
        marker_matrix[np.asarray(idx, dtype=np.int64)][:, [gene_to_idx[gene] for gene in TREG_DIAGNOSTIC_GENES]],
        columns=TREG_DIAGNOSTIC_GENES,
    )
    plot_df = pd.concat([plot_df.reset_index(drop=True), sampled_expr.reset_index(drop=True)], axis=1)
    make_feature_plot_panel(
        plot_df,
        TREG_DIAGNOSTIC_GENES,
        FIGURE_DIR / "plus6_umap_treg_cd8_diagnostic.png",
        "plus6 Treg/CD8 diagnostic feature plots",
    )
    save_categorical_umap_with_legend(
        plot_df,
        "simple_annotation_plus6",
        FIGURE_DIR / "plus6_umap_by_simple_annotation_corrected.png",
        "plus6 UMAP by corrected simple annotation",
    )

    selected_clusters = sorted(
        set(change_df.loc[change_df["changed"], "leiden"].astype(str))
        | set(cluster_summary.loc[cluster_summary["simple_annotation_plus6"] == "Treg", "leiden"].astype(str))
        | set(change_df.loc[change_df["new_label"] == "CD8_T", "leiden"].astype(str))
    )
    if not selected_clusters:
        return
    leiden = adata.obs["leiden"].astype(str).to_numpy()
    dot_gene_idx = [gene_to_idx[gene] for gene in TREG_CLUSTER_DOTPLOT_GENES]
    dot_matrix = marker_matrix[:, dot_gene_idx]
    mean_expr_rows = []
    frac_expr_rows = []
    for cluster in selected_clusters:
        mask = leiden == cluster
        cluster_matrix = dot_matrix[mask, :]
        mean_expr_rows.append(cluster_matrix.mean(axis=0))
        frac_expr_rows.append((cluster_matrix > 0).mean(axis=0))
    mean_expr = pd.DataFrame(mean_expr_rows, index=selected_clusters, columns=TREG_CLUSTER_DOTPLOT_GENES)
    frac_expr = pd.DataFrame(frac_expr_rows, index=selected_clusters, columns=TREG_CLUSTER_DOTPLOT_GENES)

    x_positions = np.arange(len(TREG_CLUSTER_DOTPLOT_GENES))
    y_positions = np.arange(len(selected_clusters))
    fig, ax = plt.subplots(figsize=(14, max(4, 0.55 * len(selected_clusters) + 2)), constrained_layout=True)
    norm = plt.Normalize(vmin=0.0, vmax=float(np.percentile(mean_expr.to_numpy().ravel(), 95.0)) if np.any(mean_expr.to_numpy() > 0) else 1.0)
    for yi, cluster in enumerate(selected_clusters):
        for xi, gene in enumerate(TREG_CLUSTER_DOTPLOT_GENES):
            ax.scatter(
                xi,
                yi,
                s=30 + 260 * float(frac_expr.loc[cluster, gene]),
                c=[plt.cm.viridis(norm(float(mean_expr.loc[cluster, gene])))],
                edgecolors="none",
            )
    ax.set_xticks(x_positions)
    ax.set_xticklabels(TREG_CLUSTER_DOTPLOT_GENES, rotation=45, ha="right")
    ax.set_yticks(y_positions)
    ax.set_yticklabels(selected_clusters)
    ax.set_xlabel("Marker gene")
    ax.set_ylabel("Leiden cluster")
    ax.set_title("plus6 Treg/CD8 cluster dotplot")
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Mean log1p expression")
    fig.savefig(FIGURE_DIR / "plus6_treg_cd8_cluster_dotplot.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def classify_clusters_legacy(adata: ad.AnnData, marker_matrix: np.ndarray, marker_genes: list[str]) -> tuple[pd.DataFrame, np.ndarray]:
    """Preserve the previous annotation logic for change auditing."""
    score_df = pd.DataFrame(
        {
            "leiden": adata.obs["leiden"].astype(str).to_numpy(),
            "source_gse_id": adata.obs["source_gse_id"].astype(str).to_numpy(),
            "phase4_trd_score": adata.obs["phase4_trd_score"].to_numpy(dtype=np.float32),
            "phase4_trab_score": adata.obs["phase4_trab_score"].to_numpy(dtype=np.float32),
            "phase4_trd_minus_trab": adata.obs["phase4_trd_minus_trab"].to_numpy(dtype=np.float32),
            "has_TRA_TRB_paired": adata.obs["has_TRA_TRB_paired"].to_numpy(dtype=bool),
            "has_TRG_TRD_paired": adata.obs["has_TRG_TRD_paired"].to_numpy(dtype=bool),
        }
    )
    gene_avail = {gene for gene in marker_genes if gene in ANNOTATION_GENES}
    gene_to_idx = {gene: idx for idx, gene in enumerate(marker_genes)}

    def mean_available(row_values: np.ndarray, genes: list[str]) -> float:
        avail = [gene_to_idx[gene] for gene in genes if gene in gene_avail]
        if not avail:
            return 0.0
        return float(np.mean(row_values[avail]))

    summaries = []
    cluster_labels = {}
    leiden = score_df["leiden"].to_numpy()
    source_gse = score_df["source_gse_id"].to_numpy()
    phase4_trd_score = score_df["phase4_trd_score"].to_numpy(dtype=np.float32)
    phase4_trab_score = score_df["phase4_trab_score"].to_numpy(dtype=np.float32)
    phase4_trd_minus_trab = score_df["phase4_trd_minus_trab"].to_numpy(dtype=np.float32)
    has_ab = score_df["has_TRA_TRB_paired"].to_numpy(dtype=bool)
    has_gd = score_df["has_TRG_TRD_paired"].to_numpy(dtype=bool)
    for cluster in sorted(pd.unique(leiden)):
        mask = leiden == cluster
        marker_med = np.median(marker_matrix[mask, :], axis=0)
        row = {
            "leiden": cluster,
            "n_cells": int(mask.sum()),
            "n_gses": int(pd.unique(source_gse[mask]).size),
            "paired_ab_fraction": float(has_ab[mask].mean()),
            "paired_gd_fraction": float(has_gd[mask].mean()),
            "phase4_trd_score_median": float(np.median(phase4_trd_score[mask])),
            "phase4_trab_score_median": float(np.median(phase4_trab_score[mask])),
            "phase4_trd_minus_trab_median": float(np.median(phase4_trd_minus_trab[mask])),
            "gdt_marker_mean": mean_available(marker_med, GD_ANNOTATION_ALLOWED_GENES),
            "treg_marker_mean": mean_available(marker_med, ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2"]),
            "nk_marker_mean": mean_available(marker_med, ["NKG7", "GNLY", "KLRD1", "FCER1G", "PRF1", "CTSW"]),
            "tcell_marker_mean": mean_available(marker_med, ["CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2"]),
            "cd4_marker_mean": mean_available(marker_med, ["CD4", "IL7R", "LTB", "MAL"]),
            "cd8_marker_mean": mean_available(marker_med, ["CD8A", "CD8B", "CCL5", "PRF1"]),
        }
        if (
            row["paired_gd_fraction"] >= 0.10
            or (row["phase4_trd_score_median"] >= 0.12 and row["phase4_trd_minus_trab_median"] >= 0.08 and row["gdt_marker_mean"] >= 0.15)
        ):
            label = "gdT_cell"
        elif row["treg_marker_mean"] >= 0.18 and row["tcell_marker_mean"] >= 0.10:
            label = "Treg"
        elif row["nk_marker_mean"] >= 0.25 and row["tcell_marker_mean"] < 0.10:
            label = "NK_cell"
        elif row["tcell_marker_mean"] >= 0.10 and row["cd8_marker_mean"] >= row["cd4_marker_mean"] + 0.03:
            label = "CD8_T"
        elif row["tcell_marker_mean"] >= 0.08 and row["cd4_marker_mean"] >= row["cd8_marker_mean"]:
            label = "CD4_T"
        else:
            label = "other"
        row["simple_annotation_plus6"] = label
        summaries.append(row)
        cluster_labels[cluster] = label

    summary_df = pd.DataFrame(summaries).sort_values(
        ["simple_annotation_plus6", "phase4_trd_minus_trab_median", "n_cells"],
        ascending=[True, False, False],
    )
    labels = pd.Series(leiden).map(cluster_labels).fillna("other").to_numpy(dtype=object)
    return summary_df, labels


def write_annotation_figures(adata: ad.AnnData) -> None:
    sns.set_theme(style="white", context="talk")
    plot_df = sample_plot_df(adata)
    for column, out_name in [
        ("simple_annotation_plus6", "plus6_umap_by_simple_annotation.png"),
        ("simple_annotation_plus6", "plus6_umap_by_simple_annotation_corrected.png"),
        ("source_gse_id", "plus6_umap_by_gse_after_annotation.png"),
        ("tissue", "plus6_umap_by_tissue_after_annotation.png"),
    ]:
        save_categorical_umap_with_legend(plot_df, column, FIGURE_DIR / out_name, f"plus6 UMAP by {column}")


def stage_annotate() -> None:
    logging.info("Stage annotate: adding simple_annotation_plus6")
    old_cluster_summary = pd.read_csv(ANNOTATION_CLUSTER_SUMMARY_CSV)
    old_cluster_labels = legacy_relabel_cluster_summary(old_cluster_summary)
    cluster_summary = relabel_cluster_summary(old_cluster_summary.drop(columns=["simple_annotation_plus6"], errors="ignore"))
    with h5py.File(PLUS6_INTEGRATED, "r") as handle:
        leiden = read_obs_values(handle, "leiden").astype(str)
        source_gse_id = read_obs_values(handle, "source_gse_id").astype(str)
        tissue_col = "tissue_corrected" if "tissue_corrected" in handle["obs"] else "tissue"
        tissue = read_obs_values(handle, tissue_col).astype(str)
        n_obs = int(handle["obs"]["_index"].shape[0])
        n_vars = int(handle["var"]["_index"].shape[0])
        umap = np.asarray(handle["obsm"]["X_umap"], dtype=np.float32)
    label_map = dict(zip(cluster_summary["leiden"].astype(str), cluster_summary["simple_annotation_plus6"].astype(str), strict=False))
    labels = pd.Series(leiden).map(label_map).fillna("other").to_numpy(dtype=object)
    append_text_obs_columns(
        PLUS6_INTEGRATED,
        {
            "simple_annotation_plus6": labels,
        },
    )
    change_df = build_annotation_change_table(cluster_summary, old_cluster_labels)
    cluster_summary.to_csv(ANNOTATION_CLUSTER_SUMMARY_CSV, index=False)
    change_df.to_csv(ANNOTATION_CHANGES_CSV, index=False)
    (
        pd.DataFrame({"source_gse_id": source_gse_id, "simple_annotation_plus6": labels})
        .groupby(["source_gse_id", "simple_annotation_plus6"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values(["source_gse_id", "n_cells"], ascending=[True, False])
        .to_csv(ANNOTATION_GSE_COUNTS_CSV, index=False)
    )
    (
        pd.DataFrame({tissue_col: tissue, "simple_annotation_plus6": labels})
        .groupby([tissue_col, "simple_annotation_plus6"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
        .sort_values([tissue_col, "n_cells"], ascending=[True, False])
        .to_csv(ANNOTATION_TISSUE_COUNTS_CSV, index=False)
    )
    changed_clusters = change_df.loc[change_df["changed"]].copy()
    treg_clusters = cluster_summary.loc[cluster_summary["simple_annotation_plus6"] == "Treg"].copy()
    treg_cd8_fail = treg_clusters.loc[treg_clusters["cd8_marker_mean"] > treg_clusters["cd4_marker_mean"]]
    treg_nk_fail = treg_clusters.loc[treg_clusters["nk_marker_mean"] >= 0.35]
    lines = [
        "# plus6 simple annotation summary",
        "",
        f"- gdT markers used in decision logic: `{', '.join(GD_ANNOTATION_ALLOWED_GENES)}`",
        "- gdT markers excluded from decision logic: `TRGC1`, `TRGC2`",
        f"- clusters labelled: `{cluster_summary.shape[0]}`",
        f"- changed clusters after Treg/CD8 correction: `{int(changed_clusters.shape[0])}`",
        f"- Treg clusters violating `cd4_marker_mean >= cd8_marker_mean`: `{int(treg_cd8_fail.shape[0])}`",
        f"- Treg clusters violating `nk_marker_mean < 0.35`: `{int(treg_nk_fail.shape[0])}`",
        "",
        "## Changed clusters",
        "",
    ]
    if changed_clusters.empty:
        lines.append("- none")
    else:
        for row in changed_clusters.itertuples(index=False):
            lines.append(
                f"- Leiden {row.leiden}: `{row.old_label}` -> `{row.new_label}` | n={row.n_cells:,} | treg={row.treg_marker_mean:.3f} | cd4={row.cd4_marker_mean:.3f} | cd8={row.cd8_marker_mean:.3f} | nk={row.nk_marker_mean:.3f}"
            )
    lines.extend(
        [
            "",
        "## Cluster labels",
        "",
        ]
    )
    for row in cluster_summary.itertuples(index=False):
        lines.append(
            f"- Leiden {row.leiden}: `{row.simple_annotation_plus6}` | n={row.n_cells:,} | paired_gd={row.paired_gd_fraction:.3f} | paired_ab={row.paired_ab_fraction:.3f} | TRD-TRAB median={row.phase4_trd_minus_trab_median:.3f} | treg={row.treg_marker_mean:.3f} | cd4={row.cd4_marker_mean:.3f} | cd8={row.cd8_marker_mean:.3f} | nk={row.nk_marker_mean:.3f}"
        )
    ANNOTATION_SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    sample_idx = choose_plot_indices(n_obs)
    plot_df = pd.DataFrame(
        {
            "umap1": umap[sample_idx, 0],
            "umap2": umap[sample_idx, 1],
            "simple_annotation_plus6": labels[sample_idx],
            "source_gse_id": source_gse_id[sample_idx],
            "tissue": tissue[sample_idx],
        }
    )
    for column, out_name in [
        ("simple_annotation_plus6", "plus6_umap_by_simple_annotation.png"),
        ("simple_annotation_plus6", "plus6_umap_by_simple_annotation_corrected.png"),
        ("source_gse_id", "plus6_umap_by_gse_after_annotation.png"),
        ("tissue", "plus6_umap_by_tissue_after_annotation.png"),
    ]:
        save_categorical_umap_with_legend(plot_df, column, FIGURE_DIR / out_name, f"plus6 UMAP by {column}")

    sample_expr = extract_log1p_gene_expression_for_sample(PLUS6_INTEGRATED, sample_idx, TREG_DIAGNOSTIC_GENES, CHUNK_SIZE)
    diagnostic_df = pd.concat([plot_df[["umap1", "umap2", "simple_annotation_plus6"]].reset_index(drop=True), sample_expr.reset_index(drop=True)], axis=1)
    make_feature_plot_panel(
        diagnostic_df,
        TREG_DIAGNOSTIC_GENES,
        FIGURE_DIR / "plus6_umap_treg_cd8_diagnostic.png",
        "plus6 Treg/CD8 diagnostic feature plots",
    )

    selected_clusters = sorted(
        set(change_df.loc[change_df["changed"], "leiden"].astype(str))
        | set(cluster_summary.loc[cluster_summary["simple_annotation_plus6"] == "Treg", "leiden"].astype(str))
        | set(change_df.loc[(change_df["changed"]) & (change_df["new_label"] == "CD8_T"), "leiden"].astype(str))
    )
    if selected_clusters:
        metric_cols = [
            "treg_marker_mean",
            "cd4_marker_mean",
            "cd8_marker_mean",
            "nk_marker_mean",
            "tcell_marker_mean",
            "paired_ab_fraction",
            "paired_gd_fraction",
        ]
        dot_df = cluster_summary.copy()
        dot_df["leiden"] = dot_df["leiden"].astype(str)
        dot_df = dot_df.set_index("leiden").reindex(selected_clusters)[metric_cols]
        x_positions = np.arange(len(metric_cols))
        y_positions = np.arange(len(selected_clusters))
        fig, ax = plt.subplots(figsize=(14, max(4, 0.55 * len(selected_clusters) + 2)), constrained_layout=True)
        vmax = float(np.percentile(dot_df.to_numpy().ravel(), 95.0)) if np.any(dot_df.to_numpy() > 0) else 1.0
        norm = plt.Normalize(vmin=0.0, vmax=vmax)
        for yi, cluster in enumerate(selected_clusters):
            for xi, metric in enumerate(metric_cols):
                ax.scatter(
                    xi,
                    yi,
                    s=160,
                    c=[plt.cm.viridis(norm(float(dot_df.loc[cluster, metric])))],
                    edgecolors="none",
                )
        ax.set_xticks(x_positions)
        ax.set_xticklabels(metric_cols, rotation=45, ha="right")
        ax.set_yticks(y_positions)
        ax.set_yticklabels(selected_clusters)
        ax.set_xlabel("Cluster summary metric")
        ax.set_ylabel("Leiden cluster")
        ax.set_title("plus6 Treg/CD8 cluster summary dotplot")
        sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
        cbar.set_label("Metric value")
        fig.savefig(FIGURE_DIR / "plus6_treg_cd8_cluster_dotplot.png", dpi=FIGURE_DPI, bbox_inches="tight")
        plt.close(fig)


def dataframe_to_html_table(df: pd.DataFrame, max_rows: int = 30) -> str:
    if df.shape[0] > max_rows:
        df = df.head(max_rows).copy()
    return df.to_html(index=False, classes="table table-sm", border=0)


def dataframe_to_markdown_fallback(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if max_rows is not None and df.shape[0] > max_rows:
        df = df.head(max_rows).copy()
    try:
        return df.to_markdown(index=False)
    except Exception:
        return "```text\n" + df.to_string(index=False) + "\n```"


def sample_plot_df(adata: ad.AnnData, sample_n: int = LEGEND_SAMPLE_LIMIT) -> pd.DataFrame:
    rng = np.random.default_rng(RANDOM_STATE)
    sample_n = min(sample_n, adata.n_obs)
    idx = np.sort(rng.choice(np.arange(adata.n_obs), size=sample_n, replace=False))
    plot_df = pd.DataFrame(
        {
            "umap1": adata.obsm["X_umap"][idx, 0],
            "umap2": adata.obsm["X_umap"][idx, 1],
            "source_gse_id": adata.obs["source_gse_id"].astype(str).to_numpy()[idx],
            "tissue": adata.obs.get("tissue_corrected", adata.obs.get("tissue", "")).astype(str).to_numpy()[idx],
            "leiden": adata.obs["leiden"].astype(str).to_numpy()[idx],
        }
    )
    if "simple_annotation_plus6" in adata.obs.columns:
        plot_df["simple_annotation_plus6"] = adata.obs["simple_annotation_plus6"].astype(str).to_numpy()[idx]
    if "phase4_trab_score" in adata.obs.columns:
        plot_df["phase4_trab_score"] = adata.obs["phase4_trab_score"].to_numpy(dtype=np.float32)[idx]
    if "phase4_trd_score" in adata.obs.columns:
        plot_df["phase4_trd_score"] = adata.obs["phase4_trd_score"].to_numpy(dtype=np.float32)[idx]
    if "phase4_trd_minus_trab" in adata.obs.columns:
        plot_df["phase4_trd_minus_trab"] = adata.obs["phase4_trd_minus_trab"].to_numpy(dtype=np.float32)[idx]
    if "has_TRA_TRB_paired" in adata.obs.columns:
        plot_df["has_TRA_TRB_paired"] = adata.obs["has_TRA_TRB_paired"].to_numpy(dtype=bool, copy=False)[idx]
    if "has_TRG_TRD_paired" in adata.obs.columns:
        plot_df["has_TRG_TRD_paired"] = adata.obs["has_TRG_TRD_paired"].to_numpy(dtype=bool, copy=False)[idx]
    return plot_df


def categorical_palette(categories: list[str]) -> list:
    n = len(categories)
    if n <= 10:
        return sns.color_palette("tab10", n)
    if n <= 20:
        return sns.color_palette("tab20", n)
    return sns.color_palette("husl", n)


def save_categorical_umap_with_legend(plot_df: pd.DataFrame, column: str, out_path: Path, title: str) -> None:
    values = plot_df[column].astype(str)
    categories = sorted(pd.unique(values))
    palette = dict(zip(categories, categorical_palette(categories), strict=False))
    fig, ax = plt.subplots(figsize=(14, 10), constrained_layout=False)
    sns.scatterplot(
        data=plot_df,
        x="umap1",
        y="umap2",
        hue=column,
        hue_order=categories,
        palette=palette,
        s=4,
        linewidth=0,
        ax=ax,
        rasterized=True,
        legend="full",
    )
    ax.set_title(title)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_aspect("equal", adjustable="box")
    legend = ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=True,
        borderaxespad=0.0,
        ncol=1 if len(categories) <= 15 else 2 if len(categories) <= 35 else 3,
        fontsize=7,
        title=column,
        title_fontsize=8,
        markerscale=2.5,
    )
    if legend is not None:
        for lh in legend.legend_handles:
            try:
                lh.set_alpha(1.0)
            except Exception:
                pass
    fig.tight_layout(rect=[0, 0, 0.82, 1])
    fig.savefig(out_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def horizontal_barplot(df: pd.DataFrame, category_col: str, count_col: str, title: str, out_path: Path, top_n: int = 30) -> None:
    plot_df = df.sort_values(count_col, ascending=False).head(top_n).copy()
    if plot_df.empty:
        return
    fig_h = max(6, 0.35 * plot_df.shape[0] + 1.5)
    fig, ax = plt.subplots(figsize=(12, fig_h), constrained_layout=True)
    sns.barplot(data=plot_df, y=category_col, x=count_col, color="#2F6690", ax=ax)
    ax.set_title(title)
    ax.set_xlabel("Cell count")
    ax.set_ylabel(category_col)
    fig.savefig(out_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def make_image_card(title: str, image_path: Path) -> str:
    rel = os.path.relpath(image_path, OUTPUT_ROOT)
    return (
        f"<section class='figure-card'><h3>{html.escape(title)}</h3>"
        f"<img src='{html.escape(rel)}' alt='{html.escape(title)}'></section>"
    )


def refresh_plus6_report_assets() -> None:
    required_paths = [
        PHASE3_SUMMARY_CSV,
        PREP_SUMMARY_CSV,
        PHASE4_SCORE_SUMMARY_CSV,
        PHASE4_GSE_SUMMARY_CSV,
        PHASE4_TISSUE_SUMMARY_CSV,
        ANNOTATION_CLUSTER_SUMMARY_CSV,
        ANNOTATION_CHANGES_CSV,
        ANNOTATION_GSE_COUNTS_CSV,
        ANNOTATION_TISSUE_COUNTS_CSV,
        FIGURE_DIR / "plus6_phase3_umap_by_gse.png",
        FIGURE_DIR / "plus6_phase3_umap_by_tissue.png",
        FIGURE_DIR / "plus6_phase3_umap_by_leiden.png",
        FIGURE_DIR / "plus6_umap_by_simple_annotation_corrected.png",
        FIGURE_DIR / "plus6_umap_treg_cd8_diagnostic.png",
        FIGURE_DIR / "plus6_treg_cd8_cluster_dotplot.png",
        FIGURE_DIR / "plus6_phase4_umap_trd_score.png",
        FIGURE_DIR / "plus6_phase4_umap_trab_score.png",
        FIGURE_DIR / "plus6_phase4_umap_trd_minus_trab.png",
        FIGURE_DIR / "plus6_tnk_marker_umap_panel.png",
        FIGURE_DIR / "plus6_phase4_score_distributions.png",
        FIGURE_DIR / "plus6_phase4_trab_vs_trd_raw.png",
        FIGURE_DIR / "plus6_phase4_trab_vs_trd_paired_tcr.png",
        FIGURE_DIR / "plus6_phase4_trd_over_trab_by_tissue_barplot.png",
        FIGURE_DIR / "plus6_phase4_trd_over_trab_by_gse_barplot.png",
        FIGURE_DIR / "plus6_phase4_trd_gt_0p1_by_tissue_barplot.png",
        FIGURE_DIR / "plus6_phase4_trd_gt_0p1_by_gse_barplot.png",
    ]
    if all(path.exists() for path in required_paths):
        logging.info("plus6 report assets already exist; skipping heavy asset refresh")
        return
    adata = ad.read_h5ad(PLUS6_INTEGRATED)
    if "X_umap" in adata.obsm and "leiden" in adata.obs.columns:
        write_phase3_figures(adata)
    if "simple_annotation_plus6" in adata.obs.columns:
        write_annotation_figures(adata)
    score_cols = set(PHASE4_SCORE_COLUMNS.values()) | set(PHASE4_SCALED_SCORE_COLUMNS.values())
    if score_cols.issubset(set(adata.obs.columns)):
        scores = {
            "tra": adata.obs[PHASE4_SCORE_COLUMNS["tra"]].to_numpy(dtype=np.float32),
            "trb": adata.obs[PHASE4_SCORE_COLUMNS["trb"]].to_numpy(dtype=np.float32),
            "trab": adata.obs[PHASE4_SCORE_COLUMNS["trab"]].to_numpy(dtype=np.float32),
            "trd": adata.obs[PHASE4_SCORE_COLUMNS["trd"]].to_numpy(dtype=np.float32),
            "trd_minus_trab": adata.obs[PHASE4_SCORE_COLUMNS["trd_minus_trab"]].to_numpy(dtype=np.float32),
            "trd_scaled": adata.obs[PHASE4_SCALED_SCORE_COLUMNS["trd_scaled"]].to_numpy(dtype=np.float32),
            "trab_scaled": adata.obs[PHASE4_SCALED_SCORE_COLUMNS["trab_scaled"]].to_numpy(dtype=np.float32),
            "trd_minus_trab_scaled": adata.obs[PHASE4_SCALED_SCORE_COLUMNS["trd_minus_trab_scaled"]].to_numpy(dtype=np.float32),
        }
        write_plus6_phase4_tables(adata, scores)
        write_plus6_phase4_figures(adata, scores)


def stage_report(send_email: bool) -> None:
    logging.info("Stage report: building plus6 HTML/PDF report")
    refresh_plus6_report_assets()
    subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "build_plus6_gdt_report_assets.py"),
        ],
        check=True,
    )
    phase3_summary = pd.read_csv(PHASE3_SUMMARY_CSV)
    prep_summary = pd.read_csv(PREP_SUMMARY_CSV)
    phase4_summary = pd.read_csv(PHASE4_SCORE_SUMMARY_CSV)
    gse_summary = pd.read_csv(PHASE4_GSE_SUMMARY_CSV)
    tissue_summary = pd.read_csv(PHASE4_TISSUE_SUMMARY_CSV)
    annotation_cluster = pd.read_csv(ANNOTATION_CLUSTER_SUMMARY_CSV)
    annotation_changes = pd.read_csv(ANNOTATION_CHANGES_CSV)
    annotation_gse = pd.read_csv(ANNOTATION_GSE_COUNTS_CSV)
    annotation_tissue = pd.read_csv(ANNOTATION_TISSUE_COUNTS_CSV)
    gdt_stats = pd.read_csv(TABLE_DIR / "plus6_gdt_candidate_statistics.csv")
    gdt_overlap = pd.read_csv(TABLE_DIR / "plus6_gdt_candidate_overlap_gt0p4.csv")
    gdt_paired_by_tissue = pd.read_csv(TABLE_DIR / "plus6_gdt_paired_gdtcr_by_tissue.csv")
    gdt_criteria_by_tissue = pd.read_csv(TABLE_DIR / "plus6_gdt_three_criteria_by_tissue.csv")

    md_lines = [
        "# plus6 integrated milestone report",
        "",
        "## Overview",
        "",
        f"- Base integrated milestone: `{BASE_INTEGRATED}`",
        f"- New datasets merged: `{len(NEW_INPUTS)}`",
        f"- Final plus6 cells: `{int(phase3_summary.loc[0, 'n_cells']):,}`",
        f"- Final plus6 GSEs: `{int(phase3_summary.loc[0, 'n_gses']):,}`",
        f"- Final plus6 Leiden clusters: `{int(phase3_summary.loc[0, 'n_leiden']):,}`",
        "",
        "## Input compatibility",
        "",
        dataframe_to_markdown_fallback(prep_summary),
        "",
        "## Phase 3",
        "",
        dataframe_to_markdown_fallback(phase3_summary),
        "",
        "## Phase 4 score summary",
        "",
        dataframe_to_markdown_fallback(phase4_summary),
        "",
        "## Top GSE composition",
        "",
        dataframe_to_markdown_fallback(gse_summary, max_rows=25),
        "",
        "## Tissue composition",
        "",
        dataframe_to_markdown_fallback(tissue_summary, max_rows=25),
        "",
        "## Annotation cluster summary",
        "",
        dataframe_to_markdown_fallback(annotation_cluster),
        "",
        "## Annotation changes",
        "",
        dataframe_to_markdown_fallback(annotation_changes, max_rows=40),
        "",
        "## Annotation by GSE",
        "",
        dataframe_to_markdown_fallback(annotation_gse, max_rows=60),
        "",
        "## Annotation by tissue",
        "",
        dataframe_to_markdown_fallback(annotation_tissue, max_rows=60),
        "",
        "## γδ-focused candidate statistics",
        "",
        dataframe_to_markdown_fallback(gdt_stats, max_rows=30),
        "",
        "## γδ-focused overlap breakdown",
        "",
        dataframe_to_markdown_fallback(gdt_overlap, max_rows=20),
        "",
        "## gdT cells with paired gdTCR by tissue",
        "",
        dataframe_to_markdown_fallback(gdt_paired_by_tissue, max_rows=40),
        "",
        "## gdT cells meeting at least one of the three criteria by tissue",
        "",
        dataframe_to_markdown_fallback(gdt_criteria_by_tissue, max_rows=40),
        "",
    ]
    PLUS6_REPORT_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    figures = [
        ("plus6 Phase 3 UMAP by GSE", FIGURE_DIR / "plus6_phase3_umap_by_gse.png"),
        ("plus6 Phase 3 UMAP by tissue", FIGURE_DIR / "plus6_phase3_umap_by_tissue.png"),
        ("plus6 Phase 3 UMAP by Leiden", FIGURE_DIR / "plus6_phase3_umap_by_leiden.png"),
        ("plus6 UMAP by simple annotation", FIGURE_DIR / "plus6_umap_by_simple_annotation.png"),
        ("plus6 UMAP by corrected simple annotation", FIGURE_DIR / "plus6_umap_by_simple_annotation_corrected.png"),
        ("plus6 Treg/CD8 diagnostic feature plots", FIGURE_DIR / "plus6_umap_treg_cd8_diagnostic.png"),
        ("plus6 Treg/CD8 cluster dotplot", FIGURE_DIR / "plus6_treg_cd8_cluster_dotplot.png"),
        ("plus6 paired TRA/TRB, paired TRG/TRD, and Sorted_gdT highlight UMAPs", FIGURE_DIR / "plus6_umap_paired_tcr_sorted_gdt.png"),
        ("plus6 paired TCR and doublet-proxy UMAPs", FIGURE_DIR / "plus6_umap_paired_tcr_doublets.png"),
        ("plus6 Sorted_gdT highlight UMAP", FIGURE_DIR / "plus6_umap_sorted_gdt_highlight.png"),
        ("plus6 UMAP by TRD score", FIGURE_DIR / "plus6_phase4_umap_trd_score.png"),
        ("plus6 UMAP by TRAB score", FIGURE_DIR / "plus6_phase4_umap_trab_score.png"),
        ("plus6 UMAP by TRD minus TRAB", FIGURE_DIR / "plus6_phase4_umap_trd_minus_trab.png"),
        ("plus6 TNK marker UMAP panel", FIGURE_DIR / "plus6_tnk_marker_umap_panel.png"),
        ("plus6 Phase 4 score distributions", FIGURE_DIR / "plus6_phase4_score_distributions.png"),
        ("plus6 Raw TRAB-versus-TRD score space", FIGURE_DIR / "plus6_phase4_trab_vs_trd_raw.png"),
        ("plus6 TRAB versus TRD colored by paired TCR", FIGURE_DIR / "plus6_phase4_trab_vs_trd_paired_tcr.png"),
        ("High-Confidence TRD-over-TRAB Candidates by Tissue", FIGURE_DIR / "plus6_phase4_trd_over_trab_by_tissue_barplot.png"),
        ("High-Confidence TRD-over-TRAB Candidates by GSE", FIGURE_DIR / "plus6_phase4_trd_over_trab_by_gse_barplot.png"),
        ("Broad TRD-Enriched Candidates by Tissue", FIGURE_DIR / "plus6_phase4_trd_gt_0p1_by_tissue_barplot.png"),
        ("Broad TRD-Enriched Candidates by GSE", FIGURE_DIR / "plus6_phase4_trd_gt_0p1_by_gse_barplot.png"),
    ]

    html_parts = [
        "<html><head><meta charset='utf-8'><title>plus6 integrated milestone report</title>",
        "<style>body{font-family:Arial,sans-serif;margin:24px;line-height:1.5} h1,h2,h3{margin-top:1.2em} .figure-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:24px} .figure-card img{width:100%;border:1px solid #ccc} table{border-collapse:collapse;font-size:12px} th,td{border:1px solid #ddd;padding:4px 6px;text-align:left}</style>",
        "</head><body>",
        "<h1>plus6 integrated milestone report</h1>",
        f"<p>Base current integrated milestone plus six new datasets, merged into a parallel lineage at <code>{html.escape(str(PLUS6_INTEGRATED))}</code>.</p>",
        "<h2>Key metrics</h2>",
        dataframe_to_html_table(phase3_summary),
        "<h2>Input compatibility</h2>",
        dataframe_to_html_table(prep_summary, max_rows=20),
        "<h2>Phase 4 score summary</h2>",
        dataframe_to_html_table(phase4_summary),
        "<h2>Annotation cluster summary</h2>",
        dataframe_to_html_table(annotation_cluster, max_rows=40),
        "<h2>Annotation changes</h2>",
        dataframe_to_html_table(annotation_changes, max_rows=40),
        "<h2>Top GSE composition</h2>",
        dataframe_to_html_table(gse_summary, max_rows=30),
        "<h2>Tissue composition</h2>",
        dataframe_to_html_table(tissue_summary, max_rows=30),
        "<h2>γδ-focused candidate statistics</h2>",
        dataframe_to_html_table(gdt_stats, max_rows=30),
        "<h2>γδ-focused overlap breakdown</h2>",
        dataframe_to_html_table(gdt_overlap, max_rows=20),
        "<h2>gdT cells with paired gdTCR by tissue</h2>",
        dataframe_to_html_table(gdt_paired_by_tissue, max_rows=40),
        "<h2>gdT cells meeting at least one of the three criteria by tissue</h2>",
        dataframe_to_html_table(gdt_criteria_by_tissue, max_rows=40),
        "<h2>Figures</h2>",
        "<div class='figure-grid'>",
    ]
    for title, path in figures:
        if path.exists():
            html_parts.append(make_image_card(title, path))
    html_parts.extend(["</div>", "</body></html>"])
    PLUS6_REPORT_HTML.write_text("\n".join(html_parts), encoding="utf-8")

    chrome_candidates = ["/usr/bin/google-chrome", "/usr/bin/google-chrome-stable"]
    chrome = next((Path(path) for path in chrome_candidates if Path(path).exists()), None)
    if chrome is None:
        raise FileNotFoundError("google-chrome not found for PDF export.")
    subprocess.run(
        [
            str(chrome),
            "--headless",
            "--disable-gpu",
            "--print-to-pdf-no-header",
            f"--print-to-pdf={PLUS6_REPORT_PDF}",
            str(PLUS6_REPORT_HTML),
        ],
        check=True,
    )

    if send_email:
        helper = Path("/home/tanlikai/.codex/skills/email-to-likai/scripts/send_email_to_likai.py")
        subprocess.run(
            [
                "python3",
                str(helper),
                "--subject",
                "[plus6] Integrated milestone profile report",
                "--report-file",
                str(PLUS6_REPORT_MD),
                "--attachment",
                str(PLUS6_REPORT_PDF),
            ],
            check=True,
        )


def main() -> None:
    args = parse_args()
    configure_logging()
    logging.info("Starting plus6 pipeline stage=%s", args.stage)

    if args.stage in {"prepare", "all"}:
        stage_prepare()
    if args.stage in {"integrate", "all"}:
        stage_integrate()
    if args.stage in {"score", "all"}:
        stage_score()
    if args.stage in {"annotate", "all"}:
        stage_annotate()
    if args.stage in {"report", "all"}:
        stage_report(send_email=args.send_email)

    logging.info("plus6 pipeline stage=%s completed", args.stage)


if __name__ == "__main__":
    main()
