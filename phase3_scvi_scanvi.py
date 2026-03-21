#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Phase 3 scVI integration and post-scVI scANVI annotation.

This script consumes `Integrated_dataset/TNK_cleaned.h5ad`, attaches harmonized
metadata, trains scVI for integration, computes a RAPIDS-backed neighborhood
graph and UMAP, then maps the cleaned query to the saved blood-immune scANVI
reference as a coarse T/NK annotation prior.

Outputs:
- `Integrated_dataset/integrated.h5ad`
- Phase 3 tables under `Integrated_dataset/tables/`
- Phase 3 markdown summary under `Integrated_dataset/logs/`
- Phase 3 PNG figures under `Integrated_dataset/figures/`
"""

from __future__ import annotations

import argparse
import gc
import logging
import os
from pathlib import Path
import shutil
import sys

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import torch

# Import torch before rapids_singlecell. In this env, the import order avoids
# a CUDA runtime symbol-resolution failure inside libc10_cuda/libcudart.
import rapids_singlecell as rsc
import scvi
import seaborn as sns


# Config
INPUT_H5AD = Path("Integrated_dataset/TNK_cleaned.h5ad")
OUTPUT_H5AD = Path("Integrated_dataset/integrated.h5ad")
OUTPUT_ROOT = Path("Integrated_dataset")
TABLE_DIR = OUTPUT_ROOT / "tables"
FIGURE_DIR = OUTPUT_ROOT / "figures"
LOG_DIR = OUTPUT_ROOT / "logs"
MODEL_DIR = OUTPUT_ROOT / "models"
PHASE3_SCVI_MODEL_DIR = MODEL_DIR / "phase3_scvi_model"
PHASE3_HVG_GENES_TXT = TABLE_DIR / "phase3_hvg_genes.txt"


def resolve_high_speed_temp_root() -> Path:
    """Pick a writable local temp root while preserving one mirrored output tree.

    The Phase 3 bottleneck is dominated by large H5AD I/O. To keep the eventual
    migration back to NFS simple, stage large files under a parallel
    `Integrated_dataset/` tree on fast local storage instead of scattering temp
    files across unrelated paths.
    """
    env_override = os.environ.get("TNK_HIGH_SPEED_TEMP_ROOT", "").strip()
    candidates = []
    if env_override:
        candidates.append(Path(env_override))
    candidates.extend(
        [
            Path("/ssd/tnk_phase3"),
            Path("/dev/shm/tnk_phase3"),
            Path("/tmp/tnk_phase3"),
        ]
    )

    for candidate in candidates:
        try:
            candidate.mkdir(parents=True, exist_ok=True)
            probe = candidate / ".write_test"
            probe.write_text("ok", encoding="utf-8")
            probe.unlink()
            return candidate
        except OSError:
            continue

    raise RuntimeError(
        "Could not find a writable high-speed temp root. "
        "Set TNK_HIGH_SPEED_TEMP_ROOT to a writable local path."
    )


HIGH_SPEED_TEMP_ROOT = resolve_high_speed_temp_root()
LOCAL_OUTPUT_ROOT = HIGH_SPEED_TEMP_ROOT / OUTPUT_ROOT.name
LOCAL_STAGE_DIR = LOCAL_OUTPUT_ROOT
TEMP_OUTPUT_H5AD = LOCAL_OUTPUT_ROOT / "integrated.phase3_tmp.h5ad"
NFS_LINK_ROOT = Path("high_speed_temp")
NFS_LINK_OUTPUT_ROOT = NFS_LINK_ROOT / OUTPUT_ROOT.name
LOCAL_FINAL_OUTPUT_H5AD = LOCAL_OUTPUT_ROOT / OUTPUT_H5AD.name

MAIN_METADATA_CSV = Path("analysis_26GSE_V4/outputs/harmonized_metadata_v4.csv")
SUPP_METADATA_CSV = Path("analysis_26GSE_V4/outputs/harmonized_metadata_supp.csv")

REFERENCE_DIR = Path(
    "/home/tanlikai/databank/owndata/fasting/raw/report_from_niuxian/models/census_scanvi_ref_v1"
)
REFERENCE_MODEL_PT = REFERENCE_DIR / "model.pt"
REFERENCE_H5AD = REFERENCE_DIR / "census_reference_subset.h5ad"

N_HVGS = 4000
SCVI_MAX_EPOCHS = 20
SCANVI_MAX_EPOCHS = 20
BATCH_SIZE = 4096
LEIDEN_RESOLUTION = 0.8
QC_SAMPLE_MAX_CELLS = 300_000
SCANVI_SUBSET_MAX_CELLS = 300_000
SCANVI_MAX_PER_STRATUM = 2_000
MAX_COUNTS_PER_CELL = 100_000

PHASE3_SUMMARY_CSV = TABLE_DIR / "phase3_integration_summary.csv"
PHASE3_METADATA_CSV = TABLE_DIR / "phase3_metadata_join_summary.csv"
PHASE3_BATCH_CSV = TABLE_DIR / "phase3_batch_key_summary.csv"
PHASE3_LEIDEN_CSV = TABLE_DIR / "phase3_leiden_cluster_summary.csv"
PHASE3_LABEL_CSV = TABLE_DIR / "phase3_scanvi_label_summary.csv"
PHASE3_MARKER_CSV = TABLE_DIR / "phase3_marker_agreement_summary.csv"
PHASE3_HVG_FILTER_CSV = TABLE_DIR / "phase3_hvg_filter_summary.csv"
PHASE3_SCANVI_SUBSET_CSV = TABLE_DIR / "phase3_scanvi_subset_summary.csv"
PHASE3_CLUSTER_LABEL_CSV = TABLE_DIR / "phase3_cluster_label_summary.csv"
PHASE3_SANITIZATION_CSV = TABLE_DIR / "phase3_input_sanitization.csv"
PHASE3_QC_SAMPLE_CSV = TABLE_DIR / "phase3_qc_sample_obs.csv.gz"
PHASE3_QC_MD = LOG_DIR / "phase3_qc_summary.md"
PHASE3_SCVI_HISTORY_CSV = TABLE_DIR / "phase3_scvi_history.csv"
PHASE3_SCANVI_HISTORY_CSV = TABLE_DIR / "phase3_scanvi_history.csv"
PHASE3_RUN_LOG = LOG_DIR / "phase3_run.log"

PHASE3_UMAP_GSE_PNG = FIGURE_DIR / "phase3_umap_by_gse.png"
PHASE3_UMAP_BATCH_PNG = FIGURE_DIR / "phase3_umap_by_batch_level.png"
PHASE3_UMAP_LEIDEN_PNG = FIGURE_DIR / "phase3_umap_by_leiden.png"
PHASE3_UMAP_LABEL_PNG = FIGURE_DIR / "phase3_umap_by_scanvi_detailed_label.png"
PHASE3_UMAP_SUPER_PNG = FIGURE_DIR / "phase3_umap_by_scanvi_superclass.png"
PHASE3_CONFIDENCE_PNG = FIGURE_DIR / "phase3_scanvi_confidence_distribution.png"
PHASE3_MARKER_PNG = FIGURE_DIR / "phase3_marker_agreement.png"
PHASE3_LEIDEN_SIZE_PNG = FIGURE_DIR / "phase3_leiden_cluster_sizes.png"

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
    "technology_simple",
    "phase1_library_label",
    "phase1_sample_label",
    "scanvi_detailed_label",
    "scanvi_tnk_superclass",
    "scanvi_reference_other_flag",
    "scanvi_transfer_method",
]

METADATA_COLUMNS = [
    "source_gse_id",
    "original_cell_id",
    "project name",
    "sampleid",
    "barcodes",
    "sample_id",
    "library_id",
    "donor_id",
    "tissue",
    "technology_simple",
]

MARKER_GENES = ["CD3D", "TRAC", "TRDC", "NKG7", "KLRD1", "FCER1G", "EPCAM"]
LABEL_TRANSFER_METHOD = "subset_scANVI_plus_latent_centroid"


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        default=INPUT_H5AD,
        help="Path to the Phase 2-cleaned milestone.",
    )
    parser.add_argument(
        "--output-h5ad",
        type=Path,
        default=OUTPUT_H5AD,
        help="Path to the integrated Phase 3 milestone.",
    )
    parser.add_argument(
        "--n-hvgs",
        type=int,
        default=N_HVGS,
        help="Number of HVGs for scVI training.",
    )
    parser.add_argument(
        "--scvi-max-epochs",
        type=int,
        default=SCVI_MAX_EPOCHS,
        help="Maximum scVI epochs.",
    )
    parser.add_argument(
        "--scanvi-max-epochs",
        type=int,
        default=SCANVI_MAX_EPOCHS,
        help="Maximum scANVI query fine-tuning epochs.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=BATCH_SIZE,
        help="Mini-batch size for scVI/scANVI.",
    )
    parser.add_argument(
        "--qc-sample-max-cells",
        type=int,
        default=QC_SAMPLE_MAX_CELLS,
        help="Maximum cells to keep in the plotting sample.",
    )
    parser.add_argument(
        "--skip-scanvi",
        action="store_true",
        help="Pause scANVI annotation and finish integration, Leiden clustering, and UMAP only.",
    )
    parser.add_argument(
        "--resume-scanvi-only",
        action="store_true",
        help="Resume from the existing integrated SSD milestone and run only the scANVI/QC update pass.",
    )
    return parser.parse_args()


def setup_logging() -> None:
    """Configure concise logging."""
    handlers = [
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(PHASE3_RUN_LOG, mode="a", encoding="utf-8"),
    ]
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
        force=True,
    )


def ensure_output_dirs() -> None:
    """Create required output directories."""
    for path in (
        TABLE_DIR,
        FIGURE_DIR,
        LOG_DIR,
        MODEL_DIR,
        LOCAL_OUTPUT_ROOT,
        LOCAL_STAGE_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)


def ensure_high_speed_symlink_view() -> None:
    """Expose the SSD-backed mirrored tree from the NFS working directory.

    This keeps the fast local execution tree discoverable from the project root
    without changing the canonical NFS output locations.
    """
    NFS_LINK_ROOT.mkdir(parents=True, exist_ok=True)
    target = LOCAL_OUTPUT_ROOT.resolve()
    link_path = NFS_LINK_OUTPUT_ROOT

    if link_path.is_symlink():
        if link_path.resolve() == target:
            return
        link_path.unlink()
    elif link_path.exists():
        raise RuntimeError(
            f"Expected {link_path} to be absent or a symlink, but a real path already exists."
        )

    link_path.symlink_to(target, target_is_directory=True)


def release_memory(context: str) -> None:
    """Force a best-effort CPU/GPU memory release at phase boundaries."""
    collected = gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    logging.info("Released memory after %s (gc_collected=%s)", context, collected)


def stage_input_locally(input_h5ad: Path) -> Path:
    """Stage the large input H5AD to local disk when enough space is available."""
    local_path = LOCAL_OUTPUT_ROOT / input_h5ad.name
    source_size = input_h5ad.stat().st_size
    source_mtime = input_h5ad.stat().st_mtime

    if local_path.exists():
        local_stat = local_path.stat()
        if local_stat.st_size == source_size and local_stat.st_mtime >= source_mtime:
            logging.info("Reusing staged local input H5AD at %s", local_path)
            return local_path

    usage = shutil.disk_usage(LOCAL_STAGE_DIR)
    required_free = int(source_size * 1.15)
    if usage.free < required_free:
        logging.warning(
            "Skipping local staging because %s has only %.1f GB free and %.1f GB is required",
            LOCAL_STAGE_DIR,
            usage.free / (1024 ** 3),
            required_free / (1024 ** 3),
        )
        return input_h5ad

    temp_path = LOCAL_OUTPUT_ROOT / f"{input_h5ad.name}.tmp"
    if temp_path.exists():
        temp_path.unlink()
    logging.info(
        "Staging %s to local disk at %s before loading (size=%.1f GB)",
        input_h5ad,
        local_path,
        source_size / (1024 ** 3),
    )
    shutil.copy2(input_h5ad, temp_path)
    os.replace(temp_path, local_path)
    logging.info("Finished local staging to %s", local_path)
    return local_path


def make_obs_write_safe(adata: ad.AnnData) -> None:
    """Normalize text-like obs columns before H5AD write.

    The final write can fail when mixed Python objects remain in string-intended
    columns such as `sampleid`. Normalize these explicitly instead of relying on
    h5py to guess object conversions.
    """
    for column in WRITE_SAFE_TEXT_COLUMNS:
        if column in adata.obs.columns:
            normalized = normalize_text(pd.Series(adata.obs[column], index=adata.obs_names))
            adata.obs[column] = pd.Series(
                np.where(normalized.isna(), "", normalized.astype(str)),
                index=adata.obs_names,
                dtype=object,
            )

    for column in adata.obs.columns:
        series = adata.obs[column]
        if pd.api.types.is_categorical_dtype(series.dtype) or pd.api.types.is_string_dtype(series.dtype):
            normalized = normalize_text(pd.Series(series, index=adata.obs_names))
            adata.obs[column] = pd.Series(
                np.where(normalized.isna(), "", normalized.astype(str)),
                index=adata.obs_names,
                dtype=object,
            )


def normalize_text(series: pd.Series) -> pd.Series:
    """Normalize text columns while preserving real missing values."""
    series = series.astype("string")
    series = series.fillna("")
    series = series.str.strip()
    series = series.mask(series.eq(""), pd.NA)
    return series


def read_metadata_subset(path: Path) -> pd.DataFrame:
    """Read only the columns needed for Phase 3 metadata attachment."""
    if not path.exists():
        logging.warning("Metadata file not found and will be skipped: %s", path)
        return pd.DataFrame(columns=METADATA_COLUMNS)

    header = pd.read_csv(path, nrows=0)
    usecols = [column for column in METADATA_COLUMNS if column in header.columns]
    frame = pd.read_csv(path, usecols=usecols, dtype="string")
    for column in usecols:
        frame[column] = normalize_text(frame[column])
    for required in ["source_gse_id", "original_cell_id"]:
        if required not in frame.columns:
            raise ValueError(f"Required metadata key column missing in {path}: {required}")
    return frame


def load_combined_metadata() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load, combine, and validate the main and supplementary metadata files."""
    main = read_metadata_subset(MAIN_METADATA_CSV)
    supp = read_metadata_subset(SUPP_METADATA_CSV)
    metadata = pd.concat([main, supp], ignore_index=True)
    metadata["metadata_key"] = (
        metadata["source_gse_id"].fillna("missing_gse")
        + "||"
        + metadata["original_cell_id"].fillna("missing_cell")
    )

    duplicate_mask = metadata["metadata_key"].duplicated(keep=False)
    duplicates = metadata.loc[duplicate_mask, ["source_gse_id", "original_cell_id"]].copy()
    if duplicate_mask.any():
        raise ValueError(
            "Combined harmonized metadata contains duplicate source_gse_id + original_cell_id keys."
        )

    metadata = metadata.drop_duplicates("metadata_key").set_index("metadata_key")
    metadata_summary = pd.DataFrame(
        [
            {
                "rows_main": int(len(main)),
                "rows_supp": int(len(supp)),
                "rows_combined": int(len(metadata)),
                "duplicate_keys": int(len(duplicates)),
            }
        ]
    )
    return metadata, metadata_summary


def attach_metadata(adata: ad.AnnData, metadata: pd.DataFrame) -> pd.DataFrame:
    """Attach harmonized metadata to the cleaned milestone by GSE + barcode."""
    if "source_gse_id" not in adata.obs.columns or "original_cell_id" not in adata.obs.columns:
        raise ValueError("The cleaned milestone must contain source_gse_id and original_cell_id.")

    adata.obs["metadata_key"] = (
        adata.obs["source_gse_id"].astype(str).fillna("missing_gse")
        + "||"
        + adata.obs["original_cell_id"].astype(str).fillna("missing_cell")
    )
    join = metadata.reindex(adata.obs["metadata_key"])

    for column in metadata.columns:
        if column in adata.obs.columns:
            continue
        adata.obs[column] = join[column].to_numpy()

    matched = join.index.notna() & join.index.isin(metadata.index)
    summary = pd.DataFrame(
        [
            {
                "cells_in_adata": int(adata.n_obs),
                "metadata_matches": int(matched.sum()),
                "metadata_unmatched": int((~matched).sum()),
            }
        ]
    )
    return summary


def build_batch_key(adata: ad.AnnData) -> None:
    """Build a robust Phase 3 batch key with a clear fallback ladder."""
    def clean_obs(column: str) -> pd.Series:
        if column not in adata.obs.columns:
            return pd.Series(pd.NA, index=adata.obs_names, dtype="string")
        return normalize_text(pd.Series(adata.obs[column], index=adata.obs_names, dtype="string"))

    source_gse = clean_obs("source_gse_id").fillna("unknown_gse")
    library_id = clean_obs("library_id")
    sampleid = clean_obs("sampleid")
    sample_id = clean_obs("sample_id")
    phase1_library = clean_obs("phase1_library_label")
    phase1_sample = clean_obs("phase1_sample_label")

    batch_level = np.full(adata.n_obs, "gse", dtype=object)
    batch_key = source_gse.to_numpy(dtype=object, copy=True)

    library_mask = library_id.notna()
    batch_key[library_mask] = (
        source_gse[library_mask] + "|library|" + library_id[library_mask]
    ).to_numpy(dtype=object)
    batch_level[library_mask] = "library_id"

    sampleid_mask = (~library_mask) & sampleid.notna()
    batch_key[sampleid_mask] = (
        source_gse[sampleid_mask] + "|sampleid|" + sampleid[sampleid_mask]
    ).to_numpy(dtype=object)
    batch_level[sampleid_mask] = "sampleid"

    sample_id_mask = (~library_mask) & (~sampleid_mask) & sample_id.notna()
    batch_key[sample_id_mask] = (
        source_gse[sample_id_mask] + "|sample|" + sample_id[sample_id_mask]
    ).to_numpy(dtype=object)
    batch_level[sample_id_mask] = "sample_id"

    phase1_library_mask = (~library_mask) & (~sampleid_mask) & (~sample_id_mask) & phase1_library.notna()
    batch_key[phase1_library_mask] = (
        source_gse[phase1_library_mask] + "|phase1_library|" + phase1_library[phase1_library_mask]
    ).to_numpy(dtype=object)
    batch_level[phase1_library_mask] = "phase1_library_label"

    phase1_sample_mask = (
        (~library_mask)
        & (~sampleid_mask)
        & (~sample_id_mask)
        & (~phase1_library_mask)
        & phase1_sample.notna()
    )
    batch_key[phase1_sample_mask] = (
        source_gse[phase1_sample_mask] + "|phase1_sample|" + phase1_sample[phase1_sample_mask]
    ).to_numpy(dtype=object)
    batch_level[phase1_sample_mask] = "phase1_sample_label"

    adata.obs["phase3_batch_key"] = pd.Categorical(batch_key)
    adata.obs["phase3_batch_level"] = pd.Categorical(batch_level)


def sanitize_phase3_input(adata: ad.AnnData) -> tuple[ad.AnnData, pd.DataFrame]:
    """Drop the tiny set of numerically unsafe cells before scVI/scANVI.

    This is intentionally narrow. It removes only:
    - cells with NaN / Inf / negative values in the sparse matrix
    - cells with an extreme library size above the configured safety threshold
    """
    if not sp.issparse(adata.X):
        matrix = sp.csr_matrix(adata.X)
        adata.X = matrix
    else:
        matrix = adata.X.tocsr()
        adata.X = matrix

    row_sums = np.asarray(matrix.sum(axis=1)).ravel().astype(np.float64, copy=False)
    invalid_rows = np.zeros(adata.n_obs, dtype=bool)

    bad_value_mask = np.isnan(matrix.data) | np.isinf(matrix.data) | (matrix.data < 0)
    if bad_value_mask.any():
        bad_positions = np.where(bad_value_mask)[0]
        bad_rows = np.searchsorted(matrix.indptr, bad_positions, side="right") - 1
        invalid_rows[np.unique(bad_rows)] = True

    extreme_rows = row_sums > MAX_COUNTS_PER_CELL
    combined_bad_rows = invalid_rows | extreme_rows

    reason = np.full(adata.n_obs, "kept", dtype=object)
    reason[extreme_rows] = "extreme_library_size"
    reason[invalid_rows] = "invalid_matrix_values"
    reason[invalid_rows & extreme_rows] = "invalid_and_extreme"

    sanitization = pd.DataFrame(
        {
            "obs_name": adata.obs_names.astype(str),
            "source_gse_id": adata.obs["source_gse_id"].astype(str).to_numpy(),
            "phase3_sanitization_reason": reason,
            "phase3_total_counts": row_sums,
        }
    )
    removed = sanitization.loc[sanitization["phase3_sanitization_reason"] != "kept"].copy()
    removed.to_csv(PHASE3_SANITIZATION_CSV, index=False)

    if combined_bad_rows.any():
        logging.warning(
            "Phase 3 sanitization removed %s cells before scVI/scANVI (%s invalid, %s extreme count > %s)",
            int(combined_bad_rows.sum()),
            int(invalid_rows.sum()),
            int(extreme_rows.sum()),
            MAX_COUNTS_PER_CELL,
        )
        adata = adata[~combined_bad_rows].copy()
    else:
        logging.info("Phase 3 sanitization found no invalid or extreme-count cells.")

    return adata, removed


def build_hvg_exclusion_frame(var_names: pd.Index) -> pd.DataFrame:
    """Flag genes that should be excluded from HVG-driven clustering/UMAP."""
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


def select_hvg_genes(adata: ad.AnnData, n_hvgs: int) -> tuple[list[str], str, pd.DataFrame]:
    """Select allowed HVGs while excluding mitochondrial, ribosomal, and noncoding-like genes."""
    hvg_batch_key = "source_gse_id" if "source_gse_id" in adata.obs.columns else "phase3_batch_key"
    exclusion = build_hvg_exclusion_frame(adata.var_names)
    exclude_mask = exclusion["exclude_from_hvg"].to_numpy(dtype=bool, copy=False)
    allowed_mask = ~exclude_mask
    candidate_hvgs = min(adata.n_vars, max(n_hvgs * 3, n_hvgs + 5000))

    logging.info(
        "Selecting %s allowed HVGs for scVI training with primary method seurat_v3 and batch key %s",
        n_hvgs,
        hvg_batch_key,
    )
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
        logging.warning(
            "seurat_v3 HVG selection failed (%s). Falling back to cell_ranger HVGs.",
            exc,
        )
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=candidate_hvgs,
            flavor="cell_ranger",
            subset=False,
        )
        hvg_method = "cell_ranger_fallback"

    if "highly_variable_rank" in adata.var.columns:
        rank = pd.to_numeric(adata.var["highly_variable_rank"], errors="coerce")
    else:
        rank = pd.Series(np.nan, index=adata.var_names, dtype=np.float64)

    if rank.notna().sum() > 0:
        score = rank.copy()
        score[~allowed_mask] = np.inf
        selected = score.nsmallest(n_hvgs).index.astype(str).tolist()
    else:
        if "variances_norm" in adata.var.columns:
            score = pd.to_numeric(adata.var["variances_norm"], errors="coerce").fillna(-np.inf)
        elif "dispersions_norm" in adata.var.columns:
            score = pd.to_numeric(adata.var["dispersions_norm"], errors="coerce").fillna(-np.inf)
        else:
            raise ValueError("No HVG ranking metric found after highly_variable_genes.")
        score[~allowed_mask] = -np.inf
        selected = score.nlargest(n_hvgs).index.astype(str).tolist()

    if len(selected) < n_hvgs:
        raise ValueError(
            f"Unable to recover {n_hvgs} allowed HVGs after excluding mitochondrial/ribosomal/noncoding-like genes."
        )

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
    summary.to_csv(PHASE3_HVG_FILTER_CSV, index=False)
    return selected, hvg_method, summary


def train_scvi(
    adata: ad.AnnData,
    hvg_genes: list[str],
    n_hvgs: int,
    max_epochs: int,
    batch_size: int,
) -> tuple[scvi.model.SCVI, ad.AnnData]:
    """Select HVGs, train scVI, and return the model with the training AnnData."""
    adata_train = adata[:, hvg_genes].copy()
    PHASE3_HVG_GENES_TXT.write_text(
        "\n".join(adata_train.var_names.astype(str).tolist()) + "\n",
        encoding="utf-8",
    )

    logging.info(
        "Training scVI on shape (%s, %s) with %s batches",
        adata_train.n_obs,
        adata_train.n_vars,
        adata.obs["phase3_batch_key"].nunique(),
    )
    scvi.model.SCVI.setup_anndata(adata_train, batch_key="phase3_batch_key")
    model = scvi.model.SCVI(
        adata_train,
        n_hidden=128,
        n_layers=2,
        n_latent=30,
        gene_likelihood="nb",
    )
    model.train(
        max_epochs=max_epochs,
        batch_size=batch_size,
        accelerator="gpu",
        devices=1,
        early_stopping=True,
        early_stopping_patience=5,
    )
    return model, adata_train


def load_or_train_scvi(
    adata: ad.AnnData,
    n_hvgs: int,
    max_epochs: int,
    batch_size: int,
) -> tuple[scvi.model.SCVI, ad.AnnData, str, pd.DataFrame]:
    """Reuse a saved scVI model when possible; otherwise train a new one."""
    hvg_genes, hvg_method, hvg_filter_summary = select_hvg_genes(adata, n_hvgs)
    canonical_model_ready = PHASE3_SCVI_MODEL_DIR.exists()
    if PHASE3_SCVI_MODEL_DIR.exists() and PHASE3_HVG_GENES_TXT.exists():
        saved_hvg_genes = [
            line.strip()
            for line in PHASE3_HVG_GENES_TXT.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
    else:
        saved_hvg_genes = []

    if canonical_model_ready and saved_hvg_genes == hvg_genes:
        logging.info("Reusing canonical scVI model from %s", PHASE3_SCVI_MODEL_DIR)
        adata_train = adata[:, hvg_genes].copy()
        scvi.model.SCVI.setup_anndata(adata_train, batch_key="phase3_batch_key")
        model = scvi.model.SCVI.load(PHASE3_SCVI_MODEL_DIR, adata=adata_train)
        return model, adata_train, "reused_saved_model", hvg_filter_summary

    if canonical_model_ready and saved_hvg_genes != hvg_genes:
        logging.info(
            "Saved scVI HVGs do not match the current exclusion-aware HVG set. Retraining scVI."
        )

    model, adata_train = train_scvi(
        adata=adata,
        hvg_genes=hvg_genes,
        n_hvgs=n_hvgs,
        max_epochs=max_epochs,
        batch_size=batch_size,
    )
    model.save(PHASE3_SCVI_MODEL_DIR, overwrite=True, save_anndata=False)
    return model, adata_train, hvg_method, hvg_filter_summary


def load_saved_scvi_for_resume(
    adata: ad.AnnData,
) -> tuple[scvi.model.SCVI, ad.AnnData, int, str, pd.DataFrame]:
    """Load the saved scVI model for a scanVI-only resume pass.

    This path assumes the integrated milestone already carries the latent space
    and Leiden labels, so only the model registry and the historical HVG
    metadata need to be reloaded for QC bookkeeping.
    """
    if not PHASE3_SCVI_MODEL_DIR.exists():
        raise FileNotFoundError(f"Saved scVI model directory not found: {PHASE3_SCVI_MODEL_DIR}")
    if not PHASE3_HVG_GENES_TXT.exists():
        raise FileNotFoundError(f"Saved Phase 3 HVG gene list not found: {PHASE3_HVG_GENES_TXT}")

    saved_hvg_genes = [
        line.strip()
        for line in PHASE3_HVG_GENES_TXT.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    missing_hvgs = [gene for gene in saved_hvg_genes if gene not in adata.var_names]
    if missing_hvgs:
        raise ValueError(
            f"The integrated milestone is missing {len(missing_hvgs)} saved scVI HVGs; "
            "cannot resume scANVI-only safely."
        )

    adata_train = adata[:, saved_hvg_genes].copy()
    scvi.model.SCVI.setup_anndata(adata_train, batch_key="phase3_batch_key")
    model = scvi.model.SCVI.load(PHASE3_SCVI_MODEL_DIR, adata=adata_train)

    if PHASE3_HVG_FILTER_CSV.exists():
        hvg_filter_summary = pd.read_csv(PHASE3_HVG_FILTER_CSV)
    else:
        _, _, hvg_filter_summary = select_hvg_genes(adata, len(saved_hvg_genes))

    hvg_method = "resume_saved_model"
    if PHASE3_SUMMARY_CSV.exists():
        try:
            previous_summary = pd.read_csv(PHASE3_SUMMARY_CSV)
            if "hvg_method" in previous_summary.columns and len(previous_summary):
                hvg_method = str(previous_summary.iloc[0]["hvg_method"])
        except Exception:
            pass

    return model, adata_train, len(saved_hvg_genes), hvg_method, hvg_filter_summary


def attach_scanvi_labels_to_qc_sample(qc_sample: ad.AnnData, adata: ad.AnnData) -> None:
    """Propagate newly created scANVI label fields from the full object to the QC sample."""
    label_columns = [
        "scanvi_detailed_label",
        "scanvi_tnk_superclass",
        "scanvi_reference_other_flag",
        "scanvi_transfer_method",
        "scanvi_label_confidence",
    ]
    for column in label_columns:
        if column in adata.obs.columns:
            qc_sample.obs[column] = adata.obs.loc[qc_sample.obs_names, column].to_numpy()


def compute_scvi_latent(adata: ad.AnnData, model: scvi.model.SCVI, adata_train: ad.AnnData) -> None:
    """Store the scVI latent representation in the full object."""
    logging.info("Extracting scVI latent representation")
    adata.obsm["X_scVI"] = model.get_latent_representation(adata_train).astype(np.float32, copy=False)


def run_rapids_embedding(adata: ad.AnnData) -> None:
    """Compute neighbors, Leiden, and UMAP on a latent-only AnnData.

    `rapids_singlecell.get.anndata_to_GPU()` moves `X` to the GPU. To avoid
    sending the full count matrix to GPU memory, this function builds a compact
    AnnData whose `X` is the scVI latent representation only, runs the RAPIDS
    graph/UMAP workflow there, then copies `leiden` and `X_umap` back into the
    full milestone object.
    """
    logging.info("Running RAPIDS neighbors/Leiden/UMAP on the latent-only AnnData")
    latent_obs = pd.DataFrame(index=adata.obs_names.copy())
    latent_adata = ad.AnnData(X=adata.obsm["X_scVI"].copy(), obs=latent_obs)
    rsc.get.anndata_to_GPU(latent_adata)
    rsc.pp.neighbors(latent_adata, n_neighbors=15, metric="euclidean")
    rsc.tl.leiden(latent_adata, resolution=LEIDEN_RESOLUTION)
    rsc.tl.umap(latent_adata)
    rsc.get.anndata_to_CPU(latent_adata)

    adata.obs["leiden"] = latent_adata.obs["leiden"].astype("string").to_numpy()
    adata.obsm["X_umap"] = np.asarray(latent_adata.obsm["X_umap"]).astype(np.float32, copy=False)
    del latent_adata
    release_memory("RAPIDS embedding")


def get_reference_var_names() -> list[str]:
    """Recover the exact scANVI training genes from the saved model."""
    payload = torch.load(REFERENCE_MODEL_PT, map_location="cpu", weights_only=False)
    return payload["var_names"].tolist()


def load_reference_model() -> scvi.model.SCANVI:
    """Load the saved scANVI reference model with the correct training gene space."""
    ref_var_names = get_reference_var_names()
    ref = ad.read_h5ad(REFERENCE_H5AD)
    missing = [gene for gene in ref_var_names if gene not in ref.var_names]
    if missing:
        raise ValueError(
            f"Reference companion H5AD is missing {len(missing)} training genes from model.pt."
        )
    ref = ref[:, ref_var_names].copy()
    return scvi.model.SCANVI.load(REFERENCE_DIR, adata=ref)


def build_reference_symbol_map() -> pd.Series:
    """Build a unique gene_symbol -> Ensembl mapping for the reference model genes."""
    ref_var_names = get_reference_var_names()
    ref = ad.read_h5ad(REFERENCE_H5AD, backed="r")
    ref_var = ref.var.loc[ref_var_names, ["ensembl_id", "gene_symbol"]].copy()
    if getattr(ref, "file", None) is not None:
        ref.file.close()

    ref_var["gene_symbol"] = normalize_text(pd.Series(ref_var["gene_symbol"], dtype="string"))
    ref_var["ensembl_id"] = normalize_text(pd.Series(ref_var["ensembl_id"], dtype="string"))
    ref_var = ref_var.dropna(subset=["gene_symbol", "ensembl_id"])
    ref_var = ref_var.reset_index(drop=True)
    ref_var = ref_var.drop_duplicates(subset=["gene_symbol"], keep="first")
    return pd.Series(
        ref_var["ensembl_id"].to_numpy(),
        index=ref_var["gene_symbol"].to_numpy(),
        dtype="string",
    )


def collapse_scanvi_label(label: str) -> str:
    """Collapse detailed reference labels into a coarse T/NK superclass."""
    value = (label or "").strip().lower()
    if not value:
        return "reference_other"
    if "natural killer" in value or value == "nk cell":
        return "NK_cell"
    t_like_terms = [
        "t cell",
        "t-cell",
        "alpha-beta t",
        "gamma-delta t",
        "gd t",
        "cytotoxic t",
        "helper t",
        "regulatory t",
        "memory t",
        "naive t",
    ]
    if any(term in value for term in t_like_terms):
        return "T_cell"
    return "reference_other"


def build_scanvi_subset_indices(
    adata: ad.AnnData,
    max_cells: int,
    max_per_stratum: int,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Select a stratified scANVI subset across GSE x Leiden strata."""
    gse_codes, gse_labels = pd.factorize(adata.obs["source_gse_id"].astype(str), sort=False)
    leiden_codes, leiden_labels = pd.factorize(adata.obs["leiden"].astype(str), sort=False)
    pairs = pd.DataFrame({"gse_code": gse_codes, "leiden_code": leiden_codes})
    counts = pairs.value_counts(sort=False).reset_index(name="cells")
    counts["target_cells"] = np.floor(max_cells * counts["cells"] / counts["cells"].sum()).astype(int)
    counts["target_cells"] = counts["target_cells"].clip(lower=1, upper=max_per_stratum)
    counts["target_cells"] = np.minimum(counts["target_cells"], counts["cells"]).astype(int)

    if counts["target_cells"].sum() > max_cells:
        scale = max_cells / counts["target_cells"].sum()
        counts["target_cells"] = np.floor(counts["target_cells"] * scale).astype(int).clip(lower=1)
        counts["target_cells"] = np.minimum(counts["target_cells"], counts["cells"]).astype(int)

    overflow = int(counts["target_cells"].sum() - max_cells)
    if overflow > 0:
        reducible = counts.sort_values("target_cells", ascending=False).index.tolist()
        for idx in reducible:
            if overflow <= 0:
                break
            if counts.at[idx, "target_cells"] > 1:
                counts.at[idx, "target_cells"] -= 1
                overflow -= 1

    order = np.lexsort((leiden_codes, gse_codes))
    gse_sorted = gse_codes[order]
    leiden_sorted = leiden_codes[order]
    change = np.where((np.diff(gse_sorted) != 0) | (np.diff(leiden_sorted) != 0))[0] + 1
    starts = np.r_[0, change]
    ends = np.r_[change, len(order)]
    target_map = {
        (int(row.gse_code), int(row.leiden_code)): int(row.target_cells)
        for row in counts.itertuples(index=False)
    }

    rng = np.random.default_rng(0)
    selected_parts: list[np.ndarray] = []
    for start, end in zip(starts, ends):
        idx = order[start:end]
        key = (int(gse_sorted[start]), int(leiden_sorted[start]))
        target = target_map.get(key, 0)
        if target <= 0:
            continue
        if len(idx) <= target:
            selected_parts.append(idx)
        else:
            selected_parts.append(np.sort(rng.choice(idx, size=target, replace=False)))

    selected_idx = np.sort(np.concatenate(selected_parts))
    summary = counts.copy()
    summary["source_gse_id"] = gse_labels.take(summary["gse_code"]).astype(str)
    summary["leiden"] = leiden_labels.take(summary["leiden_code"]).astype(str)
    summary["selection_fraction"] = summary["target_cells"] / summary["cells"]
    summary = summary[
        ["source_gse_id", "leiden", "cells", "target_cells", "selection_fraction"]
    ].sort_values(["target_cells", "cells"], ascending=[False, False]).reset_index(drop=True)
    summary.to_csv(PHASE3_SCANVI_SUBSET_CSV, index=False)
    return selected_idx, summary


def transfer_labels_by_centroid(
    full_latent: np.ndarray,
    subset_idx: np.ndarray,
    subset_labels: np.ndarray,
    subset_confidence: np.ndarray,
    full_obs_names: pd.Index,
) -> tuple[pd.Series, np.ndarray, pd.Series]:
    """Transfer scANVI labels from a subset to all cells by nearest latent centroid."""
    subset_idx = np.asarray(subset_idx, dtype=np.int64)
    subset_labels = np.asarray(subset_labels, dtype=object)
    subset_confidence = np.asarray(subset_confidence, dtype=np.float32)
    subset_latent = np.asarray(full_latent[subset_idx], dtype=np.float32)

    if subset_latent.shape[0] != subset_labels.shape[0] or subset_labels.shape[0] != subset_confidence.shape[0]:
        raise ValueError(
            "scANVI subset transfer inputs must have matching lengths: "
            f"subset_latent={subset_latent.shape[0]}, subset_labels={subset_labels.shape[0]}, "
            f"subset_confidence={subset_confidence.shape[0]}"
        )

    valid_mask = pd.Series(subset_labels, dtype="string").notna().to_numpy()
    subset_idx_valid = subset_idx[valid_mask]
    subset_latent = subset_latent[valid_mask]
    valid_labels = pd.Series(subset_labels[valid_mask], dtype="string").reset_index(drop=True)
    valid_confidence = np.asarray(subset_confidence[valid_mask], dtype=np.float32)

    if subset_latent.shape[0] != valid_labels.shape[0] or valid_labels.shape[0] != valid_confidence.shape[0]:
        raise ValueError(
            "Filtered scANVI subset transfer inputs must have matching lengths: "
            f"subset_latent={subset_latent.shape[0]}, valid_labels={valid_labels.shape[0]}, "
            f"valid_confidence={valid_confidence.shape[0]}"
        )

    label_frame = pd.DataFrame(
        {
            "label": valid_labels,
            "confidence": valid_confidence,
        }
    )
    unique_labels = label_frame["label"].drop_duplicates().astype(str).tolist()
    if not unique_labels:
        raise ValueError("scANVI subset labeling produced no valid labels.")

    centroids = []
    for label in unique_labels:
        mask = label_frame["label"].astype(str).eq(label).to_numpy()
        weights = np.clip(label_frame.loc[mask, "confidence"].to_numpy(dtype=np.float32), 1e-4, None)
        centroid = np.average(subset_latent[mask], axis=0, weights=weights)
        centroids.append(np.asarray(centroid, dtype=np.float32))
    centroids = np.vstack(centroids).astype(np.float32, copy=False)
    centroid_norm = np.sum(centroids * centroids, axis=1)

    assigned = np.empty(full_latent.shape[0], dtype=object)
    assigned_confidence = np.zeros(full_latent.shape[0], dtype=np.float32)
    transfer_method = np.full(full_latent.shape[0], "latent_centroid_transfer", dtype=object)
    chunk_size = 250_000

    for start in range(0, full_latent.shape[0], chunk_size):
        end = min(start + chunk_size, full_latent.shape[0])
        chunk = np.asarray(full_latent[start:end], dtype=np.float32)
        chunk_norm = np.sum(chunk * chunk, axis=1, keepdims=True)
        distances = chunk_norm + centroid_norm[None, :] - 2.0 * np.matmul(chunk, centroids.T)
        distances = np.maximum(distances, 0.0)
        best = distances.argmin(axis=1)
        assigned[start:end] = np.asarray(unique_labels, dtype=object)[best]

        if distances.shape[1] == 1:
            assigned_confidence[start:end] = 1.0
        else:
            top2 = np.partition(distances, kth=1, axis=1)[:, :2]
            nearest = top2.min(axis=1)
            second = top2.max(axis=1)
            confidence = 1.0 - (nearest / np.clip(second, 1e-6, None))
            assigned_confidence[start:end] = confidence.astype(np.float32, copy=False)

    assigned[subset_idx_valid] = label_frame["label"].astype(str).to_numpy()
    assigned_confidence[subset_idx_valid] = label_frame["confidence"].to_numpy(dtype=np.float32, copy=False)
    transfer_method[subset_idx_valid] = "direct_subset_scANVI"
    label_series = pd.Series(assigned, index=full_obs_names, dtype="string")
    method_series = pd.Series(transfer_method, index=full_obs_names, dtype="string")
    return label_series, assigned_confidence, method_series


def run_scanvi_mapping(
    adata: ad.AnnData,
    max_epochs: int,
    batch_size: int,
) -> tuple[scvi.model.SCANVI, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Map a stratified subset with scANVI, then transfer labels to all cells in latent space."""
    logging.info("Loading saved scANVI reference model")
    reference_model = load_reference_model()
    reference_var_names = pd.Index(reference_model.adata.var_names)
    reference_symbol_map = build_reference_symbol_map()
    query_symbols = pd.Index(adata.var_names.astype(str))
    matched_symbols = query_symbols.intersection(reference_symbol_map.index)
    mapped_ensembl = pd.Index(reference_symbol_map.loc[matched_symbols].astype(str))
    common_var_names = mapped_ensembl.intersection(reference_var_names)
    logging.info(
        "Preparing scANVI query with %s / %s reference genes present after symbol-to-Ensembl mapping",
        len(common_var_names),
        len(reference_var_names),
    )
    if len(common_var_names) == 0:
        raise ValueError(
            "No reference var names found after mapping query gene symbols to the reference Ensembl IDs."
        )

    matched_symbols = [symbol for symbol in matched_symbols if reference_symbol_map[symbol] in common_var_names]
    subset_idx, subset_summary = build_scanvi_subset_indices(
        adata=adata,
        max_cells=SCANVI_SUBSET_MAX_CELLS,
        max_per_stratum=SCANVI_MAX_PER_STRATUM,
    )
    logging.info(
        "Running scANVI on a stratified subset of %s cells across %s GSE x Leiden strata",
        len(subset_idx),
        len(subset_summary),
    )

    query_subset = adata[subset_idx, matched_symbols].copy()
    query_subset.var["query_gene_symbol"] = query_subset.var_names.astype(str)
    query_subset.var["ensembl_id"] = reference_symbol_map.loc[matched_symbols].astype(str).to_numpy()
    query_subset.var_names = pd.Index(query_subset.var["ensembl_id"].astype(str))
    if "batch" not in query_subset.obs.columns:
        query_subset.obs["batch"] = pd.Categorical(np.repeat("query_batch", query_subset.n_obs))
    duplicate_mask = query_subset.var_names.duplicated()
    if duplicate_mask.any():
        query_subset = query_subset[:, ~duplicate_mask].copy()

    scvi.model.SCANVI.prepare_query_anndata(query_subset, reference_model)
    query_model = scvi.model.SCANVI.load_query_data(query_subset, reference_model)
    query_model.train(
        max_epochs=max_epochs,
        batch_size=batch_size,
        accelerator="gpu",
        devices=1,
        plan_kwargs={"weight_decay": 0.0},
        early_stopping=True,
        early_stopping_patience=5,
    )

    subset_detailed = query_model.predict()
    subset_soft = query_model.predict(soft=True)
    subset_confidence = np.asarray(subset_soft.max(axis=1), dtype=np.float32)
    detailed_series, full_confidence, transfer_method = transfer_labels_by_centroid(
        full_latent=np.asarray(adata.obsm["X_scVI"], dtype=np.float32),
        subset_idx=subset_idx,
        subset_labels=np.asarray(subset_detailed, dtype=object),
        subset_confidence=subset_confidence,
        full_obs_names=adata.obs_names,
    )
    superclass = detailed_series.map(collapse_scanvi_label).fillna("reference_other")

    adata.obs["scanvi_detailed_label"] = pd.Series(
        detailed_series.fillna("").astype(str).to_numpy(),
        index=adata.obs_names,
        dtype=object,
    )
    adata.obs["scanvi_label_confidence"] = full_confidence.astype(np.float32, copy=False)
    adata.obs["scanvi_tnk_superclass"] = pd.Series(
        superclass.fillna("").astype(str).to_numpy(),
        index=adata.obs_names,
        dtype=object,
    )
    adata.obs["scanvi_reference_other_flag"] = pd.Series(
        np.where(superclass.eq("reference_other"), "reference_other", "tnk_like"),
        index=adata.obs_names,
        dtype=object,
    )
    adata.obs["scanvi_transfer_method"] = pd.Series(
        transfer_method.fillna("").astype(str).to_numpy(),
        index=adata.obs_names,
        dtype=object,
    )

    label_summary = (
        pd.DataFrame(
            {
                "scanvi_detailed_label": detailed_series,
                "scanvi_tnk_superclass": superclass,
                "scanvi_label_confidence": full_confidence,
            }
        )
        .groupby(["scanvi_detailed_label", "scanvi_tnk_superclass"], dropna=False)
        .agg(
            cells=("scanvi_label_confidence", "size"),
            mean_confidence=("scanvi_label_confidence", "mean"),
        )
        .reset_index()
        .sort_values(["cells", "scanvi_detailed_label"], ascending=[False, True])
        .reset_index(drop=True)
    )
    cluster_summary = (
        adata.obs.assign(
            scanvi_detailed_label=adata.obs["scanvi_detailed_label"].astype(str),
            scanvi_tnk_superclass=adata.obs["scanvi_tnk_superclass"].astype(str),
            leiden=adata.obs["leiden"].astype(str),
            scanvi_transfer_method=adata.obs["scanvi_transfer_method"].astype(str),
        )
        .groupby(["leiden", "scanvi_detailed_label", "scanvi_tnk_superclass"], dropna=False)
        .size()
        .reset_index(name="cells")
        .sort_values(["leiden", "cells"], ascending=[True, False])
    )
    cluster_top = cluster_summary.drop_duplicates("leiden", keep="first").rename(
        columns={
            "scanvi_detailed_label": "dominant_scanvi_detailed_label",
            "scanvi_tnk_superclass": "dominant_scanvi_tnk_superclass",
            "cells": "dominant_label_cells",
        }
    )
    cluster_sizes = (
        adata.obs["leiden"].astype(str).value_counts().rename_axis("leiden").reset_index(name="cluster_cells")
    )
    cluster_top = cluster_top.merge(cluster_sizes, on="leiden", how="left")
    cluster_top["dominant_label_fraction"] = (
        cluster_top["dominant_label_cells"] / cluster_top["cluster_cells"]
    )
    cluster_top.to_csv(PHASE3_CLUSTER_LABEL_CSV, index=False)

    del subset_soft
    del subset_detailed
    del superclass
    del query_subset
    del reference_model
    del reference_symbol_map
    release_memory("scANVI mapping outputs")
    return query_model, label_summary, subset_summary, cluster_top


def build_qc_sample(adata: ad.AnnData, max_cells: int) -> ad.AnnData:
    """Create a representative QC plotting sample with proportional GSE sampling."""
    if adata.n_obs <= max_cells:
        sample = adata.copy()
        sample.obs["phase3_qc_sampled"] = True
        return sample

    gse_counts = adata.obs["source_gse_id"].astype(str).value_counts()
    target_counts = np.maximum(
        1,
        np.floor(max_cells * (gse_counts / gse_counts.sum())).astype(int),
    )
    selected: list[np.ndarray] = []
    for gse_id, target in target_counts.items():
        idx = np.where(adata.obs["source_gse_id"].astype(str).to_numpy() == gse_id)[0]
        if len(idx) <= target:
            selected.append(idx)
        else:
            selected.append(np.random.default_rng(0).choice(idx, size=target, replace=False))

    selected_idx = np.sort(np.concatenate(selected))
    if len(selected_idx) > max_cells:
        selected_idx = np.sort(
            np.random.default_rng(1).choice(selected_idx, size=max_cells, replace=False)
        )

    sample = adata[selected_idx].copy()
    sample.obs["phase3_qc_sampled"] = True
    return sample


def collapse_top_categories(series: pd.Series, top_n: int = 15, other_label: str = "other") -> pd.Series:
    """Collapse rare categories for cleaner UMAP legends."""
    series = series.astype(str)
    top = series.value_counts().head(top_n).index
    collapsed = np.where(series.isin(top), series, other_label)
    return pd.Series(collapsed, index=series.index, dtype="string")


def plot_umap(sample: ad.AnnData, color: str, out_path: Path, title: str, categorical: bool = True) -> None:
    """Plot a QC-sample UMAP to PNG."""
    fig = sc.pl.umap(
        sample,
        color=color,
        show=False,
        return_fig=True,
        frameon=False,
        title=title,
    )
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_leiden_summary(adata: ad.AnnData) -> pd.DataFrame:
    """Summarize Leiden cluster sizes."""
    summary = (
        adata.obs["leiden"]
        .astype(str)
        .value_counts()
        .rename_axis("leiden")
        .reset_index(name="cells")
        .sort_values("cells", ascending=False)
        .reset_index(drop=True)
    )
    summary["fraction"] = summary["cells"] / float(adata.n_obs)
    summary.to_csv(PHASE3_LEIDEN_CSV, index=False)
    return summary


def plot_leiden_sizes(leiden_summary: pd.DataFrame) -> None:
    """Plot the Leiden cluster size distribution."""
    top = leiden_summary.head(40).copy()
    fig, ax = plt.subplots(figsize=(12, 5.6))
    sns.barplot(data=top, x="leiden", y="cells", ax=ax, color="#3a718d")
    ax.set_title("Phase 3 Leiden Cluster Sizes")
    ax.set_xlabel("Leiden cluster")
    ax.set_ylabel("Cells")
    ax.tick_params(axis="x", rotation=90)
    fig.tight_layout()
    fig.savefig(PHASE3_LEIDEN_SIZE_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_history_table(history: object, out_path: Path) -> None:
    """Persist scvi-tools history objects whether they are DataFrames or dict-like."""
    if history is None:
        pd.DataFrame(columns=["metric", "step", "value"]).to_csv(out_path, index=False)
        return

    if isinstance(history, pd.DataFrame):
        history.to_csv(out_path, index=True)
        return

    if isinstance(history, dict):
        frames: list[pd.DataFrame] = []
        for metric, values in history.items():
            if isinstance(values, pd.DataFrame):
                frame = values.copy().reset_index().rename(columns={"index": "step"})
                if "value" not in frame.columns and frame.shape[1] >= 2:
                    value_col = [col for col in frame.columns if col != "step"][0]
                    frame = frame.rename(columns={value_col: "value"})
                frame.insert(0, "metric", str(metric))
                frames.append(frame)
                continue

            if isinstance(values, pd.Series):
                series = values.reset_index()
                series.columns = ["step", "value"]
                series.insert(0, "metric", str(metric))
                frames.append(series)
                continue

            array = np.asarray(values)
            if array.ndim == 0:
                frame = pd.DataFrame({"metric": [str(metric)], "step": [0], "value": [array.item()]})
            elif array.ndim == 1:
                frame = pd.DataFrame(
                    {"metric": str(metric), "step": np.arange(array.shape[0]), "value": array}
                )
            else:
                frame = pd.DataFrame(array)
                frame.insert(0, "metric", str(metric))
            frames.append(frame)

        if frames:
            pd.concat(frames, ignore_index=True).to_csv(out_path, index=False)
        else:
            pd.DataFrame(columns=["metric", "step", "value"]).to_csv(out_path, index=False)
        return

    pd.DataFrame({"history_repr": [repr(history)]}).to_csv(out_path, index=False)


def build_marker_agreement(sample: ad.AnnData) -> pd.DataFrame:
    """Summarize key marker expression across scANVI coarse classes."""
    group_labels = sample.obs["scanvi_tnk_superclass"].astype(str).to_numpy()
    result_rows: list[dict[str, object]] = []
    total_counts = np.asarray(sample.X.sum(axis=1)).ravel().astype(np.float64, copy=False)

    for gene in MARKER_GENES:
        if gene not in sample.var_names:
            continue
        gene_idx = sample.var_names.get_loc(gene)
        values = np.asarray(sample.X[:, gene_idx].toarray()).ravel() if sp.issparse(sample.X) else sample.X[:, gene_idx]
        normalized = np.log1p(1e4 * values / np.clip(total_counts, 1.0, None))
        frame = pd.DataFrame(
            {
                "scanvi_tnk_superclass": group_labels,
                "marker_gene": gene,
                "log1p_cpm": normalized,
            }
        )
        summary = (
            frame.groupby(["scanvi_tnk_superclass", "marker_gene"], dropna=False)["log1p_cpm"]
            .mean()
            .reset_index()
        )
        result_rows.append(summary)

    if not result_rows:
        return pd.DataFrame(columns=["scanvi_tnk_superclass", "marker_gene", "log1p_cpm"])
    return pd.concat(result_rows, ignore_index=True)


def plot_marker_agreement(marker_summary: pd.DataFrame) -> None:
    """Plot a marker agreement heatmap."""
    pivot = marker_summary.pivot(
        index="scanvi_tnk_superclass",
        columns="marker_gene",
        values="log1p_cpm",
    ).fillna(0)
    fig, ax = plt.subplots(figsize=(8.5, 3.8))
    sns.heatmap(pivot, cmap="viridis", ax=ax)
    ax.set_title("Phase 3 Marker Agreement By scANVI Superclass")
    ax.set_xlabel("Marker gene")
    ax.set_ylabel("scANVI superclass")
    fig.tight_layout()
    fig.savefig(PHASE3_MARKER_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_qc_tables(
    adata: ad.AnnData,
    metadata_summary: pd.DataFrame,
    metadata_join_summary: pd.DataFrame,
    leiden_summary: pd.DataFrame,
    label_summary: pd.DataFrame,
    subset_summary: pd.DataFrame,
    cluster_summary: pd.DataFrame,
    marker_summary: pd.DataFrame,
    sanitization_removed: pd.DataFrame,
    hvg_filter_summary: pd.DataFrame,
    scvi_model: scvi.model.SCVI,
    scanvi_model: scvi.model.SCANVI | None,
    n_hvgs: int,
    hvg_method: str,
    scanvi_executed: bool,
) -> pd.DataFrame:
    """Write Phase 3 tables and return the top-level summary table."""
    batch_summary = (
        adata.obs["phase3_batch_level"]
        .astype(str)
        .value_counts()
        .rename_axis("phase3_batch_level")
        .reset_index(name="cells")
    )
    batch_summary.to_csv(PHASE3_BATCH_CSV, index=False)
    leiden_summary.to_csv(PHASE3_LEIDEN_CSV, index=False)
    metadata_summary.to_csv(PHASE3_METADATA_CSV, index=False)
    label_summary.to_csv(PHASE3_LABEL_CSV, index=False)
    marker_summary.to_csv(PHASE3_MARKER_CSV, index=False)
    subset_summary.to_csv(PHASE3_SCANVI_SUBSET_CSV, index=False)
    cluster_summary.to_csv(PHASE3_CLUSTER_LABEL_CSV, index=False)

    save_history_table(scvi_model.history, PHASE3_SCVI_HISTORY_CSV)
    if scanvi_model is not None:
        save_history_table(scanvi_model.history, PHASE3_SCANVI_HISTORY_CSV)

    summary = pd.DataFrame(
        [
            {
                "cells": int(adata.n_obs),
                "genes": int(adata.n_vars),
                "leiden_clusters": int(adata.obs["leiden"].astype(str).nunique()),
                "hvgs": int(n_hvgs),
                "hvg_method": hvg_method,
                "hvg_excluded_genes_total": int(hvg_filter_summary.iloc[0]["genes_excluded_total"]),
                "hvg_excluded_mitochondrial": int(hvg_filter_summary.iloc[0]["genes_excluded_mitochondrial"]),
                "hvg_excluded_ribosomal": int(hvg_filter_summary.iloc[0]["genes_excluded_ribosomal"]),
                "hvg_excluded_noncoding_like": int(hvg_filter_summary.iloc[0]["genes_excluded_noncoding_like"]),
                "phase3_sanitized_cells_removed": int(len(sanitization_removed)),
                "unique_phase3_batches": int(adata.obs["phase3_batch_key"].nunique()),
                "unique_gses": int(adata.obs["source_gse_id"].astype(str).nunique()),
                "scanvi_executed": bool(scanvi_executed),
                "scanvi_subset_cells": int(subset_summary["target_cells"].sum()) if len(subset_summary) else 0,
                "scanvi_subset_strata": int(len(subset_summary)),
                "label_transfer_method": LABEL_TRANSFER_METHOD if scanvi_executed else "not_run",
                "scanvi_detailed_labels": int(adata.obs["scanvi_detailed_label"].astype(str).nunique()) if "scanvi_detailed_label" in adata.obs.columns else 0,
                "scanvi_superclasses": int(adata.obs["scanvi_tnk_superclass"].astype(str).nunique()) if "scanvi_tnk_superclass" in adata.obs.columns else 0,
                "reference_other_cells": int(
                    (adata.obs["scanvi_tnk_superclass"].astype(str) == "reference_other").sum()
                ) if "scanvi_tnk_superclass" in adata.obs.columns else 0,
                "mean_scanvi_confidence": float(
                    pd.to_numeric(adata.obs["scanvi_label_confidence"], errors="coerce").mean()
                ) if "scanvi_label_confidence" in adata.obs.columns else np.nan,
                "metadata_matches": int(metadata_join_summary.iloc[0]["metadata_matches"]),
                "metadata_unmatched": int(metadata_join_summary.iloc[0]["metadata_unmatched"]),
            }
        ]
    )
    summary.to_csv(PHASE3_SUMMARY_CSV, index=False)
    return summary


def write_qc_summary(summary: pd.DataFrame, label_summary: pd.DataFrame, scanvi_executed: bool) -> None:
    """Write the Phase 3 markdown QC summary."""
    row = summary.iloc[0]
    top_labels = label_summary.head(10)
    lines = [
        "# Phase 3 QC Summary",
        "",
        "## Scope",
        f"- Input milestone: `{INPUT_H5AD}`",
        f"- Canonical output milestone target: `{OUTPUT_H5AD}`",
        f"- Current validated high-speed output path: `{LOCAL_FINAL_OUTPUT_H5AD}`",
        f"- Execution env: `rapids_sc_py310`",
        "- RAPIDS import note: `torch` was imported before `rapids_singlecell` to avoid the CUDA symbol-resolution issue in this env",
        "- Figures were written in PNG format only",
        "",
        "## Integration summary",
        f"- Cells integrated: {int(row['cells']):,}",
        f"- Genes retained in the integrated milestone: {int(row['genes']):,}",
        f"- Leiden clusters: {int(row['leiden_clusters']):,}",
        f"- HVGs used for scVI: {int(row['hvgs']):,}",
        f"- HVG method: {row['hvg_method']}",
        f"- HVG exclusions before clustering/UMAP: {int(row['hvg_excluded_genes_total']):,} total "
        f"(mitochondrial={int(row['hvg_excluded_mitochondrial']):,}, ribosomal={int(row['hvg_excluded_ribosomal']):,}, "
        f"noncoding_like={int(row['hvg_excluded_noncoding_like']):,})",
        f"- Phase 3 preflight cells removed for numerical safety: {int(row['phase3_sanitized_cells_removed']):,}",
        f"- Unique Phase 3 batches: {int(row['unique_phase3_batches']):,}",
        f"- Unique GSEs represented: {int(row['unique_gses']):,}",
        "",
    ]
    if scanvi_executed:
        lines.extend(
            [
                "## scANVI annotation summary",
                f"- scANVI subset cells: {int(row['scanvi_subset_cells']):,}",
                f"- scANVI subset strata: {int(row['scanvi_subset_strata']):,}",
                f"- Label transfer method: `{row['label_transfer_method']}`",
                f"- Detailed labels observed: {int(row['scanvi_detailed_labels']):,}",
                f"- Coarse superclasses observed: {int(row['scanvi_superclasses']):,}",
                f"- `reference_other` cells: {int(row['reference_other_cells']):,}",
                f"- Mean scANVI confidence: {row['mean_scanvi_confidence']:.4f}",
                "",
            ]
        )
    else:
        lines.extend(
            [
                "## scANVI annotation summary",
                "- scANVI was paused by user instruction for this run.",
                "- `integrated.h5ad` currently contains scVI integration, Leiden clustering, and UMAP only.",
                "",
            ]
        )

    lines.extend(
        [
        "## Metadata attachment",
        f"- Metadata-matched cells: {int(row['metadata_matches']):,}",
        f"- Metadata-unmatched cells: {int(row['metadata_unmatched']):,}",
        "",
            ]
    )
    if scanvi_executed:
        lines.append("## Top detailed labels")
        for _, label_row in top_labels.iterrows():
            lines.append(
                f"- {label_row['scanvi_detailed_label']}: {int(label_row['cells']):,} cells "
                f"(superclass={label_row['scanvi_tnk_superclass']}, "
                f"mean_confidence={label_row['mean_confidence']:.4f})"
            )

    lines.extend(
        [
            "",
            "## Outputs",
            f"- Summary table: `{PHASE3_SUMMARY_CSV}`",
            f"- Metadata join summary: `{PHASE3_METADATA_CSV}`",
            f"- Batch summary: `{PHASE3_BATCH_CSV}`",
            f"- Leiden summary: `{PHASE3_LEIDEN_CSV}`",
            f"- Label summary: `{PHASE3_LABEL_CSV}`",
            f"- Marker agreement summary: `{PHASE3_MARKER_CSV}`",
            f"- HVG filter summary: `{PHASE3_HVG_FILTER_CSV}`",
            f"- scANVI subset summary: `{PHASE3_SCANVI_SUBSET_CSV}`",
            f"- Cluster label summary: `{PHASE3_CLUSTER_LABEL_CSV}`",
            f"- Input sanitization table: `{PHASE3_SANITIZATION_CSV}`",
            f"- Figures: `{PHASE3_UMAP_GSE_PNG}`, `{PHASE3_UMAP_BATCH_PNG}`, `{PHASE3_UMAP_LEIDEN_PNG}`, `{PHASE3_LEIDEN_SIZE_PNG}`"
            + (f", `{PHASE3_UMAP_LABEL_PNG}`, `{PHASE3_UMAP_SUPER_PNG}`, `{PHASE3_CONFIDENCE_PNG}`, `{PHASE3_MARKER_PNG}`" if scanvi_executed else ""),
            "",
            "## QC conclusion",
            "- The integrated milestone currently carries scVI integration, Leiden clustering, and UMAP results.",
            "- If scANVI is resumed later, its labels should be interpreted as a coarse blood-immune prior rather than final tissue-state truth.",
            "",
        ]
    )
    PHASE3_QC_MD.write_text("\n".join(lines), encoding="utf-8")


def finalize_phase3_output(temp_h5ad: Path, output_h5ad: Path) -> Path:
    """Finalize the validated Phase 3 output without auto-migrating back to NFS.

    During high-speed execution, keep the validated H5AD in the mirrored local
    tree and migrate it back to the canonical NFS path only on explicit user
    instruction.
    """
    local_output_h5ad = LOCAL_OUTPUT_ROOT / output_h5ad.name
    logging.info(
        "Keeping validated Phase 3 output in the mirrored high-speed tree at %s; "
        "canonical NFS path %s will be updated only on explicit user instruction",
        local_output_h5ad,
        output_h5ad,
    )
    if temp_h5ad.resolve() == local_output_h5ad.resolve():
        logging.info("Validated Phase 3 output is already at the final mirrored SSD path: %s", local_output_h5ad)
        return local_output_h5ad
    if local_output_h5ad.exists():
        local_output_h5ad.unlink()
    os.replace(temp_h5ad, local_output_h5ad)
    return local_output_h5ad


def main() -> None:
    """Run Phase 3 integration and annotation."""
    args = parse_args()
    setup_logging()
    ensure_output_dirs()
    ensure_high_speed_symlink_view()
    logging.info(
        "Using high-speed temp root %s with mirrored local output tree %s",
        HIGH_SPEED_TEMP_ROOT,
        LOCAL_OUTPUT_ROOT,
    )

    torch.set_float32_matmul_precision("medium")
    np.random.seed(0)

    metadata, metadata_summary = load_combined_metadata()
    if args.resume_scanvi_only:
        resume_path = LOCAL_FINAL_OUTPUT_H5AD if LOCAL_FINAL_OUTPUT_H5AD.exists() else args.output_h5ad
        logging.info("Resuming scANVI-only from existing integrated milestone at %s", resume_path)
        adata = ad.read_h5ad(resume_path)
        if "X_scVI" not in adata.obsm:
            raise ValueError("Resume scanVI-only mode requires `X_scVI` in `adata.obsm`.")
        if "X_umap" not in adata.obsm:
            raise ValueError("Resume scanVI-only mode requires `X_umap` in `adata.obsm`.")
        if "leiden" not in adata.obs.columns:
            raise ValueError("Resume scanVI-only mode requires `leiden` in `adata.obs`.")

        metadata_join_summary = attach_metadata(adata, metadata)
        build_batch_key(adata)
        scvi_model, adata_train, effective_n_hvgs, hvg_method, hvg_filter_summary = load_saved_scvi_for_resume(adata)
        del adata_train
        sanitization_removed = (
            pd.read_csv(PHASE3_SANITIZATION_CSV) if PHASE3_SANITIZATION_CSV.exists()
            else pd.DataFrame(columns=["obs_name", "source_gse_id", "phase3_sanitization_reason", "phase3_total_counts"])
        )
        leiden_summary = build_leiden_summary(adata)
    else:
        staged_input_h5ad = stage_input_locally(args.input_h5ad)
        logging.info("Loading cleaned milestone from %s", staged_input_h5ad)
        adata = ad.read_h5ad(staged_input_h5ad)
        adata, sanitization_removed = sanitize_phase3_input(adata)
        metadata_join_summary = attach_metadata(adata, metadata)
        build_batch_key(adata)

        scvi_model, adata_train, hvg_method, hvg_filter_summary = load_or_train_scvi(
            adata=adata,
            n_hvgs=args.n_hvgs,
            max_epochs=args.scvi_max_epochs,
            batch_size=args.batch_size,
        )
        effective_n_hvgs = args.n_hvgs
        compute_scvi_latent(adata=adata, model=scvi_model, adata_train=adata_train)
        del adata_train
        release_memory("scVI latent extraction")
        run_rapids_embedding(adata)
        leiden_summary = build_leiden_summary(adata)

    qc_sample = build_qc_sample(adata, max_cells=args.qc_sample_max_cells)
    qc_sample.obs["plot_gse_top"] = collapse_top_categories(qc_sample.obs["source_gse_id"], top_n=18)
    qc_sample.obs["plot_phase3_batch_level"] = collapse_top_categories(
        qc_sample.obs["phase3_batch_level"], top_n=12
    )
    qc_sample.obs["plot_leiden_top"] = collapse_top_categories(qc_sample.obs["leiden"], top_n=24)

    scanvi_model: scvi.model.SCANVI | None = None
    if args.skip_scanvi:
        label_summary = pd.DataFrame(
            columns=["scanvi_detailed_label", "scanvi_tnk_superclass", "cells", "mean_confidence"]
        )
        subset_summary = pd.DataFrame(
            columns=["source_gse_id", "leiden", "cells", "target_cells", "selection_fraction"]
        )
        cluster_summary = leiden_summary.rename(
            columns={
                "cells": "cluster_cells",
            }
        ).assign(
            dominant_scanvi_detailed_label="not_run",
            dominant_scanvi_tnk_superclass="not_run",
            dominant_label_cells=0,
            dominant_label_fraction=0.0,
        )[
            [
                "leiden",
                "dominant_scanvi_detailed_label",
                "dominant_scanvi_tnk_superclass",
                "dominant_label_cells",
                "cluster_cells",
                "dominant_label_fraction",
            ]
        ]
    else:
        scanvi_model, label_summary, subset_summary, cluster_summary = run_scanvi_mapping(
            adata=adata,
            max_epochs=args.scanvi_max_epochs,
            batch_size=args.batch_size,
        )
        attach_scanvi_labels_to_qc_sample(qc_sample, adata)
        qc_sample.obs["plot_scanvi_label_top"] = collapse_top_categories(
            qc_sample.obs["scanvi_detailed_label"], top_n=18
        )

    qc_sample.obs.to_csv(PHASE3_QC_SAMPLE_CSV, index=True, compression="gzip")

    plot_umap(
        qc_sample,
        color="plot_gse_top",
        out_path=PHASE3_UMAP_GSE_PNG,
        title="Phase 3 QC Sample UMAP By GSE",
    )
    plot_umap(
        qc_sample,
        color="plot_phase3_batch_level",
        out_path=PHASE3_UMAP_BATCH_PNG,
        title="Phase 3 QC Sample UMAP By Batch Level",
    )
    plot_umap(
        qc_sample,
        color="plot_leiden_top",
        out_path=PHASE3_UMAP_LEIDEN_PNG,
        title="Phase 3 QC Sample UMAP By Leiden Cluster",
    )
    plot_leiden_sizes(leiden_summary)

    if args.skip_scanvi:
        marker_summary = pd.DataFrame(columns=["scanvi_tnk_superclass", "marker_gene", "log1p_cpm"])
    else:
        plot_umap(
            qc_sample,
            color="plot_scanvi_label_top",
            out_path=PHASE3_UMAP_LABEL_PNG,
            title="Phase 3 QC Sample UMAP By scANVI Detailed Label",
        )
        plot_umap(
            qc_sample,
            color="scanvi_tnk_superclass",
            out_path=PHASE3_UMAP_SUPER_PNG,
            title="Phase 3 QC Sample UMAP By scANVI T/NK Superclass",
        )

        fig, ax = plt.subplots(figsize=(7.2, 4.8))
        ax.hist(
            pd.to_numeric(adata.obs["scanvi_label_confidence"], errors="coerce").fillna(0).to_numpy(),
            bins=60,
            color="#3a718d",
            alpha=0.9,
        )
        ax.set_xlabel("scANVI label confidence")
        ax.set_ylabel("Cells")
        ax.set_title("Phase 3 scANVI Confidence Distribution")
        fig.tight_layout()
        fig.savefig(PHASE3_CONFIDENCE_PNG, dpi=300, bbox_inches="tight")
        plt.close(fig)

        marker_summary = build_marker_agreement(qc_sample)
        plot_marker_agreement(marker_summary)

    summary = write_qc_tables(
        adata=adata,
        metadata_summary=metadata_summary,
        metadata_join_summary=metadata_join_summary,
        leiden_summary=leiden_summary,
        label_summary=label_summary,
        subset_summary=subset_summary,
        cluster_summary=cluster_summary,
        marker_summary=marker_summary,
        sanitization_removed=sanitization_removed,
        hvg_filter_summary=hvg_filter_summary,
        scvi_model=scvi_model,
        scanvi_model=scanvi_model,
        n_hvgs=effective_n_hvgs,
        hvg_method=hvg_method,
        scanvi_executed=not args.skip_scanvi,
    )
    write_qc_summary(summary=summary, label_summary=label_summary, scanvi_executed=not args.skip_scanvi)

    if TEMP_OUTPUT_H5AD.exists():
        TEMP_OUTPUT_H5AD.unlink()
    make_obs_write_safe(adata)
    adata.write_h5ad(TEMP_OUTPUT_H5AD, convert_strings_to_categoricals=False)
    final_output_h5ad = finalize_phase3_output(TEMP_OUTPUT_H5AD, args.output_h5ad)

    logging.info(
        "Phase 3 complete: cells=%s genes=%s batches=%s final_output=%s",
        int(summary.iloc[0]["cells"]),
        int(summary.iloc[0]["genes"]),
        int(summary.iloc[0]["unique_phase3_batches"]),
        final_output_h5ad,
    )


if __name__ == "__main__":
    try:
        main()
    except Exception:
        logging.exception("Phase 3 failed")
        raise
