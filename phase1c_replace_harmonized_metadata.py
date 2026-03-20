#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Phase 1c merged metadata backup and replacement.

This step exports the current merged candidate `adata.obs`, validates the
required join keys, joins the harmonized metadata by barcode plus GSE, writes a
backup of the previous metadata target, and replaces the canonical CSV only
after validation passes.
"""

from __future__ import annotations

import argparse
import logging
import os
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


# Config
INPUT_H5AD = Path("Integrated_dataset/TNK_candidates.h5ad")
HARMONIZED_METADATA = Path("analysis_26GSE_V4/outputs/harmonized_metadata_v4.csv")
BACKUP_METADATA = Path("analysis_26GSE_V4/outputs/metadata.csv.bk")
OUTPUT_ROOT = Path("Integrated_dataset")
TABLE_DIR = OUTPUT_ROOT / "tables"
LOG_DIR = OUTPUT_ROOT / "logs"

OBS_EXPORT_CSV = TABLE_DIR / "phase1c_merged_obs_export.csv.gz"
SUMMARY_CSV = TABLE_DIR / "phase1c_metadata_replacement_summary.csv"
JOIN_BY_GSE_CSV = TABLE_DIR / "phase1c_metadata_join_by_gse.csv"
QC_LOG_MD = LOG_DIR / "phase1c_metadata_replacement.md"

TEMP_METADATA = Path("analysis_26GSE_V4/outputs/harmonized_metadata_v4.csv.tmp")

REQUIRED_OBS_COLUMNS = [
    "original_cell_id",
    "source_gse_id",
    "phase1_sample_label",
]

OBS_JOIN_KEY = ["project name", "sampleid", "barcodes"]
RAW_JOIN_LEFT = ["source_gse_id", "original_cell_id"]
RAW_JOIN_RIGHT = ["gse_id", "cell_id"]
PSEUDO_MISSING_STRINGS = {"na", "nan", "none", "null", "<na>"}


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        default=INPUT_H5AD,
        help="Filtered merged candidate H5AD used as the driver for Phase 1c.",
    )
    parser.add_argument(
        "--metadata-csv",
        type=Path,
        default=HARMONIZED_METADATA,
        help="Existing harmonized metadata CSV to subset and replace.",
    )
    return parser.parse_args()


def setup_logging() -> None:
    """Configure logging."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )


def ensure_output_dirs() -> None:
    """Create Phase 1c output directories."""
    for path in (TABLE_DIR, LOG_DIR):
        path.mkdir(parents=True, exist_ok=True)


def normalize_string_series(series: pd.Series) -> pd.Series:
    """Normalize string-like metadata fields and collapse pseudo-missing tokens.

    The harmonized metadata contains some literal placeholder strings such as
    `NA`. These should be treated as missing values, not as real sample IDs.
    """
    normalized = series.fillna("").astype(str).str.strip()
    return normalized.mask(normalized.str.lower().isin(PSEUDO_MISSING_STRINGS), "")


def load_candidate_obs(input_h5ad: Path) -> pd.DataFrame:
    """Load the obs columns required for Phase 1c from the merged candidate."""
    logging.info("Loading candidate obs from %s", input_h5ad)
    adata = ad.read_h5ad(input_h5ad, backed="r")
    missing = [col for col in REQUIRED_OBS_COLUMNS if col not in adata.obs.columns]
    if missing:
        if getattr(adata, "file", None) is not None:
            adata.file.close()
        raise ValueError(f"Missing required obs columns for Phase 1c: {missing}")

    obs = adata.obs.copy()
    if getattr(adata, "file", None) is not None:
        adata.file.close()
    return obs


def build_obs_export(obs: pd.DataFrame) -> pd.DataFrame:
    """Build the exported Phase 1c obs table with explicit join columns."""
    export = obs.copy()
    export["project name"] = normalize_string_series(export["source_gse_id"])
    export["sampleid"] = normalize_string_series(export["phase1_sample_label"])
    export["barcodes"] = normalize_string_series(export["original_cell_id"])
    export.insert(0, "obs_name", export.index.astype(str))
    return export


def validate_required_columns(df: pd.DataFrame, columns: list[str], label: str) -> None:
    """Fail loudly if required columns are missing."""
    missing = [col for col in columns if col not in df.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def validate_unique_key(df: pd.DataFrame, key_cols: list[str], label: str) -> None:
    """Fail loudly if a join key is not unique."""
    duplicate_n = int(df.duplicated(key_cols).sum())
    if duplicate_n:
        raise ValueError(f"{label} has {duplicate_n} duplicated rows for key {key_cols}")


def load_harmonized_metadata(metadata_csv: Path) -> pd.DataFrame:
    """Load the existing harmonized metadata with stable string join columns."""
    logging.info("Loading harmonized metadata from %s", metadata_csv)
    metadata = pd.read_csv(
        metadata_csv,
        dtype={
            "cell_id": str,
            "gse_id": str,
            "sample_id": str,
        },
        low_memory=False,
    )
    metadata["cell_id"] = normalize_string_series(metadata["cell_id"])
    metadata["gse_id"] = normalize_string_series(metadata["gse_id"])
    if "sample_id" not in metadata.columns:
        raise ValueError("Existing harmonized metadata does not contain `sample_id`.")
    metadata["sample_id"] = normalize_string_series(metadata["sample_id"])
    return metadata


def join_metadata(obs_export: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Join the current obs export to harmonized metadata by GSE and barcode."""
    logging.info("Joining obs export to harmonized metadata by GSE + barcode")
    validate_required_columns(obs_export, RAW_JOIN_LEFT, "obs export")
    validate_required_columns(metadata, RAW_JOIN_RIGHT, "harmonized metadata")
    validate_unique_key(obs_export, RAW_JOIN_LEFT, "obs export")
    validate_unique_key(metadata, RAW_JOIN_RIGHT, "harmonized metadata")

    merged = obs_export.merge(
        metadata,
        left_on=RAW_JOIN_LEFT,
        right_on=RAW_JOIN_RIGHT,
        how="left",
        indicator=True,
        suffixes=("_obs", ""),
    )

    match_counts = merged["_merge"].value_counts()
    left_only = int(match_counts.get("left_only", 0))
    if left_only:
        missing_preview = merged.loc[
            merged["_merge"] != "both",
            ["source_gse_id", "original_cell_id"],
        ].head(20)
        raise ValueError(
            "Phase 1c join failed: "
            f"{left_only} candidate rows did not match harmonized metadata.\n"
            f"Preview:\n{missing_preview.to_string(index=False)}"
        )

    merged = merged.drop(columns=["_merge"])
    merged["sampleid"] = np.where(
        normalize_string_series(merged["sample_id"]) != "",
        normalize_string_series(merged["sample_id"]),
        normalize_string_series(merged["sampleid"]),
    )
    return merged


def reorder_columns(merged: pd.DataFrame) -> pd.DataFrame:
    """Move the explicit Phase 1c join columns to the front."""
    preferred = [
        "project name",
        "sampleid",
        "barcodes",
        "obs_name",
        "gse_id",
        "sample_id",
        "cell_id",
        "source_gse_id",
        "original_cell_id",
        "phase1_sample_label",
    ]
    ordered = [col for col in preferred if col in merged.columns]
    remainder = [col for col in merged.columns if col not in ordered]
    return merged[ordered + remainder]


def build_join_by_gse_table(merged: pd.DataFrame) -> pd.DataFrame:
    """Summarize joined metadata rows by GSE."""
    summary = (
        merged.groupby("source_gse_id")
        .size()
        .rename("rows_after_join")
        .reset_index()
        .sort_values(["rows_after_join", "source_gse_id"], ascending=[False, True])
        .reset_index(drop=True)
    )
    return summary


def write_obs_export(obs_export: pd.DataFrame) -> None:
    """Write the exported obs table for record keeping."""
    logging.info("Writing Phase 1c obs export to %s", OBS_EXPORT_CSV)
    obs_export.to_csv(OBS_EXPORT_CSV, index=False, compression="gzip")


def backup_existing_metadata(metadata_csv: Path) -> None:
    """Write the required backup copy before replacing the target."""
    if BACKUP_METADATA.exists():
        logging.info("Backup already exists at %s; preserving it", BACKUP_METADATA)
        return
    logging.info("Backing up %s to %s", metadata_csv, BACKUP_METADATA)
    shutil.copy2(metadata_csv, BACKUP_METADATA)


def write_replacement_metadata(merged: pd.DataFrame, metadata_csv: Path) -> None:
    """Write the replacement metadata to a temp file and atomically replace."""
    if TEMP_METADATA.exists():
        TEMP_METADATA.unlink()

    logging.info("Writing replacement metadata to temporary path %s", TEMP_METADATA)
    merged.to_csv(TEMP_METADATA, index=False)
    os.replace(TEMP_METADATA, metadata_csv)


def validate_written_metadata(
    metadata_csv: Path,
    expected_rows: int,
) -> dict[str, int]:
    """Re-read the replacement metadata and validate row counts and uniqueness."""
    logging.info("Validating written metadata at %s", metadata_csv)
    check = pd.read_csv(
        metadata_csv,
        dtype={
            "project name": str,
            "sampleid": str,
            "barcodes": str,
        },
        low_memory=False,
        usecols=["project name", "sampleid", "barcodes"],
    )
    check["project name"] = normalize_string_series(check["project name"])
    check["sampleid"] = normalize_string_series(check["sampleid"])
    check["barcodes"] = normalize_string_series(check["barcodes"])

    row_count = len(check)
    duplicate_n = int(check.duplicated(OBS_JOIN_KEY).sum())
    blank_project = int((check["project name"] == "").sum())
    blank_barcodes = int((check["barcodes"] == "").sum())
    blank_sampleid = int((check["sampleid"] == "").sum())

    if row_count != expected_rows:
        raise ValueError(
            f"Replacement metadata row count mismatch: expected {expected_rows}, got {row_count}"
        )
    if duplicate_n:
        raise ValueError(
            f"Replacement metadata has {duplicate_n} duplicated rows for key {OBS_JOIN_KEY}"
        )
    if blank_project or blank_barcodes:
        raise ValueError(
            "Replacement metadata has blank required join values: "
            f"project name={blank_project}, barcodes={blank_barcodes}"
        )

    return {
        "row_count": row_count,
        "duplicate_n": duplicate_n,
        "blank_project": blank_project,
        "blank_barcodes": blank_barcodes,
        "blank_sampleid": blank_sampleid,
    }


def count_reference_rows(csv_path: Path, reference_column: str) -> int:
    """Count rows in a CSV using a lightweight single-column read."""
    return int(pd.read_csv(csv_path, usecols=[reference_column]).shape[0])


def write_outputs(
    obs_export: pd.DataFrame,
    merged: pd.DataFrame,
    join_by_gse: pd.DataFrame,
    validation_stats: dict[str, int],
    metadata_csv: Path,
    metadata_rows_before: int,
) -> None:
    """Write Phase 1c summary outputs."""
    sampleid_blank_after_fill = int((merged["sampleid"].fillna("").astype(str) == "").sum())
    summary = pd.DataFrame(
        [
            {
                "candidate_rows": int(len(obs_export)),
                "replacement_rows": int(len(merged)),
                "metadata_rows_before": int(metadata_rows_before),
                "rows_removed_vs_previous_metadata": int(metadata_rows_before - len(merged)),
                "blank_sampleid_after_fill": sampleid_blank_after_fill,
                "written_blank_sampleid": int(validation_stats["blank_sampleid"]),
                "backup_path": str(BACKUP_METADATA),
                "target_path": str(metadata_csv),
            }
        ]
    )
    summary.to_csv(SUMMARY_CSV, index=False)
    join_by_gse.to_csv(JOIN_BY_GSE_CSV, index=False)

    lines = [
        "# Phase 1c Metadata Replacement",
        "",
        "## Scope",
        f"- Input candidate milestone: `{INPUT_H5AD}`",
        f"- Replaced metadata target: `{metadata_csv}`",
        f"- Backup file: `{BACKUP_METADATA}`",
        "- Join executed by `source_gse_id + original_cell_id` to `gse_id + cell_id` with explicit string typing.",
        "- Required output join columns were written as `project name`, `sampleid`, and `barcodes`.",
        "",
        "## Validation",
        f"- Candidate rows: {len(obs_export):,}",
        f"- Replacement rows: {len(merged):,}",
        f"- Previous metadata rows: {metadata_rows_before:,}",
        f"- Rows removed vs previous metadata: {metadata_rows_before - len(merged):,}",
        f"- Duplicate replacement join keys: {validation_stats['duplicate_n']:,}",
        f"- Blank `project name` after write: {validation_stats['blank_project']:,}",
        f"- Blank `barcodes` after write: {validation_stats['blank_barcodes']:,}",
        f"- Blank `sampleid` after write: {validation_stats['blank_sampleid']:,}",
        "",
        "## Outputs",
        f"- Obs export: `{OBS_EXPORT_CSV}`",
        f"- Summary table: `{SUMMARY_CSV}`",
        f"- GSE join summary: `{JOIN_BY_GSE_CSV}`",
        "",
        "## Conclusion",
        "- Phase 1c completed successfully.",
        "- The harmonized metadata target now corresponds exactly to the current filtered merged candidate object.",
        "- Phase 2 must still wait for explicit user approval.",
        "",
    ]
    QC_LOG_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Run Phase 1c metadata backup and replacement."""
    args = parse_args()
    setup_logging()
    ensure_output_dirs()

    obs = load_candidate_obs(args.input_h5ad)
    obs_export = build_obs_export(obs)
    validate_required_columns(obs_export, OBS_JOIN_KEY, "obs export")
    validate_unique_key(obs_export, OBS_JOIN_KEY, "obs export")
    write_obs_export(obs_export)

    metadata = load_harmonized_metadata(args.metadata_csv)
    metadata_rows_before = len(metadata)
    merged = join_metadata(obs_export, metadata)
    merged = reorder_columns(merged)
    validate_required_columns(merged, OBS_JOIN_KEY, "replacement metadata")
    validate_unique_key(merged, OBS_JOIN_KEY, "replacement metadata")

    join_by_gse = build_join_by_gse_table(merged)

    backup_existing_metadata(args.metadata_csv)
    write_replacement_metadata(merged, args.metadata_csv)
    validation_stats = validate_written_metadata(args.metadata_csv, expected_rows=len(obs_export))
    reference_rows = count_reference_rows(BACKUP_METADATA, reference_column="cell_id")
    write_outputs(
        obs_export=obs_export,
        merged=merged,
        join_by_gse=join_by_gse,
        validation_stats=validation_stats,
        metadata_csv=args.metadata_csv,
        metadata_rows_before=reference_rows,
    )

    logging.info(
        "Phase 1c completed: replacement rows=%s, blank sampleid after write=%s",
        len(merged),
        validation_stats["blank_sampleid"],
    )


if __name__ == "__main__":
    main()
