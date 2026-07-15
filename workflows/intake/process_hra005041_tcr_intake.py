#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Standardize HRA005041 metadata and integrate productive ab/gd TCR fields.

This intake step treats `newdata/HRA005041/HRA005041_T_cells_subset.h5ad` as the
source of truth and writes a project-format per-dataset H5AD under
`downloads/per_gse_h5ad_with_metadata/`.

Key rules:
- join TCR by `sample_id + barcode_core`, never barcode alone
- use productive alpha-beta chains only
- use productive-like gamma-delta chains only, approximated from MiXCR clone
  rows that have a valid TRG/TRD chain call plus amino-acid and nucleotide CDR3
- keep metadata fields aligned with the project harmonized schema
"""

from __future__ import annotations

# TNK_WORKFLOW_BOOTSTRAP
import sys as _tnk_sys
from pathlib import Path as _TnkPath

_TNK_PROJECT_ROOT = _TnkPath(__file__).resolve().parents[2]
for _tnk_path in (
    _TNK_PROJECT_ROOT,
    _TNK_PROJECT_ROOT / "src",
    _TNK_PROJECT_ROOT / "workflows" / "integration",
    _TNK_PROJECT_ROOT / "workflows" / "intake",
    _TNK_PROJECT_ROOT / "workflows" / "metadata",
    _TNK_PROJECT_ROOT / "workflows" / "analysis",
    _TNK_PROJECT_ROOT / "workflows" / "gdtai",
    _TNK_PROJECT_ROOT / "workflows" / "gdt_atlas",
    _TNK_PROJECT_ROOT / "workflows" / "reporting",
    _TNK_PROJECT_ROOT / "workflows" / "maintenance",
):
    _tnk_value = str(_tnk_path)
    if _tnk_value not in _tnk_sys.path:
        _tnk_sys.path.insert(0, _tnk_value)


import argparse
import csv
import logging
import os
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
INPUT_H5AD = PROJECT_ROOT / "newdata" / "HRA005041" / "HRA005041_T_cells_subset.h5ad"
AB_DIR = PROJECT_ROOT / "newdata" / "HRA005041" / "abTCR"
GD_DIR = PROJECT_ROOT / "newdata" / "HRA005041" / "gdTCR"
OUTPUT_H5AD = PROJECT_ROOT / "downloads" / "per_gse_h5ad_with_metadata" / "HRA005041_T_cells_subset.h5ad"
TABLE_DIR = PROJECT_ROOT / "Integrated_dataset" / "tables" / "HRA005041_tcr_intake"
LOG_DIR = PROJECT_ROOT / "Integrated_dataset" / "logs"
SUMMARY_CSV = TABLE_DIR / "HRA005041_tcr_join_summary.csv"
SAMPLE_SUMMARY_CSV = TABLE_DIR / "HRA005041_tcr_sample_summary.csv"
UNMATCHED_CSV = TABLE_DIR / "HRA005041_tcr_unmatched_summary.csv"
QC_MD = LOG_DIR / "HRA005041_tcr_intake.md"
LOG_FILE = LOG_DIR / "HRA005041_tcr_intake.log"

DATASET_ID = "HRA005041"
RAW_ASSAY = "10x 5'"
ASSAY_TYPE = "10x 5' gene expression + VDJ"
SAMPLE_TYPE = "T_cell_subset"
ENRICHMENT = "T_cell_subset"

AB_REQUIRED_COLS = [
    "barcode",
    "is_cell",
    "high_confidence",
    "chain",
    "v_gene",
    "d_gene",
    "j_gene",
    "c_gene",
    "productive",
    "cdr3",
    "cdr3_nt",
    "reads",
    "umis",
    "raw_clonotype_id",
    "cellbarcode",
    "Samplename",
]

CANONICAL_TCR_COLS = [
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


def reshape_chain_rows_wide(
    collapsed: pd.DataFrame,
    *,
    field_map: dict[str, str],
    chains: list[str],
) -> pd.DataFrame:
    """Convert one-row-per-chain data into one-row-per-cell wide format."""
    join_cols = ["sample_id", "barcode_core"]
    wide: pd.DataFrame | None = None
    for chain in chains:
        subset = collapsed.loc[collapsed["chain"].eq(chain), join_cols + list(field_map)].copy()
        if subset.empty:
            continue
        subset = subset.rename(columns={src: f"{chain}_{dst}" for src, dst in field_map.items()})
        subset = subset.drop(columns=["chain"], errors="ignore")
        ensure_unique_key(subset, join_cols, f"{chain} collapsed chain table")
        wide = subset if wide is None else wide.merge(subset, on=join_cols, how="outer", validate="one_to_one")

    if wide is None:
        return pd.DataFrame(columns=join_cols)
    ensure_unique_key(wide, join_cols, "wide chain table")
    return wide


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, default=INPUT_H5AD)
    parser.add_argument("--ab-dir", type=Path, default=AB_DIR)
    parser.add_argument("--gd-dir", type=Path, default=GD_DIR)
    parser.add_argument("--output-h5ad", type=Path, default=OUTPUT_H5AD)
    return parser.parse_args()


def configure_logging() -> None:
    """Configure logging to file and stdout."""
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(LOG_FILE, mode="w", encoding="utf-8"),
            logging.StreamHandler(),
        ],
        force=True,
    )


def clean_text(value: object) -> str:
    """Normalize missing-like values to an empty string."""
    if pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"", "na", "nan", "none", "null", "<na>"} else text


def normalize_barcode_core(value: object) -> str:
    """Return the canonical 10x barcode core without suffixes or sample labels."""
    text = clean_text(value).upper()
    if not text:
        return ""
    text = text.split("_")[0]
    return text.split("-")[0]


def normalize_bool_text(value: object) -> bool:
    """Normalize yes/true-like text to bool."""
    return clean_text(value).lower() in {"true", "t", "yes", "y", "1"}


def normalize_string_series(series: pd.Series) -> pd.Series:
    """Normalize a series into plain object strings."""
    return series.map(clean_text).astype(object)


def first_nonblank_series(df: pd.DataFrame, columns: list[str], default: str = "") -> pd.Series:
    """Return the first non-blank value across several columns."""
    out = pd.Series([default] * len(df), index=df.index, dtype=object)
    for column in columns:
        if column not in df.columns:
            continue
        values = normalize_string_series(df[column])
        mask = (out == "") & (values != "")
        out.loc[mask] = values.loc[mask]
    return out


def parse_gene_hit(value: object) -> str:
    """Extract the first gene symbol from Cell Ranger or MiXCR hit text."""
    text = clean_text(value)
    if not text:
        return ""
    first = re.split(r"[;,]", text)[0].strip()
    match = re.match(r"([A-Za-z0-9-]+)\*", first)
    if match:
        return match.group(1)
    match = re.match(r"([A-Za-z0-9-]+)\(", first)
    if match:
        return match.group(1)
    return re.sub(r"\(.*$", "", first).strip()


def ensure_unique_key(df: pd.DataFrame, key_cols: list[str], label: str) -> None:
    """Fail loudly on duplicate join keys."""
    dup_n = int(df.duplicated(key_cols).sum())
    if dup_n:
        raise ValueError(f"{label} has {dup_n} duplicated rows for key {key_cols}")


def load_productive_ab_tcr(ab_dir: Path) -> tuple[pd.DataFrame, dict[str, int]]:
    """Load productive alpha-beta chains and collapse to one row per cell and chain."""
    parts: list[pd.DataFrame] = []
    file_n = 0
    raw_rows = 0
    kept_rows = 0

    for path in sorted(ab_dir.glob("*.TCR.csv")):
        file_n += 1
        frame = pd.read_csv(path, usecols=AB_REQUIRED_COLS, low_memory=False)
        raw_rows += len(frame)
        frame["sample_id"] = normalize_string_series(frame["Samplename"])
        frame["barcode_core"] = frame["barcode"].map(normalize_barcode_core)
        keep_mask = (
            frame["chain"].astype(str).isin(["TRA", "TRB"])
            & frame["is_cell"].map(normalize_bool_text)
            & frame["high_confidence"].map(normalize_bool_text)
            & frame["productive"].map(normalize_bool_text)
            & frame["sample_id"].ne("")
            & frame["barcode_core"].ne("")
            & frame["cdr3"].map(clean_text).ne("")
        )
        frame = frame.loc[keep_mask].copy()
        if frame.empty:
            continue

        kept_rows += len(frame)
        frame["v_gene"] = frame["v_gene"].map(parse_gene_hit)
        frame["d_gene"] = frame["d_gene"].map(parse_gene_hit)
        frame["j_gene"] = frame["j_gene"].map(parse_gene_hit)
        frame["c_gene"] = frame["c_gene"].map(parse_gene_hit)
        frame["cdr3"] = normalize_string_series(frame["cdr3"])
        frame["cdr3_nt"] = normalize_string_series(frame["cdr3_nt"])
        frame["raw_clonotype_id"] = normalize_string_series(frame["raw_clonotype_id"])
        frame["umis"] = pd.to_numeric(frame["umis"], errors="coerce").fillna(0).astype(int)
        frame["reads"] = pd.to_numeric(frame["reads"], errors="coerce").fillna(0).astype(int)

        frame = frame.sort_values(
            ["sample_id", "barcode_core", "chain", "umis", "reads"],
            ascending=[True, True, True, False, False],
        )
        frame = frame.drop_duplicates(["sample_id", "barcode_core", "chain"], keep="first")
        parts.append(frame)

    if not parts:
        raise ValueError(f"No productive alpha-beta TCR rows found under {ab_dir}")

    collapsed = pd.concat(parts, ignore_index=True)
    wide = reshape_chain_rows_wide(
        collapsed,
        field_map={
            "cdr3": "cdr3",
            "v_gene": "v",
            "d_gene": "d",
            "j_gene": "j",
            "cdr3_nt": "cdr3_nt",
            "raw_clonotype_id": "clone_id",
            "umis": "umis",
            "reads": "reads",
            "c_gene": "c_gene",
        },
        chains=["TRA", "TRB"],
    )
    ensure_unique_key(wide, ["sample_id", "barcode_core"], "productive alpha-beta TCR table")
    return wide, {"ab_file_n": file_n, "ab_raw_rows": raw_rows, "ab_kept_rows": kept_rows, "ab_cells": len(wide)}


def load_productive_like_gd_tcr(gd_dir: Path) -> tuple[pd.DataFrame, dict[str, int]]:
    """Load productive-like gamma-delta chains from MiXCR clone tables."""
    parts: list[pd.DataFrame] = []
    file_n = 0
    raw_rows = 0
    kept_rows = 0

    for path in sorted(gd_dir.glob("*_result.clones.tsv")):
        file_n += 1
        sample_id = re.sub(r"-T_result\.clones\.tsv$", "", path.name)
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        raw_rows += len(frame)
        frame["sample_id"] = sample_id
        frame["chain"] = frame["topChains"].map(clean_text)
        frame["barcode_core"] = frame["tagValueCELL"].map(normalize_barcode_core)
        frame["aaSeqCDR3"] = normalize_string_series(frame["aaSeqCDR3"])
        frame["nSeqCDR3"] = normalize_string_series(frame["nSeqCDR3"])
        keep_mask = (
            frame["chain"].isin(["TRG", "TRD"])
            & frame["barcode_core"].ne("")
            & frame["aaSeqCDR3"].ne("")
            & frame["nSeqCDR3"].ne("")
            & ~frame["aaSeqCDR3"].str.contains(r"[*_]", regex=True)
        )
        frame = frame.loc[keep_mask].copy()
        if frame.empty:
            continue

        kept_rows += len(frame)
        frame["v_gene"] = frame["allVHitsWithScore"].map(parse_gene_hit)
        frame["d_gene"] = frame["allDHitsWithScore"].map(parse_gene_hit)
        frame["j_gene"] = frame["allJHitsWithScore"].map(parse_gene_hit)
        frame["c_gene"] = frame["allCHitsWithScore"].map(parse_gene_hit)
        frame["clone_id"] = normalize_string_series(frame["cloneId"])
        frame["uniqueMoleculeCount"] = pd.to_numeric(frame["uniqueMoleculeCount"], errors="coerce").fillna(0).astype(int)
        frame["readCount"] = pd.to_numeric(frame["readCount"], errors="coerce").fillna(0).astype(int)
        frame = frame.sort_values(
            ["sample_id", "barcode_core", "chain", "uniqueMoleculeCount", "readCount"],
            ascending=[True, True, True, False, False],
        )
        frame = frame.drop_duplicates(["sample_id", "barcode_core", "chain"], keep="first")
        parts.append(frame)

    if not parts:
        return pd.DataFrame(columns=["sample_id", "barcode_core"]), {
            "gd_file_n": file_n,
            "gd_raw_rows": raw_rows,
            "gd_kept_rows": 0,
            "gd_cells": 0,
        }

    collapsed = pd.concat(parts, ignore_index=True)
    wide = reshape_chain_rows_wide(
        collapsed,
        field_map={
            "aaSeqCDR3": "cdr3",
            "v_gene": "v",
            "d_gene": "d",
            "j_gene": "j",
            "nSeqCDR3": "cdr3_nt",
            "clone_id": "clone_id",
            "uniqueMoleculeCount": "umis",
            "readCount": "reads",
            "c_gene": "c_gene",
        },
        chains=["TRG", "TRD"],
    )
    ensure_unique_key(wide, ["sample_id", "barcode_core"], "productive-like gamma-delta TCR table")
    return wide, {"gd_file_n": file_n, "gd_raw_rows": raw_rows, "gd_kept_rows": kept_rows, "gd_cells": len(wide)}


def build_harmonized_obs(obs: pd.DataFrame, input_h5ad: Path) -> pd.DataFrame:
    """Add canonical project metadata fields to the source H5AD obs."""
    out = obs.copy()
    out["raw_obs_name"] = out.index.astype(str)
    out["sample_id"] = normalize_string_series(out["Sample_name"])
    out["sampleid"] = out["sample_id"].copy()
    out["library_id"] = out["sample_id"].copy()
    out["cell_id"] = normalize_string_series(out["cellbarcode"])
    out["original_cell_id"] = out["cell_id"].copy()
    out["barcodes"] = out["cell_id"].copy()
    out["barcode"] = out["cell_id"].map(normalize_barcode_core)
    out["barcode_core"] = out["barcode"].copy()
    out["gse_id"] = DATASET_ID
    out["source_gse_id"] = DATASET_ID
    out["project name"] = DATASET_ID
    out["obs_name"] = out["cell_id"].copy()
    out["source_h5ad"] = str(input_h5ad.resolve())
    out["source_root"] = str(input_h5ad.parent.resolve())
    out["sample_type"] = SAMPLE_TYPE
    out["donor_patient"] = normalize_string_series(out["Donor"])
    out["donor_id"] = normalize_string_series(out["Donor"])
    out["age"] = normalize_string_series(out["Age"])
    out["sex"] = ""
    out["tissue"] = first_nonblank_series(out, ["Tissue_merged", "Tissue"]).str.lower()
    out["condition"] = normalize_string_series(out["Stage"])
    out["treatment"] = ""
    out["enrichment_strategy"] = ENRICHMENT
    out["assay_type"] = ASSAY_TYPE
    out["technology_simple"] = RAW_ASSAY
    out["tcr_vdj_flag"] = "yes"
    out["tcr_availability"] = "ab_and_gd_external"
    out["original_cell_annotation"] = first_nonblank_series(
        out,
        ["tnk_detailed_annotation", "tnk_major_cell", "Major_cluster", "Cluster", "Recluster"],
    )
    out["input_population"] = SAMPLE_TYPE
    out["tcr_chain_mode"] = "ab_and_gd"
    return out


def merge_tcr_tables(obs: pd.DataFrame, ab_tcr: pd.DataFrame, gd_tcr: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Merge sample-aware TCR tables into the H5AD obs dataframe."""
    join_cols = ["sample_id", "barcode_core"]
    ensure_unique_key(obs, join_cols, "HRA005041 H5AD obs")
    ab_tcr = ab_tcr.copy()
    gd_tcr = gd_tcr.copy()
    merged_tcr = ab_tcr.merge(gd_tcr, on=join_cols, how="outer")
    ensure_unique_key(merged_tcr, join_cols, "merged TCR table")

    obs_with_tcr = obs.merge(merged_tcr, on=join_cols, how="left", validate="one_to_one", indicator=True)
    unmatched = (
        merged_tcr.merge(obs[join_cols], on=join_cols, how="left", indicator=True)
        .loc[lambda d: d["_merge"] == "left_only", join_cols]
        .copy()
    )
    obs_with_tcr = obs_with_tcr.drop(columns=["_merge"])

    for col in CANONICAL_TCR_COLS:
        if col not in obs_with_tcr.columns:
            obs_with_tcr[col] = ""

    string_cols = [col for col in CANONICAL_TCR_COLS if not col.endswith(("_umis", "_reads"))]
    for col in string_cols:
        obs_with_tcr[col] = normalize_string_series(obs_with_tcr[col])

    int_cols = [col for col in CANONICAL_TCR_COLS if col.endswith(("_umis", "_reads"))]
    for col in int_cols:
        obs_with_tcr[col] = pd.to_numeric(obs_with_tcr[col], errors="coerce").fillna(0).astype(np.int32)

    obs_with_tcr["has_any_ab_tcr"] = (
        obs_with_tcr["TRA_cdr3"].ne("") | obs_with_tcr["TRB_cdr3"].ne("")
    )
    obs_with_tcr["has_any_gd_tcr"] = (
        obs_with_tcr["TRG_cdr3"].ne("") | obs_with_tcr["TRD_cdr3"].ne("")
    )
    obs_with_tcr["has_TRA_TRB_paired"] = (
        obs_with_tcr["TRA_cdr3"].ne("") & obs_with_tcr["TRB_cdr3"].ne("")
    )
    obs_with_tcr["has_TRG_TRD_paired"] = (
        obs_with_tcr["TRG_cdr3"].ne("") & obs_with_tcr["TRD_cdr3"].ne("")
    )
    obs_with_tcr["TCRseq"] = np.where(
        obs_with_tcr["has_any_ab_tcr"] | obs_with_tcr["has_any_gd_tcr"],
        "yes",
        "no",
    )

    return obs_with_tcr, unmatched


def build_sample_summary(obs: pd.DataFrame) -> pd.DataFrame:
    """Build one sample-level TCR coverage summary."""
    rows = []
    grouped = obs.groupby("sample_id", sort=True)
    for sample_id, group in grouped:
        rows.append(
            {
                "sample_id": sample_id,
                "n_cells": int(len(group)),
                "any_ab_tcr_cells": int(group["has_any_ab_tcr"].sum()),
                "paired_ab_cells": int(group["has_TRA_TRB_paired"].sum()),
                "any_gd_tcr_cells": int(group["has_any_gd_tcr"].sum()),
                "paired_gd_cells": int(group["has_TRG_TRD_paired"].sum()),
                "any_tcr_cells": int(group["TCRseq"].eq("yes").sum()),
            }
        )
    sample_df = pd.DataFrame(rows)
    for numer, denom, out_col in [
        ("any_ab_tcr_cells", "n_cells", "any_ab_fraction"),
        ("paired_ab_cells", "n_cells", "paired_ab_fraction"),
        ("any_gd_tcr_cells", "n_cells", "any_gd_fraction"),
        ("paired_gd_cells", "n_cells", "paired_gd_fraction"),
        ("any_tcr_cells", "n_cells", "any_tcr_fraction"),
    ]:
        sample_df[out_col] = np.where(sample_df[denom] > 0, sample_df[numer] / sample_df[denom], 0.0)
    return sample_df


def write_summary(
    *,
    obs: pd.DataFrame,
    sample_df: pd.DataFrame,
    unmatched: pd.DataFrame,
    ab_stats: dict[str, int],
    gd_stats: dict[str, int],
) -> None:
    """Write QC summaries to CSV and markdown."""
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    summary_rows = [
        {"metric": "n_cells", "value": int(len(obs))},
        {"metric": "n_samples", "value": int(obs["sample_id"].nunique())},
        {"metric": "ab_file_n", "value": int(ab_stats["ab_file_n"])},
        {"metric": "ab_raw_rows", "value": int(ab_stats["ab_raw_rows"])},
        {"metric": "ab_kept_productive_rows", "value": int(ab_stats["ab_kept_rows"])},
        {"metric": "ab_unique_cells", "value": int(ab_stats["ab_cells"])},
        {"metric": "gd_file_n", "value": int(gd_stats["gd_file_n"])},
        {"metric": "gd_raw_rows", "value": int(gd_stats["gd_raw_rows"])},
        {"metric": "gd_kept_productive_like_rows", "value": int(gd_stats["gd_kept_rows"])},
        {"metric": "gd_unique_cells", "value": int(gd_stats["gd_cells"])},
        {"metric": "cells_with_any_ab_tcr", "value": int(obs["has_any_ab_tcr"].sum())},
        {"metric": "cells_with_paired_ab_tcr", "value": int(obs["has_TRA_TRB_paired"].sum())},
        {"metric": "cells_with_any_gd_tcr", "value": int(obs["has_any_gd_tcr"].sum())},
        {"metric": "cells_with_paired_gd_tcr", "value": int(obs["has_TRG_TRD_paired"].sum())},
        {"metric": "cells_with_any_tcr", "value": int(obs["TCRseq"].eq("yes").sum())},
        {"metric": "unmatched_tcr_rows", "value": int(len(unmatched))},
    ]
    pd.DataFrame(summary_rows).to_csv(SUMMARY_CSV, index=False)
    sample_df.to_csv(SAMPLE_SUMMARY_CSV, index=False)
    unmatched.groupby("sample_id").size().rename("unmatched_tcr_rows").reset_index().to_csv(UNMATCHED_CSV, index=False)

    lines = [
        "# HRA005041 TCR Intake Summary",
        "",
        f"- Input H5AD: `{INPUT_H5AD}`",
        f"- Output H5AD: `{OUTPUT_H5AD}`",
        f"- Cells: `{len(obs):,}`",
        f"- Samples: `{obs['sample_id'].nunique():,}`",
        f"- Productive abTCR files: `{ab_stats['ab_file_n']}`",
        f"- Productive-like gdTCR files: `{gd_stats['gd_file_n']}`",
        "",
        "## TCR coverage",
        "",
        f"- Cells with any alpha-beta TCR: `{int(obs['has_any_ab_tcr'].sum()):,}`",
        f"- Cells with paired TRA/TRB: `{int(obs['has_TRA_TRB_paired'].sum()):,}`",
        f"- Cells with any gamma-delta TCR: `{int(obs['has_any_gd_tcr'].sum()):,}`",
        f"- Cells with paired TRG/TRD: `{int(obs['has_TRG_TRD_paired'].sum()):,}`",
        f"- Cells with any productive TCR evidence: `{int(obs['TCRseq'].eq('yes').sum()):,}`",
        f"- Unmatched TCR rows after sample-aware join: `{len(unmatched):,}`",
        "",
        "## Join rule",
        "",
        "- TCR was joined by `sample_id + barcode_core`, never barcode alone.",
        "- Alpha-beta rows were filtered to `is_cell=true`, `high_confidence=true`, and `productive=true`.",
        "- Gamma-delta rows were filtered to `topChains in {TRG, TRD}` with non-empty amino-acid and nucleotide CDR3 and no stop/underscore characters in `aaSeqCDR3`.",
        "",
        "## Outputs",
        "",
        f"- `{SUMMARY_CSV}`",
        f"- `{SAMPLE_SUMMARY_CSV}`",
        f"- `{UNMATCHED_CSV}`",
    ]
    QC_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_output_h5ad(adata: ad.AnnData, output_h5ad: Path) -> None:
    """Write the processed H5AD atomically."""
    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    legacy_tmp = output_h5ad.with_suffix(output_h5ad.suffix + ".tmp")
    if legacy_tmp.exists():
        legacy_tmp.unlink()
    tmp_path = output_h5ad.parent / f".{output_h5ad.name}.{os.getpid()}.tmp"
    if tmp_path.exists():
        tmp_path.unlink()
    adata.write_h5ad(tmp_path, convert_strings_to_categoricals=False)
    os.replace(tmp_path, output_h5ad)


def main() -> None:
    """Run HRA005041 productive TCR intake."""
    args = parse_args()
    configure_logging()
    logging.info("Loading source H5AD from %s", args.input_h5ad)
    adata = ad.read_h5ad(args.input_h5ad)
    if "Sample_name" not in adata.obs.columns or "cellbarcode" not in adata.obs.columns:
        raise KeyError("Source H5AD is missing required obs columns `Sample_name` and/or `cellbarcode`.")

    obs = build_harmonized_obs(adata.obs.copy(), args.input_h5ad)
    if obs["cell_id"].duplicated().any():
        dup_n = int(obs["cell_id"].duplicated().sum())
        raise ValueError(f"Source H5AD has {dup_n} duplicated `cell_id` values after normalization.")

    logging.info("Loading productive alpha-beta TCR tables from %s", args.ab_dir)
    ab_tcr, ab_stats = load_productive_ab_tcr(args.ab_dir)
    logging.info(
        "Loaded productive alpha-beta TCR: %s files, %s kept rows, %s unique cells",
        ab_stats["ab_file_n"],
        ab_stats["ab_kept_rows"],
        ab_stats["ab_cells"],
    )
    logging.info("Loading productive-like gamma-delta TCR tables from %s", args.gd_dir)
    gd_tcr, gd_stats = load_productive_like_gd_tcr(args.gd_dir)
    logging.info(
        "Loaded productive-like gamma-delta TCR: %s files, %s kept rows, %s unique cells",
        gd_stats["gd_file_n"],
        gd_stats["gd_kept_rows"],
        gd_stats["gd_cells"],
    )

    logging.info("Merging TCR tables into H5AD obs")
    merged_obs, unmatched = merge_tcr_tables(obs, ab_tcr, gd_tcr)
    logging.info("Merged TCR into obs: %s unmatched TCR rows", len(unmatched))

    merged_obs = merged_obs.copy()
    merged_obs.index = pd.Index(merged_obs["cell_id"].astype(str), name=None)
    adata.obs = merged_obs
    adata.obs.index.name = None
    adata.uns["dataset_accession"] = DATASET_ID
    adata.uns["tcr_join_rule"] = "sample_id + barcode_core"
    adata.uns["tcr_sources"] = {
        "ab_dir": str(args.ab_dir.resolve()),
        "gd_dir": str(args.gd_dir.resolve()),
    }

    sample_df = build_sample_summary(merged_obs)
    write_summary(obs=merged_obs, sample_df=sample_df, unmatched=unmatched, ab_stats=ab_stats, gd_stats=gd_stats)
    logging.info("Wrote join summary tables and markdown QC summary")

    logging.info("Writing processed H5AD to %s", args.output_h5ad)
    write_output_h5ad(adata, args.output_h5ad)
    logging.info("Wrote %s, %s, %s, and %s", args.output_h5ad, SUMMARY_CSV, SAMPLE_SUMMARY_CSV, QC_MD)


if __name__ == "__main__":
    main()
