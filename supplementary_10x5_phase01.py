#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Build and audit the supplementary six-dataset 10x 5' intake lane.

This helper creates per-GSE H5AD files for the approved supplementary datasets,
writes a supplementary registry, runs a separate Phase 0 audit package, and can
run a separate Phase 1 extraction package after explicit QC approval.

The default stop point is Phase 0 to respect the project QC gate.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import re
import shutil
import subprocess
import tarfile
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread

import phase0_dataset_audit as phase0
import phase1_extract_tnk_candidates as phase1


PROJECT_ROOT = Path(__file__).resolve().parent
DOWNLOADS_DIR = PROJECT_ROOT / "downloads"
PER_GSE_OUT_DIR = DOWNLOADS_DIR / "per_gse_h5ad_with_metadata"
R_HELPER = PROJECT_ROOT / "supplementary_export_rds_payloads.R"
PYTHON_BIN = Path("/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python")

OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
SUPP_TABLE_DIR = OUTPUT_ROOT / "tables" / "supplementary_10x5"
SUPP_LOG_DIR = OUTPUT_ROOT / "logs" / "supplementary_10x5"
SUPP_FIG_DIR = OUTPUT_ROOT / "figures" / "supplementary_10x5"
SUPP_TMP_DIR = OUTPUT_ROOT / "_tmp_supplementary_payloads"
SUPP_PHASE1_TMP_DIR = OUTPUT_ROOT / "_tmp_phase1_candidates_supp"

SUPP_REGISTRY_CSV = SUPP_TABLE_DIR / "h5ad_supplementary_10x5.csv"
SUPP_BUILD_SUMMARY_CSV = SUPP_TABLE_DIR / "supplementary_h5ad_build_summary.csv"
SUPP_BASE_METADATA_CSV = SUPP_TABLE_DIR / "supplementary_metadata_all_cells.csv.gz"
SUPP_METADATA_CANDIDATE_CSV = PROJECT_ROOT / "analysis_26GSE_V4" / "outputs" / "harmonized_metadata_supp.csv"
SUPP_CANDIDATES_H5AD = OUTPUT_ROOT / "TNK_candidates_supp.h5ad"

SUPP_GSES = [
    "GSE179994",
    "GSE235863",
    "GSE240865",
    "GSE287301",
    "GSE234069",
    "GSE287541",
]

BASE_METADATA_COLS = [
    "gse_id",
    "cell_id",
    "source_h5ad",
    "source_root",
    "library_id",
    "sample_id",
    "donor_id",
    "age",
    "sex",
    "tissue",
    "condition",
    "treatment",
    "enrichment_strategy",
    "assay_type",
    "tcr_availability",
    "original_cell_annotation",
    "sample_type",
    "donor_patient",
    "technology_simple",
    "tcr_vdj_flag",
    "TCRseq",
    "TRA_cdr3",
    "TRA_v",
    "TRA_d",
    "TRA_j",
    "TRA_cdr3_nt",
    "TRA_clone_id",
    "TRA_umis",
    "TRA_reads",
    "TRB_cdr3",
    "TRB_v",
    "TRB_d",
    "TRB_j",
    "TRB_cdr3_nt",
    "TRB_clone_id",
    "TRB_umis",
    "TRB_reads",
    "barcode",
]

TCR_COLS = [
    "TRA_cdr3",
    "TRA_v",
    "TRA_d",
    "TRA_j",
    "TRA_cdr3_nt",
    "TRA_clone_id",
    "TRA_umis",
    "TRA_reads",
    "TRB_cdr3",
    "TRB_v",
    "TRB_d",
    "TRB_j",
    "TRB_cdr3_nt",
    "TRB_clone_id",
    "TRB_umis",
    "TRB_reads",
]

HARMONIZED_HEADER_SOURCE = PROJECT_ROOT / "analysis_26GSE_V4" / "outputs" / "harmonized_metadata_v4.csv"


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gse",
        nargs="+",
        choices=SUPP_GSES,
        help="Optional subset of supplementary GSEs to build.",
    )
    parser.add_argument(
        "--stop-after",
        choices=["phase0", "phase1"],
        default="phase0",
        help="Run through supplementary Phase 0 only, or continue to supplementary Phase 1 after approval.",
    )
    parser.add_argument(
        "--force-rebuild",
        action="store_true",
        help="Rebuild per-GSE supplementary H5AD files even if they already exist.",
    )
    parser.add_argument(
        "--build-only",
        action="store_true",
        help="Only build per-GSE H5AD files and the supplementary registry; do not run Phase 0 or Phase 1.",
    )
    return parser.parse_args()


def log(message: str) -> None:
    """Print a readable progress message."""
    print(message, flush=True)


def ensure_dirs() -> None:
    """Create required supplementary directories."""
    for path in [PER_GSE_OUT_DIR, SUPP_TABLE_DIR, SUPP_LOG_DIR, SUPP_FIG_DIR, SUPP_TMP_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def clean_text(value: Any) -> str:
    """Normalize missing-like values to an empty string."""
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "na", "none", "null", "<na>"}:
        return ""
    return text


def normalize_yes_no(value: Any) -> str:
    """Normalize truthy metadata to yes/no text."""
    text = clean_text(value).lower()
    if text in {"yes", "y", "true", "1"}:
        return "yes"
    if text in {"no", "n", "false", "0"}:
        return "no"
    return ""


def normalize_barcode_core(value: Any) -> str:
    """Return the 10x sequence core before the first dash."""
    text = clean_text(value).upper()
    if not text:
        return ""
    match = re.search(r"([ACGTN]+)-\d+$", text)
    if match:
        return match.group(1)
    match = re.search(r"([ACGTN]+)-\d+", text)
    if match:
        return match.group(1)
    return text.split("-")[0]


def extract_10x_barcode(value: Any) -> str:
    """Extract a full 10x barcode with lane suffix if present."""
    text = clean_text(value)
    match = re.search(r"([ACGTN]+-\d+)", text, flags=re.IGNORECASE)
    return match.group(1).upper() if match else text


def normalize_string_series(series: pd.Series) -> pd.Series:
    """Return a plain object string series with blank missing values."""
    return series.map(clean_text).astype(object)


def obs_series(obs: pd.DataFrame, column: str, default: str = "") -> pd.Series:
    """Return one obs column or a same-length blank series."""
    if column in obs.columns:
        return normalize_string_series(obs[column])
    return pd.Series([default] * len(obs), index=obs.index, dtype=object)


def first_nonblank_series(df: pd.DataFrame, columns: list[str]) -> pd.Series:
    """Return the first non-blank value across several columns."""
    result = pd.Series([""] * len(df), index=df.index, dtype=object)
    for column in columns:
        if column not in df.columns:
            continue
        values = normalize_string_series(df[column])
        mask = (result == "") & (values != "")
        result.loc[mask] = values.loc[mask]
    return result


def normalize_phase1_sampleid(base_sample: pd.Series, phase1_sample: pd.Series) -> pd.Series:
    """Prefer harmonized sample_id and fall back to the Phase 1 sample label."""
    out = normalize_string_series(base_sample)
    fallback = normalize_string_series(phase1_sample)
    mask = (out == "") & (fallback != "")
    out.loc[mask] = fallback.loc[mask]
    return out


def parse_geo_series_matrix_samples(gse_id: str) -> pd.DataFrame:
    """Parse GEO series-matrix sample metadata into one dataframe."""
    path = DOWNLOADS_DIR / gse_id / "matrix" / f"{gse_id}_series_matrix.txt.gz"
    grouped: dict[str, list[list[str]]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            if raw_line.startswith("!sample_table_begin") or raw_line.startswith('"ID_REF"'):
                break
            if not raw_line.startswith("!Sample_"):
                continue
            parts = next(csv.reader([raw_line.rstrip("\n")], delimiter="\t", quotechar='"'))
            key = parts[0].replace("!Sample_", "").lower()
            grouped.setdefault(key, []).append(parts[1:])

    sample_n = max(len(values[0]) for values in grouped.values())
    rows = [dict() for _ in range(sample_n)]

    for key, repeats in grouped.items():
        if len(repeats) == 1:
            values = repeats[0]
            for idx, value in enumerate(values):
                rows[idx][key] = value
        else:
            for repeat_idx, values in enumerate(repeats, start=1):
                column = f"{key}_{repeat_idx}"
                for idx, value in enumerate(values):
                    rows[idx][column] = value

    sample_df = pd.DataFrame(rows)
    if "geo_accession" not in sample_df.columns:
        raise ValueError(f"{gse_id} series matrix did not provide !Sample_geo_accession")

    for column in [c for c in sample_df.columns if c.startswith("characteristics_ch1")]:
        values = sample_df[column].map(clean_text)
        for idx, value in values.items():
            if ":" not in value:
                continue
            key, parsed_value = value.split(":", 1)
            parsed_key = re.sub(r"[^a-z0-9]+", "_", key.strip().lower()).strip("_")
            if not parsed_key:
                continue
            existing = clean_text(sample_df.at[idx, parsed_key]) if parsed_key in sample_df.columns else ""
            if not existing:
                sample_df.at[idx, parsed_key] = parsed_value.strip()

    sample_df["gse_id"] = gse_id
    return sample_df


def write_h5ad_atomic(adata: ad.AnnData, output_path: Path) -> None:
    """Write a dataset H5AD atomically."""
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    if tmp_path.exists():
        tmp_path.unlink()
    adata.write_h5ad(tmp_path)
    tmp_path.replace(output_path)


def sanitize_dataframe_for_h5ad(frame: pd.DataFrame) -> pd.DataFrame:
    """Convert mixed dataframe columns to H5AD-safe dtypes."""
    clean = frame.copy()
    clean.index = clean.index.map(str)

    for column in clean.columns:
        series = clean[column]
        if isinstance(series.dtype, pd.CategoricalDtype):
            series = series.astype(object)

        if pd.api.types.is_object_dtype(series) or pd.api.types.is_string_dtype(series):
            non_missing = series[series.notna()]
            bool_like = non_missing.map(lambda value: isinstance(value, (bool, np.bool_))).all() if len(non_missing) else False
            if bool_like:
                clean[column] = series.fillna(False).astype(bool)
            else:
                clean[column] = series.map(clean_text).astype(object)
    return clean


def run_subprocess(command: list[str]) -> None:
    """Run a subprocess and fail loudly on non-zero exit."""
    log(f"Running: {' '.join(command)}")
    subprocess.run(command, check=True, cwd=PROJECT_ROOT)


def maybe_backup_existing(output_path: Path) -> None:
    """Create a lightweight backup before overwriting an existing per-GSE H5AD."""
    if not output_path.exists():
        return
    backup_path = output_path.with_name(output_path.name + ".bak_before_supplementary_rebuild")
    if backup_path.exists():
        return
    shutil.copy2(output_path, backup_path)


def load_payload(prefix: Path) -> tuple[sp.csr_matrix, pd.DataFrame, pd.DataFrame]:
    """Load a sparse payload exported by the R helper."""
    matrix = mmread(str(prefix) + "_counts.mtx").tocsr().transpose().tocsr()
    genes = pd.read_csv(str(prefix) + "_genes.tsv.gz", sep="\t", header=None, names=["gene_id", "gene_name"])
    metadata_path = Path(str(prefix) + "_metadata.csv.gz")
    metadata = pd.read_csv(metadata_path) if metadata_path.exists() else pd.DataFrame()
    return matrix, genes, metadata


def read_triplet_with_feature_filter(prefix: Path) -> ad.AnnData:
    """Read one 10x triplet and keep only gene-expression rows."""
    def read_table(path: Path) -> pd.DataFrame:
        try:
            return pd.read_csv(path, sep="\t", header=None)
        except Exception:
            return pd.read_csv(path, sep="\t", header=None, compression=None)

    barcodes = read_table(Path(str(prefix) + "_barcodes.tsv.gz"))
    features = read_table(Path(str(prefix) + "_features.tsv.gz"))
    matrix = mmread(str(prefix) + "_matrix.mtx.gz").tocsr().transpose().tocsr()

    if features.shape[1] >= 3:
        keep = features.iloc[:, 2].astype(str) == "Gene Expression"
        features = features.loc[keep].reset_index(drop=True)
        matrix = matrix[:, keep.to_numpy()]

    obs = pd.DataFrame(index=barcodes.iloc[:, 0].astype(str))
    obs["barcode"] = obs.index.astype(str)
    var = pd.DataFrame(index=features.iloc[:, 1].astype(str))
    var.index.name = None
    var["gene_ids"] = features.iloc[:, 0].astype(str).to_numpy()
    if features.shape[1] >= 3:
        var["feature_types"] = features.iloc[:, 2].astype(str).to_numpy()

    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    adata.var_names_make_unique()
    return adata


def read_tarred_feature_bc_matrix(path: Path) -> tuple[sp.csr_matrix, pd.DataFrame, pd.DataFrame]:
    """Read a tarred 10x feature-barcode matrix and keep only gene-expression rows."""
    with tarfile.open(path, "r:gz") as tar:
        members = {member.name: member for member in tar.getmembers() if member.isfile()}
        matrix_name = next(name for name in members if name.endswith("matrix.mtx.gz"))
        barcodes_name = next(name for name in members if name.endswith("barcodes.tsv.gz"))
        features_name = next(name for name in members if name.endswith("features.tsv.gz"))

        with tar.extractfile(members[matrix_name]) as handle:
            matrix = mmread(gzip.open(handle, "rt")).tocsr().transpose().tocsr()
        with tar.extractfile(members[barcodes_name]) as handle:
            barcodes = pd.read_csv(gzip.open(handle, "rt"), sep="\t", header=None)
        with tar.extractfile(members[features_name]) as handle:
            features = pd.read_csv(gzip.open(handle, "rt"), sep="\t", header=None)

    if features.shape[1] >= 3:
        keep = features.iloc[:, 2].astype(str) == "Gene Expression"
        features = features.loc[keep].reset_index(drop=True)
        matrix = matrix[:, keep.to_numpy()]
    return matrix, barcodes, features


def find_tar_member(tar: tarfile.TarFile, suffix: str) -> tarfile.TarInfo:
    """Return the preferred real data member, skipping AppleDouble sidecars."""
    candidates = [
        member
        for member in tar.getmembers()
        if member.isfile()
        and member.name.endswith(suffix)
        and not Path(member.name).name.startswith("._")
    ]
    if not candidates:
        raise FileNotFoundError(f"No tar member matched suffix `{suffix}`")
    return candidates[0]


def parse_cdr3_pairs(value: Any) -> dict[str, str]:
    """Parse a cdr3s_aa string like 'TRB:...;TRA:...' into first-chain values."""
    out = {"TRA_cdr3": "", "TRB_cdr3": ""}
    text = clean_text(value)
    if not text:
        return out
    for token in text.split(";"):
        token = token.strip()
        if ":" not in token:
            continue
        chain, cdr3 = token.split(":", 1)
        chain = chain.strip().upper()
        cdr3 = cdr3.strip()
        if chain == "TRA" and not out["TRA_cdr3"]:
            out["TRA_cdr3"] = cdr3
        if chain == "TRB" and not out["TRB_cdr3"]:
            out["TRB_cdr3"] = cdr3
    return out


def empty_tcr_frame(index: pd.Index) -> pd.DataFrame:
    """Return an empty TCR dataframe on one index."""
    frame = pd.DataFrame(index=index)
    for column in TCR_COLS:
        frame[column] = ""
    frame["TCRseq"] = "no"
    return frame


def first_nonempty_value(series: pd.Series) -> str:
    """Return the first non-empty normalized value from a series."""
    values = normalize_string_series(series)
    non_empty = values[values != ""]
    return non_empty.iloc[0] if len(non_empty) else ""


def aggregate_tcr_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse duplicate TCR rows by sample plus barcode core."""
    if df.empty:
        return df
    keep_cols = ["sample_id", "barcode_core"] + TCR_COLS
    work = df[keep_cols].copy()
    numeric_cols = ["TRA_umis", "TRA_reads", "TRB_umis", "TRB_reads"]
    string_cols = [column for column in TCR_COLS if column not in numeric_cols]

    work["sample_id"] = normalize_string_series(work["sample_id"])
    work["barcode_core"] = normalize_string_series(work["barcode_core"])
    for column in string_cols:
        work[column] = normalize_string_series(work[column])
    for column in numeric_cols:
        work[column] = pd.to_numeric(work[column], errors="coerce").fillna(0).astype(int)

    work = work[(work["sample_id"] != "") & (work["barcode_core"] != "")].copy()
    if work.empty:
        return pd.DataFrame(columns=["sample_id", "barcode_core"] + TCR_COLS)

    agg_spec: dict[str, Any] = {column: "max" for column in numeric_cols}
    for column in string_cols:
        agg_spec[column] = first_nonempty_value
    return work.groupby(["sample_id", "barcode_core"], as_index=False, sort=False).agg(agg_spec)


def merge_tcr_into_obs(obs: pd.DataFrame, tcr_df: pd.DataFrame) -> pd.DataFrame:
    """Merge standardized TCR columns into one obs dataframe."""
    merged = obs.copy()
    merged["barcode_core"] = merged["barcode"].map(normalize_barcode_core)
    merged["sample_id"] = normalize_string_series(merged["sample_id"])

    if tcr_df.empty:
        merged = merged.join(empty_tcr_frame(merged.index))
        return merged.drop(columns=["barcode_core"])

    tcr_use = aggregate_tcr_rows(tcr_df.copy())
    joined = (
        merged.reset_index(names="obs_index")
        .merge(
            tcr_use,
            on=["sample_id", "barcode_core"],
            how="left",
            suffixes=("", "_tcr"),
        )
        .set_index("obs_index")
    )

    numeric_cols = ["TRA_umis", "TRA_reads", "TRB_umis", "TRB_reads"]
    string_cols = [column for column in TCR_COLS if column not in numeric_cols]
    for column in TCR_COLS:
        if column not in joined.columns:
            joined[column] = 0 if column in numeric_cols else ""
    for column in string_cols:
        joined[column] = normalize_string_series(joined[column]).astype(object)
    for column in numeric_cols:
        joined[column] = pd.to_numeric(joined[column], errors="coerce").fillna(0).astype(int)

    tcr_mask = np.zeros(len(joined), dtype=bool)
    for column in string_cols:
        values = joined[column].astype(str)
        tcr_mask |= (values != "") & (values != "0")
    for column in numeric_cols:
        tcr_mask |= joined[column].to_numpy() != 0
    joined["TCRseq"] = np.where(tcr_mask, "yes", "no")
    return joined.drop(columns=["barcode_core"])


def standardize_tcr_from_alpha_beta(
    df: pd.DataFrame,
    sample_column: str,
    barcode_column: str,
    tra_prefix: str,
    trb_prefix: str,
    clone_column: str | None = None,
) -> pd.DataFrame:
    """Standardize wide alpha/beta TCR tables to the shared schema."""
    out = pd.DataFrame()
    out["sample_id"] = df[sample_column].map(clean_text)
    out["barcode_core"] = df[barcode_column].map(extract_10x_barcode).map(normalize_barcode_core)

    def pick_series(base: str, first_suffix: str, second_suffix: str | None = None) -> pd.Series:
        candidates = [f"{base}{first_suffix}"]
        if second_suffix:
            candidates.append(f"{base}{second_suffix}")
        result = pd.Series([""] * len(df), index=df.index, dtype=object)
        for column in candidates:
            if column not in df.columns:
                continue
            values = df[column].map(clean_text)
            mask = (result == "") & (values != "")
            result.loc[mask] = values.loc[mask]
        return result

    out["TRA_cdr3"] = pick_series("", f"{tra_prefix}cdr3", None)
    out["TRA_v"] = pick_series("", f"{tra_prefix}v_gene", None)
    out["TRA_d"] = pick_series("", f"{tra_prefix}d_gene", None)
    out["TRA_j"] = pick_series("", f"{tra_prefix}j_gene", None)
    out["TRA_cdr3_nt"] = pick_series("", f"{tra_prefix}cdr3_nt", None)
    out["TRA_clone_id"] = df[clone_column].map(clean_text) if clone_column and clone_column in df.columns else ""
    out["TRA_umis"] = pd.to_numeric(df.get(f"{tra_prefix}expr", df.get(f"{tra_prefix}umis", 0)), errors="coerce").fillna(0).astype(int)
    out["TRA_reads"] = pd.to_numeric(df.get(f"{tra_prefix}reads", df.get(f"{tra_prefix}expr", 0)), errors="coerce").fillna(0).astype(int)

    out["TRB_cdr3"] = pick_series("", f"{trb_prefix}cdr3", None)
    out["TRB_v"] = pick_series("", f"{trb_prefix}v_gene", None)
    out["TRB_d"] = pick_series("", f"{trb_prefix}d_gene", None)
    out["TRB_j"] = pick_series("", f"{trb_prefix}j_gene", None)
    out["TRB_cdr3_nt"] = pick_series("", f"{trb_prefix}cdr3_nt", None)
    out["TRB_clone_id"] = df[clone_column].map(clean_text) if clone_column and clone_column in df.columns else ""
    out["TRB_umis"] = pd.to_numeric(df.get(f"{trb_prefix}expr", df.get(f"{trb_prefix}umis", 0)), errors="coerce").fillna(0).astype(int)
    out["TRB_reads"] = pd.to_numeric(df.get(f"{trb_prefix}reads", df.get(f"{trb_prefix}expr", 0)), errors="coerce").fillna(0).astype(int)

    return out


def standardize_tcr_from_contig(df: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    """Standardize 10x filtered_contig_annotations tables."""
    if "barcode" not in df.columns or "chain" not in df.columns:
        return pd.DataFrame(columns=["sample_id", "barcode_core"] + TCR_COLS)

    work = pd.DataFrame(
        {
            "sample_id": sample_id,
            "barcode_core": df["barcode"].map(extract_10x_barcode).map(normalize_barcode_core),
            "chain": df["chain"].astype(str).str.upper(),
            "v_gene": df.get("v_gene", "").map(clean_text),
            "d_gene": df.get("d_gene", "").map(clean_text),
            "j_gene": df.get("j_gene", "").map(clean_text),
            "cdr3": df.get("cdr3", "").map(clean_text),
            "cdr3_nt": df.get("cdr3_nt", "").map(clean_text),
            "clone_id": df.get("raw_clonotype_id", "").map(clean_text),
            "umis": pd.to_numeric(df.get("umis", 0), errors="coerce").fillna(0).astype(int),
            "reads": pd.to_numeric(df.get("reads", 0), errors="coerce").fillna(0).astype(int),
        }
    )
    work = work[work["barcode_core"] != ""].copy()
    work = work[work["chain"].isin(["TRA", "TRB"])].copy()
    if work.empty:
        return pd.DataFrame(columns=["sample_id", "barcode_core"] + TCR_COLS)

    agg = (
        work.groupby(["sample_id", "barcode_core", "chain"], as_index=False, sort=False)
        .agg(
            {
                "cdr3": first_nonempty_value,
                "v_gene": first_nonempty_value,
                "d_gene": first_nonempty_value,
                "j_gene": first_nonempty_value,
                "cdr3_nt": first_nonempty_value,
                "clone_id": first_nonempty_value,
                "umis": "max",
                "reads": "max",
            }
        )
    )

    base = work[["sample_id", "barcode_core"]].drop_duplicates().copy()
    outputs = [base]
    rename_map = {
        "cdr3": "cdr3",
        "v_gene": "v",
        "d_gene": "d",
        "j_gene": "j",
        "cdr3_nt": "cdr3_nt",
        "clone_id": "clone_id",
        "umis": "umis",
        "reads": "reads",
    }
    for chain in ["TRA", "TRB"]:
        subset = agg.loc[agg["chain"] == chain, ["sample_id", "barcode_core", *rename_map.keys()]].copy()
        subset = subset.rename(columns={column: f"{chain}_{suffix}" for column, suffix in rename_map.items()})
        outputs.append(subset)

    out = outputs[0]
    for subset in outputs[1:]:
        out = out.merge(subset, on=["sample_id", "barcode_core"], how="left")

    numeric_cols = ["TRA_umis", "TRA_reads", "TRB_umis", "TRB_reads"]
    string_cols = [column for column in TCR_COLS if column not in numeric_cols]
    for column in string_cols:
        if column not in out.columns:
            out[column] = ""
        out[column] = normalize_string_series(out[column]).astype(object)
    for column in numeric_cols:
        if column not in out.columns:
            out[column] = 0
        out[column] = pd.to_numeric(out[column], errors="coerce").fillna(0).astype(int)
    return out[["sample_id", "barcode_core"] + TCR_COLS]


def parse_235863_tcr(path: Path) -> pd.DataFrame:
    """Parse the HCC paired-TCR CSV schema into the shared layout."""
    df = pd.read_csv(path, low_memory=False)
    out = pd.DataFrame()
    out["sample_id"] = df["sample"].map(clean_text)
    out["barcode_core"] = df["barcode"].map(extract_10x_barcode).map(normalize_barcode_core)
    out["TRA_cdr3"] = df.get("cdr3.A", "").map(clean_text)
    out["TRA_v"] = df.get("v_gene.A", "").map(clean_text)
    out["TRA_d"] = ""
    out["TRA_j"] = df.get("j_gene.A", "").map(clean_text)
    out["TRA_cdr3_nt"] = df.get("cdr3_nt.A", "").map(clean_text)
    out["TRA_clone_id"] = df.get("clone.id", "").map(clean_text)
    out["TRA_umis"] = 0
    out["TRA_reads"] = 0
    out["TRB_cdr3"] = df.get("cdr3.B", "").map(clean_text)
    out["TRB_v"] = df.get("v_gene.B", "").map(clean_text)
    out["TRB_d"] = ""
    out["TRB_j"] = df.get("j_gene.B", "").map(clean_text)
    out["TRB_cdr3_nt"] = df.get("cdr3_nt.B", "").map(clean_text)
    out["TRB_clone_id"] = df.get("clone.id", "").map(clean_text)
    out["TRB_umis"] = 0
    out["TRB_reads"] = 0
    return out


def parse_179994_tcr(path: Path) -> pd.DataFrame:
    """Parse the lung-cancer alpha/beta TCR table into the shared layout."""
    df = pd.read_csv(path, sep="\t", low_memory=False)
    out = pd.DataFrame()
    out["sample_id"] = df["sample"].map(clean_text)
    out["barcode_core"] = df["CellName"].map(extract_10x_barcode).map(normalize_barcode_core)
    out["TRA_cdr3"] = first_nonblank_series(df, ["CDR3(Alpha1)", "CDR3(Alpha2)"])
    out["TRA_v"] = first_nonblank_series(df, ["V_gene(Alpha1)", "V_gene(Alpha2)"])
    out["TRA_d"] = ""
    out["TRA_j"] = first_nonblank_series(df, ["J_gene(Alpha1)", "J_gene(Alpha2)"])
    out["TRA_cdr3_nt"] = first_nonblank_series(df, ["CDR3_nt(Alpha1)", "CDR3_nt(Alpha2)"])
    out["TRA_clone_id"] = df.get("clone.id", "").map(clean_text)
    out["TRA_umis"] = pd.to_numeric(first_nonblank_series(df, ["nUMI(Alpha1)", "nUMI(Alpha2)"]), errors="coerce").fillna(0).astype(int)
    out["TRA_reads"] = pd.to_numeric(first_nonblank_series(df, ["nRead(Alpha1)", "nRead(Alpha2)"]), errors="coerce").fillna(0).astype(int)
    out["TRB_cdr3"] = first_nonblank_series(df, ["CDR3(Beta1)", "CDR3(Beta2)"])
    out["TRB_v"] = first_nonblank_series(df, ["V_gene(Beta1)", "V_gene(Beta2)"])
    out["TRB_d"] = ""
    out["TRB_j"] = first_nonblank_series(df, ["J_gene(Beta1)", "J_gene(Beta2)"])
    out["TRB_cdr3_nt"] = first_nonblank_series(df, ["CDR3_nt(Beta1)", "CDR3_nt(Beta2)"])
    out["TRB_clone_id"] = df.get("clone.id", "").map(clean_text)
    out["TRB_umis"] = pd.to_numeric(first_nonblank_series(df, ["nUMI(Beta1)", "nUMI(Beta2)"]), errors="coerce").fillna(0).astype(int)
    out["TRB_reads"] = pd.to_numeric(first_nonblank_series(df, ["nRead(Beta1)", "nRead(Beta2)"]), errors="coerce").fillna(0).astype(int)
    return out


def set_cell_ids(adata: ad.AnnData, library_series: pd.Series) -> None:
    """Set unique obs_names from library/sample plus barcode."""
    barcode = normalize_string_series(adata.obs["barcode"])
    library = normalize_string_series(library_series)
    adata.obs_names = pd.Index(
        [
            f"{lib}:{bc}" if lib else bc
            for lib, bc in zip(library.tolist(), barcode.tolist())
        ]
    )


def finalize_per_gse_adata(
    adata: ad.AnnData,
    gse_id: str,
    output_path: Path,
    source_note: str,
    matrix_file: str,
    platform_id: str,
    technology_simple: str,
) -> dict[str, Any]:
    """Finalize one per-GSE H5AD and write it to the shared output directory."""
    adata.obs.index.name = None
    adata.var.index.name = None
    adata.var_names = pd.Index(adata.var_names.astype(str))
    adata.var_names_make_unique()
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(np.asarray(adata.X))
    elif not sp.isspmatrix_csr(adata.X):
        adata.X = adata.X.tocsr()

    obs = adata.obs.copy()
    obs["gse"] = gse_id
    obs["GSE_ID"] = gse_id
    obs["project_id"] = gse_id
    obs["dataset"] = "0"
    obs["source"] = source_note
    obs["matrix_file"] = matrix_file
    obs["platform_id"] = platform_id
    obs["technology_simple"] = technology_simple
    tcrseq_series = obs_series(obs, "TCRseq")
    obs["sample_type"] = obs_series(obs, "sample_type")
    obs["donor_patient"] = obs_series(obs, "donor_patient")
    obs["tcr_vdj_flag"] = np.where(tcrseq_series.str.lower() == "yes", "yes", obs_series(obs, "tcr_vdj_flag"))
    obs["tcr_availability"] = np.where(tcrseq_series.str.lower() == "yes", "yes", obs_series(obs, "tcr_availability"))

    for column in ["sample_id", "library_id", "barcode"]:
        if column not in obs.columns:
            raise ValueError(f"{gse_id} is missing required obs column `{column}` before final write")

    set_cell_ids(adata, normalize_string_series(obs["library_id"]).replace("", np.nan).fillna(normalize_string_series(obs["sample_id"])))
    obs.index = adata.obs_names.copy()
    adata.obs = sanitize_dataframe_for_h5ad(obs)
    adata.var = sanitize_dataframe_for_h5ad(adata.var)
    adata.uns = {}
    adata.obsm = {}
    adata.varm = {}
    adata.obsp = {}
    adata.varp = {}
    adata.raw = None

    maybe_backup_existing(output_path)
    write_h5ad_atomic(adata, output_path)

    obs = adata.obs.copy()
    summary = {
        "gse_id": gse_id,
        "h5ad_path": str(output_path),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "samples": int(normalize_string_series(obs["sample_id"]).replace("", np.nan).nunique(dropna=True)),
        "libraries": int(normalize_string_series(obs["library_id"]).replace("", np.nan).nunique(dropna=True)),
        "tcr_yes_cells": int((normalize_string_series(obs["TCRseq"]).str.lower() == "yes").sum()),
    }
    return summary


def series_patient_to_short(patient_text: str) -> str:
    """Convert series-matrix patient labels like P001 to P1."""
    text = clean_text(patient_text).upper()
    match = re.fullmatch(r"P0*([0-9]+)", text)
    if match:
        return f"P{int(match.group(1))}"
    return text


def build_179994_gsm_map(series_df: pd.DataFrame) -> dict[str, str]:
    """Map GSE179994 sample IDs to GSM accessions."""
    mapping: dict[str, str] = {}
    for _, row in series_df.iterrows():
        patient = series_patient_to_short(row.get("patient_id", ""))
        title = clean_text(row.get("title", ""))
        condition = clean_text(row.get("condition", "")).lower()
        if not patient:
            continue
        if "post" in condition or ".post." in title.lower():
            match = re.search(r"\.post\.(\d+)$", title.lower())
            suffix = int(match.group(1)) if match else 1
            sample_id = f"{patient}.post.{suffix}"
        else:
            sample_id = f"{patient}.pre"
        mapping[sample_id] = clean_text(row.get("geo_accession", ""))
    return mapping


def build_gse179994(temp_dir: Path) -> tuple[Path, dict[str, Any]]:
    """Build the supplementary H5AD for GSE179994."""
    gse_id = "GSE179994"
    log(f"Building {gse_id}")
    prefix = temp_dir / gse_id / gse_id
    prefix.parent.mkdir(parents=True, exist_ok=True)
    run_subprocess(
        [
            "Rscript",
            str(R_HELPER),
            "--mode",
            "dgcmatrix",
            "--input",
            str(DOWNLOADS_DIR / gse_id / "suppl" / "GSE179994_all.Tcell.rawCounts.rds.gz"),
            "--out-prefix",
            str(prefix),
        ]
    )
    matrix, genes, _ = load_payload(prefix)
    metadata = pd.read_csv(DOWNLOADS_DIR / gse_id / "suppl" / "GSE179994_Tcell.metadata.tsv.gz", sep="\t", low_memory=False)
    metadata = metadata.drop_duplicates("cellid", keep="first").set_index("cellid")

    cell_ids = pd.Index(pd.read_csv(str(prefix) + "_barcodes.tsv.gz", sep="\t", header=None)[0].astype(str))
    keep_cells = cell_ids.isin(metadata.index)
    matrix = matrix[np.asarray(keep_cells), :]
    cell_ids = cell_ids[keep_cells]
    obs = metadata.loc[cell_ids].copy()
    obs.index = cell_ids
    obs["barcode"] = obs.index.to_series().map(extract_10x_barcode).to_numpy()

    series_df = parse_geo_series_matrix_samples(gse_id)
    sample_to_gsm = build_179994_gsm_map(series_df)

    obs["sample_id"] = normalize_string_series(obs["sample"])
    obs["sample"] = obs["sample_id"]
    obs["library_id"] = obs["sample_id"].map(lambda value: sample_to_gsm.get(value, value))
    obs["donor_id"] = normalize_string_series(obs["patient"])
    obs["donor_patient"] = obs["donor_id"]
    obs["sample_type"] = "tumor biopsy"
    obs["tissue"] = "tumor biopsy"
    obs["condition"] = np.where(obs["sample_id"].str.contains(".post.", regex=False), "On-treatment", "Pre-treatment")
    obs["treatment"] = "PD-1 blockade plus chemotherapy"
    obs["age"] = ""
    obs["sex"] = ""
    obs["enrichment_strategy"] = ""
    obs["assay_type"] = "10x Chromium Single cell 5' + VDJ"
    obs["original_cell_annotation"] = first_nonblank_series(obs, ["celltype", "cluster"])

    tcr = parse_179994_tcr(DOWNLOADS_DIR / gse_id / "suppl" / "GSE179994_all.scTCR.tsv.gz")
    obs = merge_tcr_into_obs(obs, tcr)

    var = pd.DataFrame(index=pd.Index(genes["gene_name"].astype(str)))
    var["gene_ids"] = genes["gene_id"].astype(str).to_numpy()
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    output_path = PER_GSE_OUT_DIR / f"{gse_id}_with_tcr.h5ad"
    summary = finalize_per_gse_adata(
        adata=adata,
        gse_id=gse_id,
        output_path=output_path,
        source_note="GSE179994_all.Tcell.rawCounts.rds.gz + GSE179994_Tcell.metadata.tsv.gz",
        matrix_file=str(DOWNLOADS_DIR / gse_id / "matrix" / f"{gse_id}_series_matrix.txt.gz"),
        platform_id="GPL24676",
        technology_simple="10x 5'",
    )
    return output_path, summary


def load_gzipped_h5ad(path: Path) -> ad.AnnData:
    """Load a gzip-compressed H5AD by expanding it to a temporary file first."""
    with tempfile.NamedTemporaryFile(suffix=".h5ad", dir=SUPP_TMP_DIR, delete=False) as handle:
        temp_path = Path(handle.name)
    try:
        with gzip.open(path, "rb") as src, temp_path.open("wb") as dst:
            shutil.copyfileobj(src, dst)
        adata = ad.read_h5ad(temp_path)
    finally:
        if temp_path.exists():
            temp_path.unlink()
    return adata


def build_gse235863(_temp_dir: Path | None = None) -> tuple[Path, dict[str, Any]]:
    """Build the supplementary H5AD for GSE235863."""
    gse_id = "GSE235863"
    log(f"Building {gse_id}")

    inputs = [
        (
            DOWNLOADS_DIR / gse_id / "suppl" / "GSE235863_five_patients_scRNAseq_cd8t_raw_counts.h5ad.gz",
            DOWNLOADS_DIR / gse_id / "suppl" / "GSE235863_five_patients_TCR.csv.gz",
        ),
        (
            DOWNLOADS_DIR / gse_id / "suppl" / "GSE235863_nine_patients_scRNAseq_cd45_raw_counts.h5ad.gz",
            DOWNLOADS_DIR / gse_id / "suppl" / "GSE235863_nine_patients_TCR.csv.gz",
        ),
    ]

    adatas: list[ad.AnnData] = []
    for gex_path, tcr_path in inputs:
        adata = load_gzipped_h5ad(gex_path)
        obs = adata.obs.copy()
        obs["barcode"] = pd.Index(adata.obs_names.astype(str))
        obs["sample_id"] = normalize_string_series(obs["sample"])
        obs["sample"] = obs["sample_id"]
        obs["library_id"] = obs["sample_id"]
        obs["donor_id"] = normalize_string_series(obs["patient"])
        obs["donor_patient"] = obs["donor_id"]
        obs["sample_type"] = obs["sample_id"]
        obs["tissue"] = normalize_string_series(obs["tissue"]).replace(
            {"P": "blood", "T": "liver tumor", "N": "adjacent normal liver tissue"}
        )
        obs["condition"] = np.select(
            [
                obs["sample_id"].astype(str).str.contains("-pre-", regex=False),
                obs["sample_id"].astype(str).str.contains("-post-", regex=False),
            ],
            ["pre-treatment", "post-treatment"],
            default="",
        )
        obs["treatment"] = "anti-PD-1 plus lenvatinib combination therapy"
        obs["age"] = ""
        obs["sex"] = ""
        obs["enrichment_strategy"] = ""
        obs["assay_type"] = "10x Chromium Single cell 5' + VDJ"
        obs["original_cell_annotation"] = first_nonblank_series(obs, ["sub_cluster", "major_cluster"])
        obs = merge_tcr_into_obs(obs, parse_235863_tcr(tcr_path))
        adata.obs = obs
        adatas.append(adata)

    combined = ad.concat(adatas, join="outer", merge="first", fill_value=0)
    output_path = PER_GSE_OUT_DIR / f"{gse_id}_with_tcr.h5ad"
    summary = finalize_per_gse_adata(
        adata=combined,
        gse_id=gse_id,
        output_path=output_path,
        source_note="five_patients_scRNAseq_cd8t_raw_counts.h5ad.gz + nine_patients_scRNAseq_cd45_raw_counts.h5ad.gz",
        matrix_file=str(DOWNLOADS_DIR / gse_id / "matrix" / f"{gse_id}_series_matrix.txt.gz"),
        platform_id="GPL24676",
        technology_simple="10x 5'",
    )
    return output_path, summary


def build_gse240865(_temp_dir: Path | None = None) -> tuple[Path, dict[str, Any]]:
    """Build the supplementary H5AD for GSE240865."""
    gse_id = "GSE240865"
    log(f"Building {gse_id}")

    sample_prefixes = sorted(
        {
            str(path).rsplit("_barcodes.tsv.gz", 1)[0]
            for path in (DOWNLOADS_DIR / gse_id / "suppl").glob("*_barcodes.tsv.gz")
        }
    )
    adatas: list[ad.AnnData] = []
    for prefix_str in sample_prefixes:
        prefix = Path(prefix_str)
        sample_id = prefix.name.replace(f"{gse_id}_", "")
        adata = read_triplet_with_feature_filter(prefix)
        obs = adata.obs.copy()
        obs["sample_id"] = sample_id
        obs["sample"] = sample_id
        obs["library_id"] = sample_id
        obs["donor_id"] = ""
        obs["donor_patient"] = ""
        obs["sample_type"] = sample_id.replace("_", " ")
        obs["tissue"] = "PBMC"
        obs["condition"] = "filarial" if "filarial" in sample_id.lower() else "control"
        obs["treatment"] = ""
        obs["age"] = ""
        obs["sex"] = ""
        obs["enrichment_strategy"] = "sorted memory CD4 T cells"
        obs["assay_type"] = "10x 5' GEX + feature barcode + VDJ"
        obs["original_cell_annotation"] = ""
        obs = merge_tcr_into_obs(
            obs,
            standardize_tcr_from_contig(
                pd.read_csv(Path(str(prefix) + "_all_contig_annotations.csv.gz"), low_memory=False),
                sample_id=sample_id,
            ),
        )
        adata.obs = obs
        adatas.append(adata)

    combined = ad.concat(adatas, join="outer", merge="first", fill_value=0)
    output_path = PER_GSE_OUT_DIR / f"{gse_id}_with_tcr.h5ad"
    summary = finalize_per_gse_adata(
        adata=combined,
        gse_id=gse_id,
        output_path=output_path,
        source_note="six per-library 10x 5' GEX feature-barcode triplets with matching all_contig_annotations",
        matrix_file=str(DOWNLOADS_DIR / gse_id / "matrix" / f"{gse_id}_series_matrix.txt.gz"),
        platform_id="GPL34284",
        technology_simple="10x 5'",
    )
    return output_path, summary


def build_gse287301(_temp_dir: Path | None = None) -> tuple[Path, dict[str, Any]]:
    """Build the supplementary H5AD for GSE287301."""
    gse_id = "GSE287301"
    log(f"Building {gse_id}")

    pool_order = [
        "chip1pool1",
        "chip1pool2",
        "chip1pool3",
        "chip1pool4",
        "chip1pool5",
        "chip1pool6",
        "chip1pool7",
        "chip1pool8",
        "chip2pool1",
        "chip2pool2",
        "chip2pool3",
        "chip2pool4",
        "chip2pool5",
        "chip2pool6",
        "chip2pool7",
        "chip2pool16",
    ]
    pool_to_gsm = {
        pool: f"GSM87434{42 + idx:02d}" if idx < 10 else f"GSM87434{42 + idx}"
        for idx, pool in enumerate(pool_order)
    }

    patient_matrix = pd.read_csv(DOWNLOADS_DIR / gse_id / "suppl" / "GSE287301_patient_matrix.txt.gz", sep="\t")
    donor_by_pool: dict[str, str] = {}
    for column in [c for c in patient_matrix.columns if c != "Unnamed: 0"]:
        pool_values = patient_matrix[column].astype(str).tolist()
        donor_by_pool[column] = ";".join(pool_values)
        donor_by_pool[f"chip1{column}"] = ";".join(pool_values)
        donor_by_pool[f"chip2{column}"] = ";".join(pool_values)

    matrix, barcodes, features = read_tarred_feature_bc_matrix(
        DOWNLOADS_DIR / gse_id / "suppl" / "GSE287301_filtered_feature_bc_matrix.tar.gz"
    )
    obs = pd.DataFrame(index=barcodes.iloc[:, 0].astype(str))
    obs["barcode"] = obs.index.astype(str)
    suffix = obs["barcode"].astype(str).str.rsplit("-", n=1).str[-1].astype(int)
    obs["sample_id"] = suffix.map(lambda value: pool_order[value - 1])
    obs["sample"] = obs["sample_id"]
    obs["library_id"] = obs["sample_id"].map(lambda value: pool_to_gsm.get(value, value))
    obs["donor_id"] = ""
    obs["donor_patient"] = obs["sample_id"].map(lambda value: donor_by_pool.get(value, ""))
    obs["sample_type"] = "tumor CD3+ T-cell pool"
    obs["tissue"] = "tumor"
    obs["condition"] = "HNSCC"
    obs["treatment"] = ""
    obs["age"] = ""
    obs["sex"] = ""
    obs["enrichment_strategy"] = "CD3-positive tumor T cells"
    obs["assay_type"] = "10x GEM-X 5' GEX + CITE + TCR"
    obs["original_cell_annotation"] = ""

    tcr_frames: list[pd.DataFrame] = []
    for tar_path in sorted((DOWNLOADS_DIR / gse_id / "suppl" / "extracted_analysis").glob("GSM*_*.tar.gz")):
        sample_id = tar_path.name.replace(".tar.gz", "").split("_", 1)[1]
        if sample_id not in pool_order:
            continue
        with tarfile.open(tar_path, "r:gz") as tar:
            member = find_tar_member(tar, "filtered_contig_annotations.csv")
            with tar.extractfile(member) as handle:
                df = pd.read_csv(handle, low_memory=False)
        tcr_frames.append(standardize_tcr_from_contig(df, sample_id=sample_id))
    tcr_df = pd.concat(tcr_frames, ignore_index=True) if tcr_frames else pd.DataFrame()
    obs = merge_tcr_into_obs(obs, tcr_df)

    var = pd.DataFrame(index=features.iloc[:, 1].astype(str))
    var["gene_ids"] = features.iloc[:, 0].astype(str).to_numpy()
    if features.shape[1] >= 3:
        var["feature_types"] = features.iloc[:, 2].astype(str).to_numpy()
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    output_path = PER_GSE_OUT_DIR / f"{gse_id}_with_tcr.h5ad"
    summary = finalize_per_gse_adata(
        adata=adata,
        gse_id=gse_id,
        output_path=output_path,
        source_note="GSE287301 filtered_feature_bc_matrix.tar.gz plus 16 per-pool filtered_contig_annotations tarballs",
        matrix_file=str(DOWNLOADS_DIR / gse_id / "matrix" / f"{gse_id}_series_matrix.txt.gz"),
        platform_id="GPL30173",
        technology_simple="10x 5'",
    )
    return output_path, summary


def build_gse234069(temp_dir: Path) -> tuple[Path, dict[str, Any]]:
    """Build the supplementary H5AD for the 5' lane of GSE234069."""
    gse_id = "GSE234069"
    log(f"Building {gse_id}")
    prefix = temp_dir / gse_id / gse_id
    prefix.parent.mkdir(parents=True, exist_ok=True)
    run_subprocess(
        [
            "Rscript",
            str(R_HELPER),
            "--mode",
            "seurat",
            "--input",
            str(DOWNLOADS_DIR / gse_id / "suppl" / "10x_5" / "GSE234069_data.5prime.n16.rds.gz"),
            "--out-prefix",
            str(prefix),
            "--assay",
            "RNA",
        ]
    )
    matrix, genes, metadata = load_payload(prefix)
    metadata = metadata.set_index("cell_id")

    series_df = parse_geo_series_matrix_samples(gse_id)
    series_df = series_df[series_df.get("seq_type", "").astype(str).str.lower() == "5prime"].copy()
    age_map = dict(zip(series_df["sample_id"], series_df["age"]))
    gsm_map = {
        clean_text(row["sample_id"]): clean_text(row["geo_accession"])
        for _, row in series_df.iterrows()
        if clean_text(row.get("title", "")).endswith("_GEX_scRNA-seq")
    }

    obs = metadata.copy()
    obs["barcode"] = obs.index.to_series().map(extract_10x_barcode).to_numpy()
    sample_id = obs_series(obs, "orig.ident")
    if (sample_id == "").all():
        sample_id = obs_series(obs, "sample")
    donor_id = sample_id.str.replace("FRZ$", "", regex=True)

    parsed_pairs = obs_series(obs, "cdr3s_aa").map(parse_cdr3_pairs)
    obs["sample_id"] = sample_id
    obs["sample"] = sample_id
    obs["library_id"] = sample_id.map(lambda value: gsm_map.get(value, value))
    obs["donor_id"] = donor_id
    obs["donor_patient"] = donor_id
    obs["sample_type"] = "CSF 5prime"
    obs["tissue"] = "Cerebrospinal fluid (CSF)"
    obs["condition"] = obs_series(obs, "group")
    obs["treatment"] = ""
    obs["age"] = donor_id.map(lambda value: clean_text(age_map.get(value, "")))
    obs["sex"] = ""
    obs["enrichment_strategy"] = ""
    obs["assay_type"] = "10x Chromium 5' GEX + VDJ"
    obs["original_cell_annotation"] = first_nonblank_series(obs, ["annotation", "annotation_L1", "annotation_L2", "annotation_L3"])
    obs["TRA_cdr3"] = parsed_pairs.map(lambda item: item["TRA_cdr3"]).astype(object)
    obs["TRB_cdr3"] = parsed_pairs.map(lambda item: item["TRB_cdr3"]).astype(object)
    obs["TRA_v"] = ""
    obs["TRA_d"] = ""
    obs["TRA_j"] = ""
    obs["TRA_cdr3_nt"] = ""
    obs["TRA_clone_id"] = obs_series(obs, "clonotype_id")
    obs["TRA_umis"] = 0
    obs["TRA_reads"] = 0
    obs["TRB_v"] = ""
    obs["TRB_d"] = ""
    obs["TRB_j"] = ""
    obs["TRB_cdr3_nt"] = ""
    obs["TRB_clone_id"] = obs_series(obs, "clonotype_id")
    obs["TRB_umis"] = 0
    obs["TRB_reads"] = 0
    obs["TCRseq"] = np.where(
        (normalize_string_series(obs["TRA_cdr3"]) != "")
        | (normalize_string_series(obs["TRB_cdr3"]) != ""),
        "yes",
        "no",
    )

    var = pd.DataFrame(index=genes["gene_name"].astype(str))
    var["gene_ids"] = genes["gene_id"].astype(str).to_numpy()
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    output_path = PER_GSE_OUT_DIR / f"{gse_id}_with_tcr.h5ad"
    summary = finalize_per_gse_adata(
        adata=adata,
        gse_id=gse_id,
        output_path=output_path,
        source_note="GSE234069_data.5prime.n16.rds.gz from suppl/10x_5 only",
        matrix_file=str(DOWNLOADS_DIR / gse_id / "matrix" / f"{gse_id}_series_matrix.txt.gz"),
        platform_id="GPL24676",
        technology_simple="10x 5'",
    )
    return output_path, summary


def build_gse287541(temp_dir: Path) -> tuple[Path, dict[str, Any]]:
    """Build the supplementary H5AD for GSE287541."""
    gse_id = "GSE287541"
    log(f"Building {gse_id}")
    prefix = temp_dir / gse_id / gse_id
    prefix.parent.mkdir(parents=True, exist_ok=True)
    run_subprocess(
        [
            "Rscript",
            str(R_HELPER),
            "--mode",
            "seurat",
            "--input",
            str(DOWNLOADS_DIR / gse_id / "suppl" / "GSE287541_all_pbmc_combined_filtered_SeuratObj.RDS.gz"),
            "--out-prefix",
            str(prefix),
            "--assay",
            "RNA",
        ]
    )
    matrix, genes, metadata = load_payload(prefix)
    metadata = metadata.set_index("cell_id")

    obs = metadata.copy()
    obs["barcode"] = obs.index.to_series().map(extract_10x_barcode).to_numpy()
    obs["sample_id"] = obs_series(obs, "PID_visit")
    obs["sample"] = obs["sample_id"]
    library_series = obs_series(obs, "orig.ident")
    library_series = library_series.where(library_series != "", obs["sample_id"])
    obs["library_id"] = library_series
    donor = obs_series(obs, "Person.CHILI.ID")
    if (donor == "").all():
        donor = obs["sample_id"].str.extract(r"^(N\d+)").fillna("").iloc[:, 0].astype(object)
    obs["donor_id"] = donor
    obs["donor_patient"] = donor
    obs["sample_type"] = "PBMC"
    obs["tissue"] = "PBMC"
    obs["condition"] = obs_series(obs, "Group")
    obs["treatment"] = obs_series(obs, "CPI.regime")
    obs["age"] = obs_series(obs, "Age")
    obs["sex"] = obs_series(obs, "Gender")
    obs["enrichment_strategy"] = ""
    obs["assay_type"] = "10x Chromium Single Cell 5' Gene Expression and V(D)J"
    obs["original_cell_annotation"] = first_nonblank_series(
        obs,
        [
            "predicted.celltype.l3",
            "predicted.celltype.l2",
            "predicted.celltype.l1",
            "seurat_clusters",
        ],
    )

    paired = obs_series(obs, "cdr3s_aa").map(parse_cdr3_pairs)
    chain = obs_series(obs, "chain").str.upper()
    obs["TRA_cdr3"] = paired.map(lambda item: item["TRA_cdr3"]).astype(object)
    obs["TRB_cdr3"] = paired.map(lambda item: item["TRB_cdr3"]).astype(object)
    cdr3_series = obs_series(obs, "cdr3")
    v_series = obs_series(obs, "v_gene")
    j_series = obs_series(obs, "j_gene")
    obs.loc[(obs["TRA_cdr3"] == "") & (chain == "TRA"), "TRA_cdr3"] = cdr3_series.loc[(obs["TRA_cdr3"] == "") & (chain == "TRA")]
    obs.loc[(obs["TRB_cdr3"] == "") & (chain == "TRB"), "TRB_cdr3"] = cdr3_series.loc[(obs["TRB_cdr3"] == "") & (chain == "TRB")]
    obs["TRA_v"] = np.where(chain == "TRA", v_series, "")
    obs["TRA_d"] = ""
    obs["TRA_j"] = np.where(chain == "TRA", j_series, "")
    obs["TRA_cdr3_nt"] = ""
    obs["TRA_clone_id"] = obs_series(obs, "clonotype_id")
    obs["TRA_umis"] = 0
    obs["TRA_reads"] = 0
    obs["TRB_v"] = np.where(chain == "TRB", v_series, "")
    obs["TRB_d"] = ""
    obs["TRB_j"] = np.where(chain == "TRB", j_series, "")
    obs["TRB_cdr3_nt"] = ""
    obs["TRB_clone_id"] = obs_series(obs, "clonotype_id")
    obs["TRB_umis"] = 0
    obs["TRB_reads"] = 0
    obs["TCRseq"] = np.where(
        (normalize_string_series(obs["TRA_cdr3"]) != "")
        | (normalize_string_series(obs["TRB_cdr3"]) != ""),
        "yes",
        "no",
    )

    var = pd.DataFrame(index=genes["gene_name"].astype(str))
    var["gene_ids"] = genes["gene_id"].astype(str).to_numpy()
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    output_path = PER_GSE_OUT_DIR / f"{gse_id}_with_tcr.h5ad"
    summary = finalize_per_gse_adata(
        adata=adata,
        gse_id=gse_id,
        output_path=output_path,
        source_note="GSE287541_all_pbmc_combined_filtered_SeuratObj.RDS.gz",
        matrix_file=str(DOWNLOADS_DIR / gse_id / "matrix" / f"{gse_id}_series_matrix.txt.gz"),
        platform_id="GPL30173",
        technology_simple="10x 5'",
    )
    return output_path, summary


def build_all_per_gse_h5ads(force_rebuild: bool, selected_gses: list[str] | None = None) -> pd.DataFrame:
    """Build all six supplementary per-GSE H5AD files and return the registry."""
    ensure_dirs()
    builder_map = {
        "GSE179994": build_gse179994,
        "GSE235863": build_gse235863,
        "GSE240865": build_gse240865,
        "GSE287301": build_gse287301,
        "GSE234069": build_gse234069,
        "GSE287541": build_gse287541,
    }
    target_gses = selected_gses or SUPP_GSES
    builders = [builder_map[gse_id] for gse_id in target_gses]
    summaries: list[dict[str, Any]] = []
    registry_rows: list[dict[str, str]] = []

    with tempfile.TemporaryDirectory(dir=SUPP_TMP_DIR) as temp_dir_name:
        temp_dir = Path(temp_dir_name)
        for builder in builders:
            gse_id = builder.__name__.replace("build_", "").upper()
            output_path = PER_GSE_OUT_DIR / f"{gse_id}_with_tcr.h5ad"
            if output_path.exists() and not force_rebuild:
                log(f"Reusing existing {gse_id} per-GSE H5AD: {output_path}")
                adata = ad.read_h5ad(output_path, backed="r")
                try:
                    summaries.append(
                        {
                            "gse_id": gse_id,
                            "h5ad_path": str(output_path),
                            "n_obs": int(adata.n_obs),
                            "n_vars": int(adata.n_vars),
                            "samples": int(normalize_string_series(adata.obs.get("sample_id", pd.Series([], dtype=object))).replace("", np.nan).nunique(dropna=True)),
                            "libraries": int(normalize_string_series(adata.obs.get("library_id", pd.Series([], dtype=object))).replace("", np.nan).nunique(dropna=True)),
                            "tcr_yes_cells": int((normalize_string_series(adata.obs.get("TCRseq", pd.Series([], dtype=object))).str.lower() == "yes").sum()),
                        }
                    )
                finally:
                    if getattr(adata, "file", None) is not None:
                        adata.file.close()
            else:
                output_path, summary = builder(temp_dir)
                summaries.append(summary)

            registry_rows.append(
                {
                    "h5ad_path": str(output_path),
                    "gse_id": gse_id,
                    "source_root": str(PER_GSE_OUT_DIR),
                }
            )

    build_summary = pd.DataFrame(summaries).sort_values("gse_id").reset_index(drop=True)
    build_summary.to_csv(SUPP_BUILD_SUMMARY_CSV, index=False)
    registry = pd.DataFrame(registry_rows).sort_values("gse_id").reset_index(drop=True)
    registry.to_csv(SUPP_REGISTRY_CSV, index=False)
    log(f"Wrote {SUPP_BUILD_SUMMARY_CSV}")
    log(f"Wrote {SUPP_REGISTRY_CSV}")
    return registry


def export_base_metadata(registry: pd.DataFrame) -> pd.DataFrame:
    """Export supplementary all-cell metadata from the six per-GSE H5AD files."""
    rows: list[pd.DataFrame] = []
    for _, record in registry.iterrows():
        gse_id = record["gse_id"]
        h5ad_path = Path(record["h5ad_path"])
        adata = ad.read_h5ad(h5ad_path, backed="r")
        try:
            obs = adata.obs.copy()
        finally:
            if getattr(adata, "file", None) is not None:
                adata.file.close()

        frame = pd.DataFrame(index=obs.index.astype(str))
        frame["gse_id"] = gse_id
        frame["cell_id"] = frame.index.astype(str)
        frame["source_h5ad"] = str(h5ad_path)
        frame["source_root"] = record["source_root"]

        for column in BASE_METADATA_COLS:
            if column in {"gse_id", "cell_id", "source_h5ad", "source_root"}:
                continue
            if column in obs.columns:
                frame[column] = normalize_string_series(obs[column])
            else:
                frame[column] = ""
        rows.append(frame.reset_index(drop=True))

    base_metadata = pd.concat(rows, ignore_index=True)
    base_metadata.to_csv(SUPP_BASE_METADATA_CSV, index=False, compression="gzip")
    log(f"Wrote {SUPP_BASE_METADATA_CSV}")
    return base_metadata


def configure_phase0_outputs() -> None:
    """Point the shared Phase 0 helper at the supplementary output package."""
    phase0.FIGURE_DIR = SUPP_FIG_DIR
    phase0.TABLE_DIR = SUPP_TABLE_DIR
    phase0.LOG_DIR = SUPP_LOG_DIR
    phase0.PHASE0_AUDIT_CSV = SUPP_TABLE_DIR / "phase0_dataset_audit.csv"
    phase0.PHASE0_CATEGORY_CSV = SUPP_TABLE_DIR / "phase0_category_summary.csv"
    phase0.PHASE0_QC_MD = SUPP_LOG_DIR / "phase0_qc_summary.md"


def configure_phase1_outputs() -> None:
    """Point the shared Phase 1 helper at the supplementary output package."""
    phase1.OUTPUT_ROOT = OUTPUT_ROOT
    phase1.FIGURE_DIR = SUPP_FIG_DIR
    phase1.TABLE_DIR = SUPP_TABLE_DIR
    phase1.LOG_DIR = SUPP_LOG_DIR
    phase1.TMP_DIR = SUPP_PHASE1_TMP_DIR
    phase1.TNK_CANDIDATES_H5AD = SUPP_CANDIDATES_H5AD
    phase1.AUDIT_CSV = SUPP_TABLE_DIR / "phase0_dataset_audit.csv"
    phase1.PHASE1_SUMMARY_CSV = SUPP_TABLE_DIR / "phase1_categoryA_selection_summary.csv"
    phase1.PHASE1_MARKER_CSV = SUPP_TABLE_DIR / "phase1_categoryA_marker_availability.csv"
    phase1.PHASE1_QC_MD = SUPP_LOG_DIR / "phase1_qc_summary.md"


def run_supplementary_phase0(registry_path: Path) -> pd.DataFrame:
    """Run supplementary Phase 0 on the six-dataset registry."""
    configure_phase0_outputs()
    phase0.ensure_output_dirs()
    audit = phase0.run_phase0_audit(registry_path)
    return audit


def run_supplementary_phase1() -> pd.DataFrame:
    """Run supplementary Phase 1 after explicit approval."""
    configure_phase1_outputs()
    phase1.ensure_output_dirs()
    phase1.main()
    qc_text = (SUPP_LOG_DIR / "phase1_qc_summary.md").read_text(encoding="utf-8")
    qc_text = qc_text.replace("TNK_candidates.h5ad", "TNK_candidates_supp.h5ad")
    (SUPP_LOG_DIR / "phase1_qc_summary.md").write_text(qc_text, encoding="utf-8")
    return pd.read_csv(SUPP_TABLE_DIR / "phase1_categoryA_selection_summary.csv")


def build_candidate_metadata(base_metadata: pd.DataFrame) -> pd.DataFrame:
    """Build a candidate-only supplementary metadata CSV aligned to the current header."""
    header = list(pd.read_csv(HARMONIZED_HEADER_SOURCE, nrows=0).columns)
    adata = ad.read_h5ad(SUPP_CANDIDATES_H5AD, backed="r")
    try:
        obs = adata.obs.copy()
    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()

    obs_export = obs.copy()
    obs_export.insert(0, "obs_name", obs_export.index.astype(str))
    obs_export["project name"] = normalize_string_series(obs_export["source_gse_id"])
    obs_export["sampleid"] = normalize_string_series(obs_export["phase1_sample_label"])
    obs_export["barcodes"] = normalize_string_series(obs_export["original_cell_id"])

    merged = obs_export.merge(
        base_metadata,
        left_on=["source_gse_id", "original_cell_id"],
        right_on=["gse_id", "cell_id"],
        how="left",
        validate="one_to_one",
    )

    if merged["cell_id"].isna().any():
        missing = merged.loc[merged["cell_id"].isna(), ["source_gse_id", "original_cell_id"]].head(20)
        raise ValueError(
            "Supplementary candidate metadata join failed for some Phase 1 cells.\n"
            f"Preview:\n{missing.to_string(index=False)}"
        )

    merged["sampleid_obs"] = normalize_string_series(merged["sampleid"])
    merged["project name_obs"] = normalize_string_series(merged["project name"])
    merged["barcodes_obs"] = normalize_string_series(merged["barcodes"])
    merged["obs_name_obs"] = normalize_string_series(merged["obs_name"])
    for column in [
        "original_cell_id",
        "source_gse_id",
        "phase1_donor_label",
        "phase1_sample_label",
        "phase1_library_label",
        "phase1_annotation_label",
        "phase1_annotation_keep",
        "phase1_marker_keep",
        "phase1_keep",
        "phase1_selection_reason",
        "phase1_t_hits",
        "phase1_gd_hits",
        "phase1_nk_hits",
        "phase1_nk_strong_hits",
        "phase1_t_score",
        "phase1_gd_score",
        "phase1_nk_score",
        "phase1_contam_score",
        "source_gse_id_concat",
    ]:
        obs_column = f"{column}_obs"
        if column in merged.columns and obs_column not in merged.columns:
            merged[obs_column] = merged[column]

    merged["sampleid"] = normalize_phase1_sampleid(merged["sample_id"], merged["sampleid"])
    merged["source_root"] = normalize_string_series(merged["source_root"])
    merged["tcr_availability"] = np.where(
        normalize_string_series(merged["TCRseq"]).str.lower() == "yes",
        "yes",
        normalize_string_series(merged["tcr_availability"]),
    )

    for column in header:
        if column not in merged.columns:
            merged[column] = ""

    final = merged[header].copy()
    final.to_csv(SUPP_METADATA_CANDIDATE_CSV, index=False)
    log(f"Wrote {SUPP_METADATA_CANDIDATE_CSV}")
    return final


def main() -> None:
    """Run the supplementary build through the requested stop point."""
    args = parse_args()
    ensure_dirs()

    registry = build_all_per_gse_h5ads(force_rebuild=args.force_rebuild, selected_gses=args.gse)
    if args.build_only:
        log("")
        log("Supplementary build-only run complete.")
        log(f"Registry rows written: {len(registry)}")
        return

    base_metadata = export_base_metadata(registry)
    audit = run_supplementary_phase0(SUPP_REGISTRY_CSV)

    log("")
    log("Supplementary Phase 0 complete.")
    log(f"Registry rows audited: {len(audit)}")
    log(f"Category counts: {audit['phase0_category'].value_counts().to_dict()}")

    if args.stop_after == "phase0":
        log("Stopping after supplementary Phase 0 per QC gate.")
        return

    phase1_summary = run_supplementary_phase1()
    build_candidate_metadata(base_metadata)
    log("")
    log("Supplementary Phase 1 complete.")
    log(f"Candidate rows: {int(phase1_summary['candidate_n_obs'].sum())}")


if __name__ == "__main__":
    main()
