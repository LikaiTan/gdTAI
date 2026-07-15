#!/usr/bin/env python3
"""Rebuild raw TCR calls for flagged per-GSE H5ADs and rerun Phase 4-style QC.

This helper is intentionally conservative:

1. discover the raw TCR source under ``downloads/<GSE>/suppl``
2. rebuild one canonical per-cell TCR table using a sample-aware join key
3. delete the old TCR fields in the target H5AD and replace them from raw TCR
4. recompute Phase 4-style TRAB/TRD continuous scores from a temporary
   normalize-total + log1p copy of count-space ``X``
5. write a TRAB-vs-TRD scatter plot colored by whether the cell has paired
   TRA/TRB CDR3 calls

The script supports dry-run mode by default. Use ``--write-h5ad`` to persist the
repaired TCR fields and Phase 4 score columns back into each per-GSE H5AD.
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
import gzip
import io
import logging
import math
import re
import shutil
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse
from sklearn.utils import sparsefuncs


PROJECT_ROOT = Path(__file__).resolve().parents[2]
REGISTRY_CSV = PROJECT_ROOT / "configs" / "datasets" / "integration_inputs.csv"
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
TABLE_DIR = OUTPUT_ROOT / "tables" / "tcr_rebuild_phase4"
FIGURE_DIR = OUTPUT_ROOT / "figures" / "tcr_rebuild_phase4"
LOG_DIR = OUTPUT_ROOT / "logs"
LOG_PATH = LOG_DIR / "tcr_rebuild_phase4.log"
SUMMARY_CSV = TABLE_DIR / "tcr_rebuild_phase4_summary.csv"
SUMMARY_MD = LOG_DIR / "tcr_rebuild_phase4_qc.md"

TARGET_GSES = [
    "GSE161918",
    "GSE168859",
    "GSE188620",
    "GSE227709",
    "GSE243572",
    "GSE243905",
    "GSE254249",
    "GSE212217",
    "GSE311112",
    "GSE308075",
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

TARGET_SUM = 10_000.0
N_BINS = 25
CTRL_SIZE = 50
RANDOM_STATE = 1
PHASE4_SCORE_COLUMNS = {
    "tra": "phase4_tra_score",
    "trb": "phase4_trb_score",
    "trab": "phase4_trab_score",
    "trd": "phase4_trd_score",
    "trd_minus_trab": "phase4_trd_minus_trab",
}
PHASE4_SCALED_SCORE_COLUMNS = {
    "trd_scaled": "phase4_trd_score_scaled",
    "trab_scaled": "phase4_trab_score_scaled",
    "trd_minus_trab_scaled": "phase4_trd_minus_trab_scaled",
}
MODULE_PATTERNS = {
    "tra": re.compile(r"^TRAC$|^TRAV|^TRAJ"),
    "trb": re.compile(r"^TRBC|^TRBV|^TRBJ"),
    "trd": re.compile(r"^TRDC$|^TRDV|^TRDJ"),
}

PAIRING_SCATTER_SAMPLE_SIZE = 80_000
FIGURE_DPI = 300
PHASE4_REPAIR_COLUMNS = [
    *PHASE4_SCORE_COLUMNS.values(),
    *PHASE4_SCALED_SCORE_COLUMNS.values(),
    "phase4_trab_minus_trd",
    "phase4_trab_minus_trd_scaled",
]
BACKUP_SUFFIX = ".bak_before_raw_tcr_repair"
RAW_CONTIG_USECOLS = lambda column: column in {
    "barcode",
    "chain",
    "productive",
    "cdr3",
    "v_gene",
    "d_gene",
    "j_gene",
    "cdr3_nt",
    "raw_clonotype_id",
    "clonotype_id",
    "clone_id",
    "umis",
    "reads",
}


@dataclass
class RawTcrBundle:
    """Canonical raw TCR payload for one dataset."""

    gse_id: str
    source_label: str
    source_files: list[str]
    tcr_df: pd.DataFrame


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gse",
        nargs="+",
        default=TARGET_GSES,
        help="Flagged GSE ids to process.",
    )
    parser.add_argument(
        "--registry",
        type=Path,
        default=REGISTRY_CSV,
        help="Registry CSV with per-GSE H5AD paths.",
    )
    parser.add_argument(
        "--write-h5ad",
        action="store_true",
        help="Persist repaired TCR and Phase 4 score columns into each target H5AD.",
    )
    parser.add_argument(
        "--scatter-sample-size",
        type=int,
        default=PAIRING_SCATTER_SAMPLE_SIZE,
        help="Maximum number of cells to draw in each scatter plot.",
    )
    return parser.parse_args()


def configure_logging() -> None:
    """Configure console and file logging."""
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    handlers = [
        logging.FileHandler(LOG_PATH, mode="a", encoding="utf-8"),
        logging.StreamHandler(),
    ]
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
        force=True,
    )


def ensure_dirs() -> None:
    """Create output directories."""
    for path in [TABLE_DIR, FIGURE_DIR, LOG_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def normalize_key(text: object) -> str:
    """Normalize sample-like text into a stable lowercase key."""
    value = clean_text(text).lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def clean_text(value: object) -> str:
    """Normalize missing-like values to an empty string."""
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "na", "none", "null", "<na>"}:
        return ""
    return text


def normalize_string_series(series: pd.Series) -> pd.Series:
    """Return a plain object string series with blank missing values."""
    return series.map(clean_text).astype(object)


def normalize_barcode_core(value: object) -> str:
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


def extract_10x_barcode(value: object) -> str:
    """Extract a full 10x barcode with lane suffix if present."""
    text = clean_text(value)
    match = re.search(r"([ACGTN]+-\d+)", text, flags=re.IGNORECASE)
    return match.group(1).upper() if match else text


def write_h5ad_atomic(adata: ad.AnnData, output_path: Path) -> None:
    """Write a dataset H5AD atomically."""
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    if tmp_path.exists():
        tmp_path.unlink()
    adata.write_h5ad(tmp_path)
    tmp_path.replace(output_path)


def first_nonempty_series_value(series: pd.Series) -> str:
    """Return the first non-empty normalized value from a series."""
    values = normalize_string_series(series)
    non_empty = values[values != ""]
    return non_empty.iloc[0] if len(non_empty) else ""


def aggregate_tcr_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse duplicate TCR rows by sample plus barcode core."""
    if df.empty:
        return pd.DataFrame(columns=["sample_id", "barcode_core", *TCR_COLS])
    keep_cols = ["sample_id", "barcode_core", *TCR_COLS]
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
        return pd.DataFrame(columns=["sample_id", "barcode_core", *TCR_COLS])
    work["has_paired_trab"] = work["TRA_cdr3"].ne("") & work["TRB_cdr3"].ne("")
    work["has_any_tcr"] = work["has_paired_trab"] | work["TRA_cdr3"].ne("") | work["TRB_cdr3"].ne("")
    work["numeric_support"] = work[numeric_cols].sum(axis=1)
    selected = (
        work.sort_values(
            ["sample_id", "barcode_core", "has_paired_trab", "has_any_tcr", "numeric_support"],
            ascending=[True, True, False, False, False],
            kind="mergesort",
        )
        .drop_duplicates(subset=["sample_id", "barcode_core"], keep="first")
        .drop(columns=["has_paired_trab", "has_any_tcr", "numeric_support"])
    )
    return selected.reset_index(drop=True)


def empty_tcr_frame(index: pd.Index) -> pd.DataFrame:
    """Return an empty TCR dataframe on one index."""
    frame = pd.DataFrame(index=index)
    for column in TCR_COLS:
        frame[column] = ""
    frame["TCRseq"] = "no"
    return frame


def merge_tcr_into_obs(obs: pd.DataFrame, tcr_df: pd.DataFrame) -> pd.DataFrame:
    """Merge a pre-aggregated TCR table into one obs dataframe."""
    merged = obs.copy()
    merged["barcode_core"] = merged["barcode"].map(normalize_barcode_core)
    merged["sample_id"] = normalize_string_series(merged["sample_id"])
    if tcr_df.empty:
        return merged.join(empty_tcr_frame(merged.index)).drop(columns=["barcode_core"])

    joined = (
        merged.reset_index(names="obs_index")
        .merge(tcr_df, on=["sample_id", "barcode_core"], how="left", suffixes=("", "_tcr"))
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


def find_module_genes(var_names: pd.Index) -> dict[str, pd.Index]:
    """Find the exact Phase 4 TRA/TRB/TRD module genes."""
    modules: dict[str, pd.Index] = {}
    for name, pattern in MODULE_PATTERNS.items():
        genes = [gene for gene in var_names if pattern.match(str(gene))]
        if len(genes) == 0:
            raise ValueError(f"Module `{name}` has no genes in the current object.")
        modules[name] = pd.Index(sorted(set(genes)), dtype="string")
    modules["trab"] = pd.Index(sorted(set(modules["tra"]).union(set(modules["trb"]))), dtype="string")
    return modules


def pick_control_genes(
    gene_list: pd.Index,
    gene_pool: pd.Index,
    gene_means: np.ndarray,
    *,
    random_state: int,
) -> pd.Index:
    """Reproduce the Scanpy/Seurat control-gene binning and sampling logic."""
    np.random.seed(random_state)
    obs_avg = pd.Series(gene_means, index=gene_pool)
    obs_avg = obs_avg[np.isfinite(obs_avg)]
    n_items = int(np.round(len(obs_avg) / (N_BINS - 1)))
    if n_items <= 0:
        raise ValueError("Invalid control-gene bin size during Phase 4 setup.")
    obs_cut = obs_avg.rank(method="min") // n_items
    control_genes = pd.Index([], dtype="string")
    for cut in np.unique(obs_cut.loc[gene_list]):
        r_genes = obs_cut[obs_cut == cut].index
        if len(r_genes) == 0:
            continue
        if CTRL_SIZE < len(r_genes):
            r_genes = r_genes.to_series().sample(CTRL_SIZE).index
        r_genes = r_genes.difference(gene_list)
        control_genes = control_genes.union(r_genes)
    if len(control_genes) == 0:
        raise RuntimeError(f"No control genes found for module with {len(gene_list)} genes.")
    return control_genes


def score_chunk(chunk: sparse.csr_matrix, gene_idx: np.ndarray, ctrl_idx: np.ndarray) -> np.ndarray:
    """Score one module on one normalized/log1p CSR matrix."""
    gene_mean = np.asarray(chunk[:, gene_idx].sum(axis=1)).ravel() / float(len(gene_idx))
    ctrl_mean = np.asarray(chunk[:, ctrl_idx].sum(axis=1)).ravel() / float(len(ctrl_idx))
    return (gene_mean - ctrl_mean).astype(np.float32, copy=False)


def minmax_scale(values: np.ndarray) -> np.ndarray:
    """Scale one score vector into the 0-1 range."""
    value_min = float(np.min(values))
    value_max = float(np.max(values))
    if math.isclose(value_min, value_max):
        return np.zeros_like(values, dtype=np.float32)
    return ((values - value_min) / (value_max - value_min)).astype(np.float32, copy=False)


def add_scaled_scores(scores: dict[str, np.ndarray]) -> tuple[dict[str, np.ndarray], dict[str, dict[str, float]]]:
    """Add scaled TRD/TRAB scores and scaled TRD-TRAB difference."""
    trd_scaled = minmax_scale(scores["trd"])
    trab_scaled = minmax_scale(scores["trab"])
    scores["trd_scaled"] = trd_scaled
    scores["trab_scaled"] = trab_scaled
    scores["trd_minus_trab_scaled"] = (trd_scaled - trab_scaled).astype(np.float32, copy=False)
    return scores, {}


def load_registry_map(registry_csv: Path) -> dict[str, Path]:
    """Load the per-GSE H5AD path map from the registry."""
    registry = pd.read_csv(registry_csv, dtype=str)
    required = {"gse_id", "h5ad_path"}
    missing = required.difference(registry.columns)
    if missing:
        raise ValueError(f"Registry is missing required columns: {sorted(missing)}")
    mapping = {
        clean_text(row.gse_id): Path(clean_text(row.h5ad_path))
        for row in registry.itertuples(index=False)
        if clean_text(row.gse_id)
    }
    return mapping


def read_csv_from_handle(handle, name: str) -> pd.DataFrame:
    """Read a CSV or CSV.GZ stream from a tarfile handle."""
    if name.lower().endswith(".gz"):
        with gzip.open(handle, "rt", encoding="utf-8", errors="ignore") as text_handle:
            return pd.read_csv(text_handle, low_memory=False, usecols=RAW_CONTIG_USECOLS)
    return pd.read_csv(handle, low_memory=False, usecols=RAW_CONTIG_USECOLS)


def truthy_series(series: pd.Series) -> pd.Series:
    """Normalize a productive-like column to boolean."""
    values = series.map(clean_text).str.lower()
    return values.isin({"true", "t", "1", "yes", "y"})


def filter_productive(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only productive contigs when the source provides that flag."""
    if "productive" not in df.columns:
        return df.copy()
    return df.loc[truthy_series(df["productive"])].copy()


def first_nonempty_value(series: pd.Series) -> str:
    """Return the first non-empty string in a series."""
    values = normalize_string_series(series)
    non_empty = values[values != ""]
    return non_empty.iloc[0] if len(non_empty) else ""


def column_or_blank(df: pd.DataFrame, column: str) -> pd.Series:
    """Return one dataframe column or a same-length blank object series."""
    if column in df.columns:
        return df[column]
    return pd.Series([""] * len(df), index=df.index, dtype=object)


def standardize_contig_table(df: pd.DataFrame, sample_keys: pd.Series | str) -> pd.DataFrame:
    """Standardize 10x contig rows into one canonical per-cell TCR table."""
    if "barcode" not in df.columns or "chain" not in df.columns:
        return pd.DataFrame(columns=["sample_id", "barcode_core", *TCR_COLS])

    work = pd.DataFrame(index=df.index)
    if isinstance(sample_keys, pd.Series):
        work["sample_id"] = sample_keys.map(normalize_key)
    else:
        work["sample_id"] = normalize_key(sample_keys)
    work["barcode_core"] = df["barcode"].map(extract_10x_barcode).map(normalize_barcode_core)
    work["chain"] = df["chain"].astype(str).str.upper()
    work["cdr3"] = column_or_blank(df, "cdr3").map(clean_text)
    work["v_gene"] = column_or_blank(df, "v_gene").map(clean_text)
    work["d_gene"] = column_or_blank(df, "d_gene").map(clean_text)
    work["j_gene"] = column_or_blank(df, "j_gene").map(clean_text)
    work["cdr3_nt"] = column_or_blank(df, "cdr3_nt").map(clean_text)
    work["clone_id"] = (
        column_or_blank(
            df,
            "raw_clonotype_id" if "raw_clonotype_id" in df.columns else "clonotype_id" if "clonotype_id" in df.columns else "clone_id",
        ).map(clean_text)
    )
    work["umis"] = pd.to_numeric(column_or_blank(df, "umis"), errors="coerce").fillna(0).astype(int)
    work["reads"] = pd.to_numeric(column_or_blank(df, "reads"), errors="coerce").fillna(0).astype(int)
    work = work.loc[(work["sample_id"] != "") & (work["barcode_core"] != "")]
    work = work.loc[work["chain"].isin(["TRA", "TRB"])].copy()
    if work.empty:
        return pd.DataFrame(columns=["sample_id", "barcode_core", *TCR_COLS])

    # Rank rows so the first row per cell/chain keeps the most informative contig.
    work["has_cdr3"] = work["cdr3"].ne("")
    work["has_cdr3_nt"] = work["cdr3_nt"].ne("")
    selected = (
        work.sort_values(
            ["sample_id", "barcode_core", "chain", "has_cdr3", "has_cdr3_nt", "umis", "reads"],
            ascending=[True, True, True, False, False, False, False],
            kind="mergesort",
        )
        .drop_duplicates(subset=["sample_id", "barcode_core", "chain"], keep="first")
        .drop(columns=["has_cdr3", "has_cdr3_nt"])
    )

    base = selected[["sample_id", "barcode_core"]].drop_duplicates().copy()
    tra = selected.loc[
        selected["chain"] == "TRA",
        ["sample_id", "barcode_core", "cdr3", "v_gene", "d_gene", "j_gene", "cdr3_nt", "clone_id", "umis", "reads"],
    ].rename(
        columns={
            "cdr3": "TRA_cdr3",
            "v_gene": "TRA_v",
            "d_gene": "TRA_d",
            "j_gene": "TRA_j",
            "cdr3_nt": "TRA_cdr3_nt",
            "clone_id": "TRA_clone_id",
            "umis": "TRA_umis",
            "reads": "TRA_reads",
        }
    )
    trb = selected.loc[
        selected["chain"] == "TRB",
        ["sample_id", "barcode_core", "cdr3", "v_gene", "d_gene", "j_gene", "cdr3_nt", "clone_id", "umis", "reads"],
    ].rename(
        columns={
            "cdr3": "TRB_cdr3",
            "v_gene": "TRB_v",
            "d_gene": "TRB_d",
            "j_gene": "TRB_j",
            "cdr3_nt": "TRB_cdr3_nt",
            "clone_id": "TRB_clone_id",
            "umis": "TRB_umis",
            "reads": "TRB_reads",
        }
    )
    out = base.merge(tra, on=["sample_id", "barcode_core"], how="left")
    out = out.merge(trb, on=["sample_id", "barcode_core"], how="left")

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
    return out[["sample_id", "barcode_core", *TCR_COLS]]


def sample_key_from_obs(gse_id: str, obs: pd.DataFrame) -> pd.Series:
    """Build one safe sample key per cell from the current H5AD obs."""
    if gse_id == "GSE161918":
        sample = obs["sample"].astype(str)
        matches = sample.str.extract(r"_(b\d+)_10xlane0*(\d+)_", flags=re.IGNORECASE)
        return matches.apply(
            lambda row: normalize_key(f"{row.iloc[0]}_10xlane{int(row.iloc[1])}")
            if pd.notna(row.iloc[0]) and pd.notna(row.iloc[1])
            else "",
            axis=1,
        )

    if gse_id == "GSE168859":
        return obs["sample_id"].map(normalize_key)

    if gse_id == "GSE254249":
        barcode = obs["barcode"].astype(str)
        return barcode.str.replace(r"^[ACGTN]+-\d+_", "", regex=True).map(normalize_key)

    sample_column = "sample_id" if "sample_id" in obs.columns and gse_id not in {"GSE243905"} else "sample"
    sample = obs[sample_column].astype(str)

    if gse_id in {"GSE188620", "GSE227709", "GSE243572", "GSE243905", "GSE212217", "GSE308075"}:
        cleaned = sample.str.replace(r"^gsm\d+_", "", regex=True)
        cleaned = cleaned.str.replace(r"_(matrix|gex|scrna)$", "", regex=True)
        return cleaned.map(normalize_key)

    if gse_id == "GSE311112":
        cleaned = sample.str.replace(r"^gsm\d+_", "", regex=True)
        extracted = cleaned.str.extract(r"(?:pt_|cll_)?(\d+_(?:baseline|3yrs|5yrs))", expand=False)
        return extracted.map(normalize_key)

    raise ValueError(f"Unsupported obs sample-key rule for {gse_id}")


def barcode_series_from_obs(obs: pd.DataFrame) -> pd.Series:
    """Return one usable barcode-like series from obs."""
    for column in ["barcode", "cell_barcode"]:
        if column in obs.columns:
            values = obs[column].astype(str).map(clean_text)
            if values.ne("").any():
                return values
    return pd.Series(obs.index.astype(str), index=obs.index, dtype=object)


def merge_barcode_series_from_obs(obs: pd.DataFrame) -> pd.Series:
    """Return one barcode string suitable for the sample-aware TCR merge."""
    barcode_like = barcode_series_from_obs(obs)
    full_barcode = barcode_like.str.extract(r"([ACGTN]+-\d+)", expand=False).fillna("")
    barcode_core = barcode_like.str.extract(r"([ACGTN]+)", expand=False).fillna("")
    return full_barcode.where(full_barcode != "", barcode_core)


def barcode_core_from_obs(obs: pd.DataFrame) -> pd.Series:
    """Extract the 10x core barcode from the current H5AD obs."""
    barcode_like = barcode_series_from_obs(obs)
    full_barcode = barcode_like.map(extract_10x_barcode)
    fallback = barcode_like.str.extract(r"([ACGTN]+)", expand=False).fillna("")
    full_barcode = full_barcode.where(full_barcode != "", fallback)
    return full_barcode.map(normalize_barcode_core)


def raw_sample_key_from_name(gse_id: str, filename: str) -> str:
    """Normalize one raw TCR filename into the same sample key used by obs."""
    name = Path(filename).name
    stem = re.sub(r"\.(csv|txt)(\.gz)?$", "", name, flags=re.IGNORECASE)
    stem = re.sub(r"\.tar(\.gz)?$", "", stem, flags=re.IGNORECASE)

    if gse_id == "GSE168859":
        match = re.search(r"_S(\d+)_filtered_contig_annotations$", stem, flags=re.IGNORECASE)
        return normalize_key(f"S{match.group(1)}") if match else ""

    if gse_id == "GSE188620":
        return normalize_key(re.sub(r"^GSM\d+_|_filtered_contig_annotations$", "", stem, flags=re.IGNORECASE))

    if gse_id == "GSE227709":
        return normalize_key(re.sub(r"^GSM\d+_|_all_contig_annotations$", "", stem, flags=re.IGNORECASE))

    if gse_id == "GSE243572":
        return normalize_key(re.sub(r"^GSM\d+_|_TCR_filtered_contig_annotations$", "", stem, flags=re.IGNORECASE))

    if gse_id == "GSE243905":
        return normalize_key(re.sub(r"^GSM\d+_|_TCR$", "", stem, flags=re.IGNORECASE))

    if gse_id == "GSE212217":
        return normalize_key(re.sub(r"^GSM\d+_|_scTCR_filtered_contig_annotations$", "", stem, flags=re.IGNORECASE))

    if gse_id == "GSE311112":
        match = re.search(r"(?:cll_|pt_)?(\d+_(?:baseline|3yrs|5yrs))", stem, flags=re.IGNORECASE)
        return normalize_key(match.group(1)) if match else ""

    if gse_id == "GSE308075":
        return normalize_key(re.sub(r"^GSM\d+_|_T_filtered_contig_annotations$", "", stem, flags=re.IGNORECASE))

    raise ValueError(f"Unsupported raw filename rule for {gse_id}")


def ensure_gse161918_quick_repair_feasible(obs: pd.DataFrame, raw_df: pd.DataFrame) -> None:
    """Fail fast if the current GSE161918 H5AD lacks the lane identity needed for a safe join."""
    obs_sample_key = sample_key_from_obs("GSE161918", obs)
    if (obs_sample_key == "").any():
        raise RuntimeError("GSE161918 obs sample labels do not encode a recoverable lane key.")
    if raw_df.empty:
        raise RuntimeError("GSE161918 raw TCR table is empty after parsing.")
    if raw_df["sample_id"].eq("").any():
        raise RuntimeError("GSE161918 raw TCR rows do not encode a recoverable lane key.")


def load_gse161918_raw_tcr(gse_dir: Path) -> RawTcrBundle:
    """Parse the direct B1/B2/B3 raw contig files for GSE161918."""
    frames: list[pd.DataFrame] = []
    source_files: list[str] = []
    for path in sorted(gse_dir.glob("GSE161918_B*TCR_filtered_contig_annotations.csv.gz")):
        match = re.search(r"_B(\d+)TCR_", path.name, flags=re.IGNORECASE)
        if match is None:
            continue
        batch_key = f"b{int(match.group(1))}"
        df = filter_productive(pd.read_csv(path, low_memory=False, usecols=RAW_CONTIG_USECOLS))
        lane = df["barcode"].astype(str).str.extract(r"-(\d+)$", expand=False)
        sample_keys = lane.map(lambda value: normalize_key(f"{batch_key}_10xlane{int(value)}"))
        frames.append(standardize_contig_table(df, sample_keys))
        source_files.append(str(path.relative_to(PROJECT_ROOT)))
    tcr_df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["sample_id", "barcode_core", *TCR_COLS])
    return RawTcrBundle("GSE161918", "direct_contig_csv", source_files, tcr_df)


def load_tarred_raw_tcr(
    gse_id: str,
    tar_path: Path,
    member_filter: Callable[[str], bool],
) -> RawTcrBundle:
    """Parse one tar archive that directly contains contig CSV members."""
    frames: list[pd.DataFrame] = []
    source_files: list[str] = []
    with tarfile.open(tar_path, "r:*") as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue
            member_name = Path(member.name).name
            if not member_filter(member_name):
                continue
            handle = tar.extractfile(member)
            if handle is None:
                continue
            df = filter_productive(read_csv_from_handle(handle, member_name))
            sample_key = raw_sample_key_from_name(gse_id, member_name)
            frames.append(standardize_contig_table(df, sample_key))
            source_files.append(f"{tar_path.relative_to(PROJECT_ROOT)}::{member_name}")
    tcr_df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["sample_id", "barcode_core", *TCR_COLS])
    return RawTcrBundle(gse_id, "tarred_contig_csv", source_files, tcr_df)


def load_gse243905_raw_tcr(gse_dir: Path) -> RawTcrBundle:
    """Parse nested per-sample TCR tarballs from the main raw archive."""
    tar_path = gse_dir / "GSE243905_RAW.tar"
    frames: list[pd.DataFrame] = []
    source_files: list[str] = []
    with tarfile.open(tar_path, "r:*") as outer_tar:
        for member in outer_tar.getmembers():
            if not member.isfile():
                continue
            outer_name = Path(member.name).name
            if not outer_name.endswith("_TCR.tar.gz"):
                continue
            sample_key = raw_sample_key_from_name("GSE243905", outer_name)
            handle = outer_tar.extractfile(member)
            if handle is None:
                continue
            raw_bytes = handle.read()
            with tarfile.open(fileobj=io.BytesIO(raw_bytes), mode="r:gz") as inner_tar:
                for inner_member in inner_tar.getmembers():
                    if not inner_member.isfile():
                        continue
                    inner_name = Path(inner_member.name).name
                    if "filtered_contig_annotations.csv" not in inner_name:
                        continue
                    inner_handle = inner_tar.extractfile(inner_member)
                    if inner_handle is None:
                        continue
                    df = filter_productive(read_csv_from_handle(inner_handle, inner_name))
                    frames.append(standardize_contig_table(df, sample_key))
                    source_files.append(
                        f"{tar_path.relative_to(PROJECT_ROOT)}::{outer_name}::{inner_name}"
                    )
    tcr_df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["sample_id", "barcode_core", *TCR_COLS])
    return RawTcrBundle("GSE243905", "nested_tcr_tarballs", source_files, tcr_df)


def load_gse254249_raw_tcr(gse_dir: Path) -> RawTcrBundle:
    """Parse the flat merged 10x VDJ table for GSE254249."""
    path = gse_dir / "GSE254249_10X_VDJ.merge.txt.gz"
    df = filter_productive(pd.read_csv(path, low_memory=False, usecols=RAW_CONTIG_USECOLS))
    sample_keys = df["barcode"].astype(str).str.replace(r"^[ACGTN]+-\d+_", "", regex=True).map(normalize_key)
    tcr_df = standardize_contig_table(df, sample_keys)
    return RawTcrBundle(
        "GSE254249",
        "flat_10x_vdj_table",
        [str(path.relative_to(PROJECT_ROOT))],
        tcr_df,
    )


def load_raw_tcr_bundle(gse_id: str) -> RawTcrBundle:
    """Discover and parse the raw TCR source for one flagged dataset."""
    gse_dir = PROJECT_ROOT / "downloads" / gse_id / "suppl"
    if not gse_dir.exists():
        raise FileNotFoundError(f"Raw supplementary directory not found for {gse_id}: {gse_dir}")

    if gse_id == "GSE161918":
        return load_gse161918_raw_tcr(gse_dir)
    if gse_id == "GSE243905":
        return load_gse243905_raw_tcr(gse_dir)
    if gse_id == "GSE254249":
        return load_gse254249_raw_tcr(gse_dir)

    tar_path = next(gse_dir.glob(f"{gse_id}_RAW.tar"))
    member_filters: dict[str, Callable[[str], bool]] = {
        "GSE168859": lambda name: name.endswith("_filtered_contig_annotations.csv.gz"),
        "GSE188620": lambda name: name.endswith("_filtered_contig_annotations.csv.gz"),
        "GSE227709": lambda name: name.endswith("_all_contig_annotations.csv.gz"),
        "GSE243572": lambda name: name.endswith("_TCR_filtered_contig_annotations.csv.gz"),
        "GSE212217": lambda name: name.endswith("_scTCR_filtered_contig_annotations.csv.gz"),
        "GSE311112": lambda name: name.endswith("_tcr_all_contig_annotations.csv.gz"),
        "GSE308075": lambda name: name.endswith("_T_filtered_contig_annotations.csv.gz"),
    }
    if gse_id not in member_filters:
        raise ValueError(f"No raw-source parser registered for {gse_id}")
    return load_tarred_raw_tcr(gse_id, tar_path, member_filters[gse_id])


def compute_tcr_flags(obs: pd.DataFrame) -> pd.DataFrame:
    """Return boolean TRA/TRB presence flags from one obs frame."""
    tra = normalize_string_series(obs.get("TRA_cdr3", pd.Series("", index=obs.index))).ne("")
    trb = normalize_string_series(obs.get("TRB_cdr3", pd.Series("", index=obs.index))).ne("")
    return pd.DataFrame(
        {
            "has_tra_cdr3": tra.to_numpy(),
            "has_trb_cdr3": trb.to_numpy(),
            "has_any_tcr": (tra | trb).to_numpy(),
            "has_paired_trab": (tra & trb).to_numpy(),
        },
        index=obs.index,
    )


def clear_old_tcr_columns(obs: pd.DataFrame) -> pd.DataFrame:
    """Drop the current TCR fields before rebuilding them from raw source."""
    drop_cols = [column for column in ["TCRseq", *TCR_COLS] if column in obs.columns]
    if drop_cols:
        return obs.drop(columns=drop_cols)
    return obs


def clear_old_phase4_columns(obs: pd.DataFrame) -> pd.DataFrame:
    """Drop previously computed per-GSE Phase 4 score columns."""
    drop_cols = [column for column in PHASE4_REPAIR_COLUMNS if column in obs.columns]
    if drop_cols:
        return obs.drop(columns=drop_cols)
    return obs


def maybe_backup_h5ad(path: Path) -> None:
    """Create one pre-repair backup before the first overwrite."""
    backup = path.with_name(path.name + BACKUP_SUFFIX)
    if path.exists() and not backup.exists():
        shutil.copy2(path, backup)


def compute_phase4_scores_on_adata(adata: ad.AnnData) -> tuple[dict[str, np.ndarray], dict[str, list[str]], dict[str, list[str]]]:
    """Compute Phase 4-style module scores on a temporary normalized/log1p copy."""
    modules = find_module_genes(pd.Index(adata.var_names.astype(str)))
    x = adata.X.copy()
    x = x.tocsr(copy=False) if sparse.issparse(x) else sparse.csr_matrix(x)
    x = x.astype(np.float32, copy=False)
    row_sums = np.asarray(x.sum(axis=1)).ravel().astype(np.float32, copy=False)
    scale = np.zeros_like(row_sums, dtype=np.float32)
    nonzero_mask = row_sums > 0
    scale[nonzero_mask] = TARGET_SUM / row_sums[nonzero_mask]
    sparsefuncs.inplace_row_scale(x, scale)
    np.log1p(x.data, out=x.data)

    var_names = pd.Index(adata.var_names.astype(str))
    gene_means = np.asarray(x.mean(axis=0)).ravel()

    scores: dict[str, np.ndarray] = {}
    control_sets: dict[str, list[str]] = {}
    module_sets = {name: [str(gene) for gene in genes] for name, genes in modules.items()}
    for module_name in ["tra", "trb", "trab", "trd"]:
        gene_list = modules[module_name]
        control_genes = pick_control_genes(
            gene_list,
            var_names,
            gene_means,
            random_state=RANDOM_STATE,
        )
        gene_idx = var_names.get_indexer(gene_list)
        ctrl_idx = var_names.get_indexer(control_genes)
        scores[module_name] = score_chunk(x, gene_idx, ctrl_idx)
        control_sets[module_name] = [str(gene) for gene in control_genes]

    scores["trd_minus_trab"] = (scores["trd"] - scores["trab"]).astype(np.float32, copy=False)
    scores, _ = add_scaled_scores(scores)
    return scores, module_sets, control_sets


def attach_phase4_scores(adata: ad.AnnData, scores: dict[str, np.ndarray]) -> None:
    """Write Phase 4 score columns into one AnnData obs."""
    score_payload = {
        **{column_name: scores[module_name] for module_name, column_name in PHASE4_SCORE_COLUMNS.items()},
        **{column_name: scores[module_name] for module_name, column_name in PHASE4_SCALED_SCORE_COLUMNS.items()},
    }
    for column, values in score_payload.items():
        adata.obs[column] = values
    adata.obs["phase4_trab_minus_trd"] = adata.obs["phase4_trab_score"] - adata.obs["phase4_trd_score"]
    adata.obs["phase4_trab_minus_trd_scaled"] = (
        adata.obs["phase4_trab_score_scaled"] - adata.obs["phase4_trd_score_scaled"]
    )


def downsample_obs_index(obs_index: pd.Index, n_cells: int) -> pd.Index:
    """Sample up to ``n_cells`` rows for plotting."""
    if len(obs_index) <= n_cells:
        return obs_index
    rng = np.random.default_rng(RANDOM_STATE)
    positions = np.sort(rng.choice(len(obs_index), size=n_cells, replace=False))
    return obs_index[positions]


def write_paired_trab_scatter(gse_id: str, obs: pd.DataFrame, sample_size: int) -> Path:
    """Write a Phase 4-style TRAB-vs-TRD scatter plot colored by paired TRA/TRB."""
    plot_index = downsample_obs_index(obs.index, sample_size)
    plot_df = obs.loc[
        plot_index,
        [
            "phase4_trab_score",
            "phase4_trd_score",
            "phase4_trab_score_scaled",
            "phase4_trd_score_scaled",
            "TRA_cdr3",
            "TRB_cdr3",
        ],
    ].copy()
    plot_df["has_paired_trab"] = (
        normalize_string_series(plot_df["TRA_cdr3"]).ne("")
        & normalize_string_series(plot_df["TRB_cdr3"]).ne("")
    )
    plot_df["paired_trab_group"] = np.where(
        plot_df["has_paired_trab"],
        "With paired TRA/TRB",
        "Without paired TRA/TRB",
    )

    palette = {
        "Without paired TRA/TRB": "#D1495B",
        "With paired TRA/TRB": "#0077B6",
    }
    fig, axes = plt.subplots(1, 2, figsize=(12, 6), constrained_layout=True)
    specs = [
        ("phase4_trab_score", "phase4_trd_score", "Raw scores"),
        ("phase4_trab_score_scaled", "phase4_trd_score_scaled", "Scaled scores"),
    ]
    for ax, (x_col, y_col, title) in zip(axes, specs):
        sns.scatterplot(
            data=plot_df,
            x=x_col,
            y=y_col,
            hue="paired_trab_group",
            hue_order=["Without paired TRA/TRB", "With paired TRA/TRB"],
            palette=palette,
            s=6,
            linewidth=0,
            alpha=0.8,
            ax=ax,
        )
        ax.set_title(f"{gse_id} {title}")
        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)
        legend = ax.get_legend()
        if legend is not None:
            legend.set_title("TRA/TRB pairing")
    out_path = FIGURE_DIR / f"{gse_id}_trab_vs_trd_paired_trab.png"
    fig.savefig(out_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    return out_path


def compute_validation_metrics(gse_id: str, obs: pd.DataFrame) -> tuple[dict[str, object], Path]:
    """Validate the rebuilt TCRs against the requested TRD/TRAB score patterns."""
    work = obs[
        [
            "phase4_trd_score",
            "phase4_trab_score",
            "TRA_cdr3",
            "TRB_cdr3",
        ]
    ].copy()
    work["has_paired_trab"] = (
        normalize_string_series(work["TRA_cdr3"]).ne("")
        & normalize_string_series(work["TRB_cdr3"]).ne("")
    )
    work["trd_bin"] = np.select(
        [
            work["phase4_trd_score"] > 0.3,
            (work["phase4_trd_score"] >= 0.1) & (work["phase4_trd_score"] <= 0.3),
        ],
        ["trd_gt_0_3", "trd_0_1_to_0_3"],
        default="trd_lt_0_1",
    )
    bin_summary = (
        work.groupby("trd_bin", dropna=False)
        .agg(
            n_cells=("trd_bin", "size"),
            paired_trab_cells=("has_paired_trab", "sum"),
            paired_trab_fraction=("has_paired_trab", "mean"),
            trab_score_mean=("phase4_trab_score", "mean"),
        )
        .reset_index()
    )
    bin_path = TABLE_DIR / f"{gse_id}_validation_bins.csv"
    bin_summary.to_csv(bin_path, index=False)

    bin_map = bin_summary.set_index("trd_bin")["paired_trab_fraction"].to_dict()
    high_frac = float(bin_map.get("trd_gt_0_3", np.nan))
    mid_frac = float(bin_map.get("trd_0_1_to_0_3", np.nan))
    low_frac = float(bin_map.get("trd_lt_0_1", np.nan))
    paired_trab_mean = float(work.loc[work["has_paired_trab"], "phase4_trab_score"].mean())
    unpaired_trab_mean = float(work.loc[~work["has_paired_trab"], "phase4_trab_score"].mean())

    paired_order_pass = (
        np.isfinite(high_frac)
        and np.isfinite(mid_frac)
        and np.isfinite(low_frac)
        and (high_frac <= mid_frac <= low_frac)
        and (high_frac < low_frac)
    )
    trab_shift_pass = (
        np.isfinite(paired_trab_mean)
        and np.isfinite(unpaired_trab_mean)
        and (paired_trab_mean > unpaired_trab_mean)
    )
    validation_pass = bool(paired_order_pass and trab_shift_pass)

    return (
        {
            "validation_pass": validation_pass,
            "validation_note": (
                "paired coverage ordered by TRD bins and paired TRA/TRB cells shift upward on TRAB"
                if validation_pass
                else "pattern check failed"
            ),
            "trd_gt_0_3_paired_fraction": high_frac,
            "trd_0_1_to_0_3_paired_fraction": mid_frac,
            "trd_lt_0_1_paired_fraction": low_frac,
            "paired_trab_score_mean": paired_trab_mean,
            "unpaired_trab_score_mean": unpaired_trab_mean,
        },
        bin_path,
    )


def summarize_repair(
    gse_id: str,
    h5ad_path: Path,
    source_bundle: RawTcrBundle,
    old_flags: pd.DataFrame,
    new_flags: pd.DataFrame,
    changed_pair_cells: int,
    duplicated_raw_key_rows: int,
    score_status: str,
    score_note: str,
    validation_metrics: dict[str, object],
    validation_table: Path | None,
    figure_path: Path | None,
) -> dict[str, object]:
    """Build one summary row for the QC table."""
    aggregated = aggregate_tcr_rows(source_bundle.tcr_df)
    return {
        "gse_id": gse_id,
        "h5ad_path": str(h5ad_path),
        "raw_source_label": source_bundle.source_label,
        "raw_source_file_n": len(source_bundle.source_files),
        "raw_rows": int(len(source_bundle.tcr_df)),
        "raw_unique_keys": int(len(aggregated)),
        "raw_duplicated_key_rows": int(duplicated_raw_key_rows),
        "cells_total": int(len(new_flags)),
        "old_any_tcr_cells": int(old_flags["has_any_tcr"].sum()),
        "old_paired_trab_cells": int(old_flags["has_paired_trab"].sum()),
        "new_any_tcr_cells": int(new_flags["has_any_tcr"].sum()),
        "new_paired_trab_cells": int(new_flags["has_paired_trab"].sum()),
        "old_any_tcr_fraction": float(old_flags["has_any_tcr"].mean()),
        "old_paired_trab_fraction": float(old_flags["has_paired_trab"].mean()),
        "new_any_tcr_fraction": float(new_flags["has_any_tcr"].mean()),
        "new_paired_trab_fraction": float(new_flags["has_paired_trab"].mean()),
        "changed_pair_cells": int(changed_pair_cells),
        "score_status": score_status,
        "score_note": score_note,
        "validation_pass": validation_metrics.get("validation_pass", False),
        "validation_note": validation_metrics.get("validation_note", ""),
        "trd_gt_0_3_paired_fraction": validation_metrics.get("trd_gt_0_3_paired_fraction", np.nan),
        "trd_0_1_to_0_3_paired_fraction": validation_metrics.get("trd_0_1_to_0_3_paired_fraction", np.nan),
        "trd_lt_0_1_paired_fraction": validation_metrics.get("trd_lt_0_1_paired_fraction", np.nan),
        "paired_trab_score_mean": validation_metrics.get("paired_trab_score_mean", np.nan),
        "unpaired_trab_score_mean": validation_metrics.get("unpaired_trab_score_mean", np.nan),
        "validation_table": str(validation_table.relative_to(PROJECT_ROOT)) if validation_table else "",
        "figure_path": str(figure_path.relative_to(PROJECT_ROOT)) if figure_path else "",
    }


def write_summary_markdown(summary_df: pd.DataFrame) -> None:
    """Write a compact markdown QC summary."""
    detail_columns = [
        "gse_id",
        "old_paired_trab_fraction",
        "new_paired_trab_fraction",
        "trd_gt_0_3_paired_fraction",
        "trd_lt_0_1_paired_fraction",
        "paired_trab_score_mean",
        "unpaired_trab_score_mean",
        "validation_pass",
        "changed_pair_cells",
        "score_status",
        "validation_note",
    ]
    lines = [
        "# Raw TCR Repair And Phase 4 QC",
        "",
        f"- GSEs processed: {len(summary_df)}",
        f"- GSEs with scoreable Phase 4 plots: {int((summary_df['score_status'] == 'scored').sum())}",
        f"- GSEs passing the requested TRD/TRAB validation pattern: {int(summary_df['validation_pass'].fillna(False).sum())}",
        f"- GSEs skipped at scoring due to missing module genes or blocking issues: {int((summary_df['score_status'] != 'scored').sum())}",
        "",
        "## Per-GSE summary",
        "",
    ]
    lines.append("| " + " | ".join(detail_columns) + " |")
    lines.append("|" + "|".join(["---"] * len(detail_columns)) + "|")
    for row in summary_df[detail_columns].itertuples(index=False):
        values = [clean_text(value) for value in row]
        lines.append("| " + " | ".join(values) + " |")
    lines.append("")
    SUMMARY_MD.write_text("\n".join(lines), encoding="utf-8")


def process_one_gse(
    gse_id: str,
    h5ad_path: Path,
    *,
    write_h5ad: bool,
    scatter_sample_size: int,
) -> dict[str, object]:
    """Repair one flagged GSE and return the QC summary row."""
    logging.info("Processing %s from %s", gse_id, h5ad_path)
    source_bundle = load_raw_tcr_bundle(gse_id)
    logging.info("%s raw source parsed: rows=%s files=%s", gse_id, len(source_bundle.tcr_df), len(source_bundle.source_files))
    duplicated_raw_key_rows = int(
        source_bundle.tcr_df.duplicated(subset=["sample_id", "barcode_core"], keep=False).sum()
    )

    adata = ad.read_h5ad(h5ad_path)
    obs = adata.obs.copy()
    if gse_id == "GSE161918":
        ensure_gse161918_quick_repair_feasible(obs, source_bundle.tcr_df)

    old_flags = compute_tcr_flags(obs)
    old_pair_key = (
        normalize_string_series(obs.get("TRA_cdr3", pd.Series("", index=obs.index)))
        + "||"
        + normalize_string_series(obs.get("TRB_cdr3", pd.Series("", index=obs.index)))
    )

    join_obs = pd.DataFrame(index=obs.index)
    join_obs["sample_id"] = sample_key_from_obs(gse_id, obs)
    join_obs["barcode"] = merge_barcode_series_from_obs(obs)

    rebuilt_tcr = aggregate_tcr_rows(source_bundle.tcr_df)
    merged_tcr = merge_tcr_into_obs(join_obs, rebuilt_tcr)
    logging.info("%s merge prepared: obs_cells=%s raw_unique_keys=%s", gse_id, adata.n_obs, len(rebuilt_tcr))

    adata.obs = clear_old_phase4_columns(clear_old_tcr_columns(obs))
    for column in [*TCR_COLS, "TCRseq"]:
        adata.obs[column] = merged_tcr[column].reindex(adata.obs.index)

    new_flags = compute_tcr_flags(adata.obs)
    new_pair_key = (
        normalize_string_series(adata.obs["TRA_cdr3"])
        + "||"
        + normalize_string_series(adata.obs["TRB_cdr3"])
    )
    changed_pair_cells = int((old_pair_key != new_pair_key).sum())

    figure_path: Path | None = None
    score_status = "scored"
    score_note = ""
    validation_metrics: dict[str, object] = {
        "validation_pass": False,
        "validation_note": "not_scored",
    }
    validation_table: Path | None = None
    try:
        logging.info("%s Phase 4 scoring started", gse_id)
        scores, module_sets, control_sets = compute_phase4_scores_on_adata(adata)
        attach_phase4_scores(adata, scores)
        adata.uns["raw_tcr_repair"] = {
            "gse_id": gse_id,
            "source_label": source_bundle.source_label,
            "source_files": source_bundle.source_files,
            "join_key": "sample_id + barcode_core",
            "phase4_mode": "temporary_normalize_total_log1p_on_count_space_X",
            "phase4_module_genes": module_sets,
            "phase4_control_genes": control_sets,
        }
        figure_path = write_paired_trab_scatter(gse_id, adata.obs, scatter_sample_size)
        validation_metrics, validation_table = compute_validation_metrics(gse_id, adata.obs)
        logging.info("%s Phase 4 scoring complete", gse_id)
    except Exception as exc:
        score_status = "unscoreable"
        score_note = str(exc)
        validation_metrics = {
            "validation_pass": False,
            "validation_note": f"unscoreable: {exc}",
        }
        logging.warning("Skipping Phase 4 plot for %s: %s", gse_id, exc)

    if write_h5ad:
        maybe_backup_h5ad(h5ad_path)
        write_h5ad_atomic(adata, h5ad_path)
        logging.info("%s H5AD updated in place", gse_id)

    return summarize_repair(
        gse_id=gse_id,
        h5ad_path=h5ad_path,
        source_bundle=source_bundle,
        old_flags=old_flags,
        new_flags=new_flags,
        changed_pair_cells=changed_pair_cells,
        duplicated_raw_key_rows=duplicated_raw_key_rows,
        score_status=score_status,
        score_note=score_note,
        validation_metrics=validation_metrics,
        validation_table=validation_table,
        figure_path=figure_path,
    )


def main() -> None:
    """Run the flagged raw-TCR repair workflow."""
    args = parse_args()
    ensure_dirs()
    configure_logging()
    registry_map = load_registry_map(args.registry)

    summary_rows: list[dict[str, object]] = []
    for gse_id in args.gse:
        if gse_id not in registry_map:
            raise KeyError(f"{gse_id} is not present in {args.registry}")
        summary_rows.append(
            process_one_gse(
                gse_id,
                registry_map[gse_id],
                write_h5ad=args.write_h5ad,
                scatter_sample_size=args.scatter_sample_size,
            )
        )

    summary_df = pd.DataFrame(summary_rows).sort_values("gse_id").reset_index(drop=True)
    summary_df.to_csv(SUMMARY_CSV, index=False)
    write_summary_markdown(summary_df)
    logging.info("Wrote summary table to %s", SUMMARY_CSV)
    logging.info("Wrote QC markdown to %s", SUMMARY_MD)


if __name__ == "__main__":
    main()
