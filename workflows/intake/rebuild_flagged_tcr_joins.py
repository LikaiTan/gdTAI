#!/usr/bin/env python3
"""Rebuild flagged TCR joins as validated, sample-aware metadata sidecars.

The workflow is deliberately read-only with respect to every H5AD. Productive
raw VDJ contigs are joined to RNA metadata only by ``sample_id + barcode_core``.
All four TCR chains are retained, including per-chain UMI/read support and
source-file provenance. Sources or samples with ambiguous identity fail closed.
"""

from __future__ import annotations

# TNK_WORKFLOW_BOOTSTRAP
import sys as _tnk_sys
from pathlib import Path as _TnkPath

_TNK_PROJECT_ROOT = _TnkPath(__file__).resolve().parents[2]
for _tnk_path in (_TNK_PROJECT_ROOT, _TNK_PROJECT_ROOT / "src"):
    _tnk_value = str(_tnk_path)
    if _tnk_value not in _tnk_sys.path:
        _tnk_sys.path.insert(0, _tnk_value)

import argparse
import gzip
import hashlib
import io
import json
import logging
import re
import tarfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import BinaryIO

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SOURCE_AUDIT = (
    PROJECT_ROOT
    / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_tcr_conflict_audit/source_join_integrity_audit.csv"
)
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
TABLE_DIR = OUTPUT_ROOT / "tables/tcr_join_rebuild"
PER_SOURCE_DIR = TABLE_DIR / "per_source"
LOG_DIR = OUTPUT_ROOT / "logs/tcr_join_rebuild"
SUMMARY_CSV = TABLE_DIR / "source_rebuild_summary.csv"
SAMPLE_AUDIT_CSV = TABLE_DIR / "sample_mapping_audit.csv"
REUSE_AUDIT_CSV = TABLE_DIR / "receptor_reuse_audit.csv"
SUMMARY_MD = LOG_DIR / "tcr_join_rebuild_summary.md"
RUN_MANIFEST = LOG_DIR / "run_manifest.json"
LOG_PATH = LOG_DIR / "tcr_join_rebuild.log"
STAGED_REPLACEMENT = TABLE_DIR / "validated_tcr_replacement_sidecar.parquet"
STAGED_REPLACEMENT_MANIFEST = TABLE_DIR / "validated_tcr_replacement_manifest.csv"

CHAINS = ("TRA", "TRB", "TRG", "TRD")
CHAIN_FIELDS = (
    "cdr3",
    "v",
    "d",
    "j",
    "cdr3_nt",
    "clone_id",
    "umis",
    "reads",
    "umi_available",
    "read_available",
    "n_productive_contigs",
    "source_file",
    "selected_contig_id",
)
STRING_FIELDS = (
    "cdr3",
    "v",
    "d",
    "j",
    "cdr3_nt",
    "clone_id",
    "source_file",
    "selected_contig_id",
)
RAW_COLUMNS = {
    "barcode",
    "is_cell",
    "contig_id",
    "high_confidence",
    "chain",
    "v_gene",
    "d_gene",
    "j_gene",
    "productive",
    "full_length",
    "cdr3",
    "cdr3_nt",
    "reads",
    "umis",
    "raw_clonotype_id",
    "clonotype_id",
    "clone_id",
}
DEFAULT_GSES = (
    "GSE125527",
    "GSE171037",
    "GSE178882",
    "GSE188620",
    "GSE211504",
    "GSE212217",
    "GSE228597",
    "GSE235863",
    "GSE243572",
    "GSE243905",
    "GSE254176",
    "GSE254249",
    "GSE287541",
    "GSE311112",
)
UNRESOLVED_DONOR_SOURCES = {
    "GSE178882": "paired libraries; donor grouping is not retained in RNA metadata",
    "GSE243572": "pooled libraries; donor identity is not cell-resolved",
}


@dataclass(frozen=True)
class SourceSpec:
    """Raw-input and sample-key contract for one source."""

    loader: str
    member_pattern: str = ""
    source_reason: str = ""


@dataclass
class RawBundle:
    """Long-form productive contigs and source-level parse accounting."""

    contigs: pd.DataFrame
    source_files: list[str] = field(default_factory=list)
    n_input_rows: int = 0
    n_productive_rows: int = 0


SOURCE_SPECS: dict[str, SourceSpec] = {
    "GSE125527": SourceSpec(
        "quarantine",
        source_reason="RNA metadata does not retain the raw library/patient identity required for a safe join",
    ),
    "GSE171037": SourceSpec("direct_tar", r"_filtered_contig_annotations\.csv\.gz$"),
    "GSE178882": SourceSpec("direct_tar", r"_filtered_contig_annotations\.csv\.gz$"),
    "GSE188620": SourceSpec("direct_tar", r"_filtered_contig_annotations\.csv\.gz$"),
    "GSE211504": SourceSpec("direct_tar", r"filtered_contig_annotations_.*\.csv\.gz$"),
    "GSE212217": SourceSpec("direct_tar", r"_scTCR_filtered_contig_annotations\.csv\.gz$"),
    "GSE228597": SourceSpec("direct_tar", r"_all_contig_annotations_tcr\.csv\.gz$"),
    "GSE235863": SourceSpec("wide_gse235863"),
    "GSE243572": SourceSpec("direct_tar", r"_TCR_filtered_contig_annotations\.csv\.gz$"),
    "GSE243905": SourceSpec("nested_gse243905"),
    "GSE254176": SourceSpec("direct_tar", r"-TCR_filtered_contig_annotations\.csv\.gz$"),
    "GSE254249": SourceSpec("flat_gse254249"),
    "GSE287541": SourceSpec(
        "quarantine",
        source_reason="No raw VDJ contig or receptor table is present in the project",
    ),
    "GSE311112": SourceSpec("direct_tar", r"_tcr_all_contig_annotations\.csv\.gz$"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gse", nargs="+", default=list(DEFAULT_GSES))
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing sidecars. H5AD files are never written.",
    )
    return parser.parse_args()


def ensure_dirs() -> None:
    for path in (TABLE_DIR, PER_SOURCE_DIR, LOG_DIR):
        path.mkdir(parents=True, exist_ok=True)


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(LOG_PATH, mode="w"), logging.StreamHandler()],
        force=True,
    )


def clean_text(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"", "nan", "none", "null", "na", "<na>"} else text


def normalize_key(value: object) -> str:
    text = clean_text(value).lower()
    return re.sub(r"[^a-z0-9]+", "_", text).strip("_")


def barcode_core(value: object) -> str:
    """Extract a 10x nucleotide barcode without lane suffixes."""
    text = clean_text(value).upper()
    matches = re.findall(r"[ACGTN]{14,}(?=-\d+|$)", text)
    if matches:
        return matches[-1]
    matches = re.findall(r"[ACGTN]{14,}", text)
    return matches[-1] if matches else ""


def truthy(series: pd.Series) -> pd.Series:
    values = series.astype("string").str.strip().str.lower()
    return values.isin({"true", "t", "1", "yes", "y"}).fillna(False)


def obs_barcode_series(obs: pd.DataFrame) -> pd.Series:
    for column in ("barcode", "cell_barcode"):
        if column not in obs.columns:
            continue
        values = obs[column].map(barcode_core)
        if values.ne("").mean() >= 0.5:
            return obs[column]
    return pd.Series(obs.index.astype(str), index=obs.index, dtype=object)


def rna_sample_keys(gse_id: str, obs: pd.DataFrame) -> pd.Series:
    """Derive a documented canonical RNA sample key for one source."""
    if gse_id == "GSE254249":
        raw = obs_barcode_series(obs).map(clean_text)
        return raw.str.replace(r"^.*?[ACGTN]{14,}-\d+_", "", regex=True).map(normalize_key)
    if gse_id == "GSE235863":
        return obs["sample"].map(normalize_key)
    if gse_id == "GSE228597":
        return obs["sample_id"].map(normalize_key)

    source_column = "sample"
    if source_column not in obs.columns:
        raise KeyError(f"{gse_id} lacks required RNA sample column `{source_column}`")
    sample = obs[source_column].map(clean_text).str.lower()
    sample = sample.str.replace(r"^gsm\d+_", "", regex=True)
    if gse_id == "GSE171037":
        return sample.map(normalize_key)
    if gse_id == "GSE178882":
        return sample.map(normalize_key)
    if gse_id == "GSE188620":
        return sample.str.replace(r"_matrix$", "", regex=True).map(normalize_key)
    if gse_id == "GSE211504":
        return sample.str.extract(r"([0-9]+\.[0-9]+)$", expand=False).map(normalize_key)
    if gse_id == "GSE212217":
        return sample.str.replace(r"_scrna$", "", regex=True).map(normalize_key)
    if gse_id == "GSE243572":
        return sample.str.replace(r"_gex$", "", regex=True).map(normalize_key)
    if gse_id == "GSE243905":
        return sample.str.replace(r"_gex$", "", regex=True).map(normalize_key)
    if gse_id == "GSE254176":
        return sample.map(normalize_key)
    if gse_id == "GSE311112":
        extracted = sample.str.extract(
            r"(?:pt_|cll_)?([0-9]+_(?:baseline|3yrs|5yrs|relapse|btk_clone))",
            expand=False,
        )
        return extracted.map(normalize_key)
    raise ValueError(f"No RNA sample-key rule for {gse_id}")


def rna_donor_keys(gse_id: str, obs: pd.DataFrame, sample_keys: pd.Series) -> pd.Series:
    """Return the best available donor key without weakening the join key."""
    if gse_id == "GSE178882":
        donor = sample_keys.str.extract(r"^([0-9]+)(?:_[0-9]+)$", expand=False)
        donor = donor.fillna(sample_keys)
    elif gse_id == "GSE188620":
        donor = sample_keys.str.extract(r"^([0-9]+_[0-9]+)", expand=False)
    elif gse_id == "GSE212217":
        donor = sample_keys.str.extract(r"^(pem[0-9]+)", expand=False)
    elif gse_id == "GSE235863":
        donor = sample_keys.str.extract(r"^(p[0-9]+)", expand=False)
    elif gse_id == "GSE243905":
        donor = sample_keys.str.extract(r"^([cp][0-9]+)", expand=False)
    elif gse_id == "GSE254249":
        donor = sample_keys.str.extract(r"^(crc[0-9]+)", expand=False)
    elif gse_id in {"GSE211504", "GSE311112"}:
        donor = sample_keys.str.extract(r"^([0-9]+)", expand=False)
    elif gse_id == "GSE254176":
        donor = sample_keys.str.extract(r"^(om_(?:hl|hs)_[0-9]+)", expand=False)
    else:
        donor = pd.Series("", index=obs.index, dtype=object)
        for column in ("donor_id", "patient", "PatientID", "patient_assignment", "donor"):
            if column not in obs.columns:
                continue
            candidate = obs[column].map(normalize_key)
            if candidate.ne("").any():
                donor = candidate
                break
    donor = pd.Series(donor.astype(object), index=obs.index).map(normalize_key)
    return donor.where(donor.ne(""), sample_keys)


def raw_sample_key_from_name(gse_id: str, filename: str) -> str:
    """Derive the same sample key from a raw VDJ member name."""
    stem = Path(filename).name
    stem = re.sub(r"\.(csv|tsv)(\.gz)?$", "", stem, flags=re.I)
    stem = re.sub(r"\.tar(\.gz)?$", "", stem, flags=re.I)
    stem = re.sub(r"^GSM\d+_", "", stem, flags=re.I)
    if gse_id == "GSE171037":
        stem = re.sub(r"_filtered_contig_annotations$", "", stem, flags=re.I)
    elif gse_id in {"GSE178882", "GSE188620"}:
        stem = re.sub(r"_filtered_contig_annotations$", "", stem, flags=re.I)
    elif gse_id == "GSE211504":
        match = re.search(r"filtered_contig_annotations_([0-9]+\.[0-9]+)$", stem, flags=re.I)
        stem = match.group(1) if match else ""
    elif gse_id == "GSE212217":
        stem = re.sub(r"_scTCR_filtered_contig_annotations$", "", stem, flags=re.I)
    elif gse_id == "GSE228597":
        stem = re.sub(r"_all_contig_annotations_tcr$", "", stem, flags=re.I)
        stem = re.sub(r"_\d+_heart$", "", stem, flags=re.I)
    elif gse_id == "GSE243572":
        stem = re.sub(r"_TCR_filtered_contig_annotations$", "", stem, flags=re.I)
    elif gse_id == "GSE243905":
        stem = re.sub(r"_TCR$", "", stem, flags=re.I)
    elif gse_id == "GSE254176":
        stem = re.sub(r"-TCR_filtered_contig_annotations$", "", stem, flags=re.I)
    elif gse_id == "GSE311112":
        match = re.search(
            r"(?:cll_|pt_)?([0-9]+_(?:baseline|3yrs|5yrs|relapse|btk_clone))",
            stem,
            flags=re.I,
        )
        stem = match.group(1) if match else ""
    else:
        raise ValueError(f"No raw filename sample-key rule for {gse_id}")
    return normalize_key(stem)


def read_contig_csv(handle: BinaryIO, name: str) -> pd.DataFrame:
    usecols = lambda column: column in RAW_COLUMNS
    if name.lower().endswith(".gz"):
        with gzip.open(handle, "rt", encoding="utf-8", errors="replace") as text:
            return pd.read_csv(text, low_memory=False, usecols=usecols)
    return pd.read_csv(handle, low_memory=False, usecols=usecols)


def standardize_contigs(
    df: pd.DataFrame,
    *,
    sample_id: str | pd.Series,
    source_file: str,
) -> tuple[pd.DataFrame, int, int]:
    """Filter and standardize one raw 10x-like contig table."""
    n_input = len(df)
    work = df.copy()
    for flag in ("is_cell", "high_confidence", "productive"):
        if flag in work.columns:
            work = work.loc[truthy(work[flag])].copy()
    n_productive = len(work)
    if "barcode" not in work.columns or "chain" not in work.columns:
        return empty_contig_frame(), n_input, n_productive

    out = pd.DataFrame(index=work.index)
    if isinstance(sample_id, pd.Series):
        out["sample_id"] = sample_id.reindex(work.index).map(normalize_key)
    else:
        out["sample_id"] = normalize_key(sample_id)
    out["barcode_core"] = work["barcode"].map(barcode_core)
    out["chain"] = work["chain"].map(clean_text).str.upper()
    for source, target in (
        ("cdr3", "cdr3"),
        ("v_gene", "v"),
        ("d_gene", "d"),
        ("j_gene", "j"),
        ("cdr3_nt", "cdr3_nt"),
    ):
        out[target] = work[source].map(clean_text) if source in work.columns else ""
    clone_column = next((c for c in ("raw_clonotype_id", "clonotype_id", "clone_id") if c in work.columns), None)
    out["clone_id"] = work[clone_column].map(clean_text) if clone_column else ""
    out["contig_id"] = work["contig_id"].map(clean_text) if "contig_id" in work.columns else ""
    out["umis"] = pd.to_numeric(work["umis"], errors="coerce") if "umis" in work.columns else np.nan
    out["reads"] = pd.to_numeric(work["reads"], errors="coerce") if "reads" in work.columns else np.nan
    out["umi_available"] = "umis" in work.columns
    out["read_available"] = "reads" in work.columns
    out["full_length"] = truthy(work["full_length"]) if "full_length" in work.columns else False
    out["source_file"] = source_file
    out = out.loc[
        out["sample_id"].ne("")
        & out["barcode_core"].ne("")
        & out["chain"].isin(CHAINS)
        & out["cdr3"].ne("")
    ].copy()
    return out, n_input, n_productive


def empty_contig_frame() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "sample_id",
            "barcode_core",
            "chain",
            "cdr3",
            "v",
            "d",
            "j",
            "cdr3_nt",
            "clone_id",
            "contig_id",
            "umis",
            "reads",
            "umi_available",
            "read_available",
            "full_length",
            "source_file",
        ]
    )


def relative_label(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path)


def load_direct_tar(gse_id: str, spec: SourceSpec) -> RawBundle:
    tar_path = PROJECT_ROOT / f"data/datasets/{gse_id}/raw/legacy_source/suppl/{gse_id}_RAW.tar"
    frames: list[pd.DataFrame] = []
    files: list[str] = []
    n_input = 0
    n_productive = 0
    pattern = re.compile(spec.member_pattern, flags=re.I)
    with tarfile.open(tar_path, "r:*") as archive:
        for member in archive.getmembers():
            name = Path(member.name).name
            if not member.isfile() or not pattern.search(name):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                continue
            source_file = f"{relative_label(tar_path)}::{name}"
            df = read_contig_csv(handle, name)
            standardized, n_all, n_prod = standardize_contigs(
                df,
                sample_id=raw_sample_key_from_name(gse_id, name),
                source_file=source_file,
            )
            frames.append(standardized)
            files.append(source_file)
            n_input += n_all
            n_productive += n_prod
    return RawBundle(
        pd.concat(frames, ignore_index=True) if frames else empty_contig_frame(),
        files,
        n_input,
        n_productive,
    )


def load_nested_gse243905() -> RawBundle:
    gse_id = "GSE243905"
    tar_path = PROJECT_ROOT / f"data/datasets/{gse_id}/raw/legacy_source/suppl/{gse_id}_RAW.tar"
    frames: list[pd.DataFrame] = []
    files: list[str] = []
    n_input = 0
    n_productive = 0
    with tarfile.open(tar_path, "r:*") as outer:
        for member in outer.getmembers():
            outer_name = Path(member.name).name
            if not member.isfile() or not outer_name.endswith("_TCR.tar.gz"):
                continue
            handle = outer.extractfile(member)
            if handle is None:
                continue
            with tarfile.open(fileobj=io.BytesIO(handle.read()), mode="r:gz") as inner:
                for inner_member in inner.getmembers():
                    inner_name = Path(inner_member.name).name
                    if not inner_member.isfile() or "filtered_contig_annotations.csv" not in inner_name:
                        continue
                    inner_handle = inner.extractfile(inner_member)
                    if inner_handle is None:
                        continue
                    source_file = f"{relative_label(tar_path)}::{outer_name}::{inner_name}"
                    df = read_contig_csv(inner_handle, inner_name)
                    standardized, n_all, n_prod = standardize_contigs(
                        df,
                        sample_id=raw_sample_key_from_name(gse_id, outer_name),
                        source_file=source_file,
                    )
                    frames.append(standardized)
                    files.append(source_file)
                    n_input += n_all
                    n_productive += n_prod
    return RawBundle(pd.concat(frames, ignore_index=True), files, n_input, n_productive)


def load_flat_gse254249() -> RawBundle:
    path = PROJECT_ROOT / "data/datasets/GSE254249/raw/legacy_source/suppl/GSE254249_10X_VDJ.merge.txt.gz"
    df = pd.read_csv(path, low_memory=False, usecols=lambda column: column in RAW_COLUMNS | {"sample"})
    sample = df["sample"].map(normalize_key)
    standardized, n_all, n_prod = standardize_contigs(
        df,
        sample_id=sample,
        source_file=relative_label(path),
    )
    return RawBundle(standardized, [relative_label(path)], n_all, n_prod)


def load_wide_gse235863() -> RawBundle:
    paths = [
        PROJECT_ROOT / "data/datasets/GSE235863/raw/legacy_source/suppl/GSE235863_five_patients_TCR.csv.gz",
        PROJECT_ROOT / "data/datasets/GSE235863/raw/legacy_source/suppl/GSE235863_nine_patients_TCR.csv.gz",
    ]
    frames: list[pd.DataFrame] = []
    n_input = 0
    for path in paths:
        df = pd.read_csv(path, low_memory=False)
        n_input += len(df)
        for suffix, chain in (("A", "TRA"), ("B", "TRB")):
            out = pd.DataFrame(index=df.index)
            out["sample_id"] = df["sample"].map(normalize_key)
            out["barcode_core"] = df["barcode"].map(barcode_core)
            out["chain"] = chain
            for source, target in (
                (f"cdr3.{suffix}", "cdr3"),
                (f"v_gene.{suffix}", "v"),
                (f"j_gene.{suffix}", "j"),
                (f"cdr3_nt.{suffix}", "cdr3_nt"),
            ):
                out[target] = df[source].map(clean_text) if source in df.columns else ""
            out["d"] = ""
            out["clone_id"] = df["clone.id"].map(clean_text) if "clone.id" in df.columns else ""
            out["contig_id"] = ""
            out["umis"] = np.nan
            out["reads"] = np.nan
            out["umi_available"] = False
            out["read_available"] = False
            out["full_length"] = True
            out["source_file"] = relative_label(path)
            out = out.loc[out["barcode_core"].ne("") & out["cdr3"].ne("")]
            frames.append(out)
    contigs = pd.concat(frames, ignore_index=True) if frames else empty_contig_frame()
    return RawBundle(contigs, [relative_label(p) for p in paths], n_input, len(contigs))


def load_raw_bundle(gse_id: str, spec: SourceSpec) -> RawBundle:
    if spec.loader == "direct_tar":
        return load_direct_tar(gse_id, spec)
    if spec.loader == "nested_gse243905":
        return load_nested_gse243905()
    if spec.loader == "flat_gse254249":
        return load_flat_gse254249()
    if spec.loader == "wide_gse235863":
        return load_wide_gse235863()
    if spec.loader == "quarantine":
        return RawBundle(empty_contig_frame())
    raise ValueError(f"Unknown loader `{spec.loader}` for {gse_id}")


def select_best_contigs(contigs: pd.DataFrame) -> pd.DataFrame:
    """Select one productive contig per sample, barcode, and chain."""
    if contigs.empty:
        return contigs.copy()
    work = contigs.copy()
    work["n_productive_contigs"] = work.groupby(
        ["sample_id", "barcode_core", "chain"], dropna=False
    )["chain"].transform("size")
    work["umi_rank"] = pd.to_numeric(work["umis"], errors="coerce").fillna(-1)
    work["read_rank"] = pd.to_numeric(work["reads"], errors="coerce").fillna(-1)
    work = work.sort_values(
        [
            "sample_id",
            "barcode_core",
            "chain",
            "umi_rank",
            "read_rank",
            "full_length",
            "contig_id",
        ],
        ascending=[True, True, True, False, False, False, True],
        kind="mergesort",
    )
    return work.drop_duplicates(["sample_id", "barcode_core", "chain"], keep="first").drop(
        columns=["umi_rank", "read_rank"]
    )


def pivot_chains(selected: pd.DataFrame) -> pd.DataFrame:
    keys = ["sample_id", "barcode_core"]
    base = selected[keys].drop_duplicates().copy() if not selected.empty else pd.DataFrame(columns=keys)
    for chain in CHAINS:
        sub = selected.loc[selected["chain"].eq(chain)].copy()
        chain_frame = sub[keys].copy()
        for field in CHAIN_FIELDS:
            source = "contig_id" if field == "selected_contig_id" else field
            chain_frame[f"{chain}_{field}"] = sub[source].to_numpy() if source in sub.columns else np.nan
        base = base.merge(chain_frame, on=keys, how="outer", validate="one_to_one")
    return base


def old_tcr_flags(obs: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(index=obs.index)
    for chain in CHAINS:
        column = f"{chain}_cdr3"
        out[f"old_has_{chain.lower()}"] = (
            obs[column].map(clean_text).ne("") if column in obs.columns else False
        )
    out["old_has_any_ab"] = out["old_has_tra"] | out["old_has_trb"]
    out["old_has_any_gd"] = out["old_has_trg"] | out["old_has_trd"]
    out["old_has_paired_ab"] = out["old_has_tra"] & out["old_has_trb"]
    out["old_has_paired_gd"] = out["old_has_trg"] & out["old_has_trd"]
    return out


def source_h5ad_path(gse_id: str) -> Path:
    return PROJECT_ROOT / f"data/datasets/{gse_id}/processed/current.h5ad"


def build_quarantine_sidecar(
    gse_id: str,
    obs: pd.DataFrame,
    reason: str,
) -> pd.DataFrame:
    sidecar = pd.DataFrame(
        {
            "source_gse_id": gse_id,
            "source_obs_name": obs.index.astype(str),
            "join_sample_id": "",
            "join_donor_id": "",
            "barcode_core": obs_barcode_series(obs).map(barcode_core).to_numpy(),
            "join_status": "source_quarantined",
            "join_reason": reason,
            "tcr_assignment_eligible": False,
            "replacement_eligible": False,
        },
        index=obs.index,
    )
    return add_empty_chain_columns(sidecar)


def add_empty_chain_columns(frame: pd.DataFrame) -> pd.DataFrame:
    for chain in CHAINS:
        for field in CHAIN_FIELDS:
            column = f"{chain}_{field}"
            if field in STRING_FIELDS:
                frame[column] = ""
            elif field in {"umi_available", "read_available"}:
                frame[column] = False
            else:
                frame[column] = pd.Series(pd.NA, index=frame.index, dtype="Int64")
    return frame


def build_sidecar(
    gse_id: str,
    obs: pd.DataFrame,
    bundle: RawBundle,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    rna = pd.DataFrame(index=obs.index)
    rna["source_obs_name"] = obs.index.astype(str)
    sample_keys = rna_sample_keys(gse_id, obs)
    rna["sample_id"] = sample_keys.to_numpy()
    rna["donor_id"] = rna_donor_keys(gse_id, obs, sample_keys).to_numpy()
    rna["barcode_core"] = obs_barcode_series(obs).map(barcode_core).to_numpy()
    valid_key = rna["sample_id"].ne("") & rna["barcode_core"].ne("")
    rna["rna_key_multiplicity"] = 0
    rna.loc[valid_key, "rna_key_multiplicity"] = (
        rna.loc[valid_key].groupby(["sample_id", "barcode_core"])["source_obs_name"].transform("size")
    )

    raw = bundle.contigs.copy()
    raw_files_per_sample = (
        raw.groupby("sample_id")["source_file"].nunique().rename("n_raw_files")
        if not raw.empty
        else pd.Series(dtype=int, name="n_raw_files")
    )
    ambiguous_samples = set(raw_files_per_sample[raw_files_per_sample.gt(1)].index)
    eligible_raw = raw.loc[~raw["sample_id"].isin(ambiguous_samples)].copy()
    selected = select_best_contigs(eligible_raw)
    wide = pivot_chains(selected)

    raw_sample_counts = raw.groupby("sample_id").size().rename("n_productive_contigs") if not raw.empty else pd.Series(dtype=int)
    rna_sample_counts = rna.groupby("sample_id").size().rename("n_rna_cells")
    sample_audit = pd.concat([rna_sample_counts, raw_sample_counts, raw_files_per_sample], axis=1).fillna(0).reset_index()
    sample_audit.insert(0, "source_gse_id", gse_id)
    sample_audit["sample_in_rna"] = sample_audit["n_rna_cells"].gt(0)
    sample_audit["sample_in_raw"] = sample_audit["n_productive_contigs"].gt(0)
    sample_audit["sample_mapping_status"] = np.select(
        [
            sample_audit["n_raw_files"].gt(1),
            sample_audit["sample_in_rna"] & sample_audit["sample_in_raw"],
            sample_audit["sample_in_raw"] & ~sample_audit["sample_in_rna"],
        ],
        ["quarantined_multiple_raw_libraries", "eligible", "raw_sample_not_in_rna"],
        default="rna_sample_without_raw_vdj",
    )

    merged = rna.reset_index(drop=True).merge(wide, on=["sample_id", "barcode_core"], how="left", validate="many_to_one")
    merged.index = obs.index
    merged.insert(0, "source_gse_id", gse_id)
    merged = merged.rename(columns={"sample_id": "join_sample_id", "donor_id": "join_donor_id"})
    for chain in CHAINS:
        for field in CHAIN_FIELDS:
            column = f"{chain}_{field}"
            if column not in merged.columns:
                merged[column] = "" if field in STRING_FIELDS else False if field in {"umi_available", "read_available"} else pd.NA
            if field in STRING_FIELDS:
                merged[column] = merged[column].map(clean_text)
            elif field in {"umis", "reads", "n_productive_contigs"}:
                merged[column] = pd.to_numeric(merged[column], errors="coerce").astype("Int64")
            elif field in {"umi_available", "read_available"}:
                merged[column] = merged[column].astype("boolean").fillna(False).astype(bool)

    in_raw_sample = merged["join_sample_id"].isin(set(raw["sample_id"]))
    raw_sample_ambiguous = merged["join_sample_id"].isin(ambiguous_samples)
    missing_key = merged["join_sample_id"].eq("") | merged["barcode_core"].eq("")
    ambiguous_rna = merged["rna_key_multiplicity"].gt(1)
    invalid_row = missing_key | ambiguous_rna | raw_sample_ambiguous
    for chain in CHAINS:
        for field in CHAIN_FIELDS:
            column = f"{chain}_{field}"
            if field in STRING_FIELDS:
                merged.loc[invalid_row, column] = ""
            elif field in {"umi_available", "read_available"}:
                merged.loc[invalid_row, column] = False
            else:
                merged.loc[invalid_row, column] = pd.NA

    has_chain = {chain: merged[f"{chain}_cdr3"].ne("") for chain in CHAINS}
    merged["has_any_ab_tcr_rebuilt"] = has_chain["TRA"] | has_chain["TRB"]
    merged["has_any_gd_tcr_rebuilt"] = has_chain["TRG"] | has_chain["TRD"]
    merged["has_TRA_TRB_paired_rebuilt"] = has_chain["TRA"] & has_chain["TRB"]
    merged["has_TRG_TRD_paired_rebuilt"] = has_chain["TRG"] & has_chain["TRD"]
    merged["TCRseq_rebuilt"] = np.where(
        merged["has_any_ab_tcr_rebuilt"] | merged["has_any_gd_tcr_rebuilt"], "yes", "no"
    )

    matched = merged["TCRseq_rebuilt"].eq("yes")
    merged["join_status"] = np.select(
        [missing_key, ambiguous_rna, raw_sample_ambiguous, ~in_raw_sample, matched],
        [
            "missing_rna_join_key",
            "ambiguous_rna_join_key",
            "raw_sample_quarantined",
            "no_raw_vdj_for_sample",
            "matched_productive_tcr",
        ],
        default="raw_sample_no_productive_tcr",
    )
    merged["join_reason"] = np.select(
        [missing_key, ambiguous_rna, raw_sample_ambiguous],
        [
            "sample_id or barcode_core is unavailable",
            "sample_id + barcode_core maps to multiple RNA rows",
            "multiple raw VDJ libraries collapse to one RNA sample key",
        ],
        default="",
    )
    merged["tcr_assignment_eligible"] = ~(missing_key | ambiguous_rna | raw_sample_ambiguous)
    # A passing source replaces all stale rows. Ambiguous rows are deliberately
    # written as blank and remain ineligible as receptor truth.
    merged["replacement_eligible"] = True

    old = old_tcr_flags(obs).reset_index(drop=True)
    old.index = merged.index
    merged = pd.concat([merged, old], axis=1)
    stats = {
        "n_ambiguous_raw_samples": len(ambiguous_samples),
        "n_ambiguous_rna_rows": int(ambiguous_rna.sum()),
        "n_missing_rna_keys": int(missing_key.sum()),
        "n_eligible_raw_cell_keys": int(len(wide)),
        "n_matched_raw_cell_keys": int(matched.sum()),
        "source_mapping_quarantined": 0,
    }
    return merged, sample_audit, stats


def receptor_reuse_metrics(gse_id: str, sidecar: pd.DataFrame) -> dict[str, object]:
    paired = sidecar.loc[
        sidecar["has_TRA_TRB_paired_rebuilt"],
        ["join_sample_id", "join_donor_id", "TRA_cdr3", "TRB_cdr3"],
    ].copy()
    if paired.empty:
        return {
            "source_gse_id": gse_id,
            "donor_key_status": UNRESOLVED_DONOR_SOURCES.get(gse_id, "resolved_or_sample_derived"),
            "n_paired_ab": 0,
            "n_unique_pairs": 0,
            "n_pairs_across_samples": 0,
            "n_cells_in_pairs_across_samples": 0,
            "max_samples_per_pair": 0,
            "n_pairs_across_donors": 0,
            "n_cells_in_pairs_across_donors": 0,
            "max_donors_per_pair": 0,
        }
    paired["pair"] = paired["TRA_cdr3"] + "|" + paired["TRB_cdr3"]
    per_pair = paired.groupby("pair").agg(
        n_cells=("pair", "size"),
        n_samples=("join_sample_id", "nunique"),
        n_donors=("join_donor_id", "nunique"),
    )
    cross_sample = per_pair["n_samples"].gt(1)
    cross_donor = per_pair["n_donors"].gt(1)
    result = {
        "source_gse_id": gse_id,
        "donor_key_status": UNRESOLVED_DONOR_SOURCES.get(gse_id, "resolved_or_sample_derived"),
        "n_paired_ab": len(paired),
        "n_unique_pairs": len(per_pair),
        "n_pairs_across_samples": int(cross_sample.sum()),
        "n_cells_in_pairs_across_samples": int(per_pair.loc[cross_sample, "n_cells"].sum()),
        "max_samples_per_pair": int(per_pair["n_samples"].max()),
        "n_pairs_across_donors": int(cross_donor.sum()),
        "n_cells_in_pairs_across_donors": int(per_pair.loc[cross_donor, "n_cells"].sum()),
        "max_donors_per_pair": int(per_pair["n_donors"].max()),
    }
    if gse_id in UNRESOLVED_DONOR_SOURCES:
        for column in (
            "n_pairs_across_donors",
            "n_cells_in_pairs_across_donors",
            "max_donors_per_pair",
        ):
            result[column] = np.nan
    return result


def author_lineage_conflicts(obs: pd.DataFrame, sidecar: pd.DataFrame) -> dict[str, object]:
    """Summarize rebuilt alpha/beta calls in author-labeled NK or non-T cells."""
    annotation_columns = [
        column
        for column in (
            "cell_type",
            "celltype",
            "cluster_name",
            "lineage",
            "major_cell_type",
            "original_cell_annotation",
            "predicted.celltype.l2",
        )
        if column in obs.columns
    ]
    if not annotation_columns:
        return {
            "author_annotation_columns": "",
            "n_author_nk_cells": 0,
            "n_rebuilt_ab_in_author_nk": 0,
            "n_author_non_t_cells": 0,
            "n_rebuilt_ab_in_author_non_t": 0,
        }
    labels = obs[annotation_columns].astype(str).agg(" | ".join, axis=1).str.lower()
    nk = labels.str.contains(r"(?:^|[^a-z])nk(?:[^a-z]|$)|natural killer", regex=True, na=False)
    non_t = labels.str.contains(
        r"(?:^|[^a-z])(?:b cell|b-cell|plasma|myeloid|monocyte|macrophage|dendritic|epithelial|fibroblast|endothelial|erythroid)(?:[^a-z]|$)",
        regex=True,
        na=False,
    )
    rebuilt_ab = sidecar["has_any_ab_tcr_rebuilt"].to_numpy()
    return {
        "author_annotation_columns": ";".join(annotation_columns),
        "n_author_nk_cells": int(nk.sum()),
        "n_rebuilt_ab_in_author_nk": int((nk.to_numpy() & rebuilt_ab).sum()),
        "n_author_non_t_cells": int(non_t.sum()),
        "n_rebuilt_ab_in_author_non_t": int((non_t.to_numpy() & rebuilt_ab).sum()),
    }


def summarize_source(
    gse_id: str,
    obs: pd.DataFrame,
    sidecar: pd.DataFrame,
    bundle: RawBundle,
    stats: dict[str, int],
    source_reason: str,
    source_stat_before: tuple[int, int],
    source_stat_after: tuple[int, int],
) -> dict[str, object]:
    new_any = sidecar["has_any_ab_tcr_rebuilt"] | sidecar["has_any_gd_tcr_rebuilt"]
    old_any = sidecar["old_has_any_ab"] | sidecar["old_has_any_gd"]
    umi_calls = pd.concat(
        [
            pd.to_numeric(sidecar.loc[sidecar[f"{chain}_cdr3"].ne(""), f"{chain}_umis"], errors="coerce")
            for chain in CHAINS
        ],
        ignore_index=True,
    )
    if SOURCE_SPECS[gse_id].loader == "quarantine":
        status = "QUARANTINED_NO_SAFE_RAW_JOIN"
    elif stats["source_mapping_quarantined"]:
        status = "QUARANTINED_INSUFFICIENT_RAW_TO_RNA_MATCH"
    elif new_any.sum() == 0:
        status = "FAIL_NO_PRODUCTIVE_TCR_MATCHES"
    elif stats["n_ambiguous_raw_samples"] or stats["n_ambiguous_rna_rows"]:
        status = "PASS_PARTIAL_WITH_QUARANTINE"
    else:
        status = "PASS_REBUILT"
    result = {
        "source_gse_id": gse_id,
        "rebuild_status": status,
        "source_reason": source_reason,
        "source_h5ad": relative_label(source_h5ad_path(gse_id)),
        "n_source_cells": len(obs),
        "n_raw_source_files": len(bundle.source_files),
        "n_raw_input_rows": bundle.n_input_rows,
        "n_raw_productive_rows_before_chain_filter": bundle.n_productive_rows,
        "n_standardized_productive_chain_contigs": int(bundle.contigs.shape[0]),
        "n_old_any_tcr": int(old_any.sum()),
        "n_rebuilt_any_tcr": int(new_any.sum()),
        "n_rebuilt_any_ab": int(sidecar["has_any_ab_tcr_rebuilt"].sum()),
        "n_rebuilt_paired_ab": int(sidecar["has_TRA_TRB_paired_rebuilt"].sum()),
        "n_rebuilt_any_gd": int(sidecar["has_any_gd_tcr_rebuilt"].sum()),
        "n_rebuilt_paired_gd": int(sidecar["has_TRG_TRD_paired_rebuilt"].sum()),
        "n_old_calls_removed": int((old_any & ~new_any).sum()),
        "n_new_calls_added": int((~old_any & new_any).sum()),
        "n_calls_agree_any_tcr": int((old_any == new_any).sum()),
        "n_rebuilt_chain_calls_with_umi": int(umi_calls.notna().sum()),
        "n_rebuilt_chain_calls_umi_eq_1": int(umi_calls.eq(1).sum()),
        "n_ambiguous_raw_samples": stats["n_ambiguous_raw_samples"],
        "n_ambiguous_rna_rows": stats["n_ambiguous_rna_rows"],
        "n_missing_rna_keys": stats["n_missing_rna_keys"],
        "n_eligible_raw_cell_keys": stats["n_eligible_raw_cell_keys"],
        "n_matched_raw_cell_keys": stats["n_matched_raw_cell_keys"],
        "raw_to_rna_match_fraction": (
            stats["n_matched_raw_cell_keys"] / stats["n_eligible_raw_cell_keys"]
            if stats["n_eligible_raw_cell_keys"]
            else np.nan
        ),
        "n_replacement_eligible_rows": int(sidecar["replacement_eligible"].sum()),
        "n_tcr_assignment_eligible_rows": int(sidecar["tcr_assignment_eligible"].sum()),
        "source_h5ad_unchanged": source_stat_before == source_stat_after,
        "source_h5ad_size": source_stat_after[0],
        "source_h5ad_mtime_ns": source_stat_after[1],
    }
    result.update(author_lineage_conflicts(obs, sidecar))
    return result


def write_markdown(summary: pd.DataFrame, sample_audit: pd.DataFrame) -> None:
    pass_count = int(summary["rebuild_status"].str.startswith("PASS").sum())
    quarantined = int(summary["rebuild_status"].str.startswith("QUARANTINED").sum())
    lines = [
        "# Sample-aware TCR join rebuild",
        "",
        "## Result",
        "",
        f"- Sources passing a full or fail-closed partial rebuild: **{pass_count}/{len(summary)}**",
        f"- Sources quarantined because no safe raw join can be reconstructed: **{quarantined}**",
        f"- Rebuilt cells with any productive TCR: **{int(summary['n_rebuilt_any_tcr'].sum()):,}**",
        f"- Rebuilt productive chain calls with observed UMI support: **{int(summary['n_rebuilt_chain_calls_with_umi'].sum()):,}**",
        "- No H5AD was modified. Outputs are metadata sidecars only.",
        "",
        "## Method",
        "",
        "Productive, cell-associated, high-confidence raw contigs were reduced to one contig per chain by highest UMI, then reads, full-length status, and contig id. TRA, TRB, TRG, and TRD were joined only by canonical sample id plus barcode core. Missing UMI/read fields remain null. Multiple raw libraries collapsing to one sample key and duplicated RNA join keys fail closed.",
        "",
        "## Source summary",
        "",
        summary[
            [
                "source_gse_id",
                "rebuild_status",
                "n_source_cells",
                "n_rebuilt_any_tcr",
                "n_rebuilt_paired_ab",
                "n_rebuilt_any_gd",
                "n_rebuilt_chain_calls_with_umi",
                "n_ambiguous_raw_samples",
            ]
        ].to_markdown(index=False),
        "",
        "## Quarantined sample mappings",
        "",
    ]
    quarantined_samples = sample_audit.loc[
        ~sample_audit["sample_mapping_status"].isin({"eligible", "rna_sample_without_raw_vdj"})
    ]
    lines.append(
        quarantined_samples.to_markdown(index=False)
        if len(quarantined_samples)
        else "No sample-level mapping was quarantined among rebuildable sources."
    )
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def write_staged_replacement(summary: pd.DataFrame) -> dict[str, object]:
    """Combine passing source sidecars without modifying any H5AD."""
    passing = summary.loc[summary["rebuild_status"].str.startswith("PASS")].copy()
    columns = [
        "source_gse_id",
        "source_obs_name",
        "join_sample_id",
        "join_donor_id",
        "barcode_core",
        "join_status",
        "join_reason",
        "tcr_assignment_eligible",
        "replacement_eligible",
        *[f"{chain}_{field}" for chain in CHAINS for field in CHAIN_FIELDS],
        "has_any_ab_tcr_rebuilt",
        "has_any_gd_tcr_rebuilt",
        "has_TRA_TRB_paired_rebuilt",
        "has_TRG_TRD_paired_rebuilt",
        "TCRseq_rebuilt",
    ]
    tmp_path = STAGED_REPLACEMENT.with_suffix(STAGED_REPLACEMENT.suffix + ".tmp")
    if tmp_path.exists():
        tmp_path.unlink()
    writer: pq.ParquetWriter | None = None
    manifest_rows: list[dict[str, object]] = []
    try:
        for row in passing.itertuples(index=False):
            path = PER_SOURCE_DIR / f"{row.source_gse_id}.parquet"
            frame = pd.read_parquet(path, columns=columns)
            frame = frame.loc[frame["replacement_eligible"]].copy()
            table = pa.Table.from_pandas(frame, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(tmp_path, table.schema, compression="zstd")
            elif table.schema != writer.schema:
                table = table.cast(writer.schema)
            writer.write_table(table)
            manifest_rows.append(
                {
                    "source_gse_id": row.source_gse_id,
                    "rebuild_status": row.rebuild_status,
                    "n_replacement_rows": len(frame),
                    "n_tcr_assignment_eligible_rows": int(frame["tcr_assignment_eligible"].sum()),
                    "n_rebuilt_any_tcr": int(frame["TCRseq_rebuilt"].eq("yes").sum()),
                }
            )
    finally:
        if writer is not None:
            writer.close()
    if writer is None:
        raise RuntimeError("No passing sources were available for the staged replacement.")
    tmp_path.replace(STAGED_REPLACEMENT)
    manifest_df = pd.DataFrame(manifest_rows)
    manifest_df.to_csv(STAGED_REPLACEMENT_MANIFEST, index=False)
    return {
        "path": relative_label(STAGED_REPLACEMENT),
        "sha256": sha256_file(STAGED_REPLACEMENT),
        "n_rows": int(manifest_df["n_replacement_rows"].sum()),
        "n_sources": len(manifest_df),
        "manifest": relative_label(STAGED_REPLACEMENT_MANIFEST),
    }


def process_source(gse_id: str, overwrite: bool) -> tuple[dict[str, object], pd.DataFrame, dict[str, object]]:
    if gse_id not in SOURCE_SPECS:
        raise KeyError(f"No source specification for {gse_id}")
    spec = SOURCE_SPECS[gse_id]
    h5ad_path = source_h5ad_path(gse_id)
    resolved = h5ad_path.resolve()
    before = (resolved.stat().st_size, resolved.stat().st_mtime_ns)
    logging.info("Reading %s from %s", gse_id, h5ad_path)
    adata = ad.read_h5ad(h5ad_path, backed="r")
    obs = adata.obs.copy()
    adata.file.close()

    bundle = load_raw_bundle(gse_id, spec)
    if spec.loader == "quarantine":
        sidecar = build_quarantine_sidecar(gse_id, obs, spec.source_reason)
        old = old_tcr_flags(obs)
        sidecar = pd.concat([sidecar, old], axis=1)
        sidecar["has_any_ab_tcr_rebuilt"] = False
        sidecar["has_any_gd_tcr_rebuilt"] = False
        sidecar["has_TRA_TRB_paired_rebuilt"] = False
        sidecar["has_TRG_TRD_paired_rebuilt"] = False
        sidecar["TCRseq_rebuilt"] = "no"
        sample_audit = pd.DataFrame(
            [{
                "source_gse_id": gse_id,
                "sample_id": "",
                "n_rna_cells": len(obs),
                "n_productive_contigs": 0,
                "n_raw_files": 0,
                "sample_in_rna": True,
                "sample_in_raw": False,
                "sample_mapping_status": "source_quarantined",
            }]
        )
        stats = {
            "n_ambiguous_raw_samples": 0,
            "n_ambiguous_rna_rows": 0,
            "n_missing_rna_keys": len(obs),
            "n_eligible_raw_cell_keys": 0,
            "n_matched_raw_cell_keys": 0,
            "source_mapping_quarantined": 1,
        }
    else:
        sidecar, sample_audit, stats = build_sidecar(gse_id, obs, bundle)
        match_fraction = (
            stats["n_matched_raw_cell_keys"] / stats["n_eligible_raw_cell_keys"]
            if stats["n_eligible_raw_cell_keys"]
            else 0.0
        )
        if stats["n_eligible_raw_cell_keys"] and match_fraction < 0.5:
            stats["source_mapping_quarantined"] = 1
            sidecar["replacement_eligible"] = False
            sidecar["tcr_assignment_eligible"] = False
            sidecar["join_reason"] = (
                "source quarantined because fewer than 50% of eligible raw VDJ cell keys map to RNA"
            )
            for chain in CHAINS:
                for field in CHAIN_FIELDS:
                    column = f"{chain}_{field}"
                    if field in STRING_FIELDS:
                        sidecar[column] = ""
                    elif field in {"umi_available", "read_available"}:
                        sidecar[column] = False
                    else:
                        sidecar[column] = pd.Series(pd.NA, index=sidecar.index, dtype="Int64")
            sidecar["has_any_ab_tcr_rebuilt"] = False
            sidecar["has_any_gd_tcr_rebuilt"] = False
            sidecar["has_TRA_TRB_paired_rebuilt"] = False
            sidecar["has_TRG_TRD_paired_rebuilt"] = False
            sidecar["TCRseq_rebuilt"] = "no"

    sidecar_path = PER_SOURCE_DIR / f"{gse_id}.parquet"
    if sidecar_path.exists() and not overwrite:
        raise FileExistsError(f"Sidecar exists; rerun with --overwrite: {sidecar_path}")
    sidecar.reset_index(drop=True).to_parquet(sidecar_path, index=False, compression="zstd")
    after = (resolved.stat().st_size, resolved.stat().st_mtime_ns)
    summary = summarize_source(gse_id, obs, sidecar, bundle, stats, spec.source_reason, before, after)
    reuse = receptor_reuse_metrics(gse_id, sidecar)
    logging.info(
        "%s: %s, rebuilt any TCR %s/%s",
        gse_id,
        summary["rebuild_status"],
        f"{summary['n_rebuilt_any_tcr']:,}",
        f"{len(obs):,}",
    )
    return summary, sample_audit, reuse


def main() -> None:
    args = parse_args()
    ensure_dirs()
    configure_logging()
    source_rows: list[dict[str, object]] = []
    sample_frames: list[pd.DataFrame] = []
    reuse_rows: list[dict[str, object]] = []
    for gse_id in args.gse:
        summary, sample_audit, reuse = process_source(gse_id, args.overwrite)
        source_rows.append(summary)
        sample_frames.append(sample_audit)
        reuse_rows.append(reuse)

    summary_df = pd.DataFrame(source_rows)
    sample_df = pd.concat(sample_frames, ignore_index=True)
    reuse_df = pd.DataFrame(reuse_rows)
    summary_df.to_csv(SUMMARY_CSV, index=False)
    sample_df.to_csv(SAMPLE_AUDIT_CSV, index=False)
    reuse_df.to_csv(REUSE_AUDIT_CSV, index=False)
    write_markdown(summary_df, sample_df)
    staged_replacement = write_staged_replacement(summary_df)
    manifest = {
        "workflow": relative_label(Path(__file__)),
        "source_audit": relative_label(SOURCE_AUDIT),
        "sources": list(args.gse),
        "h5ad_write_performed": False,
        "sidecar_format": "parquet-zstd",
        "join_key": ["sample_id", "barcode_core"],
        "chains": list(CHAINS),
        "outputs": {
            "source_summary": relative_label(SUMMARY_CSV),
            "sample_audit": relative_label(SAMPLE_AUDIT_CSV),
            "receptor_reuse_audit": relative_label(REUSE_AUDIT_CSV),
            "summary_markdown": relative_label(SUMMARY_MD),
            "per_source_sidecars": relative_label(PER_SOURCE_DIR),
            "staged_replacement": relative_label(STAGED_REPLACEMENT),
            "staged_replacement_manifest": relative_label(STAGED_REPLACEMENT_MANIFEST),
        },
        "staged_replacement": staged_replacement,
    }
    RUN_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    logging.info("Wrote %s", SUMMARY_CSV)


if __name__ == "__main__":
    main()
