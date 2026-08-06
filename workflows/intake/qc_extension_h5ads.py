#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Run the standalone Phase 0 review gate for built extension H5ADs.

The workflow is deliberately read-only. It validates each built cohort without
merging objects, streams sparse count data without materializing a dense
expression matrix, and writes review artifacts only after all inputs have been
audited.
"""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any

import anndata as ad
import h5py
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = (
    PROJECT_ROOT / "data/interim/extension_intake/built_h5ads_manifest.csv"
)
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
REVIEW_SUBDIR = "extension_intake"

CHAINS = ("TRA", "TRB", "TRG", "TRD")
STRING_SUFFIXES = ("cdr3", "v", "d", "j", "cdr3_nt", "clone_id")
NUMERIC_SUFFIXES = ("umis", "reads")
TCR_COLUMNS = tuple(
    f"{chain}_{suffix}"
    for chain in CHAINS
    for suffix in (*STRING_SUFFIXES, *NUMERIC_SUFFIXES)
)
TCR_FLAGS = (
    "has_TRA",
    "has_TRB",
    "has_TRG",
    "has_TRD",
    "has_TRA_TRB_paired",
    "has_TRG_TRD_paired",
    "has_any_ab_tcr",
    "has_any_gd_tcr",
)
REQUIRED_HARMONIZED_METADATA = (
    "sample_id",
    "library_id",
    "donor_id",
    "tissue",
    "diagnosis",
    "technology_simple",
    "source_accession",
    "barcode",
    "barcode_core",
    "tcr_schema_provenance",
)

COHORT_ALIASES = ("cohort_id", "dataset_id", "gse_id")
PATH_ALIASES = (
    "h5ad_path",
    "output_h5ad",
    "output_path",
    "built_h5ad",
    "path",
)
STATUS_ALIASES = ("build_status", "status")
CELL_COUNT_ALIASES = ("n_cells", "cells", "n_obs")
SUCCESS_STATUSES = {
    "built",
    "complete",
    "completed",
    "ok",
    "pass",
    "passed",
    "ready",
    "success",
    "true",
    "1",
}
TRUE_VALUES = {"true", "t", "yes", "y", "1"}
FALSE_VALUES = {"false", "f", "no", "n", "0"}
MISSING_TEXT = {"", "na", "nan", "none", "null", "<na>", "n/a"}
MAX_SPARSE_VALUES_PER_BLOCK = 2_000_000
MAX_MAJOR_AXIS_PER_BLOCK = 50_000


class Phase0QCError(RuntimeError):
    """Raised when the Phase 0 gate cannot safely audit an input."""


@dataclass(frozen=True)
class ManifestEntry:
    cohort_id: str
    h5ad_path: Path
    expected_cells: int | None
    manifest_row: int


@dataclass
class SparseMatrixReview:
    encoding: str
    dtype: str
    nnz: int
    explicit_zero_entries: int
    count_min: float
    count_max: float
    count_sum: float
    total_counts: np.ndarray
    detected_genes: np.ndarray
    mitochondrial_counts: np.ndarray


def clean_scalar(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    return "" if text.casefold() in MISSING_TEXT else text


def clean_text_series(series: pd.Series) -> pd.Series:
    text = series.astype("string").fillna("").str.strip()
    return text.mask(text.str.casefold().isin(MISSING_TEXT), "")


def nonempty_text(series: pd.Series) -> np.ndarray:
    return clean_text_series(series).ne("").to_numpy(dtype=bool)


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(block_size):
            digest.update(chunk)
    return digest.hexdigest()


def first_present(columns: Sequence[str], aliases: Sequence[str]) -> str | None:
    normalized = {str(column).strip().casefold(): column for column in columns}
    for alias in aliases:
        if alias.casefold() in normalized:
            return normalized[alias.casefold()]
    return None


def resolve_input_path(raw_path: str, manifest_path: Path) -> Path:
    candidate = Path(raw_path).expanduser()
    if candidate.is_absolute():
        return candidate.resolve()
    choices = (manifest_path.parent / candidate, PROJECT_ROOT / candidate)
    for choice in choices:
        if choice.exists():
            return choice.resolve()
    return choices[0].resolve()


def load_manifest(manifest_path: Path) -> tuple[list[ManifestEntry], pd.DataFrame]:
    if not manifest_path.is_file():
        raise Phase0QCError(f"Built-H5AD manifest does not exist: {manifest_path}")
    frame = pd.read_csv(
        manifest_path,
        dtype=str,
        keep_default_na=False,
        na_filter=False,
    )
    if frame.empty:
        raise Phase0QCError(f"Built-H5AD manifest is empty: {manifest_path}")

    cohort_column = first_present(frame.columns, COHORT_ALIASES)
    path_column = first_present(frame.columns, PATH_ALIASES)
    if cohort_column is None or path_column is None:
        raise Phase0QCError(
            "Built-H5AD manifest must contain a cohort identifier "
            f"({COHORT_ALIASES}) and H5AD path ({PATH_ALIASES})"
        )
    status_column = first_present(frame.columns, STATUS_ALIASES)
    count_column = first_present(frame.columns, CELL_COUNT_ALIASES)

    entries: list[ManifestEntry] = []
    audit_rows: list[dict[str, Any]] = []
    for offset, row in frame.iterrows():
        manifest_row = int(offset) + 2
        cohort_id = clean_scalar(row[cohort_column])
        raw_path = clean_scalar(row[path_column])
        status = clean_scalar(row[status_column]).casefold() if status_column else "built"
        if not cohort_id:
            raise Phase0QCError(f"Manifest row {manifest_row} has a blank cohort ID")
        if not raw_path:
            raise Phase0QCError(f"Manifest row {manifest_row} has a blank H5AD path")
        if status not in SUCCESS_STATUSES:
            raise Phase0QCError(
                f"Manifest row {manifest_row} for {cohort_id} is not built: {status!r}"
            )
        path = resolve_input_path(raw_path, manifest_path)
        if path.suffix.casefold() != ".h5ad":
            raise Phase0QCError(
                f"Manifest row {manifest_row} does not reference an H5AD: {path}"
            )
        if not path.is_file():
            raise Phase0QCError(
                f"Manifest row {manifest_row} references a missing H5AD: {path}"
            )

        expected_cells: int | None = None
        if count_column and clean_scalar(row[count_column]):
            try:
                expected_cells = int(clean_scalar(row[count_column]))
            except ValueError as exc:
                raise Phase0QCError(
                    f"Manifest row {manifest_row} has invalid expected cell count "
                    f"{row[count_column]!r}"
                ) from exc
            if expected_cells < 0:
                raise Phase0QCError(
                    f"Manifest row {manifest_row} has negative expected cell count"
                )

        entries.append(
            ManifestEntry(
                cohort_id=cohort_id,
                h5ad_path=path,
                expected_cells=expected_cells,
                manifest_row=manifest_row,
            )
        )
        audit_rows.append(
            {
                "manifest_row": manifest_row,
                "cohort_id": cohort_id,
                "h5ad_path": str(path),
                "status": status,
                "expected_cells": expected_cells,
                "path_exists": True,
            }
        )

    cohort_ids = pd.Series([entry.cohort_id for entry in entries])
    duplicate_cohorts = sorted(cohort_ids[cohort_ids.duplicated(keep=False)].unique())
    if duplicate_cohorts:
        raise Phase0QCError(
            "Manifest contains more than one built H5AD for cohort(s): "
            + ", ".join(duplicate_cohorts)
        )
    paths = pd.Series([str(entry.h5ad_path) for entry in entries])
    duplicate_paths = sorted(paths[paths.duplicated(keep=False)].unique())
    if duplicate_paths:
        raise Phase0QCError(
            "Manifest reuses built H5AD path(s): " + ", ".join(duplicate_paths)
        )
    return entries, pd.DataFrame(audit_rows)


def decode_h5_attr(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.ndarray) and value.size == 1:
        return decode_h5_attr(value.item())
    return str(value)


def major_axis_blocks(indptr: np.ndarray) -> list[tuple[int, int]]:
    """Bound both major-axis size and sparse-value count for each read."""
    blocks: list[tuple[int, int]] = []
    n_major = len(indptr) - 1
    start = 0
    while start < n_major:
        max_end = min(n_major, start + MAX_MAJOR_AXIS_PER_BLOCK)
        target = int(indptr[start]) + MAX_SPARSE_VALUES_PER_BLOCK
        by_values = int(np.searchsorted(indptr, target, side="right") - 1)
        end = min(max_end, max(start + 1, by_values))
        blocks.append((start, end))
        start = end
    return blocks


def validate_sparse_values(values: np.ndarray, label: str) -> None:
    if not np.issubdtype(values.dtype, np.number) or np.issubdtype(
        values.dtype, np.bool_
    ):
        raise Phase0QCError(f"{label} has a non-numeric count dtype: {values.dtype}")
    if values.size == 0:
        return
    if not np.isfinite(values).all():
        raise Phase0QCError(f"{label} contains non-finite count values")
    if np.any(values < 0):
        raise Phase0QCError(f"{label} contains negative count values")
    if not np.issubdtype(values.dtype, np.integer) and not np.equal(
        values, np.floor(values)
    ).all():
        raise Phase0QCError(f"{label} contains non-integer count values")


def validate_indptr(indptr: np.ndarray, nnz: int, expected_length: int, label: str) -> None:
    if len(indptr) != expected_length:
        raise Phase0QCError(
            f"{label} indptr length is {len(indptr)}, expected {expected_length}"
        )
    if not np.issubdtype(indptr.dtype, np.integer):
        raise Phase0QCError(f"{label} indptr is not integer typed")
    if int(indptr[0]) != 0 or int(indptr[-1]) != nnz:
        raise Phase0QCError(
            f"{label} indptr endpoints are inconsistent with nnz={nnz}"
        )
    if np.any(np.diff(indptr) < 0):
        raise Phase0QCError(f"{label} indptr is not monotonic")


def inspect_sparse_x(
    path: Path,
    expected_shape: tuple[int, int],
    mitochondrial_gene_mask: np.ndarray,
) -> SparseMatrixReview:
    """Validate and summarize H5AD X using sparse on-disk arrays only."""
    with h5py.File(path, "r") as handle:
        if "X" not in handle:
            raise Phase0QCError(f"{path} has no X matrix")
        group = handle["X"]
        if not isinstance(group, h5py.Group):
            raise Phase0QCError(f"{path} stores dense X; sparse raw counts are required")
        encoding = decode_h5_attr(group.attrs.get("encoding-type", ""))
        if encoding not in {"csr_matrix", "csc_matrix"}:
            raise Phase0QCError(
                f"{path} X encoding is {encoding!r}; expected csr_matrix or csc_matrix"
            )
        missing = [name for name in ("data", "indices", "indptr") if name not in group]
        if missing:
            raise Phase0QCError(f"{path} sparse X lacks arrays: {missing}")
        raw_shape = group.attrs.get("shape")
        if raw_shape is None:
            raise Phase0QCError(f"{path} sparse X lacks its shape attribute")
        shape = tuple(int(value) for value in np.asarray(raw_shape).tolist())
        if shape != expected_shape:
            raise Phase0QCError(
                f"{path} sparse X shape {shape} does not match AnnData {expected_shape}"
            )

        data_ds = group["data"]
        indices_ds = group["indices"]
        indptr_ds = group["indptr"]
        count_dtype = str(data_ds.dtype)
        if len(data_ds) != len(indices_ds):
            raise Phase0QCError(f"{path} sparse X data/indices lengths differ")
        if not np.issubdtype(indices_ds.dtype, np.integer):
            raise Phase0QCError(f"{path} sparse X indices are not integer typed")

        nnz = int(len(data_ds))
        n_obs, n_vars = shape
        major_size = n_obs if encoding == "csr_matrix" else n_vars
        minor_size = n_vars if encoding == "csr_matrix" else n_obs
        indptr = np.asarray(indptr_ds[:])
        validate_indptr(indptr, nnz, major_size + 1, str(path))

        total_counts = np.zeros(n_obs, dtype=np.float64)
        detected_genes = np.zeros(n_obs, dtype=np.int64)
        mitochondrial_counts = np.zeros(n_obs, dtype=np.float64)
        explicit_zeros = 0
        count_min = math.inf
        count_max = -math.inf
        count_sum = 0.0

        for major_start, major_end in major_axis_blocks(indptr):
            value_start = int(indptr[major_start])
            value_end = int(indptr[major_end])
            values = np.asarray(data_ds[value_start:value_end])
            indices = np.asarray(indices_ds[value_start:value_end])
            validate_sparse_values(values, f"{path} X")
            if indices.size and (np.any(indices < 0) or np.any(indices >= minor_size)):
                raise Phase0QCError(f"{path} sparse X contains out-of-range indices")

            if values.size:
                explicit_zeros += int(np.count_nonzero(values == 0))
                count_min = min(count_min, float(values.min()))
                count_max = max(count_max, float(values.max()))
                count_sum += float(values.sum(dtype=np.float64))

            local_indptr = indptr[major_start : major_end + 1] - value_start
            cumulative = np.empty(values.size + 1, dtype=np.float64)
            cumulative[0] = 0.0
            np.cumsum(values, dtype=np.float64, out=cumulative[1:])
            positive = values > 0
            positive_cumulative = np.empty(values.size + 1, dtype=np.int64)
            positive_cumulative[0] = 0
            np.cumsum(positive, dtype=np.int64, out=positive_cumulative[1:])

            if encoding == "csr_matrix":
                total_counts[major_start:major_end] = (
                    cumulative[local_indptr[1:]] - cumulative[local_indptr[:-1]]
                )
                detected_genes[major_start:major_end] = (
                    positive_cumulative[local_indptr[1:]]
                    - positive_cumulative[local_indptr[:-1]]
                )
                mt_values = values * mitochondrial_gene_mask[indices]
                mt_cumulative = np.empty(values.size + 1, dtype=np.float64)
                mt_cumulative[0] = 0.0
                np.cumsum(mt_values, dtype=np.float64, out=mt_cumulative[1:])
                mitochondrial_counts[major_start:major_end] = (
                    mt_cumulative[local_indptr[1:]]
                    - mt_cumulative[local_indptr[:-1]]
                )
            else:
                if values.size:
                    np.add.at(total_counts, indices, values)
                    np.add.at(detected_genes, indices, positive.astype(np.int64))
                    column_ids = np.repeat(
                        np.arange(major_start, major_end, dtype=np.int64),
                        np.diff(local_indptr),
                    )
                    mt_values = values * mitochondrial_gene_mask[column_ids]
                    np.add.at(mitochondrial_counts, indices, mt_values)

    return SparseMatrixReview(
        encoding=encoding,
        dtype=count_dtype,
        nnz=nnz,
        explicit_zero_entries=explicit_zeros,
        count_min=0.0 if math.isinf(count_min) else count_min,
        count_max=0.0 if math.isinf(count_max) else count_max,
        count_sum=count_sum,
        total_counts=total_counts,
        detected_genes=detected_genes,
        mitochondrial_counts=mitochondrial_counts,
    )


def parse_boolean(series: pd.Series, field: str) -> np.ndarray:
    if pd.api.types.is_bool_dtype(series.dtype):
        if series.isna().any():
            raise Phase0QCError(f"{field} contains missing boolean values")
        return series.to_numpy(dtype=bool)
    values = clean_text_series(series).str.casefold()
    invalid = ~(values.isin(TRUE_VALUES) | values.isin(FALSE_VALUES))
    if invalid.any():
        examples = sorted(values.loc[invalid].unique().tolist())[:5]
        raise Phase0QCError(f"{field} has invalid boolean values: {examples}")
    return values.isin(TRUE_VALUES).to_numpy(dtype=bool)


def normalize_barcode_core(value: Any) -> str:
    text = clean_scalar(value).upper()
    embedded = re.search(r"([ACGTN]{8,}(?:-\d+)?)", text)
    full = embedded.group(1) if embedded else text
    sequence = re.match(r"([ACGTN]{8,})", full)
    if sequence:
        return sequence.group(1)
    return full.split("-")[0].split(".")[0].split(":")[-1]


def json_from_uns(value: Any, label: str) -> Any:
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    if isinstance(value, np.ndarray) and value.size == 1:
        value = value.item()
    if isinstance(value, str):
        try:
            return json.loads(value)
        except json.JSONDecodeError as exc:
            raise Phase0QCError(f"{label} is not valid JSON") from exc
    return value


def issue(
    cohort_id: str,
    severity: str,
    check: str,
    message: str,
    library_id: str = "",
) -> dict[str, str]:
    return {
        "cohort_id": cohort_id,
        "library_id": library_id,
        "severity": severity,
        "check": check,
        "message": message,
    }


def validate_required_obs(obs: pd.DataFrame, cohort_id: str) -> dict[str, pd.Series]:
    missing = [column for column in REQUIRED_HARMONIZED_METADATA if column not in obs]
    missing.extend(column for column in (*TCR_COLUMNS, *TCR_FLAGS) if column not in obs)
    if missing:
        raise Phase0QCError(f"Missing required obs columns: {sorted(set(missing))}")

    cleaned: dict[str, pd.Series] = {}
    for column in REQUIRED_HARMONIZED_METADATA:
        cleaned[column] = clean_text_series(obs[column])
        blanks = int(cleaned[column].eq("").sum())
        if blanks:
            raise Phase0QCError(f"{column} contains {blanks} blank value(s)")

    technologies = set(cleaned["technology_simple"].unique())
    if technologies != {"10x 5'"}:
        raise Phase0QCError(
            f"technology_simple must be exactly 10x 5'; observed {sorted(technologies)}"
        )
    if "cohort_id" in obs:
        observed = set(clean_text_series(obs["cohort_id"]).unique())
        if observed != {cohort_id}:
            raise Phase0QCError(
                f"obs cohort_id values {sorted(observed)} do not match manifest {cohort_id}"
            )
    return cleaned


def validate_join_keys(cleaned: Mapping[str, pd.Series]) -> None:
    key = pd.MultiIndex.from_arrays(
        [cleaned["sample_id"].to_numpy(), cleaned["barcode_core"].to_numpy()]
    )
    duplicated = key.duplicated(keep=False)
    if duplicated.any():
        raise Phase0QCError(
            "sample_id+barcode_core is not unique; "
            f"{int(duplicated.sum())} cell row(s) participate in duplicate keys"
        )
    expected_core = cleaned["barcode"].map(normalize_barcode_core)
    mismatch = expected_core.ne(cleaned["barcode_core"])
    if mismatch.any():
        raise Phase0QCError(
            f"barcode_core disagrees with normalized barcode in {int(mismatch.sum())} cell(s)"
        )


def validate_tcr_logic(obs: pd.DataFrame) -> dict[str, np.ndarray]:
    flags = {column: parse_boolean(obs[column], column) for column in TCR_FLAGS}
    chain_present: dict[str, np.ndarray] = {}
    for chain in CHAINS:
        cdr3_present = nonempty_text(obs[f"{chain}_cdr3"])
        chain_present[chain] = cdr3_present
        mismatch = flags[f"has_{chain}"] != cdr3_present
        if mismatch.any():
            raise Phase0QCError(
                f"has_{chain} disagrees with {chain}_cdr3 in {int(mismatch.sum())} cell(s)"
            )
        for suffix in STRING_SUFFIXES[1:]:
            orphan = nonempty_text(obs[f"{chain}_{suffix}"]) & ~cdr3_present
            if orphan.any():
                raise Phase0QCError(
                    f"{chain}_{suffix} is populated without productive {chain}_cdr3 "
                    f"in {int(orphan.sum())} cell(s)"
                )
        for suffix in NUMERIC_SUFFIXES:
            column = f"{chain}_{suffix}"
            raw = clean_text_series(obs[column])
            numeric = pd.to_numeric(raw.mask(raw.eq(""), "0"), errors="coerce")
            invalid = numeric.isna() | ~np.isfinite(numeric) | (numeric < 0)
            invalid |= numeric.ne(np.floor(numeric))
            if invalid.any():
                raise Phase0QCError(
                    f"{column} contains {int(invalid.sum())} invalid count value(s)"
                )
            orphan = numeric.to_numpy(dtype=np.float64) > 0
            orphan &= ~cdr3_present
            if orphan.any():
                raise Phase0QCError(
                    f"{column} is positive without productive {chain}_cdr3 in "
                    f"{int(orphan.sum())} cell(s)"
                )

    expected = {
        "has_TRA_TRB_paired": chain_present["TRA"] & chain_present["TRB"],
        "has_TRG_TRD_paired": chain_present["TRG"] & chain_present["TRD"],
        "has_any_ab_tcr": chain_present["TRA"] | chain_present["TRB"],
        "has_any_gd_tcr": chain_present["TRG"] | chain_present["TRD"],
    }
    for field, values in expected.items():
        mismatch = flags[field] != values
        if mismatch.any():
            raise Phase0QCError(
                f"{field} is logically inconsistent in {int(mismatch.sum())} cell(s)"
            )

    if "TCRseq" in obs:
        expected_tcrseq = expected["has_any_ab_tcr"] | expected["has_any_gd_tcr"]
        observed = clean_text_series(obs["TCRseq"]).str.casefold()
        invalid = ~observed.isin({"yes", "no"})
        mismatch = observed.eq("yes").to_numpy(dtype=bool) != expected_tcrseq
        if invalid.any() or mismatch.any():
            raise Phase0QCError(
                "TCRseq is inconsistent with canonical chain evidence in "
                f"{int((invalid.to_numpy() | mismatch).sum())} cell(s)"
            )
    return flags


def validate_library_metadata(cleaned: Mapping[str, pd.Series]) -> None:
    frame = pd.DataFrame(cleaned)
    for field in (
        "sample_id",
        "donor_id",
        "tissue",
        "diagnosis",
        "technology_simple",
        "source_accession",
        "tcr_schema_provenance",
    ):
        cardinality = frame.groupby("library_id", observed=True)[field].nunique()
        inconsistent = cardinality[cardinality != 1]
        if not inconsistent.empty:
            raise Phase0QCError(
                f"library_id does not map one-to-one to {field}: "
                + ", ".join(inconsistent.index.astype(str)[:5])
            )
    libraries_per_sample = frame.groupby("sample_id", observed=True)["library_id"].nunique()
    ambiguous_samples = libraries_per_sample[libraries_per_sample != 1]
    if not ambiguous_samples.empty:
        raise Phase0QCError(
            "sample_id does not map one-to-one to library_id: "
            + ", ".join(ambiguous_samples.index.astype(str)[:5])
        )


def parse_join_summaries(
    uns: Mapping[str, Any],
    cleaned: Mapping[str, pd.Series],
    flags: Mapping[str, np.ndarray],
) -> tuple[pd.DataFrame, str]:
    provenance = " ".join(
        sorted(cleaned["tcr_schema_provenance"].str.casefold().unique())
    )
    requires_external_join = "productive_contig" in provenance or "external_join" in provenance
    raw = uns.get("tcr_join_summaries")
    if raw is None:
        if requires_external_join:
            raise Phase0QCError(
                "Productive-contig provenance requires uns['tcr_join_summaries']"
            )
        return pd.DataFrame(), "embedded_or_not_applicable"

    payload = json_from_uns(raw, "uns['tcr_join_summaries']")
    if isinstance(payload, Mapping):
        payload = [payload]
    if not isinstance(payload, list) or not payload:
        raise Phase0QCError("uns['tcr_join_summaries'] must be a non-empty list")
    frame = pd.DataFrame(payload)
    required = {"sample_id", "expression_cells", "matched_cells", "unmatched_tcr_cells"}
    if not required.issubset(frame.columns):
        raise Phase0QCError(
            "TCR join summaries lack fields: " + ", ".join(sorted(required - set(frame.columns)))
        )
    frame["sample_id"] = clean_text_series(frame["sample_id"])
    if frame["sample_id"].eq("").any() or frame["sample_id"].duplicated().any():
        raise Phase0QCError("TCR join summaries have blank or duplicate sample_id values")
    for column in ("expression_cells", "matched_cells", "unmatched_tcr_cells"):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
        invalid = frame[column].isna() | (frame[column] < 0) | frame[column].ne(
            np.floor(frame[column])
        )
        if invalid.any():
            raise Phase0QCError(f"TCR join summary {column} contains invalid values")
        frame[column] = frame[column].astype(np.int64)
    if (frame["matched_cells"] > frame["expression_cells"]).any():
        raise Phase0QCError("TCR join matched_cells exceeds expression_cells")

    sample_series = cleaned["sample_id"]
    observed_cells = sample_series.value_counts(sort=False)
    tcr_positive = flags["has_any_ab_tcr"] | flags["has_any_gd_tcr"]
    observed_positive = pd.Series(tcr_positive).groupby(
        sample_series.reset_index(drop=True), observed=True
    ).sum()
    if set(frame["sample_id"]) != set(observed_cells.index):
        raise Phase0QCError(
            "TCR join summaries do not cover exactly the H5AD sample_id values"
        )
    by_sample = frame.set_index("sample_id")
    if not by_sample["expression_cells"].reindex(observed_cells.index).eq(
        observed_cells.astype(np.int64)
    ).all():
        raise Phase0QCError("TCR join expression_cells disagrees with H5AD sample counts")
    if not by_sample["matched_cells"].reindex(observed_positive.index).eq(
        observed_positive.astype(np.int64)
    ).all():
        raise Phase0QCError(
            "TCR join matched_cells disagrees with canonical TCR-positive cells"
        )
    frame["join_coverage"] = frame["matched_cells"] / frame["expression_cells"].replace(
        0, np.nan
    )
    return frame, "sample_id+barcode_core"


def quantile(values: pd.Series, q: float) -> float:
    return float(values.quantile(q)) if len(values) else math.nan


def make_library_rows(
    cohort_id: str,
    cleaned: Mapping[str, pd.Series],
    flags: Mapping[str, np.ndarray],
    matrix: SparseMatrixReview,
    join_frame: pd.DataFrame,
    join_mode: str,
) -> list[dict[str, Any]]:
    total_counts = matrix.total_counts
    pct_mt = np.divide(
        matrix.mitochondrial_counts * 100.0,
        total_counts,
        out=np.zeros_like(total_counts),
        where=total_counts > 0,
    )
    cell_frame = pd.DataFrame(
        {
            "library_id": cleaned["library_id"].to_numpy(),
            "sample_id": cleaned["sample_id"].to_numpy(),
            "donor_id": cleaned["donor_id"].to_numpy(),
            "tissue": cleaned["tissue"].to_numpy(),
            "diagnosis": cleaned["diagnosis"].to_numpy(),
            "source_accession": cleaned["source_accession"].to_numpy(),
            "tcr_schema_provenance": cleaned["tcr_schema_provenance"].to_numpy(),
            "total_counts": total_counts,
            "detected_genes": matrix.detected_genes,
            "pct_counts_mt": pct_mt,
        }
    )
    for field, values in flags.items():
        cell_frame[field] = values
    cell_frame["tcr_positive"] = flags["has_any_ab_tcr"] | flags["has_any_gd_tcr"]

    join_by_sample = join_frame.set_index("sample_id") if not join_frame.empty else None
    rows: list[dict[str, Any]] = []
    for library_id, group in cell_frame.groupby("library_id", sort=True, observed=True):
        sample_id = str(group["sample_id"].iloc[0])
        row: dict[str, Any] = {
            "cohort_id": cohort_id,
            "library_id": str(library_id),
            "sample_id": sample_id,
            "donor_id": str(group["donor_id"].iloc[0]),
            "tissue": str(group["tissue"].iloc[0]),
            "diagnosis": str(group["diagnosis"].iloc[0]),
            "source_accession": str(group["source_accession"].iloc[0]),
            "tcr_schema_provenance": str(group["tcr_schema_provenance"].iloc[0]),
            "n_cells": int(len(group)),
            "median_total_counts": float(group["total_counts"].median()),
            "q05_total_counts": quantile(group["total_counts"], 0.05),
            "q95_total_counts": quantile(group["total_counts"], 0.95),
            "median_detected_genes": float(group["detected_genes"].median()),
            "q05_detected_genes": quantile(group["detected_genes"], 0.05),
            "q95_detected_genes": quantile(group["detected_genes"], 0.95),
            "median_pct_counts_mt": float(group["pct_counts_mt"].median()),
            "q95_pct_counts_mt": quantile(group["pct_counts_mt"], 0.95),
            "zero_count_cells": int(group["total_counts"].eq(0).sum()),
            "cells_lt_200_genes": int(group["detected_genes"].lt(200).sum()),
            "cells_gt_20pct_mt": int(group["pct_counts_mt"].gt(20).sum()),
            "tcr_positive_cells": int(group["tcr_positive"].sum()),
            "has_TRA_cells": int(group["has_TRA"].sum()),
            "has_TRB_cells": int(group["has_TRB"].sum()),
            "has_TRG_cells": int(group["has_TRG"].sum()),
            "has_TRD_cells": int(group["has_TRD"].sum()),
            "paired_TRA_TRB_cells": int(group["has_TRA_TRB_paired"].sum()),
            "paired_TRG_TRD_cells": int(group["has_TRG_TRD_paired"].sum()),
            "tcr_positive_rate": float(group["tcr_positive"].mean()),
            "paired_TRA_TRB_rate": float(group["has_TRA_TRB_paired"].mean()),
            "paired_TRG_TRD_rate": float(group["has_TRG_TRD_paired"].mean()),
            "join_mode": join_mode,
            "join_expression_cells": math.nan,
            "join_matched_cells": math.nan,
            "unmatched_tcr_cells": math.nan,
            "join_coverage": math.nan,
        }
        if join_by_sample is not None:
            join = join_by_sample.loc[sample_id]
            row.update(
                {
                    "join_expression_cells": int(join["expression_cells"]),
                    "join_matched_cells": int(join["matched_cells"]),
                    "unmatched_tcr_cells": int(join["unmatched_tcr_cells"]),
                    "join_coverage": float(join["join_coverage"]),
                }
            )
        rows.append(row)
    return rows


def cohort_from_libraries(
    entry: ManifestEntry,
    adata: ad.AnnData,
    matrix: SparseMatrixReview,
    library_rows: list[dict[str, Any]],
    flags: Mapping[str, np.ndarray],
    n_mt_genes: int,
    join_mode: str,
    stat: Any,
) -> dict[str, Any]:
    total_counts = matrix.total_counts
    pct_mt = np.divide(
        matrix.mitochondrial_counts * 100.0,
        total_counts,
        out=np.zeros_like(total_counts),
        where=total_counts > 0,
    )
    tcr_positive = flags["has_any_ab_tcr"] | flags["has_any_gd_tcr"]
    join_expression = sum(
        int(row["join_expression_cells"])
        for row in library_rows
        if not pd.isna(row["join_expression_cells"])
    )
    join_matched = sum(
        int(row["join_matched_cells"])
        for row in library_rows
        if not pd.isna(row["join_matched_cells"])
    )
    unmatched_tcr = sum(
        int(row["unmatched_tcr_cells"])
        for row in library_rows
        if not pd.isna(row["unmatched_tcr_cells"])
    )
    return {
        "cohort_id": entry.cohort_id,
        "h5ad_path": str(entry.h5ad_path),
        "passed_structural_qc": True,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_libraries": len(library_rows),
        "n_samples": int(adata.obs["sample_id"].nunique()),
        "n_donors": int(adata.obs["donor_id"].nunique()),
        "sparse_encoding": matrix.encoding,
        "count_dtype": matrix.dtype,
        "matrix_nnz": matrix.nnz,
        "explicit_zero_entries": matrix.explicit_zero_entries,
        "count_min": matrix.count_min,
        "count_max": matrix.count_max,
        "total_count_sum": matrix.count_sum,
        "median_total_counts": float(np.median(total_counts)),
        "median_detected_genes": float(np.median(matrix.detected_genes)),
        "median_pct_counts_mt": float(np.median(pct_mt)),
        "n_mitochondrial_genes": n_mt_genes,
        "zero_count_cells": int(np.count_nonzero(total_counts == 0)),
        "cells_lt_200_genes": int(np.count_nonzero(matrix.detected_genes < 200)),
        "cells_gt_20pct_mt": int(np.count_nonzero(pct_mt > 20)),
        "tcr_positive_cells": int(tcr_positive.sum()),
        "has_TRA_cells": int(flags["has_TRA"].sum()),
        "has_TRB_cells": int(flags["has_TRB"].sum()),
        "has_TRG_cells": int(flags["has_TRG"].sum()),
        "has_TRD_cells": int(flags["has_TRD"].sum()),
        "paired_TRA_TRB_cells": int(flags["has_TRA_TRB_paired"].sum()),
        "paired_TRG_TRD_cells": int(flags["has_TRG_TRD_paired"].sum()),
        "tcr_positive_rate": float(tcr_positive.mean()),
        "paired_TRA_TRB_rate": float(flags["has_TRA_TRB_paired"].mean()),
        "paired_TRG_TRD_rate": float(flags["has_TRG_TRD_paired"].mean()),
        "join_mode": join_mode,
        "join_expression_cells": join_expression if join_expression else math.nan,
        "join_matched_cells": join_matched if join_expression else math.nan,
        "unmatched_tcr_cells": unmatched_tcr if join_expression else math.nan,
        "join_coverage": join_matched / join_expression if join_expression else math.nan,
        "source_size_bytes": int(stat.st_size),
        "source_mtime_ns": int(stat.st_mtime_ns),
    }


def audit_h5ad(
    entry: ManifestEntry,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, str]]]:
    before = entry.h5ad_path.stat()
    issues: list[dict[str, str]] = []
    adata = ad.read_h5ad(entry.h5ad_path, backed="r")
    try:
        if adata.n_obs <= 0 or adata.n_vars <= 0:
            raise Phase0QCError("H5AD must contain at least one cell and one gene")
        if not adata.obs_names.is_unique:
            raise Phase0QCError("obs_names are not unique")
        if not adata.var_names.is_unique:
            raise Phase0QCError("var_names are not unique")
        if clean_text_series(pd.Series(adata.obs_names)).eq("").any():
            raise Phase0QCError("obs_names contain blank values")
        if clean_text_series(pd.Series(adata.var_names)).eq("").any():
            raise Phase0QCError("var_names contain blank values")

        obs = adata.obs
        cleaned = validate_required_obs(obs, entry.cohort_id)
        validate_join_keys(cleaned)
        validate_library_metadata(cleaned)
        flags = validate_tcr_logic(obs)

        var_upper = pd.Index(adata.var_names.astype(str)).str.upper()
        mt_mask = np.asarray(var_upper.str.startswith("MT-"), dtype=bool)
        matrix = inspect_sparse_x(
            entry.h5ad_path,
            (int(adata.n_obs), int(adata.n_vars)),
            mt_mask,
        )
        if entry.expected_cells is not None and entry.expected_cells != adata.n_obs:
            raise Phase0QCError(
                f"Manifest expected {entry.expected_cells} cells but H5AD has {adata.n_obs}"
            )
        join_frame, join_mode = parse_join_summaries(adata.uns, cleaned, flags)
        library_rows = make_library_rows(
            entry.cohort_id,
            cleaned,
            flags,
            matrix,
            join_frame,
            join_mode,
        )
        summary = cohort_from_libraries(
            entry,
            adata,
            matrix,
            library_rows,
            flags,
            int(mt_mask.sum()),
            join_mode,
            before,
        )

        if matrix.explicit_zero_entries:
            issues.append(
                issue(
                    entry.cohort_id,
                    "WARNING",
                    "sparse_explicit_zeros",
                    f"Sparse X stores {matrix.explicit_zero_entries} explicit zero entries",
                )
            )
        if not mt_mask.any():
            issues.append(
                issue(
                    entry.cohort_id,
                    "WARNING",
                    "mitochondrial_gene_detection",
                    "No MT- gene symbols were available; mitochondrial QC is zero/unavailable",
                )
            )
        if summary["zero_count_cells"]:
            issues.append(
                issue(
                    entry.cohort_id,
                    "WARNING",
                    "zero_count_cells",
                    f"{summary['zero_count_cells']} cells have zero total counts",
                )
            )
        if summary["tcr_positive_cells"] == 0:
            issues.append(
                issue(
                    entry.cohort_id,
                    "WARNING",
                    "tcr_yield",
                    "No cell has canonical productive TCR evidence",
                )
            )
        if summary["unmatched_tcr_cells"] and not pd.isna(
            summary["unmatched_tcr_cells"]
        ):
            issues.append(
                issue(
                    entry.cohort_id,
                    "WARNING",
                    "unmatched_tcr_cells",
                    f"{int(summary['unmatched_tcr_cells'])} productive TCR key(s) did not match GEX",
                )
            )
    finally:
        adata.file.close()

    after = entry.h5ad_path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise Phase0QCError("Source H5AD size or modification time changed during read-only QC")
    return summary, library_rows, issues


def empty_failure_summary(entry: ManifestEntry, message: str) -> dict[str, Any]:
    stat = entry.h5ad_path.stat()
    return {
        "cohort_id": entry.cohort_id,
        "h5ad_path": str(entry.h5ad_path),
        "passed_structural_qc": False,
        "failure_reason": message,
        "source_size_bytes": int(stat.st_size),
        "source_mtime_ns": int(stat.st_mtime_ns),
    }


def records_for_json(frame: pd.DataFrame) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for record in frame.to_dict("records"):
        cleaned: dict[str, Any] = {}
        for key, value in record.items():
            if value is None or (isinstance(value, float) and math.isnan(value)):
                cleaned[key] = None
            elif isinstance(value, np.generic):
                cleaned[key] = value.item()
            else:
                cleaned[key] = value
        records.append(cleaned)
    return records


def markdown_table(frame: pd.DataFrame, columns: Sequence[str]) -> str:
    if frame.empty:
        return "_No rows._"
    available = [column for column in columns if column in frame]
    header = "| " + " | ".join(available) + " |"
    separator = "| " + " | ".join("---" for _ in available) + " |"
    rows = [header, separator]
    for _, record in frame[available].iterrows():
        values: list[str] = []
        for value in record:
            if pd.isna(value):
                text = "NA"
            elif isinstance(value, float):
                text = f"{value:.4g}"
            else:
                text = str(value)
            values.append(text.replace("|", "\\|"))
        rows.append("| " + " | ".join(values) + " |")
    return "\n".join(rows)


def save_figure(fig: Figure, path: Path) -> None:
    FigureCanvasAgg(fig)
    fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white")


def placeholder_figure(title: str, message: str, path: Path) -> None:
    fig = Figure(figsize=(8, 3.5), constrained_layout=True)
    axis = fig.subplots()
    axis.axis("off")
    axis.text(0.5, 0.62, title, ha="center", va="center", fontsize=15, weight="bold")
    axis.text(0.5, 0.38, message, ha="center", va="center", fontsize=11)
    save_figure(fig, path)


def plot_cohort_overview(frame: pd.DataFrame, path: Path) -> None:
    valid = frame.loc[frame.get("passed_structural_qc", False).astype(bool)].copy()
    if valid.empty:
        placeholder_figure(
            "Extension Phase 0 cohort review",
            "No cohort passed structural QC.",
            path,
        )
        return
    valid = valid.sort_values("n_cells", ascending=True)
    labels = valid["cohort_id"].astype(str).tolist()
    y = np.arange(len(valid))
    height = max(5.5, 0.48 * len(valid) + 2.5)
    fig = Figure(figsize=(12, height), constrained_layout=True)
    axes = fig.subplots(2, 2)
    panels = (
        ("n_cells", "Cells", "#277da1", False),
        ("median_detected_genes", "Median detected genes", "#43aa8b", False),
        ("tcr_positive_rate", "Productive TCR-positive cells (%)", "#f8961e", True),
        ("paired_TRG_TRD_rate", "Paired TRG/TRD cells (%)", "#9b5de5", True),
    )
    for axis, (column, title, color, percent) in zip(axes.ravel(), panels, strict=True):
        values = pd.to_numeric(valid[column], errors="coerce").fillna(0).to_numpy()
        if percent:
            values = values * 100.0
        axis.barh(y, values, color=color, height=0.68)
        axis.set_yticks(y, labels if axis in axes[:, 0] else [""] * len(labels))
        axis.set_title(title, fontsize=11)
        axis.grid(axis="x", color="#d9d9d9", linewidth=0.7)
        axis.set_axisbelow(True)
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Extension cohort standalone Phase 0 review", fontsize=15, weight="bold")
    save_figure(fig, path)


def plot_library_qc(frame: pd.DataFrame, path: Path) -> None:
    if frame.empty:
        placeholder_figure(
            "Extension Phase 0 library QC",
            "No library-level metrics were available.",
            path,
        )
        return
    frame = frame.copy().sort_values(["cohort_id", "library_id"], ascending=False)
    labels = (frame["cohort_id"].astype(str) + " | " + frame["library_id"].astype(str)).tolist()
    y = np.arange(len(frame))
    height = min(24, max(6, 0.35 * len(frame) + 3))
    fig = Figure(figsize=(13, height), constrained_layout=True)
    axes = fig.subplots(2, 2)
    panels = (
        ("median_total_counts", "Median total counts", "#277da1"),
        ("median_detected_genes", "Median detected genes", "#43aa8b"),
        ("median_pct_counts_mt", "Median mitochondrial counts (%)", "#f94144"),
        ("tcr_positive_rate", "Productive TCR-positive cells (%)", "#f8961e"),
    )
    for axis, (column, title, color) in zip(axes.ravel(), panels, strict=True):
        values = pd.to_numeric(frame[column], errors="coerce").fillna(0).to_numpy()
        if column.endswith("_rate"):
            values = values * 100.0
        axis.scatter(values, y, color=color, s=24, zorder=3)
        axis.set_yticks(y, labels if axis in axes[:, 0] else [""] * len(labels))
        axis.set_title(title, fontsize=11)
        axis.grid(axis="x", color="#d9d9d9", linewidth=0.7)
        axis.set_axisbelow(True)
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Extension cohort library-level QC", fontsize=15, weight="bold")
    save_figure(fig, path)


def write_review_artifacts(
    output_root: Path,
    manifest_path: Path,
    manifest_audit: pd.DataFrame,
    cohort_frame: pd.DataFrame,
    library_frame: pd.DataFrame,
    issue_frame: pd.DataFrame,
) -> dict[str, str]:
    table_dir = output_root / "tables" / REVIEW_SUBDIR
    log_dir = output_root / "logs" / REVIEW_SUBDIR
    figure_dir = output_root / "figures" / REVIEW_SUBDIR
    for directory in (table_dir, log_dir, figure_dir):
        directory.mkdir(parents=True, exist_ok=True)

    artifacts = {
        "manifest_audit_csv": str(table_dir / "extension_phase0_manifest_audit.csv"),
        "cohort_summary_csv": str(table_dir / "extension_phase0_cohort_summary.csv"),
        "library_summary_csv": str(table_dir / "extension_phase0_library_summary.csv"),
        "issues_csv": str(table_dir / "extension_phase0_validation_issues.csv"),
        "summary_json": str(log_dir / "extension_phase0_qc_summary.json"),
        "summary_markdown": str(log_dir / "extension_phase0_qc_summary.md"),
        "cohort_overview_png": str(figure_dir / "extension_phase0_cohort_overview.png"),
        "library_qc_png": str(figure_dir / "extension_phase0_library_qc.png"),
    }
    manifest_audit.to_csv(artifacts["manifest_audit_csv"], index=False)
    cohort_frame.to_csv(artifacts["cohort_summary_csv"], index=False)
    library_frame.to_csv(artifacts["library_summary_csv"], index=False)
    issue_frame.to_csv(artifacts["issues_csv"], index=False)

    errors = int(issue_frame["severity"].eq("ERROR").sum()) if not issue_frame.empty else 0
    warnings = (
        int(issue_frame["severity"].eq("WARNING").sum()) if not issue_frame.empty else 0
    )
    gate_status = "FAIL" if errors else "PASS_REVIEW_REQUIRED"
    generated_at = datetime.now(timezone.utc).isoformat()
    cell_count_values = (
        pd.to_numeric(cohort_frame["n_cells"], errors="coerce").fillna(0)
        if "n_cells" in cohort_frame
        else pd.Series(dtype=np.int64)
    )
    payload = {
        "schema_version": 1,
        "generated_at_utc": generated_at,
        "manifest_path": str(manifest_path),
        "manifest_sha256": sha256_file(manifest_path),
        "gate_status": gate_status,
        "review_required": True,
        "merge_approved": False,
        "source_h5ads_opened_read_only": True,
        "source_h5ads_mutated": False,
        "cohort_count": int(len(cohort_frame)),
        "passing_cohort_count": int(
            cohort_frame.get("passed_structural_qc", pd.Series(dtype=bool)).fillna(False).sum()
        ),
        "library_count": int(len(library_frame)),
        "cell_count": int(cell_count_values.sum()),
        "error_count": errors,
        "warning_count": warnings,
        "artifacts": artifacts,
        "cohorts": records_for_json(cohort_frame),
        "issues": records_for_json(issue_frame),
    }
    Path(artifacts["summary_json"]).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    lines = [
        "# Extension H5AD Standalone Phase 0 QC",
        "",
        f"- Generated (UTC): `{generated_at}`",
        f"- Manifest: `{manifest_path}`",
        f"- Gate status: **{gate_status}**",
        f"- Cohorts passing structural QC: `{payload['passing_cohort_count']}/{payload['cohort_count']}`",
        f"- Audited cells: `{payload['cell_count']:,}`",
        f"- Structural errors: `{errors}`",
        f"- Review warnings: `{warnings}`",
        "- Merge approved: **No**",
        "",
        "> A structural PASS means only that the standalone files are eligible for user review. "
        "It does not authorize merging, integration, model training, or source-file mutation.",
        "",
        "## Cohort Summary",
        "",
        markdown_table(
            cohort_frame,
            (
                "cohort_id",
                "passed_structural_qc",
                "n_cells",
                "n_libraries",
                "n_donors",
                "median_detected_genes",
                "median_pct_counts_mt",
                "tcr_positive_cells",
                "paired_TRA_TRB_cells",
                "paired_TRG_TRD_cells",
                "join_coverage",
            ),
        ),
        "",
        "## Validation Issues",
        "",
        markdown_table(
            issue_frame,
            ("severity", "cohort_id", "library_id", "check", "message"),
        ),
        "",
        "## Gate Semantics",
        "",
        "- Every H5AD was opened with read-only backing; sparse `X` values were streamed in bounded chunks.",
        "- Structural errors fail the gate and produce a nonzero process exit code.",
        "- Biological QC observations are warnings for manual review and never trigger automatic filtering.",
        "- The workflow never merges cohorts and never writes to a source H5AD.",
    ]
    Path(artifacts["summary_markdown"]).write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    plot_cohort_overview(cohort_frame, Path(artifacts["cohort_overview_png"]))
    plot_library_qc(library_frame, Path(artifacts["library_qc_png"]))
    return artifacts


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run read-only standalone Phase 0 QC on built extension H5ADs."
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=DEFAULT_MANIFEST,
        help=f"Built-H5AD manifest (default: {DEFAULT_MANIFEST})",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Integrated_dataset-compatible output root containing tables/logs/figures",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate the manifest and source paths, print the plan, and write nothing",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest_path = args.manifest.expanduser().resolve()
    output_root = args.output_root.expanduser().resolve()
    try:
        entries, manifest_audit = load_manifest(manifest_path)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    if args.dry_run:
        print(
            json.dumps(
                {
                    "dry_run": True,
                    "manifest": str(manifest_path),
                    "output_root": str(output_root),
                    "cohort_count": len(entries),
                    "inputs": [
                        {"cohort_id": entry.cohort_id, "h5ad_path": str(entry.h5ad_path)}
                        for entry in entries
                    ],
                    "writes_planned": False,
                    "merge_planned": False,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0

    cohort_rows: list[dict[str, Any]] = []
    library_rows: list[dict[str, Any]] = []
    issues: list[dict[str, str]] = []
    for entry in entries:
        print(f"Auditing {entry.cohort_id}: {entry.h5ad_path}", flush=True)
        try:
            cohort, libraries, cohort_issues = audit_h5ad(entry)
            cohort_rows.append(cohort)
            library_rows.extend(libraries)
            issues.extend(cohort_issues)
        except Exception as exc:
            message = str(exc)
            cohort_rows.append(empty_failure_summary(entry, message))
            issues.append(issue(entry.cohort_id, "ERROR", "structural_qc", message))

    cohort_frame = pd.DataFrame(cohort_rows).sort_values("cohort_id").reset_index(drop=True)
    library_frame = pd.DataFrame(library_rows)
    if not library_frame.empty:
        library_frame = library_frame.sort_values(["cohort_id", "library_id"]).reset_index(
            drop=True
        )
    issue_frame = pd.DataFrame(
        issues,
        columns=("cohort_id", "library_id", "severity", "check", "message"),
    )
    if not issue_frame.empty:
        issue_frame = issue_frame.sort_values(
            ["severity", "cohort_id", "library_id", "check"]
        ).reset_index(drop=True)

    artifacts = write_review_artifacts(
        output_root,
        manifest_path,
        manifest_audit,
        cohort_frame,
        library_frame,
        issue_frame,
    )
    error_count = int(issue_frame["severity"].eq("ERROR").sum()) if not issue_frame.empty else 0
    print(f"Wrote Phase 0 review: {artifacts['summary_markdown']}")
    if error_count:
        print(f"FAIL: {error_count} structural error(s); merge remains blocked", file=sys.stderr)
        return 2
    print("PASS_REVIEW_REQUIRED: structural QC passed; merge remains unapproved")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
