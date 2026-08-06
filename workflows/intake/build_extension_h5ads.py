#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Build provenance-preserving extension H5ADs from immutable staged sources."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import tempfile
import uuid
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

from validate_extension_cohorts import (
    DEFAULT_CANONICAL_REGISTRY,
    DEFAULT_COHORTS,
    DEFAULT_LIBRARIES,
    ManifestError,
    assert_extension_output_root,
    parse_bool,
    split_tokens,
    validate_built_h5ad,
    validate_extension_manifests,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_EXTENSION_ROOT = PROJECT_ROOT / "data/interim/extension_intake"
DEFAULT_STAGE_ROOT = DEFAULT_EXTENSION_ROOT / "staged"
DEFAULT_OUTPUT_ROOT = DEFAULT_EXTENSION_ROOT / "built_h5ads"
DEFAULT_BUILT_MANIFEST = DEFAULT_EXTENSION_ROOT / "built_h5ads_manifest.csv"
DEFAULT_R_EXPORTER = Path(__file__).with_name("export_gse159251_rds.R")

CHAINS = ("TRA", "TRB", "TRG", "TRD")
CHAIN_SUFFIXES = ("cdr3", "v", "d", "j", "cdr3_nt", "clone_id", "umis", "reads")
STRING_SUFFIXES = ("cdr3", "v", "d", "j", "cdr3_nt", "clone_id")
NUMERIC_SUFFIXES = ("umis", "reads")
TCR_COLUMNS = tuple(f"{chain}_{suffix}" for chain in CHAINS for suffix in CHAIN_SUFFIXES)
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
REQUIRED_CELL_METADATA = ("sample_id", "library_id", "donor_id", "tissue", "diagnosis")
SOURCE_METASHEET_REQUIRED = ("gse", "gsm", "sample_title", "modality")
TRUE_VALUES = {"true", "t", "yes", "y", "1"}
FALSE_VALUES = {"false", "f", "no", "n", "0", ""}


class IntakeBuildError(RuntimeError):
    """Raised when an extension build would lose or misattribute data."""


def clean_text(value: Any) -> str:
    """Normalize scalar metadata without turning missing values into text."""
    if value is None or (not isinstance(value, (list, dict)) and pd.isna(value)):
        return ""
    text = str(value).strip()
    return "" if text.casefold() in {"", "na", "nan", "none", "null", "<na>", "n/a"} else text


def clean_series(series: pd.Series) -> pd.Series:
    return series.map(clean_text).astype(object)


def normalize_bool(value: Any, *, field: str = "boolean") -> bool:
    text = clean_text(value).casefold()
    if text in TRUE_VALUES:
        return True
    if text in FALSE_VALUES:
        return False
    raise IntakeBuildError(f"Unrecognized {field} value: {value!r}")


def extract_full_barcode(value: Any) -> str:
    """Extract a full 10x barcode when it is embedded in a composite cell ID."""
    text = clean_text(value).upper()
    match = re.search(r"([ACGTN]{8,}(?:-\d+)?)", text)
    return match.group(1) if match else text


def barcode_core(value: Any) -> str:
    """Return the upper-case 10x sequence before lane/sample suffixes."""
    full = extract_full_barcode(value)
    match = re.match(r"([ACGTN]{8,})", full)
    if match:
        return match.group(1)
    return full.split("-")[0].split(".")[0].split(":")[-1]


def ensure_unique_key(frame: pd.DataFrame, columns: Sequence[str], label: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise IntakeBuildError(f"{label} lacks key column(s): {missing}")
    blank = frame[list(columns)].apply(lambda col: clean_series(col).eq("")).any(axis=1)
    if blank.any():
        raise IntakeBuildError(f"{label} has {int(blank.sum())} blank join-key row(s)")
    duplicated = frame.duplicated(list(columns), keep=False)
    if duplicated.any():
        examples = frame.loc[duplicated, list(columns)].head(5).to_dict("records")
        raise IntakeBuildError(
            f"{label} has {int(duplicated.sum())} rows with duplicate key {list(columns)}; "
            f"examples={examples}"
        )


def validate_sparse_integer_counts(matrix: Any, label: str) -> sparse.csr_matrix:
    """Require nonnegative, finite, integer-valued sparse counts without densifying."""
    if not sparse.issparse(matrix):
        raise IntakeBuildError(f"{label} is dense; sparse raw counts are required")
    result = matrix.tocsr(copy=False)
    values = result.data
    if values.size:
        if not np.isfinite(values).all():
            raise IntakeBuildError(f"{label} contains non-finite values")
        if np.any(values < 0):
            raise IntakeBuildError(f"{label} contains negative values")
        if not np.issubdtype(values.dtype, np.integer) and not np.equal(
            values, np.floor(values)
        ).all():
            raise IntakeBuildError(f"{label} is not integer-valued raw counts")
    return result


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", header=None, compression="infer", dtype=str)


def _pick_one(candidates: Iterable[Path], label: str) -> Path:
    unique = sorted({path.resolve() for path in candidates if path.is_file()}, key=str)
    if len(unique) != 1:
        raise IntakeBuildError(f"Expected one {label}, found {len(unique)}: {[str(p) for p in unique]}")
    return unique[0]


def _matrix_bundle_files(path: Path) -> tuple[Path, Path, Path]:
    if path.is_dir():
        matrix = _pick_one(
            list(path.glob("*matrix.mtx*")) + list(path.glob("*counts.mtx*")),
            "Matrix Market count file",
        )
        directory = path
    else:
        matrix = path.resolve()
        directory = matrix.parent
    if not matrix.is_file():
        raise FileNotFoundError(matrix)

    prefix = re.sub(r"(?:matrix|counts)\.mtx(?:\.gz)?$", "", matrix.name, flags=re.IGNORECASE)
    barcode_candidates = list(directory.glob("*barcodes.tsv*")) + list(directory.glob("*barcode.tsv*"))
    feature_candidates = (
        list(directory.glob("*features.tsv*"))
        + list(directory.glob("*feature.tsv*"))
        + list(directory.glob("*genes.tsv*"))
        + list(directory.glob("*gene.tsv*"))
    )
    prefixed_barcodes = [candidate for candidate in barcode_candidates if candidate.name.startswith(prefix)]
    prefixed_features = [candidate for candidate in feature_candidates if candidate.name.startswith(prefix)]
    barcodes = _pick_one(prefixed_barcodes or barcode_candidates, "barcode file")
    features = _pick_one(prefixed_features or feature_candidates, "feature file")
    return matrix, barcodes, features


def read_mtx_counts(path: Path) -> ad.AnnData:
    """Read one genes-by-cells or cells-by-genes Matrix Market bundle sparsely."""
    matrix_path, barcode_path, feature_path = _matrix_bundle_files(path)
    matrix = mmread(matrix_path)
    if not sparse.issparse(matrix):
        raise IntakeBuildError(f"Matrix Market source unexpectedly decoded densely: {matrix_path}")
    barcodes = _read_tsv(barcode_path).iloc[:, 0].map(clean_text).tolist()
    features = _read_tsv(feature_path)
    feature_ids = features.iloc[:, 0].map(clean_text)
    feature_names = features.iloc[:, 1].map(clean_text) if features.shape[1] > 1 else feature_ids

    if matrix.shape == (len(features), len(barcodes)):
        matrix = matrix.T
    elif matrix.shape != (len(barcodes), len(features)):
        raise IntakeBuildError(
            f"Count shape {matrix.shape} does not match {len(barcodes)} barcodes and {len(features)} features"
        )
    matrix = validate_sparse_integer_counts(matrix, str(matrix_path))
    obs = pd.DataFrame(index=pd.Index(barcodes, name="cell_id"))
    var = pd.DataFrame(
        {"feature_id": feature_ids.to_numpy(dtype=object)},
        index=pd.Index(feature_names, name="gene_symbol"),
    )
    if features.shape[1] > 2:
        var["feature_type"] = features.iloc[:, 2].map(clean_text).to_numpy(dtype=object)
    result = ad.AnnData(X=matrix, obs=obs, var=var)
    result.var_names_make_unique()
    result.uns["count_source"] = str(matrix_path)
    result.uns["count_matrix_state"] = "raw_integer_counts"
    return result


def _counts_from_h5ad(source: ad.AnnData, path: Path) -> tuple[sparse.csr_matrix, pd.DataFrame]:
    if "counts" in source.layers:
        matrix = source.layers["counts"]
        var = source.var.copy()
        location = "layers/counts"
    elif source.raw is not None:
        matrix = source.raw.X
        var = source.raw.var.copy()
        location = "raw/X"
    else:
        matrix = source.X
        var = source.var.copy()
        location = "X"
    counts = validate_sparse_integer_counts(matrix, f"{path}:{location}")
    return counts, var


def read_expression_source(path: Path) -> ad.AnnData:
    """Read a staged MTX, 10x H5, or H5AD source and retain sparse counts in X."""
    path = path.expanduser().resolve()
    if path.is_dir() or ".mtx" in path.name.casefold():
        return read_mtx_counts(path)
    if path.suffix.casefold() == ".h5ad":
        source = ad.read_h5ad(path)
        counts, var = _counts_from_h5ad(source, path)
        result = ad.AnnData(X=counts.copy(), obs=source.obs.copy(), var=var)
        result.var_names_make_unique()
        result.uns["count_source"] = str(path)
        result.uns["count_matrix_state"] = "raw_integer_counts"
        return result
    if path.suffix.casefold() in {".h5", ".hdf5"}:
        import scanpy as sc

        result = sc.read_10x_h5(path, gex_only=True)
        result.X = validate_sparse_integer_counts(result.X, str(path))
        result.var_names_make_unique()
        result.uns["count_source"] = str(path)
        result.uns["count_matrix_state"] = "raw_integer_counts"
        return result
    raise IntakeBuildError(f"Unsupported expression source: {path}")


def read_table(path: Path, **kwargs: Any) -> pd.DataFrame:
    name = path.name.casefold()
    if name.endswith((".xlsx", ".xls")):
        return pd.read_excel(path, **kwargs)
    separator = "\t" if name.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")) else ","
    return pd.read_csv(path, sep=separator, compression="infer", low_memory=False, **kwargs)


def empty_tcr_frame(index: pd.Index | None = None) -> pd.DataFrame:
    result = pd.DataFrame(index=index)
    for chain in CHAINS:
        for suffix in STRING_SUFFIXES:
            result[f"{chain}_{suffix}"] = ""
        for suffix in NUMERIC_SUFFIXES:
            result[f"{chain}_{suffix}"] = 0
    return result


def _first_column(frame: pd.DataFrame, aliases: Sequence[str], *, required: bool = False) -> pd.Series:
    by_normalized = {re.sub(r"[^a-z0-9]+", "", str(column).casefold()): column for column in frame.columns}
    for alias in aliases:
        key = re.sub(r"[^a-z0-9]+", "", alias.casefold())
        if key in by_normalized:
            return frame[by_normalized[key]]
    if required:
        raise IntakeBuildError(f"None of required columns {list(aliases)} found in {list(frame.columns)}")
    return pd.Series([""] * len(frame), index=frame.index, dtype=object)


def collapse_productive_contigs(frame: pd.DataFrame, sample_id: str | None = None) -> pd.DataFrame:
    """Select one productive contig per sample/barcode/chain and pivot canonically."""
    if frame.empty:
        return pd.DataFrame(columns=["sample_id", "barcode_core", *TCR_COLUMNS])

    sample = _first_column(frame, ("sample_id", "sample", "library_id"))
    if sample.map(clean_text).eq("").all():
        if not sample_id:
            raise IntakeBuildError("Contig table has no sample_id and no sample override")
        sample = pd.Series([sample_id] * len(frame), index=frame.index)
    barcode = _first_column(frame, ("barcode", "cell_barcode", "barcode_core"), required=True)
    chain = _first_column(frame, ("chain", "locus"), required=True).map(clean_text).str.upper()
    productive_raw = _first_column(frame, ("productive", "is_productive"), required=True)
    productive = productive_raw.map(lambda value: normalize_bool(value, field="productive"))
    is_cell_raw = _first_column(frame, ("is_cell",))
    high_conf_raw = _first_column(frame, ("high_confidence",))
    is_cell = is_cell_raw.map(lambda value: normalize_bool(value, field="is_cell"))
    high_conf = high_conf_raw.map(lambda value: normalize_bool(value, field="high_confidence"))

    work = pd.DataFrame(
        {
            "sample_id": sample.map(clean_text),
            "barcode_core": barcode.map(barcode_core),
            "chain": chain,
            "cdr3": _first_column(frame, ("cdr3", "junction_aa")).map(clean_text),
            "v": _first_column(frame, ("v_gene", "v_call")).map(clean_text),
            "d": _first_column(frame, ("d_gene", "d_call")).map(clean_text),
            "j": _first_column(frame, ("j_gene", "j_call")).map(clean_text),
            "cdr3_nt": _first_column(frame, ("cdr3_nt", "junction")).map(clean_text),
            "clone_id": _first_column(
                frame, ("raw_clonotype_id", "clone_id", "clonotype_id")
            ).map(clean_text),
            "umis": pd.to_numeric(
                _first_column(frame, ("umis", "umi_count", "duplicate_count")), errors="coerce"
            ).fillna(0).astype(np.int64),
            "reads": pd.to_numeric(
                _first_column(frame, ("reads", "read_count", "consensus_count")), errors="coerce"
            ).fillna(0).astype(np.int64),
            "productive": productive,
            "is_cell": is_cell if not is_cell_raw.map(clean_text).eq("").all() else True,
            "high_confidence": high_conf if not high_conf_raw.map(clean_text).eq("").all() else True,
        }
    )
    work = work.loc[
        work["productive"]
        & work["is_cell"]
        & work["high_confidence"]
        & work["chain"].isin(CHAINS)
        & work["sample_id"].ne("")
        & work["barcode_core"].ne("")
        & work["cdr3"].ne("")
    ].copy()
    if work.empty:
        return pd.DataFrame(columns=["sample_id", "barcode_core", *TCR_COLUMNS])

    exact_columns = [
        "sample_id",
        "barcode_core",
        "chain",
        "cdr3",
        "v",
        "d",
        "j",
        "cdr3_nt",
        "clone_id",
        "umis",
        "reads",
    ]
    exact_duplicates = work.duplicated(exact_columns, keep=False)
    if exact_duplicates.any():
        raise IntakeBuildError(
            f"Contig table contains {int(exact_duplicates.sum())} exact duplicate productive contigs"
        )

    work = work.sort_values(
        ["sample_id", "barcode_core", "chain", "umis", "reads", "cdr3"],
        ascending=[True, True, True, False, False, True],
        kind="mergesort",
    )
    selected = work.groupby(["sample_id", "barcode_core", "chain"], sort=False).head(1)
    keys = selected[["sample_id", "barcode_core"]].drop_duplicates()
    result = keys.copy()
    for chain_name in CHAINS:
        subset = selected.loc[selected["chain"] == chain_name, [
            "sample_id", "barcode_core", *CHAIN_SUFFIXES
        ]].copy()
        subset = subset.rename(
            columns={suffix: f"{chain_name}_{suffix}" for suffix in CHAIN_SUFFIXES}
        )
        result = result.merge(subset, on=["sample_id", "barcode_core"], how="left", validate="one_to_one")
    for chain_name in CHAINS:
        for suffix in STRING_SUFFIXES:
            column = f"{chain_name}_{suffix}"
            result[column] = clean_series(result.get(column, pd.Series("", index=result.index)))
        for suffix in NUMERIC_SUFFIXES:
            column = f"{chain_name}_{suffix}"
            result[column] = pd.to_numeric(result.get(column, 0), errors="coerce").fillna(0).astype(np.int64)
    ensure_unique_key(result, ("sample_id", "barcode_core"), "productive TCR table")
    return result[["sample_id", "barcode_core", *TCR_COLUMNS]]


def standardize_tcr_table(path: Path, sample_id: str | None = None) -> pd.DataFrame:
    """Read either a canonical wide TCR table or a productive contig table."""
    frame = read_table(path)
    canonical_present = any(column in frame.columns for column in TCR_COLUMNS)
    if not canonical_present:
        return collapse_productive_contigs(frame, sample_id)

    result = pd.DataFrame(index=frame.index)
    sample = _first_column(frame, ("sample_id", "sample", "library_id"))
    if sample.map(clean_text).eq("").all():
        if not sample_id:
            raise IntakeBuildError("Canonical TCR table lacks sample_id")
        sample = pd.Series([sample_id] * len(frame), index=frame.index)
    result["sample_id"] = sample.map(clean_text)
    result["barcode_core"] = _first_column(
        frame, ("barcode_core", "barcode", "cell_barcode"), required=True
    ).map(barcode_core)
    for chain in CHAINS:
        for suffix in STRING_SUFFIXES:
            column = f"{chain}_{suffix}"
            result[column] = clean_series(frame[column]) if column in frame else ""
        for suffix in NUMERIC_SUFFIXES:
            column = f"{chain}_{suffix}"
            result[column] = (
                pd.to_numeric(frame[column], errors="coerce").fillna(0).astype(np.int64)
                if column in frame
                else 0
            )
    ensure_unique_key(result, ("sample_id", "barcode_core"), "canonical TCR table")
    return result[["sample_id", "barcode_core", *TCR_COLUMNS]]


def ensure_canonical_tcr(obs: pd.DataFrame) -> pd.DataFrame:
    result = obs.copy()
    for chain in CHAINS:
        for suffix in STRING_SUFFIXES:
            column = f"{chain}_{suffix}"
            result[column] = clean_series(result[column]) if column in result else ""
        for suffix in NUMERIC_SUFFIXES:
            column = f"{chain}_{suffix}"
            result[column] = (
                pd.to_numeric(result[column], errors="coerce").fillna(0).astype(np.int64)
                if column in result
                else 0
            )
        result[f"has_{chain}"] = result[f"{chain}_cdr3"].ne("")
    result["has_TRA_TRB_paired"] = result["has_TRA"] & result["has_TRB"]
    result["has_TRG_TRD_paired"] = result["has_TRG"] & result["has_TRD"]
    result["has_any_ab_tcr"] = result["has_TRA"] | result["has_TRB"]
    result["has_any_gd_tcr"] = result["has_TRG"] | result["has_TRD"]
    result["TCRseq"] = np.where(result["has_any_ab_tcr"] | result["has_any_gd_tcr"], "yes", "no")
    return result


def join_tcr_by_sample_barcode(adata: ad.AnnData, tcr: pd.DataFrame) -> dict[str, int]:
    """Left-join a unique canonical TCR table without dropping expression cells."""
    keys = adata.obs[["sample_id", "barcode_core"]].copy()
    keys["sample_id"] = clean_series(keys["sample_id"])
    keys["barcode_core"] = keys["barcode_core"].map(barcode_core)
    ensure_unique_key(keys, ("sample_id", "barcode_core"), "expression cells")
    ensure_unique_key(tcr, ("sample_id", "barcode_core"), "TCR join table")

    left = keys.reset_index(names="_obs_name")
    joined = left.merge(tcr, on=["sample_id", "barcode_core"], how="left", validate="one_to_one", indicator=True)
    if len(joined) != adata.n_obs:
        raise IntakeBuildError("TCR join changed the expression row count")
    joined = joined.set_index("_obs_name").reindex(adata.obs_names)
    for column in TCR_COLUMNS:
        adata.obs[column] = joined[column].to_numpy() if column in joined else ""
    matched = int(joined["_merge"].eq("both").sum())
    unmatched_tcr = int(
        tcr.merge(keys.drop_duplicates(), on=["sample_id", "barcode_core"], how="left", indicator=True)[
            "_merge"
        ].eq("left_only").sum()
    )
    adata.obs = ensure_canonical_tcr(adata.obs)
    return {"expression_cells": adata.n_obs, "matched_cells": matched, "unmatched_tcr_cells": unmatched_tcr}


def load_stage_manifest(stage_dir: Path) -> dict[str, Any]:
    path = stage_dir / "stage_manifest.json"
    if not path.is_file():
        raise FileNotFoundError(f"Missing stage manifest: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("cohort_id") != stage_dir.name:
        raise IntakeBuildError(
            f"Stage manifest cohort {payload.get('cohort_id')!r} does not match directory {stage_dir.name!r}"
        )
    if payload.get("join_key") != "sample_id+barcode_core":
        raise IntakeBuildError(f"Unsafe or missing stage join rule in {path}")
    return payload


def staged_sources(stage_dir: Path, role: str | None = None) -> list[Path]:
    payload = load_stage_manifest(stage_dir)
    paths: list[Path] = []
    for record in payload.get("sources", []):
        if role and record.get("source_role") != role:
            continue
        path = stage_dir / str(record.get("staged_path", ""))
        if not path.is_file():
            raise FileNotFoundError(f"Staged source is missing: {path}")
        paths.append(path)
    return sorted(set(paths), key=str)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def write_built_manifest(output_root: Path, built: Sequence[Mapping[str, Any]]) -> Path:
    """Atomically publish the immutable inputs consumed by the Phase 0 QC gate."""
    manifest_path = output_root.parent / DEFAULT_BUILT_MANIFEST.name
    if manifest_path.exists():
        raise FileExistsError(f"Refusing to replace built-H5AD manifest: {manifest_path}")
    rows: list[dict[str, Any]] = []
    for item in built:
        h5ad_path = Path(str(item["h5ad"])).resolve()
        qc_path = Path(str(item["standalone_qc"])).resolve()
        stat = h5ad_path.stat()
        rows.append(
            {
                "cohort_id": str(item["cohort_id"]),
                "h5ad_path": str(h5ad_path),
                "build_status": "built",
                "n_cells": int(item["cells"]),
                "n_genes": int(item["genes"]),
                "standalone_qc_path": str(qc_path),
                "size_bytes": int(stat.st_size),
                "mtime_ns": int(stat.st_mtime_ns),
                "sha256": _sha256(h5ad_path),
            }
        )
    temporary = manifest_path.with_name(f".{manifest_path.name}.{uuid.uuid4().hex}.tmp")
    try:
        pd.DataFrame(rows).sort_values("cohort_id").to_csv(temporary, index=False)
        temporary.replace(manifest_path)
    finally:
        temporary.unlink(missing_ok=True)
    return manifest_path


def validated_existing_build(
    output_path: Path,
    qc_path: Path,
    cohort: Mapping[str, str],
) -> dict[str, Any]:
    """Admit an interrupted-run output only after serialized and QC validation."""
    if not output_path.is_file() or not qc_path.is_file():
        raise IntakeBuildError(
            f"Resume requires both H5AD and standalone QC for {cohort['cohort_id']}"
        )
    issues = validate_built_h5ad(output_path, cohort)
    errors = [issue for issue in issues if issue.severity == "error"]
    if errors:
        details = "; ".join(f"{issue.code}: {issue.message}" for issue in errors)
        raise IntakeBuildError(f"Existing build failed validation for {cohort['cohort_id']}: {details}")
    summary = json.loads(qc_path.read_text(encoding="utf-8"))
    if not summary.get("passed"):
        raise IntakeBuildError(f"Existing standalone QC did not pass for {cohort['cohort_id']}")
    if str(summary.get("cohort_id", "")) != str(cohort["cohort_id"]):
        raise IntakeBuildError(f"Existing standalone QC cohort mismatch for {cohort['cohort_id']}")
    return {
        "cohort_id": str(cohort["cohort_id"]),
        "h5ad": str(output_path),
        "standalone_qc": str(qc_path),
        "cells": int(summary["cells"]),
        "genes": int(summary["genes"]),
        "resumed_existing": True,
    }


def read_source_metasheet(stage_dir: Path) -> pd.DataFrame:
    """Read and verify the accession-filtered shared source metadata snapshot."""
    path = stage_dir / "sample_metasheet.csv"
    if not path.is_file():
        raise FileNotFoundError(f"Filtered shared metasheet snapshot is missing: {path}")
    payload = load_stage_manifest(stage_dir)
    record = payload.get("metasheet") or {}
    expected_hash = clean_text(record.get("staged_sha256"))
    if expected_hash and _sha256(path) != expected_hash:
        raise IntakeBuildError(f"Filtered shared metasheet checksum mismatch: {path}")
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = [column for column in SOURCE_METASHEET_REQUIRED if column not in frame]
    if missing:
        raise IntakeBuildError(f"{path} lacks required shared metadata columns: {missing}")
    for column in frame.columns:
        frame[column] = clean_series(frame[column])
    ensure_unique_key(frame, ("gsm",), "filtered shared source metasheet")
    return frame


def biological_sample_key(title: Any) -> str:
    """Normalize paired GEX/VDJ titles to one biological-sample key."""
    text = clean_text(title)
    if ";" in text and re.search(r"sc(?:rna|vdj|tcr)", text, flags=re.IGNORECASE):
        text = text.rsplit(";", 1)[-1].strip()
    text = re.sub(r"\s*\[VDJ sequencing\]\s*$", "", text, flags=re.IGNORECASE)
    text = re.sub(r"\s*\(TCR\)\s*$", "", text, flags=re.IGNORECASE)
    text = re.sub(r"_(?:GEX|VDJ|TCR)$", "", text, flags=re.IGNORECASE)
    text = re.sub(r"\s+", "_", text.strip())
    text = re.sub(r"[^A-Za-z0-9_.+-]+", "_", text).strip("_")
    if not text:
        raise IntakeBuildError(f"Cannot derive biological sample key from title {title!r}")
    return text


def _primary_sources_by_gsm(stage_dir: Path, role: str) -> dict[str, Path]:
    sources: dict[str, Path] = {}
    for path in staged_sources(stage_dir, role):
        name = path.name.casefold()
        if role == "gex":
            primary = ".mtx" in name and any(token in name for token in ("matrix", "counts"))
            primary = primary or name.endswith((".h5", ".hdf5", ".h5ad"))
        elif role == "vdj":
            primary = "filtered_contig_annotations" in name or "contig_annotations" in name
        elif role == "combined":
            primary = name.endswith((".rds", ".rds.gz", ".rds.gz.gz", ".h5ad"))
        else:
            primary = True
        if not primary:
            continue
        match = re.search(r"GSM\d+", path.name, flags=re.IGNORECASE)
        if not match:
            raise IntakeBuildError(f"Staged {role} source lacks a GSM prefix: {path}")
        gsm = match.group(0).upper()
        if gsm in sources:
            raise IntakeBuildError(f"Multiple staged {role} primary sources for {gsm}")
        sources[gsm] = path
    return sources


def _first_metadata_value(row: Mapping[str, Any], columns: Sequence[str]) -> str:
    for column in columns:
        value = clean_text(row.get(column, ""))
        if value:
            return value
    return ""


def _normalize_tissue(raw_value: str, cohort: Mapping[str, str]) -> str:
    if cohort["cohort_id"] == "GSE294273_GSE294274":
        return "lymph_node_metastasis"
    text = clean_text(raw_value).casefold()
    if not text:
        return cohort["default_tissue"]
    if "blood" in text:
        return "blood"
    if "lymph" in text and "node" in text:
        return "lymph_node"
    if "subcu" in text and "met" in text:
        return "subcutaneous_metastasis"
    if "paracancer" in text or "adjacent" in text:
        return "adjacent_tissue"
    if "renal" in text and "tumor" in text:
        return "renal_tumor"
    if "breast" in text and "tumor" in text:
        return "breast_tumor"
    if "tumor" in text:
        return "tumor"
    return re.sub(r"[^a-z0-9]+", "_", text).strip("_")


def _derive_donor(
    row: Mapping[str, Any], sample_key: str, source_path: Path
) -> str:
    explicit = _first_metadata_value(row, ("char:donor", "char:patient id"))
    if explicit:
        return explicit
    search_text = " ".join((sample_key, source_path.name))
    match = re.search(r"(?:BC|GU|HM|LC|K|E)\d+", search_text, flags=re.IGNORECASE)
    if match:
        return match.group(0).upper()
    if re.fullmatch(r"\d+[A-Za-z]", sample_key):
        return f"patient_{sample_key[:-1]}"
    raise IntakeBuildError(
        f"Cannot derive donor for sample {sample_key!r} from shared metadata/source name"
    )


def _sample_metadata(
    row: Mapping[str, Any], cohort: Mapping[str, str], source_path: Path
) -> dict[str, str]:
    sample_id = biological_sample_key(row.get("sample_title", ""))
    raw_tissue = _first_metadata_value(
        row,
        ("char:tissue", "char:tissue type", "char:resident tissue", "char:site", "source"),
    )
    return {
        "sample_id": sample_id,
        "library_id": f"{clean_text(row.get('gsm'))}_{sample_id}",
        "donor_id": _derive_donor(row, sample_id, source_path),
        "tissue": _normalize_tissue(raw_tissue, cohort),
        "diagnosis": cohort["default_diagnosis"],
    }


def derive_paired_sample_specs(
    stage_dir: Path, cohort: Mapping[str, str]
) -> list[dict[str, Any]]:
    """Pair GEX/TCR rows by normalized biological title and GSM-prefixed sources."""
    metasheet = read_source_metasheet(stage_dir)
    gex_rows = metasheet.loc[metasheet["modality"].str.upper().eq("GEX")]
    tcr_rows = metasheet.loc[metasheet["modality"].str.upper().isin({"TCR", "VDJ"})]
    if gex_rows.empty or tcr_rows.empty:
        raise IntakeBuildError(f"{cohort['cohort_id']} requires both GEX and TCR metadata rows")
    tcr_by_key: dict[str, Mapping[str, Any]] = {}
    for _, row in tcr_rows.iterrows():
        key = biological_sample_key(row["sample_title"])
        if key in tcr_by_key:
            raise IntakeBuildError(f"Duplicate TCR metadata pairing key {key!r}")
        tcr_by_key[key] = row
    gex_sources = _primary_sources_by_gsm(stage_dir, "gex")
    tcr_sources = _primary_sources_by_gsm(stage_dir, "vdj")
    specs: list[dict[str, Any]] = []
    for _, gex_row in gex_rows.sort_values(["gse", "gsm"]).iterrows():
        key = biological_sample_key(gex_row["sample_title"])
        tcr_row = tcr_by_key.pop(key, None)
        if tcr_row is None:
            raise IntakeBuildError(f"No paired TCR metadata row for biological sample {key!r}")
        gex_gsm = gex_row["gsm"].upper()
        tcr_gsm = clean_text(tcr_row["gsm"]).upper()
        if gex_gsm not in gex_sources or tcr_gsm not in tcr_sources:
            raise IntakeBuildError(
                f"Missing staged source for paired sample {key}: GEX={gex_gsm}, TCR={tcr_gsm}"
            )
        spec = _sample_metadata(gex_row, cohort, gex_sources[gex_gsm])
        spec.update(
            {
                "gex_path": gex_sources[gex_gsm],
                "tcr_path": tcr_sources[tcr_gsm],
                "gex_gsm": gex_gsm,
                "tcr_gsm": tcr_gsm,
            }
        )
        specs.append(spec)
    if tcr_by_key:
        raise IntakeBuildError(f"Unpaired TCR metadata rows remain: {sorted(tcr_by_key)}")
    ensure_unique_key(pd.DataFrame(specs), ("sample_id", "library_id"), "derived paired samples")
    return specs


def derive_combined_sample_specs(
    stage_dir: Path, cohort: Mapping[str, str]
) -> list[dict[str, Any]]:
    """Map one GSM-prefixed combined source to its exact shared metadata row."""
    metasheet = read_source_metasheet(stage_dir).set_index("gsm", drop=False)
    sources = _primary_sources_by_gsm(stage_dir, "combined")
    specs: list[dict[str, Any]] = []
    for gsm, source in sorted(sources.items()):
        if gsm not in metasheet.index:
            raise IntakeBuildError(f"No shared metadata row for staged combined source {gsm}")
        row = metasheet.loc[gsm]
        if isinstance(row, pd.DataFrame):
            raise IntakeBuildError(f"Duplicate shared metadata rows for {gsm}")
        spec = _sample_metadata(row, cohort, source)
        spec.update({"gex_path": source, "gex_gsm": gsm})
        specs.append(spec)
    if not specs:
        raise IntakeBuildError(f"No staged combined sources for {cohort['cohort_id']}")
    ensure_unique_key(pd.DataFrame(specs), ("sample_id", "library_id"), "derived combined samples")
    return specs


def attach_sample_metadata(
    adata: ad.AnnData,
    sample: Mapping[str, Any],
    cohort: Mapping[str, str],
) -> None:
    metadata = {
        "sample_id": clean_text(sample.get("sample_id")),
        "library_id": clean_text(sample.get("library_id")),
        "donor_id": clean_text(sample.get("donor_id") or sample.get("donor")),
        "tissue": clean_text(sample.get("tissue")) or cohort["default_tissue"],
        "diagnosis": clean_text(sample.get("diagnosis")) or cohort["default_diagnosis"],
    }
    unresolved = [field for field, value in metadata.items() if not value]
    if unresolved:
        raise IntakeBuildError(
            f"Unresolved metadata for cohort {cohort['cohort_id']}: {unresolved}"
        )
    if "barcode" in adata.obs:
        full_barcodes = clean_series(adata.obs["barcode"])
    else:
        full_barcodes = pd.Series(adata.obs_names.astype(str), index=adata.obs_names).map(
            extract_full_barcode
        )
    cores = full_barcodes.map(barcode_core)
    if cores.eq("").any() or cores.duplicated(keep=False).any():
        duplicate_count = int(cores.duplicated(keep=False).sum())
        raise IntakeBuildError(
            f"Sample {metadata['sample_id']} has blank or duplicate barcode_core values "
            f"(duplicate rows={duplicate_count})"
        )
    adata.obs["barcode"] = full_barcodes.to_numpy(dtype=object)
    adata.obs["barcode_core"] = cores.to_numpy(dtype=object)
    for field, value in metadata.items():
        adata.obs[field] = value
    adata.obs["donor"] = metadata["donor_id"]
    adata.obs["technology_simple"] = "10x 5'"
    adata.obs["source_accession"] = cohort["primary_accession"]
    adata.obs_names = pd.Index(
        [f"{metadata['library_id']}:{core}" for core in cores], name="cell_id"
    )


def standalone_sample_qc(adata: ad.AnnData, label: str) -> dict[str, Any]:
    """Validate one sample completely before it is eligible for concatenation."""
    adata.X = validate_sparse_integer_counts(adata.X, f"{label} X")
    if not adata.obs_names.is_unique:
        raise IntakeBuildError(f"{label} has non-unique obs_names")
    missing_columns = [column for column in REQUIRED_CELL_METADATA if column not in adata.obs]
    if missing_columns:
        raise IntakeBuildError(f"{label} lacks cell metadata columns {missing_columns}")
    for column in REQUIRED_CELL_METADATA:
        values = clean_series(adata.obs[column])
        if values.eq("").any():
            raise IntakeBuildError(f"{label} has blank {column} values")
        adata.obs[column] = values
    adata.obs = ensure_canonical_tcr(adata.obs)
    missing_tcr = [column for column in (*TCR_COLUMNS, *TCR_FLAGS) if column not in adata.obs]
    if missing_tcr:
        raise IntakeBuildError(f"{label} lacks canonical TCR fields/flags: {missing_tcr}")
    ensure_unique_key(adata.obs, ("sample_id", "barcode_core"), f"{label} join keys")
    values = adata.X.data
    report = {
        "label": label,
        "passed": True,
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "nnz": int(adata.X.nnz),
        "count_min": float(values.min()) if values.size else 0.0,
        "count_max": float(values.max()) if values.size else 0.0,
        "tcr_positive_cells": int(
            (adata.obs["has_any_ab_tcr"] | adata.obs["has_any_gd_tcr"]).sum()
        ),
        "join_key": "sample_id+barcode_core",
        "obs_names_unique": True,
        "sparse_integer_x": True,
    }
    adata.uns["standalone_qc_passed"] = True
    adata.uns["standalone_qc"] = json.dumps(report, sort_keys=True)
    return report


def concatenate_samples(samples: Sequence[ad.AnnData], cohort_id: str) -> ad.AnnData:
    if not samples:
        raise IntakeBuildError(f"No sample matrices were built for {cohort_id}")
    failed = [
        index
        for index, sample in enumerate(samples)
        if not bool(sample.uns.get("standalone_qc_passed", False))
    ]
    if failed:
        raise IntakeBuildError(
            f"{cohort_id} sample(s) lack passing standalone QC before merge: {failed}"
        )
    result = ad.concat(
        list(samples),
        axis=0,
        join="outer",
        merge="first",
        uns_merge="first",
        index_unique=None,
        fill_value=0,
    )
    result.X = validate_sparse_integer_counts(result.X, f"{cohort_id} concatenated X")
    if not result.obs_names.is_unique:
        raise IntakeBuildError(f"Concatenated {cohort_id} cell IDs are not unique")
    ensure_unique_key(result.obs, ("sample_id", "barcode_core"), f"{cohort_id} expression cells")
    return result


def build_tenx_gex_vdj(stage_dir: Path, cohort: Mapping[str, str]) -> ad.AnnData:
    specs = derive_paired_sample_specs(stage_dir, cohort)
    samples: list[ad.AnnData] = []
    join_summaries: list[dict[str, Any]] = []
    qc_reports: list[dict[str, Any]] = []
    for spec in specs:
        sample = read_expression_source(Path(spec["gex_path"]))
        attach_sample_metadata(sample, spec, cohort)
        tcr_path = Path(spec["tcr_path"])
        tcr = standardize_tcr_table(tcr_path, spec["sample_id"])
        tcr = tcr.loc[tcr["sample_id"].eq(spec["sample_id"])].copy()
        summary = join_tcr_by_sample_barcode(sample, tcr)
        summary.update({"sample_id": spec["sample_id"], "tcr_source": str(tcr_path)})
        join_summaries.append(summary)
        sample.obs["tcr_schema_provenance"] = "productive_contigs"
        qc_reports.append(standalone_sample_qc(sample, f"{cohort['cohort_id']}:{spec['sample_id']}"))
        samples.append(sample)
    result = concatenate_samples(samples, cohort["cohort_id"])
    result.uns["tcr_join_summaries"] = json.dumps(join_summaries, sort_keys=True)
    result.uns["standalone_sample_qc"] = json.dumps(qc_reports, sort_keys=True)
    return result


def parse_embedded_tra_trb(values: pd.Series) -> pd.DataFrame:
    """Parse GSE159251's partial TRA_CDR3|TRB_CDR3 representation."""
    result = empty_tcr_frame(values.index)
    for index, raw_value in values.items():
        text = clean_text(raw_value)
        if text.casefold() in {"notcr", "no_tcr", "no tcr", "unassigned"}:
            text = ""
        if not text:
            continue
        if text.count("|") != 1:
            raise IntakeBuildError(
                f"GSE159251 embedded TCR must contain exactly one '|': index={index!r}, value={text!r}"
            )
        tra, trb = (clean_text(part) for part in text.split("|", 1))
        result.at[index, "TRA_cdr3"] = tra
        result.at[index, "TRB_cdr3"] = trb
    return result


def _read_exported_gse159251(directory: Path) -> tuple[ad.AnnData, pd.DataFrame]:
    sample = read_mtx_counts(directory)
    metadata_path = _pick_one(
        list(directory.glob("metadata.csv")) + list(directory.glob("metadata.csv.gz")),
        "GSE159251 exported metadata",
    )
    metadata = pd.read_csv(metadata_path, compression="infer", dtype=str, keep_default_na=False)
    barcode_values = _first_column(metadata, ("barcode", "cell_barcode"), required=True).map(
        extract_full_barcode
    )
    metadata = metadata.copy()
    metadata["_barcode"] = barcode_values
    ensure_unique_key(metadata, ("_barcode",), "GSE159251 exported metadata")
    metadata = metadata.set_index("_barcode").reindex(sample.obs_names.map(extract_full_barcode))
    if metadata.index.hasnans or metadata.isna().all(axis=1).any():
        raise IntakeBuildError("GSE159251 exported metadata does not cover every count-matrix barcode")
    metadata.index = sample.obs_names
    return sample, metadata


def _export_gse159251_rds(
    source: Path,
    destination: Path,
    r_exporter: Path,
    tcr_column: str,
) -> None:
    command = [
        "Rscript",
        str(r_exporter),
        "--input",
        str(source),
        "--output-dir",
        str(destination),
    ]
    if tcr_column:
        command.extend(["--tcr-column", tcr_column])
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    if completed.returncode:
        raise IntakeBuildError(
            f"GSE159251 RDS export failed for {source}:\n{completed.stdout}\n{completed.stderr}"
        )


def build_gse159251(
    stage_dir: Path,
    cohort: Mapping[str, str],
    work_root: Path,
    r_exporter: Path,
) -> ad.AnnData:
    specs = derive_combined_sample_specs(stage_dir, cohort)
    samples: list[ad.AnnData] = []
    qc_reports: list[dict[str, Any]] = []
    with tempfile.TemporaryDirectory(prefix="gse159251_export_", dir=work_root) as temporary:
        temporary_root = Path(temporary)
        for row_number, spec in enumerate(specs, start=1):
            source = Path(spec["gex_path"])
            if source.is_dir():
                exported = source
            else:
                exported = temporary_root / f"sample_{row_number:03d}"
                exported.mkdir()
                _export_gse159251_rds(
                    source,
                    exported,
                    r_exporter,
                    "",
                )
            sample, metadata = _read_exported_gse159251(exported)
            attach_sample_metadata(sample, spec, cohort)
            tcr_values = _first_column(
                metadata,
                (
                    "embedded_tcr",
                    "TCR",
                    "TRA_CDR3|TRB_CDR3",
                ),
                required=True,
            )
            partial = parse_embedded_tra_trb(tcr_values)
            for column in TCR_COLUMNS:
                sample.obs[column] = partial[column].to_numpy()
            sample.obs = ensure_canonical_tcr(sample.obs)
            sample.obs["tcr_schema_provenance"] = "GSE159251:partial_embedded_TRA_CDR3|TRB_CDR3"
            qc_reports.append(
                standalone_sample_qc(sample, f"{cohort['cohort_id']}:{spec['sample_id']}")
            )
            samples.append(sample)
    result = concatenate_samples(samples, cohort["cohort_id"])
    result.uns["tcr_schema"] = "partial_embedded_tra_trb_cdr3"
    result.uns["standalone_sample_qc"] = json.dumps(qc_reports, sort_keys=True)
    result.uns["tcr_unavailable_fields"] = ";".join(
        f"{chain}_{suffix}"
        for chain in ("TRA", "TRB")
        for suffix in CHAIN_SUFFIXES
        if suffix != "cdr3"
    )
    return result


def map_embedded_airr_tra_trb(obs: pd.DataFrame) -> pd.DataFrame:
    """Map productive AIRR chain slots in GSE315928 to canonical TRA/TRB fields."""
    if any(column in obs for column in ("TRA_cdr3", "TRB_cdr3")):
        canonical = ensure_canonical_tcr(obs)[list(TCR_COLUMNS)].copy()
        if canonical["TRA_cdr3"].ne("").any() or canonical["TRB_cdr3"].ne("").any():
            return canonical

    result = empty_tcr_frame(obs.index)
    prefixes = ("IR_VJ_1", "IR_VJ_2", "IR_VDJ_1", "IR_VDJ_2")
    required_slot_columns = [f"{prefix}_locus" for prefix in prefixes]
    if not any(column in obs for column in required_slot_columns):
        raise IntakeBuildError(
            "Embedded AIRR source has neither canonical TRA/TRB fields nor IR_VJ/IR_VDJ slots"
        )

    mapped_cells = 0
    for index, row in obs.iterrows():
        by_chain: dict[str, list[dict[str, Any]]] = {"TRA": [], "TRB": []}
        for prefix in prefixes:
            locus = clean_text(row.get(f"{prefix}_locus", "")).upper()
            if locus not in by_chain:
                continue
            if not normalize_bool(row.get(f"{prefix}_productive", False), field="AIRR productive"):
                continue
            cdr3 = clean_text(row.get(f"{prefix}_junction_aa", ""))
            if not cdr3:
                continue
            duplicate_count = pd.to_numeric(
                pd.Series([row.get(f"{prefix}_duplicate_count", 0)]), errors="coerce"
            ).fillna(0).astype(np.int64).iloc[0]
            consensus_count = pd.to_numeric(
                pd.Series([row.get(f"{prefix}_consensus_count", 0)]), errors="coerce"
            ).fillna(0).astype(np.int64).iloc[0]
            by_chain[locus].append(
                {
                    "cdr3": cdr3,
                    "v": clean_text(row.get(f"{prefix}_v_call", "")),
                    "d": clean_text(row.get(f"{prefix}_d_call", "")),
                    "j": clean_text(row.get(f"{prefix}_j_call", "")),
                    "cdr3_nt": clean_text(row.get(f"{prefix}_junction", "")),
                    "clone_id": clean_text(row.get("clonotype", "")),
                    "umis": int(duplicate_count),
                    "reads": int(consensus_count),
                }
            )
        cell_has_chain = False
        for chain, candidates in by_chain.items():
            if not candidates:
                continue
            candidates.sort(
                key=lambda item: (-item["umis"], -item["reads"], item["cdr3"])
            )
            selected = candidates[0]
            for suffix in CHAIN_SUFFIXES:
                result.at[index, f"{chain}_{suffix}"] = selected[suffix]
            cell_has_chain = True
        mapped_cells += int(cell_has_chain)
    if mapped_cells == 0:
        raise IntakeBuildError("No productive embedded AIRR TRA/TRB chains were mapped")
    return result


def read_gse315928_publisher_metadata(stage_dir: Path) -> pd.DataFrame:
    sources = [
        path
        for path in staged_sources(stage_dir, "sample_metadata")
        if path.suffix.casefold() in {".xlsx", ".xls"}
    ]
    path = _pick_one(sources, "GSE315928 publisher sample XLSX")
    frame = pd.read_excel(path, sheet_name="Metadata", header=6, dtype=str).fillna("")
    required = {"*title", "**tissue", "Sample ID"}
    if not required.issubset(frame.columns):
        raise IntakeBuildError(f"GSE315928 publisher XLSX lacks columns {sorted(required)}")
    processed_columns = [
        column for column in frame.columns if str(column).startswith("processed data file")
    ]
    if not processed_columns:
        raise IntakeBuildError("GSE315928 publisher XLSX lacks processed-data filenames")
    processed = pd.Series([""] * len(frame), index=frame.index, dtype=object)
    for column in processed_columns:
        values = clean_series(frame[column])
        processed = processed.mask(processed.eq("") & values.ne(""), values)
    result = pd.DataFrame(
        {
            "sample_title": clean_series(frame["*title"]),
            "donor_id": clean_series(frame["Sample ID"]),
            "tissue": clean_series(frame["**tissue"]),
            "processed_file": processed,
        }
    )
    result = result.loc[result["sample_title"].ne("") & result["processed_file"].ne("")].copy()
    ensure_unique_key(result, ("sample_title",), "GSE315928 publisher metadata")
    ensure_unique_key(result, ("processed_file",), "GSE315928 publisher processed files")
    return result


def build_gse315928(stage_dir: Path, cohort: Mapping[str, str]) -> ad.AnnData:
    shared = read_source_metasheet(stage_dir).set_index("gsm", drop=False)
    publisher = read_gse315928_publisher_metadata(stage_dir).set_index("sample_title", drop=False)
    sources = _primary_sources_by_gsm(stage_dir, "combined")
    samples: list[ad.AnnData] = []
    qc_reports: list[dict[str, Any]] = []
    for gsm, source_path in sorted(sources.items()):
        if gsm not in shared.index:
            raise IntakeBuildError(f"GSE315928 shared metadata lacks {gsm}")
        shared_row = shared.loc[gsm]
        sample_title = clean_text(shared_row["sample_title"])
        if sample_title not in publisher.index:
            raise IntakeBuildError(f"GSE315928 publisher XLSX lacks sample {sample_title}")
        publisher_row = publisher.loc[sample_title]
        expected_file = clean_text(publisher_row["processed_file"])
        if expected_file not in source_path.name:
            raise IntakeBuildError(
                f"GSE315928 source mismatch for {gsm}: expected {expected_file}, found {source_path.name}"
            )
        spec = _sample_metadata(shared_row, cohort, source_path)
        spec["donor_id"] = clean_text(publisher_row["donor_id"])
        spec["tissue"] = _normalize_tissue(clean_text(publisher_row["tissue"]), cohort)
        sample = read_expression_source(source_path)
        mapped_tcr = map_embedded_airr_tra_trb(sample.obs)
        attach_sample_metadata(sample, spec, cohort)
        for column in TCR_COLUMNS:
            sample.obs[column] = mapped_tcr[column].to_numpy()
        sample.obs = ensure_canonical_tcr(sample.obs)
        sample.obs["tcr_schema_provenance"] = "GSE315928:embedded_productive_AIRR_TRA_TRB"
        qc_reports.append(
            standalone_sample_qc(sample, f"{cohort['cohort_id']}:{spec['sample_id']}")
        )
        samples.append(sample)
    if len(samples) != len(shared):
        raise IntakeBuildError(
            f"GSE315928 source/shared metadata count mismatch: {len(samples)} vs {len(shared)}"
        )
    result = concatenate_samples(samples, cohort["cohort_id"])
    result.uns["standalone_sample_qc"] = json.dumps(qc_reports, sort_keys=True)
    result.uns["tcr_schema"] = "embedded_airr_tra_trb"
    return result


def finalize_extension_adata(
    adata: ad.AnnData,
    cohort: Mapping[str, str],
    stage_manifest: Mapping[str, Any],
) -> dict[str, Any]:
    """Apply and validate the complete extension schema before serialization."""
    adata.X = validate_sparse_integer_counts(adata.X, f"{cohort['cohort_id']} final X")
    adata.obs = ensure_canonical_tcr(adata.obs)
    if not adata.obs_names.is_unique:
        raise IntakeBuildError(f"{cohort['cohort_id']} final obs_names are not unique")
    for column in REQUIRED_CELL_METADATA:
        if column not in adata.obs:
            raise IntakeBuildError(f"Final H5AD lacks required metadata column {column}")
        adata.obs[column] = clean_series(adata.obs[column])
        if adata.obs[column].eq("").any():
            raise IntakeBuildError(f"Final H5AD has blank {column} values")
    ensure_unique_key(
        adata.obs,
        ("sample_id", "barcode_core"),
        f"{cohort['cohort_id']} final expression cells",
    )
    missing_tcr = [column for column in (*TCR_COLUMNS, *TCR_FLAGS) if column not in adata.obs]
    if missing_tcr:
        raise IntakeBuildError(f"Final H5AD lacks canonical TCR fields/flags: {missing_tcr}")
    qc_payload = json.loads(clean_text(adata.uns.get("standalone_sample_qc", "[]")) or "[]")
    if not qc_payload or any(not item.get("passed") for item in qc_payload):
        raise IntakeBuildError("Final H5AD lacks complete passing standalone sample QC")

    if cohort["cohort_id"] == "GSE294273_GSE294274":
        if not adata.obs["tissue"].eq("lymph_node_metastasis").all():
            raise IntakeBuildError("GSE294273/GSE294274 tissue must be lymph_node_metastasis")
        if not adata.obs["diagnosis"].eq("melanoma").all():
            raise IntakeBuildError("GSE294273/GSE294274 diagnosis must be melanoma")

    adata.obs["intake_role"] = cohort["intake_role"]
    adata.obs["holdout_type"] = cohort["holdout_type"]
    adata.obs["sealed_holdout"] = parse_bool(cohort["sealed"])
    adata.uns["extension_intake_schema_version"] = 1
    adata.uns["cohort_id"] = cohort["cohort_id"]
    adata.uns["primary_accession"] = cohort["primary_accession"]
    adata.uns["gex_accessions"] = split_tokens(cohort["gex_accessions"])
    adata.uns["vdj_accessions"] = split_tokens(cohort["vdj_accessions"])
    adata.uns["count_matrix_location"] = "X"
    adata.uns["count_matrix_state"] = "raw_sparse_integer_counts"
    adata.uns["count_provenance"] = cohort["count_provenance"]
    adata.uns["tcr_join_rule"] = "sample_id+barcode_core"
    adata.uns["tcr_schema"] = cohort["tcr_schema"]
    adata.uns["canonical_tcr_flags"] = list(TCR_FLAGS)
    adata.uns["intake_role"] = cohort["intake_role"]
    adata.uns["holdout_type"] = cohort["holdout_type"]
    adata.uns["sealed_holdout"] = parse_bool(cohort["sealed"])
    adata.uns["stage_manifest"] = json.dumps(stage_manifest, sort_keys=True)
    adata.var_names_make_unique()
    return {
        "cohort_id": cohort["cohort_id"],
        "passed": True,
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "samples": int(adata.obs["sample_id"].nunique()),
        "libraries": int(adata.obs["library_id"].nunique()),
        "donors": int(adata.obs["donor_id"].nunique()),
        "nnz": int(adata.X.nnz),
        "tcr_positive_cells": int(
            (adata.obs["has_any_ab_tcr"] | adata.obs["has_any_gd_tcr"]).sum()
        ),
        "sealed_holdout": parse_bool(cohort["sealed"]),
        "join_key": "sample_id+barcode_core",
        "standalone_sample_qc": qc_payload,
    }


def build_one_cohort(
    stage_dir: Path,
    cohort: Mapping[str, str],
    work_root: Path,
    r_exporter: Path,
) -> tuple[ad.AnnData, dict[str, Any]]:
    if not parse_bool(cohort["build_enabled"]):
        raise IntakeBuildError(
            f"Cohort {cohort['cohort_id']} is build-disabled: {cohort['exclusion_reason']}"
        )
    stage_manifest = load_stage_manifest(stage_dir)
    adapter = cohort["builder_adapter"]
    if adapter == "tenx_gex_vdj":
        adata = build_tenx_gex_vdj(stage_dir, cohort)
    elif adapter == "gse159251_double_gzip_rds":
        adata = build_gse159251(stage_dir, cohort, work_root, r_exporter)
    elif adapter == "gse315928_embedded_airr":
        adata = build_gse315928(stage_dir, cohort)
    elif adapter == "provenance_blocked":
        raise IntakeBuildError(
            "GSE169246 is provenance-blocked; bootstrap H5AD and mmc3/mmc5 are forbidden"
        )
    else:
        raise IntakeBuildError(f"Unsupported extension builder adapter {adapter!r}")
    summary = finalize_extension_adata(adata, cohort, stage_manifest)
    return adata, summary


def build_extension_cohorts(
    stage_root: Path,
    output_root: Path,
    cohorts: Sequence[dict[str, str]],
    dry_run: bool = False,
    r_exporter: Path = DEFAULT_R_EXPORTER,
    resume_existing: bool = False,
) -> dict[str, Any]:
    """Plan or build all selected, enabled extension cohorts."""
    stage_root = stage_root.expanduser().resolve()
    output_root = assert_extension_output_root(output_root)
    blocked = [
        row["cohort_id"]
        for row in cohorts
        if row["intake_role"] in {"blocked_provenance", "excluded_duplicate"}
    ]
    if blocked and len(blocked) == len(cohorts):
        raise IntakeBuildError(f"Selected cohort(s) are build-disabled: {', '.join(sorted(blocked))}")
    enabled = [row for row in cohorts if parse_bool(row["build_enabled"])]
    plans: list[dict[str, Any]] = []
    for cohort in enabled:
        stage_dir = stage_root / cohort["cohort_id"]
        output_path = output_root / f"{cohort['cohort_id']}.h5ad"
        plans.append(
            {
                "cohort_id": cohort["cohort_id"],
                "adapter": cohort["builder_adapter"],
                "stage_dir": str(stage_dir),
                "stage_ready": (stage_dir / "stage_manifest.json").is_file(),
                "output_path": str(output_path),
                "output_exists": output_path.exists(),
            }
        )
    summary: dict[str, Any] = {
        "dry_run": dry_run,
        "stage_root": str(stage_root),
        "output_root": str(output_root),
        "blocked_or_excluded": blocked,
        "resume_existing": resume_existing,
        "plans": plans,
        "built": [],
    }
    if dry_run:
        return summary
    if not r_exporter.is_file():
        raise FileNotFoundError(f"GSE159251 R exporter is missing: {r_exporter}")
    missing_stages = [plan["cohort_id"] for plan in plans if not plan["stage_ready"]]
    if missing_stages:
        raise FileNotFoundError(f"Missing staged cohort(s): {missing_stages}")
    existing_outputs = [plan["output_path"] for plan in plans if plan["output_exists"]]
    if existing_outputs and not resume_existing:
        raise FileExistsError(f"Refusing to replace extension H5AD(s): {existing_outputs}")
    manifest_path = output_root.parent / DEFAULT_BUILT_MANIFEST.name
    if manifest_path.exists():
        raise FileExistsError(f"Refusing to replace built-H5AD manifest: {manifest_path}")

    output_root.mkdir(parents=True, exist_ok=True)
    built: list[dict[str, Any]] = []
    for cohort in enabled:
        cohort_id = cohort["cohort_id"]
        output_path = output_root / f"{cohort_id}.h5ad"
        qc_path = output_root / f"{cohort_id}.standalone_qc.json"
        if output_path.exists() or qc_path.exists():
            if not resume_existing:
                raise FileExistsError(f"Refusing existing extension output for {cohort_id}")
            built.append(validated_existing_build(output_path, qc_path, cohort))
            continue
        temporary_h5ad = output_root / f".{cohort_id}.{uuid.uuid4().hex}.h5ad"
        temporary_qc = output_root / f".{cohort_id}.{uuid.uuid4().hex}.qc.json"
        try:
            adata, cohort_summary = build_one_cohort(
                stage_root / cohort_id, cohort, output_root, r_exporter
            )
            temporary_qc.write_text(
                json.dumps(cohort_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
            adata.write_h5ad(temporary_h5ad, compression="gzip")
            temporary_h5ad.replace(output_path)
            temporary_qc.replace(qc_path)
            built.append(
                {
                    "cohort_id": cohort_id,
                    "h5ad": str(output_path),
                    "standalone_qc": str(qc_path),
                    "cells": cohort_summary["cells"],
                    "genes": cohort_summary["genes"],
                }
            )
        finally:
            temporary_h5ad.unlink(missing_ok=True)
            temporary_qc.unlink(missing_ok=True)
    summary["built"] = built
    summary["built_manifest"] = str(write_built_manifest(output_root, built))
    return summary


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage-root", type=Path, default=DEFAULT_STAGE_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--cohorts-manifest", type=Path, default=DEFAULT_COHORTS)
    parser.add_argument("--libraries-manifest", type=Path, default=DEFAULT_LIBRARIES)
    parser.add_argument("--canonical-registry", type=Path, default=DEFAULT_CANONICAL_REGISTRY)
    parser.add_argument(
        "--shared-metasheet",
        type=Path,
        default=PROJECT_ROOT / "data/compat/downloads/new_candidata_data/sample_metasheet.csv",
    )
    parser.add_argument("--r-exporter", type=Path, default=DEFAULT_R_EXPORTER)
    parser.add_argument("--cohort", action="append", default=[], help="Build one cohort ID; repeatable.")
    parser.add_argument(
        "--resume-existing",
        action="store_true",
        help="Reuse only existing cohort outputs that pass serialized and standalone-QC validation.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Report build readiness without reading or writing H5ADs.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    report, cohorts, _ = validate_extension_manifests(
        args.cohorts_manifest,
        args.libraries_manifest,
        args.canonical_registry,
        args.cohort,
        args.shared_metasheet,
    )
    if not report.passed:
        print(json.dumps(report.to_dict(), indent=2, sort_keys=True), file=sys.stderr)
        return 1
    try:
        summary = build_extension_cohorts(
            args.stage_root,
            args.output_root,
            cohorts,
            args.dry_run,
            args.r_exporter,
            args.resume_existing,
        )
    except (FileNotFoundError, FileExistsError, IntakeBuildError, ManifestError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
