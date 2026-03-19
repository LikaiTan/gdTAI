#!/usr/bin/env python3
"""Repair an H5AD by restoring count-like data from adata.raw onto current var."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


# Config
PRESERVE_UNS_KEYS = {"Channels", "genome", "modality", "uid"}


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Restore integer-like counts from adata.raw onto the current var space, "
            "drop stale derived embeddings, and replace the input H5AD atomically."
        )
    )
    parser.add_argument("input_h5ad", help="Path to the H5AD to repair in place.")
    parser.add_argument(
        "--summary-path",
        default="",
        help="Optional markdown path for a short repair report.",
    )
    return parser.parse_args()


def integer_like(values: np.ndarray) -> bool:
    """Return whether the provided numeric vector is integer-like."""
    if values.size == 0:
        return True
    return bool(np.allclose(values, np.round(values), atol=1e-6))


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


def load_repaired_counts(source: ad.AnnData) -> tuple[sp.csr_matrix, dict]:
    """Extract raw counts for the current var space from adata.raw."""
    if source.raw is None:
        raise ValueError("adata.raw is missing; there is no raw source to repair from")

    raw_var_index = pd.Index(source.raw.var_names.astype(str))
    current_var_names = pd.Index(source.var_names.astype(str))
    raw_positions = raw_var_index.get_indexer(current_var_names)

    if (raw_positions < 0).any():
        missing = current_var_names[raw_positions < 0][:10].tolist()
        raise ValueError(f"Current var contains genes absent from adata.raw: {missing}")

    counts = source.raw.X[:, raw_positions]
    if not sp.issparse(counts):
        counts = sp.csr_matrix(counts)
    elif not sp.isspmatrix_csr(counts):
        counts = counts.tocsr()

    max_value = float(counts.data.max()) if counts.nnz else 0.0
    min_value = float(counts.data.min()) if counts.nnz else 0.0
    is_integer_like = integer_like(counts.data)
    if is_integer_like and min_value >= np.iinfo(np.int32).min and max_value <= np.iinfo(np.int32).max:
        counts.data = counts.data.astype(np.int32, copy=False)

    metadata = {
        "raw_n_vars": int(source.raw.n_vars),
        "current_n_vars": int(source.n_vars),
        "nnz": int(counts.nnz),
        "min_value": min_value,
        "max_value": max_value,
        "integer_like": is_integer_like,
    }
    return counts, metadata


def build_repaired_adata(source: ad.AnnData, counts: sp.csr_matrix, counts_meta: dict) -> ad.AnnData:
    """Construct a repaired AnnData with counts in X and stale derived data removed."""
    repaired = ad.AnnData(
        X=counts,
        obs=sanitize_dataframe(source.obs),
        var=sanitize_dataframe(source.var),
    )

    kept_uns = {}
    removed_uns = []
    for key in source.uns.keys():
        if key in PRESERVE_UNS_KEYS:
            kept_uns[key] = source.uns[key]
        else:
            removed_uns.append(key)
    repaired.uns = kept_uns
    repaired.uns["phase0_raw_rescue"] = {
        "source": "adata.raw",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "kept_var_space_n": int(source.n_vars),
        "original_raw_n_vars": counts_meta["raw_n_vars"],
        "counts_nnz": counts_meta["nnz"],
        "counts_dtype": str(counts.dtype),
        "counts_integer_like": bool(counts_meta["integer_like"]),
        "removed_obsm_keys": list(source.obsm.keys()),
        "removed_varm_keys": list(source.varm.keys()),
        "removed_uns_keys": removed_uns,
    }
    return repaired


def validate_repair(repaired: ad.AnnData, original_n_obs: int, original_var_names: pd.Index) -> None:
    """Fail loudly if the repaired object does not match expectations."""
    if repaired.n_obs != original_n_obs:
        raise ValueError(f"Cell count changed unexpectedly: {repaired.n_obs} != {original_n_obs}")
    if not repaired.var_names.equals(original_var_names):
        raise ValueError("Var names changed unexpectedly during raw rescue")
    if repaired.raw is not None:
        raise ValueError("Repaired object should not retain the old raw slot")
    if len(list(repaired.obsm.keys())) != 0:
        raise ValueError(f"Repaired object should not retain obsm entries: {list(repaired.obsm.keys())}")
    if len(list(repaired.varm.keys())) != 0:
        raise ValueError(f"Repaired object should not retain varm entries: {list(repaired.varm.keys())}")

    data = repaired.X.data if sp.issparse(repaired.X) else np.asarray(repaired.X).ravel()
    if not integer_like(np.asarray(data)):
        raise ValueError("Repaired X is not integer-like")


def validate_repair_backed(repaired: ad.AnnData, original_n_obs: int, original_var_names: pd.Index) -> None:
    """Validate the written H5AD using backed reads only."""
    if repaired.n_obs != original_n_obs:
        raise ValueError(f"Cell count changed unexpectedly after write: {repaired.n_obs} != {original_n_obs}")
    if not pd.Index(repaired.var_names.astype(str)).equals(original_var_names):
        raise ValueError("Var names changed unexpectedly after write")
    if repaired.raw is not None:
        raise ValueError("Written H5AD unexpectedly retained a raw slot")
    if len(list(repaired.obsm.keys())) != 0:
        raise ValueError(f"Written H5AD unexpectedly retained obsm entries: {list(repaired.obsm.keys())}")
    if len(list(repaired.varm.keys())) != 0:
        raise ValueError(f"Written H5AD unexpectedly retained varm entries: {list(repaired.varm.keys())}")

    sample = repaired.X[:1000, :100]
    if sp.issparse(sample):
        data = sample.data
    else:
        data = np.asarray(sample).ravel()
    if not integer_like(np.asarray(data)):
        raise ValueError("Written H5AD X sample is not integer-like")


def write_summary(summary_path: Path, input_path: Path, repaired: ad.AnnData) -> None:
    """Write a short markdown summary of the repair."""
    rescue = repaired.uns["phase0_raw_rescue"]
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("# H5AD Raw Rescue Summary\n\n")
        handle.write(f"- Input path: `{input_path}`\n")
        handle.write("- Repair source: `adata.raw`\n")
        handle.write(f"- Repaired cells: {repaired.n_obs}\n")
        handle.write(f"- Repaired genes: {repaired.n_vars}\n")
        handle.write(f"- Counts nnz: {rescue['counts_nnz']}\n")
        handle.write(f"- Counts dtype: `{rescue['counts_dtype']}`\n")
        handle.write(f"- Counts integer-like: {rescue['counts_integer_like']}\n")
        handle.write(f"- Removed obsm keys: `{', '.join(rescue['removed_obsm_keys'])}`\n")
        handle.write(f"- Removed varm keys: `{', '.join(rescue['removed_varm_keys'])}`\n")
        handle.write(f"- Removed uns keys: `{', '.join(rescue['removed_uns_keys'])}`\n")


def main() -> None:
    """Repair the requested H5AD in place via an atomic temp-file replacement."""
    args = parse_args()
    input_path = Path(args.input_h5ad).resolve()
    temp_path = input_path.with_suffix(input_path.suffix + ".tmp")
    summary_path = Path(args.summary_path).resolve() if args.summary_path else None

    source = ad.read_h5ad(input_path, backed="r")
    original_n_obs = int(source.n_obs)
    original_var_names = pd.Index(source.var_names.astype(str))

    try:
        counts, counts_meta = load_repaired_counts(source)
        repaired = build_repaired_adata(source, counts, counts_meta)
    finally:
        if getattr(source, "file", None) is not None:
            source.file.close()

    validate_repair(repaired, original_n_obs, original_var_names)

    if temp_path.exists():
        temp_path.unlink()
    repaired.write_h5ad(temp_path)
    temp_loaded = ad.read_h5ad(temp_path, backed="r")
    try:
        validate_repair_backed(temp_loaded, original_n_obs, original_var_names)
    finally:
        if getattr(temp_loaded, "file", None) is not None:
            temp_loaded.file.close()

    temp_path.replace(input_path)

    if summary_path is not None:
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        write_summary(summary_path, input_path, repaired)

    print(f"Repaired {input_path}")
    print(f"n_obs={repaired.n_obs}")
    print(f"n_vars={repaired.n_vars}")
    print(f"nnz={repaired.X.nnz if sp.issparse(repaired.X) else int(np.count_nonzero(repaired.X))}")


if __name__ == "__main__":
    main()
