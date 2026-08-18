#!/usr/bin/env python3
"""Create and validate a metadata-only full-atlas H5AD candidate.

This is an explicit, checksum-bound transaction. It never changes the source
atlas or the canonical symlink and it never writes TCR chain calls.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import time
from collections import Counter
from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pyarrow.compute as pc
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[2]
SOURCE = Path("/ssd/tnk_phase3/Integrated_dataset/full_atlas/integrated_full_atlas.h5ad")
OVERLAY = ROOT / "Integrated_dataset/tables/metadata_harmonization/full_atlas_v2/metadata_overlay_v2.parquet"
CANDIDATE = Path(
    "/ssd/tnk_phase3/Integrated_dataset/full_atlas/metadata_corrected/integrated_full_atlas.h5ad"
)
LOG_DIR = ROOT / "Integrated_dataset/logs/metadata_harmonization/full_atlas_v2"

SOURCE_SHA256 = "f5fc491a70f12adeeda5764cb116a59bb441460285905f297bbc2ff691559802"
OVERLAY_SHA256 = "4da4eea32e9d275de790775e2a0f59d4f6553d72756e9c8c6935f35bb398984f"
EXPECTED_CELLS = 5_933_312
EXPECTED_SOURCES = 40
MIN_FREE_AFTER_COPY = 150 * 1024**3

KEY_COLUMNS = ("source_gse_id", "original_cell_id")
ADD_COLUMNS = (
    "source_accession_harmonized_v2",
    "tissue_harmonized_v2",
    "specimen_context_harmonized_v2",
    "tumor_type_harmonized_v2",
    "sample_id_harmonized_v2",
    "library_id_harmonized_v2",
    "donor_id_harmonized_v2",
    "sample_identity_rule_v2",
    "tcr_library_join_id_v2",
    "metadata_rule_id_v2",
    "metadata_evidence_url_v2",
    "metadata_evidence_level_v2",
    "metadata_status_v2",
)
FORBIDDEN_CONTEXT = {"cd4", "cd8", "nk", "treg", "naive", "proliferating"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--overlay", type=Path, default=OVERLAY)
    parser.add_argument("--candidate", type=Path, default=CANDIDATE)
    parser.add_argument("--log-dir", type=Path, default=LOG_DIR)
    return parser.parse_args()


def sha256_file(path: Path, chunk_size: int = 64 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def file_state(path: Path, include_sha: bool = False) -> dict[str, Any]:
    stat = path.stat()
    state: dict[str, Any] = {
        "path": str(path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "inode": stat.st_ino,
    }
    if include_sha:
        state["sha256"] = sha256_file(path)
    return state


def decode(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [value.decode() if isinstance(value, bytes) else str(value) for value in values],
        dtype=object,
    )


def obs_values(obs: h5py.Group, name: str, positions: np.ndarray | slice = slice(None)) -> np.ndarray:
    node = obs[name]
    if isinstance(node, h5py.Group) and {"categories", "codes"}.issubset(node):
        categories = decode(node["categories"][:])
        codes = node["codes"][positions]
        result = np.full(len(codes), "", dtype=object)
        present = codes >= 0
        result[present] = categories[codes[present]]
        return result
    if isinstance(node, h5py.Dataset):
        return decode(node[positions])
    raise ValueError(f"Unsupported obs encoding: {name}")


def categorical_parts(obs: h5py.Group, name: str) -> tuple[np.ndarray, h5py.Dataset]:
    node = obs[name]
    if not isinstance(node, h5py.Group) or not {"categories", "codes"}.issubset(node):
        raise ValueError(f"Expected categorical obs column: {name}")
    return decode(node["categories"][:]), node["codes"]


def canonical_categorical_hash(categories: list[str] | np.ndarray, codes: np.ndarray) -> str:
    digest = hashlib.sha256()
    for value in categories:
        encoded = str(value).encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "little"))
        digest.update(encoded)
    digest.update(np.asarray(codes, dtype="<i4").tobytes())
    return digest.hexdigest()


def sampled_dataset_signature(dataset: h5py.Dataset) -> dict[str, Any]:
    signature: dict[str, Any] = {
        "shape": list(dataset.shape),
        "dtype": str(dataset.dtype),
        "chunks": list(dataset.chunks) if dataset.chunks else None,
        "compression": dataset.compression,
    }
    if dataset.size == 0:
        signature["sample_sha256"] = hashlib.sha256(b"").hexdigest()
        return signature
    first_axis = dataset.shape[0] if dataset.shape else 1
    points = sorted({0, first_axis // 4, first_axis // 2, (3 * first_axis) // 4, first_axis - 1})
    digest = hashlib.sha256()
    for point in points:
        value = dataset[point] if dataset.shape else dataset[()]
        array = np.asarray(value)
        if array.dtype.kind in {"O", "S", "U"}:
            for item in array.reshape(-1):
                encoded = item if isinstance(item, bytes) else str(item).encode("utf-8")
                digest.update(len(encoded).to_bytes(8, "little"))
                digest.update(encoded)
        else:
            digest.update(np.ascontiguousarray(array).tobytes())
    signature["sample_sha256"] = digest.hexdigest()
    return signature


def structural_signature(handle: h5py.File, old_obs_columns: list[str]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for root_name in ("X", "obsm", "var"):
        root = handle[root_name]
        datasets: dict[str, Any] = {}
        root.visititems(
            lambda name, node: datasets.__setitem__(name, sampled_dataset_signature(node))
            if isinstance(node, h5py.Dataset)
            else None
        )
        result[root_name] = datasets
    obs_result: dict[str, Any] = {}
    obs = handle["obs"]
    for column in old_obs_columns:
        node = obs[column]
        if isinstance(node, h5py.Dataset):
            obs_result[column] = {"": sampled_dataset_signature(node)}
        else:
            obs_result[column] = {
                name: sampled_dataset_signature(child)
                for name, child in node.items()
                if isinstance(child, h5py.Dataset)
            }
    result["obs_existing"] = obs_result
    return result


def code_dtype(n_categories: int) -> np.dtype:
    if n_categories <= np.iinfo(np.int8).max:
        return np.dtype("int8")
    if n_categories <= np.iinfo(np.int16).max:
        return np.dtype("int16")
    return np.dtype("int32")


def write_categorical(obs: h5py.Group, name: str, categories: list[str], codes: np.ndarray) -> None:
    group = obs.create_group(name)
    group.attrs["encoding-type"] = "categorical"
    group.attrs["encoding-version"] = "0.2.0"
    group.attrs["ordered"] = False
    strings = np.asarray(categories, dtype=object)
    category_ds = group.create_dataset(
        "categories",
        data=strings,
        dtype=h5py.string_dtype(encoding="utf-8"),
        compression="gzip",
    )
    category_ds.attrs["encoding-type"] = "string-array"
    category_ds.attrs["encoding-version"] = "0.2.0"
    code_ds = group.create_dataset(
        "codes",
        data=codes.astype(code_dtype(len(categories)), copy=False),
        compression="gzip",
        chunks=(min(len(codes), 262_144),),
    )
    code_ds.attrs["encoding-type"] = "array"
    code_ds.attrs["encoding-version"] = "0.2.0"


def build_reordered_codes(source: Path, overlay: Path) -> tuple[dict[str, list[str]], dict[str, np.ndarray]]:
    table = pq.read_table(overlay, columns=[*KEY_COLUMNS, *ADD_COLUMNS])
    if table.num_rows != EXPECTED_CELLS:
        raise RuntimeError(f"Unexpected overlay row count: {table.num_rows:,}")
    categories: dict[str, list[str]] = {}
    value_to_code: dict[str, dict[str, int]] = {}
    codes: dict[str, np.ndarray] = {}
    for name in ADD_COLUMNS:
        values = ["" if value is None else str(value) for value in pc.unique(table[name]).to_pylist()]
        categories[name] = sorted(set(values))
        value_to_code[name] = {value: code for code, value in enumerate(categories[name])}
        codes[name] = np.full(EXPECTED_CELLS, -1, dtype=np.int32)

    with h5py.File(source, "r") as handle:
        obs = handle["obs"]
        source_categories, source_codes = categorical_parts(obs, "source_gse_id")
        original_categories, original_codes = categorical_parts(obs, "original_cell_id")
        if len(source_codes) != EXPECTED_CELLS:
            raise RuntimeError(f"Unexpected source atlas row count: {len(source_codes):,}")
        if len(source_categories) != EXPECTED_SOURCES:
            raise RuntimeError(f"Unexpected source count: {len(source_categories)}")

        for source_code, source_name in enumerate(source_categories):
            positions = np.flatnonzero(source_codes[:] == source_code)
            atlas_ids = original_categories[original_codes[positions]]
            source_table = table.filter(pc.equal(table["source_gse_id"], source_name))
            if source_table.num_rows != len(positions):
                raise RuntimeError(
                    f"Row mismatch for {source_name}: atlas={len(positions):,}, overlay={source_table.num_rows:,}"
                )
            overlay_ids = np.asarray(source_table["original_cell_id"].to_pylist(), dtype=object)
            if len(set(overlay_ids)) != len(overlay_ids):
                raise RuntimeError(f"Duplicate overlay keys within {source_name}")
            lookup = {cell_id: index for index, cell_id in enumerate(overlay_ids)}
            order = np.fromiter((lookup.get(cell_id, -1) for cell_id in atlas_ids), dtype=np.int64)
            if np.any(order < 0) or len(set(atlas_ids)) != len(atlas_ids):
                raise RuntimeError(f"Non-bijective source/cell key mapping for {source_name}")
            for name in ADD_COLUMNS:
                source_values = source_table[name].to_pylist()
                mapping = value_to_code[name]
                source_value_codes = np.fromiter(
                    (mapping["" if value is None else str(value)] for value in source_values),
                    dtype=np.int32,
                    count=len(source_values),
                )
                codes[name][positions] = source_value_codes[order]
            print(f"keyed {source_name}: {len(positions):,} cells", flush=True)
    if any(np.any(column_codes < 0) for column_codes in codes.values()):
        raise RuntimeError("At least one additive metadata code was not populated")
    return categories, codes


def validate_semantics(categories: dict[str, list[str]], codes: dict[str, np.ndarray]) -> dict[str, bool]:
    def count(column: str, value: str) -> int:
        try:
            target = categories[column].index(value)
        except ValueError:
            return 0
        return int(np.count_nonzero(codes[column] == target))

    sample_values = set(categories["sample_id_harmonized_v2"])
    library_values = set(categories["library_id_harmonized_v2"])
    tissue_values = set(categories["tissue_harmonized_v2"])
    context_values = set(categories["specimen_context_harmonized_v2"])
    return {
        "sample_nonblank": "" not in sample_values,
        "library_nonblank": "" not in library_values,
        "sample_unresolved_exact": count("sample_id_harmonized_v2", "unresolved") == 4_611,
        "tissue_unresolved_exact": count("tissue_harmonized_v2", "unresolved") == 39_920,
        "context_unresolved_exact": count("specimen_context_harmonized_v2", "unresolved") == 39_920,
        "no_cell_type_tissue_values": not ({value.casefold() for value in tissue_values} & FORBIDDEN_CONTEXT),
        "no_cell_type_context_values": not ({value.casefold() for value in context_values} & FORBIDDEN_CONTEXT),
    }


def main() -> None:
    args = parse_args()
    source = args.source.resolve()
    overlay = args.overlay.resolve()
    candidate = args.candidate.resolve()
    partial = candidate.with_suffix(candidate.suffix + ".partial")
    log_dir = args.log_dir.resolve()
    log_dir.mkdir(parents=True, exist_ok=True)
    candidate.parent.mkdir(parents=True, exist_ok=True)
    if candidate.exists() or partial.exists():
        raise RuntimeError(f"Refusing to overwrite an existing candidate: {candidate} or {partial}")

    started = time.time()
    source_before = file_state(source, include_sha=True)
    overlay_sha = sha256_file(overlay)
    if source_before["sha256"] != SOURCE_SHA256:
        raise RuntimeError("Source atlas SHA-256 does not match the frozen write plan")
    if overlay_sha != OVERLAY_SHA256:
        raise RuntimeError("Metadata overlay SHA-256 does not match the frozen write plan")
    free = shutil.disk_usage(candidate.parent).free
    if free - source_before["size"] < MIN_FREE_AFTER_COPY:
        raise RuntimeError("Insufficient SSD free space for candidate plus 150-GiB reserve")

    with h5py.File(source, "r") as handle:
        old_columns = [str(value.decode() if isinstance(value, bytes) else value) for value in handle["obs"].attrs["column-order"]]
        if set(ADD_COLUMNS) & set(old_columns):
            raise RuntimeError("At least one additive column already exists in the source atlas")
        source_signature = structural_signature(handle, old_columns)

    categories, codes = build_reordered_codes(source, overlay)
    semantic_checks = validate_semantics(categories, codes)
    if not all(semantic_checks.values()):
        raise RuntimeError(f"Pre-write semantic validation failed: {semantic_checks}")
    expected_hashes = {
        name: canonical_categorical_hash(categories[name], codes[name]) for name in ADD_COLUMNS
    }
    expected_counts = {
        name: dict(Counter(categories[name][code] for code in codes[name])) for name in ADD_COLUMNS
    }

    subprocess.run(
        ["cp", "--reflink=auto", "--sparse=always", "--preserve=mode,timestamps", str(source), str(partial)],
        check=True,
    )
    with h5py.File(partial, "r+") as handle:
        obs = handle["obs"]
        for name in ADD_COLUMNS:
            write_categorical(obs, name, categories[name], codes[name])
        string_dtype = h5py.string_dtype(encoding="utf-8")
        del obs.attrs["column-order"]
        obs.attrs.create(
            "column-order",
            np.asarray([*old_columns, *ADD_COLUMNS], dtype=string_dtype),
        )
        handle.flush()

    checks: dict[str, bool] = dict(semantic_checks)
    observed_hashes: dict[str, str] = {}
    observed_counts: dict[str, dict[str, int]] = {}
    with h5py.File(partial, "r") as handle:
        observed_columns = [str(value.decode() if isinstance(value, bytes) else value) for value in handle["obs"].attrs["column-order"]]
        checks["only_additive_obs_columns"] = observed_columns == [*old_columns, *ADD_COLUMNS]
        checks["row_count"] = len(handle["obs"]["source_gse_id"]["codes"]) == EXPECTED_CELLS
        checks["existing_structure_sampled_identity"] = structural_signature(handle, old_columns) == source_signature
        checks["x_shape"] = tuple(handle["X"].attrs["shape"]) == (EXPECTED_CELLS, 27_413)
        checks["x_nnz"] = len(handle["X"]["data"]) == 9_127_088_723
        checks["scvi_shape"] = handle["obsm"]["X_scVI"].shape == (EXPECTED_CELLS, 30)
        checks["umap_shape"] = handle["obsm"]["X_umap"].shape == (EXPECTED_CELLS, 2)
        for name in ADD_COLUMNS:
            observed_categories, observed_code_ds = categorical_parts(handle["obs"], name)
            observed_codes = observed_code_ds[:].astype(np.int32)
            observed_hashes[name] = canonical_categorical_hash(observed_categories, observed_codes)
            observed_counts[name] = dict(Counter(observed_categories[code] for code in observed_codes))
        checks["additive_ordered_hashes"] = observed_hashes == expected_hashes
        checks["additive_value_counts"] = observed_counts == expected_counts

    backed = ad.read_h5ad(partial, backed="r")
    checks["anndata_backed_open"] = backed.n_obs == EXPECTED_CELLS and backed.n_vars == 27_413
    backed.file.close()
    source_after = file_state(source, include_sha=True)
    checks["source_state_and_sha_unchanged"] = source_after == source_before
    checks["no_tcr_chain_calls_applied"] = True
    status = "PASS_METADATA_H5AD_CANDIDATE" if all(checks.values()) else "FAIL_METADATA_H5AD_CANDIDATE"

    manifest = {
        "status": status,
        "source_atlas": source_before,
        "overlay": {"path": str(overlay), "sha256": overlay_sha, "rows": EXPECTED_CELLS},
        "candidate": str(candidate),
        "partial": str(partial),
        "n_cells": EXPECTED_CELLS,
        "n_added_columns": len(ADD_COLUMNS),
        "added_columns": list(ADD_COLUMNS),
        "checks": checks,
        "ordered_value_hashes": observed_hashes,
        "elapsed_seconds": round(time.time() - started, 3),
        "canonical_symlink_changed": False,
        "tcr_chain_calls_applied": False,
    }
    manifest_path = log_dir / "metadata_h5ad_candidate_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if status != "PASS_METADATA_H5AD_CANDIDATE":
        print(json.dumps(manifest, indent=2, sort_keys=True))
        raise SystemExit(1)
    os.replace(partial, candidate)
    manifest["validated_candidate_exists"] = candidate.exists()
    manifest["candidate_size_bytes"] = candidate.stat().st_size
    manifest["candidate_sha256"] = sha256_file(candidate)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (log_dir / "metadata_h5ad_candidate_summary.md").write_text(
        "# Full-atlas metadata candidate\n\n"
        f"- Status: `{status}`\n"
        f"- Candidate: `{candidate}`\n"
        f"- Candidate SHA-256: `{manifest['candidate_sha256']}`\n"
        f"- Cells: {EXPECTED_CELLS:,}\n"
        f"- Additive metadata columns: {len(ADD_COLUMNS)}\n"
        "- Canonical atlas changed: no\n"
        "- TCR chain calls applied: no\n"
        "- Next gate: validated TCR sidecar transaction using this candidate as input\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
