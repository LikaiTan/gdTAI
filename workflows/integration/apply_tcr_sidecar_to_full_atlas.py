#!/usr/bin/env python3
"""Apply the validated TCR sidecar to a separate full-atlas H5AD candidate."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[2]
SOURCE = Path(
    "/ssd/tnk_phase3/Integrated_dataset/full_atlas/metadata_corrected/integrated_full_atlas.h5ad"
)
SIDECAR = ROOT / "Integrated_dataset/tables/tcr_join_rebuild/validated_tcr_replacement_sidecar.parquet"
CANDIDATE = Path(
    "/ssd/tnk_phase3/Integrated_dataset/full_atlas/tcr_corrected/integrated_full_atlas.h5ad"
)
TABLE_DIR = ROOT / "Integrated_dataset/tables/tcr_sidecar_application"
LOG_DIR = ROOT / "Integrated_dataset/logs/tcr_sidecar_application"

SOURCE_SHA256 = "7f1c5e1cac1074a8e2863703bc1862e225defc5ba1a3adbaabd3f6e023d5871c"
SIDECAR_SHA256 = "3114e70719301d693ae1a2bc2c63bac6c8bd57e3e8ac73a88c24320eaabfc2f0"
EXPECTED_CELLS = 5_933_312
EXPECTED_GENES = 27_413
EXPECTED_NNZ = 9_127_088_723
EXPECTED_SIDECAR_ROWS = 3_041_871
EXPECTED_AFFECTED_ROWS = 2_155_409
EXPECTED_AFFECTED_SOURCES = 14
EXPECTED_INELIGIBLE_ATLAS_ROWS = 109
EXPECTED_WHOLE_ANY_AB = 2_270_138
EXPECTED_WHOLE_PAIRED_AB = 1_938_158
MIN_FREE_AFTER_COPY = 150 * 1024**3
EXPECTED_CANONICAL_TARGET = Path(
    "/ssd/tnk_phase3/Integrated_dataset/full_atlas/integrated_full_atlas.h5ad"
)

CHAINS = ("TRA", "TRB", "TRG", "TRD")
CHAIN_STRING_SUFFIXES = ("cdr3", "cdr3_nt", "v", "d", "j", "clone_id")
CHAIN_NUMERIC_SUFFIXES = ("umis", "reads")
CANONICAL_CHAIN_COLUMNS = tuple(
    f"{chain}_{suffix}"
    for chain in CHAINS
    for suffix in (*CHAIN_STRING_SUFFIXES, *CHAIN_NUMERIC_SUFFIXES, "c_gene")
)
CANONICAL_BOOL_MAP = {
    "has_TRA": None,
    "has_TRB": None,
    "has_TRG": None,
    "has_TRD": None,
    "has_TRA_TRB_paired": "has_TRA_TRB_paired_rebuilt",
    "has_TRG_TRD_paired": "has_TRG_TRD_paired_rebuilt",
    "has_any_ab_tcr": "has_any_ab_tcr_rebuilt",
    "has_any_gd_tcr": "has_any_gd_tcr_rebuilt",
}
REPLACED_CATEGORICAL = (*CANONICAL_CHAIN_COLUMNS, "TCRseq", "tcr_schema_provenance")
REPLACED_COLUMNS = (*REPLACED_CATEGORICAL, *CANONICAL_BOOL_MAP)

GENERAL_PROVENANCE = {
    "tcr_join_status_rebuilt_v2": "join_status",
    "tcr_join_reason_rebuilt_v2": "join_reason",
    "tcr_join_donor_id_rebuilt_v2": "join_donor_id",
    "tcr_barcode_core_rebuilt_v2": "barcode_core",
}
GENERAL_BOOL_PROVENANCE = {
    "tcr_sidecar_covered_v2": None,
    "tcr_assignment_eligible_rebuilt_v2": "tcr_assignment_eligible",
    "tcr_replacement_eligible_rebuilt_v2": "replacement_eligible",
}
CHAIN_PROVENANCE = tuple(
    f"{chain}_{suffix}_rebuilt_v2"
    for chain in CHAINS
    for suffix in (
        "umi_available",
        "read_available",
        "n_productive_contigs",
        "source_file",
        "selected_contig_id",
    )
)
ADD_COLUMNS = (*GENERAL_PROVENANCE, *GENERAL_BOOL_PROVENANCE, *CHAIN_PROVENANCE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--sidecar", type=Path, default=SIDECAR)
    parser.add_argument("--candidate", type=Path, default=CANDIDATE)
    parser.add_argument("--table-dir", type=Path, default=TABLE_DIR)
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
    result: dict[str, Any] = {
        "path": str(path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "inode": stat.st_ino,
    }
    if include_sha:
        result["sha256"] = sha256_file(path)
    return result


def decode(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [value.decode() if isinstance(value, bytes) else str(value) for value in values],
        dtype=object,
    )


def categorical_values(obs: h5py.Group, name: str) -> np.ndarray:
    node = obs[name]
    if not isinstance(node, h5py.Group) or not {"categories", "codes"}.issubset(node):
        raise ValueError(f"Expected categorical obs column: {name}")
    categories = decode(node["categories"][:])
    codes = node["codes"][:]
    values = np.full(len(codes), "", dtype=object)
    present = codes >= 0
    values[present] = categories[codes[present]]
    return values


def categorical_values_at(obs: h5py.Group, name: str, positions: np.ndarray) -> np.ndarray:
    node = obs[name]
    categories = decode(node["categories"][:])
    codes = node["codes"][positions]
    values = np.full(len(codes), "", dtype=object)
    present = codes >= 0
    values[present] = categories[codes[present]]
    return values


def code_dtype(n_categories: int) -> np.dtype:
    if n_categories <= np.iinfo(np.int8).max:
        return np.dtype("int8")
    if n_categories <= np.iinfo(np.int16).max:
        return np.dtype("int16")
    return np.dtype("int32")


def write_categorical(obs: h5py.Group, name: str, values: np.ndarray) -> None:
    values = np.asarray(values, dtype=object)
    categories = sorted(set(str(value) for value in values))
    mapping = {value: code for code, value in enumerate(categories)}
    codes = np.fromiter((mapping[str(value)] for value in values), dtype=np.int32, count=len(values))
    if name in obs:
        del obs[name]
    group = obs.create_group(name)
    group.attrs["encoding-type"] = "categorical"
    group.attrs["encoding-version"] = "0.2.0"
    group.attrs["ordered"] = False
    cat_ds = group.create_dataset(
        "categories",
        data=np.asarray(categories, dtype=object),
        dtype=h5py.string_dtype("utf-8"),
        compression="gzip",
    )
    cat_ds.attrs["encoding-type"] = "string-array"
    cat_ds.attrs["encoding-version"] = "0.2.0"
    code_ds = group.create_dataset(
        "codes",
        data=codes.astype(code_dtype(len(categories)), copy=False),
        compression="gzip",
        chunks=(min(len(codes), 262_144),),
    )
    code_ds.attrs["encoding-type"] = "array"
    code_ds.attrs["encoding-version"] = "0.2.0"


def write_boolean(obs: h5py.Group, name: str, values: np.ndarray) -> None:
    if name in obs:
        del obs[name]
    dataset = obs.create_dataset(
        name,
        data=np.asarray(values, dtype=bool),
        compression="gzip",
        chunks=(min(len(values), 262_144),),
    )
    dataset.attrs["encoding-type"] = "array"
    dataset.attrs["encoding-version"] = "0.2.0"


def value_hash(values: np.ndarray) -> str:
    hashes = pd.util.hash_array(np.asarray(values, dtype=object), categorize=True)
    return hashlib.sha256(np.asarray(hashes, dtype="<u8").tobytes()).hexdigest()


def dataset_signature(dataset: h5py.Dataset) -> dict[str, Any]:
    result: dict[str, Any] = {
        "shape": list(dataset.shape),
        "dtype": str(dataset.dtype),
        "chunks": list(dataset.chunks) if dataset.chunks else None,
        "compression": dataset.compression,
    }
    if dataset.size == 0:
        result["sample_sha256"] = hashlib.sha256(b"").hexdigest()
        return result
    n = dataset.shape[0] if dataset.shape else 1
    points = sorted({0, n // 4, n // 2, (3 * n) // 4, n - 1})
    digest = hashlib.sha256()
    for point in points:
        value = np.asarray(dataset[point] if dataset.shape else dataset[()])
        if value.dtype.kind in {"O", "S", "U"}:
            for item in value.reshape(-1):
                encoded = item if isinstance(item, bytes) else str(item).encode()
                digest.update(len(encoded).to_bytes(8, "little"))
                digest.update(encoded)
        else:
            digest.update(np.ascontiguousarray(value).tobytes())
    result["sample_sha256"] = digest.hexdigest()
    return result


def protected_signature(handle: h5py.File, protected_obs: list[str]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for root_name in ("X", "obsm", "var"):
        datasets: dict[str, Any] = {}
        handle[root_name].visititems(
            lambda name, node: datasets.__setitem__(name, dataset_signature(node))
            if isinstance(node, h5py.Dataset)
            else None
        )
        result[root_name] = datasets
    obs_result: dict[str, Any] = {}
    obs = handle["obs"]
    for name in protected_obs:
        node = obs[name]
        if isinstance(node, h5py.Dataset):
            obs_result[name] = {"": dataset_signature(node)}
        else:
            obs_result[name] = {
                child_name: dataset_signature(child)
                for child_name, child in node.items()
                if isinstance(child, h5py.Dataset)
            }
    result["obs"] = obs_result
    return result


def sidecar_string_values(array: Any) -> np.ndarray:
    return np.asarray(["" if value is None else str(value) for value in array.to_pylist()], dtype=object)


def build_keyed_overlays(source: Path, sidecar: Path) -> tuple[dict[str, tuple[np.ndarray, Any]], int]:
    table = pq.read_table(sidecar)
    if table.num_rows != EXPECTED_SIDECAR_ROWS:
        raise RuntimeError(f"Unexpected sidecar row count: {table.num_rows:,}")
    sources = sorted(str(value) for value in pc.unique(table["source_gse_id"]).to_pylist())
    if len(sources) != EXPECTED_AFFECTED_SOURCES:
        raise RuntimeError(f"Unexpected sidecar source count: {len(sources)}")
    overlays: dict[str, tuple[np.ndarray, Any]] = {}
    total = 0
    with h5py.File(source, "r") as handle:
        obs = handle["obs"]
        source_values = categorical_values(obs, "source_gse_id")
        original_values = categorical_values(obs, "original_cell_id")
        for source_name in sources:
            positions = np.flatnonzero(source_values == source_name)
            source_table = table.filter(pc.equal(table["source_gse_id"], source_name))
            side_ids = sidecar_string_values(source_table["source_obs_name"])
            if len(set(side_ids)) != len(side_ids):
                raise RuntimeError(f"Duplicate sidecar source/cell keys: {source_name}")
            lookup = {cell_id: index for index, cell_id in enumerate(side_ids)}
            order = np.fromiter((lookup.get(cell_id, -1) for cell_id in original_values[positions]), dtype=np.int64)
            if np.any(order < 0) or len(positions) != len(set(original_values[positions])):
                raise RuntimeError(f"Incomplete or non-unique atlas key mapping: {source_name}")
            ordered = source_table.take(order)
            overlays[source_name] = (positions, ordered)
            total += len(positions)
            print(f"keyed {source_name}: {len(positions):,} atlas cells", flush=True)
    if total != EXPECTED_AFFECTED_ROWS:
        raise RuntimeError(f"Unexpected affected atlas rows: {total:,}")
    return overlays, total


def replacement_strings(side_table: Any, canonical_name: str) -> np.ndarray:
    if canonical_name == "TCRseq":
        return sidecar_string_values(side_table["TCRseq_rebuilt"])
    if canonical_name == "tcr_schema_provenance":
        return np.full(side_table.num_rows, "validated_tcr_replacement_sidecar_v2", dtype=object)
    chain, suffix = canonical_name.split("_", 1)
    if suffix == "c_gene":
        return np.full(side_table.num_rows, "", dtype=object)
    values = side_table[canonical_name].to_pylist()
    if suffix in CHAIN_NUMERIC_SUFFIXES:
        return np.asarray(["" if value is None else str(int(value)) for value in values], dtype=object)
    return np.asarray(["" if value is None else str(value) for value in values], dtype=object)


def replacement_bool(side_table: Any, canonical_name: str) -> np.ndarray:
    source_name = CANONICAL_BOOL_MAP[canonical_name]
    if source_name is not None:
        return np.asarray(side_table[source_name].to_numpy(zero_copy_only=False), dtype=bool)
    chain = canonical_name.removeprefix("has_")
    return sidecar_string_values(side_table[f"{chain}_cdr3"]) != ""


def provenance_values(side_table: Any, output_name: str) -> tuple[np.ndarray, bool]:
    if output_name in GENERAL_PROVENANCE:
        return sidecar_string_values(side_table[GENERAL_PROVENANCE[output_name]]), False
    if output_name in GENERAL_BOOL_PROVENANCE:
        source_name = GENERAL_BOOL_PROVENANCE[output_name]
        if source_name is None:
            return np.ones(side_table.num_rows, dtype=bool), True
        return np.asarray(side_table[source_name].to_numpy(zero_copy_only=False), dtype=bool), True
    base = output_name.removesuffix("_rebuilt_v2")
    if base.endswith(("_umi_available", "_read_available")):
        return np.asarray(side_table[base].to_numpy(zero_copy_only=False), dtype=bool), True
    values = side_table[base].to_pylist()
    if base.endswith("_n_productive_contigs"):
        return np.asarray(["" if value is None else str(int(value)) for value in values], dtype=object), False
    return np.asarray(["" if value is None else str(value) for value in values], dtype=object), False


def main() -> None:
    args = parse_args()
    source = args.source.resolve()
    sidecar = args.sidecar.resolve()
    candidate = args.candidate.resolve()
    partial = candidate.with_suffix(candidate.suffix + ".partial")
    table_dir = args.table_dir.resolve()
    log_dir = args.log_dir.resolve()
    table_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    candidate.parent.mkdir(parents=True, exist_ok=True)
    if candidate.exists() or partial.exists():
        raise RuntimeError(f"Refusing to overwrite existing candidate or partial: {candidate}")

    started = time.time()
    source_before = file_state(source, include_sha=True)
    sidecar_sha = sha256_file(sidecar)
    if source_before["sha256"] != SOURCE_SHA256:
        raise RuntimeError("Metadata-corrected source SHA-256 does not match the frozen input")
    if sidecar_sha != SIDECAR_SHA256:
        raise RuntimeError("TCR sidecar SHA-256 does not match the validated input")
    if shutil.disk_usage(candidate.parent).free - source_before["size"] < MIN_FREE_AFTER_COPY:
        raise RuntimeError("Insufficient SSD free space for candidate plus 150-GiB reserve")

    overlays, affected_rows = build_keyed_overlays(source, sidecar)
    affected_positions = np.sort(
        np.concatenate([positions for positions, _ in overlays.values()])
    )
    unaffected_mask = np.ones(EXPECTED_CELLS, dtype=bool)
    unaffected_mask[affected_positions] = False

    with h5py.File(source, "r") as handle:
        obs = handle["obs"]
        old_columns = [str(v.decode() if isinstance(v, bytes) else v) for v in obs.attrs["column-order"]]
        missing = set(REPLACED_COLUMNS) - set(old_columns)
        if missing:
            raise RuntimeError(f"Missing canonical TCR columns: {sorted(missing)}")
        if set(ADD_COLUMNS) & set(old_columns):
            raise RuntimeError("At least one rebuilt provenance column already exists")
        protected_obs = [name for name in old_columns if name not in REPLACED_COLUMNS]
        source_signature = protected_signature(handle, protected_obs)
        old_any_ab_total = int(np.asarray(obs["has_any_ab_tcr"][:], dtype=bool).sum())
        old_paired_ab_total = int(np.asarray(obs["has_TRA_TRB_paired"][:], dtype=bool).sum())

    source_rows: list[dict[str, Any]] = []
    with h5py.File(source, "r") as handle:
        obs = handle["obs"]
        for source_name, (positions, side_table) in overlays.items():
            old_any = np.asarray(obs["has_any_ab_tcr"][positions], dtype=bool)
            old_paired = np.asarray(obs["has_TRA_TRB_paired"][positions], dtype=bool)
            new_any = replacement_bool(side_table, "has_any_ab_tcr")
            new_paired = replacement_bool(side_table, "has_TRA_TRB_paired")
            eligible = np.asarray(side_table["tcr_assignment_eligible"].to_numpy(zero_copy_only=False), dtype=bool)
            source_rows.append(
                {
                    "source_gse_id": source_name,
                    "n_atlas_rows": len(positions),
                    "n_old_any_ab": int(old_any.sum()),
                    "n_rebuilt_any_ab": int(new_any.sum()),
                    "n_old_paired_ab": int(old_paired.sum()),
                    "n_rebuilt_paired_ab": int(new_paired.sum()),
                    "n_old_only_any_ab": int((old_any & ~new_any).sum()),
                    "n_rebuilt_only_any_ab": int((~old_any & new_any).sum()),
                    "n_assignment_ineligible": int((~eligible).sum()),
                }
            )
    source_summary = pd.DataFrame(source_rows)
    source_summary.to_csv(table_dir / "tcr_replacement_by_source.csv", index=False)
    expected_whole_any_ab = old_any_ab_total - int(source_summary["n_old_any_ab"].sum()) + int(source_summary["n_rebuilt_any_ab"].sum())
    expected_whole_paired_ab = old_paired_ab_total - int(source_summary["n_old_paired_ab"].sum()) + int(source_summary["n_rebuilt_paired_ab"].sum())
    if expected_whole_any_ab != EXPECTED_WHOLE_ANY_AB or expected_whole_paired_ab != EXPECTED_WHOLE_PAIRED_AB:
        raise RuntimeError("Pre-write whole-atlas TCR totals disagree with frozen overlay audit")

    subprocess.run(
        ["cp", "--reflink=auto", "--sparse=always", "--preserve=mode,timestamps", str(source), str(partial)],
        check=True,
    )
    unchanged_hashes: dict[str, str] = {}
    with h5py.File(partial, "r+") as handle:
        obs = handle["obs"]
        for name in REPLACED_CATEGORICAL:
            values = categorical_values(obs, name)
            unchanged_hashes[name] = value_hash(values[unaffected_mask])
            for positions, side_table in overlays.values():
                values[positions] = replacement_strings(side_table, name)
            write_categorical(obs, name, values)
            print(f"replaced {name}", flush=True)
        for name in CANONICAL_BOOL_MAP:
            values = np.asarray(obs[name][:], dtype=bool)
            unchanged_hashes[name] = hashlib.sha256(values[unaffected_mask].tobytes()).hexdigest()
            for positions, side_table in overlays.values():
                values[positions] = replacement_bool(side_table, name)
            write_boolean(obs, name, values)
            print(f"replaced {name}", flush=True)
        for name in ADD_COLUMNS:
            is_bool = name in GENERAL_BOOL_PROVENANCE or name.endswith(("_umi_available_rebuilt_v2", "_read_available_rebuilt_v2"))
            values = np.zeros(EXPECTED_CELLS, dtype=bool) if is_bool else np.full(EXPECTED_CELLS, "", dtype=object)
            for positions, side_table in overlays.values():
                side_values, side_is_bool = provenance_values(side_table, name)
                if side_is_bool != is_bool:
                    raise RuntimeError(f"Provenance type mismatch for {name}")
                values[positions] = side_values
            if is_bool:
                write_boolean(obs, name, values)
            else:
                write_categorical(obs, name, values)
            print(f"added {name}", flush=True)
        del obs.attrs["column-order"]
        obs.attrs.create(
            "column-order",
            np.asarray([*old_columns, *ADD_COLUMNS], dtype=h5py.string_dtype("utf-8")),
        )
        handle.flush()

    checks: dict[str, bool] = {}
    with h5py.File(partial, "r") as handle:
        obs = handle["obs"]
        observed_columns = [str(v.decode() if isinstance(v, bytes) else v) for v in obs.attrs["column-order"]]
        checks["only_expected_provenance_columns_added"] = observed_columns == [*old_columns, *ADD_COLUMNS]
        checks["protected_structure_sampled_identity"] = protected_signature(handle, protected_obs) == source_signature
        checks["x_shape"] = tuple(handle["X"].attrs["shape"]) == (EXPECTED_CELLS, EXPECTED_GENES)
        checks["x_nnz"] = len(handle["X"]["data"]) == EXPECTED_NNZ
        checks["scvi_shape"] = handle["obsm"]["X_scVI"].shape == (EXPECTED_CELLS, 30)
        checks["umap_shape"] = handle["obsm"]["X_umap"].shape == (EXPECTED_CELLS, 2)
        checks["whole_any_ab_total"] = int(np.asarray(obs["has_any_ab_tcr"][:], dtype=bool).sum()) == EXPECTED_WHOLE_ANY_AB
        checks["whole_paired_ab_total"] = int(np.asarray(obs["has_TRA_TRB_paired"][:], dtype=bool).sum()) == EXPECTED_WHOLE_PAIRED_AB
        covered = np.asarray(obs["tcr_sidecar_covered_v2"][:], dtype=bool)
        eligible = np.asarray(obs["tcr_assignment_eligible_rebuilt_v2"][:], dtype=bool)
        checks["sidecar_coverage_exact"] = int(covered.sum()) == EXPECTED_AFFECTED_ROWS
        checks["ineligible_atlas_rows_exact"] = int((covered & ~eligible).sum()) == EXPECTED_INELIGIBLE_ATLAS_ROWS
        unchanged_ok = True
        for name in REPLACED_CATEGORICAL:
            values = categorical_values(obs, name)
            unchanged_ok &= value_hash(values[unaffected_mask]) == unchanged_hashes[name]
        for name in CANONICAL_BOOL_MAP:
            values = np.asarray(obs[name][:], dtype=bool)
            unchanged_ok &= hashlib.sha256(values[unaffected_mask].tobytes()).hexdigest() == unchanged_hashes[name]
        checks["unaffected_tcr_values_exact"] = bool(unchanged_ok)

        affected_exact = True
        for positions, side_table in overlays.values():
            for name in REPLACED_CATEGORICAL:
                affected_exact &= np.array_equal(categorical_values_at(obs, name, positions), replacement_strings(side_table, name))
            for name in CANONICAL_BOOL_MAP:
                affected_exact &= np.array_equal(np.asarray(obs[name][positions], dtype=bool), replacement_bool(side_table, name))
        checks["affected_tcr_values_exact"] = bool(affected_exact)

        ineligible_positions = []
        for positions, side_table in overlays.values():
            eligible_side = np.asarray(side_table["tcr_assignment_eligible"].to_numpy(zero_copy_only=False), dtype=bool)
            ineligible_positions.append(positions[~eligible_side])
        ineligible_positions_array = np.concatenate(ineligible_positions)
        cleared = True
        for chain in CHAINS:
            cleared &= bool(np.all(categorical_values_at(obs, f"{chain}_cdr3", ineligible_positions_array) == ""))
            cleared &= bool(np.all(~np.asarray(obs[f"has_{chain}"][ineligible_positions_array], dtype=bool)))
        checks["ineligible_rows_fail_closed"] = cleared
        checks["canonical_flag_consistency_affected"] = all(
            np.array_equal(
                np.asarray(obs[f"has_{chain}"][affected_positions], dtype=bool),
                categorical_values_at(obs, f"{chain}_cdr3", affected_positions) != "",
            )
            for chain in CHAINS
        )

    backed = ad.read_h5ad(partial, backed="r")
    checks["anndata_backed_open"] = backed.shape == (EXPECTED_CELLS, EXPECTED_GENES)
    backed.file.close()
    source_after = file_state(source, include_sha=True)
    checks["metadata_source_state_and_sha_unchanged"] = source_after == source_before
    checks["canonical_symlink_unchanged"] = (
        ROOT / "Integrated_dataset/integrated_full_atlas.h5ad"
    ).resolve() == EXPECTED_CANONICAL_TARGET
    status = "PASS_TCR_H5AD_CANDIDATE" if all(checks.values()) else "FAIL_TCR_H5AD_CANDIDATE"
    manifest: dict[str, Any] = {
        "status": status,
        "source_metadata_candidate": source_before,
        "sidecar": {"path": str(sidecar), "sha256": sidecar_sha, "rows": EXPECTED_SIDECAR_ROWS},
        "candidate": str(candidate),
        "n_cells": EXPECTED_CELLS,
        "n_affected_sources": len(overlays),
        "n_affected_atlas_rows": affected_rows,
        "n_assignment_ineligible_atlas_rows": EXPECTED_INELIGIBLE_ATLAS_ROWS,
        "n_whole_atlas_any_ab_tcr": EXPECTED_WHOLE_ANY_AB,
        "n_whole_atlas_paired_tra_trb": EXPECTED_WHOLE_PAIRED_AB,
        "n_replaced_columns": len(REPLACED_COLUMNS),
        "n_added_provenance_columns": len(ADD_COLUMNS),
        "checks": checks,
        "canonical_symlink_switched": False,
        "elapsed_seconds": round(time.time() - started, 3),
    }
    manifest_path = log_dir / "tcr_h5ad_candidate_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if status != "PASS_TCR_H5AD_CANDIDATE":
        print(json.dumps(manifest, indent=2, sort_keys=True))
        raise SystemExit(1)
    os.replace(partial, candidate)
    manifest["candidate_size_bytes"] = candidate.stat().st_size
    manifest["candidate_sha256"] = sha256_file(candidate)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (log_dir / "tcr_h5ad_candidate_summary.md").write_text(
        "# Full-atlas TCR-corrected candidate\n\n"
        f"- Status: `{status}`\n"
        f"- Candidate: `{candidate}`\n"
        f"- Candidate SHA-256: `{manifest['candidate_sha256']}`\n"
        f"- Affected atlas rows: {affected_rows:,} across {len(overlays)} sources\n"
        f"- Whole-atlas productive alpha-beta cells: {EXPECTED_WHOLE_ANY_AB:,}\n"
        f"- Whole-atlas paired TRA/TRB cells: {EXPECTED_WHOLE_PAIRED_AB:,}\n"
        "- Canonical atlas changed: no\n"
        "- Next action: regenerate TCR-derived truth/control and boundary audits\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
