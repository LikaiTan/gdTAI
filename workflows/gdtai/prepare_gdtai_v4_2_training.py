#!/usr/bin/env python3
"""Freeze V4.2 features/folds and build the labeled-cell expression cache."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TRUTH_DIR = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_ground_truth"
OUT_DIR = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_training"
CACHE_DIR = ROOT / "Integrated_dataset/cache/gdT_prediction/gdtai_v4_2_training"
LOG_DIR = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_training"
ATLAS = Path("/ssd/tnk_phase3/Integrated_dataset/full_atlas/tcr_corrected/integrated_full_atlas.h5ad")
ATLAS_SHA = "d32c9d2bdb955b12e1eafbed8322f8cb965cf3a225191e612b53f3d3783480d5"
STAGE1_GENES = {
    "BCL11B", "CCL5", "CCR7", "CD2", "CD247", "CD3D", "CD3E", "CD3G", "CD4", "CD5", "CD6",
    "CD8A", "CD8B", "CTLA4", "CTSW", "EOMES", "FCER1G", "FCGR3A", "FOXP3", "GNLY", "GZMB",
    "IKZF2", "IL2RA", "IL7R", "KLRB1", "KLRC1", "KLRD1", "KLRF1", "LAT", "LCK", "LEF1", "LTB",
    "MAL", "NKG7", "PRF1", "RUNX3", "SELL", "TBX21", "TCF7", "THEMIS", "TIGIT", "TRAT1", "TYROBP",
    "ZNF683", "TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2",
}
POSITIVE_SOURCES = ["HRA005041", "GSE144469", "MalteGDT"]
NK_SOURCES = ["GSE125527", "GSE228597", "GSE212217", "GSE234069", "GSE235863", "GSE243013", "GSE287541"]


def decode(value: object) -> str:
    return value.decode() if isinstance(value, bytes) else str(value)


def read_var_names(handle: h5py.File) -> list[str]:
    var = handle["var"]
    key = decode(var.attrs.get("_index", "_index"))
    node = var[key]
    return [decode(value) for value in node[:]]


def matrix_group(handle: h5py.File, key: str) -> h5py.Group:
    node = handle[key]
    if not isinstance(node, h5py.Group) or not {"data", "indices", "indptr"}.issubset(node):
        raise TypeError("Expected CSR matrix")
    return node


def extract_csr_rows(group: h5py.Group, rows: np.ndarray, source_to_output: np.ndarray,
                     output_feature_count: int, row_chunk: int) -> np.ndarray:
    requested = np.asarray(rows, dtype=np.int64)
    output = np.zeros((len(requested), output_feature_count), dtype=np.float32)
    indptr = np.asarray(group["indptr"][:], dtype=np.int64)
    request_to_output = np.full(len(indptr) - 1, -1, dtype=np.int64)
    request_to_output[requested] = np.arange(len(requested))
    totals = np.full(len(requested), np.nan, dtype=np.float64)
    for start in range(0, len(indptr) - 1, row_chunk):
        stop = min(start + row_chunk, len(indptr) - 1)
        requested_output = request_to_output[start:stop]
        if not (requested_output >= 0).any():
            continue
        data_start, data_stop = int(indptr[start]), int(indptr[stop])
        values = np.asarray(group["data"][data_start:data_stop], dtype=np.float64)
        columns = np.asarray(group["indices"][data_start:data_stop], dtype=np.int64)
        offsets = indptr[start:stop + 1] - data_start
        cumulative = np.concatenate([[0.0], np.cumsum(values, dtype=np.float64)])
        row_totals = cumulative[offsets[1:]] - cumulative[offsets[:-1]]
        selected_local = np.flatnonzero(requested_output >= 0)
        totals[requested_output[selected_local]] = row_totals[selected_local]
        mapped = source_to_output[columns]
        selected_positions = np.flatnonzero(mapped >= 0)
        if selected_positions.size:
            local_rows = np.searchsorted(offsets[1:], selected_positions, side="right")
            cache_rows = requested_output[local_rows]
            keep = cache_rows >= 0
            np.add.at(output, (cache_rows[keep], mapped[selected_positions[keep]]), values[selected_positions[keep]])
    if not np.isfinite(totals).all() or (totals <= 0).any():
        raise RuntimeError("Invalid library sizes in requested rows")
    return np.log1p(output * (10_000.0 / totals[:, None])).astype(np.float32)


def sha256_file(path: Path, chunk: int = 64 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def assign_balanced_groups(frame: pd.DataFrame, folds: int = 3) -> dict[str, int]:
    """Greedily balance whole groups within source/truth strata."""
    assignments: dict[str, int] = {}
    grouped = frame.groupby(["source_gse_id", "truth_class", "group_id"], observed=True).size().reset_index(name="n")
    for (_, _), stratum in grouped.groupby(["source_gse_id", "truth_class"], observed=True):
        loads = np.zeros(folds, dtype=np.int64)
        for row in stratum.sort_values(["n", "group_id"], ascending=[False, True]).itertuples(index=False):
            fold = int(np.argmin(loads))
            assignments[str(row.group_id)] = fold
            loads[fold] += int(row.n)
    return assignments


def freeze_manifests() -> dict:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    labels_path = TRUTH_DIR / "v4_2_label_manifest.parquet"
    labels = pd.read_parquet(labels_path)
    old_features = pd.read_csv(ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_preflight/feature_manifest.csv")
    features = old_features.sort_values("feature_index").copy()
    features["stage1"] = features["gene"].isin(STAGE1_GENES)
    features["stage2"] = True
    if len(features) != 197 or int(features.stage1.sum()) != 50:
        raise RuntimeError(f"Unexpected feature dimensions: {len(features)} total, {features.stage1.sum()} Stage 1")
    features.to_csv(OUT_DIR / "feature_manifest.csv", index=False)

    development = labels[labels.allow_fit].copy()
    assignments = assign_balanced_groups(development)
    development["inner_fold"] = development.group_id.map(assignments).astype(int)
    splits = (development.groupby(["source_gse_id", "truth_class", "group_id", "inner_fold"], observed=True)
              .size().reset_index(name="n_cells"))
    splits.to_csv(OUT_DIR / "grouped_split_manifest.csv", index=False)
    outer = pd.DataFrame(
        [{"outer_fold_id": f"positive_transfer::{source}", "heldout_axis": "positive_source", "heldout_source": source} for source in POSITIVE_SOURCES]
        + [{"outer_fold_id": f"nk_transfer::{source}", "heldout_axis": "nk_source", "heldout_source": source} for source in NK_SOURCES]
    )
    outer.to_csv(OUT_DIR / "outer_fold_manifest.csv", index=False)
    balance = (development.groupby(["inner_fold", "source_gse_id", "truth_class"], observed=True).size()
               .reset_index(name="n_cells"))
    balance.to_csv(OUT_DIR / "inner_fold_balance.csv", index=False)
    lockbox = labels[~labels.allow_fit].copy()
    lockbox.to_parquet(OUT_DIR / "locked_evaluation_manifest.parquet", index=False, compression="zstd")

    checks = pd.DataFrame([
        {"check": "feature_count_le_300", "pass": len(features) <= 300, "detail": str(len(features))},
        {"check": "stage1_feature_count", "pass": int(features.stage1.sum()) == 50, "detail": str(int(features.stage1.sum()))},
        {"check": "no_score_features", "pass": not features.gene.str.contains("score", case=False).any(), "detail": "individual genes only"},
        {"check": "three_positive_outer_folds", "pass": int((outer.heldout_axis == "positive_source").sum()) == 3, "detail": ",".join(POSITIVE_SOURCES)},
        {"check": "seven_nk_transfer_folds", "pass": int((outer.heldout_axis == "nk_source").sum()) == 7, "detail": ",".join(NK_SOURCES)},
        {"check": "lockbox_has_no_fit", "pass": not lockbox.allow_fit.any(), "detail": str(len(lockbox))},
        {"check": "source_namespaced_groups", "pass": development.apply(lambda r: str(r.group_id).startswith(str(r.source_gse_id) + "::"), axis=1).all(), "detail": "source::donor/sample/library"},
    ])
    checks.to_csv(OUT_DIR / "precache_checks.csv", index=False)
    if not checks["pass"].all():
        raise RuntimeError("Feature/fold freeze failed")
    result = {
        "status": "PASS_FEATURE_AND_FOLD_FREEZE",
        "label_manifest": str(labels_path),
        "label_manifest_sha256": sha256_file(labels_path),
        "feature_manifest": str(OUT_DIR / "feature_manifest.csv"),
        "feature_manifest_sha256": sha256_file(OUT_DIR / "feature_manifest.csv"),
        "split_manifest": str(OUT_DIR / "grouped_split_manifest.csv"),
        "split_manifest_sha256": sha256_file(OUT_DIR / "grouped_split_manifest.csv"),
        "n_labels": len(labels), "n_development": len(development), "n_lockbox": len(lockbox),
        "n_features": len(features), "n_stage1_features": int(features.stage1.sum()),
    }
    (LOG_DIR / "freeze_summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def build_cache(row_chunk: int, force: bool) -> dict:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    manifest_path = CACHE_DIR / "feature_cache_manifest.json"
    if manifest_path.exists() and not force:
        return json.loads(manifest_path.read_text())
    labels_path = TRUTH_DIR / "v4_2_label_manifest.parquet"
    labels = pd.read_parquet(labels_path)
    features = pd.read_csv(OUT_DIR / "feature_manifest.csv").sort_values("feature_index")
    feature_names = features.gene.astype(str).tolist()
    rows = labels.atlas_row.to_numpy(np.int64)
    order = np.argsort(rows, kind="mergesort")
    if np.unique(rows).size != len(rows):
        raise RuntimeError("Cache requires one manifest row per atlas cell")
    temporary = CACHE_DIR / "gene_features.tmp.npy"
    final = CACHE_DIR / "gene_features.npy"
    if temporary.exists():
        temporary.unlink()
    output = np.lib.format.open_memmap(temporary, mode="w+", dtype=np.float32, shape=(len(labels), len(features)))
    with h5py.File(ATLAS, "r") as handle:
        var_names = read_var_names(handle)
        lookup = {gene: index for index, gene in enumerate(var_names)}
        missing = [gene for gene in feature_names if gene not in lookup]
        if missing:
            raise RuntimeError(f"Missing frozen features: {missing}")
        mapping = np.full(len(var_names), -1, dtype=np.int32)
        for output_index, gene in enumerate(feature_names):
            mapping[lookup[gene]] = output_index
        values = extract_csr_rows(matrix_group(handle, "X"), rows[order], mapping, len(features), row_chunk)
    output[order] = values
    output.flush(); del output, values
    os.replace(temporary, final)
    cache = np.load(final, mmap_mode="r")
    coverage = pd.DataFrame({"gene": feature_names, "detected_cells": np.asarray((cache > 0).sum(axis=0)).ravel()})
    coverage["detected_fraction"] = coverage.detected_cells / len(labels)
    coverage.to_csv(OUT_DIR / "feature_detection.csv", index=False)
    summary = {
        "status": "PASS_FEATURE_CACHE",
        "atlas_path": str(ATLAS), "atlas_sha256": ATLAS_SHA,
        "label_manifest_sha256": sha256_file(labels_path),
        "feature_manifest_sha256": sha256_file(OUT_DIR / "feature_manifest.csv"),
        "matrix_path": str(final), "matrix_sha256": sha256_file(final),
        "shape": list(cache.shape), "dtype": str(cache.dtype), "normalization": "log1p(CP10K)",
        "h5ad_modified": False,
    }
    manifest_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["freeze", "cache", "all"], default="all")
    parser.add_argument("--row-chunk", type=int, default=20_000)
    parser.add_argument("--force-cache", action="store_true")
    args = parser.parse_args()
    result = {}
    if args.stage in {"freeze", "all"}:
        result["freeze"] = freeze_manifests()
    if args.stage in {"cache", "all"}:
        result["cache"] = build_cache(args.row_chunk, args.force_cache)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
