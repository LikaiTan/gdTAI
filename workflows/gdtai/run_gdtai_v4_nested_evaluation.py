#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Run the precommitted gdTAI V4 nested cross-dataset evaluation.

Project H5ADs are opened read-only. Project-data cache construction and model
fitting require a checksum-bound Step-2 approval record created only after the
Step-1 PASS package has been reviewed. ``--stage validate`` is deliberately
non-fitting and can be used before that approval.
"""

# ruff: noqa: E402

from __future__ import annotations

# TNK_WORKFLOW_BOOTSTRAP
import sys as _tnk_sys
from pathlib import Path as _TnkPath

_TNK_PROJECT_ROOT = _TnkPath(__file__).resolve().parents[2]
for _tnk_path in (_TNK_PROJECT_ROOT, _TNK_PROJECT_ROOT / "src", _TNK_PROJECT_ROOT / "workflows" / "gdtai"):
    _tnk_value = str(_tnk_path)
    if _tnk_value not in _tnk_sys.path:
        _tnk_sys.path.insert(0, _tnk_value)

import argparse
import hashlib
import html
import json
import logging
import math
import os
import pickle
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import h5py
from joblib import Parallel, delayed, parallel_config
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import expit
from sklearn.metrics import average_precision_score

from gdtai_v4_nested_core import (
    DERIVED_FEATURE_ORDER,
    FittedBinaryEstimator,
    PlattCalibrator,
    ThresholdDecision,
    apply_two_stage_call,
    balanced_training_weights,
    binary_metrics,
    bootstrap_interval,
    candidate_id,
    derive_features,
    exclusion_flags,
    fit_binary_estimator,
    fit_platt_calibrator,
    grouped_confusion_counts,
    paired_hierarchical_bootstrap_f1_difference,
    parameter_grid,
    select_stage1_threshold,
    select_stage2_threshold,
    stable_seed,
)
from compare_frozen_gdtai_profiles import (  # noqa: E402
    CANONICAL_OBS_COLUMNS as EXTENSION_OBS_COLUMNS,
    FeatureSpec as ExtensionFeatureSpec,
    normalize_annotation as normalize_extension_annotation,
    read_obs as read_extension_obs,
    truth_frame as extension_truth_frame,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = PROJECT_ROOT / "configs/models/gdtai/v4_nested_evaluation.json"
PREFLIGHT_CONFIG = PROJECT_ROOT / "configs/models/gdtai/v4_preflight.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--stage",
        choices=("validate", "cache", "evaluate", "report", "all"),
        default="validate",
        help="validate never fits; all other stages require the Step-2 approval record",
    )
    parser.add_argument("--outer-fold", action="append", default=[], help="Optional held-out source filter")
    parser.add_argument("--candidate-jobs", type=int, default=1)
    parser.add_argument("--matrix-row-chunk", type=int, default=20000)
    parser.add_argument("--skip-source-hashes", action="store_true", help="Development only; evaluate refuses this")
    parser.add_argument("--force-cache", action="store_true")
    parser.add_argument("--no-pdf", action="store_true")
    parser.add_argument("--log-level", default="INFO")
    return parser.parse_args()


def resolve(path: str | Path) -> Path:
    value = Path(path)
    return value if value.is_absolute() else PROJECT_ROOT / value


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(value: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def sha256_file(path: Path, chunk_bytes: int = 64 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_bytes):
            digest.update(chunk)
    return digest.hexdigest()


def decode(value: Any) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def read_strings(dataset: h5py.Dataset) -> np.ndarray:
    try:
        return np.asarray(dataset.asstr()[:], dtype=object)
    except Exception:
        return np.asarray([decode(value) for value in dataset[:]], dtype=object)


def read_axis_column(group: h5py.Group, column: str) -> np.ndarray:
    obj = group[column]
    if isinstance(obj, h5py.Group):
        encoding = decode(obj.attrs.get("encoding-type", ""))
        if encoding == "categorical":
            categories = read_strings(obj["categories"])
            codes = np.asarray(obj["codes"][:], dtype=np.int64)
            output = np.full(codes.shape, "", dtype=object)
            valid = codes >= 0
            output[valid] = categories[codes[valid]]
            return output
        if encoding in {"nullable-boolean", "nullable-integer"}:
            values = np.asarray(obj["values"][:])
            mask = np.asarray(obj["mask"][:], dtype=bool)
            output = values.astype(object)
            output[mask] = np.nan
            return output
        raise TypeError(f"Unsupported encoded H5AD column: {column} ({encoding})")
    values = obj[:]
    return read_strings(obj) if values.dtype.kind in {"S", "O", "U"} else np.asarray(values)


def axis_index_key(group: h5py.Group) -> str:
    value = group.attrs.get("_index")
    if value is not None:
        return decode(value)
    if "_index" in group:
        return "_index"
    raise KeyError("H5AD axis has no index")


def read_var_names(handle: h5py.File) -> list[str]:
    return [str(value) for value in read_axis_column(handle["var"], axis_index_key(handle["var"]))]


def matrix_group(handle: h5py.File, matrix_key: str) -> h5py.Group:
    obj: h5py.Group | h5py.Dataset = handle
    for token in matrix_key.strip("/").split("/"):
        obj = obj[token]  # type: ignore[index]
    if not isinstance(obj, h5py.Group) or not {"data", "indices", "indptr"}.issubset(obj.keys()):
        raise TypeError(f"{matrix_key} must be a CSR sparse matrix")
    return obj


def setup_logging(log_dir: Path, level: str) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    handlers = [
        logging.FileHandler(log_dir / "nested_evaluation.log", mode="a", encoding="utf-8"),
        logging.StreamHandler(),
    ]
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
        force=True,
    )


def validate_approval(config: Mapping[str, Any]) -> dict[str, Any]:
    path = resolve(config["step2_approval_file"])
    if not path.exists():
        raise PermissionError(
            f"Step 2 is not authorized. Missing checksum-bound approval record: {path}"
        )
    approval = load_json(path)
    required = {
        "approved": True,
        "protocol_version": str(config["protocol_version"]),
        "cell_manifest_sha256": str(config["preflight"]["cell_manifest_sha256"]),
        "split_manifest_sha256": str(config["preflight"]["split_manifest_sha256"]),
        "cd4_treg_rule_sha256": str(config["preflight"]["cd4_treg_rule_sha256"]),
        "evaluation_config_sha256": sha256_file(DEFAULT_CONFIG),
        "runner_sha256": sha256_file(Path(__file__)),
        "core_sha256": sha256_file(Path(__file__).with_name("gdtai_v4_nested_core.py")),
    }
    for key, expected in required.items():
        if approval.get(key) != expected:
            raise PermissionError(f"Step-2 approval field {key!r} does not match the frozen gate")
    if not str(approval.get("approved_by", "")).strip() or not str(approval.get("approved_at", "")).strip():
        raise PermissionError("Step-2 approval must identify approver and approval time")
    return approval


def validate_contract(config: Mapping[str, Any], verify_source_hashes: bool) -> dict[str, Any]:
    preflight = config["preflight"]
    summary_path = resolve(preflight["summary"])
    summary = load_json(summary_path)
    failures: list[str] = []
    if summary.get("overall_status") != "PASS":
        failures.append("Step-1 preflight status is not PASS")
    if summary.get("protocol_version") != config["protocol_version"]:
        failures.append("Protocol version disagrees with Step-1 summary")
    if summary.get("training_started") is not False:
        failures.append("Step-1 summary unexpectedly records training")
    if int(summary.get("n_preflight_failures", -1)) != 0 or int(summary.get("n_preflight_warnings", -1)) != 0:
        failures.append("Step-1 summary contains failures or warnings")
    digest_paths = {
        "cell_manifest_sha256": resolve(preflight["cell_manifest"]),
        "split_manifest_sha256": resolve(preflight["split_manifest"]),
        "feature_manifest_sha256": resolve(preflight["feature_manifest"]),
    }
    actual_digests: dict[str, str] = {}
    for key, path in digest_paths.items():
        actual = sha256_file(path)
        actual_digests[key] = actual
        if actual != str(preflight[key]):
            failures.append(f"Checksum mismatch: {path}")
    preflight_config = load_json(PREFLIGHT_CONFIG)
    frozen_rule_payload = {
        "cd4_helper_rule": preflight_config["cd4_helper_rule"],
        "treg_rule": preflight_config["treg_rule"],
    }
    actual_rule_hash = hashlib.sha256(
        json.dumps(frozen_rule_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    if summary.get("cd4_treg_rule_sha256") != preflight["cd4_treg_rule_sha256"] or actual_rule_hash != preflight["cd4_treg_rule_sha256"]:
        failures.append("CD4/Treg rule hash disagrees with Step-1 summary")
    expected_cd4 = {
        "marker": "CD4",
        "marker_min": float(preflight_config["cd4_helper_rule"]["cd4_min"]),
        "support_genes": preflight_config["cd4_helper_genes"][1:],
        "support_detected_min": int(preflight_config["cd4_helper_rule"]["support_detected_min"]),
        "program_genes": preflight_config["cd4_helper_genes"],
        "program_mean_min": float(preflight_config["cd4_helper_rule"]["program_mean_min"]),
    }
    expected_treg = {
        "marker": "FOXP3",
        "marker_min": float(preflight_config["treg_rule"]["foxp3_min"]),
        "support_genes": preflight_config["treg_genes"][1:],
        "support_detected_min": int(preflight_config["treg_rule"]["support_detected_min"]),
        "program_genes": preflight_config["treg_genes"],
        "program_mean_min": float(preflight_config["treg_rule"]["program_mean_min"]),
    }
    if config["cd4_helper_rule"] != expected_cd4 or config["treg_rule"] != expected_treg:
        failures.append("Expanded Step-2 CD4/Treg rule is not semantically identical to Step 1")

    features = pd.read_csv(resolve(preflight["feature_manifest"]))
    if features.shape[0] != int(config["feature_policy"]["expected_gene_count"]):
        failures.append("Frozen feature count changed")
    if features["gene"].duplicated().any():
        failures.append("Frozen feature manifest contains duplicates")
    forbidden = [gene for gene in features["gene"].astype(str) if gene.startswith(("MS4A", "CD79", "LST1", "LYZ"))]
    if forbidden:
        failures.append(f"Forbidden B/myeloid features found: {forbidden}")

    cells = pd.read_csv(
        resolve(preflight["cell_manifest"]),
        usecols=["cell_id", "source_gse_id", "truth_class", "stage1_role", "stage2_role", "group_key"],
    )
    splits = pd.read_csv(resolve(preflight["split_manifest"]))
    if cells["cell_id"].duplicated().any():
        failures.append("Cell-label manifest has duplicate cell IDs")
    if set(cells.loc[cells["stage2_role"].ne("none"), "truth_class"]) & {
        "gdT_silver",
        "sorted_sensitivity",
        "dual_or_ambiguous",
    }:
        failures.append("Sensitivity or ambiguous truth entered Stage 2")
    leakage = splits.groupby(["outer_fold_id", "stage", "group_key"])["inner_fold"].nunique().max()
    if int(leakage) != 1:
        failures.append("A protected group crosses inner folds")
    if set(splits["inner_fold"].unique()) != {0, 1, 2}:
        failures.append("Frozen inner-fold IDs changed")

    source_hash_rows: list[dict[str, Any]] = []
    if verify_source_hashes:
        inputs = pd.read_csv(resolve(preflight["input_manifest"]))
        for row in inputs.itertuples(index=False):
            path = Path(row.path)
            actual = sha256_file(path)
            passed = actual == str(row.sha256)
            source_hash_rows.append(
                {"input_id": row.input_id, "path": str(path), "expected_sha256": row.sha256, "actual_sha256": actual, "pass": passed}
            )
            if not passed:
                failures.append(f"Source input checksum changed: {row.input_id}")
    return {
        "status": "PASS" if not failures else "FAIL",
        "failures": failures,
        "protocol_version": config["protocol_version"],
        "n_cells": int(cells.shape[0]),
        "n_features": int(features.shape[0]),
        "n_split_groups": int(splits.shape[0]),
        "actual_digests": actual_digests,
        "source_hashes_verified": bool(verify_source_hashes),
        "source_hash_rows": source_hash_rows,
    }


def extract_csr_rows(
    group: h5py.Group,
    rows: np.ndarray,
    source_to_output: np.ndarray,
    output_feature_count: int,
    matrix_state: str,
    target_sum: float,
    row_chunk: int,
) -> np.ndarray:
    """Extract selected genes with bounded contiguous CSR row-block scans."""
    requested = np.asarray(rows, dtype=np.int64)
    if requested.size and (np.diff(requested) < 0).any():
        raise ValueError("CSR rows must be sorted")
    if output_feature_count <= 0:
        raise ValueError("Output feature count must be positive")
    mapped_columns = source_to_output[source_to_output >= 0]
    if mapped_columns.size and int(mapped_columns.max()) >= output_feature_count:
        raise ValueError("Source-to-output mapping exceeds the requested feature width")
    output = np.zeros((requested.size, output_feature_count), dtype=np.float32)
    indptr = np.asarray(group["indptr"][:], dtype=np.int64)
    indices_ds = group["indices"]
    data_ds = group["data"]
    n_obs = indptr.size - 1
    if requested.size and (requested[0] < 0 or requested[-1] >= n_obs):
        raise IndexError("Requested CSR row is outside matrix dimensions")
    request_to_output = np.full(n_obs, -1, dtype=np.int64)
    request_to_output[requested] = np.arange(requested.size, dtype=np.int64)
    totals = np.full(requested.size, np.nan, dtype=np.float64)
    for block_start in range(0, n_obs, row_chunk):
        block_stop = min(block_start + row_chunk, n_obs)
        requested_output = request_to_output[block_start:block_stop]
        if not (requested_output >= 0).any():
            continue
        data_start = int(indptr[block_start])
        data_stop = int(indptr[block_stop])
        values = np.asarray(data_ds[data_start:data_stop], dtype=np.float64)
        source_columns = np.asarray(indices_ds[data_start:data_stop], dtype=np.int64)
        offsets = indptr[block_start : block_stop + 1] - data_start
        requested_local_rows = np.flatnonzero(requested_output >= 0)
        requested_cache_rows = requested_output[requested_local_rows]
        if matrix_state == "raw_counts":
            cumulative = np.empty(values.size + 1, dtype=np.float64)
            cumulative[0] = 0.0
            np.cumsum(values, dtype=np.float64, out=cumulative[1:])
            row_totals = cumulative[offsets[1:]] - cumulative[offsets[:-1]]
            totals[requested_cache_rows] = row_totals[requested_local_rows]
        elif matrix_state != "log1p_cp10k_registered":
            raise ValueError(f"Unsupported matrix state: {matrix_state}")
        mapped_columns = source_to_output[source_columns]
        selected_positions = np.flatnonzero(mapped_columns >= 0)
        if selected_positions.size:
            selected_local_rows = np.searchsorted(offsets[1:], selected_positions, side="right")
            selected_cache_rows = requested_output[selected_local_rows]
            selected_requested = selected_cache_rows >= 0
            if selected_requested.any():
                np.add.at(
                    output,
                    (
                        selected_cache_rows[selected_requested],
                        mapped_columns[selected_positions[selected_requested]],
                    ),
                    values[selected_positions[selected_requested]].astype(np.float32),
                )
    if matrix_state == "raw_counts":
        if not np.isfinite(totals).all() or (totals <= 0).any():
            bad = int((~np.isfinite(totals) | (totals <= 0)).sum())
            raise ValueError(f"{bad} requested raw rows have invalid library sizes")
        output = np.log1p(output * (target_sum / totals[:, None])).astype(np.float32)
    return output


def apply_frozen_exclusions(
    cells: pd.DataFrame,
    gene_matrix: np.ndarray,
    feature_names: Sequence[str],
    config: Mapping[str, Any],
) -> dict[str, Any]:
    """Reproduce the preflight scope, then apply frozen exclusions to all rows."""
    cd4, treg, union = exclusion_flags(
        gene_matrix,
        feature_names,
        config["cd4_helper_rule"],
        config["treg_rule"],
    )
    preflight_scope = (
        cells["stage1_role"].eq("nk_negative")
        | cells["truth_class"].isin(
            ["gdT_primary", "gdT_silver", "sorted_gdT_weak", "sorted_sensitivity"]
        )
    ).to_numpy(dtype=bool)
    comparisons = (
        ("cd4_helper_exclusion", cd4),
        ("treg_exclusion", treg),
        ("exclusion_union", union),
    )
    for column, recomputed in comparisons:
        frozen = cells[column].to_numpy(dtype=bool)
        if not np.array_equal(recomputed[preflight_scope], frozen[preflight_scope]):
            mismatch = int((recomputed[preflight_scope] != frozen[preflight_scope]).sum())
            raise RuntimeError(
                f"Recomputed {column} disagrees with the frozen preflight scope for {mismatch} cells"
            )
        cells[column] = recomputed
    return {
        "preflight_scope_cells": int(preflight_scope.sum()),
        "preflight_scope_reproduction_pass": True,
        "all_cell_cd4_helper_excluded": int(cd4.sum()),
        "all_cell_treg_excluded": int(treg.sum()),
        "all_cell_union_excluded": int(union.sum()),
    }


def extract_feature_cache(config: Mapping[str, Any], row_chunk: int, force: bool) -> dict[str, Any]:
    cache_dir = resolve(config["outputs"]["cache_dir"])
    cache_dir.mkdir(parents=True, exist_ok=True)
    summary_path = cache_dir / "feature_cache_manifest.json"
    if summary_path.exists() and not force:
        summary = load_json(summary_path)
        expected = config["preflight"]["cell_manifest_sha256"]
        if summary.get("cell_manifest_sha256") == expected:
            logging.info("Reusing checksum-compatible feature cache")
            return summary
        raise RuntimeError("Existing feature cache is incompatible; use --force-cache after review")

    cells = pd.read_csv(resolve(config["preflight"]["cell_manifest"]))
    features = pd.read_csv(resolve(config["preflight"]["feature_manifest"])).sort_values("feature_index")
    feature_names = features["gene"].astype(str).tolist()
    inputs = pd.read_csv(resolve(config["preflight"]["input_manifest"]))
    preflight_config = load_json(PREFLIGHT_CONFIG)
    matrix_states = preflight_config["matrix_states"]
    matrix_path = cache_dir / "gene_features.npy"
    temporary = cache_dir / "gene_features.npy.tmp"
    if temporary.exists():
        temporary.unlink()
    cache = np.lib.format.open_memmap(
        temporary,
        mode="w+",
        dtype=np.float32,
        shape=(cells.shape[0], len(feature_names)),
    )
    cache[:] = 0
    extraction_rows: list[dict[str, Any]] = []
    for input_id in sorted(cells["expression_input_id"].unique()):
        input_rows = inputs[inputs["input_id"].eq(input_id)]
        if input_rows.shape[0] != 1:
            raise RuntimeError(f"Input manifest does not uniquely define {input_id}")
        row = input_rows.iloc[0]
        subset_index = cells.index[cells["expression_input_id"].eq(input_id)].to_numpy(dtype=np.int64)
        expression_rows = cells.loc[subset_index, "expression_row"].to_numpy(dtype=np.int64)
        order = np.argsort(expression_rows, kind="mergesort")
        sorted_expression_rows = expression_rows[order]
        sorted_manifest_rows = subset_index[order]
        if np.unique(sorted_expression_rows).size != sorted_expression_rows.size:
            raise RuntimeError(f"Duplicate expression rows in {input_id}")
        with h5py.File(Path(row.path), "r") as handle:
            var_names = read_var_names(handle)
            lookup = {gene: index for index, gene in enumerate(var_names)}
            if len(lookup) != len(var_names):
                raise RuntimeError(f"Duplicate var names in {input_id}")
            mapping = np.full(len(var_names), -1, dtype=np.int32)
            present = 0
            for output_index, gene in enumerate(feature_names):
                if gene in lookup:
                    mapping[lookup[gene]] = output_index
                    present += 1
            group = matrix_group(handle, str(row.matrix_key))
            n_obs = int(group["indptr"].shape[0] - 1)
            if sorted_expression_rows.size and int(sorted_expression_rows[-1]) >= n_obs:
                raise IndexError(f"Expression row exceeds {input_id} dimensions")
            values = extract_csr_rows(
                group,
                sorted_expression_rows,
                mapping,
                len(feature_names),
                str(matrix_states[input_id]),
                float(config["normalization"]["target_sum"]),
                row_chunk,
            )
        cache[sorted_manifest_rows] = values
        extraction_rows.append(
            {
                "input_id": input_id,
                "n_rows": int(sorted_manifest_rows.size),
                "features_present": int(present),
                "features_expected": int(len(feature_names)),
                "coverage": float(present / len(feature_names)),
                "matrix_state": matrix_states[input_id],
            }
        )
        logging.info("Cached %s rows from %s", f"{sorted_manifest_rows.size:,}", input_id)
    cache.flush()
    del cache
    os.replace(temporary, matrix_path)
    gene_matrix = np.load(matrix_path, mmap_mode="r")
    derived, derived_names = derive_features(
        gene_matrix,
        feature_names,
        config["feature_policy"]["family_prefixes"],
        config["feature_policy"]["derived_panels"],
    )
    derived_path = cache_dir / "derived_features.npy"
    np.save(derived_path, derived, allow_pickle=False)
    legacy_score_path = cache_dir / "legacy_trd_minus_trab.npy"
    legacy_score = extract_legacy_score(cells, config)
    np.save(legacy_score_path, legacy_score, allow_pickle=False)
    exclusion_audit = apply_frozen_exclusions(cells, gene_matrix, feature_names, config)
    extraction_table = resolve(config["outputs"]["table_dir"]) / "feature_cache_inputs.csv"
    extraction_table.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(extraction_rows).to_csv(extraction_table, index=False)
    summary = {
        "protocol_version": config["protocol_version"],
        "cell_manifest_sha256": config["preflight"]["cell_manifest_sha256"],
        "feature_manifest_sha256": config["preflight"]["feature_manifest_sha256"],
        "n_cells": int(cells.shape[0]),
        "n_gene_features": int(len(feature_names)),
        "n_derived_features": int(len(derived_names)),
        "gene_feature_path": str(matrix_path),
        "gene_feature_sha256": sha256_file(matrix_path),
        "derived_feature_path": str(derived_path),
        "derived_feature_sha256": sha256_file(derived_path),
        "derived_feature_names": derived_names,
        "legacy_score_path": str(legacy_score_path),
        "legacy_score_sha256": sha256_file(legacy_score_path),
        "legacy_score_primary_finite": bool(
            np.isfinite(legacy_score[cells["truth_class"].isin(["gdT_primary", "abT_primary"]).to_numpy(dtype=bool)]).all()
        ),
        "exclusion_audit": exclusion_audit,
    }
    write_json(summary, summary_path)
    return summary


def extract_legacy_score(cells: pd.DataFrame, config: Mapping[str, Any]) -> np.ndarray:
    """Read the frozen Phase-4 TRD-minus-TRAB comparator without recomputation."""
    inputs = pd.read_csv(resolve(config["preflight"]["input_manifest"]))
    output = np.full(cells.shape[0], np.nan, dtype=np.float32)
    for input_id in sorted(cells["expression_input_id"].unique()):
        if input_id == "gse144469":
            continue
        input_row = inputs[inputs["input_id"].eq(input_id)]
        if input_row.shape[0] != 1:
            raise RuntimeError(f"Legacy comparator input is not unique: {input_id}")
        path = Path(input_row.iloc[0]["path"])
        subset = cells.index[cells["expression_input_id"].eq(input_id)].to_numpy(dtype=np.int64)
        expression_rows = cells.loc[subset, "expression_row"].to_numpy(dtype=np.int64)
        with h5py.File(path, "r") as handle:
            if "phase4_trd_minus_trab" not in handle["obs"]:
                logging.warning("No legacy score in sensitivity input %s", input_id)
                continue
            values = np.asarray(read_axis_column(handle["obs"], "phase4_trd_minus_trab"), dtype=np.float32)
        output[subset] = values[expression_rows]

    gse_rows = cells.index[cells["expression_input_id"].eq("gse144469")].to_numpy(dtype=np.int64)
    legacy_input = inputs[inputs["input_id"].eq("legacy_metadata_reference")]
    if legacy_input.shape[0] != 1:
        raise RuntimeError("Legacy metadata reference is not uniquely registered")
    with h5py.File(Path(legacy_input.iloc[0]["path"]), "r") as handle:
        source = np.asarray(read_axis_column(handle["obs"], "source_gse_id"), dtype=object)
        source_rows = np.flatnonzero(source == "GSE144469")
        ids = np.asarray(read_axis_column(handle["obs"], axis_index_key(handle["obs"])), dtype=object)[source_rows]
        score = np.asarray(read_axis_column(handle["obs"], "phase4_trd_minus_trab"), dtype=np.float32)[source_rows]
    if pd.Index(ids).duplicated().any():
        raise RuntimeError("Legacy GSE144469 score IDs are duplicated")
    mapping = pd.Series(score, index=pd.Index(ids))
    joined = mapping.reindex(cells.loc[gse_rows, "cell_id"].astype(str)).to_numpy(dtype=np.float32)
    if not np.isfinite(joined).all():
        raise RuntimeError("Legacy score did not join every GSE144469 manifest cell")
    output[gse_rows] = joined
    primary = cells["truth_class"].isin(["gdT_primary", "abT_primary"]).to_numpy(dtype=bool)
    if not np.isfinite(output[primary]).all():
        raise RuntimeError("Legacy comparator is missing primary truth scores")
    return output


@dataclass
class CrossfitResult:
    candidate_id: str
    family: str
    parameters: dict[str, Any]
    train_rows: np.ndarray
    calibrated_oof: np.ndarray
    calibrated_control: np.ndarray
    calibrator: PlattCalibrator
    fold_models: list[FittedBinaryEstimator]
    retained_feature_counts: list[int]
    retained_gene_counts: list[int]


def inner_assignments(
    cells: pd.DataFrame,
    split_manifest: pd.DataFrame,
    outer_fold_id: str,
    stage: str,
    train_rows: np.ndarray,
) -> np.ndarray:
    frozen = split_manifest[
        split_manifest["outer_fold_id"].eq(outer_fold_id) & split_manifest["stage"].eq(stage)
    ]
    mapping = frozen.set_index("group_key")["inner_fold"]
    assignments = cells.loc[train_rows, "group_key"].map(mapping)
    if assignments.isna().any():
        missing = cells.loc[train_rows[assignments.isna().to_numpy()], "group_key"].drop_duplicates().head(10).tolist()
        raise RuntimeError(f"Frozen {stage} assignments missing groups: {missing}")
    result = assignments.to_numpy(dtype=np.int8)
    if set(np.unique(result)) != {0, 1, 2}:
        raise RuntimeError(f"Frozen {stage} assignments do not contain all three folds")
    return result


def crossfit_candidate(
    *,
    matrix: np.ndarray,
    feature_names: Sequence[str],
    gene_feature_count: int,
    train_rows: np.ndarray,
    assignments: np.ndarray,
    labels: np.ndarray,
    sources: np.ndarray,
    reliability: np.ndarray,
    control_rows: np.ndarray,
    family: str,
    parameters: Mapping[str, Any],
    config: Mapping[str, Any],
    seed_tokens: Sequence[object],
) -> CrossfitResult:
    train_rows = np.asarray(train_rows, dtype=np.int64)
    assignments = np.asarray(assignments, dtype=np.int8)
    y = np.asarray(labels, dtype=np.int8)
    source = np.asarray(sources, dtype=object)
    rel = np.asarray(reliability, dtype=np.float32)
    if not (train_rows.shape == assignments.shape == y.shape == source.shape == rel.shape):
        raise ValueError("Cross-fitting arrays have inconsistent shapes")
    model_id = candidate_id(family, parameters)
    raw_oof = np.full(train_rows.size, np.nan, dtype=np.float64)
    raw_control = np.zeros(control_rows.size, dtype=np.float64)
    fold_models: list[FittedBinaryEstimator] = []
    retained_features: list[int] = []
    retained_genes: list[int] = []
    for fold in (0, 1, 2):
        training = assignments != fold
        validation = assignments == fold
        fit_rows = train_rows[training]
        validation_rows = train_rows[validation]
        weights = balanced_training_weights(y[training], source[training], rel[training])
        model = fit_binary_estimator(
            np.asarray(matrix[fit_rows]),
            y[training],
            weights,
            feature_names,
            gene_feature_count,
            family,
            parameters,
            float(config["feature_policy"]["minimum_detection_fraction"]),
            int(config["feature_policy"]["maximum_retained_genes"]),
            stable_seed(int(config["random_seed"]), *seed_tokens, model_id, fold),
        )
        raw_oof[validation] = model.predict_probability(np.asarray(matrix[validation_rows]))
        if control_rows.size:
            raw_control += model.predict_probability(np.asarray(matrix[control_rows])) / 3.0
        fold_models.append(model)
        retained_features.append(len(model.selected_feature_names))
        retained_genes.append(int((model.selected_columns < gene_feature_count).sum()))
    if not np.isfinite(raw_oof).all():
        raise RuntimeError(f"Incomplete OOF predictions for {model_id}")
    outer_weights = balanced_training_weights(y, source, rel)
    calibrator = fit_platt_calibrator(
        raw_oof,
        y,
        outer_weights,
        stable_seed(int(config["random_seed"]), *seed_tokens, model_id, "calibration"),
    )
    return CrossfitResult(
        candidate_id=model_id,
        family=family,
        parameters=dict(parameters),
        train_rows=train_rows,
        calibrated_oof=calibrator.predict(raw_oof),
        calibrated_control=calibrator.predict(raw_control) if control_rows.size else np.asarray([], dtype=np.float64),
        calibrator=calibrator,
        fold_models=fold_models,
        retained_feature_counts=retained_features,
        retained_gene_counts=retained_genes,
    )


def fit_outer_model(
    *,
    matrix: np.ndarray,
    feature_names: Sequence[str],
    gene_feature_count: int,
    train_rows: np.ndarray,
    labels: np.ndarray,
    sources: np.ndarray,
    reliability: np.ndarray,
    family: str,
    parameters: Mapping[str, Any],
    config: Mapping[str, Any],
    seed_tokens: Sequence[object],
) -> FittedBinaryEstimator:
    weights = balanced_training_weights(labels, sources, reliability)
    return fit_binary_estimator(
        np.asarray(matrix[train_rows]),
        labels,
        weights,
        feature_names,
        gene_feature_count,
        family,
        parameters,
        float(config["feature_policy"]["minimum_detection_fraction"]),
        int(config["feature_policy"]["maximum_retained_genes"]),
        stable_seed(int(config["random_seed"]), *seed_tokens, "outer_fit"),
    )


def add_fixed_parameters(grid: list[dict[str, Any]], fixed: Mapping[str, Any]) -> list[dict[str, Any]]:
    return [{**row, **fixed} for row in grid]


def stage1_candidates(config: Mapping[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    spec = config["models"]["stage1_elastic_net"]
    grid = parameter_grid({"C": spec["C"], "l1_ratio": spec["l1_ratio"]})
    fixed = {"max_iter": spec["max_iter"], "tolerance": spec["tolerance"]}
    return [("elastic_net", parameters) for parameters in add_fixed_parameters(grid, fixed)]


def stage2_candidates(config: Mapping[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    elastic = config["models"]["stage2_elastic_net"]
    elastic_grid = add_fixed_parameters(
        parameter_grid({"C": elastic["C"], "l1_ratio": elastic["l1_ratio"]}),
        {"max_iter": elastic["max_iter"], "tolerance": elastic["tolerance"]},
    )
    hgb = config["models"]["stage2_hist_gradient_boosting"]
    hgb_grid = add_fixed_parameters(
        parameter_grid(
            {
                "learning_rate": hgb["learning_rate"],
                "max_leaf_nodes": hgb["max_leaf_nodes"],
                "min_samples_leaf": hgb["min_samples_leaf"],
                "l2_regularization": hgb["l2_regularization"],
            }
        ),
        {"max_iter": hgb["max_iter"], "early_stopping": hgb["early_stopping"]},
    )
    return (
        [("elastic_net", parameters) for parameters in elastic_grid],
        [("hist_gradient_boosting", parameters) for parameters in hgb_grid],
    )


def flatten_stage2_candidates(config: Mapping[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    elastic, hgb = stage2_candidates(config)
    return [*elastic, *hgb]


def select_best_stage1(rows: list[dict[str, Any]]) -> dict[str, Any]:
    eligible = [row for row in rows if row["decision"].passed]
    if not eligible:
        raise RuntimeError("No Stage-1 candidate satisfies the frozen recall guardrails")
    return sorted(
        eligible,
        key=lambda row: (
            -float(row["decision"].objective_value),
            -float(row["decision"].macro_recall),
            float(row["mean_retained_genes"]),
            str(row["candidate_id"]),
        ),
    )[0]


def select_best_stage2(rows: list[dict[str, Any]]) -> dict[str, Any] | None:
    eligible = [row for row in rows if row["balanced"].passed]
    if not eligible:
        return None
    return sorted(
        eligible,
        key=lambda row: (
            -float(row["balanced"].objective_value),
            float(row["balanced"].strict_nk_fpr),
            float(row["balanced"].paired_abt_fpr),
            float(row["mean_retained_genes"]),
            str(row["candidate_id"]),
        ),
    )[0]


def decision_record(decision: ThresholdDecision) -> dict[str, Any]:
    return asdict(decision)


def candidate_record(
    result: CrossfitResult,
    decisions: Mapping[str, ThresholdDecision],
    heldout_source: str,
    stage: str,
) -> dict[str, Any]:
    output: dict[str, Any] = {
        "heldout_source": heldout_source,
        "stage": stage,
        "candidate_id": result.candidate_id,
        "family": result.family,
        "parameters_json": json.dumps(result.parameters, sort_keys=True),
        "mean_retained_features": float(np.mean(result.retained_feature_counts)),
        "min_retained_features": int(min(result.retained_feature_counts)),
        "max_retained_features": int(max(result.retained_feature_counts)),
        "mean_retained_genes": float(np.mean(result.retained_gene_counts)),
        "calibrator_intercept": result.calibrator.intercept,
        "calibrator_slope": result.calibrator.slope,
    }
    for mode, decision in decisions.items():
        for key, value in decision_record(decision).items():
            output[f"{mode}_{key}"] = value
    return output


def evaluate_stage1_outer(
    *,
    cells: pd.DataFrame,
    splits: pd.DataFrame,
    gene_matrix: np.ndarray,
    feature_names: list[str],
    heldout_source: str,
    outer_fold_id: str,
    config: Mapping[str, Any],
) -> tuple[dict[str, Any], CrossfitResult, FittedBinaryEstimator, np.ndarray, pd.DataFrame]:
    stage1_lookup = [feature_names.index(gene) for gene in config["feature_policy"]["stage1_genes"]]
    matrix = np.asarray(gene_matrix[:, stage1_lookup], dtype=np.float32)
    names = [feature_names[index] for index in stage1_lookup]
    eligible = cells["stage1_role"].isin(["t_positive", "nk_negative"])
    train_rows = cells.index[eligible & cells["source_gse_id"].ne(heldout_source)].to_numpy(dtype=np.int64)
    assignments = inner_assignments(cells, splits, outer_fold_id, "stage1", train_rows)
    labels = cells.loc[train_rows, "stage1_role"].eq("t_positive").to_numpy(dtype=np.int8)
    sources = cells.loc[train_rows, "source_gse_id"].to_numpy(dtype=object)
    reliability = cells.loc[train_rows, "stage1_weight"].to_numpy(dtype=np.float32)
    primary_training_sources = [source for source in config_primary_sources(config) if source != heldout_source]
    guard = config["stage1_guardrails"]

    def fit_one(family: str, parameters: dict[str, Any]) -> tuple[dict[str, Any], CrossfitResult]:
        result = crossfit_candidate(
            matrix=matrix,
            feature_names=names,
            gene_feature_count=len(names),
            train_rows=train_rows,
            assignments=assignments,
            labels=labels,
            sources=sources,
            reliability=reliability,
            control_rows=np.asarray([], dtype=np.int64),
            family=family,
            parameters=parameters,
            config=config,
            seed_tokens=(outer_fold_id, "stage1"),
        )
        local = cells.loc[train_rows]
        decision = select_stage1_threshold(
            result.calibrated_oof,
            sources,
            local["truth_class"].eq("gdT_primary").to_numpy(dtype=bool),
            local["has_any_ab_tcr"].to_numpy(dtype=bool) & local["stage1_role"].eq("t_positive").to_numpy(dtype=bool),
            local["stage1_role"].eq("nk_negative").to_numpy(dtype=bool),
            primary_training_sources,
            float(guard["gdt_recall_per_source"]),
            float(guard["abt_recall_per_source"]),
        )
        row = {
            "candidate_id": result.candidate_id,
            "family": family,
            "parameters": parameters,
            "decision": decision,
            "mean_retained_genes": float(np.mean(result.retained_gene_counts)),
        }
        return row, result

    candidates = stage1_candidates(config)
    jobs = int(config.get("_candidate_jobs", 1))
    with parallel_config(backend="loky", inner_max_num_threads=1):
        fitted = Parallel(n_jobs=jobs)(delayed(fit_one)(family, parameters) for family, parameters in candidates)
    results = [row for row, _ in fitted]
    retained_results = {result.candidate_id: result for _, result in fitted}
    best = select_best_stage1(results)
    selected = retained_results[best["candidate_id"]]
    final_model = fit_outer_model(
        matrix=matrix,
        feature_names=names,
        gene_feature_count=len(names),
        train_rows=train_rows,
        labels=labels,
        sources=sources,
        reliability=reliability,
        family=selected.family,
        parameters=selected.parameters,
        config=config,
        seed_tokens=(outer_fold_id, "stage1", selected.candidate_id),
    )
    final_probability = selected.calibrator.predict(final_model.predict_probability(matrix))
    crossfit_probability = final_probability.copy()
    crossfit_probability[train_rows] = selected.calibrated_oof
    table = pd.DataFrame(
        [
            candidate_record(retained_results[row["candidate_id"]], {"stage1": row["decision"]}, heldout_source, "stage1")
            for row in results
        ]
    )
    return best, selected, final_model, crossfit_probability, table


def config_primary_sources(config: Mapping[str, Any]) -> list[str]:
    preflight_config = load_json(PREFLIGHT_CONFIG)
    return [str(value) for value in preflight_config["primary_sources"]]


def evaluate_stage2_candidate_grid(
    *,
    cells: pd.DataFrame,
    matrix: np.ndarray,
    feature_names: list[str],
    gene_feature_count: int,
    train_rows: np.ndarray,
    assignments: np.ndarray,
    control_rows: np.ndarray,
    stage1_probability: np.ndarray,
    stage1_threshold: float,
    heldout_source: str,
    outer_fold_id: str,
    config: Mapping[str, Any],
    candidates: Sequence[tuple[str, dict[str, Any]]],
    stage_name: str,
) -> tuple[list[dict[str, Any]], dict[str, CrossfitResult]]:
    labels = cells.loc[train_rows, "stage2_role"].eq("positive").to_numpy(dtype=np.int8)
    sources = cells.loc[train_rows, "source_gse_id"].to_numpy(dtype=object)
    reliability = cells.loc[train_rows, "truth_reliability"].to_numpy(dtype=np.float32)
    primary_eval_mask = cells.loc[train_rows, "truth_class"].isin(["gdT_primary", "abT_primary"]).to_numpy(dtype=bool)
    eval_rows = train_rows[primary_eval_mask]
    control_exclusion = cells.loc[control_rows, "exclusion_union"].to_numpy(dtype=bool)
    def fit_one(family: str, parameters: dict[str, Any]) -> tuple[dict[str, Any], CrossfitResult]:
        result = crossfit_candidate(
            matrix=matrix,
            feature_names=feature_names,
            gene_feature_count=gene_feature_count,
            train_rows=train_rows,
            assignments=assignments,
            labels=labels,
            sources=sources,
            reliability=reliability,
            control_rows=control_rows,
            family=family,
            parameters=parameters,
            config=config,
            seed_tokens=(outer_fold_id, stage_name),
        )
        eval_probability = result.calibrated_oof[primary_eval_mask]
        eval_frame = cells.loc[eval_rows]
        decisions: dict[str, ThresholdDecision] = {}
        for mode, mode_spec in config["operating_modes"].items():
            decisions[mode] = select_stage2_threshold(
                eval_probability,
                eval_frame["truth_class"].eq("gdT_primary").to_numpy(dtype=np.int8),
                eval_frame["source_gse_id"].to_numpy(dtype=object),
                eval_frame["truth_class"].eq("abT_primary").to_numpy(dtype=bool),
                result.calibrated_control,
                stage1_probability[eval_rows],
                stage1_probability[control_rows],
                stage1_threshold,
                eval_frame["exclusion_union"].to_numpy(dtype=bool),
                control_exclusion,
                mode,
                mode_spec,
            )
        record = {
            "candidate_id": result.candidate_id,
            "family": family,
            "parameters": parameters,
            "balanced": decisions["balanced"],
            "high_purity": decisions["high_purity"],
            "mean_retained_genes": float(np.mean(result.retained_gene_counts)),
        }
        return record, result

    jobs = int(config.get("_candidate_jobs", 1))
    with parallel_config(backend="loky", inner_max_num_threads=1):
        fitted = Parallel(n_jobs=jobs)(delayed(fit_one)(family, parameters) for family, parameters in candidates)
    records = [record for record, _ in fitted]
    results = {result.candidate_id: result for _, result in fitted}
    return records, results


def stage2_candidate_table(
    records: list[dict[str, Any]],
    results: Mapping[str, CrossfitResult],
    heldout_source: str,
    stage_name: str,
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            candidate_record(
                results[row["candidate_id"]],
                {"balanced": row["balanced"], "high_purity": row["high_purity"]},
                heldout_source,
                stage_name,
            )
            for row in records
        ]
    )


def evaluate_outer_predictions(
    *,
    cells: pd.DataFrame,
    heldout_source: str,
    model_id: str,
    mode: str,
    stage1_probability: np.ndarray,
    stage2_probability: np.ndarray,
    stage1_threshold: float,
    stage2_threshold: float,
) -> tuple[dict[str, Any], pd.DataFrame]:
    primary_mask = cells["source_gse_id"].eq(heldout_source) & cells["truth_class"].isin(["gdT_primary", "abT_primary"])
    rows = cells.index[primary_mask].to_numpy(dtype=np.int64)
    prediction = apply_two_stage_call(
        stage1_probability[rows],
        stage2_probability[rows],
        stage1_threshold,
        stage2_threshold,
        cells.loc[rows, "exclusion_union"].to_numpy(dtype=bool),
    )
    labels = cells.loc[rows, "truth_class"].eq("gdT_primary").to_numpy(dtype=np.int8)
    metrics = binary_metrics(labels, prediction, stage2_probability[rows])
    metrics.update(
        {
            "heldout_source": heldout_source,
            "model_id": model_id,
            "mode": mode,
            "evaluation": "primary_gold",
            "stage1_threshold": stage1_threshold,
            "stage2_threshold": stage2_threshold,
        }
    )
    predictions = pd.DataFrame(
        {
            "cell_id": cells.loc[rows, "cell_id"].to_numpy(dtype=object),
            "source_gse_id": heldout_source,
            "group_key": cells.loc[rows, "group_key"].to_numpy(dtype=object),
            "truth_class": cells.loc[rows, "truth_class"].to_numpy(dtype=object),
            "stage1_probability": stage1_probability[rows],
            "stage2_probability": stage2_probability[rows],
            "exclusion_union": cells.loc[rows, "exclusion_union"].to_numpy(dtype=bool),
            "prediction": prediction,
            "model_id": model_id,
            "mode": mode,
        }
    )
    return metrics, predictions


def heldout_nk_control_metrics(
    *,
    cells: pd.DataFrame,
    heldout_source: str,
    model_id: str,
    mode: str,
    stage1_probability: np.ndarray,
    stage2_probability: np.ndarray,
    stage1_threshold: float,
    stage2_threshold: float,
) -> dict[str, Any]:
    mask = cells["source_gse_id"].eq(heldout_source) & cells["stage1_role"].eq("nk_negative")
    rows = cells.index[mask].to_numpy(dtype=np.int64)
    prediction = apply_two_stage_call(
        stage1_probability[rows],
        stage2_probability[rows],
        stage1_threshold,
        stage2_threshold,
        cells.loc[rows, "exclusion_union"].to_numpy(dtype=bool),
    )
    total = int(rows.size)
    predicted = int(prediction.sum())
    return {
        "heldout_source": heldout_source,
        "model_id": model_id,
        "mode": mode,
        "evaluation": "strict_nk",
        "n_cells": total,
        "predicted_positive": predicted,
        "fpr": float(predicted / total) if total else math.nan,
        "stage1_threshold": stage1_threshold,
        "stage2_threshold": stage2_threshold,
    }


def heldout_stage1_guardrails(
    cells: pd.DataFrame,
    heldout_source: str,
    stage1_probability: np.ndarray,
    threshold: float,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    source = cells["source_gse_id"].eq(heldout_source)
    gdt = source & cells["truth_class"].eq("gdT_primary")
    abt = source & cells["stage1_role"].eq("t_positive") & cells["has_any_ab_tcr"].astype(bool)
    nk = source & cells["stage1_role"].eq("nk_negative")

    def rate(mask: pd.Series) -> float:
        rows = cells.index[mask].to_numpy(dtype=np.int64)
        return float((stage1_probability[rows] >= threshold).mean()) if rows.size else math.nan

    gdt_recall = rate(gdt)
    abt_recall = rate(abt)
    nk_call_rate = rate(nk)
    guard = config["stage1_guardrails"]
    return {
        "heldout_source": heldout_source,
        "stage1_threshold": threshold,
        "n_gdt": int(gdt.sum()),
        "gdt_recall": gdt_recall,
        "required_gdt_recall": float(guard["gdt_recall_per_source"]),
        "n_abt": int(abt.sum()),
        "abt_recall": abt_recall,
        "required_abt_recall": float(guard["abt_recall_per_source"]),
        "n_strict_nk": int(nk.sum()),
        "strict_nk_call_rate": nk_call_rate,
        "pass": bool(
            np.isfinite(gdt_recall)
            and gdt_recall >= float(guard["gdt_recall_per_source"])
            and np.isfinite(abt_recall)
            and abt_recall >= float(guard["abt_recall_per_source"])
        ),
    }


def evaluate_single_stage_comparator(
    *,
    cells: pd.DataFrame,
    splits: pd.DataFrame,
    matrix: np.ndarray,
    feature_names: list[str],
    gene_feature_count: int,
    heldout_source: str,
    outer_fold_id: str,
    config: Mapping[str, Any],
    model_id: str,
    candidates: Sequence[tuple[str, dict[str, Any]]],
) -> dict[str, Any]:
    eligible = cells["stage2_role"].isin(["positive", "negative"])
    train_rows = cells.index[eligible & cells["source_gse_id"].ne(heldout_source)].to_numpy(dtype=np.int64)
    assignments = inner_assignments(cells, splits, outer_fold_id, "stage2", train_rows)
    control_rows = cells.index[
        cells["stage1_role"].eq("nk_negative") & cells["source_gse_id"].ne(heldout_source)
    ].to_numpy(dtype=np.int64)
    ones = np.ones(cells.shape[0], dtype=np.float64)
    records, results = evaluate_stage2_candidate_grid(
        cells=cells,
        matrix=matrix,
        feature_names=feature_names,
        gene_feature_count=gene_feature_count,
        train_rows=train_rows,
        assignments=assignments,
        control_rows=control_rows,
        stage1_probability=ones,
        stage1_threshold=0.0,
        heldout_source=heldout_source,
        outer_fold_id=outer_fold_id,
        config=config,
        candidates=candidates,
        stage_name=model_id,
    )
    selected_record = select_best_stage2(records)
    output: dict[str, Any] = {
        "candidate_table": stage2_candidate_table(records, results, heldout_source, model_id),
        "selected": None,
        "metrics": [],
        "predictions": [],
        "control_metrics": [],
    }
    if selected_record is None:
        return output
    selected = results[selected_record["candidate_id"]]
    labels = cells.loc[train_rows, "stage2_role"].eq("positive").to_numpy(dtype=np.int8)
    sources = cells.loc[train_rows, "source_gse_id"].to_numpy(dtype=object)
    reliability = cells.loc[train_rows, "truth_reliability"].to_numpy(dtype=np.float32)
    final_model = fit_outer_model(
        matrix=matrix,
        feature_names=feature_names,
        gene_feature_count=gene_feature_count,
        train_rows=train_rows,
        labels=labels,
        sources=sources,
        reliability=reliability,
        family=selected.family,
        parameters=selected.parameters,
        config=config,
        seed_tokens=(outer_fold_id, model_id, selected.candidate_id),
    )
    probability = selected.calibrator.predict(final_model.predict_probability(matrix))
    output["selected"] = {
        "candidate_id": selected.candidate_id,
        "family": selected.family,
        "parameters": selected.parameters,
        "balanced": decision_record(selected_record["balanced"]),
        "high_purity": decision_record(selected_record["high_purity"]),
        "model": final_model,
        "calibrator": selected.calibrator,
    }
    for mode in ("balanced", "high_purity"):
        decision = selected_record[mode]
        if not decision.passed:
            continue
        metrics, predictions = evaluate_outer_predictions(
            cells=cells,
            heldout_source=heldout_source,
            model_id=model_id,
            mode=mode,
            stage1_probability=ones,
            stage2_probability=probability,
            stage1_threshold=0.0,
            stage2_threshold=float(decision.threshold),
        )
        output["metrics"].append(metrics)
        output["predictions"].append(predictions)
        output["control_metrics"].append(
            heldout_nk_control_metrics(
                cells=cells,
                heldout_source=heldout_source,
                model_id=model_id,
                mode=mode,
                stage1_probability=ones,
                stage2_probability=probability,
                stage1_threshold=0.0,
                stage2_threshold=float(decision.threshold),
            )
        )
    return output


def evaluate_legacy_comparator(
    *,
    cells: pd.DataFrame,
    legacy_score: np.ndarray,
    heldout_source: str,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    probability = expit(np.asarray(legacy_score, dtype=np.float64))
    train_primary = cells["source_gse_id"].ne(heldout_source) & cells["truth_class"].isin(["gdT_primary", "abT_primary"])
    train_rows = cells.index[train_primary].to_numpy(dtype=np.int64)
    control_rows = cells.index[
        cells["stage1_role"].eq("nk_negative") & cells["source_gse_id"].ne(heldout_source)
    ].to_numpy(dtype=np.int64)
    if not np.isfinite(probability[np.concatenate([train_rows, control_rows])]).all():
        raise RuntimeError("Legacy comparator lacks a training or strict-NK score")
    train = cells.loc[train_rows]
    decisions: dict[str, ThresholdDecision] = {}
    ones_primary = np.ones(train_rows.size, dtype=np.float64)
    ones_nk = np.ones(control_rows.size, dtype=np.float64)
    for mode, mode_spec in config["operating_modes"].items():
        decisions[mode] = select_stage2_threshold(
            probability[train_rows],
            train["truth_class"].eq("gdT_primary").to_numpy(dtype=np.int8),
            train["source_gse_id"].to_numpy(dtype=object),
            train["truth_class"].eq("abT_primary").to_numpy(dtype=bool),
            probability[control_rows],
            ones_primary,
            ones_nk,
            0.0,
            train["exclusion_union"].to_numpy(dtype=bool),
            cells.loc[control_rows, "exclusion_union"].to_numpy(dtype=bool),
            mode,
            mode_spec,
        )
    output: dict[str, Any] = {
        "selected": {
            mode: {
                **decision_record(decision),
                "raw_trd_minus_trab_threshold": float(np.log(decision.threshold / (1.0 - decision.threshold)))
                if decision.passed and 0 < decision.threshold < 1
                else math.nan,
            }
            for mode, decision in decisions.items()
        },
        "metrics": [],
        "predictions": [],
        "control_metrics": [],
    }
    ones = np.ones(cells.shape[0], dtype=np.float64)
    for mode, decision in decisions.items():
        if not decision.passed:
            continue
        metrics, predictions = evaluate_outer_predictions(
            cells=cells,
            heldout_source=heldout_source,
            model_id="legacy_trd_minus_trab",
            mode=mode,
            stage1_probability=ones,
            stage2_probability=probability,
            stage1_threshold=0.0,
            stage2_threshold=float(decision.threshold),
        )
        metrics["raw_trd_minus_trab_threshold"] = output["selected"][mode]["raw_trd_minus_trab_threshold"]
        output["metrics"].append(metrics)
        output["predictions"].append(predictions)
        output["control_metrics"].append(
            heldout_nk_control_metrics(
                cells=cells,
                heldout_source=heldout_source,
                model_id="legacy_trd_minus_trab",
                mode=mode,
                stage1_probability=ones,
                stage2_probability=probability,
                stage1_threshold=0.0,
                stage2_threshold=float(decision.threshold),
            )
        )
    return output


def selected_model_feature_importance(
    *,
    model: FittedBinaryEstimator,
    matrix: np.ndarray,
    labels: np.ndarray,
    feature_names: Sequence[str],
    heldout_source: str,
    config: Mapping[str, Any],
) -> pd.DataFrame:
    selected_names = list(model.selected_feature_names)
    if model.family == "elastic_net":
        coefficients = np.asarray(model.estimator.coef_[0], dtype=np.float64)
        importance = np.abs(coefficients)
        direction = np.sign(coefficients).astype(np.int8)
        method = "standardized_elastic_net_coefficient"
    elif model.family == "hist_gradient_boosting":
        settings = config["feature_stability"]
        rng = np.random.default_rng(
            stable_seed(int(config["random_seed"]), heldout_source, "permutation_importance")
        )
        n_max = int(settings["permutation_max_cells"])
        if matrix.shape[0] > n_max:
            sampled = np.sort(rng.choice(matrix.shape[0], size=n_max, replace=False))
        else:
            sampled = np.arange(matrix.shape[0], dtype=np.int64)
        x = np.asarray(matrix[sampled], dtype=np.float32)
        y = np.asarray(labels, dtype=np.int8)[sampled]
        baseline_probability = model.predict_probability(x)
        baseline = average_precision_score(y, baseline_probability)
        importance = np.zeros(len(selected_names), dtype=np.float64)
        direction = np.zeros(len(selected_names), dtype=np.int8)
        repeats = int(settings["permutation_repeats"])
        selected_global_columns = model.selected_columns
        for local_index, global_column in enumerate(selected_global_columns):
            losses: list[float] = []
            for _ in range(repeats):
                permuted = x.copy()
                permuted[:, global_column] = permuted[rng.permutation(permuted.shape[0]), global_column]
                losses.append(baseline - average_precision_score(y, model.predict_probability(permuted)))
            importance[local_index] = max(0.0, float(np.mean(losses)))
            values = x[:, global_column]
            if np.ptp(values) > 0:
                correlation = np.corrcoef(values, baseline_probability)[0, 1]
                direction[local_index] = 0 if not np.isfinite(correlation) else int(np.sign(correlation))
        method = "heldout_permutation_pr_auc_drop"
    else:
        raise ValueError(f"Unsupported selected V4 family for importance: {model.family}")
    total = float(importance.sum())
    normalized = importance / total if total > 0 else np.zeros_like(importance)
    table = pd.DataFrame(
        {
            "heldout_source": heldout_source,
            "family": model.family,
            "method": method,
            "feature": selected_names,
            "importance": importance,
            "normalized_importance": normalized,
            "direction": direction,
        }
    )
    table["rank"] = table["importance"].rank(method="first", ascending=False).astype(int)
    return table.sort_values("rank", kind="mergesort").reset_index(drop=True)


def run_outer_fold(
    *,
    cells: pd.DataFrame,
    splits: pd.DataFrame,
    gene_matrix: np.ndarray,
    derived_matrix: np.ndarray,
    legacy_score: np.ndarray,
    feature_names: list[str],
    heldout_source: str,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    outer_index = config_primary_sources(config).index(heldout_source)
    outer_fold_id = f"outer_{outer_index}_{heldout_source}"
    logging.info("Starting outer fold %s", outer_fold_id)
    stage1_best, stage1_crossfit, stage1_final, p1, stage1_table = evaluate_stage1_outer(
        cells=cells,
        splits=splits,
        gene_matrix=gene_matrix,
        feature_names=feature_names,
        heldout_source=heldout_source,
        outer_fold_id=outer_fold_id,
        config=config,
    )
    stage1_threshold = float(stage1_best["decision"].threshold)
    stage2_eligible = cells["stage2_role"].isin(["positive", "negative"])
    train_rows = cells.index[stage2_eligible & cells["source_gse_id"].ne(heldout_source)].to_numpy(dtype=np.int64)
    assignments = inner_assignments(cells, splits, outer_fold_id, "stage2", train_rows)
    control_rows = cells.index[
        cells["stage1_role"].eq("nk_negative") & cells["source_gse_id"].ne(heldout_source)
    ].to_numpy(dtype=np.int64)
    stage2_matrix = np.column_stack([gene_matrix, derived_matrix, p1]).astype(np.float32, copy=False)
    stage2_names = [*feature_names, *DERIVED_FEATURE_ORDER, "stage1_probability"]
    records, results = evaluate_stage2_candidate_grid(
        cells=cells,
        matrix=stage2_matrix,
        feature_names=stage2_names,
        gene_feature_count=len(feature_names),
        train_rows=train_rows,
        assignments=assignments,
        control_rows=control_rows,
        stage1_probability=p1,
        stage1_threshold=stage1_threshold,
        heldout_source=heldout_source,
        outer_fold_id=outer_fold_id,
        config=config,
        candidates=flatten_stage2_candidates(config),
        stage_name="stage2_v4",
    )
    selected_record = select_best_stage2(records)
    result: dict[str, Any] = {
        "outer_fold_id": outer_fold_id,
        "heldout_source": heldout_source,
        "stage1_candidate_table": stage1_table,
        "stage2_candidate_table": stage2_candidate_table(records, results, heldout_source, "stage2_v4"),
        "stage1_selected": {
            "candidate_id": stage1_crossfit.candidate_id,
            "parameters": stage1_crossfit.parameters,
            "threshold": stage1_threshold,
        },
        "stage1_final_model": stage1_final,
        "stage1_calibrator": stage1_crossfit.calibrator,
        "stage2_selected": None,
        "comparator_candidate_tables": {},
        "comparator_selected": {},
        "metrics": [],
        "predictions": [],
        "control_metrics": [],
        "feature_importance": pd.DataFrame(),
        "stage1_heldout_guardrails": heldout_stage1_guardrails(
            cells,
            heldout_source,
            p1,
            stage1_threshold,
            config,
        ),
    }
    if selected_record is None:
        logging.error("No V4 Stage-2 candidate passed balanced guardrails in %s", outer_fold_id)
    else:
        selected = results[selected_record["candidate_id"]]
        labels = cells.loc[train_rows, "stage2_role"].eq("positive").to_numpy(dtype=np.int8)
        sources = cells.loc[train_rows, "source_gse_id"].to_numpy(dtype=object)
        reliability = cells.loc[train_rows, "truth_reliability"].to_numpy(dtype=np.float32)
        final_model = fit_outer_model(
            matrix=stage2_matrix,
            feature_names=stage2_names,
            gene_feature_count=len(feature_names),
            train_rows=train_rows,
            labels=labels,
            sources=sources,
            reliability=reliability,
            family=selected.family,
            parameters=selected.parameters,
            config=config,
            seed_tokens=(outer_fold_id, "stage2_v4", selected.candidate_id),
        )
        p2 = selected.calibrator.predict(final_model.predict_probability(stage2_matrix))
        result["stage2_selected"] = {
            "candidate_id": selected.candidate_id,
            "family": selected.family,
            "parameters": selected.parameters,
            "balanced": decision_record(selected_record["balanced"]),
            "high_purity": decision_record(selected_record["high_purity"]),
        }
        result["stage2_final_model"] = final_model
        result["stage2_calibrator"] = selected.calibrator
        heldout_primary_rows = cells.index[
            cells["source_gse_id"].eq(heldout_source)
            & cells["truth_class"].isin(["gdT_primary", "abT_primary"])
        ].to_numpy(dtype=np.int64)
        result["feature_importance"] = selected_model_feature_importance(
            model=final_model,
            matrix=np.asarray(stage2_matrix[heldout_primary_rows]),
            labels=cells.loc[heldout_primary_rows, "truth_class"].eq("gdT_primary").to_numpy(dtype=np.int8),
            feature_names=stage2_names,
            heldout_source=heldout_source,
            config=config,
        )
        for mode in ("balanced", "high_purity"):
            decision = selected_record[mode]
            if not decision.passed:
                continue
            metrics, predictions = evaluate_outer_predictions(
                cells=cells,
                heldout_source=heldout_source,
                model_id="v4_nested_selected",
                mode=mode,
                stage1_probability=p1,
                stage2_probability=p2,
                stage1_threshold=stage1_threshold,
                stage2_threshold=float(decision.threshold),
            )
            result["metrics"].append(metrics)
            result["predictions"].append(predictions)
            result["control_metrics"].append(
                heldout_nk_control_metrics(
                    cells=cells,
                    heldout_source=heldout_source,
                    model_id="v4_nested_selected",
                    mode=mode,
                    stage1_probability=p1,
                    stage2_probability=p2,
                    stage1_threshold=stage1_threshold,
                    stage2_threshold=float(decision.threshold),
                )
            )

        # Mandatory individual-gene-only ablation: selected family and
        # hyperparameters, but no deterministic summaries or Stage-1 feature.
        ablation = evaluate_single_stage_comparator(
            cells=cells,
            splits=splits,
            matrix=np.asarray(gene_matrix),
            feature_names=feature_names,
            gene_feature_count=len(feature_names),
            heldout_source=heldout_source,
            outer_fold_id=outer_fold_id,
            config=config,
            model_id="v4_individual_gene_ablation",
            candidates=[(selected.family, selected.parameters)],
        )
        result["comparator_candidate_tables"]["v4_individual_gene_ablation"] = ablation["candidate_table"]
        result["comparator_selected"]["v4_individual_gene_ablation"] = ablation["selected"]
        result["metrics"].extend(ablation["metrics"])
        result["predictions"].extend(ablation["predictions"])
        result["control_metrics"].extend(ablation["control_metrics"])

    elastic_candidates, _ = stage2_candidates(config)
    compact_indices = [feature_names.index(gene) for gene in config["models"]["compact_logistic_genes"]]
    compact = evaluate_single_stage_comparator(
        cells=cells,
        splits=splits,
        matrix=np.asarray(gene_matrix[:, compact_indices]),
        feature_names=[feature_names[index] for index in compact_indices],
        gene_feature_count=len(compact_indices),
        heldout_source=heldout_source,
        outer_fold_id=outer_fold_id,
        config=config,
        model_id="compact_7gene_logistic",
        candidates=elastic_candidates,
    )
    result["comparator_candidate_tables"]["compact_7gene_logistic"] = compact["candidate_table"]
    result["comparator_selected"]["compact_7gene_logistic"] = compact["selected"]
    result["metrics"].extend(compact["metrics"])
    result["predictions"].extend(compact["predictions"])
    result["control_metrics"].extend(compact["control_metrics"])

    tcr_indices = [
        index
        for index, gene in enumerate(feature_names)
        if gene != "TRAT1" and gene.startswith(("TRA", "TRB", "TRG", "TRD"))
    ]
    if len(tcr_indices) != 153:
        raise RuntimeError(f"TCR comparator expected 153 genes, found {len(tcr_indices)}")
    tcr = evaluate_single_stage_comparator(
        cells=cells,
        splits=splits,
        matrix=np.asarray(gene_matrix[:, tcr_indices]),
        feature_names=[feature_names[index] for index in tcr_indices],
        gene_feature_count=len(tcr_indices),
        heldout_source=heldout_source,
        outer_fold_id=outer_fold_id,
        config=config,
        model_id="v2_like_tcr_logistic",
        candidates=elastic_candidates,
    )
    result["comparator_candidate_tables"]["v2_like_tcr_logistic"] = tcr["candidate_table"]
    result["comparator_selected"]["v2_like_tcr_logistic"] = tcr["selected"]
    result["metrics"].extend(tcr["metrics"])
    result["predictions"].extend(tcr["predictions"])
    result["control_metrics"].extend(tcr["control_metrics"])

    legacy = evaluate_legacy_comparator(
        cells=cells,
        legacy_score=legacy_score,
        heldout_source=heldout_source,
        config=config,
    )
    result["comparator_selected"]["legacy_trd_minus_trab"] = legacy["selected"]
    result["metrics"].extend(legacy["metrics"])
    result["predictions"].extend(legacy["predictions"])
    result["control_metrics"].extend(legacy["control_metrics"])
    logging.info("Completed outer fold %s", outer_fold_id)
    return result


def save_outer_result(result: Mapping[str, Any], config: Mapping[str, Any]) -> None:
    table_dir = resolve(config["outputs"]["table_dir"])
    cache_dir = resolve(config["outputs"]["cache_dir"])
    table_dir.mkdir(parents=True, exist_ok=True)
    fold = str(result["outer_fold_id"])
    result["stage1_candidate_table"].to_csv(table_dir / f"{fold}_stage1_candidates.csv", index=False)
    result["stage2_candidate_table"].to_csv(table_dir / f"{fold}_stage2_candidates.csv", index=False)
    for model_id, table in result["comparator_candidate_tables"].items():
        table.to_csv(table_dir / f"{fold}_{model_id}_candidates.csv", index=False)
    if result["metrics"]:
        pd.DataFrame(result["metrics"]).to_csv(table_dir / f"{fold}_metrics.csv", index=False)
    if result["control_metrics"]:
        pd.DataFrame(result["control_metrics"]).to_csv(table_dir / f"{fold}_negative_controls.csv", index=False)
    if not result["feature_importance"].empty:
        result["feature_importance"].to_csv(table_dir / f"{fold}_feature_importance.csv", index=False)
    pd.DataFrame([result["stage1_heldout_guardrails"]]).to_csv(
        table_dir / f"{fold}_stage1_heldout_guardrails.csv", index=False
    )
    if result["predictions"]:
        pd.concat(result["predictions"], ignore_index=True).to_csv(
            table_dir / f"{fold}_predictions.csv.gz", index=False, compression={"method": "gzip", "mtime": 0}
        )
    model_payload = {
        key: result[key]
        for key in (
            "outer_fold_id",
            "heldout_source",
            "stage1_selected",
            "stage1_final_model",
            "stage1_calibrator",
            "stage2_selected",
            "stage2_final_model",
            "stage2_calibrator",
            "comparator_selected",
        )
        if key in result
    }
    model_path = cache_dir / f"{fold}_selected_models.pkl"
    with model_path.open("wb") as handle:
        pickle.dump(model_payload, handle, protocol=pickle.HIGHEST_PROTOCOL)
    summary = {
        "outer_fold_id": fold,
        "heldout_source": result["heldout_source"],
        "stage1_selected": result["stage1_selected"],
        "stage2_selected": result["stage2_selected"],
        "comparator_selected": {
            model_id: {
                key: value
                for key, value in selected.items()
                if key not in {"model", "calibrator"}
            }
            if selected is not None
            else None
            for model_id, selected in result["comparator_selected"].items()
        },
        "model_path": str(model_path),
        "model_sha256": sha256_file(model_path),
        "metrics": result["metrics"],
        "control_metrics": result["control_metrics"],
        "stage1_heldout_guardrails": result["stage1_heldout_guardrails"],
    }
    write_json(summary, table_dir / f"{fold}_summary.json")


def fold_bootstrap_counts(result: Mapping[str, Any]) -> dict[str, pd.DataFrame]:
    balanced = {
        frame["model_id"].iloc[0]: frame
        for frame in result["predictions"]
        if not frame.empty and frame["mode"].iloc[0] == "balanced"
    }
    if "v4_nested_selected" not in balanced:
        return {}
    candidate = balanced["v4_nested_selected"].sort_values("cell_id", kind="mergesort").reset_index(drop=True)
    output: dict[str, pd.DataFrame] = {}
    for comparator_id in ("legacy_trd_minus_trab", "compact_7gene_logistic", "v2_like_tcr_logistic"):
        if comparator_id not in balanced:
            continue
        comparator = balanced[comparator_id].sort_values("cell_id", kind="mergesort").reset_index(drop=True)
        if not candidate["cell_id"].equals(comparator["cell_id"]):
            raise RuntimeError(f"Paired bootstrap cells disagree for {comparator_id}")
        labels = candidate["truth_class"].eq("gdT_primary").to_numpy(dtype=np.int8)
        output[comparator_id] = grouped_confusion_counts(
            labels,
            {
                "v4": candidate["prediction"].to_numpy(dtype=bool),
                "comparator": comparator["prediction"].to_numpy(dtype=bool),
            },
            candidate["source_gse_id"].to_numpy(dtype=object),
            candidate["group_key"].to_numpy(dtype=object),
        )
    return output


def aggregate_nested_results(
    metrics: list[dict[str, Any]],
    controls: list[dict[str, Any]],
    bootstrap_counts: Mapping[str, list[pd.DataFrame]],
    extension_controls: pd.DataFrame,
    feature_importance: pd.DataFrame,
    stage1_guardrails: pd.DataFrame,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    table_dir = resolve(config["outputs"]["table_dir"])
    log_dir = resolve(config["outputs"]["log_dir"])
    table_dir.mkdir(parents=True, exist_ok=True)
    metrics_frame = pd.DataFrame(metrics)
    controls_frame = pd.DataFrame(controls)
    metrics_frame.to_csv(table_dir / "outer_fold_metrics.csv", index=False)
    controls_frame.to_csv(table_dir / "outer_fold_negative_controls.csv", index=False)
    stage1_guardrails.to_csv(table_dir / "stage1_heldout_guardrails.csv", index=False)
    macro_columns = [
        "precision",
        "recall",
        "specificity",
        "f1",
        "f0_5",
        "balanced_accuracy",
        "mcc",
        "roc_auc",
        "pr_auc",
        "brier_score",
        "calibration_intercept",
        "calibration_slope",
        "ece",
    ]
    macro = (
        metrics_frame.groupby(["model_id", "mode"], observed=True)[macro_columns]
        .mean()
        .reset_index()
    )
    fold_counts = (
        metrics_frame.groupby(["model_id", "mode"], observed=True)["heldout_source"]
        .nunique()
        .rename("n_outer_sources")
        .reset_index()
    )
    macro = macro.merge(fold_counts, on=["model_id", "mode"], how="left")
    macro.to_csv(table_dir / "dataset_macro_metrics.csv", index=False)
    fair = {"legacy_trd_minus_trab", "compact_7gene_logistic", "v2_like_tcr_logistic"}
    balanced_macro = macro[macro["mode"].eq("balanced") & macro["model_id"].isin(fair)]
    complete_balanced = balanced_macro[balanced_macro["n_outer_sources"].eq(len(config_primary_sources(config)))]
    strongest = None if complete_balanced.empty else str(complete_balanced.sort_values("f1", ascending=False).iloc[0]["model_id"])
    bootstrap_rows: list[dict[str, Any]] = []
    bootstrap_values: dict[str, np.ndarray] = {}
    for comparator_id, frames in bootstrap_counts.items():
        if len(frames) != len(config_primary_sources(config)):
            continue
        counts = pd.concat(frames, ignore_index=True)
        values = paired_hierarchical_bootstrap_f1_difference(
            counts,
            "v4",
            "comparator",
            int(config["bootstrap_replicates"]),
            stable_seed(int(config["random_seed"]), "paired_bootstrap", comparator_id),
        )
        low, high = bootstrap_interval(values)
        bootstrap_values[comparator_id] = values
        bootstrap_rows.append(
            {
                "candidate": "v4_nested_selected",
                "comparator": comparator_id,
                "replicates": int(values.size),
                "mean_f1_difference": float(values.mean()),
                "median_f1_difference": float(np.median(values)),
                "ci95_low": low,
                "ci95_high": high,
                "fraction_difference_gt_0": float((values > 0).mean()),
            }
        )
    bootstrap_frame = pd.DataFrame(bootstrap_rows)
    bootstrap_frame.to_csv(table_dir / "paired_hierarchical_bootstrap_f1.csv", index=False)

    v4_balanced = metrics_frame[
        metrics_frame["model_id"].eq("v4_nested_selected") & metrics_frame["mode"].eq("balanced")
    ]
    v4_macro_row = macro[
        macro["model_id"].eq("v4_nested_selected") & macro["mode"].eq("balanced")
    ]
    strongest_row = complete_balanced[complete_balanced["model_id"].eq(strongest)] if strongest else pd.DataFrame()
    observed_difference = (
        float(v4_macro_row.iloc[0]["f1"] - strongest_row.iloc[0]["f1"])
        if not v4_macro_row.empty and not strongest_row.empty
        else math.nan
    )
    strongest_bootstrap = bootstrap_frame[bootstrap_frame["comparator"].eq(strongest)] if strongest else pd.DataFrame()
    ci_low = float(strongest_bootstrap.iloc[0]["ci95_low"]) if not strongest_bootstrap.empty else math.nan
    source_recall_pass = bool(
        v4_balanced["heldout_source"].nunique() == len(config_primary_sources(config))
        and (v4_balanced["recall"] >= 0.70).all()
    )
    paired_fpr = v4_balanced["fp"] / v4_balanced["n_negative"]
    paired_guard_pass = bool(
        v4_balanced["heldout_source"].nunique() == len(config_primary_sources(config))
        and (paired_fpr <= float(config["operating_modes"]["balanced"]["maximum_paired_abt_fpr"])).all()
    )
    v4_nk = controls_frame[
        controls_frame["model_id"].eq("v4_nested_selected") & controls_frame["mode"].eq("balanced")
    ]
    nk_guard_pass = bool(
        v4_nk["heldout_source"].nunique() == len(config_primary_sources(config))
        and (v4_nk["fpr"] <= float(config["operating_modes"]["balanced"]["maximum_strict_nk_fpr"])).all()
    )
    stability_table, stability_summary = aggregate_feature_stability(feature_importance, config)
    stability_table.to_csv(table_dir / "feature_stability.csv", index=False)
    guardrails = pd.DataFrame(
        [
            {
                "condition": "stage1_heldout_recall_guardrails",
                "pass": bool(
                    stage1_guardrails.shape[0] == len(config_primary_sources(config))
                    and stage1_guardrails["pass"].astype(bool).all()
                ),
                "observed": float(
                    min(stage1_guardrails["gdt_recall"].min(), stage1_guardrails["abt_recall"].min())
                )
                if not stage1_guardrails.empty
                else math.nan,
                "required": "gdT >=0.99 and productive-abT >=0.98 in every held-out source",
            },
            {
                "condition": "macro_f1_improvement_at_least_0p01",
                "pass": bool(np.isfinite(observed_difference) and observed_difference >= 0.01),
                "observed": observed_difference,
                "required": ">=0.01",
            },
            {
                "condition": "paired_bootstrap_ci_low_gt_0",
                "pass": bool(np.isfinite(ci_low) and ci_low > 0),
                "observed": ci_low,
                "required": ">0",
            },
            {
                "condition": "recall_each_primary_source_at_least_0p70",
                "pass": source_recall_pass,
                "observed": float(v4_balanced["recall"].min()) if not v4_balanced.empty else math.nan,
                "required": ">=0.70",
            },
            {
                "condition": "paired_abt_fpr_each_primary_source_at_most_0p002",
                "pass": paired_guard_pass,
                "observed": float(paired_fpr.max()) if not paired_fpr.empty else math.nan,
                "required": "<=0.002",
            },
            {
                "condition": "strict_nk_fpr_each_primary_source_at_most_0p01",
                "pass": nk_guard_pass,
                "observed": float(v4_nk["fpr"].max()) if not v4_nk.empty else math.nan,
                "required": "<=0.01",
            },
            {
                "condition": "extension_negative_control_guardrails",
                "pass": extension_guardrails_pass(extension_controls, config),
                "observed": float(extension_controls["fpr"].max()) if not extension_controls.empty else math.nan,
                "required": "paired-abT <=0.002 and strict-NK <=0.01 in every eligible cohort/outer model",
            },
            {
                "condition": "feature_stability",
                "pass": bool(stability_summary["pass"]),
                "observed": float(stability_summary["stable_fraction"]),
                "required": "stable fraction >=0.70 and maximum single-dataset share <=0.50",
            },
            {
                "condition": "release_semantic_checksum_consistency",
                "pass": False,
                "observed": math.nan,
                "required": "Step 3 only",
            },
        ]
    )
    guardrails.to_csv(table_dir / "promotion_guardrails_step2.csv", index=False)
    summary = {
        "strongest_fair_comparator": strongest,
        "observed_macro_f1_difference": observed_difference,
        "paired_bootstrap_ci95_low": ci_low,
        "all_current_step2_guardrails_pass": bool(guardrails.iloc[:-1]["pass"].all()),
        "promotion_authorized": False,
        "promotion_reason": "Step 2 cannot authorize promotion; extension, stability, release, and supervision gates remain.",
        "table_paths": {
            "outer_metrics": str(table_dir / "outer_fold_metrics.csv"),
            "outer_negative_controls": str(table_dir / "outer_fold_negative_controls.csv"),
            "stage1_heldout_guardrails": str(table_dir / "stage1_heldout_guardrails.csv"),
            "dataset_macro_metrics": str(table_dir / "dataset_macro_metrics.csv"),
            "bootstrap": str(table_dir / "paired_hierarchical_bootstrap_f1.csv"),
            "promotion_guardrails": str(table_dir / "promotion_guardrails_step2.csv"),
            "feature_stability": str(table_dir / "feature_stability.csv"),
        },
        "feature_stability": stability_summary,
    }
    write_json(summary, log_dir / "nested_evaluation_aggregate.json")
    return summary


def aggregate_feature_stability(
    feature_importance: pd.DataFrame,
    config: Mapping[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    settings = config["feature_stability"]
    if feature_importance.empty:
        return pd.DataFrame(), {
            "pass": False,
            "stable_fraction": math.nan,
            "maximum_single_dataset_share": math.nan,
            "reason": "No selected V4 feature-importance tables",
        }
    rows: list[dict[str, Any]] = []
    top_k = int(settings["top_k"])
    for feature, frame in feature_importance.groupby("feature", sort=True):
        total = float(frame["normalized_importance"].sum())
        max_share = float(frame["normalized_importance"].max() / total) if total > 0 else math.nan
        direction_counts = frame.loc[frame["direction"].ne(0), "direction"].value_counts()
        direction_agreement = int(direction_counts.max()) if not direction_counts.empty else 0
        top_fold_count = int((frame["rank"] <= top_k).sum())
        rows.append(
            {
                "feature": feature,
                "aggregate_normalized_importance": total,
                "outer_fold_count": int(frame["heldout_source"].nunique()),
                "top_k_fold_count": top_fold_count,
                "direction_agreement_folds": direction_agreement,
                "directionally_stable": bool(top_fold_count >= 2 and direction_agreement >= 2),
                "maximum_single_dataset_share": max_share,
            }
        )
    table = pd.DataFrame(rows).sort_values(
        ["aggregate_normalized_importance", "feature"], ascending=[False, True], kind="mergesort"
    )
    table["aggregate_rank"] = np.arange(1, table.shape[0] + 1)
    top = table.head(top_k)
    stable_fraction = float(top["directionally_stable"].mean()) if not top.empty else math.nan
    max_share = float(top["maximum_single_dataset_share"].max()) if not top.empty else math.nan
    passed = bool(
        np.isfinite(stable_fraction)
        and stable_fraction >= float(settings["minimum_stable_fraction"])
        and np.isfinite(max_share)
        and max_share <= float(settings["maximum_single_dataset_share"])
    )
    return table, {
        "pass": passed,
        "top_k": top_k,
        "stable_fraction": stable_fraction,
        "minimum_stable_fraction": float(settings["minimum_stable_fraction"]),
        "maximum_single_dataset_share": max_share,
        "maximum_allowed_single_dataset_share": float(settings["maximum_single_dataset_share"]),
    }


def extension_guardrails_pass(extension_controls: pd.DataFrame, config: Mapping[str, Any]) -> bool:
    if extension_controls.empty:
        return False
    expected_cohorts = 7
    expected_folds = len(config_primary_sources(config))
    balanced = extension_controls[extension_controls["mode"].eq("balanced")]
    paired = balanced[balanced["stratum"].eq("paired_abT")]
    nk = balanced[balanced["stratum"].eq("strict_NK")]
    paired_complete = paired[["dataset_id", "outer_fold_id"]].drop_duplicates().shape[0] == expected_cohorts * expected_folds
    nk_complete = nk[["dataset_id", "outer_fold_id"]].drop_duplicates().shape[0] == expected_cohorts * expected_folds
    return bool(
        paired_complete
        and nk_complete
        and (paired["fpr"] <= float(config["operating_modes"]["balanced"]["maximum_paired_abt_fpr"])).all()
        and (nk["fpr"] <= float(config["operating_modes"]["balanced"]["maximum_strict_nk_fpr"])).all()
    )


def extension_feature_matrix(
    path: Path,
    feature_names: Sequence[str],
    config: Mapping[str, Any],
    row_chunk: int = 20000,
) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray, np.ndarray, float, list[str]]:
    with h5py.File(path, "r") as handle:
        var_names = read_var_names(handle)
        lookup = {gene: index for index, gene in enumerate(var_names)}
        mapping = np.full(len(var_names), -1, dtype=np.int32)
        for output_index, gene in enumerate(feature_names):
            if gene in lookup:
                mapping[lookup[gene]] = output_index
        group = matrix_group(handle, "X")
        n_obs = int(group["indptr"].shape[0] - 1)
        matrix = extract_csr_rows(
            group,
            np.arange(n_obs, dtype=np.int64),
            mapping,
            len(feature_names),
            "raw_counts",
            float(config["normalization"]["target_sum"]),
            row_chunk,
        )
        obs = {column: read_extension_obs(handle, column) for column in EXTENSION_OBS_COLUMNS}
    available = set(var_names)
    missing_critical = sorted(set(config["feature_policy"]["critical_genes"]) - available)
    coverage = float(sum(gene in available for gene in feature_names) / len(feature_names))
    spec = ExtensionFeatureSpec(
        gene_names=list(feature_names),
        gene_indices=np.arange(len(feature_names), dtype=np.int32),
        gene_feature_names=[f"{gene}_log1p_cp10k" for gene in feature_names],
        engineered_feature_names=[],
        model_feature_names=[f"{gene}_log1p_cp10k" for gene in feature_names],
        gene_to_col={gene: index for index, gene in enumerate(feature_names)},
        engineered_to_col={},
    )
    annotation, annotation_source = normalize_extension_annotation(obs, matrix, spec)
    return matrix, obs, annotation, annotation_source, coverage, missing_critical


def score_extension_negative_controls(
    outer_results: Sequence[Mapping[str, Any]],
    feature_names: list[str],
    config: Mapping[str, Any],
) -> pd.DataFrame:
    manifest = pd.read_csv(PROJECT_ROOT / "data/interim/extension_intake/tnk_filtered_h5ads_manifest.csv")
    eligible = manifest[manifest["cohort_id"].ne("GSE169246")].copy()
    rows: list[dict[str, Any]] = []
    stage1_indices = [feature_names.index(gene) for gene in config["feature_policy"]["stage1_genes"]]
    for manifest_row in eligible.itertuples(index=False):
        dataset_id = str(manifest_row.cohort_id)
        path = Path(manifest_row.output_h5ad)
        logging.info("Scoring V4 extension negative controls: %s", dataset_id)
        genes, obs, annotation, annotation_source, coverage, missing_critical = extension_feature_matrix(
            path, feature_names, config
        )
        if missing_critical or coverage < float(config["feature_policy"]["minimum_inference_coverage"]):
            raise RuntimeError(
                f"Extension {dataset_id} is ineligible: coverage={coverage:.3f}, missing critical={missing_critical}"
            )
        truth = extension_truth_frame(obs, annotation, annotation_source)
        paired_mask = truth["paired_abT"].to_numpy(dtype=bool)
        nk_mask = truth["strict_NK"].to_numpy(dtype=bool)
        derived, _ = derive_features(
            genes,
            feature_names,
            config["feature_policy"]["family_prefixes"],
            config["feature_policy"]["derived_panels"],
        )
        _, _, exclusion = exclusion_flags(
            genes,
            feature_names,
            config["cd4_helper_rule"],
            config["treg_rule"],
        )
        stage1_matrix = genes[:, stage1_indices]
        for outer in outer_results:
            if outer.get("stage2_selected") is None:
                continue
            p1 = outer["stage1_calibrator"].predict(
                outer["stage1_final_model"].predict_probability(stage1_matrix)
            )
            stage2_matrix = np.column_stack([genes, derived, p1]).astype(np.float32, copy=False)
            p2 = outer["stage2_calibrator"].predict(
                outer["stage2_final_model"].predict_probability(stage2_matrix)
            )
            stage1_threshold = float(outer["stage1_selected"]["threshold"])
            for mode in ("balanced", "high_purity"):
                decision = outer["stage2_selected"][mode]
                if not bool(decision["passed"]):
                    continue
                call = apply_two_stage_call(
                    p1,
                    p2,
                    stage1_threshold,
                    float(decision["threshold"]),
                    exclusion,
                )
                for stratum, mask in (("paired_abT", paired_mask), ("strict_NK", nk_mask)):
                    total = int(mask.sum())
                    predicted = int((call & mask).sum())
                    rows.append(
                        {
                            "dataset_id": dataset_id,
                            "outer_fold_id": outer["outer_fold_id"],
                            "heldout_source": outer["heldout_source"],
                            "mode": mode,
                            "stratum": stratum,
                            "n_cells": total,
                            "predicted_positive": predicted,
                            "fpr": float(predicted / total) if total else math.nan,
                            "feature_coverage": coverage,
                            "missing_critical_genes": ";".join(missing_critical),
                        }
                    )
        del genes, derived
    output = pd.DataFrame(rows)
    table_path = resolve(config["outputs"]["table_dir"]) / "extension_negative_controls.csv"
    output.to_csv(table_path, index=False)
    return output


def report_table(frame: pd.DataFrame, columns: Sequence[str] | None = None, digits: int = 4) -> str:
    shown = frame.copy() if columns is None else frame[[column for column in columns if column in frame]].copy()
    float_columns = shown.select_dtypes(include=["float", "float32", "float64"]).columns
    shown[float_columns] = shown[float_columns].round(digits)
    return shown.to_html(index=False, border=0, classes="dataframe")


def build_report_figures(config: Mapping[str, Any]) -> list[tuple[str, Path]]:
    table_dir = resolve(config["outputs"]["table_dir"])
    figure_dir = resolve(config["outputs"]["figure_dir"])
    static_dir = resolve(config["outputs"]["static_dir"])
    assets = static_dir / "assets"
    figure_dir.mkdir(parents=True, exist_ok=True)
    assets.mkdir(parents=True, exist_ok=True)
    metrics = pd.read_csv(table_dir / "outer_fold_metrics.csv")
    controls = pd.read_csv(table_dir / "outer_fold_negative_controls.csv")
    bootstrap = pd.read_csv(table_dir / "paired_hierarchical_bootstrap_f1.csv")
    stability = pd.read_csv(table_dir / "feature_stability.csv")
    figures: list[tuple[str, Path]] = []

    balanced = metrics[metrics["mode"].eq("balanced")]
    model_order = [
        "legacy_trd_minus_trab",
        "compact_7gene_logistic",
        "v2_like_tcr_logistic",
        "v4_individual_gene_ablation",
        "v4_nested_selected",
    ]
    sources = config_primary_sources(config)
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2), constrained_layout=True)
    for axis, metric, label in zip(axes, ["f1", "recall", "precision"], ["F1", "Recall", "Precision"], strict=True):
        pivot = balanced.pivot(index="heldout_source", columns="model_id", values=metric).reindex(index=sources, columns=model_order)
        pivot.plot(kind="bar", ax=axis, width=0.82)
        axis.set_ylim(0, 1.02)
        axis.set_xlabel("")
        axis.set_ylabel(label)
        axis.tick_params(axis="x", rotation=20)
        axis.grid(axis="y", alpha=0.25)
        if axis is not axes[-1] and axis.get_legend() is not None:
            axis.get_legend().remove()
    axes[-1].legend(title="Model", fontsize=7, title_fontsize=8, loc="lower left")
    path = figure_dir / "balanced_metrics_by_outer_source.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    figures.append(("Balanced-mode held-out metrics", path))

    fig, axis = plt.subplots(figsize=(9.5, 4.8), constrained_layout=True)
    control = controls[controls["mode"].eq("balanced")].copy()
    control["label"] = control["heldout_source"].astype(str) + " | " + control["model_id"].astype(str)
    axis.barh(control["label"], np.maximum(control["fpr"], 1e-6), color="#3b7f78")
    axis.axvline(float(config["operating_modes"]["balanced"]["maximum_strict_nk_fpr"]), color="#b33a3a", linestyle="--", label="Frozen NK limit")
    axis.set_xscale("log")
    axis.set_xlabel("Strict-NK false-positive rate (log scale)")
    axis.legend()
    axis.grid(axis="x", alpha=0.25)
    path = figure_dir / "heldout_strict_nk_fpr.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    figures.append(("Held-out strict-NK false-positive rates", path))

    if not bootstrap.empty:
        fig, axis = plt.subplots(figsize=(8.0, 3.8), constrained_layout=True)
        y = np.arange(bootstrap.shape[0])
        axis.errorbar(
            bootstrap["mean_f1_difference"],
            y,
            xerr=np.vstack(
                [
                    bootstrap["mean_f1_difference"] - bootstrap["ci95_low"],
                    bootstrap["ci95_high"] - bootstrap["mean_f1_difference"],
                ]
            ),
            fmt="o",
            color="#176b5f",
            capsize=4,
        )
        axis.axvline(0, color="#333333", linewidth=1)
        axis.axvline(0.01, color="#b33a3a", linestyle="--", linewidth=1, label="Required point difference")
        axis.set_yticks(y, bootstrap["comparator"])
        axis.set_xlabel("Dataset/group hierarchical-bootstrap F1 difference (V4 - comparator)")
        axis.legend(fontsize=8)
        axis.grid(axis="x", alpha=0.25)
        path = figure_dir / "paired_bootstrap_f1_differences.png"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        figures.append(("Paired hierarchical-bootstrap F1 differences", path))

    if not stability.empty:
        top = stability.head(30).sort_values("aggregate_normalized_importance")
        fig, axis = plt.subplots(figsize=(8.5, 7.0), constrained_layout=True)
        colors = np.where(top["directionally_stable"], "#176b5f", "#9aa0a6")
        axis.barh(top["feature"], top["aggregate_normalized_importance"], color=colors)
        axis.set_xlabel("Aggregate normalized importance")
        axis.grid(axis="x", alpha=0.25)
        path = figure_dir / "feature_stability_top30.png"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        figures.append(("Selected V4 feature stability", path))

    for _, source in figures:
        shutil.copy2(source, assets / source.name)
    return [(title, assets / path.name) for title, path in figures]


def render_nested_report(config: Mapping[str, Any], no_pdf: bool) -> dict[str, str]:
    table_dir = resolve(config["outputs"]["table_dir"])
    log_dir = resolve(config["outputs"]["log_dir"])
    static_dir = resolve(config["outputs"]["static_dir"])
    static_dir.mkdir(parents=True, exist_ok=True)
    aggregate = load_json(log_dir / "nested_evaluation_aggregate.json")
    metrics = pd.read_csv(table_dir / "outer_fold_metrics.csv")
    macro = pd.read_csv(table_dir / "dataset_macro_metrics.csv")
    controls = pd.read_csv(table_dir / "outer_fold_negative_controls.csv")
    stage1_guardrails = pd.read_csv(table_dir / "stage1_heldout_guardrails.csv")
    extension = pd.read_csv(table_dir / "extension_negative_controls.csv")
    bootstrap = pd.read_csv(table_dir / "paired_hierarchical_bootstrap_f1.csv")
    guardrails = pd.read_csv(table_dir / "promotion_guardrails_step2.csv")
    stability = pd.read_csv(table_dir / "feature_stability.csv")
    figures = build_report_figures(config)
    figure_html = "".join(
        f"<figure><img src='assets/{html.escape(path.name)}' alt='{html.escape(title)}'><figcaption>{html.escape(title)}</figcaption></figure>"
        for title, path in figures
    )
    decision = "Eligible for Step 3 review" if bool(guardrails.iloc[:-1]["pass"].all()) else "Not eligible for promotion under the frozen Step 2 evidence"
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>gdTAI V4 Nested Evaluation</title>
<style>
@page {{ size: A4; margin: 13mm 12mm 15mm; }}
* {{ box-sizing: border-box; }}
body {{ margin: 0; color: #1f2933; font: 10pt/1.42 Arial, sans-serif; background: #f3f6f5; }}
main {{ max-width: 1180px; margin: 0 auto; background: white; padding: 28px 34px 60px; }}
h1 {{ color: #124f47; font-size: 25pt; margin: 0 0 8px; border-bottom: 4px solid #176b5f; padding-bottom: 10px; }}
h2 {{ color: #124f47; font-size: 16pt; margin: 26px 0 8px; border-bottom: 1px solid #c9d8d5; padding-bottom: 5px; }}
h3 {{ color: #176b5f; font-size: 12pt; margin: 18px 0 6px; }}
.status {{ border-left: 5px solid #176b5f; background: #eef6f4; padding: 12px 15px; margin: 16px 0; }}
.warning {{ border-left-color: #b33a3a; background: #fff2f0; }}
.grid {{ display: grid; grid-template-columns: repeat(2, minmax(0,1fr)); gap: 16px; }}
figure {{ margin: 10px 0 18px; break-inside: avoid; }}
figure img {{ width: 100%; height: auto; border: 1px solid #d5dfdd; }}
figcaption {{ font-size: 9pt; color: #52606d; margin-top: 4px; }}
.table-wrap {{ overflow-x: auto; margin: 8px 0 18px; }}
table {{ border-collapse: collapse; width: 100%; font-size: 7.8pt; table-layout: auto; }}
th {{ background: #176b5f; color: white; text-align: left; }}
th,td {{ padding: 4px 5px; border: 1px solid #cbd5d3; vertical-align: top; overflow-wrap: anywhere; }}
tr:nth-child(even) td {{ background: #f3f7f6; }}
code {{ font-size: 8.5pt; }}
@media print {{ body {{ background: white; }} main {{ max-width: none; padding: 0; }} .grid {{ display: block; }} figure {{ page-break-inside: avoid; }} h2 {{ page-break-after: avoid; }} }}
</style></head><body><main>
<h1>gdTAI V4 Nested Cross-Dataset Evaluation</h1>
<p><strong>Protocol:</strong> v{html.escape(str(config['protocol_version']))}<br>
<strong>Design:</strong> nested leave-one-dataset-out with frozen grouped inner folds<br>
<strong>Promotion authorized:</strong> No; Step 2 is a supervision gate</p>
<div class="status {'warning' if not bool(guardrails.iloc[:-1]['pass'].all()) else ''}"><strong>Decision:</strong> {html.escape(decision)}.<br>
Strongest fair comparator: <code>{html.escape(str(aggregate.get('strongest_fair_comparator')))}</code>.
Observed macro-F1 difference: {aggregate.get('observed_macro_f1_difference', math.nan):.4f}; paired-bootstrap lower 95% bound: {aggregate.get('paired_bootstrap_ci95_low', math.nan):.4f}.</div>
<h2>Methodological Contract</h2>
<p>Primary labels use productive paired TRG/TRD without TRA/TRB versus productive paired TRA/TRB without TRG/TRD. Silver, sorted, dual, and ambiguous cells never enter fitting, calibration, threshold selection, or primary performance. Single-chain alpha-beta cells have reliability 0.5 and are training-only. All preprocessing, feature filtering, Stage-1 probabilities, calibration, model selection, and thresholds are learned within the frozen outer-training partition.</p>
<p>V4 uses no legacy TRD/TRAB scores. Stage 1 is a soft T-lineage elastic-net model. Stage 2 compares elastic-net and histogram gradient boosting using 197 frozen genes, eight deterministic summaries, and cross-fitted Stage-1 probability. Calls then apply the immutable CD4-helper/Treg exclusions.</p>
<h2>Promotion Guardrails</h2><div class="table-wrap">{report_table(guardrails)}</div>
<h3>Stage-1 held-out recall</h3><div class="table-wrap">{report_table(stage1_guardrails)}</div>
<h2>Primary Performance</h2><div class="table-wrap">{report_table(macro, ['model_id','mode','n_outer_sources','f1','f0_5','precision','recall','specificity','mcc','pr_auc','roc_auc','brier_score','ece'])}</div>
<div class="grid">{figure_html}</div>
<h3>Per-source metrics</h3><div class="table-wrap">{report_table(metrics, ['heldout_source','model_id','mode','n_positive','n_negative','tp','fp','fn','precision','recall','specificity','f1','mcc','pr_auc','roc_auc'])}</div>
<h2>False-Positive Controls</h2>
<h3>Held-out primary-source NK</h3><div class="table-wrap">{report_table(controls)}</div>
<h3>Frozen extension cohorts</h3><div class="table-wrap">{report_table(extension)}</div>
<h2>Paired Hierarchical Bootstrap</h2><div class="table-wrap">{report_table(bootstrap)}</div>
<h2>Feature Stability</h2><p>Top-30 recurrence, direction agreement, and per-dataset importance shares are evaluated only after held-out predictions are frozen.</p><div class="table-wrap">{report_table(stability.head(40))}</div>
<h2>Interpretation Boundary</h2><p>All positive cohorts have influenced prior development or review. These are nested development estimates, not independent prospective validation. Whole-atlas inference and final fitting remain blocked until this report is reviewed and Step 3 is explicitly authorized.</p>
</main></body></html>"""
    html_path = static_dir / "index.html"
    html_path.write_text(document, encoding="utf-8")
    summary_md = log_dir / "gdtai_v4_nested_evaluation_summary.md"
    summary_md.write_text(
        "# gdTAI V4 Nested Evaluation\n\n"
        f"- Protocol: v{config['protocol_version']}\n"
        f"- Decision: {decision}\n"
        f"- Strongest fair comparator: `{aggregate.get('strongest_fair_comparator')}`\n"
        f"- Macro-F1 difference: `{aggregate.get('observed_macro_f1_difference', math.nan):.6f}`\n"
        f"- Bootstrap 95% lower bound: `{aggregate.get('paired_bootstrap_ci95_low', math.nan):.6f}`\n"
        "- Promotion authorized: No\n",
        encoding="utf-8",
    )
    pdf_path = static_dir / "gdtai_v4_nested_evaluation_report.pdf"
    if not no_pdf:
        try:
            from weasyprint import HTML

            HTML(filename=str(html_path)).write_pdf(str(pdf_path))
        except Exception as error:
            raise RuntimeError(f"PDF rendering failed: {error}") from error
    return {"html": str(html_path), "pdf": str(pdf_path) if pdf_path.exists() else "", "summary": str(summary_md)}


def run_nested_evaluation(config: Mapping[str, Any], heldout_filter: Sequence[str]) -> dict[str, Any]:
    cache_summary = load_json(resolve(config["outputs"]["cache_dir"]) / "feature_cache_manifest.json")
    if cache_summary.get("cell_manifest_sha256") != config["preflight"]["cell_manifest_sha256"]:
        raise RuntimeError("Feature cache does not match the frozen cell manifest")
    cells = pd.read_csv(resolve(config["preflight"]["cell_manifest"]))
    splits = pd.read_csv(resolve(config["preflight"]["split_manifest"]))
    features = pd.read_csv(resolve(config["preflight"]["feature_manifest"])).sort_values("feature_index")
    feature_names = features["gene"].astype(str).tolist()
    gene_matrix = np.load(cache_summary["gene_feature_path"], mmap_mode="r")
    derived_matrix = np.load(cache_summary["derived_feature_path"], mmap_mode="r")
    legacy_score = np.load(cache_summary["legacy_score_path"], mmap_mode="r")
    if gene_matrix.shape != (cells.shape[0], len(feature_names)):
        raise RuntimeError("Cached gene matrix dimensions changed")
    exclusion_audit = apply_frozen_exclusions(cells, gene_matrix, feature_names, config)
    if exclusion_audit != cache_summary.get("exclusion_audit"):
        raise RuntimeError("Full-cell exclusion audit changed after cache construction")
    sources = config_primary_sources(config)
    if heldout_filter:
        unknown = sorted(set(heldout_filter) - set(sources))
        if unknown:
            raise ValueError(f"Unknown held-out sources: {unknown}")
        sources = [source for source in sources if source in set(heldout_filter)]
    summaries: list[dict[str, Any]] = []
    all_metrics: list[dict[str, Any]] = []
    all_controls: list[dict[str, Any]] = []
    bootstrap_counts: dict[str, list[pd.DataFrame]] = {
        "legacy_trd_minus_trab": [],
        "compact_7gene_logistic": [],
        "v2_like_tcr_logistic": [],
    }
    outer_results: list[dict[str, Any]] = []
    feature_importance_frames: list[pd.DataFrame] = []
    stage1_guardrail_rows: list[dict[str, Any]] = []
    for heldout_source in sources:
        result = run_outer_fold(
            cells=cells,
            splits=splits,
            gene_matrix=gene_matrix,
            derived_matrix=derived_matrix,
            legacy_score=legacy_score,
            feature_names=feature_names,
            heldout_source=heldout_source,
            config=config,
        )
        save_outer_result(result, config)
        outer_results.append(result)
        if not result["feature_importance"].empty:
            feature_importance_frames.append(result["feature_importance"])
        stage1_guardrail_rows.append(result["stage1_heldout_guardrails"])
        all_metrics.extend(result["metrics"])
        all_controls.extend(result["control_metrics"])
        for comparator_id, counts in fold_bootstrap_counts(result).items():
            bootstrap_counts[comparator_id].append(counts)
        summaries.append(
            {
                "outer_fold_id": result["outer_fold_id"],
                "heldout_source": heldout_source,
                "stage2_selected": result["stage2_selected"],
                "metrics": result["metrics"],
            }
        )
    output = {
        "protocol_version": config["protocol_version"],
        "cell_manifest_sha256": config["preflight"]["cell_manifest_sha256"],
        "split_manifest_sha256": config["preflight"]["split_manifest_sha256"],
        "outer_folds": summaries,
        "complete": len(summaries) == len(config_primary_sources(config)),
    }
    if output["complete"]:
        extension_controls = score_extension_negative_controls(
            outer_results,
            feature_names,
            config,
        )
        output["aggregate"] = aggregate_nested_results(
            all_metrics,
            all_controls,
            bootstrap_counts,
            extension_controls,
            pd.concat(feature_importance_frames, ignore_index=True)
            if feature_importance_frames
            else pd.DataFrame(),
            pd.DataFrame(stage1_guardrail_rows),
            config,
        )
    write_json(output, resolve(config["outputs"]["log_dir"]) / "nested_evaluation_summary.json")
    return output


def main() -> None:
    args = parse_args()
    config_path = args.config.resolve()
    config = load_json(config_path)
    if not 1 <= int(args.candidate_jobs) <= 80:
        raise ValueError("--candidate-jobs must be between 1 and 80")
    log_dir = resolve(config["outputs"]["log_dir"])
    setup_logging(log_dir, args.log_level)
    logging.info("gdTAI V4 protocol %s, requested stage %s", config["protocol_version"], args.stage)
    verify_hashes = not args.skip_source_hashes
    if args.stage in {"evaluate", "all"} and not verify_hashes:
        raise ValueError("Evaluation cannot skip source-file hashes")
    validation = validate_contract(config, verify_source_hashes=verify_hashes)
    write_json(validation, log_dir / "step2_contract_validation.json")
    if validation["status"] != "PASS":
        raise RuntimeError("Step-2 contract validation failed: " + "; ".join(validation["failures"]))
    logging.info("Frozen Step-2 contract validation passed")
    if args.stage == "validate":
        logging.info("Validation-only stage complete; no model was fitted")
        return
    approval = validate_approval(config)
    logging.info("Step-2 approval accepted: %s at %s", approval["approved_by"], approval["approved_at"])
    runtime_config = dict(config)
    runtime_config["_candidate_jobs"] = int(args.candidate_jobs)
    if args.stage in {"cache", "all"}:
        extract_feature_cache(runtime_config, args.matrix_row_chunk, args.force_cache)
    if args.stage in {"evaluate", "all"}:
        run_nested_evaluation(runtime_config, args.outer_fold)
    if args.stage in {"report", "all"}:
        render_nested_report(runtime_config, args.no_pdf)


if __name__ == "__main__":
    main()
