#!/usr/bin/env python3
"""Core, testable utilities for the gdTAI V4.2 NK-reference sidecar."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any, Iterable, Sequence

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
BOOL_ROLE_COLUMNS = [
    "include_in_integration_fit",
    "include_in_cluster_label_design",
    "allow_pseudolabel_training",
    "allow_model_tuning",
    "allow_locked_evaluation",
]
MISSING_TEXT = {"", "nan", "none", "na", "unknown", "<na>"}
TCR_PATTERN = re.compile(r"^TR[ABGD](?:V|J|D|C)|^TRAC$|^TRBC|^TRGC|^TRDC$")
OTHER_EXCLUDED_PATTERN = re.compile(r"^(?:MT-|RPS|RPL|IG[HKL])")


def resolve(path: str | Path) -> Path:
    value = Path(path)
    return value if value.is_absolute() else ROOT / value


def read_json(path: str | Path) -> dict[str, Any]:
    return json.loads(resolve(path).read_text(encoding="utf-8"))


def sha256_file(path: str | Path, chunk_size: int = 32 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with resolve(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def decode(value: Any) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def axis_names(handle: h5py.File, axis: str) -> np.ndarray:
    group = handle[axis]
    index_name = decode(group.attrs["_index"])
    return np.asarray([decode(value) for value in group[index_name][:]], dtype=object)


def obs_values(handle: h5py.File, column: str) -> np.ndarray:
    obj = handle["obs"][column]
    if isinstance(obj, h5py.Group) and "categories" in obj and "codes" in obj:
        categories = np.asarray([decode(value) for value in obj["categories"][:]], dtype=object)
        codes = np.asarray(obj["codes"][:], dtype=np.int64)
        values = np.full(codes.size, "", dtype=object)
        valid = codes >= 0
        values[valid] = categories[codes[valid]]
        return values
    raw = np.asarray(obj[:])
    if raw.dtype.kind in {"S", "O", "U"}:
        return np.asarray([decode(value) for value in raw], dtype=object)
    return raw


def bool_values(values: np.ndarray | pd.Series) -> np.ndarray:
    array = np.asarray(values)
    if array.dtype.kind == "b":
        return array.astype(bool)
    if array.dtype.kind in {"i", "u", "f"}:
        return np.nan_to_num(array.astype(float), nan=0.0) != 0
    text = pd.Series(array, dtype="string").fillna("").str.strip().str.lower()
    return text.isin(["true", "1", "yes", "y", "t"]).to_numpy(bool)


def load_roles(config: dict[str, Any]) -> pd.DataFrame:
    frozen = pd.read_csv(resolve(config["preflight"]["frozen_role_contract"]))
    for column in BOOL_ROLE_COLUMNS:
        frozen[column] = bool_values(frozen[column])
    return frozen


def role_contract_sha256(frozen_roles: pd.DataFrame) -> str:
    contract = frozen_roles.drop(columns=["path"]).to_dict(orient="records")
    return canonical_sha256(contract)


def validate_role_separation(roles: pd.DataFrame) -> list[str]:
    errors: list[str] = []
    if roles["cohort_id"].duplicated().any():
        errors.append("cohort roles are not unique")
    locked = roles[roles["allow_locked_evaluation"]]
    leak_columns = [
        "include_in_integration_fit",
        "include_in_cluster_label_design",
        "allow_pseudolabel_training",
        "allow_model_tuning",
    ]
    if locked[leak_columns].any(axis=None):
        errors.append("a locked cohort has a development permission")
    development = roles[roles["include_in_integration_fit"]]
    if development.empty:
        errors.append("no integration-development cohort is defined")
    if development["allow_locked_evaluation"].any():
        errors.append("a development cohort is marked as locked evaluation")
    return errors


def validate_preflight_approval(config: dict[str, Any]) -> dict[str, Any]:
    preflight = config["preflight"]
    approval = read_json(preflight["approval"])
    summary = read_json(preflight["summary"])
    roles = load_roles(config)
    expected = {
        "protocol_version": read_json(preflight["config"])["protocol_version"],
        "config_sha256": sha256_file(preflight["config"]),
        "roles_sha256": sha256_file(preflight["roles"]),
        "cohort_role_contract_sha256": role_contract_sha256(roles),
    }
    errors = validate_role_separation(roles)
    if approval.get("approved") is not True:
        errors.append("preflight implementation approval is not active")
    for key, value in expected.items():
        if approval.get(key) != value:
            errors.append(f"preflight approval mismatch for {key}")
        if summary.get(key) != value:
            errors.append(f"preflight summary mismatch for {key}")
    if summary.get("result") != "PASS_REVIEW_REQUIRED":
        errors.append("preflight result is not PASS_REVIEW_REQUIRED")
    if errors:
        raise RuntimeError("; ".join(errors))
    return {"approval": approval, "summary": summary, "roles": roles, **expected}


def validate_project_data_approval(
    config_path: str | Path,
    config: dict[str, Any],
    runner_path: str | Path,
    core_path: str | Path,
) -> dict[str, Any]:
    approval_path = resolve(config["execution_approval"])
    if not approval_path.exists():
        raise PermissionError(f"Project-data integration approval is absent: {approval_path}")
    approval = read_json(approval_path)
    expected = {
        "protocol_version": config["protocol_version"],
        "execution_config_sha256": sha256_file(config_path),
        "runner_sha256": sha256_file(runner_path),
        "core_sha256": sha256_file(core_path),
    }
    if config.get("current_atlas_recovery", {}).get("active"):
        validate_recovery_preflight(config)
        expected["recovery_contract_sha256"] = recovery_contract_sha256(config)
    errors = []
    if approval.get("approved") is not True:
        errors.append("project-data integration approval is not active")
    for key, value in expected.items():
        if approval.get(key) != value:
            errors.append(f"project-data approval mismatch for {key}")
    if errors:
        raise PermissionError("; ".join(errors))
    return approval


def apply_current_atlas_recovery(config: dict[str, Any], roles: pd.DataFrame) -> pd.DataFrame:
    """Replace a missing integrated input with its frozen raw-count precursor."""
    recovery = config.get("current_atlas_recovery", {})
    if not recovery.get("active"):
        return roles.copy()
    result = roles.copy()
    current_mask = result["cohort_id"].eq("current_atlas")
    if current_mask.sum() != 1:
        raise RuntimeError("Exactly one current_atlas role is required")
    original = Path(result.loc[current_mask, "path"].iloc[0])
    if recovery.get("require_original_missing") and original.exists():
        raise RuntimeError("Recovery is forbidden while the original integrated H5AD exists")
    source = resolve(recovery["source_h5ad"])
    if not source.exists():
        raise FileNotFoundError(f"Recovery source is missing: {source}")
    if int(result.loc[current_mask, "expected_cells"].iloc[0]) != int(recovery["expected_effective_cells"]):
        raise RuntimeError("Recovery effective-cell count differs from the frozen role contract")
    result.loc[current_mask, "path"] = str(source)
    result.loc[current_mask, "expected_sha256"] = str(recovery["source_sha256"])
    result.loc[current_mask, "recovery_active"] = True
    result["recovery_active"] = result["recovery_active"].fillna(False).astype(bool)
    return result


def recovery_contract_payload(config: dict[str, Any]) -> dict[str, Any]:
    recovery = config["current_atlas_recovery"]
    return {
        "source_h5ad": str(resolve(recovery["source_h5ad"])),
        "source_sha256": recovery["source_sha256"],
        "expected_raw_cells": int(recovery["expected_raw_cells"]),
        "expected_effective_cells": int(recovery["expected_effective_cells"]),
        "expected_genes": int(recovery["expected_genes"]),
        "row_exclusion_manifest": str(resolve(recovery["row_exclusion_manifest"])),
        "row_exclusion_manifest_sha256": recovery["row_exclusion_manifest_sha256"],
        "expected_intersecting_exclusions": int(recovery["expected_intersecting_exclusions"]),
        "metadata_sources": [
            {"path": str(resolve(item["path"])), "sha256": item["sha256"]}
            for item in recovery["metadata_sources"]
        ],
        "metadata_columns": list(recovery["metadata_columns"]),
        "expected_metadata_audit": recovery["expected_metadata_audit"],
    }


def recovery_contract_sha256(config: dict[str, Any]) -> str:
    return canonical_sha256(recovery_contract_payload(config))


def validate_recovery_preflight(config: dict[str, Any]) -> dict[str, Any]:
    recovery = config.get("current_atlas_recovery", {})
    if not recovery.get("active"):
        return {"active": False}
    summary_path = resolve(recovery["recovery_preflight_summary"])
    if not summary_path.exists():
        raise RuntimeError(f"Recovery preflight summary is absent: {summary_path}")
    summary = read_json(summary_path)
    expected = recovery_contract_sha256(config)
    if summary.get("result") != "PASS_REVIEW_REQUIRED":
        raise RuntimeError("Recovery preflight did not pass")
    if summary.get("recovery_contract_sha256") != expected:
        raise RuntimeError("Recovery preflight summary does not match the active recovery contract")
    return summary


def is_hvg_excluded(gene: str) -> bool:
    name = str(gene).upper()
    return bool(TCR_PATTERN.match(name) or OTHER_EXCLUDED_PATTERN.match(name))


def common_eligible_genes(gene_sets: Sequence[set[str]]) -> list[str]:
    if not gene_sets:
        raise ValueError("At least one gene set is required")
    common = set.intersection(*gene_sets)
    return sorted(gene for gene in common if not is_hvg_excluded(gene))


def normalize_text(values: pd.Series) -> pd.Series:
    text = values.astype("string").fillna("").str.strip()
    return text.mask(text.str.lower().isin(MISSING_TEXT), "")


def make_integration_batch(obs: pd.DataFrame, hierarchy: Sequence[str]) -> pd.DataFrame:
    if "source_gse_id" not in obs:
        raise KeyError("source_gse_id is required for integration batching")
    source = normalize_text(obs["source_gse_id"])
    if source.eq("").any():
        raise ValueError("source_gse_id contains missing values")
    batch = pd.Series("", index=obs.index, dtype="string")
    level = pd.Series("", index=obs.index, dtype="string")
    for column in hierarchy:
        if column not in obs:
            continue
        values = normalize_text(obs[column])
        take = batch.eq("") & values.ne("")
        batch.loc[take] = source.loc[take] + "::" + column + "::" + values.loc[take]
        level.loc[take] = column
    if batch.eq("").any():
        take = batch.eq("")
        batch.loc[take] = source.loc[take] + "::source_gse_id::" + source.loc[take]
        level.loc[take] = "source_gse_id"
    return pd.DataFrame({"integration_batch": batch, "integration_batch_level": level}, index=obs.index)


def source_balanced_sample_indices(
    sources: Sequence[str], cap_per_source: int, seed: int
) -> np.ndarray:
    if cap_per_source <= 0:
        raise ValueError("cap_per_source must be positive")
    frame = pd.DataFrame({"source": pd.Series(sources, dtype="string"), "row": np.arange(len(sources))})
    if frame["source"].isna().any() or frame["source"].str.strip().eq("").any():
        raise ValueError("source labels must be complete")
    selected: list[np.ndarray] = []
    for offset, (_, group) in enumerate(frame.groupby("source", sort=True, observed=True)):
        rows = group["row"].to_numpy(np.int64)
        if rows.size > cap_per_source:
            rng = np.random.default_rng(seed + offset * 104729)
            rows = np.sort(rng.choice(rows, size=cap_per_source, replace=False))
        selected.append(rows)
    return np.sort(np.concatenate(selected)) if selected else np.asarray([], dtype=np.int64)


def h5ad_obs_frame(path: str | Path, columns: Sequence[str]) -> pd.DataFrame:
    with h5py.File(resolve(path), "r") as handle:
        index = pd.Index(axis_names(handle, "obs"), name="original_cell_id")
        data: dict[str, Any] = {}
        for column in columns:
            if column in handle["obs"]:
                data[column] = obs_values(handle, column)
            else:
                data[column] = np.full(index.size, "", dtype=object)
    return pd.DataFrame(data, index=index)


def h5ad_var_names(path: str | Path) -> np.ndarray:
    with h5py.File(resolve(path), "r") as handle:
        return axis_names(handle, "var")


def recovery_excluded_cell_ids(config: dict[str, Any]) -> set[str]:
    recovery = config["current_atlas_recovery"]
    manifest = pd.read_csv(resolve(recovery["row_exclusion_manifest"]), usecols=["obs_name"])
    return set(manifest["obs_name"].astype(str))


def load_recovery_metadata(config: dict[str, Any]) -> pd.DataFrame:
    recovery = config["current_atlas_recovery"]
    usecols = ["source_gse_id", "original_cell_id", *recovery["metadata_columns"]]
    frames = []
    for item in recovery["metadata_sources"]:
        path = resolve(item["path"])
        header = pd.read_csv(path, nrows=0)
        missing = sorted(set(usecols) - set(header.columns))
        if missing:
            raise KeyError(f"Recovery metadata source {path} lacks columns: {missing}")
        frames.append(pd.read_csv(path, usecols=usecols, dtype="string", low_memory=False))
    metadata = pd.concat(frames, ignore_index=True, copy=False)
    for column in usecols:
        metadata[column] = normalize_text(metadata[column])
    metadata["metadata_key"] = (
        metadata["source_gse_id"].fillna("missing_gse")
        + "||"
        + metadata["original_cell_id"].fillna("missing_cell")
    )
    if metadata["metadata_key"].duplicated().any():
        raise RuntimeError("Recovery metadata contains duplicate source-plus-cell keys")
    return metadata.set_index("metadata_key")


def attach_recovery_metadata(frame: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    required = ["source_gse_id", "source_original_cell_id"]
    missing = [column for column in required if column not in frame]
    if missing:
        raise KeyError(f"Recovery source metadata lacks columns: {missing}")
    metadata = load_recovery_metadata(config)
    key = (
        normalize_text(frame["source_gse_id"]).replace("", "missing_gse")
        + "||"
        + normalize_text(frame["source_original_cell_id"]).replace("", "missing_cell")
    )
    joined = metadata.reindex(key)
    result = frame.copy()
    for column in config["current_atlas_recovery"]["metadata_columns"]:
        result[column] = joined[column].to_numpy()
    return result


def _matrix_group(handle: h5py.File, key: str) -> h5py.Group:
    value: Any = handle
    for component in key.split("/"):
        value = value[component]
    if not isinstance(value, h5py.Group):
        raise TypeError(f"{key} is not an H5AD sparse-matrix group")
    return value


def read_sparse_rows_genes(
    path: str | Path,
    matrix_key: str,
    genes: Sequence[str],
    rows: np.ndarray | None = None,
    row_chunk_size: int = 25000,
    dtype: np.dtype = np.float32,
) -> sparse.csr_matrix:
    """Read selected H5AD sparse rows/genes without materializing all genes."""
    rows = None if rows is None else np.asarray(rows, dtype=np.int64)
    if rows is not None and (rows.size and (np.any(np.diff(rows) < 0) or rows[0] < 0)):
        raise ValueError("rows must be sorted nonnegative indices")
    with h5py.File(resolve(path), "r") as handle:
        var_names = axis_names(handle, "var")
        lookup = {gene: index for index, gene in enumerate(var_names)}
        missing = [gene for gene in genes if gene not in lookup]
        if missing:
            raise KeyError(f"Missing requested genes: {missing[:10]}")
        columns = np.asarray([lookup[gene] for gene in genes], dtype=np.int64)
        group = _matrix_group(handle, matrix_key)
        encoding = decode(group.attrs.get("encoding-type", ""))
        shape = tuple(int(value) for value in group.attrs["shape"])
        if rows is not None and rows.size and rows[-1] >= shape[0]:
            raise IndexError("row index exceeds matrix shape")

        blocks: list[sparse.csr_matrix] = []
        if encoding == "csr_matrix":
            for start in range(0, shape[0], row_chunk_size):
                end = min(shape[0], start + row_chunk_size)
                if rows is None:
                    local_rows: slice | np.ndarray = slice(None)
                else:
                    lo = np.searchsorted(rows, start, side="left")
                    hi = np.searchsorted(rows, end, side="left")
                    if lo == hi:
                        continue
                    local_rows = rows[lo:hi] - start
                pointers = np.asarray(group["indptr"][start : end + 1], dtype=np.int64)
                data_start, data_end = int(pointers[0]), int(pointers[-1])
                block = sparse.csr_matrix(
                    (
                        np.asarray(group["data"][data_start:data_end], dtype=dtype),
                        np.asarray(group["indices"][data_start:data_end], dtype=np.int64),
                        pointers - data_start,
                    ),
                    shape=(end - start, shape[1]),
                )
                blocks.append(block[local_rows, :][:, columns].tocsr())
        elif encoding == "csc_matrix":
            full = sparse.csc_matrix(
                (
                    np.asarray(group["data"][:], dtype=dtype),
                    np.asarray(group["indices"][:], dtype=np.int64),
                    np.asarray(group["indptr"][:], dtype=np.int64),
                ),
                shape=shape,
            )
            selected = full[:, columns]
            selected = selected if rows is None else selected[rows, :]
            blocks.append(selected.tocsr())
        else:
            raise TypeError(f"Unsupported sparse encoding: {encoding}")
    if not blocks:
        return sparse.csr_matrix((0, len(genes)), dtype=dtype)
    result = sparse.vstack(blocks, format="csr", dtype=dtype)
    result.sort_indices()
    return result


def mixed_boundary_mask(
    partitions: np.ndarray,
    primary_nk_anchor: np.ndarray,
    productive_tcr_anchor: np.ndarray,
    minimum_mixed_runs: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    labels = np.asarray(partitions)
    if labels.ndim != 2:
        raise ValueError("partitions must be a cells-by-runs matrix")
    nk = np.asarray(primary_nk_anchor, dtype=bool)
    tcr = np.asarray(productive_tcr_anchor, dtype=bool)
    if labels.shape[0] != nk.size or nk.size != tcr.size:
        raise ValueError("partition and anchor lengths differ")
    mixed_count = np.zeros(labels.shape[0], dtype=np.int16)
    for run in range(labels.shape[1]):
        values = labels[:, run].astype(np.int64)
        valid = values >= 0
        n_clusters = int(values[valid].max()) + 1 if valid.any() else 0
        nk_count = np.bincount(values[nk & valid], minlength=n_clusters)
        tcr_count = np.bincount(values[tcr & valid], minlength=n_clusters)
        mixed = (nk_count > 0) & (tcr_count > 0)
        mixed_count[valid] += mixed[values[valid]].astype(np.int16)
    return mixed_count >= minimum_mixed_runs, mixed_count


def validate_pseudolabel_contract(contract: dict[str, Any]) -> None:
    if contract.get("marker_thresholds_may_define_truth") is not False:
        raise ValueError("marker thresholds must not define pseudo-label truth")
    if contract.get("may_control_validation_guardrails") is not False:
        raise ValueError("pseudo-labels must not control validation guardrails")
    if not 0.5 <= float(contract["minimum_run_agreement"]) <= 1.0:
        raise ValueError("minimum_run_agreement is outside [0.5, 1]")
    if not 0.5 <= float(contract["minimum_anchor_nk_purity"]) <= 1.0:
        raise ValueError("minimum_anchor_nk_purity is outside [0.5, 1]")
    if not 0 <= float(contract["maximum_productive_tcr_contamination"]) <= 0.5:
        raise ValueError("maximum_productive_tcr_contamination is outside [0, 0.5]")
    if int(contract["minimum_development_sources"]) < 2:
        raise ValueError("minimum_development_sources must be at least two")


def cluster_consensus_votes(
    partitions: np.ndarray,
    primary_nk_anchor: np.ndarray,
    productive_tcr_anchor: np.ndarray,
    candidate_eligible: np.ndarray,
    candidate_sources: Sequence[str],
    contract: dict[str, Any],
    run_names: Sequence[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Evaluate independent-anchor cluster rules and return cell votes/stats."""
    validate_pseudolabel_contract(contract)
    labels = np.asarray(partitions, dtype=np.int64)
    n_cells, n_runs = labels.shape
    nk = np.asarray(primary_nk_anchor, dtype=bool)
    tcr = np.asarray(productive_tcr_anchor, dtype=bool)
    candidate = np.asarray(candidate_eligible, dtype=bool)
    sources = pd.Series(candidate_sources, dtype="string").fillna("").to_numpy(str)
    if any(values.size != n_cells for values in [nk, tcr, candidate, sources]):
        raise ValueError("cell-level consensus inputs have different lengths")
    run_names = list(run_names or [f"run_{index}" for index in range(n_runs)])
    if len(run_names) != n_runs:
        raise ValueError("run_names length differs from partition count")

    source_levels = np.sort(np.unique(sources[candidate]))
    source_lookup = {value: index for index, value in enumerate(source_levels)}
    source_code = np.full(n_cells, -1, dtype=np.int32)
    source_code[candidate] = np.asarray([source_lookup[value] for value in sources[candidate]], dtype=np.int32)
    votes = np.zeros(n_cells, dtype=np.int16)
    available = np.zeros(n_cells, dtype=np.int16)
    rows: list[dict[str, Any]] = []

    for run_index, run_name in enumerate(run_names):
        values = labels[:, run_index]
        valid = values >= 0
        if not valid.any():
            continue
        n_clusters = int(values[valid].max()) + 1
        total = np.bincount(values[valid], minlength=n_clusters)
        nk_count = np.bincount(values[nk & valid], minlength=n_clusters)
        tcr_count = np.bincount(values[tcr & valid], minlength=n_clusters)
        anchor_count = nk_count + tcr_count
        purity = np.divide(nk_count, anchor_count, out=np.zeros(n_clusters, dtype=float), where=anchor_count > 0)
        contamination = np.divide(tcr_count, total, out=np.ones(n_clusters, dtype=float), where=total > 0)

        source_counts = np.zeros((n_clusters, len(source_levels)), dtype=np.int64)
        candidate_valid = candidate & valid
        if candidate_valid.any() and source_levels.size:
            flat = values[candidate_valid] * len(source_levels) + source_code[candidate_valid]
            source_counts = np.bincount(flat, minlength=n_clusters * len(source_levels)).reshape(n_clusters, len(source_levels))
        candidate_count = source_counts.sum(axis=1)
        source_n = (source_counts > 0).sum(axis=1)
        max_source = source_counts.max(axis=1) if source_levels.size else np.zeros(n_clusters, dtype=np.int64)
        max_fraction = np.divide(max_source, candidate_count, out=np.ones(n_clusters, dtype=float), where=candidate_count > 0)

        qualifies = (
            (anchor_count >= int(contract["minimum_independent_anchors_per_cluster"]))
            & (nk_count >= int(contract["minimum_primary_nk_anchors_per_cluster"]))
            & (purity >= float(contract["minimum_anchor_nk_purity"]))
            & (contamination <= float(contract["maximum_productive_tcr_contamination"]))
            & (source_n >= int(contract["minimum_development_sources"]))
            & (max_fraction <= float(contract["maximum_single_source_fraction"]))
        )
        available[candidate_valid] += 1
        votes[candidate_valid] += qualifies[values[candidate_valid]].astype(np.int16)
        for cluster in range(n_clusters):
            rows.append(
                {
                    "run": run_name,
                    "cluster": cluster,
                    "n_cells": int(total[cluster]),
                    "n_primary_nk_anchors": int(nk_count[cluster]),
                    "n_productive_tcr_anchors": int(tcr_count[cluster]),
                    "anchor_nk_purity": float(purity[cluster]),
                    "productive_tcr_contamination": float(contamination[cluster]),
                    "n_candidate_cells": int(candidate_count[cluster]),
                    "n_candidate_sources": int(source_n[cluster]),
                    "maximum_candidate_source_fraction": float(max_fraction[cluster]),
                    "qualifies_as_nk": bool(qualifies[cluster]),
                }
            )

    agreement = np.divide(votes, available, out=np.zeros(n_cells, dtype=float), where=available > 0)
    selected = candidate & (available > 0) & (agreement >= float(contract["minimum_run_agreement"]))
    cell = pd.DataFrame(
        {
            "candidate_eligible": candidate,
            "qualifying_votes": votes,
            "available_runs": available,
            "run_agreement": agreement,
            "selected_pseudo_nk": selected,
        }
    )
    return cell, pd.DataFrame(rows)


def atomic_write_json(path: str | Path, value: Any) -> None:
    output = resolve(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".partial")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(output)


def input_file_state(roles: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in roles.itertuples(index=False):
        path = Path(row.path)
        stat = path.stat()
        rows.append(
            {
                "cohort_id": row.cohort_id,
                "path": str(path),
                "size_bytes": int(stat.st_size),
                "mtime_ns": int(stat.st_mtime_ns),
                "expected_sha256": row.expected_sha256,
            }
        )
    return pd.DataFrame(rows)


def ensure_no_locked_cohorts(frame: pd.DataFrame, roles: pd.DataFrame) -> None:
    locked = set(roles.loc[roles["allow_locked_evaluation"], "cohort_id"])
    observed = set(frame["input_cohort_id"].astype(str).unique())
    overlap = sorted(locked & observed)
    if overlap:
        raise RuntimeError(f"Locked cohorts leaked into development data: {overlap}")


def iter_parameter_grid(resolutions: Iterable[float], seeds: Iterable[int]) -> list[tuple[float, int]]:
    return [(float(resolution), int(seed)) for resolution in resolutions for seed in seeds]
