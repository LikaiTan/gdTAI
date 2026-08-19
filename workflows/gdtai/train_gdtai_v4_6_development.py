#!/usr/bin/env python3
"""Train gdTAI V4.6 with expanded negative diversity and receptor ablation."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import time
from pathlib import Path
from typing import Callable

import cupy as cp
import numpy as np
import pandas as pd

import train_gdtai_v4_4_dual_mode as v44
import train_gdtai_v4_5_development as v45
from train_gdtai_v4_2_nested import sha256_file, source_balanced_weights


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "configs/models/gdtai/v4_6_negative_diversity_receptor_ablation.json"
FIT_CLASSES = ["gdT_gold", "abT_gold", "single_ab_support", "nk_tier1", "nk_tier2"]
ORIGINAL_COMPOSE = v44.compose_stage2


def resolve(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def stable_hash(value: str) -> int:
    return int.from_bytes(hashlib.sha256(value.encode()).digest()[:8], "little")


def load_protocol(path: Path) -> tuple[dict, dict, dict]:
    protocol = json.loads(path.read_text())
    base = json.loads(resolve(protocol["base_config"]).read_text())
    v45_protocol = json.loads(resolve(protocol["v4_5_config"]).read_text())
    base["protocol_version"] = protocol["protocol_version"]
    base["output_table_dir"] = protocol["output_table_dir"]
    base["output_model_dir"] = protocol["output_model_dir"]
    base["output_log_dir"] = protocol["output_log_dir"]
    return protocol, base, v45_protocol


def nk_label_rows(frame: pd.DataFrame) -> pd.DataFrame:
    output = pd.DataFrame({
        "cell_id": frame.cell_id.astype(str),
        "atlas_row": frame.atlas_row.astype(np.int64),
        "truth_class": "nk_tier1",
        "label": 0,
        "stage": "extension_author_nk",
        "reliability_weight": 1.0,
        "label_source": "author_nk_common_benchmark",
        "source_gse_id": frame.source_gse_id.astype(str),
        "donor_id": frame.donor_id.astype(str),
        "sample_id": frame.sample_id.astype(str),
        "library_id": frame.sample_id.astype(str),
        "cohort_role": "locked_v4_6_preparation",
        "allow_fit": False,
        "allow_threshold_selection": False,
        "has_any_ab_tcr": False,
        "has_any_gd_tcr": False,
        "dual_receptor": False,
        "group_id": frame.group_id.astype(str),
        "inner_fold": np.nan,
    })
    return output


def prepare_augmented_cache(protocol: dict, base: dict) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, dict]:
    labels, features, base_matrix, base_cache = v44.load_inputs(base)
    benchmark = pd.read_parquet(resolve(protocol["common_benchmark_predictions"]))
    nk = benchmark[
        benchmark.truth_class.eq("nk_lockbox")
        & benchmark.source_gse_id.isin(protocol["additional_author_nk_sources"])
    ].copy()
    if nk.empty or nk.cell_id.duplicated().any():
        raise RuntimeError("V4.6 author-NK augmentation is empty or duplicated")
    if labels.cell_id.isin(nk.cell_id).any():
        raise RuntimeError("Author-NK augmentation overlaps the frozen label manifest")
    append_labels = nk_label_rows(nk)
    augmented_labels = pd.concat([labels, append_labels], ignore_index=True, sort=False)

    cache_dir = resolve(protocol["augmented_cache_dir"])
    cache_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = cache_dir / "gene_features.npy"
    labels_path = cache_dir / "augmented_label_manifest.parquet"
    manifest_path = cache_dir / "cache_manifest.json"
    common_matrix_path = resolve(protocol["common_benchmark_cache"])
    common_feature_path = resolve(protocol["common_benchmark_feature_manifest"])
    common_names = pd.read_csv(common_feature_path).sort_values("feature_index").gene.astype(str).tolist()
    common_lookup = {gene: index for index, gene in enumerate(common_names)}
    feature_names = features.gene.astype(str).tolist()
    missing = sorted(set(feature_names) - set(common_lookup))
    required = (
        set(base["stage1"]["feature_genes"])
        | set(features.loc[features.feature_class.eq("TCR"), "gene"].astype(str))
        | {"CD3D", "CD3E", "CD3G", "CD4", "FOXP3"}
    )
    for exclusion in base["frozen_exclusions"].values():
        if isinstance(exclusion, dict):
            required.update(exclusion.get("support_genes", []))
            required.update(exclusion.get("program_genes", []))
            if exclusion.get("marker"):
                required.add(exclusion["marker"])
    missing_required = sorted(set(missing) & required)
    if missing_required:
        raise RuntimeError(f"Common benchmark lacks required V4.6 genes: {missing_required}")
    present_target = np.asarray(
        [index for index, gene in enumerate(feature_names) if gene in common_lookup], dtype=np.int64
    )
    common_columns = np.asarray(
        [common_lookup[feature_names[index]] for index in present_target], dtype=np.int64
    )
    common_matrix = np.load(common_matrix_path, mmap_mode="r")
    nk_positions = benchmark.index.get_indexer(nk.index)
    if (nk_positions < 0).any():
        raise RuntimeError("Could not map author-NK rows to the common feature cache")

    expected = {
        "protocol_version": protocol["protocol_version"],
        "base_matrix_sha256": base_cache["matrix_sha256"],
        "common_matrix_sha256": sha256_file(common_matrix_path),
        "common_predictions_sha256": sha256_file(resolve(protocol["common_benchmark_predictions"])),
        "n_base_rows": len(labels),
        "n_author_nk_rows": len(nk),
        "shape": [len(augmented_labels), len(features)],
        "augmentation_unavailable_unused_genes": missing,
    }
    reusable = False
    if manifest_path.exists() and matrix_path.exists() and labels_path.exists():
        previous = json.loads(manifest_path.read_text())
        reusable = all(previous.get(key) == value for key, value in expected.items())
        reusable &= tuple(np.load(matrix_path, mmap_mode="r").shape) == tuple(expected["shape"])
        reusable &= previous.get("augmented_labels_sha256") == sha256_file(labels_path)
    if not reusable:
        temporary = cache_dir / "gene_features.tmp.npy"
        if temporary.exists():
            temporary.unlink()
        target = np.lib.format.open_memmap(
            temporary, mode="w+", dtype=np.float32, shape=tuple(expected["shape"])
        )
        chunk = 100_000
        for start in range(0, len(labels), chunk):
            stop = min(start + chunk, len(labels))
            target[start:stop] = np.asarray(base_matrix[start:stop], dtype=np.float32)
        offset = len(labels)
        target[offset:] = 0.0
        target[offset:, present_target] = np.asarray(
            common_matrix[nk_positions][:, common_columns], dtype=np.float32
        )
        target.flush()
        del target
        os.replace(temporary, matrix_path)
        augmented_labels.to_parquet(labels_path, index=False, compression="zstd")
        manifest = dict(expected)
        manifest.update({
            "status": "PASS_V4_6_AUGMENTED_CACHE",
            "matrix_path": str(matrix_path),
            "matrix_sha256": sha256_file(matrix_path),
            "augmented_labels_sha256": sha256_file(labels_path),
            "normalization": base_cache["normalization"],
        })
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    else:
        augmented_labels = pd.read_parquet(labels_path)
    manifest = json.loads(manifest_path.read_text())
    if manifest["matrix_sha256"] != sha256_file(matrix_path):
        raise RuntimeError("V4.6 augmented cache checksum mismatch")
    matrix = np.load(matrix_path, mmap_mode="r")
    return augmented_labels, features, matrix, manifest


def assign_negative_source_groups(labels: pd.DataFrame, source: str) -> dict[str, str]:
    counts = labels.loc[labels.source_gse_id.eq(source)].groupby("group_id", observed=True).size()
    if counts.empty:
        raise RuntimeError(f"Promoted negative source is absent: {source}")
    total = float(counts.sum())
    targets = {"fit": 0.60 * total, "threshold_validation": 0.20 * total, "calibration": 0.20 * total}
    loads = {key: 0.0 for key in targets}
    assignments: dict[str, str] = {}
    order = sorted(counts.index.astype(str), key=lambda group: (-int(counts.loc[group]), stable_hash(f"v4.6::{source}::{group}")))
    priority = {"fit": 0, "threshold_validation": 1, "calibration": 2}
    for group in order:
        partition = min(
            targets,
            key=lambda name: (loads[name] / max(targets[name], 1.0), priority[name]),
        )
        assignments[group] = partition
        loads[partition] += float(counts.loc[group])
    return assignments


def prepare_partitions(labels: pd.DataFrame, base: dict, v45_protocol: dict, protocol: dict) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    labels, partition = v45.prepare_labels_and_partitions(labels, base, v45_protocol)
    assignment_rows = []
    for source in protocol["promoted_negative_sources"]:
        source_rows = labels.source_gse_id.eq(source)
        labels.loc[source_rows, "allow_fit"] = True
        labels.loc[source_rows, "allow_threshold_selection"] = True
        labels.loc[source_rows, "cohort_role"] = "development_v4_6_consumed_negative"
        assignments = assign_negative_source_groups(labels, source)
        mapped = labels.loc[source_rows, "group_id"].astype(str).map(assignments)
        if mapped.isna().any():
            raise RuntimeError(f"Incomplete V4.6 partition assignment: {source}")
        partition.loc[source_rows] = mapped
        for group, value in assignments.items():
            assignment_rows.append({"source_gse_id": source, "group_id": group, "partition": value})

    development = labels.allow_fit.to_numpy(bool)
    if partition.loc[development].eq("locked_test").any():
        raise RuntimeError("A V4.6 development row remains locked")
    grouped = pd.DataFrame({"group_id": labels.group_id, "partition": partition})
    if grouped.groupby("group_id", observed=True).partition.nunique().max() != 1:
        raise RuntimeError("A biological group crosses V4.6 partitions")
    if labels.loc[development, "truth_class"].astype(str).str.contains("silver", case=False).any():
        raise RuntimeError("Silver positives entered V4.6 development")
    return labels, partition, pd.DataFrame(assignment_rows)


def install_weight_policy(protocol: dict) -> None:
    source_multipliers = {
        source: float(spec["source_weight"])
        for source, spec in protocol["promoted_sorted_positive_sources"].items()
    }
    truth_multipliers = {str(key): float(value) for key, value in protocol["truth_weight_multipliers"].items()}

    def weighted(frame: pd.DataFrame, rows: np.ndarray, stage: int) -> np.ndarray:
        weights = source_balanced_weights(frame, rows, stage).astype(np.float64)
        local = frame.iloc[rows]
        weights *= local.source_gse_id.map(source_multipliers).fillna(1.0).to_numpy(float)
        weights *= local.truth_class.map(truth_multipliers).fillna(1.0).to_numpy(float)
        return (weights / weights.mean()).astype(np.float32)

    v44.source_balanced_weights = weighted
    v44.choose_stage1_threshold = v45.fast_stage1_threshold


def architecture_composer(architecture: str) -> Callable[[np.ndarray, list[str]], tuple[np.ndarray, list[str]]]:
    if architecture == "full_symmetric":
        return ORIGINAL_COMPOSE

    def compose(values: np.ndarray, names: list[str]) -> tuple[np.ndarray, list[str]]:
        full, full_names = ORIGINAL_COMPOSE(values, names)
        lookup = {name: index for index, name in enumerate(full_names)}
        dropped_engineered = {
            "eng_TRGC_TRGV_joint", "eng_gd_receptor_sum", "eng_gd_minus_ab_receptor_sum"
        }
        keep = []
        for index, name in enumerate(full_names):
            standalone_trg = not name.startswith("eng_") and name.startswith("TRG")
            standalone_trg_engineered = name.startswith("eng_TRGV") or name.startswith("eng_TRGC")
            if standalone_trg or standalone_trg_engineered or name in dropped_engineered:
                continue
            keep.append(index)
        output = [full[:, keep]]
        output_names = [full_names[index] for index in keep]
        trdc = values[:, names.index("TRDC")]
        delta_sum = trdc + full[:, lookup["eng_TRDV_all_sum"]] + full[:, lookup["eng_TRDJ_sum"]]
        ab_sum = (
            full[:, lookup["eng_ab_v_sum"]] + full[:, lookup["eng_ab_j_sum"]]
            + full[:, lookup["eng_ab_constant_sum"]]
        )
        engineered = [delta_sum, delta_sum - ab_sum]
        engineered_names = ["eng_delta_receptor_sum", "eng_delta_minus_ab_receptor_sum"]
        if architecture == "conditional_gamma":
            delta_support = np.maximum(trdc, full[:, lookup["eng_TRDV_all_max"]])
            gamma_max = np.maximum(full[:, lookup["eng_TRGV_max"]], full[:, lookup["eng_TRGC_max"]])
            gamma_sum = full[:, lookup["eng_TRGV_sum"]] + full[:, lookup["eng_TRGC_sum"]]
            conditional_max = np.minimum(gamma_max, delta_support)
            conditional_sum = gamma_sum * (delta_support > 0)
            engineered.extend([conditional_max, conditional_sum, delta_sum + conditional_sum])
            engineered_names.extend([
                "eng_gamma_given_delta_max", "eng_gamma_given_delta_sum",
                "eng_conditional_gd_receptor_sum",
            ])
        elif architecture != "delta_dominant":
            raise ValueError(f"Unknown V4.6 architecture: {architecture}")
        output.extend([np.asarray(value, dtype=np.float32)[:, None] for value in engineered])
        return np.column_stack(output).astype(np.float32), output_names + engineered_names

    return compose


def validate_architectures(matrix: np.ndarray, features: pd.DataFrame, base: dict,
                           protocol: dict) -> pd.DataFrame:
    """Freeze and audit the effective feature contract before any model fit."""
    all_names = features.gene.astype(str).tolist()
    config = candidate_config(base)
    _, stage2_columns = v45.feature_columns(features, base, config["stage2"]["control_genes"])
    sample = np.asarray(matrix[: min(32, len(matrix)), :][:, stage2_columns], dtype=np.float32)
    base_names = [all_names[index] for index in stage2_columns]
    rows = []
    maximum = int(protocol["scientific_guardrails"]["maximum_effective_features"])
    for candidate, spec in protocol["feature_candidates"].items():
        architecture = str(spec["architecture"])
        values, names = architecture_composer(architecture)(sample, base_names)
        if values.shape != (len(sample), len(names)) or len(names) != len(set(names)):
            raise RuntimeError(f"Invalid effective feature contract: {candidate}")
        if len(names) > maximum:
            raise RuntimeError(f"Feature limit exceeded before fitting: {candidate}")
        standalone_trg = [name for name in names if not name.startswith("eng_") and name.startswith("TRG")]
        forbidden_gamma = standalone_trg + [
            name for name in names if name.startswith(("eng_TRGV", "eng_TRGC"))
        ]
        if architecture in {"delta_dominant", "conditional_gamma"} and forbidden_gamma:
            raise RuntimeError(f"Standalone gamma evidence remains in {candidate}: {forbidden_gamma}")
        conditional = [name for name in names if name.startswith("eng_gamma_given_delta")]
        if architecture == "conditional_gamma" and len(conditional) != 2:
            raise RuntimeError("Conditional-gamma terms are incomplete")
        if architecture != "conditional_gamma" and conditional:
            raise RuntimeError(f"Conditional-gamma terms leaked into {candidate}")
        rows.append({
            "candidate": candidate,
            "architecture": architecture,
            "n_effective_features": len(names),
            "n_standalone_trg_features": len(standalone_trg),
            "n_forbidden_gamma_features": len(forbidden_gamma),
            "n_conditional_gamma_features": len(conditional),
            "feature_names_sha256": hashlib.sha256("\n".join(names).encode()).hexdigest(),
        })
    return pd.DataFrame(rows)


def candidate_config(base: dict) -> dict:
    config = copy.deepcopy(base)
    config["stage2"]["control_genes"] = ["CD3D", "CD3E", "CD3G", "CD4", "FOXP3"]
    return config


def evaluate_negative_source(bundle: dict, labels: pd.DataFrame, matrix: np.ndarray,
                             all_names: list[str], stage1_columns: np.ndarray,
                             stage2_columns: np.ndarray, source: str) -> pd.DataFrame:
    rows = np.flatnonzero(
        labels.allow_fit.to_numpy() & labels.source_gse_id.eq(source).to_numpy()
        & labels.truth_class.isin(FIT_CLASSES).to_numpy()
    )
    stage1 = bundle["stage1_calibrator"].apply(
        v44.predict_xgb(bundle["stage1_model"], np.asarray(matrix[rows][:, stage1_columns]))
    )
    stage2_values, effective_names = v44.compose_stage2(
        np.asarray(matrix[rows][:, stage2_columns]), bundle["stage2_base_feature_names"]
    )
    if effective_names != bundle["stage2_effective_feature_names"]:
        raise RuntimeError("Candidate feature order changed during negative holdout evaluation")
    stage2 = bundle["stage2_calibrator"].apply(v44.predict_xgb(bundle["stage2_model"], stage2_values))
    excluded = v44.exclusion_for_rows(matrix, rows, all_names)
    score = stage2.copy()
    score[(stage1 < bundle["stage1_threshold"]) | excluded] = -1.0
    frame = labels.iloc[rows].reset_index(drop=True)
    output = []
    for mode, payload in bundle["operating_modes"].items():
        calls = score >= float(payload["threshold"])
        ab = frame.truth_class.eq("abT_gold").to_numpy()
        nk1 = frame.truth_class.eq("nk_tier1").to_numpy()
        nk2 = frame.truth_class.eq("nk_tier2").to_numpy()
        negative = ~frame.truth_class.eq("gdT_gold").to_numpy()
        output.append({
            "heldout_source": source,
            "mode": mode,
            "n_negative": int(negative.sum()),
            "combined_negative_fpr": float(calls[negative].mean()),
            "abt_fpr": float(calls[ab].mean()) if ab.any() else np.nan,
            "tier1_nk_fpr": float(calls[nk1].mean()) if nk1.any() else np.nan,
            "tier2_nk_fpr": float(calls[nk2].mean()) if nk2.any() else np.nan,
            "threshold": float(payload["threshold"]),
            "threshold_selected_without_heldout_source": True,
        })
    return pd.DataFrame(output)


def run_candidate(candidate: str, architecture: str, config: dict, labels: pd.DataFrame,
                  matrix: np.ndarray, all_names: list[str], stage1_columns: np.ndarray,
                  stage2_columns: np.ndarray, partition: pd.Series,
                  protocol: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    positive_tables, negative_tables, validation_tables = [], [], []
    sources = [
        *(('positive', source) for source in protocol["candidate_selection"]["positive_source_holdouts"]),
        *(('negative', source) for source in protocol["candidate_selection"]["negative_source_holdouts"]),
    ]
    v44.compose_stage2 = architecture_composer(architecture)
    for index, (role, source) in enumerate(sources, start=1):
        print(f"[{candidate}] {role} holdout {index}/{len(sources)}: {source}", flush=True)
        bundle = v44.train_bundle(
            labels, matrix, all_names, stage1_columns, stage2_columns,
            partition, config, excluded_source=source,
        )
        if role == "positive":
            result = v44.evaluate_source_holdout(
                bundle, labels, matrix, all_names, stage1_columns, stage2_columns, source
            )
            result["candidate"] = candidate
            positive_tables.append(result)
        else:
            result = evaluate_negative_source(
                bundle, labels, matrix, all_names, stage1_columns, stage2_columns, source
            )
            result["candidate"] = candidate
            negative_tables.append(result)
        validation = bundle["validation_overall"].copy()
        validation["candidate"] = candidate
        validation["excluded_source"] = source
        validation["holdout_role"] = role
        validation_tables.append(validation)
        del bundle
        cp.get_default_memory_pool().free_all_blocks()
    return (
        pd.concat(positive_tables, ignore_index=True),
        pd.concat(negative_tables, ignore_index=True),
        pd.concat(validation_tables, ignore_index=True),
    )


def prevalence_metrics(recall: float, fpr: float, prevalence: float) -> tuple[float, float]:
    precision = recall * prevalence / (recall * prevalence + fpr * (1 - prevalence))
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return precision, f1


def summarize_candidates(positive: pd.DataFrame, negative: pd.DataFrame, protocol: dict) -> pd.DataFrame:
    prevalence = float(protocol["candidate_selection"]["target_prevalence"])
    rows = []
    for candidate in sorted(positive.candidate.unique()):
        item = {"candidate": candidate}
        for mode in ("highest_f1", "high_purity"):
            pos = positive[positive.candidate.eq(candidate) & positive["mode"].eq(mode)]
            neg = negative[negative.candidate.eq(candidate) & negative["mode"].eq(mode)]
            recall = float(pos.gdt_recall.mean())
            fpr = float(neg.combined_negative_fpr.mean())
            precision, f1 = prevalence_metrics(recall, fpr, prevalence)
            item.update({
                f"{mode}_macro_positive_recall": recall,
                f"{mode}_macro_negative_fpr": fpr,
                f"{mode}_precision_at_1pct": precision,
                f"{mode}_f1_at_1pct": f1,
                f"{mode}_minimum_positive_recall": float(pos.gdt_recall.min()),
                f"{mode}_maximum_negative_fpr": float(neg.combined_negative_fpr.max()),
            })
        rows.append(item)
    summary = pd.DataFrame(rows).sort_values(
        ["highest_f1_f1_at_1pct", "high_purity_f1_at_1pct",
         "highest_f1_macro_positive_recall", "highest_f1_macro_negative_fpr"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)
    summary["selected"] = False
    summary.loc[0, "selected"] = True
    summary["selection_objective"] = protocol["candidate_selection"]["objective"]
    return summary


def finalize_contract(summary: dict, protocol: dict, cache_manifest: dict,
                      assignments: pd.DataFrame, config_path: Path) -> dict:
    model_dir = resolve(protocol["output_model_dir"])
    log_dir = resolve(protocol["output_log_dir"])
    contract_path = model_dir / "model_contract.json"
    contract = json.loads(contract_path.read_text())
    input_hashes = dict(contract.get("input_hashes", {}))
    input_hashes.pop("v4_5_config", None)
    input_hashes.update({
        "v4_6_config": sha256_file(config_path),
        "augmented_labels": cache_manifest["augmented_labels_sha256"],
    })
    contract.update({
        "status": "DEVELOPMENT_FROZEN_NOT_PROMOTABLE_NO_NEW_TEST",
        "iteration": "V4.6",
        "promoted_consumed_negative_sources": protocol["promoted_negative_sources"],
        "additional_author_nk_sources": protocol["additional_author_nk_sources"],
        "candidate_selection_protocol": protocol["candidate_selection"],
        "balf_used_for_fit_calibration_selection_or_threshold": False,
        "augmented_cache": {
            "matrix_sha256": cache_manifest["matrix_sha256"],
            "labels_sha256": cache_manifest["augmented_labels_sha256"],
            "n_author_nk_rows": cache_manifest["n_author_nk_rows"],
        },
        "negative_partition_assignments_sha256": sha256_file(
            resolve(protocol["output_table_dir"]) / "promoted_negative_group_partitions.csv"
        ),
        "failure_audit_sha256": sha256_file(
            ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_5_failure_audit/summary.json"
        ),
        "v4_6_config_sha256": sha256_file(config_path),
        "input_hashes": input_hashes,
    })
    contract_path.write_text(json.dumps(v44.json_safe(contract), indent=2, sort_keys=True) + "\n")
    summary.update({
        "model_contract_sha256": sha256_file(contract_path),
        "status": "PASS_V4_6_DEVELOPMENT_FROZEN_NOT_PROMOTABLE",
        "new_untouched_test_scored": False,
        "model_promoted": False,
    })
    (log_dir / "training_summary.json").write_text(
        json.dumps(v44.json_safe(summary), indent=2, sort_keys=True) + "\n"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    started = time.monotonic()
    config_path = resolve(args.config)
    protocol, base, v45_protocol = load_protocol(config_path)
    try:
        cuda_devices = int(cp.cuda.runtime.getDeviceCount())
    except cp.cuda.runtime.CUDARuntimeError:
        cuda_devices = 0
    if not args.preflight_only and (
        cuda_devices < 1 or base["xgboost_fixed"].get("device") != "cuda"
    ):
        raise RuntimeError("V4.6 requires CUDA XGBoost")

    labels, features, matrix, cache_manifest = prepare_augmented_cache(protocol, base)
    labels, partition, assignments = prepare_partitions(labels, base, v45_protocol, protocol)
    install_weight_policy(protocol)
    table_dir = resolve(protocol["output_table_dir"])
    table_dir.mkdir(parents=True, exist_ok=True)
    assignments.to_csv(table_dir / "promoted_negative_group_partitions.csv", index=False)
    architecture_audit = validate_architectures(matrix, features, base, protocol)
    architecture_audit.to_csv(table_dir / "feature_architecture_preflight.csv", index=False)
    counts = labels[["source_gse_id", "truth_class"]].copy()
    counts["partition"] = partition
    counts.groupby(["partition", "source_gse_id", "truth_class"], observed=True).size().reset_index(
        name="n_cells"
    ).to_csv(table_dir / "development_partition_counts.csv", index=False)
    if args.preflight_only:
        preflight = {
            "status": "PASS_V4_6_PREFLIGHT_ONLY",
            "partition_counts": partition.value_counts().to_dict(),
            "development_truth_counts": labels.loc[labels.allow_fit, "truth_class"].value_counts().to_dict(),
            "development_sources": int(labels.loc[labels.allow_fit, "source_gse_id"].nunique()),
            "augmented_cache_shape": list(matrix.shape),
            "author_nk_added": int(cache_manifest["n_author_nk_rows"]),
            "feature_architectures": architecture_audit.to_dict("records"),
            "config_sha256": sha256_file(config_path),
            "augmented_cache_sha256": cache_manifest["matrix_sha256"],
            "balf_used": False,
            "model_fitted": False,
            "cuda_devices_visible_at_preflight": cuda_devices,
        }
        log_dir = resolve(protocol["output_log_dir"])
        log_dir.mkdir(parents=True, exist_ok=True)
        (log_dir / "preflight_summary.json").write_text(
            json.dumps(v44.json_safe(preflight), indent=2, sort_keys=True) + "\n"
        )
        print(json.dumps(v44.json_safe(preflight), indent=2))
        return

    all_names = features.gene.astype(str).tolist()
    config = candidate_config(base)
    stage1_columns, stage2_columns = v45.feature_columns(
        features, base, config["stage2"]["control_genes"]
    )
    positive_tables, negative_tables, validation_tables = [], [], []
    for candidate, spec in protocol["feature_candidates"].items():
        positive, negative, validation = run_candidate(
            candidate, spec["architecture"], config, labels, matrix, all_names,
            stage1_columns, stage2_columns, partition, protocol,
        )
        positive_tables.append(positive)
        negative_tables.append(negative)
        validation_tables.append(validation)
    positive = pd.concat(positive_tables, ignore_index=True)
    negative = pd.concat(negative_tables, ignore_index=True)
    outer_validation = pd.concat(validation_tables, ignore_index=True)
    candidate_summary = summarize_candidates(positive, negative, protocol)
    selected = str(candidate_summary.loc[candidate_summary.selected, "candidate"].iloc[0])
    architecture = protocol["feature_candidates"][selected]["architecture"]
    print(candidate_summary.to_string(index=False), flush=True)

    v44.compose_stage2 = architecture_composer(architecture)
    print(f"Fitting exact V4.6 development scorer: {selected}", flush=True)
    final_bundle = v44.train_bundle(
        labels, matrix, all_names, stage1_columns, stage2_columns, partition, config
    )
    if len(final_bundle["stage2_effective_feature_names"]) > int(
        protocol["scientific_guardrails"]["maximum_effective_features"]
    ):
        raise RuntimeError("V4.6 exceeds the frozen feature limit")
    combined_holdout = pd.concat([
        positive.assign(holdout_role="positive"),
        negative.assign(holdout_role="negative"),
    ], ignore_index=True, sort=False)
    summary = v45.save_final(
        final_bundle, selected, candidate_summary, combined_holdout, outer_validation,
        labels, partition, protocol, base, config_path, cache_manifest, started,
    )
    manifest_path = table_dir / "development_partition_manifest.parquet"
    manifest = pd.read_parquet(manifest_path)
    if "v4_5_partition" not in manifest or "v4_6_partition" in manifest:
        raise RuntimeError("Unexpected inherited partition manifest schema")
    manifest = manifest.rename(columns={"v4_5_partition": "v4_6_partition"})
    manifest.to_parquet(manifest_path, index=False, compression="zstd")
    manifest.groupby(
        ["v4_6_partition", "source_gse_id", "truth_class"], observed=True
    ).size().reset_index(name="n_cells").to_csv(
        table_dir / "development_partition_counts.csv", index=False
    )
    summary = finalize_contract(summary, protocol, cache_manifest, assignments, config_path)
    print(json.dumps(v44.json_safe(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
