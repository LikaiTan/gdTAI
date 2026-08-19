#!/usr/bin/env python3
"""Run the precommitted gdTAI V4.5 development-only feature ablation."""

from __future__ import annotations

import argparse
import copy
import json
import math
import time
from pathlib import Path
from typing import Any

import cupy as cp
import numpy as np
import pandas as pd

import train_gdtai_v4_4_dual_mode as v44
from train_gdtai_v4_2_nested import sha256_file, source_balanced_weights


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "configs/models/gdtai/v4_5_positive_diversity_ablation.json"


def resolve(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def json_safe(value: Any) -> Any:
    return v44.json_safe(value)


def load_protocol(path: Path) -> tuple[dict, dict]:
    protocol = json.loads(path.read_text())
    base_path = resolve(protocol["base_config"])
    base = json.loads(base_path.read_text())
    base["protocol_version"] = protocol["protocol_version"]
    base["output_table_dir"] = protocol["output_table_dir"]
    base["output_model_dir"] = protocol["output_model_dir"]
    base["output_log_dir"] = protocol["output_log_dir"]
    return protocol, base


def prepare_labels_and_partitions(
    labels: pd.DataFrame, base: dict, protocol: dict
) -> tuple[pd.DataFrame, pd.Series]:
    labels = labels.copy()
    partition = v44.assign_development_partitions(labels, base)
    promoted = protocol["promoted_sorted_positive_sources"]
    for source, spec in promoted.items():
        source_rows = labels.source_gse_id.eq(source)
        if not source_rows.any():
            raise RuntimeError(f"Promoted positive source is absent: {source}")
        if not labels.loc[source_rows, "truth_class"].eq("gdT_gold").all():
            raise RuntimeError(f"Promoted source contains non-gdT truth: {source}")
        labels.loc[source_rows, "allow_fit"] = True
        labels.loc[source_rows, "allow_threshold_selection"] = True
        labels.loc[source_rows, "cohort_role"] = "development_v4_5_consumed_sorted"
        assignments = spec["group_partitions"]
        observed = set(labels.loc[source_rows, "group_id"].astype(str))
        if observed != set(assignments):
            raise RuntimeError(
                f"Frozen group map mismatch for {source}: missing={sorted(observed-set(assignments))}, "
                f"extra={sorted(set(assignments)-observed)}"
            )
        partition.loc[source_rows] = labels.loc[source_rows, "group_id"].map(assignments)

    development = labels.allow_fit.to_numpy(bool)
    if partition.loc[development].eq("locked_test").any():
        raise RuntimeError("A V4.5 development cell remains in locked_test")
    if labels.loc[development, "truth_class"].astype(str).str.contains("silver", case=False).any():
        raise RuntimeError("Silver positives entered V4.5 development")
    grouped = pd.DataFrame({"group_id": labels.group_id, "partition": partition})
    if grouped.groupby("group_id", observed=True).partition.nunique().max() != 1:
        raise RuntimeError("A biological group crosses V4.5 partitions")
    return labels, partition


def install_weight_policy(labels: pd.DataFrame, protocol: dict) -> None:
    source_multipliers = {
        source: float(spec["source_weight"])
        for source, spec in protocol["promoted_sorted_positive_sources"].items()
    }
    truth_multipliers = {
        str(key): float(value) for key, value in protocol["truth_weight_multipliers"].items()
    }

    def weighted(frame: pd.DataFrame, rows: np.ndarray, stage: int) -> np.ndarray:
        weights = source_balanced_weights(frame, rows, stage).astype(np.float64)
        local = frame.iloc[rows]
        weights *= local.source_gse_id.map(source_multipliers).fillna(1.0).to_numpy(float)
        weights *= local.truth_class.map(truth_multipliers).fillna(1.0).to_numpy(float)
        return (weights / weights.mean()).astype(np.float32)

    v44.source_balanced_weights = weighted


def feature_columns(
    features: pd.DataFrame, base: dict, context_genes: list[str]
) -> tuple[np.ndarray, np.ndarray]:
    names = features.gene.astype(str)
    stage1_names = set(base["stage1"]["feature_genes"])
    stage2_context = set(context_genes)
    stage1 = np.flatnonzero(names.isin(stage1_names).to_numpy())
    stage2 = np.flatnonzero(
        features.feature_class.eq("TCR").to_numpy() | names.isin(stage2_context).to_numpy()
    )
    missing = sorted((stage1_names - set(names.iloc[stage1])) | (stage2_context - set(names.iloc[stage2])))
    if missing:
        raise RuntimeError(f"Candidate features missing from cache: {missing}")
    return stage1, stage2


def candidate_config(base: dict, context_genes: list[str]) -> dict:
    config = copy.deepcopy(base)
    config["stage2"]["control_genes"] = list(context_genes)
    return config


def finite_mean(values: pd.Series) -> float:
    values = pd.to_numeric(values, errors="coerce")
    return float(values[np.isfinite(values)].mean()) if np.isfinite(values).any() else math.nan


def run_candidate_holdouts(
    candidate: str,
    config: dict,
    labels: pd.DataFrame,
    matrix: np.ndarray,
    all_names: list[str],
    stage1_columns: np.ndarray,
    stage2_columns: np.ndarray,
    partition: pd.Series,
    sources: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    heldout_tables: list[pd.DataFrame] = []
    validation_tables: list[pd.DataFrame] = []
    for index, source in enumerate(sources, start=1):
        print(f"[{candidate}] source holdout {index}/{len(sources)}: {source}", flush=True)
        bundle = v44.train_bundle(
            labels, matrix, all_names, stage1_columns, stage2_columns,
            partition, config, excluded_source=source,
        )
        heldout = v44.evaluate_source_holdout(
            bundle, labels, matrix, all_names, stage1_columns, stage2_columns, source
        )
        heldout["candidate"] = candidate
        heldout_tables.append(heldout)
        local_validation = bundle["validation_overall"].copy()
        local_validation["candidate"] = candidate
        local_validation["excluded_source"] = source
        validation_tables.append(local_validation)
        del bundle
        cp.get_default_memory_pool().free_all_blocks()
    return pd.concat(heldout_tables, ignore_index=True), pd.concat(validation_tables, ignore_index=True)


def summarize_candidates(heldout: pd.DataFrame, protocol: dict) -> pd.DataFrame:
    primary_mode = protocol["candidate_selection"]["primary_mode"]
    selected = heldout.loc[heldout["mode"].eq(primary_mode)].copy()
    rows = []
    for candidate, frame in selected.groupby("candidate", observed=True):
        rows.append({
            "candidate": candidate,
            "macro_source_holdout_f1": finite_mean(frame["f1"]),
            "macro_source_holdout_recall": finite_mean(frame["gdt_recall"]),
            "macro_source_holdout_abt_fpr": finite_mean(frame["abt_fpr"]),
            "minimum_source_holdout_recall": float(frame.gdt_recall.min()),
            "n_positive_source_holdouts": int(len(frame)),
        })
    summary = pd.DataFrame(rows).sort_values(
        ["macro_source_holdout_f1", "macro_source_holdout_recall", "macro_source_holdout_abt_fpr"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    summary["selected"] = False
    summary.loc[0, "selected"] = True
    summary["selection_objective"] = protocol["candidate_selection"]["objective"]
    return summary


def save_final(
    bundle: dict,
    candidate: str,
    candidate_summary: pd.DataFrame,
    heldout: pd.DataFrame,
    outer_validation: pd.DataFrame,
    labels: pd.DataFrame,
    partition: pd.Series,
    protocol: dict,
    base: dict,
    config_path: Path,
    cache_manifest: dict,
    started: float,
) -> dict:
    table_dir = resolve(protocol["output_table_dir"])
    model_dir = resolve(protocol["output_model_dir"])
    log_dir = resolve(protocol["output_log_dir"])
    for path in (table_dir, model_dir, log_dir):
        path.mkdir(parents=True, exist_ok=True)

    candidate_summary.to_csv(table_dir / "candidate_selection_summary.csv", index=False)
    heldout.to_csv(table_dir / "candidate_source_holdout_metrics.csv", index=False)
    outer_validation.to_csv(table_dir / "candidate_outer_validation_metrics.csv", index=False)
    bundle["frontier"].to_csv(table_dir / "threshold_validation_frontier.csv", index=False)
    bundle["validation_overall"].to_csv(table_dir / "threshold_validation_metrics.csv", index=False)
    bundle["validation_per_source"].to_csv(table_dir / "threshold_validation_per_source.csv", index=False)

    manifest = labels[[
        "cell_id", "atlas_row", "source_gse_id", "truth_class", "group_id", "cohort_role"
    ]].copy()
    manifest["v4_5_partition"] = partition
    manifest.to_parquet(table_dir / "development_partition_manifest.parquet", index=False, compression="zstd")
    manifest.groupby(
        ["v4_5_partition", "source_gse_id", "truth_class"], observed=True
    ).size().reset_index(name="n_cells").to_csv(table_dir / "development_partition_counts.csv", index=False)

    rows = bundle["threshold_rows"]
    predictions = labels.iloc[rows][
        ["cell_id", "atlas_row", "source_gse_id", "truth_class", "group_id"]
    ].copy()
    predictions["stage1_probability"] = bundle["threshold_stage1_score"]
    predictions["stage2_probability"] = bundle["threshold_stage2_score"]
    predictions["effective_score"] = bundle["threshold_effective_score"]
    predictions["cd4_or_treg_excluded"] = bundle["threshold_excluded"]
    low = float(bundle["operating_modes"]["highest_f1"]["threshold"])
    high = float(bundle["operating_modes"]["high_purity"]["threshold"])
    if high < low:
        raise RuntimeError("High-purity threshold is below highest-F1 threshold")
    predictions["predicted_highest_f1"] = predictions.effective_score.ge(low)
    predictions["predicted_high_purity"] = predictions.effective_score.ge(high)
    predictions["decision_band"] = np.select(
        [predictions.effective_score.ge(high), predictions.effective_score.ge(low)],
        ["high_confidence_gdt", "uncertain_gdt"], default="negative"
    )
    predictions.to_parquet(
        table_dir / "threshold_validation_predictions.parquet", index=False, compression="zstd"
    )
    band_summary = predictions.groupby(
        ["decision_band", "truth_class"], observed=True
    ).size().reset_index(name="n_cells")
    band_summary["fraction_within_truth"] = band_summary.n_cells / band_summary.groupby(
        "truth_class", observed=True
    ).n_cells.transform("sum")
    band_summary.to_csv(table_dir / "uncertainty_band_summary.csv", index=False)

    feature_contract = pd.DataFrame({"feature_name": bundle["stage2_effective_feature_names"]})
    feature_contract["feature_type"] = np.where(
        feature_contract.feature_name.str.startswith("eng_"),
        "engineered_receptor_aggregate", "individual_gene"
    )
    feature_contract.to_csv(table_dir / "stage2_feature_contract.csv", index=False)

    stage1_path = model_dir / "stage1_t_lineage_gate.ubj"
    stage2_path = model_dir / "stage2_gdt_classifier.ubj"
    bundle["stage1_model"].save_model(stage1_path)
    bundle["stage2_model"].save_model(stage2_path)
    calibration_path = model_dir / "platt_calibration.json"
    calibration_path.write_text(json.dumps({
        "stage1": bundle["stage1_calibrator"].as_dict(),
        "stage2": bundle["stage2_calibrator"].as_dict(),
    }, indent=2, sort_keys=True) + "\n")

    contract = {
        "protocol_version": protocol["protocol_version"],
        "status": "DEVELOPMENT_FROZEN_NOT_PROMOTABLE_NO_NEW_TEST",
        "selected_candidate": candidate,
        "normalization": "log1p(raw_counts_per_10000)",
        "development_frozen": True,
        "new_untouched_test_scored": False,
        "model_promoted": False,
        "consumed_sorted_sources_reclassified_as_development": sorted(
            protocol["promoted_sorted_positive_sources"]
        ),
        "silver_used": False,
        "operating_modes": bundle["operating_modes"],
        "uncertainty_policy": protocol["uncertainty_policy"],
        "stage1_threshold": bundle["stage1_threshold"],
        "stage1_feature_names": bundle["stage1_feature_names"],
        "stage2_base_feature_names": bundle["stage2_base_feature_names"],
        "stage2_effective_feature_names": bundle["stage2_effective_feature_names"],
        "n_effective_features": len(bundle["stage2_effective_feature_names"]),
        "hard_trdv_expression_cutoff": False,
        "hard_single_nk_gene_exclusion": False,
        "frozen_cd4_treg_exclusions": base["frozen_exclusions"],
        "candidate_selection": candidate_summary.to_dict("records"),
        "source_holdout_metrics": heldout.to_dict("records"),
        "threshold_validation_metrics": bundle["validation_overall"].to_dict("records"),
        "scientific_guardrails": protocol["scientific_guardrails"],
        "input_hashes": {
            "v4_5_config": sha256_file(config_path),
            "labels": sha256_file(resolve(base["label_manifest"])),
            "features": sha256_file(resolve(base["feature_manifest"])),
            "feature_cache": cache_manifest["matrix_sha256"],
        },
        "artifact_hashes": {
            "stage1_model": sha256_file(stage1_path),
            "stage2_model": sha256_file(stage2_path),
            "platt_calibration": sha256_file(calibration_path),
        },
    }
    contract_path = model_dir / "model_contract.json"
    contract_path.write_text(json.dumps(json_safe(contract), indent=2, sort_keys=True) + "\n")
    summary = {
        "status": "PASS_DEVELOPMENT_FROZEN_NOT_PROMOTABLE",
        "runtime_seconds": time.monotonic() - started,
        "selected_candidate": candidate,
        "model_contract": str(contract_path),
        "model_contract_sha256": sha256_file(contract_path),
        "operating_modes": bundle["operating_modes"],
        "new_untouched_test_scored": False,
        "model_promoted": False,
    }
    (log_dir / "training_summary.json").write_text(
        json.dumps(json_safe(summary), indent=2, sort_keys=True) + "\n"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    started = time.monotonic()
    config_path = resolve(args.config)
    protocol, base = load_protocol(config_path)
    if cp.cuda.runtime.getDeviceCount() < 1:
        raise RuntimeError("No CUDA device is available")
    if base["xgboost_fixed"].get("device") != "cuda":
        raise RuntimeError("V4.5 requires GPU XGBoost")

    labels, features, matrix, cache_manifest = v44.load_inputs(base)
    labels, partition = prepare_labels_and_partitions(labels, base, protocol)
    install_weight_policy(labels, protocol)
    table_dir = resolve(protocol["output_table_dir"])
    table_dir.mkdir(parents=True, exist_ok=True)
    preview = labels[["source_gse_id", "truth_class"]].copy()
    preview["partition"] = partition
    preview.groupby(
        ["partition", "source_gse_id", "truth_class"], observed=True
    ).size().reset_index(name="n_cells").to_csv(
        table_dir / "development_partition_counts.csv", index=False
    )
    if args.preflight_only:
        print(json.dumps({
            "status": "PASS_PREFLIGHT_ONLY",
            "partition_counts": partition.value_counts().to_dict(),
            "development_positive_sources": sorted(
                labels.loc[labels.allow_fit & labels.truth_class.eq("gdT_gold"), "source_gse_id"].unique()
            ),
        }, indent=2))
        return

    all_names = features.gene.astype(str).tolist()
    source_holdouts = list(protocol["candidate_selection"]["positive_source_holdouts"])
    all_heldout, all_outer_validation = [], []
    candidate_columns: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for candidate, spec in protocol["feature_candidates"].items():
        config = candidate_config(base, spec["stage2_context_genes"])
        stage1_columns, stage2_columns = feature_columns(
            features, base, spec["stage2_context_genes"]
        )
        candidate_columns[candidate] = (stage1_columns, stage2_columns)
        heldout, outer_validation = run_candidate_holdouts(
            candidate, config, labels, matrix, all_names,
            stage1_columns, stage2_columns, partition, source_holdouts,
        )
        all_heldout.append(heldout)
        all_outer_validation.append(outer_validation)
    heldout = pd.concat(all_heldout, ignore_index=True)
    outer_validation = pd.concat(all_outer_validation, ignore_index=True)
    candidate_summary = summarize_candidates(heldout, protocol)
    selected_candidate = str(candidate_summary.loc[candidate_summary.selected, "candidate"].iloc[0])
    print(candidate_summary.to_string(index=False), flush=True)

    selected_spec = protocol["feature_candidates"][selected_candidate]
    selected_config = candidate_config(base, selected_spec["stage2_context_genes"])
    stage1_columns, stage2_columns = candidate_columns[selected_candidate]
    print(f"Fitting exact final development scorer: {selected_candidate}", flush=True)
    final_bundle = v44.train_bundle(
        labels, matrix, all_names, stage1_columns, stage2_columns, partition, selected_config
    )
    if len(final_bundle["stage2_effective_feature_names"]) > int(
        protocol["scientific_guardrails"]["maximum_effective_features"]
    ):
        raise RuntimeError("Selected V4.5 candidate exceeds the precommitted feature limit")
    summary = save_final(
        final_bundle, selected_candidate, candidate_summary, heldout, outer_validation,
        labels, partition, protocol, base, config_path, cache_manifest, started,
    )
    print(json.dumps(json_safe(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
