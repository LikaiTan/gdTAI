#!/usr/bin/env python3
"""Compare frozen gdTAI V4.5 with V2/V3 on the consumed common benchmark."""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from evaluate_gdtai_v4_4_reused_lockbox import (
    calibration_from,
    load_model,
    overall_metrics,
    overlap_table,
    per_source_metrics,
    predict_booster,
)
from train_gdtai_v4_2_nested import exclusion_flags, sha256_file
from train_gdtai_v4_4_dual_mode import compose_stage2, json_safe


ROOT = Path(__file__).resolve().parents[2]
MODEL_DIR = ROOT / "Integrated_dataset/models/gdT_prediction/gdtai_v4_5_development"
LOCKBOX_DIR = ROOT / "Integrated_dataset/cache/gdT_prediction/gdtai_v4_3_common_lockbox"
INPUT_PREDICTIONS = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_3_common_lockbox/common_lockbox_predictions.parquet"
INPUT_SUMMARY = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_3_common_lockbox/evaluation_summary.json"
OUTPUT_TABLES = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_5_vs_v2_v3_consumed"
OUTPUT_LOGS = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_5_vs_v2_v3_consumed"
TRAINING_EXPOSED_SOURCES = {"GDT_2020AUG_woCOV", "GDTlung2023july_7p"}


def definitions() -> dict[str, tuple[str, str]]:
    return {
        "v4_5_highest_f1": ("v4_5_effective_score", "v4_5_highest_f1"),
        "v4_5_high_purity": ("v4_5_effective_score", "v4_5_high_purity"),
        "v3_balanced": ("v3_balanced_score", "v3_balanced"),
        "v2_high_f1": ("v2_score", "v2_high_f1"),
        "v2_high_purity": ("v2_score", "v2_high_purity"),
    }


def load_and_verify() -> tuple[pd.DataFrame, np.ndarray, list[str], dict, dict]:
    contract_path = MODEL_DIR / "model_contract.json"
    contract = json.loads(contract_path.read_text())
    if contract["status"] != "DEVELOPMENT_FROZEN_NOT_PROMOTABLE_NO_NEW_TEST":
        raise RuntimeError("V4.5 contract is not development-frozen")
    if contract["model_promoted"] or contract["new_untouched_test_scored"]:
        raise RuntimeError("V4.5 contract unexpectedly records promotion or new-test scoring")
    if contract["selected_candidate"] != "receptor_context_baseline":
        raise RuntimeError("Unexpected V4.5 selected candidate")

    original = json.loads(INPUT_SUMMARY.read_text())
    if sha256_file(INPUT_PREDICTIONS) != original["predictions_sha256"]:
        raise RuntimeError("Frozen V2/V3 common-benchmark predictions changed")
    cache_summary = json.loads((LOCKBOX_DIR / "lockbox_feature_cache_summary.json").read_text())
    cache_path = LOCKBOX_DIR / "lockbox_gene_features.npy"
    manifest_path = LOCKBOX_DIR / "lockbox_feature_manifest.csv"
    if sha256_file(cache_path) != cache_summary["feature_cache_sha256"]:
        raise RuntimeError("Common-benchmark feature cache checksum mismatch")
    if sha256_file(manifest_path) != cache_summary["feature_manifest_sha256"]:
        raise RuntimeError("Common-benchmark feature manifest checksum mismatch")

    frame = pd.read_parquet(INPUT_PREDICTIONS).reset_index(drop=True)
    matrix = np.load(cache_path, mmap_mode="r")
    names = pd.read_csv(manifest_path).sort_values("feature_index").gene.astype(str).tolist()
    if len(frame) != len(matrix) or len(frame) != int(cache_summary["n_cells"]):
        raise RuntimeError("Common-benchmark row count mismatch")
    if set(frame.source_gse_id[frame.truth_class.eq("gdT_gold")].unique()) != {
        "BALF_BLOOD_COPD", "GDT_2020AUG_woCOV", "GDTlung2023july_7p"
    }:
        raise RuntimeError("Positive common-benchmark sources changed")
    return frame, matrix, names, contract, cache_summary


def score_v45(frame: pd.DataFrame, matrix: np.ndarray, names: list[str], contract: dict) -> pd.DataFrame:
    lookup = {gene: index for index, gene in enumerate(names)}
    required = set(contract["stage1_feature_names"]) | set(contract["stage2_base_feature_names"])
    missing = sorted(required - set(lookup))
    if missing:
        raise RuntimeError(f"Common benchmark lacks V4.5 features: {missing}")
    stage1_columns = np.asarray([lookup[gene] for gene in contract["stage1_feature_names"]], dtype=np.int64)
    stage2_columns = np.asarray([lookup[gene] for gene in contract["stage2_base_feature_names"]], dtype=np.int64)
    calibration = json.loads((MODEL_DIR / "platt_calibration.json").read_text())

    stage1_model = load_model(MODEL_DIR / "stage1_t_lineage_gate.ubj")
    stage2_model = load_model(MODEL_DIR / "stage2_gdt_classifier.ubj")
    stage1_probability = calibration_from(calibration["stage1"]).apply(
        predict_booster(stage1_model, np.asarray(matrix[:, stage1_columns]))
    )
    stage2_values, effective_names = compose_stage2(
        np.asarray(matrix[:, stage2_columns]), contract["stage2_base_feature_names"]
    )
    if effective_names != contract["stage2_effective_feature_names"]:
        raise RuntimeError("V4.5 engineered feature order differs from its frozen contract")
    stage2_probability = calibration_from(calibration["stage2"]).apply(
        predict_booster(stage2_model, stage2_values)
    )
    excluded = exclusion_flags(np.asarray(matrix), names)[2]
    effective_score = stage2_probability.copy()
    effective_score[(stage1_probability < float(contract["stage1_threshold"])) | excluded] = -1.0

    output = frame.copy()
    output["v4_5_stage1_probability"] = stage1_probability
    output["v4_5_stage2_probability"] = stage2_probability
    output["v4_5_effective_score"] = effective_score
    output["v4_5_cd4_treg_excluded"] = excluded
    for mode, payload in contract["operating_modes"].items():
        output[f"v4_5_{mode}"] = effective_score >= float(payload["threshold"])
    output["v4_5_exposure"] = np.where(
        output.source_gse_id.isin(TRAINING_EXPOSED_SOURCES),
        "training_exposed_consumed", "unexposed_but_previously_consumed",
    )
    return output


def stratified_metrics(frame: pd.DataFrame, model_definitions: dict[str, tuple[str, str]]) -> pd.DataFrame:
    strata = {
        "all_consumed_benchmark": np.ones(len(frame), dtype=bool),
        "v4_5_unexposed_consumed": ~frame.source_gse_id.isin(TRAINING_EXPOSED_SOURCES).to_numpy(),
        "balf_unexposed_consumed": frame.source_gse_id.eq("BALF_BLOOD_COPD").to_numpy(),
        "v4_5_training_exposed_sorted": frame.source_gse_id.isin(TRAINING_EXPOSED_SOURCES).to_numpy(),
    }
    tables = []
    for stratum, mask in strata.items():
        local = overall_metrics(frame.loc[mask].reset_index(drop=True), model_definitions)
        local.insert(0, "stratum", stratum)
        local["v4_5_external_evidence"] = False
        local["interpretation"] = (
            "training-exposed descriptive result" if stratum == "v4_5_training_exposed_sorted"
            else "post-hoc consumed diagnostic"
        )
        tables.append(local)
    return pd.concat(tables, ignore_index=True)


def source_overlap(frame: pd.DataFrame, model_definitions: dict[str, tuple[str, str]]) -> pd.DataFrame:
    comparators = ["v3_balanced", "v2_high_f1", "v2_high_purity"]
    output = []
    for mode in ("v4_5_highest_f1", "v4_5_high_purity"):
        left = frame[model_definitions[mode][1]].to_numpy(bool)
        for comparator in comparators:
            right = frame[model_definitions[comparator][1]].to_numpy(bool)
            for (source, truth), indexes in frame.groupby(
                ["source_gse_id", "truth_class"], observed=True
            ).groups.items():
                idx = np.asarray(indexes, dtype=np.int64)
                a, b = left[idx], right[idx]
                union = int((a | b).sum())
                output.append({
                    "v4_5_mode": mode,
                    "comparator": comparator,
                    "source_gse_id": source,
                    "truth_class": truth,
                    "v4_5_exposure": frame.iloc[idx[0]].v4_5_exposure,
                    "n_cells": int(len(idx)),
                    "both_positive": int((a & b).sum()),
                    "v4_5_only": int((a & ~b).sum()),
                    "comparator_only": int((~a & b).sum()),
                    "both_negative": int((~a & ~b).sum()),
                    "jaccard_positive": float((a & b).sum() / union) if union else np.nan,
                })
    return pd.DataFrame(output)


def main() -> None:
    started = time.monotonic()
    frame, matrix, names, contract, cache_summary = load_and_verify()
    frame = score_v45(frame, matrix, names, contract)
    model_definitions = definitions()
    overall = stratified_metrics(frame, model_definitions)
    per_source = per_source_metrics(frame, model_definitions)
    exposure = frame[["source_gse_id", "v4_5_exposure"]].drop_duplicates()
    per_source = per_source.merge(exposure, on="source_gse_id", how="left", validate="many_to_one")
    overlap = overlap_table(frame, model_definitions)
    detailed_overlap = source_overlap(frame, model_definitions)

    OUTPUT_TABLES.mkdir(parents=True, exist_ok=True)
    OUTPUT_LOGS.mkdir(parents=True, exist_ok=True)
    predictions_path = OUTPUT_TABLES / "consumed_benchmark_predictions.parquet"
    frame.to_parquet(predictions_path, index=False, compression="zstd")
    overall.to_csv(OUTPUT_TABLES / "stratified_overall_metrics.csv", index=False)
    per_source.to_csv(OUTPUT_TABLES / "per_source_metrics.csv", index=False)
    overlap.to_csv(OUTPUT_TABLES / "overall_prediction_overlap.csv", index=False)
    detailed_overlap.to_csv(OUTPUT_TABLES / "per_source_truth_prediction_overlap.csv", index=False)
    exposure_counts = frame.groupby(
        ["v4_5_exposure", "source_gse_id", "truth_class"], observed=True
    ).size().reset_index(name="n_cells")
    exposure_counts.to_csv(OUTPUT_TABLES / "benchmark_exposure_manifest.csv", index=False)

    integrity = pd.DataFrame([
        {"check": "v4_5_thresholds_unchanged", "pass": True,
         "detail": json.dumps({key: value["threshold"] for key, value in contract["operating_modes"].items()}, sort_keys=True)},
        {"check": "v2_v3_predictions_checksum", "pass": True, "detail": sha256_file(INPUT_PREDICTIONS)},
        {"check": "common_feature_cache_checksum", "pass": True, "detail": cache_summary["feature_cache_sha256"]},
        {"check": "training_exposed_sources_labeled", "pass": True, "detail": ",".join(sorted(TRAINING_EXPOSED_SOURCES))},
        {"check": "eligible_for_promotion", "pass": True, "detail": "false; consumed diagnostic only"},
    ])
    integrity.to_csv(OUTPUT_TABLES / "implementation_integrity.csv", index=False)
    summary = {
        "status": "PASS_CONSUMED_DIAGNOSTIC_ONLY",
        "runtime_seconds": time.monotonic() - started,
        "n_cells": len(frame),
        "v4_5_model_contract_sha256": sha256_file(MODEL_DIR / "model_contract.json"),
        "input_predictions_sha256": sha256_file(INPUT_PREDICTIONS),
        "output_predictions_sha256": sha256_file(predictions_path),
        "thresholds_retuned": False,
        "features_or_models_changed": False,
        "training_exposed_sources": sorted(TRAINING_EXPOSED_SOURCES),
        "new_untouched_test_scored": False,
        "eligible_for_promotion": False,
        "model_promoted": False,
        "stratified_metrics": overall.to_dict("records"),
    }
    (OUTPUT_LOGS / "evaluation_summary.json").write_text(
        json.dumps(json_safe(summary), indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(json_safe(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
