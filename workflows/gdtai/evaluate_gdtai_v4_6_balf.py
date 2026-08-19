#!/usr/bin/env python3
"""Score frozen gdTAI V4.6 on BALF and compare with historical V2/V3 calls."""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from evaluate_gdtai_v4_4_reused_lockbox import (
    calibration_from,
    load_model,
    overall_metrics,
    per_source_metrics,
    predict_booster,
)
from train_gdtai_v4_2_nested import exclusion_flags, sha256_file
from train_gdtai_v4_4_dual_mode import json_safe
from train_gdtai_v4_6_development import architecture_composer


ROOT = Path(__file__).resolve().parents[2]
MODEL_DIR = ROOT / "Integrated_dataset/models/gdT_prediction/gdtai_v4_6_development"
CACHE_DIR = ROOT / "Integrated_dataset/cache/gdT_prediction/gdtai_v4_3_common_lockbox"
INPUT = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_3_common_lockbox/common_lockbox_predictions.parquet"
INPUT_SUMMARY = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_3_common_lockbox/evaluation_summary.json"
OUTPUT_TABLES = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_6_balf"
OUTPUT_LOGS = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_6_balf"
BOOTSTRAP_ITERATIONS = 5_000
BOOTSTRAP_SEED = 20260820


def definitions() -> dict[str, tuple[str, str]]:
    return {
        "v4_6_highest_f1": ("v4_6_effective_score", "v4_6_highest_f1"),
        "v4_6_high_purity": ("v4_6_effective_score", "v4_6_high_purity"),
        "v3_balanced_historical": ("v3_balanced_score", "v3_balanced"),
        "v2_high_f1_historical": ("v2_score", "v2_high_f1"),
        "v2_high_purity_historical": ("v2_score", "v2_high_purity"),
    }


def load_inputs() -> tuple[pd.DataFrame, np.ndarray, list[str], dict, dict]:
    contract_path = MODEL_DIR / "model_contract.json"
    contract = json.loads(contract_path.read_text())
    if contract.get("iteration") != "V4.6" or contract.get("selected_candidate") != "conditional_gamma":
        raise RuntimeError("Unexpected V4.6 frozen model contract")
    if contract.get("model_promoted") or contract.get("new_untouched_test_scored"):
        raise RuntimeError("V4.6 contract already records test scoring or promotion")
    if contract.get("balf_used_for_fit_calibration_selection_or_threshold") is not False:
        raise RuntimeError("BALF isolation is not explicit in the V4.6 contract")

    original = json.loads(INPUT_SUMMARY.read_text())
    if sha256_file(INPUT) != original["predictions_sha256"]:
        raise RuntimeError("Frozen common-benchmark predictions changed")
    cache_summary = json.loads((CACHE_DIR / "lockbox_feature_cache_summary.json").read_text())
    matrix_path = CACHE_DIR / "lockbox_gene_features.npy"
    feature_path = CACHE_DIR / "lockbox_feature_manifest.csv"
    if sha256_file(matrix_path) != cache_summary["feature_cache_sha256"]:
        raise RuntimeError("Common feature cache checksum mismatch")
    if sha256_file(feature_path) != cache_summary["feature_manifest_sha256"]:
        raise RuntimeError("Common feature manifest checksum mismatch")
    frame = pd.read_parquet(INPUT).reset_index(drop=True)
    matrix = np.load(matrix_path, mmap_mode="r")
    names = pd.read_csv(feature_path).sort_values("feature_index").gene.astype(str).tolist()
    if len(frame) != len(matrix):
        raise RuntimeError("Common prediction/cache row mismatch")
    return frame, matrix, names, contract, cache_summary


def score_v46(frame: pd.DataFrame, matrix: np.ndarray, names: list[str], contract: dict) -> pd.DataFrame:
    lookup = {gene: index for index, gene in enumerate(names)}
    required = set(contract["stage1_feature_names"]) | set(contract["stage2_base_feature_names"])
    missing = sorted(required - set(lookup))
    if missing:
        raise RuntimeError(f"Common benchmark lacks required V4.6 features: {missing}")
    stage1_columns = np.asarray([lookup[gene] for gene in contract["stage1_feature_names"]])
    stage2_columns = np.asarray([lookup[gene] for gene in contract["stage2_base_feature_names"]])
    calibration = json.loads((MODEL_DIR / "platt_calibration.json").read_text())
    stage1_model = load_model(MODEL_DIR / "stage1_t_lineage_gate.ubj")
    stage2_model = load_model(MODEL_DIR / "stage2_gdt_classifier.ubj")
    stage1_probability = calibration_from(calibration["stage1"]).apply(
        predict_booster(stage1_model, np.asarray(matrix[:, stage1_columns]))
    )
    composer = architecture_composer("conditional_gamma")
    stage2_values, effective_names = composer(
        np.asarray(matrix[:, stage2_columns]), contract["stage2_base_feature_names"]
    )
    if effective_names != contract["stage2_effective_feature_names"]:
        raise RuntimeError("V4.6 effective feature order differs from the frozen contract")
    stage2_probability = calibration_from(calibration["stage2"]).apply(
        predict_booster(stage2_model, stage2_values)
    )
    excluded = exclusion_flags(np.asarray(matrix), names)[2]
    effective = stage2_probability.copy()
    effective[(stage1_probability < float(contract["stage1_threshold"])) | excluded] = -1.0
    output = frame.copy()
    output["v4_6_stage1_probability"] = stage1_probability
    output["v4_6_stage2_probability"] = stage2_probability
    output["v4_6_effective_score"] = effective
    output["v4_6_cd4_treg_excluded"] = excluded
    for mode, payload in contract["operating_modes"].items():
        output[f"v4_6_{mode}"] = effective >= float(payload["threshold"])
    return output


def truth_stratified_group_bootstrap(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    models = ["v4_6_highest_f1", "v3_balanced", "v2_high_f1"]
    local = frame[["truth_class", "source_gse_id", "donor_id", "group_id", *models]].copy()
    donor = local.donor_id.fillna("").astype(str)
    local["unit"] = np.where(
        donor.str.strip().ne(""), local.source_gse_id.astype(str) + "::" + donor,
        local.group_id.astype(str),
    )
    units = local.groupby(["truth_class", "unit"], as_index=False, observed=True).agg(
        n_cells=("unit", "size"), **{model: (model, "sum") for model in models}
    )
    strata = [group.reset_index(drop=True) for _, group in units.groupby("truth_class", sort=True)]
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    rows = []
    for iteration in range(BOOTSTRAP_ITERATIONS):
        sample = pd.concat(
            [stratum.iloc[rng.integers(0, len(stratum), size=len(stratum))] for stratum in strata],
            ignore_index=True,
        )
        counts = {
            truth: {
                "n": int(group.n_cells.sum()),
                **{model: int(group[model].sum()) for model in models},
            }
            for truth, group in sample.groupby("truth_class", observed=True)
        }
        metrics = {}
        for model in models:
            tp = counts["gdT_gold"][model]
            fn = counts["gdT_gold"]["n"] - tp
            fp = counts["abT_gold"][model] + counts["nk_lockbox"][model]
            precision = tp / (tp + fp) if tp + fp else 0.0
            recall = tp / (tp + fn) if tp + fn else 0.0
            metrics[model] = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        for comparator in ("v3_balanced", "v2_high_f1"):
            rows.append({
                "iteration": iteration,
                "comparator": comparator,
                "delta_f1": metrics["v4_6_highest_f1"] - metrics[comparator],
            })
    draws = pd.DataFrame(rows)
    summary = draws.groupby("comparator", observed=True).delta_f1.agg(
        median="median",
        ci_low=lambda values: values.quantile(0.025),
        ci_high=lambda values: values.quantile(0.975),
        probability_gt_zero=lambda values: (values > 0).mean(),
    ).reset_index()
    summary.insert(1, "iterations", BOOTSTRAP_ITERATIONS)
    summary["bootstrap_unit"] = "truth-stratified donor/group"
    return draws, summary


def paired_error_tests(frame: pd.DataFrame) -> pd.DataFrame:
    truth = frame.truth_class.eq("gdT_gold").to_numpy()
    v46 = frame.v4_6_highest_f1.to_numpy(bool) == truth
    rows = []
    for comparator, column in (("v3_balanced", "v3_balanced"), ("v2_high_f1", "v2_high_f1")):
        old = frame[column].to_numpy(bool) == truth
        v46_only = int((v46 & ~old).sum())
        old_only = int((~v46 & old).sum())
        discordant = v46_only + old_only
        rows.append({
            "comparator": comparator,
            "v4_6_correct_comparator_wrong": v46_only,
            "v4_6_wrong_comparator_correct": old_only,
            "discordant": discordant,
            "exact_mcnemar_p": float(binomtest(v46_only, discordant, 0.5).pvalue) if discordant else 1.0,
            "interpretation": "descriptive cell-level paired test; donor correlation ignored",
        })
    return pd.DataFrame(rows)


def main() -> None:
    started = time.monotonic()
    frame, matrix, names, contract, cache_summary = load_inputs()
    frame = score_v46(frame, matrix, names, contract)
    bcore = frame.source_gse_id.eq("BALF_BLOOD_COPD")
    if not bcore.any():
        raise RuntimeError("BALF is absent from the frozen common benchmark")
    bframe = frame.loc[bcore].reset_index(drop=True)
    model_definitions = definitions()
    overall = overall_metrics(bframe, model_definitions)
    overall.insert(0, "stratum", "BALF_untouched_for_V4_6")
    overall["comparison_status"] = np.where(
        overall.model.str.startswith("v4_6"),
        "single_untouched_V4_6_test",
        "historical_comparator_previously_consumed",
    )
    per_source = per_source_metrics(bframe, model_definitions)

    overlap_rows = []
    for v46_mode in ("v4_6_highest_f1", "v4_6_high_purity"):
        left = bframe[model_definitions[v46_mode][1]].to_numpy(bool)
        for comparator in ("v3_balanced_historical", "v2_high_f1_historical", "v2_high_purity_historical"):
            right = bframe[model_definitions[comparator][1]].to_numpy(bool)
            for truth, indexes in bframe.groupby("truth_class", observed=True).groups.items():
                idx = np.asarray(indexes, dtype=np.int64)
                a, b = left[idx], right[idx]
                overlap_rows.append({
                    "v4_6_mode": v46_mode,
                    "comparator": comparator,
                    "truth_class": truth,
                    "n_cells": len(idx),
                    "both_positive": int((a & b).sum()),
                    "v4_6_only": int((a & ~b).sum()),
                    "comparator_only": int((~a & b).sum()),
                    "both_negative": int((~a & ~b).sum()),
                })
    overlap = pd.DataFrame(overlap_rows)
    bootstrap_draws, bootstrap_summary = truth_stratified_group_bootstrap(bframe)
    mcnemar = paired_error_tests(bframe)

    OUTPUT_TABLES.mkdir(parents=True, exist_ok=True)
    OUTPUT_LOGS.mkdir(parents=True, exist_ok=True)
    predictions_path = OUTPUT_TABLES / "balf_predictions.parquet"
    bframe.to_parquet(predictions_path, index=False, compression="zstd")
    overall.to_csv(OUTPUT_TABLES / "overall_metrics.csv", index=False)
    per_source.to_csv(OUTPUT_TABLES / "per_source_metrics.csv", index=False)
    overlap.to_csv(OUTPUT_TABLES / "prediction_overlap.csv", index=False)
    bootstrap_draws.to_parquet(OUTPUT_TABLES / "balf_group_bootstrap_draws.parquet", index=False, compression="zstd")
    bootstrap_summary.to_csv(OUTPUT_TABLES / "balf_group_bootstrap_summary.csv", index=False)
    mcnemar.to_csv(OUTPUT_TABLES / "balf_paired_error_tests.csv", index=False)
    integrity = pd.DataFrame([
        {"check": "v4_6_contract_frozen_before_balf", "pass": True, "detail": sha256_file(MODEL_DIR / "model_contract.json")},
        {"check": "thresholds_retuned_on_balf", "pass": True, "detail": "false"},
        {"check": "features_or_models_changed_after_balf", "pass": True, "detail": "false"},
        {"check": "historical_comparators_marked_consumed", "pass": True, "detail": "true"},
        {"check": "eligible_for_promotion", "pass": True, "detail": "false; no untouched confirmatory test remains"},
    ])
    integrity.to_csv(OUTPUT_TABLES / "implementation_integrity.csv", index=False)
    summary = {
        "status": "PASS_V4_6_BALF_SINGLE_UNTOUCHED_DIAGNOSTIC_NOT_PROMOTABLE",
        "runtime_seconds": time.monotonic() - started,
        "n_balf_cells": len(bframe),
        "truth_counts": bframe.truth_class.value_counts().to_dict(),
        "v4_6_model_contract_sha256_before_scoring": sha256_file(MODEL_DIR / "model_contract.json"),
        "common_cache_sha256": cache_summary["feature_cache_sha256"],
        "output_predictions_sha256": sha256_file(predictions_path),
        "thresholds_retuned": False,
        "features_or_models_changed": False,
        "historical_comparators_previously_consumed": True,
        "model_promoted": False,
        "overall_metrics": overall.to_dict("records"),
        "group_bootstrap": bootstrap_summary.to_dict("records"),
        "paired_error_tests": mcnemar.to_dict("records"),
    }
    (OUTPUT_LOGS / "evaluation_summary.json").write_text(
        json.dumps(json_safe(summary), indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(json_safe(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
