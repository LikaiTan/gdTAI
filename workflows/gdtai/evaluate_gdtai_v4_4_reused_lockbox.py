#!/usr/bin/env python3
"""Evaluate frozen gdTAI V4.4 modes on the consumed V4.3 lockbox."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb

from train_gdtai_v4_2_nested import exclusion_flags, sha256_file
from train_gdtai_v4_4_dual_mode import (
    PlattCalibration,
    binary_metrics,
    compose_stage2,
    json_safe,
    resolve,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MODEL_DIR = ROOT / "Integrated_dataset/models/gdT_prediction/gdtai_v4_4_dual_mode"
LOCKBOX_DIR = ROOT / "Integrated_dataset/cache/gdT_prediction/gdtai_v4_3_common_lockbox"
INPUT_PREDICTIONS = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_3_common_lockbox/common_lockbox_predictions.parquet"
OUTPUT_TABLE_DIR = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_4_reused_lockbox"
OUTPUT_LOG_DIR = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_4_reused_lockbox"


def load_model(path: Path) -> xgb.Booster:
    model = xgb.Booster()
    model.load_model(path)
    model.set_param({"device": "cpu", "nthread": 16})
    return model


def predict_booster(model: xgb.Booster, values: np.ndarray, chunk: int = 100_000) -> np.ndarray:
    output = []
    for start in range(0, len(values), chunk):
        block = xgb.DMatrix(np.asarray(values[start:start + chunk], dtype=np.float32))
        output.append(model.predict(block))
    return np.concatenate(output).astype(np.float32) if output else np.empty(0, dtype=np.float32)


def calibration_from(payload: dict) -> PlattCalibration:
    return PlattCalibration(float(payload["coefficient"]), float(payload["intercept"]))


def prevalence_precision(recall: float, fpr: float, prevalence: float) -> float:
    numerator = recall * prevalence
    denominator = numerator + fpr * (1 - prevalence)
    return numerator / denominator if denominator else 0.0


def overall_metrics(frame: pd.DataFrame, definitions: dict[str, tuple[str, str]]) -> pd.DataFrame:
    primary = frame.truth_class.isin(["gdT_gold", "abT_gold", "nk_lockbox"]).to_numpy()
    y = frame.truth_class.eq("gdT_gold").to_numpy(np.int8)
    ab = frame.truth_class.eq("abT_gold").to_numpy()
    nk = frame.truth_class.eq("nk_lockbox").to_numpy()
    output = []
    for model_name, (score_column, call_column) in definitions.items():
        score = frame[score_column].to_numpy(np.float32)
        calls = frame[call_column].to_numpy(bool)
        metrics = binary_metrics(y[primary], score[primary], calls[primary])
        metrics.update({
            "model": model_name,
            "predicted_positive": int(calls.sum()),
            "abt_fpr": float(calls[ab].mean()),
            "author_nk_fpr": float(calls[nk].mean()),
        })
        total_negative_fpr = float(calls[ab | nk].mean())
        for prevalence in (0.001, 0.005, 0.01, 0.03, 0.05):
            metrics[f"precision_at_{prevalence:g}_prevalence"] = prevalence_precision(
                float(metrics["recall"]), total_negative_fpr, prevalence
            )
        output.append(metrics)
    return pd.DataFrame(output)


def per_source_metrics(frame: pd.DataFrame, definitions: dict[str, tuple[str, str]]) -> pd.DataFrame:
    output = []
    for model_name, (_, call_column) in definitions.items():
        calls = frame[call_column].to_numpy(bool)
        for source, indexes in frame.groupby("source_gse_id", observed=True).groups.items():
            idx = np.asarray(indexes, dtype=np.int64)
            truth = frame.iloc[idx].truth_class.to_numpy()
            local_calls = calls[idx]
            gd = truth == "gdT_gold"
            ab = truth == "abT_gold"
            nk = truth == "nk_lockbox"
            output.append({
                "model": model_name,
                "source_gse_id": source,
                "n_cells": int(len(idx)),
                "n_gdt_gold": int(gd.sum()),
                "n_abt_gold": int(ab.sum()),
                "n_author_nk": int(nk.sum()),
                "predicted_positive": int(local_calls.sum()),
                "gdt_recall": float(local_calls[gd].mean()) if gd.any() else np.nan,
                "abt_fpr": float(local_calls[ab].mean()) if ab.any() else np.nan,
                "author_nk_fpr": float(local_calls[nk].mean()) if nk.any() else np.nan,
            })
    return pd.DataFrame(output)


def overlap_table(frame: pd.DataFrame, definitions: dict[str, tuple[str, str]]) -> pd.DataFrame:
    output = []
    names = list(definitions)
    for left_index, left in enumerate(names):
        left_calls = frame[definitions[left][1]].to_numpy(bool)
        for right in names[left_index + 1:]:
            right_calls = frame[definitions[right][1]].to_numpy(bool)
            output.append({
                "model_a": left, "model_b": right,
                "both_positive": int((left_calls & right_calls).sum()),
                "a_only": int((left_calls & ~right_calls).sum()),
                "b_only": int((~left_calls & right_calls).sum()),
                "both_negative": int((~left_calls & ~right_calls).sum()),
                "jaccard_positive": float((left_calls & right_calls).sum() / max(1, (left_calls | right_calls).sum())),
            })
    return pd.DataFrame(output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-dir", type=Path, default=DEFAULT_MODEL_DIR)
    args = parser.parse_args()
    started = time.monotonic()
    model_dir = resolve(args.model_dir)
    contract_path = model_dir / "model_contract.json"
    contract = json.loads(contract_path.read_text())
    if contract["thresholds_selected_from_test"] is not False:
        raise RuntimeError("V4.4 thresholds were not development-only")
    if contract["test_scored"] is not False:
        raise RuntimeError("Expected an unscored frozen development contract")
    if contract["test_policy"]["consumed_v4_3_lockbox_is_diagnostic_only"] is not True:
        raise RuntimeError("The frozen contract does not mark this lockbox diagnostic-only")

    feature_manifest_path = LOCKBOX_DIR / "lockbox_feature_manifest.csv"
    cache_path = LOCKBOX_DIR / "lockbox_gene_features.npy"
    cache_summary = json.loads((LOCKBOX_DIR / "lockbox_feature_cache_summary.json").read_text())
    if sha256_file(cache_path) != cache_summary["feature_cache_sha256"]:
        raise RuntimeError("Common-lockbox feature cache checksum mismatch")
    if sha256_file(feature_manifest_path) != cache_summary["feature_manifest_sha256"]:
        raise RuntimeError("Common-lockbox feature manifest checksum mismatch")
    feature_manifest = pd.read_csv(feature_manifest_path).sort_values("feature_index")
    lockbox_names = feature_manifest.gene.astype(str).tolist()
    lookup = {gene: index for index, gene in enumerate(lockbox_names)}
    required = set(contract["stage1_feature_names"]) | set(contract["stage2_base_feature_names"])
    missing = sorted(required - set(lookup))
    if missing:
        raise RuntimeError(f"Common lockbox lacks V4.4 features: {missing}")

    frame = pd.read_parquet(INPUT_PREDICTIONS).reset_index(drop=True)
    matrix = np.load(cache_path, mmap_mode="r")
    if len(frame) != len(matrix) or len(frame) != int(cache_summary["n_cells"]):
        raise RuntimeError("Common-lockbox row count mismatch")
    stage1_columns = np.asarray([lookup[gene] for gene in contract["stage1_feature_names"]], dtype=np.int64)
    stage2_columns = np.asarray([lookup[gene] for gene in contract["stage2_base_feature_names"]], dtype=np.int64)
    calibration = json.loads((model_dir / "platt_calibration.json").read_text())

    stage1_model = load_model(model_dir / "stage1_t_lineage_gate.ubj")
    stage2_model = load_model(model_dir / "stage2_receptor_classifier.ubj")
    stage1_probability = calibration_from(calibration["stage1"]).apply(
        predict_booster(stage1_model, np.asarray(matrix[:, stage1_columns]))
    )
    stage2_values, effective_names = compose_stage2(
        np.asarray(matrix[:, stage2_columns]), contract["stage2_base_feature_names"]
    )
    if effective_names != contract["stage2_effective_feature_names"]:
        raise RuntimeError("Stage-2 engineered feature order differs from the frozen contract")
    stage2_probability = calibration_from(calibration["stage2"]).apply(
        predict_booster(stage2_model, stage2_values)
    )
    excluded = exclusion_flags(np.asarray(matrix), lockbox_names)[2]
    effective_score = stage2_probability.copy()
    effective_score[(stage1_probability < float(contract["stage1_threshold"])) | excluded] = -1.0
    frame["v4_4_stage1_probability"] = stage1_probability
    frame["v4_4_stage2_probability"] = stage2_probability
    frame["v4_4_effective_score"] = effective_score
    frame["v4_4_cd4_treg_excluded"] = excluded
    for mode, payload in contract["operating_modes"].items():
        frame[f"v4_4_{mode}"] = effective_score >= float(payload["threshold"])

    definitions = {
        "v4_4_highest_f1": ("v4_4_effective_score", "v4_4_highest_f1"),
        "v4_4_high_purity": ("v4_4_effective_score", "v4_4_high_purity"),
        "v4_3_balanced": ("v4_3_balanced_score", "v4_3_balanced"),
        "v3_balanced": ("v3_balanced_score", "v3_balanced"),
        "v2_high_f1": ("v2_score", "v2_high_f1"),
        "v2_high_purity": ("v2_score", "v2_high_purity"),
    }
    overall = overall_metrics(frame, definitions)
    per_source = per_source_metrics(frame, definitions)
    overlap = overlap_table(frame, definitions)

    OUTPUT_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_LOG_DIR.mkdir(parents=True, exist_ok=True)
    predictions_path = OUTPUT_TABLE_DIR / "reused_lockbox_predictions.parquet"
    frame.to_parquet(predictions_path, index=False, compression="zstd")
    overall.to_csv(OUTPUT_TABLE_DIR / "overall_metrics.csv", index=False)
    per_source.to_csv(OUTPUT_TABLE_DIR / "per_source_metrics.csv", index=False)
    overlap.to_csv(OUTPUT_TABLE_DIR / "prediction_overlap.csv", index=False)
    integrity = pd.DataFrame([
        {"check": "thresholds_selected_from_test", "pass": not contract["thresholds_selected_from_test"], "detail": "false"},
        {"check": "lockbox_is_diagnostic_only", "pass": contract["test_policy"]["consumed_v4_3_lockbox_is_diagnostic_only"], "detail": "consumed V4.3 common lockbox"},
        {"check": "feature_cache_checksum", "pass": True, "detail": cache_summary["feature_cache_sha256"]},
        {"check": "thresholds_unchanged", "pass": True, "detail": json.dumps({key: value["threshold"] for key, value in contract["operating_modes"].items()}, sort_keys=True)},
    ])
    integrity.to_csv(OUTPUT_TABLE_DIR / "implementation_integrity.csv", index=False)
    summary = {
        "status": "PASS_REUSED_LOCKBOX_DIAGNOSTIC",
        "runtime_seconds": time.monotonic() - started,
        "n_cells": len(frame),
        "lockbox_sha256": cache_summary["lockbox_sha256"],
        "model_contract_sha256": sha256_file(contract_path),
        "predictions_sha256": sha256_file(predictions_path),
        "thresholds_retuned": False,
        "thresholds_selected_from_test": False,
        "lockbox_reused_after_v4_3_consumption": True,
        "eligible_for_promotion_review": False,
        "model_promoted": False,
        "overall_metrics": overall.to_dict("records"),
    }
    (OUTPUT_LOG_DIR / "evaluation_summary.json").write_text(
        json.dumps(json_safe(summary), indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(json_safe(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
