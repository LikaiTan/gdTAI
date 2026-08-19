#!/usr/bin/env python3
"""Freeze the gdTAI V4.3 model and thresholds using development cells only."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from train_gdtai_v4_2_nested import (
    audit_rows,
    candidate_id,
    choose_profile_threshold,
    choose_stage1_threshold,
    compose_stage2_features,
    exclusion_flags,
    fit_xgb,
    predict,
    profile_metrics,
    receptor_rescue_flags,
    resolve,
    sampled_rows,
    sha256_file,
    source_balanced_weights,
    stable_seed,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "configs/models/gdtai/v4_3_rescue_training.json"
FIT_CLASSES = ["gdT_gold", "abT_gold", "single_ab_support", "nk_tier1", "nk_tier2"]
STAGE2_CLASSES = ["gdT_gold", "abT_gold", "single_ab_support"]


def effective_columns(features: pd.DataFrame, config: dict) -> tuple[np.ndarray, np.ndarray, bool]:
    stage1_excluded = set(config.get("stage1_exclude_genes", []))
    stage1 = np.flatnonzero(
        features.stage1.to_numpy(bool)
        & ~features.gene.astype(str).isin(stage1_excluded).to_numpy()
    )
    policy = config.get("stage2_feature_policy", {})
    if policy:
        controls = set(policy.get("control_genes", []))
        stage2 = np.flatnonzero(
            features.feature_class.eq("TCR").to_numpy()
            | features.gene.astype(str).isin(controls).to_numpy()
        )
    else:
        stage2 = np.flatnonzero(features.stage2.to_numpy(bool))
    return stage1, stage2, bool(config.get("stage2_include_stage1_probability", True))


def fit_stage1_oof(
    labels: pd.DataFrame,
    audit: np.ndarray,
    x: np.ndarray,
    columns: np.ndarray,
    params: dict,
    config: dict,
) -> tuple[np.ndarray, float, dict]:
    frame = labels.iloc[audit].reset_index(drop=True)
    folds = sorted(frame.inner_fold.dropna().astype(int).unique())
    oof = np.full(len(audit), np.nan, dtype=np.float32)
    model_id = candidate_id(params)
    for fold in folds:
        fit_local = np.flatnonzero(
            frame.inner_fold.ne(fold).to_numpy()
            & frame.truth_class.isin(FIT_CLASSES).to_numpy()
        )
        fit_frame = frame.iloc[fit_local]
        t_local = np.flatnonzero(~fit_frame.truth_class.str.startswith("nk_").to_numpy())
        nk_local = np.flatnonzero(fit_frame.truth_class.str.startswith("nk_").to_numpy())
        t_selected = sampled_rows(
            fit_frame.reset_index(drop=True),
            t_local,
            config["maximum_stage1_t_cells"],
            stable_seed("v4.3-final", model_id, fold, "stage1-t"),
        )
        train_local = fit_local[np.concatenate([t_selected, nk_local])]
        train_rows = audit[train_local]
        validation_local = np.flatnonzero(frame.inner_fold.eq(fold).to_numpy())
        validation_rows = audit[validation_local]
        model = fit_xgb(
            np.asarray(x[train_rows][:, columns]),
            (~labels.iloc[train_rows].truth_class.str.startswith("nk_")).to_numpy(np.int8),
            source_balanced_weights(labels, train_rows, 1),
            params,
            config["fixed"],
            stable_seed("v4.3-final", model_id, fold, "stage1"),
            (
                np.asarray(x[validation_rows][:, columns]),
                (~frame.iloc[validation_local].truth_class.str.startswith("nk_")).to_numpy(np.int8),
            ),
        )
        oof[validation_local] = predict(model, np.asarray(x[validation_rows][:, columns]))
    if not np.isfinite(oof).all():
        raise RuntimeError("Stage 1 grouped OOF predictions are incomplete")
    threshold, details = choose_stage1_threshold(frame, oof, config["stage1"])
    if details.get("eligible", True) is not True:
        raise RuntimeError("Final Stage 1 candidate failed the frozen development gate")
    return oof, float(threshold), details


def fit_stage2_oof(
    labels: pd.DataFrame,
    audit: np.ndarray,
    x: np.ndarray,
    columns: np.ndarray,
    stage1_oof: np.ndarray,
    include_stage1_probability: bool,
    params: dict,
    config: dict,
) -> np.ndarray:
    frame = labels.iloc[audit].reset_index(drop=True)
    folds = sorted(frame.inner_fold.dropna().astype(int).unique())
    oof = np.full(len(audit), np.nan, dtype=np.float32)
    model_id = candidate_id(params)
    for fold in folds:
        fit_local = np.flatnonzero(
            frame.inner_fold.ne(fold).to_numpy()
            & frame.truth_class.isin(STAGE2_CLASSES).to_numpy()
        )
        fit_frame = frame.iloc[fit_local]
        positive = np.flatnonzero(fit_frame.truth_class.eq("gdT_gold").to_numpy())
        negative = np.flatnonzero(~fit_frame.truth_class.eq("gdT_gold").to_numpy())
        negative_selected = sampled_rows(
            fit_frame.reset_index(drop=True),
            negative,
            config["maximum_stage2_negative_cells"],
            stable_seed("v4.3-final", model_id, fold, "stage2-negative"),
        )
        train_local = fit_local[np.concatenate([positive, negative_selected])]
        train_rows = audit[train_local]
        validation_local = np.flatnonzero(frame.inner_fold.eq(fold).to_numpy())
        validation_rows = audit[validation_local]
        validation_t = validation_local[frame.iloc[validation_local].truth_class.isin(STAGE2_CLASSES).to_numpy()]
        model = fit_xgb(
            compose_stage2_features(x, train_rows, columns, stage1_oof[train_local], include_stage1_probability),
            labels.iloc[train_rows].truth_class.eq("gdT_gold").to_numpy(np.int8),
            source_balanced_weights(labels, train_rows, 2),
            params,
            config["fixed"],
            stable_seed("v4.3-final", model_id, fold, "stage2"),
            (
                compose_stage2_features(
                    x,
                    audit[validation_t],
                    columns,
                    stage1_oof[validation_t],
                    include_stage1_probability,
                ),
                frame.iloc[validation_t].truth_class.eq("gdT_gold").to_numpy(np.int8),
            ),
        )
        oof[validation_local] = predict(
            model,
            compose_stage2_features(
                x,
                validation_rows,
                columns,
                stage1_oof[validation_local],
                include_stage1_probability,
            ),
        )
    if not np.isfinite(oof).all():
        raise RuntimeError("Stage 2 grouped OOF predictions are incomplete")
    return oof


def fit_final_models(
    labels: pd.DataFrame,
    x: np.ndarray,
    stage1_columns: np.ndarray,
    stage2_columns: np.ndarray,
    include_stage1_probability: bool,
    params: dict,
    config: dict,
    model_dir: Path,
) -> tuple[Path, Path, int, int]:
    development = np.flatnonzero(labels.allow_fit.to_numpy() & labels.truth_class.isin(FIT_CLASSES).to_numpy())
    frame = labels.iloc[development]
    t_local = np.flatnonzero(~frame.truth_class.str.startswith("nk_").to_numpy())
    nk_local = np.flatnonzero(frame.truth_class.str.startswith("nk_").to_numpy())
    t_selected = sampled_rows(
        frame.reset_index(drop=True),
        t_local,
        config["maximum_stage1_t_cells"],
        stable_seed("v4.3-final-fit", "stage1-t"),
    )
    stage1_rows = development[np.concatenate([t_selected, nk_local])]
    stage1_model = fit_xgb(
        np.asarray(x[stage1_rows][:, stage1_columns]),
        (~labels.iloc[stage1_rows].truth_class.str.startswith("nk_")).to_numpy(np.int8),
        source_balanced_weights(labels, stage1_rows, 1),
        params,
        config["fixed"],
        stable_seed("v4.3-final-fit", "stage1"),
        None,
    )

    stage2_pool = np.flatnonzero(labels.allow_fit.to_numpy() & labels.truth_class.isin(STAGE2_CLASSES).to_numpy())
    stage2_frame = labels.iloc[stage2_pool]
    positive = np.flatnonzero(stage2_frame.truth_class.eq("gdT_gold").to_numpy())
    negative = np.flatnonzero(~stage2_frame.truth_class.eq("gdT_gold").to_numpy())
    negative_selected = sampled_rows(
        stage2_frame.reset_index(drop=True),
        negative,
        config["maximum_stage2_negative_cells"],
        stable_seed("v4.3-final-fit", "stage2-negative"),
    )
    stage2_rows = stage2_pool[np.concatenate([positive, negative_selected])]
    stage2_model = fit_xgb(
        compose_stage2_features(
            x,
            stage2_rows,
            stage2_columns,
            np.zeros(len(stage2_rows), dtype=np.float32),
            include_stage1_probability,
        ),
        labels.iloc[stage2_rows].truth_class.eq("gdT_gold").to_numpy(np.int8),
        source_balanced_weights(labels, stage2_rows, 2),
        params,
        config["fixed"],
        stable_seed("v4.3-final-fit", "stage2"),
        None,
    )
    model_dir.mkdir(parents=True, exist_ok=True)
    stage1_path = model_dir / "stage1_t_lineage_gate.ubj"
    stage2_path = model_dir / "stage2_receptor_classifier.ubj"
    stage1_model.save_model(stage1_path)
    stage2_model.save_model(stage2_path)
    return stage1_path, stage2_path, len(stage1_rows), len(stage2_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    args = parser.parse_args()
    started = time.monotonic()
    config_path = resolve(args.config)
    config = json.loads(config_path.read_text())
    table_dir = resolve(config["output_table_dir"]) / "final_development"
    model_dir = resolve(config["output_model_dir"]) / "final_development"
    log_dir = resolve(config["output_log_dir"]) / "final_development"
    for path in (table_dir, model_dir, log_dir):
        path.mkdir(parents=True, exist_ok=True)

    labels_path = resolve(config["label_manifest"])
    features_path = resolve(config["feature_manifest"])
    matrix_path = resolve(config["cache_matrix"])
    split_path = resolve(config["split_manifest"])
    cache_manifest_path = resolve(config["cache_manifest"])
    cache_manifest = json.loads(cache_manifest_path.read_text())
    if cache_manifest["label_manifest_sha256"] != sha256_file(labels_path):
        raise RuntimeError("Label manifest differs from the frozen feature cache")
    if cache_manifest["feature_manifest_sha256"] != sha256_file(features_path):
        raise RuntimeError("Feature manifest differs from the frozen feature cache")

    labels = pd.read_parquet(labels_path).reset_index(drop=True)
    split = pd.read_csv(split_path)
    group_to_fold = split.drop_duplicates("group_id").set_index("group_id")["inner_fold"]
    labels["inner_fold"] = labels.group_id.map(group_to_fold)
    if labels.loc[labels.allow_fit, "inner_fold"].isna().any():
        raise RuntimeError("Development fold assignments are incomplete")
    features = pd.read_csv(features_path).sort_values("feature_index").reset_index(drop=True)
    names = features.gene.astype(str).tolist()
    stage1_columns, stage2_columns, include_stage1_probability = effective_columns(features, config)
    if include_stage1_probability:
        raise RuntimeError("V4.3 receptor-isolated protocol must not propagate Stage 1 probability")
    banned = set(config.get("stage1_exclude_genes", []))
    if banned & set(features.iloc[stage1_columns].gene.astype(str)):
        raise RuntimeError("A banned shared cytotoxic/NK gene entered Stage 1")
    if banned & set(features.iloc[stage2_columns].gene.astype(str)):
        raise RuntimeError("A banned shared cytotoxic/NK gene entered Stage 2")

    x = np.load(matrix_path, mmap_mode="r")
    params = dict(config["final_candidate"])
    profile_name = str(config["final_profile"])
    development_rows = np.flatnonzero(labels.allow_fit.to_numpy())
    audit = audit_rows(labels, development_rows, config, stable_seed("v4.3-final", "audit"))
    if not labels.iloc[audit].allow_fit.all():
        raise RuntimeError("Lockbox row entered the development audit")
    stage1_oof, stage1_threshold, stage1_details = fit_stage1_oof(
        labels, audit, x, stage1_columns, params, config
    )
    stage2_oof = fit_stage2_oof(
        labels,
        audit,
        x,
        stage2_columns,
        stage1_oof,
        include_stage1_probability,
        params,
        config,
    )
    audit_frame = labels.iloc[audit].reset_index(drop=True)
    excluded = exclusion_flags(np.asarray(x[audit]), names)[2]
    rescue = receptor_rescue_flags(
        np.asarray(x[audit]), names, config.get("receptor_rescue", {}).get(profile_name)
    )
    threshold_scores = stage2_oof.copy()
    threshold_scores[rescue] = np.finfo(np.float32).max
    stage2_threshold, threshold_metrics = choose_profile_threshold(
        audit_frame,
        stage1_oof,
        stage1_threshold,
        threshold_scores,
        excluded,
        config["profiles"][profile_name],
    )
    if threshold_metrics.get("eligible") is not True:
        raise RuntimeError("Final balanced profile failed the frozen grouped-development gate")
    calls = (
        (stage1_oof >= stage1_threshold)
        & ((stage2_oof >= stage2_threshold) | rescue)
        & ~excluded
    )
    score = stage1_oof * stage2_oof
    score[rescue & (stage1_oof >= stage1_threshold) & ~excluded] = 1.0
    development_metrics = profile_metrics(audit_frame, calls, score)

    oof_table = audit_frame[
        ["cell_id", "atlas_row", "truth_class", "source_gse_id", "donor_id", "sample_id", "group_id"]
    ].copy()
    oof_table["stage1_score"] = stage1_oof
    oof_table["stage2_score"] = stage2_oof
    oof_table["receptor_rescue"] = rescue
    oof_table["cd4_or_treg_excluded"] = excluded
    oof_table["predicted_gdT"] = calls
    oof_path = table_dir / "grouped_oof_predictions.parquet"
    oof_table.to_parquet(oof_path, index=False, compression="zstd")

    stage1_path, stage2_path, n_stage1_fit, n_stage2_fit = fit_final_models(
        labels,
        x,
        stage1_columns,
        stage2_columns,
        include_stage1_probability,
        params,
        config,
        model_dir,
    )
    contract = {
        "protocol_version": config["protocol_version"],
        "status": "DEVELOPMENT_FROZEN_NOT_PROMOTED",
        "development_frozen": True,
        "lockbox_scored": False,
        "model_promoted": False,
        "normalization": "log1p(raw_counts_per_10000)",
        "candidate_id": candidate_id(params),
        "candidate_parameters": params,
        "profile": profile_name,
        "stage1_threshold": stage1_threshold,
        "stage2_threshold": float(stage2_threshold),
        "receptor_rescue": config.get("receptor_rescue", {}).get(profile_name, {}),
        "stage1_feature_names": features.iloc[stage1_columns].gene.astype(str).tolist(),
        "stage2_feature_names": features.iloc[stage2_columns].gene.astype(str).tolist(),
        "stage2_includes_stage1_probability": include_stage1_probability,
        "n_stage1_fit_rows": n_stage1_fit,
        "n_stage2_fit_rows": n_stage2_fit,
        "n_grouped_oof_rows": len(audit),
        "stage1_oof_gate": stage1_details,
        "stage2_oof_gate": threshold_metrics,
        "development_oof_metrics": development_metrics,
        "input_hashes": {
            "config": sha256_file(config_path),
            "labels": sha256_file(labels_path),
            "features": sha256_file(features_path),
            "feature_cache": cache_manifest["matrix_sha256"],
            "split_manifest": sha256_file(split_path),
        },
        "artifact_hashes": {
            "stage1_model": sha256_file(stage1_path),
            "stage2_model": sha256_file(stage2_path),
            "grouped_oof_predictions": sha256_file(oof_path),
        },
        "common_lockbox_decision": config["common_lockbox_decision"],
    }
    contract_path = model_dir / "model_contract.json"
    contract_path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n")
    summary = {
        "status": "PASS_DEVELOPMENT_FROZEN",
        "runtime_seconds": time.monotonic() - started,
        "model_contract": str(contract_path),
        "model_contract_sha256": sha256_file(contract_path),
        "development_frozen": True,
        "lockbox_scored": False,
        "model_promoted": False,
        "stage1_threshold": stage1_threshold,
        "stage2_threshold": float(stage2_threshold),
        "development_oof_metrics": development_metrics,
    }
    (log_dir / "final_fit_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
