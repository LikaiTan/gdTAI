#!/usr/bin/env python3
"""Development-only gdTAI V4.2 recovery with a source-regularized GPU logistic Stage 2."""

from __future__ import annotations

import argparse
import json
import math
import pickle
import sys
import time
from pathlib import Path

import cupy as cp
import numpy as np
import pandas as pd
from cuml.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from workflows.gdtai.train_gdtai_v4_2_nested import (
    POSITIVE_SOURCES,
    audit_rows,
    choose_profile_threshold,
    choose_stage1_threshold,
    exclusion_flags,
    fit_xgb,
    predict,
    profile_metrics,
    resolve,
    sampled_rows,
    sha256_file,
    source_balanced_weights,
    stable_seed,
)


DEFAULT_CONFIG = ROOT / "configs/models/gdtai/v4_2_recovery.json"


def panel_columns(feature_names: list[str], panel: str) -> np.ndarray:
    if panel == "all_frozen":
        return np.arange(len(feature_names), dtype=np.int64)
    controls = {"CD3D", "CD3E", "CD3G", "CD4", "FOXP3"}
    selected = [i for i, gene in enumerate(feature_names) if gene.startswith(("TRA", "TRB", "TRG", "TRD")) or gene in controls]
    return np.asarray(selected, dtype=np.int64)


def panel_values(x: np.ndarray, rows: np.ndarray, columns: np.ndarray, panel: str) -> np.ndarray:
    values = np.asarray(x[rows][:, columns])
    if panel == "receptor_detection":
        values = (values > 0).astype(np.float32)
    return values


def standardize_fit(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean = values.mean(axis=0, dtype=np.float64).astype(np.float32)
    scale = values.std(axis=0, dtype=np.float64).astype(np.float32)
    scale[scale < 1e-6] = 1.0
    return ((values - mean) / scale).astype(np.float32), mean, scale


def standardize_apply(values: np.ndarray, mean: np.ndarray, scale: np.ndarray) -> np.ndarray:
    return ((values - mean) / scale).astype(np.float32)


def fit_logistic(values: np.ndarray, y: np.ndarray, weight: np.ndarray, c_value: float,
                 max_iter: int, tol: float) -> LogisticRegression:
    model = LogisticRegression(C=c_value, penalty="l2", max_iter=max_iter, tol=tol, output_type="cupy")
    model.fit(cp.asarray(values), cp.asarray(y), sample_weight=cp.asarray(weight))
    return model


def predict_logistic(model: LogisticRegression, values: np.ndarray, chunk: int = 250_000) -> np.ndarray:
    output = []
    for start in range(0, len(values), chunk):
        probability = model.predict_proba(cp.asarray(values[start:start + chunk]))[:, 1]
        output.append(cp.asnumpy(probability))
    return np.concatenate(output).astype(np.float32) if output else np.asarray([], dtype=np.float32)


def selected_candidate(table: pd.DataFrame, profile: str, beta: float) -> tuple[pd.Series, dict] | None:
    choices = []
    for _, row in table.iterrows():
        payload = json.loads(row.profiles_json)[profile]
        metrics = payload["metrics"]
        if metrics.get("eligible") is True:
            objective = metrics["f1"] if beta == 1 else metrics["f0_5"]
            choices.append((objective, metrics.get("pr_auc", -math.inf), row, payload))
    if not choices:
        return None
    _, _, row, payload = max(choices, key=lambda item: (item[0], item[1]))
    return row, payload


def save_logistic(path: Path, model: LogisticRegression, mean: np.ndarray, scale: np.ndarray,
                  feature_names: list[str], metadata: dict) -> None:
    payload = {
        "model": model,
        "mean": mean,
        "scale": scale,
        "feature_names": feature_names,
        "metadata": metadata,
        "diagnostic_only": True,
        "lockbox_scored": False,
        "model_promoted": False,
    }
    with path.open("wb") as handle:
        pickle.dump(payload, handle, protocol=pickle.HIGHEST_PROTOCOL)


def run_outer(labels: pd.DataFrame, x: np.ndarray, names: list[str], stage1_columns: np.ndarray,
              heldout: str, base: dict, recovery: dict, model_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    development = labels.allow_fit.to_numpy()
    outer_test_rows = np.flatnonzero(development & labels.source_gse_id.eq(heldout).to_numpy())
    outer_train = np.flatnonzero(development & ~labels.source_gse_id.eq(heldout).to_numpy())
    audit = audit_rows(labels, outer_train, base, stable_seed(heldout, "recovery_audit"))
    audit_frame = labels.iloc[audit].reset_index(drop=True)
    inner = audit_frame.inner_fold.to_numpy(np.int8)
    excluded_audit = exclusion_flags(np.asarray(x[audit]), names)[2]
    params = recovery["stage1_candidate"]

    stage1_oof = np.full(len(audit), np.nan, dtype=np.float32)
    for fold in (0, 1, 2):
        fit_local = np.flatnonzero((inner != fold) & audit_frame.truth_class.isin(
            ["gdT_gold", "abT_gold", "single_ab_support", "nk_tier1", "nk_tier2"]
        ).to_numpy())
        fit_frame = audit_frame.iloc[fit_local]
        t_local = np.flatnonzero(~fit_frame.truth_class.str.startswith("nk_").to_numpy())
        nk_local = np.flatnonzero(fit_frame.truth_class.str.startswith("nk_").to_numpy())
        t_selected = sampled_rows(fit_frame.reset_index(drop=True), t_local, base["maximum_stage1_t_cells"], stable_seed(heldout, fold, "recovery_s1"))
        train_local = fit_local[np.concatenate([t_selected, nk_local])]
        validation = np.flatnonzero(inner == fold)
        model = fit_xgb(
            np.asarray(x[audit[train_local]][:, stage1_columns]),
            (~audit_frame.iloc[train_local].truth_class.str.startswith("nk_")).to_numpy(np.int8),
            source_balanced_weights(audit_frame, train_local, 1), params, base["fixed"],
            stable_seed(heldout, fold, "recovery_s1_fit"),
            (np.asarray(x[audit[validation]][:, stage1_columns]),
             (~audit_frame.iloc[validation].truth_class.str.startswith("nk_")).to_numpy(np.int8)),
        )
        stage1_oof[validation] = predict(model, np.asarray(x[audit[validation]][:, stage1_columns]))
    stage1_threshold, stage1_details = choose_stage1_threshold(audit_frame, stage1_oof, base["stage1"])
    if stage1_details.get("eligible", True) is not True:
        raise RuntimeError(f"Recovery Stage 1 failed frozen gates for {heldout}")

    candidate_rows = []
    candidate_oof = {}
    for panel in recovery["stage2_panels"]:
        columns = panel_columns(names, panel)
        for c_value in recovery["stage2_C"]:
            candidate_id = f"{panel}__C-{c_value:g}"
            print(f"[{heldout}] Stage 2 {candidate_id}", flush=True)
            oof = np.full(len(audit), np.nan, dtype=np.float32)
            for fold in (0, 1, 2):
                train_local = np.flatnonzero((inner != fold) & audit_frame.truth_class.isin(
                    ["gdT_gold", "abT_gold", "single_ab_support"]
                ).to_numpy())
                train_frame = audit_frame.iloc[train_local]
                positive = np.flatnonzero(train_frame.truth_class.eq("gdT_gold").to_numpy())
                negative = np.flatnonzero(~train_frame.truth_class.eq("gdT_gold").to_numpy())
                negative = sampled_rows(train_frame.reset_index(drop=True), negative, base["maximum_stage2_negative_cells"], stable_seed(heldout, candidate_id, fold, "recovery_neg"))
                chosen_local = train_local[np.concatenate([positive, negative])]
                train_values, mean, scale = standardize_fit(panel_values(x, audit[chosen_local], columns, panel))
                y = audit_frame.iloc[chosen_local].truth_class.eq("gdT_gold").to_numpy(np.int32)
                weight = source_balanced_weights(audit_frame, chosen_local, 2)
                model = fit_logistic(train_values, y, weight, c_value, recovery["stage2_max_iter"], recovery["stage2_tol"])
                validation = np.flatnonzero(inner == fold)
                validation_values = standardize_apply(panel_values(x, audit[validation], columns, panel), mean, scale)
                oof[validation] = predict_logistic(model, validation_values)
            profiles = {}
            for profile, spec in base["profiles"].items():
                threshold, metrics = choose_profile_threshold(
                    audit_frame, stage1_oof, stage1_threshold, oof, excluded_audit, spec
                )
                profiles[profile] = {"threshold": threshold, "metrics": metrics}
            candidate_oof[candidate_id] = oof
            candidate_rows.append({
                "heldout_source": heldout,
                "candidate_id": candidate_id,
                "panel": panel,
                "C": c_value,
                "n_features": len(columns),
                "complete_oof": bool(np.isfinite(oof).all()),
                "profiles_json": json.dumps(profiles, sort_keys=True),
            })
    candidate_table = pd.DataFrame(candidate_rows)

    outer_frame = labels.iloc[outer_train]
    t_local = np.flatnonzero(~outer_frame.truth_class.str.startswith("nk_").to_numpy())
    nk_local = np.flatnonzero(outer_frame.truth_class.str.startswith("nk_").to_numpy())
    t_selected = sampled_rows(outer_frame.reset_index(drop=True), t_local, base["maximum_stage1_t_cells"], stable_seed(heldout, "recovery_final_s1_t"))
    final_stage1_rows = outer_train[np.concatenate([t_selected, nk_local])]
    final_stage1 = fit_xgb(
        np.asarray(x[final_stage1_rows][:, stage1_columns]),
        (~labels.iloc[final_stage1_rows].truth_class.str.startswith("nk_")).to_numpy(np.int8),
        source_balanced_weights(labels, final_stage1_rows, 1), params, base["fixed"],
        stable_seed(heldout, "recovery_final_s1"), None,
    )
    outer_stage1 = predict(final_stage1, np.asarray(x[outer_test_rows][:, stage1_columns]))
    outer_excluded = exclusion_flags(np.asarray(x[outer_test_rows]), names)[2]
    outer_test_frame = labels.iloc[outer_test_rows].reset_index(drop=True)
    outer_model_dir = model_dir / f"outer_{heldout}"
    outer_model_dir.mkdir(parents=True, exist_ok=True)
    final_stage1.save_model(outer_model_dir / "stage1.ubj")

    outer_results = []
    fitted = {}
    profile_contracts = {}
    for profile, profile_spec in base["profiles"].items():
        selection = selected_candidate(candidate_table, profile, profile_spec["objective_beta"])
        if selection is None:
            outer_results.append({"heldout_source": heldout, "profile": profile, "eligible_inner": False})
            continue
        selected_row, payload = selection
        candidate_id = selected_row.candidate_id
        panel, c_value = selected_row.panel, float(selected_row.C)
        columns = panel_columns(names, panel)
        if candidate_id not in fitted:
            train_local = np.flatnonzero(audit_frame.truth_class.isin(["gdT_gold", "abT_gold", "single_ab_support"]).to_numpy())
            train_frame = audit_frame.iloc[train_local]
            positive = np.flatnonzero(train_frame.truth_class.eq("gdT_gold").to_numpy())
            negative = np.flatnonzero(~train_frame.truth_class.eq("gdT_gold").to_numpy())
            negative = sampled_rows(train_frame.reset_index(drop=True), negative, base["maximum_stage2_negative_cells"], stable_seed(heldout, candidate_id, "recovery_final_neg"))
            chosen = train_local[np.concatenate([positive, negative])]
            train_values, mean, scale = standardize_fit(panel_values(x, audit[chosen], columns, panel))
            model = fit_logistic(
                train_values, audit_frame.iloc[chosen].truth_class.eq("gdT_gold").to_numpy(np.int32),
                source_balanced_weights(audit_frame, chosen, 2), c_value,
                recovery["stage2_max_iter"], recovery["stage2_tol"],
            )
            fitted[candidate_id] = (model, mean, scale, columns)
            save_logistic(
                outer_model_dir / f"stage2_{candidate_id}.pkl", model, mean, scale,
                [names[i] for i in columns], {"heldout_source": heldout, "candidate_id": candidate_id},
            )
        model, mean, scale, columns = fitted[candidate_id]
        outer_stage2 = predict_logistic(model, standardize_apply(panel_values(x, outer_test_rows, columns, panel), mean, scale))
        threshold = float(payload["threshold"])
        calls = (outer_stage1 >= stage1_threshold) & (outer_stage2 >= threshold) & ~outer_excluded
        metrics = profile_metrics(outer_test_frame, calls, outer_stage1 * outer_stage2)
        outer_results.append({
            "heldout_source": heldout, "profile": profile, "eligible_inner": True,
            "stage1_threshold": stage1_threshold, "stage2_candidate_id": candidate_id,
            "stage2_threshold": threshold,
            **{key: value for key, value in metrics.items() if not isinstance(value, dict)},
            "recall_by_source_json": json.dumps(metrics["recall_by_source"], sort_keys=True),
            "tier1_nk_fpr_by_source_json": json.dumps(metrics["tier1_nk_fpr_by_source"], sort_keys=True),
            "tier2_nk_fpr_by_source_json": json.dumps(metrics["tier2_nk_fpr_by_source"], sort_keys=True),
        })
        profile_contracts[profile] = {
            "stage2_candidate_id": candidate_id,
            "stage2_threshold": threshold,
            "eligible_inner": True,
        }
    contract = {
        "heldout_source": heldout,
        "stage1_threshold": stage1_threshold,
        "stage1_details": stage1_details,
        "stage1_feature_names": [names[i] for i in stage1_columns],
        "profiles": profile_contracts,
        "diagnostic_only": True,
        "lockbox_scored": False,
        "model_promoted": False,
    }
    (outer_model_dir / "contract.json").write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n")
    return candidate_table, pd.DataFrame(outer_results)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--heldout-source", action="append", default=[])
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--panel", choices=["receptor_detection", "receptor_controls", "all_frozen"])
    parser.add_argument("--C", type=float, dest="c_value")
    args = parser.parse_args()
    recovery = json.loads(resolve(args.config).read_text())
    if args.panel:
        recovery["stage2_panels"] = [args.panel]
    if args.c_value is not None:
        recovery["stage2_C"] = [args.c_value]
    if args.smoke:
        recovery["stage2_panels"] = recovery["stage2_panels"][:1]
        recovery["stage2_C"] = recovery["stage2_C"][:1]
    base = json.loads(resolve(recovery["base_training_config"]).read_text())
    table_dir, model_dir, log_dir = map(resolve, [recovery["output_table_dir"], recovery["output_model_dir"], recovery["output_log_dir"]])
    for path in (table_dir, model_dir, log_dir):
        path.mkdir(parents=True, exist_ok=True)
    cache_manifest = json.loads(resolve(base["cache_manifest"]).read_text())
    label_path, feature_path, matrix_path = map(resolve, [base["label_manifest"], base["feature_manifest"], base["cache_matrix"]])
    if cache_manifest["label_manifest_sha256"] != sha256_file(label_path) or cache_manifest["feature_manifest_sha256"] != sha256_file(feature_path):
        raise RuntimeError("Recovery inputs disagree with the frozen cache contract")
    labels = pd.read_parquet(label_path).reset_index(drop=True)
    split = pd.read_csv(resolve(base["split_manifest"]))
    labels["inner_fold"] = labels.group_id.map(split.drop_duplicates("group_id").set_index("group_id")["inner_fold"])
    if labels.loc[labels.allow_fit, "inner_fold"].isna().any():
        raise RuntimeError("Frozen inner-fold assignment is incomplete")
    features = pd.read_csv(feature_path).sort_values("feature_index")
    names = features.gene.astype(str).tolist()
    stage1_columns = np.flatnonzero(features.stage1.to_numpy(bool))
    x = np.load(matrix_path, mmap_mode="r")
    heldouts = args.heldout_source or POSITIVE_SOURCES
    started = time.monotonic()
    candidates, outer = [], []
    for heldout in heldouts:
        print(f"Starting recovery outer fold: {heldout}", flush=True)
        candidate_table, outer_table = run_outer(labels, x, names, stage1_columns, heldout, base, recovery, model_dir)
        candidates.append(candidate_table)
        outer.append(outer_table)
        pd.concat(candidates).to_csv(table_dir / "nested_candidates.csv", index=False)
        pd.concat(outer).to_csv(table_dir / "nested_outer_metrics.csv", index=False)
    result = {
        "status": "PASS_RECOVERY_NESTED_RUN",
        "heldout_sources": heldouts,
        "runtime_seconds": time.monotonic() - started,
        "lockbox_scored": False,
        "model_promoted": False,
    }
    (log_dir / "nested_summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
