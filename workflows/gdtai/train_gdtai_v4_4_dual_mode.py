#!/usr/bin/env python3
"""Train and freeze gdTAI V4.4 with validation-only dual operating modes."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import cupy as cp
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

from train_gdtai_v4_2_nested import (
    exclusion_flags,
    sampled_rows,
    sha256_file,
    source_balanced_weights,
    stable_seed,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "configs/models/gdtai/v4_4_dual_mode_training.json"
FIT_CLASSES = ["gdT_gold", "abT_gold", "single_ab_support", "nk_tier1", "nk_tier2"]
STAGE2_AB_CLASSES = ["abT_gold", "single_ab_support"]


def resolve(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def group_tiebreak(group_id: str) -> int:
    digest = hashlib.sha256(f"gdtai-v4.4::{group_id}".encode()).digest()
    return int.from_bytes(digest[:8], "little")


def balanced_holdout_partitions(labels: pd.DataFrame, heldout: np.ndarray) -> dict[str, int]:
    """Balance whole groups across calibration/threshold sets within each source."""
    frame = labels.loc[heldout, ["source_gse_id", "truth_class", "group_id"]].copy()
    counts = frame.groupby(
        ["source_gse_id", "group_id", "truth_class"], observed=True
    ).size().unstack(fill_value=0)
    assignments: dict[str, int] = {}
    for source, local in counts.groupby(level="source_gse_id", observed=True):
        local = local.droplevel("source_gse_id")
        totals = local.sum(axis=0).to_numpy(float)
        scale = np.where(totals > 0, totals, 1.0)
        order = sorted(
            local.index.astype(str),
            key=lambda group: (
                -float(np.linalg.norm(local.loc[group].to_numpy(float) / scale)),
                group_tiebreak(group),
            ),
        )
        loads = np.zeros((2, local.shape[1]), dtype=float)
        for group in order:
            vector = local.loc[group].to_numpy(float)
            costs = []
            for target in (0, 1):
                candidate = loads.copy()
                candidate[target] += vector
                costs.append(float(np.square((candidate[0] - candidate[1]) / scale).sum()))
            if math.isclose(costs[0], costs[1], rel_tol=0.0, abs_tol=1e-14):
                target = group_tiebreak(f"{source}::{group}") & 1
            else:
                target = int(np.argmin(costs))
            assignments[str(group)] = target
            loads[target] += vector
    return assignments


def load_inputs(config: dict) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, dict]:
    labels_path = resolve(config["label_manifest"])
    features_path = resolve(config["feature_manifest"])
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
    matrix = np.load(resolve(config["cache_matrix"]), mmap_mode="r")
    if tuple(matrix.shape) != (len(labels), len(features)):
        raise RuntimeError(f"Feature-cache shape mismatch: {matrix.shape}")
    return labels, features, matrix, cache_manifest


def assign_development_partitions(labels: pd.DataFrame, config: dict) -> pd.Series:
    spec = config["development_partition"]
    fit_folds = set(map(int, spec["fit_inner_folds"]))
    holdout_fold = int(spec["heldout_inner_fold"])
    partition = pd.Series("locked_test", index=labels.index, dtype="object")
    development = labels.allow_fit.to_numpy(bool)
    fit = development & labels.inner_fold.isin(fit_folds).to_numpy()
    heldout = development & labels.inner_fold.eq(holdout_fold).to_numpy()
    partition.loc[fit] = "fit"
    assignments = balanced_holdout_partitions(labels, heldout)
    parity = labels.loc[heldout, "group_id"].astype(str).map(assignments)
    if parity.isna().any():
        raise RuntimeError("A held-out development group lacks a calibration/threshold assignment")
    partition.loc[parity.index[parity.eq(int(spec["calibration_partition"]))]] = "calibration"
    partition.loc[parity.index[parity.eq(int(spec["threshold_validation_partition"]))]] = "threshold_validation"
    if (partition.loc[development] == "locked_test").any():
        raise RuntimeError("A development row lacks a fit/calibration/threshold partition")
    for name in ("fit", "calibration", "threshold_validation"):
        if not (partition == name).any():
            raise RuntimeError(f"Empty development partition: {name}")
    grouped = pd.DataFrame({"group_id": labels.group_id, "partition": partition})
    if grouped.groupby("group_id", observed=True).partition.nunique().max() != 1:
        raise RuntimeError("A biological group crosses development partitions")
    return partition


def fit_xgb_model(x_values: np.ndarray, y: np.ndarray, weight: np.ndarray, params: dict,
                  fixed: dict, seed: int) -> xgb.XGBClassifier:
    model = xgb.XGBClassifier(**(dict(fixed) | dict(params) | {"random_state": seed}))
    model.fit(cp.asarray(np.asarray(x_values, dtype=np.float32)), cp.asarray(y),
              sample_weight=cp.asarray(weight), verbose=False)
    return model


def predict_xgb(model: xgb.XGBClassifier, values: np.ndarray, chunk: int = 250_000) -> np.ndarray:
    output: list[np.ndarray] = []
    for start in range(0, len(values), chunk):
        block = cp.asarray(np.asarray(values[start:start + chunk], dtype=np.float32))
        score = model.predict_proba(block)[:, 1]
        output.append(cp.asnumpy(score) if isinstance(score, cp.ndarray) else np.asarray(score))
    return np.concatenate(output).astype(np.float32) if output else np.empty(0, dtype=np.float32)


@dataclass(frozen=True)
class PlattCalibration:
    coefficient: float
    intercept: float

    def apply(self, score: np.ndarray) -> np.ndarray:
        clipped = np.clip(np.asarray(score, dtype=np.float64), 1e-7, 1 - 1e-7)
        logit = np.log(clipped / (1 - clipped))
        linear = np.clip(self.coefficient * logit + self.intercept, -40, 40)
        return (1.0 / (1.0 + np.exp(-linear))).astype(np.float32)

    def as_dict(self) -> dict[str, float]:
        return {"coefficient": self.coefficient, "intercept": self.intercept}


def fit_platt(score: np.ndarray, y: np.ndarray, weight: np.ndarray, seed: int) -> PlattCalibration:
    clipped = np.clip(np.asarray(score, dtype=np.float64), 1e-7, 1 - 1e-7)
    logit = np.log(clipped / (1 - clipped)).reshape(-1, 1)
    model = LogisticRegression(solver="lbfgs", max_iter=500, random_state=seed)
    model.fit(logit, y, sample_weight=weight)
    return PlattCalibration(float(model.coef_[0, 0]), float(model.intercept_[0]))


def feature_columns(features: pd.DataFrame, config: dict) -> tuple[np.ndarray, np.ndarray]:
    names = features.gene.astype(str)
    stage1_names = set(config["stage1"]["feature_genes"])
    stage2_controls = set(config["stage2"]["control_genes"])
    stage1 = np.flatnonzero(names.isin(stage1_names).to_numpy())
    stage2 = np.flatnonzero(features.feature_class.eq("TCR").to_numpy() | names.isin(stage2_controls).to_numpy())
    missing_stage1 = sorted(stage1_names - set(names.iloc[stage1]))
    missing_controls = sorted(stage2_controls - set(names.iloc[stage2]))
    if missing_stage1 or missing_controls:
        raise RuntimeError(f"Missing frozen features: stage1={missing_stage1}, stage2={missing_controls}")
    banned = {"KLRD1", "NKG7", "GNLY", "FCER1G", "TYROBP", "FCGR3A"}
    effective = set(names.iloc[stage1]) | set(names.iloc[stage2])
    if banned & effective:
        raise RuntimeError(f"Banned shared NK/cytotoxic genes entered V4.4: {sorted(banned & effective)}")
    return stage1, stage2


def aggregate(values: np.ndarray, indexes: list[int], operation: str) -> np.ndarray:
    if not indexes:
        return np.zeros(len(values), dtype=np.float32)
    block = values[:, indexes]
    if operation == "max":
        return np.max(block, axis=1).astype(np.float32)
    if operation == "sum":
        return np.sum(block, axis=1, dtype=np.float32)
    if operation == "count":
        return np.sum(block > 0, axis=1, dtype=np.float32)
    raise ValueError(operation)


def compose_stage2(values: np.ndarray, names: list[str]) -> tuple[np.ndarray, list[str]]:
    base = np.asarray(values, dtype=np.float32)
    lookup = {gene: index for index, gene in enumerate(names)}
    groups = {
        "TRDV_dedicated": [i for i, g in enumerate(names) if g.startswith("TRDV")],
        "TRDV_all": [i for i, g in enumerate(names) if g.startswith("TRDV") or (g.startswith("TRAV") and "DV" in g)],
        "TRDJ": [i for i, g in enumerate(names) if g.startswith("TRDJ")],
        "TRGV": [i for i, g in enumerate(names) if g.startswith("TRGV")],
        "TRGC": [i for i, g in enumerate(names) if g.startswith("TRGC")],
        "TRAV_alpha": [i for i, g in enumerate(names) if g.startswith("TRAV") and "DV" not in g],
        "TRAJ": [i for i, g in enumerate(names) if g.startswith("TRAJ")],
        "TRBV": [i for i, g in enumerate(names) if g.startswith("TRBV")],
        "TRBJ": [i for i, g in enumerate(names) if g.startswith("TRBJ")],
        "TRBC": [i for i, g in enumerate(names) if g.startswith("TRBC")],
    }
    engineered: list[np.ndarray] = []
    engineered_names: list[str] = []

    def add(name: str, vector: np.ndarray) -> None:
        engineered_names.append(name)
        engineered.append(np.asarray(vector, dtype=np.float32))

    trdc = base[:, lookup["TRDC"]]
    trac = base[:, lookup["TRAC"]]
    add("eng_TRDC_detected", trdc > 0)
    add("eng_TRAC_detected", trac > 0)
    cached: dict[tuple[str, str], np.ndarray] = {}
    for group, indexes in groups.items():
        for operation in ("max", "sum", "count"):
            value = aggregate(base, indexes, operation)
            cached[(group, operation)] = value
            add(f"eng_{group}_{operation}", value)
    add("eng_TRDV_dedicated_detected", cached[("TRDV_dedicated", "max")] > 0)
    add("eng_TRDV_all_detected", cached[("TRDV_all", "max")] > 0)
    add("eng_TRGV_detected", cached[("TRGV", "max")] > 0)
    add("eng_TRDC_TRDV_joint", np.minimum(trdc, cached[("TRDV_all", "max")]))
    add("eng_TRGC_TRGV_joint", np.minimum(cached[("TRGC", "max")], cached[("TRGV", "max")]))
    gd_sum = trdc + cached[("TRDV_all", "sum")] + cached[("TRDJ", "sum")] + cached[("TRGC", "sum")] + cached[("TRGV", "sum")]
    ab_v_sum = cached[("TRAV_alpha", "sum")] + cached[("TRBV", "sum")]
    ab_j_sum = cached[("TRAJ", "sum")] + cached[("TRBJ", "sum")]
    ab_constant_sum = trac + cached[("TRBC", "sum")]
    add("eng_gd_receptor_sum", gd_sum)
    add("eng_ab_v_sum", ab_v_sum)
    add("eng_ab_j_sum", ab_j_sum)
    add("eng_ab_constant_sum", ab_constant_sum)
    add("eng_gd_minus_ab_receptor_sum", gd_sum - ab_v_sum - ab_j_sum - ab_constant_sum)
    output = np.column_stack([base, *engineered]).astype(np.float32)
    return output, names + engineered_names


def make_stage1_rows(labels: pd.DataFrame, candidate_rows: np.ndarray, maximum_t: int, seed: int) -> np.ndarray:
    frame = labels.iloc[candidate_rows]
    t_local = np.flatnonzero(~frame.truth_class.str.startswith("nk_").to_numpy())
    nk_local = np.flatnonzero(frame.truth_class.str.startswith("nk_").to_numpy())
    t_selected = sampled_rows(frame.reset_index(drop=True), t_local, maximum_t, seed)
    return candidate_rows[np.unique(np.concatenate([t_selected, nk_local]))]


def make_stage2_rows(labels: pd.DataFrame, candidate_rows: np.ndarray, maximum_ab: int,
                     seed: int) -> np.ndarray:
    frame = labels.iloc[candidate_rows]
    positive = np.flatnonzero(frame.truth_class.eq("gdT_gold").to_numpy())
    ab = np.flatnonzero(frame.truth_class.isin(STAGE2_AB_CLASSES).to_numpy())
    nk = np.flatnonzero(frame.truth_class.isin(["nk_tier1", "nk_tier2"]).to_numpy())
    ab_selected = sampled_rows(frame.reset_index(drop=True), ab, maximum_ab, seed)
    return candidate_rows[np.unique(np.concatenate([positive, ab_selected, nk]))]


def source_rates(frame: pd.DataFrame, calls: np.ndarray, truth: str) -> dict[str, float]:
    truth_mask = frame.truth_class.eq(truth).to_numpy()
    return {
        str(source): float(calls[truth_mask & frame.source_gse_id.eq(source).to_numpy()].mean())
        for source in sorted(frame.loc[truth_mask, "source_gse_id"].unique())
    }


def choose_stage1_threshold(frame: pd.DataFrame, score: np.ndarray, spec: dict) -> tuple[float, dict]:
    thresholds = np.unique(np.quantile(score, np.linspace(0, 1, 4001)))
    best: tuple[float, dict] | None = None
    for threshold in thresholds[::-1]:
        calls = score >= threshold
        gd = source_rates(frame, calls, "gdT_gold")
        ab = source_rates(frame, calls, "abT_gold")
        eligible = (
            (not gd or min(gd.values()) >= float(spec["minimum_gdt_recall_per_source"]))
            and (not ab or min(ab.values()) >= float(spec["minimum_abt_recall_per_source"]))
        )
        if eligible:
            best = (float(threshold), {
                "eligible": True,
                "gdt_recall_by_source": gd,
                "abt_recall_by_source": ab,
                "tier1_nk_passage_by_source": source_rates(frame, calls, "nk_tier1"),
            })
            break
    if best is None:
        raise RuntimeError("No Stage-1 validation threshold satisfies the frozen T-cell recall rules")
    return best


def fbeta(tp: int, fp: int, fn: int, beta: float) -> float:
    b2 = beta * beta
    denominator = (1 + b2) * tp + b2 * fn + fp
    return (1 + b2) * tp / denominator if denominator else 0.0


def threshold_frontier(y: np.ndarray, score: np.ndarray) -> pd.DataFrame:
    order = np.argsort(-score, kind="mergesort")
    sorted_score = score[order]
    sorted_y = y[order].astype(np.int8)
    cumulative_tp = np.cumsum(sorted_y, dtype=np.int64)
    cumulative_fp = np.cumsum(1 - sorted_y, dtype=np.int64)
    endpoints = np.r_[np.flatnonzero(sorted_score[:-1] != sorted_score[1:]), len(sorted_score) - 1]
    tp = cumulative_tp[endpoints]
    fp = cumulative_fp[endpoints]
    positive = int(sorted_y.sum())
    negative = int(len(sorted_y) - positive)
    fn = positive - tp
    tn = negative - fp
    precision = np.divide(tp, tp + fp, out=np.zeros(len(tp), dtype=float), where=(tp + fp) > 0)
    recall = np.divide(tp, positive, out=np.zeros(len(tp), dtype=float), where=positive > 0)
    specificity = np.divide(tn, negative, out=np.zeros(len(tp), dtype=float), where=negative > 0)
    f1 = np.asarray([fbeta(int(a), int(b), int(c), 1.0) for a, b, c in zip(tp, fp, fn)])
    f0_5 = np.asarray([fbeta(int(a), int(b), int(c), 0.5) for a, b, c in zip(tp, fp, fn)])
    return pd.DataFrame({
        "threshold": sorted_score[endpoints], "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "precision": precision, "recall": recall, "specificity": specificity,
        "f1": f1, "f0_5": f0_5,
    })


def select_modes(frontier: pd.DataFrame, config: dict) -> dict[str, dict]:
    selected: dict[str, dict] = {}
    for mode, spec in config["threshold_selection"]["operating_modes"].items():
        if spec.get("hard_fpr_constraint") is not False:
            raise RuntimeError("V4.4 operating modes must not impose an FPR threshold")
        objective = str(spec["objective"])
        ranked = frontier.sort_values(
            [objective, "precision", "recall", "threshold"],
            ascending=[False, False, False, False],
        )
        selected[mode] = json_safe(ranked.iloc[0].to_dict())
        selected[mode]["objective"] = objective
        selected[mode]["threshold_source"] = "grouped_development_threshold_validation"
    return selected


def binary_metrics(y: np.ndarray, score: np.ndarray, calls: np.ndarray) -> dict:
    tp = int((calls & (y == 1)).sum())
    fp = int((calls & (y == 0)).sum())
    fn = int((~calls & (y == 1)).sum())
    tn = int((~calls & (y == 0)).sum())
    valid_auc = np.unique(y).size == 2
    return {
        "n": int(len(y)), "positive": int((y == 1).sum()), "negative": int((y == 0).sum()),
        "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "precision": tp / (tp + fp) if tp + fp else 0.0,
        "recall": tp / (tp + fn) if tp + fn else 0.0,
        "specificity": tn / (tn + fp) if tn + fp else 0.0,
        "f1": fbeta(tp, fp, fn, 1.0), "f0_5": fbeta(tp, fp, fn, 0.5),
        "mcc": float(matthews_corrcoef(y, calls)) if valid_auc else math.nan,
        "roc_auc": float(roc_auc_score(y, score)) if valid_auc else math.nan,
        "pr_auc": float(average_precision_score(y, score)) if valid_auc else math.nan,
    }


def exclusion_for_rows(matrix: np.ndarray, rows: np.ndarray, all_names: list[str]) -> np.ndarray:
    return exclusion_flags(np.asarray(matrix[rows]), all_names)[2]


def validation_metrics(labels: pd.DataFrame, rows: np.ndarray, score: np.ndarray,
                       selected: dict[str, dict]) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = labels.iloc[rows].reset_index(drop=True)
    primary = frame.truth_class.isin(["gdT_gold", "abT_gold", "nk_tier1"]).to_numpy()
    y = frame.truth_class.eq("gdT_gold").to_numpy(np.int8)
    overall, source_rows = [], []
    for mode, payload in selected.items():
        calls = score >= float(payload["threshold"])
        metrics = binary_metrics(y[primary], score[primary], calls[primary])
        metrics.update({"mode": mode, "threshold": float(payload["threshold"]), "evaluation_set": "threshold_validation_primary"})
        overall.append(metrics)
        for source, local in frame.groupby("source_gse_id", observed=True).groups.items():
            indexes = np.asarray(local, dtype=np.int64)
            source_primary = frame.iloc[indexes].truth_class.isin(["gdT_gold", "abT_gold", "nk_tier1"]).to_numpy()
            if not source_primary.any():
                continue
            local_y = y[indexes][source_primary]
            local_score = score[indexes][source_primary]
            local_calls = calls[indexes][source_primary]
            item = binary_metrics(local_y, local_score, local_calls)
            item.update({"mode": mode, "source_gse_id": source, "threshold": float(payload["threshold"])})
            source_rows.append(item)
    return pd.DataFrame(overall), pd.DataFrame(source_rows)


def train_bundle(labels: pd.DataFrame, matrix: np.ndarray, all_names: list[str],
                 stage1_columns: np.ndarray, stage2_columns: np.ndarray,
                 partition: pd.Series, config: dict, excluded_source: str | None = None) -> dict:
    allowed_source = np.ones(len(labels), dtype=bool)
    if excluded_source is not None:
        allowed_source &= labels.source_gse_id.ne(excluded_source).to_numpy()
    fit_pool = np.flatnonzero(
        partition.eq("fit").to_numpy() & allowed_source & labels.truth_class.isin(FIT_CLASSES).to_numpy()
    )
    calibration_rows = np.flatnonzero(
        partition.eq("calibration").to_numpy() & allowed_source & labels.truth_class.isin(FIT_CLASSES).to_numpy()
    )
    threshold_rows = np.flatnonzero(
        partition.eq("threshold_validation").to_numpy() & allowed_source & labels.truth_class.isin(FIT_CLASSES).to_numpy()
    )
    if not len(fit_pool) or not len(calibration_rows) or not len(threshold_rows):
        raise RuntimeError("An empty development partition entered model training")
    calibration_classes = labels.iloc[calibration_rows].truth_class.eq("gdT_gold")
    threshold_classes = labels.iloc[threshold_rows].truth_class.eq("gdT_gold")
    if calibration_classes.nunique() != 2:
        raise RuntimeError("Stage-2 calibration requires both gdT and non-gdT development cells")
    if threshold_classes.nunique() != 2:
        raise RuntimeError("Operating-threshold validation requires both gdT and non-gdT development cells")
    seed_tag = excluded_source or "final"

    stage1_rows = make_stage1_rows(
        labels, fit_pool, int(config["stage1"]["maximum_t_cells"]),
        stable_seed("v4.4", seed_tag, "stage1-rows"),
    )
    stage1_model = fit_xgb_model(
        np.asarray(matrix[stage1_rows][:, stage1_columns]),
        (~labels.iloc[stage1_rows].truth_class.str.startswith("nk_")).to_numpy(np.int8),
        source_balanced_weights(labels, stage1_rows, 1),
        config["stage1"]["xgboost"], config["xgboost_fixed"],
        stable_seed("v4.4", seed_tag, "stage1-fit"),
    )
    calibration_stage1_raw = predict_xgb(stage1_model, np.asarray(matrix[calibration_rows][:, stage1_columns]))
    stage1_calibrator = fit_platt(
        calibration_stage1_raw,
        (~labels.iloc[calibration_rows].truth_class.str.startswith("nk_")).to_numpy(np.int8),
        source_balanced_weights(labels, calibration_rows, 1),
        stable_seed("v4.4", seed_tag, "stage1-calibration"),
    )
    threshold_stage1 = stage1_calibrator.apply(
        predict_xgb(stage1_model, np.asarray(matrix[threshold_rows][:, stage1_columns]))
    )
    stage1_threshold, stage1_details = choose_stage1_threshold(
        labels.iloc[threshold_rows].reset_index(drop=True), threshold_stage1, config["stage1"]
    )

    stage2_rows = make_stage2_rows(
        labels, fit_pool, int(config["stage2"]["maximum_ab_negative_cells"]),
        stable_seed("v4.4", seed_tag, "stage2-rows"),
    )
    stage2_names = [all_names[index] for index in stage2_columns]
    stage2_train, effective_names = compose_stage2(np.asarray(matrix[stage2_rows][:, stage2_columns]), stage2_names)
    if len(effective_names) > 300:
        raise RuntimeError(f"V4.4 exceeds the 300-feature limit: {len(effective_names)}")
    stage2_model = fit_xgb_model(
        stage2_train,
        labels.iloc[stage2_rows].truth_class.eq("gdT_gold").to_numpy(np.int8),
        source_balanced_weights(labels, stage2_rows, 2),
        config["stage2"]["xgboost"], config["xgboost_fixed"],
        stable_seed("v4.4", seed_tag, "stage2-fit"),
    )
    calibration_stage2, _ = compose_stage2(np.asarray(matrix[calibration_rows][:, stage2_columns]), stage2_names)
    calibration_stage2_raw = predict_xgb(stage2_model, calibration_stage2)
    stage2_calibrator = fit_platt(
        calibration_stage2_raw,
        labels.iloc[calibration_rows].truth_class.eq("gdT_gold").to_numpy(np.int8),
        source_balanced_weights(labels, calibration_rows, 2),
        stable_seed("v4.4", seed_tag, "stage2-calibration"),
    )
    threshold_stage2, _ = compose_stage2(np.asarray(matrix[threshold_rows][:, stage2_columns]), stage2_names)
    threshold_stage2_score = stage2_calibrator.apply(predict_xgb(stage2_model, threshold_stage2))
    excluded = exclusion_for_rows(matrix, threshold_rows, all_names)
    eligible = (threshold_stage1 >= stage1_threshold) & ~excluded
    effective_score = threshold_stage2_score.copy()
    effective_score[~eligible] = -1.0
    threshold_spec = config["threshold_selection"]
    primary = labels.iloc[threshold_rows].truth_class.isin(
        threshold_spec["positive_truth_classes"] + threshold_spec["negative_truth_classes"]
    ).to_numpy()
    y = labels.iloc[threshold_rows].truth_class.eq("gdT_gold").to_numpy(np.int8)
    frontier = threshold_frontier(y[primary], effective_score[primary])
    selected = select_modes(frontier, config)
    overall, per_source = validation_metrics(labels, threshold_rows, effective_score, selected)
    return {
        "stage1_model": stage1_model,
        "stage2_model": stage2_model,
        "stage1_calibrator": stage1_calibrator,
        "stage2_calibrator": stage2_calibrator,
        "stage1_threshold": stage1_threshold,
        "stage1_details": stage1_details,
        "operating_modes": selected,
        "frontier": frontier,
        "validation_overall": overall,
        "validation_per_source": per_source,
        "threshold_rows": threshold_rows,
        "threshold_stage1_score": threshold_stage1,
        "threshold_stage2_score": threshold_stage2_score,
        "threshold_effective_score": effective_score,
        "threshold_excluded": excluded,
        "stage1_fit_rows": stage1_rows,
        "stage2_fit_rows": stage2_rows,
        "stage1_feature_names": [all_names[index] for index in stage1_columns],
        "stage2_base_feature_names": stage2_names,
        "stage2_effective_feature_names": effective_names,
        "calibration_rows": calibration_rows,
        "fit_pool": fit_pool,
    }


def evaluate_source_holdout(bundle: dict, labels: pd.DataFrame, matrix: np.ndarray,
                            all_names: list[str], stage1_columns: np.ndarray,
                            stage2_columns: np.ndarray, source: str) -> pd.DataFrame:
    rows = np.flatnonzero(labels.allow_fit.to_numpy() & labels.source_gse_id.eq(source).to_numpy()
                         & labels.truth_class.isin(FIT_CLASSES).to_numpy())
    stage1 = bundle["stage1_calibrator"].apply(
        predict_xgb(bundle["stage1_model"], np.asarray(matrix[rows][:, stage1_columns]))
    )
    stage2_values, _ = compose_stage2(
        np.asarray(matrix[rows][:, stage2_columns]), bundle["stage2_base_feature_names"]
    )
    stage2 = bundle["stage2_calibrator"].apply(predict_xgb(bundle["stage2_model"], stage2_values))
    excluded = exclusion_for_rows(matrix, rows, all_names)
    score = stage2.copy()
    score[(stage1 < bundle["stage1_threshold"]) | excluded] = -1.0
    frame = labels.iloc[rows].reset_index(drop=True)
    primary = frame.truth_class.isin(["gdT_gold", "abT_gold", "nk_tier1"]).to_numpy()
    y = frame.truth_class.eq("gdT_gold").to_numpy(np.int8)
    output = []
    for mode, payload in bundle["operating_modes"].items():
        calls = score >= float(payload["threshold"])
        metrics = binary_metrics(y[primary], score[primary], calls[primary]) if primary.any() else {}
        metrics.update({
            "heldout_source": source, "mode": mode, "threshold": float(payload["threshold"]),
            "gdt_recall": float(calls[frame.truth_class.eq("gdT_gold").to_numpy()].mean()) if frame.truth_class.eq("gdT_gold").any() else math.nan,
            "abt_fpr": float(calls[frame.truth_class.eq("abT_gold").to_numpy()].mean()) if frame.truth_class.eq("abT_gold").any() else math.nan,
            "tier1_nk_fpr": float(calls[frame.truth_class.eq("nk_tier1").to_numpy()].mean()) if frame.truth_class.eq("nk_tier1").any() else math.nan,
            "threshold_selected_without_heldout_source": True,
        })
        output.append(metrics)
    return pd.DataFrame(output)


def save_bundle(bundle: dict, labels: pd.DataFrame, features: pd.DataFrame, partition: pd.Series,
                config: dict, config_path: Path, cache_manifest: dict,
                source_holdout: pd.DataFrame, started: float) -> dict:
    table_dir = resolve(config["output_table_dir"])
    model_dir = resolve(config["output_model_dir"])
    log_dir = resolve(config["output_log_dir"])
    for path in (table_dir, model_dir, log_dir):
        path.mkdir(parents=True, exist_ok=True)

    partition_table = labels[["cell_id", "atlas_row", "source_gse_id", "truth_class", "group_id", "allow_fit", "cohort_role"]].copy()
    partition_table["inner_fold"] = labels.inner_fold
    partition_table["v4_4_partition"] = partition
    partition_path = table_dir / "development_partition_manifest.parquet"
    partition_table.to_parquet(partition_path, index=False, compression="zstd")

    feature_contract = pd.DataFrame({"feature_name": bundle["stage2_effective_feature_names"]})
    feature_contract["feature_type"] = np.where(feature_contract.feature_name.str.startswith("eng_"), "engineered_receptor_aggregate", "individual_gene")
    feature_path = table_dir / "stage2_feature_contract.csv"
    feature_contract.to_csv(feature_path, index=False)
    bundle["frontier"].to_csv(table_dir / "threshold_validation_frontier.csv", index=False)
    bundle["validation_overall"].to_csv(table_dir / "threshold_validation_metrics.csv", index=False)
    bundle["validation_per_source"].to_csv(table_dir / "threshold_validation_per_source.csv", index=False)
    source_holdout.to_csv(table_dir / "source_holdout_metrics.csv", index=False)

    threshold_rows = bundle["threshold_rows"]
    predictions = labels.iloc[threshold_rows][
        ["cell_id", "atlas_row", "source_gse_id", "truth_class", "group_id"]
    ].copy()
    predictions["stage1_probability"] = bundle["threshold_stage1_score"]
    predictions["stage2_probability"] = bundle["threshold_stage2_score"]
    predictions["effective_score"] = bundle["threshold_effective_score"]
    predictions["cd4_or_treg_excluded"] = bundle["threshold_excluded"]
    for mode, payload in bundle["operating_modes"].items():
        predictions[f"predicted_{mode}"] = predictions.effective_score.to_numpy() >= float(payload["threshold"])
    predictions_path = table_dir / "threshold_validation_predictions.parquet"
    predictions.to_parquet(predictions_path, index=False, compression="zstd")

    stage1_path = model_dir / "stage1_t_lineage_gate.ubj"
    stage2_path = model_dir / "stage2_receptor_classifier.ubj"
    bundle["stage1_model"].save_model(stage1_path)
    bundle["stage2_model"].save_model(stage2_path)
    calibration_path = model_dir / "platt_calibration.json"
    calibration_path.write_text(json.dumps({
        "stage1": bundle["stage1_calibrator"].as_dict(),
        "stage2": bundle["stage2_calibrator"].as_dict(),
    }, indent=2, sort_keys=True) + "\n")

    partition_counts = partition_table.groupby(["v4_4_partition", "source_gse_id", "truth_class"], observed=True).size().reset_index(name="n_cells")
    partition_counts.to_csv(table_dir / "development_partition_counts.csv", index=False)
    threshold_provenance = partition_counts[partition_counts.v4_4_partition.eq("threshold_validation")].copy()
    threshold_provenance["used_for_threshold"] = threshold_provenance.truth_class.isin(
        config["threshold_selection"]["positive_truth_classes"] + config["threshold_selection"]["negative_truth_classes"]
    )
    threshold_provenance.to_csv(table_dir / "threshold_provenance.csv", index=False)

    contract = {
        "protocol_version": config["protocol_version"],
        "status": "DEVELOPMENT_FROZEN_NOT_PROMOTED",
        "normalization": "log1p(raw_counts_per_10000)",
        "development_frozen": True,
        "test_scored": False,
        "model_promoted": False,
        "thresholds_selected_from_test": False,
        "threshold_selection_partition": "threshold_validation",
        "threshold_selection_truth": config["threshold_selection"],
        "stage1_threshold": bundle["stage1_threshold"],
        "stage1_threshold_details": bundle["stage1_details"],
        "operating_modes": bundle["operating_modes"],
        "stage1_feature_names": bundle["stage1_feature_names"],
        "stage2_base_feature_names": bundle["stage2_base_feature_names"],
        "stage2_effective_feature_names": bundle["stage2_effective_feature_names"],
        "stage2_includes_stage1_probability": False,
        "shared_nk_cytotoxic_features_used": False,
        "nk_negatives_in_stage2": ["nk_tier1", "nk_tier2"],
        "hard_trdv_expression_cutoff": False,
        "frozen_exclusions": config["frozen_exclusions"],
        "test_policy": config["test_policy"],
        "n_stage1_fit_rows": len(bundle["stage1_fit_rows"]),
        "n_stage2_fit_rows": len(bundle["stage2_fit_rows"]),
        "n_calibration_rows": len(bundle["calibration_rows"]),
        "n_threshold_validation_rows": len(bundle["threshold_rows"]),
        "threshold_validation_metrics": bundle["validation_overall"].to_dict("records"),
        "source_holdout_metrics": source_holdout.to_dict("records"),
        "input_hashes": {
            "config": sha256_file(config_path),
            "labels": sha256_file(resolve(config["label_manifest"])),
            "features": sha256_file(resolve(config["feature_manifest"])),
            "split_manifest": sha256_file(resolve(config["split_manifest"])),
            "feature_cache": cache_manifest["matrix_sha256"],
        },
        "artifact_hashes": {
            "stage1_model": sha256_file(stage1_path),
            "stage2_model": sha256_file(stage2_path),
            "platt_calibration": sha256_file(calibration_path),
            "development_partition_manifest": sha256_file(partition_path),
            "stage2_feature_contract": sha256_file(feature_path),
            "threshold_validation_predictions": sha256_file(predictions_path),
        },
    }
    contract_path = model_dir / "model_contract.json"
    contract_path.write_text(json.dumps(json_safe(contract), indent=2, sort_keys=True) + "\n")
    summary = {
        "status": "PASS_DEVELOPMENT_FROZEN",
        "runtime_seconds": time.monotonic() - started,
        "model_contract": str(contract_path),
        "model_contract_sha256": sha256_file(contract_path),
        "thresholds_selected_from_test": False,
        "test_scored": False,
        "model_promoted": False,
        "operating_modes": bundle["operating_modes"],
    }
    (log_dir / "training_summary.json").write_text(json.dumps(json_safe(summary), indent=2, sort_keys=True) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--skip-source-holdout", action="store_true")
    parser.add_argument("--partition-only", action="store_true")
    args = parser.parse_args()
    started = time.monotonic()
    config_path = resolve(args.config)
    config = json.loads(config_path.read_text())
    if config["test_policy"]["thresholds_may_use_test_cells"] is not False:
        raise RuntimeError("The V4.4 contract must forbid test-set threshold selection")
    if config["xgboost_fixed"].get("device") != "cuda":
        raise RuntimeError("The canonical V4.4 run requires GPU XGBoost")
    if cp.cuda.runtime.getDeviceCount() < 1:
        raise RuntimeError("No CUDA device is available")
    labels, features, matrix, cache_manifest = load_inputs(config)
    partition = assign_development_partitions(labels, config)
    table_dir = resolve(config["output_table_dir"])
    table_dir.mkdir(parents=True, exist_ok=True)
    partition_preview = labels[["source_gse_id", "truth_class"]].copy()
    partition_preview["partition"] = partition
    partition_preview.groupby(["partition", "source_gse_id", "truth_class"], observed=True).size().reset_index(name="n_cells").to_csv(
        table_dir / "development_partition_counts.csv", index=False
    )
    if args.partition_only:
        print(json.dumps({"status": "PASS_PARTITION_ONLY", "counts": partition.value_counts().to_dict()}, indent=2))
        return

    all_names = features.gene.astype(str).tolist()
    stage1_columns, stage2_columns = feature_columns(features, config)
    final_bundle = train_bundle(
        labels, matrix, all_names, stage1_columns, stage2_columns, partition, config
    )
    holdout_tables = []
    if not args.skip_source_holdout:
        for source in ("HRA005041", "GSE144469", "MalteGDT"):
            outer_partition = partition.copy()
            if source == "GSE144469":
                # HRA fold 2 has one indivisible positive group. Move the independent
                # Malte group from fitting to calibration so the GSE144469 outer run
                # retains disjoint positive evidence for fitting, calibration, and
                # threshold validation without splitting a biological group.
                malte_fit = labels.source_gse_id.eq("MalteGDT") & outer_partition.eq("fit")
                outer_partition.loc[malte_fit] = "calibration"
            outer_bundle = train_bundle(
                labels, matrix, all_names, stage1_columns, stage2_columns, outer_partition, config,
                excluded_source=source,
            )
            holdout_tables.append(
                evaluate_source_holdout(
                    outer_bundle, labels, matrix, all_names, stage1_columns, stage2_columns, source
                )
            )
            del outer_bundle
            cp.get_default_memory_pool().free_all_blocks()
    source_holdout = pd.concat(holdout_tables, ignore_index=True) if holdout_tables else pd.DataFrame()
    summary = save_bundle(
        final_bundle, labels, features, partition, config, config_path, cache_manifest,
        source_holdout, started,
    )
    print(json.dumps(json_safe(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
