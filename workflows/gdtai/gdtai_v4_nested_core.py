"""Pure modeling utilities for the precommitted gdTAI V4 evaluation.

This module contains no project-path or H5AD mutation logic. The command-line
workflow owns data access and checkpointing; these functions keep weighting,
fold-local filtering, calibration, operating-point selection, and statistics
small enough to test on synthetic data before Step 2 is authorized.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from dataclasses import dataclass
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.special import expit, logit
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.preprocessing import StandardScaler


DERIVED_FEATURE_ORDER = (
    "TRA_detected",
    "TRB_detected",
    "TRG_detected",
    "TRD_detected",
    "T_lineage_mean",
    "NK_context_mean",
    "CD4_helper_mean",
    "Treg_mean",
)


class GuardrailFailure(RuntimeError):
    """Raised when no operating threshold satisfies a frozen mode."""


def canonical_json_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def stable_seed(base_seed: int, *tokens: object) -> int:
    payload = "|".join([str(base_seed), *(str(token) for token in tokens)])
    digest = hashlib.sha256(payload.encode("utf-8")).digest()
    return int.from_bytes(digest[:4], "little", signed=False)


def parameter_grid(spec: Mapping[str, Sequence[Any]]) -> list[dict[str, Any]]:
    keys = sorted(spec)
    return [dict(zip(keys, values, strict=True)) for values in itertools.product(*(spec[key] for key in keys))]


def candidate_id(family: str, parameters: Mapping[str, Any]) -> str:
    suffix = "__".join(f"{key}={parameters[key]}" for key in sorted(parameters))
    return f"{family}__{suffix}" if suffix else family


def family_columns(feature_names: Sequence[str], prefixes: Sequence[str]) -> np.ndarray:
    return np.asarray(
        [index for index, gene in enumerate(feature_names) if any(gene.startswith(prefix) for prefix in prefixes)],
        dtype=np.int64,
    )


def named_columns(feature_names: Sequence[str], names: Sequence[str]) -> np.ndarray:
    lookup = {name: index for index, name in enumerate(feature_names)}
    missing = sorted(set(names) - set(lookup))
    if missing:
        raise ValueError(f"Missing frozen features: {', '.join(missing)}")
    return np.asarray([lookup[name] for name in names], dtype=np.int64)


def derive_features(
    gene_matrix: np.ndarray,
    feature_names: Sequence[str],
    family_prefixes: Mapping[str, Sequence[str]],
    panels: Mapping[str, Sequence[str]],
) -> tuple[np.ndarray, list[str]]:
    """Build the eight allowed deterministic features from log1p(CP10K)."""
    x = np.asarray(gene_matrix, dtype=np.float32)
    if x.ndim != 2 or x.shape[1] != len(feature_names):
        raise ValueError("Gene matrix shape does not match the frozen feature order")
    columns: list[np.ndarray] = []
    output_names: list[str] = []
    for name in DERIVED_FEATURE_ORDER[:4]:
        if name not in family_prefixes:
            raise ValueError(f"Missing family-prefix definition for {name}")
        indices = family_columns(feature_names, family_prefixes[name])
        if indices.size == 0:
            raise ValueError(f"No frozen genes match the {name} family")
        columns.append((x[:, indices] > 0).sum(axis=1, dtype=np.int32).astype(np.float32))
        output_names.append(name)
    for name in DERIVED_FEATURE_ORDER[4:]:
        if name not in panels:
            raise ValueError(f"Missing panel definition for {name}")
        indices = named_columns(feature_names, panels[name])
        columns.append(x[:, indices].mean(axis=1, dtype=np.float32))
        output_names.append(name)
    return np.column_stack(columns).astype(np.float32, copy=False), output_names


def exclusion_flags(
    gene_matrix: np.ndarray,
    feature_names: Sequence[str],
    cd4_rule: Mapping[str, Any],
    treg_rule: Mapping[str, Any],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply the immutable post-model CD4-helper and Treg exclusions."""
    x = np.asarray(gene_matrix, dtype=np.float32)
    lookup = {name: index for index, name in enumerate(feature_names)}

    def rule_mask(rule: Mapping[str, Any]) -> np.ndarray:
        required = [rule["marker"], *rule["support_genes"], *rule["program_genes"]]
        missing = sorted(set(required) - set(lookup))
        if missing:
            raise ValueError(f"Exclusion rule is missing genes: {', '.join(missing)}")
        marker = x[:, lookup[str(rule["marker"])]] > float(rule["marker_min"])
        support_idx = np.asarray([lookup[str(gene)] for gene in rule["support_genes"]], dtype=np.int64)
        program_idx = np.asarray([lookup[str(gene)] for gene in rule["program_genes"]], dtype=np.int64)
        support = (x[:, support_idx] > 0).sum(axis=1) >= int(rule["support_detected_min"])
        program = x[:, program_idx].mean(axis=1) >= float(rule["program_mean_min"])
        return marker & support & program

    cd4 = rule_mask(cd4_rule)
    treg = rule_mask(treg_rule)
    return cd4, treg, cd4 | treg


def inference_coverage(
    available_genes: Iterable[str],
    retained_genes: Sequence[str],
    critical_genes: Sequence[str],
    minimum_coverage: float,
) -> dict[str, Any]:
    available = set(available_genes)
    retained = list(retained_genes)
    missing_critical = sorted(set(critical_genes) - available)
    present = sum(gene in available for gene in retained)
    coverage = present / len(retained) if retained else 0.0
    abstain = bool(missing_critical or coverage < minimum_coverage)
    return {
        "retained_gene_count": len(retained),
        "present_gene_count": present,
        "coverage": coverage,
        "missing_critical_genes": missing_critical,
        "abstain": abstain,
    }


def fold_feature_mask(
    matrix: np.ndarray,
    gene_feature_count: int,
    minimum_detection_fraction: float,
    maximum_retained_genes: int,
) -> np.ndarray:
    """Return a fold-local mask for genes and deterministic non-gene features."""
    x = np.asarray(matrix)
    if x.ndim != 2 or not 0 <= gene_feature_count <= x.shape[1]:
        raise ValueError("Invalid feature matrix or gene feature count")
    if gene_feature_count > maximum_retained_genes:
        raise ValueError("Candidate gene set exceeds the frozen maximum")
    keep = np.zeros(x.shape[1], dtype=bool)
    if gene_feature_count:
        genes = x[:, :gene_feature_count]
        minimum_detected = max(1, int(math.ceil(x.shape[0] * minimum_detection_fraction)))
        detected = np.count_nonzero(genes > 0, axis=0)
        finite = np.isfinite(genes).all(axis=0)
        nonconstant = np.ptp(genes, axis=0) > 0
        keep[:gene_feature_count] = finite & nonconstant & (detected >= minimum_detected)
    if gene_feature_count < x.shape[1]:
        derived = x[:, gene_feature_count:]
        keep[gene_feature_count:] = np.isfinite(derived).all(axis=0) & (np.ptp(derived, axis=0) > 0)
    if int(keep[:gene_feature_count].sum()) > maximum_retained_genes:
        raise RuntimeError("Fold-local filter retained too many genes")
    if not keep.any():
        raise RuntimeError("Fold-local filtering removed every candidate feature")
    return keep


def balanced_training_weights(
    labels: np.ndarray,
    sources: np.ndarray,
    reliability: np.ndarray,
) -> np.ndarray:
    """Apply the frozen reliability x dataset x within-dataset class weights."""
    y = np.asarray(labels, dtype=np.int8)
    source = np.asarray(sources, dtype=object)
    rel = np.asarray(reliability, dtype=np.float64)
    if not (y.shape == source.shape == rel.shape):
        raise ValueError("Labels, sources, and reliability must have identical shapes")
    if y.size == 0 or set(np.unique(y)) - {0, 1}:
        raise ValueError("Training labels must be nonempty binary values")
    if not np.isfinite(rel).all() or (rel < 0).any():
        raise ValueError("Reliability weights must be finite and nonnegative")
    frame = pd.DataFrame({"source": source, "label": y})
    cell_counts = frame.groupby(["source", "label"], observed=True).size().to_dict()
    class_counts = frame.groupby("label", observed=True).size().to_dict()
    source_counts = frame.groupby("source", observed=True).size().to_dict()
    sources_per_class = frame.groupby("label", observed=True)["source"].nunique().to_dict()
    classes_per_source = frame.groupby("source", observed=True)["label"].nunique().to_dict()
    weights = np.empty(y.size, dtype=np.float64)
    for index, (src, label) in enumerate(zip(source, y, strict=True)):
        pair_count = float(cell_counts[(src, int(label))])
        dataset_within_class = float(class_counts[int(label)]) / (
            float(sources_per_class[int(label)]) * pair_count
        )
        class_within_dataset = float(source_counts[src]) / (
            float(classes_per_source[src]) * pair_count
        )
        weights[index] = rel[index] * dataset_within_class * class_within_dataset
    positive = weights > 0
    if not positive.any():
        raise ValueError("All training weights are zero")
    weights /= weights[positive].mean()
    return weights.astype(np.float32)


@dataclass
class FittedBinaryEstimator:
    family: str
    parameters: dict[str, Any]
    selected_columns: np.ndarray
    selected_feature_names: list[str]
    scaler: StandardScaler | None
    estimator: Any
    n_iter: int
    converged: bool
    convergence_applicable: bool

    def predict_probability(self, matrix: np.ndarray) -> np.ndarray:
        x = np.asarray(matrix)[:, self.selected_columns]
        if self.scaler is not None:
            x = self.scaler.transform(x)
        return self.estimator.predict_proba(x)[:, 1].astype(np.float64)


def fit_binary_estimator(
    matrix: np.ndarray,
    labels: np.ndarray,
    sample_weight: np.ndarray,
    feature_names: Sequence[str],
    gene_feature_count: int,
    family: str,
    parameters: Mapping[str, Any],
    minimum_detection_fraction: float,
    maximum_retained_genes: int,
    seed: int,
) -> FittedBinaryEstimator:
    x = np.asarray(matrix, dtype=np.float32)
    y = np.asarray(labels, dtype=np.int8)
    weight = np.asarray(sample_weight, dtype=np.float32)
    if len(feature_names) != x.shape[1]:
        raise ValueError("Feature names do not match matrix columns")
    if np.unique(y).size != 2:
        raise ValueError("Estimator training requires both classes")
    mask = fold_feature_mask(x, gene_feature_count, minimum_detection_fraction, maximum_retained_genes)
    selected = np.flatnonzero(mask)
    selected_x = x[:, selected]
    scaler: StandardScaler | None = None
    if family in {"elastic_net", "compact_logistic", "tcr_logistic"}:
        scaler = StandardScaler()
        selected_x = scaler.fit_transform(selected_x, sample_weight=weight)
        estimator = LogisticRegression(
            C=float(parameters["C"]),
            l1_ratio=float(parameters["l1_ratio"]),
            penalty="elasticnet",
            solver="saga",
            max_iter=int(parameters["max_iter"]),
            tol=float(parameters["tolerance"]),
            random_state=int(seed),
            n_jobs=1,
        )
    elif family == "hist_gradient_boosting":
        estimator = HistGradientBoostingClassifier(
            learning_rate=float(parameters["learning_rate"]),
            max_leaf_nodes=int(parameters["max_leaf_nodes"]),
            min_samples_leaf=int(parameters["min_samples_leaf"]),
            l2_regularization=float(parameters["l2_regularization"]),
            max_iter=int(parameters["max_iter"]),
            early_stopping=bool(parameters["early_stopping"]),
            random_state=int(seed),
        )
    else:
        raise ValueError(f"Unknown model family: {family}")
    estimator.fit(selected_x, y, sample_weight=weight)
    if family in {"elastic_net", "compact_logistic", "tcr_logistic"}:
        n_iter = int(np.max(estimator.n_iter_))
        converged = n_iter < int(parameters["max_iter"])
        convergence_applicable = True
    else:
        n_iter = int(estimator.n_iter_)
        converged = True
        convergence_applicable = False
    return FittedBinaryEstimator(
        family=family,
        parameters=dict(parameters),
        selected_columns=selected,
        selected_feature_names=[str(feature_names[index]) for index in selected],
        scaler=scaler,
        estimator=estimator,
        n_iter=n_iter,
        converged=converged,
        convergence_applicable=convergence_applicable,
    )


@dataclass(frozen=True)
class PlattCalibrator:
    intercept: float
    slope: float

    def predict(self, raw_probability: np.ndarray) -> np.ndarray:
        raw = np.clip(np.asarray(raw_probability, dtype=np.float64), 1e-7, 1 - 1e-7)
        return expit(self.intercept + self.slope * logit(raw))


def fit_platt_calibrator(
    raw_probability: np.ndarray,
    labels: np.ndarray,
    sample_weight: np.ndarray,
    seed: int,
) -> PlattCalibrator:
    raw = np.clip(np.asarray(raw_probability, dtype=np.float64), 1e-7, 1 - 1e-7)
    y = np.asarray(labels, dtype=np.int8)
    weight = np.asarray(sample_weight, dtype=np.float64)
    if np.unique(y).size != 2:
        raise ValueError("Calibration requires both classes")
    model = LogisticRegression(
        penalty=None,
        solver="lbfgs",
        max_iter=1000,
        tol=1e-10,
        random_state=int(seed),
    )
    model.fit(logit(raw).reshape(-1, 1), y, sample_weight=weight)
    return PlattCalibrator(intercept=float(model.intercept_[0]), slope=float(model.coef_[0, 0]))


def threshold_candidates(*score_arrays: np.ndarray) -> np.ndarray:
    values = [np.asarray(scores, dtype=np.float64) for scores in score_arrays if np.asarray(scores).size]
    if not values:
        return np.asarray([0.0, 1.0], dtype=np.float64)
    combined = np.concatenate([*values, np.asarray([0.0, 1.0], dtype=np.float64)])
    combined = combined[np.isfinite(combined)]
    return np.unique(np.clip(combined, 0.0, 1.0))


def count_greater_equal(values: np.ndarray, thresholds: np.ndarray) -> np.ndarray:
    sorted_values = np.sort(np.asarray(values, dtype=np.float64))
    return sorted_values.size - np.searchsorted(sorted_values, thresholds, side="left")


def f_beta_from_counts(tp: np.ndarray, fp: np.ndarray, fn: np.ndarray, beta: float) -> np.ndarray:
    beta2 = beta * beta
    denominator = (1 + beta2) * tp + beta2 * fn + fp
    return np.divide((1 + beta2) * tp, denominator, out=np.zeros_like(tp, dtype=float), where=denominator > 0)


@dataclass(frozen=True)
class ThresholdDecision:
    mode: str
    passed: bool
    threshold: float
    objective_value: float
    macro_precision: float
    macro_recall: float
    macro_f1: float
    macro_f0_5: float
    paired_abt_fpr: float
    strict_nk_fpr: float
    candidate_count: int
    valid_candidate_count: int
    candidate_sha256: str
    reason: str


def select_stage1_threshold(
    calibrated_probability: np.ndarray,
    sources: np.ndarray,
    gdt_mask: np.ndarray,
    abt_mask: np.ndarray,
    strict_nk_mask: np.ndarray,
    primary_sources: Sequence[str],
    minimum_gdt_recall: float,
    minimum_abt_recall: float,
) -> ThresholdDecision:
    score = np.asarray(calibrated_probability, dtype=np.float64)
    source = np.asarray(sources, dtype=object)
    gdt = np.asarray(gdt_mask, dtype=bool)
    abt = np.asarray(abt_mask, dtype=bool)
    nk = np.asarray(strict_nk_mask, dtype=bool)
    thresholds = threshold_candidates(score)
    valid = np.ones(thresholds.size, dtype=bool)
    recalls_gdt: list[np.ndarray] = []
    recalls_abt: list[np.ndarray] = []
    for primary_source in primary_sources:
        source_gdt = score[gdt & (source == primary_source)]
        source_abt = score[abt & (source == primary_source)]
        if source_gdt.size == 0 or source_abt.size == 0:
            return ThresholdDecision(
                "stage1",
                False,
                math.nan,
                math.nan,
                math.nan,
                math.nan,
                math.nan,
                math.nan,
                math.nan,
                math.nan,
                int(thresholds.size),
                0,
                hashlib.sha256(thresholds.tobytes()).hexdigest(),
                f"Source {primary_source} lacks a Stage-1 recall stratum",
            )
        gdt_recall = count_greater_equal(source_gdt, thresholds) / source_gdt.size
        abt_recall = count_greater_equal(source_abt, thresholds) / source_abt.size
        valid &= (gdt_recall >= minimum_gdt_recall) & (abt_recall >= minimum_abt_recall)
        recalls_gdt.append(gdt_recall)
        recalls_abt.append(abt_recall)
    nk_total = int(nk.sum())
    nk_fpr = count_greater_equal(score[nk], thresholds) / nk_total if nk_total else np.full(thresholds.size, np.nan)
    valid_indices = np.flatnonzero(valid)
    digest = hashlib.sha256(thresholds.tobytes()).hexdigest()
    if valid_indices.size == 0:
        return ThresholdDecision(
            "stage1",
            False,
            math.nan,
            math.nan,
            math.nan,
            math.nan,
            math.nan,
            math.nan,
            math.nan,
            math.nan,
            int(thresholds.size),
            0,
            digest,
            "No threshold satisfies every per-source Stage-1 recall guardrail",
        )
    # Highest feasible threshold is deterministic and minimizes NK calls for a
    # monotone score while retaining the required T-cell recalls.
    best = int(valid_indices[-1])
    macro_gdt = float(np.mean([values[best] for values in recalls_gdt]))
    macro_abt = float(np.mean([values[best] for values in recalls_abt]))
    return ThresholdDecision(
        "stage1",
        True,
        float(thresholds[best]),
        float(1.0 - nk_fpr[best]) if nk_total else math.nan,
        math.nan,
        float((macro_gdt + macro_abt) / 2),
        math.nan,
        math.nan,
        math.nan,
        float(nk_fpr[best]) if nk_total else math.nan,
        int(thresholds.size),
        int(valid_indices.size),
        digest,
        "PASS",
    )


def select_stage2_threshold(
    calibrated_primary_probability: np.ndarray,
    labels: np.ndarray,
    sources: np.ndarray,
    paired_abt_mask: np.ndarray,
    calibrated_nk_probability: np.ndarray,
    stage1_primary_probability: np.ndarray,
    stage1_nk_probability: np.ndarray,
    stage1_threshold: float,
    primary_exclusion: np.ndarray,
    nk_exclusion: np.ndarray,
    mode: str,
    mode_spec: Mapping[str, Any],
) -> ThresholdDecision:
    score = np.asarray(calibrated_primary_probability, dtype=np.float64)
    y = np.asarray(labels, dtype=np.int8)
    source = np.asarray(sources, dtype=object)
    paired = np.asarray(paired_abt_mask, dtype=bool)
    nk_score = np.asarray(calibrated_nk_probability, dtype=np.float64)
    p1 = np.asarray(stage1_primary_probability, dtype=np.float64)
    p1_nk = np.asarray(stage1_nk_probability, dtype=np.float64)
    excluded = np.asarray(primary_exclusion, dtype=bool)
    nk_excluded = np.asarray(nk_exclusion, dtype=bool)
    if not (score.shape == y.shape == source.shape == paired.shape == p1.shape == excluded.shape):
        raise ValueError("Primary threshold arrays have inconsistent shapes")
    if not (nk_score.shape == p1_nk.shape == nk_excluded.shape):
        raise ValueError("NK threshold arrays have inconsistent shapes")
    primary_effective = np.where((p1 >= stage1_threshold) & (~excluded), score, -np.inf)
    nk_effective = np.where((p1_nk >= stage1_threshold) & (~nk_excluded), nk_score, -np.inf)
    thresholds = threshold_candidates(score, nk_score)
    source_values = sorted(pd.unique(source).tolist())
    precision_by_source: list[np.ndarray] = []
    recall_by_source: list[np.ndarray] = []
    f1_by_source: list[np.ndarray] = []
    f05_by_source: list[np.ndarray] = []
    for source_value in source_values:
        source_mask = source == source_value
        pos_values = primary_effective[source_mask & (y == 1)]
        neg_values = primary_effective[source_mask & (y == 0)]
        if pos_values.size == 0 or neg_values.size == 0:
            raise ValueError(f"Source {source_value} lacks a Stage-2 class")
        tp = count_greater_equal(pos_values, thresholds).astype(float)
        fp = count_greater_equal(neg_values, thresholds).astype(float)
        fn = pos_values.size - tp
        precision = np.divide(tp, tp + fp, out=np.zeros_like(tp), where=(tp + fp) > 0)
        recall = tp / pos_values.size
        precision_by_source.append(precision)
        recall_by_source.append(recall)
        f1_by_source.append(f_beta_from_counts(tp, fp, fn, beta=1.0))
        f05_by_source.append(f_beta_from_counts(tp, fp, fn, beta=0.5))
    macro_precision = np.mean(np.vstack(precision_by_source), axis=0)
    macro_recall = np.mean(np.vstack(recall_by_source), axis=0)
    macro_f1 = np.mean(np.vstack(f1_by_source), axis=0)
    macro_f05 = np.mean(np.vstack(f05_by_source), axis=0)
    paired_total = int(paired.sum())
    paired_fpr = (
        count_greater_equal(primary_effective[paired], thresholds) / paired_total
        if paired_total
        else np.full(thresholds.size, np.nan)
    )
    nk_total = int(nk_effective.size)
    nk_fpr = (
        count_greater_equal(nk_effective, thresholds) / nk_total
        if nk_total
        else np.full(thresholds.size, np.nan)
    )
    valid = (
        (macro_recall >= float(mode_spec["minimum_macro_recall"]))
        & (paired_fpr <= float(mode_spec["maximum_paired_abt_fpr"]))
        & (nk_fpr <= float(mode_spec["maximum_strict_nk_fpr"]))
    )
    valid_indices = np.flatnonzero(valid)
    digest = hashlib.sha256(thresholds.tobytes()).hexdigest()
    if valid_indices.size == 0:
        return ThresholdDecision(
            mode,
            False,
            math.nan,
            math.nan,
            math.nan,
            math.nan,
            math.nan,
            math.nan,
            math.nan,
            math.nan,
            int(thresholds.size),
            0,
            digest,
            "No threshold satisfies the frozen operating-mode guardrails",
        )
    objective_name = str(mode_spec["objective"])
    objective = macro_f1 if objective_name == "f1" else macro_f05
    valid_objective = objective[valid_indices]
    best_value = float(valid_objective.max())
    tied = valid_indices[np.isclose(valid_objective, best_value, rtol=0.0, atol=1e-12)]
    # Apply the frozen threshold tie order: strict-NK FPR, paired-abT FPR,
    # then a conservative higher threshold when predictions are equivalent.
    best = int(
        sorted(
            tied.tolist(),
            key=lambda index: (float(nk_fpr[index]), float(paired_fpr[index]), -float(thresholds[index])),
        )[0]
    )
    return ThresholdDecision(
        mode,
        True,
        float(thresholds[best]),
        float(objective[best]),
        float(macro_precision[best]),
        float(macro_recall[best]),
        float(macro_f1[best]),
        float(macro_f05[best]),
        float(paired_fpr[best]),
        float(nk_fpr[best]),
        int(thresholds.size),
        int(valid_indices.size),
        digest,
        "PASS",
    )


def apply_two_stage_call(
    stage1_probability: np.ndarray,
    stage2_probability: np.ndarray,
    stage1_threshold: float,
    stage2_threshold: float,
    exclusion: np.ndarray,
) -> np.ndarray:
    return (
        (np.asarray(stage1_probability) >= stage1_threshold)
        & (np.asarray(stage2_probability) >= stage2_threshold)
        & (~np.asarray(exclusion, dtype=bool))
    )


def vdj_rescue_mask(
    rna_call: np.ndarray,
    exclusion: np.ndarray,
    paired_gdt: np.ndarray,
    any_abt: np.ndarray,
) -> np.ndarray:
    """Return optional rescue calls without changing the RNA metric frame."""
    call = np.asarray(rna_call, dtype=bool)
    return (~call) & np.asarray(exclusion, dtype=bool) & np.asarray(paired_gdt, dtype=bool) & (~np.asarray(any_abt, dtype=bool))


def binary_metrics(labels: np.ndarray, prediction: np.ndarray, probability: np.ndarray) -> dict[str, Any]:
    y = np.asarray(labels, dtype=np.int8)
    pred = np.asarray(prediction, dtype=bool)
    score = np.asarray(probability, dtype=np.float64)
    tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel().astype(int)
    specificity = tn / (tn + fp) if tn + fp else math.nan
    output: dict[str, Any] = {
        "n_cells": int(y.size),
        "n_positive": int((y == 1).sum()),
        "n_negative": int((y == 0).sum()),
        "predicted_positive": int(pred.sum()),
        "tp": int(tp),
        "fp": int(fp),
        "tn": int(tn),
        "fn": int(fn),
        "precision": float(precision_score(y, pred, zero_division=0)),
        "recall": float(recall_score(y, pred, zero_division=0)),
        "specificity": float(specificity),
        "f1": float(f1_score(y, pred, zero_division=0)),
        "f0_5": float(f_beta_from_counts(np.asarray([tp]), np.asarray([fp]), np.asarray([fn]), 0.5)[0]),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)) if np.unique(y).size == 2 else math.nan,
        "mcc": float(matthews_corrcoef(y, pred)) if np.unique(y).size == 2 else math.nan,
        "roc_auc": float(roc_auc_score(y, score)) if np.unique(y).size == 2 else math.nan,
        "pr_auc": float(average_precision_score(y, score)) if np.unique(y).size == 2 else math.nan,
        "brier_score": float(brier_score_loss(y, score)),
    }
    output.update(calibration_metrics(y, score))
    return output


def calibration_metrics(labels: np.ndarray, probability: np.ndarray, bins: int = 10) -> dict[str, float]:
    y = np.asarray(labels, dtype=np.int8)
    score = np.clip(np.asarray(probability, dtype=np.float64), 1e-7, 1 - 1e-7)
    if np.unique(y).size != 2:
        return {"calibration_intercept": math.nan, "calibration_slope": math.nan, "ece": math.nan}
    model = LogisticRegression(penalty=None, solver="lbfgs", max_iter=1000, tol=1e-10)
    model.fit(logit(score).reshape(-1, 1), y)
    edges = np.linspace(0.0, 1.0, bins + 1)
    assignments = np.clip(np.digitize(score, edges[1:-1], right=False), 0, bins - 1)
    ece = 0.0
    for index in range(bins):
        mask = assignments == index
        if mask.any():
            ece += float(mask.mean()) * abs(float(score[mask].mean()) - float(y[mask].mean()))
    return {
        "calibration_intercept": float(model.intercept_[0]),
        "calibration_slope": float(model.coef_[0, 0]),
        "ece": float(ece),
    }


def prevalence_ppv_table(sensitivity: float, specificity: float, prevalences: Sequence[float]) -> pd.DataFrame:
    rows: list[dict[str, float]] = []
    fpr = 1.0 - specificity
    for prevalence in prevalences:
        numerator = sensitivity * prevalence
        denominator = numerator + fpr * (1.0 - prevalence)
        ppv = numerator / denominator if denominator > 0 else math.nan
        rows.append({"prevalence": float(prevalence), "ppv": float(ppv), "fdr": float(1.0 - ppv)})
    return pd.DataFrame(rows)


def grouped_confusion_counts(
    labels: np.ndarray,
    predictions: Mapping[str, np.ndarray],
    sources: np.ndarray,
    groups: np.ndarray,
) -> pd.DataFrame:
    frame = pd.DataFrame({"source": sources, "group": groups, "y": np.asarray(labels, dtype=np.int8)})
    for name, prediction in predictions.items():
        pred = np.asarray(prediction, dtype=bool)
        frame[f"{name}_tp"] = (frame["y"].to_numpy() == 1) & pred
        frame[f"{name}_fp"] = (frame["y"].to_numpy() == 0) & pred
        frame[f"{name}_fn"] = (frame["y"].to_numpy() == 1) & (~pred)
    aggregate_columns = [column for column in frame if column.endswith(("_tp", "_fp", "_fn"))]
    return frame.groupby(["source", "group"], observed=True, sort=True)[aggregate_columns].sum().reset_index()


def paired_hierarchical_bootstrap_f1_difference(
    grouped_counts: pd.DataFrame,
    candidate: str,
    comparator: str,
    replicates: int,
    seed: int,
) -> np.ndarray:
    """Resample sources, then groups, and return macro-F1 differences."""
    rng = np.random.default_rng(seed)
    source_values = sorted(grouped_counts["source"].astype(str).unique().tolist())
    by_source = {
        source: grouped_counts[grouped_counts["source"].astype(str).eq(source)].reset_index(drop=True)
        for source in source_values
    }
    output = np.empty(replicates, dtype=np.float64)

    def f1_from_rows(frame: pd.DataFrame, indices: np.ndarray, prefix: str) -> float:
        tp = float(frame.loc[indices, f"{prefix}_tp"].sum())
        fp = float(frame.loc[indices, f"{prefix}_fp"].sum())
        fn = float(frame.loc[indices, f"{prefix}_fn"].sum())
        denominator = 2 * tp + fp + fn
        return 2 * tp / denominator if denominator > 0 else 0.0

    for replicate in range(replicates):
        sampled_sources = rng.choice(source_values, size=len(source_values), replace=True)
        differences: list[float] = []
        for source in sampled_sources:
            frame = by_source[str(source)]
            indices = rng.integers(0, frame.shape[0], size=frame.shape[0])
            differences.append(f1_from_rows(frame, indices, candidate) - f1_from_rows(frame, indices, comparator))
        output[replicate] = float(np.mean(differences))
    return output


def bootstrap_interval(values: np.ndarray, confidence: float = 0.95) -> tuple[float, float]:
    alpha = (1.0 - confidence) / 2.0
    return tuple(float(value) for value in np.quantile(np.asarray(values, dtype=float), [alpha, 1.0 - alpha]))
