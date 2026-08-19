#!/usr/bin/env python3
"""Train the precommitted gdTAI V4.2 two-stage GPU XGBoost cascade."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import time
from pathlib import Path

import cupy as cp
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "configs/models/gdtai/v4_2_training.json"
POSITIVE_SOURCES = ["HRA005041", "GSE144469", "MalteGDT"]
NK_SOURCES = ["GSE125527", "GSE228597", "GSE212217", "GSE234069", "GSE235863", "GSE243013", "GSE287541"]


def resolve(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def sha256_file(path: Path, chunk: int = 64 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def stable_seed(*parts: object) -> int:
    payload = "::".join(map(str, parts)).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:4], "little") & 0x7FFFFFFF


def sampled_rows(frame: pd.DataFrame, candidates: np.ndarray, maximum: int, seed: int) -> np.ndarray:
    candidates = np.asarray(candidates, dtype=np.int64)
    if len(candidates) <= maximum:
        return candidates
    local = frame.iloc[candidates]
    chosen = []
    groups = list(local.groupby("source_gse_id", observed=True).groups.items())
    per_source = max(1, maximum // max(1, len(groups)))
    for source, indexes in groups:
        positions = frame.index.get_indexer(np.asarray(indexes))
        rng = np.random.default_rng(stable_seed(seed, source))
        chosen.extend(rng.choice(positions, size=min(per_source, len(positions)), replace=False))
    chosen = np.unique(np.asarray(chosen, dtype=np.int64))
    if len(chosen) < maximum:
        remaining = np.setdiff1d(candidates, chosen, assume_unique=False)
        rng = np.random.default_rng(seed)
        chosen = np.concatenate([chosen, rng.choice(remaining, size=min(maximum - len(chosen), len(remaining)), replace=False)])
    return np.sort(chosen[:maximum])


def source_balanced_weights(frame: pd.DataFrame, rows: np.ndarray, stage: int) -> np.ndarray:
    local = frame.iloc[rows]
    if stage == 1:
        classes = np.where(local.truth_class.str.startswith("nk_"), "NK", "T")
    else:
        classes = np.where(local.truth_class.eq("gdT_gold"), "gdT", "non_gdT")
    weights = local.reliability_weight.to_numpy(float)
    key = pd.DataFrame({"class": classes, "source": local.source_gse_id.to_numpy(), "base": weights})
    totals = key.groupby(["class", "source"], observed=True).base.transform("sum").to_numpy(float)
    source_counts = key.groupby("class", observed=True).source.transform("nunique").to_numpy(float)
    weights = np.divide(weights, totals * source_counts, out=np.zeros_like(weights), where=(totals * source_counts) > 0)
    # Equalize the two top-level classes.
    class_totals = pd.Series(weights).groupby(classes).transform("sum").to_numpy(float)
    weights = np.divide(weights, class_totals, out=np.zeros_like(weights), where=class_totals > 0)
    if stage == 1:
        tier1 = local.truth_class.eq("nk_tier1").to_numpy()
        tier2 = local.truth_class.eq("nk_tier2").to_numpy()
        if tier1.any() and tier2.any() and weights[tier2].sum() > weights[tier1].sum():
            weights[tier2] *= weights[tier1].sum() / weights[tier2].sum()
    return (weights / weights.mean()).astype(np.float32)


def fit_xgb(x: np.ndarray, y: np.ndarray, weight: np.ndarray, params: dict, fixed: dict,
            seed: int, eval_data: tuple[np.ndarray, np.ndarray] | None = None) -> xgb.XGBClassifier:
    spec = dict(fixed) | dict(params) | {"random_state": seed}
    if eval_data is None:
        spec.pop("early_stopping_rounds", None)
    model = xgb.XGBClassifier(**spec)
    kwargs = {"sample_weight": cp.asarray(weight), "verbose": False}
    if eval_data is not None:
        kwargs["eval_set"] = [(cp.asarray(eval_data[0]), cp.asarray(eval_data[1]))]
    model.fit(cp.asarray(x), cp.asarray(y), **kwargs)
    return model


def predict(model: xgb.XGBClassifier, x: np.ndarray, chunk: int = 250_000) -> np.ndarray:
    output = []
    for start in range(0, len(x), chunk):
        value = model.predict_proba(cp.asarray(np.asarray(x[start:start + chunk], dtype=np.float32)))[:, 1]
        output.append(cp.asnumpy(value) if isinstance(value, cp.ndarray) else np.asarray(value))
    return np.concatenate(output).astype(np.float32) if output else np.asarray([], dtype=np.float32)


def exclusion_flags(x: np.ndarray, names: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lookup = {name: i for i, name in enumerate(names)}
    def values(gene: str) -> np.ndarray:
        return np.asarray(x[:, lookup[gene]])
    cd4_genes = ["CD4", "IL7R", "CCR7", "SELL", "LTB", "MAL"]
    treg_genes = ["FOXP3", "IL2RA", "CTLA4", "TIGIT"]
    cd4 = (values("CD4") > 0) & (np.column_stack([values(g) > 0 for g in cd4_genes[1:]]).sum(1) >= 2) & (np.column_stack([values(g) for g in cd4_genes]).mean(1) >= 0.5)
    treg = (values("FOXP3") > 0) & (np.column_stack([values(g) > 0 for g in treg_genes[1:]]).sum(1) >= 1) & (np.column_stack([values(g) for g in treg_genes]).mean(1) >= 0.5)
    return cd4, treg, cd4 | treg


def receptor_rescue_flags(x: np.ndarray, names: list[str], thresholds: dict[str, float] | None) -> np.ndarray:
    if not thresholds:
        return np.zeros(len(x), dtype=bool)
    lookup = {name: index for index, name in enumerate(names)}
    missing = sorted(set(thresholds) - set(lookup))
    if missing:
        raise RuntimeError(f"Receptor-rescue genes are missing from the frozen feature cache: {missing}")
    return np.logical_or.reduce([
        np.asarray(x[:, lookup[gene]]) >= float(threshold)
        for gene, threshold in thresholds.items()
    ])


def fbeta(tp: int, fp: int, fn: int, beta: float) -> float:
    b2 = beta * beta
    denominator = (1 + b2) * tp + b2 * fn + fp
    return (1 + b2) * tp / denominator if denominator else 0.0


def source_rates(frame: pd.DataFrame, calls: np.ndarray, truth: str) -> dict[str, float]:
    mask = frame.truth_class.eq(truth).to_numpy()
    return {source: float(calls[mask & frame.source_gse_id.eq(source).to_numpy()].mean()) for source in sorted(frame.loc[mask, "source_gse_id"].unique())}


def sorted_scores_by_source(frame: pd.DataFrame, score: np.ndarray, truth: str,
                            base: np.ndarray | None = None) -> tuple[dict[str, np.ndarray], dict[str, int]]:
    truth_mask = frame.truth_class.eq(truth).to_numpy()
    if base is None:
        base = np.ones(len(frame), dtype=bool)
    arrays, totals = {}, {}
    sources = frame.source_gse_id.to_numpy()
    for source in sorted(set(sources[truth_mask])):
        group = truth_mask & (sources == source)
        totals[source] = int(group.sum())
        arrays[source] = np.sort(score[group & base])
    return arrays, totals


def rates_from_sorted(arrays: dict[str, np.ndarray], totals: dict[str, int], threshold: float) -> dict[str, float]:
    return {
        source: float((len(values) - np.searchsorted(values, threshold, side="left")) / totals[source])
        for source, values in arrays.items()
    }


def choose_stage1_threshold(frame: pd.DataFrame, score: np.ndarray, spec: dict) -> tuple[float, dict]:
    thresholds = np.unique(np.quantile(score, np.linspace(0, 1, 2001)))
    gd_arrays, gd_totals = sorted_scores_by_source(frame, score, "gdT_gold")
    ab_arrays, ab_totals = sorted_scores_by_source(frame, score, "abT_gold")
    nk_arrays, nk_totals = sorted_scores_by_source(frame, score, "nk_tier1")
    best = None
    for threshold in thresholds[::-1]:
        gd = rates_from_sorted(gd_arrays, gd_totals, threshold)
        ab = rates_from_sorted(ab_arrays, ab_totals, threshold)
        nk1 = rates_from_sorted(nk_arrays, nk_totals, threshold)
        eligible = (not gd or min(gd.values()) >= spec["minimum_gdt_recall_per_source"]) and (not ab or min(ab.values()) >= spec["minimum_abt_recall_per_source"]) and (not nk1 or max(nk1.values()) <= spec["maximum_tier1_nk_passage_per_source"])
        if eligible:
            best = (float(threshold), {"gdt_recall": gd, "abt_recall": ab, "tier1_nk_passage": nk1})
            break
    if best is None:
        # Record the maximum-recall fallback as ineligible; downstream profiles will fail closed.
        threshold = float(score.min())
        calls = score >= threshold
        best = (threshold, {"gdt_recall": source_rates(frame, calls, "gdT_gold"), "abt_recall": source_rates(frame, calls, "abT_gold"), "tier1_nk_passage": source_rates(frame, calls, "nk_tier1"), "eligible": False})
    return best


def profile_metrics(frame: pd.DataFrame, calls: np.ndarray, score: np.ndarray) -> dict:
    positive = frame.truth_class.eq("gdT_gold").to_numpy()
    negative = frame.truth_class.isin(["abT_gold", "nk_tier1", "nk_tier2"]).to_numpy()
    tp, fn = int((calls & positive).sum()), int((~calls & positive).sum())
    fp, tn = int((calls & negative).sum()), int((~calls & negative).sum())
    recall_by_source = source_rates(frame, calls, "gdT_gold")
    tier1 = source_rates(frame, calls, "nk_tier1")
    tier2 = source_rates(frame, calls, "nk_tier2")
    ab_mask = frame.truth_class.eq("abT_gold").to_numpy()
    y = positive[positive | negative].astype(int)
    s = score[positive | negative]
    auc_valid = np.unique(y).size == 2
    return {
        "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "precision": tp / (tp + fp) if tp + fp else 0.0,
        "recall": tp / (tp + fn) if tp + fn else 0.0,
        "specificity": tn / (tn + fp) if tn + fp else 0.0,
        "f1": fbeta(tp, fp, fn, 1.0), "f0_5": fbeta(tp, fp, fn, 0.5),
        "mcc": float(matthews_corrcoef(positive[positive | negative], calls[positive | negative])),
        "roc_auc": float(roc_auc_score(y, s)) if auc_valid else math.nan,
        "pr_auc": float(average_precision_score(y, s)) if positive.any() else math.nan,
        "abt_fpr": float(calls[ab_mask].mean()) if ab_mask.any() else math.nan,
        "recall_by_source": recall_by_source, "tier1_nk_fpr_by_source": tier1, "tier2_nk_fpr_by_source": tier2,
        "macro_recall": float(np.mean(list(recall_by_source.values()))) if recall_by_source else math.nan,
        "tier2_macro_nk_fpr": float(np.mean(list(tier2.values()))) if tier2 else 0.0,
    }


def choose_profile_threshold(frame: pd.DataFrame, stage1_score: np.ndarray, stage1_threshold: float,
                             stage2_score: np.ndarray, excluded: np.ndarray, spec: dict) -> tuple[float, dict]:
    thresholds = np.unique(np.quantile(stage2_score, np.linspace(0, 1, 3001)))
    base = (stage1_score >= stage1_threshold) & ~excluded
    truth_arrays = {}
    truth_totals = {}
    for truth in ("gdT_gold", "abT_gold", "nk_tier1", "nk_tier2"):
        truth_arrays[truth], truth_totals[truth] = sorted_scores_by_source(frame, stage2_score, truth, base)
    positive_total = sum(truth_totals["gdT_gold"].values())
    negative_total = sum(sum(truth_totals[name].values()) for name in ("abT_gold", "nk_tier1", "nk_tier2"))
    best = None
    near_miss = None
    for threshold in thresholds:
        recalls = rates_from_sorted(truth_arrays["gdT_gold"], truth_totals["gdT_gold"], threshold)
        ab = rates_from_sorted(truth_arrays["abT_gold"], truth_totals["abT_gold"], threshold)
        tier1 = rates_from_sorted(truth_arrays["nk_tier1"], truth_totals["nk_tier1"], threshold)
        tier2 = rates_from_sorted(truth_arrays["nk_tier2"], truth_totals["nk_tier2"], threshold)
        tp = int(round(sum(recalls[s] * truth_totals["gdT_gold"][s] for s in recalls)))
        fp = int(round(sum(ab[s] * truth_totals["abT_gold"][s] for s in ab) + sum(tier1[s] * truth_totals["nk_tier1"][s] for s in tier1) + sum(tier2[s] * truth_totals["nk_tier2"][s] for s in tier2)))
        fn, tn = positive_total - tp, negative_total - fp
        metrics = {
            "f1": fbeta(tp, fp, fn, 1.0), "f0_5": fbeta(tp, fp, fn, 0.5),
            "macro_recall": float(np.mean(list(recalls.values()))) if recalls else math.nan,
            "abt_fpr": float(sum(ab[s] * truth_totals["abT_gold"][s] for s in ab) / max(1, sum(truth_totals["abT_gold"].values()))),
            "tier2_macro_nk_fpr": float(np.mean(list(tier2.values()))) if tier2 else 0.0,
            "recall_by_source": recalls, "tier1_nk_fpr_by_source": tier1, "tier2_nk_fpr_by_source": tier2,
        }
        eligible = (
            metrics["macro_recall"] >= spec["minimum_macro_recall"]
            and (not recalls or min(recalls.values()) >= spec["minimum_per_source_recall"])
            and metrics["abt_fpr"] <= spec["maximum_abt_fpr"]
            and (not tier1 or max(tier1.values()) <= spec["maximum_tier1_nk_fpr"])
            and metrics["tier2_macro_nk_fpr"] <= spec["maximum_tier2_macro_nk_fpr"]
        )
        objective = metrics["f1"] if spec["objective_beta"] == 1 else metrics["f0_5"]
        violations = {
            "macro_recall_deficit": max(0.0, spec["minimum_macro_recall"] - metrics["macro_recall"]),
            "minimum_source_recall_deficit": max(0.0, spec["minimum_per_source_recall"] - min(recalls.values(), default=1.0)),
            "abt_fpr_excess": max(0.0, metrics["abt_fpr"] - spec["maximum_abt_fpr"]),
            "tier1_nk_fpr_excess": max(0.0, max(tier1.values(), default=0.0) - spec["maximum_tier1_nk_fpr"]),
            "tier2_macro_nk_fpr_excess": max(0.0, metrics["tier2_macro_nk_fpr"] - spec["maximum_tier2_macro_nk_fpr"]),
        }
        normalized_violation = (
            violations["macro_recall_deficit"] / max(spec["minimum_macro_recall"], 1e-12)
            + violations["minimum_source_recall_deficit"] / max(spec["minimum_per_source_recall"], 1e-12)
            + violations["abt_fpr_excess"] / max(spec["maximum_abt_fpr"], 1e-12)
            + violations["tier1_nk_fpr_excess"] / max(spec["maximum_tier1_nk_fpr"], 1e-12)
            + violations["tier2_macro_nk_fpr_excess"] / max(spec["maximum_tier2_macro_nk_fpr"], 1e-12)
        )
        if near_miss is None or (normalized_violation, -objective) < (near_miss[0], -near_miss[1]):
            near_miss = (normalized_violation, objective, float(threshold), metrics, violations)
        if eligible and (best is None or objective > best[2]):
            best = (float(threshold), metrics, objective)
    if best is None:
        payload = {"eligible": False}
        if near_miss is not None:
            payload.update({
                "near_miss_threshold": near_miss[2],
                "near_miss_metrics": near_miss[3],
                "near_miss_violations": near_miss[4],
                "near_miss_normalized_violation": near_miss[0],
            })
        return math.nan, payload
    threshold = best[0]
    calls = base & (stage2_score >= threshold)
    metrics = profile_metrics(frame, calls, stage1_score * stage2_score)
    metrics["eligible"] = True
    return threshold, metrics


def grid(config: dict) -> list[dict]:
    values = config["grid"]
    return [dict(zip(values, combination)) for combination in itertools.product(*(values[key] for key in values))]


def audit_rows(labels: pd.DataFrame, outer_train: np.ndarray, config: dict, seed: int) -> np.ndarray:
    local = labels.iloc[outer_train]
    positive_nk = local.truth_class.isin(["gdT_gold", "nk_tier1", "nk_tier2"]).to_numpy()
    ab_local = np.flatnonzero(local.truth_class.isin(["abT_gold", "single_ab_support"]).to_numpy())
    ab_selected = sampled_rows(local.reset_index(drop=True), ab_local, config["maximum_audit_ab_cells"], seed)
    return outer_train[np.unique(np.concatenate([np.flatnonzero(positive_nk), ab_selected]))]


def candidate_id(params: dict) -> str:
    return "__".join(f"{key}-{params[key]}" for key in sorted(params))


def compose_stage2_features(x: np.ndarray, rows: np.ndarray, columns: np.ndarray,
                            stage1_score: np.ndarray, include_stage1_probability: bool) -> np.ndarray:
    values = np.asarray(x[rows][:, columns], dtype=np.float32)
    if include_stage1_probability:
        values = np.column_stack([values, stage1_score]).astype(np.float32)
    return values


def selected_profile_candidate(stage2_table: pd.DataFrame, profile: str, beta: float) -> tuple[pd.Series, dict] | None:
    choices = []
    for _, row in stage2_table.iterrows():
        payload = json.loads(row.profiles_json)[profile]
        metrics = payload["metrics"]
        if metrics.get("eligible") is True:
            objective = metrics["f1"] if beta == 1 else metrics["f0_5"]
            choices.append((objective, row, payload))
    if not choices:
        return None
    _, row, payload = max(choices, key=lambda value: value[0])
    return row, payload


def run_outer(labels: pd.DataFrame, x: np.ndarray, feature_names: list[str], stage1_columns: np.ndarray,
              stage2_columns: np.ndarray, include_stage1_probability: bool,
              heldout: str, config: dict, candidates: list[dict], smoke: bool, model_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    development = labels.allow_fit.to_numpy()
    outer_test = development & labels.source_gse_id.eq(heldout).to_numpy()
    outer_train = np.flatnonzero(development & ~labels.source_gse_id.eq(heldout).to_numpy())
    audit = audit_rows(labels, outer_train, config, stable_seed(heldout, "audit"))
    audit_frame = labels.iloc[audit].reset_index(drop=True)
    inner = audit_frame.inner_fold.to_numpy(np.int8)
    excluded_audit = exclusion_flags(np.asarray(x[audit]), feature_names)[2]
    stage1_oof_by_candidate: dict[str, np.ndarray] = {}
    stage1_rows = []
    folds = [0, 1, 2]
    for params in candidates:
        model_id = candidate_id(params)
        oof = np.full(len(audit), np.nan, dtype=np.float32)
        for fold in folds:
            fit_candidates = audit[(inner != fold) & audit_frame.truth_class.isin(["gdT_gold", "abT_gold", "single_ab_support", "nk_tier1", "nk_tier2"]).to_numpy()]
            fit_frame = labels.iloc[fit_candidates]
            t_local = np.flatnonzero(~fit_frame.truth_class.str.startswith("nk_").to_numpy())
            nk_local = np.flatnonzero(fit_frame.truth_class.str.startswith("nk_").to_numpy())
            t_selected = sampled_rows(fit_frame.reset_index(drop=True), t_local, config["maximum_stage1_t_cells"], stable_seed(heldout, model_id, fold, "t"))
            train_rows = fit_candidates[np.concatenate([t_selected, nk_local])]
            y = (~labels.iloc[train_rows].truth_class.str.startswith("nk_")).to_numpy(np.int8)
            w = source_balanced_weights(labels, train_rows, 1)
            validation = np.flatnonzero(inner == fold)
            model = fit_xgb(np.asarray(x[train_rows][:, stage1_columns]), y, w, params, config["fixed"], stable_seed(heldout, model_id, fold, 1), (np.asarray(x[audit[validation]][:, stage1_columns]), (~audit_frame.iloc[validation].truth_class.str.startswith("nk_")).to_numpy(np.int8)))
            oof[validation] = predict(model, np.asarray(x[audit[validation]][:, stage1_columns]))
        decision = choose_stage1_threshold(audit_frame, oof, config["stage1"])
        stage1_oof_by_candidate[model_id] = oof
        stage1_rows.append({
            "heldout_source": heldout,
            "candidate_id": model_id,
            "threshold": decision[0],
            "eligible": decision[1].get("eligible", True) is True,
            "details_json": json.dumps(decision[1], sort_keys=True),
            "complete_oof": bool(np.isfinite(oof).all()),
        })
    stage1_table = pd.DataFrame(stage1_rows)
    eligible_stage1 = stage1_table[stage1_table.complete_oof & stage1_table.eligible]
    if eligible_stage1.empty:
        raise RuntimeError(f"No Stage-1 candidate passed the frozen source-wise gates for {heldout}")
    selected_stage1 = eligible_stage1.sort_values("threshold", ascending=False).iloc[0]
    stage1_params = next(p for p in candidates if candidate_id(p) == selected_stage1.candidate_id)
    s1_oof = stage1_oof_by_candidate[selected_stage1.candidate_id]
    s1_threshold = float(selected_stage1.threshold)

    stage2_rows = []
    for params in candidates:
        model_id = candidate_id(params)
        oof = np.full(len(audit), np.nan, dtype=np.float32)
        for fold in folds:
            train_local = np.flatnonzero((inner != fold) & audit_frame.truth_class.isin(["gdT_gold", "abT_gold", "single_ab_support"]).to_numpy())
            train_frame = audit_frame.iloc[train_local]
            pos_local = np.flatnonzero(train_frame.truth_class.eq("gdT_gold").to_numpy())
            neg_local = np.flatnonzero(~train_frame.truth_class.eq("gdT_gold").to_numpy())
            neg_selected = sampled_rows(train_frame.reset_index(drop=True), neg_local, config["maximum_stage2_negative_cells"], stable_seed(heldout, model_id, fold, "neg"))
            chosen_local = np.concatenate([pos_local, neg_selected])
            train_rows = audit[train_local[chosen_local]]
            z = compose_stage2_features(
                x, train_rows, stage2_columns, s1_oof[train_local[chosen_local]],
                include_stage1_probability,
            )
            y = labels.iloc[train_rows].truth_class.eq("gdT_gold").to_numpy(np.int8)
            w = source_balanced_weights(labels, train_rows, 2)
            validation = np.flatnonzero(inner == fold)
            validation_t = validation[audit_frame.iloc[validation].truth_class.isin(["gdT_gold", "abT_gold", "single_ab_support"]).to_numpy()]
            zv_t = compose_stage2_features(
                x, audit[validation_t], stage2_columns, s1_oof[validation_t],
                include_stage1_probability,
            )
            model = fit_xgb(z, y, w, params, config["fixed"], stable_seed(heldout, model_id, fold, 2), (zv_t, audit_frame.iloc[validation_t].truth_class.eq("gdT_gold").to_numpy(np.int8)))
            zv = compose_stage2_features(
                x, audit[validation], stage2_columns, s1_oof[validation],
                include_stage1_probability,
            )
            oof[validation] = predict(model, zv)
        valid = np.isfinite(oof)
        profiles = {}
        for name, spec in config["profiles"].items():
            profile_score = oof[valid].copy()
            rescue = receptor_rescue_flags(
                np.asarray(x[audit[valid]]), feature_names,
                config.get("receptor_rescue", {}).get(name),
            )
            profile_score[rescue] = np.finfo(np.float32).max
            threshold, metrics = choose_profile_threshold(
                audit_frame[valid].reset_index(drop=True), s1_oof[valid], s1_threshold,
                profile_score, excluded_audit[valid], spec,
            )
            profiles[name] = {"threshold": threshold, "metrics": metrics}
        stage2_rows.append({"heldout_source": heldout, "candidate_id": model_id, "complete_oof": bool(valid.all()), "profiles_json": json.dumps(profiles, sort_keys=True)})
    stage2_table = pd.DataFrame(stage2_rows)
    # Fit the selected Stage-1 model on all outer-development rows.
    outer_frame = labels.iloc[outer_train]
    t_local = np.flatnonzero(~outer_frame.truth_class.str.startswith("nk_").to_numpy())
    nk_local = np.flatnonzero(outer_frame.truth_class.str.startswith("nk_").to_numpy())
    t_selected = sampled_rows(outer_frame.reset_index(drop=True), t_local, config["maximum_stage1_t_cells"], stable_seed(heldout, "final_s1_t"))
    final_stage1_rows = outer_train[np.concatenate([t_selected, nk_local])]
    final_stage1_y = (~labels.iloc[final_stage1_rows].truth_class.str.startswith("nk_")).to_numpy(np.int8)
    final_stage1 = fit_xgb(
        np.asarray(x[final_stage1_rows][:, stage1_columns]), final_stage1_y,
        source_balanced_weights(labels, final_stage1_rows, 1), stage1_params, config["fixed"],
        stable_seed(heldout, "final_s1"), None,
    )
    outer_test_rows = np.flatnonzero(outer_test)
    outer_test_frame = labels.iloc[outer_test_rows].reset_index(drop=True)
    outer_s1 = predict(final_stage1, np.asarray(x[outer_test_rows][:, stage1_columns]))
    outer_excluded = exclusion_flags(np.asarray(x[outer_test_rows]), feature_names)[2]
    outer_results = []
    trained_stage2: dict[str, xgb.XGBClassifier] = {}
    outer_model_dir = model_dir / f"outer_{heldout}"
    outer_model_dir.mkdir(parents=True, exist_ok=True)
    final_stage1.save_model(outer_model_dir / "stage1.ubj")
    profile_contracts = {}
    for profile, profile_spec in config["profiles"].items():
        selection = selected_profile_candidate(stage2_table, profile, profile_spec["objective_beta"])
        if selection is None:
            outer_results.append({"heldout_source": heldout, "profile": profile, "eligible_inner": False})
            continue
        selected_row, payload = selection
        selected_id = selected_row.candidate_id
        params = next(value for value in candidates if candidate_id(value) == selected_id)
        if selected_id not in trained_stage2:
            train_local = np.flatnonzero(audit_frame.truth_class.isin(["gdT_gold", "abT_gold", "single_ab_support"]).to_numpy())
            train_frame = audit_frame.iloc[train_local]
            pos_local = np.flatnonzero(train_frame.truth_class.eq("gdT_gold").to_numpy())
            neg_local = np.flatnonzero(~train_frame.truth_class.eq("gdT_gold").to_numpy())
            neg_selected = sampled_rows(train_frame.reset_index(drop=True), neg_local, config["maximum_stage2_negative_cells"], stable_seed(heldout, selected_id, "final_neg"))
            chosen = train_local[np.concatenate([pos_local, neg_selected])]
            train_rows = audit[chosen]
            z = compose_stage2_features(
                x, train_rows, stage2_columns, s1_oof[chosen],
                include_stage1_probability,
            )
            trained_stage2[selected_id] = fit_xgb(
                z, labels.iloc[train_rows].truth_class.eq("gdT_gold").to_numpy(np.int8),
                source_balanced_weights(labels, train_rows, 2), params, config["fixed"],
                stable_seed(heldout, selected_id, "final_s2"), None,
            )
            trained_stage2[selected_id].save_model(outer_model_dir / f"stage2_{selected_id}.ubj")
        outer_s2 = predict(
            trained_stage2[selected_id],
            compose_stage2_features(
                x, outer_test_rows, stage2_columns, outer_s1,
                include_stage1_probability,
            ),
        )
        threshold = float(payload["threshold"])
        outer_rescue = receptor_rescue_flags(
            np.asarray(x[outer_test_rows]), feature_names,
            config.get("receptor_rescue", {}).get(profile),
        )
        calls = (outer_s1 >= s1_threshold) & ((outer_s2 >= threshold) | outer_rescue) & ~outer_excluded
        evaluation_score = outer_s1 * outer_s2
        evaluation_score[outer_rescue & (outer_s1 >= s1_threshold) & ~outer_excluded] = 1.0
        metrics = profile_metrics(outer_test_frame, calls, evaluation_score)
        outer_results.append({
            "heldout_source": heldout, "profile": profile, "eligible_inner": True,
            "stage1_candidate_id": selected_stage1.candidate_id, "stage1_threshold": s1_threshold,
            "stage2_candidate_id": selected_id, "stage2_threshold": threshold,
            **{key: value for key, value in metrics.items() if not isinstance(value, dict)},
            "recall_by_source_json": json.dumps(metrics["recall_by_source"], sort_keys=True),
            "tier1_nk_fpr_by_source_json": json.dumps(metrics["tier1_nk_fpr_by_source"], sort_keys=True),
            "tier2_nk_fpr_by_source_json": json.dumps(metrics["tier2_nk_fpr_by_source"], sort_keys=True),
        })
        profile_contracts[profile] = {
            "stage2_candidate_id": selected_id,
            "stage2_threshold": threshold,
            "eligible_inner": True,
            "receptor_rescue": config.get("receptor_rescue", {}).get(profile, {}),
        }
    contract = {
        "heldout_source": heldout, "stage1_candidate_id": selected_stage1.candidate_id,
        "stage1_threshold": s1_threshold, "feature_names": feature_names,
        "stage1_feature_names": [feature_names[i] for i in stage1_columns],
        "stage2_feature_names": [feature_names[i] for i in stage2_columns],
        "stage2_includes_stage1_probability": include_stage1_probability,
        "profiles": profile_contracts,
        "diagnostic_only": True,
        "lockbox_scored": False,
        "model_promoted": False,
    }
    (outer_model_dir / "contract.json").write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n")
    return stage1_table, stage2_table, pd.DataFrame(outer_results)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--feature-contract-only", action="store_true")
    parser.add_argument("--heldout-source", action="append", default=[])
    args = parser.parse_args()
    config = json.loads(resolve(args.config).read_text())
    table_dir, model_dir, log_dir = map(resolve, [config["output_table_dir"], config["output_model_dir"], config["output_log_dir"]])
    for path in (table_dir, model_dir, log_dir): path.mkdir(parents=True, exist_ok=True)
    cache_manifest = json.loads(resolve(config["cache_manifest"]).read_text())
    labels_path, features_path, matrix_path = map(resolve, [config["label_manifest"], config["feature_manifest"], config["cache_matrix"]])
    if cache_manifest["label_manifest_sha256"] != sha256_file(labels_path) or cache_manifest["feature_manifest_sha256"] != sha256_file(features_path):
        raise RuntimeError("Training inputs disagree with frozen cache contract")
    labels = pd.read_parquet(labels_path).reset_index(drop=True)
    split = pd.read_csv(resolve(config["split_manifest"]))
    group_to_fold = split.drop_duplicates("group_id").set_index("group_id")["inner_fold"]
    labels["inner_fold"] = labels.group_id.map(group_to_fold)
    if labels.loc[labels.allow_fit, "inner_fold"].isna().any():
        raise RuntimeError("Frozen inner-fold assignment is incomplete")
    features = pd.read_csv(features_path).sort_values("feature_index")
    names = features.gene.astype(str).tolist()
    stage1_excluded = set(config.get("stage1_exclude_genes", []))
    stage1_columns = np.flatnonzero(
        features.stage1.to_numpy(bool)
        & ~features.gene.astype(str).isin(stage1_excluded).to_numpy()
    )
    stage2_policy = config.get("stage2_feature_policy", {})
    if stage2_policy:
        controls = set(stage2_policy.get("control_genes", []))
        stage2_columns = np.flatnonzero(
            features.feature_class.eq("TCR").to_numpy()
            | features.gene.astype(str).isin(controls).to_numpy()
        )
    else:
        stage2_columns = np.flatnonzero(features.stage2.to_numpy(bool))
    include_stage1_probability = bool(config.get("stage2_include_stage1_probability", True))
    if len(stage2_columns) == 0:
        raise RuntimeError("Stage-2 feature policy selected zero genes")
    effective_features = features[["feature_index", "gene", "feature_class"]].copy()
    effective_features["used_by_stage1_t_nk_gate"] = False
    effective_features["used_by_stage2_gdt_classifier"] = False
    effective_features.loc[stage1_columns, "used_by_stage1_t_nk_gate"] = True
    effective_features.loc[stage2_columns, "used_by_stage2_gdt_classifier"] = True
    effective_features["stage2_role"] = np.where(
        effective_features.used_by_stage2_gdt_classifier,
        np.where(effective_features.feature_class.eq("TCR"), "individual_TCR_gene", "lineage_control"),
        "excluded",
    )
    effective_features.to_csv(table_dir / "effective_feature_contract.csv", index=False)
    feature_contract = {
        "protocol_version": config.get("protocol_version", "unknown"),
        "n_cached_genes": int(len(effective_features)),
        "n_stage1_genes": int(len(stage1_columns)),
        "n_stage2_genes": int(len(stage2_columns)),
        "stage1_excluded_genes": sorted(stage1_excluded),
        "stage2_includes_stage1_probability": include_stage1_probability,
        "stage2_feature_policy": stage2_policy,
        "cytotoxic_or_nk_context_genes_are_stage2_evidence": False,
    }
    (log_dir / "effective_feature_contract.json").write_text(
        json.dumps(feature_contract, indent=2, sort_keys=True) + "\n"
    )
    if args.feature_contract_only:
        print(json.dumps(feature_contract, indent=2, sort_keys=True))
        return
    x = np.load(matrix_path, mmap_mode="r")
    candidates = grid(config)
    if args.smoke: candidates = candidates[:1]
    heldouts = args.heldout_source or POSITIVE_SOURCES
    started = time.monotonic()
    all_stage1, all_stage2, all_outer = [], [], []
    for heldout in heldouts:
        s1, s2, outer = run_outer(
            labels, x, names, stage1_columns, stage2_columns,
            include_stage1_probability, heldout, config, candidates,
            args.smoke, model_dir,
        )
        all_stage1.append(s1); all_stage2.append(s2); all_outer.append(outer)
        pd.concat(all_stage1).to_csv(table_dir / ("smoke_stage1_candidates.csv" if args.smoke else "nested_stage1_candidates.csv"), index=False)
        pd.concat(all_stage2).to_csv(table_dir / ("smoke_stage2_candidates.csv" if args.smoke else "nested_stage2_candidates.csv"), index=False)
        pd.concat(all_outer).to_csv(table_dir / ("smoke_outer_metrics.csv" if args.smoke else "nested_outer_metrics.csv"), index=False)
    result = {"status": "PASS_SMOKE" if args.smoke else "PASS_NESTED_DEVELOPMENT", "heldout_sources": heldouts, "n_candidates": len(candidates), "runtime_seconds": time.monotonic() - started, "lockbox_scored": False, "model_promoted": False}
    (log_dir / ("smoke_summary.json" if args.smoke else "nested_summary.json")).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
