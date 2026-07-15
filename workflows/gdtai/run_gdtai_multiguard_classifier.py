#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Train and audit a gdTAI multi-guard candidate model.

This lane extends the NK-guard attempt by adding alpha-beta TCR guard negatives
into both training and threshold selection. The H5AD is opened read-only. The
current shared gdTAI model is not overwritten unless every acceptance gate passes.
"""

from __future__ import annotations

# TNK_WORKFLOW_BOOTSTRAP
import sys as _tnk_sys
from pathlib import Path as _TnkPath

_TNK_PROJECT_ROOT = _TnkPath(__file__).resolve().parents[2]
for _tnk_path in (
    _TNK_PROJECT_ROOT,
    _TNK_PROJECT_ROOT / "src",
    _TNK_PROJECT_ROOT / "workflows" / "integration",
    _TNK_PROJECT_ROOT / "workflows" / "intake",
    _TNK_PROJECT_ROOT / "workflows" / "metadata",
    _TNK_PROJECT_ROOT / "workflows" / "analysis",
    _TNK_PROJECT_ROOT / "workflows" / "gdtai",
    _TNK_PROJECT_ROOT / "workflows" / "gdt_atlas",
    _TNK_PROJECT_ROOT / "workflows" / "reporting",
    _TNK_PROJECT_ROOT / "workflows" / "maintenance",
):
    _tnk_value = str(_tnk_path)
    if _tnk_value not in _tnk_sys.path:
        _tnk_sys.path.insert(0, _tnk_value)


import argparse
import html
import json
import logging
import pickle
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

from run_gdt_prediction_package_evaluation import (
    DEFAULT_INPUT_H5AD,
    clean_group_values,
    dataframe_to_html,
    dataframe_to_markdown,
    normalize_strings,
    read_bool_obs,
    read_float_obs,
    read_obs_column,
    read_string_dataset,
)
from run_gdt_deg_tcr_classifier_training import metric_dict, safe_div
from run_gdt_gse144469_holdout_tcrgene_classifier import (
    EXTRA_TRAB_HOLDOUT_SOURCE,
    GDT2020_HOLDOUT_TISSUE,
    GDT2020_SOURCE,
    HOLDOUT_SOURCE,
    MAX_FEATURE_GENES,
    MAX_NEGATIVE_TRAIN,
    PENALTY_CONTROL_GENES,
    RANDOM_SEED,
    SUBOPTIMAL_SOURCE,
    TCR_PREFIXES,
    FeatureSpec,
    bounded_scatter_limits,
    build_obs_metadata,
    extract_gene_matrix,
    load_original_trd_trab_threshold,
    read_nonempty_if_present,
    summarize_full_predictions,
    tcr_family,
    tcr_priority,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
STATIC_DIR = PROJECT_ROOT / "gdT_prediction"
OUT_PREFIX = "gdtai_multiguard"
TABLE_DIR = OUTPUT_ROOT / "tables" / "gdT_prediction" / OUT_PREFIX
FIGURE_DIR = OUTPUT_ROOT / "figures" / "gdT_prediction" / OUT_PREFIX
LOG_DIR = OUTPUT_ROOT / "logs" / "gdT_prediction" / OUT_PREFIX
MODEL_DIR = OUTPUT_ROOT / "models" / "gdT_prediction_classifier" / OUT_PREFIX
REPORT_MD = LOG_DIR / "gdtai_multiguard_report.md"
REPORT_HTML = STATIC_DIR / "gdtai_multiguard_report.html"
REPORT_PDF = STATIC_DIR / "gdtai_multiguard_report.pdf"
RUN_LOG = LOG_DIR / "run.log"
CURRENT_MODEL_PATH = OUTPUT_ROOT / "models" / "gdT_prediction_classifier" / "gse144469_holdout_tcrgene" / "selected_model.pkl"
CURRENT_TABLE_DIR = OUTPUT_ROOT / "tables" / "gdT_prediction" / "gse144469_holdout_tcrgene"

NK_STRESS_SOURCES = ("GSE243572", "GSE243905")
NK_CONTROL_GENES = [
    "NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "KLRD1", "KLRF1", "FCGR3A", "NCAM1", "TYROBP", "FCER1G", "CST7",
]
MAX_NK_GUARD_TRAIN = 150_000
MAX_TCRAB_GUARD_TRAIN = 180_000
NK_TUNE_FPR_MAX = 0.0175
TCRAB_TUNE_FPR_MAX = 0.0020
PRIMARY_TUNE_F1_FRACTION = 0.95
FULL_APPLY_CHUNK = 50_000
SCATTER_SAMPLE_CELLS = 250_000

ACCEPT_PRIMARY_F1_MIN = 0.887
ACCEPT_PRIMARY_RECALL_MIN = 0.84
ACCEPT_PRIMARY_PRECISION_MIN = 0.94
ACCEPT_NK_STRESS_REDUCTION_MIN = 0.30
ACCEPT_FULL_NK_FP_MAX = 4_500
ACCEPT_FULL_PAIRED_TCRAB_FP_MAX = 2_600
ACCEPT_FULL_KNOWN_FP_FRACTION_MAX = 0.025


@dataclass
class SplitBundle:
    train_idx: np.ndarray
    tune_idx: np.ndarray
    validation_idx: np.ndarray
    nk_train_idx: np.ndarray
    nk_tune_idx: np.ndarray
    nk_stress_idx: np.ndarray
    tcrab_train_idx: np.ndarray
    tcrab_tune_idx: np.ndarray
    tcrab_stress_idx: np.ndarray
    excluded_idx: np.ndarray
    split_overall: pd.DataFrame
    split_by_source: pd.DataFrame
    annotation: np.ndarray
    nk_guard_mask: np.ndarray
    tcrab_guard_mask: np.ndarray


@dataclass
class Candidate:
    model: str
    feature_set: str
    spec: FeatureSpec
    model_object: Any
    threshold: float
    threshold_valid: bool
    threshold_reason: str
    primary_tune_metrics: dict[str, Any]
    nk_tune_fpr: float
    nk_tune_predicted: int
    tcrab_tune_fpr: float
    tcrab_tune_predicted: int
    current_primary_tune_f1: float
    primary_validation_metrics: dict[str, Any]
    nk_stress_fpr: float
    nk_stress_predicted: int
    current_nk_stress_fpr: float
    current_nk_stress_predicted: int
    tcrab_stress_fpr: float
    tcrab_stress_predicted: int
    current_tcrab_stress_fpr: float
    current_tcrab_stress_predicted: int
    notes: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train gdTAI multi-guard candidate models.")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--current-model-pkl", type=Path, default=CURRENT_MODEL_PATH)
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    parser.add_argument("--max-negative-train", type=int, default=MAX_NEGATIVE_TRAIN)
    parser.add_argument("--max-nk-guard-train", type=int, default=MAX_NK_GUARD_TRAIN)
    parser.add_argument("--max-tcrab-guard-train", type=int, default=MAX_TCRAB_GUARD_TRAIN)
    parser.add_argument("--scatter-sample-cells", type=int, default=SCATTER_SAMPLE_CELLS)
    parser.add_argument("--threshold-only", action="store_true", help="Only retune the current gdTAI model threshold against NK and TCRAB guard sets.")
    parser.add_argument("--skip-xgb", action="store_true", help="Skip the XGBoost candidate for a faster, lower-overfitting logistic-only run.")
    parser.add_argument("--no-pdf", action="store_true")
    return parser.parse_args()


def ensure_dirs() -> None:
    for path in [TABLE_DIR, FIGURE_DIR, LOG_DIR, MODEL_DIR, STATIC_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def setup_logging() -> None:
    ensure_dirs()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(RUN_LOG, mode="w", encoding="utf-8"), logging.StreamHandler()],
        force=True,
    )


def multiguard_family(gene: str) -> str:
    if gene in NK_CONTROL_GENES:
        return "NK_control"
    return tcr_family(gene)


def build_feature_spec(handle: h5py.File, *, include_nk_controls: bool, name: str) -> FeatureSpec:
    var_names = pd.Index(read_string_dataset(handle["var"]["_index"]), dtype="string").astype(str).tolist()
    gene_to_idx = {gene: idx for idx, gene in enumerate(var_names)}
    tcr_genes = sorted([gene for gene in var_names if gene.upper().startswith(TCR_PREFIXES)], key=tcr_priority)
    controls = [gene for gene in PENALTY_CONTROL_GENES if gene in gene_to_idx]
    nk_controls = [gene for gene in NK_CONTROL_GENES if gene in gene_to_idx] if include_nk_controls else []
    selected: list[str] = []
    for gene in [*tcr_genes, *controls, *nk_controls]:
        if gene not in selected:
            selected.append(gene)
    if len(selected) > MAX_FEATURE_GENES:
        selected = selected[:MAX_FEATURE_GENES]
    name_to_col = {gene: idx for idx, gene in enumerate(selected)}
    spec = FeatureSpec(
        feature_names=[f"{gene}_log1p_cp10k" for gene in selected],
        gene_names=selected,
        gene_indices=np.asarray([gene_to_idx[gene] for gene in selected], dtype=np.int32),
        tcr_gene_names=[gene for gene in selected if gene.upper().startswith(TCR_PREFIXES)],
        penalty_gene_names=[gene for gene in selected if gene in controls or gene in nk_controls],
        family_by_feature=[multiguard_family(gene) for gene in selected],
        cd3_cols=[name_to_col[gene] for gene in ["CD3D", "CD3E", "CD3G"] if gene in name_to_col],
        cd4_col=name_to_col.get("CD4"),
        foxp3_col=name_to_col.get("FOXP3"),
    )
    pd.DataFrame(
        {
            "feature": spec.feature_names,
            "gene": spec.gene_names,
            "family": spec.family_by_feature,
            "feature_index": np.arange(len(spec.feature_names), dtype=int),
        }
    ).to_csv(TABLE_DIR / f"feature_manifest_{name}.csv", index=False)
    return spec


def split_by_source(rows_idx: np.ndarray, source: np.ndarray, label: str, rng: np.random.Generator, rows: list[dict[str, Any]]) -> tuple[np.ndarray, np.ndarray]:
    train_parts: list[np.ndarray] = []
    tune_parts: list[np.ndarray] = []
    if rows_idx.size == 0:
        return np.asarray([], dtype=np.int64), np.asarray([], dtype=np.int64)
    df = pd.DataFrame({"obs_index": rows_idx, "source_gse_id": source[rows_idx]})
    for src, group in df.groupby("source_gse_id", sort=True):
        idx = group["obs_index"].to_numpy(dtype=np.int64)
        rng.shuffle(idx)
        n_tune = 0 if idx.size < 5 else max(1, int(round(idx.size * 0.20)))
        if n_tune:
            tune_parts.append(idx[:n_tune])
            train_parts.append(idx[n_tune:])
        else:
            train_parts.append(idx)
        rows.append({"source_gse_id": src, "label": label, "n_cells": int(idx.size), "train": int(idx.size - n_tune), "tune": int(n_tune)})
    train = np.concatenate(train_parts).astype(np.int64) if train_parts else np.asarray([], dtype=np.int64)
    tune = np.concatenate(tune_parts).astype(np.int64) if tune_parts else np.asarray([], dtype=np.int64)
    return train, tune


def build_splits(handle: h5py.File, seed: int) -> tuple[Any, SplitBundle]:
    obs = build_obs_metadata(handle)
    n_obs = obs.source.size
    annotation = clean_group_values(read_obs_column(handle, "simple_annotation_plus6")) if "simple_annotation_plus6" in handle["obs"] else np.full(n_obs, "unknown", dtype=object)
    annotation_s = pd.Series(annotation, copy=False).astype(str)
    is_nk_annotation = annotation_s.str.upper().eq("NK_CELL").to_numpy(dtype=bool)
    primary = (obs.class_code == 1) | (obs.class_code == 2)

    validation_primary_idx = np.flatnonzero(primary & (obs.source == HOLDOUT_SOURCE)).astype(np.int64)
    validation_gdt2020_idx = np.flatnonzero(
        primary
        & (obs.source == GDT2020_SOURCE)
        & (obs.class_code == 2)
        & (pd.Series(obs.tissue, copy=False).astype(str).str.lower().to_numpy() == GDT2020_HOLDOUT_TISSUE.lower())
    ).astype(np.int64)
    validation_extra_trab_idx = np.flatnonzero(
        (obs.source == EXTRA_TRAB_HOLDOUT_SOURCE)
        & obs.has_TRA_TRB_paired
        & (~obs.corrected_has_any_gd_tcr)
    ).astype(np.int64)
    validation_idx = np.unique(np.concatenate([validation_primary_idx, validation_gdt2020_idx, validation_extra_trab_idx])).astype(np.int64)
    validation_mask = np.zeros(n_obs, dtype=bool)
    validation_mask[validation_idx] = True
    excluded_idx = np.flatnonzero(primary & (obs.source == SUBOPTIMAL_SOURCE)).astype(np.int64)

    train_tune_pool = np.flatnonzero(
        primary
        & (~validation_mask)
        & (obs.source != HOLDOUT_SOURCE)
        & (obs.source != SUBOPTIMAL_SOURCE)
        & (obs.source != EXTRA_TRAB_HOLDOUT_SOURCE)
    ).astype(np.int64)
    if validation_primary_idx.size == 0 or validation_gdt2020_idx.size == 0 or validation_extra_trab_idx.size == 0:
        raise RuntimeError("Primary validation cohorts are incomplete.")
    validation_y = (obs.class_code[validation_idx] == 2).astype(np.int8)
    if validation_y.sum() == 0 or (validation_y == 0).sum() == 0:
        raise RuntimeError("Combined primary validation must contain positives and negatives.")

    tra_nonempty = read_nonempty_if_present(handle, "TRA_cdr3", n_obs)
    trb_nonempty = read_nonempty_if_present(handle, "TRB_cdr3", n_obs)
    ab_evidence = obs.has_TRA_TRB_paired | obs.has_any_ab_tcr | tra_nonempty | trb_nonempty

    nk_guard_mask = (
        is_nk_annotation
        & (obs.class_code != 2)
        & (~obs.silver_mask)
        & (~obs.corrected_has_any_gd_tcr)
    )
    tcrab_guard_mask = (
        ab_evidence
        & (obs.class_code != 2)
        & (~obs.silver_mask)
        & (~obs.corrected_has_any_gd_tcr)
    )

    nk_stress_idx = np.flatnonzero(nk_guard_mask & np.isin(obs.source, NK_STRESS_SOURCES)).astype(np.int64)
    nk_train_tune_pool = np.flatnonzero(
        nk_guard_mask
        & (~np.isin(obs.source, NK_STRESS_SOURCES))
        & (obs.source != HOLDOUT_SOURCE)
        & (obs.source != EXTRA_TRAB_HOLDOUT_SOURCE)
        & (obs.source != GDT2020_SOURCE)
        & (obs.source != SUBOPTIMAL_SOURCE)
        & (~validation_mask)
    ).astype(np.int64)
    tcrab_stress_idx = np.flatnonzero(tcrab_guard_mask & (obs.source == EXTRA_TRAB_HOLDOUT_SOURCE)).astype(np.int64)
    tcrab_train_tune_pool = np.flatnonzero(
        tcrab_guard_mask
        & (obs.source != HOLDOUT_SOURCE)
        & (obs.source != EXTRA_TRAB_HOLDOUT_SOURCE)
        & (obs.source != GDT2020_SOURCE)
        & (obs.source != SUBOPTIMAL_SOURCE)
        & (~np.isin(obs.source, NK_STRESS_SOURCES))
        & (~validation_mask)
    ).astype(np.int64)
    if nk_stress_idx.size == 0 or nk_train_tune_pool.size == 0:
        raise RuntimeError("NK-guard train/tune or stress validation set is empty.")
    if tcrab_stress_idx.size == 0 or tcrab_train_tune_pool.size == 0:
        raise RuntimeError("TCRAB-guard train/tune or stress validation set is empty.")

    rng = np.random.default_rng(seed)
    train_parts: list[np.ndarray] = []
    tune_parts: list[np.ndarray] = []
    rows: list[dict[str, Any]] = []
    split_df = pd.DataFrame(
        {
            "obs_index": train_tune_pool,
            "source_gse_id": obs.source[train_tune_pool],
            "label": np.where(obs.class_code[train_tune_pool] == 2, "gdT_gold", "abT_gold"),
        }
    )
    for (source, label), group in split_df.groupby(["source_gse_id", "label"], sort=True):
        idx = group["obs_index"].to_numpy(dtype=np.int64)
        rng.shuffle(idx)
        n_tune = 0 if idx.size < 5 else max(1, int(round(idx.size * 0.20)))
        if n_tune:
            tune_parts.append(idx[:n_tune])
            train_parts.append(idx[n_tune:])
        else:
            train_parts.append(idx)
        rows.append({"source_gse_id": source, "label": label, "n_cells": int(idx.size), "train": int(idx.size - n_tune), "tune": int(n_tune)})
    train_idx = np.concatenate(train_parts).astype(np.int64)
    tune_idx = np.concatenate(tune_parts).astype(np.int64)

    nk_train_idx, nk_tune_idx = split_by_source(nk_train_tune_pool, obs.source, "NK_guard_negative", rng, rows)
    tcrab_train_idx, tcrab_tune_idx = split_by_source(tcrab_train_tune_pool, obs.source, "TCRAB_guard_negative", rng, rows)

    split_by_source_df = pd.DataFrame(rows)
    split_by_source_df.to_csv(TABLE_DIR / "split_by_source_label.csv", index=False)
    overall_rows = []
    for split, idx in [
        ("train_primary", train_idx),
        ("tune_primary", tune_idx),
        ("train_NK_guard", nk_train_idx),
        ("tune_NK_guard", nk_tune_idx),
        ("train_TCRAB_guard_single_or_paired", tcrab_train_idx),
        ("tune_TCRAB_guard_single_or_paired", tcrab_tune_idx),
        ("validation_primary_combined", validation_idx),
        (f"validation_NK_stress_{'_'.join(NK_STRESS_SOURCES)}", nk_stress_idx),
        (f"validation_TCRAB_stress_{EXTRA_TRAB_HOLDOUT_SOURCE}", tcrab_stress_idx),
        (f"sensitivity_excluded_{SUBOPTIMAL_SOURCE}", excluded_idx),
    ]:
        labels = obs.class_code[idx]
        overall_rows.append(
            {
                "split": split,
                "n_cells": int(idx.size),
                "gdT_gold": int((labels == 2).sum()),
                "abT_gold": int((labels == 1).sum()),
                "NK_guard_negative": int(nk_guard_mask[idx].sum()),
                "TCRAB_guard_negative": int(tcrab_guard_mask[idx].sum()),
                "gdT_prevalence": safe_div(int((labels == 2).sum()), int(idx.size)),
            }
        )
    split_overall = pd.DataFrame(overall_rows)
    split_overall.to_csv(TABLE_DIR / "split_overall.csv", index=False)
    conflicts = {
        "nk_guard_conflict_with_gdT_evidence": int((nk_guard_mask & ((obs.class_code == 2) | obs.silver_mask | obs.corrected_has_any_gd_tcr)).sum()),
        "tcrab_guard_conflict_with_gdT_evidence": int((tcrab_guard_mask & ((obs.class_code == 2) | obs.silver_mask | obs.corrected_has_any_gd_tcr)).sum()),
        "nk_and_tcrab_guard_overlap": int((nk_guard_mask & tcrab_guard_mask).sum()),
    }
    pd.DataFrame([conflicts]).to_csv(TABLE_DIR / "guard_label_conflict_audit.csv", index=False)
    if conflicts["nk_guard_conflict_with_gdT_evidence"] or conflicts["tcrab_guard_conflict_with_gdT_evidence"]:
        raise RuntimeError(f"Guard label conflict detected: {conflicts}")
    return obs, SplitBundle(
        train_idx=train_idx,
        tune_idx=tune_idx,
        validation_idx=validation_idx,
        nk_train_idx=nk_train_idx,
        nk_tune_idx=nk_tune_idx,
        nk_stress_idx=nk_stress_idx,
        tcrab_train_idx=tcrab_train_idx,
        tcrab_tune_idx=tcrab_tune_idx,
        tcrab_stress_idx=tcrab_stress_idx,
        excluded_idx=excluded_idx,
        split_overall=split_overall,
        split_by_source=split_by_source_df,
        annotation=annotation,
        nk_guard_mask=nk_guard_mask,
        tcrab_guard_mask=tcrab_guard_mask,
    )


def local_positions(eval_rows: np.ndarray, target_rows: np.ndarray) -> np.ndarray:
    lookup = pd.Series(np.arange(eval_rows.size, dtype=np.int64), index=eval_rows)
    return lookup.loc[target_rows].to_numpy(dtype=np.int64)


def sample_train_positions(
    pos_train_primary: np.ndarray,
    pos_train_nk: np.ndarray,
    pos_train_tcrab: np.ndarray,
    y_eval: np.ndarray,
    is_nk_eval: np.ndarray,
    is_tcrab_eval: np.ndarray,
    *,
    max_negative: int,
    max_nk: int,
    max_tcrab: int,
    seed: int,
) -> tuple[np.ndarray, pd.DataFrame, np.ndarray]:
    rng = np.random.default_rng(seed)
    positives = pos_train_primary[y_eval[pos_train_primary] == 1]
    primary_neg = pos_train_primary[y_eval[pos_train_primary] == 0]
    nk_neg = pos_train_nk[is_nk_eval[pos_train_nk]]
    tcrab_neg = pos_train_tcrab[is_tcrab_eval[pos_train_tcrab]]
    if primary_neg.size > max_negative:
        primary_neg = rng.choice(primary_neg, size=max_negative, replace=False)
    if nk_neg.size > max_nk:
        nk_neg = rng.choice(nk_neg, size=max_nk, replace=False)
    if tcrab_neg.size > max_tcrab:
        tcrab_neg = rng.choice(tcrab_neg, size=max_tcrab, replace=False)
    raw = np.concatenate([positives, primary_neg, nk_neg, tcrab_neg]).astype(np.int64)
    out = np.unique(raw).astype(np.int64)
    rng.shuffle(out)
    rows = [
        {"sample_class": "gdT_positive", "n_cells": int(y_eval[out].sum())},
        {"sample_class": "primary_abT_negative", "n_cells": int(((y_eval[out] == 0) & (~is_nk_eval[out]) & (~is_tcrab_eval[out])).sum())},
        {"sample_class": "NK_guard_negative", "n_cells": int(is_nk_eval[out].sum())},
        {"sample_class": "TCRAB_guard_negative_single_or_paired", "n_cells": int(is_tcrab_eval[out].sum())},
        {"sample_class": "NK_and_TCRAB_guard_overlap", "n_cells": int((is_nk_eval[out] & is_tcrab_eval[out]).sum())},
        {"sample_class": "unique_training_cells", "n_cells": int(out.size)},
    ]
    weights = np.ones(out.size, dtype=np.float32)
    weights[is_nk_eval[out]] *= 1.5
    weights[is_tcrab_eval[out]] *= 2.0
    return out, pd.DataFrame(rows), weights


def score_model(model: Any, x: np.ndarray) -> np.ndarray:
    return model.predict_proba(x)[:, 1].astype(np.float32)


def choose_multiguard_threshold(
    primary_y: np.ndarray,
    primary_score: np.ndarray,
    nk_score: np.ndarray,
    tcrab_score: np.ndarray,
    *,
    current_primary_tune_f1: float,
) -> tuple[float, dict[str, Any], float, int, float, int, bool, str]:
    fpr, tpr, thresholds = roc_curve(primary_y, primary_score)
    keep = np.isfinite(thresholds)
    thresholds = thresholds[keep]
    fpr = fpr[keep]
    tpr = tpr[keep]
    if thresholds.size == 0:
        threshold = float(np.nanmax(primary_score) + 1e-6)
        metrics = metric_dict(primary_y, (primary_score >= threshold).astype(np.int8), primary_score)
        return threshold, metrics, float("nan"), 0, float("nan"), 0, False, "no_finite_thresholds"

    n_pos = int(primary_y.sum())
    n_neg = int((primary_y == 0).sum())
    tp = np.rint(tpr * n_pos).astype(np.int64)
    fp = np.rint(fpr * n_neg).astype(np.int64)
    fn = n_pos - tp
    tn = n_neg - fp
    pred_pos = tp + fp
    precision = np.divide(tp, pred_pos, out=np.zeros_like(tp, dtype=float), where=pred_pos > 0)
    recall = np.divide(tp, tp + fn, out=np.zeros_like(tp, dtype=float), where=(tp + fn) > 0)
    specificity = np.divide(tn, tn + fp, out=np.zeros_like(tn, dtype=float), where=(tn + fp) > 0)
    f1 = np.divide(2.0 * precision * recall, precision + recall, out=np.zeros_like(precision), where=(precision + recall) > 0)

    nk_sorted = np.sort(nk_score)
    tcrab_sorted = np.sort(tcrab_score)
    nk_predicted = nk_score.size - np.searchsorted(nk_sorted, thresholds, side="left")
    tcrab_predicted = tcrab_score.size - np.searchsorted(tcrab_sorted, thresholds, side="left")
    nk_fpr = nk_predicted / max(nk_score.size, 1)
    tcrab_fpr = tcrab_predicted / max(tcrab_score.size, 1)

    min_f1 = PRIMARY_TUNE_F1_FRACTION * current_primary_tune_f1
    nonempty = pred_pos > 0
    valid = nonempty & (f1 >= min_f1) & (nk_fpr <= NK_TUNE_FPR_MAX) & (tcrab_fpr <= TCRAB_TUNE_FPR_MAX)
    if valid.any():
        pool = np.flatnonzero(valid)
        best = max(pool, key=lambda i: (f1[i], -tcrab_fpr[i], -nk_fpr[i], precision[i], recall[i], specificity[i], thresholds[i]))
        reason = "valid_multi_guard_threshold"
        is_valid = True
    else:
        pool = np.flatnonzero(nonempty)
        if pool.size == 0:
            threshold = float(np.nanmax(primary_score) + 1e-6)
            metrics = metric_dict(primary_y, (primary_score >= threshold).astype(np.int8), primary_score)
            return threshold, metrics, float("nan"), 0, float("nan"), 0, False, "no_nonempty_primary_predictions"
        violation = np.maximum(0.0, nk_fpr - NK_TUNE_FPR_MAX) / NK_TUNE_FPR_MAX + np.maximum(0.0, tcrab_fpr - TCRAB_TUNE_FPR_MAX) / TCRAB_TUNE_FPR_MAX
        best = max(pool, key=lambda i: (-violation[i], f1[i], -tcrab_fpr[i], -nk_fpr[i], precision[i], recall[i], thresholds[i]))
        reason = "fallback_no_threshold_met_primary_NK_and_TCRAB_constraints"
        is_valid = False

    threshold = float(thresholds[best])
    metrics = metric_dict(primary_y, (primary_score >= threshold).astype(np.int8), primary_score)
    return threshold, metrics, float(nk_fpr[best]), int(nk_predicted[best]), float(tcrab_fpr[best]), int(tcrab_predicted[best]), is_valid, reason

def train_candidates(
    x_base: np.ndarray,
    x_nk: np.ndarray,
    y_eval: np.ndarray,
    is_nk_guard_eval: np.ndarray,
    is_tcrab_guard_eval: np.ndarray,
    positions: dict[str, np.ndarray],
    current_payload: dict[str, Any],
    base_spec: FeatureSpec,
    nk_spec: FeatureSpec,
    args: argparse.Namespace,
) -> list[Candidate]:
    primary_tune = positions["primary_tune"]
    nk_tune = positions["nk_tune"]
    tcrab_tune = positions["tcrab_tune"]
    primary_validation = positions["primary_validation"]
    nk_stress = positions["nk_stress"]
    tcrab_stress = positions["tcrab_stress"]
    primary_train = positions["primary_train"]
    nk_train = positions["nk_train"]
    tcrab_train = positions["tcrab_train"]
    current_model = current_payload["model_object"]
    current_threshold = float(current_payload["threshold"])
    current_tune_score = score_model(current_model, x_base[primary_tune])
    current_tune_metrics = metric_dict(y_eval[primary_tune], (current_tune_score >= current_threshold).astype(np.int8), current_tune_score)
    current_primary_tune_f1 = float(current_tune_metrics["f1"])
    current_nk_stress_score = score_model(current_model, x_base[nk_stress])
    current_tcrab_stress_score = score_model(current_model, x_base[tcrab_stress])
    current_nk_predicted = int((current_nk_stress_score >= current_threshold).sum())
    current_nk_fpr = safe_div(current_nk_predicted, current_nk_stress_score.size)
    current_tcrab_predicted = int((current_tcrab_stress_score >= current_threshold).sum())
    current_tcrab_fpr = safe_div(current_tcrab_predicted, current_tcrab_stress_score.size)

    train_sample, train_comp, train_weights = sample_train_positions(
        primary_train,
        nk_train,
        tcrab_train,
        y_eval,
        is_nk_guard_eval,
        is_tcrab_guard_eval,
        max_negative=args.max_negative_train,
        max_nk=args.max_nk_guard_train,
        max_tcrab=args.max_tcrab_guard_train,
        seed=args.seed,
    )
    train_comp.to_csv(TABLE_DIR / "training_sample_composition.csv", index=False)
    logging.info("Multi-guard training sample: %s unique cells", f"{train_sample.size:,}")

    candidates: list[Candidate] = []

    def add_candidate(name: str, feature_set: str, spec: FeatureSpec, model_object: Any, x: np.ndarray, notes: str) -> None:
        tune_score = score_model(model_object, x[primary_tune])
        nk_tune_score = score_model(model_object, x[nk_tune])
        tcrab_tune_score = score_model(model_object, x[tcrab_tune])
        threshold, tune_metrics, nk_tune_fpr, nk_tune_predicted, tcrab_tune_fpr, tcrab_tune_predicted, valid, reason = choose_multiguard_threshold(
            y_eval[primary_tune], tune_score, nk_tune_score, tcrab_tune_score, current_primary_tune_f1=current_primary_tune_f1
        )
        validation_score = score_model(model_object, x[primary_validation])
        validation_pred = (validation_score >= threshold).astype(np.int8)
        validation_metrics = metric_dict(y_eval[primary_validation], validation_pred, validation_score)
        nk_stress_score = score_model(model_object, x[nk_stress])
        tcrab_stress_score = score_model(model_object, x[tcrab_stress])
        nk_stress_predicted = int((nk_stress_score >= threshold).sum())
        tcrab_stress_predicted = int((tcrab_stress_score >= threshold).sum())
        candidates.append(
            Candidate(
                model=name,
                feature_set=feature_set,
                spec=spec,
                model_object=model_object,
                threshold=float(threshold),
                threshold_valid=valid,
                threshold_reason=reason,
                primary_tune_metrics=tune_metrics,
                nk_tune_fpr=float(nk_tune_fpr),
                nk_tune_predicted=int(nk_tune_predicted),
                tcrab_tune_fpr=float(tcrab_tune_fpr),
                tcrab_tune_predicted=int(tcrab_tune_predicted),
                current_primary_tune_f1=current_primary_tune_f1,
                primary_validation_metrics=validation_metrics,
                nk_stress_fpr=float(safe_div(nk_stress_predicted, nk_stress_score.size)),
                nk_stress_predicted=int(nk_stress_predicted),
                current_nk_stress_fpr=float(current_nk_fpr),
                current_nk_stress_predicted=int(current_nk_predicted),
                tcrab_stress_fpr=float(safe_div(tcrab_stress_predicted, tcrab_stress_score.size)),
                tcrab_stress_predicted=int(tcrab_stress_predicted),
                current_tcrab_stress_fpr=float(current_tcrab_fpr),
                current_tcrab_stress_predicted=int(current_tcrab_predicted),
                notes=notes,
            )
        )

    add_candidate(
        "current_model_threshold_only_multi_guarded",
        "base_TCR_controls",
        base_spec,
        current_model,
        x_base,
        "Current shared logistic model with a retuned threshold selected against both NK and TCRAB guard false positives.",
    )
    if args.threshold_only:
        return candidates

    logging.info("Fitting logistic_TCR_controls_multi_guard_negatives")
    logistic_base = make_pipeline(
        StandardScaler(),
        LogisticRegression(solver="lbfgs", class_weight="balanced", max_iter=500, random_state=args.seed),
    )
    logistic_base.fit(x_base[train_sample], y_eval[train_sample])
    add_candidate(
        "logistic_TCR_controls_multi_guard_negatives",
        "base_TCR_controls",
        base_spec,
        logistic_base,
        x_base,
        "Balanced logistic model using TCR plus FOXP3/CD4/CD3 controls and multi-guard negatives.",
    )

    logging.info("Fitting logistic_TCR_NKmarker_controls_multi_guard_negatives")
    logistic_nk = make_pipeline(
        StandardScaler(),
        LogisticRegression(solver="lbfgs", class_weight="balanced", max_iter=500, random_state=args.seed),
    )
    logistic_nk.fit(x_nk[train_sample], y_eval[train_sample])
    add_candidate(
        "logistic_TCR_NKmarker_controls_multi_guard_negatives",
        "TCR_NKmarker_controls",
        nk_spec,
        logistic_nk,
        x_nk,
        "Balanced logistic model using TCR, FOXP3/CD4/CD3, fixed NK marker controls, and multi-guard negatives.",
    )

    logging.info("Fitting logistic_TCR_NKmarker_controls_multi_guard_weighted")
    logistic_weighted = make_pipeline(
        StandardScaler(),
        LogisticRegression(solver="lbfgs", class_weight="balanced", max_iter=500, random_state=args.seed),
    )
    logistic_weighted.fit(x_nk[train_sample], y_eval[train_sample], logisticregression__sample_weight=train_weights)
    add_candidate(
        "logistic_TCR_NKmarker_controls_multi_guard_weighted",
        "TCR_NKmarker_controls",
        nk_spec,
        logistic_weighted,
        x_nk,
        "Weighted logistic model that upweights NK and TCRAB guard negatives during fitting.",
    )

    if args.skip_xgb:
        return candidates

    logging.info("Fitting xgb_TCR_NKmarker_controls_multi_guard_weighted")
    n_pos = int(y_eval[train_sample].sum())
    n_neg = int((y_eval[train_sample] == 0).sum())
    xgb = XGBClassifier(
        n_estimators=420,
        max_depth=4,
        learning_rate=0.045,
        subsample=0.85,
        colsample_bytree=0.85,
        min_child_weight=2.0,
        reg_lambda=1.5,
        objective="binary:logistic",
        eval_metric="logloss",
        tree_method="hist",
        n_jobs=32,
        random_state=args.seed,
        scale_pos_weight=max(n_neg / max(n_pos, 1), 1.0),
    )
    xgb.fit(x_nk[train_sample], y_eval[train_sample], sample_weight=train_weights)
    add_candidate(
        "xgb_TCR_NKmarker_controls_multi_guard_weighted",
        "TCR_NKmarker_controls",
        nk_spec,
        xgb,
        x_nk,
        "Weighted XGBoost model using TCR, FOXP3/CD4/CD3, fixed NK marker controls, and multi-guard negatives.",
    )
    return candidates


def candidate_table(candidates: list[Candidate]) -> pd.DataFrame:
    rows = []
    for candidate in candidates:
        m = candidate.primary_validation_metrics
        nk_reduction = 1.0 - safe_div(candidate.nk_stress_fpr, candidate.current_nk_stress_fpr)
        tcrab_reduction = 1.0 - safe_div(candidate.tcrab_stress_fpr, candidate.current_tcrab_stress_fpr)
        rows.append(
            {
                "model": candidate.model,
                "feature_set": candidate.feature_set,
                "threshold": candidate.threshold,
                "threshold_valid": candidate.threshold_valid,
                "threshold_reason": candidate.threshold_reason,
                "current_primary_tune_f1": candidate.current_primary_tune_f1,
                "primary_tune_f1": candidate.primary_tune_metrics["f1"],
                "primary_tune_precision": candidate.primary_tune_metrics["precision"],
                "primary_tune_recall": candidate.primary_tune_metrics["recall"],
                "nk_tune_fpr": candidate.nk_tune_fpr,
                "nk_tune_predicted": candidate.nk_tune_predicted,
                "tcrab_tune_fpr": candidate.tcrab_tune_fpr,
                "tcrab_tune_predicted": candidate.tcrab_tune_predicted,
                "primary_validation_f1": m["f1"],
                "primary_validation_precision": m["precision"],
                "primary_validation_recall": m["recall"],
                "primary_validation_specificity": m["specificity"],
                "primary_validation_tp": m["tp"],
                "primary_validation_fp": m["fp"],
                "primary_validation_tn": m["tn"],
                "primary_validation_fn": m["fn"],
                "nk_stress_fpr": candidate.nk_stress_fpr,
                "nk_stress_predicted": candidate.nk_stress_predicted,
                "current_nk_stress_fpr": candidate.current_nk_stress_fpr,
                "current_nk_stress_predicted": candidate.current_nk_stress_predicted,
                "nk_stress_fpr_reduction_vs_current": nk_reduction,
                "tcrab_stress_fpr": candidate.tcrab_stress_fpr,
                "tcrab_stress_predicted": candidate.tcrab_stress_predicted,
                "current_tcrab_stress_fpr": candidate.current_tcrab_stress_fpr,
                "current_tcrab_stress_predicted": candidate.current_tcrab_stress_predicted,
                "tcrab_stress_fpr_reduction_vs_current": tcrab_reduction,
                "notes": candidate.notes,
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["threshold_valid", "primary_tune_f1", "tcrab_tune_fpr", "nk_tune_fpr"], ascending=[False, False, True, True]
    ).reset_index(drop=True)


def select_tune_candidate(candidates: list[Candidate]) -> Candidate:
    valid = [c for c in candidates if c.threshold_valid]
    pool = valid if valid else candidates
    return max(pool, key=lambda c: (bool(c.threshold_valid), c.primary_tune_metrics["f1"], -c.tcrab_tune_fpr, -c.nk_tune_fpr, c.primary_tune_metrics["precision"]))


def full_context(handle: h5py.File) -> dict[str, np.ndarray]:
    n_obs = int(handle["obs"]["_index"].shape[0])
    source = clean_group_values(read_obs_column(handle, "source_gse_id"))
    tissue_key = "tissue_corrected" if "tissue_corrected" in handle["obs"] else "tissue"
    tissue = clean_group_values(read_obs_column(handle, tissue_key))
    annotation = clean_group_values(read_obs_column(handle, "simple_annotation_plus6")) if "simple_annotation_plus6" in handle["obs"] else np.full(n_obs, "unknown", dtype=object)
    library_id = normalize_strings(read_obs_column(handle, "library_id"))
    sample_id = normalize_strings(read_obs_column(handle, "sample_id"))
    use_library = pd.Series(library_id, copy=False).astype(str)
    use_sample = pd.Series(sample_id, copy=False).astype(str)
    sublibrary = use_library.mask(use_library.str.lower().isin({"", "na", "nan", "none", "null"}), use_sample).fillna("")
    sublibrary = sublibrary.mask(sublibrary.str.lower().isin({"", "na", "nan", "none", "null"}), "unknown").astype(str).to_numpy(dtype=object)
    has_TRA_TRB_paired = read_bool_obs(handle, "has_TRA_TRB_paired")
    has_any_ab_tcr = read_bool_obs(handle, "has_any_ab_tcr") if "has_any_ab_tcr" in handle["obs"] else has_TRA_TRB_paired.copy()
    tra_nonempty = read_nonempty_if_present(handle, "TRA_cdr3", n_obs)
    trb_nonempty = read_nonempty_if_present(handle, "TRB_cdr3", n_obs)
    ab_evidence = has_TRA_TRB_paired | has_any_ab_tcr | tra_nonempty | trb_nonempty
    keys = pd.Series(source, copy=False).astype(str) + "||" + pd.Series(sublibrary, copy=False).astype(str)
    tcrab_seq_keys = set(keys[ab_evidence])
    is_tcrab_seq_library = keys.isin(tcrab_seq_keys).to_numpy(dtype=bool)
    is_nk = pd.Series(annotation, copy=False).astype(str).str.upper().str.contains("NK", regex=False).to_numpy(dtype=bool)
    return {
        "source": source,
        "tissue": tissue,
        "annotation": annotation,
        "sublibrary": sublibrary,
        "has_TRA_TRB_paired": has_TRA_TRB_paired,
        "has_any_ab_evidence": ab_evidence,
        "is_tcrab_seq_library": is_tcrab_seq_library,
        "is_nk": is_nk,
    }


def apply_candidate_full(handle: h5py.File, candidate: Candidate) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, np.ndarray]:
    n_obs = int(handle["obs"]["_index"].shape[0])
    pred = np.zeros(n_obs, dtype=bool)
    for start in range(0, n_obs, FULL_APPLY_CHUNK):
        end = min(start + FULL_APPLY_CHUNK, n_obs)
        rows = np.arange(start, end, dtype=np.int64)
        x = extract_gene_matrix(handle, rows, candidate.spec, label=f"multiguard_full_{start}_{end}")
        score = score_model(candidate.model_object, x)
        pred[start:end] = score >= candidate.threshold
        if start and start % 500_000 == 0:
            logging.info("Applied multi-guard candidate to %s / %s cells", f"{start:,}", f"{n_obs:,}")
    ctx = full_context(handle)
    outputs = summarize_full_predictions(
        strategy=f"multiguard_candidate_{candidate.model}",
        threshold=float(candidate.threshold),
        pred=pred,
        source=ctx["source"],
        tissue=ctx["tissue"],
        annotation=ctx["annotation"],
        sublibrary=ctx["sublibrary"],
        is_tcrab_seq_library=ctx["is_tcrab_seq_library"],
        has_TRA_TRB_paired=ctx["has_TRA_TRB_paired"],
        is_nk=ctx["is_nk"],
    )
    overall = outputs[0].copy()
    overall["predicted_any_TCRAB_evidence_FP"] = int((pred & ctx["has_any_ab_evidence"]).sum())
    overall["any_TCRAB_evidence_FP_fraction_of_predictions"] = safe_div(int((pred & ctx["has_any_ab_evidence"]).sum()), int(pred.sum()))
    return (overall, *outputs[1:], pred)


def combine_with_current_outputs(candidate_outputs: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    current_paths = {
        "overall": CURRENT_TABLE_DIR / "selected_model_full_dataset_prediction_overall.csv",
        "source": CURRENT_TABLE_DIR / "selected_model_full_dataset_prediction_by_source.csv",
        "tissue": CURRENT_TABLE_DIR / "selected_model_full_dataset_prediction_by_tissue.csv",
        "fp_overall": CURRENT_TABLE_DIR / "selected_model_full_dataset_false_positive_overall.csv",
        "fp_source": CURRENT_TABLE_DIR / "selected_model_full_dataset_false_positive_by_source.csv",
        "annotation": CURRENT_TABLE_DIR / "selected_model_full_dataset_prediction_by_annotation.csv",
    }
    current = {key: pd.read_csv(path) for key, path in current_paths.items()}
    cand_overall, cand_source, cand_tissue, cand_fp_source, _cand_sublib, cand_annotation = candidate_outputs
    cand_fp_overall = cand_overall[[col for col in current["fp_overall"].columns if col in cand_overall.columns]].copy()
    overall = pd.concat([cand_overall, current["overall"]], ignore_index=True, sort=False)
    source = pd.concat([cand_source, current["source"]], ignore_index=True, sort=False)
    tissue = pd.concat([cand_tissue, current["tissue"]], ignore_index=True, sort=False)
    fp_overall = pd.concat([cand_fp_overall, current["fp_overall"]], ignore_index=True, sort=False)
    fp_source = pd.concat([cand_fp_source, current["fp_source"]], ignore_index=True, sort=False)
    annotation = pd.concat([cand_annotation, current["annotation"]], ignore_index=True, sort=False)
    overall.to_csv(TABLE_DIR / "full_dataset_prediction_overall.csv", index=False)
    source.to_csv(TABLE_DIR / "full_dataset_prediction_by_source.csv", index=False)
    tissue.to_csv(TABLE_DIR / "full_dataset_prediction_by_tissue.csv", index=False)
    fp_overall.to_csv(TABLE_DIR / "full_dataset_false_positive_overall.csv", index=False)
    fp_source.to_csv(TABLE_DIR / "full_dataset_false_positive_by_source.csv", index=False)
    annotation.to_csv(TABLE_DIR / "full_dataset_prediction_by_annotation.csv", index=False)
    return overall, source, tissue, fp_overall, fp_source, annotation


def evaluate_acceptance(candidate: Candidate, fp_overall: pd.DataFrame) -> pd.DataFrame:
    cand_row = fp_overall.loc[fp_overall["strategy"].astype(str).str.startswith("multiguard_candidate_")].iloc[0]
    m = candidate.primary_validation_metrics
    nk_reduction = 1.0 - safe_div(candidate.nk_stress_fpr, candidate.current_nk_stress_fpr)
    tcrab_reduction = 1.0 - safe_div(candidate.tcrab_stress_fpr, candidate.current_tcrab_stress_fpr)
    checks = {
        "threshold_valid": bool(candidate.threshold_valid),
        "primary_validation_f1_ge_0p887": m["f1"] >= ACCEPT_PRIMARY_F1_MIN,
        "primary_validation_recall_ge_0p84": m["recall"] >= ACCEPT_PRIMARY_RECALL_MIN,
        "primary_validation_precision_ge_0p94": m["precision"] >= ACCEPT_PRIMARY_PRECISION_MIN,
        "nk_stress_fpr_reduction_ge_30pct": nk_reduction >= ACCEPT_NK_STRESS_REDUCTION_MIN,
        "full_nk_fp_le_4500": int(cand_row["predicted_NK_FP"]) <= ACCEPT_FULL_NK_FP_MAX,
        "full_paired_tcrab_fp_le_2600": int(cand_row["predicted_paired_TCRAB_FP"]) <= ACCEPT_FULL_PAIRED_TCRAB_FP_MAX,
        "full_known_fp_fraction_le_2p5pct": float(cand_row["known_FP_fraction_of_predictions"]) <= ACCEPT_FULL_KNOWN_FP_FRACTION_MAX,
    }
    row = {
        "candidate_model": candidate.model,
        "accepted_for_promotion": all(checks.values()),
        "primary_validation_f1": m["f1"],
        "primary_validation_precision": m["precision"],
        "primary_validation_recall": m["recall"],
        "primary_validation_specificity": m["specificity"],
        "nk_stress_fpr": candidate.nk_stress_fpr,
        "current_nk_stress_fpr": candidate.current_nk_stress_fpr,
        "nk_stress_fpr_reduction_vs_current": nk_reduction,
        "tcrab_stress_fpr": candidate.tcrab_stress_fpr,
        "current_tcrab_stress_fpr": candidate.current_tcrab_stress_fpr,
        "tcrab_stress_fpr_reduction_vs_current": tcrab_reduction,
        "full_predicted_putative_gdT": int(cand_row["predicted_putative_gdT"]),
        "full_predicted_NK_FP": int(cand_row["predicted_NK_FP"]),
        "full_predicted_paired_TCRAB_FP": int(cand_row["predicted_paired_TCRAB_FP"]),
        "full_known_FP_fraction": float(cand_row["known_FP_fraction_of_predictions"]),
        **checks,
    }
    out = pd.DataFrame([row])
    out.to_csv(TABLE_DIR / "promotion_acceptance_gates.csv", index=False)
    return out


def save_model(candidate: Candidate, accepted: bool) -> tuple[Path, Path | None]:
    payload = {
        "model": candidate.model,
        "threshold": candidate.threshold,
        "feature_names": candidate.spec.feature_names,
        "gene_names": candidate.spec.gene_names,
        "model_object": candidate.model_object,
        "notes": candidate.notes,
        "feature_set": candidate.feature_set,
        "nk_guard": True,
        "tcrab_guard": True,
        "multi_guard": True,
        "accepted_for_promotion": bool(accepted),
    }
    best_path = MODEL_DIR / "best_candidate_model.pkl"
    with best_path.open("wb") as out:
        pickle.dump(payload, out)
    selected_path = None
    if accepted:
        selected_path = MODEL_DIR / "selected_model.pkl"
        shutil.copyfile(best_path, selected_path)
    return best_path, selected_path


def plot_candidate_summary(candidate_df: pd.DataFrame, acceptance: pd.DataFrame, fp_overall: pd.DataFrame) -> list[Path]:
    paths: list[Path] = []
    fig, ax = plt.subplots(figsize=(10.2, 5.0), constrained_layout=True)
    plot = candidate_df.copy()
    x = np.arange(plot.shape[0])
    ax.bar(x - 0.25, plot["primary_validation_f1"], width=0.25, label="validation F1", color="#0f766e")
    ax.bar(x, 1.0 - plot["nk_stress_fpr"], width=0.25, label="1 - NK stress FPR", color="#1d4ed8")
    ax.bar(x + 0.25, 1.0 - plot["tcrab_stress_fpr"], width=0.25, label="1 - TCRAB stress FPR", color="#7c3aed")
    ax.set_xticks(x)
    ax.set_xticklabels(plot["model"], rotation=25, ha="right", fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_title("Multi-guard candidate validation summary")
    ax.legend(frameon=False)
    path = FIGURE_DIR / "candidate_validation_summary.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(7.8, 4.8), constrained_layout=True)
    comp = fp_overall.copy()
    comp["label"] = comp["strategy"].map(lambda v: "Multi-guard" if str(v).startswith("multiguard_candidate_") else ("Current gdTAI" if str(v).startswith("selected_model_") else "Original TRD-TRAB"))
    comp = comp.drop_duplicates("label")
    x = np.arange(comp.shape[0])
    ax.bar(x - 0.18, comp["predicted_NK_FP"], width=0.36, label="NK FP", color="#0f766e")
    ax.bar(x + 0.18, comp["predicted_paired_TCRAB_FP"], width=0.36, label="paired TCRAB FP", color="#b45309")
    ax.set_xticks(x)
    ax.set_xticklabels(comp["label"])
    ax.set_ylabel("False-positive count")
    ax.set_title("Full-atlas known false positives")
    ax.legend(frameon=False)
    path = FIGURE_DIR / "full_atlas_known_fp_comparison.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(5.0, 4.6), constrained_layout=True)
    accepted = bool(acceptance["accepted_for_promotion"].iloc[0])
    labels = ["Accepted", "Not accepted"]
    values = [1 if accepted else 0, 0 if accepted else 1]
    ax.bar(labels, values, color=["#0f766e", "#be123c"])
    ax.set_ylim(0, 1.2)
    ax.set_title("Promotion gate result")
    path = FIGURE_DIR / "promotion_gate_result.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)
    return paths


def plot_confusion_matrix(candidate: Candidate) -> Path:
    m = candidate.primary_validation_metrics
    cm = np.asarray([[m["tn"], m["fp"]], [m["fn"], m["tp"]]], dtype=int)
    fig, ax = plt.subplots(figsize=(4.8, 4.2), constrained_layout=True)
    im = ax.imshow(cm, cmap="Blues")
    ax.set_xticks([0, 1], labels=["Pred abT", "Pred gdT"])
    ax.set_yticks([0, 1], labels=["True abT", "True gdT"])
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f"{cm[i, j]:,}", ha="center", va="center", fontsize=13)
    ax.set_title(candidate.model)
    fig.colorbar(im, ax=ax, shrink=0.8)
    path = FIGURE_DIR / "selected_candidate_confusion_matrix.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_trd_trab_scatter(handle: h5py.File, candidate: Candidate, candidate_pred_full: np.ndarray, current_payload: dict[str, Any], base_spec: FeatureSpec, seed: int, sample_cells: int) -> Path:
    n_obs = int(handle["obs"]["_index"].shape[0])
    rng = np.random.default_rng(seed)
    random_idx = rng.choice(n_obs, size=min(sample_cells, n_obs), replace=False)
    candidate_idx = np.flatnonzero(candidate_pred_full)
    if candidate_idx.size > 250_000:
        candidate_idx = rng.choice(candidate_idx, size=250_000, replace=False)
    sample_idx = np.unique(np.concatenate([random_idx, candidate_idx])).astype(np.int64)
    logging.info("Building multi-guard TRD-vs-TRAB scatter sample with %s cells", f"{sample_idx.size:,}")
    x_candidate = extract_gene_matrix(handle, sample_idx, candidate.spec, label="multiguard_trd_trab_scatter")
    candidate_score = score_model(candidate.model_object, x_candidate)
    candidate_pred = candidate_score >= candidate.threshold
    if candidate.spec.gene_names[: len(base_spec.gene_names)] == base_spec.gene_names:
        x_current = x_candidate[:, : len(base_spec.gene_names)]
    else:
        x_current = extract_gene_matrix(handle, sample_idx, base_spec, label="current_trd_trab_scatter")
    current_score = score_model(current_payload["model_object"], x_current)
    current_pred = current_score >= float(current_payload["threshold"])
    trd = read_float_obs(handle, "phase4_trd_score")[sample_idx]
    trab = read_float_obs(handle, "phase4_trab_score")[sample_idx]
    category = np.full(sample_idx.size, "neither", dtype=object)
    category[current_pred & (~candidate_pred)] = "current only"
    category[candidate_pred & (~current_pred)] = "multi-guard only"
    category[current_pred & candidate_pred] = "both"
    df = pd.DataFrame({"phase4_trab_score": trab, "phase4_trd_score": trd, "category": category})
    df.to_csv(TABLE_DIR / "trd_vs_trab_current_vs_multiguard_scatter_sample.csv.gz", index=False)
    fig, ax = plt.subplots(figsize=(7.2, 6.0), constrained_layout=True)
    xlim, ylim = bounded_scatter_limits(trab, trd)
    colors = {"neither": "#c4c9ce", "current only": "#d97904", "multi-guard only": "#2b7c9b", "both": "#7a3db8"}
    for cat in ["neither", "current only", "multi-guard only", "both"]:
        idx = np.flatnonzero(category == cat)
        if idx.size == 0:
            continue
        max_points = 140_000 if cat == "neither" else 90_000
        if idx.size > max_points:
            idx = rng.choice(idx, size=max_points, replace=False)
        ax.scatter(trab[idx], trd[idx], s=2.0 if cat == "neither" else 4.0, alpha=0.18 if cat == "neither" else 0.7, linewidths=0, rasterized=True, color=colors[cat], label=f"{cat} (n={(category == cat).sum():,})")
    original_threshold = load_original_trd_trab_threshold()
    line_x = np.linspace(xlim[0], xlim[1], 100)
    ax.plot(line_x, line_x + original_threshold, color="#111111", lw=1.0, ls="--", label="original TRD-TRAB cutoff")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("Phase 4 TRAB score")
    ax.set_ylabel("Phase 4 TRD score")
    ax.set_title("TRD vs TRAB score space: current gdTAI vs multi-guard")
    ax.legend(frameon=False, fontsize=7)
    path = FIGURE_DIR / "trd_vs_trab_current_vs_multiguard.png"
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def table_html(df: pd.DataFrame, max_rows: int = 30) -> str:
    if df.empty:
        return "<p><em>No rows.</em></p>"
    return dataframe_to_html(df.head(max_rows))


def write_report(candidate: Candidate, candidate_df: pd.DataFrame, acceptance: pd.DataFrame, fp_overall: pd.DataFrame, fp_source: pd.DataFrame, split_overall: pd.DataFrame, figures: list[Path], best_path: Path, selected_path: Path | None) -> None:
    accepted = bool(acceptance["accepted_for_promotion"].iloc[0])
    cand_full = fp_overall.loc[fp_overall["strategy"].astype(str).str.startswith("multiguard_candidate_")].iloc[0]
    current_full = fp_overall.loc[fp_overall["strategy"].astype(str).str.startswith("selected_model_")].iloc[0]
    nk_delta = int(cand_full["predicted_NK_FP"] - current_full["predicted_NK_FP"])
    tcrab_delta = int(cand_full["predicted_paired_TCRAB_FP"] - current_full["predicted_paired_TCRAB_FP"])
    summary_lines = [
        "# gdTAI Multi-Guard Candidate Report",
        "",
        f"- Candidate selected by tune rule: `{candidate.model}`",
        f"- Promotion accepted: `{accepted}`",
        f"- Candidate threshold: `{candidate.threshold:.8f}`",
        f"- Primary validation F1: `{candidate.primary_validation_metrics['f1']:.4f}`",
        f"- Primary validation precision/recall: `{candidate.primary_validation_metrics['precision']:.4f}` / `{candidate.primary_validation_metrics['recall']:.4f}`",
        f"- NK-stress FPR reduction vs current: `{1.0 - safe_div(candidate.nk_stress_fpr, candidate.current_nk_stress_fpr):.2%}`",
        f"- TCRAB-stress FPR reduction vs current: `{1.0 - safe_div(candidate.tcrab_stress_fpr, candidate.current_tcrab_stress_fpr):.2%}`",
        f"- Full-atlas NK FP: `{int(cand_full['predicted_NK_FP']):,}` vs current `{int(current_full['predicted_NK_FP']):,}` (delta `{nk_delta:,}`)",
        f"- Full-atlas paired-TCRAB FP: `{int(cand_full['predicted_paired_TCRAB_FP']):,}` vs current `{int(current_full['predicted_paired_TCRAB_FP']):,}` (delta `{tcrab_delta:,}`)",
        f"- Full-atlas predicted gdT: `{int(cand_full['predicted_putative_gdT']):,}`",
        f"- Best candidate model: `{best_path}`",
        f"- Promoted selected model: `{selected_path if selected_path else 'not written because acceptance gates failed'}`",
        "",
        "## Interpretation And Decision",
        "",
        "- This iteration directly uses both NK guard negatives and single-or-paired TCRAB guard negatives during training and threshold selection.",
        "- The candidate is promoted only if it keeps primary validation performance, NK FP, paired-TCRAB FP, and known-FP fraction inside the preset gates.",
        "- If promoted, it is still a separate multi-guard model artifact and does not overwrite the earlier shared model in-place.",
        "",
        "## Split Summary",
        dataframe_to_markdown(split_overall),
        "",
        "## Candidate Metrics",
        dataframe_to_markdown(candidate_df),
        "",
        "## Promotion Gates",
        dataframe_to_markdown(acceptance),
        "",
        "## Full-Atlas FP Audit",
        dataframe_to_markdown(fp_overall),
        "",
        "## Top Known FP Sources",
        dataframe_to_markdown(fp_source.sort_values(["predicted_NK_FP", "predicted_paired_TCRAB_FP"], ascending=False).head(30)),
        "",
    ]
    REPORT_MD.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    asset_dir = STATIC_DIR / "assets" / OUT_PREFIX
    asset_dir.mkdir(parents=True, exist_ok=True)
    fig_blocks = []
    for path in figures:
        target = asset_dir / path.name
        if path.exists() and path.resolve() != target.resolve():
            shutil.copyfile(path, target)
        fig_blocks.append(f"<section class='figure'><h3>{html.escape(path.stem.replace('_', ' '))}</h3><img src='assets/{OUT_PREFIX}/{html.escape(path.name)}'></section>")
    fig_html = "\n".join(fig_blocks)
    status_color = "#0f766e" if accepted else "#be123c"
    css = """
    @page{size:A4 landscape;margin:8mm} body{font-family:Arial,Helvetica,sans-serif;margin:18px;color:#1d252c;line-height:1.45;background:#f5f6f7} main{max-width:1320px;margin:auto} section{background:white;border:1px solid #d9dee4;padding:16px;margin:12px 0;break-inside:avoid} h1{font-size:29px;margin:0 0 8px} h2{font-size:20px} table{border-collapse:collapse;width:100%;font-size:9px;table-layout:fixed} th,td{border:1px solid #d9dee4;padding:3px 4px;text-align:left;vertical-align:top;overflow-wrap:anywhere} th{background:#eef1f4}.grid{display:grid;grid-template-columns:repeat(5,1fr);gap:10px}.metric{background:white;border:1px solid #d9dee4;padding:12px}.value{font-size:22px;font-weight:bold}.figure img{max-width:100%;height:auto;border:1px solid #d9dee4;background:white}.status{font-weight:bold;color:__STATUS_COLOR__} code{background:#eef1f4;padding:1px 4px;border-radius:3px}
    """.replace("__STATUS_COLOR__", status_color)
    html_doc = f"""<!doctype html><html><head><meta charset='utf-8'><title>gdTAI Multi-Guard Report</title><style>{css}</style></head><body><main>
    <section><h1>gdTAI Multi-Guard Candidate Report</h1><p class='status'>Promotion accepted: {accepted}</p><p>This report trains candidates with both NK guard negatives and single-or-paired TCRAB guard negatives. The H5AD is read-only.</p></section>
    <section><h2>Headline Metrics</h2><div class='grid'>
      <div class='metric'><div>Model</div><div class='value'>{html.escape(candidate.model)}</div></div>
      <div class='metric'><div>Primary F1</div><div class='value'>{candidate.primary_validation_metrics['f1']:.3f}</div></div>
      <div class='metric'><div>Full NK FP</div><div class='value'>{int(cand_full['predicted_NK_FP']):,}</div></div>
      <div class='metric'><div>Paired TCRAB FP</div><div class='value'>{int(cand_full['predicted_paired_TCRAB_FP']):,}</div></div>
      <div class='metric'><div>Predicted gdT</div><div class='value'>{int(cand_full['predicted_putative_gdT']):,}</div></div>
    </div></section>
    <section><h2>Interpretation And Decision</h2><p>The threshold selection now jointly constrains primary F1, NK-tune FPR, and single-or-paired TCRAB-tune FPR. Promotion still depends on held-out validation and full-atlas false-positive gates.</p></section>
    <section><h2>Promotion Gates</h2>{table_html(acceptance)}</section>
    <section><h2>Candidate Metrics</h2>{table_html(candidate_df, max_rows=20)}</section>
    <section><h2>Full-Atlas FP Audit</h2>{table_html(fp_overall, max_rows=10)}</section>
    <section><h2>Top Known FP Sources</h2>{table_html(fp_source.sort_values(['predicted_NK_FP','predicted_paired_TCRAB_FP'], ascending=False), max_rows=30)}</section>
    <section><h2>Split Summary</h2>{table_html(split_overall)}</section>
    <section><h2>Figures</h2>{fig_html}</section>
    </main></body></html>"""
    REPORT_HTML.write_text(html_doc, encoding="utf-8")


def render_pdf(no_pdf: bool) -> None:
    if no_pdf:
        return
    subprocess.run(
        [
            "google-chrome",
            "--headless",
            "--disable-gpu",
            "--no-sandbox",
            "--print-to-pdf-no-header",
            f"--print-to-pdf={REPORT_PDF.resolve()}",
            REPORT_HTML.resolve().as_uri(),
        ],
        check=True,
    )


def write_protocol(candidate: Candidate, accepted: bool, best_path: Path, selected_path: Path | None) -> Path:
    path = MODEL_DIR / "USAGE_PROTOCOL.md"
    text = f"""# gdTAI Multi-Guard Candidate Usage Protocol

This model was trained to control NK and TCRAB false positives jointly. NK annotation and TCRAB metadata are used for labels/audits, not inference features.

- Candidate model: `{candidate.model}`
- Feature set: `{candidate.feature_set}`
- Threshold: `{candidate.threshold:.16g}`
- Accepted for promotion: `{accepted}`
- Best candidate pickle: `{best_path}`
- Selected model pickle: `{selected_path if selected_path else 'not written; acceptance gates failed'}`

Use `workflows/gdtai/predict_with_selected_gdt_model.py --model-pkl <pickle>` for same-server testing. The wrapper reads the saved gene list and threshold from the pickle.

Required input assumptions are unchanged from the current gdTAI model: human gene symbols in `var/_index` and count-like `X` for `log1p(counts per 10,000)` reconstruction.
"""
    path.write_text(text, encoding="utf-8")
    return path


def main() -> None:
    args = parse_args()
    ensure_dirs()
    setup_logging()
    input_h5ad = args.input_h5ad.resolve()
    stat_before = input_h5ad.stat()
    with args.current_model_pkl.open("rb") as fh:
        current_payload = pickle.load(fh)
    logging.info("Input H5AD: %s", input_h5ad)
    logging.info("Current model: %s", args.current_model_pkl)
    with h5py.File(input_h5ad, "r") as handle:
        obs, splits = build_splits(handle, args.seed)
        base_spec = build_feature_spec(handle, include_nk_controls=False, name="base_TCR_controls")
        nk_spec = build_feature_spec(handle, include_nk_controls=True, name="TCR_NKmarker_controls")
        if [str(g) for g in current_payload["gene_names"]] != base_spec.gene_names:
            raise RuntimeError("Current model gene list does not match base feature spec; refusing to compare mismatched features.")
        eval_rows = np.unique(
            np.concatenate([
                splits.train_idx,
                splits.tune_idx,
                splits.validation_idx,
                splits.nk_train_idx,
                splits.nk_tune_idx,
                splits.nk_stress_idx,
                splits.tcrab_train_idx,
                splits.tcrab_tune_idx,
                splits.tcrab_stress_idx,
            ])
        ).astype(np.int64)
        positions = {
            "primary_train": local_positions(eval_rows, splits.train_idx),
            "primary_tune": local_positions(eval_rows, splits.tune_idx),
            "primary_validation": local_positions(eval_rows, splits.validation_idx),
            "nk_train": local_positions(eval_rows, splits.nk_train_idx),
            "nk_tune": local_positions(eval_rows, splits.nk_tune_idx),
            "nk_stress": local_positions(eval_rows, splits.nk_stress_idx),
            "tcrab_train": local_positions(eval_rows, splits.tcrab_train_idx),
            "tcrab_tune": local_positions(eval_rows, splits.tcrab_tune_idx),
            "tcrab_stress": local_positions(eval_rows, splits.tcrab_stress_idx),
        }
        y_eval = np.zeros(eval_rows.size, dtype=np.int8)
        y_eval[obs.class_code[eval_rows] == 2] = 1
        is_nk_eval = splits.nk_guard_mask[eval_rows]
        is_tcrab_eval = splits.tcrab_guard_mask[eval_rows]
        logging.info("Extracting multi-guard feature matrix for %s rows", f"{eval_rows.size:,}")
        x_nk = extract_gene_matrix(handle, eval_rows, nk_spec, label="multiguard_train_tune_validation")
        x_base = x_nk[:, : len(base_spec.gene_names)]
        candidates = train_candidates(x_base, x_nk, y_eval, is_nk_eval, is_tcrab_eval, positions, current_payload, base_spec, nk_spec, args)
        candidate_df = candidate_table(candidates)
        candidate_df.to_csv(TABLE_DIR / "candidate_tune_validation_metrics.csv", index=False)
        selected_candidate = select_tune_candidate(candidates)
        logging.info("Selected tune candidate: %s", selected_candidate.model)
        full_outputs = apply_candidate_full(handle, selected_candidate)
        candidate_full_tuple = full_outputs[:6]
        candidate_pred_full = full_outputs[6]
        _full_overall, _full_source, _full_tissue, fp_overall, fp_source, _full_annotation = combine_with_current_outputs(candidate_full_tuple)
        acceptance = evaluate_acceptance(selected_candidate, fp_overall)
        best_path, selected_path = save_model(selected_candidate, bool(acceptance["accepted_for_promotion"].iloc[0]))
        protocol_path = write_protocol(selected_candidate, bool(acceptance["accepted_for_promotion"].iloc[0]), best_path, selected_path)
        figures = plot_candidate_summary(candidate_df, acceptance, fp_overall)
        figures.append(plot_confusion_matrix(selected_candidate))
        figures.append(plot_trd_trab_scatter(handle, selected_candidate, candidate_pred_full, current_payload, base_spec, args.seed, args.scatter_sample_cells))
    stat_after = input_h5ad.stat()
    if (stat_before.st_size != stat_after.st_size) or (stat_before.st_mtime_ns != stat_after.st_mtime_ns):
        raise RuntimeError("Input H5AD changed during multi-guard read-only training/audit.")
    write_report(selected_candidate, candidate_df, acceptance, fp_overall, fp_source, splits.split_overall, figures, best_path, selected_path)
    render_pdf(args.no_pdf)
    (LOG_DIR / "gdtai_multiguard_summary.json").write_text(
        json.dumps(
            {
                "input_h5ad": str(input_h5ad),
                "candidate_model": selected_candidate.model,
                "threshold": selected_candidate.threshold,
                "accepted_for_promotion": bool(acceptance["accepted_for_promotion"].iloc[0]),
                "report_html": str(REPORT_HTML),
                "report_pdf": str(REPORT_PDF),
                "best_candidate_model": str(best_path),
                "selected_model": None if selected_path is None else str(selected_path),
                "usage_protocol": str(protocol_path),
            },
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )
    logging.info("Saved multi-guard report: %s", REPORT_HTML)
    logging.info("Saved multi-guard PDF: %s", REPORT_PDF)


if __name__ == "__main__":
    main()
