#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""GSE144469-holdout gdT classifier using individual TCR genes only.

This is a stricter overfitting audit than the primary classifier report:

- all GSE144469 gold-label cells are held out for validation
- Phase 4 TRD/TRAB score columns are not read and not used
- features are individual TCR alpha/beta/gamma/delta genes, plus explicit
  FOXP3/CD4/CD3 penalty-control genes
- TCR metadata is used only to define labels, not as model features
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
import math
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
from sklearn.metrics import (
    average_precision_score,
    confusion_matrix,
    fbeta_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

from run_gdt_prediction_package_evaluation import (
    DEFAULT_INPUT_H5AD,
    TABLE_DIR,
    FIGURE_DIR,
    LOG_DIR,
    STATIC_ASSET_DIR,
    STATIC_DIR,
    build_corrected_tcr_evidence,
    build_sublibrary_labels,
    clean_group_values,
    dataframe_to_html,
    dataframe_to_markdown,
    make_truth_labels,
    normalize_strings,
    read_bool_obs,
    read_float_obs,
    read_nonempty_string_mask,
    read_obs_column,
    read_string_dataset,
)
from run_gdt_deg_tcr_classifier_training import (
    PREVALENCE_SCENARIOS,
    json_ready,
    metric_dict,
    ppv_at_prevalence,
    safe_div,
    wilson_ci,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
OUT_PREFIX = "gse144469_holdout_tcrgene"
OUT_TABLE_DIR = TABLE_DIR / OUT_PREFIX
OUT_FIGURE_DIR = FIGURE_DIR / OUT_PREFIX
OUT_LOG_DIR = LOG_DIR / OUT_PREFIX
OUT_MODEL_DIR = OUTPUT_ROOT / "models" / "gdT_prediction_classifier" / OUT_PREFIX

REPORT_MD = OUT_LOG_DIR / f"{OUT_PREFIX}_report.md"
RUN_LOG = OUT_LOG_DIR / f"{OUT_PREFIX}.log"
REPORT_HTML = STATIC_DIR / f"{OUT_PREFIX}_report.html"
REPORT_PDF = STATIC_DIR / f"{OUT_PREFIX}_report.pdf"
SUMMARY_JSON = OUT_LOG_DIR / f"{OUT_PREFIX}_summary.json"
MODEL_PKL = OUT_MODEL_DIR / "selected_model.pkl"

TARGET_SUM = 10_000.0
RANDOM_SEED = 20260501
HOLDOUT_SOURCE = "GSE144469"
EXTRA_TRAB_HOLDOUT_SOURCE = "GSE254249"
GDT2020_SOURCE = "GDT_2020AUG_woCOV"
GDT2020_HOLDOUT_TISSUE = "cord blood"
SUBOPTIMAL_SOURCE = "GDTlung2023july_7p"
MAX_FEATURE_GENES = 300
MAX_NEGATIVE_TRAIN = 250_000
MIN_SPECIFICITY = 0.0
F_BETA = 1.0
FULL_APPLY_CHUNK = 50_000

TCR_PREFIXES = (
    "TRAV",
    "TRAJ",
    "TRAC",
    "TRBV",
    "TRBJ",
    "TRBC",
    "TRGV",
    "TRGJ",
    "TRGC",
    "TRDV",
    "TRDD",
    "TRDJ",
    "TRDC",
)
PENALTY_CONTROL_GENES = ["FOXP3", "CD4", "CD3D", "CD3E", "CD3G"]


@dataclass
class ObsMetadata:
    source: np.ndarray
    tissue: np.ndarray
    class_code: np.ndarray
    silver_mask: np.ndarray
    has_TRA_TRB_paired: np.ndarray
    corrected_has_any_gd_tcr: np.ndarray
    has_any_ab_tcr: np.ndarray
    sublibrary_summary: pd.DataFrame
    conflict_df: pd.DataFrame
    tcr_evidence_audit: pd.DataFrame


@dataclass
class FeatureSpec:
    feature_names: list[str]
    gene_names: list[str]
    gene_indices: np.ndarray
    tcr_gene_names: list[str]
    penalty_gene_names: list[str]
    family_by_feature: list[str]
    cd3_cols: list[int]
    cd4_col: int | None
    foxp3_col: int | None


@dataclass
class Result:
    model: str
    score_name: str
    threshold: float
    tune_metrics: dict[str, Any]
    validation_metrics: dict[str, Any]
    validation_score: np.ndarray
    validation_pred: np.ndarray
    tune_score: np.ndarray
    tune_pred: np.ndarray
    silver_recall: float | None
    model_object: Any | None = None
    notes: str = ""


@dataclass
class FullApplyOutputs:
    overall: pd.DataFrame
    source: pd.DataFrame
    tissue: pd.DataFrame
    fp_overall: pd.DataFrame
    fp_source: pd.DataFrame
    fp_sublibrary: pd.DataFrame
    fp_annotation: pd.DataFrame


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="GSE144469-holdout TCR-gene gdT classifier.")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--holdout-source", default=HOLDOUT_SOURCE)
    parser.add_argument("--extra-trab-holdout-source", default=EXTRA_TRAB_HOLDOUT_SOURCE)
    parser.add_argument("--gdt2020-holdout-tissue", default=GDT2020_HOLDOUT_TISSUE)
    parser.add_argument("--exclude-suboptimal-source", default=SUBOPTIMAL_SOURCE)
    parser.add_argument("--max-feature-genes", type=int, default=MAX_FEATURE_GENES)
    parser.add_argument("--max-negative-train", type=int, default=MAX_NEGATIVE_TRAIN)
    parser.add_argument(
        "--min-specificity",
        type=float,
        default=MIN_SPECIFICITY,
        help="Optional specificity floor for tune-split threshold selection. Default 0.0 selects the best F1 threshold.",
    )
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    parser.add_argument("--no-full-apply", action="store_true")
    parser.add_argument(
        "--apply-best-unaccepted-full",
        action="store_true",
        help="Exploratorily apply the best learned model even if no model passed the holdout acceptance rule.",
    )
    parser.add_argument(
        "--reuse-existing-full-outputs",
        action="store_true",
        help="Reuse existing exploratory full-dataset output CSVs when regenerating report/figures.",
    )
    parser.add_argument("--scatter-sample-cells", type=int, default=250_000)
    parser.add_argument("--no-pdf", action="store_true")
    return parser.parse_args()


def ensure_dirs() -> None:
    for path in [OUT_TABLE_DIR, OUT_FIGURE_DIR, OUT_LOG_DIR, OUT_MODEL_DIR, STATIC_DIR, STATIC_ASSET_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def setup_logging() -> None:
    ensure_dirs()
    handlers = [logging.FileHandler(RUN_LOG, mode="w", encoding="utf-8"), logging.StreamHandler()]
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
        force=True,
    )


def required_label_columns_present(handle: h5py.File) -> None:
    required = [
        "source_gse_id",
        "library_id",
        "sample_id",
        "has_TRA_TRB_paired",
        "has_TRG_TRD_paired",
        "has_any_ab_tcr",
        "has_any_gd_tcr",
        "TRG_cdr3",
        "TRD_cdr3",
    ]
    present = set(handle["obs"].keys())
    missing = [key for key in required if key not in present]
    if missing:
        raise KeyError(f"Missing required label columns: {missing}")


def read_nonempty_if_present(handle: h5py.File, key: str, n_obs: int) -> np.ndarray:
    if key not in handle["obs"]:
        return np.zeros(n_obs, dtype=bool)
    return read_nonempty_string_mask(handle, key)


def build_obs_metadata(handle: h5py.File) -> ObsMetadata:
    required_label_columns_present(handle)
    source = clean_group_values(read_obs_column(handle, "source_gse_id"))
    tissue_key = "tissue_corrected" if "tissue_corrected" in handle["obs"] else "tissue"
    tissue = clean_group_values(read_obs_column(handle, tissue_key))
    library_id = normalize_strings(read_obs_column(handle, "library_id"))
    sample_id = normalize_strings(read_obs_column(handle, "sample_id"))
    has_TRA_TRB_paired = read_bool_obs(handle, "has_TRA_TRB_paired")
    has_TRG_TRD_paired = read_bool_obs(handle, "has_TRG_TRD_paired")
    has_any_ab_tcr = read_bool_obs(handle, "has_any_ab_tcr")
    has_any_gd_tcr = read_bool_obs(handle, "has_any_gd_tcr")
    trg_nonempty = read_nonempty_string_mask(handle, "TRG_cdr3")
    trd_nonempty = read_nonempty_string_mask(handle, "TRD_cdr3")
    corrected_any_gd, corrected_paired_gd, audit = build_corrected_tcr_evidence(
        source,
        has_TRA_TRB_paired,
        has_TRG_TRD_paired,
        has_any_ab_tcr,
        has_any_gd_tcr,
        trg_nonempty,
        trd_nonempty,
    )
    is_gdtcr_sublibrary, sublibrary_summary = build_sublibrary_labels(
        source,
        library_id,
        sample_id,
        corrected_any_gd,
        trg_nonempty,
        trd_nonempty,
    )
    class_code, conflict_df = make_truth_labels(
        source,
        is_gdtcr_sublibrary,
        has_TRA_TRB_paired,
        corrected_paired_gd,
        has_any_ab_tcr,
        corrected_any_gd,
        trd_nonempty,
    )
    silver_mask = class_code == 3
    return ObsMetadata(
        source=source,
        tissue=tissue,
        class_code=class_code,
        silver_mask=silver_mask,
        has_TRA_TRB_paired=has_TRA_TRB_paired,
        corrected_has_any_gd_tcr=corrected_any_gd,
        has_any_ab_tcr=has_any_ab_tcr,
        sublibrary_summary=sublibrary_summary,
        conflict_df=conflict_df,
        tcr_evidence_audit=audit,
    )


def tcr_family(gene: str) -> str:
    upper = gene.upper()
    if upper.startswith(("TRAV", "TRAJ", "TRAC")):
        return "TRA"
    if upper.startswith(("TRBV", "TRBJ", "TRBC")):
        return "TRB"
    if upper.startswith(("TRGV", "TRGJ", "TRGC")):
        return "TRG"
    if upper.startswith(("TRDV", "TRDD", "TRDJ", "TRDC")):
        return "TRD"
    if upper in {"CD3D", "CD3E", "CD3G"}:
        return "CD3_control"
    if upper in {"FOXP3", "CD4"}:
        return "death_penalty_control"
    return "other"


def tcr_priority(gene: str) -> tuple[int, str]:
    upper = gene.upper()
    if upper.startswith(("TRDC", "TRDV", "TRDD", "TRDJ")):
        return (0, upper)
    if upper.startswith(("TRGC", "TRGV", "TRGJ")):
        return (1, upper)
    if upper.startswith(("TRAC", "TRAV", "TRAJ")):
        return (2, upper)
    if upper.startswith(("TRBC", "TRBV", "TRBJ")):
        return (3, upper)
    return (4, upper)


def build_feature_spec(handle: h5py.File, max_feature_genes: int) -> FeatureSpec:
    var_names = pd.Index(read_string_dataset(handle["var"]["_index"]), dtype="string")
    gene_to_idx = {str(gene): idx for idx, gene in enumerate(var_names.astype(str))}
    all_genes = var_names.astype(str).tolist()
    tcr_genes = [gene for gene in all_genes if gene.upper().startswith(TCR_PREFIXES)]
    tcr_genes = sorted(tcr_genes, key=tcr_priority)
    penalty_genes = [gene for gene in PENALTY_CONTROL_GENES if gene in gene_to_idx]
    selected = []
    for gene in tcr_genes + penalty_genes:
        if gene not in selected:
            selected.append(gene)
    if len(selected) > max_feature_genes:
        selected = selected[:max_feature_genes]
    feature_names = [f"{gene}_log1p_cp10k" for gene in selected]
    gene_indices = np.asarray([gene_to_idx[gene] for gene in selected], dtype=np.int32)
    family_by_feature = [tcr_family(gene) for gene in selected]
    name_to_col = {gene: i for i, gene in enumerate(selected)}
    cd3_cols = [name_to_col[gene] for gene in ["CD3D", "CD3E", "CD3G"] if gene in name_to_col]
    spec = FeatureSpec(
        feature_names=feature_names,
        gene_names=selected,
        gene_indices=gene_indices,
        tcr_gene_names=[gene for gene in selected if gene.upper().startswith(TCR_PREFIXES)],
        penalty_gene_names=[gene for gene in selected if gene in penalty_genes],
        family_by_feature=family_by_feature,
        cd3_cols=cd3_cols,
        cd4_col=name_to_col.get("CD4"),
        foxp3_col=name_to_col.get("FOXP3"),
    )
    feature_df = pd.DataFrame(
        {
            "feature": spec.feature_names,
            "gene": spec.gene_names,
            "family": spec.family_by_feature,
            "feature_index": np.arange(len(spec.feature_names), dtype=int),
        }
    )
    feature_df.to_csv(OUT_TABLE_DIR / "feature_manifest.csv", index=False)
    return spec


def extract_gene_matrix(
    handle: h5py.File,
    rows: np.ndarray,
    spec: FeatureSpec,
    *,
    label: str,
) -> np.ndarray:
    rows = np.asarray(rows, dtype=np.int64)
    matrix = np.zeros((rows.size, len(spec.feature_names)), dtype=np.float32)
    gene_to_col = {int(gene_idx): col for col, gene_idx in enumerate(spec.gene_indices)}
    selected_set = set(gene_to_col)
    x_group = handle["X"]
    indptr_ds = x_group["indptr"]
    indices_ds = x_group["indices"]
    data_ds = x_group["data"]
    for out_row, obs_idx in enumerate(rows):
        start = int(indptr_ds[obs_idx])
        end = int(indptr_ds[obs_idx + 1])
        if end <= start:
            continue
        row_indices = indices_ds[start:end].astype(np.int32, copy=False)
        row_data = data_ds[start:end].astype(np.float32, copy=False)
        row_sum = float(np.sum(row_data, dtype=np.float64))
        if row_sum <= 0:
            continue
        row_values = np.log1p(row_data * (TARGET_SUM / row_sum)).astype(np.float32, copy=False)
        present = np.isin(row_indices, list(selected_set), assume_unique=False)
        if present.any():
            for gene_idx, value in zip(row_indices[present], row_values[present]):
                matrix[out_row, gene_to_col[int(gene_idx)]] = value
        if out_row and out_row % 50_000 == 0:
            logging.info("Extracted %s features for %s / %s rows", label, f"{out_row:,}", f"{rows.size:,}")
    return matrix


def make_train_tune_validation_splits(
    obs: ObsMetadata,
    *,
    holdout_source: str,
    extra_trab_holdout_source: str,
    gdt2020_holdout_tissue: str,
    exclude_suboptimal_source: str,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, pd.DataFrame]:
    primary = (obs.class_code == 1) | (obs.class_code == 2)
    validation_primary_idx = np.flatnonzero(primary & (obs.source == holdout_source)).astype(np.int64)
    validation_gdt2020_idx = np.flatnonzero(
        primary
        & (obs.source == GDT2020_SOURCE)
        & (obs.class_code == 2)
        & (pd.Series(obs.tissue, copy=False).astype(str).str.lower().to_numpy() == gdt2020_holdout_tissue.lower())
    ).astype(np.int64)
    validation_extra_trab_idx = np.flatnonzero(
        (obs.source == extra_trab_holdout_source)
        & obs.has_TRA_TRB_paired
        & (~obs.corrected_has_any_gd_tcr)
    ).astype(np.int64)
    validation_idx = np.unique(
        np.concatenate([validation_primary_idx, validation_gdt2020_idx, validation_extra_trab_idx])
    ).astype(np.int64)
    validation_mask = np.zeros(obs.source.shape[0], dtype=bool)
    validation_mask[validation_idx] = True
    excluded_idx = np.flatnonzero(primary & (obs.source == exclude_suboptimal_source)).astype(np.int64)
    train_tune_pool = np.flatnonzero(
        primary
        & (~validation_mask)
        & (obs.source != holdout_source)
        & (obs.source != exclude_suboptimal_source)
        & (obs.source != extra_trab_holdout_source)
    ).astype(np.int64)
    if validation_primary_idx.size == 0:
        raise RuntimeError(f"No primary gold cells found in holdout source {holdout_source}")
    if validation_gdt2020_idx.size == 0:
        raise RuntimeError(f"No {GDT2020_SOURCE} gdT_gold cells found with tissue `{gdt2020_holdout_tissue}`")
    if validation_extra_trab_idx.size == 0:
        raise RuntimeError(f"No paired TRA/TRB no-gdTCR cells found in extra TCRAB holdout source {extra_trab_holdout_source}")
    validation_y = (obs.class_code[validation_idx] == 2).astype(np.int8)
    if validation_y.sum() == 0 or (validation_y == 0).sum() == 0:
        raise RuntimeError("Combined validation cohort must contain both gdT-positive and abT-negative cells.")
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
        if idx.size < 5:
            train_parts.append(idx)
            n_tune = 0
        else:
            n_tune = max(1, int(round(idx.size * 0.20)))
            tune_parts.append(idx[:n_tune])
            train_parts.append(idx[n_tune:])
        rows.append(
            {
                "source_gse_id": source,
                "label": label,
                "n_cells": int(idx.size),
                "train": int(idx.size - n_tune),
                "tune": int(n_tune),
            }
        )
    train_idx = np.concatenate(train_parts).astype(np.int64)
    tune_idx = np.concatenate(tune_parts).astype(np.int64)
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_TABLE_DIR / "train_tune_split_by_source_label.csv", index=False)
    overall = []
    for split, idx in [
        ("train", train_idx),
        ("tune", tune_idx),
        (f"validation_primary_{holdout_source}", validation_primary_idx),
        (f"validation_gdT_{GDT2020_SOURCE}_{gdt2020_holdout_tissue.replace(' ', '_')}", validation_gdt2020_idx),
        (f"validation_abT_{extra_trab_holdout_source}_paired_TCRAB_no_gdTCR", validation_extra_trab_idx),
        ("validation_combined", validation_idx),
        (f"sensitivity_excluded_{exclude_suboptimal_source}", excluded_idx),
    ]:
        labels = obs.class_code[idx]
        if split.startswith("validation_abT_") or split == "validation_combined":
            gdT_gold = int((labels == 2).sum())
            abT_gold = int(idx.size - gdT_gold)
        else:
            gdT_gold = int((labels == 2).sum())
            abT_gold = int((labels == 1).sum())
        overall.append(
            {
                "split": split,
                "n_cells": int(idx.size),
                "gdT_gold": gdT_gold,
                "abT_gold": abT_gold,
                "gdT_prevalence": safe_div(gdT_gold, int(idx.size)),
            }
        )
    pd.DataFrame(overall).to_csv(OUT_TABLE_DIR / "split_overall.csv", index=False)
    return train_idx, tune_idx, validation_idx, excluded_idx, summary


def local_positions(eval_rows: np.ndarray, target_rows: np.ndarray) -> np.ndarray:
    lookup = pd.Series(np.arange(eval_rows.size, dtype=np.int64), index=eval_rows)
    return lookup.loc[target_rows].to_numpy(dtype=np.int64)


def choose_threshold(
    y_true: np.ndarray,
    score: np.ndarray,
    *,
    min_specificity: float,
    beta: float,
    min_threshold: float | None = None,
) -> tuple[float, dict[str, Any]]:
    fpr, tpr, thresholds = roc_curve(y_true, score)
    n_pos = int(y_true.sum())
    n_neg = int((y_true == 0).sum())
    tp = np.rint(tpr * n_pos).astype(int)
    fp = np.rint(fpr * n_neg).astype(int)
    fn = n_pos - tp
    tn = n_neg - fp
    specificity = np.divide(tn, tn + fp, out=np.zeros_like(tn, dtype=float), where=(tn + fp) > 0)
    recall = np.divide(tp, tp + fn, out=np.zeros_like(tp, dtype=float), where=(tp + fn) > 0)
    precision = np.divide(tp, tp + fp, out=np.zeros_like(tp, dtype=float), where=(tp + fp) > 0)
    beta2 = beta * beta
    fbeta = np.divide(
        (1.0 + beta2) * precision * recall,
        beta2 * precision + recall,
        out=np.zeros_like(precision),
        where=(beta2 * precision + recall) > 0,
    )
    valid = np.isfinite(thresholds) & (specificity >= min_specificity) & ((tp + fp) > 0)
    if min_threshold is not None:
        valid &= thresholds >= min_threshold
    if not valid.any():
        valid = np.isfinite(thresholds) & ((tp + fp) > 0)
        if min_threshold is not None:
            valid &= thresholds >= min_threshold
    if not valid.any():
        threshold = float(np.nanmax(score) + 1e-6)
        pred = (score >= threshold).astype(np.int8)
        metrics = metric_dict(y_true, pred, score)
        metrics["threshold_fallback_used"] = True
        return threshold, metrics
    valid_idx = np.flatnonzero(valid)
    order = np.lexsort((thresholds[valid_idx], recall[valid_idx], precision[valid_idx], fbeta[valid_idx]))
    idx = int(valid_idx[order[-1]])
    threshold = float(thresholds[idx])
    pred = (score >= threshold).astype(np.int8)
    metrics = metric_dict(y_true, pred, score)
    metrics["threshold_fallback_used"] = False
    return threshold, metrics


def sample_training_positions(pos_train: np.ndarray, y_eval: np.ndarray, max_negative: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    positives = pos_train[y_eval[pos_train] == 1]
    negatives = pos_train[y_eval[pos_train] == 0]
    n_neg = min(negatives.size, max(max_negative, positives.size * 4))
    if negatives.size > n_neg:
        negatives = rng.choice(negatives, size=n_neg, replace=False)
    sampled = np.concatenate([positives, negatives])
    rng.shuffle(sampled)
    return sampled.astype(np.int64)


def apply_death_penalties(score: np.ndarray, x: np.ndarray, spec: FeatureSpec) -> np.ndarray:
    out = score.astype(np.float32, copy=True)
    penalty = np.zeros(out.size, dtype=bool)
    if spec.foxp3_col is not None:
        penalty |= x[:, spec.foxp3_col] > 0.25
    if spec.cd4_col is not None:
        penalty |= x[:, spec.cd4_col] > 0.75
    if spec.cd3_cols:
        penalty |= x[:, spec.cd3_cols].sum(axis=1) < 0.25
    out[penalty] = np.minimum(out[penalty], 0.03)
    return out


def train_models(
    x_eval: np.ndarray,
    y_eval: np.ndarray,
    pos_train: np.ndarray,
    pos_tune: np.ndarray,
    pos_validation: np.ndarray,
    pos_silver: np.ndarray,
    spec: FeatureSpec,
    *,
    max_negative_train: int,
    min_specificity: float,
    seed: int,
) -> list[Result]:
    y_tune = y_eval[pos_tune]
    y_validation = y_eval[pos_validation]
    train_sample = sample_training_positions(pos_train, y_eval, max_negative_train, seed)
    logging.info(
        "Training sample: %s cells, positives=%s, negatives=%s",
        f"{train_sample.size:,}",
        f"{int(y_eval[train_sample].sum()):,}",
        f"{int((y_eval[train_sample] == 0).sum()):,}",
    )
    results: list[Result] = []

    family = np.asarray(spec.family_by_feature, dtype=object)
    trd_cols = np.flatnonzero(family == "TRD")
    trg_cols = np.flatnonzero(family == "TRG")
    tra_cols = np.flatnonzero(family == "TRA")
    trb_cols = np.flatnonzero(family == "TRB")
    baseline_specs = [
        ("baseline_individual_TRD_gene_sum", x_eval[:, trd_cols].sum(axis=1) if trd_cols.size else np.zeros(x_eval.shape[0])),
        (
            "baseline_individual_gd_minus_ab_gene_sum",
            (x_eval[:, np.concatenate([trd_cols, trg_cols])].sum(axis=1) if (trd_cols.size + trg_cols.size) else 0)
            - (x_eval[:, np.concatenate([tra_cols, trb_cols])].sum(axis=1) if (tra_cols.size + trb_cols.size) else 0),
        ),
    ]
    for model_name, full_score in baseline_specs:
        tune_score = np.asarray(full_score[pos_tune], dtype=np.float32)
        validation_score = np.asarray(full_score[pos_validation], dtype=np.float32)
        silver_score = np.asarray(full_score[pos_silver], dtype=np.float32) if pos_silver.size else None
        threshold, tune_metrics = choose_threshold(y_tune, tune_score, min_specificity=min_specificity, beta=F_BETA)
        tune_pred = (tune_score >= threshold).astype(np.int8)
        validation_pred = (validation_score >= threshold).astype(np.int8)
        validation_metrics = metric_dict(y_validation, validation_pred, validation_score)
        silver_recall = float(np.mean(silver_score >= threshold)) if silver_score is not None and silver_score.size else None
        results.append(
            Result(
                model=model_name,
                score_name="individual_gene_heuristic",
                threshold=threshold,
                tune_metrics=tune_metrics,
                validation_metrics=validation_metrics,
                validation_score=validation_score,
                validation_pred=validation_pred,
                tune_score=tune_score,
                tune_pred=tune_pred,
                silver_recall=silver_recall,
                notes="No Phase 4 scores; threshold selected on non-holdout tune split.",
            )
        )

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
        random_state=seed,
        scale_pos_weight=max(n_neg / max(n_pos, 1), 1.0),
    )
    xgb.fit(x_eval[train_sample], y_eval[train_sample])
    for model_name, use_penalty in [
        ("xgb_individual_TCRABGD_genes", False),
        ("xgb_individual_TCRABGD_genes_with_FOXP3_CD4_lowCD3_penalty", True),
    ]:
        tune_score = xgb.predict_proba(x_eval[pos_tune])[:, 1].astype(np.float32)
        validation_score = xgb.predict_proba(x_eval[pos_validation])[:, 1].astype(np.float32)
        silver_score = xgb.predict_proba(x_eval[pos_silver])[:, 1].astype(np.float32) if pos_silver.size else None
        if use_penalty:
            tune_score = apply_death_penalties(tune_score, x_eval[pos_tune], spec)
            validation_score = apply_death_penalties(validation_score, x_eval[pos_validation], spec)
            silver_score = apply_death_penalties(silver_score, x_eval[pos_silver], spec) if silver_score is not None else None
        threshold, tune_metrics = choose_threshold(
            y_tune,
            tune_score,
            min_specificity=min_specificity,
            beta=F_BETA,
        )
        tune_pred = (tune_score >= threshold).astype(np.int8)
        validation_pred = (validation_score >= threshold).astype(np.int8)
        validation_metrics = metric_dict(y_validation, validation_pred, validation_score)
        silver_recall = float(np.mean(silver_score >= threshold)) if silver_score is not None and silver_score.size else None
        results.append(
            Result(
                model=model_name,
                score_name="predicted_probability",
                threshold=threshold,
                tune_metrics=tune_metrics,
                validation_metrics=validation_metrics,
                validation_score=validation_score,
                validation_pred=validation_pred,
                tune_score=tune_score,
                tune_pred=tune_pred,
                silver_recall=silver_recall,
                model_object=xgb,
                notes="Individual TCR A/B/G/D genes only; no Phase 4 scores or TCR metadata features."
                + (" FOXP3/CD4/low-CD3 probability death penalties applied." if use_penalty else ""),
            )
        )

    logistic = make_pipeline(
        StandardScaler(),
        LogisticRegression(solver="lbfgs", class_weight="balanced", max_iter=1000, random_state=seed),
    )
    logistic.fit(x_eval[train_sample], y_eval[train_sample])
    tune_score = logistic.predict_proba(x_eval[pos_tune])[:, 1].astype(np.float32)
    validation_score = logistic.predict_proba(x_eval[pos_validation])[:, 1].astype(np.float32)
    silver_score = logistic.predict_proba(x_eval[pos_silver])[:, 1].astype(np.float32) if pos_silver.size else None
    threshold, tune_metrics = choose_threshold(y_tune, tune_score, min_specificity=min_specificity, beta=F_BETA)
    tune_pred = (tune_score >= threshold).astype(np.int8)
    validation_pred = (validation_score >= threshold).astype(np.int8)
    validation_metrics = metric_dict(y_validation, validation_pred, validation_score)
    silver_recall = float(np.mean(silver_score >= threshold)) if silver_score is not None and silver_score.size else None
    results.append(
        Result(
            model="logistic_individual_TCRABGD_genes",
            score_name="predicted_probability",
            threshold=threshold,
            tune_metrics=tune_metrics,
            validation_metrics=validation_metrics,
            validation_score=validation_score,
            validation_pred=validation_pred,
            tune_score=tune_score,
            tune_pred=tune_pred,
            silver_recall=silver_recall,
            model_object=logistic,
            notes="Logistic fallback using individual TCR A/B/G/D genes only.",
        )
    )
    return results


def results_to_frame(results: list[Result], split: str) -> pd.DataFrame:
    rows = []
    for result in results:
        metrics = result.tune_metrics if split == "tune" else result.validation_metrics
        row = {
            "model": result.model,
            "score_name": result.score_name,
            "threshold": result.threshold,
            "silver_recall_at_threshold": result.silver_recall,
            "notes": result.notes,
        }
        row.update(metrics)
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["f1", "precision", "recall"], ascending=False).reset_index(drop=True)


def select_model(results: list[Result]) -> tuple[Result | None, Result]:
    baselines = [result for result in results if result.model.startswith("baseline_")]
    best_baseline = max(
        baselines,
        key=lambda item: (
            item.validation_metrics["f1"],
            item.validation_metrics["precision"],
            item.validation_metrics["recall"],
        ),
    )
    candidates = [result for result in results if not result.model.startswith("baseline_")]
    accepted = []
    rows = []
    for result in candidates:
        metrics = result.validation_metrics
        delta_f1 = float(metrics["f1"] - best_baseline.validation_metrics["f1"])
        delta_f05 = float(metrics["f0.5"] - best_baseline.validation_metrics["f0.5"])
        accepted_flag = delta_f1 > 0
        rows.append(
            {
                "model": result.model,
                "best_individual_gene_baseline": best_baseline.model,
                "validation_f1": metrics["f1"],
                "baseline_f1": best_baseline.validation_metrics["f1"],
                "delta_f1": delta_f1,
                "validation_f0.5": metrics["f0.5"],
                "baseline_f0.5": best_baseline.validation_metrics["f0.5"],
                "delta_f0.5": delta_f05,
                "validation_precision": metrics["precision"],
                "validation_recall": metrics["recall"],
                "validation_specificity": metrics["specificity"],
                "accepted": accepted_flag,
            }
        )
        if accepted_flag:
            accepted.append(result)
    acceptance = pd.DataFrame(rows)
    acceptance.to_csv(OUT_TABLE_DIR / "acceptance_vs_individual_gene_baseline.csv", index=False)
    selected = None
    if accepted:
        selected = max(
            accepted,
            key=lambda item: (
                item.validation_metrics["f1"],
                item.validation_metrics["precision"],
                item.validation_metrics["recall"],
            ),
        )
    return selected, best_baseline


def choose_exploratory_model(results: list[Result]) -> Result | None:
    candidates = [
        result
        for result in results
        if not result.model.startswith("baseline_") and result.model_object is not None
    ]
    if not candidates:
        return None
    return max(
        candidates,
        key=lambda item: (
            item.validation_metrics["f1"],
            item.validation_metrics["precision"],
            item.validation_metrics["recall"],
            item.validation_metrics["specificity"],
        ),
    )


def prevalence_table(results: list[Result]) -> pd.DataFrame:
    rows = []
    for result in results:
        m = result.validation_metrics
        tp, fp, tn, fn = int(m["tp"]), int(m["fp"]), int(m["tn"]), int(m["fn"])
        rec_low, _rec_high = wilson_ci(tp, tp + fn)
        spec_low, _spec_high = wilson_ci(tn, tn + fp)
        scenarios = sorted(set(PREVALENCE_SCENARIOS + [safe_div(tp + fn, tp + fp + tn + fn)]))
        for prev in scenarios:
            rows.append(
                {
                    "model": result.model,
                    "prevalence": prev,
                    "prevalence_percent": prev * 100.0,
                    "validation_recall": m["recall"],
                    "validation_specificity": m["specificity"],
                    "ppv_observed_at_prevalence": ppv_at_prevalence(m["recall"], m["specificity"], prev),
                    "ppv_conservative_wilson_at_prevalence": ppv_at_prevalence(rec_low, spec_low, prev),
                    "expected_false_positive_per_million_observed": (1.0 - m["specificity"]) * (1.0 - prev) * 1_000_000.0,
                    "expected_false_positive_per_million_conservative": (1.0 - spec_low) * (1.0 - prev) * 1_000_000.0,
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(OUT_TABLE_DIR / "prevalence_aware_ppv_scenarios.csv", index=False)
    return out


def write_validation_strata(selected: Result, validation_idx: np.ndarray, y_validation: np.ndarray, obs: ObsMetadata) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows_source = []
    rows_tissue = []
    for key, values, rows in [
        ("source_gse_id", obs.source[validation_idx], rows_source),
        ("tissue", obs.tissue[validation_idx], rows_tissue),
    ]:
        for value in sorted(pd.Series(values).astype(str).unique()):
            mask = values.astype(str) == value
            metrics = metric_dict(y_validation[mask], selected.validation_pred[mask], selected.validation_score[mask])
            rows.append({key: value, **metrics})
    source_df = pd.DataFrame(rows_source)
    tissue_df = pd.DataFrame(rows_tissue)
    source_df.to_csv(OUT_TABLE_DIR / "selected_validation_metrics_by_source.csv", index=False)
    tissue_df.to_csv(OUT_TABLE_DIR / "selected_validation_metrics_by_tissue.csv", index=False)
    return source_df, tissue_df


def sublibrary_ids(library_id: np.ndarray, sample_id: np.ndarray) -> np.ndarray:
    use_library = pd.Series(library_id, copy=False).astype("string").fillna("").str.strip()
    use_sample = pd.Series(sample_id, copy=False).astype("string").fillna("").str.strip()
    sublibrary = use_library.mask(use_library.str.lower().isin({"", "na", "nan", "none", "null"}), use_sample).fillna("")
    sublibrary = sublibrary.mask(sublibrary.str.lower().isin({"", "na", "nan", "none", "null"}), "unknown").astype(str)
    return sublibrary.to_numpy(dtype=object)


def load_original_trd_trab_threshold() -> float:
    selected_json = LOG_DIR / "gdT_prediction_selected_model.json"
    if selected_json.exists():
        info = json.loads(selected_json.read_text(encoding="utf-8"))
        if "selected_threshold_on_score" in info:
            return float(info["selected_threshold_on_score"])
        if "selected_raw_cutoff_for_original_score" in info:
            return float(info["selected_raw_cutoff_for_original_score"])
    selected_csv = TABLE_DIR / "selected_model_summary.csv"
    if selected_csv.exists():
        df = pd.read_csv(selected_csv)
        if "selected_threshold_on_score" in df and not df.empty:
            return float(df["selected_threshold_on_score"].iloc[0])
    return 0.37491393089294434


def add_prediction_ratios(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    out["predicted_fraction"] = out["predicted_putative_gdT"] / out["total_cells"].replace(0, np.nan)
    out["paired_TCRAB_FP_ratio_among_paired_TCRAB"] = (
        out["predicted_paired_TCRAB_FP"] / out["paired_TCRAB_cells"].replace(0, np.nan)
    )
    out["paired_TCRAB_FP_fraction_of_predictions"] = (
        out["predicted_paired_TCRAB_FP"] / out["predicted_putative_gdT"].replace(0, np.nan)
    )
    out["NK_FP_ratio_among_NK"] = out["predicted_NK_FP"] / out["NK_cells"].replace(0, np.nan)
    out["NK_FP_fraction_of_predictions"] = out["predicted_NK_FP"] / out["predicted_putative_gdT"].replace(0, np.nan)
    out["known_FP_fraction_of_predictions"] = (
        out["predicted_paired_TCRAB_or_NK_FP"] / out["predicted_putative_gdT"].replace(0, np.nan)
    )
    return out.replace([np.inf, -np.inf], np.nan)


def read_existing_full_outputs(output_prefix: str) -> FullApplyOutputs:
    paths = {
        "overall": OUT_TABLE_DIR / f"{output_prefix}_full_dataset_prediction_overall.csv",
        "source": OUT_TABLE_DIR / f"{output_prefix}_full_dataset_prediction_by_source.csv",
        "tissue": OUT_TABLE_DIR / f"{output_prefix}_full_dataset_prediction_by_tissue.csv",
        "fp_overall": OUT_TABLE_DIR / f"{output_prefix}_full_dataset_false_positive_overall.csv",
        "fp_source": OUT_TABLE_DIR / f"{output_prefix}_full_dataset_false_positive_by_source.csv",
        "fp_sublibrary": OUT_TABLE_DIR / f"{output_prefix}_full_dataset_false_positive_by_sublibrary.csv",
        "fp_annotation": OUT_TABLE_DIR / f"{output_prefix}_full_dataset_prediction_by_annotation.csv",
    }
    if not paths["overall"].exists():
        empty = pd.DataFrame()
        return FullApplyOutputs(empty, empty, empty, empty, empty, empty, empty)
    return FullApplyOutputs(**{key: pd.read_csv(path) if path.exists() else pd.DataFrame() for key, path in paths.items()})


def summarize_full_predictions(
    *,
    strategy: str,
    threshold: float,
    pred: np.ndarray,
    source: np.ndarray,
    tissue: np.ndarray,
    annotation: np.ndarray,
    sublibrary: np.ndarray,
    is_tcrab_seq_library: np.ndarray,
    has_TRA_TRB_paired: np.ndarray,
    is_nk: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    total_cells = int(pred.size)
    predicted = pred.astype(bool, copy=False)
    paired_fp = predicted & has_TRA_TRB_paired
    nk_fp = predicted & is_nk
    known_fp = paired_fp | nk_fp
    base = pd.DataFrame(
        {
            "source_gse_id": source,
            "tissue": tissue,
            "annotation": annotation,
            "sublibrary_id": sublibrary,
            "is_tcrab_seq_library": is_tcrab_seq_library,
            "has_TRA_TRB_paired": has_TRA_TRB_paired,
            "is_NK_cell": is_nk,
            "predicted_putative_gdT": predicted,
            "predicted_paired_TCRAB_FP": paired_fp,
            "predicted_NK_FP": nk_fp,
            "predicted_paired_TCRAB_or_NK_FP": known_fp,
        }
    )
    overall = pd.DataFrame(
        [
            {
                "strategy": strategy,
                "threshold": threshold,
                "total_cells": total_cells,
                "predicted_putative_gdT": int(predicted.sum()),
                "tcrab_seq_library_cells": int(is_tcrab_seq_library.sum()),
                "predicted_in_tcrab_seq_libraries": int((predicted & is_tcrab_seq_library).sum()),
                "paired_TCRAB_cells": int(has_TRA_TRB_paired.sum()),
                "predicted_paired_TCRAB_FP": int(paired_fp.sum()),
                "NK_cells": int(is_nk.sum()),
                "predicted_NK_FP": int(nk_fp.sum()),
                "predicted_paired_TCRAB_or_NK_FP": int(known_fp.sum()),
            }
        ]
    )
    overall = add_prediction_ratios(overall)

    def summarize_group(group_cols: list[str]) -> pd.DataFrame:
        grouped = (
            base.groupby(group_cols, dropna=False, as_index=False)
            .agg(
                total_cells=("predicted_putative_gdT", "size"),
                predicted_putative_gdT=("predicted_putative_gdT", "sum"),
                tcrab_seq_library_cells=("is_tcrab_seq_library", "sum"),
                paired_TCRAB_cells=("has_TRA_TRB_paired", "sum"),
                predicted_paired_TCRAB_FP=("predicted_paired_TCRAB_FP", "sum"),
                NK_cells=("is_NK_cell", "sum"),
                predicted_NK_FP=("predicted_NK_FP", "sum"),
                predicted_paired_TCRAB_or_NK_FP=("predicted_paired_TCRAB_or_NK_FP", "sum"),
            )
        )
        grouped.insert(0, "strategy", strategy)
        grouped.insert(1, "threshold", threshold)
        return add_prediction_ratios(grouped)

    source_df = summarize_group(["source_gse_id"]).sort_values(
        ["predicted_putative_gdT", "predicted_paired_TCRAB_or_NK_FP"], ascending=[False, False]
    )
    tissue_df = summarize_group(["tissue"]).sort_values(
        ["predicted_putative_gdT", "predicted_paired_TCRAB_or_NK_FP"], ascending=[False, False]
    )
    annotation_df = summarize_group(["annotation"]).sort_values("predicted_putative_gdT", ascending=False)
    sublibrary_df = summarize_group(["source_gse_id", "sublibrary_id"])
    sublibrary_df = sublibrary_df[
        (sublibrary_df["predicted_putative_gdT"] > 0) | (sublibrary_df["tcrab_seq_library_cells"] > 0)
    ].sort_values(["predicted_paired_TCRAB_FP", "predicted_putative_gdT"], ascending=[False, False])
    fp_source_df = source_df[
        (source_df["predicted_paired_TCRAB_or_NK_FP"] > 0) | (source_df["paired_TCRAB_cells"] > 0) | (source_df["NK_cells"] > 0)
    ].copy()
    return overall, source_df, tissue_df, fp_source_df, sublibrary_df, annotation_df


def apply_full_dataset(
    handle: h5py.File,
    selected: Result,
    spec: FeatureSpec,
    *,
    output_prefix: str,
    exploratory: bool,
) -> FullApplyOutputs:
    if selected.model_object is None:
        empty = pd.DataFrame()
        return FullApplyOutputs(empty, empty, empty, empty, empty, empty, empty)
    n_obs = int(handle["obs"]["_index"].shape[0])
    source = clean_group_values(read_obs_column(handle, "source_gse_id"))
    tissue_key = "tissue_corrected" if "tissue_corrected" in handle["obs"] else "tissue"
    tissue = clean_group_values(read_obs_column(handle, tissue_key))
    annotation_key = "simple_annotation_plus6" if "simple_annotation_plus6" in handle["obs"] else None
    annotation = clean_group_values(read_obs_column(handle, annotation_key)) if annotation_key else np.full(n_obs, "unknown", dtype=object)
    library_id = normalize_strings(read_obs_column(handle, "library_id"))
    sample_id = normalize_strings(read_obs_column(handle, "sample_id"))
    sublibrary = sublibrary_ids(library_id, sample_id)
    has_TRA_TRB_paired = read_bool_obs(handle, "has_TRA_TRB_paired")
    has_any_ab_tcr = read_bool_obs(handle, "has_any_ab_tcr") if "has_any_ab_tcr" in handle["obs"] else has_TRA_TRB_paired.copy()
    tra_nonempty = read_nonempty_if_present(handle, "TRA_cdr3", n_obs)
    trb_nonempty = read_nonempty_if_present(handle, "TRB_cdr3", n_obs)
    ab_evidence = has_TRA_TRB_paired | has_any_ab_tcr | tra_nonempty | trb_nonempty
    keys = pd.Series(source, copy=False).astype(str) + "||" + pd.Series(sublibrary, copy=False).astype(str)
    tcrab_seq_keys = set(keys[ab_evidence])
    is_tcrab_seq_library = keys.isin(tcrab_seq_keys).to_numpy(dtype=bool)
    is_nk = pd.Series(annotation, copy=False).astype(str).str.upper().str.contains("NK", regex=False).to_numpy(dtype=bool)

    model_pred = np.zeros(n_obs, dtype=bool)
    for start in range(0, n_obs, FULL_APPLY_CHUNK):
        end = min(start + FULL_APPLY_CHUNK, n_obs)
        rows = np.arange(start, end, dtype=np.int64)
        x = extract_gene_matrix(handle, rows, spec, label=f"full_{start}_{end}")
        score = selected.model_object.predict_proba(x)[:, 1].astype(np.float32)
        if "penalty" in selected.model:
            score = apply_death_penalties(score, x, spec)
        model_pred[start:end] = score >= selected.threshold
        if start and start % 500_000 == 0:
            logging.info("Applied %s model to %s / %s cells", output_prefix, f"{start:,}", f"{n_obs:,}")

    original_threshold = load_original_trd_trab_threshold()
    original_score = read_float_obs(handle, "phase4_trd_minus_trab")
    original_pred = original_score >= original_threshold
    strategy_name = f"{output_prefix}_{selected.model}"
    if exploratory:
        strategy_name += "_UNACCEPTED"
    summaries = [
        summarize_full_predictions(
            strategy=strategy_name,
            threshold=float(selected.threshold),
            pred=model_pred,
            source=source,
            tissue=tissue,
            annotation=annotation,
            sublibrary=sublibrary,
            is_tcrab_seq_library=is_tcrab_seq_library,
            has_TRA_TRB_paired=has_TRA_TRB_paired,
            is_nk=is_nk,
        ),
        summarize_full_predictions(
            strategy="original_TRD_minus_TRAB_score_strategy",
            threshold=float(original_threshold),
            pred=original_pred,
            source=source,
            tissue=tissue,
            annotation=annotation,
            sublibrary=sublibrary,
            is_tcrab_seq_library=is_tcrab_seq_library,
            has_TRA_TRB_paired=has_TRA_TRB_paired,
            is_nk=is_nk,
        ),
    ]
    overall = pd.concat([item[0] for item in summaries], ignore_index=True)
    source_df = pd.concat([item[1] for item in summaries], ignore_index=True)
    tissue_df = pd.concat([item[2] for item in summaries], ignore_index=True)
    fp_source_df = pd.concat([item[3] for item in summaries], ignore_index=True)
    sublibrary_df = pd.concat([item[4] for item in summaries], ignore_index=True)
    annotation_df = pd.concat([item[5] for item in summaries], ignore_index=True)
    fp_overall = overall[
        [
            "strategy",
            "threshold",
            "total_cells",
            "predicted_putative_gdT",
            "predicted_fraction",
            "tcrab_seq_library_cells",
            "predicted_in_tcrab_seq_libraries",
            "paired_TCRAB_cells",
            "predicted_paired_TCRAB_FP",
            "paired_TCRAB_FP_ratio_among_paired_TCRAB",
            "paired_TCRAB_FP_fraction_of_predictions",
            "NK_cells",
            "predicted_NK_FP",
            "NK_FP_ratio_among_NK",
            "NK_FP_fraction_of_predictions",
            "predicted_paired_TCRAB_or_NK_FP",
            "known_FP_fraction_of_predictions",
        ]
    ].copy()
    overall.to_csv(OUT_TABLE_DIR / f"{output_prefix}_full_dataset_prediction_overall.csv", index=False)
    source_df.to_csv(OUT_TABLE_DIR / f"{output_prefix}_full_dataset_prediction_by_source.csv", index=False)
    tissue_df.to_csv(OUT_TABLE_DIR / f"{output_prefix}_full_dataset_prediction_by_tissue.csv", index=False)
    fp_overall.to_csv(OUT_TABLE_DIR / f"{output_prefix}_full_dataset_false_positive_overall.csv", index=False)
    fp_source_df.to_csv(OUT_TABLE_DIR / f"{output_prefix}_full_dataset_false_positive_by_source.csv", index=False)
    sublibrary_df.to_csv(OUT_TABLE_DIR / f"{output_prefix}_full_dataset_false_positive_by_sublibrary.csv", index=False)
    annotation_df.to_csv(OUT_TABLE_DIR / f"{output_prefix}_full_dataset_prediction_by_annotation.csv", index=False)
    return FullApplyOutputs(
        overall=overall,
        source=source_df,
        tissue=tissue_df,
        fp_overall=fp_overall,
        fp_source=fp_source_df,
        fp_sublibrary=sublibrary_df,
        fp_annotation=annotation_df,
    )


def bounded_scatter_limits(x: np.ndarray, y: np.ndarray) -> tuple[tuple[float, float], tuple[float, float]]:
    finite = np.isfinite(x) & np.isfinite(y)
    if not finite.any():
        return (-1.0, 1.0), (-1.0, 1.0)
    x_lo, x_hi = np.nanpercentile(x[finite], [0.2, 99.8])
    y_lo, y_hi = np.nanpercentile(y[finite], [0.2, 99.8])
    pad_x = max((x_hi - x_lo) * 0.08, 0.05)
    pad_y = max((y_hi - y_lo) * 0.08, 0.05)
    return (float(x_lo - pad_x), float(x_hi + pad_x)), (float(y_lo - pad_y), float(y_hi + pad_y))


def plot_category_scatter(
    df: pd.DataFrame,
    *,
    category_col: str,
    category_order: list[str],
    colors: dict[str, str],
    title: str,
    path: Path,
    original_threshold: float,
    seed: int = RANDOM_SEED,
) -> None:
    rng = np.random.default_rng(seed)
    fig, ax = plt.subplots(figsize=(7.2, 6.0), constrained_layout=True)
    x = df["phase4_trab_score"].to_numpy(dtype=float)
    y = df["phase4_trd_score"].to_numpy(dtype=float)
    xlim, ylim = bounded_scatter_limits(x, y)
    for category in category_order:
        idx = np.flatnonzero(df[category_col].astype(str).to_numpy() == category)
        if idx.size == 0:
            continue
        max_points = 140_000 if category == "neither" else 90_000
        if idx.size > max_points:
            idx = rng.choice(idx, size=max_points, replace=False)
        alpha = 0.18 if category == "neither" else 0.72
        size = 2.0 if category == "neither" else 4.0
        ax.scatter(
            x[idx],
            y[idx],
            s=size,
            alpha=alpha,
            linewidths=0,
            rasterized=True,
            color=colors.get(category, "#777777"),
            label=f"{category} (n={int((df[category_col].astype(str) == category).sum()):,})",
        )
    line_x = np.linspace(xlim[0], xlim[1], 100)
    ax.plot(line_x, line_x + original_threshold, color="#111111", lw=1.0, ls="--", label="original TRD-TRAB cutoff")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("Phase 4 TRAB score")
    ax.set_ylabel("Phase 4 TRD score")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=7, loc="best")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_trd_trab_prediction_scatters(
    handle: h5py.File,
    model: Result | None,
    spec: FeatureSpec,
    *,
    sample_cells: int,
    seed: int,
) -> list[Path]:
    if model is None or model.model_object is None:
        return []
    n_obs = int(handle["obs"]["_index"].shape[0])
    trd = read_float_obs(handle, "phase4_trd_score")
    trab = read_float_obs(handle, "phase4_trab_score")
    trd_minus_trab = read_float_obs(handle, "phase4_trd_minus_trab")
    original_threshold = load_original_trd_trab_threshold()
    original_pred_full = trd_minus_trab >= original_threshold
    rng = np.random.default_rng(seed)
    random_idx = rng.choice(n_obs, size=min(sample_cells, n_obs), replace=False)
    original_idx = np.flatnonzero(original_pred_full)
    sample_idx = np.unique(np.concatenate([random_idx, original_idx])).astype(np.int64)
    logging.info("Building TRD-vs-TRAB scatter sample with %s cells", f"{sample_idx.size:,}")
    x_sample = extract_gene_matrix(handle, sample_idx, spec, label="trd_trab_scatter_sample")
    model_score = model.model_object.predict_proba(x_sample)[:, 1].astype(np.float32)
    if "penalty" in model.model:
        model_score = apply_death_penalties(model_score, x_sample, spec)
    model_pred = model_score >= model.threshold
    original_pred = original_pred_full[sample_idx]
    annotation_key = "simple_annotation_plus6" if "simple_annotation_plus6" in handle["obs"] else None
    annotation = clean_group_values(read_obs_column(handle, annotation_key))[sample_idx] if annotation_key else np.full(sample_idx.size, "unknown", dtype=object)
    source = clean_group_values(read_obs_column(handle, "source_gse_id"))[sample_idx]
    has_TRA_TRB_paired = read_bool_obs(handle, "has_TRA_TRB_paired")[sample_idx]
    is_nk = pd.Series(annotation, copy=False).astype(str).str.upper().str.contains("NK", regex=False).to_numpy(dtype=bool)
    agreement = np.full(sample_idx.size, "neither", dtype=object)
    agreement[original_pred & (~model_pred)] = "original score only"
    agreement[model_pred & (~original_pred)] = "TCR-gene only"
    agreement[model_pred & original_pred] = "both methods"
    fp_status = np.full(sample_idx.size, "not TCR-gene predicted", dtype=object)
    fp_status[model_pred] = "TCR-gene predicted other"
    fp_status[model_pred & has_TRA_TRB_paired] = "TCR-gene predicted paired TCRAB"
    fp_status[model_pred & is_nk] = "TCR-gene predicted NK"
    fp_status[model_pred & has_TRA_TRB_paired & is_nk] = "TCR-gene predicted paired TCRAB+NK"
    scatter_df = pd.DataFrame(
        {
            "obs_index": sample_idx,
            "source_gse_id": source,
            "annotation": annotation,
            "phase4_trab_score": trab[sample_idx],
            "phase4_trd_score": trd[sample_idx],
            "phase4_trd_minus_trab": trd_minus_trab[sample_idx],
            "original_trd_minus_trab_predicted": original_pred,
            "tcr_gene_model_score": model_score,
            "tcr_gene_model_predicted": model_pred,
            "method_agreement": agreement,
            "has_TRA_TRB_paired": has_TRA_TRB_paired,
            "is_NK_cell": is_nk,
            "known_fp_status": fp_status,
        }
    )
    scatter_df.to_csv(OUT_TABLE_DIR / "trd_vs_trab_prediction_scatter_sample.csv.gz", index=False)
    method_path = OUT_FIGURE_DIR / "trd_vs_trab_prediction_method_agreement.png"
    fp_path = OUT_FIGURE_DIR / "trd_vs_trab_tcrgene_known_fp_status.png"
    plot_category_scatter(
        scatter_df,
        category_col="method_agreement",
        category_order=["neither", "original score only", "TCR-gene only", "both methods"],
        colors={
            "neither": "#b8c0c8",
            "original score only": "#d97904",
            "TCR-gene only": "#2b7c9b",
            "both methods": "#7a3db8",
        },
        title="TRD vs TRAB score space: method agreement",
        path=method_path,
        original_threshold=original_threshold,
        seed=seed,
    )
    plot_category_scatter(
        scatter_df,
        category_col="known_fp_status",
        category_order=[
            "not TCR-gene predicted",
            "TCR-gene predicted other",
            "TCR-gene predicted paired TCRAB",
            "TCR-gene predicted NK",
            "TCR-gene predicted paired TCRAB+NK",
        ],
        colors={
            "not TCR-gene predicted": "#c4c9ce",
            "TCR-gene predicted other": "#2b7c9b",
            "TCR-gene predicted paired TCRAB": "#b51f2a",
            "TCR-gene predicted NK": "#5a8f28",
            "TCR-gene predicted paired TCRAB+NK": "#111111",
        },
        title="TRD vs TRAB score space: TCR-gene predictions and known FP proxies",
        path=fp_path,
        original_threshold=original_threshold,
        seed=seed,
    )
    return [method_path, fp_path]


def plot_curves(results: list[Result], y_validation: np.ndarray) -> tuple[Path, Path]:
    roc_path = OUT_FIGURE_DIR / "validation_roc.png"
    pr_path = OUT_FIGURE_DIR / "validation_pr.png"
    fig, ax = plt.subplots(figsize=(6.2, 5.2), constrained_layout=True)
    for result in results:
        fpr, tpr, _ = roc_curve(y_validation, result.validation_score)
        ax.plot(fpr, tpr, lw=1.3, label=f"{result.model} ({result.validation_metrics['roc_auc']:.3f})")
    ax.plot([0, 1], [0, 1], ls="--", lw=0.8, color="#777")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("Multi-cohort holdout ROC")
    ax.legend(frameon=False, fontsize=6)
    fig.savefig(roc_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 5.2), constrained_layout=True)
    for result in results:
        precision, recall, _ = precision_recall_curve(y_validation, result.validation_score)
        ax.plot(recall, precision, lw=1.3, label=f"{result.model} ({result.validation_metrics['pr_auc']:.3f})")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Multi-cohort holdout precision-recall")
    ax.legend(frameon=False, fontsize=6)
    fig.savefig(pr_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return roc_path, pr_path


def plot_confusion(selected: Result, y_validation: np.ndarray) -> Path | None:
    if selected is None:
        return None
    path = OUT_FIGURE_DIR / "selected_validation_confusion_matrix.png"
    cm = confusion_matrix(y_validation, selected.validation_pred, labels=[0, 1])
    fig, ax = plt.subplots(figsize=(4.8, 4.2), constrained_layout=True)
    im = ax.imshow(cm, cmap="Blues")
    ax.set_xticks([0, 1], labels=["Pred abT", "Pred gdT"])
    ax.set_yticks([0, 1], labels=["True abT", "True gdT"])
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f"{cm[i, j]:,}", ha="center", va="center", fontsize=13)
    ax.set_title(selected.model)
    fig.colorbar(im, ax=ax, shrink=0.8)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_feature_importance(selected: Result, spec: FeatureSpec) -> Path | None:
    if selected is None or selected.model_object is None:
        return None
    path = OUT_FIGURE_DIR / "selected_feature_importance.png"
    names = np.asarray(spec.gene_names, dtype=object)
    if hasattr(selected.model_object, "feature_importances_"):
        values = np.asarray(selected.model_object.feature_importances_, dtype=float)
        order = np.argsort(values)[::-1][:40]
        plot_values = values[order]
        colors = "#4a6fa5"
        xlabel = "XGBoost importance"
        title = "Selected individual-gene feature importance"
        add_zero = False
    elif hasattr(selected.model_object, "named_steps") and "logisticregression" in selected.model_object.named_steps:
        coef = np.asarray(selected.model_object.named_steps["logisticregression"].coef_[0], dtype=float)
        order = np.argsort(np.abs(coef))[::-1][:40]
        plot_values = coef[order]
        colors = np.where(plot_values >= 0, "#0f766e", "#be123c")
        xlabel = "Standardized logistic coefficient"
        title = "Selected individual-gene logistic coefficients"
        add_zero = True
    else:
        return None
    fig, ax = plt.subplots(figsize=(8, max(5, 0.22 * order.size + 1.5)), constrained_layout=True)
    y = np.arange(order.size)
    ax.barh(y, plot_values, color=colors)
    if add_zero:
        ax.axvline(0, color="#334155", linewidth=0.8)
    ax.set_yticks(y, labels=names[order])
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def copy_assets(paths: list[Path | None]) -> list[str]:
    copied = []
    for path in paths:
        if path is None or not path.exists():
            continue
        target = STATIC_ASSET_DIR / path.name
        shutil.copyfile(path, target)
        copied.append(path.name)
    return copied


def render_pdf(no_pdf: bool) -> None:
    if no_pdf:
        return
    chrome_candidates = [
        Path("/usr/bin/google-chrome"),
        Path("/usr/bin/google-chrome-stable"),
        shutil.which("google-chrome"),
        shutil.which("google-chrome-stable"),
    ]
    chrome = next((Path(path) for path in chrome_candidates if path and Path(path).exists()), None)
    if chrome is None:
        logging.warning("Skipping PDF export because google-chrome was not found")
        return
    subprocess.run(
        [
            str(chrome),
            "--headless",
            "--disable-gpu",
            "--no-sandbox",
            "--print-to-pdf-no-header",
            f"--print-to-pdf={REPORT_PDF}",
            str(REPORT_HTML),
        ],
        check=True,
    )


def compact_numeric(df: pd.DataFrame, digits: int = 4) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_float_dtype(out[col]):
            out[col] = out[col].map(lambda value: "" if pd.isna(value) else round(float(value), digits))
    return out


def pick_columns(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    if df.empty:
        return df
    return df[[col for col in columns if col in df.columns]].copy()


def write_reports(
    input_h5ad: Path,
    obs: ObsMetadata,
    spec: FeatureSpec,
    tune_df: pd.DataFrame,
    validation_df: pd.DataFrame,
    acceptance_df: pd.DataFrame,
    prevalence_df: pd.DataFrame,
    selected: Result | None,
    best_baseline: Result,
    source_metrics: pd.DataFrame,
    tissue_metrics: pd.DataFrame,
    full_overall: pd.DataFrame,
    full_source: pd.DataFrame,
    full_fp_overall: pd.DataFrame,
    full_fp_source: pd.DataFrame,
    full_fp_annotation: pd.DataFrame,
    copied_assets: list[str],
    *,
    holdout_source: str,
    extra_trab_holdout_source: str,
    gdt2020_holdout_tissue: str,
    exclude_suboptimal_source: str,
    model_path: Path | None,
    full_application_note: str,
) -> None:
    selected_name = selected.model if selected is not None else "none"
    selected_prev = prevalence_df.loc[prevalence_df["model"] == selected_name].copy()
    selected_prev_display = selected_prev[
        [
            "prevalence_percent",
            "ppv_observed_at_prevalence",
            "ppv_conservative_wilson_at_prevalence",
            "expected_false_positive_per_million_observed",
            "expected_false_positive_per_million_conservative",
        ]
    ] if not selected_prev.empty else pd.DataFrame()
    split_overall = pd.read_csv(OUT_TABLE_DIR / "split_overall.csv")
    feature_summary = (
        pd.DataFrame({"family": spec.family_by_feature})
        .value_counts("family")
        .rename_axis("feature_family")
        .reset_index(name="n_features")
    )
    feature_examples = (
        pd.DataFrame(
            {
                "feature_family": spec.family_by_feature,
                "gene": spec.gene_names,
                "feature": spec.feature_names,
            }
        )
        .groupby("feature_family", sort=True)
        .head(12)
        .reset_index(drop=True)
    )
    validation_display = compact_numeric(
        pick_columns(
            validation_df,
            [
                "model",
                "threshold",
                "n_cells",
                "n_positive",
                "n_negative",
                "predicted_positive",
                "tp",
                "fp",
                "tn",
                "fn",
                "precision",
                "recall",
                "specificity",
                "f1",
                "f0.5",
                "roc_auc",
                "pr_auc",
            ],
        )
    )
    acceptance_display = compact_numeric(acceptance_df)
    full_overall_display = compact_numeric(
        pick_columns(
            full_overall,
            ["strategy", "threshold", "total_cells", "predicted_putative_gdT", "predicted_fraction"],
        )
    )
    fp_overall_display = compact_numeric(
        pick_columns(
            full_fp_overall,
            [
                "strategy",
                "predicted_putative_gdT",
                "predicted_fraction",
                "predicted_paired_TCRAB_FP",
                "paired_TCRAB_FP_ratio_among_paired_TCRAB",
                "paired_TCRAB_FP_fraction_of_predictions",
                "predicted_NK_FP",
                "NK_FP_ratio_among_NK",
                "NK_FP_fraction_of_predictions",
                "predicted_paired_TCRAB_or_NK_FP",
                "known_FP_fraction_of_predictions",
            ],
        )
    )
    source_display = compact_numeric(
        pick_columns(
            full_source.head(25) if not full_source.empty else full_source,
            ["strategy", "source_gse_id", "total_cells", "predicted_putative_gdT", "predicted_fraction"],
        )
    )
    fp_source_display = compact_numeric(
        pick_columns(
            full_fp_source.head(25) if not full_fp_source.empty else full_fp_source,
            [
                "strategy",
                "source_gse_id",
                "predicted_putative_gdT",
                "predicted_paired_TCRAB_FP",
                "predicted_NK_FP",
                "predicted_paired_TCRAB_or_NK_FP",
                "known_FP_fraction_of_predictions",
            ],
        )
    )
    annotation_display = compact_numeric(
        pick_columns(
            full_fp_annotation,
            [
                "strategy",
                "annotation",
                "total_cells",
                "predicted_putative_gdT",
                "predicted_fraction",
                "predicted_paired_TCRAB_FP",
                "predicted_NK_FP",
                "known_FP_fraction_of_predictions",
            ],
        )
    )
    learned_validation = validation_df.loc[~validation_df["model"].astype(str).str.startswith("baseline_")].copy()
    best_learned = learned_validation.iloc[0] if not learned_validation.empty else None
    overfit_lines = [
        "The current evidence reduces, but does not eliminate, overfitting risk.",
        f"All `{holdout_source}` primary gold cells, `{GDT2020_SOURCE}` `{gdt2020_holdout_tissue}` gdT cells, and paired-TCRAB/no-gdTCR cells from `{extra_trab_holdout_source}` were held out from both training and threshold tuning.",
        "Phase 4 TRD/TRAB score columns were not used as model features.",
        "Remaining risks are source-specific TCR-gene capture, differences in TCR library chemistry, and labels that are themselves derived from TCR metadata.",
    ]
    if best_learned is not None:
        overfit_lines.append(
            f"The best learned model on the strict holdout had precision `{float(best_learned['precision']):.3f}`, "
            f"recall `{float(best_learned['recall']):.3f}`, and specificity `{float(best_learned['specificity']):.3f}` "
            f"({int(best_learned['fp']):,} FP among {int(best_learned['n_negative']):,} combined holdout negative cells)."
        )
    if not full_fp_overall.empty:
        first_fp = full_fp_overall.iloc[0]
        overfit_lines.append(
            f"On the whole-dataset proxy audit, paired-TCRAB FP was `{int(first_fp['predicted_paired_TCRAB_FP']):,}` "
            f"and NK FP was `{int(first_fp['predicted_NK_FP']):,}`; their union was "
            f"`{float(first_fp['known_FP_fraction_of_predictions']):.2%}` of predictions."
        )
    overfit_lines.append(
        "Interpretation: this is strong enough to justify targeted external validation, but not strong enough to claim overfitting is ruled out."
    )
    selected_block = []
    if selected is not None:
        selected_block = [
            f"- Selected model: `{selected.model}`",
            f"- Threshold: `{selected.threshold}`",
            f"- Holdout precision: `{selected.validation_metrics['precision']:.4f}`",
            f"- Holdout recall: `{selected.validation_metrics['recall']:.4f}`",
            f"- Holdout specificity: `{selected.validation_metrics['specificity']:.4f}`",
            f"- Holdout F1: `{selected.validation_metrics['f1']:.4f}`",
            f"- Holdout F0.5: `{selected.validation_metrics['f0.5']:.4f}`",
            f"- Holdout ROC-AUC: `{selected.validation_metrics['roc_auc']:.4f}`",
            f"- Holdout PR-AUC: `{selected.validation_metrics['pr_auc']:.4f}`",
            f"- GSE144469 silver recall: `{selected.silver_recall}`",
            f"- Model artifact: `{model_path}`",
        ]
    else:
        selected_block = [
            "- No model passed the acceptance rule on the GSE144469 holdout.",
            "- No production model artifact was published.",
        ]
    full_block = []
    if not full_overall.empty:
        full_block = [
            f"- {full_application_note}",
            "- Full-dataset application includes the strict TCR-gene classifier and the original `TRD_minus_TRAB` score strategy for comparison.",
            "- Paired `TRA/TRB` cells predicted as gdT are counted as false positives in TCRAB-sequenced libraries.",
            "- `simple_annotation_plus6 == NK_cell` cells predicted as gdT are counted as false positives.",
        ]
    else:
        full_block = ["- Full-dataset application was skipped because no model passed or `--no-full-apply` was used."]

    lines = [
        "# Multi-Cohort Holdout Individual TCR-Gene gdT Classifier",
        "",
        "## Scope",
        "",
        f"- Input H5AD: `{input_h5ad}`",
        f"- Primary holdout validation source: `{holdout_source}`. All primary gold cells from this source are excluded from training and threshold tuning.",
        f"- Additional abT validation source: `{extra_trab_holdout_source}` paired-TCRAB/no-gdTCR cells.",
        f"- Additional gdT validation cohort: `{GDT2020_SOURCE}` cells with `tissue_corrected == {gdt2020_holdout_tissue!r}`.",
        f"- Suboptimal sorted source excluded from train/tune: `{exclude_suboptimal_source}`.",
        "- Phase 4 `phase4_trd_score`, `phase4_trab_score`, and `phase4_trd_minus_trab` columns are not read and not used.",
        "- TCR metadata is used only to define gold/silver labels. It is not used as model features.",
        "- Feature matrix uses individual TCR alpha/beta/gamma/delta genes plus FOXP3, CD4, and CD3D/CD3E/CD3G penalty-control genes.",
        f"- Feature gene cap: `{MAX_FEATURE_GENES}`; selected feature genes: `{len(spec.gene_names)}`.",
        "",
        "## Split Summary",
        "",
        dataframe_to_markdown(split_overall),
        "",
        "## Training and Test Protocol",
        "",
        "- Primary labels were built from project TCR metadata. `gdT_gold` and `abT_gold` were the only classes used for model fitting and primary validation; `gdT_silver` was sensitivity-only.",
        f"- Validation was multi-cohort-held-out: every primary gold cell from `{holdout_source}`, every `{GDT2020_SOURCE}` `{gdt2020_holdout_tissue}` gdT cell, and paired-TCRAB/no-gdTCR cells from `{extra_trab_holdout_source}` were excluded from training and threshold tuning.",
        f"- `{exclude_suboptimal_source}` was excluded from train/tune because its library quality was flagged as suboptimal.",
        "- The remaining primary gold cells were split 80/20 into train and tune within each source and label stratum.",
        f"- Training used all positive train cells and at most `{MAX_NEGATIVE_TRAIN:,}` randomly sampled negative train cells to control the severe class imbalance.",
        "- Candidate models were two individual-gene heuristic baselines, XGBoost, XGBoost plus FOXP3/CD4/low-CD3 death penalties, and balanced logistic regression.",
        "- Thresholds were selected only on the tune split by maximizing F1, without the previous high-specificity constraint.",
        "- The GSE144469 holdout was used once for final validation and acceptance.",
        "- A learned model was accepted if holdout F1 exceeded the best individual-gene heuristic baseline under the same F1 threshold-selection rule.",
        "- The whole-dataset application uses the selected F1-optimized model and compares it to the original `TRD_minus_TRAB` score strategy.",
        "",
        "## Feature Summary",
        "",
        "Feature selection was deterministic and did not use model performance:",
        "",
        "- Candidate genes were read from `var/_index` in the input H5AD.",
        "- TCR genes were selected by gene-symbol prefixes: `TRAV`, `TRAJ`, `TRAC`, `TRBV`, `TRBJ`, `TRBC`, `TRGV`, `TRGJ`, `TRGC`, `TRDV`, `TRDD`, `TRDJ`, and `TRDC`.",
        "- Genes were ordered before applying the cap as: delta (`TRDC/TRDV/TRDD/TRDJ`), gamma (`TRGC/TRGV/TRGJ`), alpha (`TRAC/TRAV/TRAJ`), then beta (`TRBC/TRBV/TRBJ`), each alphabetically within priority group.",
        "- `FOXP3`, `CD4`, `CD3D`, `CD3E`, and `CD3G` were appended when present as penalty/control genes, not as TCR-score features.",
        f"- The cap was `{MAX_FEATURE_GENES}` genes. This object had `{len(spec.tcr_gene_names)}` TCR genes plus `{len(spec.penalty_gene_names)}` control genes, so no selected feature was dropped by the cap.",
        "- Every feature value is per-cell `log1p(counts per 10,000)` computed directly from the sparse `X` matrix.",
        "",
        dataframe_to_markdown(feature_summary),
        "",
        "Feature examples from each family:",
        "",
        dataframe_to_markdown(feature_examples, max_rows=80),
        "",
        "## Selected Model",
        "",
        *selected_block,
        f"- Best individual-gene heuristic baseline: `{best_baseline.model}` with holdout F1 `{best_baseline.validation_metrics['f1']:.4f}`.",
        "",
        "## Holdout Validation Metrics",
        "",
        dataframe_to_markdown(validation_display, max_rows=20),
        "",
        "## Acceptance",
        "",
        dataframe_to_markdown(acceptance_display, max_rows=20),
        "",
        "## Overfitting Risk Interpretation",
        "",
        *[f"- {line}" for line in overfit_lines],
        "",
        "## Prevalence-Aware PPV",
        "",
        dataframe_to_markdown(selected_prev_display, max_rows=20),
        "",
        "## Validation Strata",
        "",
        "By source:",
        "",
        dataframe_to_markdown(source_metrics, max_rows=20),
        "",
        "By tissue:",
        "",
        dataframe_to_markdown(tissue_metrics, max_rows=30),
        "",
        "## Full Dataset Application",
        "",
        *full_block,
        "",
        dataframe_to_markdown(full_overall_display, max_rows=20),
        "",
        "Known false-positive audit:",
        "",
        dataframe_to_markdown(fp_overall_display, max_rows=20),
        "",
        "Top full-dataset predicted sources:",
        "",
        dataframe_to_markdown(source_display, max_rows=30),
        "",
        "False-positive audit by source:",
        "",
        dataframe_to_markdown(fp_source_display, max_rows=30),
        "",
        "Prediction by annotation:",
        "",
        dataframe_to_markdown(annotation_display, max_rows=40),
        "",
        "## TRD-vs-TRAB Score Scatter",
        "",
        "- Scatter plots use all original `TRD_minus_TRAB` positive cells plus a reproducible random background sample.",
        "- The classifier predictions in these figures are computed from individual TCR-gene expression only; TRD/TRAB scores are plotted for interpretation, not used as classifier features.",
        "- Dashed line marks the original conservative `TRD_minus_TRAB` cutoff.",
        "",
        "## Outputs",
        "",
        f"- Validation metrics: `{OUT_TABLE_DIR / 'validation_metrics.csv'}`",
        f"- Acceptance table: `{OUT_TABLE_DIR / 'acceptance_vs_individual_gene_baseline.csv'}`",
        f"- Feature manifest: `{OUT_TABLE_DIR / 'feature_manifest.csv'}`",
        f"- Prevalence table: `{OUT_TABLE_DIR / 'prevalence_aware_ppv_scenarios.csv'}`",
        f"- Full application comparison: `{OUT_TABLE_DIR / 'selected_model_full_dataset_prediction_overall.csv'}`",
        f"- Full false-positive audit: `{OUT_TABLE_DIR / 'selected_model_full_dataset_false_positive_overall.csv'}`",
        f"- TRD-vs-TRAB scatter sample: `{OUT_TABLE_DIR / 'trd_vs_trab_prediction_scatter_sample.csv.gz'}`",
        f"- TRD-vs-TRAB method scatter: `{OUT_FIGURE_DIR / 'trd_vs_trab_prediction_method_agreement.png'}`",
        f"- TRD-vs-TRAB FP scatter: `{OUT_FIGURE_DIR / 'trd_vs_trab_tcrgene_known_fp_status.png'}`",
        f"- HTML report: `{REPORT_HTML}`",
        f"- PDF report: `{REPORT_PDF}`",
        "",
    ]
    REPORT_MD.write_text("\n".join(lines), encoding="utf-8")

    figure_titles = {
        "validation_roc.png": "Multi-cohort holdout ROC",
        "validation_pr.png": "Multi-cohort holdout precision-recall",
        "selected_validation_confusion_matrix.png": "Selected confusion matrix",
        "selected_feature_importance.png": "Selected feature importance",
        "trd_vs_trab_prediction_method_agreement.png": "TRD-vs-TRAB score space: method agreement",
        "trd_vs_trab_tcrgene_known_fp_status.png": "TRD-vs-TRAB score space: known FP proxies",
    }
    figure_html = "\n".join(
        f"<section class='figure'><h3>{html.escape(figure_titles.get(name, name))}</h3>"
        f"<img src='assets/{html.escape(name)}' alt='{html.escape(figure_titles.get(name, name))}'></section>"
        for name in copied_assets
    )
    css = """
    :root{--bg:#f5f6f7;--paper:#fff;--ink:#1d252c;--muted:#63707c;--line:#d9dee4;--accent:#8f2d2d}
    @page{size:A4 landscape;margin:8mm}
    *{box-sizing:border-box} body{margin:0;background:var(--bg);color:var(--ink);font-family:Arial,Helvetica,sans-serif;line-height:1.5}
    main{width:min(1320px,calc(100vw - 28px));margin:18px auto 38px} section{background:var(--paper);border:1px solid var(--line);padding:18px;margin-top:14px;break-inside:avoid}.hero{border-top:6px solid var(--accent)}
    h1{font-size:30px;line-height:1.15;margin:0 0 10px;letter-spacing:0} h2{font-size:21px;margin:0 0 12px} h3{font-size:15px;margin:0 0 8px} p{margin:0 0 12px;color:var(--muted)} code{background:#eef1f4;padding:1px 5px;border-radius:4px}
    .table-wrap{width:100%;overflow-x:auto;margin:8px 0 14px}
    table{border-collapse:collapse;width:100%;font-size:9px;line-height:1.25;table-layout:fixed} th,td{border:1px solid var(--line);padding:3px 4px;text-align:left;vertical-align:top;overflow-wrap:anywhere} th{background:#eef1f4}
    .figures{display:grid;grid-template-columns:repeat(auto-fit,minmax(480px,1fr));gap:14px}.figure{padding:10px;background:#fbfcfd}.figure img{width:100%;height:auto;border:1px solid var(--line);background:white}
    ul{margin:6px 0 0 18px;padding:0}.note{color:var(--muted)}
    """
    method_steps = [
        f"Hold out every `{holdout_source}` primary gold cell for validation.",
        f"Hold out `{GDT2020_SOURCE}` `{gdt2020_holdout_tissue}` gdT cells as an independent sorted gdT-positive cohort.",
        f"Hold out `{extra_trab_holdout_source}` paired-TCRAB/no-gdTCR cells as an independent abT false-positive challenge cohort.",
        f"Exclude `{exclude_suboptimal_source}` from train/tune as suboptimal-quality sensitivity source.",
        "Split remaining primary gold cells 80/20 into train and tune within each source and label.",
        f"Train with all positives and up to {MAX_NEGATIVE_TRAIN:,} negatives.",
        "Select thresholds on tune by maximizing F1, without the previous high-specificity constraint.",
        "Accept a learned model if holdout F1 beats the best individual-gene heuristic baseline.",
    ]
    method_html = "<ul>" + "".join(f"<li>{html.escape(step)}</li>" for step in method_steps) + "</ul>"
    overfit_html = "<ul>" + "".join(f"<li>{html.escape(line)}</li>" for line in overfit_lines) + "</ul>"
    html_parts = [
        "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'>",
        "<title>Multi-cohort holdout TCR-gene gdT classifier</title>",
        f"<style>{css}</style></head><body><main>",
        "<section class='hero'><h1>Multi-Cohort Holdout TCR-Gene gdT Classifier</h1>",
        f"<p><code>{html.escape(holdout_source)}</code>, <code>{html.escape(extra_trab_holdout_source)}</code> paired-TCRAB cells, and <code>{html.escape(GDT2020_SOURCE)}</code> {html.escape(gdt2020_holdout_tissue)} gdT cells are held out. Phase 4 TRD/TRAB scores are not used as model features.</p></section>",
        f"<section><h2>Scope</h2><div class='table-wrap'>{dataframe_to_html(split_overall)}</div></section>",
        f"<section><h2>Training and Test Protocol</h2>{method_html}<p class='note'>The selected operating point maximizes F1 on the tune split, so recall is allowed to increase when it improves the precision-recall balance.</p></section>",
        "<section><h2>Feature Selection</h2>"
        "<p>Features were selected deterministically from <code>var/_index</code>: all individual TCR genes matching "
        "<code>TRAV/TRAJ/TRAC</code>, <code>TRBV/TRBJ/TRBC</code>, <code>TRGV/TRGJ/TRGC</code>, and "
        "<code>TRDV/TRDD/TRDJ/TRDC</code>, ordered as delta, gamma, alpha, beta before the 300-gene cap. "
        "<code>FOXP3</code>, <code>CD4</code>, and <code>CD3D/CD3E/CD3G</code> were appended as penalty/control genes. "
        "This run used 182 TCR genes plus 5 controls; no feature was dropped by the cap. Values are per-cell "
        "<code>log1p(counts per 10,000)</code> from sparse <code>X</code>.</p>"
        f"<div class='table-wrap'>{dataframe_to_html(feature_summary)}</div><div class='table-wrap'>{dataframe_to_html(feature_examples, max_rows=80)}</div></section>",
        f"<section><h2>Holdout Validation</h2><p>Accepted production model: <code>{html.escape(selected_name)}</code></p><div class='table-wrap'>{dataframe_to_html(validation_display, max_rows=20)}</div></section>",
        f"<section><h2>Acceptance and Overfitting Risk</h2><div class='table-wrap'>{dataframe_to_html(acceptance_display, max_rows=20)}</div>{overfit_html}</section>",
        f"<section><h2>Figures</h2><div class='figures'>{figure_html}</div></section>",
        f"<section><h2>Prevalence-Aware PPV</h2><div class='table-wrap'>{dataframe_to_html(selected_prev_display, max_rows=20)}</div></section>",
        f"<section><h2>Full Dataset Application</h2><p>{html.escape(full_application_note)}</p><div class='table-wrap'>{dataframe_to_html(full_overall_display if not full_overall_display.empty else pd.DataFrame({'status': [full_block[0]]}))}</div></section>",
        f"<section><h2>False-Positive Audit</h2><div class='table-wrap'>{dataframe_to_html(fp_overall_display, max_rows=20)}</div><div class='table-wrap'>{dataframe_to_html(fp_source_display, max_rows=30)}</div></section>",
        f"<section><h2>Prediction by Source and Annotation</h2><div class='table-wrap'>{dataframe_to_html(source_display, max_rows=30)}</div><div class='table-wrap'>{dataframe_to_html(annotation_display, max_rows=40)}</div></section>",
        "</main></body></html>",
    ]
    REPORT_HTML.write_text("\n".join(html_parts), encoding="utf-8")


def main() -> None:
    args = parse_args()
    setup_logging()
    input_h5ad = args.input_h5ad.resolve()
    stat_before = input_h5ad.stat()
    with h5py.File(input_h5ad, "r") as handle:
        n_obs = int(handle["obs"]["_index"].shape[0])
        logging.info("Reading labels for %s cells", f"{n_obs:,}")
        obs = build_obs_metadata(handle)
        spec = build_feature_spec(handle, args.max_feature_genes)
        train_idx, tune_idx, validation_idx, excluded_idx, _split_summary = make_train_tune_validation_splits(
            obs,
            holdout_source=args.holdout_source,
            extra_trab_holdout_source=args.extra_trab_holdout_source,
            gdt2020_holdout_tissue=args.gdt2020_holdout_tissue,
            exclude_suboptimal_source=args.exclude_suboptimal_source,
            seed=args.seed,
        )
        silver_idx = np.flatnonzero(obs.silver_mask & (obs.source == args.holdout_source)).astype(np.int64)
        eval_rows = np.concatenate([train_idx, tune_idx, validation_idx, silver_idx, excluded_idx])
        eval_rows = np.unique(eval_rows)
        y_eval = (obs.class_code[eval_rows] == 2).astype(np.int8)
        pos_train = local_positions(eval_rows, train_idx)
        pos_tune = local_positions(eval_rows, tune_idx)
        pos_validation = local_positions(eval_rows, validation_idx)
        pos_silver = local_positions(eval_rows, silver_idx) if silver_idx.size else np.array([], dtype=np.int64)
        logging.info("Extracting individual TCR-gene feature matrix for %s rows", f"{eval_rows.size:,}")
        x_eval = extract_gene_matrix(handle, eval_rows, spec, label="train_tune_holdout")
        results = train_models(
            x_eval,
            y_eval,
            pos_train,
            pos_tune,
            pos_validation,
            pos_silver,
            spec,
            max_negative_train=args.max_negative_train,
            min_specificity=args.min_specificity,
            seed=args.seed,
        )
        tune_df = results_to_frame(results, "tune")
        validation_df = results_to_frame(results, "validation")
        tune_df.to_csv(OUT_TABLE_DIR / "tune_metrics.csv", index=False)
        validation_df.to_csv(OUT_TABLE_DIR / "validation_metrics.csv", index=False)
        selected, best_baseline = select_model(results)
        exploratory_model = choose_exploratory_model(results)
        acceptance_df = pd.read_csv(OUT_TABLE_DIR / "acceptance_vs_individual_gene_baseline.csv")
        prevalence_df = prevalence_table(results)
        if selected is not None:
            source_metrics, tissue_metrics = write_validation_strata(
                selected,
                validation_idx,
                (obs.class_code[validation_idx] == 2).astype(np.int8),
                obs,
            )
            with MODEL_PKL.open("wb") as out:
                pickle.dump(
                    {
                        "model": selected.model,
                        "threshold": selected.threshold,
                        "feature_names": spec.feature_names,
                        "gene_names": spec.gene_names,
                        "model_object": selected.model_object,
                        "notes": selected.notes,
                    },
                    out,
                )
            model_path: Path | None = MODEL_PKL
        else:
            source_metrics = pd.DataFrame()
            tissue_metrics = pd.DataFrame()
            model_path = None
        empty = pd.DataFrame()
        full_outputs = FullApplyOutputs(empty, empty, empty, empty, empty, empty, empty)
        full_application_note = "Full-dataset application was skipped."
        applied_full_model: Result | None = selected
        output_prefix = "selected_model"
        exploratory = False
        scatter_paths: list[Path] = []
        if selected is None and args.apply_best_unaccepted_full:
            applied_full_model = exploratory_model
            output_prefix = "exploratory_best_unaccepted"
            exploratory = True
        if applied_full_model is not None and not args.no_full_apply:
            if exploratory:
                full_application_note = (
                    f"Exploratory application used the best learned but unaccepted holdout model "
                    f"`{applied_full_model.model}`. This is not a production classifier."
                )
            else:
                full_application_note = f"Accepted holdout model `{applied_full_model.model}` was applied to the full 5-million-cell atlas."
            logging.info("Applying %s to full 5-million-cell atlas", applied_full_model.model)
            full_outputs = apply_full_dataset(
                handle,
                applied_full_model,
                spec,
                output_prefix=output_prefix,
                exploratory=exploratory,
            )
        elif applied_full_model is not None and args.reuse_existing_full_outputs:
            full_outputs = read_existing_full_outputs(output_prefix)
            if not full_outputs.overall.empty:
                if exploratory:
                    full_application_note = (
                        f"Reused existing exploratory full-dataset output for the best learned but unaccepted model "
                        f"`{applied_full_model.model}`. This is not a production classifier."
                    )
                else:
                    full_application_note = f"Reused existing full-dataset output for `{applied_full_model.model}`."
        if applied_full_model is not None:
            scatter_paths = plot_trd_trab_prediction_scatters(
                handle,
                applied_full_model,
                spec,
                sample_cells=args.scatter_sample_cells,
                seed=args.seed,
            )

    stat_after = input_h5ad.stat()
    if (stat_before.st_size != stat_after.st_size) or (stat_before.st_mtime_ns != stat_after.st_mtime_ns):
        raise RuntimeError("Input H5AD changed during read-only holdout classifier analysis.")

    y_validation = (obs.class_code[validation_idx] == 2).astype(np.int8)
    roc_path, pr_path = plot_curves(results, y_validation)
    confusion_path = plot_confusion(selected, y_validation)
    importance_path = plot_feature_importance(selected, spec)
    copied_assets = copy_assets([roc_path, pr_path, confusion_path, importance_path, *scatter_paths])
    write_reports(
        input_h5ad,
        obs,
        spec,
        tune_df,
        validation_df,
        acceptance_df,
        prevalence_df,
        selected,
        best_baseline,
        source_metrics,
        tissue_metrics,
        full_outputs.overall,
        full_outputs.source,
        full_outputs.fp_overall,
        full_outputs.fp_source,
        full_outputs.fp_annotation,
        copied_assets,
        holdout_source=args.holdout_source,
        extra_trab_holdout_source=args.extra_trab_holdout_source,
        gdt2020_holdout_tissue=args.gdt2020_holdout_tissue,
        exclude_suboptimal_source=args.exclude_suboptimal_source,
        model_path=model_path,
        full_application_note=full_application_note,
    )
    render_pdf(args.no_pdf)
    SUMMARY_JSON.write_text(
        json.dumps(
            json_ready(
                {
                    "input_h5ad": str(input_h5ad),
                    "holdout_source": args.holdout_source,
                    "extra_trab_holdout_source": args.extra_trab_holdout_source,
                    "gdt2020_holdout_tissue": args.gdt2020_holdout_tissue,
                    "selected_model": None if selected is None else selected.model,
                    "exploratory_full_model": None if applied_full_model is None else applied_full_model.model,
                    "exploratory_full_application": bool(exploratory and applied_full_model is not None and not args.no_full_apply),
                    "report_md": str(REPORT_MD),
                    "report_html": str(REPORT_HTML),
                    "report_pdf": str(REPORT_PDF),
                    "model_path": None if model_path is None else str(model_path),
                }
            ),
            indent=2,
        ),
        encoding="utf-8",
    )
    logging.info("Saved holdout report: %s", REPORT_MD)
    logging.info("Saved holdout HTML: %s", REPORT_HTML)
    if REPORT_PDF.exists():
        logging.info("Saved holdout PDF: %s", REPORT_PDF)


if __name__ == "__main__":
    main()
