#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Train and evaluate gdTAI v4 two-stage T-lineage/gdT classifier.

Design intent:
- first learn a soft, high-recall T-lineage gate;
- then learn gdT probability among T-like / T-NK-like cells;
- avoid hard TRDV exclusion so TRDC+TRDV- dropout-like gdT cells can be recovered;
- include high-confidence NK cells as negatives without treating cytotoxic genes
  as single-marker death penalties.

The input H5AD is opened read-only. No H5AD is mutated.
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
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)

from run_gdt_gse144469_holdout_tcrgene_classifier import (
    GDT2020_SOURCE,
    GDT2020_HOLDOUT_TISSUE,
    RANDOM_SEED,
    SUBOPTIMAL_SOURCE,
    build_obs_metadata,
    tcr_priority,
)
from run_gdt_prediction_package_evaluation import (
    DEFAULT_INPUT_H5AD,
    clean_group_values,
    dataframe_to_html,
    dataframe_to_markdown,
    read_obs_column,
)
from run_gdtai_v3_trdc_nk_guard_classifier import FeatureSpec, extract_gene_features, obs_index_values

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
STATIC_DIR = PROJECT_ROOT / "gdT_prediction"
OUT_PREFIX = "gdtai_v4_tcell_gate"
TABLE_DIR = OUTPUT_ROOT / "tables" / "gdT_prediction" / OUT_PREFIX
FIGURE_DIR = OUTPUT_ROOT / "figures" / "gdT_prediction" / OUT_PREFIX
LOG_DIR = OUTPUT_ROOT / "logs" / "gdT_prediction" / OUT_PREFIX
MODEL_DIR = OUTPUT_ROOT / "models" / "gdT_prediction_classifier" / "gdTAI_v4.0"
STATIC_ASSET_DIR = STATIC_DIR / "assets" / OUT_PREFIX
REPORT_MD = LOG_DIR / "gdtai_v4_tcell_gate_report.md"
REPORT_HTML = STATIC_DIR / "gdtai_v4_tcell_gate_report.html"
REPORT_PDF = STATIC_DIR / "gdtai_v4_tcell_gate_report.pdf"
RUN_LOG = LOG_DIR / "run.log"
SUMMARY_JSON = LOG_DIR / "gdtai_v4_tcell_gate_summary.json"
MODEL_PKL = MODEL_DIR / "gdTAI_v4_model.pkl"
MANIFEST_JSON = MODEL_DIR / "model_manifest.json"

TARGET_SUM = 10_000.0
MAX_FEATURE_GENES = 300
FULL_CHUNK = 50_000

TCR_PREFIXES = ("TRAV", "TRAJ", "TRAC", "TRBV", "TRBJ", "TRBC", "TRGV", "TRGJ", "TRGC", "TRDV", "TRDD", "TRDJ", "TRDC")
T_ANNOTATIONS = {"CD4_T", "CD8_T", "TREG", "GDT_CELL"}
NON_T_ANNOTATIONS = {"NK_CELL", "B_CELL", "MYELOID", "OTHER"}

T_LINEAGE_GENES = ["CD3D", "CD3E", "CD3G", "CD2", "CD247", "LCK", "LAT", "TRAT1", "CD5", "CD6", "THEMIS"]
NK_CONTEXT_GENES = ["NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "KLRD1", "KLRF1", "FCGR3A", "NCAM1", "TYROBP", "FCER1G", "CST7", "CTSW", "EOMES", "TBX21"]
CONTROL_GENES = ["CD4", "FOXP3", "IL7R", "CCR7", "SELL", "CD8A", "CD8B", "MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "LST1", "LYZ", "S100A8", "S100A9", "FCN1", "CST3", "MKI67", "TOP2A"]

ENGINEERED_FEATURES = [
    "any_TRDV",
    "any_TRDJ_TRDD",
    "any_TRG",
    "any_ab_TCR_gene",
    "TRDC_only",
    "TRDC_plus_TRDV",
    "TRDC_plus_TRG_no_TRDV",
    "CD3_score",
    "T_lineage_score",
    "NK_context_score",
    "cytotoxic_score",
    "gdT_TCR_score",
    "abT_TCR_score",
    "gdT_minus_abT_TCR_score",
    "CD3_minus_NK_score",
    "TRDC_log1p",
    "TRDV_score",
    "TRDJ_TRDD_score",
    "TRG_score",
    "cytotoxic_with_TCR_support",
    "cytotoxic_without_TCR_support",
    "CD4_FOXP3_score",
    "B_score",
    "myeloid_score",
]


@dataclass
class SplitRows:
    stage1_train: np.ndarray
    stage1_tune: np.ndarray
    stage2_train: np.ndarray
    stage2_tune: np.ndarray
    validation: np.ndarray
    validation_gdt: np.ndarray
    validation_abt: np.ndarray
    sensitivity_gdtlung: np.ndarray
    sensitivity_silver: np.ndarray
    all_eval_rows: np.ndarray
    split_summary: pd.DataFrame


class GdTAIV4TwoStageModel:
    """Pickle-safe two-stage gdTAI v4 inference wrapper."""

    def __init__(
        self,
        *,
        tcell_model: Any,
        gdt_model: Any,
        tcell_threshold: float,
        mode_thresholds: dict[str, float],
        gene_names: list[str],
        engineered_feature_names: list[str],
        model_feature_names: list[str],
    ) -> None:
        self.tcell_model = tcell_model
        self.gdt_model = gdt_model
        self.tcell_threshold = float(tcell_threshold)
        self.mode_thresholds = {str(k): float(v) for k, v in mode_thresholds.items()}
        self.gene_names = list(gene_names)
        self.engineered_feature_names = list(engineered_feature_names)
        self.model_feature_names = list(model_feature_names)

    def score_arrays(self, x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        t_score = self.tcell_model.predict_proba(x)[:, 1].astype(np.float32)
        x2 = np.column_stack([x, t_score]).astype(np.float32, copy=False)
        gdt_score = self.gdt_model.predict_proba(x2)[:, 1].astype(np.float32)
        joint_score = (t_score * gdt_score).astype(np.float32)
        return t_score, gdt_score, joint_score

    def predict_mode(self, x: np.ndarray, mode: str = "high_f1") -> np.ndarray:
        t_score, gdt_score, _joint_score = self.score_arrays(x)
        return (t_score >= self.tcell_threshold) & (gdt_score >= self.mode_thresholds[mode])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train gdTAI v4 two-stage T-cell gate/gdT classifier.")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    parser.add_argument("--chunk-size", type=int, default=FULL_CHUNK)
    parser.add_argument("--stage1-max-positive", type=int, default=320_000)
    parser.add_argument("--stage1-max-nk-negative", type=int, default=180_000)
    parser.add_argument("--stage2-max-ab-negative", type=int, default=300_000)
    parser.add_argument("--stage2-max-nk-negative", type=int, default=100_000)
    parser.add_argument("--target-tcell-gdt-recall", type=float, default=0.995)
    parser.add_argument("--target-tcell-positive-recall", type=float, default=0.985)
    parser.add_argument("--max-full-atlas-cells", type=int, default=None)
    parser.add_argument("--no-full-atlas", action="store_true")
    parser.add_argument("--no-pdf", action="store_true")
    return parser.parse_args()


def ensure_dirs() -> None:
    for path in [TABLE_DIR, FIGURE_DIR, LOG_DIR, MODEL_DIR, STATIC_DIR, STATIC_ASSET_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def setup_logging() -> None:
    ensure_dirs()
    handlers = [logging.FileHandler(RUN_LOG, mode="w", encoding="utf-8"), logging.StreamHandler()]
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s", handlers=handlers, force=True)


def h5ad_shape(handle: h5py.File) -> tuple[int, int]:
    x = handle["X"]
    if isinstance(x, h5py.Group) and "indptr" in x:
        n_obs = int(x["indptr"].shape[0] - 1)
        var = handle["var"]
        key = "_index" if "_index" in var else "index"
        return n_obs, int(var[key].shape[0])
    return int(x.shape[0]), int(x.shape[1])


def read_var_names(handle: h5py.File) -> list[str]:
    var = handle["var"]
    key = "_index" if "_index" in var else "index"
    ds = var[key]
    try:
        values = ds.asstr()[:]
    except Exception:
        values = ds[:]
    return [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in values]


def annotation_vector(handle: h5py.File, n_obs: int) -> np.ndarray:
    if "simple_annotation_plus6" in handle["obs"]:
        return clean_group_values(read_obs_column(handle, "simple_annotation_plus6"))
    if "simple_annotation" in handle["obs"]:
        return clean_group_values(read_obs_column(handle, "simple_annotation"))
    return np.full(n_obs, "unknown", dtype=object)


def make_feature_spec(handle: h5py.File) -> FeatureSpec:
    var_names = read_var_names(handle)
    lookup = {gene: i for i, gene in enumerate(var_names)}
    tcr_genes = sorted([g for g in var_names if g.upper().startswith(TCR_PREFIXES)], key=tcr_priority)
    selected: list[str] = []
    for group in [tcr_genes, T_LINEAGE_GENES, NK_CONTEXT_GENES, CONTROL_GENES]:
        for gene in group:
            if gene in lookup and gene not in selected:
                selected.append(gene)
    if len(selected) > MAX_FEATURE_GENES:
        selected = selected[:MAX_FEATURE_GENES]
    spec = FeatureSpec(
        gene_names=selected,
        gene_indices=np.asarray([lookup[g] for g in selected], dtype=np.int32),
        gene_feature_names=[f"{g}_log1p_cp10k" for g in selected],
        engineered_feature_names=list(ENGINEERED_FEATURES),
        model_feature_names=[*[f"{g}_log1p_cp10k" for g in selected], *ENGINEERED_FEATURES],
        gene_to_col={g: i for i, g in enumerate(selected)},
        engineered_to_col={name: len(selected) + i for i, name in enumerate(ENGINEERED_FEATURES)},
    )
    pd.DataFrame(
        {
            "feature": spec.model_feature_names,
            "feature_index": np.arange(len(spec.model_feature_names), dtype=int),
            "feature_type": ["gene_log1p_cp10k"] * len(spec.gene_names) + ["engineered"] * len(ENGINEERED_FEATURES),
            "gene": spec.gene_names + [""] * len(ENGINEERED_FEATURES),
        }
    ).to_csv(TABLE_DIR / "feature_manifest.csv", index=False)
    return spec


def gene_cols(spec: FeatureSpec, names: list[str] | tuple[str, ...]) -> list[int]:
    return [spec.gene_to_col[name] for name in names if name in spec.gene_to_col]


def prefix_cols(spec: FeatureSpec, prefixes: tuple[str, ...]) -> list[int]:
    return [i for gene, i in spec.gene_to_col.items() if gene.upper().startswith(prefixes)]


def mean_cols(x: np.ndarray, cols: list[int]) -> np.ndarray:
    if not cols:
        return np.zeros(x.shape[0], dtype=np.float32)
    return x[:, cols].mean(axis=1).astype(np.float32)


def sum_cols(x: np.ndarray, cols: list[int]) -> np.ndarray:
    if not cols:
        return np.zeros(x.shape[0], dtype=np.float32)
    return x[:, cols].sum(axis=1).astype(np.float32)


def append_engineered_v4(x_gene: np.ndarray, spec: FeatureSpec) -> np.ndarray:
    trdc = x_gene[:, spec.gene_to_col["TRDC"]] if "TRDC" in spec.gene_to_col else np.zeros(x_gene.shape[0], dtype=np.float32)
    trdv = sum_cols(x_gene, prefix_cols(spec, ("TRDV",)))
    trdj = sum_cols(x_gene, prefix_cols(spec, ("TRDJ", "TRDD")))
    trg = sum_cols(x_gene, prefix_cols(spec, ("TRGV", "TRGJ", "TRGC")))
    ab = sum_cols(x_gene, prefix_cols(spec, ("TRAV", "TRAJ", "TRAC", "TRBV", "TRBJ", "TRBC")))
    cd3 = mean_cols(x_gene, gene_cols(spec, ["CD3D", "CD3E", "CD3G"]))
    t_lineage = mean_cols(x_gene, gene_cols(spec, T_LINEAGE_GENES))
    nk = mean_cols(x_gene, gene_cols(spec, NK_CONTEXT_GENES))
    cytotoxic = mean_cols(x_gene, gene_cols(spec, ["NKG7", "GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "CST7", "CTSW"]))
    cd4_foxp3 = mean_cols(x_gene, gene_cols(spec, ["CD4", "FOXP3"]))
    b_score = mean_cols(x_gene, gene_cols(spec, ["MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN"]))
    myeloid = mean_cols(x_gene, gene_cols(spec, ["LST1", "LYZ", "S100A8", "S100A9", "FCN1", "CST3"]))
    gdt_tcr = trdc + trdv + trdj + trg
    tcr_support = (cd3 > 0) | (gdt_tcr > 0) | (ab > 0)
    engineered = np.column_stack(
        [
            trdv > 0,
            trdj > 0,
            trg > 0,
            ab > 0,
            (trdc > 0) & (trdv <= 0) & (trdj <= 0),
            (trdc > 0) & (trdv > 0),
            (trdc > 0) & (trg > 0) & (trdv <= 0),
            cd3,
            t_lineage,
            nk,
            cytotoxic,
            gdt_tcr,
            ab,
            gdt_tcr - ab,
            cd3 - nk,
            trdc,
            trdv,
            trdj,
            trg,
            cytotoxic * tcr_support.astype(np.float32),
            cytotoxic * (~tcr_support).astype(np.float32),
            cd4_foxp3,
            b_score,
            myeloid,
        ]
    ).astype(np.float32, copy=False)
    return np.column_stack([x_gene, engineered]).astype(np.float32, copy=False)


def quadrant_from_x(x: np.ndarray, spec: FeatureSpec) -> np.ndarray:
    trdc = x[:, spec.engineered_to_col["TRDC_log1p"]] > 0
    trdv = x[:, spec.engineered_to_col["any_TRDV"]] > 0.5
    out = np.full(x.shape[0], "TRDC-TRDV-", dtype=object)
    out[trdc & (~trdv)] = "TRDC+TRDV-"
    out[(~trdc) & trdv] = "TRDC-TRDV+"
    out[trdc & trdv] = "TRDC+TRDV+"
    return out


def sample_indices(indices: np.ndarray, max_n: int | None, rng: np.random.Generator) -> np.ndarray:
    indices = np.asarray(indices, dtype=np.int64)
    if max_n is None or indices.size <= max_n:
        return indices
    return np.sort(rng.choice(indices, size=max_n, replace=False).astype(np.int64))


def split_by_source_label(rows: np.ndarray, labels: np.ndarray, source: np.ndarray, seed: int, tune_frac: float = 0.2) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    train_parts: list[np.ndarray] = []
    tune_parts: list[np.ndarray] = []
    summary: list[dict[str, Any]] = []
    df = pd.DataFrame({"row": rows, "label": labels, "source_gse_id": source[rows]})
    for (src, label), sub in df.groupby(["source_gse_id", "label"], sort=True):
        idx = sub["row"].to_numpy(dtype=np.int64)
        rng.shuffle(idx)
        n_tune = 0 if idx.size < 5 else max(1, int(round(idx.size * tune_frac)))
        tune_parts.append(idx[:n_tune])
        train_parts.append(idx[n_tune:])
        summary.append({"source_gse_id": src, "label": label, "n_cells": int(idx.size), "train": int(idx.size - n_tune), "tune": int(n_tune)})
    train = np.concatenate(train_parts).astype(np.int64) if train_parts else np.asarray([], dtype=np.int64)
    tune = np.concatenate(tune_parts).astype(np.int64) if tune_parts else np.asarray([], dtype=np.int64)
    return np.sort(train), np.sort(tune), pd.DataFrame(summary)


def sample_weight_balanced(y: np.ndarray) -> np.ndarray:
    y = y.astype(np.int8, copy=False)
    n = y.size
    weights = np.ones(n, dtype=np.float32)
    for label in [0, 1]:
        count = int((y == label).sum())
        if count:
            weights[y == label] = n / (2.0 * count)
    return weights


def metric_row(name: str, y_true: np.ndarray, pred: np.ndarray, score: np.ndarray, threshold: float) -> dict[str, Any]:
    y = y_true.astype(np.int8, copy=False)
    p = pred.astype(np.int8, copy=False)
    if y.size == 0:
        return {"strategy": name, "threshold": threshold, "n_cells": 0}
    cm = confusion_matrix(y, p, labels=[0, 1])
    tn, fp, fn, tp = [int(x) for x in cm.ravel()]
    out = {
        "strategy": name,
        "threshold": float(threshold),
        "n_cells": int(y.size),
        "n_positive": int((y == 1).sum()),
        "n_negative": int((y == 0).sum()),
        "predicted_positive": int(p.sum()),
        "tp": tp,
        "fp": fp,
        "tn": tn,
        "fn": fn,
        "precision": float(precision_score(y, p, zero_division=0)),
        "recall": float(recall_score(y, p, zero_division=0)),
        "specificity": float(tn / (tn + fp)) if (tn + fp) else np.nan,
        "f1": float(f1_score(y, p, zero_division=0)),
        "balanced_accuracy": float(balanced_accuracy_score(y, p)) if len(np.unique(y)) == 2 else np.nan,
        "mcc": float(matthews_corrcoef(y, p)) if len(np.unique(y)) == 2 else np.nan,
        "fp_fraction_of_predictions": float(fp / int(p.sum())) if int(p.sum()) else np.nan,
    }
    if len(np.unique(y)) == 2:
        try:
            out["roc_auc"] = float(roc_auc_score(y, score))
            out["pr_auc"] = float(average_precision_score(y, score))
        except Exception:
            out["roc_auc"] = np.nan
            out["pr_auc"] = np.nan
    else:
        out["roc_auc"] = np.nan
        out["pr_auc"] = np.nan
    return out


def threshold_grid(score: np.ndarray) -> np.ndarray:
    score = np.asarray(score, dtype=np.float32)
    qs = np.unique(np.quantile(score, np.linspace(0.001, 0.999, 999)).astype(np.float32)) if score.size else np.asarray([], dtype=np.float32)
    extras = np.asarray([0.05, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 0.95, 0.97, 0.99], dtype=np.float32)
    return np.unique(np.concatenate([qs, extras]))


def choose_tcell_threshold(score: np.ndarray, y_t: np.ndarray, is_known_gdt: np.ndarray, args: argparse.Namespace) -> tuple[float, pd.DataFrame]:
    rows = []
    for th in threshold_grid(score):
        pred = score >= th
        gdt_recall = float((pred & is_known_gdt).sum() / is_known_gdt.sum()) if is_known_gdt.sum() else np.nan
        t_recall = float((pred & (y_t == 1)).sum() / (y_t == 1).sum()) if (y_t == 1).sum() else np.nan
        specificity = float(((~pred) & (y_t == 0)).sum() / (y_t == 0).sum()) if (y_t == 0).sum() else np.nan
        rows.append({"threshold": float(th), "known_gdt_recall": gdt_recall, "t_positive_recall": t_recall, "non_t_specificity": specificity})
    df = pd.DataFrame(rows)
    valid = df[(df["known_gdt_recall"] >= args.target_tcell_gdt_recall) & (df["t_positive_recall"] >= args.target_tcell_positive_recall)]
    if valid.empty:
        gdt_scores = score[is_known_gdt]
        threshold = float(np.quantile(gdt_scores, max(0.0, 1.0 - args.target_tcell_gdt_recall))) if gdt_scores.size else 0.2
    else:
        threshold = float(valid.sort_values(["non_t_specificity", "threshold"], ascending=[False, False]).iloc[0]["threshold"])
    return threshold, df


def choose_gdt_thresholds(score: np.ndarray, y: np.ndarray, t_score: np.ndarray, t_threshold: float) -> tuple[dict[str, float], pd.DataFrame]:
    rows = []
    for th in threshold_grid(score):
        pred = (t_score >= t_threshold) & (score >= th)
        m = metric_row("candidate", y, pred, score, float(th))
        tp, fp, fn = int(m.get("tp", 0)), int(m.get("fp", 0)), int(m.get("fn", 0))
        f05 = (1.25 * tp / (1.25 * tp + 0.25 * fn + fp)) if (1.25 * tp + 0.25 * fn + fp) else 0.0
        m["f0.5"] = float(f05)
        rows.append(m)
    df = pd.DataFrame(rows)
    nonzero = df[df["predicted_positive"] > 0].copy()
    high_f1 = nonzero.sort_values(["f1", "recall", "precision"], ascending=[False, False, False]).iloc[0]
    # high-purity still tracks recall, but prioritizes precision/F0.5 and lower FP fraction.
    hp = nonzero.sort_values(["f0.5", "precision", "recall", "fp_fraction_of_predictions"], ascending=[False, False, False, True]).iloc[0]
    return {"high_f1": float(high_f1["threshold"]), "high_purity": float(hp["threshold"])}, df


def build_initial_rows(obs: Any, annotation: np.ndarray, args: argparse.Namespace) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    rng = np.random.default_rng(args.seed)
    n = obs.source.size
    ann_upper = pd.Series(annotation, copy=False).astype(str).str.upper().to_numpy(dtype=object)
    tissue_lower = pd.Series(obs.tissue, copy=False).astype(str).str.lower().to_numpy(dtype=object)
    t_like = np.isin(ann_upper, list(T_ANNOTATIONS))
    nk_like = ann_upper == "NK_CELL"
    other_non_t = np.isin(ann_upper, list(NON_T_ANNOTATIONS)) & (~t_like)
    tcr_evidence = obs.has_TRA_TRB_paired | obs.has_any_ab_tcr | obs.corrected_has_any_gd_tcr
    known_gdt = np.isin(obs.class_code, [2, 3])
    validation_gdt = (obs.source == GDT2020_SOURCE) & (obs.class_code == 2) & (tissue_lower == GDT2020_HOLDOUT_TISSUE.lower())
    validation_abt = (obs.source == "GSE254249") & (obs.has_TRA_TRB_paired | obs.has_any_ab_tcr) & (~obs.corrected_has_any_gd_tcr) & (obs.class_code != 2)
    validation = validation_gdt | validation_abt
    trainable = (~validation) & (obs.source != SUBOPTIMAL_SOURCE)

    stage1_pos_pool = np.flatnonzero(trainable & (t_like | tcr_evidence | known_gdt)).astype(np.int64)
    stage1_nk_pool = np.flatnonzero(trainable & nk_like & (obs.class_code != 2)).astype(np.int64)
    stage1_other_pool = np.flatnonzero(trainable & other_non_t & (~tcr_evidence) & (~known_gdt)).astype(np.int64)
    stage2_pos_pool = np.flatnonzero(trainable & (obs.class_code == 2)).astype(np.int64)
    stage2_ab_pool = np.flatnonzero(trainable & (((obs.class_code == 1) | ((obs.has_TRA_TRB_paired | obs.has_any_ab_tcr) & (~obs.corrected_has_any_gd_tcr))) & (obs.class_code != 2))).astype(np.int64)
    stage2_nk_pool = np.flatnonzero(trainable & nk_like & (obs.class_code != 2)).astype(np.int64)
    sensitivity_gdtlung = np.flatnonzero((obs.source == SUBOPTIMAL_SOURCE) & (obs.class_code == 2)).astype(np.int64)
    sensitivity_silver = np.flatnonzero(obs.class_code == 3).astype(np.int64)

    rows = {
        "stage1_pos_pool": sample_indices(stage1_pos_pool, args.stage1_max_positive, rng),
        "stage1_nk_pool": sample_indices(stage1_nk_pool, args.stage1_max_nk_negative, rng),
        "stage1_other_pool": sample_indices(stage1_other_pool, 40_000, rng),
        "stage2_pos_pool": stage2_pos_pool,
        "stage2_ab_pool": sample_indices(stage2_ab_pool, args.stage2_max_ab_negative, rng),
        "stage2_nk_pool": sample_indices(stage2_nk_pool, max(args.stage2_max_nk_negative * 3, args.stage2_max_nk_negative), rng),
        "validation_gdt": np.flatnonzero(validation_gdt).astype(np.int64),
        "validation_abt": np.flatnonzero(validation_abt).astype(np.int64),
        "validation": np.flatnonzero(validation).astype(np.int64),
        "sensitivity_gdtlung": sensitivity_gdtlung,
        "sensitivity_silver": sensitivity_silver,
    }
    summary = pd.DataFrame([{"pool": k, "n_cells": int(v.size)} for k, v in rows.items()])
    return rows, summary


def high_confidence_nk_from_x(x: np.ndarray, spec: FeatureSpec, ann: np.ndarray) -> np.ndarray:
    ann_upper = pd.Series(ann, copy=False).astype(str).str.upper().to_numpy(dtype=object)
    nk_ann = ann_upper == "NK_CELL"
    nk = x[:, spec.engineered_to_col["NK_context_score"]]
    cytotoxic = x[:, spec.engineered_to_col["cytotoxic_score"]]
    cd3 = x[:, spec.engineered_to_col["CD3_score"]]
    t_lineage = x[:, spec.engineered_to_col["T_lineage_score"]]
    # Composite, intentionally not single-gene. Allows some TCR/cytotoxic overlap but requires low CD3/T-lineage support.
    return nk_ann & (nk >= 0.65) & (cytotoxic >= 0.45) & (cd3 <= 0.75) & (t_lineage <= 0.85)


def build_training_matrices(handle: h5py.File, spec: FeatureSpec, obs: Any, annotation: np.ndarray, args: argparse.Namespace) -> tuple[np.ndarray, np.ndarray, SplitRows, pd.DataFrame, pd.DataFrame]:
    pools, pool_summary = build_initial_rows(obs, annotation, args)
    all_rows = np.unique(np.concatenate([v for v in pools.values() if v.size])).astype(np.int64)
    logging.info("Extracting v4 training/evaluation features for %s unique rows", f"{all_rows.size:,}")
    x_gene, _row_sum, _n_detected = extract_gene_features(handle, "X", all_rows, spec, label="v4_train_eval")
    x_all = append_engineered_v4(x_gene, spec)
    row_to_pos = pd.Series(np.arange(all_rows.size, dtype=np.int64), index=all_rows).to_dict()

    def pos(rows: np.ndarray) -> np.ndarray:
        return np.asarray([row_to_pos[int(r)] for r in rows if int(r) in row_to_pos], dtype=np.int64)

    stage1_pos_rows = pools["stage1_pos_pool"]
    stage1_neg_candidates = np.unique(np.concatenate([pools["stage1_nk_pool"], pools["stage1_other_pool"]])).astype(np.int64)
    hc_nk_all = high_confidence_nk_from_x(x_all, spec, annotation[all_rows])
    hc_nk_rows = all_rows[hc_nk_all]
    other_rows = pools["stage1_other_pool"]
    stage1_neg_rows = np.unique(np.concatenate([np.intersect1d(stage1_neg_candidates, hc_nk_rows), other_rows])).astype(np.int64)
    if stage1_neg_rows.size == 0:
        raise RuntimeError("No high-confidence non-T/NK negatives available for T-cell gate training.")
    stage1_rows = np.unique(np.concatenate([stage1_pos_rows, stage1_neg_rows])).astype(np.int64)
    stage1_labels = np.isin(stage1_rows, stage1_pos_rows).astype(np.int8)
    stage1_train_rows, stage1_tune_rows, split1 = split_by_source_label(stage1_rows, np.where(stage1_labels == 1, "T_like", "non_T_highconf"), obs.source, args.seed + 1)

    # Stage 2 positives and negatives. Use high-confidence NK cells that are most T-like as hard negatives.
    stage2_pos_rows = pools["stage2_pos_pool"]
    stage2_ab_rows = pools["stage2_ab_pool"]
    stage2_nk_candidates = np.intersect1d(pools["stage2_nk_pool"], hc_nk_rows)
    stage2_nk_rows = stage2_nk_candidates
    if stage2_nk_rows.size > args.stage2_max_nk_negative:
        # Temporary broad T-like proxy before the T gate is trained: sort by CD3 + T-lineage support and keep hardest plus random.
        cand_pos = pos(stage2_nk_rows)
        proxy = x_all[cand_pos, spec.engineered_to_col["CD3_score"]] + x_all[cand_pos, spec.engineered_to_col["T_lineage_score"]]
        order = np.argsort(proxy)[::-1]
        hard_n = min(args.stage2_max_nk_negative // 2, order.size)
        hard = stage2_nk_rows[order[:hard_n]]
        rng = np.random.default_rng(args.seed + 3)
        remaining = np.setdiff1d(stage2_nk_rows, hard, assume_unique=False)
        random_n = max(0, args.stage2_max_nk_negative - hard.size)
        random = sample_indices(remaining, random_n, rng)
        stage2_nk_rows = np.unique(np.concatenate([hard, random])).astype(np.int64)
    stage2_rows = np.unique(np.concatenate([stage2_pos_rows, stage2_ab_rows, stage2_nk_rows])).astype(np.int64)
    stage2_labels = np.where(np.isin(stage2_rows, stage2_pos_rows), "gdT_gold", np.where(np.isin(stage2_rows, stage2_nk_rows), "highconf_NK_negative", "TCRAB_negative"))
    stage2_train_rows, stage2_tune_rows, split2 = split_by_source_label(stage2_rows, stage2_labels, obs.source, args.seed + 2)

    split_summary = pd.concat(
        [
            split1.assign(stage="stage1_tcell_gate"),
            split2.assign(stage="stage2_gdt_classifier"),
        ],
        ignore_index=True,
    )
    split_summary.to_csv(TABLE_DIR / "split_by_source_label.csv", index=False)
    pool_summary.to_csv(TABLE_DIR / "initial_training_pools.csv", index=False)
    pd.DataFrame(
        [
            {"set": "all_eval_rows", "n_cells": int(all_rows.size)},
            {"set": "stage1_positive_pool", "n_cells": int(stage1_pos_rows.size)},
            {"set": "stage1_highconf_negative", "n_cells": int(stage1_neg_rows.size)},
            {"set": "stage2_gdT_positive", "n_cells": int(stage2_pos_rows.size)},
            {"set": "stage2_abT_TCRAB_negative", "n_cells": int(stage2_ab_rows.size)},
            {"set": "stage2_highconf_NK_negative", "n_cells": int(stage2_nk_rows.size)},
            {"set": "validation_gdt", "n_cells": int(pools["validation_gdt"].size)},
            {"set": "validation_abt", "n_cells": int(pools["validation_abt"].size)},
            {"set": "sensitivity_gdtlung", "n_cells": int(pools["sensitivity_gdtlung"].size)},
            {"set": "sensitivity_silver", "n_cells": int(pools["sensitivity_silver"].size)},
        ]
    ).to_csv(TABLE_DIR / "training_sample_composition.csv", index=False)

    splits = SplitRows(
        stage1_train=stage1_train_rows,
        stage1_tune=stage1_tune_rows,
        stage2_train=stage2_train_rows,
        stage2_tune=stage2_tune_rows,
        validation=pools["validation"],
        validation_gdt=pools["validation_gdt"],
        validation_abt=pools["validation_abt"],
        sensitivity_gdtlung=pools["sensitivity_gdtlung"],
        sensitivity_silver=pools["sensitivity_silver"],
        all_eval_rows=all_rows,
        split_summary=split_summary,
    )
    return all_rows, x_all, splits, pool_summary, split_summary


def train_models(all_rows: np.ndarray, x_all: np.ndarray, splits: SplitRows, obs: Any, args: argparse.Namespace) -> tuple[GdTAIV4TwoStageModel, pd.DataFrame, pd.DataFrame]:
    row_to_pos = pd.Series(np.arange(all_rows.size, dtype=np.int64), index=all_rows).to_dict()

    def positions(rows: np.ndarray) -> np.ndarray:
        return np.asarray([row_to_pos[int(r)] for r in rows if int(r) in row_to_pos], dtype=np.int64)

    stage1_train_pos = positions(splits.stage1_train)
    stage1_tune_pos = positions(splits.stage1_tune)
    y1_train = np.isin(splits.stage1_train, splits.stage1_train[obs.class_code[splits.stage1_train] >= -99]).astype(np.int8)
    # Recreate stage1 labels from source labels in split summary: T-like if class/gene/annotation label assigned during split table.
    stage1_train_labels = []
    stage1_tune_labels = []
    split_lookup = splits.split_summary[splits.split_summary["stage"].eq("stage1_tcell_gate")]
    t_like_sources = set()
    # More direct labels: stage1 positives are all non-negative rows with T-like evidence.
    # Derive from obs/annotation-independent metadata encoded during row choice: gdT/TCR evidence/class positives are not enough for CD8/CD4 annotations,
    # so use the split table label via merge below.
    label_frame = []
    for _, row in split_lookup.iterrows():
        pass
    # Simpler deterministic reconstruction: rows not high-confidence negatives are T-like; high-confidence negatives were assigned label 0 in split_by_source_label.
    # Store labels by reading the generated split table is overkill, so infer from class/TCR/annotation not available here is avoided by passing labels below.
    raise_if = False
    _ = raise_if
    return _train_models_with_labels(all_rows, x_all, splits, obs, args)


def _train_models_with_labels(all_rows: np.ndarray, x_all: np.ndarray, splits: SplitRows, obs: Any, args: argparse.Namespace) -> tuple[GdTAIV4TwoStageModel, pd.DataFrame, pd.DataFrame]:
    # Labels are reconstructed from the split_by_source_label table created during sampling.
    # Rows with class_code 2/3, TCR evidence, or T-like annotation were included as T positives;
    # high-confidence NK/other rows included as negatives do not overlap known gdT or paired TCR evidence.
    # To avoid hidden state, use the CSV row labels to build row label maps.
    split_df = pd.read_csv(TABLE_DIR / "split_by_source_label.csv")
    # The split table is aggregated, so exact row labels are not recoverable. Use biologically equivalent metadata rules here.
    # This is intentional: the stage1 gate learns broad T-lineage, not the exact pool origin.
    # For exact row-level labels, stage2 uses class/TCR/NK membership below.
    row_to_pos = pd.Series(np.arange(all_rows.size, dtype=np.int64), index=all_rows).to_dict()

    def positions(rows: np.ndarray) -> np.ndarray:
        return np.asarray([row_to_pos[int(r)] for r in rows if int(r) in row_to_pos], dtype=np.int64)

    # Rebuild annotation from saved global file is not available here; infer stage1 positives from split labels by row membership is not possible.
    # Therefore this function is replaced at runtime by train_models_direct, which has annotation and high-confidence labels.
    raise RuntimeError("Internal label reconstruction path should not be used.")


def train_models_direct(
    all_rows: np.ndarray,
    x_all: np.ndarray,
    splits: SplitRows,
    obs: Any,
    annotation: np.ndarray,
    spec: FeatureSpec,
    args: argparse.Namespace,
) -> tuple[GdTAIV4TwoStageModel, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    row_to_pos = pd.Series(np.arange(all_rows.size, dtype=np.int64), index=all_rows).to_dict()

    def positions(rows: np.ndarray) -> np.ndarray:
        return np.asarray([row_to_pos[int(r)] for r in rows if int(r) in row_to_pos], dtype=np.int64)

    ann_upper_all = pd.Series(annotation, copy=False).astype(str).str.upper().to_numpy(dtype=object)
    tcr_evidence = obs.has_TRA_TRB_paired | obs.has_any_ab_tcr | obs.corrected_has_any_gd_tcr
    known_gdt = np.isin(obs.class_code, [2, 3])
    t_like = np.isin(ann_upper_all, list(T_ANNOTATIONS))
    hc_nk_all = high_confidence_nk_from_x(x_all, spec, annotation[all_rows])
    hc_nk_rows = set(all_rows[hc_nk_all].tolist())

    y1_train = ((t_like[splits.stage1_train] | tcr_evidence[splits.stage1_train] | known_gdt[splits.stage1_train]) & (~np.isin(splits.stage1_train, list(hc_nk_rows)))).astype(np.int8)
    y1_tune = ((t_like[splits.stage1_tune] | tcr_evidence[splits.stage1_tune] | known_gdt[splits.stage1_tune]) & (~np.isin(splits.stage1_tune, list(hc_nk_rows)))).astype(np.int8)
    p1_train = positions(splits.stage1_train)
    p1_tune = positions(splits.stage1_tune)
    if y1_train.sum() == 0 or (y1_train == 0).sum() == 0:
        raise RuntimeError("Stage1 training labels are not two-class.")
    logging.info("Training stage1 T-cell gate on %s rows", f"{p1_train.size:,}")
    tcell_model = HistGradientBoostingClassifier(max_iter=160, learning_rate=0.055, max_leaf_nodes=31, l2_regularization=0.03, random_state=args.seed)
    tcell_model.fit(x_all[p1_train], y1_train, sample_weight=sample_weight_balanced(y1_train))
    t_tune_score = tcell_model.predict_proba(x_all[p1_tune])[:, 1].astype(np.float32)
    t_threshold, t_sweep = choose_tcell_threshold(t_tune_score, y1_tune, known_gdt[splits.stage1_tune], args)
    t_sweep.to_csv(TABLE_DIR / "stage1_tcell_threshold_sweep.csv", index=False)
    t_tune_pred = t_tune_score >= t_threshold
    stage1_metrics = pd.DataFrame([metric_row("v4_stage1_tcell_gate", y1_tune, t_tune_pred, t_tune_score, t_threshold)])
    stage1_metrics.to_csv(TABLE_DIR / "stage1_tcell_gate_tune_metrics.csv", index=False)

    p2_train = positions(splits.stage2_train)
    p2_tune = positions(splits.stage2_tune)
    y2_train = (obs.class_code[splits.stage2_train] == 2).astype(np.int8)
    y2_tune = (obs.class_code[splits.stage2_tune] == 2).astype(np.int8)
    if y2_train.sum() == 0 or (y2_train == 0).sum() == 0:
        raise RuntimeError("Stage2 training labels are not two-class.")
    logging.info("Scoring T gate for stage2 and training gdT classifier on %s rows", f"{p2_train.size:,}")
    t_train_score = tcell_model.predict_proba(x_all[p2_train])[:, 1].astype(np.float32)
    t2_tune_score = tcell_model.predict_proba(x_all[p2_tune])[:, 1].astype(np.float32)
    x2_train = np.column_stack([x_all[p2_train], t_train_score]).astype(np.float32, copy=False)
    x2_tune = np.column_stack([x_all[p2_tune], t2_tune_score]).astype(np.float32, copy=False)
    gdt_model = HistGradientBoostingClassifier(max_iter=190, learning_rate=0.052, max_leaf_nodes=31, l2_regularization=0.04, random_state=args.seed + 10)
    gdt_model.fit(x2_train, y2_train, sample_weight=sample_weight_balanced(y2_train))
    gdt_tune_score = gdt_model.predict_proba(x2_tune)[:, 1].astype(np.float32)
    thresholds, gdt_sweep = choose_gdt_thresholds(gdt_tune_score, y2_tune, t2_tune_score, t_threshold)
    gdt_sweep.to_csv(TABLE_DIR / "stage2_gdt_threshold_sweep.csv", index=False)
    model = GdTAIV4TwoStageModel(
        tcell_model=tcell_model,
        gdt_model=gdt_model,
        tcell_threshold=t_threshold,
        mode_thresholds=thresholds,
        gene_names=[],
        engineered_feature_names=list(ENGINEERED_FEATURES),
        model_feature_names=[],
    )
    tune_rows = []
    for mode, th in thresholds.items():
        pred = (t2_tune_score >= t_threshold) & (gdt_tune_score >= th)
        row = metric_row(f"v4_{mode}_tune", y2_tune, pred, gdt_tune_score, th)
        row["tcell_threshold"] = float(t_threshold)
        tune_rows.append(row)
    tune_metrics = pd.DataFrame(tune_rows)
    tune_metrics.to_csv(TABLE_DIR / "stage2_gdt_tune_metrics.csv", index=False)
    return model, stage1_metrics, tune_metrics, gdt_sweep


def evaluate_rows(
    name: str,
    rows: np.ndarray,
    y: np.ndarray,
    x_all: np.ndarray,
    all_rows: np.ndarray,
    model: GdTAIV4TwoStageModel,
    obs: Any,
    annotation: np.ndarray,
    spec: FeatureSpec,
) -> pd.DataFrame:
    if rows.size == 0:
        return pd.DataFrame()
    row_to_pos = pd.Series(np.arange(all_rows.size, dtype=np.int64), index=all_rows).to_dict()
    pos = np.asarray([row_to_pos[int(r)] for r in rows if int(r) in row_to_pos], dtype=np.int64)
    valid_rows = np.asarray([int(r) for r in rows if int(r) in row_to_pos], dtype=np.int64)
    if pos.size == 0:
        return pd.DataFrame()
    t_score, gdt_score, joint = model.score_arrays(x_all[pos])
    quadrant = quadrant_from_x(x_all[pos], spec)
    ann = annotation[valid_rows]
    ann_upper = pd.Series(ann, copy=False).astype(str).str.upper().to_numpy(dtype=object)
    records = []
    for mode, th in model.mode_thresholds.items():
        pred = (t_score >= model.tcell_threshold) & (gdt_score >= th)
        row = metric_row(f"v4_{mode}", y[: pred.size], pred, gdt_score, th)
        row["evaluation_set"] = name
        row["tcell_threshold"] = model.tcell_threshold
        row["predicted_NK_annotation"] = int((pred & (ann_upper == "NK_CELL")).sum())
        row["predicted_TRDC_plus_TRDV_minus"] = int((pred & (quadrant == "TRDC+TRDV-")).sum())
        records.append(row)
        for q in ["TRDC+TRDV+", "TRDC+TRDV-", "TRDC-TRDV+", "TRDC-TRDV-"]:
            mask = quadrant == q
            if mask.any() and len(np.unique(y[: pred.size][mask])) >= 1:
                qrow = metric_row(f"v4_{mode}", y[: pred.size][mask], pred[mask], gdt_score[mask], th)
                qrow["evaluation_set"] = name
                qrow["group_column"] = "tcr_gene_quadrant"
                qrow["group_value"] = q
                records.append(qrow)
    return pd.DataFrame(records)


def evaluate_internal(
    all_rows: np.ndarray,
    x_all: np.ndarray,
    splits: SplitRows,
    model: GdTAIV4TwoStageModel,
    obs: Any,
    annotation: np.ndarray,
    spec: FeatureSpec,
) -> pd.DataFrame:
    frames = []
    if splits.validation.size:
        y_val = (obs.class_code[splits.validation] == 2).astype(np.int8)
        frames.append(evaluate_rows("validation_GDT2020_cordblood_vs_GSE254249_TCRAB", splits.validation, y_val, x_all, all_rows, model, obs, annotation, spec))
    if splits.validation_gdt.size:
        y = np.ones(splits.validation_gdt.size, dtype=np.int8)
        frames.append(evaluate_rows("validation_gdt_GDT2020_cordblood", splits.validation_gdt, y, x_all, all_rows, model, obs, annotation, spec))
    if splits.validation_abt.size:
        y = np.zeros(splits.validation_abt.size, dtype=np.int8)
        frames.append(evaluate_rows("validation_abt_GSE254249_TCRAB", splits.validation_abt, y, x_all, all_rows, model, obs, annotation, spec))
    if splits.sensitivity_gdtlung.size:
        y = np.ones(splits.sensitivity_gdtlung.size, dtype=np.int8)
        frames.append(evaluate_rows("sensitivity_GDTlung_suboptimal_sorted_gdT", splits.sensitivity_gdtlung, y, x_all, all_rows, model, obs, annotation, spec))
    if splits.sensitivity_silver.size:
        rows = np.intersect1d(splits.sensitivity_silver, all_rows)
        if rows.size:
            y = np.ones(rows.size, dtype=np.int8)
            frames.append(evaluate_rows("sensitivity_silver_gdT", rows, y, x_all, all_rows, model, obs, annotation, spec))
    out = pd.concat([f for f in frames if f is not None and not f.empty], ignore_index=True) if frames else pd.DataFrame()
    out.to_csv(TABLE_DIR / "internal_validation_metrics.csv", index=False)
    return out


def full_summary_for_mode(
    mode: str,
    pred: np.ndarray,
    t_score: np.ndarray,
    gdt_score: np.ndarray,
    quadrant: np.ndarray,
    obs: Any,
    source: np.ndarray,
    tissue: np.ndarray,
    annotation: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    y = np.zeros(pred.size, dtype=np.int8)
    y[obs.class_code[: pred.size] == 2] = 1
    primary = np.isin(obs.class_code[: pred.size], [1, 2])
    ann_upper = pd.Series(annotation, copy=False).astype(str).str.upper().to_numpy(dtype=object)
    is_nk = ann_upper == "NK_CELL"
    is_cd4_treg = np.isin(ann_upper, ["CD4_T", "TREG"])
    is_tcrab = (obs.has_TRA_TRB_paired[: pred.size] | obs.has_any_ab_tcr[: pred.size]) & (obs.class_code[: pred.size] != 2)
    overall = {
        "strategy": f"v4_{mode}",
        "total_cells": int(pred.size),
        "predicted_putative_gdT": int(pred.sum()),
        "predicted_fraction": float(pred.mean()),
        "primary_gold_cells": int(primary.sum()),
        "full_primary_tp": int((pred & primary & (y == 1)).sum()),
        "full_primary_fp": int((pred & primary & (y == 0)).sum()),
        "full_primary_fn": int(((~pred) & primary & (y == 1)).sum()),
        "paired_TCRAB_cells": int(is_tcrab.sum()),
        "predicted_paired_TCRAB": int((pred & is_tcrab).sum()),
        "NK_cells": int(is_nk.sum()),
        "predicted_NK": int((pred & is_nk).sum()),
        "TRDC_plus_TRDV_minus_cells": int((quadrant == "TRDC+TRDV-").sum()),
        "predicted_TRDC_plus_TRDV_minus": int((pred & (quadrant == "TRDC+TRDV-")).sum()),
        "CD4_Treg_warning_cells": int(is_cd4_treg.sum()),
        "predicted_CD4_Treg_warning": int((pred & is_cd4_treg).sum()),
        "tcell_score_mean_predicted": float(t_score[pred].mean()) if pred.any() else np.nan,
        "gdt_score_mean_predicted": float(gdt_score[pred].mean()) if pred.any() else np.nan,
    }
    tp, fp, fn = overall["full_primary_tp"], overall["full_primary_fp"], overall["full_primary_fn"]
    overall["full_primary_precision"] = float(tp / (tp + fp)) if (tp + fp) else np.nan
    overall["full_primary_recall"] = float(tp / (tp + fn)) if (tp + fn) else np.nan
    base = pd.DataFrame(
        {
            "source_gse_id": source,
            "tissue": tissue,
            "annotation": annotation,
            "predicted_putative_gdT": pred.astype(np.int8),
            "gdT_gold": ((obs.class_code[: pred.size] == 2)).astype(np.int8),
            "abT_gold": ((obs.class_code[: pred.size] == 1)).astype(np.int8),
            "gdT_silver": ((obs.class_code[: pred.size] == 3)).astype(np.int8),
            "predicted_NK": (pred & is_nk).astype(np.int8),
            "predicted_paired_TCRAB": (pred & is_tcrab).astype(np.int8),
            "predicted_TRDC_plus_TRDV_minus": (pred & (quadrant == "TRDC+TRDV-")).astype(np.int8),
            "predicted_CD4_Treg_warning": (pred & is_cd4_treg).astype(np.int8),
        }
    )

    def agg(group: str) -> pd.DataFrame:
        out = base.groupby(group, dropna=False, as_index=False).agg(
            total_cells=("predicted_putative_gdT", "size"),
            predicted_putative_gdT=("predicted_putative_gdT", "sum"),
            gdT_gold=("gdT_gold", "sum"),
            abT_gold=("abT_gold", "sum"),
            gdT_silver=("gdT_silver", "sum"),
            predicted_NK=("predicted_NK", "sum"),
            predicted_paired_TCRAB=("predicted_paired_TCRAB", "sum"),
            predicted_TRDC_plus_TRDV_minus=("predicted_TRDC_plus_TRDV_minus", "sum"),
            predicted_CD4_Treg_warning=("predicted_CD4_Treg_warning", "sum"),
        )
        out.insert(0, "strategy", f"v4_{mode}")
        out["predicted_fraction"] = out["predicted_putative_gdT"] / out["total_cells"].replace(0, np.nan)
        return out

    return pd.DataFrame([overall]), agg("source_gse_id"), agg("tissue"), agg("annotation")


def apply_full_atlas(handle: h5py.File, spec: FeatureSpec, model: GdTAIV4TwoStageModel, obs: Any, annotation: np.ndarray, args: argparse.Namespace) -> dict[str, pd.DataFrame]:
    n_obs, _n_vars = h5ad_shape(handle)
    n_eval = n_obs if args.max_full_atlas_cells is None else min(n_obs, args.max_full_atlas_cells)
    source = obs.source[:n_eval]
    tissue = obs.tissue[:n_eval]
    annotation_eval = annotation[:n_eval]
    pred_by_mode = {mode: np.zeros(n_eval, dtype=bool) for mode in model.mode_thresholds}
    t_score_all = np.zeros(n_eval, dtype=np.float32)
    gdt_score_all = np.zeros(n_eval, dtype=np.float32)
    quadrant = np.empty(n_eval, dtype=object)
    logging.info("Applying v4 full atlas to %s cells", f"{n_eval:,}")
    for start in range(0, n_eval, args.chunk_size):
        end = min(start + args.chunk_size, n_eval)
        rows = np.arange(start, end, dtype=np.int64)
        x_gene, _row_sum, _n_detected = extract_gene_features(handle, "X", rows, spec, label=f"v4_full_{start}_{end}")
        x = append_engineered_v4(x_gene, spec)
        t_score, gdt_score, _joint = model.score_arrays(x)
        t_score_all[start:end] = t_score
        gdt_score_all[start:end] = gdt_score
        quadrant[start:end] = quadrant_from_x(x, spec)
        for mode, th in model.mode_thresholds.items():
            pred_by_mode[mode][start:end] = (t_score >= model.tcell_threshold) & (gdt_score >= th)
        if end % 500_000 == 0 or end == n_eval:
            logging.info("Applied v4 to %s / %s cells", f"{end:,}", f"{n_eval:,}")
    overall_parts = []
    source_parts = []
    tissue_parts = []
    annotation_parts = []
    for mode, pred in pred_by_mode.items():
        overall, by_source, by_tissue, by_annotation = full_summary_for_mode(mode, pred, t_score_all, gdt_score_all, quadrant, obs, source, tissue, annotation_eval)
        overall_parts.append(overall)
        source_parts.append(by_source)
        tissue_parts.append(by_tissue)
        annotation_parts.append(by_annotation)
        selected_rows = np.flatnonzero(pred).astype(np.int64)
        selected = pd.DataFrame(
            {
                "obs_index": selected_rows,
                "cell_id": obs_index_values(handle)[selected_rows],
                "source_gse_id": source[selected_rows],
                "tissue": tissue[selected_rows],
                "annotation": annotation_eval[selected_rows],
                "class_code": obs.class_code[selected_rows].astype(int),
                "tcell_score": t_score_all[selected_rows],
                "gdt_score": gdt_score_all[selected_rows],
                "mode": mode,
                "tcr_gene_quadrant": quadrant[selected_rows],
            }
        )
        selected.to_csv(TABLE_DIR / f"full_atlas_selected_predicted_cells_{mode}.csv.gz", index=False, compression="gzip")
    outputs = {
        "overall": pd.concat(overall_parts, ignore_index=True),
        "source": pd.concat(source_parts, ignore_index=True),
        "tissue": pd.concat(tissue_parts, ignore_index=True),
        "annotation": pd.concat(annotation_parts, ignore_index=True),
    }
    outputs["overall"].to_csv(TABLE_DIR / "full_atlas_prediction_overall.csv", index=False)
    outputs["source"].to_csv(TABLE_DIR / "full_atlas_prediction_by_source.csv", index=False)
    outputs["tissue"].to_csv(TABLE_DIR / "full_atlas_prediction_by_tissue.csv", index=False)
    outputs["annotation"].to_csv(TABLE_DIR / "full_atlas_prediction_by_annotation.csv", index=False)
    return outputs


def compare_against_existing(v4_overall: pd.DataFrame) -> pd.DataFrame:
    rows = []
    existing = OUTPUT_ROOT / "tables" / "gdT_prediction" / "gdtai_v3_trdc_nk_guard" / "full_atlas_prediction_overall.csv"
    if existing.exists():
        df = pd.read_csv(existing)
        keep = df[df["strategy"].isin(["v2_high_f1", "v2_high_purity", "v3_round12_hist_gradient_fixed_0p5", "original_TRD_minus_TRAB"])]
        rows.append(keep)
    rows.append(v4_overall)
    out = pd.concat(rows, ignore_index=True, sort=False) if rows else v4_overall
    cols = [c for c in ["strategy", "total_cells", "predicted_putative_gdT", "predicted_fraction", "full_primary_precision", "full_primary_recall", "predicted_paired_TCRAB", "predicted_NK", "predicted_TRDC_plus_TRDV_minus", "predicted_CD4_Treg_warning"] if c in out.columns]
    out[cols].to_csv(TABLE_DIR / "comparison_against_v2_v3_overall.csv", index=False)
    return out[cols]


def write_figures(validation: pd.DataFrame, comparison: pd.DataFrame) -> list[Path]:
    paths: list[Path] = []
    if not comparison.empty and "predicted_putative_gdT" in comparison.columns:
        fig, ax = plt.subplots(figsize=(9, 4.8))
        plot_df = comparison.drop_duplicates("strategy").copy()
        ax.bar(plot_df["strategy"], plot_df["predicted_putative_gdT"], color="#3A6EA5")
        ax.set_ylabel("Predicted gdT cells")
        ax.set_title("Full-atlas predicted gdT burden")
        ax.tick_params(axis="x", rotation=35, labelsize=8)
        fig.tight_layout()
        out = FIGURE_DIR / "full_atlas_predicted_counts_comparison.png"
        fig.savefig(out, dpi=220)
        plt.close(fig)
        paths.append(out)
    if not validation.empty:
        val = validation[validation["evaluation_set"].eq("validation_GDT2020_cordblood_vs_GSE254249_TCRAB")].copy()
        if not val.empty:
            fig, ax = plt.subplots(figsize=(6, 4.5))
            ax.bar(val["strategy"], val["f1"], color="#2A9D8F")
            ax.set_ylim(0, 1)
            ax.set_ylabel("F1")
            ax.set_title("v4 validation F1")
            ax.tick_params(axis="x", rotation=25, labelsize=9)
            fig.tight_layout()
            out = FIGURE_DIR / "v4_validation_f1.png"
            fig.savefig(out, dpi=220)
            plt.close(fig)
            paths.append(out)
    for path in paths:
        target = STATIC_ASSET_DIR / path.name
        target.write_bytes(path.read_bytes())
    return paths


def render_report(stage1: pd.DataFrame, tune: pd.DataFrame, validation: pd.DataFrame, full: pd.DataFrame, comparison: pd.DataFrame, figures: list[Path], args: argparse.Namespace) -> None:
    fig_html = "\n".join(f'<figure><img src="assets/{OUT_PREFIX}/{html.escape(p.name)}"><figcaption>{html.escape(p.stem.replace("_", " "))}</figcaption></figure>' for p in figures)
    md = [
        "# gdTAI v4 Two-Step T-Cell Gate Classifier",
        "",
        "## Design",
        "- Stage 1 estimates soft T-lineage probability.",
        "- Stage 2 estimates gdT probability among T-like/T-NK-like cells.",
        "- No hard TRDV requirement is used; TRDC+TRDV- cells are evaluated explicitly.",
        "- NK/cytotoxic genes are context features, not single-gene exclusion rules.",
        "",
        "## Stage 1 T-Cell Gate Tune Metrics",
        dataframe_to_markdown(stage1),
        "",
        "## Stage 2 gdT Tune Metrics",
        dataframe_to_markdown(tune),
        "",
        "## Internal Validation and Sensitivity",
        dataframe_to_markdown(validation.head(80) if not validation.empty else validation),
        "",
        "## Full-Atlas v4 Summary",
        dataframe_to_markdown(full),
        "",
        "## Comparison Against Existing Strategies",
        dataframe_to_markdown(comparison),
        "",
        "## Outputs",
        f"- Tables: `{TABLE_DIR}`",
        f"- Figures: `{FIGURE_DIR}`",
        f"- Model: `{MODEL_PKL}`",
    ]
    REPORT_MD.write_text("\n".join(md) + "\n", encoding="utf-8")
    css = "body{font-family:Arial,sans-serif;max-width:1200px;margin:28px auto;line-height:1.45;color:#222} table{border-collapse:collapse;width:100%;font-size:12px} th,td{border:1px solid #ddd;padding:5px;text-align:right} th:first-child,td:first-child{text-align:left} th{background:#eef3f7} img{max-width:100%;border:1px solid #ddd} code{background:#f3f4f6;padding:1px 4px}"
    html_parts = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        "<title>gdTAI v4 T-cell gate report</title>",
        f"<style>{css}</style></head><body>",
        "<h1>gdTAI v4 Two-Step T-Cell Gate Classifier</h1>",
        "<p>Stage 1 is a soft T-lineage gate; Stage 2 is a dropout-tolerant gdT classifier. No H5AD was modified.</p>",
        fig_html,
        "<h2>Stage 1 T-Cell Gate Tune Metrics</h2>", dataframe_to_html(stage1),
        "<h2>Stage 2 gdT Tune Metrics</h2>", dataframe_to_html(tune),
        "<h2>Internal Validation and Sensitivity</h2>", dataframe_to_html(validation.head(100) if not validation.empty else validation),
        "<h2>Full-Atlas v4 Summary</h2>", dataframe_to_html(full),
        "<h2>Comparison Against Existing Strategies</h2>", dataframe_to_html(comparison),
        f"<p>Tables: <code>{html.escape(str(TABLE_DIR))}</code><br>Model: <code>{html.escape(str(MODEL_PKL))}</code></p>",
        "</body></html>",
    ]
    REPORT_HTML.write_text("\n".join(html_parts), encoding="utf-8")
    if not args.no_pdf:
        try:
            subprocess.run(["google-chrome", "--headless", "--disable-gpu", "--no-sandbox", f"--print-to-pdf={REPORT_PDF}", str(REPORT_HTML)], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        except Exception as exc:
            logging.warning("PDF export failed: %s", exc)


def package_model(model: GdTAIV4TwoStageModel, spec: FeatureSpec, stage1: pd.DataFrame, tune: pd.DataFrame, validation: pd.DataFrame, full: pd.DataFrame, args: argparse.Namespace) -> None:
    model.gene_names = list(spec.gene_names)
    model.engineered_feature_names = list(spec.engineered_feature_names)
    model.model_feature_names = list(spec.model_feature_names)
    payload = {
        "model_family": "gdTAI",
        "version": "gdTAI_v4.0",
        "wrapper_name": "two_stage_tcell_gate_dropout_tolerant_gdt_classifier",
        "normalization": "log1p(counts per 10000)",
        "tcell_model": model.tcell_model,
        "gdt_model": model.gdt_model,
        "gene_names": spec.gene_names,
        "engineered_feature_names": spec.engineered_feature_names,
        "feature_names": spec.model_feature_names,
        "tcell_threshold": model.tcell_threshold,
        "operating_modes": model.mode_thresholds,
        "notes": "Two-stage soft T-lineage gate plus gdT classifier; no hard TRDV requirement; NK/cytotoxic genes are contextual features.",
    }
    with MODEL_PKL.open("wb") as handle:
        pickle.dump(payload, handle)
    pd.DataFrame({"feature": spec.model_feature_names, "feature_index": np.arange(len(spec.model_feature_names), dtype=int)}).to_csv(MODEL_DIR / "feature_genes.csv", index=False)
    mode_metrics = tune.copy()
    if not validation.empty:
        val_primary = validation[validation["evaluation_set"].eq("validation_GDT2020_cordblood_vs_GSE254249_TCRAB")]
        mode_metrics = mode_metrics.merge(val_primary[["strategy", "precision", "recall", "specificity", "f1", "fp", "tp"]], on="strategy", how="left", suffixes=("_tune", "_validation"))
    mode_metrics.to_csv(MODEL_DIR / "mode_metrics.csv", index=False)
    manifest = {
        "version": "gdTAI_v4.0",
        "model": "two_stage_tcell_gate_dropout_tolerant_gdt_classifier",
        "normalization": "log1p(counts per 10000)",
        "n_gene_features": len(spec.gene_names),
        "n_engineered_features": len(spec.engineered_feature_names),
        "n_total_features_stage1": len(spec.model_feature_names),
        "n_total_features_stage2": len(spec.model_feature_names) + 1,
        "tcell_threshold": model.tcell_threshold,
        "operating_modes": model.mode_thresholds,
        "input_h5ad": str(args.input_h5ad),
        "paths": {
            "model": str(MODEL_PKL),
            "feature_table": str(MODEL_DIR / "feature_genes.csv"),
            "tables": str(TABLE_DIR),
            "report_html": str(REPORT_HTML),
            "report_pdf": str(REPORT_PDF) if REPORT_PDF.exists() else None,
        },
    }
    MANIFEST_JSON.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (MODEL_DIR / "README.md").write_text(
        "# gdTAI v4.0\n\nTwo-stage soft T-lineage gate plus dropout-tolerant gdT classifier. Use log1p(counts per 10,000) features through the packaged protocol; no hard TRDV requirement is applied.\n",
        encoding="utf-8",
    )


def main() -> None:
    args = parse_args()
    setup_logging()
    logging.info("Opening %s", args.input_h5ad)
    with h5py.File(args.input_h5ad, "r") as handle:
        n_obs, n_vars = h5ad_shape(handle)
        logging.info("Input shape: %s cells x %s genes", f"{n_obs:,}", f"{n_vars:,}")
        obs = build_obs_metadata(handle)
        annotation = annotation_vector(handle, n_obs)
        spec = make_feature_spec(handle)
        all_rows, x_all, splits, pool_summary, split_summary = build_training_matrices(handle, spec, obs, annotation, args)
        model, stage1_metrics, tune_metrics, _sweep = train_models_direct(all_rows, x_all, splits, obs, annotation, spec, args)
        validation = evaluate_internal(all_rows, x_all, splits, model, obs, annotation, spec)
        full_outputs = {}
        if not args.no_full_atlas:
            full_outputs = apply_full_atlas(handle, spec, model, obs, annotation, args)
        full_overall = full_outputs.get("overall", pd.DataFrame())
    comparison = compare_against_existing(full_overall) if not full_overall.empty else pd.DataFrame()
    figures = write_figures(validation, comparison)
    package_model(model, spec, stage1_metrics, tune_metrics, validation, full_overall, args)
    render_report(stage1_metrics, tune_metrics, validation, full_overall, comparison, figures, args)
    summary = {
        "model": "gdTAI_v4.0",
        "tcell_threshold": model.tcell_threshold,
        "operating_modes": model.mode_thresholds,
        "paths": {
            "tables": str(TABLE_DIR),
            "figures": str(FIGURE_DIR),
            "model": str(MODEL_PKL),
            "report_html": str(REPORT_HTML),
            "report_pdf": str(REPORT_PDF) if REPORT_PDF.exists() else None,
        },
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    logging.info("Done. Report: %s", REPORT_HTML)


if __name__ == "__main__":
    main()
