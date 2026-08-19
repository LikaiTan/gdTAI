#!/usr/bin/env python3
"""Diagnose the frozen gdTAI V4.3 recall failure without retuning it."""

from __future__ import annotations

import hashlib
import html
import json
import math
import os
import pickle
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-gdtai-v43-recall")

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xgboost as xgb
from scipy.stats import binomtest, fisher_exact, mannwhitneyu
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_confint


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "workflows/gdtai") not in sys.path:
    sys.path.insert(0, str(ROOT / "workflows/gdtai"))

from run_gdtai_v4_2_ground_truth import axis_node_values  # noqa: E402


ATLAS = Path("/ssd/tnk_phase3/Integrated_dataset/full_atlas/tcr_corrected/integrated_full_atlas.h5ad")
PREDICTIONS = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_3_common_lockbox/common_lockbox_predictions.parquet"
LOCKBOX_MATRIX = ROOT / "Integrated_dataset/cache/gdT_prediction/gdtai_v4_3_common_lockbox/lockbox_gene_features.npy"
LOCKBOX_FEATURES = ROOT / "Integrated_dataset/cache/gdT_prediction/gdtai_v4_3_common_lockbox/lockbox_feature_manifest.csv"
LABELS = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_ground_truth/v4_2_label_manifest.parquet"
DEVELOPMENT_MATRIX = ROOT / "Integrated_dataset/cache/gdT_prediction/gdtai_v4_2_training/gene_features.npy"
DEVELOPMENT_FEATURES = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_training/feature_manifest.csv"
OOF = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_3_receptor_isolated/final_development/grouped_oof_predictions.parquet"
V4_DIR = ROOT / "Integrated_dataset/models/gdT_prediction/gdtai_v4_3_receptor_isolated/final_development"
V4_CONTRACT = V4_DIR / "model_contract.json"
V2_MODEL = ROOT / "Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v2.0/gdTAI_v2_model.pkl"
V3_MODEL = ROOT / "Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/gdTAI_v3_model.pkl"

TABLE_DIR = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_3_recall_failure"
FIGURE_DIR = ROOT / "Integrated_dataset/figures/gdT_prediction/gdtai_v4_3_recall_failure"
LOG_DIR = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_3_recall_failure"
REPORT_DIR = ROOT / "gdT_prediction/gdtai_v4_3_recall_failure"
ASSET_DIR = REPORT_DIR / "assets"

STAGE1_THRESHOLD = 0.01939413930755109
FROZEN_THRESHOLD = 0.9940980704625443
COLOR_V4 = "#cc5a3d"
COLOR_V3 = "#2878b5"
COLOR_V2 = "#2f8f6b"
COLOR_FN = "#d97706"
COLOR_TP = "#2f8f6b"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_pickle(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        payload = pickle.load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected dictionary payload: {path}")
    return payload


def decode_subset(node: h5py.Group | h5py.Dataset, rows: np.ndarray) -> np.ndarray:
    if isinstance(node, h5py.Group) and {"categories", "codes"}.issubset(node):
        categories = axis_node_values(node["categories"])
        codes = np.asarray(node["codes"][rows], dtype=np.int64)
        output = np.full(len(rows), "", dtype=object)
        present = codes >= 0
        output[present] = categories[codes[present]]
        return output
    values = np.asarray(node[rows])
    if values.dtype.kind in {"O", "S", "U"}:
        return np.asarray([
            value.decode("utf-8", errors="replace") if isinstance(value, (bytes, np.bytes_)) else str(value)
            for value in values
        ], dtype=object)
    return values


def bh(values: Iterable[float]) -> np.ndarray:
    values = np.asarray(list(values), dtype=float)
    return multipletests(values, method="fdr_bh")[1] if len(values) else values


def ci(successes: int, total: int) -> tuple[float, float]:
    if total == 0:
        return math.nan, math.nan
    low, high = proportion_confint(successes, total, method="wilson")
    return float(low), float(high)


def full_v4_call(frame: pd.DataFrame, threshold: float) -> np.ndarray:
    return (
        (frame.v4_3_stage1_score.to_numpy(float) >= STAGE1_THRESHOLD)
        & (
            (frame.v4_3_stage2_score.to_numpy(float) >= threshold)
            | frame.v4_3_receptor_rescue.to_numpy(bool)
        )
        & ~frame.v4_3_cd4_treg_excluded.to_numpy(bool)
    )


def full_oof_call(frame: pd.DataFrame, threshold: float) -> np.ndarray:
    return (
        (frame.stage1_score.to_numpy(float) >= STAGE1_THRESHOLD)
        & ((frame.stage2_score.to_numpy(float) >= threshold) | frame.receptor_rescue.to_numpy(bool))
        & ~frame.cd4_or_treg_excluded.to_numpy(bool)
    )


def binary_metrics(frame: pd.DataFrame, calls: np.ndarray) -> dict[str, float | int]:
    positive = frame.truth_class.eq("gdT_gold").to_numpy()
    ab = frame.truth_class.eq("abT_gold").to_numpy()
    nk = frame.truth_class.astype(str).str.startswith("nk_").to_numpy()
    negative = ab | nk
    tp = int((calls & positive).sum())
    fn = int((~calls & positive).sum())
    fp = int((calls & negative).sum())
    tn = int((~calls & negative).sum())
    recall = tp / (tp + fn) if tp + fn else math.nan
    precision = tp / (tp + fp) if tp + fp else math.nan
    return {
        "tp": tp,
        "fn": fn,
        "fp": fp,
        "tn": tn,
        "recall": recall,
        "precision": precision,
        "f1": 2 * precision * recall / (precision + recall) if precision + recall else math.nan,
        "abt_fpr": float(calls[ab].mean()) if ab.any() else math.nan,
        "nk_fpr": float(calls[nk].mean()) if nk.any() else math.nan,
    }


def diagnostic_best_threshold(frame: pd.DataFrame) -> tuple[float, dict[str, float | int]]:
    eligible = (
        (frame.v4_3_stage1_score.to_numpy(float) >= STAGE1_THRESHOLD)
        & ~frame.v4_3_cd4_treg_excluded.to_numpy(bool)
    )
    score = np.where(eligible, frame.v4_3_stage2_score.to_numpy(float), -1.0)
    score[eligible & frame.v4_3_receptor_rescue.to_numpy(bool)] = 2.0
    y = frame.truth_class.eq("gdT_gold").to_numpy(np.int8)
    precision, recall, thresholds = precision_recall_curve(y, score)
    f1 = np.divide(2 * precision * recall, precision + recall, out=np.zeros_like(precision), where=(precision + recall) > 0)
    index = int(np.nanargmax(f1[:-1]))
    threshold = float(thresholds[index])
    return threshold, binary_metrics(frame, full_v4_call(frame, threshold))


def receptor_states(matrix: np.ndarray, genes: list[str]) -> pd.DataFrame:
    lookup = {gene: index for index, gene in enumerate(genes)}
    trdv = [lookup[gene] for gene in genes if gene.startswith("TRDV")]
    trgv = [lookup[gene] for gene in genes if gene.startswith("TRGV")]
    alpha = [lookup[gene] for gene in genes if gene.startswith(("TRAV", "TRAJ"))]
    beta = [lookup[gene] for gene in genes if gene.startswith(("TRBV", "TRBJ"))]
    trdc = matrix[:, lookup["TRDC"]] > 0
    any_trdv = (matrix[:, trdv] > 0).any(axis=1)
    any_trgv = (matrix[:, trgv] > 0).any(axis=1)
    trdv1 = matrix[:, lookup["TRDV1"]] > 0
    trdv2 = matrix[:, lookup["TRDV2"]] > 0
    output = pd.DataFrame({
        "TRDC_detected": trdc,
        "any_TRDV_detected": any_trdv,
        "any_TRGV_detected": any_trgv,
        "any_alpha_VJ_detected": (matrix[:, alpha] > 0).any(axis=1),
        "any_beta_VJ_detected": (matrix[:, beta] > 0).any(axis=1),
        "n_TRD_genes_detected": (matrix[:, [lookup[g] for g in genes if g.startswith("TRD")]] > 0).sum(axis=1),
        "n_TRG_genes_detected": (matrix[:, [lookup[g] for g in genes if g.startswith("TRG")]] > 0).sum(axis=1),
        "n_AB_VJ_genes_detected": (matrix[:, alpha + beta] > 0).sum(axis=1),
    })
    output["TRDC_TRDV_state"] = np.select(
        [trdc & any_trdv, trdc & ~any_trdv, ~trdc & any_trdv, ~trdc & ~any_trdv],
        ["TRDC+ TRDV+", "TRDC+ TRDV-", "TRDC- TRDV+", "TRDC- TRDV-"],
        default="unknown",
    )
    output["TRDV1_TRDV2_state"] = np.select(
        [trdv1 & trdv2, trdv1 & ~trdv2, ~trdv1 & trdv2, ~trdv1 & ~trdv2],
        ["TRDV1+ TRDV2+", "TRDV1+ TRDV2-", "TRDV1- TRDV2+", "TRDV1- TRDV2-"],
        default="unknown",
    )
    return output


def group_recall(frame: pd.DataFrame, column: str) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for value, group in frame.groupby(column, dropna=False, sort=False):
        for model, call_column in [
            ("V4.3 frozen", "v4_3_balanced"),
            ("V3 balanced", "v3_balanced"),
            ("V2 high-F1", "v2_high_f1"),
            ("V2 high-purity", "v2_high_purity"),
        ]:
            successes = int(group[call_column].sum())
            low, high = ci(successes, len(group))
            rows.append({
                column: value,
                "model": model,
                "n_cells": len(group),
                "detected": successes,
                "recall": successes / len(group),
                "ci_low": low,
                "ci_high": high,
            })
    return pd.DataFrame(rows)


def score_quantiles(frame: pd.DataFrame, role: str) -> pd.DataFrame:
    rows = []
    score_columns = {
        "V4.3 Stage 2": "v4_3_stage2_score",
        "V3 balanced": "v3_balanced_score",
        "V2 base": "v2_score",
    }
    for source, group in frame.groupby("source_gse_id", sort=True):
        for model, column in score_columns.items():
            values = group[column].to_numpy(float)
            values = values[np.isfinite(values)]
            if not len(values):
                continue
            quantiles = np.quantile(values, [0.1, 0.25, 0.5, 0.75, 0.9])
            rows.append({
                "role": role,
                "source_gse_id": source,
                "model": model,
                "n_cells": len(group),
                "q10": quantiles[0],
                "q25": quantiles[1],
                "median": quantiles[2],
                "q75": quantiles[3],
                "q90": quantiles[4],
            })
    return pd.DataFrame(rows)


def final_model_on_oof(oof: pd.DataFrame, contract: dict[str, Any]) -> tuple[np.ndarray, pd.DataFrame]:
    labels = pd.read_parquet(LABELS, columns=["cell_id"])
    positions = pd.Index(labels.cell_id).get_indexer(oof.cell_id)
    if (positions < 0).any():
        raise RuntimeError("OOF rows do not map to the frozen feature cache")
    features = pd.read_csv(DEVELOPMENT_FEATURES).sort_values("feature_index")
    lookup = {gene: index for index, gene in enumerate(features.gene.astype(str))}
    columns = np.asarray([lookup[gene] for gene in contract["stage2_feature_names"]], dtype=np.int64)
    matrix = np.load(DEVELOPMENT_MATRIX, mmap_mode="r")
    booster = xgb.Booster()
    booster.load_model(V4_DIR / "stage2_receptor_classifier.ubj")
    booster.set_param({"device": "cpu", "nthread": 16})
    output = []
    for start in range(0, len(positions), 100_000):
        local = np.asarray(matrix[positions[start:start + 100_000]][:, columns], dtype=np.float32)
        output.append(booster.predict(xgb.DMatrix(local)))
    score = np.concatenate(output).astype(np.float32)
    rows = []
    for truth, group_index in oof.groupby("truth_class", sort=False).groups.items():
        index = np.asarray(group_index, dtype=np.int64)
        old = oof.stage2_score.to_numpy(float)[index]
        new = score[index]
        rows.append({
            "truth_class": truth,
            "n_cells": len(index),
            "oof_median": float(np.median(old)),
            "final_model_median": float(np.median(new)),
            "pearson_correlation": float(np.corrcoef(old, new)[0, 1]),
            "oof_above_frozen_threshold": float((old >= FROZEN_THRESHOLD).mean()),
            "final_model_above_frozen_threshold": float((new >= FROZEN_THRESHOLD).mean()),
            "final_model_tree_count": booster.num_boosted_rounds(),
        })
    return score, pd.DataFrame(rows)


def feature_statistics(frame: pd.DataFrame, matrix: np.ndarray, genes: list[str]) -> pd.DataFrame:
    selected = [
        "TRDC", "TRDV1", "TRDV2", "TRDV3", "TRGC1", "TRGC2", "TRGV3", "TRGV4",
        "TRGV5", "TRGV9", "TRBC1", "TRBC2", "CD3D", "CD3E", "CD3G",
    ]
    lookup = {gene: index for index, gene in enumerate(genes)}
    tp = frame.v4_3_balanced.to_numpy(bool)
    fn = ~tp
    rows = []
    for gene in selected:
        values = matrix[:, lookup[gene]].astype(float)
        a = values[tp]
        b = values[fn]
        table = [[int((a > 0).sum()), int((a == 0).sum())], [int((b > 0).sum()), int((b == 0).sum())]]
        fisher_p = fisher_exact(table)[1]
        mw = mannwhitneyu(a, b, alternative="two-sided")
        rank_biserial = 2 * float(mw.statistic) / (len(a) * len(b)) - 1
        rows.append({
            "gene": gene,
            "v4_tp_detected_fraction": float((a > 0).mean()),
            "v4_fn_detected_fraction": float((b > 0).mean()),
            "detected_fraction_difference_tp_minus_fn": float((a > 0).mean() - (b > 0).mean()),
            "v4_tp_median_expression": float(np.median(a)),
            "v4_fn_median_expression": float(np.median(b)),
            "rank_biserial_tp_minus_fn": rank_biserial,
            "fisher_p": float(fisher_p),
            "mannwhitney_p": float(mw.pvalue),
        })
    table = pd.DataFrame(rows)
    table["fisher_fdr"] = bh(table.fisher_p)
    table["mannwhitney_fdr"] = bh(table.mannwhitney_p)
    return table


def qc_statistics(frame: pd.DataFrame) -> pd.DataFrame:
    tp = frame.v4_3_balanced.to_numpy(bool)
    rows = []
    for column in ["n_genes_by_counts", "total_counts", "pct_counts_mt", "TRD_umis", "TRG_umis"]:
        values = pd.to_numeric(frame[column], errors="coerce").to_numpy(float)
        a = values[tp & np.isfinite(values)]
        b = values[~tp & np.isfinite(values)]
        test = mannwhitneyu(a, b, alternative="two-sided")
        rank_biserial = 2 * float(test.statistic) / (len(a) * len(b)) - 1
        rows.append({
            "variable": column,
            "v4_tp_n": len(a),
            "v4_fn_n": len(b),
            "v4_tp_median": float(np.median(a)),
            "v4_fn_median": float(np.median(b)),
            "rank_biserial_tp_minus_fn": rank_biserial,
            "mannwhitney_p": float(test.pvalue),
        })
    table = pd.DataFrame(rows)
    table["mannwhitney_fdr"] = bh(table.mannwhitney_p)
    return table


def stage2_shap(frame: pd.DataFrame, matrix: np.ndarray, genes: list[str], contract: dict[str, Any]) -> pd.DataFrame:
    lookup = {gene: index for index, gene in enumerate(genes)}
    feature_names = contract["stage2_feature_names"]
    columns = [lookup[gene] for gene in feature_names]
    local = np.asarray(matrix[:, columns], dtype=np.float32)
    booster = xgb.Booster()
    booster.load_model(V4_DIR / "stage2_receptor_classifier.ubj")
    contributions = booster.predict(xgb.DMatrix(local), pred_contribs=True)
    reconstructed = 1 / (1 + np.exp(-contributions.sum(axis=1)))
    if np.max(np.abs(reconstructed - frame.v4_3_stage2_score.to_numpy(float))) > 1e-5:
        raise RuntimeError("SHAP contributions do not reconstruct frozen V4.3 scores")
    rows = []
    groups = {
        "V4_TP": frame.v4_3_balanced.to_numpy(bool),
        "V4_FN_V2_TP": (~frame.v4_3_balanced & frame.v2_high_f1).to_numpy(bool),
        "V4_FN_V2_FN": (~frame.v4_3_balanced & ~frame.v2_high_f1).to_numpy(bool),
    }
    for group_name, mask in groups.items():
        for index, gene in enumerate(feature_names):
            rows.append({
                "failure_group": group_name,
                "gene": gene,
                "n_cells": int(mask.sum()),
                "mean_shap": float(contributions[mask, index].mean()),
                "median_shap": float(np.median(contributions[mask, index])),
                "detected_fraction": float((local[mask, index] > 0).mean()),
                "mean_expression": float(local[mask, index].mean()),
            })
    return pd.DataFrame(rows)


def v2_coefficients(v2: dict[str, Any], v4_genes: set[str]) -> pd.DataFrame:
    pipe = v2["base_model"]["model_object"]
    scaler = pipe.named_steps["standardscaler"]
    model = pipe.named_steps["logisticregression"]
    genes = np.asarray(v2["base_model"]["gene_names"], dtype=object)
    return pd.DataFrame({
        "gene": genes,
        "standardized_coefficient": model.coef_[0],
        "coefficient_per_log1p_cp10k": model.coef_[0] / scaler.scale_,
        "training_mean": scaler.mean_,
        "training_sd": scaler.scale_,
        "present_in_v4_stage2": [gene in v4_genes for gene in genes],
    })


def make_figures(
    balf: pd.DataFrame,
    balf_matrix: np.ndarray,
    genes: list[str],
    receptor_recall: pd.DataFrame,
    feature_stats: pd.DataFrame,
    shap_table: pd.DataFrame,
    oof: pd.DataFrame,
    per_source_ab: pd.DataFrame,
    diagnostic_threshold: float,
) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10, "figure.dpi": 140})

    positive = balf.truth_class.eq("gdT_gold")
    ab = balf.truth_class.eq("abT_gold")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    for axis, column, title, threshold, color in [
        (axes[0], "v4_3_stage2_score", "V4.3 Stage-2 score", FROZEN_THRESHOLD, COLOR_V4),
        (axes[1], "v2_score", "V2 linear score", float(balf.v2_high_f1_threshold.iloc[0]), COLOR_V2),
    ]:
        for mask, label, local_color in [(positive, "gdT gold", color), (ab, "paired alpha-beta", "#6b7280")]:
            values = np.clip(balf.loc[mask, column].to_numpy(float), 0, 1 - 1e-7)
            transformed = -np.log10(1 - values + 1e-7)
            axis.hist(transformed, bins=70, density=True, histtype="step", linewidth=1.7, color=local_color, label=label)
        axis.axvline(-math.log10(1 - threshold + 1e-7), color="#111827", linestyle="--", linewidth=1.2, label="frozen cutoff")
        axis.set(title=title, xlabel="-log10(1 - model score)", ylabel="Density")
        axis.legend(frameon=False)
    fig.suptitle("BALF_BLOOD_COPD: score calibration differs despite strong ranking")
    save_figure(fig, "balf_score_distributions.png")

    states = ["TRDC+ TRDV+", "TRDC+ TRDV-", "TRDC- TRDV+", "TRDC- TRDV-"]
    models = ["V4.3 frozen", "V3 balanced", "V2 high-F1"]
    colors = [COLOR_V4, COLOR_V3, COLOR_V2]
    fig, axis = plt.subplots(figsize=(9.5, 4.5), constrained_layout=True)
    x = np.arange(len(states))
    width = 0.24
    for offset, model, color in zip([-width, 0, width], models, colors):
        subset = receptor_recall[receptor_recall.model.eq(model)].set_index("TRDC_TRDV_state").reindex(states)
        axis.bar(x + offset, subset.recall, width, color=color, label=model)
    counts = receptor_recall[receptor_recall.model.eq("V4.3 frozen")].set_index("TRDC_TRDV_state").reindex(states).n_cells
    axis.set_xticks(x, [f"{state}\n(n={int(count)})" for state, count in zip(states, counts)])
    axis.set_ylim(0, 1.05)
    axis.set_ylabel("Recall")
    axis.set_title("BALF gold recall by transcriptomic TRDC/TRDV detection state")
    axis.legend(frameon=False, ncol=3, loc="upper right")
    save_figure(fig, "balf_recall_by_receptor_state.png")

    fig, axis = plt.subplots(figsize=(6.2, 5.4), constrained_layout=True)
    groups = np.where(balf.v4_3_balanced, "V4 TP", np.where(balf.v2_high_f1, "V4 FN / V2 TP", "V4 FN / V2 FN"))
    palette = {"V4 TP": COLOR_TP, "V4 FN / V2 TP": COLOR_FN, "V4 FN / V2 FN": "#7c3aed"}
    for group in ["V4 FN / V2 FN", "V4 FN / V2 TP", "V4 TP"]:
        mask = groups == group
        axis.scatter(
            balf.loc[mask, "v2_score"], balf.loc[mask, "v4_3_stage2_score"],
            s=13, alpha=0.55, edgecolors="none", color=palette[group], label=f"{group} (n={mask.sum()})",
        )
    axis.axhline(FROZEN_THRESHOLD, color="#111827", linestyle="--", linewidth=1)
    axis.axvline(float(balf.v2_high_f1_threshold.iloc[0]), color="#111827", linestyle=":", linewidth=1)
    axis.set(xlabel="V2 base score", ylabel="V4.3 Stage-2 score", title="Most V4.3 BALF misses remain high-scoring V2 positives", xlim=(-0.02, 1.02), ylim=(-0.02, 1.02))
    axis.legend(frameon=False, fontsize=8, loc="lower left")
    save_figure(fig, "balf_v4_vs_v2_scores.png")

    top = feature_stats.reindex(feature_stats.detected_fraction_difference_tp_minus_fn.abs().sort_values(ascending=False).index).head(12)
    fig, axis = plt.subplots(figsize=(8.4, 5.0), constrained_layout=True)
    y = np.arange(len(top))
    axis.barh(y + 0.18, top.v4_tp_detected_fraction, 0.36, color=COLOR_TP, label="V4 TP")
    axis.barh(y - 0.18, top.v4_fn_detected_fraction, 0.36, color=COLOR_FN, label="V4 FN")
    axis.set_yticks(y, top.gene)
    axis.invert_yaxis()
    axis.set(xlabel="Detected fraction", title="Receptor-expression differences among BALF gold cells", xlim=(0, 1))
    axis.legend(frameon=False)
    save_figure(fig, "balf_feature_detection_tp_fn.png")

    pivot = shap_table.pivot(index="gene", columns="failure_group", values="mean_shap").fillna(0)
    pivot["TP_minus_FN"] = pivot["V4_TP"] - pivot["V4_FN_V2_TP"]
    selected = pivot.reindex(pivot.TP_minus_FN.abs().sort_values(ascending=False).head(15).index).sort_values("TP_minus_FN")
    fig, axis = plt.subplots(figsize=(8.7, 5.2), constrained_layout=True)
    axis.barh(selected.index, selected.TP_minus_FN, color=np.where(selected.TP_minus_FN >= 0, COLOR_TP, COLOR_FN))
    axis.axvline(0, color="#111827", linewidth=0.8)
    axis.set(xlabel="Mean SHAP difference: V4 TP - V4 FN/V2 TP", title="V4.3 decision drivers separating recovered and missed BALF gold cells")
    save_figure(fig, "balf_shap_tp_fn_difference.png")

    thresholds = np.unique(np.r_[np.linspace(0.05, 0.95, 45), np.linspace(0.955, 0.9995, 70), diagnostic_threshold, FROZEN_THRESHOLD])
    records = []
    for threshold in thresholds:
        for name, frame, caller in [("Development OOF", oof, full_oof_call), ("BALF external", balf, full_v4_call)]:
            calls = caller(frame, float(threshold))
            positive_mask = frame.truth_class.eq("gdT_gold").to_numpy()
            ab_mask = frame.truth_class.eq("abT_gold").to_numpy()
            records.append({
                "threshold": threshold,
                "x": -math.log10(1 - threshold + 1e-7),
                "cohort": name,
                "recall": float(calls[positive_mask].mean()),
                "abt_fpr": float(calls[ab_mask].mean()),
            })
    frontier = pd.DataFrame(records)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)
    for cohort, group in frontier.groupby("cohort"):
        color = COLOR_V4 if cohort == "BALF external" else COLOR_V3
        axes[0].plot(group.x, group.recall, color=color, linewidth=1.8, label=cohort)
        axes[1].plot(group.x, group.abt_fpr * 100, color=color, linewidth=1.8, label=cohort)
    frozen_x = -math.log10(1 - FROZEN_THRESHOLD + 1e-7)
    for axis in axes:
        axis.axvline(frozen_x, color="#111827", linestyle="--", linewidth=1.1)
        ticks = [0.5, 0.8, 0.9, 0.95, 0.99, FROZEN_THRESHOLD, 0.999]
        axis.set_xticks([-math.log10(1 - value + 1e-7) for value in ticks], [f"{value:.3f}" for value in ticks], rotation=35)
        axis.set_xlabel("V4.3 Stage-2 threshold")
    axes[0].set(ylabel="gdT-gold recall", title="Recall transfer")
    axes[1].set(ylabel="Paired-alpha-beta FPR (%)", title="Specificity transfer")
    axes[1].axhline(0.2, color="#9ca3af", linestyle=":", linewidth=1, label="development limit")
    axes[0].legend(frameon=False)
    axes[1].legend(frameon=False)
    fig.suptitle("The frozen threshold sits in a source-sensitive saturated score region")
    save_figure(fig, "threshold_transfer.png")

    lookup = {gene: index for index, gene in enumerate(genes)}
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 4.0), constrained_layout=True)
    group_masks = [balf.v4_3_balanced.to_numpy(bool), (~balf.v4_3_balanced & balf.v2_high_f1).to_numpy(bool)]
    labels = ["V4 TP", "V4 FN / V2 TP"]
    for axis, (column, title) in zip(axes, [("TRDC", "TRDC"), ("TRDV2", "TRDV2"), ("TRGV9", "TRGV9")]):
        values = [balf_matrix[mask, lookup[column]] for mask in group_masks]
        axis.boxplot(values, tick_labels=labels, showfliers=False, patch_artist=True, boxprops={"facecolor": "#dbeafe"})
        axis.set(title=title, ylabel="log1p CP10K")
        axis.tick_params(axis="x", rotation=20)
    fig.suptitle("BALF V4 misses retain receptor evidence but have weaker TRDC / asymmetric V-gene support")
    save_figure(fig, "balf_receptor_expression_boxes.png")

    local = per_source_ab.sort_values("v4_minus_v3_fpr", ascending=True)
    fig, axis = plt.subplots(figsize=(8.8, 5.1), constrained_layout=True)
    colors = np.where(local.v4_minus_v3_fpr > 0, COLOR_V4, COLOR_TP)
    axis.barh(local.source_gse_id, local.v4_minus_v3_fpr * 100, color=colors)
    axis.axvline(0, color="#111827", linewidth=0.9)
    axis.set(
        xlabel="Paired-alpha-beta FPR difference, V4.3 - V3 (percentage points)",
        title="V4.3 has higher paired-alpha-beta FPR in most new datasets",
    )
    for index, row in enumerate(local.itertuples(index=False)):
        difference = row.v4_minus_v3_fpr * 100
        if difference < 0:
            axis.text(difference / 2, index, f"{difference:+.3f}", va="center", ha="center", fontsize=8, color="white")
        else:
            axis.text(difference + 0.012, index, f"{difference:+.3f}", va="center", ha="left", fontsize=8)
    save_figure(fig, "per_source_ab_fpr_difference.png")


def save_figure(fig: plt.Figure, name: str) -> None:
    canonical = FIGURE_DIR / name
    report = ASSET_DIR / name
    fig.savefig(canonical, dpi=220, bbox_inches="tight", facecolor="white")
    fig.savefig(report, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fmt(value: Any, digits: int = 3) -> str:
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "NA"
    if isinstance(value, (int, np.integer)):
        return f"{int(value):,}"
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.{digits}f}"
    return str(value)


def dataframe_html(frame: pd.DataFrame, columns: list[str] | None = None, percent: set[str] | None = None, limit: int | None = None) -> str:
    local = frame.copy()
    if columns is not None:
        local = local[columns]
    if limit is not None:
        local = local.head(limit)
    percent = percent or set()
    for column in local.columns:
        if column in percent:
            local[column] = local[column].map(lambda value: "NA" if pd.isna(value) else f"{100 * float(value):.2f}%")
        elif pd.api.types.is_float_dtype(local[column]):
            local[column] = local[column].map(lambda value: fmt(value, 4))
        elif pd.api.types.is_integer_dtype(local[column]):
            local[column] = local[column].map(lambda value: f"{int(value):,}")
    return local.to_html(index=False, border=0, escape=True, classes="data-table")


def render_report(
    exposure: pd.DataFrame,
    headline: pd.DataFrame,
    operating: pd.DataFrame,
    receptor_recall: pd.DataFrame,
    trdv_recall: pd.DataFrame,
    failure_summary: pd.DataFrame,
    feature_stats: pd.DataFrame,
    qc_stats: pd.DataFrame,
    threshold_drivers: pd.DataFrame,
    score_shift: pd.DataFrame,
    contract_comparison: pd.DataFrame,
    per_source_ab: pd.DataFrame,
    diagnostic_threshold: float,
    summary: dict[str, Any],
) -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    key_receptor = receptor_recall[receptor_recall.model.isin(["V4.3 frozen", "V3 balanced", "V2 high-F1"])]
    top_features = feature_stats.sort_values("detected_fraction_difference_tp_minus_fn", ascending=False).head(10)
    driver = threshold_drivers[
        threshold_drivers.threshold_label.isin(["BALF diagnostic best-F1", "0.8006 diagnostic", "Frozen 0.9941"])
    ]
    style = """
    @page { size: A4 landscape; margin: 12mm; }
    * { box-sizing: border-box; }
    body { margin: 0; color: #172033; background: #f5f7fa; font-family: Arial, Helvetica, sans-serif; font-size: 10pt; line-height: 1.42; }
    main { max-width: 1240px; margin: 0 auto; background: white; }
    header { padding: 28px 34px 24px; background: #172033; color: white; border-bottom: 7px solid #cc5a3d; }
    h1 { margin: 0 0 8px; font-size: 28px; letter-spacing: 0; }
    header p { margin: 0; color: #dbe3ee; max-width: 980px; }
    section { padding: 24px 34px; border-bottom: 1px solid #dfe5ec; break-inside: avoid; }
    h2 { margin: 0 0 12px; font-size: 19px; color: #172033; }
    h3 { margin: 14px 0 7px; font-size: 13px; }
    p { margin: 7px 0; }
    .cards { display: grid; grid-template-columns: repeat(4, 1fr); gap: 10px; margin: 14px 0; }
    .card { border: 1px solid #d7dee7; border-top: 4px solid #2878b5; padding: 11px; background: #fff; }
    .card strong { display: block; font-size: 20px; margin: 3px 0; }
    .callout { padding: 12px 14px; border-left: 5px solid #cc5a3d; background: #fff4ee; margin: 12px 0; }
    .good { border-left-color: #2f8f6b; background: #effaf5; }
    .warning { border-left-color: #d97706; background: #fff8e8; }
    .grid2 { display: grid; grid-template-columns: 1fr 1fr; gap: 18px; align-items: start; }
    img { display: block; width: 100%; height: auto; }
    figure { margin: 8px 0; break-inside: avoid; }
    figcaption { color: #526071; font-size: 8.8pt; margin-top: 5px; }
    .data-table { width: 100%; border-collapse: collapse; font-size: 8.4pt; margin: 9px 0; table-layout: auto; break-inside: avoid; }
    .data-table th { background: #e9eef4; text-align: left; padding: 5px 6px; border: 1px solid #cfd8e3; }
    .data-table td { padding: 4px 6px; border: 1px solid #dbe2ea; vertical-align: top; }
    code { font-family: Consolas, monospace; background: #eef2f6; padding: 1px 3px; }
    ul { margin: 7px 0 7px 20px; padding: 0; }
    li { margin: 4px 0; }
    .small { font-size: 8.5pt; color: #526071; }
    .page-break { break-before: page; }
    @media (max-width: 800px) { .cards, .grid2 { grid-template-columns: 1fr; } section, header { padding-left: 18px; padding-right: 18px; } }
    """
    body = f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>gdTAI V4.3 recall failure investigation</title><style>{style}</style></head><body><main>
    <header><h1>Why gdTAI V4.3 lost recall</h1><p>Frozen, read-only investigation of BALF_BLOOD_COPD and the common external lockbox. No model was fitted, no cutoff was replaced, and the consumed lockbox remains unavailable for future tuning.</p></header>
    <section>
      <h2>Conclusion</h2>
      <div class="callout"><strong>The primary failure is threshold transfer, not inability to rank BALF gdT cells.</strong> V4.3 separates BALF receptor-gold positives from paired-alpha-beta controls extremely well, but its XGBoost probabilities are saturated and not calibrated across cohorts. The development safety constraints forced a cutoff of <code>0.994098</code>, just above the BALF positive median of <code>{summary['balf_v4_positive_median']:.6f}</code>.</div>
      <div class="cards">
        <div class="card"><span>BALF V4.3 recall</span><strong>{100 * summary['frozen_balf_recall']:.1f}%</strong><span>at the frozen threshold</span></div>
        <div class="card"><span>BALF Stage-2 ROC-AUC</span><strong>{summary['balf_v4_auc']:.4f}</strong><span>better ranking than V2/V3</span></div>
        <div class="card"><span>V4 misses recovered by V2</span><strong>{summary['v4_fn_v2_tp']:,}</strong><span>of {summary['v4_fn']:,} V4 false negatives</span></div>
        <div class="card"><span>New datasets with worse AB FPR</span><strong>{summary['new_sources_v4_worse_ab_fpr']}/{summary['new_sources_compared']}</strong><span>versus V3</span></div>
      </div>
      <p>V4.3 therefore failed as a deployable frozen classifier, even though its continuous BALF ranking is strong. A cutoff chosen after seeing BALF would recover performance but would be lockbox overfitting and is reported only as a diagnostic.</p>
      <div class="callout warning"><strong>V4.3 has no overall deployment advantage.</strong> It lowers author-NK FPR, but paired-alpha-beta FPR is higher than V3 in seven of eight new datasets, while recall is substantially lower. The gain is specific to one negative endpoint, not general specificity.</div>
    </section>
    <section>
      <h2>Was BALF used by V2 or V3?</h2>
      <div class="callout good"><strong>Your correction is valid for coefficient fitting.</strong> BALF was not used to fit V2 coefficients, and V3 contains the exact same V2 scaler and logistic coefficients. However, BALF was repeatedly evaluated and used later in V3 Round-12-versus-Round-14 promotion selection, so V3's BALF estimate is not untouched.</div>
      {dataframe_html(exposure)}
      <p>Numerical equality checks on the V2 and V3 base models found exact equality for scaler means, scaler scales, coefficients, and intercept. V3 differs through a conditional NK/TRDC gate and a threshold change from V2's <code>0.906359</code> to <code>0.936</code>. On BALF positives, that gate changed only {summary['v3_v2_score_changed']} of 852 scores.</p>
      {dataframe_html(contract_comparison)}
    </section>
    <section class="page-break">
      <h2>Failure localization</h2>
      {dataframe_html(headline, percent={'recall','precision','f1','abt_fpr','nk_fpr'})}
      <div class="grid2">
        <figure><img src="assets/balf_score_distributions.png"><figcaption>Scores are shown near the saturated upper tail. V4 positives sit close to, but often just below, the frozen cutoff.</figcaption></figure>
        <figure><img src="assets/balf_v4_vs_v2_scores.png"><figcaption>Most V4 false negatives remain high-confidence V2 positives.</figcaption></figure>
      </div>
      {dataframe_html(failure_summary, percent={'fraction_of_balf_gold'})}
      <ul>
        <li>Stage 1 passes {100 * summary['balf_stage1_pass']:.2f}% of BALF gold cells.</li>
        <li>CD4/Treg exclusions affect {100 * summary['balf_excluded']:.2f}%.</li>
        <li>Only {100 * summary['balf_raw_stage2_pass']:.2f}% pass the frozen Stage-2 cutoff directly; receptor rescue raises final recall to {100 * summary['frozen_balf_recall']:.2f}%.</li>
      </ul>
      <h3>Specificity by independent dataset</h3>
      <figure><img src="assets/per_source_ab_fpr_difference.png"><figcaption>Red indicates a higher V4.3 paired-alpha-beta false-positive rate; green indicates improvement. BALF is shown separately from the eight extension datasets.</figcaption></figure>
      {dataframe_html(per_source_ab, columns=['source_gse_id','n_abt_gold','v4_3_abt_fpr','v3_abt_fpr','v4_minus_v3_fpr','v4_only_fp','v3_only_fp','mcnemar_fdr'], percent={'v4_3_abt_fpr','v3_abt_fpr','v4_minus_v3_fpr'})}
      <p>Among the eight extension datasets, only <code>GSE159251</code> has lower paired-alpha-beta FPR under V4.3; {summary['new_sources_v4_significantly_worse_ab_fpr']} of the seven adverse differences remain significant after paired exact McNemar testing and BH correction. BALF also improves, but BALF is the cohort where V4.3 loses more than half of its gold positives at the frozen cutoff.</p>
    </section>
    <section>
      <h2>Threshold transfer failure</h2>
      <figure><img src="assets/threshold_transfer.png"><figcaption>The same threshold has materially different sensitivity and specificity behavior in development and BALF.</figcaption></figure>
      {dataframe_html(driver, percent={'recall','precision','f1','abt_fpr','nk_fpr','tier2_macro_nk_fpr'})}
      <div class="callout warning"><strong>Diagnostic only:</strong> the BALF best-F1 threshold is <code>{diagnostic_threshold:.6f}</code>. It gives high BALF performance, but it was found after opening BALF and must never become the released threshold or be used as evidence of generalization.</div>
      <h3>Why development selected 0.9941</h3>
      <p>At a Stage-2 threshold near 0.8006, development alpha-beta FPR remained above its 0.2% limit and tier-2 NK macro-FPR remained above its 2% limit. The threshold rose to 0.9941 to suppress high-scoring negative tails from heterogeneous development sources. BALF negatives were easier, so this development-driven cutoff became unnecessarily strict there.</p>
      <h3>OOF-to-final score mismatch</h3>
      {dataframe_html(score_shift, percent={'oof_above_frozen_threshold','final_model_above_frozen_threshold'})}
      <p>The OOF threshold was learned from fold models fitted with early stopping. The release candidate was then refitted as a single 800-tree model without early stopping, but no held-out probability calibration mapped that final model onto the OOF score scale. This is a model-construction defect even before considering biological source shift.</p>
    </section>
    <section class="page-break">
      <h2>Receptor dropout and V-gene asymmetry</h2>
      <div class="grid2">
        <figure><img src="assets/balf_recall_by_receptor_state.png"><figcaption>V4 depends strongly on simultaneous transcriptomic TRDC and TRDV detection.</figcaption></figure>
        <figure><img src="assets/balf_receptor_expression_boxes.png"><figcaption>Missed cells retain substantial receptor evidence but show weaker TRDC and different V-gene combinations.</figcaption></figure>
      </div>
      {dataframe_html(key_receptor, percent={'recall','ci_low','ci_high'})}
      <h3>Sequenced TRD V-gene calls</h3>
      {dataframe_html(trdv_recall[trdv_recall.model.isin(['V4.3 frozen','V2 high-F1'])], percent={'recall','ci_low','ci_high'})}
      <p>TRDV2 is the dominant productive V call in BALF. V4 recovers only {100 * summary['trdv2_v4_recall']:.1f}% of these cells, while the V2 linear model recovers {100 * summary['trdv2_v2_recall']:.1f}%. SHAP analysis shows that absent TRDV1 contributes a strong negative term even when TRDV2 and TRGV9 are present. The trees learned source-specific conjunctions rather than treating alternate gdT receptor programs as interchangeable evidence.</p>
      <figure><img src="assets/balf_shap_tp_fn_difference.png"><figcaption>Positive values support recovery in V4 TPs; negative values are stronger in V4 misses that V2 recovers.</figcaption></figure>
    </section>
    <section>
      <h2>Are the false negatives low-quality or weak truth?</h2>
      {dataframe_html(qc_stats, columns=['variable','v4_tp_n','v4_fn_n','v4_tp_median','v4_fn_median','rank_biserial_tp_minus_fn','mannwhitney_fdr'])}
      <p>No. V4 false negatives do not have fewer genes or fewer total counts; their mitochondrial fractions are similar. Cells missed by V4 but recovered by V2 have median productive TRD UMI support of {summary['v4_fn_v2_tp_trd_umi_median']:.1f}, versus {summary['v4_tp_trd_umi_median']:.1f} among V4 TPs. These are not one-UMI ambient receptor calls. The main expression difference is partial receptor dropout/asymmetry, particularly lower TRDC and absent TRDV1 in genuine TRDV2 cells.</p>
      <figure><img src="assets/balf_feature_detection_tp_fn.png"><figcaption>Detection-rate contrasts are tested with Fisher exact tests; expression contrasts use Mann-Whitney tests with BH correction.</figcaption></figure>
      {dataframe_html(top_features, columns=['gene','v4_tp_detected_fraction','v4_fn_detected_fraction','detected_fraction_difference_tp_minus_fn','v4_tp_median_expression','v4_fn_median_expression','fisher_fdr','mannwhitney_fdr'], percent={'v4_tp_detected_fraction','v4_fn_detected_fraction','detected_fraction_difference_tp_minus_fn'})}
    </section>
    <section class="page-break">
      <h2>Why V2 transfers better</h2>
      <ol>
        <li><strong>Smoother decision boundary.</strong> V2 is a standardized linear logistic model. It adds evidence across TRDV/TRGV genes and absence of many alpha/beta V genes. V4's shallow trees still form nonlinear receptor combinations and sharp magnitude splits.</li>
        <li><strong>Broader positive source used by V2.</strong> V2 fitted HRA005041, MalteGDT, and GDT_2020AUG_woCOV positives. V4 used GSE144469, HRA005041, and MalteGDT; it deliberately locked GDT_2020 out. This likely reduced exposure to the Vdelta2-rich program that resembles BALF, although BALF itself remained unseen.</li>
        <li><strong>Architectural mismatch.</strong> V4 Stage 2 was trained on gdT versus alpha-beta cells, not NK cells. NK limits were nevertheless imposed during threshold selection after Stage 1 allowed a substantial NK tail through. The cutoff, rather than the classifier, absorbed that conflict.</li>
        <li><strong>Calibration mismatch.</strong> Fold-model OOF scores selected the cutoff, while a structurally different final 800-tree fit produced deployment scores. No source-held-out Platt/isotonic calibration or fold ensemble preserved the OOF scale.</li>
      </ol>
      <p>The 29 TCR genes present in V2 but absent from V4 are not the explanation. They account for only {100 * summary['v2_only_abs_coefficient_share']:.2f}% of V2's absolute standardized coefficient mass. Setting all 29 to zero changed BALF V2 recall from {100 * summary['v2_original_recall']:.2f}% to {100 * summary['v2_ablated_recall']:.2f}%.</p>
    </section>
    <section>
      <h2>Implications</h2>
      <div class="callout"><strong>V4.3 remains non-promotable.</strong> Its BALF ranking suggests the receptor feature concept is viable, but the observed performance cannot be repaired by changing the threshold after lockbox review.</div>
      <p>A scientifically valid next version needs a new untouched evaluation cohort and should precommit:</p>
      <ul>
        <li>a fold-model ensemble or a separate calibration split so the selected operating point applies to the deployed predictor;</li>
        <li>source-held-out probability calibration and threshold stability as an explicit gate;</li>
        <li>alternate receptor-program handling that treats TRDV1 and TRDV2 as valid substitutes and is robust to TRDC dropout;</li>
        <li>NK negatives in Stage-2 training, or a stronger transcriptomic T-cell gate, instead of forcing Stage 2 to solve unseen NK tails through an extreme threshold;</li>
        <li>validation of recall by productive TRD V call, receptor-dropout state, source, donor, and tissue.</li>
      </ul>
      <p class="small">Input checksums and exact output paths are recorded in <code>Integrated_dataset/logs/gdT_prediction/gdtai_v4_3_recall_failure/summary.json</code>. No H5AD was modified.</p>
    </section>
    </main></body></html>"""
    index = REPORT_DIR / "index.html"
    index.write_text(body)
    pdf = REPORT_DIR / "gdtai_v4_3_recall_failure_report.pdf"
    result = subprocess.run(["/home/anaconda3/bin/weasyprint", str(index), str(pdf)], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"WeasyPrint failed: {result.stderr}")


def main() -> None:
    started = time.monotonic()
    for path in [TABLE_DIR, FIGURE_DIR, LOG_DIR, REPORT_DIR, ASSET_DIR]:
        path.mkdir(parents=True, exist_ok=True)

    required = [ATLAS, PREDICTIONS, LOCKBOX_MATRIX, LOCKBOX_FEATURES, LABELS, DEVELOPMENT_MATRIX, DEVELOPMENT_FEATURES, OOF, V4_CONTRACT, V2_MODEL, V3_MODEL]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing required inputs: {missing}")

    predictions = pd.read_parquet(PREDICTIONS)
    balf = predictions[predictions.source_gse_id.eq("BALF_BLOOD_COPD")].reset_index(drop=True)
    balf_gold = balf[balf.truth_class.eq("gdT_gold")].reset_index(drop=True)
    all_genes = pd.read_csv(LOCKBOX_FEATURES).sort_values("feature_index").gene.astype(str).tolist()
    lockbox_matrix = np.load(LOCKBOX_MATRIX, mmap_mode="r")
    balf_rows = predictions.index[predictions.source_gse_id.eq("BALF_BLOOD_COPD")].to_numpy(np.int64)
    balf_gold_rows = predictions.index[
        predictions.source_gse_id.eq("BALF_BLOOD_COPD") & predictions.truth_class.eq("gdT_gold")
    ].to_numpy(np.int64)
    balf_matrix = np.asarray(lockbox_matrix[balf_rows])
    balf_gold_matrix = np.asarray(lockbox_matrix[balf_gold_rows])
    states = receptor_states(balf_gold_matrix, all_genes)
    balf_gold = pd.concat([balf_gold, states], axis=1)

    metadata_columns = [
        "TRD_v", "TRG_v", "TRD_umis", "TRG_umis", "n_genes_by_counts", "total_counts",
        "pct_counts_mt", "source_cell_type", "tissue_harmonized_v2", "specimen_context_harmonized_v2",
        "technology_simple", "library_id", "sample_id", "donor_id",
    ]
    atlas_rows = balf_gold.atlas_row.to_numpy(np.int64)
    with h5py.File(ATLAS, "r") as handle:
        for column in metadata_columns:
            balf_gold[column] = decode_subset(handle["obs"][column], atlas_rows)
    for column in ["TRD_umis", "TRG_umis", "n_genes_by_counts", "total_counts", "pct_counts_mt"]:
        balf_gold[column] = pd.to_numeric(balf_gold[column], errors="coerce")

    balf_gold["failure_group"] = np.where(
        balf_gold.v4_3_balanced,
        "V4_TP",
        np.where(balf_gold.v2_high_f1, "V4_FN_V2_TP", "V4_FN_V2_FN"),
    )
    balf_gold.to_parquet(TABLE_DIR / "balf_gold_cell_diagnostics.parquet", index=False, compression="zstd")

    v4_contract = json.loads(V4_CONTRACT.read_text())
    v2 = load_pickle(V2_MODEL)
    v3 = load_pickle(V3_MODEL)
    v2_pipe = v2["base_model"]["model_object"]
    v3_base = v3["model_object"].base_model
    equality = {
        "scaler_mean": np.array_equal(v2_pipe.named_steps["standardscaler"].mean_, v3_base.named_steps["standardscaler"].mean_),
        "scaler_scale": np.array_equal(v2_pipe.named_steps["standardscaler"].scale_, v3_base.named_steps["standardscaler"].scale_),
        "coefficient": np.array_equal(v2_pipe.named_steps["logisticregression"].coef_, v3_base.named_steps["logisticregression"].coef_),
        "intercept": np.array_equal(v2_pipe.named_steps["logisticregression"].intercept_, v3_base.named_steps["logisticregression"].intercept_),
    }
    if not all(equality.values()):
        raise RuntimeError("V2 and V3 base models are not numerically identical")

    exposure = pd.DataFrame([
        {
            "model": "V2 base",
            "coefficient_fit_sources": "HRA005041; MalteGDT; GDT_2020AUG_woCOV",
            "BALF_in_coefficient_fit": False,
            "BALF_in_threshold_selection": False,
            "BALF_in_promotion_selection": False,
            "interpretation": "BALF was fully unseen by V2 development.",
        },
        {
            "model": "V3 balanced",
            "coefficient_fit_sources": "Exact V2 base coefficients",
            "BALF_in_coefficient_fit": False,
            "BALF_in_threshold_selection": False,
            "BALF_in_promotion_selection": True,
            "interpretation": "BALF influenced later Round 12/14 promotion, not fitted weights.",
        },
        {
            "model": "V4.3",
            "coefficient_fit_sources": "GSE144469; HRA005041; MalteGDT (positive truth)",
            "BALF_in_coefficient_fit": False,
            "BALF_in_threshold_selection": False,
            "BALF_in_promotion_selection": False,
            "interpretation": "BALF was a once-opened V4 lockbox.",
        },
    ])
    exposure.to_csv(TABLE_DIR / "model_balf_exposure_audit.csv", index=False)

    contract_comparison = pd.DataFrame([
        {"property": "Base algorithm", "V2": "Standardized logistic", "V3": "Exact V2 logistic + conditional gate", "V4.3": "Two XGBoost stages"},
        {"property": "Stage-2 fitted genes", "V2": "187", "V3": "V2's 187 for base score", "V4.3": "158"},
        {"property": "NK/cytotoxic genes in base score", "V2": "No", "V3": "No; only conditional gate uses NK score", "V4.3": "No"},
        {"property": "NK cells in gdT-vs-ab Stage-2 fitting", "V2": "No explicit NK class", "V3": "Conditional post-gate", "V4.3": "No"},
        {"property": "Probability calibration", "V2": "Implicit logistic", "V3": "Implicit logistic", "V4.3": "None"},
        {"property": "Deployed model matches threshold-generating model", "V2": "Yes", "V3": "Yes", "V4.3": "No: OOF fold models vs final 800-tree refit"},
    ])
    contract_comparison.to_csv(TABLE_DIR / "model_contract_comparison.csv", index=False)

    headline_rows = []
    for model, column in [
        ("V4.3 frozen", "v4_3_balanced"),
        ("V3 balanced", "v3_balanced"),
        ("V2 high-F1", "v2_high_f1"),
        ("V2 high-purity", "v2_high_purity"),
    ]:
        headline_rows.append({"model": model, **binary_metrics(balf, balf[column].to_numpy(bool))})
    headline = pd.DataFrame(headline_rows)
    headline.to_csv(TABLE_DIR / "balf_headline_metrics.csv", index=False)

    ab_frame = predictions[predictions.truth_class.eq("abT_gold")]
    per_source_ab = ab_frame.groupby("source_gse_id", as_index=False).agg(
        n_abt_gold=("cell_id", "size"),
        v4_3_abt_fpr=("v4_3_balanced", "mean"),
        v3_abt_fpr=("v3_balanced", "mean"),
        v2_high_f1_abt_fpr=("v2_high_f1", "mean"),
        v2_high_purity_abt_fpr=("v2_high_purity", "mean"),
    )
    per_source_ab["v4_minus_v3_fpr"] = per_source_ab.v4_3_abt_fpr - per_source_ab.v3_abt_fpr
    per_source_ab["v4_better_than_v3"] = per_source_ab.v4_3_abt_fpr < per_source_ab.v3_abt_fpr
    discordance = []
    for source, group in ab_frame.groupby("source_gse_id", sort=True):
        v4_call = group.v4_3_balanced.to_numpy(bool)
        v3_call = group.v3_balanced.to_numpy(bool)
        v4_only = int((v4_call & ~v3_call).sum())
        v3_only = int((~v4_call & v3_call).sum())
        discordance.append({
            "source_gse_id": source,
            "v4_only_fp": v4_only,
            "v3_only_fp": v3_only,
            "mcnemar_exact_p": float(binomtest(v4_only, v4_only + v3_only, 0.5).pvalue) if v4_only + v3_only else 1.0,
        })
    discordance = pd.DataFrame(discordance)
    discordance["mcnemar_fdr"] = bh(discordance.mcnemar_exact_p)
    per_source_ab = per_source_ab.merge(discordance, on="source_gse_id", how="left", validate="one_to_one")
    per_source_ab["cohort_type"] = np.where(
        per_source_ab.source_gse_id.eq("BALF_BLOOD_COPD"), "BALF reused benchmark", "new extension dataset"
    )
    per_source_ab.to_csv(TABLE_DIR / "paired_ab_fpr_by_source_v4_vs_v3.csv", index=False)

    diagnostic_threshold, diagnostic_metrics = diagnostic_best_threshold(balf)
    operating_rows = []
    oof = pd.read_parquet(OOF).reset_index(drop=True)
    threshold_specs = [
        ("BALF diagnostic best-F1", diagnostic_threshold),
        ("V3-matching diagnostic 0.5624", 0.5624107122421265),
        ("0.8006 diagnostic", 0.8005564212799072),
        ("Frozen 0.9941", FROZEN_THRESHOLD),
    ]
    for label, threshold in threshold_specs:
        for cohort, frame, caller in [("BALF", balf, full_v4_call), ("Development OOF", oof, full_oof_call)]:
            metrics = binary_metrics(frame, caller(frame, threshold))
            if cohort == "Development OOF":
                nk2 = frame.truth_class.eq("nk_tier2").to_numpy()
                local_calls = caller(frame, threshold)
                tier2 = frame.loc[nk2, ["source_gse_id"]].copy()
                tier2["call"] = local_calls[nk2]
                metrics["tier2_macro_nk_fpr"] = tier2.groupby("source_gse_id").call.mean().mean()
            else:
                metrics["tier2_macro_nk_fpr"] = math.nan
            operating_rows.append({"threshold_label": label, "threshold": threshold, "cohort": cohort, **metrics})
    operating = pd.DataFrame(operating_rows)
    operating.to_csv(TABLE_DIR / "threshold_transfer_diagnostics.csv", index=False)

    receptor_recall = group_recall(balf_gold, "TRDC_TRDV_state")
    receptor_recall.to_csv(TABLE_DIR / "balf_recall_by_trdc_trdv_state.csv", index=False)
    trdv_recall = group_recall(balf_gold, "TRD_v")
    trdv_recall.to_csv(TABLE_DIR / "balf_recall_by_productive_trd_v.csv", index=False)
    trgv_recall = group_recall(balf_gold, "TRG_v")
    trgv_recall.to_csv(TABLE_DIR / "balf_recall_by_productive_trg_v.csv", index=False)
    donor_recall = group_recall(balf_gold, "donor_id")
    donor_recall.to_csv(TABLE_DIR / "balf_recall_by_donor.csv", index=False)
    sample_recall = group_recall(balf_gold, "sample_id")
    sample_recall.to_csv(TABLE_DIR / "balf_recall_by_sample.csv", index=False)

    failure_summary = balf_gold.groupby("failure_group", as_index=False).agg(
        n_cells=("cell_id", "size"),
        median_stage1_score=("v4_3_stage1_score", "median"),
        median_stage2_score=("v4_3_stage2_score", "median"),
        median_v2_score=("v2_score", "median"),
        median_n_genes=("n_genes_by_counts", "median"),
        median_total_counts=("total_counts", "median"),
        median_TRD_umis=("TRD_umis", "median"),
        median_TRG_umis=("TRG_umis", "median"),
    )
    failure_summary["fraction_of_balf_gold"] = failure_summary.n_cells / len(balf_gold)
    failure_summary.to_csv(TABLE_DIR / "balf_failure_group_summary.csv", index=False)

    feature_stats = feature_statistics(balf_gold, balf_gold_matrix, all_genes)
    feature_stats.to_csv(TABLE_DIR / "balf_feature_tp_fn_statistics.csv", index=False)
    qc_stats = qc_statistics(balf_gold)
    qc_stats.to_csv(TABLE_DIR / "balf_qc_tp_fn_statistics.csv", index=False)

    shap_table = stage2_shap(balf_gold, balf_gold_matrix, all_genes, v4_contract)
    shap_table.to_csv(TABLE_DIR / "balf_v4_stage2_shap_by_failure_group.csv", index=False)
    coefficient_table = v2_coefficients(v2, set(v4_contract["stage2_feature_names"]))
    coefficient_table.to_csv(TABLE_DIR / "v2_logistic_coefficients.csv", index=False)

    v2_genes = v2["base_model"]["gene_names"]
    lookup = {gene: index for index, gene in enumerate(all_genes)}
    v2_matrix = np.asarray(balf_gold_matrix[:, [lookup[gene] for gene in v2_genes]], dtype=np.float32)
    original_v2_score = v2_pipe.predict_proba(v2_matrix)[:, 1]
    v2_only = np.asarray([index for index, gene in enumerate(v2_genes) if gene not in set(v4_contract["stage2_feature_names"])], dtype=np.int64)
    ablated = v2_matrix.copy()
    ablated[:, v2_only] = 0
    ablated_v2_score = v2_pipe.predict_proba(ablated)[:, 1]
    v2_threshold = float(v2["operating_modes"]["high_f1"]["threshold"])
    feature_contract = pd.DataFrame([
        {"comparison": "V2 genes", "n_genes": len(v2_genes), "detail": "182 TCR genes + CD3D/E/G, CD4, FOXP3"},
        {"comparison": "V4.3 Stage-2 genes", "n_genes": len(v4_contract["stage2_feature_names"]), "detail": "153 TCR genes + same five controls"},
        {"comparison": "V2-only genes", "n_genes": len(v2_only), "detail": "; ".join(np.asarray(v2_genes)[v2_only])},
    ])
    feature_contract.to_csv(TABLE_DIR / "feature_contract_counts.csv", index=False)

    source_scores = score_quantiles(
        predictions[predictions.truth_class.eq("gdT_gold")], "External lockbox"
    )
    oof_score = oof[oof.truth_class.eq("gdT_gold")].copy()
    oof_score = oof_score.rename(columns={"stage2_score": "v4_3_stage2_score"})
    oof_score["v3_balanced_score"] = np.nan
    oof_score["v2_score"] = np.nan
    source_scores = pd.concat([source_scores, score_quantiles(oof_score, "Development OOF")], ignore_index=True)
    source_scores.to_csv(TABLE_DIR / "positive_score_quantiles_by_source.csv", index=False)

    final_score, score_shift = final_model_on_oof(oof, v4_contract)
    score_shift.to_csv(TABLE_DIR / "oof_vs_final_model_score_shift.csv", index=False)

    threshold_drivers = operating.copy()
    threshold_drivers.to_csv(TABLE_DIR / "development_threshold_drivers.csv", index=False)

    balf_primary = balf[balf.truth_class.isin(["gdT_gold", "abT_gold"])].copy()
    y_balf = balf_primary.truth_class.eq("gdT_gold").to_numpy(np.int8)
    ranking = []
    for model, column in [("V4.3 Stage 2", "v4_3_stage2_score"), ("V3 balanced", "v3_balanced_score"), ("V2 base", "v2_score")]:
        score = balf_primary[column].to_numpy(float)
        ranking.append({"model": model, "roc_auc": roc_auc_score(y_balf, score), "average_precision": average_precision_score(y_balf, score)})
    ranking = pd.DataFrame(ranking)
    ranking.to_csv(TABLE_DIR / "balf_ranking_metrics.csv", index=False)

    make_figures(
        balf, balf_matrix, all_genes, receptor_recall, feature_stats, shap_table, oof,
        per_source_ab, diagnostic_threshold,
    )

    trdv2 = trdv_recall[trdv_recall.TRD_v.astype(str).str.upper().eq("TRDV2")].set_index("model")
    failure_lookup = failure_summary.set_index("failure_group")
    v2_abs = coefficient_table.standardized_coefficient.abs()
    new_source_ab = per_source_ab[per_source_ab.cohort_type.eq("new extension dataset")]
    summary = {
        "status": "PASS_V4_3_RECALL_FAILURE_INVESTIGATION",
        "scientific_conclusion": "threshold_transfer_and_nonlinear_receptor_calibration_failure",
        "model_promoted": False,
        "threshold_retuned": False,
        "diagnostic_threshold_not_deployable": diagnostic_threshold,
        "frozen_balf_recall": float(headline.set_index("model").loc["V4.3 frozen", "recall"]),
        "balf_v4_auc": float(ranking.set_index("model").loc["V4.3 Stage 2", "roc_auc"]),
        "balf_v4_ap": float(ranking.set_index("model").loc["V4.3 Stage 2", "average_precision"]),
        "balf_v4_positive_median": float(balf_gold.v4_3_stage2_score.median()),
        "balf_stage1_pass": float((balf_gold.v4_3_stage1_score >= STAGE1_THRESHOLD).mean()),
        "balf_excluded": float(balf_gold.v4_3_cd4_treg_excluded.mean()),
        "balf_raw_stage2_pass": float((
            (balf_gold.v4_3_stage1_score >= STAGE1_THRESHOLD)
            & (balf_gold.v4_3_stage2_score >= FROZEN_THRESHOLD)
            & ~balf_gold.v4_3_cd4_treg_excluded
        ).mean()),
        "v4_fn": int((~balf_gold.v4_3_balanced).sum()),
        "v4_fn_v2_tp": int((~balf_gold.v4_3_balanced & balf_gold.v2_high_f1).sum()),
        "trdv2_v4_recall": float(trdv2.loc["V4.3 frozen", "recall"]),
        "trdv2_v2_recall": float(trdv2.loc["V2 high-F1", "recall"]),
        "v4_fn_v2_tp_trd_umi_median": float(failure_lookup.loc["V4_FN_V2_TP", "median_TRD_umis"]),
        "v4_tp_trd_umi_median": float(failure_lookup.loc["V4_TP", "median_TRD_umis"]),
        "v3_v2_score_changed": int((np.abs(balf_gold.v3_balanced_score - balf_gold.v2_score) > 1e-8).sum()),
        "v2_only_gene_count": int(len(v2_only)),
        "v2_only_abs_coefficient_share": float(v2_abs[~coefficient_table.present_in_v4_stage2].sum() / v2_abs.sum()),
        "v2_original_recall": float((original_v2_score >= v2_threshold).mean()),
        "v2_ablated_recall": float((ablated_v2_score >= v2_threshold).mean()),
        "new_sources_compared": int(len(new_source_ab)),
        "new_sources_v4_worse_ab_fpr": int((new_source_ab.v4_minus_v3_fpr > 0).sum()),
        "new_sources_v4_better_ab_fpr": int((new_source_ab.v4_minus_v3_fpr < 0).sum()),
        "new_sources_v4_significantly_worse_ab_fpr": int(((new_source_ab.v4_minus_v3_fpr > 0) & (new_source_ab.mcnemar_fdr < 0.05)).sum()),
        "v2_v3_base_exact_equality": equality,
        "input_hashes": {
            "predictions": sha256(PREDICTIONS),
            "lockbox_matrix": sha256(LOCKBOX_MATRIX),
            "v4_contract": sha256(V4_CONTRACT),
            "v2_model": sha256(V2_MODEL),
            "v3_model": sha256(V3_MODEL),
            "oof_predictions": sha256(OOF),
        },
        "runtime_seconds": time.monotonic() - started,
        "no_h5ad_modified": True,
    }
    render_report(
        exposure, headline, operating, receptor_recall, trdv_recall, failure_summary,
        feature_stats, qc_stats, threshold_drivers, score_shift, contract_comparison, per_source_ab,
        diagnostic_threshold, summary,
    )
    summary["html_report"] = str(REPORT_DIR / "index.html")
    summary["pdf_report"] = str(REPORT_DIR / "gdtai_v4_3_recall_failure_report.pdf")
    summary["runtime_seconds"] = time.monotonic() - started
    (LOG_DIR / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    (LOG_DIR / "conclusion.md").write_text(
        "# gdTAI V4.3 recall-failure conclusion\n\n"
        "V4.3 ranks BALF_BLOOD_COPD receptor-gold cells well but fails at its frozen operating threshold. "
        "The development-derived XGBoost threshold does not transfer to the final refit or to new source distributions. "
        "Receptor dropout and TRDV1/TRDV2 asymmetry aggravate the loss, but BALF false negatives retain strong TCR UMI support and are not generally low-quality cells. "
        "V2 and the V3 base coefficients were not trained on BALF; V3 did reuse BALF for later promotion selection. "
        "No diagnostic threshold in this analysis is eligible for deployment.\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
