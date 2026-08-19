#!/usr/bin/env python3
"""Render the final gdTAI V4.3 common-lockbox evaluation and conclusion."""

from __future__ import annotations

import html
import json
import subprocess
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score, roc_curve


ROOT = Path(__file__).resolve().parents[2]
TABLE_ROOT = ROOT / "Integrated_dataset/tables/gdT_prediction"
FIGURE_DIR = ROOT / "Integrated_dataset/figures/gdT_prediction/gdtai_v4_3_final_evaluation"
LOG_DIR = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_3_final_evaluation"
STATIC_DIR = ROOT / "gdT_prediction/gdtai_v4_3_final_evaluation"
PREDICTIONS = TABLE_ROOT / "gdtai_v4_3_common_lockbox/common_lockbox_predictions.parquet"
OVERALL = TABLE_ROOT / "gdtai_v4_3_common_lockbox/overall_metrics.csv"
PER_SOURCE = TABLE_ROOT / "gdtai_v4_3_common_lockbox/per_source_metrics.csv"
BOOTSTRAP = TABLE_ROOT / "gdtai_v4_3_common_lockbox/cluster_bootstrap_summary.csv"
DECISION = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_3_common_lockbox/promotion_decision.json"
V4_CONTRACT = ROOT / "Integrated_dataset/models/gdT_prediction/gdtai_v4_3_receptor_isolated/final_development/model_contract.json"
LOCKBOX = TABLE_ROOT / "gdtai_v4_3_lockbox/lockbox_manifest.parquet"
OLD_LABELS = TABLE_ROOT / "gdtai_v4_2_ground_truth/v4_2_label_manifest.parquet"
OLD_FEATURES = TABLE_ROOT / "gdtai_v4_2_training/feature_manifest.csv"
OLD_CACHE = ROOT / "Integrated_dataset/cache/gdT_prediction/gdtai_v4_2_training/gene_features.npy"
NEW_FEATURES = ROOT / "Integrated_dataset/cache/gdT_prediction/gdtai_v4_3_common_lockbox/lockbox_feature_manifest.csv"
NEW_CACHE = ROOT / "Integrated_dataset/cache/gdT_prediction/gdtai_v4_3_common_lockbox/lockbox_gene_features.npy"
MODEL_LABELS = {
    "v4_3_balanced": "V4.3 frozen",
    "v3_balanced": "V3 balanced",
    "v2_high_f1": "V2 high-F1",
    "v2_high_purity": "V2 high-purity",
}
COLORS = {
    "v4_3_balanced": "#b42318",
    "v3_balanced": "#176b87",
    "v2_high_f1": "#4d7c0f",
    "v2_high_purity": "#6b7280",
}


def setup() -> None:
    for path in (FIGURE_DIR, LOG_DIR, STATIC_DIR / "assets"):
        path.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "figure.dpi": 150,
        "savefig.dpi": 240,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def pct(value: float) -> str:
    return f"{100 * value:.2f}%"


def savefig(fig: plt.Figure, name: str) -> Path:
    path = FIGURE_DIR / name
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    target = STATIC_DIR / "assets" / name
    target.write_bytes(path.read_bytes())
    return path


def implementation_integrity(v4: dict) -> pd.DataFrame:
    labels = pd.read_parquet(OLD_LABELS, columns=["cell_id"])
    lockbox = pd.read_parquet(LOCKBOX, columns=["cell_id"])
    old = np.load(OLD_CACHE, mmap_mode="r")
    new = np.load(NEW_CACHE, mmap_mode="r")
    old_genes = pd.read_csv(OLD_FEATURES).sort_values("feature_index").gene.astype(str).tolist()
    new_genes = pd.read_csv(NEW_FEATURES).sort_values("feature_index").gene.astype(str).tolist()
    used = []
    for gene in [*v4["stage1_feature_names"], *v4["stage2_feature_names"]]:
        if gene not in used:
            used.append(gene)
    old_columns = np.asarray([old_genes.index(gene) for gene in used])
    new_columns = np.asarray([new_genes.index(gene) for gene in used])
    lookup = pd.Series(np.arange(len(labels), dtype=np.int64), index=labels.cell_id)
    mapped = lookup.reindex(lockbox.cell_id).to_numpy(float)
    present = np.isfinite(mapped)
    old_rows = mapped[present].astype(np.int64)
    new_rows = np.flatnonzero(present)
    maximum = 0.0
    differences = 0
    compared = 0
    for start in range(0, len(old_rows), 20_000):
        first = np.asarray(old[old_rows[start:start + 20_000]][:, old_columns])
        second = np.asarray(new[new_rows[start:start + 20_000]][:, new_columns])
        delta = np.abs(first - second)
        maximum = max(maximum, float(delta.max()))
        differences += int((delta > 1e-6).sum())
        compared += delta.size
    result = pd.DataFrame([
        {"check": "lockbox cells matching original frozen cache", "value": int(present.sum()), "result": "PASS"},
        {"check": "new author-NK cells requiring direct extraction", "value": int((~present).sum()), "result": "EXPECTED"},
        {"check": "effective V4 genes compared", "value": len(used), "result": "PASS"},
        {"check": "feature values compared", "value": compared, "result": "PASS"},
        {"check": "values differing by >1e-6", "value": differences, "result": "PASS" if differences == 0 else "FAIL"},
        {"check": "maximum absolute feature difference", "value": maximum, "result": "PASS" if maximum == 0 else "FAIL"},
    ])
    result.to_csv(TABLE_ROOT / "gdtai_v4_3_common_lockbox/implementation_integrity.csv", index=False)
    return result


def stage_decomposition(predictions: pd.DataFrame, v4: dict) -> pd.DataFrame:
    positive = predictions[predictions.truth_class.eq("gdT_gold")].copy()
    positive["stage1_pass"] = positive.v4_3_stage1_score >= float(v4["stage1_threshold"])
    positive["not_excluded"] = ~positive.v4_3_cd4_treg_excluded
    positive["base_eligible"] = positive.stage1_pass & positive.not_excluded
    positive["raw_stage2_pass"] = positive.v4_3_stage2_score >= float(v4["stage2_threshold"])
    rows = []
    for source, group in positive.groupby("source_gse_id", sort=True):
        rows.append({
            "source_gse_id": source,
            "n_gdT_gold": len(group),
            "stage1_pass": float(group.stage1_pass.mean()),
            "not_cd4_treg_excluded": float(group.not_excluded.mean()),
            "raw_stage2_pass": float((group.base_eligible & group.raw_stage2_pass).mean()),
            "receptor_rescue": float((group.base_eligible & group.v4_3_receptor_rescue).mean()),
            "final_recall": float(group.v4_3_balanced.mean()),
            "v3_recall": float(group.v3_balanced.mean()),
            "v2_high_f1_recall": float(group.v2_high_f1.mean()),
            "v3_recovered_v4_missed": int((~group.v4_3_balanced & group.v3_balanced).sum()),
        })
    result = pd.DataFrame(rows)
    result.to_csv(TABLE_ROOT / "gdtai_v4_3_common_lockbox/stage_decomposition.csv", index=False)
    return result


def ranking_metrics(predictions: pd.DataFrame) -> pd.DataFrame:
    truth = predictions.truth_class.eq("gdT_gold").astype(int)
    rows = []
    for model, score in [
        ("v4_3_balanced", "v4_3_stage2_score"),
        ("v3_balanced", "v3_balanced_score"),
        ("v2_high_f1", "v2_score"),
    ]:
        rows.append({
            "model": model,
            "roc_auc": roc_auc_score(truth, predictions[score]),
            "average_precision": average_precision_score(truth, predictions[score]),
        })
    result = pd.DataFrame(rows)
    result.to_csv(TABLE_ROOT / "gdtai_v4_3_common_lockbox/ranking_metrics.csv", index=False)
    return result


def threshold_frontier(predictions: pd.DataFrame, v4: dict, overall: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    truth = predictions.truth_class.eq("gdT_gold").to_numpy()
    abt = predictions.truth_class.eq("abT_gold").to_numpy()
    nk = predictions.truth_class.eq("nk_lockbox").to_numpy()
    score = predictions.v4_3_stage2_score.to_numpy()
    base = (
        (predictions.v4_3_stage1_score.to_numpy() >= float(v4["stage1_threshold"]))
        & ~predictions.v4_3_cd4_treg_excluded.to_numpy()
    )
    fixed = base & predictions.v4_3_receptor_rescue.to_numpy()
    variable = base & ~fixed
    index = np.flatnonzero(variable)
    order = index[np.argsort(-score[index], kind="mergesort")]
    tp = int((fixed & truth).sum()) + np.cumsum(truth[order])
    ab_fp = int((fixed & abt).sum()) + np.cumsum(abt[order])
    nk_fp = int((fixed & nk).sum()) + np.cumsum(nk[order])
    fp = ab_fp + nk_fp
    recall = tp / truth.sum()
    precision = tp / (tp + fp)
    f1 = 2 * precision * recall / (precision + recall)
    boundary = np.r_[score[order][:-1] != score[order][1:], True]
    valid = np.flatnonzero(boundary)
    frontier = pd.DataFrame({
        "threshold": score[order][valid],
        "f1": f1[valid],
        "recall": recall[valid],
        "precision": precision[valid],
        "abt_fpr": ab_fp[valid] / abt.sum(),
        "author_nk_fpr": nk_fp[valid] / nk.sum(),
    })
    best = frontier.loc[frontier.f1.idxmax()]
    v3_recall = float(overall.set_index("model").loc["v3_balanced", "recall"])
    diagnostics = [{"diagnostic": "posthoc_maximum_f1", **best.to_dict()}]
    for label, target in [("minimum_70pct_recall", 0.70), ("match_v3_recall", v3_recall), ("minimum_80pct_recall", 0.80)]:
        eligible = frontier[frontier.recall >= target]
        if len(eligible):
            selected = eligible.sort_values(["abt_fpr", "author_nk_fpr", "threshold"], ascending=[True, True, False]).iloc[0]
            diagnostics.append({"diagnostic": label, **selected.to_dict()})
    frozen = overall.set_index("model").loc["v4_3_balanced"]
    diagnostics.append({
        "diagnostic": "frozen_operating_point",
        "threshold": float(v4["stage2_threshold"]),
        "f1": float(frozen.f1),
        "recall": float(frozen.recall),
        "precision": float(frozen.precision),
        "abt_fpr": float(frozen.abt_fpr),
        "author_nk_fpr": float(frozen.author_nk_fpr),
    })
    diagnostic_table = pd.DataFrame(diagnostics)
    diagnostic_table.to_csv(TABLE_ROOT / "gdtai_v4_3_common_lockbox/posthoc_threshold_diagnostics.csv", index=False)
    sampled = frontier.iloc[np.unique(np.linspace(0, len(frontier) - 1, min(2500, len(frontier))).astype(int))]
    sampled.to_csv(TABLE_ROOT / "gdtai_v4_3_common_lockbox/posthoc_threshold_frontier_sample.csv", index=False)
    return diagnostic_table, sampled


def feature_detection(v4: dict) -> pd.DataFrame:
    panel = ["TRDC", "TRDV1", "TRDV2", "TRDV3", "TRGC1", "TRGC2", "TRGV3", "TRGV9", "CD3D", "CD3E", "CD3G"]
    labels = pd.read_parquet(OLD_LABELS)
    old = np.load(OLD_CACHE, mmap_mode="r")
    old_genes = pd.read_csv(OLD_FEATURES).sort_values("feature_index").gene.astype(str).tolist()
    old_columns = [old_genes.index(gene) for gene in panel]
    rows = []
    development = labels[labels.allow_fit & labels.truth_class.eq("gdT_gold")]
    for source, indexes in development.groupby("source_gse_id").groups.items():
        matrix = np.asarray(old[np.asarray(indexes, dtype=np.int64)][:, old_columns])
        for column, gene in enumerate(panel):
            rows.append({"source_gse_id": source, "role": "development", "gene": gene, "n_cells": len(matrix), "detected_fraction": float((matrix[:, column] > 0).mean()), "median_expression": float(np.median(matrix[:, column]))})
    lockbox = pd.read_parquet(LOCKBOX)
    new = np.load(NEW_CACHE, mmap_mode="r")
    new_genes = pd.read_csv(NEW_FEATURES).sort_values("feature_index").gene.astype(str).tolist()
    new_columns = [new_genes.index(gene) for gene in panel]
    locked_positive = lockbox[lockbox.truth_class.eq("gdT_gold")]
    for source, indexes in locked_positive.groupby("source_gse_id").groups.items():
        matrix = np.asarray(new[np.asarray(indexes, dtype=np.int64)][:, new_columns])
        for column, gene in enumerate(panel):
            rows.append({"source_gse_id": source, "role": "lockbox", "gene": gene, "n_cells": len(matrix), "detected_fraction": float((matrix[:, column] > 0).mean()), "median_expression": float(np.median(matrix[:, column]))})
    result = pd.DataFrame(rows)
    result.to_csv(TABLE_ROOT / "gdtai_v4_3_common_lockbox/positive_source_feature_detection.csv", index=False)
    return result


def model_gain(v4: dict) -> pd.DataFrame:
    booster = xgb.Booster()
    booster.load_model(ROOT / "Integrated_dataset/models/gdT_prediction/gdtai_v4_3_receptor_isolated/final_development/stage2_receptor_classifier.ubj")
    rows = []
    for key, gain in booster.get_score(importance_type="gain").items():
        index = int(key[1:])
        rows.append({"gene": v4["stage2_feature_names"][index], "gain": float(gain)})
    result = pd.DataFrame(rows).sort_values("gain", ascending=False)
    result.to_csv(TABLE_ROOT / "gdtai_v4_3_common_lockbox/stage2_feature_gain.csv", index=False)
    return result


def make_figures(predictions: pd.DataFrame, overall: pd.DataFrame, source: pd.DataFrame, bootstrap: pd.DataFrame,
                 stage: pd.DataFrame, ranking: pd.DataFrame, diagnostics: pd.DataFrame,
                 frontier: pd.DataFrame, detection: pd.DataFrame, gain: pd.DataFrame) -> None:
    order = list(MODEL_LABELS)
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.5), constrained_layout=True)
    for ax, metric, title in zip(axes, ["recall", "precision", "f1"], ["Recall", "Precision", "F1"]):
        values = overall.set_index("model").loc[order, metric]
        ax.bar(range(len(order)), values, color=[COLORS[key] for key in order])
        ax.set_xticks(range(len(order)), [MODEL_LABELS[key] for key in order], rotation=25, ha="right")
        ax.set_ylim(0, 1.02); ax.set_title(title)
        for i, value in enumerate(values): ax.text(i, value + 0.018, f"{value:.3f}", ha="center", fontsize=8)
    savefig(fig, "overall_performance.png")

    positive = source[source.truth_class.eq("gdT_gold")].copy()
    sources = positive.source_gse_id.drop_duplicates().tolist()
    fig, ax = plt.subplots(figsize=(9.5, 4.2), constrained_layout=True)
    width = 0.19; x = np.arange(len(sources))
    for offset, model in enumerate(order):
        values = positive[positive.model.eq(model)].set_index("source_gse_id").loc[sources, "positive_rate"]
        ax.bar(x + (offset - 1.5) * width, values, width, label=MODEL_LABELS[model], color=COLORS[model])
    ax.set_xticks(x, sources); ax.set_ylim(0, 1.02); ax.set_ylabel("Recall"); ax.set_title("External gdT-gold recall by source")
    ax.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.18))
    savefig(fig, "positive_source_recall.png")

    negatives = source[source.truth_class.eq("abT_gold")].copy()
    sources = negatives.source_gse_id.drop_duplicates().tolist()
    fig, ax = plt.subplots(figsize=(11, 4.5), constrained_layout=True)
    width = 0.19; x = np.arange(len(sources))
    for offset, model in enumerate(order):
        values = negatives[negatives.model.eq(model)].set_index("source_gse_id").loc[sources, "positive_rate"]
        ax.bar(x + (offset - 1.5) * width, 100 * values, width, label=MODEL_LABELS[model], color=COLORS[model])
    ax.set_xticks(x, sources, rotation=30, ha="right"); ax.set_ylabel("False-positive rate (%)"); ax.set_title("Paired alpha-beta TCR false positives by source")
    ax.legend(frameon=False, ncol=4)
    savefig(fig, "abt_fpr_by_source.png")

    nk = source[source.truth_class.eq("nk_lockbox")].copy()
    fig, ax = plt.subplots(figsize=(7.8, 4), constrained_layout=True)
    sources = nk.source_gse_id.drop_duplicates().tolist(); x = np.arange(len(sources)); width = 0.19
    for offset, model in enumerate(order):
        values = nk[nk.model.eq(model)].set_index("source_gse_id").loc[sources, "positive_rate"]
        ax.bar(x + (offset - 1.5) * width, 100 * values, width, label=MODEL_LABELS[model], color=COLORS[model])
    ax.set_xticks(x, sources); ax.set_ylabel("False-positive rate (%)"); ax.set_title("Independent author-NK false positives"); ax.legend(frameon=False, ncol=2)
    savefig(fig, "nk_fpr_by_source.png")

    fig, ax = plt.subplots(figsize=(9.2, 4.2), constrained_layout=True)
    columns = ["stage1_pass", "not_cd4_treg_excluded", "raw_stage2_pass", "receptor_rescue", "final_recall"]
    labels = ["Stage 1 pass", "Not excluded", "Raw Stage 2", "Receptor rescue", "Final"]
    x = np.arange(len(stage)); width = 0.16
    for offset, (column, label) in enumerate(zip(columns, labels)):
        ax.bar(x + (offset - 2) * width, stage[column], width, label=label)
    ax.set_xticks(x, stage.source_gse_id); ax.set_ylim(0, 1.05); ax.set_ylabel("Fraction of gdT gold"); ax.set_title("V4.3 recall loss occurs in Stage 2"); ax.legend(frameon=False, ncol=3)
    savefig(fig, "stage_decomposition.png")

    truth = predictions.truth_class.eq("gdT_gold").astype(int)
    curves = [("v4_3_balanced", "v4_3_stage2_score"), ("v3_balanced", "v3_balanced_score"), ("v2_high_f1", "v2_score")]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4), constrained_layout=True)
    for model, column in curves:
        fpr, tpr, _ = roc_curve(truth, predictions[column]); precision, recall, _ = precision_recall_curve(truth, predictions[column])
        values = ranking.set_index("model").loc[model]
        axes[0].plot(fpr, tpr, color=COLORS[model], label=f"{MODEL_LABELS[model]} ({values.roc_auc:.3f})")
        axes[1].plot(recall, precision, color=COLORS[model], label=f"{MODEL_LABELS[model]} ({values.average_precision:.3f})")
    axes[0].plot([0, 1], [0, 1], "--", color="#9ca3af", lw=0.8); axes[0].set(xlabel="FPR", ylabel="TPR", title="ROC")
    axes[1].set(xlabel="Recall", ylabel="Precision", title="Precision-recall")
    for ax in axes: ax.legend(frameon=False)
    savefig(fig, "ranking_curves.png")

    fig, ax = plt.subplots(figsize=(7.5, 4.2), constrained_layout=True)
    ax.plot(frontier.recall, frontier.f1, color=COLORS["v4_3_balanced"], lw=1.5, label="V4.3 post-hoc frontier")
    v3 = overall.set_index("model").loc["v3_balanced"]
    frozen = overall.set_index("model").loc["v4_3_balanced"]
    ax.scatter([frozen.recall], [frozen.f1], color="#111827", marker="x", s=70, label="V4.3 frozen")
    ax.scatter([v3.recall], [v3.f1], color=COLORS["v3_balanced"], s=50, label="V3 frozen")
    ax.set(xlabel="Recall", ylabel="F1", title="Even the post-hoc V4.3 threshold frontier does not exceed V3")
    ax.set_xlim(0.35, 1.0); ax.set_ylim(0.55, 0.9); ax.legend(frameon=False)
    savefig(fig, "threshold_frontier.png")

    pivot = detection.pivot(index="source_gse_id", columns="gene", values="detected_fraction")
    role = detection.drop_duplicates("source_gse_id").set_index("source_gse_id").role
    pivot = pivot.loc[[*role[role.eq("development")].index, *role[role.eq("lockbox")].index]]
    fig, ax = plt.subplots(figsize=(10, 4.6), constrained_layout=True)
    image = ax.imshow(100 * pivot, vmin=0, vmax=100, cmap="viridis", aspect="auto")
    ax.set_xticks(np.arange(len(pivot.columns)), pivot.columns, rotation=45, ha="right")
    labels = [f"{source} ({role[source]})" for source in pivot.index]
    ax.set_yticks(np.arange(len(pivot.index)), labels); ax.set_title("Gene detection in development and locked gdT-gold sources")
    fig.colorbar(image, ax=ax, label="Detected cells (%)", fraction=0.03)
    savefig(fig, "positive_feature_detection.png")

    top = gain.head(15).sort_values("gain")
    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    ax.barh(top.gene, top.gain, color="#7c3aed"); ax.set_xlabel("XGBoost gain"); ax.set_title("V4.3 Stage 2 is dominated by receptor magnitude")
    savefig(fig, "stage2_feature_gain.png")

    v3_boot = bootstrap.set_index("older_model").loc["v3_balanced"]
    fig, ax = plt.subplots(figsize=(6.8, 3.4), constrained_layout=True)
    metrics = ["f1", "recall", "abt_fpr", "author_nk_fpr"]
    med = [v3_boot[f"delta_{m}_median"] for m in metrics]
    low = [v3_boot[f"delta_{m}_ci_low"] for m in metrics]
    high = [v3_boot[f"delta_{m}_ci_high"] for m in metrics]
    x = np.arange(len(metrics)); ax.errorbar(x, med, yerr=[np.array(med) - low, np.array(high) - med], fmt="o", color=COLORS["v4_3_balanced"], capsize=4)
    ax.axhline(0, color="#111827", lw=0.8); ax.set_xticks(x, ["F1", "Recall", "AB FPR", "NK FPR"]); ax.set_ylabel("V4.3 minus V3"); ax.set_title("Donor-cluster bootstrap differences")
    savefig(fig, "bootstrap_deltas_vs_v3.png")


def table_html(frame: pd.DataFrame, percent_columns: list[str] | None = None, digits: int = 3) -> str:
    display = frame.copy()
    for column in percent_columns or []:
        if column in display:
            display[column] = display[column].map(lambda value: pct(float(value)))
    for column in display.select_dtypes(include=["float"]).columns:
        if column not in (percent_columns or []):
            display[column] = display[column].map(lambda value: f"{value:.{digits}f}")
    return display.to_html(index=False, escape=True, border=0, classes="data-table")


def render_report(overall: pd.DataFrame, source: pd.DataFrame, bootstrap: pd.DataFrame, decision: dict,
                  integrity: pd.DataFrame, stage: pd.DataFrame, ranking: pd.DataFrame,
                  diagnostics: pd.DataFrame, gain: pd.DataFrame, v4: dict) -> tuple[Path, Path]:
    headline = overall.copy()
    headline["model"] = headline.model.map(MODEL_LABELS)
    headline = headline[["model", "predicted_positive", "precision", "recall", "f1", "abt_fpr", "author_nk_fpr", "precision_at_1pct_prevalence"]]
    positive = source[source.truth_class.eq("gdT_gold")].pivot(index="source_gse_id", columns="model", values="positive_rate").reset_index()
    positive = positive[["source_gse_id", *MODEL_LABELS]].rename(columns=MODEL_LABELS)
    positive.columns.name = None
    negative = source[source.truth_class.eq("abT_gold")].pivot(index="source_gse_id", columns="model", values="positive_rate").reset_index()
    negative = negative[["source_gse_id", *MODEL_LABELS]].rename(columns=MODEL_LABELS)
    negative.columns.name = None
    nk = source[source.truth_class.eq("nk_lockbox")].pivot(index="source_gse_id", columns="model", values="positive_rate").reset_index()
    nk = nk[["source_gse_id", *MODEL_LABELS]].rename(columns=MODEL_LABELS)
    nk.columns.name = None
    bootstrap_display = bootstrap[["older_model", "delta_f1_median", "delta_f1_ci_low", "delta_f1_ci_high", "probability_delta_f1_gt_0", "delta_recall_median", "delta_recall_ci_low", "delta_recall_ci_high"]].copy()
    bootstrap_display["older_model"] = bootstrap_display.older_model.map(MODEL_LABELS)
    diagnostic_display = diagnostics.copy()
    stage_display = stage[["source_gse_id", "n_gdT_gold", "stage1_pass", "not_cd4_treg_excluded", "raw_stage2_pass", "receptor_rescue", "final_recall", "v3_recall", "v3_recovered_v4_missed"]]
    ranking_display = ranking.copy(); ranking_display["model"] = ranking_display.model.map(MODEL_LABELS)
    iteration = pd.DataFrame([
        {"iteration": "V4.1 GPU", "change": "Two-stage T/NK architecture", "result": "Stopped at Gate C; Stage 2 never ran"},
        {"iteration": "V4.2", "change": "Repaired TCR joins, expression-independent truth, seven-source NK negatives", "result": "Failed development transfer; Malte recall 36.5%"},
        {"iteration": "V4.3 initial", "change": "Receptor-isolated Stage 2; removed Stage 1 probability", "result": "Passed development, but KLRD1 remained in Stage 1"},
        {"iteration": "V4.3 final", "change": "Removed KLRD1/shared cytotoxic genes from both stages", "result": "Passed development; failed locked external superiority"},
    ])
    css = """
    :root{--ink:#17212b;--muted:#5e6b76;--line:#d8dee5;--red:#b42318;--blue:#176b87;--pale:#f4f7f8}
    *{box-sizing:border-box}body{font-family:Arial,Helvetica,sans-serif;color:var(--ink);margin:0;background:#fff;line-height:1.42}
    main{max-width:1180px;margin:0 auto;padding:28px 34px 60px}header{border-bottom:4px solid var(--red);padding-bottom:16px;margin-bottom:22px}
    h1{font-size:30px;margin:0 0 6px;letter-spacing:0}h2{font-size:20px;margin:28px 0 10px;border-bottom:1px solid var(--line);padding-bottom:5px}h3{font-size:15px;margin:18px 0 7px}
    p,li{font-size:13px}.subtitle{color:var(--muted);font-size:15px}.decision{border-left:6px solid var(--red);background:#fff3f1;padding:14px 17px;margin:18px 0;font-size:14px}
    .method{border-left:5px solid var(--blue);background:#eff7fa;padding:12px 16px;margin:15px 0}.grid{display:grid;grid-template-columns:1fr 1fr;gap:16px}.metric-grid{display:grid;grid-template-columns:repeat(4,1fr);gap:9px;margin:14px 0}
    .metric{border:1px solid var(--line);padding:10px;background:var(--pale)}.metric b{font-size:20px;display:block}.metric span{font-size:11px;color:var(--muted)}
    figure{margin:10px 0 18px;break-inside:avoid}figure img{width:100%;border:1px solid var(--line)}figcaption{font-size:10px;color:var(--muted);margin-top:4px}
    .table-wrap{overflow-x:auto;margin:8px 0 18px}.data-table{width:100%;border-collapse:collapse;font-size:10px;table-layout:auto}.data-table th,.data-table td{border:1px solid var(--line);padding:4px 5px;vertical-align:top;overflow-wrap:anywhere}.data-table th{background:#e9eef1;text-align:left}
    code{font-family:Consolas,monospace;font-size:11px;background:#eef1f3;padding:1px 4px}.note{font-size:11px;color:var(--muted)}.page-break{break-before:page}
    @media print{@page{size:A4 landscape;margin:9mm}main{max-width:none;padding:0}h2{break-after:avoid}.grid{gap:10px}.data-table{font-size:8px}.metric-grid{grid-template-columns:repeat(4,1fr)}figure img{max-height:158mm;object-fit:contain}.page-break{break-before:page}}
    """
    asset = lambda name: f"assets/{name}"
    report = f"""<!doctype html><html><head><meta charset='utf-8'><title>gdTAI V4 Final Evaluation</title><style>{css}</style></head><body><main>
    <header><h1>gdTAI V4 Final Evaluation</h1><div class='subtitle'>Development-frozen training, common external lockbox comparison, and final no-promotion decision</div></header>
    <div class='decision'><b>Decision: do not promote gdTAI V4.3.</b> V4.3 met its false-positive limits but failed external recall and did not outperform V3. The current V3 balanced model remains the default. No V4 artifact was released, pushed, or applied to the 5.9-million-cell atlas.</div>
    <div class='metric-grid'><div class='metric'><b>{pct(float(overall.set_index('model').loc['v4_3_balanced','recall']))}</b><span>V4.3 external gdT recall</span></div><div class='metric'><b>{pct(float(overall.set_index('model').loc['v4_3_balanced','abt_fpr']))}</b><span>V4.3 paired-alpha-beta FPR</span></div><div class='metric'><b>{pct(float(overall.set_index('model').loc['v4_3_balanced','author_nk_fpr']))}</b><span>V4.3 author-NK FPR</span></div><div class='metric'><b>{float(overall.set_index('model').loc['v4_3_balanced','f1']):.3f}</b><span>V4.3 external F1</span></div></div>
    <h2>1. Evaluation Design</h2><p>The final V4.3 candidate, all thresholds, receptor-rescue rules, and the feature contract were frozen using development cells only. The external lockbox was checksum-bound before scoring and contained <b>335,479 cells</b>: 41,931 gdT gold, 285,524 paired-alpha-beta gold, and 8,024 independent author-NK negatives. V4.3, V3, and both V2 modes were recomputed from the same raw-count matrix and identical cells.</p>
    <div class='method'><b>Ground truth.</b> gdT positives were sorted gdT cohorts or strict paired productive TRG/TRD cells without alpha-beta receptor evidence. Alpha-beta negatives had productive alpha/beta receptor evidence without gamma/delta evidence. Author-NK negatives had NK annotation and no productive TCR. Where chain UMI support existed, one-UMI calls were not accepted as gold; datasets without UMI fields retained published productive calls.</div>
    <p><b>V4.3 architecture:</b> a 36-gene high-recall T-lineage support gate followed by a 158-gene receptor classifier. Stage 2 uses 153 individual TCR genes plus <code>CD3D/E/G</code>, <code>CD4</code>, and <code>FOXP3</code>. <code>KLRD1</code>, <code>NKG7</code>, <code>GNLY</code>, <code>FCER1G</code>, <code>TYROBP</code>, and the shared cytotoxic program are absent from both stages. Stage 1 probability is not supplied to Stage 2.</p>
    <h3>Implementation integrity</h3><div class='table-wrap'>{table_html(integrity)}</div>
    <h2>2. Headline Performance</h2><div class='table-wrap'>{table_html(headline, ['precision','recall','f1','abt_fpr','author_nk_fpr','precision_at_1pct_prevalence'])}</div>
    <figure><img src='{asset('overall_performance.png')}'><figcaption>V4.3 gained NK specificity relative to V3 but lost too much gdT recall. Precision at 1% prevalence is included because the atlas is strongly imbalanced.</figcaption></figure>
    <div class='grid'><figure><img src='{asset('ranking_curves.png')}'><figcaption>V4.3 has lower ROC-AUC and average precision than V3 and V2; this is a ranking deficit, not only threshold calibration.</figcaption></figure><figure><img src='{asset('bootstrap_deltas_vs_v3.png')}'><figcaption>Donor-cluster bootstrap. The V4.3-minus-V3 F1 and recall intervals remain entirely below zero.</figcaption></figure></div>
    <h3>Ranking metrics</h3><div class='table-wrap'>{table_html(ranking_display)}</div>
    <h3>Donor-cluster bootstrap</h3><div class='table-wrap'>{table_html(bootstrap_display)}</div>
    <div class='page-break'></div><h2>3. Dataset-by-Dataset Results</h2><figure><img src='{asset('positive_source_recall.png')}'><figcaption>V4.3 under-recovers every locked positive source. GDTlung is retained as a documented low-quality stress cohort, but BALF_BLOOD_COPD and GDT_2020 also show large deficits.</figcaption></figure>
    <div class='table-wrap'>{table_html(positive, list(MODEL_LABELS.values()))}</div>
    <div class='grid'><figure><img src='{asset('abt_fpr_by_source.png')}'><figcaption>Paired-alpha-beta false-positive rates vary by external source. V4.3 is not uniformly more specific than V3.</figcaption></figure><figure><img src='{asset('nk_fpr_by_source.png')}'><figcaption>V4.3 reduces author-NK false positives relative to V3 and V2, but purity does not compensate for the recall collapse.</figcaption></figure></div>
    <h3>Paired-alpha-beta false-positive rates</h3><div class='table-wrap'>{table_html(negative, list(MODEL_LABELS.values()))}</div><h3>Author-NK false-positive rates</h3><div class='table-wrap'>{table_html(nk, list(MODEL_LABELS.values()))}</div>
    <div class='page-break'></div><h2>4. Why V4.3 Failed</h2><figure><img src='{asset('stage_decomposition.png')}'><figcaption>Stage 1 passes 98-100% of locked gdT gold and CD4/Treg exclusions affect no more than 1.4%. Most false negatives arise from the Stage 2 receptor classifier.</figcaption></figure>
    <div class='table-wrap'>{table_html(stage_display, ['stage1_pass','not_cd4_treg_excluded','raw_stage2_pass','receptor_rescue','final_recall','v3_recall'])}</div>
    <div class='grid'><figure><img src='{asset('stage2_feature_gain.png')}'><figcaption>TRDC and TRDV magnitude dominate the nonlinear Stage 2 model despite the broader individual-receptor panel.</figcaption></figure><figure><img src='{asset('positive_feature_detection.png')}'><figcaption>Positive sources differ substantially in receptor detection. GDTlung has especially weak TRDC/CD3 support, while BALF and GDT_2020 also differ from development cohorts.</figcaption></figure></div>
    <p>The final Stage 2 threshold was <code>{float(v4['stage2_threshold']):.6f}</code>, selected from grouped development OOF predictions. Lowering it after seeing the lockbox would be leakage. A diagnostic-only threshold sweep nevertheless shows whether calibration alone could solve the problem.</p>
    <figure><img src='{asset('threshold_frontier.png')}'><figcaption>The best post-hoc V4.3 F1 is still below frozen V3. Matching V3 recall requires substantially more alpha-beta false positives, so threshold adjustment cannot establish superiority.</figcaption></figure>
    <div class='table-wrap'>{table_html(diagnostic_display, ['f1','recall','precision','abt_fpr','author_nk_fpr'])}</div>
    <h2>5. Iteration Record</h2><div class='table-wrap'>{table_html(iteration)}</div>
    <p>V4 development addressed the major known risks in sequence: weak NK labels, unsafe TCR joins, source leakage, shared cytotoxic genes, continuous Stage 1 leakage, and source-heldout evaluation. V4.3 is the first iteration to reach a clean common lockbox comparison, and it fails there.</p>
    <h2>6. Final Scientific Conclusion</h2><div class='decision'><b>With the current datasets, no defensible V4 improvement over V3 was demonstrated.</b> Training another candidate after inspecting this lockbox would convert the external test into a tuning set. The ranking deficit also shows that a threshold-only revision is insufficient. Therefore V4 development should stop here rather than optimize against the consumed validation cohorts.</div>
    <p>This conclusion is scoped to the available evidence, not a claim that gdT classification can never improve. A future V4 attempt requires genuinely new independent positive and NK/alpha-beta cohorts, or a prespecified domain-generalization design with a new untouched final test. Until then:</p><ul><li>retain V3 balanced as the current default;</li><li>retain V2 high-purity as a conservative annotation-assisted fallback;</li><li>keep all V4.3 artifacts diagnostic and unpromoted;</li><li>do not use lockbox-derived thresholds or feature changes in a release model.</li></ul>
    <p class='note'>Core hashes: V4 contract <code>{html.escape(str(json.load(open(ROOT / 'Integrated_dataset/logs/gdT_prediction/gdtai_v4_3_common_lockbox/evaluation_summary.json'))['v4_contract_sha256']))}</code>; lockbox <code>{html.escape(str(json.load(open(ROOT / 'Integrated_dataset/logs/gdT_prediction/gdtai_v4_3_common_lockbox/evaluation_summary.json'))['lockbox_sha256']))}</code>. Model promotion status: false.</p>
    </main></body></html>"""
    html_path = STATIC_DIR / "index.html"
    pdf_path = STATIC_DIR / "gdtai_v4_3_final_evaluation_report.pdf"
    html_path.write_text(report, encoding="utf-8")
    try:
        with tempfile.TemporaryDirectory(prefix="gdtai-v4-report-") as profile:
            subprocess.run([
                "google-chrome", "--headless", "--no-sandbox", "--disable-gpu", "--disable-crash-reporter",
                f"--user-data-dir={profile}", "--no-pdf-header-footer",
                f"--print-to-pdf={pdf_path}", html_path.resolve().as_uri(),
            ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=240)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        subprocess.run(
            ["/home/anaconda3/bin/weasyprint", str(html_path), str(pdf_path)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=240,
        )
    return html_path, pdf_path


def main() -> None:
    setup()
    predictions = pd.read_parquet(PREDICTIONS)
    overall = pd.read_csv(OVERALL)
    source = pd.read_csv(PER_SOURCE)
    bootstrap = pd.read_csv(BOOTSTRAP)
    decision = json.loads(DECISION.read_text())
    v4 = json.loads(V4_CONTRACT.read_text())
    integrity = implementation_integrity(v4)
    stage = stage_decomposition(predictions, v4)
    ranking = ranking_metrics(predictions)
    diagnostics, frontier = threshold_frontier(predictions, v4, overall)
    detection = feature_detection(v4)
    gain = model_gain(v4)
    make_figures(predictions, overall, source, bootstrap, stage, ranking, diagnostics, frontier, detection, gain)
    html_path, pdf_path = render_report(overall, source, bootstrap, decision, integrity, stage, ranking, diagnostics, gain, v4)
    conclusion = f"""# gdTAI V4 final conclusion

- decision: `{decision['status']}`
- eligible for promotion review: `{decision['eligible_for_promotion_review']}`
- V4.3 recall / precision / F1: `{overall.set_index('model').loc['v4_3_balanced','recall']:.4f}` / `{overall.set_index('model').loc['v4_3_balanced','precision']:.4f}` / `{overall.set_index('model').loc['v4_3_balanced','f1']:.4f}`
- V3 recall / precision / F1: `{overall.set_index('model').loc['v3_balanced','recall']:.4f}` / `{overall.set_index('model').loc['v3_balanced','precision']:.4f}` / `{overall.set_index('model').loc['v3_balanced','f1']:.4f}`
- V4.3 paired-alpha-beta FPR: `{overall.set_index('model').loc['v4_3_balanced','abt_fpr']:.6f}`
- V4.3 author-NK FPR: `{overall.set_index('model').loc['v4_3_balanced','author_nk_fpr']:.6f}`
- post-hoc maximum V4.3 F1 (diagnostic only): `{diagnostics.set_index('diagnostic').loc['posthoc_maximum_f1','f1']:.4f}`
- lockbox retuning: `forbidden`
- model promoted: `false`
- HTML: `{html_path}`
- PDF: `{pdf_path}`

V4.3 does not surpass V3. The effective feature matrix is exact, Stage 1 and exclusions are not the cause, and even the post-hoc V4 threshold frontier does not exceed frozen V3. With the final lockbox now consumed, further model selection on these cohorts would be overfitting. V3 remains the current default.
"""
    summary_path = LOG_DIR / "gdtai_v4_final_conclusion.md"
    summary_path.write_text(conclusion, encoding="utf-8")
    result = {
        "status": "PASS_FINAL_REPORT_NO_PROMOTION",
        "html": str(html_path),
        "pdf": str(pdf_path),
        "summary": str(summary_path),
        "model_promoted": False,
        "v4_surpasses_older": False,
        "further_lockbox_tuning_allowed": False,
    }
    (LOG_DIR / "report_summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
