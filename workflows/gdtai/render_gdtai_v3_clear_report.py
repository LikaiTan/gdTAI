#!/usr/bin/env python3
"""Render a reader-facing gdTAI v3 decision report from saved result tables."""

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
import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT_PREFIX = "gdtai_v3_trdc_nk_guard"
TABLE_DIR = ROOT / "Integrated_dataset" / "tables" / "gdT_prediction" / OUT_PREFIX
LOG_DIR = ROOT / "Integrated_dataset" / "logs" / "gdT_prediction" / OUT_PREFIX
MODEL_DIR = ROOT / "Integrated_dataset" / "models" / "gdT_prediction_classifier" / OUT_PREFIX
STATIC_DIR = ROOT / "gdT_prediction"
ASSET_DIR = STATIC_DIR / "assets" / OUT_PREFIX

SUMMARY_JSON = LOG_DIR / "gdtai_v3_trdc_nk_guard_summary.json"
REPORT_MD = LOG_DIR / "gdtai_v3_trdc_nk_guard_report.md"
REPORT_HTML = STATIC_DIR / "gdtai_v3_trdc_nk_guard_report.html"
REPORT_PDF = STATIC_DIR / "gdtai_v3_trdc_nk_guard_report.pdf"


SHORT_NAMES = {
    "v2_high_f1": "v2 high-F1",
    "v2_high_purity": "v2 high-purity",
    "original_TRD_minus_TRAB": "TRD-TRAB score",
    "v3_round03_hist_gradient_target_recall80_fp05": "v3 hist-gradient",
    "v3_round08_hist_gradient_fixed_0p64": "v3 hist-gradient relaxed",
    "v3_round14_v2_score_trdc_gate_fixed_0p936": "v3 selected",
}


def fmt_int(value: float | int) -> str:
    return f"{int(round(float(value))):,}"


def fmt_pct(value: float, digits: int = 1) -> str:
    return f"{100.0 * float(value):.{digits}f}%"


def fmt_dec(value: float, digits: int = 3) -> str:
    return f"{float(value):.{digits}f}"


def read_csv(name: str) -> pd.DataFrame:
    path = TABLE_DIR / name
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def one_row(df: pd.DataFrame, column: str, value: str) -> pd.Series:
    match = df.loc[df[column].astype(str) == value]
    if match.empty:
        raise KeyError(f"{value!r} not found in column {column!r}")
    return match.iloc[0]


def strategy_label(strategy: str) -> str:
    return SHORT_NAMES.get(strategy, strategy.replace("_", " "))


def rel_asset(path: Path) -> str:
    return path.relative_to(STATIC_DIR).as_posix()


def small_markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    out = ["| " + " | ".join(headers) + " |"]
    out.append("| " + " | ".join(["---"] * len(headers)) + " |")
    for row in rows:
        out.append("| " + " | ".join(row) + " |")
    return "\n".join(out)


def html_table(headers: list[str], rows: list[list[str]]) -> str:
    head = "".join(f"<th>{html.escape(h)}</th>" for h in headers)
    body_rows = []
    for row in rows:
        body_rows.append("<tr>" + "".join(f"<td>{html.escape(str(cell))}</td>" for cell in row) + "</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def save_metric_bars(metrics: pd.DataFrame, selected: str) -> Path:
    strategies = ["v2_high_purity", selected, "v2_high_f1", "original_TRD_minus_TRAB"]
    cols = ["precision", "recall", "f1"]
    sub = metrics.loc[metrics["strategy"].isin(strategies)].copy()
    sub["strategy"] = pd.Categorical(sub["strategy"], strategies, ordered=True)
    sub = sub.sort_values("strategy")
    x = np.arange(len(cols))
    width = 0.18
    colors = ["#059669", "#2563eb", "#f59e0b", "#6b7280"]

    fig, ax = plt.subplots(figsize=(8.6, 4.7))
    for i, (_, row) in enumerate(sub.iterrows()):
        values = [float(row[col]) for col in cols]
        offset = (i - 1.5) * width
        bars = ax.bar(x + offset, values, width, label=strategy_label(str(row["strategy"])), color=colors[i])
        for bar, value in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width() / 2, value + 0.008, fmt_dec(value, 2), ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0.74, 1.0)
    ax.set_xticks(x, ["Precision", "Recall", "F1"])
    ax.set_ylabel("External primary metric")
    ax.set_title("External primary performance")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, ncol=2, loc="lower left")
    fig.tight_layout()
    out = ASSET_DIR / "decision_external_metric_bars.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def save_fp_group_bars(fp: pd.DataFrame, selected: str) -> Path:
    groups = [
        "NK_cells",
        "paired_TCRAB_cells",
        "TRDC_plus_TRDV_minus_cells",
        "CD4_Treg_warning_cells",
        "B_myeloid_cells",
    ]
    labels = ["NK", "paired TCRAB", "TRDC+/TRDV-", "CD4/Treg warn", "B/myeloid"]
    strategies = ["v2_high_purity", selected, "original_TRD_minus_TRAB"]
    pivot = fp.pivot_table(index="false_positive_group", columns="strategy", values="false_positive_cells", aggfunc="sum").reindex(groups)
    y = np.arange(len(groups))
    width = 0.24
    colors = ["#059669", "#2563eb", "#6b7280"]

    fig, ax = plt.subplots(figsize=(8.6, 4.9))
    for i, strategy in enumerate(strategies):
        values = pivot[strategy].fillna(0).to_numpy(dtype=float)
        bars = ax.barh(y + (i - 1) * width, values, width, label=strategy_label(strategy), color=colors[i])
        for bar, value in zip(bars, values):
            ax.text(value + 1.5, bar.get_y() + bar.get_height() / 2, fmt_int(value), va="center", fontsize=8)
    ax.set_yticks(y, labels)
    ax.set_xlabel("External false-positive cells")
    ax.set_title("Where false positives remain")
    ax.grid(axis="x", alpha=0.25)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    out = ASSET_DIR / "decision_external_fp_groups.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def save_recall_group_bars(recall: pd.DataFrame, selected: str) -> Path:
    groups = [
        "all_real_gdt",
        "TRDV_positive_gdT",
        "TRDC_plus_TRDV_minus_gdT",
        "cytotoxic_NK_marker_high_gdT",
        "Vd2_gdT_cells",
    ]
    labels = ["all real gdT", "TRDV+ gdT", "TRDC+/TRDV- gdT", "cytotoxic gdT", "Vd2 gdT"]
    strategies = ["v2_high_purity", selected, "v2_high_f1"]
    pivot = recall.pivot_table(index="recall_group", columns="strategy", values="recall", aggfunc="first").reindex(groups)
    y = np.arange(len(groups))
    width = 0.24
    colors = ["#059669", "#2563eb", "#f59e0b"]

    fig, ax = plt.subplots(figsize=(8.6, 4.9))
    for i, strategy in enumerate(strategies):
        values = pivot[strategy].fillna(0).to_numpy(dtype=float)
        bars = ax.barh(y + (i - 1) * width, values, width, label=strategy_label(strategy), color=colors[i])
        for bar, value in zip(bars, values):
            ax.text(value + 0.01, bar.get_y() + bar.get_height() / 2, fmt_pct(value, 0), va="center", fontsize=8)
    ax.set_xlim(0.35, 1.02)
    ax.set_yticks(y, labels)
    ax.set_xlabel("External recall")
    ax.set_title("Recall by gdT subgroup")
    ax.grid(axis="x", alpha=0.25)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    out = ASSET_DIR / "decision_external_recall_groups.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def save_full_atlas_tradeoff(target: pd.DataFrame, selected: str, recall_target: float, fp_target: float) -> Path:
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    x = 100.0 * target["estimated_fp_fraction_of_predictions"].to_numpy(dtype=float)
    y = 100.0 * target["full_primary_recall"].to_numpy(dtype=float)
    met = target["target_met_full_atlas"].astype(bool).to_numpy()
    ax.scatter(x[~met], y[~met], s=55, color="#9ca3af", label="candidate failed target")
    ax.scatter(x[met], y[met], s=70, color="#059669", label="candidate met target")
    selected_mask = target["strategy"].astype(str).eq(selected).to_numpy()
    ax.scatter(x[selected_mask], y[selected_mask], s=150, facecolors="none", edgecolors="#dc2626", linewidths=2.5, label="selected v3")
    ax.axvline(100.0 * fp_target, color="#dc2626", linestyle="--", linewidth=1.2)
    ax.axhline(100.0 * recall_target, color="#2563eb", linestyle="--", linewidth=1.2)
    for _, row in target.loc[selected_mask | target["target_met_full_atlas"].astype(bool)].iterrows():
        label = strategy_label(str(row["strategy"]))
        if str(row["strategy"]) == selected:
            label = "round14 selected"
        elif "round15" in str(row["strategy"]):
            label = "round15"
        elif "round16" in str(row["strategy"]):
            label = "round16"
        else:
            continue
        ax.annotate(
            label,
            (100.0 * float(row["estimated_fp_fraction_of_predictions"]), 100.0 * float(row["full_primary_recall"])),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )
    ax.set_xlabel("Estimated false positives / predictions")
    ax.set_ylabel("Full-atlas primary-gold recall")
    ax.set_title("Full-atlas target trade-off")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    out = ASSET_DIR / "decision_full_atlas_tradeoff.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def save_gate_status(acceptance: pd.Series) -> Path:
    gate_cols = [
        ("external_NK_FP_lower_than_v2_high_purity", "NK FP lower than v2 HP"),
        ("external_TRDC_plus_TRDV_minus_burden_lower_than_v2_high_purity", "TRDC+/TRDV- burden lower"),
        ("external_F1_not_worse_than_v2_high_purity_by_gt_0p01", "F1 not worse by >0.01"),
        ("cytotoxic_gdT_recall_not_substantially_degraded", "cytotoxic recall preserved"),
        ("paired_TCRAB_FP_not_above_v2_high_purity", "paired TCRAB FP not higher"),
    ]
    labels = [label for _, label in gate_cols]
    values = [bool(acceptance[col]) for col, _ in gate_cols]
    colors = ["#059669" if value else "#dc2626" for value in values]
    y = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(8.6, 3.7))
    ax.barh(y, [1] * len(labels), color=colors)
    for yi, value in zip(y, values):
        ax.text(0.5, yi, "PASS" if value else "FAIL", ha="center", va="center", color="white", weight="bold")
    ax.set_yticks(y, labels)
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.invert_yaxis()
    ax.set_title("Promotion gate checklist")
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    out = ASSET_DIR / "decision_promotion_gates.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def save_confusion_matrix(selected_metrics: pd.Series) -> Path:
    matrix = np.array(
        [
            [int(selected_metrics["tp"]), int(selected_metrics["fn"])],
            [int(selected_metrics["fp"]), int(selected_metrics["tn"])],
        ],
        dtype=float,
    )
    labels = [["TP", "FN"], ["FP", "TN"]]
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    ax.imshow(np.log10(matrix + 1), cmap="Blues")
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f"{labels[i][j]}\n{fmt_int(matrix[i, j])}", ha="center", va="center", fontsize=13, weight="bold")
    ax.set_xticks([0, 1], ["predicted gdT", "predicted non-gdT"])
    ax.set_yticks([0, 1], ["true gdT", "true negative"])
    ax.set_title("External confusion matrix for selected v3")
    fig.tight_layout()
    out = ASSET_DIR / "decision_external_confusion_matrix.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def save_figures(
    metrics: pd.DataFrame,
    fp: pd.DataFrame,
    recall: pd.DataFrame,
    target: pd.DataFrame,
    acceptance: pd.Series,
    selected: str,
    recall_target: float,
    fp_target: float,
) -> list[Path]:
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    selected_metrics = one_row(metrics, "strategy", selected)
    return [
        save_metric_bars(metrics, selected),
        save_full_atlas_tradeoff(target, selected, recall_target, fp_target),
        save_gate_status(acceptance),
        save_fp_group_bars(fp, selected),
        save_recall_group_bars(recall, selected),
        save_confusion_matrix(selected_metrics),
    ]


def build_report(no_pdf: bool) -> None:
    with SUMMARY_JSON.open("r", encoding="utf-8") as handle:
        summary = json.load(handle)
    manifest_path = MODEL_DIR / "model_manifest.json"
    with manifest_path.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)

    selected = str(summary["selected_candidate"])
    metrics = read_csv("external_primary_metrics.csv")
    fp = read_csv("external_false_positive_groups.csv")
    recall = read_csv("external_recall_groups.csv")
    acceptance_df = read_csv("promotion_acceptance_gates.csv")
    target = read_csv("full_atlas_target_selection.csv")
    overall = read_csv("full_atlas_prediction_overall.csv")
    selected_strategy = read_csv("full_atlas_selected_strategy.csv")
    tune = read_csv("internal_tune_candidate_metrics.csv")
    split = read_csv("split_overall.csv")
    hard_filter = read_csv("hard_negative_expression_filter_summary.csv")
    weak_filter = read_csv("weak_trdc_hard_negative_filter_summary.csv")
    hard_by_source = read_csv("hard_negative_expression_filter_by_source.csv")

    selected_metrics = one_row(metrics, "strategy", selected)
    v2hp_metrics = one_row(metrics, "strategy", "v2_high_purity")
    v2hf_metrics = one_row(metrics, "strategy", "v2_high_f1")
    trd_metrics = one_row(metrics, "strategy", "original_TRD_minus_TRAB")
    selected_full = one_row(overall, "strategy", selected)
    v2hp_full = one_row(overall, "strategy", "v2_high_purity")
    acceptance = acceptance_df.iloc[0]
    selection = selected_strategy.iloc[0]

    recall_target = float(selection["target_recall"]) + float(selection["target_recall_margin"])
    fp_target = float(selection["target_max_estimated_fp_fraction"])

    fp_pivot = fp.pivot_table(index="strategy", columns="false_positive_group", values="false_positive_cells", aggfunc="sum").fillna(0)
    recall_pivot = recall.pivot_table(index="strategy", columns="recall_group", values="recall", aggfunc="first").fillna(np.nan)

    nk_sources = hard_by_source.loc[hard_by_source["retained_NK_TRDCpos_TRDVneg"] > 0, "source_gse_id"].nunique()
    train_nk = int(hard_filter.loc[hard_filter["split"] == "train_hard_negative", "retained_NK_TRDCpos_TRDVneg"].iloc[0])
    tune_nk = int(hard_filter.loc[hard_filter["split"] == "tune_hard_negative", "retained_NK_TRDCpos_TRDVneg"].iloc[0])
    combined_nk = int(hard_filter["retained_NK_TRDCpos_TRDVneg"].sum())
    excluded_nk = int(hard_filter["excluded_NK_not_TRDCpos_TRDVneg"].sum())
    excluded_overlap = int(hard_filter["excluded_NK_TCRAB_overlap"].sum())
    weak_overlap_excluded = int(weak_filter["weak_trdc_excluded_NK_TCRAB_overlap"].iloc[0])

    validation_row = one_row(split, "split", "validation_combined_non_GSE144469")
    gse144469_row = one_row(split, "split", "GSE144469_primary_in_train_tune")

    target_met = target.loc[target["target_met_full_atlas"].astype(bool)].copy()
    selected_target = one_row(target, "strategy", selected)
    best_external_v3 = metrics.loc[metrics["strategy"].str.startswith("v3_")].sort_values("f1", ascending=False).iloc[0]
    previous_internal = str(selection["previous_internal_selected_strategy"])
    previous_target = one_row(target, "strategy", previous_internal)

    figures = save_figures(metrics, fp, recall, target, acceptance, selected, recall_target, fp_target)

    gate_rows = [
        [
            "NK false positives lower than v2 high-purity",
            "PASS" if bool(acceptance["external_NK_FP_lower_than_v2_high_purity"]) else "FAIL",
            f"v3 {fmt_int(acceptance['external_NK_FP'])} vs v2 high-purity {fmt_int(acceptance['v2_high_purity_NK_FP'])}",
        ],
        [
            "TRDC+/TRDV- false-positive burden lower",
            "PASS" if bool(acceptance["external_TRDC_plus_TRDV_minus_burden_lower_than_v2_high_purity"]) else "FAIL",
            f"v3 {fmt_int(acceptance['external_TRDC_plus_TRDV_minus_FP'])} vs v2 high-purity {fmt_int(acceptance['v2_high_purity_TRDC_plus_TRDV_minus_FP'])}",
        ],
        [
            "External F1 not worse by more than 0.01",
            "PASS" if bool(acceptance["external_F1_not_worse_than_v2_high_purity_by_gt_0p01"]) else "FAIL",
            f"v3 {fmt_dec(acceptance['external_f1'], 4)} vs v2 high-purity {fmt_dec(acceptance['v2_high_purity_external_f1'], 4)}",
        ],
        [
            "Cytotoxic gdT recall preserved",
            "PASS" if bool(acceptance["cytotoxic_gdT_recall_not_substantially_degraded"]) else "FAIL",
            f"v3 {fmt_pct(recall_pivot.loc[selected, 'cytotoxic_NK_marker_high_gdT'])} vs v2 high-purity {fmt_pct(recall_pivot.loc['v2_high_purity', 'cytotoxic_NK_marker_high_gdT'])}",
        ],
        [
            "Paired-TCRAB false positives not higher",
            "PASS" if bool(acceptance["paired_TCRAB_FP_not_above_v2_high_purity"]) else "FAIL",
            f"v3 {fmt_int(acceptance['external_paired_TCRAB_FP'])} vs v2 high-purity {fmt_int(acceptance['v2_high_purity_paired_TCRAB_FP'])}",
        ],
    ]

    metric_rows = [
        ["External precision", fmt_pct(selected_metrics["precision"], 1), fmt_pct(v2hp_metrics["precision"], 1), "v3 higher"],
        ["External recall", fmt_pct(selected_metrics["recall"], 1), fmt_pct(v2hp_metrics["recall"], 1), "v3 slightly lower"],
        ["External F1", fmt_dec(selected_metrics["f1"], 4), fmt_dec(v2hp_metrics["f1"], 4), "near-tie"],
        ["External FP / predictions", fmt_pct(selected_metrics["fp_fraction_of_predictions"], 1), fmt_pct(v2hp_metrics["fp_fraction_of_predictions"], 1), "v3 lower"],
        ["External NK FP", fmt_int(acceptance["external_NK_FP"]), fmt_int(acceptance["v2_high_purity_NK_FP"]), "v3 worse"],
    ]

    full_rows = [
        ["Predicted putative gdT", fmt_int(selected_full["predicted_putative_gdT"]), fmt_int(v2hp_full["predicted_putative_gdT"])],
        ["Primary-gold recall", fmt_pct(selected_full["full_primary_recall"], 2), fmt_pct(v2hp_full["full_primary_recall"], 2)],
        ["Primary-gold precision", fmt_pct(selected_full["full_primary_precision"], 2), fmt_pct(v2hp_full["full_primary_precision"], 2)],
        ["Estimated total abT FP", fmt_int(selected_full["estimated_total_abT_fp"]), fmt_int(v2hp_full["estimated_total_abT_fp"])],
        ["Estimated FP / predictions", fmt_pct(selected_full["estimated_fp_fraction_of_predictions"], 2), fmt_pct(v2hp_full["estimated_fp_fraction_of_predictions"], 2)],
        ["Predicted paired TCRAB", fmt_int(selected_full["predicted_paired_TCRAB"]), fmt_int(v2hp_full["predicted_paired_TCRAB"])],
        ["Predicted TRDC+/TRDV-", fmt_int(selected_full["predicted_TRDC_plus_TRDV_minus"]), fmt_int(v2hp_full["predicted_TRDC_plus_TRDV_minus"])],
    ]

    md_sections = [
        "# gdTAI v3 TRDC/NK Guard Decision Report",
        "",
        "## Executive decision",
        "",
        "**Decision: do not promote this v3 candidate as the default release model yet.**",
        "",
        (
            f"The selected candidate `{selected}` met the user target on the full atlas "
            f"({fmt_pct(selected_full['full_primary_recall'], 2)} primary-gold recall and "
            f"{fmt_pct(selected_full['estimated_fp_fraction_of_predictions'], 2)} estimated FP/prediction), "
            f"but it failed the explicit promotion gate for NK false positives on the independent external test "
            f"({fmt_int(acceptance['external_NK_FP'])} NK FP vs {fmt_int(acceptance['v2_high_purity_NK_FP'])} for v2 high-purity)."
        ),
        "",
        "## What changed",
        "",
        (
            f"- `GSE144469` was moved into train/tune: {fmt_int(gse144469_row['n_cells'])} cells, "
            f"including {fmt_int(gse144469_row['gdT_gold'])} gdT positives and {fmt_int(gse144469_row['abT_gold'])} abT negatives."
        ),
        (
            f"- Other validation cohorts were retained: {fmt_int(validation_row['n_cells'])} held-out non-GSE144469 cells."
        ),
        "- The independent external H5AD was used only after model and threshold selection.",
        "- External inference used `layers[\"counts\"]` with log1p(counts per 10,000), not normalized `X`.",
        (
            f"- Hard NK negatives followed the requested rule: {fmt_int(combined_nk)} NK `TRDC+TRDV-` cells "
            f"from {nk_sources} datasets were retained; {fmt_int(excluded_nk)} NK cells not matching that expression rule were excluded. "
            f"Also, {fmt_int(weak_overlap_excluded)} expression-passing weak-TRDC NK+TCRAB-overlap cells were excluded."
        ),
        "",
        "**Step conclusion:** the evaluation design is now aligned with the requested split policy, and the new external H5AD remains independent.",
        "",
        "## Algorithm",
        "",
        "- The selected model is not XGBoost.",
        "- It is a conditional gate around the gdTAI v2 transcriptome score.",
        "- It uses count-derived gene features plus engineered TCR/CD3/NK features.",
        "- It does not require TCR-seq metadata at inference time.",
        "- The operating threshold is `0.936`.",
        "",
        "**Step conclusion:** this is a conservative v2-score-plus-guard candidate, not a fully new tree model.",
        "",
        "## Candidate selection",
        "",
        (
            f"The initial internal winner was `{previous_internal}`, but it had only "
            f"{fmt_pct(previous_target['full_primary_recall'], 2)} full-atlas primary-gold recall. "
            f"The final selected round had {fmt_pct(selected_target['full_primary_recall'], 2)} recall and "
            f"{fmt_pct(selected_target['estimated_fp_fraction_of_predictions'], 2)} estimated FP/prediction, passing the "
            f"{fmt_pct(recall_target, 1)} recall-with-margin and {fmt_pct(fp_target, 1)} estimated-FP targets."
        ),
        "",
        (
            f"External testing was not used for model selection. For example, `{best_external_v3['strategy']}` had the best v3 external F1 "
            f"({fmt_dec(best_external_v3['f1'], 4)}), but it was not selected because the final decision rule prioritized the full-atlas target."
        ),
        "",
        "**Step conclusion:** round14 is the selected candidate because it is the first practical point satisfying the atlas-level recall and FP estimate target.",
        "",
        "## External performance",
        "",
        small_markdown_table(["Metric", "v3 selected", "v2 high-purity", "Interpretation"], metric_rows),
        "",
        (
            f"Compared with v2 high-purity, v3 reduced total external false positives from "
            f"{fmt_int(v2hp_metrics['fp'])} to {fmt_int(selected_metrics['fp'])} and reduced paired-TCRAB false positives "
            f"from {fmt_int(acceptance['v2_high_purity_paired_TCRAB_FP'])} to {fmt_int(acceptance['external_paired_TCRAB_FP'])}. "
            f"However, NK false positives increased from {fmt_int(acceptance['v2_high_purity_NK_FP'])} to {fmt_int(acceptance['external_NK_FP'])}."
        ),
        "",
        "**Step conclusion:** the candidate is competitive overall, but it does not solve the NK-specific failure mode strongly enough.",
        "",
        "## Full-atlas application",
        "",
        small_markdown_table(["Full-atlas metric", "v3 selected", "v2 high-purity"], full_rows),
        "",
        (
            "The full-atlas recall is measured on primary gold labels only. The estimated false-positive count extrapolates observed "
            "paired-TCRAB false-positive behavior from TCR-seq sources to sources without TCR-seq, so it is useful but still an estimate."
        ),
        "",
        "**Step conclusion:** the candidate meets the numerical full-atlas target, but this does not override a direct external NK regression.",
        "",
        "## Promotion gate",
        "",
        small_markdown_table(["Gate", "Status", "Evidence"], gate_rows),
        "",
        "**Final conclusion:** keep this model as a packaged v3 candidate, not as promoted `gdTAI_v3.0`. The reason is not low overall performance; the reason is that the named failure mode, NK contamination, is worse than v2 high-purity on the independent external test.",
        "",
        "## Detailed artifact paths",
        "",
        f"- Model package: `{MODEL_DIR.relative_to(ROOT)}`",
        f"- Full result tables: `{TABLE_DIR.relative_to(ROOT)}`",
        f"- Figures: `{ASSET_DIR.relative_to(ROOT)}`",
        f"- HTML report: `{REPORT_HTML.relative_to(ROOT)}`",
        f"- PDF report: `{REPORT_PDF.relative_to(ROOT)}`",
    ]
    REPORT_MD.write_text("\n".join(str(x) for x in md_sections) + "\n", encoding="utf-8")

    cards = [
        ("Decision", "Do not promote", "The v3 candidate remains packaged but should not replace v2 high-purity yet.", "bad"),
        ("External F1", fmt_dec(selected_metrics["f1"], 3), f"v2 high-purity: {fmt_dec(v2hp_metrics['f1'], 3)}", "neutral"),
        ("External recall", fmt_pct(selected_metrics["recall"], 1), "above 0.8 target, slightly below v2 high-purity", "neutral"),
        ("Full-atlas recall", fmt_pct(selected_full["full_primary_recall"], 2), "primary-gold labels only", "good"),
        ("Estimated atlas FP/pred", fmt_pct(selected_full["estimated_fp_fraction_of_predictions"], 2), "below 5% target", "good"),
        ("External NK FP", f"{fmt_int(acceptance['external_NK_FP'])} vs {fmt_int(acceptance['v2_high_purity_NK_FP'])}", "failed promotion gate", "bad"),
    ]
    card_html = "".join(
        f"<div class='card {kind}'><div class='card-title'>{html.escape(title)}</div><div class='card-value'>{html.escape(value)}</div><div class='card-note'>{html.escape(note)}</div></div>"
        for title, value, note, kind in cards
    )
    figure_html = "".join(
        f"<figure><img src='{html.escape(rel_asset(path))}'><figcaption>{html.escape(path.stem.replace('decision_', '').replace('_', ' '))}</figcaption></figure>"
        for path in figures
    )

    css = """
    @page{size:A4;margin:12mm}
    body{font-family:Arial,Helvetica,sans-serif;margin:0;color:#18212b;background:#f5f7fa;line-height:1.48}
    main{max-width:1180px;margin:0 auto;padding:22px}
    section{background:white;border:1px solid #dce3eb;border-radius:8px;margin:14px 0;padding:18px;break-inside:avoid}
    h1{font-size:28px;margin:0 0 8px}h2{font-size:20px;margin:0 0 10px}h3{font-size:15px;margin:16px 0 8px}
    p{margin:8px 0}.muted{color:#5f6b7a}.bad-text{color:#b91c1c;font-weight:700}.good-text{color:#047857;font-weight:700}
    .cards{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:10px;margin-top:12px}
    .card{border:1px solid #dce3eb;border-left-width:5px;border-radius:7px;padding:12px;background:#fbfdff}
    .card.good{border-left-color:#059669}.card.bad{border-left-color:#dc2626}.card.neutral{border-left-color:#2563eb}
    .card-title{font-size:12px;text-transform:uppercase;color:#5f6b7a;letter-spacing:.02em}.card-value{font-size:24px;font-weight:700;margin:4px 0}.card-note{font-size:13px;color:#4b5563}
    .conclusion{border-left:5px solid #2563eb;background:#eff6ff;padding:10px 12px;border-radius:6px;margin-top:12px}
    .warning{border-left-color:#dc2626;background:#fef2f2}
    .flow{display:grid;grid-template-columns:repeat(5,1fr);gap:8px;margin:12px 0}
    .flow div{border:1px solid #dce3eb;border-radius:7px;padding:10px;background:#f8fafc;text-align:center;font-size:13px}
    table{border-collapse:collapse;width:100%;font-size:13px;margin-top:8px}th,td{border:1px solid #dce3eb;padding:7px;text-align:left;vertical-align:top}th{background:#eef2f7}
    .figures{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:12px}
    figure{margin:0;border:1px solid #dce3eb;border-radius:7px;padding:8px;background:#fbfdff;break-inside:avoid}
    img{width:100%;height:auto}figcaption{font-size:12px;color:#5f6b7a;margin-top:4px}
    code{background:#eef2f7;padding:1px 4px;border-radius:4px}
    @media(max-width:850px){main{padding:10px}.cards,.figures,.flow{grid-template-columns:1fr}}
    """
    html_doc = f"""<!doctype html>
<html>
<head><meta charset='utf-8'><title>gdTAI v3 TRDC/NK Guard Decision Report</title><style>{css}</style></head>
<body><main>
<section>
  <h1>gdTAI v3 TRDC/NK Guard Decision Report</h1>
  <p class='bad-text'>Decision: do not promote this candidate as the default gdTAI v3.0 model yet.</p>
  <p>The selected candidate meets the requested atlas-level recall and estimated false-positive target, but it fails the explicit NK false-positive promotion gate on the independent external test.</p>
  <div class='cards'>{card_html}</div>
</section>

<section>
  <h2>1. What changed in this iteration</h2>
  <div class='flow'>
    <div>GSE144469 moved into train/tune</div>
    <div>Other validation cohorts retained</div>
    <div>NK hard negatives restricted to TRDC+/TRDV-</div>
    <div>Independent external H5AD held out</div>
    <div>Full 5.13M-cell atlas applied</div>
  </div>
  <p><code>GSE144469</code> contributed {fmt_int(gse144469_row['n_cells'])} train/tune cells ({fmt_int(gse144469_row['gdT_gold'])} gdT positives and {fmt_int(gse144469_row['abT_gold'])} abT negatives). The non-GSE144469 validation block remained held out with {fmt_int(validation_row['n_cells'])} cells.</p>
  <p>Hard NK negatives follow the requested rule: {fmt_int(combined_nk)} NK <code>TRDC+TRDV-</code> cells from {nk_sources} datasets were retained ({fmt_int(train_nk)} train, {fmt_int(tune_nk)} tune). NK cells outside this expression rule were excluded from NK hard negatives ({fmt_int(excluded_nk)} cells). NK+TCRAB-overlap cells were not used as explicit NK hard negatives; {fmt_int(weak_overlap_excluded)} expression-passing weak-TRDC overlap cells were also removed.</p>
  <div class='conclusion'>Step conclusion: the split design now matches the requested policy, and the external H5AD remains an independent final test.</div>
</section>

<section>
  <h2>2. What the algorithm is</h2>
  <p>The selected candidate is <code>{html.escape(selected)}</code> at threshold <code>0.936</code>. It is <b>not XGBoost</b>. It is a conditional guard around the gdTAI v2 transcriptome score, with engineered features for TCR gene evidence, CD3 strength, NK-marker strength, and TRDC-only risk.</p>
  <p>External inference requires raw counts from <code>layers["counts"]</code> and computes <code>log1p(counts per 10,000)</code>. It does not require TCR-seq metadata at inference time.</p>
  <div class='conclusion'>Step conclusion: this model is a conservative v2-score-plus-guard candidate, designed to reduce specific false-positive modes without making NK markers automatic exclusion genes.</div>
</section>

<section>
  <h2>3. How the candidate was selected</h2>
  <p>The initial internal winner was <code>{html.escape(previous_internal)}</code>, but it reached only {fmt_pct(previous_target['full_primary_recall'], 2)} full-atlas primary-gold recall. The final selected round reached {fmt_pct(selected_target['full_primary_recall'], 2)} recall and {fmt_pct(selected_target['estimated_fp_fraction_of_predictions'], 2)} estimated FP/prediction, passing the {fmt_pct(recall_target, 1)} recall-with-margin and {fmt_pct(fp_target, 1)} estimated-FP targets.</p>
  <p>External test results were not used for model selection. The best external-F1 v3 row was <code>{html.escape(str(best_external_v3['strategy']))}</code> with F1 {fmt_dec(best_external_v3['f1'], 4)}, but that was evaluated only after the selection decision.</p>
  <div class='figures'><figure><img src='{html.escape(rel_asset(figures[1]))}'><figcaption>Full-atlas recall versus estimated FP fraction. The selected point passes the target box.</figcaption></figure><figure><img src='{html.escape(rel_asset(figures[0]))}'><figcaption>External primary metrics after selection.</figcaption></figure></div>
  <div class='conclusion'>Step conclusion: round14 was selected because it is the first practical operating point satisfying the atlas-level target, not because it looked best on the external test.</div>
</section>

<section>
  <h2>4. External test result</h2>
  {html_table(["Metric", "v3 selected", "v2 high-purity", "Interpretation"], metric_rows)}
  <p>v3 reduced total external false positives from {fmt_int(v2hp_metrics['fp'])} to {fmt_int(selected_metrics['fp'])}. It also reduced paired-TCRAB false positives from {fmt_int(acceptance['v2_high_purity_paired_TCRAB_FP'])} to {fmt_int(acceptance['external_paired_TCRAB_FP'])} and TRDC+/TRDV- false positives from {fmt_int(acceptance['v2_high_purity_TRDC_plus_TRDV_minus_FP'])} to {fmt_int(acceptance['external_TRDC_plus_TRDV_minus_FP'])}.</p>
  <p class='bad-text'>The problem is NK specificity: NK false positives increased from {fmt_int(acceptance['v2_high_purity_NK_FP'])} in v2 high-purity to {fmt_int(acceptance['external_NK_FP'])} in selected v3.</p>
  <div class='figures'><figure><img src='{html.escape(rel_asset(figures[3]))}'><figcaption>False-positive stress groups.</figcaption></figure><figure><img src='{html.escape(rel_asset(figures[4]))}'><figcaption>Recall stress groups.</figcaption></figure></div>
  <div class='conclusion warning'>Step conclusion: v3 is competitive overall, but it does not convincingly fix the NK false-positive failure mode.</div>
</section>

<section>
  <h2>5. Full-atlas application</h2>
  {html_table(["Full-atlas metric", "v3 selected", "v2 high-purity"], full_rows)}
  <p>The full-atlas recall is measured on primary gold labels only: {fmt_int(selected_full['full_primary_tp'])} / ({fmt_int(selected_full['full_primary_tp'])} + {fmt_int(selected_full['full_primary_fn'])}) = {fmt_pct(selected_full['full_primary_recall'], 2)}. This is not a claim of absolute biological recall for every unlabeled gdT cell.</p>
  <p>The estimated abT FP count uses observed paired-TCRAB false-positive behavior in TCR-seq sources and extrapolates it to datasets without TCR-seq. It is the best current estimate, but it is still an estimate.</p>
  <div class='figures'><figure><img src='{html.escape(rel_asset(figures[5]))}'><figcaption>External confusion matrix for selected v3.</figcaption></figure><figure><img src='{html.escape(rel_asset(figures[2]))}'><figcaption>Promotion gate checklist.</figcaption></figure></div>
  <div class='conclusion'>Step conclusion: the candidate meets the full-atlas numeric target, but the direct external NK regression prevents promotion.</div>
</section>

<section>
  <h2>6. Why this should not be promoted yet</h2>
  {html_table(["Gate", "Status", "Evidence"], gate_rows)}
  <p>The promotion rule was intentionally stricter than overall F1. The previous concern was NK-like contamination, especially around weak TRD evidence. A release model should not improve some aggregate metrics while worsening the named failure mode.</p>
  <p>Therefore the correct status is: <b>package and keep as a v3 candidate</b>, but do not write it as the promoted <code>gdTAI_v3.0</code> default model.</p>
  <div class='conclusion warning'>Final conclusion: not promoted because external NK FP is worse than v2 high-purity, while external F1 is essentially unchanged and external recall is slightly lower than v2 high-purity.</div>
</section>

<section>
  <h2>7. Detailed artifacts</h2>
  <p>Large tables are intentionally not embedded in this report. They remain available for audit:</p>
  <ul>
    <li>Model package: <code>{html.escape(str(MODEL_DIR.relative_to(ROOT)))}</code></li>
    <li>Full result tables: <code>{html.escape(str(TABLE_DIR.relative_to(ROOT)))}</code></li>
    <li>Figures: <code>{html.escape(str(ASSET_DIR.relative_to(ROOT)))}</code></li>
    <li>Markdown report: <code>{html.escape(str(REPORT_MD.relative_to(ROOT)))}</code></li>
    <li>PDF report: <code>{html.escape(str(REPORT_PDF.relative_to(ROOT)))}</code></li>
  </ul>
</section>
</main></body></html>
"""
    REPORT_HTML.write_text(html_doc, encoding="utf-8")

    if not no_pdf:
        chrome = shutil.which("google-chrome")
        if chrome is None:
            raise RuntimeError("google-chrome not found; rerun with --no-pdf or install Chrome.")
        subprocess.run(
            [
                chrome,
                "--headless",
                "--disable-gpu",
                "--no-sandbox",
                "--print-to-pdf-no-header",
                f"--print-to-pdf={REPORT_PDF.resolve()}",
                REPORT_HTML.resolve().as_uri(),
            ],
            check=True,
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-pdf", action="store_true", help="Write markdown and HTML only.")
    args = parser.parse_args()
    build_report(no_pdf=args.no_pdf)
    print(f"Wrote {REPORT_HTML.relative_to(ROOT)}")
    if not args.no_pdf:
        print(f"Wrote {REPORT_PDF.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
