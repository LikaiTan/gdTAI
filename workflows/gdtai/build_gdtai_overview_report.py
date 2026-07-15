#!/usr/bin/env python3
"""Build a self-contained gdTAI overview and performance report.

The report is intentionally a presentation layer over existing gdTAI outputs.
It reads canonical CSVs and PNG figures, creates a few overview graphics, and
writes a standalone HTML file that can be sent as one attachment.
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


import base64
import hashlib
import html
import math
import textwrap
from datetime import datetime
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch


ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = ROOT / "Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene"
FAILURE_TABLE_DIR = ROOT / "Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene_failure_modes"
FIGURE_DIR = ROOT / "Integrated_dataset/figures/gdT_prediction/gse144469_holdout_tcrgene"
FAILURE_FIGURE_DIR = ROOT / "Integrated_dataset/figures/gdT_prediction/gse144469_holdout_tcrgene_failure_modes"
OVERVIEW_FIGURE_DIR = ROOT / "Integrated_dataset/figures/gdT_prediction/gdTAI_overview"
LOG_DIR = ROOT / "Integrated_dataset/logs/gdT_prediction"
STATIC_DIR = ROOT / "gdT_prediction"
MODEL_PATH = ROOT / "Integrated_dataset/models/gdT_prediction_classifier/gse144469_holdout_tcrgene/selected_model.pkl"
USAGE_PROTOCOL = ROOT / "Integrated_dataset/models/gdT_prediction_classifier/gse144469_holdout_tcrgene/USAGE_PROTOCOL.md"
ATLAS_UMAP = ROOT / "Integrated_dataset/figures/plus6/plus6_umap_by_corrected_simple_annotation_100x120mm.png"

SELECTED_MODEL = "xgb_individual_TCRABGD_genes_with_FOXP3_CD4_lowCD3_penalty"
ORIGINAL_STRATEGY = "original_TRD_minus_TRAB_score_strategy"


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def fmt_int(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return f"{int(round(float(value))):,}"


def fmt_num(value: object, digits: int = 3) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return f"{float(value):.{digits}f}"


def fmt_pct(value: object, digits: int = 1) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return f"{float(value) * 100:.{digits}f}%"


def short_strategy(value: str) -> str:
    if value == SELECTED_MODEL or value.startswith("selected_model_"):
        return "gdTAI"
    if value == ORIGINAL_STRATEGY:
        return "Original TRD-TRAB"
    if value == "baseline_individual_TRD_gene_sum":
        return "TRD-gene baseline"
    if value == "baseline_individual_gd_minus_ab_gene_sum":
        return "gd-minus-ab baseline"
    if value == "xgb_individual_TCRABGD_genes":
        return "XGBoost TCR genes"
    if value == SELECTED_MODEL:
        return "gdTAI"
    if value == "logistic_individual_TCRABGD_genes":
        return "Logistic TCR genes"
    return value.replace("_", " ")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def image_data_uri(path: Path) -> str:
    data = path.read_bytes()
    encoded = base64.b64encode(data).decode("ascii")
    suffix = path.suffix.lower().lstrip(".") or "png"
    if suffix == "jpg":
        suffix = "jpeg"
    return f"data:image/{suffix};base64,{encoded}"


def image_block(path: Path, caption: str, css_class: str = "") -> str:
    if not path.exists():
        return f'<div class="missing">Missing figure: <code>{html.escape(rel(path))}</code></div>'
    cls = f' class="figure {css_class}"' if css_class else ' class="figure"'
    return (
        f"<figure{cls}>"
        f'<img src="{image_data_uri(path)}" alt="{html.escape(caption)}">'
        f"<figcaption>{html.escape(caption)}</figcaption>"
        f"</figure>"
    )


def format_table(
    df: pd.DataFrame,
    columns: Iterable[str] | None = None,
    rename: dict[str, str] | None = None,
    percent_columns: Iterable[str] = (),
    integer_columns: Iterable[str] = (),
    float_columns: Iterable[str] = (),
    max_rows: int | None = None,
) -> str:
    out = df.copy()
    if columns is not None:
        out = out[list(columns)]
    if max_rows is not None:
        out = out.head(max_rows)
    for col in percent_columns:
        if col in out.columns:
            out[col] = out[col].map(lambda x: fmt_pct(x, 2))
    for col in integer_columns:
        if col in out.columns:
            out[col] = out[col].map(fmt_int)
    for col in float_columns:
        if col in out.columns:
            out[col] = out[col].map(lambda x: fmt_num(x, 3))
    if rename:
        out = out.rename(columns=rename)
    return out.to_html(index=False, escape=True, classes="data-table", border=0)


def markdown_table(
    df: pd.DataFrame,
    columns: Iterable[str] | None = None,
    rename: dict[str, str] | None = None,
    percent_columns: Iterable[str] = (),
    integer_columns: Iterable[str] = (),
    float_columns: Iterable[str] = (),
    max_rows: int | None = None,
) -> str:
    out = df.copy()
    if columns is not None:
        out = out[list(columns)]
    if max_rows is not None:
        out = out.head(max_rows)
    for col in percent_columns:
        if col in out.columns:
            out[col] = out[col].map(lambda x: fmt_pct(x, 2))
    for col in integer_columns:
        if col in out.columns:
            out[col] = out[col].map(fmt_int)
    for col in float_columns:
        if col in out.columns:
            out[col] = out[col].map(lambda x: fmt_num(x, 3))
    if rename:
        out = out.rename(columns=rename)
    return dataframe_to_markdown(out)


def dataframe_to_markdown(df: pd.DataFrame) -> str:
    """Render a compact pipe table without pandas optional tabulate dependency."""
    if df.empty:
        return "_No rows._"
    str_df = df.fillna("").astype(str)
    columns = list(str_df.columns)
    widths = {
        col: max(len(str(col)), *(len(value) for value in str_df[col].tolist()))
        for col in columns
    }
    header = "| " + " | ".join(str(col).ljust(widths[col]) for col in columns) + " |"
    sep = "| " + " | ".join("-" * widths[col] for col in columns) + " |"
    rows = [
        "| " + " | ".join(str(row[col]).ljust(widths[col]) for col in columns) + " |"
        for _, row in str_df.iterrows()
    ]
    return "\n".join([header, sep, *rows])


def save_headline_figure(selected: pd.Series, full_selected: pd.Series, fp_selected: pd.Series) -> Path:
    out = OVERVIEW_FIGURE_DIR / "gdtai_headline_metrics.png"
    fig, ax = plt.subplots(figsize=(11.5, 4.4), dpi=180)
    ax.axis("off")
    cards = [
        ("Precision", fmt_pct(selected["precision"], 1), f"TP {fmt_int(selected['tp'])}, FP {fmt_int(selected['fp'])}"),
        ("Specificity", fmt_pct(selected["specificity"], 2), f"TN {fmt_int(selected['tn'])} / negatives {fmt_int(selected['n_negative'])}"),
        ("Recall", fmt_pct(selected["recall"], 1), f"TP {fmt_int(selected['tp'])} / positives {fmt_int(selected['n_positive'])}"),
        ("F1", fmt_num(selected["f1"], 3), "selected precision-recall balance"),
        ("Recovered in atlas", fmt_int(full_selected["predicted_putative_gdT"]), "5-million-cell TNK atlas"),
        ("Known FP fraction", fmt_pct(fp_selected["known_FP_fraction_of_predictions"], 2), "paired TCRAB or NK audit"),
    ]
    colors = ["#0f766e", "#1d4ed8", "#b45309", "#7c3aed", "#047857", "#be123c"]
    for idx, (title, value, note) in enumerate(cards):
        row, col = divmod(idx, 3)
        x = 0.04 + col * 0.32
        y = 0.56 - row * 0.44
        rect = FancyBboxPatch(
            (x, y),
            0.28,
            0.32,
            boxstyle="round,pad=0.018,rounding_size=0.018",
            facecolor="#ffffff",
            edgecolor="#d7dee8",
            linewidth=1.2,
            transform=ax.transAxes,
        )
        ax.add_patch(rect)
        ax.text(x + 0.025, y + 0.225, title, transform=ax.transAxes, fontsize=10.5, color="#475569")
        ax.text(x + 0.025, y + 0.118, value, transform=ax.transAxes, fontsize=23, fontweight="bold", color=colors[idx])
        ax.text(x + 0.025, y + 0.047, note, transform=ax.transAxes, fontsize=9.2, color="#64748b")
    fig.suptitle("gdTAI performance at the selected F1-optimized threshold", fontsize=15, fontweight="bold", y=0.98)
    fig.savefig(out, bbox_inches="tight", facecolor="#f8fafc")
    plt.close(fig)
    return out


def save_strategy_comparison(full_overall: pd.DataFrame) -> Path:
    out = OVERVIEW_FIGURE_DIR / "gdtai_vs_original_score_strategy.png"
    df = full_overall.copy()
    df["label"] = df["strategy"].map(short_strategy)
    df = df.set_index("label")
    order = ["gdTAI", "Original TRD-TRAB"]
    df = df.loc[[label for label in order if label in df.index]]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7), dpi=180, gridspec_kw={"width_ratios": [1.5, 1]})
    count_metrics = [
        ("predicted_putative_gdT", "Recovered gdT calls"),
        ("predicted_paired_TCRAB_FP", "Paired TCRAB FP"),
        ("predicted_NK_FP", "NK FP"),
        ("predicted_paired_TCRAB_or_NK_FP", "Known FP union"),
    ]
    x = np.arange(len(count_metrics))
    width = 0.34
    colors = {"gdTAI": "#0f766e", "Original TRD-TRAB": "#94a3b8"}
    for offset, label in [(-width / 2, "gdTAI"), (width / 2, "Original TRD-TRAB")]:
        if label in df.index:
            axes[0].bar(
                x + offset,
                [df.loc[label, metric] for metric, _ in count_metrics],
                width,
                label=label,
                color=colors[label],
            )
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([label for _, label in count_metrics], rotation=25, ha="right")
    axes[0].set_ylabel("Cells")
    axes[0].set_yscale("log")
    axes[0].grid(axis="y", linestyle=":", alpha=0.35)
    axes[0].legend(frameon=False)
    axes[0].set_title("Full-atlas calls and obvious FP counts")

    metric = "known_FP_fraction_of_predictions"
    vals = [df.loc[label, metric] * 100 for label in df.index]
    axes[1].bar(df.index.tolist(), vals, color=[colors.get(label, "#64748b") for label in df.index])
    axes[1].set_ylabel("Percent of predicted gdT calls")
    axes[1].set_title("Known FP fraction")
    axes[1].grid(axis="y", linestyle=":", alpha=0.35)
    for i, val in enumerate(vals):
        axes[1].text(i, val + 0.12, f"{val:.2f}%", ha="center", va="bottom", fontsize=10)
    fig.suptitle("gdTAI reduces paired-TCRAB and NK false positives versus the original score rule", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


def save_feature_family_figure(feature_manifest: pd.DataFrame) -> Path:
    out = OVERVIEW_FIGURE_DIR / "gdtai_feature_family_counts.png"
    counts = (
        feature_manifest["family"]
        .value_counts()
        .rename_axis("family")
        .reset_index(name="n_features")
        .sort_values("n_features", ascending=True)
    )
    colors = ["#64748b", "#64748b", "#0f766e", "#0f766e", "#1d4ed8", "#1d4ed8"][-len(counts) :]
    fig, ax = plt.subplots(figsize=(7.4, 4.2), dpi=180)
    ax.barh(counts["family"], counts["n_features"], color=colors)
    ax.set_xlabel("Number of features")
    ax.set_ylabel("")
    ax.set_title("gdTAI feature families", fontweight="bold")
    ax.grid(axis="x", linestyle=":", alpha=0.35)
    for y, val in enumerate(counts["n_features"]):
        ax.text(val + 1.0, y, str(int(val)), va="center", fontsize=10)
    ax.set_xlim(0, max(counts["n_features"]) * 1.18)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


def save_workflow_figure() -> Path:
    out = OVERVIEW_FIGURE_DIR / "gdtai_training_workflow.png"
    fig, ax = plt.subplots(figsize=(12, 3.2), dpi=180)
    ax.axis("off")
    boxes = [
        ("Gold labels", "sorted gdT + TCR evidence\nsilver cells kept sensitivity-only"),
        ("Leakage control", "hold out all GSE144469,\nGSE254249 abT, cord-blood gdT"),
        ("Features", "187 individual TCR genes\nplus FOXP3/CD4/CD3 controls"),
        ("Tune threshold", "maximize F1\non tune split"),
        ("Validate/apply", "combined holdout metrics\nthen full 5M-cell atlas"),
    ]
    x_positions = np.linspace(0.03, 0.79, len(boxes))
    for idx, ((title, body), x) in enumerate(zip(boxes, x_positions)):
        rect = FancyBboxPatch(
            (x, 0.28),
            0.18,
            0.48,
            boxstyle="round,pad=0.014,rounding_size=0.02",
            facecolor="#ffffff",
            edgecolor="#cbd5e1",
            linewidth=1.1,
            transform=ax.transAxes,
        )
        ax.add_patch(rect)
        ax.text(x + 0.015, 0.62, title, transform=ax.transAxes, fontsize=10.5, fontweight="bold", color="#0f172a")
        ax.text(x + 0.015, 0.43, body, transform=ax.transAxes, fontsize=8.8, color="#475569", linespacing=1.25)
        if idx < len(boxes) - 1:
            ax.annotate(
                "",
                xy=(x + 0.197, 0.52),
                xytext=(x + 0.184, 0.52),
                xycoords=ax.transAxes,
                arrowprops={"arrowstyle": "->", "lw": 1.6, "color": "#64748b"},
            )
    fig.suptitle("How gdTAI was built and evaluated", fontsize=14, fontweight="bold", y=0.95)
    fig.savefig(out, bbox_inches="tight", facecolor="#f8fafc")
    plt.close(fig)
    return out


def make_css() -> str:
    return """
    :root {
      --ink: #102033;
      --muted: #536477;
      --line: #d8e0ea;
      --panel: #ffffff;
      --bg: #f5f7fb;
      --teal: #0f766e;
      --blue: #1d4ed8;
      --amber: #b45309;
      --rose: #be123c;
    }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Arial, sans-serif;
      color: var(--ink);
      background: var(--bg);
      line-height: 1.55;
    }
    .hero {
      background: linear-gradient(135deg, #102033 0%, #123d48 52%, #0f766e 100%);
      color: #fff;
      padding: 46px 34px 34px;
    }
    .hero-inner { max-width: 1160px; margin: 0 auto; }
    .eyebrow { text-transform: uppercase; letter-spacing: .12em; font-size: 12px; opacity: .82; }
    h1 { margin: 10px 0 12px; font-size: 42px; line-height: 1.05; letter-spacing: 0; }
    .subtitle { max-width: 920px; font-size: 18px; opacity: .92; }
    .meta { display: flex; flex-wrap: wrap; gap: 10px; margin-top: 22px; }
    .pill {
      border: 1px solid rgba(255,255,255,.26);
      background: rgba(255,255,255,.10);
      padding: 7px 11px;
      border-radius: 999px;
      font-size: 13px;
    }
    main { max-width: 1160px; margin: 0 auto; padding: 26px 24px 64px; }
    section {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 10px;
      padding: 24px;
      margin: 18px 0;
      box-shadow: 0 8px 24px rgba(15, 23, 42, .04);
    }
    h2 { margin: 0 0 14px; font-size: 24px; letter-spacing: 0; }
    h3 { margin: 22px 0 10px; font-size: 17px; letter-spacing: 0; }
    p { margin: 10px 0; }
    .lead { font-size: 17px; color: #26374a; }
    .grid { display: grid; gap: 16px; }
    .grid.two { grid-template-columns: repeat(2, minmax(0, 1fr)); }
    .grid.three { grid-template-columns: repeat(3, minmax(0, 1fr)); }
    .metric {
      border: 1px solid var(--line);
      border-radius: 9px;
      padding: 15px;
      background: #fbfdff;
    }
    .metric .label { color: var(--muted); font-size: 13px; }
    .metric .value { font-size: 27px; font-weight: 750; margin-top: 4px; }
    .metric .note { color: var(--muted); font-size: 12.5px; margin-top: 4px; }
    .value.teal { color: var(--teal); }
    .value.blue { color: var(--blue); }
    .value.amber { color: var(--amber); }
    .value.rose { color: var(--rose); }
    ul { padding-left: 20px; margin: 10px 0; }
    li { margin: 5px 0; }
    code {
      background: #eef3f8;
      border: 1px solid #dce5ee;
      border-radius: 4px;
      padding: 1px 5px;
      font-size: .92em;
    }
    pre {
      background: #0f172a;
      color: #e2e8f0;
      padding: 14px;
      border-radius: 8px;
      overflow-x: auto;
      font-size: 13px;
    }
    .figure { margin: 12px 0 6px; }
    .figure img {
      width: 100%;
      height: auto;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: #fff;
      display: block;
    }
    .figure.small img { max-height: 520px; object-fit: contain; }
    figcaption { color: var(--muted); font-size: 13px; margin-top: 7px; }
    .table-wrap { overflow-x: auto; margin: 12px 0; border: 1px solid var(--line); border-radius: 8px; }
    table.data-table { border-collapse: collapse; width: 100%; min-width: 720px; font-size: 13px; }
    table.data-table th {
      background: #edf3f8;
      color: #1e293b;
      text-align: left;
      padding: 9px 10px;
      border-bottom: 1px solid var(--line);
      white-space: nowrap;
    }
    table.data-table td {
      padding: 8px 10px;
      border-bottom: 1px solid #e8eef5;
      vertical-align: top;
      white-space: nowrap;
    }
    .callout {
      border-left: 4px solid var(--teal);
      background: #ecfdf5;
      padding: 12px 14px;
      border-radius: 6px;
      color: #12372f;
    }
    .warning {
      border-left: 4px solid var(--amber);
      background: #fff7ed;
      padding: 12px 14px;
      border-radius: 6px;
      color: #3b2a16;
    }
    .missing {
      border: 1px dashed #cbd5e1;
      padding: 14px;
      border-radius: 8px;
      color: var(--muted);
      background: #f8fafc;
    }
    .captioned-list strong { color: #0f172a; }
    footer { color: var(--muted); font-size: 12px; margin-top: 22px; }
    @media (max-width: 820px) {
      h1 { font-size: 34px; }
      main { padding: 18px 14px 42px; }
      section { padding: 18px; }
      .grid.two, .grid.three { grid-template-columns: 1fr; }
    }
    """


def build_report() -> tuple[Path, Path]:
    global SELECTED_MODEL
    OVERVIEW_FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    STATIC_DIR.mkdir(parents=True, exist_ok=True)

    validation = read_csv(TABLE_DIR / "validation_metrics.csv")
    split = read_csv(TABLE_DIR / "split_overall.csv")
    acceptance = read_csv(TABLE_DIR / "acceptance_vs_individual_gene_baseline.csv")
    feature_manifest = read_csv(TABLE_DIR / "feature_manifest.csv")
    by_source = read_csv(TABLE_DIR / "selected_validation_metrics_by_source.csv")
    full_overall = read_csv(TABLE_DIR / "selected_model_full_dataset_prediction_overall.csv")
    full_by_source = read_csv(TABLE_DIR / "selected_model_full_dataset_prediction_by_source.csv")
    full_by_tissue = read_csv(TABLE_DIR / "selected_model_full_dataset_prediction_by_tissue.csv")
    full_by_annotation = read_csv(TABLE_DIR / "selected_model_full_dataset_prediction_by_annotation.csv")
    fp_overall = read_csv(TABLE_DIR / "selected_model_full_dataset_false_positive_overall.csv")
    ppv = read_csv(TABLE_DIR / "prevalence_aware_ppv_scenarios.csv")
    failure_summary = read_csv(FAILURE_TABLE_DIR / "failure_mode_group_summary.csv")
    failure_flags = read_csv(FAILURE_TABLE_DIR / "failure_mode_flag_enrichment.csv")

    full_selected = full_overall.loc[full_overall["strategy"].str.startswith("selected_model_")].iloc[0]
    full_original = full_overall.loc[full_overall["strategy"] == ORIGINAL_STRATEGY].iloc[0]
    fp_selected = fp_overall.loc[fp_overall["strategy"].str.startswith("selected_model_")].iloc[0]
    fp_original = fp_overall.loc[fp_overall["strategy"] == ORIGINAL_STRATEGY].iloc[0]
    SELECTED_MODEL = str(full_selected["strategy"]).replace("selected_model_", "", 1)
    selected = validation.loc[validation["model"] == SELECTED_MODEL].iloc[0]
    selected_acceptance = acceptance.loc[acceptance["model"] == SELECTED_MODEL].iloc[0]

    headline_fig = save_headline_figure(selected, full_selected, fp_selected)
    strategy_fig = save_strategy_comparison(full_overall)
    feature_fig = save_feature_family_figure(feature_manifest)
    workflow_fig = save_workflow_figure()

    feature_counts = (
        feature_manifest.groupby("family", as_index=False)
        .size()
        .rename(columns={"size": "n_features"})
        .sort_values(["n_features", "family"], ascending=[False, True])
    )
    top_features = feature_manifest.groupby("family").head(8).copy()

    strategy_compare = full_overall.copy()
    strategy_compare["strategy_short"] = strategy_compare["strategy"].map(short_strategy)
    source_top = (
        full_by_source.loc[full_by_source["strategy"].str.startswith("selected_model_")]
        .sort_values("predicted_putative_gdT", ascending=False)
        .head(12)
        .copy()
    )
    tissue_top = (
        full_by_tissue.loc[full_by_tissue["strategy"].str.startswith("selected_model_")]
        .sort_values("predicted_putative_gdT", ascending=False)
        .head(10)
        .copy()
    )
    annotation_compare = full_by_annotation.copy()
    annotation_compare["strategy_short"] = annotation_compare["strategy"].map(short_strategy)
    annotation_compare = annotation_compare.sort_values(["annotation", "strategy_short"])

    ppv_selected = ppv.loc[ppv["model"] == SELECTED_MODEL].copy()
    if ppv_selected.empty:
        ppv_selected = ppv.loc[ppv["model"].str.contains("FOXP3", regex=False)].copy()

    fp_tcrab_reduction = 1 - (fp_selected["predicted_paired_TCRAB_FP"] / fp_original["predicted_paired_TCRAB_FP"])
    fp_nk_reduction = 1 - (fp_selected["predicted_NK_FP"] / fp_original["predicted_NK_FP"])
    known_fp_reduction = 1 - (
        fp_selected["known_FP_fraction_of_predictions"] / fp_original["known_FP_fraction_of_predictions"]
    )
    model_sha = sha256_file(MODEL_PATH) if MODEL_PATH.exists() else "missing"
    model_size = MODEL_PATH.stat().st_size if MODEL_PATH.exists() else 0
    report_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    metrics_cards = "\n".join(
        [
            f'<div class="metric"><div class="label">Holdout precision</div><div class="value teal">{fmt_pct(selected["precision"], 1)}</div><div class="note">TP {fmt_int(selected["tp"])}, FP {fmt_int(selected["fp"])}</div></div>',
            f'<div class="metric"><div class="label">Holdout specificity</div><div class="value blue">{fmt_pct(selected["specificity"], 2)}</div><div class="note">TN {fmt_int(selected["tn"])} / negatives {fmt_int(selected["n_negative"])}</div></div>',
            f'<div class="metric"><div class="label">Holdout recall</div><div class="value amber">{fmt_pct(selected["recall"], 1)}</div><div class="note">FN {fmt_int(selected["fn"])}</div></div>',
            f'<div class="metric"><div class="label">PR-AUC</div><div class="value teal">{fmt_num(selected["pr_auc"], 3)}</div><div class="note">Class-imbalance aware summary</div></div>',
            f'<div class="metric"><div class="label">Full-atlas recovery</div><div class="value blue">{fmt_int(full_selected["predicted_putative_gdT"])}</div><div class="note">{fmt_pct(full_selected["predicted_fraction"], 2)} of 5,128,904 cells</div></div>',
            f'<div class="metric"><div class="label">Known FP fraction</div><div class="value rose">{fmt_pct(fp_selected["known_FP_fraction_of_predictions"], 2)}</div><div class="note">paired TCRAB or NK audit</div></div>',
        ]
    )

    validation_table = format_table(
        validation.assign(model_short=validation["model"].map(short_strategy)),
        columns=[
            "model_short",
            "threshold",
            "n_cells",
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
        rename={
            "model_short": "model",
            "n_cells": "cells",
            "predicted_positive": "predicted gdT",
            "f1": "F1",
            "f0.5": "F0.5",
            "roc_auc": "ROC-AUC",
            "pr_auc": "PR-AUC",
        },
        integer_columns=["n_cells", "predicted_positive", "tp", "fp", "tn", "fn"],
        float_columns=["threshold", "precision", "recall", "specificity", "f1", "f0.5", "roc_auc", "pr_auc"],
    )
    split_table = format_table(
        split,
        integer_columns=["n_cells", "gdT_gold", "abT_gold"],
        percent_columns=["gdT_prevalence"],
    )
    acceptance_table = format_table(
        acceptance.assign(model_short=acceptance["model"].map(short_strategy)),
        columns=[
            "model_short",
            "best_individual_gene_baseline",
            "validation_f1",
            "baseline_f1",
            "delta_f1",
            "validation_f0.5",
            "delta_f0.5",
            "validation_precision",
            "validation_recall",
            "validation_specificity",
            "accepted",
        ],
        rename={
            "model_short": "model",
            "best_individual_gene_baseline": "baseline",
            "validation_f1": "validation F1",
            "baseline_f1": "baseline F1",
            "delta_f1": "delta F1",
            "validation_f0.5": "validation F0.5",
            "delta_f0.5": "delta F0.5",
        },
        float_columns=[
            "validation_f1",
            "baseline_f1",
            "delta_f1",
            "validation_f0.5",
            "delta_f0.5",
            "validation_precision",
            "validation_recall",
            "validation_specificity",
        ],
    )
    feature_counts_table = format_table(feature_counts, integer_columns=["n_features"])
    top_features_table = format_table(top_features, max_rows=48, integer_columns=["feature_index"])
    by_source_table = format_table(
        by_source,
        integer_columns=["n_cells", "n_positive", "n_negative", "predicted_positive", "tp", "fp", "tn", "fn"],
        float_columns=["precision", "recall", "specificity", "f1", "f0.5", "roc_auc", "pr_auc"],
    )
    strategy_table = format_table(
        strategy_compare,
        columns=[
            "strategy_short",
            "threshold",
            "total_cells",
            "predicted_putative_gdT",
            "predicted_fraction",
            "predicted_paired_TCRAB_FP",
            "paired_TCRAB_FP_fraction_of_predictions",
            "predicted_NK_FP",
            "NK_FP_fraction_of_predictions",
            "predicted_paired_TCRAB_or_NK_FP",
            "known_FP_fraction_of_predictions",
        ],
        rename={
            "strategy_short": "strategy",
            "total_cells": "cells",
            "predicted_putative_gdT": "predicted gdT",
            "predicted_fraction": "predicted fraction",
            "predicted_paired_TCRAB_FP": "paired TCRAB FP",
            "paired_TCRAB_FP_fraction_of_predictions": "paired TCRAB FP / predictions",
            "predicted_NK_FP": "NK FP",
            "NK_FP_fraction_of_predictions": "NK FP / predictions",
            "predicted_paired_TCRAB_or_NK_FP": "known FP union",
            "known_FP_fraction_of_predictions": "known FP / predictions",
        },
        integer_columns=[
            "total_cells",
            "predicted_putative_gdT",
            "predicted_paired_TCRAB_FP",
            "predicted_NK_FP",
            "predicted_paired_TCRAB_or_NK_FP",
        ],
        percent_columns=[
            "predicted_fraction",
            "paired_TCRAB_FP_fraction_of_predictions",
            "NK_FP_fraction_of_predictions",
            "known_FP_fraction_of_predictions",
        ],
        float_columns=["threshold"],
    )
    source_top_table = format_table(
        source_top,
        columns=[
            "source_gse_id",
            "total_cells",
            "predicted_putative_gdT",
            "predicted_fraction",
            "predicted_paired_TCRAB_FP",
            "predicted_NK_FP",
            "known_FP_fraction_of_predictions",
        ],
        rename={
            "source_gse_id": "source",
            "total_cells": "cells",
            "predicted_putative_gdT": "predicted gdT",
            "predicted_fraction": "predicted fraction",
            "predicted_paired_TCRAB_FP": "paired TCRAB FP",
            "predicted_NK_FP": "NK FP",
            "known_FP_fraction_of_predictions": "known FP / predictions",
        },
        integer_columns=["total_cells", "predicted_putative_gdT", "predicted_paired_TCRAB_FP", "predicted_NK_FP"],
        percent_columns=["predicted_fraction", "known_FP_fraction_of_predictions"],
    )
    tissue_top_table = format_table(
        tissue_top,
        columns=["tissue", "total_cells", "predicted_putative_gdT", "predicted_fraction", "known_FP_fraction_of_predictions"],
        rename={
            "total_cells": "cells",
            "predicted_putative_gdT": "predicted gdT",
            "predicted_fraction": "predicted fraction",
            "known_FP_fraction_of_predictions": "known FP / predictions",
        },
        integer_columns=["total_cells", "predicted_putative_gdT"],
        percent_columns=["predicted_fraction", "known_FP_fraction_of_predictions"],
    )
    annotation_table = format_table(
        annotation_compare,
        columns=[
            "strategy_short",
            "annotation",
            "total_cells",
            "predicted_putative_gdT",
            "predicted_fraction",
            "predicted_paired_TCRAB_FP",
            "predicted_NK_FP",
            "known_FP_fraction_of_predictions",
        ],
        rename={
            "strategy_short": "strategy",
            "total_cells": "cells",
            "predicted_putative_gdT": "predicted gdT",
            "predicted_fraction": "predicted fraction",
            "predicted_paired_TCRAB_FP": "paired TCRAB FP",
            "predicted_NK_FP": "NK FP",
            "known_FP_fraction_of_predictions": "known FP / predictions",
        },
        integer_columns=["total_cells", "predicted_putative_gdT", "predicted_paired_TCRAB_FP", "predicted_NK_FP"],
        percent_columns=["predicted_fraction", "known_FP_fraction_of_predictions"],
    )
    ppv_table = format_table(
        ppv_selected,
        columns=[
            "prevalence_percent",
            "ppv_observed_at_prevalence",
            "ppv_conservative_wilson_at_prevalence",
            "expected_false_positive_per_million_observed",
            "expected_false_positive_per_million_conservative",
        ],
        rename={
            "prevalence_percent": "true gdT prevalence",
            "ppv_observed_at_prevalence": "PPV from observed specificity",
            "ppv_conservative_wilson_at_prevalence": "PPV conservative Wilson",
            "expected_false_positive_per_million_observed": "FP / million, observed",
            "expected_false_positive_per_million_conservative": "FP / million, conservative",
        },
        percent_columns=["ppv_observed_at_prevalence", "ppv_conservative_wilson_at_prevalence"],
        float_columns=[
            "prevalence_percent",
            "expected_false_positive_per_million_observed",
            "expected_false_positive_per_million_conservative",
        ],
    )
    failure_summary_table = format_table(
        failure_summary,
        columns=[
            "outcome",
            "n_cells",
            "score_median",
            "source_low_quality_proxy_fraction",
            "multi_lineage_doublet_proxy_fraction",
            "tcr_dual_gene_proxy_fraction",
            "nk_like_proxy_fraction",
            "any_death_penalty_proxy_fraction",
        ],
        rename={
            "score_median": "median score",
            "source_low_quality_proxy_fraction": "low-quality proxy",
            "multi_lineage_doublet_proxy_fraction": "multi-lineage proxy",
            "tcr_dual_gene_proxy_fraction": "TCR dual-gene proxy",
            "nk_like_proxy_fraction": "NK-like proxy",
            "any_death_penalty_proxy_fraction": "death penalty proxy",
        },
        integer_columns=["n_cells"],
        percent_columns=[
            "source_low_quality_proxy_fraction",
            "multi_lineage_doublet_proxy_fraction",
            "tcr_dual_gene_proxy_fraction",
            "nk_like_proxy_fraction",
            "any_death_penalty_proxy_fraction",
        ],
        float_columns=["score_median"],
    )
    failure_flag_table = format_table(
        failure_flags.loc[failure_flags["contrast"].isin(["FP_vs_TN", "FN_vs_TP"])].copy(),
        columns=["contrast", "flag", "fraction_a", "fraction_b", "fraction_diff_a_minus_b", "q_value"],
        rename={
            "fraction_a": "failure fraction",
            "fraction_b": "comparison fraction",
            "fraction_diff_a_minus_b": "difference",
        },
        percent_columns=["fraction_a", "fraction_b", "fraction_diff_a_minus_b"],
        float_columns=["q_value"],
        max_rows=22,
    )

    html_body = f"""
    <header class="hero">
      <div class="hero-inner">
        <div class="eyebrow">TNK atlas classifier report</div>
        <h1>gdTAI: F1-optimized recovery of gamma-delta T cells</h1>
        <p class="subtitle">
          gdTAI is the current gdT AI classifier for the 5-million-cell TNK atlas we built. This report explains
          what the model uses, how it was trained and tested, and how well it performs, with
          performance, recovery, and false-positive audit as the main focus.
        </p>
        <div class="meta">
          <span class="pill">Generated {html.escape(report_time)}</span>
          <span class="pill">Model threshold {fmt_num(selected["threshold"], 6)}</span>
          <span class="pill">187 features, below 300-gene cap</span>
          <span class="pill">No TRD/TRAB score inputs</span>
        </div>
      </div>
    </header>

    <main>
      <section id="summary">
        <h2>Executive Summary</h2>
        <p class="lead">
          The selected gdTAI model is an F1-optimized TCR-gene expression classifier. On the
          combined held-out validation set it reached <strong>{fmt_pct(selected["precision"], 1)}
          precision</strong>, <strong>{fmt_pct(selected["specificity"], 2)} specificity</strong>,
          and <strong>{fmt_pct(selected["recall"], 1)} recall</strong>. It substantially reduced
          obvious false positives compared with the original <code>TRD - TRAB</code> score rule.
        </p>
        <div class="grid three">{metrics_cards}</div>
        <div class="callout">
          Main conclusion: gdTAI is practical as a gdT recovery tool for the TNK atlas when the
          operating point should balance precision and recall. Its false-positive behavior remains
          explicitly audited against paired TCRAB and NK compartments, and it still needs independent
          external testing before being treated as a general-purpose biological truth classifier.
        </div>
        {image_block(headline_fig, "Headline validation and full-atlas performance metrics.")}
      </section>

      <section id="what">
        <h2>What Is gdTAI?</h2>
        <p>
          gdTAI is a classifier that calls putative gamma-delta T cells from expression of
          individual TCR genes and a few penalty/control genes. It was built to improve on
          simple score thresholds while balancing precision and recall in a very imbalanced
          5.1-million-cell atlas.
        </p>
        <ul>
          <li><strong>Used by the model:</strong> individual TCRA, TCRB, TCRG, and TCRD gene expression; FOXP3; CD4; CD3D/CD3E/CD3G.</li>
          <li><strong>Not used by the model:</strong> Phase 4 <code>TRD_score</code>, <code>TRAB_score</code>, <code>TRD_minus_TRAB</code>, CDR3 metadata, productive TCR metadata, dataset ID, tissue, sample ID, or cluster labels.</li>
          <li><strong>Feature scale:</strong> per-cell <code>log1p(counts per 10,000)</code> rebuilt from the H5AD expression matrix.</li>
          <li><strong>Operating point:</strong> threshold selected on the tune split by maximizing F1, without the previous high-specificity constraint.</li>
        </ul>
        {image_block(workflow_fig, "Build and evaluation workflow for gdTAI.")}
      </section>

      <section id="data">
        <h2>Data And Validation Design</h2>
        <p>
          Gold labels were created from sorted gdT datasets and TCR evidence. Silver gdT cells were
          kept for sensitivity checks only. To reduce leakage risk, entire validation cohorts were
          held out before training and threshold tuning.
        </p>
        <ul>
          <li>All primary gold cells from <code>GSE144469</code> were held out.</li>
          <li>Paired-TCRAB/no-gdTCR cells from <code>GSE254249</code> were held out as alpha-beta T negatives.</li>
          <li><code>GDT_2020AUG_woCOV</code> cord-blood gdT cells were held out as an additional positive cohort.</li>
          <li><code>GDTlung2023july_7p</code> was excluded from train/tune because library quality was flagged as suboptimal.</li>
        </ul>
        <div class="table-wrap">{split_table}</div>
        <div class="grid two">
          {image_block(ATLAS_UMAP, "5-million-cell TNK atlas context: corrected simple annotation UMAP.", "small")}
          {image_block(ROOT / "Integrated_dataset/figures/gdT_prediction/umap_predicted_putative_gdt.png", "Whole-atlas UMAP highlighting putative gdT calls from the package-style score report.", "small")}
        </div>
      </section>

      <section id="features">
        <h2>Feature Selection</h2>
        <p>
          Feature selection was deterministic and did not use validation performance. Candidate
          genes were read from <code>var/_index</code>, prioritized as delta, gamma, alpha, then
          beta TCR genes, and capped at 300 genes. The current object contained 182 TCR genes plus
          5 penalty/control genes, so no selected feature was dropped by the cap.
        </p>
        <div class="grid two">
          {image_block(feature_fig, "Feature-family counts in the selected model.")}
          <div>
            <h3>Feature Families</h3>
            <div class="table-wrap">{feature_counts_table}</div>
            <h3>Representative Features</h3>
            <div class="table-wrap">{top_features_table}</div>
          </div>
        </div>
      </section>

      <section id="validation">
        <h2>How Good Is It?</h2>
        <p class="lead">
          The selected model made {fmt_int(selected["predicted_positive"])} positive calls in
          {fmt_int(selected["n_cells"])} validation cells. The confusion matrix was TP
          {fmt_int(selected["tp"])}, FP {fmt_int(selected["fp"])}, TN {fmt_int(selected["tn"])},
          and FN {fmt_int(selected["fn"])}.
        </p>
        <div class="grid two">
          {image_block(FIGURE_DIR / "validation_roc.png", "ROC curves on the held-out validation set.")}
          {image_block(FIGURE_DIR / "validation_pr.png", "Precision-recall curves on the held-out validation set.")}
        </div>
        <div class="grid two">
          {image_block(FIGURE_DIR / "selected_validation_confusion_matrix.png", "Selected gdTAI confusion matrix.")}
          {image_block(FIGURE_DIR / "selected_feature_importance.png", "Selected gdTAI feature importance.")}
        </div>
        <h3>Validation Metrics</h3>
        <div class="table-wrap">{validation_table}</div>
        <h3>Acceptance Against Individual-Gene Baselines</h3>
        <p>
          The selected gdTAI model improved F1 by
          <strong>{fmt_num(selected_acceptance["delta_f1"], 3)}</strong> over the best
          individual-gene heuristic baseline and met the F1 acceptance criterion.
        </p>
        <div class="table-wrap">{acceptance_table}</div>
        <h3>Validation By Held-Out Source</h3>
        <div class="table-wrap">{by_source_table}</div>
      </section>

      <section id="comparison">
        <h2>Comparison With The Original TRD-TRAB Score Strategy</h2>
        <p>
          On the full 5-million-cell atlas, gdTAI recovered {fmt_int(full_selected["predicted_putative_gdT"])}
          putative gdT cells versus {fmt_int(full_original["predicted_putative_gdT"])} from the original
          score rule. Because this version uses an F1-optimized threshold, the table below should be read
          as a precision-recall and false-positive tradeoff: paired-TCRAB FP relative change
          {fmt_pct(fp_tcrab_reduction, 1)}, NK FP relative change {fmt_pct(fp_nk_reduction, 1)}, and
          known-FP fraction relative change {fmt_pct(known_fp_reduction, 1)}.
        </p>
        {image_block(strategy_fig, "Full-atlas gdTAI versus original TRD-TRAB score strategy.")}
        <div class="table-wrap">{strategy_table}</div>
        <div class="grid two">
          {image_block(FIGURE_DIR / "trd_vs_trab_prediction_method_agreement.png", "TRD versus TRAB score scatter colored by agreement between gdTAI and the original score strategy.")}
          {image_block(FIGURE_DIR / "trd_vs_trab_tcrgene_known_fp_status.png", "TRD versus TRAB score scatter highlighting known false-positive compartments.")}
        </div>
      </section>

      <section id="full-atlas">
        <h2>Full 5-Million-Cell Atlas Application</h2>
        <p>
          Applying the selected threshold to all {fmt_int(full_selected["total_cells"])} cells in the
          5-million-cell atlas recovered {fmt_int(full_selected["predicted_putative_gdT"])} putative gdT cells.
          The tables below show where most recovered cells occur and where known false-positive
          audits remain important.
        </p>
        <h3>Top Sources By gdTAI Calls</h3>
        <div class="table-wrap">{source_top_table}</div>
        <h3>Top Tissues By gdTAI Calls</h3>
        <div class="table-wrap">{tissue_top_table}</div>
        <h3>Prediction Counts By Broad Annotation</h3>
        <div class="table-wrap">{annotation_table}</div>
      </section>

      <section id="prevalence">
        <h2>Class Imbalance And Prevalence-Aware PPV</h2>
        <p>
          gdT cells are rare, so the F1 score alone is not enough to judge usefulness. At low true
          prevalence, even a small false-positive rate can dominate the predicted-positive pool.
          This table shows the expected PPV and false positives per million under plausible prevalence
          values for the selected F1-optimized operating point.
        </p>
        <div class="table-wrap">{ppv_table}</div>
      </section>

      <section id="failure">
        <h2>Failure Modes And Overfitting Risk</h2>
        <p>
          Explicit doublet and QC labels were not available in the full atlas object, so the audit used
          expression-derived proxies. Failed positives were more likely to have FOXP3/CD4/low-CD3
          penalty-control proxy signals, lower gdT TCR-gene signal, and low-quality proxies. False positives were strongly
          enriched for dual TCR-gene and NK-like expression proxies rather than simple global
          low-quality signatures.
        </p>
        <div class="warning">
          Overfitting is reduced but not ruled out. Dataset-level holdouts, no TRD/TRAB score inputs,
          and no metadata features reduce leakage. Remaining risk comes from source-specific TCR
          capture, TCR chemistry differences, and labels that are partly derived from TCR metadata.
        </div>
        <div class="grid two">
          {image_block(FAILURE_FIGURE_DIR / "failure_mode_flag_fractions.png", "Failure-mode proxy flag fractions.")}
          {image_block(FAILURE_FIGURE_DIR / "failure_mode_marker_boxplots.png", "Failure-mode marker expression shifts.")}
        </div>
        <h3>Outcome-Level Failure Summary</h3>
        <div class="table-wrap">{failure_summary_table}</div>
        <h3>Top Proxy Enrichments</h3>
        <div class="table-wrap">{failure_flag_table}</div>
      </section>

      <section id="use">
        <h2>How Others Should Test gdTAI</h2>
        <p>
          The current model is small enough to share locally. The trusted pickle is
          {fmt_int(model_size)} bytes and has SHA256 checksum <code>{html.escape(model_sha)}</code>.
          The input H5AD is opened read-only by the wrapper.
        </p>
        <pre><code>cd {html.escape(str(ROOT))}

/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \\
  workflows/gdtai/predict_with_selected_gdt_model.py \\
  --input-h5ad /path/to/test_dataset.h5ad \\
  --model-pkl {html.escape(rel(MODEL_PATH))} \\
  --dataset-id TEST_DATASET \\
  --chunk-size 50000</code></pre>
        <ul>
          <li>Model: <code>{html.escape(rel(MODEL_PATH))}</code></li>
          <li>Usage protocol: <code>{html.escape(rel(USAGE_PROTOCOL))}</code></li>
          <li>Main output: <code>Integrated_dataset/tables/gdT_prediction/external_tests/TEST_DATASET/gdt_predictions.csv.gz</code></li>
          <li>Static output: <code>gdT_prediction/external_tests/TEST_DATASET/index.html</code></li>
        </ul>
      </section>

      <section id="sources">
        <h2>Source Files Used</h2>
        <ul class="captioned-list">
          <li><strong>Validation metrics:</strong> <code>{html.escape(rel(TABLE_DIR / "validation_metrics.csv"))}</code></li>
          <li><strong>Full-atlas application:</strong> <code>{html.escape(rel(TABLE_DIR / "selected_model_full_dataset_prediction_overall.csv"))}</code></li>
          <li><strong>Known FP audit:</strong> <code>{html.escape(rel(TABLE_DIR / "selected_model_full_dataset_false_positive_overall.csv"))}</code></li>
          <li><strong>Feature manifest:</strong> <code>{html.escape(rel(TABLE_DIR / "feature_manifest.csv"))}</code></li>
          <li><strong>Failure-mode audit:</strong> <code>{html.escape(rel(FAILURE_TABLE_DIR / "failure_mode_group_summary.csv"))}</code></li>
        </ul>
      </section>

      <footer>
        Generated by <code>build_gdtai_overview_report.py</code>. This report does not read or
        modify H5AD files; it summarizes existing gdTAI tables and figures.
      </footer>
    </main>
    """

    html_out = STATIC_DIR / "gdTAI_overview_report.html"
    html_doc = f"<!doctype html><html><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"><title>gdTAI Overview and Performance Report</title><style>{make_css()}</style></head><body>{html_body}</body></html>"
    html_out.write_text(html_doc, encoding="utf-8")

    md_out = LOG_DIR / "gdTAI_overview_report.md"
    md_text = f"""# gdTAI Overview And Performance Report

Generated: {report_time}

## Executive Summary

gdTAI is the current gdT AI classifier for the 5-million-cell TNK atlas we built. It uses
individual TCRA/TCRB/TCRG/TCRD gene expression plus FOXP3, CD4, and CD3D/CD3E/CD3G
penalty/control genes. It does not use Phase 4 TRD/TRAB scores, TCR metadata,
dataset ID, tissue, sample ID, or cluster labels as prediction features.

Main holdout performance for the selected gdTAI model:

- Precision: {fmt_pct(selected["precision"], 1)}
- Recall: {fmt_pct(selected["recall"], 1)}
- Specificity: {fmt_pct(selected["specificity"], 2)}
- F1: {fmt_num(selected["f1"], 3)}
- F0.5: {fmt_num(selected["f0.5"], 3)}
- MCC: {fmt_num(selected["mcc"], 3)}
- ROC-AUC: {fmt_num(selected["roc_auc"], 3)}
- PR-AUC: {fmt_num(selected["pr_auc"], 3)}
- Confusion matrix: TP {fmt_int(selected["tp"])}, FP {fmt_int(selected["fp"])}, TN {fmt_int(selected["tn"])}, FN {fmt_int(selected["fn"])}

Full 5-million-cell atlas application:

- Total cells: {fmt_int(full_selected["total_cells"])}
- gdTAI predicted gdT cells: {fmt_int(full_selected["predicted_putative_gdT"])} ({fmt_pct(full_selected["predicted_fraction"], 2)})
- Original TRD-TRAB strategy predicted gdT cells: {fmt_int(full_original["predicted_putative_gdT"])} ({fmt_pct(full_original["predicted_fraction"], 2)})
- gdTAI known FP fraction: {fmt_pct(fp_selected["known_FP_fraction_of_predictions"], 2)}
- Original TRD-TRAB known FP fraction: {fmt_pct(fp_original["known_FP_fraction_of_predictions"], 2)}

## Validation Splits

{markdown_table(split, integer_columns=["n_cells", "gdT_gold", "abT_gold"], percent_columns=["gdT_prevalence"])}

## Validation Metrics

{markdown_table(validation.assign(model_short=validation["model"].map(short_strategy)), columns=["model_short", "threshold", "n_cells", "predicted_positive", "tp", "fp", "tn", "fn", "precision", "recall", "specificity", "f1", "f0.5", "roc_auc", "pr_auc"], rename={"model_short": "model"}, integer_columns=["n_cells", "predicted_positive", "tp", "fp", "tn", "fn"], float_columns=["threshold", "precision", "recall", "specificity", "f1", "f0.5", "roc_auc", "pr_auc"])}

## gdTAI Versus Original Score Strategy

{markdown_table(strategy_compare, columns=["strategy_short", "threshold", "total_cells", "predicted_putative_gdT", "predicted_fraction", "predicted_paired_TCRAB_FP", "predicted_NK_FP", "predicted_paired_TCRAB_or_NK_FP", "known_FP_fraction_of_predictions"], rename={"strategy_short": "strategy"}, integer_columns=["total_cells", "predicted_putative_gdT", "predicted_paired_TCRAB_FP", "predicted_NK_FP", "predicted_paired_TCRAB_or_NK_FP"], percent_columns=["predicted_fraction", "known_FP_fraction_of_predictions"], float_columns=["threshold"])}

## Feature Families

{markdown_table(feature_counts, integer_columns=["n_features"])}

## Failure-Mode Summary

{markdown_table(failure_summary, columns=["outcome", "n_cells", "score_median", "source_low_quality_proxy_fraction", "multi_lineage_doublet_proxy_fraction", "tcr_dual_gene_proxy_fraction", "nk_like_proxy_fraction", "any_death_penalty_proxy_fraction"], integer_columns=["n_cells"], percent_columns=["source_low_quality_proxy_fraction", "multi_lineage_doublet_proxy_fraction", "tcr_dual_gene_proxy_fraction", "nk_like_proxy_fraction", "any_death_penalty_proxy_fraction"], float_columns=["score_median"])}

## Outputs

- HTML report: `{rel(html_out)}`
- Markdown summary: `{rel(md_out)}`
- Overview figures: `{rel(OVERVIEW_FIGURE_DIR)}`
- Model: `{rel(MODEL_PATH)}`
- Usage protocol: `{rel(USAGE_PROTOCOL)}`

## Interpretation

gdTAI is practical for F1-optimized gdT recovery in the TNK atlas when the goal
is to balance precision and recall while still auditing paired-TCRAB and NK false
positives. Dataset-level holdouts, no score inputs, and no metadata features
reduce leakage risk, but external validation remains the next important test.
"""
    md_out.write_text(md_text, encoding="utf-8")
    return html_out, md_out


def main() -> None:
    html_out, md_out = build_report()
    print(f"Wrote {rel(html_out)}")
    print(f"Wrote {rel(md_out)}")


if __name__ == "__main__":
    main()
