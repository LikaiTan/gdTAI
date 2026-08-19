#!/usr/bin/env python3
"""Render the frozen V4.5 versus V2/V3 consumed-benchmark report."""

from __future__ import annotations

import argparse
import html
import json
import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter


ROOT = Path(__file__).resolve().parents[2]
TABLES = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_5_vs_v2_v3_consumed"
LOGS = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_5_vs_v2_v3_consumed"
FIGURES = ROOT / "Integrated_dataset/figures/gdT_prediction/gdtai_v4_5_vs_v2_v3_consumed"
REPORT = ROOT / "gdT_prediction/gdtai_v4_5_vs_v2_v3_consumed"
ORDER = ["v4_5_highest_f1", "v4_5_high_purity", "v3_balanced", "v2_high_f1", "v2_high_purity"]
LABELS = {
    "v4_5_highest_f1": "V4.5 highest F1",
    "v4_5_high_purity": "V4.5 high purity",
    "v3_balanced": "V3 balanced",
    "v2_high_f1": "V2 highest F1",
    "v2_high_purity": "V2 high purity",
}
COLORS = {
    "v4_5_highest_f1": "#16817a",
    "v4_5_high_purity": "#9b3f52",
    "v3_balanced": "#176b87",
    "v2_high_f1": "#d18b32",
    "v2_high_purity": "#725a7a",
}


def savefig(fig: plt.Figure, name: str) -> Path:
    FIGURES.mkdir(parents=True, exist_ok=True)
    path = FIGURES / name
    fig.savefig(path, dpi=260, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def table_html(frame: pd.DataFrame, columns: list[str], labels: dict[str, str] | None = None,
               percentages: set[str] | None = None, decimals: set[str] | None = None) -> str:
    labels, percentages, decimals = labels or {}, percentages or set(), decimals or set()
    output = ["<table><thead><tr>"]
    output.extend(f"<th>{html.escape(labels.get(column, column))}</th>" for column in columns)
    output.append("</tr></thead><tbody>")
    for _, row in frame.iterrows():
        output.append("<tr>")
        for column in columns:
            value = row.get(column)
            if pd.isna(value):
                rendered = "-"
            elif column in percentages:
                rendered = f"{100 * float(value):.3f}%"
            elif column in decimals:
                rendered = f"{float(value):.4f}"
            elif isinstance(value, (int, np.integer)):
                rendered = f"{int(value):,}"
            elif isinstance(value, (float, np.floating)):
                rendered = f"{value:.3f}"
            else:
                rendered = str(value)
            output.append(f"<td>{html.escape(rendered)}</td>")
        output.append("</tr>")
    output.append("</tbody></table>")
    return "".join(output)


def plot_performance(metrics: pd.DataFrame, stratum: str, name: str, title: str) -> Path:
    view = metrics[metrics.stratum.eq(stratum)].set_index("model").loc[ORDER]
    x = np.arange(len(ORDER))
    fig, ax = plt.subplots(figsize=(10.2, 4.6))
    width = 0.24
    for index, metric in enumerate(["precision", "recall", "f1"]):
        values = view[metric].to_numpy()
        bars = ax.bar(x + (index - 1) * width, values, width,
                      color=["#5b8e7d", "#d18b32", "#176b87"][index], label=metric.title())
        if metric == "f1":
            ax.bar_label(bars, labels=[f"{value:.3f}" for value in values], fontsize=7.5, padding=2)
    ax.set_xticks(x, [LABELS[model] for model in ORDER], rotation=14, ha="right")
    ax.set_ylim(0, 1.03)
    ax.set_title(title, weight="bold")
    ax.legend(frameon=False, ncol=3)
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, name)


def plot_fpr(metrics: pd.DataFrame) -> Path:
    view = metrics[metrics.stratum.eq("v4_5_unexposed_consumed")].set_index("model").loc[ORDER]
    x = np.arange(len(ORDER))
    fig, ax = plt.subplots(figsize=(10.0, 4.7))
    width = 0.36
    ax.bar(x - width / 2, view.abt_fpr, width, color="#176b87", label="Alpha-beta FPR")
    ax.bar(x + width / 2, view.author_nk_fpr, width, color="#9b3f52", label="Author-NK FPR")
    ax.set_yscale("log")
    ax.set_xticks(x, [LABELS[model] for model in ORDER], rotation=14, ha="right")
    ax.set_ylabel("False-positive rate (log scale)")
    ax.set_title("V4.5-unexposed consumed negative controls", weight="bold")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.18)
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, "unexposed_false_positive_rates.png")


def plot_source_recall(per_source: pd.DataFrame) -> Path:
    sources = ["BALF_BLOOD_COPD", "GDT_2020AUG_woCOV", "GDTlung2023july_7p"]
    view = per_source[per_source.model.isin(ORDER) & per_source.source_gse_id.isin(sources)]
    pivot = view.pivot(index="source_gse_id", columns="model", values="gdt_recall").reindex(sources)
    x = np.arange(len(sources))
    fig, ax = plt.subplots(figsize=(10.2, 4.8))
    width = 0.16
    for index, model in enumerate(ORDER):
        values = pivot[model].to_numpy()
        ax.bar(x + (index - 2) * width, values, width, color=COLORS[model], label=LABELS[model])
    ax.set_xticks(x, ["BALF (unexposed, consumed)", "GDT2020 (training-exposed)", "GDTlung (training-exposed)"])
    ax.set_ylim(0, 1.03)
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.set_ylabel("gdT-gold recall")
    ax.set_title("Recall by positive benchmark source", weight="bold")
    ax.legend(frameon=False, ncol=3, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, "positive_source_recall.png")


def plot_negative_sources(per_source: pd.DataFrame) -> Path:
    view = per_source[per_source.n_abt_gold.gt(0) & per_source.model.isin(ORDER)].copy()
    pivot = view.pivot(index="source_gse_id", columns="model", values="abt_fpr")
    sources = pivot.max(axis=1).sort_values(ascending=False).index
    pivot = pivot.reindex(sources)
    fig, ax = plt.subplots(figsize=(10.5, 5.3))
    x = np.arange(len(pivot))
    width = 0.16
    for index, model in enumerate(ORDER):
        ax.bar(x + (index - 2) * width, pivot[model], width, color=COLORS[model], label=LABELS[model])
    ax.set_yscale("log")
    ax.set_xticks(x, pivot.index, rotation=25, ha="right")
    ax.set_ylabel("Alpha-beta FPR (log scale)")
    ax.set_title("False positives in each alpha-beta benchmark source", weight="bold")
    ax.legend(frameon=False, ncol=3, fontsize=8)
    ax.grid(axis="y", alpha=0.18)
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, "per_source_abt_fpr.png")


def plot_overlap(overlap: pd.DataFrame) -> Path:
    view = overlap[
        overlap.v4_5_mode.eq("v4_5_highest_f1")
        & overlap.comparator.isin(["v3_balanced", "v2_high_purity"])
        & overlap.truth_class.eq("gdT_gold")
    ].copy()
    view["comparison"] = view.comparator.map(LABELS)
    pivot = view.pivot(index="source_gse_id", columns="comparison", values="jaccard_positive").reindex(
        ["BALF_BLOOD_COPD", "GDT_2020AUG_woCOV", "GDTlung2023july_7p"]
    )
    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    x = np.arange(len(pivot))
    width = 0.36
    for index, comparator in enumerate(["V3 balanced", "V2 high purity"]):
        values = pivot[comparator].to_numpy()
        bars = ax.bar(x + (index - 0.5) * width, values, width,
                      color=[COLORS["v3_balanced"], COLORS["v2_high_purity"]][index], label=comparator)
        ax.bar_label(bars, labels=[f"{value:.3f}" for value in values], fontsize=8, padding=2)
    ax.set_xticks(x, ["BALF", "GDT2020", "GDTlung"])
    ax.set_ylim(0, 1.03)
    ax.set_ylabel("Positive-call Jaccard")
    ax.set_title("V4.5 highest-F1 prediction overlap", weight="bold")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, "positive_prediction_overlap.png")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-pdf", action="store_true")
    args = parser.parse_args()
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 12, "axes.labelsize": 10})
    metrics = pd.read_csv(TABLES / "stratified_overall_metrics.csv")
    per_source = pd.read_csv(TABLES / "per_source_metrics.csv")
    overlap = pd.read_csv(TABLES / "per_source_truth_prediction_overlap.csv")
    exposure = pd.read_csv(TABLES / "benchmark_exposure_manifest.csv")
    integrity = pd.read_csv(TABLES / "implementation_integrity.csv")
    summary = json.loads((LOGS / "evaluation_summary.json").read_text())

    figures = {
        "pooled": plot_performance(metrics, "all_consumed_benchmark", "pooled_performance.png",
                                   "Pooled consumed benchmark: biased by V4.5 training exposure"),
        "unexposed": plot_performance(metrics, "v4_5_unexposed_consumed", "unexposed_performance.png",
                                      "V4.5-unexposed but previously consumed benchmark"),
        "fpr": plot_fpr(metrics),
        "recall": plot_source_recall(per_source),
        "negative": plot_negative_sources(per_source),
        "overlap": plot_overlap(overlap),
    }
    REPORT.mkdir(parents=True, exist_ok=True)
    assets = REPORT / "assets"
    assets.mkdir(parents=True, exist_ok=True)
    for source in figures.values():
        shutil.copy2(source, assets / source.name)

    def metric_view(stratum: str) -> pd.DataFrame:
        view = metrics[metrics.stratum.eq(stratum)].set_index("model").loc[ORDER].reset_index()
        view["model"] = view.model.map(LABELS)
        return view[["model", "predicted_positive", "precision", "recall", "f1", "abt_fpr", "author_nk_fpr"]]

    pooled = metric_view("all_consumed_benchmark")
    unexposed = metric_view("v4_5_unexposed_consumed")
    balf = metric_view("balf_unexposed_consumed")
    exposed = metric_view("v4_5_training_exposed_sorted")[["model", "predicted_positive", "recall", "f1"]]
    positive_sources = per_source[per_source.n_gdt_gold.gt(0) & per_source.model.isin(ORDER)][[
        "source_gse_id", "v4_5_exposure", "model", "n_gdt_gold", "gdt_recall"
    ]].copy()
    positive_sources["model"] = positive_sources.model.map(LABELS)
    negative_sources = per_source[
        per_source.n_abt_gold.gt(0) & per_source.model.isin(["v4_5_highest_f1", "v4_5_high_purity", "v3_balanced", "v2_high_purity"])
    ][["source_gse_id", "model", "n_abt_gold", "abt_fpr", "n_author_nk", "author_nk_fpr"]].copy()
    negative_sources["model"] = negative_sources.model.map(LABELS)
    exposure_view = exposure.groupby(["v4_5_exposure", "source_gse_id"], observed=True).n_cells.sum().reset_index()

    css = """
    :root{--ink:#17272d;--muted:#52656c;--line:#ccd8dc;--teal:#16817a;--red:#9b3f52;--blue:#176b87}
    *{box-sizing:border-box}body{margin:0;background:#f3f6f7;color:var(--ink);font-family:Arial,Helvetica,sans-serif;line-height:1.45}
    header{background:#17272d;color:white;padding:34px 42px 28px}header h1{font-size:28px;margin:0 0 7px;letter-spacing:0}header p{margin:0;color:#dbe6e8}
    main{max-width:1180px;margin:auto;background:white}section{padding:27px 40px;border-bottom:1px solid var(--line)}
    h2{font-size:20px;margin:0 0 13px;color:#174453;letter-spacing:0}h3{font-size:15px;margin:17px 0 7px;letter-spacing:0}p{margin:8px 0}
    .verdict{border-left:5px solid var(--red);background:#fff6f7;padding:12px 15px}.note{border-left-color:#d18b32;background:#fff9ef}.ok{border-left-color:var(--teal);background:#f2fbfa}
    .grid{display:grid;grid-template-columns:1fr 1fr;gap:22px;align-items:start}figure{margin:12px 0 5px}figure img{width:100%;height:auto;display:block}
    figcaption,.small{font-size:9.3pt;color:var(--muted)}table{width:100%;border-collapse:collapse;font-size:8.6pt;margin:10px 0 15px;table-layout:fixed}
    th{background:#e8f0f2;color:#174453;text-align:left}th,td{border:1px solid var(--line);padding:5px 6px;vertical-align:top;overflow-wrap:anywhere}tr:nth-child(even) td{background:#f8fafb}
    code{font-family:Consolas,monospace;font-size:.9em}ul{padding-left:19px;margin:8px 0}
    @media(max-width:780px){.grid{grid-template-columns:1fr}section{padding:21px}header{padding:27px 21px}}
    @media print{@page{size:A4 landscape;margin:9mm}body{background:white;font-size:9pt}header{padding:17px 23px}header h1{font-size:21px}main{max-width:none}
      section{padding:13px 20px}h2{font-size:15px}h3{font-size:11.5px}figure{break-inside:avoid}figure img{max-height:121mm;object-fit:contain}
      table{font-size:6.8pt;table-layout:fixed}th,td{padding:3.4px}tr{break-inside:avoid}.page{break-before:page}.grid{gap:12px}}
    """
    body = f"""<!doctype html><html><head><meta charset='utf-8'><title>gdTAI V4.5 versus V2/V3</title><style>{css}</style></head><body>
    <header><h1>gdTAI V4.5 versus V2/V3</h1><p>Frozen post-hoc comparison on the consumed 335,479-cell common benchmark | 19 August 2026</p></header><main>
    <section><h2>Executive conclusion</h2><div class='verdict'><strong>Do not promote V4.5.</strong> Its pooled highest-F1 result is inflated by 41,079 training-exposed sorted positives. On V4.5-unexposed consumed cells, highest-F1 has excessive alpha-beta FPR, while high-purity does not clearly outperform V3 balanced or V2 high-purity.</div>
    <p>V4.5 highest-F1 pooled F1 is 0.875, but its unexposed alpha-beta FPR is <strong>0.923%</strong> versus <strong>0.146%</strong> for V3 and <strong>0.099%</strong> for V2 high-purity. V4.5 high-purity reduces alpha-beta FPR to <strong>0.189%</strong>, but unexposed F1 is 0.663 versus 0.665 for V3 and 0.753 for V2 high-purity.</p></section>
    <section><h2>Validity boundary</h2><div class='grid'><div><div class='verdict note'><strong>Not a new external test.</strong> Every benchmark cohort was previously inspected. GDT2020 and GDTlung were additionally used by V4.5 development.</div><ul><li>Training-exposed for V4.5: GDT2020 and GDTlung.</li><li>Unexposed to V4.5 but consumed previously: BALF and eight extension sources.</li><li>No threshold, feature, exclusion, or calibrator changed.</li><li>No result in this report is eligible for promotion.</li></ul></div><div>{table_html(exposure_view, list(exposure_view.columns), {'v4_5_exposure':'V4.5 exposure','source_gse_id':'Source','n_cells':'Cells'})}</div></div></section>
    <section class='page'><h2>Pooled benchmark</h2><div class='verdict note'>Pooled metrics mix training-exposed and unexposed cells. They describe prediction behavior but cannot rank generalization.</div><figure><img src='assets/{figures['pooled'].name}'><figcaption>V4.5 highest-F1 benefits from its training exposure to 41,079 of 41,931 benchmark positives.</figcaption></figure>
    {table_html(pooled, list(pooled.columns), {'model':'Model','predicted_positive':'Predicted','abt_fpr':'alpha-beta FPR','author_nk_fpr':'NK FPR'}, percentages={'precision','recall','abt_fpr','author_nk_fpr'}, decimals={'f1'})}</section>
    <section class='page'><h2>V4.5-unexposed consumed cells</h2><div class='grid'><figure><img src='assets/{figures['unexposed'].name}'><figcaption>This subset contains BALF as the only positive source plus all alpha-beta/NK negative extension controls.</figcaption></figure><figure><img src='assets/{figures['fpr'].name}'><figcaption>Highest-F1 is not sufficiently specific for the very low gdT prevalence. High-purity approaches V3/V2 specificity but does not dominate them.</figcaption></figure></div>
    {table_html(unexposed, list(unexposed.columns), {'model':'Model','predicted_positive':'Predicted','abt_fpr':'alpha-beta FPR','author_nk_fpr':'NK FPR'}, percentages={'precision','recall','abt_fpr','author_nk_fpr'}, decimals={'f1'})}</section>
    <section class='page'><h2>BALF positive benchmark</h2><p>BALF was not fitted by V4.5, but it was already inspected during earlier model rounds. It is diagnostic, not external.</p>
    {table_html(balf, list(balf.columns), {'model':'Model','predicted_positive':'Predicted','abt_fpr':'alpha-beta FPR','author_nk_fpr':'NK FPR'}, percentages={'precision','recall','abt_fpr','author_nk_fpr'}, decimals={'f1'})}
    <p>V4.5 highest-F1 has the best BALF F1 (0.949) and recall 94.01%, with 0.129% alpha-beta FPR and zero observed author-NK false positives. V4.5 high-purity reaches 99.21% precision but recalls 88.73%.</p></section>
    <section class='page'><h2>Positive sources and exposure effect</h2><figure><img src='assets/{figures['recall'].name}'><figcaption>Only BALF is V4.5-unexposed. Higher V4.5 recall in GDT2020/GDTlung is expected after those cohorts entered development.</figcaption></figure>
    {table_html(positive_sources, list(positive_sources.columns), {'source_gse_id':'Source','v4_5_exposure':'V4.5 exposure','model':'Model','n_gdt_gold':'gdT gold','gdt_recall':'Recall'}, percentages={'gdt_recall'})}
    <h3>Training-exposed pooled recall</h3>{table_html(exposed, list(exposed.columns), {'model':'Model','predicted_positive':'Predicted','recall':'Recall'}, percentages={'recall'}, decimals={'f1'})}</section>
    <section class='page'><h2>Per-dataset false positives</h2><figure><img src='assets/{figures['negative'].name}'><figcaption>GSE159251 dominates V4.5 false positives: 6.91% for highest-F1 and 1.27% for high-purity, versus 0.78% for V3 and 0.51% for V2 high-purity.</figcaption></figure>
    {table_html(negative_sources, list(negative_sources.columns), {'source_gse_id':'Source','model':'Model','n_abt_gold':'alpha-beta gold','abt_fpr':'alpha-beta FPR','n_author_nk':'Author NK','author_nk_fpr':'NK FPR'}, percentages={'abt_fpr','author_nk_fpr'})}</section>
    <section class='page'><h2>Prediction overlap</h2><figure><img src='assets/{figures['overlap'].name}'><figcaption>V4.5 agrees strongly with V2/V3 on BALF and GDT2020. GDTlung overlap is lower, showing that V4.5 gains there come from a distinct training-exposed subset.</figcaption></figure>
    <p>For BALF gdT gold, V4.5 highest-F1 finds 46 cells missed by V3 but misses 11 found by V3. In GDTlung it finds 3,088 cells beyond V3 while missing 351 V3 calls; this difference cannot be interpreted as external improvement because GDTlung trained V4.5.</p></section>
    <section class='page'><h2>Model decision</h2><ul><li><strong>V3 balanced remains the default.</strong> It retains the strongest defensible balance across the consumed benchmark without V4.5's new training exposure.</li><li><strong>V2 high-purity remains the conservative fallback.</strong> It has the lowest unexposed alpha-beta FPR and the best unexposed consumed F1.</li><li><strong>V4.5 highest-F1 is not deployment-ready.</strong> Its GSE159251 false-positive rate is unacceptable for a rare-cell classifier.</li><li><strong>V4.5 high-purity is promising on BALF</strong> but lacks a reproducible advantage over V3/V2 across negative sources.</li><li>No V4.5 threshold should be changed using this report. A new untouched positive and negative benchmark is required.</li></ul></section>
    <section class='page'><h2>Integrity and reproducibility</h2>{table_html(integrity, ['check','pass','detail'], {'check':'Check','pass':'Pass','detail':'Evidence'})}
    <div class='verdict ok'><strong>Frozen diagnostic passed.</strong> No fitting, calibration, feature selection, or threshold selection occurred; no H5AD was modified and no model was promoted or pushed.</div>
    <p class='small'>Evaluator: <code>workflows/gdtai/evaluate_gdtai_v4_5_vs_v2_v3_consumed.py</code><br>V4.5 contract SHA-256: <code>{summary['v4_5_model_contract_sha256']}</code><br>Input V2/V3 predictions SHA-256: <code>{summary['input_predictions_sha256']}</code><br>Output predictions SHA-256: <code>{summary['output_predictions_sha256']}</code></p></section>
    </main></body></html>"""
    html_path = REPORT / "index.html"
    html_path.write_text(body)
    pdf_path = REPORT / "gdtai_v4_5_vs_v2_v3_consumed_report.pdf"
    if not args.skip_pdf:
        subprocess.run([
            "google-chrome", "--headless", "--no-sandbox", "--disable-gpu",
            "--run-all-compositor-stages-before-draw", "--virtual-time-budget=3000",
            f"--print-to-pdf={pdf_path}", html_path.resolve().as_uri(),
        ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elif not pdf_path.exists():
        raise FileNotFoundError(f"--skip-pdf requested but PDF is absent: {pdf_path}")

    (LOGS / "comparison_summary.md").write_text(f"""# gdTAI V4.5 versus V2/V3

## Decision

V4.5 is not promoted. Pooled improvement is confounded by training exposure.

- V4.5 highest-F1 unexposed alpha-beta FPR: `0.009232`
- V4.5 high-purity unexposed F1: `0.6629`
- V3 balanced unexposed F1: `0.6649`
- V2 high-purity unexposed F1: `0.7527`
- HTML: `{html_path}`
- PDF: `{pdf_path}`
""")
    print(json.dumps({"status": "PASS_REPORT", "html": str(html_path), "pdf": str(pdf_path),
                      "pdf_bytes": pdf_path.stat().st_size, "figures": len(figures)}, indent=2))


if __name__ == "__main__":
    main()
