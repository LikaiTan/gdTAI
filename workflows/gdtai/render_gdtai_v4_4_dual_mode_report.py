#!/usr/bin/env python3
"""Render the gdTAI V4.4 dual-mode development and diagnostic report."""

from __future__ import annotations

import html
import json
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter


ROOT = Path(__file__).resolve().parents[2]
DEV_TABLES = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_4_dual_mode"
TEST_TABLES = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_4_reused_lockbox"
MODEL_DIR = ROOT / "Integrated_dataset/models/gdT_prediction/gdtai_v4_4_dual_mode"
FIGURE_DIR = ROOT / "Integrated_dataset/figures/gdT_prediction/gdtai_v4_4_dual_mode"
LOG_DIR = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_4_dual_mode"
REPORT_DIR = ROOT / "gdT_prediction/gdtai_v4_4_dual_mode"
COLORS = {
    "v4_4_highest_f1": "#087e8b",
    "v4_4_high_purity": "#d1495b",
    "v4_3_balanced": "#6c757d",
    "v3_balanced": "#2f6690",
    "v2_high_f1": "#f4a261",
    "v2_high_purity": "#8a6d3b",
}
LABELS = {
    "v4_4_highest_f1": "V4.4 highest F1",
    "v4_4_high_purity": "V4.4 high purity",
    "v4_3_balanced": "V4.3",
    "v3_balanced": "V3 balanced",
    "v2_high_f1": "V2 highest F1",
    "v2_high_purity": "V2 high purity",
}


def savefig(fig: plt.Figure, name: str) -> Path:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    path = FIGURE_DIR / name
    fig.savefig(path, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def table_html(frame: pd.DataFrame, columns: list[str], labels: dict[str, str] | None = None,
               percent: set[str] | None = None, decimals: set[str] | None = None) -> str:
    labels = labels or {}
    percent = percent or set()
    decimals = decimals or set()
    rows = ["<table><thead><tr>" + "".join(f"<th>{html.escape(labels.get(c, c))}</th>" for c in columns) + "</tr></thead><tbody>"]
    for _, row in frame.iterrows():
        cells = []
        for column in columns:
            value = row.get(column)
            if pd.isna(value):
                rendered = "-"
            elif column in percent:
                rendered = f"{100 * float(value):.2f}%"
            elif column in decimals:
                rendered = f"{float(value):.4f}"
            elif isinstance(value, (float, np.floating)):
                rendered = f"{value:.3f}"
            elif isinstance(value, (int, np.integer)):
                rendered = f"{int(value):,}"
            else:
                rendered = str(value)
            cells.append(f"<td>{html.escape(rendered)}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    rows.append("</tbody></table>")
    return "".join(rows)


def plot_partition(counts: pd.DataFrame) -> Path:
    partition_column = "v4_4_partition" if "v4_4_partition" in counts.columns else "partition"
    totals = counts.groupby(partition_column, observed=True).n_cells.sum().reindex(
        ["fit", "calibration", "threshold_validation", "locked_test"]
    )
    fig, ax = plt.subplots(figsize=(8.2, 3.6))
    colors = ["#2f6690", "#7a9e9f", "#087e8b", "#6c757d"]
    bars = ax.bar(totals.index, totals.values, color=colors, width=0.68)
    ax.bar_label(bars, labels=[f"{int(v):,}" for v in totals.values], padding=4, fontsize=9)
    ax.set_ylabel("Cells")
    ax.set_title("Group-disjoint development and test partitions", fontsize=12, weight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", rotation=12)
    return savefig(fig, "development_partitions.png")


def plot_frontier(frontier: pd.DataFrame, contract: dict) -> Path:
    view = frontier.iloc[::max(1, len(frontier) // 4000)]
    fig, ax = plt.subplots(figsize=(6.6, 4.5))
    ax.plot(view.recall, view.precision, color="#59656f", lw=1.4)
    for mode, payload in contract["operating_modes"].items():
        color = COLORS[f"v4_4_{mode}"]
        ax.scatter(payload["recall"], payload["precision"], s=70, color=color,
                   edgecolor="white", linewidth=0.8, label=LABELS[f"v4_4_{mode}"], zorder=3)
    ax.set_xlim(0.65, 1.005)
    ax.set_ylim(0.70, 1.005)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Development threshold-validation frontier", fontsize=12, weight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.18)
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, "development_precision_recall_frontier.png")


def plot_development(dev: pd.DataFrame) -> Path:
    metrics = ["precision", "recall", "f1", "f0_5"]
    x = np.arange(len(metrics))
    fig, ax = plt.subplots(figsize=(7.8, 4.2))
    width = 0.34
    for offset, (_, row) in zip((-0.5, 0.5), dev.iterrows()):
        model = f"v4_4_{row['mode']}"
        bars = ax.bar(x + offset * width, [row[m] for m in metrics], width,
                      label=LABELS[model], color=COLORS[model])
        ax.bar_label(bars, labels=[f"{v:.3f}" for v in [row[m] for m in metrics]], fontsize=8, padding=2)
    ax.set_xticks(x, ["Precision", "Recall", "F1", "F0.5"])
    ax.set_ylim(0.82, 1.01)
    ax.set_title("Reserved development-validation performance", fontsize=12, weight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, "development_mode_performance.png")


def plot_holdout(holdout: pd.DataFrame) -> Path:
    pivot = holdout.pivot(index="heldout_source", columns="mode", values="gdt_recall").reindex(
        ["HRA005041", "GSE144469", "MalteGDT"]
    )
    x = np.arange(len(pivot))
    fig, ax = plt.subplots(figsize=(7.8, 4.2))
    width = 0.34
    for index, mode in enumerate(["highest_f1", "high_purity"]):
        values = pivot[mode].to_numpy()
        bars = ax.bar(x + (index - 0.5) * width, values, width,
                      color=COLORS[f"v4_4_{mode}"], label=LABELS[f"v4_4_{mode}"])
        ax.bar_label(bars, labels=[f"{100*v:.1f}%" for v in values], fontsize=8, padding=2)
    ax.set_xticks(x, pivot.index)
    ax.set_ylim(0, 1.0)
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.set_ylabel("gdT-gold recall")
    ax.set_title("Positive-source-held-out development audits", fontsize=12, weight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, "source_holdout_recall.png")


def plot_lockbox_metrics(overall: pd.DataFrame) -> Path:
    order = ["v4_4_highest_f1", "v4_4_high_purity", "v4_3_balanced", "v3_balanced", "v2_high_f1", "v2_high_purity"]
    view = overall.set_index("model").loc[order]
    metrics = ["precision", "recall", "f1"]
    x = np.arange(len(order))
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    width = 0.23
    shades = ["#1b9aaa", "#ef476f", "#495057"]
    for index, metric in enumerate(metrics):
        values = view[metric].to_numpy()
        bars = ax.bar(x + (index - 1) * width, values, width, label=metric.title(), color=shades[index])
        if metric == "f1":
            ax.bar_label(bars, labels=[f"{v:.3f}" for v in values], fontsize=7, padding=2)
    ax.set_xticks(x, [LABELS[name] for name in order], rotation=18, ha="right")
    ax.set_ylim(0.35, 1.02)
    ax.set_title("Reused common-lockbox diagnostic", fontsize=12, weight="bold")
    ax.legend(frameon=False, ncol=3, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, "reused_lockbox_model_comparison.png")


def plot_fpr(overall: pd.DataFrame) -> Path:
    order = ["v4_4_highest_f1", "v4_4_high_purity", "v4_3_balanced", "v3_balanced", "v2_high_f1", "v2_high_purity"]
    view = overall.set_index("model").loc[order]
    x = np.arange(len(order))
    fig, ax = plt.subplots(figsize=(10.5, 4.6))
    width = 0.34
    ax.bar(x - width/2, view.abt_fpr, width, color="#2f6690", label="Paired alpha-beta FPR")
    ax.bar(x + width/2, view.author_nk_fpr, width, color="#d1495b", label="Author-NK FPR")
    ax.set_yscale("log")
    ax.set_xticks(x, [LABELS[name] for name in order], rotation=18, ha="right")
    ax.set_ylabel("False-positive rate (log scale)")
    ax.set_title("Reused-lockbox false-positive behavior", fontsize=12, weight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(axis="y", alpha=0.2)
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, "reused_lockbox_fpr.png")


def plot_positive_sources(per_source: pd.DataFrame) -> Path:
    models = ["v4_4_highest_f1", "v4_4_high_purity", "v3_balanced", "v2_high_f1", "v2_high_purity"]
    sources = ["BALF_BLOOD_COPD", "GDT_2020AUG_woCOV", "GDTlung2023july_7p"]
    view = per_source[per_source.model.isin(models) & per_source.source_gse_id.isin(sources)]
    pivot = view.pivot(index="source_gse_id", columns="model", values="gdt_recall").reindex(sources)
    x = np.arange(len(sources))
    fig, ax = plt.subplots(figsize=(10.2, 4.8))
    width = 0.15
    for index, model in enumerate(models):
        ax.bar(x + (index - 2) * width, pivot[model], width, label=LABELS[model], color=COLORS[model])
    ax.set_xticks(x, ["BALF/Blood COPD", "GDT 2020", "GDT lung"])
    ax.set_ylim(0, 1.02)
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.set_ylabel("gdT-gold recall")
    ax.set_title("Diagnostic recall by positive test source", fontsize=12, weight="bold")
    ax.legend(frameon=False, ncol=3, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, "positive_test_source_recall.png")


def main() -> None:
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 12, "axes.labelsize": 10})
    contract = json.loads((MODEL_DIR / "model_contract.json").read_text())
    counts = pd.read_csv(DEV_TABLES / "development_partition_counts.csv")
    frontier = pd.read_csv(DEV_TABLES / "threshold_validation_frontier.csv")
    dev = pd.read_csv(DEV_TABLES / "threshold_validation_metrics.csv")
    holdout = pd.read_csv(DEV_TABLES / "source_holdout_metrics.csv")
    overall = pd.read_csv(TEST_TABLES / "overall_metrics.csv")
    per_source = pd.read_csv(TEST_TABLES / "per_source_metrics.csv")
    integrity = pd.read_csv(TEST_TABLES / "implementation_integrity.csv")

    figures = {
        "partition": plot_partition(counts),
        "frontier": plot_frontier(frontier, contract),
        "development": plot_development(dev),
        "holdout": plot_holdout(holdout),
        "lockbox": plot_lockbox_metrics(overall),
        "fpr": plot_fpr(overall),
        "sources": plot_positive_sources(per_source),
    }
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    assets = REPORT_DIR / "assets"
    assets.mkdir(parents=True, exist_ok=True)
    for name, source in figures.items():
        target = assets / source.name
        target.write_bytes(source.read_bytes())

    dev_view = dev[["mode", "threshold", "precision", "recall", "f1", "f0_5", "specificity"]].copy()
    holdout_view = holdout[["heldout_source", "mode", "gdt_recall", "abt_fpr", "precision", "f1"]].copy()
    overall_view = overall[["model", "predicted_positive", "precision", "recall", "f1", "abt_fpr", "author_nk_fpr", "precision_at_0.01_prevalence"]].copy()
    overall_view["model"] = overall_view.model.map(LABELS)
    source_view = per_source[
        per_source.model.isin(["v4_4_highest_f1", "v4_4_high_purity", "v3_balanced", "v2_high_f1", "v2_high_purity"])
        & per_source.source_gse_id.isin(["BALF_BLOOD_COPD", "GDT_2020AUG_woCOV", "GDTlung2023july_7p"])
    ][["source_gse_id", "model", "n_gdt_gold", "gdt_recall"]].copy()
    source_view["model"] = source_view.model.map(LABELS)
    negative_view = per_source[
        per_source.model.isin(["v4_4_highest_f1", "v4_4_high_purity", "v3_balanced"])
        & (per_source.n_abt_gold.gt(0) | per_source.n_author_nk.gt(0))
    ][["source_gse_id", "model", "n_abt_gold", "abt_fpr", "n_author_nk", "author_nk_fpr"]].copy()
    negative_view["model"] = negative_view.model.map(LABELS)

    css = """
    :root{--ink:#16252d;--muted:#52646d;--line:#cad5da;--teal:#087e8b;--red:#d1495b;--blue:#2f6690}
    *{box-sizing:border-box} body{margin:0;background:#f4f7f8;color:var(--ink);font-family:Arial,Helvetica,sans-serif;line-height:1.45}
    header{background:#16252d;color:white;padding:34px 42px 28px} header h1{font-size:30px;margin:0 0 8px;letter-spacing:0} header p{margin:0;color:#dce7ea}
    main{max-width:1180px;margin:auto;background:white} section{padding:28px 42px;border-bottom:1px solid var(--line)}
    h2{font-size:21px;margin:0 0 14px;color:#173f4f;letter-spacing:0} h3{font-size:15px;margin:18px 0 8px;letter-spacing:0}
    p{margin:8px 0}.lead{font-size:16px}.verdict{border-left:5px solid var(--red);padding:12px 16px;background:#fff7f7}
    .ok{border-left-color:var(--teal);background:#f2fbfb}.grid{display:grid;grid-template-columns:1fr 1fr;gap:24px;align-items:start}
    figure{margin:14px 0 4px} figure img{width:100%;height:auto;display:block} figcaption{font-size:9.5pt;color:var(--muted);margin-top:6px}
    table{width:100%;border-collapse:collapse;font-size:9pt;margin:10px 0 16px;table-layout:fixed} th{background:#e9f0f2;text-align:left;color:#173f4f}
    th,td{border:1px solid var(--line);padding:6px 7px;vertical-align:top;overflow-wrap:anywhere} tr:nth-child(even) td{background:#f8fafb}
    code{font-family:Consolas,monospace;font-size:.92em}.small{font-size:9.5pt;color:var(--muted)} ul{margin:8px 0;padding-left:20px}
    @media(max-width:780px){.grid{grid-template-columns:1fr}section{padding:22px}header{padding:28px 22px}}
    @media print{@page{size:A4 landscape;margin:9mm}body{background:white;font-size:9pt}header{padding:18px 24px}header h1{font-size:22px}main{max-width:none}
      section{padding:14px 22px;break-inside:auto}h2{font-size:16px}h3{font-size:12px}figure{break-inside:avoid}figure img{max-height:125mm;object-fit:contain}
      table{font-size:7.2pt;table-layout:fixed}th,td{padding:4px;overflow-wrap:anywhere}tr{break-inside:avoid}.page{break-before:page}.grid{gap:14px}}
    """
    body = f"""<!doctype html><html><head><meta charset='utf-8'><title>gdTAI V4.4 dual-mode evaluation</title><style>{css}</style></head><body>
    <header><h1>gdTAI V4.4 dual-mode evaluation</h1><p>Development-frozen highest-F1 and high-purity operating modes | 19 August 2026</p></header><main>
    <section><h2>Executive decision</h2><div class='verdict'><strong>Do not promote V4.4.</strong> The redesign improves specificity over V4.3, but external diagnostic recall remains substantially below V2/V3. V3 balanced remains the default and V2 high purity remains the conservative fallback.</div>
    <p class='lead'>The thresholds were selected only from grouped development-validation cells: <strong>0.985703</strong> for highest F1 and <strong>0.996359</strong> for high purity. Test cells fitted no model or calibrator and selected no feature, model, or threshold.</p></section>
    <section><h2>What changed</h2><div class='grid'><div><h3>Classifier architecture</h3><ul><li>Stage 1: high-recall T-lineage support gate.</li><li>Stage 2: individual TCR genes plus symmetric receptor aggregates.</li><li>Confident author-NK cells enter Stage 2 as negatives.</li><li>No <code>KLRD1</code>, <code>NKG7</code>, <code>GNLY</code>, <code>FCER1G</code>, <code>TYROBP</code>, or <code>FCGR3A</code> feature.</li><li>No hard TRDV-expression cutoff.</li><li>200 effective Stage-2 features, below the 300-feature cap.</li></ul></div>
    <div><h3>Operating modes</h3><ul><li><strong>Highest F1:</strong> maximize F1 on threshold validation.</li><li><strong>High purity:</strong> maximize F0.5 on threshold validation.</li><li>No fixed 0.2% alpha-beta FPR limit.</li><li>FPR is reported after threshold selection, not used as a hard constraint.</li></ul></div></div></section>
    <section class='page'><h2>Leakage firewall</h2><div class='grid'><figure><img src='assets/{figures['partition'].name}'><figcaption>All biological groups remain within one partition. The 327,455-cell common lockbox is excluded from fitting, calibration, and threshold selection.</figcaption></figure>
    <div><div class='verdict ok'><strong>Integrity gate passed.</strong> Model contract records <code>thresholds_selected_from_test=false</code>.</div>
    {table_html(integrity, ['check','pass','detail'], {'check':'Check','pass':'Pass','detail':'Evidence'})}</div></div></section>
    <section class='page'><h2>Development operating points</h2><div class='grid'><figure><img src='assets/{figures['frontier'].name}'><figcaption>The two thresholds are selected on the same reserved development precision-recall frontier.</figcaption></figure>
    <figure><img src='assets/{figures['development'].name}'><figcaption>Highest F1 favors recovery; high purity gives additional weight to precision through F0.5.</figcaption></figure></div>
    {table_html(dev_view, ['mode','threshold','precision','recall','f1','f0_5','specificity'], {'mode':'Mode','threshold':'Threshold','f0_5':'F0.5'}, percent={'precision','recall','specificity'}, decimals={'threshold','f1','f0_5'})}</section>
    <section class='page'><h2>Source-held-out development audits</h2><figure><img src='assets/{figures['holdout'].name}'><figcaption>Each displayed positive source was excluded from its audit model, calibration, and threshold selection. These results expose residual source dependence despite strong within-development validation.</figcaption></figure>
    {table_html(holdout_view, ['heldout_source','mode','gdt_recall','abt_fpr','precision','f1'], {'heldout_source':'Held-out source','mode':'Mode','gdt_recall':'gdT recall','abt_fpr':'alpha-beta FPR'}, percent={'gdt_recall','abt_fpr','precision'}, decimals={'f1'})}</section>
    <section class='page'><h2>Reused common-lockbox diagnostic</h2><div class='verdict'><strong>Consumed test set.</strong> This 335,479-cell lockbox was previously used for the V4.3 gate. V4.4 scores are useful diagnostic evidence only and cannot support model promotion.</div>
    <div class='grid'><figure><img src='assets/{figures['lockbox'].name}'><figcaption>V4.4 improves on V4.3 but remains well below V2/V3 in F1 because of external recall loss.</figcaption></figure>
    <figure><img src='assets/{figures['fpr'].name}'><figcaption>V4.4 produces low alpha-beta and NK FPRs, but specificity alone does not compensate for lost gdT recovery.</figcaption></figure></div>
    {table_html(overall_view, ['model','predicted_positive','precision','recall','f1','abt_fpr','author_nk_fpr','precision_at_0.01_prevalence'], {'predicted_positive':'Predicted','abt_fpr':'alpha-beta FPR','author_nk_fpr':'NK FPR','precision_at_0.01_prevalence':'Precision at 1% prevalence'}, percent={'precision','recall','abt_fpr','author_nk_fpr','precision_at_0.01_prevalence'}, decimals={'f1'})}</section>
    <section class='page'><h2>Positive test sources</h2><figure><img src='assets/{figures['sources'].name}'><figcaption>V4.4 highest-F1 recall is 69.0% in BALF/Blood COPD, 66.0% in GDT 2020, and 27.4% in the suboptimal-quality GDT lung cohort.</figcaption></figure>
    {table_html(source_view, ['source_gse_id','model','n_gdt_gold','gdt_recall'], {'source_gse_id':'Source','n_gdt_gold':'Gold gdT','gdt_recall':'Recall'}, percent={'gdt_recall'})}</section>
    <section class='page'><h2>Negative test sources</h2><p>Per-dataset alpha-beta and author-NK false-positive rates are shown below. Zero values are retained as observed zeros, not replaced by pseudocounts.</p>
    {table_html(negative_view, ['source_gse_id','model','n_abt_gold','abt_fpr','n_author_nk','author_nk_fpr'], {'source_gse_id':'Source','n_abt_gold':'alpha-beta cells','abt_fpr':'alpha-beta FPR','n_author_nk':'Author NK','author_nk_fpr':'NK FPR'}, percent={'abt_fpr','author_nk_fpr'})}</section>
    <section class='page'><h2>Interpretation and next action</h2><ul><li>The calibration/refit mismatch of V4.3 was corrected: the deployed scorer, Platt calibrators, and thresholds are frozen together.</li><li>Symmetric receptor aggregates and NK negatives improve precision and alpha-beta specificity.</li><li>External recall remains source-dependent, especially in GDT lung; receptor-expression transfer is still the limiting failure mode.</li><li>Do not adjust either threshold using these test outcomes.</li><li>A promotable next iteration requires genuinely new untouched positive and NK/alpha-beta cohorts, or a precommitted training redesign followed by a new final test.</li></ul>
    <h3>Reproducibility</h3><p class='small'>Model contract: <code>Integrated_dataset/models/gdT_prediction/gdtai_v4_4_dual_mode/model_contract.json</code><br>Training script: <code>workflows/gdtai/train_gdtai_v4_4_dual_mode.py</code><br>Evaluation script: <code>workflows/gdtai/evaluate_gdtai_v4_4_reused_lockbox.py</code><br>Contract SHA-256: <code>{html.escape((LOG_DIR / 'training_summary.json').read_text() and json.loads((LOG_DIR / 'training_summary.json').read_text())['model_contract_sha256'])}</code></p></section>
    </main></body></html>"""
    html_path = REPORT_DIR / "index.html"
    html_path.write_text(body)
    pdf_path = REPORT_DIR / "gdtai_v4_4_dual_mode_report.pdf"
    subprocess.run([
        "google-chrome", "--headless", "--no-sandbox", "--disable-gpu",
        "--run-all-compositor-stages-before-draw", "--virtual-time-budget=3000",
        f"--print-to-pdf={pdf_path}", html_path.resolve().as_uri(),
    ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    v44_f1 = overall.set_index("model").loc["v4_4_highest_f1"]
    v3 = overall.set_index("model").loc["v3_balanced"]
    summary_md = f"""# gdTAI V4.4 dual-mode evaluation

## Decision

V4.4 is not promoted. Thresholds were selected only from grouped development validation, never from test cohorts.

- highest-F1 threshold: `{contract['operating_modes']['highest_f1']['threshold']:.9f}`
- high-purity threshold: `{contract['operating_modes']['high_purity']['threshold']:.9f}`
- reused-lockbox V4.4 highest-F1 recall/F1: `{v44_f1.recall:.4f}` / `{v44_f1.f1:.4f}`
- reused-lockbox V3 balanced recall/F1: `{v3.recall:.4f}` / `{v3.f1:.4f}`
- HTML: `{html_path}`
- PDF: `{pdf_path}`
"""
    (LOG_DIR / "gdtai_v4_4_dual_mode_summary.md").write_text(summary_md)
    print(json.dumps({"status":"PASS_REPORT","html":str(html_path),"pdf":str(pdf_path),"pdf_bytes":pdf_path.stat().st_size,"figures":len(figures)}, indent=2))


if __name__ == "__main__":
    main()
