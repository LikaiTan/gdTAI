#!/usr/bin/env python3
"""Render the gdTAI V4.6 development and BALF diagnostic report."""

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
DEV_TABLES = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_6_development"
BALF_TABLES = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_6_balf"
FIGURES = ROOT / "Integrated_dataset/figures/gdT_prediction/gdtai_v4_6_development"
REPORT = ROOT / "gdT_prediction/gdtai_v4_6_development"
MODEL_CONTRACT = ROOT / "Integrated_dataset/models/gdT_prediction/gdtai_v4_6_development/model_contract.json"
LABELS = {
    "v4_6_highest_f1": "V4.6 highest F1",
    "v4_6_high_purity": "V4.6 high purity",
    "v3_balanced_historical": "V3 balanced*",
    "v2_high_f1_historical": "V2 highest F1*",
    "v2_high_purity_historical": "V2 high purity*",
}
COLORS = ["#16817a", "#9b3f52", "#176b87", "#d18b32", "#725a7a"]


def savefig(fig: plt.Figure, name: str) -> Path:
    FIGURES.mkdir(parents=True, exist_ok=True)
    path = FIGURES / name
    fig.savefig(path, dpi=260, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def table_html(frame: pd.DataFrame, columns: list[str], *, percentages: set[str] = set(),
               decimals: set[str] = set(), labels: dict[str, str] | None = None) -> str:
    labels = labels or {}
    out = ["<div class='table-wrap'><table><thead><tr>"]
    out.extend(f"<th>{html.escape(labels.get(column, column))}</th>" for column in columns)
    out.append("</tr></thead><tbody>")
    for _, row in frame.iterrows():
        out.append("<tr>")
        for column in columns:
            value = row.get(column)
            if pd.isna(value):
                rendered = "-"
            elif column in percentages:
                rendered = f"{100 * float(value):.3f}%"
            elif column in decimals:
                rendered = f"{float(value):.4f}"
            elif isinstance(value, (int, np.integer)) or (isinstance(value, float) and value.is_integer()):
                rendered = f"{int(value):,}"
            else:
                rendered = str(value)
            out.append(f"<td>{html.escape(rendered)}</td>")
        out.append("</tr>")
    out.append("</tbody></table></div>")
    return "".join(out)


def architecture_figure(summary: pd.DataFrame) -> Path:
    labels = summary.candidate.str.replace("_", " ").str.title()
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.7))
    specs = [
        ("highest_f1_macro_positive_recall", "Macro positive recall", False),
        ("highest_f1_macro_negative_fpr", "Macro negative FPR", True),
        ("highest_f1_f1_at_1pct", "F1 at 1% prevalence", False),
    ]
    for ax, (column, title, percent) in zip(axes, specs):
        bars = ax.bar(labels, summary[column], color=["#16817a" if value else "#9aaeb5" for value in summary.selected])
        ax.set_title(title, weight="bold")
        ax.tick_params(axis="x", rotation=18)
        ax.spines[["top", "right"]].set_visible(False)
        if percent:
            ax.yaxis.set_major_formatter(PercentFormatter(1))
        else:
            ax.set_ylim(0, 1)
        ax.bar_label(bars, labels=[f"{100*v:.2f}%" if percent else f"{v:.3f}" for v in summary[column]], fontsize=8)
    fig.suptitle("Precommitted source-heldout architecture selection", fontsize=13, weight="bold")
    fig.tight_layout()
    return savefig(fig, "architecture_selection.png")


def holdout_figure(holdout: pd.DataFrame) -> Path:
    selected = holdout[holdout.candidate.eq("conditional_gamma")]
    positive = selected[selected.holdout_role.eq("positive")].pivot(
        index="heldout_source", columns="mode", values="gdt_recall"
    )
    negative = selected[selected.holdout_role.eq("negative")].pivot(
        index="heldout_source", columns="mode", values="combined_negative_fpr"
    ).sort_values("highest_f1", ascending=False)
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), gridspec_kw={"width_ratios": [0.82, 1.45]})
    for ax, frame, title, ylabel in (
        (axes[0], positive, "Whole-positive-source recall", "gdT-gold recall"),
        (axes[1], negative, "Whole-negative-source false positives", "Combined negative FPR"),
    ):
        x = np.arange(len(frame))
        width = 0.36
        for index, mode in enumerate(("highest_f1", "high_purity")):
            ax.bar(x + (index - 0.5) * width, frame[mode], width,
                   color=["#16817a", "#9b3f52"][index], label=mode.replace("_", " ").title())
        ax.set_xticks(x, frame.index, rotation=28, ha="right")
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.set_title(title, weight="bold")
        ax.set_ylabel(ylabel)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    return savefig(fig, "source_holdout_performance.png")


def balf_performance_figure(metrics: pd.DataFrame) -> Path:
    order = list(LABELS)
    view = metrics.set_index("model").loc[order]
    x = np.arange(len(view))
    fig, ax = plt.subplots(figsize=(10.8, 4.5))
    width = 0.23
    for index, metric in enumerate(("precision", "recall", "f1")):
        bars = ax.bar(x + (index - 1) * width, view[metric], width,
                      color=["#5b8e7d", "#d18b32", "#176b87"][index], label=metric.title())
        if metric == "f1":
            ax.bar_label(bars, labels=[f"{value:.3f}" for value in view[metric]], fontsize=8, padding=2)
    ax.set_xticks(x, [LABELS[item] for item in order], rotation=13, ha="right")
    ax.set_ylim(0.82, 1.01)
    ax.set_title("BALF_BLOOD_COPD: frozen operating-point performance", weight="bold")
    ax.legend(frameon=False, ncol=3)
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, "balf_performance.png")


def balf_score_figure(predictions: pd.DataFrame, contract: dict) -> Path:
    fig, ax = plt.subplots(figsize=(9.3, 4.6))
    bins = np.linspace(0, 1, 81)
    for truth, label, color in (
        ("gdT_gold", "gdT gold", "#16817a"),
        ("abT_gold", "alpha-beta gold", "#176b87"),
        ("nk_lockbox", "author NK", "#9b3f52"),
    ):
        values = predictions.loc[predictions.truth_class.eq(truth), "v4_6_stage2_probability"]
        ax.hist(values, bins=bins, density=True, histtype="step", linewidth=2, label=f"{label} (n={len(values):,})", color=color)
    for mode, color in (("highest_f1", "#16817a"), ("high_purity", "#9b3f52")):
        ax.axvline(contract["operating_modes"][mode]["threshold"], color=color, linestyle="--", linewidth=1.5,
                   label=f"{mode.replace('_', ' ')} threshold")
    ax.set_yscale("log")
    ax.set_xlabel("Calibrated V4.6 Stage-2 probability")
    ax.set_ylabel("Density (log scale)")
    ax.set_title("Untouched BALF score distributions", weight="bold")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, "balf_score_distributions.png")


def bootstrap_figure(bootstrap: pd.DataFrame) -> Path:
    display = bootstrap.copy()
    display["label"] = display.comparator.map({"v3_balanced": "V3 balanced*", "v2_high_f1": "V2 highest F1*"})
    y = np.arange(len(display))
    fig, ax = plt.subplots(figsize=(7.6, 3.3))
    ax.errorbar(display["median"], y,
                xerr=[display["median"] - display["ci_low"], display["ci_high"] - display["median"]],
                fmt="o", color="#16817a", capsize=5, linewidth=2)
    ax.axvline(0, color="#52656c", linewidth=1)
    ax.set_yticks(y, display.label)
    ax.set_xlabel("V4.6 highest-F1 minus comparator F1")
    ax.set_title("Truth-stratified donor/group bootstrap", weight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, "balf_bootstrap_f1_difference.png")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-pdf", action="store_true")
    args = parser.parse_args()
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10})
    candidate = pd.read_csv(DEV_TABLES / "candidate_selection_summary.csv")
    holdout = pd.read_csv(DEV_TABLES / "candidate_source_holdout_metrics.csv")
    validation = pd.read_csv(DEV_TABLES / "threshold_validation_metrics.csv")
    architecture = pd.read_csv(DEV_TABLES / "feature_architecture_preflight.csv")
    metrics = pd.read_csv(BALF_TABLES / "overall_metrics.csv")
    predictions = pd.read_parquet(BALF_TABLES / "balf_predictions.parquet")
    bootstrap = pd.read_csv(BALF_TABLES / "balf_group_bootstrap_summary.csv")
    mcnemar = pd.read_csv(BALF_TABLES / "balf_paired_error_tests.csv")
    integrity = pd.read_csv(BALF_TABLES / "implementation_integrity.csv")
    contract = json.loads(MODEL_CONTRACT.read_text())

    figures = [
        architecture_figure(candidate), holdout_figure(holdout), balf_performance_figure(metrics),
        balf_score_figure(predictions, contract), bootstrap_figure(bootstrap),
    ]
    REPORT.mkdir(parents=True, exist_ok=True)
    assets = REPORT / "assets"
    assets.mkdir(parents=True, exist_ok=True)
    for figure in figures:
        shutil.copy2(figure, assets / figure.name)

    metric_view = metrics.copy()
    metric_view["model"] = metric_view.model.map(LABELS)
    selected_positive = holdout[
        holdout.candidate.eq("conditional_gamma") & holdout.holdout_role.eq("positive")
    ][["heldout_source", "mode", "gdt_recall"]]
    selected_negative = holdout[
        holdout.candidate.eq("conditional_gamma") & holdout.holdout_role.eq("negative")
    ][["heldout_source", "mode", "n_negative", "combined_negative_fpr", "abt_fpr", "tier1_nk_fpr"]]
    selected_negative = selected_negative.sort_values(["mode", "combined_negative_fpr"], ascending=[True, False])

    css = """
    :root{--ink:#17272d;--muted:#52656c;--line:#cad6da;--teal:#16817a;--red:#9b3f52;--blue:#176b87}
    *{box-sizing:border-box}body{margin:0;background:#eef3f4;color:var(--ink);font-family:Arial,Helvetica,sans-serif;line-height:1.42}
    header{background:#17272d;color:white;padding:32px 42px}header h1{font-size:28px;margin:0 0 7px;letter-spacing:0}header p{margin:0;color:#d8e5e8}
    main{max-width:1180px;margin:auto;background:white}section{padding:25px 38px;border-bottom:1px solid var(--line)}
    h2{font-size:20px;color:#174453;margin:0 0 12px;letter-spacing:0}h3{font-size:15px;margin:16px 0 7px;letter-spacing:0}p{margin:7px 0}
    .callout{padding:12px 15px;border-left:5px solid var(--teal);background:#f0faf8}.caution{border-left-color:#d18b32;background:#fff8ec}.stop{border-left-color:var(--red);background:#fff4f6}
    .grid{display:grid;grid-template-columns:1fr 1fr;gap:20px;align-items:start}.kpis{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin:13px 0}
    .kpi{border-top:4px solid var(--teal);background:#f4f8f9;padding:10px}.kpi b{display:block;font-size:19px;color:#174453}.kpi span{font-size:9pt;color:var(--muted)}
    figure{margin:10px 0}figure img{width:100%;height:auto;display:block}figcaption,.small{font-size:9pt;color:var(--muted)}
    .table-wrap{overflow-x:auto}table{width:100%;border-collapse:collapse;font-size:8.4pt;margin:9px 0 14px;table-layout:fixed}th{background:#e7f0f2;color:#174453;text-align:left}
    th,td{border:1px solid var(--line);padding:5px 6px;vertical-align:top;overflow-wrap:anywhere}tr:nth-child(even) td{background:#f8fafb}code{font-family:Consolas,monospace;font-size:.9em}
    ul{margin:7px 0;padding-left:19px}@media(max-width:780px){.grid,.kpis{grid-template-columns:1fr}section{padding:20px}header{padding:25px 20px}}
    @media print{@page{size:A4 landscape;margin:8mm}body{background:white;font-size:8.8pt}main{max-width:none}header{padding:16px 22px}header h1{font-size:20px}
      section{padding:12px 18px}h2{font-size:15px}h3{font-size:11px}.kpi b{font-size:15px}figure{break-inside:avoid}figure img{max-height:116mm;object-fit:contain}
      table{font-size:7.1pt;table-layout:auto}tr{break-inside:avoid}.page-break{break-before:page}.table-wrap{overflow:visible}
      .kpis{display:flex}.kpi{flex:1}.compact img{max-height:72mm}}
    """
    body = f"""<!doctype html><html><head><meta charset='utf-8'><title>gdTAI V4.6 development report</title><style>{css}</style></head><body>
    <header><h1>gdTAI V4.6: negative-diversity and receptor-ablation iteration</h1><p>Frozen source-heldout development, dual operating modes, and one untouched BALF diagnostic | 20 August 2026</p></header><main>
    <section><h2>Decision</h2><div class='callout'><b>V4.6 highest-F1 is the first V4 iteration to improve F1 over V2/V3 on BALF while sharply reducing the V4.5 GSE159251 failure.</b> The gain is donor-bootstrap supported, but V4.6 is not promoted because BALF is the only new positive test and is now consumed.</div>
    <div class='kpis'><div class='kpi'><b>0.9544</b><span>BALF V4.6 highest-F1 F1</span></div><div class='kpi'><b>0.029%</b><span>BALF alpha-beta FPR</span></div><div class='kpi'><b>0 / 254</b><span>BALF author-NK false positives</span></div><div class='kpi'><b>0.92%</b><span>GSE159251 source-heldout FPR</span></div></div>
    <div class='callout caution'><b>Operational status:</b> V3 balanced remains the default and V2 high-purity remains the conservative fallback. V4.6 is a frozen development candidate awaiting a genuinely new confirmatory cohort.</div></section>

    <section><h2>Why this iteration</h2><p>V4.5 produced 2,240 highest-F1 alpha-beta false positives in GSE159251. Only 13.35% expressed a dedicated TRDV gene, while broad TRGV and alpha/beta V aggregates drove many errors. V4.6 therefore expanded source-diverse negatives and tested whether gamma evidence should be removed or conditioned on delta support.</p>
    <ul><li>No hard TRDV cutoff: dropout-tolerant gdT cells remain recoverable.</li><li>No single NK/cytotoxic gene veto: cytotoxic gdT biology is not penalized directly.</li><li>Silver positives are excluded. Consumed cohorts are development evidence, never relabeled as external tests.</li><li>BALF did not fit a model/calibrator or select architecture/threshold.</li></ul></section>

    <section class='page-break'><h2>Architecture selection</h2><figure><img src='assets/architecture_selection.png'><figcaption>Candidate choice uses five whole-positive-source and ten whole-negative-source exclusions. F1 is reconstructed from macro recall and macro FPR at 1% gdT prevalence.</figcaption></figure>
    {table_html(candidate, ['candidate','highest_f1_macro_positive_recall','highest_f1_macro_negative_fpr','highest_f1_precision_at_1pct','highest_f1_f1_at_1pct','highest_f1_minimum_positive_recall','highest_f1_maximum_negative_fpr','selected'], percentages={'highest_f1_macro_positive_recall','highest_f1_macro_negative_fpr','highest_f1_precision_at_1pct','highest_f1_minimum_positive_recall','highest_f1_maximum_negative_fpr'}, decimals={'highest_f1_f1_at_1pct'})}
    <p><b>Selected:</b> conditional gamma. It retains gamma support only when TRDC or dedicated TRDV support is detected. Its 182 effective features include zero standalone TRG genes and zero unconditioned TRGV/TRGC aggregates.</p>
    {table_html(architecture, ['candidate','architecture','n_effective_features','n_standalone_trg_features','n_forbidden_gamma_features','feature_names_sha256'])}</section>

    <section class='page-break'><h2>Cross-source development performance</h2><figure><img src='assets/source_holdout_performance.png'><figcaption>Every point comes from a model trained without the named source. GDTlung remains the dominant positive domain gap.</figcaption></figure>
    <h3>Positive-source recall</h3>{table_html(selected_positive, ['heldout_source','mode','gdt_recall'], percentages={'gdt_recall'})}
    <h3>Negative-source false-positive rates</h3>{table_html(selected_negative, ['heldout_source','mode','n_negative','combined_negative_fpr','abt_fpr','tier1_nk_fpr'], percentages={'combined_negative_fpr','abt_fpr','tier1_nk_fpr'})}
    <div class='callout caution'><b>Residual weakness:</b> GDTlung highest-F1 recall is 36.61%, consistent with its known suboptimal library quality. GSE228597 is the worst negative holdout at 1.05% FPR; GSE159251 is reduced to 0.92% but not fully below historical V3/V2 controls.</div></section>

    <section class='page-break'><h2>Frozen dual operating modes</h2><p>Both thresholds were selected once from the group-disjoint development threshold-validation partition after architecture selection.</p>
    {table_html(validation, ['mode','threshold','tp','fp','tn','fn','precision','recall','specificity','f1','f0_5'], percentages={'precision','recall','specificity'}, decimals={'threshold','f1','f0_5'})}
    <p><b>Highest F1:</b> threshold 0.894354. <b>High purity:</b> threshold 0.971830. Scores between them form the uncertainty band; no test result altered either threshold.</p></section>

    <section class='page-break'><h2>Untouched BALF diagnostic</h2><div class='grid'><figure><img src='assets/balf_performance.png'><figcaption>*V2/V3 are historical comparators whose prior project use means they are not independent tests here.</figcaption></figure><figure><img src='assets/balf_score_distributions.png'><figcaption>Distributions use the frozen calibrated scorer. The Stage-1 gate and CD4/Treg exclusions are applied to final calls.</figcaption></figure></div>
    {table_html(metric_view, ['model','positive','negative','tp','fp','fn','precision','recall','f1','abt_fpr','author_nk_fpr','roc_auc','pr_auc','comparison_status'], percentages={'precision','recall','abt_fpr','author_nk_fpr'}, decimals={'f1','roc_auc','pr_auc'})}
    <p>V4.6 highest-F1 recovered 785/852 gdT-gold cells with 8 alpha-beta false positives and no author-NK false positives. Its ranking AUC is lower than V2/V3, so the advantage is at the frozen operating point rather than across every possible threshold.</p></section>

    <section class='page-break'><h2>Paired uncertainty</h2><div class='grid'><figure class='compact'><img src='assets/balf_bootstrap_f1_difference.png'><figcaption>5,000 truth-stratified donor/group resamples; intervals are paired V4.6 minus comparator F1.</figcaption></figure><div>{table_html(bootstrap, ['comparator','iterations','median','ci_low','ci_high','probability_gt_zero','bootstrap_unit'], decimals={'median','ci_low','ci_high','probability_gt_zero'})}</div></div>
    <h3>Descriptive paired cell-error tests</h3>{table_html(mcnemar, ['comparator','v4_6_correct_comparator_wrong','v4_6_wrong_comparator_correct','discordant','exact_mcnemar_p','interpretation'], decimals={'exact_mcnemar_p'})}
    <p>The bootstrap is the primary uncertainty analysis. McNemar tests are descriptive because cell-level independence is false within donors.</p></section>

    <section class='page-break'><h2>Interpretation and next gate</h2><div class='callout stop'><b>Do not promote or push V4.6 yet.</b> One positive external cohort is insufficient to rule out another source-specific failure, and BALF cannot be reused for future tuning after this report.</div>
    <h3>What improved</h3><ul><li>V4.5 GSE159251 highest-F1 FPR fell from 6.91% to 0.92% under honest whole-source exclusion.</li><li>BALF highest-F1 improved both recall and precision versus historical V3, with a positive donor-bootstrap F1 interval.</li><li>Author-NK FPR was zero in BALF without hard exclusion of NKG7/GNLY/KLRD1-like cytotoxic programs.</li></ul>
    <h3>What remains unresolved</h3><ul><li>GDTlung recall remains poor and cannot be solved by threshold transfer alone.</li><li>GSE159251 and GSE228597 still show source-specific alpha-beta false positives.</li><li>V4.6 high-purity F1/recall does not beat V2 high-purity on BALF.</li><li>V4.6 has lower BALF ROC-AUC/PR-AUC than V2/V3 despite a better frozen highest-F1 operating point.</li></ul>
    <h3>Required next evidence</h3><p>Freeze V4.6 unchanged and evaluate it once on a genuinely new cohort containing receptor-confirmed gdT positives, alpha-beta T cells, and independently annotated NK cells. Promotion should require dataset-level superiority or non-inferiority in both modes, not another adaptive search on consumed cohorts.</p></section>

    <section><h2>Reproducibility and integrity</h2>{table_html(integrity, ['check','pass','detail'])}<p class='small'>Model contract: <code>{html.escape(str(MODEL_CONTRACT.relative_to(ROOT)))}</code><br>Contract SHA-256: <code>{html.escape(json.loads((ROOT/'Integrated_dataset/logs/gdT_prediction/gdtai_v4_6_balf/evaluation_summary.json').read_text())['v4_6_model_contract_sha256_before_scoring'])}</code><br>Normalization: log1p raw counts normalized to 10,000 counts per cell. No H5AD was mutated. No V4 artifact was promoted, atlas-applied, pushed, or served on port 8000.</p></section>
    </main></body></html>"""
    html_path = REPORT / "index.html"
    html_path.write_text(body)
    pdf_path = REPORT / "gdtai_v4_6_development_report.pdf"
    if not args.skip_pdf:
        subprocess.run([
            "google-chrome", "--headless", "--no-sandbox", "--disable-gpu",
            "--disable-dev-shm-usage", "--user-data-dir=/tmp/gdtai-v4-6-chrome-profile",
            "--run-all-compositor-stages-before-draw", f"--print-to-pdf={pdf_path}",
            html_path.resolve().as_uri(),
        ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elif not pdf_path.exists():
        raise FileNotFoundError(f"--skip-pdf requested but PDF is absent: {pdf_path}")
    print(json.dumps({
        "status": "PASS_V4_6_REPORT_RENDERED", "html": str(html_path), "pdf": str(pdf_path),
        "pdf_bytes": pdf_path.stat().st_size, "figures": len(figures),
    }, indent=2))


if __name__ == "__main__":
    main()
