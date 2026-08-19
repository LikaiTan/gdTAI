#!/usr/bin/env python3
"""Render the gdTAI V4.5 development-only ablation report."""

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
TABLES = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_5_development"
MODEL_DIR = ROOT / "Integrated_dataset/models/gdT_prediction/gdtai_v4_5_development"
FIGURES = ROOT / "Integrated_dataset/figures/gdT_prediction/gdtai_v4_5_development"
LOGS = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_5_development"
REPORT = ROOT / "gdT_prediction/gdtai_v4_5_development"
COLORS = {"baseline": "#176b87", "cytotoxic": "#c05a47", "f1": "#16817a", "purity": "#9b3f52"}


def savefig(fig: plt.Figure, name: str) -> Path:
    FIGURES.mkdir(parents=True, exist_ok=True)
    path = FIGURES / name
    fig.savefig(path, dpi=260, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def table_html(
    frame: pd.DataFrame,
    columns: list[str],
    labels: dict[str, str] | None = None,
    percentages: set[str] | None = None,
    decimals: set[str] | None = None,
) -> str:
    labels = labels or {}
    percentages = percentages or set()
    decimals = decimals or set()
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
                rendered = f"{100 * float(value):.2f}%"
            elif column in decimals:
                rendered = f"{float(value):.4f}"
            elif isinstance(value, (float, np.floating)):
                rendered = f"{value:.3f}"
            elif isinstance(value, (int, np.integer)):
                rendered = f"{int(value):,}"
            else:
                rendered = str(value)
            output.append(f"<td>{html.escape(rendered)}</td>")
        output.append("</tr>")
    output.append("</tbody></table>")
    return "".join(output)


def plot_partitions(counts: pd.DataFrame) -> Path:
    totals = counts.groupby("v4_5_partition", observed=True).n_cells.sum().reindex(
        ["fit", "calibration", "threshold_validation", "locked_test"]
    )
    fig, ax = plt.subplots(figsize=(8.4, 3.8))
    bars = ax.bar(totals.index, totals, color=["#176b87", "#5b8e7d", "#16817a", "#7a8589"])
    ax.bar_label(bars, labels=[f"{int(value):,}" for value in totals], padding=4, fontsize=9)
    ax.set_ylabel("Labeled/audit cells")
    ax.set_title("Whole-group V4.5 partitions", weight="bold")
    ax.tick_params(axis="x", rotation=10)
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, "development_partitions.png")


def plot_candidate_summary(candidates: pd.DataFrame) -> Path:
    view = candidates.set_index("candidate").loc[
        ["receptor_context_baseline", "restored_cytotoxic_context"]
    ]
    labels = ["Receptor context", "Restored cytotoxic context"]
    metrics = ["macro_source_holdout_f1", "macro_source_holdout_recall"]
    x = np.arange(2)
    fig, ax = plt.subplots(figsize=(7.6, 4.3))
    width = 0.34
    for index, metric in enumerate(metrics):
        values = view[metric].to_numpy()
        bars = ax.bar(x + (index - 0.5) * width, values, width,
                      color=["#176b87", "#c05a47"][index],
                      label=["Macro F1", "Macro recall"][index])
        ax.bar_label(bars, labels=[f"{value:.3f}" for value in values], fontsize=8, padding=3)
    ax.set_xticks(x, labels)
    ax.set_ylim(0.78, 0.88)
    ax.set_title("Precommitted whole-source candidate selection", weight="bold")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, "candidate_macro_comparison.png")


def plot_holdouts(heldout: pd.DataFrame) -> Path:
    view = heldout.loc[heldout["mode"].eq("highest_f1")]
    pivot = view.pivot(index="heldout_source", columns="candidate", values="gdt_recall").reindex([
        "HRA005041", "GSE144469", "MalteGDT", "GDT_2020AUG_woCOV", "GDTlung2023july_7p"
    ])
    labels = ["HRA005041", "GSE144469", "Malte", "GDT2020", "GDTlung"]
    x = np.arange(len(pivot))
    fig, ax = plt.subplots(figsize=(9.4, 4.6))
    width = 0.36
    for index, candidate in enumerate(["receptor_context_baseline", "restored_cytotoxic_context"]):
        values = pivot[candidate].to_numpy()
        bars = ax.bar(x + (index - 0.5) * width, values, width,
                      color=[COLORS["baseline"], COLORS["cytotoxic"]][index],
                      label=["Receptor context", "Restored cytotoxic context"][index])
        ax.bar_label(bars, labels=[f"{100*value:.1f}%" for value in values], fontsize=7.5, padding=2)
    ax.set_xticks(x, labels)
    ax.set_ylim(0, 1.08)
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.set_ylabel("gdT-gold recall")
    ax.set_title("Five whole-source development holdouts", weight="bold")
    ax.legend(frameon=False, loc="lower left")
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, "candidate_source_holdout_recall.png")


def plot_frontier(frontier: pd.DataFrame, contract: dict) -> Path:
    step = max(1, len(frontier) // 5000)
    view = frontier.iloc[::step]
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    ax.plot(view.recall, view.precision, color="#5c6970", lw=1.4)
    for mode, color, label in [
        ("highest_f1", COLORS["f1"], "Highest F1"),
        ("high_purity", COLORS["purity"], "High purity"),
    ]:
        payload = contract["operating_modes"][mode]
        ax.scatter(payload["recall"], payload["precision"], s=72, color=color,
                   edgecolor="white", linewidth=0.8, label=label, zorder=3)
    ax.set_xlim(0.70, 1.005)
    ax.set_ylim(0.75, 1.005)
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Disjoint threshold-validation frontier", weight="bold")
    ax.grid(alpha=0.16)
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, "threshold_validation_frontier.png")


def plot_modes(metrics: pd.DataFrame) -> Path:
    view = metrics.set_index("mode").loc[["highest_f1", "high_purity"]]
    values = view[["precision", "recall", "f1", "f0_5"]]
    x = np.arange(4)
    fig, ax = plt.subplots(figsize=(7.8, 4.3))
    width = 0.34
    for index, (mode, color, label) in enumerate([
        ("highest_f1", COLORS["f1"], "Highest F1"),
        ("high_purity", COLORS["purity"], "High purity"),
    ]):
        row = values.loc[mode].to_numpy()
        bars = ax.bar(x + (index - 0.5) * width, row, width, color=color, label=label)
        ax.bar_label(bars, labels=[f"{value:.3f}" for value in row], fontsize=8, padding=2)
    ax.set_xticks(x, ["Precision", "Recall", "F1", "F0.5"])
    ax.set_ylim(0.75, 1.02)
    ax.set_title("Frozen development operating modes", weight="bold")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, "operating_mode_metrics.png")


def plot_uncertainty(bands: pd.DataFrame) -> Path:
    classes = ["gdT_gold", "abT_gold", "nk_tier1", "nk_tier2"]
    band_order = ["negative", "uncertain_gdt", "high_confidence_gdt"]
    pivot = bands[bands.truth_class.isin(classes)].pivot(
        index="truth_class", columns="decision_band", values="fraction_within_truth"
    ).fillna(0).reindex(classes)[band_order]
    fig, ax = plt.subplots(figsize=(8.3, 4.4))
    left = np.zeros(len(pivot))
    for band, color, label in zip(
        band_order, ["#cdd6d9", "#e0a458", COLORS["f1"]],
        ["Negative", "Uncertain gdT", "High-confidence gdT"],
    ):
        values = pivot[band].to_numpy()
        ax.barh(pivot.index, values, left=left, color=color, label=label)
        left += values
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(PercentFormatter(1))
    ax.set_xlabel("Fraction within truth class")
    ax.set_title("Decision-band composition", weight="bold")
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.32))
    ax.spines[["top", "right"]].set_visible(False)
    return savefig(fig, "uncertainty_band_composition.png")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-pdf", action="store_true")
    args = parser.parse_args()
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 12, "axes.labelsize": 10})
    contract = json.loads((MODEL_DIR / "model_contract.json").read_text())
    training = json.loads((LOGS / "training_summary.json").read_text())
    counts = pd.read_csv(TABLES / "development_partition_counts.csv")
    candidates = pd.read_csv(TABLES / "candidate_selection_summary.csv")
    heldout = pd.read_csv(TABLES / "candidate_source_holdout_metrics.csv")
    frontier = pd.read_csv(TABLES / "threshold_validation_frontier.csv")
    metrics = pd.read_csv(TABLES / "threshold_validation_metrics.csv")
    bands = pd.read_csv(TABLES / "uncertainty_band_summary.csv")
    figures = {
        "partitions": plot_partitions(counts),
        "candidates": plot_candidate_summary(candidates),
        "holdouts": plot_holdouts(heldout),
        "frontier": plot_frontier(frontier, contract),
        "modes": plot_modes(metrics),
        "bands": plot_uncertainty(bands),
    }

    REPORT.mkdir(parents=True, exist_ok=True)
    assets = REPORT / "assets"
    assets.mkdir(parents=True, exist_ok=True)
    for source in figures.values():
        shutil.copy2(source, assets / source.name)

    selected = candidates.loc[candidates.selected].iloc[0]
    holdout_view = heldout.loc[heldout["mode"].eq("highest_f1"), [
        "candidate", "heldout_source", "gdt_recall", "precision", "f1", "abt_fpr"
    ]].copy()
    holdout_view["candidate"] = holdout_view.candidate.map({
        "receptor_context_baseline": "Receptor context",
        "restored_cytotoxic_context": "Restored cytotoxic context",
    })
    candidate_view = candidates[[
        "candidate", "macro_source_holdout_f1", "macro_source_holdout_recall",
        "minimum_source_holdout_recall", "macro_source_holdout_abt_fpr", "selected"
    ]].copy()
    candidate_view["candidate"] = candidate_view.candidate.map({
        "receptor_context_baseline": "Receptor context",
        "restored_cytotoxic_context": "Restored cytotoxic context",
    })
    metrics_view = metrics[["mode", "threshold", "precision", "recall", "f1", "f0_5", "specificity"]]
    band_view = bands[bands.truth_class.isin(["gdT_gold", "abT_gold", "nk_tier1", "nk_tier2"])][[
        "decision_band", "truth_class", "n_cells", "fraction_within_truth"
    ]]
    positive_counts = counts[counts.truth_class.eq("gdT_gold")][[
        "v4_5_partition", "source_gse_id", "n_cells"
    ]]

    css = """
    :root{--ink:#17272d;--muted:#52656c;--line:#ccd8dc;--teal:#16817a;--red:#9b3f52;--blue:#176b87}
    *{box-sizing:border-box}body{margin:0;background:#f3f6f7;color:var(--ink);font-family:Arial,Helvetica,sans-serif;line-height:1.45}
    header{background:#17272d;color:white;padding:34px 42px 28px}header h1{font-size:29px;margin:0 0 7px;letter-spacing:0}header p{margin:0;color:#dbe6e8}
    main{max-width:1160px;margin:auto;background:white}section{padding:27px 40px;border-bottom:1px solid var(--line)}
    h2{font-size:20px;margin:0 0 13px;color:#174453;letter-spacing:0}h3{font-size:15px;margin:17px 0 7px;letter-spacing:0}p{margin:8px 0}
    .verdict{border-left:5px solid var(--red);background:#fff6f7;padding:12px 15px}.ok{border-left-color:var(--teal);background:#f2fbfa}
    .grid{display:grid;grid-template-columns:1fr 1fr;gap:22px;align-items:start}figure{margin:12px 0 5px}figure img{width:100%;height:auto;display:block}
    figcaption,.small{font-size:9.3pt;color:var(--muted)}table{width:100%;border-collapse:collapse;font-size:8.8pt;margin:10px 0 15px;table-layout:fixed}
    th{background:#e8f0f2;color:#174453;text-align:left}th,td{border:1px solid var(--line);padding:5px 6px;vertical-align:top;overflow-wrap:anywhere}tr:nth-child(even) td{background:#f8fafb}
    code{font-family:Consolas,monospace;font-size:.9em}ul{padding-left:19px;margin:8px 0}
    @media(max-width:780px){.grid{grid-template-columns:1fr}section{padding:21px}header{padding:27px 21px}}
    @media print{@page{size:A4 landscape;margin:9mm}body{background:white;font-size:9pt}header{padding:17px 23px}header h1{font-size:21px}main{max-width:none}
      section{padding:13px 20px}h2{font-size:15px}h3{font-size:11.5px}figure{break-inside:avoid}figure img{max-height:121mm;object-fit:contain}
      table{font-size:7pt;table-layout:fixed}th,td{padding:3.5px}tr{break-inside:avoid}.page{break-before:page}.grid{gap:12px}}
    """
    body = f"""<!doctype html><html><head><meta charset='utf-8'><title>gdTAI V4.5 development report</title><style>{css}</style></head><body>
    <header><h1>gdTAI V4.5 development report</h1><p>Positive-domain expansion, cytotoxic-feature ablation, dual operating modes | 19 August 2026</p></header><main>
    <section><h2>Executive decision</h2><div class='verdict'><strong>Development succeeded; promotion is blocked.</strong> Receptor context outperformed the restored cytotoxic panel under the frozen five-source criterion, but no genuinely new untouched test cohort exists. V3 balanced remains the deployed default.</div>
    <p>The selected candidate reached macro source-heldout F1 <strong>{selected.macro_source_holdout_f1:.3f}</strong> and macro recall <strong>{selected.macro_source_holdout_recall:.3f}</strong>. GDTlung recall remained only <strong>40.5%</strong>, so broader positive training reduced but did not eliminate source dependence.</p></section>
    <section><h2>What V4.5 tested</h2><div class='grid'><div><h3>Scientific changes</h3><ul><li>Added sorted GDT2020 and GDTlung positives to development with source weights 0.60 and 0.35.</li><li>Kept silver positives out of fitting, calibration, candidate selection, and thresholds.</li><li>Compared 200-feature receptor context against a bounded cytotoxic-context alternative.</li><li>Used curated tier-1 and tier-2 NK negatives; no single NK gene is a hard exclusion.</li><li>Retained the frozen CD4/Treg rules and no hard TRDV-expression cutoff.</li></ul></div>
    <div><h3>Critical validity boundary</h3><p>GDT2020 and GDTlung came from the already-consumed common lockbox. They are now development data and can never again count as external validation. BALF and extension alpha-beta cohorts were not used by V4.5, but their previous inspection also means the project still needs a new final test.</p></div></div></section>
    <section class='page'><h2>Ground truth and partitions</h2><div class='grid'><figure><img src='assets/{figures['partitions'].name}'><figcaption>Every biological group stays in exactly one partition. Locked rows do not fit models or calibrators and cannot select features or thresholds.</figcaption></figure><div>{table_html(positive_counts, ['v4_5_partition','source_gse_id','n_cells'], {'v4_5_partition':'Partition','source_gse_id':'Positive source','n_cells':'gdT gold'})}</div></div>
    <div class='verdict ok'><strong>Preflight passed.</strong> Fit/calibration/threshold-validation contain 47,460/3,783/13,639 gdT-gold cells. No silver cell entered development.</div></section>
    <section class='page'><h2>Feature ablation</h2><div class='grid'><figure><img src='assets/{figures['candidates'].name}'><figcaption>Candidate selection used macro F1 across five independently refitted whole-source holdouts.</figcaption></figure><figure><img src='assets/{figures['holdouts'].name}'><figcaption>Restored shared cytotoxic genes did not improve aggregate transfer and reduced the weakest-source recall.</figcaption></figure></div>
    {table_html(candidate_view, list(candidate_view.columns), {'candidate':'Candidate','macro_source_holdout_f1':'Macro F1','macro_source_holdout_recall':'Macro recall','minimum_source_holdout_recall':'Minimum recall','macro_source_holdout_abt_fpr':'Macro alpha-beta FPR','selected':'Selected'}, percentages={'macro_source_holdout_recall','minimum_source_holdout_recall','macro_source_holdout_abt_fpr'}, decimals={'macro_source_holdout_f1'})}</section>
    <section class='page'><h2>Whole-source failure modes</h2>{table_html(holdout_view, list(holdout_view.columns), {'candidate':'Candidate','heldout_source':'Held-out source','gdt_recall':'gdT recall','abt_fpr':'alpha-beta FPR'}, percentages={'gdt_recall','precision','abt_fpr'}, decimals={'f1'})}
    <p>HRA precision is lower because its held-out set contains 506,985 alpha-beta gold cells, whereas the three sorted-only held-out sources contain no negatives and therefore cannot estimate precision or FPR. The sorted-source F1 values should be read primarily as recall summaries.</p></section>
    <section class='page'><h2>Frozen operating modes</h2><div class='grid'><figure><img src='assets/{figures['frontier'].name}'><figcaption>Both operating points come only from the disjoint development threshold-validation partition.</figcaption></figure><figure><img src='assets/{figures['modes'].name}'><figcaption>Highest F1 prioritizes recovery; high purity maximizes F0.5 without a hard FPR constraint.</figcaption></figure></div>
    {table_html(metrics_view, list(metrics_view.columns), {'mode':'Mode','f0_5':'F0.5'}, percentages={'precision','recall','specificity'}, decimals={'threshold','f1','f0_5'})}</section>
    <section class='page'><h2>Uncertainty band</h2><div class='grid'><figure><img src='assets/{figures['bands'].name}'><figcaption>Scores between {contract['operating_modes']['highest_f1']['threshold']:.6f} and {contract['operating_modes']['high_purity']['threshold']:.6f} are labeled uncertain, not automatically correct.</figcaption></figure><div>{table_html(band_view, list(band_view.columns), {'decision_band':'Band','truth_class':'Truth','n_cells':'Cells','fraction_within_truth':'Within truth'}, percentages={'fraction_within_truth'})}</div></div>
    <p>The uncertainty interval contains 1,031 gdT-gold and 560 alpha-beta-gold validation cells. It supports conservative review workflows but does not repair low-scoring GDTlung positives below the highest-F1 threshold.</p></section>
    <section class='page'><h2>Interpretation</h2><ul><li>The V4.3 threshold-transfer defect is addressed because each held-out model and the final deployment scorer are calibrated directly, and final thresholds are derived from that exact scorer.</li><li>Positive-domain expansion markedly improves held-out HRA, GSE144469, Malte, and GDT2020 recall.</li><li>The cytotoxic panel is not promoted: NK labels were sufficient to keep its FPR similar, but the added genes provided no reproducible macro-F1 gain.</li><li>GDTlung remains the dominant domain failure, consistent with its known suboptimal library quality and weak receptor-expression transfer.</li><li>V4.5 development-validation F1 is not directly comparable to V4.4's easier threshold partition because V4.5 adds difficult consumed sorted cohorts.</li><li>Do not score the consumed common lockbox again for model selection. Promotion requires a new untouched positive plus alpha-beta/NK test.</li></ul></section>
    <section class='page'><h2>Reproducibility and status</h2><div class='verdict'><strong>Status: development frozen, not promotable.</strong> No H5AD was modified, no whole-atlas inference was run, and no V4 artifact was pushed to GitHub.</div>
    <p class='small'>Training script: <code>workflows/gdtai/train_gdtai_v4_5_development.py</code><br>Configuration: <code>configs/models/gdtai/v4_5_positive_diversity_ablation.json</code><br>Model contract: <code>Integrated_dataset/models/gdT_prediction/gdtai_v4_5_development/model_contract.json</code><br>Contract SHA-256: <code>{training['model_contract_sha256']}</code><br>Runtime: {training['runtime_seconds']:.1f} seconds after exact vectorization of the unchanged Stage-1 threshold grid.</p></section>
    </main></body></html>"""
    html_path = REPORT / "index.html"
    html_path.write_text(body)
    pdf_path = REPORT / "gdtai_v4_5_development_report.pdf"
    if not args.skip_pdf:
        subprocess.run([
            "google-chrome", "--headless", "--no-sandbox", "--disable-gpu",
            "--run-all-compositor-stages-before-draw", "--virtual-time-budget=3000",
            f"--print-to-pdf={pdf_path}", html_path.resolve().as_uri(),
        ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elif not pdf_path.exists():
        raise FileNotFoundError(f"--skip-pdf requested but PDF is absent: {pdf_path}")

    summary_md = f"""# gdTAI V4.5 development report

## Decision

V4.5 development passed, but the model is not promotable because there is no new untouched final test.

- selected candidate: `{contract['selected_candidate']}`
- macro source-heldout F1: `{selected.macro_source_holdout_f1:.4f}`
- macro source-heldout recall: `{selected.macro_source_holdout_recall:.4f}`
- highest-F1 threshold: `{contract['operating_modes']['highest_f1']['threshold']:.9f}`
- high-purity threshold: `{contract['operating_modes']['high_purity']['threshold']:.9f}`
- model contract SHA-256: `{training['model_contract_sha256']}`
- HTML: `{html_path}`
- PDF: `{pdf_path}`
"""
    (LOGS / "gdtai_v4_5_development_summary.md").write_text(summary_md)
    print(json.dumps({
        "status": "PASS_REPORT", "html": str(html_path), "pdf": str(pdf_path),
        "pdf_bytes": pdf_path.stat().st_size, "figures": len(figures),
    }, indent=2))


if __name__ == "__main__":
    main()
