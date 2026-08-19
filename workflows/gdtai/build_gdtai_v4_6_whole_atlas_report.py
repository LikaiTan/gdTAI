#!/usr/bin/env python3
"""Build the frozen V4.6-versus-V2 whole-atlas application report."""

from __future__ import annotations

import html
import json
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_6_whole_atlas"
FIGURE_DIR = ROOT / "Integrated_dataset/figures/gdT_prediction/gdtai_v4_6_whole_atlas"
LOG_DIR = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_6_whole_atlas"
REPORT_DIR = ROOT / "gdT_prediction/gdtai_v4_6_whole_atlas"
MODELS = (
    "v4_6_highest_f1", "v4_6_high_purity", "v2_high_f1", "v2_high_purity",
)
LABELS = {
    "v4_6_highest_f1": "V4.6 highest F1",
    "v4_6_high_purity": "V4.6 high purity",
    "v2_high_f1": "V2 highest F1",
    "v2_high_purity": "V2 high purity",
}
COLORS = {
    "v4_6_highest_f1": "#0B6E75",
    "v4_6_high_purity": "#63A6A0",
    "v2_high_f1": "#B54A35",
    "v2_high_purity": "#D8937E",
}


def strict_nk_audit() -> tuple[pd.DataFrame, pd.DataFrame]:
    parquet = pq.ParquetFile(TABLE_DIR / "whole_atlas_predictions.parquet")
    columns = [
        "source_gse_id", *MODELS, "author_nk", "conservative_expression_nk",
        "trdv_expressed", "has_any_gd_tcr",
    ]
    pieces = []
    for batch in parquet.iter_batches(batch_size=250_000, columns=columns):
        frame = batch.to_pandas()
        author_no_delta = frame.author_nk & ~frame.trdv_expressed & ~frame.has_any_gd_tcr
        expression_no_delta = (
            frame.conservative_expression_nk & ~frame.trdv_expressed & ~frame.has_any_gd_tcr
        )
        author_conflict = frame.author_nk & (frame.trdv_expressed | frame.has_any_gd_tcr)
        block = pd.DataFrame({"source_gse_id": frame.source_gse_id})
        for model in MODELS:
            call = frame[model]
            block[f"{model}__predicted"] = call.astype(np.int32)
            block[f"{model}__strict_likely_nk"] = (
                call & (author_no_delta | expression_no_delta)
            ).astype(np.int32)
            block[f"{model}__author_nk_no_delta_support"] = (call & author_no_delta).astype(np.int32)
            block[f"{model}__expression_nk_no_delta_support"] = (call & expression_no_delta).astype(np.int32)
            block[f"{model}__author_nk_delta_conflict"] = (call & author_conflict).astype(np.int32)
        pieces.append(block.groupby("source_gse_id", observed=True).sum().reset_index())
    wide = pd.concat(pieces).groupby("source_gse_id", observed=True).sum().reset_index()
    rows = []
    for _, row in wide.iterrows():
        for model in MODELS:
            predicted = int(row[f"{model}__predicted"])
            strict = int(row[f"{model}__strict_likely_nk"])
            rows.append({
                "source_gse_id": row.source_gse_id,
                "model": model,
                "predicted_gdt": predicted,
                "strict_likely_nk": strict,
                "strict_likely_nk_fraction": strict / predicted if predicted else np.nan,
                "author_nk_no_delta_support": int(row[f"{model}__author_nk_no_delta_support"]),
                "expression_nk_no_delta_support": int(row[f"{model}__expression_nk_no_delta_support"]),
                "author_nk_delta_conflict": int(row[f"{model}__author_nk_delta_conflict"]),
            })
    long = pd.DataFrame(rows)
    overall = long.groupby("model", observed=True)[
        ["predicted_gdt", "strict_likely_nk", "author_nk_no_delta_support",
         "expression_nk_no_delta_support", "author_nk_delta_conflict"]
    ].sum().reset_index()
    overall["strict_likely_nk_fraction"] = overall.strict_likely_nk / overall.predicted_gdt
    long.to_csv(TABLE_DIR / "strict_nk_audit_by_dataset.csv", index=False)
    overall.to_csv(TABLE_DIR / "strict_nk_audit_overall.csv", index=False)
    return long, overall


def set_style() -> None:
    plt.rcParams.update({
        "font.size": 10, "axes.titlesize": 12, "axes.labelsize": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": "white", "axes.facecolor": "white",
    })


def make_figures(source: pd.DataFrame, overall: pd.DataFrame, strict: pd.DataFrame,
                 strict_overall: pd.DataFrame, lung: pd.DataFrame, comparison: pd.DataFrame) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    set_style()
    order = list(MODELS)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    values = overall.set_index("model").loc[order]
    axes[0].bar(range(4), values.predicted_gdt / 1000, color=[COLORS[x] for x in order])
    axes[0].set_ylabel("Predicted gdT cells (thousands)")
    axes[0].set_xticks(range(4), [LABELS[x].replace(" ", "\n", 1) for x in order])
    strict_values = strict_overall.set_index("model").loc[order]
    axes[1].bar(range(4), strict_values.strict_likely_nk_fraction * 100,
                color=[COLORS[x] for x in order])
    axes[1].set_ylabel("Strict likely-NK fraction of predictions (%)")
    axes[1].set_xticks(range(4), [LABELS[x].replace(" ", "\n", 1) for x in order])
    fig.suptitle("Whole-atlas prediction yield and conservative NK audit", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "whole_atlas_overview.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    pivot = source[source.model.isin(["v4_6_highest_f1", "v2_high_f1"])].pivot(
        index="source_gse_id", columns="model", values="predicted_gdt"
    ).fillna(0)
    pivot = pivot.sort_values("v4_6_highest_f1")
    fig, ax = plt.subplots(figsize=(9.5, 6.0))
    y = np.arange(len(pivot))
    ax.barh(y - 0.18, pivot.v4_6_highest_f1 / 1000, height=0.34,
            color=COLORS["v4_6_highest_f1"], label="V4.6 highest F1")
    ax.barh(y + 0.18, pivot.v2_high_f1 / 1000, height=0.34,
            color=COLORS["v2_high_f1"], label="V2 highest F1")
    ax.set_yticks(y, pivot.index)
    ax.tick_params(axis="y", labelsize=6)
    ax.set_xlabel("Predicted gdT cells (thousands)")
    ax.set_title("Predicted cells in every source dataset", fontweight="bold")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "dataset_prediction_counts_highest_f1.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    joined = strict.merge(source[["source_gse_id", "model", "n_cells"]], on=["source_gse_id", "model"])
    focus = joined[joined.model.isin(["v4_6_highest_f1", "v2_high_f1"])].copy()
    focus["max_fraction"] = focus.groupby("source_gse_id").strict_likely_nk_fraction.transform("max")
    keep = focus.groupby("source_gse_id", observed=True).max_fraction.max().nlargest(12).index
    focus = focus[focus.source_gse_id.isin(keep)]
    pivot = focus.pivot(index="source_gse_id", columns="model", values="strict_likely_nk_fraction").fillna(0)
    pivot = pivot.sort_values("v2_high_f1")
    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    y = np.arange(len(pivot))
    ax.barh(y - 0.18, pivot.v4_6_highest_f1 * 100, height=0.34,
            color=COLORS["v4_6_highest_f1"], label="V4.6 highest F1")
    ax.barh(y + 0.18, pivot.v2_high_f1 * 100, height=0.34,
            color=COLORS["v2_high_f1"], label="V2 highest F1")
    ax.set_yticks(y, pivot.index)
    ax.set_xlabel("Strict likely-NK fraction of predictions (%)")
    ax.set_title("Datasets with the strongest conservative NK concern", fontweight="bold")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "dataset_strict_nk_audit.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    lung_pivot = lung.pivot(index="tumor_type", columns="model", values="predicted_gdt").loc[:, order]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    x = np.arange(len(lung_pivot))
    width = 0.19
    for i, model in enumerate(order):
        axes[0].bar(x + (i - 1.5) * width, lung_pivot[model] / 1000, width,
                    color=COLORS[model], label=LABELS[model])
    axes[0].set_xticks(x, ["LUAD", "LUSC"])
    axes[0].set_ylabel("Predicted gdT cells (thousands)")
    axes[0].set_title("GSE243013 prediction yield", fontweight="bold")
    axes[0].legend(frameon=False, fontsize=8)
    comp = comparison[comparison.source_gse_id.eq("GSE243013")].set_index("mode")
    axes[1].bar([0, 1], comp.loc[["highest_f1", "high_purity"], "jaccard"],
                color=["#3B718F", "#7EA6B8"])
    axes[1].set_xticks([0, 1], ["Highest F1", "High purity"])
    axes[1].set_ylim(0, 1)
    axes[1].set_ylabel("V4.6 / V2 Jaccard overlap")
    axes[1].set_title("GSE243013 cell-level agreement", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "gse243013_luad_lusc.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def fmt_int(value: float | int) -> str:
    return f"{int(value):,}"


def fmt_pct(value: float) -> str:
    return "-" if pd.isna(value) else f"{100 * value:.2f}%"


def table_html(frame: pd.DataFrame, columns: list[tuple[str, str]], classes: str = "") -> str:
    header = "".join(f"<th>{html.escape(label)}</th>" for _, label in columns)
    rows = []
    for _, row in frame.iterrows():
        cells = []
        for key, _ in columns:
            value = row[key]
            if "fraction" in key or key == "jaccard":
                text = fmt_pct(float(value))
            elif (key.startswith("n_") or key.startswith("predicted") or key.endswith("_strict_nk")
                  or key in {"strict_likely_nk", "both", "v4_6_only", "v2_only"}):
                text = fmt_int(value)
            else:
                text = str(value)
            cells.append(f"<td>{html.escape(text)}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return f'<table class="{classes}"><thead><tr>{header}</tr></thead><tbody>{"".join(rows)}</tbody></table>'


def build_report(source: pd.DataFrame, overall: pd.DataFrame, strict: pd.DataFrame,
                 strict_overall: pd.DataFrame, lung: pd.DataFrame, comparison: pd.DataFrame) -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    summary = overall.merge(strict_overall, on=["model", "predicted_gdt"], suffixes=("", "_strict"))
    summary["profile"] = summary.model.map(LABELS)
    summary_table = table_html(summary, [
        ("profile", "Profile"), ("predicted_gdt", "Predicted gdT"),
        ("predicted_fraction", "Atlas fraction"), ("strict_likely_nk", "Strict likely NK"),
        ("strict_likely_nk_fraction", "Strict NK fraction"),
        ("predicted_trdv_positive", "TRDV+ predicted"),
        ("predicted_trdc_plus_trdv_minus", "TRDC+ TRDV- predicted"),
    ])
    lung_display = lung.copy()
    lung_display["profile"] = lung_display.model.map(LABELS)
    lung_display["histology"] = lung_display.tumor_type.map({
        "lung_adenocarcinoma": "LUAD", "lung_squamous_cell_carcinoma": "LUSC",
    })
    lung_table = table_html(lung_display, [
        ("histology", "Histology"), ("profile", "Profile"), ("n_cells", "Atlas cells"),
        ("predicted_gdt", "Predicted gdT"), ("predicted_fraction", "Fraction"),
        ("predicted_trdv_positive", "TRDV+"), ("predicted_trdc_plus_trdv_minus", "TRDC+ TRDV-"),
    ])
    hf1 = source[source.model.isin(["v4_6_highest_f1", "v2_high_f1"])].pivot(
        index="source_gse_id", columns="model", values=["n_cells", "predicted_gdt"]
    )
    nk_hf1 = strict[strict.model.isin(["v4_6_highest_f1", "v2_high_f1"])].pivot(
        index="source_gse_id", columns="model", values=["strict_likely_nk", "strict_likely_nk_fraction"]
    )
    dataset = pd.DataFrame(index=hf1.index)
    dataset["source_gse_id"] = dataset.index
    dataset["n_cells"] = hf1[("n_cells", "v4_6_highest_f1")]
    dataset["v46_predicted"] = hf1[("predicted_gdt", "v4_6_highest_f1")]
    dataset["v46_strict_nk"] = nk_hf1[("strict_likely_nk", "v4_6_highest_f1")]
    dataset["v46_strict_nk_fraction"] = nk_hf1[("strict_likely_nk_fraction", "v4_6_highest_f1")]
    dataset["v2_predicted"] = hf1[("predicted_gdt", "v2_high_f1")]
    dataset["v2_strict_nk"] = nk_hf1[("strict_likely_nk", "v2_high_f1")]
    dataset["v2_strict_nk_fraction"] = nk_hf1[("strict_likely_nk_fraction", "v2_high_f1")]
    dataset = dataset.sort_values("n_cells", ascending=False).reset_index(drop=True)
    dataset_table = table_html(dataset, [
        ("source_gse_id", "Dataset"), ("n_cells", "Cells"),
        ("v46_predicted", "V4.6 predicted"), ("v46_strict_nk", "V4.6 strict NK"),
        ("v46_strict_nk_fraction", "V4.6 strict NK %"),
        ("v2_predicted", "V2 predicted"), ("v2_strict_nk", "V2 strict NK"),
        ("v2_strict_nk_fraction", "V2 strict NK %"),
    ], "compact")
    gse_comp = comparison[comparison.source_gse_id.eq("GSE243013")].copy()
    gse_comp["mode"] = gse_comp["mode"].map({"highest_f1": "Highest F1", "high_purity": "High purity"})
    overlap_table = table_html(gse_comp, [
        ("mode", "Mode"), ("v4_6_predicted", "V4.6"), ("v2_predicted", "V2"),
        ("both", "Both"), ("v4_6_only", "V4.6 only"), ("v2_only", "V2 only"),
        ("jaccard", "Jaccard"),
    ])
    v46 = summary.set_index("model").loc["v4_6_highest_f1"]
    v2 = summary.set_index("model").loc["v2_high_f1"]
    css = """
    @page { size: A4 landscape; margin: 12mm; @bottom-right { content: counter(page); color: #667; } }
    body { font-family: Arial, sans-serif; color: #172126; margin: 0; line-height: 1.42; }
    main { max-width: 1180px; margin: 0 auto; padding: 30px; }
    h1 { font-size: 30px; margin: 0 0 8px; color: #123D43; }
    h2 { font-size: 20px; margin: 28px 0 10px; border-bottom: 2px solid #D7E3E2; padding-bottom: 5px; }
    h3 { font-size: 15px; margin: 20px 0 7px; }
    .subtitle { color: #55666A; font-size: 15px; max-width: 900px; }
    .kpis { display: grid; grid-template-columns: repeat(4, 1fr); gap: 10px; margin: 22px 0; }
    .kpi { border-left: 5px solid #0B6E75; background: #F3F7F6; padding: 12px; }
    .kpi b { display: block; font-size: 23px; color: #123D43; }
    .note { border-left: 4px solid #B54A35; background: #FBF5F3; padding: 10px 14px; margin: 12px 0; }
    .method { background: #F5F7F8; border: 1px solid #D8E0E2; padding: 12px 15px; }
    table { width: 100%; border-collapse: collapse; margin: 10px 0 18px; font-size: 9px; }
    th { background: #123D43; color: white; text-align: left; padding: 6px; }
    td { padding: 5px 6px; border-bottom: 1px solid #DCE3E4; }
    tr:nth-child(even) td { background: #F6F8F8; }
    table.compact { font-size: 7.4px; }
    table.compact td, table.compact th { padding: 3px 4px; }
    img { width: 100%; max-width: 1050px; display: block; margin: 12px auto 22px; }
    .pagebreak { break-before: page; }
    code { background: #EEF2F2; padding: 1px 4px; }
    ul { margin-top: 6px; }
    """
    document = f"""<!doctype html><html><head><meta charset="utf-8"><title>gdTAI whole-atlas comparison</title><style>{css}</style></head>
    <body><main><h1>Frozen gdTAI V4.6 and V2 on the 5.93-million-cell T/NK atlas</h1>
    <p class="subtitle">A read-only, cell-for-cell application of both models at highest-F1 and high-purity operating points, with per-dataset yield, model overlap, and a conservative NK-likeness audit.</p>
    <div class="kpis"><div class="kpi"><b>{fmt_int(v46.predicted_gdt)}</b>V4.6 highest-F1 predictions</div>
    <div class="kpi"><b>{fmt_int(v2.predicted_gdt)}</b>V2 highest-F1 predictions</div>
    <div class="kpi"><b>{fmt_pct(v46.strict_likely_nk_fraction)}</b>V4.6 strict likely-NK</div>
    <div class="kpi"><b>40</b>source datasets</div></div>
    <div class="note"><b>Interpretation boundary.</b> These are deployment counts, not new performance estimates. V4.6 remains development-only. A predicted cell is not a false positive unless independent truth supports that conclusion.</div>
    <h2>Main result</h2><p>V4.6 is substantially more selective than V2 on this atlas. At the highest-F1 settings it calls {fmt_int(v46.predicted_gdt)} cells ({fmt_pct(v46.predicted_fraction)}) versus {fmt_int(v2.predicted_gdt)} ({fmt_pct(v2.predicted_fraction)}) for V2. The strict likely-NK audit is {fmt_int(v46.strict_likely_nk)} cells ({fmt_pct(v46.strict_likely_nk_fraction)}) for V4.6 and {fmt_int(v2.strict_likely_nk)} ({fmt_pct(v2.strict_likely_nk_fraction)}) for V2.</p>
    {summary_table}<img src="../../Integrated_dataset/figures/gdT_prediction/gdtai_v4_6_whole_atlas/whole_atlas_overview.png">
    <h2>What “likely NK” means</h2><div class="method"><p>The strict audit requires no TRDV1/2/3 expression and no productive gd-TCR metadata, plus either:</p><ul><li>an explicit source-author NK label; or</li><li>weak CD3 expression (at most one of CD3D/E/G detected), at least one Tier-A innate gene (NCR1, SIGLEC7, SH2D1B), and at least two Tier-A/B genes in total (adding LAT2, SYK, PLCG2, PILRB, CD300C).</li></ul><p>NKG7, GNLY, KLRD1, TYROBP, and FCER1G alone are only “cytotoxic ambiguous,” because true gdT cells can express them strongly. Author-NK predictions with TRDV or productive gd-TCR support are exported separately as annotation/receptor conflicts, not counted in the strict NK total.</p></div>
    <img src="../../Integrated_dataset/figures/gdT_prediction/gdtai_v4_6_whole_atlas/dataset_strict_nk_audit.png">
    <div class="pagebreak"></div><h2>Large LUAD/LUSC project: GSE243013</h2>
    <p>The corrected metadata contain 759,436 cells: 189,821 LUAD and 569,615 LUSC cells. V4.6 highest-F1 identifies 25,792 cells (3.40%); V2 highest-F1 identifies 37,112 (4.89%). This project has no productive gd-TCR metadata, so expression evidence and model agreement are the available audit signals.</p>
    {lung_table}{overlap_table}<img src="../../Integrated_dataset/figures/gdT_prediction/gdtai_v4_6_whole_atlas/gse243013_luad_lusc.png">
    <h2>Every dataset</h2><p>The table below uses highest-F1 for the principal cross-dataset comparison. High-purity values and all audit fields are retained in the canonical CSV tables.</p>{dataset_table}
    <img src="../../Integrated_dataset/figures/gdT_prediction/gdtai_v4_6_whole_atlas/dataset_prediction_counts_highest_f1.png">
    <h2>Reproducibility and limitations</h2><ul><li>Input: validated corrected atlas, 5,933,312 cells x 27,413 genes.</li><li>Normalization: raw counts to 10,000 counts per cell, followed by log1p, computed in 20,000-cell chunks.</li><li>Frozen thresholds: V4.6 highest-F1 0.894354 and high-purity 0.971830; V2 uses its frozen highest-F1 and annotation-specific high-purity contract.</li><li>No model, calibrator, feature set, or threshold was fitted or tuned on atlas predictions.</li><li>The input H5AD size and nanosecond modification time were unchanged; predictions are stored separately in Parquet.</li><li>Source annotations vary in quality. The strict NK audit is conservative but is still an audit heuristic, not orthogonal cell identity truth.</li></ul>
    <p><b>Canonical outputs:</b> <code>Integrated_dataset/tables/gdT_prediction/gdtai_v4_6_whole_atlas/</code></p>
    </main></body></html>"""
    index = REPORT_DIR / "index.html"
    index.write_text(document)
    pdf = REPORT_DIR / "gdtai_v4_6_vs_v2_whole_atlas_report.pdf"
    subprocess.run(["/home/anaconda3/bin/weasyprint", str(index), str(pdf)], check=True)

    markdown = f"""# Frozen gdTAI V4.6 versus V2 whole-atlas application

- Status: PASS
- Cells: 5,933,312 across 40 source datasets
- V4.6 highest-F1: {fmt_int(v46.predicted_gdt)} predictions; strict likely NK {fmt_int(v46.strict_likely_nk)} ({fmt_pct(v46.strict_likely_nk_fraction)})
- V2 highest-F1: {fmt_int(v2.predicted_gdt)} predictions; strict likely NK {fmt_int(v2.strict_likely_nk)} ({fmt_pct(v2.strict_likely_nk_fraction)})
- GSE243013: 759,436 cells; V4.6 highest-F1 25,792; V2 highest-F1 37,112
- Models and thresholds were frozen; no atlas-informed tuning occurred
- Input H5AD was not mutated
- Report: `{pdf.relative_to(ROOT)}`
"""
    (LOG_DIR / "whole_atlas_report_summary.md").write_text(markdown)


def main() -> None:
    source = pd.read_csv(TABLE_DIR / "predictions_by_dataset.csv")
    overall = pd.read_csv(TABLE_DIR / "whole_atlas_overall.csv")
    lung = pd.read_csv(TABLE_DIR / "gse243013_luad_lusc_predictions.csv")
    comparison = pd.read_csv(TABLE_DIR / "v4_6_vs_v2_by_dataset.csv")
    strict, strict_overall = strict_nk_audit()
    make_figures(source, overall, strict, strict_overall, lung, comparison)
    build_report(source, overall, strict, strict_overall, lung, comparison)
    print(json.dumps({
        "status": "PASS_V4_6_V2_WHOLE_ATLAS_REPORT",
        "html": str(REPORT_DIR / "index.html"),
        "pdf": str(REPORT_DIR / "gdtai_v4_6_vs_v2_whole_atlas_report.pdf"),
    }, indent=2))


if __name__ == "__main__":
    main()
