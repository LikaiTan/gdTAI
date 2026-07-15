#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Failure-mode audit for gdTAI v3 TRDC/NK guard predictions.

The script consumes prediction tables from
``run_gdtai_v3_trdc_nk_guard_classifier.py`` and writes read-only audit
summaries stratified by TRDC/TRDV quadrant, labels, NK/TCRAB/CD4-like stress
groups, TCR-gene evidence, CD3/NK scores, QC metrics, source/sample/tissue, and
candidate strategy.
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
import shutil
import subprocess
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from run_gdt_deg_tcr_classifier_training import safe_div
from run_gdt_prediction_package_evaluation import dataframe_to_html, dataframe_to_markdown


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
STATIC_DIR = PROJECT_ROOT / "gdT_prediction"
OUT_PREFIX = "gdtai_v3_trdc_nk_guard"
TABLE_DIR = OUTPUT_ROOT / "tables" / "gdT_prediction" / OUT_PREFIX
FIGURE_DIR = OUTPUT_ROOT / "figures" / "gdT_prediction" / OUT_PREFIX
LOG_DIR = OUTPUT_ROOT / "logs" / "gdT_prediction" / OUT_PREFIX
REPORT_MD = LOG_DIR / "gdtai_v3_trdc_nk_failure_audit.md"
REPORT_HTML = STATIC_DIR / "gdtai_v3_trdc_nk_failure_audit.html"
REPORT_PDF = STATIC_DIR / "gdtai_v3_trdc_nk_failure_audit.pdf"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit gdTAI v3 TRDC/NK failure modes.")
    parser.add_argument("--external-predictions", type=Path, default=TABLE_DIR / "external_predictions_wide.csv.gz")
    parser.add_argument("--internal-predictions", type=Path, default=TABLE_DIR / "internal_tune_predictions_wide.csv.gz")
    parser.add_argument("--no-pdf", action="store_true")
    return parser.parse_args()


def ensure_dirs() -> None:
    for path in [TABLE_DIR, FIGURE_DIR, LOG_DIR, STATIC_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def infer_strategies(df: pd.DataFrame) -> list[str]:
    out = []
    for col in df.columns:
        if not col.endswith("_pred"):
            continue
        strategy = col[: -len("_pred")]
        if f"{strategy}_score" in df.columns:
            out.append(strategy)
    return out


def bool_array(series: pd.Series) -> np.ndarray:
    if series.dtype == bool:
        return series.to_numpy(dtype=bool)
    return series.astype("string").fillna("").str.strip().str.lower().isin({"true", "1", "yes", "t"}).to_numpy(dtype=bool)


def summarize_binary(y: np.ndarray, pred: np.ndarray) -> dict[str, Any]:
    primary = y >= 0
    yy = y[primary]
    pp = pred[primary]
    tp = int(((pp == 1) & (yy == 1)).sum())
    fp = int(((pp == 1) & (yy == 0)).sum())
    tn = int(((pp == 0) & (yy == 0)).sum())
    fn = int(((pp == 0) & (yy == 1)).sum())
    precision = safe_div(tp, tp + fp)
    recall = safe_div(tp, tp + fn)
    specificity = safe_div(tn, tn + fp)
    f1 = safe_div(2 * tp, 2 * tp + fp + fn)
    return {
        "primary_eval_cells": int(primary.sum()),
        "positives": int((yy == 1).sum()),
        "negatives": int((yy == 0).sum()),
        "predicted": int(pp.sum()),
        "tp": tp,
        "fp": fp,
        "tn": tn,
        "fn": fn,
        "precision": precision,
        "recall": recall,
        "specificity": specificity,
        "f1": f1,
    }


def numeric_summary(group: pd.DataFrame, col: str) -> dict[str, float]:
    if col not in group:
        return {f"{col}_median": np.nan, f"{col}_mean": np.nan}
    values = pd.to_numeric(group[col], errors="coerce")
    return {f"{col}_median": float(values.median()) if values.notna().any() else np.nan, f"{col}_mean": float(values.mean()) if values.notna().any() else np.nan}


def summarize_groups(df: pd.DataFrame, strategies: list[str], *, dataset: str, group_cols: list[str]) -> pd.DataFrame:
    rows = []
    y = pd.to_numeric(df.get("y_true", pd.Series([-1] * len(df))), errors="coerce").fillna(-1).to_numpy(dtype=np.int8)
    for strategy in strategies:
        pred = bool_array(df[f"{strategy}_pred"])
        for group_col in group_cols:
            if group_col not in df:
                continue
            for value, group in df.groupby(group_col, dropna=False, sort=True):
                idx = group.index.to_numpy()
                row = {
                    "dataset": dataset,
                    "strategy": strategy,
                    "group_column": group_col,
                    "group_value": str(value),
                    "n_cells": int(group.shape[0]),
                }
                row.update(summarize_binary(y[idx], pred[idx].astype(np.int8)))
                for label_col in ["real_gdt", "real_abt_tcr_strict", "tcr_paired_ab", "CD4_T_like_warning", "is_NK_guard", "is_TCRAB_guard"]:
                    if label_col in group:
                        row[f"{label_col}_cells"] = int(bool_array(group[label_col]).sum())
                if "cell_type" in group:
                    row["NK_annotation_cells"] = int(group["cell_type"].astype(str).eq("NK").sum())
                    row["B_myeloid_annotation_cells"] = int(group["cell_type"].astype(str).isin(["B_cell", "Inflammatory_myeloid"]).sum())
                if "annotation" in group:
                    row["NK_annotation_cells"] = int(group["annotation"].astype(str).str.upper().eq("NK_CELL").sum())
                for flag in ["any_TRDV", "any_TRDJ", "any_TRG", "TRDC_only"]:
                    if flag in group:
                        row[f"{flag}_cells"] = int(bool_array(group[flag]).sum())
                for num_col in [
                    "CD3_score",
                    "NK_score",
                    "gdT_TCR_score",
                    "abT_TCR_score",
                    "NK_minus_CD3_score",
                    "phase4_trd_score",
                    "phase4_trab_score",
                    "phase4_trd_minus_trab",
                    "n_counts",
                    "n_genes",
                    "pct_mito",
                ]:
                    row.update(numeric_summary(group, num_col))
                if f"{strategy}_score" in group:
                    row.update(numeric_summary(group, f"{strategy}_score"))
                rows.append(row)
    return pd.DataFrame(rows)


def trdc_comparison(df: pd.DataFrame, strategies: list[str], *, dataset: str) -> pd.DataFrame:
    rows = []
    if "tcr_gene_quadrant" not in df:
        return pd.DataFrame()
    y = pd.to_numeric(df.get("y_true", pd.Series([-1] * len(df))), errors="coerce").fillna(-1).to_numpy(dtype=np.int8)
    for strategy in strategies:
        pred = bool_array(df[f"{strategy}_pred"])
        score = pd.to_numeric(df[f"{strategy}_score"], errors="coerce").to_numpy(dtype=float)
        for quadrant, group in df.groupby("tcr_gene_quadrant", dropna=False, sort=True):
            idx = group.index.to_numpy()
            primary = y[idx] >= 0
            fp = pred[idx] & primary & (y[idx] == 0)
            fn = (~pred[idx]) & primary & (y[idx] == 1)
            rows.append(
                {
                    "dataset": dataset,
                    "strategy": strategy,
                    "tcr_gene_quadrant": str(quadrant),
                    "n_cells": int(group.shape[0]),
                    "primary_eval_cells": int(primary.sum()),
                    "predicted": int(pred[idx].sum()),
                    "fp": int(fp.sum()),
                    "fn": int(fn.sum()),
                    "median_score": float(np.nanmedian(score[idx])) if idx.size else np.nan,
                    "real_gdt_cells": int(bool_array(group["real_gdt"]).sum()) if "real_gdt" in group else int((y[idx] == 1).sum()),
                    "real_abt_cells": int(bool_array(group["real_abt_tcr_strict"]).sum()) if "real_abt_tcr_strict" in group else int((y[idx] == 0).sum()),
                    "unlabeled_or_stress_cells": int((~primary).sum()),
                    "NK_cells": int(group["cell_type"].astype(str).eq("NK").sum()) if "cell_type" in group else int(group.get("annotation", pd.Series([], dtype=str)).astype(str).str.upper().eq("NK_CELL").sum()),
                    "paired_TCRAB_cells": int(bool_array(group["tcr_paired_ab"]).sum()) if "tcr_paired_ab" in group else int(bool_array(group["is_TCRAB_guard"]).sum()) if "is_TCRAB_guard" in group else 0,
                    "CD3_score_median": float(pd.to_numeric(group.get("CD3_score", pd.Series(dtype=float)), errors="coerce").median()) if "CD3_score" in group else np.nan,
                    "NK_score_median": float(pd.to_numeric(group.get("NK_score", pd.Series(dtype=float)), errors="coerce").median()) if "NK_score" in group else np.nan,
                    "TRDJ_evidence_cells": int(bool_array(group["any_TRDJ"]).sum()) if "any_TRDJ" in group else 0,
                    "TRG_evidence_cells": int(bool_array(group["any_TRG"]).sum()) if "any_TRG" in group else 0,
                }
            )
    return pd.DataFrame(rows)


def strategy_delta_table(external: pd.DataFrame, strategies: list[str]) -> pd.DataFrame:
    if "v2_high_purity" not in strategies:
        return pd.DataFrame()
    rows = []
    baseline = bool_array(external["v2_high_purity_pred"])
    groups = {
        "all_cells": np.ones(external.shape[0], dtype=bool),
        "NK_cells": external["cell_type"].astype(str).eq("NK").to_numpy(dtype=bool) if "cell_type" in external else np.zeros(external.shape[0], dtype=bool),
        "paired_TCRAB_cells": bool_array(external["tcr_paired_ab"]) if "tcr_paired_ab" in external else np.zeros(external.shape[0], dtype=bool),
        "TRDC_plus_TRDV_minus": external["tcr_gene_quadrant"].astype(str).eq("TRDC+TRDV-").to_numpy(dtype=bool) if "tcr_gene_quadrant" in external else np.zeros(external.shape[0], dtype=bool),
        "CD4_Treg_warning": bool_array(external["CD4_T_like_warning"]) if "CD4_T_like_warning" in external else np.zeros(external.shape[0], dtype=bool),
        "B_myeloid": external["cell_type"].astype(str).isin(["B_cell", "Inflammatory_myeloid"]).to_numpy(dtype=bool) if "cell_type" in external else np.zeros(external.shape[0], dtype=bool),
    }
    for strategy in strategies:
        pred = bool_array(external[f"{strategy}_pred"])
        for group, mask in groups.items():
            rows.append(
                {
                    "strategy": strategy,
                    "comparison_baseline": "v2_high_purity",
                    "group": group,
                    "group_cells": int(mask.sum()),
                    "predicted_cells": int((pred & mask).sum()),
                    "baseline_predicted_cells": int((baseline & mask).sum()),
                    "delta_vs_baseline": int((pred & mask).sum()) - int((baseline & mask).sum()),
                }
            )
    return pd.DataFrame(rows)


def plot_quadrant_summary(quadrant: pd.DataFrame) -> list[Path]:
    paths: list[Path] = []
    if quadrant.empty:
        return paths
    external = quadrant[quadrant["dataset"] == "external"].copy()
    if external.empty:
        return paths
    pivot = external.pivot_table(index="strategy", columns="tcr_gene_quadrant", values="predicted", aggfunc="sum").fillna(0)
    fig, ax = plt.subplots(figsize=(10.2, 5.0), constrained_layout=True)
    pivot.plot(kind="bar", ax=ax, width=0.82)
    ax.set_ylabel("predicted cells")
    ax.set_title("External predictions by TRDC/TRDV quadrant")
    ax.tick_params(axis="x", labelrotation=25)
    ax.legend(frameon=False, fontsize=8)
    path = FIGURE_DIR / "v3_failure_audit_trdc_trdv_predictions.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    fp_pivot = external.pivot_table(index="strategy", columns="tcr_gene_quadrant", values="fp", aggfunc="sum").fillna(0)
    fig, ax = plt.subplots(figsize=(10.2, 5.0), constrained_layout=True)
    fp_pivot.plot(kind="bar", ax=ax, width=0.82)
    ax.set_ylabel("false-positive cells")
    ax.set_title("External false positives by TRDC/TRDV quadrant")
    ax.tick_params(axis="x", labelrotation=25)
    ax.legend(frameon=False, fontsize=8)
    path = FIGURE_DIR / "v3_failure_audit_trdc_trdv_false_positives.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)
    return paths


def write_report(
    *,
    external_path: Path,
    internal_path: Path,
    quadrant: pd.DataFrame,
    groups: pd.DataFrame,
    deltas: pd.DataFrame,
    figures: list[Path],
) -> None:
    top_external = quadrant[quadrant["dataset"] == "external"].sort_values(["strategy", "predicted"], ascending=[True, False]) if not quadrant.empty else quadrant
    top_groups = groups[(groups["dataset"] == "external") & (groups["group_column"].isin(["cell_type", "tcr_strict_cell_label", "tcr_gene_quadrant"]))] if not groups.empty else groups
    lines = [
        "# gdTAI v3 TRDC/NK Failure Audit",
        "",
        f"- External predictions: `{external_path}`",
        f"- Internal tune predictions: `{internal_path}`",
        "",
        "## TRDC/TRDV Quadrant Summary",
        dataframe_to_markdown(top_external),
        "",
        "## Strategy Delta Versus v2 High-Purity",
        dataframe_to_markdown(deltas),
        "",
        "## External Label And Stress Groups",
        dataframe_to_markdown(top_groups.head(120) if not top_groups.empty else top_groups),
        "",
        "## Outputs",
        f"- Tables: `{TABLE_DIR}`",
        f"- Figures: `{FIGURE_DIR}`",
        f"- HTML report: `{REPORT_HTML}`",
        f"- PDF report: `{REPORT_PDF}`",
    ]
    REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    asset_dir = STATIC_DIR / "assets" / OUT_PREFIX
    asset_dir.mkdir(parents=True, exist_ok=True)
    fig_html = []
    for fig in figures:
        target = asset_dir / fig.name
        if fig.exists() and fig.resolve() != target.resolve():
            shutil.copyfile(fig, target)
        fig_html.append(f"<figure><img src='assets/{OUT_PREFIX}/{html.escape(fig.name)}'><figcaption>{html.escape(fig.stem.replace('_', ' '))}</figcaption></figure>")
    css = """
    @page{size:A4 landscape;margin:8mm}body{font-family:Arial,Helvetica,sans-serif;margin:18px;color:#20272e;background:#f5f6f7;line-height:1.45}main{max-width:1320px;margin:auto}section{background:white;border:1px solid #d9dee4;padding:14px;margin:12px 0;break-inside:avoid}h1{font-size:28px;margin:0 0 8px}h2{font-size:19px}table{border-collapse:collapse;width:100%;font-size:8.5px;table-layout:fixed}th,td{border:1px solid #d9dee4;padding:3px 4px;text-align:left;vertical-align:top;overflow-wrap:anywhere}th{background:#eef1f4}img{max-width:100%;border:1px solid #d9dee4}code{background:#eef1f4;padding:1px 4px;border-radius:3px}
    """
    html_doc = f"""<!doctype html><html><head><meta charset='utf-8'><title>gdTAI v3 Failure Audit</title><style>{css}</style></head><body><main>
    <section><h1>gdTAI v3 TRDC/NK Failure Audit</h1><p>Predictions are stratified by TRDC/TRDV quadrant, gold/stress labels, NK and paired-TCRAB annotations, CD3/NK scores, and QC metrics.</p></section>
    <section><h2>TRDC/TRDV Quadrant Summary</h2>{dataframe_to_html(top_external)}</section>
    <section><h2>Delta Versus v2 High-Purity</h2>{dataframe_to_html(deltas)}</section>
    <section><h2>External Label And Stress Groups</h2>{dataframe_to_html(top_groups.head(120) if not top_groups.empty else top_groups)}</section>
    <section><h2>Figures</h2>{''.join(fig_html)}</section>
    </main></body></html>"""
    REPORT_HTML.write_text(html_doc, encoding="utf-8")


def render_pdf(no_pdf: bool) -> None:
    if no_pdf:
        return
    subprocess.run(
        [
            "google-chrome",
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
    args = parse_args()
    ensure_dirs()
    if not args.external_predictions.exists():
        raise FileNotFoundError(f"Missing external prediction table: {args.external_predictions}")
    external = pd.read_csv(args.external_predictions)
    external_strategies = infer_strategies(external)
    if not external_strategies:
        raise RuntimeError("No strategy prediction columns found in external prediction table.")

    group_tables = []
    quadrant_tables = []
    external_groups = summarize_groups(
        external,
        external_strategies,
        dataset="external",
        group_cols=["tcr_gene_quadrant", "cell_type", "tcr_strict_cell_label", "tissue", "sample_id"],
    )
    group_tables.append(external_groups)
    quadrant_tables.append(trdc_comparison(external, external_strategies, dataset="external"))
    deltas = strategy_delta_table(external, external_strategies)
    deltas.to_csv(TABLE_DIR / "v3_failure_audit_strategy_delta_vs_v2_high_purity.csv", index=False)

    if args.internal_predictions.exists():
        internal = pd.read_csv(args.internal_predictions)
        internal_strategies = infer_strategies(internal)
        if internal_strategies:
            group_tables.append(
                summarize_groups(
                    internal,
                    internal_strategies,
                    dataset="internal_tune",
                    group_cols=["tcr_gene_quadrant", "source_gse_id", "annotation", "tissue"],
                )
            )
            quadrant_tables.append(trdc_comparison(internal, internal_strategies, dataset="internal_tune"))

    groups = pd.concat(group_tables, ignore_index=True, sort=False)
    quadrant = pd.concat(quadrant_tables, ignore_index=True, sort=False)
    groups.to_csv(TABLE_DIR / "v3_failure_audit_group_summary.csv", index=False)
    quadrant.to_csv(TABLE_DIR / "v3_failure_audit_trdc_trdv_quadrant_summary.csv", index=False)
    figures = plot_quadrant_summary(quadrant)
    write_report(
        external_path=args.external_predictions,
        internal_path=args.internal_predictions,
        quadrant=quadrant,
        groups=groups,
        deltas=deltas,
        figures=figures,
    )
    render_pdf(args.no_pdf)


if __name__ == "__main__":
    main()
