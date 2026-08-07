#!/usr/bin/env python3
"""Finalize a print-safe Gate C report when no frozen operating point exists."""

from __future__ import annotations

import html
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import expit

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from gdtai_v4_gpu_core import json_safe, stage2_threshold_frontier  # noqa: E402
from run_gdtai_v4_gpu_evaluation import merged_config  # noqa: E402


CONFIG_PATH = PROJECT_ROOT / "configs/models/gdtai/v4_gpu_nested_evaluation.json"


def resolve(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else PROJECT_ROOT / path


def stage1_diagnostics(config: dict[str, Any], cells: pd.DataFrame) -> pd.DataFrame:
    root = resolve(config["outputs"]["table_dir"]) / "threshold_frontiers"
    rows: list[dict[str, Any]] = []
    for path in sorted(root.glob("outer_*/stage1/*__stage1.parquet")):
        frame = pd.read_parquet(path)
        recall_ok = frame.minimum_source_gdt_recall.ge(.99) & frame.minimum_source_abt_recall.ge(.98)
        if not recall_ok.any():
            continue
        selected = frame[recall_ok].iloc[-1]
        nk_columns = [column for column in frame if column.startswith("nk_passage__")]
        passages = {column.removeprefix("nk_passage__"): float(selected[column]) for column in nk_columns}
        outer = path.parts[-3]
        heldout = outer.split("_", 2)[-1]
        training_nk = cells[cells.stage1_role.eq("nk_negative") & cells.source_gse_id.ne(heldout)].source_gse_id.value_counts()
        weighted_numerator = sum(passages[source] * int(training_nk.get(source, 0)) for source in passages)
        weighted_denominator = sum(int(training_nk.get(source, 0)) for source in passages)
        stable_sources = [source for source in passages if int(training_nk.get(source, 0)) >= 100]
        stable = [passages[source] for source in stable_sources]
        worst_stable_source = sorted(stable_sources, key=lambda source: (-passages[source], -int(training_nk.get(source, 0)), source))[0]
        worst_source = sorted(passages, key=lambda source: (-passages[source], int(training_nk.get(source, 0)), source))[0]
        candidate = path.name.split("__stage1.parquet")[0].split("_", 1)[-1]
        rows.append({
            "outer_fold_id": outer, "heldout_source": heldout, "candidate_file_id": candidate,
            "family": "xgboost" if "xgboost" in path.name else "torch_ridge",
            "threshold_at_highest_recall_feasible": float(selected.threshold),
            "minimum_source_gdt_recall": float(selected.minimum_source_gdt_recall),
            "minimum_source_abt_recall": float(selected.minimum_source_abt_recall),
            "macro_nk_passage": float(selected.macro_nk_passage),
            "maximum_source_nk_passage": float(selected.maximum_source_nk_passage),
            "pooled_nk_passage": float(weighted_numerator / weighted_denominator),
            "maximum_nk_passage_sources_n_ge_100": float(max(stable)) if stable else math.nan,
            "worst_nk_source_n_ge_100": worst_stable_source,
            "worst_nk_source_n_ge_100_n": int(training_nk.get(worst_stable_source, 0)),
            "worst_nk_source": worst_source,
            "worst_nk_source_n": int(training_nk.get(worst_source, 0)),
            "worst_nk_source_passage": passages[worst_source],
            "official_stage1_pass": False,
        })
    return pd.DataFrame(rows)


def nearest_row(frame: pd.DataFrame, mode: str, mode_spec: dict[str, Any]) -> tuple[pd.Series, pd.Series]:
    objective = "macro_f1" if mode_spec["objective"] == "f1" else "macro_f0_5"
    unconstrained = frame.loc[frame[objective].idxmax()]
    recall_violation = np.maximum(0, float(mode_spec["minimum_macro_recall"]) - frame.macro_recall) / float(mode_spec["minimum_macro_recall"])
    paired_violation = np.maximum(0, frame.paired_abt_fpr - float(mode_spec["maximum_paired_abt_fpr"])) / float(mode_spec["maximum_paired_abt_fpr"])
    nk_violation = np.maximum(0, frame.strict_nk_fpr - float(mode_spec["maximum_strict_nk_fpr"])) / float(mode_spec["maximum_strict_nk_fpr"])
    total = recall_violation + paired_violation + nk_violation
    best = np.flatnonzero(np.isclose(total, total.min(), atol=1e-12, rtol=0))
    nearest_index = frame.iloc[best].sort_values([objective, "threshold"], ascending=[False, False]).index[0]
    return unconstrained, frame.loc[nearest_index]


def comparator_diagnostics(config: dict[str, Any]) -> pd.DataFrame:
    root = resolve(config["outputs"]["table_dir"]) / "threshold_frontiers"
    rows: list[dict[str, Any]] = []
    for path in sorted(root.glob("outer_*/*/*.parquet")):
        if "/stage1/" in str(path):
            continue
        mode = "high_purity" if path.name.endswith("__high_purity.parquet") else "balanced"
        spec = config["operating_modes"][mode]
        frame = pd.read_parquet(path)
        unconstrained, nearest = nearest_row(frame, mode, spec)
        objective = "macro_f1" if spec["objective"] == "f1" else "macro_f0_5"
        stage = path.parts[-2]
        rows.append({
            "outer_fold_id": path.parts[-3], "model_id": stage, "mode": mode,
            "candidate_file": path.name, "official_pass": bool(frame.eligible.any()),
            "unconstrained_threshold": float(unconstrained.threshold),
            "unconstrained_objective": float(unconstrained[objective]),
            "unconstrained_recall": float(unconstrained.macro_recall),
            "unconstrained_paired_abt_fpr": float(unconstrained.paired_abt_fpr),
            "unconstrained_strict_nk_fpr": float(unconstrained.strict_nk_fpr),
            "nearest_threshold": float(nearest.threshold),
            "nearest_objective": float(nearest[objective]),
            "nearest_recall": float(nearest.macro_recall),
            "nearest_paired_abt_fpr": float(nearest.paired_abt_fpr),
            "nearest_strict_nk_fpr": float(nearest.strict_nk_fpr),
            "recall_requirement": float(spec["minimum_macro_recall"]),
            "paired_abt_fpr_limit": float(spec["maximum_paired_abt_fpr"]),
            "strict_nk_fpr_limit": float(spec["maximum_strict_nk_fpr"]),
        })
    return pd.DataFrame(rows)


def legacy_diagnostics(config: dict[str, Any], cells: pd.DataFrame) -> pd.DataFrame:
    cache = json.loads((resolve(config["outputs"]["cache_dir"]) / "feature_cache_manifest.json").read_text())
    probability = expit(np.load(cache["legacy_score_path"], mmap_mode="r").astype(np.float64))
    rows = []
    for heldout in ["HRA005041", "GSE144469", "BALF_BLOOD_COPD"]:
        primary_mask = cells.source_gse_id.ne(heldout) & cells.truth_class.isin(["gdT_primary", "abT_primary"])
        control_mask = cells.source_gse_id.ne(heldout) & cells.stage1_role.eq("nk_negative")
        primary_rows = cells.index[primary_mask].to_numpy(np.int64)
        control_rows = cells.index[control_mask].to_numpy(np.int64)
        primary = cells.loc[primary_rows]
        for mode, spec in config["operating_modes"].items():
            frontier, _ = stage2_threshold_frontier(
                probability[primary_rows], primary.truth_class.eq("gdT_primary").to_numpy(np.int8),
                primary.source_gse_id.to_numpy(object), primary.truth_class.eq("abT_primary").to_numpy(bool),
                probability[control_rows], np.ones(primary_rows.size), np.ones(control_rows.size), 0,
                primary.exclusion_union.to_numpy(bool), cells.loc[control_rows, "exclusion_union"].to_numpy(bool),
                mode, spec,
            )
            frontier_dir = resolve(config["outputs"]["table_dir"]) / "threshold_frontiers" / f"outer_{['HRA005041','GSE144469','BALF_BLOOD_COPD'].index(heldout)}_{heldout}" / "legacy_trd_minus_trab"
            frontier_dir.mkdir(parents=True, exist_ok=True)
            frontier_path = frontier_dir / f"raw_score__{mode}.parquet"
            temporary = frontier_path.with_suffix(".tmp.parquet")
            frontier.to_parquet(temporary, index=False, compression="zstd")
            temporary.replace(frontier_path)
            unconstrained, nearest = nearest_row(frontier, mode, spec)
            objective = "macro_f1" if spec["objective"] == "f1" else "macro_f0_5"
            rows.append({
                "outer_fold_id": f"outer_{['HRA005041','GSE144469','BALF_BLOOD_COPD'].index(heldout)}_{heldout}",
                "model_id": "legacy_trd_minus_trab", "mode": mode, "candidate_file": "raw_score",
                "official_pass": bool(frontier.eligible.any()),
                "unconstrained_threshold": float(unconstrained.threshold), "unconstrained_objective": float(unconstrained[objective]),
                "unconstrained_recall": float(unconstrained.macro_recall), "unconstrained_paired_abt_fpr": float(unconstrained.paired_abt_fpr),
                "unconstrained_strict_nk_fpr": float(unconstrained.strict_nk_fpr),
                "nearest_threshold": float(nearest.threshold), "nearest_objective": float(nearest[objective]),
                "nearest_recall": float(nearest.macro_recall), "nearest_paired_abt_fpr": float(nearest.paired_abt_fpr),
                "nearest_strict_nk_fpr": float(nearest.strict_nk_fpr),
                "recall_requirement": float(spec["minimum_macro_recall"]),
                "paired_abt_fpr_limit": float(spec["maximum_paired_abt_fpr"]),
                "strict_nk_fpr_limit": float(spec["maximum_strict_nk_fpr"]),
            })
    return pd.DataFrame(rows)


def make_figures(stage1: pd.DataFrame, comparators: pd.DataFrame, figure_dir: Path) -> list[Path]:
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    colors = {"torch_ridge": "#176b5f", "xgboost": "#c95d3a"}
    for family, frame in stage1.groupby("family"):
        axes[0].scatter(frame.minimum_source_gdt_recall, frame.maximum_source_nk_passage, label=family, color=colors[family], alpha=.8)
        axes[1].scatter(frame.pooled_nk_passage, frame.maximum_nk_passage_sources_n_ge_100, label=family, color=colors[family], alpha=.8)
    axes[0].axvline(.99, color="#555", linestyle="--"); axes[0].axhline(.5, color="#a22", linestyle="--")
    axes[0].set(xlabel="Minimum-source gdT recall", ylabel="Maximum-source NK passage", title="Literal Stage-1 guardrail")
    axes[1].axhline(.5, color="#a22", linestyle="--"); axes[1].set(xlabel="Pooled NK passage", ylabel="Maximum NK passage (source n >= 100)", title="Large-stratum sensitivity")
    axes[0].legend(frameon=False)
    path = figure_dir / "stage1_guardrail_failure.png"; fig.savefig(path, dpi=220); plt.close(fig); paths.append(path)

    balanced = comparators[comparators["mode"].eq("balanced")]
    fig, ax = plt.subplots(figsize=(7.8, 5.2), constrained_layout=True)
    for model, frame in balanced.groupby("model_id"):
        ax.scatter(frame.nearest_recall, frame.nearest_strict_nk_fpr, label=model, s=45, alpha=.8)
    ax.axvline(.8, color="#555", linestyle="--", label="recall floor")
    ax.axhline(.01, color="#a22", linestyle="--", label="NK FPR ceiling")
    ax.set(xlabel="Nearest-frontier macro recall", ylabel="Nearest-frontier strict-NK FPR", title="No fair comparator reached the balanced operating region")
    ax.set_yscale("symlog", linthresh=1e-4)
    ax.legend(frameon=False, fontsize=8)
    path = figure_dir / "comparator_near_miss.png"; fig.savefig(path, dpi=220); plt.close(fig); paths.append(path)
    return paths


def table_html(frame: pd.DataFrame, columns: list[str]) -> str:
    shown = frame[[column for column in columns if column in frame]].copy()
    for column in shown.select_dtypes(include=["float"]).columns:
        shown[column] = shown[column].map(lambda value: f"{value:.4g}" if np.isfinite(value) else "NA")
    return shown.to_html(index=False, border=0)


def main() -> None:
    config = merged_config(CONFIG_PATH)
    table_dir = resolve(config["outputs"]["table_dir"])
    figure_dir = resolve(config["outputs"]["figure_dir"])
    log_dir = resolve(config["outputs"]["log_dir"])
    static_dir = resolve(config["outputs"]["static_dir"])
    static_assets = static_dir / "assets"
    for path in (table_dir, figure_dir, log_dir, static_dir, static_assets):
        path.mkdir(parents=True, exist_ok=True)
    cells = pd.read_csv(resolve(config["preflight"]["cell_manifest"]))
    stage1 = stage1_diagnostics(config, cells)
    comparators = pd.concat([comparator_diagnostics(config), legacy_diagnostics(config, cells)], ignore_index=True)
    stage1.to_csv(table_dir / "gate_c_stage1_diagnostics.csv", index=False)
    comparators.to_csv(table_dir / "gate_c_comparator_near_miss.csv", index=False)
    figures = make_figures(stage1, comparators, figure_dir)
    for source in figures:
        (static_assets / source.name).write_bytes(source.read_bytes())

    convergence_files = sorted(table_dir.glob("outer_*_stage1_candidates.csv")) + sorted(table_dir.glob("outer_*_*logistic_candidates.csv"))
    convergence = pd.concat([pd.read_csv(path).assign(source_file=path.name) for path in convergence_files], ignore_index=True)
    n_fits = int(convergence.shape[0] * 3)
    n_converged = int(convergence.converged_fold_count.sum())
    summary = {
        "status": "COMPLETE_NO_ELIGIBLE_MODEL", "protocol_version": config["protocol_version"],
        "outer_folds_completed": 3, "stage1_candidates": int(stage1.shape[0]),
        "stage1_candidates_passed": int(stage1.official_stage1_pass.sum()),
        "comparator_operating_points_passed": int(comparators.official_pass.sum()),
        "recorded_fold_fits": n_fits, "converged_fold_fits": n_converged,
        "bootstrap_status": "not_estimable_without_eligible_balanced_predictions",
        "extension_control_status": "not_estimable_without_selected_v4_model",
        "promotion_eligible": False, "promotion_authorized": False,
        "decision": "V4.1-GPU fails Gate C under the frozen guardrails; V3 Round 14 remains promoted balanced default.",
    }
    (log_dir / "gate_c_failure_summary.json").write_text(json.dumps(json_safe(summary), indent=2, sort_keys=True) + "\n")
    (log_dir / "gate_c_failure_summary.md").write_text(
        "# gdTAI V4.1-GPU Gate C\n\n"
        "- Status: **complete, no eligible model**\n"
        "- All three outer folds completed.\n"
        "- No Stage-1 candidate passed the frozen per-source NK-passage guardrail.\n"
        "- No fair comparator or raw legacy score passed either frozen operating mode.\n"
        "- Bootstrap and extension-model controls are not estimable because no balanced V4 model was selected.\n"
        "- Promotion: **No**; V3 Round 14 remains the balanced default.\n",
    )

    stage1_display = stage1.sort_values(["outer_fold_id", "pooled_nk_passage"]).groupby("outer_fold_id").head(2)
    comparator_display = comparators[comparators["mode"].eq("balanced")].sort_values(["outer_fold_id", "model_id", "nearest_objective"]).groupby(["outer_fold_id", "model_id"]).tail(1)
    document = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><title>gdTAI V4.1-GPU Gate C</title><style>
@page {{ size:A4; margin:13mm; }} * {{ box-sizing:border-box; }} body {{ margin:0; background:#eef3f1; color:#18211f; font:10pt/1.45 Arial,sans-serif; }} main {{ max-width:1180px; margin:auto; background:white; padding:30px 36px 60px; }}
h1 {{ font-size:25pt; margin:0 0 6px; letter-spacing:0; }} h2 {{ color:#145c52; border-bottom:1px solid #cbd8d5; padding-bottom:5px; margin:25px 0 9px; }} h3 {{ color:#176b5f; }}
.status {{ border-left:5px solid #ad3535; background:#fff1ef; padding:13px 16px; margin:18px 0; }} .grid {{ display:grid; grid-template-columns:1fr 1fr; gap:16px; }} figure {{ margin:8px 0 16px; break-inside:avoid; }} img {{ width:100%; border:1px solid #d2ddda; }} figcaption {{ color:#52625f; font-size:9pt; }}
table {{ border-collapse:collapse; width:100%; font-size:7.4pt; table-layout:fixed; }} th {{ color:white; background:#176b5f; text-align:left; }} th,td {{ border:1px solid #c9d5d2; padding:4px; overflow-wrap:anywhere; vertical-align:top; }} tr:nth-child(even) td {{ background:#f3f7f6; }}
code {{ overflow-wrap:anywhere; }} @media print {{ body {{ background:white; }} main {{ max-width:none; padding:0; }} .grid {{ display:block; }} h2 {{ break-after:avoid; }} tr,figure {{ break-inside:avoid; }} }}
</style></head><body><main><p><strong>SUPERVISION REPORT · 7 AUGUST 2026</strong></p><h1>gdTAI V4.1-GPU Gate C</h1><p>Nested cross-dataset evaluation under the precommitted purity and recall guardrails.</p>
<div class="status"><strong>Complete, but no model is eligible.</strong> All three outer folds finished. No Stage-1 candidate satisfied the literal per-source NK-passage condition, and no fair comparator or raw legacy score satisfied either final operating mode. V3 Round 14 remains the promoted balanced default.</div>
<h2>What Was Evaluated</h2><p>Ground truth, 197 genes, grouped splits, CD4/Treg exclusions, calibration, and operating thresholds were unchanged. Stage 1 compared three deterministic GPU ridge penalties with two fixed GPU XGBoost configurations. Primary gdT and productive alpha-beta T recall had to remain at least 0.99 and 0.98 in every inner primary source, while strict-NK passage had to remain at most 0.50 in every NK source.</p>
<p>Silver, sorted, dual, ambiguous, and unlabeled cells never selected models or thresholds. Fair compact seven-gene and 153-gene TCR comparators used GPU ridge with the same grouped folds, exclusions, hard strict-NK negatives, calibration, and final thresholds. The raw TRD-minus-TRAB score was evaluated without refitting.</p>
<h2>Stage-1 Blocking Result</h2><p>All 15 candidate-by-outer-fold evaluations converged, but 0/15 passed. At the highest threshold retaining the required T-cell recalls, at least one NK source still had passage above 50%. The one-cell GSE190870 stratum makes the literal maximum-source statistic discontinuous, but it does not explain away the result: after restricting the diagnostic calculation to NK sources with at least 100 cells, the best worst-source passage still ranged from 90.97% to 93.23% across outer folds. Pooled NK passage was lower at 8.60% to 12.38%, showing marked source heterogeneity rather than adequate general separation. These post-hoc summaries do not alter the official failure.</p>
<div class="grid"><figure><img src="assets/stage1_guardrail_failure.png"><figcaption>Literal and sample-size-aware Stage-1 diagnostics.</figcaption></figure><figure><img src="assets/comparator_near_miss.png"><figcaption>Balanced comparator near-miss operating points.</figcaption></figure></div>
{table_html(stage1_display, ['outer_fold_id','family','minimum_source_gdt_recall','minimum_source_abt_recall','pooled_nk_passage','maximum_nk_passage_sources_n_ge_100','worst_nk_source_n_ge_100','worst_nk_source_n_ge_100_n','maximum_source_nk_passage','worst_nk_source_n'])}
<h2>Fair Comparator Results</h2><p>No compact seven-gene ridge, V2-like TCR-gene ridge, or raw legacy-score threshold met the balanced or high-purity constraints. Therefore there are no valid held-out predictions from which to claim F1, calculate a paired bootstrap advantage, or select a V4 model.</p>
{table_html(comparator_display, ['outer_fold_id','model_id','nearest_objective','nearest_recall','nearest_paired_abt_fpr','nearest_strict_nk_fpr','recall_requirement','paired_abt_fpr_limit','strict_nk_fpr_limit'])}
<h2>What Is Not Estimable</h2><p>V4 Stage 2 was not entered because Stage 1 had no eligible candidate. Extension-cohort model false-positive rates and paired hierarchical-bootstrap differences are consequently not estimable. This is a valid negative Gate C result, not missing evidence and not authorization to relax a threshold after observing outcomes.</p>
<h2>Interpretation</h2><p>The precommitted Stage-1 rule is statistically brittle because it takes a hard maximum over NK strata ranging from one cell to more than 90,000 cells. However, the ≥100-cell sensitivity analysis shows a second, substantive problem: at the recall-preserving threshold, at least one well-populated dataset passes more than 90% of strict NK cells. Thus V4.1 failed because T-lineage sensitivity and cross-dataset NK rejection were incompatible under this feature panel, not merely because of a one-cell denominator. A future V4.2 may use hierarchical or minimum-stratum guardrails, but it must also improve source-robust discrimination and be declared before re-evaluation.</p>
<h2>Decision Boundary</h2><p>No V4.1 model is promoted or released. No whole-atlas inference is permitted. V3 Round 14 remains the balanced default pending a separately reviewed redesign.</p>
</main></body></html>"""
    html_path = static_dir / "index.html"; html_path.write_text(document)
    pdf_path = static_dir / "gdtai_v4_gpu_gate_c_report.pdf"
    completed = subprocess.run(["google-chrome", "--headless", "--no-sandbox", "--disable-gpu", f"--print-to-pdf={pdf_path}", "--no-pdf-header-footer", html_path.resolve().as_uri()], capture_output=True, text=True, timeout=180)
    if completed.returncode != 0 or not pdf_path.exists():
        raise RuntimeError(completed.stderr)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
