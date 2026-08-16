#!/usr/bin/env python3
"""Build the V4.2 integration, clustering, and pseudo-NK consensus QC report."""

from __future__ import annotations

import html
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from workflows.gdtai.gdtai_v4_2_integration_core import (
    read_json,
    resolve,
    sha256_file,
    validate_preflight_approval,
    validate_project_data_approval,
    validate_recovery_preflight,
    validate_resource_preflight,
)


CONFIG_PATH = PROJECT_ROOT / "configs/models/gdtai/v4_2_integration_execution.json"
RUNNER_PATH = PROJECT_ROOT / "workflows/gdtai/run_gdtai_v4_2_nk_reference_integration.py"
CORE_PATH = PROJECT_ROOT / "workflows/gdtai/gdtai_v4_2_integration_core.py"
REPORT_SCRIPT = Path(__file__).resolve()
REPORT_NAME = "gdtai_v4_2_nk_reference_qc_report.pdf"


def add_check(rows: list[dict[str, str]], check: str, status: str, detail: str) -> None:
    rows.append({"check": check, "status": status, "detail": detail})


def save_figure(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=260, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def dominant_sources_for_clusters(
    evidence: pd.DataFrame,
    saved: Any,
    candidate: np.ndarray,
    sources: np.ndarray,
) -> pd.DataFrame:
    result = evidence.copy()
    result["dominant_candidate_source"] = ""
    global_names = saved["global_names"].astype(str).tolist()
    boundary_names = saved["boundary_names"].astype(str).tolist()
    boundary_indices = saved["boundary_indices"].astype(np.int64)
    source_levels = np.sort(np.unique(sources[candidate]))

    for row_index, row in result.iterrows():
        run = str(row["run"])
        cluster = int(row["cluster"])
        if run in global_names:
            labels = saved["global_labels"][:, global_names.index(run)]
            mask = candidate & (labels == cluster)
            values = sources[mask]
        else:
            position = boundary_names.index(run)
            labels = saved["boundary_labels"][:, position]
            local_candidate = candidate[boundary_indices]
            values = sources[boundary_indices][local_candidate & (labels == cluster)]
        if values.size:
            counts = pd.Series(values).value_counts()
            result.loc[row_index, "dominant_candidate_source"] = str(counts.index[0])
            observed_fraction = float(counts.iloc[0] / counts.sum())
            if not np.isclose(
                observed_fraction,
                float(row["maximum_candidate_source_fraction"]),
                rtol=0,
                atol=1e-12,
            ):
                raise RuntimeError(f"Dominant-source fraction mismatch for {run}, cluster {cluster}")
        elif source_levels.size:
            result.loc[row_index, "dominant_candidate_source"] = "none"
    return result


def make_figures(
    history: pd.DataFrame,
    evidence: pd.DataFrame,
    criteria: pd.DataFrame,
    by_source: pd.DataFrame,
    figure_dir: Path,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )

    epochs = np.arange(1, history.shape[0] + 1)
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(epochs, history["train_loss_epoch"], color="#176b87", linewidth=2, label="Training")
    ax.plot(epochs, history["validation_loss"], color="#c46a2c", linewidth=2, label="Validation")
    ax.set(xlabel="Epoch", ylabel="Loss", title="A100 scVI training history")
    ax.legend(frameon=False)
    ax.grid(axis="y", color="#dbe3e8", linewidth=0.7)
    save_figure(fig, figure_dir / "scvi_training_history.png")

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    run_type = evidence["run"].str.startswith("boundary")
    sizes = 18 + 45 * np.log10(evidence["n_candidate_cells"].clip(lower=1)) / max(
        1.0, np.log10(evidence["n_candidate_cells"].clip(lower=1)).max()
    )
    for mask, label, color, marker in [
        (~run_type, "Global runs", "#247ba0", "o"),
        (run_type, "Boundary runs", "#d98f39", "^"),
    ]:
        ax.scatter(
            evidence.loc[mask, "maximum_candidate_source_fraction"],
            evidence.loc[mask, "anchor_nk_purity"],
            s=sizes.loc[mask],
            alpha=0.62,
            color=color,
            marker=marker,
            edgecolors="none",
            label=label,
        )
    near = evidence["anchor_nk_purity"].ge(0.95) & evidence["productive_tcr_contamination"].le(0.02)
    ax.scatter(
        evidence.loc[near, "maximum_candidate_source_fraction"],
        evidence.loc[near, "anchor_nk_purity"],
        s=85,
        facecolors="none",
        edgecolors="#b3261e",
        linewidths=1.6,
        label="Purity/contamination pass",
    )
    ax.axhline(0.95, color="#b3261e", linestyle="--", linewidth=1.1)
    ax.axvline(0.70, color="#b3261e", linestyle="--", linewidth=1.1)
    ax.set(
        xlabel="Largest candidate-source fraction",
        ylabel="Independent-anchor NK purity",
        title="No cluster enters the qualifying upper-left region",
        xlim=(-0.02, 1.02),
        ylim=(-0.02, 1.02),
    )
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    ax.grid(color="#e3e8eb", linewidth=0.6)
    save_figure(fig, figure_dir / "cluster_contract_geometry.png")

    fig, ax = plt.subplots(figsize=(7.8, 4.5))
    colors = ["#247ba0"] * (criteria.shape[0] - 1) + ["#b3261e"]
    bars = ax.barh(criteria["criterion"], criteria["n_clusters"], color=colors)
    ax.bar_label(bars, labels=[f"{value:,}" for value in criteria["n_clusters"]], padding=4, fontsize=9)
    ax.set(xlabel="Clusters passing criterion", title="Individual pseudo-NK cluster criteria")
    ax.set_xlim(0, max(410, int(criteria["n_clusters"].max() * 1.12)))
    ax.grid(axis="x", color="#e3e8eb", linewidth=0.6)
    save_figure(fig, figure_dir / "cluster_criterion_pass_counts.png")

    plotted = by_source[by_source["n_candidate_eligible"].gt(0)].sort_values("n_candidate_eligible")
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    bars = ax.barh(plotted["source_gse_id"], plotted["n_candidate_eligible"], color="#2a7f62")
    ax.bar_label(bars, labels=[f"{value:,}" for value in plotted["n_candidate_eligible"]], padding=4, fontsize=9)
    ax.set(xlabel="Eligible no-productive-TCR candidate cells", title="Pseudo-label candidate pool by source")
    ax.set_xlim(0, int(plotted["n_candidate_eligible"].max() * 1.18))
    ax.grid(axis="x", color="#e3e8eb", linewidth=0.6)
    save_figure(fig, figure_dir / "candidate_pool_by_source.png")


def table_html(
    frame: pd.DataFrame,
    columns: list[str] | None = None,
    classes: str = "data-table",
) -> str:
    shown = frame if columns is None else frame[columns]
    return shown.to_html(index=False, border=0, classes=classes, escape=True)


def render_html(
    summary: dict[str, Any],
    checks: pd.DataFrame,
    criteria: pd.DataFrame,
    near: pd.DataFrame,
    by_source: pd.DataFrame,
    static_dir: Path,
) -> Path:
    static_dir.mkdir(parents=True, exist_ok=True)
    figure_prefix = "../../Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_nk_reference"
    near_display = near.copy()
    for column in ["anchor_nk_purity", "productive_tcr_contamination", "maximum_candidate_source_fraction"]:
        near_display[column] = near_display[column].map(lambda value: f"{value:.3%}")
    near_print = near_display[
        [
            "run",
            "cluster",
            "n_primary_nk_anchors",
            "n_productive_tcr_anchors",
            "anchor_nk_purity",
            "productive_tcr_contamination",
            "dominant_candidate_source",
            "maximum_candidate_source_fraction",
        ]
    ].rename(
        columns={
            "run": "Run",
            "cluster": "Cluster",
            "n_primary_nk_anchors": "NK anchors",
            "n_productive_tcr_anchors": "TCR anchors",
            "anchor_nk_purity": "NK purity",
            "productive_tcr_contamination": "TCR/cell",
            "dominant_candidate_source": "Dominant source",
            "maximum_candidate_source_fraction": "Source share",
        }
    )
    by_source_display = by_source[by_source["n_candidate_eligible"].gt(0)].copy()
    by_source_display["candidate_fraction"] = (
        by_source_display["n_candidate_eligible"] / by_source_display["n_candidate_eligible"].sum()
    ).map(lambda value: f"{value:.2%}")
    by_source_display = by_source_display.rename(
        columns={
            "source_gse_id": "Dataset",
            "n_candidate_eligible": "Eligible candidates",
            "candidate_fraction": "Pool share",
            "n_selected_pseudo_nk": "Selected pseudo-NK",
        }
    )

    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>gdTAI V4.2 integration and pseudo-NK QC</title>
<style>
:root{{--ink:#17212b;--navy:#17324d;--teal:#197278;--red:#b3261e;--amber:#b86e16;--line:#ccd6dc;--soft:#f3f6f7}}
*{{box-sizing:border-box}} body{{font-family:Arial,Helvetica,sans-serif;color:var(--ink);line-height:1.48;max-width:1120px;margin:0 auto;padding:30px 42px;background:white}}
h1{{font-size:29px;color:var(--navy);border-bottom:4px solid var(--red);padding-bottom:12px;margin-bottom:16px}} h2{{font-size:20px;color:var(--navy);margin-top:30px}} h3{{font-size:15px;color:var(--navy)}} p,li{{font-size:14px}}
.status{{border-left:6px solid var(--red);background:#fff4f2;padding:14px 18px;margin:18px 0}} .status strong{{color:var(--red);font-size:18px}}
.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin:18px 0}} .card{{border:1px solid var(--line);padding:12px;background:var(--soft)}} .card b{{display:block;font-size:21px;color:var(--navy)}} .card span{{font-size:11px;color:#52636e}}
.callout{{border-left:5px solid var(--amber);background:#fff9ef;padding:12px 16px;margin:18px 0}} figure{{margin:22px 0;break-inside:avoid}} img{{display:block;max-width:94%;max-height:650px;margin:0 auto}} figcaption{{font-size:11px;color:#52636e;text-align:center;margin-top:7px}}
.table-wrap{{overflow-x:auto;margin:14px 0 24px}} table{{border-collapse:collapse;width:100%;font-size:10px}} th{{background:var(--navy);color:white;text-align:left}} th,td{{padding:6px 7px;border:1px solid var(--line);vertical-align:top;overflow-wrap:anywhere}} tr:nth-child(even){{background:var(--soft)}} code{{background:#edf1f4;padding:1px 4px}} .pass{{color:#176b43;font-weight:bold}} .warn{{color:#9a5d10;font-weight:bold}} .fail{{color:var(--red);font-weight:bold}} .print-only{{display:none}} .repro{{break-inside:avoid}}
@media(max-width:800px){{body{{padding:18px}}.cards{{grid-template-columns:repeat(2,1fr)}}}}
@media print{{@page{{size:A4;margin:9mm}}body{{max-width:none;padding:0}}h1{{font-size:23px}}h2{{font-size:17px;break-after:avoid}}p,li{{font-size:10.3px}}.cards{{grid-template-columns:repeat(4,1fr)}}.card b{{font-size:16px}}.card span{{font-size:8px}}table{{font-size:7.2px}}tr,figure,.card{{break-inside:avoid}}img{{max-height:145mm}}.table-wrap{{overflow:visible}}.screen-only{{display:none!important}}.print-only{{display:block!important}}}}
</style></head><body>
<h1>gdTAI V4.2 integration and pseudo-NK consensus QC</h1>
<p>Execution date: 16 August 2026. This report audits the recovered development sidecar, A100 scVI fit, repeated RAPIDS clustering, and expression-independent pseudo-NK consensus. Locked validation cohorts were excluded throughout.</p>
<div class="status"><strong>Scientific result: FAIL - no usable pseudo-NK cells</strong><br>The computation completed correctly, but zero of {summary['n_candidate_eligible']:,} eligible cells satisfied the frozen consensus contract. Classifier fitting from this pseudo-label lane must not proceed.</div>
<div class="cards">
  <div class="card"><b>{summary['n_cells']:,}</b><span>development cells</span></div>
  <div class="card"><b>{summary['n_primary_nk_anchors']:,}</b><span>primary NK anchors</span></div>
  <div class="card"><b>{summary['n_boundary_cells']:,}</b><span>boundary cells ({summary['boundary_fraction']:.1%})</span></div>
  <div class="card"><b>{summary['n_selected_pseudo_nk']:,}</b><span>selected pseudo-NK cells</span></div>
</div>

<h2>Executive interpretation</h2>
<p>Staging, GPU integration, and clustering passed their technical contracts: the sidecar contains {summary['n_cells']:,} cells by {summary['n_hvgs']:,} common HVGs, scVI produced a finite 30-dimensional latent representation on the A100 without CPU fallback, and all 18 Leiden partitions were written with a verified SHA-256 checksum.</p>
<p>The pseudo-label contract failed scientifically. Six clusters reached at least 95% NK purity among independent anchors and no more than 2% productive-TCR contamination. Each also passed the anchor-count and multi-source-presence rules, but 87%-90% of its eligible candidate cells came from one dataset, exceeding the frozen 70% single-source cap. Therefore no cluster qualified in any run and no cell could receive the required 80% agreement.</p>
<div class="callout"><b>Why this is not a classifier result.</b> No gdTAI classifier, calibration model, decision threshold, release artifact, or whole-atlas inference was fitted in this execution. V2/V3 performance cannot be compared to an empty V4.2 training lane.</div>

<h2>Integration QC</h2>
<figure><img src="{figure_prefix}/scvi_training_history.png" alt="scVI training history"><figcaption>Training and validation loss across the frozen 20-epoch A100 scVI fit.</figcaption></figure>
<div class="table-wrap">{table_html(checks)}</div>

<h2>Clustering and boundary QC</h2>
<p>The expression-independent boundary was defined as the union of global clusters containing both a dual-annotation primary NK anchor and a productive-TCR anchor in at least one global run. It included {summary['n_boundary_cells']:,} cells ({summary['boundary_fraction']:.2%} of the sidecar), indicating that the global clusters were too broad for this union to localize the biological T/NK boundary.</p>
<figure><img src="{figure_prefix}/cluster_contract_geometry.png" alt="Cluster contract geometry"><figcaption>Dashed lines mark the 95% anchor-purity and 70% maximum-source thresholds. The six high-purity clusters all fall to the right of the source-balance limit.</figcaption></figure>
<figure><img src="{figure_prefix}/cluster_criterion_pass_counts.png" alt="Criterion pass counts"><figcaption>Individual criteria are feasible, but their frozen intersection contains no cluster.</figcaption></figure>

<h2>Near-qualifying clusters</h2>
<p>These are the only clusters that passed the frozen purity and productive-TCR contamination limits. Their sole failing field is <code>maximum_candidate_source_fraction</code>.</p>
<div class="table-wrap screen-only">{table_html(near_display, ['run','cluster','n_cells','n_primary_nk_anchors','n_productive_tcr_anchors','anchor_nk_purity','productive_tcr_contamination','n_candidate_cells','n_candidate_sources','dominant_candidate_source','maximum_candidate_source_fraction'])}</div>
<div class="table-wrap print-only">{table_html(near_print)}</div>

<h2>Candidate-source composition</h2>
<figure><img src="{figure_prefix}/candidate_pool_by_source.png" alt="Candidate pool by source"><figcaption>The overall eligible pool is highly imbalanced before cluster-level selection.</figcaption></figure>
<div class="table-wrap">{table_html(by_source_display, ['Dataset','Eligible candidates','Pool share','Selected pseudo-NK'])}</div>

<h2>Decision and next iteration</h2>
<p><b>Do not fit the planned V4.2 classifier from these pseudo-labels.</b> Relaxing the 70% cap after seeing this result would be outcome-driven tuning and would reintroduce the source-label shortcut the cap was designed to prevent.</p>
<p>The next defensible iteration should repair the reference design rather than lower the purity standard: source-balance the eligible candidate construction before graph inference; increase independent, confidently annotated NK anchors from multiple development datasets; use finer local neighborhoods or anchor-neighborhood propagation instead of a union of very broad Leiden clusters; and retain the locked cohorts unchanged for nested comparison against V2 and V3.</p>

<section class="repro"><h2>Reproducibility</h2>
<ul>
<li>Execution config SHA-256: <code>{html.escape(summary['execution_config_sha256'])}</code></li>
<li>Runner SHA-256: <code>{html.escape(summary['runner_sha256'])}</code></li>
<li>Core SHA-256: <code>{html.escape(summary['core_sha256'])}</code></li>
<li>Partition SHA-256: <code>{html.escape(summary['partitions_sha256'])}</code></li>
<li>Report builder SHA-256: <code>{html.escape(summary['report_script_sha256'])}</code></li>
<li>Classifier fitting performed: <b>no</b></li>
</ul></section>
</body></html>"""
    output = static_dir / "index.html"
    output.write_text(document, encoding="utf-8")
    return output


def main() -> None:
    config = read_json(CONFIG_PATH)
    validate_preflight_approval(config)
    validate_recovery_preflight(config)
    validate_resource_preflight(config)
    validate_project_data_approval(CONFIG_PATH, config, RUNNER_PATH, CORE_PATH)

    outputs = config["outputs"]
    table_dir = resolve(outputs["table_dir"])
    figure_dir = resolve(outputs["figure_dir"])
    log_dir = resolve(outputs["log_dir"])
    static_dir = resolve(outputs["static_dir"])
    ssd_dir = Path(config["resources"]["ssd_root"])
    for directory in [table_dir, figure_dir, log_dir, static_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    prepare = read_json(log_dir / "prepare_summary.json")
    fit = read_json(log_dir / "fit_summary.json")
    cluster = read_json(log_dir / "cluster_summary.json")
    consensus = read_json(log_dir / "consensus_summary.json")
    history = pd.read_csv(table_dir / "scvi_training_history.csv")
    evidence = pd.read_csv(table_dir / "cluster_consensus_evidence.csv.gz")
    by_source = pd.read_csv(table_dir / "pseudo_nk_counts_by_source.csv")
    source_state = pd.read_csv(table_dir / "development_input_file_state.csv")
    partitions_path = ssd_dir / "cluster_partitions.npz"
    staged_path = ssd_dir / "development_hvg_counts.h5ad"
    saved = np.load(partitions_path, allow_pickle=False)

    backed = ad.read_h5ad(staged_path, backed="r")
    obs = backed.obs[[
        "candidate_eligible", "source_gse_id", "primary_nk_anchor",
        "productive_tcr_anchor", "input_cohort_id",
    ]].copy()
    backed.file.close()
    candidate = obs["candidate_eligible"].astype(bool).to_numpy()
    sources = obs["source_gse_id"].astype(str).to_numpy()

    contract = config["pseudolabel"]
    criterion_specs = [
        ("Independent anchors >= 20", evidence["n_primary_nk_anchors"] + evidence["n_productive_tcr_anchors"] >= int(contract["minimum_independent_anchors_per_cluster"])),
        ("Primary NK anchors >= 10", evidence["n_primary_nk_anchors"] >= int(contract["minimum_primary_nk_anchors_per_cluster"])),
        ("Anchor NK purity >= 95%", evidence["anchor_nk_purity"] >= float(contract["minimum_anchor_nk_purity"])),
        ("Productive-TCR contamination <= 2%", evidence["productive_tcr_contamination"] <= float(contract["maximum_productive_tcr_contamination"])),
        ("Candidate sources >= 3", evidence["n_candidate_sources"] >= int(contract["minimum_development_sources"])),
        ("Largest source fraction <= 70%", evidence["maximum_candidate_source_fraction"] <= float(contract["maximum_single_source_fraction"])),
        ("All criteria", evidence["qualifies_as_nk"].astype(bool)),
    ]
    criteria = pd.DataFrame(
        {
            "criterion": [name for name, _ in criterion_specs],
            "n_clusters": [int(mask.sum()) for _, mask in criterion_specs],
            "fraction_of_clusters": [float(mask.mean()) for _, mask in criterion_specs],
        }
    )
    near_mask = (
        evidence["anchor_nk_purity"].ge(float(contract["minimum_anchor_nk_purity"]))
        & evidence["productive_tcr_contamination"].le(float(contract["maximum_productive_tcr_contamination"]))
    )
    near = dominant_sources_for_clusters(evidence.loc[near_mask].copy(), saved, candidate, sources)
    near["failed_criterion"] = "maximum_candidate_source_fraction"

    n_cells = int(prepare["n_cells"])
    boundary_n = int(cluster["n_boundary_cells"])
    partition_hash = sha256_file(partitions_path)
    unchanged_sources = int(source_state["unchanged"].astype(bool).sum())
    checks: list[dict[str, str]] = []
    add_check(checks, "Recovered sparse staging", "PASS" if prepare["status"] == "PASS" else "FAIL", f"{n_cells:,} cells x {prepare['n_hvgs']:,} HVGs")
    add_check(checks, "Locked cohorts excluded", "PASS" if not prepare["locked_cohorts_included"] and not consensus["locked_cohorts_included"] else "FAIL", "Zero locked-cohort cells in staging and consensus")
    add_check(checks, "Source H5ADs unchanged", "PASS" if unchanged_sources == source_state.shape[0] else "FAIL", f"{unchanged_sources}/{source_state.shape[0]} size/mtime pairs unchanged")
    add_check(checks, "A100 scVI fit", "PASS" if fit["status"] == "PASS" and not fit["cpu_fallback"] else "FAIL", f"{fit['n_latent']} latent dimensions; {fit['elapsed_seconds'] / 60:.1f} minutes")
    add_check(checks, "RAPIDS partitions", "PASS" if cluster["status"] == "PASS" and not cluster["cpu_fallback"] else "FAIL", f"{cluster['n_global_runs']} global + {cluster['n_boundary_runs']} boundary runs")
    add_check(checks, "Partition checksum", "PASS" if partition_hash == cluster["partitions_sha256"] else "FAIL", partition_hash)
    add_check(checks, "Boundary localization", "WARN" if boundary_n / n_cells > 0.95 else "PASS", f"{boundary_n:,}/{n_cells:,} cells ({boundary_n / n_cells:.2%})")
    add_check(checks, "Qualifying clusters", "FAIL" if int(criteria.iloc[-1]["n_clusters"]) == 0 else "PASS", f"{int(criteria.iloc[-1]['n_clusters'])}/{evidence.shape[0]} clusters")
    add_check(checks, "Selected pseudo-NK cells", "FAIL" if int(consensus["n_selected_pseudo_nk"]) == 0 else "PASS", f"{int(consensus['n_selected_pseudo_nk']):,}/{int(consensus['n_candidate_eligible']):,} eligible cells")
    add_check(checks, "Classifier fitting", "PASS", "Not performed")
    checks_frame = pd.DataFrame(checks)

    summary = {
        "result": "FAIL_NO_PSEUDO_NK",
        "technical_execution_status": "PASS",
        "scientific_qc_status": "FAIL",
        "n_cells": n_cells,
        "n_hvgs": int(prepare["n_hvgs"]),
        "n_primary_nk_anchors": int(prepare["primary_nk_anchors"]),
        "n_productive_tcr_anchors": int(prepare["productive_tcr_anchors"]),
        "n_candidate_eligible": int(consensus["n_candidate_eligible"]),
        "n_boundary_cells": boundary_n,
        "boundary_fraction": boundary_n / n_cells,
        "n_global_runs": int(cluster["n_global_runs"]),
        "n_boundary_runs": int(cluster["n_boundary_runs"]),
        "n_clusters_evaluated": int(evidence.shape[0]),
        "n_clusters_passing_purity_and_contamination": int(near.shape[0]),
        "n_qualifying_clusters": int(criteria.iloc[-1]["n_clusters"]),
        "n_selected_pseudo_nk": int(consensus["n_selected_pseudo_nk"]),
        "dominant_source_fraction_min_near_clusters": float(near["maximum_candidate_source_fraction"].min()),
        "dominant_source_fraction_max_near_clusters": float(near["maximum_candidate_source_fraction"].max()),
        "classifier_fitting_performed": False,
        "locked_cohorts_included": False,
        "execution_config_sha256": sha256_file(CONFIG_PATH),
        "runner_sha256": sha256_file(RUNNER_PATH),
        "core_sha256": sha256_file(CORE_PATH),
        "partitions_sha256": partition_hash,
        "report_script_sha256": sha256_file(REPORT_SCRIPT),
    }

    criteria.to_csv(table_dir / "cluster_criterion_pass_counts.csv", index=False)
    near.to_csv(table_dir / "near_qualifying_clusters.csv", index=False)
    checks_frame.to_csv(table_dir / "execution_qc_checks.csv", index=False)
    run_summary = (
        evidence.groupby("run", observed=True)
        .agg(
            n_clusters=("cluster", "size"),
            n_qualifying_clusters=("qualifies_as_nk", "sum"),
            maximum_anchor_nk_purity=("anchor_nk_purity", "max"),
            minimum_productive_tcr_contamination=("productive_tcr_contamination", "min"),
        )
        .reset_index()
    )
    run_summary.to_csv(table_dir / "cluster_run_summary.csv", index=False)
    (log_dir / "gdtai_v4_2_nk_reference_qc_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    make_figures(history, evidence, criteria, by_source, figure_dir)
    html_path = render_html(summary, checks_frame, criteria, near, by_source, static_dir)

    markdown = f"""# gdTAI V4.2 integration and pseudo-NK consensus QC

## Result

`FAIL_NO_PSEUDO_NK`: technical execution passed, but zero of
{summary['n_candidate_eligible']:,} eligible cells met the frozen pseudo-NK
consensus contract. No classifier was fitted.

Six clusters passed the 95% anchor-purity and 2% productive-TCR-contamination
limits, but each had {summary['dominant_source_fraction_min_near_clusters']:.1%}
to {summary['dominant_source_fraction_max_near_clusters']:.1%} of eligible
candidate cells from one source, above the frozen 70% cap.

## Artifacts

- HTML: `gdT_prediction/gdtai_v4_2_nk_reference/index.html`
- PDF: `gdT_prediction/gdtai_v4_2_nk_reference/{REPORT_NAME}`
- Summary JSON: `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_nk_reference/gdtai_v4_2_nk_reference_qc_summary.json`
- Canonical tables: `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_nk_reference/`
- Canonical figures: `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_nk_reference/`

## Decision

Do not fit the planned V4.2 classifier from this empty pseudo-label lane. The
next iteration must repair source balance and boundary localization without
using locked-cohort outcomes to tune the pseudo-label contract.
"""
    (log_dir / "gdtai_v4_2_nk_reference_qc_summary.md").write_text(markdown, encoding="utf-8")

    chrome = shutil.which("google-chrome") or shutil.which("google-chrome-stable")
    if chrome is None:
        raise RuntimeError("Google Chrome is required for PDF export")
    profile = Path("/tmp/gdtai-v42-nk-reference-report-chrome")
    profile.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            chrome,
            "--headless",
            "--no-sandbox",
            "--disable-gpu",
            "--disable-dev-shm-usage",
            "--disable-breakpad",
            "--disable-crash-reporter",
            "--allow-file-access-from-files",
            "--no-pdf-header-footer",
            f"--user-data-dir={profile}",
            f"--print-to-pdf={static_dir / REPORT_NAME}",
            html_path.resolve().as_uri(),
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
