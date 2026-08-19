#!/usr/bin/env python3
"""Render the gdTAI V4.2 nested-development gate report."""

from __future__ import annotations

import html
import json
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xgboost as xgb

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from workflows.gdtai.train_gdtai_v4_2_nested import exclusion_flags


TABLE = ROOT / "Integrated_dataset/tables/gdT_prediction"
LOG = ROOT / "Integrated_dataset/logs/gdT_prediction"
FIGURE = ROOT / "Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_development_gate"
STATIC = ROOT / "gdT_prediction/gdtai_v4_2_development_gate"
MODEL = ROOT / "Integrated_dataset/models/gdT_prediction/gdtai_v4_2_training"


def pct(value: float) -> str:
    return "NA" if pd.isna(value) else f"{100 * value:.2f}%"


def table_html(frame: pd.DataFrame, formats: dict[str, callable] | None = None) -> str:
    local = frame.copy()
    for column, formatter in (formats or {}).items():
        if column in local:
            local[column] = local[column].map(formatter)
    return local.to_html(index=False, border=0, classes="data", escape=True)


def stage2_eligibility(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    rows = []
    for _, row in frame.iterrows():
        payload = json.loads(row.profiles_json)
        for profile in ("balanced", "high_purity"):
            rows.append({
                "heldout_source": row.heldout_source,
                "profile": profile,
                "eligible": payload[profile]["metrics"].get("eligible", False) is True,
            })
    return pd.DataFrame(rows).groupby(["heldout_source", "profile"], as_index=False).eligible.sum().rename(columns={"eligible": "eligible_candidates"})


def exclusion_cost() -> pd.DataFrame:
    labels = pd.read_parquet(TABLE / "gdtai_v4_2_ground_truth/v4_2_label_manifest.parquet")
    features = pd.read_csv(TABLE / "gdtai_v4_2_training/feature_manifest.csv").sort_values("feature_index")
    matrix = np.load(ROOT / "Integrated_dataset/cache/gdT_prediction/gdtai_v4_2_training/gene_features.npy", mmap_mode="r")
    rows = np.flatnonzero(labels.truth_class.eq("gdT_gold").to_numpy())
    cd4, treg, excluded = exclusion_flags(np.asarray(matrix[rows]), features.gene.astype(str).tolist())
    audit = labels.iloc[rows][["source_gse_id"]].copy()
    audit["cd4_rule"] = cd4
    audit["treg_rule"] = treg
    audit["excluded"] = excluded
    result = audit.groupby("source_gse_id", as_index=False).agg(
        gold_gdT=("excluded", "size"), cd4_excluded=("cd4_rule", "sum"),
        treg_excluded=("treg_rule", "sum"), total_excluded=("excluded", "sum"),
        recall_cost=("excluded", "mean"),
    )
    result.to_csv(TABLE / "gdtai_v4_2_training/cd4_treg_exclusion_cost.csv", index=False)
    return result


def feature_gain() -> pd.DataFrame:
    names = pd.read_csv(TABLE / "gdtai_v4_2_training/feature_manifest.csv").sort_values("feature_index").gene.astype(str).tolist()
    rows = []
    for heldout in ("HRA005041", "MalteGDT"):
        contract = json.loads((MODEL / f"outer_{heldout}/contract.json").read_text())
        for profile, profile_spec in contract.get("profiles", {}).items():
            path = MODEL / f"outer_{heldout}/stage2_{profile_spec['stage2_candidate_id']}.ubj"
            booster = xgb.Booster()
            booster.load_model(path)
            for key, gain in booster.get_score(importance_type="gain").items():
                index = int(key[1:])
                rows.append({
                    "heldout_source": heldout,
                    "profile": profile,
                    "feature": names[index] if index < len(names) else "stage1_probability",
                    "gain": gain,
                })
    frame = pd.DataFrame(rows)
    top = frame.sort_values(["heldout_source", "profile", "gain"], ascending=[True, True, False]).groupby(["heldout_source", "profile"]).head(12)
    top.to_csv(TABLE / "gdtai_v4_2_training/outer_model_top_gain.csv", index=False)
    return top


def make_figures(outer: pd.DataFrame, eligibility: pd.DataFrame, exclusions: pd.DataFrame,
                 nk: pd.DataFrame, recovery: pd.DataFrame) -> list[Path]:
    FIGURE.mkdir(parents=True, exist_ok=True)
    paths = []
    palette = {"primary": "#24577A", "recovery": "#B04A3A", "guard": "#444444"}

    balanced = outer[outer.profile.eq("balanced")].set_index("heldout_source")
    labels = ["HRA005041", "GSE144469", "MalteGDT", "Malte recovery"]
    values = [balanced.recall.get("HRA005041", np.nan), np.nan, balanced.recall.get("MalteGDT", np.nan), recovery.loc[recovery.profile.eq("balanced"), "recall"].iloc[0]]
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    bars = ax.bar(labels, np.nan_to_num(values), color=[palette["primary"], "#B8BEC3", palette["primary"], palette["recovery"]])
    ax.axhline(0.70, color=palette["guard"], linestyle="--", linewidth=1.4, label="Frozen per-source floor (70%)")
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, max(0.025, 0 if pd.isna(value) else value + 0.025), "No eligible\ninner profile" if pd.isna(value) else pct(value), ha="center", va="bottom", fontsize=9)
    ax.set_ylim(0, 1.08); ax.set_ylabel("Held-out gdT recall"); ax.set_title("Source-held-out transfer is the blocking failure")
    ax.legend(frameon=False, loc="lower left"); ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); path = FIGURE / "heldout_recall.png"; fig.savefig(path, dpi=240); plt.close(fig); paths.append(path)

    pivot = eligibility.pivot(index="heldout_source", columns="profile", values="eligible_candidates").fillna(0)
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    pivot[[c for c in ["balanced", "high_purity"] if c in pivot]].plot.bar(ax=ax, color=["#2A7F62", "#D08A24"])
    ax.set_ylabel("Eligible candidates out of 12"); ax.set_xlabel(""); ax.set_title("Frozen inner-profile eligibility")
    ax.legend(title="Profile", frameon=False); ax.spines[["top", "right"]].set_visible(False); ax.tick_params(axis="x", rotation=0)
    fig.tight_layout(); path = FIGURE / "candidate_eligibility.png"; fig.savefig(path, dpi=240); plt.close(fig); paths.append(path)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), gridspec_kw={"width_ratios": [1.15, 1]})
    ordered = exclusions.sort_values("recall_cost", ascending=False)
    axes[0].barh(ordered.source_gse_id, 100 * ordered.recall_cost, color="#5D7A8C")
    axes[0].set_xlabel("Gold gdT removed (%)"); axes[0].set_title("CD4/Treg exclusion recall cost")
    axes[0].spines[["top", "right"]].set_visible(False)
    nk2 = nk.sort_values("n_eligible")
    axes[1].barh(nk2.source_gse_id, nk2.n_eligible, color=np.where(nk2.tier.eq("tier1"), "#7B3F61", "#77966D"))
    axes[1].set_xlabel("Eligible NK negatives"); axes[1].set_title("Seven independent NK sources")
    axes[1].spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); path = FIGURE / "exclusions_and_nk_truth.png"; fig.savefig(path, dpi=240); plt.close(fig); paths.append(path)
    return paths


def main() -> None:
    STATIC.mkdir(parents=True, exist_ok=True)
    outer = pd.read_csv(TABLE / "gdtai_v4_2_training/nested_outer_metrics.csv")
    stage1 = pd.read_csv(TABLE / "gdtai_v4_2_training/nested_stage1_candidates.csv")
    eligibility = stage2_eligibility(TABLE / "gdtai_v4_2_training/nested_stage2_candidates.csv")
    recovery = pd.read_csv(TABLE / "gdtai_v4_2_recovery/nested_outer_metrics.csv")
    exclusions = exclusion_cost()
    gains = feature_gain()
    nk = pd.read_csv(TABLE / "gdtai_v4_2_ground_truth/nk_author_curation_audit.csv")
    ground = json.loads((LOG / "gdtai_v4_2_ground_truth/summary.json").read_text())
    nested = json.loads((LOG / "gdtai_v4_2_training/nested_summary.json").read_text())
    figures = make_figures(outer, eligibility, exclusions, nk, recovery)
    for path in figures:
        target = STATIC / path.name
        target.write_bytes(path.read_bytes())

    primary_display = outer[["heldout_source", "profile", "eligible_inner", "precision", "recall", "abt_fpr", "f1", "pr_auc"]].copy()
    recovery_display = recovery[["heldout_source", "profile", "eligible_inner", "stage2_candidate_id", "precision", "recall", "f1"]].copy()
    stage1_summary = stage1.groupby("heldout_source", as_index=False).agg(candidates=("candidate_id", "size"), eligible=("eligible", "sum"), threshold_min=("threshold", "min"), threshold_max=("threshold", "max"))
    historical = pd.concat([
        pd.read_csv(ROOT / "Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v2.0/mode_metrics.csv").assign(version="V2").rename(columns={"mode": "profile", "validation_precision": "precision", "validation_recall": "recall", "validation_f1": "f1"})[["version", "profile", "precision", "recall", "f1"]],
        pd.read_csv(ROOT / "Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/mode_metrics.csv").query("evaluation == 'independent_external_primary'").assign(version="V3 round14", profile="fixed 0.936")[["version", "profile", "precision", "recall", "f1"]],
    ], ignore_index=True)

    css = """
    @page { size: A4 landscape; margin: 11mm; }
    body { font-family: Arial, sans-serif; color:#1f2933; margin:0 auto; max-width:1180px; line-height:1.43; font-size:14px; }
    h1 { color:#173B57; font-size:30px; margin-bottom:4px; } h2 { color:#24577A; border-bottom:2px solid #D9E3EA; padding-bottom:5px; margin-top:23px; }
    h3 { color:#315D72; } .subtitle { color:#596B78; margin-top:0; } .decision { background:#FFF0EC; border-left:6px solid #B04A3A; padding:14px 18px; }
    .ok { background:#EDF7F2; border-left:5px solid #2A7F62; padding:10px 14px; } .grid { display:grid; grid-template-columns:repeat(4,1fr); gap:10px; }
    .metric { border:1px solid #CAD5DC; padding:10px; } .metric b { display:block; font-size:22px; color:#173B57; }
    table.data { border-collapse:collapse; width:100%; font-size:11px; table-layout:auto; } table.data th { background:#244F68; color:white; }
    table.data th, table.data td { padding:5px 7px; border:1px solid #D5DEE3; text-align:right; white-space:normal; overflow-wrap:anywhere; }
    table.data th:first-child, table.data td:first-child { text-align:left; } img { width:100%; max-width:1000px; display:block; margin:12px auto; }
    .note { color:#566773; font-size:12px; } code { background:#EDF1F4; padding:1px 4px; } .page { break-before:page; }
    ol { margin-top:4px; margin-bottom:4px; } li { margin-top:1px; margin-bottom:1px; }
    @media print { .screen-only { display:none; } }
    """
    document = f"""<!doctype html><html><head><meta charset='utf-8'><title>gdTAI V4.2 development gate</title><style>{css}</style></head><body>
    <h1>gdTAI V4.2: Development Gate Report</h1><p class='subtitle'>Expression-independent truth, multi-source NK controls, nested source transfer, and recovery diagnostics | 19 August 2026</p>
    <div class='decision'><b>Decision: do not promote V4.2.</b> The two-stage XGBoost cascade fails precommitted source-transfer requirements. No lockbox was scored, no release model was fitted, and no whole-atlas inference was run.</div>
    <div class='grid'><div class='metric'><b>{ground['n_manifest_rows']:,}</b>Labeled/audit cells</div><div class='metric'><b>{ground['n_nk_sources']}</b>Independent NK sources</div><div class='metric'><b>{nested['runtime_seconds']/60:.1f} min</b>Full A100 nested run</div><div class='metric'><b>0</b>Lockbox cells scored</div></div>

    <h2>1. What Was Frozen</h2><p>The analysis used the checksum-bound TCR-corrected 5.93-million-cell T/NK atlas (<code>{html.escape(ground['atlas_sha256'])}</code>) without modifying an H5AD. Ground truth was based on sorted gdT origin or productive paired receptors, not expression scores or model annotation.</p>
    <ul><li>When chain-level UMI information exists, both defining paired chains require UMI &ge;2; UMI=1 calls are excluded.</li><li>When a source truly lacks UMI information, published productive paired calls remain eligible.</li><li>Silver and dual-receptor cells have zero fitting and threshold-selection weight.</li><li>All eight extension cohorts, BALF_BLOOD_COPD, GDT_2020, and GDTlung remain locked.</li></ul>
    <h3>NK truth</h3>{table_html(nk[['source_gse_id','tier','n_author_nk','n_tcr_conflicts','n_eligible']])}
    <img src='exclusions_and_nk_truth.png'><p class='note'>The frozen CD4/Treg rules cost at most 1.47% recall in any gold-positive source and do not explain the Malte failure.</p>

    <h2 class='page'>2. Training And Evaluation</h2><p>Stage 1 separates transcriptomic T from curated NK. Stage 2 separates gdT from alpha-beta T using 197 individual genes, including TCR V/J/constant genes and a small T/NK context panel. Aggregate TRD/TRAB scores, source IDs, scVI coordinates, UMAP, and TCR metadata were not model features.</p>
    <p>Three source-held-out outer folds were used: HRA005041, GSE144469, and MalteGDT. Within each outer training set, group-preserving inner folds selected thresholds under frozen per-source recall and false-positive constraints. Twelve shallow GPU XGBoost candidates were evaluated.</p>
    <h3>Stage 1</h3>{table_html(stage1_summary, {'threshold_min':lambda x:f'{x:.4f}','threshold_max':lambda x:f'{x:.4f}'})}
    <h3>Stage 2 candidate eligibility</h3>{table_html(eligibility)}<img src='candidate_eligibility.png'>

    <h2 class='page'>3. Primary Results</h2>{table_html(primary_display, {'precision':pct,'recall':pct,'abt_fpr':pct,'f1':lambda x:'NA' if pd.isna(x) else f'{x:.3f}','pr_auc':lambda x:'NA' if pd.isna(x) else f'{x:.3f}'})}
    <img src='heldout_recall.png'><div class='ok'><b>HRA transfer is encouraging but insufficient.</b> Balanced recall is 95.19%; however, alpha-beta FPR is 0.359%, above the 0.2% inner target. GSE144 has no eligible inner profile. Malte recall collapses to 36.47%.</div>

    <h2>4. Failure Decomposition</h2><p>Malte Stage 1 passage is 99.79%, and only 1/7,800 Malte gold cells is removed by the fixed CD4/Treg rules. The loss therefore occurs in Stage 2. Feature-gain audits show dependence on TRDC/TRDV genes and, in some folds, the cross-fitted Stage 1 probability. This is consistent with nonlinear cohort-specific calibration rather than failure to recognize T lineage.</p>
    {table_html(gains.sort_values(['heldout_source','profile','gain'],ascending=[True,True,False]).groupby(['heldout_source','profile']).head(8), {'gain':lambda x:f'{x:,.1f}'})}

    <h2 class='page'>5. Recovery Attempts</h2><p>Two bounded development-only probes retained the same labels, groups, and gates while replacing Stage 2 with GPU L2 logistic regression and removing the stacked Stage 1 probability. Continuous receptor/control expression reached 46.72% held-out Malte recall. Detected/not-detected receptor features reached 41.18% balanced recall and 15.38% high-purity recall.</p>
    {table_html(recovery_display, {'precision':pct,'recall':pct,'f1':lambda x:'NA' if pd.isna(x) else f'{x:.3f}'})}
    <p>Neither recovery clears the frozen 70% per-source floor. Further Malte-directed hyperparameter search would convert the development cohort into a tuning target and increase overfitting risk, so it was stopped.</p>

    <h2>6. Relationship To V2 And V3</h2>{table_html(historical, {'precision':pct,'recall':pct,'f1':lambda x:f'{x:.3f}'})}<p class='note'><b>These historical numbers are descriptive only.</b> V2/V3 used different training sets, truth versions, validation cohorts, and in V2 high-purity mode an annotation-specific threshold. They cannot be interpreted as a fair head-to-head estimate against the new nested folds. The current result instead demonstrates why the older headline metrics do not rule out source overfitting.</p>

    <h2>7. Conclusion And Next Iteration</h2><ol><li>Do not promote, release-fit, or apply V4.2 to the whole atlas.</li><li>Preserve the untouched broad lockbox for a future precommitted iteration.</li><li>Do not continue adaptive searches on Malte. A new protocol should precommit source-invariant feature transformations or domain-generalization objectives before reusing these folds.</li><li>Retain V2/V3 as historical comparators, not as evidence that V4.2 should be relaxed.</li></ol>
    <p class='note screen-only'>Canonical tables: <code>Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_training/</code>. Recovery diagnostics: <code>Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_recovery/</code>. No H5AD was mutated.</p>
    </body></html>"""
    html_path = STATIC / "index.html"
    html_path.write_text(document)
    pdf_path = STATIC / "gdtai_v4_2_development_gate_report.pdf"
    subprocess.run([
        "google-chrome", "--headless", "--no-sandbox", "--disable-gpu",
        f"--print-to-pdf={pdf_path}", f"file://{html_path}",
    ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    summary = {
        "status": "FAIL_DEVELOPMENT_TRANSFER",
        "html": str(html_path.relative_to(ROOT)),
        "pdf": str(pdf_path.relative_to(ROOT)),
        "lockbox_scored": False,
        "model_promoted": False,
    }
    (LOG / "gdtai_v4_2_training/development_gate_report_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
