#!/usr/bin/env python3
"""Build a read-only audit of gdTAI training and evaluation provenance.

The workflow reads frozen model artifacts and existing evaluation tables. It does
not read or mutate the atlas H5AD. The generated report separates empirical
re-analysis from methodological findings and recommendations.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    f1_score,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)


ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_methodology_audit"
FIGURE_DIR = ROOT / "Integrated_dataset/figures/gdT_prediction/gdtai_methodology_audit"
LOG_DIR = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_methodology_audit"
STATIC_DIR = ROOT / "gdT_prediction/gdtai_methodology_audit"
STATIC_FIGURE_DIR = STATIC_DIR / "figures"

REGISTRY_PATH = ROOT / "configs/models/gdtai/model_registry.csv"
V2_MODEL_PATH = ROOT / "Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v2.0/gdTAI_v2_model.pkl"
V3_MODEL_PATH = ROOT / "Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/gdTAI_v3_model.pkl"
V3_MANIFEST_PATH = ROOT / "Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/model_manifest.json"
V3_RUN_SUMMARY_PATH = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v3_trdc_nk_guard/gdtai_v3_trdc_nk_guard_summary.json"
V3_TRAINING_COMPOSITION_PATH = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v3_trdc_nk_guard/training_sample_composition.csv"
V2_SPLIT_PATH = ROOT / "Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene/train_tune_split_by_source_label.csv"
V2_VALIDATION_PATH = ROOT / "Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene/validation_metrics.csv"
V3_GUARDRAIL_PATH = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v3_round12_vs_round14/promotion_guardrails.csv"
EXTERNAL_PREDICTIONS_PATH = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v3_trdc_nk_guard/external_predictions_wide.csv.gz"
EXTENSION_SCREEN_PATH = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_extension_evaluation/screen_cohort_summary.csv"

REPORT_HTML = STATIC_DIR / "index.html"
REPORT_PDF = STATIC_DIR / "gdtai_methodology_audit_report.pdf"
SUMMARY_MD = LOG_DIR / "gdtai_methodology_audit_summary.md"


@dataclass(frozen=True)
class ModelSpec:
    label: str
    score_column: str
    prediction_column: str
    operating_mode: str


MODEL_SPECS = (
    ModelSpec("V2 high-F1", "v2_high_f1_score", "v2_high_f1_pred", "fixed 0.9064"),
    ModelSpec("V2 high-purity", "v2_high_purity_score", "v2_high_purity_pred", "annotation-specific"),
    ModelSpec(
        "V3 Round 12",
        "v3_round12_hist_gradient_fixed_0p5_score",
        "v3_round12_hist_gradient_fixed_0p5_pred",
        "fixed 0.5",
    ),
    ModelSpec(
        "V3 Round 14",
        "v3_round20_v2_score_trdc_gate_fixed_0p936_score",
        "v3_round20_v2_score_trdc_gate_fixed_0p936_pred",
        "fixed 0.936",
    ),
)

EXTENSION_PROFILE_LABELS = {
    "v2_high_f1": "V2 high-F1",
    "v2_high_purity": "V2 high-purity",
    "v3_round12_high_purity": "V3 Round 12",
    "v3_round14_balanced": "V3 Round 14",
}

COLORS = {
    "V2 high-F1": "#7a6f63",
    "V2 high-purity": "#2f7d6d",
    "V3 Round 12": "#d58b28",
    "V3 Round 14": "#b4433c",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-pdf", action="store_true", help="Skip Chrome PDF export.")
    return parser.parse_args()


def ensure_directories() -> None:
    for path in (TABLE_DIR, FIGURE_DIR, LOG_DIR, STATIC_DIR, STATIC_FIGURE_DIR):
        path.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def as_bool(series: pd.Series) -> np.ndarray:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False).to_numpy(dtype=bool)
    return (
        series.astype("string")
        .str.strip()
        .str.lower()
        .isin({"true", "1", "yes", "t"})
        .to_numpy(dtype=bool)
    )


def metric_row(label: str, y: np.ndarray, score: np.ndarray, pred: np.ndarray, mode: str) -> dict[str, Any]:
    y = y.astype(np.int8, copy=False)
    pred = pred.astype(bool, copy=False)
    tp = int((pred & (y == 1)).sum())
    fp = int((pred & (y == 0)).sum())
    tn = int(((~pred) & (y == 0)).sum())
    fn = int(((~pred) & (y == 1)).sum())
    return {
        "model": label,
        "operating_mode": mode,
        "n_cells": int(y.size),
        "n_positive": int((y == 1).sum()),
        "n_negative": int((y == 0).sum()),
        "tp": tp,
        "fp": fp,
        "tn": tn,
        "fn": fn,
        "precision": float(precision_score(y, pred, zero_division=0)),
        "recall": float(recall_score(y, pred, zero_division=0)),
        "specificity": float(tn / (tn + fp)) if (tn + fp) else np.nan,
        "f1": float(f1_score(y, pred, zero_division=0)),
        "balanced_accuracy": float(balanced_accuracy_score(y, pred)),
        "mcc": float(matthews_corrcoef(y, pred)),
        "roc_auc": float(roc_auc_score(y, score)),
        "pr_auc": float(average_precision_score(y, score)),
    }


def audit_model_artifacts() -> pd.DataFrame:
    registry = pd.read_csv(REGISTRY_PATH)
    rows: list[dict[str, Any]] = []
    for _, row in registry.iterrows():
        artifact = ROOT / str(row["artifact_path"])
        actual = sha256(artifact) if artifact.exists() else "missing"
        rows.append(
            {
                "model_id": row["model_id"],
                "version": row["version"],
                "status": row["status"],
                "artifact_path": row["artifact_path"],
                "registry_sha256": row["sha256"],
                "actual_sha256": actual,
                "checksum_match": bool(actual == row["sha256"]),
                "size_bytes": int(artifact.stat().st_size) if artifact.exists() else np.nan,
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(TABLE_DIR / "artifact_integrity.csv", index=False)
    if not result["checksum_match"].all():
        raise RuntimeError("At least one registered model artifact failed checksum verification.")
    return result


def inspect_v2_coefficients() -> tuple[pd.DataFrame, pd.DataFrame]:
    payload = joblib.load(V2_MODEL_PATH)
    base = payload["base_model"]
    pipeline = base["model_object"]
    scaler = pipeline.named_steps["standardscaler"]
    logistic = pipeline.named_steps["logisticregression"]
    features = np.asarray(base["feature_names"], dtype=object)
    coefficients = logistic.coef_[0].astype(float)
    table = pd.DataFrame(
        {
            "feature": features,
            "standardized_coefficient": coefficients,
            "absolute_coefficient": np.abs(coefficients),
            "scaler_mean": scaler.mean_,
            "scaler_scale": scaler.scale_,
        }
    ).sort_values("absolute_coefficient", ascending=False)
    table.to_csv(TABLE_DIR / "v2_standardized_coefficients.csv", index=False)
    summary = pd.DataFrame(
        [
            {
                "n_features": int(features.size),
                "n_nonzero_coefficients": int((coefficients != 0).sum()),
                "minimum_coefficient": float(coefficients.min()),
                "maximum_coefficient": float(coefficients.max()),
                "maximum_absolute_coefficient": float(np.abs(coefficients).max()),
                "l1_norm": float(np.abs(coefficients).sum()),
                "uses_standard_scaler": True,
                "interpretation": "Magnitude is not intrinsically excessive after standardization; split and selection design dominate overfitting risk.",
            }
        ]
    )
    summary.to_csv(TABLE_DIR / "v2_coefficient_summary.csv", index=False)
    return table, summary


def audit_v3_semantic_provenance() -> pd.DataFrame:
    workflow_dir = str(ROOT / "workflows/gdtai")
    if workflow_dir not in sys.path:
        sys.path.insert(0, workflow_dir)
    payload = joblib.load(V3_MODEL_PATH)
    manifest = json.loads(V3_MANIFEST_PATH.read_text(encoding="utf-8"))
    run_summary = json.loads(V3_RUN_SUMMARY_PATH.read_text(encoding="utf-8"))
    composition = pd.read_csv(V3_TRAINING_COMPOSITION_PATH).set_index("sample_class")["n_cells"]
    row = {
        "manifest_model_name": manifest.get("model_name"),
        "manifest_status": manifest.get("status"),
        "manifest_accepted_for_promotion": manifest.get("accepted_for_promotion"),
        "payload_model": payload.get("model"),
        "payload_version": payload.get("version"),
        "payload_accepted_for_promotion": payload.get("accepted_for_promotion"),
        "completed_run_selected_candidate": run_summary.get("selected_candidate"),
        "completed_run_accepted_for_promotion": run_summary.get("accepted_for_promotion"),
        "recorded_training_negatives": int(composition.get("primary_or_hard_negative", 0)),
        "methodology_documented_training_negatives": 100_000,
        "semantic_consistency": False,
        "interpretation": "The checksum-pinned artifact is identifiable, but manifest, embedded payload, completed-run summary, and methodology do not describe one coherent build record.",
    }
    result = pd.DataFrame([row])
    result.to_csv(TABLE_DIR / "v3_semantic_provenance_audit.csv", index=False)
    return result


def corrected_external_benchmark() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    required = {
        "tcr_strict_cell_label",
        "real_abt_tcr_strict",
        "library_id",
        *(spec.score_column for spec in MODEL_SPECS),
        *(spec.prediction_column for spec in MODEL_SPECS),
    }
    data = pd.read_csv(EXTERNAL_PREDICTIONS_PATH, usecols=sorted(required), low_memory=False)
    strict_label = data["tcr_strict_cell_label"].astype("string").fillna("unknown")
    pure_positive = strict_label.isin(["real_gdT", "cd4_like_gdTCR_excluded"]).to_numpy(dtype=bool)
    pure_negative = as_bool(data["real_abt_tcr_strict"])
    primary = pure_positive | pure_negative
    y = pure_positive[primary].astype(np.int8)

    overall_rows: list[dict[str, Any]] = []
    library_rows: list[dict[str, Any]] = []
    subgroup_rows: list[dict[str, Any]] = []
    for spec in MODEL_SPECS:
        score = pd.to_numeric(data.loc[primary, spec.score_column], errors="coerce").to_numpy(dtype=float)
        pred_all = as_bool(data[spec.prediction_column])
        pred = pred_all[primary]
        overall_rows.append(metric_row(spec.label, y, score, pred, spec.operating_mode))
        for library, group_index in data.loc[primary].groupby("library_id", dropna=False).groups.items():
            positions = data.index.get_indexer(group_index)
            group_y = pure_positive[positions].astype(np.int8)
            group_score = pd.to_numeric(data.loc[group_index, spec.score_column], errors="coerce").to_numpy(dtype=float)
            group_pred = pred_all[positions]
            row = metric_row(spec.label, group_y, group_score, group_pred, spec.operating_mode)
            row["library_id"] = str(library)
            library_rows.append(row)
        for subgroup, subgroup_mask in {
            "TCR-only gdT retained by original label": strict_label.eq("real_gdT").to_numpy(dtype=bool),
            "TCR-only gdT excluded by CD4-like expression rule": strict_label.eq("cd4_like_gdTCR_excluded").to_numpy(dtype=bool),
            "All paired productive gdT/no productive abT": pure_positive,
        }.items():
            n = int(subgroup_mask.sum())
            subgroup_rows.append(
                {
                    "model": spec.label,
                    "subgroup": subgroup,
                    "n_positive": n,
                    "true_positive": int((pred_all & subgroup_mask).sum()),
                    "false_negative": int(((~pred_all) & subgroup_mask).sum()),
                    "recall": float(pred_all[subgroup_mask].mean()) if n else np.nan,
                }
            )

    overall = pd.DataFrame(overall_rows)
    per_library = pd.DataFrame(library_rows)
    subgroup = pd.DataFrame(subgroup_rows)
    overall.to_csv(TABLE_DIR / "corrected_external_pure_tcr_metrics.csv", index=False)
    per_library.to_csv(TABLE_DIR / "corrected_external_pure_tcr_metrics_by_library.csv", index=False)
    subgroup.to_csv(TABLE_DIR / "corrected_external_positive_subgroup_recall.csv", index=False)

    if int(pure_positive.sum()) != 1046 or int(pure_negative.sum()) != 33117:
        raise RuntimeError("Corrected TCR-only external cohort counts changed; audit assumptions must be reviewed.")
    return overall, per_library, subgroup


def prevalence_scenarios(metrics: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for _, metric in metrics.iterrows():
        sensitivity = float(metric["recall"])
        specificity = float(metric["specificity"])
        for prevalence in (0.001, 0.0025, 0.005, 0.01, 0.02, 0.05):
            tp_rate = sensitivity * prevalence
            fp_rate = (1.0 - specificity) * (1.0 - prevalence)
            ppv = tp_rate / (tp_rate + fp_rate)
            rows.append(
                {
                    "model": metric["model"],
                    "prevalence": prevalence,
                    "sensitivity_basis": sensitivity,
                    "specificity_basis": specificity,
                    "expected_true_positive_per_million": int(round(tp_rate * 1_000_000)),
                    "expected_false_positive_per_million": int(round(fp_rate * 1_000_000)),
                    "expected_predicted_positive_per_million": int(round((tp_rate + fp_rate) * 1_000_000)),
                    "expected_ppv": ppv,
                    "expected_fdr": 1.0 - ppv,
                    "basis": "Corrected BALF_BLOOD_COPD TCR-only benchmark; reused evaluation cohort, not an untouched external test.",
                }
            )
    result = pd.DataFrame(rows)
    result.to_csv(TABLE_DIR / "prevalence_ppv_scenarios.csv", index=False)
    return result


def extension_negative_summary() -> tuple[pd.DataFrame, pd.DataFrame]:
    data = pd.read_csv(EXTENSION_SCREEN_PATH)
    data = data[data["profile_id"].isin(EXTENSION_PROFILE_LABELS)].copy()
    data["model"] = data["profile_id"].map(EXTENSION_PROFILE_LABELS)
    data["strict_NK_n"] = pd.to_numeric(data["strict_NK_n"], errors="coerce").fillna(0).astype(int)
    data["strict_NK_fp"] = pd.to_numeric(data["strict_NK_fp"], errors="coerce").fillna(0).astype(int)
    data["abT_n"] = pd.to_numeric(data["abT_n"], errors="coerce").fillna(0).astype(int)
    data["abT_fp"] = pd.to_numeric(data["abT_fp"], errors="coerce").fillna(0).astype(int)
    data["paired_abT_n"] = pd.to_numeric(data["paired_abT_n"], errors="coerce").fillna(0).astype(int)
    data["paired_abT_fp"] = pd.to_numeric(data["paired_abT_fp"], errors="coerce").fillna(0).astype(int)
    pooled = (
        data.groupby("model", sort=False)
        .agg(
            n_datasets=("dataset_id", "nunique"),
            total_cells=("total_cells", "sum"),
            predicted_gdT=("predicted_gdT", "sum"),
            abT_n=("abT_n", "sum"),
            abT_fp=("abT_fp", "sum"),
            paired_abT_n=("paired_abT_n", "sum"),
            paired_abT_fp=("paired_abT_fp", "sum"),
            strict_NK_n=("strict_NK_n", "sum"),
            strict_NK_fp=("strict_NK_fp", "sum"),
        )
        .reset_index()
    )
    pooled["predicted_rate"] = pooled["predicted_gdT"] / pooled["total_cells"]
    pooled["abT_fpr"] = pooled["abT_fp"] / pooled["abT_n"]
    pooled["paired_abT_fpr"] = pooled["paired_abT_fp"] / pooled["paired_abT_n"]
    pooled["strict_NK_fpr"] = pooled["strict_NK_fp"] / pooled["strict_NK_n"]
    data.to_csv(TABLE_DIR / "extension_negative_challenge_by_dataset.csv", index=False)
    pooled.to_csv(TABLE_DIR / "extension_negative_challenge_pooled.csv", index=False)
    return data, pooled


def build_issue_register() -> pd.DataFrame:
    rows = [
        {
            "issue_id": "AUD-01",
            "severity": "Critical",
            "stage": "V2 model selection",
            "finding": "The combined cohort named validation was used to accept and select the V2 algorithm.",
            "evidence": "workflows/gdtai/run_gdt_gse144469_holdout_tcrgene_classifier.py:784-831; gdTAI_v2.0/METHODOLOGY.md:304-320",
            "impact": "Reported V2 validation F1 is a selection estimate, not an untouched final-test estimate.",
            "certainty": "Confirmed from code and output tables",
            "recommended_action": "Re-estimate with nested grouped outer folds; reserve model choice and thresholding for inner folds.",
        },
        {
            "issue_id": "AUD-02",
            "severity": "Critical",
            "stage": "V3 promotion",
            "finding": "The cohort called independent external was included in Round 12 versus Round 14 guardrails and mean-F1 promotion.",
            "evidence": "workflows/gdtai/compare_gdtai_v3_round12_vs_round14.py:850-910; gdTAI_v3.0/METHODOLOGY.md:13",
            "impact": "The cohort is no longer independent of model selection, so its metrics are optimistic for generalization claims.",
            "certainty": "Confirmed from code, manifest, and promotion JSON",
            "recommended_action": "Rename it a reused cross-study benchmark and stop using it for tuning or promotion.",
        },
        {
            "issue_id": "AUD-03",
            "severity": "High",
            "stage": "Train/tune splitting",
            "finding": "Within each source-label stratum, individual cells were shuffled into train and tune splits.",
            "evidence": "run_gdt_gse144469_holdout_tcrgene_classifier.py:486-505; run_gdtai_v3_trdc_nk_guard_classifier.py:635-653",
            "impact": "Donor, library, and clonotype-correlated expression can appear in both splits and inflate tuning performance.",
            "certainty": "Confirmed from code",
            "recommended_action": "Group by dataset then donor/library, and by clonotype where TCR clonotypes are available.",
        },
        {
            "issue_id": "AUD-04",
            "severity": "High",
            "stage": "External ground truth",
            "finding": "The external gdT label excludes paired productive gdT/no-abT cells using CD4 annotation or CD4 expression, while CD4 is a model feature; external negatives also include expression annotations.",
            "evidence": "/home/tanlikai/databank/owndata/singlecell/data/phase4_tcell_focused_analysis.py:292-311; run_gdtai_v3_trdc_nk_guard_classifier.py:1676-1705",
            "impact": "Feature-label circularity can make the benchmark easier and hides performance on CD4-like TCR-confirmed gdT cells.",
            "certainty": "Confirmed from label-building code",
            "recommended_action": "Use TCR-only labels for primary metrics and report CD4-like, NK-like, and low-CD3 cells as stress strata.",
        },
        {
            "issue_id": "AUD-05",
            "severity": "High",
            "stage": "Validation composition",
            "finding": "Large validation components are positive-only cord blood and negative-only GSE254249, linking class to cohort and biology.",
            "evidence": "gdTAI_v2.0/METHODOLOGY.md:155-182; split_overall.csv",
            "impact": "A pooled cell-level metric partly measures source/platform separation and is dominated by artificial class composition.",
            "certainty": "Confirmed from split counts",
            "recommended_action": "Report within-dataset contrasts, dataset-macro metrics, and leave-one-dataset-out performance.",
        },
        {
            "issue_id": "AUD-06",
            "severity": "High",
            "stage": "Adaptive iteration",
            "finding": "Many candidate rounds, thresholds, gates, and post-hoc comparisons repeatedly inspected the same atlas, held-out, and external cohorts.",
            "evidence": "external_predictions_wide.csv.gz contains candidate rounds 02-38; Round 12/14 promotion reuses all three evaluation frames",
            "impact": "Repeated adaptive inspection creates selection overfitting even when model fitting never reads the external labels directly.",
            "certainty": "Confirmed from retained prediction schema and reports",
            "recommended_action": "Pre-register a compact candidate grid and evaluation rule before the next run; keep all decisions inside nested resampling.",
        },
        {
            "issue_id": "AUD-07",
            "severity": "High",
            "stage": "Deployment prevalence",
            "finding": "Training/evaluation prevalence is higher and more controlled than deployment prevalence in many T/NK cohorts.",
            "evidence": "V2 validation prevalence 6.58%; corrected external prevalence 3.06%; extension cohorts contain no positive truth",
            "impact": "Precision can fall sharply at 0.1-1% prevalence despite specificity near 99.9%.",
            "certainty": "Confirmed; deployment prevalence varies",
            "recommended_action": "Publish PPV/FDR curves at realistic prevalence and select operating modes by explicit error costs.",
        },
        {
            "issue_id": "AUD-08",
            "severity": "High",
            "stage": "Positive labels",
            "finding": "All cells in sorted gdT datasets are treated as gold positives, including possible contaminants, doublets, and low-quality cells.",
            "evidence": "gdTAI_v2.0/METHODOLOGY.md ground-truth definition and sorted-source counts",
            "impact": "Label noise can reward source signatures and penalize biologically valid but weak-expression gdT cells.",
            "certainty": "Design risk; contamination rate is not directly measured",
            "recommended_action": "Run source-exclusion sensitivity analyses and robust/noisy-label training; do not silently relabel by model features.",
        },
        {
            "issue_id": "AUD-09",
            "severity": "High",
            "stage": "Feature design",
            "finding": "Individual TCR V/J genes can encode donor repertoire, clonotype, tissue, and chemistry; stability was not established for the released V2/Round14 coefficients.",
            "evidence": "gdTAI_v2.0/METHODOLOGY.md:202-220; v2_standardized_coefficients.csv; bootstrap table supports a different elastic-net candidate",
            "impact": "Good within-source discrimination may not transfer when repertoire or gene detection changes.",
            "certainty": "Biologically plausible and incompletely tested",
            "recommended_action": "Compare family-level sums, detection indicators, and individual V/J genes using outer-fold stability selection and missing-panel stress tests.",
        },
        {
            "issue_id": "AUD-10",
            "severity": "High",
            "stage": "Full-atlas FP estimation",
            "finding": "Hidden abT false positives in non-TCR-sequenced sources are extrapolated from observed paired-abT errors in TCR-sequenced sources.",
            "evidence": "workflows/gdtai/run_gdtai_v3_trdc_nk_guard_classifier.py:2091-2099",
            "impact": "The estimate assumes exchangeability across studies and does not capture NK or other non-abT errors under domain shift.",
            "certainty": "Confirmed assumption in code",
            "recommended_action": "Present this as a sensitivity range stratified by dataset/technology, not a single estimated FP fraction.",
        },
        {
            "issue_id": "AUD-11",
            "severity": "High",
            "stage": "Domain-shift challenge",
            "finding": "Frozen-model extension cohorts show materially higher and dataset-variable NK/abT false-positive rates than the reused external benchmark.",
            "evidence": "gdtai_extension_evaluation/screen_cohort_summary.csv and extension_negative_challenge_pooled.csv",
            "impact": "One pooled external specificity does not describe cross-study deployment risk.",
            "certainty": "Confirmed on negative-control strata; these cohorts lack positive truth",
            "recommended_action": "Use dataset-level stress panels and macro/maximum FPR guardrails, while keeping sensitivity evaluation separate.",
        },
        {
            "issue_id": "AUD-12",
            "severity": "Moderate",
            "stage": "Reproducibility governance",
            "finding": "Model checksums are registered, but frozen manifests do not fully pin input H5AD hashes, split cell IDs, Git commit, package lock, and fold-specific feature selection.",
            "evidence": "configs/models/gdtai/model_registry.csv and gdTAI_v2.0/gdTAI_v3.0 model_manifest.json",
            "impact": "The artifact is verifiable, but exact training reconstruction and leakage auditing are incomplete.",
            "certainty": "Confirmed from manifests",
            "recommended_action": "Add immutable run manifests with data hashes, environment lock, split IDs, fold definitions, seed, commit, and output checksums.",
        },
        {
            "issue_id": "AUD-13",
            "severity": "High",
            "stage": "Full-atlas candidate selection",
            "finding": "V3 candidate selection used full-atlas recall, F1, TRDC/TRDV composition, and estimated FP, even though the full atlas includes training cells.",
            "evidence": "workflows/gdtai/run_gdtai_v3_trdc_nk_guard_classifier.py:2188-2255",
            "impact": "The atlas-wide performance and composition targets are development criteria, not independent evidence of generalization.",
            "certainty": "Confirmed independently from code",
            "recommended_action": "Move every selection criterion into nested outer-fold out-of-fold predictions and treat full-atlas application as post-selection inference only.",
        },
        {
            "issue_id": "AUD-14",
            "severity": "High",
            "stage": "Release provenance",
            "finding": "The Round 14 manifest says promoted, its pickle embeds accepted_for_promotion=False and candidate versioning, the completed-run summary selects Round 12, and training-negative counts disagree with methodology.",
            "evidence": "gdTAI_v3.0/model_manifest.json; gdTAI_v3_model.pkl; gdtai_v3_trdc_nk_guard_summary.json; training_sample_composition.csv",
            "impact": "The artifact checksum is sound, but the promoted build cannot be reconstructed from one internally coherent run record.",
            "certainty": "Confirmed independently and rechecked by main agent",
            "recommended_action": "Create one immutable release bundle and semantic-consistency test spanning pickle fields, manifest, registry, run summary, command, data hashes, and composition.",
        },
        {
            "issue_id": "AUD-15",
            "severity": "High",
            "stage": "Extension holdout governance",
            "finding": "The extension config names sealed holdouts, but screen-manifest loading forces sealed_holdout=False; minimum model-gene completeness is also configured as zero.",
            "evidence": "configs/models/gdtai/extension_evaluation.json; workflows/gdtai/compare_frozen_gdtai_profiles.py:190-203",
            "impact": "Screened cohorts are no longer prospective holdouts, and incomplete gene panels can produce misleadingly confident calls.",
            "certainty": "Confirmed independently from config and code",
            "recommended_action": "Relabel screened cohorts as development stress sets; enforce a gene-completeness floor and abstain when critical genes or sufficient model coverage are absent.",
        },
        {
            "issue_id": "AUD-16",
            "severity": "Moderate",
            "stage": "Statistical uncertainty",
            "finding": "Canonical comparisons emphasize pooled cell metrics and cell-level paired tests without donor/dataset-cluster confidence intervals.",
            "evidence": "gdtai_v3_round12_vs_round14/paired_prediction_tests.csv and per-dataset comparison tables",
            "impact": "Large cell counts overstate effective sample size; model ranking may change across donors or datasets.",
            "certainty": "Confirmed reporting gap",
            "recommended_action": "Report donor- and dataset-macro metrics with clustered bootstrap confidence intervals and worst-dataset guardrails.",
        },
        {
            "issue_id": "AUD-17",
            "severity": "Moderate",
            "stage": "Calibration",
            "finding": "Round 14 is a gated probability-like score rather than a calibrated probability, while its threshold was chosen under enriched development prevalence.",
            "evidence": "gdTAI_v3.0/USAGE_PROTOCOL.md and promotion workflow",
            "impact": "A score of 0.936 cannot be interpreted as a 93.6% probability, and fixed thresholds may not transfer across prevalence or chemistry.",
            "certainty": "Confirmed documentation and design limitation",
            "recommended_action": "Cross-fit calibration inside grouped folds and publish prevalence-specific PPV/FDR plus an abstention state.",
        },
    ]
    table = pd.DataFrame(rows)
    order = pd.Categorical(table["severity"], categories=["Critical", "High", "Moderate", "Low"], ordered=True)
    table = table.assign(_order=order).sort_values(["_order", "issue_id"]).drop(columns="_order")
    table.to_csv(TABLE_DIR / "issue_register.csv", index=False)
    return table


def build_recommendation_register() -> pd.DataFrame:
    rows = [
        ("P0", "Correct claims and freeze comparators", "Treat V2, Round 12, and Round 14 as legacy comparators. Rename BALF_BLOOD_COPD a reused cross-study benchmark.", "Documentation-only; prevents unsupported independent-validation claims."),
        ("P0", "Lock the analysis specification", "Before rerunning, freeze labels, outer/inner groups, candidate families, feature sets, metrics, guardrails, and promotion rule in a committed config.", "Stops further adaptive evaluation leakage."),
        ("P0", "Rebuild expression-independent labels", "Use paired productive gd/no productive ab as positives and paired/single productive ab/no gd as separately reported negatives. Keep transcriptomic flags only as stress strata.", "Removes CD4/NK feature-label circularity."),
        ("P1", "Nested leave-one-dataset-out evaluation", "Outer folds hold out complete datasets. Inner folds group donor/library and clonotype when available; preprocessing, feature selection, calibration, and thresholds occur inside inner training only.", "Best available honest generalization estimate without new data."),
        ("P1", "Dataset-balanced learning", "Cap or weight each dataset and class so HRA005041 does not dominate. Compare matched within-dataset positives and negatives wherever possible.", "Reduces source and prevalence confounding."),
        ("P1", "Feature-family ablation and stability", "Compare regularized logistic models using family sums/detection, individual V/J genes, transcriptome lineage features, and combinations. Retain features stable across outer folds.", "Tests whether repertoire-specific V/J signals are portable."),
        ("P1", "Soft two-stage T-lineage model", "Generate an out-of-fold high-sensitivity T-lineage probability, then combine it with the gd-versus-ab score. Do not hard-exclude NK-marker-high or TRDV-dropout cells.", "Addresses NK FPs while preserving cytotoxic and dropout-prone gdT cells."),
        ("P1", "Prevalence-aware calibration and abstention", "Calibrate only within inner folds; publish PPV/FDR at 0.1%, 0.25%, 0.5%, 1%, 2%, and 5% prevalence; abstain when critical genes or sufficient panel coverage are missing.", "Makes operating modes interpretable in rare-cell deployment and prevents confident calls from incomplete panels."),
        ("P1", "Failure and robustness panel", "Report strict NK, NKT/MAIT, CD8, CD4/Treg, low-CD3, TRDC+TRDV-, doublet/QC proxies, missing genes, 3-prime versus 5-prime, tissue, disease, and source strata.", "Surfaces clinically relevant failure modes hidden by pooled metrics."),
        ("P2", "Uncertainty and model simplification", "Use dataset bootstrap confidence intervals and compare the simplest elastic-net/logistic model against gates and trees. Promote complexity only for consistent macro-fold gains.", "Controls variance and supports portable deployment."),
        ("P2", "Leakage diagnostics", "Run label permutation within dataset/donor and predict source from candidate features. Large residual performance or source predictability indicates confounding.", "Directly tests source leakage using current data."),
        ("P2", "Immutable semantic release check", "Require pickle fields, registry, manifest, run summary, training composition, command, Git commit, environment, and input hashes to agree before promotion.", "Makes the selected build exactly reconstructable and catches stale mixed-run metadata."),
    ]
    table = pd.DataFrame(rows, columns=["priority", "action", "implementation", "purpose"])
    table.to_csv(TABLE_DIR / "recommended_actions.csv", index=False)
    return table


def plot_corrected_performance(metrics: pd.DataFrame) -> Path:
    measures = ["precision", "recall", "f1", "specificity"]
    labels = ["Precision", "Recall", "F1", "Specificity"]
    x = np.arange(len(measures))
    width = 0.19
    fig, ax = plt.subplots(figsize=(10.8, 5.4), constrained_layout=True)
    for i, (_, row) in enumerate(metrics.iterrows()):
        values = [float(row[m]) for m in measures]
        ax.bar(x + (i - 1.5) * width, values, width=width, label=row["model"], color=COLORS[row["model"]])
    ax.set_xticks(x, labels)
    ax.set_ylim(0.78, 1.005)
    ax.set_ylabel("Metric")
    ax.set_title("Corrected TCR-only benchmark: operating-point performance")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, ncol=2, loc="lower right")
    ax.text(
        0.01,
        0.01,
        "1,046 paired productive gdT/no productive abT positives; 33,117 strict paired abT/no gdT negatives.\nThis is a reused benchmark, not an untouched external test.",
        transform=ax.transAxes,
        fontsize=8.5,
        color="#4c4a46",
    )
    path = FIGURE_DIR / "corrected_tcr_only_performance.png"
    fig.savefig(path, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def plot_prevalence_ppv(scenarios: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(9.4, 5.4), constrained_layout=True)
    for model, group in scenarios.groupby("model", sort=False):
        ax.plot(
            group["prevalence"] * 100,
            group["expected_ppv"] * 100,
            marker="o",
            linewidth=2.2,
            label=model,
            color=COLORS[model],
        )
    ax.set_xscale("log")
    ax.set_xticks([0.1, 0.25, 0.5, 1, 2, 5], ["0.1", "0.25", "0.5", "1", "2", "5"])
    ax.set_xlabel("True gdT prevalence (%)")
    ax.set_ylabel("Expected positive predictive value (%)")
    ax.set_ylim(0, 100)
    ax.set_title("Purity depends strongly on deployment prevalence")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, ncol=2)
    path = FIGURE_DIR / "prevalence_ppv_scenarios.png"
    fig.savefig(path, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def plot_extension_challenge(data: pd.DataFrame, pooled: pd.DataFrame) -> Path:
    filtered = data[data["strict_NK_n"] >= 100].copy()
    datasets = sorted(filtered["dataset_id"].unique())
    fig, ax = plt.subplots(figsize=(11.6, 5.8), constrained_layout=True)
    x = np.arange(len(datasets))
    width = 0.19
    for i, model in enumerate(COLORS):
        group = filtered[filtered["model"] == model].set_index("dataset_id")
        values = [100 * float(group.loc[d, "strict_NK_fpr"]) if d in group.index else np.nan for d in datasets]
        ax.bar(x + (i - 1.5) * width, values, width=width, color=COLORS[model], label=model)
    ax.set_xticks(x, datasets, rotation=28, ha="right")
    ax.set_ylabel("Strict NK false-positive rate (%)")
    ax.set_title("Frozen-model NK challenge varies substantially by dataset")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, ncol=2)
    ax.text(
        0.01,
        0.98,
        "Datasets with fewer than 100 strict-NK controls are omitted from bars. These cohorts have no positive truth.",
        transform=ax.transAxes,
        va="top",
        fontsize=8.5,
        color="#4c4a46",
    )
    path = FIGURE_DIR / "extension_strict_nk_challenge.png"
    fig.savefig(path, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def plot_evaluation_reuse_map() -> Path:
    fig, ax = plt.subplots(figsize=(11.3, 4.8), constrained_layout=True)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4)
    ax.axis("off")
    boxes = [
        (0.35, 2.45, 2.15, 0.85, "Train/tune cells", "Random within source-label\nV2 and V3"),
        (3.05, 2.45, 2.15, 0.85, "Named validation", "Compared candidate algorithms\nand selected V2"),
        (5.75, 2.45, 2.15, 0.85, "BALF/BLOOD/COPD", "Repeated candidate evaluation\nplus Round 12/14 promotion"),
        (8.45, 2.45, 1.25, 0.85, "Release", "Round 14"),
    ]
    for i, (x, y, w, h, title, subtitle) in enumerate(boxes):
        color = "#f4ece4" if i == 0 else "#f8e2df"
        edge = "#8a7562" if i == 0 else "#b4433c"
        ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=color, edgecolor=edge, linewidth=1.6))
        ax.text(x + w / 2, y + 0.57, title, ha="center", va="center", fontsize=10.5, weight="bold", color="#272522")
        ax.text(x + w / 2, y + 0.25, subtitle, ha="center", va="center", fontsize=8.2, color="#4c4a46")
        if i < len(boxes) - 1:
            next_x = boxes[i + 1][0]
            ax.annotate("", xy=(next_x - 0.08, y + h / 2), xytext=(x + w + 0.08, y + h / 2), arrowprops={"arrowstyle": "->", "color": "#6f6a64", "lw": 1.4})
    ax.text(5.0, 1.45, "Selection estimates", ha="center", fontsize=10.5, weight="bold", color="#b4433c")
    ax.text(5.0, 1.08, "The same evaluation evidence influenced algorithm or release choice.", ha="center", fontsize=9.3, color="#4c4a46")
    ax.text(5.0, 0.55, "Remedy: nested dataset-level outer folds; donor/library/clonotype-grouped inner folds; precommitted promotion rule.", ha="center", fontsize=9.5, weight="bold", color="#2f7d6d")
    path = FIGURE_DIR / "evaluation_reuse_map.png"
    fig.savefig(path, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def fmt_int(value: Any) -> str:
    return f"{int(value):,}"


def fmt_pct(value: Any, digits: int = 2) -> str:
    if pd.isna(value):
        return "NA"
    return f"{100 * float(value):.{digits}f}%"


def table_html(frame: pd.DataFrame, columns: list[tuple[str, str, str]], classes: str = "") -> str:
    head = "".join(f"<th>{html.escape(label)}</th>" for _, label, _ in columns)
    body_rows = []
    for _, row in frame.iterrows():
        cells = []
        for key, _, kind in columns:
            value = row[key]
            if kind == "int":
                rendered = fmt_int(value)
            elif kind == "pct":
                rendered = fmt_pct(value)
            elif kind == "float3":
                rendered = f"{float(value):.3f}"
            else:
                rendered = html.escape(str(value))
            cells.append(f"<td>{rendered}</td>")
        body_rows.append("<tr>" + "".join(cells) + "</tr>")
    return f'<div class="table-wrap"><table class="{classes}"><thead><tr>{head}</tr></thead><tbody>{"".join(body_rows)}</tbody></table></div>'


def render_html(
    issues: pd.DataFrame,
    recommendations: pd.DataFrame,
    corrected: pd.DataFrame,
    per_library: pd.DataFrame,
    subgroup: pd.DataFrame,
    scenarios: pd.DataFrame,
    extension_pooled: pd.DataFrame,
    coefficient_table: pd.DataFrame,
    coefficient_summary: pd.DataFrame,
    artifact_integrity: pd.DataFrame,
    v3_provenance: pd.DataFrame,
) -> None:
    round12 = corrected.set_index("model").loc["V3 Round 12"]
    round14 = corrected.set_index("model").loc["V3 Round 14"]
    scenario_one_pct = scenarios[np.isclose(scenarios["prevalence"], 0.01)].copy()
    top_coefficients = pd.concat(
        [
            coefficient_table.nsmallest(8, "standardized_coefficient"),
            coefficient_table.nlargest(8, "standardized_coefficient"),
        ]
    ).drop_duplicates("feature")
    top_coefficients = top_coefficients.sort_values("standardized_coefficient")
    issue_rows = []
    for _, row in issues.iterrows():
        issue_rows.append(
            f"<tr><td><span class='severity {row['severity'].lower()}'>{html.escape(row['severity'])}</span></td>"
            f"<td><strong>{html.escape(row['issue_id'])}</strong><br>{html.escape(row['stage'])}</td>"
            f"<td>{html.escape(row['finding'])}</td><td>{html.escape(row['impact'])}</td>"
            f"<td>{html.escape(row['recommended_action'])}</td></tr>"
        )
    recommendation_rows = "".join(
        f"<tr><td><span class='priority'>{html.escape(row['priority'])}</span></td><td><strong>{html.escape(row['action'])}</strong></td><td>{html.escape(row['implementation'])}</td><td>{html.escape(row['purpose'])}</td></tr>"
        for _, row in recommendations.iterrows()
    )
    artifact_rows = "".join(
        f"<tr><td>{html.escape(str(row['model_id']))}</td><td>{html.escape(str(row['status']))}</td><td>{fmt_int(row['size_bytes'])}</td><td><code>{html.escape(str(row['actual_sha256'])[:16])}...</code></td><td>{'PASS' if row['checksum_match'] else 'FAIL'}</td></tr>"
        for _, row in artifact_integrity.iterrows()
    )
    html_text = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>gdTAI training and evaluation audit</title>
<style>
:root {{ --ink:#272522; --muted:#69645e; --paper:#ffffff; --line:#d9d4cd; --soft:#f5f2ed; --red:#b4433c; --amber:#d58b28; --green:#2f7d6d; --taupe:#7a6f63; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:#ece9e4; color:var(--ink); font-family:Arial,Helvetica,sans-serif; line-height:1.42; letter-spacing:0; }}
main {{ max-width:1180px; margin:0 auto; background:var(--paper); min-height:100vh; }}
header {{ padding:42px 54px 34px; border-top:8px solid var(--red); border-bottom:1px solid var(--line); }}
h1 {{ margin:0 0 8px; font-size:34px; letter-spacing:0; }}
.subtitle {{ margin:0; color:var(--muted); font-size:16px; max-width:900px; }}
.meta {{ margin-top:18px; color:var(--muted); font-size:12px; }}
section {{ padding:32px 54px; border-bottom:1px solid var(--line); break-inside:auto; }}
h2 {{ margin:0 0 18px; font-size:23px; letter-spacing:0; }}
h3 {{ margin:24px 0 10px; font-size:16px; letter-spacing:0; }}
p {{ margin:8px 0 14px; }}
.lead {{ font-size:17px; max-width:960px; }}
.facts {{ display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:12px; margin:20px 0; }}
.fact {{ border:1px solid var(--line); border-top:4px solid var(--taupe); padding:15px; min-height:112px; }}
.fact strong {{ display:block; font-size:25px; margin-bottom:6px; }}
.fact span {{ color:var(--muted); font-size:12px; }}
.callout {{ border-left:5px solid var(--red); background:#fbf3f1; padding:16px 18px; margin:18px 0; }}
.callout.green {{ border-color:var(--green); background:#eef7f4; }}
.figure {{ margin:18px 0 8px; text-align:center; break-inside:avoid; }}
.figure img {{ width:100%; max-width:1030px; height:auto; }}
.caption {{ color:var(--muted); font-size:11px; text-align:left; margin-top:4px; }}
.table-wrap {{ overflow-x:auto; margin:12px 0 18px; }}
table {{ border-collapse:collapse; width:100%; font-size:11.5px; }}
th {{ background:#eeeae4; text-align:left; padding:8px 7px; border:1px solid var(--line); vertical-align:bottom; }}
td {{ padding:7px; border:1px solid var(--line); vertical-align:top; }}
tbody tr:nth-child(even) {{ background:#faf9f7; }}
.issue-table {{ font-size:10.5px; table-layout:fixed; }}
.issue-table th:nth-child(1) {{ width:8%; }} .issue-table th:nth-child(2) {{ width:12%; }} .issue-table th:nth-child(3) {{ width:25%; }} .issue-table th:nth-child(4) {{ width:24%; }} .issue-table th:nth-child(5) {{ width:31%; }}
.severity,.priority {{ display:inline-block; padding:2px 6px; border-radius:3px; font-size:9px; font-weight:bold; text-transform:uppercase; }}
.critical {{ background:#f3d2cf; color:#7b211c; }} .high {{ background:#f7e3c1; color:#71420d; }} .moderate {{ background:#e3e0dc; color:#47433e; }}
.priority {{ background:#dcece7; color:#205d50; }}
code {{ font-family:Consolas,monospace; font-size:10px; }}
ul {{ margin:8px 0 14px 20px; padding:0; }} li {{ margin:6px 0; }}
.small {{ color:var(--muted); font-size:11px; }}
footer {{ padding:24px 54px; color:var(--muted); font-size:11px; }}
@media (max-width:760px) {{ header,section,footer {{ padding-left:22px; padding-right:22px; }} .facts {{ grid-template-columns:1fr 1fr; }} h1 {{ font-size:28px; }} }}
@media print {{
  @page {{ size:A4 landscape; margin:10mm; }}
  body {{ background:white; font-size:9.5pt; }} main {{ max-width:none; }}
  header {{ padding:14mm 10mm 8mm; }} section {{ padding:8mm 10mm; }} footer {{ padding:6mm 10mm; }}
  h1 {{ font-size:24pt; }} h2 {{ font-size:17pt; }} h3 {{ font-size:12pt; }}
  .facts {{ gap:5mm; }} .fact {{ min-height:24mm; padding:4mm; }}
  .figure {{ margin:3mm 0 1.5mm; }} .figure img {{ width:auto; max-width:100%; max-height:75mm; object-fit:contain; }}
  table {{ font-size:7.2pt; }} .issue-table {{ font-size:6.2pt; }} th,td {{ padding:1.3mm; }}
  .print-omit {{ display:none; }}
  .page-break {{ break-before:page; }} .table-wrap {{ overflow:visible; }}
}}
</style>
</head>
<body><main>
<header>
  <h1>gdTAI training and evaluation audit</h1>
  <p class="subtitle">Independent methodology review, reconciled against frozen artifacts and executable repository code. Scope: data preparation, labels, feature construction, algorithms, splitting, tuning, evaluation, deployment claims, and reproducibility.</p>
  <p class="meta">Audit date: {datetime.now().astimezone().strftime('%Y-%m-%d %Z')} | Read-only analysis: no H5AD was opened or modified | Atlas reference: 5.13 million cells</p>
</header>

<section>
  <h2>Executive conclusion</h2>
  <p class="lead"><strong>gdTAI has real discriminatory signal, but the published internal numbers do not yet rule out selection overfitting.</strong> The largest problems are evaluation reuse, cell-level train/tune splits, source-class confounding, and a partly transcriptome-informed external label. These are more consequential than the size of the fitted coefficients.</p>
  <div class="facts">
    <div class="fact"><strong>{int((issues['severity'] == 'Critical').sum())}</strong><span>critical design findings requiring claim correction</span></div>
    <div class="fact"><strong>{fmt_pct(round12['f1'])}</strong><span>Round 12 F1 on corrected TCR-only benchmark</span></div>
    <div class="fact"><strong>{fmt_pct(round14['f1'])}</strong><span>Round 14 F1 on corrected TCR-only benchmark</span></div>
    <div class="fact"><strong>{float(coefficient_summary.iloc[0]['maximum_absolute_coefficient']):.2f}</strong><span>maximum absolute V2 coefficient after standardization</span></div>
  </div>
  <div class="callout"><strong>Claim boundary.</strong> BALF_BLOOD_COPD was used in Round 12 versus Round 14 promotion. It should be described as a reused cross-study benchmark, not an untouched independent external test. Current releases remain suitable as research comparators while evaluation is rebuilt.</div>
  <div class="callout green"><strong>Best next move without new data.</strong> Do not start with another ad hoc retraining round. First run nested leave-one-dataset-out evaluation with donor/library/clonotype grouping, expression-independent TCR labels, dataset-balanced weights, and all feature selection/calibration inside each training fold.</div>
</section>

<section>
  <h2>Independent reviewer result</h2>
  <p>A separately instantiated read-only reviewer examined the repository without editing files. Its conclusions were then reconciled against the main-agent provenance inventory.</p>
  <ul>
    <li>It independently confirmed both critical defects: V2 validation-based algorithm selection and BALF_BLOOD_COPD reuse in V3 promotion.</li>
    <li>It independently reproduced the corrected TCR-only result: Round 12 F1 0.8929 and Round 14 F1 0.8922, with Round 14 gaining 10 TP but adding 14 strict-abT FP.</li>
    <li>It identified the mixed-run V3 release record, full-atlas candidate selection, extension holdout unsealing, missing gene-completeness guardrails, and absence of grouped uncertainty as additional confirmed issues.</li>
    <li>It agreed that no currently inspected cohort can support a new claim of untouched external validation.</li>
  </ul>
</section>

<section>
  <h2>Where evaluation independence was lost</h2>
  <div class="figure"><img src="figures/evaluation_reuse_map.png" alt="Evaluation reuse map"><div class="caption">Code-level provenance shows selection reuse. This is selection leakage, not evidence that external labels were directly used to fit model coefficients.</div></div>
  <ul>
    <li>V2 thresholds were tuned on the tune split, but the final algorithm was selected by F1 on the cohort named validation.</li>
    <li>V3 Round 14 wraps the V2 logistic score with a deterministic TRDC/NK gate; it is not a newly fitted 210-gene classifier.</li>
    <li>Round 12 versus Round 14 promotion used atlas, held-out, and BALF_BLOOD_COPD F1 plus external specificity guardrails.</li>
    <li>Train/tune rows were randomized at cell level within source-label groups, allowing donor/library/clonotype correlations across both sets.</li>
  </ul>
</section>

<section class="page-break">
  <h2>Corrected TCR-only re-evaluation</h2>
  <p>The existing BALF_BLOOD_COPD predictions were re-scored without changing model outputs. Primary positives are all cells labeled by the source workflow as paired productive TRG/TRD with no productive TRA/TRB, including 98 CD4-like cells previously excluded by expression. Primary negatives are strict paired productive TRA/TRB with no productive TRG/TRD. NK, B, and myeloid annotations are excluded from this primary table and retained conceptually as stress tests.</p>
  <div class="figure"><img src="figures/corrected_tcr_only_performance.png" alt="Corrected TCR-only performance"><div class="caption">Round 12 and Round 14 have nearly identical F1. Round 12 trades slightly lower recall for higher precision and fewer strict abT false positives.</div></div>
  {table_html(corrected, [('model','Model','text'),('n_positive','TCR gdT','int'),('n_negative','TCR abT','int'),('tp','TP','int'),('fp','FP','int'),('fn','FN','int'),('precision','Precision','pct'),('recall','Recall','pct'),('specificity','Specificity','pct'),('f1','F1','pct'),('pr_auc','PR-AUC','float3')])}
  <p><strong>Interpretation:</strong> Round 12 F1 is {round12['f1']:.4f}; Round 14 F1 is {round14['f1']:.4f}. The difference is too small, and the cohort too reused, to justify a general superiority claim. Round 14 gains {100*(round14['recall']-round12['recall']):.2f} percentage points of recall while adding {int(round14['fp']-round12['fp'])} strict-abT false positives.</p>
  <h3>Recall of TCR-confirmed positives excluded by the old CD4-like rule</h3>
  {table_html(subgroup[subgroup['subgroup'].str.contains('excluded')], [('model','Model','text'),('n_positive','N','int'),('true_positive','Detected','int'),('false_negative','Missed','int'),('recall','Recall','pct')])}
  <h3>Library-level stability</h3>
  {table_html(per_library[per_library['model'].isin(['V2 high-purity','V3 Round 12','V3 Round 14'])], [('library_id','Library','text'),('model','Model','text'),('n_positive','Positive','int'),('n_negative','Negative','int'),('fp','FP','int'),('precision','Precision','pct'),('recall','Recall','pct'),('f1','F1','pct')])}
</section>

<section class="page-break">
  <h2>Rare-cell deployment changes purity</h2>
  <p>Sensitivity and specificity were taken from the corrected TCR-only benchmark and projected to plausible deployment prevalence. These are scenarios, not empirical estimates for every dataset.</p>
  <div class="figure"><img src="figures/prevalence_ppv_scenarios.png" alt="Prevalence-aware PPV"><div class="caption">Even small specificity differences matter when gdT prevalence is below 1%.</div></div>
  <h3>Expected results at 1% true prevalence per million cells</h3>
  {table_html(scenario_one_pct, [('model','Model','text'),('expected_true_positive_per_million','Expected TP','int'),('expected_false_positive_per_million','Expected FP','int'),('expected_predicted_positive_per_million','Predicted positive','int'),('expected_ppv','Expected PPV','pct'),('expected_fdr','Expected FDR','pct')])}
</section>

<section>
  <h2>Cross-dataset negative-control challenge</h2>
  <p>The frozen models were recently screened on eight extension cohorts containing 758,135 T/NK cells. They provide useful abT and NK specificity stress tests but no gdT-positive truth, so they cannot support sensitivity or final model promotion.</p>
  <div class="figure"><img src="figures/extension_strict_nk_challenge.png" alt="Extension strict NK challenge"><div class="caption">Dataset-level variation is more informative than a single pooled specificity. Small strict-NK strata are omitted from the chart but retained in canonical tables.</div></div>
  {table_html(extension_pooled, [('model','Model','text'),('n_datasets','Datasets','int'),('total_cells','Cells','int'),('abT_n','abT controls','int'),('abT_fp','abT FP','int'),('abT_fpr','abT FPR','pct'),('strict_NK_n','Strict NK','int'),('strict_NK_fp','NK FP','int'),('strict_NK_fpr','NK FPR','pct')])}
</section>

<section class="page-break">
  <h2>Issue register</h2>
  <p>Findings are ordered by severity. “Confirmed” means the behavior is directly visible in code or retained artifacts; design risks are explicitly labeled as such.</p>
  <div class="table-wrap"><table class="issue-table"><thead><tr><th>Severity</th><th>ID / stage</th><th>Finding</th><th>Consequence</th><th>Required action</th></tr></thead><tbody>{''.join(issue_rows)}</tbody></table></div>
  <p class="small print-omit">Full evidence paths, certainty labels, and recommended actions are in <code>issue_register.csv</code>.</p>
</section>

<section class="page-break">
  <h2>What can be done without essentially more data</h2>
  <p>The existing data are sufficient for a much more honest and informative evaluation. The constraint is not sample count alone; it is how independent units, labels, feature selection, and adaptive decisions are handled.</p>
  <div class="table-wrap"><table><thead><tr><th>Priority</th><th>Action</th><th>Implementation</th><th>Why it matters</th></tr></thead><tbody>{recommendation_rows}</tbody></table></div>
  <h3>Proposed promotion rule for the next iteration</h3>
  <ul>
    <li>Primary objective: maximize dataset-macro F1 across outer folds, not pooled cell F1.</li>
    <li>Safety guardrails: report median and worst-dataset strict NK FPR, paired-abT FPR, and recall in every positive source.</li>
    <li>Reject a more complex model unless its outer-fold improvement is consistent and confidence intervals exclude a practically negligible gain.</li>
    <li>Freeze two operating points from out-of-fold predictions: highest F1 and high purity under a prespecified FPR/FDR target.</li>
    <li>After model choice, fit once on all eligible data. Do not present an internal resampling result as a new independent external validation.</li>
  </ul>
</section>

<section>
  <h2>Coefficient and feature assessment</h2>
  <p>The V2 release uses <code>StandardScaler</code> followed by class-balanced logistic regression. Its standardized coefficients range from {coefficient_summary.iloc[0]['minimum_coefficient']:.3f} to {coefficient_summary.iloc[0]['maximum_coefficient']:.3f}; the maximum absolute value is {coefficient_summary.iloc[0]['maximum_absolute_coefficient']:.3f}. That magnitude is not, by itself, evidence of instability. The stronger overfitting concerns are split leakage, source confounding, repeated model selection, and unmeasured feature stability across held-out datasets.</p>
  {table_html(top_coefficients, [('feature','Feature','text'),('standardized_coefficient','Standardized coefficient','float3'),('absolute_coefficient','Absolute value','float3')])}
  <p class="small">Positive TRBC/TRAC coefficients and large individual V/J effects should be treated as prompts for outer-fold stability analysis, not as biological conclusions.</p>
</section>

<section>
  <h2>Reproducibility and scope</h2>
  <p>All registered model artifacts passed SHA256 verification during this audit.</p>
  <div class="table-wrap"><table><thead><tr><th>Model</th><th>Status</th><th>Bytes</th><th>SHA256 prefix</th><th>Check</th></tr></thead><tbody>{artifact_rows}</tbody></table></div>
  <h3>V3 semantic provenance mismatch</h3>
  <ul>
    <li>Manifest model/status: <code>{html.escape(str(v3_provenance.iloc[0]['manifest_model_name']))}</code> / <code>{html.escape(str(v3_provenance.iloc[0]['manifest_status']))}</code>.</li>
    <li>Embedded pickle model/version/acceptance: <code>{html.escape(str(v3_provenance.iloc[0]['payload_model']))}</code> / <code>{html.escape(str(v3_provenance.iloc[0]['payload_version']))}</code> / <code>{html.escape(str(v3_provenance.iloc[0]['payload_accepted_for_promotion']))}</code>.</li>
    <li>Completed-run summary selected <code>{html.escape(str(v3_provenance.iloc[0]['completed_run_selected_candidate']))}</code>, while the canonical path now contains Round 14.</li>
    <li>Recorded training negatives: {fmt_int(v3_provenance.iloc[0]['recorded_training_negatives'])}; methodology statement: {fmt_int(v3_provenance.iloc[0]['methodology_documented_training_negatives'])}.</li>
  </ul>
  <p><strong>Still missing for exact reconstruction:</strong> immutable input-data hashes, exact split cell IDs, grouped fold assignments, the Git commit for every historical run, a package/environment lock, and fold-local feature-selection records.</p>
  <p><strong>Audit limitation:</strong> no unused cohort remains after repeated iteration. The corrected TCR-only analysis improves label independence but does not restore test independence. Nested grouped evaluation is therefore the strongest defensible next estimate available without new data.</p>
</section>

<footer>Canonical tables: Integrated_dataset/tables/gdT_prediction/gdtai_methodology_audit/ &nbsp;|&nbsp; Canonical figures: Integrated_dataset/figures/gdT_prediction/gdtai_methodology_audit/ &nbsp;|&nbsp; Summary: Integrated_dataset/logs/gdT_prediction/gdtai_methodology_audit/</footer>
</main></body></html>"""
    REPORT_HTML.write_text(html_text, encoding="utf-8")


def render_summary(
    issues: pd.DataFrame,
    corrected: pd.DataFrame,
    extension_pooled: pd.DataFrame,
    coefficient_summary: pd.DataFrame,
    v3_provenance: pd.DataFrame,
) -> None:
    corrected_index = corrected.set_index("model")
    extension_index = extension_pooled.set_index("model")
    summary = f"""# gdTAI Training And Evaluation Audit Summary

Date: {datetime.now().astimezone().strftime('%Y-%m-%d %Z')}

## Conclusion

gdTAI contains useful biological signal, but current performance estimates do not rule out selection overfitting. The highest-priority issues are:

1. V2 algorithm selection used the cohort named validation.
2. V3 Round 12 versus Round 14 promotion used the cohort described as independent external.
3. Train/tune splitting was at cell level within source-label strata rather than donor/library/clonotype groups.
4. The external positive label used CD4 annotation/expression exclusion while CD4 is a model feature.
5. Positive/negative labels and held-out components are strongly source-confounded.
6. The V3 pickle, manifest, completed-run summary, and documented training composition do not describe one coherent promoted build record.

## Corrected TCR-only benchmark

The existing external predictions were re-evaluated using 1,046 paired productive gdT/no productive abT positives and 33,117 strict paired productive abT/no productive gdT negatives.

| Model | Precision | Recall | Specificity | F1 | FP |
| --- | ---: | ---: | ---: | ---: | ---: |
| V2 high-purity | {corrected_index.loc['V2 high-purity','precision']:.4f} | {corrected_index.loc['V2 high-purity','recall']:.4f} | {corrected_index.loc['V2 high-purity','specificity']:.4f} | {corrected_index.loc['V2 high-purity','f1']:.4f} | {int(corrected_index.loc['V2 high-purity','fp'])} |
| V3 Round 12 | {corrected_index.loc['V3 Round 12','precision']:.4f} | {corrected_index.loc['V3 Round 12','recall']:.4f} | {corrected_index.loc['V3 Round 12','specificity']:.4f} | {corrected_index.loc['V3 Round 12','f1']:.4f} | {int(corrected_index.loc['V3 Round 12','fp'])} |
| V3 Round 14 | {corrected_index.loc['V3 Round 14','precision']:.4f} | {corrected_index.loc['V3 Round 14','recall']:.4f} | {corrected_index.loc['V3 Round 14','specificity']:.4f} | {corrected_index.loc['V3 Round 14','f1']:.4f} | {int(corrected_index.loc['V3 Round 14','fp'])} |

This benchmark is still reused evidence, not an untouched external test. Round 12 and Round 14 are effectively tied on F1; Round 12 is purer and Round 14 is slightly more sensitive.

## Coefficients

The V2 standardized logistic coefficients have maximum absolute magnitude {coefficient_summary.iloc[0]['maximum_absolute_coefficient']:.3f}. Coefficient magnitude is not the primary overfitting concern.

The V3 payload embeds `accepted_for_promotion={v3_provenance.iloc[0]['payload_accepted_for_promotion']}` while the manifest says promoted; the retained completed-run summary selects `{v3_provenance.iloc[0]['completed_run_selected_candidate']}`. This is a semantic provenance defect, not a checksum failure.

## Extension challenge

- V2 high-purity pooled strict-NK FPR: {extension_index.loc['V2 high-purity','strict_NK_fpr']:.4%}
- V3 Round 12 pooled strict-NK FPR: {extension_index.loc['V3 Round 12','strict_NK_fpr']:.4%}
- V3 Round 14 pooled strict-NK FPR: {extension_index.loc['V3 Round 14','strict_NK_fpr']:.4%}

These eight extension cohorts contain no positive truth and therefore cannot select the final model.

## Recommended next action

Run a precommitted nested leave-one-dataset-out re-evaluation. Use donor/library/clonotype-grouped inner folds; pure TCR labels independent of model features; dataset-balanced weights; fold-local feature selection/calibration; macro-dataset metrics; realistic-prevalence PPV/FDR; and explicit NK/abT/dropout stress panels. Keep V2, Round 12, and Round 14 frozen as comparators.

## Outputs

- HTML: `gdT_prediction/gdtai_methodology_audit/index.html`
- PDF: `gdT_prediction/gdtai_methodology_audit/gdtai_methodology_audit_report.pdf`
- Independent reviewer record: `Integrated_dataset/logs/gdT_prediction/gdtai_methodology_audit/independent_reviewer_record.md`
- Issue register: `Integrated_dataset/tables/gdT_prediction/gdtai_methodology_audit/issue_register.csv`
- Recommended actions: `Integrated_dataset/tables/gdT_prediction/gdtai_methodology_audit/recommended_actions.csv`
"""
    SUMMARY_MD.write_text(summary, encoding="utf-8")


def copy_figures(paths: list[Path]) -> None:
    for path in paths:
        shutil.copy2(path, STATIC_FIGURE_DIR / path.name)


def export_pdf() -> None:
    chrome = shutil.which("google-chrome") or shutil.which("google-chrome-stable") or shutil.which("chromium")
    if not chrome:
        raise FileNotFoundError("Chrome/Chromium is required for PDF export.")
    profile_dir = Path("/tmp/gdtai-methodology-audit-chrome-profile")
    profile_dir.mkdir(parents=True, exist_ok=True)
    temporary_pdf = Path("/tmp") / f"gdtai-methodology-audit-render-{os.getpid()}.pdf"
    command = [
        chrome,
        "--headless",
        "--disable-gpu",
        "--disable-dev-shm-usage",
        "--disable-breakpad",
        "--disable-crash-reporter",
        "--no-sandbox",
        "--allow-file-access-from-files",
        "--no-pdf-header-footer",
        f"--user-data-dir={profile_dir}",
        f"--print-to-pdf={temporary_pdf}",
        REPORT_HTML.resolve().as_uri(),
    ]
    last_error: subprocess.CalledProcessError | None = None
    for _ in range(2):
        try:
            subprocess.run(
                command,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            last_error = None
            break
        except subprocess.CalledProcessError as exc:
            last_error = exc
    if last_error is not None:
        raise RuntimeError(f"Chrome PDF export failed twice: {last_error.stderr}") from last_error
    shutil.copy2(temporary_pdf, REPORT_PDF)
    if not REPORT_PDF.exists() or REPORT_PDF.stat().st_size < 10_000:
        raise RuntimeError("PDF export did not produce a plausible report file.")


def write_run_manifest(outputs: list[Path], artifacts: pd.DataFrame) -> None:
    manifest = {
        "generated_at": datetime.now().astimezone().isoformat(),
        "workflow": str(Path(__file__).resolve().relative_to(ROOT)),
        "read_only_no_h5ad": True,
        "inputs": {
            str(path.relative_to(ROOT)): {"sha256": sha256(path), "size_bytes": path.stat().st_size}
            for path in (
                REGISTRY_PATH,
                V2_MODEL_PATH,
                V3_MODEL_PATH,
                V3_MANIFEST_PATH,
                V3_RUN_SUMMARY_PATH,
                V3_TRAINING_COMPOSITION_PATH,
                V2_SPLIT_PATH,
                V2_VALIDATION_PATH,
                V3_GUARDRAIL_PATH,
                EXTERNAL_PREDICTIONS_PATH,
                EXTENSION_SCREEN_PATH,
            )
        },
        "registered_artifacts_verified": artifacts[["model_id", "actual_sha256", "checksum_match"]].to_dict(orient="records"),
        "outputs": [str(path.relative_to(ROOT)) for path in outputs if path.exists()],
    }
    (LOG_DIR / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main() -> None:
    args = parse_args()
    ensure_directories()
    artifacts = audit_model_artifacts()
    coefficient_table, coefficient_summary = inspect_v2_coefficients()
    v3_provenance = audit_v3_semantic_provenance()
    corrected, per_library, subgroup = corrected_external_benchmark()
    scenarios = prevalence_scenarios(corrected)
    extension_data, extension_pooled = extension_negative_summary()
    issues = build_issue_register()
    recommendations = build_recommendation_register()
    figures = [
        plot_corrected_performance(corrected),
        plot_prevalence_ppv(scenarios),
        plot_extension_challenge(extension_data, extension_pooled),
        plot_evaluation_reuse_map(),
    ]
    copy_figures(figures)
    render_html(
        issues,
        recommendations,
        corrected,
        per_library,
        subgroup,
        scenarios,
        extension_pooled,
        coefficient_table,
        coefficient_summary,
        artifacts,
        v3_provenance,
    )
    render_summary(issues, corrected, extension_pooled, coefficient_summary, v3_provenance)
    if not args.no_pdf:
        export_pdf()
    outputs = [
        REPORT_HTML,
        REPORT_PDF,
        SUMMARY_MD,
        TABLE_DIR / "issue_register.csv",
        TABLE_DIR / "recommended_actions.csv",
        TABLE_DIR / "corrected_external_pure_tcr_metrics.csv",
        TABLE_DIR / "prevalence_ppv_scenarios.csv",
        *figures,
    ]
    write_run_manifest(outputs, artifacts)
    print(f"Wrote {REPORT_HTML.relative_to(ROOT)}")
    if not args.no_pdf and REPORT_PDF.exists():
        print(f"Wrote {REPORT_PDF.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
