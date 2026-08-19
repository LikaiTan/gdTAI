#!/usr/bin/env python3
"""Score frozen gdTAI V4.3, V3, and V2 on one common external lockbox."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import pickle
import sys
import time
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
import xgboost as xgb


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "workflows/gdtai") not in sys.path:
    sys.path.insert(0, str(ROOT / "workflows/gdtai"))

from run_gdtai_v4_2_ground_truth import axis_node_values  # noqa: E402
from prepare_gdtai_v4_2_training import extract_csr_rows, matrix_group, read_var_names  # noqa: E402
from run_gdtai_v3_trdc_nk_guard_classifier import (  # noqa: E402
    ENGINEERED_FEATURES,
    FeatureSpec,
    append_engineered_features,
)
from train_gdtai_v4_2_nested import (  # noqa: E402
    exclusion_flags,
    receptor_rescue_flags,
    sha256_file,
)


ATLAS = Path("/ssd/tnk_phase3/Integrated_dataset/full_atlas/tcr_corrected/integrated_full_atlas.h5ad")
LOCKBOX = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_3_lockbox/lockbox_manifest.parquet"
LOCKBOX_FREEZE = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_3_lockbox/freeze_summary.json"
V4_DIR = ROOT / "Integrated_dataset/models/gdT_prediction/gdtai_v4_3_receptor_isolated/final_development"
V4_CONTRACT = V4_DIR / "model_contract.json"
V2_MODEL = ROOT / "Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v2.0/gdTAI_v2_model.pkl"
V3_MODEL = ROOT / "Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/gdTAI_v3_model.pkl"
CONFIG = ROOT / "configs/models/gdtai/v4_3_rescue_training.json"
CACHE_DIR = ROOT / "Integrated_dataset/cache/gdT_prediction/gdtai_v4_3_common_lockbox"
TABLE_DIR = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_3_common_lockbox"
LOG_DIR = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_3_common_lockbox"
FEATURE_CACHE = CACHE_DIR / "lockbox_gene_features.npy"
FEATURE_MANIFEST = CACHE_DIR / "lockbox_feature_manifest.csv"
FEATURE_SUMMARY = CACHE_DIR / "lockbox_feature_cache_summary.json"
FINAL_SUMMARY = LOG_DIR / "evaluation_summary.json"
MODEL_ORDER = ["v4_3_balanced", "v3_balanced", "v2_high_f1", "v2_high_purity"]


def load_pickle(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        value = pickle.load(handle)
    if not isinstance(value, dict):
        raise TypeError(f"Expected dict payload in {path}")
    return value


def lockbox_contract() -> tuple[pd.DataFrame, dict, dict, dict[str, Any], dict[str, Any]]:
    freeze = json.loads(LOCKBOX_FREEZE.read_text())
    if freeze["manifest_sha256"] != sha256_file(LOCKBOX):
        raise RuntimeError("Lockbox manifest differs from its pre-score freeze")
    if freeze["model_scored"] or freeze["threshold_selected"]:
        raise RuntimeError("Lockbox freeze does not have pre-score status")
    v4 = json.loads(V4_CONTRACT.read_text())
    if not v4["development_frozen"] or v4["lockbox_scored"] or v4["model_promoted"]:
        raise RuntimeError("V4 final-development contract is not eligible for first lockbox scoring")
    config = json.loads(CONFIG.read_text())
    v2 = load_pickle(V2_MODEL)
    v3 = load_pickle(V3_MODEL)
    lockbox = pd.read_parquet(LOCKBOX).reset_index(drop=True)
    if lockbox.cell_id.duplicated().any() or lockbox.atlas_row.duplicated().any():
        raise RuntimeError("Lockbox rows are not unique")
    if lockbox.allow_fit.any() or lockbox.allow_threshold_selection.any():
        raise RuntimeError("Lockbox contains a fitting or threshold-selection row")
    return lockbox, freeze, v4, v2, v3


def union_genes(v4: dict, v2: dict, v3: dict) -> list[str]:
    ordered: list[str] = []
    for gene in [
        *v4["stage1_feature_names"],
        *v4["stage2_feature_names"],
        *v2["base_model"]["gene_names"],
        *v3["gene_names"],
    ]:
        value = str(gene)
        if value not in ordered:
            ordered.append(value)
    return ordered


def prepare_features(row_chunk: int) -> dict:
    lockbox, freeze, v4, v2, v3 = lockbox_contract()
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    genes = union_genes(v4, v2, v3)
    with h5py.File(ATLAS, "r") as handle:
        var_names = read_var_names(handle)
        lookup = {gene: index for index, gene in enumerate(var_names)}
        missing = sorted(set(genes) - set(lookup))
        if missing:
            raise RuntimeError(f"Common-lockbox genes missing from atlas: {missing}")
        mapping = np.full(len(var_names), -1, dtype=np.int32)
        for output_index, gene in enumerate(genes):
            mapping[lookup[gene]] = output_index
        matrix = extract_csr_rows(
            matrix_group(handle, "X"),
            lockbox.atlas_row.to_numpy(np.int64),
            mapping,
            len(genes),
            row_chunk,
        )
    temporary = FEATURE_CACHE.with_suffix(".tmp.npy")
    if temporary.exists():
        temporary.unlink()
    np.save(temporary, matrix)
    os.replace(temporary, FEATURE_CACHE)
    manifest = pd.DataFrame({"feature_index": np.arange(len(genes)), "gene": genes})
    manifest.to_csv(FEATURE_MANIFEST, index=False)
    summary = {
        "status": "PASS_COMMON_LOCKBOX_FEATURE_CACHE",
        "model_scored": False,
        "n_cells": len(lockbox),
        "n_genes": len(genes),
        "normalization": "log1p(raw_counts_per_10000)",
        "atlas_sha256": "d32c9d2bdb955b12e1eafbed8322f8cb965cf3a225191e612b53f3d3783480d5",
        "lockbox_sha256": freeze["manifest_sha256"],
        "v4_contract_sha256": sha256_file(V4_CONTRACT),
        "v2_model_sha256": sha256_file(V2_MODEL),
        "v3_model_sha256": sha256_file(V3_MODEL),
        "feature_manifest_sha256": sha256_file(FEATURE_MANIFEST),
        "feature_cache_sha256": sha256_file(FEATURE_CACHE),
    }
    FEATURE_SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def verify_feature_cache(lockbox: pd.DataFrame, freeze: dict, v4: dict, v2: dict, v3: dict) -> tuple[np.ndarray, list[str]]:
    summary = json.loads(FEATURE_SUMMARY.read_text())
    expected = {
        "lockbox_sha256": freeze["manifest_sha256"],
        "v4_contract_sha256": sha256_file(V4_CONTRACT),
        "v2_model_sha256": sha256_file(V2_MODEL),
        "v3_model_sha256": sha256_file(V3_MODEL),
        "feature_manifest_sha256": sha256_file(FEATURE_MANIFEST),
        "feature_cache_sha256": sha256_file(FEATURE_CACHE),
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            raise RuntimeError(f"Feature-cache contract mismatch: {key}")
    genes = pd.read_csv(FEATURE_MANIFEST).sort_values("feature_index").gene.astype(str).tolist()
    if genes != union_genes(v4, v2, v3):
        raise RuntimeError("Feature-cache gene order differs from model union")
    matrix = np.load(FEATURE_CACHE, mmap_mode="r")
    if matrix.shape != (len(lockbox), len(genes)) or not np.isfinite(matrix).all():
        raise RuntimeError("Common-lockbox feature matrix is invalid")
    return matrix, genes


def subset_annotation(rows: np.ndarray) -> np.ndarray:
    rows = np.asarray(rows, dtype=np.int64)
    with h5py.File(ATLAS, "r") as handle:
        node = handle["obs"]["source_cell_type"]
        if not isinstance(node, h5py.Group) or not {"categories", "codes"}.issubset(node):
            raise TypeError("Expected categorical source_cell_type in the frozen atlas")
        categories = axis_node_values(node["categories"]).astype(str)
        codes = np.asarray(node["codes"][:], dtype=np.int64)[rows]
        values = np.full(len(rows), "", dtype=object)
        present = codes >= 0
        values[present] = categories[codes[present]]
    normalized = pd.Series(values).str.upper().str.strip()
    out = np.full(len(values), "OTHER", dtype=object)
    out[normalized.str.startswith("NK").to_numpy()] = "NK"
    out[normalized.str.contains("GDT|GAMMA|DELTA", regex=True).to_numpy()] = "GDT"
    out[normalized.str.startswith("CD8").to_numpy()] = "CD8_T"
    out[normalized.str.startswith("CD4").to_numpy()] = "CD4_T"
    out[normalized.str.contains("TREG|REGULATORY").to_numpy()] = "TREG"
    return out


def v2_thresholds(payload: dict, mode: str, annotation: np.ndarray) -> np.ndarray:
    specification = payload["operating_modes"][mode]
    if mode == "high_f1":
        return np.full(len(annotation), float(specification["threshold"]), dtype=np.float32)
    values = specification["annotation_thresholds"]
    output = np.full(len(annotation), float(values["other_threshold"]), dtype=np.float32)
    mapping = {
        "GDT": "gdt_threshold",
        "CD8_T": "cd8_threshold",
        "CD4_T": "cd4_threshold",
        "TREG": "treg_threshold",
        "NK": "nk_threshold",
    }
    for annotation_name, threshold_name in mapping.items():
        mask = annotation == annotation_name
        value = values[threshold_name]
        output[mask] = np.inf if not np.isfinite(value) else float(value)
    return output


def v3_feature_spec(v3: dict) -> FeatureSpec:
    genes = [str(gene) for gene in v3["gene_names"]]
    engineered = [str(value) for value in v3.get("engineered_feature_names", ENGINEERED_FEATURES)]
    return FeatureSpec(
        gene_names=genes,
        gene_indices=np.arange(len(genes), dtype=np.int32),
        gene_feature_names=[f"{gene}_log1p_cp10k" for gene in genes],
        engineered_feature_names=engineered,
        model_feature_names=[*[f"{gene}_log1p_cp10k" for gene in genes], *engineered],
        gene_to_col={gene: index for index, gene in enumerate(genes)},
        engineered_to_col={name: len(genes) + index for index, name in enumerate(engineered)},
    )


def v3_columns(spec: FeatureSpec, v3: dict) -> np.ndarray:
    output = []
    suffix = "_log1p_cp10k"
    for feature in map(str, v3["feature_names"]):
        if feature.endswith(suffix):
            output.append(spec.gene_to_col[feature[: -len(suffix)]])
        else:
            output.append(spec.engineered_to_col[feature])
    return np.asarray(output, dtype=np.int64)


def predict_xgb_cpu(model: xgb.XGBClassifier, matrix: np.ndarray, chunk: int = 100_000) -> np.ndarray:
    output = []
    model.set_params(device="cpu")
    for start in range(0, len(matrix), chunk):
        output.append(
            model.predict_proba(np.asarray(matrix[start:start + chunk], dtype=np.float32))[:, 1]
        )
    return np.concatenate(output).astype(np.float32)


def score_models(matrix: np.ndarray, genes: list[str], lockbox: pd.DataFrame, v4: dict, v2: dict, v3: dict) -> pd.DataFrame:
    lookup = {gene: index for index, gene in enumerate(genes)}
    stage1_columns = np.asarray([lookup[gene] for gene in v4["stage1_feature_names"]], dtype=np.int64)
    stage2_columns = np.asarray([lookup[gene] for gene in v4["stage2_feature_names"]], dtype=np.int64)
    stage1 = xgb.XGBClassifier()
    stage1.load_model(V4_DIR / "stage1_t_lineage_gate.ubj")
    stage2 = xgb.XGBClassifier()
    stage2.load_model(V4_DIR / "stage2_receptor_classifier.ubj")
    stage1_score = predict_xgb_cpu(stage1, np.asarray(matrix[:, stage1_columns]))
    stage2_score = predict_xgb_cpu(stage2, np.asarray(matrix[:, stage2_columns]))
    excluded = exclusion_flags(np.asarray(matrix), genes)[2]
    rescue = receptor_rescue_flags(np.asarray(matrix), genes, v4["receptor_rescue"])
    v4_call = (
        (stage1_score >= float(v4["stage1_threshold"]))
        & ((stage2_score >= float(v4["stage2_threshold"])) | rescue)
        & ~excluded
    )

    v2_columns = np.asarray([lookup[str(gene)] for gene in v2["base_model"]["gene_names"]], dtype=np.int64)
    v2_score = v2["base_model"]["model_object"].predict_proba(np.asarray(matrix[:, v2_columns]))[:, 1].astype(np.float32)
    annotation = subset_annotation(lockbox.atlas_row.to_numpy(np.int64))
    v2_high_f1_threshold = v2_thresholds(v2, "high_f1", annotation)
    v2_high_purity_threshold = v2_thresholds(v2, "high_purity", annotation)

    v3_gene_columns = np.asarray([lookup[str(gene)] for gene in v3["gene_names"]], dtype=np.int64)
    spec = v3_feature_spec(v3)
    engineered = append_engineered_features(np.asarray(matrix[:, v3_gene_columns]), spec)
    v3_score = v3["model_object"].predict_proba(engineered[:, v3_columns(spec, v3)])[:, 1].astype(np.float32)
    v3_call = v3_score >= float(v3["threshold"])

    output = lockbox[
        ["cell_id", "atlas_row", "truth_class", "source_gse_id", "donor_id", "sample_id", "group_id"]
    ].copy()
    output["annotation_for_v2"] = annotation
    output["v4_3_stage1_score"] = stage1_score
    output["v4_3_stage2_score"] = stage2_score
    output["v4_3_receptor_rescue"] = rescue
    output["v4_3_cd4_treg_excluded"] = excluded
    output["v4_3_balanced_score"] = stage1_score * stage2_score
    output["v4_3_balanced"] = v4_call
    output["v3_balanced_score"] = v3_score
    output["v3_balanced"] = v3_call
    output["v2_score"] = v2_score
    output["v2_high_f1_threshold"] = v2_high_f1_threshold
    output["v2_high_f1"] = v2_score >= v2_high_f1_threshold
    output["v2_high_purity_threshold"] = v2_high_purity_threshold
    output["v2_high_purity"] = v2_score >= v2_high_purity_threshold
    return output


def confusion_metrics(frame: pd.DataFrame, model: str) -> dict[str, float | int | str]:
    call = frame[model].to_numpy(bool)
    positive = frame.truth_class.eq("gdT_gold").to_numpy()
    abt = frame.truth_class.eq("abT_gold").to_numpy()
    nk = frame.truth_class.eq("nk_lockbox").to_numpy()
    negative = abt | nk
    tp = int((call & positive).sum())
    fn = int((~call & positive).sum())
    fp = int((call & negative).sum())
    tn = int((~call & negative).sum())
    recall = tp / (tp + fn)
    precision = tp / (tp + fp) if tp + fp else 0.0
    specificity = tn / (tn + fp)
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    abt_fpr = float(call[abt].mean())
    nk_fpr = float(call[nk].mean())
    result: dict[str, float | int | str] = {
        "model": model,
        "tp": tp,
        "fp": fp,
        "tn": tn,
        "fn": fn,
        "predicted_positive": int(call.sum()),
        "precision": precision,
        "recall": recall,
        "specificity": specificity,
        "f1": f1,
        "abt_fpr": abt_fpr,
        "author_nk_fpr": nk_fpr,
    }
    combined_fpr = 1.0 - specificity
    for prevalence in (0.01, 0.03, 0.05):
        adjusted_precision = (recall * prevalence) / (
            recall * prevalence + combined_fpr * (1.0 - prevalence)
        ) if recall > 0 or combined_fpr > 0 else 0.0
        result[f"precision_at_{int(prevalence * 100)}pct_prevalence"] = adjusted_precision
    return result


def per_source_metrics(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (truth, source), group in frame.groupby(["truth_class", "source_gse_id"], sort=True):
        for model in MODEL_ORDER:
            calls = group[model].to_numpy(bool)
            rows.append({
                "truth_class": truth,
                "source_gse_id": source,
                "model": model,
                "n_cells": len(group),
                "predicted_positive": int(calls.sum()),
                "positive_rate": float(calls.mean()),
                "metric_name": "recall" if truth == "gdT_gold" else "false_positive_rate",
            })
    return pd.DataFrame(rows)


def overlap_table(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for older in MODEL_ORDER[1:]:
        v4 = frame.v4_3_balanced.to_numpy(bool)
        old = frame[older].to_numpy(bool)
        for truth, mask in frame.groupby("truth_class").groups.items():
            idx = np.asarray(mask, dtype=np.int64)
            rows.append({
                "older_model": older,
                "truth_class": truth,
                "n_cells": len(idx),
                "both": int((v4[idx] & old[idx]).sum()),
                "v4_only": int((v4[idx] & ~old[idx]).sum()),
                "older_only": int((~v4[idx] & old[idx]).sum()),
                "neither": int((~v4[idx] & ~old[idx]).sum()),
            })
    return pd.DataFrame(rows)


def bootstrap_units(frame: pd.DataFrame) -> pd.DataFrame:
    donor = frame.donor_id.fillna("").astype(str)
    fallback = frame.group_id.fillna("").astype(str)
    unit = np.where(donor.str.strip().ne(""), frame.source_gse_id.astype(str) + "::" + donor, fallback)
    local = frame[["truth_class", "source_gse_id"]].copy()
    local["bootstrap_unit"] = unit
    for model in MODEL_ORDER:
        local[model] = frame[model].astype(np.int64)
    aggregated = local.groupby(
        ["truth_class", "source_gse_id", "bootstrap_unit"], as_index=False, observed=True
    ).agg(n_cells=("bootstrap_unit", "size"), **{model: (model, "sum") for model in MODEL_ORDER})
    return aggregated


def counts_to_metrics(counts: dict[str, dict[str, int]], model: str) -> dict[str, float]:
    tp = counts["gdT_gold"][model]
    fn = counts["gdT_gold"]["n_cells"] - tp
    ab_fp = counts["abT_gold"][model]
    nk_fp = counts["nk_lockbox"][model]
    fp = ab_fp + nk_fp
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    return {
        "f1": 2 * precision * recall / (precision + recall) if precision + recall else 0.0,
        "recall": recall,
        "abt_fpr": ab_fp / counts["abT_gold"]["n_cells"],
        "author_nk_fpr": nk_fp / counts["nk_lockbox"]["n_cells"],
    }


def cluster_bootstrap(frame: pd.DataFrame, iterations: int, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    units = bootstrap_units(frame)
    strata = [group.reset_index(drop=True) for _, group in units.groupby(["truth_class", "source_gse_id"], sort=True)]
    rng = np.random.default_rng(seed)
    draws = []
    for iteration in range(iterations):
        sampled = []
        for stratum in strata:
            indexes = rng.integers(0, len(stratum), size=len(stratum))
            sampled.append(stratum.iloc[indexes])
        sample = pd.concat(sampled, ignore_index=True)
        counts: dict[str, dict[str, int]] = {}
        for truth, group in sample.groupby("truth_class", sort=False):
            counts[str(truth)] = {
                "n_cells": int(group.n_cells.sum()),
                **{model: int(group[model].sum()) for model in MODEL_ORDER},
            }
        metrics = {model: counts_to_metrics(counts, model) for model in MODEL_ORDER}
        for older in MODEL_ORDER[1:]:
            draws.append({
                "iteration": iteration,
                "older_model": older,
                **{f"delta_{metric}": metrics["v4_3_balanced"][metric] - metrics[older][metric]
                   for metric in ("f1", "recall", "abt_fpr", "author_nk_fpr")},
            })
    draw_table = pd.DataFrame(draws)
    summary_rows = []
    for older, group in draw_table.groupby("older_model", sort=False):
        row: dict[str, Any] = {"older_model": older, "iterations": iterations}
        for metric in ("f1", "recall", "abt_fpr", "author_nk_fpr"):
            values = group[f"delta_{metric}"].to_numpy(float)
            row[f"delta_{metric}_median"] = float(np.median(values))
            row[f"delta_{metric}_ci_low"] = float(np.quantile(values, 0.025))
            row[f"delta_{metric}_ci_high"] = float(np.quantile(values, 0.975))
            row[f"probability_delta_{metric}_gt_0"] = float((values > 0).mean())
        summary_rows.append(row)
    return draw_table, pd.DataFrame(summary_rows)


def promotion_decision(overall: pd.DataFrame, per_source: pd.DataFrame, bootstrap: pd.DataFrame, config: dict) -> dict:
    rule = config["common_lockbox_decision"]
    metrics = overall.set_index("model")
    positive = per_source[per_source.truth_class.eq("gdT_gold")]
    pivot = positive.pivot(index="source_gse_id", columns="model", values="positive_rate")
    older = ["v3_balanced", "v2_high_f1", "v2_high_purity"]
    deficits = pivot[older].max(axis=1) - pivot["v4_3_balanced"]
    v3_bootstrap = bootstrap.set_index("older_model").loc["v3_balanced"]
    checks = {
        "v4_abt_fpr_within_limit": bool(metrics.loc["v4_3_balanced", "abt_fpr"] <= rule["maximum_v4_abt_fpr"]),
        "v4_author_nk_fpr_within_limit": bool(metrics.loc["v4_3_balanced", "author_nk_fpr"] <= rule["maximum_v4_author_nk_fpr"]),
        "no_source_recall_deficit_over_limit": bool(deficits.max() <= rule["maximum_source_recall_deficit_vs_best_older"]),
        "cluster_bootstrap_f1_improvement_vs_v3": bool(
            v3_bootstrap.probability_delta_f1_gt_0 >= rule["minimum_cluster_bootstrap_probability_f1_improvement_vs_v3"]
        ),
        "observed_f1_above_v3": bool(metrics.loc["v4_3_balanced", "f1"] > metrics.loc["v3_balanced", "f1"]),
    }
    passed = all(checks.values())
    return {
        "status": "PASS_V4_SURPASSES_OLDER" if passed else "FAIL_V4_NOT_SUPERIOR",
        "eligible_for_promotion_review": passed,
        "model_promoted": False,
        "no_lockbox_retuning": True,
        "checks": checks,
        "maximum_observed_source_recall_deficit_vs_best_older": float(deficits.max()),
        "source_recall_deficits": {str(key): float(value) for key, value in deficits.items()},
    }


def score(iterations: int, seed: int) -> dict:
    if FINAL_SUMMARY.exists():
        raise RuntimeError(f"Common lockbox was already scored: {FINAL_SUMMARY}")
    lockbox, freeze, v4, v2, v3 = lockbox_contract()
    matrix, genes = verify_feature_cache(lockbox, freeze, v4, v2, v3)
    config = json.loads(CONFIG.read_text())
    started = time.monotonic()
    predictions = score_models(matrix, genes, lockbox, v4, v2, v3)
    overall = pd.DataFrame([confusion_metrics(predictions, model) for model in MODEL_ORDER])
    source = per_source_metrics(predictions)
    overlap = overlap_table(predictions)
    draws, bootstrap = cluster_bootstrap(predictions, iterations, seed)
    decision = promotion_decision(overall, source, bootstrap, config)

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    predictions_path = TABLE_DIR / "common_lockbox_predictions.parquet"
    predictions.to_parquet(predictions_path, index=False, compression="zstd")
    overall.to_csv(TABLE_DIR / "overall_metrics.csv", index=False)
    source.to_csv(TABLE_DIR / "per_source_metrics.csv", index=False)
    overlap.to_csv(TABLE_DIR / "prediction_overlap.csv", index=False)
    draws.to_parquet(TABLE_DIR / "cluster_bootstrap_draws.parquet", index=False, compression="zstd")
    bootstrap.to_csv(TABLE_DIR / "cluster_bootstrap_summary.csv", index=False)
    (LOG_DIR / "promotion_decision.json").write_text(json.dumps(decision, indent=2, sort_keys=True) + "\n")
    summary = {
        "status": "PASS_COMMON_LOCKBOX_EVALUATION",
        "decision": decision,
        "runtime_seconds": time.monotonic() - started,
        "lockbox_scored": True,
        "thresholds_retuned": False,
        "model_promoted": False,
        "n_lockbox_cells": len(lockbox),
        "lockbox_sha256": freeze["manifest_sha256"],
        "v4_contract_sha256": sha256_file(V4_CONTRACT),
        "predictions_sha256": sha256_file(predictions_path),
        "overall_metrics": overall.to_dict(orient="records"),
    }
    temporary = FINAL_SUMMARY.with_suffix(".tmp.json")
    temporary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, FINAL_SUMMARY)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["prepare", "score", "all"], default="all")
    parser.add_argument("--row-chunk", type=int, default=10_000)
    parser.add_argument("--bootstrap-iterations", type=int, default=2_000)
    parser.add_argument("--seed", type=int, default=20260819)
    args = parser.parse_args()
    result = {}
    if args.stage in {"prepare", "all"}:
        result["prepare"] = prepare_features(args.row_chunk)
        print(json.dumps(result["prepare"], indent=2, sort_keys=True))
    if args.stage in {"score", "all"}:
        result["score"] = score(args.bootstrap_iterations, args.seed)


if __name__ == "__main__":
    main()
