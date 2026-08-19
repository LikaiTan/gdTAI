#!/usr/bin/env python3
"""Audit the dominant GSE159251 false-positive failure of frozen gdTAI V4.5."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb

from train_gdtai_v4_2_nested import sha256_file
from train_gdtai_v4_4_dual_mode import compose_stage2, json_safe


ROOT = Path(__file__).resolve().parents[2]
PREDICTIONS = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_5_vs_v2_v3_consumed/consumed_benchmark_predictions.parquet"
CACHE_DIR = ROOT / "Integrated_dataset/cache/gdT_prediction/gdtai_v4_3_common_lockbox"
MODEL_DIR = ROOT / "Integrated_dataset/models/gdT_prediction/gdtai_v4_5_development"
TABLES = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_5_failure_audit"
LOGS = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_5_failure_audit"
SOURCE = "GSE159251"
AUDIT_GENES = [
    "TRDC", "TRGC1", "TRGC2", "TRDV1", "TRDV2", "TRDV3",
    "TRAV14DV4", "TRAV23DV6", "TRAV29DV5", "TRAV38-2DV8",
    "TRGV1", "TRGV2", "TRGV3", "TRGV4", "TRGV5", "TRGV7", "TRGV8", "TRGV9",
    "TRAC", "TRBC1", "TRBC2", "CD3D", "CD3E", "CD3G", "CD4", "FOXP3",
    "KLRD1", "NKG7", "GNLY", "FCER1G", "TYROBP",
]


def main() -> None:
    contract = json.loads((MODEL_DIR / "model_contract.json").read_text())
    cache_summary = json.loads((CACHE_DIR / "lockbox_feature_cache_summary.json").read_text())
    cache_path = CACHE_DIR / "lockbox_gene_features.npy"
    manifest_path = CACHE_DIR / "lockbox_feature_manifest.csv"
    if sha256_file(cache_path) != cache_summary["feature_cache_sha256"]:
        raise RuntimeError("Benchmark feature cache checksum mismatch")
    frame = pd.read_parquet(PREDICTIONS).reset_index(drop=True)
    matrix = np.load(cache_path, mmap_mode="r")
    names = pd.read_csv(manifest_path).sort_values("feature_index").gene.astype(str).tolist()
    lookup = {gene: index for index, gene in enumerate(names)}
    source_mask = frame.source_gse_id.eq(SOURCE).to_numpy()
    if int(source_mask.sum()) != 32420 or not frame.loc[source_mask, "truth_class"].eq("abT_gold").all():
        raise RuntimeError("GSE159251 benchmark membership changed")

    groups = {
        "v4_5_highest_f1_fp": source_mask & frame.v4_5_highest_f1.to_numpy(bool),
        "v4_5_high_purity_fp": source_mask & frame.v4_5_high_purity.to_numpy(bool),
        "v4_5_highest_f1_tn": source_mask & ~frame.v4_5_highest_f1.to_numpy(bool),
        "v3_balanced_fp": source_mask & frame.v3_balanced.to_numpy(bool),
        "v2_high_purity_fp": source_mask & frame.v2_high_purity.to_numpy(bool),
    }
    group_summary = pd.DataFrame([
        {"group": group, "n_cells": int(mask.sum()), "fraction_of_source": float(mask.mean() / source_mask.mean())}
        for group, mask in groups.items()
    ])

    present_genes = [gene for gene in AUDIT_GENES if gene in lookup]
    expression_rows = []
    for group, mask in groups.items():
        values = np.asarray(matrix[np.flatnonzero(mask)][:, [lookup[gene] for gene in present_genes]])
        for index, gene in enumerate(present_genes):
            expression_rows.append({
                "group": group, "gene": gene, "mean_log1p_cp10k": float(values[:, index].mean()),
                "detected_fraction": float((values[:, index] > 0).mean()),
            })
    expression = pd.DataFrame(expression_rows)

    overlap_rows = []
    v45 = groups["v4_5_highest_f1_fp"]
    for comparator, call_column in [
        ("v3_balanced", "v3_balanced"),
        ("v2_high_f1", "v2_high_f1"),
        ("v2_high_purity", "v2_high_purity"),
    ]:
        other = source_mask & frame[call_column].to_numpy(bool)
        overlap_rows.append({
            "comparator": comparator,
            "both_positive": int((v45 & other).sum()),
            "v4_5_only": int((v45 & ~other).sum()),
            "comparator_only": int((~v45 & other).sum()),
            "both_negative": int((source_mask & ~v45 & ~other).sum()),
        })
    overlap = pd.DataFrame(overlap_rows)

    source_rows = np.flatnonzero(source_mask)
    stage2_columns = [lookup[gene] for gene in contract["stage2_base_feature_names"]]
    stage2, effective_names = compose_stage2(
        np.asarray(matrix[source_rows][:, stage2_columns]), contract["stage2_base_feature_names"]
    )
    if effective_names != contract["stage2_effective_feature_names"]:
        raise RuntimeError("V4.5 effective feature order changed")
    booster = xgb.Booster()
    booster.load_model(MODEL_DIR / "stage2_gdt_classifier.ubj")
    contributions = booster.predict(xgb.DMatrix(stage2), pred_contribs=True)[:, :-1]
    local = frame.iloc[source_rows].reset_index(drop=True)
    contribution_groups = {
        "v4_5_highest_f1_fp": local.v4_5_highest_f1.to_numpy(bool),
        "v4_5_high_purity_fp": local.v4_5_high_purity.to_numpy(bool),
        "v4_5_highest_f1_tn": ~local.v4_5_highest_f1.to_numpy(bool),
    }
    contribution_rows = []
    for group, mask in contribution_groups.items():
        mean_values = contributions[mask].mean(axis=0)
        for feature, value in zip(effective_names, mean_values):
            contribution_rows.append({"group": group, "feature": feature, "mean_contribution": float(value)})
    contribution = pd.DataFrame(contribution_rows)

    gain_rows = []
    for key, value in booster.get_score(importance_type="gain").items():
        index = int(key[1:])
        gain_rows.append({"feature": effective_names[index], "feature_index": index, "gain": float(value)})
    gain = pd.DataFrame(gain_rows).sort_values("gain", ascending=False)

    TABLES.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    group_summary.to_csv(TABLES / "group_summary.csv", index=False)
    expression.to_csv(TABLES / "marker_expression_by_group.csv", index=False)
    overlap.to_csv(TABLES / "prediction_overlap.csv", index=False)
    contribution.to_csv(TABLES / "mean_stage2_contributions.csv", index=False)
    gain.to_csv(TABLES / "stage2_gain_importance.csv", index=False)

    hf1 = expression[expression.group.eq("v4_5_highest_f1_fp")].set_index("gene")
    hp = expression[expression.group.eq("v4_5_high_purity_fp")].set_index("gene")
    hf1_rows = np.flatnonzero(groups["v4_5_highest_f1_fp"])
    trdv_columns = [lookup[gene] for gene in ["TRDV1", "TRDV2", "TRDV3"] if gene in lookup]
    trgv_columns = [
        lookup[gene] for gene in ["TRGV1", "TRGV2", "TRGV3", "TRGV4", "TRGV5", "TRGV7", "TRGV8", "TRGV9"]
        if gene in lookup
    ]
    hf1_trdv_any = np.asarray(matrix[hf1_rows][:, trdv_columns]).max(axis=1) > 0
    hf1_trgv_any = np.asarray(matrix[hf1_rows][:, trgv_columns]).max(axis=1) > 0
    summary = {
        "status": "PASS_V4_5_FAILURE_AUDIT",
        "source": SOURCE,
        "n_abt_gold": int(source_mask.sum()),
        "n_v4_5_highest_f1_fp": int(groups["v4_5_highest_f1_fp"].sum()),
        "n_v4_5_high_purity_fp": int(groups["v4_5_high_purity_fp"].sum()),
        "v4_5_highest_f1_fp_fraction": float(groups["v4_5_highest_f1_fp"].sum() / source_mask.sum()),
        "v4_5_high_purity_fp_fraction": float(groups["v4_5_high_purity_fp"].sum() / source_mask.sum()),
        "highest_f1_fp_trdv_any_detected_fraction": float(hf1_trdv_any.mean()),
        "highest_f1_fp_trgv_any_detected_fraction": float(hf1_trgv_any.mean()),
        "high_purity_fp_trdv1_detected_fraction": float(hp.loc["TRDV1", "detected_fraction"]),
        "interpretation": (
            "V4.5-only false positives are dominated by alpha-beta cells with broad TRGV and alpha/beta V evidence, "
            "while most lack dedicated TRDV expression; standalone gamma and symmetric receptor aggregates require ablation."
        ),
        "input_hashes": {
            "predictions": sha256_file(PREDICTIONS),
            "feature_cache": cache_summary["feature_cache_sha256"],
            "model_contract": sha256_file(MODEL_DIR / "model_contract.json"),
        },
    }
    (LOGS / "summary.json").write_text(json.dumps(json_safe(summary), indent=2, sort_keys=True) + "\n")
    print(json.dumps(json_safe(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
