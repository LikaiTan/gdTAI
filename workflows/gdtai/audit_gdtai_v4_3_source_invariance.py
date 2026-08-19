#!/usr/bin/env python3
"""Audit source-invariant receptor contrasts on frozen gdTAI development labels."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from workflows.gdtai.train_gdtai_v4_2_nested import (  # noqa: E402
    choose_profile_threshold,
    exclusion_flags,
    profile_metrics,
)


TABLE_DIR = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_3_source_invariance"
LOG_DIR = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_3_source_invariance"
LABEL_PATH = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_ground_truth/v4_2_label_manifest.parquet"
FEATURE_PATH = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_training/feature_manifest.csv"
CACHE_PATH = ROOT / "Integrated_dataset/cache/gdT_prediction/gdtai_v4_2_training/gene_features.npy"
HELDOUTS = ["HRA005041", "GSE144469", "MalteGDT"]

PROFILES = {
    "balanced": {
        "objective_beta": 1.0,
        "minimum_macro_recall": 0.80,
        "minimum_per_source_recall": 0.70,
        "maximum_abt_fpr": 0.002,
        "maximum_tier1_nk_fpr": 0.01,
        "maximum_tier2_macro_nk_fpr": 0.02,
    },
    "high_purity": {
        "objective_beta": 0.5,
        "minimum_macro_recall": 0.70,
        "minimum_per_source_recall": 0.60,
        "maximum_abt_fpr": 0.001,
        "maximum_tier1_nk_fpr": 0.005,
        "maximum_tier2_macro_nk_fpr": 0.01,
    },
}


def family_max(x: np.ndarray, names: list[str], prefixes: tuple[str, ...]) -> np.ndarray:
    output = np.zeros(x.shape[0], dtype=np.float32)
    columns = [index for index, gene in enumerate(names) if gene.startswith(prefixes)]
    for position, column in enumerate(columns, 1):
        np.maximum(output, np.asarray(x[:, column]), out=output)
        if position % 25 == 0:
            print(f"  family max {prefixes}: {position}/{len(columns)}", flush=True)
    return output


def family_sum(x: np.ndarray, names: list[str], prefixes: tuple[str, ...]) -> np.ndarray:
    output = np.zeros(x.shape[0], dtype=np.float32)
    columns = [index for index, gene in enumerate(names) if gene.startswith(prefixes)]
    for column in columns:
        output += np.asarray(x[:, column])
    return output


def chunked_exclusions(x: np.ndarray, names: list[str], chunk: int = 250_000) -> np.ndarray:
    output = np.zeros(x.shape[0], dtype=bool)
    for start in range(0, x.shape[0], chunk):
        stop = min(x.shape[0], start + chunk)
        output[start:stop] = exclusion_flags(np.asarray(x[start:stop]), names)[2]
    return output


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    labels = pd.read_parquet(LABEL_PATH).reset_index(drop=True)
    features = pd.read_csv(FEATURE_PATH).sort_values("feature_index")
    names = features.gene.astype(str).tolist()
    lookup = {gene: index for index, gene in enumerate(names)}
    x = np.load(CACHE_PATH, mmap_mode="r")
    if x.shape != (len(labels), len(names)):
        raise RuntimeError("Feature cache shape does not match frozen labels/features")

    print("Calculating receptor-family summaries", flush=True)
    delta_v = family_max(x, names, ("TRDV",))
    gamma_v = family_max(x, names, ("TRGV",))
    alpha_v = family_max(x, names, ("TRAV",))
    beta_v = family_max(x, names, ("TRBV",))
    ab_v = np.maximum(alpha_v, beta_v)
    delta_v_sum = family_sum(x, names, ("TRDV",))
    ab_v_sum = family_sum(x, names, ("TRAV", "TRBV"))
    trdc = np.asarray(x[:, lookup["TRDC"]])
    trgc = np.maximum(np.asarray(x[:, lookup["TRGC1"]]), np.asarray(x[:, lookup["TRGC2"]]))
    ab_constant = np.maximum.reduce([
        np.asarray(x[:, lookup["TRAC"]]),
        np.asarray(x[:, lookup["TRBC1"]]),
        np.asarray(x[:, lookup["TRBC2"]]),
    ])
    scores = {
        "delta_v_max": delta_v,
        "delta_v_minus_ab_v": delta_v - ab_v,
        "delta_v_minus_ab_constant": delta_v - ab_constant,
        "delta_v_plus_trdc_minus_ab_v": delta_v + 0.25 * trdc - 0.25 * ab_v,
        "delta_v_plus_trgc_minus_ab_v": delta_v + 0.25 * trgc - 0.25 * ab_v,
        "delta_v_plus_constants_minus_ab_v": delta_v + 0.20 * trdc + 0.10 * trgc - 0.25 * ab_v,
        "delta_v_fraction_of_v_max": delta_v / (delta_v + ab_v + 1e-4),
        "delta_v_sum_minus_scaled_ab_v_sum": delta_v_sum - 0.10 * ab_v_sum,
        "delta_v_detected_minus_ab_v_detected": (delta_v > 0).astype(np.float32) - 0.25 * (ab_v > 0),
    }
    score_manifest = pd.DataFrame({
        "score_id": list(scores),
        "definition": [
            "max(TRDV1,TRDV2,TRDV3)",
            "max(TRDV)-max(TRAV,TRBV)",
            "max(TRDV)-max(TRAC,TRBC1,TRBC2)",
            "max(TRDV)+0.25*TRDC-0.25*max(TRAV,TRBV)",
            "max(TRDV)+0.25*max(TRGC1,TRGC2)-0.25*max(TRAV,TRBV)",
            "max(TRDV)+0.20*TRDC+0.10*max(TRGC1,TRGC2)-0.25*max(TRAV,TRBV)",
            "max(TRDV)/(max(TRDV)+max(TRAV,TRBV)+1e-4)",
            "sum(TRDV)-0.10*sum(TRAV,TRBV)",
            "I(any TRDV)-0.25*I(any TRAV/TRBV)",
        ],
    })
    score_manifest.to_csv(TABLE_DIR / "score_manifest.csv", index=False)

    print("Calculating frozen exclusion flags", flush=True)
    excluded = chunked_exclusions(x, names)
    development = labels.allow_fit.to_numpy()
    rows = []
    for heldout in HELDOUTS:
        print(f"Screening held-out source: {heldout}", flush=True)
        train = development & ~labels.source_gse_id.eq(heldout).to_numpy()
        test = development & labels.source_gse_id.eq(heldout).to_numpy()
        train_frame = labels.loc[train].reset_index(drop=True)
        test_frame = labels.loc[test].reset_index(drop=True)
        for score_id, score in scores.items():
            for profile, spec in PROFILES.items():
                threshold, inner = choose_profile_threshold(
                    train_frame, np.ones(train.sum(), dtype=np.float32), 0.0,
                    score[train], excluded[train], spec,
                )
                record = {
                    "heldout_source": heldout,
                    "score_id": score_id,
                    "profile": profile,
                    "inner_eligible": inner.get("eligible") is True,
                    "threshold": threshold,
                }
                if inner.get("eligible") is True:
                    calls = (score[test] >= threshold) & ~excluded[test]
                    outer = profile_metrics(test_frame, calls, score[test])
                    record.update({
                        "inner_f1": inner["f1"],
                        "inner_macro_recall": inner["macro_recall"],
                        "inner_abt_fpr": inner["abt_fpr"],
                        "inner_max_tier1_nk_fpr": max(inner["tier1_nk_fpr_by_source"].values(), default=0.0),
                        "inner_tier2_macro_nk_fpr": inner["tier2_macro_nk_fpr"],
                        "outer_precision": outer["precision"],
                        "outer_recall": outer["recall"],
                        "outer_abt_fpr": outer["abt_fpr"],
                        "outer_f1": outer["f1"],
                    })
                rows.append(record)
    results = pd.DataFrame(rows)
    results.to_csv(TABLE_DIR / "source_heldout_score_screen.csv", index=False)
    summary = (
        results[results.inner_eligible]
        .groupby(["score_id", "profile"], as_index=False)
        .agg(
            eligible_folds=("heldout_source", "nunique"),
            minimum_outer_recall=("outer_recall", "min"),
            mean_outer_recall=("outer_recall", "mean"),
            maximum_outer_abt_fpr=("outer_abt_fpr", "max"),
            mean_outer_f1=("outer_f1", "mean"),
        )
        .sort_values(["eligible_folds", "minimum_outer_recall"], ascending=False)
    )
    summary.to_csv(TABLE_DIR / "source_heldout_score_summary.csv", index=False)
    payload = {
        "status": "PASS_SOURCE_INVARIANCE_AUDIT",
        "n_scores": len(scores),
        "n_outer_folds": len(HELDOUTS),
        "lockbox_scored": False,
        "model_promoted": False,
    }
    (LOG_DIR / "summary.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(summary.to_string(index=False), flush=True)
    print(json.dumps(payload, indent=2), flush=True)


if __name__ == "__main__":
    main()
