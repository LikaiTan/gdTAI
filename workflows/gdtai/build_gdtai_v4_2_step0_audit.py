#!/usr/bin/env python3
"""Build the read-only gdTAI V4.2 label audit and OOF counterfactual report.

This script never fits a model. It reuses checksum-bound V4.1 out-of-fold
probabilities to isolate the effect of NK-control provenance on the failed
Stage-1 gate.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import math
import subprocess
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "configs/models/gdtai/v4_2_precommit.json"


def resolve(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def sha256_file(path: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def array_sha256(values: np.ndarray) -> str:
    return hashlib.sha256(memoryview(np.ascontiguousarray(values))).hexdigest()


def threshold_candidates(values: np.ndarray) -> np.ndarray:
    score = np.asarray(values, dtype=np.float64)
    score = score[np.isfinite(score)]
    return np.unique(np.clip(np.concatenate([score, np.asarray([0.0, 1.0])]), 0.0, 1.0))


def count_greater_equal(values: np.ndarray, thresholds: np.ndarray) -> np.ndarray:
    ordered = np.sort(np.asarray(values, dtype=np.float64))
    return ordered.size - np.searchsorted(ordered, thresholds, side="left")


def wilson_interval(successes: int, total: int, z: float) -> tuple[float, float]:
    if total == 0:
        return math.nan, math.nan
    p = successes / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt((p * (1 - p) + z * z / (4 * total)) / total) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def markdown_table(frame: pd.DataFrame, columns: list[str], formats: dict[str, str] | None = None) -> str:
    formats = formats or {}
    selected = frame.loc[:, columns].copy()
    for column, template in formats.items():
        selected[column] = selected[column].map(
            lambda value: "NA" if pd.isna(value) else template.format(value)
        )
    return selected.to_markdown(index=False)


def read_inputs(config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, list[dict[str, Any]]]:
    inputs = config["inputs"]
    expected = config["expected_sha256"]
    checks: list[dict[str, Any]] = []
    paths = {
        "cell_manifest": resolve(inputs["cell_manifest"]),
        "feature_manifest": resolve(inputs["feature_manifest"]),
        "gene_feature_cache": resolve(inputs["gene_feature_cache"]),
    }
    for key, path in paths.items():
        observed = sha256_file(path)
        wanted = expected[key]
        checks.append(
            {
                "check": f"sha256_{key}",
                "status": "PASS" if observed == wanted else "FAIL",
                "detail": observed,
            }
        )
        if observed != wanted:
            raise RuntimeError(f"Checksum mismatch for {path}: {observed} != {wanted}")
    cells = pd.read_csv(paths["cell_manifest"], low_memory=False)
    features = pd.read_csv(paths["feature_manifest"])
    matrix = np.load(paths["gene_feature_cache"], mmap_mode="r")
    if matrix.shape != (cells.shape[0], features.shape[0]):
        raise RuntimeError(f"Feature cache shape {matrix.shape} does not match manifests")
    checks.append({"check": "feature_cache_shape", "status": "PASS", "detail": str(matrix.shape)})
    return cells, features, matrix, checks


def expression_audit(
    cells: pd.DataFrame,
    features: pd.DataFrame,
    matrix: np.ndarray,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    policy = config["expression_audit"]
    genes = policy["cd3_genes"] + policy["abt_constant_genes"] + policy["gdt_constant_genes"]
    feature_index = features.set_index("gene")["feature_index"].astype(int).to_dict()
    missing = [gene for gene in genes if gene not in feature_index]
    if missing:
        raise RuntimeError(f"Expression-audit genes are absent: {missing}")

    nk_rows = cells.index[cells["stage1_role"].eq("nk_negative")].to_numpy(dtype=np.int64)
    values = np.asarray(matrix[nk_rows][:, [feature_index[gene] for gene in genes]], dtype=np.float32)
    detected = values > 0
    offsets = {gene: idx for idx, gene in enumerate(genes)}
    cd3 = detected[:, [offsets[gene] for gene in policy["cd3_genes"]]].sum(axis=1)
    ab = detected[:, [offsets[gene] for gene in policy["abt_constant_genes"]]].any(axis=1)
    gd = detected[:, [offsets[gene] for gene in policy["gdt_constant_genes"]]].any(axis=1)
    cd3_coherent = cd3 >= int(policy["minimum_cd3_detected"])
    coherent_t = cd3_coherent & (ab | gd)

    flags = pd.DataFrame(
        {
            "manifest_row": nk_rows,
            "cd3_two_of_three": cd3_coherent,
            "ab_constant_detected": ab,
            "gd_constant_detected": gd,
            "expression_coherent_t": coherent_t,
        }
    ).set_index("manifest_row")
    for gene in genes:
        flags[f"detected_{gene}"] = detected[:, offsets[gene]]

    nk = cells.loc[nk_rows, ["source_gse_id", "nk_annotation_strength"]].copy()
    nk = nk.join(flags)
    source_rows: list[dict[str, Any]] = []
    strength_rows: list[dict[str, Any]] = []
    for source, frame in nk.groupby("source_gse_id", sort=True):
        source_rows.append(summarize_expression_group(source, "all", frame, genes))
        for strength, subset in frame.groupby("nk_annotation_strength", sort=True):
            source_rows.append(summarize_expression_group(source, f"strength_{strength:g}", subset, genes))
    for strength, frame in nk.groupby("nk_annotation_strength", sort=True):
        strength_rows.append(summarize_expression_group("ALL", f"strength_{strength:g}", frame, genes))
    strength_rows.append(summarize_expression_group("ALL", "all", nk, genes))
    return pd.DataFrame(source_rows), pd.DataFrame(strength_rows), flags["expression_coherent_t"].reindex(cells.index, fill_value=False).to_numpy(bool)


def summarize_expression_group(source: str, subset: str, frame: pd.DataFrame, genes: list[str]) -> dict[str, Any]:
    row: dict[str, Any] = {
        "source_gse_id": source,
        "subset": subset,
        "n_cells": int(frame.shape[0]),
        "cd3_two_of_three_fraction": float(frame["cd3_two_of_three"].mean()),
        "ab_constant_fraction": float(frame["ab_constant_detected"].mean()),
        "gd_constant_fraction": float(frame["gd_constant_detected"].mean()),
        "expression_coherent_t_fraction": float(frame["expression_coherent_t"].mean()),
        "ab_cd3_coherent_fraction": float((frame["ab_constant_detected"] & frame["cd3_two_of_three"]).mean()),
        "gd_cd3_coherent_fraction": float((frame["gd_constant_detected"] & frame["cd3_two_of_three"]).mean()),
    }
    for gene in genes:
        row[f"detected_{gene}_fraction"] = float(frame[f"detected_{gene}"].mean())
    return row


def annotation_provenance(cells: pd.DataFrame) -> pd.DataFrame:
    nk = cells[cells["stage1_role"].eq("nk_negative")]
    counts = (
        nk.groupby(["source_gse_id", "nk_annotation_strength"], observed=True)
        .size()
        .unstack(fill_value=0)
        .rename(columns={0.5: "single_annotation_n", 1.0: "dual_annotation_n"})
        .reset_index()
    )
    for column in ["single_annotation_n", "dual_annotation_n"]:
        if column not in counts:
            counts[column] = 0
    counts["n_cells"] = counts["single_annotation_n"] + counts["dual_annotation_n"]
    counts["dual_annotation_fraction"] = counts["dual_annotation_n"] / counts["n_cells"]
    return counts.sort_values(["n_cells", "source_gse_id"], ascending=[False, True]).reset_index(drop=True)


def highest_recall_feasible_threshold(
    probability: np.ndarray,
    local: pd.DataFrame,
    primary_sources: list[str],
    heldout: str,
    policy: dict[str, Any],
) -> tuple[float, float, float]:
    thresholds = threshold_candidates(probability)
    source = local["source_gse_id"].to_numpy(dtype=object)
    gdt = local["truth_class"].eq("gdT_primary").to_numpy(bool)
    abt = (local["stage1_role"].eq("t_positive") & local["has_any_ab_tcr"].astype(bool)).to_numpy(bool)
    valid = np.ones(thresholds.size, dtype=bool)
    gdt_recalls: list[np.ndarray] = []
    abt_recalls: list[np.ndarray] = []
    for value in primary_sources:
        if value == heldout:
            continue
        gd_values = probability[gdt & (source == value)]
        ab_values = probability[abt & (source == value)]
        if gd_values.size == 0 or ab_values.size == 0:
            raise RuntimeError(f"Missing primary Stage-1 stratum for {value}")
        gd_recall = count_greater_equal(gd_values, thresholds) / gd_values.size
        ab_recall = count_greater_equal(ab_values, thresholds) / ab_values.size
        gdt_recalls.append(gd_recall)
        abt_recalls.append(ab_recall)
        valid &= gd_recall >= float(policy["gdt_recall_per_source"])
        valid &= ab_recall >= float(policy["abt_recall_per_source"])
    eligible = np.flatnonzero(valid)
    if eligible.size == 0:
        raise RuntimeError(f"No recall-feasible Stage-1 threshold for {heldout}")
    index = int(eligible[-1])
    return float(thresholds[index]), float(np.min(np.vstack(gdt_recalls)[:, index])), float(np.min(np.vstack(abt_recalls)[:, index]))


def counterfactual_audit(
    cells: pd.DataFrame,
    coherent_t: np.ndarray,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, Any]]]:
    checkpoint_root = resolve(config["inputs"]["v4_1_checkpoints"])
    policy = config["stage1_counterfactual"]
    primary_sources = list(config["primary_t_sources"])
    minimum_n = int(config["nk_truth"]["minimum_source_size_for_gate"])
    z = float(policy["wilson_confidence_z"])
    records: list[dict[str, Any]] = []
    source_records: list[dict[str, Any]] = []
    checks: list[dict[str, Any]] = []
    checkpoint_paths = sorted(checkpoint_root.glob("outer_*/stage1/*/checkpoint.json"))
    if not checkpoint_paths:
        raise RuntimeError(f"No Stage-1 checkpoints under {checkpoint_root}")

    for checkpoint_path in checkpoint_paths:
        metadata = json.loads(checkpoint_path.read_text())
        outer_fold_id = metadata["contract"]["outer_fold_id"]
        heldout = outer_fold_id.split("_", 2)[2]
        eligible = cells["stage1_role"].isin(["t_positive", "nk_negative"]) & cells["source_gse_id"].ne(heldout)
        train_rows = cells.index[eligible].to_numpy(dtype=np.int64)
        observed_hash = array_sha256(train_rows)
        expected_hash = metadata["contract"]["train_rows_sha256"]
        checks.append(
            {
                "check": f"checkpoint_rows::{outer_fold_id}::{metadata['candidate_id']}",
                "status": "PASS" if observed_hash == expected_hash else "FAIL",
                "detail": observed_hash,
            }
        )
        if observed_hash != expected_hash:
            raise RuntimeError(f"Train-row hash mismatch for {checkpoint_path}")
        with np.load(checkpoint_path.with_name("probabilities.npz")) as archive:
            probability = np.asarray(archive["calibrated_oof"], dtype=np.float64)
        if probability.size != train_rows.size or not np.isfinite(probability).all():
            raise RuntimeError(f"Invalid OOF probability array for {checkpoint_path}")
        local = cells.loc[train_rows].copy()
        threshold, minimum_gdt, minimum_abt = highest_recall_feasible_threshold(
            probability, local, primary_sources, heldout, policy
        )
        source = local["source_gse_id"].to_numpy(dtype=object)
        nk = local["stage1_role"].eq("nk_negative").to_numpy(bool)
        strength = local["nk_annotation_strength"].to_numpy(dtype=float)
        local_coherent = coherent_t[train_rows]
        variants = {
            "all_v4_1_controls": nk,
            "dual_annotation_primary": nk & (strength == float(config["nk_truth"]["primary_strength"])),
            "single_annotation_stress": nk & (strength == float(config["nk_truth"]["stress_strength"])),
            "all_minus_Tcoherent_diagnostic": nk & (~local_coherent),
            "single_minus_Tcoherent_diagnostic": nk & (strength == float(config["nk_truth"]["stress_strength"])) & (~local_coherent),
        }
        for variant, mask in variants.items():
            subset_rows: list[dict[str, Any]] = []
            for value in sorted(pd.unique(source[mask]).tolist()):
                values = probability[mask & (source == value)]
                passed = int((values >= threshold).sum())
                lower, upper = wilson_interval(passed, int(values.size), z)
                item = {
                    "outer_fold_id": outer_fold_id,
                    "heldout_source": heldout,
                    "candidate_id": metadata["candidate_id"],
                    "family": metadata["family"],
                    "variant": variant,
                    "source_gse_id": value,
                    "n_cells": int(values.size),
                    "n_passed": passed,
                    "passage": passed / values.size,
                    "passage_wilson_low": lower,
                    "passage_wilson_high": upper,
                    "eligible_source_n_ge_minimum": bool(values.size >= minimum_n),
                }
                source_records.append(item)
                subset_rows.append(item)
            if not subset_rows:
                continue
            subset = pd.DataFrame(subset_rows)
            large = subset[subset["eligible_source_n_ge_minimum"]]
            worst_all = subset.sort_values(["passage", "n_cells", "source_gse_id"], ascending=[False, False, True]).iloc[0]
            worst_large = None if large.empty else large.sort_values(["passage", "n_cells", "source_gse_id"], ascending=[False, False, True]).iloc[0]
            total = int(subset["n_cells"].sum())
            passed = int(subset["n_passed"].sum())
            records.append(
                {
                    "outer_fold_id": outer_fold_id,
                    "heldout_source": heldout,
                    "candidate_id": metadata["candidate_id"],
                    "family": metadata["family"],
                    "parameters": json.dumps(metadata["parameters"], sort_keys=True),
                    "threshold": threshold,
                    "minimum_source_gdt_recall": minimum_gdt,
                    "minimum_source_abt_recall": minimum_abt,
                    "variant": variant,
                    "n_cells": total,
                    "pooled_passage": passed / total,
                    "macro_source_passage": float(subset["passage"].mean()),
                    "maximum_source_passage": float(worst_all["passage"]),
                    "worst_source": str(worst_all["source_gse_id"]),
                    "worst_source_n": int(worst_all["n_cells"]),
                    "maximum_source_passage_n_ge_minimum": math.nan if worst_large is None else float(worst_large["passage"]),
                    "maximum_source_wilson_high_n_ge_minimum": math.nan if worst_large is None else float(large["passage_wilson_high"].max()),
                    "worst_source_n_ge_minimum": "" if worst_large is None else str(worst_large["source_gse_id"]),
                    "worst_source_n_ge_minimum_n": 0 if worst_large is None else int(worst_large["n_cells"]),
                    "legacy_all_source_cap_pass": bool(worst_all["passage"] <= float(policy["maximum_nk_passage_per_source"])),
                    "sample_size_aware_cap_pass": bool(
                        worst_large is not None
                        and worst_large["passage"] <= float(policy["maximum_nk_passage_per_source"])
                    ),
                    "wilson_upper_cap_pass": bool(
                        not large.empty
                        and float(large["passage_wilson_high"].max()) <= float(policy["maximum_nk_passage_per_source"])
                    ),
                }
            )
    return pd.DataFrame(records), pd.DataFrame(source_records), checks


def build_figures(provenance: pd.DataFrame, expression: pd.DataFrame, counterfactual: pd.DataFrame, figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 11, "axes.labelsize": 9})

    view = provenance.sort_values("n_cells", ascending=True)
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.barh(view["source_gse_id"], view["single_annotation_n"], color="#8b9eb7", label="Single annotation")
    ax.barh(
        view["source_gse_id"], view["dual_annotation_n"], left=view["single_annotation_n"],
        color="#14866d", label="Dual annotation agreement",
    )
    ax.set_xscale("symlog", linthresh=100)
    ax.set_xlabel("NK-control cells (symlog scale)")
    ax.set_title("V4.1 NK controls were dominated by single-annotation labels")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(figure_dir / "nk_control_provenance_by_source.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    view = expression[(expression["subset"] == "all") & (expression["n_cells"] >= 100)].copy()
    view = view.sort_values("expression_coherent_t_fraction", ascending=True)
    colors = np.where(view["expression_coherent_t_fraction"] >= 0.5, "#b33a3a", "#3a78a8")
    fig, ax = plt.subplots(figsize=(9, 6.8))
    ax.barh(view["source_gse_id"], 100 * view["expression_coherent_t_fraction"], color=colors)
    ax.axvline(50, color="#555555", lw=1, ls="--")
    ax.set_xlim(0, 100)
    ax.set_xlabel("Cells with CD3 2-of-3 and any TCR constant gene (%)")
    ax.set_title("Several single-annotation NK sources are transcriptomically T-like")
    fig.tight_layout()
    fig.savefig(figure_dir / "nk_control_tlineage_coherence.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    variants = ["all_v4_1_controls", "dual_annotation_primary", "single_annotation_stress"]
    view = counterfactual[counterfactual["variant"].isin(variants)].copy()
    labels = {
        "all_v4_1_controls": "All V4.1 controls",
        "dual_annotation_primary": "Dual-annotation primary",
        "single_annotation_stress": "Single-annotation stress",
    }
    fig, ax = plt.subplots(figsize=(9, 5.2))
    palette = {"all_v4_1_controls": "#9a3e3e", "dual_annotation_primary": "#14866d", "single_annotation_stress": "#c58b2a"}
    xmap = {variant: index for index, variant in enumerate(variants)}
    offsets = {fold: offset for fold, offset in zip(sorted(view["outer_fold_id"].unique()), [-0.12, 0.0, 0.12])}
    for _, row in view.iterrows():
        ax.scatter(
            xmap[row["variant"]] + offsets[row["outer_fold_id"]],
            100 * row["maximum_source_passage_n_ge_minimum"],
            color=palette[row["variant"]], alpha=0.72, s=28, edgecolor="white", linewidth=0.35,
        )
    medians = view.groupby("variant")["maximum_source_passage_n_ge_minimum"].median()
    for variant, median in medians.items():
        ax.hlines(100 * median, xmap[variant] - 0.25, xmap[variant] + 0.25, color="#111111", lw=2)
    ax.axhline(50, color="#555555", lw=1, ls="--", label="Frozen 50% cap")
    ax.set_xticks(range(len(variants)), [labels[value] for value in variants])
    ax.set_ylabel("Worst source passage, n >= 100 (%)")
    ax.set_title("Saved OOF predictions: NK-truth provenance changes Gate C")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figure_dir / "stage1_counterfactual_feasibility.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def render_report(
    cells: pd.DataFrame,
    provenance: pd.DataFrame,
    expression: pd.DataFrame,
    counterfactual: pd.DataFrame,
    checks: pd.DataFrame,
    config: dict[str, Any],
    log_dir: Path,
    static_dir: Path,
) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    static_dir.mkdir(parents=True, exist_ok=True)
    nk = cells[cells["stage1_role"].eq("nk_negative")]
    dual_n = int((nk["nk_annotation_strength"] == 1.0).sum())
    single_n = int((nk["nk_annotation_strength"] == 0.5).sum())
    dual = counterfactual[counterfactual["variant"].eq("dual_annotation_primary")]
    original = counterfactual[counterfactual["variant"].eq("all_v4_1_controls")]
    worst_expression = expression[(expression["subset"] == "all") & (expression["n_cells"] >= 100)].nlargest(6, "expression_coherent_t_fraction")
    compact = dual.sort_values(["heldout_source", "pooled_passage"]).groupby("heldout_source", as_index=False).first()

    report = f"""# gdTAI V4.2 Step 0: why V4.1 failed and how to repair it

## Decision summary

V4.1 did **not** demonstrate that a two-step gdTAI model is intrinsically worse than V2 or V3. It stopped before Stage 2 because the Stage-1 eligibility gate treated a heterogeneous, mostly weak NK-control pool as equally definitive ground truth.

- V4.1 used **{len(nk):,}** strict-NK controls.
- Only **{dual_n:,} ({dual_n / len(nk):.1%})** had dual scVI/author NK agreement.
- **{single_n:,} ({single_n / len(nk):.1%})** had only the scVI NK label after conflict filtering.
- Under the unchanged T-recall thresholds, the worst source with at least 100 cells passed at **{100 * original['maximum_source_passage_n_ge_minimum'].min():.1f}% to {100 * original['maximum_source_passage_n_ge_minimum'].max():.1f}%** for all V4.1 controls.
- The same saved out-of-fold predictions passed only **{100 * dual['maximum_source_passage_n_ge_minimum'].min():.1f}% to {100 * dual['maximum_source_passage_n_ge_minimum'].max():.1f}%** of dual-annotation NK controls. All 15 candidate/fold combinations satisfy the frozen 50% cap under that stronger truth definition.

This is a **label-only counterfactual**, not V4.2 performance. The V4.1 models were still trained using the weak controls, and no model was fitted in this audit.

## Why V4.1 looked weak

1. **Stage 2 never ran.** There is no V4.1 final F1, precision, or recall to compare directly with V2/V3.
2. **The Stage-1 failure was source-specific.** Large single-annotation sources such as GSE144469 and GSE287301 were highly T-like by CD3/TCR-constant expression, while genuine-looking BALF NK controls also showed substantial passage in one outer fold. This combines label uncertainty with real cross-dataset calibration shift.
3. **The gate was brittle.** It used the literal maximum source passage, including a one-cell source, and did not distinguish primary labels from stress labels.
4. **Historical V2/V3 figures are optimistic for model selection.** Their validation cohorts were reused during development. V4.2 must compare algorithms on identical grouped, fold-local splits and report the frozen historical artifacts separately.

![NK-control provenance](../../Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_step0/nk_control_provenance_by_source.png)

## NK-label audit

The expression audit is used only to diagnose label conflict. It does not create, remove, or relabel training truth.

{markdown_table(worst_expression, ['source_gse_id', 'n_cells', 'cd3_two_of_three_fraction', 'expression_coherent_t_fraction', 'ab_cd3_coherent_fraction', 'gd_cd3_coherent_fraction'], {'cd3_two_of_three_fraction': '{:.1%}', 'expression_coherent_t_fraction': '{:.1%}', 'ab_cd3_coherent_fraction': '{:.1%}', 'gd_cd3_coherent_fraction': '{:.1%}'})}

![T-lineage coherence](../../Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_step0/nk_control_tlineage_coherence.png)

## Saved-OOF counterfactual

For each of the 15 completed V4.1 Stage-1 candidate/fold combinations, the audit reconstructed the exact checksum-bound training rows, loaded calibrated OOF probabilities, and selected the highest threshold satisfying per-source gdT recall >= 0.99 and productive alpha-beta T-cell recall >= 0.98. Only the NK evaluation mask changed.

{markdown_table(compact, ['heldout_source', 'family', 'threshold', 'minimum_source_gdt_recall', 'minimum_source_abt_recall', 'pooled_passage', 'maximum_source_passage_n_ge_minimum', 'worst_source_n_ge_minimum'], {'threshold': '{:.5f}', 'minimum_source_gdt_recall': '{:.2%}', 'minimum_source_abt_recall': '{:.2%}', 'pooled_passage': '{:.2%}', 'maximum_source_passage_n_ge_minimum': '{:.2%}'})}

![Counterfactual feasibility](../../Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_step0/stage1_counterfactual_feasibility.png)

## Frozen V4.2 repair

The detailed precommit is in `docs/GDTAI_V4_2_PRECOMMITTED_PLAN.md`. Its key changes are:

- Primary Stage-1 NK truth requires dual scVI/author NK agreement, no productive TCR, and no doublet flag.
- Single-annotation NK cells receive zero fitting and selection weight; they remain a required stress test.
- CD3/TCR expression coherence remains diagnostic only and cannot define truth.
- Stage 1 remains a high-recall T-cell gate. It cannot use a hard cytotoxic-gene veto, and Stage 2 has no hard TRDC/TRDV requirement.
- Stage 2 uses individual genes and the cross-fitted Stage-1 probability, never TRD/TRAB scores, source IDs, or batch IDs; total features remain below 300.
- Source-size-aware constraints use strata with at least 100 controls. Smaller strata and Wilson intervals are reported separately rather than controlling a literal maximum.
- V2-like and V3-like algorithms are retrained fold-locally on the same labels and splits. Frozen V2/V3 artifacts remain descriptive because their development data overlap the benchmark.
- Promotion requires a clinically meaningful F1 improvement under grouped resampling while retaining the frozen paired-alpha-beta and NK false-positive guardrails.

## QC and supervision gate

{markdown_table(checks, ['check', 'status', 'detail'])}

**Result: {('PASS' if checks['status'].eq('PASS').all() else 'FAIL')}.** This authorizes review of the V4.2 protocol only. It does not authorize fitting. The next action requires explicit supervision approval.

## Reproducible outputs

- Full tables: `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_step0/`
- Machine summary: `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_step0/summary.json`
- Frozen protocol: `docs/GDTAI_V4_2_PRECOMMITTED_PLAN.md`
- Reader HTML: `gdT_prediction/gdtai_v4_2_step0/index.html`
- PDF: `gdT_prediction/gdtai_v4_2_step0/gdtai_v4_2_step0_report.pdf`
"""
    md_path = log_dir / "gdtai_v4_2_step0_summary.md"
    md_path.write_text(report, encoding="utf-8")

    css = """
body{font-family:Arial,Helvetica,sans-serif;color:#17212b;line-height:1.48;max-width:1120px;margin:0 auto;padding:30px 42px;background:#fff}
h1{font-size:30px;color:#17324d;border-bottom:4px solid #14866d;padding-bottom:12px}h2{font-size:21px;color:#17324d;margin-top:32px}
p,li{font-size:14px}table{border-collapse:collapse;width:100%;font-size:11px;margin:16px 0 24px}th{background:#17324d;color:#fff;text-align:left}th,td{padding:6px 7px;border:1px solid #cdd6df;vertical-align:top}tr:nth-child(even){background:#f3f6f8}
img{display:block;max-width:96%;max-height:690px;margin:22px auto}code{background:#edf1f4;padding:1px 4px;border-radius:2px}strong{color:#102f43}
@media print{body{max-width:none;padding:8mm 10mm}h1{font-size:24px}h2{font-size:18px;page-break-after:avoid}table,img{page-break-inside:avoid}p,li{font-size:11px}table{font-size:8.5px}img{max-height:170mm}}
""".strip()
    (static_dir / "print.css").write_text(css + "\n", encoding="utf-8")
    html_path = static_dir / "index.html"
    subprocess.run(
        ["pandoc", str(md_path), "--standalone", "--metadata", "title=gdTAI V4.2 Step 0", "--css", "print.css", "-o", str(html_path)],
        cwd=ROOT,
        check=True,
    )
    pdf_path = static_dir / "gdtai_v4_2_step0_report.pdf"
    subprocess.run(
        [
            "google-chrome", "--headless", "--no-sandbox", "--disable-gpu",
            "--disable-dev-shm-usage", "--disable-breakpad", "--disable-crash-reporter",
            "--allow-file-access-from-files", "--no-pdf-header-footer",
            "--user-data-dir=/tmp/gdtai-v42-step0-chrome-profile",
            f"--print-to-pdf={pdf_path}", html_path.resolve().as_uri(),
        ],
        cwd=ROOT,
        check=True,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = json.loads(resolve(args.config).read_text())
    outputs = config["outputs"]
    table_dir = resolve(outputs["table_dir"])
    figure_dir = resolve(outputs["figure_dir"])
    log_dir = resolve(outputs["log_dir"])
    static_dir = resolve(outputs["static_dir"])
    for path in [table_dir, figure_dir, log_dir, static_dir]:
        path.mkdir(parents=True, exist_ok=True)

    cells, features, matrix, checks = read_inputs(config)
    provenance = annotation_provenance(cells)
    expression, expression_strength, coherent_t = expression_audit(cells, features, matrix, config)
    counterfactual, source_counterfactual, checkpoint_checks = counterfactual_audit(cells, coherent_t, config)
    checks.extend(checkpoint_checks)

    dual = counterfactual[counterfactual["variant"].eq("dual_annotation_primary")]
    checks.extend(
        [
            {
                "check": "stage1_checkpoints_complete",
                "status": "PASS" if dual.shape[0] == 15 else "FAIL",
                "detail": f"{dual.shape[0]}/15 candidate-fold combinations",
            },
            {
                "check": "dual_annotation_nk_sources",
                "status": "PASS" if int((provenance["dual_annotation_n"] > 0).sum()) >= 2 else "FAIL",
                "detail": f"{int((provenance['dual_annotation_n'] > 0).sum())} sources",
            },
            {
                "check": "dual_annotation_counterfactual_gate",
                "status": "PASS" if dual["wilson_upper_cap_pass"].all() else "FAIL",
                "detail": f"{int(dual['wilson_upper_cap_pass'].sum())}/{dual.shape[0]} pass Wilson-aware 50% cap",
            },
            {
                "check": "expression_used_for_truth",
                "status": "PASS" if config["nk_truth"]["expression_coherence_is_diagnostic_only"] else "FAIL",
                "detail": "diagnostic only; no labels changed",
            },
            {
                "check": "model_fitting_performed",
                "status": "PASS",
                "detail": "none; saved V4.1 OOF probabilities only",
            },
        ]
    )
    checks_frame = pd.DataFrame(checks)

    provenance.to_csv(table_dir / "nk_annotation_provenance_by_source.csv", index=False)
    expression.to_csv(table_dir / "nk_expression_coherence_by_source.csv", index=False)
    expression_strength.to_csv(table_dir / "nk_expression_coherence_by_strength.csv", index=False)
    counterfactual.to_csv(table_dir / "stage1_counterfactual_feasibility.csv", index=False)
    source_counterfactual.to_csv(table_dir / "stage1_counterfactual_source_passage.csv", index=False)
    checks_frame.to_csv(table_dir / "step0_checks.csv", index=False)
    build_figures(provenance, expression, counterfactual, figure_dir)

    nk = cells[cells["stage1_role"].eq("nk_negative")]
    summary = {
        "protocol_version": config["protocol_version"],
        "result": "PASS" if checks_frame["status"].eq("PASS").all() else "FAIL",
        "model_fitting_performed": False,
        "n_v4_1_nk_controls": int(nk.shape[0]),
        "n_single_annotation_nk": int((nk["nk_annotation_strength"] == 0.5).sum()),
        "n_dual_annotation_nk": int((nk["nk_annotation_strength"] == 1.0).sum()),
        "n_dual_annotation_sources": int((provenance["dual_annotation_n"] > 0).sum()),
        "n_checkpoint_candidate_folds": int(dual.shape[0]),
        "n_dual_counterfactual_wilson_gate_pass": int(dual["wilson_upper_cap_pass"].sum()),
        "dual_counterfactual_max_source_passage_range": [
            float(dual["maximum_source_passage_n_ge_minimum"].min()),
            float(dual["maximum_source_passage_n_ge_minimum"].max()),
        ],
        "all_control_max_source_passage_range": [
            float(counterfactual.loc[counterfactual["variant"].eq("all_v4_1_controls"), "maximum_source_passage_n_ge_minimum"].min()),
            float(counterfactual.loc[counterfactual["variant"].eq("all_v4_1_controls"), "maximum_source_passage_n_ge_minimum"].max()),
        ],
        "next_action": "Obtain explicit supervision approval before implementing or fitting V4.2.",
    }
    (log_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    render_report(cells, provenance, expression, counterfactual, checks_frame, config, log_dir, static_dir)
    if not checks_frame["status"].eq("PASS").all():
        raise SystemExit("V4.2 Step 0 audit failed")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
