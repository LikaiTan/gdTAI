#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Evaluate a transcriptome-gated cascade around the frozen gdTAI model.

This script is read-only with respect to the input H5AD. It does not retrain the
selected gdTAI TCR-gene model. It only tunes wrapper thresholds on the existing
non-holdout tune split, then evaluates the frozen-model cascade on the strict
multi-cohort validation split and the full 5 million cell atlas.
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
import json
import logging
import pickle
import shutil
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix

from run_gdt_prediction_package_evaluation import (
    DEFAULT_INPUT_H5AD,
    clean_group_values,
    dataframe_to_html,
    dataframe_to_markdown,
    normalize_strings,
    read_bool_obs,
    read_float_obs,
    read_obs_column,
    read_string_dataset,
)
from run_gdt_deg_tcr_classifier_training import json_ready, metric_dict, safe_div
from run_gdt_gse144469_holdout_tcrgene_classifier import (
    EXTRA_TRAB_HOLDOUT_SOURCE,
    GDT2020_HOLDOUT_TISSUE,
    GDT2020_SOURCE,
    HOLDOUT_SOURCE,
    RANDOM_SEED,
    SUBOPTIMAL_SOURCE,
    TARGET_SUM,
    TCR_PREFIXES,
    FeatureSpec,
    build_obs_metadata,
    extract_gene_matrix,
    load_original_trd_trab_threshold,
    read_nonempty_if_present,
    sublibrary_ids,
    summarize_full_predictions,
    tcr_family,
    tcr_priority,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
STATIC_DIR = PROJECT_ROOT / "gdT_prediction"
OUT_PREFIX = "gdtai_transcriptome_gate"

TABLE_DIR = OUTPUT_ROOT / "tables" / "gdT_prediction" / OUT_PREFIX
FIGURE_DIR = OUTPUT_ROOT / "figures" / "gdT_prediction" / OUT_PREFIX
LOG_DIR = OUTPUT_ROOT / "logs" / "gdT_prediction" / OUT_PREFIX
MODEL_DIR = OUTPUT_ROOT / "models" / "gdT_prediction_classifier" / OUT_PREFIX
STATIC_ASSET_DIR = STATIC_DIR / "assets" / OUT_PREFIX

REPORT_MD = LOG_DIR / "gdtai_transcriptome_gate_report.md"
REPORT_HTML = STATIC_DIR / "gdtai_transcriptome_gate_report.html"
REPORT_PDF = STATIC_DIR / "gdtai_transcriptome_gate_report.pdf"
RUN_LOG = LOG_DIR / "run.log"
SUMMARY_JSON = LOG_DIR / "gdtai_transcriptome_gate_summary.json"
WRAPPER_PKL = MODEL_DIR / "selected_transcriptome_gate_wrapper.pkl"

CURRENT_MODEL_PATH = (
    OUTPUT_ROOT / "models" / "gdT_prediction_classifier" / "gse144469_holdout_tcrgene" / "selected_model.pkl"
)
CURRENT_TABLE_DIR = OUTPUT_ROOT / "tables" / "gdT_prediction" / "gse144469_holdout_tcrgene"
NKGUARD_TABLE_DIR = OUTPUT_ROOT / "tables" / "gdT_prediction" / "gdtai_nkguard"
MULTIGUARD_TABLE_DIR = OUTPUT_ROOT / "tables" / "gdT_prediction" / "gdtai_multiguard"

FULL_APPLY_CHUNK = 50_000
SCATTER_SAMPLE_CELLS = 250_000
UMAP_BACKGROUND_CELLS = 100_000
UMAP_POSITIVE_CELLS = 60_000

HIGH_CONF_ANNOTATIONS = {"GDT_CELL", "CD8_T"}
CD4_ANNOTATION = "CD4_T"
RESTRICTED_ANNOTATIONS = {"TREG", "NK_CELL", "OTHER"}

T_CELL_MARKERS = ["CD3D", "CD3E", "CD3G", "CD247", "LCK"]
NK_MARKERS = ["NKG7", "GNLY", "PRF1", "KLRD1", "KLRF1", "FCGR3A", "NCR1", "TYROBP", "FCER1G", "CST7"]
DEATH_PENALTY_MARKERS = ["FOXP3", "IL2RA", "CTLA4", "CD4"]
NON_T_AUDIT_MARKERS = [
    "MS4A1",
    "CD79A",
    "CD14",
    "LYZ",
    "FCGR3A",
    "EPCAM",
    "KRT19",
    "COL1A1",
    "PECAM1",
]
MARKER_GENES = sorted(set(T_CELL_MARKERS + NK_MARKERS + DEATH_PENALTY_MARKERS + NON_T_AUDIT_MARKERS))

GATE_CODE_TO_NAME = {
    0: "below_threshold",
    1: "high_conf_pass",
    2: "cd4_strict_pass",
    3: "restricted_rescue_pass",
    4: "blocked_cd4_or_treg_penalty",
    5: "blocked_strong_nk_or_low_cd3",
    6: "blocked_by_transcriptome_gate",
}

ACCEPT_RECALL_MIN = 0.84
ACCEPT_VALIDATION_DELTA_RECALL_MIN = -0.01
ACCEPT_VALIDATION_DELTA_F1_MIN = -0.01
ACCEPT_FULL_NK_FP_REDUCTION_MIN = 0.50
ACCEPT_FULL_TCRAB_FP_MAX_VS_CURRENT = 1.00
ACCEPT_FULL_KNOWN_FP_FRACTION_DELTA_MAX = -1e-12


@dataclass
class LoadedModel:
    model_name: str
    model_object: Any
    threshold: float
    gene_names: list[str]
    feature_names: list[str]
    notes: str


@dataclass
class FeatureBundle:
    spec: FeatureSpec
    model_cols: np.ndarray
    marker_cols: dict[str, int]
    missing_marker_genes: list[str]


@dataclass(frozen=True)
class CascadeParams:
    high_conf_threshold: float
    cd4_threshold: float
    rescue_threshold: float
    cd3_min: float
    nk_max: float
    foxp3_max: float
    cd4_rescue_max: float


@dataclass
class SplitBundle:
    train_idx: np.ndarray
    tune_idx: np.ndarray
    validation_idx: np.ndarray
    excluded_idx: np.ndarray
    split_overall: pd.DataFrame
    split_by_source: pd.DataFrame


@dataclass
class ScoreBlock:
    rows: np.ndarray
    y_true: np.ndarray
    annotation: np.ndarray
    score: np.ndarray
    cd3_score: np.ndarray
    nk_score: np.ndarray
    foxp3_expr: np.ndarray
    cd4_expr: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate transcriptome-gated gdTAI cascade.")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--current-model-pkl", type=Path, default=CURRENT_MODEL_PATH)
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    parser.add_argument("--chunk-size", type=int, default=FULL_APPLY_CHUNK)
    parser.add_argument("--scatter-sample-cells", type=int, default=SCATTER_SAMPLE_CELLS)
    parser.add_argument("--umap-background-cells", type=int, default=UMAP_BACKGROUND_CELLS)
    parser.add_argument("--umap-positive-cells", type=int, default=UMAP_POSITIVE_CELLS)
    parser.add_argument("--skip-full-apply", action="store_true")
    parser.add_argument("--no-pdf", action="store_true")
    return parser.parse_args()


def ensure_dirs() -> None:
    for path in [TABLE_DIR, FIGURE_DIR, LOG_DIR, MODEL_DIR, STATIC_DIR, STATIC_ASSET_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def setup_logging() -> None:
    ensure_dirs()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(RUN_LOG, mode="w", encoding="utf-8"), logging.StreamHandler()],
        force=True,
    )


def load_current_model(path: Path) -> LoadedModel:
    with path.open("rb") as handle:
        payload = pickle.load(handle)
    required = ["model", "model_object", "threshold", "gene_names", "feature_names"]
    missing = [key for key in required if key not in payload]
    if missing:
        raise KeyError(f"Selected gdTAI model pickle is missing keys: {missing}")
    return LoadedModel(
        model_name=str(payload["model"]),
        model_object=payload["model_object"],
        threshold=float(payload["threshold"]),
        gene_names=[str(gene) for gene in payload["gene_names"]],
        feature_names=[str(feature) for feature in payload["feature_names"]],
        notes=str(payload.get("notes", "")),
    )


def build_feature_bundle(handle: h5py.File, model: LoadedModel) -> FeatureBundle:
    var_names = pd.Index(read_string_dataset(handle["var"]["_index"]), dtype="string").astype(str).tolist()
    gene_to_idx = {gene: idx for idx, gene in enumerate(var_names)}
    missing_model = [gene for gene in model.gene_names if gene not in gene_to_idx]
    if missing_model:
        raise KeyError(f"Input H5AD is missing selected model genes: {missing_model[:20]}")

    selected: list[str] = []
    for gene in [*model.gene_names, *MARKER_GENES]:
        if gene in gene_to_idx and gene not in selected:
            selected.append(gene)
    name_to_col = {gene: idx for idx, gene in enumerate(selected)}
    marker_cols = {gene: name_to_col[gene] for gene in MARKER_GENES if gene in name_to_col}
    missing_markers = [gene for gene in MARKER_GENES if gene not in marker_cols]
    model_cols = np.asarray([name_to_col[gene] for gene in model.gene_names], dtype=np.int32)

    def feature_family(gene: str) -> str:
        if gene in T_CELL_MARKERS:
            return "T_cell_marker"
        if gene in NK_MARKERS:
            return "NK_marker"
        if gene in DEATH_PENALTY_MARKERS:
            return "death_penalty_marker"
        if gene in NON_T_AUDIT_MARKERS:
            return "non_T_audit_marker"
        return tcr_family(gene)

    cd3_cols = [name_to_col[gene] for gene in ["CD3D", "CD3E", "CD3G"] if gene in name_to_col]
    spec = FeatureSpec(
        feature_names=[f"{gene}_log1p_cp10k" for gene in selected],
        gene_names=selected,
        gene_indices=np.asarray([gene_to_idx[gene] for gene in selected], dtype=np.int32),
        tcr_gene_names=[gene for gene in selected if gene.upper().startswith(TCR_PREFIXES)],
        penalty_gene_names=[gene for gene in selected if gene in DEATH_PENALTY_MARKERS or gene in T_CELL_MARKERS],
        family_by_feature=[feature_family(gene) for gene in selected],
        cd3_cols=cd3_cols,
        cd4_col=name_to_col.get("CD4"),
        foxp3_col=name_to_col.get("FOXP3"),
    )
    manifest = pd.DataFrame(
        {
            "feature_index": np.arange(len(selected), dtype=int),
            "gene": selected,
            "feature": spec.feature_names,
            "family": spec.family_by_feature,
            "used_by_frozen_model": [gene in set(model.gene_names) for gene in selected],
            "used_by_transcriptome_gate": [gene in set(MARKER_GENES) for gene in selected],
        }
    )
    manifest.to_csv(TABLE_DIR / "feature_manifest.csv", index=False)
    pd.DataFrame({"missing_marker_gene": missing_markers}).to_csv(TABLE_DIR / "missing_marker_genes.csv", index=False)
    return FeatureBundle(spec=spec, model_cols=model_cols, marker_cols=marker_cols, missing_marker_genes=missing_markers)


def build_splits(obs: Any, seed: int) -> SplitBundle:
    primary = (obs.class_code == 1) | (obs.class_code == 2)
    validation_primary_idx = np.flatnonzero(primary & (obs.source == HOLDOUT_SOURCE)).astype(np.int64)
    validation_gdt2020_idx = np.flatnonzero(
        primary
        & (obs.source == GDT2020_SOURCE)
        & (obs.class_code == 2)
        & (pd.Series(obs.tissue, copy=False).astype(str).str.lower().to_numpy() == GDT2020_HOLDOUT_TISSUE.lower())
    ).astype(np.int64)
    validation_extra_trab_idx = np.flatnonzero(
        (obs.source == EXTRA_TRAB_HOLDOUT_SOURCE)
        & obs.has_TRA_TRB_paired
        & (~obs.corrected_has_any_gd_tcr)
    ).astype(np.int64)
    validation_idx = np.unique(
        np.concatenate([validation_primary_idx, validation_gdt2020_idx, validation_extra_trab_idx])
    ).astype(np.int64)
    validation_mask = np.zeros(obs.source.shape[0], dtype=bool)
    validation_mask[validation_idx] = True
    excluded_idx = np.flatnonzero(primary & (obs.source == SUBOPTIMAL_SOURCE)).astype(np.int64)
    train_tune_pool = np.flatnonzero(
        primary
        & (~validation_mask)
        & (obs.source != HOLDOUT_SOURCE)
        & (obs.source != SUBOPTIMAL_SOURCE)
        & (obs.source != EXTRA_TRAB_HOLDOUT_SOURCE)
    ).astype(np.int64)
    if validation_primary_idx.size == 0:
        raise RuntimeError(f"No primary gold cells found in holdout source {HOLDOUT_SOURCE}")
    if validation_gdt2020_idx.size == 0:
        raise RuntimeError(f"No {GDT2020_SOURCE} gdT cells found with tissue `{GDT2020_HOLDOUT_TISSUE}`")
    if validation_extra_trab_idx.size == 0:
        raise RuntimeError(f"No paired TCRAB/no-gdTCR cells found in {EXTRA_TRAB_HOLDOUT_SOURCE}")

    rng = np.random.default_rng(seed)
    train_parts: list[np.ndarray] = []
    tune_parts: list[np.ndarray] = []
    by_source_rows: list[dict[str, Any]] = []
    pool_df = pd.DataFrame(
        {
            "obs_index": train_tune_pool,
            "source_gse_id": obs.source[train_tune_pool],
            "label": np.where(obs.class_code[train_tune_pool] == 2, "gdT_gold", "abT_gold"),
        }
    )
    for (source, label), group in pool_df.groupby(["source_gse_id", "label"], sort=True):
        idx = group["obs_index"].to_numpy(dtype=np.int64)
        rng.shuffle(idx)
        n_tune = 0 if idx.size < 5 else max(1, int(round(idx.size * 0.20)))
        if n_tune:
            tune_parts.append(idx[:n_tune])
            train_parts.append(idx[n_tune:])
        else:
            train_parts.append(idx)
        by_source_rows.append(
            {
                "source_gse_id": source,
                "label": label,
                "n_cells": int(idx.size),
                "train": int(idx.size - n_tune),
                "tune": int(n_tune),
            }
        )
    train_idx = np.concatenate(train_parts).astype(np.int64) if train_parts else np.asarray([], dtype=np.int64)
    tune_idx = np.concatenate(tune_parts).astype(np.int64) if tune_parts else np.asarray([], dtype=np.int64)

    overall_rows: list[dict[str, Any]] = []
    for split, idx in [
        ("train", train_idx),
        ("tune", tune_idx),
        (f"validation_primary_{HOLDOUT_SOURCE}", validation_primary_idx),
        (f"validation_gdT_{GDT2020_SOURCE}_{GDT2020_HOLDOUT_TISSUE.replace(' ', '_')}", validation_gdt2020_idx),
        (f"validation_abT_{EXTRA_TRAB_HOLDOUT_SOURCE}_paired_TCRAB_no_gdTCR", validation_extra_trab_idx),
        ("validation_combined", validation_idx),
        (f"sensitivity_excluded_{SUBOPTIMAL_SOURCE}", excluded_idx),
    ]:
        labels = obs.class_code[idx]
        if split.startswith("validation_abT_") or split == "validation_combined":
            gdT_gold = int((labels == 2).sum())
            abT_gold = int(idx.size - gdT_gold)
        else:
            gdT_gold = int((labels == 2).sum())
            abT_gold = int((labels == 1).sum())
        overall_rows.append(
            {
                "split": split,
                "n_cells": int(idx.size),
                "gdT_gold": gdT_gold,
                "abT_gold": abT_gold,
                "gdT_prevalence": safe_div(gdT_gold, int(idx.size)),
            }
        )
    split_overall = pd.DataFrame(overall_rows)
    split_by_source = pd.DataFrame(by_source_rows)
    split_overall.to_csv(TABLE_DIR / "split_overall.csv", index=False)
    split_by_source.to_csv(TABLE_DIR / "split_by_source_label.csv", index=False)
    return SplitBundle(
        train_idx=train_idx,
        tune_idx=tune_idx,
        validation_idx=validation_idx,
        excluded_idx=excluded_idx,
        split_overall=split_overall,
        split_by_source=split_by_source,
    )


def normalize_annotation(annotation: np.ndarray) -> np.ndarray:
    return pd.Series(annotation, copy=False).astype(str).str.upper().str.strip().to_numpy(dtype=object)


def sum_cols(x: np.ndarray, cols: list[int]) -> np.ndarray:
    if not cols:
        return np.zeros(x.shape[0], dtype=np.float32)
    return x[:, cols].sum(axis=1).astype(np.float32, copy=False)


def mean_cols(x: np.ndarray, cols: list[int]) -> np.ndarray:
    if not cols:
        return np.zeros(x.shape[0], dtype=np.float32)
    return x[:, cols].mean(axis=1).astype(np.float32, copy=False)


def marker_scores(x: np.ndarray, bundle: FeatureBundle) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    cd3_cols = [bundle.marker_cols[gene] for gene in T_CELL_MARKERS if gene in bundle.marker_cols]
    nk_cols = [bundle.marker_cols[gene] for gene in NK_MARKERS if gene in bundle.marker_cols]
    cd3_score = sum_cols(x, cd3_cols)
    nk_score = mean_cols(x, nk_cols)
    foxp3_expr = x[:, bundle.marker_cols["FOXP3"]].astype(np.float32, copy=False) if "FOXP3" in bundle.marker_cols else np.zeros(x.shape[0], dtype=np.float32)
    cd4_expr = x[:, bundle.marker_cols["CD4"]].astype(np.float32, copy=False) if "CD4" in bundle.marker_cols else np.zeros(x.shape[0], dtype=np.float32)
    return cd3_score, nk_score, foxp3_expr, cd4_expr


def predict_model_score(model: LoadedModel, x_model: np.ndarray) -> np.ndarray:
    return model.model_object.predict_proba(x_model)[:, 1].astype(np.float32)


def score_rows(
    handle: h5py.File,
    rows: np.ndarray,
    annotation: np.ndarray,
    y_true: np.ndarray,
    model: LoadedModel,
    bundle: FeatureBundle,
    *,
    label: str,
) -> ScoreBlock:
    x = extract_gene_matrix(handle, rows, bundle.spec, label=label)
    score = predict_model_score(model, x[:, bundle.model_cols])
    cd3_score, nk_score, foxp3_expr, cd4_expr = marker_scores(x, bundle)
    return ScoreBlock(
        rows=rows,
        y_true=y_true,
        annotation=annotation[rows],
        score=score,
        cd3_score=cd3_score,
        nk_score=nk_score,
        foxp3_expr=foxp3_expr,
        cd4_expr=cd4_expr,
    )


def apply_cascade(block: ScoreBlock, params: CascadeParams) -> tuple[np.ndarray, np.ndarray]:
    ann = normalize_annotation(block.annotation)
    high_conf = np.isin(ann, list(HIGH_CONF_ANNOTATIONS))
    cd4 = ann == CD4_ANNOTATION
    restricted = np.isin(ann, list(RESTRICTED_ANNOTATIONS)) | (~high_conf & ~cd4)
    low_cd3 = block.cd3_score < params.cd3_min
    strong_nk_or_low_cd3 = (block.nk_score > params.nk_max) | low_cd3
    cd4_or_treg_penalty = (block.foxp3_expr > params.foxp3_max) | (restricted & (block.cd4_expr > params.cd4_rescue_max))

    pred = np.zeros(block.score.size, dtype=bool)
    code = np.zeros(block.score.size, dtype=np.uint8)
    high_pass = high_conf & (block.score >= params.high_conf_threshold)
    cd4_pass = cd4 & (block.score >= params.cd4_threshold) & (block.foxp3_expr <= params.foxp3_max)
    rescue_pass = (
        restricted
        & (block.score >= params.rescue_threshold)
        & (~strong_nk_or_low_cd3)
        & (~cd4_or_treg_penalty)
    )
    pred[high_pass | cd4_pass | rescue_pass] = True
    code[high_pass] = 1
    code[cd4_pass] = 2
    code[rescue_pass] = 3

    base_positive = block.score >= params.high_conf_threshold
    code[(~pred) & base_positive & cd4_or_treg_penalty] = 4
    code[(~pred) & base_positive & strong_nk_or_low_cd3] = 5
    code[(~pred) & base_positive & (code == 0)] = 6
    return pred, code


def metrics_for_strategy(name: str, y_true: np.ndarray, score: np.ndarray, pred: np.ndarray, threshold: float) -> dict[str, Any]:
    out = {
        "strategy": name,
        "threshold": float(threshold),
    }
    out.update(metric_dict(y_true, pred.astype(np.int8), score))
    return out


def tune_cascade(block: ScoreBlock, model: LoadedModel) -> tuple[CascadeParams, pd.DataFrame, pd.DataFrame]:
    current_pred = block.score >= model.threshold
    current_metrics = metrics_for_strategy("current_gdTAI", block.y_true, block.score, current_pred, model.threshold)
    positive_nk = block.nk_score[block.y_true == 1]
    nk_reference = float(np.nanpercentile(positive_nk, 95)) if positive_nk.size else 1.0
    nk_candidates = sorted(set([max(nk_reference * 0.75, 0.25), max(nk_reference, 0.5), np.inf]))
    cd3_candidates = [0.0, 0.25, 0.5]
    high_candidates = sorted(set([model.threshold, max(model.threshold, 0.93), 0.95, 0.97]))
    cd4_candidates = sorted(set([model.threshold, max(model.threshold, 0.95), 0.97, 0.985]))
    rescue_candidates = sorted(set([max(model.threshold, 0.97), 0.985, 0.995]))
    foxp3_candidates = [0.25, 0.5]
    cd4_rescue_candidates = [0.75, 1.25]

    recall_floor = max(ACCEPT_RECALL_MIN, float(current_metrics["recall"]) - 0.02)
    f1_floor = float(current_metrics["f1"]) - 0.01
    rows: list[dict[str, Any]] = []
    for high in high_candidates:
        for cd4 in cd4_candidates:
            if cd4 < high:
                continue
            for rescue in rescue_candidates:
                if rescue < cd4:
                    continue
                for cd3_min in cd3_candidates:
                    for nk_max in nk_candidates:
                        for foxp3_max in foxp3_candidates:
                            for cd4_rescue_max in cd4_rescue_candidates:
                                params = CascadeParams(
                                    high_conf_threshold=float(high),
                                    cd4_threshold=float(cd4),
                                    rescue_threshold=float(rescue),
                                    cd3_min=float(cd3_min),
                                    nk_max=float(nk_max),
                                    foxp3_max=float(foxp3_max),
                                    cd4_rescue_max=float(cd4_rescue_max),
                                )
                                pred, code = apply_cascade(block, params)
                                metrics = metrics_for_strategy("cascade_tune_candidate", block.y_true, block.score, pred, high)
                                row = {
                                    **asdict(params),
                                    **{f"tune_{k}": v for k, v in metrics.items() if k not in {"strategy", "threshold"}},
                                    "current_tune_f1": current_metrics["f1"],
                                    "current_tune_recall": current_metrics["recall"],
                                    "current_tune_precision": current_metrics["precision"],
                                    "passes_recall_floor": metrics["recall"] >= recall_floor,
                                    "passes_f1_floor": metrics["f1"] >= f1_floor,
                                    "n_high_conf_pass": int((code == 1).sum()),
                                    "n_cd4_strict_pass": int((code == 2).sum()),
                                    "n_restricted_rescue_pass": int((code == 3).sum()),
                                    "n_blocked_cd4_or_treg_penalty": int((code == 4).sum()),
                                    "n_blocked_strong_nk_or_low_cd3": int((code == 5).sum()),
                                    "n_blocked_by_transcriptome_gate": int((code == 6).sum()),
                                }
                                rows.append(row)
    grid = pd.DataFrame(rows)
    valid = grid["passes_recall_floor"] & grid["passes_f1_floor"]
    if valid.any():
        ranked = grid.loc[valid].sort_values(
            ["tune_fp", "tune_precision", "tune_f1", "tune_recall"],
            ascending=[True, False, False, False],
        )
        reason = "recall_and_f1_floor_then_lowest_tune_fp"
    else:
        ranked = grid.sort_values(["tune_f1", "tune_precision", "tune_recall"], ascending=[False, False, False])
        reason = "fallback_best_tune_f1_no_floor_candidate"
    chosen = ranked.iloc[0]
    params = CascadeParams(
        high_conf_threshold=float(chosen["high_conf_threshold"]),
        cd4_threshold=float(chosen["cd4_threshold"]),
        rescue_threshold=float(chosen["rescue_threshold"]),
        cd3_min=float(chosen["cd3_min"]),
        nk_max=float(chosen["nk_max"]),
        foxp3_max=float(chosen["foxp3_max"]),
        cd4_rescue_max=float(chosen["cd4_rescue_max"]),
    )
    grid.insert(0, "selection_reason", "")
    grid.loc[grid.index == chosen.name, "selection_reason"] = reason
    grid.sort_values(
        ["passes_recall_floor", "passes_f1_floor", "tune_fp", "tune_precision", "tune_f1"],
        ascending=[False, False, True, False, False],
    ).to_csv(TABLE_DIR / "cascade_threshold_tuning_grid.csv", index=False)
    current_df = pd.DataFrame([current_metrics])
    current_df.to_csv(TABLE_DIR / "current_model_tune_metrics.csv", index=False)
    return params, grid, current_df


def evaluate_blocks(
    tune_block: ScoreBlock,
    validation_block: ScoreBlock,
    params: CascadeParams,
    model: LoadedModel,
    handle: h5py.File,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    current_tune_pred = tune_block.score >= model.threshold
    cascade_tune_pred, tune_code = apply_cascade(tune_block, params)
    current_valid_pred = validation_block.score >= model.threshold
    cascade_valid_pred, valid_code = apply_cascade(validation_block, params)

    original_threshold = load_original_trd_trab_threshold()
    trd_minus = read_float_obs(handle, "phase4_trd_minus_trab")
    original_tune_score = trd_minus[tune_block.rows]
    original_valid_score = trd_minus[validation_block.rows]
    original_tune_pred = original_tune_score >= original_threshold
    original_valid_pred = original_valid_score >= original_threshold

    tune_rows = [
        metrics_for_strategy("current_gdTAI", tune_block.y_true, tune_block.score, current_tune_pred, model.threshold),
        metrics_for_strategy("transcriptome_gated_cascade", tune_block.y_true, tune_block.score, cascade_tune_pred, params.high_conf_threshold),
        metrics_for_strategy("original_TRD_minus_TRAB_score_strategy", tune_block.y_true, original_tune_score, original_tune_pred, original_threshold),
    ]
    validation_rows = [
        metrics_for_strategy("current_gdTAI", validation_block.y_true, validation_block.score, current_valid_pred, model.threshold),
        metrics_for_strategy("transcriptome_gated_cascade", validation_block.y_true, validation_block.score, cascade_valid_pred, params.high_conf_threshold),
        metrics_for_strategy("original_TRD_minus_TRAB_score_strategy", validation_block.y_true, original_valid_score, original_valid_pred, original_threshold),
    ]
    tune_df = pd.DataFrame(tune_rows)
    validation_df = pd.DataFrame(validation_rows)
    tune_df.to_csv(TABLE_DIR / "tune_metrics.csv", index=False)
    validation_df.to_csv(TABLE_DIR / "validation_metrics.csv", index=False)

    gate_summary = pd.DataFrame(
        {
            "split": np.repeat(["tune", "validation"], [tune_code.size, valid_code.size]),
            "gate_code": np.concatenate([tune_code, valid_code]).astype(int),
            "gate_status": [GATE_CODE_TO_NAME[int(code)] for code in np.concatenate([tune_code, valid_code])],
            "label": np.concatenate([
                np.where(tune_block.y_true == 1, "gdT_gold", "abT_gold"),
                np.where(validation_block.y_true == 1, "gdT_gold", "abT_gold"),
            ]),
        }
    )
    (
        gate_summary.groupby(["split", "label", "gate_code", "gate_status"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .to_csv(TABLE_DIR / "gate_status_by_split_label.csv", index=False)
    )
    return tune_df, validation_df


def build_full_metadata(handle: h5py.File) -> dict[str, np.ndarray]:
    n_obs = int(handle["obs"]["_index"].shape[0])
    source = clean_group_values(read_obs_column(handle, "source_gse_id"))
    tissue_key = "tissue_corrected" if "tissue_corrected" in handle["obs"] else "tissue"
    tissue = clean_group_values(read_obs_column(handle, tissue_key))
    annotation = clean_group_values(read_obs_column(handle, "simple_annotation_plus6"))
    library_id = normalize_strings(read_obs_column(handle, "library_id"))
    sample_id = normalize_strings(read_obs_column(handle, "sample_id"))
    sublibrary = sublibrary_ids(library_id, sample_id)
    has_TRA_TRB_paired = read_bool_obs(handle, "has_TRA_TRB_paired")
    has_any_ab_tcr = read_bool_obs(handle, "has_any_ab_tcr") if "has_any_ab_tcr" in handle["obs"] else has_TRA_TRB_paired.copy()
    tra_nonempty = read_nonempty_if_present(handle, "TRA_cdr3", n_obs)
    trb_nonempty = read_nonempty_if_present(handle, "TRB_cdr3", n_obs)
    ab_evidence = has_TRA_TRB_paired | has_any_ab_tcr | tra_nonempty | trb_nonempty
    keys = pd.Series(source, copy=False).astype(str) + "||" + pd.Series(sublibrary, copy=False).astype(str)
    tcrab_seq_keys = set(keys[ab_evidence])
    is_tcrab_seq_library = keys.isin(tcrab_seq_keys).to_numpy(dtype=bool)
    is_nk = pd.Series(annotation, copy=False).astype(str).str.upper().str.contains("NK", regex=False).to_numpy(dtype=bool)
    return {
        "source": source,
        "tissue": tissue,
        "annotation": annotation,
        "sublibrary": sublibrary,
        "is_tcrab_seq_library": is_tcrab_seq_library,
        "has_TRA_TRB_paired": has_TRA_TRB_paired,
        "is_nk": is_nk,
    }


def apply_full_atlas(
    handle: h5py.File,
    model: LoadedModel,
    bundle: FeatureBundle,
    params: CascadeParams,
    *,
    chunk_size: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, np.ndarray]]:
    n_obs = int(handle["obs"]["_index"].shape[0])
    meta = build_full_metadata(handle)
    current_score = np.zeros(n_obs, dtype=np.float32)
    current_pred = np.zeros(n_obs, dtype=bool)
    cascade_pred = np.zeros(n_obs, dtype=bool)
    gate_code = np.zeros(n_obs, dtype=np.uint8)

    for start in range(0, n_obs, chunk_size):
        end = min(start + chunk_size, n_obs)
        rows = np.arange(start, end, dtype=np.int64)
        x = extract_gene_matrix(handle, rows, bundle.spec, label=f"{OUT_PREFIX}_full_{start}_{end}")
        score = predict_model_score(model, x[:, bundle.model_cols])
        cd3_score, nk_score, foxp3_expr, cd4_expr = marker_scores(x, bundle)
        block = ScoreBlock(
            rows=rows,
            y_true=np.zeros(rows.size, dtype=np.int8),
            annotation=meta["annotation"][rows],
            score=score,
            cd3_score=cd3_score,
            nk_score=nk_score,
            foxp3_expr=foxp3_expr,
            cd4_expr=cd4_expr,
        )
        pred, code = apply_cascade(block, params)
        current_score[start:end] = score
        current_pred[start:end] = score >= model.threshold
        cascade_pred[start:end] = pred
        gate_code[start:end] = code
        if start and start % 500_000 == 0:
            logging.info("Applied transcriptome-gated cascade to %s / %s cells", f"{start:,}", f"{n_obs:,}")

    original_threshold = load_original_trd_trab_threshold()
    original_score = read_float_obs(handle, "phase4_trd_minus_trab")
    original_pred = original_score >= original_threshold
    summaries = []
    for strategy, threshold, pred in [
        ("current_gdTAI", model.threshold, current_pred),
        ("transcriptome_gated_cascade", params.high_conf_threshold, cascade_pred),
        ("original_TRD_minus_TRAB_score_strategy", original_threshold, original_pred),
    ]:
        summaries.append(
            summarize_full_predictions(
                strategy=strategy,
                threshold=float(threshold),
                pred=pred,
                source=meta["source"],
                tissue=meta["tissue"],
                annotation=meta["annotation"],
                sublibrary=meta["sublibrary"],
                is_tcrab_seq_library=meta["is_tcrab_seq_library"],
                has_TRA_TRB_paired=meta["has_TRA_TRB_paired"],
                is_nk=meta["is_nk"],
            )
        )
    overall = pd.concat([item[0] for item in summaries], ignore_index=True)
    source_df = pd.concat([item[1] for item in summaries], ignore_index=True)
    tissue_df = pd.concat([item[2] for item in summaries], ignore_index=True)
    fp_source_df = pd.concat([item[3] for item in summaries], ignore_index=True)
    sublibrary_df = pd.concat([item[4] for item in summaries], ignore_index=True)
    annotation_df = pd.concat([item[5] for item in summaries], ignore_index=True)
    fp_overall = overall[
        [
            "strategy",
            "threshold",
            "total_cells",
            "predicted_putative_gdT",
            "predicted_fraction",
            "tcrab_seq_library_cells",
            "predicted_in_tcrab_seq_libraries",
            "paired_TCRAB_cells",
            "predicted_paired_TCRAB_FP",
            "paired_TCRAB_FP_ratio_among_paired_TCRAB",
            "paired_TCRAB_FP_fraction_of_predictions",
            "NK_cells",
            "predicted_NK_FP",
            "NK_FP_ratio_among_NK",
            "NK_FP_fraction_of_predictions",
            "predicted_paired_TCRAB_or_NK_FP",
            "known_FP_fraction_of_predictions",
        ]
    ].copy()

    gate_full = pd.DataFrame(
        {
            "gate_code": gate_code.astype(int),
            "gate_status": [GATE_CODE_TO_NAME[int(code)] for code in gate_code],
            "current_gdTAI_predicted": current_pred,
            "cascade_predicted": cascade_pred,
            "annotation": meta["annotation"],
        }
    )
    (
        gate_full.groupby(["annotation", "gate_code", "gate_status", "current_gdTAI_predicted", "cascade_predicted"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .to_csv(TABLE_DIR / "full_atlas_gate_status_by_annotation.csv", index=False)
    )
    overall.to_csv(TABLE_DIR / "full_dataset_prediction_overall.csv", index=False)
    source_df.to_csv(TABLE_DIR / "full_dataset_prediction_by_source.csv", index=False)
    tissue_df.to_csv(TABLE_DIR / "full_dataset_prediction_by_tissue.csv", index=False)
    fp_overall.to_csv(TABLE_DIR / "full_dataset_false_positive_overall.csv", index=False)
    fp_source_df.to_csv(TABLE_DIR / "full_dataset_false_positive_by_source.csv", index=False)
    sublibrary_df.to_csv(TABLE_DIR / "full_dataset_false_positive_by_sublibrary.csv", index=False)
    annotation_df.to_csv(TABLE_DIR / "full_dataset_prediction_by_annotation.csv", index=False)
    arrays = {
        "current_score": current_score,
        "current_pred": current_pred,
        "cascade_pred": cascade_pred,
        "gate_code": gate_code,
        "original_pred": original_pred,
    }
    np.savez(
        TABLE_DIR / "full_prediction_arrays.npz",
        current_score=current_score,
        current_pred=current_pred,
        cascade_pred=cascade_pred,
        gate_code=gate_code,
        original_pred=original_pred,
    )
    logging.info("Saved full prediction arrays for report reuse")
    return overall, source_df, tissue_df, fp_overall, fp_source_df, annotation_df, arrays


def compact_numeric(df: pd.DataFrame, digits: int = 4) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_float_dtype(out[col]):
            out[col] = out[col].map(lambda value: "" if pd.isna(value) else round(float(value), digits))
    return out


def pick_columns(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    if df.empty:
        return df
    return df[[col for col in columns if col in df.columns]].copy()


def add_guard_comparisons(full_fp_overall: pd.DataFrame) -> pd.DataFrame:
    frames = [full_fp_overall.copy()]
    for label, path in [
        ("nkguard", NKGUARD_TABLE_DIR / "full_dataset_false_positive_overall.csv"),
        ("multiguard", MULTIGUARD_TABLE_DIR / "full_dataset_false_positive_overall.csv"),
    ]:
        if not path.exists():
            continue
        df = pd.read_csv(path)
        df = df.loc[~df["strategy"].astype(str).isin({"selected_model_logistic_individual_TCRABGD_genes", "original_TRD_minus_TRAB_score_strategy"})].copy()
        df["comparison_source"] = label
        frames.append(df)
    out = pd.concat(frames, ignore_index=True, sort=False)
    out.to_csv(TABLE_DIR / "full_dataset_false_positive_comparison_all_strategies.csv", index=False)
    return out


def evaluate_acceptance(validation_df: pd.DataFrame, full_fp_overall: pd.DataFrame) -> tuple[bool, pd.DataFrame]:
    current_valid = validation_df.loc[validation_df["strategy"] == "current_gdTAI"].iloc[0]
    cascade_valid = validation_df.loc[validation_df["strategy"] == "transcriptome_gated_cascade"].iloc[0]
    current_full = full_fp_overall.loc[full_fp_overall["strategy"] == "current_gdTAI"].iloc[0]
    cascade_full = full_fp_overall.loc[full_fp_overall["strategy"] == "transcriptome_gated_cascade"].iloc[0]
    checks = [
        {
            "gate": "validation_recall_min",
            "observed": float(cascade_valid["recall"]),
            "required": ACCEPT_RECALL_MIN,
            "passed": float(cascade_valid["recall"]) >= ACCEPT_RECALL_MIN,
        },
        {
            "gate": "validation_recall_delta_vs_current",
            "observed": float(cascade_valid["recall"] - current_valid["recall"]),
            "required": ACCEPT_VALIDATION_DELTA_RECALL_MIN,
            "passed": float(cascade_valid["recall"] - current_valid["recall"]) >= ACCEPT_VALIDATION_DELTA_RECALL_MIN,
        },
        {
            "gate": "validation_f1_delta_vs_current",
            "observed": float(cascade_valid["f1"] - current_valid["f1"]),
            "required": ACCEPT_VALIDATION_DELTA_F1_MIN,
            "passed": float(cascade_valid["f1"] - current_valid["f1"]) >= ACCEPT_VALIDATION_DELTA_F1_MIN,
        },
        {
            "gate": "full_nk_fp_reduction_vs_current",
            "observed": 1.0 - safe_div(float(cascade_full["predicted_NK_FP"]), float(current_full["predicted_NK_FP"])),
            "required": ACCEPT_FULL_NK_FP_REDUCTION_MIN,
            "passed": 1.0 - safe_div(float(cascade_full["predicted_NK_FP"]), float(current_full["predicted_NK_FP"])) >= ACCEPT_FULL_NK_FP_REDUCTION_MIN,
        },
        {
            "gate": "full_paired_tcrab_fp_no_worse_than_current",
            "observed": safe_div(float(cascade_full["predicted_paired_TCRAB_FP"]), float(current_full["predicted_paired_TCRAB_FP"])),
            "required": ACCEPT_FULL_TCRAB_FP_MAX_VS_CURRENT,
            "passed": float(cascade_full["predicted_paired_TCRAB_FP"]) <= float(current_full["predicted_paired_TCRAB_FP"]),
        },
        {
            "gate": "full_known_fp_fraction_lower_than_current",
            "observed": float(cascade_full["known_FP_fraction_of_predictions"] - current_full["known_FP_fraction_of_predictions"]),
            "required": ACCEPT_FULL_KNOWN_FP_FRACTION_DELTA_MAX,
            "passed": float(cascade_full["known_FP_fraction_of_predictions"]) < float(current_full["known_FP_fraction_of_predictions"]),
        },
    ]
    out = pd.DataFrame(checks)
    out.to_csv(TABLE_DIR / "promotion_acceptance_gates.csv", index=False)
    return bool(out["passed"].all()), out


def save_wrapper_if_accepted(accepted: bool, model_path: Path, model: LoadedModel, params: CascadeParams, bundle: FeatureBundle) -> Path | None:
    if not accepted:
        return None
    payload = {
        "wrapper_name": "gdTAI_transcriptome_gated_cascade",
        "base_model_path": str(model_path),
        "base_model": {
            "model": model.model_name,
            "model_object": model.model_object,
            "threshold": model.threshold,
            "gene_names": model.gene_names,
            "feature_names": model.feature_names,
            "notes": model.notes,
        },
        "cascade_params": asdict(params),
        "annotation_policy": {
            "high_confidence_pass": sorted(HIGH_CONF_ANNOTATIONS),
            "conditional_cd4": CD4_ANNOTATION,
            "rescue_only_or_restricted": sorted(RESTRICTED_ANNOTATIONS),
            "gate_status_codes": GATE_CODE_TO_NAME,
        },
        "marker_genes": {
            "t_cell_markers": T_CELL_MARKERS,
            "nk_markers": NK_MARKERS,
            "death_penalty_markers": DEATH_PENALTY_MARKERS,
            "missing_marker_genes": bundle.missing_marker_genes,
        },
    }
    with WRAPPER_PKL.open("wb") as handle:
        pickle.dump(payload, handle)
    return WRAPPER_PKL


def plot_validation_metric_bars(validation_df: pd.DataFrame) -> Path:
    path = FIGURE_DIR / "validation_metric_comparison.png"
    plot_df = validation_df.loc[
        validation_df["strategy"].isin(["current_gdTAI", "transcriptome_gated_cascade", "original_TRD_minus_TRAB_score_strategy"])
    ].copy()
    metrics = ["precision", "recall", "specificity", "f1"]
    fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    x = np.arange(len(metrics))
    width = 0.24
    colors = {
        "current_gdTAI": "#2563eb",
        "transcriptome_gated_cascade": "#0f766e",
        "original_TRD_minus_TRAB_score_strategy": "#d97706",
    }
    for i, (_, row) in enumerate(plot_df.iterrows()):
        ax.bar(x + (i - 1) * width, [row[m] for m in metrics], width=width, label=row["strategy"], color=colors.get(row["strategy"], "#64748b"))
    ax.set_xticks(x, labels=metrics)
    ax.set_ylim(0, 1.04)
    ax.set_ylabel("Metric value")
    ax.set_title("Strict validation performance")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_confusion_pair(validation_block: ScoreBlock, params: CascadeParams, model: LoadedModel) -> Path:
    path = FIGURE_DIR / "validation_confusion_current_vs_cascade.png"
    current_pred = validation_block.score >= model.threshold
    cascade_pred, _ = apply_cascade(validation_block, params)
    cms = [
        ("Current gdTAI", confusion_matrix(validation_block.y_true, current_pred.astype(int), labels=[0, 1])),
        ("Transcriptome-gated cascade", confusion_matrix(validation_block.y_true, cascade_pred.astype(int), labels=[0, 1])),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.8), constrained_layout=True)
    vmax = max(int(cm.max()) for _, cm in cms)
    for ax, (title, cm) in zip(axes, cms):
        ax.imshow(cm, cmap="Blues", vmin=0, vmax=vmax)
        ax.set_xticks([0, 1], labels=["Pred abT", "Pred gdT"])
        ax.set_yticks([0, 1], labels=["True abT", "True gdT"])
        ax.set_title(title)
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f"{cm[i, j]:,}", ha="center", va="center", fontsize=11)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_full_fp_comparison(comparison_df: pd.DataFrame) -> Path:
    path = FIGURE_DIR / "full_atlas_known_fp_comparison.png"
    keep = comparison_df.drop_duplicates("strategy").copy()
    keep = keep.sort_values("predicted_paired_TCRAB_or_NK_FP", ascending=False)
    fig, ax = plt.subplots(figsize=(8.8, max(4.0, 0.45 * len(keep) + 1.2)), constrained_layout=True)
    y = np.arange(len(keep))
    ax.barh(y - 0.12, keep["predicted_NK_FP"], height=0.24, color="#4f46e5", label="NK FP")
    ax.barh(y + 0.12, keep["predicted_paired_TCRAB_FP"], height=0.24, color="#dc2626", label="paired TCRAB FP")
    ax.set_yticks(y, labels=keep["strategy"])
    ax.invert_yaxis()
    ax.set_xlabel("Full-atlas known false-positive proxy count")
    ax.set_title("Known FP proxy comparison")
    ax.legend(frameon=False)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def bounded_limits(x: np.ndarray, y: np.ndarray) -> tuple[tuple[float, float], tuple[float, float]]:
    finite = np.isfinite(x) & np.isfinite(y)
    if not finite.any():
        return (-1.0, 1.0), (-1.0, 1.0)
    x_lo, x_hi = np.nanpercentile(x[finite], [0.2, 99.8])
    y_lo, y_hi = np.nanpercentile(y[finite], [0.2, 99.8])
    pad_x = max((x_hi - x_lo) * 0.08, 0.05)
    pad_y = max((y_hi - y_lo) * 0.08, 0.05)
    return (float(x_lo - pad_x), float(x_hi + pad_x)), (float(y_lo - pad_y), float(y_hi + pad_y))


def plot_trd_trab_scatters(
    handle: h5py.File,
    arrays: dict[str, np.ndarray],
    *,
    sample_cells: int,
    seed: int,
) -> list[Path]:
    logging.info("Rendering TRD-vs-TRAB scatter figures")
    n_obs = arrays["current_pred"].size
    rng = np.random.default_rng(seed)
    random_idx = rng.choice(n_obs, size=min(sample_cells, n_obs), replace=False)
    current_idx = np.flatnonzero(arrays["current_pred"])
    cascade_idx = np.flatnonzero(arrays["cascade_pred"])
    original_idx = np.flatnonzero(arrays["original_pred"])
    for name, idx in [("current", current_idx), ("cascade", cascade_idx), ("original", original_idx)]:
        if idx.size > sample_cells:
            idx = rng.choice(idx, size=sample_cells, replace=False)
        if name == "current":
            current_idx = idx
        elif name == "cascade":
            cascade_idx = idx
        else:
            original_idx = idx
    sample_idx = np.unique(np.concatenate([random_idx, current_idx, cascade_idx, original_idx])).astype(np.int64)
    trab = read_float_obs(handle, "phase4_trab_score")[sample_idx]
    trd = read_float_obs(handle, "phase4_trd_score")[sample_idx]
    source = clean_group_values(read_obs_column(handle, "source_gse_id"))[sample_idx]
    annotation = clean_group_values(read_obs_column(handle, "simple_annotation_plus6"))[sample_idx]
    current_pred = arrays["current_pred"][sample_idx]
    cascade_pred = arrays["cascade_pred"][sample_idx]
    original_pred = arrays["original_pred"][sample_idx]
    gate_status = np.asarray([GATE_CODE_TO_NAME[int(code)] for code in arrays["gate_code"][sample_idx]], dtype=object)
    agreement = np.full(sample_idx.size, "neither", dtype=object)
    agreement[current_pred & (~cascade_pred)] = "current only"
    agreement[cascade_pred & (~current_pred)] = "cascade only"
    agreement[current_pred & cascade_pred] = "both current and cascade"
    scatter_df = pd.DataFrame(
        {
            "obs_index": sample_idx,
            "source_gse_id": source,
            "annotation": annotation,
            "phase4_trab_score": trab,
            "phase4_trd_score": trd,
            "current_gdTAI_score": arrays["current_score"][sample_idx],
            "current_gdTAI_predicted": current_pred,
            "cascade_predicted": cascade_pred,
            "original_trd_minus_trab_predicted": original_pred,
            "current_vs_cascade": agreement,
            "gate_status": gate_status,
        }
    )
    scatter_df.to_csv(TABLE_DIR / "trd_vs_trab_current_vs_cascade_scatter_sample.csv.gz", index=False)

    paths: list[Path] = []
    xlim, ylim = bounded_limits(trab, trd)
    for col, order, colors, title, name in [
        (
            "current_vs_cascade",
            ["neither", "current only", "cascade only", "both current and cascade"],
            {
                "neither": "#cbd5e1",
                "current only": "#2563eb",
                "cascade only": "#0f766e",
                "both current and cascade": "#7c3aed",
            },
            "TRD vs TRAB: current gdTAI and transcriptome-gated cascade",
            "trd_vs_trab_current_vs_cascade.png",
        ),
        (
            "gate_status",
            [
                "below_threshold",
                "high_conf_pass",
                "cd4_strict_pass",
                "restricted_rescue_pass",
                "blocked_cd4_or_treg_penalty",
                "blocked_strong_nk_or_low_cd3",
                "blocked_by_transcriptome_gate",
            ],
            {
                "below_threshold": "#cbd5e1",
                "high_conf_pass": "#0f766e",
                "cd4_strict_pass": "#22c55e",
                "restricted_rescue_pass": "#f59e0b",
                "blocked_cd4_or_treg_penalty": "#dc2626",
                "blocked_strong_nk_or_low_cd3": "#4f46e5",
                "blocked_by_transcriptome_gate": "#64748b",
            },
            "TRD vs TRAB: transcriptome gate status",
            "trd_vs_trab_cascade_gate_status.png",
        ),
    ]:
        path = FIGURE_DIR / name
        fig, ax = plt.subplots(figsize=(7.4, 6.0), constrained_layout=True)
        for category in order:
            mask = scatter_df[col].astype(str).to_numpy() == category
            if not mask.any():
                continue
            alpha = 0.18 if category in {"neither", "below_threshold"} else 0.70
            size = 2.0 if category in {"neither", "below_threshold"} else 4.5
            ax.scatter(
                trab[mask],
                trd[mask],
                s=size,
                linewidths=0,
                alpha=alpha,
                rasterized=True,
                color=colors.get(category, "#64748b"),
                label=f"{category} (n={int(mask.sum()):,})",
            )
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_xlabel("Phase 4 TRAB score")
        ax.set_ylabel("Phase 4 TRD score")
        ax.set_title(title)
        ax.legend(frameon=False, fontsize=7, loc="best")
        fig.savefig(path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
    return paths


def plot_umap(handle: h5py.File, arrays: dict[str, np.ndarray], *, seed: int, background_cells: int, positive_cells: int) -> Path | None:
    logging.info("Rendering UMAP prediction figure with background=%s positive=%s", f"{background_cells:,}", f"{positive_cells:,}")
    if "obsm" not in handle or "X_umap" not in handle["obsm"]:
        return None
    path = FIGURE_DIR / "umap_transcriptome_gated_cascade_predictions.png"
    n_obs = arrays["cascade_pred"].size
    rng = np.random.default_rng(seed)
    bg_idx = rng.choice(n_obs, size=min(background_cells, n_obs), replace=False)
    pos_idx = np.flatnonzero(arrays["cascade_pred"])
    current_only_idx = np.flatnonzero(arrays["current_pred"] & (~arrays["cascade_pred"]))
    if pos_idx.size > positive_cells:
        pos_idx = rng.choice(pos_idx, size=positive_cells, replace=False)
    if current_only_idx.size > positive_cells // 2:
        current_only_idx = rng.choice(current_only_idx, size=positive_cells // 2, replace=False)
    umap = handle["obsm"]["X_umap"]
    bg = umap[np.sort(bg_idx), :]
    pos = umap[np.sort(pos_idx), :] if pos_idx.size else np.empty((0, 2))
    current_only = umap[np.sort(current_only_idx), :] if current_only_idx.size else np.empty((0, 2))
    fig, ax = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    ax.scatter(bg[:, 0], bg[:, 1], s=0.35, color="#d1d5db", alpha=0.18, linewidths=0, rasterized=True, label="atlas background")
    if current_only.size:
        ax.scatter(current_only[:, 0], current_only[:, 1], s=0.7, color="#2563eb", alpha=0.35, linewidths=0, rasterized=True, label="current only")
    if pos.size:
        ax.scatter(pos[:, 0], pos[:, 1], s=0.9, color="#0f766e", alpha=0.55, linewidths=0, rasterized=True, label="cascade gdT")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("5 million cell atlas: transcriptome-gated gdTAI")
    ax.legend(frameon=False, fontsize=8, loc="best")
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    logging.info("Saved UMAP prediction figure: %s", path)
    return path


def copy_assets(paths: list[Path | None]) -> list[str]:
    copied: list[str] = []
    for path in paths:
        if path is None or not path.exists():
            continue
        target = STATIC_ASSET_DIR / path.name
        shutil.copyfile(path, target)
        copied.append(path.name)
    return copied


def render_pdf(no_pdf: bool) -> None:
    if no_pdf:
        return
    chrome_candidates = [
        Path("/usr/bin/google-chrome"),
        Path("/usr/bin/google-chrome-stable"),
        shutil.which("google-chrome"),
        shutil.which("google-chrome-stable"),
    ]
    chrome = next((Path(path) for path in chrome_candidates if path and Path(path).exists()), None)
    if chrome is None:
        logging.warning("Skipping PDF export because google-chrome was not found")
        return
    subprocess.run(
        [
            str(chrome),
            "--headless",
            "--disable-gpu",
            "--no-sandbox",
            "--print-to-pdf-no-header",
            f"--print-to-pdf={REPORT_PDF}",
            str(REPORT_HTML),
        ],
        check=True,
    )


def write_protocol(model: LoadedModel, params: CascadeParams, accepted: bool, wrapper_path: Path | None) -> Path:
    path = MODEL_DIR / "gdtai_transcriptome_gate_protocol.md"
    lines = [
        "# gdTAI Transcriptome-Gated Cascade Protocol",
        "",
        "## Purpose",
        "",
        "This wrapper first uses the scVI-derived transcriptome annotation to reject or rescue cells, then applies the frozen individual TCR-gene gdTAI model. It is a post-processing cascade, not a newly trained classifier.",
        "",
        "## Required Inputs",
        "",
        "- Human gene symbols in `var/_index`.",
        "- Count-like sparse `X` so the wrapper can reconstruct `log1p(counts per 10,000)` features.",
        "- `obs['simple_annotation_plus6']` or an equivalent simple transcriptome annotation mapped to `CD8_T`, `CD4_T`, `Treg`, `NK_cell`, `gdT_cell`, and `other`.",
        "",
        "## Frozen Model",
        "",
        f"- Base model: `{model.model_name}`",
        f"- Base model threshold: `{model.threshold}`",
        f"- Selected TCR-gene features: `{len(model.gene_names)}`",
        "",
        "## Cascade Parameters",
        "",
        dataframe_to_markdown(pd.DataFrame([asdict(params)])),
        "",
        "## Decision Policy",
        "",
        "- `gdT_cell` and `CD8_T` annotations pass with the high-confidence threshold.",
        "- `CD4_T` cells require the stricter CD4 threshold and low FOXP3.",
        "- `Treg`, `NK_cell`, `other`, and unknown annotations are rescue-only: they require the rescue threshold, sufficient CD3/T-lineage marker support, low FOXP3/CD4 rescue penalty, and no strong NK-marker veto.",
        "- Cells blocked by the gate should be interpreted as low-confidence/no-call rather than re-labeled as biologically impossible gdT cells.",
        "",
        "## Status",
        "",
        f"- Promotion accepted: `{accepted}`",
        f"- Wrapper artifact: `{wrapper_path}`" if wrapper_path is not None else "- Wrapper artifact: not written because promotion was not accepted.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_reports(
    input_h5ad: Path,
    model: LoadedModel,
    params: CascadeParams,
    tune_df: pd.DataFrame,
    validation_df: pd.DataFrame,
    full_overall: pd.DataFrame,
    full_fp_overall: pd.DataFrame,
    full_source: pd.DataFrame,
    full_annotation: pd.DataFrame,
    comparison_df: pd.DataFrame,
    acceptance_df: pd.DataFrame,
    split_overall: pd.DataFrame,
    copied_assets: list[str],
    accepted: bool,
    wrapper_path: Path | None,
) -> None:
    tune_display = compact_numeric(
        pick_columns(
            tune_df,
            ["strategy", "threshold", "n_cells", "n_positive", "n_negative", "predicted_positive", "tp", "fp", "tn", "fn", "precision", "recall", "specificity", "f1", "roc_auc", "pr_auc"],
        )
    )
    valid_display = compact_numeric(
        pick_columns(
            validation_df,
            ["strategy", "threshold", "n_cells", "n_positive", "n_negative", "predicted_positive", "tp", "fp", "tn", "fn", "precision", "recall", "specificity", "f1", "roc_auc", "pr_auc"],
        )
    )
    fp_display = compact_numeric(
        pick_columns(
            comparison_df,
            [
                "strategy",
                "threshold",
                "predicted_putative_gdT",
                "predicted_fraction",
                "predicted_paired_TCRAB_FP",
                "paired_TCRAB_FP_fraction_of_predictions",
                "predicted_NK_FP",
                "NK_FP_fraction_of_predictions",
                "predicted_paired_TCRAB_or_NK_FP",
                "known_FP_fraction_of_predictions",
            ],
        )
    )
    source_subset = (
        full_source.loc[full_source["strategy"] == "transcriptome_gated_cascade"].head(30)
        if "strategy" in full_source.columns
        else pd.DataFrame()
    )
    source_display = compact_numeric(
        pick_columns(
            source_subset,
            ["source_gse_id", "total_cells", "predicted_putative_gdT", "predicted_fraction", "predicted_paired_TCRAB_FP", "predicted_NK_FP", "known_FP_fraction_of_predictions"],
        )
    )
    annotation_subset = (
        full_annotation.loc[full_annotation["strategy"] == "transcriptome_gated_cascade"]
        if "strategy" in full_annotation.columns
        else pd.DataFrame()
    )
    annotation_display = compact_numeric(
        pick_columns(
            annotation_subset,
            ["annotation", "total_cells", "predicted_putative_gdT", "predicted_fraction", "predicted_paired_TCRAB_FP", "predicted_NK_FP", "known_FP_fraction_of_predictions"],
        )
    )
    acceptance_display = compact_numeric(acceptance_df)
    param_display = compact_numeric(pd.DataFrame([asdict(params)]))

    current_fp = full_fp_overall.loc[full_fp_overall["strategy"] == "current_gdTAI"].iloc[0] if not full_fp_overall.empty else None
    cascade_fp = full_fp_overall.loc[full_fp_overall["strategy"] == "transcriptome_gated_cascade"].iloc[0] if not full_fp_overall.empty else None
    decision = "accepted and written as a wrapper artifact" if accepted else "not promoted; kept as a candidate report"
    if current_fp is not None and cascade_fp is not None:
        decision_detail = (
            f"The cascade predicts {int(cascade_fp['predicted_putative_gdT']):,} putative gdT cells, "
            f"with {int(cascade_fp['predicted_NK_FP']):,} NK false-positive proxies and "
            f"{int(cascade_fp['predicted_paired_TCRAB_FP']):,} paired-TCRAB false-positive proxies. "
            f"Current gdTAI had {int(current_fp['predicted_putative_gdT']):,} predictions, "
            f"{int(current_fp['predicted_NK_FP']):,} NK FP proxies, and "
            f"{int(current_fp['predicted_paired_TCRAB_FP']):,} paired-TCRAB FP proxies."
        )
    else:
        decision_detail = "Full-atlas application was skipped."

    md_lines = [
        "# gdTAI Transcriptome-Gated Cascade Report",
        "",
        "## Executive Summary",
        "",
        f"- Input atlas: `{input_h5ad}`",
        "- This is a read-only wrapper around the frozen individual TCR-gene gdTAI model.",
        "- Phase 4 TRD/TRAB scores are used only for comparison figures, not for cascade prediction decisions.",
        f"- Promotion decision: `{decision}`.",
        f"- {decision_detail}",
        f"- Wrapper artifact: `{wrapper_path}`" if wrapper_path is not None else "- Wrapper artifact: not written.",
        "",
        "## Why This Iteration",
        "",
        "Most current known false positives are transcriptome-inconsistent cells, especially NK-annotated cells. A hard T-cell-only filter would reduce NK false positives but could discard true gdT cells with NK-like or low-quality transcriptomes. This cascade keeps a rescue path and treats contradictory cells as low-confidence/no-call.",
        "",
        "## Training And Test Protocol",
        "",
        "- No model weights were retrained.",
        "- The frozen gdTAI model uses individual TCR A/B/G/D genes plus the existing penalty-control genes from the selected model.",
        "- Wrapper thresholds were tuned only on the existing non-holdout tune split.",
        f"- Strict validation holds out all `{HOLDOUT_SOURCE}` primary gold cells, paired-TCRAB/no-gdTCR cells from `{EXTRA_TRAB_HOLDOUT_SOURCE}`, and `{GDT2020_SOURCE}` `{GDT2020_HOLDOUT_TISSUE}` gdT cells.",
        f"- `{SUBOPTIMAL_SOURCE}` remains excluded from train/tune sensitivity because its library quality was flagged as suboptimal.",
        "- The cascade gate uses the scVI-derived `simple_annotation_plus6` layer plus small marker panels for CD3/T-lineage support, NK risk, and FOXP3/CD4 death-penalty risk.",
        "",
        "## Split Summary",
        "",
        dataframe_to_markdown(compact_numeric(split_overall)),
        "",
        "## Cascade Parameters",
        "",
        dataframe_to_markdown(param_display),
        "",
        "## Tune Metrics",
        "",
        dataframe_to_markdown(tune_display),
        "",
        "## Strict Validation Metrics",
        "",
        dataframe_to_markdown(valid_display),
        "",
        "## Acceptance Gates",
        "",
        dataframe_to_markdown(acceptance_display),
        "",
        "## Full Atlas Known False-Positive Proxy Comparison",
        "",
        dataframe_to_markdown(fp_display),
        "",
        "## Cascade Predictions By Annotation",
        "",
        dataframe_to_markdown(annotation_display),
        "",
        "## Top Cascade Prediction Sources",
        "",
        dataframe_to_markdown(source_display),
        "",
        "## Output Files",
        "",
        f"- Tables: `{TABLE_DIR}`",
        f"- Figures: `{FIGURE_DIR}`",
        f"- HTML: `{REPORT_HTML}`",
        f"- PDF: `{REPORT_PDF}`",
        f"- Protocol: `{MODEL_DIR / 'gdtai_transcriptome_gate_protocol.md'}`",
    ]
    REPORT_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    figures_html = []
    for asset in copied_assets:
        title = asset.rsplit(".", 1)[0].replace("_", " ")
        figures_html.append(f"<section class='figure'><h3>{html.escape(title)}</h3><img src='assets/{OUT_PREFIX}/{html.escape(asset)}'></section>")
    css = """
    body{font-family:Arial,Helvetica,sans-serif;color:#111827;background:#f8fafc;margin:0}
    main{max-width:1180px;margin:0 auto;padding:28px 28px 64px;background:#fff}
    h1{font-size:30px;margin:0 0 8px} h2{font-size:21px;margin-top:28px;border-top:1px solid #e5e7eb;padding-top:18px}
    h3{font-size:15px;margin:14px 0 8px}.status{font-weight:700;color:#0f766e}
    table{border-collapse:collapse;width:100%;font-size:12px;margin:10px 0 18px;table-layout:auto}
    th,td{border:1px solid #e5e7eb;padding:5px 7px;text-align:left;vertical-align:top}
    th{background:#f1f5f9}.table-wrap{overflow-x:auto;margin-bottom:12px}
    img{max-width:100%;height:auto;border:1px solid #e5e7eb;background:#fff}
    code{background:#f1f5f9;padding:1px 4px;border-radius:4px}.note{color:#475569}
    @media print{main{max-width:none;padding:12mm} table{font-size:9px} h1{font-size:22px} h2{font-size:16px} .table-wrap{overflow:visible}}
    """
    html_doc = f"""<!doctype html><html><head><meta charset='utf-8'><title>gdTAI Transcriptome-Gated Cascade</title><style>{css}</style></head><body><main>
    <h1>gdTAI Transcriptome-Gated Cascade</h1>
    <p class='status'>Promotion decision: {html.escape(decision)}</p>
    <p>{html.escape(decision_detail)}</p>
    <h2>Protocol</h2>
    <p>This report evaluates a two-step transcriptome-aware wrapper around the frozen gdTAI TCR-gene model. It does not retrain the model and does not modify any H5AD.</p>
    <p>Prediction order: transcriptome annotation gate, marker-based rescue/veto checks, then the frozen individual TCR-gene model threshold. TRD/TRAB scores are shown only for interpretation.</p>
    <h2>Split Summary</h2><div class='table-wrap'>{dataframe_to_html(compact_numeric(split_overall))}</div>
    <h2>Cascade Parameters</h2><div class='table-wrap'>{dataframe_to_html(param_display)}</div>
    <h2>Tune Metrics</h2><div class='table-wrap'>{dataframe_to_html(tune_display)}</div>
    <h2>Strict Validation Metrics</h2><div class='table-wrap'>{dataframe_to_html(valid_display)}</div>
    <h2>Acceptance Gates</h2><div class='table-wrap'>{dataframe_to_html(acceptance_display)}</div>
    <h2>Full Atlas Known FP Proxy Comparison</h2><div class='table-wrap'>{dataframe_to_html(fp_display)}</div>
    <h2>Cascade Predictions By Annotation</h2><div class='table-wrap'>{dataframe_to_html(annotation_display)}</div>
    <h2>Top Cascade Prediction Sources</h2><div class='table-wrap'>{dataframe_to_html(source_display)}</div>
    <h2>Figures</h2>{''.join(figures_html)}
    <h2>Artifacts</h2>
    <p>Tables: <code>{html.escape(str(TABLE_DIR))}</code><br>Figures: <code>{html.escape(str(FIGURE_DIR))}</code><br>Markdown: <code>{html.escape(str(REPORT_MD))}</code><br>Wrapper: <code>{html.escape(str(wrapper_path)) if wrapper_path is not None else 'not written'}</code></p>
    </main></body></html>"""
    REPORT_HTML.write_text(html_doc, encoding="utf-8")


def write_summary_json(
    input_h5ad: Path,
    model_path: Path,
    model: LoadedModel,
    params: CascadeParams,
    accepted: bool,
    wrapper_path: Path | None,
    validation_df: pd.DataFrame,
    full_fp_overall: pd.DataFrame,
) -> None:
    summary = {
        "input_h5ad": str(input_h5ad),
        "current_model_path": str(model_path),
        "current_model": model.model_name,
        "current_threshold": model.threshold,
        "cascade_params": asdict(params),
        "promotion_accepted": accepted,
        "wrapper_path": str(wrapper_path) if wrapper_path else None,
        "validation_metrics": validation_df.to_dict(orient="records"),
        "full_false_positive_overall": full_fp_overall.to_dict(orient="records") if not full_fp_overall.empty else [],
        "report_html": str(REPORT_HTML),
        "report_pdf": str(REPORT_PDF),
    }
    SUMMARY_JSON.write_text(json.dumps(json_ready(summary), indent=2), encoding="utf-8")


def main() -> None:
    args = parse_args()
    ensure_dirs()
    setup_logging()
    input_h5ad = args.input_h5ad.resolve()
    model_path = args.current_model_pkl.resolve()
    before_mtime_ns = input_h5ad.stat().st_mtime_ns
    logging.info("Input H5AD: %s", input_h5ad)
    logging.info("Current model: %s", model_path)

    model = load_current_model(model_path)
    logging.info("Loaded frozen model %s with threshold %.6f and %s genes", model.model_name, model.threshold, len(model.gene_names))
    with h5py.File(input_h5ad, "r") as handle:
        obs = build_obs_metadata(handle)
        annotation = clean_group_values(read_obs_column(handle, "simple_annotation_plus6"))
        splits = build_splits(obs, args.seed)
        bundle = build_feature_bundle(handle, model)

        tune_y = (obs.class_code[splits.tune_idx] == 2).astype(np.int8)
        valid_y = (obs.class_code[splits.validation_idx] == 2).astype(np.int8)
        tune_block = score_rows(handle, splits.tune_idx, annotation, tune_y, model, bundle, label="cascade_tune")
        validation_block = score_rows(handle, splits.validation_idx, annotation, valid_y, model, bundle, label="cascade_validation")
        params, tuning_grid, _current_tune_df = tune_cascade(tune_block, model)
        logging.info("Selected cascade params: %s", params)
        pd.DataFrame([asdict(params)]).to_csv(TABLE_DIR / "selected_cascade_parameters.csv", index=False)
        tune_df, validation_df = evaluate_blocks(tune_block, validation_block, params, model, handle)

        full_overall = pd.DataFrame()
        full_source = pd.DataFrame()
        full_fp_overall = pd.DataFrame()
        full_fp_source = pd.DataFrame()
        full_annotation = pd.DataFrame()
        arrays: dict[str, np.ndarray] = {}
        if not args.skip_full_apply:
            full_overall, full_source, full_tissue, full_fp_overall, full_fp_source, full_annotation, arrays = apply_full_atlas(
                handle,
                model,
                bundle,
                params,
                chunk_size=args.chunk_size,
            )
        else:
            logging.warning("Skipping full-atlas application by request")

        if full_fp_overall.empty:
            accepted = False
            acceptance_df = pd.DataFrame()
            comparison_df = pd.DataFrame()
            wrapper_path = None
            copied_assets = copy_assets([plot_validation_metric_bars(validation_df), plot_confusion_pair(validation_block, params, model)])
        else:
            comparison_df = add_guard_comparisons(full_fp_overall)
            accepted, acceptance_df = evaluate_acceptance(validation_df, full_fp_overall)
            wrapper_path = save_wrapper_if_accepted(accepted, model_path, model, params, bundle)
            figure_paths = [
                plot_validation_metric_bars(validation_df),
                plot_confusion_pair(validation_block, params, model),
                plot_full_fp_comparison(comparison_df),
            ]
            if arrays:
                figure_paths.extend(plot_trd_trab_scatters(handle, arrays, sample_cells=args.scatter_sample_cells, seed=args.seed))
                figure_paths.append(
                    plot_umap(
                        handle,
                        arrays,
                        seed=args.seed,
                        background_cells=args.umap_background_cells,
                        positive_cells=args.umap_positive_cells,
                    )
                )
            copied_assets = copy_assets(figure_paths)

        write_protocol(model, params, accepted, wrapper_path)
        write_reports(
            input_h5ad=input_h5ad,
            model=model,
            params=params,
            tune_df=tune_df,
            validation_df=validation_df,
            full_overall=full_overall,
            full_fp_overall=full_fp_overall,
            full_source=full_source,
            full_annotation=full_annotation,
            comparison_df=comparison_df,
            acceptance_df=acceptance_df,
            split_overall=splits.split_overall,
            copied_assets=copied_assets,
            accepted=accepted,
            wrapper_path=wrapper_path,
        )
        write_summary_json(input_h5ad, model_path, model, params, accepted, wrapper_path, validation_df, full_fp_overall)

    render_pdf(args.no_pdf)
    after_mtime_ns = input_h5ad.stat().st_mtime_ns
    if after_mtime_ns != before_mtime_ns:
        raise RuntimeError("Input H5AD changed during read-only transcriptome-gated cascade analysis.")
    logging.info("Saved report: %s", REPORT_HTML)
    logging.info("Saved summary: %s", SUMMARY_JSON)


if __name__ == "__main__":
    main()
