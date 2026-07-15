#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Annotation-specific threshold cascade around the frozen gdTAI model.

This read-only candidate keeps the selected individual TCR-gene gdTAI model
frozen, then tunes separate thresholds for scVI-derived simple annotations. It
is designed to preserve CD8/gdT recovery while applying strict thresholds only
to high-risk transcriptome compartments such as NK, Treg, and other.
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

import run_gdtai_transcriptome_gate_cascade as base
from run_gdt_prediction_package_evaluation import (
    DEFAULT_INPUT_H5AD,
    clean_group_values,
    dataframe_to_html,
    dataframe_to_markdown,
    read_float_obs,
    read_obs_column,
)
from run_gdt_deg_tcr_classifier_training import json_ready, metric_dict, safe_div


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
STATIC_DIR = PROJECT_ROOT / "gdT_prediction"
OUT_PREFIX = "gdtai_annotation_specific_cascade"
TABLE_DIR = OUTPUT_ROOT / "tables" / "gdT_prediction" / OUT_PREFIX
FIGURE_DIR = OUTPUT_ROOT / "figures" / "gdT_prediction" / OUT_PREFIX
LOG_DIR = OUTPUT_ROOT / "logs" / "gdT_prediction" / OUT_PREFIX
MODEL_DIR = OUTPUT_ROOT / "models" / "gdT_prediction_classifier" / OUT_PREFIX
STATIC_ASSET_DIR = STATIC_DIR / "assets" / OUT_PREFIX
REPORT_MD = LOG_DIR / "gdtai_annotation_specific_cascade_report.md"
REPORT_HTML = STATIC_DIR / "gdtai_annotation_specific_cascade_report.html"
REPORT_PDF = STATIC_DIR / "gdtai_annotation_specific_cascade_report.pdf"
SUMMARY_JSON = LOG_DIR / "gdtai_annotation_specific_cascade_summary.json"
RUN_LOG = LOG_DIR / "run.log"
WRAPPER_PKL = MODEL_DIR / "selected_annotation_specific_wrapper.pkl"
PROTOCOL_MD = MODEL_DIR / "gdtai_annotation_specific_cascade_protocol.md"

DEFAULT_SCORE_CACHE = OUTPUT_ROOT / "tables" / "gdT_prediction" / "gdtai_transcriptome_gate" / "full_prediction_arrays.npz"
CURRENT_MODEL_PATH = base.CURRENT_MODEL_PATH
PREVIOUS_GLOBAL_TABLE_DIR = OUTPUT_ROOT / "tables" / "gdT_prediction" / "gdtai_transcriptome_gate"
NKGUARD_TABLE_DIR = OUTPUT_ROOT / "tables" / "gdT_prediction" / "gdtai_nkguard"
MULTIGUARD_TABLE_DIR = OUTPUT_ROOT / "tables" / "gdT_prediction" / "gdtai_multiguard"

FULL_APPLY_CHUNK = 50_000
SCATTER_SAMPLE_CELLS = 250_000
UMAP_BACKGROUND_CELLS = 50_000
UMAP_POSITIVE_CELLS = 20_000
RANDOM_SEED = base.RANDOM_SEED

ACCEPT_RECALL_MIN = 0.86
ACCEPT_VALIDATION_DELTA_RECALL_MIN = -0.01
ACCEPT_VALIDATION_DELTA_F1_MIN = -0.005
ACCEPT_FULL_NK_FP_REDUCTION_MIN = 0.50
ACCEPT_FULL_TCRAB_FP_MAX_VS_CURRENT = 1.00
ACCEPT_FULL_KNOWN_FP_FRACTION_DELTA_MAX = -1e-12

ANNOTATION_ORDER = ["gdT_cell", "CD8_T", "CD4_T", "Treg", "NK_cell", "other"]
ANNOTATION_KEYS = ["GDT_CELL", "CD8_T", "CD4_T", "TREG", "NK_CELL"]


@dataclass(frozen=True)
class AnnotationThresholds:
    gdt_threshold: float
    cd8_threshold: float
    cd4_threshold: float
    treg_threshold: float
    nk_threshold: float
    other_threshold: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate annotation-specific gdTAI cascade.")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--current-model-pkl", type=Path, default=CURRENT_MODEL_PATH)
    parser.add_argument("--score-cache", type=Path, default=DEFAULT_SCORE_CACHE)
    parser.add_argument("--recompute-score-cache", action="store_true")
    parser.add_argument("--chunk-size", type=int, default=FULL_APPLY_CHUNK)
    parser.add_argument("--scatter-sample-cells", type=int, default=SCATTER_SAMPLE_CELLS)
    parser.add_argument("--umap-background-cells", type=int, default=UMAP_BACKGROUND_CELLS)
    parser.add_argument("--umap-positive-cells", type=int, default=UMAP_POSITIVE_CELLS)
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    parser.add_argument("--no-pdf", action="store_true")
    return parser.parse_args()


def configure_base_paths() -> None:
    base.OUT_PREFIX = OUT_PREFIX
    base.TABLE_DIR = TABLE_DIR
    base.FIGURE_DIR = FIGURE_DIR
    base.LOG_DIR = LOG_DIR
    base.MODEL_DIR = MODEL_DIR
    base.STATIC_ASSET_DIR = STATIC_ASSET_DIR
    base.REPORT_MD = REPORT_MD
    base.REPORT_HTML = REPORT_HTML
    base.REPORT_PDF = REPORT_PDF
    base.SUMMARY_JSON = SUMMARY_JSON


def ensure_dirs() -> None:
    configure_base_paths()
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


def normalize_annotation(annotation: np.ndarray) -> np.ndarray:
    return pd.Series(annotation, copy=False).astype(str).str.upper().str.strip().to_numpy(dtype=object)


def predict_annotation_specific(score: np.ndarray, annotation: np.ndarray, params: AnnotationThresholds) -> np.ndarray:
    ann = normalize_annotation(annotation)
    pred = np.zeros(score.size, dtype=bool)
    pred[(ann == "GDT_CELL") & (score >= params.gdt_threshold)] = True
    pred[(ann == "CD8_T") & (score >= params.cd8_threshold)] = True
    pred[(ann == "CD4_T") & (score >= params.cd4_threshold)] = True
    pred[(ann == "TREG") & (score >= params.treg_threshold)] = True
    pred[(ann == "NK_CELL") & (score >= params.nk_threshold)] = True
    other = ~np.isin(ann, ANNOTATION_KEYS)
    pred[other & (score >= params.other_threshold)] = True
    return pred


def gate_status(score: np.ndarray, annotation: np.ndarray, pred: np.ndarray, current_threshold: float) -> np.ndarray:
    ann = normalize_annotation(annotation)
    out = np.full(score.size, "below_current_threshold", dtype=object)
    current_positive = score >= current_threshold
    out[current_positive & (~pred)] = "blocked_or_stricter_threshold"
    out[pred & (ann == "GDT_CELL")] = "gdT_annotation_pass"
    out[pred & (ann == "CD8_T")] = "CD8_T_pass"
    out[pred & (ann == "CD4_T")] = "CD4_T_pass"
    out[pred & (ann == "TREG")] = "Treg_rescue_pass"
    out[pred & (ann == "NK_CELL")] = "NK_rescue_pass"
    out[pred & (~np.isin(ann, ANNOTATION_KEYS))] = "other_rescue_pass"
    return out


def metrics_row(strategy: str, y_true: np.ndarray, score: np.ndarray, pred: np.ndarray, threshold: float | str) -> dict[str, Any]:
    row = {"strategy": strategy, "threshold": threshold}
    row.update(metric_dict(y_true, pred.astype(np.int8), score))
    return row


def threshold_value_for_csv(value: float) -> str | float:
    return "disabled" if np.isinf(value) else float(value)


def params_for_csv(params: AnnotationThresholds) -> dict[str, Any]:
    return {key: threshold_value_for_csv(value) for key, value in asdict(params).items()}


def threshold_label(params: AnnotationThresholds) -> str:
    return "; ".join(f"{key.replace('_threshold','')}={threshold_value_for_csv(value)}" for key, value in asdict(params).items())


def load_or_compute_scores(
    handle: h5py.File,
    model: base.LoadedModel,
    score_cache: Path,
    *,
    recompute: bool,
    chunk_size: int,
) -> dict[str, np.ndarray]:
    n_obs = int(handle["obs"]["_index"].shape[0])
    if score_cache.exists() and not recompute:
        logging.info("Loading existing full-score cache: %s", score_cache)
        with np.load(score_cache) as arrays:
            score = np.asarray(arrays["current_score"], dtype=np.float32)
            original_pred = np.asarray(arrays["original_pred"], dtype=bool) if "original_pred" in arrays else None
        if score.shape[0] != n_obs:
            raise RuntimeError(f"Score cache length {score.shape[0]} does not match H5AD n_obs {n_obs}")
        if original_pred is None:
            original_pred = read_float_obs(handle, "phase4_trd_minus_trab") >= base.load_original_trd_trab_threshold()
        current_pred = score >= model.threshold
        return {"current_score": score, "current_pred": current_pred, "original_pred": original_pred}

    logging.info("Computing full current gdTAI scores because cache is missing or recompute was requested")
    bundle = base.build_feature_bundle(handle, model)
    score = np.zeros(n_obs, dtype=np.float32)
    for start in range(0, n_obs, chunk_size):
        end = min(start + chunk_size, n_obs)
        rows = np.arange(start, end, dtype=np.int64)
        x = base.extract_gene_matrix(handle, rows, bundle.spec, label=f"{OUT_PREFIX}_score_{start}_{end}")
        score[start:end] = base.predict_model_score(model, x[:, bundle.model_cols])
        if start and start % 500_000 == 0:
            logging.info("Computed frozen gdTAI score for %s / %s cells", f"{start:,}", f"{n_obs:,}")
    current_pred = score >= model.threshold
    original_pred = read_float_obs(handle, "phase4_trd_minus_trab") >= base.load_original_trd_trab_threshold()
    np.savez(TABLE_DIR / "full_current_score_cache.npz", current_score=score, current_pred=current_pred, original_pred=original_pred)
    return {"current_score": score, "current_pred": current_pred, "original_pred": original_pred}


def tune_thresholds(score: np.ndarray, annotation: np.ndarray, y_true: np.ndarray, current_threshold: float) -> tuple[AnnotationThresholds, pd.DataFrame, pd.DataFrame]:
    current_pred = score >= current_threshold
    current_metrics = metrics_row("current_gdTAI", y_true, score, current_pred, current_threshold)
    recall_floor = max(ACCEPT_RECALL_MIN, float(current_metrics["recall"]) - 0.01)
    f1_floor = float(current_metrics["f1"]) - 0.005
    gdt_candidates = [current_threshold]
    cd8_candidates = [current_threshold, 0.93, 0.95]
    cd4_candidates = [current_threshold, 0.93, 0.95, 0.97]
    restricted_candidates = [0.97, 0.985, 0.995, np.inf]
    rows: list[dict[str, Any]] = []
    for gdt in gdt_candidates:
        for cd8 in cd8_candidates:
            for cd4 in cd4_candidates:
                for treg in restricted_candidates:
                    for nk in restricted_candidates:
                        for other in restricted_candidates:
                            params = AnnotationThresholds(
                                gdt_threshold=float(gdt),
                                cd8_threshold=float(cd8),
                                cd4_threshold=float(cd4),
                                treg_threshold=float(treg),
                                nk_threshold=float(nk),
                                other_threshold=float(other),
                            )
                            pred = predict_annotation_specific(score, annotation, params)
                            metrics = metrics_row("annotation_specific_tune_candidate", y_true, score, pred, threshold_label(params))
                            row = params_for_csv(params)
                            row.update({f"tune_{key}": value for key, value in metrics.items() if key not in {"strategy", "threshold"}})
                            row["current_tune_precision"] = current_metrics["precision"]
                            row["current_tune_recall"] = current_metrics["recall"]
                            row["current_tune_f1"] = current_metrics["f1"]
                            row["passes_recall_floor"] = metrics["recall"] >= recall_floor
                            row["passes_f1_floor"] = metrics["f1"] >= f1_floor
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
        reason = "fallback_best_tune_f1"
    chosen = ranked.iloc[0]

    def parse_threshold(value: Any) -> float:
        return float("inf") if str(value) == "disabled" else float(value)

    params = AnnotationThresholds(
        gdt_threshold=parse_threshold(chosen["gdt_threshold"]),
        cd8_threshold=parse_threshold(chosen["cd8_threshold"]),
        cd4_threshold=parse_threshold(chosen["cd4_threshold"]),
        treg_threshold=parse_threshold(chosen["treg_threshold"]),
        nk_threshold=parse_threshold(chosen["nk_threshold"]),
        other_threshold=parse_threshold(chosen["other_threshold"]),
    )
    grid.insert(0, "selection_reason", "")
    grid.loc[grid.index == chosen.name, "selection_reason"] = reason
    grid.sort_values(
        ["passes_recall_floor", "passes_f1_floor", "tune_fp", "tune_precision", "tune_f1"],
        ascending=[False, False, True, False, False],
    ).to_csv(TABLE_DIR / "annotation_threshold_tuning_grid.csv", index=False)
    current_df = pd.DataFrame([current_metrics])
    current_df.to_csv(TABLE_DIR / "current_model_tune_metrics.csv", index=False)
    return params, grid, current_df


def evaluate_splits(
    score: np.ndarray,
    annotation: np.ndarray,
    obs: Any,
    splits: base.SplitBundle,
    params: AnnotationThresholds,
    current_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows_out: list[dict[str, Any]] = []
    gate_rows: list[pd.DataFrame] = []
    original_threshold = base.load_original_trd_trab_threshold()
    original_score_all = None
    for split_name, idx in [("tune", splits.tune_idx), ("validation", splits.validation_idx)]:
        y = (obs.class_code[idx] == 2).astype(np.int8)
        split_score = score[idx]
        split_ann = annotation[idx]
        current_pred = split_score >= current_threshold
        candidate_pred = predict_annotation_specific(split_score, split_ann, params)
        rows_out.append(metrics_row(f"{split_name}:current_gdTAI", y, split_score, current_pred, current_threshold))
        rows_out.append(metrics_row(f"{split_name}:annotation_specific_cascade", y, split_score, candidate_pred, threshold_label(params)))
        status = gate_status(split_score, split_ann, candidate_pred, current_threshold)
        gate_rows.append(pd.DataFrame({"split": split_name, "label": np.where(y == 1, "gdT_gold", "abT_gold"), "annotation": split_ann, "gate_status": status}))
    metrics_df = pd.DataFrame(rows_out)
    tune_df = metrics_df.loc[metrics_df["strategy"].str.startswith("tune:")].copy()
    validation_df = metrics_df.loc[metrics_df["strategy"].str.startswith("validation:")].copy()
    tune_df["strategy"] = tune_df["strategy"].str.replace("tune:", "", regex=False)
    validation_df["strategy"] = validation_df["strategy"].str.replace("validation:", "", regex=False)
    tune_df.to_csv(TABLE_DIR / "tune_metrics.csv", index=False)
    validation_df.to_csv(TABLE_DIR / "validation_metrics.csv", index=False)
    gate_df = pd.concat(gate_rows, ignore_index=True)
    gate_summary = gate_df.groupby(["split", "label", "annotation", "gate_status"], dropna=False).size().reset_index(name="n_cells")
    gate_summary.to_csv(TABLE_DIR / "gate_status_by_split_label_annotation.csv", index=False)
    return tune_df, validation_df


def summarize_full(
    handle: h5py.File,
    arrays: dict[str, np.ndarray],
    annotation: np.ndarray,
    params: AnnotationThresholds,
    current_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, np.ndarray]]:
    meta = base.build_full_metadata(handle)
    score = arrays["current_score"]
    current_pred = score >= current_threshold
    candidate_pred = predict_annotation_specific(score, annotation, params)
    original_threshold = base.load_original_trd_trab_threshold()
    original_pred = arrays.get("original_pred")
    if original_pred is None:
        original_pred = read_float_obs(handle, "phase4_trd_minus_trab") >= original_threshold
    summaries = []
    for strategy, threshold, pred in [
        ("current_gdTAI", float(current_threshold), current_pred),
        ("annotation_specific_cascade", threshold_label(params), candidate_pred),
        ("original_TRD_minus_TRAB_score_strategy", float(original_threshold), original_pred),
    ]:
        summaries.append(
            base.summarize_full_predictions(
                strategy=strategy,
                threshold=threshold,
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
    fp_overall = overall[[
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
    ]].copy()
    status = gate_status(score, annotation, candidate_pred, current_threshold)
    pd.DataFrame({
        "annotation": annotation,
        "gate_status": status,
        "current_gdTAI_predicted": current_pred,
        "annotation_specific_predicted": candidate_pred,
    }).groupby(["annotation", "gate_status", "current_gdTAI_predicted", "annotation_specific_predicted"], dropna=False).size().reset_index(name="n_cells").to_csv(TABLE_DIR / "full_atlas_gate_status_by_annotation.csv", index=False)
    overall.to_csv(TABLE_DIR / "full_dataset_prediction_overall.csv", index=False)
    source_df.to_csv(TABLE_DIR / "full_dataset_prediction_by_source.csv", index=False)
    tissue_df.to_csv(TABLE_DIR / "full_dataset_prediction_by_tissue.csv", index=False)
    fp_overall.to_csv(TABLE_DIR / "full_dataset_false_positive_overall.csv", index=False)
    fp_source_df.to_csv(TABLE_DIR / "full_dataset_false_positive_by_source.csv", index=False)
    sublibrary_df.to_csv(TABLE_DIR / "full_dataset_false_positive_by_sublibrary.csv", index=False)
    annotation_df.to_csv(TABLE_DIR / "full_dataset_prediction_by_annotation.csv", index=False)
    out_arrays = {
        "current_score": score,
        "current_pred": current_pred,
        "annotation_specific_pred": candidate_pred,
        "original_pred": original_pred,
    }
    np.savez(TABLE_DIR / "full_prediction_arrays.npz", **out_arrays)
    return overall, source_df, tissue_df, fp_overall, fp_source_df, annotation_df, out_arrays


def add_comparison_tables(fp_overall: pd.DataFrame) -> pd.DataFrame:
    frames = [fp_overall.copy()]
    prev = PREVIOUS_GLOBAL_TABLE_DIR / "full_dataset_false_positive_overall.csv"
    if prev.exists():
        df = pd.read_csv(prev)
        df = df.loc[df["strategy"].astype(str).eq("transcriptome_gated_cascade")].copy()
        if not df.empty:
            df["strategy"] = "previous_global_transcriptome_gate"
            frames.append(df)
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


def evaluate_acceptance(validation_df: pd.DataFrame, fp_overall: pd.DataFrame) -> tuple[bool, pd.DataFrame]:
    current_valid = validation_df.loc[validation_df["strategy"] == "current_gdTAI"].iloc[0]
    candidate_valid = validation_df.loc[validation_df["strategy"] == "annotation_specific_cascade"].iloc[0]
    current_full = fp_overall.loc[fp_overall["strategy"] == "current_gdTAI"].iloc[0]
    candidate_full = fp_overall.loc[fp_overall["strategy"] == "annotation_specific_cascade"].iloc[0]
    checks = [
        {
            "gate": "validation_recall_min",
            "observed": float(candidate_valid["recall"]),
            "required": ACCEPT_RECALL_MIN,
            "passed": float(candidate_valid["recall"]) >= ACCEPT_RECALL_MIN,
        },
        {
            "gate": "validation_recall_delta_vs_current",
            "observed": float(candidate_valid["recall"] - current_valid["recall"]),
            "required": ACCEPT_VALIDATION_DELTA_RECALL_MIN,
            "passed": float(candidate_valid["recall"] - current_valid["recall"]) >= ACCEPT_VALIDATION_DELTA_RECALL_MIN,
        },
        {
            "gate": "validation_f1_delta_vs_current",
            "observed": float(candidate_valid["f1"] - current_valid["f1"]),
            "required": ACCEPT_VALIDATION_DELTA_F1_MIN,
            "passed": float(candidate_valid["f1"] - current_valid["f1"]) >= ACCEPT_VALIDATION_DELTA_F1_MIN,
        },
        {
            "gate": "full_nk_fp_reduction_vs_current",
            "observed": 1.0 - safe_div(float(candidate_full["predicted_NK_FP"]), float(current_full["predicted_NK_FP"])),
            "required": ACCEPT_FULL_NK_FP_REDUCTION_MIN,
            "passed": 1.0 - safe_div(float(candidate_full["predicted_NK_FP"]), float(current_full["predicted_NK_FP"])) >= ACCEPT_FULL_NK_FP_REDUCTION_MIN,
        },
        {
            "gate": "full_paired_tcrab_fp_no_worse_than_current",
            "observed": safe_div(float(candidate_full["predicted_paired_TCRAB_FP"]), float(current_full["predicted_paired_TCRAB_FP"])),
            "required": ACCEPT_FULL_TCRAB_FP_MAX_VS_CURRENT,
            "passed": float(candidate_full["predicted_paired_TCRAB_FP"]) <= float(current_full["predicted_paired_TCRAB_FP"]),
        },
        {
            "gate": "full_known_fp_fraction_lower_than_current",
            "observed": float(candidate_full["known_FP_fraction_of_predictions"] - current_full["known_FP_fraction_of_predictions"]),
            "required": ACCEPT_FULL_KNOWN_FP_FRACTION_DELTA_MAX,
            "passed": float(candidate_full["known_FP_fraction_of_predictions"]) < float(current_full["known_FP_fraction_of_predictions"]),
        },
    ]
    out = pd.DataFrame(checks)
    out.to_csv(TABLE_DIR / "promotion_acceptance_gates.csv", index=False)
    return bool(out["passed"].all()), out


def save_wrapper_if_accepted(accepted: bool, model: base.LoadedModel, model_path: Path, params: AnnotationThresholds) -> Path | None:
    if not accepted:
        return None
    payload = {
        "wrapper_name": "gdTAI_annotation_specific_cascade",
        "base_model_path": str(model_path),
        "base_model": {
            "model": model.model_name,
            "model_object": model.model_object,
            "threshold": model.threshold,
            "gene_names": model.gene_names,
            "feature_names": model.feature_names,
            "notes": model.notes,
        },
        "annotation_thresholds": asdict(params),
        "annotation_policy": {
            "gdT_cell": "current threshold",
            "CD8_T": "tuned near current threshold to preserve recall",
            "CD4_T": "tuned stricter compartment threshold",
            "Treg": "rescue-only or disabled by threshold",
            "NK_cell": "strict rescue-only threshold",
            "other": "strict rescue-only threshold",
        },
    }
    with WRAPPER_PKL.open("wb") as handle:
        pickle.dump(payload, handle)
    return WRAPPER_PKL


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


def plot_validation_bars(validation_df: pd.DataFrame) -> Path:
    path = FIGURE_DIR / "validation_metric_comparison.png"
    metrics = ["precision", "recall", "specificity", "f1"]
    colors = {"current_gdTAI": "#2563eb", "annotation_specific_cascade": "#0f766e"}
    fig, ax = plt.subplots(figsize=(6.8, 4.3), constrained_layout=True)
    x = np.arange(len(metrics))
    for i, (_, row) in enumerate(validation_df.iterrows()):
        ax.bar(x + (i - 0.5) * 0.28, [row[m] for m in metrics], width=0.28, color=colors.get(row["strategy"], "#64748b"), label=row["strategy"])
    ax.set_xticks(x, labels=metrics)
    ax.set_ylim(0, 1.04)
    ax.set_ylabel("Metric value")
    ax.set_title("Strict validation performance")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_confusion(validation_df: pd.DataFrame, score: np.ndarray, annotation: np.ndarray, y_true: np.ndarray, params: AnnotationThresholds, current_threshold: float) -> Path:
    path = FIGURE_DIR / "validation_confusion_current_vs_annotation_specific.png"
    current_pred = score >= current_threshold
    candidate_pred = predict_annotation_specific(score, annotation, params)
    panels = [
        ("Current gdTAI", confusion_matrix(y_true, current_pred.astype(int), labels=[0, 1])),
        ("Annotation-specific cascade", confusion_matrix(y_true, candidate_pred.astype(int), labels=[0, 1])),
    ]
    vmax = max(int(cm.max()) for _, cm in panels)
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.8), constrained_layout=True)
    for ax, (title, cm) in zip(axes, panels):
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


def plot_fp_comparison(comparison_df: pd.DataFrame) -> Path:
    path = FIGURE_DIR / "full_atlas_known_fp_comparison.png"
    keep = comparison_df.drop_duplicates("strategy").copy()
    keep = keep.sort_values("predicted_paired_TCRAB_or_NK_FP", ascending=False)
    fig, ax = plt.subplots(figsize=(8.8, max(4.0, 0.42 * len(keep) + 1.2)), constrained_layout=True)
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


def plot_annotation_counts(annotation_df: pd.DataFrame) -> Path:
    path = FIGURE_DIR / "prediction_counts_by_annotation.png"
    df = annotation_df.loc[annotation_df["strategy"].isin(["current_gdTAI", "annotation_specific_cascade"])].copy()
    piv = df.pivot_table(index="annotation", columns="strategy", values="predicted_putative_gdT", aggfunc="sum").fillna(0)
    piv = piv.reindex([a for a in ANNOTATION_ORDER if a in piv.index] + [a for a in piv.index if a not in ANNOTATION_ORDER])
    fig, ax = plt.subplots(figsize=(7.5, 4.5), constrained_layout=True)
    x = np.arange(len(piv.index))
    width = 0.34
    ax.bar(x - width / 2, piv.get("current_gdTAI", pd.Series(0, index=piv.index)), width=width, color="#2563eb", label="current")
    ax.bar(x + width / 2, piv.get("annotation_specific_cascade", pd.Series(0, index=piv.index)), width=width, color="#0f766e", label="annotation-specific")
    ax.set_xticks(x, labels=piv.index, rotation=25, ha="right")
    ax.set_ylabel("Predicted gdT cells")
    ax.set_title("Full-atlas predictions by scVI simple annotation")
    ax.legend(frameon=False)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def bounded_limits(x: np.ndarray, y: np.ndarray) -> tuple[tuple[float, float], tuple[float, float]]:
    finite = np.isfinite(x) & np.isfinite(y)
    if not finite.any():
        return (-1, 1), (-1, 1)
    x_lo, x_hi = np.nanpercentile(x[finite], [0.2, 99.8])
    y_lo, y_hi = np.nanpercentile(y[finite], [0.2, 99.8])
    return (float(x_lo - max((x_hi - x_lo) * 0.08, 0.05)), float(x_hi + max((x_hi - x_lo) * 0.08, 0.05))), (float(y_lo - max((y_hi - y_lo) * 0.08, 0.05)), float(y_hi + max((y_hi - y_lo) * 0.08, 0.05)))


def plot_trd_trab(handle: h5py.File, arrays: dict[str, np.ndarray], annotation: np.ndarray, *, sample_cells: int, seed: int) -> Path:
    path = FIGURE_DIR / "trd_vs_trab_current_vs_annotation_specific.png"
    n_obs = arrays["current_score"].shape[0]
    rng = np.random.default_rng(seed)
    random_idx = rng.choice(n_obs, size=min(sample_cells, n_obs), replace=False)
    current_idx = np.flatnonzero(arrays["current_pred"])
    candidate_idx = np.flatnonzero(arrays["annotation_specific_pred"])
    if current_idx.size > sample_cells:
        current_idx = rng.choice(current_idx, size=sample_cells, replace=False)
    if candidate_idx.size > sample_cells:
        candidate_idx = rng.choice(candidate_idx, size=sample_cells, replace=False)
    sample_idx = np.unique(np.concatenate([random_idx, current_idx, candidate_idx])).astype(np.int64)
    sample_idx.sort()
    trab = read_float_obs(handle, "phase4_trab_score")[sample_idx]
    trd = read_float_obs(handle, "phase4_trd_score")[sample_idx]
    current = arrays["current_pred"][sample_idx]
    candidate = arrays["annotation_specific_pred"][sample_idx]
    category = np.full(sample_idx.size, "neither", dtype=object)
    category[current & (~candidate)] = "current only"
    category[candidate & (~current)] = "annotation-specific only"
    category[current & candidate] = "both"
    scatter_df = pd.DataFrame({
        "obs_index": sample_idx,
        "annotation": annotation[sample_idx],
        "phase4_trab_score": trab,
        "phase4_trd_score": trd,
        "current_gdTAI_score": arrays["current_score"][sample_idx],
        "current_gdTAI_predicted": current,
        "annotation_specific_predicted": candidate,
        "method_category": category,
    })
    scatter_df.to_csv(TABLE_DIR / "trd_vs_trab_current_vs_annotation_specific_scatter_sample.csv.gz", index=False)
    xlim, ylim = bounded_limits(trab, trd)
    order = ["neither", "current only", "annotation-specific only", "both"]
    colors = {"neither": "#cbd5e1", "current only": "#2563eb", "annotation-specific only": "#f59e0b", "both": "#0f766e"}
    fig, ax = plt.subplots(figsize=(7.4, 6.0), constrained_layout=True)
    for cat in order:
        mask = category == cat
        if not mask.any():
            continue
        ax.scatter(trab[mask], trd[mask], s=2.0 if cat == "neither" else 4.0, alpha=0.18 if cat == "neither" else 0.70, linewidths=0, rasterized=True, color=colors[cat], label=f"{cat} (n={int(mask.sum()):,})")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("Phase 4 TRAB score")
    ax.set_ylabel("Phase 4 TRD score")
    ax.set_title("TRD vs TRAB: current vs annotation-specific gdTAI")
    ax.legend(frameon=False, fontsize=7, loc="best")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_umap(handle: h5py.File, arrays: dict[str, np.ndarray], *, seed: int, background_cells: int, positive_cells: int) -> Path | None:
    if "obsm" not in handle or "X_umap" not in handle["obsm"]:
        return None
    path = FIGURE_DIR / "umap_annotation_specific_cascade_predictions.png"
    n_obs = arrays["current_score"].shape[0]
    rng = np.random.default_rng(seed)
    bg_idx = rng.choice(n_obs, size=min(background_cells, n_obs), replace=False)
    candidate_idx = np.flatnonzero(arrays["annotation_specific_pred"])
    current_only_idx = np.flatnonzero(arrays["current_pred"] & (~arrays["annotation_specific_pred"]))
    if candidate_idx.size > positive_cells:
        candidate_idx = rng.choice(candidate_idx, size=positive_cells, replace=False)
    if current_only_idx.size > positive_cells // 2:
        current_only_idx = rng.choice(current_only_idx, size=positive_cells // 2, replace=False)
    umap = handle["obsm"]["X_umap"]
    bg = umap[np.sort(bg_idx), :]
    candidate = umap[np.sort(candidate_idx), :] if candidate_idx.size else np.empty((0, 2))
    current_only = umap[np.sort(current_only_idx), :] if current_only_idx.size else np.empty((0, 2))
    fig, ax = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    ax.scatter(bg[:, 0], bg[:, 1], s=0.35, color="#d1d5db", alpha=0.18, linewidths=0, rasterized=True, label="atlas background")
    if current_only.size:
        ax.scatter(current_only[:, 0], current_only[:, 1], s=0.7, color="#2563eb", alpha=0.32, linewidths=0, rasterized=True, label="current only")
    if candidate.size:
        ax.scatter(candidate[:, 0], candidate[:, 1], s=0.9, color="#0f766e", alpha=0.55, linewidths=0, rasterized=True, label="annotation-specific gdT")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("5 million cell atlas: annotation-specific gdTAI")
    ax.legend(frameon=False, fontsize=8, loc="best")
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
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
    chrome_candidates = [Path("/usr/bin/google-chrome"), Path("/usr/bin/google-chrome-stable"), shutil.which("google-chrome"), shutil.which("google-chrome-stable")]
    chrome = next((Path(path) for path in chrome_candidates if path and Path(path).exists()), None)
    if chrome is None:
        logging.warning("Skipping PDF export because google-chrome was not found")
        return
    subprocess.run([
        str(chrome),
        "--headless",
        "--disable-gpu",
        "--no-sandbox",
        "--print-to-pdf-no-header",
        f"--print-to-pdf={REPORT_PDF}",
        str(REPORT_HTML),
    ], check=True)


def write_protocol(model: base.LoadedModel, params: AnnotationThresholds, accepted: bool, wrapper_path: Path | None) -> None:
    lines = [
        "# gdTAI Annotation-Specific Cascade Protocol",
        "",
        "This wrapper keeps the selected individual TCR-gene gdTAI model frozen and applies separate score thresholds by `simple_annotation_plus6`.",
        "",
        "## Frozen Model",
        "",
        f"- Base model: `{model.model_name}`",
        f"- Base threshold: `{model.threshold}`",
        f"- Base feature genes: `{len(model.gene_names)}`",
        "",
        "## Annotation Thresholds",
        "",
        dataframe_to_markdown(pd.DataFrame([params_for_csv(params)])),
        "",
        "## Policy",
        "",
        "- `gdT_cell` and `CD8_T` keep recovery-oriented thresholds.",
        "- `CD4_T` uses a stricter tuned threshold.",
        "- `Treg`, `NK_cell`, and `other` are rescue-only or disabled by high thresholds.",
        "- This wrapper does not use Phase 4 TRD/TRAB scores for prediction decisions.",
        "",
        "## Status",
        "",
        f"- Promotion accepted: `{accepted}`",
        f"- Wrapper artifact: `{wrapper_path}`" if wrapper_path else "- Wrapper artifact: not written.",
    ]
    PROTOCOL_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_reports(
    input_h5ad: Path,
    model: base.LoadedModel,
    params: AnnotationThresholds,
    tune_df: pd.DataFrame,
    validation_df: pd.DataFrame,
    split_overall: pd.DataFrame,
    fp_overall: pd.DataFrame,
    comparison_df: pd.DataFrame,
    source_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    acceptance_df: pd.DataFrame,
    copied_assets: list[str],
    accepted: bool,
    wrapper_path: Path | None,
) -> None:
    current_fp = fp_overall.loc[fp_overall["strategy"] == "current_gdTAI"].iloc[0]
    candidate_fp = fp_overall.loc[fp_overall["strategy"] == "annotation_specific_cascade"].iloc[0]
    current_valid = validation_df.loc[validation_df["strategy"] == "current_gdTAI"].iloc[0]
    candidate_valid = validation_df.loc[validation_df["strategy"] == "annotation_specific_cascade"].iloc[0]
    decision = "accepted and written as a wrapper artifact" if accepted else "not promoted; kept as a candidate report"
    decision_detail = (
        f"Validation recall {candidate_valid['recall']:.4f} vs current {current_valid['recall']:.4f}; "
        f"F1 {candidate_valid['f1']:.4f} vs current {current_valid['f1']:.4f}. "
        f"Full atlas predictions {int(candidate_fp['predicted_putative_gdT']):,} vs current {int(current_fp['predicted_putative_gdT']):,}; "
        f"NK FP {int(candidate_fp['predicted_NK_FP']):,} vs {int(current_fp['predicted_NK_FP']):,}; "
        f"paired TCRAB FP {int(candidate_fp['predicted_paired_TCRAB_FP']):,} vs {int(current_fp['predicted_paired_TCRAB_FP']):,}."
    )
    tune_display = compact_numeric(pick_columns(tune_df, ["strategy", "threshold", "n_cells", "n_positive", "n_negative", "predicted_positive", "tp", "fp", "tn", "fn", "precision", "recall", "specificity", "f1", "roc_auc", "pr_auc"]))
    valid_display = compact_numeric(pick_columns(validation_df, ["strategy", "threshold", "n_cells", "n_positive", "n_negative", "predicted_positive", "tp", "fp", "tn", "fn", "precision", "recall", "specificity", "f1", "roc_auc", "pr_auc"]))
    fp_display = compact_numeric(pick_columns(comparison_df, ["strategy", "threshold", "predicted_putative_gdT", "predicted_fraction", "predicted_paired_TCRAB_FP", "paired_TCRAB_FP_fraction_of_predictions", "predicted_NK_FP", "NK_FP_fraction_of_predictions", "predicted_paired_TCRAB_or_NK_FP", "known_FP_fraction_of_predictions"]))
    source_display = compact_numeric(pick_columns(source_df.loc[source_df["strategy"] == "annotation_specific_cascade"].head(30), ["source_gse_id", "total_cells", "predicted_putative_gdT", "predicted_fraction", "predicted_paired_TCRAB_FP", "predicted_NK_FP", "known_FP_fraction_of_predictions"]))
    annotation_display = compact_numeric(pick_columns(annotation_df.loc[annotation_df["strategy"] == "annotation_specific_cascade"], ["annotation", "total_cells", "predicted_putative_gdT", "predicted_fraction", "predicted_paired_TCRAB_FP", "predicted_NK_FP", "known_FP_fraction_of_predictions"]))
    param_display = pd.DataFrame([params_for_csv(params)])
    acceptance_display = compact_numeric(acceptance_df)
    md_lines = [
        "# gdTAI Annotation-Specific Cascade Report",
        "",
        "## Executive Summary",
        "",
        f"- Input atlas: `{input_h5ad}`",
        "- Frozen selected gdTAI model; no model retraining.",
        "- Per-annotation thresholds tuned on non-holdout tune cells only.",
        "- Phase 4 TRD/TRAB scores are used only for figures/comparison, not prediction decisions.",
        f"- Promotion decision: `{decision}`.",
        f"- {decision_detail}",
        f"- Wrapper artifact: `{wrapper_path}`" if wrapper_path else "- Wrapper artifact: not written.",
        "",
        "## Split Summary",
        "",
        dataframe_to_markdown(compact_numeric(split_overall)),
        "",
        "## Annotation Thresholds",
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
        "## Full Atlas Known FP Proxy Comparison",
        "",
        dataframe_to_markdown(fp_display),
        "",
        "## Annotation-Specific Predictions By Annotation",
        "",
        dataframe_to_markdown(annotation_display),
        "",
        "## Top Prediction Sources",
        "",
        dataframe_to_markdown(source_display),
        "",
        "## Output Files",
        "",
        f"- Tables: `{TABLE_DIR}`",
        f"- Figures: `{FIGURE_DIR}`",
        f"- HTML: `{REPORT_HTML}`",
        f"- PDF: `{REPORT_PDF}`",
        f"- Protocol: `{PROTOCOL_MD}`",
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
    th,td{border:1px solid #e5e7eb;padding:5px 7px;text-align:left;vertical-align:top} th{background:#f1f5f9}
    .table-wrap{overflow-x:auto;margin-bottom:12px} img{max-width:100%;height:auto;border:1px solid #e5e7eb;background:#fff}
    code{background:#f1f5f9;padding:1px 4px;border-radius:4px}
    @media print{main{max-width:none;padding:12mm} table{font-size:9px} h1{font-size:22px} h2{font-size:16px} .table-wrap{overflow:visible}}
    """
    html_doc = f"""<!doctype html><html><head><meta charset='utf-8'><title>gdTAI Annotation-Specific Cascade</title><style>{css}</style></head><body><main>
    <h1>gdTAI Annotation-Specific Cascade</h1>
    <p class='status'>Promotion decision: {html.escape(decision)}</p>
    <p>{html.escape(decision_detail)}</p>
    <h2>Protocol</h2><p>This wrapper keeps the frozen selected gdTAI model and changes only score thresholds by scVI simple annotation. It does not modify any H5AD.</p>
    <h2>Split Summary</h2><div class='table-wrap'>{dataframe_to_html(compact_numeric(split_overall))}</div>
    <h2>Annotation Thresholds</h2><div class='table-wrap'>{dataframe_to_html(param_display)}</div>
    <h2>Tune Metrics</h2><div class='table-wrap'>{dataframe_to_html(tune_display)}</div>
    <h2>Strict Validation Metrics</h2><div class='table-wrap'>{dataframe_to_html(valid_display)}</div>
    <h2>Acceptance Gates</h2><div class='table-wrap'>{dataframe_to_html(acceptance_display)}</div>
    <h2>Full Atlas Known FP Proxy Comparison</h2><div class='table-wrap'>{dataframe_to_html(fp_display)}</div>
    <h2>Predictions By Annotation</h2><div class='table-wrap'>{dataframe_to_html(annotation_display)}</div>
    <h2>Top Prediction Sources</h2><div class='table-wrap'>{dataframe_to_html(source_display)}</div>
    <h2>Figures</h2>{''.join(figures_html)}
    <h2>Artifacts</h2><p>Tables: <code>{html.escape(str(TABLE_DIR))}</code><br>Figures: <code>{html.escape(str(FIGURE_DIR))}</code><br>Markdown: <code>{html.escape(str(REPORT_MD))}</code><br>Wrapper: <code>{html.escape(str(wrapper_path)) if wrapper_path else 'not written'}</code></p>
    </main></body></html>"""
    REPORT_HTML.write_text(html_doc, encoding="utf-8")


def write_summary_json(input_h5ad: Path, model_path: Path, model: base.LoadedModel, params: AnnotationThresholds, accepted: bool, wrapper_path: Path | None, validation_df: pd.DataFrame, fp_overall: pd.DataFrame) -> None:
    summary = {
        "input_h5ad": str(input_h5ad),
        "current_model_path": str(model_path),
        "current_model": model.model_name,
        "current_threshold": model.threshold,
        "annotation_thresholds": params_for_csv(params),
        "promotion_accepted": accepted,
        "wrapper_path": str(wrapper_path) if wrapper_path else None,
        "validation_metrics": validation_df.to_dict(orient="records"),
        "full_false_positive_overall": fp_overall.to_dict(orient="records"),
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
    model = base.load_current_model(model_path)
    logging.info("Input H5AD: %s", input_h5ad)
    logging.info("Loaded frozen model %s threshold %.6f", model.model_name, model.threshold)
    with h5py.File(input_h5ad, "r") as handle:
        obs = base.build_obs_metadata(handle)
        splits = base.build_splits(obs, args.seed)
        annotation = clean_group_values(read_obs_column(handle, "simple_annotation_plus6"))
        arrays = load_or_compute_scores(handle, model, args.score_cache.resolve(), recompute=args.recompute_score_cache, chunk_size=args.chunk_size)
        score = arrays["current_score"]
        tune_y = (obs.class_code[splits.tune_idx] == 2).astype(np.int8)
        params, _grid, _current_tune = tune_thresholds(score[splits.tune_idx], annotation[splits.tune_idx], tune_y, model.threshold)
        logging.info("Selected annotation thresholds: %s", params)
        pd.DataFrame([params_for_csv(params)]).to_csv(TABLE_DIR / "selected_annotation_thresholds.csv", index=False)
        tune_df, validation_df = evaluate_splits(score, annotation, obs, splits, params, model.threshold)
        full_overall, source_df, tissue_df, fp_overall, fp_source_df, annotation_df, full_arrays = summarize_full(handle, arrays, annotation, params, model.threshold)
        comparison_df = add_comparison_tables(fp_overall)
        accepted, acceptance_df = evaluate_acceptance(validation_df, fp_overall)
        wrapper_path = save_wrapper_if_accepted(accepted, model, model_path, params)
        valid_idx = splits.validation_idx
        valid_y = (obs.class_code[valid_idx] == 2).astype(np.int8)
        figures = [
            plot_validation_bars(validation_df),
            plot_confusion(validation_df, score[valid_idx], annotation[valid_idx], valid_y, params, model.threshold),
            plot_fp_comparison(comparison_df),
            plot_annotation_counts(annotation_df),
            plot_trd_trab(handle, full_arrays, annotation, sample_cells=args.scatter_sample_cells, seed=args.seed),
            plot_umap(handle, full_arrays, seed=args.seed, background_cells=args.umap_background_cells, positive_cells=args.umap_positive_cells),
        ]
        copied_assets = copy_assets(figures)
        write_protocol(model, params, accepted, wrapper_path)
        write_reports(
            input_h5ad=input_h5ad,
            model=model,
            params=params,
            tune_df=tune_df,
            validation_df=validation_df,
            split_overall=splits.split_overall,
            fp_overall=fp_overall,
            comparison_df=comparison_df,
            source_df=source_df,
            annotation_df=annotation_df,
            acceptance_df=acceptance_df,
            copied_assets=copied_assets,
            accepted=accepted,
            wrapper_path=wrapper_path,
        )
        write_summary_json(input_h5ad, model_path, model, params, accepted, wrapper_path, validation_df, fp_overall)
    render_pdf(args.no_pdf)
    if input_h5ad.stat().st_mtime_ns != before_mtime_ns:
        raise RuntimeError("Input H5AD changed during read-only annotation-specific cascade analysis.")
    logging.info("Saved report: %s", REPORT_HTML)
    logging.info("Saved summary: %s", SUMMARY_JSON)


if __name__ == "__main__":
    main()
