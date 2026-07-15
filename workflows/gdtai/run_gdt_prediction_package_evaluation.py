#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Evaluate package-style gdT score classifiers on the plus6 milestone.

This script is intentionally read-only with respect to the input H5AD. It reads
metadata and Phase 4 scores directly from HDF5, builds project-specific
gold/silver labels from TCR metadata, evaluates the shared package score rules,
and renders static report outputs.
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
import math
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    precision_recall_curve,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT_H5AD = PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "integrated_plus6.h5ad"
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
TABLE_DIR = OUTPUT_ROOT / "tables" / "gdT_prediction"
FIGURE_DIR = OUTPUT_ROOT / "figures" / "gdT_prediction"
LOG_DIR = OUTPUT_ROOT / "logs" / "gdT_prediction"
STATIC_DIR = PROJECT_ROOT / "gdT_prediction"
STATIC_ASSET_DIR = STATIC_DIR / "assets"

SUMMARY_MD = LOG_DIR / "gdT_prediction_summary.md"
RUN_LOG = LOG_DIR / "gdT_prediction_run.log"
SELECTED_JSON = LOG_DIR / "gdT_prediction_selected_model.json"
REPORT_HTML = STATIC_DIR / "index.html"
REPORT_PDF = STATIC_DIR / "gdT_prediction_report.pdf"

SORTED_GDT_SOURCES = {"GDTlung2023july_7p", "MalteGDT", "GDT_2020AUG_woCOV"}
TCR_EVALUATION_SOURCES = {"HRA005041", "GSE144469"}
GDT_BOOLEAN_REPAIR_SOURCES = {"GSE144469"}
INVALID_STRINGS = {"", "nan", "none", "na", "n/a", "<na>", "null"}
TRD_MINUS_TRAB_LEGACY_THRESHOLD = 0.4
FIGURE_DPI = 300
CONSERVATIVE_MIN_SPECIFICITY = 0.995
CONSERVATIVE_F_BETA = 0.5
OPERATING_POINT_SPECIFICITIES = (0.99, 0.995, 0.9975)

REQUIRED_COLUMNS = [
    "source_gse_id",
    "library_id",
    "sample_id",
    "has_TRA_TRB_paired",
    "has_TRG_TRD_paired",
    "has_any_ab_tcr",
    "has_any_gd_tcr",
    "TRG_cdr3",
    "TRD_cdr3",
    "phase4_trd_score",
    "phase4_trab_score",
    "phase4_trd_minus_trab",
]


@dataclass
class ScoreModel:
    model: str
    score_label: str
    score: np.ndarray
    threshold: float | None = None
    youden_threshold: float | None = None
    raw_cutoff: float | None = None
    y_pred_primary: np.ndarray | None = None
    y_pred_full: np.ndarray | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate package-style gdT prediction thresholds on integrated_plus6.h5ad."
    )
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--cluster-col", default="leiden")
    parser.add_argument("--n-top-clusters", type=int, default=3)
    parser.add_argument("--scatter-sample", type=int, default=160_000)
    parser.add_argument("--umap-background-sample", type=int, default=350_000)
    parser.add_argument("--umap-target-sample", type=int, default=250_000)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--no-pdf", action="store_true", help="Skip headless Chrome PDF export.")
    return parser.parse_args()


def setup_logging() -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    handlers = [
        logging.FileHandler(RUN_LOG, mode="w", encoding="utf-8"),
        logging.StreamHandler(),
    ]
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
        force=True,
    )


def ensure_dirs() -> None:
    for path in [TABLE_DIR, FIGURE_DIR, LOG_DIR, STATIC_DIR, STATIC_ASSET_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def decode_attr(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def read_string_dataset(dataset: h5py.Dataset) -> np.ndarray:
    try:
        return np.asarray(dataset.asstr()[:], dtype=object)
    except Exception:
        values = dataset[:]
        return np.asarray(
            [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values],
            dtype=object,
        )


def read_obs_column(handle: h5py.File, column: str) -> np.ndarray:
    obj = handle["obs"][column]
    if isinstance(obj, h5py.Group):
        encoding_type = decode_attr(obj.attrs.get("encoding-type", ""))
        if encoding_type == "categorical":
            categories = read_string_dataset(obj["categories"])
            codes = obj["codes"][:]
            out = np.full(codes.shape, "", dtype=object)
            valid = codes >= 0
            out[valid] = categories[codes[valid]]
            return out
        if "values" in obj and "mask" in obj:
            values = obj["values"][:]
            mask = obj["mask"][:].astype(bool)
            if values.dtype.kind in {"S", "O", "U"}:
                out = np.asarray(
                    [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values],
                    dtype=object,
                )
                out[mask] = ""
                return out
            out = np.asarray(values)
            if mask.any():
                out = out.astype(object)
                out[mask] = np.nan
            return out
        raise TypeError(f"Unsupported H5AD obs group encoding for column `{column}`: {encoding_type}")

    values = obj[:]
    if values.dtype.kind in {"S", "O", "U"}:
        return read_string_dataset(obj)
    return np.asarray(values)


def normalize_strings(values: np.ndarray, *, invalid_to: str = "") -> np.ndarray:
    series = pd.Series(values, copy=False).astype("string").fillna("").str.strip()
    invalid = series.str.lower().isin(INVALID_STRINGS)
    if invalid_to:
        series = series.mask(invalid, invalid_to)
    else:
        series = series.mask(invalid, "")
    return series.astype(object).to_numpy()


def clean_group_values(values: np.ndarray) -> np.ndarray:
    return normalize_strings(values, invalid_to="unknown")


def read_bool_obs(handle: h5py.File, column: str) -> np.ndarray:
    values = read_obs_column(handle, column)
    if np.issubdtype(values.dtype, np.bool_):
        return values.astype(bool, copy=False)
    if np.issubdtype(values.dtype, np.number):
        return values.astype(float) != 0
    lowered = pd.Series(values, copy=False).astype("string").fillna("").str.strip().str.lower()
    return lowered.isin({"true", "1", "yes", "y", "t"}).to_numpy(dtype=bool)


def read_nonempty_string_mask(handle: h5py.File, column: str) -> np.ndarray:
    obj = handle["obs"][column]
    if isinstance(obj, h5py.Group) and decode_attr(obj.attrs.get("encoding-type", "")) == "categorical":
        categories = read_string_dataset(obj["categories"])
        category_nonempty = (
            pd.Series(categories, copy=False).astype("string").fillna("").str.strip().str.lower().isin(INVALID_STRINGS)
            == False
        ).to_numpy(dtype=bool)
        codes = obj["codes"][:]
        out = codes >= 0
        clipped = np.clip(codes, 0, len(category_nonempty) - 1)
        out &= category_nonempty[clipped]
        return out
    values = normalize_strings(read_obs_column(handle, column))
    return pd.Series(values, copy=False).astype("string").str.lower().isin(INVALID_STRINGS).to_numpy(dtype=bool) == False


def read_float_obs(handle: h5py.File, column: str) -> np.ndarray:
    return np.asarray(read_obs_column(handle, column), dtype=np.float32)


def required_columns_present(handle: h5py.File) -> None:
    present = set(handle["obs"].keys())
    missing = [column for column in REQUIRED_COLUMNS if column not in present]
    if missing:
        raise KeyError(f"Missing required obs columns: {missing}")


def safe_div(numerator: float, denominator: float) -> float:
    if denominator == 0:
        return float("nan")
    return float(numerator / denominator)


def best_youden_threshold(y_true: np.ndarray, score: np.ndarray) -> dict[str, float]:
    fpr, tpr, thresholds = roc_curve(y_true, score)
    youden = tpr - fpr
    finite = np.isfinite(thresholds)
    if finite.any():
        finite_indices = np.flatnonzero(finite)
        idx = int(finite_indices[np.argmax(youden[finite])])
    else:
        idx = int(np.argmax(youden))
    return {
        "threshold": float(thresholds[idx]),
        "sensitivity": float(tpr[idx]),
        "specificity": float(1.0 - fpr[idx]),
        "youden_j": float(youden[idx]),
    }


def threshold_metrics_from_curve(
    y_true: np.ndarray,
    score: np.ndarray,
    *,
    min_specificity: float,
    beta: float,
) -> dict[str, float | int]:
    fpr, tpr, thresholds = roc_curve(y_true, score)
    n_positive = int(y_true.sum())
    n_negative = int((y_true == 0).sum())
    tp = np.rint(tpr * n_positive).astype(int)
    fp = np.rint(fpr * n_negative).astype(int)
    fn = n_positive - tp
    tn = n_negative - fp
    specificity = np.divide(tn, tn + fp, out=np.full_like(tn, np.nan, dtype=float), where=(tn + fp) > 0)
    recall = np.divide(tp, tp + fn, out=np.zeros_like(tp, dtype=float), where=(tp + fn) > 0)
    precision = np.divide(tp, tp + fp, out=np.zeros_like(tp, dtype=float), where=(tp + fp) > 0)
    beta2 = beta * beta
    fbeta = np.divide(
        (1.0 + beta2) * precision * recall,
        (beta2 * precision) + recall,
        out=np.zeros_like(precision, dtype=float),
        where=((beta2 * precision) + recall) > 0,
    )
    valid = np.isfinite(thresholds) & (specificity >= min_specificity) & ((tp + fp) > 0)
    if not valid.any():
        valid = np.isfinite(thresholds) & ((tp + fp) > 0)
    if not valid.any():
        raise RuntimeError("No finite threshold with predicted positives was available.")

    valid_indices = np.flatnonzero(valid)
    order = np.lexsort(
        (
            thresholds[valid_indices],
            recall[valid_indices],
            precision[valid_indices],
            fbeta[valid_indices],
        )
    )
    idx = int(valid_indices[order[-1]])
    return {
        "threshold": float(thresholds[idx]),
        "precision": float(precision[idx]),
        "recall": float(recall[idx]),
        "specificity": float(specificity[idx]),
        "f_beta": float(fbeta[idx]),
        "tp": int(tp[idx]),
        "fp": int(fp[idx]),
        "tn": int(tn[idx]),
        "fn": int(fn[idx]),
    }


def binary_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    score: np.ndarray | None = None,
) -> dict[str, float | int]:
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    has_both_true_classes = np.unique(y_true).size == 2
    out: dict[str, float | int] = {
        "n_cells": int(len(y_true)),
        "n_positive": int(y_true.sum()),
        "n_negative": int((y_true == 0).sum()),
        "predicted_positive": int(y_pred.sum()),
        "tp": int(tp),
        "fp": int(fp),
        "tn": int(tn),
        "fn": int(fn),
        "precision": float(precision_score(y_true, y_pred, zero_division=0)),
        "recall": float(recall_score(y_true, y_pred, zero_division=0)),
        "specificity": safe_div(tn, tn + fp),
        "f1": float(f1_score(y_true, y_pred, zero_division=0)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)) if has_both_true_classes else float("nan"),
        "mcc": float(matthews_corrcoef(y_true, y_pred)) if has_both_true_classes else float("nan"),
    }
    if score is not None and has_both_true_classes:
        out["roc_auc"] = float(roc_auc_score(y_true, score))
        out["pr_auc"] = float(average_precision_score(y_true, score))
    else:
        out["roc_auc"] = float("nan")
        out["pr_auc"] = float("nan")
    return out


def dataframe_to_markdown(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if max_rows is not None:
        df = df.head(max_rows)
    if df.empty:
        return "_No rows._"
    columns = [str(column) for column in df.columns]
    rows = [["" if pd.isna(value) else str(value) for value in row] for row in df.itertuples(index=False, name=None)]
    widths = [len(column) for column in columns]
    for row in rows:
        for idx, value in enumerate(row):
            widths[idx] = max(widths[idx], len(value))

    def format_row(row_values: list[str]) -> str:
        return "| " + " | ".join(value.ljust(widths[idx]) for idx, value in enumerate(row_values)) + " |"

    header = format_row(columns)
    separator = "| " + " | ".join("-" * width for width in widths) + " |"
    body = [format_row(row) for row in rows]
    return "\n".join([header, separator, *body])


def dataframe_to_html(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if max_rows is not None:
        df = df.head(max_rows)
    if df.empty:
        return "<p class='empty-table'>No rows.</p>"
    return df.to_html(index=False, escape=True, border=0, classes="dataframe")


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, list):
        return [json_ready(v) for v in value]
    if isinstance(value, tuple):
        return [json_ready(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        if math.isnan(float(value)):
            return None
        return float(value)
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value, float) and math.isnan(value):
        return None
    return value


def build_corrected_tcr_evidence(
    source: np.ndarray,
    has_TRA_TRB_paired: np.ndarray,
    has_TRG_TRD_paired_raw: np.ndarray,
    has_any_ab_tcr: np.ndarray,
    has_any_gd_tcr_raw: np.ndarray,
    trg_nonempty: np.ndarray,
    trd_nonempty: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """Return report-local gdTCR booleans after cleaning known missing-token artifacts.

    GSE144469 carries literal "NA" in missing TRD_cdr3 cells. The plus6
    harmonizer interpreted those strings as nonempty, which inflated
    has_any_gd_tcr to all cells. The integrated object no longer carries the
    original has_TRG/has_TRD columns, but for this dataset the raw
    has_TRG_TRD_paired column preserves the original TRG-evidence mask because
    TRD was inflated to all cells. We therefore use it as a TRG proxy and
    require clean TRD_cdr3 evidence for paired gdTCR.
    """
    corrected_any_gd = has_any_gd_tcr_raw.copy()
    corrected_paired_gd = has_TRG_TRD_paired_raw.copy()

    for source_id in GDT_BOOLEAN_REPAIR_SOURCES:
        mask = source == source_id
        if not mask.any():
            continue
        trg_proxy = has_TRG_TRD_paired_raw[mask] | trg_nonempty[mask]
        corrected_any_gd[mask] = trg_proxy | trd_nonempty[mask]
        corrected_paired_gd[mask] = trg_proxy & trd_nonempty[mask]

    audit_rows: list[dict[str, Any]] = []
    for source_id in sorted(TCR_EVALUATION_SOURCES):
        mask = source == source_id
        if not mask.any():
            continue
        audit_rows.append(
            {
                "source_gse_id": source_id,
                "n_cells": int(mask.sum()),
                "raw_has_any_gd_tcr": int(has_any_gd_tcr_raw[mask].sum()),
                "corrected_has_any_gd_tcr": int(corrected_any_gd[mask].sum()),
                "raw_has_TRG_TRD_paired": int(has_TRG_TRD_paired_raw[mask].sum()),
                "corrected_has_TRG_TRD_paired": int(corrected_paired_gd[mask].sum()),
                "clean_TRG_cdr3_nonempty": int(trg_nonempty[mask].sum()),
                "clean_TRD_cdr3_nonempty": int(trd_nonempty[mask].sum()),
                "has_any_ab_tcr": int(has_any_ab_tcr[mask].sum()),
                "has_TRA_TRB_paired": int(has_TRA_TRB_paired[mask].sum()),
                "repair_applied": bool(source_id in GDT_BOOLEAN_REPAIR_SOURCES),
            }
        )
    audit_df = pd.DataFrame(audit_rows)
    audit_df.to_csv(TABLE_DIR / "tcr_evidence_correction_audit.csv", index=False)
    return corrected_any_gd, corrected_paired_gd, audit_df


def build_sublibrary_labels(
    source: np.ndarray,
    library_id: np.ndarray,
    sample_id: np.ndarray,
    corrected_has_any_gd_tcr: np.ndarray,
    trg_nonempty: np.ndarray,
    trd_nonempty: np.ndarray,
) -> tuple[np.ndarray, pd.DataFrame]:
    use_library = pd.Series(library_id, copy=False).astype("string").fillna("").str.strip()
    use_sample = pd.Series(sample_id, copy=False).astype("string").fillna("").str.strip()
    sublibrary = use_library.mask(use_library.str.lower().isin(INVALID_STRINGS), use_sample).fillna("")
    sublibrary = sublibrary.mask(sublibrary.str.lower().isin(INVALID_STRINGS), "unknown").astype(str)

    target_source = np.isin(source, list(TCR_EVALUATION_SOURCES))
    evidence = corrected_has_any_gd_tcr | trg_nonempty | trd_nonempty
    target_df = pd.DataFrame(
        {
            "source_gse_id": source[target_source],
            "sublibrary_id": sublibrary.to_numpy(dtype=object)[target_source],
            "has_productive_gd_tcr_evidence": evidence[target_source],
        }
    )
    sublibrary_summary = (
        target_df.groupby(["source_gse_id", "sublibrary_id"], dropna=False, as_index=False)
        .agg(
            n_cells=("has_productive_gd_tcr_evidence", "size"),
            is_gdtcr_sequenced=("has_productive_gd_tcr_evidence", "any"),
        )
        .sort_values(["source_gse_id", "is_gdtcr_sequenced", "n_cells"], ascending=[True, False, False])
        .reset_index(drop=True)
    )

    target_keys = pd.Series(source[target_source], copy=False).astype(str) + "||" + target_df["sublibrary_id"].astype(str)
    sequenced_keys = set(target_keys[target_df["has_productive_gd_tcr_evidence"].to_numpy(dtype=bool)])
    all_keys = pd.Series(source, copy=False).astype(str) + "||" + sublibrary.astype(str)
    is_gdtcr_sequenced = target_source & all_keys.isin(sequenced_keys).to_numpy(dtype=bool)
    return is_gdtcr_sequenced, sublibrary_summary


def make_truth_labels(
    source: np.ndarray,
    is_gdtcr_sequenced_sublibrary: np.ndarray,
    has_TRA_TRB_paired: np.ndarray,
    corrected_has_TRG_TRD_paired: np.ndarray,
    has_any_ab_tcr: np.ndarray,
    corrected_has_any_gd_tcr: np.ndarray,
    trd_nonempty: np.ndarray,
) -> tuple[np.ndarray, pd.DataFrame]:
    sorted_gdt = np.isin(source, list(SORTED_GDT_SOURCES))
    tcr_eval_source = np.isin(source, list(TCR_EVALUATION_SOURCES))

    gdt_gold_rule_sorted = sorted_gdt
    gdt_gold_rule_tcr = (
        tcr_eval_source & is_gdtcr_sequenced_sublibrary & corrected_has_TRG_TRD_paired & (~has_any_ab_tcr)
    )
    gdt_gold = gdt_gold_rule_sorted | gdt_gold_rule_tcr

    abt_gold = tcr_eval_source & has_TRA_TRB_paired & (~corrected_has_any_gd_tcr)
    gdt_silver = (
        tcr_eval_source
        & is_gdtcr_sequenced_sublibrary
        & trd_nonempty
        & (~corrected_has_TRG_TRD_paired)
        & (~has_any_ab_tcr)
    )

    conflicts = gdt_gold & abt_gold
    conflict_df = pd.DataFrame(
        [
            {"conflict_type": "gdT_gold_and_abT_gold", "n_cells": int(conflicts.sum())},
            {"conflict_type": "gdT_gold_and_gdT_silver", "n_cells": int((gdt_gold & gdt_silver).sum())},
            {"conflict_type": "abT_gold_and_gdT_silver", "n_cells": int((abt_gold & gdt_silver).sum())},
        ]
    )

    class_code = np.zeros(source.shape[0], dtype=np.int8)
    class_code[abt_gold] = 1
    class_code[gdt_gold] = 2
    class_code[gdt_silver & (~gdt_gold) & (~abt_gold)] = 3
    return class_code, conflict_df


def categorical_truth(class_code: np.ndarray) -> pd.Categorical:
    categories = ["unlabeled_or_ambiguous", "abT_gold", "gdT_gold", "gdT_silver"]
    return pd.Categorical.from_codes(class_code.astype(int), categories=categories)


def write_ground_truth_tables(
    class_code: np.ndarray,
    source: np.ndarray,
    tissue: np.ndarray,
    sublibrary_summary: pd.DataFrame,
    conflict_df: pd.DataFrame,
    tcr_evidence_audit: pd.DataFrame,
) -> dict[str, pd.DataFrame]:
    truth = categorical_truth(class_code)
    class_counts = (
        pd.Series(truth, name="ground_truth_class")
        .value_counts(dropna=False)
        .rename_axis("ground_truth_class")
        .reset_index(name="n_cells")
    )
    class_counts["fraction_of_all_cells"] = class_counts["n_cells"] / int(class_code.size)

    by_source = (
        pd.DataFrame({"source_gse_id": source, "ground_truth_class": truth})
        .groupby(["source_gse_id", "ground_truth_class"], observed=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["source_gse_id", "ground_truth_class"])
    )
    by_tissue = (
        pd.DataFrame({"tissue": tissue, "ground_truth_class": truth})
        .groupby(["tissue", "ground_truth_class"], observed=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["tissue", "ground_truth_class"])
    )
    sublibrary_counts = (
        sublibrary_summary.groupby("source_gse_id", as_index=False)
        .agg(
            n_sublibraries=("sublibrary_id", "nunique"),
            n_gdtcr_sequenced_sublibraries=("is_gdtcr_sequenced", "sum"),
            n_cells=("n_cells", "sum"),
        )
        .sort_values("source_gse_id")
    )

    outputs = {
        "class_counts": class_counts,
        "by_source": by_source,
        "by_tissue": by_tissue,
        "sublibrary_summary": sublibrary_summary,
        "sublibrary_counts": sublibrary_counts,
        "conflicts": conflict_df,
        "tcr_evidence_audit": tcr_evidence_audit,
    }
    class_counts.to_csv(TABLE_DIR / "ground_truth_class_counts.csv", index=False)
    by_source.to_csv(TABLE_DIR / "ground_truth_by_source_gse.csv", index=False)
    by_tissue.to_csv(TABLE_DIR / "ground_truth_by_tissue.csv", index=False)
    sublibrary_summary.to_csv(TABLE_DIR / "gdtcr_sublibrary_summary.csv", index=False)
    sublibrary_counts.to_csv(TABLE_DIR / "gdtcr_sublibrary_counts.csv", index=False)
    conflict_df.to_csv(TABLE_DIR / "ground_truth_conflicts.csv", index=False)
    return outputs


def evaluate_score_models(
    models: list[ScoreModel],
    y_true: np.ndarray,
    primary_mask: np.ndarray,
) -> tuple[pd.DataFrame, ScoreModel]:
    rows: list[dict[str, Any]] = []
    y_primary = y_true[primary_mask]
    for model in models:
        score_primary = model.score[primary_mask]
        youden_info = best_youden_threshold(y_primary, score_primary)
        conservative_info = threshold_metrics_from_curve(
            y_primary,
            score_primary,
            min_specificity=CONSERVATIVE_MIN_SPECIFICITY,
            beta=CONSERVATIVE_F_BETA,
        )
        model.youden_threshold = float(youden_info["threshold"])
        model.threshold = float(conservative_info["threshold"])
        model.y_pred_primary = (score_primary >= model.threshold).astype(np.int8)
        model.y_pred_full = model.score >= model.threshold
        if model.model == "TRAB_score_low":
            model.raw_cutoff = -model.threshold
        else:
            model.raw_cutoff = model.threshold
        metrics = binary_metrics(y_primary, model.y_pred_primary, score_primary)
        rows.append(
            {
                "model": model.model,
                "score_used": model.score_label,
                "threshold_policy": (
                    f"specificity>={CONSERVATIVE_MIN_SPECIFICITY:g}; "
                    f"maximize_F{CONSERVATIVE_F_BETA:g}"
                ),
                "youden_threshold_on_score": youden_info["threshold"],
                "youden_sensitivity": youden_info["sensitivity"],
                "youden_specificity": youden_info["specificity"],
                "youden_j": youden_info["youden_j"],
                "selected_threshold_on_score": model.threshold,
                "raw_cutoff_for_original_score": model.raw_cutoff,
                "selected_specificity_constraint": CONSERVATIVE_MIN_SPECIFICITY,
                "selected_f_beta_beta": CONSERVATIVE_F_BETA,
                "selected_f_beta": conservative_info["f_beta"],
                **metrics,
                "full_dataset_predicted_positive": int(model.y_pred_full.sum()),
                "full_dataset_predicted_fraction": float(model.y_pred_full.mean()),
            }
        )

    threshold_df = pd.DataFrame(rows)
    threshold_df.to_csv(TABLE_DIR / "model_thresholds.csv", index=False)
    metrics_df = threshold_df.copy()
    metrics_df.to_csv(TABLE_DIR / "classifier_metrics_primary_gold.csv", index=False)

    ranked = threshold_df.sort_values(
        ["selected_f_beta", "precision", "pr_auc", "recall", "model"],
        ascending=[False, False, False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    best_name = str(ranked.loc[0, "model"])
    best_model = next(model for model in models if model.model == best_name)
    return threshold_df, best_model


def build_operating_point_grid(
    models: list[ScoreModel],
    y_true: np.ndarray,
    primary_mask: np.ndarray,
    min_specificities: tuple[float, ...],
) -> pd.DataFrame:
    y_primary = y_true[primary_mask]
    rows: list[dict[str, Any]] = []
    for model in models:
        score_primary = model.score[primary_mask]
        for min_specificity in min_specificities:
            info = threshold_metrics_from_curve(
                y_primary,
                score_primary,
                min_specificity=min_specificity,
                beta=CONSERVATIVE_F_BETA,
            )
            rows.append(
                {
                    "model": model.model,
                    "score_used": model.score_label,
                    "min_specificity": min_specificity,
                    "threshold_on_score": info["threshold"],
                    "raw_cutoff_for_original_score": -info["threshold"] if model.model == "TRAB_score_low" else info["threshold"],
                    "precision": info["precision"],
                    "recall": info["recall"],
                    "specificity": info["specificity"],
                    f"f{CONSERVATIVE_F_BETA:g}": info["f_beta"],
                    "tp": info["tp"],
                    "fp": info["fp"],
                    "tn": info["tn"],
                    "fn": info["fn"],
                    "full_dataset_predicted_positive": int((model.score >= float(info["threshold"])).sum()),
                    "full_dataset_predicted_fraction": float((model.score >= float(info["threshold"])).mean()),
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "conservative_operating_points.csv", index=False)
    return out


def evaluate_dual_and_cluster_rules(
    y_true: np.ndarray,
    primary_mask: np.ndarray,
    trd_score: np.ndarray,
    trab_score: np.ndarray,
    trd_minus_trab: np.ndarray,
    trd_threshold: float,
    trab_raw_cutoff: float,
    cluster_values: np.ndarray | None,
    n_top_clusters: int,
) -> tuple[pd.DataFrame, list[str]]:
    y_primary = y_true[primary_mask]
    rows: list[dict[str, Any]] = []

    dual_full = (trd_score >= trd_threshold) & (trab_score <= trab_raw_cutoff)
    dual_primary = dual_full[primary_mask].astype(np.int8)
    rows.append(
        {
            "model": "Dual_TRD_high_AND_TRAB_low",
            "score_used": f"phase4_trd_score >= {trd_threshold:.6g} AND phase4_trab_score <= {trab_raw_cutoff:.6g}",
            **binary_metrics(y_primary, dual_primary, None),
            "full_dataset_predicted_positive": int(dual_full.sum()),
            "full_dataset_predicted_fraction": float(dual_full.mean()),
        }
    )

    top_clusters: list[str] = []
    if cluster_values is not None:
        cluster_df = pd.DataFrame(
            {
                "cluster": pd.Series(cluster_values, copy=False).astype(str),
                "trd_minus_trab": trd_minus_trab,
            }
        )
        cluster_summary = (
            cluster_df.groupby("cluster", as_index=False)
            .agg(median_trd_minus_trab=("trd_minus_trab", "median"), n_cells=("trd_minus_trab", "size"))
            .sort_values(["median_trd_minus_trab", "n_cells"], ascending=[False, False])
        )
        cluster_summary.to_csv(TABLE_DIR / "cluster_trd_minus_trab_summary.csv", index=False)
        top_clusters = cluster_summary["cluster"].head(n_top_clusters).astype(str).tolist()
        cluster_full = pd.Series(cluster_values, copy=False).astype(str).isin(top_clusters).to_numpy(dtype=bool) | (
            trd_score >= trd_threshold
        )
        cluster_primary = cluster_full[primary_mask].astype(np.int8)
        rows.append(
            {
                "model": "Secondary_cluster_top_TRD_minus_TRAB_or_TRD_high",
                "score_used": f"top {n_top_clusters} clusters by median TRD-TRAB OR phase4_trd_score >= {trd_threshold:.6g}",
                **binary_metrics(y_primary, cluster_primary, None),
                "full_dataset_predicted_positive": int(cluster_full.sum()),
                "full_dataset_predicted_fraction": float(cluster_full.mean()),
            }
        )

    df = pd.DataFrame(rows)
    df.to_csv(TABLE_DIR / "secondary_rule_metrics_primary_gold.csv", index=False)
    return df, top_clusters


def evaluate_selected_by_group(
    group_values: np.ndarray,
    group_name: str,
    primary_mask: np.ndarray,
    y_true: np.ndarray,
    y_pred_full: np.ndarray,
    score: np.ndarray,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    groups = pd.Series(group_values[primary_mask], copy=False).astype(str)
    y = y_true[primary_mask]
    pred = y_pred_full[primary_mask].astype(np.int8)
    score_primary = score[primary_mask]
    for group, idx in groups.groupby(groups, sort=True).groups.items():
        idx_array = np.asarray(list(idx), dtype=int)
        y_group = y[idx_array]
        pred_group = pred[idx_array]
        score_group = score_primary[idx_array]
        metrics = binary_metrics(y_group, pred_group, score_group if np.unique(y_group).size == 2 else None)
        rows.append({group_name: group, **metrics})
    return pd.DataFrame(rows).sort_values(["n_cells", group_name], ascending=[False, True]).reset_index(drop=True)


def selected_counts_by_group(group_values: np.ndarray, group_name: str, y_pred_full: np.ndarray) -> pd.DataFrame:
    df = pd.DataFrame({group_name: group_values, "predicted_putative_gdT": y_pred_full.astype(bool)})
    out = (
        df.groupby(group_name, dropna=False, as_index=False)
        .agg(total_cells=("predicted_putative_gdT", "size"), predicted_putative_gdT=("predicted_putative_gdT", "sum"))
        .sort_values(["predicted_putative_gdT", "total_cells"], ascending=[False, False])
        .reset_index(drop=True)
    )
    out["predicted_fraction"] = out["predicted_putative_gdT"] / out["total_cells"]
    return out


def build_overlap_table(
    y_pred_full: np.ndarray,
    sorted_gdt: np.ndarray,
    has_TRG_TRD_paired: np.ndarray,
    has_any_ab_tcr: np.ndarray,
    trd_minus_trab: np.ndarray,
    simple_annotation: np.ndarray | None,
) -> pd.DataFrame:
    legacy_score = trd_minus_trab > TRD_MINUS_TRAB_LEGACY_THRESHOLD
    criteria: list[tuple[str, np.ndarray]] = [
        ("Sorted_gdT_true", sorted_gdt),
        ("has_TRG_TRD_paired", has_TRG_TRD_paired),
        ("has_TRG_TRD_paired_no_any_ab_tcr", has_TRG_TRD_paired & (~has_any_ab_tcr)),
        ("TRD_minus_TRAB_gt_0p4", legacy_score),
    ]
    if simple_annotation is not None:
        not_nk = pd.Series(simple_annotation, copy=False).astype(str).to_numpy() != "NK_cell"
        criteria.append(("TRD_minus_TRAB_gt_0p4_no_any_ab_tcr_not_NK", legacy_score & (~has_any_ab_tcr) & not_nk))

    rows = []
    for name, mask in criteria:
        rows.append(
            {
                "criterion": name,
                "criterion_cells": int(mask.sum()),
                "selected_predicted_putative_gdT_overlap": int((mask & y_pred_full).sum()),
                "fraction_of_criterion_predicted_putative_gdT": safe_div(float((mask & y_pred_full).sum()), float(mask.sum())),
                "fraction_of_selected_predictions_in_criterion": safe_div(
                    float((mask & y_pred_full).sum()), float(y_pred_full.sum())
                ),
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "overlap_with_reference_criteria.csv", index=False)
    return out


def sample_indices(mask: np.ndarray, max_n: int, rng: np.random.Generator) -> np.ndarray:
    indices = np.flatnonzero(mask)
    if len(indices) <= max_n:
        return indices
    return np.sort(rng.choice(indices, size=max_n, replace=False))


def plot_roc_pr(models: list[ScoreModel], primary_mask: np.ndarray, y_true: np.ndarray) -> None:
    y = y_true[primary_mask]

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for model in models:
        score = model.score[primary_mask]
        fpr, tpr, _ = roc_curve(y, score)
        auc_val = roc_auc_score(y, score)
        ax.plot(fpr, tpr, lw=2, label=f"{model.model} (AUC={auc_val:.3f})")
    ax.plot([0, 1], [0, 1], color="0.7", lw=1, linestyle="--")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("ROC curves on primary gold labels")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.2)
    fig.savefig(FIGURE_DIR / "roc_curves.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    prevalence = y.mean()
    for model in models:
        score = model.score[primary_mask]
        precision, recall, _ = precision_recall_curve(y, score)
        ap_val = average_precision_score(y, score)
        ax.plot(recall, precision, lw=2, label=f"{model.model} (AP={ap_val:.3f})")
    ax.axhline(prevalence, color="0.7", lw=1, linestyle="--", label=f"prevalence={prevalence:.3f}")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Precision-recall curves on primary gold labels")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.2)
    fig.savefig(FIGURE_DIR / "precision_recall_curves.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_score_distributions(models: list[ScoreModel], class_code: np.ndarray, rng: np.random.Generator) -> None:
    labels = [
        (2, "gdT_gold", "#b51f2a"),
        (1, "abT_gold", "#2d6f9f"),
        (3, "gdT_silver", "#5b8c3b"),
    ]
    fig, axes = plt.subplots(1, len(models), figsize=(18, 5), constrained_layout=True)
    for ax, model in zip(axes, models):
        for code, label, color in labels:
            idx = sample_indices(class_code == code, 200_000, rng)
            if len(idx) == 0:
                continue
            ax.hist(model.score[idx], bins=80, density=True, histtype="step", lw=2, color=color, label=f"{label} (n={len(idx):,})")
        if model.threshold is not None:
            ax.axvline(model.threshold, color="black", lw=1.5, linestyle="--", label="selected conservative cutoff")
        ax.set_title(model.model)
        ax.set_xlabel(model.score_label)
        ax.set_ylabel("Density")
        ax.grid(alpha=0.15)
        ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIGURE_DIR / "score_distributions_by_truth.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_confusion(y_true: np.ndarray, y_pred: np.ndarray, model_name: str) -> None:
    matrix = confusion_matrix(y_true, y_pred, labels=[1, 0])
    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    im = ax.imshow(matrix, cmap="Blues")
    ax.set_xticks([0, 1], labels=["Pred gdT", "Pred abT"])
    ax.set_yticks([0, 1], labels=["True gdT", "True abT"])
    ax.set_title(f"Confusion matrix: {model_name}")
    for (row, col), value in np.ndenumerate(matrix):
        ax.text(col, row, f"{int(value):,}", ha="center", va="center", color="black", fontsize=12)
    fig.colorbar(im, ax=ax, shrink=0.85)
    fig.savefig(FIGURE_DIR / "selected_confusion_matrix.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_trab_vs_trd_scatter(
    trd_score: np.ndarray,
    trab_score: np.ndarray,
    class_code: np.ndarray,
    selected_pred: np.ndarray,
    rng: np.random.Generator,
    max_sample: int,
) -> None:
    primary_or_silver = class_code > 0
    background = (class_code == 0) & (~selected_pred)
    unlabeled_pred = (class_code == 0) & selected_pred
    idx = np.unique(
        np.concatenate(
            [
                sample_indices(background, max_sample, rng),
                sample_indices(unlabeled_pred, max_sample // 2, rng),
                np.flatnonzero(primary_or_silver),
            ]
        )
    )

    code = class_code[idx]
    pred = selected_pred[idx]
    categories = [
        ("unlabeled/background", (code == 0) & (~pred), "#d8d8d8", 4, 0.25),
        ("unlabeled predicted gdT", (code == 0) & pred, "#444444", 5, 0.45),
        ("abT_gold TN", (code == 1) & (~pred), "#2d6f9f", 10, 0.65),
        ("abT_gold FP", (code == 1) & pred, "#f28e2b", 14, 0.9),
        ("gdT_gold TP", (code == 2) & pred, "#b51f2a", 14, 0.9),
        ("gdT_gold FN", (code == 2) & (~pred), "#7a5195", 18, 0.9),
        ("gdT_silver predicted", (code == 3) & pred, "#5b8c3b", 12, 0.8),
        ("gdT_silver not predicted", (code == 3) & (~pred), "#9fbf6a", 12, 0.65),
    ]

    fig, ax = plt.subplots(figsize=(7.5, 6.2))
    for label, mask, color, size, alpha in categories:
        if not mask.any():
            continue
        sub = idx[mask]
        ax.scatter(trab_score[sub], trd_score[sub], s=size, c=color, alpha=alpha, linewidths=0, rasterized=True, label=label)
    ax.set_xlabel("phase4_trab_score")
    ax.set_ylabel("phase4_trd_score")
    ax.set_title("TRAB versus TRD score space by truth and selected prediction")
    ax.legend(frameon=True, fontsize=7, loc="best")
    ax.grid(alpha=0.15)
    fig.savefig(FIGURE_DIR / "trab_vs_trd_truth_prediction_scatter.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_umap_predictions(
    input_h5ad: Path,
    selected_pred: np.ndarray,
    rng: np.random.Generator,
    background_n: int,
    target_n: int,
) -> None:
    with h5py.File(input_h5ad, "r") as handle:
        umap = np.asarray(handle["obsm"]["X_umap"][:], dtype=np.float32)
    background_idx = sample_indices(~selected_pred, background_n, rng)
    target_idx = sample_indices(selected_pred, target_n, rng)

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    ax.scatter(
        umap[background_idx, 0],
        umap[background_idx, 1],
        s=1.5,
        c="#d0d0d0",
        alpha=0.35,
        linewidths=0,
        rasterized=True,
        label=f"other cells sample (n={len(background_idx):,})",
    )
    ax.scatter(
        umap[target_idx, 0],
        umap[target_idx, 1],
        s=2.2,
        c="#b51f2a",
        alpha=0.75,
        linewidths=0,
        rasterized=True,
        label=f"predicted putative gdT sample (n={len(target_idx):,})",
    )
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("Whole plus6 UMAP: selected predicted putative gdT cells")
    ax.legend(frameon=True, fontsize=8, loc="best")
    fig.savefig(FIGURE_DIR / "umap_predicted_putative_gdt.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def plot_count_bars(source_counts: pd.DataFrame, tissue_counts: pd.DataFrame, overlap_df: pd.DataFrame) -> None:
    for df, group_name, filename, title, max_rows in [
        (source_counts, "source_gse_id", "predicted_counts_by_source_gse.png", "Predicted putative gdT cells by source", 25),
        (tissue_counts, "tissue", "predicted_counts_by_tissue.png", "Predicted putative gdT cells by tissue", 25),
    ]:
        plot_df = df.head(max_rows).iloc[::-1]
        fig, ax = plt.subplots(figsize=(8, max(5, 0.28 * len(plot_df) + 1.5)))
        ax.barh(plot_df[group_name].astype(str), plot_df["predicted_putative_gdT"], color="#2d6f9f")
        ax.set_xlabel("Predicted putative gdT cells")
        ax.set_title(title)
        ax.grid(axis="x", alpha=0.2)
        fig.savefig(FIGURE_DIR / filename, dpi=FIGURE_DPI, bbox_inches="tight")
        plt.close(fig)

    plot_df = overlap_df.iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, max(4.5, 0.45 * len(plot_df) + 1.5)))
    y = np.arange(len(plot_df))
    ax.barh(y - 0.18, plot_df["criterion_cells"], height=0.35, color="#b8b8b8", label="criterion total")
    ax.barh(
        y + 0.18,
        plot_df["selected_predicted_putative_gdT_overlap"],
        height=0.35,
        color="#b51f2a",
        label="overlap with selected predictions",
    )
    ax.set_yticks(y, labels=plot_df["criterion"])
    ax.set_xlabel("Cells")
    ax.set_title("Overlap with sorted, paired gdTCR, and legacy score criteria")
    ax.legend(frameon=False)
    ax.grid(axis="x", alpha=0.2)
    fig.savefig(FIGURE_DIR / "reference_overlap_barplot.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def write_summary_markdown(
    input_h5ad: Path,
    n_obs: int,
    best_model: ScoreModel,
    threshold_df: pd.DataFrame,
    operating_point_df: pd.DataFrame,
    secondary_df: pd.DataFrame,
    truth_tables: dict[str, pd.DataFrame],
    source_counts: pd.DataFrame,
    tissue_counts: pd.DataFrame,
    overlap_df: pd.DataFrame,
    silver_recall: float,
    top_clusters: list[str],
) -> None:
    selected_row = threshold_df.loc[threshold_df["model"] == best_model.model].iloc[0]
    lines = [
        "# gdT prediction package-style evaluation",
        "",
        f"- Input H5AD: `{input_h5ad}`",
        f"- Whole plus6 cells: `{n_obs:,}`",
        f"- Selected single-score model: `{best_model.model}`",
        f"- Selected conservative cutoff on score: `{float(best_model.threshold):.6g}`",
        f"- Threshold policy: specificity >= `{CONSERVATIVE_MIN_SPECIFICITY:.3f}` on the primary gold set, maximizing `F{CONSERVATIVE_F_BETA:g}`",
        f"- Youden cutoff for selected model, retained only as a diagnostic: `{float(best_model.youden_threshold):.6g}`",
        f"- Raw cutoff for original score: `{float(best_model.raw_cutoff):.6g}`",
        f"- Full-dataset predicted putative gdT cells: `{int(selected_row['full_dataset_predicted_positive']):,}`",
        f"- Full-dataset predicted putative gdT fraction: `{float(selected_row['full_dataset_predicted_fraction']):.4%}`",
        f"- Silver gdT sensitivity recall under selected cutoff: `{silver_recall:.4f}`",
        "",
        "## Ground Truth Definition",
        "",
        "Primary gdT-gold positives were assigned from all cells in sorted gdT datasets "
        "`GDTlung2023july_7p`, `MalteGDT`, and `GDT_2020AUG_woCOV`, plus cells from "
        "`HRA005041` and `GSE144469` that are in gdTCR-sequenced sublibraries with "
        "corrected `has_TRG_TRD_paired == True` and `has_any_ab_tcr == False`.",
        "",
        "Primary abT-gold negatives were assigned from `HRA005041` and `GSE144469` cells "
        "with `has_TRA_TRB_paired == True` and corrected `has_any_gd_tcr == False`.",
        "",
        "Sensitivity-only gdT-silver cells were assigned from `HRA005041` and `GSE144469` "
        "gdTCR-sequenced sublibraries when clean productive `TRD_cdr3` was present, "
        "corrected paired `TRG/TRD` was absent, and `has_any_ab_tcr == False`. Silver "
        "cells were excluded from threshold fitting and primary performance metrics.",
        "",
        "A sublibrary used `library_id` when available, otherwise `sample_id`. A sublibrary "
        "was marked gdTCR-sequenced when it contained any productive gdTCR evidence: "
        "corrected `has_any_gd_tcr == True` or clean nonempty `TRG_cdr3` / `TRD_cdr3`.",
        "",
        "Report-local TCR evidence correction treats literal missing tokens such as `NA`, "
        "`nan`, `None`, and empty strings as missing. For `GSE144469`, raw `has_any_gd_tcr` "
        "was inflated by `NA` in missing `TRD_cdr3`, so the report recalculates any gdTCR "
        "evidence from raw TRG-proxy evidence plus clean `TRD_cdr3`, and requires both for "
        "corrected paired gdTCR. The H5AD itself is not modified.",
        "",
        "## Ground Truth Audit",
        "",
        dataframe_to_markdown(truth_tables["class_counts"]),
        "",
        "### TCR evidence correction audit",
        "",
        dataframe_to_markdown(truth_tables["tcr_evidence_audit"]),
        "",
        "### gdTCR-sequenced sublibraries",
        "",
        dataframe_to_markdown(truth_tables["sublibrary_counts"]),
        "",
        "### Overlap conflicts",
        "",
        dataframe_to_markdown(truth_tables["conflicts"]),
        "",
        "## Primary Performance",
        "",
        "The final operating point is precision-prioritized. Youden thresholds are kept in the table for comparison, but final predictions use the selected conservative threshold.",
        "",
        dataframe_to_markdown(threshold_df),
        "",
        "## Conservative Operating Points",
        "",
        dataframe_to_markdown(operating_point_df),
        "",
        "## Secondary Rules",
        "",
        dataframe_to_markdown(secondary_df),
        "",
        f"Top clusters used by the secondary cluster rule: `{', '.join(top_clusters) if top_clusters else 'not available'}`",
        "",
        "## Whole-Dataset Prediction Counts",
        "",
        "### Top source GSEs",
        "",
        dataframe_to_markdown(source_counts, max_rows=25),
        "",
        "### Top tissues",
        "",
        dataframe_to_markdown(tissue_counts, max_rows=25),
        "",
        "## Reference Criteria Overlap",
        "",
        dataframe_to_markdown(overlap_df),
        "",
    ]
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_image_section(title: str, filename: str) -> str:
    return (
        f"<section class='figure-card'><h3>{html.escape(title)}</h3>"
        f"<img src='assets/{html.escape(filename)}' alt='{html.escape(title)}' class='zoomable-figure' "
        f"data-title='{html.escape(title)}'></section>"
    )


def copy_static_assets() -> list[str]:
    filenames = [
        "roc_curves.png",
        "precision_recall_curves.png",
        "score_distributions_by_truth.png",
        "selected_confusion_matrix.png",
        "trab_vs_trd_truth_prediction_scatter.png",
        "umap_predicted_putative_gdt.png",
        "predicted_counts_by_source_gse.png",
        "predicted_counts_by_tissue.png",
        "reference_overlap_barplot.png",
    ]
    copied = []
    for filename in filenames:
        src = FIGURE_DIR / filename
        if src.exists():
            shutil.copyfile(src, STATIC_ASSET_DIR / filename)
            copied.append(filename)
    return copied


def render_html_report(
    input_h5ad: Path,
    n_obs: int,
    best_model: ScoreModel,
    threshold_df: pd.DataFrame,
    operating_point_df: pd.DataFrame,
    secondary_df: pd.DataFrame,
    truth_tables: dict[str, pd.DataFrame],
    source_metrics: pd.DataFrame,
    tissue_metrics: pd.DataFrame,
    source_counts: pd.DataFrame,
    tissue_counts: pd.DataFrame,
    overlap_df: pd.DataFrame,
    prevalence_df: pd.DataFrame,
    silver_recall: float,
    copied_assets: list[str],
) -> None:
    selected_row = threshold_df.loc[threshold_df["model"] == best_model.model].iloc[0]
    asset_titles = {
        "roc_curves.png": "ROC curves",
        "precision_recall_curves.png": "Precision-recall curves",
        "score_distributions_by_truth.png": "Score distributions by truth class",
        "selected_confusion_matrix.png": "Selected cutoff confusion matrix",
        "trab_vs_trd_truth_prediction_scatter.png": "TRAB versus TRD truth/prediction scatter",
        "umap_predicted_putative_gdt.png": "UMAP highlight of predicted putative gdT cells",
        "predicted_counts_by_source_gse.png": "Predicted putative gdT counts by source GSE",
        "predicted_counts_by_tissue.png": "Predicted putative gdT counts by tissue",
        "reference_overlap_barplot.png": "Overlap with sorted, paired gdTCR, and legacy criteria",
    }

    figure_html = "\n".join(make_image_section(asset_titles[name], name) for name in copied_assets)
    css = """
    :root{--bg:#f6f7f8;--paper:#ffffff;--ink:#1e252b;--muted:#5d6975;--line:#d9dee3;--accent:#b51f2a;--accent2:#2d6f9f}
    *{box-sizing:border-box} body{margin:0;background:var(--bg);color:var(--ink);font-family:Arial,Helvetica,sans-serif;line-height:1.55}
    .page{width:min(1220px,calc(100vw - 32px));margin:20px auto 42px}
    .hero{background:var(--paper);border:1px solid var(--line);border-top:6px solid var(--accent);padding:28px 32px;margin-bottom:18px}
    h1{margin:0 0 10px;font-size:36px;letter-spacing:0;line-height:1.12} h2{margin:0 0 12px;font-size:24px} h3{margin:18px 0 8px;font-size:17px}
    p{margin:0 0 12px;color:var(--muted)} code{background:#eef1f4;padding:1px 5px;border-radius:4px}
    .metric-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:10px;margin-top:18px}
    .metric{border:1px solid var(--line);background:#fbfcfd;padding:12px 14px}
    .metric .value{display:block;font-size:24px;font-weight:700;color:var(--ink)} .metric .label{display:block;font-size:12px;text-transform:uppercase;color:var(--muted)}
    .section{background:var(--paper);border:1px solid var(--line);padding:24px;margin-top:18px}
    .figure-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(430px,1fr));gap:18px}
    .figure-card{border:1px solid var(--line);background:#fbfcfd;padding:14px}
    .figure-card img{width:100%;height:auto;border:1px solid var(--line);background:white;cursor:zoom-in}
    table{border-collapse:collapse;width:100%;font-size:12px;margin:8px 0 16px} th,td{border:1px solid var(--line);padding:5px 7px;text-align:left;vertical-align:top}
    th{background:#eef1f4}.empty-table{font-style:italic}
    .lightbox{position:fixed;inset:0;display:none;align-items:center;justify-content:center;padding:24px;background:rgba(0,0,0,0.84);z-index:20}
    .lightbox.open{display:flex}.lightbox img{max-width:96vw;max-height:88vh;background:white}.lightbox-caption{color:white;margin-top:8px}
    .lightbox-close{position:absolute;top:16px;right:18px;width:38px;height:38px;border:0;background:white;color:#111;font-size:26px;cursor:pointer}
    """

    html_parts = [
        "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'>",
        "<title>gdT prediction package-style evaluation</title>",
        f"<style>{css}</style></head><body><main class='page'>",
        "<section class='hero'>",
        "<h1>gdT Prediction Package-Style Evaluation</h1>",
        (
            "<p>Read-only evaluation of Phase 4 score classifiers on "
            f"<code>{html.escape(str(input_h5ad))}</code>, using project-defined TCR metadata labels.</p>"
        ),
        "<div class='metric-grid'>",
        f"<div class='metric'><span class='value'>{n_obs:,}</span><span class='label'>Plus6 cells</span></div>",
        f"<div class='metric'><span class='value'>{html.escape(best_model.model)}</span><span class='label'>Selected model</span></div>",
        f"<div class='metric'><span class='value'>{float(best_model.threshold):.4g}</span><span class='label'>Conservative cutoff</span></div>",
        f"<div class='metric'><span class='value'>{int(selected_row['full_dataset_predicted_positive']):,}</span><span class='label'>Predicted putative gdT</span></div>",
        f"<div class='metric'><span class='value'>{float(selected_row['full_dataset_predicted_fraction']):.3%}</span><span class='label'>Predicted fraction</span></div>",
        f"<div class='metric'><span class='value'>{silver_recall:.3f}</span><span class='label'>Silver recall</span></div>",
        "</div></section>",
        "<section class='section'><h2>Ground Truth Definition</h2>",
        "<p><strong>gdT-gold positives</strong>: all cells from <code>GDTlung2023july_7p</code>, <code>MalteGDT</code>, and <code>GDT_2020AUG_woCOV</code>; plus <code>HRA005041</code> and <code>GSE144469</code> cells in gdTCR-sequenced sublibraries with corrected <code>has_TRG_TRD_paired == True</code> and <code>has_any_ab_tcr == False</code>.</p>",
        "<p><strong>abT-gold negatives</strong>: <code>HRA005041</code> and <code>GSE144469</code> cells with <code>has_TRA_TRB_paired == True</code> and corrected <code>has_any_gd_tcr == False</code>.</p>",
        "<p><strong>gdT-silver sensitivity cells</strong>: <code>HRA005041</code> and <code>GSE144469</code> cells in gdTCR-sequenced sublibraries with clean productive <code>TRD_cdr3</code>, no corrected paired <code>TRG/TRD</code>, and <code>has_any_ab_tcr == False</code>. Silver cells are excluded from threshold fitting and primary metrics.</p>",
        "<p>Sublibraries use <code>library_id</code> when available, otherwise <code>sample_id</code>. A sublibrary is gdTCR-sequenced when any cell has corrected <code>has_any_gd_tcr == True</code> or clean nonempty <code>TRG_cdr3</code>/<code>TRD_cdr3</code>.</p>",
        "<p>Report-local TCR evidence correction treats literal missing tokens such as <code>NA</code>, <code>nan</code>, <code>None</code>, and empty strings as missing. For <code>GSE144469</code>, raw <code>has_any_gd_tcr</code> was inflated by <code>NA</code> in missing <code>TRD_cdr3</code>, so the report recalculates any gdTCR evidence from raw TRG-proxy evidence plus clean <code>TRD_cdr3</code>, and requires both for corrected paired gdTCR. The H5AD itself is not modified.</p>",
        "<h3>Class Counts</h3>",
        dataframe_to_html(truth_tables["class_counts"]),
        "<h3>TCR Evidence Correction Audit</h3>",
        dataframe_to_html(truth_tables["tcr_evidence_audit"]),
        "<h3>gdTCR-Sequenced Sublibraries</h3>",
        dataframe_to_html(truth_tables["sublibrary_counts"]),
        "<h3>Overlap Conflicts</h3>",
        dataframe_to_html(truth_tables["conflicts"]),
        "</section>",
        "<section class='section'><h2>Performance Results</h2>",
        (
            f"<p>Final predictions use a precision-prioritized conservative cutoff: specificity &ge; "
            f"{CONSERVATIVE_MIN_SPECIFICITY:.3f} on the primary gold set, maximizing F{CONSERVATIVE_F_BETA:g}. "
            "Youden thresholds are retained in the table as diagnostics only.</p>"
        ),
        dataframe_to_html(threshold_df),
        "<h3>Conservative Operating Points</h3>",
        dataframe_to_html(operating_point_df),
        "<h3>Secondary Rules</h3>",
        dataframe_to_html(secondary_df),
        "<h3>Prevalence-Aware Counts</h3>",
        dataframe_to_html(prevalence_df),
        "</section>",
        "<section class='section'><h2>Stratified Selected-Model Metrics</h2>",
        "<h3>By Source GSE</h3>",
        dataframe_to_html(source_metrics, max_rows=40),
        "<h3>By Tissue</h3>",
        dataframe_to_html(tissue_metrics, max_rows=40),
        "</section>",
        "<section class='section'><h2>Whole-Dataset Predictions</h2>",
        "<h3>Counts by Source GSE</h3>",
        dataframe_to_html(source_counts, max_rows=40),
        "<h3>Counts by Tissue</h3>",
        dataframe_to_html(tissue_counts, max_rows=40),
        "<h3>Reference Criteria Overlap</h3>",
        dataframe_to_html(overlap_df),
        "</section>",
        f"<section class='section'><h2>Figures</h2><div class='figure-grid'>{figure_html}</div></section>",
        "<div class='lightbox' id='lightbox'><button class='lightbox-close' id='lightbox-close' aria-label='Close'>&times;</button><div><img id='lightbox-image' alt='Expanded figure'><div class='lightbox-caption' id='lightbox-caption'></div></div></div>",
        "<script>const lb=document.getElementById('lightbox'),img=document.getElementById('lightbox-image'),cap=document.getElementById('lightbox-caption');document.querySelectorAll('.zoomable-figure').forEach(el=>el.addEventListener('click',()=>{img.src=el.src;cap.textContent=el.dataset.title||el.alt||'';lb.classList.add('open');}));function closeLb(){lb.classList.remove('open');img.removeAttribute('src');}document.getElementById('lightbox-close').addEventListener('click',closeLb);lb.addEventListener('click',e=>{if(e.target===lb)closeLb();});document.addEventListener('keydown',e=>{if(e.key==='Escape')closeLb();});</script>",
        "</main></body></html>",
    ]
    REPORT_HTML.write_text("\n".join(html_parts), encoding="utf-8")

    if "Ground Truth Definition" not in REPORT_HTML.read_text(encoding="utf-8"):
        raise RuntimeError("Rendered HTML is missing the Ground Truth Definition section.")


def render_pdf() -> None:
    chrome_candidates = [
        Path("/usr/bin/google-chrome"),
        Path("/usr/bin/google-chrome-stable"),
        shutil.which("google-chrome"),
        shutil.which("google-chrome-stable"),
    ]
    chrome = next((Path(path) for path in chrome_candidates if path and Path(path).exists()), None)
    if chrome is None:
        raise FileNotFoundError("google-chrome not found for PDF export.")
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


def main() -> None:
    args = parse_args()
    ensure_dirs()
    setup_logging()
    rng = np.random.default_rng(args.seed)
    input_h5ad = args.input_h5ad.resolve()
    if not input_h5ad.exists():
        raise FileNotFoundError(input_h5ad)
    stat_before = input_h5ad.stat()

    logging.info("Reading plus6 metadata and Phase 4 scores from %s", input_h5ad)
    with h5py.File(input_h5ad, "r") as handle:
        required_columns_present(handle)
        n_obs = int(handle["obs"]["_index"].shape[0])
        source = clean_group_values(read_obs_column(handle, "source_gse_id"))
        library_id = normalize_strings(read_obs_column(handle, "library_id"))
        sample_id = normalize_strings(read_obs_column(handle, "sample_id"))
        tissue_column = "tissue_corrected" if "tissue_corrected" in handle["obs"] else "tissue"
        tissue = clean_group_values(read_obs_column(handle, tissue_column))
        has_TRA_TRB_paired = read_bool_obs(handle, "has_TRA_TRB_paired")
        has_TRG_TRD_paired_raw = read_bool_obs(handle, "has_TRG_TRD_paired")
        has_any_ab_tcr = read_bool_obs(handle, "has_any_ab_tcr")
        has_any_gd_tcr_raw = read_bool_obs(handle, "has_any_gd_tcr")
        trg_nonempty = read_nonempty_string_mask(handle, "TRG_cdr3")
        trd_nonempty = read_nonempty_string_mask(handle, "TRD_cdr3")
        trd_score = read_float_obs(handle, "phase4_trd_score")
        trab_score = read_float_obs(handle, "phase4_trab_score")
        trd_minus_trab = read_float_obs(handle, "phase4_trd_minus_trab")
        sorted_gdt = read_bool_obs(handle, "Sorted_gdT") if "Sorted_gdT" in handle["obs"] else np.isin(source, list(SORTED_GDT_SOURCES))
        simple_annotation = (
            clean_group_values(read_obs_column(handle, "simple_annotation_plus6"))
            if "simple_annotation_plus6" in handle["obs"]
            else None
        )
        cluster_values = clean_group_values(read_obs_column(handle, args.cluster_col)) if args.cluster_col in handle["obs"] else None

    logging.info("Building ground-truth labels")
    corrected_has_any_gd_tcr, corrected_has_TRG_TRD_paired, tcr_evidence_audit = build_corrected_tcr_evidence(
        source=source,
        has_TRA_TRB_paired=has_TRA_TRB_paired,
        has_TRG_TRD_paired_raw=has_TRG_TRD_paired_raw,
        has_any_ab_tcr=has_any_ab_tcr,
        has_any_gd_tcr_raw=has_any_gd_tcr_raw,
        trg_nonempty=trg_nonempty,
        trd_nonempty=trd_nonempty,
    )
    is_gdtcr_sequenced_sublibrary, sublibrary_summary = build_sublibrary_labels(
        source=source,
        library_id=library_id,
        sample_id=sample_id,
        corrected_has_any_gd_tcr=corrected_has_any_gd_tcr,
        trg_nonempty=trg_nonempty,
        trd_nonempty=trd_nonempty,
    )
    class_code, conflict_df = make_truth_labels(
        source=source,
        is_gdtcr_sequenced_sublibrary=is_gdtcr_sequenced_sublibrary,
        has_TRA_TRB_paired=has_TRA_TRB_paired,
        corrected_has_TRG_TRD_paired=corrected_has_TRG_TRD_paired,
        has_any_ab_tcr=has_any_ab_tcr,
        corrected_has_any_gd_tcr=corrected_has_any_gd_tcr,
        trd_nonempty=trd_nonempty,
    )
    truth_tables = write_ground_truth_tables(
        class_code, source, tissue, sublibrary_summary, conflict_df, tcr_evidence_audit
    )

    primary_mask = (class_code == 1) | (class_code == 2)
    silver_mask = class_code == 3
    y_true_full = (class_code == 2).astype(np.int8)
    y_true_primary = y_true_full[primary_mask]
    if int((class_code == 2).sum()) == 0 or int((class_code == 1).sum()) == 0:
        raise RuntimeError("Primary positives and negatives must both be nonempty.")
    if int(conflict_df["n_cells"].sum()) != 0:
        raise RuntimeError("Ground-truth rules produced overlapping labels; see ground_truth_conflicts.csv.")
    if bool((silver_mask & primary_mask).any()):
        raise RuntimeError("Silver cells overlapped primary threshold-fitting labels.")

    logging.info(
        "Primary labels: %s gdT_gold, %s abT_gold; silver sensitivity cells: %s",
        int((class_code == 2).sum()),
        int((class_code == 1).sum()),
        int(silver_mask.sum()),
    )

    score_models = [
        ScoreModel("TRD_score_high", "phase4_trd_score", trd_score),
        ScoreModel("TRAB_score_low", "-phase4_trab_score", -trab_score),
        ScoreModel("TRD_minus_TRAB_high", "phase4_trd_minus_trab", trd_minus_trab),
    ]
    logging.info("Selecting conservative high-specificity thresholds on primary gold labels")
    threshold_df, best_model = evaluate_score_models(score_models, y_true_full, primary_mask)
    operating_point_df = build_operating_point_grid(
        score_models,
        y_true_full,
        primary_mask,
        OPERATING_POINT_SPECIFICITIES,
    )

    trd_model = next(model for model in score_models if model.model == "TRD_score_high")
    trab_low_model = next(model for model in score_models if model.model == "TRAB_score_low")
    secondary_df, top_clusters = evaluate_dual_and_cluster_rules(
        y_true=y_true_full,
        primary_mask=primary_mask,
        trd_score=trd_score,
        trab_score=trab_score,
        trd_minus_trab=trd_minus_trab,
        trd_threshold=float(trd_model.threshold),
        trab_raw_cutoff=float(trab_low_model.raw_cutoff),
        cluster_values=cluster_values,
        n_top_clusters=args.n_top_clusters,
    )

    selected_pred = np.asarray(best_model.y_pred_full, dtype=bool)
    selected_count = int(selected_pred.sum())
    if selected_count != int(np.asarray(best_model.score >= float(best_model.threshold), dtype=bool).sum()):
        raise RuntimeError("Whole-dataset prediction count does not match selected cutoff mask.")

    selected_primary_pred = selected_pred[primary_mask].astype(np.int8)
    selected_primary_score = best_model.score[primary_mask]
    silver_recall = float(selected_pred[silver_mask].mean()) if silver_mask.any() else float("nan")

    source_metrics = evaluate_selected_by_group(
        source, "source_gse_id", primary_mask, y_true_full, selected_pred, best_model.score
    )
    tissue_metrics = evaluate_selected_by_group(tissue, "tissue", primary_mask, y_true_full, selected_pred, best_model.score)
    source_metrics.to_csv(TABLE_DIR / "selected_model_metrics_by_source_gse.csv", index=False)
    tissue_metrics.to_csv(TABLE_DIR / "selected_model_metrics_by_tissue.csv", index=False)

    source_counts = selected_counts_by_group(source, "source_gse_id", selected_pred)
    tissue_counts = selected_counts_by_group(tissue, "tissue", selected_pred)
    source_counts.to_csv(TABLE_DIR / "prediction_counts_full_by_source_gse.csv", index=False)
    tissue_counts.to_csv(TABLE_DIR / "prediction_counts_full_by_tissue.csv", index=False)

    prevalence_rows = []
    primary_prevalence = float(y_true_primary.mean())
    for model in score_models:
        metrics = binary_metrics(y_true_primary, model.y_pred_primary, model.score[primary_mask])
        pred_n = int(np.asarray(model.y_pred_full, dtype=bool).sum())
        precision = float(metrics["precision"])
        prevalence_rows.append(
            {
                "model": model.model,
                "full_predicted_positive": pred_n,
                "full_predicted_fraction": pred_n / n_obs,
                "primary_positive_prevalence": primary_prevalence,
                "primary_precision": precision,
                "primary_recall": float(metrics["recall"]),
                "primary_specificity": float(metrics["specificity"]),
                "precision_adjusted_estimated_true_positive_if_precision_transfers": pred_n * precision,
                "precision_adjusted_estimated_false_positive_if_precision_transfers": pred_n * (1.0 - precision),
            }
        )
    prevalence_df = pd.DataFrame(prevalence_rows)
    prevalence_df.to_csv(TABLE_DIR / "prediction_prevalence_adjusted_counts.csv", index=False)

    overlap_df = build_overlap_table(
        selected_pred,
        sorted_gdt,
        corrected_has_TRG_TRD_paired,
        has_any_ab_tcr,
        trd_minus_trab,
        simple_annotation,
    )

    selected_metrics = binary_metrics(y_true_primary, selected_primary_pred, selected_primary_score)
    selected_summary = {
        "input_h5ad": str(input_h5ad),
        "n_obs": n_obs,
        "tissue_column": tissue_column,
        "threshold_policy": (
            f"specificity>={CONSERVATIVE_MIN_SPECIFICITY:g}; "
            f"maximize_F{CONSERVATIVE_F_BETA:g}"
        ),
        "selected_model": best_model.model,
        "selected_score_label": best_model.score_label,
        "youden_threshold_on_score_for_selected_model": best_model.youden_threshold,
        "selected_threshold_on_score": best_model.threshold,
        "selected_raw_cutoff_for_original_score": best_model.raw_cutoff,
        "selected_specificity_constraint": CONSERVATIVE_MIN_SPECIFICITY,
        "selected_f_beta_beta": CONSERVATIVE_F_BETA,
        "selected_full_dataset_predicted_positive": selected_count,
        "selected_full_dataset_predicted_fraction": selected_count / n_obs,
        "selected_primary_metrics": selected_metrics,
        "silver_gdT_recall_under_selected_cutoff": silver_recall,
        "primary_label_counts": {
            "gdT_gold": int((class_code == 2).sum()),
            "abT_gold": int((class_code == 1).sum()),
            "gdT_silver_sensitivity_only": int((class_code == 3).sum()),
        },
        "tcr_evidence_correction_audit": tcr_evidence_audit.to_dict(orient="records"),
        "secondary_cluster_top_clusters": top_clusters,
        "input_h5ad_size_bytes_before": stat_before.st_size,
        "input_h5ad_mtime_ns_before": stat_before.st_mtime_ns,
    }

    logging.info("Rendering figures")
    plot_roc_pr(score_models, primary_mask, y_true_full)
    plot_score_distributions(score_models, class_code, rng)
    plot_confusion(y_true_primary, selected_primary_pred, best_model.model)
    plot_trab_vs_trd_scatter(trd_score, trab_score, class_code, selected_pred, rng, args.scatter_sample)
    plot_umap_predictions(input_h5ad, selected_pred, rng, args.umap_background_sample, args.umap_target_sample)
    plot_count_bars(source_counts, tissue_counts, overlap_df)

    logging.info("Writing summaries and static report")
    write_summary_markdown(
        input_h5ad=input_h5ad,
        n_obs=n_obs,
        best_model=best_model,
        threshold_df=threshold_df,
        operating_point_df=operating_point_df,
        secondary_df=secondary_df,
        truth_tables=truth_tables,
        source_counts=source_counts,
        tissue_counts=tissue_counts,
        overlap_df=overlap_df,
        silver_recall=silver_recall,
        top_clusters=top_clusters,
    )
    copied_assets = copy_static_assets()
    render_html_report(
        input_h5ad=input_h5ad,
        n_obs=n_obs,
        best_model=best_model,
        threshold_df=threshold_df,
        operating_point_df=operating_point_df,
        secondary_df=secondary_df,
        truth_tables=truth_tables,
        source_metrics=source_metrics,
        tissue_metrics=tissue_metrics,
        source_counts=source_counts,
        tissue_counts=tissue_counts,
        overlap_df=overlap_df,
        prevalence_df=prevalence_df,
        silver_recall=silver_recall,
        copied_assets=copied_assets,
    )
    if not args.no_pdf:
        render_pdf()

    stat_after = input_h5ad.stat()
    selected_summary["input_h5ad_size_bytes_after"] = stat_after.st_size
    selected_summary["input_h5ad_mtime_ns_after"] = stat_after.st_mtime_ns
    if (stat_after.st_size != stat_before.st_size) or (stat_after.st_mtime_ns != stat_before.st_mtime_ns):
        raise RuntimeError("Input H5AD file size or mtime changed during read-only evaluation.")

    with SELECTED_JSON.open("w", encoding="utf-8") as handle:
        json.dump(json_ready(selected_summary), handle, indent=2)
    pd.DataFrame([json_ready(selected_summary)]).to_csv(TABLE_DIR / "selected_model_summary.csv", index=False)

    logging.info("Saved HTML report: %s", REPORT_HTML)
    if REPORT_PDF.exists():
        logging.info("Saved PDF report: %s", REPORT_PDF)
    logging.info("Selected model: %s; full predictions: %s", best_model.model, f"{selected_count:,}")


if __name__ == "__main__":
    main()
