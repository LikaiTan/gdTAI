#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Audit TP versus FN validation-positive cells for the promoted gdTAI wrapper.

The script is read-only with respect to the input H5AD. It reuses the accepted
annotation-specific gdTAI prediction array, then compares strict held-out
gdT-positive TP and FN cells across TCR metadata, source composition, targeted
marker expression, and a sparse whole-transcriptome screening summary.
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

try:
    from scipy import stats
except Exception:  # pragma: no cover - scipy is available in the canonical env
    stats = None

import run_gdtai_transcriptome_gate_cascade as base
from run_gdt_prediction_package_evaluation import (
    DEFAULT_INPUT_H5AD,
    clean_group_values,
    dataframe_to_html,
    dataframe_to_markdown,
    normalize_strings,
    read_bool_obs,
    read_float_obs,
    read_nonempty_string_mask,
    read_obs_column,
    read_string_dataset,
)
from run_gdt_gse144469_holdout_tcrgene_classifier import (
    TARGET_SUM,
    TCR_PREFIXES,
    FeatureSpec,
    extract_gene_matrix,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
STATIC_DIR = PROJECT_ROOT / "gdT_prediction"
OUT_PREFIX = "gdtai_annotation_specific_tp_fn_audit"
TABLE_DIR = OUTPUT_ROOT / "tables" / "gdT_prediction" / OUT_PREFIX
FIGURE_DIR = OUTPUT_ROOT / "figures" / "gdT_prediction" / OUT_PREFIX
LOG_DIR = OUTPUT_ROOT / "logs" / "gdT_prediction" / OUT_PREFIX
STATIC_ASSET_DIR = STATIC_DIR / "assets" / OUT_PREFIX
REPORT_MD = LOG_DIR / "gdtai_annotation_specific_tp_fn_audit_report.md"
REPORT_HTML = STATIC_DIR / "gdtai_annotation_specific_tp_fn_audit.html"
REPORT_PDF = STATIC_DIR / "gdtai_annotation_specific_tp_fn_audit.pdf"
SUMMARY_JSON = LOG_DIR / "gdtai_annotation_specific_tp_fn_audit_summary.json"
RUN_LOG = LOG_DIR / "run.log"

CASCADE_TABLE_DIR = OUTPUT_ROOT / "tables" / "gdT_prediction" / "gdtai_annotation_specific_cascade"
DEFAULT_PREDICTION_ARRAYS = CASCADE_TABLE_DIR / "full_prediction_arrays.npz"
DEFAULT_MODEL_PATH = base.CURRENT_MODEL_PATH
RANDOM_SEED = 20260501

T_CELL_MARKERS = ["CD3D", "CD3E", "CD3G", "CD247", "TRAC", "LCK", "LAT", "ZAP70"]
GDT_TCR_MARKERS = ["TRDC", "TRDV1", "TRDV2", "TRGV9", "TRGC1", "TRGC2", "TRDJ1", "TRDJ4"]
AB_TCR_MARKERS = ["TRAC", "TRBC1", "TRBC2", "TRAV1-2", "TRBV20-1"]
CYTOTOXIC_NK_MARKERS = [
    "NKG7",
    "GNLY",
    "PRF1",
    "GZMB",
    "GZMA",
    "KLRD1",
    "KLRF1",
    "NCR1",
    "TYROBP",
    "FCER1G",
    "FCGR3A",
    "CST7",
]
PENALTY_MARKERS = ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT", "CD4"]
LINEAGE_CONTAMINATION_MARKERS = ["MS4A1", "CD79A", "CD14", "LYZ", "EPCAM", "KRT19", "COL1A1", "PECAM1"]
TARGETED_MARKER_GENES = sorted(
    set(T_CELL_MARKERS + GDT_TCR_MARKERS + AB_TCR_MARKERS + CYTOTOXIC_NK_MARKERS + PENALTY_MARKERS + LINEAGE_CONTAMINATION_MARKERS)
)


@dataclass
class ValidationPositiveSet:
    rows: np.ndarray
    component: np.ndarray
    outcome: np.ndarray
    is_fn: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit TP versus FN cells for the promoted annotation-specific gdTAI wrapper.")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--prediction-arrays", type=Path, default=DEFAULT_PREDICTION_ARRAYS)
    parser.add_argument("--model-pkl", type=Path, default=DEFAULT_MODEL_PATH)
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    parser.add_argument("--umap-background-sample", type=int, default=160_000)
    parser.add_argument("--no-pdf", action="store_true")
    parser.add_argument("--skip-transcriptome-screen", action="store_true")
    return parser.parse_args()


def ensure_dirs() -> None:
    for path in [TABLE_DIR, FIGURE_DIR, LOG_DIR, STATIC_DIR, STATIC_ASSET_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def setup_logging() -> None:
    ensure_dirs()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(RUN_LOG, mode="w", encoding="utf-8"), logging.StreamHandler()],
        force=True,
    )


def safe_div(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if denominator else float("nan")


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    vals = p[valid]
    order = np.argsort(vals)
    ranked = vals[order]
    n = ranked.size
    adj_sorted = ranked * n / np.arange(1, n + 1, dtype=float)
    adj_sorted = np.minimum.accumulate(adj_sorted[::-1])[::-1]
    adj_sorted = np.clip(adj_sorted, 0, 1)
    valid_idx = np.flatnonzero(valid)
    out[valid_idx[order]] = adj_sorted
    return out


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, list):
        return [json_ready(v) for v in value]
    if isinstance(value, tuple):
        return [json_ready(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if not np.isfinite(float(value)) else float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if pd.isna(value) if not isinstance(value, (str, bytes)) else False:
        return None
    return value


def read_optional_bool(handle: h5py.File, column: str, n_obs: int) -> np.ndarray:
    if column not in handle["obs"]:
        return np.zeros(n_obs, dtype=bool)
    return read_bool_obs(handle, column)


def read_optional_nonempty(handle: h5py.File, column: str, n_obs: int) -> np.ndarray:
    if column not in handle["obs"]:
        return np.zeros(n_obs, dtype=bool)
    return read_nonempty_string_mask(handle, column)


def read_optional_group(handle: h5py.File, column: str, n_obs: int, *, default: str = "unknown") -> np.ndarray:
    if column not in handle["obs"]:
        return np.full(n_obs, default, dtype=object)
    return clean_group_values(read_obs_column(handle, column))


def read_optional_string(handle: h5py.File, column: str, n_obs: int) -> np.ndarray:
    if column not in handle["obs"]:
        return np.full(n_obs, "", dtype=object)
    return normalize_strings(read_obs_column(handle, column))


def sublibrary_ids(library_id: np.ndarray, sample_id: np.ndarray) -> np.ndarray:
    library = pd.Series(library_id, copy=False).astype("string").fillna("").str.strip()
    sample = pd.Series(sample_id, copy=False).astype("string").fillna("").str.strip()
    out = library.mask(library == "", sample)
    out = out.mask(out == "", "unknown")
    return out.astype(object).to_numpy()


def validation_positive_set(obs: Any, pred: np.ndarray) -> ValidationPositiveSet:
    primary = (obs.class_code == 1) | (obs.class_code == 2)
    validation_primary_idx = np.flatnonzero(primary & (obs.source == base.HOLDOUT_SOURCE)).astype(np.int64)
    validation_gdt2020_idx = np.flatnonzero(
        primary
        & (obs.source == base.GDT2020_SOURCE)
        & (obs.class_code == 2)
        & (pd.Series(obs.tissue, copy=False).astype(str).str.lower().to_numpy() == base.GDT2020_HOLDOUT_TISSUE.lower())
    ).astype(np.int64)
    validation_extra_trab_idx = np.flatnonzero(
        (obs.source == base.EXTRA_TRAB_HOLDOUT_SOURCE)
        & obs.has_TRA_TRB_paired
        & (~obs.corrected_has_any_gd_tcr)
    ).astype(np.int64)
    validation_idx = np.unique(np.concatenate([validation_primary_idx, validation_gdt2020_idx, validation_extra_trab_idx])).astype(np.int64)
    positive_rows = validation_idx[obs.class_code[validation_idx] == 2]
    if positive_rows.size == 0:
        raise RuntimeError("No validation gdT-positive cells were found.")

    component = np.full(positive_rows.size, "validation_gdT_positive_other", dtype=object)
    source = obs.source[positive_rows]
    tissue = pd.Series(obs.tissue[positive_rows], copy=False).astype(str).str.lower().to_numpy()
    component[source == base.HOLDOUT_SOURCE] = f"validation_primary_{base.HOLDOUT_SOURCE}_gdT_gold"
    component[
        (source == base.GDT2020_SOURCE) & (tissue == base.GDT2020_HOLDOUT_TISSUE.lower())
    ] = f"validation_gdT_{base.GDT2020_SOURCE}_{base.GDT2020_HOLDOUT_TISSUE.replace(' ', '_')}"
    outcome = np.where(pred[positive_rows], "TP", "FN")
    is_fn = outcome == "FN"
    if int((outcome == "TP").sum()) == 0 or int(is_fn.sum()) == 0:
        raise RuntimeError("The validation-positive audit requires both TP and FN cells.")
    return ValidationPositiveSet(rows=positive_rows, component=component, outcome=outcome, is_fn=is_fn)


def tcr_family(gene: str) -> str:
    upper = gene.upper()
    if upper.startswith(("TRAV", "TRAJ", "TRAC")):
        return "TRA"
    if upper.startswith(("TRBV", "TRBJ", "TRBC")):
        return "TRB"
    if upper.startswith(("TRGV", "TRGJ", "TRGC")):
        return "TRG"
    if upper.startswith(("TRDV", "TRDD", "TRDJ", "TRDC")):
        return "TRD"
    return "non_TCR"


def build_feature_spec(handle: h5py.File, genes: list[str], model_genes: set[str]) -> tuple[FeatureSpec, list[str]]:
    var_names = pd.Index(read_string_dataset(handle["var"]["_index"]), dtype="string").astype(str)
    gene_to_idx = {gene: idx for idx, gene in enumerate(var_names)}
    selected: list[str] = []
    missing: list[str] = []
    for gene in genes:
        if gene in gene_to_idx:
            if gene not in selected:
                selected.append(gene)
        else:
            missing.append(gene)
    spec = FeatureSpec(
        feature_names=[f"{gene}_log1p_cp10k" for gene in selected],
        gene_names=selected,
        gene_indices=np.asarray([gene_to_idx[gene] for gene in selected], dtype=np.int32),
        tcr_gene_names=[gene for gene in selected if gene.upper().startswith(TCR_PREFIXES)],
        penalty_gene_names=[gene for gene in selected if gene in PENALTY_MARKERS],
        family_by_feature=[tcr_family(gene) for gene in selected],
        cd3_cols=[i for i, gene in enumerate(selected) if gene in {"CD3D", "CD3E", "CD3G"}],
        cd4_col=selected.index("CD4") if "CD4" in selected else None,
        foxp3_col=selected.index("FOXP3") if "FOXP3" in selected else None,
    )
    manifest = pd.DataFrame(
        {
            "feature_index": np.arange(len(selected), dtype=int),
            "gene": selected,
            "feature": spec.feature_names,
            "family": spec.family_by_feature,
            "used_by_promoted_base_model": [gene in model_genes for gene in selected],
            "targeted_marker_panel": [gene in set(TARGETED_MARKER_GENES) for gene in selected],
        }
    )
    manifest.to_csv(TABLE_DIR / "targeted_feature_manifest.csv", index=False)
    pd.DataFrame({"missing_gene": sorted(set(missing))}).to_csv(TABLE_DIR / "targeted_missing_genes.csv", index=False)
    return spec, missing


def mean_available(x: np.ndarray, spec: FeatureSpec, genes: list[str]) -> np.ndarray:
    cols = [spec.gene_names.index(gene) for gene in genes if gene in spec.gene_names]
    if not cols:
        return np.zeros(x.shape[0], dtype=np.float32)
    return x[:, cols].mean(axis=1).astype(np.float32, copy=False)


def add_signature_scores(cell_df: pd.DataFrame, x: np.ndarray, spec: FeatureSpec) -> None:
    signatures = {
        "expr_t_cell_marker_mean": T_CELL_MARKERS,
        "expr_gdt_tcr_marker_mean": GDT_TCR_MARKERS,
        "expr_ab_tcr_marker_mean": AB_TCR_MARKERS,
        "expr_cytotoxic_nk_marker_mean": CYTOTOXIC_NK_MARKERS,
        "expr_penalty_marker_mean": PENALTY_MARKERS,
        "expr_lineage_contamination_marker_mean": LINEAGE_CONTAMINATION_MARKERS,
    }
    for name, genes in signatures.items():
        cell_df[name] = mean_available(x, spec, genes)
    for gene in ["TRDC", "TRDV1", "TRDV2", "TRGV9", "TRGC1", "TRGC2", "CD3D", "CD3E", "CD3G", "CD4", "FOXP3", "NKG7", "GNLY", "KLRD1"]:
        if gene in spec.gene_names:
            cell_df[f"expr_{gene}"] = x[:, spec.gene_names.index(gene)]


def categorical_summary(cell_df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for column in columns:
        if column not in cell_df:
            continue
        values = cell_df[column].astype(str).fillna("unknown")
        for level, idx in values.groupby(values, sort=True).groups.items():
            mask = cell_df.index.isin(idx)
            n_total = int(mask.sum())
            n_tp = int(((cell_df["outcome"] == "TP") & mask).sum())
            n_fn = int(((cell_df["outcome"] == "FN") & mask).sum())
            rows.append(
                {
                    "variable": column,
                    "level": level,
                    "n_cells": n_total,
                    "tp": n_tp,
                    "fn": n_fn,
                    "fn_rate": safe_div(n_fn, n_total),
                    "fraction_of_all_fn": safe_div(n_fn, int((cell_df["outcome"] == "FN").sum())),
                    "fraction_of_all_tp": safe_div(n_tp, int((cell_df["outcome"] == "TP").sum())),
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["variable", "fn_rate", "n_cells"], ascending=[True, False, False])
    return out


def bool_summary(cell_df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for column in columns:
        if column not in cell_df:
            continue
        present = cell_df[column].astype(bool).to_numpy()
        is_fn = (cell_df["outcome"].to_numpy() == "FN")
        is_tp = ~is_fn
        fn_present = int((present & is_fn).sum())
        tp_present = int((present & is_tp).sum())
        fn_absent = int(((~present) & is_fn).sum())
        tp_absent = int(((~present) & is_tp).sum())
        fn_frac = safe_div(fn_present, int(is_fn.sum()))
        tp_frac = safe_div(tp_present, int(is_tp.sum()))
        odds_ratio = float("nan")
        p_value = float("nan")
        if stats is not None:
            try:
                odds_ratio, p_value = stats.fisher_exact([[fn_present, fn_absent], [tp_present, tp_absent]])
            except Exception:
                pass
        rows.append(
            {
                "variable": column,
                "fn_present": fn_present,
                "fn_absent": fn_absent,
                "tp_present": tp_present,
                "tp_absent": tp_absent,
                "fn_fraction_present": fn_frac,
                "tp_fraction_present": tp_frac,
                "fn_minus_tp_fraction": fn_frac - tp_frac,
                "fn_vs_tp_odds_ratio": odds_ratio,
                "p_value": p_value,
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_value"] = bh_adjust(out["p_value"].to_numpy(dtype=float))
        out = out.sort_values(["q_value", "fn_minus_tp_fraction"], ascending=[True, False])
    return out


def numeric_summary(cell_df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    is_fn = cell_df["outcome"].to_numpy() == "FN"
    for column in columns:
        if column not in cell_df:
            continue
        vals = pd.to_numeric(cell_df[column], errors="coerce").to_numpy(dtype=float)
        fn_vals = vals[is_fn & np.isfinite(vals)]
        tp_vals = vals[(~is_fn) & np.isfinite(vals)]
        p_value = float("nan")
        if stats is not None and fn_vals.size and tp_vals.size:
            try:
                p_value = float(stats.mannwhitneyu(fn_vals, tp_vals, alternative="two-sided").pvalue)
            except Exception:
                p_value = float("nan")
        rows.append(
            {
                "variable": column,
                "tp_n": int(tp_vals.size),
                "fn_n": int(fn_vals.size),
                "tp_mean": float(np.mean(tp_vals)) if tp_vals.size else float("nan"),
                "fn_mean": float(np.mean(fn_vals)) if fn_vals.size else float("nan"),
                "fn_minus_tp_mean": (float(np.mean(fn_vals)) - float(np.mean(tp_vals))) if fn_vals.size and tp_vals.size else float("nan"),
                "tp_median": float(np.median(tp_vals)) if tp_vals.size else float("nan"),
                "fn_median": float(np.median(fn_vals)) if fn_vals.size else float("nan"),
                "fn_minus_tp_median": (float(np.median(fn_vals)) - float(np.median(tp_vals))) if fn_vals.size and tp_vals.size else float("nan"),
                "tp_q25": float(np.quantile(tp_vals, 0.25)) if tp_vals.size else float("nan"),
                "tp_q75": float(np.quantile(tp_vals, 0.75)) if tp_vals.size else float("nan"),
                "fn_q25": float(np.quantile(fn_vals, 0.25)) if fn_vals.size else float("nan"),
                "fn_q75": float(np.quantile(fn_vals, 0.75)) if fn_vals.size else float("nan"),
                "p_value": p_value,
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_value"] = bh_adjust(out["p_value"].to_numpy(dtype=float))
        out = out.sort_values(["q_value", "fn_minus_tp_median"], ascending=[True, False])
    return out


def summarize_targeted_expression(x: np.ndarray, spec: FeatureSpec, outcome: np.ndarray) -> pd.DataFrame:
    is_fn = outcome == "FN"
    rows: list[dict[str, Any]] = []
    for col, gene in enumerate(spec.gene_names):
        fn_vals = x[is_fn, col]
        tp_vals = x[~is_fn, col]
        p_value = float("nan")
        if stats is not None and fn_vals.size and tp_vals.size:
            try:
                p_value = float(stats.mannwhitneyu(fn_vals, tp_vals, alternative="two-sided").pvalue)
            except Exception:
                p_value = float("nan")
        tp_mean = float(tp_vals.mean()) if tp_vals.size else float("nan")
        fn_mean = float(fn_vals.mean()) if fn_vals.size else float("nan")
        rows.append(
            {
                "gene": gene,
                "family": tcr_family(gene),
                "tp_mean_log1p_cp10k": tp_mean,
                "fn_mean_log1p_cp10k": fn_mean,
                "fn_minus_tp_mean": fn_mean - tp_mean,
                "tp_pct_expressing": safe_div(float((tp_vals > 0).sum()), float(tp_vals.size)),
                "fn_pct_expressing": safe_div(float((fn_vals > 0).sum()), float(fn_vals.size)),
                "fn_minus_tp_pct_expressing": safe_div(float((fn_vals > 0).sum()), float(fn_vals.size)) - safe_div(float((tp_vals > 0).sum()), float(tp_vals.size)),
                "p_value": p_value,
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_value"] = bh_adjust(out["p_value"].to_numpy(dtype=float))
        out = out.sort_values(["q_value", "fn_minus_tp_mean"], ascending=[True, False])
    return out


def screen_transcriptome(handle: h5py.File, rows: np.ndarray, outcome: np.ndarray) -> pd.DataFrame:
    var_names = pd.Index(read_string_dataset(handle["var"]["_index"]), dtype="string").astype(str)
    n_genes = len(var_names)
    is_fn = outcome == "FN"
    n_fn = int(is_fn.sum())
    n_tp = int((~is_fn).sum())
    sums_fn = np.zeros(n_genes, dtype=np.float64)
    sums_tp = np.zeros(n_genes, dtype=np.float64)
    sumsq_fn = np.zeros(n_genes, dtype=np.float64)
    sumsq_tp = np.zeros(n_genes, dtype=np.float64)
    nnz_fn = np.zeros(n_genes, dtype=np.int32)
    nnz_tp = np.zeros(n_genes, dtype=np.int32)
    x_group = handle["X"]
    indptr = x_group["indptr"]
    indices = x_group["indices"]
    data = x_group["data"]
    for local_pos, obs_idx in enumerate(rows):
        start = int(indptr[obs_idx])
        end = int(indptr[obs_idx + 1])
        if end <= start:
            continue
        idx = indices[start:end].astype(np.int64, copy=False)
        vals = data[start:end].astype(np.float32, copy=False)
        row_sum = float(np.sum(vals, dtype=np.float64))
        if row_sum <= 0:
            continue
        norm = np.log1p(vals * (TARGET_SUM / row_sum)).astype(np.float64, copy=False)
        if is_fn[local_pos]:
            np.add.at(sums_fn, idx, norm)
            np.add.at(sumsq_fn, idx, norm * norm)
            np.add.at(nnz_fn, idx, 1)
        else:
            np.add.at(sums_tp, idx, norm)
            np.add.at(sumsq_tp, idx, norm * norm)
            np.add.at(nnz_tp, idx, 1)
        if local_pos and local_pos % 2500 == 0:
            logging.info("Screened transcriptome for %s / %s validation positives", f"{local_pos:,}", f"{rows.size:,}")

    mean_fn = sums_fn / max(n_fn, 1)
    mean_tp = sums_tp / max(n_tp, 1)
    var_fn = np.maximum((sumsq_fn - n_fn * mean_fn * mean_fn) / max(n_fn - 1, 1), 0)
    var_tp = np.maximum((sumsq_tp - n_tp * mean_tp * mean_tp) / max(n_tp - 1, 1), 0)
    p_value = np.full(n_genes, np.nan, dtype=float)
    if stats is not None:
        _, p_value = stats.ttest_ind_from_stats(
            mean1=mean_fn,
            std1=np.sqrt(var_fn),
            nobs1=n_fn,
            mean2=mean_tp,
            std2=np.sqrt(var_tp),
            nobs2=n_tp,
            equal_var=False,
        )
    upper = var_names.str.upper()
    is_tcr = np.asarray([gene.startswith(TCR_PREFIXES) for gene in upper], dtype=bool)
    out = pd.DataFrame(
        {
            "gene": var_names.to_numpy(dtype=object),
            "is_tcr_gene": is_tcr,
            "tp_mean_log1p_cp10k": mean_tp,
            "fn_mean_log1p_cp10k": mean_fn,
            "fn_minus_tp_mean": mean_fn - mean_tp,
            "tp_pct_expressing": nnz_tp / max(n_tp, 1),
            "fn_pct_expressing": nnz_fn / max(n_fn, 1),
            "fn_minus_tp_pct_expressing": (nnz_fn / max(n_fn, 1)) - (nnz_tp / max(n_tp, 1)),
            "log2_ratio_mean_log1p_cp10k": np.log2((mean_fn + 0.01) / (mean_tp + 0.01)),
            "p_value": p_value,
        }
    )
    out["q_value"] = bh_adjust(out["p_value"].to_numpy(dtype=float))
    out["abs_mean_delta"] = out["fn_minus_tp_mean"].abs()
    return out.sort_values(["q_value", "abs_mean_delta"], ascending=[True, False])


def top_nontcr_tables(screen_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if screen_df.empty:
        empty = pd.DataFrame()
        return empty, empty, empty
    nontcr = screen_df.loc[~screen_df["is_tcr_gene"]].copy()
    tcr = screen_df.loc[screen_df["is_tcr_gene"]].copy()
    fn_enriched = nontcr.sort_values(["fn_minus_tp_mean", "q_value"], ascending=[False, True]).head(40)
    tp_enriched = nontcr.sort_values(["fn_minus_tp_mean", "q_value"], ascending=[True, True]).head(40)
    tcr_top = tcr.sort_values(["q_value", "abs_mean_delta"], ascending=[True, False]).head(80)
    return fn_enriched, tp_enriched, tcr_top


def plot_source_fn_rate(source_summary: pd.DataFrame) -> Path:
    path = FIGURE_DIR / "source_fn_rate_validation_positive.png"
    df = source_summary.sort_values("n_cells", ascending=False).head(20)
    fig, ax = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
    colors = ["#dc2626" if value > 0.2 else "#2563eb" for value in df["fn_rate"]]
    ax.barh(df["source_gse_id"].astype(str), df["fn_rate"], color=colors)
    ax.invert_yaxis()
    ax.set_xlabel("FN fraction among validation gdT positives")
    ax.set_ylabel("Source")
    ax.set_title("Source composition of missed validation gdT cells")
    for i, row in enumerate(df.itertuples(index=False)):
        ax.text(row.fn_rate + 0.006, i, f"{int(row.fn)}/{int(row.n_cells)}", va="center", fontsize=8)
    ax.set_xlim(0, max(0.05, min(1.0, float(df["fn_rate"].max()) * 1.25)))
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_bool_tcr_summary(summary: pd.DataFrame) -> Path:
    path = FIGURE_DIR / "tcr_metadata_fraction_tp_fn.png"
    if summary.empty:
        return path
    df = summary.sort_values("fn_minus_tp_fraction", key=lambda s: s.abs(), ascending=False).head(18).iloc[::-1]
    y = np.arange(df.shape[0])
    fig, ax = plt.subplots(figsize=(7.8, 5.2), constrained_layout=True)
    ax.barh(y - 0.18, df["tp_fraction_present"], height=0.34, color="#2563eb", label="TP")
    ax.barh(y + 0.18, df["fn_fraction_present"], height=0.34, color="#dc2626", label="FN")
    ax.set_yticks(y, df["variable"])
    ax.set_xlabel("Fraction present")
    ax.set_title("TCR metadata differences")
    ax.legend(frameon=False, loc="lower right")
    ax.set_xlim(0, 1)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_score_distributions(cell_df: pd.DataFrame) -> Path:
    path = FIGURE_DIR / "score_distributions_tp_fn.png"
    cols = [
        ("current_score", "gdTAI score"),
        ("phase4_trd_score", "TRD score"),
        ("phase4_trab_score", "TRAB score"),
        ("phase4_trd_minus_trab", "TRD minus TRAB"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.2), constrained_layout=True)
    for ax, (col, title) in zip(axes.flat, cols):
        data = [cell_df.loc[cell_df["outcome"] == outcome, col].dropna().to_numpy(dtype=float) for outcome in ["TP", "FN"]]
        parts = ax.violinplot(data, positions=[0, 1], showmeans=False, showmedians=True, widths=0.75)
        for body, color in zip(parts["bodies"], ["#2563eb", "#dc2626"]):
            body.set_facecolor(color)
            body.set_alpha(0.35)
            body.set_edgecolor("none")
        parts["cmedians"].set_color("#111827")
        ax.set_xticks([0, 1], ["TP", "FN"])
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.18)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_trd_vs_trab(cell_df: pd.DataFrame) -> Path:
    path = FIGURE_DIR / "trd_vs_trab_tp_fn_validation_positive.png"
    fig, ax = plt.subplots(figsize=(6.4, 5.4), constrained_layout=True)
    for outcome, color, alpha, size in [("TP", "#2563eb", 0.25, 6), ("FN", "#dc2626", 0.55, 9)]:
        df = cell_df.loc[cell_df["outcome"] == outcome]
        ax.scatter(df["phase4_trab_score"], df["phase4_trd_score"], s=size, c=color, alpha=alpha, linewidths=0, rasterized=True, label=outcome)
    ax.set_xlabel("TRAB score")
    ax.set_ylabel("TRD score")
    ax.set_title("Validation gdT positives: TRD versus TRAB score")
    ax.legend(frameon=False, markerscale=2)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_marker_heatmap(summary: pd.DataFrame) -> Path:
    path = FIGURE_DIR / "targeted_marker_expression_tp_fn.png"
    if summary.empty:
        return path
    markers = [gene for gene in TARGETED_MARKER_GENES if gene in set(summary["gene"])]
    df = summary.set_index("gene").loc[markers].copy()
    df = df.reindex(df["fn_minus_tp_mean"].abs().sort_values(ascending=False).head(30).index)
    mat = df[["tp_mean_log1p_cp10k", "fn_mean_log1p_cp10k"]].to_numpy()
    fig, ax = plt.subplots(figsize=(5.2, max(4.0, 0.21 * df.shape[0] + 1.5)), constrained_layout=True)
    im = ax.imshow(mat, aspect="auto", cmap="viridis")
    ax.set_xticks([0, 1], ["TP", "FN"])
    ax.set_yticks(np.arange(df.shape[0]), df.index.tolist(), fontsize=8)
    ax.set_title("Targeted marker expression")
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("mean log1p CP10K")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_top_transcriptome(screen_df: pd.DataFrame) -> Path:
    path = FIGURE_DIR / "top_nontcr_transcriptome_differences.png"
    if screen_df.empty:
        return path
    nontcr = screen_df.loc[~screen_df["is_tcr_gene"]].copy()
    fn_top = nontcr.sort_values(["fn_minus_tp_mean", "q_value"], ascending=[False, True]).head(12)
    tp_top = nontcr.sort_values(["fn_minus_tp_mean", "q_value"], ascending=[True, True]).head(12)
    df = pd.concat([fn_top, tp_top], ignore_index=True)
    if df.empty:
        return path
    df = df.drop_duplicates("gene").set_index("gene")
    mat = df[["tp_mean_log1p_cp10k", "fn_mean_log1p_cp10k"]].to_numpy()
    fig, ax = plt.subplots(figsize=(5.4, max(4.0, 0.24 * df.shape[0] + 1.5)), constrained_layout=True)
    im = ax.imshow(mat, aspect="auto", cmap="magma")
    ax.set_xticks([0, 1], ["TP", "FN"])
    ax.set_yticks(np.arange(df.shape[0]), df.index.tolist(), fontsize=8)
    ax.set_title("Top non-TCR transcriptome differences")
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("mean log1p CP10K")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_umap(handle: h5py.File, cell_df: pd.DataFrame, rows: np.ndarray, rng: np.random.Generator, background_sample: int) -> Path:
    path = FIGURE_DIR / "umap_validation_positive_tp_fn.png"
    if "obsm" not in handle or "X_umap" not in handle["obsm"]:
        return path
    n_obs = int(handle["obs"]["_index"].shape[0])
    bg_size = min(background_sample, n_obs)
    bg_idx = rng.choice(np.arange(n_obs, dtype=np.int64), size=bg_size, replace=False)
    # Read the two-column UMAP contiguously; random HDF5 point selection is very slow here.
    umap_all = np.asarray(handle["obsm"]["X_umap"][:], dtype=np.float32)
    bg = umap_all[bg_idx, :]
    pos = umap_all[rows, :]
    fig, ax = plt.subplots(figsize=(6.7, 5.8), constrained_layout=True)
    ax.scatter(bg[:, 0], bg[:, 1], s=0.25, c="#d1d5db", alpha=0.15, linewidths=0, rasterized=True, label="atlas background")
    for outcome, color, size, alpha in [("TP", "#2563eb", 5, 0.25), ("FN", "#dc2626", 8, 0.65)]:
        mask = cell_df["outcome"].to_numpy() == outcome
        ax.scatter(pos[mask, 0], pos[mask, 1], s=size, c=color, alpha=alpha, linewidths=0, rasterized=True, label=outcome)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("Validation gdT positives on atlas UMAP")
    ax.legend(frameon=False, markerscale=2, loc="best")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def compact_numeric(df: pd.DataFrame, digits: int = 4) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.copy()
    for col in out.columns:
        if pd.api.types.is_float_dtype(out[col]):
            out[col] = out[col].map(lambda value: "" if pd.isna(value) else round(float(value), digits))
    return out


def copy_assets(paths: list[Path]) -> list[str]:
    copied: list[str] = []
    for path in paths:
        if path.exists():
            target = STATIC_ASSET_DIR / path.name
            shutil.copyfile(path, target)
            copied.append(path.name)
    return copied


def render_pdf(no_pdf: bool) -> None:
    if no_pdf:
        return
    candidates = [Path("/usr/bin/google-chrome"), Path("/usr/bin/google-chrome-stable"), shutil.which("google-chrome"), shutil.which("google-chrome-stable")]
    chrome = next((Path(item) for item in candidates if item and Path(item).exists()), None)
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


def summarize_group(cell_df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    out = (
        cell_df.groupby(group_cols + ["outcome"], dropna=False)
        .size()
        .unstack(fill_value=0)
        .reset_index()
        .rename_axis(None, axis=1)
    )
    for col in ["TP", "FN"]:
        if col not in out:
            out[col] = 0
    out["n_cells"] = out["TP"] + out["FN"]
    out["tp"] = out["TP"].astype(int)
    out["fn"] = out["FN"].astype(int)
    out["fn_rate"] = out["fn"] / out["n_cells"].replace(0, np.nan)
    out = out.drop(columns=["TP", "FN"])
    return out.sort_values(["fn_rate", "n_cells"], ascending=[False, False])


def build_summary_text(
    cell_df: pd.DataFrame,
    source_summary: pd.DataFrame,
    bool_df: pd.DataFrame,
    numeric_df: pd.DataFrame,
    targeted_df: pd.DataFrame,
    screen_df: pd.DataFrame,
) -> list[str]:
    n_tp = int((cell_df["outcome"] == "TP").sum())
    n_fn = int((cell_df["outcome"] == "FN").sum())
    lines = [
        f"Validation-positive audit size: {n_tp + n_fn:,} gdT_gold cells, with {n_tp:,} TP and {n_fn:,} FN under the promoted annotation-specific wrapper.",
    ]
    if not source_summary.empty:
        top_source = source_summary.sort_values(["fn_rate", "n_cells"], ascending=[False, False]).iloc[0]
        lines.append(
            f"Highest source-level FN fraction: {top_source['source_gse_id']} ({int(top_source['fn'])}/{int(top_source['n_cells'])}, {top_source['fn_rate']:.1%})."
        )
    if not bool_df.empty:
        tcr_bool = bool_df.loc[~bool_df["variable"].isin(["current_pred", "original_trd_minus_trab_pred"])].copy()
        if not tcr_bool.empty:
            top_bool = tcr_bool.sort_values("fn_minus_tp_fraction", key=lambda s: s.abs(), ascending=False).iloc[0]
            direction = "more common in FN" if top_bool["fn_minus_tp_fraction"] > 0 else "less common in FN"
            lines.append(
                f"Strongest TCR metadata shift: {top_bool['variable']} is {direction} by {abs(top_bool['fn_minus_tp_fraction']):.1%}."
            )
    signature_rows = numeric_df[numeric_df["variable"].str.startswith("expr_", na=False)] if not numeric_df.empty else pd.DataFrame()
    if not signature_rows.empty:
        top_sig = signature_rows.sort_values("fn_minus_tp_mean", key=lambda s: s.abs(), ascending=False).iloc[0]
        direction = "higher in FN" if top_sig["fn_minus_tp_mean"] > 0 else "lower in FN"
        lines.append(
            f"Largest targeted transcriptome signature shift: {top_sig['variable']} is {direction} (mean delta {top_sig['fn_minus_tp_mean']:.3f})."
        )
    if not targeted_df.empty:
        top_tcr = targeted_df[targeted_df["family"].isin(["TRD", "TRG", "TRA", "TRB"])].sort_values(
            "fn_minus_tp_mean", key=lambda s: s.abs(), ascending=False
        )
        if not top_tcr.empty:
            row = top_tcr.iloc[0]
            direction = "higher in FN" if row["fn_minus_tp_mean"] > 0 else "lower in FN"
            lines.append(f"Largest targeted TCR-gene expression shift: {row['gene']} is {direction} (mean delta {row['fn_minus_tp_mean']:.3f}).")
    if not screen_df.empty:
        nontcr = screen_df.loc[~screen_df["is_tcr_gene"]].copy()
        if not nontcr.empty:
            row = nontcr.sort_values("fn_minus_tp_mean", key=lambda s: s.abs(), ascending=False).iloc[0]
            direction = "higher in FN" if row["fn_minus_tp_mean"] > 0 else "lower in FN"
            lines.append(f"Top non-TCR transcriptome screen difference by mean delta: {row['gene']} is {direction} (mean delta {row['fn_minus_tp_mean']:.3f}).")
    return lines


def write_reports(
    input_h5ad: Path,
    prediction_arrays: Path,
    model_path: Path,
    cell_df: pd.DataFrame,
    source_summary: pd.DataFrame,
    tissue_summary: pd.DataFrame,
    annotation_summary: pd.DataFrame,
    bool_df: pd.DataFrame,
    numeric_df: pd.DataFrame,
    targeted_df: pd.DataFrame,
    tcr_top_df: pd.DataFrame,
    fn_top_df: pd.DataFrame,
    tp_top_df: pd.DataFrame,
    summary_lines: list[str],
    copied_assets: list[str],
    h5ad_mtime_unchanged: bool,
) -> None:
    source_display = compact_numeric(source_summary.head(20))
    tissue_display = compact_numeric(tissue_summary.head(20))
    annotation_display = compact_numeric(annotation_summary.head(30))
    tcr_display = compact_numeric(bool_df.head(30))
    numeric_display = compact_numeric(numeric_df.head(35))
    targeted_display = compact_numeric(targeted_df.head(35))
    tcr_gene_display = compact_numeric(tcr_top_df.head(35))
    fn_gene_display = compact_numeric(fn_top_df.head(25))
    tp_gene_display = compact_numeric(tp_top_df.head(25))

    md_lines = [
        "# gdTAI TP versus FN Audit",
        "",
        "## Definition",
        "",
        "- TP: strict held-out `gdT_gold` validation cells predicted positive by the promoted annotation-specific gdTAI wrapper.",
        "- FN: strict held-out `gdT_gold` validation cells missed by the promoted annotation-specific gdTAI wrapper.",
        f"- Positive validation cells come from `{base.HOLDOUT_SOURCE}` primary gold positives and `{base.GDT2020_SOURCE}` `{base.GDT2020_HOLDOUT_TISSUE}` gdT positives.",
        "- The audit does not include validation negatives except indirectly through the already-promoted model decision.",
        "- TCR metadata follows the harmonized project TCR fields; transcriptome values are log1p CP10K extracted from count-space `X`.",
        "- The input H5AD is opened read-only and was not modified.",
        "",
        "## Summary",
        "",
        *[f"- {line}" for line in summary_lines],
        f"- H5AD mtime unchanged after audit: `{h5ad_mtime_unchanged}`.",
        "",
        "## Source Differences",
        "",
        dataframe_to_markdown(source_display),
        "",
        "## Tissue Differences",
        "",
        dataframe_to_markdown(tissue_display),
        "",
        "## Annotation Differences",
        "",
        dataframe_to_markdown(annotation_display),
        "",
        "## TCR Metadata Differences",
        "",
        dataframe_to_markdown(tcr_display),
        "",
        "## Numeric Score And Marker Differences",
        "",
        dataframe_to_markdown(numeric_display),
        "",
        "## Targeted TCR And Marker Expression",
        "",
        dataframe_to_markdown(targeted_display),
        "",
        "## Top TCR Gene Expression Differences",
        "",
        dataframe_to_markdown(tcr_gene_display),
        "",
        "## FN-Enriched Non-TCR Transcriptome Screen",
        "",
        dataframe_to_markdown(fn_gene_display),
        "",
        "## TP-Enriched Non-TCR Transcriptome Screen",
        "",
        dataframe_to_markdown(tp_gene_display),
        "",
        "## Files",
        "",
        f"- Input H5AD: `{input_h5ad}`",
        f"- Prediction arrays: `{prediction_arrays}`",
        f"- Base model: `{model_path}`",
        f"- Tables: `{TABLE_DIR}`",
        f"- Figures: `{FIGURE_DIR}`",
        f"- HTML: `{REPORT_HTML}`",
        f"- PDF: `{REPORT_PDF}`",
    ]
    REPORT_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    figure_html = []
    for asset in copied_assets:
        title = asset.rsplit(".", 1)[0].replace("_", " ")
        figure_html.append(f"<section class='figure'><h3>{html.escape(title)}</h3><img src='assets/{OUT_PREFIX}/{html.escape(asset)}'></section>")
    summary_html = "".join(f"<li>{html.escape(line)}</li>" for line in summary_lines)
    css = """
    body{font-family:Arial,Helvetica,sans-serif;color:#111827;background:#f8fafc;margin:0}
    main{max-width:1180px;margin:0 auto;padding:28px 28px 64px;background:#fff}
    h1{font-size:30px;margin:0 0 8px} h2{font-size:20px;margin-top:28px;border-top:1px solid #e5e7eb;padding-top:18px}
    h3{font-size:15px;margin:14px 0 8px}.note{color:#475569}.ok{font-weight:700;color:#0f766e}
    table{border-collapse:collapse;width:100%;font-size:11px;margin:10px 0 18px;table-layout:auto}
    th,td{border:1px solid #e5e7eb;padding:5px 7px;text-align:left;vertical-align:top} th{background:#f1f5f9}
    .table-wrap{overflow-x:auto;margin-bottom:12px} img{max-width:100%;height:auto;border:1px solid #e5e7eb;background:#fff}
    code{background:#f1f5f9;padding:1px 4px;border-radius:4px}
    @media print{main{max-width:none;padding:11mm} table{font-size:8px} h1{font-size:22px} h2{font-size:15px}.table-wrap{overflow:visible}}
    """
    html_doc = f"""<!doctype html><html><head><meta charset='utf-8'><title>gdTAI TP versus FN Audit</title><style>{css}</style></head><body><main>
    <h1>gdTAI TP versus FN Audit</h1>
    <p class='note'>TP and FN are defined only among strict held-out gdT_gold validation positives. No H5AD was modified.</p>
    <h2>Definition</h2>
    <p>TP are validation gdT_gold cells predicted positive by the promoted annotation-specific gdTAI wrapper; FN are validation gdT_gold cells missed by that wrapper. Positive validation cells come from <code>{html.escape(base.HOLDOUT_SOURCE)}</code> and <code>{html.escape(base.GDT2020_SOURCE)} {html.escape(base.GDT2020_HOLDOUT_TISSUE)}</code>. TCR metadata uses harmonized project fields, and transcriptome values are log1p CP10K from count-space X.</p>
    <h2>Summary</h2><ul>{summary_html}<li>H5AD mtime unchanged after audit: <code>{str(h5ad_mtime_unchanged)}</code></li></ul>
    <h2>Source Differences</h2><div class='table-wrap'>{dataframe_to_html(source_display)}</div>
    <h2>Tissue Differences</h2><div class='table-wrap'>{dataframe_to_html(tissue_display)}</div>
    <h2>Annotation Differences</h2><div class='table-wrap'>{dataframe_to_html(annotation_display)}</div>
    <h2>TCR Metadata Differences</h2><div class='table-wrap'>{dataframe_to_html(tcr_display)}</div>
    <h2>Numeric Score And Marker Differences</h2><div class='table-wrap'>{dataframe_to_html(numeric_display)}</div>
    <h2>Targeted TCR And Marker Expression</h2><div class='table-wrap'>{dataframe_to_html(targeted_display)}</div>
    <h2>Top TCR Gene Expression Differences</h2><div class='table-wrap'>{dataframe_to_html(tcr_gene_display)}</div>
    <h2>FN-Enriched Non-TCR Transcriptome Screen</h2><div class='table-wrap'>{dataframe_to_html(fn_gene_display)}</div>
    <h2>TP-Enriched Non-TCR Transcriptome Screen</h2><div class='table-wrap'>{dataframe_to_html(tp_gene_display)}</div>
    <h2>Figures</h2>{''.join(figure_html)}
    <h2>Artifacts</h2><p>Tables: <code>{html.escape(str(TABLE_DIR))}</code><br>Figures: <code>{html.escape(str(FIGURE_DIR))}</code><br>Markdown: <code>{html.escape(str(REPORT_MD))}</code></p>
    </main></body></html>"""
    REPORT_HTML.write_text(html_doc, encoding="utf-8")


def main() -> None:
    args = parse_args()
    ensure_dirs()
    setup_logging()
    input_h5ad = args.input_h5ad.resolve()
    prediction_arrays = args.prediction_arrays.resolve()
    model_path = args.model_pkl.resolve()
    rng = np.random.default_rng(args.seed)
    before_mtime = input_h5ad.stat().st_mtime_ns
    logging.info("Input H5AD: %s", input_h5ad)
    logging.info("Prediction arrays: %s", prediction_arrays)

    with h5py.File(input_h5ad, "r") as handle:
        n_obs = int(handle["obs"]["_index"].shape[0])
        obs = base.build_obs_metadata(handle)
        model = base.load_current_model(model_path)
        with np.load(prediction_arrays) as arrays:
            pred = np.asarray(arrays["annotation_specific_pred"], dtype=bool)
            current_score = np.asarray(arrays["current_score"], dtype=np.float32)
            current_pred = np.asarray(arrays["current_pred"], dtype=bool)
            original_pred = np.asarray(arrays["original_pred"], dtype=bool)
        if pred.shape[0] != n_obs:
            raise RuntimeError(f"Prediction array length {pred.shape[0]} does not match H5AD n_obs {n_obs}")

        validation = validation_positive_set(obs, pred)
        rows = validation.rows
        logging.info("Validation positives: %s TP, %s FN", f"{int((validation.outcome == 'TP').sum()):,}", f"{int((validation.outcome == 'FN').sum()):,}")

        obs_index = read_string_dataset(handle["obs"]["_index"])
        source = clean_group_values(read_obs_column(handle, "source_gse_id"))
        tissue_key = "tissue_corrected" if "tissue_corrected" in handle["obs"] else "tissue"
        tissue = clean_group_values(read_obs_column(handle, tissue_key))
        annotation = read_optional_group(handle, "simple_annotation_plus6", n_obs)
        library_id = read_optional_string(handle, "library_id", n_obs)
        sample_id = read_optional_string(handle, "sample_id", n_obs)
        sublibrary = sublibrary_ids(library_id, sample_id)

        phase4_trd = read_float_obs(handle, "phase4_trd_score")
        phase4_trab = read_float_obs(handle, "phase4_trab_score")
        phase4_delta = read_float_obs(handle, "phase4_trd_minus_trab")

        bool_cols = [
            "TCRseq",
            "Sorted_gdT",
            "has_TRA_TRB_paired",
            "has_TRG_TRD_paired",
            "has_any_ab_tcr",
            "has_any_gd_tcr",
        ]
        cdr3_cols = ["TRA_cdr3", "TRB_cdr3", "TRG_cdr3", "TRD_cdr3"]
        vj_cols = ["TRA_v", "TRB_v", "TRG_v", "TRD_v", "TRA_j", "TRB_j", "TRG_j", "TRD_j"]

        cell_df = pd.DataFrame(
            {
                "obs_row": rows,
                "obs_index": obs_index[rows],
                "outcome": validation.outcome,
                "validation_component": validation.component,
                "source_gse_id": source[rows],
                "tissue": tissue[rows],
                "simple_annotation_plus6": annotation[rows],
                "library_id": library_id[rows],
                "sample_id": sample_id[rows],
                "sublibrary": sublibrary[rows],
                "current_score": current_score[rows],
                "current_pred": current_pred[rows],
                "annotation_specific_pred": pred[rows],
                "original_trd_minus_trab_pred": original_pred[rows],
                "phase4_trd_score": phase4_trd[rows],
                "phase4_trab_score": phase4_trab[rows],
                "phase4_trd_minus_trab": phase4_delta[rows],
                "corrected_has_any_gd_tcr": obs.corrected_has_any_gd_tcr[rows],
            }
        )
        for col in bool_cols:
            cell_df[col] = read_optional_bool(handle, col, n_obs)[rows]
        for col in cdr3_cols:
            cell_df[f"{col}_nonempty"] = read_optional_nonempty(handle, col, n_obs)[rows]
        for col in vj_cols:
            cell_df[col] = read_optional_group(handle, col, n_obs)[rows]

        target_genes = list(dict.fromkeys([*model.gene_names, *TARGETED_MARKER_GENES]))
        spec, _ = build_feature_spec(handle, target_genes, set(model.gene_names))
        x_target = extract_gene_matrix(handle, rows, spec, label=OUT_PREFIX)
        add_signature_scores(cell_df, x_target, spec)

        cell_df.to_csv(TABLE_DIR / "validation_positive_tp_fn_cell_table.csv.gz", index=False)

        source_summary = summarize_group(cell_df, ["source_gse_id"])
        tissue_summary = summarize_group(cell_df, ["tissue"])
        annotation_summary = summarize_group(cell_df, ["simple_annotation_plus6"])
        component_summary = summarize_group(cell_df, ["validation_component"])
        source_tissue_annotation = summarize_group(cell_df, ["source_gse_id", "tissue", "simple_annotation_plus6"])
        source_summary.to_csv(TABLE_DIR / "source_tp_fn_summary.csv", index=False)
        tissue_summary.to_csv(TABLE_DIR / "tissue_tp_fn_summary.csv", index=False)
        annotation_summary.to_csv(TABLE_DIR / "annotation_tp_fn_summary.csv", index=False)
        component_summary.to_csv(TABLE_DIR / "validation_component_tp_fn_summary.csv", index=False)
        source_tissue_annotation.to_csv(TABLE_DIR / "source_tissue_annotation_tp_fn_summary.csv", index=False)

        bool_summary_cols = [
            "TCRseq",
            "Sorted_gdT",
            "has_TRA_TRB_paired",
            "has_TRG_TRD_paired",
            "has_any_ab_tcr",
            "has_any_gd_tcr",
            "corrected_has_any_gd_tcr",
            "TRA_cdr3_nonempty",
            "TRB_cdr3_nonempty",
            "TRG_cdr3_nonempty",
            "TRD_cdr3_nonempty",
            "current_pred",
            "original_trd_minus_trab_pred",
        ]
        bool_df = bool_summary(cell_df, bool_summary_cols)
        bool_df.to_csv(TABLE_DIR / "tcr_metadata_bool_tp_fn_summary.csv", index=False)

        categorical_cols = [
            "validation_component",
            "source_gse_id",
            "tissue",
            "simple_annotation_plus6",
            "TRD_v",
            "TRG_v",
            "TRD_j",
            "TRG_j",
            "TRA_v",
            "TRB_v",
        ]
        cat_df = categorical_summary(cell_df, categorical_cols)
        cat_df.to_csv(TABLE_DIR / "categorical_tp_fn_summary.csv", index=False)

        numeric_cols = [
            "current_score",
            "phase4_trd_score",
            "phase4_trab_score",
            "phase4_trd_minus_trab",
            "expr_t_cell_marker_mean",
            "expr_gdt_tcr_marker_mean",
            "expr_ab_tcr_marker_mean",
            "expr_cytotoxic_nk_marker_mean",
            "expr_penalty_marker_mean",
            "expr_lineage_contamination_marker_mean",
        ] + [col for col in cell_df.columns if col.startswith("expr_") and col not in {
            "expr_t_cell_marker_mean",
            "expr_gdt_tcr_marker_mean",
            "expr_ab_tcr_marker_mean",
            "expr_cytotoxic_nk_marker_mean",
            "expr_penalty_marker_mean",
            "expr_lineage_contamination_marker_mean",
        }]
        numeric_df = numeric_summary(cell_df, numeric_cols)
        numeric_df.to_csv(TABLE_DIR / "numeric_score_marker_tp_fn_summary.csv", index=False)

        targeted_df = summarize_targeted_expression(x_target, spec, validation.outcome)
        targeted_df.to_csv(TABLE_DIR / "targeted_tcr_marker_expression_tp_fn_summary.csv", index=False)

        if args.skip_transcriptome_screen:
            transcriptome_df = pd.DataFrame()
            fn_top_df = pd.DataFrame()
            tp_top_df = pd.DataFrame()
            tcr_top_df = targeted_df.loc[targeted_df["family"].isin(["TRD", "TRG", "TRA", "TRB"])].copy()
        else:
            transcriptome_df = screen_transcriptome(handle, rows, validation.outcome)
            transcriptome_df.to_csv(TABLE_DIR / "transcriptome_tp_fn_gene_screen_all.csv.gz", index=False)
            fn_top_df, tp_top_df, tcr_top_df = top_nontcr_tables(transcriptome_df)
            fn_top_df.to_csv(TABLE_DIR / "transcriptome_fn_enriched_nontcr_top_genes.csv", index=False)
            tp_top_df.to_csv(TABLE_DIR / "transcriptome_tp_enriched_nontcr_top_genes.csv", index=False)
            tcr_top_df.to_csv(TABLE_DIR / "transcriptome_tcr_gene_tp_fn_top.csv", index=False)

        figures = [
            plot_source_fn_rate(source_summary),
            plot_bool_tcr_summary(bool_df),
            plot_score_distributions(cell_df),
            plot_trd_vs_trab(cell_df),
            plot_marker_heatmap(targeted_df),
            plot_top_transcriptome(transcriptome_df),
            plot_umap(handle, cell_df, rows, rng, args.umap_background_sample),
        ]
        copied_assets = copy_assets(figures)

    after_mtime = input_h5ad.stat().st_mtime_ns
    h5ad_mtime_unchanged = before_mtime == after_mtime
    summary_lines = build_summary_text(cell_df, source_summary, bool_df, numeric_df, targeted_df, transcriptome_df)
    write_reports(
        input_h5ad=input_h5ad,
        prediction_arrays=prediction_arrays,
        model_path=model_path,
        cell_df=cell_df,
        source_summary=source_summary,
        tissue_summary=tissue_summary,
        annotation_summary=annotation_summary,
        bool_df=bool_df,
        numeric_df=numeric_df,
        targeted_df=targeted_df,
        tcr_top_df=tcr_top_df,
        fn_top_df=fn_top_df,
        tp_top_df=tp_top_df,
        summary_lines=summary_lines,
        copied_assets=copied_assets,
        h5ad_mtime_unchanged=h5ad_mtime_unchanged,
    )
    summary = {
        "input_h5ad": str(input_h5ad),
        "prediction_arrays": str(prediction_arrays),
        "model_path": str(model_path),
        "n_validation_positive": int(cell_df.shape[0]),
        "tp": int((cell_df["outcome"] == "TP").sum()),
        "fn": int((cell_df["outcome"] == "FN").sum()),
        "h5ad_mtime_unchanged": h5ad_mtime_unchanged,
        "report_html": str(REPORT_HTML),
        "report_pdf": str(REPORT_PDF),
        "tables": str(TABLE_DIR),
        "figures": str(FIGURE_DIR),
        "summary_lines": summary_lines,
    }
    SUMMARY_JSON.write_text(json.dumps(json_ready(summary), indent=2), encoding="utf-8")
    render_pdf(args.no_pdf)
    if not h5ad_mtime_unchanged:
        raise RuntimeError("Input H5AD mtime changed during read-only audit.")
    logging.info("Saved TP/FN audit report: %s", REPORT_HTML)


if __name__ == "__main__":
    main()
