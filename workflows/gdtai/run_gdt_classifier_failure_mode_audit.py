#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Audit failure modes of the selected gdT classifier.

The script is read-only with respect to the input H5AD. It reuses the same
multi-cohort validation definition as
`run_gdt_gse144469_holdout_tcrgene_classifier.py`, reapplies the saved selected
model, and asks whether FP/FN validation cells are enriched for expression-based
low-quality, high-complexity/doublet, or lineage-contamination proxies.
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
import logging
import math
import pickle
import subprocess
from pathlib import Path
from typing import Any

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

from run_gdt_gse144469_holdout_tcrgene_classifier import (
    DEFAULT_INPUT_H5AD,
    EXTRA_TRAB_HOLDOUT_SOURCE,
    GDT2020_HOLDOUT_TISSUE,
    GDT2020_SOURCE,
    HOLDOUT_SOURCE,
    MODEL_PKL,
    RANDOM_SEED,
    SUBOPTIMAL_SOURCE,
    TARGET_SUM,
    TCR_PREFIXES,
    FeatureSpec,
    apply_death_penalties,
    build_obs_metadata,
    tcr_family,
)
from run_gdt_prediction_package_evaluation import read_string_dataset


OUT_ROOT = _TNK_PROJECT_ROOT / "Integrated_dataset"
TABLE_DIR = OUT_ROOT / "tables" / "gdT_prediction" / "gse144469_holdout_tcrgene_failure_modes"
FIGURE_DIR = OUT_ROOT / "figures" / "gdT_prediction" / "gse144469_holdout_tcrgene_failure_modes"
LOG_DIR = OUT_ROOT / "logs" / "gdT_prediction" / "gse144469_holdout_tcrgene_failure_modes"
STATIC_DIR = _TNK_PROJECT_ROOT / "gdT_prediction"
REPORT_MD = LOG_DIR / "gse144469_holdout_tcrgene_failure_mode_audit.md"
REPORT_HTML = STATIC_DIR / "gse144469_holdout_tcrgene_failure_mode_audit.html"
REPORT_PDF = STATIC_DIR / "gse144469_holdout_tcrgene_failure_mode_audit.pdf"
RUN_LOG = LOG_DIR / "run.log"


QC_CONTINUOUS = [
    "total_x_counts",
    "n_genes_by_x",
    "pct_mito_x",
    "pct_ribo_x",
    "cd3_score",
    "cd4_score",
    "foxp3_score",
    "cd8_score",
    "tra_trb_score",
    "trg_trd_score",
    "nk_score",
    "myeloid_score",
    "b_cell_score",
    "epithelial_score",
    "erythroid_score",
    "non_t_lineage_score",
]

QC_FLAGS = [
    "source_high_complexity_q95",
    "source_low_complexity_q05",
    "source_low_quality_proxy",
    "source_high_mito_q95",
    "multi_lineage_doublet_proxy",
    "tcr_dual_gene_proxy",
    "nk_like_proxy",
    "cd4_penalty_proxy",
    "foxp3_penalty_proxy",
    "low_cd3_proxy",
    "any_death_penalty_proxy",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit gdT classifier validation failure modes.")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--model-pkl", type=Path, default=MODEL_PKL)
    parser.add_argument("--holdout-source", default=HOLDOUT_SOURCE)
    parser.add_argument("--extra-trab-holdout-source", default=EXTRA_TRAB_HOLDOUT_SOURCE)
    parser.add_argument("--gdt2020-holdout-tissue", default=GDT2020_HOLDOUT_TISSUE)
    parser.add_argument("--exclude-suboptimal-source", default=SUBOPTIMAL_SOURCE)
    parser.add_argument("--plot-sample-cells", type=int, default=60_000)
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    parser.add_argument("--no-pdf", action="store_true")
    return parser.parse_args()


def ensure_dirs() -> None:
    for path in [TABLE_DIR, FIGURE_DIR, LOG_DIR, STATIC_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def setup_logging() -> None:
    ensure_dirs()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(RUN_LOG, mode="w", encoding="utf-8"), logging.StreamHandler()],
        force=True,
    )


def bh_adjust(p_values: list[float]) -> list[float]:
    p = np.asarray([1.0 if (x is None or not np.isfinite(x)) else float(x) for x in p_values], dtype=float)
    n = p.size
    if n == 0:
        return []
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.empty(n, dtype=float)
    running = 1.0
    for i in range(n - 1, -1, -1):
        running = min(running, ranked[i] * n / float(i + 1))
        adjusted[order[i]] = running
    return adjusted.tolist()


def safe_div(num: float, den: float) -> float:
    return float(num / den) if den else float("nan")


def make_validation_indices_no_write(
    obs: Any,
    *,
    holdout_source: str,
    extra_trab_holdout_source: str,
    gdt2020_holdout_tissue: str,
) -> tuple[np.ndarray, pd.DataFrame]:
    primary = (obs.class_code == 1) | (obs.class_code == 2)
    validation_primary_idx = np.flatnonzero(primary & (obs.source == holdout_source)).astype(np.int64)
    validation_gdt2020_idx = np.flatnonzero(
        primary
        & (obs.source == GDT2020_SOURCE)
        & (obs.class_code == 2)
        & (pd.Series(obs.tissue, copy=False).astype(str).str.lower().to_numpy() == gdt2020_holdout_tissue.lower())
    ).astype(np.int64)
    validation_extra_trab_idx = np.flatnonzero(
        (obs.source == extra_trab_holdout_source)
        & obs.has_TRA_TRB_paired
        & (~obs.corrected_has_any_gd_tcr)
    ).astype(np.int64)
    validation_idx = np.unique(
        np.concatenate([validation_primary_idx, validation_gdt2020_idx, validation_extra_trab_idx])
    ).astype(np.int64)
    parts = [
        (f"validation_primary_{holdout_source}", validation_primary_idx),
        (f"validation_gdT_{GDT2020_SOURCE}_{gdt2020_holdout_tissue.replace(' ', '_')}", validation_gdt2020_idx),
        (f"validation_abT_{extra_trab_holdout_source}_paired_TCRAB_no_gdTCR", validation_extra_trab_idx),
        ("validation_combined", validation_idx),
    ]
    rows = []
    for split, idx in parts:
        labels = obs.class_code[idx]
        if split.startswith("validation_abT_") or split == "validation_combined":
            gdt = int((labels == 2).sum())
            abt = int(idx.size - gdt)
        else:
            gdt = int((labels == 2).sum())
            abt = int((labels == 1).sum())
        rows.append(
            {
                "split": split,
                "n_cells": int(idx.size),
                "gdT_gold": gdt,
                "abT_gold": abt,
                "gdT_prevalence": safe_div(gdt, int(idx.size)),
            }
        )
    return validation_idx, pd.DataFrame(rows)


def feature_spec_from_model(handle: h5py.File, model_payload: dict[str, Any]) -> FeatureSpec:
    var_names = pd.Index(read_string_dataset(handle["var"]["_index"]), dtype="string").astype(str).tolist()
    gene_to_idx = {gene: idx for idx, gene in enumerate(var_names)}
    genes = [str(gene) for gene in model_payload["gene_names"]]
    missing = [gene for gene in genes if gene not in gene_to_idx]
    if missing:
        raise KeyError(f"Model genes missing from H5AD var: {missing[:10]}")
    family = [tcr_family(gene) for gene in genes]
    name_to_col = {gene: i for i, gene in enumerate(genes)}
    return FeatureSpec(
        feature_names=[str(x) for x in model_payload["feature_names"]],
        gene_names=genes,
        gene_indices=np.asarray([gene_to_idx[gene] for gene in genes], dtype=np.int32),
        tcr_gene_names=[gene for gene in genes if gene.upper().startswith(TCR_PREFIXES)],
        penalty_gene_names=[gene for gene in genes if gene in {"FOXP3", "CD4", "CD3D", "CD3E", "CD3G"}],
        family_by_feature=family,
        cd3_cols=[name_to_col[gene] for gene in ["CD3D", "CD3E", "CD3G"] if gene in name_to_col],
        cd4_col=name_to_col.get("CD4"),
        foxp3_col=name_to_col.get("FOXP3"),
    )


def audit_gene_panels(model_genes: list[str]) -> dict[str, list[str]]:
    model_upper = {gene.upper(): gene for gene in model_genes}
    panels = {
        "cd3": ["CD3D", "CD3E", "CD3G"],
        "cd4": ["CD4"],
        "foxp3": ["FOXP3"],
        "cd8": ["CD8A", "CD8B"],
        "nk": ["NKG7", "GNLY", "PRF1", "GZMB", "KLRD1", "FCGR3A", "NCAM1", "TYROBP", "FCER1G"],
        "myeloid": ["LST1", "LYZ", "S100A8", "S100A9", "FCN1", "MS4A7", "CST3"],
        "b_cell": ["MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN"],
        "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT14"],
        "erythroid": ["HBB", "HBA1", "HBA2", "ALAS2"],
        "proliferation": ["MKI67", "TOP2A"],
        "tra_trb": [gene for gene in model_genes if gene.upper().startswith(("TRAV", "TRAJ", "TRAC", "TRBV", "TRBJ", "TRBC"))],
        "trg_trd": [gene for gene in model_genes if gene.upper().startswith(("TRGV", "TRGJ", "TRGC", "TRDV", "TRDD", "TRDJ", "TRDC"))],
    }
    cleaned = {}
    for key, genes in panels.items():
        deduped = []
        for gene in genes:
            actual = model_upper.get(gene.upper(), gene)
            if actual not in deduped:
                deduped.append(actual)
        cleaned[key] = deduped
    return cleaned


def extract_validation_matrix_and_qc(
    handle: h5py.File,
    rows: np.ndarray,
    gene_names: list[str],
) -> tuple[np.ndarray, pd.DataFrame, dict[str, list[str]]]:
    var_names = pd.Index(read_string_dataset(handle["var"]["_index"]), dtype="string").astype(str).tolist()
    var_upper = np.asarray([gene.upper() for gene in var_names], dtype=object)
    gene_to_idx = {gene: idx for idx, gene in enumerate(var_names)}
    panel_defs = audit_gene_panels(gene_names)
    requested = []
    for gene in gene_names:
        if gene in gene_to_idx and gene not in requested:
            requested.append(gene)
    for genes in panel_defs.values():
        for gene in genes:
            if gene in gene_to_idx and gene not in requested:
                requested.append(gene)
    requested_idx = np.asarray([gene_to_idx[gene] for gene in requested], dtype=np.int64)
    requested_to_col = {int(idx): col for col, idx in enumerate(requested_idx)}
    model_cols = [requested.index(gene) for gene in gene_names]

    matrix_all = np.zeros((rows.size, len(requested)), dtype=np.float32)
    total_counts = np.zeros(rows.size, dtype=np.float64)
    n_genes = np.zeros(rows.size, dtype=np.int32)
    mt_counts = np.zeros(rows.size, dtype=np.float64)
    ribo_counts = np.zeros(rows.size, dtype=np.float64)

    is_mt = np.asarray([name.startswith("MT-") for name in var_upper], dtype=bool)
    is_ribo = np.asarray([name.startswith("RPS") or name.startswith("RPL") for name in var_upper], dtype=bool)
    requested_sorted = np.asarray(sorted(requested_to_col), dtype=np.int64)
    x_group = handle["X"]
    indptr = x_group["indptr"]
    indices = x_group["indices"]
    data = x_group["data"]

    for out_row, obs_idx in enumerate(rows):
        start = int(indptr[obs_idx])
        end = int(indptr[obs_idx + 1])
        if end <= start:
            continue
        row_indices = indices[start:end].astype(np.int64, copy=False)
        row_data = data[start:end].astype(np.float32, copy=False)
        row_sum = float(np.sum(row_data, dtype=np.float64))
        total_counts[out_row] = row_sum
        n_genes[out_row] = int(np.count_nonzero(row_data > 0))
        mt_counts[out_row] = float(np.sum(row_data[is_mt[row_indices]], dtype=np.float64))
        ribo_counts[out_row] = float(np.sum(row_data[is_ribo[row_indices]], dtype=np.float64))
        if row_sum <= 0:
            continue
        present = np.isin(row_indices, requested_sorted, assume_unique=False)
        if present.any():
            row_values = np.log1p(row_data[present] * (TARGET_SUM / row_sum)).astype(np.float32, copy=False)
            for gene_idx, value in zip(row_indices[present], row_values):
                matrix_all[out_row, requested_to_col[int(gene_idx)]] = value
        if out_row and out_row % 25_000 == 0:
            logging.info("Extracted QC/gene matrix for %s / %s validation rows", f"{out_row:,}", f"{rows.size:,}")

    qc = pd.DataFrame(
        {
            "total_x_counts": total_counts,
            "n_genes_by_x": n_genes,
            "pct_mito_x": np.divide(mt_counts, total_counts, out=np.zeros_like(mt_counts), where=total_counts > 0) * 100.0,
            "pct_ribo_x": np.divide(ribo_counts, total_counts, out=np.zeros_like(ribo_counts), where=total_counts > 0) * 100.0,
        }
    )
    panel_to_cols = {
        key: [requested.index(gene) for gene in genes if gene in requested]
        for key, genes in panel_defs.items()
    }
    for key, cols in panel_to_cols.items():
        qc[f"{key}_score"] = matrix_all[:, cols].sum(axis=1) if cols else 0.0
    qc["non_t_lineage_score"] = (
        qc["myeloid_score"] + qc["b_cell_score"] + qc["epithelial_score"] + qc["erythroid_score"]
    )
    return matrix_all[:, model_cols].astype(np.float32, copy=False), qc, panel_defs


def add_flags(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["cd4_penalty_proxy"] = out["cd4_score"] > 0.75
    out["foxp3_penalty_proxy"] = out["foxp3_score"] > 0.25
    out["low_cd3_proxy"] = out["cd3_score"] < 0.25
    out["any_death_penalty_proxy"] = (
        out["cd4_penalty_proxy"] | out["foxp3_penalty_proxy"] | out["low_cd3_proxy"]
    )
    out["multi_lineage_doublet_proxy"] = (
        (out["cd3_score"] > 0.25)
        & (
            (out["myeloid_score"] > 0.5)
            | (out["b_cell_score"] > 0.5)
            | (out["epithelial_score"] > 0.5)
            | (out["erythroid_score"] > 0.5)
        )
    )
    out["tcr_dual_gene_proxy"] = (out["tra_trb_score"] > 0.5) & (out["trg_trd_score"] > 0.5)
    out["nk_like_proxy"] = out["nk_score"] > 1.0

    for col in ["total_x_counts", "n_genes_by_x", "pct_mito_x"]:
        q05 = out.groupby("source_gse_id")[col].transform(lambda x: x.quantile(0.05))
        q95 = out.groupby("source_gse_id")[col].transform(lambda x: x.quantile(0.95))
        out[f"{col}_source_q05"] = q05
        out[f"{col}_source_q95"] = q95
    out["source_high_complexity_q95"] = (
        (out["total_x_counts"] > out["total_x_counts_source_q95"])
        | (out["n_genes_by_x"] > out["n_genes_by_x_source_q95"])
    )
    out["source_low_complexity_q05"] = (
        (out["total_x_counts"] < out["total_x_counts_source_q05"])
        | (out["n_genes_by_x"] < out["n_genes_by_x_source_q05"])
    )
    out["source_high_mito_q95"] = out["pct_mito_x"] > out["pct_mito_x_source_q95"]
    out["source_low_quality_proxy"] = out["source_low_complexity_q05"] | out["source_high_mito_q95"]
    return out


def summarize_groups(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name, group in df.groupby("outcome", sort=True):
        row: dict[str, Any] = {
            "outcome": name,
            "n_cells": int(group.shape[0]),
            "fraction_of_validation": safe_div(group.shape[0], df.shape[0]),
            "score_median": float(group["selected_model_score"].median()),
        }
        for col in QC_CONTINUOUS:
            row[f"{col}_median"] = float(group[col].median())
            row[f"{col}_iqr"] = float(group[col].quantile(0.75) - group[col].quantile(0.25))
        for col in QC_FLAGS:
            row[f"{col}_fraction"] = float(group[col].mean()) if group.shape[0] else float("nan")
        rows.append(row)
    return pd.DataFrame(rows)


def compare_continuous(df: pd.DataFrame, contrast: str, a_name: str, b_name: str, cols: list[str]) -> pd.DataFrame:
    a = df[df["outcome"] == a_name]
    b = df[df["outcome"] == b_name]
    rows = []
    for col in cols:
        av = a[col].replace([np.inf, -np.inf], np.nan).dropna().to_numpy(dtype=float)
        bv = b[col].replace([np.inf, -np.inf], np.nan).dropna().to_numpy(dtype=float)
        p = float("nan")
        if av.size and bv.size and (np.unique(np.concatenate([av, bv])).size > 1):
            p = float(mannwhitneyu(av, bv, alternative="two-sided").pvalue)
        rows.append(
            {
                "contrast": contrast,
                "metric": col,
                "group_a": a_name,
                "group_b": b_name,
                "n_a": int(av.size),
                "n_b": int(bv.size),
                "median_a": float(np.median(av)) if av.size else float("nan"),
                "median_b": float(np.median(bv)) if bv.size else float("nan"),
                "median_diff_a_minus_b": (float(np.median(av) - np.median(bv)) if av.size and bv.size else float("nan")),
                "mannwhitney_p": p,
            }
        )
    out = pd.DataFrame(rows)
    out["q_value"] = bh_adjust(out["mannwhitney_p"].tolist())
    return out


def compare_flags(df: pd.DataFrame, contrast: str, a_name: str, b_name: str, flags: list[str]) -> pd.DataFrame:
    a = df[df["outcome"] == a_name]
    b = df[df["outcome"] == b_name]
    rows = []
    for col in flags:
        a_true = int(a[col].sum())
        a_false = int(a.shape[0] - a_true)
        b_true = int(b[col].sum())
        b_false = int(b.shape[0] - b_true)
        p = float("nan")
        odds = float("nan")
        if a.shape[0] and b.shape[0]:
            odds, p = fisher_exact([[a_true, a_false], [b_true, b_false]], alternative="two-sided")
        rows.append(
            {
                "contrast": contrast,
                "flag": col,
                "group_a": a_name,
                "group_b": b_name,
                "n_a": int(a.shape[0]),
                "n_b": int(b.shape[0]),
                "fraction_a": safe_div(a_true, int(a.shape[0])),
                "fraction_b": safe_div(b_true, int(b.shape[0])),
                "fraction_diff_a_minus_b": safe_div(a_true, int(a.shape[0])) - safe_div(b_true, int(b.shape[0])),
                "odds_ratio": float(odds) if np.isfinite(odds) else odds,
                "fisher_p": float(p) if np.isfinite(p) else p,
            }
        )
    out = pd.DataFrame(rows)
    out["q_value"] = bh_adjust(out["fisher_p"].tolist())
    return out


def compare_within_sources(df: pd.DataFrame, a_name: str, b_name: str, flags: list[str], cols: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    flag_rows = []
    cont_rows = []
    for source, source_df in df.groupby("source_gse_id", sort=True):
        if not ((source_df["outcome"] == a_name).any() and (source_df["outcome"] == b_name).any()):
            continue
        flag_rows.append(compare_flags(source_df, f"{a_name}_vs_{b_name}_within_{source}", a_name, b_name, flags))
        cont_rows.append(compare_continuous(source_df, f"{a_name}_vs_{b_name}_within_{source}", a_name, b_name, cols))
    flag_df = pd.concat(flag_rows, ignore_index=True) if flag_rows else pd.DataFrame()
    cont_df = pd.concat(cont_rows, ignore_index=True) if cont_rows else pd.DataFrame()
    return flag_df, cont_df


def plot_qc_boxplots(df: pd.DataFrame, rng: np.random.Generator) -> Path:
    plot_cols = ["total_x_counts", "n_genes_by_x", "pct_mito_x", "pct_ribo_x"]
    plot_df = df.copy()
    plot_df["log10_total_x_counts"] = np.log10(plot_df["total_x_counts"].clip(lower=1))
    plot_df["log10_n_genes_by_x"] = np.log10(plot_df["n_genes_by_x"].clip(lower=1))
    rename = {
        "log10_total_x_counts": "log10 total X counts",
        "log10_n_genes_by_x": "log10 detected genes",
        "pct_mito_x": "% mitochondrial",
        "pct_ribo_x": "% ribosomal",
    }
    plot_vars = ["log10_total_x_counts", "log10_n_genes_by_x", "pct_mito_x", "pct_ribo_x"]
    sample_idx = np.arange(plot_df.shape[0])
    if plot_df.shape[0] > 60_000:
        sample_idx = rng.choice(sample_idx, size=60_000, replace=False)
    sample = plot_df.iloc[sample_idx]
    outcomes = ["TN", "FP", "TP", "FN"]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for ax, col in zip(axes.ravel(), plot_vars):
        data = [sample.loc[sample["outcome"] == outcome, col].dropna().to_numpy() for outcome in outcomes]
        ax.boxplot(data, tick_labels=outcomes, showfliers=False)
        ax.set_title(rename[col])
        ax.set_xlabel("validation outcome")
    path = FIGURE_DIR / "failure_mode_qc_boxplots.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_marker_boxplots(df: pd.DataFrame, rng: np.random.Generator) -> Path:
    plot_vars = ["cd3_score", "trg_trd_score", "tra_trb_score", "nk_score", "non_t_lineage_score", "cd4_score"]
    sample_idx = np.arange(df.shape[0])
    if df.shape[0] > 60_000:
        sample_idx = rng.choice(sample_idx, size=60_000, replace=False)
    sample = df.iloc[sample_idx]
    outcomes = ["TN", "FP", "TP", "FN"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    for ax, col in zip(axes.ravel(), plot_vars):
        data = [sample.loc[sample["outcome"] == outcome, col].dropna().to_numpy() for outcome in outcomes]
        ax.boxplot(data, tick_labels=outcomes, showfliers=False)
        ax.set_title(col)
        ax.set_xlabel("validation outcome")
    path = FIGURE_DIR / "failure_mode_marker_boxplots.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_flag_fractions(summary: pd.DataFrame) -> Path:
    flags = [
        "source_high_complexity_q95_fraction",
        "source_low_quality_proxy_fraction",
        "multi_lineage_doublet_proxy_fraction",
        "tcr_dual_gene_proxy_fraction",
        "nk_like_proxy_fraction",
        "any_death_penalty_proxy_fraction",
    ]
    labels = [flag.replace("_fraction", "").replace("_", " ") for flag in flags]
    outcomes = ["TN", "FP", "TP", "FN"]
    values = []
    for outcome in outcomes:
        row = summary.loc[summary["outcome"] == outcome]
        values.append([float(row[flag].iloc[0]) if not row.empty else 0.0 for flag in flags])
    arr = np.asarray(values)
    x = np.arange(len(labels))
    width = 0.18
    fig, ax = plt.subplots(figsize=(14, 6), constrained_layout=True)
    for i, outcome in enumerate(outcomes):
        ax.bar(x + (i - 1.5) * width, arr[i], width=width, label=outcome)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha="right")
    ax.set_ylabel("fraction of cells")
    ax.set_title("Failure-mode proxy fractions by validation outcome")
    ax.legend(frameon=False, ncol=4)
    path = FIGURE_DIR / "failure_mode_flag_fractions.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def df_to_markdown(df: pd.DataFrame, max_rows: int = 30) -> str:
    if df.empty:
        return "_No rows._"
    view = df.head(max_rows).copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            view[col] = view[col].map(lambda x: "" if pd.isna(x) else f"{x:.4g}")
    view = view.fillna("").astype(str)
    headers = list(view.columns)
    rows = view.values.tolist()

    def escape_cell(value: str) -> str:
        return value.replace("|", "\\|").replace("\n", " ")

    out = ["| " + " | ".join(escape_cell(x) for x in headers) + " |"]
    out.append("| " + " | ".join("---" for _ in headers) + " |")
    for row in rows:
        out.append("| " + " | ".join(escape_cell(x) for x in row) + " |")
    return "\n".join(out)


def df_to_html(df: pd.DataFrame, max_rows: int = 30) -> str:
    if df.empty:
        return "<p><em>No rows.</em></p>"
    view = df.head(max_rows).copy()
    return view.to_html(index=False, escape=True, classes="dataframe")


def write_report(
    *,
    validation_summary: pd.DataFrame,
    group_summary: pd.DataFrame,
    flag_tests: pd.DataFrame,
    continuous_tests: pd.DataFrame,
    source_summary: pd.DataFrame,
    figures: list[Path],
    threshold: float,
    model_notes: str,
    no_pdf: bool,
) -> None:
    key_flags = flag_tests.sort_values(["contrast", "q_value", "flag"]).copy()
    key_cont = continuous_tests.sort_values(["contrast", "q_value", "metric"]).copy()
    group_display_cols = ["outcome", "n_cells", "fraction_of_validation", "score_median"]
    for flag in [
        "source_high_complexity_q95_fraction",
        "source_low_quality_proxy_fraction",
        "multi_lineage_doublet_proxy_fraction",
        "tcr_dual_gene_proxy_fraction",
        "nk_like_proxy_fraction",
        "any_death_penalty_proxy_fraction",
    ]:
        if flag in group_summary:
            group_display_cols.append(flag)
    group_display = group_summary[group_display_cols]
    lines = [
        "# gdT classifier failure-mode audit",
        "",
        f"- Selected model threshold: `{threshold:.8f}`",
        f"- Model notes: {model_notes}",
        "- Explicit doublet/QC obs columns were not present in `integrated_plus6.h5ad`; this audit therefore uses expression-derived proxies.",
        "- Source-specific high/low complexity flags use each validation source's 5th/95th percentiles to reduce source-level QC confounding.",
        "",
        "## Validation cohorts",
        df_to_markdown(validation_summary),
        "",
        "## Outcome-level summary",
        df_to_markdown(group_display, max_rows=20),
        "",
        "## Flag enrichment",
        df_to_markdown(key_flags, max_rows=40),
        "",
        "## Continuous metric shifts",
        df_to_markdown(key_cont, max_rows=40),
        "",
        "## Source and outcome counts",
        df_to_markdown(source_summary, max_rows=80),
        "",
        "## Output files",
        f"- Per-cell audit table: `{TABLE_DIR / 'validation_failure_mode_cell_metrics.csv.gz'}`",
        f"- Group summary: `{TABLE_DIR / 'failure_mode_group_summary.csv'}`",
        f"- Flag enrichment tests: `{TABLE_DIR / 'failure_mode_flag_enrichment.csv'}`",
        f"- Continuous tests: `{TABLE_DIR / 'failure_mode_continuous_tests.csv'}`",
    ]
    for fig in figures:
        lines.append(f"- Figure: `{fig}`")
    REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    style = """
    body{font-family:Arial,Helvetica,sans-serif;margin:24px;color:#1d252c;line-height:1.5}
    h1{font-size:30px} h2{font-size:21px;margin-top:28px}
    table{border-collapse:collapse;width:100%;font-size:12px;margin:10px 0 18px}
    th,td{border:1px solid #d7dde3;padding:5px 7px;text-align:left;vertical-align:top}
    th{background:#eef1f4}.figure{margin:18px 0}.figure img{max-width:100%;border:1px solid #d7dde3}
    code{background:#eef1f4;padding:1px 4px;border-radius:3px}
    @media print{body{margin:12mm} table{font-size:9px} h1{font-size:22px} h2{font-size:16px}}
    """
    figure_html = "\n".join(
        f"<div class='figure'><h3>{html.escape(fig.stem.replace('_', ' '))}</h3>"
        f"<img src='../{html.escape(str(fig))}'></div>"
        for fig in figures
    )
    html_text = f"""<!doctype html><html><head><meta charset='utf-8'>
<title>gdT classifier failure-mode audit</title><style>{style}</style></head><body>
<h1>gdT classifier failure-mode audit</h1>
<p><strong>Threshold:</strong> <code>{threshold:.8f}</code></p>
<p>{html.escape(model_notes)}</p>
<p>Explicit doublet/QC obs columns were not present in <code>integrated_plus6.h5ad</code>; this audit uses expression-derived proxies and source-specific QC quantiles.</p>
<h2>Validation cohorts</h2>{df_to_html(validation_summary)}
<h2>Outcome-level summary</h2>{df_to_html(group_display, max_rows=20)}
<h2>Flag enrichment</h2>{df_to_html(key_flags, max_rows=40)}
<h2>Continuous metric shifts</h2>{df_to_html(key_cont, max_rows=40)}
<h2>Source and outcome counts</h2>{df_to_html(source_summary, max_rows=80)}
<h2>Figures</h2>{figure_html}
</body></html>"""
    REPORT_HTML.write_text(html_text, encoding="utf-8")
    if not no_pdf:
        try:
            subprocess.run(
                [
                    "google-chrome",
                    "--headless",
                    "--disable-gpu",
                    f"--print-to-pdf={REPORT_PDF.resolve()}",
                    REPORT_HTML.resolve().as_uri(),
                ],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
        except Exception as exc:
            logging.warning("PDF rendering failed: %s", exc)


def main() -> None:
    args = parse_args()
    setup_logging()
    rng = np.random.default_rng(args.seed)
    logging.info("Input H5AD: %s", args.input_h5ad)
    logging.info("Model pickle: %s", args.model_pkl)
    with args.model_pkl.open("rb") as fh:
        model_payload = pickle.load(fh)
    threshold = float(model_payload["threshold"])
    model = model_payload["model_object"]
    model_name = str(model_payload.get("model", "selected_model"))
    with h5py.File(args.input_h5ad, "r") as handle:
        obs = build_obs_metadata(handle)
        validation_idx, validation_summary = make_validation_indices_no_write(
            obs,
            holdout_source=args.holdout_source,
            extra_trab_holdout_source=args.extra_trab_holdout_source,
            gdt2020_holdout_tissue=args.gdt2020_holdout_tissue,
        )
        spec = feature_spec_from_model(handle, model_payload)
        x_validation, qc_df, panel_defs = extract_validation_matrix_and_qc(handle, validation_idx, spec.gene_names)

    score = model.predict_proba(x_validation)[:, 1].astype(np.float32)
    if "penalty" in model_name or "death penalties" in str(model_payload.get("notes", "")).lower():
        score = apply_death_penalties(score, x_validation, spec)
    pred = score >= threshold
    truth_gdt = obs.class_code[validation_idx] == 2
    outcome = np.full(validation_idx.size, "TN", dtype=object)
    outcome[(~truth_gdt) & pred] = "FP"
    outcome[truth_gdt & pred] = "TP"
    outcome[truth_gdt & (~pred)] = "FN"

    cell_df = pd.DataFrame(
        {
            "obs_index": validation_idx,
            "source_gse_id": obs.source[validation_idx].astype(str),
            "tissue": obs.tissue[validation_idx].astype(str),
            "truth_label": np.where(truth_gdt, "gdT_gold", "abT_gold"),
            "selected_model_score": score,
            "selected_model_prediction": np.where(pred, "predicted_gdT", "predicted_non_gdT"),
            "outcome": outcome,
        }
    )
    cell_df = pd.concat([cell_df, qc_df], axis=1)
    cell_df = add_flags(cell_df)
    cell_df.to_csv(TABLE_DIR / "validation_failure_mode_cell_metrics.csv.gz", index=False)

    validation_summary.to_csv(TABLE_DIR / "failure_mode_validation_cohorts.csv", index=False)
    group_summary = summarize_groups(cell_df)
    group_summary.to_csv(TABLE_DIR / "failure_mode_group_summary.csv", index=False)
    source_summary = (
        cell_df.groupby(["source_gse_id", "truth_label", "outcome"], sort=True)
        .size()
        .reset_index(name="n_cells")
    )
    source_summary.to_csv(TABLE_DIR / "failure_mode_source_outcome_counts.csv", index=False)

    flag_tests = pd.concat(
        [
            compare_flags(cell_df, "FP_vs_TN", "FP", "TN", QC_FLAGS),
            compare_flags(cell_df, "FN_vs_TP", "FN", "TP", QC_FLAGS),
            compare_flags(cell_df, "all_failures_vs_correct", "failure", "correct", QC_FLAGS)
            if False
            else pd.DataFrame(),
        ],
        ignore_index=True,
    )
    cell_df_for_failure = cell_df.copy()
    cell_df_for_failure["outcome"] = np.where(cell_df_for_failure["outcome"].isin(["FP", "FN"]), "failure", "correct")
    flag_tests = pd.concat(
        [flag_tests, compare_flags(cell_df_for_failure, "all_failures_vs_correct", "failure", "correct", QC_FLAGS)],
        ignore_index=True,
    )
    source_flag_fp, source_cont_fp = compare_within_sources(cell_df, "FP", "TN", QC_FLAGS, QC_CONTINUOUS)
    source_flag_fn, source_cont_fn = compare_within_sources(cell_df, "FN", "TP", QC_FLAGS, QC_CONTINUOUS)
    flag_tests = pd.concat([flag_tests, source_flag_fp, source_flag_fn], ignore_index=True)
    flag_tests.to_csv(TABLE_DIR / "failure_mode_flag_enrichment.csv", index=False)

    continuous_tests = pd.concat(
        [
            compare_continuous(cell_df, "FP_vs_TN", "FP", "TN", QC_CONTINUOUS),
            compare_continuous(cell_df, "FN_vs_TP", "FN", "TP", QC_CONTINUOUS),
            compare_continuous(cell_df_for_failure, "all_failures_vs_correct", "failure", "correct", QC_CONTINUOUS),
            source_cont_fp,
            source_cont_fn,
        ],
        ignore_index=True,
    )
    continuous_tests.to_csv(TABLE_DIR / "failure_mode_continuous_tests.csv", index=False)

    panel_rows = []
    for panel, genes in panel_defs.items():
        panel_rows.append({"panel": panel, "available_genes": ", ".join(genes), "n_available_genes": len(genes)})
    pd.DataFrame(panel_rows).to_csv(TABLE_DIR / "failure_mode_marker_panels.csv", index=False)

    figures = [
        plot_qc_boxplots(cell_df, rng),
        plot_marker_boxplots(cell_df, rng),
        plot_flag_fractions(group_summary),
    ]
    write_report(
        validation_summary=validation_summary,
        group_summary=group_summary,
        flag_tests=flag_tests,
        continuous_tests=continuous_tests,
        source_summary=source_summary,
        figures=figures,
        threshold=threshold,
        model_notes=str(model_payload.get("notes", "")),
        no_pdf=args.no_pdf,
    )
    logging.info("Wrote failure-mode report: %s", REPORT_MD)


if __name__ == "__main__":
    main()
