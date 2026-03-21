#!/usr/bin/env python3
"""Phase 4 gdT module scoring on the integrated milestone.

This helper adapts the shared `gdt_tcr_module_sharing_package_full` logic to the
project-scale integrated AnnData object. The package script is evaluation-first
and expects ground-truth labels; this helper keeps the package-faithful TRA/TRB/TRD
module definitions, but computes only continuous per-cell scores and writes them
back into the canonical integrated H5AD.

Design constraints for this project:
- large H5AD files live on SSD
- tables, logs, and PNG figures live on NFS
- score all cells
- use temporary normalize+log1p for scoring only
- do not create a new milestone H5AD
- do not write hard gdT/abT calls in Phase 4
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import re
import time
from pathlib import Path
from typing import Iterator

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata._io.specs import write_elem
from scipy import sparse
from sklearn.utils import sparsefuncs


PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
FIGURE_DIR = OUTPUT_ROOT / "figures"
TABLE_DIR = OUTPUT_ROOT / "tables"
LOG_DIR = OUTPUT_ROOT / "logs"
HIGH_SPEED_LINK = PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "integrated.h5ad"
SSD_INTEGRATED = Path("/ssd/tnk_phase3/Integrated_dataset/integrated.h5ad")
PACKAGE_ZIP = PROJECT_ROOT / "gdt_tcr_module_sharing_package_full.zip"
PHASE4_LOG = LOG_DIR / "phase4_gdt_module_scoring.log"
PHASE4_QC_MD = LOG_DIR / "phase4_qc_summary.md"

TARGET_SUM = 10_000.0
CHUNK_SIZE = 50_000
N_BINS = 25
CTRL_SIZE = 50
RANDOM_STATE = 1
PLOT_SAMPLE_SIZE = 200_000
SCATTER_PLOT_SAMPLE_SIZE = 80_000
TOP_CELL_N = 5_000
FIGURE_DPI = 300

PHASE4_SCORE_COLUMNS = {
    "tra": "phase4_tra_score",
    "trb": "phase4_trb_score",
    "trab": "phase4_trab_score",
    "trd": "phase4_trd_score",
    "trd_minus_trab": "phase4_trd_minus_trab",
}

PHASE4_SCALED_SCORE_COLUMNS = {
    "trd_scaled": "phase4_trd_score_scaled",
    "trab_scaled": "phase4_trab_score_scaled",
    "trd_minus_trab_scaled": "phase4_trd_minus_trab_scaled",
}

MODULE_PATTERNS = {
    "tra": re.compile(r"^TRAC$|^TRAV|^TRAJ"),
    "trb": re.compile(r"^TRBC|^TRBV|^TRBJ"),
    "trd": re.compile(r"^TRDC$|^TRDV|^TRDJ"),
}

MARKER_GENES = ["TRDC", "TRGC1", "TRGC2", "TRAC", "TRBC1", "TRBC2"]
SCATTER_COLOR_GENES = ["TRDC", "TRDV1", "TRDV2", "NCAM1", "FOXP3", "CD4", "CD8A", "CD8B"]


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description="Run Phase 4 gdT module scoring.")
    parser.add_argument(
        "--integrated-h5ad",
        type=Path,
        default=None,
        help="Optional override for the canonical integrated H5AD.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=CHUNK_SIZE,
        help="Number of cells processed per sparse chunk.",
    )
    parser.add_argument(
        "--plot-sample-size",
        type=int,
        default=PLOT_SAMPLE_SIZE,
        help="Number of cells to sample for UMAP and distribution plots.",
    )
    parser.add_argument(
        "--top-cell-n",
        type=int,
        default=TOP_CELL_N,
        help="Number of top cells to export by phase4_trd_minus_trab.",
    )
    return parser.parse_args()


def configure_logging() -> None:
    """Configure console and file logging."""
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    handlers = [
        logging.FileHandler(PHASE4_LOG, mode="a", encoding="utf-8"),
        logging.StreamHandler(),
    ]
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
        force=True,
    )


def resolve_integrated_h5ad(arg_path: Path | None) -> Path:
    """Resolve the canonical large-file Phase 4 input path."""
    candidates = [arg_path, HIGH_SPEED_LINK, SSD_INTEGRATED, OUTPUT_ROOT / "integrated.h5ad"]
    for candidate in candidates:
        if candidate is not None and candidate.exists():
            return candidate.resolve()
    raise FileNotFoundError("Could not resolve the canonical integrated.h5ad path.")


def read_string_dataset(dataset: h5py.Dataset) -> np.ndarray:
    """Read an HDF5 string dataset into a plain numpy string array."""
    values = dataset[:]
    return np.asarray(
        [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values],
        dtype=object,
    )


def ensure_output_dirs() -> None:
    """Ensure NFS-side Phase 4 output directories exist."""
    for path in [FIGURE_DIR, TABLE_DIR, LOG_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def build_csr_chunk(x_group: h5py.Group, start: int, end: int, n_vars: int) -> sparse.csr_matrix:
    """Build one CSR chunk from the on-disk AnnData `X` group."""
    indptr = x_group["indptr"][start : end + 1].astype(np.int64, copy=False)
    offset = int(indptr[0])
    indptr = indptr - offset
    nnz = int(indptr[-1])
    data = x_group["data"][offset : offset + nnz].astype(np.float32, copy=False)
    indices = x_group["indices"][offset : offset + nnz].astype(np.int32, copy=False)
    return sparse.csr_matrix((data, indices, indptr), shape=(end - start, n_vars))


def iter_csr_chunks(x_group: h5py.Group, n_obs: int, n_vars: int, chunk_size: int) -> Iterator[tuple[int, int, sparse.csr_matrix]]:
    """Yield CSR chunks over the cell axis."""
    for start in range(0, n_obs, chunk_size):
        end = min(start + chunk_size, n_obs)
        yield start, end, build_csr_chunk(x_group, start, end, n_vars)


def normalize_log1p_chunk(chunk: sparse.csr_matrix, target_sum: float) -> sparse.csr_matrix:
    """Normalize a CSR chunk to fixed depth and log1p in place."""
    row_sums = np.asarray(chunk.sum(axis=1)).ravel().astype(np.float32, copy=False)
    scale = np.zeros_like(row_sums, dtype=np.float32)
    nonzero_mask = row_sums > 0
    scale[nonzero_mask] = target_sum / row_sums[nonzero_mask]
    sparsefuncs.inplace_row_scale(chunk, scale)
    np.log1p(chunk.data, out=chunk.data)
    return chunk


def find_module_genes(var_names: pd.Index) -> dict[str, pd.Index]:
    """Find the exact package-faithful TRA/TRB/TRD module genes."""
    modules: dict[str, pd.Index] = {}
    for name, pattern in MODULE_PATTERNS.items():
        genes = [gene for gene in var_names if pattern.match(str(gene))]
        if len(genes) == 0:
            raise ValueError(f"Module `{name}` has no genes in the integrated object.")
        modules[name] = pd.Index(sorted(set(genes)), dtype="string")
    modules["trab"] = pd.Index(sorted(set(modules["tra"]).union(set(modules["trb"]))), dtype="string")
    return modules


def compute_gene_means(integrated_h5ad: Path, n_obs: int, n_vars: int, chunk_size: int) -> np.ndarray:
    """Compute mean log-normalized expression per gene in sparse chunks."""
    logging.info("Computing global mean log-normalized expression for control-gene selection")
    gene_sums = np.zeros(n_vars, dtype=np.float64)
    with h5py.File(integrated_h5ad, "r") as handle:
        x_group = handle["X"]
        for start, end, chunk in iter_csr_chunks(x_group, n_obs, n_vars, chunk_size):
            chunk = normalize_log1p_chunk(chunk, TARGET_SUM)
            gene_sums += np.asarray(chunk.sum(axis=0)).ravel()
            if start == 0 or end == n_obs or (start // chunk_size) % 10 == 0:
                logging.info("Gene-mean pass: processed cells %s-%s / %s", start, end, n_obs)
    return gene_sums / float(n_obs)


def pick_control_genes(gene_list: pd.Index, gene_pool: pd.Index, gene_means: np.ndarray, *, random_state: int) -> pd.Index:
    """Reproduce the Scanpy/Seurat control-gene binning and sampling logic."""
    np.random.seed(random_state)
    obs_avg = pd.Series(gene_means, index=gene_pool)
    obs_avg = obs_avg[np.isfinite(obs_avg)]
    n_items = int(np.round(len(obs_avg) / (N_BINS - 1)))
    if n_items <= 0:
        raise ValueError("Invalid control-gene bin size during Phase 4 setup.")
    obs_cut = obs_avg.rank(method="min") // n_items
    control_genes = pd.Index([], dtype="string")
    for cut in np.unique(obs_cut.loc[gene_list]):
        r_genes = obs_cut[obs_cut == cut].index
        if len(r_genes) == 0:
            continue
        if CTRL_SIZE < len(r_genes):
            r_genes = r_genes.to_series().sample(CTRL_SIZE).index
        r_genes = r_genes.difference(gene_list)
        control_genes = control_genes.union(r_genes)
    if len(control_genes) == 0:
        raise RuntimeError(f"No control genes found for module with {len(gene_list)} genes.")
    return control_genes


def score_chunk(chunk: sparse.csr_matrix, gene_idx: np.ndarray, ctrl_idx: np.ndarray) -> np.ndarray:
    """Score one module on one normalized/log1p CSR chunk."""
    gene_mean = np.asarray(chunk[:, gene_idx].sum(axis=1)).ravel() / float(len(gene_idx))
    ctrl_mean = np.asarray(chunk[:, ctrl_idx].sum(axis=1)).ravel() / float(len(ctrl_idx))
    return (gene_mean - ctrl_mean).astype(np.float32, copy=False)


def minmax_scale(values: np.ndarray) -> tuple[np.ndarray, dict[str, float]]:
    """Scale one score vector to the 0-1 range with global min-max scaling."""
    value_min = float(np.min(values))
    value_max = float(np.max(values))
    if math.isclose(value_min, value_max):
        scaled = np.zeros_like(values, dtype=np.float32)
    else:
        scaled = ((values - value_min) / (value_max - value_min)).astype(np.float32, copy=False)
    return scaled, {"min": value_min, "max": value_max}


def add_scaled_scores(scores: dict[str, np.ndarray]) -> tuple[dict[str, np.ndarray], dict[str, dict[str, float]]]:
    """Add scaled TRD/TRAB scores and scaled TRD-TRAB difference."""
    trd_scaled, trd_stats = minmax_scale(scores["trd"])
    trab_scaled, trab_stats = minmax_scale(scores["trab"])
    scores["trd_scaled"] = trd_scaled
    scores["trab_scaled"] = trab_scaled
    scores["trd_minus_trab_scaled"] = (trd_scaled - trab_scaled).astype(np.float32, copy=False)
    return scores, {"trd_scaled": trd_stats, "trab_scaled": trab_stats}


def update_marker_detection(
    raw_chunk: sparse.csr_matrix,
    marker_idx: np.ndarray,
    cluster_codes: np.ndarray,
    cluster_counts: np.ndarray,
    marker_detection_counts: np.ndarray,
) -> None:
    """Accumulate cluster-level marker detection fractions from raw counts."""
    ones_chunk = raw_chunk[:, marker_idx].copy()
    ones_chunk.data = np.ones_like(ones_chunk.data, dtype=np.float32)
    for cluster_code in np.unique(cluster_codes):
        mask = cluster_codes == cluster_code
        cluster_counts[cluster_code] += int(mask.sum())
        cluster_matrix = ones_chunk[mask]
        marker_detection_counts[cluster_code] += np.asarray(cluster_matrix.sum(axis=0)).ravel().astype(np.int64)


def compute_scores(
    integrated_h5ad: Path,
    *,
    n_obs: int,
    n_vars: int,
    chunk_size: int,
    module_gene_idx: dict[str, np.ndarray],
    module_ctrl_idx: dict[str, np.ndarray],
    leiden_codes: np.ndarray,
    marker_idx: np.ndarray,
) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
    """Compute all Phase 4 scores and cluster-level marker detection counts."""
    logging.info("Computing Phase 4 module scores in sparse chunks")
    scores = {
        "tra": np.zeros(n_obs, dtype=np.float32),
        "trb": np.zeros(n_obs, dtype=np.float32),
        "trab": np.zeros(n_obs, dtype=np.float32),
        "trd": np.zeros(n_obs, dtype=np.float32),
    }
    n_clusters = int(leiden_codes.max()) + 1
    cluster_counts = np.zeros(n_clusters, dtype=np.int64)
    marker_detection_counts = np.zeros((n_clusters, len(marker_idx)), dtype=np.int64)

    with h5py.File(integrated_h5ad, "r") as handle:
        x_group = handle["X"]
        for start, end, raw_chunk in iter_csr_chunks(x_group, n_obs, n_vars, chunk_size):
            chunk_cluster_codes = leiden_codes[start:end]
            update_marker_detection(raw_chunk, marker_idx, chunk_cluster_codes, cluster_counts, marker_detection_counts)
            score_chunk_matrix = normalize_log1p_chunk(raw_chunk, TARGET_SUM)
            for module_name in ("tra", "trb", "trab", "trd"):
                scores[module_name][start:end] = score_chunk(
                    score_chunk_matrix,
                    module_gene_idx[module_name],
                    module_ctrl_idx[module_name],
                )
            if start == 0 or end == n_obs or (start // chunk_size) % 10 == 0:
                logging.info("Scoring pass: processed cells %s-%s / %s", start, end, n_obs)

    scores["trd_minus_trab"] = (scores["trd"] - scores["trab"]).astype(np.float32, copy=False)
    return scores, cluster_counts, marker_detection_counts


def read_obs_strings(handle: h5py.File, column_name: str) -> np.ndarray:
    """Read one string obs column from the integrated H5AD."""
    return read_string_dataset(handle["obs"][column_name])


def read_umap(handle: h5py.File) -> np.ndarray:
    """Read the Phase 3 UMAP embedding."""
    return handle["obsm"]["X_umap"][:].astype(np.float32, copy=False)


def append_obs_columns_in_place(h5ad_path: Path, score_columns: dict[str, np.ndarray], uns_payload: dict) -> None:
    """Append Phase 4 columns to obs and Phase 4 metadata to uns in place."""
    logging.info("Writing Phase 4 score columns back into %s", h5ad_path)
    with h5py.File(h5ad_path, "r+") as handle:
        obs_group = handle["obs"]
        column_order = list(obs_group.attrs["column-order"])
        final_names = list(score_columns.keys())

        for final_name, values in score_columns.items():
            tmp_name = f"__{final_name}_tmp"
            for existing in (tmp_name, final_name):
                if existing in obs_group:
                    del obs_group[existing]
            write_elem(obs_group, tmp_name, np.asarray(values, dtype=np.float32))
            obs_group.move(tmp_name, final_name)
            if final_name not in column_order:
                column_order.append(final_name)

        obs_group.attrs["column-order"] = np.asarray(column_order, dtype=object)

        uns_group = handle["uns"]
        if "phase4_gdt_module" in uns_group:
            del uns_group["phase4_gdt_module"]
        write_elem(uns_group, "phase4_gdt_module", uns_payload)


def select_plot_sample(leiden_labels: np.ndarray, sample_size: int, random_state: int) -> np.ndarray:
    """Sample cells for Phase 4 plots while preserving small clusters."""
    rng = np.random.default_rng(random_state)
    leiden_series = pd.Series(leiden_labels, dtype="string")
    cluster_counts = leiden_series.value_counts(sort=False)
    total_cells = int(cluster_counts.sum())
    target = min(sample_size, total_cells)
    sampled_idx: list[np.ndarray] = []
    for cluster, count in cluster_counts.items():
        cluster_idx = np.flatnonzero(leiden_series.to_numpy() == cluster)
        if count <= 0:
            continue
        cluster_target = int(round(target * (count / total_cells)))
        cluster_target = max(500, cluster_target)
        cluster_target = min(cluster_target, int(count))
        sampled_idx.append(rng.choice(cluster_idx, size=cluster_target, replace=False))
    sample = np.concatenate(sampled_idx)
    if len(sample) > target:
        sample = rng.choice(sample, size=target, replace=False)
    return np.sort(sample)


def load_selected_strings(dataset: h5py.Dataset, indices: np.ndarray) -> np.ndarray:
    """Load selected rows from a string dataset in original order."""
    sorted_idx = np.sort(indices)
    selected = dataset[sorted_idx]
    sorted_vals = np.asarray(
        [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in selected],
        dtype=object,
    )
    position_map = {idx: pos for pos, idx in enumerate(sorted_idx.tolist())}
    return np.asarray([sorted_vals[position_map[idx]] for idx in indices.tolist()], dtype=object)


def downsample_indices(indices: np.ndarray, target_size: int, random_state: int) -> np.ndarray:
    """Downsample a sorted index array without replacement."""
    if indices.size <= target_size:
        return indices
    rng = np.random.default_rng(random_state)
    sampled = rng.choice(indices, size=target_size, replace=False)
    return np.sort(sampled.astype(np.int64, copy=False))


def extract_log1p_gene_expression_for_sample(
    integrated_h5ad: Path,
    sample_idx: np.ndarray,
    gene_names: list[str],
    chunk_size: int,
) -> pd.DataFrame:
    """Extract temporary normalize+log1p expression values for selected genes on sampled cells."""
    sample_idx = np.asarray(sample_idx, dtype=np.int64)
    with h5py.File(integrated_h5ad, "r") as handle:
        var_names = pd.Index(read_string_dataset(handle["var"]["_index"]), dtype="string")
        gene_idx = var_names.get_indexer(pd.Index(gene_names, dtype="string"))
        if np.any(gene_idx < 0):
            missing = [gene for gene, idx in zip(gene_names, gene_idx.tolist()) if idx < 0]
            raise ValueError(f"Missing genes for scatter panel: {missing}")

        x_group = handle["X"]
        n_obs = int(x_group["indptr"].shape[0] - 1)
        n_vars = int(handle["var"]["_index"].shape[0])
        expr_values = np.zeros((sample_idx.size, len(gene_names)), dtype=np.float32)

        for start, end, raw_chunk in iter_csr_chunks(x_group, n_obs, n_vars, chunk_size):
            lo = int(np.searchsorted(sample_idx, start, side="left"))
            hi = int(np.searchsorted(sample_idx, end, side="left"))
            if lo == hi:
                continue
            local_idx = sample_idx[lo:hi] - start
            selected_chunk = raw_chunk[local_idx].copy()
            selected_chunk = normalize_log1p_chunk(selected_chunk, TARGET_SUM)
            expr_values[lo:hi, :] = selected_chunk[:, gene_idx].toarray().astype(np.float32, copy=False)

    expr_df = pd.DataFrame(expr_values, columns=gene_names)
    return expr_df


def write_tables(
    integrated_h5ad: Path,
    scores: dict[str, np.ndarray],
    leiden_labels: np.ndarray,
    leiden_categories: np.ndarray,
    source_gse_id: np.ndarray,
    cluster_counts: np.ndarray,
    marker_detection_counts: np.ndarray,
    marker_genes: list[str],
    module_gene_sets: dict[str, pd.Index],
    top_cell_n: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Write Phase 4 summary tables to NFS."""
    logging.info("Writing Phase 4 tables")

    overall_rows = []
    score_column_map = {**PHASE4_SCORE_COLUMNS, **PHASE4_SCALED_SCORE_COLUMNS}
    for module_name, column_name in score_column_map.items():
        values = scores[module_name]
        overall_rows.append(
            {
                "score_name": column_name,
                "n_cells": values.size,
                "mean": float(np.mean(values)),
                "std": float(np.std(values)),
                "min": float(np.min(values)),
                "p01": float(np.quantile(values, 0.01)),
                "p05": float(np.quantile(values, 0.05)),
                "median": float(np.quantile(values, 0.50)),
                "p95": float(np.quantile(values, 0.95)),
                "p99": float(np.quantile(values, 0.99)),
                "max": float(np.max(values)),
            }
        )
    overall_df = pd.DataFrame(overall_rows)
    overall_df.to_csv(TABLE_DIR / "phase4_score_summary.csv", index=False)

    membership_rows = []
    for module_name, genes in module_gene_sets.items():
        for gene in genes:
            membership_rows.append({"module_name": module_name, "gene_symbol": str(gene)})
    membership_df = pd.DataFrame(membership_rows)
    membership_df.to_csv(TABLE_DIR / "phase4_module_gene_membership.csv", index=False)

    summary_df = pd.DataFrame(
        {
            "leiden": pd.Categorical(leiden_labels),
            "source_gse_id": pd.Categorical(source_gse_id),
            **{column_name: scores[module_name] for module_name, column_name in score_column_map.items()},
        }
    )

    leiden_summary = summary_df.groupby("leiden", observed=True).agg(
        n_cells=("leiden", "size"),
        phase4_tra_score_mean=("phase4_tra_score", "mean"),
        phase4_trb_score_mean=("phase4_trb_score", "mean"),
        phase4_trab_score_mean=("phase4_trab_score", "mean"),
        phase4_trd_score_mean=("phase4_trd_score", "mean"),
        phase4_trd_minus_trab_mean=("phase4_trd_minus_trab", "mean"),
        phase4_trd_minus_trab_median=("phase4_trd_minus_trab", "median"),
        phase4_trd_score_scaled_mean=("phase4_trd_score_scaled", "mean"),
        phase4_trab_score_scaled_mean=("phase4_trab_score_scaled", "mean"),
        phase4_trd_minus_trab_scaled_mean=("phase4_trd_minus_trab_scaled", "mean"),
        phase4_trd_minus_trab_scaled_median=("phase4_trd_minus_trab_scaled", "median"),
    ).reset_index()
    marker_fraction_df = pd.DataFrame({"leiden": leiden_categories.astype(object)})
    for marker_idx, marker_gene in enumerate(marker_genes):
        marker_fraction_df[f"{marker_gene}_detected_fraction"] = (
            marker_detection_counts[:, marker_idx] / np.maximum(cluster_counts, 1)
        )
    leiden_summary = leiden_summary.merge(marker_fraction_df, on="leiden", how="left", validate="one_to_one")
    leiden_summary = leiden_summary.sort_values("phase4_trd_minus_trab_median", ascending=False)
    leiden_summary.to_csv(TABLE_DIR / "phase4_leiden_score_summary.csv", index=False)

    gse_summary = summary_df.groupby("source_gse_id", observed=True).agg(
        n_cells=("source_gse_id", "size"),
        phase4_tra_score_mean=("phase4_tra_score", "mean"),
        phase4_trb_score_mean=("phase4_trb_score", "mean"),
        phase4_trab_score_mean=("phase4_trab_score", "mean"),
        phase4_trd_score_mean=("phase4_trd_score", "mean"),
        phase4_trd_minus_trab_mean=("phase4_trd_minus_trab", "mean"),
        phase4_trd_minus_trab_median=("phase4_trd_minus_trab", "median"),
        phase4_trd_score_scaled_mean=("phase4_trd_score_scaled", "mean"),
        phase4_trab_score_scaled_mean=("phase4_trab_score_scaled", "mean"),
        phase4_trd_minus_trab_scaled_mean=("phase4_trd_minus_trab_scaled", "mean"),
        phase4_trd_minus_trab_scaled_median=("phase4_trd_minus_trab_scaled", "median"),
    ).reset_index()
    gse_summary = gse_summary.sort_values("n_cells", ascending=False)
    gse_summary.to_csv(TABLE_DIR / "phase4_gse_score_summary.csv", index=False)

    top_idx = np.argpartition(scores["trd_minus_trab"], -top_cell_n)[-top_cell_n:]
    top_idx = top_idx[np.argsort(scores["trd_minus_trab"][top_idx])[::-1]]
    with h5py.File(integrated_h5ad, "r") as handle:
        top_df = pd.DataFrame(
            {
                "cell_id": load_selected_strings(handle["obs"]["_index"], top_idx),
                "original_cell_id": load_selected_strings(handle["obs"]["original_cell_id"], top_idx),
                "source_gse_id": load_selected_strings(handle["obs"]["source_gse_id"], top_idx),
                "sampleid": load_selected_strings(handle["obs"]["sampleid"], top_idx),
                "project_name": load_selected_strings(handle["obs"]["project name"], top_idx),
                "barcodes": load_selected_strings(handle["obs"]["barcodes"], top_idx),
                "leiden": load_selected_strings(handle["obs"]["leiden"], top_idx),
                "phase4_tra_score": scores["tra"][top_idx],
                "phase4_trb_score": scores["trb"][top_idx],
                "phase4_trab_score": scores["trab"][top_idx],
                "phase4_trd_score": scores["trd"][top_idx],
                "phase4_trd_minus_trab": scores["trd_minus_trab"][top_idx],
                "phase4_trab_score_scaled": scores["trab_scaled"][top_idx],
                "phase4_trd_score_scaled": scores["trd_scaled"][top_idx],
                "phase4_trd_minus_trab_scaled": scores["trd_minus_trab_scaled"][top_idx],
            }
        )
    top_df.to_csv(TABLE_DIR / "phase4_top_cells_by_trd_minus_trab.csv", index=False)
    return leiden_summary, gse_summary


def write_figures(
    *,
    sample_df: pd.DataFrame,
    leiden_summary: pd.DataFrame,
    gse_summary: pd.DataFrame,
) -> None:
    """Generate Phase 4 PNG figures on NFS."""
    logging.info("Writing Phase 4 PNG figures")
    sns.set_theme(style="whitegrid", context="talk")

    dist_fig, dist_axes = plt.subplots(2, 3, figsize=(18, 10), constrained_layout=True)
    dist_specs = [
        ("phase4_trd_score", "TRD score"),
        ("phase4_trab_score", "TRAB score"),
        ("phase4_trd_minus_trab", "TRD - TRAB"),
        ("phase4_trd_score_scaled", "TRD score scaled 0-1"),
        ("phase4_trab_score_scaled", "TRAB score scaled 0-1"),
        ("phase4_trd_minus_trab_scaled", "Scaled TRD - TRAB"),
    ]
    for ax, (column, title) in zip(dist_axes.flatten(), dist_specs):
        sns.histplot(sample_df[column], bins=100, ax=ax, color="#2F6690")
        ax.set_title(title)
        ax.set_xlabel(column)
    dist_fig.savefig(FIGURE_DIR / "phase4_score_distributions.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(dist_fig)

    umap_fig, umap_axes = plt.subplots(2, 3, figsize=(18, 10), constrained_layout=True)
    for ax, (column, title) in zip(umap_axes.flatten(), dist_specs):
        scatter = ax.scatter(
            sample_df["umap1"],
            sample_df["umap2"],
            c=sample_df[column],
            cmap="viridis",
            s=3,
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(title)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        umap_fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
    umap_fig.savefig(FIGURE_DIR / "phase4_umap_score_overlays.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(umap_fig)

    heatmap_cols = [
        "phase4_tra_score_mean",
        "phase4_trb_score_mean",
        "phase4_trab_score_mean",
        "phase4_trd_score_mean",
        "phase4_trd_minus_trab_median",
        "phase4_trd_score_scaled_mean",
        "phase4_trab_score_scaled_mean",
        "phase4_trd_minus_trab_scaled_median",
    ]
    leiden_heatmap = leiden_summary.set_index("leiden")[heatmap_cols]
    heatmap_fig, heatmap_ax = plt.subplots(figsize=(10, 10), constrained_layout=True)
    sns.heatmap(leiden_heatmap, cmap="vlag", center=0, ax=heatmap_ax)
    heatmap_ax.set_title("Phase 4 Score Summary by Leiden Cluster")
    heatmap_fig.savefig(FIGURE_DIR / "phase4_leiden_score_summary.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(heatmap_fig)

    gse_plot_df = gse_summary.sort_values("n_cells", ascending=False).head(20).copy()
    gse_fig, gse_ax = plt.subplots(figsize=(12, 7), constrained_layout=True)
    sns.barplot(
        data=gse_plot_df,
        y="source_gse_id",
        x="phase4_trd_minus_trab_median",
        color="#4C956C",
        ax=gse_ax,
    )
    gse_ax.set_title("Top GSEs by Cell Count: Median TRD - TRAB")
    gse_ax.set_xlabel("Median phase4_trd_minus_trab")
    gse_ax.set_ylabel("source_gse_id")
    gse_fig.savefig(FIGURE_DIR / "phase4_gse_score_summary.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(gse_fig)

    marker_cols = [
        "phase4_trd_score_mean",
        "phase4_trab_score_mean",
        "phase4_trd_minus_trab_median",
        "TRDC_detected_fraction",
        "TRGC1_detected_fraction",
        "TRGC2_detected_fraction",
        "TRAC_detected_fraction",
        "TRBC1_detected_fraction",
        "TRBC2_detected_fraction",
    ]
    marker_heatmap = leiden_summary.set_index("leiden")[marker_cols]
    marker_fig, marker_ax = plt.subplots(figsize=(12, 10), constrained_layout=True)
    sns.heatmap(marker_heatmap, cmap="rocket", ax=marker_ax)
    marker_ax.set_title("Phase 4 Scores and TCR Marker Detection by Leiden Cluster")
    marker_fig.savefig(FIGURE_DIR / "phase4_marker_score_comparison.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(marker_fig)


def write_trab_trd_scatter_panel(sample_df: pd.DataFrame) -> None:
    """Write one figure containing raw-score and scaled-score TRAB-vs-TRD scatter panels."""
    logging.info("Writing Phase 4 TRAB-vs-TRD scatter panel")
    raw_color_specs = [("phase4_trab_minus_trd", "TRAB - TRD")] + [(gene, gene) for gene in SCATTER_COLOR_GENES]
    scaled_color_specs = [("phase4_trab_minus_trd_scaled", "Scaled TRAB - TRD")] + [
        (gene, gene) for gene in SCATTER_COLOR_GENES
    ]

    fig, axes = plt.subplots(2, len(raw_color_specs), figsize=(4 * len(raw_color_specs), 8), constrained_layout=True)

    def draw_row(row_axes: np.ndarray, x_col: str, y_col: str, color_specs: list[tuple[str, str]], row_title: str) -> None:
        for ax, (color_col, title) in zip(row_axes, color_specs):
            cmap = "coolwarm" if "TRAB - TRD" in title else "viridis"
            scatter = ax.scatter(
                sample_df[x_col],
                sample_df[y_col],
                c=sample_df[color_col],
                cmap=cmap,
                s=3,
                linewidths=0,
                rasterized=True,
            )
            ax.set_title(title)
            ax.set_xlabel(x_col)
            ax.set_ylabel(y_col)
            fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
        row_axes[0].annotate(
            row_title,
            xy=(-0.38, 0.5),
            xycoords="axes fraction",
            rotation=90,
            va="center",
            ha="center",
            fontsize=16,
            fontweight="bold",
        )

    draw_row(axes[0], "phase4_trab_score", "phase4_trd_score", raw_color_specs, "Raw scores")
    draw_row(
        axes[1],
        "phase4_trab_score_scaled",
        "phase4_trd_score_scaled",
        scaled_color_specs,
        "Scaled scores",
    )
    fig.savefig(FIGURE_DIR / "phase4_trab_vs_trd_scatter_panel.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def write_qc_summary(
    integrated_h5ad: Path,
    *,
    module_gene_sets: dict[str, pd.Index],
    module_control_sets: dict[str, pd.Index],
    overall_score_summary: pd.DataFrame,
    leiden_summary: pd.DataFrame,
    gse_summary: pd.DataFrame,
) -> None:
    """Write the Phase 4 QC summary markdown."""
    top_leiden = leiden_summary.head(5)[
        ["leiden", "phase4_trd_minus_trab_median", "phase4_trd_minus_trab_scaled_median", "n_cells"]
    ]
    top_gse = gse_summary.sort_values("phase4_trd_minus_trab_median", ascending=False).head(5)[
        ["source_gse_id", "phase4_trd_minus_trab_median", "phase4_trd_minus_trab_scaled_median", "n_cells"]
    ]
    lines = [
        "# Phase 4 QC Summary",
        "",
        "## Scope",
        f"- Input and updated output milestone: `{integrated_h5ad}`",
        f"- Package source: `{PACKAGE_ZIP}`",
        "- Scoring mode: exact package TRA/TRB/TRD modules on a temporary normalized/log1p copy of count-space `X`",
        "- Canonical result type: continuous scores only; no hard gdT/abT call was written in Phase 4",
        "- Derived scaled outputs: min-max scaled `phase4_trd_score_scaled` and `phase4_trab_score_scaled` in the 0-1 range, plus `phase4_trd_minus_trab_scaled` as their difference",
        "- scANVI labels remain reference-only and were not used to define Phase 4 outputs",
        "",
        "## Module sizes",
    ]
    for module_name in ("tra", "trb", "trab", "trd"):
        lines.append(
            f"- `{module_name}` genes: {len(module_gene_sets[module_name])}; control genes used: {len(module_control_sets[module_name])}"
        )
    lines.extend(
        [
            "",
            "## Score summary",
        ]
    )
    for _, row in overall_score_summary.iterrows():
        lines.append(
            f"- `{row['score_name']}`: mean={row['mean']:.4f}, median={row['median']:.4f}, p95={row['p95']:.4f}, max={row['max']:.4f}"
        )
    lines.extend(["", "## Top Leiden clusters by median TRD - TRAB"])
    for _, row in top_leiden.iterrows():
        lines.append(
            f"- Leiden `{row['leiden']}`: raw_median={row['phase4_trd_minus_trab_median']:.4f}, scaled_median={row['phase4_trd_minus_trab_scaled_median']:.4f}, n_cells={int(row['n_cells'])}"
        )
    lines.extend(["", "## Top GSEs by median TRD - TRAB"])
    for _, row in top_gse.iterrows():
        lines.append(
            f"- `{row['source_gse_id']}`: raw_median={row['phase4_trd_minus_trab_median']:.4f}, scaled_median={row['phase4_trd_minus_trab_scaled_median']:.4f}, n_cells={int(row['n_cells'])}"
        )
    lines.extend(
        [
            "",
            "## Outputs",
            "- Tables: `phase4_score_summary.csv`, `phase4_module_gene_membership.csv`, `phase4_leiden_score_summary.csv`, `phase4_gse_score_summary.csv`, `phase4_top_cells_by_trd_minus_trab.csv`",
            "- Figures: `phase4_score_distributions.png`, `phase4_umap_score_overlays.png`, `phase4_leiden_score_summary.png`, `phase4_gse_score_summary.png`, `phase4_marker_score_comparison.png`, `phase4_trab_vs_trd_scatter_panel.png`",
            "",
            "## QC conclusion",
            "- Phase 4 wrote continuous TRA/TRB/TRAB/TRD module scores and `TRD - TRAB` back into the canonical integrated milestone.",
            "- Phase 4 also wrote min-max scaled `TRD` and `TRAB` scores in the 0-1 range plus the scaled `TRD - TRAB` difference.",
            "- The shared package was used as the method source of truth for module definitions, but its ground-truth evaluation branch was intentionally not used here.",
            "- Phase 4 results should be interpreted jointly with Leiden structure, marker genes, and the scVI embedding.",
        ]
    )
    PHASE4_QC_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """Run the full Phase 4 scoring workflow."""
    args = parse_args()
    ensure_output_dirs()
    configure_logging()

    start_time = time.time()
    integrated_h5ad = resolve_integrated_h5ad(args.integrated_h5ad)
    logging.info("Using integrated H5AD at %s", integrated_h5ad)

    with h5py.File(integrated_h5ad, "r") as handle:
        obs_group = handle["obs"]
        x_group = handle["X"]
        obsm_group = handle["obsm"]
        if "X_scVI" not in obsm_group or "X_umap" not in obsm_group:
            raise ValueError("Phase 4 requires `X_scVI` and `X_umap` from Phase 3.")
        required_obs = {"leiden", "source_gse_id"}
        missing_obs = [col for col in required_obs if col not in obs_group]
        if missing_obs:
            raise ValueError(f"Phase 4 requires missing obs columns: {missing_obs}")
        n_obs = int(x_group["indptr"].shape[0] - 1)
        n_vars = int(handle["var"]["_index"].shape[0])
        var_names = pd.Index(read_string_dataset(handle["var"]["_index"]), dtype="string")
        leiden_labels = read_obs_strings(handle, "leiden")
        source_gse_id = read_obs_strings(handle, "source_gse_id")

    modules = find_module_genes(var_names)
    for module_name, genes in modules.items():
        logging.info("Module %s contains %s genes", module_name, len(genes))

    gene_means = compute_gene_means(integrated_h5ad, n_obs, n_vars, args.chunk_size)
    gene_pool = pd.Index(var_names, dtype="string")
    module_controls = {
        module_name: pick_control_genes(gene_list, gene_pool, gene_means, random_state=RANDOM_STATE)
        for module_name, gene_list in modules.items()
    }

    module_gene_idx = {name: var_names.get_indexer(genes).astype(np.int32) for name, genes in modules.items()}
    module_ctrl_idx = {
        name: var_names.get_indexer(ctrl_genes).astype(np.int32) for name, ctrl_genes in module_controls.items()
    }
    marker_idx = var_names.get_indexer(pd.Index(MARKER_GENES, dtype="string"))
    if np.any(marker_idx < 0):
        missing_markers = [gene for gene, idx in zip(MARKER_GENES, marker_idx.tolist()) if idx < 0]
        raise ValueError(f"Missing required marker genes for Phase 4 QC: {missing_markers}")

    leiden_codes, leiden_categories = pd.factorize(leiden_labels, sort=True)
    scores, cluster_counts, marker_detection_counts = compute_scores(
        integrated_h5ad,
        n_obs=n_obs,
        n_vars=n_vars,
        chunk_size=args.chunk_size,
        module_gene_idx=module_gene_idx,
        module_ctrl_idx=module_ctrl_idx,
        leiden_codes=leiden_codes.astype(np.int32, copy=False),
        marker_idx=marker_idx.astype(np.int32, copy=False),
    )
    scores, scaling_stats = add_scaled_scores(scores)

    uns_payload = {
        "package_source": str(PACKAGE_ZIP),
        "integrated_h5ad": str(integrated_h5ad),
        "scoring_mode": "temporary_normalize_total_log1p_on_count_space_X",
        "continuous_only": True,
        "random_state": RANDOM_STATE,
        "target_sum": TARGET_SUM,
        "ctrl_size": CTRL_SIZE,
        "n_bins": N_BINS,
        "module_genes": {name: [str(gene) for gene in genes] for name, genes in modules.items()},
        "control_genes": {name: [str(gene) for gene in genes] for name, genes in module_controls.items()},
        "score_columns": {**PHASE4_SCORE_COLUMNS, **PHASE4_SCALED_SCORE_COLUMNS},
        "scaled_score_columns": PHASE4_SCALED_SCORE_COLUMNS,
        "scaled_score_ranges": scaling_stats,
        "scANVI_usage": "reference_only_not_used_for_phase4_calls",
    }
    append_obs_columns_in_place(
        integrated_h5ad,
        {
            **{column_name: scores[module_name] for module_name, column_name in PHASE4_SCORE_COLUMNS.items()},
            **{column_name: scores[module_name] for module_name, column_name in PHASE4_SCALED_SCORE_COLUMNS.items()},
        },
        uns_payload,
    )

    leiden_summary, gse_summary = write_tables(
        integrated_h5ad,
        scores,
        leiden_labels,
        leiden_categories.astype(object),
        source_gse_id,
        cluster_counts,
        marker_detection_counts,
        MARKER_GENES,
        modules,
        args.top_cell_n,
    )

    with h5py.File(integrated_h5ad, "r") as handle:
        umap = read_umap(handle)
    sample_idx = select_plot_sample(leiden_labels, args.plot_sample_size, RANDOM_STATE)
    sample_df = pd.DataFrame(
        {
            "umap1": umap[sample_idx, 0],
            "umap2": umap[sample_idx, 1],
            "leiden": leiden_labels[sample_idx],
            "source_gse_id": source_gse_id[sample_idx],
            "phase4_trd_score": scores["trd"][sample_idx],
            "phase4_trab_score": scores["trab"][sample_idx],
            "phase4_trd_minus_trab": scores["trd_minus_trab"][sample_idx],
            "phase4_trd_score_scaled": scores["trd_scaled"][sample_idx],
            "phase4_trab_score_scaled": scores["trab_scaled"][sample_idx],
            "phase4_trd_minus_trab_scaled": scores["trd_minus_trab_scaled"][sample_idx],
        }
    )
    sample_df["phase4_trab_minus_trd"] = sample_df["phase4_trab_score"] - sample_df["phase4_trd_score"]
    sample_df["phase4_trab_minus_trd_scaled"] = (
        sample_df["phase4_trab_score_scaled"] - sample_df["phase4_trd_score_scaled"]
    )
    scatter_idx = downsample_indices(sample_idx, SCATTER_PLOT_SAMPLE_SIZE, RANDOM_STATE)
    scatter_selector = np.searchsorted(sample_idx, scatter_idx)
    scatter_df = sample_df.iloc[scatter_selector].reset_index(drop=True)
    scatter_df = pd.concat(
        [
            scatter_df,
            extract_log1p_gene_expression_for_sample(
                integrated_h5ad,
                scatter_idx,
                SCATTER_COLOR_GENES,
                args.chunk_size,
            ),
        ],
        axis=1,
    )
    write_figures(sample_df=sample_df, leiden_summary=leiden_summary, gse_summary=gse_summary)
    write_trab_trd_scatter_panel(scatter_df)

    overall_score_summary = pd.read_csv(TABLE_DIR / "phase4_score_summary.csv")
    write_qc_summary(
        integrated_h5ad,
        module_gene_sets=modules,
        module_control_sets=module_controls,
        overall_score_summary=overall_score_summary,
        leiden_summary=leiden_summary,
        gse_summary=gse_summary,
    )

    elapsed = time.time() - start_time
    logging.info(
        "Phase 4 complete: cells=%s genes=%s elapsed_minutes=%.2f output=%s",
        n_obs,
        n_vars,
        elapsed / 60.0,
        integrated_h5ad,
    )


if __name__ == "__main__":
    main()
