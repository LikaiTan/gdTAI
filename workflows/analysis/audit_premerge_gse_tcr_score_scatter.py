#!/usr/bin/env python3
"""Audit pre-merge per-GSE H5ADs for suspicious paired TRA/TRB coverage.

This helper works on the individual GSE H5ADs before merge. For each selected
dataset, it:

1. finds a scoreable RNA matrix source
2. computes raw TRAB and TRD scores with the Phase 4 module definitions
3. colors cells by paired TRA/TRB versus all other cells
4. writes one per-GSE scatter PNG
5. consolidates suspicious-dataset summaries

Large datasets are capped at 30,000 randomly sampled cells to keep the audit
bounded and directly comparable across GSEs.
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
import gc
import hashlib
import json
import logging
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import anndata as ad
import h5py
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


matplotlib.use("Agg")


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
TABLE_ROOT = OUTPUT_ROOT / "tables" / "premerge_tcr_audit"
FIGURE_ROOT = OUTPUT_ROOT / "figures" / "premerge_tcr_audit"
LOG_ROOT = OUTPUT_ROOT / "logs"
REGISTRY_CSV = PROJECT_ROOT / "configs" / "datasets" / "integration_inputs.csv"
SCANPY_PROJECT_ROOT = PROJECT_ROOT / "analysis_26GSE_V4" / "scanpy_projects"
PER_GSE_METADATA_ROOT = PROJECT_ROOT / "downloads" / "per_gse_h5ad_with_metadata"

AUDIT_LOG = LOG_ROOT / "premerge_tcr_audit.log"
AUDIT_SUMMARY_MD = LOG_ROOT / "premerge_tcr_audit_summary.md"
MANIFEST_CSV = TABLE_ROOT / "premerge_tcr_audit_manifest.csv"
BATCH_CSV = TABLE_ROOT / "premerge_tcr_audit_batches.csv"
MATRIX_STATE_CSV = TABLE_ROOT / "premerge_tcr_audit_gse_matrix_state.csv"
SCORE_SUMMARY_CSV = TABLE_ROOT / "premerge_tcr_audit_gse_score_summary.csv"
PAIRED_COVERAGE_CSV = TABLE_ROOT / "premerge_tcr_audit_gse_paired_coverage.csv"
SUSPICIOUS_CSV = TABLE_ROOT / "premerge_tcr_audit_suspicious_ranked.csv"

SAMPLED_SCORE_DIR = TABLE_ROOT / "sampled_scores"
BATCH_RESULT_DIR = TABLE_ROOT / "batch_results"

MAX_FULL_CELLS = 60_000
MAX_SAMPLED_CELLS = 30_000
TARGET_SUM = 10_000.0
N_BINS = 25
CTRL_SIZE = 50
RANDOM_STATE = 1
FIGURE_DPI = 300
POINT_SIZE = 4.0
ALPHA = 0.55

TCR_METADATA_COLUMNS = [
    "TRA_cdr3",
    "TRB_cdr3",
    "TRG_cdr3",
    "TRD_cdr3",
    "TRA_cdr3_nt",
    "TRB_cdr3_nt",
    "TRG_cdr3_nt",
    "TRD_cdr3_nt",
    "TCRseq",
]
PAIR_REQUIRED_COLUMNS = ["TRA_cdr3", "TRB_cdr3"]

MODULE_PATTERNS = {
    "tra": re.compile(r"^TRAC$|^TRAV|^TRAJ"),
    "trb": re.compile(r"^TRBC|^TRBV|^TRBJ"),
    "trd": re.compile(r"^TRDC$|^TRDV|^TRDJ"),
}

PATH_PRIORITY = {
    "per_gse_h5ad_with_metadata": 0,
    "scanpy_project_output": 1,
    "registry_path": 2,
    "other": 99,
}


@dataclass(frozen=True)
class MatrixChoice:
    """Describe which matrix source should be scored for one dataset."""

    source_name: str
    normalization_mode: str
    x_state: str
    counts_state: str
    raw_state: str
    x_integer_like_fraction: float
    counts_integer_like_fraction: float
    raw_integer_like_fraction: float
    x_has_negative: bool
    counts_has_negative: bool
    raw_has_negative: bool
    x_max_sample: float
    counts_max_sample: float
    raw_max_sample: float


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description="Audit pre-merge per-GSE TCR score patterns.")
    parser.add_argument("--prepare-batches", action="store_true", help="Discover eligible H5ADs and write the batch manifest.")
    parser.add_argument("--run-batch", type=int, default=None, help="Run one prepared batch id.")
    parser.add_argument("--consolidate", action="store_true", help="Combine completed batch outputs into final tables.")
    parser.add_argument("--n-batches", type=int, default=3, help="Number of balanced batches for the audit.")
    return parser.parse_args()


def configure_logging() -> None:
    """Configure file and console logging."""
    LOG_ROOT.mkdir(parents=True, exist_ok=True)
    handlers = [
        logging.FileHandler(AUDIT_LOG, mode="a", encoding="utf-8"),
        logging.StreamHandler(),
    ]
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
        force=True,
    )


def ensure_output_dirs() -> None:
    """Create audit output directories."""
    for path in [TABLE_ROOT, FIGURE_ROOT, SAMPLED_SCORE_DIR, BATCH_RESULT_DIR, LOG_ROOT]:
        path.mkdir(parents=True, exist_ok=True)


def stable_seed(label: str) -> int:
    """Build a stable integer seed from one string label."""
    digest = hashlib.sha256(label.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], byteorder="big", signed=False) % (2**32 - 1)


def extract_gse_id(path: Path) -> str | None:
    """Extract the GSE id from a path or filename."""
    match = re.search(r"(GSE\d+)", str(path))
    return match.group(1) if match else None


def classify_path_source(path: Path) -> str:
    """Assign a human-readable path source label."""
    path_str = str(path)
    if "downloads/per_gse_h5ad_with_metadata" in path_str:
        return "per_gse_h5ad_with_metadata"
    if "analysis_26GSE_V4/scanpy_projects" in path_str:
        return "scanpy_project_output"
    if path == Path(path_str) and REGISTRY_CSV.exists():
        return "registry_path"
    return "other"


def registry_candidates() -> list[Path]:
    """Collect registry-listed H5AD paths when available."""
    if not REGISTRY_CSV.exists():
        return []
    registry = pd.read_csv(REGISTRY_CSV, dtype=str)
    if "h5ad_path" not in registry.columns:
        return []
    paths: list[Path] = []
    for raw_path in registry["h5ad_path"].dropna().astype(str):
        candidate = Path(raw_path)
        if not candidate.is_absolute():
            candidate = (PROJECT_ROOT / candidate).resolve()
        if candidate.exists():
            paths.append(candidate)
    return paths


def discover_candidate_paths() -> dict[str, list[Path]]:
    """Discover per-GSE H5AD candidates from the known project locations."""
    grouped: dict[str, list[Path]] = {}
    candidate_paths: list[Path] = []
    candidate_paths.extend(registry_candidates())
    if PER_GSE_METADATA_ROOT.exists():
        candidate_paths.extend(sorted(PER_GSE_METADATA_ROOT.glob("*.h5ad")))
    if SCANPY_PROJECT_ROOT.exists():
        candidate_paths.extend(sorted(SCANPY_PROJECT_ROOT.glob("GSE*/outputs/*.h5ad")))

    seen: set[Path] = set()
    for path in candidate_paths:
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        gse_id = extract_gse_id(resolved)
        if gse_id is None:
            continue
        grouped.setdefault(gse_id, []).append(resolved)
    return grouped


def obs_column_names(path: Path) -> list[str]:
    """Return the obs column names stored in one H5AD."""
    with h5py.File(path, "r") as handle:
        if "obs" not in handle:
            return []
        return sorted(handle["obs"].keys())


def has_tcr_metadata(obs_columns: Iterable[str]) -> bool:
    """Return True when the H5AD exposes any TCR-related metadata columns."""
    return any(column in obs_columns for column in TCR_METADATA_COLUMNS)


def matrix_shape(path: Path) -> tuple[int, int]:
    """Read `n_obs` and `n_vars` without loading the matrix."""
    with h5py.File(path, "r") as handle:
        obs_n = int(handle["obs"].attrs["_index"][0]) if False else len(handle["obs"].get(next(iter(handle["obs"].keys())), []))
        var_n = len(handle["var"].get(next(iter(handle["var"].keys())), []))
        if "X" in handle:
            x_node = handle["X"]
            if isinstance(x_node, h5py.Group) and "shape" in x_node.attrs:
                n_obs, n_vars = map(int, x_node.attrs["shape"])
                return n_obs, n_vars
            if isinstance(x_node, h5py.Dataset):
                return tuple(map(int, x_node.shape))
        if "_n_obs" in handle.attrs and "_n_vars" in handle.attrs:
            return int(handle.attrs["_n_obs"]), int(handle.attrs["_n_vars"])
    return obs_n, var_n


def read_dataset_values(node: h5py.Group | h5py.Dataset, sample_size: int = 200_000) -> np.ndarray:
    """Sample numeric values from one matrix-like HDF5 node."""
    if isinstance(node, h5py.Group):
        if "data" not in node:
            return np.asarray([], dtype=np.float32)
        data = node["data"]
        if data.shape[0] == 0:
            return np.asarray([], dtype=np.float32)
        step = max(1, int(math.ceil(data.shape[0] / sample_size)))
        return np.asarray(data[::step], dtype=np.float32)
    if node.size == 0:
        return np.asarray([], dtype=np.float32)
    if node.ndim == 2:
        row_count = min(node.shape[0], 256)
        col_count = min(node.shape[1], 256)
        return np.asarray(node[:row_count, :col_count], dtype=np.float32).ravel()
    step = max(1, int(math.ceil(node.shape[0] / sample_size)))
    return np.asarray(node[::step], dtype=np.float32)


def summarize_numeric_values(values: np.ndarray) -> tuple[str, float, bool, float]:
    """Classify sampled matrix values into raw-like, normalized-like, or scaled-like."""
    if values.size == 0:
        return "empty", 0.0, False, 0.0
    finite_values = values[np.isfinite(values)]
    if finite_values.size == 0:
        return "empty", 0.0, False, 0.0
    has_negative = bool(np.any(finite_values < 0))
    integer_like_fraction = float(np.mean(np.isclose(finite_values, np.round(finite_values), atol=1e-3)))
    max_sample = float(np.max(finite_values))
    if has_negative:
        return "scaled_like", integer_like_fraction, True, max_sample
    if integer_like_fraction >= 0.95 or max_sample > 30.0:
        return "raw_like", integer_like_fraction, False, max_sample
    return "normalized_like", integer_like_fraction, False, max_sample


def inspect_matrix_choice(path: Path) -> MatrixChoice:
    """Inspect X, layers/counts, and raw/X and pick the best scoring source."""
    x_state = counts_state = raw_state = "missing"
    x_integer_like_fraction = counts_integer_like_fraction = raw_integer_like_fraction = 0.0
    x_has_negative = counts_has_negative = raw_has_negative = False
    x_max_sample = counts_max_sample = raw_max_sample = 0.0

    with h5py.File(path, "r") as handle:
        if "X" in handle:
            x_state, x_integer_like_fraction, x_has_negative, x_max_sample = summarize_numeric_values(
                read_dataset_values(handle["X"])
            )
        if "layers" in handle and "counts" in handle["layers"]:
            counts_state, counts_integer_like_fraction, counts_has_negative, counts_max_sample = summarize_numeric_values(
                read_dataset_values(handle["layers"]["counts"])
            )
        if "raw" in handle and "X" in handle["raw"]:
            raw_state, raw_integer_like_fraction, raw_has_negative, raw_max_sample = summarize_numeric_values(
                read_dataset_values(handle["raw"]["X"])
            )

    if x_state == "normalized_like":
        return MatrixChoice(
            source_name="X",
            normalization_mode="already_normalized",
            x_state=x_state,
            counts_state=counts_state,
            raw_state=raw_state,
            x_integer_like_fraction=x_integer_like_fraction,
            counts_integer_like_fraction=counts_integer_like_fraction,
            raw_integer_like_fraction=raw_integer_like_fraction,
            x_has_negative=x_has_negative,
            counts_has_negative=counts_has_negative,
            raw_has_negative=raw_has_negative,
            x_max_sample=x_max_sample,
            counts_max_sample=counts_max_sample,
            raw_max_sample=raw_max_sample,
        )
    if x_state == "raw_like":
        return MatrixChoice(
            source_name="X",
            normalization_mode="normalize_from_counts",
            x_state=x_state,
            counts_state=counts_state,
            raw_state=raw_state,
            x_integer_like_fraction=x_integer_like_fraction,
            counts_integer_like_fraction=counts_integer_like_fraction,
            raw_integer_like_fraction=raw_integer_like_fraction,
            x_has_negative=x_has_negative,
            counts_has_negative=counts_has_negative,
            raw_has_negative=raw_has_negative,
            x_max_sample=x_max_sample,
            counts_max_sample=counts_max_sample,
            raw_max_sample=raw_max_sample,
        )
    if counts_state == "raw_like":
        return MatrixChoice(
            source_name="layers/counts",
            normalization_mode="normalize_from_counts",
            x_state=x_state,
            counts_state=counts_state,
            raw_state=raw_state,
            x_integer_like_fraction=x_integer_like_fraction,
            counts_integer_like_fraction=counts_integer_like_fraction,
            raw_integer_like_fraction=raw_integer_like_fraction,
            x_has_negative=x_has_negative,
            counts_has_negative=counts_has_negative,
            raw_has_negative=raw_has_negative,
            x_max_sample=x_max_sample,
            counts_max_sample=counts_max_sample,
            raw_max_sample=raw_max_sample,
        )
    if raw_state == "raw_like":
        return MatrixChoice(
            source_name="raw/X",
            normalization_mode="normalize_from_counts",
            x_state=x_state,
            counts_state=counts_state,
            raw_state=raw_state,
            x_integer_like_fraction=x_integer_like_fraction,
            counts_integer_like_fraction=counts_integer_like_fraction,
            raw_integer_like_fraction=raw_integer_like_fraction,
            x_has_negative=x_has_negative,
            counts_has_negative=counts_has_negative,
            raw_has_negative=raw_has_negative,
            x_max_sample=x_max_sample,
            counts_max_sample=counts_max_sample,
            raw_max_sample=raw_max_sample,
        )
    raise RuntimeError(
        f"No scoreable matrix source found. X={x_state}, counts={counts_state}, raw={raw_state}"
    )


def choose_best_h5ad(candidates: list[Path]) -> tuple[Path, list[str], str]:
    """Choose the best candidate H5AD path for one GSE."""
    scored_candidates: list[tuple[int, str, Path, list[str]]] = []
    for path in candidates:
        columns = obs_column_names(path)
        if not has_tcr_metadata(columns):
            continue
        source_label = classify_path_source(path)
        scored_candidates.append((PATH_PRIORITY[source_label], str(path), path, columns))
    if not scored_candidates:
        raise FileNotFoundError("No H5AD candidate with TCR-related metadata columns was found.")
    _, _, best_path, best_columns = sorted(scored_candidates)[0]
    return best_path, best_columns, classify_path_source(best_path)


def build_manifest() -> pd.DataFrame:
    """Discover and summarize all eligible pre-merge GSE H5ADs."""
    rows: list[dict[str, object]] = []
    for gse_id, candidates in sorted(discover_candidate_paths().items()):
        try:
            path, columns, source_label = choose_best_h5ad(candidates)
            n_obs, n_vars = matrix_shape(path)
            rows.append(
                {
                    "gse_id": gse_id,
                    "h5ad_path": str(path),
                    "source_label": source_label,
                    "n_total_cells": n_obs,
                    "n_vars": n_vars,
                    "obs_columns": json.dumps(columns),
                    "has_pair_columns": all(column in columns for column in PAIR_REQUIRED_COLUMNS),
                }
            )
        except Exception as exc:  # noqa: BLE001
            rows.append(
                {
                    "gse_id": gse_id,
                    "h5ad_path": "",
                    "source_label": "unresolved",
                    "n_total_cells": 0,
                    "n_vars": 0,
                    "obs_columns": "[]",
                    "has_pair_columns": False,
                    "discovery_error": str(exc),
                }
            )
    manifest = pd.DataFrame(rows).sort_values(["source_label", "gse_id"]).reset_index(drop=True)
    manifest.to_csv(MANIFEST_CSV, index=False)
    return manifest


def assign_balanced_batches(manifest: pd.DataFrame, n_batches: int) -> pd.DataFrame:
    """Assign GSEs to balanced batches using greedy cell-count bin packing."""
    eligible = manifest.loc[manifest["h5ad_path"].astype(str).ne("")].copy()
    eligible["batch_id"] = 0
    running_totals = {batch_id: 0 for batch_id in range(1, n_batches + 1)}
    for row in eligible.sort_values("n_total_cells", ascending=False).itertuples(index=False):
        batch_id = min(running_totals, key=running_totals.get)
        eligible.loc[eligible["gse_id"] == row.gse_id, "batch_id"] = batch_id
        running_totals[batch_id] += int(row.n_total_cells)
    manifest = manifest.merge(eligible[["gse_id", "batch_id"]], on="gse_id", how="left")
    manifest["batch_id"] = manifest["batch_id"].fillna(0).astype(int)
    manifest.to_csv(BATCH_CSV, index=False)
    return manifest


def read_batch_manifest() -> pd.DataFrame:
    """Load the prepared batch manifest."""
    if not BATCH_CSV.exists():
        raise FileNotFoundError("Batch manifest is missing. Run --prepare-batches first.")
    return pd.read_csv(BATCH_CSV)


def sample_indices(n_obs: int, gse_id: str) -> np.ndarray:
    """Sample cell indices with the agreed large-dataset cap."""
    if n_obs <= MAX_FULL_CELLS:
        return np.arange(n_obs, dtype=np.int64)
    rng = np.random.default_rng(stable_seed(gse_id) + RANDOM_STATE)
    return np.sort(rng.choice(n_obs, size=MAX_SAMPLED_CELLS, replace=False).astype(np.int64))


def load_subset(path: Path, sample_index: np.ndarray) -> ad.AnnData:
    """Read only the sampled subset from one backed H5AD."""
    backed = ad.read_h5ad(path, backed="r")
    try:
        subset = backed[sample_index, :].to_memory()
    finally:
        if getattr(backed, "isbacked", False):
            backed.file.close()
    return subset


def prepare_score_adata(subset: ad.AnnData, matrix_choice: MatrixChoice) -> ad.AnnData:
    """Build a temporary AnnData object in the correct expression space for scoring."""
    if matrix_choice.source_name == "X":
        matrix = subset.X.copy()
        var = subset.var.copy()
    elif matrix_choice.source_name == "layers/counts":
        if "counts" not in subset.layers:
            raise KeyError("Expected layers['counts'] but it is missing in the sampled subset.")
        matrix = subset.layers["counts"].copy()
        var = subset.var.copy()
    elif matrix_choice.source_name == "raw/X":
        if subset.raw is None:
            raise KeyError("Expected .raw.X but raw is missing in the sampled subset.")
        matrix = subset.raw.X.copy()
        var = subset.raw.var.copy()
    else:
        raise ValueError(f"Unsupported source name: {matrix_choice.source_name}")

    if sparse.issparse(matrix):
        matrix = matrix.tocsr().astype(np.float32)
    else:
        matrix = np.asarray(matrix, dtype=np.float32)
    score_adata = ad.AnnData(X=matrix, obs=subset.obs.copy(), var=var.copy())
    if matrix_choice.normalization_mode == "normalize_from_counts":
        sc.pp.normalize_total(score_adata, target_sum=TARGET_SUM, inplace=True)
        sc.pp.log1p(score_adata)
    return score_adata


def find_module_genes(var_names: pd.Index) -> dict[str, pd.Index]:
    """Find the Phase 4 module genes in one per-GSE var index."""
    modules: dict[str, pd.Index] = {}
    for name, pattern in MODULE_PATTERNS.items():
        genes = [gene for gene in var_names if pattern.match(str(gene))]
        modules[name] = pd.Index(sorted(set(genes)), dtype="string")
    modules["trab"] = pd.Index(sorted(set(modules["tra"]).union(set(modules["trb"]))), dtype="string")
    return modules


def pick_control_genes(gene_list: pd.Index, gene_pool: pd.Index, gene_means: np.ndarray, *, random_state: int) -> pd.Index:
    """Reproduce the Phase 4 control-gene binning logic."""
    rng = np.random.default_rng(random_state)
    obs_avg = pd.Series(gene_means, index=gene_pool)
    obs_avg = obs_avg[np.isfinite(obs_avg)]
    n_items = int(np.round(len(obs_avg) / (N_BINS - 1)))
    if n_items <= 0:
        raise ValueError("Invalid control-gene bin size.")
    obs_cut = obs_avg.rank(method="min") // n_items
    control_genes = pd.Index([], dtype="string")
    for cut in np.unique(obs_cut.loc[gene_list]):
        r_genes = obs_cut[obs_cut == cut].index
        if len(r_genes) == 0:
            continue
        if CTRL_SIZE < len(r_genes):
            sampled_positions = rng.choice(len(r_genes), size=CTRL_SIZE, replace=False)
            r_genes = r_genes.take(sampled_positions)
        r_genes = r_genes.difference(gene_list)
        control_genes = control_genes.union(r_genes)
    if len(control_genes) == 0:
        raise RuntimeError("No control genes found for the requested module.")
    return control_genes


def gene_means_from_matrix(matrix: sparse.spmatrix | np.ndarray) -> np.ndarray:
    """Compute per-gene mean expression on the scored matrix."""
    if sparse.issparse(matrix):
        return np.asarray(matrix.mean(axis=0)).ravel().astype(np.float32)
    return np.asarray(matrix.mean(axis=0), dtype=np.float32).ravel()


def score_matrix_columns(matrix: sparse.spmatrix | np.ndarray, gene_idx: np.ndarray, ctrl_idx: np.ndarray) -> np.ndarray:
    """Compute one score vector as module mean minus control-gene mean."""
    if sparse.issparse(matrix):
        gene_mean = np.asarray(matrix[:, gene_idx].sum(axis=1)).ravel() / float(len(gene_idx))
        ctrl_mean = np.asarray(matrix[:, ctrl_idx].sum(axis=1)).ravel() / float(len(ctrl_idx))
        return (gene_mean - ctrl_mean).astype(np.float32)
    gene_mean = np.asarray(matrix[:, gene_idx], dtype=np.float32).mean(axis=1)
    ctrl_mean = np.asarray(matrix[:, ctrl_idx], dtype=np.float32).mean(axis=1)
    return (gene_mean - ctrl_mean).astype(np.float32)


def compute_scores(score_adata: ad.AnnData, gse_id: str) -> tuple[pd.DataFrame, dict[str, int]]:
    """Compute raw TRA, TRB, TRAB, and TRD scores for one sampled dataset."""
    modules = find_module_genes(score_adata.var_names.astype(str))
    if len(modules["trab"]) == 0 or len(modules["trd"]) == 0:
        raise RuntimeError(
            f"Missing required module genes. TRAB={len(modules['trab'])}, TRD={len(modules['trd'])}"
        )

    matrix = score_adata.X
    gene_means = gene_means_from_matrix(matrix)
    scores: dict[str, np.ndarray] = {}

    for module_name in ["tra", "trb", "trab", "trd"]:
        gene_list = modules[module_name]
        if len(gene_list) == 0:
            scores[module_name] = np.full(score_adata.n_obs, np.nan, dtype=np.float32)
            continue
        controls = pick_control_genes(
            gene_list,
            score_adata.var_names,
            gene_means,
            random_state=stable_seed(f"{gse_id}|{module_name}") + RANDOM_STATE,
        )
        gene_idx = score_adata.var_names.get_indexer(gene_list)
        ctrl_idx = score_adata.var_names.get_indexer(controls)
        scores[module_name] = score_matrix_columns(matrix, gene_idx, ctrl_idx)

    frame = pd.DataFrame(
        {
            "phase4_tra_score": scores["tra"],
            "phase4_trb_score": scores["trb"],
            "phase4_trab_score": scores["trab"],
            "phase4_trd_score": scores["trd"],
        }
    )
    frame["phase4_trd_minus_trab"] = frame["phase4_trd_score"] - frame["phase4_trab_score"]
    module_sizes = {f"{name}_module_genes": int(len(genes)) for name, genes in modules.items()}
    return frame, module_sizes


def strip_series(series: pd.Series) -> pd.Series:
    """Normalize string-like obs columns for robust emptiness checks."""
    string_series = series.astype("string")
    return string_series.fillna("").str.strip()


def pairing_groups(obs: pd.DataFrame) -> tuple[pd.Series, bool]:
    """Build the paired-versus-other labels from obs metadata."""
    if not all(column in obs.columns for column in PAIR_REQUIRED_COLUMNS):
        return pd.Series(np.repeat("Other cells", obs.shape[0]), index=obs.index, dtype=object), False
    paired = strip_series(obs["TRA_cdr3"]).ne("") & strip_series(obs["TRB_cdr3"]).ne("")
    labels = pd.Series(np.where(paired, "Paired TRA/TRB", "Other cells"), index=obs.index, dtype=object)
    return labels, True


def write_scatter(gse_id: str, frame: pd.DataFrame, figure_path: Path, n_total_cells: int, sampled_mode: str) -> None:
    """Write the per-GSE TRAB-vs-TRD scatter PNG."""
    palette = {"Other cells": "#0077B6", "Paired TRA/TRB": "#D1495B"}
    fig, ax = plt.subplots(figsize=(6.6, 5.6))
    for group_name in ["Other cells", "Paired TRA/TRB"]:
        group_df = frame.loc[frame["tcr_pairing_group"] == group_name]
        if group_df.empty:
            continue
        ax.scatter(
            group_df["phase4_trab_score"],
            group_df["phase4_trd_score"],
            s=POINT_SIZE,
            linewidths=0,
            alpha=ALPHA,
            rasterized=True,
            color=palette[group_name],
            label=group_name,
        )
    paired_fraction = float(np.mean(frame["tcr_pairing_group"].eq("Paired TRA/TRB"))) if not frame.empty else 0.0
    ax.set_title(
        f"{gse_id} pre-merge raw TRAB-vs-TRD\n"
        f"{sampled_mode}; plotted={frame.shape[0]:,}; total={n_total_cells:,}; paired={paired_fraction:.1%}",
        fontsize=13,
    )
    ax.set_xlabel("Raw TRAB score", fontsize=11)
    ax.set_ylabel("Raw TRD score", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, title="Pairing context", fontsize=9, title_fontsize=10)
    fig.tight_layout()
    fig.savefig(figure_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def analyze_one_gse(row: pd.Series) -> tuple[dict[str, object], dict[str, object], dict[str, object]]:
    """Run the per-GSE audit for one manifest row."""
    gse_id = str(row["gse_id"])
    path = Path(str(row["h5ad_path"]))
    figure_path = FIGURE_ROOT / f"{gse_id}_trab_vs_trd_paired_context.png"
    sampled_score_path = SAMPLED_SCORE_DIR / f"{gse_id}_sampled_scores.csv.gz"

    logging.info("Auditing %s from %s", gse_id, path)
    matrix_choice = inspect_matrix_choice(path)

    backed = ad.read_h5ad(path, backed="r")
    try:
        n_total_cells = int(backed.n_obs)
        all_pairing_labels, pairing_complete = pairing_groups(backed.obs)
        full_paired_fraction = float(np.mean(all_pairing_labels.eq("Paired TRA/TRB")))
        sample_index = sample_indices(n_total_cells, gse_id)
        sampled_mode = "full" if n_total_cells <= MAX_FULL_CELLS else f"random_{MAX_SAMPLED_CELLS}"
        subset = backed[sample_index, :].to_memory()
    finally:
        if getattr(backed, "isbacked", False):
            backed.file.close()

    subset_pairing_labels, _ = pairing_groups(subset.obs)
    score_adata = prepare_score_adata(subset, matrix_choice)
    score_frame, module_sizes = compute_scores(score_adata, gse_id)
    score_frame["tcr_pairing_group"] = subset_pairing_labels.to_numpy()
    score_frame.to_csv(sampled_score_path, index=False, compression="gzip")
    write_scatter(gse_id, score_frame, figure_path, n_total_cells, sampled_mode)

    paired_mask = score_frame["tcr_pairing_group"].eq("Paired TRA/TRB")
    trd_gt_0p1 = score_frame["phase4_trd_score"] > 0.1
    trd_gt_0p5 = score_frame["phase4_trd_score"] > 0.5

    matrix_row = {
        "gse_id": gse_id,
        "h5ad_path": str(path),
        "source_label": row["source_label"],
        "source_name": matrix_choice.source_name,
        "normalization_mode": matrix_choice.normalization_mode,
        "x_state": matrix_choice.x_state,
        "counts_state": matrix_choice.counts_state,
        "raw_state": matrix_choice.raw_state,
        "x_integer_like_fraction": matrix_choice.x_integer_like_fraction,
        "counts_integer_like_fraction": matrix_choice.counts_integer_like_fraction,
        "raw_integer_like_fraction": matrix_choice.raw_integer_like_fraction,
        "x_has_negative": matrix_choice.x_has_negative,
        "counts_has_negative": matrix_choice.counts_has_negative,
        "raw_has_negative": matrix_choice.raw_has_negative,
        "x_max_sample": matrix_choice.x_max_sample,
        "counts_max_sample": matrix_choice.counts_max_sample,
        "raw_max_sample": matrix_choice.raw_max_sample,
        "n_total_cells": n_total_cells,
        "n_scored_cells": int(score_frame.shape[0]),
        "sampling_mode": sampled_mode,
        "pairing_complete": pairing_complete,
        "figure_path": str(figure_path),
        "sampled_score_path": str(sampled_score_path),
        **module_sizes,
        "status": "success",
    }
    score_row = {
        "gse_id": gse_id,
        "n_total_cells": n_total_cells,
        "n_scored_cells": int(score_frame.shape[0]),
        "sampling_mode": sampled_mode,
        "phase4_trab_score_mean": float(score_frame["phase4_trab_score"].mean()),
        "phase4_trd_score_mean": float(score_frame["phase4_trd_score"].mean()),
        "phase4_trd_minus_trab_mean": float(score_frame["phase4_trd_minus_trab"].mean()),
        "phase4_trab_score_median": float(score_frame["phase4_trab_score"].median()),
        "phase4_trd_score_median": float(score_frame["phase4_trd_score"].median()),
        "phase4_trd_minus_trab_median": float(score_frame["phase4_trd_minus_trab"].median()),
        "figure_path": str(figure_path),
        "sampled_score_path": str(sampled_score_path),
        **module_sizes,
        "status": "success",
    }
    coverage_row = {
        "gse_id": gse_id,
        "n_total_cells": n_total_cells,
        "n_scored_cells": int(score_frame.shape[0]),
        "pairing_complete": pairing_complete,
        "full_paired_fraction": full_paired_fraction,
        "sampled_paired_fraction": float(np.mean(paired_mask)),
        "paired_fraction_trd_gt_0p1": float(np.mean(paired_mask.loc[trd_gt_0p1])) if trd_gt_0p1.any() else np.nan,
        "paired_fraction_trd_gt_0p5": float(np.mean(paired_mask.loc[trd_gt_0p5])) if trd_gt_0p5.any() else np.nan,
        "n_trd_gt_0p1": int(trd_gt_0p1.sum()),
        "n_trd_gt_0p5": int(trd_gt_0p5.sum()),
        "n_paired_cells_sampled": int(paired_mask.sum()),
        "n_other_cells_sampled": int((~paired_mask).sum()),
        "status": "success",
    }

    del subset
    del score_adata
    gc.collect()
    return matrix_row, score_row, coverage_row


def batch_result_paths(batch_id: int) -> dict[str, Path]:
    """Return the output paths for one batch id."""
    return {
        "matrix": BATCH_RESULT_DIR / f"batch_{batch_id}_matrix_state.csv",
        "score": BATCH_RESULT_DIR / f"batch_{batch_id}_score_summary.csv",
        "coverage": BATCH_RESULT_DIR / f"batch_{batch_id}_paired_coverage.csv",
    }


def run_batch(batch_id: int) -> None:
    """Run one prepared batch."""
    manifest = read_batch_manifest()
    batch_rows = manifest.loc[manifest["batch_id"] == batch_id].copy()
    if batch_rows.empty:
        raise ValueError(f"Batch {batch_id} has no assigned GSEs.")

    matrix_rows: list[dict[str, object]] = []
    score_rows: list[dict[str, object]] = []
    coverage_rows: list[dict[str, object]] = []

    for row in batch_rows.itertuples(index=False):
        row_series = pd.Series(row._asdict())
        try:
            matrix_row, score_row, coverage_row = analyze_one_gse(row_series)
        except Exception as exc:  # noqa: BLE001
            logging.exception("Failed to audit %s", row_series["gse_id"])
            matrix_row = {
                "gse_id": row_series["gse_id"],
                "h5ad_path": row_series["h5ad_path"],
                "source_label": row_series["source_label"],
                "status": "error",
                "error": str(exc),
            }
            score_row = {
                "gse_id": row_series["gse_id"],
                "status": "error",
                "error": str(exc),
            }
            coverage_row = {
                "gse_id": row_series["gse_id"],
                "status": "error",
                "error": str(exc),
            }
        matrix_rows.append(matrix_row)
        score_rows.append(score_row)
        coverage_rows.append(coverage_row)

    output_paths = batch_result_paths(batch_id)
    pd.DataFrame(matrix_rows).sort_values("gse_id").to_csv(output_paths["matrix"], index=False)
    pd.DataFrame(score_rows).sort_values("gse_id").to_csv(output_paths["score"], index=False)
    pd.DataFrame(coverage_rows).sort_values("gse_id").to_csv(output_paths["coverage"], index=False)
    logging.info("Finished batch %s with %s datasets", batch_id, batch_rows.shape[0])


def consolidate_results() -> None:
    """Combine all batch outputs into canonical audit tables and one summary note."""
    manifest = read_batch_manifest()
    matrix_frames: list[pd.DataFrame] = []
    score_frames: list[pd.DataFrame] = []
    coverage_frames: list[pd.DataFrame] = []
    for batch_id in sorted(manifest["batch_id"].dropna().astype(int).unique()):
        if batch_id <= 0:
            continue
        paths = batch_result_paths(batch_id)
        if paths["matrix"].exists():
            matrix_frames.append(pd.read_csv(paths["matrix"]))
        if paths["score"].exists():
            score_frames.append(pd.read_csv(paths["score"]))
        if paths["coverage"].exists():
            coverage_frames.append(pd.read_csv(paths["coverage"]))

    if not matrix_frames or not score_frames or not coverage_frames:
        raise RuntimeError("At least one batch output family is missing; cannot consolidate.")

    matrix_df = pd.concat(matrix_frames, ignore_index=True).sort_values("gse_id")
    score_df = pd.concat(score_frames, ignore_index=True).sort_values("gse_id")
    coverage_df = pd.concat(coverage_frames, ignore_index=True).sort_values("gse_id")

    matrix_df.to_csv(MATRIX_STATE_CSV, index=False)
    score_df.to_csv(SCORE_SUMMARY_CSV, index=False)
    coverage_df.to_csv(PAIRED_COVERAGE_CSV, index=False)

    suspicious = score_df.merge(coverage_df, on=["gse_id", "status"], how="outer", suffixes=("_score", "_coverage"))
    suspicious = suspicious.merge(
        manifest[["gse_id", "h5ad_path", "source_label", "n_total_cells", "batch_id"]],
        on="gse_id",
        how="left",
    )
    suspicious["suspicion_score"] = (
        suspicious["full_paired_fraction"].fillna(0.0)
        + 2.0 * suspicious["paired_fraction_trd_gt_0p1"].fillna(0.0)
        + 5.0 * suspicious["paired_fraction_trd_gt_0p5"].fillna(0.0)
    )
    suspicious = suspicious.sort_values(["status", "suspicion_score"], ascending=[True, False])
    suspicious.to_csv(SUSPICIOUS_CSV, index=False)

    success_count = int(np.sum(matrix_df["status"].eq("success")))
    error_count = int(np.sum(matrix_df["status"].eq("error")))
    top_suspicious = suspicious.loc[suspicious["status"] == "success"].head(10)
    with AUDIT_SUMMARY_MD.open("w", encoding="utf-8") as handle:
        handle.write("# Pre-merge TCR audit summary\n\n")
        handle.write(f"- Eligible GSE H5ADs discovered: `{manifest.loc[manifest['batch_id'] > 0].shape[0]}`\n")
        handle.write(f"- Successful per-GSE audits: `{success_count}`\n")
        handle.write(f"- Failed per-GSE audits: `{error_count}`\n")
        handle.write("- Large datasets were scored on random 30,000-cell subsets; smaller datasets were kept full.\n")
        handle.write("- Raw TRAB-vs-TRD scatter plots are colored red for `Paired TRA/TRB` and blue for `Other cells`.\n\n")
        handle.write("## Top suspicious datasets\n\n")
        if top_suspicious.empty:
            handle.write("No successful datasets were available for suspicious ranking.\n")
        else:
            for row in top_suspicious.itertuples(index=False):
                handle.write(
                    f"- `{row.gse_id}`: suspicion score `{row.suspicion_score:.3f}`, "
                    f"full paired `{row.full_paired_fraction:.1%}`, "
                    f"paired in `TRD > 0.1` `{row.paired_fraction_trd_gt_0p1:.1%}`, "
                    f"paired in `TRD > 0.5` `{row.paired_fraction_trd_gt_0p5:.1%}`\n"
                )
        if error_count:
            handle.write("\n## Failed datasets\n\n")
            for row in matrix_df.loc[matrix_df["status"] == "error"].itertuples(index=False):
                handle.write(f"- `{row.gse_id}`: `{getattr(row, 'error', 'unknown error')}`\n")


def main() -> None:
    """CLI entrypoint."""
    args = parse_args()
    configure_logging()
    ensure_output_dirs()

    selected_modes = [args.prepare_batches, args.run_batch is not None, args.consolidate]
    if sum(bool(mode) for mode in selected_modes) != 1:
        raise SystemExit("Choose exactly one of --prepare-batches, --run-batch, or --consolidate.")

    if args.prepare_batches:
        logging.info("Preparing pre-merge TCR audit manifest and %s balanced batches", args.n_batches)
        manifest = build_manifest()
        manifest = assign_balanced_batches(manifest, args.n_batches)
        logging.info(
            "Prepared %s eligible GSEs across %s batches",
            int(np.sum(manifest["batch_id"] > 0)),
            args.n_batches,
        )
        return
    if args.run_batch is not None:
        logging.info("Running pre-merge TCR audit batch %s", args.run_batch)
        run_batch(args.run_batch)
        return
    logging.info("Consolidating completed pre-merge TCR audit batches")
    consolidate_results()


if __name__ == "__main__":
    main()
