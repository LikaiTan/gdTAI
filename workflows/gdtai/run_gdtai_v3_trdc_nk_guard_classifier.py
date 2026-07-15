#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Train and externally test gdTAI v3 TRDC/NK guard candidates.

This downstream lane revises the previous GSE144469-holdout setup:

- GSE144469 primary gold cells are eligible for atlas train/tune splits.
- The independent external H5AD is never used for fitting, feature selection,
  threshold tuning, or model selection.
- External features are rebuilt only from ``layers["counts"]`` as
  log1p(counts per 10,000); the script fails if that layer is missing.
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
import gzip
import html
import json
import logging
import math
import pickle
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier, HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.metrics import brier_score_loss, confusion_matrix
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

try:
    from lightgbm import LGBMClassifier

    LIGHTGBM_AVAILABLE = True
except Exception:
    LGBMClassifier = None
    LIGHTGBM_AVAILABLE = False

from run_gdt_deg_tcr_classifier_training import metric_dict, safe_div
from run_gdt_gse144469_holdout_tcrgene_classifier import (
    EXTRA_TRAB_HOLDOUT_SOURCE,
    GDT2020_HOLDOUT_TISSUE,
    GDT2020_SOURCE,
    RANDOM_SEED,
    SUBOPTIMAL_SOURCE,
    TCR_PREFIXES,
    build_obs_metadata,
    load_original_trd_trab_threshold,
    tcr_priority,
)
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


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
STATIC_DIR = PROJECT_ROOT / "gdT_prediction"
OUT_PREFIX = "gdtai_v3_trdc_nk_guard"
TABLE_DIR = OUTPUT_ROOT / "tables" / "gdT_prediction" / OUT_PREFIX
FIGURE_DIR = OUTPUT_ROOT / "figures" / "gdT_prediction" / OUT_PREFIX
LOG_DIR = OUTPUT_ROOT / "logs" / "gdT_prediction" / OUT_PREFIX
MODEL_DIR = OUTPUT_ROOT / "models" / "gdT_prediction_classifier" / OUT_PREFIX
PROMOTED_DIR = OUTPUT_ROOT / "models" / "gdT_prediction_classifier" / "gdTAI_v3.0"
REPORT_MD = LOG_DIR / "gdtai_v3_trdc_nk_guard_report.md"
REPORT_HTML = STATIC_DIR / "gdtai_v3_trdc_nk_guard_report.html"
REPORT_PDF = STATIC_DIR / "gdtai_v3_trdc_nk_guard_report.pdf"
RUN_LOG = LOG_DIR / "run.log"
SUMMARY_JSON = LOG_DIR / "gdtai_v3_trdc_nk_guard_summary.json"

EXTERNAL_H5AD = Path("/home/tanlikai/databank/owndata/singlecell/data/results/phase4_final_annotated.h5ad")
V2_WRAPPER = OUTPUT_ROOT / "models" / "gdT_prediction_classifier" / "gdTAI_v2.0" / "gdTAI_v2_model.pkl"

TARGET_SUM = 10_000.0
MAX_FEATURE_GENES = 340
FULL_CHUNK = 50_000
DEFAULT_MAX_PRIMARY_NEGATIVE_TRAIN = 260_000
DEFAULT_MAX_HARD_NEGATIVE_TRAIN = 220_000
DEFAULT_MAX_WEAK_TRDC_PREFILTER = 180_000
DEFAULT_BOOTSTRAP_ITERATIONS = 12

NK_CONTROL_GENES = [
    "NKG7",
    "GNLY",
    "PRF1",
    "GZMB",
    "GZMH",
    "KLRD1",
    "KLRF1",
    "FCGR3A",
    "NCAM1",
    "TYROBP",
    "FCER1G",
    "CST7",
]
CD3_GENES = ["CD3D", "CD3E", "CD3G"]
PENALTY_CONTROL_GENES = ["FOXP3", "CD4", *CD3_GENES]
B_CELL_MARKERS = ["MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN"]
MYELOID_MARKERS = ["LYZ", "LST1", "S100A8", "S100A9", "FCN1", "CST3"]
EXTRA_CONTROL_GENES = [*NK_CONTROL_GENES, *B_CELL_MARKERS, *MYELOID_MARKERS]

ENGINEERED_FEATURES = [
    "any_TRDV",
    "any_TRDJ",
    "any_TRG",
    "any_ab_TCR_gene",
    "TRDC_only",
    "TRDC_plus_TRDV",
    "TRDC_plus_TRDJ",
    "CD3_score",
    "NK_score",
    "gdT_TCR_score",
    "abT_TCR_score",
    "NK_minus_CD3_score",
    "TRDC_log1p",
    "TRDV_score",
    "TRDJ_score",
    "TRG_score",
]

V2_BASELINE_STRATEGIES = ["v2_high_f1", "v2_high_purity"]
REFERENCE_STRATEGIES = [*V2_BASELINE_STRATEGIES, "original_TRD_minus_TRAB"]


@dataclass
class FeatureSpec:
    gene_names: list[str]
    gene_indices: np.ndarray
    gene_feature_names: list[str]
    engineered_feature_names: list[str]
    model_feature_names: list[str]
    gene_to_col: dict[str, int]
    engineered_to_col: dict[str, int]


@dataclass
class SplitBundle:
    train_idx: np.ndarray
    tune_idx: np.ndarray
    hard_negative_train_idx: np.ndarray
    hard_negative_tune_idx: np.ndarray
    weak_trdc_prefilter_idx: np.ndarray
    validation_idx: np.ndarray
    validation_gdt2020_idx: np.ndarray
    validation_extra_trab_idx: np.ndarray
    split_overall: pd.DataFrame
    split_by_source: pd.DataFrame
    annotation: np.ndarray
    nk_guard_mask: np.ndarray
    tcrab_guard_mask: np.ndarray


@dataclass
class CandidateResult:
    name: str
    kind: str
    model_object: Any
    threshold: float
    tune_metrics: dict[str, Any]
    tune_score: np.ndarray
    tune_pred: np.ndarray
    notes: str
    promotable: bool = True
    iteration_round: int = 0
    model_family: str = ""
    threshold_policy: str = "f1"


class PlattCalibratedModel:
    """Small pickle-safe Platt calibrator around a fitted sklearn model."""

    def __init__(self, base_model: Any, calibrator: LogisticRegression):
        self.base_model = base_model
        self.calibrator = calibrator

    @staticmethod
    def _logit(prob: np.ndarray) -> np.ndarray:
        clipped = np.clip(prob.astype(np.float64), 1e-6, 1.0 - 1e-6)
        return np.log(clipped / (1.0 - clipped)).reshape(-1, 1)

    def predict_proba(self, x: np.ndarray) -> np.ndarray:
        raw = self.base_model.predict_proba(x)[:, 1]
        calibrated = self.calibrator.predict_proba(self._logit(raw))[:, 1].astype(np.float32)
        return np.column_stack([1.0 - calibrated, calibrated])


class ConditionalGateModel:
    """Wrap another probabilistic model with conditional TRDC-only risk gates."""

    def __init__(
        self,
        base_model: Any,
        *,
        base_cols: list[int] | None,
        engineered_to_col: dict[str, int],
        mode: str,
    ):
        self.base_model = base_model
        self.base_cols = base_cols
        self.engineered_to_col = engineered_to_col
        self.mode = mode

    def _col(self, name: str) -> int:
        return int(self.engineered_to_col[name])

    def predict_proba(self, x: np.ndarray) -> np.ndarray:
        model_x = x[:, self.base_cols] if self.base_cols is not None else x
        score = self.base_model.predict_proba(model_x)[:, 1].astype(np.float32)
        trdc_only = x[:, self._col("TRDC_only")] > 0.5
        any_trdv = x[:, self._col("any_TRDV")] > 0.5
        any_trdj = x[:, self._col("any_TRDJ")] > 0.5
        nk_minus_cd3 = x[:, self._col("NK_minus_CD3_score")]
        gd_tcr = x[:, self._col("gdT_TCR_score")]
        ab_tcr = x[:, self._col("abT_TCR_score")]
        cd3 = x[:, self._col("CD3_score")]
        risk = trdc_only & (~any_trdv) & (~any_trdj) & (nk_minus_cd3 > 0.35)
        if self.mode == "two_stage_tcell_gate":
            risk |= (cd3 < 0.20) & (nk_minus_cd3 > 0.25) & (gd_tcr < 0.55)
        risk |= trdc_only & (ab_tcr > gd_tcr + 0.40) & (nk_minus_cd3 > 0.15)
        gated = score.copy()
        gated[risk] = np.minimum(gated[risk], 0.02)
        return np.column_stack([1.0 - gated, gated])


class StrictTRDVGateModel:
    """Post-gate a probabilistic model by requiring TRDV expression for TRDC+ cells."""

    def __init__(
        self,
        base_model: Any,
        *,
        base_cols: list[int] | None,
        engineered_to_col: dict[str, int],
        block_score: float = 0.001,
    ):
        self.base_model = base_model
        self.base_cols = base_cols
        self.engineered_to_col = engineered_to_col
        self.block_score = float(block_score)

    def _col(self, name: str) -> int:
        return int(self.engineered_to_col[name])

    def predict_proba(self, x: np.ndarray) -> np.ndarray:
        model_x = x[:, self.base_cols] if self.base_cols is not None else x
        score = self.base_model.predict_proba(model_x)[:, 1].astype(np.float32)
        trdc = x[:, self._col("TRDC_log1p")] > 0
        trdv = x[:, self._col("any_TRDV")] > 0.5
        block = trdc & (~trdv)
        gated = score.copy()
        gated[block] = np.minimum(gated[block], self.block_score)
        return np.column_stack([1.0 - gated, gated])


class TRDVRescuePostGateModel:
    """Use a low-FP base model, then rescue only TRDV+ cells with strong v2 support."""

    def __init__(
        self,
        base_model: Any,
        rescue_model: Any,
        *,
        base_cols: list[int] | None,
        rescue_cols: list[int] | None,
        engineered_to_col: dict[str, int],
        rescue_threshold: float,
        rescue_floor: float = 0.97,
        block_score: float = 0.001,
    ):
        self.base_model = base_model
        self.rescue_model = rescue_model
        self.base_cols = base_cols
        self.rescue_cols = rescue_cols
        self.engineered_to_col = engineered_to_col
        self.rescue_threshold = float(rescue_threshold)
        self.rescue_floor = float(rescue_floor)
        self.block_score = float(block_score)

    def _col(self, name: str) -> int:
        return int(self.engineered_to_col[name])

    def predict_proba(self, x: np.ndarray) -> np.ndarray:
        base_x = x[:, self.base_cols] if self.base_cols is not None else x
        rescue_x = x[:, self.rescue_cols] if self.rescue_cols is not None else x
        score = self.base_model.predict_proba(base_x)[:, 1].astype(np.float32)
        rescue_score = self.rescue_model.predict_proba(rescue_x)[:, 1].astype(np.float32)
        trdc = x[:, self._col("TRDC_log1p")] > 0
        trdv = x[:, self._col("any_TRDV")] > 0.5
        rescue = trdv & (rescue_score >= self.rescue_threshold)
        gated = score.copy()
        gated[rescue] = np.maximum(gated[rescue], self.rescue_floor)
        block = trdc & (~trdv)
        gated[block] = np.minimum(gated[block], self.block_score)
        return np.column_stack([1.0 - gated, gated])


# Saved model payloads should be loadable by importing this script as a module,
# even when the training run executed it as __main__.
sys.modules.setdefault("run_gdtai_v3_trdc_nk_guard_classifier", sys.modules[__name__])
PlattCalibratedModel.__module__ = "run_gdtai_v3_trdc_nk_guard_classifier"
ConditionalGateModel.__module__ = "run_gdtai_v3_trdc_nk_guard_classifier"
StrictTRDVGateModel.__module__ = "run_gdtai_v3_trdc_nk_guard_classifier"
TRDVRescuePostGateModel.__module__ = "run_gdtai_v3_trdc_nk_guard_classifier"

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train and test gdTAI v3 TRDC/NK guard candidates.")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--external-h5ad", type=Path, default=EXTERNAL_H5AD)
    parser.add_argument("--v2-model-pkl", type=Path, default=V2_WRAPPER)
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    parser.add_argument("--chunk-size", type=int, default=FULL_CHUNK)
    parser.add_argument("--max-primary-negative-train", type=int, default=DEFAULT_MAX_PRIMARY_NEGATIVE_TRAIN)
    parser.add_argument("--max-hard-negative-train", type=int, default=DEFAULT_MAX_HARD_NEGATIVE_TRAIN)
    parser.add_argument("--max-weak-trdc-prefilter", type=int, default=DEFAULT_MAX_WEAK_TRDC_PREFILTER)
    parser.add_argument("--bootstrap-iterations", type=int, default=DEFAULT_BOOTSTRAP_ITERATIONS)
    parser.add_argument("--bootstrap-sample-cells", type=int, default=60_000)
    parser.add_argument("--target-recall", type=float, default=0.80, help="Target recall for target-aware candidate thresholds.")
    parser.add_argument("--target-recall-margin", type=float, default=0.005, help="Small full-atlas recall cushion used when choosing among target-satisfying candidates.")
    parser.add_argument("--max-estimated-fp-fraction", type=float, default=0.05, help="Maximum estimated FP fraction among predictions.")
    parser.add_argument("--fp-comparator-strategy", default="v2_high_purity", help="Full-atlas strategy whose estimated FP fraction is the not-worse-than comparator.")
    parser.add_argument("--max-trdc-trdv-minus-fraction", type=float, default=0.05, help="Maximum TRDC+TRDV- fraction among full-atlas predicted gdT cells.")
    parser.add_argument("--min-iteration-rounds", type=int, default=5, help="Evaluate at least this many candidate rounds unless the target is reached earlier.")
    parser.add_argument("--skip-lightgbm", action="store_true", help="Skip the LightGBM candidate family.")
    parser.add_argument("--skip-extra-trees", action="store_true", help="Skip the ExtraTrees candidate family if --use-extra-trees is set.")
    parser.add_argument("--use-extra-trees", action="store_true", help="Opt in to the slower ExtraTrees candidate family.")
    parser.add_argument("--skip-hist-gradient", action="store_true", help="Skip the sklearn histogram-gradient candidate family.")
    parser.add_argument("--max-external-cells", type=int, default=None, help="Smoke-test cap only; default evaluates all external cells.")
    parser.add_argument("--extra-trab-validation-source", default=EXTRA_TRAB_HOLDOUT_SOURCE)
    parser.add_argument("--gdt2020-validation-tissue", default=GDT2020_HOLDOUT_TISSUE)
    parser.add_argument("--no-full-atlas", action="store_true", help="Skip final whole-atlas inference; intended only for smoke tests.")
    parser.add_argument("--max-full-atlas-cells", type=int, default=None, help="Smoke-test cap only; default applies to all input H5AD cells.")
    parser.add_argument("--skip-bootstrap", action="store_true")
    parser.add_argument("--no-pdf", action="store_true")
    return parser.parse_args()


def ensure_dirs() -> None:
    for path in [TABLE_DIR, FIGURE_DIR, LOG_DIR, MODEL_DIR, STATIC_DIR, PROMOTED_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def setup_logging() -> None:
    ensure_dirs()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(RUN_LOG, mode="w", encoding="utf-8"), logging.StreamHandler()],
        force=True,
    )


def load_pickle(path: Path) -> Any:
    with path.open("rb") as handle:
        return pickle.load(handle)


def valid_bool(values: np.ndarray | pd.Series) -> np.ndarray:
    series = pd.Series(values, copy=False).astype("string").fillna("").str.strip().str.lower()
    return series.isin({"true", "1", "yes", "y", "t"}).to_numpy(dtype=bool)


def n_obs_from_handle(handle: h5py.File) -> int:
    if "_index" in handle["obs"]:
        return int(handle["obs"]["_index"].shape[0])
    if "X" in handle:
        return h5ad_shape(handle, "X")[0]
    if "layers" in handle and "counts" in handle["layers"]:
        return h5ad_shape(handle, "layers/counts")[0]
    raise KeyError("Cannot infer n_obs from H5AD; missing obs[_index], X, and layers[counts].")


def obs_index_values(handle: h5py.File) -> np.ndarray:
    n_obs = n_obs_from_handle(handle)
    for key in ["_index", "cell_id", "barcode", "barcode_original"]:
        if key in handle["obs"]:
            return read_string_like(handle, key)
    return np.asarray([str(i) for i in range(n_obs)], dtype=object)


def read_bool_like(handle: h5py.File, key: str, *, default: bool = False) -> np.ndarray:
    n_obs = n_obs_from_handle(handle)
    if key not in handle["obs"]:
        return np.full(n_obs, default, dtype=bool)
    values = read_obs_column(handle, key)
    if np.asarray(values).dtype == bool:
        return np.asarray(values, dtype=bool)
    return valid_bool(values)


def read_string_like(handle: h5py.File, key: str, *, default: str = "") -> np.ndarray:
    n_obs = n_obs_from_handle(handle)
    if key not in handle["obs"]:
        return np.full(n_obs, default, dtype=object)
    return pd.Series(read_obs_column(handle, key), copy=False).astype("string").fillna(default).astype(str).to_numpy(dtype=object)


def h5ad_shape(handle: h5py.File, matrix_path: str) -> tuple[int, int]:
    obj = get_h5_object(handle, matrix_path)
    if isinstance(obj, h5py.Dataset):
        shape = obj.shape
    else:
        shape = obj.attrs.get("shape")
    if shape is None:
        raise TypeError(f"Cannot infer shape for matrix path `{matrix_path}`")
    return int(shape[0]), int(shape[1])


def get_h5_object(handle: h5py.File, matrix_path: str) -> h5py.Dataset | h5py.Group:
    obj: h5py.Dataset | h5py.Group = handle
    for part in matrix_path.split("/"):
        if not part:
            continue
        obj = obj[part]
    return obj


def matrix_encoding(handle: h5py.File, matrix_path: str) -> str:
    obj = get_h5_object(handle, matrix_path)
    if isinstance(obj, h5py.Dataset):
        return "dense"
    encoding = obj.attrs.get("encoding-type", "")
    if isinstance(encoding, bytes):
        encoding = encoding.decode("utf-8")
    if not encoding and {"data", "indices", "indptr"}.issubset(obj.keys()):
        encoding = "csr_matrix"
    return str(encoding)


def build_feature_spec(handle: h5py.File, max_feature_genes: int = MAX_FEATURE_GENES) -> FeatureSpec:
    var_names = pd.Index(read_string_dataset(handle["var"]["_index"]), dtype="string").astype(str).tolist()
    gene_to_idx = {gene: idx for idx, gene in enumerate(var_names)}
    tcr_genes = sorted([gene for gene in var_names if gene.upper().startswith(TCR_PREFIXES)], key=tcr_priority)
    controls = [gene for gene in [*PENALTY_CONTROL_GENES, *EXTRA_CONTROL_GENES] if gene in gene_to_idx]
    selected: list[str] = []
    for gene in [*tcr_genes, *controls]:
        if gene not in selected:
            selected.append(gene)
    if len(selected) > max_feature_genes:
        selected = selected[:max_feature_genes]
    gene_feature_names = [f"{gene}_log1p_cp10k" for gene in selected]
    engineered_to_col = {name: len(gene_feature_names) + i for i, name in enumerate(ENGINEERED_FEATURES)}
    spec = FeatureSpec(
        gene_names=selected,
        gene_indices=np.asarray([gene_to_idx[gene] for gene in selected], dtype=np.int32),
        gene_feature_names=gene_feature_names,
        engineered_feature_names=list(ENGINEERED_FEATURES),
        model_feature_names=[*gene_feature_names, *ENGINEERED_FEATURES],
        gene_to_col={gene: i for i, gene in enumerate(selected)},
        engineered_to_col=engineered_to_col,
    )
    pd.DataFrame(
        {
            "feature": spec.model_feature_names,
            "feature_type": ["gene_log1p_cp10k"] * len(spec.gene_feature_names) + ["engineered"] * len(spec.engineered_feature_names),
            "gene": spec.gene_names + [""] * len(spec.engineered_feature_names),
            "feature_index": np.arange(len(spec.model_feature_names), dtype=int),
        }
    ).to_csv(TABLE_DIR / "feature_manifest.csv", index=False)
    return spec


def selected_gene_mapping(var_names: list[str], gene_names: list[str]) -> pd.DataFrame:
    lookup = {gene: i for i, gene in enumerate(var_names)}
    rows = []
    for i, gene in enumerate(gene_names):
        idx = lookup.get(gene)
        rows.append({"feature_index": i, "gene": gene, "available_in_h5ad": idx is not None, "h5ad_gene_index": "" if idx is None else idx})
    return pd.DataFrame(rows)


def extract_gene_features(
    handle: h5py.File,
    matrix_path: str,
    rows: np.ndarray,
    spec: FeatureSpec,
    *,
    label: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = np.asarray(rows, dtype=np.int64)
    features = np.zeros((rows.size, len(spec.gene_names)), dtype=np.float32)
    row_sum = np.zeros(rows.size, dtype=np.float32)
    n_detected = np.zeros(rows.size, dtype=np.int32)
    matrix = get_h5_object(handle, matrix_path)
    encoding = matrix_encoding(handle, matrix_path)
    if encoding == "dense":
        selected = spec.gene_indices.astype(np.int64)
        for start in range(0, rows.size, 5_000):
            end = min(start + 5_000, rows.size)
            dense = matrix[rows[start:end], :]
            totals = dense.sum(axis=1).astype(np.float32)
            row_sum[start:end] = totals
            n_detected[start:end] = (dense > 0).sum(axis=1).astype(np.int32)
            if selected.size:
                values = dense[:, selected].astype(np.float32, copy=False)
                with np.errstate(divide="ignore", invalid="ignore"):
                    transformed = np.log1p(values * (TARGET_SUM / totals[:, None]))
                transformed[~np.isfinite(transformed)] = 0.0
                features[start:end, :] = transformed.astype(np.float32, copy=False)
        return features, row_sum, n_detected

    if encoding != "csr_matrix":
        raise TypeError(f"Matrix path `{matrix_path}` must be CSR or dense; found `{encoding}`.")
    gene_to_col = {int(gene_idx): col for col, gene_idx in enumerate(spec.gene_indices)}
    selected = np.asarray(sorted(gene_to_col), dtype=np.int64)
    indptr = matrix["indptr"]
    indices = matrix["indices"]
    data = matrix["data"]
    for out_i, obs_i in enumerate(rows):
        start = int(indptr[obs_i])
        end = int(indptr[obs_i + 1])
        if end <= start:
            continue
        row_indices = indices[start:end].astype(np.int64, copy=False)
        row_data = data[start:end].astype(np.float32, copy=False)
        total = float(np.sum(row_data, dtype=np.float64))
        row_sum[out_i] = total
        n_detected[out_i] = int(np.count_nonzero(row_data > 0))
        if total <= 0 or selected.size == 0:
            continue
        present = np.isin(row_indices, selected, assume_unique=False)
        if present.any():
            values = np.log1p(row_data[present] * (TARGET_SUM / total)).astype(np.float32, copy=False)
            for gene_idx, value in zip(row_indices[present], values):
                features[out_i, gene_to_col[int(gene_idx)]] = value
        if out_i and out_i % 50_000 == 0:
            logging.info("Extracted %s features for %s / %s rows", label, f"{out_i:,}", f"{rows.size:,}")
    return features, row_sum, n_detected


def gene_cols(spec: FeatureSpec, names: list[str] | tuple[str, ...]) -> list[int]:
    return [spec.gene_to_col[name] for name in names if name in spec.gene_to_col]


def prefix_cols(spec: FeatureSpec, prefixes: tuple[str, ...]) -> list[int]:
    return [i for gene, i in spec.gene_to_col.items() if gene.upper().startswith(prefixes)]


def score_cols(x_gene: np.ndarray, cols: list[int]) -> np.ndarray:
    if not cols:
        return np.zeros(x_gene.shape[0], dtype=np.float32)
    return x_gene[:, cols].mean(axis=1).astype(np.float32)


def sum_cols(x_gene: np.ndarray, cols: list[int]) -> np.ndarray:
    if not cols:
        return np.zeros(x_gene.shape[0], dtype=np.float32)
    return x_gene[:, cols].sum(axis=1).astype(np.float32)


def append_engineered_features(x_gene: np.ndarray, spec: FeatureSpec) -> np.ndarray:
    trdc = x_gene[:, spec.gene_to_col["TRDC"]] if "TRDC" in spec.gene_to_col else np.zeros(x_gene.shape[0], dtype=np.float32)
    trdv = sum_cols(x_gene, prefix_cols(spec, ("TRDV",)))
    trdj = sum_cols(x_gene, prefix_cols(spec, ("TRDJ", "TRDD")))
    trg = sum_cols(x_gene, prefix_cols(spec, ("TRGV", "TRGJ", "TRGC")))
    ab = sum_cols(x_gene, prefix_cols(spec, ("TRAV", "TRAJ", "TRAC", "TRBV", "TRBJ", "TRBC")))
    cd3 = score_cols(x_gene, gene_cols(spec, CD3_GENES))
    nk = score_cols(x_gene, gene_cols(spec, NK_CONTROL_GENES))
    gdt_tcr = trdc + trdv + trdj + trg
    engineered = np.column_stack(
        [
            trdv > 0,
            trdj > 0,
            trg > 0,
            ab > 0,
            (trdc > 0) & (trdv <= 0) & (trdj <= 0),
            (trdc > 0) & (trdv > 0),
            (trdc > 0) & (trdj > 0),
            cd3,
            nk,
            gdt_tcr,
            ab,
            nk - cd3,
            trdc,
            trdv,
            trdj,
            trg,
        ]
    ).astype(np.float32, copy=False)
    return np.column_stack([x_gene, engineered]).astype(np.float32, copy=False)


def trdc_trdv_quadrant(x: np.ndarray, spec: FeatureSpec) -> np.ndarray:
    trdc = x[:, spec.engineered_to_col["TRDC_log1p"]] > 0
    trdv = x[:, spec.engineered_to_col["any_TRDV"]] > 0.5
    out = np.full(x.shape[0], "TRDC-TRDV-", dtype=object)
    out[trdc & (~trdv)] = "TRDC+TRDV-"
    out[(~trdc) & trdv] = "TRDC-TRDV+"
    out[trdc & trdv] = "TRDC+TRDV+"
    return out


def split_by_source_label(indices: np.ndarray, labels: np.ndarray, source: np.ndarray, seed: int) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    rows: list[dict[str, Any]] = []
    train_parts: list[np.ndarray] = []
    tune_parts: list[np.ndarray] = []
    df = pd.DataFrame({"obs_index": indices, "source_gse_id": source[indices], "label": labels[indices]})
    for (src, label), group in df.groupby(["source_gse_id", "label"], sort=True):
        idx = group["obs_index"].to_numpy(dtype=np.int64)
        rng.shuffle(idx)
        n_tune = 0 if idx.size < 5 else max(1, int(round(idx.size * 0.20)))
        if n_tune:
            tune_parts.append(idx[:n_tune])
            train_parts.append(idx[n_tune:])
        else:
            train_parts.append(idx)
        rows.append({"source_gse_id": src, "label": label, "n_cells": int(idx.size), "train": int(idx.size - n_tune), "tune": int(n_tune)})
    train = np.concatenate(train_parts).astype(np.int64) if train_parts else np.asarray([], dtype=np.int64)
    tune = np.concatenate(tune_parts).astype(np.int64) if tune_parts else np.asarray([], dtype=np.int64)
    return train, tune, pd.DataFrame(rows)


def make_splits(handle: h5py.File, args: argparse.Namespace) -> tuple[Any, SplitBundle]:
    obs = build_obs_metadata(handle)
    n_obs = obs.source.size
    annotation = clean_group_values(read_obs_column(handle, "simple_annotation_plus6")) if "simple_annotation_plus6" in handle["obs"] else np.full(n_obs, "unknown", dtype=object)
    annotation_upper = pd.Series(annotation, copy=False).astype(str).str.upper().to_numpy(dtype=object)
    nk_annotation_mask = annotation_upper == "NK_CELL"
    ab_tcr_mask = obs.has_TRA_TRB_paired | obs.has_any_ab_tcr
    tissue_lower = pd.Series(obs.tissue, copy=False).astype(str).str.lower().to_numpy(dtype=object)
    validation_gdt2020_idx = np.flatnonzero(
        (obs.source == GDT2020_SOURCE)
        & (obs.class_code == 2)
        & (tissue_lower == str(args.gdt2020_validation_tissue).lower())
        & (obs.source != SUBOPTIMAL_SOURCE)
    ).astype(np.int64)
    validation_extra_trab_idx = np.flatnonzero(
        (obs.source == args.extra_trab_validation_source)
        & obs.has_TRA_TRB_paired
        & (~obs.corrected_has_any_gd_tcr)
        & (~nk_annotation_mask)
        & (obs.source != SUBOPTIMAL_SOURCE)
    ).astype(np.int64)
    validation_idx = np.unique(np.concatenate([validation_gdt2020_idx, validation_extra_trab_idx])).astype(np.int64)
    validation_mask = np.zeros(n_obs, dtype=bool)
    validation_mask[validation_idx] = True
    validation_source_mask = obs.source == args.extra_trab_validation_source
    if validation_gdt2020_idx.size == 0:
        raise RuntimeError(f"No {GDT2020_SOURCE} gdT validation cells found with tissue `{args.gdt2020_validation_tissue}`")
    if validation_extra_trab_idx.size == 0:
        raise RuntimeError(f"No paired TRA/TRB no-gdTCR validation cells found in `{args.extra_trab_validation_source}` after NK-overlap exclusion")

    excluded_nk_tcrab_negative_mask = (
        nk_annotation_mask
        & ab_tcr_mask
        & (obs.class_code != 2)
        & (~obs.corrected_has_any_gd_tcr)
        & (obs.source != SUBOPTIMAL_SOURCE)
    )
    primary = (
        ((obs.class_code == 2) | ((obs.class_code == 1) & (~excluded_nk_tcrab_negative_mask)))
        & (~validation_mask)
        & (~validation_source_mask)
        & (obs.source != SUBOPTIMAL_SOURCE)
    )
    primary_idx = np.flatnonzero(primary).astype(np.int64)
    primary_label = np.where(obs.class_code == 2, "gdT_gold", "abT_gold")
    train_idx, tune_idx, primary_split = split_by_source_label(primary_idx, primary_label, obs.source, args.seed)

    nk_guard_mask = (
        nk_annotation_mask
        & (obs.class_code != 2)
        & (~obs.silver_mask)
        & (~obs.corrected_has_any_gd_tcr)
    )
    tcrab_guard_mask = (
        ab_tcr_mask
        & (obs.class_code != 2)
        & (~obs.corrected_has_any_gd_tcr)
    )
    nk_tcrab_overlap_mask = nk_guard_mask & tcrab_guard_mask
    nk_non_tcrab_guard_mask = nk_guard_mask & (~tcrab_guard_mask)
    tcrab_non_nk_guard_mask = tcrab_guard_mask & (~nk_guard_mask)
    # Primary abT/gdT rows are already split by source and label. Extra guard
    # negatives are drawn only from non-primary cells to avoid train/tune leakage
    # when an abT row also satisfies the paired-TCRAB guard definition.
    # User rule: do not use NK+TCRAB-overlap cells as hard negatives. NK hard
    # negatives are expression-filtered to TRDC+TRDV- after feature extraction.
    hard_pool = np.flatnonzero(
        (nk_non_tcrab_guard_mask | tcrab_non_nk_guard_mask)
        & (~primary)
        & (~validation_mask)
        & (~validation_source_mask)
        & (obs.source != SUBOPTIMAL_SOURCE)
    ).astype(np.int64)
    hard_labels = np.full(n_obs, "hard_negative", dtype=object)
    hard_labels[nk_non_tcrab_guard_mask] = "NK_TRDC_TRDV_expression_prefilter"
    hard_labels[tcrab_non_nk_guard_mask] = "TCRAB_guard_negative"
    hard_train, hard_tune, hard_split = split_by_source_label(hard_pool, hard_labels, obs.source, args.seed + 1)
    rng = np.random.default_rng(args.seed + 2)
    max_hard_train_prefeature = max(args.max_hard_negative_train * 2, args.max_hard_negative_train)
    max_hard_tune_prefeature = max(args.max_hard_negative_train // 2, 50_000)
    if hard_train.size > max_hard_train_prefeature:
        hard_train = rng.choice(hard_train, size=max_hard_train_prefeature, replace=False).astype(np.int64)
    if hard_tune.size > max_hard_tune_prefeature:
        hard_tune = rng.choice(hard_tune, size=max_hard_tune_prefeature, replace=False).astype(np.int64)


    phase4_trd = read_float_obs(handle, "phase4_trd_score") if "phase4_trd_score" in handle["obs"] else np.zeros(n_obs, dtype=np.float32)
    phase4_trab = read_float_obs(handle, "phase4_trab_score") if "phase4_trab_score" in handle["obs"] else np.zeros(n_obs, dtype=np.float32)
    weak_pool = np.flatnonzero(
        (obs.class_code == 0)
        & (~obs.corrected_has_any_gd_tcr)
        & (phase4_trd > np.nanpercentile(phase4_trd[np.isfinite(phase4_trd)], 75))
        & (phase4_trab >= np.nanpercentile(phase4_trab[np.isfinite(phase4_trab)], 40))
        & (~validation_mask)
        & (~validation_source_mask)
        & (obs.source != SUBOPTIMAL_SOURCE)
    ).astype(np.int64)
    if weak_pool.size > args.max_weak_trdc_prefilter:
        weak_pool = rng.choice(weak_pool, size=args.max_weak_trdc_prefilter, replace=False).astype(np.int64)

    split_by_source = pd.concat([primary_split, hard_split], ignore_index=True)
    split_by_source.to_csv(TABLE_DIR / "split_by_source_label.csv", index=False)
    overall_rows = []
    for name, idx in [
        ("train_primary", train_idx),
        ("tune_primary", tune_idx),
        ("train_hard_negative_prefilter", hard_train),
        ("tune_hard_negative_prefilter", hard_tune),
        ("weak_TRDC_prefilter", weak_pool),
        (f"validation_gdT_{GDT2020_SOURCE}_{str(args.gdt2020_validation_tissue).replace(' ', '_')}", validation_gdt2020_idx),
        (f"validation_abT_{args.extra_trab_validation_source}_paired_TCRAB_no_gdTCR", validation_extra_trab_idx),
        ("validation_combined_non_GSE144469", validation_idx),
        ("GSE144469_primary_in_train_tune", primary_idx[obs.source[primary_idx] == "GSE144469"]),
        ("excluded_NK_TCRAB_overlap_negative", np.flatnonzero(excluded_nk_tcrab_negative_mask).astype(np.int64)),
        (f"sensitivity_excluded_{SUBOPTIMAL_SOURCE}", np.flatnonzero(((obs.class_code == 1) | (obs.class_code == 2)) & (obs.source == SUBOPTIMAL_SOURCE)).astype(np.int64)),
    ]:
        labels = obs.class_code[idx]
        overall_rows.append(
            {
                "split": name,
                "n_cells": int(idx.size),
                "gdT_gold": int((labels == 2).sum()),
                "abT_gold": int((labels == 1).sum()),
                "silver": int((labels == 3).sum()),
                "unlabeled": int((labels == 0).sum()),
                "NK_guard_negative": int(nk_guard_mask[idx].sum()),
                "TCRAB_guard_negative": int(tcrab_guard_mask[idx].sum()),
                "NK_TCRAB_overlap_in_split": int(nk_tcrab_overlap_mask[idx].sum()),
            }
        )
    split_overall = pd.DataFrame(overall_rows)
    split_overall.to_csv(TABLE_DIR / "split_overall.csv", index=False)
    return obs, SplitBundle(
        train_idx=train_idx,
        tune_idx=tune_idx,
        hard_negative_train_idx=hard_train,
        hard_negative_tune_idx=hard_tune,
        weak_trdc_prefilter_idx=weak_pool,
        validation_idx=validation_idx,
        validation_gdt2020_idx=validation_gdt2020_idx,
        validation_extra_trab_idx=validation_extra_trab_idx,
        split_overall=split_overall,
        split_by_source=split_by_source,
        annotation=annotation,
        nk_guard_mask=nk_guard_mask,
        tcrab_guard_mask=tcrab_guard_mask,
    )


def local_positions(eval_rows: np.ndarray, target_rows: np.ndarray) -> np.ndarray:
    lookup = pd.Series(np.arange(eval_rows.size, dtype=np.int64), index=eval_rows)
    return lookup.loc[target_rows].to_numpy(dtype=np.int64)


def filter_hard_negative_positions(
    eval_rows: np.ndarray,
    x_eval: np.ndarray,
    spec: FeatureSpec,
    split: SplitBundle,
    hard_positions: np.ndarray,
    split_name: str,
) -> tuple[np.ndarray, dict[str, Any]]:
    if hard_positions.size == 0:
        return hard_positions, {
            "split": split_name,
            "input_hard_prefilter": 0,
            "input_NK_prefilter": 0,
            "input_TCRAB_prefilter": 0,
            "input_NK_TCRAB_overlap": 0,
            "retained_total": 0,
            "retained_NK_TRDCpos_TRDVneg": 0,
            "retained_TCRAB_only": 0,
            "excluded_NK_not_TRDCpos_TRDVneg": 0,
            "excluded_NK_TCRAB_overlap": 0,
            "excluded_other": 0,
        }
    rows = eval_rows[hard_positions]
    is_nk = split.nk_guard_mask[rows]
    is_tcrab = split.tcrab_guard_mask[rows]
    nk_overlap = is_nk & is_tcrab
    trdc_expr = x_eval[hard_positions, spec.engineered_to_col["TRDC_log1p"]] > 0
    trdv_expr = x_eval[hard_positions, spec.engineered_to_col["any_TRDV"]] > 0.5
    nk_trdc_trdv_minus = is_nk & (~is_tcrab) & trdc_expr & (~trdv_expr)
    tcrab_only = is_tcrab & (~is_nk)
    keep = nk_trdc_trdv_minus | tcrab_only
    stats = {
        "split": split_name,
        "input_hard_prefilter": int(hard_positions.size),
        "input_NK_prefilter": int(is_nk.sum()),
        "input_TCRAB_prefilter": int(is_tcrab.sum()),
        "input_NK_TCRAB_overlap": int(nk_overlap.sum()),
        "retained_total": int(keep.sum()),
        "retained_NK_TRDCpos_TRDVneg": int(nk_trdc_trdv_minus.sum()),
        "retained_TCRAB_only": int(tcrab_only.sum()),
        "excluded_NK_not_TRDCpos_TRDVneg": int((is_nk & (~is_tcrab) & (~nk_trdc_trdv_minus)).sum()),
        "excluded_NK_TCRAB_overlap": int(nk_overlap.sum()),
        "excluded_other": int(((~keep) & (~is_nk) & (~is_tcrab)).sum()),
    }
    return hard_positions[keep].astype(np.int64), stats


def hard_negative_filter_by_source(
    eval_rows: np.ndarray,
    obs: Any,
    split: SplitBundle,
    train_positions: np.ndarray,
    tune_positions: np.ndarray,
) -> pd.DataFrame:
    rows_out: list[dict[str, Any]] = []
    for split_name, positions in [("train_hard_negative", train_positions), ("tune_hard_negative", tune_positions)]:
        rows = eval_rows[positions]
        if rows.size == 0:
            continue
        retained_nk = split.nk_guard_mask[rows] & (~split.tcrab_guard_mask[rows])
        retained_tcrab = split.tcrab_guard_mask[rows] & (~split.nk_guard_mask[rows])
        frame = pd.DataFrame(
            {
                "source_gse_id": obs.source[rows],
                "retained_NK_TRDCpos_TRDVneg": retained_nk.astype(np.int8),
                "retained_TCRAB_only": retained_tcrab.astype(np.int8),
                "retained_total": np.ones(rows.size, dtype=np.int8),
            }
        )
        grouped = frame.groupby("source_gse_id", sort=True).sum(numeric_only=True).reset_index()
        grouped.insert(0, "split", split_name)
        rows_out.extend(grouped.to_dict(orient="records"))
    columns = ["split", "source_gse_id", "retained_total", "retained_NK_TRDCpos_TRDVneg", "retained_TCRAB_only"]
    return pd.DataFrame(rows_out, columns=columns)


def candidate_thresholds(score: np.ndarray) -> np.ndarray:
    score = np.asarray(score, dtype=np.float32)
    finite = np.isfinite(score)
    if not finite.all():
        score = np.where(finite, score, 0.0).astype(np.float32)
    thresholds = np.unique(score)
    if thresholds.size > 1200:
        thresholds = np.unique(np.quantile(score, np.linspace(0.001, 0.999, 1200)).astype(np.float32))
    extras = np.asarray([float(np.nanmin(score)) - 1e-6, float(np.nanmax(score)) + 1e-6], dtype=np.float32)
    return np.unique(np.concatenate([thresholds, extras]))


def threshold_counts(y_true: np.ndarray, score: np.ndarray, threshold: float) -> dict[str, Any]:
    y = y_true.astype(np.int8, copy=False)
    pred = score >= threshold
    predicted = int(pred.sum())
    tp = int((pred & (y == 1)).sum())
    fp = predicted - tp
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    fn = n_pos - tp
    tn = n_neg - fp
    precision = safe_div(tp, tp + fp)
    recall = safe_div(tp, tp + fn)
    specificity = safe_div(tn, tn + fp)
    f1 = safe_div(2 * tp, 2 * tp + fp + fn)
    fp_fraction = safe_div(fp, predicted)
    return {
        "predicted_positive": predicted,
        "tp": tp,
        "fp": fp,
        "tn": tn,
        "fn": fn,
        "precision": 0.0 if not np.isfinite(precision) else float(precision),
        "recall": 0.0 if not np.isfinite(recall) else float(recall),
        "specificity": 0.0 if not np.isfinite(specificity) else float(specificity),
        "f1": 0.0 if not np.isfinite(f1) else float(f1),
        "fp_fraction_of_predictions": 1.0 if not np.isfinite(fp_fraction) else float(fp_fraction),
    }


def threshold_metric_dict(y_true: np.ndarray, score: np.ndarray, threshold: float) -> dict[str, Any]:
    y = y_true.astype(np.int8, copy=False)
    finite_score = np.asarray(score, dtype=np.float32)
    finite = np.isfinite(finite_score)
    if not finite.all():
        finite_score = np.where(finite, finite_score, 0.0).astype(np.float32)
    pred = (finite_score >= threshold).astype(np.int8)
    metrics = metric_dict(y, pred, finite_score)
    metrics["fp_fraction_of_predictions"] = safe_div(int(metrics.get("fp", 0)), int(metrics.get("predicted_positive", 0)))
    metrics["threshold"] = float(threshold)
    return metrics


def choose_threshold(y_true: np.ndarray, score: np.ndarray) -> tuple[float, dict[str, Any]]:
    y = y_true.astype(np.int8, copy=False)
    score = np.asarray(score, dtype=np.float32)
    finite = np.isfinite(score)
    if not finite.all():
        score = np.where(finite, score, 0.0).astype(np.float32)
    best_key: tuple[float, float, float, float, float] | None = None
    best_threshold: float | None = None
    for threshold in candidate_thresholds(score):
        metrics = threshold_counts(y, score, float(threshold))
        predicted = int(metrics.get("predicted_positive", 0))
        if predicted == 0:
            continue
        key = (
            float(metrics.get("f1", 0.0)),
            float(metrics.get("precision", 0.0)),
            float(metrics.get("recall", 0.0)),
            float(metrics.get("specificity", 0.0)),
            -float(predicted),
        )
        if best_key is None or key > best_key:
            best_key = key
            best_threshold = float(threshold)
    if best_threshold is None:
        best_threshold = float(np.nanmax(score) + 1e-6)
    metrics = threshold_metric_dict(y, score, best_threshold)
    metrics["threshold_policy"] = "max_f1"
    return best_threshold, metrics


def choose_target_threshold(
    y_true: np.ndarray,
    score: np.ndarray,
    *,
    target_recall: float,
    max_fp_fraction: float,
) -> tuple[float, dict[str, Any]]:
    y = y_true.astype(np.int8, copy=False)
    score = np.asarray(score, dtype=np.float32)
    finite = np.isfinite(score)
    if not finite.all():
        score = np.where(finite, score, 0.0).astype(np.float32)
    best_key: tuple[float, ...] | None = None
    best_threshold: float | None = None
    for threshold in candidate_thresholds(score):
        metrics = threshold_counts(y, score, float(threshold))
        predicted = int(metrics.get("predicted_positive", 0))
        if predicted == 0:
            continue
        recall = float(metrics.get("recall", 0.0))
        fp_fraction = float(metrics.get("fp_fraction_of_predictions", 1.0))
        precision = float(metrics.get("precision", 0.0))
        f1 = float(metrics.get("f1", 0.0))
        target_met = recall >= target_recall and fp_fraction <= max_fp_fraction
        if target_met:
            key = (1.0, recall, -fp_fraction, precision, f1, -float(predicted))
        else:
            recall_gap = max(0.0, target_recall - recall)
            fp_excess = max(0.0, fp_fraction - max_fp_fraction)
            key = (0.0, -recall_gap, -fp_excess, f1, precision, recall, -float(predicted))
        if best_key is None or key > best_key:
            best_key = key
            best_threshold = float(threshold)
    if best_threshold is None:
        best_threshold = float(np.nanmax(score) + 1e-6)
    best_metrics = threshold_metric_dict(y, score, best_threshold)
    best_metrics["threshold_policy"] = "target_recall_fp_fraction"
    best_metrics["target_recall"] = float(target_recall)
    best_metrics["target_max_fp_fraction"] = float(max_fp_fraction)
    best_metrics["target_met_on_tune"] = bool(
        float(best_metrics.get("recall", 0.0)) >= target_recall
        and float(best_metrics.get("fp_fraction_of_predictions", 1.0)) <= max_fp_fraction
    )
    return best_threshold, best_metrics

def sample_train_positions(
    pos_primary_train: np.ndarray,
    pos_hard_train: np.ndarray,
    y_eval: np.ndarray,
    x_eval: np.ndarray,
    spec: FeatureSpec,
    args: argparse.Namespace,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(args.seed + 3)
    positives = pos_primary_train[y_eval[pos_primary_train] == 1]
    primary_neg = pos_primary_train[y_eval[pos_primary_train] == 0]
    hard_neg = pos_hard_train[y_eval[pos_hard_train] == 0]
    if primary_neg.size > args.max_primary_negative_train:
        primary_neg = rng.choice(primary_neg, size=args.max_primary_negative_train, replace=False)
    if hard_neg.size > args.max_hard_negative_train:
        hard_neg = rng.choice(hard_neg, size=args.max_hard_negative_train, replace=False)
    out = np.unique(np.concatenate([positives, primary_neg, hard_neg])).astype(np.int64)
    rng.shuffle(out)
    weights = np.ones(out.size, dtype=np.float32)
    nk_score = x_eval[out, spec.engineered_to_col["NK_score"]]
    gd_score = x_eval[out, spec.engineered_to_col["gdT_TCR_score"]]
    cd3_score = x_eval[out, spec.engineered_to_col["CD3_score"]]
    is_pos = y_eval[out] == 1
    cytotoxic_gdt = is_pos & (nk_score >= np.nanpercentile(nk_score[is_pos], 70) if is_pos.any() else False) & (gd_score > 0)
    weights[cytotoxic_gdt] = 1.6
    hard_mask = np.isin(out, hard_neg)
    weights[hard_mask] = np.maximum(weights[hard_mask], 1.8)
    weak_trdc = (x_eval[out, spec.engineered_to_col["TRDC_only"]] > 0.5) & (x_eval[out, spec.engineered_to_col["NK_minus_CD3_score"]] > 0.25) & (cd3_score < 0.5)
    weights[weak_trdc & (~is_pos)] = np.maximum(weights[weak_trdc & (~is_pos)], 2.0)
    return out, weights


def fit_platt_calibrated_model(x_train: np.ndarray, y_train: np.ndarray, sample_weight: np.ndarray, seed: int) -> PlattCalibratedModel:
    rng = np.random.default_rng(seed)
    idx = np.arange(x_train.shape[0])
    rng.shuffle(idx)
    n_cal = max(2000, int(round(idx.size * 0.18)))
    n_cal = min(n_cal, max(1000, idx.size // 3))
    cal_idx = idx[:n_cal]
    fit_idx = idx[n_cal:]
    if fit_idx.size == 0:
        fit_idx = idx
        cal_idx = idx
    base = make_pipeline(
        StandardScaler(),
        SGDClassifier(
            loss="log_loss",
            penalty="elasticnet",
            l1_ratio=0.35,
            alpha=1e-4,
            class_weight="balanced",
            max_iter=15,
            tol=1e-3,
            n_jobs=32,
            random_state=seed,
        ),
    )
    base.fit(x_train[fit_idx], y_train[fit_idx], sgdclassifier__sample_weight=sample_weight[fit_idx])
    raw = base.predict_proba(x_train[cal_idx])[:, 1]
    logits = PlattCalibratedModel._logit(raw)
    calibrator = LogisticRegression(solver="lbfgs", max_iter=300, random_state=seed)
    calibrator.fit(logits, y_train[cal_idx], sample_weight=sample_weight[cal_idx])
    return PlattCalibratedModel(base, calibrator)


def platt_split_indices(n_rows: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    idx = np.arange(n_rows)
    rng.shuffle(idx)
    n_cal = max(2000, int(round(idx.size * 0.18)))
    n_cal = min(n_cal, max(1000, idx.size // 3))
    cal_idx = idx[:n_cal]
    fit_idx = idx[n_cal:]
    if fit_idx.size == 0:
        fit_idx = idx
        cal_idx = idx
    return fit_idx, cal_idx


def fit_platt_calibrated_base_model(
    base_model: Any,
    x_train: np.ndarray,
    y_train: np.ndarray,
    sample_weight: np.ndarray,
    seed: int,
) -> PlattCalibratedModel:
    fit_idx, cal_idx = platt_split_indices(x_train.shape[0], seed)
    base_model.fit(x_train[fit_idx], y_train[fit_idx], sample_weight=sample_weight[fit_idx])
    raw = base_model.predict_proba(x_train[cal_idx])[:, 1]
    calibrator = LogisticRegression(solver="lbfgs", max_iter=300, random_state=seed)
    calibrator.fit(PlattCalibratedModel._logit(raw), y_train[cal_idx], sample_weight=sample_weight[cal_idx])
    return PlattCalibratedModel(base_model, calibrator)


def fit_hist_gradient_model(x_train: np.ndarray, y_train: np.ndarray, sample_weight: np.ndarray, seed: int) -> PlattCalibratedModel:
    base = HistGradientBoostingClassifier(
        loss="log_loss",
        learning_rate=0.055,
        max_iter=140,
        max_leaf_nodes=23,
        min_samples_leaf=35,
        l2_regularization=0.15,
        early_stopping=True,
        validation_fraction=0.12,
        random_state=seed,
    )
    return fit_platt_calibrated_base_model(base, x_train, y_train, sample_weight, seed + 31)


def fit_lightgbm_model(x_train: np.ndarray, y_train: np.ndarray, sample_weight: np.ndarray, seed: int) -> PlattCalibratedModel | None:
    if not LIGHTGBM_AVAILABLE or LGBMClassifier is None:
        return None
    base = LGBMClassifier(
        objective="binary",
        n_estimators=420,
        learning_rate=0.035,
        num_leaves=31,
        min_child_samples=35,
        subsample=0.85,
        subsample_freq=1,
        colsample_bytree=0.85,
        reg_alpha=0.05,
        reg_lambda=3.0,
        random_state=seed,
        n_jobs=32,
        verbosity=-1,
    )
    return fit_platt_calibrated_base_model(base, x_train, y_train, sample_weight, seed + 41)


def fit_extra_trees_model(x_train: np.ndarray, y_train: np.ndarray, sample_weight: np.ndarray, seed: int) -> PlattCalibratedModel:
    base = ExtraTreesClassifier(
        n_estimators=80,
        max_depth=10,
        min_samples_leaf=20,
        max_features="sqrt",
        class_weight="balanced_subsample",
        n_jobs=16,
        random_state=seed,
    )
    return fit_platt_calibrated_base_model(base, x_train, y_train, sample_weight, seed + 51)


def append_candidate(
    candidates: list[CandidateResult],
    *,
    name: str,
    kind: str,
    model_object: Any,
    score: np.ndarray,
    y_tune: np.ndarray,
    notes: str,
    iteration_round: int,
    model_family: str,
    threshold_policy: str,
    args: argparse.Namespace,
    promotable: bool = True,
) -> None:
    if threshold_policy == "target_recall_fp_fraction":
        threshold, tune_metrics = choose_target_threshold(
            y_tune,
            score,
            target_recall=args.target_recall,
            max_fp_fraction=args.max_estimated_fp_fraction,
        )
    elif threshold_policy == "max_f1":
        threshold, tune_metrics = choose_threshold(y_tune, score)
    else:
        raise ValueError(f"Unknown threshold policy: {threshold_policy}")
    tune_metrics["iteration_round"] = int(iteration_round)
    tune_metrics["model_family"] = model_family
    tune_metrics["threshold_policy"] = threshold_policy
    candidates.append(
        CandidateResult(
            name=name,
            kind=kind,
            model_object=model_object,
            threshold=threshold,
            tune_metrics=tune_metrics,
            tune_score=score,
            tune_pred=(score >= threshold).astype(np.int8),
            notes=notes,
            promotable=promotable,
            iteration_round=iteration_round,
            model_family=model_family,
            threshold_policy=threshold_policy,
        )
    )


def append_fixed_threshold_candidate(
    candidates: list[CandidateResult],
    *,
    name: str,
    kind: str,
    model_object: Any,
    score: np.ndarray,
    y_tune: np.ndarray,
    threshold: float,
    notes: str,
    iteration_round: int,
    model_family: str,
    args: argparse.Namespace,
    promotable: bool = True,
) -> None:
    tune_metrics = threshold_metric_dict(y_tune, score, float(threshold))
    tune_metrics["iteration_round"] = int(iteration_round)
    tune_metrics["model_family"] = model_family
    tune_metrics["threshold_policy"] = "fixed_recall_rescue"
    tune_metrics["target_recall"] = float(args.target_recall)
    tune_metrics["target_max_fp_fraction"] = float(args.max_estimated_fp_fraction)
    tune_metrics["target_met_on_tune"] = bool(
        float(tune_metrics.get("recall", 0.0)) >= args.target_recall
        and float(tune_metrics.get("fp_fraction_of_predictions", 1.0)) <= args.max_estimated_fp_fraction
    )
    candidates.append(
        CandidateResult(
            name=name,
            kind=kind,
            model_object=model_object,
            threshold=float(threshold),
            tune_metrics=tune_metrics,
            tune_score=score,
            tune_pred=(score >= float(threshold)).astype(np.int8),
            notes=notes,
            promotable=promotable,
            iteration_round=iteration_round,
            model_family=model_family,
            threshold_policy="fixed_recall_rescue",
        )
    )


def train_candidates(
    x_eval: np.ndarray,
    y_eval: np.ndarray,
    positions: dict[str, np.ndarray],
    spec: FeatureSpec,
    args: argparse.Namespace,
    v2_payload: dict[str, Any],
) -> list[CandidateResult]:
    train_sample, train_weights = sample_train_positions(
        positions["primary_train"],
        positions["hard_train"],
        y_eval,
        x_eval,
        spec,
        args,
    )
    pd.DataFrame(
        [
            {"sample_class": "gdT_positive", "n_cells": int(y_eval[train_sample].sum())},
            {"sample_class": "primary_or_hard_negative", "n_cells": int((y_eval[train_sample] == 0).sum())},
            {"sample_class": "total", "n_cells": int(train_sample.size)},
        ]
    ).to_csv(TABLE_DIR / "training_sample_composition.csv", index=False)
    tune_pos = np.unique(np.concatenate([positions["primary_tune"], positions["hard_tune"]])).astype(np.int64)
    y_tune = y_eval[tune_pos]
    x_train = x_eval[train_sample]
    y_train = y_eval[train_sample]
    x_tune = x_eval[tune_pos]
    candidates: list[CandidateResult] = []
    target_tag = f"recall{int(round(args.target_recall * 100)):02d}_fp{int(round(args.max_estimated_fp_fraction * 100)):02d}"
    hist: PlattCalibratedModel | None = None
    hist_score: np.ndarray | None = None

    logging.info("Round 1: fitting elastic-net calibrated logistic v3 model")
    elastic = fit_platt_calibrated_model(x_train, y_train, train_weights, args.seed)
    elastic_score = elastic.predict_proba(x_tune)[:, 1].astype(np.float32)
    append_candidate(
        candidates,
        name="v3_elasticnet_conditional_features",
        kind="learned",
        model_object=elastic,
        score=elastic_score,
        y_tune=y_tune,
        notes="Round 1. Elastic-net logistic model on TCR genes, NK/CD3 controls, and engineered conditional features; Platt-calibrated on training holdback; threshold maximizes tune F1.",
        iteration_round=1,
        model_family="elasticnet_logistic",
        threshold_policy="max_f1",
        args=args,
    )
    append_candidate(
        candidates,
        name=f"v3_round02_elasticnet_target_{target_tag}",
        kind="learned_target",
        model_object=elastic,
        score=elastic_score,
        y_tune=y_tune,
        notes="Round 2. Same calibrated elastic-net score, but threshold selected for target recall and FP-fraction constraints.",
        iteration_round=2,
        model_family="elasticnet_logistic",
        threshold_policy="target_recall_fp_fraction",
        args=args,
    )

    if not args.skip_hist_gradient:
        logging.info("Round 3: fitting sklearn histogram-gradient v3 model")
        hist = fit_hist_gradient_model(x_train, y_train, train_weights, args.seed)
        hist_score = hist.predict_proba(x_tune)[:, 1].astype(np.float32)
        append_candidate(
            candidates,
            name=f"v3_round03_hist_gradient_target_{target_tag}",
            kind="learned_target",
            model_object=hist,
            score=hist_score,
            y_tune=y_tune,
            notes="Round 3. Non-XGBoost sklearn histogram-gradient model with engineered conditional features and target-aware threshold.",
            iteration_round=3,
            model_family="sklearn_hist_gradient",
            threshold_policy="target_recall_fp_fraction",
            args=args,
        )
        hist_fixed_thresholds = [0.80, 0.72, 0.64, 0.58, 0.55, 0.52, 0.50, 0.48, 0.45]
        for rescue_round, rescue_threshold in enumerate(hist_fixed_thresholds, start=6):
            append_fixed_threshold_candidate(
                candidates,
                name=f"v3_round{rescue_round:02d}_hist_gradient_fixed_{str(rescue_threshold).replace('.', 'p')}",
                kind="learned_fixed_threshold",
                model_object=hist,
                score=hist_score,
                y_tune=y_tune,
                threshold=rescue_threshold,
                notes=f"Round {rescue_round}. Histogram-gradient fixed threshold {rescue_threshold}; tests full-atlas recall rescue while preserving the v2 high-purity FP comparator and TRDC+TRDV- gate.",
                iteration_round=rescue_round,
                model_family="sklearn_hist_gradient",
                args=args,
            )

    if not args.skip_lightgbm:
        lightgbm = fit_lightgbm_model(x_train, y_train, train_weights, args.seed)
        if lightgbm is not None:
            logging.info("Round 4: fitted LightGBM v3 model")
            lightgbm_score = lightgbm.predict_proba(x_tune)[:, 1].astype(np.float32)
            append_candidate(
                candidates,
                name=f"v3_round04_lightgbm_target_{target_tag}",
                kind="learned_target",
                model_object=lightgbm,
                score=lightgbm_score,
                y_tune=y_tune,
                notes="Round 4. LightGBM model with engineered conditional features and target-aware threshold; not XGBoost.",
                iteration_round=4,
                model_family="lightgbm",
                threshold_policy="target_recall_fp_fraction",
                args=args,
            )
        else:
            logging.warning("LightGBM is not available; skipping round 4 LightGBM candidate")

    if args.use_extra_trees and not args.skip_extra_trees:
        round_id = 5 if any(c.iteration_round == 4 for c in candidates) else 4
        logging.info("Round %s: fitting ExtraTrees v3 model", round_id)
        extra = fit_extra_trees_model(x_train, y_train, train_weights, args.seed)
        extra_score = extra.predict_proba(x_tune)[:, 1].astype(np.float32)
        append_candidate(
            candidates,
            name=f"v3_round{round_id:02d}_extra_trees_target_{target_tag}",
            kind="learned_target",
            model_object=extra,
            score=extra_score,
            y_tune=y_tune,
            notes=f"Round {round_id}. ExtraTrees ensemble with engineered conditional features and target-aware threshold; not XGBoost.",
            iteration_round=round_id,
            model_family="extra_trees",
            threshold_policy="target_recall_fp_fraction",
            args=args,
        )

    round_id = max([c.iteration_round for c in candidates], default=3) + 1
    two_stage_target = ConditionalGateModel(
        elastic,
        base_cols=None,
        engineered_to_col=spec.engineered_to_col,
        mode="two_stage_tcell_gate",
    )
    two_stage_target_score = two_stage_target.predict_proba(x_tune)[:, 1].astype(np.float32)
    append_candidate(
        candidates,
        name=f"v3_round{round_id:02d}_two_stage_target_{target_tag}",
        kind="two_stage_target",
        model_object=two_stage_target,
        score=two_stage_target_score,
        y_tune=y_tune,
        notes=f"Round {round_id}. Elastic-net score with transcriptome T-cell gate and target-aware threshold; no additional tree model fit.",
        iteration_round=round_id,
        model_family="elasticnet_two_stage_gate",
        threshold_policy="target_recall_fp_fraction",
        args=args,
    )

    base = v2_payload["base_model"]
    v2_genes = [str(gene) for gene in base["gene_names"]]
    missing_v2 = [gene for gene in v2_genes if gene not in spec.gene_to_col]
    if missing_v2:
        raise KeyError(f"v2 model genes missing from v3 feature spec: {missing_v2[:10]}")
    v2_cols = [spec.gene_to_col[gene] for gene in v2_genes]
    v2_gate = ConditionalGateModel(
        base["model_object"],
        base_cols=v2_cols,
        engineered_to_col=spec.engineered_to_col,
        mode="trdc_only_nk_like_gate",
    )
    v2_gate_score = v2_gate.predict_proba(x_tune)[:, 1].astype(np.float32)
    round_id = max([c.iteration_round for c in candidates], default=4) + 1
    append_candidate(
        candidates,
        name=f"v3_round{round_id:02d}_v2_score_trdc_gate_target_{target_tag}",
        kind="postgate_target",
        model_object=v2_gate,
        score=v2_gate_score,
        y_tune=y_tune,
        notes=f"Round {round_id}. Frozen v2 score with conditional TRDC-only/NK-like post-gate and target-aware threshold; no TCR-seq metadata features.",
        iteration_round=round_id,
        model_family="v2_score_conditional_gate",
        threshold_policy="target_recall_fp_fraction",
        args=args,
        promotable=True,
    )
    v2_high_f1_threshold = float(v2_payload["operating_modes"]["high_f1"]["threshold"])
    for rescue_threshold in [0.95, 0.944, 0.94, 0.936, 0.93, v2_high_f1_threshold]:
        rescue_round = max([c.iteration_round for c in candidates], default=round_id) + 1
        threshold_label = f"{rescue_threshold:.3f}".replace('.', 'p')
        append_fixed_threshold_candidate(
            candidates,
            name=f"v3_round{rescue_round:02d}_v2_score_trdc_gate_fixed_{threshold_label}",
            kind="postgate_fixed_threshold",
            model_object=v2_gate,
            score=v2_gate_score,
            y_tune=y_tune,
            threshold=rescue_threshold,
            notes=f"Round {rescue_round}. Frozen v2 score with conditional TRDC-only/NK-like post-gate and fixed recall-rescue threshold {rescue_threshold:.6f}; no TCR-seq metadata features.",
            iteration_round=rescue_round,
            model_family="v2_score_conditional_gate",
            args=args,
            promotable=True,
        )
    append_candidate(
        candidates,
        name="v2_score_trdc_only_nk_like_postgate",
        kind="postgate",
        model_object=v2_gate,
        score=v2_gate_score,
        y_tune=y_tune,
        notes="Reference. Frozen v2 score with a conditional TRDC-only/NK-like post-gate; no TCR-seq metadata features; threshold maximizes tune F1.",
        iteration_round=round_id + 1,
        model_family="v2_score_conditional_gate",
        threshold_policy="max_f1",
        args=args,
        promotable=False,
    )

    v2_trdv_required = StrictTRDVGateModel(
        base["model_object"],
        base_cols=v2_cols,
        engineered_to_col=spec.engineered_to_col,
    )
    v2_trdv_required_score = v2_trdv_required.predict_proba(x_tune)[:, 1].astype(np.float32)
    strict_round = max([c.iteration_round for c in candidates], default=round_id) + 1
    append_candidate(
        candidates,
        name=f"v3_round{strict_round:02d}_v2_score_trdv_required_target_{target_tag}",
        kind="postgate_trdv_required_target",
        model_object=v2_trdv_required,
        score=v2_trdv_required_score,
        y_tune=y_tune,
        notes=f"Round {strict_round}. Frozen v2 score with an explicit TRDV-required gate for all TRDC+ cells; target-aware threshold.",
        iteration_round=strict_round,
        model_family="v2_score_trdv_required_gate",
        threshold_policy="target_recall_fp_fraction",
        args=args,
        promotable=True,
    )
    for strict_threshold in [0.99, 0.97, 0.95, 0.936, v2_high_f1_threshold]:
        strict_round = max([c.iteration_round for c in candidates], default=strict_round) + 1
        threshold_label = f"{strict_threshold:.3f}".replace('.', 'p')
        append_fixed_threshold_candidate(
            candidates,
            name=f"v3_round{strict_round:02d}_v2_score_trdv_required_fixed_{threshold_label}",
            kind="postgate_trdv_required_fixed_threshold",
            model_object=v2_trdv_required,
            score=v2_trdv_required_score,
            y_tune=y_tune,
            threshold=strict_threshold,
            notes=f"Round {strict_round}. Frozen v2 score with all TRDC+TRDV- cells blocked and fixed threshold {strict_threshold:.6f}.",
            iteration_round=strict_round,
            model_family="v2_score_trdv_required_gate",
            args=args,
            promotable=True,
        )

    if hist is not None and hist_score is not None:
        hist_trdv_required = StrictTRDVGateModel(
            hist,
            base_cols=None,
            engineered_to_col=spec.engineered_to_col,
        )
        hist_trdv_required_score = hist_trdv_required.predict_proba(x_tune)[:, 1].astype(np.float32)
        strict_round = max([c.iteration_round for c in candidates], default=strict_round) + 1
        append_candidate(
            candidates,
            name=f"v3_round{strict_round:02d}_hist_gradient_trdv_required_target_{target_tag}",
            kind="hist_gradient_trdv_required_target",
            model_object=hist_trdv_required,
            score=hist_trdv_required_score,
            y_tune=y_tune,
            notes=f"Round {strict_round}. Low-FP histogram-gradient score with all TRDC+TRDV- cells blocked; target-aware threshold.",
            iteration_round=strict_round,
            model_family="hist_gradient_trdv_required_gate",
            threshold_policy="target_recall_fp_fraction",
            args=args,
            promotable=True,
        )
        for strict_hist_threshold in [0.58, 0.55, 0.52, 0.50, 0.48, 0.45]:
            strict_round = max([c.iteration_round for c in candidates], default=strict_round) + 1
            threshold_label = f"{strict_hist_threshold:.3f}".replace('.', 'p')
            append_fixed_threshold_candidate(
                candidates,
                name=f"v3_round{strict_round:02d}_hist_gradient_trdv_required_fixed_{threshold_label}",
                kind="hist_gradient_trdv_required_fixed_threshold",
                model_object=hist_trdv_required,
                score=hist_trdv_required_score,
                y_tune=y_tune,
                threshold=strict_hist_threshold,
                notes=f"Round {strict_round}. Low-FP histogram-gradient score with all TRDC+TRDV- cells blocked and fixed threshold {strict_hist_threshold:.6f}.",
                iteration_round=strict_round,
                model_family="hist_gradient_trdv_required_gate",
                args=args,
                promotable=True,
            )
        for rescue_threshold in [0.99, 0.97, 0.95]:
            strict_round = max([c.iteration_round for c in candidates], default=strict_round) + 1
            rescue_label = f"{rescue_threshold:.3f}".replace('.', 'p')
            hybrid = TRDVRescuePostGateModel(
                hist,
                base["model_object"],
                base_cols=None,
                rescue_cols=v2_cols,
                engineered_to_col=spec.engineered_to_col,
                rescue_threshold=rescue_threshold,
                rescue_floor=0.97,
            )
            hybrid_score = hybrid.predict_proba(x_tune)[:, 1].astype(np.float32)
            append_candidate(
                candidates,
                name=f"v3_round{strict_round:02d}_hist_gradient_trdv_v2rescue_{rescue_label}_target_{target_tag}",
                kind="hist_gradient_trdv_v2_rescue_target",
                model_object=hybrid,
                score=hybrid_score,
                y_tune=y_tune,
                notes=f"Round {strict_round}. Low-FP histogram-gradient score, all TRDC+TRDV- cells blocked, and TRDV+ cells rescued when frozen v2 score is >= {rescue_threshold:.3f}.",
                iteration_round=strict_round,
                model_family="hist_gradient_trdv_rescue_gate",
                threshold_policy="target_recall_fp_fraction",
                args=args,
                promotable=True,
            )

    two_stage = ConditionalGateModel(
        elastic,
        base_cols=None,
        engineered_to_col=spec.engineered_to_col,
        mode="two_stage_tcell_gate",
    )
    two_stage_score = two_stage.predict_proba(x_tune)[:, 1].astype(np.float32)
    append_candidate(
        candidates,
        name="v3_two_stage_tcell_gate_then_tcr_classifier",
        kind="two_stage",
        model_object=two_stage,
        score=two_stage_score,
        y_tune=y_tune,
        notes="Reference. Elastic-net v3 score with a transcriptome T-cell gate that only blocks low-CD3/TRDC-only/NK-like cells; threshold maximizes tune F1.",
        iteration_round=round_id + 2,
        model_family="elasticnet_two_stage_gate",
        threshold_policy="max_f1",
        args=args,
    )

    if len({c.iteration_round for c in candidates if c.promotable}) < args.min_iteration_rounds:
        logging.warning(
            "Only %s promotable iteration rounds are available; requested at least %s.",
            len({c.iteration_round for c in candidates if c.promotable}),
            args.min_iteration_rounds,
        )
    return candidates

def candidate_metrics_frame(candidates: list[CandidateResult]) -> pd.DataFrame:
    rows = []
    for candidate in candidates:
        row = {
            "strategy": candidate.name,
            "iteration_round": candidate.iteration_round,
            "model_family": candidate.model_family,
            "threshold_policy": candidate.threshold_policy,
            "kind": candidate.kind,
            "threshold": candidate.threshold,
            "promotable": candidate.promotable,
            "notes": candidate.notes,
        }
        row.update(candidate.tune_metrics)
        rows.append(row)
    out = pd.DataFrame(rows).sort_values(["target_met_on_tune", "recall", "fp_fraction_of_predictions", "f1"], ascending=[False, False, True, False], na_position="last")
    out.to_csv(TABLE_DIR / "internal_tune_candidate_metrics.csv", index=False)
    return out


def select_internal_candidate(candidates: list[CandidateResult], args: argparse.Namespace) -> CandidateResult:
    promotable = [candidate for candidate in candidates if candidate.promotable]
    pool = promotable if promotable else candidates

    def key(candidate: CandidateResult) -> tuple[float, ...]:
        metrics = candidate.tune_metrics
        recall = float(metrics.get("recall", 0.0))
        fp_fraction = float(metrics.get("fp_fraction_of_predictions", 1.0))
        precision = float(metrics.get("precision", 0.0))
        f1 = float(metrics.get("f1", 0.0))
        target_met = recall >= args.target_recall and fp_fraction <= args.max_estimated_fp_fraction
        if target_met:
            return (1.0, recall, -fp_fraction, precision, f1, -float(metrics.get("predicted_positive", 0)))
        recall_gap = max(0.0, args.target_recall - recall)
        fp_excess = max(0.0, fp_fraction - args.max_estimated_fp_fraction)
        return (0.0, -recall_gap, -fp_excess, f1, precision, recall, -float(metrics.get("predicted_positive", 0)))

    return max(pool, key=key)


def external_annotation_for_v2(cell_type: np.ndarray) -> np.ndarray:
    normalized = pd.Series(cell_type, copy=False).astype(str).str.upper().str.strip()
    out = np.full(normalized.size, "OTHER", dtype=object)
    out[normalized.eq("VD2_GDT") | normalized.eq("GDT_CELL") | normalized.eq("GDT")] = "GDT_CELL"
    out[normalized.eq("CD8_T")] = "CD8_T"
    out[normalized.eq("CD4_T")] = "CD4_T"
    out[normalized.eq("TREG")] = "TREG"
    out[normalized.eq("NK") | normalized.eq("NK_CELL")] = "NK_CELL"
    return out


def v2_threshold_vector(v2_payload: dict[str, Any], mode: str, cell_type: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mode_payload = v2_payload["operating_modes"][mode]
    n = cell_type.size
    if mode == "high_f1":
        return np.full(n, float(mode_payload["threshold"]), dtype=np.float32), np.full(n, "all_cells", dtype=object)
    thresholds = mode_payload["annotation_thresholds"]
    ann = external_annotation_for_v2(cell_type)
    out = np.full(n, float(thresholds["other_threshold"]), dtype=np.float32)
    label = np.full(n, "other", dtype=object)
    mapping = {
        "GDT_CELL": ("gdt_threshold", "gdT_cell"),
        "CD8_T": ("cd8_threshold", "CD8_T"),
        "CD4_T": ("cd4_threshold", "CD4_T"),
        "TREG": ("treg_threshold", "Treg"),
        "NK_CELL": ("nk_threshold", "NK_cell"),
    }
    for ann_key, (threshold_key, display) in mapping.items():
        mask = ann == ann_key
        out[mask] = float(thresholds[threshold_key])
        label[mask] = display
    return out, label


def external_truth(handle: h5py.File, rows: np.ndarray) -> pd.DataFrame:
    real_gdt = read_bool_like(handle, "real_gdt")[rows]
    real_abt = read_bool_like(handle, "real_abt_tcr_strict")[rows]
    cell_type = read_string_like(handle, "cell_type", default="unknown")[rows]
    strict_label = read_string_like(handle, "tcr_strict_cell_label", default="unknown")[rows]
    neg_type = pd.Series(cell_type, copy=False).isin(["NK", "B_cell", "Inflammatory_myeloid"]).to_numpy(dtype=bool)
    primary_negative = (real_abt | neg_type) & (~real_gdt)
    primary_mask = real_gdt | primary_negative
    y = np.full(rows.size, -1, dtype=np.int8)
    y[primary_negative] = 0
    y[real_gdt] = 1
    cd4_like = read_bool_like(handle, "cd4_t_like")[rows] if "cd4_t_like" in handle["obs"] else pd.Series(cell_type).astype(str).eq("CD4_T").to_numpy(dtype=bool)
    out = pd.DataFrame(
        {
            "external_row": rows,
            "cell_id": obs_index_values(handle)[rows],
            "primary_eval": primary_mask,
            "y_true": y,
            "label_group": np.where(real_gdt, "real_gdT_positive", np.where(primary_negative, "primary_negative", "stress_ambiguous")),
            "real_gdt": real_gdt,
            "real_abt_tcr_strict": real_abt,
            "cell_type": cell_type,
            "tcr_strict_cell_label": strict_label,
            "tcr_paired_ab": read_bool_like(handle, "tcr_paired_ab")[rows],
            "tcr_paired_gd": read_bool_like(handle, "tcr_paired_gd")[rows],
            "tcr_has_ab": read_bool_like(handle, "tcr_has_ab")[rows],
            "tcr_has_gd": read_bool_like(handle, "tcr_has_gd")[rows],
            "has_TRD": read_bool_like(handle, "has_TRD")[rows],
            "has_TRG": read_bool_like(handle, "has_TRG")[rows],
            "CD4_T_like_warning": cd4_like,
            "tissue": read_string_like(handle, "tissue", default="unknown")[rows],
            "sample_id": read_string_like(handle, "sample_id", default="unknown")[rows],
            "library_id": read_string_like(handle, "library_id", default="unknown")[rows],
            "phase4_trd_score": read_float_obs(handle, "phase4_trd_score")[rows] if "phase4_trd_score" in handle["obs"] else np.nan,
            "phase4_trab_score": read_float_obs(handle, "phase4_trab_score")[rows] if "phase4_trab_score" in handle["obs"] else np.nan,
            "phase4_trd_minus_trab": read_float_obs(handle, "phase4_trd_minus_trab")[rows] if "phase4_trd_minus_trab" in handle["obs"] else np.nan,
            "n_counts": read_float_obs(handle, "n_counts")[rows] if "n_counts" in handle["obs"] else np.nan,
            "n_genes": read_float_obs(handle, "n_genes")[rows] if "n_genes" in handle["obs"] else np.nan,
            "pct_mito": read_float_obs(handle, "pct_mito")[rows] if "pct_mito" in handle["obs"] else np.nan,
        }
    )
    return out


def safe_metric_row(strategy: str, y_true: np.ndarray, score: np.ndarray, pred: np.ndarray, threshold: float | str) -> dict[str, Any]:
    row = {"strategy": strategy, "threshold": threshold}
    try:
        row.update(metric_dict(y_true.astype(np.int8), pred.astype(np.int8), score.astype(np.float32)))
    except Exception as exc:  # sklearn metrics can fail on single-class strata.
        tp = int(((pred == 1) & (y_true == 1)).sum())
        fp = int(((pred == 1) & (y_true == 0)).sum())
        tn = int(((pred == 0) & (y_true == 0)).sum())
        fn = int(((pred == 0) & (y_true == 1)).sum())
        row.update(
            {
                "n_cells": int(y_true.size),
                "n_positive": int((y_true == 1).sum()),
                "n_negative": int((y_true == 0).sum()),
                "predicted_positive": int(pred.sum()),
                "tp": tp,
                "fp": fp,
                "tn": tn,
                "fn": fn,
                "precision": safe_div(tp, tp + fp),
                "recall": safe_div(tp, tp + fn),
                "specificity": safe_div(tn, tn + fp),
                "f1": safe_div(2 * tp, 2 * tp + fp + fn),
                "mcc": np.nan,
                "roc_auc": np.nan,
                "pr_auc": np.nan,
                "metric_warning": str(exc),
            }
        )
    row["fp_fraction_of_predictions"] = safe_div(int(row.get("fp", 0)), int(row.get("predicted_positive", 0)))
    try:
        row["brier_score"] = float(brier_score_loss(y_true.astype(np.int8), np.clip(score, 0, 1)))
    except Exception:
        row["brier_score"] = np.nan
    return row


def reliability_table(df: pd.DataFrame, strategies: list[str]) -> pd.DataFrame:
    rows = []
    primary = df["primary_eval"].to_numpy(dtype=bool)
    y = df.loc[primary, "y_true"].to_numpy(dtype=np.int8)
    for strategy in strategies:
        score = df.loc[primary, f"{strategy}_score"].to_numpy(dtype=np.float32)
        bins = np.linspace(0.0, 1.0, 11)
        bin_id = np.digitize(score, bins, right=False) - 1
        bin_id = np.clip(bin_id, 0, 9)
        for i in range(10):
            mask = bin_id == i
            if not mask.any():
                continue
            rows.append(
                {
                    "strategy": strategy,
                    "probability_bin": f"{bins[i]:.1f}-{bins[i + 1]:.1f}",
                    "n_cells": int(mask.sum()),
                    "mean_score": float(score[mask].mean()),
                    "observed_positive_fraction": float(y[mask].mean()),
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "external_calibration_bins.csv", index=False)
    return out


def evaluate_external_predictions(df: pd.DataFrame, strategies: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    primary = df["primary_eval"].to_numpy(dtype=bool)
    y = df.loc[primary, "y_true"].to_numpy(dtype=np.int8)
    metric_rows = []
    for strategy in strategies:
        score = df.loc[primary, f"{strategy}_score"].to_numpy(dtype=np.float32)
        pred = df.loc[primary, f"{strategy}_pred"].to_numpy(dtype=bool)
        threshold = df[f"{strategy}_threshold"].iloc[0] if f"{strategy}_threshold" in df else np.nan
        metric_rows.append(safe_metric_row(strategy, y, score, pred, threshold))
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(TABLE_DIR / "external_primary_metrics.csv", index=False)

    fp_rows = []
    recall_rows = []
    group_rows = []
    for strategy in strategies:
        pred = df[f"{strategy}_pred"].to_numpy(dtype=bool)
        score = df[f"{strategy}_score"].to_numpy(dtype=np.float32)
        fp = pred & primary & (df["y_true"].to_numpy(dtype=np.int8) == 0)
        fp_groups = {
            "NK_cells": pd.Series(df["cell_type"]).astype(str).eq("NK").to_numpy(dtype=bool),
            "paired_TCRAB_cells": df["tcr_paired_ab"].to_numpy(dtype=bool),
            "TRDC_plus_TRDV_minus_cells": pd.Series(df["tcr_gene_quadrant"]).astype(str).eq("TRDC+TRDV-").to_numpy(dtype=bool),
            "CD4_Treg_warning_cells": df["CD4_T_like_warning"].to_numpy(dtype=bool) | pd.Series(df["cell_type"]).astype(str).eq("CD4_T").to_numpy(dtype=bool),
            "B_myeloid_cells": pd.Series(df["cell_type"]).astype(str).isin(["B_cell", "Inflammatory_myeloid"]).to_numpy(dtype=bool),
        }
        for group_name, group_mask in fp_groups.items():
            denom = primary & (df["y_true"].to_numpy(dtype=np.int8) == 0) & group_mask
            fp_rows.append(
                {
                    "strategy": strategy,
                    "false_positive_group": group_name,
                    "group_cells": int(denom.sum()),
                    "false_positive_cells": int((fp & group_mask).sum()),
                    "false_positive_rate": safe_div(int((fp & group_mask).sum()), int(denom.sum())),
                }
            )
        real_gdt = df["real_gdt"].to_numpy(dtype=bool)
        recall_groups = {
            "all_real_gdt": real_gdt,
            "TRDV_positive_gdT": real_gdt & (df["any_TRDV"].to_numpy(dtype=bool)),
            "TRDC_plus_TRDV_minus_gdT": real_gdt & pd.Series(df["tcr_gene_quadrant"]).astype(str).eq("TRDC+TRDV-").to_numpy(dtype=bool),
            "cytotoxic_NK_marker_high_gdT": real_gdt & (df["NK_score"].to_numpy(dtype=float) >= np.nanpercentile(df.loc[real_gdt, "NK_score"], 70) if real_gdt.any() else False),
            "Vd2_gdT_cells": real_gdt & pd.Series(df["cell_type"]).astype(str).eq("Vd2_gdT").to_numpy(dtype=bool),
        }
        for group_name, group_mask in recall_groups.items():
            recall_rows.append(
                {
                    "strategy": strategy,
                    "recall_group": group_name,
                    "group_cells": int(group_mask.sum()),
                    "true_positive_cells": int((pred & group_mask).sum()),
                    "recall": safe_div(int((pred & group_mask).sum()), int(group_mask.sum())),
                }
            )
        for group_col in ["cell_type", "tcr_strict_cell_label", "tcr_gene_quadrant", "tissue", "sample_id"]:
            for value, group in df.groupby(group_col, dropna=False, sort=True):
                mask = group.index.to_numpy()
                group_rows.append(
                    {
                        "strategy": strategy,
                        "group_column": group_col,
                        "group_value": str(value),
                        "n_cells": int(group.shape[0]),
                        "primary_eval_cells": int(group["primary_eval"].sum()),
                        "real_gdt_cells": int(group["real_gdt"].sum()),
                        "predicted_cells": int(pred[mask].sum()),
                        "mean_score": float(np.nanmean(score[mask])) if mask.size else np.nan,
                    }
                )
    fp_df = pd.DataFrame(fp_rows)
    recall_df = pd.DataFrame(recall_rows)
    group_df = pd.DataFrame(group_rows)
    fp_df.to_csv(TABLE_DIR / "external_false_positive_groups.csv", index=False)
    recall_df.to_csv(TABLE_DIR / "external_recall_groups.csv", index=False)
    group_df.to_csv(TABLE_DIR / "external_prediction_by_group.csv", index=False)
    reliability_table(df, strategies)
    return metrics, fp_df, recall_df, group_df


def plot_external_summary(metrics: pd.DataFrame, fp: pd.DataFrame, recall: pd.DataFrame) -> list[Path]:
    paths: list[Path] = []
    fig, ax = plt.subplots(figsize=(9.4, 4.8), constrained_layout=True)
    plot = metrics.copy()
    x = np.arange(plot.shape[0])
    ax.bar(x - 0.24, plot["precision"], width=0.24, label="precision", color="#0f766e")
    ax.bar(x, plot["recall"], width=0.24, label="recall", color="#1d4ed8")
    ax.bar(x + 0.24, plot["f1"], width=0.24, label="F1", color="#7c3aed")
    ax.set_ylim(0, 1.02)
    ax.set_xticks(x)
    ax.set_xticklabels(plot["strategy"], rotation=25, ha="right", fontsize=8)
    ax.set_title("External primary metrics")
    ax.legend(frameon=False)
    path = FIGURE_DIR / "external_primary_metrics.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(10.2, 4.8), constrained_layout=True)
    pivot = fp.pivot_table(index="strategy", columns="false_positive_group", values="false_positive_cells", aggfunc="sum").fillna(0)
    pivot.plot(kind="bar", stacked=False, ax=ax, width=0.82)
    ax.set_ylabel("false-positive cells")
    ax.set_title("External false positives by stress group")
    ax.tick_params(axis="x", labelrotation=25)
    ax.legend(frameon=False, fontsize=7)
    path = FIGURE_DIR / "external_false_positive_groups.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(10.2, 4.8), constrained_layout=True)
    pivot = recall.pivot_table(index="strategy", columns="recall_group", values="recall", aggfunc="first").fillna(0)
    pivot.plot(kind="bar", ax=ax, width=0.82)
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("recall")
    ax.set_title("External gdT recall by subgroup")
    ax.tick_params(axis="x", labelrotation=25)
    ax.legend(frameon=False, fontsize=7)
    path = FIGURE_DIR / "external_recall_groups.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    paths.append(path)
    return paths


def write_internal_tune_predictions(
    x_eval: np.ndarray,
    y_eval: np.ndarray,
    eval_rows: np.ndarray,
    positions: dict[str, np.ndarray],
    candidates: list[CandidateResult],
    spec: FeatureSpec,
    obs: Any,
    split: SplitBundle,
) -> pd.DataFrame:
    tune_pos = np.unique(np.concatenate([positions["primary_tune"], positions["hard_tune"]])).astype(np.int64)
    rows = eval_rows[tune_pos]
    df = pd.DataFrame(
        {
            "obs_index": rows,
            "source_gse_id": obs.source[rows],
            "tissue": obs.tissue[rows],
            "annotation": split.annotation[rows],
            "y_true": y_eval[tune_pos].astype(np.int8),
            "class_code": obs.class_code[rows].astype(int),
            "is_NK_guard": split.nk_guard_mask[rows],
            "is_TCRAB_guard": split.tcrab_guard_mask[rows],
            "tcr_gene_quadrant": trdc_trdv_quadrant(x_eval[tune_pos], spec),
            "CD3_score": x_eval[tune_pos, spec.engineered_to_col["CD3_score"]],
            "NK_score": x_eval[tune_pos, spec.engineered_to_col["NK_score"]],
            "gdT_TCR_score": x_eval[tune_pos, spec.engineered_to_col["gdT_TCR_score"]],
        }
    )
    for candidate in candidates:
        score = candidate.model_object.predict_proba(x_eval[tune_pos])[:, 1].astype(np.float32)
        df[f"{candidate.name}_score"] = score
        df[f"{candidate.name}_pred"] = score >= candidate.threshold
    path = TABLE_DIR / "internal_tune_predictions_wide.csv.gz"
    df.to_csv(path, index=False, compression="gzip")
    source_rows = []
    for candidate in candidates:
        for source, group in df.groupby("source_gse_id", sort=True):
            y = group["y_true"].to_numpy(dtype=np.int8)
            pred = group[f"{candidate.name}_pred"].to_numpy(dtype=bool)
            score = group[f"{candidate.name}_score"].to_numpy(dtype=np.float32)
            if np.unique(y).size < 2:
                continue
            row = safe_metric_row(candidate.name, y, score, pred, candidate.threshold)
            row["source_gse_id"] = source
            source_rows.append(row)
    source_df = pd.DataFrame(source_rows)
    source_df.to_csv(TABLE_DIR / "internal_tune_leave_one_source_diagnostics.csv", index=False)
    return df



def validation_component_labels(rows: np.ndarray, split: SplitBundle) -> np.ndarray:
    out = np.full(rows.size, "unknown_validation", dtype=object)
    if split.validation_gdt2020_idx.size:
        out[np.isin(rows, split.validation_gdt2020_idx)] = f"gdT_{GDT2020_SOURCE}_{GDT2020_HOLDOUT_TISSUE.replace(' ', '_')}"
    if split.validation_extra_trab_idx.size:
        out[np.isin(rows, split.validation_extra_trab_idx)] = f"abT_{EXTRA_TRAB_HOLDOUT_SOURCE}_paired_TCRAB_no_gdTCR"
    return out


def write_atlas_validation_predictions(
    x_eval: np.ndarray,
    y_eval: np.ndarray,
    eval_rows: np.ndarray,
    positions: dict[str, np.ndarray],
    candidates: list[CandidateResult],
    spec: FeatureSpec,
    obs: Any,
    split: SplitBundle,
    v2_payload: dict[str, Any],
    phase4_trd_minus_trab_eval: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    val_pos = positions.get("validation", np.asarray([], dtype=np.int64))
    rows = eval_rows[val_pos]
    if rows.size == 0:
        empty = pd.DataFrame()
        empty.to_csv(TABLE_DIR / "atlas_validation_predictions_wide.csv.gz", index=False, compression="gzip")
        empty.to_csv(TABLE_DIR / "atlas_validation_metrics.csv", index=False)
        empty.to_csv(TABLE_DIR / "atlas_validation_metrics_by_group.csv", index=False)
        return empty, empty, empty
    df = pd.DataFrame(
        {
            "obs_index": rows,
            "source_gse_id": obs.source[rows],
            "tissue": obs.tissue[rows],
            "annotation": split.annotation[rows],
            "validation_component": validation_component_labels(rows, split),
            "y_true": y_eval[val_pos].astype(np.int8),
            "class_code": obs.class_code[rows].astype(int),
            "is_NK_guard": split.nk_guard_mask[rows],
            "is_TCRAB_guard": split.tcrab_guard_mask[rows],
            "tcr_gene_quadrant": trdc_trdv_quadrant(x_eval[val_pos], spec),
            "CD3_score": x_eval[val_pos, spec.engineered_to_col["CD3_score"]],
            "NK_score": x_eval[val_pos, spec.engineered_to_col["NK_score"]],
            "gdT_TCR_score": x_eval[val_pos, spec.engineered_to_col["gdT_TCR_score"]],
        }
    )
    strategies: list[str] = []
    base = v2_payload["base_model"]
    v2_cols = [spec.gene_to_col[str(gene)] for gene in base["gene_names"]]
    v2_score = base["model_object"].predict_proba(x_eval[val_pos][:, v2_cols])[:, 1].astype(np.float32)
    for mode in ["high_f1", "high_purity"]:
        mode_key = f"v2_{mode}"
        threshold, threshold_label = v2_threshold_vector(v2_payload, mode, df["annotation"].to_numpy(dtype=object))
        df[f"{mode_key}_score"] = v2_score
        df[f"{mode_key}_threshold"] = threshold
        df[f"{mode_key}_threshold_label"] = threshold_label
        df[f"{mode_key}_pred"] = v2_score >= threshold
        strategies.append(mode_key)
    original_threshold = load_original_trd_trab_threshold()
    original_score = phase4_trd_minus_trab_eval[val_pos].astype(np.float32)
    df["original_TRD_minus_TRAB_score"] = original_score
    df["original_TRD_minus_TRAB_threshold"] = original_threshold
    df["original_TRD_minus_TRAB_pred"] = original_score >= original_threshold
    strategies.append("original_TRD_minus_TRAB")
    validation_score_cache: dict[int, np.ndarray] = {}
    for candidate in candidates:
        cache_key = id(candidate.model_object)
        if cache_key not in validation_score_cache:
            validation_score_cache[cache_key] = candidate.model_object.predict_proba(x_eval[val_pos])[:, 1].astype(np.float32)
        score = validation_score_cache[cache_key]
        df[f"{candidate.name}_score"] = score
        df[f"{candidate.name}_threshold"] = candidate.threshold
        df[f"{candidate.name}_pred"] = score >= candidate.threshold
        strategies.append(candidate.name)
    df.to_csv(TABLE_DIR / "atlas_validation_predictions_wide.csv.gz", index=False, compression="gzip")

    y = df["y_true"].to_numpy(dtype=np.int8)
    metric_rows = []
    group_rows = []
    for strategy in strategies:
        score = df[f"{strategy}_score"].to_numpy(dtype=np.float32)
        pred = df[f"{strategy}_pred"].to_numpy(dtype=bool)
        threshold = df[f"{strategy}_threshold"].iloc[0] if f"{strategy}_threshold" in df else np.nan
        metric_rows.append(safe_metric_row(strategy, y, score, pred, threshold))
        for group_col in ["validation_component", "source_gse_id", "tissue", "tcr_gene_quadrant"]:
            for value, group in df.groupby(group_col, sort=True):
                idx = group.index.to_numpy(dtype=np.int64)
                row = safe_metric_row(strategy, y[idx], score[idx], pred[idx], threshold)
                row["group_column"] = group_col
                row["group_value"] = value
                group_rows.append(row)
    metrics = pd.DataFrame(metric_rows)
    groups = pd.DataFrame(group_rows)
    metrics.to_csv(TABLE_DIR / "atlas_validation_metrics.csv", index=False)
    groups.to_csv(TABLE_DIR / "atlas_validation_metrics_by_group.csv", index=False)
    return df, metrics, groups


def quadrant_from_flags(trdc: np.ndarray, trdv: np.ndarray) -> np.ndarray:
    out = np.full(trdc.size, "TRDC-TRDV-", dtype=object)
    out[trdc & (~trdv)] = "TRDC+TRDV-"
    out[(~trdc) & trdv] = "TRDC-TRDV+"
    out[trdc & trdv] = "TRDC+TRDV+"
    return out


def full_group_summary(
    *,
    strategy: str,
    threshold: float | str,
    pred: np.ndarray,
    score: np.ndarray,
    y: np.ndarray,
    primary_mask: np.ndarray,
    source: np.ndarray,
    tissue: np.ndarray,
    annotation: np.ndarray,
    quadrant: np.ndarray,
    is_nk: np.ndarray,
    is_tcrab: np.ndarray,
    is_cd4_treg: np.ndarray,
    is_b_myeloid: np.ndarray,
    source_has_tcrseq: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pred_bool = pred.astype(bool, copy=False)
    primary_pred = pred_bool[primary_mask]
    primary_score = score[primary_mask].astype(np.float32, copy=False)
    primary_y = y[primary_mask]
    metric = safe_metric_row(strategy, primary_y, primary_score, primary_pred, threshold)
    gold_gdt = primary_mask & (y == 1)
    predicted_count = int(pred_bool.sum())
    observed_paired_abt_fp = int((pred_bool & is_tcrab).sum())
    tcrseq_pred_non_gold = pred_bool & source_has_tcrseq & (~gold_gdt)
    non_tcrseq_pred_non_gold = pred_bool & (~source_has_tcrseq) & (~gold_gdt)
    tcrseq_basis_count = int(tcrseq_pred_non_gold.sum())
    non_tcrseq_basis_count = int(non_tcrseq_pred_non_gold.sum())
    observed_abt_fp_rate = safe_div(observed_paired_abt_fp, tcrseq_basis_count)
    estimated_hidden_abt_fp = float(observed_abt_fp_rate * non_tcrseq_basis_count) if tcrseq_basis_count else np.nan
    estimated_total_abt_fp = float(observed_paired_abt_fp + estimated_hidden_abt_fp) if np.isfinite(estimated_hidden_abt_fp) else np.nan
    estimated_fp_fraction = safe_div(estimated_total_abt_fp, predicted_count) if np.isfinite(estimated_total_abt_fp) else np.nan
    overall = {
        "strategy": strategy,
        "threshold": threshold,
        "total_cells": int(pred_bool.size),
        "predicted_putative_gdT": int(pred_bool.sum()),
        "predicted_fraction": safe_div(predicted_count, int(pred_bool.size)),
        "primary_gold_cells": int(primary_mask.sum()),
        "tcrseq_source_cells": int(source_has_tcrseq.sum()),
        "non_tcrseq_source_cells": int((~source_has_tcrseq).sum()),
        "predicted_in_tcrseq_sources": int((pred_bool & source_has_tcrseq).sum()),
        "predicted_without_tcrseq_sources": int((pred_bool & (~source_has_tcrseq)).sum()),
        "predicted_tcrseq_non_gold_basis": tcrseq_basis_count,
        "predicted_non_tcrseq_non_gold_basis": non_tcrseq_basis_count,
        "observed_paired_TCRAB_abT_fp": observed_paired_abt_fp,
        "observed_abT_fp_rate_in_tcrseq_non_gold_predictions": float(observed_abt_fp_rate),
        "estimated_hidden_abT_fp_without_tcrseq": estimated_hidden_abt_fp,
        "estimated_total_abT_fp": estimated_total_abt_fp,
        "estimated_fp_fraction_of_predictions": estimated_fp_fraction,
        "full_primary_tp": int(metric.get("tp", 0)),
        "full_primary_fp": int(metric.get("fp", 0)),
        "full_primary_tn": int(metric.get("tn", 0)),
        "full_primary_fn": int(metric.get("fn", 0)),
        "full_primary_precision": float(metric.get("precision", np.nan)),
        "full_primary_recall": float(metric.get("recall", np.nan)),
        "full_primary_specificity": float(metric.get("specificity", np.nan)),
        "full_primary_f1": float(metric.get("f1", np.nan)),
        "full_primary_mcc": float(metric.get("mcc", np.nan)),
        "paired_TCRAB_cells": int(is_tcrab.sum()),
        "predicted_paired_TCRAB": int((pred_bool & is_tcrab).sum()),
        "NK_cells": int(is_nk.sum()),
        "predicted_NK": int((pred_bool & is_nk).sum()),
        "TRDC_plus_TRDV_minus_cells": int((quadrant == "TRDC+TRDV-").sum()),
        "predicted_TRDC_plus_TRDV_minus": int((pred_bool & (quadrant == "TRDC+TRDV-")).sum()),
        "CD4_Treg_warning_cells": int(is_cd4_treg.sum()),
        "predicted_CD4_Treg_warning": int((pred_bool & is_cd4_treg).sum()),
        "B_myeloid_cells": int(is_b_myeloid.sum()),
        "predicted_B_myeloid": int((pred_bool & is_b_myeloid).sum()),
    }
    full_overall = pd.DataFrame([overall])

    base = pd.DataFrame(
        {
            "source_gse_id": source,
            "tissue": tissue,
            "annotation": annotation,
            "predicted_putative_gdT": pred_bool.astype(np.int8),
            "primary_gold": primary_mask.astype(np.int8),
            "gdT_gold": ((y == 1) & primary_mask).astype(np.int8),
            "abT_gold": ((y == 0) & primary_mask).astype(np.int8),
            "paired_TCRAB": is_tcrab.astype(np.int8),
            "predicted_paired_TCRAB": (pred_bool & is_tcrab).astype(np.int8),
            "NK_cell": is_nk.astype(np.int8),
            "predicted_NK": (pred_bool & is_nk).astype(np.int8),
            "TRDC_plus_TRDV_minus": (quadrant == "TRDC+TRDV-").astype(np.int8),
            "predicted_TRDC_plus_TRDV_minus": (pred_bool & (quadrant == "TRDC+TRDV-")).astype(np.int8),
            "CD4_Treg_warning": is_cd4_treg.astype(np.int8),
            "predicted_CD4_Treg_warning": (pred_bool & is_cd4_treg).astype(np.int8),
            "B_myeloid": is_b_myeloid.astype(np.int8),
            "predicted_B_myeloid": (pred_bool & is_b_myeloid).astype(np.int8),
        }
    )

    def summarize(group_col: str) -> pd.DataFrame:
        out = base.groupby(group_col, dropna=False, as_index=False).agg(
            total_cells=("predicted_putative_gdT", "size"),
            predicted_putative_gdT=("predicted_putative_gdT", "sum"),
            primary_gold_cells=("primary_gold", "sum"),
            gdT_gold=("gdT_gold", "sum"),
            abT_gold=("abT_gold", "sum"),
            paired_TCRAB_cells=("paired_TCRAB", "sum"),
            predicted_paired_TCRAB=("predicted_paired_TCRAB", "sum"),
            NK_cells=("NK_cell", "sum"),
            predicted_NK=("predicted_NK", "sum"),
            TRDC_plus_TRDV_minus_cells=("TRDC_plus_TRDV_minus", "sum"),
            predicted_TRDC_plus_TRDV_minus=("predicted_TRDC_plus_TRDV_minus", "sum"),
            CD4_Treg_warning_cells=("CD4_Treg_warning", "sum"),
            predicted_CD4_Treg_warning=("predicted_CD4_Treg_warning", "sum"),
            B_myeloid_cells=("B_myeloid", "sum"),
            predicted_B_myeloid=("predicted_B_myeloid", "sum"),
        )
        out.insert(0, "strategy", strategy)
        out.insert(1, "threshold", threshold)
        out["predicted_fraction"] = out["predicted_putative_gdT"] / out["total_cells"].replace(0, np.nan)
        return out.replace([np.inf, -np.inf], np.nan)

    return full_overall, summarize("source_gse_id"), summarize("tissue"), summarize("annotation")


def select_candidate_from_full_atlas(
    candidates: list[CandidateResult],
    full_overall: pd.DataFrame,
    args: argparse.Namespace,
) -> CandidateResult:
    candidate_by_name = {candidate.name: candidate for candidate in candidates if candidate.promotable}
    rows = full_overall.loc[full_overall["strategy"].isin(candidate_by_name)].copy()
    if rows.empty:
        return select_internal_candidate(candidates, args)

    comparator_strategy = str(getattr(args, "fp_comparator_strategy", "v2_high_purity"))
    comparator = full_overall.loc[full_overall["strategy"].astype(str) == comparator_strategy]
    comparator_fp = np.nan
    if not comparator.empty:
        comparator_fp = float(comparator["estimated_fp_fraction_of_predictions"].iloc[0])
    effective_max_fp = float(args.max_estimated_fp_fraction)
    if np.isfinite(comparator_fp):
        effective_max_fp = min(effective_max_fp, comparator_fp)

    rows["target_recall"] = float(args.target_recall)
    rows["target_recall_margin"] = float(args.target_recall_margin)
    rows["target_max_estimated_fp_fraction"] = float(args.max_estimated_fp_fraction)
    rows["fp_comparator_strategy"] = comparator_strategy
    rows["fp_comparator_estimated_fp_fraction"] = comparator_fp
    rows["effective_max_estimated_fp_fraction"] = effective_max_fp
    rows["target_max_TRDC_plus_TRDV_minus_fraction"] = float(args.max_trdc_trdv_minus_fraction)
    rows["predicted_TRDC_plus_TRDV_minus_fraction"] = (
        rows["predicted_TRDC_plus_TRDV_minus"].astype(float) / rows["predicted_putative_gdT"].replace(0, np.nan).astype(float)
    ).fillna(1.0)
    rows["target_recall_met_full_atlas"] = rows["full_primary_recall"].astype(float) > args.target_recall
    rows["target_estimated_fp_fraction_met_full_atlas"] = rows["estimated_fp_fraction_of_predictions"].astype(float) <= effective_max_fp
    rows["target_TRDC_plus_TRDV_minus_fraction_met_full_atlas"] = rows["predicted_TRDC_plus_TRDV_minus_fraction"].astype(float) < args.max_trdc_trdv_minus_fraction
    rows["target_met_full_atlas"] = (
        rows["target_recall_met_full_atlas"]
        & rows["target_estimated_fp_fraction_met_full_atlas"]
        & rows["target_TRDC_plus_TRDV_minus_fraction_met_full_atlas"]
    )
    rows["iteration_round"] = rows["strategy"].map(lambda name: candidate_by_name[str(name)].iteration_round)
    rows["model_family"] = rows["strategy"].map(lambda name: candidate_by_name[str(name)].model_family)
    rows["threshold_policy"] = rows["strategy"].map(lambda name: candidate_by_name[str(name)].threshold_policy)
    rows.to_csv(TABLE_DIR / "full_atlas_target_selection.csv", index=False)

    def key(row: pd.Series) -> tuple[float, ...]:
        recall = float(row.get("full_primary_recall", 0.0))
        estimated_fp = float(row.get("estimated_fp_fraction_of_predictions", np.nan))
        if not np.isfinite(estimated_fp):
            estimated_fp = 1.0
        trdc_fraction = float(row.get("predicted_TRDC_plus_TRDV_minus_fraction", np.nan))
        if not np.isfinite(trdc_fraction):
            trdc_fraction = 1.0
        f1 = float(row.get("full_primary_f1", 0.0))
        precision = float(row.get("full_primary_precision", 0.0))
        predicted = float(row.get("predicted_putative_gdT", 0.0))
        target_met = (
            recall > args.target_recall
            and estimated_fp <= effective_max_fp
            and trdc_fraction < args.max_trdc_trdv_minus_fraction
        )
        recall_margin_met = recall >= args.target_recall + args.target_recall_margin
        if target_met:
            return (1.0, float(recall_margin_met), -estimated_fp, -trdc_fraction, recall, f1, precision, -predicted)
        recall_gap = max(0.0, args.target_recall - recall)
        fp_excess = max(0.0, estimated_fp - effective_max_fp)
        trdc_excess = max(0.0, trdc_fraction - args.max_trdc_trdv_minus_fraction)
        return (0.0, -recall_gap, -fp_excess, -trdc_excess, f1, precision, recall, -predicted)

    selected_row = max([row for _, row in rows.iterrows()], key=key)
    return candidate_by_name[str(selected_row["strategy"])]


def apply_full_atlas(
    handle: h5py.File,
    spec: FeatureSpec,
    candidates: list[CandidateResult],
    selected: CandidateResult,
    v2_payload: dict[str, Any],
    obs: Any,
    split: SplitBundle,
    args: argparse.Namespace,
) -> dict[str, pd.DataFrame]:
    n_obs = obs.source.size
    n_eval = n_obs if args.max_full_atlas_cells is None else min(n_obs, int(args.max_full_atlas_cells))
    rows_all = np.arange(n_eval, dtype=np.int64)
    source = obs.source[rows_all]
    tissue = obs.tissue[rows_all]
    annotation = split.annotation[rows_all]
    annotation_upper = pd.Series(annotation, copy=False).astype(str).str.upper().to_numpy(dtype=object)
    y = np.zeros(n_eval, dtype=np.int8)
    y[obs.class_code[rows_all] == 2] = 1
    primary_mask = np.isin(obs.class_code[rows_all], [1, 2])
    is_nk = np.char.find(annotation_upper.astype(str), "NK") >= 0
    is_tcrab = (obs.has_TRA_TRB_paired[rows_all] | obs.has_any_ab_tcr[rows_all]) & (obs.class_code[rows_all] != 2)
    is_cd4_treg = np.char.find(annotation_upper.astype(str), "CD4") >= 0
    is_cd4_treg |= np.char.find(annotation_upper.astype(str), "TREG") >= 0
    is_b_myeloid = np.char.find(annotation_upper.astype(str), "B_CELL") >= 0
    is_b_myeloid |= np.char.find(annotation_upper.astype(str), "MYELOID") >= 0
    tcr_evidence = obs.has_TRA_TRB_paired[rows_all] | obs.has_any_ab_tcr[rows_all] | obs.corrected_has_any_gd_tcr[rows_all]
    source_has_tcrseq = (
        pd.DataFrame({"source_gse_id": source, "tcr_evidence": tcr_evidence.astype(np.int8)})
        .groupby("source_gse_id", sort=False)["tcr_evidence"]
        .transform("max")
        .to_numpy(dtype=bool)
    )
    trdc_flag = np.zeros(n_eval, dtype=bool)
    trdv_flag = np.zeros(n_eval, dtype=bool)

    strategies = ["v2_high_f1", "v2_high_purity", "original_TRD_minus_TRAB", *[c.name for c in candidates]]
    scores = {strategy: np.zeros(n_eval, dtype=np.float32) for strategy in strategies}
    preds = {strategy: np.zeros(n_eval, dtype=bool) for strategy in strategies}
    thresholds: dict[str, float | str] = {"original_TRD_minus_TRAB": float(load_original_trd_trab_threshold())}
    base = v2_payload["base_model"]
    v2_cols = [spec.gene_to_col[str(gene)] for gene in base["gene_names"]]
    original_full = read_float_obs(handle, "phase4_trd_minus_trab")[:n_eval].astype(np.float32) if "phase4_trd_minus_trab" in handle["obs"] else np.zeros(n_eval, dtype=np.float32)
    scores["original_TRD_minus_TRAB"] = original_full
    preds["original_TRD_minus_TRAB"] = original_full >= float(thresholds["original_TRD_minus_TRAB"])

    for start in range(0, n_eval, args.chunk_size):
        end = min(start + args.chunk_size, n_eval)
        rows = np.arange(start, end, dtype=np.int64)
        x_gene, _row_sum, _n_detected = extract_gene_features(handle, "X", rows, spec, label=f"full_atlas_{start}_{end}")
        x = append_engineered_features(x_gene, spec)
        trdc_flag[start:end] = x[:, spec.engineered_to_col["TRDC_log1p"]] > 0
        trdv_flag[start:end] = x[:, spec.engineered_to_col["any_TRDV"]] > 0.5
        v2_score = base["model_object"].predict_proba(x[:, v2_cols])[:, 1].astype(np.float32)
        for mode in ["high_f1", "high_purity"]:
            key = f"v2_{mode}"
            threshold_vec, _threshold_label = v2_threshold_vector(v2_payload, mode, annotation[start:end])
            scores[key][start:end] = v2_score
            preds[key][start:end] = v2_score >= threshold_vec
            thresholds[key] = "annotation_specific" if mode == "high_purity" else float(v2_payload["operating_modes"][mode]["threshold"])
        candidate_score_cache: dict[int, np.ndarray] = {}
        for candidate in candidates:
            cache_key = id(candidate.model_object)
            if cache_key not in candidate_score_cache:
                candidate_score_cache[cache_key] = candidate.model_object.predict_proba(x)[:, 1].astype(np.float32)
            score = candidate_score_cache[cache_key]
            scores[candidate.name][start:end] = score
            preds[candidate.name][start:end] = score >= candidate.threshold
            thresholds[candidate.name] = float(candidate.threshold)
        if end % 500_000 == 0 or end == n_eval:
            logging.info("Applied v3 strategies to %s / %s full-atlas cells", f"{end:,}", f"{n_eval:,}")

    quadrant = quadrant_from_flags(trdc_flag, trdv_flag)
    overall_parts = []
    source_parts = []
    tissue_parts = []
    annotation_parts = []
    for strategy in strategies:
        overall, by_source, by_tissue, by_annotation = full_group_summary(
            strategy=strategy,
            threshold=thresholds.get(strategy, np.nan),
            pred=preds[strategy],
            score=scores[strategy],
            y=y,
            primary_mask=primary_mask,
            source=source,
            tissue=tissue,
            annotation=annotation,
            quadrant=quadrant,
            is_nk=is_nk,
            is_tcrab=is_tcrab,
            is_cd4_treg=is_cd4_treg,
            is_b_myeloid=is_b_myeloid,
            source_has_tcrseq=source_has_tcrseq,
        )
        overall_parts.append(overall)
        source_parts.append(by_source)
        tissue_parts.append(by_tissue)
        annotation_parts.append(by_annotation)
    full_overall = pd.concat(overall_parts, ignore_index=True)
    full_source = pd.concat(source_parts, ignore_index=True)
    full_tissue = pd.concat(tissue_parts, ignore_index=True)
    full_annotation = pd.concat(annotation_parts, ignore_index=True)
    full_overall.to_csv(TABLE_DIR / "full_atlas_prediction_overall.csv", index=False)
    full_source.to_csv(TABLE_DIR / "full_atlas_prediction_by_source.csv", index=False)
    full_tissue.to_csv(TABLE_DIR / "full_atlas_prediction_by_tissue.csv", index=False)
    full_annotation.to_csv(TABLE_DIR / "full_atlas_prediction_by_annotation.csv", index=False)

    selected_for_full = select_candidate_from_full_atlas(candidates, full_overall, args)
    if selected_for_full.name != selected.name:
        logging.info("Full-atlas target selector changed selected candidate from %s to %s", selected.name, selected_for_full.name)
    pd.DataFrame(
        [
            {
                "selected_strategy": selected_for_full.name,
                "previous_internal_selected_strategy": selected.name,
                "target_recall": float(args.target_recall),
                "target_recall_margin": float(args.target_recall_margin),
                "target_max_estimated_fp_fraction": float(args.max_estimated_fp_fraction),
                "fp_comparator_strategy": str(getattr(args, "fp_comparator_strategy", "v2_high_purity")),
                "target_max_TRDC_plus_TRDV_minus_fraction": float(args.max_trdc_trdv_minus_fraction),
            }
        ]
    ).to_csv(TABLE_DIR / "full_atlas_selected_strategy.csv", index=False)

    selected_pred = preds[selected_for_full.name]
    selected_rows = rows_all[selected_pred]
    selected_cells = pd.DataFrame(
        {
            "obs_index": selected_rows,
            "cell_id": obs_index_values(handle)[selected_rows],
            "source_gse_id": source[selected_pred],
            "tissue": tissue[selected_pred],
            "annotation": annotation[selected_pred],
            "class_code": obs.class_code[selected_rows].astype(int),
            "score": scores[selected_for_full.name][selected_pred],
            "threshold": float(selected_for_full.threshold),
            "tcr_gene_quadrant": quadrant[selected_pred],
            "is_NK_annotation": is_nk[selected_pred],
            "is_TCRAB_no_gdT_gold": is_tcrab[selected_pred],
        }
    )
    selected_cells.to_csv(TABLE_DIR / "full_atlas_selected_predicted_cells.csv.gz", index=False, compression="gzip")
    return {
        "overall": full_overall,
        "source": full_source,
        "tissue": full_tissue,
        "annotation": full_annotation,
        "selected_cells": selected_cells,
        "selected_strategy": pd.DataFrame([{"strategy": selected_for_full.name}]),
        "target_selection": pd.read_csv(TABLE_DIR / "full_atlas_target_selection.csv"),
    }

def bootstrap_coefficients(
    x_train: np.ndarray,
    y_train: np.ndarray,
    sample_weight: np.ndarray,
    feature_names: list[str],
    *,
    seed: int,
    iterations: int,
    sample_cells: int,
) -> pd.DataFrame:
    if iterations <= 0:
        out = pd.DataFrame()
        out.to_csv(TABLE_DIR / "elasticnet_bootstrap_coefficient_stability.csv", index=False)
        return out
    rng = np.random.default_rng(seed + 9)
    rows = []
    n = x_train.shape[0]
    size = min(sample_cells, n)
    for i in range(iterations):
        idx = rng.choice(n, size=size, replace=True)
        model = make_pipeline(
            StandardScaler(),
            SGDClassifier(
                loss="log_loss",
                penalty="elasticnet",
                l1_ratio=0.35,
                alpha=1e-4,
                class_weight="balanced",
                max_iter=10,
                tol=1e-3,
                n_jobs=16,
                random_state=seed + 100 + i,
            ),
        )
        model.fit(x_train[idx], y_train[idx], sgdclassifier__sample_weight=sample_weight[idx])
        coefs = model.named_steps["sgdclassifier"].coef_[0]
        top = np.argsort(np.abs(coefs))[-40:]
        for j in top:
            rows.append({"bootstrap_iteration": i, "feature": feature_names[j], "coefficient": float(coefs[j]), "abs_coefficient": float(abs(coefs[j]))})
        logging.info("Bootstrap coefficient fit %s / %s complete", i + 1, iterations)
    df = pd.DataFrame(rows)
    if not df.empty:
        summary = (
            df.groupby("feature", as_index=False)
            .agg(
                n_selected=("coefficient", "size"),
                mean_coefficient=("coefficient", "mean"),
                sd_coefficient=("coefficient", "std"),
                positive_fraction=("coefficient", lambda s: float((s > 0).mean())),
                mean_abs_coefficient=("abs_coefficient", "mean"),
            )
            .sort_values(["n_selected", "mean_abs_coefficient"], ascending=False)
        )
    else:
        summary = pd.DataFrame()
    summary.to_csv(TABLE_DIR / "elasticnet_bootstrap_coefficient_stability.csv", index=False)
    return summary


def save_best_model(candidate: CandidateResult, spec: FeatureSpec, external_acceptance: pd.DataFrame) -> tuple[Path, Path | None]:
    accepted = bool(external_acceptance["accepted_for_promotion"].iloc[0])
    payload = {
        "model": candidate.name,
        "threshold": float(candidate.threshold),
        "feature_names": spec.model_feature_names,
        "gene_names": spec.gene_names,
        "engineered_feature_names": spec.engineered_feature_names,
        "model_object": candidate.model_object,
        "notes": candidate.notes,
        "version": "gdTAI_v3_candidate",
        "accepted_for_promotion": accepted,
        "external_acceptance": external_acceptance.to_dict(orient="records")[0],
        "requires_counts_layer_for_external_h5ad": True,
    }
    best_path = MODEL_DIR / "best_candidate_model.pkl"
    with best_path.open("wb") as out:
        pickle.dump(payload, out)
    selected_path = None
    if accepted:
        selected_path = PROMOTED_DIR / "gdTAI_v3_model.pkl"
        with selected_path.open("wb") as out:
            promoted = dict(payload)
            promoted["version"] = "gdTAI_v3.0"
            pickle.dump(promoted, out)
        pd.DataFrame({"feature": spec.model_feature_names, "feature_index": np.arange(len(spec.model_feature_names))}).to_csv(PROMOTED_DIR / "feature_genes.csv", index=False)
        external_acceptance.to_csv(PROMOTED_DIR / "mode_metrics.csv", index=False)
        (PROMOTED_DIR / "README.md").write_text(
            "# gdTAI v3.0\n\n"
            "Promoted TRDC/NK-guard gdT classifier. External H5AD inference must use raw counts and rebuild log1p CP10K features.\n",
            encoding="utf-8",
        )
        (PROMOTED_DIR / "METHODOLOGY.md").write_text(
            "# gdTAI v3.0 Methodology\n\n"
            "- GSE144469 is included in atlas train/tune splits.\n"
            "- The independent external H5AD is used only for final testing.\n"
            "- Features are TCR genes, CD3/NK lineage controls, and conditional TRDC/TRDV/TRDJ features.\n"
            "- NK markers are conditional risk features, not standalone exclusion rules.\n"
            "- NK+TCRAB-overlap cells are excluded from negative training/tuning labels; explicit NK hard negatives require expression TRDC+TRDV-.\n",
            encoding="utf-8",
        )
        examples = PROMOTED_DIR / "examples"
        examples.mkdir(exist_ok=True)
        (examples / "run_external_test.sh").write_text(
            "#!/usr/bin/env bash\n"
            "set -euo pipefail\n"
            "/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python workflows/gdtai/run_gdtai_v3_trdc_nk_guard_classifier.py --external-h5ad \"$1\"\n",
            encoding="utf-8",
        )
    return best_path, selected_path


def promotion_acceptance(
    selected: CandidateResult,
    external_metrics: pd.DataFrame,
    fp: pd.DataFrame,
    recall: pd.DataFrame,
    full_atlas_target_selection: pd.DataFrame | None = None,
) -> pd.DataFrame:
    v3 = external_metrics.loc[external_metrics["strategy"] == selected.name].iloc[0]
    v2hp = external_metrics.loc[external_metrics["strategy"] == "v2_high_purity"].iloc[0]
    fp_pivot = fp.pivot_table(index="strategy", columns="false_positive_group", values="false_positive_cells", aggfunc="sum").fillna(0)
    rec_pivot = recall.pivot_table(index="strategy", columns="recall_group", values="recall", aggfunc="first").fillna(np.nan)

    full_row = pd.Series(dtype=object)
    if full_atlas_target_selection is not None and not full_atlas_target_selection.empty:
        full_hits = full_atlas_target_selection.loc[full_atlas_target_selection["strategy"].astype(str) == selected.name]
        if not full_hits.empty:
            full_row = full_hits.iloc[0]
    full_recall = float(full_row.get("full_primary_recall", np.nan)) if not full_row.empty else np.nan
    full_estimated_fp = float(full_row.get("estimated_fp_fraction_of_predictions", np.nan)) if not full_row.empty else np.nan
    full_effective_max_fp = float(full_row.get("effective_max_estimated_fp_fraction", np.nan)) if not full_row.empty else np.nan
    full_trdc_fraction = float(full_row.get("predicted_TRDC_plus_TRDV_minus_fraction", np.nan)) if not full_row.empty else np.nan
    full_target_met = bool(full_row.get("target_met_full_atlas", False)) if not full_row.empty else False

    checks = {
        "full_atlas_active_goal_gates_met": full_target_met,
        "external_NK_FP_not_above_v2_high_purity": float(fp_pivot.loc[selected.name, "NK_cells"]) <= float(fp_pivot.loc["v2_high_purity", "NK_cells"]),
        "external_TRDC_plus_TRDV_minus_burden_lower_than_v2_high_purity": float(fp_pivot.loc[selected.name, "TRDC_plus_TRDV_minus_cells"]) < float(fp_pivot.loc["v2_high_purity", "TRDC_plus_TRDV_minus_cells"]),
        "external_F1_not_worse_than_v2_high_purity_by_gt_0p01": float(v3["f1"]) >= float(v2hp["f1"]) - 0.01,
        "cytotoxic_gdT_recall_not_substantially_degraded": float(rec_pivot.loc[selected.name, "cytotoxic_NK_marker_high_gdT"]) >= float(rec_pivot.loc["v2_high_purity", "cytotoxic_NK_marker_high_gdT"]) - 0.03,
        "paired_TCRAB_FP_not_above_v2_high_purity": float(fp_pivot.loc[selected.name, "paired_TCRAB_cells"]) <= float(fp_pivot.loc["v2_high_purity", "paired_TCRAB_cells"]),
    }
    row = {
        "candidate_model": selected.name,
        "accepted_for_promotion": all(checks.values()),
        "full_atlas_recall": full_recall,
        "full_atlas_estimated_fp_fraction": full_estimated_fp,
        "full_atlas_effective_max_estimated_fp_fraction": full_effective_max_fp,
        "full_atlas_TRDC_plus_TRDV_minus_fraction": full_trdc_fraction,
        "external_f1": float(v3["f1"]),
        "v2_high_purity_external_f1": float(v2hp["f1"]),
        "external_precision": float(v3["precision"]),
        "external_recall": float(v3["recall"]),
        "external_specificity": float(v3["specificity"]),
        "external_NK_FP": int(fp_pivot.loc[selected.name, "NK_cells"]),
        "v2_high_purity_NK_FP": int(fp_pivot.loc["v2_high_purity", "NK_cells"]),
        "external_TRDC_plus_TRDV_minus_FP": int(fp_pivot.loc[selected.name, "TRDC_plus_TRDV_minus_cells"]),
        "v2_high_purity_TRDC_plus_TRDV_minus_FP": int(fp_pivot.loc["v2_high_purity", "TRDC_plus_TRDV_minus_cells"]),
        "external_paired_TCRAB_FP": int(fp_pivot.loc[selected.name, "paired_TCRAB_cells"]),
        "v2_high_purity_paired_TCRAB_FP": int(fp_pivot.loc["v2_high_purity", "paired_TCRAB_cells"]),
        **checks,
    }
    out = pd.DataFrame([row])
    out.to_csv(TABLE_DIR / "promotion_acceptance_gates.csv", index=False)
    return out


def write_report(
    *,
    selected: CandidateResult,
    candidate_metrics: pd.DataFrame,
    external_metrics: pd.DataFrame,
    fp: pd.DataFrame,
    recall: pd.DataFrame,
    acceptance: pd.DataFrame,
    split_overall: pd.DataFrame,
    hard_negative_filter_summary: pd.DataFrame,
    atlas_validation_metrics: pd.DataFrame,
    atlas_validation_groups: pd.DataFrame,
    full_atlas_overall: pd.DataFrame,
    full_atlas_target_selection: pd.DataFrame,
    coefficient_stability: pd.DataFrame,
    figures: list[Path],
    best_model: Path,
    promoted_model: Path | None,
) -> None:
    accepted = bool(acceptance["accepted_for_promotion"].iloc[0])
    selected_external = external_metrics.loc[external_metrics["strategy"] == selected.name].iloc[0]
    lines = [
        "# gdTAI v3 TRDC/NK Guard Report",
        "",
        f"- Selected internal candidate: `{selected.name}`",
        f"- Promotion accepted: `{accepted}`",
        f"- Internal threshold: `{selected.threshold:.8f}`",
        f"- External precision/recall/F1: `{selected_external['precision']:.4f}` / `{selected_external['recall']:.4f}` / `{selected_external['f1']:.4f}`",
        f"- Best candidate model: `{best_model}`",
        f"- Promoted gdTAI v3.0 model: `{promoted_model if promoted_model else 'not written because gates failed'}`",
        "",
        "## Scope",
        "",
        "- GSE144469 is included in atlas train/tune splits and is no longer a holdout.",
        "- Non-GSE144469 validation cohorts remain held out: GSE254249 paired-TCRAB/no-gdTCR negatives and GDT_2020AUG_woCOV cord-blood gdT positives.",
        "- The final selected v3 candidate is also applied to the full atlas input H5AD for whole-dataset evaluation.",
        "- The independent external H5AD is used only after model fitting and threshold selection.",
        "- External features are computed from `layers[\"counts\"]`; normalized `X` is not used.",
        "- NK+TCRAB-overlap cells are excluded from negative train/tune labels; explicit NK hard negatives must be expression `TRDC+TRDV-`.",
        "",
        "## Split Summary",
        dataframe_to_markdown(split_overall),
        "",
        "## Hard Negative Expression Filter",
        dataframe_to_markdown(hard_negative_filter_summary),
        "",
        "## Internal Tune Candidates",
        dataframe_to_markdown(candidate_metrics),
        "",
        "## Atlas Validation Metrics",
        dataframe_to_markdown(atlas_validation_metrics),
        "",
        "## Atlas Validation Metrics By Group",
        dataframe_to_markdown(atlas_validation_groups.head(80) if not atlas_validation_groups.empty else atlas_validation_groups),
        "",
        "## Full-Atlas Application Summary",
        dataframe_to_markdown(full_atlas_overall),
        "",
        "## Full-Atlas Target Selection",
        dataframe_to_markdown(full_atlas_target_selection),
        "",
        "## External Primary Metrics",
        dataframe_to_markdown(external_metrics),
        "",
        "## Promotion Gates",
        dataframe_to_markdown(acceptance),
        "",
        "## External FP Stress Groups",
        dataframe_to_markdown(fp),
        "",
        "## External Recall Stress Groups",
        dataframe_to_markdown(recall),
        "",
        "## Coefficient Stability",
        dataframe_to_markdown(coefficient_stability.head(40) if not coefficient_stability.empty else coefficient_stability),
        "",
        "## Outputs",
        f"- Tables: `{TABLE_DIR}`",
        f"- Figures: `{FIGURE_DIR}`",
        f"- Logs: `{LOG_DIR}`",
        f"- HTML report: `{REPORT_HTML}`",
        f"- PDF report: `{REPORT_PDF}`",
    ]
    REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    asset_dir = STATIC_DIR / "assets" / OUT_PREFIX
    asset_dir.mkdir(parents=True, exist_ok=True)
    fig_html = []
    for fig in figures:
        target = asset_dir / fig.name
        if fig.exists() and fig.resolve() != target.resolve():
            shutil.copyfile(fig, target)
        fig_html.append(f"<figure><img src='assets/{OUT_PREFIX}/{html.escape(fig.name)}'><figcaption>{html.escape(fig.stem.replace('_', ' '))}</figcaption></figure>")
    css = """
    @page{size:A4 landscape;margin:8mm}body{font-family:Arial,Helvetica,sans-serif;margin:18px;color:#20272e;background:#f5f6f7;line-height:1.45}main{max-width:1320px;margin:auto}section{background:white;border:1px solid #d9dee4;padding:14px;margin:12px 0;break-inside:avoid}h1{font-size:28px;margin:0 0 8px}h2{font-size:19px}.grid{display:grid;grid-template-columns:repeat(5,1fr);gap:10px}.metric{border:1px solid #d9dee4;padding:10px;background:white}.value{font-size:22px;font-weight:bold}table{border-collapse:collapse;width:100%;font-size:9px;table-layout:fixed}th,td{border:1px solid #d9dee4;padding:3px 4px;text-align:left;vertical-align:top;overflow-wrap:anywhere}th{background:#eef1f4}img{max-width:100%;border:1px solid #d9dee4}.status{font-weight:bold;color:__COLOR__}code{background:#eef1f4;padding:1px 4px;border-radius:3px}
    """.replace("__COLOR__", "#0f766e" if accepted else "#be123c")
    html_doc = f"""<!doctype html><html><head><meta charset='utf-8'><title>gdTAI v3 TRDC/NK Guard</title><style>{css}</style></head><body><main>
    <section><h1>gdTAI v3 TRDC/NK Guard Report</h1><p class='status'>Promotion accepted: {accepted}</p><p>GSE144469 is train/tune eligible. The external H5AD is held independent and read from <code>layers["counts"]</code>.</p></section>
    <section><h2>Headline Metrics</h2><div class='grid'>
      <div class='metric'><div>Selected</div><div class='value'>{html.escape(selected.name)}</div></div>
      <div class='metric'><div>External F1</div><div class='value'>{selected_external['f1']:.3f}</div></div>
      <div class='metric'><div>Precision</div><div class='value'>{selected_external['precision']:.3f}</div></div>
      <div class='metric'><div>Recall</div><div class='value'>{selected_external['recall']:.3f}</div></div>
      <div class='metric'><div>Specificity</div><div class='value'>{selected_external['specificity']:.3f}</div></div>
    </div></section>
    <section><h2>Promotion Gates</h2>{dataframe_to_html(acceptance)}</section>
    <section><h2>External Primary Metrics</h2>{dataframe_to_html(external_metrics)}</section>
    <section><h2>External FP Stress Groups</h2>{dataframe_to_html(fp)}</section>
    <section><h2>External Recall Groups</h2>{dataframe_to_html(recall)}</section>
    <section><h2>Internal Tune Candidates</h2>{dataframe_to_html(candidate_metrics)}</section>
    <section><h2>Atlas Validation Metrics</h2>{dataframe_to_html(atlas_validation_metrics)}</section>
    <section><h2>Full-Atlas Application Summary</h2>{dataframe_to_html(full_atlas_overall)}</section>
    <section><h2>Full-Atlas Target Selection</h2>{dataframe_to_html(full_atlas_target_selection)}</section>
    <section><h2>Split Summary</h2>{dataframe_to_html(split_overall)}</section>
    <section><h2>Hard Negative Expression Filter</h2>{dataframe_to_html(hard_negative_filter_summary)}</section>
    <section><h2>Figures</h2>{''.join(fig_html)}</section>
    </main></body></html>"""
    REPORT_HTML.write_text(html_doc, encoding="utf-8")


def render_pdf(no_pdf: bool) -> None:
    if no_pdf:
        return
    subprocess.run(
        [
            "google-chrome",
            "--headless",
            "--disable-gpu",
            "--no-sandbox",
            "--print-to-pdf-no-header",
            f"--print-to-pdf={REPORT_PDF.resolve()}",
            REPORT_HTML.resolve().as_uri(),
        ],
        check=True,
    )


def main() -> None:
    args = parse_args()
    setup_logging()
    logging.info("Atlas training H5AD: %s", args.input_h5ad)
    logging.info("Independent external H5AD: %s", args.external_h5ad)
    v2_payload = load_pickle(args.v2_model_pkl)
    if "base_model" not in v2_payload:
        raise KeyError("Expected gdTAI v2 wrapper payload with `base_model`.")

    with h5py.File(args.input_h5ad, "r") as handle:
        obs, split = make_splits(handle, args)
        spec = build_feature_spec(handle)
        eval_rows = np.unique(
            np.concatenate(
                [
                    split.train_idx,
                    split.tune_idx,
                    split.hard_negative_train_idx,
                    split.hard_negative_tune_idx,
                    split.weak_trdc_prefilter_idx,
                    split.validation_idx,
                ]
            )
        ).astype(np.int64)
        logging.info("Extracting atlas v3 features for %s train/tune/prefilter rows", f"{eval_rows.size:,}")
        x_gene, _row_sum, _n_detected = extract_gene_features(handle, "X", eval_rows, spec, label="atlas_train_tune")
        x_eval = append_engineered_features(x_gene, spec)
        weak_pos = local_positions(eval_rows, split.weak_trdc_prefilter_idx) if split.weak_trdc_prefilter_idx.size else np.asarray([], dtype=np.int64)
        weak_keep_prefilter = weak_pos[
            (x_eval[weak_pos, spec.engineered_to_col["TRDC_only"]] > 0.5)
            & (x_eval[weak_pos, spec.engineered_to_col["any_TRDV"]] < 0.5)
            & (x_eval[weak_pos, spec.engineered_to_col["any_TRDJ"]] < 0.5)
            & (x_eval[weak_pos, spec.engineered_to_col["NK_minus_CD3_score"]] > 0.20)
        ]
        weak_rows = eval_rows[weak_keep_prefilter]
        weak_nk_tcrab_overlap = split.nk_guard_mask[weak_rows] & split.tcrab_guard_mask[weak_rows]
        weak_keep = weak_keep_prefilter[~weak_nk_tcrab_overlap]
        y_eval = np.zeros(eval_rows.size, dtype=np.int8)
        y_eval[obs.class_code[eval_rows] == 2] = 1
        hard_train_prefilter_pos = local_positions(eval_rows, split.hard_negative_train_idx)
        hard_tune_prefilter_pos = local_positions(eval_rows, split.hard_negative_tune_idx)
        hard_train_filtered, hard_train_filter_stats = filter_hard_negative_positions(
            eval_rows, x_eval, spec, split, hard_train_prefilter_pos, "train_hard_negative"
        )
        hard_tune_filtered, hard_tune_filter_stats = filter_hard_negative_positions(
            eval_rows, x_eval, spec, split, hard_tune_prefilter_pos, "tune_hard_negative"
        )
        hard_negative_filter_summary = pd.DataFrame([hard_train_filter_stats, hard_tune_filter_stats])
        hard_negative_filter_summary.to_csv(TABLE_DIR / "hard_negative_expression_filter_summary.csv", index=False)
        hard_negative_filter_by_source(
            eval_rows, obs, split, hard_train_filtered, hard_tune_filtered
        ).to_csv(TABLE_DIR / "hard_negative_expression_filter_by_source.csv", index=False)
        positions = {
            "primary_train": local_positions(eval_rows, split.train_idx),
            "primary_tune": local_positions(eval_rows, split.tune_idx),
            "hard_train": np.unique(np.concatenate([hard_train_filtered, weak_keep])).astype(np.int64),
            "hard_tune": hard_tune_filtered,
            "validation": local_positions(eval_rows, split.validation_idx),
        }
        pd.DataFrame(
            [
                {
                    "weak_trdc_prefilter_cells": int(weak_pos.size),
                    "weak_trdc_expression_prefilter_kept": int(weak_keep_prefilter.size),
                    "weak_trdc_excluded_NK_TCRAB_overlap": int(weak_nk_tcrab_overlap.sum()),
                    "weak_trdc_kept_as_hard_negative": int(weak_keep.size),
                }
            ]
        ).to_csv(TABLE_DIR / "weak_trdc_hard_negative_filter_summary.csv", index=False)
        candidates = train_candidates(x_eval, y_eval, positions, spec, args, v2_payload)
        candidate_metrics = candidate_metrics_frame(candidates)
        selected = select_internal_candidate(candidates, args)
        logging.info("Selected internal candidate: %s", selected.name)
        internal_tune_df = write_internal_tune_predictions(x_eval, y_eval, eval_rows, positions, candidates, spec, obs, split)
        phase4_minus_eval = read_float_obs(handle, "phase4_trd_minus_trab")[eval_rows] if "phase4_trd_minus_trab" in handle["obs"] else np.zeros(eval_rows.size, dtype=np.float32)
        atlas_validation_df, atlas_validation_metrics, atlas_validation_groups = write_atlas_validation_predictions(
            x_eval, y_eval, eval_rows, positions, candidates, spec, obs, split, v2_payload, phase4_minus_eval
        )
        train_sample, train_weights = sample_train_positions(positions["primary_train"], positions["hard_train"], y_eval, x_eval, spec, args)
        coefficient_stability = pd.DataFrame()
        if not args.skip_bootstrap:
            coefficient_stability = bootstrap_coefficients(
                x_eval[train_sample],
                y_eval[train_sample],
                train_weights,
                spec.model_feature_names,
                seed=args.seed,
                iterations=args.bootstrap_iterations,
                sample_cells=args.bootstrap_sample_cells,
            )
        full_atlas_outputs = {}
        if not args.no_full_atlas:
            full_atlas_outputs = apply_full_atlas(handle, spec, candidates, selected, v2_payload, obs, split, args)
            selected_strategy_df = full_atlas_outputs.get("selected_strategy", pd.DataFrame())
            if not selected_strategy_df.empty:
                selected_name = str(selected_strategy_df["strategy"].iloc[0])
                candidate_lookup = {candidate.name: candidate for candidate in candidates}
                if selected_name in candidate_lookup:
                    selected = candidate_lookup[selected_name]
                    logging.info("Using full-atlas target-selected candidate for external test and saved model: %s", selected.name)

    with h5py.File(args.external_h5ad, "r") as handle:
        if "layers" not in handle or "counts" not in handle["layers"]:
            raise RuntimeError("External H5AD must contain layers['counts']; refusing to use normalized/log X.")
        var_names = pd.Index(read_string_dataset(handle["var"]["_index"]), dtype="string").astype(str).tolist()
        availability = selected_gene_mapping(var_names, spec.gene_names)
        missing = availability.loc[~availability["available_in_h5ad"], "gene"].astype(str).tolist()
        if missing:
            raise KeyError(f"External H5AD is missing v3 feature genes: {missing[:20]}")
        external_gene_to_idx = {gene: idx for idx, gene in enumerate(var_names)}
        external_spec = FeatureSpec(
            gene_names=spec.gene_names,
            gene_indices=np.asarray([external_gene_to_idx[gene] for gene in spec.gene_names], dtype=np.int32),
            gene_feature_names=spec.gene_feature_names,
            engineered_feature_names=spec.engineered_feature_names,
            model_feature_names=spec.model_feature_names,
            gene_to_col=spec.gene_to_col,
            engineered_to_col=spec.engineered_to_col,
        )
        availability.to_csv(TABLE_DIR / "external_feature_gene_availability.csv", index=False)
        n_external = h5ad_shape(handle, "layers/counts")[0]
        rows = np.arange(n_external, dtype=np.int64)
        if args.max_external_cells is not None:
            rows = rows[: int(args.max_external_cells)]
            logging.warning("Smoke mode: evaluating only first %s external cells", rows.size)
        logging.info("Extracting external features from layers['counts'] for %s cells", f"{rows.size:,}")
        x_gene_ext, row_sum, n_detected = extract_gene_features(handle, "layers/counts", rows, external_spec, label="external_counts")
        x_ext = append_engineered_features(x_gene_ext, external_spec)
        external_df = external_truth(handle, rows)
        external_df["row_sum_counts_layer"] = row_sum
        external_df["n_detected_genes_counts_layer"] = n_detected
        for feature in ["any_TRDV", "any_TRDJ", "any_TRG", "TRDC_only", "CD3_score", "NK_score", "gdT_TCR_score", "abT_TCR_score", "NK_minus_CD3_score"]:
            external_df[feature] = x_ext[:, external_spec.engineered_to_col[feature]]
        external_df["tcr_gene_quadrant"] = trdc_trdv_quadrant(x_ext, external_spec)

        strategies: list[str] = []
        base = v2_payload["base_model"]
        v2_cols = [external_spec.gene_to_col[str(gene)] for gene in base["gene_names"]]
        v2_score = base["model_object"].predict_proba(x_ext[:, v2_cols])[:, 1].astype(np.float32)
        for mode in ["high_f1", "high_purity"]:
            mode_key = f"v2_{mode}"
            threshold, threshold_label = v2_threshold_vector(v2_payload, mode, external_df["cell_type"].to_numpy(dtype=object))
            external_df[f"{mode_key}_score"] = v2_score
            external_df[f"{mode_key}_threshold"] = threshold
            external_df[f"{mode_key}_threshold_label"] = threshold_label
            external_df[f"{mode_key}_pred"] = v2_score >= threshold
            strategies.append(mode_key)

        original_threshold = load_original_trd_trab_threshold()
        original_score = external_df["phase4_trd_minus_trab"].to_numpy(dtype=np.float32)
        external_df["original_TRD_minus_TRAB_score"] = original_score
        external_df["original_TRD_minus_TRAB_threshold"] = original_threshold
        external_df["original_TRD_minus_TRAB_pred"] = original_score >= original_threshold
        strategies.append("original_TRD_minus_TRAB")

        external_score_cache: dict[int, np.ndarray] = {}
        for candidate in candidates:
            cache_key = id(candidate.model_object)
            if cache_key not in external_score_cache:
                external_score_cache[cache_key] = candidate.model_object.predict_proba(x_ext)[:, 1].astype(np.float32)
            score = external_score_cache[cache_key]
            external_df[f"{candidate.name}_score"] = score
            external_df[f"{candidate.name}_threshold"] = candidate.threshold
            external_df[f"{candidate.name}_pred"] = score >= candidate.threshold
            strategies.append(candidate.name)

    prediction_path = TABLE_DIR / "external_predictions_wide.csv.gz"
    external_df.to_csv(prediction_path, index=False, compression="gzip")
    external_metrics, fp, recall, group_df = evaluate_external_predictions(external_df, strategies)
    acceptance = promotion_acceptance(selected, external_metrics, fp, recall, full_atlas_outputs.get("target_selection", pd.DataFrame()))
    figures = plot_external_summary(external_metrics, fp, recall)
    best_model, promoted_model = save_best_model(selected, spec, acceptance)
    write_report(
        selected=selected,
        candidate_metrics=candidate_metrics,
        external_metrics=external_metrics,
        fp=fp,
        recall=recall,
        acceptance=acceptance,
        split_overall=split.split_overall,
        hard_negative_filter_summary=hard_negative_filter_summary,
        atlas_validation_metrics=atlas_validation_metrics,
        atlas_validation_groups=atlas_validation_groups,
        full_atlas_overall=full_atlas_outputs.get("overall", pd.DataFrame()),
        full_atlas_target_selection=full_atlas_outputs.get("target_selection", pd.DataFrame()),
        coefficient_stability=coefficient_stability,
        figures=figures,
        best_model=best_model,
        promoted_model=promoted_model,
    )
    render_pdf(args.no_pdf)
    SUMMARY_JSON.write_text(
        json.dumps(
            {
                "atlas_input_h5ad": str(args.input_h5ad),
                "external_h5ad": str(args.external_h5ad),
                "selected_candidate": selected.name,
                "selected_threshold": float(selected.threshold),
                "accepted_for_promotion": bool(acceptance["accepted_for_promotion"].iloc[0]),
                "external_predictions": str(prediction_path),
                "atlas_validation_metrics": str(TABLE_DIR / "atlas_validation_metrics.csv"),
                "full_atlas_prediction_overall": str(TABLE_DIR / "full_atlas_prediction_overall.csv") if full_atlas_outputs else None,
                "full_atlas_prediction_by_source": str(TABLE_DIR / "full_atlas_prediction_by_source.csv") if full_atlas_outputs else None,
                "full_atlas_target_selection": str(TABLE_DIR / "full_atlas_target_selection.csv") if full_atlas_outputs else None,
                "full_atlas_selected_strategy": str(TABLE_DIR / "full_atlas_selected_strategy.csv") if full_atlas_outputs else None,
                "full_atlas_selected_predicted_cells": str(TABLE_DIR / "full_atlas_selected_predicted_cells.csv.gz") if full_atlas_outputs else None,
                "report_html": str(REPORT_HTML),
                "report_pdf": str(REPORT_PDF),
                "best_candidate_model": str(best_model),
                "promoted_model": None if promoted_model is None else str(promoted_model),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    logging.info("Saved v3 report: %s", REPORT_HTML)


if __name__ == "__main__":
    main()
