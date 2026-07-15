#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Compare local gdTAI v3.0 latest predictions against gdTAI v2.0 high-purity.

This script is read-only with respect to H5AD files. It applies both models to
all cells in the integrated atlas, saves cell-level overlap for the union of
predictions, and aggregates the differences by dataset/source, tissue, and
simple annotation.
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

from run_gdt_gse144469_holdout_tcrgene_classifier import build_obs_metadata
from run_gdt_prediction_package_evaluation import (
    DEFAULT_INPUT_H5AD,
    clean_group_values,
    dataframe_to_html,
    dataframe_to_markdown,
    read_bool_obs,
    read_float_obs,
    read_obs_column,
)
from run_gdtai_v3_trdc_nk_guard_classifier import (
    ENGINEERED_FEATURES,
    FeatureSpec,
    append_engineered_features,
    extract_gene_features,
    obs_index_values,
    quadrant_from_flags,
    safe_div,
    v2_threshold_vector,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
STATIC_DIR = PROJECT_ROOT / "gdT_prediction"
OUT_PREFIX = "gdtai_v3_vs_v2_high_purity"
TABLE_DIR = OUTPUT_ROOT / "tables" / "gdT_prediction" / OUT_PREFIX
FIGURE_DIR = OUTPUT_ROOT / "figures" / "gdT_prediction" / OUT_PREFIX
LOG_DIR = OUTPUT_ROOT / "logs" / "gdT_prediction" / OUT_PREFIX
STATIC_ASSET_DIR = STATIC_DIR / "assets" / OUT_PREFIX
REPORT_MD = LOG_DIR / "gdtai_v3_vs_v2_high_purity_report.md"
REPORT_HTML = STATIC_DIR / "gdtai_v3_vs_v2_high_purity_report.html"
REPORT_PDF = STATIC_DIR / "gdtai_v3_vs_v2_high_purity_report.pdf"
RUN_LOG = LOG_DIR / "run.log"
SUMMARY_JSON = LOG_DIR / "gdtai_v3_vs_v2_high_purity_summary.json"

V3_MODEL = OUTPUT_ROOT / "models" / "gdT_prediction_classifier" / "gdTAI_v3.0" / "gdTAI_v3_model.pkl"
V3_MANIFEST = OUTPUT_ROOT / "models" / "gdT_prediction_classifier" / "gdTAI_v3.0" / "model_manifest.json"
V2_MODEL = OUTPUT_ROOT / "models" / "gdT_prediction_classifier" / "gdTAI_v2.0" / "gdTAI_v2_model.pkl"
V2_MANIFEST = OUTPUT_ROOT / "models" / "gdT_prediction_classifier" / "gdTAI_v2.0" / "model_manifest.json"

CLASS_LABELS = {
    0: "unlabeled_or_ambiguous",
    1: "abT_gold",
    2: "gdT_gold",
    3: "gdT_silver",
}
OVERLAP_ORDER = ["both", "latest_only", "v2_high_purity_only", "neither"]
UNION_ORDER = ["both", "latest_only", "v2_high_purity_only"]


@dataclass
class PredictionBundle:
    latest_score: np.ndarray
    latest_pred: np.ndarray
    v2_score: np.ndarray
    v2_threshold: np.ndarray
    v2_threshold_label: np.ndarray
    v2_pred: np.ndarray
    trdc_flag: np.ndarray
    trdv_flag: np.ndarray
    quadrant: np.ndarray


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare gdTAI latest local model with gdTAI v2 high-purity on the full atlas.")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--latest-model-pkl", type=Path, default=V3_MODEL)
    parser.add_argument("--v2-model-pkl", type=Path, default=V2_MODEL)
    parser.add_argument("--chunk-size", type=int, default=50_000)
    parser.add_argument("--score-scatter-sample", type=int, default=180_000)
    parser.add_argument("--top-source-figure-n", type=int, default=35)
    parser.add_argument("--no-pdf", action="store_true")
    return parser.parse_args()


def ensure_dirs() -> None:
    for path in [TABLE_DIR, FIGURE_DIR, LOG_DIR, STATIC_DIR, STATIC_ASSET_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def setup_logging() -> None:
    ensure_dirs()
    handlers = [logging.FileHandler(RUN_LOG, mode="w", encoding="utf-8"), logging.StreamHandler()]
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s", handlers=handlers, force=True)


def load_pickle(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        payload = pickle.load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected dict payload in {path}; found {type(payload)!r}")
    return payload


def h5ad_shape(handle: h5py.File) -> tuple[int, int]:
    if "X" in handle:
        x = handle["X"]
        if isinstance(x, h5py.Group) and "indptr" in x:
            n_obs = int(x["indptr"].shape[0] - 1)
            n_vars = int(handle["var"].get("_index", handle["var"].get("index")).shape[0])
            return n_obs, n_vars
        if hasattr(x, "shape"):
            return int(x.shape[0]), int(x.shape[1])
    raise RuntimeError("Could not determine H5AD shape from X.")


def read_var_names(handle: h5py.File) -> list[str]:
    var = handle["var"]
    key = "_index" if "_index" in var else "index"
    ds = var[key]
    try:
        values = ds.asstr()[:]
    except Exception:
        values = ds[:]
    return [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in values]


def build_feature_spec(handle: h5py.File, latest_payload: dict[str, Any], v2_payload: dict[str, Any]) -> FeatureSpec:
    var_names = read_var_names(handle)
    lookup = {gene: i for i, gene in enumerate(var_names)}
    latest_genes = [str(g) for g in latest_payload["gene_names"]]
    v2_genes = [str(g) for g in v2_payload["base_model"]["gene_names"]]
    selected: list[str] = []
    for gene in [*latest_genes, *v2_genes]:
        if gene not in selected:
            selected.append(gene)
    missing = [gene for gene in selected if gene not in lookup]
    if missing:
        raise KeyError(f"Model feature genes are missing from H5AD var: {missing[:20]} (n={len(missing)})")
    engineered = [str(x) for x in latest_payload.get("engineered_feature_names", ENGINEERED_FEATURES)]
    if engineered != list(ENGINEERED_FEATURES):
        logging.warning("Latest model engineered feature list differs from script constant; using payload order.")
    gene_feature_names = [f"{gene}_log1p_cp10k" for gene in selected]
    engineered_to_col = {name: len(selected) + i for i, name in enumerate(engineered)}
    return FeatureSpec(
        gene_names=selected,
        gene_indices=np.asarray([lookup[gene] for gene in selected], dtype=np.int32),
        gene_feature_names=gene_feature_names,
        engineered_feature_names=engineered,
        model_feature_names=[*gene_feature_names, *engineered],
        gene_to_col={gene: i for i, gene in enumerate(selected)},
        engineered_to_col=engineered_to_col,
    )


def column_indices_for_payload(spec: FeatureSpec, payload: dict[str, Any]) -> list[int]:
    gene_names = [str(g) for g in payload["gene_names"]]
    feature_names = [str(f) for f in payload.get("feature_names", [])]
    if not feature_names:
        return [spec.gene_to_col[g] for g in gene_names]
    cols: list[int] = []
    for feature in feature_names:
        if feature.endswith("_log1p_cp10k"):
            gene = feature[: -len("_log1p_cp10k")]
            cols.append(spec.gene_to_col[gene])
        else:
            cols.append(spec.engineered_to_col[feature])
    return cols


def annotation_vector(handle: h5py.File, n_obs: int) -> np.ndarray:
    if "simple_annotation_plus6" in handle["obs"]:
        return clean_group_values(read_obs_column(handle, "simple_annotation_plus6"))
    if "simple_annotation" in handle["obs"]:
        return clean_group_values(read_obs_column(handle, "simple_annotation"))
    return np.full(n_obs, "unknown", dtype=object)


def optional_float_obs(handle: h5py.File, key: str, n_obs: int) -> np.ndarray:
    if key in handle["obs"]:
        return read_float_obs(handle, key).astype(np.float32, copy=False)
    return np.full(n_obs, np.nan, dtype=np.float32)


def optional_bool_obs(handle: h5py.File, key: str, n_obs: int) -> np.ndarray:
    if key in handle["obs"]:
        return read_bool_obs(handle, key)
    return np.zeros(n_obs, dtype=bool)


def apply_models(handle: h5py.File, spec: FeatureSpec, latest_payload: dict[str, Any], v2_payload: dict[str, Any], annotation: np.ndarray, chunk_size: int) -> PredictionBundle:
    n_obs, _ = h5ad_shape(handle)
    latest_score = np.zeros(n_obs, dtype=np.float32)
    v2_score = np.zeros(n_obs, dtype=np.float32)
    v2_threshold = np.zeros(n_obs, dtype=np.float32)
    v2_threshold_label = np.empty(n_obs, dtype=object)
    latest_pred = np.zeros(n_obs, dtype=bool)
    v2_pred = np.zeros(n_obs, dtype=bool)
    trdc_flag = np.zeros(n_obs, dtype=bool)
    trdv_flag = np.zeros(n_obs, dtype=bool)

    latest_cols = column_indices_for_payload(spec, latest_payload)
    v2_cols = [spec.gene_to_col[str(gene)] for gene in v2_payload["base_model"]["gene_names"]]
    latest_model = latest_payload["model_object"]
    latest_threshold = float(latest_payload["threshold"])
    v2_model = v2_payload["base_model"]["model_object"]

    logging.info("Applying latest model threshold=%s and v2 high-purity thresholds to %s cells", latest_threshold, f"{n_obs:,}")
    for start in range(0, n_obs, chunk_size):
        end = min(start + chunk_size, n_obs)
        rows = np.arange(start, end, dtype=np.int64)
        x_gene, _row_sum, _n_detected = extract_gene_features(handle, "X", rows, spec, label=f"compare_{start}_{end}")
        x_all = append_engineered_features(x_gene, spec)
        trdc_flag[start:end] = x_all[:, spec.engineered_to_col["TRDC_log1p"]] > 0
        trdv_flag[start:end] = x_all[:, spec.engineered_to_col["any_TRDV"]] > 0.5
        local_latest_score = latest_model.predict_proba(x_all[:, latest_cols])[:, 1].astype(np.float32)
        local_v2_score = v2_model.predict_proba(x_all[:, v2_cols])[:, 1].astype(np.float32)
        threshold_vec, threshold_label = v2_threshold_vector(v2_payload, "high_purity", annotation[start:end])
        latest_score[start:end] = local_latest_score
        v2_score[start:end] = local_v2_score
        v2_threshold[start:end] = threshold_vec.astype(np.float32, copy=False)
        v2_threshold_label[start:end] = threshold_label
        latest_pred[start:end] = local_latest_score >= latest_threshold
        v2_pred[start:end] = local_v2_score >= threshold_vec
        if end % 500_000 == 0 or end == n_obs:
            logging.info("Compared models on %s / %s cells", f"{end:,}", f"{n_obs:,}")
    quadrant = quadrant_from_flags(trdc_flag, trdv_flag)
    return PredictionBundle(latest_score, latest_pred, v2_score, v2_threshold, v2_threshold_label, v2_pred, trdc_flag, trdv_flag, quadrant)


def overlap_category(latest_pred: np.ndarray, v2_pred: np.ndarray) -> np.ndarray:
    out = np.full(latest_pred.size, "neither", dtype=object)
    out[latest_pred & v2_pred] = "both"
    out[latest_pred & (~v2_pred)] = "latest_only"
    out[(~latest_pred) & v2_pred] = "v2_high_purity_only"
    return out


def class_flags(class_code: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "gdT_gold": class_code == 2,
        "abT_gold": class_code == 1,
        "gdT_silver": class_code == 3,
        "unlabeled_or_ambiguous": class_code == 0,
    }


def summarize_global(category: np.ndarray, bundle: PredictionBundle, class_code: np.ndarray, is_nk: np.ndarray, is_tcrab: np.ndarray, is_cd4_treg: np.ndarray, is_b_myeloid: np.ndarray) -> pd.DataFrame:
    n = category.size
    latest = bundle.latest_pred
    v2 = bundle.v2_pred
    both = latest & v2
    union = latest | v2
    rows = [
        {
            "comparison": "latest_vs_v2_high_purity",
            "total_cells": int(n),
            "latest_predicted": int(latest.sum()),
            "v2_high_purity_predicted": int(v2.sum()),
            "intersection_both": int(both.sum()),
            "union_predicted_by_either": int(union.sum()),
            "latest_only": int((latest & (~v2)).sum()),
            "v2_high_purity_only": int(((~latest) & v2).sum()),
            "neither": int((~union).sum()),
            "jaccard_index": safe_div(int(both.sum()), int(union.sum())),
            "fraction_of_latest_also_v2": safe_div(int(both.sum()), int(latest.sum())),
            "fraction_of_v2_also_latest": safe_div(int(both.sum()), int(v2.sum())),
            "latest_minus_v2_predicted": int(latest.sum()) - int(v2.sum()),
            "latest_only_fraction_of_latest": safe_div(int((latest & (~v2)).sum()), int(latest.sum())),
            "v2_only_fraction_of_v2": safe_div(int(((~latest) & v2).sum()), int(v2.sum())),
            "latest_predicted_NK": int((latest & is_nk).sum()),
            "v2_predicted_NK": int((v2 & is_nk).sum()),
            "both_predicted_NK": int((both & is_nk).sum()),
            "latest_only_NK": int((latest & (~v2) & is_nk).sum()),
            "v2_only_NK": int(((~latest) & v2 & is_nk).sum()),
            "latest_predicted_paired_or_any_TCRAB": int((latest & is_tcrab).sum()),
            "v2_predicted_paired_or_any_TCRAB": int((v2 & is_tcrab).sum()),
            "both_paired_or_any_TCRAB": int((both & is_tcrab).sum()),
            "latest_only_paired_or_any_TCRAB": int((latest & (~v2) & is_tcrab).sum()),
            "v2_only_paired_or_any_TCRAB": int(((~latest) & v2 & is_tcrab).sum()),
            "latest_predicted_CD4_Treg_warning": int((latest & is_cd4_treg).sum()),
            "v2_predicted_CD4_Treg_warning": int((v2 & is_cd4_treg).sum()),
            "latest_predicted_B_myeloid_warning": int((latest & is_b_myeloid).sum()),
            "v2_predicted_B_myeloid_warning": int((v2 & is_b_myeloid).sum()),
        }
    ]
    for label, mask in class_flags(class_code).items():
        rows[0][f"latest_{label}"] = int((latest & mask).sum())
        rows[0][f"v2_{label}"] = int((v2 & mask).sum())
        rows[0][f"both_{label}"] = int((both & mask).sum())
        rows[0][f"latest_only_{label}"] = int((latest & (~v2) & mask).sum())
        rows[0][f"v2_only_{label}"] = int(((~latest) & v2 & mask).sum())
    return pd.DataFrame(rows)


def aggregate_by_group(group_name: str, values: np.ndarray, category: np.ndarray, bundle: PredictionBundle, class_code: np.ndarray, is_nk: np.ndarray, is_tcrab: np.ndarray, is_cd4_treg: np.ndarray, is_b_myeloid: np.ndarray, quadrant: np.ndarray) -> pd.DataFrame:
    base = pd.DataFrame(
        {
            group_name: values,
            "latest_predicted": bundle.latest_pred.astype(np.int8),
            "v2_high_purity_predicted": bundle.v2_pred.astype(np.int8),
            "both": ((bundle.latest_pred) & (bundle.v2_pred)).astype(np.int8),
            "latest_only": ((bundle.latest_pred) & (~bundle.v2_pred)).astype(np.int8),
            "v2_high_purity_only": ((~bundle.latest_pred) & (bundle.v2_pred)).astype(np.int8),
            "neither": ((~bundle.latest_pred) & (~bundle.v2_pred)).astype(np.int8),
            "union_predicted_by_either": ((bundle.latest_pred) | (bundle.v2_pred)).astype(np.int8),
            "gdT_gold": (class_code == 2).astype(np.int8),
            "abT_gold": (class_code == 1).astype(np.int8),
            "gdT_silver": (class_code == 3).astype(np.int8),
            "unlabeled_or_ambiguous": (class_code == 0).astype(np.int8),
            "latest_gdT_gold": (bundle.latest_pred & (class_code == 2)).astype(np.int8),
            "v2_gdT_gold": (bundle.v2_pred & (class_code == 2)).astype(np.int8),
            "both_gdT_gold": (bundle.latest_pred & bundle.v2_pred & (class_code == 2)).astype(np.int8),
            "latest_only_gdT_gold": (bundle.latest_pred & (~bundle.v2_pred) & (class_code == 2)).astype(np.int8),
            "v2_only_gdT_gold": ((~bundle.latest_pred) & bundle.v2_pred & (class_code == 2)).astype(np.int8),
            "latest_gdT_silver": (bundle.latest_pred & (class_code == 3)).astype(np.int8),
            "v2_gdT_silver": (bundle.v2_pred & (class_code == 3)).astype(np.int8),
            "latest_abT_gold": (bundle.latest_pred & (class_code == 1)).astype(np.int8),
            "v2_abT_gold": (bundle.v2_pred & (class_code == 1)).astype(np.int8),
            "latest_NK": (bundle.latest_pred & is_nk).astype(np.int8),
            "v2_NK": (bundle.v2_pred & is_nk).astype(np.int8),
            "both_NK": (bundle.latest_pred & bundle.v2_pred & is_nk).astype(np.int8),
            "latest_only_NK": (bundle.latest_pred & (~bundle.v2_pred) & is_nk).astype(np.int8),
            "v2_only_NK": ((~bundle.latest_pred) & bundle.v2_pred & is_nk).astype(np.int8),
            "latest_TCRAB": (bundle.latest_pred & is_tcrab).astype(np.int8),
            "v2_TCRAB": (bundle.v2_pred & is_tcrab).astype(np.int8),
            "latest_only_TCRAB": (bundle.latest_pred & (~bundle.v2_pred) & is_tcrab).astype(np.int8),
            "v2_only_TCRAB": ((~bundle.latest_pred) & bundle.v2_pred & is_tcrab).astype(np.int8),
            "latest_CD4_Treg_warning": (bundle.latest_pred & is_cd4_treg).astype(np.int8),
            "v2_CD4_Treg_warning": (bundle.v2_pred & is_cd4_treg).astype(np.int8),
            "latest_B_myeloid_warning": (bundle.latest_pred & is_b_myeloid).astype(np.int8),
            "v2_B_myeloid_warning": (bundle.v2_pred & is_b_myeloid).astype(np.int8),
            "latest_TRDC_plus_TRDV_minus": (bundle.latest_pred & (quadrant == "TRDC+TRDV-")).astype(np.int8),
            "v2_TRDC_plus_TRDV_minus": (bundle.v2_pred & (quadrant == "TRDC+TRDV-")).astype(np.int8),
            "latest_only_TRDC_plus_TRDV_minus": (bundle.latest_pred & (~bundle.v2_pred) & (quadrant == "TRDC+TRDV-")).astype(np.int8),
            "v2_only_TRDC_plus_TRDV_minus": ((~bundle.latest_pred) & bundle.v2_pred & (quadrant == "TRDC+TRDV-")).astype(np.int8),
        }
    )
    grouped = base.groupby(group_name, dropna=False, sort=True).sum(numeric_only=True).reset_index()
    counts = base.groupby(group_name, dropna=False, sort=True).size().rename("total_cells").reset_index()
    out = counts.merge(grouped, on=group_name, how="left")
    out["latest_predicted_fraction"] = out["latest_predicted"] / out["total_cells"].replace(0, np.nan)
    out["v2_high_purity_predicted_fraction"] = out["v2_high_purity_predicted"] / out["total_cells"].replace(0, np.nan)
    out["latest_minus_v2_predicted"] = out["latest_predicted"] - out["v2_high_purity_predicted"]
    out["jaccard_index"] = out["both"] / out["union_predicted_by_either"].replace(0, np.nan)
    out["fraction_of_latest_also_v2"] = out["both"] / out["latest_predicted"].replace(0, np.nan)
    out["fraction_of_v2_also_latest"] = out["both"] / out["v2_high_purity_predicted"].replace(0, np.nan)
    out["latest_only_fraction_of_latest"] = out["latest_only"] / out["latest_predicted"].replace(0, np.nan)
    out["v2_only_fraction_of_v2"] = out["v2_high_purity_only"] / out["v2_high_purity_predicted"].replace(0, np.nan)
    out["latest_NK_fraction_of_latest"] = out["latest_NK"] / out["latest_predicted"].replace(0, np.nan)
    out["v2_NK_fraction_of_v2"] = out["v2_NK"] / out["v2_high_purity_predicted"].replace(0, np.nan)
    out["latest_TCRAB_fraction_of_latest"] = out["latest_TCRAB"] / out["latest_predicted"].replace(0, np.nan)
    out["v2_TCRAB_fraction_of_v2"] = out["v2_TCRAB"] / out["v2_high_purity_predicted"].replace(0, np.nan)
    out["latest_TRDC_plus_TRDV_minus_fraction_of_latest"] = out["latest_TRDC_plus_TRDV_minus"] / out["latest_predicted"].replace(0, np.nan)
    out["v2_TRDC_plus_TRDV_minus_fraction_of_v2"] = out["v2_TRDC_plus_TRDV_minus"] / out["v2_high_purity_predicted"].replace(0, np.nan)
    return out.replace([np.inf, -np.inf], np.nan)


def aggregate_overlap_by_class(category: np.ndarray, class_code: np.ndarray, is_nk: np.ndarray, is_tcrab: np.ndarray, quadrant: np.ndarray) -> pd.DataFrame:
    rows = []
    for cat in UNION_ORDER:
        mask = category == cat
        row = {"overlap_category": cat, "n_cells": int(mask.sum())}
        for code, label in CLASS_LABELS.items():
            row[label] = int((mask & (class_code == code)).sum())
            row[f"{label}_fraction"] = safe_div(row[label], row["n_cells"])
        row["NK_annotation"] = int((mask & is_nk).sum())
        row["NK_annotation_fraction"] = safe_div(row["NK_annotation"], row["n_cells"])
        row["paired_or_any_TCRAB"] = int((mask & is_tcrab).sum())
        row["paired_or_any_TCRAB_fraction"] = safe_div(row["paired_or_any_TCRAB"], row["n_cells"])
        for q in ["TRDC+TRDV+", "TRDC+TRDV-", "TRDC-TRDV+", "TRDC-TRDV-"]:
            row[q] = int((mask & (quadrant == q)).sum())
            row[f"{q}_fraction"] = safe_div(row[q], row["n_cells"])
        rows.append(row)
    return pd.DataFrame(rows)


def write_union_cells(
    path: Path,
    handle: h5py.File,
    obs: Any,
    source: np.ndarray,
    tissue: np.ndarray,
    annotation: np.ndarray,
    category: np.ndarray,
    bundle: PredictionBundle,
    is_nk: np.ndarray,
    is_tcrab: np.ndarray,
    is_cd4_treg: np.ndarray,
    is_b_myeloid: np.ndarray,
    phase4_trd: np.ndarray,
    phase4_trab: np.ndarray,
    phase4_delta: np.ndarray,
) -> pd.DataFrame:
    union = bundle.latest_pred | bundle.v2_pred
    idx = np.flatnonzero(union).astype(np.int64)
    cell_ids = obs_index_values(handle)[idx]
    df = pd.DataFrame(
        {
            "obs_index": idx,
            "cell_id": cell_ids,
            "source_gse_id": source[idx],
            "tissue": tissue[idx],
            "annotation": annotation[idx],
            "ground_truth_class_code": obs.class_code[idx].astype(int),
            "ground_truth_class": [CLASS_LABELS.get(int(x), "unknown") for x in obs.class_code[idx]],
            "overlap_category": category[idx],
            "latest_pred": bundle.latest_pred[idx].astype(int),
            "latest_score": bundle.latest_score[idx],
            "latest_threshold": 0.5,
            "v2_high_purity_pred": bundle.v2_pred[idx].astype(int),
            "v2_high_purity_score": bundle.v2_score[idx],
            "v2_high_purity_threshold": bundle.v2_threshold[idx],
            "v2_high_purity_threshold_label": bundle.v2_threshold_label[idx],
            "tcr_gene_quadrant": bundle.quadrant[idx],
            "is_NK_annotation": is_nk[idx],
            "is_paired_or_any_TCRAB_non_gdT_gold": is_tcrab[idx],
            "is_CD4_Treg_warning": is_cd4_treg[idx],
            "is_B_myeloid_warning": is_b_myeloid[idx],
            "has_TRA_TRB_paired": obs.has_TRA_TRB_paired[idx],
            "has_any_ab_tcr": obs.has_any_ab_tcr[idx],
            "corrected_has_any_gd_tcr": obs.corrected_has_any_gd_tcr[idx],
            "phase4_trd_score": phase4_trd[idx],
            "phase4_trab_score": phase4_trab[idx],
            "phase4_trd_minus_trab": phase4_delta[idx],
        }
    )
    df.to_csv(path, index=False, compression="gzip")
    return df


def copy_figures_to_static() -> None:
    for png in FIGURE_DIR.glob("*.png"):
        target = STATIC_ASSET_DIR / png.name
        target.write_bytes(png.read_bytes())


def plot_overlap_bar(overall: pd.DataFrame) -> Path:
    row = overall.iloc[0]
    labels = ["Both", "Latest only", "v2 high-purity only"]
    values = [int(row["intersection_both"]), int(row["latest_only"]), int(row["v2_high_purity_only"])]
    colors = ["#3A6EA5", "#2A9D8F", "#E76F51"]
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    bars = ax.bar(labels, values, color=colors)
    ax.set_ylabel("Cells")
    ax.set_title("Prediction Overlap: Latest gdTAI vs v2 High-Purity")
    ax.tick_params(axis="x", rotation=10)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{val:,}", ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    out = FIGURE_DIR / "prediction_overlap_counts.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def plot_source_delta(source_df: pd.DataFrame, top_n: int) -> Path:
    df = source_df.copy()
    df["difference_abs"] = df["latest_minus_v2_predicted"].abs()
    df = df.sort_values("difference_abs", ascending=False).head(top_n).sort_values("latest_minus_v2_predicted")
    colors = np.where(df["latest_minus_v2_predicted"] >= 0, "#2A9D8F", "#E76F51")
    fig_h = max(6.0, 0.22 * len(df) + 1.5)
    fig, ax = plt.subplots(figsize=(9.2, fig_h))
    ax.barh(df["source_gse_id"].astype(str), df["latest_minus_v2_predicted"], color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Latest predicted cells minus v2 high-purity predicted cells")
    ax.set_ylabel("Dataset")
    ax.set_title("Largest Per-Dataset Prediction Count Differences")
    fig.tight_layout()
    out = FIGURE_DIR / "per_dataset_prediction_delta_top.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def plot_source_jaccard(source_df: pd.DataFrame, top_n: int) -> Path:
    df = source_df.loc[source_df["union_predicted_by_either"] > 0].copy()
    df = df.sort_values("union_predicted_by_either", ascending=False).head(top_n).sort_values("jaccard_index")
    fig_h = max(6.0, 0.22 * len(df) + 1.5)
    fig, ax = plt.subplots(figsize=(9.2, fig_h))
    ax.barh(df["source_gse_id"].astype(str), df["jaccard_index"], color="#3A6EA5")
    ax.set_xlim(0, 1)
    ax.set_xlabel("Jaccard index among cells predicted by either model")
    ax.set_ylabel("Dataset")
    ax.set_title("Per-Dataset Agreement for Datasets with Most Predicted Cells")
    fig.tight_layout()
    out = FIGURE_DIR / "per_dataset_jaccard_top_predicted.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def plot_score_scatter(union_df: pd.DataFrame, sample_n: int) -> Path:
    if union_df.empty:
        raise RuntimeError("No union predictions available for scatter plot.")
    rng = np.random.default_rng(7)
    if union_df.shape[0] > sample_n:
        df = union_df.iloc[rng.choice(union_df.shape[0], size=sample_n, replace=False)].copy()
    else:
        df = union_df.copy()
    palette = {"both": "#3A6EA5", "latest_only": "#2A9D8F", "v2_high_purity_only": "#E76F51"}
    fig, ax = plt.subplots(figsize=(6.6, 5.8))
    for cat in UNION_ORDER:
        sub = df.loc[df["overlap_category"] == cat]
        if sub.empty:
            continue
        ax.scatter(sub["v2_high_purity_score"], sub["latest_score"], s=3, alpha=0.25, linewidths=0, c=palette[cat], label=f"{cat} ({sub.shape[0]:,} sampled)")
    ax.axhline(0.5, color="black", linestyle="--", linewidth=0.8, label="latest threshold 0.5")
    ax.set_xlabel("v2 base-model probability")
    ax.set_ylabel("latest gdTAI probability")
    ax.set_title("Score Relationship Among Cells Predicted by Either Model")
    ax.legend(markerscale=3, frameon=False, fontsize=8)
    fig.tight_layout()
    out = FIGURE_DIR / "score_scatter_union_predictions.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


def top_table(df: pd.DataFrame, sort_col: str, n: int, columns: list[str]) -> pd.DataFrame:
    return df.sort_values(sort_col, ascending=False).head(n)[columns]


def render_report(
    args: argparse.Namespace,
    latest_payload: dict[str, Any],
    v2_payload: dict[str, Any],
    overall: pd.DataFrame,
    source_df: pd.DataFrame,
    tissue_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    category_df: pd.DataFrame,
    figures: list[Path],
) -> None:
    model_name = str(latest_payload.get("model", latest_payload.get("notes", "latest local gdTAI")))
    latest_version = str(latest_payload.get("version", "unknown"))
    v2_version = str(v2_payload.get("version", "unknown"))
    source_cols = [
        "source_gse_id",
        "total_cells",
        "latest_predicted",
        "v2_high_purity_predicted",
        "both",
        "latest_only",
        "v2_high_purity_only",
        "latest_minus_v2_predicted",
        "jaccard_index",
        "latest_NK",
        "v2_NK",
        "latest_TCRAB",
        "v2_TCRAB",
        "latest_TRDC_plus_TRDV_minus",
        "v2_TRDC_plus_TRDV_minus",
    ]
    top_delta = top_table(source_df.assign(abs_delta=source_df["latest_minus_v2_predicted"].abs()), "abs_delta", 25, source_cols)
    top_latest_only = top_table(source_df, "latest_only", 25, source_cols)
    top_v2_only = top_table(source_df, "v2_high_purity_only", 25, source_cols)
    top_union = top_table(source_df, "union_predicted_by_either", 25, source_cols)
    tissue_cols = ["tissue", "total_cells", "latest_predicted", "v2_high_purity_predicted", "both", "latest_only", "v2_high_purity_only", "latest_minus_v2_predicted", "jaccard_index"]
    annotation_cols = ["annotation", "total_cells", "latest_predicted", "v2_high_purity_predicted", "both", "latest_only", "v2_high_purity_only", "latest_minus_v2_predicted", "jaccard_index"]
    md_parts = [
        "# gdTAI Latest vs v2 High-Purity Prediction Comparison",
        "",
        "## Models Compared",
        f"- Latest local model: `{model_name}` (`{latest_version}`), threshold `{float(latest_payload['threshold']):.6g}`.",
        f"- v2 comparator: `gdTAI_v{v2_version}` high-purity mode, annotation-specific thresholds from the v2 wrapper.",
        f"- Input H5AD: `{args.input_h5ad}`.",
        "- H5AD mutation: none; all predictions are computed from read-only HDF5 access.",
        "",
        "## Whole-Atlas Overlap",
        dataframe_to_markdown(overall),
        "",
        "## Overlap Category Composition",
        dataframe_to_markdown(category_df),
        "",
        "## Per-Dataset Differences: Largest Absolute Count Changes",
        dataframe_to_markdown(top_delta),
        "",
        "## Per-Dataset Differences: Latest-Only Enrichment",
        dataframe_to_markdown(top_latest_only),
        "",
        "## Per-Dataset Differences: v2-Only Enrichment",
        dataframe_to_markdown(top_v2_only),
        "",
        "## Per-Dataset Agreement: Largest Prediction Sets",
        dataframe_to_markdown(top_union),
        "",
        "## Per-Tissue Summary",
        dataframe_to_markdown(tissue_df.sort_values("union_predicted_by_either", ascending=False).head(30)[tissue_cols]),
        "",
        "## Per-Annotation Summary",
        dataframe_to_markdown(annotation_df.sort_values("union_predicted_by_either", ascending=False).head(30)[annotation_cols]),
        "",
        "## Output Files",
        f"- Tables: `{TABLE_DIR}`",
        f"- Figures: `{FIGURE_DIR}`",
        f"- Union predicted cell table: `{TABLE_DIR / 'union_predicted_cells_latest_vs_v2_high_purity.csv.gz'}`",
    ]
    REPORT_MD.write_text("\n".join(md_parts) + "\n", encoding="utf-8")

    fig_html = "\n".join(
        f'<figure><img src="assets/{OUT_PREFIX}/{html.escape(fig.name)}" alt="{html.escape(fig.stem)}"><figcaption>{html.escape(fig.stem.replace("_", " "))}</figcaption></figure>'
        for fig in figures
    )
    css = """
    body { font-family: Arial, sans-serif; margin: 28px auto; max-width: 1200px; color: #222; line-height: 1.45; }
    h1, h2 { color: #15324a; }
    table { border-collapse: collapse; width: 100%; font-size: 12px; margin: 12px 0 24px 0; }
    th, td { border: 1px solid #d8dee4; padding: 5px 7px; text-align: right; }
    th:first-child, td:first-child { text-align: left; }
    th { background: #eef3f7; }
    figure { margin: 22px 0; page-break-inside: avoid; }
    img { max-width: 100%; height: auto; border: 1px solid #ddd; }
    code { background: #f3f4f6; padding: 1px 4px; }
    .note { background: #f8fafc; border-left: 4px solid #3A6EA5; padding: 10px 12px; }
    """
    html_parts = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        "<title>gdTAI Latest vs v2 High-Purity</title>",
        f"<style>{css}</style>",
        "</head><body>",
        "<h1>gdTAI Latest vs v2 High-Purity Prediction Comparison</h1>",
        "<section class='note'>",
        f"<p><b>Latest local model:</b> {html.escape(model_name)} ({html.escape(latest_version)}), threshold {float(latest_payload['threshold']):.6g}.<br>",
        f"<b>Comparator:</b> gdTAI_v{html.escape(v2_version)} high-purity with annotation-specific thresholds.<br>",
        f"<b>Input:</b> {html.escape(str(args.input_h5ad))}. H5AD mutation: none.</p>",
        "</section>",
        "<h2>Figures</h2>",
        fig_html,
        "<h2>Whole-Atlas Overlap</h2>",
        dataframe_to_html(overall),
        "<h2>Overlap Category Composition</h2>",
        dataframe_to_html(category_df),
        "<h2>Per-Dataset Differences: Largest Absolute Count Changes</h2>",
        dataframe_to_html(top_delta),
        "<h2>Per-Dataset Differences: Latest-Only Enrichment</h2>",
        dataframe_to_html(top_latest_only),
        "<h2>Per-Dataset Differences: v2-Only Enrichment</h2>",
        dataframe_to_html(top_v2_only),
        "<h2>Per-Dataset Agreement: Largest Prediction Sets</h2>",
        dataframe_to_html(top_union),
        "<h2>Per-Tissue Summary</h2>",
        dataframe_to_html(tissue_df.sort_values("union_predicted_by_either", ascending=False).head(40)[tissue_cols]),
        "<h2>Per-Annotation Summary</h2>",
        dataframe_to_html(annotation_df.sort_values("union_predicted_by_either", ascending=False).head(40)[annotation_cols]),
        "<h2>Output Files</h2>",
        f"<p>Tables: <code>{html.escape(str(TABLE_DIR))}</code><br>Figures: <code>{html.escape(str(FIGURE_DIR))}</code></p>",
        "</body></html>",
    ]
    REPORT_HTML.write_text("\n".join(html_parts), encoding="utf-8")
    if not args.no_pdf:
        chrome = "/usr/bin/google-chrome"
        try:
            subprocess.run(
                [chrome, "--headless", "--disable-gpu", "--no-sandbox", f"--print-to-pdf={REPORT_PDF}", str(REPORT_HTML)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
        except Exception as exc:
            logging.warning("PDF export failed: %s", exc)


def main() -> None:
    args = parse_args()
    setup_logging()
    logging.info("Loading model payloads")
    latest_payload = load_pickle(args.latest_model_pkl)
    v2_payload = load_pickle(args.v2_model_pkl)
    manifest = {}
    if V3_MANIFEST.exists():
        manifest["latest_manifest"] = json.loads(V3_MANIFEST.read_text(encoding="utf-8"))
    if V2_MANIFEST.exists():
        manifest["v2_manifest"] = json.loads(V2_MANIFEST.read_text(encoding="utf-8"))

    with h5py.File(args.input_h5ad, "r") as handle:
        n_obs, n_vars = h5ad_shape(handle)
        logging.info("Input H5AD shape: %s cells x %s genes", f"{n_obs:,}", f"{n_vars:,}")
        obs = build_obs_metadata(handle)
        annotation = annotation_vector(handle, n_obs)
        source = obs.source
        tissue = obs.tissue
        annotation_upper = pd.Series(annotation, copy=False).astype(str).str.upper().to_numpy(dtype=object)
        is_nk = np.char.find(annotation_upper.astype(str), "NK") >= 0
        is_tcrab = (obs.has_TRA_TRB_paired | obs.has_any_ab_tcr) & (obs.class_code != 2)
        is_cd4_treg = (np.char.find(annotation_upper.astype(str), "CD4") >= 0) | (np.char.find(annotation_upper.astype(str), "TREG") >= 0)
        is_b_myeloid = (np.char.find(annotation_upper.astype(str), "B_CELL") >= 0) | (np.char.find(annotation_upper.astype(str), "MYELOID") >= 0)
        phase4_trd = optional_float_obs(handle, "phase4_trd_score", n_obs)
        phase4_trab = optional_float_obs(handle, "phase4_trab_score", n_obs)
        phase4_delta = optional_float_obs(handle, "phase4_trd_minus_trab", n_obs)
        spec = build_feature_spec(handle, latest_payload, v2_payload)
        pd.DataFrame(
            {
                "feature_index": np.arange(len(spec.model_feature_names), dtype=int),
                "feature": spec.model_feature_names,
                "feature_type": ["gene_log1p_cp10k"] * len(spec.gene_feature_names) + ["engineered"] * len(spec.engineered_feature_names),
                "gene": spec.gene_names + [""] * len(spec.engineered_feature_names),
            }
        ).to_csv(TABLE_DIR / "comparison_feature_manifest.csv", index=False)
        bundle = apply_models(handle, spec, latest_payload, v2_payload, annotation, args.chunk_size)
        category = overlap_category(bundle.latest_pred, bundle.v2_pred)
        overall = summarize_global(category, bundle, obs.class_code, is_nk, is_tcrab, is_cd4_treg, is_b_myeloid)
        source_df = aggregate_by_group("source_gse_id", source, category, bundle, obs.class_code, is_nk, is_tcrab, is_cd4_treg, is_b_myeloid, bundle.quadrant)
        tissue_df = aggregate_by_group("tissue", tissue, category, bundle, obs.class_code, is_nk, is_tcrab, is_cd4_treg, is_b_myeloid, bundle.quadrant)
        annotation_df = aggregate_by_group("annotation", annotation, category, bundle, obs.class_code, is_nk, is_tcrab, is_cd4_treg, is_b_myeloid, bundle.quadrant)
        category_df = aggregate_overlap_by_class(category, obs.class_code, is_nk, is_tcrab, bundle.quadrant)
        union_df = write_union_cells(
            TABLE_DIR / "union_predicted_cells_latest_vs_v2_high_purity.csv.gz",
            handle,
            obs,
            source,
            tissue,
            annotation,
            category,
            bundle,
            is_nk,
            is_tcrab,
            is_cd4_treg,
            is_b_myeloid,
            phase4_trd,
            phase4_trab,
            phase4_delta,
        )

    overall.to_csv(TABLE_DIR / "whole_atlas_overlap_summary.csv", index=False)
    source_df.to_csv(TABLE_DIR / "per_source_overlap_summary.csv", index=False)
    tissue_df.to_csv(TABLE_DIR / "per_tissue_overlap_summary.csv", index=False)
    annotation_df.to_csv(TABLE_DIR / "per_annotation_overlap_summary.csv", index=False)
    category_df.to_csv(TABLE_DIR / "overlap_category_composition.csv", index=False)

    figures = [
        plot_overlap_bar(overall),
        plot_source_delta(source_df, args.top_source_figure_n),
        plot_source_jaccard(source_df, args.top_source_figure_n),
        plot_score_scatter(union_df, args.score_scatter_sample),
    ]
    copy_figures_to_static()
    render_report(args, latest_payload, v2_payload, overall, source_df, tissue_df, annotation_df, category_df, figures)

    summary = {
        "input_h5ad": str(args.input_h5ad),
        "latest_model_pkl": str(args.latest_model_pkl),
        "v2_model_pkl": str(args.v2_model_pkl),
        "latest_model_name": str(latest_payload.get("model", "latest")),
        "latest_threshold": float(latest_payload["threshold"]),
        "whole_atlas_overlap": overall.iloc[0].to_dict(),
        "paths": {
            "tables": str(TABLE_DIR),
            "figures": str(FIGURE_DIR),
            "report_md": str(REPORT_MD),
            "report_html": str(REPORT_HTML),
            "report_pdf": str(REPORT_PDF) if REPORT_PDF.exists() else None,
        },
        **manifest,
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    logging.info("Done. Wrote report to %s", REPORT_HTML)


if __name__ == "__main__":
    main()
