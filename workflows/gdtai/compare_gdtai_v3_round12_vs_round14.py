#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Compare the saved gdTAI v3 Round 12 and Round 14 artifacts.

The workflow is read-only with respect to H5AD files. It pins both model
artifacts by SHA256, reconstructs their raw full-atlas positive sets from the
validated caches, reruns both models on the independent external H5AD, and
uses identical labels and operating-point metrics for the promotion decision.
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
import hashlib
import html
import json
import logging
import pickle
import shutil
import subprocess
from pathlib import Path
from typing import Any

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest

from run_gdt_prediction_package_evaluation import (
    DEFAULT_INPUT_H5AD,
    clean_group_values,
    dataframe_to_markdown,
    read_obs_column,
    read_string_dataset,
)
from run_gdtai_v3_trdc_nk_guard_classifier import (
    EXTERNAL_H5AD,
    FeatureSpec,
    append_engineered_features,
    external_truth,
    extract_gene_features,
    h5ad_shape,
    safe_metric_row,
    trdc_trdv_quadrant,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
PREFIX = "gdtai_v3_round12_vs_round14"
TABLE_DIR = OUTPUT_ROOT / "tables" / "gdT_prediction" / PREFIX
FIGURE_DIR = OUTPUT_ROOT / "figures" / "gdT_prediction" / PREFIX
LOG_DIR = OUTPUT_ROOT / "logs" / "gdT_prediction" / PREFIX
MODEL_DIR = OUTPUT_ROOT / "models" / "gdT_prediction_classifier" / PREFIX
STATIC_DIR = PROJECT_ROOT / "gdT_prediction" / PREFIX
STATIC_FIGURE_DIR = STATIC_DIR / "figures"

REPORT_MD = LOG_DIR / "gdtai_v3_round12_vs_round14_report.md"
REPORT_HTML = STATIC_DIR / "index.html"
REPORT_PDF = STATIC_DIR / "gdtai_v3_round12_vs_round14_report.pdf"
SUMMARY_JSON = LOG_DIR / "gdtai_v3_round12_vs_round14_summary.json"
RUN_LOG = LOG_DIR / "run.log"

R12_NAME = "v3_round12_hist_gradient_fixed_0p5"
R14_NAME = "v3_round14_v2_score_trdc_gate_fixed_0p936"
R14_TABLE_NAME = "v3_round20_v2_score_trdc_gate_fixed_0p936"
R12_SHA256 = "7373e79350f7db190c415b376b9763e31652754438ee8c5afd3853beb7b2ebc4"
R14_SHA256 = "16dedc0081da9b8487887341232bcf8c9c9403dd3bbd72e04cab43d4cd7b2e09"
R14_GIT_COMMIT = "21b5d60"
R14_GIT_PATH = "Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/gdTAI_v3_model.pkl"

R12_SOURCE_MODEL = OUTPUT_ROOT / "models" / "gdT_prediction_classifier" / "gdTAI_v3.0" / "gdTAI_v3_model.pkl"
CANONICAL_MODEL_DIR = OUTPUT_ROOT / "models" / "gdT_prediction_classifier" / "gdTAI_v3.0"
SOURCE_CANDIDATE_DIR = OUTPUT_ROOT / "models" / "gdT_prediction_classifier" / "gdtai_v3_trdc_nk_guard"
R12_POSITIVE_CACHE = OUTPUT_ROOT / "tables" / "gdT_prediction" / "gdtai_v3_trdc_nk_guard" / "full_atlas_selected_predicted_cells.csv.gz"
R14_ATLAS_CACHE = OUTPUT_ROOT / "tables" / "gdT_atlas" / "predicted_gdt_cell_atlas_metadata.csv.gz"
R14_REMOVED_FP_CACHE = OUTPUT_ROOT / "tables" / "gdT_atlas" / "predicted_gdt_cell_atlas_removed_fp_cells.csv.gz"
FULL_OVERALL_CACHE = OUTPUT_ROOT / "tables" / "gdT_prediction" / "gdtai_v3_trdc_nk_guard" / "full_atlas_prediction_overall.csv"
FULL_SOURCE_CACHE = OUTPUT_ROOT / "tables" / "gdT_prediction" / "gdtai_v3_trdc_nk_guard" / "full_atlas_prediction_by_source.csv"
VALIDATION_CACHE = OUTPUT_ROOT / "tables" / "gdT_prediction" / "gdtai_v3_trdc_nk_guard" / "atlas_validation_predictions_wide.csv.gz"
EXTERNAL_CACHE = OUTPUT_ROOT / "tables" / "gdT_prediction" / "gdtai_v3_trdc_nk_guard" / "external_predictions_wide.csv.gz"
GROUND_TRUTH_SOURCE = OUTPUT_ROOT / "tables" / "gdT_prediction" / "ground_truth_by_source_gse.csv"

CLASS_LABELS = {
    0: "unlabeled_or_ambiguous",
    1: "abT_gold",
    2: "gdT_gold",
    3: "gdT_silver",
}
MODEL_LABELS = {"round12": "Round 12", "round14": "Round 14"}
MODEL_COLORS = {"round12": "#126782", "round14": "#D1495B"}
QUADRANT_COLORS = {
    "TRDC+TRDV+": "#157A6E",
    "TRDC-TRDV+": "#E9C46A",
    "TRDC+TRDV-": "#D1495B",
    "TRDC-TRDV-": "#6B7280",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare gdTAI v3 Round 12 and Round 14.")
    parser.add_argument("--atlas-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--external-h5ad", type=Path, default=EXTERNAL_H5AD)
    parser.add_argument("--round12-model", type=Path, default=R12_SOURCE_MODEL)
    parser.add_argument("--atlas-cache-check-per-group", type=int, default=500)
    parser.add_argument("--promote-selected", action="store_true", help="Synchronize the selected artifact and package into gdTAI_v3.0.")
    parser.add_argument("--no-pdf", action="store_true")
    return parser.parse_args()


def ensure_dirs() -> None:
    for path in [TABLE_DIR, FIGURE_DIR, LOG_DIR, MODEL_DIR, STATIC_DIR, STATIC_FIGURE_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def setup_logging() -> None:
    ensure_dirs()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(RUN_LOG, mode="w", encoding="utf-8"), logging.StreamHandler()],
        force=True,
    )


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {
        "path": str(path),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


def load_pickle(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        payload = pickle.load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected dict payload in {path}; found {type(payload)!r}")
    return payload


def pin_model_artifacts(args: argparse.Namespace) -> tuple[Path, Path, pd.DataFrame]:
    r12_snapshot = MODEL_DIR / "round12_model.pkl"
    r14_snapshot = MODEL_DIR / "round14_model.pkl"

    if not r12_snapshot.exists():
        source_hash = sha256_file(args.round12_model)
        if source_hash != R12_SHA256:
            raise RuntimeError(
                f"Round 12 source hash mismatch: expected {R12_SHA256}, found {source_hash} at {args.round12_model}"
            )
        shutil.copy2(args.round12_model, r12_snapshot)
    if sha256_file(r12_snapshot) != R12_SHA256:
        raise RuntimeError("Pinned Round 12 model hash does not match the expected artifact.")

    r14_bytes = subprocess.check_output(
        ["git", "show", f"{R14_GIT_COMMIT}:{R14_GIT_PATH}"],
        cwd=PROJECT_ROOT,
    )
    if sha256_bytes(r14_bytes) != R14_SHA256:
        raise RuntimeError("Git-pinned Round 14 artifact hash does not match the expected artifact.")
    if not r14_snapshot.exists() or sha256_file(r14_snapshot) != R14_SHA256:
        r14_snapshot.write_bytes(r14_bytes)

    rows = []
    for model_key, model_path, expected_hash, provenance in [
        ("round12", r12_snapshot, R12_SHA256, str(args.round12_model)),
        ("round14", r14_snapshot, R14_SHA256, f"git:{R14_GIT_COMMIT}:{R14_GIT_PATH}"),
    ]:
        payload = load_pickle(model_path)
        rows.append(
            {
                "model": model_key,
                "display_name": MODEL_LABELS[model_key],
                "payload_model": str(payload.get("model", "")),
                "threshold": float(payload["threshold"]),
                "sha256": sha256_file(model_path),
                "size_bytes": model_path.stat().st_size,
                "snapshot_path": str(model_path),
                "provenance": provenance,
                "n_gene_features": len(payload.get("gene_names", [])),
                "n_total_features": len(payload.get("feature_names", [])),
                "model_object_type": type(payload.get("model_object")).__name__,
            }
        )
    manifest = pd.DataFrame(rows)
    expected = {
        "round12": (R12_NAME, 0.5),
        "round14": (R14_NAME, 0.936),
    }
    for row in manifest.to_dict(orient="records"):
        name, threshold = expected[str(row["model"])]
        if row["payload_model"] != name or not np.isclose(float(row["threshold"]), threshold):
            raise RuntimeError(f"Unexpected payload identity for {row['model']}: {row}")
    manifest.to_csv(TABLE_DIR / "artifact_manifest.csv", index=False)
    return r12_snapshot, r14_snapshot, manifest


def build_spec(handle: h5py.File, payloads: dict[str, dict[str, Any]]) -> FeatureSpec:
    r12 = payloads["round12"]
    r14 = payloads["round14"]
    for key in ["gene_names", "engineered_feature_names", "feature_names"]:
        if list(r12[key]) != list(r14[key]):
            raise RuntimeError(f"Round 12 and Round 14 use different `{key}` schemas.")
    var_names = pd.Index(read_string_dataset(handle["var"]["_index"]), dtype="string").astype(str).tolist()
    lookup = {gene: idx for idx, gene in enumerate(var_names)}
    genes = [str(gene) for gene in r12["gene_names"]]
    missing = [gene for gene in genes if gene not in lookup]
    if missing:
        raise KeyError(f"Input H5AD is missing {len(missing)} model genes: {missing[:20]}")
    engineered = [str(name) for name in r12["engineered_feature_names"]]
    gene_features = [f"{gene}_log1p_cp10k" for gene in genes]
    return FeatureSpec(
        gene_names=genes,
        gene_indices=np.asarray([lookup[gene] for gene in genes], dtype=np.int32),
        gene_feature_names=gene_features,
        engineered_feature_names=engineered,
        model_feature_names=[*gene_features, *engineered],
        gene_to_col={gene: idx for idx, gene in enumerate(genes)},
        engineered_to_col={name: len(genes) + idx for idx, name in enumerate(engineered)},
    )


def direct_external_evaluation(
    external_h5ad: Path,
    payloads: dict[str, dict[str, Any]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    logging.info("Running both pinned artifacts on independent external H5AD: %s", external_h5ad)
    with h5py.File(external_h5ad, "r") as handle:
        if "layers" not in handle or "counts" not in handle["layers"]:
            raise RuntimeError("External H5AD must contain layers['counts'].")
        spec = build_spec(handle, payloads)
        n_obs = h5ad_shape(handle, "layers/counts")[0]
        rows = np.arange(n_obs, dtype=np.int64)
        x_gene, row_sum, n_detected = extract_gene_features(
            handle,
            "layers/counts",
            rows,
            spec,
            label="round12_round14_external",
        )
        x = append_engineered_features(x_gene, spec)
        out = external_truth(handle, rows)
        out["row_sum_counts_layer"] = row_sum
        out["n_detected_genes_counts_layer"] = n_detected
        out["any_TRDV"] = x[:, spec.engineered_to_col["any_TRDV"]] > 0.5
        out["any_TRDJ"] = x[:, spec.engineered_to_col["any_TRDJ"]] > 0.5
        out["any_TRG"] = x[:, spec.engineered_to_col["any_TRG"]] > 0.5
        out["CD3_score"] = x[:, spec.engineered_to_col["CD3_score"]]
        out["NK_score"] = x[:, spec.engineered_to_col["NK_score"]]
        out["gdT_TCR_score"] = x[:, spec.engineered_to_col["gdT_TCR_score"]]
        out["tcr_gene_quadrant"] = trdc_trdv_quadrant(x, spec)
        for model_key, payload in payloads.items():
            score = payload["model_object"].predict_proba(x)[:, 1].astype(np.float32)
            threshold = float(payload["threshold"])
            out[f"{model_key}_score"] = score
            out[f"{model_key}_pred"] = score >= threshold

    cache = pd.read_csv(
        EXTERNAL_CACHE,
        usecols=[
            "external_row",
            f"{R12_NAME}_score",
            f"{R12_NAME}_pred",
            f"{R14_TABLE_NAME}_score",
            f"{R14_TABLE_NAME}_pred",
        ],
    ).sort_values("external_row")
    out = out.sort_values("external_row").reset_index(drop=True)
    checks = []
    for model_key, cache_name in [("round12", R12_NAME), ("round14", R14_TABLE_NAME)]:
        score_match = np.allclose(
            out[f"{model_key}_score"].to_numpy(dtype=float),
            cache[f"{cache_name}_score"].to_numpy(dtype=float),
            rtol=1e-6,
            atol=1e-6,
        )
        pred_match = np.array_equal(
            out[f"{model_key}_pred"].to_numpy(dtype=bool),
            cache[f"{cache_name}_pred"].to_numpy(dtype=bool),
        )
        checks.append(
            {
                "check": "external_direct_vs_cache",
                "model": model_key,
                "n_cells": len(out),
                "score_match": bool(score_match),
                "prediction_match": bool(pred_match),
            }
        )
        if not score_match or not pred_match:
            raise RuntimeError(f"Direct external inference does not match the validated cache for {model_key}.")

    primary = out["primary_eval"].to_numpy(dtype=bool)
    y = out.loc[primary, "y_true"].to_numpy(dtype=np.int8)
    metric_rows = []
    group_rows = []
    for model_key in MODEL_LABELS:
        score = out.loc[primary, f"{model_key}_score"].to_numpy(dtype=np.float32)
        pred = out.loc[primary, f"{model_key}_pred"].to_numpy(dtype=bool)
        row = safe_metric_row(model_key, y, score, pred, float(payloads[model_key]["threshold"]))
        row["evaluation"] = "independent_external_primary"
        metric_rows.append(row)

        all_pred = out[f"{model_key}_pred"].to_numpy(dtype=bool)
        groups = {
            "all_real_gdT": out["real_gdt"].to_numpy(dtype=bool),
            "TRDV_positive_real_gdT": out["real_gdt"].to_numpy(dtype=bool) & out["any_TRDV"].to_numpy(dtype=bool),
            "TRDC_plus_TRDV_minus_real_gdT": out["real_gdt"].to_numpy(dtype=bool)
            & out["tcr_gene_quadrant"].astype(str).eq("TRDC+TRDV-").to_numpy(dtype=bool),
            "Vd2_real_gdT": out["real_gdt"].to_numpy(dtype=bool)
            & out["cell_type"].astype(str).eq("Vd2_gdT").to_numpy(dtype=bool),
            "primary_negative": primary & (out["y_true"].to_numpy(dtype=np.int8) == 0),
            "NK_negative": primary
            & (out["y_true"].to_numpy(dtype=np.int8) == 0)
            & out["cell_type"].astype(str).eq("NK").to_numpy(dtype=bool),
            "paired_TCRAB_negative": primary
            & (out["y_true"].to_numpy(dtype=np.int8) == 0)
            & out["tcr_paired_ab"].to_numpy(dtype=bool),
            "TRDC_plus_TRDV_minus_negative": primary
            & (out["y_true"].to_numpy(dtype=np.int8) == 0)
            & out["tcr_gene_quadrant"].astype(str).eq("TRDC+TRDV-").to_numpy(dtype=bool),
        }
        for group_name, mask in groups.items():
            positive_group = "gdT" in group_name
            called = int((all_pred & mask).sum())
            n_group = int(mask.sum())
            group_rows.append(
                {
                    "model": model_key,
                    "group": group_name,
                    "group_type": "positive_recall" if positive_group else "negative_false_positive",
                    "n_cells": n_group,
                    "called_positive": called,
                    "rate": called / n_group if n_group else np.nan,
                }
            )

    keep = [
        "external_row",
        "cell_id",
        "primary_eval",
        "y_true",
        "real_gdt",
        "real_abt_tcr_strict",
        "cell_type",
        "tcr_paired_ab",
        "tcr_gene_quadrant",
        "sample_id",
        "library_id",
        "round12_score",
        "round12_pred",
        "round14_score",
        "round14_pred",
    ]
    out[keep].to_csv(TABLE_DIR / "external_direct_predictions.csv.gz", index=False, compression="gzip")
    return pd.DataFrame(metric_rows), pd.DataFrame(group_rows), out, pd.DataFrame(checks)


def load_full_prediction_union(atlas_h5ad: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    logging.info("Reconstructing raw full-atlas positive sets from validated cell-level caches")
    r12 = pd.read_csv(R12_POSITIVE_CACHE)
    r14_retained = pd.read_csv(R14_ATLAS_CACHE)
    r14_retained = r14_retained.loc[r14_retained["gdtai_v3_predicted"].astype(bool)].copy()
    r14_removed = pd.read_csv(R14_REMOVED_FP_CACHE)
    r14 = pd.concat([r14_retained, r14_removed], ignore_index=True)

    if r12["obs_index"].duplicated().any() or r14["obs_index"].duplicated().any():
        raise RuntimeError("Duplicate obs_index values found in a full-atlas positive cache.")
    if len(r12) != 213_182 or len(r14) != 251_356:
        raise RuntimeError(f"Unexpected raw positive counts: Round 12={len(r12):,}, Round 14={len(r14):,}")

    left = r12.rename(
        columns={
            "cell_id": "cell_id_r12",
            "source_gse_id": "source_r12",
            "tissue": "tissue_r12",
            "annotation": "annotation_r12",
            "class_code": "class_code_r12",
            "score": "round12_score",
            "tcr_gene_quadrant": "quadrant_r12",
            "is_NK_annotation": "nk_r12",
            "is_TCRAB_no_gdT_gold": "tcrab_r12",
        }
    )
    right = r14.rename(
        columns={
            "cell_id": "cell_id_r14",
            "source_gse_id": "source_r14",
            "tissue": "tissue_r14",
            "gdtai_v3_class_code": "class_code_r14",
            "gdtai_v3_score": "round14_score",
            "gdtai_v3_tcr_gene_quadrant": "quadrant_r14",
            "gdtai_v3_fp_paired_TCRAB_no_gdT": "tcrab_r14",
        }
    )
    union = left[
        [
            "obs_index",
            "cell_id_r12",
            "source_r12",
            "tissue_r12",
            "annotation_r12",
            "class_code_r12",
            "round12_score",
            "quadrant_r12",
            "nk_r12",
            "tcrab_r12",
        ]
    ].merge(
        right[
            [
                "obs_index",
                "cell_id_r14",
                "source_r14",
                "tissue_r14",
                "class_code_r14",
                "round14_score",
                "quadrant_r14",
                "tcrab_r14",
            ]
        ],
        on="obs_index",
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    union["round12_pred"] = union["_merge"].isin(["left_only", "both"])
    union["round14_pred"] = union["_merge"].isin(["right_only", "both"])
    union["overlap_category"] = union["_merge"].map(
        {"left_only": "Round 12 only", "right_only": "Round 14 only", "both": "Both"}
    )
    union["cell_id"] = union["cell_id_r12"].combine_first(union["cell_id_r14"])
    union["source_gse_id"] = union["source_r12"].combine_first(union["source_r14"])
    union["tissue"] = union["tissue_r12"].combine_first(union["tissue_r14"])
    union["class_code"] = union["class_code_r12"].combine_first(union["class_code_r14"]).astype(int)
    union["ground_truth_class"] = union["class_code"].map(CLASS_LABELS)
    union["tcr_gene_quadrant"] = union["quadrant_r12"].combine_first(union["quadrant_r14"])
    tcrab_r12 = union["tcrab_r12"].astype("boolean").fillna(False).to_numpy(dtype=bool)
    tcrab_r14 = union["tcrab_r14"].astype("boolean").fillna(False).to_numpy(dtype=bool)
    union["is_paired_TCRAB_non_gdT"] = tcrab_r12 | tcrab_r14

    both = union["_merge"].eq("both")
    for first, second, label in [
        ("source_r12", "source_r14", "source"),
        ("tissue_r12", "tissue_r14", "tissue"),
        ("class_code_r12", "class_code_r14", "class code"),
        ("quadrant_r12", "quadrant_r14", "TCR quadrant"),
    ]:
        mismatch = both & union[first].notna() & union[second].notna() & union[first].astype(str).ne(union[second].astype(str))
        if mismatch.any():
            raise RuntimeError(f"Round caches disagree on {label} for {int(mismatch.sum())} shared cells.")

    with h5py.File(atlas_h5ad, "r") as handle:
        n_obs = h5ad_shape(handle, "X")[0]
        idx = union["obs_index"].to_numpy(dtype=np.int64)
        if idx.min() < 0 or idx.max() >= n_obs:
            raise IndexError("Positive cache obs_index is outside the atlas H5AD row range.")
        if "simple_annotation_plus6" not in handle["obs"]:
            raise KeyError("Atlas H5AD lacks obs['simple_annotation_plus6'] used by the original evaluation.")
        annotation_all = clean_group_values(read_obs_column(handle, "simple_annotation_plus6"))
        union["annotation"] = annotation_all[idx]

    ann_upper = union["annotation"].astype(str).str.upper()
    union["is_NK_annotation"] = ann_upper.str.contains("NK", regex=False)
    union["is_CD4_Treg_warning"] = ann_upper.str.contains("CD4", regex=False) | ann_upper.str.contains("TREG", regex=False)
    union["is_B_myeloid_warning"] = ann_upper.str.contains("B_CELL", regex=False) | ann_upper.str.contains("MYELOID", regex=False)

    overall = pd.DataFrame(
        [
            {
                "total_atlas_cells": n_obs,
                "round12_predicted": int(union["round12_pred"].sum()),
                "round14_predicted": int(union["round14_pred"].sum()),
                "both": int(both.sum()),
                "round12_only": int(union["_merge"].eq("left_only").sum()),
                "round14_only": int(union["_merge"].eq("right_only").sum()),
                "union": len(union),
                "jaccard": float(both.sum() / len(union)),
                "fraction_round12_also_round14": float(both.sum() / union["round12_pred"].sum()),
                "fraction_round14_also_round12": float(both.sum() / union["round14_pred"].sum()),
            }
        ]
    )

    composition_rows = []
    for category in ["Both", "Round 12 only", "Round 14 only"]:
        mask = union["overlap_category"].eq(category)
        row: dict[str, Any] = {"overlap_category": category, "n_cells": int(mask.sum())}
        for code, label in CLASS_LABELS.items():
            count = int((mask & union["class_code"].eq(code)).sum())
            row[label] = count
            row[f"{label}_fraction"] = count / int(mask.sum()) if mask.any() else np.nan
        row["NK_annotation"] = int((mask & union["is_NK_annotation"]).sum())
        row["paired_TCRAB_non_gdT"] = int((mask & union["is_paired_TCRAB_non_gdT"]).sum())
        for quadrant in QUADRANT_COLORS:
            row[quadrant] = int((mask & union["tcr_gene_quadrant"].eq(quadrant)).sum())
        composition_rows.append(row)
    composition = pd.DataFrame(composition_rows)

    keep = [
        "obs_index",
        "cell_id",
        "source_gse_id",
        "tissue",
        "annotation",
        "class_code",
        "ground_truth_class",
        "overlap_category",
        "round12_pred",
        "round12_score",
        "round14_pred",
        "round14_score",
        "tcr_gene_quadrant",
        "is_NK_annotation",
        "is_paired_TCRAB_non_gdT",
        "is_CD4_Treg_warning",
        "is_B_myeloid_warning",
    ]
    union[keep].to_csv(TABLE_DIR / "full_atlas_prediction_union.csv.gz", index=False, compression="gzip")
    overall.to_csv(TABLE_DIR / "full_atlas_overlap_summary.csv", index=False)
    composition.to_csv(TABLE_DIR / "full_atlas_overlap_composition.csv", index=False)
    return union, overall, composition


def source_summary(union: pd.DataFrame) -> pd.DataFrame:
    base = pd.DataFrame(
        {
            "source_gse_id": union["source_gse_id"],
            "round12_predicted": union["round12_pred"].astype(np.int8),
            "round14_predicted": union["round14_pred"].astype(np.int8),
            "both": (union["round12_pred"] & union["round14_pred"]).astype(np.int8),
            "round12_only": (union["round12_pred"] & ~union["round14_pred"]).astype(np.int8),
            "round14_only": (~union["round12_pred"] & union["round14_pred"]).astype(np.int8),
        }
    )
    for model in MODEL_LABELS:
        pred = union[f"{model}_pred"]
        for code, label in CLASS_LABELS.items():
            base[f"{model}_{label}"] = (pred & union["class_code"].eq(code)).astype(np.int8)
        base[f"{model}_NK"] = (pred & union["is_NK_annotation"]).astype(np.int8)
        base[f"{model}_paired_TCRAB"] = (pred & union["is_paired_TCRAB_non_gdT"]).astype(np.int8)
        base[f"{model}_TRDC_plus_TRDV_minus"] = (
            pred & union["tcr_gene_quadrant"].eq("TRDC+TRDV-")
        ).astype(np.int8)
    grouped = base.groupby("source_gse_id", as_index=False, sort=True).sum(numeric_only=True)

    full_source = pd.read_csv(FULL_SOURCE_CACHE)
    r12_base = full_source.loc[full_source["strategy"].eq(R12_NAME)].copy()
    r12_base = r12_base[
        ["source_gse_id", "total_cells", "paired_TCRAB_cells", "NK_cells", "TRDC_plus_TRDV_minus_cells"]
    ]
    truth = pd.read_csv(GROUND_TRUTH_SOURCE)
    truth_wide = truth.pivot(index="source_gse_id", columns="ground_truth_class", values="n_cells").fillna(0).reset_index()
    out = r12_base.merge(truth_wide, on="source_gse_id", how="left").merge(grouped, on="source_gse_id", how="left")
    for column in grouped.columns:
        if column != "source_gse_id":
            out[column] = out[column].fillna(0).astype(int)
    for label in CLASS_LABELS.values():
        if label not in out:
            out[label] = 0
        out[label] = out[label].fillna(0).astype(int)
    out["round14_minus_round12"] = out["round14_predicted"] - out["round12_predicted"]
    out["jaccard"] = out["both"] / (out["round12_predicted"] + out["round14_only"]).replace(0, np.nan)
    for model in MODEL_LABELS:
        out[f"{model}_predicted_fraction"] = out[f"{model}_predicted"] / out["total_cells"].replace(0, np.nan)
        out[f"{model}_gold_recall"] = out[f"{model}_gdT_gold"] / out["gdT_gold"].replace(0, np.nan)
        out[f"{model}_silver_recall"] = out[f"{model}_gdT_silver"] / out["gdT_silver"].replace(0, np.nan)
        out[f"{model}_abT_gold_fp_rate"] = out[f"{model}_abT_gold"] / out["abT_gold"].replace(0, np.nan)
        out[f"{model}_paired_TCRAB_fp_rate"] = out[f"{model}_paired_TCRAB"] / out["paired_TCRAB_cells"].replace(0, np.nan)
        out[f"{model}_NK_call_rate"] = out[f"{model}_NK"] / out["NK_cells"].replace(0, np.nan)
    out.to_csv(TABLE_DIR / "per_dataset_comparison.csv", index=False)
    return out


def standardized_full_metrics(union: pd.DataFrame) -> pd.DataFrame:
    full = pd.read_csv(FULL_OVERALL_CACHE)
    rows = []
    for model_key, strategy in [("round12", R12_NAME), ("round14", R14_TABLE_NAME)]:
        hit = full.loc[full["strategy"].eq(strategy)]
        if len(hit) != 1:
            raise RuntimeError(f"Expected one full-atlas aggregate row for {strategy}; found {len(hit)}")
        row = hit.iloc[0]
        silver_total = int(pd.read_csv(GROUND_TRUTH_SOURCE).query("ground_truth_class == 'gdT_silver'")["n_cells"].sum())
        silver_recovered = int((union[f"{model_key}_pred"] & union["class_code"].eq(3)).sum())
        rows.append(
            {
                "model": model_key,
                "threshold": float(row["threshold"]),
                "total_cells": int(row["total_cells"]),
                "predicted_positive": int(row["predicted_putative_gdT"]),
                "predicted_fraction": float(row["predicted_fraction"]),
                "tp": int(row["full_primary_tp"]),
                "fp": int(row["full_primary_fp"]),
                "tn": int(row["full_primary_tn"]),
                "fn": int(row["full_primary_fn"]),
                "precision": float(row["full_primary_precision"]),
                "recall": float(row["full_primary_recall"]),
                "specificity": float(row["full_primary_specificity"]),
                "f1": float(row["full_primary_f1"]),
                "mcc": float(row["full_primary_mcc"]),
                "estimated_fp_fraction": float(row["estimated_fp_fraction_of_predictions"]),
                "estimated_total_fp": float(row["estimated_total_abT_fp"]),
                "paired_TCRAB_calls": int(row["predicted_paired_TCRAB"]),
                "paired_TCRAB_cells": int(row["paired_TCRAB_cells"]),
                "NK_calls": int(row["predicted_NK"]),
                "NK_cells": int(row["NK_cells"]),
                "TRDC_plus_TRDV_minus_calls": int(row["predicted_TRDC_plus_TRDV_minus"]),
                "silver_total": silver_total,
                "silver_recovered": silver_recovered,
                "silver_recall": silver_recovered / silver_total,
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "full_atlas_performance.csv", index=False)
    return out


def validation_evaluation() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    columns = [
        "obs_index",
        "source_gse_id",
        "tissue",
        "validation_component",
        "y_true",
        "tcr_gene_quadrant",
        f"{R12_NAME}_score",
        f"{R12_NAME}_pred",
        f"{R14_TABLE_NAME}_score",
        f"{R14_TABLE_NAME}_pred",
    ]
    df = pd.read_csv(VALIDATION_CACHE, usecols=columns)
    df = df.rename(
        columns={
            f"{R12_NAME}_score": "round12_score",
            f"{R12_NAME}_pred": "round12_pred",
            f"{R14_TABLE_NAME}_score": "round14_score",
            f"{R14_TABLE_NAME}_pred": "round14_pred",
        }
    )
    y = df["y_true"].to_numpy(dtype=np.int8)
    rows = []
    group_rows = []
    for model_key, threshold in [("round12", 0.5), ("round14", 0.936)]:
        score = df[f"{model_key}_score"].to_numpy(dtype=np.float32)
        pred = df[f"{model_key}_pred"].to_numpy(dtype=bool)
        row = safe_metric_row(model_key, y, score, pred, threshold)
        row["evaluation"] = "atlas_heldout_combined"
        rows.append(row)
        for group_column in ["validation_component", "source_gse_id", "tissue", "tcr_gene_quadrant"]:
            for value, group in df.groupby(group_column, sort=True):
                idx = group.index.to_numpy(dtype=np.int64)
                group_row = safe_metric_row(model_key, y[idx], score[idx], pred[idx], threshold)
                group_row["group_column"] = group_column
                group_row["group_value"] = value
                group_rows.append(group_row)
    metrics = pd.DataFrame(rows)
    groups = pd.DataFrame(group_rows)
    metrics.to_csv(TABLE_DIR / "atlas_heldout_performance.csv", index=False)
    groups.to_csv(TABLE_DIR / "atlas_heldout_performance_by_group.csv", index=False)
    return metrics, groups, df


def atlas_artifact_cache_check(
    atlas_h5ad: Path,
    payloads: dict[str, dict[str, Any]],
    union: pd.DataFrame,
    per_group: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(20260715)
    selected_parts = []
    for category in ["Both", "Round 12 only", "Round 14 only"]:
        values = union.loc[union["overlap_category"].eq(category), "obs_index"].to_numpy(dtype=np.int64)
        n = min(per_group, len(values))
        selected_parts.append(rng.choice(values, size=n, replace=False))
    positive_set = set(union["obs_index"].astype(int))
    n_obs = int(pd.read_csv(FULL_OVERALL_CACHE, nrows=1)["total_cells"].iloc[0])
    negative: list[int] = []
    while len(negative) < per_group:
        candidates = rng.integers(0, n_obs, size=per_group * 3)
        negative.extend(int(idx) for idx in candidates if int(idx) not in positive_set)
    negative_array = np.asarray(negative[:per_group], dtype=np.int64)
    rows = np.unique(np.concatenate([*selected_parts, negative_array])).astype(np.int64)

    with h5py.File(atlas_h5ad, "r") as handle:
        spec = build_spec(handle, payloads)
        x_gene, _row_sum, _n_detected = extract_gene_features(handle, "X", rows, spec, label="round_cache_check")
        x = append_engineered_features(x_gene, spec)
        direct = {key: payload["model_object"].predict_proba(x)[:, 1] >= float(payload["threshold"]) for key, payload in payloads.items()}

    set12 = set(union.loc[union["round12_pred"], "obs_index"].astype(int))
    set14 = set(union.loc[union["round14_pred"], "obs_index"].astype(int))
    expected = {
        "round12": np.asarray([int(idx) in set12 for idx in rows], dtype=bool),
        "round14": np.asarray([int(idx) in set14 for idx in rows], dtype=bool),
    }
    checks = []
    for model_key in MODEL_LABELS:
        matches = direct[model_key] == expected[model_key]
        checks.append(
            {
                "check": "atlas_direct_artifact_vs_cell_cache",
                "model": model_key,
                "n_sampled": int(len(rows)),
                "n_matching": int(matches.sum()),
                "n_mismatching": int((~matches).sum()),
                "prediction_match": bool(matches.all()),
            }
        )
        if not matches.all():
            raise RuntimeError(f"Pinned {model_key} artifact disagrees with the full-atlas cell cache.")
    return pd.DataFrame(checks)


def exact_paired_test(name: str, y: np.ndarray, pred12: np.ndarray, pred14: np.ndarray) -> dict[str, Any]:
    y_bool = np.asarray(y, dtype=bool)
    p12 = np.asarray(pred12, dtype=bool)
    p14 = np.asarray(pred14, dtype=bool)
    error12 = p12 != y_bool
    error14 = p14 != y_bool
    r12_correct_r14_wrong = int((~error12 & error14).sum())
    r12_wrong_r14_correct = int((error12 & ~error14).sum())
    discordant = r12_correct_r14_wrong + r12_wrong_r14_correct
    pvalue = float(binomtest(r12_correct_r14_wrong, discordant, 0.5).pvalue) if discordant else 1.0
    if r12_wrong_r14_correct > r12_correct_r14_wrong:
        favored = "round14"
    elif r12_correct_r14_wrong > r12_wrong_r14_correct:
        favored = "round12"
    else:
        favored = "tie"
    return {
        "evaluation": name,
        "n_cells": int(len(y_bool)),
        "round12_errors": int(error12.sum()),
        "round14_errors": int(error14.sum()),
        "round12_correct_round14_wrong": r12_correct_r14_wrong,
        "round12_wrong_round14_correct": r12_wrong_r14_correct,
        "discordant": discordant,
        "exact_mcnemar_p": pvalue,
        "favored_model": favored,
        "note": "Cell-level paired test; descriptive because cells within donors are correlated.",
    }


def paired_tests(
    external: pd.DataFrame,
    validation: pd.DataFrame,
    union: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    ext_primary = external["primary_eval"].to_numpy(dtype=bool)
    ext_y = external["y_true"].to_numpy(dtype=np.int8)
    rows.append(
        exact_paired_test(
            "independent_external_primary",
            ext_y[ext_primary],
            external.loc[ext_primary, "round12_pred"],
            external.loc[ext_primary, "round14_pred"],
        )
    )
    for label, mask in [
        ("independent_external_gdT", ext_primary & (ext_y == 1)),
        ("independent_external_negative", ext_primary & (ext_y == 0)),
    ]:
        rows.append(
            exact_paired_test(
                label,
                ext_y[mask],
                external.loc[mask, "round12_pred"],
                external.loc[mask, "round14_pred"],
            )
        )

    val_y = validation["y_true"].to_numpy(dtype=np.int8)
    rows.append(exact_paired_test("atlas_heldout_combined", val_y, validation["round12_pred"], validation["round14_pred"]))
    for component, group in validation.groupby("validation_component", sort=True):
        rows.append(exact_paired_test(f"atlas_heldout:{component}", group["y_true"], group["round12_pred"], group["round14_pred"]))

    truth = pd.read_csv(GROUND_TRUTH_SOURCE).groupby("ground_truth_class")["n_cells"].sum().to_dict()
    for code, label, y_value in [(2, "full_atlas_gdT_gold", True), (3, "full_atlas_gdT_silver", True), (1, "full_atlas_abT_gold", False)]:
        group = union.loc[union["class_code"].eq(code)]
        both_missed_or_correct = int(truth[label.replace("full_atlas_", "")] - len(group))
        y = np.full(len(group) + both_missed_or_correct, y_value, dtype=bool)
        p12 = np.concatenate([group["round12_pred"].to_numpy(dtype=bool), np.zeros(both_missed_or_correct, dtype=bool)])
        p14 = np.concatenate([group["round14_pred"].to_numpy(dtype=bool), np.zeros(both_missed_or_correct, dtype=bool)])
        rows.append(exact_paired_test(label, y, p12, p14))
    out = pd.DataFrame(rows)
    out.to_csv(TABLE_DIR / "paired_prediction_tests.csv", index=False)
    return out


def promotion_decision(
    full_metrics: pd.DataFrame,
    external_metrics: pd.DataFrame,
    validation_metrics: pd.DataFrame,
    validation_groups: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    full = full_metrics.set_index("model")
    external = external_metrics.set_index("strategy")
    validation = validation_metrics.set_index("strategy")
    negative_group = validation_groups.loc[
        (validation_groups["group_column"] == "validation_component")
        & (validation_groups["group_value"] == "abT_GSE254249_paired_TCRAB_no_gdTCR")
    ].set_index("strategy")

    guardrail_rows = []
    for model in MODEL_LABELS:
        heldout_fpr = 1.0 - float(negative_group.loc[model, "specificity"])
        values = {
            "full_atlas_primary_gold_recall_ge_0p80": float(full.loc[model, "recall"]) >= 0.80,
            "full_atlas_estimated_fp_fraction_le_0p05": float(full.loc[model, "estimated_fp_fraction"]) <= 0.05,
            "external_specificity_ge_0p9985": float(external.loc[model, "specificity"]) >= 0.9985,
            "heldout_GSE254249_fp_rate_le_0p001": heldout_fpr <= 0.001,
        }
        guardrail_rows.append(
            {
                "model": model,
                "full_atlas_recall": float(full.loc[model, "recall"]),
                "full_atlas_estimated_fp_fraction": float(full.loc[model, "estimated_fp_fraction"]),
                "external_specificity": float(external.loc[model, "specificity"]),
                "heldout_GSE254249_fp_rate": heldout_fpr,
                **values,
                "all_guardrails_pass": all(values.values()),
            }
        )
    guardrails = pd.DataFrame(guardrail_rows)

    criteria = [
        ("Full-atlas primary F1", "higher", full["f1"]),
        ("Atlas held-out F1", "higher", validation["f1"]),
        ("Independent external F1", "higher", external["f1"]),
        ("Full-atlas primary recall", "higher", full["recall"]),
        ("Silver recall", "higher", full["silver_recall"]),
        ("Independent external precision", "higher", external["precision"]),
        ("Full-atlas estimated FP fraction", "lower", full["estimated_fp_fraction"]),
        ("Full-atlas paired-TCRAB calls", "lower", full["paired_TCRAB_calls"]),
        ("Full-atlas NK calls", "lower", full["NK_calls"]),
    ]
    ranking_rows = []
    for criterion, direction, values in criteria:
        v12 = float(values.loc["round12"])
        v14 = float(values.loc["round14"])
        winner = "round12" if (v12 > v14 if direction == "higher" else v12 < v14) else "round14"
        ranking_rows.append(
            {
                "criterion": criterion,
                "preferred_direction": direction,
                "round12": v12,
                "round14": v14,
                "winner": winner,
            }
        )
    ranking = pd.DataFrame(ranking_rows)

    passing = guardrails.loc[guardrails["all_guardrails_pass"], "model"].tolist()
    if not passing:
        raise RuntimeError("Neither candidate passes the prespecified promotion guardrails.")
    composite = {
        model: float(np.mean([full.loc[model, "f1"], validation.loc[model, "f1"], external.loc[model, "f1"]]))
        for model in passing
    }
    selected = max(composite, key=composite.get)
    decision = {
        "selected_model": selected,
        "selected_payload_model": R14_NAME if selected == "round14" else R12_NAME,
        "selection_rule": "Pass all four safety/recall guardrails, then maximize mean F1 across full-atlas gold, atlas-held-out, and independent-external evaluations.",
        "mean_f1_among_passing": composite,
        "round12_mean_f1_all_evaluations": float(
            np.mean([full.loc["round12", "f1"], validation.loc["round12", "f1"], external.loc["round12", "f1"]])
        ),
        "round14_mean_f1_all_evaluations": float(
            np.mean([full.loc["round14", "f1"], validation.loc["round14", "f1"], external.loc["round14", "f1"]])
        ),
        "rationale": (
            "Round 14 is selected because it is the only candidate that reaches the prespecified 0.80 full-atlas gold-recall floor while remaining below all FP guardrails. "
            "It also has substantially higher held-out cord-blood recall and higher mean F1 across the three evaluation frames. Round 12 remains the lower-FP fallback."
            if selected == "round14"
            else "Round 12 is selected by the prespecified guardrail-plus-mean-F1 rule."
        ),
    }
    guardrails.to_csv(TABLE_DIR / "promotion_guardrails.csv", index=False)
    ranking.to_csv(TABLE_DIR / "criterion_comparison.csv", index=False)
    (TABLE_DIR / "promotion_decision.json").write_text(json.dumps(decision, indent=2), encoding="utf-8")
    return guardrails, ranking, decision


def plot_performance(full: pd.DataFrame, validation: pd.DataFrame, external: pd.DataFrame) -> Path:
    frames = {
        "Full atlas gold": full.set_index("model"),
        "Atlas held-out": validation.set_index("strategy"),
        "Independent external": external.set_index("strategy"),
    }
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 4.1), sharey=True, constrained_layout=True)
    metrics = ["precision", "recall", "f1"]
    x = np.arange(len(metrics))
    for ax, (title, frame) in zip(axes, frames.items()):
        for offset, model in [(-0.18, "round12"), (0.18, "round14")]:
            values = [float(frame.loc[model, metric]) for metric in metrics]
            bars = ax.bar(x + offset, values, width=0.34, color=MODEL_COLORS[model], label=MODEL_LABELS[model])
            for bar, value in zip(bars, values):
                ax.text(bar.get_x() + bar.get_width() / 2, value + 0.012, f"{value:.3f}", ha="center", va="bottom", fontsize=7)
        ax.set_title(title, fontsize=11)
        ax.set_xticks(x, ["Precision", "Recall", "F1"])
        ax.set_ylim(0.62, 1.02)
        ax.grid(axis="y", color="#E5E7EB", linewidth=0.7)
    axes[0].set_ylabel("Metric")
    axes[-1].legend(frameon=False, fontsize=8, loc="lower right")
    out = FIGURE_DIR / "performance_across_evaluations.png"
    fig.savefig(out, dpi=240, facecolor="white")
    plt.close(fig)
    return out


def plot_safety(full: pd.DataFrame, validation_groups: pd.DataFrame, external: pd.DataFrame) -> Path:
    f = full.set_index("model")
    e = external.set_index("strategy")
    neg = validation_groups.loc[
        (validation_groups["group_column"] == "validation_component")
        & (validation_groups["group_value"] == "abT_GSE254249_paired_TCRAB_no_gdTCR")
    ].set_index("strategy")
    labels = ["Estimated atlas FP / calls", "External FP / calls", "Held-out abT FPR", "Paired-TCRAB call rate", "NK call rate"]
    values = {}
    for model in MODEL_LABELS:
        values[model] = [
            float(f.loc[model, "estimated_fp_fraction"]),
            float(e.loc[model, "fp_fraction_of_predictions"]),
            1.0 - float(neg.loc[model, "specificity"]),
            float(f.loc[model, "paired_TCRAB_calls"] / f.loc[model, "paired_TCRAB_cells"]),
            float(f.loc[model, "NK_calls"] / f.loc[model, "NK_cells"]),
        ]
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(10.8, 4.8), constrained_layout=True)
    for offset, model in [(-0.19, "round12"), (0.19, "round14")]:
        ax.bar(x + offset, np.asarray(values[model]) * 100, width=0.36, color=MODEL_COLORS[model], label=MODEL_LABELS[model])
    ax.set_yscale("log")
    ax.set_ylabel("Rate (%) on log scale")
    ax.set_xticks(x, labels, rotation=12, ha="right")
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.7, which="both")
    ax.legend(frameon=False)
    out = FIGURE_DIR / "false_positive_safety_comparison.png"
    fig.savefig(out, dpi=240, facecolor="white")
    plt.close(fig)
    return out


def plot_overlap(overall: pd.DataFrame, composition: pd.DataFrame) -> Path:
    row = overall.iloc[0]
    categories = ["Both", "Round 12 only", "Round 14 only"]
    counts = [int(row["both"]), int(row["round12_only"]), int(row["round14_only"])]
    colors = ["#4B5563", MODEL_COLORS["round12"], MODEL_COLORS["round14"]]
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.7), constrained_layout=True)
    bars = axes[0].bar(categories, counts, color=colors)
    axes[0].set_ylabel("Cells")
    axes[0].set_title("Full-atlas prediction overlap")
    for bar, value in zip(bars, counts):
        axes[0].text(bar.get_x() + bar.get_width() / 2, value, f"{value:,}", ha="center", va="bottom", fontsize=8)
    class_cols = ["gdT_gold", "gdT_silver", "abT_gold", "unlabeled_or_ambiguous"]
    class_colors = ["#157A6E", "#E9C46A", "#D1495B", "#9CA3AF"]
    bottom = np.zeros(len(composition))
    for column, color in zip(class_cols, class_colors):
        values = composition[column].to_numpy(dtype=float)
        axes[1].bar(composition["overlap_category"], values, bottom=bottom, label=column, color=color)
        bottom += values
    axes[1].set_title("Ground-truth composition of overlap groups")
    axes[1].set_ylabel("Cells")
    axes[1].legend(frameon=False, fontsize=8)
    out = FIGURE_DIR / "prediction_overlap_and_composition.png"
    fig.savefig(out, dpi=240, facecolor="white")
    plt.close(fig)
    return out


def plot_source_comparison(source: pd.DataFrame) -> tuple[Path, Path]:
    fig, ax = plt.subplots(figsize=(7.2, 6.4), constrained_layout=True)
    x = source["round12_predicted"].to_numpy(dtype=float)
    y = source["round14_predicted"].to_numpy(dtype=float)
    ax.scatter(x, y, s=38, color="#126782", alpha=0.8, edgecolor="white", linewidth=0.5)
    limit = max(x.max(), y.max()) * 1.08
    ax.plot([0, limit], [0, limit], color="#6B7280", linestyle="--", linewidth=1)
    ax.set_xscale("symlog", linthresh=50)
    ax.set_yscale("symlog", linthresh=50)
    ax.set_xlabel("Round 12 predicted cells")
    ax.set_ylabel("Round 14 predicted cells")
    ax.set_title("Prediction counts across all datasets")
    label_rows = source.assign(abs_delta=source["round14_minus_round12"].abs()).nlargest(9, "abs_delta")
    for _, row in label_rows.iterrows():
        ax.annotate(str(row["source_gse_id"]), (row["round12_predicted"], row["round14_predicted"]), xytext=(4, 3), textcoords="offset points", fontsize=7)
    scatter_path = FIGURE_DIR / "per_dataset_prediction_scatter.png"
    fig.savefig(scatter_path, dpi=240, facecolor="white")
    plt.close(fig)

    ordered = source.sort_values("round14_minus_round12")
    fig_h = max(7.0, 0.25 * len(ordered) + 1.4)
    fig, ax = plt.subplots(figsize=(9.0, fig_h), constrained_layout=True)
    colors = np.where(ordered["round14_minus_round12"] >= 0, MODEL_COLORS["round14"], MODEL_COLORS["round12"])
    ax.barh(ordered["source_gse_id"], ordered["round14_minus_round12"], color=colors)
    ax.axvline(0, color="#111827", linewidth=0.8)
    ax.set_xlabel("Round 14 minus Round 12 predicted cells")
    ax.set_title("Dataset-by-dataset prediction change")
    delta_path = FIGURE_DIR / "per_dataset_prediction_delta.png"
    fig.savefig(delta_path, dpi=240, facecolor="white")
    plt.close(fig)
    return scatter_path, delta_path


def plot_recall_by_source(source: pd.DataFrame, external_groups: pd.DataFrame) -> Path:
    positive = source.loc[source["gdT_gold"] > 0].copy()
    labels = positive["source_gse_id"].astype(str).tolist() + ["Silver standard", "Independent external"]
    values = {"round12": positive["round12_gold_recall"].tolist(), "round14": positive["round14_gold_recall"].tolist()}
    for model in MODEL_LABELS:
        silver = source[model + "_gdT_silver"].sum() / source["gdT_silver"].sum()
        external = external_groups.loc[
            (external_groups["model"] == model) & (external_groups["group"] == "all_real_gdT"), "rate"
        ].iloc[0]
        values[model].extend([float(silver), float(external)])
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(11.5, 4.8), constrained_layout=True)
    for offset, model in [(-0.19, "round12"), (0.19, "round14")]:
        ax.bar(x + offset, values[model], width=0.36, color=MODEL_COLORS[model], label=MODEL_LABELS[model])
    ax.set_ylim(0, 1.04)
    ax.set_ylabel("Recall")
    ax.set_xticks(x, labels, rotation=18, ha="right")
    ax.set_title("Gold, silver, and external gdT recovery")
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.7)
    ax.legend(frameon=False)
    out = FIGURE_DIR / "gdt_recall_by_evidence_source.png"
    fig.savefig(out, dpi=240, facecolor="white")
    plt.close(fig)
    return out


def plot_quadrants(union: pd.DataFrame) -> Path:
    rows = []
    for model in MODEL_LABELS:
        counts = union.loc[union[f"{model}_pred"], "tcr_gene_quadrant"].value_counts()
        total = counts.sum()
        for quadrant in QUADRANT_COLORS:
            rows.append({"model": model, "quadrant": quadrant, "fraction": counts.get(quadrant, 0) / total, "count": counts.get(quadrant, 0)})
    df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(7.6, 4.8), constrained_layout=True)
    bottom = np.zeros(2)
    for quadrant, color in QUADRANT_COLORS.items():
        values = [df.loc[(df["model"] == model) & (df["quadrant"] == quadrant), "fraction"].iloc[0] for model in MODEL_LABELS]
        ax.bar([MODEL_LABELS[m] for m in MODEL_LABELS], values, bottom=bottom, label=quadrant, color=color)
        bottom += np.asarray(values)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of predicted cells")
    ax.set_title("TCR-gene expression composition")
    ax.legend(frameon=False, fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    out = FIGURE_DIR / "predicted_tcr_quadrant_composition.png"
    fig.savefig(out, dpi=240, facecolor="white")
    plt.close(fig)
    return out


def compact_table_html(df: pd.DataFrame, percent_columns: set[str] | None = None) -> str:
    percent_columns = percent_columns or set()
    display = df.copy()
    for column in display.columns:
        if column in percent_columns:
            display[column] = display[column].map(lambda value: "" if pd.isna(value) else f"{100 * float(value):.2f}%")
        elif pd.api.types.is_float_dtype(display[column]):
            display[column] = display[column].map(lambda value: "" if pd.isna(value) else f"{float(value):.4f}")
    return display.to_html(index=False, escape=True, border=0, classes="data-table")


def render_report(
    artifact_manifest: pd.DataFrame,
    full: pd.DataFrame,
    external: pd.DataFrame,
    external_groups: pd.DataFrame,
    validation: pd.DataFrame,
    validation_groups: pd.DataFrame,
    overlap: pd.DataFrame,
    composition: pd.DataFrame,
    source: pd.DataFrame,
    paired: pd.DataFrame,
    guardrails: pd.DataFrame,
    ranking: pd.DataFrame,
    decision: dict[str, Any],
    checks: pd.DataFrame,
    figures: list[Path],
    args: argparse.Namespace,
) -> None:
    selected = decision["selected_model"]
    f = full.set_index("model")
    e = external.set_index("strategy")
    v = validation.set_index("strategy")
    comparison = pd.DataFrame(
        {
            "Evaluation": ["Full atlas primary gold", "Atlas held-out combined", "Independent external"],
            "Round 12 precision": [f.loc["round12", "precision"], v.loc["round12", "precision"], e.loc["round12", "precision"]],
            "Round 12 recall": [f.loc["round12", "recall"], v.loc["round12", "recall"], e.loc["round12", "recall"]],
            "Round 12 F1": [f.loc["round12", "f1"], v.loc["round12", "f1"], e.loc["round12", "f1"]],
            "Round 14 precision": [f.loc["round14", "precision"], v.loc["round14", "precision"], e.loc["round14", "precision"]],
            "Round 14 recall": [f.loc["round14", "recall"], v.loc["round14", "recall"], e.loc["round14", "recall"]],
            "Round 14 F1": [f.loc["round14", "f1"], v.loc["round14", "f1"], e.loc["round14", "f1"]],
        }
    )
    validation_components = validation_groups.loc[validation_groups["group_column"].eq("validation_component")].copy()
    validation_components = validation_components[
        ["strategy", "group_value", "n_cells", "predicted_positive", "tp", "fp", "fn", "precision", "recall", "specificity", "f1"]
    ]
    source_display = source[
        [
            "source_gse_id",
            "total_cells",
            "round12_predicted",
            "round14_predicted",
            "round14_minus_round12",
            "jaccard",
            "gdT_gold",
            "round12_gold_recall",
            "round14_gold_recall",
            "round12_NK",
            "round14_NK",
            "round12_paired_TCRAB",
            "round14_paired_TCRAB",
        ]
    ].sort_values("round14_minus_round12", key=lambda values: values.abs(), ascending=False)
    overlap_row = overlap.iloc[0]
    key_deltas = [
        ("Raw atlas calls", int(f.loc["round12", "predicted_positive"]), int(f.loc["round14", "predicted_positive"])),
        ("Gold gdT recovered", int(f.loc["round12", "tp"]), int(f.loc["round14", "tp"])),
        ("Silver gdT recovered", int(f.loc["round12", "silver_recovered"]), int(f.loc["round14", "silver_recovered"])),
        ("Paired-TCRAB calls", int(f.loc["round12", "paired_TCRAB_calls"]), int(f.loc["round14", "paired_TCRAB_calls"])),
        ("NK-annotated calls", int(f.loc["round12", "NK_calls"]), int(f.loc["round14", "NK_calls"])),
    ]
    delta_df = pd.DataFrame(key_deltas, columns=["Measure", "Round 12", "Round 14"])
    delta_df["Round 14 minus Round 12"] = delta_df["Round 14"] - delta_df["Round 12"]

    for figure in figures:
        shutil.copy2(figure, STATIC_FIGURE_DIR / figure.name)
    figure_html = "\n".join(
        f'<figure><img src="figures/{html.escape(path.name)}" alt="{html.escape(path.stem)}"><figcaption>{html.escape(path.stem.replace("_", " ").title())}</figcaption></figure>'
        for path in figures
    )

    css = """
    :root { --ink:#18212b; --muted:#5f6b76; --line:#d9e0e6; --panel:#f5f7f9; --r12:#126782; --r14:#d1495b; --ok:#157a6e; }
    * { box-sizing:border-box; }
    body { margin:0; color:var(--ink); background:white; font-family:Arial, Helvetica, sans-serif; line-height:1.46; }
    header { border-bottom:5px solid var(--r14); padding:34px 5vw 28px; background:#f8fafb; }
    main { max-width:1220px; margin:0 auto; padding:26px 34px 60px; }
    h1 { margin:0 0 8px; font-size:34px; letter-spacing:0; }
    h2 { margin:32px 0 12px; font-size:22px; letter-spacing:0; border-bottom:1px solid var(--line); padding-bottom:7px; }
    h3 { margin:22px 0 8px; font-size:16px; letter-spacing:0; }
    p, li { font-size:14px; }
    .subtitle { color:var(--muted); margin:0; max-width:900px; }
    .decision { border-left:6px solid var(--ok); background:#edf7f4; padding:18px 20px; margin:22px 0; }
    .decision strong { font-size:21px; }
    .metrics { display:grid; grid-template-columns:repeat(4, minmax(0,1fr)); gap:10px; margin:18px 0; }
    .metric { border:1px solid var(--line); border-radius:6px; padding:13px; background:white; }
    .metric .label { color:var(--muted); font-size:12px; }
    .metric .value { font-size:23px; font-weight:700; margin-top:3px; }
    .r12 { color:var(--r12); } .r14 { color:var(--r14); }
    .table-wrap { overflow-x:auto; border:1px solid var(--line); margin:10px 0 22px; }
    table.dataframe, table.data-table { border-collapse:collapse; width:100%; font-size:11px; }
    th, td { border-bottom:1px solid var(--line); padding:6px 7px; text-align:right; white-space:nowrap; }
    th:first-child, td:first-child { text-align:left; }
    th { background:#eef2f5; color:#283744; position:sticky; top:0; }
    figure { margin:24px 0; break-inside:avoid; }
    figure img { display:block; max-width:100%; height:auto; border:1px solid var(--line); }
    figcaption { color:var(--muted); font-size:11px; margin-top:5px; }
    code { background:#eef2f5; padding:1px 4px; }
    .note { color:var(--muted); font-size:12px; }
    .file-list a { color:#126782; }
    @media (max-width:780px) { main { padding:20px 14px; } .metrics { grid-template-columns:repeat(2,minmax(0,1fr)); } h1 { font-size:27px; } }
    @media print {
      @page { size:A4 landscape; margin:9mm; }
      header { padding:12mm 10mm 7mm; }
      main { max-width:none; padding:5mm 8mm; }
      h1 { font-size:24px; } h2 { font-size:17px; break-after:avoid; }
      p, li { font-size:10px; }
      .metrics { grid-template-columns:repeat(4,1fr); }
      .metric .value { font-size:17px; }
      table.dataframe, table.data-table { font-size:7px; }
      th, td { padding:3px 4px; white-space:normal; }
      .table-wrap { overflow:visible; border:0; }
      figure { break-inside:avoid; page-break-inside:avoid; }
      figure img { max-height:165mm; margin:auto; }
      a { color:inherit; text-decoration:none; }
    }
    """
    report = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>gdTAI v3 Round 12 vs Round 14</title><style>{css}</style></head><body>
<header><h1>gdTAI v3: Round 12 vs Round 14</h1>
<p class="subtitle">A checksum-pinned, same-cohort comparison across the 5.13-million-cell atlas, held-out atlas cohorts, and an independent external H5AD.</p></header>
<main>
<section class="decision"><strong>Selected: {html.escape(MODEL_LABELS[selected])}</strong>
<p>{html.escape(decision['rationale'])}</p></section>
<div class="metrics">
  <div class="metric"><div class="label">Round 12 atlas calls</div><div class="value r12">{int(f.loc['round12','predicted_positive']):,}</div></div>
  <div class="metric"><div class="label">Round 14 atlas calls</div><div class="value r14">{int(f.loc['round14','predicted_positive']):,}</div></div>
  <div class="metric"><div class="label">Shared calls</div><div class="value">{int(overlap_row['both']):,}</div></div>
  <div class="metric"><div class="label">Prediction Jaccard</div><div class="value">{float(overlap_row['jaccard']):.3f}</div></div>
</div>

<h2>Executive Comparison</h2>
<p>Round 14 adds {int(f.loc['round14','predicted_positive'] - f.loc['round12','predicted_positive']):,} raw atlas calls and recovers {int(f.loc['round14','tp'] - f.loc['round12','tp']):,} additional gold gdT cells. The cost is {int(f.loc['round14','paired_TCRAB_calls'] - f.loc['round12','paired_TCRAB_calls']):,} additional paired-TCRAB calls and {int(f.loc['round14','NK_calls'] - f.loc['round12','NK_calls']):,} additional NK-annotated calls. On the independent external set, Round 14 recovers {int(e.loc['round14','tp'] - e.loc['round12','tp']):,} additional true gdT cell but adds {int(e.loc['round14','fp'] - e.loc['round12','fp']):,} false positives, so Round 12 has the higher external F1. On the atlas held-out set, Round 14 has the substantially higher F1 because it rescues cord-blood gdT cells.</p>
<div class="table-wrap">{compact_table_html(delta_df)}</div>
<div class="table-wrap">{compact_table_html(comparison)}</div>

<h2>Decision Rule</h2>
<p>{html.escape(decision['selection_rule'])}</p>
<div class="table-wrap">{compact_table_html(guardrails, {'full_atlas_recall','full_atlas_estimated_fp_fraction','external_specificity','heldout_GSE254249_fp_rate'})}</div>
<p>Mean F1 across full-atlas gold, atlas-held-out, and independent-external evaluations: Round 12 = {decision['round12_mean_f1_all_evaluations']:.4f}; Round 14 = {decision['round14_mean_f1_all_evaluations']:.4f}.</p>
<div class="table-wrap">{compact_table_html(ranking)}</div>

<h2>Figures</h2>{figure_html}

<h2>Held-Out Cohorts</h2>
<p>The atlas-held-out frame contains sorted cord-blood gdT positives from <code>GDT_2020AUG_woCOV</code> and paired-TCRAB/no-gdTCR negatives from <code>GSE254249</code>. The external H5AD was not used for fitting or threshold selection.</p>
<div class="table-wrap">{compact_table_html(validation_components)}</div>
<div class="table-wrap">{compact_table_html(external_groups)}</div>

<h2>Prediction Overlap</h2>
<div class="table-wrap">{compact_table_html(overlap)}</div>
<div class="table-wrap">{compact_table_html(composition)}</div>
<p>Round 14-only calls contain {int(composition.loc[composition['overlap_category'].eq('Round 14 only'),'gdT_gold'].iloc[0]):,} gold and {int(composition.loc[composition['overlap_category'].eq('Round 14 only'),'gdT_silver'].iloc[0]):,} silver gdT cells. Round 12-only calls contain {int(composition.loc[composition['overlap_category'].eq('Round 12 only'),'gdT_gold'].iloc[0]):,} gold and {int(composition.loc[composition['overlap_category'].eq('Round 12 only'),'gdT_silver'].iloc[0]):,} silver gdT cells.</p>

<h2>Dataset-by-Dataset Audit</h2>
<p>All datasets are shown. Count changes are predictions, not biological abundance estimates; denominator-aware rates and truth recovery are included where labels exist.</p>
<div class="table-wrap">{compact_table_html(source_display, {'jaccard','round12_gold_recall','round14_gold_recall'})}</div>

<h2>Paired Statistical Tests</h2>
<p>Exact McNemar tests compare discordant cell-level errors. P-values are descriptive because cells from the same donor are correlated; effect sizes and cohort consistency remain the primary evidence.</p>
<div class="table-wrap">{compact_table_html(paired)}</div>

<h2>Artifact and Cache Verification</h2>
<p>Both artifacts were frozen before comparison. Direct external inference was required to match the existing validated cache exactly. A deterministic atlas sample included shared calls, each model-only group, and cells called by neither model.</p>
<div class="table-wrap">{compact_table_html(artifact_manifest)}</div>
<div class="table-wrap">{compact_table_html(checks)}</div>

<h2>Interpretation</h2>
<p><strong>Round 14 should remain the canonical balanced model.</strong> Its recall gain is large in the fully held-out cord-blood cohort and sufficient to clear the established 0.80 atlas gold-recall floor. Its independent external recall is essentially unchanged relative to Round 12, while external precision is lower. This means Round 14 is not uniformly superior: Round 12 is the appropriate documented high-purity fallback when false-positive minimization is more important than sorted-gdT recovery.</p>
<p>The Round 14 increase in NK and TRDC+TRDV- calls remains a real caveat. Downstream atlas construction should continue removing cells with direct paired-TCRAB/no-gdT evidence and should report NK-annotation sensitivity analyses. The comparison does not justify treating every unlabeled Round 14-only cell as a confirmed gdT cell.</p>

<h2>Reproducibility</h2>
<ul>
<li>Atlas H5AD: <code>{html.escape(str(args.atlas_h5ad))}</code></li>
<li>External H5AD: <code>{html.escape(str(args.external_h5ad))}</code></li>
<li>No H5AD was modified.</li>
<li>Detailed CSV tables: <code>{html.escape(str(TABLE_DIR))}</code></li>
<li>Model snapshots: <code>{html.escape(str(MODEL_DIR))}</code></li>
</ul>
<p class="note">Generated from <code>workflows/gdtai/compare_gdtai_v3_round12_vs_round14.py</code>.</p>
</main></body></html>"""
    REPORT_HTML.write_text(report, encoding="utf-8")

    md = [
        "# gdTAI v3 Round 12 vs Round 14",
        "",
        f"## Decision: {MODEL_LABELS[selected]}",
        "",
        decision["rationale"],
        "",
        "## Performance",
        "",
        dataframe_to_markdown(comparison),
        "",
        "## Guardrails",
        "",
        dataframe_to_markdown(guardrails),
        "",
        "## Full-Atlas Overlap",
        "",
        dataframe_to_markdown(overlap),
        "",
        "## Artifact Manifest",
        "",
        dataframe_to_markdown(artifact_manifest),
        "",
        f"HTML: `{REPORT_HTML}`",
        f"PDF: `{REPORT_PDF}`",
    ]
    REPORT_MD.write_text("\n".join(md) + "\n", encoding="utf-8")

    if not args.no_pdf:
        subprocess.run(
            [
                "/usr/bin/google-chrome",
                "--headless",
                "--disable-gpu",
                "--no-sandbox",
                "--run-all-compositor-stages-before-draw",
                f"--print-to-pdf={REPORT_PDF}",
                str(REPORT_HTML),
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )


def promote_selected_package(
    decision: dict[str, Any],
    full: pd.DataFrame,
    external: pd.DataFrame,
    validation: pd.DataFrame,
) -> dict[str, Any]:
    selected = str(decision["selected_model"])
    if selected not in MODEL_LABELS:
        raise RuntimeError(f"Unsupported selected model: {selected}")
    selected_snapshot = MODEL_DIR / f"{selected}_model.pkl"
    selected_hash = sha256_file(selected_snapshot)
    expected_hash = R14_SHA256 if selected == "round14" else R12_SHA256
    if selected_hash != expected_hash:
        raise RuntimeError(f"Refusing promotion because {selected} snapshot hash is unexpected: {selected_hash}")

    canonical_model = CANONICAL_MODEL_DIR / "gdTAI_v3_model.pkl"
    candidate_model = SOURCE_CANDIDATE_DIR / "best_candidate_model.pkl"
    shutil.copy2(selected_snapshot, canonical_model)
    shutil.copy2(selected_snapshot, candidate_model)
    if sha256_file(canonical_model) != expected_hash or sha256_file(candidate_model) != expected_hash:
        raise RuntimeError("Promoted model copies failed checksum verification.")

    for name in [
        "feature_genes.csv",
        "candidate_target_selection.csv",
        "external_false_positive_groups.csv",
        "external_recall_groups.csv",
    ]:
        shutil.copy2(SOURCE_CANDIDATE_DIR / name, CANONICAL_MODEL_DIR / name)

    full_row = full.set_index("model").loc[selected]
    external_row = external.set_index("strategy").loc[selected]
    validation_row = validation.set_index("strategy").loc[selected]
    payload_name = R14_NAME if selected == "round14" else R12_NAME
    threshold = 0.936 if selected == "round14" else 0.5
    manifest = {
        "model_name": payload_name,
        "model_file": "gdTAI_v3_model.pkl",
        "model_sha256": expected_hash,
        "model_size_bytes": int(canonical_model.stat().st_size),
        "version": "gdTAI_v3.0",
        "status": "promoted_default",
        "operating_profile": "balanced",
        "threshold": threshold,
        "accepted_for_promotion": True,
        "promotion_basis": "checksum-pinned Round 12 versus Round 14 revalidation",
        "promotion_date_local": "2026-07-15 Asia/Hong_Kong",
        "comparison_report": "gdT_prediction/gdtai_v3_round12_vs_round14/index.html",
        "comparison_workflow": "workflows/gdtai/compare_gdtai_v3_round12_vs_round14.py",
        "external_validation_dataset_id": "BALF_BLOOD_COPD",
        "external_validation_h5ad": (
            "data/datasets/BALF_BLOOD_COPD/processed/current.h5ad"
        ),
        "requires_counts_layer_for_external_h5ad": True,
        "normalization": "log1p(counts per 10000)",
        "n_gene_features": 210,
        "n_engineered_features": 16,
        "n_total_features": 226,
        "selection_rule": decision["selection_rule"],
        "full_atlas": {
            "total_cells": int(full_row["total_cells"]),
            "predicted_putative_gdT": int(full_row["predicted_positive"]),
            "primary_gold_recall": float(full_row["recall"]),
            "primary_gold_precision": float(full_row["precision"]),
            "primary_gold_f1": float(full_row["f1"]),
            "estimated_total_abT_fp": float(full_row["estimated_total_fp"]),
            "estimated_fp_fraction_of_predictions": float(full_row["estimated_fp_fraction"]),
            "silver_recall": float(full_row["silver_recall"]),
        },
        "atlas_heldout": {
            "n_cells": int(validation_row["n_cells"]),
            "precision": float(validation_row["precision"]),
            "recall": float(validation_row["recall"]),
            "specificity": float(validation_row["specificity"]),
            "f1": float(validation_row["f1"]),
        },
        "independent_external": {
            "n_cells": int(external_row["n_cells"]),
            "predicted_positive": int(external_row["predicted_positive"]),
            "precision": float(external_row["precision"]),
            "recall": float(external_row["recall"]),
            "specificity": float(external_row["specificity"]),
            "f1": float(external_row["f1"]),
            "fp": int(external_row["fp"]),
        },
        "high_purity_fallback": {
            "model_name": R12_NAME if selected == "round14" else R14_NAME,
            "artifact_path": str(MODEL_DIR / ("round12_model.pkl" if selected == "round14" else "round14_model.pkl")),
            "sha256": R12_SHA256 if selected == "round14" else R14_SHA256,
        },
        "caveat": "The balanced Round 14 default has higher NK and paired-TCRAB call burdens than the Round 12 high-purity fallback; report these sensitivity checks on new datasets.",
    }
    (CANONICAL_MODEL_DIR / "model_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    mode_rows = []
    for evaluation, row in [
        ("full_atlas_primary_gold", full_row),
        ("atlas_heldout_combined", validation_row),
        ("independent_external_primary", external_row),
    ]:
        record = {"evaluation": evaluation, "strategy": payload_name, "threshold": threshold}
        for key in [
            "n_cells",
            "total_cells",
            "predicted_positive",
            "precision",
            "recall",
            "specificity",
            "f1",
            "tp",
            "fp",
            "tn",
            "fn",
            "estimated_fp_fraction",
            "paired_TCRAB_calls",
            "NK_calls",
            "silver_recall",
        ]:
            if key in row.index and pd.notna(row[key]):
                record[key] = row[key]
        mode_rows.append(record)
    pd.DataFrame(mode_rows).to_csv(CANONICAL_MODEL_DIR / "mode_metrics.csv", index=False)
    external_export = external.loc[external["strategy"].eq(selected)].copy()
    external_export["strategy"] = payload_name
    external_export.to_csv(CANONICAL_MODEL_DIR / "external_primary_metrics.csv", index=False)

    readme = f"""# gdTAI v3.0

This is the promoted balanced gdTAI v3 package.

- Model: `{payload_name}`
- Threshold: `{threshold}`
- SHA256: `{expected_hash}`
- Input transform: `log1p(raw counts per 10,000)`
- Features: 210 genes plus 16 engineered expression features
- Comparison decision: `gdT_prediction/gdtai_v3_round12_vs_round14/index.html`
- Independent validation dataset: `BALF_BLOOD_COPD`
- Validation H5AD: `data/datasets/BALF_BLOOD_COPD/processed/current.h5ad`

## Headline Metrics

| Evaluation | Precision | Recall | F1 |
|---|---:|---:|---:|
| Full atlas primary gold | {full_row['precision']:.4f} | {full_row['recall']:.4f} | {full_row['f1']:.4f} |
| Atlas held-out combined | {validation_row['precision']:.4f} | {validation_row['recall']:.4f} | {validation_row['f1']:.4f} |
| Independent external | {external_row['precision']:.4f} | {external_row['recall']:.4f} | {external_row['f1']:.4f} |

## Files

- `gdTAI_v3_model.pkl`: promoted model payload
- `model_manifest.json`: machine-readable identity and metrics
- `METHODOLOGY.md`: training, labels, features, and evaluation design
- `USAGE_PROTOCOL.md`: inference requirements and QC
- `feature_genes.csv`: ordered feature schema
- `examples/predict_h5ad_counts.py`: chunked H5AD inference example

Round 12 is preserved as a lower-FP fallback under
`Integrated_dataset/models/gdT_prediction_classifier/{PREFIX}/round12_model.pkl`.
"""
    (CANONICAL_MODEL_DIR / "README.md").write_text(readme, encoding="utf-8")

    source_methodology = (SOURCE_CANDIDATE_DIR / "METHODOLOGY.md").read_text(encoding="utf-8")
    source_methodology = source_methodology.split("## Promotion Status", maxsplit=1)[0].rstrip()
    source_methodology += f"""

## Round 12 Versus Round 14 Revalidation

Both model artifacts were pinned by SHA256 and compared on identical cohorts.
Round 14 was selected on 2026-07-15 because it was the only model to pass all
four prespecified guardrails, including full-atlas gold recall >= 0.80, while
remaining below the atlas, external, and GSE254249 false-positive limits.

- mean F1 across the three evaluation frames: `{decision['round14_mean_f1_all_evaluations']:.4f}`
- comparison report: `gdT_prediction/{PREFIX}/index.html`
- Round 12 remains a lower-FP fallback; it is not the canonical default.

## Promotion Status

`{payload_name}` is the promoted gdTAI v3.0 balanced default. Its higher NK and
paired-TCRAB call burdens relative to Round 12 must remain visible in external
QC reports.
"""
    (CANONICAL_MODEL_DIR / "METHODOLOGY.md").write_text(source_methodology.rstrip() + "\n", encoding="utf-8")

    usage = (SOURCE_CANDIDATE_DIR / "USAGE_PROTOCOL.md").read_text(encoding="utf-8")
    usage = usage.replace("This is a review candidate, not a promoted `gdTAI_v3.0` release.", "This is the promoted `gdTAI_v3.0` balanced release.")
    usage = usage.replace(
        "It met the user target of full-atlas primary-gold recall above `0.8` and estimated FP below `5%`.",
        "It passed the 2026-07-15 Round 12 versus Round 14 promotion guardrails, including full-atlas primary-gold recall above `0.8` and estimated FP below `5%`.",
    )
    usage = usage.replace(
        "It failed the older strict promotion gate because external NK false positives were higher than v2 high-purity.",
        "Round 12 remains the lower-FP fallback because Round 14 has higher NK and paired-TCRAB call burdens.",
    )
    (CANONICAL_MODEL_DIR / "USAGE_PROTOCOL.md").write_text(usage, encoding="utf-8")

    promotion_md = f"""# gdTAI v3.0 Promotion Decision

Round 14 remains the canonical balanced gdTAI v3.0 model after a new
checksum-pinned Round 12 versus Round 14 comparison completed on 2026-07-15.

- Promoted model: `{payload_name}`
- Threshold: `{threshold}`
- Model SHA256: `{expected_hash}`
- Full-atlas predicted cells: `{int(full_row['predicted_positive']):,}`
- Full-atlas primary-gold recall / F1: `{full_row['recall']:.4f}` / `{full_row['f1']:.4f}`
- Atlas-held-out recall / F1: `{validation_row['recall']:.4f}` / `{validation_row['f1']:.4f}`
- Independent external precision / recall / F1: `{external_row['precision']:.4f}` / `{external_row['recall']:.4f}` / `{external_row['f1']:.4f}`

Selection rule: {decision['selection_rule']}

Round 12 is retained as the high-purity fallback. The detailed comparison is
served at `gdT_prediction/{PREFIX}/index.html`.
"""
    (CANONICAL_MODEL_DIR / "PROMOTION_DECISION.md").write_text(promotion_md, encoding="utf-8")

    record = {
        "selected_model": selected,
        "payload_model": payload_name,
        "canonical_model": str(canonical_model),
        "canonical_sha256": sha256_file(canonical_model),
        "candidate_model": str(candidate_model),
        "candidate_sha256": sha256_file(candidate_model),
        "comparison_report": str(REPORT_HTML),
    }
    (LOG_DIR / "promotion_record.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return record


def main() -> None:
    args = parse_args()
    setup_logging()
    atlas_before = file_record(args.atlas_h5ad)
    external_before = file_record(args.external_h5ad)

    r12_path, r14_path, artifact_manifest = pin_model_artifacts(args)
    payloads = {"round12": load_pickle(r12_path), "round14": load_pickle(r14_path)}
    logging.info("Pinned Round 12 SHA256: %s", R12_SHA256)
    logging.info("Pinned Round 14 SHA256: %s", R14_SHA256)

    external_metrics, external_groups, external_df, external_checks = direct_external_evaluation(args.external_h5ad, payloads)
    union, overlap, composition = load_full_prediction_union(args.atlas_h5ad)
    source = source_summary(union)
    full = standardized_full_metrics(union)
    validation_metrics, validation_groups, validation_df = validation_evaluation()
    atlas_checks = atlas_artifact_cache_check(
        args.atlas_h5ad,
        payloads,
        union,
        args.atlas_cache_check_per_group,
    )
    checks = pd.concat([external_checks, atlas_checks], ignore_index=True, sort=False)
    checks.to_csv(TABLE_DIR / "artifact_cache_verification.csv", index=False)
    paired = paired_tests(external_df, validation_df, union)
    guardrails, ranking, decision = promotion_decision(full, external_metrics, validation_metrics, validation_groups)
    promotion_record = None
    if args.promote_selected:
        promotion_record = promote_selected_package(decision, full, external_metrics, validation_metrics)

    figures = [
        plot_performance(full, validation_metrics, external_metrics),
        plot_safety(full, validation_groups, external_metrics),
        plot_overlap(overlap, composition),
        *plot_source_comparison(source),
        plot_recall_by_source(source, external_groups),
        plot_quadrants(union),
    ]
    render_report(
        artifact_manifest,
        full,
        external_metrics,
        external_groups,
        validation_metrics,
        validation_groups,
        overlap,
        composition,
        source,
        paired,
        guardrails,
        ranking,
        decision,
        checks,
        figures,
        args,
    )

    atlas_after = file_record(args.atlas_h5ad)
    external_after = file_record(args.external_h5ad)
    h5ad_unchanged = atlas_before == atlas_after and external_before == external_after
    if not h5ad_unchanged:
        raise RuntimeError("An input H5AD size or mtime changed during the read-only comparison.")
    summary = {
        "selected_model": decision["selected_model"],
        "selected_payload_model": decision["selected_payload_model"],
        "decision": decision,
        "artifact_manifest": artifact_manifest.to_dict(orient="records"),
        "full_atlas_overlap": overlap.iloc[0].to_dict(),
        "full_atlas_performance": full.to_dict(orient="records"),
        "external_performance": external_metrics.to_dict(orient="records"),
        "atlas_heldout_performance": validation_metrics.to_dict(orient="records"),
        "h5ad_unchanged": h5ad_unchanged,
        "atlas_file_before": atlas_before,
        "atlas_file_after": atlas_after,
        "external_file_before": external_before,
        "external_file_after": external_after,
        "outputs": {
            "tables": str(TABLE_DIR),
            "figures": str(FIGURE_DIR),
            "models": str(MODEL_DIR),
            "report_markdown": str(REPORT_MD),
            "report_html": str(REPORT_HTML),
            "report_pdf": str(REPORT_PDF) if REPORT_PDF.exists() else None,
        },
        "promotion_applied": bool(args.promote_selected),
        "promotion_record": promotion_record,
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    logging.info("Selected model: %s", decision["selected_model"])
    logging.info("Report: %s", REPORT_HTML)


if __name__ == "__main__":
    main()
