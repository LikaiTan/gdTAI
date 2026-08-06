#!/usr/bin/env python3
"""Apply released gdTAI v2 and v3 models to all deposited GSE305372 T cells."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import html
import json
import os
import pickle
import shutil
import sys
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/gdtai_gse305372_all_t_matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import mmread

PROJECT_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_DIR = PROJECT_ROOT / "workflows" / "gdtai"
if str(WORKFLOW_DIR) not in sys.path:
    sys.path.insert(0, str(WORKFLOW_DIR))

from run_gdtai_v3_trdc_nk_guard_classifier import (  # noqa: E402
    FeatureSpec,
    append_engineered_features,
    trdc_trdv_quadrant,
)

DATASET_ID = "GSE305372"
DATASET_ROOT = PROJECT_ROOT / "data" / "datasets" / DATASET_ID
DEFAULT_PAYLOAD_DIR = DATASET_ROOT / "interim" / "all_t_model_payload"
DEFAULT_V2_MODEL = (
    PROJECT_ROOT
    / "Integrated_dataset"
    / "models"
    / "gdT_prediction_classifier"
    / "gdTAI_v2.0"
    / "gdTAI_v2_model.pkl"
)
DEFAULT_V3_MODEL = (
    PROJECT_ROOT
    / "Integrated_dataset"
    / "models"
    / "gdT_prediction_classifier"
    / "gdTAI_v3.0"
    / "gdTAI_v3_model.pkl"
)
DEFAULT_OUTPUT_DIR = (
    PROJECT_ROOT / "Integrated_dataset" / "tables" / "gdT_prediction" / "GSE305372_all_T"
)
DEFAULT_FIGURE_DIR = (
    PROJECT_ROOT / "Integrated_dataset" / "figures" / "gdT_prediction" / "GSE305372_all_T"
)
DEFAULT_LOG_DIR = (
    PROJECT_ROOT / "Integrated_dataset" / "logs" / "gdT_prediction" / "GSE305372_all_T"
)
DEFAULT_REPORT_DIR = PROJECT_ROOT / "reports" / "GSE305372_gdtai_all_t"
MODEL_REGISTRY = PROJECT_ROOT / "configs" / "models" / "gdtai" / "model_registry.csv"
DICE_TCR_RESOURCE = (
    DATASET_ROOT
    / "raw"
    / "source"
    / "GSE305372_DICE-TCR_BL-CD4-ALL_ProcessedObject.csv.gz"
)

OBJECTS = {
    "lung_CD4": {
        "prefix": "LG_CD4",
        "compartment": "lung",
        "lineage": "CD4",
        "cluster_column": "RNA_snn_res.0.2",
        "cluster_labels": {
            "0": "TRM",
            "1": "TCM",
            "2": "CD4-CTL",
            "3": "TREG",
            "4": "TFH",
            "5": "Prolif",
        },
    },
    "lung_CD8": {
        "prefix": "LG_CD8",
        "compartment": "lung",
        "lineage": "CD8",
        "cluster_column": "RNA_snn_res.0.2",
        "cluster_labels": {
            "0": "TRM",
            "1": "TEM",
            "2": "GZMKhi",
            "3": "TCM",
            "4": "NKG2Cpos",
            "5": "MAIT",
            "6": "Prolif",
        },
    },
    "lymph_node_CD4": {
        "prefix": "LN_CD4",
        "compartment": "lymph_node",
        "lineage": "CD4",
        "cluster_column": "RNA_snn_res.0.3",
        "cluster_labels": {
            "0": "TN",
            "1": "TI",
            "2": "TCM",
            "3": "TREG",
            "4": "TREG",
            "5": "IFNR",
            "6": "TFH",
            "7": "TRM",
            "8": "TFH",
        },
    },
    "lymph_node_CD8": {
        "prefix": "LN_CD8",
        "compartment": "lymph_node",
        "lineage": "CD8",
        "cluster_column": "RNA_snn_res.0.3",
        "cluster_labels": {
            "0": "GZMKhi",
            "1": "TI",
            "2": "TCM",
            "3": "TN",
            "4": "IFNR",
            "5": "TRM",
            "6": "TACT",
            "7": "MAIT",
            "8": "TEM",
        },
    },
}

MODEL_COLUMNS = {
    "v2_high_f1": ("gdtai_v2_score", "gdtai_v2_high_f1_predicted"),
    "v2_high_purity": ("gdtai_v2_score", "gdtai_v2_high_purity_predicted"),
    "v3_round14": ("gdtai_v3_score", "gdtai_v3_predicted"),
}

MANUSCRIPT_REPORTED_COUNTS = {
    "lung_CD4": 370_233,
    "lung_CD8": 227_800,
    "lymph_node_CD4": 198_890,
    "lymph_node_CD8": 129_737,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--payload-dir", type=Path, default=DEFAULT_PAYLOAD_DIR)
    parser.add_argument("--v2-model-pkl", type=Path, default=DEFAULT_V2_MODEL)
    parser.add_argument("--v3-model-pkl", type=Path, default=DEFAULT_V3_MODEL)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--log-dir", type=Path, default=DEFAULT_LOG_DIR)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--chunk-size", type=int, default=50_000)
    parser.add_argument("--max-plot-points", type=int, default=120_000)
    parser.add_argument("--seed", type=int, default=20260730)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def inspect_dice_tcr_resource(path: Path, model_genes: list[str]) -> dict[str, Any]:
    if not path.is_file():
        return {"path": str(path.relative_to(PROJECT_ROOT)), "available": False}
    columns = pd.read_csv(path, nrows=0).columns.astype(str).tolist()
    with gzip.open(path, "rb") as handle:
        data_rows = sum(1 for _ in handle) - 1
    return {
        "path": str(path.relative_to(PROJECT_ROOT)),
        "available": True,
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "data_rows": int(data_rows),
        "n_columns": len(columns),
        "columns": columns,
        "model_genes_present_as_columns": sorted(set(columns).intersection(model_genes)),
        "eligible_for_gdtai": False,
        "ineligibility_reason": (
            "TCR/metadata table only; no raw transcriptome matrix or whole-library "
            "count denominator"
        ),
    }


def build_v3_spec(payload: dict[str, Any]) -> FeatureSpec:
    genes = [str(value) for value in payload["gene_names"]]
    engineered = [str(value) for value in payload["engineered_feature_names"]]
    model_features = [str(value) for value in payload["feature_names"]]
    gene_features = [f"{gene}_log1p_cp10k" for gene in genes]
    if [*gene_features, *engineered] != model_features:
        raise RuntimeError("The gdTAI v3 feature order is internally inconsistent")
    return FeatureSpec(
        gene_names=genes,
        gene_indices=np.arange(len(genes), dtype=np.int32),
        gene_feature_names=gene_features,
        engineered_feature_names=engineered,
        model_feature_names=model_features,
        gene_to_col={gene: index for index, gene in enumerate(genes)},
        engineered_to_col={
            name: len(genes) + index for index, name in enumerate(engineered)
        },
    )


def is_present(values: pd.Series) -> np.ndarray:
    strings = values.astype("string").fillna("").str.strip().str.lower()
    return (~strings.isin(["", "na", "nan", "none", "null"])).to_numpy(dtype=bool)


def add_tcr_metadata_fields(metadata: pd.DataFrame) -> pd.DataFrame:
    has_tra = np.zeros(len(metadata), dtype=bool)
    has_trb = np.zeros(len(metadata), dtype=bool)
    for column in ("cdr3a.aa.seq", "cdr3a.nt.seq", "tra.v", "TRA_cdr3", "TRA_v_gene"):
        if column in metadata:
            has_tra |= is_present(metadata[column])
    for column in ("cdr3b.aa.seq", "cdr3b.nt.seq", "trb.v", "TRB_cdr3", "TRB_v_gene"):
        if column in metadata:
            has_trb |= is_present(metadata[column])
    metadata["has_TRA_evidence"] = has_tra
    metadata["has_TRB_evidence"] = has_trb
    metadata["has_paired_TRA_TRB_evidence"] = has_tra & has_trb
    return metadata


def lineage_high_purity_thresholds(
    lineages: np.ndarray,
    author_clusters: np.ndarray,
    v2_payload: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    config = v2_payload["operating_modes"]["high_purity"]["annotation_thresholds"]
    thresholds = np.full(len(lineages), float(config["other_threshold"]), dtype=np.float32)
    threshold_labels = np.full(len(lineages), "other", dtype=object)
    thresholds[np.asarray(lineages) == "CD4"] = float(config["cd4_threshold"])
    threshold_labels[np.asarray(lineages) == "CD4"] = "CD4_T"
    thresholds[np.asarray(lineages) == "CD8"] = float(config["cd8_threshold"])
    threshold_labels[np.asarray(lineages) == "CD8"] = "CD8_T"
    treg = pd.Series(author_clusters, copy=False).astype(str).str.upper().eq("TREG").to_numpy()
    thresholds[treg] = np.inf
    threshold_labels[treg] = "Treg_disabled"
    return thresholds, threshold_labels


def read_object_payload(
    *,
    payload_dir: Path,
    object_id: str,
    config: dict[str, Any],
    v2_payload: dict[str, Any],
    v3_payload: dict[str, Any],
    v3_spec: FeatureSpec,
    chunk_size: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    prefix = payload_dir / config["prefix"]
    matrix_path = Path(f"{prefix}_model_gene_counts.mtx.gz")
    genes_path = Path(f"{prefix}_model_genes.tsv")
    metadata_path = Path(f"{prefix}_cell_metadata.csv.gz")
    export_manifest_path = Path(f"{prefix}_export_manifest.csv")
    for path in (matrix_path, genes_path, metadata_path, export_manifest_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    genes = genes_path.read_text(encoding="utf-8").splitlines()
    if genes != v3_spec.gene_names:
        raise RuntimeError(f"Model-gene order mismatch in {genes_path}")
    with gzip.open(matrix_path, "rb") as handle:
        gene_by_cell = mmread(handle).tocsr().astype(np.float32)
    metadata = pd.read_csv(metadata_path, low_memory=False)
    if gene_by_cell.shape != (len(v3_spec.gene_names), len(metadata)):
        raise RuntimeError(
            f"Payload shape mismatch for {object_id}: {gene_by_cell.shape} versus "
            f"{len(v3_spec.gene_names)} genes x {len(metadata)} cells"
        )
    if metadata["cell_id"].duplicated().any():
        raise RuntimeError(f"Duplicated cell IDs in {metadata_path}")
    if not metadata["author_lineage"].astype(str).eq(config["lineage"]).all():
        raise RuntimeError(f"Lineage metadata mismatch for {object_id}")
    if not metadata["source_compartment"].astype(str).eq(config["compartment"]).all():
        raise RuntimeError(f"Compartment metadata mismatch for {object_id}")

    totals = pd.to_numeric(metadata["row_sum_counts_layer"], errors="raise").to_numpy(
        dtype=np.float32
    )
    if np.any(~np.isfinite(totals)) or np.any(totals <= 0):
        raise RuntimeError(f"Non-positive library sizes in {metadata_path}")

    v2_base = v2_payload["base_model"]
    v2_genes = [str(gene) for gene in v2_base["gene_names"]]
    v2_indices = np.asarray([v3_spec.gene_to_col[gene] for gene in v2_genes], dtype=np.int32)
    v2_model = v2_base["model_object"]
    v2_high_f1_threshold = float(
        v2_payload["operating_modes"]["high_f1"]["threshold"]
    )
    cluster_column = str(config["cluster_column"])
    if cluster_column not in metadata:
        raise KeyError(f"Published author cluster column is missing: {cluster_column}")
    cluster_numeric = pd.to_numeric(metadata[cluster_column], errors="coerce")
    cluster_id = metadata[cluster_column].astype("string").fillna("unassigned").astype(str)
    integral = cluster_numeric.notna() & (cluster_numeric % 1 == 0)
    cluster_id.loc[integral] = cluster_numeric.loc[integral].astype(int).astype(str)
    author_cluster = cluster_id.map(config["cluster_labels"]).fillna("unassigned")
    v2_high_purity_threshold, v2_high_purity_threshold_annotation = (
        lineage_high_purity_thresholds(
            metadata["author_lineage"].astype(str).to_numpy(dtype=object),
            author_cluster.to_numpy(dtype=object),
            v2_payload,
        )
    )
    v3_model = v3_payload["model_object"]
    v3_threshold = float(v3_payload["threshold"])

    n_cells = len(metadata)
    v2_score = np.empty(n_cells, dtype=np.float32)
    v3_score = np.empty(n_cells, dtype=np.float32)
    quadrants = np.empty(n_cells, dtype=object)
    trdc = np.empty(n_cells, dtype=np.float32)
    trdv = np.empty(n_cells, dtype=np.float32)
    trg = np.empty(n_cells, dtype=np.float32)
    trab = np.empty(n_cells, dtype=np.float32)
    cd3_detected_count = np.empty(n_cells, dtype=np.int8)
    cd3_indices = np.asarray(
        [v3_spec.gene_to_col[gene] for gene in ("CD3D", "CD3E", "CD3G")],
        dtype=np.int32,
    )

    for start in range(0, n_cells, chunk_size):
        end = min(start + chunk_size, n_cells)
        values = gene_by_cell[:, start:end].T.toarray().astype(np.float32, copy=False)
        with np.errstate(divide="ignore", invalid="ignore"):
            x_gene = np.log1p(values * (10_000.0 / totals[start:end, None]))
        x_gene[~np.isfinite(x_gene)] = 0.0
        v2_score[start:end] = v2_model.predict_proba(x_gene[:, v2_indices])[:, 1]
        x_v3 = append_engineered_features(x_gene, v3_spec)
        v3_score[start:end] = v3_model.predict_proba(x_v3)[:, 1]
        quadrants[start:end] = trdc_trdv_quadrant(x_v3, v3_spec)
        trdc[start:end] = x_v3[:, v3_spec.engineered_to_col["TRDC_log1p"]]
        trdv[start:end] = x_v3[:, v3_spec.engineered_to_col["TRDV_score"]]
        trg[start:end] = x_v3[:, v3_spec.engineered_to_col["TRG_score"]]
        trab[start:end] = x_v3[:, v3_spec.engineered_to_col["abT_TCR_score"]]
        cd3_detected_count[start:end] = (x_gene[:, cd3_indices] > 0).sum(axis=1)

    metadata = add_tcr_metadata_fields(metadata)
    metadata["author_cluster_source_column"] = cluster_column
    metadata["author_cluster_id"] = cluster_id
    metadata["author_cluster"] = author_cluster
    metadata["object_id"] = object_id
    metadata["tissue"] = (
        "lung"
        if config["compartment"] == "lung"
        else "lung-associated lymph node"
    )
    metadata["gdtai_v2_score"] = v2_score
    metadata["gdtai_v2_high_f1_threshold"] = v2_high_f1_threshold
    metadata["gdtai_v2_high_f1_predicted"] = v2_score >= v2_high_f1_threshold
    metadata["gdtai_v2_high_purity_threshold"] = v2_high_purity_threshold
    metadata["gdtai_v2_high_purity_threshold_annotation"] = (
        v2_high_purity_threshold_annotation
    )
    metadata["gdtai_v2_high_purity_predicted"] = v2_score >= v2_high_purity_threshold
    metadata["gdtai_v3_score"] = v3_score
    metadata["gdtai_v3_threshold"] = v3_threshold
    metadata["gdtai_v3_predicted"] = v3_score >= v3_threshold
    metadata["tcr_gene_quadrant"] = quadrants
    metadata["TRDC_log1p_cp10k"] = trdc
    metadata["TRDV_sum_log1p_cp10k"] = trdv
    metadata["TRG_sum_log1p_cp10k"] = trg
    metadata["TRAB_sum_log1p_cp10k"] = trab
    metadata["CD3_detected_gene_count"] = cd3_detected_count
    metadata["CD3_any"] = cd3_detected_count >= 1
    metadata["CD3_at_least_2"] = cd3_detected_count >= 2
    metadata["CD3_all_3"] = cd3_detected_count == 3
    metadata["TRDC_positive"] = trdc > 0
    metadata["TRDV_positive"] = trdv > 0
    metadata["TRG_positive"] = trg > 0
    metadata["CD3_any_TRDC_positive"] = metadata["CD3_any"] & metadata["TRDC_positive"]
    metadata["CD3_any_TRDC_or_TRDV_positive"] = metadata["CD3_any"] & (
        metadata["TRDC_positive"] | metadata["TRDV_positive"]
    )

    export_manifest = pd.read_csv(export_manifest_path).iloc[0].to_dict()
    export_manifest.update(
        object_id=object_id,
        payload_matrix=str(matrix_path.relative_to(PROJECT_ROOT)),
        payload_cells=int(n_cells),
        payload_matrix_nnz=int(gene_by_cell.nnz),
        payload_matrix_sha256=sha256_file(matrix_path),
        payload_genes_sha256=sha256_file(genes_path),
        payload_metadata_sha256=sha256_file(metadata_path),
        payload_export_manifest_sha256=sha256_file(export_manifest_path),
        cluster_column=cluster_column,
        cluster_labels=config["cluster_labels"],
    )
    return metadata, export_manifest


def wilson_interval(
    successes: int, total: int, z: float = 1.959963984540054
) -> tuple[float, float]:
    if total == 0:
        return float("nan"), float("nan")
    fraction = successes / total
    denominator = 1 + z * z / total
    center = (fraction + z * z / (2 * total)) / denominator
    margin = z * np.sqrt(
        fraction * (1 - fraction) / total + z * z / (4 * total * total)
    ) / denominator
    return max(0.0, center - margin), min(1.0, center + margin)


def summarize_models(data: pd.DataFrame, group_columns: list[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    groups: Any
    if group_columns:
        groups = data.groupby(group_columns, dropna=False, observed=True, sort=True)
    else:
        groups = [((), data)]
    for key, frame in groups:
        keys = key if isinstance(key, tuple) else (key,)
        if not group_columns:
            keys = ()
        for model_name, (score_column, prediction_column) in MODEL_COLUMNS.items():
            n_cells = len(frame)
            n_predicted = int(frame[prediction_column].sum())
            low, high = wilson_interval(n_predicted, n_cells)
            row = dict(zip(group_columns, keys))
            row.update(
                model=model_name,
                n_cells=n_cells,
                n_predicted=n_predicted,
                predicted_fraction=n_predicted / n_cells if n_cells else np.nan,
                predicted_fraction_ci_low=low,
                predicted_fraction_ci_high=high,
                median_score=float(frame[score_column].median()),
                p95_score=float(frame[score_column].quantile(0.95)),
            )
            rows.append(row)
    return pd.DataFrame(rows)


def make_tables(data: pd.DataFrame, output_dir: Path) -> dict[str, pd.DataFrame]:
    overall = summarize_models(data, [])
    by_object = summarize_models(
        data, ["object_id", "source_compartment", "author_lineage", "tissue"]
    )
    by_donor = summarize_models(
        data, ["object_id", "donor.id.tag"]
        if "donor.id.tag" in data
        else ["object_id"]
    )
    by_cluster = summarize_models(
        data,
        [
            "object_id",
            "author_cluster_source_column",
            "author_cluster_id",
            "author_cluster",
        ],
    )
    by_quadrant = summarize_models(data, ["object_id", "tcr_gene_quadrant"])
    by_tcr_evidence = summarize_models(
        data,
        ["object_id", "has_TRA_evidence", "has_TRB_evidence", "has_paired_TRA_TRB_evidence"],
    )
    by_cite_tag = (
        summarize_models(data, ["object_id", "cite.cell.type.tag"])
        if "cite.cell.type.tag" in data
        else pd.DataFrame()
    )

    scope_rows: list[dict[str, Any]] = []
    for object_id in OBJECTS:
        frame = data.loc[data["object_id"].eq(object_id)]
        deposited = len(frame)
        mapped = int(frame["author_cluster"].ne("unassigned").sum())
        manuscript = MANUSCRIPT_REPORTED_COUNTS[object_id]
        if object_id == "lung_CD4":
            note = "189 deposited cells are in unmapped author cluster IDs 6-8"
        elif object_id.startswith("lymph_node"):
            note = "GEO object is a partial lymph-node deposit"
        else:
            note = "deposited object matches the manuscript count"
        scope_rows.append(
            {
                "object_id": object_id,
                "deposited_cells_analyzed": deposited,
                "mapped_author_cluster_cells": mapped,
                "manuscript_reported_cells": manuscript,
                "deposited_minus_manuscript": deposited - manuscript,
                "deposited_coverage_fraction": deposited / manuscript,
                "scope_note": note,
            }
        )
    scope_audit = pd.DataFrame(scope_rows)

    overlap_columns = [prediction for _, prediction in MODEL_COLUMNS.values()]
    overlap = (
        data.groupby(overlap_columns, observed=True, sort=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    overlap["fraction_of_all_cells"] = overlap["n_cells"] / len(data)
    overlap["prediction_pattern"] = overlap.apply(
        lambda row: "+".join(
            model
            for model, (_, prediction) in MODEL_COLUMNS.items()
            if bool(row[prediction])
        )
        or "none",
        axis=1,
    )
    overlap = overlap.sort_values("n_cells", ascending=False).reset_index(drop=True)

    proxy_definitions = {
        "CD3+TRDC+": data["CD3_any_TRDC_positive"].to_numpy(dtype=bool),
        "CD3+(TRDC_or_TRDV)+": data["CD3_any_TRDC_or_TRDV_positive"].to_numpy(dtype=bool),
        "TRDC+TRDV+": data["tcr_gene_quadrant"].eq("TRDC+TRDV+").to_numpy(),
        "TRDC+TRDV-": data["tcr_gene_quadrant"].eq("TRDC+TRDV-").to_numpy(),
        "TRDC-TRDV+": data["tcr_gene_quadrant"].eq("TRDC-TRDV+").to_numpy(),
    }
    proxy_rows: list[dict[str, Any]] = []
    for proxy_name, mask in proxy_definitions.items():
        for object_id in ["all", *OBJECTS]:
            object_mask = (
                np.ones(len(data), dtype=bool)
                if object_id == "all"
                else data["object_id"].eq(object_id).to_numpy()
            )
            reference = mask & object_mask
            for model_name, (_, prediction_column) in MODEL_COLUMNS.items():
                prediction = data[prediction_column].to_numpy(dtype=bool)
                n_reference = int(reference.sum())
                n_captured = int((reference & prediction).sum())
                proxy_rows.append(
                    {
                        "object_id": object_id,
                        "expression_proxy": proxy_name,
                        "model": model_name,
                        "n_reference_cells": n_reference,
                        "n_model_positive": n_captured,
                        "model_capture_fraction": n_captured / n_reference
                        if n_reference
                        else np.nan,
                    }
                )
    proxy_capture = pd.DataFrame(proxy_rows)

    pairwise_rows: list[dict[str, Any]] = []
    models = list(MODEL_COLUMNS)
    for index, model_a in enumerate(models):
        pred_a = data[MODEL_COLUMNS[model_a][1]].to_numpy(dtype=bool)
        for model_b in models[index + 1 :]:
            pred_b = data[MODEL_COLUMNS[model_b][1]].to_numpy(dtype=bool)
            for object_id in ["all", *OBJECTS]:
                scope = (
                    np.ones(len(data), dtype=bool)
                    if object_id == "all"
                    else data["object_id"].eq(object_id).to_numpy()
                )
                union = scope & (pred_a | pred_b)
                intersection = scope & pred_a & pred_b
                pairwise_rows.append(
                    {
                        "object_id": object_id,
                        "model_a": model_a,
                        "model_b": model_b,
                        "a_positive": int((scope & pred_a).sum()),
                        "b_positive": int((scope & pred_b).sum()),
                        "both_positive": int(intersection.sum()),
                        "a_only": int((scope & pred_a & ~pred_b).sum()),
                        "b_only": int((scope & ~pred_a & pred_b).sum()),
                        "jaccard": int(intersection.sum()) / int(union.sum())
                        if union.any()
                        else np.nan,
                    }
                )
    pairwise = pd.DataFrame(pairwise_rows)

    tables = {
        "overall": overall,
        "by_object": by_object,
        "by_donor": by_donor,
        "by_cluster": by_cluster,
        "by_quadrant": by_quadrant,
        "by_tcr_evidence": by_tcr_evidence,
        "by_cite_tag": by_cite_tag,
        "scope_audit": scope_audit,
        "prediction_overlap": overlap,
        "expression_proxy_capture": proxy_capture,
        "pairwise_overlap": pairwise,
    }
    for name, table in tables.items():
        table.to_csv(output_dir / f"{name}.csv", index=False)
    return tables


def sampled_indices(
    data: pd.DataFrame, prediction_column: str, maximum: int, seed: int
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    positive = np.flatnonzero(data[prediction_column].to_numpy(dtype=bool))
    negative = np.flatnonzero(~data[prediction_column].to_numpy(dtype=bool))
    positive_keep = positive
    if len(positive_keep) > maximum // 2:
        positive_keep = rng.choice(positive_keep, maximum // 2, replace=False)
    remaining = max(0, maximum - len(positive_keep))
    negative_keep = negative
    if len(negative_keep) > remaining:
        negative_keep = rng.choice(negative_keep, remaining, replace=False)
    return np.concatenate([negative_keep, positive_keep])


def plot_umaps(
    data: pd.DataFrame, figure_dir: Path, max_points: int, seed: int
) -> Path:
    output = figure_dir / "all_t_model_prediction_umaps.png"
    fig, axes = plt.subplots(3, 4, figsize=(16, 11.5), constrained_layout=True)
    for row, (model_name, (_, prediction_column)) in enumerate(MODEL_COLUMNS.items()):
        for column, object_id in enumerate(OBJECTS):
            axis = axes[row, column]
            frame = data.loc[data["object_id"].eq(object_id)].reset_index(drop=True)
            if not {"UMAP_1", "UMAP_2"}.issubset(frame.columns):
                axis.text(0.5, 0.5, "UMAP unavailable", ha="center", va="center")
                axis.set_axis_off()
                continue
            indices = sampled_indices(
                frame, prediction_column, max_points, seed + row * 31 + column
            )
            sampled = frame.iloc[indices]
            predicted = sampled[prediction_column].to_numpy(dtype=bool)
            axis.scatter(
                sampled.loc[~predicted, "UMAP_1"],
                sampled.loc[~predicted, "UMAP_2"],
                s=0.6,
                c="#C9CED3",
                alpha=0.25,
                linewidths=0,
                rasterized=True,
            )
            axis.scatter(
                sampled.loc[predicted, "UMAP_1"],
                sampled.loc[predicted, "UMAP_2"],
                s=3.2,
                c="#B33A3A",
                alpha=0.8,
                linewidths=0,
                rasterized=True,
            )
            total_positive = int(frame[prediction_column].sum())
            axis.set_title(
                f"{object_id.replace('_', ' ')}\n{total_positive:,}/{len(frame):,}",
                fontsize=10,
            )
            if column == 0:
                axis.set_ylabel(f"{model_name}\nUMAP 2")
            else:
                axis.set_ylabel("")
            if row == 2:
                axis.set_xlabel("UMAP 1")
            else:
                axis.set_xlabel("")
            axis.set_xticks([])
            axis.set_yticks([])
            axis.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.suptitle("GSE305372: gdTAI predictions in all deposited T-cell objects", fontsize=15)
    fig.savefig(output, dpi=300, facecolor="white")
    plt.close(fig)
    return output


def plot_object_fractions(by_object: pd.DataFrame, figure_dir: Path) -> Path:
    output = figure_dir / "model_call_fraction_by_object.png"
    pivot = by_object.pivot(index="object_id", columns="model", values="predicted_fraction")
    pivot = pivot.reindex(index=list(OBJECTS), columns=list(MODEL_COLUMNS))
    axis = pivot.plot(
        kind="bar",
        figsize=(10.5, 5.5),
        color=["#3F6F8F", "#6B9A79", "#B33A3A"],
        width=0.78,
    )
    axis.set_ylabel("Predicted fraction")
    axis.set_xlabel("")
    axis.set_xticklabels([value.replace("_", " ") for value in pivot.index], rotation=0)
    axis.legend(title="Model", frameon=False)
    axis.spines[["top", "right"]].set_visible(False)
    axis.figure.tight_layout()
    axis.figure.savefig(output, dpi=300, facecolor="white")
    plt.close(axis.figure)
    return output


def plot_overlap(overlap: pd.DataFrame, figure_dir: Path) -> Path:
    output = figure_dir / "three_model_prediction_overlap.png"
    frame = overlap.sort_values("n_cells").copy()
    colors = np.where(frame["prediction_pattern"].eq("none"), "#9AA2A9", "#B85C38")
    fig, axis = plt.subplots(
        figsize=(10, max(4.5, 0.48 * len(frame))), constrained_layout=True
    )
    y = np.arange(len(frame))
    axis.barh(y, frame["n_cells"], color=colors)
    axis.set_yticks(y, frame["prediction_pattern"])
    axis.set_xlabel("Cells")
    axis.set_xscale("log")
    axis.spines[["top", "right"]].set_visible(False)
    for position, value in enumerate(frame["n_cells"]):
        axis.text(value, position, f"  {int(value):,}", va="center", fontsize=8)
    fig.savefig(output, dpi=300, facecolor="white")
    plt.close(fig)
    return output


def plot_score_scatter(data: pd.DataFrame, figure_dir: Path, seed: int) -> Path:
    output = figure_dir / "v2_vs_v3_score_scatter.png"
    rng = np.random.default_rng(seed)
    size = min(180_000, len(data))
    indices = rng.choice(len(data), size=size, replace=False)
    sampled = data.iloc[indices]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), constrained_layout=True)
    for axis, lineage in zip(axes, ("CD4", "CD8")):
        frame = sampled.loc[sampled["author_lineage"].eq(lineage)]
        artist = axis.hexbin(
            frame["gdtai_v2_score"],
            frame["gdtai_v3_score"],
            gridsize=75,
            mincnt=1,
            bins="log",
            cmap="viridis",
        )
        axis.axhline(0.936, color="#B33A3A", linestyle="--", linewidth=1)
        axis.axvline(
            0.97 if lineage == "CD4" else 0.93,
            color="#6B9A79",
            linestyle="--",
            linewidth=1,
        )
        axis.set_title(f"{lineage} objects")
        axis.set_xlabel("gdTAI v2 score")
        axis.set_ylabel("gdTAI v3 score")
        axis.spines[["top", "right"]].set_visible(False)
        fig.colorbar(artist, ax=axis, label="log10(cell count)")
    fig.savefig(output, dpi=300, facecolor="white")
    plt.close(fig)
    return output


def plot_proxy_capture(proxy_capture: pd.DataFrame, figure_dir: Path) -> Path:
    output = figure_dir / "expression_proxy_capture_by_model.png"
    frame = proxy_capture.loc[proxy_capture["object_id"].eq("all")]
    pivot = frame.pivot(
        index="expression_proxy", columns="model", values="model_capture_fraction"
    ).reindex(columns=list(MODEL_COLUMNS))
    axis = pivot.plot(
        kind="bar",
        figsize=(11, 5.8),
        color=["#3F6F8F", "#6B9A79", "#B33A3A"],
        width=0.78,
    )
    axis.set_ylim(0, 1)
    axis.set_ylabel("Model-positive fraction within expression proxy")
    axis.set_xlabel("")
    axis.tick_params(axis="x", rotation=20)
    axis.legend(title="Model", frameon=False)
    axis.spines[["top", "right"]].set_visible(False)
    axis.figure.tight_layout()
    axis.figure.savefig(output, dpi=300, facecolor="white")
    plt.close(axis.figure)
    return output


def format_table(table: pd.DataFrame, columns: list[str], max_rows: int = 60) -> str:
    selected = [column for column in columns if column in table]
    view = table.loc[:, selected].head(max_rows).copy()
    for column in view.columns:
        if "fraction" in column or column in {"median_score", "p95_score", "jaccard"}:
            view[column] = pd.to_numeric(view[column], errors="coerce").map(
                lambda value: "" if pd.isna(value) else f"{value:.4f}"
            )
    return view.fillna("").to_html(
        index=False, border=0, classes="data-table", escape=True
    )


def render_report(
    *,
    data: pd.DataFrame,
    tables: dict[str, pd.DataFrame],
    figures: list[Path],
    output_dir: Path,
    report_dir: Path,
    manifest: dict[str, Any],
) -> Path:
    report_dir.mkdir(parents=True, exist_ok=True)
    report_figures = report_dir / "figures"
    report_figures.mkdir(parents=True, exist_ok=True)
    for source in figures:
        shutil.copy2(source, report_figures / source.name)

    counts = tables["overall"].set_index("model")
    v2_f1 = int(counts.loc["v2_high_f1", "n_predicted"])
    v2_purity = int(counts.loc["v2_high_purity", "n_predicted"])
    v3 = int(counts.loc["v3_round14", "n_predicted"])
    cd3_trdc = int(data["CD3_any_TRDC_positive"].sum())
    broad_proxy = int(data["CD3_any_TRDC_or_TRDV_positive"].sum())
    paired = int(data["has_paired_TRA_TRB_evidence"].sum())
    paired_rows = []
    for model, (_, prediction_column) in MODEL_COLUMNS.items():
        paired_predictions = int(
            (data[prediction_column] & data["has_paired_TRA_TRB_evidence"]).sum()
        )
        paired_rows.append(
            {
                "model": model,
                "paired_TRA_TRB_cells": paired,
                "model_positive_paired_TRA_TRB": paired_predictions,
                "conflict_screening_fraction": paired_predictions / paired if paired else np.nan,
            }
        )
    paired_table = pd.DataFrame(paired_rows)
    paired_table.to_csv(output_dir / "paired_trab_conflict_screening.csv", index=False)

    report = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GSE305372 all-T gdTAI v2/v3 comparison</title>
<style>
:root {{ --ink:#20262c; --navy:#17324d; --line:#d6dee5; --soft:#f3f6f8; --red:#b33a3a; --green:#527b61; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; font-family:Arial,Helvetica,sans-serif; color:var(--ink); line-height:1.48; background:white; }}
header {{ background:var(--navy); color:white; padding:30px max(5vw,28px); }}
header h1 {{ margin:0 0 7px; font-size:30px; letter-spacing:0; }}
main {{ max-width:1200px; margin:0 auto; padding:28px; }}
h2 {{ margin-top:34px; border-bottom:2px solid var(--line); padding-bottom:7px; }}
h3 {{ margin-top:25px; }}
.metrics {{ display:grid; grid-template-columns:repeat(4,minmax(155px,1fr)); gap:12px; margin:20px 0; }}
.metric {{ border:1px solid var(--line); border-radius:6px; padding:14px; min-height:86px; }}
.metric strong {{ display:block; color:var(--navy); font-size:24px; }}
.note {{ border-left:4px solid #d17a22; background:#fff7ed; padding:12px 16px; margin:18px 0; }}
.method {{ border-left:4px solid var(--green); background:#f1f7f3; padding:12px 16px; margin:18px 0; }}
figure {{ margin:24px 0; break-inside:avoid; }}
figure img {{ max-width:100%; height:auto; border:1px solid var(--line); }}
.data-table {{ border-collapse:collapse; width:100%; font-size:12px; display:block; overflow-x:auto; }}
.data-table th,.data-table td {{ border:1px solid var(--line); padding:6px 8px; text-align:left; white-space:nowrap; }}
.data-table th {{ background:var(--soft); }}
code {{ background:var(--soft); padding:2px 4px; }}
@media(max-width:800px) {{ .metrics {{ grid-template-columns:1fr 1fr; }} }}
@page {{ size:A4 landscape; margin:10mm; }}
@media print {{
 body {{ font-size:8.5pt; }} main {{ max-width:none; padding:0; }} header {{ padding:8mm; }}
 h2 {{ break-after:avoid; }} .metrics {{ grid-template-columns:repeat(4,1fr); gap:5px; }}
 .metric {{ min-height:0; padding:7px; }} .metric strong {{ font-size:15pt; }}
 figure {{ break-inside:avoid; page-break-inside:avoid; }}
 figure img {{ display:block; width:auto; max-width:100%; max-height:165mm; margin:0 auto; object-fit:contain; }}
 .data-table {{ display:table; width:100%; font-size:6.5pt; overflow:visible; }}
 .data-table thead {{ display:table-header-group; }} .data-table tr {{ break-inside:avoid; }}
 .data-table th,.data-table td {{ padding:2px 3px; white-space:normal; overflow-wrap:anywhere; }}
}}
</style>
</head>
<body>
<header>
 <h1>gdTAI v2 and v3 in GSE305372</h1>
 <div>Corrected inference on every cell in all four deposited author-finalized CD4/CD8 T-cell objects</div>
</header>
<main>
<div class="metrics">
 <div class="metric"><strong>{len(data):,}</strong>transcriptome-eligible deposited T cells</div>
 <div class="metric"><strong>{v2_f1:,}</strong>V2 high-F1 calls</div>
 <div class="metric"><strong>{v2_purity:,}</strong>V2 high-purity calls</div>
 <div class="metric"><strong>{v3:,}</strong>V3 Round 14 calls</div>
 <div class="metric"><strong>{cd3_trdc:,}</strong>CD3+TRDC+ expression proxy</div>
 <div class="metric"><strong>{broad_proxy:,}</strong>CD3+(TRDC or TRDV)+ proxy</div>
 <div class="metric"><strong>{paired:,}</strong>paired TRA/TRB evidence</div>
 <div class="metric"><strong>0</strong>CITE-based exclusions</div>
</div>

<div class="method"><strong>Corrected inclusion rule.</strong> Every cell in the authors' processed lung CD4, lung CD8, lymph-node CD4, and lymph-node CD8 Seurat objects was analyzed. <code>cite.cell.type.tag</code> is retained only as audit metadata and never used to select cells.</div>
<div class="note"><strong>Inference, not ground-truth evaluation.</strong> The processed objects contain TRA/TRB metadata but no matched TRG/TRD clonotype calls. Consequently, model-positive paired-TRA/TRB cells are conflict-screening candidates rather than proven false positives, and expression-proxy capture is not biological recall.</div>

<h2>Denominator audit</h2>
<p>The deposited lung CD4 object has 370,422 cells, whereas the manuscript reports 370,233. The difference is exactly 189 cells in cluster IDs 6-8, which are not present in the authors' finalized phenotype map; no model called any of those 189 cells. They remain included here under the all-deposited-cells rule. Both lymph-node objects are partial deposits relative to the manuscript cohort.</p>
{format_table(tables['scope_audit'], ['object_id','deposited_cells_analyzed','mapped_author_cluster_cells','manuscript_reported_cells','deposited_minus_manuscript','deposited_coverage_fraction','scope_note'])}
<p>The separate DICE-TCR blood CD4 CSV contains {manifest['dice_tcr_blood_cd4_resource'].get('data_rows', 0):,} TCR/metadata rows and {manifest['dice_tcr_blood_cd4_resource'].get('n_columns', 0)} columns, but no raw transcriptome matrix or whole-library count denominator. It is therefore not eligible for gdTAI expression inference and is not part of the 690,715-cell denominator.</p>

<h2>Models and normalization</h2>
<p>Both releases were applied to the same raw RNA counts. Each model gene was transformed as <code>log1p(raw count * 10,000 / whole-transcriptome raw-count total)</code>. V2 high-F1 uses its fixed threshold <code>{manifest['v2_high_f1_threshold']:.6f}</code>. V2 high-purity uses the released annotation-specific policy: <code>0.97</code> for CD4, <code>0.93</code> for CD8, and no calls from author-defined Treg clusters. V3 is the checksum-pinned promoted Round 14 model at <code>{manifest['v3_threshold']:.3f}</code>.</p>
<p>V2 SHA256: <code>{manifest['v2_model_sha256']}</code><br>V3 SHA256: <code>{manifest['v3_model_sha256']}</code></p>

<h2>Overall and object-level calls</h2>
{format_table(tables['overall'], ['model','n_cells','n_predicted','predicted_fraction','predicted_fraction_ci_low','predicted_fraction_ci_high','median_score','p95_score'])}
{format_table(tables['by_object'], ['object_id','source_compartment','author_lineage','model','n_cells','n_predicted','predicted_fraction','predicted_fraction_ci_low','predicted_fraction_ci_high'])}
<figure><img src="figures/model_call_fraction_by_object.png" alt="Model-positive fractions by deposited T-cell object"></figure>
<figure><img src="figures/all_t_model_prediction_umaps.png" alt="UMAPs showing model predictions in all deposited T-cell objects"></figure>

<h3>Author phenotype clusters</h3>
<p>Cluster resolutions and labels follow the authors' public figure script: <code>RNA_snn_res.0.2</code> for lung and <code>RNA_snn_res.0.3</code> for lymph node. These labels are also the source for the V2 high-purity Treg exclusion.</p>
{format_table(tables['by_cluster'], ['object_id','author_cluster_id','author_cluster','model','n_cells','n_predicted','predicted_fraction'], 120)}

<h2>Model agreement</h2>
{format_table(tables['pairwise_overlap'].loc[tables['pairwise_overlap']['object_id'].eq('all')], ['model_a','model_b','a_positive','b_positive','both_positive','a_only','b_only','jaccard'])}
{format_table(tables['prediction_overlap'], ['prediction_pattern','n_cells','fraction_of_all_cells'])}
<figure><img src="figures/three_model_prediction_overlap.png" alt="Prediction overlap among V2 and V3"></figure>
<figure><img src="figures/v2_vs_v3_score_scatter.png" alt="V2 versus V3 score scatter"></figure>

<h2>TCR-expression audit</h2>
<p><code>CD3+</code> means at least one raw UMI in <code>CD3D</code>, <code>CD3E</code>, or <code>CD3G</code>. Gene-expression categories are dropout-sensitive proxies and are not used as ground truth.</p>
{format_table(tables['expression_proxy_capture'].loc[tables['expression_proxy_capture']['object_id'].eq('all')], ['expression_proxy','model','n_reference_cells','n_model_positive','model_capture_fraction'])}
<figure><img src="figures/expression_proxy_capture_by_model.png" alt="Expression proxy capture by model"></figure>
{format_table(tables['by_quadrant'], ['object_id','tcr_gene_quadrant','model','n_cells','n_predicted','predicted_fraction'], 80)}

<h2>TRA/TRB conflict screening</h2>
{format_table(paired_table, ['model','paired_TRA_TRB_cells','model_positive_paired_TRA_TRB','conflict_screening_fraction'])}

<h2>CITE-tag audit</h2>
<p>CITE tags are shown only to expose which cells were previously lost. They do not determine inclusion or prediction thresholds in this corrected analysis.</p>
{format_table(tables['by_cite_tag'], ['object_id','cite.cell.type.tag','model','n_cells','n_predicted','predicted_fraction'], 100)}

<h2>Scope limitation</h2>
<p>The GEO deposit contains the complete manuscript lung objects but only a subset of the manuscript lymph-node cohort. This report therefore describes all <em>deposited</em> T cells, not every lymph-node cell reported in the article.</p>
<p>Source: <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305372">GEO GSE305372</a>. Study: <a href="https://www.nature.com/articles/s41590-026-02592-6">Fajardo-Rosas et al., Nature Immunology (2026)</a>.</p>
<p>Workflow: <code>workflows/intake/export_gse305372_all_t_model_payload.R</code> and <code>workflows/gdtai/run_gdtai_gse305372_all_t.py</code>.</p>
</main>
</body>
</html>
"""
    output = report_dir / "index.html"
    output.write_text(report, encoding="utf-8")
    return output


def prediction_export_columns(data: pd.DataFrame) -> list[str]:
    preferred = [
        "cell_id",
        "cell_id_original",
        "object_id",
        "source_compartment",
        "tissue",
        "author_lineage",
        "donor.id.tag",
        "sample.id.tag",
        "library.id.tag",
        "cite.cell.type.tag",
        "author_cluster_source_column",
        "author_cluster_id",
        "author_cluster",
        "UMAP_1",
        "UMAP_2",
        "has_TRA_evidence",
        "has_TRB_evidence",
        "has_paired_TRA_TRB_evidence",
        "CD3_detected_gene_count",
        "CD3_any",
        "CD3_at_least_2",
        "CD3_all_3",
        "TRDC_positive",
        "TRDV_positive",
        "TRG_positive",
        "CD3_any_TRDC_positive",
        "CD3_any_TRDC_or_TRDV_positive",
        "tcr_gene_quadrant",
        "TRDC_log1p_cp10k",
        "TRDV_sum_log1p_cp10k",
        "TRG_sum_log1p_cp10k",
        "TRAB_sum_log1p_cp10k",
        "gdtai_v2_score",
        "gdtai_v2_high_f1_threshold",
        "gdtai_v2_high_f1_predicted",
        "gdtai_v2_high_purity_threshold",
        "gdtai_v2_high_purity_threshold_annotation",
        "gdtai_v2_high_purity_predicted",
        "gdtai_v3_score",
        "gdtai_v3_threshold",
        "gdtai_v3_predicted",
        "row_sum_counts_layer",
        "n_detected_genes_counts_layer",
    ]
    return [column for column in preferred if column in data]


def main() -> None:
    args = parse_args()
    for directory in (args.output_dir, args.figure_dir, args.log_dir, args.report_dir):
        directory.mkdir(parents=True, exist_ok=True)
    prediction_path = args.output_dir / "all_t_v2_v3_predictions.csv.gz"
    if prediction_path.exists() and not args.overwrite:
        raise FileExistsError(f"Output exists: {prediction_path}; pass --overwrite")

    registry = pd.read_csv(MODEL_REGISTRY)
    expected_v2 = registry.loc[registry["model_id"].eq("gdtai_v2")].iloc[0]
    expected_v3 = registry.loc[
        registry["model_id"].eq("gdtai_v3") & registry["status"].eq("promoted_default")
    ].iloc[0]
    v2_sha = sha256_file(args.v2_model_pkl)
    v3_sha = sha256_file(args.v3_model_pkl)
    if v2_sha != str(expected_v2["sha256"]):
        raise RuntimeError("gdTAI v2 checksum does not match the model registry")
    if v3_sha != str(expected_v3["sha256"]):
        raise RuntimeError("gdTAI v3 checksum does not match the model registry")

    with args.v2_model_pkl.open("rb") as handle:
        v2_payload = pickle.load(handle)
    with args.v3_model_pkl.open("rb") as handle:
        v3_payload = pickle.load(handle)
    v3_spec = build_v3_spec(v3_payload)
    v2_genes = [str(gene) for gene in v2_payload["base_model"]["gene_names"]]
    if not set(v2_genes).issubset(v3_spec.gene_names):
        raise RuntimeError("The shared V3 payload does not contain every V2 model gene")
    dice_resource = inspect_dice_tcr_resource(DICE_TCR_RESOURCE, v3_spec.gene_names)

    frames: list[pd.DataFrame] = []
    source_manifests: list[dict[str, Any]] = []
    for object_id, config in OBJECTS.items():
        frame, source_manifest = read_object_payload(
            payload_dir=args.payload_dir,
            object_id=object_id,
            config=config,
            v2_payload=v2_payload,
            v3_payload=v3_payload,
            v3_spec=v3_spec,
            chunk_size=args.chunk_size,
        )
        frames.append(frame)
        source_manifests.append(source_manifest)
    data = pd.concat(frames, ignore_index=True, sort=False)

    v2_high_f1_threshold = float(
        v2_payload["operating_modes"]["high_f1"]["threshold"]
    )
    v3_threshold = float(v3_payload["threshold"])
    validation_checks = {
        "combined_cell_ids_unique": bool(~data["cell_id"].duplicated().any()),
        "four_author_objects_present": set(data["object_id"].unique()) == set(OBJECTS),
        "all_export_manifests_are_unfiltered": all(
            item.get("inclusion_rule") == "all_cells_in_author_finalized_T_cell_object"
            for item in source_manifests
        ),
        "combined_cells_equal_exported_cells": len(data)
        == sum(int(item["exported_cells"]) for item in source_manifests),
        "v2_scores_finite_and_bounded": bool(
            np.isfinite(data["gdtai_v2_score"]).all()
            and data["gdtai_v2_score"].between(0, 1, inclusive="both").all()
        ),
        "v3_scores_finite_and_bounded": bool(
            np.isfinite(data["gdtai_v3_score"]).all()
            and data["gdtai_v3_score"].between(0, 1, inclusive="both").all()
        ),
        "v2_high_f1_mask_exact": bool(
            np.array_equal(
                data["gdtai_v2_high_f1_predicted"].to_numpy(dtype=bool),
                data["gdtai_v2_score"].to_numpy() >= v2_high_f1_threshold,
            )
        ),
        "v2_high_purity_mask_exact": bool(
            np.array_equal(
                data["gdtai_v2_high_purity_predicted"].to_numpy(dtype=bool),
                data["gdtai_v2_score"].to_numpy()
                >= data["gdtai_v2_high_purity_threshold"].to_numpy(),
            )
        ),
        "v2_high_purity_treg_disabled": bool(
            np.isinf(
                data.loc[
                    data["author_cluster"].eq("TREG"),
                    "gdtai_v2_high_purity_threshold",
                ].to_numpy()
            ).all()
            and not data.loc[
                data["author_cluster"].eq("TREG"),
                "gdtai_v2_high_purity_predicted",
            ].any()
        ),
        "v3_mask_exact": bool(
            np.array_equal(
                data["gdtai_v3_predicted"].to_numpy(dtype=bool),
                data["gdtai_v3_score"].to_numpy() >= v3_threshold,
            )
        ),
        "v2_gene_set_is_v3_subset": set(v2_genes).issubset(v3_spec.gene_names),
        "v3_gene_count_is_210": len(v3_spec.gene_names) == 210,
        "dice_tcr_resource_correctly_ineligible": bool(
            dice_resource.get("available")
            and not dice_resource.get("eligible_for_gdtai")
            and not dice_resource.get("model_genes_present_as_columns")
        ),
        "lung_cd4_mapped_clusters_match_manuscript_count": int(
            (
                data["object_id"].eq("lung_CD4")
                & data["author_cluster"].ne("unassigned")
            ).sum()
        )
        == MANUSCRIPT_REPORTED_COUNTS["lung_CD4"],
        "lung_cd8_deposit_matches_manuscript_count": int(
            data["object_id"].eq("lung_CD8").sum()
        )
        == MANUSCRIPT_REPORTED_COUNTS["lung_CD8"],
    }
    failed = [name for name, passed in validation_checks.items() if not passed]
    if failed:
        raise RuntimeError(f"Validation failed: {', '.join(failed)}")

    export_columns = prediction_export_columns(data)
    data.loc[:, export_columns].to_csv(
        prediction_path, index=False, compression="gzip"
    )
    any_positive = np.zeros(len(data), dtype=bool)
    for _, prediction_column in MODEL_COLUMNS.values():
        any_positive |= data[prediction_column].to_numpy(dtype=bool)
    data.loc[any_positive, export_columns].to_csv(
        args.output_dir / "all_t_any_model_positive_cells.csv.gz",
        index=False,
        compression="gzip",
    )

    tables = make_tables(data, args.output_dir)
    figures = [
        plot_umaps(data, args.figure_dir, args.max_plot_points, args.seed),
        plot_object_fractions(tables["by_object"], args.figure_dir),
        plot_overlap(tables["prediction_overlap"], args.figure_dir),
        plot_score_scatter(data, args.figure_dir, args.seed),
        plot_proxy_capture(tables["expression_proxy_capture"], args.figure_dir),
    ]

    manifest = {
        "dataset_id": DATASET_ID,
        "analysis": "all_deposited_author_finalized_T_cells_gdtai_v2_v3_inference",
        "inclusion_rule": "all cells in four deposited CD4/CD8 Seurat objects",
        "cite_filter_used": False,
        "normalization": "log1p(raw_count * 10000 / whole_transcriptome_total_count)",
        "total_cells": int(len(data)),
        "v2_model_path": str(args.v2_model_pkl.relative_to(PROJECT_ROOT)),
        "v2_model_sha256": v2_sha,
        "v2_high_f1_threshold": v2_high_f1_threshold,
        "v2_high_purity_thresholds": {
            "CD4": 0.97,
            "CD8": 0.93,
            "Treg": "disabled",
        },
        "v3_model_path": str(args.v3_model_pkl.relative_to(PROJECT_ROOT)),
        "v3_model_sha256": v3_sha,
        "v3_threshold": v3_threshold,
        "prediction_counts": {
            model: int(data[prediction_column].sum())
            for model, (_, prediction_column) in MODEL_COLUMNS.items()
        },
        "expression_proxy_counts": {
            "CD3_any_TRDC_positive": int(data["CD3_any_TRDC_positive"].sum()),
            "CD3_any_TRDC_or_TRDV_positive": int(
                data["CD3_any_TRDC_or_TRDV_positive"].sum()
            ),
        },
        "source_payloads": source_manifests,
        "dice_tcr_blood_cd4_resource": dice_resource,
        "validation_checks": validation_checks,
        "ground_truth_available": False,
        "tcr_metadata_limit": (
            "TRA/TRB fields are available; matched TRG/TRD clonotype fields are absent"
        ),
    }
    manifest_path = args.log_dir / "gdtai_v2_v3_gse305372_all_t_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (args.log_dir / "validation_checks.json").write_text(
        json.dumps(validation_checks, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    report = render_report(
        data=data,
        tables=tables,
        figures=figures,
        output_dir=args.output_dir,
        report_dir=args.report_dir,
        manifest=manifest,
    )
    summary_lines = [
        "# GSE305372 all-T gdTAI V2/V3 summary",
        "",
        f"- cells: {len(data):,}",
        "- inclusion: every cell in all four deposited CD4/CD8 Seurat objects",
        "- CITE filter: none",
        f"- V2 high-F1: {int(data['gdtai_v2_high_f1_predicted'].sum()):,}",
        f"- V2 high-purity: {int(data['gdtai_v2_high_purity_predicted'].sum()):,}",
        f"- V3 Round 14: {int(data['gdtai_v3_predicted'].sum()):,}",
        f"- CD3+TRDC+ expression proxy: {int(data['CD3_any_TRDC_positive'].sum()):,}",
        "- evaluation status: inference only; matched gamma-delta TCR ground truth unavailable",
        f"- report: `{report.relative_to(PROJECT_ROOT)}`",
    ]
    (args.log_dir / "gdtai_v2_v3_gse305372_all_t_summary.md").write_text(
        "\n".join(summary_lines) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "cells": len(data),
                "prediction_counts": manifest["prediction_counts"],
                "expression_proxy_counts": manifest["expression_proxy_counts"],
                "report": str(report),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
