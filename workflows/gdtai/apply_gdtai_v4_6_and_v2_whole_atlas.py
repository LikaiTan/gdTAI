#!/usr/bin/env python3
"""Apply frozen gdTAI V4.6 and V2 modes to the corrected 5.93M-cell atlas."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import time
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
for value in (ROOT, ROOT / "src", ROOT / "workflows/gdtai"):
    if str(value) not in sys.path:
        sys.path.insert(0, str(value))

import compare_frozen_gdtai_profiles as v2_eval
from evaluate_gdtai_v4_4_reused_lockbox import calibration_from, load_model, predict_booster
from run_gdt_prediction_package_evaluation import read_obs_column, read_string_dataset
from train_gdtai_v4_2_nested import exclusion_flags
from train_gdtai_v4_6_development import architecture_composer


DEFAULT_CONFIG = ROOT / "configs/models/gdtai/v4_6_whole_atlas_inference.json"
MODEL_KEYS = ("v4_6_highest_f1", "v4_6_high_purity", "v2_high_f1", "v2_high_purity")
FLAG_KEYS = (
    "likely_nk", "author_nk", "conservative_expression_nk", "shared_cytotoxic_ambiguous",
    "paired_ab_without_gd", "productive_gd", "trdc_positive", "trdv_positive", "trdc_plus_trdv_minus",
)
OBS_COLUMNS = (
    "cell_id", "source_gse_id", "source_accession_harmonized_v2", "sample_id_harmonized_v2",
    "donor_id_harmonized_v2", "tissue_harmonized_v2", "specimen_context_harmonized_v2",
    "tumor_type_harmonized_v2", "source_cell_type", "source_cell_type_level1",
    "has_TRA_TRB_paired", "has_any_ab_tcr", "has_any_gd_tcr", "has_TRG_TRD_paired",
)


def resolve(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def as_bool(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values)
    if array.dtype.kind == "b":
        return array.astype(bool, copy=False)
    return pd.Series(array, copy=False).astype("string").fillna("").str.strip().str.lower().isin(
        {"true", "1", "yes", "y", "t"}
    ).to_numpy(bool)


def read_csr_chunk(matrix: h5py.Group, start: int, stop: int, n_vars: int) -> sparse.csr_matrix:
    absolute = np.asarray(matrix["indptr"][start:stop + 1], dtype=np.int64)
    first, last = int(absolute[0]), int(absolute[-1])
    return sparse.csr_matrix(
        (
            np.asarray(matrix["data"][first:last]),
            np.asarray(matrix["indices"][first:last]),
            absolute - first,
        ),
        shape=(stop - start, n_vars),
    )


def author_nk_mask(primary: np.ndarray, secondary: np.ndarray, exclusions: list[str]) -> np.ndarray:
    first = pd.Series(primary, copy=False).astype("string").fillna("").str.strip()
    second = pd.Series(secondary, copy=False).astype("string").fillna("").str.strip()
    text = first.where(first.ne(""), second).str.upper()
    result = text.str.contains(r"(?:^|[^A-Z])NK(?:[^A-Z]|$)", regex=True).to_numpy(bool)
    for token in exclusions:
        result &= ~text.str.contains(str(token).upper(), regex=False).to_numpy(bool)
    return result


def feature_values(matrix: np.ndarray, lookup: dict[str, int], genes: list[str]) -> np.ndarray:
    return np.column_stack([matrix[:, lookup[gene]] for gene in genes]).astype(np.float32, copy=False)


def nk_audit_flags(x: np.ndarray, lookup: dict[str, int], config: dict, author_nk: np.ndarray) -> dict[str, np.ndarray]:
    audit = config["nk_audit"]
    tier_a = feature_values(x, lookup, audit["tier_a_genes"]) > 0
    tier_b = feature_values(x, lookup, audit["tier_b_genes"]) > 0
    shared = feature_values(x, lookup, audit["shared_cytotoxic_genes"]) > 0
    cd3 = feature_values(x, lookup, audit["t_lineage_genes"]) > 0
    trdv = feature_values(x, lookup, audit["dedicated_delta_genes"]) > 0
    tier_a_count = tier_a.sum(axis=1)
    innate_count = tier_a_count + tier_b.sum(axis=1)
    trdv_positive = trdv.any(axis=1)
    conservative = (
        ~trdv_positive
        & (cd3.sum(axis=1) <= 1)
        & (tier_a_count >= 1)
        & (innate_count >= 2)
    )
    likely = author_nk | conservative
    return {
        "likely_nk": likely,
        "author_nk": author_nk,
        "conservative_expression_nk": conservative,
        "shared_cytotoxic_ambiguous": (shared.sum(axis=1) >= 3) & ~likely,
        "trdv_positive": trdv_positive,
        "tier_a_detected_count": tier_a_count.astype(np.int8),
        "innate_detected_count": innate_count.astype(np.int8),
        "cd3_detected_count": cd3.sum(axis=1).astype(np.int8),
        "shared_cytotoxic_detected_count": shared.sum(axis=1).astype(np.int8),
    }


def prepare_models(config: dict) -> dict[str, Any]:
    model_dir = resolve(config["v4_6_model_dir"])
    contract_path = model_dir / "model_contract.json"
    if sha256_file(contract_path) != config["v4_6_contract_sha256"]:
        raise RuntimeError("Frozen V4.6 contract checksum mismatch")
    contract = json.loads(contract_path.read_text())
    if contract.get("selected_candidate") != "conditional_gamma" or contract.get("model_promoted"):
        raise RuntimeError("Unexpected V4.6 contract state")
    calibration = json.loads((model_dir / "platt_calibration.json").read_text())
    v2_config = json.loads((ROOT / "configs/models/gdtai/extension_evaluation.json").read_text())
    profiles = [profile for profile in v2_eval.load_profiles(v2_config) if profile.profile_id.startswith("v2_")]
    if {profile.profile_id for profile in profiles} != {"v2_high_f1", "v2_high_purity"}:
        raise RuntimeError("Frozen V2 profile set is incomplete")
    if sha256_file(resolve(config["v2_model"])) != config["v2_model_sha256"]:
        raise RuntimeError("Frozen V2 model checksum mismatch")
    return {
        "contract": contract,
        "stage1_model": load_model(model_dir / "stage1_t_lineage_gate.ubj"),
        "stage2_model": load_model(model_dir / "stage2_gdt_classifier.ubj"),
        "stage1_calibration": calibration_from(calibration["stage1"]),
        "stage2_calibration": calibration_from(calibration["stage2"]),
        "stage2_composer": architecture_composer("conditional_gamma"),
        "v2_profiles": profiles,
    }


def model_feature_union(models: dict, config: dict) -> list[str]:
    contract = models["contract"]
    genes: list[str] = []
    sources = [contract["stage1_feature_names"], contract["stage2_base_feature_names"]]
    v2_genes, _ = v2_eval.feature_schema(models["v2_profiles"])
    sources.append(v2_genes)
    audit = config["nk_audit"]
    sources.extend([
        audit["tier_a_genes"], audit["tier_b_genes"], audit["shared_cytotoxic_genes"],
        audit["t_lineage_genes"], audit["dedicated_delta_genes"],
        ["TRDC", "TRAC", "TRBC1", "TRBC2"],
    ])
    for source in sources:
        for gene in source:
            if gene not in genes:
                genes.append(str(gene))
    return genes


def preflight(config_path: Path, config: dict, models: dict) -> dict[str, Any]:
    input_path = resolve(config["input_h5ad"])
    manifest = json.loads(resolve(config["input_manifest"]).read_text())
    stat = input_path.stat()
    if Path(manifest["candidate"]).resolve() != input_path.resolve():
        raise RuntimeError("Input path differs from the validated TCR-corrected candidate")
    if manifest["candidate_sha256"] != config["expected_input_sha256"]:
        raise RuntimeError("Configured input checksum differs from its validation manifest")
    if int(manifest["candidate_size_bytes"]) != stat.st_size:
        raise RuntimeError("Input H5AD size differs from its validation manifest")
    genes = model_feature_union(models, config)
    with h5py.File(input_path, "r") as handle:
        matrix = handle["X"]
        shape = tuple(int(value) for value in matrix.attrs["shape"])
        if shape[0] != int(config["expected_cells"]) or matrix.attrs.get("encoding-type") != "csr_matrix":
            raise RuntimeError(f"Unexpected input matrix contract: {shape}")
        var_names = read_string_dataset(handle["var"]["gene"]).astype(str)
        missing = sorted(set(genes) - set(var_names))
        if missing:
            raise RuntimeError(f"Whole atlas lacks required inference/audit genes: {missing}")
        missing_obs = sorted(set(OBS_COLUMNS) - set(handle["obs"].keys()))
        if missing_obs:
            raise RuntimeError(f"Whole atlas lacks required metadata columns: {missing_obs}")
        source = read_obs_column(handle, "source_gse_id").astype(str)
        large = source == config["large_luad_lusc_source"]
        if int(large.sum()) != 759436:
            raise RuntimeError("GSE243013 cell count changed")
    output = {
        "status": "PASS_V4_6_V2_WHOLE_ATLAS_PREFLIGHT",
        "protocol_version": config["protocol_version"],
        "config_sha256": sha256_file(config_path),
        "input_path": str(input_path),
        "input_manifest_sha256": config["expected_input_sha256"],
        "input_size_bytes": stat.st_size,
        "input_mtime_ns": stat.st_mtime_ns,
        "n_cells": shape[0],
        "n_genes": shape[1],
        "n_required_union_genes": len(genes),
        "large_luad_lusc_source": config["large_luad_lusc_source"],
        "large_luad_lusc_cells": int(large.sum()),
        "models_or_thresholds_retuned": False,
        "h5ad_mutated": False,
    }
    log_dir = resolve(config["output_log_dir"])
    log_dir.mkdir(parents=True, exist_ok=True)
    (log_dir / "preflight_summary.json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    return output


def smoke_test(config: dict, models: dict, n_cells: int = 2000) -> dict[str, Any]:
    input_path = resolve(config["input_h5ad"])
    union_genes = model_feature_union(models, config)
    with h5py.File(input_path, "r") as handle:
        matrix = handle["X"]
        _, n_vars = (int(value) for value in matrix.attrs["shape"])
        var_names = read_string_dataset(handle["var"]["gene"]).astype(str)
        var_lookup = {gene: index for index, gene in enumerate(var_names)}
        lookup = {gene: index for index, gene in enumerate(union_genes)}
        counts = read_csr_chunk(matrix, 0, n_cells, n_vars)
        totals = np.asarray(counts.sum(axis=1), dtype=np.float64).ravel()
        x = counts[:, np.asarray([var_lookup[gene] for gene in union_genes])].toarray().astype(np.float32)
        x *= (float(config["normalization_target_sum"]) / totals).astype(np.float32)[:, None]
        np.log1p(x, out=x)
        v46_score, v46_hf1, v46_hp, _ = v4_scores(x, lookup, models)
        v2_score, v2_calls, annotation = v2_scores(x, union_genes, models)
    for name, values in (("V4.6", v46_score), ("V2", v2_score)):
        if not np.isfinite(values).all():
            raise RuntimeError(f"{name} smoke-test scores contain non-finite values")
    return {
        "status": "PASS_V4_6_V2_SMOKE_TEST",
        "n_cells": n_cells,
        "v4_6_highest_f1": int(v46_hf1.sum()),
        "v4_6_high_purity": int(v46_hp.sum()),
        "v2_high_f1": int(v2_calls["v2_high_f1"].sum()),
        "v2_high_purity": int(v2_calls["v2_high_purity"].sum()),
        "v2_expression_annotations": pd.Series(annotation).value_counts().to_dict(),
    }


def v4_scores(x: np.ndarray, lookup: dict[str, int], models: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    contract = models["contract"]
    stage1_names = contract["stage1_feature_names"]
    stage2_names = contract["stage2_base_feature_names"]
    stage1 = models["stage1_calibration"].apply(
        predict_booster(models["stage1_model"], feature_values(x, lookup, stage1_names))
    )
    stage2_base = feature_values(x, lookup, stage2_names)
    stage2_values, effective_names = models["stage2_composer"](stage2_base, stage2_names)
    if effective_names != contract["stage2_effective_feature_names"]:
        raise RuntimeError("V4.6 effective feature order changed during atlas inference")
    stage2 = models["stage2_calibration"].apply(
        predict_booster(models["stage2_model"], stage2_values)
    )
    excluded = exclusion_flags(x, list(lookup))[2]
    effective = stage2.copy()
    effective[(stage1 < float(contract["stage1_threshold"])) | excluded] = -1.0
    high_f1 = effective >= float(contract["operating_modes"]["highest_f1"]["threshold"])
    high_purity = effective >= float(contract["operating_modes"]["high_purity"]["threshold"])
    return effective, high_f1, high_purity, excluded


def v2_scores(x: np.ndarray, genes: list[str], models: dict) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray]:
    spec = v2_eval.FeatureSpec(
        gene_names=genes,
        gene_indices=np.arange(len(genes), dtype=np.int32),
        gene_feature_names=[f"{gene}_log1p_cp10k" for gene in genes],
        engineered_feature_names=[],
        model_feature_names=[f"{gene}_log1p_cp10k" for gene in genes],
        gene_to_col={gene: index for index, gene in enumerate(genes)},
        engineered_to_col={},
    )
    annotation, _ = v2_eval.normalize_annotation({}, x, spec)
    predictions = v2_eval.predict_profiles(models["v2_profiles"], x, spec, annotation)
    first = predictions["v2_high_f1"][0]
    second = predictions["v2_high_purity"][0]
    if not np.allclose(first, second, rtol=0, atol=0):
        raise RuntimeError("V2 operating modes unexpectedly changed the underlying score")
    return first, {key: value[1] for key, value in predictions.items()}, annotation


def add_group_summary(parts: list[pd.DataFrame], keys: list[str], frame: pd.DataFrame) -> None:
    numeric = [column for column in frame.columns if column not in keys]
    parts.append(frame.groupby(keys, observed=True, dropna=False)[numeric].sum().reset_index())


def finalize_summary(parts: list[pd.DataFrame], keys: list[str]) -> pd.DataFrame:
    merged = pd.concat(parts, ignore_index=True)
    numeric = [column for column in merged.columns if column not in keys]
    return merged.groupby(keys, observed=True, dropna=False)[numeric].sum().reset_index()


def long_model_summary(wide: pd.DataFrame, group_keys: list[str]) -> pd.DataFrame:
    rows = []
    for _, row in wide.iterrows():
        for model in MODEL_KEYS:
            predicted = int(row[model])
            item = {key: row[key] for key in group_keys}
            item.update({
                "model": model,
                "n_cells": int(row["n_cells"]),
                "predicted_gdt": predicted,
                "predicted_fraction": predicted / row["n_cells"] if row["n_cells"] else np.nan,
            })
            for flag in FLAG_KEYS:
                count = int(row[f"{model}__{flag}"])
                item[f"predicted_{flag}"] = count
                item[f"{flag}_fraction_of_predictions"] = count / predicted if predicted else np.nan
            rows.append(item)
    return pd.DataFrame(rows)


def run(config_path: Path, config: dict, models: dict, preflight_result: dict) -> dict[str, Any]:
    started = time.monotonic()
    input_path = resolve(config["input_h5ad"])
    table_dir = resolve(config["output_table_dir"])
    log_dir = resolve(config["output_log_dir"])
    table_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    prediction_path = table_dir / "whole_atlas_predictions.parquet"
    temporary = table_dir / ".whole_atlas_predictions.tmp.parquet"
    temporary.unlink(missing_ok=True)
    union_genes = model_feature_union(models, config)
    source_parts: list[pd.DataFrame] = []
    source_tissue_parts: list[pd.DataFrame] = []
    lung_parts: list[pd.DataFrame] = []
    writer: pq.ParquetWriter | None = None
    before = input_path.stat()
    try:
        with h5py.File(input_path, "r") as handle:
            matrix = handle["X"]
            n_obs, n_vars = (int(value) for value in matrix.attrs["shape"])
            var_names = read_string_dataset(handle["var"]["gene"]).astype(str)
            var_lookup = {gene: index for index, gene in enumerate(var_names)}
            gene_indices = np.asarray([var_lookup[gene] for gene in union_genes], dtype=np.int64)
            lookup = {gene: index for index, gene in enumerate(union_genes)}
            obs = {column: read_obs_column(handle, column) for column in OBS_COLUMNS}
            chunk_size = int(config["chunk_size"])
            for start in range(0, n_obs, chunk_size):
                stop = min(start + chunk_size, n_obs)
                counts = read_csr_chunk(matrix, start, stop, n_vars)
                totals = np.asarray(counts.sum(axis=1), dtype=np.float64).ravel()
                if np.any(totals <= 0):
                    raise RuntimeError(f"Zero-library cells in rows {start}:{stop}")
                x = counts[:, gene_indices].toarray().astype(np.float32, copy=False)
                x *= (float(config["normalization_target_sum"]) / totals).astype(np.float32)[:, None]
                np.log1p(x, out=x)

                v46_score, v46_hf1, v46_hp, excluded = v4_scores(x, lookup, models)
                v2_score, v2_calls, v2_annotation = v2_scores(x, union_genes, models)
                source_primary = obs["source_cell_type"][start:stop]
                source_secondary = obs["source_cell_type_level1"][start:stop]
                author_nk = author_nk_mask(source_primary, source_secondary, config["nk_audit"]["source_annotation_excludes"])
                nk = nk_audit_flags(x, lookup, config, author_nk)
                paired_ab = as_bool(obs["has_TRA_TRB_paired"][start:stop])
                any_gd = as_bool(obs["has_any_gd_tcr"][start:stop])
                trdc = x[:, lookup["TRDC"]] > 0
                flags = {
                    **{key: nk[key] for key in ("likely_nk", "author_nk", "conservative_expression_nk", "shared_cytotoxic_ambiguous", "trdv_positive")},
                    "paired_ab_without_gd": paired_ab & ~any_gd,
                    "productive_gd": any_gd,
                    "trdc_positive": trdc,
                    "trdc_plus_trdv_minus": trdc & ~nk["trdv_positive"],
                }
                calls = {
                    "v4_6_highest_f1": v46_hf1,
                    "v4_6_high_purity": v46_hp,
                    "v2_high_f1": v2_calls["v2_high_f1"],
                    "v2_high_purity": v2_calls["v2_high_purity"],
                }
                data: dict[str, Any] = {
                    "cell_id": obs["cell_id"][start:stop].astype(str),
                    "source_gse_id": obs["source_gse_id"][start:stop].astype(str),
                    "source_accession": obs["source_accession_harmonized_v2"][start:stop].astype(str),
                    "sample_id": obs["sample_id_harmonized_v2"][start:stop].astype(str),
                    "donor_id": obs["donor_id_harmonized_v2"][start:stop].astype(str),
                    "tissue": obs["tissue_harmonized_v2"][start:stop].astype(str),
                    "specimen_context": obs["specimen_context_harmonized_v2"][start:stop].astype(str),
                    "tumor_type": obs["tumor_type_harmonized_v2"][start:stop].astype(str),
                    "source_cell_type": source_primary.astype(str),
                    "v4_6_score": v46_score.astype(np.float32),
                    "v2_score": v2_score.astype(np.float32),
                    "v2_expression_annotation": v2_annotation.astype(str),
                    "v4_6_cd4_treg_excluded": excluded,
                    "likely_nk": flags["likely_nk"],
                    "author_nk": flags["author_nk"],
                    "conservative_expression_nk": flags["conservative_expression_nk"],
                    "shared_cytotoxic_ambiguous": flags["shared_cytotoxic_ambiguous"],
                    "tier_a_nk_detected_count": nk["tier_a_detected_count"],
                    "innate_nk_detected_count": nk["innate_detected_count"],
                    "cd3_detected_count": nk["cd3_detected_count"],
                    "shared_cytotoxic_detected_count": nk["shared_cytotoxic_detected_count"],
                    "has_paired_ab_tcr": paired_ab,
                    "has_any_gd_tcr": any_gd,
                    "trdc_expressed": trdc,
                    "trdv_expressed": nk["trdv_positive"],
                }
                data.update(calls)
                cell_frame = pd.DataFrame(data)
                arrow = pa.Table.from_pandas(cell_frame, preserve_index=False)
                if writer is None:
                    writer = pq.ParquetWriter(temporary, arrow.schema, compression="zstd")
                writer.write_table(arrow, row_group_size=len(cell_frame))

                agg = pd.DataFrame({
                    "source_gse_id": data["source_gse_id"],
                    "tissue": data["tissue"],
                    "specimen_context": data["specimen_context"],
                    "tumor_type": data["tumor_type"],
                    "n_cells": np.ones(stop - start, dtype=np.int32),
                })
                for model, call in calls.items():
                    agg[model] = call.astype(np.int32)
                    for flag in FLAG_KEYS:
                        agg[f"{model}__{flag}"] = (call & flags[flag]).astype(np.int32)
                for label, left, right in (
                    ("highest_f1", v46_hf1, v2_calls["v2_high_f1"]),
                    ("high_purity", v46_hp, v2_calls["v2_high_purity"]),
                ):
                    agg[f"{label}__both"] = (left & right).astype(np.int32)
                    agg[f"{label}__v4_only"] = (left & ~right).astype(np.int32)
                    agg[f"{label}__v2_only"] = (~left & right).astype(np.int32)
                add_group_summary(source_parts, ["source_gse_id"], agg.drop(columns=["tissue", "specimen_context", "tumor_type"]))
                add_group_summary(source_tissue_parts, ["source_gse_id", "tissue", "specimen_context"], agg.drop(columns=["tumor_type"]))
                lung = agg.source_gse_id.eq(config["large_luad_lusc_source"])
                if lung.any():
                    add_group_summary(lung_parts, ["source_gse_id", "tumor_type", "specimen_context"], agg.loc[lung].drop(columns=["tissue"]))
                print(f"Scored {stop:,} / {n_obs:,} cells", flush=True)
    finally:
        if writer is not None:
            writer.close()
    os.replace(temporary, prediction_path)

    source_wide = finalize_summary(source_parts, ["source_gse_id"])
    source_tissue_wide = finalize_summary(source_tissue_parts, ["source_gse_id", "tissue", "specimen_context"])
    lung_wide = finalize_summary(lung_parts, ["source_gse_id", "tumor_type", "specimen_context"])
    source_long = long_model_summary(source_wide, ["source_gse_id"])
    source_tissue_long = long_model_summary(source_tissue_wide, ["source_gse_id", "tissue", "specimen_context"])
    lung_long = long_model_summary(lung_wide, ["source_gse_id", "tumor_type", "specimen_context"])
    source_long.to_csv(table_dir / "predictions_by_dataset.csv", index=False)
    source_tissue_long.to_csv(table_dir / "predictions_by_dataset_tissue_context.csv", index=False)
    lung_long.to_csv(table_dir / "gse243013_luad_lusc_predictions.csv", index=False)
    source_wide.to_parquet(table_dir / "source_aggregate_wide.parquet", index=False)

    comparison_rows = []
    for _, row in source_wide.iterrows():
        for mode, left, right in (
            ("highest_f1", "v4_6_highest_f1", "v2_high_f1"),
            ("high_purity", "v4_6_high_purity", "v2_high_purity"),
        ):
            both = int(row[f"{mode}__both"])
            v4_only = int(row[f"{mode}__v4_only"])
            v2_only = int(row[f"{mode}__v2_only"])
            union = both + v4_only + v2_only
            comparison_rows.append({
                "source_gse_id": row.source_gse_id,
                "mode": mode,
                "n_cells": int(row.n_cells),
                "v4_6_predicted": int(row[left]),
                "v2_predicted": int(row[right]),
                "both": both,
                "v4_6_only": v4_only,
                "v2_only": v2_only,
                "jaccard": both / union if union else np.nan,
            })
    comparison = pd.DataFrame(comparison_rows)
    comparison.to_csv(table_dir / "v4_6_vs_v2_by_dataset.csv", index=False)

    overall_wide = source_wide.drop(columns=["source_gse_id"]).sum(numeric_only=True)
    overall_long = long_model_summary(pd.DataFrame([overall_wide]), [])
    overall_long.to_csv(table_dir / "whole_atlas_overall.csv", index=False)
    after = input_path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError("Input H5AD size or mtime changed during read-only inference")
    if int(source_wide.n_cells.sum()) != int(config["expected_cells"]):
        raise RuntimeError("Per-source totals do not reproduce the whole atlas")
    parquet_rows = pq.ParquetFile(prediction_path).metadata.num_rows
    if parquet_rows != int(config["expected_cells"]):
        raise RuntimeError("Per-cell Parquet row count differs from atlas cells")
    summary = {
        "status": "PASS_V4_6_V2_WHOLE_ATLAS_INFERENCE",
        "runtime_seconds": time.monotonic() - started,
        "protocol_version": config["protocol_version"],
        "config_sha256": sha256_file(config_path),
        "input_candidate_sha256": config["expected_input_sha256"],
        "input_size_mtime_unchanged": True,
        "n_cells": int(config["expected_cells"]),
        "n_sources": int(source_wide.source_gse_id.nunique()),
        "prediction_parquet": str(prediction_path),
        "prediction_parquet_rows": parquet_rows,
        "prediction_parquet_sha256": sha256_file(prediction_path),
        "models_or_thresholds_retuned": False,
        "h5ad_mutated": False,
        "overall": overall_long.to_dict("records"),
        "gse243013": lung_long.to_dict("records"),
        "preflight": preflight_result,
    }
    (log_dir / "inference_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--preflight-only", action="store_true")
    parser.add_argument("--smoke-test", action="store_true")
    args = parser.parse_args()
    config_path = resolve(args.config)
    config = json.loads(config_path.read_text())
    models = prepare_models(config)
    result = preflight(config_path, config, models)
    if args.smoke_test:
        result = smoke_test(config, models)
    elif not args.preflight_only:
        result = run(config_path, config, models, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
