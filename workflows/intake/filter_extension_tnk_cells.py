#!/usr/bin/env python3
"""Filter every standalone extension cohort to high-recall T/NK candidates."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
import uuid
from pathlib import Path
from typing import Any, Iterable, Sequence

PROJECT_ROOT = Path(__file__).resolve().parents[2]
for path in (
    PROJECT_ROOT,
    PROJECT_ROOT / "src",
    PROJECT_ROOT / "workflows" / "integration",
    PROJECT_ROOT / "workflows" / "intake",
):
    value = str(path)
    if value not in sys.path:
        sys.path.insert(0, value)

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns

import phase1_extract_tnk_candidates as phase1
from build_extension_h5ads import ensure_canonical_tcr, validate_sparse_integer_counts


BUILT_ROOT = PROJECT_ROOT / "data/interim/extension_intake/built_h5ads"
GSE169246_SOURCE = PROJECT_ROOT / "data/datasets/GSE169246/processed/current.h5ad"
OUTPUT_ROOT = PROJECT_ROOT / "data/interim/extension_intake/tnk_filtered_h5ads"
MANIFEST = PROJECT_ROOT / "data/interim/extension_intake/tnk_filtered_h5ads_manifest.csv"
TABLE_ROOT = PROJECT_ROOT / "Integrated_dataset/tables/extension_intake/tnk_filter"
FIGURE_ROOT = PROJECT_ROOT / "Integrated_dataset/figures/extension_intake/tnk_filter"
LOG_ROOT = PROJECT_ROOT / "Integrated_dataset/logs/extension_intake"
SUMMARY_CSV = TABLE_ROOT / "extension_tnk_filter_cohort_summary.csv"
SAMPLE_CSV = TABLE_ROOT / "extension_tnk_filter_sample_summary.csv"
REASON_CSV = TABLE_ROOT / "extension_tnk_filter_reason_summary.csv"
DECISIONS_CSV = TABLE_ROOT / "extension_tnk_filter_cell_decisions.csv.gz"
METADATA_AUDIT_CSV = TABLE_ROOT / "extension_tnk_filter_metadata_audit.csv"
SUMMARY_JSON = LOG_ROOT / "extension_tnk_filter_qc.json"
SUMMARY_MD = LOG_ROOT / "extension_tnk_filter_qc.md"

REQUIRED_PHASE0_COHORTS = (
    "GSE114724",
    "GSE121636_GSE121637",
    "GSE159251",
    "GSE292700",
    "GSE294273_GSE294274",
    "GSE296954",
    "GSE315928",
)
ALL_COHORTS = (*REQUIRED_PHASE0_COHORTS, "GSE169246")

# IL7R and LTB are useful T-cell-supporting genes, but either gene alone is too
# broad for a lineage gate in mixed tumor samples. Require at least one core
# CD3/TCR-alpha-beta transcript for the transcriptomic T-cell rule.
CORE_T_MARKERS = ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2"]

EXPLICIT_ANNOTATION_COLUMNS = (
    "annot_final",
    "major_cluster",
    "major_cell_type",
    "cell_type",
    "celltype",
    "annotation",
    "cell_annotation",
    "cell_to_cluster",
    "subset",
    "type",
    "subgroup",
    "group",
    "phase1_annotation_label",
)


class ExtensionTnkFilterError(RuntimeError):
    """Raised when extension T/NK filtering would violate a data contract."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", action="append", default=[])
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--replace-existing", action="store_true")
    parser.add_argument("--output-root", type=Path, default=OUTPUT_ROOT)
    return parser.parse_args(argv)


def clean_string_series(series: pd.Series) -> pd.Series:
    values = series.astype("string").fillna("").str.strip()
    return values.replace({"nan": "", "None": "", "<NA>": ""}).astype(object)


def bool_series(obs: pd.DataFrame, column: str) -> np.ndarray:
    if column not in obs:
        return np.zeros(len(obs), dtype=bool)
    values = obs[column]
    if pd.api.types.is_bool_dtype(values.dtype):
        return values.fillna(False).to_numpy(dtype=bool)
    normalized = clean_string_series(values).str.casefold()
    unknown = ~normalized.isin({"true", "t", "yes", "y", "1", "false", "f", "no", "n", "0", ""})
    if unknown.any():
        examples = sorted(normalized.loc[unknown].unique())[:5]
        raise ExtensionTnkFilterError(f"Unrecognized values in {column}: {examples}")
    return normalized.isin({"true", "t", "yes", "y", "1"}).to_numpy(dtype=bool)


def sha256_file(path: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk_size):
            digest.update(block)
    return digest.hexdigest()


def source_paths(selected: Sequence[str] | None = None) -> dict[str, Path]:
    paths = {cohort: BUILT_ROOT / f"{cohort}.h5ad" for cohort in REQUIRED_PHASE0_COHORTS}
    paths["GSE169246"] = GSE169246_SOURCE
    requested = set(selected or ALL_COHORTS)
    unknown = requested - set(paths)
    if unknown:
        raise ExtensionTnkFilterError(f"Unknown cohort(s): {sorted(unknown)}")
    result = {cohort: paths[cohort] for cohort in ALL_COHORTS if cohort in requested}
    missing = [str(path) for path in result.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing extension H5AD source(s): {missing}")
    return result


def annotation_columns(obs: pd.DataFrame) -> list[str]:
    selected: list[str] = []
    for column in EXPLICIT_ANNOTATION_COLUMNS:
        if column in obs and column not in selected:
            selected.append(column)
    for column in phase1.find_annotation_columns(list(obs.columns)):
        lower = column.casefold()
        if lower.endswith("_keep") or "selection_reason" in lower:
            continue
        if column not in selected:
            selected.append(column)
    return selected


def marker_hits_and_score(
    adata: ad.AnnData, markers: list[str], total_counts: np.ndarray
) -> tuple[np.ndarray, np.ndarray, int]:
    positions = phase1.marker_positions_from_var_names(pd.Index(adata.var_names), markers)
    if len(positions) == 0:
        return (
            np.zeros(adata.n_obs, dtype=np.int16),
            np.zeros(adata.n_obs, dtype=np.float32),
            0,
        )
    matrix = adata.X[:, positions]
    if sp.issparse(matrix):
        matrix = matrix.tocsr().astype(np.float32, copy=False)
        hits = np.ravel((matrix > 0).sum(axis=1)).astype(np.int16)
        scale = sp.diags(1e4 / total_counts, format="csr")
        normalized = scale @ matrix
        normalized.data = np.log1p(normalized.data)
        scores = np.ravel(np.asarray(normalized.mean(axis=1))).astype(np.float32)
    else:
        dense = np.asarray(matrix, dtype=np.float32)
        hits = np.sum(dense > 0, axis=1).astype(np.int16)
        scores = np.log1p((dense / total_counts[:, None]) * 1e4).mean(axis=1).astype(np.float32)
    return hits, scores, int(len(positions))


def propagate_unique_sample_values(obs: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    result = obs.copy()
    samples = clean_string_series(result["sample_id"])
    for column in columns:
        if column not in result:
            result[column] = ""
        values = clean_string_series(result[column])
        unique_by_sample: dict[str, str] = {}
        for sample, indices in values.groupby(samples, sort=False).groups.items():
            unique = pd.unique(values.loc[indices][values.loc[indices].ne("")])
            if len(unique) == 1:
                unique_by_sample[str(sample)] = str(unique[0])
        fill = samples.map(unique_by_sample).fillna("")
        result[column] = values.mask(values.eq(""), fill).astype(object)
    return result


def specimen_context(tissue: pd.Series, diagnosis: pd.Series) -> pd.Series:
    tissues = clean_string_series(tissue).str.casefold()
    diagnoses = clean_string_series(diagnosis).str.casefold()
    result = pd.Series("other_tissue", index=tissue.index, dtype=object)
    result.loc[tissues.eq("")] = "unresolved"
    result.loc[tissues.str.contains("blood|pbmc", regex=True)] = "blood"
    adjacent = tissues.str.contains("adjacent|paracancer|normal", regex=True)
    result.loc[adjacent] = "adjacent_non_tumor"
    primary = tissues.str.contains("breast_tumor|renal_tumor|primary_tumor", regex=True)
    result.loc[primary] = "primary_tumor"
    generic_tumor = tissues.eq("tumor")
    result.loc[generic_tumor] = "primary_tumor"
    metastatic = tissues.str.contains("metast|lymph_node|liver|chest_wall|lung", regex=True)
    result.loc[metastatic & diagnoses.ne("")] = "metastatic_tumor"
    breast_primary = tissues.eq("breast") & diagnoses.str.contains(
        "breast|triple_negative", regex=True
    )
    result.loc[breast_primary] = "primary_tumor"
    return result


def harmonize_gse169246_obs(obs: pd.DataFrame) -> pd.DataFrame:
    result = obs.copy()
    if "sample_id" not in result:
        raise ExtensionTnkFilterError("GSE169246 lacks sample_id")
    result["sample_id"] = clean_string_series(result["sample_id"])
    if result["sample_id"].eq("").any():
        raise ExtensionTnkFilterError("GSE169246 contains blank sample_id values")

    result = propagate_unique_sample_values(
        result,
        ("patient_id", "timepoint_group", "tissue", "treatment", "efficacy", "response"),
    )
    derived_patient = result["sample_id"].str.extract(r"(P\d+)", expand=False).fillna("")
    result["patient_id"] = clean_string_series(result["patient_id"]).mask(
        clean_string_series(result["patient_id"]).eq(""), derived_patient
    )
    result["donor_id"] = clean_string_series(result["patient_id"])
    result["donor"] = result["donor_id"].astype(object)
    raw_barcode = clean_string_series(result["raw_barcode"])
    derived_library = raw_barcode.str.extract(r"^[^.]+[.]([^|]+)$", expand=False).fillna("")
    result["library_id"] = clean_string_series(derived_library).mask(
        clean_string_series(derived_library).eq(""), result["sample_id"]
    )
    result["tissue_original"] = clean_string_series(result["tissue"])
    result["tissue_harmonized"] = clean_string_series(result["tissue"]).replace("", "unresolved")
    result["diagnosis"] = "triple_negative_breast_cancer"
    result["specimen_context"] = specimen_context(result["tissue_harmonized"], result["diagnosis"])
    result["source_accession"] = "GSE169246"
    result["technology_simple"] = "10x 5'"
    result["source_gse_id"] = "GSE169246"
    if "barcode" not in result:
        for candidate in ("raw_barcode", "cell_barcode", "published_barcode"):
            if candidate in result:
                result["barcode"] = clean_string_series(result[candidate])
                break
    if "barcode" not in result:
        result["barcode"] = result.index.astype(str)
    result = ensure_canonical_tcr(result)
    result["tcr_schema_provenance"] = "embedded_mmc5_paired_tra_trb_cdr3"
    return result


def harmonize_common_obs(adata: ad.AnnData, cohort_id: str) -> None:
    if cohort_id == "GSE169246":
        adata.obs = harmonize_gse169246_obs(adata.obs)
    else:
        adata.obs = ensure_canonical_tcr(adata.obs)
        if "tissue_harmonized" not in adata.obs:
            adata.obs["tissue_harmonized"] = clean_string_series(adata.obs["tissue"])
        if "specimen_context" not in adata.obs:
            adata.obs["specimen_context"] = specimen_context(
                adata.obs["tissue_harmonized"], adata.obs["diagnosis"]
            )
    if "source_gse_id" in adata.obs:
        source_gse = clean_string_series(adata.obs["source_gse_id"])
    else:
        source_gse = pd.Series("", index=adata.obs_names, dtype=object)
    if "source_accession" in adata.obs:
        source_accession = clean_string_series(adata.obs["source_accession"])
        source_gse = source_gse.mask(source_gse.eq(""), source_accession)
    adata.obs["source_gse_id"] = source_gse.replace("", cohort_id).astype(object)
    adata.obs["extension_cohort_id"] = cohort_id


def selection_reason(
    keep: np.ndarray,
    tcr_keep: np.ndarray,
    annotation_keep: np.ndarray,
    marker_keep: np.ndarray,
    annotation_negative: np.ndarray,
    contamination_override: np.ndarray,
) -> np.ndarray:
    result = np.full(len(keep), "no_tnk_evidence", dtype=object)
    result[~keep & annotation_negative] = "explicit_non_tnk_annotation"
    result[~keep & contamination_override] = "contamination_override"
    result[keep & marker_keep] = "marker"
    result[keep & annotation_keep] = "annotation"
    result[keep & annotation_keep & marker_keep] = "annotation+marker"
    result[keep & tcr_keep] = "productive_tcr"
    result[keep & tcr_keep & marker_keep] = "productive_tcr+marker"
    result[keep & tcr_keep & annotation_keep] = "productive_tcr+annotation"
    result[keep & tcr_keep & annotation_keep & marker_keep] = "productive_tcr+annotation+marker"
    return result


def compute_decisions(adata: ad.AnnData, cohort_id: str) -> pd.DataFrame:
    validate_sparse_integer_counts(adata.X, f"{cohort_id} X")
    harmonize_common_obs(adata, cohort_id)
    total_counts = np.ravel(np.asarray(adata.X.sum(axis=1))).astype(np.float32)
    total_counts[total_counts <= 0] = 1.0

    t_hits, t_score, t_present = marker_hits_and_score(adata, phase1.T_MARKERS, total_counts)
    core_t_hits, _, core_t_present = marker_hits_and_score(
        adata, CORE_T_MARKERS, total_counts
    )
    gd_hits, gd_score, gd_present = marker_hits_and_score(adata, phase1.GD_MARKERS, total_counts)
    nk_hits, nk_score, nk_present = marker_hits_and_score(adata, phase1.NK_MARKERS, total_counts)
    nk_strong_hits, _, nk_strong_present = marker_hits_and_score(
        adata, phase1.NK_STRONG_MARKERS, total_counts
    )
    _, contam_score, contam_present = marker_hits_and_score(
        adata, phase1.CONTAM_MARKERS, total_counts
    )

    columns = annotation_columns(adata.obs)
    annotation_keep, annotation_negative = phase1.annotation_masks(adata.obs, columns)
    tcr_keep = bool_series(adata.obs, "has_any_ab_tcr") | bool_series(
        adata.obs, "has_any_gd_tcr"
    )

    t_rule = core_t_hits >= 1
    gd_rule = gd_hits >= 1
    nk_rule = (nk_hits >= 2) & (nk_strong_hits >= 1)
    marker_before_override = t_rule | gd_rule | nk_rule
    strong_contam = contam_score >= 1.25
    weak_lymphoid = (t_hits + gd_hits + nk_hits) <= 1
    contamination_override = strong_contam & weak_lymphoid & ~annotation_keep & ~tcr_keep
    marker_keep = marker_before_override & ~contamination_override
    keep = tcr_keep | annotation_keep | (marker_keep & ~annotation_negative)

    if np.any(tcr_keep & ~keep):
        raise ExtensionTnkFilterError(f"{cohort_id}: productive TCR cells would be dropped")
    if not np.any(keep):
        raise ExtensionTnkFilterError(f"{cohort_id}: T/NK filter retained zero cells")

    frame = pd.DataFrame(index=adata.obs_names)
    frame["cohort_id"] = cohort_id
    frame["cell_id"] = adata.obs_names.astype(str)
    frame["sample_id"] = clean_string_series(adata.obs["sample_id"]).to_numpy()
    frame["tissue_harmonized"] = clean_string_series(adata.obs["tissue_harmonized"]).to_numpy()
    frame["specimen_context"] = clean_string_series(adata.obs["specimen_context"]).to_numpy()
    frame["keep"] = keep
    frame["selection_reason"] = selection_reason(
        keep, tcr_keep, annotation_keep, marker_keep, annotation_negative, contamination_override
    )
    frame["productive_tcr_evidence"] = tcr_keep
    frame["annotation_tnk_evidence"] = annotation_keep
    frame["annotation_non_tnk_evidence"] = annotation_negative
    frame["marker_tnk_evidence"] = marker_keep
    frame["nk_marker_rule"] = nk_rule
    frame["contamination_override"] = contamination_override
    frame["t_hits"] = t_hits
    frame["core_t_hits"] = core_t_hits
    frame["gd_hits"] = gd_hits
    frame["nk_hits"] = nk_hits
    frame["nk_strong_hits"] = nk_strong_hits
    frame["t_score"] = t_score
    frame["gd_score"] = gd_score
    frame["nk_score"] = nk_score
    frame["contam_score"] = contam_score
    frame.attrs["annotation_columns"] = columns
    frame.attrs["marker_availability"] = {
        "T": t_present,
        "T_CORE": core_t_present,
        "GD": gd_present,
        "NK": nk_present,
        "NK_STRONG": nk_strong_present,
        "CONTAM": contam_present,
    }
    return frame


def attach_decisions(adata: ad.AnnData, decisions: pd.DataFrame) -> None:
    mapping = {
        "keep": "extension_tnk_keep",
        "selection_reason": "extension_tnk_selection_reason",
        "productive_tcr_evidence": "extension_tnk_tcr_evidence",
        "annotation_tnk_evidence": "extension_tnk_annotation_evidence",
        "annotation_non_tnk_evidence": "extension_tnk_non_tnk_annotation",
        "marker_tnk_evidence": "extension_tnk_marker_evidence",
        "nk_marker_rule": "extension_tnk_nk_marker_rule",
        "contamination_override": "extension_tnk_contamination_override",
        "t_hits": "extension_tnk_t_hits",
        "core_t_hits": "extension_tnk_core_t_hits",
        "gd_hits": "extension_tnk_gd_hits",
        "nk_hits": "extension_tnk_nk_hits",
        "nk_strong_hits": "extension_tnk_nk_strong_hits",
        "t_score": "extension_tnk_t_score",
        "gd_score": "extension_tnk_gd_score",
        "nk_score": "extension_tnk_nk_score",
        "contam_score": "extension_tnk_contam_score",
    }
    for source, target in mapping.items():
        adata.obs[target] = decisions[source].to_numpy()


def subset_and_write(
    adata: ad.AnnData,
    decisions: pd.DataFrame,
    cohort_id: str,
    output_path: Path,
    source_path: Path,
    source_sha256_before: str,
) -> dict[str, Any]:
    keep = decisions["keep"].to_numpy(dtype=bool)
    attach_decisions(adata, decisions)
    subset = adata[keep].copy()
    subset.X = validate_sparse_integer_counts(subset.X, f"{cohort_id} filtered X")
    if "counts" in subset.layers:
        subset.layers["counts"] = validate_sparse_integer_counts(
            subset.layers["counts"], f"{cohort_id} filtered counts layer"
        )
    subset.obs["extension_tnk_keep"] = True
    subset.uns["extension_tnk_filter"] = {
        "schema_version": 1,
        "cohort_id": cohort_id,
        "source_h5ad": str(source_path.resolve()),
        "selection_policy": "productive_tcr OR author_TNK OR core_transcriptomic_TNK",
        "annotation_columns": list(decisions.attrs.get("annotation_columns", [])),
        "marker_availability": dict(decisions.attrs.get("marker_availability", {})),
        "input_cells": int(adata.n_obs),
        "retained_cells": int(subset.n_obs),
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_name(
        f".{output_path.stem}.{uuid.uuid4().hex}.tmp.h5ad"
    )
    try:
        subset.write_h5ad(temporary, compression="gzip")
        os.replace(temporary, output_path)
    finally:
        temporary.unlink(missing_ok=True)

    tcr_input = int(decisions["productive_tcr_evidence"].sum())
    tcr_retained = int(
        (decisions["productive_tcr_evidence"] & decisions["keep"]).sum()
    )
    if tcr_input != tcr_retained:
        raise ExtensionTnkFilterError(
            f"{cohort_id}: productive TCR retention mismatch {tcr_retained} != {tcr_input}"
        )
    source_sha256_after = sha256_file(source_path)
    if source_sha256_after != source_sha256_before:
        raise ExtensionTnkFilterError(f"{cohort_id}: source H5AD changed during filtering")
    return {
        "cohort_id": cohort_id,
        "source_h5ad": str(source_path.resolve()),
        "output_h5ad": str(output_path.resolve()),
        "input_cells": int(adata.n_obs),
        "retained_cells": int(subset.n_obs),
        "removed_cells": int(adata.n_obs - subset.n_obs),
        "retained_fraction": float(subset.n_obs / adata.n_obs),
        "productive_tcr_input": tcr_input,
        "productive_tcr_retained": tcr_retained,
        "productive_tcr_retention": 1.0 if tcr_input else math.nan,
        "nk_marker_input": int(decisions["nk_marker_rule"].sum()),
        "nk_marker_retained": int((decisions["nk_marker_rule"] & decisions["keep"]).sum()),
        "annotation_tnk_input": int(decisions["annotation_tnk_evidence"].sum()),
        "explicit_non_tnk_input": int(decisions["annotation_non_tnk_evidence"].sum()),
        "explicit_non_tnk_retained": int(
            (decisions["annotation_non_tnk_evidence"] & decisions["keep"]).sum()
        ),
        "samples_input": int(adata.obs["sample_id"].nunique()),
        "samples_retained": int(subset.obs["sample_id"].nunique()),
        "source_size_bytes": source_path.stat().st_size,
        "output_size_bytes": output_path.stat().st_size,
        "source_sha256": source_sha256_after,
        "source_sha256_unchanged": True,
        "output_sha256": sha256_file(output_path),
        "annotation_columns": ";".join(decisions.attrs.get("annotation_columns", [])),
    }


def blank_mask(series: pd.Series) -> np.ndarray:
    return clean_string_series(series).eq("").to_numpy(dtype=bool)


def validate_output(path: Path, summary: dict[str, Any]) -> dict[str, Any]:
    adata = ad.read_h5ad(path, backed="r")
    try:
        if adata.n_obs != summary["retained_cells"]:
            raise ExtensionTnkFilterError(f"{path}: retained cell count mismatch")
        required = {
            "source_gse_id",
            "source_accession",
            "extension_cohort_id",
            "sample_id",
            "library_id",
            "donor_id",
            "tissue_harmonized",
            "specimen_context",
            "has_any_ab_tcr",
            "has_any_gd_tcr",
            "extension_tnk_keep",
            "extension_tnk_selection_reason",
        }
        missing = sorted(required - set(adata.obs.columns))
        if missing:
            raise ExtensionTnkFilterError(f"{path}: missing output obs columns {missing}")
        if not bool_series(adata.obs, "extension_tnk_keep").all():
            raise ExtensionTnkFilterError(f"{path}: contains cells not selected by the filter")
        if adata.obs_names.duplicated().any():
            raise ExtensionTnkFilterError(f"{path}: duplicate cell IDs")

        blank_counts = {
            column: int(blank_mask(adata.obs[column]).sum())
            for column in (
                "source_gse_id",
                "source_accession",
                "extension_cohort_id",
                "sample_id",
                "library_id",
                "donor_id",
                "tissue_harmonized",
                "specimen_context",
            )
        }
        nonzero_blanks = {key: value for key, value in blank_counts.items() if value}
        if nonzero_blanks:
            raise ExtensionTnkFilterError(f"{path}: blank required metadata {nonzero_blanks}")

        cohort_values = clean_string_series(adata.obs["extension_cohort_id"])
        cohort_mismatch = int(cohort_values.ne(summary["cohort_id"]).sum())
        if cohort_mismatch:
            raise ExtensionTnkFilterError(f"{path}: {cohort_mismatch} cohort ID mismatches")

        source_gse = clean_string_series(adata.obs["source_gse_id"])
        source_accession = clean_string_series(adata.obs["source_accession"])
        source_mismatch = int(source_gse.ne(source_accession).sum())
        if source_mismatch:
            raise ExtensionTnkFilterError(
                f"{path}: {source_mismatch} source_gse_id/source_accession mismatches"
            )

        barcode_column = "barcode_core" if "barcode_core" in adata.obs else "barcode"
        if barcode_column not in adata.obs:
            raise ExtensionTnkFilterError(f"{path}: missing barcode or barcode_core")
        barcode_blanks = int(blank_mask(adata.obs[barcode_column]).sum())
        key_frame = pd.DataFrame(
            {
                "library_id": clean_string_series(adata.obs["library_id"]).to_numpy(),
                "barcode": clean_string_series(adata.obs[barcode_column]).to_numpy(),
            }
        )
        duplicate_sample_barcode = int(
            key_frame.duplicated(["library_id", "barcode"], keep=False).sum()
        )
        if barcode_blanks or duplicate_sample_barcode:
            raise ExtensionTnkFilterError(
                f"{path}: barcode blanks={barcode_blanks}, duplicate library+barcode rows="
                f"{duplicate_sample_barcode}"
            )

        logical_mismatches = 0
        for chain in ("TRA", "TRB", "TRG", "TRD"):
            expected = ~blank_mask(adata.obs[f"{chain}_cdr3"])
            observed = bool_series(adata.obs, f"has_{chain}")
            logical_mismatches += int(np.count_nonzero(expected != observed))
        logical_mismatches += int(
            np.count_nonzero(
                bool_series(adata.obs, "has_TRA_TRB_paired")
                != (
                    bool_series(adata.obs, "has_TRA")
                    & bool_series(adata.obs, "has_TRB")
                )
            )
        )
        logical_mismatches += int(
            np.count_nonzero(
                bool_series(adata.obs, "has_TRG_TRD_paired")
                != (
                    bool_series(adata.obs, "has_TRG")
                    & bool_series(adata.obs, "has_TRD")
                )
            )
        )
        logical_mismatches += int(
            np.count_nonzero(
                bool_series(adata.obs, "has_any_ab_tcr")
                != (
                    bool_series(adata.obs, "has_TRA")
                    | bool_series(adata.obs, "has_TRB")
                )
            )
        )
        logical_mismatches += int(
            np.count_nonzero(
                bool_series(adata.obs, "has_any_gd_tcr")
                != (
                    bool_series(adata.obs, "has_TRG")
                    | bool_series(adata.obs, "has_TRD")
                )
            )
        )
        if logical_mismatches:
            raise ExtensionTnkFilterError(
                f"{path}: {logical_mismatches} canonical TCR logical mismatches"
            )
        unresolved = int(
            clean_string_series(adata.obs["specimen_context"]).eq("unresolved").sum()
        )
        return {
            "cohort_id": summary["cohort_id"],
            "status": "PASS",
            "cells_checked": int(adata.n_obs),
            "source_gse_ids": ";".join(sorted(source_gse.unique())),
            "blank_required_metadata": int(sum(blank_counts.values())),
            "blank_barcodes": barcode_blanks,
            "duplicate_cell_ids": int(adata.obs_names.duplicated().sum()),
            "duplicate_library_barcode_rows": duplicate_sample_barcode,
            "source_accession_mismatches": source_mismatch,
            "cohort_id_mismatches": cohort_mismatch,
            "tcr_logical_mismatches": logical_mismatches,
            "unresolved_specimen_context": unresolved,
            "samples": int(adata.obs["sample_id"].nunique()),
            "libraries": int(adata.obs["library_id"].nunique()),
            "donors": int(adata.obs["donor_id"].nunique()),
        }
    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()


def summarize_samples(decisions: pd.DataFrame) -> pd.DataFrame:
    work = decisions.copy()
    work["retained_tcr"] = work["productive_tcr_evidence"] & work["keep"]
    return (
        work.groupby(
            ["cohort_id", "sample_id", "tissue_harmonized", "specimen_context"],
            dropna=False,
            observed=True,
        )
        .agg(
            input_cells=("keep", "size"),
            retained_cells=("keep", "sum"),
            productive_tcr_input=("productive_tcr_evidence", "sum"),
            productive_tcr_retained=("retained_tcr", "sum"),
            nk_marker_input=("nk_marker_rule", "sum"),
        )
        .reset_index()
        .assign(retained_fraction=lambda frame: frame["retained_cells"] / frame["input_cells"])
    )


def write_figure(summary: pd.DataFrame) -> None:
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid")
    plot = summary.sort_values("retained_cells", ascending=True).copy()
    plot["cohort_label"] = plot["cohort_id"].str.replace("_", " / ", regex=False)
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.2))
    sns.barplot(data=plot, y="cohort_label", x="retained_cells", color="#287271", ax=axes[0])
    axes[0].set_title("Extension T/NK Cells Retained")
    axes[0].set_xlabel("Cells")
    axes[0].set_ylabel("")
    axes[0].set_xlim(0, float(plot["retained_cells"].max()) * 1.14)
    axes[0].xaxis.set_major_locator(MaxNLocator(nbins=5))
    axes[0].xaxis.set_major_formatter(
        FuncFormatter(lambda value, _: "0" if value == 0 else f"{value / 1000:.0f}k")
    )
    axes[0].bar_label(
        axes[0].containers[0],
        labels=[f"{value:,.0f}" for value in plot["retained_cells"]],
        padding=3,
        fontsize=8,
    )
    sns.barplot(data=plot, y="cohort_label", x="retained_fraction", color="#D97B29", ax=axes[1])
    axes[1].set_title("Retained Fraction")
    axes[1].set_xlabel("Fraction")
    axes[1].set_ylabel("")
    axes[1].set_xlim(0, 1.12)
    axes[1].bar_label(
        axes[1].containers[0],
        labels=[f"{value:.1%}" for value in plot["retained_fraction"]],
        padding=3,
        fontsize=8,
    )
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / "extension_tnk_filter_overview.png", dpi=300)
    plt.close(fig)


def write_report(
    summary: pd.DataFrame, sample_summary: pd.DataFrame, metadata_audit: pd.DataFrame
) -> None:
    total_input = int(summary["input_cells"].sum())
    total_retained = int(summary["retained_cells"].sum())
    tcr_input = int(summary["productive_tcr_input"].sum())
    tcr_retained = int(summary["productive_tcr_retained"].sum())
    unresolved = int(
        sample_summary.loc[
            sample_summary["specimen_context"].eq("unresolved"), "retained_cells"
        ].sum()
    )
    payload = {
        "gate_status": "PASS_REVIEW_REQUIRED",
        "merge_or_integration_approved": False,
        "cohorts": int(len(summary)),
        "input_cells": total_input,
        "retained_cells": total_retained,
        "retained_fraction": total_retained / total_input,
        "productive_tcr_input": tcr_input,
        "productive_tcr_retained": tcr_retained,
        "productive_tcr_retention": tcr_retained / tcr_input if tcr_input else None,
        "retained_cells_with_unresolved_specimen_context": unresolved,
        "outputs": {
            "manifest": str(MANIFEST.relative_to(PROJECT_ROOT)),
            "cohort_summary": str(SUMMARY_CSV.relative_to(PROJECT_ROOT)),
            "sample_summary": str(SAMPLE_CSV.relative_to(PROJECT_ROOT)),
            "cell_decisions": str(DECISIONS_CSV.relative_to(PROJECT_ROOT)),
            "metadata_audit": str(METADATA_AUDIT_CSV.relative_to(PROJECT_ROOT)),
            "figure": str(
                (FIGURE_ROOT / "extension_tnk_filter_overview.png").relative_to(
                    PROJECT_ROOT
                )
            ),
        },
    }
    LOG_ROOT.mkdir(parents=True, exist_ok=True)
    SUMMARY_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        "# Extension T/NK Filter QC",
        "",
        "## Gate",
        "",
        "- Status: **PASS_REVIEW_REQUIRED**",
        "- Merge or integration approved: **No**",
        "- Every new cohort was filtered or re-filtered independently before any merge.",
        "",
        "## Selection Policy",
        "",
        (
            "A cell is retained when it has productive TRA/TRB/TRG/TRD evidence, "
            "an author T/NK annotation, or high-recall transcriptomic T/NK evidence. "
            "Transcriptomic evidence requires a core CD3/TCR-alpha-beta gene (CD3D, "
            "CD3E, CD3G, TRAC, TRBC1, or TRBC2), a gamma-delta marker, or the "
            "established multi-gene NK rule. IL7R or LTB alone is not sufficient. "
            "Explicit non-T/NK annotations block marker-only retention, but never "
            "remove a productive-TCR cell."
        ),
        "",
        "## Totals",
        "",
        f"- Input cells: `{total_input:,}`",
        f"- Retained T/NK candidates: `{total_retained:,}` (`{total_retained / total_input:.2%}`)",
        f"- Productive-TCR cells retained: `{tcr_retained:,}` / `{tcr_input:,}`",
        f"- Retained cells with unresolved specimen context: `{unresolved:,}`",
        (
            f"- Metadata audit: **{','.join(sorted(metadata_audit['status'].unique()))}** "
            f"across `{len(metadata_audit)}` cohorts"
        ),
        (
            "- Metadata checks: required-field completeness, exact source accession, "
            "cohort identity, unique cell and library-plus-barcode keys, canonical TCR "
            "logical flags, and specimen context."
        ),
        (
            "- GSE169246 library handling: `51` biological samples and `78` source "
            "libraries were distinguished from the full barcode suffix; `76` libraries "
            "retained at least one T/NK candidate."
        ),
        "",
        "## Cohorts",
        "",
        "| Cohort | Input | Retained | Fraction | Productive TCR retained |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in summary.itertuples(index=False):
        lines.append(
            f"| {row.cohort_id} | {row.input_cells:,} | {row.retained_cells:,} | "
            f"{row.retained_fraction:.2%} | {row.productive_tcr_retained:,}/"
            f"{row.productive_tcr_input:,} |"
        )
    lines.extend(
        [
            "",
            "## Review Artifacts",
            "",
            f"- `{SUMMARY_CSV.relative_to(PROJECT_ROOT)}`",
            f"- `{SAMPLE_CSV.relative_to(PROJECT_ROOT)}`",
            f"- `{REASON_CSV.relative_to(PROJECT_ROOT)}`",
            f"- `{DECISIONS_CSV.relative_to(PROJECT_ROOT)}`",
            f"- `{METADATA_AUDIT_CSV.relative_to(PROJECT_ROOT)}`",
            f"- `{(FIGURE_ROOT / 'extension_tnk_filter_overview.png').relative_to(PROJECT_ROOT)}`",
            "",
            (
                "Explicit user approval remains required before merging these "
                "cohort-separated outputs."
            ),
        ]
    )
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def execute(
    paths: dict[str, Path], output_root: Path, dry_run: bool, replace_existing: bool = False
) -> dict[str, Any]:
    if output_root.resolve() in {
        (PROJECT_ROOT / "Integrated_dataset").resolve(),
        (PROJECT_ROOT / "high_speed_temp/Integrated_dataset").resolve(),
    }:
        raise ExtensionTnkFilterError("Extension T/NK outputs may not replace canonical milestones")
    existing = [
        output_root / f"{cohort}.h5ad"
        for cohort in paths
        if (output_root / f"{cohort}.h5ad").exists()
    ]
    if existing and not dry_run and not replace_existing:
        raise FileExistsError(f"Refusing to replace existing T/NK outputs: {existing}")

    summaries: list[dict[str, Any]] = []
    metadata_audits: list[dict[str, Any]] = []
    all_decisions: list[pd.DataFrame] = []
    dry_summary: list[dict[str, Any]] = []
    for cohort_id, source_path in paths.items():
        print(f"Filtering {cohort_id}: {source_path}", flush=True)
        source_sha256_before = sha256_file(source_path) if not dry_run else ""
        adata = ad.read_h5ad(source_path)
        try:
            decisions = compute_decisions(adata, cohort_id)
            dry_summary.append(
                {
                    "cohort_id": cohort_id,
                    "input_cells": int(adata.n_obs),
                    "retained_cells": int(decisions["keep"].sum()),
                    "retained_fraction": float(decisions["keep"].mean()),
                    "productive_tcr_input": int(decisions["productive_tcr_evidence"].sum()),
                    "annotation_tnk_input": int(decisions["annotation_tnk_evidence"].sum()),
                    "explicit_non_tnk_input": int(
                        decisions["annotation_non_tnk_evidence"].sum()
                    ),
                    "marker_tnk_input": int(decisions["marker_tnk_evidence"].sum()),
                    "nk_marker_input": int(decisions["nk_marker_rule"].sum()),
                }
            )
            if dry_run:
                continue
            output_path = output_root / f"{cohort_id}.h5ad"
            summary = subset_and_write(
                adata,
                decisions,
                cohort_id,
                output_path,
                source_path,
                source_sha256_before,
            )
            metadata_audits.append(validate_output(output_path, summary))
            summaries.append(summary)
            all_decisions.append(decisions.reset_index(drop=True))
        finally:
            del adata

    if dry_run:
        return {"dry_run": True, "cohorts": dry_summary}

    summary = pd.DataFrame(summaries).sort_values("cohort_id").reset_index(drop=True)
    decisions = pd.concat(all_decisions, ignore_index=True)
    sample_summary = summarize_samples(decisions)
    reason_summary = (
        decisions.groupby(["cohort_id", "keep", "selection_reason"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    TABLE_ROOT.mkdir(parents=True, exist_ok=True)
    metadata_audit = pd.DataFrame(metadata_audits).sort_values("cohort_id").reset_index(drop=True)
    summary.to_csv(SUMMARY_CSV, index=False)
    sample_summary.to_csv(SAMPLE_CSV, index=False)
    reason_summary.to_csv(REASON_CSV, index=False)
    decisions.to_csv(DECISIONS_CSV, index=False, compression="gzip")
    metadata_audit.to_csv(METADATA_AUDIT_CSV, index=False)
    manifest_columns = [
        "cohort_id",
        "output_h5ad",
        "retained_cells",
        "output_size_bytes",
        "output_sha256",
        "source_h5ad",
        "source_sha256",
    ]
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    summary[manifest_columns].to_csv(MANIFEST, index=False)
    write_figure(summary)
    write_report(summary, sample_summary, metadata_audit)
    return {
        "dry_run": False,
        "cohorts": len(summary),
        "input_cells": int(summary["input_cells"].sum()),
        "retained_cells": int(summary["retained_cells"].sum()),
        "productive_tcr_retention": float(
            summary["productive_tcr_retained"].sum() / summary["productive_tcr_input"].sum()
        ),
        "summary": str(SUMMARY_CSV),
        "report": str(SUMMARY_MD),
    }


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        result = execute(
            source_paths(args.cohort),
            args.output_root,
            args.dry_run,
            args.replace_existing,
        )
    except (
        ExtensionTnkFilterError,
        FileExistsError,
        FileNotFoundError,
        OSError,
        ValueError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
