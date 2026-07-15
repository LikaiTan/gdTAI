#!/usr/bin/env python3
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


import re
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd


MISSING_TEXT = {"", "nan", "na", "none", "null"}
CANONICAL_FIELDS = [
    "sample_id",
    "donor_id",
    "age",
    "sex",
    "tissue",
    "condition",
    "treatment",
    "enrichment_strategy",
    "assay_type",
    "tcr_availability",
    "original_cell_annotation",
]
RAW_LIBRARY_FIELDS = [
    "library_id",
    "sample_type",
    "donor_patient",
    "technology_simple",
    "tcr_vdj_flag",
]
GENERIC_CANDIDATES: dict[str, list[str]] = {
    "sample_id": [
        "sample_id",
        "sample",
        "Sample",
        "sample_name",
        "sample_meta",
        "sampleID",
        "orig.ident",
        "orig_ident",
        "sample_ID",
        "sample_label",
        "sample_alt_id",
        "sample_bsi_id",
        "library_name",
        "Run",
    ],
    "donor_id": [
        "donor_id",
        "donor",
        "Donor",
        "patient",
        "patient_id",
        "Subject",
        "subject",
        "subject_id",
        "patient_assignment",
        "patient_code",
        "patients",
        "patientid",
        "PatientID",
        "Patient_id",
        "donor_publicationID",
        "individual",
        "host",
        "pid",
    ],
    "age": [
        "age",
        "Age",
        "sample_age",
        "patient_age",
        "donor_age",
        "age_in_years",
        "age_years",
    ],
    "sex": [
        "sex",
        "Sex",
        "gender",
        "Gender",
        "sample_sex",
        "patient_gender",
        "donor_sex",
        "male_female",
    ],
    "tissue": [
        "tissue",
        "site",
        "organ",
        "source",
        "source_name_ch1",
        "tissue_type",
        "tissue_assignment",
        "location",
        "anatomical_site",
        "origin",
    ],
    "condition": [
        "condition",
        "disease",
        "Diagnosis",
        "status",
        "study_group",
        "phenotype",
        "clinical_status",
        "severity",
        "response",
        "disease_state",
    ],
    "treatment": [
        "treatment",
        "therapy",
        "intervention",
        "treatment_status",
        "drug",
        "anti.PD1_therapy",
        "chemotherapy",
        "targeted_therapy",
        "prior.therapy",
    ],
    "enrichment_strategy": [
        "sorting",
        "isolation",
        "purification",
        "enrichment",
        "cell_selection",
        "facs",
        "macs",
        "cell_markers",
    ],
    "assay_type": [
        "assay_type",
        "platform",
        "technology",
        "assay",
        "sequencing_method",
        "protocol",
        "chemistry",
    ],
    "tcr_availability": [
        "tcr_availability",
        "has_tcr",
        "TCRseq",
        "vdj",
        "tcr_seq",
        "tcr_data",
        "withTCR",
    ],
    "original_cell_annotation": [
        "original_cell_annotation",
        "cell_type",
        "celltype",
        "annotation",
        "author_annotation",
        "major_lineage",
        "published_celltype",
        "published_meta_celltype",
        "harmonized_celltype",
        "T_new_name",
        "major_cell_type",
        "sub_cell_type",
        "lineage",
        "ident",
        "cluster",
        "seurat_clusters",
    ],
}


def normalize_key(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(text).strip().lower())


def relative_or_str(path: Path, base: Path) -> str:
    try:
        return str(path.relative_to(base))
    except ValueError:
        return str(path)


def clean_scalar(value: Any) -> Any:
    if pd.isna(value):
        return pd.NA
    text = str(value).strip()
    if text.lower() in MISSING_TEXT:
        return pd.NA
    return text


def clean_series(series: pd.Series) -> pd.Series:
    cleaned = series.map(clean_scalar)
    return cleaned.astype("string")


def parse_bool(value: Any) -> Any:
    if pd.isna(value):
        return pd.NA
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float, np.integer, np.floating)):
        if value == 1:
            return True
        if value == 0:
            return False
    text = str(value).strip().lower()
    if text in {"true", "t", "yes", "y", "1", "present", "available", "positive", "with tcr/vdj"}:
        return True
    if text in {"false", "f", "no", "n", "0", "absent", "none", "negative", "without tcr/vdj"}:
        return False
    return pd.NA


def bool_series(series: pd.Series) -> pd.Series:
    parsed = series.map(parse_bool)
    return parsed.astype("boolean")


def find_matching_columns(obs_columns: list[str], candidates: list[str]) -> list[str]:
    norm_to_cols: dict[str, list[str]] = {}
    for col in obs_columns:
        norm_to_cols.setdefault(normalize_key(col), []).append(col)

    matched: list[str] = []
    for candidate in candidates:
        key = normalize_key(candidate)
        for col in norm_to_cols.get(key, []):
            if col not in matched:
                matched.append(col)
    return matched


def coalesce_columns(obs: pd.DataFrame, columns: list[str]) -> pd.Series:
    result = pd.Series(pd.NA, index=obs.index, dtype="string")
    for col in columns:
        if col not in obs.columns:
            continue
        cleaned = clean_series(obs[col])
        mask = result.isna() & cleaned.notna()
        if mask.any():
            result.loc[mask] = cleaned.loc[mask]
    return result


def extract_gsm(value: Any) -> str:
    if pd.isna(value):
        return ""
    m = re.search(r"(GSM\d+)", str(value), flags=re.IGNORECASE)
    return m.group(1).upper() if m else ""


def best_library_series(obs: pd.DataFrame) -> pd.Series:
    candidates = [
        "library_id",
        "gsm",
        "sample",
        "sample_id",
        "sample_name",
        "sample_meta",
        "Sample",
        "orig.ident",
        "orig_ident",
    ]
    best = pd.Series([""] * len(obs), index=obs.index, dtype=object)
    best_n = 0
    for col in candidates:
        if col not in obs.columns:
            continue
        values = obs[col].astype(str).map(extract_gsm)
        n = int((values != "").sum())
        if n > best_n:
            best = values
            best_n = n
    if best_n == 0:
        values = pd.Series(obs.index.astype(str), index=obs.index).map(extract_gsm)
        if int((values != "").sum()) > 0:
            best = values
    return clean_series(best)


def normalize_identifier_series(series: pd.Series) -> pd.Series:
    return clean_series(series)


def clean_donor_patient_value(value: Any) -> Any:
    text = clean_scalar(value)
    if pd.isna(text):
        return pd.NA
    assert isinstance(text, str)
    prefix_match = re.match(r"^(subject|patient|donor|individual|host)\s*:\s*(.+)$", text, flags=re.IGNORECASE)
    if prefix_match:
        cleaned = prefix_match.group(2).strip()
        return cleaned if cleaned else pd.NA
    if " " in text and not re.search(r"[A-Za-z]+\d+|\d+[A-Za-z]+|[_-]", text):
        return pd.NA
    return text


def clean_donor_patient_series(series: pd.Series) -> pd.Series:
    return series.map(clean_donor_patient_value).astype("string")


def parse_cell_prefix(series: pd.Series, pattern: str) -> pd.Series:
    extracted = series.astype(str).str.extract(pattern, expand=False)
    return clean_series(extracted)


def first_token_before_dash(series: pd.Series) -> pd.Series:
    out = series.astype(str).str.extract(r"^([^-]+)", expand=False)
    return clean_series(out)


def parse_regex_token(series: pd.Series, pattern: str) -> pd.Series:
    extracted = series.astype(str).str.extract(pattern, flags=re.IGNORECASE, expand=False)
    return clean_series(extracted)


def parse_sample_type_tissue(series: pd.Series) -> pd.Series:
    result = pd.Series(pd.NA, index=series.index, dtype="string")
    texts = clean_series(series)
    lower = texts.fillna("").str.lower()

    mapping = [
        (r"blood cell culture", "Blood cell culture"),
        (r"peripheral blood mononuclear cells|^pbmc$|pbmcs|human peripheral blood", "PBMC"),
        (r"blood immune cell", "Blood"),
        (r"heart tissue|heart", "Heart"),
        (r"rectum", "Rectum"),
        (r"colon", "Colon"),
        (r"lung", "Lung"),
        (r"tumor", "Tumor"),
        (r"blood", "Blood"),
    ]
    for pattern, label in mapping:
        mask = result.isna() & lower.str.contains(pattern, regex=True, na=False)
        if mask.any():
            result.loc[mask] = label
    return result


def parse_sample_type_enrichment(series: pd.Series) -> pd.Series:
    result = pd.Series(pd.NA, index=series.index, dtype="string")
    texts = clean_series(series)
    lower = texts.fillna("").str.lower()

    mapping = [
        (r"gdtcr|γδ|gd tcr", "gdTCR+ enriched"),
        (r"cd3\+", "CD3+ enriched"),
        (r"cd4\+", "CD4+ enriched"),
        (r"cd8\+", "CD8+ enriched"),
    ]
    for pattern, label in mapping:
        mask = result.isna() & lower.str.contains(pattern, regex=True, na=False)
        if mask.any():
            result.loc[mask] = label
    return result


def find_tcr_detail_columns(obs_columns: list[str]) -> list[str]:
    patterns = (
        "tcr",
        "cdr3",
        "clonotype",
        "vdj",
        "tra_",
        "trb_",
        "trav",
        "trbv",
        "traj",
        "trbj",
    )
    out: list[str] = []
    for col in obs_columns:
        lower = col.lower()
        if any(pattern in lower for pattern in patterns):
            out.append(col)
    return out


def note_source(source_map: dict[str, list[str]], field: str, note: str) -> None:
    notes = source_map.setdefault(field, [])
    if note and note not in notes:
        notes.append(note)


def fill_missing(out: pd.DataFrame, field: str, values: pd.Series, source_map: dict[str, list[str]], note: str) -> None:
    cleaned = normalize_identifier_series(values)
    mask = out[field].isna() & cleaned.notna()
    if mask.any():
        out.loc[mask, field] = cleaned.loc[mask]
        note_source(source_map, field, note)


def replace_field(out: pd.DataFrame, field: str, values: pd.Series, source_map: dict[str, list[str]], note: str) -> None:
    cleaned = normalize_identifier_series(values)
    out[field] = cleaned
    note_source(source_map, field, note)


def fill_missing_bool(out: pd.DataFrame, field: str, values: pd.Series, source_map: dict[str, list[str]], note: str) -> None:
    cleaned = bool_series(values)
    mask = out[field].isna() & cleaned.notna()
    if mask.any():
        out.loc[mask, field] = cleaned.loc[mask]
        note_source(source_map, field, note)


def replace_bool(out: pd.DataFrame, field: str, values: pd.Series, source_map: dict[str, list[str]], note: str) -> None:
    out[field] = bool_series(values)
    note_source(source_map, field, note)


def load_h5ad_registry(workspace_root: Path) -> pd.DataFrame:
    registry = pd.read_csv(workspace_root / "configs" / "datasets" / "integration_inputs.csv")
    registry["gse_id"] = registry["gse_id"].astype(str).str.upper()
    return registry


def load_sample_info(workspace_root: Path) -> pd.DataFrame:
    candidates = [
        workspace_root / "configs" / "datasets" / "sample_information_final_full.csv",
        workspace_root / "analysis_26GSE_V4" / "reports" / "sample_information_final_full.csv",
    ]
    info_path = next((path for path in candidates if path.exists()), None)
    if info_path is None:
        raise FileNotFoundError("sample_information_final_full.csv not found")

    info = pd.read_csv(info_path, low_memory=False)
    info = info.rename(columns={"gse": "gse_id"})
    info["gse_id"] = info["gse_id"].astype(str).str.upper()
    info["library_id"] = info["library_id"].astype(str).str.upper()
    for col in ["sample_type", "donor_patient", "technology_simple", "tcr_vdj_flag"]:
        info[col] = clean_series(info[col])
    keep_cols = ["gse_id", "library_id", "sample_type", "donor_patient", "technology_simple", "tcr_vdj_flag"]
    info = info[keep_cols].drop_duplicates(subset=["gse_id", "library_id"], keep="first").copy()
    return info


def apply_dataset_overrides(gse_id: str, obs: pd.DataFrame, out: pd.DataFrame, source_map: dict[str, list[str]]) -> None:
    cell_ids = pd.Series(out["cell_id"], index=out.index, dtype="string")

    if gse_id == "GSE161918":
        sample = coalesce_columns(obs, ["sample_id", "sample_name", "Sample"])
        replace_field(out, "sample_id", sample, source_map, "dataset_override:obs.sample_id|obs.sample_name|obs.Sample")
        donor = coalesce_columns(obs, ["Subject", "Donor"])
        replace_field(out, "donor_id", donor, source_map, "dataset_override:obs.Subject|obs.Donor")
        parsed_donor = parse_regex_token(sample, r"((?:HGR|AA)\d+)")
        fill_missing(out, "donor_id", parsed_donor, source_map, "dataset_override:parse_donor_from_sample_fields")
        fill_missing(out, "age", coalesce_columns(obs, ["Age"]), source_map, "dataset_override:obs.Age")
        fill_missing(out, "condition", coalesce_columns(obs, ["condition"]), source_map, "dataset_override:obs.condition")

    elif gse_id == "GSE212217":
        replace_field(out, "sample_id", coalesce_columns(obs, ["sample", "orig.ident"]), source_map, "dataset_override:obs.sample|obs.orig.ident")
        replace_field(out, "donor_id", coalesce_columns(obs, ["patient"]), source_map, "dataset_override:obs.patient")
        parsed_donor = parse_regex_token(out["sample_id"], r"pem(\d+)")
        fill_missing(out, "donor_id", parsed_donor, source_map, "dataset_override:parse_patient_from_sample_id")

    elif gse_id == "GSE168859":
        replace_field(out, "sample_id", coalesce_columns(obs, ["sample_id"]), source_map, "dataset_override:obs.sample_id")
        replace_field(out, "donor_id", clean_donor_patient_series(out["donor_patient"]), source_map, "dataset_override:clean(donor_patient)")

    elif gse_id == "GSE168163":
        parsed = parse_cell_prefix(cell_ids, r"^(S\d+)")
        replace_field(out, "sample_id", parsed, source_map, "dataset_override:cell_id_prefix")
        replace_field(out, "donor_id", parsed, source_map, "dataset_override:sample_id_as_donor")

    elif gse_id == "GSE240361":
        sample = coalesce_columns(obs, ["sample_id"])
        replace_field(out, "sample_id", sample, source_map, "dataset_override:obs.sample_id")
        replace_field(out, "donor_id", sample, source_map, "dataset_override:sample_id_as_donor")

    elif gse_id == "GSE241783":
        sample = coalesce_columns(obs, ["tcr_sample_inferred"])
        replace_field(out, "sample_id", sample, source_map, "dataset_override:obs.tcr_sample_inferred")
        fill_missing(out, "sample_id", out["library_id"], source_map, "dataset_override:library_id_fallback")
        replace_field(out, "donor_id", out["sample_id"], source_map, "dataset_override:sample_id_as_donor")

    elif gse_id == "GSE228597":
        sample = coalesce_columns(obs, ["sample_id"])
        replace_field(out, "sample_id", sample, source_map, "dataset_override:obs.sample_id")
        replace_field(out, "donor_id", sample, source_map, "dataset_override:sample_id_as_donor")

    elif gse_id == "GSE243013":
        sample = coalesce_columns(obs, ["sampleID"])
        replace_field(out, "sample_id", sample, source_map, "dataset_override:obs.sampleID")
        fill_missing(out, "sample_id", first_token_before_dash(cell_ids), source_map, "dataset_override:cell_id_prefix")
        replace_field(out, "donor_id", out["sample_id"], source_map, "dataset_override:sample_id_as_donor")

    elif gse_id in {"GSE155222", "GSE155223"}:
        replace_field(out, "donor_id", pd.Series(pd.NA, index=out.index, dtype="string"), source_map, "dataset_override:clear_generic_donor")
        fill_missing(out, "sample_id", out["library_id"], source_map, "dataset_override:library_id_fallback")


def build_harmonized_tables(root_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    workspace_root = root_dir.parent
    registry = load_h5ad_registry(workspace_root)
    sample_info = load_sample_info(workspace_root)

    harmonized_parts: list[pd.DataFrame] = []
    mapping_rows: list[dict[str, Any]] = []
    join_rows: list[dict[str, Any]] = []

    for reg_row in registry.itertuples(index=False):
        gse_id = str(reg_row.gse_id).upper()
        h5ad_path = Path(reg_row.h5ad_path)
        source_root = Path(reg_row.source_root)

        adata = ad.read_h5ad(h5ad_path, backed="r")
        obs = adata.obs.copy()
        obs_columns = list(obs.columns)

        out = pd.DataFrame(index=obs.index)
        out["cell_id"] = obs.index.astype("string")
        out["gse_id"] = gse_id
        out["source_h5ad"] = relative_or_str(h5ad_path, workspace_root)
        out["source_root"] = relative_or_str(source_root, workspace_root)

        source_map: dict[str, list[str]] = {"gse_id": ["registry:h5ad.csv"]}

        join_library_id = best_library_series(obs)
        info_sub = sample_info[sample_info["gse_id"] == gse_id].copy()
        info_rows = len(info_sub)
        info_map = info_sub.set_index("library_id").to_dict(orient="index") if info_rows else {}
        matched = join_library_id.map(lambda x: x in info_map if pd.notna(x) else False).astype("boolean")

        out["sample_info_join_library_id"] = join_library_id
        out["sample_info_matched"] = matched

        note_source(source_map, "library_id", "sample_info_join:best_library_series")
        out["library_id"] = join_library_id

        for raw_field in ["sample_type", "donor_patient", "technology_simple", "tcr_vdj_flag"]:
            if info_rows:
                mapped = join_library_id.map(lambda x: info_map.get(x, {}).get(raw_field, pd.NA) if pd.notna(x) else pd.NA)
                out[raw_field] = normalize_identifier_series(mapped)
            else:
                out[raw_field] = pd.Series(pd.NA, index=out.index, dtype="string")

        fill_missing(out, "library_id", coalesce_columns(obs, ["library_id", "gsm"]), source_map, "obs:library_id|gsm")
        fill_missing(out, "sample_type", coalesce_columns(obs, ["sample_type"]), source_map, "obs:sample_type")
        fill_missing(out, "donor_patient", coalesce_columns(obs, ["donor_patient"]), source_map, "obs:donor_patient")
        fill_missing(out, "technology_simple", coalesce_columns(obs, ["technology_simple"]), source_map, "obs:technology_simple")
        fill_missing(out, "tcr_vdj_flag", coalesce_columns(obs, ["tcr_vdj_flag"]), source_map, "obs:tcr_vdj_flag")
        if info_rows:
            for raw_field in RAW_LIBRARY_FIELDS[1:]:
                note_source(source_map, raw_field, "sample_info_join")

        for field in CANONICAL_FIELDS:
            dtype = "boolean" if field == "tcr_availability" else "string"
            out[field] = pd.Series(pd.NA, index=out.index, dtype=dtype)

        for field, candidates in GENERIC_CANDIDATES.items():
            matched_cols = find_matching_columns(obs_columns, candidates)
            if not matched_cols:
                continue
            if field == "tcr_availability":
                fill_missing_bool(out, field, coalesce_columns(obs, matched_cols), source_map, f"obs:{'|'.join(matched_cols)}")
            else:
                fill_missing(out, field, coalesce_columns(obs, matched_cols), source_map, f"obs:{'|'.join(matched_cols)}")

        fill_missing(out, "assay_type", out["technology_simple"], source_map, "sample_info:technology_simple")
        fill_missing(out, "tissue", parse_sample_type_tissue(out["sample_type"]), source_map, "parsed:sample_type->tissue")
        fill_missing(out, "enrichment_strategy", parse_sample_type_enrichment(out["sample_type"]), source_map, "parsed:sample_type->enrichment")
        fill_missing(out, "donor_id", clean_donor_patient_series(out["donor_patient"]), source_map, "clean:donor_patient")

        if "TCRseq" in obs.columns:
            fill_missing_bool(out, "tcr_availability", obs["TCRseq"], source_map, "obs:TCRseq")
        tcr_cols = [c for c in obs.columns if "tcr_" in c.lower() or c.lower().startswith("tra_") or c.lower().startswith("trb_")]
        if tcr_cols:
            has_any_tcr = obs[tcr_cols].notna().any(axis=1)
            fill_missing_bool(out, "tcr_availability", has_any_tcr, source_map, "derived:tcr_chain_presence")
        fill_missing_bool(out, "tcr_availability", out["tcr_vdj_flag"], source_map, "sample_info:tcr_vdj_flag")

        apply_dataset_overrides(gse_id, obs, out, source_map)

        fill_missing(out, "sample_id", out["library_id"], source_map, "fallback:library_id_as_sample")
        fill_missing(out, "donor_id", out["sample_id"], source_map, "fallback:sample_id_as_donor")

        out["tcr_availability"] = bool_series(out["tcr_availability"])

        tcr_detail_cols = find_tcr_detail_columns(obs_columns)
        for col in tcr_detail_cols:
            out[col] = clean_series(obs[col])

        harmonized_parts.append(out.reset_index(drop=True))

        for field in ["gse_id"] + RAW_LIBRARY_FIELDS + ["sample_info_join_library_id", "sample_info_matched"] + CANONICAL_FIELDS:
            mapping_rows.append(
                {
                    "gse_id": gse_id,
                    "source_h5ad": relative_or_str(h5ad_path, workspace_root),
                    "canonical_field": field,
                    "source_column": ";".join(source_map.get(field, [])),
                }
            )

        cells_with_lib = int(out["sample_info_join_library_id"].notna().sum())
        cells_matched = int(out["sample_info_matched"].fillna(False).sum())
        join_rows.append(
            {
                "gse_id": gse_id,
                "source_h5ad": relative_or_str(h5ad_path, workspace_root),
                "cells": int(adata.n_obs),
                "info_rows_for_gse": info_rows,
                "library_ids_in_info": int(info_sub["library_id"].nunique()) if info_rows else 0,
                "cells_with_library_id": cells_with_lib,
                "cells_matched": cells_matched,
                "match_fraction": float(cells_matched / adata.n_obs) if adata.n_obs else 0.0,
            }
        )

        adata.file.close()

    harmonized = pd.concat(harmonized_parts, axis=0, ignore_index=True)
    mapping_df = pd.DataFrame(mapping_rows)
    join_df = pd.DataFrame(join_rows).sort_values(["gse_id", "source_h5ad"]).reset_index(drop=True)

    first_columns = [
        "cell_id",
        "gse_id",
        "source_h5ad",
        "source_root",
        "library_id",
        "sample_type",
        "donor_patient",
        "technology_simple",
        "tcr_vdj_flag",
        "sample_info_join_library_id",
        "sample_info_matched",
    ] + CANONICAL_FIELDS
    other_columns = [col for col in harmonized.columns if col not in first_columns]
    harmonized = harmonized[first_columns + other_columns]
    return harmonized, mapping_df, join_df


def downsample_per_h5ad(df: pd.DataFrame, n_per_h5ad: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    sampled_parts: list[pd.DataFrame] = []
    for source_h5ad, part in df.groupby("source_h5ad", sort=True):
        n_take = min(n_per_h5ad, len(part))
        take_idx = rng.choice(part.index.to_numpy(), size=n_take, replace=False)
        sampled_parts.append(df.loc[take_idx])
    return pd.concat(sampled_parts, axis=0, ignore_index=True)


def main() -> None:
    root_dir = Path(__file__).resolve().parents[2]
    outputs_dir = root_dir / "outputs"
    reports_dir = root_dir / "reports"
    outputs_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    harmonized, mapping_df, join_df = build_harmonized_tables(root_dir)
    downsampled = downsample_per_h5ad(harmonized, n_per_h5ad=50, seed=42)

    harmonized_path = outputs_dir / "harmonized_metadata_v4.csv"
    downsampled_path = outputs_dir / "total_meta_downsample_v4_50_per_h5ad.csv"
    mapping_path = reports_dir / "metadata_harmonization_column_mapping_v4.csv"
    join_report_path = reports_dir / "sample_information_merge_report_v4.csv"

    harmonized.to_csv(harmonized_path, index=False)
    downsampled.to_csv(downsampled_path, index=False)
    mapping_df.to_csv(mapping_path, index=False)
    join_df.to_csv(join_report_path, index=False)

    print(f"Wrote: {harmonized_path}")
    print(f"Wrote: {downsampled_path}")
    print(f"Wrote: {mapping_path}")
    print(f"Wrote: {join_report_path}")
    print(f"Total rows in harmonized metadata: {len(harmonized)}")
    print(f"Total rows in downsampled metadata: {len(downsampled)}")
    print(f"Source H5ADs in downsampled metadata: {downsampled['source_h5ad'].nunique()}")


if __name__ == "__main__":
    main()
