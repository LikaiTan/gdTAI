#!/usr/bin/env python3
"""Audit and apply harmonized healthy/disease status correction.

This workflow mirrors the tissue-correction pattern:
1. Sample up to N cells per GSE from the large integrated milestone and write
   a compact audit package for review.
2. Apply a deterministic rule set to the full integrated object and create
   `obs["disease_status_corrected"]` without touching the expression matrix.

The final output column is intentionally coarse:
- healthy
- disease
- unknown
"""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import pandas as pd
from anndata.io import write_elem


REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_INTEGRATED = REPO_ROOT / "high_speed_temp/Integrated_dataset/integrated.h5ad"
DEFAULT_RULES = REPO_ROOT / "disease_status_correction_rules.json"
DEFAULT_METADATA = [
    REPO_ROOT / "analysis_26GSE_V4/outputs/harmonized_metadata_v4.csv",
    REPO_ROOT / "analysis_26GSE_V4/outputs/harmonized_metadata_supp.csv",
]
DEFAULT_TABLE_DIR = REPO_ROOT / "Integrated_dataset/tables/disease_status_correction"
OBS_COLUMNS = [
    "metadata_key",
    "original_cell_id",
    "source_gse_id",
    "project name",
    "sampleid",
    "barcodes",
    "sample_id",
    "library_id",
    "technology_simple",
]
METADATA_COLUMNS = [
    "project name",
    "sampleid",
    "barcodes",
    "obs_name",
    "source_gse_id",
    "original_cell_id",
    "sample_id",
    "library_id",
    "technology_simple",
    "condition",
    "sample_type",
    "donor_patient",
]


@dataclass
class Resolution:
    """One resolved disease-status value plus provenance."""

    value: str
    source: str
    reason: str


def configure_logging(verbose: bool) -> None:
    """Configure concise console logging."""
    level = logging.INFO if verbose else logging.WARNING
    logging.basicConfig(level=level, format="%(asctime)s | %(levelname)s | %(message)s")


def load_rules(path: Path) -> dict:
    """Load the JSON rule configuration."""
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def clean_text(value: object) -> str:
    """Normalize missing-like values to empty strings."""
    if value is None:
        return ""
    if isinstance(value, float) and np.isnan(value):
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "null", "na"}:
        return ""
    return text


def normalize_match_text(text: str) -> str:
    """Create a comparison-friendly text representation."""
    text = clean_text(text).lower()
    text = re.sub(r"([a-z])([0-9])", r"\1 \2", text)
    text = re.sub(r"([0-9])([a-z])", r"\1 \2", text)
    text = re.sub(r"[_/|+-]+", " ", text)
    text = re.sub(r"[\(\)\[\],;:]+", " ", text)
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def looks_suspicious(text: str, rules: dict) -> bool:
    """Flag obvious non-status placeholders such as matrix file names."""
    normalized = clean_text(text).lower()
    return any(token.lower() in normalized for token in rules["ignore_if_contains"])


def read_string_vector(dataset: h5py.Dataset, indices: np.ndarray | None = None) -> np.ndarray:
    """Read a string-array dataset as a NumPy object array."""
    if indices is None:
        values = dataset.asstr()[:]
    else:
        values = dataset.asstr()[indices]
    return np.asarray(values, dtype=object)


def load_obs_frame(h5ad_path: Path, columns: Iterable[str], indices: np.ndarray | None = None) -> pd.DataFrame:
    """Read selected obs columns from the integrated H5AD without loading X."""
    with h5py.File(h5ad_path, "r") as handle:
        obs = handle["obs"]
        index_values = read_string_vector(obs["_index"], indices)
        data = {"obs_name": index_values}
        for column in columns:
            data[column] = read_string_vector(obs[column], indices)
    frame = pd.DataFrame(data)
    for column in frame.columns:
        frame[column] = frame[column].map(clean_text)
    return frame


def build_metadata_key(source_gse_id: str, original_cell_id: str, barcodes: str) -> str:
    """Build the canonical metadata key used in the integrated milestone."""
    barcode = clean_text(original_cell_id) or clean_text(barcodes)
    gse = clean_text(source_gse_id)
    if not gse or not barcode:
        return ""
    return f"{gse}||{barcode}"


def load_harmonized_metadata(paths: list[Path]) -> pd.DataFrame:
    """Load the metadata columns needed for disease-status correction."""
    frames: list[pd.DataFrame] = []
    for path in paths:
        if not path.exists():
            logging.warning("Skipping missing metadata file: %s", path)
            continue
        logging.info("Loading metadata: %s", path)
        frame = pd.read_csv(
            path,
            usecols=METADATA_COLUMNS,
            dtype=str,
            keep_default_na=False,
            low_memory=False,
        )
        for column in frame.columns:
            frame[column] = frame[column].map(clean_text)
        frame["metadata_key"] = [
            build_metadata_key(gse, original, barcode)
            for gse, original, barcode in zip(
                frame["source_gse_id"], frame["original_cell_id"], frame["barcodes"], strict=False
            )
        ]
        frames.append(frame)
    if not frames:
        raise FileNotFoundError("No harmonized metadata files could be loaded.")
    metadata = pd.concat(frames, ignore_index=True)
    metadata = metadata[metadata["metadata_key"].ne("")]
    metadata = metadata.drop_duplicates(subset="metadata_key", keep="first")
    keep_columns = [
        "metadata_key",
        "source_gse_id",
        "project name",
        "sampleid",
        "barcodes",
        "obs_name",
        "sample_id",
        "library_id",
        "technology_simple",
        "condition",
        "sample_type",
        "donor_patient",
    ]
    return metadata[keep_columns].copy()


def sample_indices_by_gse(gse_values: np.ndarray, sample_per_gse: int, seed: int) -> np.ndarray:
    """Sample up to N cells per GSE."""
    rng = np.random.default_rng(seed)
    series = pd.Series(gse_values)
    sampled: list[np.ndarray] = []
    for _, index in series.groupby(series).groups.items():
        index_array = np.asarray(index, dtype=np.int64)
        if len(index_array) <= sample_per_gse:
            sampled.append(index_array)
        else:
            sampled.append(np.sort(rng.choice(index_array, size=sample_per_gse, replace=False)))
    sampled_indices = np.concatenate(sampled)
    sampled_indices.sort()
    return sampled_indices


def parse_series_matrix(path: Path) -> dict[str, dict[str, str]]:
    """Parse one GEO series_matrix file into lookup maps keyed by GSM and title."""
    fields: dict[str, list[str]] = {}
    characteristic_rows: list[list[str]] = []
    try:
        handle = gzip.open(path, "rt", errors="ignore") if path.suffix == ".gz" else open(path, "rt", errors="ignore")
        with handle as stream:
            for raw_line in stream:
                if not raw_line.startswith("!Sample_"):
                    continue
                parts = raw_line.rstrip("\n").split("\t")
                key = parts[0][8:]
                values = [part.strip().strip('"') for part in parts[1:]]
                if key == "characteristics_ch1":
                    characteristic_rows.append(values)
                else:
                    fields[key] = values
    except gzip.BadGzipFile:
        with open(path, "rt", errors="ignore") as stream:
            for raw_line in stream:
                if not raw_line.startswith("!Sample_"):
                    continue
                parts = raw_line.rstrip("\n").split("\t")
                key = parts[0][8:]
                values = [part.strip().strip('"') for part in parts[1:]]
                if key == "characteristics_ch1":
                    characteristic_rows.append(values)
                else:
                    fields[key] = values

    sample_count = 0
    if fields:
        sample_count = max(len(values) for values in fields.values())
    if characteristic_rows:
        sample_count = max(sample_count, max(len(values) for values in characteristic_rows))

    records = [
        {"geo_accession": "", "title": "", "source_name": "", "characteristics": []}
        for _ in range(sample_count)
    ]
    for key, values in fields.items():
        base = key.replace("_ch1", "")
        if base == "source_name":
            base = "source_name"
        for idx, value in enumerate(values):
            if idx < sample_count and base in {"geo_accession", "title", "source_name"}:
                records[idx][base] = clean_text(value)
    for values in characteristic_rows:
        for idx, value in enumerate(values):
            if idx < sample_count and clean_text(value):
                records[idx]["characteristics"].append(clean_text(value))

    by_gsm: dict[str, str] = {}
    by_title: dict[str, str] = {}
    for record in records:
        bundle_parts = [record["title"], record["source_name"], *record["characteristics"]]
        bundle = " | ".join(part for part in bundle_parts if clean_text(part))
        gsm = clean_text(record["geo_accession"])
        title = normalize_match_text(record["title"])
        if gsm and bundle:
            by_gsm[gsm] = bundle
        if title and bundle:
            by_title[title] = bundle
    return {"by_gsm": by_gsm, "by_title": by_title}


def discover_series_matrix_paths(gse_ids: Iterable[str]) -> dict[str, list[Path]]:
    """Auto-discover series_matrix files in common local GEO layouts."""
    found: dict[str, list[Path]] = {}
    for gse_id in sorted({clean_text(gse) for gse in gse_ids if clean_text(gse)}):
        base_dir = REPO_ROOT / "downloads" / gse_id
        if not base_dir.exists():
            continue
        paths = sorted(base_dir.glob("*series_matrix.txt.gz"))
        paths += sorted((base_dir / "matrix").glob("*series_matrix.txt.gz"))
        paths += sorted((base_dir / "metadata").glob("*series_matrix.txt"))
        paths += sorted((base_dir / "metadata").glob("*series_matrix.txt.gz"))
        deduped: list[Path] = []
        seen: set[Path] = set()
        for path in paths:
            resolved = path.resolve()
            if resolved not in seen:
                deduped.append(path)
                seen.add(resolved)
        paths = deduped
        if paths:
            found[gse_id] = paths
    return found


def build_series_matrix_maps(gse_ids: Iterable[str]) -> dict[str, dict[str, dict[str, str]]]:
    """Build GEO series_matrix lookup maps for all discoverable GSEs."""
    maps: dict[str, dict[str, dict[str, str]]] = {}
    for gse_id, paths in discover_series_matrix_paths(gse_ids).items():
        merged = {"by_gsm": {}, "by_title": {}}
        for path in paths:
            logging.info("Parsing series_matrix for %s", gse_id)
            parsed = parse_series_matrix(path)
            merged["by_gsm"].update(parsed["by_gsm"])
            merged["by_title"].update(parsed["by_title"])
        if merged["by_gsm"] or merged["by_title"]:
            maps[gse_id] = merged
    return maps


def resolve_series_bundle(row: pd.Series, series_maps: dict[str, dict[str, dict[str, str]]]) -> str:
    """Resolve the best-matching series_matrix bundle for one row."""
    gse_id = clean_text(row.get("source_gse_id"))
    if gse_id not in series_maps:
        return ""
    mapping = series_maps[gse_id]
    for field in ("library_id", "sample_id", "sampleid", "project name"):
        value = clean_text(row.get(field))
        if not value:
            continue
        gsm_match = re.search(r"GSM\d+", value, flags=re.I)
        if gsm_match:
            gsm = gsm_match.group(0).upper()
            if gsm in mapping["by_gsm"]:
                return mapping["by_gsm"][gsm]
        normalized = normalize_match_text(value)
        if normalized in mapping["by_title"]:
            return mapping["by_title"][normalized]
        for title_key, bundle in mapping["by_title"].items():
            if normalized and (normalized in title_key or title_key in normalized):
                return bundle
    return ""


def apply_regex_rules(text_values: dict[str, str], ruleset: list[dict]) -> Resolution | None:
    """Apply ordered regex rules to the selected fields."""
    for rule in ruleset:
        haystack = " || ".join(
            clean_text(text_values.get(field)) for field in rule["fields"] if clean_text(text_values.get(field))
        )
        if haystack and re.search(rule["pattern"], haystack, flags=re.IGNORECASE):
            return Resolution(rule["value"], "rule", rule["reason"])
    return None


def apply_global_rules(text_values: dict[str, str], rules: dict) -> Resolution | None:
    """Apply the global exact and regex mappings."""
    ordered_fields = [
        "condition_meta",
        "sample_type",
        "series_bundle",
        "project_name",
        "donor_patient",
        "sample_id",
        "sampleid",
        "library_id",
    ]
    for field in ordered_fields:
        value = clean_text(text_values.get(field))
        if not value:
            continue
        if field in {"condition_meta", "sample_type", "series_bundle"} and looks_suspicious(value, rules):
            continue
        normalized = normalize_match_text(value)
        if normalized in rules["global_exact_map"]:
            return Resolution(rules["global_exact_map"][normalized], field, f"global exact match on {field}")
        for rule in rules["global_regex_rules"]:
            if re.search(rule["pattern"], normalized, flags=re.IGNORECASE):
                return Resolution(rule["value"], field, rule["reason"])
    return None


def resolve_disease_status(
    row: pd.Series,
    rules: dict,
    series_maps: dict[str, dict[str, dict[str, str]]],
) -> Resolution:
    """Resolve one harmonized healthy/disease label for one cell."""
    text_values = {
        "obs_name": clean_text(row.get("obs_name")),
        "condition_meta": clean_text(row.get("condition_meta")),
        "sample_type": clean_text(row.get("sample_type")),
        "library_id": clean_text(row.get("library_id")),
        "sample_id": clean_text(row.get("sample_id")),
        "sampleid": clean_text(row.get("sampleid")),
        "project_name": clean_text(row.get("project name")),
        "donor_patient": clean_text(row.get("donor_patient")),
    }
    text_values["series_bundle"] = resolve_series_bundle(row, series_maps)

    gse_id = clean_text(row.get("source_gse_id"))
    if gse_id in rules["gse_specific_rules"]:
        result = apply_regex_rules(text_values, rules["gse_specific_rules"][gse_id])
        if result is not None:
            return Resolution(result.value, result.source, f"{gse_id}: {result.reason}")

    result = apply_global_rules(text_values, rules)
    if result is not None:
        return result

    if gse_id in rules.get("gse_default_map", {}):
        return Resolution(rules["gse_default_map"][gse_id], "gse_default", f"{gse_id}: dataset default")

    return Resolution(rules["unknown_label"], "fallback", "no matching disease-status rule")


def attach_metadata(obs: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Align harmonized metadata to the integrated obs rows by metadata_key."""
    aligned = metadata.set_index("metadata_key").reindex(obs["metadata_key"])
    merged = obs.copy()
    merged["project_name_meta"] = aligned["project name"].to_numpy(dtype=object)
    merged["sampleid_meta"] = aligned["sampleid"].to_numpy(dtype=object)
    merged["sample_id_meta"] = aligned["sample_id"].to_numpy(dtype=object)
    merged["library_id_meta"] = aligned["library_id"].to_numpy(dtype=object)
    merged["technology_simple_meta"] = aligned["technology_simple"].to_numpy(dtype=object)
    merged["condition_meta"] = aligned["condition"].to_numpy(dtype=object)
    merged["sample_type"] = aligned["sample_type"].to_numpy(dtype=object)
    merged["donor_patient"] = aligned["donor_patient"].to_numpy(dtype=object)

    merged["project name"] = [
        clean_text(left) or clean_text(right)
        for left, right in zip(merged["project name"], merged["project_name_meta"], strict=False)
    ]
    merged["sampleid"] = [
        clean_text(left) or clean_text(right)
        for left, right in zip(merged["sampleid"], merged["sampleid_meta"], strict=False)
    ]
    merged["sample_id"] = [
        clean_text(left) or clean_text(right)
        for left, right in zip(merged["sample_id"], merged["sample_id_meta"], strict=False)
    ]
    merged["library_id"] = [
        clean_text(left) or clean_text(right)
        for left, right in zip(merged["library_id"], merged["library_id_meta"], strict=False)
    ]
    merged["technology_simple"] = [
        clean_text(left) or clean_text(right)
        for left, right in zip(merged["technology_simple"], merged["technology_simple_meta"], strict=False)
    ]
    return merged.drop(columns=["project_name_meta", "sampleid_meta", "sample_id_meta", "library_id_meta", "technology_simple_meta"])


def build_series_bundle_column(
    frame: pd.DataFrame,
    series_maps: dict[str, dict[str, dict[str, str]]],
) -> pd.Series:
    """Build series-matrix bundles using unique field values, not per-cell rows."""
    bundle = pd.Series("", index=frame.index, dtype=object)
    if not series_maps:
        return bundle
    mask = frame["source_gse_id"].isin(series_maps)
    if not mask.any():
        return bundle
    for gse_id in sorted(frame.loc[mask, "source_gse_id"].unique()):
        subset_mask = frame["source_gse_id"].eq(gse_id)
        subset = frame.loc[subset_mask, ["library_id", "sample_id", "sampleid", "project name"]].copy()
        mapping = series_maps[gse_id]
        resolved = pd.Series("", index=subset.index, dtype=object)
        for field in ["library_id", "sample_id", "sampleid", "project name"]:
            values = subset[field].map(clean_text)
            unresolved_index = resolved.index[resolved.eq("")]
            if unresolved_index.empty:
                break
            unique_values = [value for value in values.loc[unresolved_index].drop_duplicates().tolist() if value]
            field_map: dict[str, str] = {}
            for value in unique_values:
                gsm_match = re.search(r"GSM\d+", value, flags=re.I)
                if gsm_match:
                    gsm = gsm_match.group(0).upper()
                    if gsm in mapping["by_gsm"]:
                        field_map[value] = mapping["by_gsm"][gsm]
                        continue
                normalized = normalize_match_text(value)
                if normalized in mapping["by_title"]:
                    field_map[value] = mapping["by_title"][normalized]
                    continue
                for title_key, text_bundle in mapping["by_title"].items():
                    if normalized and (normalized in title_key or title_key in normalized):
                        field_map[value] = text_bundle
                        break
            resolved.loc[unresolved_index] = values.loc[unresolved_index].map(field_map).fillna("")
        bundle.loc[subset_mask] = resolved
    return bundle


def add_disease_status_vectorized(
    frame: pd.DataFrame,
    rules: dict,
    series_maps: dict[str, dict[str, dict[str, str]]],
) -> pd.DataFrame:
    """Vectorized disease-status correction for the full integrated object."""
    frame = frame.copy()
    for field in ["obs_name", "condition_meta", "sample_type", "library_id", "sample_id", "sampleid", "project name", "donor_patient"]:
        frame[field] = frame[field].map(clean_text)
    if series_maps:
        frame["series_bundle"] = build_series_bundle_column(frame, series_maps)
    else:
        frame["series_bundle"] = ""

    normalized_fields: dict[str, pd.Series] = {}
    ordered_fields = ["condition_meta", "sample_type", "series_bundle", "project name", "donor_patient", "sample_id", "sampleid", "library_id"]
    for field in ordered_fields:
        series = frame[field].astype(str)
        if field in {"condition_meta", "sample_type", "series_bundle"}:
            suspicious = series.map(lambda value: looks_suspicious(value, rules))
            series = series.mask(suspicious, "")
        normalized_fields[field] = series.map(normalize_match_text)
    normalized_fields["project_name"] = normalized_fields["project name"]

    corrected = pd.Series("", index=frame.index, dtype=object)
    source = pd.Series("", index=frame.index, dtype=object)
    reason = pd.Series("", index=frame.index, dtype=object)

    for gse_id, ruleset in rules["gse_specific_rules"].items():
        remaining_mask = frame["source_gse_id"].eq(gse_id) & corrected.eq("")
        if not remaining_mask.any():
            continue
        for rule in ruleset:
            current_index = corrected.index[remaining_mask]
            if current_index.empty:
                break
            combined = pd.Series("", index=current_index, dtype=object)
            for field in rule["fields"]:
                combined = combined.str.cat(
                    normalized_fields[field].loc[current_index].fillna(""),
                    sep=" || ",
                    na_rep="",
                )
            matches = combined.str.contains(rule["pattern"], case=False, regex=True, na=False)
            assign_index = combined.index[matches]
            if len(assign_index) == 0:
                continue
            corrected.loc[assign_index] = rule["value"]
            source.loc[assign_index] = "rule"
            reason.loc[assign_index] = f"{gse_id}: {rule['reason']}"
            remaining_mask = frame["source_gse_id"].eq(gse_id) & corrected.eq("")

    for field in ordered_fields:
        remaining_index = corrected.index[corrected.eq("")]
        if remaining_index.empty:
            break
        series = normalized_fields[field].loc[remaining_index]
        mapped = series.map(rules["global_exact_map"]).fillna("")
        assign_index = mapped.index[mapped.ne("")]
        if len(assign_index) > 0:
            corrected.loc[assign_index] = mapped.loc[assign_index]
            source.loc[assign_index] = field
            reason.loc[assign_index] = f"global exact match on {field}"
        for global_rule in rules["global_regex_rules"]:
            remaining_index = corrected.index[corrected.eq("")]
            if remaining_index.empty:
                break
            series = normalized_fields[field].loc[remaining_index]
            matches = series.str.contains(global_rule["pattern"], case=False, regex=True, na=False)
            assign_index = series.index[matches]
            if len(assign_index) == 0:
                continue
            corrected.loc[assign_index] = global_rule["value"]
            source.loc[assign_index] = field
            reason.loc[assign_index] = global_rule["reason"]

    for gse_id, default_value in rules.get("gse_default_map", {}).items():
        assign_index = corrected.index[corrected.eq("") & frame["source_gse_id"].eq(gse_id)]
        if len(assign_index) == 0:
            continue
        corrected.loc[assign_index] = default_value
        source.loc[assign_index] = "gse_default"
        reason.loc[assign_index] = f"{gse_id}: dataset default"

    frame["disease_status_corrected"] = corrected.mask(corrected.eq(""), rules["unknown_label"])
    frame["disease_status_correction_source"] = source.mask(source.eq(""), "fallback")
    frame["disease_status_correction_reason"] = reason.mask(reason.eq(""), "no matching disease-status rule")
    return frame.drop(columns=["series_bundle"])


def add_disease_status(
    frame: pd.DataFrame,
    rules: dict,
    series_maps: dict[str, dict[str, dict[str, str]]],
) -> pd.DataFrame:
    """Row-wise disease-status correction for downsampled audit tables."""
    resolutions = [resolve_disease_status(row, rules, series_maps) for _, row in frame.iterrows()]
    frame = frame.copy()
    frame["disease_status_corrected"] = [resolution.value for resolution in resolutions]
    frame["disease_status_correction_source"] = [resolution.source for resolution in resolutions]
    frame["disease_status_correction_reason"] = [resolution.reason for resolution in resolutions]
    return frame


def write_markdown_report(path: Path, title: str, lines: list[str]) -> None:
    """Write a simple English markdown report."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"# {title}\n\n")
        for line in lines:
            handle.write(f"{line}\n")


def ensure_string_array_column(h5ad_path: Path, column_name: str, values: np.ndarray) -> None:
    """Write or replace one string-array obs column in place."""
    with h5py.File(h5ad_path, "r+") as handle:
        obs = handle["obs"]
        if column_name in obs:
            del obs[column_name]
        write_elem(obs, column_name, np.asarray(values, dtype=object))
        column_order = [str(item) for item in obs.attrs["column-order"]]
        if column_name not in column_order:
            column_order.append(column_name)
            del obs.attrs["column-order"]
            obs.attrs["column-order"] = np.asarray(column_order, dtype=h5py.string_dtype(encoding="utf-8"))


def run_audit(args: argparse.Namespace) -> None:
    """Run the bounded 500-cells-per-GSE audit workflow."""
    rules = load_rules(Path(args.rules))
    metadata = load_harmonized_metadata([Path(path) for path in args.metadata])
    with h5py.File(args.integrated, "r") as handle:
        source_gse = read_string_vector(handle["obs"]["source_gse_id"])
    series_maps = build_series_matrix_maps(source_gse)
    sampled_indices = sample_indices_by_gse(source_gse, rules["sample_per_gse"], rules["random_seed"])

    sampled_obs = load_obs_frame(Path(args.integrated), OBS_COLUMNS, sampled_indices)
    sampled = attach_metadata(sampled_obs, metadata)
    sampled = add_disease_status(sampled, rules, series_maps)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    sampled.to_csv(output_dir / "disease_status_correction_downsample_500_per_gse.csv", index=False)

    by_gse = (
        sampled.groupby("source_gse_id")
        .agg(
            sampled_cells=("obs_name", "size"),
            unique_condition=("condition_meta", lambda x: x.nunique(dropna=True)),
            unique_sample_type=("sample_type", lambda x: x.nunique(dropna=True)),
            unique_corrected=("disease_status_corrected", lambda x: x.nunique(dropna=True)),
            unknown_cells=("disease_status_corrected", lambda x: int((x == rules["unknown_label"]).sum())),
        )
        .reset_index()
        .sort_values("source_gse_id")
    )
    by_gse.to_csv(output_dir / "disease_status_correction_downsample_by_gse.csv", index=False)

    distinct_rows: list[pd.DataFrame] = []
    for column in ["condition_meta", "sample_type", "disease_status_corrected"]:
        value_counts = (
            sampled.assign(_value=sampled[column].map(clean_text))
            .loc[lambda df: df["_value"].ne("")]
            .groupby(["source_gse_id", "_value"])
            .size()
            .reset_index(name="cell_n")
            .rename(columns={"_value": "value"})
        )
        value_counts["column_name"] = column
        distinct_rows.append(value_counts)
    pd.concat(distinct_rows, ignore_index=True).to_csv(
        output_dir / "disease_status_correction_distinct_values.csv", index=False
    )

    report_lines = [
        f"- Integrated input: `{args.integrated}`",
        f"- Sampled cells: `{len(sampled):,}`",
        f"- GSE count: `{sampled['source_gse_id'].nunique():,}`",
        f"- Unknown cells after correction: `{int((sampled['disease_status_corrected'] == rules['unknown_label']).sum()):,}`",
        "",
        "## Top corrected values",
    ]
    for value, count in sampled["disease_status_corrected"].value_counts().head(20).items():
        report_lines.append(f"- `{value}`: `{count:,}`")
    report_lines.extend(["", "## Per-GSE unknown counts"])
    unknown_by_gse = by_gse.loc[by_gse["unknown_cells"] > 0]
    if unknown_by_gse.empty:
        report_lines.append("- none")
    else:
        for _, row in unknown_by_gse.iterrows():
            report_lines.append(f"- `{row['source_gse_id']}`: `{int(row['unknown_cells']):,}` unknown cells in sample")
    write_markdown_report(
        output_dir / "disease_status_correction_audit.md",
        "Disease Status Correction Audit",
        report_lines,
    )


def run_apply(args: argparse.Namespace) -> None:
    """Apply disease-status correction to the full integrated H5AD and export validation outputs."""
    rules = load_rules(Path(args.rules))
    metadata = load_harmonized_metadata([Path(path) for path in args.metadata])
    obs = load_obs_frame(Path(args.integrated), OBS_COLUMNS)
    series_maps = build_series_matrix_maps(obs["source_gse_id"].unique())
    full = attach_metadata(obs, metadata)
    full = add_disease_status_vectorized(full, rules, series_maps)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    by_gse = (
        full.groupby("source_gse_id")
        .agg(
            cell_n=("obs_name", "size"),
            original_blank_condition=("condition_meta", lambda x: int(pd.Series(x).map(clean_text).eq("").sum())),
            corrected_unknown=("disease_status_corrected", lambda x: int((x == rules["unknown_label"]).sum())),
            unique_corrected=("disease_status_corrected", lambda x: x.nunique(dropna=True)),
        )
        .reset_index()
        .sort_values("source_gse_id")
    )
    by_gse.to_csv(output_dir / "disease_status_corrected_by_gse.csv", index=False)

    value_counts = (
        full["disease_status_corrected"]
        .value_counts(dropna=False)
        .rename_axis("disease_status_corrected")
        .reset_index(name="cell_n")
    )
    value_counts.to_csv(output_dir / "disease_status_corrected_value_counts.csv", index=False)

    reason_counts = (
        full["disease_status_correction_reason"]
        .value_counts(dropna=False)
        .rename_axis("reason")
        .reset_index(name="cell_n")
    )
    reason_counts.to_csv(output_dir / "disease_status_correction_reason_counts.csv", index=False)

    full[
        [
            "obs_name",
            "metadata_key",
            "source_gse_id",
            "project name",
            "sampleid",
            "sample_id",
            "library_id",
            "condition_meta",
            "sample_type",
            "donor_patient",
            "disease_status_corrected",
            "disease_status_correction_reason",
        ]
    ].to_csv(output_dir / "disease_status_corrected_column_export.csv.gz", index=False)

    full.loc[full["disease_status_corrected"].eq(rules["unknown_label"]), [
        "obs_name",
        "source_gse_id",
        "condition_meta",
        "sample_type",
        "donor_patient",
        "project name",
        "library_id",
        "sample_id",
        "sampleid",
    ]].head(5000).to_csv(output_dir / "disease_status_corrected_unknown_examples.csv", index=False)

    if args.write_h5ad:
        ensure_string_array_column(
            Path(args.integrated),
            "disease_status_corrected",
            full["disease_status_corrected"].to_numpy(dtype=object),
        )

    report_lines = [
        f"- Integrated input: `{args.integrated}`",
        f"- Cells processed: `{len(full):,}`",
        f"- Distinct corrected values: `{full['disease_status_corrected'].nunique():,}`",
        f"- Unknown cells: `{int(full['disease_status_corrected'].eq(rules['unknown_label']).sum()):,}`",
        f"- Wrote back into H5AD: `{args.write_h5ad}`",
        f"- Column export: `{output_dir / 'disease_status_corrected_column_export.csv.gz'}`",
        "",
        "## Top corrected values",
    ]
    for _, row in value_counts.iterrows():
        report_lines.append(f"- `{row['disease_status_corrected']}`: `{int(row['cell_n']):,}`")
    write_markdown_report(
        output_dir / "disease_status_correction_apply.md",
        "Disease Status Correction Apply Summary",
        report_lines,
    )


def build_parser() -> argparse.ArgumentParser:
    """Create the CLI parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--integrated", default=str(DEFAULT_INTEGRATED), help="Path to the integrated H5AD.")
    parser.add_argument("--rules", default=str(DEFAULT_RULES), help="Path to the disease-status correction rules JSON.")
    parser.add_argument(
        "--metadata",
        nargs="+",
        default=[str(path) for path in DEFAULT_METADATA],
        help="One or more harmonized metadata CSV files.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_TABLE_DIR),
        help="Directory for audit or validation tables.",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable INFO logging.")

    subparsers = parser.add_subparsers(dest="command", required=True)

    audit = subparsers.add_parser("audit", help="Run the 500-cells-per-GSE audit workflow.")
    audit.set_defaults(func=run_audit)

    apply = subparsers.add_parser("apply", help="Apply disease-status correction to the full integrated object.")
    apply.add_argument(
        "--write-h5ad",
        action="store_true",
        help="Write disease_status_corrected back into the integrated H5AD after export and validation.",
    )
    apply.set_defaults(func=run_apply)
    return parser


def main() -> None:
    """CLI entrypoint."""
    parser = build_parser()
    args = parser.parse_args()
    configure_logging(args.verbose)
    args.func(args)


if __name__ == "__main__":
    main()
