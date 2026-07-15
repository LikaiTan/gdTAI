#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Build a metadata-first gdT atlas from gdTAI v2 high-purity predictions.

This script is read-only with respect to the input H5AD. It uses the
NK-optimized gdTAI v2 prediction mask (`annotation_specific_pred`) and creates
report-only harmonized metadata columns for tissue, disease, and available age.
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
import math
import re
import subprocess
from pathlib import Path
from typing import Any, Iterable

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_H5AD = PROJECT_ROOT / "high_speed_temp/Integrated_dataset/integrated_plus6.h5ad"
DEFAULT_PREDICTIONS = (
    PROJECT_ROOT
    / "Integrated_dataset/tables/gdT_prediction/gdtai_annotation_specific_cascade/full_prediction_arrays.npz"
)
DEFAULT_RULES = PROJECT_ROOT / "configs" / "gdt_atlas" / "metadata_rules.json"
DEFAULT_METADATA = [
    PROJECT_ROOT / "analysis_26GSE_V4/outputs/harmonized_metadata_v4.csv",
    PROJECT_ROOT / "analysis_26GSE_V4/outputs/harmonized_metadata_supp.csv",
]
TABLE_DIR = PROJECT_ROOT / "Integrated_dataset/tables/gdT_atlas"
FIGURE_DIR = PROJECT_ROOT / "Integrated_dataset/figures/gdT_atlas"
LOG_DIR = PROJECT_ROOT / "Integrated_dataset/logs/gdT_atlas"
STATIC_DIR = PROJECT_ROOT / "gdT_atlas"
REPORT_HTML = STATIC_DIR / "index.html"
REPORT_PDF = STATIC_DIR / "gdT_atlas_report.pdf"
SUMMARY_MD = LOG_DIR / "gdtai_nk_optimized_gdt_atlas_summary.md"
RUN_LOG = LOG_DIR / "run.log"

EXPECTED_N_OBS = 5_128_904
EXPECTED_HIGH_PURITY_COUNT = 359_857
PREDICTION_KEY = "annotation_specific_pred"
INVALID_STRINGS = {"", "nan", "none", "na", "n/a", "<na>", "null"}
TARGET_OBS_COLUMNS = [
    "source_gse_id",
    "tissue",
    "tissue_corrected",
    "condition",
    "simple_annotation_plus6",
    "library_id",
    "sample_id",
    "sampleid",
    "donor_id",
    "has_TRA_TRB_paired",
    "has_any_ab_tcr",
    "has_any_gd_tcr",
]
METADATA_COLUMNS = [
    "source_gse_id",
    "sample_id",
    "sampleid",
    "library_id",
    "donor_id",
    "age",
    "sex",
    "tissue",
    "condition",
    "sample_type",
    "donor_patient",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build gdTAI high-purity predicted gdT atlas report.")
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--prediction-npz", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--rules", type=Path, default=DEFAULT_RULES)
    parser.add_argument("--metadata", type=Path, action="append", default=DEFAULT_METADATA)
    parser.add_argument("--sample-per-source", type=int, default=500)
    parser.add_argument("--umap-background", type=int, default=300_000)
    parser.add_argument("--umap-positive", type=int, default=120_000)
    parser.add_argument("--seed", type=int, default=20260617)
    parser.add_argument("--no-pdf", action="store_true")
    return parser.parse_args()


def ensure_dirs() -> None:
    for path in [TABLE_DIR, FIGURE_DIR, LOG_DIR, STATIC_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def setup_logging() -> None:
    ensure_dirs()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(RUN_LOG, mode="w", encoding="utf-8"), logging.StreamHandler()],
        force=True,
    )


def clean_text(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    text = str(value).strip()
    if text.lower() in INVALID_STRINGS:
        return ""
    return text


def normalize_text(value: object) -> str:
    text = clean_text(value).lower()
    text = re.sub(r"([a-z])([0-9])", r"\1 \2", text)
    text = re.sub(r"([0-9])([a-z])", r"\1 \2", text)
    text = re.sub(r"[_/|+-]+", " ", text)
    text = re.sub(r"[\(\)\[\],;:]+", " ", text)
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def decode_attr(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def read_string_dataset(dataset: h5py.Dataset) -> np.ndarray:
    try:
        return np.asarray(dataset.asstr()[:], dtype=object)
    except Exception:
        values = dataset[:]
        return np.asarray([v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in values], dtype=object)


def read_obs_column(handle: h5py.File, column: str) -> np.ndarray:
    obs = handle["obs"]
    if column not in obs:
        n = int(obs["_index"].shape[0])
        return np.full(n, "", dtype=object)
    obj = obs[column]
    if isinstance(obj, h5py.Group):
        encoding_type = decode_attr(obj.attrs.get("encoding-type", ""))
        if encoding_type == "categorical":
            categories = read_string_dataset(obj["categories"])
            codes = obj["codes"][:]
            out = np.full(codes.shape, "", dtype=object)
            valid = codes >= 0
            out[valid] = categories[codes[valid]]
            return out
        if "values" in obj and "mask" in obj:
            values = obj["values"][:]
            mask = obj["mask"][:].astype(bool)
            if values.dtype.kind in {"S", "O", "U"}:
                out = np.asarray([v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in values], dtype=object)
                out[mask] = ""
                return out
            out = np.asarray(values)
            if mask.any():
                out = out.astype(object)
                out[mask] = np.nan
            return out
        raise TypeError(f"Unsupported obs group encoding for `{column}`: {encoding_type}")
    values = obj[:]
    if values.dtype.kind in {"S", "O", "U"}:
        return read_string_dataset(obj)
    return np.asarray(values)


def read_bool_obs(handle: h5py.File, column: str) -> np.ndarray:
    if column not in handle["obs"]:
        return np.zeros(int(handle["obs"]["_index"].shape[0]), dtype=bool)
    values = read_obs_column(handle, column)
    if np.issubdtype(values.dtype, np.bool_):
        return values.astype(bool, copy=False)
    if np.issubdtype(values.dtype, np.number):
        return values.astype(float) != 0
    lowered = pd.Series(values, copy=False).astype("string").fillna("").str.strip().str.lower()
    return lowered.isin({"true", "1", "yes", "y", "t"}).to_numpy(dtype=bool)


def clean_series(values: np.ndarray | pd.Series) -> pd.Series:
    return pd.Series(values, copy=False).astype("string").fillna("").str.strip().mask(
        lambda s: s.str.lower().isin(INVALID_STRINGS), ""
    )


def load_h5ad_obs(path: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    stat_before = path.stat()
    with h5py.File(path, "r") as handle:
        n_obs = int(handle["obs"]["_index"].shape[0])
        obs_index = read_string_dataset(handle["obs"]["_index"])
        data: dict[str, Any] = {"obs_name": obs_index}
        for column in TARGET_OBS_COLUMNS:
            if column.startswith("has_"):
                data[column] = read_bool_obs(handle, column)
            else:
                data[column] = read_obs_column(handle, column)
        has_umap = "obsm" in handle and "X_umap" in handle["obsm"]
    df = pd.DataFrame(data)
    for column in df.columns:
        if df[column].dtype == object or str(df[column].dtype).startswith("string"):
            df[column] = clean_series(df[column]).astype(str)
    info = {
        "input_h5ad": str(path),
        "n_obs": n_obs,
        "has_umap": has_umap,
        "size_before": stat_before.st_size,
        "mtime_ns_before": stat_before.st_mtime_ns,
    }
    return df, info


def load_predictions(path: Path, n_obs: int) -> np.ndarray:
    with np.load(path) as arrays:
        if PREDICTION_KEY not in arrays:
            raise KeyError(f"Prediction key `{PREDICTION_KEY}` not found in {path}")
        pred = np.asarray(arrays[PREDICTION_KEY], dtype=bool)
    if pred.shape[0] != n_obs:
        raise RuntimeError(f"Prediction length {pred.shape[0]:,} does not match H5AD n_obs {n_obs:,}")
    if n_obs != EXPECTED_N_OBS:
        raise RuntimeError(f"Unexpected H5AD n_obs {n_obs:,}; expected {EXPECTED_N_OBS:,}")
    if int(pred.sum()) != EXPECTED_HIGH_PURITY_COUNT:
        raise RuntimeError(
            f"Unexpected NK-optimized prediction count {int(pred.sum()):,}; "
            f"expected {EXPECTED_HIGH_PURITY_COUNT:,}"
        )
    return pred


def mode_nonempty(series: pd.Series) -> str:
    s = series.astype("string").fillna("").str.strip()
    s = s[~s.str.lower().isin(INVALID_STRINGS)]
    if s.empty:
        return ""
    return str(s.value_counts(sort=True).index[0])


def load_metadata_lookup(paths: Iterable[Path]) -> dict[str, pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    for path in paths:
        if not path.exists():
            logging.warning("Metadata file missing: %s", path)
            continue
        logging.info("Loading metadata lookup columns from %s", path)
        frame = pd.read_csv(
            path,
            usecols=lambda col: col in METADATA_COLUMNS,
            dtype=str,
            keep_default_na=False,
            low_memory=False,
        )
        for column in frame.columns:
            frame[column] = clean_series(frame[column]).astype(str)
        if "source_gse_id" not in frame:
            continue
        frames.append(frame)
    if not frames:
        return {}
    metadata = pd.concat(frames, ignore_index=True, sort=False).fillna("")
    lookups: dict[str, pd.DataFrame] = {}
    value_columns = ["age", "sex", "tissue", "condition", "sample_type", "donor_patient"]
    for key_column in ["sample_id", "library_id", "sampleid"]:
        if key_column not in metadata.columns:
            continue
        tmp = metadata.loc[metadata[key_column].astype(str).ne(""), ["source_gse_id", key_column, *value_columns]].copy()
        if tmp.empty:
            continue
        tmp["join_key"] = tmp["source_gse_id"].astype(str) + "||" + tmp[key_column].astype(str)
        grouped = tmp.groupby("join_key", sort=False)[value_columns].agg(mode_nonempty)
        lookups[key_column] = grouped
        logging.info("Built %s metadata lookup with %s keys", key_column, f"{grouped.shape[0]:,}")
    return lookups


def apply_metadata_lookup(df: pd.DataFrame, lookups: dict[str, pd.DataFrame]) -> pd.DataFrame:
    for column in ["age", "sex", "tissue_meta", "condition_meta", "sample_type_meta", "donor_patient_meta"]:
        df[column] = ""
    lookup_to_output = {
        "age": "age",
        "sex": "sex",
        "tissue": "tissue_meta",
        "condition": "condition_meta",
        "sample_type": "sample_type_meta",
        "donor_patient": "donor_patient_meta",
    }
    for key_column in ["sample_id", "library_id", "sampleid"]:
        if key_column not in lookups or key_column not in df:
            continue
        join_key = df["source_gse_id"].astype(str) + "||" + df[key_column].astype(str)
        lookup = lookups[key_column]
        for source_column, out_column in lookup_to_output.items():
            if source_column not in lookup:
                continue
            missing = df[out_column].astype(str).eq("")
            if not missing.any():
                continue
            mapped = join_key[missing].map(lookup[source_column]).fillna("")
            df.loc[missing, out_column] = mapped.to_numpy(dtype=object)
    return df


def first_nonempty(*values: object) -> str:
    for value in values:
        text = clean_text(value)
        if text:
            return text
    return ""


def resolve_by_rules(text: str, exact_map: dict[str, str], regex_rules: list[dict[str, str]], default: str) -> str:
    norm = normalize_text(text)
    if norm in exact_map:
        return exact_map[norm]
    for rule in regex_rules:
        if re.search(rule["pattern"], norm, flags=re.IGNORECASE):
            return rule["value"]
    return default


def normalize_series(series: pd.Series) -> pd.Series:
    out = clean_series(series).str.lower()
    out = out.str.replace(r"([a-z])([0-9])", r"\1 \2", regex=True)
    out = out.str.replace(r"([0-9])([a-z])", r"\1 \2", regex=True)
    out = out.str.replace(r"[_/|+-]+", " ", regex=True)
    out = out.str.replace(r"[\(\)\[\],;:]+", " ", regex=True)
    out = out.str.replace(r"\s+", " ", regex=True).str.strip()
    return out


def resolve_series_by_rules(
    series: pd.Series,
    exact_map: dict[str, str],
    regex_rules: list[dict[str, str]],
    default: str,
) -> pd.Series:
    norm = normalize_series(series)
    out = norm.map(exact_map).fillna(default).astype(object)
    for rule in regex_rules:
        unresolved = pd.Series(out, index=series.index).eq(default)
        if not unresolved.any():
            break
        mask = unresolved & norm.str.contains(rule["pattern"], regex=True, case=False, na=False)
        if mask.any():
            out[mask] = rule["value"]
    return pd.Series(out, index=series.index).astype(str)


def resolve_from_columns(
    df: pd.DataFrame,
    columns: list[str],
    exact_map: dict[str, str],
    regex_rules: list[dict[str, str]],
    default: str,
) -> pd.Series:
    out = pd.Series(default, index=df.index, dtype=object)
    for column in columns:
        resolved = resolve_series_by_rules(df[column], exact_map, regex_rules, default)
        mask = out.eq(default) & resolved.ne(default)
        if mask.any():
            out.loc[mask] = resolved.loc[mask].to_numpy(dtype=object)
    return out.astype(str)


def concatenate_clean_columns(df: pd.DataFrame, columns: list[str]) -> pd.Series:
    out = clean_series(df[columns[0]]).astype(str)
    for column in columns[1:]:
        out = out + " " + clean_series(df[column]).astype(str)
    return out

def combination_key(df: pd.DataFrame, columns: list[str], sep: str = "\x1f") -> pd.Series:
    out = clean_series(df[columns[0]]).astype(str)
    for column in columns[1:]:
        out = out + sep + clean_series(df[column]).astype(str)
    return out


def resolve_from_columns_unique(
    df: pd.DataFrame,
    columns: list[str],
    exact_map: dict[str, str],
    regex_rules: list[dict[str, str]],
    default: str,
) -> pd.Series:
    keys = combination_key(df, columns)
    unique_mask = ~keys.duplicated()
    unique_df = df.loc[unique_mask, columns].copy()
    unique_keys = keys.loc[unique_mask]
    resolved = resolve_from_columns(unique_df, columns, exact_map, regex_rules, default)
    lookup = pd.Series(resolved.to_numpy(dtype=object), index=unique_keys.to_numpy(dtype=object))
    return keys.map(lookup).fillna(default).astype(str)


def resolve_context_unique(df: pd.DataFrame, columns: list[str], rules: dict[str, Any]) -> pd.Series:
    keys = combination_key(df, columns)
    unique_mask = ~keys.duplicated()
    unique_text = concatenate_clean_columns(df.loc[unique_mask, columns], columns)
    unique_resolved = resolve_series_by_rules(
        unique_text,
        {},
        rules["specimen_context_rules"],
        "unspecified_context",
    )
    lookup = pd.Series(unique_resolved.to_numpy(dtype=object), index=keys.loc[unique_mask].to_numpy(dtype=object))
    return keys.map(lookup).fillna("unspecified_context").astype(str)

def parse_age_value(value: str) -> float | None:
    text = normalize_text(value)
    if not text:
        return None
    if re.search(r"fetus|fetal", text):
        return -1.0
    if re.search(r"newborn|neonate|cord blood", text):
        return 0.0
    match = re.search(r"(\d+(?:\.\d+)?)", text)
    if not match:
        return None
    age = float(match.group(1))
    if re.search(r"week|wk", text):
        age = age / 52.0
    if re.search(r"month|mo", text):
        age = age / 12.0
    if age < 0 or age > 120:
        return None
    return age


def age_group(value: str, rules: dict[str, Any]) -> str:
    text = normalize_text(value)
    if re.search(r"fetus|fetal", text):
        return "fetal"
    if re.search(r"newborn|neonate|cord blood", text):
        return "newborn"
    age = parse_age_value(value)
    if age is None:
        return "unknown"
    for bin_def in rules["age_bins"]:
        if float(bin_def["min"]) <= age <= float(bin_def["max"]):
            return str(bin_def["label"])
    return "unknown"


def add_harmonized_metadata(df: pd.DataFrame, rules: dict[str, Any]) -> pd.DataFrame:
    logging.info("Resolving tissue_site_gdt_atlas")
    df["tissue_site_gdt_atlas"] = resolve_from_columns_unique(
        df,
        ["tissue_corrected", "tissue", "tissue_meta", "sample_type_meta", "condition"],
        rules["tissue_exact_map"],
        rules["tissue_regex_rules"],
        rules["unknown_label"],
    )
    logging.info("Resolving specimen_context_gdt_atlas")
    df["specimen_context_gdt_atlas"] = resolve_context_unique(
        df,
        ["tissue_corrected", "tissue", "tissue_meta", "sample_type_meta", "condition", "condition_meta"],
        rules,
    )
    blood_context = df["tissue_site_gdt_atlas"].isin({"blood", "cord_blood"}) & df[
        "specimen_context_gdt_atlas"
    ].eq("unspecified_context")
    df.loc[blood_context, "specimen_context_gdt_atlas"] = "blood"
    logging.info("Resolving disease_status_gdt_atlas")
    df["disease_status_gdt_atlas"] = resolve_from_columns_unique(
        df,
        ["condition", "condition_meta", "sample_type_meta", "donor_patient_meta", "tissue_corrected", "tissue"],
        rules["disease_exact_map"],
        rules["disease_regex_rules"],
        rules["unknown_label"],
    )
    default_mask = df["disease_status_gdt_atlas"].eq(rules["unknown_label"])
    if default_mask.any():
        default_map = rules.get("gse_default_disease", {})
        mapped = df.loc[default_mask, "source_gse_id"].map(default_map).fillna(rules["unknown_label"])
        df.loc[default_mask, "disease_status_gdt_atlas"] = mapped.to_numpy(dtype=object)
    detail = clean_series(df["condition"]).astype(str)
    for column in ["condition_meta", "sample_type_meta"]:
        missing = detail.eq("")
        if missing.any():
            detail.loc[missing] = clean_series(df.loc[missing, column]).astype(str).to_numpy(dtype=object)
    df["condition_detail_gdt_atlas"] = detail.mask(detail.eq(""), "unknown").str.replace(r"\s+", " ", regex=True)
    age_values = clean_series(df["age"]).astype(str)
    unique_age = pd.unique(age_values)
    age_year_map = {value: parse_age_value(value) for value in unique_age}
    age_group_map = {value: age_group(value, rules) for value in unique_age}
    df["age_years_gdt_atlas"] = age_values.map(age_year_map)
    df["age_group_gdt_atlas"] = age_values.map(age_group_map).fillna("unknown")
    df["sex_gdt_atlas"] = clean_series(df["sex"]).replace("", "unknown").astype(str)
    return df


def summarize_by(df: pd.DataFrame, pred_col: str, group_cols: list[str], path: Path) -> pd.DataFrame:
    grouped = (
        df.groupby(group_cols, dropna=False)
        .agg(
            total_cells=(pred_col, "size"),
            predicted_gdT=(pred_col, "sum"),
            paired_TCRAB_cells=("has_TRA_TRB_paired", "sum"),
            predicted_paired_TCRAB_FP=("predicted_paired_TCRAB_FP", "sum"),
            NK_cells=("is_NK_cell", "sum"),
            predicted_NK_FP=("predicted_NK_FP", "sum"),
        )
        .reset_index()
    )
    grouped["predicted_gdT_fraction"] = grouped["predicted_gdT"] / grouped["total_cells"].replace(0, np.nan)
    grouped["paired_TCRAB_FP_fraction_of_predictions"] = grouped["predicted_paired_TCRAB_FP"] / grouped[
        "predicted_gdT"
    ].replace(0, np.nan)
    grouped["NK_FP_fraction_of_predictions"] = grouped["predicted_NK_FP"] / grouped["predicted_gdT"].replace(0, np.nan)
    grouped = grouped.replace([np.inf, -np.inf], np.nan)
    grouped = grouped.sort_values(["predicted_gdT", "total_cells"], ascending=[False, False])
    grouped.to_csv(path, index=False)
    return grouped


def metadata_completeness(df: pd.DataFrame) -> pd.DataFrame:
    fields = {
        "raw_tissue_any": df[["tissue", "tissue_corrected", "tissue_meta", "sample_type_meta"]]
        .apply(lambda row: any(clean_text(v) for v in row), axis=1),
        "harmonized_tissue_known": df["tissue_site_gdt_atlas"].ne("unknown"),
        "disease_known": df["disease_status_gdt_atlas"].ne("unknown"),
        "condition_detail_known": df["condition_detail_gdt_atlas"].ne("unknown"),
        "age_known": df["age_group_gdt_atlas"].ne("unknown"),
        "sex_known": df["sex_gdt_atlas"].ne("unknown"),
    }
    base = pd.DataFrame({"source_gse_id": df["source_gse_id"], **fields})
    rows = []
    for source, group in base.groupby("source_gse_id", dropna=False):
        row = {"source_gse_id": source, "n_cells": int(group.shape[0])}
        for field in fields:
            row[f"{field}_n"] = int(group[field].sum())
            row[f"{field}_fraction"] = float(group[field].mean())
        rows.append(row)
    return pd.DataFrame(rows).sort_values("n_cells", ascending=False)


def sample_audit(df: pd.DataFrame, sample_per_source: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    parts = []
    for _, group in df.groupby("source_gse_id", dropna=False):
        n = min(sample_per_source, group.shape[0])
        idx = rng.choice(group.index.to_numpy(), size=n, replace=False) if group.shape[0] > n else group.index.to_numpy()
        parts.append(df.loc[idx])
    return pd.concat(parts, ignore_index=True)


def unresolved_examples(df: pd.DataFrame) -> pd.DataFrame:
    mask = (
        df["tissue_site_gdt_atlas"].eq("unknown")
        | df["disease_status_gdt_atlas"].eq("unknown")
        | df["age_group_gdt_atlas"].eq("unknown")
    )
    cols = [
        "source_gse_id",
        "tissue",
        "tissue_corrected",
        "tissue_meta",
        "sample_type_meta",
        "condition",
        "condition_meta",
        "age",
        "tissue_site_gdt_atlas",
        "disease_status_gdt_atlas",
        "age_group_gdt_atlas",
    ]
    out = (
        df.loc[mask, cols]
        .drop_duplicates()
        .groupby(cols, dropna=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values("n_cells", ascending=False)
        .head(250)
    )
    return out


def plot_bar(df: pd.DataFrame, label_col: str, value_col: str, path: Path, title: str, xlabel: str, top_n: int = 25) -> Path:
    plot = df.sort_values(value_col, ascending=False).head(top_n).copy()
    plot[label_col] = plot[label_col].astype(str)
    height = max(4.0, 0.28 * plot.shape[0] + 1.2)
    fig, ax = plt.subplots(figsize=(8.5, height), constrained_layout=True)
    y = np.arange(plot.shape[0])
    ax.barh(y, plot[value_col].to_numpy(dtype=float), color="#0f766e")
    ax.set_yticks(y, labels=plot[label_col])
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    for i, value in enumerate(plot[value_col].to_numpy(dtype=int)):
        ax.text(value, i, f" {value:,}", va="center", fontsize=8)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_heatmap(cross: pd.DataFrame, path: Path, title: str) -> Path:
    if cross.empty:
        return path
    mat = cross.pivot_table(
        index="tissue_site_gdt_atlas",
        columns="disease_status_gdt_atlas",
        values="predicted_gdT",
        aggfunc="sum",
        fill_value=0,
    )
    mat = mat.loc[mat.sum(axis=1).sort_values(ascending=False).head(20).index]
    fig, ax = plt.subplots(figsize=(7.5, max(4.0, 0.3 * mat.shape[0] + 1.8)), constrained_layout=True)
    values = np.log10(mat.to_numpy(dtype=float) + 1.0)
    im = ax.imshow(values, cmap="YlGnBu", aspect="auto")
    ax.set_xticks(np.arange(mat.shape[1]), labels=mat.columns, rotation=30, ha="right")
    ax.set_yticks(np.arange(mat.shape[0]), labels=mat.index)
    ax.set_title(title)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f"{int(mat.iat[i, j]):,}", ha="center", va="center", fontsize=7)
    fig.colorbar(im, ax=ax, label="log10(predicted gdT + 1)")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def read_umap(handle: h5py.File) -> np.ndarray | None:
    if "obsm" not in handle or "X_umap" not in handle["obsm"]:
        return None
    obj = handle["obsm"]["X_umap"]
    if isinstance(obj, h5py.Dataset):
        arr = obj[:]
    elif isinstance(obj, h5py.Group) and "data" in obj:
        arr = obj["data"][:]
    else:
        return None
    arr = np.asarray(arr)
    if arr.ndim != 2 or arr.shape[1] < 2:
        return None
    return arr[:, :2].astype(np.float32, copy=False)


def plot_umap(h5ad_path: Path, pred: np.ndarray, path: Path, background_n: int, positive_n: int, seed: int) -> Path | None:
    with h5py.File(h5ad_path, "r") as handle:
        umap = read_umap(handle)
    if umap is None:
        logging.warning("No X_umap found; skipping UMAP plot")
        return None
    rng = np.random.default_rng(seed)
    positive = np.flatnonzero(pred)
    negative = np.flatnonzero(~pred)
    pos_idx = positive if positive.size <= positive_n else rng.choice(positive, size=positive_n, replace=False)
    neg_idx = negative if negative.size <= background_n else rng.choice(negative, size=background_n, replace=False)
    fig, ax = plt.subplots(figsize=(7.2, 6.4), constrained_layout=True)
    ax.scatter(umap[neg_idx, 0], umap[neg_idx, 1], s=1, c="#c2c8cf", alpha=0.18, linewidths=0, rasterized=True)
    ax.scatter(umap[pos_idx, 0], umap[pos_idx, 1], s=2, c="#b51f2a", alpha=0.70, linewidths=0, rasterized=True)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("gdTAI v2 high-purity predicted gdT cells")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def dataframe_to_markdown(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if max_rows is not None:
        df = df.head(max_rows)
    if df.empty:
        return "_No rows._"
    view = df.copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            view[col] = view[col].map(lambda x: "" if pd.isna(x) else f"{x:.6g}")
    view = view.fillna("").astype(str)
    columns = [str(c) for c in view.columns]
    rows = view.values.tolist()
    out = ["| " + " | ".join(c.replace("|", "\\|") for c in columns) + " |"]
    out.append("| " + " | ".join("---" for _ in columns) + " |")
    for row in rows:
        out.append("| " + " | ".join(str(v).replace("|", "\\|").replace("\n", " ") for v in row) + " |")
    return "\n".join(out)


def dataframe_to_html(df: pd.DataFrame, max_rows: int = 30) -> str:
    if df.empty:
        return "<p>No rows.</p>"
    view = df.head(max_rows).copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            view[col] = view[col].map(lambda x: "" if pd.isna(x) else f"{x:.4g}")
    return view.to_html(index=False, escape=True, classes="dataframe")


def rel(path: Path) -> str:
    return html.escape(str(path.relative_to(PROJECT_ROOT)))


def write_report(
    *,
    overall: pd.DataFrame,
    source: pd.DataFrame,
    tissue: pd.DataFrame,
    disease: pd.DataFrame,
    age: pd.DataFrame,
    completeness: pd.DataFrame,
    unresolved: pd.DataFrame,
    figures: list[Path],
    h5ad_info: dict[str, Any],
    mtime_after: int,
    size_after: int,
) -> None:
    unchanged = h5ad_info["mtime_ns_before"] == mtime_after and h5ad_info["size_before"] == size_after
    lines = [
        "# gdTAI NK-Optimized gdT Atlas",
        "",
        f"- Input H5AD: `{h5ad_info['input_h5ad']}`",
        f"- Cells in atlas: `{h5ad_info['n_obs']:,}`",
        f"- Prediction mode: `gdTAI v2 high_purity / NK-optimized`",
        f"- Predicted gdT cells: `{int(overall['predicted_gdT'].iloc[0]):,}`",
        f"- Source H5AD unchanged: `{unchanged}`",
        "",
        "## Metadata Policy",
        "",
        "- Blood-like labels (`PBMC`, `peripheral blood`, `Blood`, `blood`, `BLD`) are collapsed to `blood`.",
        "- `cord blood` is kept separate as `cord_blood`.",
        "- Disease is reported as coarse `healthy`, `disease`, or `unknown`, plus detailed condition when available.",
        "- Age is reported only where joined metadata provide age or explicit fetal/newborn wording; missing age is not imputed.",
        "",
        "## Overall",
        dataframe_to_markdown(overall),
        "",
        "## Outputs",
        f"- Tables: `{TABLE_DIR}`",
        f"- Figures: `{FIGURE_DIR}`",
        f"- HTML: `{REPORT_HTML}`",
        f"- PDF: `{REPORT_PDF}`",
    ]
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    figure_html = "\n".join(
        f"<figure><img src='../{rel(fig)}'><figcaption>{html.escape(fig.name)}</figcaption></figure>" for fig in figures
    )
    css = """
    @page{size:A4 landscape;margin:8mm}
    body{font-family:Arial,Helvetica,sans-serif;margin:24px auto;max-width:1320px;color:#1f2933;line-height:1.45}
    h1{font-size:30px;margin:0 0 8px} h2{font-size:21px;margin-top:24px}
    .grid{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:16px 0}
    .metric{border:1px solid #d8dee5;background:#f8fafc;padding:12px}.value{font-size:24px;font-weight:bold}
    table{border-collapse:collapse;width:100%;font-size:10px;table-layout:fixed}
    th,td{border:1px solid #d8dee5;padding:4px 5px;text-align:left;vertical-align:top;overflow-wrap:anywhere}
    th{background:#eef2f6} img{max-width:100%;height:auto;border:1px solid #d8dee5}
    figure{break-inside:avoid;margin:18px 0}.note{background:#fff7ed;border-left:4px solid #f97316;padding:10px}
    """
    html_text = f"""<!doctype html><html><head><meta charset='utf-8'>
<title>gdTAI NK-Optimized gdT Atlas</title><style>{css}</style></head><body>
<h1>gdTAI NK-Optimized gdT Atlas</h1>
<p>This report uses gdTAI v2 <code>high_purity</code> predictions from the 5 million cell T/NK atlas and harmonizes metadata without modifying the source H5AD.</p>
<div class='grid'>
  <div class='metric'><div>Total atlas cells</div><div class='value'>{int(overall['total_cells'].iloc[0]):,}</div></div>
  <div class='metric'><div>Predicted gdT cells</div><div class='value'>{int(overall['predicted_gdT'].iloc[0]):,}</div></div>
  <div class='metric'><div>Predicted fraction</div><div class='value'>{float(overall['predicted_gdT_fraction'].iloc[0]):.2%}</div></div>
  <div class='metric'><div>H5AD unchanged</div><div class='value'>{unchanged}</div></div>
</div>
<section class='note'><strong>Metadata policy:</strong> blood-like labels are collapsed to <code>blood</code>, <code>cord blood</code> remains separate, disease is shown at coarse and detailed levels, and age is shown only where available.</section>
<h2>Figures</h2>{figure_html}
<h2>Overall</h2>{dataframe_to_html(overall)}
<h2>Top Sources</h2>{dataframe_to_html(source, 35)}
<h2>Tissue Site</h2>{dataframe_to_html(tissue, 35)}
<h2>Disease Status</h2>{dataframe_to_html(disease, 20)}
<h2>Age Group</h2>{dataframe_to_html(age, 20)}
<h2>Metadata Completeness By Source</h2>{dataframe_to_html(completeness, 35)}
<h2>Unresolved Metadata Examples</h2>{dataframe_to_html(unresolved, 50)}
</body></html>"""
    REPORT_HTML.write_text(html_text, encoding="utf-8")


def render_pdf(skip: bool) -> bool:
    if skip:
        return False
    chrome = "google-chrome"
    try:
        subprocess.run(
            [
                chrome,
                "--headless",
                "--disable-gpu",
                "--no-sandbox",
                f"--print-to-pdf={REPORT_PDF}",
                str(REPORT_HTML.resolve()),
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        return REPORT_PDF.exists()
    except Exception as exc:
        logging.warning("PDF render failed: %s", exc)
        return False


def main() -> None:
    args = parse_args()
    setup_logging()
    rules = json.loads(args.rules.read_text(encoding="utf-8"))

    logging.info("Loading H5AD obs metadata")
    df, h5ad_info = load_h5ad_obs(args.input_h5ad)
    pred = load_predictions(args.prediction_npz, int(h5ad_info["n_obs"]))
    df["predicted_gdT_high_purity"] = pred
    df["is_NK_cell"] = df["simple_annotation_plus6"].astype(str).str.upper().str.contains("NK", regex=False)
    df["predicted_paired_TCRAB_FP"] = df["predicted_gdT_high_purity"] & df["has_TRA_TRB_paired"].astype(bool)
    df["predicted_NK_FP"] = df["predicted_gdT_high_purity"] & df["is_NK_cell"].astype(bool)

    logging.info("Loading and applying upstream metadata lookups")
    lookups = load_metadata_lookup(args.metadata)
    df = apply_metadata_lookup(df, lookups)
    df = add_harmonized_metadata(df, rules)

    logging.info("Writing metadata audit tables")
    audit_cols = [
        "obs_name",
        "source_gse_id",
        "sample_id",
        "library_id",
        "tissue",
        "tissue_corrected",
        "tissue_meta",
        "sample_type_meta",
        "condition",
        "condition_meta",
        "age",
        "sex",
        "simple_annotation_plus6",
        "tissue_site_gdt_atlas",
        "specimen_context_gdt_atlas",
        "disease_status_gdt_atlas",
        "condition_detail_gdt_atlas",
        "age_group_gdt_atlas",
        "predicted_gdT_high_purity",
    ]
    sample_audit(df[audit_cols], args.sample_per_source, args.seed).to_csv(
        TABLE_DIR / "metadata_audit_sample_by_source.csv.gz", index=False
    )
    pred_df = df.loc[df["predicted_gdT_high_purity"], audit_cols + ["has_TRA_TRB_paired", "is_NK_cell"]].copy()
    pred_df.to_csv(TABLE_DIR / "predicted_gdt_high_purity_metadata.csv.gz", index=False)

    completeness = metadata_completeness(df)
    completeness.to_csv(TABLE_DIR / "metadata_completeness_by_source.csv", index=False)
    unresolved = unresolved_examples(df)
    unresolved.to_csv(TABLE_DIR / "unresolved_metadata_examples.csv", index=False)
    raw_tissue_counts = (
        df.groupby(["tissue", "tissue_corrected", "tissue_site_gdt_atlas"], dropna=False)
        .agg(total_cells=("predicted_gdT_high_purity", "size"), predicted_gdT=("predicted_gdT_high_purity", "sum"))
        .reset_index()
        .sort_values(["predicted_gdT", "total_cells"], ascending=False)
    )
    raw_tissue_counts.to_csv(TABLE_DIR / "raw_to_harmonized_tissue_counts.csv", index=False)

    logging.info("Writing atlas summaries")
    overall = pd.DataFrame(
        [
            {
                "strategy": "gdTAI_v2_high_purity_NK_optimized",
                "total_cells": int(df.shape[0]),
                "predicted_gdT": int(df["predicted_gdT_high_purity"].sum()),
                "predicted_gdT_fraction": float(df["predicted_gdT_high_purity"].mean()),
                "paired_TCRAB_cells": int(df["has_TRA_TRB_paired"].sum()),
                "predicted_paired_TCRAB_FP": int(df["predicted_paired_TCRAB_FP"].sum()),
                "NK_cells": int(df["is_NK_cell"].sum()),
                "predicted_NK_FP": int(df["predicted_NK_FP"].sum()),
                "predicted_paired_TCRAB_or_NK_FP": int(
                    (df["predicted_paired_TCRAB_FP"] | df["predicted_NK_FP"]).sum()
                ),
            }
        ]
    )
    overall["known_FP_fraction_of_predictions"] = overall["predicted_paired_TCRAB_or_NK_FP"] / overall[
        "predicted_gdT"
    ].replace(0, np.nan)
    overall.to_csv(TABLE_DIR / "atlas_overall_summary.csv", index=False)

    source = summarize_by(df, "predicted_gdT_high_purity", ["source_gse_id"], TABLE_DIR / "atlas_by_source.csv")
    tissue = summarize_by(
        df, "predicted_gdT_high_purity", ["tissue_site_gdt_atlas"], TABLE_DIR / "atlas_by_tissue_site.csv"
    )
    specimen = summarize_by(
        df,
        "predicted_gdT_high_purity",
        ["specimen_context_gdt_atlas"],
        TABLE_DIR / "atlas_by_specimen_context.csv",
    )
    disease = summarize_by(
        df, "predicted_gdT_high_purity", ["disease_status_gdt_atlas"], TABLE_DIR / "atlas_by_disease_status.csv"
    )
    condition = summarize_by(
        df, "predicted_gdT_high_purity", ["condition_detail_gdt_atlas"], TABLE_DIR / "atlas_by_condition_detail.csv"
    )
    age = summarize_by(df, "predicted_gdT_high_purity", ["age_group_gdt_atlas"], TABLE_DIR / "atlas_by_age_group.csv")
    sex = summarize_by(df, "predicted_gdT_high_purity", ["sex_gdt_atlas"], TABLE_DIR / "atlas_by_sex.csv")
    _ = specimen, condition, sex

    cross_tissue_disease = summarize_by(
        df,
        "predicted_gdT_high_purity",
        ["tissue_site_gdt_atlas", "disease_status_gdt_atlas"],
        TABLE_DIR / "atlas_tissue_by_disease.csv",
    )
    summarize_by(
        df,
        "predicted_gdT_high_purity",
        ["tissue_site_gdt_atlas", "age_group_gdt_atlas"],
        TABLE_DIR / "atlas_tissue_by_age_group.csv",
    )
    summarize_by(
        df,
        "predicted_gdT_high_purity",
        ["source_gse_id", "tissue_site_gdt_atlas", "disease_status_gdt_atlas"],
        TABLE_DIR / "atlas_source_by_tissue_disease.csv",
    )

    logging.info("Plotting figures")
    figures: list[Path] = []
    umap_path = plot_umap(
        args.input_h5ad,
        pred,
        FIGURE_DIR / "gdtai_high_purity_predicted_gdt_umap.png",
        args.umap_background,
        args.umap_positive,
        args.seed,
    )
    if umap_path is not None:
        figures.append(umap_path)
    figures.append(plot_bar(source, "source_gse_id", "predicted_gdT", FIGURE_DIR / "predicted_gdt_by_source.png", "Predicted gdT cells by source", "predicted gdT cells"))
    figures.append(plot_bar(tissue, "tissue_site_gdt_atlas", "predicted_gdT", FIGURE_DIR / "predicted_gdt_by_tissue_site.png", "Predicted gdT cells by harmonized tissue site", "predicted gdT cells"))
    figures.append(plot_bar(disease, "disease_status_gdt_atlas", "predicted_gdT", FIGURE_DIR / "predicted_gdt_by_disease_status.png", "Predicted gdT cells by disease status", "predicted gdT cells"))
    figures.append(plot_bar(age, "age_group_gdt_atlas", "predicted_gdT", FIGURE_DIR / "predicted_gdt_by_age_group.png", "Predicted gdT cells by age group where available", "predicted gdT cells"))
    figures.append(plot_heatmap(cross_tissue_disease, FIGURE_DIR / "predicted_gdt_tissue_by_disease_heatmap.png", "Predicted gdT cells: tissue by disease"))

    stat_after = args.input_h5ad.stat()
    write_report(
        overall=overall,
        source=source,
        tissue=tissue,
        disease=disease,
        age=age,
        completeness=completeness,
        unresolved=unresolved,
        figures=figures,
        h5ad_info=h5ad_info,
        mtime_after=stat_after.st_mtime_ns,
        size_after=stat_after.st_size,
    )
    pdf_ok = render_pdf(args.no_pdf)
    manifest = {
        **h5ad_info,
        "prediction_npz": str(args.prediction_npz),
        "prediction_key": PREDICTION_KEY,
        "predicted_gdT_high_purity": int(pred.sum()),
        "size_after": stat_after.st_size,
        "mtime_ns_after": stat_after.st_mtime_ns,
        "h5ad_unchanged": bool(h5ad_info["size_before"] == stat_after.st_size and h5ad_info["mtime_ns_before"] == stat_after.st_mtime_ns),
        "html": str(REPORT_HTML),
        "pdf": str(REPORT_PDF),
        "pdf_rendered": bool(pdf_ok),
    }
    (LOG_DIR / "gdtai_nk_optimized_gdt_atlas_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    logging.info("Wrote report: %s", REPORT_HTML)
    logging.info("H5AD unchanged: %s", manifest["h5ad_unchanged"])


if __name__ == "__main__":
    main()
