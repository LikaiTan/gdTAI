#!/usr/bin/env python3
"""Run the precommitted gdTAI V4 Step 1 preflight without fitting a model."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import logging
import math
import os
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.model_selection import StratifiedGroupKFold


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "configs/models/gdtai/v4_preflight.json"
TABLE_DIR = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_preflight"
FIGURE_DIR = ROOT / "Integrated_dataset/figures/gdT_prediction/gdtai_v4_preflight"
LOG_DIR = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_preflight"
STATIC_DIR = ROOT / "gdT_prediction/gdtai_v4_preflight"
STATIC_ASSET_DIR = STATIC_DIR / "assets"
ARCHIVE_DIR = ROOT / "archive/retired_experiments/gdTAI_v4_experimental_precommit_20260807"

INVALID_STRINGS = {"", "na", "nan", "none", "null", "missing", "unknown"}
PRIMARY_GDT = "gdT_primary"
PRIMARY_ABT = "abT_primary"
SINGLE_ABT = "single_abT_weak"
SORTED_WEAK = "sorted_gdT_weak"
SORTED_SENSITIVITY = "sorted_sensitivity"
SILVER = "gdT_silver"
DUAL = "dual_or_ambiguous"


@dataclass(frozen=True)
class InputSpec:
    input_id: str
    path: Path
    matrix_key: str | None
    role: str
    source_id: str
    expected_sha256: str = ""
    hard_feature_gate: bool = True
    matrix_state: str = "raw_counts"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--skip-full-hashes", action="store_true")
    parser.add_argument("--log-level", default="INFO")
    return parser.parse_args()


def ensure_dirs() -> None:
    for path in (TABLE_DIR, FIGURE_DIR, LOG_DIR, STATIC_DIR, STATIC_ASSET_DIR):
        path.mkdir(parents=True, exist_ok=True)


def setup_logging(level: str) -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    handlers: list[logging.Handler] = [
        logging.StreamHandler(),
        logging.FileHandler(LOG_DIR / "gdtai_v4_preflight.log", mode="w"),
    ]
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=handlers,
        force=True,
    )


def resolve_path(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def decode(value: Any) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def read_string_dataset(dataset: h5py.Dataset) -> np.ndarray:
    try:
        return np.asarray(dataset.asstr()[:], dtype=object)
    except Exception:
        return np.asarray([decode(value) for value in dataset[:]], dtype=object)


def read_axis_column(group: h5py.Group, column: str) -> np.ndarray:
    obj = group[column]
    if isinstance(obj, h5py.Group):
        encoding = decode(obj.attrs.get("encoding-type", ""))
        if encoding == "categorical":
            categories = read_string_dataset(obj["categories"])
            codes = np.asarray(obj["codes"][:], dtype=np.int64)
            out = np.full(codes.shape, "", dtype=object)
            valid = codes >= 0
            out[valid] = categories[codes[valid]]
            return out
        if "values" in obj and "mask" in obj:
            values = obj["values"][:]
            mask = np.asarray(obj["mask"][:], dtype=bool)
            if values.dtype.kind in {"S", "O", "U"}:
                out = np.asarray([decode(value) for value in values], dtype=object)
                out[mask] = ""
                return out
            out = np.asarray(values)
            if mask.any():
                out = out.astype(object)
                out[mask] = np.nan
            return out
        raise TypeError(f"Unsupported encoded column {column}: {encoding}")
    values = obj[:]
    if values.dtype.kind in {"S", "O", "U"}:
        return read_string_dataset(obj)
    return np.asarray(values)


def axis_index_key(group: h5py.Group) -> str:
    key = group.attrs.get("_index")
    if key is not None:
        return decode(key)
    if "_index" in group:
        return "_index"
    raise KeyError("H5AD axis has no index key")


def read_obs(handle: h5py.File, column: str) -> np.ndarray:
    return read_axis_column(handle["obs"], column)


def read_obs_index(handle: h5py.File) -> np.ndarray:
    return read_axis_column(handle["obs"], axis_index_key(handle["obs"]))


def read_var_names(handle: h5py.File) -> list[str]:
    return [str(value) for value in read_axis_column(handle["var"], axis_index_key(handle["var"]))]


def clean_strings(values: np.ndarray, default: str = "") -> np.ndarray:
    series = pd.Series(values, copy=False).astype("string").fillna("").str.strip()
    invalid = series.str.lower().isin(INVALID_STRINGS)
    return series.mask(invalid, default).astype(str).to_numpy(dtype=object)


def read_bool(handle: h5py.File, column: str, default: bool = False) -> np.ndarray:
    if column not in handle["obs"]:
        return np.full(len(read_obs_index(handle)), default, dtype=bool)
    values = read_obs(handle, column)
    if np.issubdtype(values.dtype, np.bool_):
        return values.astype(bool, copy=False)
    if np.issubdtype(values.dtype, np.number):
        numeric = np.asarray(values, dtype=float)
        return np.isfinite(numeric) & (numeric != 0)
    lowered = pd.Series(values, copy=False).astype("string").fillna("").str.strip().str.lower()
    return lowered.isin({"true", "1", "yes", "y", "t"}).to_numpy(dtype=bool)


def nonempty(values: np.ndarray) -> np.ndarray:
    cleaned = clean_strings(values)
    return cleaned != ""


def get_matrix_group(handle: h5py.File, matrix_key: str) -> h5py.Group:
    obj: h5py.Group | h5py.Dataset = handle
    for token in matrix_key.strip("/").split("/"):
        obj = obj[token]  # type: ignore[index]
    if not isinstance(obj, h5py.Group):
        raise TypeError(f"Matrix {matrix_key} is dense; sparse CSR is required")
    return obj


def file_sha256(path: Path, chunk_bytes: int = 64 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    total = path.stat().st_size
    processed = 0
    last_log = time.monotonic()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_bytes)
            if not chunk:
                break
            digest.update(chunk)
            processed += len(chunk)
            if time.monotonic() - last_log >= 20:
                logging.info("Hashing %s: %.1f%%", path.name, 100.0 * processed / max(total, 1))
                last_log = time.monotonic()
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def canonical_json_sha256(value: Any) -> str:
    return sha256_text(json.dumps(value, sort_keys=True, separators=(",", ":")))


def parse_protocol_tcr_genes(protocol_path: Path) -> list[str]:
    text = protocol_path.read_text(encoding="utf-8")
    appendix = text.split("## Appendix A.", 1)[1].split("## Appendix B.", 1)[0]
    return sorted(set(re.findall(r"\b(?:TRA|TRB|TRG|TRD)[A-Z0-9-]+\b", appendix)))


def deterministic_gzip_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=6, mtime=0) as gz:
            with io.TextIOWrapper(gz, encoding="utf-8", newline="") as text:
                df.to_csv(text, index=False, lineterminator="\n")
    os.replace(temp, path)


def write_json(value: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temp, path)


def valid_group_value(value: str) -> bool:
    return value.strip().lower() not in INVALID_STRINGS


def clone_signatures(arrays: list[np.ndarray], rows: np.ndarray) -> np.ndarray:
    out = np.full(rows.size, "", dtype=object)
    for pos, row in enumerate(rows):
        chains = sorted({clean_strings(np.asarray([array[row]]))[0] for array in arrays} - {""})
        if chains:
            out[pos] = hashlib.sha1("|".join(chains).encode("utf-8")).hexdigest()[:20]
    return out


def build_group_keys(
    source: np.ndarray,
    donor: np.ndarray,
    sample: np.ndarray,
    library: np.ndarray,
    clonotype: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    groups = np.empty(source.size, dtype=object)
    levels = np.empty(source.size, dtype=object)
    for idx in range(source.size):
        src = str(source[idx])
        if valid_group_value(str(donor[idx])):
            level, value = "donor", str(donor[idx])
        elif valid_group_value(str(clonotype[idx])):
            level, value = "clonotype", str(clonotype[idx])
        elif valid_group_value(str(sample[idx])):
            level, value = "sample", str(sample[idx])
        elif valid_group_value(str(library[idx])):
            level, value = "library", str(library[idx])
        else:
            raise ValueError(f"No valid group key for source {src} at local row {idx}")
        groups[idx] = f"{src}::{level}::{value}"
        levels[idx] = level
    return groups, levels


def make_labeled_frame(
    *,
    input_id: str,
    source_id: str,
    cell_id: np.ndarray,
    donor: np.ndarray,
    sample: np.ndarray,
    library: np.ndarray,
    any_ab: np.ndarray,
    any_gd: np.ndarray,
    pair_ab: np.ndarray,
    pair_gd: np.ndarray,
    trd_evidence: np.ndarray,
    doublet: np.ndarray,
    nk_mask: np.ndarray,
    nk_strength: np.ndarray,
    clonotype_arrays: list[np.ndarray],
    expression_rows: np.ndarray | None = None,
) -> pd.DataFrame:
    n_cells = cell_id.size
    rows = np.arange(n_cells, dtype=np.int64) if expression_rows is None else expression_rows.astype(np.int64)
    primary_gdt = pair_gd & (~any_ab) & (~doublet)
    primary_abt = pair_ab & (~any_gd) & (~doublet)
    dual = any_ab & any_gd
    single_abt = any_ab & (~pair_ab) & (~any_gd) & (~doublet)
    silver = trd_evidence & (~pair_gd) & (~any_ab) & (~doublet)
    t_positive = (any_ab | any_gd) & (~silver) & (~dual) & (~doublet)
    nk_negative = nk_mask & (~any_ab) & (~any_gd) & (~doublet)

    truth = np.full(n_cells, "unlabeled", dtype=object)
    truth[dual] = DUAL
    truth[single_abt] = SINGLE_ABT
    truth[silver] = SILVER
    truth[primary_abt] = PRIMARY_ABT
    truth[primary_gdt] = PRIMARY_GDT

    stage2 = np.full(n_cells, "none", dtype=object)
    stage2[primary_gdt] = "positive"
    stage2[primary_abt] = "negative"
    stage2[single_abt] = "negative"
    reliability = np.zeros(n_cells, dtype=np.float32)
    reliability[primary_gdt | primary_abt] = 1.0
    reliability[single_abt] = 0.5

    stage1 = np.full(n_cells, "none", dtype=object)
    stage1[t_positive] = "t_positive"
    stage1[nk_negative] = "nk_negative"
    stage1_weight = np.zeros(n_cells, dtype=np.float32)
    stage1_weight[t_positive] = 1.0
    stage1_weight[nk_negative] = nk_strength[nk_negative]

    include = (stage1 != "none") | (stage2 != "none") | (truth == SILVER) | (truth == DUAL)
    selected = np.flatnonzero(include)
    missing_donor = np.asarray([not valid_group_value(str(value)) for value in donor[selected]], dtype=bool)
    clone = np.full(selected.size, "", dtype=object)
    if missing_donor.any() and clonotype_arrays:
        clone[missing_donor] = clone_signatures(clonotype_arrays, selected[missing_donor])
    groups, levels = build_group_keys(
        np.full(selected.size, source_id, dtype=object),
        clean_strings(donor[selected]),
        clean_strings(sample[selected]),
        clean_strings(library[selected]),
        clone,
    )
    return pd.DataFrame(
        {
            "cell_id": clean_strings(cell_id[selected]),
            "source_gse_id": source_id,
            "expression_input_id": input_id,
            "expression_row": rows[selected],
            "truth_class": truth[selected],
            "truth_reliability": reliability[selected],
            "stage1_role": stage1[selected],
            "stage1_weight": stage1_weight[selected],
            "stage2_role": stage2[selected],
            "group_key": groups,
            "group_level": levels,
            "clonotype_key": clone,
            "has_any_ab_tcr": any_ab[selected],
            "has_any_gd_tcr": any_gd[selected],
            "has_paired_ab_tcr": pair_ab[selected],
            "has_paired_gd_tcr": pair_gd[selected],
            "doublet_flag": doublet[selected],
            "nk_annotation_strength": nk_strength[selected],
            "nk_sampling_stratum": "",
            "cd4_helper_exclusion": False,
            "treg_exclusion": False,
            "exclusion_union": False,
        }
    )


def load_hra_frame(path: Path) -> pd.DataFrame:
    logging.info("Loading HRA005041 labels")
    with h5py.File(path, "r") as handle:
        ids = clean_strings(read_obs(handle, "cell_id"))
        donor = clean_strings(read_obs(handle, "donor_id"))
        sample = clean_strings(read_obs(handle, "sample_id"))
        library = clean_strings(read_obs(handle, "library_id"))
        any_ab = read_bool(handle, "has_any_ab_tcr")
        any_gd = read_bool(handle, "has_any_gd_tcr")
        pair_ab = read_bool(handle, "has_TRA_TRB_paired")
        pair_gd = read_bool(handle, "has_TRG_TRD_paired")
        trd = nonempty(read_obs(handle, "TRD_cdr3"))
        major = pd.Series(read_obs(handle, "Major_cluster"), copy=False).astype("string").fillna("")
        detail = pd.Series(read_obs(handle, "original_cell_annotation"), copy=False).astype("string").fillna("")
        nk_primary = major.str.fullmatch("NK", case=False).to_numpy(dtype=bool)
        nk_secondary = detail.str.contains(r"(?:^|[^A-Za-z])NK(?:[^A-Za-z]|$)", case=False, regex=True).to_numpy(dtype=bool)
        nk_strength = np.where(nk_primary & nk_secondary, 1.0, 0.5).astype(np.float32)
        nk_mask = nk_primary
        doublet = np.zeros(ids.size, dtype=bool)
        clonotypes = [read_obs(handle, column) for column in ("TRA_cdr3", "TRB_cdr3", "TRG_cdr3", "TRD_cdr3")]
    return make_labeled_frame(
        input_id="hra005041",
        source_id="HRA005041",
        cell_id=ids,
        donor=donor,
        sample=sample,
        library=library,
        any_ab=any_ab,
        any_gd=any_gd,
        pair_ab=pair_ab,
        pair_gd=pair_gd,
        trd_evidence=trd,
        doublet=doublet,
        nk_mask=nk_mask,
        nk_strength=nk_strength,
        clonotype_arrays=clonotypes,
    )


def load_legacy_reference(path: Path) -> dict[str, np.ndarray]:
    logging.info("Loading legacy metadata reference needed for TCR and GSE144469 joins")
    with h5py.File(path, "r") as handle:
        return {
            "cell_id": clean_strings(read_obs_index(handle)),
            "source": clean_strings(read_obs(handle, "source_gse_id")),
            "sample": clean_strings(read_obs(handle, "sample_id")),
            "library": clean_strings(read_obs(handle, "library_id")),
            "barcode": clean_strings(read_obs(handle, "barcodes")),
            "donor": clean_strings(read_obs(handle, "donor_id")),
            "any_ab": read_bool(handle, "has_any_ab_tcr"),
            "any_gd": read_bool(handle, "has_any_gd_tcr"),
        }


def load_gse144469_frame(path: Path, legacy: dict[str, np.ndarray]) -> tuple[pd.DataFrame, pd.DataFrame]:
    logging.info("Loading and joining GSE144469")
    legacy_mask = legacy["source"] == "GSE144469"
    legacy_rows = np.flatnonzero(legacy_mask)

    def id_srr_barcode(value: str) -> tuple[str, str]:
        local_id = str(value).removeprefix("GSE144469__")
        parts = local_id.split("_", 1)
        if len(parts) != 2 or not re.fullmatch(r"SRR\d+", parts[0]):
            raise ValueError(f"Invalid GSE144469 cell identifier: {value}")
        return parts[0], parts[1]

    canonical_parts = [id_srr_barcode(value) for value in legacy["cell_id"][legacy_rows]]
    canonical_keys = pd.Index([f"{srr}||{barcode}" for srr, barcode in canonical_parts])
    with h5py.File(path, "r") as handle:
        raw_ids = clean_strings(read_obs_index(handle))
        sample = clean_strings(read_obs(handle, "sample"))
        source_parts = [id_srr_barcode(value) for value in raw_ids]
        id_srr = np.asarray([value[0] for value in source_parts], dtype=object)
        raw_barcode = np.asarray([value[1] for value in source_parts], dtype=object)
        source_keys = pd.Index([f"{srr}||{barcode}" for srr, barcode in source_parts])
        mapped_local = canonical_keys.get_indexer(source_keys)
        mapped_ok = mapped_local >= 0
        if not mapped_ok.all():
            raise RuntimeError(f"GSE144469 join mapped {int(mapped_ok.sum())}/{mapped_ok.size} cells")
        mapped = legacy_rows[mapped_local]
        has_tra = read_bool(handle, "has_TRA")
        has_trb = read_bool(handle, "has_TRB")
        has_trg = read_bool(handle, "has_TRG")
        has_trd = read_bool(handle, "has_TRD")
        any_ab = has_tra | has_trb
        any_gd = has_trg | has_trd
        pair_ab = has_tra & has_trb
        pair_gd = has_trg & has_trd
        cell_type = pd.Series(read_obs(handle, "cell_type"), copy=False).astype("string").fillna("")
        nk_mask = cell_type.str.fullmatch("NK", case=False).to_numpy(dtype=bool)
        nk_strength = np.full(raw_ids.size, 0.5, dtype=np.float32)
        doublet = np.zeros(raw_ids.size, dtype=bool)

    join_audit = pd.DataFrame(
        [
            {
                "source_gse_id": "GSE144469",
                "source_cells": int(raw_ids.size),
                "source_unique_keys": int(source_keys.nunique()),
                "canonical_cells": int(legacy_rows.size),
                "canonical_unique_keys": int(canonical_keys.nunique()),
                "mapped_cells": int(mapped_ok.sum()),
                "unmapped_cells": int((~mapped_ok).sum()),
                "mapping_fraction": float(mapped_ok.mean()),
                "source_duplicate_keys": int(source_keys.duplicated(keep=False).sum()),
                "canonical_duplicate_keys": int(canonical_keys.duplicated(keep=False).sum()),
                "obs_sample_disagrees_with_id_srr": int((sample != id_srr).sum()),
                "obs_sample_barcode_unique_keys": int(
                    pd.Index(pd.Series(sample).astype(str) + "||" + pd.Series(raw_barcode).astype(str)).nunique()
                ),
            }
        ]
    )
    frame = make_labeled_frame(
        input_id="gse144469",
        source_id="GSE144469",
        cell_id=legacy["cell_id"][mapped],
        donor=legacy["donor"][mapped],
        sample=legacy["sample"][mapped],
        library=legacy["library"][mapped],
        any_ab=any_ab,
        any_gd=any_gd,
        pair_ab=pair_ab,
        pair_gd=pair_gd,
        trd_evidence=has_trd,
        doublet=doublet,
        nk_mask=nk_mask,
        nk_strength=nk_strength,
        clonotype_arrays=[],
        expression_rows=np.arange(raw_ids.size, dtype=np.int64),
    )
    return frame, join_audit


def load_balf_frame(path: Path) -> pd.DataFrame:
    logging.info("Loading BALF_BLOOD_COPD labels")
    with h5py.File(path, "r") as handle:
        ids = clean_strings(read_obs_index(handle))
        donor = clean_strings(read_obs(handle, "donor_id"))
        sample = clean_strings(read_obs(handle, "sample_id"))
        library = clean_strings(read_obs(handle, "library_id"))
        tra = read_bool(handle, "productive_TRA")
        trb = read_bool(handle, "productive_TRB")
        trg = read_bool(handle, "productive_TRG")
        trd = read_bool(handle, "productive_TRD")
        any_ab, any_gd = tra | trb, trg | trd
        pair_ab = read_bool(handle, "paired_TRA_TRB_productive")
        pair_gd = read_bool(handle, "paired_TRG_TRD_productive")
        level1 = pd.Series(read_obs(handle, "cell_type_level1"), copy=False).astype("string").fillna("")
        in_scope = level1.isin(["T_cell", "NK"]).to_numpy(dtype=bool)
        nk_mask = (level1 == "NK").to_numpy(dtype=bool)
        singlet = read_bool(handle, "demux_biological_singlet", default=True)
        doublet = ~singlet
        any_ab &= in_scope
        any_gd &= in_scope
        pair_ab &= in_scope
        pair_gd &= in_scope
        nk_mask &= in_scope
        nk_strength = np.full(ids.size, 0.5, dtype=np.float32)
        clonotypes = [read_obs(handle, column) for column in ("TRA_cdr3_aa", "TRB_cdr3_aa", "TRG_cdr3_aa", "TRD_cdr3_aa")]
    return make_labeled_frame(
        input_id="balf_blood_copd",
        source_id="BALF_BLOOD_COPD",
        cell_id=ids,
        donor=donor,
        sample=sample,
        library=library,
        any_ab=any_ab,
        any_gd=any_gd,
        pair_ab=pair_ab,
        pair_gd=pair_gd,
        trd_evidence=trd & in_scope,
        doublet=doublet,
        nk_mask=nk_mask,
        nk_strength=nk_strength,
        clonotype_arrays=clonotypes,
    )


def load_sorted_frame(path: Path, input_id: str, source_id: str) -> pd.DataFrame:
    logging.info("Loading sorted source %s", source_id)
    with h5py.File(path, "r") as handle:
        source_ids = clean_strings(read_obs(handle, "cell_id") if "cell_id" in handle["obs"] else read_obs_index(handle))
        ids = np.asarray([f"{source_id}__{value}" for value in source_ids], dtype=object)
        donor = clean_strings(read_obs(handle, "donor_id")) if "donor_id" in handle["obs"] else np.full(ids.size, "", dtype=object)
        sample = clean_strings(read_obs(handle, "sample_id"))
        library = clean_strings(read_obs(handle, "library_id"))
        any_ab = read_bool(handle, "has_any_ab_tcr")
        any_gd = read_bool(handle, "has_any_gd_tcr")
        pair_ab = read_bool(handle, "has_TRA_TRB_paired")
        pair_gd = read_bool(handle, "has_TRG_TRD_paired")
        clonotype_arrays = [read_obs(handle, column) for column in ("TRA_cdr3", "TRB_cdr3", "TRG_cdr3", "TRD_cdr3") if column in handle["obs"]]
    selected = np.arange(ids.size, dtype=np.int64)
    missing_donor = np.asarray([not valid_group_value(str(value)) for value in donor], dtype=bool)
    clone = np.full(ids.size, "", dtype=object)
    if missing_donor.any() and clonotype_arrays:
        clone[missing_donor] = clone_signatures(clonotype_arrays, selected[missing_donor])
    groups, levels = build_group_keys(
        np.full(ids.size, source_id, dtype=object), donor, sample, library, clone
    )
    return pd.DataFrame(
        {
            "cell_id": ids,
            "source_gse_id": source_id,
            "expression_input_id": input_id,
            "expression_row": selected,
            "truth_class": SORTED_SENSITIVITY,
            "truth_reliability": 0.0,
            "stage1_role": "none",
            "stage1_weight": 0.0,
            "stage2_role": "none",
            "group_key": groups,
            "group_level": levels,
            "clonotype_key": clone,
            "has_any_ab_tcr": any_ab,
            "has_any_gd_tcr": any_gd,
            "has_paired_ab_tcr": pair_ab,
            "has_paired_gd_tcr": pair_gd,
            "doublet_flag": False,
            "nk_annotation_strength": 0.0,
            "nk_sampling_stratum": "",
            "cd4_helper_exclusion": False,
            "treg_exclusion": False,
            "exclusion_union": False,
        }
    )


def load_atlas_supplemental_nk(path: Path, legacy: dict[str, np.ndarray], primary_sources: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    logging.info("Building strict supplemental NK pool from the current atlas")
    with h5py.File(path, "r") as handle:
        ids = clean_strings(read_obs_index(handle))
        source = clean_strings(read_obs(handle, "source_gse_id"))
        donor = clean_strings(read_obs(handle, "donor_id"))
        sample = clean_strings(read_obs(handle, "sample_id"))
        library = clean_strings(read_obs(handle, "library_id"))
        scanvi = clean_strings(read_obs(handle, "scanvi_tnk_superclass"))
        author = clean_strings(read_obs(handle, "phase1_annotation_label"))

    candidate = (scanvi == "NK_cell") & (~np.isin(source, list(primary_sources)))
    candidate_rows = np.flatnonzero(candidate)
    author_series = pd.Series(author[candidate_rows]).astype("string").fillna("").str.upper()
    author_nk = author_series.str.contains(r"(?:^|[^A-Z])NK(?:[^A-Z]|$)", regex=True).to_numpy(dtype=bool)
    author_t = author_series.str.contains(r"CD4|CD8|TREG|TCRGD|GDT|MAIT|T[ _-]?CELL", regex=True).to_numpy(dtype=bool)

    legacy_index = pd.Index(legacy["cell_id"])
    mapped = legacy_index.get_indexer(ids[candidate_rows])
    mapped_ok = mapped >= 0
    no_tcr = np.zeros(candidate_rows.size, dtype=bool)
    no_tcr[mapped_ok] = (~legacy["any_ab"][mapped[mapped_ok]]) & (~legacy["any_gd"][mapped[mapped_ok]])
    annotation_ok = ~(author_t & (~author_nk))
    groupable = np.asarray(
        [
            valid_group_value(str(donor[row]))
            or valid_group_value(str(sample[row]))
            or valid_group_value(str(library[row]))
            for row in candidate_rows
        ],
        dtype=bool,
    )
    keep = mapped_ok & no_tcr & annotation_ok & groupable
    rows = candidate_rows[keep]
    strength = np.where(author_nk[keep], 1.0, 0.5).astype(np.float32)
    clone = np.full(rows.size, "", dtype=object)
    groups, levels = build_group_keys(source[rows], donor[rows], sample[rows], library[rows], clone)
    frame = pd.DataFrame(
        {
            "cell_id": ids[rows],
            "source_gse_id": source[rows],
            "expression_input_id": "current_atlas",
            "expression_row": rows,
            "truth_class": "unlabeled",
            "truth_reliability": 0.0,
            "stage1_role": "nk_negative",
            "stage1_weight": strength,
            "stage2_role": "none",
            "group_key": groups,
            "group_level": levels,
            "clonotype_key": clone,
            "has_any_ab_tcr": False,
            "has_any_gd_tcr": False,
            "has_paired_ab_tcr": False,
            "has_paired_gd_tcr": False,
            "doublet_flag": False,
            "nk_annotation_strength": strength,
            "nk_sampling_stratum": "",
            "cd4_helper_exclusion": False,
            "treg_exclusion": False,
            "exclusion_union": False,
        }
    )
    audit = pd.DataFrame(
        [
            {"stage": "scanvi_NK_candidates_nonprimary", "n_cells": int(candidate_rows.size)},
            {"stage": "mapped_to_legacy_TCR_metadata", "n_cells": int(mapped_ok.sum())},
            {"stage": "mapped_no_productive_TCR", "n_cells": int((mapped_ok & no_tcr).sum())},
            {"stage": "author_annotation_not_T_conflict", "n_cells": int((mapped_ok & no_tcr & annotation_ok).sum())},
            {"stage": "valid_donor_sample_or_library_group", "n_cells": int((mapped_ok & no_tcr & annotation_ok & groupable).sum())},
            {"stage": "strict_supplemental_NK_final", "n_cells": int(frame.shape[0])},
            {"stage": "consensus_NK_final", "n_cells": int((strength == 1.0).sum())},
            {"stage": "single_annotation_NK_final", "n_cells": int((strength == 0.5).sum())},
        ]
    )
    return frame, audit


def scan_sparse_matrix(
    path: Path,
    matrix_key: str,
    selected_rows: np.ndarray,
    marker_genes: list[str],
    data_chunk: int,
    row_chunk: int,
    matrix_state: str,
    transformed_target: float,
    transformed_tolerance: float,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    logging.info("Scanning expression matrix %s [%s; %s]", path, matrix_key, matrix_state)
    selected_rows = np.asarray(sorted(set(map(int, selected_rows))), dtype=np.int64)
    with h5py.File(path, "r") as handle:
        matrix = get_matrix_group(handle, matrix_key)
        encoding = decode(matrix.attrs.get("encoding-type", ""))
        if encoding != "csr_matrix":
            raise TypeError(f"{path}:{matrix_key} uses {encoding}; CSR is required")
        shape = tuple(int(value) for value in matrix.attrs["shape"])
        indptr = np.asarray(matrix["indptr"][:], dtype=np.int64)
        data_ds = matrix["data"]
        indices_ds = matrix["indices"]
        data_dtype = str(data_ds.dtype)
        if indptr.size != shape[0] + 1 or int(indptr[-1]) != int(data_ds.shape[0]):
            raise ValueError(f"Invalid CSR structure in {path}:{matrix_key}")
        if selected_rows.size and (selected_rows.min() < 0 or selected_rows.max() >= shape[0]):
            raise IndexError(f"Selected expression row is outside {path}:{matrix_key}")
        var_names = read_var_names(handle)
        var_lookup = {gene: idx for idx, gene in enumerate(var_names)}
        present_markers = [gene for gene in marker_genes if gene in var_lookup]
        marker_columns = np.asarray([var_lookup[gene] for gene in present_markers], dtype=np.int64)
        output_columns = {gene: idx for idx, gene in enumerate(marker_genes)}
        marker_output = np.zeros((selected_rows.size, len(marker_genes)), dtype=np.float32)
        totals_output = np.zeros(selected_rows.size, dtype=np.float64)

        nonfinite = 0
        negative = 0
        fractional = 0
        max_fractional_distance = 0.0
        transformed_rows_checked = 0
        transformed_empty_rows = 0
        transformed_rows_outside_tolerance = 0
        transformed_sum_min = float("inf")
        transformed_sum_max = float("-inf")
        transformed_max_abs_deviation = 0.0
        nnz = int(data_ds.shape[0])
        processed_nnz = 0
        last_log = time.monotonic()

        for row_start in range(0, shape[0], row_chunk):
            row_end = min(shape[0], row_start + row_chunk)
            data_start, data_end = int(indptr[row_start]), int(indptr[row_end])
            data = np.asarray(data_ds[data_start:data_end])
            indices = np.asarray(indices_ds[data_start:data_end], dtype=np.int64)
            local_indptr = indptr[row_start : row_end + 1] - data_start
            finite = np.isfinite(data)
            nonfinite += int((~finite).sum())
            negative += int((data[finite] < 0).sum())
            if finite.any():
                distance = np.abs(data[finite] - np.rint(data[finite]))
                fractional += int((distance > 1e-6).sum())
                max_fractional_distance = max(max_fractional_distance, float(distance.max(initial=0.0)))

            if matrix_state == "log1p_cp10k_registered":
                row_nnz = np.diff(local_indptr)
                nonempty_rows = np.flatnonzero(row_nnz > 0)
                transformed_rows_checked += int(row_nnz.size)
                transformed_empty_rows += int((row_nnz == 0).sum())
                if nonempty_rows.size:
                    transformed_data = np.expm1(data.astype(np.float64, copy=False))
                    row_sums = np.add.reduceat(transformed_data, local_indptr[:-1][nonempty_rows])
                    deviations = np.abs(row_sums - transformed_target)
                    transformed_rows_outside_tolerance += int((deviations > transformed_tolerance).sum())
                    transformed_sum_min = min(transformed_sum_min, float(row_sums.min()))
                    transformed_sum_max = max(transformed_sum_max, float(row_sums.max()))
                    transformed_max_abs_deviation = max(
                        transformed_max_abs_deviation,
                        float(deviations.max(initial=0.0)),
                    )

            left = int(np.searchsorted(selected_rows, row_start, side="left"))
            right = int(np.searchsorted(selected_rows, row_end, side="left"))
            if right > left:
                local_rows = selected_rows[left:right] - row_start
                chunk_matrix = sparse.csr_matrix((data, indices, local_indptr), shape=(row_end - row_start, shape[1]))
                selected_matrix = chunk_matrix[local_rows]
                totals_output[left:right] = np.asarray(selected_matrix.sum(axis=1)).ravel()
                if present_markers:
                    values = selected_matrix[:, marker_columns].toarray().astype(np.float32, copy=False)
                    for local_col, gene in enumerate(present_markers):
                        marker_output[left:right, output_columns[gene]] = values[:, local_col]

            processed_nnz += data.size
            if time.monotonic() - last_log >= 20:
                logging.info("Matrix scan %s: %.1f%% of nonzeros", path.name, 100.0 * processed_nnz / max(nnz, 1))
                last_log = time.monotonic()

    raw_count_pass = bool(nonfinite == 0 and negative == 0 and fractional == 0)
    transformed_cp10k_pass = bool(
        matrix_state == "log1p_cp10k_registered"
        and nonfinite == 0
        and negative == 0
        and transformed_empty_rows == 0
        and transformed_rows_outside_tolerance == 0
        and transformed_rows_checked == shape[0]
    )
    expression_contract_pass = raw_count_pass if matrix_state == "raw_counts" else transformed_cp10k_pass
    audit = pd.DataFrame(
        [
            {
                "path": str(path),
                "matrix_key": matrix_key,
                "configured_matrix_state": matrix_state,
                "n_obs": shape[0],
                "n_vars": shape[1],
                "nnz": nnz,
                "dtype": data_dtype,
                "nonfinite_values": nonfinite,
                "negative_values": negative,
                "fractional_values_gt_1e_6": fractional,
                "max_fractional_distance": max_fractional_distance,
                "raw_count_pass": raw_count_pass,
                "transformed_rows_checked": transformed_rows_checked,
                "transformed_empty_rows": transformed_empty_rows,
                "transformed_target": transformed_target if transformed_rows_checked else np.nan,
                "transformed_tolerance": transformed_tolerance if transformed_rows_checked else np.nan,
                "transformed_sum_min": transformed_sum_min if transformed_rows_checked else np.nan,
                "transformed_sum_max": transformed_sum_max if transformed_rows_checked else np.nan,
                "transformed_max_abs_deviation": transformed_max_abs_deviation if transformed_rows_checked else np.nan,
                "transformed_rows_outside_tolerance": transformed_rows_outside_tolerance,
                "transformed_cp10k_pass": transformed_cp10k_pass,
                "expression_contract_pass": expression_contract_pass,
                "selected_rows_extracted": int(selected_rows.size),
            }
        ]
    )
    return audit, marker_output, totals_output


def apply_expression_audits(
    cells: pd.DataFrame,
    input_specs: dict[str, InputSpec],
    marker_genes: list[str],
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    marker_frame = pd.DataFrame(index=cells.index, columns=marker_genes, dtype=np.float32)
    total_counts = pd.Series(np.nan, index=cells.index, dtype=float)
    matrix_audits: list[pd.DataFrame] = []
    request_mask = (
        cells["stage1_role"].eq("nk_negative")
        | cells["truth_class"].isin([PRIMARY_GDT, SILVER, SORTED_WEAK, SORTED_SENSITIVITY])
    )
    transformed_contract = config["transformed_input_contract"]
    for input_id, spec in input_specs.items():
        if spec.matrix_key is None:
            continue
        subset = cells.index[request_mask & cells["expression_input_id"].eq(input_id)]
        row_values = cells.loc[subset, "expression_row"].to_numpy(dtype=np.int64)
        order = np.argsort(row_values, kind="mergesort")
        sorted_indices = subset.to_numpy()[order]
        sorted_rows = row_values[order]
        audit, marker_values, totals = scan_sparse_matrix(
            spec.path,
            spec.matrix_key,
            sorted_rows,
            marker_genes,
            int(config["matrix_scan_data_chunk"]),
            int(config["matrix_scan_row_chunk"]),
            spec.matrix_state,
            float(transformed_contract["target_sum_expm1"]),
            float(transformed_contract["absolute_tolerance"]),
        )
        if spec.matrix_state == "raw_counts":
            if (totals <= 0).any():
                raise RuntimeError(f"Requested cells in {input_id} contain non-positive raw total counts")
            marker_values = np.log1p(marker_values * (10000.0 / totals[:, None])).astype(np.float32)
        elif spec.matrix_state != "log1p_cp10k_registered":
            raise ValueError(f"Unsupported matrix state for {input_id}: {spec.matrix_state}")
        audit.insert(0, "input_id", input_id)
        matrix_audits.append(audit)
        marker_frame.loc[sorted_indices, marker_genes] = marker_values
        total_counts.loc[sorted_indices] = totals

    requested = request_mask.to_numpy(dtype=bool)
    if marker_frame.loc[requested].isna().any().any() or total_counts.loc[requested].isna().any():
        raise RuntimeError("Expression extraction did not cover every requested cell")
    cd4_genes = config["cd4_helper_genes"]
    treg_genes = config["treg_genes"]
    t_genes = config["t_lineage_genes"]
    cd4_rule = config["cd4_helper_rule"]
    treg_rule = config["treg_rule"]
    marker_requested = marker_frame.loc[requested]
    cd4_support = cd4_genes[1:]
    cd4_flag = (
        marker_requested["CD4"].to_numpy() > float(cd4_rule["cd4_min"])
    ) & (
        (marker_requested[cd4_support].to_numpy() > 0).sum(axis=1) >= int(cd4_rule["support_detected_min"])
    ) & (
        marker_requested[cd4_genes].mean(axis=1).to_numpy() >= float(cd4_rule["program_mean_min"])
    )
    treg_support = treg_genes[1:]
    treg_flag = (
        marker_requested["FOXP3"].to_numpy() > float(treg_rule["foxp3_min"])
    ) & (
        (marker_requested[treg_support].to_numpy() > 0).sum(axis=1) >= int(treg_rule["support_detected_min"])
    ) & (
        marker_requested[treg_genes].mean(axis=1).to_numpy() >= float(treg_rule["program_mean_min"])
    )
    requested_indices = cells.index[requested]
    cells.loc[requested_indices, "cd4_helper_exclusion"] = cd4_flag
    cells.loc[requested_indices, "treg_exclusion"] = treg_flag
    cells.loc[requested_indices, "exclusion_union"] = cd4_flag | treg_flag

    cells["t_lineage_mean"] = np.nan
    cells["trdc_expression"] = np.nan
    cells["trdv_detected_count"] = np.nan
    cells.loc[requested_indices, "t_lineage_mean"] = marker_requested[t_genes].mean(axis=1).to_numpy()
    cells.loc[requested_indices, "trdc_expression"] = marker_requested["TRDC"].to_numpy()
    cells.loc[requested_indices, "trdv_detected_count"] = (
        marker_requested[config["trdv_genes"]].to_numpy() > 0
    ).sum(axis=1)
    return cells, pd.concat(matrix_audits, ignore_index=True)


def assign_nk_strata(cells: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    nk_rows = cells.index[cells["stage1_role"].eq("nk_negative")]
    n_total = nk_rows.size
    fractions = config["nk_sampling"]
    n_gdt = int(round(n_total * float(fractions["gdt_like_fraction"])))
    n_tlike = int(round(n_total * float(fractions["t_like_fraction"])))
    candidate = cells.loc[nk_rows]
    gd_pool = candidate.index[(candidate["trdc_expression"] > 0) & (candidate["trdv_detected_count"] == 0)]
    gd_ranked = (
        cells.loc[gd_pool, ["trdc_expression", "t_lineage_mean", "cell_id"]]
        .sort_values(["trdc_expression", "t_lineage_mean", "cell_id"], ascending=[False, False, True], kind="mergesort")
        .index
    )
    gd_selected = gd_ranked[:n_gdt]
    if gd_selected.size < n_gdt:
        raise RuntimeError(f"TRDC+TRDV- NK pool has {gd_selected.size} cells; {n_gdt} required")
    remaining = nk_rows.difference(gd_selected, sort=False)
    t_ranked = (
        cells.loc[remaining, ["t_lineage_mean", "trdc_expression", "cell_id"]]
        .sort_values(["t_lineage_mean", "trdc_expression", "cell_id"], ascending=[False, False, True], kind="mergesort")
        .index
    )
    t_selected = t_ranked[:n_tlike]
    cells.loc[nk_rows, "nk_sampling_stratum"] = "representative"
    cells.loc[t_selected, "nk_sampling_stratum"] = "t_like_hard"
    cells.loc[gd_selected, "nk_sampling_stratum"] = "gdt_like_TRDCpos_TRDVneg"
    return cells


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total == 0:
        return float("nan"), float("nan")
    p = successes / total
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt((p * (1 - p) + z * z / (4 * total)) / total) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def recall_cost_tables(cells: pd.DataFrame, config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    populations = [PRIMARY_GDT, SILVER, SORTED_WEAK, SORTED_SENSITIVITY]
    for truth in populations:
        truth_cells = cells[cells["truth_class"].eq(truth)]
        for source, frame in truth_cells.groupby("source_gse_id", sort=True):
            total = int(frame.shape[0])
            cd4 = frame["cd4_helper_exclusion"].to_numpy(dtype=bool)
            treg = frame["treg_exclusion"].to_numpy(dtype=bool)
            union = cd4 | treg
            kept = total - int(union.sum())
            lower, upper = wilson_interval(kept, total)
            ceiling = kept / total if total else float("nan")
            floor = float(config["recall_guards"]["per_source_floor"]) if truth == PRIMARY_GDT else float("nan")
            rows.append(
                {
                    "population": truth,
                    "source_gse_id": source,
                    "n_cells": total,
                    "cd4_helper_only": int((cd4 & ~treg).sum()),
                    "treg_only": int((treg & ~cd4).sum()),
                    "overlap": int((cd4 & treg).sum()),
                    "union_excluded": int(union.sum()),
                    "union_excluded_fraction": float(union.mean()) if total else float("nan"),
                    "recall_ceiling": ceiling,
                    "recall_ceiling_wilson_low": lower,
                    "recall_ceiling_wilson_high": upper,
                    "applicable_floor": floor,
                    "margin_above_floor": ceiling - floor if truth == PRIMARY_GDT else float("nan"),
                }
            )
    summary = pd.DataFrame(rows)
    primary = summary[summary["population"].eq(PRIMARY_GDT)].copy()
    macro_ceiling = float(primary["recall_ceiling"].mean())
    macro_floor = float(config["recall_guards"]["macro_floor"])
    macro_row = {
        "population": PRIMARY_GDT,
        "source_gse_id": "SOURCE_MACRO",
        "n_cells": int(primary["n_cells"].sum()),
        "cd4_helper_only": int(primary["cd4_helper_only"].sum()),
        "treg_only": int(primary["treg_only"].sum()),
        "overlap": int(primary["overlap"].sum()),
        "union_excluded": int(primary["union_excluded"].sum()),
        "union_excluded_fraction": float(primary["union_excluded"].sum() / max(primary["n_cells"].sum(), 1)),
        "recall_ceiling": macro_ceiling,
        "recall_ceiling_wilson_low": float(primary["recall_ceiling_wilson_low"].mean()),
        "recall_ceiling_wilson_high": float(primary["recall_ceiling_wilson_high"].mean()),
        "applicable_floor": macro_floor,
        "margin_above_floor": macro_ceiling - macro_floor,
    }
    summary = pd.concat([summary, pd.DataFrame([macro_row])], ignore_index=True)

    group_rows = (
        cells[cells["truth_class"].isin(populations)]
        .groupby(["truth_class", "source_gse_id", "group_key"], observed=True, sort=True)
        .agg(
            n_cells=("cell_id", "size"),
            cd4_helper_excluded=("cd4_helper_exclusion", "sum"),
            treg_excluded=("treg_exclusion", "sum"),
            union_excluded=("exclusion_union", "sum"),
        )
        .reset_index()
    )
    group_rows["union_excluded_fraction"] = group_rows["union_excluded"] / group_rows["n_cells"]
    return summary, group_rows


def split_manifest(cells: pd.DataFrame, config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    primary_sources = list(config["primary_sources"])
    outer_rows: list[dict[str, Any]] = []
    group_rows: list[pd.DataFrame] = []
    balance_rows: list[dict[str, Any]] = []
    for outer_idx, heldout in enumerate(primary_sources):
        for stage in ("stage1", "stage2"):
            role_column = f"{stage}_role"
            if stage == "stage1":
                eligible = cells[role_column].isin(["t_positive", "nk_negative"])
                y = cells[role_column].eq("t_positive").astype(np.int8)
            else:
                eligible = cells[role_column].isin(["positive", "negative"])
                y = cells[role_column].eq("positive").astype(np.int8)
            heldout_mask = cells["source_gse_id"].eq(heldout)
            train_mask = eligible & (~heldout_mask)
            eval_mask = eligible & heldout_mask
            train = cells.loc[train_mask, ["group_key", "cell_id", "source_gse_id"]]
            y_train = y.loc[train_mask].to_numpy(dtype=np.int8)
            groups = train["group_key"].to_numpy(dtype=object)
            splitter = StratifiedGroupKFold(
                n_splits=int(config["inner_folds"]),
                shuffle=True,
                random_state=int(config["random_seed"]) + outer_idx * 10 + (1 if stage == "stage2" else 0),
            )
            assignment = np.full(train.shape[0], -1, dtype=np.int8)
            for fold, (_, val_idx) in enumerate(splitter.split(np.zeros(train.shape[0]), y_train, groups)):
                assignment[val_idx] = fold
            if (assignment < 0).any():
                raise RuntimeError(f"Incomplete inner assignment for {heldout} {stage}")
            assignment_frame = pd.DataFrame(
                {"group_key": groups, "inner_fold": assignment, "y": y_train}
            )
            leakage = assignment_frame.groupby("group_key")["inner_fold"].nunique().max()
            if int(leakage) != 1:
                raise RuntimeError(f"Group leakage in {heldout} {stage}")
            aggregated = (
                assignment_frame.groupby(["group_key", "inner_fold"], sort=True)
                .agg(n_cells=("y", "size"), n_positive=("y", "sum"))
                .reset_index()
            )
            aggregated["n_negative"] = aggregated["n_cells"] - aggregated["n_positive"]
            aggregated.insert(0, "stage", stage)
            aggregated.insert(0, "heldout_source", heldout)
            aggregated.insert(0, "outer_fold_id", f"outer_{outer_idx}_{heldout}")
            group_rows.append(aggregated)

            for fold in range(int(config["inner_folds"])):
                fold_mask = assignment == fold
                balance_rows.append(
                    {
                        "outer_fold_id": f"outer_{outer_idx}_{heldout}",
                        "heldout_source": heldout,
                        "stage": stage,
                        "inner_fold": fold,
                        "n_cells": int(fold_mask.sum()),
                        "n_positive": int(y_train[fold_mask].sum()),
                        "n_negative": int(fold_mask.sum() - y_train[fold_mask].sum()),
                        "n_groups": int(pd.Series(groups[fold_mask]).nunique()),
                    }
                )
            outer_rows.append(
                {
                    "outer_fold_id": f"outer_{outer_idx}_{heldout}",
                    "heldout_source": heldout,
                    "stage": stage,
                    "outer_train_cells": int(train_mask.sum()),
                    "outer_train_positive": int(y.loc[train_mask].sum()),
                    "outer_eval_cells": int(eval_mask.sum()),
                    "outer_eval_positive": int(y.loc[eval_mask].sum()),
                    "outer_train_groups": int(pd.Series(groups).nunique()),
                }
            )
    return pd.DataFrame(outer_rows), pd.concat(group_rows, ignore_index=True), pd.DataFrame(balance_rows)


def feature_coverage(inputs: Iterable[InputSpec], feature_genes: list[str], critical: set[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for spec in inputs:
        if spec.matrix_key is None:
            continue
        with h5py.File(spec.path, "r") as handle:
            var = set(read_var_names(handle))
            n_obs = len(read_obs_index(handle))
        missing = sorted(set(feature_genes) - var)
        missing_critical = sorted(critical - var)
        present = len(feature_genes) - len(missing)
        rows.append(
            {
                "input_id": spec.input_id,
                "source_gse_id": spec.source_id,
                "role": spec.role,
                "n_obs": n_obs,
                "n_vars": len(var),
                "features_present": present,
                "features_expected": len(feature_genes),
                "feature_coverage": present / len(feature_genes),
                "missing_critical_count": len(missing_critical),
                "missing_critical_genes": ";".join(missing_critical),
                "missing_feature_count": len(missing),
                "missing_features": ";".join(missing),
                "hard_feature_gate": spec.hard_feature_gate,
            }
        )
    return pd.DataFrame(rows)


def input_manifest(inputs: list[InputSpec], skip_hashes: bool) -> tuple[pd.DataFrame, dict[str, tuple[int, int]]]:
    rows: list[dict[str, Any]] = []
    initial_stats: dict[str, tuple[int, int]] = {}
    for spec in inputs:
        stat = spec.path.stat()
        initial_stats[spec.input_id] = (stat.st_size, stat.st_mtime_ns)
        actual_hash = "SKIPPED" if skip_hashes else file_sha256(spec.path)
        rows.append(
            {
                "input_id": spec.input_id,
                "source_gse_id": spec.source_id,
                "role": spec.role,
                "path": str(spec.path),
                "matrix_key": spec.matrix_key or "",
                "size_bytes": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
                "sha256": actual_hash,
                "expected_sha256": spec.expected_sha256,
                "expected_hash_match": bool(not spec.expected_sha256 or actual_hash == spec.expected_sha256),
            }
        )
    return pd.DataFrame(rows), initial_stats


def verify_file_state(inputs: list[InputSpec], initial: dict[str, tuple[int, int]]) -> pd.DataFrame:
    rows = []
    for spec in inputs:
        stat = spec.path.stat()
        before_size, before_mtime = initial[spec.input_id]
        rows.append(
            {
                "input_id": spec.input_id,
                "size_before": before_size,
                "size_after": stat.st_size,
                "mtime_ns_before": before_mtime,
                "mtime_ns_after": stat.st_mtime_ns,
                "unchanged": bool(before_size == stat.st_size and before_mtime == stat.st_mtime_ns),
            }
        )
    return pd.DataFrame(rows)


def add_check(rows: list[dict[str, str]], check_id: str, status: str, observed: Any, required: Any, details: str) -> None:
    rows.append(
        {
            "check_id": check_id,
            "status": status,
            "observed": str(observed),
            "required": str(required),
            "details": details,
        }
    )


def build_checks(
    *,
    config: dict[str, Any],
    features: list[str],
    cells: pd.DataFrame,
    coverage: pd.DataFrame,
    matrix_audit: pd.DataFrame,
    join_audit: pd.DataFrame,
    recall: pd.DataFrame,
    file_state: pd.DataFrame,
    input_table: pd.DataFrame,
    outer: pd.DataFrame,
    group_splits: pd.DataFrame,
    nk_audit: pd.DataFrame,
) -> pd.DataFrame:
    checks: list[dict[str, str]] = []
    add_check(checks, "protocol_version", "PASS" if config["protocol_version"] == "1.2" else "FAIL", config["protocol_version"], "1.2", "Frozen protocol version")
    add_check(checks, "feature_count", "PASS" if len(features) == int(config["expected_total_gene_count"]) else "FAIL", len(features), config["expected_total_gene_count"], "Frozen individual-gene universe")
    hard_coverage = coverage[coverage["hard_feature_gate"]]
    failed_coverage = hard_coverage[(hard_coverage["feature_coverage"] < float(config["minimum_feature_coverage"])) | (hard_coverage["missing_critical_count"] > 0)]
    add_check(checks, "training_feature_coverage", "PASS" if failed_coverage.empty else "FAIL", ";".join(failed_coverage["input_id"].astype(str)), f">={config['minimum_feature_coverage']} and all critical genes", "Hard-gated development inputs")
    contract_fail = matrix_audit[~matrix_audit["expression_contract_pass"]]
    transformed = matrix_audit[matrix_audit["configured_matrix_state"].eq("log1p_cp10k_registered")]
    permitted_transformed = str(config["transformed_input_contract"]["permitted_input_id"])
    transformed_scope_ok = transformed["input_id"].tolist() == [permitted_transformed]
    expression_ok = contract_fail.empty and transformed_scope_ok
    add_check(
        checks,
        "expression_input_contract",
        "PASS" if expression_ok else "FAIL",
        ";".join(contract_fail["input_id"].astype(str)) if not contract_fail.empty else transformed["input_id"].tolist(),
        "raw counts, except registered HRA005041 passing full log1p(CP10K) audit",
        "Full sparse data-array and inverse-library-size scan",
    )
    hash_fail = input_table[~input_table["expected_hash_match"]]
    add_check(checks, "registered_hashes", "PASS" if hash_fail.empty else "FAIL", ";".join(hash_fail["input_id"].astype(str)), "all available expected hashes match", "Full SHA-256 unless explicitly skipped")
    add_check(checks, "input_file_state", "PASS" if file_state["unchanged"].all() else "FAIL", int(file_state["unchanged"].sum()), len(file_state), "Size and mtime unchanged after read-only workflow")
    join_ok = bool(join_audit.iloc[0]["mapped_cells"] == 107068 and join_audit.iloc[0]["source_duplicate_keys"] == 0 and join_audit.iloc[0]["canonical_duplicate_keys"] == 0)
    add_check(checks, "gse144469_join", "PASS" if join_ok else "FAIL", int(join_audit.iloc[0]["mapped_cells"]), "107068 unique SRR + barcode mappings", "Raw expression to canonical metadata")
    primary_counts = cells[cells["truth_class"].isin([PRIMARY_GDT, PRIMARY_ABT])].groupby(["source_gse_id", "truth_class"]).size().unstack(fill_value=0)
    primary_ok = all(source in primary_counts.index and primary_counts.loc[source].get(PRIMARY_GDT, 0) > 0 and primary_counts.loc[source].get(PRIMARY_ABT, 0) > 0 for source in config["primary_sources"])
    add_check(checks, "primary_labels_nonempty", "PASS" if primary_ok else "FAIL", primary_counts.to_dict(), "both classes in each primary source", "Expression-independent productive TCR rules")
    conflict_count = int(((cells["truth_class"] == PRIMARY_GDT) & (cells["truth_class"] == PRIMARY_ABT)).sum())
    add_check(checks, "primary_label_overlap", "PASS" if conflict_count == 0 else "FAIL", conflict_count, 0, "Mutually exclusive primary labels")
    sensitivity = cells["truth_class"].isin([SILVER, SORTED_WEAK, SORTED_SENSITIVITY])
    sensitivity_train = int(
        (sensitivity & ((cells["stage1_role"] != "none") | (cells["stage2_role"] != "none"))).sum()
    )
    add_check(checks, "sensitivity_excluded_from_training", "PASS" if sensitivity_train == 0 else "FAIL", sensitivity_train, 0, "Silver and all sorted cohorts are sensitivity-only")
    primary_recall = recall[recall["population"].eq(PRIMARY_GDT)]
    macro = primary_recall[primary_recall["source_gse_id"].eq("SOURCE_MACRO")].iloc[0]
    per_source = primary_recall[~primary_recall["source_gse_id"].eq("SOURCE_MACRO")]
    recall_fail = float(macro["recall_ceiling"]) < float(config["recall_guards"]["macro_floor"]) or (per_source["recall_ceiling"] < float(config["recall_guards"]["per_source_floor"])).any()
    warning = (float(macro["margin_above_floor"]) < float(config["recall_guards"]["warning_margin"])) or (per_source["margin_above_floor"] < float(config["recall_guards"]["warning_margin"])).any()
    add_check(checks, "cd4_treg_recall_ceiling", "FAIL" if recall_fail else ("WARN" if warning else "PASS"), round(float(macro["recall_ceiling"]), 6), f"macro >= {config['recall_guards']['macro_floor']}; each source >= {config['recall_guards']['per_source_floor']}", "Fixed post-model exclusions before Step 2")
    leakage = group_splits.groupby(["outer_fold_id", "stage", "group_key"])["inner_fold"].nunique().max()
    add_check(checks, "group_leakage", "PASS" if int(leakage) == 1 else "FAIL", int(leakage), 1, "No group appears in multiple inner folds")
    outer_eval_ok = (outer["outer_eval_cells"] > 0).all() and (outer["outer_train_cells"] > 0).all()
    add_check(checks, "outer_fold_nonempty", "PASS" if outer_eval_ok else "FAIL", int(outer_eval_ok), True, "Every stage has train and held-out cells")
    strata = cells[cells["stage1_role"].eq("nk_negative")]["nk_sampling_stratum"].value_counts(normalize=True)
    expected = config["nk_sampling"]
    strata_ok = all(abs(float(strata.get(name, 0.0)) - float(expected[key])) <= 1e-5 for name, key in [("representative", "representative_fraction"), ("t_like_hard", "t_like_fraction"), ("gdt_like_TRDCpos_TRDVneg", "gdt_like_fraction")])
    add_check(checks, "nk_stratum_composition", "PASS" if strata_ok else "FAIL", strata.round(6).to_dict(), expected, "All strict NK negatives assigned once")
    final_nk = int(nk_audit.loc[nk_audit["stage"].eq("strict_supplemental_NK_final"), "n_cells"].sum())
    add_check(checks, "supplemental_nk_pool", "PASS" if final_nk > 0 else "FAIL", final_nk, ">0", "Mapped, no-productive-TCR atlas NK negatives")
    archive_ok = (ARCHIVE_DIR / "gdTAI_v4_model.pkl").exists() and not (ROOT / "Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v4.0").exists()
    add_check(checks, "legacy_v4_archived", "PASS" if archive_ok else "FAIL", archive_ok, True, "No overwrite of the earlier experimental artifact")
    return pd.DataFrame(checks)


def dataframe_markdown(df: pd.DataFrame, columns: list[str] | None = None, max_rows: int | None = None) -> str:
    view = df if columns is None else df[columns]
    if max_rows is not None:
        view = view.head(max_rows)
    if view.empty:
        return "_No rows._"

    def format_value(value: Any) -> str:
        if pd.isna(value):
            return ""
        if isinstance(value, (float, np.floating)):
            value = f"{float(value):.4f}"
        return str(value).replace("|", "\\|").replace("\n", "<br>")

    headers = [format_value(column) for column in view.columns]
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    lines.extend(
        "| " + " | ".join(format_value(value) for value in row) + " |"
        for row in view.itertuples(index=False, name=None)
    )
    return "\n".join(lines)


def make_figures(coverage: pd.DataFrame, recall: pd.DataFrame, label_counts: pd.DataFrame) -> None:
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 12, "axes.labelsize": 10})
    primary = recall[(recall["population"] == PRIMARY_GDT) & (recall["source_gse_id"] != "SOURCE_MACRO")]
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    colors = ["#2a9d8f" if value >= 0.7 else "#c0392b" for value in primary["recall_ceiling"]]
    ax.bar(primary["source_gse_id"], primary["recall_ceiling"], color=colors)
    ax.axhline(0.7, color="#c0392b", linestyle="--", linewidth=1.2, label="Per-source floor 0.70")
    ax.axhline(0.8, color="#264653", linestyle=":", linewidth=1.2, label="Macro target 0.80")
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("Exclusion-imposed recall ceiling")
    ax.set_title("Fixed CD4/Treg Exclusion Feasibility on Primary gdT Cells")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "cd4_treg_recall_ceiling_by_source.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    plot_coverage = coverage.sort_values(["hard_feature_gate", "feature_coverage", "input_id"], ascending=[False, True, True])
    fig_height = max(5.0, 0.32 * plot_coverage.shape[0])
    fig, ax = plt.subplots(figsize=(9.2, fig_height))
    colors = np.where(plot_coverage["feature_coverage"] >= 0.9, "#2a9d8f", "#d1495b")
    ax.barh(plot_coverage["input_id"], plot_coverage["feature_coverage"], color=colors)
    ax.axvline(0.9, color="#264653", linestyle="--", linewidth=1.2)
    ax.set_xlim(0, 1.01)
    ax.set_xlabel("Frozen 197-gene coverage")
    ax.set_title("V4 Input Feature Coverage")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "feature_coverage_by_input.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    primary_counts = label_counts[label_counts["truth_class"].isin([PRIMARY_GDT, PRIMARY_ABT, SINGLE_ABT, SILVER])]
    pivot = primary_counts.pivot(index="source_gse_id", columns="truth_class", values="n_cells").fillna(0)
    fig, ax = plt.subplots(figsize=(9.0, 5.2))
    pivot.plot(kind="bar", stacked=True, ax=ax, color=["#457b9d", "#e76f51", "#f4a261", "#8ab17d"][: pivot.shape[1]])
    ax.set_yscale("symlog", linthresh=10)
    ax.set_ylabel("Cells (symlog scale)")
    ax.set_title("Preflight Label Composition")
    ax.legend(title="Truth class", frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "label_composition_by_source.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def render_report(
    *,
    checks: pd.DataFrame,
    coverage: pd.DataFrame,
    recall: pd.DataFrame,
    label_counts: pd.DataFrame,
    join_audit: pd.DataFrame,
    outer: pd.DataFrame,
    input_table: pd.DataFrame,
    matrix_audit: pd.DataFrame,
    summary: dict[str, Any],
) -> tuple[Path, Path, Path]:
    status = summary["overall_status"]
    report_md = LOG_DIR / "gdtai_v4_preflight_summary.md"
    failures = checks[checks["status"].eq("FAIL")]
    warnings = checks[checks["status"].eq("WARN")]
    primary_recall = recall[recall["population"].eq(PRIMARY_GDT)]
    coverage_view = coverage[
        ["input_id", "role", "features_present", "features_expected", "feature_coverage", "missing_critical_count", "hard_feature_gate"]
    ]
    input_view = input_table[["input_id", "role", "size_bytes", "sha256", "expected_hash_match"]].copy()
    input_view["sha256"] = input_view["sha256"].str.slice(0, 16)
    text = f"""# gdTAI V4 Step 1 Preflight and Split Freeze

**Overall status:** {status}<br>
**Protocol:** v{summary['protocol_version']}<br>
**Training started:** No<br>
**Step 2 authorized:** No; supervision approval is still required

## Decision

This workflow performed only the precommitted read-only input, truth-label,
feature, exclusion-feasibility, and grouped-split checks. It did not fit,
calibrate, threshold, compare, or promote a classifier.

Failures block Step 2. Warnings require explicit supervision. The frozen
CD4/Treg thresholds were evaluated exactly as written and were not tuned.

### Blocking failures

{dataframe_markdown(failures, ['check_id', 'observed', 'required', 'details'])}

### Warnings

{dataframe_markdown(warnings, ['check_id', 'observed', 'required', 'details'])}

## All preflight checks

{dataframe_markdown(checks, ['check_id', 'status', 'observed', 'required', 'details'])}

## CD4/Treg recall cost

The recall ceiling is `1 - union_excluded / gdT_primary`. It is the maximum
post-exclusion recall even for a perfect pre-exclusion classifier.

![Recall ceiling](assets/cd4_treg_recall_ceiling_by_source.png)

{dataframe_markdown(primary_recall, ['source_gse_id', 'n_cells', 'cd4_helper_only', 'treg_only', 'overlap', 'union_excluded', 'recall_ceiling', 'applicable_floor', 'margin_above_floor'])}

## Frozen feature coverage

![Feature coverage](assets/feature_coverage_by_input.png)

{dataframe_markdown(coverage_view)}

## Expression input audit

Every stored sparse value was scanned. Raw inputs must be finite, nonnegative,
and integer-like within `1e-6`.

The registered HRA005041 exception must additionally reconstruct a per-cell
library sum of 10,000 from `expm1(X)` within the frozen absolute tolerance.

{dataframe_markdown(matrix_audit, ['input_id', 'matrix_key', 'configured_matrix_state', 'n_obs', 'nnz', 'raw_count_pass', 'transformed_max_abs_deviation', 'transformed_rows_outside_tolerance', 'expression_contract_pass'])}

## Ground-truth audit

![Label composition](assets/label_composition_by_source.png)

{dataframe_markdown(label_counts, max_rows=30)}

## GSE144469 join

{dataframe_markdown(join_audit)}

## Frozen outer folds

{dataframe_markdown(outer)}

## Input identity

{dataframe_markdown(input_view)}

## Canonical artifacts

- Cell-label manifest: `Integrated_dataset/tables/gdT_prediction/gdtai_v4_preflight/cell_label_manifest.csv.gz`
- Inner group assignments: `Integrated_dataset/tables/gdT_prediction/gdtai_v4_preflight/inner_group_split_manifest.csv.gz`
- Full checks: `Integrated_dataset/tables/gdT_prediction/gdtai_v4_preflight/preflight_checks.csv`
- Machine summary: `Integrated_dataset/logs/gdT_prediction/gdtai_v4_preflight/preflight_summary.json`
- Reader HTML: `gdT_prediction/gdtai_v4_preflight/index.html`
- PDF: `gdT_prediction/gdtai_v4_preflight/gdtai_v4_preflight_report.pdf`

## Gate

Stop before Step 2. No model fitting is permitted until the user reviews this
package and explicitly approves the next gate.
"""
    report_md.write_text(text, encoding="utf-8")
    for figure in FIGURE_DIR.glob("*.png"):
        shutil.copy2(figure, STATIC_ASSET_DIR / figure.name)
    html_path = STATIC_DIR / "index.html"
    pdf_path = STATIC_DIR / "gdtai_v4_preflight_report.pdf"
    subprocess.run(
        [
            "pandoc",
            str(report_md),
            "--standalone",
            "--toc",
            "--toc-depth=2",
            "--metadata",
            "title=gdTAI V4 Step 1 Preflight and Split Freeze",
            "--css=../gdtai_v4_precommit/print.css",
            f"--output={html_path}",
        ],
        cwd=ROOT,
        check=True,
    )
    subprocess.run(
        ["weasyprint", str(html_path), str(pdf_path)],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return report_md, html_path, pdf_path


def make_input_specs(config: dict[str, Any]) -> list[InputSpec]:
    paths = config["paths"]
    keys = config["matrix_keys"]
    states = config["matrix_states"]
    expected = config.get("expected_sha256", {})
    specs = [
        InputSpec("current_atlas", resolve_path(paths["current_atlas"]), keys["current_atlas"], "stage1_NK_expression_and_future_atlas_inference", "MULTI_ATLAS", matrix_state=states["current_atlas"]),
        InputSpec("legacy_metadata_reference", resolve_path(paths["legacy_metadata_reference"]), None, "read_only_TCR_and_GSE144469_metadata_reference", "MULTI_LEGACY", hard_feature_gate=False),
        InputSpec("hra005041", resolve_path(paths["hra005041"]), keys["hra005041"], "primary_development", "HRA005041", matrix_state=states["hra005041"]),
        InputSpec("gse144469", resolve_path(paths["gse144469"]), keys["gse144469"], "primary_development", "GSE144469", expected.get("gse144469", ""), matrix_state=states["gse144469"]),
        InputSpec("balf_blood_copd", resolve_path(paths["balf_blood_copd"]), keys["balf_blood_copd"], "primary_development_reused_benchmark", "BALF_BLOOD_COPD", matrix_state=states["balf_blood_copd"]),
        InputSpec("gdt_2020aug_wocov", resolve_path(paths["gdt_2020aug_wocov"]), keys["gdt_2020aug_wocov"], "sorted_sensitivity_only", "GDT_2020AUG_woCOV", hard_feature_gate=False, matrix_state=states["gdt_2020aug_wocov"]),
        InputSpec("maltegdt", resolve_path(paths["maltegdt"]), keys["maltegdt"], "sorted_sensitivity_only", "MalteGDT", hard_feature_gate=False, matrix_state=states["maltegdt"]),
        InputSpec("gdtlung", resolve_path(paths["gdtlung"]), keys["gdtlung"], "sorted_sensitivity_only", "GDTlung2023july_7p", hard_feature_gate=False, matrix_state=states["gdtlung"]),
    ]
    extension = pd.read_csv(resolve_path(paths["extension_manifest"]))
    for row in extension.itertuples(index=False):
        cohort = str(row.cohort_id)
        specs.append(
            InputSpec(
                f"extension_{cohort}",
                Path(row.output_h5ad),
                "X",
                "reduced_feature_schema_sensitivity" if cohort == "GSE169246" else "frozen_negative_stress",
                cohort,
                str(row.output_sha256),
                hard_feature_gate=False,
            )
        )
    return specs


def artifact_checksums(paths: Iterable[Path], output: Path) -> None:
    rows = []
    for path in sorted(set(paths), key=lambda value: str(value)):
        if path.exists() and path.is_file():
            rows.append(f"{file_sha256(path)}  {path.relative_to(ROOT)}")
    output.write_text("\n".join(rows) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    ensure_dirs()
    setup_logging(args.log_level)
    config = json.loads(resolve_path(args.config).read_text(encoding="utf-8"))
    protocol = resolve_path(config["protocol_path"])
    tcr_genes = parse_protocol_tcr_genes(protocol)
    non_tcr_genes = list(config["non_tcr_genes"])
    feature_genes = sorted(set(tcr_genes + non_tcr_genes))
    if len(tcr_genes) != int(config["expected_tcr_gene_count"]):
        raise RuntimeError(f"Protocol TCR gene count is {len(tcr_genes)}, expected {config['expected_tcr_gene_count']}")
    if len(set(non_tcr_genes)) != int(config["expected_non_tcr_gene_count"]):
        raise RuntimeError("Non-TCR feature count does not match the precommit")
    if len(feature_genes) != int(config["expected_total_gene_count"]):
        raise RuntimeError("Total feature count does not match the precommit")

    input_specs = make_input_specs(config)
    input_by_id = {spec.input_id: spec for spec in input_specs}
    for spec in input_specs:
        if not spec.path.exists():
            raise FileNotFoundError(spec.path)
    input_table, initial_stats = input_manifest(input_specs, args.skip_full_hashes)
    coverage = feature_coverage(input_specs, feature_genes, set(config["critical_genes"]))

    legacy = load_legacy_reference(input_by_id["legacy_metadata_reference"].path)
    hra = load_hra_frame(input_by_id["hra005041"].path)
    gse, join_audit = load_gse144469_frame(input_by_id["gse144469"].path, legacy)
    balf = load_balf_frame(input_by_id["balf_blood_copd"].path)
    gdt2020 = load_sorted_frame(input_by_id["gdt_2020aug_wocov"].path, "gdt_2020aug_wocov", "GDT_2020AUG_woCOV")
    malte = load_sorted_frame(input_by_id["maltegdt"].path, "maltegdt", "MalteGDT")
    lung = load_sorted_frame(input_by_id["gdtlung"].path, "gdtlung", "GDTlung2023july_7p")
    atlas_nk, nk_audit = load_atlas_supplemental_nk(
        input_by_id["current_atlas"].path, legacy, set(config["primary_sources"])
    )
    cells = pd.concat([hra, gse, balf, gdt2020, malte, lung, atlas_nk], ignore_index=True)
    del legacy
    if cells["cell_id"].duplicated().any():
        duplicates = cells.loc[cells["cell_id"].duplicated(keep=False), ["cell_id", "source_gse_id"]].head(20)
        raise RuntimeError(f"Development manifest has duplicate cell IDs:\n{duplicates}")
    cells = cells.sort_values(["source_gse_id", "cell_id"], kind="mergesort").reset_index(drop=True)

    marker_genes = sorted(
        set(config["t_lineage_genes"] + config["cd4_helper_genes"] + config["treg_genes"] + ["TRDC"] + config["trdv_genes"])
    )
    cells, matrix_audit = apply_expression_audits(cells, input_by_id, marker_genes, config)
    cells = assign_nk_strata(cells, config)
    recall_table, recall_group_table = recall_cost_tables(cells, config)
    label_counts = (
        cells.groupby(["source_gse_id", "truth_class"], observed=True, sort=True)
        .size()
        .reset_index(name="n_cells")
    )
    outer, group_splits, split_balance = split_manifest(cells, config)
    outer_repeat, group_repeat, balance_repeat = split_manifest(cells, config)
    if not outer.equals(outer_repeat) or not group_splits.equals(group_repeat) or not split_balance.equals(balance_repeat):
        raise RuntimeError("Grouped split construction is not deterministic")
    file_state = verify_file_state(input_specs, initial_stats)
    checks = build_checks(
        config=config,
        features=feature_genes,
        cells=cells,
        coverage=coverage,
        matrix_audit=matrix_audit,
        join_audit=join_audit,
        recall=recall_table,
        file_state=file_state,
        input_table=input_table,
        outer=outer,
        group_splits=group_splits,
        nk_audit=nk_audit,
    )
    overall_status = "FAIL" if checks["status"].eq("FAIL").any() else ("WARN" if checks["status"].eq("WARN").any() else "PASS")

    feature_manifest = pd.DataFrame(
        {
            "feature_index": np.arange(len(feature_genes), dtype=int),
            "gene": feature_genes,
            "feature_class": ["TCR" if gene in set(tcr_genes) else "context" for gene in feature_genes],
        }
    )
    rules = {"cd4_helper_rule": config["cd4_helper_rule"], "treg_rule": config["treg_rule"]}
    feature_manifest.to_csv(TABLE_DIR / "feature_manifest.csv", index=False)
    coverage.to_csv(TABLE_DIR / "feature_coverage.csv", index=False)
    input_table.to_csv(TABLE_DIR / "input_manifest.csv", index=False)
    matrix_audit.to_csv(TABLE_DIR / "raw_count_matrix_audit.csv", index=False)
    file_state.to_csv(TABLE_DIR / "input_file_state_after_preflight.csv", index=False)
    join_audit.to_csv(TABLE_DIR / "gse144469_join_audit.csv", index=False)
    label_counts.to_csv(TABLE_DIR / "label_counts_by_source.csv", index=False)
    recall_table.to_csv(TABLE_DIR / "cd4_treg_recall_cost_by_source.csv", index=False)
    recall_group_table.to_csv(TABLE_DIR / "cd4_treg_recall_cost_by_group.csv", index=False)
    nk_audit.to_csv(TABLE_DIR / "supplemental_nk_pool_audit.csv", index=False)
    outer.to_csv(TABLE_DIR / "outer_fold_manifest.csv", index=False)
    split_balance.to_csv(TABLE_DIR / "inner_split_balance.csv", index=False)
    checks.to_csv(TABLE_DIR / "preflight_checks.csv", index=False)
    deterministic_gzip_csv(group_splits, TABLE_DIR / "inner_group_split_manifest.csv.gz")

    cell_columns = [
        "cell_id", "source_gse_id", "expression_input_id", "expression_row", "truth_class",
        "truth_reliability", "stage1_role", "stage1_weight", "stage2_role", "group_key",
        "group_level", "clonotype_key", "has_any_ab_tcr", "has_any_gd_tcr",
        "has_paired_ab_tcr", "has_paired_gd_tcr", "doublet_flag", "nk_annotation_strength",
        "nk_sampling_stratum", "cd4_helper_exclusion", "treg_exclusion", "exclusion_union",
    ]
    deterministic_gzip_csv(cells[cell_columns], TABLE_DIR / "cell_label_manifest.csv.gz")

    split_hash = file_sha256(TABLE_DIR / "inner_group_split_manifest.csv.gz")
    label_hash = file_sha256(TABLE_DIR / "cell_label_manifest.csv.gz")
    rules_hash = canonical_json_sha256(rules)
    summary = {
        "protocol_version": config["protocol_version"],
        "overall_status": overall_status,
        "step2_authorized": False,
        "training_started": False,
        "random_seed": config["random_seed"],
        "feature_count": len(feature_genes),
        "feature_manifest_sha256": file_sha256(TABLE_DIR / "feature_manifest.csv"),
        "cd4_treg_rule_sha256": rules_hash,
        "cell_label_manifest_sha256": label_hash,
        "inner_group_split_manifest_sha256": split_hash,
        "n_manifest_cells": int(cells.shape[0]),
        "n_preflight_failures": int(checks["status"].eq("FAIL").sum()),
        "n_preflight_warnings": int(checks["status"].eq("WARN").sum()),
        "failed_checks": checks.loc[checks["status"].eq("FAIL"), "check_id"].tolist(),
        "warning_checks": checks.loc[checks["status"].eq("WARN"), "check_id"].tolist(),
    }
    write_json(summary, LOG_DIR / "preflight_summary.json")
    make_figures(coverage, recall_table, label_counts)
    report_md, html_path, pdf_path = render_report(
        checks=checks,
        coverage=coverage,
        recall=recall_table,
        label_counts=label_counts,
        join_audit=join_audit,
        outer=outer,
        input_table=input_table,
        matrix_audit=matrix_audit,
        summary=summary,
    )
    checksum_paths = list(TABLE_DIR.glob("*")) + list(FIGURE_DIR.glob("*.png")) + [
        LOG_DIR / "preflight_summary.json", report_md, html_path, pdf_path,
    ]
    artifact_checksums(checksum_paths, LOG_DIR / "artifact_checksums.sha256")
    logging.info("Preflight complete: status=%s, failures=%s, warnings=%s", overall_status, summary["n_preflight_failures"], summary["n_preflight_warnings"])
    logging.info("Step 2 remains blocked pending supervision")


if __name__ == "__main__":
    main()
