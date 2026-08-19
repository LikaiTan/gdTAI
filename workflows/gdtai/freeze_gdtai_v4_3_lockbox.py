#!/usr/bin/env python3
"""Freeze gdTAI V4.3 external gold and author-NK lockbox cells before scoring."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from workflows.gdtai.run_gdtai_v4_2_ground_truth import axis_node_values, decode  # noqa: E402


ATLAS = Path("/ssd/tnk_phase3/Integrated_dataset/full_atlas/tcr_corrected/integrated_full_atlas.h5ad")
LABELS = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_ground_truth/v4_2_label_manifest.parquet"
OUT_DIR = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_3_lockbox"
LOG_DIR = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_3_lockbox"
NK_RULES = {
    "BALF_BLOOD_COPD": "NK",
    "GSE169246": "NKcell",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(64 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def exact_mask(obs: h5py.Group, name: str, value: str) -> np.ndarray:
    node = obs[name]
    if isinstance(node, h5py.Group) and {"categories", "codes"}.issubset(node):
        categories = axis_node_values(node["categories"]).astype(str)
        matches = np.flatnonzero(categories == value)
        if len(matches) != 1:
            return np.zeros(len(node["codes"]), dtype=bool)
        return np.asarray(node["codes"][:], dtype=np.int64) == matches[0]
    values = axis_node_values(node).astype(str)
    return values == value


def subset_values(obs: h5py.Group, name: str, positions: np.ndarray) -> np.ndarray:
    node = obs[name]
    if isinstance(node, h5py.Group) and {"categories", "codes"}.issubset(node):
        categories = axis_node_values(node["categories"])
        codes = np.asarray(node["codes"][positions], dtype=np.int64)
        output = np.full(len(positions), "", dtype=object)
        present = codes >= 0
        output[present] = categories[codes[present]]
        return output
    if isinstance(node, h5py.Group) and {"values", "mask"}.issubset(node):
        values = subset_values(node, "values", positions)
        mask = np.asarray(node["mask"][positions], dtype=bool)
        values = values.astype(object)
        values[mask] = ""
        return values
    values = np.asarray(node[positions])
    return decode(values) if values.dtype.kind in {"O", "S", "U"} else values


def subset_bool(obs: h5py.Group, name: str, positions: np.ndarray) -> np.ndarray:
    values = subset_values(obs, name, positions)
    if values.dtype.kind == "b":
        return values.astype(bool)
    if values.dtype.kind in "iuf":
        numeric = values.astype(float)
        return np.isfinite(numeric) & (numeric != 0)
    return pd.Series(values).astype(str).str.lower().isin({"true", "1", "yes", "y", "t"}).to_numpy()


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    labels = pd.read_parquet(LABELS)
    receptor_lockbox = labels.loc[
        ~labels.allow_fit & labels.truth_class.isin(["gdT_gold", "abT_gold"])
    ].copy()
    receptor_lockbox["lockbox_label_source"] = receptor_lockbox.label_source

    with h5py.File(ATLAS, "r") as handle:
        obs = handle["obs"]
        index_name = obs.attrs.get("_index", "_index")
        if isinstance(index_name, bytes):
            index_name = index_name.decode()
        nk_rows = []
        audit_rows = []
        for dataset, annotation in NK_RULES.items():
            positions = np.flatnonzero(
                exact_mask(obs, "source_gse_id", dataset)
                & exact_mask(obs, "source_cell_type", annotation)
            )
            any_tcr = np.logical_or.reduce([
                subset_bool(obs, f"has_{chain}", positions)
                for chain in ("TRA", "TRB", "TRG", "TRD")
            ])
            eligible_positions = positions[~any_tcr]
            cell_id = subset_values(obs, str(index_name), eligible_positions)
            donor = subset_values(obs, "donor_id_harmonized_v2", eligible_positions)
            sample = subset_values(obs, "sample_id_harmonized_v2", eligible_positions)
            library = subset_values(obs, "library_id_harmonized_v2", eligible_positions)
            nk_rows.append(pd.DataFrame({
            "cell_id": cell_id,
            "atlas_row": eligible_positions,
            "truth_class": "nk_lockbox",
            "label": 0,
            "stage": "lockbox",
            "reliability_weight": 1.0,
            "label_source": "author_NK_no_productive_TCR",
            "source_gse_id": dataset,
            "donor_id": donor,
            "sample_id": sample,
            "library_id": library,
            "cohort_role": "locked_external",
            "allow_fit": False,
            "allow_threshold_selection": False,
            "has_any_ab_tcr": False,
            "has_any_gd_tcr": False,
            "dual_receptor": False,
            "group_id": [f"{dataset}::{value}" for value in sample],
            "lockbox_label_source": f"source_cell_type={annotation}; no productive TCR",
            }))
            audit_rows.append({
            "source_gse_id": dataset,
            "author_annotation": annotation,
                "author_nk": int(len(positions)),
                "productive_tcr_conflicts": int(any_tcr.sum()),
                "eligible_nk_lockbox": int((~any_tcr).sum()),
            })
    nk = pd.concat(nk_rows, ignore_index=True)
    overlap = set(receptor_lockbox.cell_id) & set(nk.cell_id)
    if overlap:
        raise RuntimeError(f"Lockbox truth conflict in {len(overlap)} cells")
    lockbox = pd.concat([receptor_lockbox, nk], ignore_index=True, sort=False)
    if lockbox.cell_id.duplicated().any():
        raise RuntimeError("Lockbox cell IDs are not unique")
    path = OUT_DIR / "lockbox_manifest.parquet"
    lockbox.to_parquet(path, index=False)
    audit = pd.DataFrame(audit_rows)
    audit.to_csv(OUT_DIR / "author_nk_lockbox_audit.csv", index=False)
    counts = lockbox.groupby(["truth_class", "source_gse_id"], as_index=False).size()
    counts.to_csv(OUT_DIR / "lockbox_counts.csv", index=False)
    result = {
        "status": "PASS_LOCKBOX_FROZEN",
        "n_lockbox_cells": len(lockbox),
        "n_gdt_gold": int(lockbox.truth_class.eq("gdT_gold").sum()),
        "n_abt_gold": int(lockbox.truth_class.eq("abT_gold").sum()),
        "n_author_nk": int(lockbox.truth_class.eq("nk_lockbox").sum()),
        "manifest_sha256": sha256(path),
        "model_scored": False,
        "threshold_selected": False,
    }
    (LOG_DIR / "freeze_summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
