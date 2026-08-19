#!/usr/bin/env python3
"""Freeze the gdTAI V4.2 truth, NK-negative, and lockbox manifests.

This is a read-only analysis. It never writes an H5AD and it fails closed when
receptor support, source-cell mapping, or cohort role is ambiguous.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import html
import io
import json
import subprocess
import tarfile
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "configs/models/gdtai/v4_2_ground_truth.json"
TABLE_DIR = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_ground_truth"
FIG_DIR = ROOT / "Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_ground_truth"
LOG_DIR = ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_ground_truth"
REPORT_DIR = ROOT / "gdT_prediction/gdtai_v4_2_ground_truth"
ASSET_DIR = REPORT_DIR / "assets"
INVALID = {"", "na", "nan", "none", "null", "missing", "unknown"}


def resolve(path: str | Path) -> Path:
    value = Path(path)
    return value if value.is_absolute() else ROOT / value


def sha256_file(path: Path, chunk: int = 64 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def decode(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [value.decode() if isinstance(value, bytes) else str(value) for value in values],
        dtype=object,
    )


def axis_values(group: h5py.Group, name: str) -> np.ndarray:
    node = group[name]
    if isinstance(node, h5py.Group) and {"categories", "codes"}.issubset(node):
        categories = axis_node_values(node["categories"])
        codes = np.asarray(node["codes"][:], dtype=np.int64)
        result = np.full(len(codes), "", dtype=object)
        present = codes >= 0
        result[present] = categories[codes[present]]
        return result
    return axis_node_values(node)


def axis_node_values(node: h5py.Group | h5py.Dataset) -> np.ndarray:
    if isinstance(node, h5py.Group) and {"values", "mask"}.issubset(node):
        values = axis_node_values(node["values"])
        mask = np.asarray(node["mask"][:], dtype=bool)
        if values.dtype.kind in {"O", "S", "U"}:
            values = values.astype(object)
            values[mask] = ""
        else:
            values = values.astype(float)
            values[mask] = np.nan
        return values
    values = np.asarray(node[:])
    return decode(values) if values.dtype.kind in {"O", "S", "U"} else values


def index_values(obs: h5py.Group) -> np.ndarray:
    key = obs.attrs.get("_index", "_index")
    if isinstance(key, bytes):
        key = key.decode()
    return axis_values(obs, str(key))


def bool_values(obs: h5py.Group, name: str, n: int) -> np.ndarray:
    if name not in obs:
        return np.zeros(n, dtype=bool)
    values = axis_values(obs, name)
    if values.dtype.kind == "b":
        return values.astype(bool)
    if values.dtype.kind in "iuf":
        numeric = values.astype(float)
        return np.isfinite(numeric) & (numeric != 0)
    return pd.Series(values).astype(str).str.lower().isin({"true", "1", "yes", "y", "t"}).to_numpy()


def numeric_values(obs: h5py.Group, name: str, n: int) -> np.ndarray:
    if name not in obs:
        return np.full(n, np.nan)
    return pd.to_numeric(pd.Series(axis_values(obs, name)), errors="coerce").to_numpy(float)


def nonempty(values: np.ndarray) -> np.ndarray:
    cleaned = pd.Series(values).astype("string").fillna("").str.strip().str.lower()
    return ~cleaned.isin(INVALID).to_numpy()


def source_positions(source: np.ndarray, name: str) -> np.ndarray:
    return np.flatnonzero(source == name)


def choose_top_chain(frame: pd.DataFrame, chain: str) -> pd.DataFrame:
    local = frame.loc[frame["chain"].eq(chain)].copy()
    if local.empty:
        return pd.DataFrame(columns=["donor_id", "barcode_core", f"{chain}_cdr3", f"{chain}_umis"])
    local["umis"] = pd.to_numeric(local["umis"], errors="coerce")
    local["reads"] = pd.to_numeric(local.get("reads", 0), errors="coerce").fillna(0)
    local = local.sort_values(
        ["donor_id", "barcode_core", "umis", "reads", "contig_id"],
        ascending=[True, True, False, False, True],
    ).drop_duplicates(["donor_id", "barcode_core"])
    return local[["donor_id", "barcode_core", "cdr3", "umis"]].rename(
        columns={"cdr3": f"{chain}_cdr3", "umis": f"{chain}_umis"}
    )


def build_gse144_overlay(config: dict, atlas_source: np.ndarray, atlas_original: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame]:
    spec = config["gse144469"]
    selected = pd.read_csv(resolve(spec["selected_inputs"]))
    selected = selected[selected["sample_key"].str.endswith("-cd3")].copy()
    selected["gsm"] = selected["path"].str.extract(r"(GSM\d+)")
    selected["donor_id"] = selected["sample_key"].str.replace("-cd3", "", regex=False).str.upper()
    runs = pd.read_csv(resolve(spec["sra_run_table"]))
    run_to_donor = runs[["Run", "GEO_Accession (exp)"]].merge(
        selected[["gsm", "donor_id"]], left_on="GEO_Accession (exp)", right_on="gsm", how="inner"
    ).set_index("Run")["donor_id"].to_dict()

    ab = pd.read_csv(resolve(spec["ab_contigs"]), compression="gzip", low_memory=False)
    ab = ab[ab["productive"].astype(str).str.lower().eq("true") & ab["chain"].isin(["TRA", "TRB"])].copy()
    ab["donor_id"] = ab["barcode"].str.rsplit("-", n=1).str[-1].str.upper()
    ab["barcode_core"] = ab["barcode"].str.rsplit("-", n=1).str[0] + "-1"
    frames = [ab]
    with tarfile.open(resolve(spec["gd_contig_tar"])) as archive:
        for member in archive.getmembers():
            if "gdTCR_filtered_contig_annotations.csv" not in member.name:
                continue
            donor = Path(member.name).name.split("-gdTCR", 1)[0].split("_", 1)[-1].upper()
            extracted = archive.extractfile(member)
            if extracted is None:
                continue
            payload = extracted.read()
            if member.name.endswith(".gz"):
                payload = gzip.decompress(payload)
            local = pd.read_csv(io.BytesIO(payload), low_memory=False)
            local = local[local["productive"].astype(str).str.lower().eq("true") & local["chain"].isin(["TRG", "TRD"])].copy()
            local["donor_id"] = donor
            local["barcode_core"] = local["barcode"].astype(str)
            frames.append(local)
    contigs = pd.concat(frames, ignore_index=True, sort=False)
    required = ["donor_id", "barcode_core", "chain", "cdr3", "umis", "reads", "contig_id"]
    contigs = contigs[required]
    cells: pd.DataFrame | None = None
    for chain in ("TRA", "TRB", "TRG", "TRD"):
        top = choose_top_chain(contigs, chain)
        cells = top if cells is None else cells.merge(top, on=["donor_id", "barcode_core"], how="outer")
    assert cells is not None

    positions = source_positions(atlas_source, "GSE144469")
    atlas = pd.DataFrame({"atlas_row": positions, "original_cell_id": atlas_original[positions]})
    atlas["run"] = atlas["original_cell_id"].str.extract(r"^(SRR\d+)_")
    atlas["donor_id"] = atlas["run"].map(run_to_donor)
    atlas["barcode_core"] = atlas["original_cell_id"].str.extract(r"_(.{16}-1)$")
    joined = atlas.merge(cells, on=["donor_id", "barcode_core"], how="left", validate="many_to_one")
    for chain in ("TRA", "TRB", "TRG", "TRD"):
        joined[f"has_{chain}"] = joined[f"{chain}_cdr3"].fillna("").ne("")
    joined["paired_ab"] = joined["has_TRA"] & joined["has_TRB"]
    joined["paired_gd"] = joined["has_TRG"] & joined["has_TRD"]
    joined["any_ab"] = joined["has_TRA"] | joined["has_TRB"]
    joined["any_gd"] = joined["has_TRG"] | joined["has_TRD"]
    joined["ab_umi2"] = joined["paired_ab"] & joined["TRA_umis"].ge(config["umi_min_per_defining_chain"]) & joined["TRB_umis"].ge(config["umi_min_per_defining_chain"])
    joined["gd_umi2"] = joined["paired_gd"] & joined["TRG_umis"].ge(config["umi_min_per_defining_chain"]) & joined["TRD_umis"].ge(config["umi_min_per_defining_chain"])
    joined["is_cd3_library"] = joined["donor_id"].notna()
    audit = pd.DataFrame([
        {"metric": "atlas_gse144469_cells", "value": len(joined)},
        {"metric": "exact_cd3_library_cells", "value": int(joined["is_cd3_library"].sum())},
        {"metric": "paired_gd_before_umi", "value": int((joined["paired_gd"] & ~joined["any_ab"]).sum())},
        {"metric": "paired_gd_both_umi_ge2", "value": int((joined["gd_umi2"] & ~joined["any_ab"]).sum())},
        {"metric": "paired_ab_before_umi", "value": int((joined["paired_ab"] & ~joined["any_gd"]).sum())},
        {"metric": "paired_ab_both_umi_ge2", "value": int((joined["ab_umi2"] & ~joined["any_gd"]).sum())},
        {"metric": "dual_receptor", "value": int((joined["any_ab"] & joined["any_gd"]).sum())},
        {"metric": "cd45_cells_with_assigned_tcr", "value": int((~joined["is_cd3_library"] & (joined["any_ab"] | joined["any_gd"])).sum())},
    ])
    return joined, audit


def source_umi_status(source: np.ndarray, paired: np.ndarray, a_available: np.ndarray, b_available: np.ndarray) -> pd.DataFrame:
    rows = []
    for name in sorted(set(source)):
        mask = (source == name) & paired
        n = int(mask.sum())
        both = a_available[mask] & b_available[mask]
        if n == 0:
            status = "JOIN_UNRESOLVED"
        elif both.all():
            status = "AVAILABLE"
        elif (~a_available[mask] & ~b_available[mask]).all():
            status = "UNAVAILABLE"
        else:
            status = "MIXED"
        rows.append({"source_gse_id": name, "n_paired_calls": n, "n_both_chain_umi_available": int(both.sum()), "umi_status": status})
    return pd.DataFrame(rows)


def qualifying_pair(source: np.ndarray, paired: np.ndarray, a_available: np.ndarray, b_available: np.ndarray,
                    a_umi: np.ndarray, b_umi: np.ndarray, status: pd.DataFrame, minimum: int) -> np.ndarray:
    lookup = status.set_index("source_gse_id")["umi_status"].to_dict()
    result = np.zeros(len(source), dtype=bool)
    for name, state in lookup.items():
        mask = (source == name) & paired
        if state == "UNAVAILABLE":
            result[mask] = True
        elif state in {"AVAILABLE", "MIXED"}:
            result[mask] = a_available[mask] & b_available[mask] & (a_umi[mask] >= minimum) & (b_umi[mask] >= minimum)
    return result


def source_annotation_mask(path: Path, column: str, mode: str, values: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        ids = index_values(obs)
        annotation = pd.Series(axis_values(obs, column)).astype(str)
        if mode == "exact":
            mask = annotation.isin(values).to_numpy()
        elif mode == "prefix":
            mask = np.logical_or.reduce([annotation.str.startswith(value).to_numpy() for value in values])
        else:
            raise ValueError(f"Unsupported match mode: {mode}")
        return ids, mask, annotation.to_numpy(), {"n_source_cells": len(ids), "n_author_nk": int(mask.sum())}


def doublet_mask(path: Path, columns: list[str]) -> tuple[np.ndarray, list[str]]:
    with h5py.File(path, "r") as handle:
        obs = handle["obs"]
        n = len(index_values(obs))
        result = np.zeros(n, dtype=bool)
        found = []
        for column in columns:
            if column not in obs:
                continue
            found.append(column)
            values = axis_values(obs, column)
            if values.dtype.kind == "b":
                result |= values.astype(bool)
            else:
                result |= pd.Series(values).astype(str).str.lower().isin({"true", "1", "yes", "doublet"}).to_numpy()
    return result, found


def frozen_tier1_rows(staged_path: Path, atlas_source: np.ndarray, atlas_original: np.ndarray,
                      atlas_source_obs: np.ndarray) -> np.ndarray:
    with h5py.File(staged_path, "r") as handle:
        obs = handle["obs"]
        primary = bool_values(obs, "primary_nk_anchor", len(index_values(obs)))
        staged_source = axis_values(obs, "source_gse_id")[primary]
        staged_original = axis_values(obs, "original_cell_id")[primary]
        staged_cohort = axis_values(obs, "input_cohort_id")[primary]
    mapped = []
    for name in sorted(set(staged_source)):
        local = staged_source == name
        atlas_pos = source_positions(atlas_source, name)
        key_values = atlas_source_obs[atlas_pos] if np.all(staged_cohort[local] == "current_atlas") else atlas_original[atlas_pos]
        lookup = dict(zip(key_values, atlas_pos))
        rows = [lookup.get(value, -1) for value in staged_original[local]]
        if min(rows, default=0) < 0:
            raise RuntimeError(f"Incomplete frozen primary NK mapping for {name}")
        mapped.extend(rows)
    return np.asarray(mapped, dtype=np.int64)


def save_figure(fig: plt.Figure, name: str) -> None:
    target = FIG_DIR / name
    fig.savefig(target, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    (ASSET_DIR / name).write_bytes(target.read_bytes())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--skip-hash", action="store_true")
    args = parser.parse_args()
    config = json.loads(resolve(args.config).read_text())
    for path in (TABLE_DIR, FIG_DIR, LOG_DIR, REPORT_DIR, ASSET_DIR):
        path.mkdir(parents=True, exist_ok=True)
    atlas_path = resolve(config["atlas_path"])
    observed_sha = config["atlas_sha256"] if args.skip_hash else sha256_file(atlas_path)
    if observed_sha != config["atlas_sha256"]:
        raise RuntimeError("Atlas checksum mismatch")

    with h5py.File(atlas_path, "r") as handle:
        obs = handle["obs"]
        cell_id = index_values(obs)
        n = len(cell_id)
        source = axis_values(obs, "source_gse_id")
        original = axis_values(obs, "original_cell_id")
        source_obs = axis_values(obs, "source_obs_name")
        donor = axis_values(obs, "donor_id_harmonized_v2")
        sample = axis_values(obs, "sample_id_harmonized_v2")
        library = axis_values(obs, "library_id_harmonized_v2")
        has = {chain: bool_values(obs, f"has_{chain}", n) for chain in ("TRA", "TRB", "TRG", "TRD")}
        paired_ab = bool_values(obs, "has_TRA_TRB_paired", n)
        paired_gd = bool_values(obs, "has_TRG_TRD_paired", n)
        any_ab = bool_values(obs, "has_any_ab_tcr", n)
        any_gd = bool_values(obs, "has_any_gd_tcr", n)
        availability = {chain: bool_values(obs, f"{chain}_umi_available_rebuilt_v2", n) for chain in ("TRA", "TRB", "TRG", "TRD")}
        umis = {chain: numeric_values(obs, f"{chain}_umis", n) for chain in ("TRA", "TRB", "TRG", "TRD")}
        # Unaffected sources retain native per-chain UMI values but do not carry
        # a rebuilt-v2 availability flag. Positive native counts are explicit
        # quantitative support; zero/null remains unavailable, not UMI=0.
        for chain in ("TRA", "TRB", "TRG", "TRD"):
            availability[chain] |= has[chain] & np.isfinite(umis[chain]) & (umis[chain] > 0)

    gse144, gse144_audit = build_gse144_overlay(config, source, original)
    gse144.to_parquet(TABLE_DIR / "gse144469_exact_tcr_overlay.parquet", index=False, compression="zstd")
    gse144_audit.to_csv(TABLE_DIR / "gse144469_join_audit.csv", index=False)
    rows144 = gse144["atlas_row"].to_numpy(np.int64)
    for chain in ("TRA", "TRB", "TRG", "TRD"):
        has[chain][rows144] = gse144[f"has_{chain}"].to_numpy(bool)
        availability[chain][rows144] = gse144[f"has_{chain}"].to_numpy(bool)
        umis[chain][rows144] = gse144[f"{chain}_umis"].to_numpy(float)
    paired_ab[rows144] = gse144["paired_ab"].to_numpy(bool)
    paired_gd[rows144] = gse144["paired_gd"].to_numpy(bool)
    any_ab[rows144] = gse144["any_ab"].to_numpy(bool)
    any_gd[rows144] = gse144["any_gd"].to_numpy(bool)

    ab_status = source_umi_status(source, paired_ab, availability["TRA"], availability["TRB"])
    gd_status = source_umi_status(source, paired_gd, availability["TRG"], availability["TRD"])
    umi_status = ab_status.merge(gd_status, on="source_gse_id", suffixes=("_ab", "_gd"))
    umi_status.to_csv(TABLE_DIR / "tcr_umi_status_by_source.csv", index=False)
    q_ab = qualifying_pair(source, paired_ab, availability["TRA"], availability["TRB"], umis["TRA"], umis["TRB"], ab_status, config["umi_min_per_defining_chain"])
    q_gd = qualifying_pair(source, paired_gd, availability["TRG"], availability["TRD"], umis["TRG"], umis["TRD"], gd_status, config["umi_min_per_defining_chain"])

    roles = pd.read_csv(resolve(config["cohort_roles"]))
    locked_sources = set(config["locked_positive_sources"]) | set(config["extension_lockbox_sources"])
    sorted_sources = set(config["sorted_positive_sources"])
    development_positive = set(config["development_positive_sources"])
    fit_allowed = ~np.isin(source, list(locked_sources))
    sorted_positive = np.isin(source, list(sorted_sources))
    receptor_positive = q_gd & ~any_ab & np.isin(source, list(development_positive | {"BALF_BLOOD_COPD"}))
    gold_gdt = sorted_positive | receptor_positive
    gold_ab = q_ab & ~any_gd & ~gold_gdt
    dual = any_ab & any_gd
    single_ab = (has["TRA"] ^ has["TRB"]) & ~any_gd & ~gold_gdt

    nk_rules = pd.read_csv(resolve(config["nk_rules"])).fillna("")
    tier1_rows = frozen_tier1_rows(resolve(config["staged_nk_path"]), source, original, source_obs)
    tier1_rows = tier1_rows[~any_ab[tier1_rows] & ~any_gd[tier1_rows]]
    nk_records = []
    nk_audit = []
    for row in nk_rules.itertuples(index=False):
        if row.tier == "tier1":
            selected_rows = tier1_rows[source[tier1_rows] == row.source_gse_id]
            nk_audit.append({"source_gse_id": row.source_gse_id, "tier": row.tier, "n_source_cells": int((source == row.source_gse_id).sum()), "n_author_nk": len(selected_rows), "n_tcr_conflicts": 0, "n_doublets": 0, "n_eligible": len(selected_rows), "doublet_columns_found": "frozen_anchor"})
        else:
            path = resolve(row.source_h5ad)
            values = str(row.match_values).split("|")
            ids, author_mask, _, counts = source_annotation_mask(path, row.annotation_column, row.match_mode, values)
            doublets, found = doublet_mask(path, [x for x in str(row.doublet_columns).split("|") if x])
            atlas_pos = source_positions(source, row.source_gse_id)
            lookup = dict(zip(original[atlas_pos], atlas_pos))
            source_selected = np.flatnonzero(author_mask)
            mapped_all = np.asarray([lookup.get(ids[i], -1) for i in source_selected], dtype=np.int64)
            retained = mapped_all >= 0
            mapped = mapped_all[retained]
            if len(set(mapped)) != len(mapped):
                raise RuntimeError(f"Nonunique author-NK mapping for {row.source_gse_id}")
            tcr_conflict = any_ab[mapped] | any_gd[mapped]
            source_doublet = doublets[source_selected][retained]
            eligible = ~(tcr_conflict | source_doublet)
            selected_rows = mapped[eligible]
            nk_audit.append({"source_gse_id": row.source_gse_id, "tier": row.tier, **counts, "n_not_retained_in_atlas": int((~retained).sum()), "n_tcr_conflicts": int(tcr_conflict.sum()), "n_doublets": int(source_doublet.sum()), "n_eligible": int(eligible.sum()), "doublet_columns_found": "|".join(found)})
        for atlas_row in selected_rows:
            nk_records.append({"atlas_row": atlas_row, "truth_class": f"nk_{row.tier}", "label": 0, "stage": "stage1", "reliability_weight": float(row.reliability_weight), "label_source": f"author_nk_{row.tier}"})
    nk_audit_df = pd.DataFrame(nk_audit)
    nk_audit_df.to_csv(TABLE_DIR / "nk_author_curation_audit.csv", index=False)

    records = []
    def add(mask: np.ndarray, truth: str, label: int, stage: str, weight: float, provenance: str) -> None:
        for atlas_row in np.flatnonzero(mask):
            records.append({"atlas_row": atlas_row, "truth_class": truth, "label": label, "stage": stage, "reliability_weight": weight, "label_source": provenance})

    add(gold_gdt, "gdT_gold", 1, "both", 1.0, "sorted_or_qualified_paired_gd")
    add(gold_ab, "abT_gold", 0, "both", 1.0, "qualified_paired_ab")
    add(single_ab & fit_allowed, "single_ab_support", 0, "stage2", config["single_ab_support_weight"], "single_productive_TRA_or_TRB")
    records.extend(nk_records)
    manifest = pd.DataFrame(records)
    if manifest.empty:
        raise RuntimeError("Empty label manifest")
    positions = manifest["atlas_row"].to_numpy(np.int64)
    manifest.insert(0, "cell_id", cell_id[positions])
    manifest["source_gse_id"] = source[positions]
    manifest["donor_id"] = donor[positions]
    manifest["sample_id"] = sample[positions]
    manifest["library_id"] = library[positions]
    manifest["cohort_role"] = np.where(np.isin(manifest["source_gse_id"], list(locked_sources)), "lockbox", "development")
    manifest["allow_fit"] = manifest["cohort_role"].eq("development")
    manifest["allow_threshold_selection"] = manifest["allow_fit"]
    manifest["has_any_ab_tcr"] = any_ab[positions]
    manifest["has_any_gd_tcr"] = any_gd[positions]
    manifest["dual_receptor"] = dual[positions]
    manifest["group_id"] = manifest["donor_id"].where(manifest["donor_id"].astype(str).str.lower().isin(INVALID) == False, manifest["sample_id"])
    manifest["group_id"] = manifest["group_id"].where(manifest["group_id"].astype(str).str.lower().isin(INVALID) == False, manifest["library_id"])
    manifest["group_id"] = manifest["source_gse_id"].astype(str) + "::" + manifest["group_id"].astype(str)

    duplicate_roles = manifest.groupby("atlas_row")["truth_class"].nunique()
    conflicting_rows = duplicate_roles[duplicate_roles > 1].index
    # Stage-1/Stage-2 reuse of the same T cell is intentional; incompatible biological labels are not.
    biological = manifest.groupby("atlas_row")["label"].nunique()
    label_conflicts = biological[biological > 1]
    if len(label_conflicts):
        raise RuntimeError(f"{len(label_conflicts)} cells have conflicting labels")
    if manifest.loc[manifest["allow_fit"], "source_gse_id"].isin(locked_sources).any():
        raise RuntimeError("Lockbox leakage into fitting manifest")
    if manifest.loc[manifest["allow_fit"], "dual_receptor"].any():
        raise RuntimeError("Dual receptor cells leaked into fitting manifest")
    nk_source_count = manifest.loc[manifest["truth_class"].isin(["nk_tier1", "nk_tier2"]), "source_gse_id"].nunique()
    if nk_source_count < config["minimum_nk_sources"]:
        raise RuntimeError("Insufficient independent NK sources")
    if not {0, 1}.issubset(set(manifest.loc[manifest["allow_fit"], "label"])):
        raise RuntimeError("Development labels do not contain both classes")

    manifest.to_parquet(TABLE_DIR / "v4_2_label_manifest.parquet", index=False, compression="zstd")
    summary = (manifest.groupby(["cohort_role", "source_gse_id", "truth_class", "stage", "label_source"], dropna=False)
               .agg(n_cells=("atlas_row", "size"), effective_weight=("reliability_weight", "sum"), n_groups=("group_id", "nunique"))
               .reset_index())
    summary.to_csv(TABLE_DIR / "truth_counts_by_source.csv", index=False)
    roles.to_csv(TABLE_DIR / "frozen_cohort_roles.csv", index=False)
    checks = pd.DataFrame([
        {"check": "atlas_checksum", "pass": observed_sha == config["atlas_sha256"], "detail": observed_sha},
        {"check": "gse144_cd45_unassigned", "pass": int(gse144_audit.loc[gse144_audit.metric.eq("cd45_cells_with_assigned_tcr"), "value"].iloc[0]) == 0, "detail": "exact donor+barcode join"},
        {"check": "no_label_conflicts", "pass": len(label_conflicts) == 0, "detail": str(len(label_conflicts))},
        {"check": "no_lockbox_fit", "pass": not manifest.loc[manifest.allow_fit, "source_gse_id"].isin(locked_sources).any(), "detail": ",".join(sorted(locked_sources))},
        {"check": "no_dual_receptor_fit", "pass": not manifest.loc[manifest.allow_fit, "dual_receptor"].any(), "detail": str(int(manifest.loc[manifest.allow_fit, "dual_receptor"].sum()))},
        {"check": "silver_excluded", "pass": not manifest.truth_class.str.contains("silver", case=False).any(), "detail": "zero fitting and threshold weight"},
        {"check": "nk_source_gate", "pass": nk_source_count >= config["minimum_nk_sources"], "detail": f"{nk_source_count} independent sources"},
        {"check": "both_fit_classes", "pass": {0, 1}.issubset(set(manifest.loc[manifest.allow_fit, "label"])), "detail": str(sorted(set(manifest.loc[manifest.allow_fit, "label"])))},
    ])
    checks.to_csv(TABLE_DIR / "step0_checks.csv", index=False)
    status = "PASS_LABEL_MANIFEST_FROZEN" if checks["pass"].all() else "FAIL"

    plot = summary.groupby(["source_gse_id", "truth_class"])["n_cells"].sum().unstack(fill_value=0)
    plot = plot.loc[plot.sum(axis=1).sort_values().index]
    fig, ax = plt.subplots(figsize=(10.5, max(5.5, len(plot) * 0.28)))
    colors = {"gdT_gold": "#007C83", "abT_gold": "#E1A23A", "single_ab_support": "#8D99AE", "nk_tier1": "#B23A48", "nk_tier2": "#5E548E"}
    left = np.zeros(len(plot))
    for column in plot.columns:
        ax.barh(plot.index, plot[column], left=left, color=colors.get(column, "#BBBBBB"), label=column)
        left += plot[column].to_numpy()
    ax.set_xscale("symlog", linthresh=100); ax.set_xlabel("Labeled cells (symlog scale)"); ax.set_title("Frozen V4.2 truth by source")
    ax.grid(axis="x", color="#E5E7EB"); ax.legend(frameon=False, ncol=2); ax.spines[["top", "right", "left"]].set_visible(False)
    save_figure(fig, "truth_by_source.png")

    nk_plot = nk_audit_df.sort_values("n_eligible")
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    ax.barh(nk_plot.source_gse_id, nk_plot.n_author_nk, color="#D7DCE2", label="Author/frozen NK")
    ax.barh(nk_plot.source_gse_id, nk_plot.n_eligible, color=np.where(nk_plot.tier.eq("tier1"), "#B23A48", "#5E548E"), label="Eligible after exclusions")
    ax.set_xlabel("Cells"); ax.set_title("NK negatives span seven independent datasets"); ax.grid(axis="x", color="#E5E7EB")
    ax.spines[["top", "right", "left"]].set_visible(False); ax.legend(frameon=False)
    save_figure(fig, "nk_curation_by_source.png")

    result = {
        "status": status,
        "atlas_path": str(atlas_path),
        "atlas_sha256": observed_sha,
        "n_atlas_cells": n,
        "n_manifest_rows": len(manifest),
        "n_unique_labeled_cells": int(manifest.atlas_row.nunique()),
        "n_development_rows": int(manifest.allow_fit.sum()),
        "n_lockbox_rows": int((~manifest.allow_fit).sum()),
        "n_nk_sources": int(nk_source_count),
        "n_tier1_nk": int((manifest.truth_class == "nk_tier1").sum()),
        "n_tier2_nk": int((manifest.truth_class == "nk_tier2").sum()),
        "gse144469": dict(zip(gse144_audit.metric, gse144_audit.value.astype(int))),
        "h5ad_modified": False,
        "model_fitted": False,
    }
    (LOG_DIR / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")

    nk_rows = "".join(f"<tr><td>{html.escape(str(r.source_gse_id))}</td><td>{r.tier}</td><td>{r.n_author_nk:,}</td><td>{r.n_tcr_conflicts:,}</td><td>{r.n_doublets:,}</td><td>{r.n_eligible:,}</td></tr>" for r in nk_audit_df.itertuples(index=False))
    report = f"""<!doctype html><html><head><meta charset='utf-8'><title>gdTAI V4.2 ground truth</title><style>
    @page{{size:A4 landscape;margin:10mm}}body{{font:14px/1.45 Arial;color:#18232B;margin:0}}main{{max-width:1180px;margin:auto}}header{{background:#163A43;color:white;padding:34px}}section{{padding:24px 32px;border-bottom:1px solid #D9E0E3;page-break-inside:avoid}}h1,h2{{margin-top:0}}h2{{color:#163A43}}.grid{{display:grid;grid-template-columns:repeat(4,1fr);gap:10px}}.card{{border-left:4px solid #00A6A6;background:#F2F7F7;padding:12px}}.big{{font-size:24px;font-weight:700}}img{{width:100%;max-height:500px;object-fit:contain}}table{{width:100%;border-collapse:collapse;font-size:11px}}th,td{{border:1px solid #D5DDE0;padding:6px;text-align:right}}th{{background:#163A43;color:white}}th:first-child,td:first-child{{text-align:left}}.note{{border-left:4px solid #E1A23A;background:#FFF8EA;padding:12px}}.page{{page-break-before:always}}</style></head><body><main>
    <header><h1>gdTAI V4.2 Ground-Truth Freeze</h1><p>Expression-independent labels, expanded NK controls, and broad external lockbox</p></header>
    <section><div class='grid'><div class='card'><div class='big'>{result['n_unique_labeled_cells']:,}</div>unique labeled cells</div><div class='card'><div class='big'>{result['n_nk_sources']}</div>independent NK sources</div><div class='card'><div class='big'>{result['n_tier1_nk']:,}</div>Tier 1 NK</div><div class='card'><div class='big'>{result['n_tier2_nk']:,}</div>Tier 2 NK</div></div><p class='note'><b>{status}</b>. No classifier was fitted and no H5AD was changed.</p></section>
    <section><h2>Ground Truth Definition</h2><p>Primary gdT truth is a sorted gdT cell or an exclusive paired productive TRG/TRD cell. Primary alpha-beta truth is exclusive paired productive TRA/TRB. When per-chain UMI support exists, both defining chains require at least two UMIs; one-UMI calls are excluded. When public/source receptor tables genuinely lack UMI fields, paired productive calls remain gold. Mixed-support sources accept only cells with support on both defining chains. Dual-receptor and silver cells have zero fitting and threshold-selection weight.</p><p>Single productive TRA or TRB calls are retained only as 0.5-weight Stage-2 support negatives. They do not replace paired alpha-beta gold.</p></section>
    <section class='page'><h2>Truth and Cohort Roles</h2><img src='assets/truth_by_source.png'><p>GDT_2020, GDTlung, BALF_BLOOD_COPD, and every new extension cohort are locked. They cannot influence features, fitting, calibration, thresholds, or model-family selection.</p></section>
    <section class='page'><h2>NK Definition</h2><img src='assets/nk_curation_by_source.png'><table><thead><tr><th>Source</th><th>Tier</th><th>Author/frozen NK</th><th>TCR conflicts</th><th>Doublets</th><th>Eligible</th></tr></thead><tbody>{nk_rows}</tbody></table><p>Tier 1 requires the frozen dual-annotation anchor and no corrected productive TCR. Tier 2 requires an original author NK label, no productive TRA/TRB/TRG/TRD, and no explicit doublet. Tier 2 is weighted 0.5, source-balanced, and capped during fitting.</p></section>
    <section class='page'><h2>GSE144469 Repair</h2><p>The raw alpha-beta aggregate and 21 public gdT-contig files were re-parsed. CD3 GEX libraries are mapped by the official SRR-to-GSM table, donor code, and exact barcode; CD45 GEX libraries receive no TCR assignment. The resulting exclusive UMI-supported counts are {result['gse144469']['paired_gd_both_umi_ge2']:,} paired gdT and {result['gse144469']['paired_ab_both_umi_ge2']:,} paired alpha-beta cells. One-UMI exclusions are retained in the audit, not in gold truth.</p></section>
    <section><h2>Next Gate</h2><p>The next step is to freeze the feature and grouped-fold manifests, build the on-disk expression cache, and run fold-local implementation tests. Lockbox data remain untouched until the final model and both operating thresholds are frozen.</p></section>
    </main></body></html>"""
    html_path = REPORT_DIR / "index.html"
    html_path.write_text(report, encoding="utf-8")
    pdf_path = REPORT_DIR / "gdtai_v4_2_ground_truth_report.pdf"
    subprocess.run(["google-chrome", "--headless", "--no-sandbox", "--disable-gpu", "--disable-crash-reporter", f"--print-to-pdf={pdf_path}", f"file://{html_path}"], check=True)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
