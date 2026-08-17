#!/usr/bin/env python3
"""Audit unexpectedly frequent productive TRA/TRB calls in the V4.2 NK boundary.

This is a read-only analysis sidecar. It tests four explanations separately:
low-support/ambient VDJ calls, inadequate transcriptomic representation,
incorrect upstream TCR joins, and a residual population compatible with true
alpha-beta T cells. No H5AD or model artifact is modified.
"""

from __future__ import annotations

import argparse
import html
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable, Sequence

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import anndata as ad  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402
from scipy import sparse, stats  # noqa: E402
from sklearn.metrics import normalized_mutual_info_score, roc_auc_score  # noqa: E402

from workflows.gdtai.gdtai_v4_2_integration_core import (  # noqa: E402
    atomic_write_json,
    h5ad_obs_frame,
    normalize_text,
    resolve,
    sha256_file,
)
from workflows.gdtai.run_gdtai_v4_2_nk_reference_integration import (  # noqa: E402
    development_roles,
)
from workflows.gdtai.run_gdtai_v4_2_tnk_reintegration import (  # noqa: E402
    load_runtime_config,
    marker_genes,
    stage_paths,
)


DEFAULT_CONFIG = PROJECT_ROOT / "configs/models/gdtai/v4_2_tnk_reintegration.json"
SLUG = "gdtai_v4_2_tcr_conflict_audit"
TABLE_DIR = PROJECT_ROOT / "Integrated_dataset/tables/gdT_prediction" / SLUG
FIGURE_DIR = PROJECT_ROOT / "Integrated_dataset/figures/gdT_prediction" / SLUG
LOG_DIR = PROJECT_ROOT / "Integrated_dataset/logs/gdT_prediction" / SLUG
STATIC_DIR = PROJECT_ROOT / "gdT_prediction" / SLUG
SCRIPT_PATH = Path(__file__).resolve()

T_MARKERS = [
    "CD247",
    "LCK",
    "TRAT1",
    "BCL11B",
    "THEMIS",
    "CD3D",
    "CD3E",
    "CD3G",
    "TRAC",
    "TRBC1",
    "TRBC2",
]
CD3_MARKERS = ["CD3D", "CD3E", "CD3G"]
T_SIGNAL_MARKERS = ["CD247", "LCK", "TRAT1", "BCL11B", "THEMIS", *CD3_MARKERS]
TCR_TRANSCRIPT_MARKERS = ["TRAC", "TRBC1", "TRBC2"]
NK_MARKERS = [
    "KLRD1",
    "NCR1",
    "FCGR3A",
    "KLRC1",
    "KLRF1",
    "S1PR5",
    "FCER1G",
    "TYROBP",
]
TCR_DETAIL_COLUMNS = [
    "sample_id",
    "library_id",
    "donor_id",
    "original_cell_annotation",
    "assay_type",
    "tcr_availability",
    "TCRseq",
    "TRA_cdr3",
    "TRB_cdr3",
    "TRD_cdr3",
    "TRA_v",
    "TRB_v",
    "TRA_j",
    "TRB_j",
    "TRA_umis",
    "TRB_umis",
    "TRA_reads",
    "TRB_reads",
]
MISSING_TEXT = {"", "nan", "none", "na", "null", "<na>"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--skip-source-rescan", action="store_true")
    return parser.parse_args()


def clean_text(values: pd.Series | Iterable[Any]) -> pd.Series:
    result = pd.Series(values, dtype="string").fillna("").str.strip()
    return result.mask(result.str.lower().isin(MISSING_TEXT), "")


def nonempty(values: pd.Series | Iterable[Any]) -> np.ndarray:
    return clean_text(values).ne("").to_numpy(bool)


def positive_numeric(values: pd.Series | Iterable[Any]) -> pd.Series:
    result = pd.to_numeric(pd.Series(values), errors="coerce")
    # Legacy importers used zero when UMI/read support was unavailable.
    return result.mask(result <= 0)


def load_csr_npz(path: Path) -> sparse.csr_matrix:
    saved = np.load(path, allow_pickle=False)
    return sparse.csr_matrix(
        (saved["data"], saved["indices"], saved["indptr"]),
        shape=tuple(saved["shape"].astype(int)),
    )


def hit_count(matrix: sparse.csr_matrix, genes: Sequence[str], selected: Sequence[str]) -> np.ndarray:
    lookup = {gene: index for index, gene in enumerate(genes)}
    columns = [lookup[gene] for gene in selected if gene in lookup]
    if not columns:
        return np.zeros(matrix.shape[0], dtype=np.int16)
    return np.asarray((matrix[:, columns] > 0).sum(axis=1)).ravel().astype(np.int16)


def cramers_v(left: Sequence[Any], right: Sequence[Any]) -> float:
    table = pd.crosstab(pd.Series(left, dtype="string"), pd.Series(right, dtype="string"))
    if table.empty or min(table.shape) < 2:
        return float("nan")
    chi2 = stats.chi2_contingency(table, correction=False)[0]
    n = table.to_numpy().sum()
    phi2 = chi2 / n
    rows, columns = table.shape
    corrected = max(0.0, phi2 - ((columns - 1) * (rows - 1)) / max(1, n - 1))
    row_corrected = rows - ((rows - 1) ** 2) / max(1, n - 1)
    col_corrected = columns - ((columns - 1) ** 2) / max(1, n - 1)
    denominator = min(col_corrected - 1, row_corrected - 1)
    return float(math.sqrt(corrected / denominator)) if denominator > 0 else float("nan")


def table_html(frame: pd.DataFrame, max_rows: int = 40) -> str:
    shown = frame.head(max_rows).copy()
    for column in shown.select_dtypes(include=["float"]).columns:
        shown[column] = shown[column].map(lambda value: "" if pd.isna(value) else f"{value:.4g}")
    suffix = "" if frame.shape[0] <= max_rows else f"<p class='small'>Showing {max_rows} of {frame.shape[0]:,} rows.</p>"
    return shown.to_html(index=False, escape=True, classes="data") + suffix


def save_figure(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def load_boundary_context(config: dict[str, Any]) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, dict[str, np.ndarray], list[str], sparse.csr_matrix]:
    paths = stage_paths(config)
    required = [
        paths["staged_h5ad"],
        paths["marker_counts"],
        paths["boundary_partitions"],
        paths["boundary_umap_sample"],
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing V4.2 boundary artifacts: {missing}")

    backed = ad.read_h5ad(paths["staged_h5ad"], backed="r")
    try:
        obs = backed.obs.copy()
    finally:
        backed.file.close()
    partition = np.load(paths["boundary_partitions"], allow_pickle=False)
    boundary_indices = partition["boundary_indices"].astype(np.int64)
    names = partition["names"].astype(str).tolist()
    review = str(config["nk_boundary"]["review_run"])
    if review not in names:
        raise KeyError(f"Boundary review run absent: {review}")
    labels = partition["labels"][:, names.index(review)].astype(np.int32)
    sample_npz = np.load(paths["boundary_umap_sample"], allow_pickle=False)
    sample = {key: np.asarray(sample_npz[key]) for key in sample_npz.files}
    marker_table = pd.read_csv(paths["table"] / "marker_panel.csv")
    genes = marker_table.sort_values("marker_index")["gene"].astype(str).tolist()
    expected = marker_genes(config)
    if genes != expected:
        raise RuntimeError("Marker NPZ gene order differs from the marker-panel contract")
    marker_matrix = load_csr_npz(paths["marker_counts"])[boundary_indices]
    if marker_matrix.shape != (boundary_indices.size, len(genes)):
        raise RuntimeError("Boundary marker matrix shape mismatch")
    return obs, boundary_indices, labels, sample, genes, marker_matrix


def extract_boundary_tcr_details(
    config: dict[str, Any],
    obs: pd.DataFrame,
    boundary_indices: np.ndarray,
) -> pd.DataFrame:
    roles = development_roles(config)
    boundary_obs = obs.iloc[boundary_indices].copy()
    result = pd.DataFrame(index=np.arange(boundary_indices.size))
    result["input_cohort_id"] = boundary_obs["input_cohort_id"].astype(str).to_numpy()
    result["source_gse_id"] = boundary_obs["source_gse_id"].astype(str).to_numpy()
    result["integration_cell_id"] = boundary_obs.index.astype(str).to_numpy()
    result["input_row"] = boundary_obs["input_row"].to_numpy(np.int64)
    result["primary_nk_anchor"] = boundary_obs["primary_nk_anchor"].astype(bool).to_numpy()
    for column in ["source_original_cell_id", *TCR_DETAIL_COLUMNS]:
        result[column] = ""

    covered = np.zeros(result.shape[0], dtype=bool)
    for role in roles.itertuples(index=False):
        local = np.flatnonzero(result["input_cohort_id"].eq(role.cohort_id).to_numpy())
        if local.size == 0:
            continue
        rows = result.iloc[local]["input_row"].to_numpy(np.int64)
        if role.cohort_id == "current_atlas":
            source_ids = h5ad_obs_frame(role.path, ["original_cell_id"])
            target = pd.DataFrame(
                {
                    "source_gse_id": result.iloc[local]["source_gse_id"].to_numpy(),
                    "source_original_cell_id": source_ids.iloc[rows]["original_cell_id"].to_numpy(),
                },
                index=local,
            )
            target["source_gse_id"] = normalize_text(target["source_gse_id"])
            target["source_original_cell_id"] = normalize_text(target["source_original_cell_id"])
            target["key"] = target["source_gse_id"] + "||" + target["source_original_cell_id"]
            if target["key"].duplicated().any():
                raise RuntimeError("Current-atlas boundary source-plus-cell keys are not unique")
            wanted = set(target["key"])
            blocks: list[pd.DataFrame] = []
            usecols = ["source_gse_id", "original_cell_id", *TCR_DETAIL_COLUMNS]
            for item in config["current_atlas_recovery"]["metadata_sources"]:
                path = resolve(item["path"])
                header = pd.read_csv(path, nrows=0)
                absent = sorted(set(usecols) - set(header.columns))
                if absent:
                    raise KeyError(f"{path} lacks boundary audit fields: {absent}")
                for chunk in pd.read_csv(
                    path,
                    usecols=usecols,
                    dtype="string",
                    chunksize=400_000,
                    low_memory=False,
                ):
                    source = normalize_text(chunk["source_gse_id"])
                    cell = normalize_text(chunk["original_cell_id"])
                    key = source + "||" + cell
                    selected = key.isin(wanted)
                    if selected.any():
                        block = chunk.loc[selected, TCR_DETAIL_COLUMNS].copy()
                        block.index = pd.Index(key.loc[selected].to_numpy())
                        blocks.append(block)
            metadata = pd.concat(blocks, axis=0, copy=False)
            if metadata.index.duplicated().any():
                raise RuntimeError("Recovery metadata has duplicate boundary keys")
            missing = pd.Index(target["key"]).difference(metadata.index)
            if not missing.empty:
                raise RuntimeError(f"Recovery metadata missed {missing.size:,} boundary cells")
            aligned = metadata.reindex(target["key"].to_numpy())
            result.loc[local, "source_original_cell_id"] = target["source_original_cell_id"].to_numpy()
        else:
            backed = ad.read_h5ad(role.path, backed="r")
            try:
                source_obs = backed.obs
                aligned = pd.DataFrame(index=np.arange(rows.size))
                for column in TCR_DETAIL_COLUMNS:
                    aligned[column] = (
                        source_obs.iloc[rows][column].to_numpy()
                        if column in source_obs
                        else np.repeat("", rows.size)
                    )
                result.loc[local, "source_original_cell_id"] = source_obs.index[rows].astype(str)
            finally:
                backed.file.close()
        for column in TCR_DETAIL_COLUMNS:
            result.loc[local, column] = aligned[column].to_numpy()
        covered[local] = True
    if not covered.all():
        raise RuntimeError(f"TCR detail extraction missed {(~covered).sum():,} boundary cells")
    return result


def add_chain_and_expression_features(
    details: pd.DataFrame,
    labels: np.ndarray,
    marker_matrix: sparse.csr_matrix,
    genes: Sequence[str],
) -> pd.DataFrame:
    result = details.copy()
    result["boundary_cluster"] = labels
    for chain in ("TRA", "TRB", "TRD"):
        result[f"productive_{chain.lower()}"] = nonempty(result[f"{chain}_cdr3"])
    result["productive_ab"] = result["productive_tra"] | result["productive_trb"]
    result["paired_ab"] = result["productive_tra"] & result["productive_trb"]
    for chain in ("TRA", "TRB"):
        result[f"{chain.lower()}_umi"] = positive_numeric(result[f"{chain}_umis"]).to_numpy()
        result[f"{chain.lower()}_reads"] = positive_numeric(result[f"{chain}_reads"]).to_numpy()

    result["any_ab_umi_known"] = result[["tra_umi", "trb_umi"]].notna().any(axis=1)
    result["both_ab_umi_known"] = result[["tra_umi", "trb_umi"]].notna().all(axis=1)
    result["any_ab_umi_ge2"] = result[["tra_umi", "trb_umi"]].ge(2).any(axis=1)
    result["both_ab_umi_ge2"] = result[["tra_umi", "trb_umi"]].ge(2).all(axis=1)
    result["one_umi_only"] = (
        result["productive_ab"]
        & result["any_ab_umi_known"]
        & result[["tra_umi", "trb_umi"]].fillna(0).max(axis=1).eq(1)
    )
    result["both_chains_one_umi"] = (
        result["paired_ab"] & result["tra_umi"].eq(1) & result["trb_umi"].eq(1)
    )
    result["any_ab_reads_ge10"] = result[["tra_reads", "trb_reads"]].ge(10).any(axis=1)

    result["t_marker_hits"] = hit_count(marker_matrix, genes, T_MARKERS)
    result["t_signal_hits"] = hit_count(marker_matrix, genes, T_SIGNAL_MARKERS)
    result["cd3_hits"] = hit_count(marker_matrix, genes, CD3_MARKERS)
    result["tcr_transcript_hits"] = hit_count(marker_matrix, genes, TCR_TRANSCRIPT_MARKERS)
    result["nk_marker_hits"] = hit_count(marker_matrix, genes, NK_MARKERS)
    result["multi_gene_t_support"] = result["t_signal_hits"].ge(2) | result["cd3_hits"].ge(2)
    result["multi_gene_nk_support"] = result["nk_marker_hits"].ge(2)
    result["expression_class"] = "low_information"
    result.loc[result["multi_gene_t_support"] & ~result["multi_gene_nk_support"], "expression_class"] = "T_like"
    result.loc[~result["multi_gene_t_support"] & result["multi_gene_nk_support"], "expression_class"] = "NK_like"
    result.loc[result["multi_gene_t_support"] & result["multi_gene_nk_support"], "expression_class"] = "mixed_T_NK"

    result["umi_support_tier"] = "no_productive_TRA_TRB"
    ab = result["productive_ab"]
    result.loc[ab & ~result["any_ab_umi_known"], "umi_support_tier"] = "UMI_unavailable"
    result.loc[ab & result["one_umi_only"], "umi_support_tier"] = "one_UMI_only"
    result.loc[ab & result["any_ab_umi_ge2"] & ~result["paired_ab"], "umi_support_tier"] = "single_chain_UMI_ge2"
    result.loc[ab & result["paired_ab"] & result["any_ab_umi_ge2"], "umi_support_tier"] = "paired_any_UMI_ge2"
    return result


def choose_group_key(obs: pd.DataFrame, candidates: Sequence[str]) -> tuple[str, pd.Series]:
    best_name = ""
    best = pd.Series(np.repeat("", obs.shape[0]), index=obs.index, dtype="string")
    best_unique = 1
    for column in candidates:
        if column not in obs:
            continue
        values = clean_text(obs[column])
        unique = int(values[values.ne("")].nunique())
        if 1 < unique < min(10_000, max(2, obs.shape[0] // 2)) and unique > best_unique:
            best_name, best, best_unique = column, values, unique
    return best_name, best


def pair_reuse_metrics(pair: pd.Series, groups: pd.Series, prefix: str) -> dict[str, Any]:
    valid = pair.ne("") & groups.ne("")
    if valid.sum() == 0 or groups.loc[valid].nunique() < 2:
        return {
            f"n_{prefix}_groups": int(groups[groups.ne("")].nunique()),
            f"n_pairs_across_{prefix}s": 0,
            f"max_{prefix}s_per_pair": 0,
            f"fraction_ab_cells_pair_across_{prefix}s": 0.0,
        }
    work = pd.DataFrame({"pair": pair.loc[valid].to_numpy(), "group": groups.loc[valid].to_numpy()})
    spread = work.groupby("pair", observed=True)["group"].nunique()
    mapped = work["pair"].map(spread)
    return {
        f"n_{prefix}_groups": int(work["group"].nunique()),
        f"n_pairs_across_{prefix}s": int(spread.ge(2).sum()),
        f"max_{prefix}s_per_pair": int(spread.max()),
        f"fraction_ab_cells_pair_across_{prefix}s": float(mapped.ge(2).mean()),
    }


def source_h5ad_audit(
    sources: Sequence[str],
    mapping_modes: dict[str, str],
    role_paths: dict[str, str],
) -> pd.DataFrame:
    registry = pd.read_csv(PROJECT_ROOT / "configs/datasets/datasets.csv", dtype="string")
    registry_paths = dict(zip(registry["dataset_id"], registry["processed_h5ad_path"], strict=False))
    records: list[dict[str, Any]] = []
    annotation_candidates = [
        "original_cell_annotation",
        "cell_type",
        "celltype",
        "major_cell_type",
        "author_cell_type",
        "predicted.celltype.l2",
        "scvi_prediction",
    ]
    sample_candidates = [
        "library_id",
        "sample_id",
        "sampleID",
        "sample",
        "_library_id_from_obs",
        "matrix_file",
        "library_index",
        "orig.ident",
    ]
    donor_candidates = ["donor_id", "donor_patient", "patient_assignment", "donor", "sampleID"]
    for source in sorted(set(sources)):
        raw_path = role_paths.get(source, registry_paths.get(source, ""))
        path = resolve(raw_path) if raw_path else Path("")
        base: dict[str, Any] = {
            "source_gse_id": source,
            "source_h5ad": str(path) if raw_path else "",
            "reported_mapping_mode": mapping_modes.get(source, "not_recorded"),
        }
        if not raw_path or not path.is_file():
            records.append({**base, "audit_status": "source_h5ad_unavailable"})
            continue
        try:
            backed = ad.read_h5ad(path, backed="r")
            try:
                obs = backed.obs.copy()
            finally:
                backed.file.close()
        except Exception as error:  # pragma: no cover - dataset-specific corruption path
            records.append({**base, "audit_status": f"source_obs_read_failed:{type(error).__name__}"})
            continue
        blank = pd.Series(np.repeat("", obs.shape[0]), index=obs.index, dtype="string")
        tra = clean_text(obs["TRA_cdr3"]) if "TRA_cdr3" in obs else blank
        trb = clean_text(obs["TRB_cdr3"]) if "TRB_cdr3" in obs else blank
        has_tra, has_trb = tra.ne(""), trb.ne("")
        has_ab = has_tra | has_trb
        pair = (tra.astype(str) + "||" + trb.astype(str)).where(has_ab, "")
        pair_counts = pair[has_ab].value_counts()
        sample_name, sample_key = choose_group_key(obs, sample_candidates)
        donor_name, donor_key = choose_group_key(obs, donor_candidates)
        record: dict[str, Any] = {
            **base,
            "n_source_cells": int(obs.shape[0]),
            "n_source_productive_ab": int(has_ab.sum()),
            "fraction_source_productive_ab": float(has_ab.mean()),
            "n_source_paired_ab": int((has_tra & has_trb).sum()),
            "n_unique_ab_pairs": int(pair_counts.size),
            "median_pair_multiplicity": float(pair_counts.median()) if not pair_counts.empty else 0.0,
            "p95_pair_multiplicity": float(pair_counts.quantile(0.95)) if not pair_counts.empty else 0.0,
            "max_pair_multiplicity": int(pair_counts.max()) if not pair_counts.empty else 0,
            "fraction_ab_cells_pair_multiplicity_ge10": float(pair[has_ab].map(pair_counts).ge(10).mean()) if has_ab.any() else 0.0,
            "sample_key_column": sample_name,
            "donor_key_column": donor_name,
        }
        record.update(pair_reuse_metrics(pair, sample_key, "sample"))
        record.update(pair_reuse_metrics(pair, donor_key, "donor"))

        lineage = pd.Series(np.repeat("", obs.shape[0]), index=obs.index, dtype="string")
        used_annotations: list[str] = []
        for column in annotation_candidates:
            if column in obs:
                values = clean_text(obs[column])
                if values.ne("").any():
                    lineage = lineage + "|" + values
                    used_annotations.append(column)
        lower = lineage.str.lower()
        nk = lower.str.contains(r"(?:^|\W)nk(?:\W|$)|natural killer", regex=True, na=False)
        clearly_non_t = lower.str.contains(
            r"b cell|b-cell|myeloid|monocyte|macroph|dendritic|(?:^|\W)dc(?:\W|$)|mast|neutro",
            regex=True,
            na=False,
        )
        record.update(
            {
                "annotation_columns": ";".join(used_annotations),
                "n_author_nk_cells": int(nk.sum()),
                "n_author_nk_with_ab": int((nk & has_ab).sum()),
                "fraction_author_nk_with_ab": float(has_ab[nk].mean()) if nk.any() else float("nan"),
                "n_author_non_t_cells": int(clearly_non_t.sum()),
                "n_author_non_t_with_ab": int((clearly_non_t & has_ab).sum()),
                "fraction_author_non_t_with_ab": float(has_ab[clearly_non_t].mean()) if clearly_non_t.any() else float("nan"),
            }
        )
        mode = record["reported_mapping_mode"]
        documented_unsafe = "barcode_only" in str(mode)
        donor_reuse = (
            record["max_donors_per_pair"] >= 3
            and record["n_pairs_across_donors"] >= 10
            and record["fraction_ab_cells_pair_across_donors"] >= 0.02
        )
        sample_reuse = (
            record["max_samples_per_pair"] >= 5
            and record["n_pairs_across_samples"] >= 50
            and record["fraction_ab_cells_pair_across_samples"] >= 0.10
        )
        non_t_conflict = (
            record["n_author_non_t_cells"] >= 100
            and record["n_author_non_t_with_ab"] >= 100
            and record["fraction_author_non_t_with_ab"] >= 0.20
        )
        if documented_unsafe:
            status = "documented_unsafe_barcode_join"
        elif donor_reuse:
            status = "highly_suspicious_cross_donor_reuse"
        elif sample_reuse:
            status = "highly_suspicious_cross_sample_reuse"
        elif non_t_conflict:
            status = "strong_author_lineage_conflict"
        else:
            status = "no_join_failure_detected"
        record["audit_status"] = status
        record["join_red_flag"] = status != "no_join_failure_detected"
        records.append(record)
    return pd.DataFrame(records)


def boundary_pair_reuse(details: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for source, group in details.groupby("source_gse_id", observed=True):
        tra = clean_text(group["TRA_cdr3"])
        trb = clean_text(group["TRB_cdr3"])
        has_ab = tra.ne("") | trb.ne("")
        pair = (tra.astype(str) + "||" + trb.astype(str)).where(has_ab, "")
        sample = clean_text(group["library_id"])
        if sample[sample.ne("")].nunique() < 2:
            sample = clean_text(group["sample_id"])
        donor = clean_text(group["donor_id"])
        pair_counts = pair[has_ab].value_counts()
        record = {
            "source_gse_id": source,
            "n_boundary_cells": int(group.shape[0]),
            "n_boundary_productive_ab": int(has_ab.sum()),
            "fraction_boundary_productive_ab": float(has_ab.mean()),
            "n_boundary_unique_pairs": int(pair_counts.size),
            "boundary_max_pair_multiplicity": int(pair_counts.max()) if not pair_counts.empty else 0,
        }
        record.update(pair_reuse_metrics(pair, sample, "boundary_sample"))
        record.update(pair_reuse_metrics(pair, donor, "boundary_donor"))
        records.append(record)
    return pd.DataFrame(records)


def summarize_boundary(details: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    source = (
        details.groupby(["source_gse_id", "join_audit_status"], observed=True)
        .agg(
            n_boundary_cells=("productive_ab", "size"),
            n_productive_ab=("productive_ab", "sum"),
            n_paired_ab=("paired_ab", "sum"),
            n_umi_observed=("any_ab_umi_known", "sum"),
            n_one_umi_only=("one_umi_only", "sum"),
            n_any_umi_ge2=("any_ab_umi_ge2", "sum"),
            n_multigene_t_support=("multi_gene_t_support", "sum"),
            n_primary_nk_anchor=("primary_nk_anchor", "sum"),
        )
        .reset_index()
    )
    source["fraction_productive_ab"] = source["n_productive_ab"] / source["n_boundary_cells"]
    source["fraction_ab_one_umi_only"] = source["n_one_umi_only"] / source["n_productive_ab"].replace(0, np.nan)
    source["fraction_ab_multigene_t_support"] = (
        details.loc[details["productive_ab"]]
        .groupby("source_gse_id", observed=True)["multi_gene_t_support"]
        .mean()
        .reindex(source["source_gse_id"])
        .to_numpy()
    )
    source = source.sort_values("n_productive_ab", ascending=False).reset_index(drop=True)

    cluster = (
        details.groupby("boundary_cluster", observed=True)
        .agg(
            n_cells=("productive_ab", "size"),
            n_productive_ab=("productive_ab", "sum"),
            n_join_red_flag_ab=("ab_from_join_red_flag_source", "sum"),
            mean_t_marker_hits=("t_marker_hits", "mean"),
            mean_cd3_hits=("cd3_hits", "mean"),
            mean_nk_marker_hits=("nk_marker_hits", "mean"),
            fraction_multigene_t_support=("multi_gene_t_support", "mean"),
            fraction_multigene_nk_support=("multi_gene_nk_support", "mean"),
        )
        .reset_index()
    )
    cluster["fraction_productive_ab"] = cluster["n_productive_ab"] / cluster["n_cells"]
    cluster["fraction_join_red_flag_ab"] = cluster["n_join_red_flag_ab"] / cluster["n_cells"]

    umi = (
        details.groupby(["source_gse_id", "umi_support_tier"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    return source, cluster, umi


def statistical_audit(details: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    strata = {
        "all_boundary": np.ones(details.shape[0], dtype=bool),
        "join_red_flag_sources": details["join_red_flag"].to_numpy(bool),
        "sources_without_join_red_flag": ~details["join_red_flag"].to_numpy(bool),
    }
    for name, mask in strata.items():
        group = details.loc[mask]
        if group.empty or group["productive_ab"].nunique() < 2:
            continue
        table = pd.crosstab(group["productive_ab"], group["multi_gene_t_support"]).reindex(
            index=[False, True], columns=[False, True], fill_value=0
        )
        odds, fisher_p = stats.fisher_exact(table.to_numpy())
        positive = group.loc[group["productive_ab"], "t_marker_hits"]
        negative = group.loc[~group["productive_ab"], "t_marker_hits"]
        mann = stats.mannwhitneyu(positive, negative, alternative="two-sided")
        auc = roc_auc_score(group["productive_ab"], group["t_marker_hits"])
        records.append(
            {
                "stratum": name,
                "n_cells": int(group.shape[0]),
                "n_productive_ab": int(group["productive_ab"].sum()),
                "fisher_odds_ratio_multigene_t_support": float(odds),
                "fisher_p": float(fisher_p),
                "median_t_hits_ab": float(positive.median()),
                "median_t_hits_no_ab": float(negative.median()),
                "mannwhitney_u": float(mann.statistic),
                "mannwhitney_p": float(mann.pvalue),
                "auc_t_marker_hits_for_ab": float(auc),
            }
        )
    return pd.DataFrame(records)


def cluster_association_audit(details: pd.DataFrame) -> pd.DataFrame:
    records = []
    for variable in ["productive_ab", "expression_class", "source_gse_id", "join_audit_status"]:
        records.append(
            {
                "variable": variable,
                "normalized_mutual_information_with_boundary_cluster": float(
                    normalized_mutual_info_score(details["boundary_cluster"].astype(str), details[variable].astype(str))
                ),
                "cramers_v_with_boundary_cluster": cramers_v(details["boundary_cluster"], details[variable]),
            }
        )
    return pd.DataFrame(records)


def plot_umap_provenance(details: pd.DataFrame, sample: dict[str, np.ndarray], path: Path) -> None:
    local = sample["sample_local_indices"].astype(np.int64)
    frame = details.iloc[local]
    xy = sample["X_umap"]
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 10.0), sharex=True, sharey=True)
    panels = [
        ("productive_ab", {False: "#d6d8da", True: "#176f9b"}, "Raw productive TRA or TRB"),
        (
            "ab_provenance",
            {
                "No TRA/TRB": "#d6d8da",
                "AB from red-flag source": "#c84b31",
                "AB without detected join red flag": "#207a62",
            },
            "TRA/TRB calls stratified by join audit",
        ),
        (
            "umi_support_tier",
            {
                "no_productive_TRA_TRB": "#d6d8da",
                "UMI_unavailable": "#8b78a8",
                "one_UMI_only": "#d99a2b",
                "single_chain_UMI_ge2": "#4f8fbd",
                "paired_any_UMI_ge2": "#217a62",
            },
            "Per-chain UMI support",
        ),
        (
            "expression_class",
            {
                "low_information": "#b8b8b8",
                "NK_like": "#d2673f",
                "mixed_T_NK": "#8d6aae",
                "T_like": "#2277a5",
            },
            "Independent transcriptomic evidence",
        ),
    ]
    for ax, (column, palette, title) in zip(axes.ravel(), panels, strict=True):
        values = frame[column].to_numpy()
        ax.scatter(xy[:, 0], xy[:, 1], s=0.25, color="#eceeef", linewidths=0, rasterized=True)
        for label, color in palette.items():
            selected = values == label
            if selected.any():
                ax.scatter(
                    xy[selected, 0],
                    xy[selected, 1],
                    s=0.65,
                    color=color,
                    label=f"{label} ({selected.sum():,})",
                    linewidths=0,
                    rasterized=True,
                )
        ax.set_title(title, fontsize=11)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.legend(loc="best", fontsize=7, frameon=False, markerscale=5)
    fig.suptitle("V4.2 boundary: productive TRA/TRB provenance is not a cell-type label", fontsize=14)
    save_figure(fig, path)


def plot_source_forensics(source: pd.DataFrame, path: Path) -> None:
    shown = source.nlargest(18, "n_boundary_productive_ab").sort_values("n_boundary_productive_ab")
    colors = np.where(shown["join_red_flag"].fillna(False), "#c84b31", "#267a65")
    fig, axes = plt.subplots(1, 2, figsize=(14.0, 7.2))
    axes[0].barh(shown["source_gse_id"], shown["n_boundary_productive_ab"], color=colors)
    axes[0].set_xlabel("Boundary cells with productive TRA or TRB")
    axes[0].set_title("Calls contributed by each source")
    axes[0].grid(axis="x", color="#e5e5e5", linewidth=0.7)
    sample_fraction = shown["fraction_ab_cells_pair_across_samples"].fillna(0)
    donor_fraction = shown["fraction_ab_cells_pair_across_donors"].fillna(0)
    y = np.arange(shown.shape[0])
    axes[1].barh(y - 0.18, sample_fraction, height=0.34, color="#356fa3", label="Across sample keys")
    axes[1].barh(y + 0.18, donor_fraction, height=0.34, color="#c87831", label="Across donor keys")
    axes[1].set_yticks(y, shown["source_gse_id"])
    axes[1].set_xlim(0, 1)
    axes[1].set_xlabel("Fraction of AB cells whose exact TRA/TRB pair recurs")
    axes[1].set_title("Cross-group exact-pair reuse")
    axes[1].legend(frameon=False)
    axes[1].grid(axis="x", color="#e5e5e5", linewidth=0.7)
    fig.suptitle("Source-level evidence for unsafe or suspicious TCR joins", fontsize=14)
    save_figure(fig, path)


def plot_umi_support(umi: pd.DataFrame, source: pd.DataFrame, path: Path) -> None:
    top = source.nlargest(16, "n_productive_ab")["source_gse_id"].tolist()
    selected = umi.loc[umi["source_gse_id"].isin(top)].copy()
    matrix = selected.pivot_table(
        index="source_gse_id", columns="umi_support_tier", values="n_cells", aggfunc="sum", fill_value=0
    ).reindex(top)
    order = [
        "no_productive_TRA_TRB",
        "UMI_unavailable",
        "one_UMI_only",
        "single_chain_UMI_ge2",
        "paired_any_UMI_ge2",
    ]
    matrix = matrix.reindex(columns=[column for column in order if column in matrix], fill_value=0)
    fractions = matrix.div(matrix.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    palette = ["#d6d8da", "#8b78a8", "#d99a2b", "#4f8fbd", "#217a62"]
    fig, ax = plt.subplots(figsize=(12.5, 7.0))
    left = np.zeros(fractions.shape[0])
    for index, column in enumerate(fractions.columns):
        values = fractions[column].to_numpy()
        ax.barh(fractions.index, values, left=left, color=palette[index], label=column)
        left += values
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of boundary cells")
    ax.set_title("UMI support is frequently unavailable; one-UMI-only calls are separable")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.22), ncol=3, frameon=False)
    ax.grid(axis="x", color="#e5e5e5", linewidth=0.7)
    save_figure(fig, path)


def plot_cluster_diagnostics(cluster: pd.DataFrame, path: Path) -> None:
    shown = cluster.sort_values("boundary_cluster").copy()
    matrix = shown.set_index("boundary_cluster")[[
        "fraction_productive_ab",
        "fraction_join_red_flag_ab",
        "fraction_multigene_t_support",
        "fraction_multigene_nk_support",
    ]]
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2), gridspec_kw={"width_ratios": [1.2, 1]})
    sns.heatmap(matrix, cmap="vlag", vmin=0, vmax=1, annot=True, fmt=".2f", ax=axes[0], cbar_kws={"label": "Cell fraction"})
    axes[0].set_title("Chain provenance and expression by boundary cluster")
    axes[0].set_ylabel("Boundary cluster")
    x = np.arange(shown.shape[0])
    axes[1].plot(x, shown["mean_t_marker_hits"], marker="o", label="T-lineage hits", color="#2774ae")
    axes[1].plot(x, shown["mean_nk_marker_hits"], marker="o", label="NK-marker hits", color="#c65c35")
    axes[1].plot(x, shown["mean_cd3_hits"], marker="o", label="CD3 hits", color="#237b63")
    axes[1].set_xticks(x, shown["boundary_cluster"])
    axes[1].set_xlabel("Boundary cluster")
    axes[1].set_ylabel("Mean number of detected markers")
    axes[1].set_title("The clustering captures transcriptomic gradients")
    axes[1].legend(frameon=False)
    axes[1].grid(axis="y", color="#e5e5e5", linewidth=0.7)
    save_figure(fig, path)


def build_report(
    details: pd.DataFrame,
    source: pd.DataFrame,
    cluster: pd.DataFrame,
    stats_table: pd.DataFrame,
    association: pd.DataFrame,
    hvg_audit: pd.DataFrame,
    summary: dict[str, Any],
) -> None:
    figure_names = [
        "boundary_umap_tcr_provenance.png",
        "source_join_forensics.png",
        "umi_support_by_source.png",
        "cluster_marker_tcr_diagnostics.png",
    ]
    figure_static = STATIC_DIR / "figures"
    figure_static.mkdir(parents=True, exist_ok=True)
    for name in figure_names:
        shutil.copy2(FIGURE_DIR / name, figure_static / name)

    css = """
    :root{--ink:#18232b;--muted:#59666f;--line:#d8dee2;--blue:#176f9b;--green:#24765f;--red:#b94732;--amber:#a96c16;--paper:#fff}
    *{box-sizing:border-box} body{margin:0;background:#eef1f2;color:var(--ink);font:15px/1.5 Arial,sans-serif}
    main{max-width:1180px;margin:0 auto;background:var(--paper);padding:34px 48px 60px;box-shadow:0 3px 20px #0001}
    h1{font-size:30px;margin:0 0 6px;letter-spacing:0} h2{font-size:21px;border-top:2px solid var(--ink);padding-top:15px;margin-top:34px} h3{font-size:17px;margin-top:24px}
    .subtitle{color:var(--muted);font-size:16px}.grid{display:grid;grid-template-columns:repeat(4,1fr);gap:10px;margin:20px 0}.card{border:1px solid var(--line);border-top:4px solid var(--blue);padding:13px;background:#fafbfb}.card.red{border-top-color:var(--red)}.card.green{border-top-color:var(--green)}.card.amber{border-top-color:var(--amber)}.big{font-size:23px;font-weight:700}.small{font-size:12px;color:var(--muted)}
    .finding{padding:14px 16px;border-left:5px solid var(--blue);background:#f1f6f8;margin:14px 0}.finding.red{border-color:var(--red);background:#fbf1ef}.finding.green{border-color:var(--green);background:#eef7f3}.finding.amber{border-color:var(--amber);background:#fbf6ec}
    figure{margin:22px 0} img{display:block;max-width:100%;height:auto} figcaption{font-size:13px;color:var(--muted);margin-top:7px}.table-wrap{margin:15px 0}table.data{border-collapse:collapse;table-layout:fixed;width:100%;font-size:11px}table.data th,table.data td{border:1px solid var(--line);padding:5px 6px;vertical-align:top;overflow-wrap:anywhere;word-break:break-word}table.data th{background:#eef2f4;text-align:left}code{background:#eef1f2;padding:1px 4px}ol li{margin:8px 0}
    @media(max-width:800px){main{padding:24px 18px}.grid{grid-template-columns:1fr 1fr}}
    @media print{@page{size:A4 landscape;margin:10mm}body{background:#fff;font-size:10pt}main{max-width:none;box-shadow:none;padding:0}.grid{grid-template-columns:repeat(4,1fr)}h2{break-after:avoid}figure,.finding,.card{break-inside:avoid}table.data{font-size:7pt}a{color:inherit;text-decoration:none}}
    """
    source_display = source[[
        "source_gse_id",
        "n_boundary_cells",
        "n_boundary_productive_ab",
        "fraction_boundary_productive_ab",
        "audit_status",
        "reported_mapping_mode",
        "sample_key_column",
        "fraction_ab_cells_pair_across_samples",
        "donor_key_column",
        "fraction_ab_cells_pair_across_donors",
        "fraction_author_nk_with_ab",
        "fraction_author_non_t_with_ab",
    ]].sort_values("n_boundary_productive_ab", ascending=False)
    source_overview = source_display[[
        "source_gse_id",
        "n_boundary_cells",
        "n_boundary_productive_ab",
        "fraction_boundary_productive_ab",
        "audit_status",
        "reported_mapping_mode",
    ]].rename(columns={
        "source_gse_id": "source",
        "n_boundary_cells": "boundary_cells",
        "n_boundary_productive_ab": "productive_AB",
        "fraction_boundary_productive_ab": "AB_fraction",
        "audit_status": "join_audit",
        "reported_mapping_mode": "recorded_mapping",
    })
    source_evidence = source_display[[
        "source_gse_id",
        "sample_key_column",
        "fraction_ab_cells_pair_across_samples",
        "donor_key_column",
        "fraction_ab_cells_pair_across_donors",
        "fraction_author_nk_with_ab",
        "fraction_author_non_t_with_ab",
    ]].rename(columns={
        "source_gse_id": "source",
        "sample_key_column": "sample_key",
        "fraction_ab_cells_pair_across_samples": "AB_exact_pair_cross_sample",
        "donor_key_column": "donor_key",
        "fraction_ab_cells_pair_across_donors": "AB_exact_pair_cross_donor",
        "fraction_author_nk_with_ab": "author_NK_with_AB",
        "fraction_author_non_t_with_ab": "author_nonT_with_AB",
    })
    html_text = f"""<!doctype html><html><head><meta charset='utf-8'><title>V4.2 productive TRA/TRB conflict audit</title><style>{css}</style></head><body><main>
    <h1>Why is productive TRA/TRB so frequent in the V4.2 NK boundary?</h1>
    <p class='subtitle'>A four-hypothesis forensic audit of {summary['n_boundary_cells']:,} cells. This report is diagnostic and does not redefine NK truth or mutate source data.</p>
    <div class='grid'>
      <div class='card'><div class='big'>{summary['fraction_productive_ab']:.1%}</div><div>raw productive TRA/TRB</div></div>
      <div class='card red'><div class='big'>{summary['fraction_ab_from_red_flag_sources']:.1%}</div><div>of AB calls from join-red-flag sources</div></div>
      <div class='card amber'><div class='big'>{summary['fraction_one_umi_only_among_umi_observed_ab']:.1%}</div><div>one-UMI-only among AB calls with observed UMI support</div></div>
      <div class='card green'><div class='big'>{summary['n_high_confidence_residual_ab']:,}</div><div>conservative residual AB-compatible cells</div></div>
    </div>
    <div class='finding red'><b>Main conclusion.</b> The high TRA/TRB rate cannot be interpreted as an inconvenient biological fact. The dominant detected problem is upstream TCR assignment provenance: legacy barcode-only joins are documented, exact TRA/TRB pairs recur across multiple sample or donor keys in several sources, and some source objects assign TRA/TRB broadly to author-labeled B/myeloid/NK cells. UMI depth cannot repair a receptor copied to the wrong cell.</div>
    <div class='finding amber'><b>Ambient reads are a minority explanation in the auditable subset.</b> Only {summary['n_one_umi_only_ab']:,} productive-AB cells are supported exclusively by one observed UMI. UMI support is unavailable for {summary['n_ab_umi_unavailable']:,} AB cells, so ambient contamination cannot be ruled out there; it also cannot explain cross-donor reuse of exact paired receptors.</div>
    <div class='finding'><b>The UMAP did not create the TCR conflict.</b> UMAP is a visualization of the 30-dimensional scVI latent representation. The 4,000-gene model forced core T/NK lineage genes and deliberately excluded TCR V/J/D genes. Clusters retain distinct CD3/T/NK expression gradients, while suspect TRA/TRB assignments span all clusters. Re-clustering alone would not correct an unsafe metadata join.</div>

    <h2>1. Scope and definitions</h2>
    <p>The boundary is the frozen union of refined parent clusters 9 and 18 from the T/NK-restricted V4.2 integration. A productive chain call means a nonempty productive-filtered CDR3 field. Zero UMI/read values are treated as unavailable because legacy import code used zero when support was absent. A one-UMI-only call requires an observed positive UMI value and no chain with two or more UMIs.</p>
    <p>A source receives a join red flag only from auditable evidence: a documented barcode-only mapping mode, extensive exact paired-receptor reuse across sample/donor keys, or strong conflict with source author lineages. This is a quarantine signal, not an estimate of how many individual calls are wrong.</p>
    <figure><img src='figures/boundary_umap_tcr_provenance.png'><figcaption>The same cells are colored four ways. Red marks productive alpha/beta calls inherited from a source with a join red flag; expression classes are descriptive multi-gene evidence, not hard truth labels.</figcaption></figure>

    <h2>2. Hypothesis 1: ambient TCR reads</h2>
    <p>Among {summary['n_ab_umi_observed']:,} productive-AB cells with at least one observed positive chain UMI, {summary['n_one_umi_only_ab']:,} ({summary['fraction_one_umi_only_among_umi_observed_ab']:.2%}) have no chain above one UMI. In contrast, {summary['n_ab_any_umi_ge2']:,} have at least one chain with two or more UMIs, and {summary['n_ab_reads_ge10']:,} have at least ten reads for a chain. This argues against ambient one-UMI calls as the main explanation in sources with quantitative support.</p>
    <p><b>Important limitation:</b> {summary['n_ab_umi_unavailable']:,} productive-AB cells have no usable UMI field. They remain unresolved for ambient contamination. Moreover, high UMI/read support validates the receptor in its originating VDJ library, but not the RNA cell to which it was joined.</p>
    <figure><img src='figures/umi_support_by_source.png'><figcaption>UMI support tiers by source. Purple indicates missing quantitative support, not one UMI.</figcaption></figure>

    <h2>3. Hypothesis 2: non-ideal HVGs, UMAP, or clustering</h2>
    <p>The refined scVI model used {summary['n_selected_hvgs']:,} genes. It forced {summary['n_forced_lineage_selected']:,}/{summary['n_forced_lineage_expected']:,} specified T/NK lineage genes into the model. TCR V/J/D genes were intentionally excluded so clonotype usage would not drive atlas geometry; constant-region and signaling genes such as <code>CD3D/E/G</code>, <code>TRAC</code>, <code>TRBC1/2</code>, <code>LCK</code>, <code>TRAT1</code>, <code>BCL11B</code>, <code>KLRD1</code>, <code>FCER1G</code>, and <code>TYROBP</code> were retained.</p>
    <p>Boundary cluster association is quantified below. The cluster-expression association should be interpreted alongside the much broader source and raw-chain associations. A better boundary-specific integration may sharpen biology after TCR repair, but it cannot make an unsafe receptor join valid.</p>
    <figure><img src='figures/cluster_marker_tcr_diagnostics.png'><figcaption>Clusters show transcriptomic gradients even though raw productive TRA/TRB calls are pervasive. The red-flag component demonstrates why raw chain overlay is not valid NK/T truth.</figcaption></figure>
    <div class='table-wrap'>{table_html(association)}</div>
    <div class='table-wrap'>{table_html(hvg_audit)}</div>

    <h2>4. Hypothesis 3: TCRs joined to the wrong cells</h2>
    <p>The project SOP states that historical barcode-only or barcode-prefix joins are unsafe because 10x barcodes are reused across samples. The legacy mapping report explicitly records such modes for affected sources. Independently, exact TRA/TRB pairs recur across multiple sample or donor keys, sometimes in a large fraction of assigned cells. Exact paired alpha/beta receptors shared at this scale across unrelated groups are not a plausible ambient-RNA pattern.</p>
    <figure><img src='figures/source_join_forensics.png'><figcaption>Red bars are sources with a documented or strong join red flag. Cross-group pair reuse is calculated from the source H5AD where suitable grouping keys are present.</figcaption></figure>
    <h3>Source overview</h3><div class='table-wrap'>{table_html(source_overview, 60)}</div>
    <h3>Exact-pair reuse and lineage conflicts</h3><div class='table-wrap'>{table_html(source_evidence, 60)}</div>

    <h2>5. Hypothesis 4: these are truly alpha-beta T cells</h2>
    <p>Some cells are likely genuine alpha-beta T cells: the boundary was deliberately broad, cytotoxic T and NK transcriptomes overlap, and productive paired receptors with multi-gene T-lineage expression are biologically coherent. However, truth can only be assessed after quarantining red-flag sources.</p>
    <p>The conservative residual count is {summary['n_high_confidence_residual_ab']:,}: paired TRA/TRB, no detected source join red flag, at least one chain with two or more observed UMIs, and multi-gene T-lineage support. A further {summary['n_compatible_but_umi_missing_residual_ab']:,} cells are compatible with alpha-beta T identity but lack quantitative UMI support. These are descriptive audit strata, not replacement ground-truth labels.</p>
    <div class='table-wrap'>{table_html(stats_table)}</div>

    <h2>6. Decision and corrective action</h2>
    <ol>
      <li>Quarantine productive TRA/TRB fields from every documented or strongly suspicious source for model truth and NK exclusion. Do not delete the original fields; add provenance-qualified replacements.</li>
      <li>Rebuild affected datasets from raw productive VDJ contigs using the canonical <code>sample_id + barcode_core</code> key, retaining UMI/read support and stopping when sample identity is ambiguous.</li>
      <li>Validate each rebuilt source by key uniqueness, before/after coverage, exact receptor reuse across donors, and author-lineage conflicts. A high UMI count alone is not a join validation.</li>
      <li>Propagate only validated receptor fields into harmonized metadata, then regenerate the T/NK subset, HVGs, scVI latent space, and boundary clustering.</li>
      <li>Keep TCR V/J/D genes out of integration HVGs unless performing a separate sensitivity analysis; use them as classifier features or overlays after join repair, not as drivers of cell-state geometry.</li>
    </ol>
    <p class='small'>Generated by <code>{html.escape(str(SCRIPT_PATH.relative_to(PROJECT_ROOT)))}</code>. Source H5AD files were read only.</p>
    </main></body></html>"""
    STATIC_DIR.mkdir(parents=True, exist_ok=True)
    (STATIC_DIR / "index.html").write_text(html_text, encoding="utf-8")

    markdown = f"""# V4.2 Productive TRA/TRB Conflict Audit

## Main result

- Boundary cells: {summary['n_boundary_cells']:,}
- Raw productive TRA/TRB: {summary['n_productive_ab']:,} ({summary['fraction_productive_ab']:.2%})
- Productive-AB calls from join-red-flag sources: {summary['n_ab_from_red_flag_sources']:,} ({summary['fraction_ab_from_red_flag_sources']:.2%} of AB calls)
- One-UMI-only among AB calls with observed UMI support: {summary['n_one_umi_only_ab']:,}/{summary['n_ab_umi_observed']:,} ({summary['fraction_one_umi_only_among_umi_observed_ab']:.2%})
- Conservative residual AB-compatible cells: {summary['n_high_confidence_residual_ab']:,}

The dominant detected explanation is unsafe or suspicious upstream TCR joining, not one-UMI ambient contamination or UMAP failure. Rebuild affected TCR assignments with `sample_id + barcode_core` before defining NK truth or resuming V4.2 model training.

## Outputs

- HTML: `{STATIC_DIR / 'index.html'}`
- PDF: `{STATIC_DIR / 'gdtai_v4_2_tcr_conflict_audit_report.pdf'}`
- Source forensics: `{TABLE_DIR / 'source_join_integrity_audit.csv'}`
- Boundary source summary: `{TABLE_DIR / 'boundary_tcr_support_by_source.csv'}`
- Cell-level compact audit: `{TABLE_DIR / 'boundary_tcr_conflict_cell_audit.csv.gz'}`
"""
    (LOG_DIR / "gdtai_v4_2_tcr_conflict_audit_summary.md").write_text(markdown, encoding="utf-8")


def render_pdf() -> None:
    html_path = STATIC_DIR / "index.html"
    pdf_path = STATIC_DIR / "gdtai_v4_2_tcr_conflict_audit_report.pdf"
    chrome = shutil.which("google-chrome") or shutil.which("chromium") or shutil.which("chromium-browser")
    if not chrome:
        raise RuntimeError("No headless Chrome executable is available for PDF rendering")
    command = [
        chrome,
        "--headless",
        "--disable-gpu",
        "--no-sandbox",
        "--allow-file-access-from-files",
        f"--print-to-pdf={pdf_path}",
        "--no-pdf-header-footer",
        html_path.resolve().as_uri(),
    ]
    completed = subprocess.run(command, capture_output=True, text=True, timeout=180, check=False)
    if completed.returncode != 0 or not pdf_path.is_file() or pdf_path.stat().st_size == 0:
        raise RuntimeError(f"PDF rendering failed: {completed.stderr[-2000:]}")


def main() -> None:
    args = parse_args()
    for directory in (TABLE_DIR, FIGURE_DIR, LOG_DIR, STATIC_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    config = load_runtime_config(args.config)
    obs, boundary_indices, labels, sample, genes, marker_matrix = load_boundary_context(config)
    details = extract_boundary_tcr_details(config, obs, boundary_indices)
    details = add_chain_and_expression_features(details, labels, marker_matrix, genes)

    mapping_path = PROJECT_ROOT / "analysis_26GSE_V4/reports/tcr_mapping_to_h5ad_report.csv"
    mapping = pd.read_csv(mapping_path, dtype="string") if mapping_path.is_file() else pd.DataFrame()
    mapping_modes = (
        dict(zip(mapping["GSE"], mapping["mapping_mode"], strict=False))
        if not mapping.empty
        else {}
    )
    role_paths: dict[str, str] = {}
    roles = development_roles(config)
    for role in roles.itertuples(index=False):
        if role.cohort_id != "current_atlas":
            role_paths[role.cohort_id] = str(role.path)
    source_audit_path = TABLE_DIR / "source_join_integrity_audit.csv"
    if args.skip_source_rescan and source_audit_path.is_file():
        source_audit = pd.read_csv(source_audit_path)
    else:
        source_audit = source_h5ad_audit(details["source_gse_id"].unique(), mapping_modes, role_paths)
        source_audit.to_csv(source_audit_path, index=False)

    boundary_reuse = boundary_pair_reuse(details)
    boundary_reuse.to_csv(TABLE_DIR / "boundary_exact_pair_reuse_by_source.csv", index=False)
    source_audit = source_audit.merge(boundary_reuse, on="source_gse_id", how="outer")
    source_status = source_audit.set_index("source_gse_id")["audit_status"].fillna("source_audit_unresolved")
    details["join_audit_status"] = details["source_gse_id"].map(source_status).fillna("source_audit_unresolved")
    details["join_red_flag"] = details["join_audit_status"].isin(
        [
            "documented_unsafe_barcode_join",
            "highly_suspicious_cross_donor_reuse",
            "highly_suspicious_cross_sample_reuse",
            "strong_author_lineage_conflict",
        ]
    )
    details["ab_from_join_red_flag_source"] = details["productive_ab"] & details["join_red_flag"]
    details["ab_provenance"] = "No TRA/TRB"
    details.loc[details["productive_ab"] & details["join_red_flag"], "ab_provenance"] = "AB from red-flag source"
    details.loc[details["productive_ab"] & ~details["join_red_flag"], "ab_provenance"] = "AB without detected join red flag"

    source_summary, cluster_summary, umi_summary = summarize_boundary(details)
    source_audit_for_merge = source_audit.drop(
        columns=[
            "n_boundary_cells",
            "n_boundary_productive_ab",
            "fraction_boundary_productive_ab",
        ],
        errors="ignore",
    )
    source_summary = source_summary.merge(
        source_audit_for_merge, on="source_gse_id", how="left"
    )
    source_summary["n_boundary_productive_ab"] = source_summary["n_productive_ab"]
    source_summary["fraction_boundary_productive_ab"] = source_summary[
        "fraction_productive_ab"
    ]
    source_summary.to_csv(TABLE_DIR / "boundary_tcr_support_by_source.csv", index=False)
    cluster_summary.to_csv(TABLE_DIR / "boundary_tcr_support_by_cluster.csv", index=False)
    umi_summary.to_csv(TABLE_DIR / "boundary_umi_support_tiers_by_source.csv", index=False)
    stats_table = statistical_audit(details)
    stats_table.to_csv(TABLE_DIR / "boundary_expression_tcr_statistical_tests.csv", index=False)
    association = cluster_association_audit(details)
    association.to_csv(TABLE_DIR / "boundary_cluster_association_metrics.csv", index=False)

    hvg = pd.read_csv(stage_paths(config)["table"] / "refined_hvg_selection.csv")
    forced = list(config["hvg"]["forced_lineage_genes"])
    hvg_audit = pd.DataFrame(
        {
            "gene": forced,
            "selected_for_refined_scvi": [
                bool(hvg.set_index("gene")["selected_for_refined_scvi"].get(gene, False)) for gene in forced
            ],
            "role": [
                "T_lineage" if gene in T_MARKERS else "NK_or_context" for gene in forced
            ],
        }
    )
    hvg_audit.to_csv(TABLE_DIR / "forced_lineage_hvg_audit.csv", index=False)

    compact_columns = [
        "integration_cell_id",
        "source_gse_id",
        "input_cohort_id",
        "source_original_cell_id",
        "boundary_cluster",
        "primary_nk_anchor",
        "productive_tra",
        "productive_trb",
        "productive_trd",
        "productive_ab",
        "paired_ab",
        "tra_umi",
        "trb_umi",
        "tra_reads",
        "trb_reads",
        "umi_support_tier",
        "t_marker_hits",
        "t_signal_hits",
        "cd3_hits",
        "tcr_transcript_hits",
        "nk_marker_hits",
        "expression_class",
        "join_audit_status",
        "ab_provenance",
    ]
    details[compact_columns].to_csv(
        TABLE_DIR / "boundary_tcr_conflict_cell_audit.csv.gz",
        index=False,
        compression="gzip",
    )

    ab = details["productive_ab"]
    observed = ab & details["any_ab_umi_known"]
    high_confidence = (
        ab
        & details["paired_ab"]
        & ~details["join_red_flag"]
        & details["any_ab_umi_ge2"]
        & details["multi_gene_t_support"]
    )
    compatible_missing = (
        ab
        & details["paired_ab"]
        & ~details["join_red_flag"]
        & ~details["any_ab_umi_known"]
        & details["multi_gene_t_support"]
    )
    summary = {
        "status": "PASS",
        "protocol": "v4.2-boundary-tcr-conflict-audit-v1",
        "n_boundary_cells": int(details.shape[0]),
        "n_productive_ab": int(ab.sum()),
        "fraction_productive_ab": float(ab.mean()),
        "n_ab_from_red_flag_sources": int(details["ab_from_join_red_flag_source"].sum()),
        "fraction_ab_from_red_flag_sources": float(details.loc[ab, "join_red_flag"].mean()),
        "n_ab_umi_observed": int(observed.sum()),
        "n_ab_umi_unavailable": int((ab & ~details["any_ab_umi_known"]).sum()),
        "n_one_umi_only_ab": int((ab & details["one_umi_only"]).sum()),
        "fraction_one_umi_only_among_umi_observed_ab": float(
            details.loc[observed, "one_umi_only"].mean()
        ) if observed.any() else float("nan"),
        "n_ab_any_umi_ge2": int((ab & details["any_ab_umi_ge2"]).sum()),
        "n_ab_reads_ge10": int((ab & details["any_ab_reads_ge10"]).sum()),
        "n_high_confidence_residual_ab": int(high_confidence.sum()),
        "n_compatible_but_umi_missing_residual_ab": int(compatible_missing.sum()),
        "n_selected_hvgs": int(hvg["selected_for_refined_scvi"].sum()),
        "n_forced_lineage_expected": len(forced),
        "n_forced_lineage_selected": int(hvg_audit["selected_for_refined_scvi"].sum()),
        "n_sources_with_join_red_flag": int(source_audit["join_red_flag"].fillna(False).sum()),
        "n_sources_audited": int(source_audit.shape[0]),
        "source_h5ads_modified": False,
        "script_sha256": sha256_file(SCRIPT_PATH),
        "config_sha256": config["_config_sha256"],
    }
    atomic_write_json(LOG_DIR / "summary.json", summary)

    plot_umap_provenance(details, sample, FIGURE_DIR / "boundary_umap_tcr_provenance.png")
    plot_source_forensics(source_summary, FIGURE_DIR / "source_join_forensics.png")
    plot_umi_support(umi_summary, source_summary, FIGURE_DIR / "umi_support_by_source.png")
    plot_cluster_diagnostics(cluster_summary, FIGURE_DIR / "cluster_marker_tcr_diagnostics.png")
    build_report(details, source_summary, cluster_summary, stats_table, association, hvg_audit, summary)
    render_pdf()
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
