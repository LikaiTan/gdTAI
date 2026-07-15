#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Build an evidence-first cluster annotation review package.

The goal is not to overwrite the milestone object.  This script produces
cluster-level evidence tables and proposed labels from independent evidence:

- curated marker-signature expression on all cells, streamed from the H5AD
- Phase 4 TRA/TRB/TRD scores
- productive TCR-chain metadata when present
- old/simple/scANVI labels as weak priors only
- dataset and tissue composition checks for cluster purity
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
from dataclasses import dataclass
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse
from sklearn.utils import sparsefuncs


TARGET_SUM = 10_000.0
CHUNK_SIZE = 50_000
FIGURE_DPI = 300


PROJECT_ROOT = _TNK_PROJECT_ROOT
DEFAULT_INPUT = PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "integrated_plus6.h5ad"
TABLE_DIR = PROJECT_ROOT / "Integrated_dataset" / "tables" / "annotation_evidence_review"
FIGURE_DIR = PROJECT_ROOT / "Integrated_dataset" / "figures" / "annotation_evidence_review"
LOG_DIR = PROJECT_ROOT / "Integrated_dataset" / "logs" / "annotation_evidence_review"


SIGNATURES: dict[str, list[str]] = {
    "pan_t": ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "LCK", "LAT"],
    "cd4_helper": ["CD4", "IL7R", "LTB", "MAL", "CCR7", "SELL", "TCF7", "LEF1"],
    "cd8_t": ["CD8A", "CD8B", "CCL5", "GZMK", "CXCR3"],
    "cytotoxic": ["NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "GZMA", "CST7", "CTSW"],
    "nk": ["NKG7", "GNLY", "KLRD1", "KLRC1", "KLRC2", "NCR1", "FCGR3A", "TYROBP", "FCER1G", "XCL1", "XCL2"],
    "treg": ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT", "TNFRSF18", "CCR8"],
    "gdt": ["TRDC", "TRDV1", "TRDV2", "TRGV9", "TRGC1", "TRGC2"],
    "mait": ["TRAV1-2", "KLRB1", "SLC4A10", "DPP4", "IL18RAP", "ZBTB16"],
    "ilc": ["KIT", "IL7R", "KLRB1", "RORA", "GATA3", "IL1R1", "AREG"],
    "cycling": ["MKI67", "TOP2A", "STMN1", "TUBA1B", "HMGB2", "PCNA"],
    "myeloid": ["LYZ", "LST1", "S100A8", "S100A9", "FCN1", "CD14", "MS4A7", "C1QA", "C1QB", "C1QC"],
    "dc": ["FCER1A", "CLEC10A", "CD1C", "LILRA4", "CLEC4C", "IRF7"],
    "b_cell": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD74"],
    "plasma": ["MZB1", "XBP1", "JCHAIN", "IGKC", "IGHG1"],
    "platelet": ["PPBP", "PF4", "NRGN", "GNG11", "GP9", "ITGA2B", "TUBB1"],
    "erythroid": ["HBB", "HBA1", "HBA2", "AHSP", "ALAS2"],
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "MUC1", "TACSTD2"],
    "endothelial": ["PECAM1", "VWF", "KDR", "ENG"],
    "fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "COL3A1"],
}


@dataclass(frozen=True)
class ObsArrays:
    leiden: np.ndarray
    source_gse_id: np.ndarray
    tissue: np.ndarray
    old_label: np.ndarray | None
    scanvi_superclass: np.ndarray | None
    scanvi_label: np.ndarray | None
    paired_ab: np.ndarray
    paired_gd: np.ndarray
    any_ab: np.ndarray
    any_gd: np.ndarray
    sorted_gdt: np.ndarray
    phase4_trd_score: np.ndarray
    phase4_trab_score: np.ndarray
    phase4_trd_minus_trab: np.ndarray


def ensure_dirs() -> None:
    for path in (TABLE_DIR, FIGURE_DIR, LOG_DIR):
        path.mkdir(parents=True, exist_ok=True)


def cluster_sort_key(value: object) -> tuple[int, str]:
    text = str(value)
    return (0, f"{int(text):06d}") if text.isdigit() else (1, text)


def decode_string_array(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values],
        dtype=object,
    )


def read_obs_column(handle: h5py.File, column: str, default: object | None = None) -> np.ndarray:
    obs = handle["obs"]
    n_obs = int(obs["_index"].shape[0])
    if column not in obs:
        if default is None:
            raise KeyError(column)
        return np.full(n_obs, default, dtype=object)
    obj = obs[column]
    if isinstance(obj, h5py.Group) and obj.attrs.get("encoding-type") == "categorical":
        categories = decode_string_array(obj["categories"][:])
        codes = obj["codes"][:]
        out = np.full(codes.shape, "", dtype=object)
        valid = codes >= 0
        out[valid] = categories[codes[valid]]
        return out
    raw = obj[:]
    if raw.dtype.kind in {"S", "O", "U"}:
        return decode_string_array(raw)
    return np.asarray(raw)


def read_numeric_obs(handle: h5py.File, column: str) -> np.ndarray:
    n_obs = int(handle["obs"]["_index"].shape[0])
    if column not in handle["obs"]:
        return np.full(n_obs, np.nan, dtype=np.float32)
    return np.asarray(read_obs_column(handle, column), dtype=np.float32)


def normalize_bool(values: np.ndarray) -> np.ndarray:
    if values.dtype == bool:
        return values.astype(bool, copy=False)
    text = np.char.lower(values.astype(str))
    return np.isin(text, ["true", "1", "yes", "y"])


def nonempty(values: np.ndarray) -> np.ndarray:
    text = np.char.lower(values.astype(str))
    return ~np.isin(text, ["", "nan", "none", "na", "null"])


def read_bool_or_derive(handle: h5py.File, bool_col: str, chain_cols: tuple[str, str]) -> np.ndarray:
    if bool_col in handle["obs"]:
        return normalize_bool(read_obs_column(handle, bool_col))
    left = nonempty(read_obs_column(handle, chain_cols[0], default=""))
    right = nonempty(read_obs_column(handle, chain_cols[1], default=""))
    return left & right


def load_obs_arrays(handle: h5py.File) -> ObsArrays:
    n_obs = int(handle["obs"]["_index"].shape[0])
    paired_ab = read_bool_or_derive(handle, "has_TRA_TRB_paired", ("TRA_cdr3", "TRB_cdr3"))
    paired_gd = read_bool_or_derive(handle, "has_TRG_TRD_paired", ("TRG_cdr3", "TRD_cdr3"))
    any_ab = normalize_bool(read_obs_column(handle, "has_any_ab_tcr", default=False))
    if not any_ab.any():
        any_ab = nonempty(read_obs_column(handle, "TRA_cdr3", default="")) | nonempty(read_obs_column(handle, "TRB_cdr3", default=""))
    any_gd = normalize_bool(read_obs_column(handle, "has_any_gd_tcr", default=False))
    if not any_gd.any():
        any_gd = nonempty(read_obs_column(handle, "TRG_cdr3", default="")) | nonempty(read_obs_column(handle, "TRD_cdr3", default=""))

    tissue_col = "tissue_corrected" if "tissue_corrected" in handle["obs"] else "tissue"
    old_label = None
    if "simple_annotation_plus6" in handle["obs"]:
        old_label = read_obs_column(handle, "simple_annotation_plus6")
    elif "simple_annotation" in handle["obs"]:
        old_label = read_obs_column(handle, "simple_annotation")

    scanvi_superclass = read_obs_column(handle, "scanvi_tnk_superclass") if "scanvi_tnk_superclass" in handle["obs"] else None
    scanvi_label = read_obs_column(handle, "scanvi_detailed_label") if "scanvi_detailed_label" in handle["obs"] else None

    return ObsArrays(
        leiden=read_obs_column(handle, "leiden").astype(str),
        source_gse_id=read_obs_column(handle, "source_gse_id", default="unknown").astype(str),
        tissue=read_obs_column(handle, tissue_col, default="unknown").astype(str),
        old_label=old_label.astype(str) if old_label is not None else None,
        scanvi_superclass=scanvi_superclass.astype(str) if scanvi_superclass is not None else None,
        scanvi_label=scanvi_label.astype(str) if scanvi_label is not None else None,
        paired_ab=paired_ab,
        paired_gd=paired_gd,
        any_ab=any_ab,
        any_gd=any_gd,
        sorted_gdt=normalize_bool(read_obs_column(handle, "Sorted_gdT", default=False)),
        phase4_trd_score=read_numeric_obs(handle, "phase4_trd_score"),
        phase4_trab_score=read_numeric_obs(handle, "phase4_trab_score"),
        phase4_trd_minus_trab=read_numeric_obs(handle, "phase4_trd_minus_trab"),
    )


def read_csr_chunk(x_group: h5py.Group, start: int, end: int, n_vars: int) -> sparse.csr_matrix:
    indptr = x_group["indptr"][start : end + 1].astype(np.int64, copy=False)
    offset = int(indptr[0])
    indptr = indptr - offset
    nnz = int(indptr[-1])
    data = x_group["data"][offset : offset + nnz].astype(np.float32, copy=False)
    indices = x_group["indices"][offset : offset + nnz].astype(np.int32, copy=False)
    return sparse.csr_matrix((data, indices, indptr), shape=(end - start, n_vars))


def normalize_log1p_chunk(chunk: sparse.csr_matrix) -> sparse.csr_matrix:
    row_sums = np.asarray(chunk.sum(axis=1)).ravel().astype(np.float32, copy=False)
    scale = np.zeros_like(row_sums, dtype=np.float32)
    positive = row_sums > 0
    scale[positive] = TARGET_SUM / row_sums[positive]
    sparsefuncs.inplace_row_scale(chunk, scale)
    np.log1p(chunk.data, out=chunk.data)
    return chunk


def value_counts_for_mask(values: np.ndarray, mask: np.ndarray) -> tuple[str, float, int]:
    if not np.any(mask):
        return "NA", 0.0, 0
    counts = pd.Series(values[mask].astype(str)).value_counts(dropna=False)
    top_value = str(counts.index[0])
    top_count = int(counts.iloc[0])
    return top_value, float(top_count / int(mask.sum())), top_count


def build_cluster_base(obs: ObsArrays) -> tuple[pd.DataFrame, dict[str, int], np.ndarray]:
    clusters = sorted(pd.unique(obs.leiden), key=cluster_sort_key)
    cluster_to_code = {cluster: idx for idx, cluster in enumerate(clusters)}
    codes = pd.Series(obs.leiden).map(cluster_to_code).to_numpy(dtype=np.int32)
    rows: list[dict[str, object]] = []
    for cluster in clusters:
        mask = obs.leiden == cluster
        n_cells = int(mask.sum())
        top_gse, top_gse_fraction, _ = value_counts_for_mask(obs.source_gse_id, mask)
        top_tissue, top_tissue_fraction, _ = value_counts_for_mask(obs.tissue, mask)
        row: dict[str, object] = {
            "leiden": cluster,
            "n_cells": n_cells,
            "n_gses": int(pd.unique(obs.source_gse_id[mask]).size),
            "top_gse": top_gse,
            "top_gse_fraction": top_gse_fraction,
            "top_tissue": top_tissue,
            "top_tissue_fraction": top_tissue_fraction,
            "paired_ab_fraction": float(obs.paired_ab[mask].mean()) if n_cells else 0.0,
            "paired_gd_fraction": float(obs.paired_gd[mask].mean()) if n_cells else 0.0,
            "any_ab_fraction": float(obs.any_ab[mask].mean()) if n_cells else 0.0,
            "any_gd_fraction": float(obs.any_gd[mask].mean()) if n_cells else 0.0,
            "paired_ab_and_gd_fraction": float((obs.paired_ab[mask] & obs.paired_gd[mask]).mean()) if n_cells else 0.0,
            "sorted_gdt_fraction": float(obs.sorted_gdt[mask].mean()) if n_cells else 0.0,
            "phase4_trd_score_median": float(np.nanmedian(obs.phase4_trd_score[mask])),
            "phase4_trab_score_median": float(np.nanmedian(obs.phase4_trab_score[mask])),
            "phase4_trd_minus_trab_median": float(np.nanmedian(obs.phase4_trd_minus_trab[mask])),
            "phase4_trd_minus_trab_p75": float(np.nanpercentile(obs.phase4_trd_minus_trab[mask], 75)),
        }
        if obs.old_label is not None:
            old_label, old_fraction, _ = value_counts_for_mask(obs.old_label, mask)
            row["old_dominant_label"] = old_label
            row["old_dominant_label_fraction"] = old_fraction
        if obs.scanvi_superclass is not None:
            superclass, superclass_fraction, _ = value_counts_for_mask(obs.scanvi_superclass, mask)
            row["scanvi_dominant_superclass"] = superclass
            row["scanvi_dominant_superclass_fraction"] = superclass_fraction
        if obs.scanvi_label is not None:
            label, label_fraction, _ = value_counts_for_mask(obs.scanvi_label, mask)
            row["scanvi_dominant_label"] = label
            row["scanvi_dominant_label_fraction"] = label_fraction
        rows.append(row)
    return pd.DataFrame(rows), cluster_to_code, codes


def robust_z(values: pd.Series) -> pd.Series:
    arr = values.to_numpy(dtype=float)
    med = np.nanmedian(arr)
    mad = np.nanmedian(np.abs(arr - med))
    if not np.isfinite(mad) or mad <= 0:
        std = np.nanstd(arr)
        scale = std if std > 0 else 1.0
        return pd.Series((arr - np.nanmean(arr)) / scale, index=values.index)
    return pd.Series(0.6745 * (arr - med) / mad, index=values.index)


def compute_marker_stats(
    h5ad_path: Path,
    obs: ObsArrays,
    cluster_to_code: dict[str, int],
    cluster_codes: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    marker_genes = sorted({gene for genes in SIGNATURES.values() for gene in genes})
    with h5py.File(h5ad_path, "r") as handle:
        var_names = pd.Index(decode_string_array(handle["var"]["_index"][:]).astype(str))
        present_genes = [gene for gene in marker_genes if gene in var_names]
        if not present_genes:
            raise ValueError("No requested marker genes were found in the input H5AD.")
        gene_idx = np.asarray([int(var_names.get_loc(gene)) for gene in present_genes], dtype=np.int32)
        n_obs = int(handle["obs"]["_index"].shape[0])
        n_vars = int(var_names.size)
        n_clusters = len(cluster_to_code)
        expr_sum = np.zeros((n_clusters, len(present_genes)), dtype=np.float64)
        detect_sum = np.zeros((n_clusters, len(present_genes)), dtype=np.int64)
        cluster_counts = np.bincount(cluster_codes, minlength=n_clusters).astype(np.int64)
        x_group = handle["X"]
        if not isinstance(x_group, h5py.Group) or not {"data", "indices", "indptr"}.issubset(x_group.keys()):
            raise ValueError("This script currently expects CSR-backed X in the H5AD.")
        for start in range(0, n_obs, CHUNK_SIZE):
            end = min(start + CHUNK_SIZE, n_obs)
            raw_chunk = read_csr_chunk(x_group, start, end, n_vars)
            normalized = normalize_log1p_chunk(raw_chunk.copy())[:, gene_idx].toarray()
            detected = raw_chunk[:, gene_idx].copy()
            detected.data[:] = 1.0
            detected_arr = detected.toarray()
            codes = cluster_codes[start:end]
            for code in np.unique(codes):
                mask = codes == code
                expr_sum[code] += normalized[mask].sum(axis=0)
                detect_sum[code] += detected_arr[mask].sum(axis=0).astype(np.int64)
            if start == 0 or end == n_obs or (start // CHUNK_SIZE) % 20 == 0:
                print(f"processed {end:,}/{n_obs:,} cells", flush=True)

    ordered_clusters = sorted(cluster_to_code, key=cluster_sort_key)
    rows: list[dict[str, object]] = []
    safe_counts = np.maximum(cluster_counts, 1)
    for cluster in ordered_clusters:
        code = cluster_to_code[cluster]
        for pos, gene in enumerate(present_genes):
            rows.append(
                {
                    "leiden": cluster,
                    "gene": gene,
                    "mean_log1p_cpm": float(expr_sum[code, pos] / safe_counts[code]),
                    "detected_fraction": float(detect_sum[code, pos] / safe_counts[code]),
                }
            )
    gene_stats = pd.DataFrame(rows)

    signature_rows: list[dict[str, object]] = []
    gene_mean = gene_stats.pivot(index="leiden", columns="gene", values="mean_log1p_cpm")
    gene_detect = gene_stats.pivot(index="leiden", columns="gene", values="detected_fraction")
    for cluster in ordered_clusters:
        row: dict[str, object] = {"leiden": cluster}
        for signature, genes in SIGNATURES.items():
            present = [gene for gene in genes if gene in gene_mean.columns]
            row[f"{signature}_genes_present"] = len(present)
            if present:
                row[f"{signature}_score"] = float(gene_mean.loc[cluster, present].mean())
                row[f"{signature}_detected_fraction"] = float(gene_detect.loc[cluster, present].mean())
            else:
                row[f"{signature}_score"] = np.nan
                row[f"{signature}_detected_fraction"] = np.nan
        signature_rows.append(row)
    signature_df = pd.DataFrame(signature_rows)
    for signature in SIGNATURES:
        signature_df[f"{signature}_z"] = robust_z(signature_df[f"{signature}_score"])
    return gene_stats, signature_df


def pick_label(row: pd.Series) -> tuple[str, str, str]:
    reasons: list[str] = []
    warnings: list[str] = []

    pan_t = float(row.get("pan_t_z", 0.0))
    pan_t_score = float(row.get("pan_t_score", 0.0))
    cd4 = float(row.get("cd4_helper_z", 0.0))
    cd8 = float(row.get("cd8_t_z", 0.0))
    cyt = float(row.get("cytotoxic_z", 0.0))
    nk = float(row.get("nk_z", 0.0))
    treg = float(row.get("treg_z", 0.0))
    gdt = float(row.get("gdt_z", 0.0))
    mait = float(row.get("mait_z", 0.0))
    ilc = float(row.get("ilc_z", 0.0))
    cycling = float(row.get("cycling_z", 0.0))
    myeloid = float(row.get("myeloid_z", 0.0))
    dc = float(row.get("dc_z", 0.0))
    b_cell = float(row.get("b_cell_z", 0.0))
    plasma = float(row.get("plasma_z", 0.0))
    platelet = float(row.get("platelet_z", 0.0))
    erythroid = float(row.get("erythroid_z", 0.0))
    epithelial = float(row.get("epithelial_z", 0.0))
    endothelial = float(row.get("endothelial_z", 0.0))
    fibroblast = float(row.get("fibroblast_z", 0.0))
    paired_ab = float(row.get("paired_ab_fraction", 0.0))
    paired_gd = float(row.get("paired_gd_fraction", 0.0))
    any_gd = float(row.get("any_gd_fraction", 0.0))
    doublet = float(row.get("paired_ab_and_gd_fraction", 0.0))
    sorted_gdt = float(row.get("sorted_gdt_fraction", 0.0))
    trd_minus = float(row.get("phase4_trd_minus_trab_median", 0.0))
    top_gse_fraction = float(row.get("top_gse_fraction", 0.0))

    if top_gse_fraction >= 0.80:
        warnings.append("dataset_dominated_cluster")
    if doublet >= 0.05:
        warnings.append("high_ab_plus_gd_doublet_proxy")

    contaminant_scores = {
        "myeloid_contaminant": (
            max(myeloid, dc),
            max(float(row.get("myeloid_score", 0.0)), float(row.get("dc_score", 0.0))),
        ),
        "b_or_plasma_contaminant": (
            max(b_cell, plasma),
            max(float(row.get("b_cell_score", 0.0)), float(row.get("plasma_score", 0.0))),
        ),
        "platelet_contaminant": (platelet, float(row.get("platelet_score", 0.0))),
        "erythroid_contaminant": (erythroid, float(row.get("erythroid_score", 0.0))),
        "epithelial_contaminant": (epithelial, float(row.get("epithelial_score", 0.0))),
        "stromal_contaminant": (
            max(endothelial, fibroblast),
            max(float(row.get("endothelial_score", 0.0)), float(row.get("fibroblast_score", 0.0))),
        ),
    }
    contaminant_label, (contaminant_z, contaminant_score) = max(contaminant_scores.items(), key=lambda item: item[1][0])
    if contaminant_z >= 2.0 and contaminant_score >= 0.30 and pan_t_score < 0.80 and paired_gd < 0.05 and sorted_gdt < 0.05:
        if paired_ab >= 0.05:
            warnings.append("TCR_metadata_conflicts_with_low_panT_expression")
        reasons.append(
            f"{contaminant_label}_z={contaminant_z:.2f}; {contaminant_label}_score={contaminant_score:.2f}; "
            f"pan_t_score={pan_t_score:.2f}; paired_ab={paired_ab:.3f}"
        )
        return contaminant_label, "high" if contaminant_z >= 1.8 else "medium", "; ".join(reasons + warnings)

    if sorted_gdt >= 0.25 or paired_gd >= 0.10:
        reasons.append(
            f"gd evidence paired_gd={paired_gd:.3f}, sorted_gdT={sorted_gdt:.3f}, TRD-TRAB={trd_minus:.3f}, gdt_z={gdt:.2f}"
        )
        confidence = "high" if paired_gd >= 0.25 or sorted_gdt >= 0.50 or trd_minus >= 0.30 else "medium"
        return "gamma_delta_T", confidence, "; ".join(reasons + warnings)
    if trd_minus >= 0.15 and gdt >= 1.0 and any_gd >= 0.03:
        reasons.append(
            f"gd-enriched evidence paired_gd={paired_gd:.3f}, any_gd={any_gd:.3f}, TRD-TRAB={trd_minus:.3f}, gdt_z={gdt:.2f}"
        )
        return "gamma_delta_enriched_T", "medium", "; ".join(reasons + warnings)

    if (pan_t_score >= 0.90 or paired_ab >= 0.10) and treg >= 1.0 and cd4 >= cd8 - 0.2 and nk < 1.0:
        reasons.append(f"treg_z={treg:.2f}; cd4_z={cd4:.2f}; cd8_z={cd8:.2f}; nk_z={nk:.2f}")
        return "Treg", "high" if treg >= 1.5 else "medium", "; ".join(reasons + warnings)

    if pan_t_score >= 0.90 and mait >= 2.0 and paired_ab >= 0.05 and paired_gd < 0.05:
        reasons.append(f"MAIT marker z={mait:.2f}; paired_gd={paired_gd:.3f}")
        return "MAIT_like_T", "medium", "; ".join(reasons + warnings)

    if nk >= 0.8 and pan_t_score < 0.90 and paired_ab < 0.10:
        reasons.append(f"nk_z={nk:.2f}; pan_t_score={pan_t_score:.2f}; paired_ab={paired_ab:.3f}")
        return "NK_cell", "high" if nk >= 1.5 else "medium", "; ".join(reasons + warnings)

    if pan_t_score >= 0.90 or paired_ab >= 0.20:
        if cd8 + cyt >= cd4 + 0.5:
            label = "CD8_cytotoxic_T" if cyt >= 0.5 or nk >= 0.5 else "CD8_T"
            reasons.append(
                f"pan_t_score={pan_t_score:.2f}; paired_ab={paired_ab:.3f}; "
                f"cd8_z={cd8:.2f}; cytotoxic_z={cyt:.2f}; nk_z={nk:.2f}"
            )
            confidence = "high" if paired_ab >= 0.20 and cd8 > cd4 else "medium"
            return label, confidence, "; ".join(reasons + warnings)
        if cd4 >= cd8 - 0.2:
            label = "CD4_T"
            reasons.append(f"pan_t_score={pan_t_score:.2f}; paired_ab={paired_ab:.3f}; cd4_z={cd4:.2f}; cd8_z={cd8:.2f}")
            confidence = "high" if paired_ab >= 0.20 and cd4 >= cd8 else "medium"
            return label, confidence, "; ".join(reasons + warnings)
        reasons.append(f"mixed T evidence pan_t_score={pan_t_score:.2f}; paired_ab={paired_ab:.3f}; cd4_z={cd4:.2f}; cd8_z={cd8:.2f}; nk_z={nk:.2f}")
        return "mixed_T", "low", "; ".join(reasons + warnings)

    if paired_ab >= 0.05:
        reasons.append(
            f"TCR metadata without enough pan-T marker support: pan_t_score={pan_t_score:.2f}; "
            f"paired_ab={paired_ab:.3f}; cd4_z={cd4:.2f}; cd8_z={cd8:.2f}"
        )
        return "TCR_metadata_mixed_low_T_marker", "low", "; ".join(reasons + warnings)

    if nk >= 0.5 or ilc >= 1.0:
        reasons.append(f"NK/ILC evidence nk_z={nk:.2f}; ilc_z={ilc:.2f}; pan_t_z={pan_t:.2f}; paired_ab={paired_ab:.3f}")
        return "NK_or_ILC_like", "low", "; ".join(reasons + warnings)

    reasons.append("no strong marker/TCR axis")
    label = "cycling_low_confidence" if cycling >= 1.2 else "unresolved"
    return label, "low", "; ".join(reasons + warnings)


def add_proposed_labels(cluster_df: pd.DataFrame) -> pd.DataFrame:
    labels = []
    confidence = []
    evidence = []
    for _, row in cluster_df.iterrows():
        label, conf, why = pick_label(row)
        labels.append(label)
        confidence.append(conf)
        evidence.append(why)
    out = cluster_df.copy()
    out["proposed_annotation"] = labels
    out["annotation_confidence"] = confidence
    out["annotation_evidence"] = evidence
    if "old_dominant_label" in out.columns:
        old_norm = out["old_dominant_label"].astype(str).str.replace("_cell", "", regex=False).str.replace("_T", "_T", regex=False)
        out["old_label_changed"] = out["old_dominant_label"].astype(str) != out["proposed_annotation"].astype(str)
    return out


def write_heatmap(proposals: pd.DataFrame, target_name: str) -> None:
    z_cols = [f"{name}_z" for name in SIGNATURES]
    heatmap = proposals.set_index("leiden")[z_cols].rename(columns={f"{name}_z": name for name in SIGNATURES})
    height = max(7, 0.35 * heatmap.shape[0])
    width = max(12, 0.55 * heatmap.shape[1])
    fig, ax = plt.subplots(figsize=(width, height), constrained_layout=True)
    sns.heatmap(heatmap, cmap="vlag", center=0, linewidths=0.2, ax=ax)
    ax.set_title(f"{target_name}: cluster marker-signature robust z-scores")
    fig.savefig(FIGURE_DIR / f"{target_name}_cluster_signature_heatmap.png", dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def write_summary(proposals: pd.DataFrame, target_name: str, input_path: Path, gene_stats: pd.DataFrame) -> None:
    changed = proposals.loc[proposals.get("old_label_changed", False) == True] if "old_label_changed" in proposals.columns else pd.DataFrame()
    low_conf = proposals.loc[proposals["annotation_confidence"] == "low"]
    lines = [
        "# Cluster Annotation Evidence Review",
        "",
        f"- Input H5AD: `{input_path}`",
        f"- Target name: `{target_name}`",
        f"- Clusters reviewed: `{proposals.shape[0]}`",
        f"- Marker genes with usable measurements: `{gene_stats['gene'].nunique()}`",
        "- Method: all-cell marker-signature mean/detection summaries, Phase 4 scores, TCR-chain fractions, and dataset/tissue composition.",
        "- Output is review-only; the H5AD annotation column was not modified.",
        "",
        "## Proposed Label Counts",
        "",
    ]
    label_counts = proposals["proposed_annotation"].value_counts().reset_index()
    label_counts.columns = ["proposed_annotation", "n_clusters"]
    lines.append("| proposed_annotation | n_clusters |")
    lines.append("|---|---:|")
    for row in label_counts.itertuples(index=False):
        lines.append(f"| {row.proposed_annotation} | {int(row.n_clusters)} |")
    lines.extend(["", "## Low-Confidence Or Ambiguous Clusters", ""])
    if low_conf.empty:
        lines.append("- none")
    else:
        for row in low_conf.sort_values("n_cells", ascending=False).itertuples(index=False):
            lines.append(
                f"- Leiden `{row.leiden}`: proposed `{row.proposed_annotation}`, n={row.n_cells:,}, evidence: {row.annotation_evidence}"
            )
    if "old_dominant_label" in proposals.columns:
        lines.extend(["", "## Changed From Previous Dominant Label", ""])
        if changed.empty:
            lines.append("- none")
        else:
            for row in changed.sort_values("n_cells", ascending=False).head(30).itertuples(index=False):
                lines.append(
                    f"- Leiden `{row.leiden}`: `{row.old_dominant_label}` -> `{row.proposed_annotation}` "
                    f"(confidence `{row.annotation_confidence}`, n={row.n_cells:,}); {row.annotation_evidence}"
                )
    lines.extend(["", "## Key Outputs", ""])
    lines.append(f"- `{TABLE_DIR / f'{target_name}_cluster_marker_gene_stats.csv'}`")
    lines.append(f"- `{TABLE_DIR / f'{target_name}_cluster_signature_scores.csv'}`")
    lines.append(f"- `{TABLE_DIR / f'{target_name}_cluster_annotation_proposals.csv'}`")
    lines.append(f"- `{FIGURE_DIR / f'{target_name}_cluster_signature_heatmap.png'}`")
    (LOG_DIR / f"{target_name}_annotation_evidence_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def target_name_from_path(path: Path) -> str:
    stem = path.name.removesuffix(".h5ad")
    return stem.replace("integrated_plus6", "plus6").replace("integrated", "current_integrated")


def run(input_h5ad: Path, target_name: str | None = None) -> None:
    ensure_dirs()
    target = target_name or target_name_from_path(input_h5ad)
    with h5py.File(input_h5ad, "r") as handle:
        obs = load_obs_arrays(handle)
    cluster_base, cluster_to_code, cluster_codes = build_cluster_base(obs)
    gene_stats, signature_df = compute_marker_stats(input_h5ad, obs, cluster_to_code, cluster_codes)
    cluster_df = cluster_base.merge(signature_df, on="leiden", how="left", validate="one_to_one")
    proposals = add_proposed_labels(cluster_df)
    proposals = proposals.sort_values(["proposed_annotation", "annotation_confidence", "n_cells"], ascending=[True, True, False])

    gene_stats.to_csv(TABLE_DIR / f"{target}_cluster_marker_gene_stats.csv", index=False)
    signature_df.to_csv(TABLE_DIR / f"{target}_cluster_signature_scores.csv", index=False)
    proposals.to_csv(TABLE_DIR / f"{target}_cluster_annotation_proposals.csv", index=False)
    write_heatmap(proposals, target)
    write_summary(proposals, target, input_h5ad, gene_stats)
    print(f"wrote {TABLE_DIR / f'{target}_cluster_annotation_proposals.csv'}")
    print(f"wrote {LOG_DIR / f'{target}_annotation_evidence_summary.md'}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--target-name", default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run(args.input_h5ad, args.target_name)


if __name__ == "__main__":
    main()
