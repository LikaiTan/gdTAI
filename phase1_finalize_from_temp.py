#!/usr/bin/env python3
"""Finalize Phase 1 outputs from existing per-dataset temp candidate files."""

from __future__ import annotations

import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

import phase1_extract_tnk_candidates as p1


# Config
PROJECT_ROOT = Path(__file__).resolve().parent
AUDIT_CSV = PROJECT_ROOT / "Integrated_dataset" / "tables" / "phase0_dataset_audit.csv"


def ensure_temp_x_sparse(temp_path: Path) -> None:
    """Rewrite a temp H5AD so X is CSR if it is currently dense."""
    adata = ad.read_h5ad(temp_path)
    changed = False
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(np.asarray(adata.X))
        changed = True
    elif not sp.isspmatrix_csr(adata.X):
        adata.X = adata.X.tocsr()
        changed = True

    if changed:
        adata.write_h5ad(temp_path)


def annotation_columns_for_dataset(h5ad_path: str) -> str:
    """Recover the annotation columns used for one original dataset."""
    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        return ";".join(p1.find_annotation_columns(list(adata.obs.columns)))
    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()


def marker_availability_from_var_names(gse_id: str, var_names: pd.Index) -> list[dict]:
    """Build marker availability rows from var names only."""
    return [
        {
            "gse_id": gse_id,
            "marker_group": "T",
            "genes_present": int(len(p1.marker_positions_from_var_names(var_names, p1.T_MARKERS))),
            "score_genes_present": int(len(p1.marker_positions_from_var_names(var_names, p1.T_MARKERS))),
        },
        {
            "gse_id": gse_id,
            "marker_group": "GD",
            "genes_present": int(len(p1.marker_positions_from_var_names(var_names, p1.GD_MARKERS))),
            "score_genes_present": int(len(p1.marker_positions_from_var_names(var_names, p1.GD_MARKERS))),
        },
        {
            "gse_id": gse_id,
            "marker_group": "NK",
            "genes_present": int(len(p1.marker_positions_from_var_names(var_names, p1.NK_MARKERS))),
            "score_genes_present": int(len(p1.marker_positions_from_var_names(var_names, p1.NK_MARKERS))),
        },
        {
            "gse_id": gse_id,
            "marker_group": "NK_STRONG",
            "genes_present": int(len(p1.marker_positions_from_var_names(var_names, p1.NK_STRONG_MARKERS))),
            "score_genes_present": int(len(p1.marker_positions_from_var_names(var_names, p1.NK_STRONG_MARKERS))),
        },
        {
            "gse_id": gse_id,
            "marker_group": "CONTAM",
            "genes_present": int(len(p1.marker_positions_from_var_names(var_names, p1.CONTAM_MARKERS))),
            "score_genes_present": int(len(p1.marker_positions_from_var_names(var_names, p1.CONTAM_MARKERS))),
        },
    ]


def summarize_temp_subset(row: pd.Series, temp_path: Path) -> tuple[dict, list[dict]]:
    """Build a Phase 1 summary row and marker availability rows from a temp subset."""
    ensure_temp_x_sparse(temp_path)
    adata = ad.read_h5ad(temp_path, backed="r")
    try:
        obs = adata.obs
        summary = {
            "gse_id": row["gse_id"],
            "input_n_obs": int(row["n_obs"]),
            "input_n_vars": int(row["n_vars"]),
            "candidate_n_obs": int(adata.n_obs),
            "candidate_fraction": float(adata.n_obs / row["n_obs"]) if row["n_obs"] else np.nan,
            "annotation_positive_n": int(obs["phase1_annotation_keep"].sum()),
            "marker_positive_n": int(obs["phase1_marker_keep"].sum()),
            "annotation_and_marker_n": int((obs["phase1_annotation_keep"] & obs["phase1_marker_keep"]).sum()),
            "annotation_only_n": int((obs["phase1_annotation_keep"] & ~obs["phase1_marker_keep"]).sum()),
            "marker_only_n": int((~obs["phase1_annotation_keep"] & obs["phase1_marker_keep"]).sum()),
            "annotation_columns": annotation_columns_for_dataset(row["h5ad_path"]),
            "mean_t_score_candidates": float(obs["phase1_t_score"].mean()) if adata.n_obs else np.nan,
            "mean_nk_score_candidates": float(obs["phase1_nk_score"].mean()) if adata.n_obs else np.nan,
            "mean_gd_score_candidates": float(obs["phase1_gd_score"].mean()) if adata.n_obs else np.nan,
            "mean_contam_score_candidates": float(obs["phase1_contam_score"].mean()) if adata.n_obs else np.nan,
        }
        marker_rows = marker_availability_from_var_names(row["gse_id"], pd.Index(adata.var_names.astype(str)))
        return summary, marker_rows
    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()


def main() -> None:
    """Finalize Phase 1 outputs from temp subset files already on disk."""
    p1.ensure_output_dirs()
    audit = pd.read_csv(AUDIT_CSV)
    audit = audit.loc[audit["phase0_category"] == "A", ["gse_id", "h5ad_path", "n_obs", "n_vars"]]
    audit = audit.sort_values("gse_id").reset_index(drop=True)

    summary_rows: list[dict] = []
    marker_rows: list[dict] = []
    temp_files: dict[str, Path] = {}

    for _, row in audit.iterrows():
        temp_path = p1.TMP_DIR / f"{row['gse_id']}_phase1_candidates.h5ad"
        if not temp_path.exists():
            raise FileNotFoundError(f"Missing temp candidate file: {temp_path}")
        print(f"Finalizing {row['gse_id']} ...", flush=True)
        summary, marker_info = summarize_temp_subset(row, temp_path)
        summary_rows.append(summary)
        marker_rows.extend(marker_info)
        temp_files[row["gse_id"]] = temp_path

    summary_df = pd.DataFrame(summary_rows).sort_values("gse_id").reset_index(drop=True)
    marker_df = pd.DataFrame(marker_rows).sort_values(["gse_id", "marker_group"]).reset_index(drop=True)

    p1.concat_temp_subsets_on_disk(temp_files)
    p1.validate_merged_candidates(summary_df)
    summary_df.to_csv(p1.PHASE1_SUMMARY_CSV, index=False)
    marker_df.to_csv(p1.PHASE1_MARKER_CSV, index=False)
    p1.write_figures(summary_df)
    p1.write_qc_summary(summary_df)
    shutil.rmtree(p1.TMP_DIR, ignore_errors=True)

    print("\nWrote:")
    print(p1.TNK_CANDIDATES_H5AD)
    print(p1.PHASE1_SUMMARY_CSV)
    print(p1.PHASE1_MARKER_CSV)
    print(p1.PHASE1_QC_MD)


if __name__ == "__main__":
    main()
