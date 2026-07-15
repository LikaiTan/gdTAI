#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Small-scale diagnostic for Phase 3 scVI + scANVI settings.

This helper tests the current scANVI query configuration on pre-integration
small libraries. It compares:

1. the current query-batch behavior used in `phase3_scvi_scanvi.py`
2. a query setup that preserves the real `phase3_batch_key`

The goal is to separate scVI/scANVI configuration problems from any issues
introduced later by the merged milestone object.
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


import json
import logging
import time
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import torch

import phase3_scvi_scanvi as phase3
import scvi


DIAG_DIR = _TNK_PROJECT_ROOT / "Integrated_dataset/diagnosis/phase3_scvi_scanvi"
OUTPUT_DIR = DIAG_DIR / "tables"
FIGURE_DIR = DIAG_DIR / "figures"
OBJECT_DIR = DIAG_DIR / "objects"
LOG_PATH = DIAG_DIR / "phase3_scvi_scanvi_diagnostic.log"
SUMMARY_CSV = OUTPUT_DIR / "diagnostic_summary.csv"
LABEL_CSV = OUTPUT_DIR / "diagnostic_label_counts.csv"
CONFIG_JSON = DIAG_DIR / "diagnostic_config.json"
DIAG_H5AD = OBJECT_DIR / "phase3_scvi_scanvi_diagnostic_subset.h5ad"
UMAP_GSE = FIGURE_DIR / "umap_by_gse.png"
UMAP_BATCH = FIGURE_DIR / "umap_by_phase3_batch_key.png"
UMAP_LEIDEN = FIGURE_DIR / "umap_by_leiden.png"
UMAP_SCANVI = FIGURE_DIR / "umap_by_scanvi_label.png"
UMAP_SCANVI_SUPER = FIGURE_DIR / "umap_by_scanvi_superclass.png"

INPUT_H5ADS = [
    _TNK_PROJECT_ROOT / "downloads/per_gse_h5ad_with_metadata/GSE234069_with_tcr.h5ad",
    _TNK_PROJECT_ROOT / "downloads/per_gse_h5ad_with_metadata/GSE240865_with_tcr.h5ad",
]
N_HVGS = 2000
SCVI_EPOCHS = 5
SCANVI_EPOCHS = 3
BATCH_SIZE = 2048
MAX_CELLS = 50_000


def setup_logging() -> None:
    """Configure logging for the diagnostic run."""
    for path in [DIAG_DIR, OUTPUT_DIR, FIGURE_DIR, OBJECT_DIR]:
        path.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(LOG_PATH, mode="w", encoding="utf-8"),
        ],
        force=True,
    )


def normalize_obs(adata: ad.AnnData, source_path: Path) -> ad.AnnData:
    """Add the minimum metadata fields needed for the Phase 3 diagnostic."""
    gse_id = source_path.name.split("_")[0]
    adata.obs["source_gse_id"] = gse_id
    adata.obs["original_cell_id"] = adata.obs_names.astype(str)
    if "sampleid" not in adata.obs.columns and "sample_id" in adata.obs.columns:
        adata.obs["sampleid"] = adata.obs["sample_id"].astype("string")
    if "barcodes" not in adata.obs.columns:
        if "barcode" in adata.obs.columns:
            adata.obs["barcodes"] = adata.obs["barcode"].astype("string")
        else:
            adata.obs["barcodes"] = adata.obs_names.astype("string")
    return adata


def load_subset() -> ad.AnnData:
    """Load a small local subset from pre-integration H5AD libraries."""
    parts: list[ad.AnnData] = []
    for path in INPUT_H5ADS:
        logging.info("Loading %s", path)
        adata = ad.read_h5ad(path)
        adata = normalize_obs(adata, path)
        parts.append(adata)
    subset = ad.concat(parts, join="outer", merge="same", label=None, index_unique=None)
    if subset.n_obs > MAX_CELLS:
        rng = np.random.default_rng(0)
        keep_idx = np.sort(rng.choice(subset.n_obs, size=MAX_CELLS, replace=False))
        subset = subset[keep_idx].copy()
    logging.info("Loaded subset with %s cells and %s genes", subset.n_obs, subset.n_vars)
    return subset


def prepare_subset() -> ad.AnnData:
    """Prepare the combined small-library subset for the diagnostic."""
    adata = load_subset()
    adata, removed = phase3.sanitize_phase3_input(adata)
    logging.info("Sanitization removed %s cells", len(removed))
    phase3.build_batch_key(adata)
    logging.info("Unique phase3 batches: %s", int(adata.obs["phase3_batch_key"].nunique()))
    return adata


def run_small_scvi(adata: ad.AnnData) -> tuple[scvi.model.SCVI, ad.AnnData]:
    """Train a short scVI model on a reduced-HVG subset."""
    hvg_genes, _, _ = phase3.select_hvg_genes(adata, N_HVGS)
    train_adata = adata[:, hvg_genes].copy()
    scvi.model.SCVI.setup_anndata(train_adata, batch_key="phase3_batch_key")
    model = scvi.model.SCVI(
        train_adata,
        n_hidden=128,
        n_layers=2,
        n_latent=30,
        gene_likelihood="nb",
    )
    logging.info("Training small scVI diagnostic model")
    model.train(
        max_epochs=SCVI_EPOCHS,
        batch_size=BATCH_SIZE,
        accelerator="gpu",
        devices=1,
        early_stopping=True,
        early_stopping_patience=3,
    )
    adata.obsm["X_scVI"] = model.get_latent_representation(train_adata).astype(np.float32, copy=False)
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=15)
    sc.tl.leiden(adata, resolution=0.8, key_added="leiden")
    sc.tl.umap(adata)
    return model, train_adata


def build_query_subset(adata: ad.AnnData, batch_mode: str) -> tuple[ad.AnnData, pd.DataFrame]:
    """Prepare the scANVI query subset for a given batch-mode policy."""
    reference_model = phase3.load_reference_model()
    reference_var_names = pd.Index(reference_model.adata.var_names)
    reference_symbol_map = phase3.build_reference_symbol_map()
    query_symbols = pd.Index(adata.var_names.astype(str))
    matched_symbols = query_symbols.intersection(reference_symbol_map.index)
    mapped_ensembl = pd.Index(reference_symbol_map.loc[matched_symbols].astype(str))
    common_var_names = mapped_ensembl.intersection(reference_var_names)
    matched_symbols = [
        symbol for symbol in matched_symbols if reference_symbol_map[symbol] in common_var_names
    ]

    subset_idx, subset_summary = phase3.build_scanvi_subset_indices(
        adata=adata,
        max_cells=min(adata.n_obs, 20000),
        max_per_stratum=1000,
    )
    query_subset = adata[subset_idx, matched_symbols].copy()
    query_subset.var["query_gene_symbol"] = query_subset.var_names.astype(str)
    query_subset.var["ensembl_id"] = reference_symbol_map.loc[matched_symbols].astype(str).to_numpy()
    query_subset.var_names = pd.Index(query_subset.var["ensembl_id"].astype(str))
    duplicate_mask = query_subset.var_names.duplicated()
    if duplicate_mask.any():
        query_subset = query_subset[:, ~duplicate_mask].copy()

    if batch_mode == "constant_query_batch":
        query_subset.obs["batch"] = pd.Categorical(np.repeat("query_batch", query_subset.n_obs))
    elif batch_mode == "phase3_batch_key":
        query_subset.obs["batch"] = pd.Categorical(query_subset.obs["phase3_batch_key"].astype(str))
    else:
        raise ValueError(f"Unknown batch_mode: {batch_mode}")

    return query_subset, subset_summary


def run_scanvi_test(adata: ad.AnnData, batch_mode: str) -> tuple[dict[str, object], pd.DataFrame]:
    """Run a short scANVI test and collect label output statistics."""
    start = time.time()
    result: dict[str, object] = {"batch_mode": batch_mode}
    label_counts = pd.DataFrame()

    try:
        reference_model = phase3.load_reference_model()
        query_subset, subset_summary = build_query_subset(adata, batch_mode=batch_mode)
        phase3_summary = {
            "query_cells": int(query_subset.n_obs),
            "query_genes": int(query_subset.n_vars),
            "query_batches": int(query_subset.obs["batch"].astype(str).nunique()),
            "subset_strata": int(len(subset_summary)),
        }

        scvi.model.SCANVI.prepare_query_anndata(query_subset, reference_model)
        query_model = scvi.model.SCANVI.load_query_data(query_subset, reference_model)
        query_model.train(
            max_epochs=SCANVI_EPOCHS,
            batch_size=BATCH_SIZE,
            accelerator="gpu",
            devices=1,
            plan_kwargs={"weight_decay": 0.0},
            early_stopping=True,
            early_stopping_patience=3,
        )
        pred = pd.Series(query_model.predict(), dtype="string")
        soft = query_model.predict(soft=True)
        conf = soft.max(axis=1).astype(np.float32, copy=False)

        label_counts = (
            pred.value_counts(dropna=False)
            .rename_axis("scanvi_label")
            .reset_index(name="cells")
            .sort_values("cells", ascending=False)
            .reset_index(drop=True)
        )
        label_counts["batch_mode"] = batch_mode

        result.update(
            {
                "status": "success",
                "mean_confidence": float(np.mean(conf)),
                "median_confidence": float(np.median(conf)),
                "labels_n": int(pred.astype(str).nunique()),
                "top_label": str(label_counts.iloc[0]["scanvi_label"]) if len(label_counts) else "",
                "top_label_cells": int(label_counts.iloc[0]["cells"]) if len(label_counts) else 0,
            }
        )
        result.update(phase3_summary)
    except Exception as exc:
        result.update({"status": "error", "error": repr(exc)})
        logging.exception("Diagnostic run failed for batch_mode=%s", batch_mode)
    finally:
        result["elapsed_sec"] = round(time.time() - start, 2)

    return result, label_counts


def collapse_scanvi_superclass(series: pd.Series) -> pd.Series:
    """Collapse detailed scANVI labels to coarse superclasses."""
    return series.astype("string").map(phase3.collapse_scanvi_label).fillna("reference_other")


def plot_umap(adata: ad.AnnData, color: str, output: Path, title: str) -> None:
    """Save one diagnostic UMAP."""
    fig = sc.pl.umap(
        adata,
        color=color,
        show=False,
        return_fig=True,
        frameon=False,
        title=title,
    )
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_h5ad_write_safe(adata: ad.AnnData) -> ad.AnnData:
    """Convert pandas string-backed columns into plain Python strings for H5AD writing."""
    safe = adata.copy()
    for frame in [safe.obs, safe.var]:
        for column in frame.columns:
            if pd.api.types.is_string_dtype(frame[column]):
                frame[column] = frame[column].astype(object)
            elif pd.api.types.is_categorical_dtype(frame[column]):
                categories = frame[column].cat.categories
                if getattr(categories, "dtype", None) is not None and pd.api.types.is_string_dtype(categories.dtype):
                    frame[column] = frame[column].astype(str).astype("category")
    return safe


def main() -> None:
    """Run the small diagnostic comparison."""
    setup_logging()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    CONFIG_JSON.write_text(
        json.dumps(
            {
                "input_h5ads": [str(path) for path in INPUT_H5ADS],
                "max_cells": MAX_CELLS,
                "n_hvgs": N_HVGS,
                "scvi_epochs": SCVI_EPOCHS,
                "scanvi_epochs": SCANVI_EPOCHS,
                "batch_size": BATCH_SIZE,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    adata = prepare_subset()
    _, _ = run_small_scvi(adata)

    rows: list[dict[str, object]] = []
    labels: list[pd.DataFrame] = []
    best_mode = "phase3_batch_key"
    for batch_mode in ["constant_query_batch", "phase3_batch_key"]:
        logging.info("Running diagnostic for batch_mode=%s", batch_mode)
        row, label_frame = run_scanvi_test(adata, batch_mode=batch_mode)
        rows.append(row)
        if len(label_frame):
            labels.append(label_frame)

    summary = pd.DataFrame(rows)
    summary.to_csv(SUMMARY_CSV, index=False)
    if labels:
        pd.concat(labels, ignore_index=True).to_csv(LABEL_CSV, index=False)
    else:
        pd.DataFrame(columns=["scanvi_label", "cells", "batch_mode"]).to_csv(LABEL_CSV, index=False)

    best_row = summary.loc[summary["batch_mode"].eq(best_mode)].iloc[0]
    if best_row["status"] != "success":
        raise RuntimeError(f"Chosen diagnostic mode failed: {best_mode}")

    logging.info("Final plotting mode: %s", best_mode)

    # Re-run once more to store cell-level labels on the object for visualization.
    reference_model = phase3.load_reference_model()
    query_subset, _ = build_query_subset(adata, batch_mode=best_mode)
    scvi.model.SCANVI.prepare_query_anndata(query_subset, reference_model)
    query_model = scvi.model.SCANVI.load_query_data(query_subset, reference_model)
    query_model.train(
        max_epochs=SCANVI_EPOCHS,
        batch_size=BATCH_SIZE,
        accelerator="gpu",
        devices=1,
        plan_kwargs={"weight_decay": 0.0},
        early_stopping=True,
        early_stopping_patience=3,
    )
    pred = pd.Series(query_model.predict(), index=query_subset.obs_names, dtype="string")
    full_labels = pd.Series("not_in_query", index=adata.obs_names, dtype="string")
    full_labels.loc[query_subset.obs_names] = pred
    adata.obs["diagnostic_scanvi_label"] = pd.Categorical(full_labels.astype(str))
    adata.obs["diagnostic_scanvi_superclass"] = pd.Categorical(
        collapse_scanvi_superclass(full_labels)
    )

    plot_umap(adata, "source_gse_id", UMAP_GSE, "Diagnostic UMAP By GSE")
    plot_umap(adata, "phase3_batch_key", UMAP_BATCH, "Diagnostic UMAP By Phase3 Batch Key")
    plot_umap(adata, "leiden", UMAP_LEIDEN, "Diagnostic UMAP By Unsupervised Leiden Cluster")
    plot_umap(adata, "diagnostic_scanvi_label", UMAP_SCANVI, "Diagnostic UMAP By scANVI Label")
    plot_umap(adata, "diagnostic_scanvi_superclass", UMAP_SCANVI_SUPER, "Diagnostic UMAP By scANVI Superclass")
    safe_adata = make_h5ad_write_safe(adata)
    safe_adata.write_h5ad(DIAG_H5AD)
    logging.info("Diagnostic summary written to %s", SUMMARY_CSV)


if __name__ == "__main__":
    main()
