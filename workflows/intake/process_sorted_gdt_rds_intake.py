#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Convert sorted gdT Seurat RDS objects into harmonized scored H5ADs."""

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
import gzip
import logging
import os
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import io as spio
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = PROJECT_ROOT / "newdata" / "Sorted_gdT"
SORTED_GDT_OUTPUT_DIR = PROJECT_ROOT / "newdata" / "Sorted_gdT"
TABLE_ROOT = PROJECT_ROOT / "Integrated_dataset" / "tables" / "Sorted_gdT"
LOG_ROOT = PROJECT_ROOT / "Integrated_dataset" / "logs" / "Sorted_gdT"
FIGURE_ROOT = PROJECT_ROOT / "Integrated_dataset" / "figures"
R_EXPORT_HELPER = PROJECT_ROOT / "workflows" / "intake" / "supplementary_export_rds_payloads.R"
PHASE4_SCORER = PROJECT_ROOT / "workflows" / "intake" / "phase4_score_single_h5ad.py"
PYTHON_BIN = Path("/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python")
RSCRIPT_BIN = "Rscript"
MAIN_LOG = LOG_ROOT / "sorted_gdt_rds_intake.log"
FIGURE_DPI = 300
RANDOM_STATE = 17

TRA_TRB_COLS = [
    "TRA_cdr3",
    "TRA_v",
    "TRA_d",
    "TRA_j",
    "TRA_cdr3_nt",
    "TRA_clone_id",
    "TRA_umis",
    "TRA_reads",
    "TRA_c_gene",
    "TRB_cdr3",
    "TRB_v",
    "TRB_d",
    "TRB_j",
    "TRB_cdr3_nt",
    "TRB_clone_id",
    "TRB_umis",
    "TRB_reads",
    "TRB_c_gene",
]

TRG_TRD_COLS = [
    "TRG_cdr3",
    "TRG_v",
    "TRG_d",
    "TRG_j",
    "TRG_cdr3_nt",
    "TRG_clone_id",
    "TRG_umis",
    "TRG_reads",
    "TRG_c_gene",
    "TRD_cdr3",
    "TRD_v",
    "TRD_d",
    "TRD_j",
    "TRD_cdr3_nt",
    "TRD_clone_id",
    "TRD_umis",
    "TRD_reads",
    "TRD_c_gene",
]


@dataclass(frozen=True)
class DatasetConfig:
    dataset_id: str
    input_rds: Path
    output_h5ad: Path
    description: str
    default_tissue: str
    sample_rule: str


DATASETS = [
    DatasetConfig(
        dataset_id="GDT_2020AUG_woCOV",
        input_rds=INPUT_DIR / "GDT_2020AUG_woCOV.rds",
        output_h5ad=SORTED_GDT_OUTPUT_DIR / "GDT_2020AUG_woCOV_sorted_gdt.h5ad",
        description="cord blood and peripheral blood sorted gdT cells",
        default_tissue="peripheral blood",
        sample_rule="orig.ident",
    ),
    DatasetConfig(
        dataset_id="GDTlung2023july_7p",
        input_rds=INPUT_DIR / "GDTlung2023july_7p.rds",
        output_h5ad=SORTED_GDT_OUTPUT_DIR / "GDTlung2023july_7p_sorted_gdt.h5ad",
        description="sorted gdT cells from lung and lymph node",
        default_tissue="lung",
        sample_rule="patient+tissue",
    ),
    DatasetConfig(
        dataset_id="MalteGDT",
        input_rds=INPUT_DIR / "MalteGDT.rds",
        output_h5ad=SORTED_GDT_OUTPUT_DIR / "MalteGDT_sorted_gdt.h5ad",
        description="sorted gdT cells from PBMC",
        default_tissue="peripheral blood",
        sample_rule="cellhash_if_available_else_orig.ident",
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset-id",
        action="append",
        default=[],
        help="Optional dataset id(s) to process. Defaults to all three.",
    )
    return parser.parse_args()


def configure_logging() -> None:
    LOG_ROOT.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(MAIN_LOG, mode="w", encoding="utf-8"),
            logging.StreamHandler(),
        ],
        force=True,
    )


def clean_text(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"", "na", "nan", "none", "null", "<na>"} else text


def normalize_string_series(series: pd.Series) -> pd.Series:
    return series.map(clean_text).astype(object)


def normalize_barcode_core(value: object) -> str:
    text = clean_text(value).upper()
    if not text:
        return ""
    text = text.split("_")[0]
    return text.split("-")[0]


def first_nonblank_series(df: pd.DataFrame, columns: list[str], default: str = "") -> pd.Series:
    if isinstance(default, pd.Series):
        out = normalize_string_series(default.reindex(df.index)).astype(object)
    elif isinstance(default, (pd.Index, np.ndarray, list, tuple)):
        out = normalize_string_series(pd.Series(default, index=df.index)).astype(object)
    else:
        out = pd.Series([default] * len(df), index=df.index, dtype=object)
    for column in columns:
        if column not in df.columns:
            continue
        values = normalize_string_series(df[column])
        mask = (out == "") & (values != "")
        out.loc[mask] = values.loc[mask]
    return out


def ensure_unique_obs_names(cell_ids: pd.Series, sample_ids: pd.Series) -> pd.Series:
    cell_ids = normalize_string_series(cell_ids)
    sample_ids = normalize_string_series(sample_ids)
    if not cell_ids.duplicated().any():
        return cell_ids
    out = cell_ids.copy()
    dup_mask = out.duplicated(keep=False)
    out.loc[dup_mask] = sample_ids.loc[dup_mask] + "::" + cell_ids.loc[dup_mask]
    if out.duplicated().any():
        dup_mask = out.duplicated(keep=False)
        out.loc[dup_mask] = out.loc[dup_mask] + "::" + pd.Series(
            np.arange(int(dup_mask.sum())), index=out.index[dup_mask], dtype=object
        ).astype(str)
    return out


def parse_gene_hit(value: object) -> str:
    text = clean_text(value)
    if not text:
        return ""
    first = str(text).split(",")[0].split(";")[0].strip()
    if first.lower() == "none":
        return ""
    return first


def normalize_chain_label(series: pd.Series) -> pd.Series:
    return normalize_string_series(series).str.upper()


def run_command(cmd: list[str], *, cwd: Path | None = None) -> None:
    logging.info("Running command: %s", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)


def export_rds_payload(dataset: DatasetConfig, tmpdir: Path) -> Path:
    prefix = tmpdir / dataset.dataset_id
    run_command(
        [
            RSCRIPT_BIN,
            str(R_EXPORT_HELPER),
            "--mode",
            "seurat",
            "--input",
            str(dataset.input_rds),
            "--out-prefix",
            str(prefix),
            "--assay",
            "RNA",
        ],
        cwd=PROJECT_ROOT,
    )
    return prefix


def read_payload(prefix: Path) -> tuple[sparse.csr_matrix, pd.Index, pd.DataFrame]:
    counts = spio.mmread(prefix.with_name(prefix.name + "_counts.mtx")).tocsr().transpose().tocsr()
    with gzip.open(prefix.with_name(prefix.name + "_genes.tsv.gz"), "rt", encoding="utf-8") as handle:
        genes = pd.read_csv(handle, sep="\t", header=None)
    with gzip.open(prefix.with_name(prefix.name + "_barcodes.tsv.gz"), "rt", encoding="utf-8") as handle:
        barcodes = pd.read_csv(handle, sep="\t", header=None)[0].astype(str)
    with gzip.open(prefix.with_name(prefix.name + "_metadata.csv.gz"), "rt", encoding="utf-8") as handle:
        meta = pd.read_csv(handle)

    meta["cell_id"] = meta["cell_id"].astype(str)
    meta = meta.set_index("cell_id", drop=False).reindex(barcodes.to_numpy(), copy=False)
    if meta["cell_id"].isna().any():
        raise ValueError(f"Metadata rows are missing after barcodes reindex for {prefix.name}")

    var_names = pd.Index(genes.iloc[:, 1].astype(str).to_numpy(), dtype="string")
    return counts, var_names, meta


def build_base_obs(meta: pd.DataFrame, dataset: DatasetConfig) -> pd.DataFrame:
    obs = meta.copy()
    obs["raw_obs_name"] = obs.index.astype(str)
    obs["project name"] = dataset.dataset_id
    obs["gse_id"] = dataset.dataset_id
    obs["source_gse_id"] = dataset.dataset_id
    obs["source_h5ad"] = str(dataset.output_h5ad.resolve())
    obs["source_root"] = str(dataset.input_rds.resolve().parent)
    obs["sample_type"] = "sorted_gdt"
    obs["input_population"] = "purified_gdt"
    obs["enrichment_strategy"] = "sorted_gdt"
    obs["technology_simple"] = "10x 5'"
    obs["assay_type"] = "10x 5' gene expression + VDJ"
    obs["tcr_vdj_flag"] = "yes"
    obs["tcr_availability"] = "gd_embedded"
    obs["tcr_chain_mode"] = "gd_only"
    obs["Sorted_gdT"] = True
    obs["sex"] = ""
    obs["treatment"] = ""
    obs["age"] = ""
    return obs


def build_dataset_specific_obs(meta: pd.DataFrame, dataset: DatasetConfig) -> pd.DataFrame:
    obs = build_base_obs(meta, dataset)
    orig_ident = normalize_string_series(obs.get("orig.ident", pd.Series("", index=obs.index)))
    barcode_source = first_nonblank_series(obs, ["barcode", "bc_backup", "cell_id"], default=obs.index.astype(str))
    barcode_source = normalize_string_series(barcode_source)

    if dataset.dataset_id == "GDT_2020AUG_woCOV":
        sample_id = orig_ident.copy()
        donor_id = orig_ident.str.extract(r"(donor\d+)", expand=False).fillna("")
        tissue = np.where(orig_ident.str.startswith("CB_"), "cord blood", "peripheral blood")
        condition = normalize_string_series(obs.get("group", pd.Series("", index=obs.index)))
        original_label = first_nonblank_series(obs, ["Cell_cluster", "group_cl", "seurat_clusters"])
        library_id = sample_id.copy()
    elif dataset.dataset_id == "GDTlung2023july_7p":
        patient = normalize_string_series(obs.get("patient", pd.Series("", index=obs.index)))
        tissue_raw = normalize_string_series(obs.get("tissue", pd.Series("", index=obs.index)))
        tissue = tissue_raw.replace({"LLN": "lymph node", "Lung": "lung"}).str.lower()
        sample_id = patient.copy()
        sample_id = sample_id.where(patient.eq("") | tissue.eq(""), sample_id + "__" + tissue)
        sample_id = sample_id.where(sample_id.ne(""), orig_ident)
        donor_id = patient.copy()
        condition = ""
        original_label = first_nonblank_series(obs, ["Cell_cluster", "Cluster_tissue", "predicted.id", "seurat_clusters"])
        library_id = orig_ident.copy()
    elif dataset.dataset_id == "MalteGDT":
        cellhash = first_nonblank_series(obs, ["cellhash", "cellhash.x", "cellhash.y"])
        sample_id = normalize_string_series(cellhash)
        sample_id = sample_id.where(sample_id.ne(""), orig_ident)
        donor_id = orig_ident.where(orig_ident.ne(""), "MalteGDT")
        tissue = dataset.default_tissue
        condition = ""
        original_label = first_nonblank_series(obs, ["Cell_cluster", "Top5paired", "seurat_clusters"])
        library_id = orig_ident.where(orig_ident.ne(""), dataset.dataset_id)
    else:
        raise ValueError(f"Unhandled dataset: {dataset.dataset_id}")

    obs["sample_id"] = normalize_string_series(pd.Series(sample_id, index=obs.index))
    obs["sampleid"] = obs["sample_id"].copy()
    obs["library_id"] = normalize_string_series(pd.Series(library_id, index=obs.index))
    obs["donor_id"] = normalize_string_series(pd.Series(donor_id, index=obs.index))
    obs["donor_patient"] = obs["donor_id"].copy()
    obs["tissue"] = normalize_string_series(pd.Series(tissue, index=obs.index)).str.lower()
    obs["condition"] = normalize_string_series(pd.Series(condition, index=obs.index))
    obs["original_cell_annotation"] = normalize_string_series(original_label)
    obs["original_cell_id"] = normalize_string_series(obs["cell_id"] if "cell_id" in obs.columns else pd.Series(obs.index.astype(str), index=obs.index))
    obs["barcodes"] = barcode_source.copy()
    obs["barcode"] = barcode_source.copy()
    obs["barcode_core"] = obs["barcodes"].map(normalize_barcode_core)
    obs["cell_id"] = ensure_unique_obs_names(obs["original_cell_id"], obs["sample_id"])
    obs["obs_name"] = obs["cell_id"].copy()
    return obs


def harmonize_embedded_gdt_tcr(obs: pd.DataFrame) -> pd.DataFrame:
    out = obs.copy()
    for col in TRA_TRB_COLS:
        out[col] = "" if not col.endswith(("_umis", "_reads")) else 0
    for col in TRG_TRD_COLS:
        out[col] = "" if not col.endswith(("_umis", "_reads")) else 0

    mapping = {
        "TRG_v": "v_gene_TRG",
        "TRG_d": "d_gene_TRG",
        "TRG_j": "j_gene_TRG",
        "TRG_c_gene": "c_gene_TRG",
        "TRG_cdr3": "cdr3_TRG",
        "TRG_cdr3_nt": "cdr3_nt_TRG",
        "TRD_v": "v_gene_TRD",
        "TRD_d": "d_gene_TRD",
        "TRD_j": "j_gene_TRD",
        "TRD_c_gene": "c_gene_TRD",
        "TRD_cdr3": "cdr3_TRD",
        "TRD_cdr3_nt": "cdr3_nt_TRD",
    }
    for target, source in mapping.items():
        if source in out.columns:
            out[target] = normalize_string_series(out[source]).replace({"None": "", "Multi": ""})

    trg_chain = normalize_chain_label(out.get("chain_TRG", pd.Series("", index=out.index)))
    trd_chain = normalize_chain_label(out.get("chain_TRD", pd.Series("", index=out.index)))
    valid_trg = trg_chain.eq("TRG") & out["TRG_cdr3"].ne("") & out["TRG_cdr3_nt"].ne("")
    valid_trd = trd_chain.eq("TRD") & out["TRD_cdr3"].ne("") & out["TRD_cdr3_nt"].ne("")

    for col in ["TRG_v", "TRG_d", "TRG_j", "TRG_cdr3", "TRG_cdr3_nt", "TRG_c_gene"]:
        out.loc[~valid_trg, col] = ""
    for col in ["TRD_v", "TRD_d", "TRD_j", "TRD_cdr3", "TRD_cdr3_nt", "TRD_c_gene"]:
        out.loc[~valid_trd, col] = ""

    paired_gdt = out["TRG_cdr3"].ne("") & out["TRD_cdr3"].ne("")
    any_gdt = out["TRG_cdr3"].ne("") | out["TRD_cdr3"].ne("")
    out["has_TRA_TRB_paired"] = False
    out["has_any_ab_tcr"] = False
    out["has_TRG_TRD_paired"] = paired_gdt.to_numpy(dtype=bool, copy=False)
    out["has_any_gd_tcr"] = any_gdt.to_numpy(dtype=bool, copy=False)
    out["TCRseq"] = np.where(any_gdt, "yes", "no")
    return out


def build_anndata(counts: sparse.csr_matrix, var_names: pd.Index, obs: pd.DataFrame) -> ad.AnnData:
    var = pd.DataFrame(index=var_names.astype(str))
    adata = ad.AnnData(X=counts, obs=obs.copy(), var=var)
    adata.var_names_make_unique()
    adata.obs.index = pd.Index(obs["cell_id"].astype(str), name=None)
    adata.obs.index.name = None
    return adata


def sanitize_obs_for_h5ad(adata: ad.AnnData) -> None:
    for column in adata.obs.columns:
        series = adata.obs[column]
        if pd.api.types.is_bool_dtype(series):
            continue
        if pd.api.types.is_numeric_dtype(series):
            continue
        adata.obs[column] = normalize_string_series(series).astype(object)


def compute_scanpy_umap(adata: ad.AnnData, dataset: DatasetConfig) -> None:
    logging.info("Computing Scanpy UMAP for %s", dataset.dataset_id)
    work = adata.copy()
    batch_key = "sample_id" if work.obs["sample_id"].astype(str).nunique() > 1 else None
    n_top_genes = int(min(4000, max(1000, work.n_vars - 1)))
    try:
        sc.pp.highly_variable_genes(
            work,
            flavor="seurat_v3",
            n_top_genes=n_top_genes,
            batch_key=batch_key,
            subset=False,
        )
    except Exception as exc:  # pragma: no cover - fallback guard
        logging.warning("seurat_v3 HVG failed for %s (%s); falling back to cell_ranger", dataset.dataset_id, exc)
        try:
            sc.pp.highly_variable_genes(
                work,
                flavor="cell_ranger",
                n_top_genes=n_top_genes,
                subset=False,
            )
        except Exception as exc2:  # pragma: no cover - fallback guard
            logging.warning(
                "cell_ranger HVG also failed for %s (%s); using all genes for UMAP",
                dataset.dataset_id,
                exc2,
            )
            work.var["highly_variable"] = True
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    hvg_mask = work.var.get("highly_variable", pd.Series(True, index=work.var_names)).to_numpy(dtype=bool, copy=False)
    if int(hvg_mask.sum()) < 200:
        hvg_mask = np.ones(work.n_vars, dtype=bool)
    work = work[:, hvg_mask].copy()
    n_comps = int(min(50, max(2, work.n_obs - 1), max(2, work.n_vars - 1)))
    sc.pp.pca(work, n_comps=n_comps, svd_solver="arpack")
    sc.pp.neighbors(work, n_neighbors=min(15, max(5, work.n_obs - 1)), n_pcs=min(30, n_comps))
    sc.tl.umap(work, random_state=RANDOM_STATE)
    sc.tl.leiden(work, key_added="leiden", flavor="igraph", n_iterations=2)
    adata.obsm["X_umap"] = work.obsm["X_umap"].astype(np.float32, copy=False)
    adata.obsm["X_pca"] = work.obsm["X_pca"].astype(np.float32, copy=False)
    adata.obs["leiden"] = work.obs["leiden"].astype(str).to_numpy()


def write_h5ad_atomic(adata: ad.AnnData, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.parent / f".{out_path.name}.{os.getpid()}.tmp"
    if tmp_path.exists():
        tmp_path.unlink()
    sanitize_obs_for_h5ad(adata)
    adata.write_h5ad(tmp_path, convert_strings_to_categoricals=False)
    os.replace(tmp_path, out_path)


def run_phase4_scoring(dataset: DatasetConfig) -> None:
    run_command(
        [
            str(PYTHON_BIN),
            str(PHASE4_SCORER),
            "--input-h5ad",
            str(dataset.output_h5ad),
            "--alias",
            dataset.dataset_id,
            "--write-back",
        ],
        cwd=PROJECT_ROOT,
    )


def plot_raw_umap_scores(adata: ad.AnnData, figure_dir: Path) -> Path:
    if "X_umap" not in adata.obsm:
        raise KeyError("Missing X_umap for UMAP score plot.")
    umap = np.asarray(adata.obsm["X_umap"])
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for ax, column, title, cmap in [
        (axes[0], "phase4_trd_score", "Raw TRD score", "magma"),
        (axes[1], "phase4_trab_score", "Raw TRAB score", "viridis"),
    ]:
        scatter = ax.scatter(
            umap[:, 0],
            umap[:, 1],
            c=adata.obs[column].to_numpy(dtype=np.float32, copy=True),
            cmap=cmap,
            s=3,
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(title)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
    out_path = figure_dir / "phase4_umap_raw_trd_trab_scores.png"
    fig.savefig(out_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_trab_trd_paired_gdt(adata: ad.AnnData, figure_dir: Path) -> Path:
    paired = adata.obs["has_TRG_TRD_paired"].to_numpy(dtype=bool, copy=False)
    trab = adata.obs["phase4_trab_score"].to_numpy(dtype=np.float32, copy=True)
    trd = adata.obs["phase4_trd_score"].to_numpy(dtype=np.float32, copy=True)
    fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
    ax.scatter(trab[~paired], trd[~paired], s=3, c="#2F6690", linewidths=0, rasterized=True, alpha=0.65, label="Other")
    ax.scatter(trab[paired], trd[paired], s=3, c="#C1121F", linewidths=0, rasterized=True, alpha=0.75, label="Paired TRG/TRD")
    ax.set_xlabel("Raw TRAB score")
    ax.set_ylabel("Raw TRD score")
    ax.set_title(f"Paired TRG/TRD (n={int(paired.sum()):,})")
    ax.legend(loc="best", frameon=True)
    out_path = figure_dir / "phase4_trab_vs_trd_paired_trg_trd.png"
    fig.savefig(out_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    return out_path


def write_dataset_qc(dataset: DatasetConfig, adata: ad.AnnData, figure_dir: Path) -> tuple[dict[str, object], Path]:
    table_dir = TABLE_ROOT / dataset.dataset_id
    table_dir.mkdir(parents=True, exist_ok=True)
    LOG_ROOT.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    summary = {
        "dataset_id": dataset.dataset_id,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_samples": int(adata.obs["sample_id"].astype(str).nunique()),
        "n_donors": int(adata.obs["donor_id"].astype(str).replace("", np.nan).dropna().nunique()),
        "n_tissues": int(adata.obs["tissue"].astype(str).replace("", np.nan).dropna().nunique()),
        "cells_with_any_gd_tcr": int(adata.obs["has_any_gd_tcr"].sum()),
        "cells_with_paired_gd_tcr": int(adata.obs["has_TRG_TRD_paired"].sum()),
        "phase4_trd_score_mean": float(adata.obs["phase4_trd_score"].mean()),
        "phase4_trab_score_mean": float(adata.obs["phase4_trab_score"].mean()),
    }
    pd.DataFrame([summary]).to_csv(table_dir / f"{dataset.dataset_id}_summary.csv", index=False)
    tissue_counts = adata.obs["tissue"].astype(str).value_counts().rename_axis("tissue").reset_index(name="n_cells")
    sample_counts = adata.obs["sample_id"].astype(str).value_counts().rename_axis("sample_id").reset_index(name="n_cells")
    tissue_counts.to_csv(table_dir / f"{dataset.dataset_id}_tissue_counts.csv", index=False)
    sample_counts.to_csv(table_dir / f"{dataset.dataset_id}_sample_counts.csv", index=False)

    qc_path = LOG_ROOT / f"{dataset.dataset_id}_intake.md"
    lines = [
        f"# {dataset.dataset_id} Sorted gdT Intake Summary",
        "",
        f"- Input RDS: `{dataset.input_rds}`",
        f"- Output H5AD: `{dataset.output_h5ad}`",
        f"- Cells: `{summary['n_cells']:,}`",
        f"- Genes: `{summary['n_genes']:,}`",
        f"- Samples: `{summary['n_samples']}`",
        f"- Donors: `{summary['n_donors']}`",
        f"- Tissues: `{summary['n_tissues']}`",
        f"- Any gdTCR: `{summary['cells_with_any_gd_tcr']:,}`",
        f"- Paired TRG/TRD: `{summary['cells_with_paired_gd_tcr']:,}`",
        f"- Mean raw TRD score: `{summary['phase4_trd_score_mean']:.4f}`",
        f"- Mean raw TRAB score: `{summary['phase4_trab_score_mean']:.4f}`",
        "",
        "## Notes",
        "",
        f"- Dataset description: `{dataset.description}`",
        f"- Sample-key rule used: `{dataset.sample_rule}`",
        "- Scanpy UMAP was recomputed from the RNA counts; Seurat reductions were not used as the canonical output.",
        "- Embedded γδTCR metadata was harmonized into canonical `TRG_*` / `TRD_*` fields.",
        "",
        "## Figures",
        "",
        f"- `{figure_dir / 'phase4_umap_raw_trd_trab_scores.png'}`",
        f"- `{figure_dir / 'phase4_trab_vs_trd_paired_trg_trd.png'}`",
    ]
    qc_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return summary, qc_path


def process_one_dataset(dataset: DatasetConfig) -> dict[str, object]:
    logging.info("Processing sorted gdT dataset %s", dataset.dataset_id)
    with tempfile.TemporaryDirectory(prefix=f"{dataset.dataset_id}_payload_") as tmpdir:
        payload_prefix = export_rds_payload(dataset, Path(tmpdir))
        counts, var_names, meta = read_payload(payload_prefix)

    obs = build_dataset_specific_obs(meta, dataset)
    obs = harmonize_embedded_gdt_tcr(obs)

    adata = build_anndata(counts, var_names, obs)
    adata.uns["sorted_gdt_intake"] = {
        "dataset_id": dataset.dataset_id,
        "input_rds": str(dataset.input_rds.resolve()),
        "sample_rule": dataset.sample_rule,
        "description": dataset.description,
        "sorted_gdt": True,
    }
    compute_scanpy_umap(adata, dataset)
    write_h5ad_atomic(adata, dataset.output_h5ad)

    run_phase4_scoring(dataset)

    scored = ad.read_h5ad(dataset.output_h5ad)
    figure_dir = FIGURE_ROOT / f"{dataset.dataset_id}_phase4"
    plot_raw_umap_scores(scored, figure_dir)
    plot_trab_trd_paired_gdt(scored, figure_dir)
    summary, qc_path = write_dataset_qc(dataset, scored, figure_dir)
    summary["qc_markdown"] = str(qc_path)
    summary["output_h5ad"] = str(dataset.output_h5ad)
    return summary


def write_combined_summary(rows: list[dict[str, object]]) -> None:
    TABLE_ROOT.mkdir(parents=True, exist_ok=True)
    LOG_ROOT.mkdir(parents=True, exist_ok=True)
    summary_df = pd.DataFrame(rows).sort_values("dataset_id")
    summary_df.to_csv(TABLE_ROOT / "sorted_gdt_dataset_summary.csv", index=False)

    lines = [
        "# Sorted gdT Intake Summary",
        "",
        f"- Datasets processed: `{len(summary_df)}`",
        "",
    ]
    for _, row in summary_df.iterrows():
        lines.extend(
            [
                f"## {row['dataset_id']}",
                "",
                f"- Cells: `{int(row['n_cells']):,}`",
                f"- Samples: `{int(row['n_samples'])}`",
                f"- Donors: `{int(row['n_donors'])}`",
                f"- Tissues: `{int(row['n_tissues'])}`",
                f"- Paired TRG/TRD: `{int(row['cells_with_paired_gd_tcr']):,}`",
                f"- Output H5AD: `{row['output_h5ad']}`",
                f"- QC note: `{row['qc_markdown']}`",
                "",
            ]
        )
    (LOG_ROOT / "sorted_gdt_dataset_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    configure_logging()
    selected = set(args.dataset_id)
    datasets = [dataset for dataset in DATASETS if not selected or dataset.dataset_id in selected]
    if not datasets:
        raise ValueError("No datasets selected for sorted gdT intake.")

    rows = [process_one_dataset(dataset) for dataset in datasets]
    write_combined_summary(rows)
    logging.info("Processed %s sorted gdT datasets", len(rows))


if __name__ == "__main__":
    main()
