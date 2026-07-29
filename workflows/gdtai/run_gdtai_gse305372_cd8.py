#!/usr/bin/env python3
"""Apply promoted gdTAI v3 to author-labeled GSE305372 CD8 T cells."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import html
import json
import os
import pickle
import sys
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/gdtai_gse305372_matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import mmread

PROJECT_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_DIR = PROJECT_ROOT / "workflows" / "gdtai"
if str(WORKFLOW_DIR) not in sys.path:
    sys.path.insert(0, str(WORKFLOW_DIR))

from run_gdtai_v3_trdc_nk_guard_classifier import (  # noqa: E402
    FeatureSpec,
    append_engineered_features,
    trdc_trdv_quadrant,
)

DATASET_ID = "GSE305372"
DATASET_ROOT = PROJECT_ROOT / "data" / "datasets" / DATASET_ID
DEFAULT_PAYLOAD_DIR = DATASET_ROOT / "interim" / "rds_export"
DEFAULT_MODEL = (
    PROJECT_ROOT
    / "Integrated_dataset"
    / "models"
    / "gdT_prediction_classifier"
    / "gdTAI_v3.0"
    / "gdTAI_v3_model.pkl"
)
DEFAULT_OUTPUT_DIR = (
    PROJECT_ROOT / "Integrated_dataset" / "tables" / "gdT_prediction" / DATASET_ID
)
DEFAULT_FIGURE_DIR = (
    PROJECT_ROOT / "Integrated_dataset" / "figures" / "gdT_prediction" / DATASET_ID
)
DEFAULT_LOG_DIR = (
    PROJECT_ROOT / "Integrated_dataset" / "logs" / "gdT_prediction" / DATASET_ID
)
DEFAULT_REPORT_DIR = PROJECT_ROOT / "reports" / "GSE305372_gdtai_cd8"
MODEL_REGISTRY = PROJECT_ROOT / "configs" / "models" / "gdtai" / "model_registry.csv"

COMPARTMENTS = {
    "lung": {
        "prefix": "LG_CD8",
        "cluster_column": "RNA_snn_res.0.2",
        "cluster_labels": {
            "0": "TRM",
            "1": "TEM",
            "2": "GZMKhi",
            "3": "TCM",
            "4": "NKG2Cpos",
            "5": "MAIT",
            "6": "Prolif",
        },
    },
    "lymph_node": {
        "prefix": "LN_CD8",
        "cluster_column": "RNA_snn_res.0.3",
        "cluster_labels": {
            "0": "GZMKhi",
            "1": "TI",
            "2": "TCM",
            "3": "TN",
            "4": "IFNR",
            "5": "TRM",
            "6": "TACT",
            "7": "MAIT",
            "8": "TEM",
        },
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--payload-dir", type=Path, default=DEFAULT_PAYLOAD_DIR)
    parser.add_argument("--model-pkl", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--log-dir", type=Path, default=DEFAULT_LOG_DIR)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--chunk-size", type=int, default=50_000)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_present(values: pd.Series) -> np.ndarray:
    strings = values.astype("string").fillna("").str.strip().str.lower()
    return (~strings.isin(["", "na", "nan", "none", "null"])).to_numpy(dtype=bool)


def build_spec(payload: dict[str, Any]) -> FeatureSpec:
    genes = [str(value) for value in payload["gene_names"]]
    engineered = [str(value) for value in payload["engineered_feature_names"]]
    model_features = [str(value) for value in payload["feature_names"]]
    gene_features = [f"{gene}_log1p_cp10k" for gene in genes]
    expected = [*gene_features, *engineered]
    if expected != model_features:
        raise RuntimeError("Packaged model feature order is internally inconsistent")
    return FeatureSpec(
        gene_names=genes,
        gene_indices=np.arange(len(genes), dtype=np.int32),
        gene_feature_names=gene_features,
        engineered_feature_names=engineered,
        model_feature_names=model_features,
        gene_to_col={gene: index for index, gene in enumerate(genes)},
        engineered_to_col={
            name: len(genes) + index for index, name in enumerate(engineered)
        },
    )


def add_author_fields(metadata: pd.DataFrame, compartment: str) -> pd.DataFrame:
    config = COMPARTMENTS[compartment]
    cluster_column = str(config["cluster_column"])
    if cluster_column not in metadata:
        raise KeyError(f"Author cluster column is missing: {cluster_column}")
    raw_cluster = metadata[cluster_column].astype("string")
    numeric_cluster = pd.to_numeric(raw_cluster, errors="coerce")
    metadata["author_cluster_id"] = raw_cluster.fillna("unassigned").astype(str)
    integral = numeric_cluster.notna() & (numeric_cluster % 1 == 0)
    metadata.loc[integral, "author_cluster_id"] = numeric_cluster.loc[
        integral
    ].astype(int).astype(str)
    metadata["author_cluster"] = metadata["author_cluster_id"].map(
        config["cluster_labels"]
    ).fillna("unassigned")
    metadata["tissue"] = "lung" if compartment == "lung" else "lung-associated lymph node"

    has_tra = np.zeros(len(metadata), dtype=bool)
    has_trb = np.zeros(len(metadata), dtype=bool)
    for column in ("cdr3a.aa.seq", "cdr3a.nt.seq", "tra.v"):
        if column in metadata:
            has_tra |= is_present(metadata[column])
    for column in ("cdr3b.aa.seq", "cdr3b.nt.seq", "trb.v"):
        if column in metadata:
            has_trb |= is_present(metadata[column])
    metadata["has_TRA_evidence"] = has_tra
    metadata["has_TRB_evidence"] = has_trb
    metadata["has_paired_TRA_TRB_evidence"] = has_tra & has_trb
    return metadata


def read_payload(
    payload_dir: Path,
    compartment: str,
    spec: FeatureSpec,
    model: Any,
    threshold: float,
    chunk_size: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    prefix = payload_dir / str(COMPARTMENTS[compartment]["prefix"])
    matrix_path = Path(f"{prefix}_model_gene_counts.mtx.gz")
    genes_path = Path(f"{prefix}_model_genes.tsv")
    metadata_path = Path(f"{prefix}_cell_metadata.csv.gz")
    manifest_path = Path(f"{prefix}_export_manifest.csv")
    for path in (matrix_path, genes_path, metadata_path, manifest_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    genes = genes_path.read_text(encoding="utf-8").splitlines()
    if genes != spec.gene_names:
        raise RuntimeError(f"Model-gene order mismatch in {genes_path}")
    with gzip.open(matrix_path, "rb") as handle:
        gene_by_cell = mmread(handle).tocsr().astype(np.float32)
    metadata = pd.read_csv(metadata_path, low_memory=False)
    if gene_by_cell.shape != (len(spec.gene_names), len(metadata)):
        raise RuntimeError(
            f"Payload shape mismatch for {compartment}: {gene_by_cell.shape} vs "
            f"{len(spec.gene_names)} genes x {len(metadata)} cells"
        )
    if metadata["cell_id"].duplicated().any():
        raise RuntimeError(f"Duplicated cell IDs in {metadata_path}")

    totals = pd.to_numeric(metadata["row_sum_counts_layer"], errors="raise").to_numpy(
        dtype=np.float32
    )
    if np.any(totals <= 0):
        raise RuntimeError(f"Non-positive raw-count library sizes in {metadata_path}")

    scores = np.empty(len(metadata), dtype=np.float32)
    quadrants = np.empty(len(metadata), dtype=object)
    trdc = np.empty(len(metadata), dtype=np.float32)
    trdv = np.empty(len(metadata), dtype=np.float32)
    trg = np.empty(len(metadata), dtype=np.float32)
    abt = np.empty(len(metadata), dtype=np.float32)
    for start in range(0, len(metadata), chunk_size):
        end = min(start + chunk_size, len(metadata))
        values = gene_by_cell[:, start:end].T.toarray().astype(np.float32, copy=False)
        with np.errstate(divide="ignore", invalid="ignore"):
            x_gene = np.log1p(values * (10_000.0 / totals[start:end, None]))
        x_gene[~np.isfinite(x_gene)] = 0.0
        x = append_engineered_features(x_gene, spec)
        scores[start:end] = model.predict_proba(x)[:, 1].astype(np.float32)
        quadrants[start:end] = trdc_trdv_quadrant(x, spec)
        trdc[start:end] = x[:, spec.engineered_to_col["TRDC_log1p"]]
        trdv[start:end] = x[:, spec.engineered_to_col["TRDV_score"]]
        trg[start:end] = x[:, spec.engineered_to_col["TRG_score"]]
        abt[start:end] = x[:, spec.engineered_to_col["abT_TCR_score"]]

    metadata = add_author_fields(metadata, compartment)
    metadata["gdtai_model_id"] = "gdtai_v3"
    metadata["gdtai_v3_score"] = scores
    metadata["gdtai_v3_threshold"] = threshold
    metadata["gdtai_v3_predicted"] = scores >= threshold
    metadata["tcr_gene_quadrant"] = quadrants
    metadata["TRDC_log1p_cp10k"] = trdc
    metadata["TRDV_sum_log1p_cp10k"] = trdv
    metadata["TRG_sum_log1p_cp10k"] = trg
    metadata["TRAB_sum_log1p_cp10k"] = abt

    manifest = pd.read_csv(manifest_path).iloc[0].to_dict()
    manifest.update(
        payload_matrix=str(matrix_path.relative_to(PROJECT_ROOT)),
        payload_cells=int(len(metadata)),
        payload_matrix_nnz=int(gene_by_cell.nnz),
        payload_matrix_sha256=sha256_file(matrix_path),
        payload_genes_sha256=sha256_file(genes_path),
        payload_metadata_sha256=sha256_file(metadata_path),
        payload_export_manifest_sha256=sha256_file(manifest_path),
    )
    return metadata, manifest


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total == 0:
        return float("nan"), float("nan")
    proportion = successes / total
    denominator = 1 + z * z / total
    center = (proportion + z * z / (2 * total)) / denominator
    margin = z * np.sqrt(
        proportion * (1 - proportion) / total + z * z / (4 * total * total)
    ) / denominator
    return max(0.0, center - margin), min(1.0, center + margin)


def summarize_groups(data: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    grouped = data.groupby(columns, dropna=False, observed=True, sort=True)
    for key, frame in grouped:
        keys = key if isinstance(key, tuple) else (key,)
        n_cells = len(frame)
        n_predicted = int(frame["gdtai_v3_predicted"].sum())
        low, high = wilson_interval(n_predicted, n_cells)
        row = dict(zip(columns, keys))
        row.update(
            n_cells=n_cells,
            n_predicted=n_predicted,
            predicted_fraction=n_predicted / n_cells if n_cells else np.nan,
            predicted_fraction_ci_low=low,
            predicted_fraction_ci_high=high,
            median_score=float(frame["gdtai_v3_score"].median()),
            p95_score=float(frame["gdtai_v3_score"].quantile(0.95)),
        )
        rows.append(row)
    return pd.DataFrame(rows)


def make_tables(data: pd.DataFrame, output_dir: Path) -> dict[str, pd.DataFrame]:
    total = summarize_groups(data.assign(scope="all_CD8A"), ["scope"])
    by_tissue = summarize_groups(data, ["source_compartment", "tissue"])
    overall = pd.concat([total, by_tissue], ignore_index=True, sort=False)
    by_cluster = summarize_groups(
        data, ["source_compartment", "author_cluster_id", "author_cluster"]
    )
    by_donor = summarize_groups(data, ["source_compartment", "donor.id.tag"])
    by_quadrant = summarize_groups(data, ["source_compartment", "tcr_gene_quadrant"])
    by_tcr = summarize_groups(
        data,
        [
            "source_compartment",
            "has_TRA_evidence",
            "has_TRB_evidence",
            "has_paired_TRA_TRB_evidence",
        ],
    )
    tables = {
        "overall": overall,
        "by_cluster": by_cluster,
        "by_donor": by_donor,
        "by_quadrant": by_quadrant,
        "by_tcr_evidence": by_tcr,
    }
    for name, table in tables.items():
        table.to_csv(output_dir / f"gdtai_v3_{name}.csv", index=False)
    return tables


def plot_umap(data: pd.DataFrame, figure_dir: Path) -> Path:
    output = figure_dir / "gdtai_v3_cd8_umap.png"
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4), constrained_layout=True)
    for axis, compartment in zip(axes, COMPARTMENTS):
        frame = data[data["source_compartment"] == compartment]
        predicted = frame["gdtai_v3_predicted"].to_numpy(dtype=bool)
        axis.scatter(
            frame.loc[~predicted, "UMAP_1"],
            frame.loc[~predicted, "UMAP_2"],
            s=1.2,
            c="#C8CDD3",
            alpha=0.35,
            linewidths=0,
            rasterized=True,
        )
        axis.scatter(
            frame.loc[predicted, "UMAP_1"],
            frame.loc[predicted, "UMAP_2"],
            s=5,
            c="#C43C39",
            alpha=0.85,
            linewidths=0,
            rasterized=True,
        )
        axis.set_title(
            f"{compartment.replace('_', ' ').title()}\n"
            f"{predicted.sum():,} / {len(frame):,} predicted"
        )
        axis.set_xlabel("UMAP 1")
        axis.set_ylabel("UMAP 2")
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle("GSE305372 author-labeled CD8 T cells: gdTAI v3 predictions", fontsize=14)
    fig.savefig(output, dpi=300, facecolor="white")
    plt.close(fig)
    return output


def plot_scores(data: pd.DataFrame, figure_dir: Path, threshold: float) -> Path:
    output = figure_dir / "gdtai_v3_cd8_score_distributions.png"
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True, sharey=True)
    for axis, compartment in zip(axes, COMPARTMENTS):
        values = data.loc[
            data["source_compartment"] == compartment, "gdtai_v3_score"
        ].to_numpy(dtype=float)
        weights = np.full(values.size, 100.0 / values.size)
        axis.hist(
            values,
            bins=np.linspace(0, 1, 81),
            weights=weights,
            color="#3B728F",
            alpha=0.9,
        )
        axis.axvline(threshold, color="#C43C39", linewidth=2, label=f"threshold {threshold:.3f}")
        axis.set_title(compartment.replace("_", " ").title())
        axis.set_xlabel("gdTAI v3 score")
        axis.set_ylabel("Cells per bin (%)")
        axis.set_yscale("log")
        axis.legend(frameon=False)
        axis.spines[["top", "right"]].set_visible(False)
    fig.savefig(output, dpi=300, facecolor="white")
    plt.close(fig)
    return output


def plot_cluster_rates(by_cluster: pd.DataFrame, figure_dir: Path) -> Path:
    output = figure_dir / "gdtai_v3_cd8_prediction_by_author_cluster.png"
    frame = by_cluster.copy()
    frame["label"] = (
        frame["source_compartment"].str.replace("_", " ").str.title()
        + " | "
        + frame["author_cluster"].astype(str)
    )
    frame = frame.sort_values(
        ["source_compartment", "predicted_fraction", "author_cluster"],
        ascending=[True, True, True],
    )
    y = np.arange(len(frame))
    fig, axis = plt.subplots(figsize=(9, max(4.5, 0.38 * len(frame))), constrained_layout=True)
    colors = np.where(frame["source_compartment"] == "lung", "#287271", "#D17A22")
    axis.barh(y, frame["predicted_fraction"], color=colors, alpha=0.9)
    axis.set_yticks(y, frame["label"])
    axis.set_xlabel("Predicted fraction")
    axis.set_xlim(left=0)
    axis.spines[["top", "right"]].set_visible(False)
    for position, (_, row) in enumerate(frame.iterrows()):
        axis.text(
            row["predicted_fraction"],
            position,
            f"  {int(row['n_predicted']):,}/{int(row['n_cells']):,}",
            va="center",
            fontsize=8,
        )
    fig.savefig(output, dpi=300, facecolor="white")
    plt.close(fig)
    return output


def format_table(table: pd.DataFrame, columns: list[str], max_rows: int = 30) -> str:
    view = table.loc[:, [column for column in columns if column in table]].head(max_rows).copy()
    for column in view.columns:
        if "fraction" in column or column in {"median_score", "p95_score"}:
            view[column] = pd.to_numeric(view[column], errors="coerce").map(
                lambda value: "" if pd.isna(value) else f"{value:.4f}"
            )
    view = view.fillna("")
    return view.to_html(index=False, border=0, classes="data-table", escape=True)


def render_report(
    data: pd.DataFrame,
    tables: dict[str, pd.DataFrame],
    figures: list[Path],
    report_dir: Path,
    manifest: dict[str, Any],
) -> Path:
    report_dir.mkdir(parents=True, exist_ok=True)
    report_figure_dir = report_dir / "figures"
    report_figure_dir.mkdir(parents=True, exist_ok=True)
    for source in figures:
        target = report_figure_dir / source.name
        target.write_bytes(source.read_bytes())

    total = tables["overall"].iloc[0]
    n_predicted = int(total["n_predicted"])
    n_cells = int(total["n_cells"])
    paired_predicted = int(
        (
            data["gdtai_v3_predicted"]
            & data["has_paired_TRA_TRB_evidence"]
        ).sum()
    )
    paired_cells = int(data["has_paired_TRA_TRB_evidence"].sum())
    paired_call_fraction = paired_predicted / paired_cells if paired_cells else np.nan
    trdv_positive_predictions = int(
        (
            data["gdtai_v3_predicted"]
            & (data["TRDV_sum_log1p_cp10k"] > 0)
        ).sum()
    )
    nk_mait = data["author_cluster"].isin(["NKG2Cpos", "MAIT"])
    nk_mait_predicted = int((data["gdtai_v3_predicted"] & nk_mait).sum())
    report = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GSE305372 CD8 gdTAI v3 inference</title>
<style>
body {{ margin: 0; font-family: Arial, sans-serif; color: #20262c; background: #fff; line-height: 1.5; }}
header {{ background: #17324d; color: white; padding: 30px max(5vw, 28px); }}
main {{ max-width: 1160px; margin: 0 auto; padding: 28px; }}
h1 {{ margin: 0 0 8px; font-size: 30px; letter-spacing: 0; }}
h2 {{ margin-top: 34px; border-bottom: 2px solid #d8e0e7; padding-bottom: 7px; }}
.summary {{ display: grid; grid-template-columns: repeat(3, minmax(170px, 1fr)); gap: 12px; margin: 20px 0; }}
.metric {{ border: 1px solid #ccd6df; border-radius: 6px; padding: 14px; }}
.metric strong {{ display: block; font-size: 24px; color: #17324d; }}
.note {{ border-left: 4px solid #d17a22; background: #fff7ed; padding: 12px 16px; }}
figure {{ margin: 24px 0; }}
figure img {{ max-width: 100%; height: auto; border: 1px solid #d8e0e7; }}
.data-table {{ border-collapse: collapse; width: 100%; font-size: 13px; display: block; overflow-x: auto; }}
.data-table th, .data-table td {{ border: 1px solid #d8e0e7; padding: 6px 8px; text-align: left; white-space: nowrap; }}
.data-table th {{ background: #eef3f6; }}
code {{ background: #eef3f6; padding: 2px 4px; }}
@media (max-width: 760px) {{ .summary {{ grid-template-columns: 1fr 1fr; }} }}
@page {{ size: A4 landscape; margin: 10mm; }}
@media print {{
  body {{ font-size: 9pt; }}
  main {{ max-width: none; padding: 0; }}
  header {{ padding: 10mm; }}
  figure {{ break-inside: avoid; }}
  .data-table {{ display: table; width: 100%; font-size: 7pt; table-layout: auto; overflow: visible; }}
  .data-table thead {{ display: table-header-group; }}
  .data-table tr {{ break-inside: avoid; }}
  .data-table th, .data-table td {{ padding: 3px 4px; white-space: normal; overflow-wrap: anywhere; }}
}}
</style>
</head>
<body>
<header>
  <h1>gdTAI v3 in GSE305372 CD8 T cells</h1>
  <div>Independent application to author-annotated lung and lung-associated lymph-node CD8 cells</div>
</header>
<main>
<div class="summary">
  <div class="metric"><strong>{n_cells:,}</strong>CD8A cells analyzed</div>
  <div class="metric"><strong>{n_predicted:,}</strong>gdTAI-positive cells</div>
  <div class="metric"><strong>{n_predicted / n_cells:.3%}</strong>predicted fraction</div>
  <div class="metric"><strong>{manifest['threshold']:.3f}</strong>registered threshold</div>
  <div class="metric"><strong>{trdv_positive_predictions:,}</strong>predictions with TRDV expression</div>
  <div class="metric"><strong>{paired_call_fraction:.3%}</strong>paired TRA/TRB conflict-screening rate</div>
</div>

<div class="note"><strong>Interpretation limit.</strong> These downloaded CD8 processed objects contain TRA/TRB metadata but no TRG/TRD clonotype fields. Therefore, this is an inference report, not a precision/recall evaluation. Model-positive cells with paired TRA/TRB evidence ({paired_predicted:,} of {paired_cells:,}; {paired_call_fraction:.3%}) are apparent alpha-beta-associated calls requiring review; they cannot be labeled definitive false positives without matched gamma-delta TCR evidence.</div>

<h2>Input and method</h2>
<p>The analysis downloaded only the two GEO CD8 processed objects and retained cells with the authors' <code>cite.cell.type.tag == CD8A</code>. It did not select cells using TRD, TRG, gdTAI score, or any model feature. Raw RNA counts for the packaged 210 genes were normalized as <code>log1p(count * 10,000 / whole-transcriptome total counts)</code>. The promoted <code>gdtai_v3</code> Round 14 artifact was applied at its fixed threshold of <code>{manifest['threshold']:.3f}</code>. This CD8A restriction means CD8-negative gamma-delta T cells are outside the analysis scope.</p>
<p>Model SHA256: <code>{html.escape(manifest['model_sha256'])}</code></p>

<h2>Overall results</h2>
{format_table(tables['overall'], ['scope', 'source_compartment', 'tissue', 'n_cells', 'n_predicted', 'predicted_fraction', 'median_score', 'p95_score'])}

<figure><img src="figures/gdtai_v3_cd8_umap.png" alt="UMAP highlighting gdTAI predictions"></figure>
<figure><img src="figures/gdtai_v3_cd8_score_distributions.png" alt="gdTAI score distributions"></figure>

<h2>Author phenotype clusters</h2>
<p>{nk_mait_predicted:,} predictions fall in the authors' NKG2C-positive or MAIT clusters. These populations receive explicit attention because cytotoxic T cells, MAIT cells, gamma-delta T cells, and NK-like states can share expression programs.</p>
{format_table(tables['by_cluster'], ['source_compartment', 'author_cluster_id', 'author_cluster', 'n_cells', 'n_predicted', 'predicted_fraction', 'predicted_fraction_ci_low', 'predicted_fraction_ci_high'], 40)}
<figure><img src="figures/gdtai_v3_cd8_prediction_by_author_cluster.png" alt="Prediction fractions by author cluster"></figure>

<h2>Donor audit</h2>
<p>Calls are shown for every available author donor identifier so that donor concentration and small-denominator effects remain visible.</p>
{format_table(tables['by_donor'], ['source_compartment', 'donor.id.tag', 'n_cells', 'n_predicted', 'predicted_fraction', 'predicted_fraction_ci_low', 'predicted_fraction_ci_high'], 80)}

<h2>TCR-gene expression quadrants</h2>
{format_table(tables['by_quadrant'], ['source_compartment', 'tcr_gene_quadrant', 'n_cells', 'n_predicted', 'predicted_fraction', 'median_score'], 30)}

<h2>TRA/TRB sequencing evidence</h2>
{format_table(tables['by_tcr_evidence'], ['source_compartment', 'has_TRA_evidence', 'has_TRB_evidence', 'has_paired_TRA_TRB_evidence', 'n_cells', 'n_predicted', 'predicted_fraction'], 30)}

<h2>Reproducibility</h2>
<p>Source: <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305372">GEO GSE305372</a>. Study: <a href="https://www.nature.com/articles/s41590-026-02592-6">Fajardo-Rosas et al., Nature Immunology (2026)</a>.</p>
<p>Workflow: <code>workflows/intake/export_gse305372_cd8_model_payload.R</code> followed by <code>workflows/gdtai/run_gdtai_gse305372_cd8.py</code>.</p>
</main>
</body>
</html>
"""
    output = report_dir / "index.html"
    output.write_text(report, encoding="utf-8")
    return output


def write_markdown_summary(
    tables: dict[str, pd.DataFrame], log_dir: Path, manifest: dict[str, Any]
) -> Path:
    overall = tables["overall"]
    rows = [
        "# GSE305372 CD8 gdTAI v3 Summary",
        "",
        f"- model: `gdtai_v3`",
        f"- model SHA256: `{manifest['model_sha256']}`",
        f"- threshold: `{manifest['threshold']:.3f}`",
        "- input subset: author `cite.cell.type.tag == CD8A`",
        "- input state: raw RNA counts normalized to `log1p(CP10K)`",
        "- evaluation status: inference only; no gamma-delta TCR ground truth is available",
        "",
        "## Counts",
        "",
    ]
    for _, row in overall.iterrows():
        label = row.get("scope")
        if pd.isna(label) or str(label).strip() == "":
            label = row.get("source_compartment", "unknown")
        rows.append(
            f"- {label}: {int(row['n_predicted']):,} / {int(row['n_cells']):,} "
            f"({float(row['predicted_fraction']):.3%})"
        )
    output = log_dir / "gdtai_v3_gse305372_cd8_summary.md"
    output.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return output


def main() -> None:
    args = parse_args()
    for directory in (args.output_dir, args.figure_dir, args.log_dir, args.report_dir):
        directory.mkdir(parents=True, exist_ok=True)
    output_predictions = args.output_dir / "gdtai_v3_all_cd8_predictions.csv.gz"
    if output_predictions.exists() and not args.overwrite:
        raise FileExistsError(f"Output exists: {output_predictions}; pass --overwrite")

    registry = pd.read_csv(MODEL_REGISTRY)
    promoted = registry[
        (registry["model_id"] == "gdtai_v3")
        & (registry["status"] == "promoted_default")
    ]
    if len(promoted) != 1:
        raise RuntimeError("Model registry must contain one promoted gdtai_v3 row")
    expected_sha = str(promoted.iloc[0]["sha256"])
    observed_sha = sha256_file(args.model_pkl)
    if observed_sha != expected_sha:
        raise RuntimeError(
            f"Promoted model checksum mismatch: {observed_sha} != {expected_sha}"
        )

    with args.model_pkl.open("rb") as handle:
        payload = pickle.load(handle)
    model = payload["model_object"]
    threshold = float(payload["threshold"])
    spec = build_spec(payload)

    feature_manifest_path = args.model_pkl.parent / "feature_genes.csv"
    if not feature_manifest_path.is_file():
        raise FileNotFoundError(feature_manifest_path)
    feature_manifest = pd.read_csv(feature_manifest_path)
    registered_genes = feature_manifest.loc[
        feature_manifest["feature_type"] == "gene_log1p_cp10k", "gene"
    ].astype(str).tolist()
    if registered_genes != spec.gene_names:
        raise RuntimeError("Model pickle and released feature manifest disagree")

    frames: list[pd.DataFrame] = []
    source_manifests: list[dict[str, Any]] = []
    for compartment in COMPARTMENTS:
        frame, source_manifest = read_payload(
            args.payload_dir,
            compartment,
            spec,
            model,
            threshold,
            args.chunk_size,
        )
        frames.append(frame)
        source_manifests.append(source_manifest)
    data = pd.concat(frames, ignore_index=True, sort=False)
    if data["cell_id"].duplicated().any():
        raise RuntimeError("Combined prediction table contains duplicated cell IDs")

    expected_predictions = data["gdtai_v3_score"].to_numpy() >= threshold
    validation_checks = {
        "cell_ids_unique": bool(~data["cell_id"].duplicated().any()),
        "all_cells_author_tagged_CD8A": bool(
            data["author_cd8_subset_value"].eq("CD8A").all()
        ),
        "scores_finite": bool(np.isfinite(data["gdtai_v3_score"]).all()),
        "scores_within_zero_one": bool(
            data["gdtai_v3_score"].between(0.0, 1.0, inclusive="both").all()
        ),
        "prediction_mask_equals_fixed_threshold": bool(
            np.array_equal(
                data["gdtai_v3_predicted"].to_numpy(dtype=bool),
                expected_predictions,
            )
        ),
        "combined_cells_equal_payload_cells": bool(
            len(data) == sum(int(item["payload_cells"]) for item in source_manifests)
        ),
        "model_gene_count_is_210": len(spec.gene_names) == 210,
    }
    failed_checks = [name for name, passed in validation_checks.items() if not passed]
    if failed_checks:
        raise RuntimeError(f"Final validation failed: {', '.join(failed_checks)}")

    data.to_csv(output_predictions, index=False, compression="gzip")
    data.loc[data["gdtai_v3_predicted"]].to_csv(
        args.output_dir / "gdtai_v3_predicted_cd8_cells.csv.gz",
        index=False,
        compression="gzip",
    )
    tables = make_tables(data, args.output_dir)
    figures = [
        plot_umap(data, args.figure_dir),
        plot_scores(data, args.figure_dir, threshold),
        plot_cluster_rates(tables["by_cluster"], args.figure_dir),
    ]

    manifest = {
        "dataset_id": DATASET_ID,
        "analysis": "author_labeled_CD8A_gdtai_v3_inference",
        "model_id": "gdtai_v3",
        "model_status": "promoted_default",
        "model_path": str(args.model_pkl.relative_to(PROJECT_ROOT)),
        "model_sha256": observed_sha,
        "threshold": threshold,
        "normalization": "log1p(raw_count * 10000 / whole_transcriptome_total_count)",
        "model_gene_count": len(spec.gene_names),
        "total_cells": int(len(data)),
        "predicted_cells": int(data["gdtai_v3_predicted"].sum()),
        "predicted_fraction": float(data["gdtai_v3_predicted"].mean()),
        "paired_TRA_TRB_cells": int(data["has_paired_TRA_TRB_evidence"].sum()),
        "paired_TRA_TRB_predicted_cells": int(
            (
                data["has_paired_TRA_TRB_evidence"]
                & data["gdtai_v3_predicted"]
            ).sum()
        ),
        "paired_TRA_TRB_conflict_screening_fraction": float(
            data.loc[
                data["has_paired_TRA_TRB_evidence"], "gdtai_v3_predicted"
            ].mean()
        ),
        "predicted_cells_with_TRDV_expression": int(
            (
                data["gdtai_v3_predicted"]
                & (data["TRDV_sum_log1p_cp10k"] > 0)
            ).sum()
        ),
        "source_payloads": source_manifests,
        "validation_checks": validation_checks,
        "ground_truth_available": False,
        "tcr_metadata_limit": (
            "TRA/TRB fields available; TRG/TRD clonotype fields absent from "
            "the downloaded CD8 processed objects"
        ),
    }
    manifest_path = args.log_dir / "gdtai_v3_gse305372_cd8_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    validation_path = args.log_dir / "gdtai_v3_gse305372_cd8_validation.json"
    validation_path.write_text(
        json.dumps(validation_checks, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_markdown_summary(tables, args.log_dir, manifest)
    report = render_report(data, tables, figures, args.report_dir, manifest)
    print(
        json.dumps(
            {
                "cells": len(data),
                "predicted": int(data["gdtai_v3_predicted"].sum()),
                "predicted_fraction": float(data["gdtai_v3_predicted"].mean()),
                "threshold": threshold,
                "report": str(report),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
