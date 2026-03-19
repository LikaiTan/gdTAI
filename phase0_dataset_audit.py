#!/usr/bin/env python3
"""Run the Phase 0 dataset audit for the T/NK integration workflow."""

from __future__ import annotations

import concurrent.futures
import math
import os
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


# Config
PROJECT_ROOT = Path(__file__).resolve().parent
REGISTRY_CSV = PROJECT_ROOT / "h5ad.csv"
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
FIGURE_DIR = OUTPUT_ROOT / "figures"
TABLE_DIR = OUTPUT_ROOT / "tables"
LOG_DIR = OUTPUT_ROOT / "logs"

PHASE0_AUDIT_CSV = TABLE_DIR / "phase0_dataset_audit.csv"
PHASE0_CATEGORY_CSV = TABLE_DIR / "phase0_category_summary.csv"
PHASE0_QC_MD = LOG_DIR / "phase0_qc_summary.md"

TARGET_TCR_GENES = [
    "TRAC",
    "TRBC1",
    "TRBC2",
    "TRDC",
    "TRGC1",
    "TRGC2",
    "CD3D",
    "CD3E",
    "CD3G",
    "NKG7",
    "KLRD1",
]
METADATA_HINTS = {
    "donor": ["donor", "patient", "subject", "individual"],
    "sample": ["sample", "orig.ident", "specimen", "biosample"],
    "library": ["library", "library_id", "batch", "lane", "channel"],
}
COUNTS_LAYER_PRIORITY = ["counts", "raw_counts", "count", "umi_counts"]


def normalize_attr(value):
    """Normalize HDF5 attrs so downstream code sees plain Python types."""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.ndarray) and value.dtype.kind == "S":
        return [v.decode("utf-8", errors="replace") for v in value.tolist()]
    return value


def sample_dense_dataset(dataset: h5py.Dataset, max_items: int = 50_000) -> np.ndarray:
    """Sample a dense dataset without reading the full matrix into memory."""
    if dataset.ndim == 0:
        return np.asarray([dataset[()]])
    if dataset.ndim == 1:
        total = dataset.shape[0]
        if total == 0:
            return np.asarray([])
        n = min(total, max_items)
        return np.asarray(dataset[:n])

    rows = min(dataset.shape[0], 256)
    cols = min(dataset.shape[1], 256)
    block = np.asarray(dataset[:rows, :cols]).ravel()
    return block[:max_items]


def sample_node_values(node) -> np.ndarray:
    """Sample numeric values from a dense or sparse matrix node."""
    if isinstance(node, h5py.Dataset):
        return sample_dense_dataset(node)

    if isinstance(node, h5py.Group) and "data" in node:
        data = node["data"]
        total = data.shape[0]
        if total == 0:
            return np.asarray([])
        n = min(total, 50_000)
        return np.asarray(data[:n])

    return np.asarray([])


def read_string_dataset(dataset: h5py.Dataset) -> pd.Index:
    """Read an HDF5 string dataset into a pandas Index."""
    data = dataset.asstr()[...] if hasattr(dataset, "asstr") else dataset[...]
    return pd.Index(np.asarray(data).astype(str))


def dataframe_index_name(group: h5py.Group) -> str:
    """Return the stored index dataset name for an H5AD obs/var dataframe group."""
    index_name = normalize_attr(group.attrs.get("_index", "_index"))
    if isinstance(index_name, list):
        return str(index_name[0])
    return str(index_name)


def inspect_obs_var_groups(handle: h5py.File, result: dict) -> None:
    """Read obs/var metadata directly from the H5AD file."""
    if "obs" in handle:
        obs_group = handle["obs"]
        obs_index = dataframe_index_name(obs_group)
        obs_columns = [key for key in obs_group.keys() if key != obs_index]
        metadata = detect_metadata_fields(obs_columns)
        result["donor_fields"] = metadata["donor"]
        result["sample_fields"] = metadata["sample"]
        result["library_fields"] = metadata["library"]
        result["obs_columns_n"] = len(obs_columns)
        if math.isnan(result["n_obs"]) and obs_index in obs_group:
            result["n_obs"] = int(obs_group[obs_index].shape[0])

    if "var" in handle:
        var_group = handle["var"]
        var_index = dataframe_index_name(var_group)
        if var_index in var_group:
            var_names = read_string_dataset(var_group[var_index])
            result["duplicate_var_names"] = int(var_names.duplicated().sum())
            result["gene_naming_style"] = detect_gene_style(var_names)
            result["n_vars"] = int(len(var_names)) if math.isnan(result["n_vars"]) else result["n_vars"]
            upper = pd.Index(var_names.str.upper())
            present = [gene for gene in TARGET_TCR_GENES if gene in upper]
            result["tcr_genes_present"] = ";".join(present)


def matrix_shape(node) -> tuple[int | float, int | float]:
    """Return a best-effort matrix shape from an H5AD X/layer node."""
    if isinstance(node, h5py.Dataset):
        shape = tuple(node.shape)
    elif isinstance(node, h5py.Group):
        shape = normalize_attr(node.attrs.get("shape", (np.nan, np.nan)))
        shape = tuple(shape) if isinstance(shape, (list, tuple, np.ndarray)) else (np.nan, np.nan)
    else:
        shape = (np.nan, np.nan)

    if len(shape) == 1:
        return int(shape[0]), 1
    if len(shape) >= 2:
        return int(shape[0]), int(shape[1])
    return np.nan, np.nan


def matrix_storage(node) -> str:
    """Describe how the matrix is stored in the H5AD file."""
    if isinstance(node, h5py.Dataset):
        return f"dense:{node.dtype}"
    if isinstance(node, h5py.Group):
        encoding = normalize_attr(node.attrs.get("encoding-type", "group"))
        return str(encoding)
    return "missing"


def matrix_nnz(node):
    """Return exact nnz for sparse matrices when available."""
    if isinstance(node, h5py.Group) and "data" in node:
        return int(node["data"].shape[0])
    return np.nan


def integer_like_fraction(values: np.ndarray) -> float:
    """Estimate how much a numeric sample behaves like raw counts."""
    if values.size == 0:
        return np.nan
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    return float(np.mean(np.isclose(values, np.round(values), atol=1e-6)))


def has_negative(values: np.ndarray) -> bool | float:
    """Detect whether a numeric sample contains negative values."""
    if values.size == 0:
        return np.nan
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    return bool(np.any(values < 0))


def max_sample(values: np.ndarray) -> float:
    """Return the maximum finite sampled value."""
    if values.size == 0:
        return np.nan
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    return float(np.max(values))


def find_counts_like_key(layer_keys: list[str]) -> str:
    """Pick the most likely counts-like layer from an H5AD layers group."""
    lowered = {key.lower(): key for key in layer_keys}
    for key in COUNTS_LAYER_PRIORITY:
        if key in lowered:
            return lowered[key]
    for key in layer_keys:
        if "count" in key.lower():
            return key
    return ""


def classify_state(
    x_negative,
    x_integer,
    has_counts_key,
    counts_integer,
    has_raw,
    raw_integer,
) -> tuple[str, str, str]:
    """Assign a matrix-state label and Phase 0 triage category."""
    if x_negative is True:
        return (
            "scaled_or_integrated_X",
            "C",
            "adata.X contains negative values, consistent with scaled or integrated data",
        )
    if not math.isnan(x_integer) and x_integer >= 0.98:
        return (
            "raw_like_X",
            "A",
            "adata.X is strongly integer-like and suitable as raw-count input",
        )
    if has_counts_key and not math.isnan(counts_integer) and counts_integer >= 0.98:
        return (
            "recoverable_counts_layer",
            "B",
            "A counts-like layer is present and integer-like",
        )
    if has_raw and not math.isnan(raw_integer) and raw_integer >= 0.98:
        return (
            "recoverable_raw_slot",
            "B",
            "adata.raw is present and appears integer-like",
        )
    if x_negative is False:
        return (
            "normalized_or_noninteger_X",
            "C",
            "adata.X is non-integer without a recoverable raw-count source",
        )
    return (
        "unreadable_or_unknown",
        "C",
        "Unable to classify a reliable raw-count source",
    )


def detect_gene_style(var_names: pd.Index) -> str:
    """Classify whether gene names look like symbols, Ensembl IDs, or mixed."""
    if len(var_names) == 0:
        return "empty"

    sample = pd.Index(var_names[: min(len(var_names), 2000)]).astype(str)
    upper = sample.str.upper()
    ensg_frac = float(np.mean(upper.str.startswith("ENSG")))
    symbol_like_frac = float(np.mean(upper.str.fullmatch(r"[A-Z0-9_.-]+", na=False)))

    if ensg_frac >= 0.7:
        return "ensembl_like"
    if symbol_like_frac >= 0.7:
        return "gene_symbol_like"
    return "mixed_or_custom"


def detect_metadata_fields(columns: list[str]) -> dict[str, str]:
    """Report candidate donor/sample/library columns from obs metadata."""
    lowered = {col.lower(): col for col in columns}
    out = {}
    for label, hints in METADATA_HINTS.items():
        matches = []
        for key, original in lowered.items():
            if any(hint in key for hint in hints):
                matches.append(original)
        out[label] = ";".join(sorted(dict.fromkeys(matches)))
    return out


def inspect_dataset(row: pd.Series) -> dict:
    """Inspect one H5AD file and return a flat audit record."""
    result = {
        "gse_id": row["gse_id"],
        "h5ad_path": row["h5ad_path"],
        "source_root": row["source_root"],
        "n_obs": np.nan,
        "n_vars": np.nan,
        "x_storage": "missing",
        "x_nnz": np.nan,
        "x_density": np.nan,
        "x_has_negative": np.nan,
        "x_integer_like_fraction": np.nan,
        "x_max_sample": np.nan,
        "has_counts_layer": False,
        "counts_layer_key": "",
        "counts_integer_like_fraction": np.nan,
        "has_raw": False,
        "raw_integer_like_fraction": np.nan,
        "state_label": "unreadable_or_unknown",
        "phase0_category": "C",
        "category_reason": "File not yet inspected",
        "duplicate_var_names": np.nan,
        "gene_naming_style": "unknown",
        "tcr_genes_present": "",
        "donor_fields": "",
        "sample_fields": "",
        "library_fields": "",
        "obs_columns_n": np.nan,
        "read_error": "",
    }

    path = Path(row["h5ad_path"])
    if not path.exists():
        result["read_error"] = "missing_file"
        result["state_label"] = "missing_file"
        result["category_reason"] = "Registry path does not exist"
        return result

    try:
        with h5py.File(path, "r") as handle:
            if "X" in handle:
                x_node = handle["X"]
                result["x_storage"] = matrix_storage(x_node)
                result["x_nnz"] = matrix_nnz(x_node)
                x_sample = sample_node_values(x_node)
                result["x_has_negative"] = has_negative(x_sample)
                result["x_integer_like_fraction"] = integer_like_fraction(x_sample)
                result["x_max_sample"] = max_sample(x_sample)

                n_obs, n_vars = matrix_shape(x_node)
                if not math.isnan(n_obs):
                    result["n_obs"] = int(n_obs)
                if not math.isnan(n_vars):
                    result["n_vars"] = int(n_vars)
                if not math.isnan(result["x_nnz"]) and result["n_obs"] and result["n_vars"]:
                    denom = int(result["n_obs"]) * int(result["n_vars"])
                    result["x_density"] = float(result["x_nnz"] / denom) if denom else np.nan

            layer_keys = list(handle["layers"].keys()) if "layers" in handle else []
            counts_key = find_counts_like_key(layer_keys)
            if counts_key:
                result["has_counts_layer"] = True
                result["counts_layer_key"] = counts_key
                counts_sample = sample_node_values(handle["layers"][counts_key])
                result["counts_integer_like_fraction"] = integer_like_fraction(counts_sample)

            if "raw" in handle and "X" in handle["raw"]:
                result["has_raw"] = True
                raw_sample = sample_node_values(handle["raw"]["X"])
                result["raw_integer_like_fraction"] = integer_like_fraction(raw_sample)

            inspect_obs_var_groups(handle, result)
    except Exception as exc:
        result["read_error"] = f"h5py:{type(exc).__name__}:{exc}"

    state, category, reason = classify_state(
        result["x_has_negative"],
        result["x_integer_like_fraction"],
        bool(result["counts_layer_key"]),
        result["counts_integer_like_fraction"],
        bool(result["has_raw"]),
        result["raw_integer_like_fraction"],
    )
    result["state_label"] = state
    result["phase0_category"] = category
    result["category_reason"] = reason

    if result["read_error"] and state == "unreadable_or_unknown":
        result["category_reason"] = "File could not be fully read for audit"

    return result


def configure_plotting() -> None:
    """Set publication-readable plotting defaults."""
    sns.set_theme(style="whitegrid")
    plt.rcParams.update(
        {
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 9,
        }
    )


def write_figures(audit: pd.DataFrame) -> None:
    """Generate the required Phase 0 PNG figures."""
    configure_plotting()

    fig1, ax1 = plt.subplots(figsize=(7, 4.5))
    cat_counts = audit["phase0_category"].value_counts().reindex(["A", "B", "C"], fill_value=0)
    cat_df = pd.DataFrame({"phase0_category": cat_counts.index, "count": cat_counts.values})
    sns.barplot(
        data=cat_df,
        x="phase0_category",
        y="count",
        hue="phase0_category",
        palette=["#2e8b57", "#d28b00", "#b22222"],
        legend=False,
        ax=ax1,
    )
    ax1.set_title("Phase 0 Dataset Categories")
    ax1.set_xlabel("Category")
    ax1.set_ylabel("Datasets")
    for idx, value in enumerate(cat_counts.values):
        ax1.text(idx, value + 0.1, str(int(value)), ha="center", va="bottom", fontsize=10)
    fig1.tight_layout()
    fig1.savefig(FIGURE_DIR / "phase0_category_distribution.png")
    plt.close(fig1)

    fig2, ax2 = plt.subplots(figsize=(10, 5.5))
    state_counts = audit["state_label"].value_counts().reset_index()
    state_counts.columns = ["state_label", "count"]
    sns.barplot(data=state_counts, x="state_label", y="count", color="#4c78a8", ax=ax2)
    ax2.set_title("Observed Matrix States Across Registered H5AD Files")
    ax2.set_xlabel("State label")
    ax2.set_ylabel("Datasets")
    ax2.tick_params(axis="x", rotation=25)
    fig2.tight_layout()
    fig2.savefig(FIGURE_DIR / "phase0_matrix_state_overview.png")
    plt.close(fig2)

    plot_df = audit.dropna(subset=["n_obs", "n_vars"]).copy()
    plot_df["log10_n_obs"] = np.log10(plot_df["n_obs"].clip(lower=1))
    plot_df["log10_n_vars"] = np.log10(plot_df["n_vars"].clip(lower=1))
    fig3, ax3 = plt.subplots(figsize=(8, 6))
    palette = {"A": "#2e8b57", "B": "#d28b00", "C": "#b22222"}
    sns.scatterplot(
        data=plot_df,
        x="log10_n_obs",
        y="log10_n_vars",
        hue="phase0_category",
        palette=palette,
        s=70,
        ax=ax3,
    )
    for _, rec in plot_df.iterrows():
        ax3.text(rec["log10_n_obs"] + 0.01, rec["log10_n_vars"] + 0.005, rec["gse_id"], fontsize=7)
    ax3.set_title("Registered Dataset Size Overview")
    ax3.set_xlabel("log10(number of cells)")
    ax3.set_ylabel("log10(number of genes)")
    fig3.tight_layout()
    fig3.savefig(FIGURE_DIR / "phase0_dataset_size_overview.png")
    plt.close(fig3)

    meta_rows = []
    for field in ["donor", "sample", "library"]:
        with_field = int((audit[f"{field}_fields"].fillna("") != "").sum())
        without_field = int((audit[f"{field}_fields"].fillna("") == "").sum())
        meta_rows.append(
            {
                "field_type": field,
                "datasets_with_field": with_field,
                "datasets_without_field": without_field,
            }
        )
    meta_df = pd.DataFrame(meta_rows)
    meta_long = meta_df.melt(id_vars="field_type", var_name="status", value_name="count")
    fig4, ax4 = plt.subplots(figsize=(7.5, 4.8))
    sns.barplot(
        data=meta_long,
        x="field_type",
        y="count",
        hue="status",
        palette=["#4c78a8", "#d3d3d3"],
        ax=ax4,
    )
    ax4.set_title("Metadata Field Coverage Across Datasets")
    ax4.set_xlabel("Metadata type")
    ax4.set_ylabel("Datasets")
    fig4.tight_layout()
    fig4.savefig(FIGURE_DIR / "phase0_metadata_completeness.png")
    plt.close(fig4)


def write_qc_summary(audit: pd.DataFrame) -> None:
    """Write a compact Phase 0 QC summary for user review."""
    cat_counts = audit["phase0_category"].value_counts().reindex(["A", "B", "C"], fill_value=0)
    state_counts = audit["state_label"].value_counts()
    readable = int((audit["read_error"].fillna("") == "").sum())
    unreadable = int((audit["read_error"].fillna("") != "").sum())
    category_b = audit.loc[audit["phase0_category"] == "B", ["gse_id", "state_label", "category_reason"]]
    category_c = audit.loc[audit["phase0_category"] == "C", ["gse_id", "state_label", "category_reason"]]

    with PHASE0_QC_MD.open("w", encoding="utf-8") as handle:
        handle.write("# Phase 0 QC Summary\n\n")
        handle.write("## Environment\n\n")
        handle.write("- Requested env from runbook: `Scanpy_gdTmodel`\n")
        handle.write(
            "- Resolved working env for this run: `/home/tanlikai/miniconda3/envs/rapids_sc_py310`\n"
        )
        handle.write("- GPU detected by `nvidia-smi`: NVIDIA A100 80GB PCIe, CUDA 12.9\n\n")

        handle.write("## Audit Scope\n\n")
        handle.write(f"- Registry rows audited: {len(audit)}\n")
        handle.write(f"- Files readable without recorded errors: {readable}\n")
        handle.write(f"- Files with read issues: {unreadable}\n\n")

        handle.write("## Category Counts\n\n")
        for category in ["A", "B", "C"]:
            handle.write(f"- Category {category}: {int(cat_counts.get(category, 0))}\n")

        handle.write("\n## Matrix State Counts\n\n")
        for label, count in state_counts.items():
            handle.write(f"- {label}: {int(count)}\n")

        handle.write("\n## Metadata Coverage\n\n")
        for field in ["donor", "sample", "library"]:
            count = int((audit[f"{field}_fields"].fillna("") != "").sum())
            handle.write(f"- Datasets with candidate {field} fields: {count}\n")

        handle.write("\n## Category B Datasets\n\n")
        if category_b.empty:
            handle.write("- None\n")
        else:
            for _, row in category_b.iterrows():
                handle.write(f"- {row['gse_id']}: {row['state_label']} | {row['category_reason']}\n")

        handle.write("\n## Category C Datasets\n\n")
        if category_c.empty:
            handle.write("- None\n")
        else:
            for _, row in category_c.iterrows():
                handle.write(f"- {row['gse_id']}: {row['state_label']} | {row['category_reason']}\n")

        handle.write("\n## QC Conclusion\n\n")
        handle.write(
            "Phase 0 audit outputs are generated and ready for user review. "
            "Do not advance to Phase 1 until the user reviews the tables and figures "
            "and explicitly approves the phase transition.\n"
        )


def ensure_output_dirs() -> None:
    """Create the canonical output structure if it does not already exist."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)


def load_registry() -> pd.DataFrame:
    """Load and validate the H5AD registry."""
    registry = pd.read_csv(REGISTRY_CSV)
    required_cols = {"h5ad_path", "gse_id", "source_root"}
    missing = required_cols - set(registry.columns)
    if missing:
        raise ValueError(f"Missing required columns in h5ad.csv: {sorted(missing)}")
    return registry


def run_phase0_audit() -> pd.DataFrame:
    """Run the dataset audit over every registry row."""
    registry = load_registry()
    workers = min(8, len(registry), max(1, os.cpu_count() or 1))
    records = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as pool:
        future_map = {
            pool.submit(inspect_dataset, row): row["gse_id"]
            for _, row in registry.iterrows()
        }
        for future in concurrent.futures.as_completed(future_map):
            gse_id = future_map[future]
            record = future.result()
            records.append(record)
            print(f"Audited {gse_id} ...", flush=True)

    audit = pd.DataFrame(records)
    audit = audit.sort_values(["phase0_category", "gse_id"]).reset_index(drop=True)
    audit.to_csv(PHASE0_AUDIT_CSV, index=False)

    summary = (
        audit.groupby(["phase0_category", "state_label"], dropna=False)
        .size()
        .reset_index(name="dataset_count")
        .sort_values(
            ["phase0_category", "dataset_count", "state_label"],
            ascending=[True, False, True],
        )
    )
    summary.to_csv(PHASE0_CATEGORY_CSV, index=False)
    write_figures(audit)
    write_qc_summary(audit)
    return audit


def main() -> None:
    """Entrypoint for the Phase 0 audit bootstrap."""
    ensure_output_dirs()
    audit = run_phase0_audit()
    print("\nWrote:")
    print(PHASE0_AUDIT_CSV)
    print(PHASE0_CATEGORY_CSV)
    print(PHASE0_QC_MD)
    for fig_name in [
        "phase0_category_distribution.png",
        "phase0_matrix_state_overview.png",
        "phase0_dataset_size_overview.png",
        "phase0_metadata_completeness.png",
    ]:
        print(FIGURE_DIR / fig_name)
    print(f"\nAudited {len(audit)} registry rows.")


if __name__ == "__main__":
    main()
