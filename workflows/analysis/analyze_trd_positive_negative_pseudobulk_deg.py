#!/usr/bin/env python3
"""Paired pseudobulk DEG for TRD-positive vs TRD-negative T cells.

The comparison uses the current repaired integrated milestone and applies the
following cell filter before pseudobulk construction:
- `scanvi_tnk_superclass == T_cell`
- paired `TRA/TRB`
- `phase4_trab_score > -0.05`
- subgroup: `phase4_trd_score > 0` versus `phase4_trd_score < 0`

Samples with fewer than 50 cells in either subgroup are excluded so the DEG is
paired at the `sampleid` level.
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


from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import false_discovery_control, wilcoxon


PROJECT_ROOT = Path(__file__).resolve().parents[2]
INTEGRATED_H5AD = PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "integrated.h5ad"
TABLE_DIR = PROJECT_ROOT / "Integrated_dataset" / "tables"
LOG_DIR = PROJECT_ROOT / "Integrated_dataset" / "logs"
FIGURE_DIR = PROJECT_ROOT / "Integrated_dataset" / "figures"

MIN_CELLS_PER_GROUP = 50
CHUNK_ROWS = 20_000
TOP_N = 50
HEATMAP_DPI = 300

COUNTS_BY_SAMPLE = TABLE_DIR / "trd_positive_negative_pseudobulk_sample_counts.csv"
DEG_TABLE = TABLE_DIR / "trd_positive_negative_pseudobulk_deg.csv"
TOP_TABLE = TABLE_DIR / "trd_positive_negative_pseudobulk_top50_up_down.csv"
HEATMAP_PATH = FIGURE_DIR / "trd_positive_negative_pseudobulk_top50_up_down_heatmap.png"
SUMMARY_MD = LOG_DIR / "trd_positive_negative_pseudobulk_deg_summary.md"


def read_string_array(group: h5py.Group, key: str) -> np.ndarray:
    """Read one H5AD string-like obs/var element."""
    obj = group[key]
    if isinstance(obj, h5py.Group) and "codes" in obj and "categories" in obj:
        codes = obj["codes"][:]
        categories = obj["categories"].asstr()[:]
        out = np.empty(len(codes), dtype=object)
        neg_mask = codes < 0
        out[neg_mask] = ""
        out[~neg_mask] = categories[codes[~neg_mask]]
        return out.astype(str)
    return obj.asstr()[:].astype(str)


def build_filtered_obs() -> tuple[pd.DataFrame, np.ndarray]:
    """Build the filtered cell-level metadata frame and keep-mask."""
    with h5py.File(INTEGRATED_H5AD, "r") as handle:
        obs = handle["obs"]
        frame = pd.DataFrame(
            {
                "sampleid": read_string_array(obs, "sampleid"),
                "source_gse_id": read_string_array(obs, "source_gse_id"),
                "scanvi_tnk_superclass": read_string_array(obs, "scanvi_tnk_superclass"),
                "TRA_cdr3": read_string_array(obs, "TRA_cdr3"),
                "TRB_cdr3": read_string_array(obs, "TRB_cdr3"),
                "phase4_trab_score": np.asarray(obs["phase4_trab_score"][:], dtype=np.float32),
                "phase4_trd_score": np.asarray(obs["phase4_trd_score"][:], dtype=np.float32),
            }
        )

    paired = frame["TRA_cdr3"].str.strip().ne("") & frame["TRB_cdr3"].str.strip().ne("")
    base = (
        frame["scanvi_tnk_superclass"].eq("T_cell")
        & paired
        & (frame["phase4_trab_score"] > -0.05)
        & frame["sampleid"].astype(str).str.strip().ne("")
    )
    frame = frame.loc[base].copy()
    frame["trd_group"] = np.where(
        frame["phase4_trd_score"] > 0,
        "TRD>0",
        np.where(frame["phase4_trd_score"] < 0, "TRD<0", "TRD=0"),
    )
    frame = frame.loc[frame["trd_group"].isin(["TRD>0", "TRD<0"])].copy()
    return frame, base.to_numpy()


def choose_paired_samples(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Keep only samples with enough cells in both subgroups."""
    counts = (
        frame.groupby(["sampleid", "source_gse_id", "trd_group"], observed=True)
        .size()
        .rename("cell_n")
        .reset_index()
    )
    pivot = (
        counts.pivot_table(
            index=["sampleid", "source_gse_id"],
            columns="trd_group",
            values="cell_n",
            fill_value=0,
        )
        .reset_index()
    )
    for column in ["TRD>0", "TRD<0"]:
        if column not in pivot.columns:
            pivot[column] = 0
    eligible = pivot.loc[(pivot["TRD>0"] >= MIN_CELLS_PER_GROUP) & (pivot["TRD<0"] >= MIN_CELLS_PER_GROUP)].copy()
    frame = frame.merge(eligible[["sampleid", "source_gse_id"]], on=["sampleid", "source_gse_id"], how="inner")
    return frame, eligible


def read_var_names() -> np.ndarray:
    """Read gene names from the integrated H5AD."""
    with h5py.File(INTEGRATED_H5AD, "r") as handle:
        var = handle["var"]
        if "_index" in var:
            return read_string_array(var, "_index")
        if "gene_symbols" in var:
            return read_string_array(var, "gene_symbols")
        raise KeyError("No usable gene-name field found in `var`.")


def build_group_index(frame: pd.DataFrame) -> tuple[np.ndarray, pd.DataFrame]:
    """Map each kept cell to a sample-pseudobulk group row."""
    sample_order = sorted(frame["sampleid"].astype(str).unique())
    group_rows = []
    for sampleid in sample_order:
        sample_df = frame.loc[frame["sampleid"] == sampleid]
        gse_id = sample_df["source_gse_id"].iloc[0]
        for trd_group in ["TRD>0", "TRD<0"]:
            group_rows.append(
                {
                    "sampleid": sampleid,
                    "source_gse_id": gse_id,
                    "trd_group": trd_group,
                    "group_key": f"{sampleid}__{trd_group}",
                }
            )
    group_df = pd.DataFrame(group_rows)
    group_to_idx = {row.group_key: idx for idx, row in enumerate(group_df.itertuples(index=False))}
    cell_group_index = (frame["sampleid"].astype(str) + "__" + frame["trd_group"].astype(str)).map(group_to_idx).to_numpy(dtype=np.int32)
    return cell_group_index, group_df


def build_pseudobulk_counts(frame: pd.DataFrame, var_names: np.ndarray) -> tuple[np.ndarray, pd.DataFrame]:
    """Construct sample-level pseudobulk count matrices from the sparse H5AD."""
    cell_indices = frame.index.to_numpy(dtype=np.int64, copy=False)
    cell_group_index, group_df = build_group_index(frame)
    n_groups = len(group_df)
    n_vars = len(var_names)
    pseudobulk = np.zeros((n_groups, n_vars), dtype=np.float64)

    with h5py.File(INTEGRATED_H5AD, "r") as handle:
        x_group = handle["X"]
        indptr = x_group["indptr"][:]
        data = x_group["data"]
        indices = x_group["indices"]

        for start in range(0, len(cell_indices), CHUNK_ROWS):
            stop = min(start + CHUNK_ROWS, len(cell_indices))
            chunk_rows = cell_indices[start:stop]
            chunk_groups = cell_group_index[start:stop]
            order = np.argsort(chunk_rows)
            sorted_rows = chunk_rows[order]
            sorted_groups = chunk_groups[order]

            chunk_indptr = np.empty(len(sorted_rows) + 1, dtype=np.int64)
            chunk_indptr[0] = 0
            chunk_data_parts = []
            chunk_index_parts = []
            nnz_total = 0
            for i, row_idx in enumerate(sorted_rows):
                row_start = indptr[row_idx]
                row_end = indptr[row_idx + 1]
                chunk_data_parts.append(data[row_start:row_end])
                chunk_index_parts.append(indices[row_start:row_end])
                nnz_total += row_end - row_start
                chunk_indptr[i + 1] = nnz_total

            if nnz_total == 0:
                continue
            chunk_matrix = sparse.csr_matrix(
                (
                    np.concatenate(chunk_data_parts).astype(np.float64, copy=False),
                    np.concatenate(chunk_index_parts).astype(np.int64, copy=False),
                    chunk_indptr,
                ),
                shape=(len(sorted_rows), n_vars),
            )
            design = sparse.csr_matrix(
                (
                    np.ones(len(sorted_groups), dtype=np.float64),
                    (np.arange(len(sorted_groups)), sorted_groups),
                ),
                shape=(len(sorted_groups), n_groups),
            )
            pseudobulk += (design.T @ chunk_matrix).toarray()

    return pseudobulk, group_df


def log_cpm(counts: np.ndarray) -> np.ndarray:
    """Convert pseudobulk counts to log1p CPM."""
    library_size = counts.sum(axis=1, keepdims=True)
    library_size[library_size == 0] = 1.0
    cpm = counts / library_size * 1_000_000.0
    return np.log1p(cpm)


def run_paired_deg(log_expr: np.ndarray, group_df: pd.DataFrame, var_names: np.ndarray) -> pd.DataFrame:
    """Run paired Wilcoxon DEG on sample-matched pseudobulks."""
    pos_idx = group_df.index[group_df["trd_group"] == "TRD>0"].to_numpy()
    neg_idx = group_df.index[group_df["trd_group"] == "TRD<0"].to_numpy()
    pos = log_expr[pos_idx]
    neg = log_expr[neg_idx]
    diff = pos - neg

    pvals = np.ones(diff.shape[1], dtype=np.float64)
    for gene_idx in range(diff.shape[1]):
        vector = diff[:, gene_idx]
        if np.allclose(vector, 0):
            pvals[gene_idx] = 1.0
            continue
        try:
            pvals[gene_idx] = wilcoxon(vector, zero_method="wilcox", alternative="two-sided", method="auto").pvalue
        except ValueError:
            pvals[gene_idx] = 1.0
    padj = false_discovery_control(pvals, method="bh")

    return pd.DataFrame(
        {
            "gene": var_names,
            "mean_log1p_cpm_trd_gt_0": pos.mean(axis=0),
            "mean_log1p_cpm_trd_lt_0": neg.mean(axis=0),
            "mean_diff_trd_gt_0_minus_lt_0": diff.mean(axis=0),
            "median_diff_trd_gt_0_minus_lt_0": np.median(diff, axis=0),
            "wilcoxon_pvalue": pvals,
            "wilcoxon_fdr": padj,
            "samples_n": pos.shape[0],
        }
    ).sort_values(["wilcoxon_fdr", "mean_diff_trd_gt_0_minus_lt_0"], ascending=[True, False]).reset_index(drop=True)


def choose_top_genes(deg_df: pd.DataFrame) -> pd.DataFrame:
    """Choose the top 50 up and top 50 down paired DEGs."""
    up = deg_df.loc[deg_df["mean_diff_trd_gt_0_minus_lt_0"] > 0].sort_values(
        ["wilcoxon_fdr", "mean_diff_trd_gt_0_minus_lt_0"],
        ascending=[True, False],
    ).head(TOP_N)
    down = deg_df.loc[deg_df["mean_diff_trd_gt_0_minus_lt_0"] < 0].sort_values(
        ["wilcoxon_fdr", "mean_diff_trd_gt_0_minus_lt_0"],
        ascending=[True, True],
    ).head(TOP_N)
    top = pd.concat([up.assign(direction="up_in_TRD_gt_0"), down.assign(direction="down_in_TRD_gt_0")], ignore_index=True)
    return top


def write_heatmap(log_expr: np.ndarray, group_df: pd.DataFrame, top_df: pd.DataFrame, var_names: np.ndarray) -> None:
    """Write a paired pseudobulk heatmap for the top up/down genes."""
    gene_to_idx = {gene: idx for idx, gene in enumerate(var_names)}
    sample_order = sorted(group_df["sampleid"].unique())

    ordered_group_idx = []
    ordered_labels = []
    for sampleid in sample_order:
        for trd_group in ["TRD>0", "TRD<0"]:
            idx = group_df.index[(group_df["sampleid"] == sampleid) & (group_df["trd_group"] == trd_group)][0]
            ordered_group_idx.append(idx)
            ordered_labels.append(f"{sampleid}\n{trd_group}")

    gene_order = top_df.loc[top_df["direction"] == "up_in_TRD_gt_0", "gene"].tolist() + top_df.loc[top_df["direction"] == "down_in_TRD_gt_0", "gene"].tolist()
    matrix = log_expr[np.asarray(ordered_group_idx)][:, [gene_to_idx[gene] for gene in gene_order]].T
    matrix = matrix.astype(np.float64, copy=False)
    matrix = (matrix - matrix.mean(axis=1, keepdims=True)) / np.maximum(matrix.std(axis=1, keepdims=True), 1e-6)

    fig_w = max(14, len(ordered_labels) * 0.16)
    fig_h = max(14, len(gene_order) * 0.12)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(matrix, aspect="auto", cmap="coolwarm", vmin=-2.5, vmax=2.5, interpolation="nearest")
    ax.set_xticks(np.arange(len(ordered_labels)))
    ax.set_xticklabels(ordered_labels, rotation=90, fontsize=6)
    ax.set_yticks(np.arange(len(gene_order)))
    ax.set_yticklabels(gene_order, fontsize=7)
    ax.set_title("Paired pseudobulk DEG heatmap: top 50 up and top 50 down genes", fontsize=14)
    ax.set_xlabel("Sample pseudobulks", fontsize=10)
    ax.set_ylabel("Genes", fontsize=10)
    plt.colorbar(im, ax=ax, fraction=0.018, pad=0.01, label="Gene-wise z-score of log1p CPM")
    fig.tight_layout()
    fig.savefig(HEATMAP_PATH, dpi=HEATMAP_DPI, bbox_inches="tight")
    plt.close(fig)


def write_summary_markdown(eligible: pd.DataFrame, top_df: pd.DataFrame) -> None:
    """Write the DEG summary markdown."""
    up_preview = ", ".join(top_df.loc[top_df["direction"] == "up_in_TRD_gt_0", "gene"].head(10))
    down_preview = ", ".join(top_df.loc[top_df["direction"] == "down_in_TRD_gt_0", "gene"].head(10))
    lines = [
        "# Paired Pseudobulk DEG: TRD Positive vs TRD Negative T Cells",
        "",
        "- Cell filter: `scanvi_tnk_superclass == T_cell`, paired `TRA/TRB`, `phase4_trab_score > -0.05`",
        "- Comparison: `phase4_trd_score > 0` versus `phase4_trd_score < 0`",
        f"- Minimum cells per sample per group: `{MIN_CELLS_PER_GROUP}`",
        f"- Eligible paired samples: `{len(eligible)}`",
        f"- Top up genes preview: `{up_preview}`",
        f"- Top down genes preview: `{down_preview}`",
        "",
        "Outputs:",
        f"- `{COUNTS_BY_SAMPLE.relative_to(PROJECT_ROOT)}`",
        f"- `{DEG_TABLE.relative_to(PROJECT_ROOT)}`",
        f"- `{TOP_TABLE.relative_to(PROJECT_ROOT)}`",
        f"- `{HEATMAP_PATH.relative_to(PROJECT_ROOT)}`",
    ]
    SUMMARY_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Run the paired pseudobulk DEG workflow."""
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    frame, _ = build_filtered_obs()
    frame, eligible = choose_paired_samples(frame)
    eligible.to_csv(COUNTS_BY_SAMPLE, index=False)

    var_names = read_var_names()
    pseudobulk_counts, group_df = build_pseudobulk_counts(frame, var_names)
    log_expr = log_cpm(pseudobulk_counts)
    deg_df = run_paired_deg(log_expr, group_df, var_names)
    deg_df.to_csv(DEG_TABLE, index=False)

    top_df = choose_top_genes(deg_df)
    top_df.to_csv(TOP_TABLE, index=False)
    write_heatmap(log_expr, group_df, top_df, var_names)
    write_summary_markdown(eligible, top_df)
    print(top_df.head(20).to_csv(index=False))


if __name__ == "__main__":
    main()
