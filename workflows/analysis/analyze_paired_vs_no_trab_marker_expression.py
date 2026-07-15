#!/usr/bin/env python3
"""Compare gdT-marker expression between paired and no-TRA/TRB cells.

This helper uses the current integrated milestone plus harmonized metadata to:
- define `Paired TRA/TRB` and `No TRA/TRB` groups from metadata, not scores
- extract log1p-normalized expression for `TRDC`, `TRGC1`, `TRDV1`, `TRDV2`
- run Wilcoxon rank-sum statistics with FDR correction
- write tables, one violin PNG, and a Markdown summary on NFS
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
import logging
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from phase4_gdt_module_scoring import TARGET_SUM, build_csr_chunk, normalize_log1p_chunk
from tissue_correction_workflow import build_metadata_key


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
FIGURE_DIR = OUTPUT_ROOT / "figures"
TABLE_DIR = OUTPUT_ROOT / "tables"
LOG_DIR = OUTPUT_ROOT / "logs"

INTEGRATED_H5AD = PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "integrated.h5ad"
METADATA_PATHS = [
    PROJECT_ROOT / "analysis_26GSE_V4" / "outputs" / "harmonized_metadata_v4.csv",
    PROJECT_ROOT / "analysis_26GSE_V4" / "outputs" / "harmonized_metadata_supp.csv",
]
GENES = ["TRDC", "TRGC1", "TRDV1", "TRDV2"]
CHUNK_SIZE = 50_000
PLOT_SAMPLE_PER_GROUP = 80_000
FIGURE_DPI = 300


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description="Compare marker expression between paired and no-TRA/TRB groups.")
    parser.add_argument(
        "--trd-min",
        type=float,
        default=None,
        help="Optional lower bound on raw phase4_trd_score before group comparison.",
    )
    parser.add_argument(
        "--output-prefix",
        type=str,
        default="paired_vs_no_trab_marker_expression",
        help="Prefix for output table, figure, and log filenames.",
    )
    return parser.parse_args()


def build_output_paths(prefix: str) -> dict[str, Path]:
    """Build one coherent output set for a given analysis prefix."""
    return {
        "stats_csv": TABLE_DIR / f"{prefix}_stats.csv",
        "gse_counts_csv": TABLE_DIR / f"{prefix}_group_counts_by_gse.csv",
        "sampleid_counts_csv": TABLE_DIR / f"{prefix}_group_counts_by_sampleid.csv",
        "plot_sample_csv": TABLE_DIR / f"{prefix}_plot_sample.csv.gz",
        "summary_md": LOG_DIR / f"{prefix}_summary.md",
        "figure_png": FIGURE_DIR / f"{prefix}_violin.png",
        "run_log": LOG_DIR / f"{prefix}.log",
    }


def configure_logging(run_log: Path) -> None:
    """Configure file and console logging."""
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(run_log, mode="a", encoding="utf-8"),
            logging.StreamHandler(),
        ],
        force=True,
    )


def ensure_output_dirs() -> None:
    """Create NFS-side output directories."""
    for path in [FIGURE_DIR, TABLE_DIR, LOG_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    """Apply Benjamini-Hochberg FDR correction."""
    p = np.asarray(p_values, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.empty(n, dtype=float)
    running = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        running = min(running, ranked[i] * n / rank)
        adjusted[i] = running
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adjusted, 0.0, 1.0)
    return out


def read_obs_strings(handle: h5py.File, column: str) -> np.ndarray:
    """Read one obs string column into a plain object array."""
    return np.asarray(
        [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in handle["obs"][column][:]],
        dtype=object,
    )


def load_group_lookup() -> pd.Series:
    """Load metadata-derived paired/no-TRA-TRB group labels keyed by metadata_key."""
    frames: list[pd.DataFrame] = []
    usecols = ["source_gse_id", "original_cell_id", "barcodes", "TRA_cdr3", "TRB_cdr3"]
    for path in METADATA_PATHS:
        if not path.exists():
            continue
        logging.info("Loading harmonized metadata: %s", path)
        frame = pd.read_csv(path, usecols=usecols, dtype=str, keep_default_na=False, low_memory=False)
        for column in usecols:
            frame[column] = frame[column].fillna("").astype(str)
        frame["metadata_key"] = [
            build_metadata_key(gse, original, barcode)
            for gse, original, barcode in zip(
                frame["source_gse_id"],
                frame["original_cell_id"],
                frame["barcodes"],
                strict=False,
            )
        ]
        frame["has_tra"] = frame["TRA_cdr3"].str.strip().ne("")
        frame["has_trb"] = frame["TRB_cdr3"].str.strip().ne("")
        frame["group_label"] = np.where(
            frame["has_tra"] & frame["has_trb"],
            "Paired TRA/TRB",
            np.where(~frame["has_tra"] & ~frame["has_trb"], "No TRA/TRB", "Single-chain only"),
        )
        frames.append(frame[["metadata_key", "group_label"]])
    if not frames:
        raise FileNotFoundError("No harmonized metadata files were found.")
    lookup = pd.concat(frames, ignore_index=True)
    lookup = lookup[lookup["metadata_key"].ne("")]
    lookup = lookup.drop_duplicates(subset=["metadata_key"], keep="first")
    lookup = lookup.set_index("metadata_key")["group_label"]
    if lookup.index.has_duplicates:
        raise ValueError("Metadata key collision detected in paired/no-TRA-TRB lookup.")
    return lookup


def load_obs_group_frame(lookup: pd.Series) -> pd.DataFrame:
    """Load join fields from the integrated milestone and assign group labels."""
    with h5py.File(INTEGRATED_H5AD, "r") as handle:
        obs_frame = pd.DataFrame(
            {
                "metadata_key": read_obs_strings(handle, "metadata_key"),
                "source_gse_id": read_obs_strings(handle, "source_gse_id"),
                "sampleid": read_obs_strings(handle, "sampleid"),
                "phase4_trd_score": handle["obs"]["phase4_trd_score"][:].astype(np.float32, copy=False),
            }
        )
    obs_frame["group_label"] = pd.Series(
        lookup.reindex(obs_frame["metadata_key"]).to_numpy(dtype=object),
        index=obs_frame.index,
        dtype="object",
    )
    obs_frame["group_label"] = obs_frame["group_label"].fillna("Single-chain only")
    allowed = {"Paired TRA/TRB", "No TRA/TRB", "Single-chain only"}
    unexpected = sorted(set(obs_frame["group_label"].dropna().unique()) - allowed)
    if unexpected:
        raise ValueError(f"Unexpected group labels after metadata join: {unexpected}")
    if not obs_frame["group_label"].eq("Paired TRA/TRB").any():
        raise ValueError("No `Paired TRA/TRB` cells were assigned.")
    if not obs_frame["group_label"].eq("No TRA/TRB").any():
        raise ValueError("No `No TRA/TRB` cells were assigned.")
    return obs_frame


def extract_gene_expression(keep_mask: np.ndarray | None = None) -> pd.DataFrame:
    """Extract all-cell log1p-normalized expression for the target genes."""
    logging.info("Extracting log1p-normalized expression for %s", ", ".join(GENES))
    with h5py.File(INTEGRATED_H5AD, "r") as handle:
        var_names = pd.Index(handle["var"]["_index"].asstr()[:], dtype="object")
        missing = [gene for gene in GENES if gene not in var_names]
        if missing:
            raise KeyError(f"Missing genes in integrated milestone: {missing}")
        gene_idx = np.asarray([int(var_names.get_loc(gene)) for gene in GENES], dtype=np.int32)
        n_obs = int(handle["obs"]["_index"].shape[0])
        n_vars = int(var_names.size)
        if keep_mask is None:
            keep_mask = np.ones(n_obs, dtype=bool)
        if keep_mask.shape[0] != n_obs:
            raise ValueError("keep_mask length does not match integrated cell count.")
        kept_n = int(np.sum(keep_mask))
        expr = np.zeros((kept_n, len(GENES)), dtype=np.float32)
        x_group = handle["X"]
        write_start = 0
        for start in range(0, n_obs, CHUNK_SIZE):
            end = min(start + CHUNK_SIZE, n_obs)
            chunk = build_csr_chunk(x_group, start, end, n_vars)
            chunk = normalize_log1p_chunk(chunk, TARGET_SUM)
            chunk_mask = keep_mask[start:end]
            if np.any(chunk_mask):
                selected = chunk[chunk_mask, :][:, gene_idx].toarray().astype(np.float32, copy=False)
                write_end = write_start + selected.shape[0]
                expr[write_start:write_end, :] = selected
                write_start = write_end
            if start == 0 or end == n_obs or (start // CHUNK_SIZE) % 10 == 0:
                logging.info("Processed expression chunk %s-%s / %s", start, end, n_obs)
    return pd.DataFrame(expr, columns=GENES)


def build_group_count_tables(obs_frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build per-GSE and per-sampleid paired/no-TRA-TRB group count summaries."""
    binary = obs_frame.loc[obs_frame["group_label"].isin(["Paired TRA/TRB", "No TRA/TRB"])].copy()
    gse = (
        binary.groupby(["source_gse_id", "group_label"], observed=True)
        .size()
        .unstack(fill_value=0)
        .reset_index()
        .sort_values("source_gse_id")
    )
    sampleid = (
        binary.assign(sampleid=binary["sampleid"].astype(str).replace("", "[blank]"))
        .groupby(["sampleid", "group_label"], observed=True)
        .size()
        .unstack(fill_value=0)
        .reset_index()
        .sort_values("sampleid")
    )
    return gse, sampleid


def compute_statistics(obs_frame: pd.DataFrame, expr: pd.DataFrame) -> pd.DataFrame:
    """Compute per-gene summaries, Wilcoxon p-values, and effect sizes."""
    paired_mask = obs_frame["group_label"].eq("Paired TRA/TRB").to_numpy()
    no_mask = obs_frame["group_label"].eq("No TRA/TRB").to_numpy()
    results: list[dict[str, float | str | int]] = []
    p_values: list[float] = []
    for gene in GENES:
        paired_values = expr.loc[paired_mask, gene].to_numpy(dtype=np.float32, copy=False)
        no_values = expr.loc[no_mask, gene].to_numpy(dtype=np.float32, copy=False)
        u_stat, p_value = stats.mannwhitneyu(
            paired_values,
            no_values,
            alternative="two-sided",
            method="asymptotic",
        )
        # Positive value means paired group tends to be higher.
        rank_biserial = 2.0 * (u_stat / (paired_values.size * no_values.size)) - 1.0
        results.append(
            {
                "gene": gene,
                "paired_n": int(paired_values.size),
                "no_trab_n": int(no_values.size),
                "paired_mean": float(np.mean(paired_values)),
                "no_trab_mean": float(np.mean(no_values)),
                "paired_median": float(np.median(paired_values)),
                "no_trab_median": float(np.median(no_values)),
                "paired_nonzero_fraction": float(np.mean(paired_values > 0)),
                "no_trab_nonzero_fraction": float(np.mean(no_values > 0)),
                "u_statistic": float(u_stat),
                "p_value": float(p_value),
                "rank_biserial": float(rank_biserial),
            }
        )
        p_values.append(float(p_value))
    stats_frame = pd.DataFrame(results)
    stats_frame["fdr_bh"] = bh_fdr(np.asarray(p_values, dtype=float))
    return stats_frame


def build_plot_sample(obs_frame: pd.DataFrame, expr: pd.DataFrame) -> pd.DataFrame:
    """Sample balanced cells for violin plotting and export the sampled values."""
    rng = np.random.default_rng(1)
    frames: list[pd.DataFrame] = []
    for group_name in ["Paired TRA/TRB", "No TRA/TRB"]:
        idx = np.flatnonzero(obs_frame["group_label"].to_numpy() == group_name)
        sample_n = min(PLOT_SAMPLE_PER_GROUP, idx.size)
        sample_idx = np.sort(rng.choice(idx, size=sample_n, replace=False))
        plot_frame = expr.iloc[sample_idx].copy()
        plot_frame["group_label"] = group_name
        frames.append(plot_frame)
    sampled = pd.concat(frames, ignore_index=True).melt(
        id_vars="group_label",
        value_vars=GENES,
        var_name="gene",
        value_name="expression",
    )
    return sampled


def write_violin_plot(plot_sample: pd.DataFrame, stats_frame: pd.DataFrame, figure_png: Path, title: str) -> None:
    """Render one four-panel violin plot with summary annotations."""
    sns.set_theme(style="whitegrid", context="talk")
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    palette = {"Paired TRA/TRB": "#D1495B", "No TRA/TRB": "#0077B6"}
    stats_lookup = stats_frame.set_index("gene")
    for ax, gene in zip(axes.flat, GENES, strict=False):
        gene_frame = plot_sample.loc[plot_sample["gene"] == gene].copy()
        sns.violinplot(
            data=gene_frame,
            x="group_label",
            y="expression",
            order=["No TRA/TRB", "Paired TRA/TRB"],
            palette=palette,
            cut=0,
            inner=None,
            linewidth=1,
            ax=ax,
        )
        sns.boxplot(
            data=gene_frame,
            x="group_label",
            y="expression",
            order=["No TRA/TRB", "Paired TRA/TRB"],
            width=0.22,
            showcaps=False,
            boxprops={"facecolor": "white", "zorder": 3},
            whiskerprops={"linewidth": 0},
            medianprops={"color": "black", "linewidth": 2},
            showfliers=False,
            ax=ax,
        )
        row = stats_lookup.loc[gene]
        ax.set_title(
            f"{gene}\nFDR={row['fdr_bh']:.2e} | paired_median={row['paired_median']:.3f} | no_median={row['no_trab_median']:.3f}",
            fontsize=12,
        )
        ax.set_xlabel("")
        ax.set_ylabel("log1p normalized expression", fontsize=11)
        ax.tick_params(axis="x", labelrotation=12)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle(title, fontsize=16)
    fig.savefig(figure_png, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def write_summary_md(
    stats_frame: pd.DataFrame,
    obs_frame: pd.DataFrame,
    summary_md: Path,
    *,
    title: str,
    trd_min: float | None,
) -> None:
    """Write a concise Markdown summary of the analysis."""
    paired_n = int(obs_frame["group_label"].eq("Paired TRA/TRB").sum())
    no_n = int(obs_frame["group_label"].eq("No TRA/TRB").sum())
    lines = [
        f"# {title}",
        "",
        f"- Integrated source: `{INTEGRATED_H5AD}`",
        f"- Raw `phase4_trd_score` filter: `{'> ' + str(trd_min) if trd_min is not None else 'none'}`",
        f"- Metadata groups: `Paired TRA/TRB` = `{paired_n:,}`, `No TRA/TRB` = `{no_n:,}`",
        "- Single-chain cells were excluded from the primary two-group comparison.",
        "",
        "## Per-gene summary",
        "",
    ]
    for _, row in stats_frame.iterrows():
        lines.append(
            f"- `{row['gene']}`: paired_median={row['paired_median']:.4f}, "
            f"no_trab_median={row['no_trab_median']:.4f}, "
            f"rank_biserial={row['rank_biserial']:.4f}, FDR={row['fdr_bh']:.2e}"
        )
    summary_md.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """Run the paired-vs-no-TRA/TRB marker comparison package."""
    args = parse_args()
    outputs = build_output_paths(args.output_prefix)
    configure_logging(outputs["run_log"])
    ensure_output_dirs()
    logging.info("Loading metadata-derived group labels")
    lookup = load_group_lookup()
    obs_frame = load_obs_group_frame(lookup)
    if args.trd_min is not None:
        keep_mask = obs_frame["phase4_trd_score"].to_numpy(dtype=np.float32, copy=False) > float(args.trd_min)
        logging.info("Applying raw phase4_trd_score > %.4f filter retained %s / %s cells", args.trd_min, int(np.sum(keep_mask)), len(obs_frame))
        obs_frame = obs_frame.loc[keep_mask].reset_index(drop=True)
    else:
        keep_mask = None
    logging.info("Loaded group labels for %s cells after any TRD filter", len(obs_frame))
    expr = extract_gene_expression(keep_mask=keep_mask)
    logging.info("Computing statistics")
    stats_frame = compute_statistics(obs_frame, expr)
    gse_counts, sampleid_counts = build_group_count_tables(obs_frame)
    plot_sample = build_plot_sample(obs_frame, expr)
    title = "Paired TRA/TRB vs No TRA/TRB Marker Expression"
    if args.trd_min is not None:
        title = f"{title} (raw TRD > {args.trd_min:g})"

    logging.info("Writing outputs")
    stats_frame.to_csv(outputs["stats_csv"], index=False)
    gse_counts.to_csv(outputs["gse_counts_csv"], index=False)
    sampleid_counts.to_csv(outputs["sampleid_counts_csv"], index=False)
    plot_sample.to_csv(outputs["plot_sample_csv"], index=False, compression="gzip")
    write_violin_plot(plot_sample, stats_frame, outputs["figure_png"], title)
    write_summary_md(stats_frame, obs_frame, outputs["summary_md"], title=title, trd_min=args.trd_min)
    logging.info(
        "Done: wrote %s, %s, %s, %s, %s, %s",
        outputs["stats_csv"],
        outputs["gse_counts_csv"],
        outputs["sampleid_counts_csv"],
        outputs["plot_sample_csv"],
        outputs["figure_png"],
        outputs["summary_md"],
    )


if __name__ == "__main__":
    main()
