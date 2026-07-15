#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Plot plus6 TNK atlas UMAP using scVI-derived broad cell-type annotations."""

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
import json
from pathlib import Path
from textwrap import fill

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


DEFAULT_H5AD = _TNK_PROJECT_ROOT / "high_speed_temp/Integrated_dataset/integrated_plus6.h5ad"
DEFAULT_ANNOTATION_KEY = "simple_annotation_plus6"
DEFAULT_OUT = _TNK_PROJECT_ROOT / "Integrated_dataset/figures/plus6/plus6_umap_by_corrected_simple_annotation_100x120mm.png"
DEFAULT_SUMMARY = _TNK_PROJECT_ROOT / "Integrated_dataset/logs/plus6/plus6_umap_by_corrected_simple_annotation_100x120mm.md"
WIDTH_MM = 100
HEIGHT_MM = 120
POINT_SIZE = 0.05
POINT_ALPHA = 0.72
DPI = 300
BROAD_ORDER = ["Other", "CD4", "CD8", "Treg", "NK", "gdT", "myeloid"]
LEGEND_ORDER = ["CD4", "CD8", "Treg", "NK", "gdT", "myeloid"]


PALETTE = {
    "CD4": "#1f77b4",
    "CD8": "#2ca02c",
    "Treg": "#9467bd",
    "NK": "#8c564b",
    "gdT": "#d62728",
    "myeloid": "#17becf",
    "Other": "#bdbdbd",
}


LABEL_DISPLAY = {
    "CD4": "CD4",
    "CD8": "CD8",
    "Treg": "Treg",
    "NK": "NK",
    "gdT": "gdT",
    "myeloid": "myeloid",
    "Other": "Other",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot plus6 TNK atlas UMAP by scVI-derived broad cell-type annotations.")
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--annotation-key", default=DEFAULT_ANNOTATION_KEY)
    parser.add_argument("--out-png", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--summary-md", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--max-cells", type=int, default=None, help="Optional reproducible downsample for debugging.")
    parser.add_argument("--seed", type=int, default=20260517)
    return parser.parse_args()


def read_string_array(obj: h5py.Dataset) -> np.ndarray:
    arr = obj[:]
    if arr.dtype.kind == "S":
        return arr.astype(str)
    if arr.dtype.kind == "O":
        return np.asarray([x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in arr], dtype=object)
    return arr.astype(str)


def read_leiden(handle: h5py.File) -> np.ndarray:
    obj = handle["obs"]["leiden"]
    if isinstance(obj, h5py.Group):
        categories = read_string_array(obj["categories"])
        codes = np.asarray(obj["codes"][:], dtype=np.int64)
        out = np.full(codes.shape[0], "NA", dtype=object)
        ok = (codes >= 0) & (codes < categories.shape[0])
        out[ok] = categories[codes[ok]]
        return out
    return read_string_array(obj)


def read_obs_column(handle: h5py.File, key: str) -> np.ndarray:
    obj = handle["obs"][key]
    if isinstance(obj, h5py.Group):
        categories = read_string_array(obj["categories"])
        codes = np.asarray(obj["codes"][:], dtype=np.int64)
        out = np.full(codes.shape[0], "", dtype=object)
        ok = (codes >= 0) & (codes < categories.shape[0])
        out[ok] = categories[codes[ok]]
        return out
    return read_string_array(obj)


def choose_indices(n_obs: int, max_cells: int | None, seed: int) -> np.ndarray | slice:
    if max_cells is None or max_cells >= n_obs:
        return slice(None)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(np.arange(n_obs), size=max_cells, replace=False))


def broad_cell_type(label: str) -> str:
    if label == "CD4_T":
        return "CD4"
    if label in {"CD8_T", "CD8_cytotoxic_T", "MAIT_like_T"}:
        return "CD8"
    if label == "Treg":
        return "Treg"
    if label == "NK_cell":
        return "NK"
    if label in {"gdT_cell", "gamma_delta_T", "gamma_delta_enriched_T"}:
        return "gdT"
    if label in {"myeloid", "myeloid_contaminant"}:
        return "myeloid"
    return "Other"


def make_label_codes(annotation: np.ndarray) -> tuple[np.ndarray, list[str], pd.DataFrame]:
    labels = pd.Series(annotation.astype(str)).map(broad_cell_type).to_numpy(dtype=object)
    counts = pd.Series(labels, dtype="string").value_counts()
    order = [label for label in BROAD_ORDER if label in set(counts.index.astype(str))]
    label_to_code = {label: i for i, label in enumerate(order)}
    codes = np.asarray([label_to_code[str(label)] for label in labels], dtype=np.int16)
    count_df = (
        counts.rename_axis("cell_type_annotation")
        .reset_index(name="n_cells")
        .assign(fraction=lambda x: x["n_cells"] / len(labels))
    )
    return codes, order, count_df


def clean_unique_count(values: np.ndarray) -> int:
    text = pd.Series(values, dtype="string").fillna("").str.strip()
    text = text[~text.str.lower().isin({"", "na", "nan", "none", "unknown", "unassigned"})]
    return int(text.nunique())


def format_int(value: int) -> str:
    return f"{int(value):,}"


def plot_umap(umap: np.ndarray, codes: np.ndarray, categories: list[str], out_png: Path, atlas_stats: dict[str, int]) -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 10,
            "axes.labelsize": 10,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "legend.title_fontsize": 10,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), dpi=DPI)
    ax = fig.add_axes([0.08, 0.30, 0.56, 0.56])
    finite = np.isfinite(umap[:, 0]) & np.isfinite(umap[:, 1])
    x = umap[:, 0]
    y = umap[:, 1]
    for code, category in enumerate(categories):
        if category == "Other":
            continue
        mask = finite & (codes == code)
        if not np.any(mask):
            continue
        ax.scatter(
            x[mask],
            y[mask],
            s=POINT_SIZE,
            c=PALETTE.get(category, "#333333"),
            alpha=POINT_ALPHA,
            linewidths=0,
            rasterized=True,
        )
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_title("TNK atlas")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    handles = []
    for category in [category for category in LEGEND_ORDER if category in categories]:
        label = fill(LABEL_DISPLAY.get(category, category.replace("_", " ")), width=14)
        handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                label=label,
                markerfacecolor=PALETTE.get(category, "#333333"),
                markeredgewidth=0,
                markersize=4.5,
            )
        )
    legend = ax.legend(
        handles=handles,
        title="Cell type",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
        frameon=False,
        handletextpad=0.35,
        handlelength=0.8,
        labelspacing=0.36,
        borderpad=0.0,
    )
    if legend is not None:
        legend._legend_box.align = "left"

    table_ax = fig.add_axes([0.06, 0.055, 0.88, 0.15])
    table_ax.axis("off")
    table = table_ax.table(
        cellText=[
            [
                format_int(atlas_stats["n_datasets"]),
                format_int(atlas_stats["n_diseases"]),
                format_int(atlas_stats["n_tissues"]),
                format_int(atlas_stats["n_cells"]),
            ]
        ],
        colLabels=["Datasets", "Diseases", "Tissues", "Cells"],
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.0, 1.22)
    for (row, _col), cell in table.get_celld().items():
        cell.set_edgecolor("#5f6b75")
        cell.set_linewidth(0.6)
        cell.set_facecolor("#f2f4f5" if row == 0 else "white")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=DPI)
    plt.close(fig)


def write_summary(
    summary_md: Path,
    *,
    h5ad: Path,
    annotation_key: str,
    out_png: Path,
    n_cells_plotted: int,
    n_cells_displayed: int,
    total_cells: int,
    count_df: pd.DataFrame,
    atlas_stats: dict[str, int],
) -> None:
    summary_md.parent.mkdir(parents=True, exist_ok=True)
    rows = [
        "| cell_type_annotation | n_cells | fraction |",
        "|---|---:|---:|",
    ]
    for row in count_df.itertuples(index=False):
        rows.append(f"| {row.cell_type_annotation} | {int(row.n_cells):,} | {float(row.fraction):.4f} |")
    payload = {
        "h5ad": str(h5ad),
        "annotation_key": annotation_key,
        "out_png": str(out_png),
        "figure_size_mm": [WIDTH_MM, HEIGHT_MM],
        "dpi": DPI,
        "text_size_pt": 10,
        "total_cells": int(total_cells),
        "cells_plotted": int(n_cells_plotted),
        "cells_displayed": int(n_cells_displayed),
            "n_datasets": int(atlas_stats["n_datasets"]),
            "n_diseases": int(atlas_stats["n_diseases"]),
            "n_tissues": int(atlas_stats["n_tissues"]),
        }
    text = "\n".join(
        [
            "# plus6 UMAP by corrected simple annotation",
            "",
            f"- Input H5AD: `{h5ad}`",
            f"- Annotation source: corrected simple annotation, `obs[\"{annotation_key}\"]`, from the plus6 H5AD",
            "- Annotation used: broad cell-type collapse of the corrected simple annotation column",
            "- Broad labels shown in the legend: `CD4`, `CD8`, `Treg`, `NK`, `gdT`, `myeloid` when present",
            "- Cells outside those broad labels are not shown on the UMAP.",
            "- TCR metadata-specific annotation proposals are not used.",
            f"- Output PNG: `{out_png}`",
            f"- Figure size: `{WIDTH_MM} x {HEIGHT_MM} mm`",
            "- Text size: `10 pt`",
            f"- Cells read: `{n_cells_plotted:,}` / `{total_cells:,}`",
            f"- Cells displayed on UMAP: `{n_cells_displayed:,}`",
            f"- Dataset count: `{atlas_stats['n_datasets']:,}`",
            f"- Disease/condition count: `{atlas_stats['n_diseases']:,}`",
            f"- Tissue type count: `{atlas_stats['n_tissues']:,}`",
            "",
            "## Cell-Type Counts",
            "",
            *rows,
            "",
            "## Metadata",
            "",
            "```json",
            json.dumps(payload, indent=2),
            "```",
        ]
    )
    summary_md.write_text(text + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    with h5py.File(args.h5ad, "r") as handle:
        n_obs = int(handle["obs"]["_index"].shape[0])
        idx = choose_indices(n_obs, args.max_cells, args.seed)
        umap = np.asarray(handle["obsm"]["X_umap"][idx], dtype=np.float32)
        annotation = read_obs_column(handle, args.annotation_key)
        source = read_obs_column(handle, "source_gse_id")
        disease_key = next(
            (key for key in ["disease_status_corrected", "disease_status", "disease", "condition"] if key in handle["obs"]),
            None,
        )
        disease = read_obs_column(handle, disease_key) if disease_key else np.asarray([""] * n_obs, dtype=object)
        tissue_key = "tissue_corrected" if "tissue_corrected" in handle["obs"] else "tissue"
        tissue = read_obs_column(handle, tissue_key)
        if not isinstance(idx, slice):
            annotation = annotation[idx]
    codes, categories, count_df = make_label_codes(annotation)
    n_cells_displayed = int(count_df.loc[count_df["cell_type_annotation"] != "Other", "n_cells"].sum())
    atlas_stats = {
        "n_datasets": clean_unique_count(source),
        "n_diseases": clean_unique_count(disease),
        "n_tissues": clean_unique_count(tissue),
        "n_cells": int(n_obs),
    }
    plot_umap(umap, codes, categories, args.out_png, atlas_stats)
    write_summary(
        args.summary_md,
        h5ad=args.h5ad,
        annotation_key=args.annotation_key,
        out_png=args.out_png,
        n_cells_plotted=int(umap.shape[0]),
        n_cells_displayed=n_cells_displayed,
        total_cells=n_obs,
        count_df=count_df,
        atlas_stats=atlas_stats,
    )
    print(args.out_png)
    print(args.summary_md)


if __name__ == "__main__":
    main()
