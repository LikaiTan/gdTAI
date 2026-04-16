#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Build plus6 γδ-focused report assets and targeted UMAP panels."""

from __future__ import annotations

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent
INPUT_H5AD = PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "integrated_plus6.h5ad"
FIGURE_DIR = PROJECT_ROOT / "Integrated_dataset" / "figures" / "plus6"
TABLE_DIR = PROJECT_ROOT / "Integrated_dataset" / "tables" / "plus6"
LOG_DIR = PROJECT_ROOT / "Integrated_dataset" / "logs" / "plus6"

UMAP_PNG = FIGURE_DIR / "plus6_umap_paired_tcr_sorted_gdt.png"
SUMMARY_CSV = TABLE_DIR / "plus6_gdt_candidate_statistics.csv"
OVERLAP_CSV = TABLE_DIR / "plus6_gdt_candidate_overlap_gt0p4.csv"
PAIRED_BY_TISSUE_CSV = TABLE_DIR / "plus6_gdt_paired_gdtcr_by_tissue.csv"
CRITERIA_BY_TISSUE_CSV = TABLE_DIR / "plus6_gdt_three_criteria_by_tissue.csv"
SUMMARY_MD = LOG_DIR / "plus6_gdt_candidate_statistics.md"

TRD_MINUS_TRAB_THRESHOLD = 0.4
SAMPLE_N = 300_000
RANDOM_SEED = 0
FIGURE_DPI = 300


def dataframe_to_markdown_fallback(df: pd.DataFrame) -> str:
    columns = [str(column) for column in df.columns]
    rows = [["" if pd.isna(value) else str(value) for value in row] for row in df.itertuples(index=False, name=None)]
    widths = [len(column) for column in columns]
    for row in rows:
        for idx, value in enumerate(row):
            widths[idx] = max(widths[idx], len(value))

    def format_row(row: list[str]) -> str:
        return "| " + " | ".join(value.ljust(widths[idx]) for idx, value in enumerate(row)) + " |"

    header = format_row(columns)
    separator = "| " + " | ".join("-" * widths[idx] for idx in range(len(widths))) + " |"
    body = [format_row(row) for row in rows]
    return "\n".join([header, separator, *body])


def read_string_dataset(dataset: h5py.Dataset) -> np.ndarray:
    values = dataset[:]
    return np.asarray([value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values], dtype=object)


def read_obs_column(handle: h5py.File, column: str) -> np.ndarray:
    obj = handle["obs"][column]
    if isinstance(obj, h5py.Group) and obj.attrs.get("encoding-type") == "categorical":
        categories = read_string_dataset(obj["categories"])
        codes = obj["codes"][:]
        out = np.full(codes.shape, "", dtype=object)
        valid = codes >= 0
        out[valid] = categories[codes[valid]]
        return out
    raw = obj[:]
    if raw.dtype.kind in {"S", "O", "U"}:
        return np.asarray([value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in raw], dtype=object)
    return np.asarray(raw)


def read_bool_obs(handle: h5py.File, column: str) -> np.ndarray:
    raw = read_obs_column(handle, column)
    if raw.dtype == bool:
        return np.asarray(raw, dtype=bool)
    lowered = np.char.lower(raw.astype(str))
    return np.asarray(lowered == "true", dtype=bool)


def choose_sample_indices(n_obs: int, sample_n: int, seed: int) -> np.ndarray:
    if n_obs <= sample_n:
        return np.arange(n_obs, dtype=np.int64)
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(np.arange(n_obs), size=sample_n, replace=False))


def clean_tissue(values: np.ndarray) -> np.ndarray:
    as_str = np.asarray(values, dtype=object).astype(str)
    as_str = np.where(np.isin(np.char.lower(as_str), ["", "nan", "none", "na"]), "unknown", as_str)
    return as_str


def build_summary_tables() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    with h5py.File(INPUT_H5AD, "r") as handle:
        n_obs = int(handle["obs"]["_index"].shape[0])
        umap_all = np.asarray(handle["obsm"]["X_umap"][:], dtype=np.float32)
        paired_ab_all = read_bool_obs(handle, "has_TRA_TRB_paired")
        paired_gd_all = read_bool_obs(handle, "has_TRG_TRD_paired")
        sorted_gdt_all = read_bool_obs(handle, "Sorted_gdT")
        any_ab_all = read_bool_obs(handle, "has_any_ab_tcr")
        simple_annotation = np.asarray(read_obs_column(handle, "simple_annotation_plus6"), dtype=object).astype(str)
        trd_minus_trab = np.asarray(handle["obs"]["phase4_trd_minus_trab"][:], dtype=np.float32)
        if "tissue_corrected" in handle["obs"]:
            tissue = clean_tissue(read_obs_column(handle, "tissue_corrected"))
            tissue_column = "tissue_corrected"
        else:
            tissue = clean_tissue(read_obs_column(handle, "tissue"))
            tissue_column = "tissue"

    doublet_all = paired_ab_all & paired_gd_all
    criterion_a = sorted_gdt_all
    criterion_b = paired_gd_all & ~doublet_all
    criterion_c = (trd_minus_trab > TRD_MINUS_TRAB_THRESHOLD) & (~any_ab_all) & (simple_annotation != "NK_cell")
    any_three = criterion_a | criterion_b | criterion_c
    gdt_annotation = simple_annotation == "gdT_cell"

    summary_df = pd.DataFrame(
        [
            {"metric": "total_cells", "value": int(n_obs)},
            {"metric": "sorted_gdT_true", "value": int(criterion_a.sum())},
            {"metric": "paired_TRG_TRD_not_doublet", "value": int(criterion_b.sum())},
            {
                "metric": f"TRD_minus_TRAB_gt_{str(TRD_MINUS_TRAB_THRESHOLD).replace('.', 'p')}_no_productive_TRA_TRB_not_NK",
                "value": int(criterion_c.sum()),
            },
            {"metric": "doublet_proxy_paired_ab_and_paired_gd", "value": int(doublet_all.sum())},
            {"metric": "union_any_of_three_criteria", "value": int(any_three.sum())},
            {"metric": "gdT_cell_annotation_total", "value": int(gdt_annotation.sum())},
            {"metric": "gdT_cell_annotation_with_paired_TRG_TRD", "value": int((gdt_annotation & paired_gd_all).sum())},
            {"metric": "gdT_cell_annotation_meeting_any_three", "value": int((gdt_annotation & any_three).sum())},
        ]
    )
    summary_df.to_csv(SUMMARY_CSV, index=False)

    overlap_df = pd.DataFrame(
        [
            {"category": "sorted_gdT_only", "n_cells": int((criterion_a & ~criterion_b & ~criterion_c).sum())},
            {"category": "paired_TRG_TRD_not_doublet_only", "n_cells": int((~criterion_a & criterion_b & ~criterion_c).sum())},
            {"category": f"TRD_minus_TRAB_gt_{str(TRD_MINUS_TRAB_THRESHOLD).replace('.', 'p')}_only", "n_cells": int((~criterion_a & ~criterion_b & criterion_c).sum())},
            {"category": "sorted_gdT_and_paired_TRG_TRD_only", "n_cells": int((criterion_a & criterion_b & ~criterion_c).sum())},
            {"category": f"sorted_gdT_and_TRD_minus_TRAB_gt_{str(TRD_MINUS_TRAB_THRESHOLD).replace('.', 'p')}_only", "n_cells": int((criterion_a & ~criterion_b & criterion_c).sum())},
            {"category": f"paired_TRG_TRD_and_TRD_minus_TRAB_gt_{str(TRD_MINUS_TRAB_THRESHOLD).replace('.', 'p')}_only", "n_cells": int((~criterion_a & criterion_b & criterion_c).sum())},
            {"category": "all_three_criteria", "n_cells": int((criterion_a & criterion_b & criterion_c).sum())},
        ]
    )
    overlap_df.to_csv(OVERLAP_CSV, index=False)

    paired_by_tissue = (
        pd.DataFrame({"tissue": tissue, "is_gdt": gdt_annotation, "paired_gd": paired_gd_all})
        .loc[lambda df: df["is_gdt"]]
        .groupby("tissue", dropna=False)
        .agg(gdt_cells=("is_gdt", "size"), gdt_with_paired_TRG_TRD=("paired_gd", "sum"))
        .reset_index()
    )
    paired_by_tissue["fraction_paired_TRG_TRD_within_gdt"] = np.where(
        paired_by_tissue["gdt_cells"] > 0,
        paired_by_tissue["gdt_with_paired_TRG_TRD"] / paired_by_tissue["gdt_cells"],
        np.nan,
    )
    paired_by_tissue = paired_by_tissue.sort_values(["gdt_with_paired_TRG_TRD", "gdt_cells"], ascending=[False, False]).reset_index(drop=True)
    paired_by_tissue.to_csv(PAIRED_BY_TISSUE_CSV, index=False)

    criteria_by_tissue = (
        pd.DataFrame({"tissue": tissue, "is_gdt": gdt_annotation, "meets_any_three": any_three})
        .loc[lambda df: df["is_gdt"]]
        .groupby("tissue", dropna=False)
        .agg(gdt_cells=("is_gdt", "size"), gdt_meeting_any_three=("meets_any_three", "sum"))
        .reset_index()
    )
    criteria_by_tissue["fraction_meeting_any_three_within_gdt"] = np.where(
        criteria_by_tissue["gdt_cells"] > 0,
        criteria_by_tissue["gdt_meeting_any_three"] / criteria_by_tissue["gdt_cells"],
        np.nan,
    )
    criteria_by_tissue = criteria_by_tissue.sort_values(["gdt_meeting_any_three", "gdt_cells"], ascending=[False, False]).reset_index(drop=True)
    criteria_by_tissue.to_csv(CRITERIA_BY_TISSUE_CSV, index=False)

    idx = choose_sample_indices(n_obs, SAMPLE_N, RANDOM_SEED)
    plot_df = pd.DataFrame(
        {
            "umap1": umap_all[idx, 0],
            "umap2": umap_all[idx, 1],
            "paired_ab": paired_ab_all[idx],
            "paired_gd": paired_gd_all[idx],
            "sorted_gdt": sorted_gdt_all[idx],
        }
    )

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
    for ax, column, title in [
        (axes[0], "paired_ab", "paired TRA/TRB"),
        (axes[1], "paired_gd", "paired TRG/TRD"),
        (axes[2], "sorted_gdt", "Sorted_gdT = True"),
    ]:
        other_df = plot_df.loc[~plot_df[column]]
        target_df = plot_df.loc[plot_df[column]]
        ax.scatter(
            other_df["umap1"],
            other_df["umap2"],
            s=2,
            c="lightgrey",
            linewidths=0,
            rasterized=True,
            alpha=0.5,
            label="other cells",
        )
        ax.scatter(
            target_df["umap1"],
            target_df["umap2"],
            s=2,
            c="#D62728",
            linewidths=0,
            rasterized=True,
            alpha=0.85,
            label=title,
        )
        ax.set_title(title)
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.legend(loc="upper right", frameon=True, fontsize=8)
    fig.savefig(UMAP_PNG, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)

    lines = [
        "# plus6 γδ-focused statistics",
        "",
        f"- tissue column used: `{tissue_column}`",
        f"- total cells: `{n_obs:,}`",
        f"- Sorted_gdT = True: `{int(criterion_a.sum()):,}`",
        f"- paired TRG/TRD but not doublet: `{int(criterion_b.sum()):,}`",
        f"- TRD - TRAB > `{TRD_MINUS_TRAB_THRESHOLD}` and no productive TRA/TRB and not NK: `{int(criterion_c.sum()):,}`",
        f"- cells meeting at least one of the three criteria: `{int(any_three.sum()):,}`",
        f"- gdT-cell annotations: `{int(gdt_annotation.sum()):,}`",
        f"- gdT-cell annotations with paired TRG/TRD: `{int((gdt_annotation & paired_gd_all).sum()):,}`",
        f"- gdT-cell annotations meeting at least one of the three criteria: `{int((gdt_annotation & any_three).sum()):,}`",
        "",
        "## Overlap breakdown",
        "",
        dataframe_to_markdown_fallback(overlap_df),
        "",
        "## gdT cells with paired gdTCR by tissue",
        "",
        dataframe_to_markdown_fallback(paired_by_tissue.head(30)),
        "",
        "## gdT cells meeting at least one of the three criteria by tissue",
        "",
        dataframe_to_markdown_fallback(criteria_by_tissue.head(30)),
        "",
    ]
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"saved_png\t{UMAP_PNG}")
    print(f"saved_csv\t{SUMMARY_CSV}")
    print(f"saved_csv\t{OVERLAP_CSV}")
    print(f"saved_csv\t{PAIRED_BY_TISSUE_CSV}")
    print(f"saved_csv\t{CRITERIA_BY_TISSUE_CSV}")
    print(f"saved_md\t{SUMMARY_MD}")


if __name__ == "__main__":
    build_summary_tables()
