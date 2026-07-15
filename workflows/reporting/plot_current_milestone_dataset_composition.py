#!/usr/bin/env python3
"""Plot per-GSE composition for the current merged milestone object.

The script prefers `Integrated_dataset/integrated.h5ad` when it exists, then
falls back to `TNK_cleaned.h5ad`, then `TNK_candidates.h5ad`.

It writes:
- one per-GSE composition table
- one cell-level source summary table
- several PNG figures
- one Markdown summary

Source comparison is resolved at sample/library/cell assignment time when
possible. No cell is ever labeled as coming from two tissues.
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
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


# Config
PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
FIGURE_DIR = OUTPUT_ROOT / "figures"
TABLE_DIR = OUTPUT_ROOT / "tables"
LOG_DIR = OUTPUT_ROOT / "logs"

METADATA_CSV = PROJECT_ROOT / "data" / "registry" / "metadata_review" / "total_meta_downsample_v3.csv"
LIBRARY_INFO_CSV = PROJECT_ROOT / "configs" / "datasets" / "sample_information_final_full.csv"
REGISTRY_CSV = PROJECT_ROOT / "configs" / "datasets" / "integration_inputs.csv"

MILESTONE_CANDIDATES = [
    OUTPUT_ROOT / "integrated.h5ad",
    OUTPUT_ROOT / "TNK_cleaned.h5ad",
    OUTPUT_ROOT / "TNK_candidates.h5ad",
]

DEFAULT_OUTPUT_PREFIX = "current_milestone"


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for optional GSE exclusion and output prefixing."""
    parser = argparse.ArgumentParser(
        description="Plot per-GSE cell/sample composition for the current merged milestone object."
    )
    parser.add_argument(
        "--exclude-gse",
        nargs="*",
        default=[],
        help="Optional list of GSE IDs to exclude from the analysis outputs.",
    )
    parser.add_argument(
        "--output-prefix",
        default=DEFAULT_OUTPUT_PREFIX,
        help="Filename prefix for written tables, figures, and summary outputs.",
    )
    return parser.parse_args()


def build_output_paths(prefix: str) -> dict[str, Path]:
    """Build one output-path bundle for the requested prefix."""
    return {
        "composition_table": TABLE_DIR / f"{prefix}_dataset_composition.csv",
        "source_summary_table": TABLE_DIR / f"{prefix}_source_group_cell_summary.csv",
        "summary_md": LOG_DIR / f"{prefix}_dataset_composition_summary.md",
        "cell_count_fig": FIGURE_DIR / f"{prefix}_gse_cell_counts.png",
        "sample_count_fig": FIGURE_DIR / f"{prefix}_gse_sample_counts.png",
        "source_group_fig": FIGURE_DIR / f"{prefix}_source_group_cell_comparison.png",
        "gse_source_split_fig": FIGURE_DIR / f"{prefix}_gse_source_split.png",
    }


def ensure_output_dirs() -> None:
    """Create canonical output directories if needed."""
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)


def configure_plotting() -> None:
    """Set publication-readable plot defaults."""
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


def select_milestone_h5ad() -> Path:
    """Return the highest-priority milestone object that currently exists."""
    for path in MILESTONE_CANDIDATES:
        if path.exists():
            return path
    raise FileNotFoundError("No current milestone H5AD found under Integrated_dataset/")


def normalize_text(value: str) -> str:
    """Normalize one metadata string for simple rule-based matching."""
    text = str(value).strip()
    return "" if text.lower() in {"", "nan", "none", "<na>"} else text


def classify_source_text(value: str) -> str:
    """Classify one tissue-like string into a broad source group."""
    text = normalize_text(value).lower()
    if not text:
        return "unknown"

    blood_like = any(token in text for token in ["pbmc", "blood", "peripheral blood", "bld"])
    tissue_like = any(
        token in text
        for token in [
            "tumor",
            "biopsy",
            "rectum",
            "rectal",
            "colon",
            "colorect",
            "nat",
            "tissue",
            "pinch",
            "mucosa",
            "bowel",
            "gut",
            "intestin",
        ]
    )

    if blood_like and not tissue_like:
        return "blood_pbmc"
    if tissue_like and not blood_like:
        return "tumor_tissue"
    return "unknown"


def collapse_source_groups(groups: list[str]) -> str:
    """Collapse several source calls into one resolved fallback group."""
    clean = {group for group in groups if group in {"blood_pbmc", "tumor_tissue"}}
    if len(clean) == 1:
        return next(iter(clean))
    return "unknown"


def load_metadata_sample_counts() -> pd.DataFrame:
    """Count unique harmonized sample IDs per GSE."""
    meta = pd.read_csv(METADATA_CSV, usecols=["dataset_id", "sample_id"])
    meta["dataset_id"] = meta["dataset_id"].astype(str).str.strip()
    meta["sample_id"] = meta["sample_id"].fillna("").astype(str).str.strip()
    summary = (
        meta.groupby("dataset_id", dropna=False)["sample_id"]
        .agg(lambda s: len({value for value in s if value}))
        .reset_index()
        .rename(columns={"dataset_id": "gse_id", "sample_id": "metadata_sample_id_n"})
    )
    return summary


def load_library_counts() -> pd.DataFrame:
    """Count unique library IDs per GSE from the library manifest."""
    libraries = pd.read_csv(LIBRARY_INFO_CSV, usecols=["gse", "library_id"])
    summary = (
        libraries.groupby("gse", dropna=False)["library_id"]
        .nunique()
        .reset_index()
        .rename(columns={"gse": "gse_id", "library_id": "library_manifest_n"})
    )
    return summary


def load_source_resolution_tables() -> tuple[dict[tuple[str, str], str], dict[tuple[str, str], str], dict[str, str], pd.DataFrame]:
    """Build sample-, library-, and dataset-level source-resolution mappings."""
    meta = pd.read_csv(METADATA_CSV, usecols=["dataset_id", "sample_id", "tissue"])
    meta["dataset_id"] = meta["dataset_id"].astype(str).str.strip()
    meta["sample_id"] = meta["sample_id"].fillna("").astype(str).str.strip()
    meta["tissue"] = meta["tissue"].fillna("").astype(str)
    meta["source_group"] = meta["tissue"].map(classify_source_text)

    sample_map: dict[tuple[str, str], str] = {}
    for row in (
        meta.groupby(["dataset_id", "sample_id"], dropna=False)["source_group"]
        .agg(lambda s: collapse_source_groups(sorted(set(s.astype(str)))))
        .reset_index()
        .itertuples(index=False)
    ):
        if normalize_text(row.sample_id):
            sample_map[(row.dataset_id, row.sample_id)] = row.source_group

    meta_dataset_fallback = (
        meta.groupby("dataset_id", dropna=False)["source_group"]
        .agg(lambda s: collapse_source_groups(sorted(set(s.astype(str)))))
        .reset_index()
        .rename(columns={"dataset_id": "gse_id", "source_group": "meta_dataset_source_group"})
    )

    tissue_values = (
        meta.groupby("dataset_id", dropna=False)["tissue"]
        .agg(lambda s: sorted({normalize_text(x) for x in s if normalize_text(x)}))
        .reset_index()
        .rename(columns={"dataset_id": "gse_id", "tissue": "tissue_values"})
    )
    tissue_values["tissue_values_joined"] = tissue_values["tissue_values"].apply(lambda values: "; ".join(values[:20]))

    library_info = pd.read_csv(LIBRARY_INFO_CSV, usecols=["gse", "library_id", "sample_type"])
    library_info["gse"] = library_info["gse"].astype(str).str.strip()
    library_info["library_id"] = library_info["library_id"].fillna("").astype(str).str.strip()
    library_info["sample_type"] = library_info["sample_type"].fillna("").astype(str)
    library_info["source_group"] = library_info["sample_type"].map(classify_source_text)

    library_map: dict[tuple[str, str], str] = {}
    for row in (
        library_info.groupby(["gse", "library_id"], dropna=False)["source_group"]
        .agg(lambda s: collapse_source_groups(sorted(set(s.astype(str)))))
        .reset_index()
        .itertuples(index=False)
    ):
        if normalize_text(row.library_id):
            library_map[(row.gse, row.library_id)] = row.source_group

    info_dataset_fallback = (
        library_info.groupby("gse", dropna=False)["source_group"]
        .agg(lambda s: collapse_source_groups(sorted(set(s.astype(str)))))
        .reset_index()
        .rename(columns={"gse": "gse_id", "source_group": "info_dataset_source_group"})
    )

    dataset_fallback_table = meta_dataset_fallback.merge(info_dataset_fallback, on="gse_id", how="outer").merge(
        tissue_values[["gse_id", "tissue_values_joined"]],
        on="gse_id",
        how="outer",
    )
    dataset_fallback_table["meta_dataset_source_group"] = dataset_fallback_table["meta_dataset_source_group"].fillna("unknown")
    dataset_fallback_table["info_dataset_source_group"] = dataset_fallback_table["info_dataset_source_group"].fillna("unknown")
    dataset_fallback_table["dataset_source_group"] = dataset_fallback_table.apply(
        lambda row: row["meta_dataset_source_group"]
        if row["meta_dataset_source_group"] in {"blood_pbmc", "tumor_tissue"}
        else row["info_dataset_source_group"],
        axis=1,
    )
    dataset_fallback_table["dataset_source_group"] = dataset_fallback_table["dataset_source_group"].fillna("unknown")
    dataset_fallback_table["tissue_values_joined"] = dataset_fallback_table["tissue_values_joined"].fillna("")

    dataset_fallback_map = dict(
        zip(dataset_fallback_table["gse_id"].astype(str), dataset_fallback_table["dataset_source_group"].astype(str))
    )
    return sample_map, library_map, dataset_fallback_map, dataset_fallback_table


def resolve_cell_source_group(
    gse_id: str,
    sample_label: str,
    library_label: str,
    sample_map: dict[tuple[str, str], str],
    library_map: dict[tuple[str, str], str],
    dataset_fallback_map: dict[str, str],
) -> str:
    """Resolve one cell into blood_pbmc, tumor_tissue, or unknown."""
    sample_label = normalize_text(sample_label)
    library_label = normalize_text(library_label)

    if sample_label:
        resolved = sample_map.get((gse_id, sample_label), "unknown")
        if resolved in {"blood_pbmc", "tumor_tissue"}:
            return resolved

    if library_label:
        resolved = library_map.get((gse_id, library_label), "unknown")
        if resolved in {"blood_pbmc", "tumor_tissue"}:
            return resolved

    return dataset_fallback_map.get(gse_id, "unknown")


def scan_current_milestone_obs(
    h5ad_path: Path,
    excluded_gse: set[str],
    sample_map: dict[tuple[str, str], str],
    library_map: dict[tuple[str, str], str],
    dataset_fallback_map: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Count per-GSE cells plus cell-level source assignments from the milestone H5AD."""
    cell_counter: Counter[str] = Counter()
    sample_labels: dict[str, set[str]] = defaultdict(set)
    library_labels: dict[str, set[str]] = defaultdict(set)
    donor_labels: dict[str, set[str]] = defaultdict(set)
    per_gse_source_counter: dict[str, Counter[str]] = defaultdict(Counter)
    global_source_counter: Counter[str] = Counter()

    with h5py.File(h5ad_path, "r") as handle:
        obs = handle["obs"]
        gse_ds = obs["source_gse_id"]
        sample_ds = obs["phase1_sample_label"] if "phase1_sample_label" in obs else None
        library_ds = obs["phase1_library_label"] if "phase1_library_label" in obs else None
        donor_ds = obs["phase1_donor_label"] if "phase1_donor_label" in obs else None

        total_rows = gse_ds.shape[0]
        chunk_size = 500_000
        for start in range(0, total_rows, chunk_size):
            stop = min(start + chunk_size, total_rows)
            gse_values = gse_ds.asstr()[start:stop]
            sample_values = sample_ds.asstr()[start:stop] if sample_ds is not None else [""] * (stop - start)
            library_values = (
                library_ds.asstr()[start:stop] if library_ds is not None else [""] * (stop - start)
            )
            donor_values = donor_ds.asstr()[start:stop] if donor_ds is not None else [""] * (stop - start)

            for gse_id, sample_label, library_label, donor_label in zip(
                gse_values,
                sample_values,
                library_values,
                donor_values,
            ):
                gse_id = normalize_text(gse_id)
                if not gse_id or gse_id in excluded_gse:
                    continue

                sample_label = normalize_text(sample_label)
                library_label = normalize_text(library_label)
                donor_label = normalize_text(donor_label)

                cell_counter[gse_id] += 1
                if sample_label:
                    sample_labels[gse_id].add(sample_label)
                if library_label:
                    library_labels[gse_id].add(library_label)
                if donor_label:
                    donor_labels[gse_id].add(donor_label)

                source_group = resolve_cell_source_group(
                    gse_id,
                    sample_label,
                    library_label,
                    sample_map,
                    library_map,
                    dataset_fallback_map,
                )
                per_gse_source_counter[gse_id][source_group] += 1
                global_source_counter[source_group] += 1

    rows = []
    for gse_id in sorted(cell_counter):
        rows.append(
            {
                "gse_id": gse_id,
                "cell_n": int(cell_counter[gse_id]),
                "milestone_sample_label_n": len(sample_labels[gse_id]),
                "milestone_library_label_n": len(library_labels[gse_id]),
                "milestone_donor_label_n": len(donor_labels[gse_id]),
                "cell_blood_pbmc_n": int(per_gse_source_counter[gse_id]["blood_pbmc"]),
                "cell_tumor_tissue_n": int(per_gse_source_counter[gse_id]["tumor_tissue"]),
                "cell_unknown_n": int(per_gse_source_counter[gse_id]["unknown"]),
            }
        )
    milestone_counts = pd.DataFrame(rows)

    source_summary = pd.DataFrame(
        [
            {"source_group": "blood_pbmc", "cell_n": int(global_source_counter["blood_pbmc"])},
            {"source_group": "tumor_tissue", "cell_n": int(global_source_counter["tumor_tissue"])},
            {"source_group": "unknown", "cell_n": int(global_source_counter["unknown"])},
        ]
    )

    gse_source_rows = []
    for gse_id in sorted(per_gse_source_counter):
        for source_group in ["blood_pbmc", "tumor_tissue", "unknown"]:
            gse_source_rows.append(
                {
                    "gse_id": gse_id,
                    "source_group": source_group,
                    "cell_n": int(per_gse_source_counter[gse_id][source_group]),
                }
            )
    gse_source_long = pd.DataFrame(gse_source_rows)
    return milestone_counts, source_summary, gse_source_long


def choose_sample_count_basis(row: pd.Series) -> tuple[int, str]:
    """Choose one practical sample-count value and record where it came from."""
    if int(row.get("metadata_sample_id_n", 0)) > 0:
        return int(row["metadata_sample_id_n"]), "metadata_sample_id"
    if int(row.get("milestone_sample_label_n", 0)) > 0:
        return int(row["milestone_sample_label_n"]), "milestone_sample_label"
    if int(row.get("library_manifest_n", 0)) > 0:
        return int(row["library_manifest_n"]), "library_manifest"
    if int(row.get("milestone_library_label_n", 0)) > 0:
        return int(row["milestone_library_label_n"]), "milestone_library_label"
    if int(row.get("milestone_donor_label_n", 0)) > 0:
        return int(row["milestone_donor_label_n"]), "milestone_donor_label"
    return 0, "missing"


def build_dataset_composition_table(
    h5ad_path: Path,
    excluded_gse: set[str],
    sample_map: dict[tuple[str, str], str],
    library_map: dict[tuple[str, str], str],
    dataset_fallback_map: dict[str, str],
    dataset_fallback_table: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Assemble per-GSE composition plus cell-level source summaries."""
    milestone_counts, source_summary, gse_source_long = scan_current_milestone_obs(
        h5ad_path,
        excluded_gse,
        sample_map,
        library_map,
        dataset_fallback_map,
    )
    metadata_sample_counts = load_metadata_sample_counts()
    library_counts = load_library_counts()
    registry = pd.read_csv(REGISTRY_CSV, usecols=["gse_id", "source_root"])

    composition = (
        registry.merge(milestone_counts, on="gse_id", how="left")
        .merge(dataset_fallback_table[["gse_id", "dataset_source_group", "tissue_values_joined"]], on="gse_id", how="left")
        .merge(metadata_sample_counts, on="gse_id", how="left")
        .merge(library_counts, on="gse_id", how="left")
    )

    for column in [
        "cell_n",
        "metadata_sample_id_n",
        "library_manifest_n",
        "milestone_sample_label_n",
        "milestone_library_label_n",
        "milestone_donor_label_n",
        "cell_blood_pbmc_n",
        "cell_tumor_tissue_n",
        "cell_unknown_n",
    ]:
        composition[column] = composition[column].fillna(0).astype(int)

    composition["gse_in_milestone_obs"] = composition["cell_n"] > 0
    composition["dataset_source_group"] = composition["dataset_source_group"].fillna("unknown")
    composition["tissue_values_joined"] = composition["tissue_values_joined"].fillna("")

    chosen = composition.apply(choose_sample_count_basis, axis=1, result_type="expand")
    composition["sample_count_for_plot"] = chosen[0].astype(int)
    composition["sample_count_basis"] = chosen[1].astype(str)

    if excluded_gse:
        composition = composition.loc[~composition["gse_id"].isin(excluded_gse)].copy()
        gse_source_long = gse_source_long.loc[~gse_source_long["gse_id"].isin(excluded_gse)].copy()

    composition = composition.sort_values("cell_n", ascending=False).reset_index(drop=True)
    source_summary = source_summary.sort_values("cell_n", ascending=False).reset_index(drop=True)
    return composition, source_summary, gse_source_long


def write_tables(composition: pd.DataFrame, source_summary: pd.DataFrame, output_paths: dict[str, Path]) -> None:
    """Write the main composition table and the source-group summary table."""
    composition.to_csv(output_paths["composition_table"], index=False)
    source_summary.to_csv(output_paths["source_summary_table"], index=False)


def write_figures(
    composition: pd.DataFrame,
    source_summary: pd.DataFrame,
    milestone_name: str,
    output_paths: dict[str, Path],
    gse_source_long: pd.DataFrame,
) -> None:
    """Write the requested GSE-level and source-comparison figures."""
    configure_plotting()

    cell_plot = composition.sort_values("cell_n", ascending=True)
    fig1, ax1 = plt.subplots(figsize=(9, 10))
    sns.barplot(data=cell_plot, x="cell_n", y="gse_id", color="#2f6c8f", ax=ax1)
    ax1.set_title(f"Cells Per GSE In {milestone_name}")
    ax1.set_xlabel("Cells")
    ax1.set_ylabel("GSE ID")
    fig1.tight_layout()
    fig1.savefig(output_paths["cell_count_fig"])
    plt.close(fig1)

    sample_plot = composition.sort_values("sample_count_for_plot", ascending=True)
    fig2, ax2 = plt.subplots(figsize=(9, 10))
    sns.barplot(
        data=sample_plot,
        x="sample_count_for_plot",
        y="gse_id",
        hue="sample_count_basis",
        dodge=False,
        palette={
            "metadata_sample_id": "#3f7cac",
            "milestone_sample_label": "#5a9367",
            "library_manifest": "#c98900",
            "milestone_library_label": "#9b59b6",
            "milestone_donor_label": "#bc6c25",
            "missing": "#bdbdbd",
        },
        ax=ax2,
    )
    ax2.set_title(f"Sample Counts Per GSE In {milestone_name}")
    ax2.set_xlabel("Sample-like groups")
    ax2.set_ylabel("GSE ID")
    ax2.legend(title="Count basis", loc="lower right")
    fig2.tight_layout()
    fig2.savefig(output_paths["sample_count_fig"])
    plt.close(fig2)

    order = ["blood_pbmc", "tumor_tissue", "unknown"]
    source_plot = source_summary.copy()
    source_plot["source_group"] = pd.Categorical(source_plot["source_group"], categories=order, ordered=True)
    source_plot = source_plot.sort_values("source_group")
    fig3, ax3 = plt.subplots(figsize=(7.5, 4.8))
    sns.barplot(
        data=source_plot,
        x="source_group",
        y="cell_n",
        hue="source_group",
        palette={
            "blood_pbmc": "#247ba0",
            "tumor_tissue": "#d1495b",
            "unknown": "#bdbdbd",
        },
        dodge=False,
        legend=False,
        ax=ax3,
    )
    ax3.set_title(f"Cells By Cell-Level Source Group In {milestone_name}")
    ax3.set_xlabel("Resolved source group")
    ax3.set_ylabel("Cells")
    for idx, value in enumerate(source_plot["cell_n"].tolist()):
        ax3.text(idx, value, f"{int(value):,}", ha="center", va="bottom", fontsize=9)
    fig3.tight_layout()
    fig3.savefig(output_paths["source_group_fig"])
    plt.close(fig3)

    if not gse_source_long.empty:
        split_plot = gse_source_long.copy()
        pivot = (
            split_plot.pivot(index="gse_id", columns="source_group", values="cell_n")
            .fillna(0)
            .reindex(columns=order, fill_value=0)
            .sort_values(order, ascending=True)
        )
        fig4, ax4 = plt.subplots(figsize=(10, 11))
        left = None
        colors = {"blood_pbmc": "#247ba0", "tumor_tissue": "#d1495b", "unknown": "#bdbdbd"}
        for source_group in order:
            values = pivot[source_group].to_numpy()
            ax4.barh(pivot.index, values, left=left, color=colors[source_group], label=source_group)
            left = values if left is None else left + values
        ax4.set_title(f"Per-GSE Cell-Level Source Split In {milestone_name}")
        ax4.set_xlabel("Cells")
        ax4.set_ylabel("GSE ID")
        ax4.legend(title="Resolved source")
        fig4.tight_layout()
        fig4.savefig(output_paths["gse_source_split_fig"])
        plt.close(fig4)


def write_summary(
    composition: pd.DataFrame,
    source_summary: pd.DataFrame,
    milestone_path: Path,
    output_paths: dict[str, Path],
    excluded_gse: list[str],
) -> None:
    """Write a concise Markdown summary of the composition analysis."""
    with output_paths["summary_md"].open("w", encoding="utf-8") as handle:
        handle.write("# Current Milestone Dataset Composition Summary\n\n")
        handle.write(f"- Milestone object used: `{milestone_path}`\n")
        handle.write(
            "- This run prefers `integrated.h5ad`, then `TNK_cleaned.h5ad`, then `TNK_candidates.h5ad`.\n"
        )
        handle.write(
            f"- Excluded GSEs for this report: `{', '.join(excluded_gse)}`\n"
            if excluded_gse
            else "- Excluded GSEs for this report: `none`\n"
        )
        handle.write(f"- GSEs listed in the current registry: {composition.shape[0]}\n")
        handle.write(
            f"- GSEs with >0 cells in the current milestone object: {int(composition['gse_in_milestone_obs'].sum())}\n"
        )
        handle.write(f"- Total cells in the current milestone object: {int(composition['cell_n'].sum()):,}\n\n")

        handle.write("## Source-group cell totals\n\n")
        for _, row in source_summary.iterrows():
            handle.write(f"- {row['source_group']}: {int(row['cell_n']):,}\n")

        handle.write("\n## Top GSEs by cell count\n\n")
        for _, row in composition.nlargest(10, "cell_n").iterrows():
            handle.write(
                f"- {row['gse_id']}: {int(row['cell_n']):,} cells; "
                f"sample_count_for_plot={int(row['sample_count_for_plot'])} "
                f"({row['sample_count_basis']}); "
                f"blood_pbmc={int(row['cell_blood_pbmc_n'])}, "
                f"tumor_tissue={int(row['cell_tumor_tissue_n'])}, "
                f"unknown={int(row['cell_unknown_n'])}\n"
            )

        handle.write("\n## Notes\n\n")
        handle.write(
            "- Cell-level source groups are resolved from harmonized sample-level `tissue` metadata when possible, "
            "then from library-level `sample_type`, then from a dataset-level fallback if the whole dataset is consistently blood-like or tissue-like.\n"
        )
        handle.write(
            "- Sample counts use the harmonized `sample_id` when available, then fall back to milestone sample labels, "
            "then library-level counts, then donor labels.\n"
        )
        handle.write(
            "- No cell is assigned to a `mixed` source group. Datasets with multiple source types are split at the sample/library level when metadata permits; otherwise unresolved cells remain `unknown`.\n"
        )
        handle.write(
            "- If `integrated.h5ad` does not yet exist, this analysis reflects the latest currently available merged milestone object rather than a true integrated latent-space object.\n"
        )


def main() -> None:
    """Run the composition analysis and write tables, figures, and summary."""
    args = parse_args()
    ensure_output_dirs()
    milestone_path = select_milestone_h5ad()
    excluded_gse = sorted({gse.strip() for gse in args.exclude_gse if gse.strip()})
    output_paths = build_output_paths(args.output_prefix)

    sample_map, library_map, dataset_fallback_map, dataset_fallback_table = load_source_resolution_tables()
    composition, source_summary, gse_source_long = build_dataset_composition_table(
        milestone_path,
        set(excluded_gse),
        sample_map,
        library_map,
        dataset_fallback_map,
        dataset_fallback_table,
    )
    write_tables(composition, source_summary, output_paths)
    write_figures(composition, source_summary, milestone_path.name, output_paths, gse_source_long)
    write_summary(composition, source_summary, milestone_path, output_paths, excluded_gse)

    print("Wrote:")
    print(output_paths["composition_table"])
    print(output_paths["source_summary_table"])
    print(output_paths["summary_md"])
    print(output_paths["cell_count_fig"])
    print(output_paths["sample_count_fig"])
    print(output_paths["source_group_fig"])
    print(output_paths["gse_source_split_fig"])


if __name__ == "__main__":
    main()
