#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Validate supplementary 10x 5' per-GSE H5AD metadata/TCR schema.

This is a small project-specific validator for the supplementary intake lane.
It checks that each H5AD follows the agreed metadata and TCR field conventions.
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
from pathlib import Path

import anndata as ad
import pandas as pd


DEFAULT_INPUT_DIR = _TNK_PROJECT_ROOT / "downloads/per_gse_h5ad_with_metadata"
DEFAULT_REPORT_DIR = _TNK_PROJECT_ROOT / "Integrated_dataset/tables/supplementary_10x5_validation"
SUPPLEMENTARY_GSES = [
    "GSE179994",
    "GSE234069",
    "GSE235863",
    "GSE240865",
    "GSE287301",
    "GSE287541",
]

REQUIRED_METADATA_COLS = [
    "library_id",
    "sample_id",
    "donor_id",
    "technology_simple",
    "barcode",
]

REQUIRED_TCR_COLS = [
    "TCRseq",
    "TRA_cdr3",
    "TRA_v",
    "TRA_d",
    "TRA_j",
    "TRA_cdr3_nt",
    "TRA_clone_id",
    "TRA_umis",
    "TRA_reads",
    "TRB_cdr3",
    "TRB_v",
    "TRB_d",
    "TRB_j",
    "TRB_cdr3_nt",
    "TRB_clone_id",
    "TRB_umis",
    "TRB_reads",
]

EXPECTED_10X5 = "10x 5'"


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "paths",
        nargs="*",
        help="Specific H5AD files to validate. Defaults to all *_with_tcr.h5ad in the supplementary output directory.",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help="Directory containing supplementary per-GSE H5AD files.",
    )
    parser.add_argument(
        "--report-dir",
        type=Path,
        default=DEFAULT_REPORT_DIR,
        help="Directory for CSV/Markdown validation reports.",
    )
    return parser.parse_args()


def clean_text(series: pd.Series) -> pd.Series:
    """Normalize text-like series to plain strings with blank missing values."""
    series = series.astype("string").fillna("").str.strip()
    return series.replace({"<NA>": "", "nan": "", "None": ""})


def infer_gse_id(path: Path) -> str:
    """Infer the GSE ID from the filename."""
    return path.name.split("_")[0]


def validate_file(path: Path) -> dict[str, object]:
    """Validate one supplementary H5AD file and return summary metrics."""
    adata = ad.read_h5ad(path, backed="r")
    obs = adata.obs
    gse_id = infer_gse_id(path)

    missing_metadata = [col for col in REQUIRED_METADATA_COLS if col not in obs.columns]
    missing_tcr = [col for col in REQUIRED_TCR_COLS if col not in obs.columns]

    technology_ok = False
    technology_unique = ""
    if "technology_simple" in obs.columns:
        tech = clean_text(pd.Series(obs["technology_simple"], index=obs.index))
        technology_unique = ";".join(sorted(v for v in tech.unique().tolist() if v))
        technology_ok = bool(len(tech)) and tech.eq(EXPECTED_10X5).all()

    tcrseq_ok = False
    tcrseq_unique = ""
    if "TCRseq" in obs.columns:
        tcrseq = clean_text(pd.Series(obs["TCRseq"], index=obs.index)).str.lower()
        allowed = {"", "yes", "no"}
        tcrseq_unique = ";".join(sorted(v for v in tcrseq.unique().tolist() if v))
        tcrseq_ok = tcrseq.isin(allowed).all()

    obs_names_unique = bool(obs.index.is_unique)

    gse240865_cell_id_ok = True
    if gse_id == "GSE240865":
        obs_names = pd.Index(obs.index.astype(str))
        gse240865_cell_id_ok = obs_names.str.contains(":").all()

    result = {
        "gse_id": gse_id,
        "path": str(path),
        "cells": int(adata.n_obs),
        "genes": int(adata.n_vars),
        "obs_names_unique": obs_names_unique,
        "missing_metadata_cols": ";".join(missing_metadata),
        "missing_tcr_cols": ";".join(missing_tcr),
        "technology_simple_unique": technology_unique,
        "technology_simple_ok": technology_ok,
        "TCRseq_unique": tcrseq_unique,
        "TCRseq_ok": tcrseq_ok,
        "gse240865_library_barcode_ids_ok": gse240865_cell_id_ok,
    }
    result["passed"] = all(
        [
            obs_names_unique,
            not missing_metadata,
            not missing_tcr,
            technology_ok,
            tcrseq_ok,
            gse240865_cell_id_ok,
        ]
    )

    if getattr(adata, "file", None) is not None:
        adata.file.close()
    return result


def resolve_paths(args: argparse.Namespace) -> list[Path]:
    """Resolve input files from explicit paths or the default directory."""
    if args.paths:
        return [Path(path).resolve() for path in args.paths]
    return [
        (args.input_dir / f"{gse_id}_with_tcr.h5ad").resolve()
        for gse_id in SUPPLEMENTARY_GSES
    ]


def write_reports(report_dir: Path, summary: pd.DataFrame) -> None:
    """Write CSV and Markdown validation reports."""
    report_dir.mkdir(parents=True, exist_ok=True)
    csv_path = report_dir / "supplementary_10x5_schema_validation.csv"
    md_path = report_dir / "supplementary_10x5_schema_validation.md"

    summary.to_csv(csv_path, index=False)

    display = summary.fillna("").astype(str)
    lines_table = [
        "| " + " | ".join(display.columns) + " |",
        "| " + " | ".join(["---"] * len(display.columns)) + " |",
    ]
    for _, row in display.iterrows():
        lines_table.append("| " + " | ".join(row.tolist()) + " |")

    lines = [
        "# Supplementary 10x 5' Schema Validation",
        "",
        f"- files_checked: {len(summary)}",
        f"- files_passed: {int(summary['passed'].sum())}",
        f"- files_failed: {int((~summary['passed']).sum())}",
        "",
        "## Per-file summary",
        "",
        *lines_table,
        "",
    ]
    md_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Run validation over the requested supplementary H5AD files."""
    args = parse_args()
    paths = resolve_paths(args)
    if not paths:
        raise FileNotFoundError("No supplementary *_with_tcr.h5ad files found to validate.")

    summary = pd.DataFrame([validate_file(path) for path in paths]).sort_values("gse_id")
    write_reports(args.report_dir, summary)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
