#!/usr/bin/env python3
"""Extract no-TCR-gene GSEs and remove them from milestone H5AD files.

This helper uses the pre-merge TRAB/TRD audit result as the source of truth for
which GSEs lack the canonical module genes required for Phase 4 scoring. It:

1. extracts those cells from the current integrated milestone into a standalone
   `No_TCR_Gene_GSEs.h5ad`
2. rewrites the candidate / cleaned / integrated milestones without those GSEs
3. validates the rewritten files before replacing the originals
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
import os
from pathlib import Path
import traceback

import anndata as ad
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
TARGET_GSES = [
    "GSE139555",
    "GSE144469",
    "GSE155222",
    "GSE155223",
    "GSE161918",
    "GSE168163",
    "GSE168859",
    "GSE179994",
    "GSE221776",
    "GSE227709",
    "GSE229858",
    "GSE232240",
    "GSE240361",
    "GSE241783",
]
DEFAULT_MILESTONES = [
    PROJECT_ROOT / "Integrated_dataset" / "TNK_candidates.h5ad",
    PROJECT_ROOT / "Integrated_dataset" / "TNK_cleaned.h5ad",
    Path("/ssd/tnk_phase3/Integrated_dataset/TNK_cleaned.h5ad"),
    PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "integrated.h5ad",
]
DEFAULT_EXTRACT_SOURCE = PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "integrated.h5ad"
DEFAULT_EXTRACT_OUTPUT = PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "No_TCR_Gene_GSEs.h5ad"
SUMMARY_CSV = PROJECT_ROOT / "Integrated_dataset" / "tables" / "no_tcr_gene_gse_removal_counts.csv"
SUMMARY_MD = PROJECT_ROOT / "Integrated_dataset" / "logs" / "no_tcr_gene_gse_removal.md"


def configure_logging() -> None:
    """Set up readable logging."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        force=True,
    )


def derive_tmp_path(output_path: Path) -> Path:
    """Build a temporary path on the same filesystem as the destination."""
    return output_path.with_suffix(output_path.suffix + ".tmp")


def extract_target_subset(source_path: Path, output_path: Path, target_gses: list[str]) -> dict[str, object]:
    """Write the standalone H5AD containing only the target GSEs."""
    if not source_path.exists():
        raise FileNotFoundError(f"Extract source H5AD not found: {source_path}")
    if output_path.exists():
        raise FileExistsError(f"Refusing to overwrite existing extract destination {output_path}")

    tmp_path = derive_tmp_path(output_path)
    if tmp_path.exists():
        tmp_path.unlink()

    logging.info("Opening extract source %s in backed mode", source_path)
    adata = ad.read_h5ad(source_path, backed="r")
    if "source_gse_id" not in adata.obs.columns:
        raise KeyError(f"{source_path} is missing required obs column 'source_gse_id'.")

    target_mask = adata.obs["source_gse_id"].astype(str).isin(target_gses).to_numpy()
    removed_n_obs = int(target_mask.sum())
    if removed_n_obs == 0:
        adata.file.close()
        raise ValueError("No target GSE cells found in the extract source.")

    before_n_obs, before_n_vars = adata.n_obs, adata.n_vars
    logging.info("Writing standalone extracted subset with %s cells to %s", removed_n_obs, output_path)
    subset_view = adata[target_mask]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        subset_view.write_h5ad(tmp_path, convert_strings_to_categoricals=False)
    finally:
        adata.file.close()

    verify = ad.read_h5ad(tmp_path, backed="r")
    after_n_obs, after_n_vars = verify.shape
    present_gses = sorted(set(verify.obs["source_gse_id"].astype(str).unique().tolist()))
    verify.file.close()
    os.replace(tmp_path, output_path)

    unexpected = sorted(set(present_gses) - set(target_gses))
    if unexpected:
        raise ValueError(f"Unexpected non-target GSEs present in extracted subset: {unexpected[:10]}")

    return {
        "path": str(source_path),
        "output_path": str(output_path),
        "before_n_obs": before_n_obs,
        "after_n_obs": after_n_obs,
        "n_vars": after_n_vars,
        "removed_n_obs": removed_n_obs,
        "written": True,
        "operation": "extract_targets",
    }


def rewrite_without_targets(path: Path, target_gses: list[str]) -> dict[str, object]:
    """Rewrite one milestone H5AD without the target GSEs."""
    if not path.exists():
        raise FileNotFoundError(f"Milestone H5AD not found: {path}")

    tmp_path = derive_tmp_path(path)
    if tmp_path.exists():
        tmp_path.unlink()

    logging.info("Opening %s in backed mode", path)
    adata = ad.read_h5ad(path, backed="r")
    if "source_gse_id" not in adata.obs.columns:
        raise KeyError(f"{path} is missing required obs column 'source_gse_id'.")

    target_mask = adata.obs["source_gse_id"].astype(str).isin(target_gses).to_numpy()
    before_n_obs, before_n_vars = adata.n_obs, adata.n_vars
    removed_n_obs = int(target_mask.sum())
    if removed_n_obs == 0:
        adata.file.close()
        logging.info("No target GSE cells present in %s; leaving file unchanged.", path)
        return {
            "path": str(path),
            "output_path": str(path),
            "before_n_obs": before_n_obs,
            "after_n_obs": before_n_obs,
            "n_vars": before_n_vars,
            "removed_n_obs": 0,
            "written": False,
            "operation": "remove_targets",
        }

    logging.info("Removing %s target cells from %s", removed_n_obs, path)
    filtered_view = adata[~target_mask]
    try:
        filtered_view.write_h5ad(tmp_path, convert_strings_to_categoricals=False)
    finally:
        adata.file.close()

    verify = ad.read_h5ad(tmp_path, backed="r")
    target_present_after = bool(verify.obs["source_gse_id"].astype(str).isin(target_gses).any())
    after_n_obs, after_n_vars = verify.shape
    verify.file.close()
    if target_present_after:
        raise ValueError(f"Target GSEs still present after rewrite of {path}")

    os.replace(tmp_path, path)
    return {
        "path": str(path),
        "output_path": str(path),
        "before_n_obs": before_n_obs,
        "after_n_obs": after_n_obs,
        "n_vars": after_n_vars,
        "removed_n_obs": removed_n_obs,
        "written": True,
        "operation": "remove_targets",
    }


def write_summary(target_gses: list[str], results: list[dict[str, object]]) -> None:
    """Write concise CSV and markdown summaries."""
    summary_df = pd.DataFrame(results)
    SUMMARY_CSV.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_MD.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(SUMMARY_CSV, index=False)

    lines = [
        "# No-TCR-Gene GSE Removal Summary",
        "",
        "Target GSEs removed from milestones:",
        "",
    ]
    lines.extend([f"- `{gse}`" for gse in target_gses])
    lines.extend(
        [
            "",
            "## Per-file results",
            "",
        ]
    )
    for row in results:
        lines.extend(
            [
                f"- `{row['operation']}` on `{row['path']}`",
                f"  - output path: `{row['output_path']}`",
                f"  - before cells: `{row['before_n_obs']}`",
                f"  - affected cells: `{row['removed_n_obs']}`",
                f"  - after cells: `{row['after_n_obs']}`",
                f"  - genes: `{row['n_vars']}`",
                f"  - written: `{row['written']}`",
            ]
        )
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--extract-source",
        type=Path,
        default=DEFAULT_EXTRACT_SOURCE,
        help="H5AD used to build the standalone No_TCR_Gene_GSEs object.",
    )
    parser.add_argument(
        "--extract-output",
        type=Path,
        default=DEFAULT_EXTRACT_OUTPUT,
        help="Output path for the standalone No_TCR_Gene_GSEs H5AD.",
    )
    parser.add_argument(
        "--h5ad",
        nargs="+",
        default=[str(path) for path in DEFAULT_MILESTONES],
        help="Milestone H5AD paths to rewrite in place.",
    )
    parser.add_argument(
        "--target-gse",
        nargs="+",
        default=TARGET_GSES,
        help="GSEs to extract and remove.",
    )
    return parser.parse_args()


def main() -> None:
    """Extract the target GSEs and rewrite the milestone H5ADs."""
    configure_logging()
    args = parse_args()
    target_gses = sorted(set(args.target_gse))
    try:
        results = [extract_target_subset(args.extract_source, args.extract_output, target_gses)]
        for path in args.h5ad:
            results.append(rewrite_without_targets(Path(path), target_gses))
    except Exception as exc:
        logging.error("No-TCR-gene GSE removal failed: %s", exc)
        logging.error(traceback.format_exc())
        raise

    write_summary(target_gses, results)
    logging.info("Wrote %s and %s", SUMMARY_CSV, SUMMARY_MD)


if __name__ == "__main__":
    main()
