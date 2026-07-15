#!/usr/bin/env python3
"""Remove metadata-confirmed 10x 3' cohorts from milestone H5AD files.

This helper uses the sample-level metadata table as the source of truth.
For the current project state, every sample recorded for GSE145926 is marked
as `technology_simple = 10x 3'`, so removing that GSE from the milestone H5ADs
is equivalent to removing the confirmed 10x 3' cells.

The safer default is to write a validated new H5AD first. Replacement of the
original file can then happen only after validation succeeds.
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
PYTHON_BIN = "/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python"
SAMPLE_INFO_CSV = PROJECT_ROOT / "configs" / "datasets" / "sample_information_final_full.csv"
DEFAULT_TARGET_GSE = "GSE145926"
DEFAULT_MILESTONES = [
    PROJECT_ROOT / "Integrated_dataset" / "TNK_candidates.h5ad",
    PROJECT_ROOT / "Integrated_dataset" / "TNK_cleaned.h5ad",
    PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "integrated.h5ad",
]
SUMMARY_CSV = PROJECT_ROOT / "Integrated_dataset" / "tables" / "remove_10x3_samples_counts.csv"
SUMMARY_MD = PROJECT_ROOT / "Integrated_dataset" / "logs" / "remove_10x3_samples.md"


def configure_logging() -> None:
    """Set up simple readable logging."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        force=True,
    )


def load_and_validate_target_sample_metadata(target_gse: str) -> pd.DataFrame:
    """Load the sample metadata and confirm the target GSE is fully 10x 3'."""
    sample_info = pd.read_csv(SAMPLE_INFO_CSV, dtype=str).fillna("")
    subset = sample_info.loc[sample_info["gse"].eq(target_gse)].copy()
    if subset.empty:
        raise ValueError(f"No sample metadata rows found for target GSE {target_gse}.")

    tech_values = sorted(v for v in subset["technology_simple"].astype(str).str.strip().unique().tolist() if v)
    if tech_values != ["10x 3'"]:
        raise ValueError(
            f"Target GSE {target_gse} is not uniformly 10x 3' in sample metadata. "
            f"Observed technology_simple values: {tech_values}"
        )

    logging.info(
        "Validated %s sample metadata rows for %s; all technology_simple values are 10x 3'.",
        len(subset),
        target_gse,
    )
    return subset


def derive_output_path(path: Path, output_path: str | None) -> Path:
    """Resolve the destination path for the filtered H5AD."""
    if output_path:
        return Path(output_path)
    return path.with_name(path.stem + ".without_10x3.h5ad")


def remove_target_gse_from_h5ad(path: Path, target_gse: str, output_path: str | None) -> dict[str, object]:
    """Write a filtered H5AD without cells from the target GSE."""
    if not path.exists():
        raise FileNotFoundError(f"Milestone H5AD not found: {path}")

    logging.info("Opening %s in backed mode", path)
    adata = ad.read_h5ad(path, backed="r")
    if "source_gse_id" not in adata.obs.columns:
        raise KeyError(f"{path} is missing required obs column 'source_gse_id'.")

    target_mask = adata.obs["source_gse_id"].astype(str).eq(target_gse).to_numpy()
    before_n_obs, before_n_vars = adata.n_obs, adata.n_vars
    removed_n_obs = int(target_mask.sum())
    if removed_n_obs == 0:
        logging.info("No %s cells present in %s; leaving file unchanged.", target_gse, path)
        adata.file.close()
        return {
            "path": str(path),
            "output_path": str(derive_output_path(path, output_path)),
            "before_n_obs": before_n_obs,
            "after_n_obs": before_n_obs,
            "n_vars": before_n_vars,
            "removed_n_obs": 0,
            "target_gse_present_after": False,
            "written": False,
        }

    out_path = derive_output_path(path, output_path)
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    logging.info("Removing %s cells from %s -> %s", removed_n_obs, path, out_path)
    filtered_view = adata[~target_mask]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        raise FileExistsError(f"Refusing to overwrite existing destination {out_path}")
    if tmp_path.exists():
        tmp_path.unlink()

    # Copying directly from a backed view writes to disk without first building
    # a full in-memory AnnData copy of these 80-119 GB milestones.
    try:
        filtered = filtered_view.copy(filename=tmp_path)
        if getattr(filtered, "isbacked", False):
            filtered.file.close()
    finally:
        adata.file.close()

    verify = ad.read_h5ad(tmp_path, backed="r")
    target_present_after = bool(verify.obs["source_gse_id"].astype(str).eq(target_gse).any())
    after_n_obs = verify.n_obs
    after_n_vars = verify.n_vars
    verify.file.close()
    os.replace(tmp_path, out_path)

    if target_present_after:
        raise ValueError(f"{target_gse} still present after write of {out_path}")

    return {
        "path": str(path),
        "output_path": str(out_path),
        "before_n_obs": before_n_obs,
        "after_n_obs": after_n_obs,
        "n_vars": after_n_vars,
        "removed_n_obs": removed_n_obs,
        "target_gse_present_after": target_present_after,
        "written": True,
    }


def write_summary(target_gse: str, sample_subset: pd.DataFrame, results: list[dict[str, object]]) -> None:
    """Persist concise validation outputs under Integrated_dataset/."""
    summary_df = pd.DataFrame(results)
    SUMMARY_CSV.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_MD.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(SUMMARY_CSV, index=False)

    libraries = ", ".join(sorted(sample_subset["library_id"].drop_duplicates().tolist()))
    lines = [
        "# 10x 3' Removal Summary",
        "",
        f"- Target GSE removed: `{target_gse}`",
        f"- Sample metadata source: `{SAMPLE_INFO_CSV.relative_to(PROJECT_ROOT)}`",
        f"- Confirmed sample metadata rows: `{len(sample_subset)}`",
        f"- Confirmed libraries: `{libraries}`",
        "",
        "## Per-file results",
        "",
    ]
    for row in results:
        lines.extend(
            [
                f"- `{row['path']}`",
                f"  - output path: `{row['output_path']}`",
                f"  - before cells: `{row['before_n_obs']}`",
                f"  - removed cells: `{row['removed_n_obs']}`",
                f"  - after cells: `{row['after_n_obs']}`",
                f"  - genes: `{row['n_vars']}`",
                f"  - written: `{row['written']}`",
            ]
        )
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-gse", default=DEFAULT_TARGET_GSE, help="Target GSE to remove.")
    parser.add_argument(
        "--h5ad",
        nargs="+",
        default=[str(path) for path in DEFAULT_MILESTONES],
        help="Milestone H5AD paths to rewrite.",
    )
    parser.add_argument(
        "--output-path",
        default=None,
        help="Optional explicit output path. Use only with a single --h5ad input.",
    )
    return parser.parse_args()


def main() -> None:
    """Validate the target metadata and rewrite the milestone H5ADs."""
    configure_logging()
    args = parse_args()
    if args.output_path and len(args.h5ad) != 1:
        raise ValueError("--output-path may only be used with a single --h5ad input.")
    sample_subset = load_and_validate_target_sample_metadata(args.target_gse)
    try:
        results = [
            remove_target_gse_from_h5ad(
                Path(path),
                args.target_gse,
                args.output_path if len(args.h5ad) == 1 else None,
            )
            for path in args.h5ad
        ]
    except Exception as exc:
        logging.error("10x 3' removal failed: %s", exc)
        logging.error(traceback.format_exc())
        raise
    write_summary(args.target_gse, sample_subset, results)
    logging.info("Wrote %s and %s", SUMMARY_CSV, SUMMARY_MD)


if __name__ == "__main__":
    main()
