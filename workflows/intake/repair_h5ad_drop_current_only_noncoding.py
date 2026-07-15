#!/usr/bin/env python3
"""Repair a processed H5AD by dropping current-only noncoding genes before count rescue.

This helper is intentionally narrow. It is meant for datasets where:
- rebuilt raw counts are available,
- obs alignment is valid,
- current var is not a subset of rebuilt var only because of a small set of
  current-only noncoding/pseudogene-style features.
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
import json
import re
from pathlib import Path

import pandas as pd
import scanpy as sc

from repair_h5ad_from_selected_inputs import (
    build_repaired,
    compare_axes,
    rebuild_counts,
)


NONCODING_RE = re.compile(
    r"(^AC\d+|^AL\d+|^AP\d+|^LINC|^MIR|^SNOR|^SNORA|^SNORD|^SCARNA|^RNU|P$|P\d+$|^AKR1C8P$)"
)


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(description="Drop current-only noncoding genes then repair from selected inputs.")
    parser.add_argument("gse_id", help="GSE accession to repair.")
    parser.add_argument(
        "--write",
        action="store_true",
        help="Write repaired H5AD in place. Default is dry-run only.",
    )
    parser.add_argument(
        "--drop-all-current-only",
        action="store_true",
        help="Drop all current-only genes instead of only noncoding-like genes.",
    )
    return parser.parse_args()


def is_noncoding_like(gene: str) -> bool:
    """Heuristic for current-only genes that are safe candidates for exclusion."""
    gene = str(gene)
    return bool(NONCODING_RE.search(gene))


def main() -> None:
    """Run the narrow current-only noncoding drop rescue."""
    args = parse_args()
    gse = args.gse_id
    project_dir = _TNK_PROJECT_ROOT / "analysis_26GSE_V4/scanpy_projects" / gse
    current_path = project_dir / "outputs" / f"{gse}_scanpy_processed.h5ad"

    current = sc.read_h5ad(current_path)
    rebuilt, _ = rebuild_counts(project_dir)
    comparison = compare_axes(current, rebuilt)

    current_only = pd.Index(current.var_names.astype(str)).difference(pd.Index(rebuilt.var_names.astype(str))).tolist()
    if args.drop_all_current_only:
        drop_genes = list(current_only)
    else:
        drop_genes = [gene for gene in current_only if is_noncoding_like(gene)]
    kept = current[:, [gene not in set(drop_genes) for gene in current.var_names]].copy()
    post = compare_axes(kept, rebuilt)

    result = {
        "gse_id": gse,
        "current_only_n": len(current_only),
        "drop_gene_n": len(drop_genes),
        "drop_genes_preview": drop_genes[:20],
        "post_obs_equal": post["obs_equal"],
        "post_obs_key_set_equal": post["obs_key_set_equal"],
        "post_var_current_subset_of_rebuilt": post["var_current_subset_of_rebuilt"],
        "post_current_n_vars": post["current_n_vars"],
        "rebuilt_n_vars": post["rebuilt_n_vars"],
    }
    print(json.dumps(result, indent=2))

    if not args.write:
        return

    if not (post["obs_equal"] or post["obs_key_set_equal"]):
        raise ValueError("Obs alignment still invalid after dropping current-only noncoding genes")
    if not post["var_current_subset_of_rebuilt"]:
        raise ValueError("Var still not a subset after dropping current-only noncoding genes")

    repaired = build_repaired(kept, rebuilt)
    temp_path = current_path.with_suffix(current_path.suffix + ".tmp")
    if temp_path.exists():
        temp_path.unlink()
    repaired.write_h5ad(temp_path)
    temp_path.replace(current_path)
    print(f"Repaired {current_path}")


if __name__ == "__main__":
    main()
