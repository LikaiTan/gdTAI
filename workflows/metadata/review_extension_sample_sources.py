#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Review disease-aware specimen context for standalone extension H5ADs."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from workflows.metadata.sample_source_refinement_workflow import (  # noqa: E402
    DEFAULT_RULES,
    file_sha256,
    review_h5ad,
    write_json,
)


DEFAULT_MANIFEST = PROJECT_ROOT / "data/interim/extension_intake/built_h5ads_manifest.csv"
DEFAULT_OUTPUT_ROOT = (
    PROJECT_ROOT / "Integrated_dataset/tables/sample_source_refinement/extensions"
)
DEFAULT_LOG = (
    PROJECT_ROOT
    / "Integrated_dataset/logs/sample_source_refinement/extension_sample_source_review.md"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--rules", type=Path, default=DEFAULT_RULES)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--log", type=Path, default=DEFAULT_LOG)
    parser.add_argument("--chunk-size", type=int, default=100_000)
    return parser.parse_args()


def resolve_h5ad_path(value: str, manifest_path: Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path.resolve()
    project_path = (PROJECT_ROOT / path).resolve()
    if project_path.exists():
        return project_path
    return (manifest_path.parent / path).resolve()


def load_inputs(manifest_path: Path) -> pd.DataFrame:
    frame = pd.read_csv(manifest_path, dtype=str, keep_default_na=False)
    cohort_column = next(
        (name for name in ("cohort_id", "dataset_id") if name in frame.columns),
        None,
    )
    path_column = next(
        (name for name in ("h5ad_path", "output_h5ad", "output_path") if name in frame.columns),
        None,
    )
    if cohort_column is None or path_column is None:
        raise ValueError("Manifest requires cohort_id/dataset_id and h5ad_path/output_h5ad")
    inputs = frame[[cohort_column, path_column]].rename(
        columns={cohort_column: "cohort_id", path_column: "h5ad_path"}
    )
    if inputs.empty or inputs["cohort_id"].eq("").any() or inputs["h5ad_path"].eq("").any():
        raise ValueError("Manifest contains blank or no extension inputs")
    if inputs["cohort_id"].duplicated().any():
        raise ValueError("Manifest contains duplicate cohort IDs")
    return inputs


def render_log(summary: pd.DataFrame, gate_pass: bool, output_root: Path) -> str:
    status = "PASS" if gate_pass else "FAIL"
    columns = [
        "cohort_id",
        "cells",
        "unresolved_cells",
        "tumor_context_violations",
        "tumor_context_gate_pass",
    ]
    table = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in summary[columns].itertuples(index=False, name=None):
        table.append("| " + " | ".join(str(value) for value in row) + " |")
    lines = [
        "# Extension Tumor-Project Specimen-Context Review",
        "",
        f"- Gate: **{status}**",
        "- Mode: read-only review; no H5AD was modified",
        f"- Review root: `{output_root}`",
        "",
        "## Cohorts",
        "",
        *table,
        "",
        "The gate requires disease-aware specimen labels for every configured tumor-project row. "
        "Original `tissue` and `diagnosis` columns remain unchanged.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    manifest_path = args.manifest.expanduser().resolve()
    rules_path = args.rules.expanduser().resolve()
    output_root = args.output_root.expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    inputs = load_inputs(manifest_path)

    rows: list[dict[str, object]] = []
    aggregate: list[pd.DataFrame] = []
    for item in inputs.itertuples(index=False):
        cohort_id = str(item.cohort_id)
        h5ad_path = resolve_h5ad_path(str(item.h5ad_path), manifest_path)
        cohort_output = output_root / cohort_id
        review_manifest_path = review_h5ad(
            h5ad_path,
            rules_path,
            cohort_output,
            chunk_size=args.chunk_size,
        )
        review_manifest = json.loads(review_manifest_path.read_text(encoding="utf-8"))
        audit_path = cohort_output / "sample_source_refined_tumor_project_audit.csv"
        audit = pd.read_csv(audit_path, dtype={"source_gse_id": str})
        if not audit.empty:
            audit.insert(0, "cohort_id", cohort_id)
            aggregate.append(audit)
        rows.append(
            {
                "cohort_id": cohort_id,
                "cells": int(review_manifest["obs_rows"]),
                "unresolved_cells": int(review_manifest["unresolved_cells"]),
                "tumor_context_violations": int(
                    review_manifest["tumor_context_violations"]
                ),
                "tumor_context_gate_pass": bool(
                    review_manifest["tumor_context_gate_pass"]
                ),
                "review_manifest": str(review_manifest_path),
                "review_manifest_sha256": file_sha256(review_manifest_path),
            }
        )

    summary = pd.DataFrame(rows)
    gate_pass = bool(summary["tumor_context_gate_pass"].all())
    summary_path = output_root / "extension_sample_source_review_summary.csv"
    summary.to_csv(summary_path, index=False)
    aggregate_path = output_root / "extension_tumor_project_context_counts.csv"
    if aggregate:
        pd.concat(aggregate, ignore_index=True).to_csv(aggregate_path, index=False)
    else:
        pd.DataFrame(
            columns=[
                "cohort_id",
                "source_gse_id",
                "sample_source_refined",
                "context_compliant",
                "cell_n",
            ]
        ).to_csv(aggregate_path, index=False)

    run_manifest = {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "mode": "review_only",
        "input_manifest": str(manifest_path),
        "input_manifest_sha256": file_sha256(manifest_path),
        "rules": str(rules_path),
        "rules_sha256": file_sha256(rules_path),
        "tumor_context_gate_pass": gate_pass,
        "tumor_context_violations": int(summary["tumor_context_violations"].sum()),
        "cohorts": rows,
        "artifacts": {
            summary_path.name: file_sha256(summary_path),
            aggregate_path.name: file_sha256(aggregate_path),
        },
    }
    run_manifest_path = output_root / "extension_sample_source_review_manifest.json"
    write_json(run_manifest_path, run_manifest)
    args.log.parent.mkdir(parents=True, exist_ok=True)
    args.log.write_text(render_log(summary, gate_pass, output_root), encoding="utf-8")
    if not gate_pass:
        raise RuntimeError(
            "Tumor-project specimen-context gate failed; inspect the review outputs"
        )


if __name__ == "__main__":
    main()
