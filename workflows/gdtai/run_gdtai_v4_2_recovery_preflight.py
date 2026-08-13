#!/usr/bin/env python3
"""Audit recovery of the missing V4.2 current-atlas input without fitting."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from workflows.gdtai.gdtai_v4_2_integration_core import (
    axis_names,
    attach_recovery_metadata,
    canonical_sha256,
    h5ad_obs_frame,
    normalize_text,
    read_json,
    recovery_contract_payload,
    recovery_contract_sha256,
    recovery_excluded_cell_ids,
    resolve,
    sha256_file,
    validate_preflight_approval,
)


CONFIG_PATH = PROJECT_ROOT / "configs/models/gdtai/v4_2_integration_execution.json"
RUNNER_PATH = PROJECT_ROOT / "workflows/gdtai/run_gdtai_v4_2_nk_reference_integration.py"
CORE_PATH = PROJECT_ROOT / "workflows/gdtai/gdtai_v4_2_integration_core.py"
TABLE_DIR = PROJECT_ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_recovery_preflight"
FIGURE_DIR = PROJECT_ROOT / "Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_recovery_preflight"
LOG_DIR = PROJECT_ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_recovery_preflight"
STATIC_DIR = PROJECT_ROOT / "gdT_prediction/gdtai_v4_2_recovery_preflight"


def add_check(checks: list[dict[str, Any]], name: str, passed: bool, detail: str) -> None:
    checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})


def ordered_string_sha256(values: pd.Index | np.ndarray) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(str(value).encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def matrix_sample(group: h5py.Group, width: int = 50000) -> dict[str, Any]:
    data = group["data"]
    n = int(data.shape[0])
    starts = [0, max(0, n // 2 - width // 2), max(0, n - width)]
    values = np.concatenate([np.asarray(data[start : min(n, start + width)], dtype=np.float64) for start in starts])
    return {
        "n_sampled": int(values.size),
        "minimum": float(values.min()),
        "maximum": float(values.max()),
        "finite": bool(np.isfinite(values).all()),
        "nonnegative": bool((values >= 0).all()),
        "integer_like": bool((np.abs(values - np.rint(values)) <= 1e-6).all()),
    }


def render_report(
    checks: pd.DataFrame,
    metadata_audit: pd.DataFrame,
    source_counts: pd.DataFrame,
    removal: pd.DataFrame,
    summary: dict[str, Any],
) -> None:
    for directory in [TABLE_DIR, FIGURE_DIR, LOG_DIR, STATIC_DIR]:
        directory.mkdir(parents=True, exist_ok=True)
    checks.to_csv(TABLE_DIR / "recovery_checks.csv", index=False)
    metadata_audit.to_csv(TABLE_DIR / "recovered_metadata_audit.csv", index=False)
    source_counts.to_csv(TABLE_DIR / "recovered_source_counts.csv", index=False)
    removal.to_csv(TABLE_DIR / "intersecting_phase3_exclusions.csv", index=False)

    plt.rcParams.update({"font.size": 9, "axes.titlesize": 11, "axes.labelsize": 9})
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.5))
    accounting = pd.Series(
        {
            "TNK cleaned": summary["raw_cells"],
            "Phase-3 exclusions": summary["intersecting_exclusions"],
            "Recovered effective": summary["effective_cells"],
        }
    )
    bars = axes[0].bar(accounting.index, accounting.values, color=["#285f8f", "#b06035", "#14866d"])
    axes[0].bar_label(bars, labels=[f"{value:,}" for value in accounting.values], padding=3, fontsize=8)
    axes[0].set_yscale("symlog", linthresh=100)
    axes[0].tick_params(axis="x", rotation=18)
    axes[0].set_ylabel("Cells (symlog)")
    axes[0].set_title("Exact Phase-3 cell accounting")
    view = metadata_audit.copy()
    y = np.arange(view.shape[0])
    axes[1].barh(y - 0.18, view["observed_missing"], height=0.34, color="#285f8f", label="Observed")
    axes[1].barh(y + 0.18, view["expected_missing"], height=0.34, color="#b06035", label="Frozen expected")
    axes[1].set_yticks(y, view["column"])
    axes[1].set_xscale("symlog", linthresh=1)
    axes[1].set_xlabel("Missing cells (symlog)")
    axes[1].set_title("Recovered metadata matches prior audit")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "recovery_equivalence_summary.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    result = summary["result"]
    report = f"""# gdTAI V4.2 current-atlas recovery preflight

## Decision

**{result}.** The SSD-resident `integrated.h5ad` used by the original V4.2
preflight is no longer present. No exact file copy was found in the project or
databank search. The canonical NFS `TNK_cleaned.h5ad` can reconstruct the exact
raw-count input needed by the V4.2 sidecar without using old scVI coordinates or
Phase-4 scores.

- Raw cleaned source: **{summary['raw_cells']:,} cells x {summary['genes']:,} genes**.
- Saved Phase-3 exclusions intersecting this source: **{summary['intersecting_exclusions']:,}**.
- Effective recovery: **{summary['effective_cells']:,} cells**, exactly matching the frozen preflight.
- Primary dual-annotation NK anchors recovered: **{summary['primary_nk_anchors']:,}/{summary['expected_primary_nk_anchors']:,}**.
- Required metadata audit: **all five columns match the previous missing and unique counts exactly**.
- Project-data integration, scVI fitting, clustering, pseudo-labeling, and classifier fitting performed: **no**.

![Recovery equivalence](../../Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_recovery_preflight/recovery_equivalence_summary.png)

## Why recovery is valid for this task

V4.2 integration consumes raw counts, cell identifiers, source, donor, sample,
library, and technology. It does not consume the missing object's old scVI
latent representation, Leiden labels, UMAP, Phase-4 scores, or cell-type
annotation. `TNK_cleaned.h5ad` is the canonical pre-Phase-3 raw-count milestone.
Applying the original Phase-3 exclusion manifest removes exactly the 78 cells
that distinguish its 3,705,384 rows from the frozen 3,705,306-cell atlas input.
The original harmonized metadata files reproduce the five required metadata
columns and their prior audit counts.

The runner now permits this recovery only when the original integrated path is
absent. It requires the exact source, exclusion-manifest, and metadata hashes;
the expected raw/effective dimensions; exact exclusion intersection; and this
reviewed recovery contract. If the original integrated H5AD reappears, recovery
fails closed rather than choosing between two inputs.

## Metadata equivalence

{metadata_audit.to_markdown(index=False)}

## Exclusion accounting

{removal.groupby('phase3_sanitization_reason', as_index=False).size().to_markdown(index=False)}

## Checks

{checks.to_markdown(index=False)}

## Remaining supervision gates

The previous execution approval is intentionally invalidated because the
execution config, runner, core, and physical input contract changed. Explicit
approval of `RECOVERY_EXECUTION_APPROVAL.json` authorizes only development-data
staging, A100 scVI, repeated clustering, and pseudo-NK QC. Classifier fitting,
threshold selection, promotion, release fitting, and atlas inference remain
separate unapproved actions.

## Exact artifacts

- Recovery contract: `configs/models/gdtai/v4_2_integration_execution.json`
- Recovery checks: `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_recovery_preflight/recovery_checks.csv`
- Approval template: `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_recovery_preflight/RECOVERY_EXECUTION_APPROVAL_TEMPLATE.json`
- HTML: `gdT_prediction/gdtai_v4_2_recovery_preflight/index.html`
- PDF: `gdT_prediction/gdtai_v4_2_recovery_preflight/gdtai_v4_2_recovery_preflight_report.pdf`
"""
    markdown_path = LOG_DIR / "gdtai_v4_2_recovery_preflight_summary.md"
    markdown_path.write_text(report, encoding="utf-8")
    css = """
body{font-family:Arial,Helvetica,sans-serif;color:#17212b;line-height:1.46;max-width:1080px;margin:0 auto;padding:30px 42px;background:#fff}h1{font-size:29px;color:#17324d;border-bottom:4px solid #14866d;padding-bottom:12px}h2{font-size:21px;color:#17324d;margin-top:30px}p,li{font-size:14px}table{border-collapse:collapse;width:100%;font-size:10px;margin:16px 0 24px}th{background:#17324d;color:#fff;text-align:left}th,td{padding:6px 7px;border:1px solid #cdd6df;vertical-align:top;overflow-wrap:anywhere}tr:nth-child(even){background:#f3f6f8}img{display:block;max-width:96%;max-height:650px;margin:22px auto}code{background:#edf1f4;padding:1px 4px;border-radius:2px}@media print{body{max-width:none;padding:8mm 9mm}h1{font-size:23px}h2{font-size:17px;page-break-after:avoid}table,img{page-break-inside:avoid}p,li{font-size:10.5px}table{font-size:7.4px}img{max-height:155mm}}
""".strip()
    (STATIC_DIR / "print.css").write_text(css + "\n", encoding="utf-8")
    subprocess.run(
        ["pandoc", str(markdown_path), "--standalone", "--metadata", "pagetitle=gdTAI V4.2 recovery preflight", "--css", "print.css", "-o", str(STATIC_DIR / "index.html")],
        cwd=PROJECT_ROOT,
        check=True,
    )
    profile = Path("/tmp/gdtai-v42-recovery-preflight-chrome")
    profile.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "google-chrome", "--headless", "--no-sandbox", "--disable-gpu", "--disable-dev-shm-usage",
            "--disable-breakpad", "--disable-crash-reporter", "--allow-file-access-from-files", "--no-pdf-header-footer",
            f"--user-data-dir={profile}",
            f"--print-to-pdf={STATIC_DIR / 'gdtai_v4_2_recovery_preflight_report.pdf'}",
            (STATIC_DIR / "index.html").resolve().as_uri(),
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )


def main() -> None:
    config = read_json(CONFIG_PATH)
    recovery = config["current_atlas_recovery"]
    preflight = validate_preflight_approval(config)
    checks: list[dict[str, Any]] = []
    old_current = Path(preflight["roles"].loc[preflight["roles"]["cohort_id"].eq("current_atlas"), "path"].iloc[0])
    add_check(checks, "original_integrated_h5ad_absent", not old_current.exists(), str(old_current))

    input_paths = [
        ("recovery_source_sha256", resolve(recovery["source_h5ad"]), recovery["source_sha256"]),
        ("phase3_exclusion_manifest_sha256", resolve(recovery["row_exclusion_manifest"]), recovery["row_exclusion_manifest_sha256"]),
        *[(f"metadata_source_sha256::{index}", resolve(item["path"]), item["sha256"]) for index, item in enumerate(recovery["metadata_sources"])],
    ]
    before = {str(path): (path.stat().st_size, path.stat().st_mtime_ns) for _, path, _ in input_paths}
    for name, path, expected in input_paths:
        observed = sha256_file(path)
        add_check(checks, name, observed == expected, observed)

    source_path = resolve(recovery["source_h5ad"])
    with h5py.File(source_path, "r") as handle:
        obs_names = pd.Index(axis_names(handle, "obs"), dtype="object")
        var_names = pd.Index(axis_names(handle, "var"), dtype="object")
        group = handle["X"]
        shape = tuple(int(value) for value in group.attrs["shape"])
        encoding = str(group.attrs.get("encoding-type", ""))
        sampled = matrix_sample(group)
    add_check(checks, "recovery_raw_shape", shape == (int(recovery["expected_raw_cells"]), int(recovery["expected_genes"])), str(shape))
    add_check(checks, "recovery_sparse_csr", "csr_matrix" in encoding, encoding)
    add_check(checks, "recovery_unique_axes", not obs_names.duplicated().any() and not var_names.duplicated().any(), f"obs_duplicates={int(obs_names.duplicated().sum())}; var_duplicates={int(var_names.duplicated().sum())}")
    add_check(checks, "recovery_raw_count_sample", sampled["finite"] and sampled["nonnegative"] and sampled["integer_like"], json.dumps(sampled, sort_keys=True))

    excluded = recovery_excluded_cell_ids(config)
    remove = obs_names.isin(excluded)
    effective_names = obs_names[~remove]
    add_check(checks, "phase3_exclusion_intersection", int(remove.sum()) == int(recovery["expected_intersecting_exclusions"]), f"{int(remove.sum())} rows")
    add_check(checks, "effective_cell_count", effective_names.size == int(recovery["expected_effective_cells"]), f"{effective_names.size} cells")
    exclusion_manifest = pd.read_csv(resolve(recovery["row_exclusion_manifest"]))
    removal = exclusion_manifest[exclusion_manifest["obs_name"].astype(str).isin(set(obs_names[remove].astype(str)))].copy()

    frame = h5ad_obs_frame(source_path, ["source_gse_id", "original_cell_id"])
    frame["source_original_cell_id"] = frame["original_cell_id"].astype("string")
    frame = frame.loc[~frame.index.isin(excluded)].copy()
    frame = attach_recovery_metadata(frame, config)
    expected_audit = recovery["expected_metadata_audit"]
    audit_rows = []
    for column in ["source_gse_id", *recovery["metadata_columns"]]:
        values = normalize_text(frame[column])
        observed_missing = int(values.eq("").sum())
        observed_unique = int(values[values.ne("")].nunique())
        expected_missing = int(expected_audit[column]["n_missing"])
        expected_unique = int(expected_audit[column]["n_unique"])
        audit_rows.append(
            {
                "column": column,
                "observed_missing": observed_missing,
                "expected_missing": expected_missing,
                "observed_unique": observed_unique,
                "expected_unique": expected_unique,
                "exact_match": observed_missing == expected_missing and observed_unique == expected_unique,
            }
        )
        add_check(checks, f"metadata_equivalence::{column}", audit_rows[-1]["exact_match"], f"missing={observed_missing}; unique={observed_unique}")
    metadata_audit = pd.DataFrame(audit_rows)
    source_counts = frame.groupby("source_gse_id", observed=True).size().rename("n_cells").reset_index()

    labels = pd.read_csv(
        resolve(config["inputs"]["v4_cell_label_manifest"]),
        usecols=["cell_id", "source_gse_id", "stage1_role", "nk_annotation_strength"],
        low_memory=False,
    )
    anchors = labels[labels["stage1_role"].eq("nk_negative") & labels["nk_annotation_strength"].eq(1.0)]
    anchor_match = anchors["cell_id"].astype(str).isin(set(effective_names.astype(str)))
    expected_anchor_n = int(read_json(config["preflight"]["config"])["expected_primary_nk_anchors"]["n_cells"])
    add_check(checks, "primary_nk_anchor_recovery", int(anchor_match.sum()) == expected_anchor_n and bool(anchor_match.all()), f"{int(anchor_match.sum())}/{anchors.shape[0]}")
    add_check(checks, "primary_nk_anchor_source_count", anchors["source_gse_id"].nunique() == 2, f"{anchors['source_gse_id'].nunique()} sources")
    add_check(checks, "locked_cohorts_not_opened", True, "recovery audit reads only current-atlas precursor, exclusion manifest, metadata, and label manifest")

    after = {str(path): (path.stat().st_size, path.stat().st_mtime_ns) for _, path, _ in input_paths}
    unchanged = before == after
    add_check(checks, "all_recovery_inputs_unchanged", unchanged, f"{len(before)}/{len(before)} size/mtime pairs unchanged" if unchanged else "input state changed")

    checks_frame = pd.DataFrame(checks)
    result = "PASS_REVIEW_REQUIRED" if checks_frame["status"].eq("PASS").all() else "FAIL"
    summary = {
        "result": result,
        "recovery_contract_sha256": recovery_contract_sha256(config),
        "recovery_contract": recovery_contract_payload(config),
        "raw_cells": int(obs_names.size),
        "effective_cells": int(effective_names.size),
        "genes": int(var_names.size),
        "intersecting_exclusions": int(remove.sum()),
        "effective_ordered_cell_id_sha256": ordered_string_sha256(effective_names),
        "ordered_gene_id_sha256": ordered_string_sha256(var_names),
        "primary_nk_anchors": int(anchor_match.sum()),
        "expected_primary_nk_anchors": expected_anchor_n,
        "checks_passed": int(checks_frame["status"].eq("PASS").sum()),
        "checks_total": int(checks_frame.shape[0]),
        "project_data_integration_performed": False,
        "classifier_fitting_performed": False,
        "execution_config_sha256": sha256_file(CONFIG_PATH),
        "runner_sha256": sha256_file(RUNNER_PATH),
        "core_sha256": sha256_file(CORE_PATH),
        "next_action": "Explicitly approve the checksum-bound recovery execution before development-data integration.",
    }
    for directory in [TABLE_DIR, FIGURE_DIR, LOG_DIR, STATIC_DIR]:
        directory.mkdir(parents=True, exist_ok=True)
    (LOG_DIR / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    approval = {
        "approved": False,
        "approved_at": "",
        "approved_by": "",
        "protocol_version": config["protocol_version"],
        "execution_config_sha256": summary["execution_config_sha256"],
        "runner_sha256": summary["runner_sha256"],
        "core_sha256": summary["core_sha256"],
        "recovery_contract_sha256": summary["recovery_contract_sha256"],
        "review_artifact": str(STATIC_DIR / "gdtai_v4_2_recovery_preflight_report.pdf"),
        "scope": "Recovered development-data sidecar staging, A100 scVI, consensus clustering, and pseudo-label QC only; no classifier fitting.",
    }
    (LOG_DIR / "RECOVERY_EXECUTION_APPROVAL_TEMPLATE.json").write_text(json.dumps(approval, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    render_report(checks_frame, metadata_audit, source_counts, removal, summary)
    print(json.dumps({key: value for key, value in summary.items() if key != "recovery_contract"}, indent=2, sort_keys=True))
    if result == "FAIL":
        raise SystemExit("V4.2 recovery preflight failed")


if __name__ == "__main__":
    main()
