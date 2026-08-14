#!/usr/bin/env python3
"""Audit the stage-specific V4.2 SSD floor after staging and scVI fitting."""

from __future__ import annotations

import json
import shutil
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
    read_json,
    recovery_contract_sha256,
    resolve,
    resource_contract_payload,
    resource_contract_sha256,
    sha256_file,
    validate_recovery_preflight,
)


CONFIG_PATH = PROJECT_ROOT / "configs/models/gdtai/v4_2_integration_execution.json"
RUNNER_PATH = PROJECT_ROOT / "workflows/gdtai/run_gdtai_v4_2_nk_reference_integration.py"
CORE_PATH = PROJECT_ROOT / "workflows/gdtai/gdtai_v4_2_integration_core.py"
TABLE_DIR = PROJECT_ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_cluster_resource_preflight"
FIGURE_DIR = PROJECT_ROOT / "Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_cluster_resource_preflight"
LOG_DIR = PROJECT_ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_cluster_resource_preflight"
STATIC_DIR = PROJECT_ROOT / "gdT_prediction/gdtai_v4_2_cluster_resource_preflight"


def check(rows: list[dict[str, Any]], name: str, passed: bool, detail: str) -> None:
    rows.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})


def render(checks: pd.DataFrame, summary: dict[str, Any]) -> None:
    for directory in [TABLE_DIR, FIGURE_DIR, LOG_DIR, STATIC_DIR]:
        directory.mkdir(parents=True, exist_ok=True)
    checks.to_csv(TABLE_DIR / "resource_amendment_checks.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    labels = ["Observed free", "Original floor", "Cluster floor", "Worst-case output"]
    values = [
        summary["observed_ssd_free_gib"],
        summary["original_minimum_ssd_free_gib"],
        summary["cluster_minimum_ssd_free_gib"],
        summary["worst_case_cluster_output_gib"],
    ]
    colors = ["#197278", "#a23e48", "#315b7d", "#d98f39"]
    bars = ax.bar(labels, values, color=colors)
    ax.bar_label(bars, labels=[f"{value:.1f}" for value in values], padding=3)
    ax.set_ylabel("GiB")
    ax.set_title("V4.2 cluster-stage storage contract")
    ax.tick_params(axis="x", rotation=15)
    fig.tight_layout()
    figure = FIGURE_DIR / "cluster_storage_contract.png"
    fig.savefig(figure, dpi=240, bbox_inches="tight")
    plt.close(fig)

    report = f"""# gdTAI V4.2 cluster-stage resource amendment

## Decision

**{summary['result']}.** Sparse staging and A100 scVI fitting completed under
the original **300 GiB** SSD floor. RAPIDS clustering did not start because an
unrelated shared-SSD BAM sort reduced free space. After that process completed
and self-cleaned, free space stabilized at **{summary['observed_ssd_free_gib']:.1f} GiB**.

The amendment retains 300 GiB for `prepare` and `fit`, and freezes **150 GiB**
for `cluster` and `consensus`. This is an operational storage guard only. It
does not change cells, HVGs, expression, scVI parameters, latent coordinates,
clustering grids, random seeds, anchors, pseudo-label rules, or validation
roles.

![Storage contract](../../Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_cluster_resource_preflight/cluster_storage_contract.png)

## Completed evidence

- Staged matrix: **{summary['n_cells']:,} cells x {summary['n_hvgs']:,} HVGs**, **{summary['matrix_nnz']:,}** nonzero counts.
- Locked cohorts in staged object: **no**.
- Primary NK anchors: **{summary['primary_nk_anchors']:,}**.
- scVI latent: **{summary['n_cells']:,} x {summary['n_latent']}**, finite, A100, no CPU fallback.
- Source H5AD size/mtime checks: **{summary['unchanged_source_files']}/{summary['source_files']} unchanged**.
- Worst-case uncompressed partition payload plus reserve: **{summary['worst_case_cluster_output_gib']:.2f} GiB**.
- Free-space margin above the amended floor and reserve: **{summary['post_reserve_margin_gib']:.1f} GiB**.

## Checks

{checks.to_markdown(index=False)}

## Supervision boundary

Activating `CLUSTER_EXECUTION_APPROVAL.json` authorizes only the already
approved RAPIDS repeated clustering and pseudo-NK consensus QC using the saved
staged H5AD and latent representation. It does not authorize re-staging,
re-fitting scVI, classifier fitting, threshold selection, promotion, release
fitting, or atlas inference.

## Artifacts

- Checks: `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_cluster_resource_preflight/resource_amendment_checks.csv`
- Approval template: `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_cluster_resource_preflight/CLUSTER_EXECUTION_APPROVAL_TEMPLATE.json`
- HTML: `gdT_prediction/gdtai_v4_2_cluster_resource_preflight/index.html`
- PDF: `gdT_prediction/gdtai_v4_2_cluster_resource_preflight/gdtai_v4_2_cluster_resource_preflight_report.pdf`
"""
    markdown = LOG_DIR / "gdtai_v4_2_cluster_resource_preflight_summary.md"
    markdown.write_text(report, encoding="utf-8")
    css = """
body{font-family:Arial,Helvetica,sans-serif;color:#17212b;line-height:1.45;max-width:1080px;margin:0 auto;padding:30px 42px}h1{font-size:28px;color:#17324d;border-bottom:4px solid #197278;padding-bottom:11px}h2{font-size:20px;color:#17324d;margin-top:28px}p,li{font-size:14px}table{border-collapse:collapse;width:100%;font-size:10px}th{background:#17324d;color:white}th,td{padding:6px;border:1px solid #ccd5dd;vertical-align:top;overflow-wrap:anywhere}th:nth-child(2),td:nth-child(2){width:8%;white-space:nowrap}tr:nth-child(even){background:#f3f6f8}img{display:block;max-width:92%;max-height:620px;margin:20px auto}code{background:#edf1f4;padding:1px 4px}@media print{body{padding:8mm 9mm}h1{font-size:23px}h2{font-size:17px;page-break-after:avoid}p,li{font-size:10.5px}table{font-size:7.5px}table,img{page-break-inside:avoid}img{max-height:145mm}}
""".strip()
    (STATIC_DIR / "print.css").write_text(css + "\n", encoding="utf-8")
    subprocess.run(
        ["pandoc", str(markdown), "--standalone", "--metadata", "pagetitle=gdTAI V4.2 cluster resource amendment", "--css", "print.css", "-o", str(STATIC_DIR / "index.html")],
        cwd=PROJECT_ROOT,
        check=True,
    )
    profile = Path("/tmp/gdtai-v42-cluster-resource-chrome")
    profile.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "google-chrome", "--headless", "--no-sandbox", "--disable-gpu",
            "--disable-dev-shm-usage", "--disable-breakpad", "--disable-crash-reporter",
            "--allow-file-access-from-files", "--no-pdf-header-footer",
            f"--user-data-dir={profile}",
            f"--print-to-pdf={STATIC_DIR / 'gdtai_v4_2_cluster_resource_preflight_report.pdf'}",
            (STATIC_DIR / "index.html").resolve().as_uri(),
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )


def main() -> None:
    config = read_json(CONFIG_PATH)
    validate_recovery_preflight(config)
    paths = {
        "staged": Path(config["resources"]["ssd_root"]) / "development_hvg_counts.h5ad",
        "latent": Path(config["resources"]["ssd_root"]) / "X_scVI.npy",
        "model": Path(config["resources"]["ssd_root"]) / "scvi_model/model.pt",
        "partitions": Path(config["resources"]["ssd_root"]) / "cluster_partitions.npz",
        "prepare": resolve(config["outputs"]["log_dir"]) / "prepare_summary.json",
        "fit": resolve(config["outputs"]["log_dir"]) / "fit_summary.json",
        "source_state": resolve(config["outputs"]["table_dir"]) / "development_input_file_state.csv",
    }
    prepare = read_json(paths["prepare"])
    fit = read_json(paths["fit"])
    rows: list[dict[str, Any]] = []

    check(rows, "prepare_passed", prepare.get("status") == "PASS", str(paths["prepare"]))
    check(rows, "fit_passed", fit.get("status") == "PASS", str(paths["fit"]))
    check(rows, "staged_h5ad_present", paths["staged"].is_file(), str(paths["staged"]))
    staged_hash = sha256_file(paths["staged"])
    check(rows, "staged_h5ad_sha256", staged_hash == prepare.get("staged_h5ad_sha256"), staged_hash)
    with h5py.File(paths["staged"], "r") as handle:
        shape = tuple(int(value) for value in handle["X"].attrs["shape"])
        encoding = str(handle["X"].attrs.get("encoding-type", ""))
        nnz = int(handle["X/data"].shape[0])
    check(rows, "staged_shape", shape == (int(prepare["n_cells"]), int(prepare["n_hvgs"])), str(shape))
    check(rows, "staged_sparse_csr", "csr_matrix" in encoding, encoding)
    check(rows, "staged_nnz", nnz == int(prepare["matrix_nnz"]), f"{nnz} nnz")
    check(rows, "locked_cohorts_absent", prepare.get("locked_cohorts_included") is False, "prepare summary false")

    check(rows, "latent_present", paths["latent"].is_file(), str(paths["latent"]))
    latent_hash = sha256_file(paths["latent"])
    check(rows, "latent_sha256", latent_hash == fit.get("latent_sha256"), latent_hash)
    latent = np.load(paths["latent"], mmap_mode="r")
    latent_finite = bool(np.isfinite(latent).all())
    check(rows, "latent_shape_and_finite", latent.shape == (int(fit["n_cells"]), int(fit["n_latent"])) and latent_finite, f"shape={latent.shape}; finite={latent_finite}")
    check(rows, "gpu_without_cpu_fallback", fit.get("gpu", {}).get("cuda_available") is True and fit.get("cpu_fallback") is False, fit.get("gpu", {}).get("device_name", ""))
    check(rows, "scvi_model_present", paths["model"].is_file(), str(paths["model"]))
    check(rows, "clustering_not_started", not paths["partitions"].exists(), str(paths["partitions"]))

    state = pd.read_csv(paths["source_state"])
    unchanged = int(state["unchanged"].astype(bool).sum())
    check(rows, "source_h5ads_unchanged", unchanged == state.shape[0], f"{unchanged}/{state.shape[0]}")

    resources = config["resources"]
    free_gib = shutil.disk_usage(resources["ssd_root"]).free / 2**30
    cluster_floor = float(resources["minimum_ssd_free_gib_by_stage"]["cluster"])
    original_floor = float(resources["minimum_ssd_free_gib"])
    n_cells = int(prepare["n_cells"])
    global_runs = len(config["clustering"]["global"]["resolutions"]) * len(config["clustering"]["global"]["seeds"])
    boundary_runs = len(config["clustering"]["boundary"]["resolutions"]) * len(config["clustering"]["boundary"]["seeds"])
    raw_partition_bytes = n_cells * (global_runs + boundary_runs) * 4 + n_cells * 16
    worst_case_gib = raw_partition_bytes / 2**30 + 2.0
    margin = free_gib - cluster_floor - worst_case_gib
    check(rows, "stage_specific_floor_preserves_original_fit_floor", float(resources["minimum_ssd_free_gib_by_stage"]["prepare"]) == original_floor and float(resources["minimum_ssd_free_gib_by_stage"]["fit"]) == original_floor, f"prepare/fit={original_floor:.1f} GiB")
    check(rows, "cluster_floor_met", free_gib >= cluster_floor, f"free={free_gib:.1f}; floor={cluster_floor:.1f} GiB")
    check(rows, "cluster_output_reserve_met", margin >= 50.0, f"post-reserve margin={margin:.1f} GiB")

    checks = pd.DataFrame(rows)
    result = "PASS_REVIEW_REQUIRED" if checks["status"].eq("PASS").all() else "FAIL"
    summary = {
        "result": result,
        "resource_contract_sha256": resource_contract_sha256(config),
        "resource_contract": resource_contract_payload(config),
        "recovery_contract_sha256": recovery_contract_sha256(config),
        "execution_config_sha256": sha256_file(CONFIG_PATH),
        "runner_sha256": sha256_file(RUNNER_PATH),
        "core_sha256": sha256_file(CORE_PATH),
        "observed_ssd_free_gib": free_gib,
        "original_minimum_ssd_free_gib": original_floor,
        "cluster_minimum_ssd_free_gib": cluster_floor,
        "worst_case_cluster_output_gib": worst_case_gib,
        "post_reserve_margin_gib": margin,
        "n_cells": int(prepare["n_cells"]),
        "n_hvgs": int(prepare["n_hvgs"]),
        "matrix_nnz": int(prepare["matrix_nnz"]),
        "primary_nk_anchors": int(prepare["primary_nk_anchors"]),
        "n_latent": int(fit["n_latent"]),
        "source_files": int(state.shape[0]),
        "unchanged_source_files": unchanged,
        "checks_passed": int(checks["status"].eq("PASS").sum()),
        "checks_total": int(checks.shape[0]),
        "clustering_performed": False,
        "classifier_fitting_performed": False,
        "scope": "Saved-latent RAPIDS clustering and pseudo-NK consensus QC only; no restaging, scVI refit, or classifier fitting.",
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
        "resource_contract_sha256": summary["resource_contract_sha256"],
        "review_artifact": str(STATIC_DIR / "gdtai_v4_2_cluster_resource_preflight_report.pdf"),
        "scope": summary["scope"],
    }
    (LOG_DIR / "CLUSTER_EXECUTION_APPROVAL_TEMPLATE.json").write_text(json.dumps(approval, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    render(checks, summary)
    print(json.dumps({key: value for key, value in summary.items() if key not in {"resource_contract"}}, indent=2, sort_keys=True))
    if result == "FAIL":
        raise SystemExit("V4.2 cluster resource-amendment preflight failed")


if __name__ == "__main__":
    main()
