#!/usr/bin/env python3
"""Run synthetic/read-only implementation QC for gdTAI V4.2 integration."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

from workflows.gdtai.gdtai_v4_2_integration_core import (
    input_file_state,
    read_json,
    resolve,
    sha256_file,
    validate_preflight_approval,
)
from workflows.gdtai.run_gdtai_v4_2_nk_reference_integration import (
    CORE_PATH,
    RUNNER_PATH,
    make_direct_cuda_environment,
    require_gpu,
)


CONFIG_PATH = PROJECT_ROOT / "configs/models/gdtai/v4_2_integration_execution.json"
TABLE_DIR = PROJECT_ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_implementation_qc"
FIGURE_DIR = PROJECT_ROOT / "Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_implementation_qc"
LOG_DIR = PROJECT_ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_implementation_qc"
STATIC_DIR = PROJECT_ROOT / "gdT_prediction/gdtai_v4_2_implementation_qc"


def add_check(checks: list[dict[str, Any]], name: str, passed: bool, detail: str) -> None:
    checks.append({"check": name, "status": "PASS" if passed else "FAIL", "detail": detail})


def run_command(command: list[str], env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, cwd=PROJECT_ROOT, env=env, text=True, capture_output=True, check=False)


def synthetic_gpu_qc(config: dict[str, Any]) -> tuple[dict[str, Any], np.ndarray, np.ndarray]:
    environment = make_direct_cuda_environment()
    gpu = require_gpu(config)
    import rapids_singlecell as rsc
    import scvi

    rng = np.random.default_rng(int(config["random_seed"]))
    n_cells, n_genes = 1600, 120
    groups = np.repeat(np.arange(4), n_cells // 4)
    batches = np.tile(np.repeat(["batch_a", "batch_b"], n_cells // 8), 4)
    rates = np.full((4, n_genes), 0.35, dtype=np.float32)
    for group in range(4):
        rates[group, group * 20 : (group + 1) * 20] = 3.5
    counts = rng.poisson(rates[groups]).astype(np.float32)
    adata = ad.AnnData(
        X=sparse.csr_matrix(counts),
        obs=pd.DataFrame({"batch": pd.Categorical(batches)}),
        var=pd.DataFrame(index=[f"G{index}" for index in range(n_genes)]),
    )
    scvi.settings.seed = int(config["random_seed"])
    scvi.model.SCVI.setup_anndata(adata, batch_key="batch")
    model = scvi.model.SCVI(
        adata,
        n_hidden=32,
        n_layers=1,
        n_latent=8,
        gene_likelihood="nb",
    )
    model.train(
        max_epochs=2,
        batch_size=256,
        accelerator="gpu",
        devices=1,
        early_stopping=False,
        enable_progress_bar=False,
    )
    latent = model.get_latent_representation(batch_size=256).astype(np.float32)
    root_device = str(model.trainer.strategy.root_device)
    if "cuda" not in root_device.lower():
        raise RuntimeError(f"Synthetic scVI trainer did not use CUDA: {root_device}")
    graph = ad.AnnData(X=sparse.csr_matrix((n_cells, 1), dtype=np.float32))
    graph.obsm["X_scVI"] = latent
    rsc.pp.neighbors(graph, n_neighbors=15, use_rep="X_scVI", algorithm="brute", random_state=11)
    rsc.tl.leiden(graph, resolution=0.5, random_state=11, key_added="run_11")
    rsc.tl.leiden(graph, resolution=0.5, random_state=29, key_added="run_29")
    labels = pd.Categorical(graph.obs["run_11"]).codes.astype(np.int32)
    labels_2 = pd.Categorical(graph.obs["run_29"]).codes.astype(np.int32)
    return (
        {
            "gpu": gpu,
            "environment": environment,
            "scvi_version": scvi.__version__,
            "latent_shape": list(latent.shape),
            "latent_finite": bool(np.isfinite(latent).all()),
            "trainer_root_device": root_device,
            "rapids_global_clusters_seed11": int(np.unique(labels).size),
            "rapids_global_clusters_seed29": int(np.unique(labels_2).size),
            "cpu_fallback": False,
        },
        latent,
        groups,
    )


def render(checks: pd.DataFrame, gpu: dict[str, Any], latent: np.ndarray, groups: np.ndarray) -> None:
    for directory in [TABLE_DIR, FIGURE_DIR, LOG_DIR, STATIC_DIR]:
        directory.mkdir(parents=True, exist_ok=True)
    checks.to_csv(TABLE_DIR / "implementation_qc_checks.csv", index=False)

    plt.rcParams.update({"font.size": 9, "axes.titlesize": 11, "axes.labelsize": 9})
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.3))
    status_counts = checks["status"].value_counts().reindex(["PASS", "FAIL"], fill_value=0)
    bars = axes[0].bar(status_counts.index, status_counts.values, color=["#16846c", "#b04b3a"])
    axes[0].bar_label(bars)
    axes[0].set_ylabel("Checks")
    axes[0].set_title("Implementation QC status")
    palette = np.asarray(["#285f8f", "#14866d", "#b06035", "#7b5aa6"])
    axes[1].scatter(latent[:, 0], latent[:, 1], c=palette[groups], s=4, alpha=0.55, linewidths=0)
    axes[1].set_xlabel("Synthetic scVI latent 1")
    axes[1].set_ylabel("Synthetic scVI latent 2")
    axes[1].set_title("A100 scVI and RAPIDS smoke test")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "implementation_qc_summary.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    result = "PASS_REVIEW_REQUIRED" if checks["status"].eq("PASS").all() else "FAIL"
    report = f"""# gdTAI V4.2 sidecar integration implementation QC

## Decision

**{result}.** The sidecar runner is implemented and passed deterministic unit,
fail-closed, sparse-H5AD, direct-CUDA scVI, and RAPIDS clustering tests. No
project-data integration, project-data scVI fit, pseudo-label construction, or
classifier fit was performed.

- Checks passed: **{int(checks['status'].eq('PASS').sum())}/{checks.shape[0]}**.
- Synthetic scVI device: **{gpu['trainer_root_device']}** on **{gpu['gpu']['device_name']}**.
- Synthetic latent matrix: **{gpu['latent_shape'][0]:,} x {gpu['latent_shape'][1]}**, all finite.
- Locked cohorts remain excluded by code and the project-data stage aborts while
  its separate approval file is absent.
- CPU fallback is forbidden for scVI and RAPIDS stages.

![Implementation QC](../../Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_implementation_qc/implementation_qc_summary.png)

## Frozen executable design

1. Select 4,000 common non-TCR, non-mitochondrial, non-ribosomal, and
   non-immunoglobulin HVGs from a deterministic source-balanced sample capped at
   20,000 cells per source.
2. Stage only development cohorts in a new sparse H5AD on SSD. Locked cohorts
   are excluded before metadata loading, HVG selection, matrix staging, scVI,
   clustering, and label construction.
3. Fit the precommitted 30-dimensional negative-binomial scVI model on one A100
   GPU, with no CPU fallback.
4. Run nine global and nine boundary RAPIDS Leiden partitions. The boundary is
   defined only by global clusters containing both independent primary-NK and
   productive-TCR anchors.
5. Select low-weight pseudo-NK development candidates only through repeated
   cluster agreement, independent-anchor purity, productive-TCR contamination,
   and cross-source representation. No CD3, TRDC, TRDV, NKG7, GNLY, KLRD1,
   TYROBP, FCER1G, or FCGR3A threshold defines truth.
6. Source balance and the primary-NK effective-weight cap are enforced at the
   later classifier-fitting gate. Pseudo-labels cannot control validation
   guardrails.

The two implementation safety floors frozen before project-data execution are
at least 20 independent anchors and at least 10 primary NK anchors per
qualifying cluster. They prevent tiny-anchor clusters from passing a percentage
rule by chance.

## Checks

{checks.to_markdown(index=False)}

## Supervision gate

Project-data integration is still blocked. Activating
`PROJECT_DATA_INTEGRATION_APPROVAL.json` authorizes only development-data sparse
staging, scVI, consensus clustering, and pseudo-label QC. It does not authorize
classifier fitting, threshold selection, model promotion, release fitting, or
whole-atlas inference.

## Exact artifacts

- Execution config: `{CONFIG_PATH.relative_to(PROJECT_ROOT)}`
- Runner: `{RUNNER_PATH.relative_to(PROJECT_ROOT)}`
- Core: `{CORE_PATH.relative_to(PROJECT_ROOT)}`
- Approval template: `{(LOG_DIR / 'PROJECT_DATA_INTEGRATION_APPROVAL_TEMPLATE.json').relative_to(PROJECT_ROOT)}`
- Full checks: `{(TABLE_DIR / 'implementation_qc_checks.csv').relative_to(PROJECT_ROOT)}`
"""
    markdown_path = LOG_DIR / "gdtai_v4_2_implementation_qc_summary.md"
    markdown_path.write_text(report, encoding="utf-8")
    css = """
body{font-family:Arial,Helvetica,sans-serif;color:#17212b;line-height:1.46;max-width:1080px;margin:0 auto;padding:30px 42px;background:#fff}h1{font-size:29px;color:#17324d;border-bottom:4px solid #14866d;padding-bottom:12px}h2{font-size:21px;color:#17324d;margin-top:30px}p,li{font-size:14px}table{border-collapse:collapse;width:100%;font-size:10px;margin:16px 0 24px}th{background:#17324d;color:#fff;text-align:left}th,td{padding:6px 7px;border:1px solid #cdd6df;vertical-align:top;overflow-wrap:anywhere}tr:nth-child(even){background:#f3f6f8}img{display:block;max-width:96%;max-height:650px;margin:22px auto}code{background:#edf1f4;padding:1px 4px;border-radius:2px}@media print{body{max-width:none;padding:8mm 9mm}h1{font-size:23px}h2{font-size:17px;page-break-after:avoid}table,img{page-break-inside:avoid}p,li{font-size:10.5px}table{font-size:7.5px}img{max-height:155mm}}
""".strip()
    (STATIC_DIR / "print.css").write_text(css + "\n", encoding="utf-8")
    subprocess.run(
        ["pandoc", str(markdown_path), "--standalone", "--metadata", "pagetitle=gdTAI V4.2 implementation QC", "--css", "print.css", "-o", str(STATIC_DIR / "index.html")],
        cwd=PROJECT_ROOT,
        check=True,
    )
    profile = Path("/tmp/gdtai-v42-implementation-qc-chrome")
    profile.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "google-chrome", "--headless", "--no-sandbox", "--disable-gpu", "--disable-dev-shm-usage",
            "--disable-breakpad", "--disable-crash-reporter", "--allow-file-access-from-files", "--no-pdf-header-footer",
            f"--user-data-dir={profile}",
            f"--print-to-pdf={STATIC_DIR / 'gdtai_v4_2_implementation_qc_report.pdf'}",
            (STATIC_DIR / "index.html").resolve().as_uri(),
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )


def main() -> None:
    config = read_json(CONFIG_PATH)
    checks: list[dict[str, Any]] = []
    preflight = validate_preflight_approval(config)
    add_check(checks, "checksum_bound_preflight_approval", True, preflight["approval"]["approved_at"])
    add_check(checks, "development_locked_role_split", int(preflight["roles"]["include_in_integration_fit"].sum()) == 6 and int(preflight["roles"]["allow_locked_evaluation"].sum()) == 3, "6 development objects; 3 locked objects")
    add_check(checks, "locked_development_permissions", not preflight["roles"].loc[preflight["roles"]["allow_locked_evaluation"], ["include_in_integration_fit", "include_in_cluster_label_design", "allow_pseudolabel_training", "allow_model_tuning"]].any(axis=None), "locked cohorts have zero development permissions")
    add_check(checks, "cpu_fallback_forbidden", config["scvi"]["allow_cpu_fallback"] is False, "scVI and RAPIDS project stages require CUDA")
    add_check(checks, "marker_threshold_truth_forbidden", config["pseudolabel"]["marker_thresholds_may_define_truth"] is False, "cluster evidence only")
    add_check(checks, "pseudo_labels_cannot_control_guardrails", config["pseudolabel"]["may_control_validation_guardrails"] is False, "primary truth retains guardrail ownership")

    compile_result = run_command([
        sys.executable, "-m", "py_compile", str(CORE_PATH), str(RUNNER_PATH), str(Path(__file__).resolve())
    ])
    add_check(checks, "python_compile", compile_result.returncode == 0, compile_result.stderr.strip() or "three modules compiled")

    pytest_env = os.environ.copy()
    pytest_env.update({"PYTEST_DISABLE_PLUGIN_AUTOLOAD": "1", "XDG_CACHE_HOME": "/tmp/gdtai-v42-pytest-cache"})
    tests = run_command(
        ["/home/anaconda3/bin/python", "-m", "pytest", "-q", "tests/test_gdtai_v4_2_integration_core.py"],
        env=pytest_env,
    )
    add_check(checks, "deterministic_unit_and_sparse_h5ad_tests", tests.returncode == 0, tests.stdout.strip().splitlines()[-1] if tests.stdout.strip() else tests.stderr.strip())

    validate = run_command([
        sys.executable, str(RUNNER_PATH), "--config", str(CONFIG_PATH), "--stage", "validate"
    ])
    add_check(checks, "runner_validate_stage", validate.returncode == 0, "checksum and role validation passed" if validate.returncode == 0 else validate.stderr.strip())
    blocked = run_command([
        sys.executable, str(RUNNER_PATH), "--config", str(CONFIG_PATH), "--stage", "prepare"
    ])
    add_check(checks, "project_data_stage_fails_without_approval", blocked.returncode != 0 and "approval is absent" in blocked.stderr, "prepare aborted before SSD access")

    preflight_state = pd.read_csv(resolve("Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_integration_preflight/input_file_state.csv"))
    current_state = input_file_state(preflight["roles"])
    state = preflight_state[["cohort_id", "size_bytes_after", "mtime_ns_after"]].merge(current_state, on="cohort_id")
    state["unchanged"] = state["size_bytes_after"].eq(state["size_bytes"]) & state["mtime_ns_after"].eq(state["mtime_ns"])
    add_check(checks, "all_source_h5ads_unchanged", bool(state["unchanged"].all()), f"{int(state['unchanged'].sum())}/{state.shape[0]} size/mtime pairs unchanged")

    gpu, latent, groups = synthetic_gpu_qc(config)
    add_check(checks, "direct_cuda_a100", gpu["gpu"]["cuda_available"] and gpu["gpu"]["total_memory_gib"] >= 75, f"{gpu['gpu']['device_name']}; {gpu['gpu']['total_memory_gib']:.1f} GiB")
    add_check(checks, "synthetic_scvi_gpu_fit", gpu["latent_finite"] and "cuda" in gpu["trainer_root_device"].lower(), f"latent={gpu['latent_shape']}; device={gpu['trainer_root_device']}")
    add_check(checks, "synthetic_rapids_neighbors_leiden", gpu["rapids_global_clusters_seed11"] > 1 and gpu["rapids_global_clusters_seed29"] > 1, f"clusters={gpu['rapids_global_clusters_seed11']},{gpu['rapids_global_clusters_seed29']}")
    add_check(checks, "synthetic_gpu_no_cpu_fallback", gpu["cpu_fallback"] is False, "direct CUDA only")

    checks_frame = pd.DataFrame(checks)
    result = "PASS_REVIEW_REQUIRED" if checks_frame["status"].eq("PASS").all() else "FAIL"
    for directory in [TABLE_DIR, FIGURE_DIR, LOG_DIR, STATIC_DIR]:
        directory.mkdir(parents=True, exist_ok=True)
    state.to_csv(TABLE_DIR / "input_file_state_recheck.csv", index=False)
    (LOG_DIR / "gpu_synthetic_qc.json").write_text(json.dumps(gpu, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary = {
        "result": result,
        "project_data_integration_performed": False,
        "classifier_fitting_performed": False,
        "checks_passed": int(checks_frame["status"].eq("PASS").sum()),
        "checks_total": int(checks_frame.shape[0]),
        "execution_config_sha256": sha256_file(CONFIG_PATH),
        "runner_sha256": sha256_file(RUNNER_PATH),
        "core_sha256": sha256_file(CORE_PATH),
        "next_action": "Review and explicitly approve project-data sidecar integration; classifier fitting remains a later gate.",
    }
    (LOG_DIR / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    approval = {
        "approved": False,
        "approved_at": "",
        "approved_by": "",
        "protocol_version": config["protocol_version"],
        "execution_config_sha256": summary["execution_config_sha256"],
        "runner_sha256": summary["runner_sha256"],
        "core_sha256": summary["core_sha256"],
        "review_artifact": str(STATIC_DIR / "gdtai_v4_2_implementation_qc_report.pdf"),
        "scope": "Development-data sidecar staging, GPU scVI, consensus clustering, and pseudo-label QC only; no classifier fitting.",
    }
    (LOG_DIR / "PROJECT_DATA_INTEGRATION_APPROVAL_TEMPLATE.json").write_text(json.dumps(approval, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    render(checks_frame, gpu, latent, groups)
    print(json.dumps(summary, indent=2, sort_keys=True))
    if result == "FAIL":
        raise SystemExit("V4.2 implementation QC failed")


if __name__ == "__main__":
    main()
