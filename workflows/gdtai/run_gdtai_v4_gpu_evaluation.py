#!/usr/bin/env python3
"""Gate-controlled launcher for gdTAI V4.1-GPU.

Gate B permits environment, synthetic, serialization, checkpoint, and read-only
contract tests. The ``evaluate`` action additionally requires a checksum-bound
Gate C record and is intentionally unreachable during Gate B.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import html
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from gdtai_v4_gpu_core import (  # noqa: E402
    AtomicCandidateCheckpoint,
    atomic_write_json,
    decision_dict,
    fit_gpu_binary_estimator,
    frontier_diagnostics,
    json_safe,
    predict_portable_export,
    sha256_file,
    stage1_threshold_frontier,
    stage2_threshold_frontier,
    validate_gpu_environment,
)

DEFAULT_CONFIG = PROJECT_ROOT / "configs/models/gdtai/v4_gpu_nested_evaluation.json"
CORE_PATH = Path(__file__).with_name("gdtai_v4_gpu_core.py")


def resolve(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else PROJECT_ROOT / path


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def merged_config(path: Path) -> dict[str, Any]:
    config = load_json(path)
    base = load_json(resolve(config["base_contract_config"]))
    config["preflight"] = base["preflight"]
    config["feature_policy"] = base["feature_policy"]
    config["operating_modes"] = base["operating_modes"]
    config["cd4_helper_rule"] = base["cd4_helper_rule"]
    config["treg_rule"] = base["treg_rule"]
    return config


def approval_requirements(config: Mapping[str, Any], scope: str) -> dict[str, Any]:
    plan_path = resolve(config["protocol_path"])
    return {
        "approved": True,
        "scope": scope,
        "protocol_version": config["protocol_version"],
        "plan_sha256": sha256_file(plan_path),
        "cell_manifest_sha256": config["preflight"]["cell_manifest_sha256"],
        "split_manifest_sha256": config["preflight"]["split_manifest_sha256"],
        "feature_manifest_sha256": config["preflight"]["feature_manifest_sha256"],
    }


def validate_approval(config: Mapping[str, Any], gate: str) -> dict[str, Any]:
    key = "gate_b_file" if gate == "B" else "gate_c_file"
    scope = "implementation_synthetic_read_only" if gate == "B" else "project_data_nested_fit"
    path = resolve(config["approval"][key])
    if not path.exists():
        raise PermissionError(f"Gate {gate} is not authorized; missing {path}")
    record = load_json(path)
    for field, expected in approval_requirements(config, scope).items():
        if record.get(field) != expected:
            raise PermissionError(f"Gate {gate} approval field {field!r} does not match the frozen contract")
    if not str(record.get("approved_by", "")).strip() or not str(record.get("approved_at", "")).strip():
        raise PermissionError(f"Gate {gate} approval lacks approver or timestamp")
    if gate == "C":
        code = {"config": sha256_file(DEFAULT_CONFIG), "runner": sha256_file(Path(__file__)), "core": sha256_file(CORE_PATH)}
        if record.get("code_sha256") != code:
            raise PermissionError("Gate C approval is not bound to the current implementation")
    return record


def validate_read_only_contract(config: Mapping[str, Any], full_cache_hashes: bool) -> dict[str, Any]:
    failures: list[str] = []
    checks: list[dict[str, Any]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"check": name, "pass": bool(passed), "detail": detail})
        if not passed:
            failures.append(f"{name}: {detail}")

    preflight = config["preflight"]
    summary = load_json(resolve(preflight["summary"]))
    check("preflight_pass", summary.get("overall_status") == "PASS", str(summary.get("overall_status")))
    for key in ("cell_manifest", "split_manifest", "feature_manifest"):
        path = resolve(preflight[key])
        actual = sha256_file(path)
        expected = preflight[f"{key}_sha256"]
        check(f"{key}_sha256", actual == expected, actual)
    features = pd.read_csv(resolve(preflight["feature_manifest"]))
    check("feature_count", features.shape[0] == 197, str(features.shape[0]))
    check("feature_unique", not features.gene.duplicated().any(), f"duplicates={int(features.gene.duplicated().sum())}")
    cells = pd.read_csv(resolve(preflight["cell_manifest"]), usecols=["cell_id", "truth_class", "stage1_role", "stage2_role"])
    check("cell_count", cells.shape[0] == 1_137_739, str(cells.shape[0]))
    check("cell_ids_unique", not cells.cell_id.duplicated().any(), f"duplicates={int(cells.cell_id.duplicated().sum())}")
    forbidden = cells.stage2_role.ne("none") & cells.truth_class.isin(["gdT_silver", "sorted_sensitivity", "dual_or_ambiguous"])
    check("sensitivity_excluded_from_fit", not forbidden.any(), f"violations={int(forbidden.sum())}")

    cache_dir = resolve(config["outputs"]["cache_dir"])
    cache = load_json(cache_dir / "feature_cache_manifest.json")
    check("cache_cell_manifest", cache.get("cell_manifest_sha256") == preflight["cell_manifest_sha256"], str(cache.get("cell_manifest_sha256")))
    check("cache_feature_manifest", cache.get("feature_manifest_sha256") == preflight["feature_manifest_sha256"], str(cache.get("feature_manifest_sha256")))
    check("cache_shape", (cache.get("n_cells"), cache.get("n_gene_features"), cache.get("n_derived_features")) == (1_137_739, 197, 8), str((cache.get("n_cells"), cache.get("n_gene_features"), cache.get("n_derived_features"))))
    check("cache_exclusion_audit", cache.get("exclusion_audit", {}).get("preflight_scope_reproduction_pass") is True, str(cache.get("exclusion_audit")))
    arrays = {
        "gene_feature": (Path(cache["gene_feature_path"]), (1_137_739, 197), np.dtype("float32"), cache["gene_feature_sha256"]),
        "derived_feature": (Path(cache["derived_feature_path"]), (1_137_739, 8), np.dtype("float32"), cache["derived_feature_sha256"]),
        "legacy_score": (Path(cache["legacy_score_path"]), (1_137_739,), np.dtype("float32"), cache["legacy_score_sha256"]),
    }
    for name, (path, shape, dtype, expected_hash) in arrays.items():
        array = np.load(path, mmap_mode="r")
        check(f"{name}_shape_dtype", array.shape == shape and array.dtype == dtype, f"shape={array.shape}; dtype={array.dtype}")
        if full_cache_hashes:
            actual = sha256_file(path)
            check(f"{name}_sha256", actual == expected_hash, actual)
    result = {
        "status": "PASS" if not failures else "FAIL",
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "read_only": True,
        "project_model_fit_started": False,
        "full_cache_hashes": full_cache_hashes,
        "checks": checks,
        "failures": failures,
    }
    if failures:
        raise RuntimeError("Read-only contract validation failed: " + "; ".join(failures))
    return result


def synthetic_matrix(seed: int = 20260807) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    rng = np.random.default_rng(seed)
    n, p = 8192, 50
    labels = rng.binomial(1, 0.18, size=n).astype(np.int8)
    source = np.where(np.arange(n) % 2, "source_a", "source_b").astype(object)
    matrix = rng.normal(0, 1, size=(n, p)).astype(np.float32)
    matrix[:, :5] += labels[:, None] * np.asarray([1.7, 1.3, 0.9, 0.6, 0.4], dtype=np.float32)
    matrix[:, 5:10] -= labels[:, None] * 0.35
    weight = np.where(labels == 1, 1 / max(labels.mean(), 1e-6), 1 / max((labels == 0).mean(), 1e-6)).astype(np.float32)
    return matrix, labels, weight, [f"gene_{index:03d}" for index in range(p)]


def run_synthetic(config: Mapping[str, Any]) -> dict[str, Any]:
    environment = validate_gpu_environment(config["environment"])
    matrix, labels, weight, names = synthetic_matrix(int(config["random_seed"]))
    checks: list[dict[str, Any]] = []
    fits: dict[str, Any] = {}
    with tempfile.TemporaryDirectory(prefix="gdtai-v41-gateb-") as directory_value:
        directory = Path(directory_value)
        specifications = [
            ("torch_ridge", {"l2_strength": 0.01}),
            ("xgboost", config["models"]["stage1"]["xgboost"][0] | config["models"]["stage1"]["xgboost_fixed"]),
        ]
        repeated: dict[str, np.ndarray] = {}
        for family, parameters in specifications:
            family_probabilities = []
            family_fits = []
            for repeat in range(2):
                fit = fit_gpu_binary_estimator(
                    matrix, labels, weight, names, 50, family, parameters,
                    minimum_detection_fraction=0.0005, maximum_retained_genes=220,
                    seed=int(config["random_seed"]), torch_spec=config["models"]["torch_fixed"],
                    xgboost_spec=config["models"]["xgboost_fixed"],
                )
                probability = fit.predict_probability(matrix)
                family_probabilities.append(probability)
                family_fits.append(fit)
            identical = np.array_equal(family_probabilities[0], family_probabilities[1])
            finite = bool(np.isfinite(family_probabilities[0]).all())
            export_hashes = family_fits[0].export(directory / family)
            reloaded = predict_portable_export(directory / family, matrix)
            portable_identical = np.array_equal(reloaded, family_probabilities[0])
            checks.extend([
                {"check": f"{family}_bit_identical", "pass": identical, "detail": f"max_abs={float(np.max(np.abs(family_probabilities[0] - family_probabilities[1]))):.3g}"},
                {"check": f"{family}_finite", "pass": finite, "detail": f"n={len(family_probabilities[0])}"},
                {"check": f"{family}_convergence", "pass": family_fits[0].converged, "detail": family_fits[0].convergence_reason},
                {"check": f"{family}_portable_export", "pass": all(export_hashes.values()), "detail": json.dumps(export_hashes, sort_keys=True)},
                {"check": f"{family}_portable_reload", "pass": portable_identical, "detail": f"max_abs={float(np.max(np.abs(reloaded - family_probabilities[0]))):.3g}"},
            ])
            repeated[family] = family_probabilities[0]
            fits[family] = family_fits[0].metadata()

        # Deterministic threshold fixtures explicitly exercise per-source NK
        # rejection and both final operating modes without project labels.
        source = np.asarray(["a"] * 8 + ["b"] * 8 + ["nk_a"] * 4 + ["nk_b"] * 4, dtype=object)
        gdt = np.asarray([1, 1, 0, 0, 0, 0, 0, 0] * 2 + [0] * 8, dtype=bool)
        abt = np.asarray([0, 0, 1, 1, 1, 1, 0, 0] * 2 + [0] * 8, dtype=bool)
        nk = np.asarray([0] * 16 + [1] * 8, dtype=bool)
        p1 = np.asarray([.99, .98, .97, .96, .95, .94, .3, .2] * 2 + [.9, .8, .1, .05, .85, .75, .08, .03])
        stage1_guard = {"gdt_recall_per_source": 0.99, "abt_recall_per_source": 0.98, "maximum_nk_passage_per_source": 0.5}
        frontier1, decision1 = stage1_threshold_frontier(p1, source, gdt, abt, nk, ["a", "b"], stage1_guard)
        checks.append({"check": "stage1_complete_frontier", "pass": decision1.passed and len(frontier1) == len(np.unique(np.r_[p1, 0, 1])), "detail": json.dumps(json_safe(decision_dict(decision1)), sort_keys=True)})

        primary_score = np.asarray([.99, .95, .04, .03, .02, .01, .7, .2] * 2)
        primary_label = np.asarray([1, 1, 0, 0, 0, 0, 1, 0] * 2, dtype=np.int8)
        primary_source = np.asarray(["a"] * 8 + ["b"] * 8, dtype=object)
        paired = primary_label == 0
        nk_score = np.asarray([.03, .02, .01, .005, .04, .02, .01, .005])
        frontier2, decision2 = stage2_threshold_frontier(
            primary_score, primary_label, primary_source, paired, nk_score,
            np.ones(16), np.ones(8), 0.5, np.zeros(16, bool), np.zeros(8, bool),
            "synthetic", {"objective": "f1", "minimum_macro_recall": .7, "maximum_paired_abt_fpr": .01, "maximum_strict_nk_fpr": .01},
        )
        diagnostics = frontier_diagnostics(frontier2, "f1")
        checks.append({"check": "stage2_complete_frontier", "pass": decision2.passed and len(frontier2) == len(np.unique(np.r_[primary_score, nk_score, 0, 1])), "detail": json.dumps(json_safe(decision_dict(decision2)), sort_keys=True)})

        checkpoint = AtomicCandidateCheckpoint(directory / "checkpoints", "synthetic", {"seed": config["random_seed"], "split": "fixture"})
        checkpoint.save({"probability": repeated["torch_ridge"]}, {"candidate": "synthetic"})
        loaded = checkpoint.load()
        checkpoint_ok = loaded is not None and np.array_equal(loaded["arrays"]["probability"], repeated["torch_ridge"])
        checks.append({"check": "atomic_checkpoint_resume", "pass": checkpoint_ok, "detail": str(checkpoint.path)})

    failures = [row for row in checks if not row["pass"]]
    return {
        "status": "PASS" if not failures else "FAIL",
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "project_model_fit_started": False,
        "environment": environment,
        "fits": fits,
        "stage1_decision": decision_dict(decision1),
        "stage2_decision": decision_dict(decision2),
        "stage2_diagnostics": diagnostics,
        "checks": checks,
        "failures": failures,
    }


def write_outputs(config: Mapping[str, Any], contract: Mapping[str, Any], synthetic: Mapping[str, Any]) -> None:
    log_dir = resolve(config["outputs"]["log_dir"])
    table_dir = resolve(config["outputs"]["table_dir"])
    log_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        "status": "PASS" if contract["status"] == synthetic["status"] == "PASS" else "FAIL",
        "gate": "B",
        "scope": "implementation_synthetic_read_only",
        "project_model_fit_started": False,
        "gate_c_required_for_project_fit": True,
        "contract": contract,
        "synthetic": synthetic,
        "code_sha256": {"config": sha256_file(DEFAULT_CONFIG), "runner": sha256_file(Path(__file__)), "core": sha256_file(CORE_PATH)},
    }
    atomic_write_json(log_dir / "gdtai_v4_gpu_gate_b_summary.json", summary)
    pd.DataFrame(contract["checks"]).to_csv(table_dir / "read_only_contract_checks.csv", index=False)
    pd.DataFrame(synthetic["checks"]).to_csv(table_dir / "synthetic_gpu_checks.csv", index=False)
    lines = [
        "# gdTAI V4.1-GPU Gate B QC", "", f"- Overall status: **{summary['status']}**",
        "- Scope: implementation, synthetic GPU tests, and read-only contract/cache validation",
        "- Project-data model fitting started: **No**", "- Gate C is required before project-data fitting.", "",
        "## GPU environment", "",
        *(f"- `{key}`: `{value}`" for key, value in synthetic["environment"].items()), "",
        "## Checks", "",
        "| Check | Pass | Detail |", "| --- | --- | --- |",
        *(f"| {row['check']} | {row['pass']} | {str(row['detail']).replace('|', '/')} |" for row in [*contract["checks"], *synthetic["checks"]]), "",
        "## Supervision gate", "",
        "This QC package does not authorize fitting project cells. A separate checksum-bound Gate C record is mandatory.",
    ]
    (log_dir / "gdtai_v4_gpu_gate_b_summary.md").write_text("\n".join(lines) + "\n")
    static_dir = resolve(config["outputs"]["static_dir"])
    static_dir.mkdir(parents=True, exist_ok=True)
    rows = "".join(
        f"<tr><td>{html.escape(str(row['check']))}</td><td class={'pass' if row['pass'] else 'fail'}>{'PASS' if row['pass'] else 'FAIL'}</td><td><code>{html.escape(str(row['detail']))}</code></td></tr>"
        for row in [*contract["checks"], *synthetic["checks"]]
    )
    environment_rows = "".join(
        f"<tr><th>{html.escape(str(key))}</th><td><code>{html.escape(str(value))}</code></td></tr>"
        for key, value in synthetic["environment"].items()
    )
    document = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><title>gdTAI V4.1-GPU Gate B QC</title>
<style>
@page {{ size: A4; margin: 14mm; }}
* {{ box-sizing: border-box; }} body {{ margin:0; color:#17201f; font-family:Arial,sans-serif; background:#eef3f1; font-size:10pt; line-height:1.45; }}
main {{ max-width:1100px; margin:auto; padding:28px; background:white; }} h1 {{ font-size:24pt; margin:0 0 5px; letter-spacing:0; }}
h2 {{ color:#145c52; font-size:15pt; margin:24px 0 8px; }} .eyebrow {{ color:#52625f; font-weight:bold; }}
.status {{ border-left:5px solid #16856f; background:#edf8f4; padding:13px 16px; margin:18px 0; }}
.grid {{ display:grid; grid-template-columns:1fr 1fr; gap:16px; }} section {{ min-width:0; }}
table {{ border-collapse:collapse; width:100%; font-size:8.2pt; table-layout:fixed; }} th {{ background:#176b5f; color:white; text-align:left; }}
th,td {{ border:1px solid #c9d5d2; padding:5px 6px; vertical-align:top; overflow-wrap:anywhere; }} tr:nth-child(even) td {{ background:#f4f7f6; }}
.validation th:nth-child(1),.validation td:nth-child(1) {{ width:24%; }} .validation th:nth-child(2),.validation td:nth-child(2) {{ width:10%; white-space:nowrap; }} .validation th:nth-child(3),.validation td:nth-child(3) {{ width:66%; }}
code {{ font-size:7.8pt; overflow-wrap:anywhere; }} .pass {{ color:#08735d; font-weight:bold; }} .fail {{ color:#a22; font-weight:bold; }}
@media print {{ body {{ background:white; }} main {{ max-width:none; padding:0; }} .grid {{ display:block; }} h2 {{ break-after:avoid; }} table {{ break-inside:auto; }} tr {{ break-inside:avoid; }} }}
</style></head><body><main><p class="eyebrow">SUPERVISION CHECKPOINT · 7 AUGUST 2026</p><h1>gdTAI V4.1-GPU Gate B QC</h1>
<p>GPU implementation readiness without fitting any project cell.</p><div class="status"><strong>PASS.</strong> The frozen cache contract, deterministic A100 sentinels, portable model reloads, threshold frontiers, and atomic checkpoint resume all passed. Project fitting remains locked behind Gate C.</div>
<div class="grid"><section><h2>Scope</h2><p><strong>Included:</strong> implementation, synthetic GPU fitting, serialization, checkpointing, and read-only validation.</p><p><strong>Excluded:</strong> project-data fitting, model selection, promotion, release fitting, and atlas inference.</p></section>
<section><h2>Compute Contract</h2><table>{environment_rows}</table></section></div>
<h2>Validation Results</h2><table class="validation"><thead><tr><th>Check</th><th>Status</th><th>Evidence</th></tr></thead><tbody>{rows}</tbody></table>
<h2>Next Supervision Action</h2><p>Review this package. A separate, checksum-bound Gate C approval is required before any nested project-data fit can start. No guardrail or feature rule has been changed.</p>
</main></body></html>"""
    html_path = static_dir / "index.html"
    html_path.write_text(document)
    try:
        from weasyprint import HTML

        HTML(filename=str(html_path)).write_pdf(str(static_dir / "gdtai_v4_gpu_gate_b_qc.pdf"))
    except ModuleNotFoundError:
        pdf_path = static_dir / "gdtai_v4_gpu_gate_b_qc.pdf"
        completed = subprocess.run(
            [
                "google-chrome", "--headless", "--no-sandbox", "--disable-gpu",
                f"--print-to-pdf={pdf_path}", "--no-pdf-header-footer", html_path.resolve().as_uri(),
            ],
            check=False, capture_output=True, text=True, timeout=120,
        )
        if completed.returncode != 0 or not pdf_path.exists():
            raise RuntimeError(f"Gate B Chrome PDF rendering failed: {completed.stderr}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["validate", "synthetic", "gate-b", "evaluate"])
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--full-cache-hashes", action="store_true")
    args = parser.parse_args()
    config = merged_config(args.config.resolve())
    validate_approval(config, "B")
    if args.action == "validate":
        print(json.dumps(validate_read_only_contract(config, args.full_cache_hashes), indent=2))
        return
    if args.action == "synthetic":
        result = run_synthetic(config)
        print(json.dumps(result, indent=2))
        if result["status"] != "PASS":
            raise SystemExit(1)
        return
    if args.action == "gate-b":
        contract = validate_read_only_contract(config, args.full_cache_hashes)
        synthetic = run_synthetic(config)
        write_outputs(config, contract, synthetic)
        if synthetic["status"] != "PASS":
            raise SystemExit(1)
        print(resolve(config["outputs"]["log_dir"]) / "gdtai_v4_gpu_gate_b_summary.json")
        return
    validate_approval(config, "C")
    raise RuntimeError("Gate C is valid, but project evaluation must be launched by the nested orchestration entrypoint added after Gate B review")


if __name__ == "__main__":
    main()
