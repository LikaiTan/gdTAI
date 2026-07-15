#!/usr/bin/env python3
"""Verify registries, migrations, links, models, and unchanged H5AD baselines."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from tnk_atlas.model_io import load_trusted_pickle  # noqa: E402
from tnk_atlas.paths import ProjectPaths  # noqa: E402
from tnk_atlas.provenance import sampled_sha256, sha256_file  # noqa: E402
from tnk_atlas.registry import read_csv, validate_dataset_registry  # noqa: E402

DEFAULT_CHECKPOINT = (
    "Integrated_dataset/logs/project_reorganization/checkpoints/"
    "pre_reorganization_20260715"
)
MIGRATION_PLANS = (
    "maintenance/reorganization/script_migration_plan.csv",
    "maintenance/reorganization/preintegration_script_migration_plan.csv",
    "maintenance/reorganization/config_data_migration_plan.csv",
    "maintenance/reorganization/config_and_r_migration_plan.csv",
    "maintenance/reorganization/legacy_archive_plan.csv",
    "maintenance/reorganization/retired_workflow_archive_plan.csv",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", default=DEFAULT_CHECKPOINT)
    parser.add_argument("--sample-h5ad-hashes", action="store_true")
    parser.add_argument("--skip-model-load", action="store_true")
    return parser.parse_args()


def check(condition: bool, message: str, failures: list[str]) -> None:
    print(f"{'PASS' if condition else 'FAIL'}  {message}")
    if not condition:
        failures.append(message)


def warn(message: str) -> None:
    print(f"WARN  {message}")


def verify_registries(failures: list[str]) -> None:
    registry_dir = PROJECT_ROOT / "configs" / "datasets"
    result = validate_dataset_registry(registry_dir / "datasets.csv", PROJECT_ROOT)
    check(result.ok, "active dataset paths and dataset IDs validate", failures)
    datasets = read_csv(registry_dir / "datasets.csv")
    libraries = read_csv(registry_dir / "libraries.csv")
    dataset_ids = {row["dataset_id"] for row in datasets}
    library_ids = {row["dataset_id"] for row in libraries}
    placeholders = [
        row for row in libraries
        if row["active"].lower() == "true"
        and row["library_id"] == "__dataset_level_pending_curation__"
    ]
    check(dataset_ids <= library_ids, "every dataset has library registry rows", failures)
    check(not placeholders, "active datasets have no placeholder libraries", failures)


def verify_migrations(failures: list[str]) -> None:
    all_rows: dict[str, list[dict[str, str]]] = {}
    next_path: dict[str, str] = {}
    for relative in MIGRATION_PLANS:
        plan = PROJECT_ROOT / relative
        if not plan.exists():
            check(False, f"migration plan exists: {relative}", failures)
            continue
        rows = read_csv(plan)
        all_rows[relative] = rows
        next_path.update({row["old_path"]: row["new_path"] for row in rows})

    def terminal_path(value: str) -> str:
        seen: set[str] = set()
        while value in next_path:
            if value in seen:
                raise RuntimeError(f"Migration cycle detected at {value}")
            seen.add(value)
            value = next_path[value]
        return value

    for relative, rows in all_rows.items():
        unresolved = []
        for row in rows:
            old = PROJECT_ROOT / row["old_path"]
            terminal = PROJECT_ROOT / terminal_path(row["new_path"])
            if old.exists() or not terminal.exists():
                unresolved.append((row["old_path"], str(terminal)))
        check(not unresolved, f"migration applied: {relative} ({len(rows)} entries)", failures)


def verify_links(failures: list[str]) -> None:
    manifest_path = PROJECT_ROOT / "data" / "registry" / "compatibility_links.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    broken: list[str] = []
    for row in manifest:
        link = PROJECT_ROOT / row["link_path"]
        target = Path(row["target_path"])
        if not target.is_absolute():
            target = PROJECT_ROOT / target
        if not link.is_symlink() or not target.exists() or link.resolve() != target.resolve():
            broken.append(row["link_path"])
    check(not broken, f"compatibility view resolves ({len(manifest)} links)", failures)


def verify_h5ads(checkpoint: Path, sample_hashes: bool, failures: list[str]) -> None:
    baseline_path = checkpoint / "h5ad_structures.json"
    baselines = json.loads(baseline_path.read_text(encoding="utf-8"))
    for row in baselines:
        path = Path(row["path"])
        same_stat = (
            path.exists()
            and path.stat().st_size == row["size_bytes"]
            and path.stat().st_mtime_ns == row["mtime_ns"]
        )
        check(same_stat, f"H5AD stat unchanged: {path.name}", failures)
        if sample_hashes and path.exists():
            check(
                sampled_sha256(path) == row["sampled_sha256"],
                f"H5AD sampled hash unchanged: {path.name}",
                failures,
            )


def verify_models(load_models: bool, failures: list[str]) -> None:
    registry_path = PROJECT_ROOT / "configs" / "models" / "gdtai" / "model_registry.csv"
    for row in read_csv(registry_path):
        artifact = PROJECT_ROOT / row["artifact_path"]
        required = row["required_in_git_checkout"].lower() == "true"
        if not artifact.exists():
            if required:
                check(False, f"required model exists: {row['model_id']}", failures)
            else:
                warn(f"optional experimental model absent: {row['model_id']}")
            continue
        observed = sha256_file(artifact)
        canonical = observed == row["sha256"]
        allowed_override = observed == row["allowed_workspace_sha256"]
        valid = canonical or allowed_override
        if allowed_override and not canonical:
            warn(f"known model workspace override: {row['model_id']} ({observed})")
        else:
            check(valid, f"model checksum: {row['model_id']}", failures)
        if load_models and valid:
            loaded = load_trusted_pickle(artifact, PROJECT_ROOT)
            check(isinstance(loaded, dict), f"model loads: {row['model_id']}", failures)


def verify_compile(failures: list[str]) -> None:
    result = subprocess.run(
        [sys.executable, "-m", "compileall", "-q", "src", "workflows", "tests"],
        cwd=PROJECT_ROOT,
        check=False,
    )
    check(result.returncode == 0, "Python sources compile", failures)


def main() -> int:
    args = parse_args()
    ProjectPaths.discover(PROJECT_ROOT)
    failures: list[str] = []
    checkpoint = Path(args.checkpoint)
    if not checkpoint.is_absolute():
        checkpoint = PROJECT_ROOT / checkpoint
    verify_registries(failures)
    verify_migrations(failures)
    verify_links(failures)
    verify_h5ads(checkpoint, args.sample_h5ad_hashes, failures)
    verify_models(not args.skip_model_load, failures)
    verify_compile(failures)
    print(f"\nverification_failures={len(failures)}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
