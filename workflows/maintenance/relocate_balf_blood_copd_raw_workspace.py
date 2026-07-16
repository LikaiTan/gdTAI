#!/usr/bin/env python3
"""Keep only the BALF_BLOOD_COPD validation H5AD inside the project."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from tnk_atlas.paths import ProjectPaths  # noqa: E402
from tnk_atlas.provenance import atomic_write_json, sampled_sha256, sha256_file  # noqa: E402
from tnk_atlas.registry import DATASET_COLUMNS, FILE_COLUMNS, read_csv  # noqa: E402
from tnk_atlas.storage import (  # noqa: E402
    apply_h5ad_only_external_layout,
    create_compatibility_link,
    rollback_h5ad_only_external_layout,
    validate_h5ad_only_external_layout,
)

DATASET_ID = "BALF_BLOOD_COPD"
DEFAULT_MIGRATION_ID = "balf_blood_copd_h5ad_only_20260716"
REGISTRY_SNAPSHOT = "pre_balf_blood_copd_h5ad_only_20260716"
BASE_MIGRATION_ID = "balf_blood_copd_20260716"
EXTERNAL_WORKSPACE = Path("/home/tanlikai/databank/owndata/singlecell")
SELECTED_RELATIVE = Path("data/results/phase4_final_annotated.h5ad")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--migration-id", default=DEFAULT_MIGRATION_ID)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("plan")
    commands.add_parser("preflight")
    apply_parser = commands.add_parser("apply")
    apply_parser.add_argument("--confirm", action="store_true")
    finalize_parser = commands.add_parser("finalize")
    finalize_parser.add_argument("--confirm", action="store_true")
    commands.add_parser("validate")
    rollback_parser = commands.add_parser("rollback")
    rollback_parser.add_argument("--confirm", action="store_true")
    return parser.parse_args()


def paths() -> ProjectPaths:
    return ProjectPaths.discover(PROJECT_ROOT)


def migration_dir(migration_id: str) -> Path:
    return PROJECT_ROOT / "data" / "registry" / "migrations" / migration_id


def project_workspace() -> Path:
    return paths().dataset(DATASET_ID) / "workspace"


def canonical_artifact() -> Path:
    return (
        paths().dataset_processed(DATASET_ID)
        / "artifacts"
        / "phase4_final_annotated.h5ad"
    )


def project_relative(path: Path) -> str:
    return path.absolute().relative_to(PROJECT_ROOT.absolute()).as_posix()


def h5ad_axis_length(handle: Any, axis: str) -> int:
    frame = handle[axis]
    index_name = frame.attrs.get("_index", "_index")
    if isinstance(index_name, bytes):
        index_name = index_name.decode("utf-8")
    return int(frame[index_name].shape[0])


def file_record(path: Path, relative: Path) -> dict[str, Any]:
    stat = path.stat()
    payload: dict[str, Any] = {
        "relative_path": relative.as_posix(),
        "device": stat.st_dev,
        "inode": stat.st_ino,
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sampled_sha256": sampled_sha256(path),
    }
    if path.suffix == ".h5ad":
        import h5py

        with h5py.File(path, "r") as handle:
            payload.update(
                top_level_keys=sorted(handle.keys()),
                layer_keys=sorted(handle["layers"].keys()),
                n_obs=h5ad_axis_length(handle, "obs"),
                n_vars=h5ad_axis_length(handle, "var"),
            )
    return payload


def validate_record(base: Path, expected: dict[str, Any]) -> str | None:
    path = base / expected["relative_path"]
    if not path.is_file():
        return f"sentinel is missing: {path}"
    observed = file_record(path, Path(expected["relative_path"]))
    for field in (
        "device",
        "inode",
        "size_bytes",
        "mtime_ns",
        "sampled_sha256",
        "top_level_keys",
        "layer_keys",
        "n_obs",
        "n_vars",
    ):
        if field in expected and observed.get(field) != expected[field]:
            return f"sentinel {field} mismatch: {path}"
    return None


def write_csv_atomic(
    path: Path,
    fieldnames: tuple[str, ...],
    rows: list[dict[str, Any]],
) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def replace_generated_link(link: Path, target: Path) -> None:
    if link.is_symlink():
        link.unlink()
    elif link.exists():
        raise FileExistsError(f"Refusing to replace non-symlink path: {link}")
    link.parent.mkdir(parents=True, exist_ok=True)
    create_compatibility_link(link, target)


def load_base_plan() -> dict[str, Any]:
    plan_path = (
        PROJECT_ROOT
        / "data"
        / "registry"
        / "migrations"
        / BASE_MIGRATION_ID
        / "plan.json"
    )
    return json.loads(plan_path.read_text(encoding="utf-8"))


def load_plan(state_dir: Path) -> dict[str, Any]:
    plan_path = state_dir / "plan.json"
    if not plan_path.exists():
        raise FileNotFoundError(f"Migration plan is missing: {plan_path}")
    return json.loads(plan_path.read_text(encoding="utf-8"))


def command_plan(state_dir: Path) -> int:
    state_dir.mkdir(parents=True, exist_ok=True)
    plan_path = state_dir / "plan.json"
    if plan_path.exists():
        raise FileExistsError(f"Migration plan already exists: {plan_path}")

    workspace = project_workspace()
    external = EXTERNAL_WORKSPACE
    artifact = canonical_artifact()
    selected = workspace / SELECTED_RELATIVE
    raw_source = paths().dataset_raw(DATASET_ID) / "legacy_source"
    current = paths().dataset_current_h5ad(DATASET_ID)
    if not workspace.is_dir() or workspace.is_symlink():
        raise FileNotFoundError(f"Project workspace is not physical: {workspace}")
    if not external.is_symlink() or external.resolve() != workspace.resolve():
        raise RuntimeError(f"External workspace is not the expected compatibility link: {external}")
    if not selected.is_file() or selected.is_symlink():
        raise FileNotFoundError(f"Selected H5AD is not physical: {selected}")
    for link, target in (
        (artifact, selected),
        (current, selected),
        (raw_source, workspace),
    ):
        if not link.is_symlink() or link.resolve() != target.resolve():
            raise RuntimeError(f"Current generated link is invalid: {link}")

    base_plan = load_base_plan()
    raw_sentinels = [
        record
        for record in base_plan["sentinels"]
        if record["relative_path"] != SELECTED_RELATIVE.as_posix()
    ]
    for record in raw_sentinels:
        error = validate_record(workspace, record)
        if error:
            raise RuntimeError(error)

    workspace_stat = workspace.stat()
    h5ad = file_record(selected, SELECTED_RELATIVE)
    snapshot = (
        PROJECT_ROOT
        / "data"
        / "registry"
        / "snapshots"
        / REGISTRY_SNAPSHOT
        / "snapshot.json"
    )
    if not snapshot.exists():
        raise FileNotFoundError(f"Registry snapshot is missing: {snapshot}")
    destination_device = external.parent.stat().st_dev
    plan = {
        "migration_id": state_dir.name,
        "dataset_id": DATASET_ID,
        "created_at": datetime.now().astimezone().isoformat(),
        "git_head": subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=PROJECT_ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip(),
        "source_layout": "full_workspace_in_project",
        "target_layout": "selected_h5ad_only_in_project",
        "external_workspace": str(external),
        "project_workspace": str(workspace),
        "selected_relative": SELECTED_RELATIVE.as_posix(),
        "canonical_artifact": str(artifact),
        "current_h5ad": str(current),
        "raw_source_link": str(raw_source),
        "workspace_device": workspace_stat.st_dev,
        "workspace_inode": workspace_stat.st_ino,
        "destination_device": destination_device,
        "same_filesystem": destination_device == workspace_stat.st_dev,
        "top_level_entries": sorted(path.name for path in workspace.iterdir()),
        "h5ad": h5ad,
        "raw_sentinels": raw_sentinels,
        "registry_snapshot": REGISTRY_SNAPSHOT,
        "base_migration_id": BASE_MIGRATION_ID,
    }
    atomic_write_json(plan_path, plan)
    shutil.copy2(paths().dataset(DATASET_ID) / "dataset.json", state_dir / "dataset_before.json")
    print(
        f"plan={plan_path} raw_sentinels={len(raw_sentinels)} "
        f"same_filesystem={plan['same_filesystem']}"
    )
    return 0


def preflight_failures(plan: dict[str, Any]) -> list[str]:
    workspace = Path(plan["project_workspace"])
    external = Path(plan["external_workspace"])
    artifact = Path(plan["canonical_artifact"])
    current = Path(plan["current_h5ad"])
    raw_source = Path(plan["raw_source_link"])
    selected = workspace / plan["selected_relative"]
    failures: list[str] = []

    if not workspace.is_dir() or workspace.is_symlink():
        failures.append(f"project workspace is not physical: {workspace}")
    else:
        stat = workspace.stat()
        if stat.st_dev != plan["workspace_device"] or stat.st_ino != plan["workspace_inode"]:
            failures.append("project workspace inode changed since planning")
        if sorted(path.name for path in workspace.iterdir()) != plan["top_level_entries"]:
            failures.append("workspace top-level entries changed since planning")
    if not external.is_symlink() or external.resolve() != workspace.resolve():
        failures.append(f"external workspace link mismatch: {external}")
    if not artifact.is_symlink() or artifact.resolve() != selected.resolve():
        failures.append(f"artifact link mismatch: {artifact}")
    if not current.is_symlink() or current.resolve() != selected.resolve():
        failures.append(f"current-H5AD link mismatch: {current}")
    if not raw_source.is_symlink() or raw_source.resolve() != workspace.resolve():
        failures.append(f"raw-source link mismatch: {raw_source}")
    h5ad_error = validate_record(workspace, plan["h5ad"])
    if h5ad_error:
        failures.append(h5ad_error)
    for record in plan["raw_sentinels"]:
        error = validate_record(workspace, record)
        if error:
            failures.append(error)
    if not plan["same_filesystem"]:
        failures.append("project and original workspace locations are not on one filesystem")
    snapshot = (
        PROJECT_ROOT
        / "data"
        / "registry"
        / "snapshots"
        / plan["registry_snapshot"]
        / "snapshot.json"
    )
    if not snapshot.exists():
        failures.append(f"registry rollback snapshot is missing: {snapshot}")
    return failures


def command_preflight(state_dir: Path) -> int:
    plan = load_plan(state_dir)
    failures = preflight_failures(plan)
    result = {
        "migration_id": state_dir.name,
        "dataset_id": DATASET_ID,
        "failure_count": len(failures),
        "failures": failures,
        "validated_at": datetime.now().astimezone().isoformat(),
    }
    atomic_write_json(state_dir / "preflight.json", result)
    print(json.dumps(result, indent=2))
    return 1 if failures else 0


def layout_errors(plan: dict[str, Any]) -> list[str]:
    return validate_h5ad_only_external_layout(
        Path(plan["external_workspace"]),
        Path(plan["project_workspace"]),
        Path(plan["selected_relative"]),
        Path(plan["canonical_artifact"]),
        expected_workspace_device=plan["workspace_device"],
        expected_workspace_inode=plan["workspace_inode"],
        expected_file_device=plan["h5ad"]["device"],
        expected_file_inode=plan["h5ad"]["inode"],
    )


def command_apply(state_dir: Path, confirmed: bool) -> int:
    if not confirmed:
        raise ValueError("Refusing physical relocation without --confirm")
    plan = load_plan(state_dir)
    external = Path(plan["external_workspace"])
    workspace = Path(plan["project_workspace"])
    artifact = Path(plan["canonical_artifact"])
    if external.is_symlink():
        failures = preflight_failures(plan)
        if failures:
            raise RuntimeError("Preflight failed:\n" + "\n".join(failures))

    status = apply_h5ad_only_external_layout(
        external,
        workspace,
        Path(plan["selected_relative"]),
        artifact,
        expected_workspace_device=plan["workspace_device"],
        expected_workspace_inode=plan["workspace_inode"],
        expected_file_device=plan["h5ad"]["device"],
        expected_file_inode=plan["h5ad"]["inode"],
    )
    replace_generated_link(Path(plan["raw_source_link"]), external)
    replace_generated_link(Path(plan["current_h5ad"]), artifact)

    failures = layout_errors(plan)
    for record in plan["raw_sentinels"]:
        error = validate_record(external, record)
        if error:
            failures.append(error)
    if failures:
        raise RuntimeError("Applied layout failed validation:\n" + "\n".join(failures))
    atomic_write_json(
        state_dir / "apply_summary.json",
        {
            "migration_id": state_dir.name,
            "dataset_id": DATASET_ID,
            "status": status,
            "external_workspace": str(external),
            "canonical_h5ad": project_relative(artifact),
            "completed_at": datetime.now().astimezone().isoformat(),
        },
    )
    print(f"{status}: raw workspace -> {external}; H5AD -> {artifact}")
    return 0


def upsert_registries(plan: dict[str, Any]) -> None:
    dataset_path = PROJECT_ROOT / "configs" / "datasets" / "datasets.csv"
    datasets = read_csv(dataset_path)
    matched = False
    for row in datasets:
        if row["dataset_id"] != DATASET_ID:
            continue
        matched = True
        row["legacy_raw_path"] = plan["external_workspace"]
        row["processed_h5ad_path"] = project_relative(Path(plan["current_h5ad"]))
        row["notes"] = (
            "Independent gdTAI validation cohort: four 10x 5' BALF/PBMC "
            "libraries from six donors. Only the selected validation H5AD is "
            "stored physically in this project; raw and interim files remain "
            "under the original local workspace. Not approved for atlas integration."
        )
    if not matched:
        raise KeyError(DATASET_ID)
    write_csv_atomic(dataset_path, DATASET_COLUMNS, datasets)

    files_path = PROJECT_ROOT / "configs" / "datasets" / "files.csv"
    files = [row for row in read_csv(files_path) if row["dataset_id"] != DATASET_ID]
    external = Path(plan["external_workspace"])
    raw_source = Path(plan["raw_source_link"])
    for record in [plan["h5ad"], *plan["raw_sentinels"]]:
        relative = Path(record["relative_path"])
        if relative == SELECTED_RELATIVE:
            canonical = Path(plan["canonical_artifact"])
            lifecycle = "processed"
            role = "independent_external_validation_h5ad"
        else:
            canonical = raw_source / relative
            lifecycle = record["lifecycle_level"]
            role = record["file_role"]
        files.append(
            {
                "dataset_id": DATASET_ID,
                "lifecycle_level": lifecycle,
                "file_role": role,
                "legacy_path": str(external / relative),
                "canonical_path": project_relative(canonical),
                "size_bytes": str(record["size_bytes"]),
                "mtime_ns": str(record["mtime_ns"]),
                "checksum_type": "sampled_sha256",
                "checksum": record["sampled_sha256"],
            }
        )
    files.sort(
        key=lambda row: (
            row["dataset_id"],
            row["lifecycle_level"],
            row["canonical_path"],
        )
    )
    write_csv_atomic(files_path, FILE_COLUMNS, files)


def write_storage_index() -> None:
    project_paths = paths()
    rows = read_csv(PROJECT_ROOT / "configs" / "datasets" / "datasets.csv")
    index_rows = []
    for row in rows:
        dataset_id = row["dataset_id"]
        current = project_paths.dataset_current_h5ad(dataset_id)
        index_rows.append(
            {
                "dataset_id": dataset_id,
                "canonical_dataset_path": project_relative(project_paths.dataset(dataset_id)),
                "raw_exists": str(project_paths.dataset_raw(dataset_id).exists()).lower(),
                "interim_exists": str(project_paths.dataset_interim(dataset_id).exists()).lower(),
                "current_h5ad": project_relative(current) if current.exists() else "",
                "legacy_raw_path": row["legacy_raw_path"],
            }
        )
    write_csv_atomic(
        PROJECT_ROOT / "data" / "registry" / "storage_index.csv",
        (
            "dataset_id",
            "canonical_dataset_path",
            "raw_exists",
            "interim_exists",
            "current_h5ad",
            "legacy_raw_path",
        ),
        index_rows,
    )


def write_dataset_json(state_dir: Path, plan: dict[str, Any]) -> None:
    project_paths = paths()
    dataset_root = project_paths.dataset(DATASET_ID)
    atomic_write_json(
        dataset_root / "dataset.json",
        {
            "dataset_id": DATASET_ID,
            "source_type": "local_external",
            "canonical_dataset_path": project_relative(dataset_root),
            "physical_raw_workspace": plan["external_workspace"],
            "raw_view": project_relative(project_paths.dataset_raw(DATASET_ID)),
            "processed_path": project_relative(project_paths.dataset_processed(DATASET_ID)),
            "current_h5ad": project_relative(Path(plan["current_h5ad"])),
            "current_h5ad_target": project_relative(Path(plan["canonical_artifact"])),
            "legacy_h5ad_path": str(
                Path(plan["external_workspace"]) / plan["selected_relative"]
            ),
            "storage_policy": "selected_h5ad_only_in_project",
            "storage_migration_id": state_dir.name,
            "integration_role": "gdTAI_independent_external_validation",
        },
    )


def command_finalize(state_dir: Path, confirmed: bool) -> int:
    if not confirmed:
        raise ValueError("Refusing registry finalization without --confirm")
    plan = load_plan(state_dir)
    failures = layout_errors(plan)
    if failures:
        raise RuntimeError("Physical layout is invalid:\n" + "\n".join(failures))
    upsert_registries(plan)
    write_dataset_json(state_dir, plan)
    write_storage_index()
    subprocess.run(
        [
            sys.executable,
            "workflows/maintenance/build_data_compatibility_view.py",
            "--apply",
            "--replace-generated-links",
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )
    atomic_write_json(
        state_dir / "finalize_summary.json",
        {
            "migration_id": state_dir.name,
            "dataset_id": DATASET_ID,
            "status": "h5ad_only_registry_cutover_complete",
            "physical_project_data": project_relative(Path(plan["canonical_artifact"])),
            "physical_raw_workspace": plan["external_workspace"],
            "completed_at": datetime.now().astimezone().isoformat(),
        },
    )
    return 0


def validation_failures(plan: dict[str, Any]) -> list[str]:
    project_paths = paths()
    external = Path(plan["external_workspace"])
    artifact = Path(plan["canonical_artifact"])
    current = Path(plan["current_h5ad"])
    raw_source = Path(plan["raw_source_link"])
    failures = layout_errors(plan)

    h5ad_error = validate_record(artifact.parent, {
        **plan["h5ad"],
        "relative_path": artifact.name,
    })
    if h5ad_error:
        failures.append(h5ad_error)
    if not raw_source.is_symlink() or raw_source.resolve() != external.resolve():
        failures.append(f"raw-source link mismatch: {raw_source}")
    if not current.is_symlink() or current.resolve() != artifact.resolve():
        failures.append(f"current-H5AD link mismatch: {current}")

    for record in plan["raw_sentinels"]:
        error = validate_record(external, record)
        if error:
            failures.append(error)
        legacy = external / record["relative_path"]
        canonical = raw_source / record["relative_path"]
        if not canonical.exists():
            failures.append(f"canonical raw view is missing: {canonical}")
        elif (legacy.stat().st_dev, legacy.stat().st_ino) != (
            canonical.stat().st_dev,
            canonical.stat().st_ino,
        ):
            failures.append(f"canonical raw view inode mismatch: {canonical}")

    dataset_root = project_paths.dataset(DATASET_ID)
    if (dataset_root / "workspace").exists() or (dataset_root / "workspace").is_symlink():
        failures.append("project-hosted workspace still exists")
    raw_children = list(project_paths.dataset_raw(DATASET_ID).iterdir())
    if raw_children != [raw_source]:
        failures.append("project raw directory contains unexpected entries")
    artifact_children = list(artifact.parent.iterdir())
    if artifact_children != [artifact]:
        failures.append("project artifact directory contains unexpected entries")

    dataset_json = json.loads((dataset_root / "dataset.json").read_text(encoding="utf-8"))
    if dataset_json.get("storage_policy") != "selected_h5ad_only_in_project":
        failures.append("dataset.json storage policy mismatch")

    dataset_rows = [
        row
        for row in read_csv(PROJECT_ROOT / "configs" / "datasets" / "datasets.csv")
        if row["dataset_id"] == DATASET_ID
    ]
    if len(dataset_rows) != 1:
        failures.append(f"dataset registry row count is {len(dataset_rows)}, expected 1")
    elif dataset_rows[0]["processed_h5ad_path"] != project_relative(current):
        failures.append("dataset registry selected-H5AD path mismatch")

    file_rows = [
        row
        for row in read_csv(PROJECT_ROOT / "configs" / "datasets" / "files.csv")
        if row["dataset_id"] == DATASET_ID
    ]
    if len(file_rows) != 1 + len(plan["raw_sentinels"]):
        failures.append("file registry row count mismatch")
    for row in file_rows:
        legacy = Path(row["legacy_path"])
        canonical = PROJECT_ROOT / row["canonical_path"]
        if not legacy.exists() or not canonical.exists():
            failures.append(f"registered path is missing: {row['canonical_path']}")
        elif (legacy.stat().st_dev, legacy.stat().st_ino) != (
            canonical.stat().st_dev,
            canonical.stat().st_ino,
        ):
            failures.append(f"legacy/canonical inode mismatch: {row['canonical_path']}")
    return failures


def command_validate(state_dir: Path) -> int:
    plan = load_plan(state_dir)
    failures = validation_failures(plan)
    result = {
        "migration_id": state_dir.name,
        "dataset_id": DATASET_ID,
        "failure_count": len(failures),
        "failures": failures,
        "validated_at": datetime.now().astimezone().isoformat(),
    }
    atomic_write_json(state_dir / "validation.json", result)
    print(json.dumps(result, indent=2))
    return 1 if failures else 0


def restore_registry_snapshot(snapshot_id: str) -> None:
    snapshot = PROJECT_ROOT / "data" / "registry" / "snapshots" / snapshot_id
    manifest = json.loads((snapshot / "snapshot.json").read_text(encoding="utf-8"))
    for name, metadata in manifest["files"].items():
        source = snapshot / name
        if sha256_file(source) != metadata["sha256"]:
            raise RuntimeError(f"Registry snapshot checksum mismatch: {source}")
        destination = PROJECT_ROOT / "configs" / "datasets" / name
        temporary = destination.with_suffix(destination.suffix + ".rollback.tmp")
        shutil.copy2(source, temporary)
        os.replace(temporary, destination)


def command_rollback(state_dir: Path, confirmed: bool) -> int:
    if not confirmed:
        raise ValueError("Refusing physical rollback without --confirm")
    plan = load_plan(state_dir)
    status = rollback_h5ad_only_external_layout(
        Path(plan["external_workspace"]),
        Path(plan["project_workspace"]),
        Path(plan["selected_relative"]),
        Path(plan["canonical_artifact"]),
        expected_workspace_device=plan["workspace_device"],
        expected_workspace_inode=plan["workspace_inode"],
        expected_file_device=plan["h5ad"]["device"],
        expected_file_inode=plan["h5ad"]["inode"],
    )
    workspace = Path(plan["project_workspace"])
    selected = workspace / plan["selected_relative"]
    replace_generated_link(Path(plan["raw_source_link"]), workspace)
    replace_generated_link(Path(plan["current_h5ad"]), selected)
    restore_registry_snapshot(plan["registry_snapshot"])
    shutil.copy2(state_dir / "dataset_before.json", paths().dataset(DATASET_ID) / "dataset.json")
    write_storage_index()
    subprocess.run(
        [
            sys.executable,
            "workflows/maintenance/build_data_compatibility_view.py",
            "--apply",
            "--replace-generated-links",
        ],
        cwd=PROJECT_ROOT,
        check=True,
    )
    atomic_write_json(
        state_dir / "rollback_summary.json",
        {
            "migration_id": state_dir.name,
            "dataset_id": DATASET_ID,
            "status": status,
            "completed_at": datetime.now().astimezone().isoformat(),
        },
    )
    print(f"{status}: raw workspace -> {workspace}")
    return 0


def main() -> int:
    args = parse_args()
    state_dir = migration_dir(args.migration_id)
    if args.command == "plan":
        return command_plan(state_dir)
    if args.command == "preflight":
        return command_preflight(state_dir)
    if args.command == "apply":
        return command_apply(state_dir, args.confirm)
    if args.command == "finalize":
        return command_finalize(state_dir, args.confirm)
    if args.command == "validate":
        return command_validate(state_dir)
    if args.command == "rollback":
        return command_rollback(state_dir, args.confirm)
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
