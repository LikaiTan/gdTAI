#!/usr/bin/env python3
"""Move the BALF/PBMC COPD validation study into canonical dataset storage."""

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
from tnk_atlas.registry import (  # noqa: E402
    DATASET_COLUMNS,
    FILE_COLUMNS,
    LIBRARY_COLUMNS,
    read_csv,
)
from tnk_atlas.storage import (  # noqa: E402
    apply_absolute_move_with_link,
    create_compatibility_link,
    rollback_absolute_move,
    validate_absolute_move,
)

DATASET_ID = "BALF_BLOOD_COPD"
DEFAULT_MIGRATION_ID = "balf_blood_copd_20260716"
REGISTRY_SNAPSHOT = "pre_balf_blood_copd_intake_20260716"
SOURCE_WORKSPACE = Path("/home/tanlikai/databank/owndata/singlecell")
SELECTED_H5AD_RELATIVE = Path("data/results/phase4_final_annotated.h5ad")

SENTINELS = (
    ("processed", "independent_external_validation_h5ad", SELECTED_H5AD_RELATIVE),
    (
        "raw",
        "filtered_count_matrix",
        Path("data/matrix/25060703_LIB1_5/filtered_feature_bc_matrix/matrix.mtx.gz"),
    ),
    (
        "raw",
        "filtered_count_matrix",
        Path("data/matrix/25060703_LIB2_5/filtered_feature_bc_matrix/matrix.mtx.gz"),
    ),
    (
        "raw",
        "filtered_count_matrix",
        Path("data/matrix/25060703_LIB3_5/filtered_feature_bc_matrix/matrix.mtx.gz"),
    ),
    (
        "raw",
        "filtered_count_matrix",
        Path("data/matrix/25060703_LIB4_5/filtered_feature_bc_matrix/matrix.mtx.gz"),
    ),
    (
        "raw",
        "raw_count_matrix",
        Path("data/raw_matrix/25060703_LIB1_5/raw_feature_bc_matrix/matrix.mtx.gz"),
    ),
    (
        "raw",
        "raw_count_matrix",
        Path("data/raw_matrix/25060703_LIB2_5/raw_feature_bc_matrix/matrix.mtx.gz"),
    ),
    (
        "raw",
        "raw_count_matrix",
        Path("data/raw_matrix/25060703_LIB3_5/raw_feature_bc_matrix/matrix.mtx.gz"),
    ),
    (
        "raw",
        "raw_count_matrix",
        Path("data/raw_matrix/25060703_LIB4_5/raw_feature_bc_matrix/matrix.mtx.gz"),
    ),
    (
        "interim",
        "harmonized_tcr_table",
        Path("data/TCR/tcr_integration_output_allchains/tcr_final_wide.tsv"),
    ),
    (
        "interim",
        "demultiplexing_assignment",
        Path("data/demultiplexing_Result/donor_ids_1_3.tsv"),
    ),
    (
        "interim",
        "demultiplexing_assignment",
        Path("data/demultiplexing_Result/donor_ids_2_4.tsv"),
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--migration-id", default=DEFAULT_MIGRATION_ID)
    parser.add_argument("--source-workspace", type=Path, default=SOURCE_WORKSPACE)
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


def migration_dir(migration_id: str) -> Path:
    return PROJECT_ROOT / "data" / "registry" / "migrations" / migration_id


def canonical_workspace() -> Path:
    return ProjectPaths.discover(PROJECT_ROOT).dataset(DATASET_ID) / "workspace"


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


def command_plan(state_dir: Path, source: Path) -> int:
    state_dir.mkdir(parents=True, exist_ok=True)
    plan_path = state_dir / "plan.json"
    if plan_path.exists():
        raise FileExistsError(f"Migration plan already exists: {plan_path}")
    if source.is_symlink() or not source.is_dir():
        raise FileNotFoundError(f"Source workspace is not a physical directory: {source}")
    destination = canonical_workspace()
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"Canonical workspace already exists: {destination}")
    stat = source.stat()
    destination_parent = destination.parent
    while not destination_parent.exists():
        destination_parent = destination_parent.parent
    records = []
    for lifecycle, role, relative in SENTINELS:
        path = source / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        record = file_record(path, relative)
        record.update(lifecycle_level=lifecycle, file_role=role)
        records.append(record)
    top_level_entries = sorted(path.name for path in source.iterdir())
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
        "source_workspace": str(source),
        "canonical_workspace": str(destination),
        "selected_h5ad_relative": SELECTED_H5AD_RELATIVE.as_posix(),
        "workspace_device": stat.st_dev,
        "workspace_inode": stat.st_ino,
        "workspace_mtime_ns": stat.st_mtime_ns,
        "destination_device": destination_parent.stat().st_dev,
        "same_filesystem": destination_parent.stat().st_dev == stat.st_dev,
        "top_level_entries": top_level_entries,
        "sentinels": records,
        "registry_snapshot": REGISTRY_SNAPSHOT,
    }
    atomic_write_json(plan_path, plan)
    print(
        f"plan={plan_path} sentinels={len(records)} "
        f"same_filesystem={plan['same_filesystem']}"
    )
    return 0


def load_plan(state_dir: Path) -> dict[str, Any]:
    path = state_dir / "plan.json"
    if not path.exists():
        raise FileNotFoundError(f"Migration plan is missing: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


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


def preflight_failures(state_dir: Path, plan: dict[str, Any]) -> list[str]:
    source = Path(plan["source_workspace"])
    destination = Path(plan["canonical_workspace"])
    failures: list[str] = []
    if source.is_symlink():
        failures.append(f"source workspace is already a symlink: {source}")
    elif not source.is_dir():
        failures.append(f"source workspace is missing: {source}")
    else:
        stat = source.stat()
        if (
            stat.st_dev != plan["workspace_device"]
            or stat.st_ino != plan["workspace_inode"]
            or stat.st_mtime_ns != plan["workspace_mtime_ns"]
        ):
            failures.append("workspace stat/inode changed since planning")
        if sorted(path.name for path in source.iterdir()) != plan["top_level_entries"]:
            failures.append("workspace top-level entries changed since planning")
        for expected in plan["sentinels"]:
            error = validate_record(source, expected)
            if error:
                failures.append(error)
    if destination.exists() or destination.is_symlink():
        failures.append(f"canonical workspace already exists: {destination}")
    if not plan["same_filesystem"]:
        failures.append("source and destination are not on the same filesystem")
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
    failures = preflight_failures(state_dir, plan)
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


def command_apply(state_dir: Path, confirmed: bool) -> int:
    if not confirmed:
        raise ValueError("Refusing physical workspace migration without --confirm")
    plan = load_plan(state_dir)
    source = Path(plan["source_workspace"])
    destination = Path(plan["canonical_workspace"])
    if not source.is_symlink():
        failures = preflight_failures(state_dir, plan)
        if failures:
            raise RuntimeError("Preflight failed:\n" + "\n".join(failures))
    status = apply_absolute_move_with_link(
        source,
        destination,
        expected_device=plan["workspace_device"],
        expected_inode=plan["workspace_inode"],
    )
    atomic_write_json(
        state_dir / "apply_summary.json",
        {
            "migration_id": state_dir.name,
            "dataset_id": DATASET_ID,
            "status": status,
            "completed_at": datetime.now().astimezone().isoformat(),
        },
    )
    print(f"{status}: {source} -> {destination}")
    return 0


def replace_generated_link(link: Path, target: Path) -> None:
    if link.is_symlink():
        link.unlink()
    elif link.exists():
        raise FileExistsError(f"Refusing to replace non-symlink generated path: {link}")
    link.parent.mkdir(parents=True, exist_ok=True)
    create_compatibility_link(link, target)


def upsert_registries(plan: dict[str, Any]) -> None:
    paths = ProjectPaths.discover(PROJECT_ROOT)
    current = paths.dataset_current_h5ad(DATASET_ID)
    datasets_path = PROJECT_ROOT / "configs" / "datasets" / "datasets.csv"
    datasets = [
        row for row in read_csv(datasets_path) if row["dataset_id"] != DATASET_ID
    ]
    datasets.append(
        {
            "dataset_id": DATASET_ID,
            "source_type": "local_external",
            "accession": DATASET_ID,
            "phase0_active": "false",
            "current_milestone_active": "false",
            "current_cell_count": "0",
            "extended_atlas_active": "false",
            "extended_atlas_cell_count": "0",
            "integration_role": "gdTAI_independent_external_validation",
            "legacy_raw_path": plan["source_workspace"],
            "processed_h5ad_path": project_relative(current),
            "status": "external_validation_only",
            "exclusion_reason": "",
            "notes": (
                "Independent gdTAI validation cohort: four 10x 5' libraries, "
                "paired BALF/PBMC samples from six donors with COPD or lung nodules. "
                "Not approved for atlas integration."
            ),
        }
    )
    datasets.sort(key=lambda row: row["dataset_id"])
    write_csv_atomic(datasets_path, DATASET_COLUMNS, datasets)

    libraries_path = PROJECT_ROOT / "configs" / "datasets" / "libraries.csv"
    libraries = [
        row for row in read_csv(libraries_path) if row["dataset_id"] != DATASET_ID
    ]
    tissue_by_library = {
        "LIB1": "BALF",
        "LIB2": "BALF",
        "LIB3": "PBMC",
        "LIB4": "PBMC",
    }
    for library_id, tissue in tissue_by_library.items():
        libraries.append(
            {
                "dataset_id": DATASET_ID,
                "library_id": library_id,
                "sample_id": library_id,
                "rna_assay": "10x 5'",
                "tcr_assay": "vdj",
                "tcr_scope": "ab_and_gd",
                "active": "false",
                "notes": (
                    f"{tissue}; multiplexed donor library. External validation only; "
                    "not active for atlas integration."
                ),
            }
        )
    libraries.sort(key=lambda row: (row["dataset_id"], row["library_id"]))
    write_csv_atomic(libraries_path, LIBRARY_COLUMNS, libraries)

    files_path = PROJECT_ROOT / "configs" / "datasets" / "files.csv"
    files = [row for row in read_csv(files_path) if row["dataset_id"] != DATASET_ID]
    source = Path(plan["source_workspace"])
    destination = Path(plan["canonical_workspace"])
    for record in plan["sentinels"]:
        relative = Path(record["relative_path"])
        files.append(
            {
                "dataset_id": DATASET_ID,
                "lifecycle_level": record["lifecycle_level"],
                "file_role": record["file_role"],
                "legacy_path": str(source / relative),
                "canonical_path": project_relative(destination / relative),
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
    paths = ProjectPaths.discover(PROJECT_ROOT)
    rows = read_csv(PROJECT_ROOT / "configs" / "datasets" / "datasets.csv")
    index_rows = []
    for row in rows:
        dataset_id = row["dataset_id"]
        current = paths.dataset_current_h5ad(dataset_id)
        index_rows.append(
            {
                "dataset_id": dataset_id,
                "canonical_dataset_path": project_relative(paths.dataset(dataset_id)),
                "raw_exists": str(paths.dataset_raw(dataset_id).exists()).lower(),
                "interim_exists": str(paths.dataset_interim(dataset_id).exists()).lower(),
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


def command_finalize(state_dir: Path, confirmed: bool) -> int:
    if not confirmed:
        raise ValueError("Refusing registry finalization without --confirm")
    plan = load_plan(state_dir)
    source = Path(plan["source_workspace"])
    workspace = Path(plan["canonical_workspace"])
    errors = validate_absolute_move(
        source,
        workspace,
        expected_device=plan["workspace_device"],
        expected_inode=plan["workspace_inode"],
    )
    for expected in plan["sentinels"]:
        error = validate_record(workspace, expected)
        if error:
            errors.append(error)
    if errors:
        raise RuntimeError("Physical migration is not valid:\n" + "\n".join(errors))

    paths = ProjectPaths.discover(PROJECT_ROOT)
    dataset_root = paths.dataset(DATASET_ID)
    selected = workspace / SELECTED_H5AD_RELATIVE
    raw_source = paths.dataset_raw(DATASET_ID) / "legacy_source"
    artifact = (
        paths.dataset_processed(DATASET_ID)
        / "artifacts"
        / "phase4_final_annotated.h5ad"
    )
    current = paths.dataset_current_h5ad(DATASET_ID)
    replace_generated_link(raw_source, workspace)
    replace_generated_link(artifact, selected)
    replace_generated_link(current, artifact)
    atomic_write_json(
        dataset_root / "dataset.json",
        {
            "dataset_id": DATASET_ID,
            "source_type": "local_external",
            "canonical_dataset_path": project_relative(dataset_root),
            "workspace_path": project_relative(workspace),
            "raw_path": project_relative(paths.dataset_raw(DATASET_ID)),
            "processed_path": project_relative(paths.dataset_processed(DATASET_ID)),
            "current_h5ad": project_relative(current),
            "current_h5ad_target": project_relative(selected),
            "legacy_workspace_path": str(source),
            "storage_migration_id": state_dir.name,
            "integration_role": "gdTAI_independent_external_validation",
        },
    )
    upsert_registries(plan)
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
            "library_count": 4,
            "sentinel_count": len(plan["sentinels"]),
            "selected_h5ad": project_relative(current),
            "completed_at": datetime.now().astimezone().isoformat(),
            "status": "canonical_registry_cutover_complete",
        },
    )
    return 0


def validation_failures(plan: dict[str, Any]) -> list[str]:
    paths = ProjectPaths.discover(PROJECT_ROOT)
    source = Path(plan["source_workspace"])
    workspace = Path(plan["canonical_workspace"])
    failures = validate_absolute_move(
        source,
        workspace,
        expected_device=plan["workspace_device"],
        expected_inode=plan["workspace_inode"],
    )
    for expected in plan["sentinels"]:
        error = validate_record(workspace, expected)
        if error:
            failures.append(error)

    artifact = (
        paths.dataset_processed(DATASET_ID)
        / "artifacts"
        / "phase4_final_annotated.h5ad"
    )
    current = paths.dataset_current_h5ad(DATASET_ID)
    selected = workspace / SELECTED_H5AD_RELATIVE
    for link, target in (
        (paths.dataset_raw(DATASET_ID) / "legacy_source", workspace),
        (artifact, selected),
        (current, artifact),
        (
            paths.data / "raw" / "local" / DATASET_ID / "legacy_source",
            workspace,
        ),
        (
            paths.data / "processed" / "per_dataset" / DATASET_ID / "current.h5ad",
            current,
        ),
    ):
        if not link.is_symlink() or link.resolve() != target.resolve():
            failures.append(f"generated link mismatch: {link}")

    dataset_rows = [
        row
        for row in read_csv(PROJECT_ROOT / "configs" / "datasets" / "datasets.csv")
        if row["dataset_id"] == DATASET_ID
    ]
    if len(dataset_rows) != 1:
        failures.append(f"dataset registry row count is {len(dataset_rows)}, expected 1")
    elif dataset_rows[0]["processed_h5ad_path"] != project_relative(current):
        failures.append("dataset registry selected H5AD path mismatch")
    library_rows = [
        row
        for row in read_csv(PROJECT_ROOT / "configs" / "datasets" / "libraries.csv")
        if row["dataset_id"] == DATASET_ID
    ]
    if {row["library_id"] for row in library_rows} != {"LIB1", "LIB2", "LIB3", "LIB4"}:
        failures.append("library registry does not contain exactly LIB1-LIB4")
    file_rows = [
        row
        for row in read_csv(PROJECT_ROOT / "configs" / "datasets" / "files.csv")
        if row["dataset_id"] == DATASET_ID
    ]
    if len(file_rows) != len(plan["sentinels"]):
        failures.append(
            f"file registry row count is {len(file_rows)}, "
            f"expected {len(plan['sentinels'])}"
        )
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


def remove_generated_dataset_paths(dataset_root: Path) -> list[str]:
    expected_links = (
        dataset_root / "raw" / "legacy_source",
        dataset_root / "processed" / "artifacts" / "phase4_final_annotated.h5ad",
        dataset_root / "processed" / "current.h5ad",
    )
    for link in expected_links:
        if link.is_symlink():
            link.unlink()
    metadata = dataset_root / "dataset.json"
    if metadata.exists():
        metadata.unlink()
    for directory in (
        dataset_root / "processed" / "artifacts",
        dataset_root / "processed",
        dataset_root / "raw",
        dataset_root / "interim",
    ):
        try:
            directory.rmdir()
        except OSError:
            pass
    leftovers = [
        project_relative(path)
        for path in dataset_root.iterdir()
        if path.name != "workspace"
        if path.exists() or path.is_symlink()
    ]
    return leftovers


def remove_lifecycle_view(path: Path, expected_links: tuple[str, ...]) -> None:
    if not path.exists() and not path.is_symlink():
        return
    if path.is_symlink():
        path.unlink()
        return
    for name in expected_links:
        link = path / name
        if link.is_symlink():
            link.unlink()
        elif link.exists():
            raise RuntimeError(f"Refusing to remove a non-symlink lifecycle path: {link}")
    unexpected = list(path.iterdir())
    if unexpected:
        raise RuntimeError(
            "Refusing to remove a nonempty lifecycle view:\n"
            + "\n".join(str(candidate) for candidate in unexpected)
        )
    path.rmdir()


def command_rollback(state_dir: Path, confirmed: bool) -> int:
    if not confirmed:
        raise ValueError("Refusing workspace rollback without --confirm")
    plan = load_plan(state_dir)
    paths = ProjectPaths.discover(PROJECT_ROOT)
    remove_lifecycle_view(
        paths.data / "raw" / "local" / DATASET_ID,
        ("legacy_source",),
    )
    remove_lifecycle_view(
        paths.data / "processed" / "per_dataset" / DATASET_ID,
        ("current.h5ad",),
    )
    dataset_root = paths.dataset(DATASET_ID)
    leftovers = remove_generated_dataset_paths(dataset_root)
    unexpected = [
        path for path in leftovers if not path.startswith(f"data/datasets/{DATASET_ID}/workspace")
    ]
    if unexpected:
        raise RuntimeError(
            "Refusing rollback because unexpected generated paths remain:\n"
            + "\n".join(unexpected)
        )
    source = Path(plan["source_workspace"])
    destination = Path(plan["canonical_workspace"])
    status = rollback_absolute_move(
        source,
        destination,
        expected_device=plan["workspace_device"],
        expected_inode=plan["workspace_inode"],
    )
    try:
        dataset_root.rmdir()
    except OSError:
        pass
    restore_registry_snapshot(plan["registry_snapshot"])
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
    print(f"{status}: {destination} -> {source}")
    return 0


def main() -> int:
    args = parse_args()
    state_dir = migration_dir(args.migration_id)
    source = args.source_workspace.expanduser().absolute()
    if args.command == "plan":
        return command_plan(state_dir, source)
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
