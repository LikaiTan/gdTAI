#!/usr/bin/env python3
"""Migrate legacy single-cell roots into dataset-centered canonical storage."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from tnk_atlas.paths import ProjectPaths  # noqa: E402
from tnk_atlas.provenance import atomic_write_json, sampled_sha256  # noqa: E402
from tnk_atlas.registry import DATASET_COLUMNS, FILE_COLUMNS, read_csv  # noqa: E402
from tnk_atlas.storage import (  # noqa: E402
    StorageMove,
    apply_move,
    create_compatibility_link,
    read_plan,
    rollback_move,
    translate_path,
    validate_move,
    write_plan,
)

DEFAULT_MIGRATION_ID = "dataset_centered_20260716"
DEFAULT_REGISTRY_SNAPSHOT = "pre_dataset_centered_migration_20260716"
GSE_PATTERN = re.compile(r"^GSE\d+$")
ROOT_ALIASES = ("downloads", "analysis_26GSE_V4", "newdata")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--migration-id", default=DEFAULT_MIGRATION_ID)
    subparsers = parser.add_subparsers(dest="command", required=True)
    plan_parser = subparsers.add_parser("plan")
    plan_parser.add_argument(
        "--refresh",
        action="store_true",
        help="Regenerate an unapplied migration plan and its baseline evidence.",
    )

    subparsers.add_parser("preflight")

    apply_parser = subparsers.add_parser("apply")
    apply_parser.add_argument("--confirm", action="store_true")

    finalize = subparsers.add_parser("finalize")
    finalize.add_argument("--confirm", action="store_true")

    subparsers.add_parser("validate")

    rollback = subparsers.add_parser("rollback")
    rollback.add_argument("--confirm", action="store_true")
    rollback.add_argument(
        "--registry-snapshot",
        default=DEFAULT_REGISTRY_SNAPSHOT,
    )
    return parser.parse_args()


def migration_dir(migration_id: str) -> Path:
    return PROJECT_ROOT / "data" / "registry" / "migrations" / migration_id


def project_relative(path: Path) -> str:
    return path.absolute().relative_to(PROJECT_ROOT.absolute()).as_posix()


def dataset_ids() -> list[str]:
    rows = read_csv(PROJECT_ROOT / "configs" / "datasets" / "datasets.csv")
    return sorted({row["dataset_id"] for row in rows}, key=lambda value: (-len(value), value))


def dataset_from_name(name: str, known_ids: list[str]) -> str | None:
    normalized = name.lstrip(".")
    for dataset_id in known_ids:
        if normalized.startswith(dataset_id):
            return dataset_id
    match = re.match(r"(GSE\d+|HRA\d+)", normalized)
    return match.group(1) if match else None


def add_move(
    moves: list[StorageMove],
    *,
    old: Path,
    new: Path,
    dataset_id: str = "",
    role: str,
    operation: str = "move_and_link",
) -> None:
    if not old.exists() or old.is_symlink():
        return
    if new.exists() or new.is_symlink():
        raise FileExistsError(f"Canonical destination already exists: {new}")
    if any(move.old_path == project_relative(old) for move in moves):
        raise ValueError(f"Duplicate source in migration plan: {old}")
    if any(move.new_path == project_relative(new) for move in moves):
        raise ValueError(f"Duplicate destination in migration plan: {new}")
    moves.append(
        StorageMove.from_path(
            sequence=len(moves) + 1,
            operation=operation,
            dataset_id=dataset_id,
            role=role,
            old_path=old,
            new_path=new,
            root=PROJECT_ROOT,
        )
    )


def build_plan() -> list[StorageMove]:
    paths = ProjectPaths.discover(PROJECT_ROOT)
    known_ids = dataset_ids()
    moves: list[StorageMove] = []
    downloads = PROJECT_ROOT / "downloads"
    analysis = PROJECT_ROOT / "analysis_26GSE_V4"
    newdata = PROJECT_ROOT / "newdata"

    # Move compatibility roots first. Child operations then create relative
    # links from the roots' final physical locations under data/compat/.
    for name in ROOT_ALIASES:
        add_move(
            moves,
            old=PROJECT_ROOT / name,
            new=paths.compatibility_data / name,
            role="top_level_legacy_root_alias",
            operation="root_alias",
        )

    for old in sorted(downloads.iterdir()):
        if old.is_dir() and GSE_PATTERN.match(old.name):
            add_move(
                moves,
                old=old,
                new=paths.dataset_raw(old.name) / "legacy_source",
                dataset_id=old.name,
                role="legacy_download_source_bundle",
            )

    scanpy_root = analysis / "scanpy_projects"
    if scanpy_root.exists():
        for old in sorted(scanpy_root.iterdir()):
            if old.is_dir() and GSE_PATTERN.match(old.name):
                add_move(
                    moves,
                    old=old,
                    new=paths.dataset_interim(old.name) / "scanpy_project",
                    dataset_id=old.name,
                    role="legacy_scanpy_project",
                )

    per_gse = downloads / "per_gse_h5ad_with_metadata"
    if per_gse.exists():
        for old in sorted(per_gse.iterdir()):
            if old.is_dir():
                destination = paths.shared_data / "legacy_processed_h5ads" / old.name
                dataset_id = ""
                role = "unclassified_processed_directory"
            else:
                dataset_id = dataset_from_name(old.name, known_ids) or ""
                if ".tmp" in old.name or old.name.endswith(".partial"):
                    bucket = "quarantine"
                    role = "incomplete_processed_artifact"
                elif ".bak" in old.name or "wrong_" in old.name:
                    bucket = "archive"
                    role = "historical_processed_artifact"
                else:
                    bucket = "artifacts"
                    role = "processed_h5ad_artifact"
                if dataset_id:
                    destination = paths.dataset_processed(dataset_id) / bucket / old.name
                else:
                    destination = paths.shared_data / "legacy_processed_h5ads" / old.name
                    role = "unclassified_processed_artifact"
            add_move(
                moves,
                old=old,
                new=destination,
                dataset_id=dataset_id,
                role=role,
            )

    hra = newdata / "HRA005041"
    add_move(
        moves,
        old=hra,
        new=paths.dataset_raw("HRA005041") / "legacy_source",
        dataset_id="HRA005041",
        role="local_source_bundle",
    )

    sorted_root = newdata / "Sorted_gdT"
    if sorted_root.exists():
        for old in sorted(sorted_root.iterdir()):
            dataset_id = dataset_from_name(old.name, known_ids)
            if dataset_id:
                if old.suffix.lower() == ".h5ad":
                    new = paths.dataset_processed(dataset_id) / "artifacts" / old.name
                    role = "sorted_gdt_processed_h5ad"
                else:
                    new = paths.dataset_raw(dataset_id) / "source" / old.name
                    role = "sorted_gdt_source"
            else:
                new = paths.shared_data / "legacy_newdata" / "Sorted_gdT" / old.name
                role = "unclassified_sorted_gdt_artifact"
            add_move(
                moves,
                old=old,
                new=new,
                dataset_id=dataset_id or "",
                role=role,
            )

    for old in sorted(newdata.iterdir()):
        if old.name in {"HRA005041", "Sorted_gdT"}:
            continue
        dataset_id = dataset_from_name(old.name, known_ids)
        if dataset_id:
            new = paths.dataset_processed(dataset_id) / "artifacts" / old.name
            role = (
                "standalone_processed_h5ad"
                if old.suffix.lower() == ".h5ad"
                else "standalone_artifact"
            )
        else:
            new = paths.shared_data / "legacy_newdata" / old.name
            role = "unclassified_newdata_artifact"
        add_move(
            moves,
            old=old,
            new=new,
            dataset_id=dataset_id or "",
            role=role,
        )

    for old in sorted(downloads.iterdir()):
        if (
            (old.is_dir() and GSE_PATTERN.match(old.name))
            or old.name == "per_gse_h5ad_with_metadata"
        ):
            continue
        if old.name.startswith(("scanpy_integration_output", "tcr_integration_output")):
            new = paths.shared_data / "legacy_integration_outputs" / old.name
            role = "legacy_integration_output"
        else:
            new = paths.shared_data / "legacy_download_root" / old.name
            role = "legacy_download_root_artifact"
        add_move(moves, old=old, new=new, role=role)

    for old in sorted(analysis.iterdir()):
        if old.name == "scanpy_projects":
            continue
        add_move(
            moves,
            old=old,
            new=paths.shared_data / "preintegration" / "analysis_26GSE_V4" / old.name,
            role="legacy_preintegration_shared_artifact",
        )

    return moves


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_json(path, payload)


def h5ad_axis_length(handle: Any, axis: str) -> int:
    frame = handle[axis]
    index_name = frame.attrs.get("_index", "_index")
    if isinstance(index_name, bytes):
        index_name = index_name.decode("utf-8")
    if index_name in frame:
        return int(frame[index_name].shape[0])
    for node in frame.values():
        if hasattr(node, "shape") and node.shape:
            return int(node.shape[0])
    raise KeyError(f"No index-like dataset found in H5AD {axis!r} group")


def h5ad_record(path: Path, dataset_id: str) -> dict[str, Any]:
    stat = path.stat()
    payload: dict[str, Any] = {
        "dataset_id": dataset_id,
        "path": project_relative(path),
        "device": stat.st_dev,
        "inode": stat.st_ino,
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sampled_sha256": sampled_sha256(path),
    }
    try:
        import h5py

        with h5py.File(path, "r") as handle:
            payload["top_level_keys"] = sorted(handle.keys())
            payload["n_obs"] = h5ad_axis_length(handle, "obs")
            payload["n_vars"] = h5ad_axis_length(handle, "var")
    except Exception as error:
        payload["read_error"] = f"{type(error).__name__}: {error}"
    return payload


def pre_migration_selected_h5ads() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for row in read_csv(PROJECT_ROOT / "configs" / "datasets" / "datasets.csv"):
        value = row["processed_h5ad_path"].strip()
        if not value:
            continue
        path = Path(value)
        if not path.is_absolute():
            path = PROJECT_ROOT / path
        if path.exists():
            records.append(h5ad_record(path, row["dataset_id"]))
    return records


def snapshot_configuration(state_dir: Path) -> None:
    snapshot = state_dir / "pre_migration_configuration"
    snapshot.mkdir(parents=True, exist_ok=True)
    for relative in (
        "configs/datasets/datasets.csv",
        "configs/datasets/files.csv",
        "configs/datasets/integration_inputs.csv",
        "data/registry/compatibility_links.json",
    ):
        source = PROJECT_ROOT / relative
        if source.exists():
            shutil.copy2(source, snapshot / source.name)


def command_plan(state_dir: Path, refresh: bool) -> int:
    state_dir.mkdir(parents=True, exist_ok=True)
    plan_path = state_dir / "move_plan.csv"
    if plan_path.exists() and not refresh:
        raise FileExistsError(f"Migration plan already exists: {plan_path}")
    if refresh and any(
        (state_dir / name).exists()
        for name in ("apply_journal.jsonl", "apply_summary.json", "finalize_summary.json")
    ):
        raise RuntimeError("Refusing to refresh a migration plan after apply has started")
    if refresh and any((PROJECT_ROOT / name).is_symlink() for name in ROOT_ALIASES):
        raise RuntimeError("Refusing to refresh a migration plan after physical cutover")
    moves = build_plan()
    write_plan(plan_path, moves)
    snapshot_configuration(state_dir)
    write_json(state_dir / "selected_h5ad_baseline.json", pre_migration_selected_h5ads())
    write_json(
        state_dir / "plan_summary.json",
        {
            "migration_id": state_dir.name,
            "created_at": datetime.now().astimezone().isoformat(),
            "git_head": subprocess.run(
                ["git", "rev-parse", "HEAD"],
                cwd=PROJECT_ROOT,
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip(),
            "operation_count": len(moves),
            "dataset_operation_count": sum(bool(move.dataset_id) for move in moves),
            "root_aliases": list(ROOT_ALIASES),
            "same_filesystem": len({move.device for move in moves}) == 1,
            "devices": sorted({move.device for move in moves}),
        },
    )
    print(f"plan={plan_path} operations={len(moves)}")
    return 0


def nearest_existing_ancestor(path: Path) -> Path:
    candidate = path.parent
    while not candidate.exists():
        if candidate == candidate.parent:
            raise FileNotFoundError(f"No existing ancestor for destination: {path}")
        candidate = candidate.parent
    return candidate


def validate_pre_migration_h5ads(
    state_dir: Path,
    failures: list[str],
) -> None:
    baseline_path = state_dir / "selected_h5ad_baseline.json"
    if not baseline_path.exists():
        failures.append(f"selected H5AD baseline is missing: {baseline_path}")
        return
    baselines = json.loads(baseline_path.read_text(encoding="utf-8"))
    for baseline in baselines:
        if baseline.get("read_error"):
            failures.append(
                f"selected H5AD baseline has read error: {baseline['path']}: "
                f"{baseline['read_error']}"
            )
            continue
        path = PROJECT_ROOT / baseline["path"]
        if not path.exists():
            failures.append(f"selected H5AD is missing before migration: {path}")
            continue
        current = h5ad_record(path, baseline["dataset_id"])
        for field in (
            "device",
            "inode",
            "size_bytes",
            "mtime_ns",
            "sampled_sha256",
            "top_level_keys",
            "n_obs",
            "n_vars",
        ):
            if current.get(field) != baseline.get(field):
                failures.append(
                    f"selected H5AD baseline mismatch for {field}: {baseline['path']}"
                )
                break


def command_preflight(state_dir: Path) -> int:
    plan_path = state_dir / "move_plan.csv"
    if not plan_path.exists():
        raise FileNotFoundError(f"Migration plan is missing: {plan_path}")
    moves = read_plan(plan_path)
    failures: list[str] = []

    expected_sequences = list(range(1, len(moves) + 1))
    if [move.sequence for move in moves] != expected_sequences:
        failures.append("move sequence is not contiguous")
    if len({move.old_path for move in moves}) != len(moves):
        failures.append("move plan contains duplicate source paths")
    if len({move.new_path for move in moves}) != len(moves):
        failures.append("move plan contains duplicate destination paths")
    if [move.old_path for move in moves[: len(ROOT_ALIASES)]] != list(ROOT_ALIASES):
        failures.append("top-level compatibility roots are not the first operations")
    if any(move.operation != "root_alias" for move in moves[: len(ROOT_ALIASES)]):
        failures.append("initial compatibility-root operations have the wrong type")
    if any(move.operation not in {"root_alias", "move_and_link"} for move in moves):
        failures.append("move plan contains an unsupported operation")
    if len({move.device for move in moves}) != 1:
        failures.append("planned sources are not all on one filesystem")

    for move in moves:
        old = PROJECT_ROOT / move.old_path
        new = PROJECT_ROOT / move.new_path
        if old.is_symlink():
            failures.append(f"source is already a symlink: {move.old_path}")
            continue
        if not old.exists():
            failures.append(f"source is missing: {move.old_path}")
            continue
        stat = old.stat()
        if stat.st_dev != move.device or stat.st_ino != move.inode:
            failures.append(f"source inode changed since planning: {move.old_path}")
        if move.expected_kind == "file":
            if not old.is_file():
                failures.append(f"planned file is not a file: {move.old_path}")
            elif stat.st_size != move.size_bytes or stat.st_mtime_ns != move.mtime_ns:
                failures.append(f"planned file stat changed: {move.old_path}")
        elif not old.is_dir():
            failures.append(f"planned directory is not a directory: {move.old_path}")
        if new.exists() or new.is_symlink():
            failures.append(f"destination already exists: {move.new_path}")
        else:
            ancestor = nearest_existing_ancestor(new)
            if ancestor.stat().st_dev != move.device:
                failures.append(
                    f"destination is not on the source filesystem: {move.new_path}"
                )

    snapshot = (
        PROJECT_ROOT
        / "data"
        / "registry"
        / "snapshots"
        / DEFAULT_REGISTRY_SNAPSHOT
        / "snapshot.json"
    )
    if not snapshot.exists():
        failures.append(f"registry rollback snapshot is missing: {snapshot}")
    validate_pre_migration_h5ads(state_dir, failures)
    result = {
        "migration_id": state_dir.name,
        "operation_count": len(moves),
        "failure_count": len(failures),
        "failures": failures,
        "validated_at": datetime.now().astimezone().isoformat(),
    }
    write_json(state_dir / "preflight.json", result)
    print(json.dumps(result, indent=2))
    return 1 if failures else 0


def migration_has_started(moves: list[StorageMove]) -> bool:
    for move in moves:
        old = PROJECT_ROOT / move.old_path
        new = PROJECT_ROOT / move.new_path
        if old.is_symlink() or (not old.exists() and new.exists()):
            return True
    return False


def append_journal(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")
        handle.flush()
        os.fsync(handle.fileno())


def command_apply(state_dir: Path, confirmed: bool) -> int:
    if not confirmed:
        raise ValueError("Refusing physical migration without --confirm")
    moves = read_plan(state_dir / "move_plan.csv")
    if not migration_has_started(moves) and command_preflight(state_dir):
        raise RuntimeError("Migration preflight failed; no paths were changed")
    journal = state_dir / "apply_journal.jsonl"
    for move in moves:
        status = apply_move(PROJECT_ROOT, move)
        append_journal(
            journal,
            {
                "sequence": move.sequence,
                "old_path": move.old_path,
                "new_path": move.new_path,
                "status": status,
                "timestamp": datetime.now().astimezone().isoformat(),
            },
        )
        print(f"{move.sequence:03d} {status:15s} {move.old_path} -> {move.new_path}")
    write_json(
        state_dir / "apply_summary.json",
        {
            "completed_at": datetime.now().astimezone().isoformat(),
            "operation_count": len(moves),
            "status": "physical_moves_complete",
        },
    )
    return 0


def translated_selected_paths(moves: list[StorageMove]) -> dict[str, Path]:
    selected: dict[str, Path] = {}
    for row in read_csv(PROJECT_ROOT / "configs" / "datasets" / "datasets.csv"):
        value = row["processed_h5ad_path"].strip()
        if value:
            registered = Path(value)
            if not registered.is_absolute():
                registered = PROJECT_ROOT / registered
            selected[row["dataset_id"]] = (
                registered.resolve()
                if registered.is_symlink()
                else translate_path(registered, PROJECT_ROOT, moves)
            )

    fallback_candidates: dict[str, list[Path]] = {}
    for move in moves:
        if move.dataset_id and move.role == "standalone_processed_h5ad":
            fallback_candidates.setdefault(move.dataset_id, []).append(
                PROJECT_ROOT / move.new_path
            )
    for dataset_id, candidates in fallback_candidates.items():
        if dataset_id not in selected and len(candidates) == 1:
            selected[dataset_id] = candidates[0]
    return selected


def replace_generated_link(link: Path, target: Path) -> None:
    if link.is_symlink():
        link.unlink()
    elif link.exists():
        raise FileExistsError(f"Refusing to replace non-symlink generated path: {link}")
    link.parent.mkdir(parents=True, exist_ok=True)
    create_compatibility_link(link, target)


def write_dataset_metadata(
    paths: ProjectPaths,
    rows: list[dict[str, str]],
    selected: dict[str, Path],
    migration_id: str,
) -> None:
    for row in rows:
        dataset_id = row["dataset_id"]
        root = paths.dataset(dataset_id)
        root.mkdir(parents=True, exist_ok=True)
        current = paths.dataset_current_h5ad(dataset_id)
        target = selected.get(dataset_id)
        if target is not None:
            if not target.exists():
                raise FileNotFoundError(f"Translated selected H5AD is missing: {target}")
            replace_generated_link(current, target)
        write_json(
            root / "dataset.json",
            {
                "dataset_id": dataset_id,
                "source_type": row["source_type"],
                "canonical_dataset_path": project_relative(root),
                "raw_path": project_relative(paths.dataset_raw(dataset_id)),
                "interim_path": project_relative(paths.dataset_interim(dataset_id)),
                "processed_path": project_relative(paths.dataset_processed(dataset_id)),
                "current_h5ad": project_relative(current) if current.exists() else "",
                "current_h5ad_target": project_relative(target) if target is not None else "",
                "legacy_raw_path": row["legacy_raw_path"],
                "storage_migration_id": migration_id,
            },
        )


def write_csv_atomic(path: Path, fieldnames: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
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


def update_dataset_registry(
    paths: ProjectPaths,
    rows: list[dict[str, str]],
    selected: dict[str, Path],
) -> None:
    for row in rows:
        dataset_id = row["dataset_id"]
        if dataset_id in selected:
            row["processed_h5ad_path"] = project_relative(paths.dataset_current_h5ad(dataset_id))
        note = row["notes"].strip()
        storage_note = "Canonical physical storage is dataset-centered under data/datasets."
        if storage_note not in note:
            row["notes"] = f"{note} {storage_note}".strip()
    write_csv_atomic(
        PROJECT_ROOT / "configs" / "datasets" / "datasets.csv",
        DATASET_COLUMNS,
        rows,
    )


def update_integration_inputs(paths: ProjectPaths) -> None:
    source = PROJECT_ROOT / "configs" / "datasets" / "integration_inputs.csv"
    rows = read_csv(source)
    for row in rows:
        dataset_id = row["gse_id"]
        current = paths.dataset_current_h5ad(dataset_id)
        if not current.exists():
            raise FileNotFoundError(f"Active integration input is missing: {current}")
        row["h5ad_path"] = str(current.absolute())
        row["source_root"] = str(paths.dataset_processed(dataset_id).absolute())
    write_csv_atomic(source, ("h5ad_path", "gse_id", "source_root"), rows)


def update_file_registry(moves: list[StorageMove]) -> None:
    source = PROJECT_ROOT / "configs" / "datasets" / "files.csv"
    rows = read_csv(source)
    known_ids = dataset_ids()
    filtered: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for row in rows:
        legacy = Path(row["legacy_path"])
        if legacy.parts[:2] == ("newdata", "Sorted_gdT"):
            inferred = dataset_from_name(legacy.name, known_ids)
            if inferred and inferred != row["dataset_id"]:
                continue
        canonical = translate_path(legacy, PROJECT_ROOT, moves)
        key = (row["dataset_id"], project_relative(canonical))
        if key in seen:
            continue
        seen.add(key)
        row["canonical_path"] = key[1]
        row["checksum_type"] = "inode_preserved_by_same_filesystem_rename"
        filtered.append(row)
    filtered.sort(key=lambda row: (row["dataset_id"], row["lifecycle_level"], row["canonical_path"]))
    write_csv_atomic(source, FILE_COLUMNS, filtered)


def write_storage_index(paths: ProjectPaths, rows: list[dict[str, str]]) -> None:
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
        raise ValueError("Refusing registry cutover without --confirm")
    moves = read_plan(state_dir / "move_plan.csv")
    errors = [error for move in moves for error in validate_move(PROJECT_ROOT, move)]
    validate_selected_h5ads(state_dir, moves, errors)
    if errors:
        raise RuntimeError("Physical migration is not valid:\n" + "\n".join(errors[:20]))
    paths = ProjectPaths.discover(PROJECT_ROOT)
    rows = read_csv(PROJECT_ROOT / "configs" / "datasets" / "datasets.csv")
    selected = translated_selected_paths(moves)
    write_dataset_metadata(paths, rows, selected, state_dir.name)
    update_dataset_registry(paths, rows, selected)
    update_integration_inputs(paths)
    update_file_registry(moves)
    write_storage_index(paths, rows)
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
    write_json(
        state_dir / "finalize_summary.json",
        {
            "completed_at": datetime.now().astimezone().isoformat(),
            "dataset_count": len(rows),
            "selected_h5ad_count": len(selected),
            "status": "canonical_registry_cutover_complete",
        },
    )
    return 0


def validate_selected_h5ads(
    state_dir: Path,
    moves: list[StorageMove],
    failures: list[str],
) -> None:
    baselines = json.loads((state_dir / "selected_h5ad_baseline.json").read_text(encoding="utf-8"))
    for row in baselines:
        canonical = translate_path(Path(row["path"]), PROJECT_ROOT, moves)
        if not canonical.exists():
            failures.append(f"selected H5AD missing: {canonical}")
            continue
        stat = canonical.stat()
        if (
            stat.st_dev != row["device"]
            or stat.st_ino != row["inode"]
            or stat.st_size != row["size_bytes"]
            or stat.st_mtime_ns != row["mtime_ns"]
        ):
            failures.append(f"selected H5AD stat/inode mismatch: {canonical}")
            continue
        if sampled_sha256(canonical) != row["sampled_sha256"]:
            failures.append(f"selected H5AD sampled hash mismatch: {canonical}")


def command_validate(state_dir: Path) -> int:
    moves = read_plan(state_dir / "move_plan.csv")
    failures = [error for move in moves for error in validate_move(PROJECT_ROOT, move)]
    validate_selected_h5ads(state_dir, moves, failures)
    paths = ProjectPaths.discover(PROJECT_ROOT)
    rows = read_csv(PROJECT_ROOT / "configs" / "datasets" / "datasets.csv")
    for row in rows:
        dataset_id = row["dataset_id"]
        if not paths.dataset(dataset_id).exists():
            failures.append(f"canonical dataset directory missing: {dataset_id}")
        value = row["processed_h5ad_path"]
        if value and not (PROJECT_ROOT / value).exists():
            failures.append(f"registered processed H5AD missing: {dataset_id}")
        legacy = row["legacy_raw_path"]
        if legacy and not (PROJECT_ROOT / legacy).exists():
            failures.append(f"legacy raw compatibility path missing: {dataset_id}")
    for name in ROOT_ALIASES:
        if not (PROJECT_ROOT / name).is_symlink():
            failures.append(f"top-level compatibility root is not a symlink: {name}")
    result = {
        "migration_id": state_dir.name,
        "operation_count": len(moves),
        "dataset_count": len(rows),
        "failure_count": len(failures),
        "failures": failures,
        "validated_at": datetime.now().astimezone().isoformat(),
    }
    write_json(state_dir / "validation.json", result)
    print(json.dumps(result, indent=2))
    return 1 if failures else 0


def restore_registry_snapshot(snapshot_id: str) -> None:
    snapshot = PROJECT_ROOT / "data" / "registry" / "snapshots" / snapshot_id
    manifest_path = snapshot / "snapshot.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Registry snapshot is missing: {snapshot}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    allowed = {
        "datasets.csv",
        "files.csv",
        "integration_inputs.csv",
        "libraries.csv",
        "sample_information_final_full.csv",
    }
    for name in manifest["files"]:
        if name not in allowed:
            continue
        source = snapshot / name
        destination = PROJECT_ROOT / "configs" / "datasets" / name
        temporary = destination.with_suffix(destination.suffix + ".rollback.tmp")
        shutil.copy2(source, temporary)
        os.replace(temporary, destination)


def remove_empty_tree(path: Path) -> list[str]:
    leftovers: list[str] = []
    if not path.exists():
        return leftovers
    for root, directories, _files in os.walk(path, topdown=False):
        directory = Path(root)
        for name in directories:
            child = directory / name
            try:
                child.rmdir()
            except OSError:
                pass
        try:
            directory.rmdir()
        except OSError:
            pass
    if path.exists():
        leftovers.extend(
            project_relative(candidate)
            for candidate in sorted(path.rglob("*"))
            if candidate.exists() or candidate.is_symlink()
        )
    return leftovers


def cleanup_finalized_layout(state_dir: Path) -> list[str]:
    paths = ProjectPaths.discover(PROJECT_ROOT)
    if (state_dir / "finalize_summary.json").exists():
        for row in read_csv(PROJECT_ROOT / "configs" / "datasets" / "datasets.csv"):
            dataset_root = paths.dataset(row["dataset_id"])
            metadata = dataset_root / "dataset.json"
            if metadata.exists():
                payload = json.loads(metadata.read_text(encoding="utf-8"))
                if payload.get("storage_migration_id") == state_dir.name:
                    metadata.unlink()
            current = paths.dataset_current_h5ad(row["dataset_id"])
            if current.is_symlink():
                current.unlink()
        storage_index = PROJECT_ROOT / "data" / "registry" / "storage_index.csv"
        if storage_index.exists():
            storage_index.unlink()

    leftovers: list[str] = []
    for root in (paths.datasets, paths.shared_data, paths.compatibility_data):
        leftovers.extend(remove_empty_tree(root))
    return leftovers


def command_rollback(state_dir: Path, confirmed: bool, registry_snapshot: str) -> int:
    if not confirmed:
        raise ValueError("Refusing rollback without --confirm")
    moves = read_plan(state_dir / "move_plan.csv")
    journal = state_dir / "rollback_journal.jsonl"
    for move in reversed(moves):
        status = rollback_move(PROJECT_ROOT, move)
        append_journal(
            journal,
            {
                "sequence": move.sequence,
                "old_path": move.old_path,
                "new_path": move.new_path,
                "status": status,
                "timestamp": datetime.now().astimezone().isoformat(),
            },
        )
        print(f"{move.sequence:03d} {status:19s} {move.new_path} -> {move.old_path}")
    leftovers = cleanup_finalized_layout(state_dir)
    restore_registry_snapshot(registry_snapshot)
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
    write_json(
        state_dir / "rollback_summary.json",
        {
            "completed_at": datetime.now().astimezone().isoformat(),
            "leftover_paths": leftovers,
            "registry_snapshot": registry_snapshot,
            "status": "rolled_back",
        },
    )
    return 0


def main() -> int:
    args = parse_args()
    state_dir = migration_dir(args.migration_id)
    if args.command == "plan":
        return command_plan(state_dir, args.refresh)
    if args.command == "preflight":
        return command_preflight(state_dir)
    if args.command == "apply":
        return command_apply(state_dir, args.confirm)
    if args.command == "finalize":
        return command_finalize(state_dir, args.confirm)
    if args.command == "validate":
        return command_validate(state_dir)
    if args.command == "rollback":
        return command_rollback(state_dir, args.confirm, args.registry_snapshot)
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
