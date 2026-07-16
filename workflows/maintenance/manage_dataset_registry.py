#!/usr/bin/env python3
"""Validate, snapshot, register, or deactivate pre-integration datasets."""

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

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from tnk_atlas.provenance import sha256_file  # noqa: E402
from tnk_atlas.registry import (  # noqa: E402
    DATASET_COLUMNS,
    LIBRARY_COLUMNS,
    read_csv,
    validate_dataset_registry,
)

REGISTRY_DIR = PROJECT_ROOT / "configs" / "datasets"
DATASETS = REGISTRY_DIR / "datasets.csv"
LIBRARIES = REGISTRY_DIR / "libraries.csv"
FILES = REGISTRY_DIR / "files.csv"
SNAPSHOTS = PROJECT_ROOT / "data" / "registry" / "snapshots"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate = subparsers.add_parser("validate")
    validate.add_argument("--strict-libraries", action="store_true")

    subparsers.add_parser("list")

    locate = subparsers.add_parser("locate")
    locate.add_argument("dataset_id")

    snapshot = subparsers.add_parser("snapshot")
    snapshot.add_argument("run_id")

    restore = subparsers.add_parser("restore")
    restore.add_argument("run_id")
    restore.add_argument(
        "--confirm",
        action="store_true",
        help="Required acknowledgement before replacing the live registry files.",
    )

    register = subparsers.add_parser("register")
    register.add_argument("dataset_id")
    register.add_argument("--source-type", choices=["GEO", "local_external"], required=True)
    register.add_argument("--raw-path", required=True)
    register.add_argument("--processed-h5ad", default="")
    register.add_argument("--integration-role", default="available_not_selected")
    register.add_argument("--activate", action="store_true")
    register.add_argument("--notes", default="")

    deactivate = subparsers.add_parser("deactivate")
    deactivate.add_argument("dataset_id")
    deactivate.add_argument("--reason", required=True)
    return parser.parse_args()


def write_csv_atomic(path: Path, fieldnames: tuple[str, ...], rows: list[dict[str, str]]) -> None:
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


def git_head() -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=PROJECT_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def snapshot(run_id: str) -> Path:
    destination = SNAPSHOTS / run_id
    if destination.exists():
        raise FileExistsError(f"Snapshot already exists: {destination}")
    destination.mkdir(parents=True)
    manifest: dict[str, object] = {
        "run_id": run_id,
        "created_at": datetime.now().astimezone().isoformat(),
        "git_head": git_head(),
        "files": {},
    }
    for path in (
        DATASETS,
        LIBRARIES,
        FILES,
        REGISTRY_DIR / "integration_inputs.csv",
        REGISTRY_DIR / "sample_information_final_full.csv",
    ):
        if not path.exists():
            continue
        target = destination / path.name
        shutil.copy2(path, target)
        manifest["files"][path.name] = {
            "sha256": sha256_file(target),
            "size_bytes": target.stat().st_size,
        }
    (destination / "snapshot.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return destination


def restore_snapshot(run_id: str) -> Path:
    source = SNAPSHOTS / run_id
    manifest_path = source / "snapshot.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Registry snapshot is incomplete: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    allowed = {
        DATASETS.name,
        LIBRARIES.name,
        FILES.name,
        "integration_inputs.csv",
        "sample_information_final_full.csv",
    }
    for name, metadata in manifest.get("files", {}).items():
        if name not in allowed:
            raise ValueError(f"Snapshot contains an unsupported registry file: {name}")
        snapshot_file = source / name
        expected = metadata["sha256"]
        observed = sha256_file(snapshot_file)
        if observed != expected:
            raise RuntimeError(
                f"Snapshot checksum mismatch for {snapshot_file}: {observed} != {expected}"
            )

    backup_id = f"before_registry_restore_{datetime.now().astimezone().strftime('%Y%m%dT%H%M%S%z')}"
    backup = snapshot(backup_id)
    for name in manifest["files"]:
        source_file = source / name
        destination = REGISTRY_DIR / name
        temporary = destination.with_suffix(destination.suffix + ".restore.tmp")
        shutil.copy2(source_file, temporary)
        os.replace(temporary, destination)
    return backup


def validate(strict_libraries: bool) -> int:
    result = validate_dataset_registry(DATASETS, PROJECT_ROOT)
    datasets = read_csv(DATASETS)
    libraries = read_csv(LIBRARIES)
    dataset_ids = {row["dataset_id"] for row in datasets}
    library_dataset_ids = {row["dataset_id"] for row in libraries}
    missing_library_rows = sorted(dataset_ids - library_dataset_ids)
    placeholders = [
        row["dataset_id"]
        for row in libraries
        if row["library_id"] == "__dataset_level_pending_curation__"
        and row["active"].lower() == "true"
    ]
    payload = {
        "ok": result.ok and not missing_library_rows and not (strict_libraries and placeholders),
        "dataset_count": result.dataset_count,
        "active_dataset_count": result.active_dataset_count,
        "missing_processed_inputs": list(result.missing_processed_inputs),
        "duplicate_dataset_ids": list(result.duplicate_dataset_ids),
        "missing_library_rows": missing_library_rows,
        "active_placeholder_library_rows": sorted(placeholders),
        "strict_libraries": strict_libraries,
    }
    print(json.dumps(payload, indent=2))
    return 0 if payload["ok"] else 1


def normalize_project_path(value: str) -> str:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    if not path.exists():
        raise FileNotFoundError(path)
    try:
        return path.resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def register_dataset(args: argparse.Namespace) -> None:
    rows = read_csv(DATASETS)
    if any(row["dataset_id"] == args.dataset_id for row in rows):
        raise ValueError(f"Dataset already exists: {args.dataset_id}")
    raw_path = normalize_project_path(args.raw_path)
    processed = normalize_project_path(args.processed_h5ad) if args.processed_h5ad else ""
    if args.activate and not processed:
        raise ValueError("An active dataset requires --processed-h5ad")
    rows.append(
        {
            "dataset_id": args.dataset_id,
            "source_type": args.source_type,
            "accession": args.dataset_id,
            "phase0_active": str(args.activate).lower(),
            "current_milestone_active": "false",
            "current_cell_count": "0",
            "extended_atlas_active": "false",
            "extended_atlas_cell_count": "0",
            "integration_role": args.integration_role,
            "legacy_raw_path": raw_path,
            "processed_h5ad_path": processed,
            "status": "selected_for_next_integration" if args.activate else "available",
            "exclusion_reason": "",
            "notes": args.notes,
        }
    )
    rows.sort(key=lambda row: row["dataset_id"])
    write_csv_atomic(DATASETS, DATASET_COLUMNS, rows)

    libraries = read_csv(LIBRARIES)
    libraries.append(
        {
            "dataset_id": args.dataset_id,
            "library_id": "__dataset_level_pending_curation__",
            "sample_id": "",
            "rna_assay": "unknown",
            "tcr_assay": "unknown",
            "tcr_scope": "unknown",
            "active": str(args.activate).lower(),
            "notes": "Must be curated before integration.",
        }
    )
    libraries.sort(key=lambda row: (row["dataset_id"], row["library_id"]))
    write_csv_atomic(LIBRARIES, LIBRARY_COLUMNS, libraries)


def deactivate_dataset(dataset_id: str, reason: str) -> None:
    rows = read_csv(DATASETS)
    matched = False
    for row in rows:
        if row["dataset_id"] == dataset_id:
            matched = True
            row["phase0_active"] = "false"
            row["status"] = "excluded_from_next_integration"
            row["exclusion_reason"] = reason
    if not matched:
        raise KeyError(dataset_id)
    write_csv_atomic(DATASETS, DATASET_COLUMNS, rows)

    libraries = read_csv(LIBRARIES)
    for row in libraries:
        if row["dataset_id"] == dataset_id:
            row["active"] = "false"
    write_csv_atomic(LIBRARIES, LIBRARY_COLUMNS, libraries)


def main() -> int:
    args = parse_args()
    if args.command == "validate":
        return validate(args.strict_libraries)
    if args.command == "list":
        for row in read_csv(DATASETS):
            print(
                "\t".join(
                    [
                        row["dataset_id"],
                        row["phase0_active"],
                        row["current_milestone_active"],
                        row["extended_atlas_active"],
                        row["status"],
                    ]
                )
            )
        return 0
    if args.command == "locate":
        matches = [
            row for row in read_csv(DATASETS) if row["dataset_id"] == args.dataset_id
        ]
        if not matches:
            raise KeyError(args.dataset_id)
        row = matches[0]
        canonical = PROJECT_ROOT / "data" / "datasets" / args.dataset_id
        current_value = row["processed_h5ad_path"]
        current = PROJECT_ROOT / current_value if current_value else None
        payload = {
            "dataset_id": args.dataset_id,
            "canonical_dataset_path": str(canonical),
            "raw_path": str(canonical / "raw"),
            "interim_path": str(canonical / "interim"),
            "processed_path": str(canonical / "processed"),
            "current_h5ad": str(current) if current else "",
            "current_h5ad_target": str(current.resolve()) if current and current.exists() else "",
            "legacy_raw_path": str(PROJECT_ROOT / row["legacy_raw_path"])
            if row["legacy_raw_path"]
            else "",
        }
        print(json.dumps(payload, indent=2))
        return 0
    if args.command == "snapshot":
        print(snapshot(args.run_id))
        return 0
    if args.command == "restore":
        if not args.confirm:
            raise ValueError("Refusing registry restore without --confirm")
        backup = restore_snapshot(args.run_id)
        print(f"restored={args.run_id} pre_restore_backup={backup}")
        return validate(strict_libraries=True)

    edit_snapshot = f"before_registry_edit_{datetime.now().astimezone().strftime('%Y%m%dT%H%M%S%z')}"
    snapshot(edit_snapshot)
    if args.command == "register":
        register_dataset(args)
    elif args.command == "deactivate":
        deactivate_dataset(args.dataset_id, args.reason)
    return validate(strict_libraries=False)


if __name__ == "__main__":
    raise SystemExit(main())
