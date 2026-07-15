#!/usr/bin/env python3
"""Create a non-destructive repository and scientific-file checkpoint."""

from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from tnk_atlas.paths import ProjectPaths  # noqa: E402
from tnk_atlas.provenance import atomic_write_json, sampled_sha256, sha256_file  # noqa: E402


SKIP_DIR_NAMES = {".git", "__pycache__", ".pytest_cache", ".mypy_cache"}
CODE_SUFFIXES = {".py", ".r", ".sh", ".md", ".toml", ".yaml", ".yml", ".json"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint-id", help="Stable checkpoint ID; defaults to a timestamp.")
    parser.add_argument(
        "--max-hash-mb",
        type=int,
        default=64,
        help="Fully hash files no larger than this threshold (default: 64 MB).",
    )
    parser.add_argument(
        "--full-hash-large-files",
        action="store_true",
        help="Compute full SHA256 hashes for large files. This can require substantial I/O.",
    )
    return parser.parse_args()


def run_git(*args: str, text: bool = True) -> str | bytes:
    result = subprocess.run(
        ["git", *args],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=text,
    )
    return result.stdout


def safe_command(command: list[str]) -> dict[str, Any]:
    try:
        result = subprocess.run(command, check=False, capture_output=True, text=True)
        return {
            "command": command,
            "returncode": result.returncode,
            "stdout": result.stdout,
            "stderr": result.stderr,
        }
    except OSError as error:
        return {"command": command, "returncode": None, "stdout": "", "stderr": str(error)}


def classify(relative: Path) -> str:
    parts = relative.parts
    name = relative.name
    suffix = relative.suffix.lower()
    if not parts:
        return "unknown"
    if parts[0] == "archive" or parts[0] in {"old_verson", "oldversion"}:
        return "historical"
    if parts[0] == "Integrated_dataset":
        if suffix == ".h5ad":
            return "scientific_milestone_or_derived_h5ad"
        if "models" in parts:
            return "model_artifact"
        return "generated_output"
    if parts[0] in {"downloads", "newdata"}:
        if suffix in {".py", ".r", ".sh"}:
            return "legacy_code_inside_data_tree"
        return "legacy_data_tree"
    if parts[0] == "analysis_26GSE_V4":
        if "outputs" in parts or "scanpy_projects" in parts:
            return "preintegration_processed_or_interim"
        return "legacy_preintegration_workflow"
    if parts[0] in {"src", "workflows", "tests", "tools"}:
        return "active_code"
    if parts[0] in {"configs", "data"}:
        return "configuration_or_registry"
    if parts[0] in {"docs", "maintenance"}:
        return "project_documentation_or_checkpoint"
    if parts[0] in {"gdT_prediction", "gdT_atlas", "reports"}:
        return "served_report"
    if len(parts) == 1 and suffix in CODE_SUFFIXES:
        return "root_active_or_unclassified"
    if name.startswith("~$") or name.startswith(".nfs") or name.endswith(".tmp"):
        return "temporary"
    return "unclassified"


def git_state_map() -> dict[str, str]:
    output = str(run_git("status", "--porcelain=v1", "--untracked-files=all"))
    states: dict[str, str] = {}
    for line in output.splitlines():
        if len(line) < 4:
            continue
        state = line[:2]
        path_text = line[3:]
        if " -> " in path_text:
            path_text = path_text.split(" -> ", 1)[1]
        states[path_text] = state
    return states


def iter_entries(checkpoint_dir: Path):
    for directory, directory_names, file_names in os.walk(PROJECT_ROOT, followlinks=False):
        base = Path(directory)
        relative_base = base.relative_to(PROJECT_ROOT)
        symlink_directories = [name for name in directory_names if (base / name).is_symlink()]
        for name in sorted(symlink_directories):
            yield base / name
        directory_names[:] = [
            name
            for name in directory_names
            if name not in SKIP_DIR_NAMES
            and not (base / name).is_symlink()
            and (base / name).resolve() != checkpoint_dir.resolve()
        ]
        parts = relative_base.parts
        if parts and parts[0] == "downloads":
            # Deep FASTQ, SRA, extracted, and Cell Ranger trees are inventoried
            # dataset-by-dataset when migrated. Keep the baseline bounded.
            if len(parts) >= 3 and parts[1].startswith("GSE"):
                directory_names[:] = []
            elif (
                len(parts) >= 2
                and not parts[1].startswith("GSE")
                and parts[1] not in {"per_gse_h5ad_with_metadata"}
            ):
                directory_names[:] = []
        if len(parts) >= 5 and parts[:2] == ("analysis_26GSE_V4", "scanpy_projects"):
            directory_names[:] = []
        for name in sorted(file_names):
            path = base / name
            if checkpoint_dir in path.parents:
                continue
            yield path


def file_record(
    path: Path,
    states: dict[str, str],
    max_hash_bytes: int,
    full_hash_large: bool,
) -> dict[str, Any]:
    relative = path.relative_to(PROJECT_ROOT)
    stat = path.lstat()
    record: dict[str, Any] = {
        "path": relative.as_posix(),
        "kind": "symlink" if path.is_symlink() else "file",
        "classification": classify(relative),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "device": stat.st_dev,
        "inode": stat.st_ino,
        "git_state": states.get(relative.as_posix(), "tracked_or_ignored"),
        "link_target": os.readlink(path) if path.is_symlink() else "",
        "checksum_type": "none",
        "checksum": "",
    }
    if path.is_symlink() or not path.is_file():
        return record
    exact_hash_classes = {
        "active_code",
        "configuration_or_registry",
        "project_documentation_or_checkpoint",
        "root_active_or_unclassified",
    }
    should_hash = (
        stat.st_size <= max_hash_bytes
        and (
            record["classification"] in exact_hash_classes
            or record["classification"] == "model_artifact"
        )
    )
    if should_hash or full_hash_large:
        record["checksum_type"] = "sha256"
        record["checksum"] = sha256_file(path)
    elif record["classification"] in {
        "legacy_data_tree",
        "preintegration_processed_or_interim",
        "scientific_milestone_or_derived_h5ad",
    }:
        record["checksum_type"] = "deferred_dataset_migration_checksum"
    else:
        record["checksum_type"] = "stat_only"
    return record


def hdf5_structure(path: Path) -> dict[str, Any]:
    payload: dict[str, Any] = {"path": str(path), "exists": path.exists()}
    if not path.exists():
        return payload
    stat = path.stat()
    payload.update(
        size_bytes=stat.st_size,
        mtime_ns=stat.st_mtime_ns,
        device=stat.st_dev,
        inode=stat.st_ino,
        sampled_sha256=sampled_sha256(path),
    )
    try:
        import h5py

        with h5py.File(path, "r") as handle:
            payload["top_level_keys"] = sorted(handle.keys())
            payload["obs_columns"] = sorted(handle["obs"].keys()) if "obs" in handle else []
            payload["var_columns"] = sorted(handle["var"].keys()) if "var" in handle else []
            payload["layer_keys"] = sorted(handle["layers"].keys()) if "layers" in handle else []
            payload["obsm_keys"] = sorted(handle["obsm"].keys()) if "obsm" in handle else []
            payload["uns_keys"] = sorted(handle["uns"].keys()) if "uns" in handle else []
            payload["datasets"] = {}
            for key in ("obs/_index", "var/_index", "X/data", "X/indices", "X/indptr"):
                if key in handle:
                    dataset = handle[key]
                    payload["datasets"][key] = {
                        "shape": list(dataset.shape),
                        "dtype": str(dataset.dtype),
                    }
            if "X" in handle:
                payload["x_encoding"] = {
                    str(key): _jsonable(value) for key, value in handle["X"].attrs.items()
                }
    except Exception as error:  # checkpoint creation must retain the error as evidence
        payload["read_error"] = f"{type(error).__name__}: {error}"
    return payload


def _jsonable(value: Any) -> Any:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if hasattr(value, "tolist"):
        return value.tolist()
    try:
        json.dumps(value)
        return value
    except TypeError:
        return repr(value)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    fieldnames = list(rows[0]) if rows else ["path"]
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def main() -> int:
    args = parse_args()
    paths = ProjectPaths.discover(PROJECT_ROOT)
    checkpoint_id = args.checkpoint_id or datetime.now().astimezone().strftime("%Y%m%dT%H%M%S%z")
    checkpoint_dir = paths.logs / "project_reorganization" / "checkpoints" / checkpoint_id
    if checkpoint_dir.exists() and any(checkpoint_dir.iterdir()):
        raise FileExistsError(f"Checkpoint already exists and is not empty: {checkpoint_dir}")
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    states = git_state_map()
    records = [
        file_record(
            path,
            states,
            max_hash_bytes=args.max_hash_mb * 1024 * 1024,
            full_hash_large=args.full_hash_large_files,
        )
        for path in iter_entries(checkpoint_dir)
    ]
    records.sort(key=lambda row: row["path"])
    write_csv(checkpoint_dir / "filesystem_manifest.csv", records)

    canonical_h5ads = [
        paths.outputs / "TNK_candidates.h5ad",
        paths.outputs / "TNK_cleaned.h5ad",
        paths.integrated_h5ad,
        paths.historical_five_million_atlas,
    ]
    atomic_write_json(
        checkpoint_dir / "h5ad_structures.json",
        [hdf5_structure(path) for path in canonical_h5ads],
    )

    (checkpoint_dir / "git_status.txt").write_text(
        str(run_git("status", "--short", "--branch", "--untracked-files=all")), encoding="utf-8"
    )
    (checkpoint_dir / "git_tracked_files.txt").write_text(
        str(run_git("ls-files", "-s")), encoding="utf-8"
    )
    (checkpoint_dir / "tracked_changes.patch").write_bytes(
        bytes(run_git("diff", "--binary", text=False))
    )
    (checkpoint_dir / "staged_changes.patch").write_bytes(
        bytes(run_git("diff", "--cached", "--binary", text=False))
    )

    git_status_lines = (checkpoint_dir / "git_status.txt").read_text(encoding="utf-8").splitlines()
    metadata = {
        "checkpoint_id": checkpoint_id,
        "created_at": datetime.now().astimezone().isoformat(),
        "project_root": str(PROJECT_ROOT),
        "git_head": str(run_git("rev-parse", "HEAD")).strip(),
        "git_branch": str(run_git("branch", "--show-current")).strip(),
        "git_status_entry_count": max(0, len(git_status_lines) - 1),
        "python": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "max_hash_mb": args.max_hash_mb,
        "full_hash_large_files": args.full_hash_large_files,
        "manifest_entry_count": len(records),
        "scientific_files_modified": False,
    }
    atomic_write_json(checkpoint_dir / "checkpoint.json", metadata)
    atomic_write_json(
        checkpoint_dir / "environment.json",
        {
            "python": metadata,
            "conda_explicit": safe_command(["conda", "list", "--explicit", "-p", sys.prefix]),
            "conda_json": safe_command(["conda", "list", "--json", "-p", sys.prefix]),
            "pip_freeze": safe_command([sys.executable, "-m", "pip", "freeze"]),
        },
    )
    latest = checkpoint_dir.parent / "LATEST"
    latest.write_text(checkpoint_id + "\n", encoding="utf-8")
    print(checkpoint_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
