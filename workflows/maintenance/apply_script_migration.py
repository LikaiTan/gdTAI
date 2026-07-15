#!/usr/bin/env python3
"""Apply or roll back the checksum-guarded root-script migration plan."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import sys
from datetime import datetime
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from tnk_atlas.provenance import sha256_file  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--plan",
        default="maintenance/reorganization/script_migration_plan.csv",
    )
    parser.add_argument("--output", help="Migration result CSV; derived from the plan by default.")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--apply", action="store_true")
    mode.add_argument("--rollback", action="store_true")
    return parser.parse_args()


def read_plan(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def content_sha256(path: Path) -> str:
    if path.is_symlink():
        return hashlib.sha256(("symlink\0" + os.readlink(path)).encode("utf-8")).hexdigest()
    if path.is_file():
        return sha256_file(path)
    if path.is_dir():
        digest = hashlib.sha256()
        for child in sorted(path.rglob("*"), key=lambda item: item.as_posix()):
            relative = child.relative_to(path).as_posix()
            if child.is_symlink():
                digest.update(f"L\0{relative}\0{os.readlink(child)}\n".encode("utf-8"))
            elif child.is_dir():
                digest.update(f"D\0{relative}\n".encode("utf-8"))
            elif child.is_file():
                digest.update(f"F\0{relative}\0{child.stat().st_size}\0".encode("utf-8"))
                digest.update(bytes.fromhex(sha256_file(child)))
        return digest.hexdigest()
    raise FileNotFoundError(path)


def move_checked(source: Path, destination: Path) -> tuple[str, str]:
    if source.exists() and destination.exists():
        raise FileExistsError(f"Both source and destination exist: {source}, {destination}")
    if not source.exists() and destination.exists():
        return "already_moved", content_sha256(destination)
    if not source.exists():
        return "missing_source", ""
    checksum = content_sha256(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    os.replace(source, destination)
    if content_sha256(destination) != checksum:
        raise RuntimeError(f"Checksum changed while moving {source} to {destination}")
    return "moved", checksum


def main() -> int:
    args = parse_args()
    rows = read_plan(PROJECT_ROOT / args.plan)
    results: list[dict[str, str]] = []
    for row in rows:
        old = PROJECT_ROOT / row["old_path"]
        new = PROJECT_ROOT / row["new_path"]
        source, destination = (new, old) if args.rollback else (old, new)
        if args.apply or args.rollback:
            status, checksum = move_checked(source, destination)
        else:
            status = "planned" if source.exists() else "source_missing"
            checksum = content_sha256(source) if source.exists() else ""
        results.append(
            {
                **row,
                "operation": "rollback" if args.rollback else "apply",
                "status": status,
                "sha256": checksum,
                "recorded_at": datetime.now().astimezone().isoformat(),
            }
        )

    plan_path = Path(args.plan)
    if args.output:
        output = PROJECT_ROOT / args.output
    elif plan_path.name == "script_migration_plan.csv":
        output = PROJECT_ROOT / "maintenance" / "reorganization" / "migration_map.csv"
    else:
        output_name = plan_path.stem.removesuffix("_plan") + "_map.csv"
        output = PROJECT_ROOT / "maintenance" / "reorganization" / output_name
    temporary = output.with_suffix(".csv.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, output)
    counts: dict[str, int] = {}
    for row in results:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    print(counts)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
