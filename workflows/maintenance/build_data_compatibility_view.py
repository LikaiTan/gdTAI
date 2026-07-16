#!/usr/bin/env python3
"""Create a reversible canonical data view using symlinks to legacy locations."""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from tnk_atlas.paths import ProjectPaths  # noqa: E402
from tnk_atlas.registry import read_csv  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", default="configs/datasets/datasets.csv")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--apply", action="store_true", help="Create missing links.")
    mode.add_argument("--verify", action="store_true", help="Verify links without modifying them.")
    parser.add_argument(
        "--replace-generated-links",
        action="store_true",
        help="Replace only existing symlinks owned by this generated view.",
    )
    return parser.parse_args()


def resolve_legacy(value: str) -> Path | None:
    if not value:
        return None
    path = Path(value)
    return path if path.is_absolute() else PROJECT_ROOT / path


def planned_links(paths: ProjectPaths, rows: list[dict[str, str]]) -> list[tuple[Path, Path, str]]:
    links: list[tuple[Path, Path, str]] = []
    for row in rows:
        dataset_id = row["dataset_id"]
        canonical_root = paths.dataset(dataset_id)
        canonical_raw = paths.dataset_raw(dataset_id)
        raw = canonical_raw if canonical_raw.exists() else resolve_legacy(row["legacy_raw_path"])
        if raw and raw.exists():
            raw_base = paths.data / "raw" / ("geo" if dataset_id.startswith("GSE") else "local") / dataset_id
            legacy_source = canonical_raw / "legacy_source"
            source = legacy_source if legacy_source.exists() else canonical_raw / "source"
            if not source.exists():
                source = raw
            links.append((raw_base / "legacy_source", source, "raw_legacy_source"))
            for legacy_name, canonical_name in (
                ("matrix", "geo_matrix"),
                ("suppl", "supplementary"),
                ("metadata", "metadata"),
                ("fastq", "fastq_legacy_mixed"),
            ):
                candidate = source / legacy_name
                if candidate.exists():
                    links.append((raw_base / canonical_name, candidate, f"raw_{canonical_name}"))
        canonical_current = paths.dataset_current_h5ad(dataset_id)
        processed = (
            canonical_current
            if canonical_current.exists()
            else resolve_legacy(row["processed_h5ad_path"])
        )
        if processed and processed.exists():
            links.append(
                (
                    paths.data / "processed" / "per_dataset" / dataset_id / "current.h5ad",
                    processed,
                    "processed_current_h5ad",
                )
            )
        canonical_scanpy = canonical_root / "interim" / "scanpy_project"
        scanpy_project = (
            canonical_scanpy
            if canonical_scanpy.exists()
            else paths.legacy_scanpy_projects / dataset_id
        )
        if scanpy_project.exists():
            links.append(
                (
                    paths.data / "interim" / dataset_id / "legacy_scanpy_project",
                    scanpy_project,
                    "interim_legacy_scanpy_project",
                )
            )
    return links


def create_link(link: Path, target: Path, replace_generated: bool) -> str:
    link.parent.mkdir(parents=True, exist_ok=True)
    relative_target = os.path.relpath(target, start=link.parent)
    if link.is_symlink():
        if link.resolve() == target.resolve() and not replace_generated:
            return "already_correct"
        if not replace_generated:
            raise FileExistsError(f"Refusing to replace mismatched symlink: {link}")
        link.unlink()
    if link.exists():
        raise FileExistsError(f"Refusing to replace existing path: {link}")
    link.symlink_to(relative_target, target_is_directory=target.is_dir())
    return "created"


def verify_link(link: Path, target: Path) -> str:
    if not link.is_symlink():
        return "missing" if not link.exists() else "not_a_symlink"
    return "correct" if link.resolve() == target.resolve() else "wrong_target"


def main() -> int:
    args = parse_args()
    paths = ProjectPaths.discover(PROJECT_ROOT)
    rows = read_csv(PROJECT_ROOT / args.registry)
    links = planned_links(paths, rows)
    manifest: list[dict[str, str]] = []
    for link, target, role in links:
        status = "planned"
        if args.apply:
            status = create_link(link, target, args.replace_generated_links)
        elif args.verify:
            status = verify_link(link, target)
        manifest.append(
            {
                "role": role,
                "link_path": link.relative_to(PROJECT_ROOT).as_posix(),
                "target_path": target.relative_to(PROJECT_ROOT).as_posix()
                if target.is_relative_to(PROJECT_ROOT)
                else str(target),
                "status": status,
            }
        )
    manifest_path = paths.data / "registry" / "compatibility_links.json"
    if args.apply:
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    mode = "apply" if args.apply else "verify" if args.verify else "dry-run"
    counts: dict[str, int] = {}
    for row in manifest:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    print(f"links={len(links)} mode={mode} statuses={json.dumps(counts, sort_keys=True)}")
    if args.verify and any(row["status"] != "correct" for row in manifest):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
