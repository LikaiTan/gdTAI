#!/usr/bin/env python3
"""Build the deterministic move plan for clearly historical project content."""

from __future__ import annotations

import csv
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTPUT = PROJECT_ROOT / "maintenance" / "reorganization" / "legacy_archive_plan.csv"
CODE_SUFFIXES = {".py", ".r", ".sh"}


def add(rows: list[dict[str, str]], old: Path, new: Path, category: str, reason: str) -> None:
    if not old.exists() and not new.exists():
        return
    rows.append(
        {
            "old_path": old.relative_to(PROJECT_ROOT).as_posix(),
            "new_path": new.relative_to(PROJECT_ROOT).as_posix(),
            "category": category,
            "reason": reason,
        }
    )


def main() -> int:
    rows: list[dict[str, str]] = []
    add(
        rows,
        PROJECT_ROOT / "old_verson",
        PROJECT_ROOT / "archive" / "legacy_pipeline" / "pre_canonical_data_intake",
        "legacy_directory",
        "Superseded pre-canonical intake code and records",
    )
    add(
        rows,
        PROJECT_ROOT / "oldversion",
        PROJECT_ROOT / "archive" / "legacy_pipeline" / "control_file_snapshots",
        "legacy_directory",
        "Superseded control-file snapshots",
    )

    download_archive = PROJECT_ROOT / "archive" / "legacy_pipeline" / "download_helpers"
    downloads = PROJECT_ROOT / "downloads"
    for path in sorted(downloads.iterdir() if downloads.exists() else []):
        if path.is_file() and path.suffix.lower() in CODE_SUFFIXES:
            add(
                rows,
                path,
                download_archive / path.name,
                "legacy_download_code",
                "Completed acquisition helper superseded by registry-driven intake",
            )

    embedded_archive = (
        PROJECT_ROOT
        / "archive"
        / "legacy_pipeline"
        / "dataset_embedded"
        / "GSE144469"
        / "fastq_scripts"
    )
    embedded = PROJECT_ROOT / "downloads" / "GSE144469" / "fastq"
    for path in sorted(embedded.iterdir() if embedded.exists() else []):
        if path.is_file() and path.suffix.lower() in CODE_SUFFIXES:
            add(
                rows,
                path,
                embedded_archive / path.name,
                "legacy_dataset_code",
                "Superseded exploratory code must not live inside a source-data directory",
            )

    for name in ("TNK_PHASES_1_3_SCRIPT.md",):
        add(
            rows,
            PROJECT_ROOT / name,
            PROJECT_ROOT / "archive" / "legacy_pipeline" / "control_files" / name,
            "legacy_control_file",
            "Superseded by TNK_PHASES_0_4_SCRIPT.md",
        )
    for name in ("h5ad_backup_before_v2_20260320.csv", "h5ad_v2.csv"):
        add(
            rows,
            PROJECT_ROOT / name,
            PROJECT_ROOT / "archive" / "legacy_pipeline" / "registries" / name,
            "legacy_registry",
            "Historical registry retained for provenance",
        )
    for path in sorted(PROJECT_ROOT.glob("~$*")):
        add(
            rows,
            path,
            PROJECT_ROOT / "archive" / "legacy_pipeline" / "quarantine" / path.name,
            "temporary_quarantine",
            "Office lock artifact quarantined rather than deleted",
        )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    temporary = OUTPUT.with_suffix(".csv.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("old_path", "new_path", "category", "reason"),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, OUTPUT)
    print(f"archive_entries={len(rows)} output={OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
