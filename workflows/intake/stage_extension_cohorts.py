#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Stage immutable extension-cohort sources from the extension manifests."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import sys
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence

from validate_extension_cohorts import (
    COHORT_COLUMNS,
    DEFAULT_CANONICAL_REGISTRY,
    DEFAULT_COHORTS,
    DEFAULT_LIBRARIES,
    LIBRARY_COLUMNS,
    ManifestError,
    assert_extension_output_root,
    cohort_accessions,
    parse_bool,
    split_tokens,
    validate_extension_manifests,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE_ROOT = PROJECT_ROOT / "data/compat/downloads/new_candidata_data"
DEFAULT_SHARED_METASHEET = DEFAULT_SOURCE_ROOT / "sample_metasheet.csv"
DEFAULT_EXTENSION_ROOT = PROJECT_ROOT / "data/interim/extension_intake"
DEFAULT_STAGE_ROOT = DEFAULT_EXTENSION_ROOT / "staged"
SHARED_METASHEET_REQUIRED_COLUMNS = ("gse", "gsm", "sample_title", "modality")
MATRIX_COMPANION_PATTERNS = (
    "barcodes.tsv",
    "barcodes.tsv.gz",
    "features.tsv",
    "features.tsv.gz",
    "genes.tsv",
    "genes.tsv.gz",
)


@dataclass(frozen=True)
class SourcePlan:
    cohort_id: str
    library_id: str
    source_role: str
    required: bool
    matches: tuple[Path, ...]


@dataclass(frozen=True)
class StagedSource:
    cohort_id: str
    library_id: str
    source_role: str
    source_path: str
    staged_path: str
    size_bytes: int
    mtime_ns: int
    sha256: str


def _is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def _matches_for_globs(source_root: Path, raw_patterns: str, stage_root: Path) -> tuple[Path, ...]:
    matches: set[Path] = set()
    stage_root = stage_root.resolve()
    for pattern in split_tokens(raw_patterns):
        for candidate in source_root.glob(pattern):
            if not candidate.is_file():
                continue
            resolved = candidate.resolve()
            if resolved == stage_root or _is_relative_to(resolved, stage_root):
                continue
            matches.add(resolved)
    return tuple(sorted(matches, key=str))


def build_stage_plan(
    source_root: Path,
    stage_root: Path,
    cohorts: Sequence[dict[str, str]],
    libraries: Sequence[dict[str, str]],
) -> list[SourcePlan]:
    """Resolve manifest source globs without changing the filesystem."""
    source_root = source_root.expanduser().resolve()
    stage_root = assert_extension_output_root(stage_root)
    enabled = {row["cohort_id"] for row in cohorts if parse_bool(row["stage_enabled"])}
    plans: list[SourcePlan] = []
    for row in libraries:
        if row["cohort_id"] not in enabled:
            continue
        matches = _matches_for_globs(source_root, row["source_glob"], stage_root)
        plans.append(
            SourcePlan(
                cohort_id=row["cohort_id"],
                library_id=row["library_id"],
                source_role=row["source_role"],
                required=parse_bool(row["required"]),
                matches=matches,
            )
        )
    return plans


def filter_shared_metasheet(
    path: Path, cohort: dict[str, str]
) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    """Filter the shared source table to the cohort's linked accessions."""
    if not path.is_file():
        raise FileNotFoundError(f"Shared source metasheet is missing: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = tuple(reader.fieldnames or ())
        missing_columns = set(SHARED_METASHEET_REQUIRED_COLUMNS) - set(fieldnames)
        if missing_columns:
            raise ManifestError(
                f"Shared source metasheet lacks columns {sorted(missing_columns)}: {path}"
            )
        accessions = cohort_accessions(cohort)
        rows = [
            {key: (value or "").strip() for key, value in row.items()}
            for row in reader
            if (row.get("gse") or "").strip().upper() in accessions
        ]
    rows.sort(key=lambda row: (row["gse"], row["modality"], row["gsm"], row["sample_title"]))
    present = {row["gse"].upper() for row in rows}
    expected = {
        accession
        for accession in accessions
        if accession.startswith("GSE")
    }
    missing_accessions = expected - present
    if missing_accessions:
        raise ManifestError(
            f"Shared source metasheet has no rows for {cohort['cohort_id']} accession(s): "
            f"{sorted(missing_accessions)}"
        )
    if not rows:
        raise ManifestError(f"Shared source metasheet has no rows for {cohort['cohort_id']}")
    duplicate_gsm = len({row["gsm"] for row in rows}) != len(rows)
    if duplicate_gsm:
        raise ManifestError(f"Shared source metasheet has duplicate GSM rows for {cohort['cohort_id']}")
    return fieldnames, rows


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _matrix_companions(path: Path) -> list[Path]:
    lowered = path.name.lower()
    if "matrix.mtx" not in lowered and "counts.mtx" not in lowered:
        return []
    companions = []
    for name in MATRIX_COMPANION_PATTERNS:
        candidate = path.parent / name
        if candidate.is_file():
            companions.append(candidate.resolve())
    if not companions:
        prefix = path.name.split("matrix.mtx", 1)[0].split("counts.mtx", 1)[0]
        for candidate in path.parent.iterdir():
            name = candidate.name.lower()
            if candidate.is_file() and candidate.name.startswith(prefix) and any(
                token in name for token in ("barcode", "feature", "gene")
            ):
                companions.append(candidate.resolve())
    return sorted(set(companions), key=str)


def _link_source(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"Refusing to replace staged source: {destination}")
    os.symlink(source, destination)


def _snapshot_csv(path: Path, columns: Sequence[str], rows: Iterable[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _stage_one_cohort(
    stage_root: Path,
    source_root: Path,
    cohort: dict[str, str],
    cohort_libraries: list[dict[str, str]],
    plans: list[SourcePlan],
    shared_metasheet: Path,
    shared_columns: Sequence[str],
    shared_rows: list[dict[str, str]],
) -> Path:
    cohort_id = cohort["cohort_id"]
    target = stage_root / cohort_id
    if target.exists() or target.is_symlink():
        raise FileExistsError(f"Refusing to replace existing stage directory: {target}")

    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{cohort_id}.", dir=target.parent))
    staged_sources: list[StagedSource] = []
    try:
        for plan in plans:
            files = set(plan.matches)
            for match in plan.matches:
                files.update(_matrix_companions(match))
            used_names: set[str] = set()
            for source in sorted(files, key=str):
                if source.name in used_names:
                    raise ManifestError(
                        f"Duplicate basename for {plan.library_id}: {source.name}; split or rename the source bundle"
                    )
                used_names.add(source.name)
                relative = Path("sources") / plan.library_id / source.name
                destination = temporary / relative
                _link_source(source, destination)
                stat = source.stat()
                staged_sources.append(
                    StagedSource(
                        cohort_id=cohort_id,
                        library_id=plan.library_id,
                        source_role=plan.source_role,
                        source_path=str(source),
                        staged_path=str(relative),
                        size_bytes=stat.st_size,
                        mtime_ns=stat.st_mtime_ns,
                        sha256=sha256_file(source),
                    )
                )

        filtered_metasheet = temporary / "sample_metasheet.csv"
        _snapshot_csv(filtered_metasheet, shared_columns, shared_rows)
        metasheet_stat = shared_metasheet.stat()
        metasheet_record: dict[str, object] = {
            "source_path": str(shared_metasheet),
            "source_size_bytes": metasheet_stat.st_size,
            "source_mtime_ns": metasheet_stat.st_mtime_ns,
            "source_sha256": sha256_file(shared_metasheet),
            "staged_path": "sample_metasheet.csv",
            "staged_sha256": sha256_file(filtered_metasheet),
            "row_count": len(shared_rows),
            "accessions": sorted({row["gse"] for row in shared_rows}),
        }

        _snapshot_csv(temporary / "cohort_manifest.csv", COHORT_COLUMNS, [cohort])
        _snapshot_csv(temporary / "library_manifest.csv", LIBRARY_COLUMNS, cohort_libraries)
        manifest = {
            "schema_version": 1,
            "cohort_id": cohort_id,
            "source_root": str(source_root),
            "join_key": "sample_id+barcode_core",
            "immutable_sources": True,
            "count_provenance": cohort["count_provenance"],
            "canonical_tcr_flags": split_tokens(cohort["canonical_tcr_flags"]),
            "metasheet": metasheet_record,
            "sources": [asdict(record) for record in staged_sources],
        }
        (temporary / "stage_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        temporary.rename(target)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return target


def stage_extension_cohorts(
    source_root: Path,
    stage_root: Path,
    cohorts: Sequence[dict[str, str]],
    libraries: Sequence[dict[str, str]],
    shared_metasheet: Path = DEFAULT_SHARED_METASHEET,
    dry_run: bool = False,
) -> dict[str, object]:
    """Plan or materialize immutable extension-cohort source stages."""
    source_root = source_root.expanduser().resolve()
    stage_root = assert_extension_output_root(stage_root)
    shared_metasheet = shared_metasheet.expanduser().resolve()
    if not source_root.is_dir():
        if not dry_run:
            raise FileNotFoundError(f"Source root does not exist: {source_root}")

    enabled_cohorts = [row for row in cohorts if parse_bool(row["stage_enabled"])]
    blocked = [row["cohort_id"] for row in cohorts if row["intake_role"] == "blocked_provenance"]
    if blocked and len(cohorts) == len(blocked):
        raise ManifestError(
            f"Provenance-blocked cohort(s) cannot be staged: {', '.join(sorted(blocked))}"
        )
    filtered_metasheets: dict[str, tuple[tuple[str, ...], list[dict[str, str]]]] = {
        cohort["cohort_id"]: filter_shared_metasheet(shared_metasheet, cohort)
        for cohort in enabled_cohorts
    }
    plans = build_stage_plan(source_root, stage_root, enabled_cohorts, libraries)
    missing = [plan for plan in plans if plan.required and not plan.matches]
    summary: dict[str, object] = {
        "dry_run": dry_run,
        "source_root": str(source_root),
        "stage_root": str(stage_root),
        "shared_metasheet": str(shared_metasheet),
        "cohorts": [row["cohort_id"] for row in enabled_cohorts],
        "metasheet_rows": {
            cohort_id: len(rows) for cohort_id, (_, rows) in filtered_metasheets.items()
        },
        "required_sources_missing": [
            {"cohort_id": plan.cohort_id, "library_id": plan.library_id} for plan in missing
        ],
        "source_matches": {
            plan.library_id: [str(path) for path in plan.matches] for plan in plans
        },
        "staged": [],
    }
    if dry_run:
        return summary
    if missing:
        details = ", ".join(f"{plan.cohort_id}/{plan.library_id}" for plan in missing)
        raise FileNotFoundError(f"Required extension sources are missing: {details}")

    libraries_by_cohort: dict[str, list[dict[str, str]]] = {}
    plans_by_cohort: dict[str, list[SourcePlan]] = {}
    for row in libraries:
        libraries_by_cohort.setdefault(row["cohort_id"], []).append(row)
    for plan in plans:
        plans_by_cohort.setdefault(plan.cohort_id, []).append(plan)

    staged: list[str] = []
    for cohort in enabled_cohorts:
        cohort_id = cohort["cohort_id"]
        shared_columns, shared_rows = filtered_metasheets[cohort_id]
        staged.append(
            str(
                _stage_one_cohort(
                    stage_root,
                    source_root,
                    cohort,
                    libraries_by_cohort.get(cohort_id, []),
                    plans_by_cohort.get(cohort_id, []),
                    shared_metasheet,
                    shared_columns,
                    shared_rows,
                )
            )
        )
    summary["staged"] = staged
    return summary


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument("--stage-root", type=Path, default=DEFAULT_STAGE_ROOT)
    parser.add_argument("--shared-metasheet", type=Path, default=DEFAULT_SHARED_METASHEET)
    parser.add_argument("--cohorts-manifest", type=Path, default=DEFAULT_COHORTS)
    parser.add_argument("--libraries-manifest", type=Path, default=DEFAULT_LIBRARIES)
    parser.add_argument("--canonical-registry", type=Path, default=DEFAULT_CANONICAL_REGISTRY)
    parser.add_argument("--cohort", action="append", default=[], help="Stage one cohort ID; repeatable.")
    parser.add_argument("--dry-run", action="store_true", help="Resolve and report inputs without writing.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    report, cohorts, libraries = validate_extension_manifests(
        args.cohorts_manifest,
        args.libraries_manifest,
        args.canonical_registry,
        args.cohort,
        args.shared_metasheet,
    )
    if not report.passed:
        print(json.dumps(report.to_dict(), indent=2, sort_keys=True), file=sys.stderr)
        return 1
    try:
        summary = stage_extension_cohorts(
            args.source_root,
            args.stage_root,
            cohorts,
            libraries,
            args.shared_metasheet,
            args.dry_run,
        )
    except (FileNotFoundError, FileExistsError, ManifestError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
