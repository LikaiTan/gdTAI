"""Dataset-registry schema and validation."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


DATASET_COLUMNS = (
    "dataset_id",
    "source_type",
    "accession",
    "phase0_active",
    "current_milestone_active",
    "current_cell_count",
    "extended_atlas_active",
    "extended_atlas_cell_count",
    "integration_role",
    "legacy_raw_path",
    "processed_h5ad_path",
    "status",
    "exclusion_reason",
    "notes",
)

LIBRARY_COLUMNS = (
    "dataset_id",
    "library_id",
    "sample_id",
    "rna_assay",
    "tcr_assay",
    "tcr_scope",
    "active",
    "notes",
)

FILE_COLUMNS = (
    "dataset_id",
    "lifecycle_level",
    "file_role",
    "legacy_path",
    "canonical_path",
    "size_bytes",
    "mtime_ns",
    "checksum_type",
    "checksum",
)


@dataclass(frozen=True)
class RegistryValidation:
    dataset_count: int
    active_dataset_count: int
    missing_processed_inputs: tuple[str, ...]
    duplicate_dataset_ids: tuple[str, ...]

    @property
    def ok(self) -> bool:
        return not self.missing_processed_inputs and not self.duplicate_dataset_ids


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def validate_dataset_registry(path: Path, project_root: Path) -> RegistryValidation:
    rows = read_csv(path)
    seen: set[str] = set()
    duplicates: set[str] = set()
    missing: list[str] = []
    active_count = 0
    for row in rows:
        dataset_id = row["dataset_id"]
        if dataset_id in seen:
            duplicates.add(dataset_id)
        seen.add(dataset_id)
        active = row["phase0_active"].strip().lower() == "true"
        if active:
            active_count += 1
            value = row["processed_h5ad_path"].strip()
            if not value:
                missing.append(f"{dataset_id}:<empty>")
            else:
                candidate = Path(value)
                if not candidate.is_absolute():
                    candidate = project_root / candidate
                if not candidate.exists():
                    missing.append(f"{dataset_id}:{value}")
    return RegistryValidation(len(rows), active_count, tuple(missing), tuple(sorted(duplicates)))
