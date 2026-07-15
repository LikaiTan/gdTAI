#!/usr/bin/env python3
"""Build deterministic dataset, library, and file registries from legacy inputs."""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from tnk_atlas.paths import ProjectPaths  # noqa: E402
from tnk_atlas.registry import DATASET_COLUMNS, FILE_COLUMNS, LIBRARY_COLUMNS  # noqa: E402


GSE_PATTERN = re.compile(r"^GSE\d+$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-registry", default="configs/datasets/integration_inputs.csv")
    parser.add_argument("--output-dir", default="configs/datasets")
    return parser.parse_args()


def relative_or_absolute(path: Path) -> str:
    try:
        return path.resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()
    except (ValueError, FileNotFoundError):
        return str(path)


def read_integration_registry(path: Path) -> dict[str, Path]:
    selected: dict[str, Path] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            dataset_id = row["gse_id"].strip()
            selected[dataset_id] = Path(row["h5ad_path"]).expanduser()
    return selected


def read_current_counts(path: Path) -> dict[str, int]:
    if not path.exists():
        return {}
    counts: dict[str, int] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            dataset_id = (row.get("gse_id") or row.get("source_gse_id") or "").strip()
            raw_count = row.get("n_cells") or row.get("cell_n") or "0"
            if dataset_id:
                counts[dataset_id] = int(float(raw_count))
    return counts


def read_extended_atlas_counts(path: Path) -> dict[str, int]:
    if not path.exists():
        return {}
    counts: dict[str, int] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            dataset_id = row.get("source_gse_id", "").strip()
            if dataset_id and dataset_id not in counts:
                counts[dataset_id] = int(float(row["total_cells"]))
    return counts


def discover_dataset_ids(
    paths: ProjectPaths,
    selected: dict[str, Path],
    extended_counts: dict[str, int],
) -> set[str]:
    dataset_ids = set(selected)
    dataset_ids.update(extended_counts)
    if paths.legacy_downloads.exists():
        dataset_ids.update(
            child.name
            for child in paths.legacy_downloads.iterdir()
            if child.is_dir() and GSE_PATTERN.match(child.name)
        )
    if paths.legacy_scanpy_projects.exists():
        dataset_ids.update(
            child.name
            for child in paths.legacy_scanpy_projects.iterdir()
            if child.is_dir() and GSE_PATTERN.match(child.name)
        )
    dataset_ids.update({"HRA005041", "GDT_2020AUG_woCOV", "GDTlung2023july_7p", "MalteGDT"})
    return dataset_ids


def discover_processed_candidates(paths: ProjectPaths) -> dict[str, list[Path]]:
    candidates: dict[str, list[Path]] = defaultdict(list)
    per_gse = paths.legacy_downloads / "per_gse_h5ad_with_metadata"
    if per_gse.exists():
        for path in sorted(per_gse.glob("*.h5ad")):
            if path.name.startswith(".") or ".bak_" in path.name:
                continue
            match = re.match(r"(GSE\d+|HRA\d+)", path.name)
            if match:
                candidates[match.group(1)].append(path)
    if paths.legacy_scanpy_projects.exists():
        for path in sorted(paths.legacy_scanpy_projects.glob("GSE*/outputs/*.h5ad")):
            if ".bak" not in path.name:
                candidates[path.parents[1].name].append(path)
    sorted_root = paths.root / "newdata" / "Sorted_gdT"
    for dataset_id in ("GDT_2020AUG_woCOV", "GDTlung2023july_7p", "MalteGDT"):
        path = sorted_root / f"{dataset_id}_sorted_gdt.h5ad"
        if path.exists():
            candidates[dataset_id].append(path)
    return candidates


def choose_processed(
    dataset_id: str,
    selected: dict[str, Path],
    candidates: dict[str, list[Path]],
) -> Path | None:
    if dataset_id in selected:
        return selected[dataset_id]
    options = candidates.get(dataset_id, [])
    return options[0] if options else None


def raw_path(paths: ProjectPaths, dataset_id: str) -> Path | None:
    if GSE_PATTERN.match(dataset_id):
        candidate = paths.legacy_downloads / dataset_id
        return candidate if candidate.exists() else None
    if dataset_id == "HRA005041":
        return paths.root / "newdata" / "HRA005041"
    if dataset_id in {"GDT_2020AUG_woCOV", "GDTlung2023july_7p", "MalteGDT"}:
        return paths.root / "newdata" / "Sorted_gdT"
    return None


def dataset_rows(
    paths: ProjectPaths,
    selected: dict[str, Path],
    current_counts: dict[str, int],
    extended_counts: dict[str, int],
    candidates: dict[str, list[Path]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for dataset_id in sorted(discover_dataset_ids(paths, selected, extended_counts)):
        processed = choose_processed(dataset_id, selected, candidates)
        phase0_active = dataset_id in selected
        current_active = dataset_id in current_counts and current_counts[dataset_id] > 0
        extended_active = dataset_id in extended_counts and extended_counts[dataset_id] > 0
        if current_active:
            role = "current_atlas_input"
            status = "active"
        elif extended_active:
            role = "extended_five_million_atlas_input"
            status = "active_extended_atlas"
        elif phase0_active:
            role = "phase0_input_not_in_current_milestone"
            status = "excluded_after_phase0"
        elif dataset_id.startswith("GDT") or dataset_id == "MalteGDT":
            role = "sorted_gdt_reference"
            status = "reference"
        elif dataset_id == "HRA005041":
            role = "external_tcr_reference"
            status = "reference"
        else:
            role = "available_not_selected"
            status = "available"
        legacy_raw = raw_path(paths, dataset_id)
        rows.append(
            {
                "dataset_id": dataset_id,
                "source_type": "GEO" if GSE_PATTERN.match(dataset_id) else "local_external",
                "accession": dataset_id,
                "phase0_active": str(phase0_active).lower(),
                "current_milestone_active": str(current_active).lower(),
                "current_cell_count": current_counts.get(dataset_id, 0),
                "extended_atlas_active": str(extended_active).lower(),
                "extended_atlas_cell_count": extended_counts.get(dataset_id, 0),
                "integration_role": role,
                "legacy_raw_path": relative_or_absolute(legacy_raw) if legacy_raw else "",
                "processed_h5ad_path": relative_or_absolute(processed) if processed else "",
                "status": status,
                "exclusion_reason": "" if current_active or not phase0_active else "not present in current milestone",
                "notes": "Generated from legacy paths; assay and cohort semantics require curated review.",
            }
        )
    return rows


def infer_file_role(path: Path) -> tuple[str, str]:
    lower = path.name.lower()
    parts = {part.lower() for part in path.parts}
    if path.suffix.lower() == ".h5ad":
        return "processed", "standardized_h5ad"
    if "matrix" in parts or "series_matrix" in lower:
        return "raw", "geo_series_matrix"
    if lower.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz", ".sra")):
        return "raw", "sequence_reads"
    if "contig" in lower or "tcr" in lower or "vdj" in lower:
        return "raw", "tcr_source"
    if "metadata" in lower or "annotation" in lower:
        return "raw", "source_metadata"
    if lower.endswith((".tar", ".tar.gz", ".tgz", ".zip", ".rds", ".rds.gz", ".rdata", ".rdata.gz")):
        return "raw", "supplementary_archive"
    if path.suffix.lower() in {".py", ".r", ".sh"}:
        return "archive", "legacy_dataset_code"
    if path.suffix.lower() in {".log", ".png", ".pdf", ".html", ".h5", ".loom"}:
        return "interim", "legacy_derived_or_analysis"
    return "raw", "unclassified_source"


def bounded_source_files(root: Path):
    if not root.exists():
        return
    for path in sorted(root.iterdir()):
        if path.is_file() or path.is_symlink():
            yield path
        elif path.is_dir():
            for child in sorted(path.iterdir()):
                if child.is_file() or child.is_symlink():
                    yield child


def canonical_file_path(dataset_id: str, level: str, role: str, path: Path) -> str:
    if level == "raw":
        return f"data/raw/geo/{dataset_id}/{role}/{path.name}"
    if level == "processed":
        return f"data/processed/per_dataset/{dataset_id}/{path.name}"
    if level == "archive":
        return f"archive/legacy_pipeline/dataset_embedded/{dataset_id}/{path.name}"
    return f"data/interim/{dataset_id}/{role}/{path.name}"


def file_rows(
    paths: ProjectPaths,
    datasets: list[dict[str, Any]],
    candidates: dict[str, list[Path]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for dataset in datasets:
        dataset_id = str(dataset["dataset_id"])
        raw = raw_path(paths, dataset_id)
        sources = list(bounded_source_files(raw)) if raw and raw.exists() else []
        sources.extend(candidates.get(dataset_id, []))
        selected_path = str(dataset["processed_h5ad_path"])
        if selected_path:
            candidate = Path(selected_path)
            if not candidate.is_absolute():
                candidate = PROJECT_ROOT / candidate
            sources.append(candidate)
        for path in sources:
            key = (dataset_id, str(path.resolve()))
            if key in seen or not path.exists():
                continue
            seen.add(key)
            level, role = infer_file_role(path)
            stat = path.stat()
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "lifecycle_level": level,
                    "file_role": role,
                    "legacy_path": relative_or_absolute(path),
                    "canonical_path": canonical_file_path(dataset_id, level, role, path),
                    "size_bytes": stat.st_size,
                    "mtime_ns": stat.st_mtime_ns,
                    "checksum_type": "deferred_until_dataset_migration",
                    "checksum": "",
                }
            )
    return sorted(rows, key=lambda row: (row["dataset_id"], row["lifecycle_level"], row["legacy_path"]))


def library_rows(
    datasets: list[dict[str, Any]],
    sample_information: Path,
) -> list[dict[str, str]]:
    active_by_dataset = {
        str(row["dataset_id"]): str(row["phase0_active"]) for row in datasets
    }
    rows: list[dict[str, str]] = []
    represented: set[str] = set()
    seen: set[tuple[str, str]] = set()
    if sample_information.exists():
        with sample_information.open(newline="", encoding="utf-8-sig") as handle:
            for source in csv.DictReader(handle):
                dataset_id = (source.get("gse") or source.get("project_id") or "").strip()
                if dataset_id not in active_by_dataset:
                    continue
                library_id = (source.get("library_id") or "").strip()
                if not library_id:
                    library_id = f"__missing_library_id_{source.get('library_index', '')}__"
                key = (dataset_id, library_id)
                if key in seen:
                    continue
                seen.add(key)
                represented.add(dataset_id)
                tcr_flag = (source.get("tcr_vdj_flag") or "").strip()
                lower_tcr = tcr_flag.lower()
                if "without" in lower_tcr or lower_tcr in {"none", "no", "false"}:
                    tcr_assay = "none"
                    tcr_scope = "none"
                elif tcr_flag:
                    tcr_assay = "vdj"
                    tcr_scope = "unspecified_ab_or_gd"
                else:
                    tcr_assay = "unknown"
                    tcr_scope = "unknown"
                notes = "; ".join(
                    value
                    for value in (
                        (source.get("sample_type") or "").strip(),
                        (source.get("donor_patient") or "").strip(),
                    )
                    if value
                )
                rows.append(
                    {
                        "dataset_id": dataset_id,
                        "library_id": library_id,
                        "sample_id": library_id,
                        "rna_assay": (source.get("technology_simple") or "unknown").strip()
                        or "unknown",
                        "tcr_assay": tcr_assay,
                        "tcr_scope": tcr_scope,
                        "active": active_by_dataset[dataset_id],
                        "notes": notes,
                    }
                )
    for dataset in datasets:
        dataset_id = str(dataset["dataset_id"])
        if dataset_id in represented:
            continue
        processed_value = str(dataset["processed_h5ad_path"])
        processed_path = Path(processed_value) if processed_value else None
        if processed_path is not None and not processed_path.is_absolute():
            processed_path = PROJECT_ROOT / processed_path
        derived_samples: list[str] = []
        derived_tcr: list[str] = []
        if processed_path is not None and processed_path.exists():
            derived_samples = h5ad_first_unique_column(
                processed_path,
                ("library_id", "sample_id", "sampleID", "sampleid", "phase1_library_label"),
            )
            derived_tcr = h5ad_first_unique_column(processed_path, ("TCRseq", "tcr_chain_mode"))
        if str(dataset["phase0_active"]).lower() == "true":
            samples = derived_samples or ["__dataset_level_no_sample_id__"]
            tcr_text = ";".join(derived_tcr)
            tcr_lower = tcr_text.lower()
            if tcr_text and not all(token in {"false", "none", "no", "0", "nan"} for token in tcr_lower.split(";")):
                tcr_assay = "vdj"
                tcr_scope = tcr_text
            else:
                tcr_assay = "none" if tcr_text else "unknown"
                tcr_scope = "none" if tcr_text else "unknown"
            for sample_id in samples:
                rows.append(
                    {
                        "dataset_id": dataset_id,
                        "library_id": sample_id,
                        "sample_id": sample_id,
                        "rna_assay": "unknown",
                        "tcr_assay": tcr_assay,
                        "tcr_scope": tcr_scope,
                        "active": "true",
                        "notes": "Derived from selected H5AD obs because the sample-information manifest had no row.",
                    }
                )
            represented.add(dataset_id)
            continue
        rows.append(
            {
                "dataset_id": dataset_id,
                "library_id": "__dataset_level_pending_curation__",
                "sample_id": "",
                "rna_assay": "unknown",
                "tcr_assay": "unknown",
                "tcr_scope": "unknown",
                "active": str(dataset["phase0_active"]),
                "notes": "No row in sample_information_final_full.csv; curate before activation.",
            }
        )
    return sorted(rows, key=lambda row: (row["dataset_id"], row["library_id"]))


def h5ad_first_unique_column(path: Path, candidates: tuple[str, ...]) -> list[str]:
    import h5py
    import numpy as np

    def decode(value: Any) -> str:
        if isinstance(value, bytes):
            return value.decode("utf-8", errors="replace")
        return str(value)

    def string_values(node: Any) -> list[str]:
        if isinstance(node, h5py.Dataset):
            return [decode(value) for value in node[:]]
        if isinstance(node, h5py.Group) and "values" in node:
            values = [decode(value) for value in node["values"][:]]
            if "mask" in node:
                mask = node["mask"][:].astype(bool)
                return [value for value, missing in zip(values, mask) if not missing]
            return values
        return []

    with h5py.File(path, "r") as handle:
        obs = handle.get("obs")
        if obs is None:
            return []
        for column in candidates:
            if column not in obs:
                continue
            node = obs[column]
            values: set[str] = set()
            if isinstance(node, h5py.Group) and "categories" in node:
                categories = string_values(node["categories"])
                if "codes" in node:
                    used_codes = np.unique(node["codes"][:])
                    values.update(
                        categories[int(code)]
                        for code in used_codes
                        if 0 <= int(code) < len(categories)
                    )
                else:
                    values.update(categories)
            elif isinstance(node, h5py.Dataset):
                chunk_size = 250_000
                for start in range(0, node.shape[0], chunk_size):
                    chunk = node[start : start + chunk_size]
                    values.update(decode(value) for value in np.unique(chunk))
            cleaned = sorted(value for value in values if value.strip() and value.lower() != "nan")
            if cleaned:
                return cleaned
    return []


def write_csv(path: Path, columns: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def main() -> int:
    args = parse_args()
    paths = ProjectPaths.discover(PROJECT_ROOT)
    input_registry = PROJECT_ROOT / args.input_registry
    output_dir = PROJECT_ROOT / args.output_dir
    selected = read_integration_registry(input_registry)
    current_counts = read_current_counts(paths.tables / "current_integrated_cells_by_gse.csv")
    extended_counts = read_extended_atlas_counts(
        paths.tables
        / "gdT_prediction"
        / "gdtai_v3_trdc_nk_guard"
        / "full_atlas_prediction_by_source.csv"
    )
    candidates = discover_processed_candidates(paths)
    datasets = dataset_rows(paths, selected, current_counts, extended_counts, candidates)
    files = file_rows(paths, datasets, candidates)
    libraries = library_rows(
        datasets,
        PROJECT_ROOT / "configs" / "datasets" / "sample_information_final_full.csv",
    )
    write_csv(output_dir / "datasets.csv", DATASET_COLUMNS, datasets)
    write_csv(output_dir / "libraries.csv", LIBRARY_COLUMNS, libraries)
    write_csv(output_dir / "files.csv", FILE_COLUMNS, files)
    print(f"datasets={len(datasets)} libraries={len(libraries)} files={len(files)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
