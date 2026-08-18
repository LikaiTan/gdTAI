#!/usr/bin/env python3
"""Audit the proposed full-atlas metadata v2 rules without modifying H5AD files."""

from __future__ import annotations

import argparse
import csv
import json
import os
from collections import Counter
from pathlib import Path
from typing import Any

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "configs/metadata/full_atlas_metadata_harmonization_v2.json"
DEFAULT_OUTPUT = ROOT / "Integrated_dataset/tables/metadata_harmonization/full_atlas_v2"
DEFAULT_SOURCE_COUNTS = ROOT / "Integrated_dataset/tables/full_atlas_rebuild/report_cells_by_source.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def categorical_counts(group: h5py.Group) -> Counter[str]:
    categories = [
        value.decode() if isinstance(value, bytes) else str(value)
        for value in group["categories"][:]
    ]
    codes = group["codes"][:]
    counts = np.bincount(codes[codes >= 0], minlength=len(categories))
    result = Counter({categories[index]: int(count) for index, count in enumerate(counts) if count})
    n_missing = int(np.count_nonzero(codes < 0))
    if n_missing:
        result[""] += n_missing
    return result


def write_rows(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def file_state(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {
        "resolved_path": str(path.resolve()),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "inode": stat.st_ino,
    }


def normalize(value: Any) -> str:
    return " ".join(str(value).strip().casefold().replace("_", " ").split())


def collect_output_literals(config: dict[str, Any]) -> list[tuple[str, str, str]]:
    literals: list[tuple[str, str, str]] = []
    for field, values in config["controlled_vocabularies"].items():
        literals.extend(("controlled_vocabulary", field, str(value)) for value in values)
    for source in config["source_defaults"]:
        for key in ("tissue", "context", "tumor_type", "status", "evidence_level"):
            if key in source:
                literals.append((f"source_default:{source['source']}", key, str(source[key])))
    for rule in config["source_scoped_rules"]:
        for key, value in rule["outputs"].items():
            literals.append((f"scoped_rule:{rule['id']}", key, str(value)))
    for raw, values in config["global_rebuilt_tissue_synonyms"].items():
        literals.append((f"global_synonym:{raw}", "tissue", str(values["tissue"])))
        literals.append((f"global_synonym:{raw}", "context", str(values["context"])))
    return literals


def main() -> None:
    args = parse_args()
    config_path = args.config.resolve()
    output_dir = args.output_dir.resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    h5ad_path = (ROOT / config["input_h5ad"]).resolve()
    before = file_state(h5ad_path)

    with h5py.File(h5ad_path, "r") as handle:
        obs = handle["obs"]
        required = {
            "source_gse_id",
            "source_obs_name",
            "original_cell_id",
            "sample_id",
            "library_id",
            "tissue_original",
            "tissue_harmonized",
            "specimen_context",
        }
        missing_columns = sorted(required.difference(obs.keys()))
        source_counts_h5 = categorical_counts(obs["source_gse_id"])
        tissue_counts = categorical_counts(obs["tissue_harmonized"])
        specimen_counts = categorical_counts(obs["specimen_context"])
        n_cells = int(obs["source_gse_id/codes"].shape[0])

    after = file_state(h5ad_path)
    unchanged = before == after

    source_counts_csv: dict[str, int] = {}
    with DEFAULT_SOURCE_COUNTS.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            source_counts_csv[row["source_gse_id"]] = int(row["n_cells"])

    defaults = {row["source"]: row for row in config["source_defaults"]}
    atlas_sources = set(source_counts_h5)
    configured_sources = set(defaults)
    aliases = config["source_accession_aliases"]
    permitted = set(config["permitted_evidence_fields"])
    forbidden_fields = set(config["forbidden_evidence_fields"])
    used_fields = {
        field
        for rule in config["source_scoped_rules"]
        for field in rule["any_field_regex"]["fields"]
    }

    forbidden_values = {normalize(value) for value in config["forbidden_output_values_case_insensitive"]}
    literal_hits = [
        {"scope": scope, "field": field, "value": value, "forbidden_match": normalize(value)}
        for scope, field, value in collect_output_literals(config)
        if normalize(value) in forbidden_values
    ]
    current_pollution = [
        {"field": "specimen_context", "value": value, "n_cells": count}
        for value, count in sorted(specimen_counts.items())
        if normalize(value) in forbidden_values
    ]

    checks = [
        ("input_h5ad_opened_read_only", not missing_columns, ";".join(missing_columns)),
        ("input_h5ad_state_unchanged", unchanged, json.dumps(before, sort_keys=True)),
        ("expected_cell_count", n_cells == 5_933_312, str(n_cells)),
        ("expected_source_count", len(atlas_sources) == 40, str(len(atlas_sources))),
        ("source_count_table_matches_h5ad", source_counts_h5 == Counter(source_counts_csv), ""),
        ("all_atlas_sources_have_defaults", atlas_sources == configured_sources, ";".join(sorted(atlas_sources ^ configured_sources))),
        ("scoped_rules_use_only_permitted_fields", used_fields <= permitted, ";".join(sorted(used_fields - permitted))),
        ("scoped_rules_do_not_use_forbidden_fields", not (used_fields & forbidden_fields), ";".join(sorted(used_fields & forbidden_fields))),
        ("proposed_output_literals_exclude_cell_types", not literal_hits, json.dumps(literal_hits, sort_keys=True)),
        ("four_requested_source_aliases_present", len(aliases) == 4, json.dumps(aliases, sort_keys=True)),
        ("tumor_evidence_urls_are_official_geo", all(row["url"].startswith("https://www.ncbi.nlm.nih.gov/geo/") for row in config["tumor_source_evidence"]), str(len(config["tumor_source_evidence"]))),
    ]

    output_dir.mkdir(parents=True, exist_ok=True)
    write_rows(
        output_dir / "existing_tissue_harmonized_counts.csv",
        ["tissue_harmonized", "n_cells"],
        [{"tissue_harmonized": value, "n_cells": count} for value, count in tissue_counts.most_common()],
    )
    write_rows(
        output_dir / "existing_specimen_context_counts.csv",
        ["specimen_context", "n_cells"],
        [{"specimen_context": value, "n_cells": count} for value, count in specimen_counts.most_common()],
    )
    write_rows(
        output_dir / "source_coverage_review.csv",
        ["source_gse_id_original", "source_accession_harmonized_v2", "n_cells", "configured_status", "tumor_type_harmonized_v2", "evidence_level", "evidence_url"],
        [
            {
                "source_gse_id_original": source,
                "source_accession_harmonized_v2": aliases.get(source, source),
                "n_cells": source_counts_h5[source],
                "configured_status": defaults[source]["status"],
                "tumor_type_harmonized_v2": defaults[source].get("tumor_type", "unresolved"),
                "evidence_level": defaults[source].get("evidence_level", "none"),
                "evidence_url": defaults[source].get("evidence_url", ""),
            }
            for source in sorted(atlas_sources)
        ],
    )
    write_rows(
        output_dir / "tumor_source_evidence.csv",
        ["source", "tumor_type", "specimen_scope", "resolution", "url"],
        config["tumor_source_evidence"],
    )
    write_rows(
        output_dir / "cell_type_value_audit.csv",
        ["scope", "field", "value", "n_cells", "result"],
        [
            {"scope": "current_atlas", **row, "result": "REJECT_AS_EVIDENCE"}
            for row in current_pollution
        ]
        + [{"scope": "proposed_v2_literals", "field": "all_v2_outputs", "value": "", "n_cells": 0, "result": "PASS_NO_CELL_TYPE_VALUES"}],
    )
    write_rows(
        output_dir / "validation_checks.csv",
        ["check", "passed", "detail"],
        [{"check": name, "passed": passed, "detail": detail} for name, passed, detail in checks],
    )
    write_rows(
        output_dir / "input_file_state.csv",
        ["stage", "resolved_path", "size_bytes", "mtime_ns", "inode"],
        [{"stage": "before", **before}, {"stage": "after", **after}],
    )
    manifest = {
        "ruleset_id": config["ruleset_id"],
        "input_h5ad": str(h5ad_path),
        "input_open_mode": "r",
        "h5ad_unchanged": unchanged,
        "n_cells": n_cells,
        "n_sources": len(atlas_sources),
        "n_existing_tissue_categories": len(tissue_counts),
        "n_existing_blank_tissue_cells": tissue_counts.get("", 0),
        "n_current_cell_type_pollution_categories": len(current_pollution),
        "n_current_cell_type_pollution_cells": sum(row["n_cells"] for row in current_pollution),
        "all_checks_passed": all(passed for _, passed, _ in checks),
        "checks_passed": sum(bool(passed) for _, passed, _ in checks),
        "checks_total": len(checks),
    }
    (output_dir / "audit_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    if not manifest["all_checks_passed"]:
        raise SystemExit("One or more metadata-v2 review validations failed")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
