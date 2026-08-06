#!/usr/bin/env python3
"""Refine sample provenance with deterministic, review-gated H5AD writeback.

The ``review`` (alias ``dry-run``) command streams selected ``obs`` columns,
exports the proposed ``sample_source_refined`` values, and writes a manifest.
It never opens the H5AD for writing. The separate ``writeback`` command
requires that manifest and an explicit review acknowledgement before it writes
one new ``obs`` column.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import logging
import os
import re
import unicodedata
from collections import Counter
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Iterator, Mapping, Sequence, TextIO

import h5py
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_H5AD = REPO_ROOT / "high_speed_temp/Integrated_dataset/integrated.h5ad"
DEFAULT_RULES = REPO_ROOT / "configs/metadata/sample_source_refinement_rules.json"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "Integrated_dataset/tables/sample_source_refinement"
DEFAULT_CHUNK_SIZE = 100_000
EXPECTED_PRECEDENCE = [
    "cell_level",
    "sample_library",
    "gse_tissue",
    "global_tissue_fallback",
    "unresolved",
]
CANONICAL_FIELDS = [
    "obs_name",
    "source_gse_id",
    "cell_id",
    "cell_source",
    "sample_id",
    "library_id",
    "patient_id",
    "sample_context",
    "tissue",
]
REVIEW_COLUMNS = [
    "row_index",
    *CANONICAL_FIELDS,
    "sample_source_refined",
    "refinement_level",
    "refinement_rule_id",
    "refinement_reason",
]
STAGING_COLUMN = "__sample_source_refined_staging__"
BACKUP_COLUMN = "__sample_source_refined_previous__"


class RefinementError(RuntimeError):
    """Raised when rules, review provenance, or H5AD safety checks fail."""


@dataclass(frozen=True)
class Refinement:
    """One refined value and its deterministic provenance."""

    value: str
    level: str
    rule_id: str
    reason: str


@dataclass
class ScanResult:
    """Digests and summaries from one streaming H5AD pass."""

    n_obs: int
    source_sha256: str
    output_values_sha256: str
    resolutions_sha256: str
    source_columns: dict[str, list[str]]
    value_counts: Counter[tuple[str, str, str]] = field(default_factory=Counter)
    gse_totals: Counter[str] = field(default_factory=Counter)
    gse_unresolved: Counter[str] = field(default_factory=Counter)
    gse_values: dict[str, set[str]] = field(default_factory=dict)
    tumor_context_counts: Counter[tuple[str, str, bool]] = field(default_factory=Counter)
    unresolved_examples: list[list[object]] = field(default_factory=list)


def configure_logging(verbose: bool) -> None:
    """Configure concise command-line logging."""
    logging.basicConfig(
        level=logging.INFO if verbose else logging.WARNING,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )


def normalize_match_text(value: object) -> str:
    """Normalize punctuation and case for deterministic matching."""
    text = clean_scalar(value).casefold()
    text = unicodedata.normalize("NFKD", text)
    text = "".join(char for char in text if not unicodedata.combining(char))
    text = re.sub(r"[_/|+\-]+", " ", text)
    text = re.sub(r"[^\w\s]", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def clean_scalar(value: object) -> str:
    """Convert one scalar to stripped text without interpreting semantics."""
    if value is None:
        return ""
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and np.isnan(value):
        return ""
    return str(value).strip()


def slugify(value: str) -> str:
    """Create a stable ASCII label for a valid global fallback value."""
    text = unicodedata.normalize("NFKD", value)
    text = text.encode("ascii", errors="ignore").decode("ascii").casefold()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_")


def file_sha256(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    """Hash a file without loading it into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def _update_digest(digest: "hashlib._Hash", values: Sequence[object]) -> None:
    payload = json.dumps(
        [clean_scalar(value) for value in values],
        ensure_ascii=False,
        separators=(",", ":"),
    ).encode("utf-8")
    digest.update(len(payload).to_bytes(8, byteorder="big", signed=False))
    digest.update(payload)


def _atomic_replace(temp_path: Path, destination: Path) -> None:
    temp_path.replace(destination)


@contextmanager
def atomic_text_writer(path: Path, *, gzip_output: bool = False) -> Iterator[TextIO]:
    """Write text atomically; gzip streams use mtime=0 for reproducibility."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temp_path.open("wb") as raw:
            if gzip_output:
                with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                    with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                        yield text
            else:
                with io.TextIOWrapper(raw, encoding="utf-8", newline="") as text:
                    yield text
        _atomic_replace(temp_path, path)
    except Exception:
        temp_path.unlink(missing_ok=True)
        raise


def write_json(path: Path, payload: Mapping[str, object]) -> None:
    """Write stable, sorted JSON atomically."""
    with atomic_text_writer(path) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, ensure_ascii=True)
        handle.write("\n")


def write_csv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    """Write a deterministic CSV atomically."""
    with atomic_text_writer(path, gzip_output=path.suffix == ".gz") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _require_nonempty_string(mapping: Mapping[str, object], key: str, context: str) -> str:
    value = mapping.get(key)
    if not isinstance(value, str) or not value.strip():
        raise RefinementError(f"{context}.{key} must be a non-empty string")
    return value.strip()


def validate_rules(rules: Mapping[str, object]) -> None:
    """Validate the bounded rule schema and reject ambiguous exact maps."""
    if rules.get("schema_version") != 1:
        raise RefinementError("Only sample-source rule schema_version=1 is supported")
    if rules.get("output_field") != "sample_source_refined":
        raise RefinementError("output_field must be sample_source_refined")
    if rules.get("precedence") != EXPECTED_PRECEDENCE:
        raise RefinementError(f"precedence must be exactly {EXPECTED_PRECEDENCE}")

    source_columns = rules.get("source_columns")
    if not isinstance(source_columns, dict):
        raise RefinementError("source_columns must be an object")
    missing_canonical = set(CANONICAL_FIELDS) - set(source_columns)
    if missing_canonical:
        raise RefinementError(f"source_columns is missing: {sorted(missing_canonical)}")
    for canonical, aliases in source_columns.items():
        if not isinstance(aliases, list) or not aliases or not all(isinstance(item, str) for item in aliases):
            raise RefinementError(f"source_columns.{canonical} must be a non-empty string list")

    seen_ids: set[str] = set()

    def register_id(rule: Mapping[str, object], context: str) -> None:
        rule_id = _require_nonempty_string(rule, "id", context)
        if rule_id in seen_ids:
            raise RefinementError(f"Duplicate rule id: {rule_id}")
        seen_ids.add(rule_id)

    for index, rule in enumerate(rules.get("cell_level_overrides", [])):
        register_id(rule, f"cell_level_overrides[{index}]")
        _require_nonempty_string(rule, "gse", f"cell_level_overrides[{index}]")
        if not rule.get("cell_id_exact"):
            raise RefinementError(f"cell_level_overrides[{index}] needs cell_id_exact")

    for index, rule in enumerate(rules.get("sample_library_rules", [])):
        register_id(rule, f"sample_library_rules[{index}]")
        _require_nonempty_string(rule, "gse", f"sample_library_rules[{index}]")
        match_fields = [
            rule.get("patient_id_regex"),
            rule.get("sample_id_exact"),
            rule.get("library_id_exact"),
        ]
        if not any(match_fields):
            raise RefinementError(f"sample_library_rules[{index}] has no matcher")
        if rule.get("patient_id_regex"):
            re.compile(str(rule["patient_id_regex"]), flags=re.IGNORECASE)

    for index, mapping in enumerate(rules.get("sample_value_maps", [])):
        register_id(mapping, f"sample_value_maps[{index}]")
        _require_nonempty_string(mapping, "gse", f"sample_value_maps[{index}]")
        key_fields = mapping.get("key_fields")
        if not isinstance(key_fields, list) or not key_fields:
            raise RefinementError(f"sample_value_maps[{index}].key_fields must not be empty")
        seen_keys: set[str] = set()
        for group_index, group in enumerate(mapping.get("value_groups", [])):
            register_id(group, f"sample_value_maps[{index}].value_groups[{group_index}]")
            _require_nonempty_string(group, "value", f"sample_value_maps[{index}].value_groups[{group_index}]")
            for key in group.get("keys", []):
                normalized = normalize_match_text(key)
                if not normalized:
                    raise RefinementError("Sample-map keys must not be missing-like")
                if normalized in seen_keys:
                    raise RefinementError(
                        f"Ambiguous key {key!r} in sample_value_maps[{index}]"
                    )
                seen_keys.add(normalized)
        unmatched = mapping.get("unmatched")
        if unmatched:
            register_id(unmatched, f"sample_value_maps[{index}].unmatched")

    for index, rule in enumerate(rules.get("gse_tissue_rules", [])):
        register_id(rule, f"gse_tissue_rules[{index}]")
        _require_nonempty_string(rule, "gse", f"gse_tissue_rules[{index}]")
        _require_nonempty_string(rule, "value", f"gse_tissue_rules[{index}]")
        if not rule.get("tissue_values"):
            raise RefinementError(f"gse_tissue_rules[{index}].tissue_values must not be empty")

    seen_tumor_gses: set[str] = set()
    forbidden_generic = {
        "blood",
        "tumor",
        "tumour",
        "adjacent_normal_tissue",
        "adjacent_tissue",
        "lymph_node",
        "lymph_node_metastasis",
        "unresolved",
    }
    for index, requirement in enumerate(rules.get("tumor_project_requirements", [])):
        context = f"tumor_project_requirements[{index}]"
        gse = _require_nonempty_string(requirement, "gse", context)
        gse_norm = normalize_match_text(gse)
        if gse_norm in seen_tumor_gses:
            raise RefinementError(f"Duplicate tumor-project requirement for {gse}")
        seen_tumor_gses.add(gse_norm)
        allowed = requirement.get("allowed_refined_values")
        if not isinstance(allowed, list) or not allowed or not all(
            isinstance(value, str) and value.strip() for value in allowed
        ):
            raise RefinementError(f"{context}.allowed_refined_values must be a non-empty string list")
        generic = sorted({normalize_match_text(value) for value in allowed} & forbidden_generic)
        if generic:
            raise RefinementError(
                f"{context} permits generic tumor-project labels: {generic}"
            )

    for pattern in rules.get("suspicious_value_patterns", []):
        re.compile(str(pattern), flags=re.IGNORECASE)


def load_rules(path: Path) -> dict[str, object]:
    """Load and validate one JSON ruleset."""
    with path.open("r", encoding="utf-8") as handle:
        rules = json.load(handle)
    validate_rules(rules)
    return rules


class SampleSourceRefiner:
    """Resolve one metadata row according to the fixed precedence contract."""

    def __init__(self, rules: Mapping[str, object]):
        validate_rules(rules)
        self.rules = rules
        self.unresolved = str(rules["unresolved_label"])
        self.missing = {normalize_match_text(value) for value in rules.get("missing_values", [])}
        self.suspicious = [
            re.compile(str(pattern), flags=re.IGNORECASE)
            for pattern in rules.get("suspicious_value_patterns", [])
        ]
        self.cell_aliases = {
            normalize_match_text(key): str(value)
            for key, value in rules.get("cell_level_value_aliases", {}).items()
        }
        self.global_aliases = {
            normalize_match_text(key): str(value)
            for key, value in rules.get("global_tissue_aliases", {}).items()
        }
        self.cell_overrides = [self._compile_cell_override(rule) for rule in rules.get("cell_level_overrides", [])]
        self.sample_rules = [self._compile_sample_rule(rule) for rule in rules.get("sample_library_rules", [])]
        self.sample_maps = [self._compile_sample_map(rule) for rule in rules.get("sample_value_maps", [])]
        self.gse_tissue_rules = [self._compile_gse_tissue_rule(rule) for rule in rules.get("gse_tissue_rules", [])]

    def clean(self, value: object) -> str:
        text = clean_scalar(value)
        return "" if normalize_match_text(text) in self.missing else text

    def usable(self, value: object) -> bool:
        text = self.clean(value)
        return bool(text) and not any(pattern.search(text) for pattern in self.suspicious)

    @staticmethod
    def _compile_cell_override(rule: Mapping[str, object]) -> dict[str, object]:
        return {
            **rule,
            "gse_norm": normalize_match_text(rule["gse"]),
            "cell_ids_norm": {normalize_match_text(value) for value in rule.get("cell_id_exact", [])},
        }

    @staticmethod
    def _compile_sample_rule(rule: Mapping[str, object]) -> dict[str, object]:
        return {
            **rule,
            "gse_norm": normalize_match_text(rule["gse"]),
            "sample_ids_norm": {normalize_match_text(value) for value in rule.get("sample_id_exact", [])},
            "library_ids_norm": {normalize_match_text(value) for value in rule.get("library_id_exact", [])},
            "patient_pattern": re.compile(str(rule["patient_id_regex"]), flags=re.IGNORECASE)
            if rule.get("patient_id_regex")
            else None,
        }

    @staticmethod
    def _compile_sample_map(mapping: Mapping[str, object]) -> dict[str, object]:
        lookup: dict[str, tuple[str, str, str]] = {}
        for group in mapping.get("value_groups", []):
            for key in group.get("keys", []):
                lookup[normalize_match_text(key)] = (
                    str(group["value"]),
                    str(group["id"]),
                    str(group.get("reason", "exact sample/library mapping")),
                )
        return {
            **mapping,
            "gse_norm": normalize_match_text(mapping["gse"]),
            "tissues_norm": {normalize_match_text(value) for value in mapping.get("tissue_values", [])},
            "lookup": lookup,
        }

    @staticmethod
    def _compile_gse_tissue_rule(rule: Mapping[str, object]) -> dict[str, object]:
        return {
            **rule,
            "gse_norm": normalize_match_text(rule["gse"]),
            "tissues_norm": {normalize_match_text(value) for value in rule.get("tissue_values", [])},
        }

    def _terminal_or_value(self, rule: Mapping[str, object], level: str) -> Refinement:
        terminal = bool(rule.get("terminal_unresolved"))
        value = self.unresolved if terminal else str(rule["value"])
        return Refinement(
            value=value,
            level=level,
            rule_id=str(rule["id"]),
            reason=str(rule.get("reason", "deterministic rule match")),
        )

    def resolve(self, row: Mapping[str, object]) -> Refinement:
        """Resolve a row using cell > sample/library > GSE+tissue > tissue."""
        cleaned = {field: self.clean(row.get(field, "")) for field in CANONICAL_FIELDS}
        gse_norm = normalize_match_text(cleaned["source_gse_id"])

        for rule in self.cell_overrides:
            if gse_norm == rule["gse_norm"] and normalize_match_text(cleaned["cell_id"]) in rule["cell_ids_norm"]:
                return self._terminal_or_value(rule, "cell_level")

        cell_source = cleaned["cell_source"]
        if self.usable(cell_source):
            normalized = normalize_match_text(cell_source)
            value = self.cell_aliases.get(normalized)
            if value is None and self.rules.get("cell_level_passthrough"):
                value = slugify(cell_source)
            if value:
                return Refinement(
                    value=value,
                    level="cell_level",
                    rule_id="cell_level_explicit_source",
                    reason="Non-missing explicit cell-level source field.",
                )

        sample_norm = normalize_match_text(cleaned["sample_id"])
        library_norm = normalize_match_text(cleaned["library_id"])
        for rule in self.sample_rules:
            if gse_norm != rule["gse_norm"]:
                continue
            matches = bool(sample_norm and sample_norm in rule["sample_ids_norm"])
            matches = matches or bool(library_norm and library_norm in rule["library_ids_norm"])
            patient_pattern = rule["patient_pattern"]
            matches = matches or bool(patient_pattern and patient_pattern.search(cleaned["patient_id"]))
            if matches:
                return self._terminal_or_value(rule, "sample_library")

        tissue_norm = normalize_match_text(cleaned["tissue"])
        for mapping in self.sample_maps:
            if gse_norm != mapping["gse_norm"]:
                continue
            if mapping["tissues_norm"] and tissue_norm not in mapping["tissues_norm"]:
                continue
            for field_name in mapping["key_fields"]:
                key = normalize_match_text(cleaned.get(str(field_name), ""))
                if key and key in mapping["lookup"]:
                    value, rule_id, reason = mapping["lookup"][key]
                    return Refinement(value, "sample_library", rule_id, reason)
            unmatched = mapping.get("unmatched")
            if unmatched and unmatched.get("terminal_unresolved"):
                return self._terminal_or_value(unmatched, "sample_library")

        for rule in self.gse_tissue_rules:
            if gse_norm == rule["gse_norm"] and tissue_norm in rule["tissues_norm"]:
                return self._terminal_or_value(rule, "gse_tissue")

        tissue = cleaned["tissue"]
        if self.usable(tissue):
            value = self.global_aliases.get(tissue_norm)
            rule_id = "global_tissue_alias"
            if value is None and self.rules.get("global_tissue_passthrough"):
                value = slugify(tissue)
                rule_id = "global_tissue_passthrough"
            if value:
                return Refinement(
                    value=value,
                    level="global_tissue_fallback",
                    rule_id=rule_id,
                    reason="Deterministic fallback from the original tissue field.",
                )

        return Refinement(
            value=self.unresolved,
            level="unresolved",
            rule_id="unresolved_no_supported_source",
            reason="No supported cell, sample/library, GSE+tissue, or tissue fallback value.",
        )


def _decode_attr(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _dataset_values(dataset: h5py.Dataset, start: int, stop: int) -> np.ndarray:
    """Read one HDF5 dataset slice as Python-like scalar objects."""
    if dataset.dtype.kind in {"O", "S", "U"}:
        try:
            return np.asarray(dataset.asstr()[start:stop], dtype=object)
        except (TypeError, ValueError):
            pass
    return np.asarray(dataset[start:stop], dtype=object)


class H5adObsReader:
    """Read selected AnnData ``obs`` columns without loading ``X``."""

    def __init__(self, path: Path, rules: Mapping[str, object], mode: str = "r"):
        self.path = path
        self.rules = rules
        self.handle = h5py.File(path, mode)
        if "obs" not in self.handle or not isinstance(self.handle["obs"], h5py.Group):
            self.handle.close()
            raise RefinementError(f"H5AD has no obs group: {path}")
        self.obs = self.handle["obs"]
        index_name = _decode_attr(self.obs.attrs.get("_index", "_index"))
        if index_name not in self.obs and "_index" in self.obs:
            index_name = "_index"
        if index_name not in self.obs:
            self.handle.close()
            raise RefinementError("The obs index dataset referenced by obs.attrs['_index'] is missing")
        self.index_name = index_name
        self.n_obs = int(self.obs[index_name].shape[0])
        self.source_columns: dict[str, list[str]] = {}
        for canonical, aliases in rules["source_columns"].items():
            resolved: list[str] = []
            for alias in aliases:
                physical_name = index_name if alias == "_index" else alias
                if physical_name in self.obs and physical_name not in resolved:
                    resolved.append(physical_name)
            self.source_columns[canonical] = resolved
        if not self.source_columns.get("obs_name"):
            self.source_columns["obs_name"] = [index_name]
        if not self.source_columns.get("source_gse_id"):
            self.handle.close()
            raise RefinementError("No source-GSE column from the configured aliases is present")
        self._categories: dict[str, np.ndarray] = {}
        self.missing = {normalize_match_text(value) for value in rules.get("missing_values", [])}

    def close(self) -> None:
        self.handle.close()

    def __enter__(self) -> "H5adObsReader":
        return self

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        self.close()

    def _read_node(self, name: str, start: int, stop: int) -> np.ndarray:
        node = self.obs[name]
        if isinstance(node, h5py.Dataset):
            return _dataset_values(node, start, stop)

        encoding = _decode_attr(node.attrs.get("encoding-type", ""))
        if encoding == "categorical":
            if name not in self._categories:
                categories = node["categories"]
                self._categories[name] = _dataset_values(categories, 0, len(categories))
            codes = np.asarray(node["codes"][start:stop], dtype=np.int64)
            output = np.full(len(codes), "", dtype=object)
            valid = codes >= 0
            if np.any(valid):
                output[valid] = self._categories[name][codes[valid]]
            return output
        if encoding in {"nullable-string-array", "nullable-integer", "nullable-boolean"}:
            values = _dataset_values(node["values"], start, stop)
            mask = np.asarray(node["mask"][start:stop], dtype=bool)
            values[mask] = ""
            return values
        raise RefinementError(f"Unsupported obs encoding for {name}: {encoding or 'unknown'}")

    def _clean_array(self, values: np.ndarray) -> np.ndarray:
        output = np.empty(len(values), dtype=object)
        for index, value in enumerate(values):
            text = clean_scalar(value)
            output[index] = "" if normalize_match_text(text) in self.missing else text
        return output

    def read_chunk(self, start: int, stop: int) -> dict[str, np.ndarray]:
        physical = {
            name
            for aliases in self.source_columns.values()
            for name in aliases
        }
        cache = {name: self._clean_array(self._read_node(name, start, stop)) for name in physical}
        length = stop - start
        result: dict[str, np.ndarray] = {}
        for canonical in CANONICAL_FIELDS:
            combined = np.full(length, "", dtype=object)
            for name in self.source_columns.get(canonical, []):
                values = cache[name]
                missing = combined == ""
                assign = missing & (values != "")
                combined[assign] = values[assign]
                if not np.any(combined == ""):
                    break
            result[canonical] = combined
        return result


def chunk_bounds(n_rows: int, chunk_size: int) -> Iterator[tuple[int, int]]:
    if chunk_size <= 0:
        raise RefinementError("chunk_size must be positive")
    for start in range(0, n_rows, chunk_size):
        yield start, min(start + chunk_size, n_rows)


RowConsumer = Callable[[int, Mapping[str, object], Refinement], None]


def scan_h5ad(
    h5ad_path: Path,
    rules: Mapping[str, object],
    *,
    chunk_size: int,
    row_consumer: RowConsumer | None = None,
    unresolved_limit: int = 0,
) -> ScanResult:
    """Stream canonical source fields, resolve rows, and compute stable digests."""
    refiner = SampleSourceRefiner(rules)
    source_digest = hashlib.sha256()
    output_digest = hashlib.sha256()
    resolution_digest = hashlib.sha256()
    value_counts: Counter[tuple[str, str, str]] = Counter()
    gse_totals: Counter[str] = Counter()
    gse_unresolved: Counter[str] = Counter()
    gse_values: dict[str, set[str]] = {}
    tumor_context_counts: Counter[tuple[str, str, bool]] = Counter()
    unresolved_examples: list[list[object]] = []
    tumor_requirements = {
        normalize_match_text(requirement["gse"]): (
            str(requirement["gse"]),
            {str(value) for value in requirement["allowed_refined_values"]},
        )
        for requirement in rules.get("tumor_project_requirements", [])
    }

    with H5adObsReader(h5ad_path, rules) as reader:
        for start, stop in chunk_bounds(reader.n_obs, chunk_size):
            chunk = reader.read_chunk(start, stop)
            for offset in range(stop - start):
                row = {field: chunk[field][offset] for field in CANONICAL_FIELDS}
                resolution = refiner.resolve(row)
                row_index = start + offset
                _update_digest(source_digest, [row_index, *(row[field] for field in CANONICAL_FIELDS)])
                _update_digest(output_digest, [resolution.value])
                _update_digest(
                    resolution_digest,
                    [resolution.value, resolution.level, resolution.rule_id, resolution.reason],
                )
                value_counts[(resolution.value, resolution.level, resolution.rule_id)] += 1
                gse = clean_scalar(row["source_gse_id"])
                gse_totals[gse] += 1
                gse_values.setdefault(gse, set()).add(resolution.value)
                tumor_requirement = tumor_requirements.get(normalize_match_text(gse))
                if tumor_requirement is not None:
                    canonical_gse, allowed = tumor_requirement
                    compliant = resolution.value in allowed
                    tumor_context_counts[(canonical_gse, resolution.value, compliant)] += 1
                if resolution.value == refiner.unresolved:
                    gse_unresolved[gse] += 1
                    if len(unresolved_examples) < unresolved_limit:
                        unresolved_examples.append(
                            [row_index, *(row[field] for field in CANONICAL_FIELDS), resolution.reason]
                        )
                if row_consumer is not None:
                    row_consumer(row_index, row, resolution)
            logging.info("Processed obs rows %s-%s", f"{start:,}", f"{stop:,}")
        source_columns = {key: list(value) for key, value in reader.source_columns.items()}
        n_obs = reader.n_obs

    return ScanResult(
        n_obs=n_obs,
        source_sha256=source_digest.hexdigest(),
        output_values_sha256=output_digest.hexdigest(),
        resolutions_sha256=resolution_digest.hexdigest(),
        source_columns=source_columns,
        value_counts=value_counts,
        gse_totals=gse_totals,
        gse_unresolved=gse_unresolved,
        gse_values=gse_values,
        tumor_context_counts=tumor_context_counts,
        unresolved_examples=unresolved_examples,
    )


def _file_identity(path: Path) -> dict[str, object]:
    stat = path.stat()
    return {
        "realpath": str(path.resolve()),
        "device": int(stat.st_dev),
        "inode": int(stat.st_ino),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


def _same_file_identity(left: Mapping[str, object], right: Mapping[str, object]) -> bool:
    keys = ["realpath", "device", "inode", "size_bytes", "mtime_ns"]
    return all(left.get(key) == right.get(key) for key in keys)


def review_h5ad(
    h5ad_path: Path,
    rules_path: Path,
    output_dir: Path,
    *,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    unresolved_limit: int = 5_000,
) -> Path:
    """Create a full dry-run review export and a writeback-gating manifest."""
    h5ad_path = h5ad_path.expanduser().resolve()
    rules_path = rules_path.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    rules = load_rules(rules_path)
    identity_before = _file_identity(h5ad_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    review_export = output_dir / "sample_source_refined_review.csv.gz"
    with atomic_text_writer(review_export, gzip_output=True) as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(REVIEW_COLUMNS)

        def consume(row_index: int, row: Mapping[str, object], resolution: Refinement) -> None:
            writer.writerow(
                [
                    row_index,
                    *(row[field] for field in CANONICAL_FIELDS),
                    resolution.value,
                    resolution.level,
                    resolution.rule_id,
                    resolution.reason,
                ]
            )

        result = scan_h5ad(
            h5ad_path,
            rules,
            chunk_size=chunk_size,
            row_consumer=consume,
            unresolved_limit=unresolved_limit,
        )

    identity_after = _file_identity(h5ad_path)
    if not _same_file_identity(identity_before, identity_after):
        review_export.unlink(missing_ok=True)
        raise RefinementError("Input H5AD changed during dry-run review generation")

    by_gse_path = output_dir / "sample_source_refined_by_gse.csv"
    by_gse_rows = []
    for gse in sorted(result.gse_totals):
        total = result.gse_totals[gse]
        unresolved = result.gse_unresolved[gse]
        by_gse_rows.append(
            [
                gse,
                total,
                total - unresolved,
                unresolved,
                f"{(total - unresolved) / total:.8f}" if total else "0.00000000",
                len(result.gse_values.get(gse, set())),
            ]
        )
    write_csv(
        by_gse_path,
        [
            "source_gse_id",
            "cell_n",
            "resolved_n",
            "unresolved_n",
            "resolved_fraction",
            "unique_refined_values",
        ],
        by_gse_rows,
    )

    value_counts_path = output_dir / "sample_source_refined_value_counts.csv"
    value_rows = [
        [value, level, rule_id, count]
        for (value, level, rule_id), count in sorted(result.value_counts.items())
    ]
    write_csv(
        value_counts_path,
        ["sample_source_refined", "refinement_level", "refinement_rule_id", "cell_n"],
        value_rows,
    )

    tumor_audit_path = output_dir / "sample_source_refined_tumor_project_audit.csv"
    tumor_rows = [
        [gse, value, compliant, count]
        for (gse, value, compliant), count in sorted(result.tumor_context_counts.items())
    ]
    write_csv(
        tumor_audit_path,
        ["source_gse_id", "sample_source_refined", "context_compliant", "cell_n"],
        tumor_rows,
    )
    tumor_violations = int(
        sum(
            count
            for (_gse, _value, compliant), count in result.tumor_context_counts.items()
            if not compliant
        )
    )

    unresolved_path = output_dir / "sample_source_refined_unresolved_examples.csv.gz"
    write_csv(
        unresolved_path,
        ["row_index", *CANONICAL_FIELDS, "refinement_reason"],
        result.unresolved_examples,
    )

    artifact_paths = [
        review_export,
        by_gse_path,
        value_counts_path,
        tumor_audit_path,
        unresolved_path,
    ]
    manifest = {
        "schema_version": 1,
        "mode": "review_only",
        "output_field": rules["output_field"],
        "ruleset_id": rules["ruleset_id"],
        "rules_path": str(rules_path),
        "rules_sha256": file_sha256(rules_path),
        "input_h5ad": identity_before,
        "obs_rows": result.n_obs,
        "source_columns": result.source_columns,
        "source_values_sha256": result.source_sha256,
        "output_values_sha256": result.output_values_sha256,
        "resolutions_sha256": result.resolutions_sha256,
        "chunk_size": chunk_size,
        "unresolved_cells": int(sum(result.gse_unresolved.values())),
        "tumor_context_gate_pass": tumor_violations == 0,
        "tumor_context_violations": tumor_violations,
        "artifacts": {
            path.name: {
                "path": str(path),
                "sha256": file_sha256(path),
            }
            for path in artifact_paths
        },
    }
    manifest_path = output_dir / "sample_source_refinement_review_manifest.json"
    write_json(manifest_path, manifest)
    logging.info("Review manifest written to %s", manifest_path)
    return manifest_path


def _load_manifest(path: Path) -> dict[str, object]:
    with path.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)
    if manifest.get("schema_version") != 1 or manifest.get("mode") != "review_only":
        raise RefinementError("Review manifest has an unsupported schema or mode")
    return manifest


def _verify_manifest_artifacts(manifest: Mapping[str, object]) -> None:
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, dict) or "sample_source_refined_review.csv.gz" not in artifacts:
        raise RefinementError("Review manifest does not reference the full review export")
    for name, metadata in artifacts.items():
        path = Path(str(metadata["path"]))
        if not path.is_file():
            raise RefinementError(f"Review artifact is missing: {path}")
        if file_sha256(path) != metadata.get("sha256"):
            raise RefinementError(f"Review artifact checksum mismatch: {name}")


def _jsonable_attr(value: object) -> object:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.ndarray):
        return [_jsonable_attr(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def h5ad_structure_sha256(handle: h5py.File, output_field: str) -> str:
    """Hash HDF5 structure outside the owned output/staging columns."""
    digest = hashlib.sha256()
    ignored = {
        f"obs/{output_field}",
        f"obs/{STAGING_COLUMN}",
        f"obs/{BACKUP_COLUMN}",
    }

    def record(name: str, node: h5py.Group | h5py.Dataset) -> None:
        if name in ignored:
            return
        object_info = h5py.h5o.get_info(node.id)
        metadata: dict[str, object] = {
            "name": name,
            "kind": "dataset" if isinstance(node, h5py.Dataset) else "group",
            "address": int(object_info.addr),
        }
        if isinstance(node, h5py.Dataset):
            metadata.update(
                {
                    "shape": list(node.shape),
                    "dtype": str(node.dtype),
                    "chunks": list(node.chunks) if node.chunks else None,
                    "compression": node.compression,
                }
            )
        attributes = {}
        for key in sorted(node.attrs):
            if name == "obs" and key == "column-order":
                continue
            attributes[key] = _jsonable_attr(node.attrs[key])
        metadata["attrs"] = attributes
        _update_digest(digest, [json.dumps(metadata, sort_keys=True, separators=(",", ":"))])

    root_attrs = {key: _jsonable_attr(handle.attrs[key]) for key in sorted(handle.attrs)}
    _update_digest(digest, [json.dumps(root_attrs, sort_keys=True, separators=(",", ":"))])
    handle.visititems(record)
    return digest.hexdigest()


def _column_order(obs: h5py.Group) -> list[str]:
    values = obs.attrs.get("column-order", [])
    return [_decode_attr(value) for value in values]


def _set_column_order(obs: h5py.Group, values: Sequence[str]) -> None:
    if "column-order" in obs.attrs:
        del obs.attrs["column-order"]
    obs.attrs["column-order"] = np.asarray(
        list(values),
        dtype=h5py.string_dtype(encoding="utf-8"),
    )


def _digest_string_dataset(dataset: h5py.Dataset, chunk_size: int) -> str:
    digest = hashlib.sha256()
    for start, stop in chunk_bounds(len(dataset), chunk_size):
        values = _dataset_values(dataset, start, stop)
        for value in values:
            _update_digest(digest, [value])
    return digest.hexdigest()


def _verify_review_for_writeback(
    h5ad_path: Path,
    rules_path: Path,
    manifest: Mapping[str, object],
    *,
    chunk_size: int,
) -> tuple[dict[str, object], ScanResult]:
    rules = load_rules(rules_path)
    if manifest.get("tumor_context_gate_pass") is not True:
        raise RefinementError(
            "Tumor-project specimen-context gate did not pass; writeback is blocked"
        )
    if manifest.get("output_field") != rules["output_field"]:
        raise RefinementError("Manifest output field does not match the rules")
    if file_sha256(rules_path) != manifest.get("rules_sha256"):
        raise RefinementError("Rules changed after the review export was created")
    current_identity = _file_identity(h5ad_path)
    if not _same_file_identity(current_identity, manifest.get("input_h5ad", {})):
        raise RefinementError("Input H5AD identity changed after review")
    _verify_manifest_artifacts(manifest)

    result = scan_h5ad(h5ad_path, rules, chunk_size=chunk_size)
    checks = {
        "obs_rows": result.n_obs,
        "source_values_sha256": result.source_sha256,
        "output_values_sha256": result.output_values_sha256,
        "resolutions_sha256": result.resolutions_sha256,
        "source_columns": result.source_columns,
        "tumor_context_violations": int(
            sum(
                count
                for (_gse, _value, compliant), count in result.tumor_context_counts.items()
                if not compliant
            )
        ),
    }
    for key, current in checks.items():
        if current != manifest.get(key):
            raise RefinementError(f"Review manifest no longer matches {key}")
    return rules, result


def writeback_h5ad(
    h5ad_path: Path,
    rules_path: Path,
    review_manifest_path: Path,
    *,
    confirm_reviewed: bool,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    replace_existing: bool = False,
) -> Path:
    """Write only the reviewed output column after all provenance checks pass."""
    if not confirm_reviewed:
        raise RefinementError("Writeback requires explicit confirm_reviewed=True")
    h5ad_path = h5ad_path.expanduser().resolve()
    rules_path = rules_path.expanduser().resolve()
    review_manifest_path = review_manifest_path.expanduser().resolve()
    manifest = _load_manifest(review_manifest_path)
    rules, reviewed = _verify_review_for_writeback(
        h5ad_path,
        rules_path,
        manifest,
        chunk_size=chunk_size,
    )
    output_field = str(rules["output_field"])
    refiner = SampleSourceRefiner(rules)

    with H5adObsReader(h5ad_path, rules, mode="r+") as reader:
        obs = reader.obs
        if output_field in obs and not replace_existing:
            raise RefinementError(
                f"obs/{output_field} already exists; use replace_existing=True to replace it"
            )
        if STAGING_COLUMN in obs or BACKUP_COLUMN in obs:
            raise RefinementError("A prior staging/backup column is present; inspect it before retrying")

        structure_before = h5ad_structure_sha256(reader.handle, output_field)
        original_order = _column_order(obs)
        had_existing = output_field in obs
        old_moved = False
        new_committed = False
        output_digest = hashlib.sha256()

        try:
            create_kwargs: dict[str, object] = {}
            if reader.n_obs:
                create_kwargs = {
                    "chunks": (min(chunk_size, reader.n_obs),),
                    "compression": "gzip",
                    "compression_opts": 4,
                }
            staging = obs.create_dataset(
                STAGING_COLUMN,
                shape=(reader.n_obs,),
                dtype=h5py.string_dtype(encoding="utf-8"),
                **create_kwargs,
            )
            staging.attrs["encoding-type"] = "string-array"
            staging.attrs["encoding-version"] = "0.2.0"

            for start, stop in chunk_bounds(reader.n_obs, chunk_size):
                chunk = reader.read_chunk(start, stop)
                values = np.empty(stop - start, dtype=object)
                for offset in range(stop - start):
                    row = {field: chunk[field][offset] for field in CANONICAL_FIELDS}
                    value = refiner.resolve(row).value
                    values[offset] = value
                    _update_digest(output_digest, [value])
                staging[start:stop] = values
                logging.info(
                    "Staged refined values for rows %s-%s",
                    f"{start:,}",
                    f"{stop:,}",
                )
            reader.handle.flush()

            if output_digest.hexdigest() != reviewed.output_values_sha256:
                raise RefinementError("Staged values differ from the reviewed values")
            if _digest_string_dataset(staging, chunk_size) != reviewed.output_values_sha256:
                raise RefinementError("Staging dataset failed readback verification")

            if had_existing:
                obs.move(output_field, BACKUP_COLUMN)
                old_moved = True
            obs.move(STAGING_COLUMN, output_field)
            new_committed = True
            updated_order = [name for name in original_order if name != output_field]
            updated_order.append(output_field)
            _set_column_order(obs, updated_order)
            reader.handle.flush()

            if _digest_string_dataset(obs[output_field], chunk_size) != reviewed.output_values_sha256:
                raise RefinementError("Committed output column failed readback verification")
            structure_after = h5ad_structure_sha256(reader.handle, output_field)
            if structure_after != structure_before:
                raise RefinementError("H5AD structure outside the owned output column changed")

            if old_moved:
                del obs[BACKUP_COLUMN]
                old_moved = False
            reader.handle.flush()
        except Exception:
            if new_committed and output_field in obs:
                del obs[output_field]
            if old_moved and BACKUP_COLUMN in obs:
                obs.move(BACKUP_COLUMN, output_field)
            if STAGING_COLUMN in obs:
                del obs[STAGING_COLUMN]
            _set_column_order(obs, original_order)
            reader.handle.flush()
            raise

    writeback_report = {
        "schema_version": 1,
        "mode": "writeback_complete",
        "input_h5ad": str(h5ad_path),
        "review_manifest": str(review_manifest_path),
        "rules_sha256": manifest["rules_sha256"],
        "output_field": output_field,
        "obs_rows": reviewed.n_obs,
        "output_values_sha256": reviewed.output_values_sha256,
        "replaced_existing_output": had_existing,
        "preserved_structure_sha256": structure_before,
    }
    report_path = review_manifest_path.parent / "sample_source_refinement_writeback.json"
    write_json(report_path, writeback_report)
    logging.info("Writeback report written to %s", report_path)
    return report_path


def _add_common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--input-h5ad", default=str(DEFAULT_H5AD), help="Input H5AD path.")
    parser.add_argument("--rules", default=str(DEFAULT_RULES), help="Refinement rules JSON path.")
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=DEFAULT_CHUNK_SIZE,
        help="Number of obs rows to process per bounded chunk.",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable progress logging.")


def build_parser() -> argparse.ArgumentParser:
    """Build the reproducible two-step CLI."""
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    review = subparsers.add_parser(
        "review",
        aliases=["dry-run"],
        help="Create non-mutating review exports and a writeback manifest.",
    )
    _add_common_arguments(review)
    review.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR), help="Review output directory.")
    review.add_argument(
        "--max-unresolved-examples",
        type=int,
        default=5_000,
        help="Maximum unresolved rows retained in the examples export.",
    )

    writeback = subparsers.add_parser(
        "writeback",
        help="Write the reviewed field after manifest validation.",
    )
    _add_common_arguments(writeback)
    writeback.add_argument("--review-manifest", required=True, help="Manifest created by review/dry-run.")
    writeback.add_argument(
        "--confirm-reviewed",
        action="store_true",
        help="Confirm that the dry-run artifacts were reviewed.",
    )
    writeback.add_argument(
        "--replace-existing",
        action="store_true",
        help="Replace an existing sample_source_refined column after full validation.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point."""
    args = build_parser().parse_args(argv)
    configure_logging(args.verbose)
    if args.command in {"review", "dry-run"}:
        manifest = review_h5ad(
            Path(args.input_h5ad),
            Path(args.rules),
            Path(args.output_dir),
            chunk_size=args.chunk_size,
            unresolved_limit=args.max_unresolved_examples,
        )
        print(manifest)
        return 0
    report = writeback_h5ad(
        Path(args.input_h5ad),
        Path(args.rules),
        Path(args.review_manifest),
        confirm_reviewed=args.confirm_reviewed,
        chunk_size=args.chunk_size,
        replace_existing=args.replace_existing,
    )
    print(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
