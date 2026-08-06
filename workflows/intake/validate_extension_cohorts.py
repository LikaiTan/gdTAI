#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Validate manifest-driven extension cohort intake without touching canonical data."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_COHORTS = PROJECT_ROOT / "configs/datasets/extension_cohorts.csv"
DEFAULT_LIBRARIES = PROJECT_ROOT / "configs/datasets/extension_libraries.csv"
DEFAULT_CANONICAL_REGISTRY = PROJECT_ROOT / "configs/datasets/datasets.csv"
DEFAULT_SHARED_METASHEET = (
    PROJECT_ROOT / "data/compat/downloads/new_candidata_data/sample_metasheet.csv"
)

COHORT_COLUMNS = (
    "cohort_id",
    "cohort_key",
    "cohort_name",
    "primary_accession",
    "gex_accessions",
    "vdj_accessions",
    "publication_ids",
    "bioproject_ids",
    "technology",
    "default_tissue",
    "default_diagnosis",
    "builder_adapter",
    "tcr_schema",
    "canonical_tcr_flags",
    "count_provenance",
    "sparse_counts",
    "preserve_raw_counts",
    "intake_role",
    "holdout_type",
    "sealed",
    "stage_enabled",
    "build_enabled",
    "duplicate_of",
    "exclusion_reason",
    "notes",
)

LIBRARY_COLUMNS = (
    "cohort_id",
    "library_id",
    "source_accession",
    "source_role",
    "modality",
    "builder_adapter",
    "source_glob",
    "required",
    "metasheet_required",
    "sample_id_rule",
    "library_id_rule",
    "donor_rule",
    "tissue_rule",
    "diagnosis_rule",
    "tcr_schema",
    "tcr_join_key",
    "count_matrix_state",
    "sparse_expected",
    "notes",
)

CANONICAL_TCR_FLAGS = (
    "has_TRA",
    "has_TRB",
    "has_TRG",
    "has_TRD",
    "has_TRA_TRB_paired",
    "has_TRG_TRD_paired",
    "has_any_ab_tcr",
    "has_any_gd_tcr",
)

REQUIRED_INCLUDED_COHORTS = {
    "GSE114724",
    "GSE121636_GSE121637",
    "GSE159251",
    "GSE169246",
    "GSE292700",
    "GSE294273_GSE294274",
    "GSE296954",
    "GSE315928",
}

TRUE_FALSE = {"true", "false"}
COUNT_ROLES = {"gex", "combined"}
ALLOWED_SOURCE_ROLES = {
    "gex",
    "vdj",
    "combined",
    "cell_metadata",
    "publication_annotation",
    "sample_metadata",
}
PROTECTED_ROOTS = (
    PROJECT_ROOT / "Integrated_dataset",
    PROJECT_ROOT / "data/datasets",
    PROJECT_ROOT / "downloads/per_gse_h5ad_with_metadata",
)


class ManifestError(ValueError):
    """Raised when a manifest cannot be parsed safely."""


@dataclass(frozen=True)
class ValidationIssue:
    severity: str
    code: str
    message: str
    cohort_id: str = ""
    library_id: str = ""


@dataclass
class ValidationReport:
    issues: list[ValidationIssue]
    cohort_count: int = 0
    included_cohort_count: int = 0
    library_count: int = 0
    h5ad_count: int = 0

    @property
    def errors(self) -> list[ValidationIssue]:
        return [issue for issue in self.issues if issue.severity == "error"]

    @property
    def warnings(self) -> list[ValidationIssue]:
        return [issue for issue in self.issues if issue.severity == "warning"]

    @property
    def passed(self) -> bool:
        return not self.errors

    def to_dict(self) -> dict[str, object]:
        return {
            "passed": self.passed,
            "cohort_count": self.cohort_count,
            "included_cohort_count": self.included_cohort_count,
            "library_count": self.library_count,
            "h5ad_count": self.h5ad_count,
            "error_count": len(self.errors),
            "warning_count": len(self.warnings),
            "issues": [asdict(issue) for issue in self.issues],
        }


def split_tokens(value: str) -> list[str]:
    """Split semicolon-delimited manifest values and remove empty tokens."""
    return [token.strip() for token in value.split(";") if token.strip()]


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized not in TRUE_FALSE:
        raise ManifestError(f"Expected true/false, received {value!r}")
    return normalized == "true"


def read_manifest(path: Path, expected_columns: Sequence[str]) -> list[dict[str, str]]:
    """Read one CSV manifest with an exact, stable header."""
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != tuple(expected_columns):
            raise ManifestError(
                f"Unexpected header in {path}: expected {list(expected_columns)}, "
                f"received {reader.fieldnames}"
            )
        rows: list[dict[str, str]] = []
        for line_number, row in enumerate(reader, start=2):
            if None in row:
                raise ManifestError(f"Extra CSV fields in {path}:{line_number}")
            cleaned = {key: (value or "").strip() for key, value in row.items()}
            if not any(cleaned.values()):
                continue
            cleaned["_line"] = str(line_number)
            rows.append(cleaned)
    return rows


def cohort_accessions(row: dict[str, str]) -> set[str]:
    values = [row["primary_accession"]]
    values.extend(split_tokens(row["gex_accessions"]))
    values.extend(split_tokens(row["vdj_accessions"]))
    return {value.upper() for value in values if value}


def _is_included(row: dict[str, str]) -> bool:
    return row["intake_role"] != "excluded_duplicate"


def _add_duplicate_issues(
    rows: Iterable[dict[str, str]],
    field_name: str,
    values_for_row,
    issues: list[ValidationIssue],
) -> None:
    owners: dict[str, str] = {}
    for row in rows:
        for raw_value in values_for_row(row):
            value = re.sub(r"\s+", " ", raw_value.strip()).casefold()
            if not value:
                continue
            previous = owners.get(value)
            if previous and previous != row["cohort_id"]:
                issues.append(
                    ValidationIssue(
                        "error",
                        f"duplicate_{field_name}",
                        f"{field_name} {raw_value!r} is shared with cohort {previous}",
                        cohort_id=row["cohort_id"],
                    )
                )
            else:
                owners[value] = row["cohort_id"]


def _canonical_identifiers(path: Path) -> set[str]:
    if not path.exists():
        raise ManifestError(f"Canonical registry not found: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"dataset_id", "accession"}
        if not required.issubset(reader.fieldnames or []):
            raise ManifestError(f"Canonical registry lacks {sorted(required)}: {path}")
        identifiers: set[str] = set()
        for row in reader:
            identifiers.update(
                value.strip().upper()
                for value in (row.get("dataset_id", ""), row.get("accession", ""))
                if value and value.strip()
            )
    return identifiers


def _validate_cohort_rows(
    rows: list[dict[str, str]],
    canonical_ids: set[str],
    issues: list[ValidationIssue],
) -> None:
    seen_ids: set[str] = set()
    included = [row for row in rows if _is_included(row)]

    for row in rows:
        cohort_id = row["cohort_id"]
        if not cohort_id:
            issues.append(ValidationIssue("error", "missing_cohort_id", "cohort_id is blank"))
            continue
        if cohort_id in seen_ids:
            issues.append(
                ValidationIssue("error", "duplicate_cohort_id", "cohort_id is repeated", cohort_id)
            )
        seen_ids.add(cohort_id)

        for field in (
            "sparse_counts",
            "preserve_raw_counts",
            "sealed",
            "stage_enabled",
            "build_enabled",
        ):
            if row[field].lower() not in TRUE_FALSE:
                issues.append(
                    ValidationIssue(
                        "error",
                        "invalid_boolean",
                        f"{field} must be true or false, received {row[field]!r}",
                        cohort_id,
                    )
                )

        if not row["cohort_key"] or not row["primary_accession"]:
            issues.append(
                ValidationIssue(
                    "error",
                    "missing_identity",
                    "cohort_key and primary_accession are required",
                    cohort_id,
                )
            )

        flags = split_tokens(row["canonical_tcr_flags"])
        if flags != list(CANONICAL_TCR_FLAGS):
            issues.append(
                ValidationIssue(
                    "error",
                    "canonical_tcr_flags",
                    "canonical_tcr_flags must contain the complete ordered TRA/TRB/TRG/TRD flag set",
                    cohort_id,
                )
            )

        if _is_included(row):
            if row["sparse_counts"] != "true" or row["preserve_raw_counts"] != "true":
                issues.append(
                    ValidationIssue(
                        "error",
                        "raw_sparse_provenance",
                        "included cohorts must preserve sparse raw integer counts",
                        cohort_id,
                    )
                )
            conflicts = cohort_accessions(row) & canonical_ids
            if conflicts:
                issues.append(
                    ValidationIssue(
                        "error",
                        "canonical_accession_duplicate",
                        f"already present in canonical registry: {sorted(conflicts)}",
                        cohort_id,
                    )
                )
        else:
            disabled = row["stage_enabled"] == "false" and row["build_enabled"] == "false"
            duplicate_of = row["duplicate_of"].upper()
            if not disabled or not duplicate_of or duplicate_of not in canonical_ids:
                issues.append(
                    ValidationIssue(
                        "error",
                        "invalid_duplicate_exclusion",
                        "excluded duplicates must disable stage/build and reference a canonical dataset",
                        cohort_id,
                    )
                )

    present_included = {row["cohort_id"] for row in included}
    missing = REQUIRED_INCLUDED_COHORTS - present_included
    extras = present_included - REQUIRED_INCLUDED_COHORTS
    if missing:
        issues.append(
            ValidationIssue(
                "error", "missing_required_cohorts", f"missing required cohorts: {sorted(missing)}"
            )
        )
    if extras:
        issues.append(
            ValidationIssue(
                "warning", "additional_cohorts", f"additional included cohorts: {sorted(extras)}"
            )
        )

    _add_duplicate_issues(included, "cohort", lambda row: [row["cohort_key"]], issues)
    _add_duplicate_issues(included, "accession", cohort_accessions, issues)
    _add_duplicate_issues(
        included, "publication", lambda row: split_tokens(row["publication_ids"]), issues
    )
    _add_duplicate_issues(
        included, "bioproject", lambda row: split_tokens(row["bioproject_ids"]), issues
    )

    policies = {
        "GSE315928": ("sealed_abt_negative_holdout", "abt_negative", "embedded_airr_tra_trb"),
        "GSE121636_GSE121637": (
            "sealed_mixed_t_nk_holdout",
            "mixed_t_nk",
            "productive_contigs",
        ),
        "GSE159251": (
            "extension_candidate",
            "none",
            "partial_embedded_tra_trb_cdr3",
        ),
    }
    by_id = {row["cohort_id"]: row for row in rows}
    for cohort_id, (role, holdout, schema) in policies.items():
        row = by_id.get(cohort_id)
        if not row:
            continue
        expected_sealed = cohort_id != "GSE159251"
        if (
            row["intake_role"] != role
            or row["holdout_type"] != holdout
            or row["tcr_schema"] != schema
            or (row["sealed"] == "true") != expected_sealed
        ):
            issues.append(
                ValidationIssue(
                    "error",
                    "cohort_policy",
                    f"required policy for {cohort_id} is not encoded exactly",
                    cohort_id,
                )
            )

    gse169246 = by_id.get("GSE169246")
    if gse169246 and (
        gse169246["intake_role"] != "blocked_provenance"
        or gse169246["builder_adapter"] != "provenance_blocked"
        or gse169246["tcr_schema"] != "unverified_blocked"
        or gse169246["stage_enabled"] != "false"
        or gse169246["build_enabled"] != "false"
    ):
        issues.append(
            ValidationIssue(
                "error",
                "gse169246_provenance_block",
                "GSE169246 must remain stage/build disabled until accession-matched provenance is proven",
                "GSE169246",
            )
        )

    gse294 = by_id.get("GSE294273_GSE294274")
    if gse294 and (
        "melanoma" not in gse294["cohort_name"].casefold()
        or gse294["default_tissue"] != "lymph_node_metastasis"
        or gse294["default_diagnosis"] != "melanoma"
        or "esophageal" in " ".join(gse294.values()).casefold()
    ):
        issues.append(
            ValidationIssue(
                "error",
                "gse294_metadata_policy",
                "GSE294273/GSE294274 must be encoded as melanoma lymph-node metastases",
                "GSE294273_GSE294274",
            )
        )

    tan_rows = [row for row in rows if row["duplicate_of"] == "GDT_2020AUG_woCOV"]
    if len(tan_rows) != 1 or "Tan et al. 2021" not in tan_rows[0].get("notes", ""):
        issues.append(
            ValidationIssue(
                "error",
                "tan_2021_exclusion",
                "exactly one disabled GDT_2020AUG_woCOV / Tan et al. 2021 exclusion is required",
            )
        )


def _validate_library_rows(
    cohorts: list[dict[str, str]],
    libraries: list[dict[str, str]],
    issues: list[ValidationIssue],
) -> None:
    cohort_by_id = {row["cohort_id"]: row for row in cohorts}
    seen_library_ids: set[str] = set()
    libraries_by_cohort: dict[str, list[dict[str, str]]] = {}

    for row in libraries:
        cohort_id = row["cohort_id"]
        library_id = row["library_id"]
        cohort = cohort_by_id.get(cohort_id)
        if not cohort:
            issues.append(
                ValidationIssue(
                    "error", "orphan_library", "cohort_id is not in extension_cohorts.csv", cohort_id, library_id
                )
            )
            continue
        libraries_by_cohort.setdefault(cohort_id, []).append(row)
        if library_id in seen_library_ids:
            issues.append(
                ValidationIssue(
                    "error", "duplicate_library_id", "library_id is repeated", cohort_id, library_id
                )
            )
        seen_library_ids.add(library_id)

        if row["source_role"] not in ALLOWED_SOURCE_ROLES:
            issues.append(
                ValidationIssue(
                    "error",
                    "invalid_source_role",
                    f"unsupported source_role {row['source_role']!r}",
                    cohort_id,
                    library_id,
                )
            )
        if row["required"] not in TRUE_FALSE or row["metasheet_required"] not in TRUE_FALSE:
            issues.append(
                ValidationIssue(
                    "error", "invalid_boolean", "required/metasheet_required must be true/false", cohort_id, library_id
                )
            )
        if row["builder_adapter"] != cohort["builder_adapter"]:
            issues.append(
                ValidationIssue(
                    "error", "adapter_mismatch", "library and cohort builder adapters differ", cohort_id, library_id
                )
            )
        if row["source_accession"].upper() not in cohort_accessions(cohort):
            issues.append(
                ValidationIssue(
                    "error",
                    "accession_not_linked",
                    "source_accession is not listed on its biological cohort",
                    cohort_id,
                    library_id,
                )
            )
        if row["tcr_join_key"] != "sample_id+barcode_core":
            issues.append(
                ValidationIssue(
                    "error", "unsafe_tcr_join_key", "TCR joins must use sample_id+barcode_core", cohort_id, library_id
                )
            )
        if row["source_role"] in COUNT_ROLES:
            if row["count_matrix_state"] != "raw_integer_counts" or row["sparse_expected"] != "true":
                issues.append(
                    ValidationIssue(
                        "error",
                        "raw_sparse_count_source",
                        "GEX/combined sources must be sparse raw integer counts",
                        cohort_id,
                        library_id,
                    )
                )

    for cohort in cohorts:
        if not _is_included(cohort):
            continue
        cohort_libraries = libraries_by_cohort.get(cohort["cohort_id"], [])
        if cohort["cohort_id"] == "GSE169246":
            forbidden = [
                row
                for row in cohort_libraries
                if any(token in row["source_glob"].casefold() for token in ("bootstrap", "mmc3", "mmc5"))
                or row["source_role"] != "gex"
                or row["required"] != "false"
            ]
            if forbidden:
                issues.append(
                    ValidationIssue(
                        "error",
                        "gse169246_forbidden_source",
                        "GSE169246 may not accept bootstrap, mmc3/mmc5, or TCR/metadata sources",
                        cohort["cohort_id"],
                    )
                )
            continue
        count_sources = [row for row in cohort_libraries if row["source_role"] in COUNT_ROLES]
        if not count_sources or not any(row["required"] == "true" for row in count_sources):
            issues.append(
                ValidationIssue(
                    "error",
                    "missing_count_source",
                    "included cohort has no required GEX/combined source",
                    cohort["cohort_id"],
                )
            )


def _is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def assert_extension_output_root(path: Path) -> Path:
    """Reject canonical/protected destinations for extension intake artifacts."""
    resolved = path.expanduser().resolve()
    for protected in PROTECTED_ROOTS:
        protected_resolved = protected.resolve()
        if resolved == protected_resolved or _is_relative_to(resolved, protected_resolved):
            raise ManifestError(f"Extension intake cannot write under protected path: {resolved}")
    return resolved


def _validate_shared_metasheet(
    cohorts: Sequence[dict[str, str]],
    path: Path,
    issues: list[ValidationIssue],
) -> None:
    required_columns = {"gse", "gsm", "sample_title", "modality"}
    try:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            if not required_columns.issubset(reader.fieldnames or []):
                raise ManifestError(
                    f"shared metasheet lacks {sorted(required_columns)}: {path}"
                )
            rows = [
                {key: (value or "").strip() for key, value in row.items()}
                for row in reader
            ]
    except (OSError, ManifestError) as exc:
        issues.append(ValidationIssue("error", "shared_metasheet", str(exc)))
        return

    gsm_values = [row["gsm"] for row in rows]
    if any(not value for value in gsm_values) or len(gsm_values) != len(set(gsm_values)):
        issues.append(
            ValidationIssue(
                "error",
                "shared_metasheet_gsm",
                "shared metasheet GSM values must be nonblank and unique",
            )
        )
    by_accession: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_accession.setdefault(row["gse"].upper(), []).append(row)
    for cohort in cohorts:
        if cohort["stage_enabled"] != "true":
            continue
        missing = {
            accession
            for accession in cohort_accessions(cohort)
            if accession.startswith("GSE") and accession not in by_accession
        }
        if missing:
            issues.append(
                ValidationIssue(
                    "error",
                    "shared_metasheet_accession",
                    f"shared metasheet lacks linked accession(s): {sorted(missing)}",
                    cohort["cohort_id"],
                )
            )

    melanoma_rows = by_accession.get("GSE294273", []) + by_accession.get("GSE294274", [])
    if melanoma_rows and any(
        "melanoma" not in row.get("cancer", "").casefold()
        or "lymph" not in " ".join(
            (row.get("source", ""), row.get("char:tissue", ""), row.get("cancer", ""))
        ).casefold()
        for row in melanoma_rows
    ):
        issues.append(
            ValidationIssue(
                "error",
                "gse294_shared_metadata",
                "shared GSE294273/GSE294274 rows are not consistently melanoma lymph-node records",
                "GSE294273_GSE294274",
            )
        )


def validate_built_h5ad(path: Path, cohort: Mapping[str, str]) -> list[ValidationIssue]:
    """Validate one serialized extension H5AD in backed/chunked mode."""
    import anndata as ad
    import h5py
    import numpy as np
    import pandas as pd

    issues: list[ValidationIssue] = []
    cohort_id = cohort["cohort_id"]
    tcr_columns = [
        f"{chain}_{suffix}"
        for chain in ("TRA", "TRB", "TRG", "TRD")
        for suffix in ("cdr3", "v", "d", "j", "cdr3_nt", "clone_id", "umis", "reads")
    ]
    required_obs = {
        "sample_id",
        "library_id",
        "donor_id",
        "tissue",
        "diagnosis",
        "barcode_core",
        *tcr_columns,
        *CANONICAL_TCR_FLAGS,
    }
    try:
        with h5py.File(path, "r") as handle:
            matrix = handle.get("X")
            encoding = ""
            if isinstance(matrix, h5py.Group):
                raw_encoding = matrix.attrs.get("encoding-type", "")
                encoding = raw_encoding.decode() if isinstance(raw_encoding, bytes) else str(raw_encoding)
            if encoding not in {"csr_matrix", "csc_matrix"} or "data" not in matrix:
                issues.append(
                    ValidationIssue(
                        "error", "h5ad_sparse_x", "serialized X is not CSR/CSC sparse", cohort_id
                    )
                )
            else:
                data = matrix["data"]
                for start in range(0, len(data), 1_000_000):
                    values = data[start : start + 1_000_000]
                    if (
                        not np.isfinite(values).all()
                        or np.any(values < 0)
                        or (not np.issubdtype(values.dtype, np.integer) and not np.equal(values, np.floor(values)).all())
                    ):
                        issues.append(
                            ValidationIssue(
                                "error",
                                "h5ad_integer_x",
                                "serialized sparse X is not finite nonnegative integer counts",
                                cohort_id,
                            )
                        )
                        break

        adata = ad.read_h5ad(path, backed="r")
        obs = adata.obs
        if not obs.index.is_unique:
            issues.append(ValidationIssue("error", "h5ad_obs_names", "obs_names are not unique", cohort_id))
        missing = sorted(required_obs - set(obs.columns))
        if missing:
            issues.append(
                ValidationIssue(
                    "error", "h5ad_schema", f"missing obs columns: {missing}", cohort_id
                )
            )
        else:
            for column in ("sample_id", "library_id", "donor_id", "tissue", "diagnosis", "barcode_core"):
                values = pd.Series(obs[column], index=obs.index).astype("string").fillna("").str.strip()
                if values.eq("").any():
                    issues.append(
                        ValidationIssue(
                            "error", "h5ad_metadata", f"blank {column} values", cohort_id
                        )
                    )
            keys = obs[["sample_id", "barcode_core"]].astype(str)
            if keys.duplicated(keep=False).any():
                issues.append(
                    ValidationIssue(
                        "error",
                        "h5ad_join_keys",
                        "duplicate sample_id+barcode_core keys",
                        cohort_id,
                    )
                )
            for chain in ("TRA", "TRB", "TRG", "TRD"):
                expected = pd.Series(obs[f"{chain}_cdr3"], index=obs.index).astype("string").fillna("").str.strip().ne("")
                observed = pd.Series(obs[f"has_{chain}"], index=obs.index).astype(bool)
                if not np.array_equal(
                    expected.to_numpy(dtype=bool), observed.to_numpy(dtype=bool)
                ):
                    issues.append(
                        ValidationIssue(
                            "error",
                            "h5ad_tcr_flags",
                            f"has_{chain} is inconsistent with {chain}_cdr3",
                            cohort_id,
                        )
                    )
        if str(adata.uns.get("count_matrix_state", "")) != "raw_sparse_integer_counts":
            issues.append(
                ValidationIssue(
                    "error", "h5ad_count_provenance", "count_matrix_state provenance is missing", cohort_id
                )
            )
        if str(adata.uns.get("tcr_join_rule", "")) != "sample_id+barcode_core":
            issues.append(
                ValidationIssue(
                    "error", "h5ad_join_rule", "TCR join rule is not sample_id+barcode_core", cohort_id
                )
            )
        if cohort_id == "GSE294273_GSE294274" and (
            not pd.Series(obs["tissue"]).astype(str).eq("lymph_node_metastasis").all()
            or not pd.Series(obs["diagnosis"]).astype(str).eq("melanoma").all()
        ):
            issues.append(
                ValidationIssue(
                    "error",
                    "gse294_h5ad_metadata",
                    "serialized GSE294273/GSE294274 metadata is not melanoma/lymph_node_metastasis",
                    cohort_id,
                )
            )
        if cohort_id == "GSE159251" and (
            "tcr_schema_provenance" not in obs
            or not pd.Series(obs["tcr_schema_provenance"]).astype(str).str.contains("partial_embedded").all()
        ):
            issues.append(
                ValidationIssue(
                    "error", "gse159251_provenance", "partial TCR provenance is missing", cohort_id
                )
            )
        if getattr(adata, "file", None) is not None:
            adata.file.close()
    except (OSError, KeyError, ValueError, TypeError) as exc:
        issues.append(ValidationIssue("error", "h5ad_read", str(exc), cohort_id))
    return issues


def validate_extension_manifests(
    cohort_path: Path = DEFAULT_COHORTS,
    library_path: Path = DEFAULT_LIBRARIES,
    canonical_registry_path: Path = DEFAULT_CANONICAL_REGISTRY,
    selected_cohorts: Sequence[str] | None = None,
    shared_metasheet_path: Path = DEFAULT_SHARED_METASHEET,
) -> tuple[ValidationReport, list[dict[str, str]], list[dict[str, str]]]:
    """Validate manifests and return the report plus selected rows."""
    issues: list[ValidationIssue] = []
    try:
        cohorts = read_manifest(cohort_path, COHORT_COLUMNS)
        libraries = read_manifest(library_path, LIBRARY_COLUMNS)
        canonical_ids = _canonical_identifiers(canonical_registry_path)
    except (OSError, ManifestError) as exc:
        report = ValidationReport(
            [ValidationIssue("error", "manifest_parse", str(exc))]
        )
        return report, [], []

    _validate_cohort_rows(cohorts, canonical_ids, issues)
    _validate_library_rows(cohorts, libraries, issues)
    _validate_shared_metasheet(cohorts, shared_metasheet_path, issues)

    selected = set(selected_cohorts or [])
    if selected:
        known = {row["cohort_id"] for row in cohorts}
        unknown = selected - known
        if unknown:
            issues.append(
                ValidationIssue(
                    "error", "unknown_cohort", f"unknown --cohort value(s): {sorted(unknown)}"
                )
            )
        cohorts = [row for row in cohorts if row["cohort_id"] in selected]
        libraries = [row for row in libraries if row["cohort_id"] in selected]

    report = ValidationReport(
        issues=issues,
        cohort_count=len(cohorts),
        included_cohort_count=sum(_is_included(row) for row in cohorts),
        library_count=len(libraries),
    )
    return report, cohorts, libraries


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohorts-manifest", type=Path, default=DEFAULT_COHORTS)
    parser.add_argument("--libraries-manifest", type=Path, default=DEFAULT_LIBRARIES)
    parser.add_argument("--canonical-registry", type=Path, default=DEFAULT_CANONICAL_REGISTRY)
    parser.add_argument("--shared-metasheet", type=Path, default=DEFAULT_SHARED_METASHEET)
    parser.add_argument("--cohort", action="append", default=[], help="Validate/select one cohort ID; repeatable.")
    parser.add_argument(
        "--h5ad-root",
        type=Path,
        help="Also require and validate <cohort_id>.h5ad files under this extension output root.",
    )
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    report, cohorts, _ = validate_extension_manifests(
        args.cohorts_manifest,
        args.libraries_manifest,
        args.canonical_registry,
        args.cohort,
        args.shared_metasheet,
    )
    if args.h5ad_root:
        root = assert_extension_output_root(args.h5ad_root)
        for cohort in cohorts:
            if cohort["build_enabled"] != "true":
                continue
            path = root / f"{cohort['cohort_id']}.h5ad"
            if not path.is_file():
                report.issues.append(
                    ValidationIssue(
                        "error", "missing_h5ad", f"built H5AD is missing: {path}", cohort["cohort_id"]
                    )
                )
                continue
            report.issues.extend(validate_built_h5ad(path, cohort))
            report.h5ad_count += 1
    if args.json:
        print(json.dumps(report.to_dict(), indent=2, sort_keys=True))
    else:
        state = "PASS" if report.passed else "FAIL"
        print(
            f"{state}: {report.cohort_count} cohorts, {report.library_count} libraries, "
            f"{len(report.errors)} errors, {len(report.warnings)} warnings"
        )
        for issue in report.issues:
            location = "/".join(value for value in (issue.cohort_id, issue.library_id) if value)
            suffix = f" [{location}]" if location else ""
            print(f"{issue.severity.upper()} {issue.code}{suffix}: {issue.message}")
    return 0 if report.passed else 1


if __name__ == "__main__":
    sys.exit(main())
