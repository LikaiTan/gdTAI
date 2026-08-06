#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Validate the bounded official-GEO reconciliation for extension cohorts.

The workflow opens extension H5AD files in backed read-only mode. It never
changes H5AD content or harmonized metadata; it only writes audit tables and a
log beneath the delegated GEO reconciliation output roots.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
from typing import Iterable
from urllib.parse import urlparse

import anndata as ad
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = (
    PROJECT_ROOT / "configs/metadata/extension_geo_metadata_reconciliation.csv"
)
COHORT_CONFIG_PATH = PROJECT_ROOT / "configs/datasets/extension_cohorts.csv"
LIBRARY_CONFIG_PATH = PROJECT_ROOT / "configs/datasets/extension_libraries.csv"
FILTERED_MANIFEST_PATH = (
    PROJECT_ROOT / "data/interim/extension_intake/tnk_filtered_h5ads_manifest.csv"
)
STAGED_ROOT = PROJECT_ROOT / "data/interim/extension_intake/staged"
TABLE_ROOT = PROJECT_ROOT / "Integrated_dataset/tables/geo_metadata_reconciliation"
LOG_ROOT = PROJECT_ROOT / "Integrated_dataset/logs/geo_metadata_reconciliation"

COHORT_ACCESSIONS = {
    "GSE114724": {"GSE114724": "combined_gex_vdj"},
    "GSE121636_GSE121637": {
        "GSE121636": "gex",
        "GSE121637": "vdj",
    },
    "GSE159251": {"GSE159251": "combined_gex_vdj_records"},
    "GSE169246": {"GSE169246": "combined_gex_vdj_plus_separate_atac"},
    "GSE292700": {"GSE292700": "gex_and_vdj_records"},
    "GSE294273_GSE294274": {
        "GSE294273": "gex",
        "GSE294274": "vdj",
    },
    "GSE296954": {"GSE296954": "gex_and_vdj_records"},
    "GSE315928": {"GSE315928": "combined_gex_vdj_records"},
}

CELL_LEVEL_ACCESSIONS = {
    "GSE114724": {"GSE114724"},
    "GSE121636_GSE121637": {"GSE121636"},
    "GSE159251": {"GSE159251"},
    "GSE169246": {"GSE169246"},
    "GSE292700": {"GSE292700"},
    "GSE294273_GSE294274": {"GSE294273"},
    "GSE296954": {"GSE296954"},
    "GSE315928": {"GSE315928"},
}

STAGED_EXPECTED = set(COHORT_ACCESSIONS) - {"GSE169246"}
ALLOWED_STATUS = {
    "resolved_verified_local",
    "resolved_geo_value_local_unresolved",
    "ambiguous_geo_partial",
    "unavailable_in_geo",
}
ALLOWED_CONFIDENCE = {"high", "medium", "low"}
ALLOWED_GEO_HOSTS = {"ftp.ncbi.nlm.nih.gov", "www.ncbi.nlm.nih.gov"}
REQUIRED_COLUMNS = (
    "cohort_id",
    "accession",
    "scope",
    "field",
    "local_column",
    "local_scope_regex",
    "local_value",
    "geo_supported_value",
    "status",
    "confidence",
    "direct_url",
    "evidence_note",
)
INVENTORY_FIELDS = (
    "sample_id",
    "library_id",
    "donor_id",
    "patient_id",
    "tissue",
    "tissue_harmonized",
    "specimen_context",
    "diagnosis",
    "source_accession",
    "source_gse_id",
    "treatment",
    "timepoint",
    "timepoint_group",
)


class AuditError(RuntimeError):
    """Raised when the reconciliation contract is not satisfied."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Validate without replacing generated audit outputs.",
    )
    return parser.parse_args()


def clean_series(series: pd.Series) -> pd.Series:
    return (
        series.astype("string")
        .fillna("")
        .str.strip()
        .replace({"nan": "", "None": "", "<NA>": ""})
    )


def display_value(value: str) -> str:
    return value if value else "<blank>"


def split_expected(value: str) -> set[str]:
    return {"" if item == "<blank>" else item for item in value.split("|")}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_join(values: Iterable[str]) -> str:
    return "|".join(sorted({display_value(value) for value in values}))


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    config = pd.read_csv(CONFIG_PATH, dtype=str, keep_default_na=False)
    cohorts = pd.read_csv(COHORT_CONFIG_PATH, dtype=str, keep_default_na=False)
    libraries = pd.read_csv(LIBRARY_CONFIG_PATH, dtype=str, keep_default_na=False)
    manifest = pd.read_csv(FILTERED_MANIFEST_PATH, dtype=str, keep_default_na=False)
    return config, cohorts, libraries, manifest


def validate_config_schema(config: pd.DataFrame) -> list[dict[str, str]]:
    issues: list[dict[str, str]] = []
    if tuple(config.columns) != REQUIRED_COLUMNS:
        issues.append(
            {
                "cohort_id": "",
                "accession": "",
                "field": "",
                "issue": "schema_mismatch",
                "detail": f"expected={list(REQUIRED_COLUMNS)} actual={list(config.columns)}",
            }
        )
        return issues

    expected_cohorts = set(COHORT_ACCESSIONS)
    actual_cohorts = set(config["cohort_id"])
    if actual_cohorts != expected_cohorts:
        issues.append(
            {
                "cohort_id": "",
                "accession": "",
                "field": "",
                "issue": "cohort_coverage_mismatch",
                "detail": f"expected={sorted(expected_cohorts)} actual={sorted(actual_cohorts)}",
            }
        )

    expected_accessions = {
        accession
        for accessions in COHORT_ACCESSIONS.values()
        for accession in accessions
    }
    actual_accessions = set(config["accession"])
    if actual_accessions != expected_accessions:
        issues.append(
            {
                "cohort_id": "",
                "accession": "",
                "field": "",
                "issue": "accession_coverage_mismatch",
                "detail": f"expected={sorted(expected_accessions)} actual={sorted(actual_accessions)}",
            }
        )

    duplicate = config.duplicated(["cohort_id", "accession", "scope", "field"])
    for row in config.loc[duplicate].itertuples(index=False):
        issues.append(
            {
                "cohort_id": row.cohort_id,
                "accession": row.accession,
                "field": row.field,
                "issue": "duplicate_reconciliation_key",
                "detail": row.scope,
            }
        )

    for row in config.itertuples(index=False):
        if row.cohort_id not in COHORT_ACCESSIONS:
            continue
        if row.accession not in COHORT_ACCESSIONS[row.cohort_id]:
            issues.append(
                {
                    "cohort_id": row.cohort_id,
                    "accession": row.accession,
                    "field": row.field,
                    "issue": "accession_not_owned_by_cohort",
                    "detail": "",
                }
            )
        if row.status not in ALLOWED_STATUS:
            issues.append(
                {
                    "cohort_id": row.cohort_id,
                    "accession": row.accession,
                    "field": row.field,
                    "issue": "invalid_status",
                    "detail": row.status,
                }
            )
        if row.confidence not in ALLOWED_CONFIDENCE:
            issues.append(
                {
                    "cohort_id": row.cohort_id,
                    "accession": row.accession,
                    "field": row.field,
                    "issue": "invalid_confidence",
                    "detail": row.confidence,
                }
            )
        parsed_url = urlparse(row.direct_url)
        if (
            parsed_url.scheme != "https"
            or parsed_url.hostname not in ALLOWED_GEO_HOSTS
        ):
            issues.append(
                {
                    "cohort_id": row.cohort_id,
                    "accession": row.accession,
                    "field": row.field,
                    "issue": "non_official_geo_direct_url",
                    "detail": row.direct_url,
                }
            )
        if not row.evidence_note.strip():
            issues.append(
                {
                    "cohort_id": row.cohort_id,
                    "accession": row.accession,
                    "field": row.field,
                    "issue": "blank_evidence_note",
                    "detail": "",
                }
            )
        if row.status == "unavailable_in_geo" and row.geo_supported_value:
            issues.append(
                {
                    "cohort_id": row.cohort_id,
                    "accession": row.accession,
                    "field": row.field,
                    "issue": "geo_unavailable_has_value",
                    "detail": row.geo_supported_value,
                }
            )
        if (
            row.status != "unavailable_in_geo"
            and not row.geo_supported_value.strip()
        ):
            issues.append(
                {
                    "cohort_id": row.cohort_id,
                    "accession": row.accession,
                    "field": row.field,
                    "issue": "geo_supported_status_has_no_value",
                    "detail": row.status,
                }
            )
    return issues


def manifest_paths(manifest: pd.DataFrame) -> dict[str, Path]:
    paths = {
        row.cohort_id: Path(row.output_h5ad)
        for row in manifest.itertuples(index=False)
    }
    if set(paths) != set(COHORT_ACCESSIONS):
        raise AuditError(
            "Filtered manifest cohort coverage differs from reconciliation scope"
        )
    return paths


def summarize_local_h5ads(
    config: pd.DataFrame, paths: dict[str, Path]
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    list[dict[str, str]],
    dict[str, set[str]],
]:
    inventory_rows: list[dict[str, object]] = []
    scoped_value_rows: list[dict[str, object]] = []
    read_only_rows: list[dict[str, object]] = []
    issues: list[dict[str, str]] = []
    observed_accessions: dict[str, set[str]] = {}

    for cohort_id in COHORT_ACCESSIONS:
        path = paths[cohort_id]
        if not path.is_file():
            raise FileNotFoundError(path)
        before = path.stat()
        adata = ad.read_h5ad(path, backed="r")
        try:
            obs = adata.obs
            for field in INVENTORY_FIELDS:
                if field not in obs:
                    continue
                values = clean_series(obs[field])
                counts = values.value_counts(sort=False).sort_index()
                for value, count in counts.items():
                    inventory_rows.append(
                        {
                            "cohort_id": cohort_id,
                            "field": field,
                            "value": display_value(str(value)),
                            "cell_count": int(count),
                        }
                    )

            source_field = (
                "source_gse_id"
                if "source_gse_id" in obs
                else "source_accession"
            )
            observed_accessions[cohort_id] = {
                value for value in clean_series(obs[source_field]) if value
            }

            cohort_rows = config.loc[config["cohort_id"].eq(cohort_id)]
            for row in cohort_rows.itertuples(index=False):
                if not row.local_column:
                    continue
                if row.local_column not in obs:
                    issues.append(
                        {
                            "cohort_id": cohort_id,
                            "accession": row.accession,
                            "field": row.field,
                            "issue": "local_column_missing",
                            "detail": row.local_column,
                        }
                    )
                    continue
                selected = pd.Series(True, index=obs.index)
                if row.local_scope_regex:
                    if "library_id" not in obs:
                        issues.append(
                            {
                                "cohort_id": cohort_id,
                                "accession": row.accession,
                                "field": row.field,
                                "issue": "scope_requires_library_id",
                                "detail": row.local_scope_regex,
                            }
                        )
                        continue
                    selected = clean_series(obs["library_id"]).str.match(
                        row.local_scope_regex, na=False
                    )
                    if not selected.any():
                        issues.append(
                            {
                                "cohort_id": cohort_id,
                                "accession": row.accession,
                                "field": row.field,
                                "issue": "local_scope_matches_no_cells",
                                "detail": row.local_scope_regex,
                            }
                        )
                        continue
                scoped_values = clean_series(obs.loc[selected, row.local_column])
                actual = set(scoped_values)
                for value, count in (
                    scoped_values.value_counts(sort=False).sort_index().items()
                ):
                    scoped_value_rows.append(
                        {
                            "cohort_id": cohort_id,
                            "accession": row.accession,
                            "scope": row.scope,
                            "field": row.field,
                            "local_column": row.local_column,
                            "local_scope_regex": row.local_scope_regex,
                            "observed_value": display_value(str(value)),
                            "cell_count": int(count),
                        }
                    )
                expected = split_expected(row.local_value)
                if actual != expected:
                    issues.append(
                        {
                            "cohort_id": cohort_id,
                            "accession": row.accession,
                            "field": row.field,
                            "issue": "local_value_mismatch",
                            "detail": (
                                f"column={row.local_column} scope={row.local_scope_regex or 'all'} "
                                f"expected={stable_join(expected)} actual={stable_join(actual)}"
                            ),
                        }
                    )
        finally:
            adata.file.close()
        after = path.stat()
        read_only_rows.append(
            {
                "cohort_id": cohort_id,
                "path": str(path),
                "size_before": before.st_size,
                "size_after": after.st_size,
                "mtime_ns_before": before.st_mtime_ns,
                "mtime_ns_after": after.st_mtime_ns,
                "unchanged_during_audit": (
                    before.st_size == after.st_size
                    and before.st_mtime_ns == after.st_mtime_ns
                ),
            }
        )

    inventory = pd.DataFrame(inventory_rows).sort_values(
        ["cohort_id", "field", "value"], kind="stable"
    )
    scoped_values = pd.DataFrame(scoped_value_rows).sort_values(
        ["cohort_id", "accession", "scope", "field", "observed_value"],
        kind="stable",
    )
    read_only = pd.DataFrame(read_only_rows)
    return inventory, scoped_values, read_only, issues, observed_accessions


def build_accession_coverage(
    config: pd.DataFrame,
    libraries: pd.DataFrame,
    observed_accessions: dict[str, set[str]],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cohort_id, accessions in COHORT_ACCESSIONS.items():
        staged_path = STAGED_ROOT / cohort_id / "sample_metasheet.csv"
        if cohort_id in STAGED_EXPECTED:
            staged = pd.read_csv(staged_path, dtype=str, keep_default_na=False)
        else:
            staged = pd.DataFrame(columns=["gse"])
        library_subset = libraries.loc[libraries["cohort_id"].eq(cohort_id)]
        for accession, role in accessions.items():
            config_rows = int(config["accession"].eq(accession).sum())
            library_rows = int(
                library_subset["source_accession"].eq(accession).sum()
            )
            staged_rows = int(staged["gse"].eq(accession).sum())
            staged_expected = cohort_id in STAGED_EXPECTED
            cell_expected = accession in CELL_LEVEL_ACCESSIONS[cohort_id]
            cell_observed = accession in observed_accessions[cohort_id]
            coverage_pass = (
                config_rows > 0
                and library_rows > 0
                and (not staged_expected or staged_rows > 0)
                and (cell_observed == cell_expected)
            )
            rows.append(
                {
                    "cohort_id": cohort_id,
                    "accession": accession,
                    "role": role,
                    "reconciliation_rows": config_rows,
                    "library_config_rows": library_rows,
                    "staged_metasheet_expected": staged_expected,
                    "staged_metasheet_rows": staged_rows,
                    "cell_level_accession_expected": cell_expected,
                    "cell_level_accession_observed": cell_observed,
                    "coverage_pass": coverage_pass,
                }
            )
    return pd.DataFrame(rows).sort_values(["cohort_id", "accession"])


def status_summary(config: pd.DataFrame) -> pd.DataFrame:
    return (
        config.groupby(["cohort_id", "status"], sort=True)
        .size()
        .rename("reconciliation_rows")
        .reset_index()
    )


def render_log(
    config: pd.DataFrame,
    summary: pd.DataFrame,
    coverage: pd.DataFrame,
    read_only: pd.DataFrame,
) -> str:
    status_counts = config["status"].value_counts().sort_index()
    summary_lines = [
        "| cohort_id | status | reconciliation_rows |",
        "|---|---|---:|",
    ]
    summary_lines.extend(
        f"| {row.cohort_id} | {row.status} | {int(row.reconciliation_rows)} |"
        for row in summary.itertuples(index=False)
    )
    lines = [
        "# Extension GEO Metadata Reconciliation Audit",
        "",
        "- Result: **PASS**",
        "- Mode: read-only H5AD metadata inspection",
        f"- Reconciliation rows: {len(config)}",
        f"- Cohorts: {config['cohort_id'].nunique()}",
        f"- Official GEO accessions: {config['accession'].nunique()}",
        f"- Accession coverage rows passing: {int(coverage['coverage_pass'].sum())}/{len(coverage)}",
        f"- H5AD files unchanged during audit: {int(read_only['unchanged_during_audit'].sum())}/{len(read_only)}",
        "",
        "## Status counts",
        "",
    ]
    for status, count in status_counts.items():
        lines.append(f"- `{status}`: {int(count)}")
    lines.extend(
        [
            "",
            "## Cohort summary",
            "",
            *summary_lines,
            "",
            "## Outputs",
            "",
            "- `Integrated_dataset/tables/geo_metadata_reconciliation/validated_reconciliation.csv`",
            "- `Integrated_dataset/tables/geo_metadata_reconciliation/reconciliation_status_summary.csv`",
            "- `Integrated_dataset/tables/geo_metadata_reconciliation/source_accession_coverage.csv`",
            "- `Integrated_dataset/tables/geo_metadata_reconciliation/local_field_inventory.csv`",
            "- `Integrated_dataset/tables/geo_metadata_reconciliation/scoped_local_value_counts.csv`",
            "- `Integrated_dataset/tables/geo_metadata_reconciliation/unresolved_or_ambiguous.csv`",
            "- `Integrated_dataset/tables/geo_metadata_reconciliation/h5ad_read_only_validation.csv`",
            "",
        ]
    )
    return "\n".join(lines)


def write_outputs(
    config: pd.DataFrame,
    inventory: pd.DataFrame,
    scoped_values: pd.DataFrame,
    read_only: pd.DataFrame,
    coverage: pd.DataFrame,
    issues: pd.DataFrame,
) -> None:
    TABLE_ROOT.mkdir(parents=True, exist_ok=True)
    LOG_ROOT.mkdir(parents=True, exist_ok=True)
    ordered = config.sort_values(
        ["cohort_id", "accession", "scope", "field"], kind="stable"
    ).reset_index(drop=True)
    summary = status_summary(ordered)
    unresolved = ordered.loc[
        ordered["status"].isin(
            {"ambiguous_geo_partial", "unavailable_in_geo"}
        )
    ].copy()

    ordered.to_csv(TABLE_ROOT / "validated_reconciliation.csv", index=False)
    summary.to_csv(TABLE_ROOT / "reconciliation_status_summary.csv", index=False)
    coverage.to_csv(TABLE_ROOT / "source_accession_coverage.csv", index=False)
    inventory.to_csv(TABLE_ROOT / "local_field_inventory.csv", index=False)
    scoped_values.to_csv(
        TABLE_ROOT / "scoped_local_value_counts.csv", index=False
    )
    unresolved.to_csv(TABLE_ROOT / "unresolved_or_ambiguous.csv", index=False)
    read_only.to_csv(TABLE_ROOT / "h5ad_read_only_validation.csv", index=False)
    issues.to_csv(TABLE_ROOT / "validation_issues.csv", index=False)

    manifest = {
        "schema_version": 1,
        "mode": "read_only",
        "config": str(CONFIG_PATH.relative_to(PROJECT_ROOT)),
        "config_sha256": sha256_file(CONFIG_PATH),
        "cohorts": len(COHORT_ACCESSIONS),
        "official_geo_accessions": sum(
            len(accessions) for accessions in COHORT_ACCESSIONS.values()
        ),
        "reconciliation_rows": len(ordered),
        "scoped_local_value_rows": len(scoped_values),
        "ambiguous_rows": int(
            ordered["status"].eq("ambiguous_geo_partial").sum()
        ),
        "geo_unavailable_rows": int(
            ordered["status"].eq("unavailable_in_geo").sum()
        ),
        "validation_issues": len(issues),
        "accession_coverage_pass": bool(coverage["coverage_pass"].all()),
        "h5ads_unchanged": bool(read_only["unchanged_during_audit"].all()),
        "h5ads_modified": False,
        "input_sha256": {
            str(path.relative_to(PROJECT_ROOT)): sha256_file(path)
            for path in (
                COHORT_CONFIG_PATH,
                LIBRARY_CONFIG_PATH,
                FILTERED_MANIFEST_PATH,
            )
        },
    }
    (TABLE_ROOT / "audit_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (LOG_ROOT / "extension_geo_metadata_reconciliation_audit.md").write_text(
        render_log(ordered, summary, coverage, read_only), encoding="utf-8"
    )


def main() -> None:
    args = parse_args()
    config, cohorts, libraries, manifest = load_inputs()
    issues = validate_config_schema(config)

    scoped_cohorts = set(
        cohorts.loc[cohorts["cohort_id"].isin(COHORT_ACCESSIONS), "cohort_id"]
    )
    if scoped_cohorts != set(COHORT_ACCESSIONS):
        issues.append(
            {
                "cohort_id": "",
                "accession": "",
                "field": "",
                "issue": "extension_cohort_config_coverage_mismatch",
                "detail": stable_join(scoped_cohorts),
            }
        )

    paths = manifest_paths(manifest)
    inventory, scoped_values, read_only, local_issues, observed = (
        summarize_local_h5ads(config, paths)
    )
    issues.extend(local_issues)
    coverage = build_accession_coverage(config, libraries, observed)
    for row in coverage.loc[~coverage["coverage_pass"]].itertuples(index=False):
        issues.append(
            {
                "cohort_id": row.cohort_id,
                "accession": row.accession,
                "field": "source_accession",
                "issue": "source_accession_coverage_failed",
                "detail": "",
            }
        )
    if not read_only["unchanged_during_audit"].all():
        issues.append(
            {
                "cohort_id": "",
                "accession": "",
                "field": "",
                "issue": "h5ad_changed_during_audit",
                "detail": "",
            }
        )

    issue_frame = pd.DataFrame(
        issues,
        columns=["cohort_id", "accession", "field", "issue", "detail"],
    )
    if not issue_frame.empty:
        raise AuditError(
            f"Reconciliation validation failed with {len(issue_frame)} issue(s):\n"
            + issue_frame.to_string(index=False)
        )

    if not args.check_only:
        write_outputs(
            config,
            inventory,
            scoped_values,
            read_only,
            coverage,
            issue_frame,
        )
    counts = config["status"].value_counts().sort_index().to_dict()
    print(
        json.dumps(
            {
                "result": "PASS",
                "check_only": args.check_only,
                "cohorts": len(COHORT_ACCESSIONS),
                "accessions": sum(
                    len(accessions) for accessions in COHORT_ACCESSIONS.values()
                ),
                "rows": len(config),
                "status_counts": counts,
                "h5ads_modified": False,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
