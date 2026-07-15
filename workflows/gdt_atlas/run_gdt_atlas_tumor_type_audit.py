#!/usr/bin/env python3
"""Export sample-aware tumor-type review tables for the predicted gdT atlas.

This audit is intentionally review-only. It never writes back to the atlas H5AD.
Cancer type is assigned only when the atlas row itself carries a tumor-like
`tissue`/`tissue_corrected` value covered by a local metadata rule. Normal,
adjacent-normal, blood, and other rows from the same GSE remain unassigned.
"""

from __future__ import annotations

# TNK_WORKFLOW_BOOTSTRAP
import sys as _tnk_sys
from pathlib import Path as _TnkPath

_TNK_PROJECT_ROOT = _TnkPath(__file__).resolve().parents[2]
for _tnk_path in (
    _TNK_PROJECT_ROOT,
    _TNK_PROJECT_ROOT / "src",
    _TNK_PROJECT_ROOT / "workflows" / "integration",
    _TNK_PROJECT_ROOT / "workflows" / "intake",
    _TNK_PROJECT_ROOT / "workflows" / "metadata",
    _TNK_PROJECT_ROOT / "workflows" / "analysis",
    _TNK_PROJECT_ROOT / "workflows" / "gdtai",
    _TNK_PROJECT_ROOT / "workflows" / "gdt_atlas",
    _TNK_PROJECT_ROOT / "workflows" / "reporting",
    _TNK_PROJECT_ROOT / "workflows" / "maintenance",
):
    _tnk_value = str(_tnk_path)
    if _tnk_value not in _tnk_sys.path:
        _tnk_sys.path.insert(0, _tnk_value)


import argparse
import json
from pathlib import Path
from typing import Iterable

import anndata as ad
import pandas as pd

DEFAULT_ATLAS = _TNK_PROJECT_ROOT / "Integrated_dataset/gdT_atlas/predicted_gdt_cell_atlas.h5ad"
DEFAULT_RULES = _TNK_PROJECT_ROOT / "configs" / "gdt_atlas" / "metadata_rules.json"
DEFAULT_TABLE_DIR = _TNK_PROJECT_ROOT / "Integrated_dataset/tables/gdT_atlas"
DEFAULT_LOG_DIR = _TNK_PROJECT_ROOT / "Integrated_dataset/logs/gdT_atlas"
TUMOR_LIKE_DEFAULTS = {"tumor", "liver_tumor", "liver tumor", "lymph_node_metastasis", "lymph node metastasis"}
NON_TUMOR_WORDS = ("blood", "normal", "adjacent", "juxta")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit broad tumor labels in the predicted gdT atlas.")
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument("--rules", type=Path, default=DEFAULT_RULES)
    parser.add_argument("--table-dir", type=Path, default=DEFAULT_TABLE_DIR)
    parser.add_argument("--log-dir", type=Path, default=DEFAULT_LOG_DIR)
    parser.add_argument("--tissue-column", default="tissue_corrected")
    parser.add_argument("--source-column", default="source_gse_id")
    return parser.parse_args()


def normalize_value(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() in {"nan", "none", "na", "n/a", "<na>"}:
        return ""
    return text


def normalize_label(value: object) -> str:
    return normalize_value(value).lower()


def existing_columns(obs: pd.DataFrame, candidates: Iterable[str]) -> list[str]:
    return [column for column in candidates if column in obs.columns]


def compact_value_counts(series: pd.Series, limit: int = 12) -> str:
    cleaned = pd.Series([normalize_value(value) for value in series], dtype="object")
    counts = cleaned.value_counts(dropna=False)
    counts = counts[counts > 0]
    if counts.empty:
        return ""
    parts = [f"{idx or '<blank>'}: {int(count)}" for idx, count in counts.head(limit).items()]
    if len(counts) > limit:
        parts.append(f"... {len(counts) - limit} more")
    return "; ".join(parts)


def unique_examples(series: pd.Series, limit: int = 8) -> str:
    seen: list[str] = []
    for raw_value in series:
        value = normalize_value(raw_value)
        if not value or value in seen:
            continue
        seen.append(value)
        if len(seen) >= limit:
            break
    return "; ".join(seen)


def read_obs(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    adata = ad.read_h5ad(path, backed="r")
    try:
        return adata.obs.copy()
    finally:
        adata.file.close()


def build_local_sample_maps(root: Path) -> dict[str, dict[str, dict[str, str]]]:
    """Build optional sample-level maps from downloaded local metadata."""
    maps: dict[str, dict[str, dict[str, str]]] = {}

    gse243013 = read_obs(root / "downloads/per_gse_h5ad_with_metadata/GSE243013_TNK.h5ad")
    if gse243013 is not None and {"sampleID", "cancer_type"}.issubset(gse243013.columns):
        sample_map: dict[str, dict[str, str]] = {}
        gse243013 = gse243013[["sampleID", "cancer_type"]].copy()
        index_samples = pd.Series(
            [normalize_value(str(value).split("-", 1)[0]) for value in gse243013.index],
            index=gse243013.index,
            dtype="object",
        )
        explicit_samples = pd.Series(
            [normalize_value(value) for value in gse243013["sampleID"]],
            index=gse243013.index,
            dtype="object",
        )
        gse243013["_sample_for_map"] = explicit_samples.where(explicit_samples != "", index_samples)
        grouped = gse243013[["_sample_for_map", "cancer_type"]].drop_duplicates()
        for sample_id, sample_df in grouped.groupby("_sample_for_map", observed=True):
            sample_id = normalize_value(sample_id)
            if not sample_id:
                continue
            values = sorted({normalize_value(v) for v in sample_df["cancer_type"] if normalize_value(v)})
            if values:
                sample_map[sample_id] = {
                    "sample_tumor_subtype": "+".join(values),
                    "sample_tumor_subtype_full": "+".join(
                        {"LUAD": "lung adenocarcinoma", "LUSC": "lung squamous cell carcinoma"}.get(v, v) for v in values
                    ),
                    "sample_level_evidence": "downloads/per_gse_h5ad_with_metadata/GSE243013_TNK.h5ad: sampleID or cell-name prefix -> cancer_type",
                }
        maps["GSE243013"] = sample_map

    gse162498 = read_obs(root / "analysis_26GSE_V4/scanpy_projects/GSE162498/outputs/GSE162498_scanpy_processed.h5ad")
    if gse162498 is not None and {"sample", "sample_info_sample_type"}.issubset(gse162498.columns):
        sample_map = {}
        grouped = gse162498[["sample", "sample_info_sample_type"]].dropna().drop_duplicates()
        for sample_id, sample_df in grouped.groupby("sample", observed=True):
            values = sorted({normalize_value(v) for v in sample_df["sample_info_sample_type"] if normalize_value(v)})
            if values:
                sample_map[normalize_value(sample_id)] = {
                    "sample_source_type": "+".join(values),
                    "sample_level_evidence": "analysis_26GSE_V4/scanpy_projects/GSE162498/outputs/GSE162498_scanpy_processed.h5ad: sample -> sample_info_sample_type",
                }
        maps["GSE162498"] = sample_map

    gse235863 = read_obs(root / "downloads/per_gse_h5ad_with_metadata/GSE235863_with_tcr.h5ad")
    if gse235863 is not None and {"sample_id", "tissue"}.issubset(gse235863.columns):
        sample_map = {}
        grouped = gse235863[["sample_id", "tissue"]].dropna().drop_duplicates()
        for sample_id, sample_df in grouped.groupby("sample_id", observed=True):
            values = sorted({normalize_value(v) for v in sample_df["tissue"] if normalize_value(v)})
            if values:
                sample_map[normalize_value(sample_id)] = {
                    "sample_source_type": "+".join(values),
                    "sample_level_evidence": "downloads/per_gse_h5ad_with_metadata/GSE235863_with_tcr.h5ad: sample_id -> tissue",
                }
        maps["GSE235863"] = sample_map

    gse287301 = read_obs(root / "downloads/per_gse_h5ad_with_metadata/GSE287301_with_tcr.h5ad")
    if gse287301 is not None and {"sample_id", "condition", "sample_type"}.issubset(gse287301.columns):
        sample_map = {}
        grouped = gse287301[["sample_id", "condition", "sample_type", "donor_patient"]].drop_duplicates()
        for sample_id, sample_df in grouped.groupby("sample_id", observed=True):
            condition = "+".join(sorted({normalize_value(v) for v in sample_df["condition"] if normalize_value(v)}))
            sample_type = "+".join(sorted({normalize_value(v) for v in sample_df["sample_type"] if normalize_value(v)}))
            donor_patient = "+".join(sorted({normalize_value(v) for v in sample_df["donor_patient"] if normalize_value(v)}))
            sample_map[normalize_value(sample_id)] = {
                "sample_source_type": sample_type,
                "sample_condition": condition,
                "sample_patient_pool": donor_patient,
                "sample_level_evidence": "downloads/per_gse_h5ad_with_metadata/GSE287301_with_tcr.h5ad: sample_id -> condition/sample_type/donor_patient",
            }
        maps["GSE287301"] = sample_map

    return maps


def source_rule_row(source: str, rules: dict[str, dict[str, object]]) -> dict[str, object]:
    rule = rules.get(source, {})
    if not rule:
        return {
            "rule_status": "unresolved",
            "tumor_type": "unknown",
            "tumor_type_full": "unknown",
            "tumor_organ": "unknown",
            "confidence": "unresolved",
            "evidence_source": "none",
            "evidence_url": "",
            "evidence_summary": "No local source/sample tumor-type rule found.",
            "notes": "Requires manual review.",
        }
    return {
        "rule_status": "mapped",
        "tumor_type": rule.get("tumor_type", "unknown"),
        "tumor_type_full": rule.get("tumor_type_full", "unknown"),
        "tumor_organ": rule.get("tumor_organ", "unknown"),
        "confidence": rule.get("confidence", "unknown"),
        "evidence_source": rule.get("evidence_source", "unknown"),
        "evidence_url": rule.get("evidence_url", ""),
        "evidence_summary": rule.get("evidence_summary", ""),
        "notes": rule.get("notes", ""),
    }


def applies_to_row(source: str, tissue: object, rules: dict[str, dict[str, object]]) -> bool:
    tissue_norm = normalize_label(tissue)
    rule = rules.get(source)
    if not rule:
        return False
    allowed = {normalize_label(v) for v in rule.get("applies_to_tissue_values", [])}
    return tissue_norm in allowed


def is_tumor_like_label(tissue: object, rules: dict[str, dict[str, object]]) -> bool:
    tissue_norm = normalize_label(tissue)
    if tissue_norm in TUMOR_LIKE_DEFAULTS:
        return True
    for rule in rules.values():
        allowed = {normalize_label(v) for v in rule.get("applies_to_tissue_values", [])}
        if tissue_norm in allowed:
            return True
    return False


def proposed_assignment(row: pd.Series, rules: dict[str, dict[str, object]]) -> dict[str, object]:
    source = row["_audit_source"]
    tissue = row["_audit_tissue_raw"]
    if not applies_to_row(source, tissue, rules):
        tissue_norm = normalize_label(tissue)
        reason = "not_tumor_tissue"
        if tissue_norm == "":
            reason = "blank_or_unknown_tissue"
        elif any(word in tissue_norm for word in NON_TUMOR_WORDS):
            reason = "explicit_non_tumor_tissue"
        return {
            "proposed_tumor_context": "not_assigned_non_tumor_or_unmatched",
            "proposed_tumor_type": "not_assigned",
            "proposed_tumor_type_full": "not_assigned",
            "proposed_tumor_organ": "not_assigned",
            "assignment_reason": reason,
            "assignment_confidence": "not_applicable",
            "assignment_evidence": "atlas row tissue is not in applies_to_tissue_values for this source",
        }
    rule = source_rule_row(source, rules)
    return {
        "proposed_tumor_context": "tumor_or_metastasis",
        "proposed_tumor_type": rule["tumor_type"],
        "proposed_tumor_type_full": rule["tumor_type_full"],
        "proposed_tumor_organ": rule["tumor_organ"],
        "assignment_reason": "atlas row tissue matched source-specific tumor-like value",
        "assignment_confidence": rule["confidence"],
        "assignment_evidence": rule["evidence_source"],
    }


def markdown_table(df: pd.DataFrame, columns: list[str]) -> str:
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join(["---"] * len(columns)) + " |"
    lines = [header, sep]
    for _, row in df[columns].iterrows():
        values = []
        for column in columns:
            text = normalize_value(row.get(column, "")).replace("|", "/")
            values.append(text)
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    root = Path.cwd()
    args.table_dir.mkdir(parents=True, exist_ok=True)
    args.log_dir.mkdir(parents=True, exist_ok=True)

    with args.rules.open() as handle:
        all_rules = json.load(handle)
    source_rules = all_rules.get("tumor_type_source_rules", {})
    sample_maps = build_local_sample_maps(root)

    adata = ad.read_h5ad(args.atlas, backed="r")
    try:
        obs = adata.obs.copy()
    finally:
        adata.file.close()

    if args.source_column not in obs.columns:
        raise KeyError(f"Required source column missing from atlas obs: {args.source_column}")
    if args.tissue_column not in obs.columns:
        fallback = "tissue" if "tissue" in obs.columns else None
        if fallback is None:
            raise KeyError(f"Required tissue column missing from atlas obs: {args.tissue_column}")
        args.tissue_column = fallback

    obs["_audit_source"] = obs[args.source_column].map(normalize_value)
    obs["_audit_tissue_raw"] = obs[args.tissue_column].map(normalize_value)
    obs["_audit_tissue_norm"] = obs[args.tissue_column].map(normalize_label)
    obs["_audit_sample"] = obs["sample_id"].map(normalize_value) if "sample_id" in obs.columns else ""
    obs["_is_tumor_like_label"] = obs["_audit_tissue_raw"].map(lambda value: is_tumor_like_label(value, source_rules))

    tumor_like_sources = sorted(obs.loc[obs["_is_tumor_like_label"], "_audit_source"].dropna().unique())
    rule_sources = sorted(source_rules.keys())
    review_sources = sorted(set(tumor_like_sources) | set(rule_sources))
    review_obs = obs.loc[obs["_audit_source"].isin(review_sources)].copy()

    assignments = review_obs.apply(lambda row: pd.Series(proposed_assignment(row, source_rules)), axis=1)
    review_obs = pd.concat([review_obs, assignments], axis=1)

    for field in ["sample_source_type", "sample_condition", "sample_patient_pool", "sample_tumor_subtype", "sample_tumor_subtype_full", "sample_level_evidence"]:
        review_obs[field] = ""
    for source, sample_map in sample_maps.items():
        mask = review_obs["_audit_source"] == source
        if not mask.any():
            continue
        for field in ["sample_source_type", "sample_condition", "sample_patient_pool", "sample_tumor_subtype", "sample_tumor_subtype_full", "sample_level_evidence"]:
            review_obs.loc[mask, field] = review_obs.loc[mask, "_audit_sample"].map(
                lambda sample: sample_map.get(sample, {}).get(field, "")
            )

    subtype_mask = (review_obs["_audit_source"] == "GSE243013") & (review_obs["proposed_tumor_context"] == "tumor_or_metastasis") & (review_obs["sample_tumor_subtype"] != "")
    review_obs.loc[subtype_mask, "proposed_tumor_type"] = review_obs.loc[subtype_mask, "sample_tumor_subtype"]
    review_obs.loc[subtype_mask, "proposed_tumor_type_full"] = review_obs.loc[subtype_mask, "sample_tumor_subtype_full"]
    review_obs.loc[subtype_mask, "assignment_reason"] = "atlas tumor row plus local sampleID cancer_type match"
    review_obs.loc[subtype_mask, "assignment_evidence"] = review_obs.loc[subtype_mask, "sample_level_evidence"]

    summary_rows: list[dict[str, object]] = []
    for (source, tissue), group in review_obs.groupby(["_audit_source", "_audit_tissue_raw"], observed=True, dropna=False):
        rule = source_rule_row(source, source_rules)
        assigned = group["proposed_tumor_context"].eq("tumor_or_metastasis")
        sample_count = int(group["_audit_sample"].loc[group["_audit_sample"] != ""].nunique())
        row = {
            "source_gse_id": source,
            "atlas_tissue_value": tissue or "<blank>",
            "cell_count": int(group.shape[0]),
            "sample_count": sample_count,
            "sample_examples": unique_examples(group["_audit_sample"]),
            "assigned_cells": int(assigned.sum()),
            "not_assigned_cells": int((~assigned).sum()),
            "proposed_tumor_type_values": compact_value_counts(group["proposed_tumor_type"]),
            "sample_source_type_values": compact_value_counts(group["sample_source_type"]),
            "sample_condition_values": compact_value_counts(group["sample_condition"]),
            "sample_tumor_subtype_values": compact_value_counts(group["sample_tumor_subtype"]),
            "rule_status": rule["rule_status"],
            "evidence_source": rule["evidence_source"],
            "evidence_summary": rule["evidence_summary"],
            "notes": rule["notes"],
        }
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows).sort_values(["source_gse_id", "assigned_cells", "cell_count"], ascending=[True, False, False])

    source_rows: list[dict[str, object]] = []
    for source, group in review_obs.groupby("_audit_source", observed=True, dropna=False):
        rule = source_rule_row(source, source_rules)
        assigned = group["proposed_tumor_context"].eq("tumor_or_metastasis")
        source_rows.append(
            {
                "source_gse_id": source,
                "atlas_cells_reviewed": int(group.shape[0]),
                "assigned_tumor_or_metastasis_cells": int(assigned.sum()),
                "not_assigned_cells_in_same_source": int((~assigned).sum()),
                "atlas_tissue_values": compact_value_counts(group[args.tissue_column]),
                "proposed_tumor_type_values": compact_value_counts(group.loc[assigned, "proposed_tumor_type"]),
                "rule_status": rule["rule_status"],
                "evidence_source": rule["evidence_source"],
                "evidence_summary": rule["evidence_summary"],
                "notes": rule["notes"],
            }
        )
    source_summary = pd.DataFrame(source_rows).sort_values(["assigned_tumor_or_metastasis_cells", "source_gse_id"], ascending=[False, True])

    export_cols = existing_columns(
        review_obs,
        [
            args.source_column,
            "sample_id",
            "sampleid",
            "library_id",
            "tissue",
            "tissue_corrected",
            "condition",
            "sample_source_type",
            "sample_condition",
            "sample_patient_pool",
            "sample_tumor_subtype",
            "sample_tumor_subtype_full",
            "proposed_tumor_context",
            "proposed_tumor_type",
            "proposed_tumor_type_full",
            "proposed_tumor_organ",
            "assignment_reason",
            "assignment_confidence",
            "assignment_evidence",
            "sample_level_evidence",
        ],
    )
    cell_review = review_obs[export_cols].copy()

    source_csv = args.table_dir / "tumor_type_source_mapping.csv"
    tissue_csv = args.table_dir / "tumor_type_by_source_tissue.csv"
    cell_csv = args.table_dir / "tumor_type_cell_annotation_review.csv.gz"
    source_summary.to_csv(source_csv, index=False)
    summary.to_csv(tissue_csv, index=False)
    cell_review.to_csv(cell_csv, index=True, index_label="cell_id", compression="gzip")

    assigned_total = int(review_obs["proposed_tumor_context"].eq("tumor_or_metastasis").sum())
    not_assigned_total = int(review_obs.shape[0] - assigned_total)
    md_path = args.log_dir / "tumor_type_source_mapping.md"
    md_lines = [
        "# Predicted gdT Atlas Tumor-Type Audit",
        "",
        "This review uses downloaded local metadata in the working directory. It does not use web lookup and does not write back to the atlas H5AD.",
        "",
        "## Main Conclusion",
        "",
        f"- Reviewed {int(review_obs.shape[0]):,} atlas cells from sources with tumor-like labels or tumor-type rules.",
        f"- Assigned tumor type to {assigned_total:,} cells whose own atlas tissue/sample context is tumor-like.",
        f"- Left {not_assigned_total:,} cells unassigned because they are blood, normal, adjacent normal, blank/unknown, or otherwise not a tumor-like row in the same accession.",
        "- Therefore, mixed accessions such as GSE162498, GSE206325, and GSE235863 are not globally labeled as tumor; only their tumor rows receive cancer type.",
        "",
        "## Source-Level Summary",
        "",
        markdown_table(
            source_summary,
            [
                "source_gse_id",
                "atlas_cells_reviewed",
                "assigned_tumor_or_metastasis_cells",
                "not_assigned_cells_in_same_source",
                "atlas_tissue_values",
                "proposed_tumor_type_values",
                "rule_status",
            ],
        ),
        "",
        "## Tissue-Level Summary",
        "",
        markdown_table(
            summary,
            [
                "source_gse_id",
                "atlas_tissue_value",
                "cell_count",
                "assigned_cells",
                "not_assigned_cells",
                "proposed_tumor_type_values",
                "sample_source_type_values",
                "sample_tumor_subtype_values",
            ],
        ),
        "",
        "## Interpretation Notes",
        "",
        "- GSE243013 is all tumor in the promoted atlas, and local obs resolves its tumor rows to LUAD or LUSC by sampleID where available.",
        "- GSE162498 contains NSCLC tumor, blood, and adjacent-normal/juxta rows. Only the `tumor` rows are assigned NSCLC.",
        "- GSE206325 contains HCC Tumor, Normal, and blank/unknown rows. Only `Tumor` rows are assigned HCC.",
        "- GSE235863 contains HBV-positive HCC liver tumor, blood, and adjacent normal liver. Only `liver_tumor` rows are assigned HBV-positive HCC.",
        "- GSE287301 atlas rows are tumor pools with local condition HNSCC; one original patient-level caveat remains because atlas rows are pool-level.",
        "- GSE190870 is included because its atlas tissue is `lymph_node_metastasis`; it is breast invasive ductal carcinoma by local metadata.",
        "",
        "## Outputs",
        "",
        f"- `{source_csv}`",
        f"- `{tissue_csv}`",
        f"- `{cell_csv}`",
    ]
    md_path.write_text("\n".join(md_lines) + "\n")

    print(f"Reviewed sources: {', '.join(review_sources)}")
    print(f"Assigned tumor/metastasis rows: {assigned_total}")
    print(f"Not assigned rows in same sources: {not_assigned_total}")
    print(f"Wrote {source_csv}")
    print(f"Wrote {tissue_csv}")
    print(f"Wrote {cell_csv}")
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
