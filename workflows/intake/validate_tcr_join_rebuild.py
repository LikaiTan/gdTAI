#!/usr/bin/env python3
"""Validate rebuilt TCR sidecars, staged replacement, and static report."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq


PROJECT_ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = PROJECT_ROOT / "Integrated_dataset/tables/tcr_join_rebuild"
LOG_DIR = PROJECT_ROOT / "Integrated_dataset/logs/tcr_join_rebuild"
SUMMARY = TABLE_DIR / "source_rebuild_summary.csv"
REPLACEMENT = TABLE_DIR / "validated_tcr_replacement_sidecar.parquet"
REPLACEMENT_MANIFEST = TABLE_DIR / "validated_tcr_replacement_manifest.csv"
RUN_MANIFEST = LOG_DIR / "run_manifest.json"
REPORT_DIR = PROJECT_ROOT / "gdT_prediction/gdtai_v4_2_tcr_join_rebuild"
OUTPUT_JSON = LOG_DIR / "final_validation.json"
OUTPUT_MD = LOG_DIR / "final_validation.md"
CHAINS = ("TRA", "TRB", "TRG", "TRD")


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def check(results: list[dict[str, object]], name: str, passed: bool, detail: str) -> None:
    results.append({"check": name, "passed": bool(passed), "detail": detail})


def main() -> None:
    summary = pd.read_csv(SUMMARY)
    replacement_manifest = pd.read_csv(REPLACEMENT_MANIFEST)
    run_manifest = json.loads(RUN_MANIFEST.read_text(encoding="utf-8"))
    results: list[dict[str, object]] = []

    check(results, "fourteen_sources", len(summary) == 14, f"observed={len(summary)}")
    n_pass = int(summary["rebuild_status"].str.startswith("PASS").sum())
    n_quarantine = int(summary["rebuild_status"].str.startswith("QUARANTINED").sum())
    check(results, "source_outcome_partition", n_pass == 11 and n_quarantine == 3, f"pass={n_pass}; quarantine={n_quarantine}")
    check(results, "recorded_h5ads_unchanged", summary["source_h5ad_unchanged"].all(), "all recorded size/mtime comparisons are true")

    sidecar_errors: list[str] = []
    current_stat_errors: list[str] = []
    for row in summary.itertuples(index=False):
        path = TABLE_DIR / "per_source" / f"{row.source_gse_id}.parquet"
        columns = [
            "source_obs_name",
            "join_status",
            "tcr_assignment_eligible",
            "replacement_eligible",
            "TCRseq_rebuilt",
            "has_TRA_TRB_paired_rebuilt",
            "has_TRG_TRD_paired_rebuilt",
            *[f"{chain}_cdr3" for chain in CHAINS],
            *[f"{chain}_umis" for chain in CHAINS],
            *[f"{chain}_umi_available" for chain in CHAINS],
        ]
        frame = pd.read_parquet(path, columns=columns)
        local: list[str] = []
        if len(frame) != int(row.n_source_cells):
            local.append("row_count")
        if not frame["source_obs_name"].is_unique:
            local.append("source_obs_name_not_unique")
        invalid = ~frame["tcr_assignment_eligible"]
        if frame.loc[invalid, "TCRseq_rebuilt"].eq("yes").any():
            local.append("ineligible_row_has_tcr")
        expected_ab = frame["TRA_cdr3"].ne("") & frame["TRB_cdr3"].ne("")
        expected_gd = frame["TRG_cdr3"].ne("") & frame["TRD_cdr3"].ne("")
        if not frame["has_TRA_TRB_paired_rebuilt"].eq(expected_ab).all():
            local.append("paired_ab_flag")
        if not frame["has_TRG_TRD_paired_rebuilt"].eq(expected_gd).all():
            local.append("paired_gd_flag")
        for chain in CHAINS:
            unavailable = ~frame[f"{chain}_umi_available"]
            if not frame.loc[unavailable, f"{chain}_umis"].isna().all():
                local.append(f"{chain}_unavailable_umi_not_null")
        if local:
            sidecar_errors.append(f"{row.source_gse_id}:{','.join(local)}")

        source_path = (PROJECT_ROOT / row.source_h5ad).resolve()
        stat = source_path.stat()
        if stat.st_size != int(row.source_h5ad_size) or stat.st_mtime_ns != int(row.source_h5ad_mtime_ns):
            current_stat_errors.append(row.source_gse_id)

    check(results, "per_source_sidecar_invariants", not sidecar_errors, "; ".join(sidecar_errors) or "all 14 sidecars pass")
    check(results, "current_h5ad_size_mtime", not current_stat_errors, ", ".join(current_stat_errors) or "all current source files match recorded size/mtime")

    parquet = pq.ParquetFile(REPLACEMENT)
    staged_rows = parquet.metadata.num_rows
    expected_rows = int(replacement_manifest["n_replacement_rows"].sum())
    required_columns = {
        *[f"{chain}_{field}" for chain in CHAINS for field in ("cdr3", "v", "d", "j", "cdr3_nt", "clone_id", "umis", "reads", "umi_available", "read_available")],
        "source_gse_id",
        "source_obs_name",
        "tcr_assignment_eligible",
        "replacement_eligible",
    }
    schema_columns = set(parquet.schema_arrow.names)
    check(results, "staged_row_count", staged_rows == expected_rows == int(run_manifest["staged_replacement"]["n_rows"]), f"rows={staged_rows:,}")
    check(results, "staged_four_chain_schema", required_columns.issubset(schema_columns), f"required={len(required_columns)}; observed={len(schema_columns)}")
    staged_sources = set(
        pd.read_parquet(REPLACEMENT, columns=["source_gse_id"])["source_gse_id"].unique()
    )
    passing_sources = set(summary.loc[summary["rebuild_status"].str.startswith("PASS"), "source_gse_id"])
    check(results, "staged_sources_only_pass", staged_sources == passing_sources, f"sources={len(staged_sources)}")
    digest = sha256_file(REPLACEMENT)
    check(results, "staged_sha256", digest == run_manifest["staged_replacement"]["sha256"], digest)
    check(results, "manifest_records_no_h5ad_write", run_manifest.get("h5ad_write_performed") is False, "h5ad_write_performed=false")

    index_html = REPORT_DIR / "index.html"
    report_pdf = REPORT_DIR / "tcr_join_rebuild_report.pdf"
    html_text = index_html.read_text(encoding="utf-8") if index_html.exists() else ""
    check(results, "report_html", index_html.exists() and "Sample-aware TCR join rebuild" in html_text, str(index_html))
    check(results, "report_documents_no_propagation", "No propagation was performed" in html_text, "explicit mutation-gate statement present")
    check(results, "report_pdf", report_pdf.exists() and report_pdf.stat().st_size > 100_000, f"bytes={report_pdf.stat().st_size if report_pdf.exists() else 0}")

    passed = sum(int(row["passed"]) for row in results)
    payload = {
        "status": "PASS" if passed == len(results) else "FAIL",
        "checks_passed": passed,
        "checks_total": len(results),
        "checks": results,
    }
    OUTPUT_JSON.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    lines = [
        "# TCR join rebuild final validation",
        "",
        f"- Status: **{payload['status']}**",
        f"- Checks passed: **{passed}/{len(results)}**",
        "",
        "| Check | Result | Detail |",
        "|---|---:|---|",
    ]
    for row in results:
        detail = str(row["detail"]).replace("|", "\\|")
        lines.append(f"| {row['check']} | {'PASS' if row['passed'] else 'FAIL'} | {detail} |")
    OUTPUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if payload["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
