#!/usr/bin/env python3
"""Build the static QC report for the sample-aware TCR join rebuild."""

from __future__ import annotations

# TNK_WORKFLOW_BOOTSTRAP
import sys as _tnk_sys
from pathlib import Path as _TnkPath

_TNK_PROJECT_ROOT = _TnkPath(__file__).resolve().parents[2]
if str(_TNK_PROJECT_ROOT) not in _tnk_sys.path:
    _tnk_sys.path.insert(0, str(_TNK_PROJECT_ROOT))

import html
import json
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = PROJECT_ROOT / "Integrated_dataset/tables/tcr_join_rebuild"
FIGURE_DIR = PROJECT_ROOT / "Integrated_dataset/figures/tcr_join_rebuild"
REPORT_DIR = PROJECT_ROOT / "gdT_prediction/gdtai_v4_2_tcr_join_rebuild"
ASSET_DIR = REPORT_DIR / "assets"
SUMMARY_CSV = TABLE_DIR / "source_rebuild_summary.csv"
REUSE_CSV = TABLE_DIR / "receptor_reuse_audit.csv"
REPLACEMENT_MANIFEST = TABLE_DIR / "validated_tcr_replacement_manifest.csv"
RUN_MANIFEST = PROJECT_ROOT / "Integrated_dataset/logs/tcr_join_rebuild/run_manifest.json"

STATUS_LABEL = {
    "PASS_REBUILT": "Pass",
    "PASS_PARTIAL_WITH_QUARANTINE": "Pass, partial quarantine",
    "QUARANTINED_NO_SAFE_RAW_JOIN": "Quarantined",
    "QUARANTINED_INSUFFICIENT_RAW_TO_RNA_MATCH": "Quarantined",
}


def ensure_dirs() -> None:
    for path in (FIGURE_DIR, REPORT_DIR, ASSET_DIR):
        path.mkdir(parents=True, exist_ok=True)


def fmt_int(value: object) -> str:
    if pd.isna(value):
        return "NA"
    return f"{int(value):,}"


def fmt_pct(value: object, digits: int = 1) -> str:
    if pd.isna(value):
        return "NA"
    return f"{100 * float(value):.{digits}f}%"


def save_figure(fig: plt.Figure, name: str) -> Path:
    path = FIGURE_DIR / name
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    shutil.copy2(path, ASSET_DIR / name)
    return path


def plot_old_vs_rebuilt(summary: pd.DataFrame) -> Path:
    work = summary.sort_values("n_old_any_tcr", ascending=True).reset_index(drop=True)
    y = np.arange(len(work))
    fig, ax = plt.subplots(figsize=(9.4, 6.6))
    ax.barh(y + 0.18, work["n_old_any_tcr"], height=0.34, color="#9CA3AF", label="Legacy call")
    colors = np.where(work["rebuild_status"].str.startswith("PASS"), "#087E8B", "#D1495B")
    ax.barh(y - 0.18, work["n_rebuilt_any_tcr"], height=0.34, color=colors, label="Rebuilt raw-VDJ call")
    ax.set_yticks(y, work["source_gse_id"])
    ax.set_xlabel("Cells with any productive TCR")
    ax.set_title("Legacy assignments versus sample-aware raw-VDJ rebuild")
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="lower right")
    ax.spines[["top", "right", "left"]].set_visible(False)
    return save_figure(fig, "old_vs_rebuilt_tcr_calls.png")


def calculate_umi_categories(summary: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for source in summary["source_gse_id"]:
        frame = pd.read_parquet(
            TABLE_DIR / "per_source" / f"{source}.parquet",
            columns=[*[f"{chain}_cdr3" for chain in ("TRA", "TRB", "TRG", "TRD")], *[f"{chain}_umis" for chain in ("TRA", "TRB", "TRG", "TRD")]],
        )
        one = 0
        ge_two = 0
        unavailable = 0
        for chain in ("TRA", "TRB", "TRG", "TRD"):
            called = frame[f"{chain}_cdr3"].ne("")
            umi = pd.to_numeric(frame[f"{chain}_umis"], errors="coerce")
            one += int((called & umi.eq(1)).sum())
            ge_two += int((called & umi.ge(2)).sum())
            unavailable += int((called & umi.isna()).sum())
        rows.append({"source_gse_id": source, "UMI = 1": one, "UMI >= 2": ge_two, "UMI unavailable": unavailable})
    return pd.DataFrame(rows)


def plot_umi_support(umi: pd.DataFrame) -> Path:
    work = umi.copy()
    work["total"] = work[["UMI = 1", "UMI >= 2", "UMI unavailable"]].sum(axis=1)
    work = work.loc[work["total"].gt(0)].sort_values("total").reset_index(drop=True)
    y = np.arange(len(work))
    fig, ax = plt.subplots(figsize=(9.4, 5.6))
    left = np.zeros(len(work))
    for column, color in (("UMI >= 2", "#087E8B"), ("UMI = 1", "#F4A261"), ("UMI unavailable", "#9CA3AF")):
        ax.barh(y, work[column], left=left, color=color, height=0.62, label=column)
        left += work[column].to_numpy()
    ax.set_yticks(y, work["source_gse_id"])
    ax.set_xlabel("Selected productive chain calls")
    ax.set_title("Per-chain UMI evidence is retained, not imputed")
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, ncol=3, loc="lower right")
    ax.spines[["top", "right", "left"]].set_visible(False)
    return save_figure(fig, "umi_support_by_source.png")


def plot_raw_key_recovery(summary: pd.DataFrame) -> Path:
    work = summary.loc[summary["n_eligible_raw_cell_keys"].gt(0)].copy()
    work = work.sort_values("raw_to_rna_match_fraction").reset_index(drop=True)
    colors = np.where(work["rebuild_status"].str.startswith("PASS"), "#087E8B", "#D1495B")
    fig, ax = plt.subplots(figsize=(9.4, 5.4))
    ax.barh(work["source_gse_id"], 100 * work["raw_to_rna_match_fraction"], color=colors, height=0.62)
    ax.axvline(50, color="#111827", linestyle="--", linewidth=1.0, label="Source quarantine floor")
    ax.set_xlim(0, 101)
    ax.set_xlabel("Eligible raw VDJ cell keys recovered in RNA metadata (%)")
    ax.set_title("Raw-to-RNA key recovery")
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="lower right")
    ax.spines[["top", "right", "left"]].set_visible(False)
    return save_figure(fig, "raw_key_recovery.png")


def plot_donor_reuse(reuse: pd.DataFrame) -> Path:
    work = reuse.loc[
        reuse["n_paired_ab"].gt(0) & reuse["n_cells_in_pairs_across_donors"].notna()
    ].copy()
    work["fraction"] = work["n_cells_in_pairs_across_donors"] / work["n_paired_ab"]
    work = work.sort_values("fraction").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(9.4, 5.4))
    ax.barh(work["source_gse_id"], 100 * work["fraction"], color="#4C78A8", height=0.62)
    ax.set_xlabel("Paired-AB cells in an exact pair seen in >1 donor (%)")
    ax.set_title("Residual exact-pair sharing after sample-aware joining")
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    return save_figure(fig, "cross_donor_receptor_reuse.png")


def table_html(headers: list[str], rows: list[list[str]], classes: str = "") -> str:
    head = "".join(f"<th>{html.escape(header)}</th>" for header in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>" for row in rows
    )
    return f'<table class="{classes}"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>'


def build_html(summary: pd.DataFrame, reuse: pd.DataFrame, replacement: pd.DataFrame, manifest: dict[str, object]) -> None:
    passing = summary["rebuild_status"].str.startswith("PASS")
    quarantined_sources = summary.loc[~passing, "source_gse_id"].tolist()
    unchanged_sources = int(summary["source_h5ad_unchanged"].sum())
    staged = manifest["staged_replacement"]
    total_chain_calls = int(summary.loc[passing, "n_rebuilt_chain_calls_with_umi"].sum())
    total_one = int(summary.loc[passing, "n_rebuilt_chain_calls_umi_eq_1"].sum())
    source_rows: list[list[str]] = []
    for row in summary.itertuples(index=False):
        status_class = "pass" if str(row.rebuild_status).startswith("PASS") else "hold"
        source_rows.append(
            [
                html.escape(row.source_gse_id),
                f'<span class="status {status_class}">{html.escape(STATUS_LABEL.get(row.rebuild_status, row.rebuild_status))}</span>',
                fmt_int(row.n_old_any_tcr),
                fmt_int(row.n_rebuilt_any_tcr),
                fmt_int(row.n_rebuilt_paired_ab),
                fmt_pct(row.raw_to_rna_match_fraction),
                fmt_pct(row.author_call_raw_concordance),
                fmt_int(row.n_rebuilt_chain_calls_with_umi),
            ]
        )

    conflict_rows: list[list[str]] = []
    merged_reuse = summary.merge(reuse, on="source_gse_id", how="left")
    for row in merged_reuse.loc[
        merged_reuse["n_paired_ab"].gt(0) | merged_reuse["n_author_nk_cells"].gt(0)
    ].itertuples(index=False):
        cross_fraction = (
            row.n_cells_in_pairs_across_donors / row.n_paired_ab if row.n_paired_ab else np.nan
        )
        nk_fraction = row.n_rebuilt_ab_in_author_nk / row.n_author_nk_cells if row.n_author_nk_cells else np.nan
        conflict_rows.append(
            [
                html.escape(row.source_gse_id),
                fmt_int(row.n_paired_ab),
                fmt_int(row.n_cells_in_pairs_across_donors),
                fmt_pct(cross_fraction, 2),
                fmt_int(row.n_author_nk_cells),
                fmt_int(row.n_rebuilt_ab_in_author_nk),
                fmt_pct(nk_fraction, 2),
            ]
        )

    replacement_rows = [
        [
            html.escape(row.source_gse_id),
            fmt_int(row.n_replacement_rows),
            fmt_int(row.n_tcr_assignment_eligible_rows),
            fmt_int(row.n_rebuilt_any_tcr),
        ]
        for row in replacement.itertuples(index=False)
    ]
    repaired_rules = {
        "GSE125527": "Published old-to-new participant map + tissue + barcode core",
        "GSE228597": "Pooled-library suffix map + barcode core",
        "GSE287541": "Round + participant visit + barcode core",
    }
    repair_rows: list[list[str]] = []
    repaired = summary.loc[summary["source_gse_id"].isin(repaired_rules)]
    for row in repaired.itertuples(index=False):
        repair_rows.append(
            [
                html.escape(row.source_gse_id),
                html.escape(repaired_rules[row.source_gse_id]),
                fmt_int(row.n_raw_source_files),
                fmt_int(row.n_matched_raw_cell_keys),
                fmt_pct(row.raw_to_rna_match_fraction),
                fmt_pct(row.author_call_raw_concordance),
                f"{float(row.author_call_enrichment_over_sample_rotation):,.1f}x"
                if pd.notna(row.author_call_enrichment_over_sample_rotation)
                else "NA",
                fmt_int(row.n_rebuilt_chain_calls_with_umi),
            ]
        )

    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>V4.2 TCR Join Rebuild QC</title>
<style>
@page {{ size: A4 landscape; margin: 11mm; }}
* {{ box-sizing: border-box; }}
body {{ margin: 0; color: #111827; font: 14px/1.45 Arial, sans-serif; background: #f5f7f8; }}
main {{ max-width: 1180px; margin: 0 auto; background: white; }}
.page {{ min-height: 185mm; padding: 28px 34px; page-break-after: always; }}
.page:last-child {{ page-break-after: auto; }}
.kicker {{ color: #087E8B; font-weight: 700; text-transform: uppercase; font-size: 12px; }}
h1 {{ margin: 7px 0 12px; font-size: 30px; letter-spacing: 0; }}
h2 {{ margin: 0 0 14px; font-size: 22px; letter-spacing: 0; }}
h3 {{ margin: 0 0 8px; font-size: 16px; letter-spacing: 0; }}
p {{ max-width: 92ch; margin: 7px 0 13px; }}
.lede {{ font-size: 17px; color: #374151; }}
.metrics {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 10px; margin: 22px 0; }}
.metric {{ border-top: 4px solid #087E8B; padding: 12px 10px; background: #f3f7f7; }}
.metric b {{ display: block; font-size: 24px; }}
.metric span {{ color: #4B5563; font-size: 12px; }}
.notice {{ border-left: 4px solid #D1495B; background: #fff6f5; padding: 10px 13px; margin: 14px 0; }}
.good {{ border-left-color: #087E8B; background: #f1f8f8; }}
.grid2 {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; align-items: start; }}
.figure {{ margin: 8px 0 0; }}
.figure img {{ width: 100%; display: block; }}
.full-figure img {{ height: 150mm; object-fit: contain; }}
.caption {{ color: #4B5563; font-size: 11px; margin-top: 5px; }}
table {{ width: 100%; border-collapse: collapse; table-layout: fixed; font-size: 10px; }}
th {{ background: #E8F0F1; text-align: left; padding: 6px 5px; overflow-wrap: anywhere; }}
td {{ padding: 5px; border-bottom: 1px solid #E5E7EB; vertical-align: top; overflow-wrap: anywhere; }}
.status {{ display: inline-block; padding: 2px 5px; font-size: 9px; font-weight: 700; }}
.status.pass {{ color: #075E66; background: #DDF1F2; }}
.status.hold {{ color: #9B2C2C; background: #FDE5E1; }}
.steps {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 10px; margin: 18px 0; }}
.step {{ border: 1px solid #D1D5DB; padding: 11px; }}
.step b {{ color: #087E8B; display: block; margin-bottom: 5px; }}
code {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 11px; }}
.small {{ font-size: 11px; color: #4B5563; }}
@media print {{ body {{ background: white; }} main {{ max-width: none; }} .page {{ padding: 0; }} }}
</style>
</head>
<body><main>
<section class="page">
<div class="kicker">gdTAI V4.2 data-integrity gate</div>
<h1>Sample-aware TCR join rebuild</h1>
<p class="lede">Productive raw VDJ contigs were rebuilt against RNA metadata using only <code>sample_id + barcode_core</code>. Per-chain UMI and read support are preserved for TRA, TRB, TRG, and TRD; missing support remains null.</p>
<div class="metrics">
  <div class="metric"><b>{int(passing.sum())}/{len(summary)}</b><span>sources pass full or partial rebuild</span></div>
  <div class="metric"><b>{fmt_int(staged['n_rows'])}</b><span>rows staged for metadata replacement</span></div>
  <div class="metric"><b>{fmt_int(summary.loc[passing, 'n_rebuilt_any_tcr'].sum())}</b><span>cells with validated productive TCR</span></div>
  <div class="metric"><b>{fmt_int(total_chain_calls)}</b><span>productive chain calls with observed UMI</span></div>
</div>
<div class="notice good"><b>QC gate:</b> the rebuild produced a checksum-bound sidecar and {unchanged_sources}/{len(summary)} source H5AD size/mtime pairs remained unchanged.</div>
<div class="notice {'good' if not quarantined_sources else ''}"><b>{'No source remains quarantined.' if not quarantined_sources else 'Sources still quarantined: ' + ', '.join(quarantined_sources) + '.'}</b> GSE125527 uses the published old-to-new patient map plus tissue; GSE228597 uses the pooled-library suffix map; GSE287541 uses Cell Ranger VDJ contigs reconstructed from all public TCR SRA runs with UMI/read support.</div>
</section>

<section class="page">
<h2>Previously quarantined sources</h2>
<p>All three repairs use an independently reconstructable sample key before barcode matching. Exact author/raw CDR3 agreement is reported as supporting evidence, not used to copy legacy calls into the rebuilt sidecar.</p>
{table_html(['Source','RNA-to-VDJ sample rule','Raw files','Matched cells','Raw-key recovery','Exact CDR3','Fold over wrong-sample rotation','Chain calls with UMI'], repair_rows)}
<div class="notice good"><b>GSE287541:</b> all 46 public TCR runs were reconstructed with Cell Ranger VDJ. The two N131V01 libraries remain distinct through their round identifiers.</div>
<div class="notice"><b>GSE125527 quantitative support:</b> the public productive receptor tables do not contain UMI/read columns. The official raw archive contains 71 TCR runs totaling 150.6 GiB compressed; support is therefore retained as unavailable rather than imputed.</div>
</section>

<section class="page">
<h2>Legacy assignments versus raw-VDJ rebuild</h2>
<div class="figure full-figure"><img src="assets/old_vs_rebuilt_tcr_calls.png" alt="Legacy versus rebuilt productive TCR counts"></div>
<div class="caption">Legacy and rebuilt counts refer to the same 14 source objects; differences reflect sample-aware raw-VDJ reconstruction, not classifier performance.</div>
</section>

<section class="page">
<h2>Rebuild method</h2>
<div class="steps">
  <div class="step"><b>1. Parse raw VDJ</b>Require cell-associated, high-confidence, productive contigs with nonempty CDR3.</div>
  <div class="step"><b>2. Select within chain</b>Highest UMI, then reads, full-length status, then stable contig id.</div>
  <div class="step"><b>3. Join safely</b>Canonical sample id plus nucleotide barcode core. Barcode-only fallback is forbidden.</div>
  <div class="step"><b>4. Fail closed</b>Duplicate RNA keys, collapsed raw libraries, or poor source recovery cannot provide TCR truth.</div>
</div>
<div class="grid2">
  <div><div class="figure"><img src="assets/raw_key_recovery.png" alt="Raw key recovery"></div><div class="caption">The primary source floor is 50% raw-key recovery. RNA-filtered multi-assay objects may pass below that floor only with at least 100 testable and 100 confirmed author calls, at least 50% exact CDR3 agreement, and at least 20-fold exact-match enrichment over a deterministic wrong-sample rotation.</div></div>
  <div><div class="figure"><img src="assets/umi_support_by_source.png" alt="UMI support"></div><div class="caption">GSE125527 and GSE235863 author receptor tables do not expose UMI/read columns, so their support remains null rather than zero. GSE125527 raw TCR recovery would require 71 SRA runs (150.6 GiB compressed) and is not needed to establish the deterministic patient+tissue join. Across UMI-auditable calls, {fmt_int(total_one)} selected chains have UMI=1.</div></div>
</div>
<p class="small">The four-chain sidecar contains CDR3 amino acid/nucleotide sequence, V/D/J calls, clone id, UMI, reads, support-availability flags, selected contig id, productive-contig multiplicity, and raw source-file provenance.</p>
</section>

<section class="page">
<h2>Source-level outcome</h2>
{table_html(['Source','Outcome','Legacy any TCR','Rebuilt any TCR','Rebuilt paired AB','Raw-key recovery','Exact CDR3 confirmation','Chain calls with UMI'], source_rows)}
<div class="notice good"><b>GSE311112 correction:</b> explicit support for relapse and BTK-clone sample names recovers 35,266 cells, compared with 9,663 under the incomplete historical parser.</div>
<p class="small">GSE235863 has 110 duplicated RNA <code>sample_id + barcode_core</code> rows. Their rebuilt receptor fields are blank and they are excluded from receptor truth, but retained in the replacement table so stale assignments can be cleared.</p>
</section>

<section class="page">
<h2>Residual biological and annotation checks</h2>
<div class="grid2">
  <div><div class="figure"><img src="assets/cross_donor_receptor_reuse.png" alt="Cross donor exact receptor reuse"></div><div class="caption">This sharing is already present in raw sample-resolved VDJ. It may include public or semi-invariant receptors; it is reported rather than used as automatic evidence of a join failure. Donor-level values are not reported for GSE178882 paired libraries or GSE243572 pooled libraries because donor identity is unresolved.</div></div>
  <div>
    <h3>Exact pair and author-label audit</h3>
    {table_html(['Source','Paired AB','Cross-donor cells','Cross-donor %','Author NK','NK + rebuilt AB','NK + AB %'], conflict_rows)}
  </div>
</div>
<div class="notice"><b>Interpretation:</b> raw productive TRA/TRB in an author-NK label is a conflict to review for NKT-like cells, doublets, or annotation error. It is not erased by transcriptomic NK labeling and is not automatically admitted as model truth.</div>
</section>

<section class="page">
<h2>Staged metadata replacement</h2>
<p>The staged sidecar contains all rows from the {int(passing.sum())} passing sources. Receptor-truth eligibility is separate from metadata-replacement eligibility, allowing unsafe legacy calls to be cleared without converting ambiguous rows into negative training truth.</p>
{table_html(['Source','Replacement rows','TCR-truth eligible','Validated productive TCR'], replacement_rows)}
<div class="notice good"><b>Artifact:</b> <code>{html.escape(str(staged['path']))}</code><br><b>SHA-256:</b> <code>{html.escape(str(staged['sha256']))}</code></div>
<div class="notice"><b>No propagation was performed.</b> Writing these fields into source or milestone H5AD metadata remains a high-risk mutation gate. The replacement should be applied atomically with backups only after explicit approval, followed by harmonized-metadata and milestone overlap checks.</div>
<p class="small">This report evaluates data integrity of TCR assignment. It does not evaluate or promote a gdTAI classifier, threshold, or atlas inference.</p>
</section>
</main></body></html>"""
    (REPORT_DIR / "index.html").write_text(document, encoding="utf-8")


def main() -> None:
    ensure_dirs()
    summary = pd.read_csv(SUMMARY_CSV)
    reuse = pd.read_csv(REUSE_CSV)
    unresolved = {
        "GSE178882": "paired libraries; donor grouping is not retained in RNA metadata",
        "GSE243572": "pooled libraries; donor identity is not cell-resolved",
    }
    if "donor_key_status" not in reuse.columns:
        reuse.insert(1, "donor_key_status", "resolved_or_sample_derived")
    for source, reason in unresolved.items():
        mask = reuse["source_gse_id"].eq(source)
        reuse.loc[mask, "donor_key_status"] = reason
        reuse.loc[
            mask,
            ["n_pairs_across_donors", "n_cells_in_pairs_across_donors", "max_donors_per_pair"],
        ] = np.nan
    reuse.to_csv(REUSE_CSV, index=False)
    replacement = pd.read_csv(REPLACEMENT_MANIFEST)
    manifest = json.loads(RUN_MANIFEST.read_text(encoding="utf-8"))
    plot_old_vs_rebuilt(summary)
    plot_umi_support(calculate_umi_categories(summary))
    plot_raw_key_recovery(summary)
    plot_donor_reuse(reuse)
    build_html(summary, reuse, replacement, manifest)
    print(REPORT_DIR / "index.html")


if __name__ == "__main__":
    main()
