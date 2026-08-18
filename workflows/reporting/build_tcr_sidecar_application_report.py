#!/usr/bin/env python3
"""Build the static HTML/PDF QC report for the atlas TCR sidecar transaction."""

from __future__ import annotations

import html
import json
import shutil
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = ROOT / "Integrated_dataset/tables/tcr_sidecar_application"
LOG_DIR = ROOT / "Integrated_dataset/logs/tcr_sidecar_application"
FIGURE_DIR = ROOT / "Integrated_dataset/figures/tcr_sidecar_application"
REPORT_DIR = ROOT / "gdT_prediction/gdtai_v4_2_tcr_sidecar_application"
ASSET_DIR = REPORT_DIR / "assets"
SUMMARY = TABLE_DIR / "tcr_replacement_by_source.csv"
MANIFEST = LOG_DIR / "tcr_h5ad_candidate_manifest.json"


def fmt_int(value: object) -> str:
    return f"{int(value):,}"


def save(fig: plt.Figure, name: str) -> None:
    path = FIGURE_DIR / name
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    shutil.copy2(path, ASSET_DIR / name)


def make_figures(frame: pd.DataFrame) -> None:
    work = frame.sort_values("n_rebuilt_any_ab").reset_index(drop=True)
    y = np.arange(len(work))
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 6.6), sharey=True)
    for ax, old, new, title in (
        (axes[0], "n_old_any_ab", "n_rebuilt_any_ab", "Any productive TRA or TRB"),
        (axes[1], "n_old_paired_ab", "n_rebuilt_paired_ab", "Paired productive TRA/TRB"),
    ):
        ax.barh(y + 0.18, work[old], height=0.34, color="#9CA3AF", label="Legacy atlas")
        ax.barh(y - 0.18, work[new], height=0.34, color="#087E8B", label="Validated rebuild")
        ax.set_title(title, fontsize=12)
        ax.set_xlabel("Cells")
        ax.grid(axis="x", color="#E5E7EB", linewidth=0.8)
        ax.set_axisbelow(True)
        ax.spines[["top", "right", "left"]].set_visible(False)
    axes[0].set_yticks(y, work["source_gse_id"])
    axes[1].legend(frameon=False, loc="lower right")
    fig.suptitle("Canonical alpha-beta TCR calls before and after sidecar application", fontsize=15)
    fig.tight_layout()
    save(fig, "old_vs_rebuilt_ab_calls.png")

    change = frame.copy()
    change["net_any_ab"] = change["n_rebuilt_any_ab"] - change["n_old_any_ab"]
    change = change.sort_values("net_any_ab").reset_index(drop=True)
    colors = np.where(change["net_any_ab"] >= 0, "#087E8B", "#D1495B")
    fig, ax = plt.subplots(figsize=(9.5, 6.0))
    ax.barh(change["source_gse_id"], change["net_any_ab"], color=colors, height=0.62)
    ax.axvline(0, color="#111827", linewidth=1)
    ax.set_xlabel("Net change in cells with productive TRA or TRB")
    ax.set_title("Sample-aware rebuilding restores missing calls and removes stale calls")
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    save(fig, "net_productive_ab_change.png")

    prevalence = frame.copy()
    prevalence["any_ab_pct"] = 100 * prevalence["n_rebuilt_any_ab"] / prevalence["n_atlas_rows"]
    prevalence["paired_ab_pct"] = 100 * prevalence["n_rebuilt_paired_ab"] / prevalence["n_atlas_rows"]
    prevalence = prevalence.sort_values("any_ab_pct").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(9.5, 6.0))
    ax.barh(prevalence["source_gse_id"], prevalence["any_ab_pct"], color="#4C78A8", height=0.62, label="Any TRA/TRB")
    ax.barh(prevalence["source_gse_id"], prevalence["paired_ab_pct"], color="#F4A261", height=0.34, label="Paired TRA/TRB")
    ax.set_xlabel("Affected atlas cells with rebuilt evidence (%)")
    ax.set_title("Rebuilt alpha-beta TCR coverage varies by dataset")
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="lower right")
    ax.spines[["top", "right", "left"]].set_visible(False)
    save(fig, "rebuilt_ab_prevalence_by_source.png")


def source_table(frame: pd.DataFrame) -> str:
    rows = []
    for row in frame.sort_values("source_gse_id").itertuples(index=False):
        rows.append(
            "<tr>"
            f"<td>{html.escape(row.source_gse_id)}</td>"
            f"<td>{fmt_int(row.n_atlas_rows)}</td>"
            f"<td>{fmt_int(row.n_old_any_ab)}</td>"
            f"<td>{fmt_int(row.n_rebuilt_any_ab)}</td>"
            f"<td>{fmt_int(row.n_old_paired_ab)}</td>"
            f"<td>{fmt_int(row.n_rebuilt_paired_ab)}</td>"
            f"<td>{fmt_int(row.n_old_only_any_ab)}</td>"
            f"<td>{fmt_int(row.n_rebuilt_only_any_ab)}</td>"
            f"<td>{fmt_int(row.n_assignment_ineligible)}</td>"
            "</tr>"
        )
    return (
        "<table><thead><tr><th>Source</th><th>Atlas cells</th><th>Old any AB</th>"
        "<th>Rebuilt any AB</th><th>Old paired</th><th>Rebuilt paired</th>"
        "<th>Old-only</th><th>Rebuilt-only</th><th>Ineligible</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table>"
    )


def build_html(frame: pd.DataFrame, manifest: dict[str, object]) -> Path:
    checks = manifest["checks"]
    check_rows = "".join(
        f"<tr><td>{html.escape(name.replace('_', ' '))}</td><td><span class='pass'>PASS</span></td></tr>"
        for name, passed in checks.items()
        if passed
    )
    source_sum = frame.sum(numeric_only=True)
    content = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Full-atlas TCR metadata correction</title>
<style>
@page {{ size: A4 landscape; margin: 10mm; }}
* {{ box-sizing: border-box; }} body {{ margin:0; color:#17202A; background:#F4F6F7; font:14px/1.45 Arial,sans-serif; }}
main {{ max-width:1180px; margin:auto; background:white; }} section {{ padding:28px 38px; border-bottom:1px solid #D5DBDB; page-break-inside:avoid; }}
.hero {{ background:#123B45; color:white; padding:42px 44px; }} h1 {{ margin:0 0 8px; font-size:30px; }} h2 {{ font-size:21px; color:#123B45; margin:0 0 14px; }} h3 {{ font-size:16px; margin:16px 0 8px; }}
.subtitle {{ font-size:17px; max-width:900px; }} .grid {{ display:grid; grid-template-columns:repeat(4,1fr); gap:12px; margin-top:22px; }}
.metric {{ border-left:4px solid #18A999; background:#F8FAFA; color:#17202A; padding:12px; }} .metric b {{ display:block; font-size:22px; }}
.note {{ background:#EAF5F3; border-left:4px solid #18A999; padding:12px 14px; margin:12px 0; }}
.warning {{ background:#FFF4E5; border-left:4px solid #F4A261; padding:12px 14px; margin:12px 0; }}
figure {{ margin:12px 0 0; }} img {{ width:100%; max-height:570px; object-fit:contain; }} figcaption {{ color:#566573; font-size:12px; }}
table {{ border-collapse:collapse; width:100%; font-size:10px; table-layout:auto; }} th {{ background:#123B45; color:white; }} th,td {{ padding:6px 7px; border:1px solid #D5DBDB; text-align:right; white-space:nowrap; }} th:first-child,td:first-child {{ text-align:left; }}
.pass {{ color:#087E8B; font-weight:bold; }} code {{ font-size:11px; overflow-wrap:anywhere; }} .page {{ page-break-before:always; }}
@media(max-width:800px) {{ .grid {{ grid-template-columns:1fr 1fr; }} section {{ padding:22px 18px; }} table {{ display:block; overflow-x:auto; }} }}
</style></head><body><main>
<header class="hero"><h1>Full-atlas TCR metadata correction</h1>
<div class="subtitle">Sample-aware validated receptor calls were applied to a separate 5.93-million-cell atlas candidate while preserving expression, integration, and the canonical rollback object.</div>
<div class="grid">
<div class="metric"><b>{fmt_int(manifest['n_cells'])}</b>atlas cells</div>
<div class="metric"><b>{fmt_int(manifest['n_affected_atlas_rows'])}</b>corrected rows</div>
<div class="metric"><b>{fmt_int(manifest['n_whole_atlas_any_ab_tcr'])}</b>productive AB</div>
<div class="metric"><b>{fmt_int(manifest['n_whole_atlas_paired_tra_trb'])}</b>paired TRA/TRB</div>
</div></header>
<section><h2>Decision summary</h2>
<div class="note"><b>{manifest['status']}</b>: all {len(checks)} transaction checks passed. The candidate was renamed from <code>.partial</code> only after exact sidecar, unaffected-value, matrix, embedding, and backed-read validation.</div>
<p><b>Candidate:</b> <code>{html.escape(manifest['candidate'])}</code><br><b>SHA-256:</b> <code>{manifest['candidate_sha256']}</code></p>
<p>The canonical symlink was not changed. The metadata-corrected input and original canonical atlas remain rollback points.</p></section>
<section><h2>What changed</h2>
<p>For 14 repaired sources, atlas <code>source_gse_id + original_cell_id</code> was joined one-to-one to sidecar <code>source_gse_id + source_obs_name</code>. Forty-six canonical chain and logical fields were replaced, including empty values that remove stale calls. Twenty-seven provenance fields record coverage, eligibility, join status, UMI/read availability, productive-contig multiplicity, source file, and selected contig.</p>
<div class="warning">Unavailable UMI/read support remains null, not zero. The 109 atlas-present ambiguous GSE235863 rows are blank and TCR-assignment-ineligible. The 886,462 sidecar rows outside the frozen T/NK atlas were not rescued.</div></section>
<section class="page"><h2>Before and after</h2><figure><img src="assets/old_vs_rebuilt_ab_calls.png"><figcaption>Legacy versus rebuilt productive and paired alpha-beta calls among affected atlas cells.</figcaption></figure></section>
<section class="page"><h2>Where calls changed</h2><figure><img src="assets/net_productive_ab_change.png"><figcaption>Positive values are recovered validated calls; GSE311112 also removes 681 stale productive-alpha-beta assignments.</figcaption></figure></section>
<section class="page"><h2>Dataset coverage</h2><figure><img src="assets/rebuilt_ab_prevalence_by_source.png"><figcaption>These percentages describe only cells retained in the frozen atlas and should not be interpreted as tissue-level biological abundance.</figcaption></figure></section>
<section class="page"><h2>Source-level audit</h2>{source_table(frame)}
<p><b>Affected-source totals:</b> {fmt_int(source_sum['n_rebuilt_any_ab'])} rebuilt productive-alpha-beta cells and {fmt_int(source_sum['n_rebuilt_paired_ab'])} rebuilt paired TRA/TRB cells.</p></section>
<section class="page"><h2>Validation matrix</h2><table><thead><tr><th>Check</th><th>Status</th></tr></thead><tbody>{check_rows}</tbody></table></section>
<section><h2>Interpretation and next boundary</h2>
<p>This transaction repairs receptor metadata; it does not retrain gdTAI, change expression scores, alter cell membership, or rerun scVI, UMAP, or Leiden. Previous NK-boundary and TCR-ground-truth summaries are now stale because they consumed legacy receptor fields.</p>
<p><b>Next:</b> regenerate TCR-derived truth/control labels and boundary/NK audits from this checksum-bound candidate. Canonical publication remains a separate supervised decision after those reports are reviewed.</p></section>
</main></body></html>"""
    path = REPORT_DIR / "index.html"
    path.write_text(content, encoding="utf-8")
    return path


def main() -> None:
    for path in (FIGURE_DIR, REPORT_DIR, ASSET_DIR):
        path.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(SUMMARY)
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if manifest["status"] != "PASS_TCR_H5AD_CANDIDATE" or not all(manifest["checks"].values()):
        raise RuntimeError("TCR candidate manifest is not fully passing")
    make_figures(frame)
    html_path = build_html(frame, manifest)
    pdf_path = REPORT_DIR / "tcr_sidecar_application_report.pdf"
    subprocess.run(
        [
            "google-chrome",
            "--headless",
            "--no-sandbox",
            "--disable-gpu",
            "--disable-crash-reporter",
            f"--print-to-pdf={pdf_path}",
            f"file://{html_path}",
        ],
        check=True,
    )
    print(json.dumps({"html": str(html_path), "pdf": str(pdf_path), "figures": 3}, indent=2))


if __name__ == "__main__":
    main()
