#!/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
"""Render the final plus6 HTML/PDF report from existing tables and figures."""

from __future__ import annotations

import html
import os
import subprocess
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
FIGURE_DIR = OUTPUT_ROOT / "figures" / "plus6"
TABLE_DIR = OUTPUT_ROOT / "tables" / "plus6"
PLUS6_REPORT_MD = OUTPUT_ROOT / "plus6_profile_report.md"
PLUS6_REPORT_HTML = OUTPUT_ROOT / "plus6_profile_report.html"
PLUS6_REPORT_PDF = OUTPUT_ROOT / "plus6_profile_report.pdf"


def dataframe_to_markdown_fallback(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if max_rows is not None:
        df = df.head(max_rows)
    columns = [str(column) for column in df.columns]
    rows = [["" if pd.isna(value) else str(value) for value in row] for row in df.itertuples(index=False, name=None)]
    widths = [len(column) for column in columns]
    for row in rows:
        for idx, value in enumerate(row):
            widths[idx] = max(widths[idx], len(value))

    def format_row(row: list[str]) -> str:
        return "| " + " | ".join(value.ljust(widths[idx]) for idx, value in enumerate(row)) + " |"

    header = format_row(columns)
    separator = "| " + " | ".join("-" * widths[idx] for idx in range(len(widths))) + " |"
    body = [format_row(row) for row in rows]
    return "\n".join([header, separator, *body])


def dataframe_to_html_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if max_rows is not None:
        df = df.head(max_rows)
    return df.to_html(index=False, escape=True, border=0, classes="dataframe")


def make_image_card(title: str, image_path: Path) -> str:
    rel = os.path.relpath(image_path, OUTPUT_ROOT)
    return (
        f"<section class='figure-card'><h3>{html.escape(title)}</h3>"
        f"<img src='{html.escape(rel)}' alt='{html.escape(title)}' class='zoomable-figure' data-title='{html.escape(title)}'>"
        f"</section>"
    )


def main() -> None:
    phase3_summary = pd.read_csv(TABLE_DIR / "plus6_phase3_summary.csv")
    prep_summary = pd.read_csv(TABLE_DIR / "plus6_input_compatibility_summary.csv")
    phase4_summary = pd.read_csv(TABLE_DIR / "plus6_phase4_score_summary.csv")
    gse_summary = pd.read_csv(TABLE_DIR / "plus6_phase4_gse_summary.csv")
    tissue_summary = pd.read_csv(TABLE_DIR / "plus6_phase4_tissue_summary.csv")
    annotation_cluster = pd.read_csv(TABLE_DIR / "plus6_simple_annotation_cluster_summary.csv")
    annotation_changes = pd.read_csv(TABLE_DIR / "plus6_simple_annotation_changes.csv")
    annotation_gse = pd.read_csv(TABLE_DIR / "plus6_simple_annotation_by_gse.csv")
    annotation_tissue = pd.read_csv(TABLE_DIR / "plus6_simple_annotation_by_tissue.csv")
    gdt_stats = pd.read_csv(TABLE_DIR / "plus6_gdt_candidate_statistics.csv")
    gdt_overlap = pd.read_csv(TABLE_DIR / "plus6_gdt_candidate_overlap_gt0p4.csv")
    gdt_paired_by_tissue = pd.read_csv(TABLE_DIR / "plus6_gdt_paired_gdtcr_by_tissue.csv")
    gdt_criteria_by_tissue = pd.read_csv(TABLE_DIR / "plus6_gdt_three_criteria_by_tissue.csv")

    md_lines = [
        "# plus6 integrated milestone report",
        "",
        "## Overview",
        "",
        f"- Final plus6 cells: `{int(phase3_summary.loc[0, 'n_cells']):,}`",
        f"- Final plus6 GSEs: `{int(phase3_summary.loc[0, 'n_gses']):,}`",
        f"- Final plus6 Leiden clusters: `{int(phase3_summary.loc[0, 'n_leiden']):,}`",
        "",
        "## Input compatibility",
        "",
        dataframe_to_markdown_fallback(prep_summary),
        "",
        "## Phase 3",
        "",
        dataframe_to_markdown_fallback(phase3_summary),
        "",
        "## Phase 4 score summary",
        "",
        dataframe_to_markdown_fallback(phase4_summary),
        "",
        "## Top GSE composition",
        "",
        dataframe_to_markdown_fallback(gse_summary, max_rows=25),
        "",
        "## Tissue composition",
        "",
        dataframe_to_markdown_fallback(tissue_summary, max_rows=25),
        "",
        "## Annotation cluster summary",
        "",
        dataframe_to_markdown_fallback(annotation_cluster),
        "",
        "## Annotation by GSE",
        "",
        dataframe_to_markdown_fallback(annotation_gse, max_rows=60),
        "",
        "## Annotation by tissue",
        "",
        dataframe_to_markdown_fallback(annotation_tissue, max_rows=60),
        "",
        "## γδ-focused candidate statistics",
        "",
        dataframe_to_markdown_fallback(gdt_stats, max_rows=30),
        "",
        "## γδ-focused overlap breakdown",
        "",
        dataframe_to_markdown_fallback(gdt_overlap, max_rows=20),
        "",
        "## gdT cells with paired gdTCR by tissue",
        "",
        dataframe_to_markdown_fallback(gdt_paired_by_tissue, max_rows=40),
        "",
        "## gdT cells meeting at least one of the three criteria by tissue",
        "",
        dataframe_to_markdown_fallback(gdt_criteria_by_tissue, max_rows=40),
        "",
    ]
    PLUS6_REPORT_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    figure_sections = [
        (
            "Phase 3 Integration",
            "The `plus6` lineage uses a fresh scVI latent space after adding six new datasets. The report keeps only the Phase 3 views that are useful for reviewing study mixing and tissue structure.",
            [
                ("plus6 Phase 3 UMAP by GSE", FIGURE_DIR / "plus6_phase3_umap_by_gse.png"),
                ("plus6 Phase 3 UMAP by tissue", FIGURE_DIR / "plus6_phase3_umap_by_tissue.png"),
            ],
        ),
        (
            "Simple Annotation",
            "Simple labels are kept as a pragmatic interpretation layer on top of the integrated manifold. The detailed Treg/CD8 diagnostic figures are intentionally omitted here.",
            [
                ("plus6 UMAP by corrected simple annotation", FIGURE_DIR / "plus6_umap_by_simple_annotation_corrected.png"),
                ("plus6 TNK marker UMAP panel", FIGURE_DIR / "plus6_tnk_marker_umap_panel.png"),
            ],
        ),
        (
            "Phase 4 Score Space",
            "Phase 4 is reported in raw score space. The report keeps the score UMAPs, the raw TRAB-versus-TRD scatter, and the threshold summaries by tissue and GSE.",
            [
                ("plus6 UMAP by TRD score", FIGURE_DIR / "plus6_phase4_umap_trd_score.png"),
                ("plus6 UMAP by TRAB score", FIGURE_DIR / "plus6_phase4_umap_trab_score.png"),
                ("plus6 UMAP by TRD minus TRAB", FIGURE_DIR / "plus6_phase4_umap_trd_minus_trab.png"),
                ("plus6 Raw TRAB-versus-TRD score space", FIGURE_DIR / "plus6_phase4_trab_vs_trd_raw.png"),
                ("plus6 TRAB versus TRD colored by paired TCR", FIGURE_DIR / "plus6_phase4_trab_vs_trd_paired_tcr.png"),
                ("plus6 Phase 4 score distributions", FIGURE_DIR / "plus6_phase4_score_distributions.png"),
                ("High-Confidence TRD-over-TRAB Candidates by Tissue", FIGURE_DIR / "plus6_phase4_trd_over_trab_by_tissue_barplot.png"),
                ("High-Confidence TRD-over-TRAB Candidates by GSE", FIGURE_DIR / "plus6_phase4_trd_over_trab_by_gse_barplot.png"),
                ("Broad TRD-Enriched Candidates by Tissue", FIGURE_DIR / "plus6_phase4_trd_gt_0p1_by_tissue_barplot.png"),
                ("Broad TRD-Enriched Candidates by GSE", FIGURE_DIR / "plus6_phase4_trd_gt_0p1_by_gse_barplot.png"),
            ],
        ),
        (
            "γδ-Focused Views",
            "This section keeps the targeted TCR and sorted-γδ overlays and the tissue-level γδ statistics used in the current review.",
            [
                ("plus6 paired TRA/TRB, paired TRG/TRD, and Sorted_gdT highlight UMAPs", FIGURE_DIR / "plus6_umap_paired_tcr_sorted_gdt.png"),
            ],
        ),
    ]

    html_parts = [
        "<!DOCTYPE html>",
        "<html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'><title>plus6 integrated milestone report</title>",
        "<style>"
        ":root{--bg:#f4efe6;--paper:#fffdf8;--ink:#1f2430;--muted:#596273;--accent:#8c2f39;--line:#d8cfc2;--shadow:0 18px 40px rgba(31,36,48,0.08);}*{box-sizing:border-box}"
        "body{margin:0;font-family:Georgia,'Times New Roman',serif;color:var(--ink);background:radial-gradient(circle at top left,rgba(140,47,57,0.10),transparent 28%),radial-gradient(circle at top right,rgba(31,92,91,0.12),transparent 26%),var(--bg);line-height:1.6}"
        ".page{width:min(1220px,calc(100vw - 32px));margin:20px auto 40px}"
        ".hero{background:linear-gradient(135deg,rgba(140,47,57,0.94),rgba(31,92,91,0.96));color:white;border-radius:24px;padding:36px 40px;box-shadow:var(--shadow)}"
        ".hero h1{margin:0 0 12px;font-size:clamp(34px,5vw,54px);line-height:1.02;letter-spacing:-0.03em}"
        ".hero p{margin:0 0 10px;max-width:920px;font-size:18px}"
        ".metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:14px;margin:18px 0 0}"
        ".metric{background:rgba(255,255,255,0.12);border:1px solid rgba(255,255,255,0.18);border-radius:18px;padding:14px 16px}"
        ".metric .value{display:block;font-size:26px;font-weight:700;letter-spacing:-0.03em}.metric .label{display:block;font-size:13px;text-transform:uppercase;letter-spacing:0.08em;opacity:0.84}"
        ".section{margin-top:28px;background:var(--paper);border:1px solid var(--line);border-radius:24px;padding:28px;box-shadow:var(--shadow)}"
        ".section h2{margin:0 0 10px;font-size:28px;letter-spacing:-0.02em}.section p{margin:0 0 12px;color:var(--muted);font-size:17px}"
        ".section code{background:rgba(31,36,48,0.06);padding:0 5px;border-radius:4px;font-size:0.95em}"
        ".figure-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:24px}"
        ".figure-card{margin-top:14px;border:1px solid var(--line);border-radius:18px;padding:16px;background:white}"
        ".figure-card img{width:100%;height:auto;border-radius:12px;border:1px solid var(--line);background:#faf7f2;cursor:zoom-in;transition:transform 0.18s ease,box-shadow 0.18s ease}"
        ".figure-card img:hover{transform:translateY(-1px);box-shadow:0 12px 24px rgba(31,36,48,0.12)}"
        "table{border-collapse:collapse;font-size:12px;width:100%} th,td{border:1px solid #ddd;padding:4px 6px;text-align:left;vertical-align:top}"
        ".lightbox{position:fixed;inset:0;display:none;align-items:center;justify-content:center;padding:24px;background:rgba(15,18,24,0.88);z-index:9999}"
        ".lightbox.open{display:flex}.lightbox-inner{width:min(96vw,1500px);max-height:96vh;display:flex;flex-direction:column;gap:10px;align-items:center}"
        ".lightbox img{max-width:100%;max-height:calc(96vh - 72px);width:auto;height:auto;border-radius:14px;background:white;box-shadow:0 22px 50px rgba(0,0,0,0.35)}"
        ".lightbox-caption{color:white;font-size:15px;text-align:center}.lightbox-close{position:absolute;top:18px;right:20px;border:0;border-radius:999px;width:40px;height:40px;font-size:28px;line-height:1;cursor:pointer;color:white;background:rgba(255,255,255,0.14)}"
        "</style></head><body><main class='page'>",
        "<section class='hero'>",
        "<h1>plus6 Integrated Milestone Report</h1>",
        "<p>This report summarizes the parallel plus6 lineage after adding six new datasets to the validated integrated T/NK base. It focuses on integration structure, simple annotation, raw Phase 4 score space, and γδ-focused review statistics.</p>",
        "<div class='metrics'>",
        f"<div class='metric'><span class='value'>{int(phase3_summary.loc[0, 'n_cells']):,}</span><span class='label'>Integrated cells</span></div>",
        f"<div class='metric'><span class='value'>{int(phase3_summary.loc[0, 'n_gses']):,}</span><span class='label'>Represented GSEs</span></div>",
        f"<div class='metric'><span class='value'>{int(phase3_summary.loc[0, 'n_leiden']):,}</span><span class='label'>Leiden clusters</span></div>",
        f"<div class='metric'><span class='value'>{int(gdt_stats.loc[gdt_stats['metric']=='sorted_gdT_true','value'].iloc[0]):,}</span><span class='label'>Sorted gdT cells</span></div>",
        f"<div class='metric'><span class='value'>{int(gdt_stats.loc[gdt_stats['metric']=='paired_TRG_TRD_not_doublet','value'].iloc[0]):,}</span><span class='label'>Paired gdTCR, not doublet</span></div>",
        f"<div class='metric'><span class='value'>{int(gdt_stats.loc[gdt_stats['metric']=='union_any_of_three_criteria','value'].iloc[0]):,}</span><span class='label'>Union of 3 gdT criteria</span></div>",
        "</div></section>",
        "<section class='section'><h2>Compatibility</h2><p>The six new inputs were standardized into the plus6 lineage before a fresh scVI integration. The table below records the compatibility status used for the merge.</p>",
        dataframe_to_html_table(prep_summary, max_rows=20),
        "</section>",
        "<section class='section'><h2>Quantitative Summary</h2><p>Phase 3, Phase 4, and simple-annotation summaries are retained as compact tables. The detailed Treg/CD8 diagnostic panels are intentionally omitted from this final report.</p>",
        "<h3>Phase 3</h3>",
        dataframe_to_html_table(phase3_summary),
        "<h3>Phase 4 Score Summary</h3>",
        dataframe_to_html_table(phase4_summary),
        "<h3>Annotation Cluster Summary</h3>",
        dataframe_to_html_table(annotation_cluster, max_rows=40),
        "<h3>Top GSE Composition</h3>",
        dataframe_to_html_table(gse_summary, max_rows=30),
        "<h3>Tissue Composition</h3>",
        dataframe_to_html_table(tissue_summary, max_rows=30),
        "</section>",
        "<section class='section'><h2>γδ Statistics</h2><p>The γδ-focused review layer includes direct counts, overlap structure at <code>TRD - TRAB &gt; 0.4</code>, and tissue-specific summaries for annotated gdT cells.</p>",
        "<h3>γδ-Focused Candidate Statistics</h3>",
        dataframe_to_html_table(gdt_stats, max_rows=30),
        "<h3>Overlap Breakdown</h3>",
        dataframe_to_html_table(gdt_overlap, max_rows=20),
        "<h3>gdT Cells With Paired gdTCR by Tissue</h3>",
        dataframe_to_html_table(gdt_paired_by_tissue, max_rows=40),
        "<h3>gdT Cells Meeting At Least One of the Three Criteria by Tissue</h3>",
        dataframe_to_html_table(gdt_criteria_by_tissue, max_rows=40),
        "</section>",
    ]
    for section_title, section_text, figures in figure_sections:
        html_parts.append(f"<section class='section'><h2>{html.escape(section_title)}</h2><p>{html.escape(section_text)}</p><div class='figure-grid'>")
        for title, path in figures:
            if path.exists():
                html_parts.append(make_image_card(title, path))
        html_parts.append("</div></section>")
    html_parts.extend(
        [
            "<div class='lightbox' id='lightbox'><button class='lightbox-close' id='lightbox-close' aria-label='Close'>&times;</button><div class='lightbox-inner'><img id='lightbox-image' alt='Expanded figure'><div class='lightbox-caption' id='lightbox-caption'></div></div></div>",
            "<script>const lightbox=document.getElementById('lightbox');const lightboxImage=document.getElementById('lightbox-image');const lightboxCaption=document.getElementById('lightbox-caption');const closeButton=document.getElementById('lightbox-close');document.querySelectorAll('.zoomable-figure').forEach(img=>{img.addEventListener('click',()=>{lightboxImage.src=img.src;lightboxCaption.textContent=img.dataset.title||img.alt||'';lightbox.classList.add('open');document.body.style.overflow='hidden';});});function closeLightbox(){lightbox.classList.remove('open');lightboxImage.removeAttribute('src');document.body.style.overflow='';}closeButton.addEventListener('click',closeLightbox);lightbox.addEventListener('click',event=>{if(event.target===lightbox){closeLightbox();}});document.addEventListener('keydown',event=>{if(event.key==='Escape'&&lightbox.classList.contains('open')){closeLightbox();}});</script>",
            "</main></body></html>",
        ]
    )
    PLUS6_REPORT_HTML.write_text("\n".join(html_parts), encoding="utf-8")

    chrome_candidates = ["/usr/bin/google-chrome", "/usr/bin/google-chrome-stable"]
    chrome = next((Path(path) for path in chrome_candidates if Path(path).exists()), None)
    if chrome is None:
        raise FileNotFoundError("google-chrome not found for PDF export.")
    subprocess.run(
        [
            str(chrome),
            "--headless",
            "--disable-gpu",
            "--print-to-pdf-no-header",
            f"--print-to-pdf={PLUS6_REPORT_PDF}",
            str(PLUS6_REPORT_HTML),
        ],
        check=True,
    )

    print(f"saved_md\t{PLUS6_REPORT_MD}")
    print(f"saved_html\t{PLUS6_REPORT_HTML}")
    print(f"saved_pdf\t{PLUS6_REPORT_PDF}")


if __name__ == "__main__":
    main()
