#!/usr/bin/env python3
"""Build the Phase 3/4 training-dataset HTML report from editable markdown.

The markdown source is narrative-first and can contain inline author comments in
double quotation marks. Those quoted comment lines are ignored during HTML
rendering so the markdown can stay review-friendly.
"""

from __future__ import annotations

import csv
import html
import re
from pathlib import Path

import h5py


PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
FIGURE_DIR = OUTPUT_ROOT / "figures"
TABLE_DIR = OUTPUT_ROOT / "tables"
REPORT_PATH = OUTPUT_ROOT / "phase3_phase4_training_dataset_report.html"
REPORT_TEXT_MD = OUTPUT_ROOT / "phase3_phase4_training_dataset_report_text.md"
INTEGRATED_H5AD = Path("/ssd/tnk_phase3/Integrated_dataset/integrated.h5ad")

FIGURE_MAP = {
    "Phase 3 UMAP by GSE": "phase3_umap_by_gse.png",
    "Phase 3 UMAP by Tissue": "phase3_umap_by_tissue.png",
    "Phase 3 UMAP by Leiden Cluster": "phase3_umap_by_leiden.png",
    "Phase 3 T/NK Marker Panel": "phase3_tnk_marker_umap_panel.png",
    "Phase 4 Raw Score UMAP Overlays": "phase4_umap_raw_score_overlays.png",
    "Raw TRAB-versus-TRD Score Space": "phase4_trab_vs_trd_raw_only.png",
    "Raw TRAB-versus-TRD Space with Paired TRA/TRB Context": "phase4_trab_vs_trd_paired_tratrb_vs_no_tcr_raw_only.png",
    "Selected Marker Genes on Raw TRAB-versus-TRD Space": "phase4_trab_vs_trd_selected_marker_panel_raw_only.png",
    "High-Confidence TRD-over-TRAB Candidates by Tissue": "phase4_trab_minus_trd_lt_neg0p6_by_tissue_barplot.png",
    "High-Confidence TRD-over-TRAB Candidates by GSE": "phase4_trab_minus_trd_lt_neg0p6_by_gse_barplot.png",
    "Broad TRD-Enriched Candidates by Tissue": "phase4_trd_gt_0p1_by_tissue_barplot.png",
    "Broad TRD-Enriched Candidates by GSE": "phase4_trd_gt_0p1_by_gse_barplot.png",
}

DEFAULT_FIGURE_CAPTIONS = {
    "Phase 3 UMAP by Tissue": (
        "Cells are colored by the harmonized tissue field `tissue_corrected`, "
        "showing how blood and tissue compartments distribute across the integrated manifold."
    ),
}


def read_total_cells_and_tissues() -> tuple[int, int]:
    """Read total cell count and tissue category count from the tissue summary."""
    total_cells = 0
    tissue_count = 0
    path = TABLE_DIR / "tissue_correction" / "tissue_corrected_value_counts.csv"
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            tissue_count += 1
            total_cells += int(row["cell_n"])
    return total_cells, tissue_count


def read_integrated_scope() -> tuple[int, int]:
    """Read the number of GSEs and Leiden clusters from the integrated object."""
    with h5py.File(INTEGRATED_H5AD, "r") as handle:
        obs = handle["obs"]
        gse_values = obs["source_gse_id"].asstr()[:]
        leiden_values = obs["leiden"].asstr()[:]
    return len(set(gse_values)), len(set(leiden_values))


def normalize_heading(title: str) -> str:
    """Normalize a markdown heading before figure lookup."""
    normalized = title.replace("`", "").strip()
    normalized = re.sub(r"\s*\d[\d,]*\s*$", "", normalized)
    normalized = re.sub(r"\s+", " ", normalized)
    return normalized.strip()


def strip_comment_lines(lines: list[str]) -> list[str]:
    """Drop inline review comments wrapped in double quotes."""
    kept: list[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('"') and stripped.endswith('"') and len(stripped) >= 2:
            continue
        kept.append(line.rstrip("\n"))
    return kept


def flush_paragraph(buffer: list[str], target: list[str]) -> None:
    """Flush a buffered paragraph into the target list."""
    if not buffer:
        return
    target.append(" ".join(part.strip() for part in buffer if part.strip()))
    buffer.clear()


def parse_report_markdown(path: Path) -> dict:
    """Parse the limited report-markdown structure into a nested dictionary."""
    lines = strip_comment_lines(path.read_text(encoding="utf-8").splitlines())
    report = {"title": "", "intro": [], "metrics": [], "sections": []}
    current_section: dict | None = None
    current_subsection: dict | None = None
    paragraph_buffer: list[str] = []

    def flush_current_paragraph() -> None:
        nonlocal paragraph_buffer
        if current_subsection is not None:
            flush_paragraph(paragraph_buffer, current_subsection["paragraphs"])
        elif current_section is not None:
            flush_paragraph(paragraph_buffer, current_section["paragraphs"])
        else:
            flush_paragraph(paragraph_buffer, report["intro"])

    for raw_line in lines + [""]:
        line = raw_line.strip()
        if not line:
            flush_current_paragraph()
            continue
        if line.startswith("# "):
            flush_current_paragraph()
            report["title"] = line[2:].strip()
            continue
        if line.startswith("## "):
            flush_current_paragraph()
            current_section = {"title": line[3:].strip(), "paragraphs": [], "subsections": []}
            report["sections"].append(current_section)
            current_subsection = None
            continue
        if line.startswith("### "):
            flush_current_paragraph()
            if current_section is None:
                continue
            current_subsection = {"title": line[4:].strip(), "paragraphs": []}
            current_section["subsections"].append(current_subsection)
            continue
        if current_section is not None and current_section["title"] == "Report Metrics" and line.startswith("- "):
            metric_line = line[2:].strip()
            if ":" in metric_line:
                label, value = metric_line.split(":", 1)
                report["metrics"].append({"label": label.strip(), "value": value.strip()})
            continue
        paragraph_buffer.append(line)

    return report


def paragraph_html(text: str) -> str:
    """Render one markdown-lite paragraph into HTML."""
    escaped = html.escape(text)
    escaped = escaped.replace("`", "<code>", 1).replace("`", "</code>", 1) if escaped.count("`") >= 2 else escaped
    while "`" in escaped:
        escaped = escaped.replace("`", "<code>", 1).replace("`", "</code>", 1)
    return f"<p>{escaped}</p>"


def figure_card(title: str, caption: str) -> str:
    """Render one figure card block for a mapped subsection title."""
    filename = FIGURE_MAP[normalize_heading(title)]
    safe_title = html.escape(normalize_heading(title))
    safe_caption = html.escape(caption)
    return f"""
    <section class="figure-card">
      <h3>{safe_title}</h3>
      <img src="figures/{filename}" alt="{safe_title}" class="zoomable-figure" data-title="{safe_title}">
      <p>{safe_caption}</p>
    </section>
    """


def render_metrics(metrics: list[dict]) -> str:
    """Render the report metrics as hero cards."""
    cards = []
    for metric in metrics:
        cards.append(
            f'<div class="metric"><span class="value">{html.escape(metric["value"])}</span>'
            f'<span class="label">{html.escape(metric["label"])}</span></div>'
        )
    return "".join(cards)


def render_section(section: dict) -> str:
    """Render one report section from parsed markdown."""
    blocks: list[str] = [f"<section class=\"section\"><h2>{html.escape(section['title'])}</h2>"]
    for paragraph in section["paragraphs"]:
        blocks.append(paragraph_html(paragraph))
    for subsection in section["subsections"]:
        blocks.append(f"<div class=\"subsection\"><h3>{html.escape(normalize_heading(subsection['title']))}</h3>")
        for paragraph in subsection["paragraphs"]:
            blocks.append(paragraph_html(paragraph))
        figure_key = normalize_heading(subsection["title"])
        if figure_key in FIGURE_MAP:
            caption = subsection["paragraphs"][-1] if subsection["paragraphs"] else DEFAULT_FIGURE_CAPTIONS.get(figure_key, figure_key)
            blocks.append(figure_card(subsection["title"], caption))
        blocks.append("</div>")
    blocks.append("</section>")
    return "\n".join(blocks)


def build_html(report: dict) -> str:
    """Build the full HTML report body from parsed markdown."""
    intro_html = "\n".join(paragraph_html(paragraph) for paragraph in report["intro"])
    section_html = "\n".join(
        render_section(section)
        for section in report["sections"]
        if section["title"] != "Report Metrics"
    )
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Phase 3-4 Training Dataset Report</title>
  <style>
    :root {{
      --bg: #f4efe6;
      --paper: #fffdf8;
      --ink: #1f2430;
      --muted: #596273;
      --accent: #8c2f39;
      --line: #d8cfc2;
      --shadow: 0 18px 40px rgba(31, 36, 48, 0.08);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: Georgia, "Times New Roman", serif;
      color: var(--ink);
      background:
        radial-gradient(circle at top left, rgba(140, 47, 57, 0.10), transparent 28%),
        radial-gradient(circle at top right, rgba(31, 92, 91, 0.12), transparent 26%),
        var(--bg);
      line-height: 1.6;
    }}
    .page {{
      width: min(1220px, calc(100vw - 32px));
      margin: 20px auto 40px;
    }}
    .hero {{
      background: linear-gradient(135deg, rgba(140, 47, 57, 0.94), rgba(31, 92, 91, 0.96));
      color: white;
      border-radius: 24px;
      padding: 36px 40px;
      box-shadow: var(--shadow);
    }}
    .hero h1 {{
      margin: 0 0 12px;
      font-size: clamp(34px, 5vw, 54px);
      line-height: 1.02;
      letter-spacing: -0.03em;
    }}
    .hero p {{
      margin: 0 0 10px;
      max-width: 920px;
      font-size: 18px;
    }}
    .metrics {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 14px;
      margin: 18px 0 0;
    }}
    .metric {{
      background: rgba(255, 255, 255, 0.12);
      border: 1px solid rgba(255, 255, 255, 0.18);
      border-radius: 18px;
      padding: 14px 16px;
    }}
    .metric .value {{
      display: block;
      font-size: 26px;
      font-weight: 700;
      letter-spacing: -0.03em;
    }}
    .metric .label {{
      display: block;
      font-size: 13px;
      text-transform: uppercase;
      letter-spacing: 0.08em;
      opacity: 0.84;
    }}
    .section {{
      margin-top: 28px;
      background: var(--paper);
      border: 1px solid var(--line);
      border-radius: 24px;
      padding: 28px;
      box-shadow: var(--shadow);
    }}
    .section h2 {{
      margin: 0 0 10px;
      font-size: 28px;
      letter-spacing: -0.02em;
    }}
    .section p {{
      margin: 0 0 12px;
      color: var(--muted);
      font-size: 17px;
    }}
    .section code {{
      background: rgba(31, 36, 48, 0.06);
      padding: 0 5px;
      border-radius: 4px;
      font-size: 0.95em;
    }}
    .subsection {{
      margin-top: 18px;
      padding-top: 12px;
      border-top: 1px solid rgba(216, 207, 194, 0.65);
    }}
    .subsection h3 {{
      margin: 0 0 10px;
      font-size: 21px;
    }}
    .figure-card {{
      margin-top: 14px;
      border: 1px solid var(--line);
      border-radius: 18px;
      padding: 16px;
      background: white;
    }}
    .figure-card img {{
      width: 100%;
      height: auto;
      border-radius: 12px;
      border: 1px solid var(--line);
      background: #faf7f2;
      cursor: zoom-in;
      transition: transform 0.18s ease, box-shadow 0.18s ease;
    }}
    .figure-card img:hover {{
      transform: translateY(-1px);
      box-shadow: 0 12px 24px rgba(31, 36, 48, 0.12);
    }}
    .figure-card p {{
      margin: 12px 0 0;
      font-size: 15px;
    }}
    .lightbox {{
      position: fixed;
      inset: 0;
      display: none;
      align-items: center;
      justify-content: center;
      padding: 24px;
      background: rgba(15, 18, 24, 0.88);
      z-index: 9999;
    }}
    .lightbox.open {{
      display: flex;
    }}
    .lightbox-inner {{
      width: min(96vw, 1500px);
      max-height: 96vh;
      display: flex;
      flex-direction: column;
      gap: 10px;
      align-items: center;
    }}
    .lightbox img {{
      max-width: 100%;
      max-height: calc(96vh - 72px);
      width: auto;
      height: auto;
      border-radius: 14px;
      background: white;
      box-shadow: 0 22px 50px rgba(0, 0, 0, 0.35);
    }}
    .lightbox-caption {{
      color: white;
      font-size: 15px;
      text-align: center;
    }}
    .lightbox-close {{
      position: absolute;
      top: 18px;
      right: 20px;
      border: 0;
      border-radius: 999px;
      width: 40px;
      height: 40px;
      font-size: 28px;
      line-height: 1;
      cursor: pointer;
      color: white;
      background: rgba(255, 255, 255, 0.14);
    }}
    .footer {{
      margin-top: 20px;
      color: var(--muted);
      font-size: 14px;
      text-align: center;
    }}
  </style>
</head>
<body>
  <main class="page">
    <section class="hero">
      <h1>{html.escape(report["title"])}</h1>
      {intro_html}
      <div class="metrics">
        {render_metrics(report["metrics"])}
      </div>
    </section>

    {section_html}

    <div class="footer">
      Report generated from the revised markdown narrative and the validated Phase 3/4 PNG outputs in <code>Integrated_dataset/figures/</code>.
    </div>

    <div class="lightbox" id="lightbox" aria-hidden="true">
      <button class="lightbox-close" id="lightbox-close" aria-label="Close figure">&times;</button>
      <div class="lightbox-inner">
        <img id="lightbox-image" alt="">
        <div class="lightbox-caption" id="lightbox-caption"></div>
      </div>
    </div>
  </main>
  <script>
    const lightbox = document.getElementById("lightbox");
    const lightboxImage = document.getElementById("lightbox-image");
    const lightboxCaption = document.getElementById("lightbox-caption");
    const lightboxClose = document.getElementById("lightbox-close");

    function closeLightbox() {{
      lightbox.classList.remove("open");
      lightbox.setAttribute("aria-hidden", "true");
      lightboxImage.removeAttribute("src");
      lightboxImage.removeAttribute("alt");
      lightboxCaption.textContent = "";
    }}

    document.querySelectorAll(".zoomable-figure").forEach((img) => {{
      img.addEventListener("click", () => {{
        lightboxImage.src = img.src;
        lightboxImage.alt = img.alt;
        lightboxCaption.textContent = img.dataset.title || img.alt;
        lightbox.classList.add("open");
        lightbox.setAttribute("aria-hidden", "false");
      }});
    }});

    lightbox.addEventListener("click", (event) => {{
      if (event.target === lightbox) {{
        closeLightbox();
      }}
    }});

    lightboxClose.addEventListener("click", closeLightbox);

    document.addEventListener("keydown", (event) => {{
      if (event.key === "Escape" && lightbox.classList.contains("open")) {{
        closeLightbox();
      }}
    }});
  </script>
</body>
</html>
"""


def main() -> None:
    """Build the HTML report from the editable markdown text file."""
    report = parse_report_markdown(REPORT_TEXT_MD)
    total_cells, tissue_count = read_total_cells_and_tissues()
    gse_count, leiden_count = read_integrated_scope()
    metrics = {metric["label"]: metric["value"] for metric in report["metrics"]}
    metrics.setdefault("Integrated cells", f"`{total_cells:,}`")
    metrics.setdefault("Represented GSEs", f"`{gse_count}`")
    metrics.setdefault("Leiden clusters", f"`{leiden_count}`")
    metrics.setdefault("Tissue groups", f"`{tissue_count}`")
    report["metrics"] = [{"label": label, "value": value} for label, value in metrics.items()]
    REPORT_PATH.write_text(build_html(report), encoding="utf-8")
    print(REPORT_PATH)


if __name__ == "__main__":
    main()
