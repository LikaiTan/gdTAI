#!/usr/bin/env python3
"""
GSE26 V4 pipeline.

Stages currently implemented:
- Discovery and metadata profiling per GSE using matrix/ + suppl/.
- Candidate expression-file inventory with filtered-first classification.
- Metadata-derived donor/tissue summary (when metadata files exist).
- Deterministic load-plan selection (filtered-first + sample-level de-dup).

This script is intended to be iteratively revised (single script policy).
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
import csv
import gzip
import json
import re
import subprocess
import shutil
import stat
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd
import h5py

ROOT = Path(__file__).resolve().parents[2]
DOWNLOADS = ROOT / "downloads"
CFG = ROOT / "configs" / "datasets" / "legacy_gse26.json"
OUT = ROOT / "analysis_26GSE_V4"

# Exclude files that are usually already-integrated/aggregate objects.
AGGREGATE_HINTS = (
    "combined",
    "integrated",
    "seuratobj",
    "seurat",
    "aggr",
    "aggregation",
)

# Keep strict expression file types for discovery.
EXPR_SUFFIXES = (
    ".h5",
    ".h5ad",
    ".mtx",
    ".mtx.gz",
)

META_SUFFIXES = (
    ".csv",
    ".csv.gz",
    ".tsv",
    ".tsv.gz",
    ".txt",
    ".txt.gz",
)

DONOR_KEYS = (
    "donor",
    "patient",
    "subject",
    "individual",
    "case",
)

TISSUE_KEYS = (
    "tissue",
    "site",
    "organ",
    "source",
    "sampletype",
    "sample_type",
    "compartment",
    "location",
)


def load_targets() -> List[str]:
    cfg = json.loads(CFG.read_text())
    return list(cfg["targets"]["gse_list"])


@dataclass
class MetadataSummary:
    donor_col: Optional[str]
    donor_count: Optional[int]
    tissue_col: Optional[str]
    tissue_values: List[str]
    meta_file: Optional[str]



def iter_files(path: Path) -> Iterable[Path]:
    if not path.exists():
        return []
    return (p for p in path.rglob("*") if p.is_file())



def is_expression_file(name: str) -> bool:
    lower = name.lower()
    return lower.endswith(EXPR_SUFFIXES)



def is_aggregate_like(path: Path) -> bool:
    n = path.name.lower()
    return any(h in n for h in AGGREGATE_HINTS)



def classify_matrix_kind(path: Path) -> str:
    n = path.name.lower()
    if "filtered" in n:
        return "filtered"
    if "raw" in n:
        return "raw"
    if path.suffix == ".h5ad":
        return "h5ad_unknown"
    return "unknown"


def modality_label(path: Path) -> str:
    n = path.name.lower()
    if "_rna_" in n or "rna_" in n or "_gex_" in n or "gex_" in n:
        return "rna"
    if "_csp_" in n or "_cite_" in n or "cite_" in n:
        return "protein"
    return "unknown"


def include_by_gse_rule(gse: str, path: Path) -> bool:
    n = path.name.lower()
    # User rule: GSE161918 should be RNA only.
    if gse == "GSE161918" and "_csp_" in n:
        return False
    return True


def source_priority(path_str: str) -> int:
    """Lower number means higher preference."""
    p = path_str.lower()
    if "/extracted_final_" in p:
        return 0
    if "/extracted_correct_" in p:
        return 1
    if "/extracted_filtered_" in p:
        return 2
    if "/extracted_v3/" in p:
        return 3
    if "/suppl/" in p:
        return 4
    return 5


def sample_key_from_path(path: str) -> str:
    """Normalize file name to a canonical sample key for deduplication."""
    pp = Path(path)
    n = pp.name.lower()
    if n in ("matrix.mtx.gz", "matrix.mtx", "matrix.mtx"):
        parent = pp.parent.name.lower()
        if parent in ("filtered_feature_bc_matrix", "raw_feature_bc_matrix", "feature_bc_matrix"):
            return pp.parents[1].name.lower()
        return parent
    for suffix in (
        ".filtered_feature_bc_matrix.h5ad",
        ".raw_feature_bc_matrix.h5ad",
        "_filtered_feature_bc_matrix.h5ad",
        "_raw_feature_bc_matrix.h5ad",
        ".filtered_feature_bc_matrix.h5",
        ".raw_feature_bc_matrix.h5",
        "_filtered_feature_bc_matrix.h5",
        "_raw_feature_bc_matrix.h5",
        ".matrix.mtx.gz",
        ".mtx.gz",
        ".h5ad",
        ".h5",
        ".mtx",
    ):
        if n.endswith(suffix):
            n = n[: -len(suffix)]
            break
    return n



def find_metadata_files(gse_dir: Path) -> List[Path]:
    suppl = gse_dir / "suppl"
    out: List[Path] = []
    if not suppl.exists():
        return out
    for p in suppl.iterdir():
        if not p.is_file():
            continue
        n = p.name.lower()
        if ("meta" in n or "metadata" in n) and n.endswith(META_SUFFIXES):
            out.append(p)
    return sorted(out)



def safe_read_table(path: Path) -> Optional[pd.DataFrame]:
    n = path.name.lower()
    sep = "\t" if n.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")) else ","
    try:
        return pd.read_csv(path, sep=sep, compression="gzip" if n.endswith(".gz") else None, low_memory=False)
    except Exception:
        return None



def choose_column(cols: Sequence[str], keys: Sequence[str]) -> Optional[str]:
    ranked: List[Tuple[int, str]] = []
    for c in cols:
        cl = c.lower().replace(" ", "").replace("-", "_")
        score = 0
        for k in keys:
            if k in cl:
                score += 1
        if score > 0:
            ranked.append((score, c))
    if not ranked:
        return None
    ranked.sort(key=lambda x: (-x[0], x[1]))
    return ranked[0][1]



def summarize_metadata(meta_files: List[Path]) -> MetadataSummary:
    for mf in meta_files:
        df = safe_read_table(mf)
        if df is None or df.empty:
            continue

        donor_col = choose_column(df.columns.tolist(), DONOR_KEYS)
        tissue_col = choose_column(df.columns.tolist(), TISSUE_KEYS)

        donor_count: Optional[int] = None
        if donor_col is not None:
            donor_count = int(df[donor_col].dropna().astype(str).nunique())

        tissue_values: List[str] = []
        if tissue_col is not None:
            vals = df[tissue_col].dropna().astype(str).str.strip()
            tissue_values = sorted(v for v in vals.unique() if v)
            tissue_values = tissue_values[:30]

        # Accept first metadata file that gives any useful signal.
        if donor_col is not None or tissue_col is not None:
            return MetadataSummary(
                donor_col=donor_col,
                donor_count=donor_count,
                tissue_col=tissue_col,
                tissue_values=tissue_values,
                meta_file=str(mf.relative_to(ROOT)),
            )

    return MetadataSummary(None, None, None, [], None)



def discover_gse(gse: str) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    gse_dir = DOWNLOADS / gse
    matrix_dir = gse_dir / "matrix"
    suppl_dir = gse_dir / "suppl"

    summary: Dict[str, object] = {
        "GSE": gse,
        "exists": gse_dir.exists(),
        "matrix_dir": matrix_dir.exists(),
        "suppl_dir": suppl_dir.exists(),
        "expected_series_matrix_files": 0,
        "candidate_expr_files": 0,
        "filtered_candidates": 0,
        "raw_candidates": 0,
        "unknown_candidates": 0,
        "aggregate_like_excluded": 0,
        "metadata_file": "",
        "donor_column": "",
        "donor_count": "",
        "tissue_column": "",
        "tissue_values": "",
    }

    details: List[Dict[str, object]] = []

    if matrix_dir.exists():
        summary["expected_series_matrix_files"] = len(list(matrix_dir.glob("*series_matrix*.txt.gz")))

    meta_summary = summarize_metadata(find_metadata_files(gse_dir))
    if meta_summary.meta_file:
        summary["metadata_file"] = meta_summary.meta_file
    if meta_summary.donor_col:
        summary["donor_column"] = meta_summary.donor_col
        summary["donor_count"] = meta_summary.donor_count if meta_summary.donor_count is not None else ""
    if meta_summary.tissue_col:
        summary["tissue_column"] = meta_summary.tissue_col
    if meta_summary.tissue_values:
        summary["tissue_values"] = "|".join(meta_summary.tissue_values)

    if not gse_dir.exists():
        return summary, details

    candidate_roots = [suppl_dir]
    for ex in ("extracted_v3", f"extracted_final_{gse}", f"extracted_correct_{gse}", f"extracted_filtered_{gse}"):
        p = suppl_dir / ex
        if p.exists():
            candidate_roots.append(p)

    seen: set[str] = set()

    for root in candidate_roots:
        for p in iter_files(root):
            rel = str(p.relative_to(ROOT))
            if rel in seen:
                continue
            seen.add(rel)

            if not is_expression_file(p.name):
                continue

            if not include_by_gse_rule(gse, p):
                details.append(
                    {
                        "GSE": gse,
                        "path": rel,
                        "kind": "excluded_rule",
                        "matrix_kind": "excluded",
                        "modality": modality_label(p),
                    }
                )
                continue

            if is_aggregate_like(p):
                summary["aggregate_like_excluded"] += 1
                details.append(
                    {
                        "GSE": gse,
                        "path": rel,
                        "kind": "excluded_aggregate",
                        "matrix_kind": "excluded",
                        "modality": modality_label(p),
                    }
                )
                continue

            mkind = classify_matrix_kind(p)
            summary["candidate_expr_files"] += 1
            if mkind == "filtered":
                summary["filtered_candidates"] += 1
            elif mkind == "raw":
                summary["raw_candidates"] += 1
            else:
                summary["unknown_candidates"] += 1

            details.append(
                {
                    "GSE": gse,
                    "path": rel,
                    "kind": "candidate",
                    "matrix_kind": mkind,
                    "modality": modality_label(p),
                }
            )

    return summary, details



def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    cols = list(rows[0].keys())
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


def read_mtx_shape(path: Path) -> Tuple[Optional[int], Optional[int]]:
    opener = gzip.open if path.name.endswith(".gz") else open
    try:
        with opener(path, "rt") as f:
            for line in f:
                if line.startswith("%"):
                    continue
                parts = line.strip().split()
                if len(parts) >= 3:
                    n_rows = int(parts[0])
                    n_cols = int(parts[1])
                    # MTX is usually genes x cells.
                    return n_cols, n_rows
                break
    except Exception:
        return None, None
    return None, None


def read_h5_shape(path: Path) -> Tuple[Optional[int], Optional[int]]:
    # 10x H5 commonly stores matrix/shape = [n_genes, n_cells]
    try:
        with h5py.File(path, "r") as h5:
            if "matrix" in h5 and "shape" in h5["matrix"]:
                shape = h5["matrix"]["shape"][()]
                if len(shape) == 2:
                    return int(shape[1]), int(shape[0])
    except Exception:
        return None, None
    return None, None


def build_plan_stats() -> None:
    plan_path = OUT / "manifests" / "gse_selected_load_plan.csv"
    if not plan_path.exists():
        raise FileNotFoundError(f"Missing load plan: {plan_path}")
    plan = pd.read_csv(plan_path)
    rows: List[Dict[str, object]] = []
    for _, r in plan.iterrows():
        p = ROOT / str(r["path"])
        n_obs = None
        n_vars = None
        if p.suffix.lower() == ".h5" or p.name.lower().endswith(".h5"):
            n_obs, n_vars = read_h5_shape(p)
        elif p.name.lower().endswith((".mtx", ".mtx.gz")):
            n_obs, n_vars = read_mtx_shape(p)
        rows.append(
            {
                "GSE": r["GSE"],
                "path": r["path"],
                "matrix_kind": r["matrix_kind"],
                "modality": r["modality"],
                "n_obs_header": n_obs if n_obs is not None else "",
                "n_vars_header": n_vars if n_vars is not None else "",
            }
        )

    write_csv(OUT / "manifests" / "gse_selected_plan_stats.csv", rows)

    ps = pd.DataFrame(rows)
    ps["n_obs_header"] = pd.to_numeric(ps["n_obs_header"], errors="coerce")
    per = (
        ps.groupby("GSE", as_index=False)
        .agg(
            selected_files=("path", "count"),
            header_cells_sum=("n_obs_header", "sum"),
            header_genes_median=("n_vars_header", "median"),
            raw_files=("matrix_kind", lambda x: int((x == "raw").sum())),
            filtered_files=("matrix_kind", lambda x: int((x == "filtered").sum())),
        )
    )
    per["likely_inflated_raw"] = (per["raw_files"] > 0) & (per["header_cells_sum"] > 1_000_000)
    per.to_csv(OUT / "manifests" / "gse_selected_plan_summary.csv", index=False)
    print("Plan stats complete")
    print(per.sort_values("header_cells_sum", ascending=False).head(10).to_string(index=False))


def scaffold_projects(targets: List[str]) -> None:
    plan_path = OUT / "manifests" / "gse_selected_load_plan.csv"
    sum_path = OUT / "manifests" / "gse_discovery_summary.csv"
    if not plan_path.exists() or not sum_path.exists():
        raise FileNotFoundError("Run --stage discovery first to generate manifests.")

    plan = pd.read_csv(plan_path)
    summary = pd.read_csv(sum_path).set_index("GSE", drop=False)
    projects_root = OUT / "scanpy_projects"
    projects_root.mkdir(parents=True, exist_ok=True)

    # Reuse existing Scanpy helper scripts as project assets.
    template_scripts = ROOT / "scanpy_tcr_project" / "scripts"
    template_readme = ROOT / "scanpy_tcr_project" / "README.md"

    created = 0
    for gse in targets:
        grp = plan[plan["GSE"] == gse].copy()
        gse_dir = projects_root / gse
        (gse_dir / "config").mkdir(parents=True, exist_ok=True)
        (gse_dir / "scripts").mkdir(parents=True, exist_ok=True)
        (gse_dir / "logs").mkdir(parents=True, exist_ok=True)
        (gse_dir / "outputs").mkdir(parents=True, exist_ok=True)
        (gse_dir / "manifests").mkdir(parents=True, exist_ok=True)

        # Keep selected file list per project.
        grp_sorted = grp.sort_values(["source_priority", "sample_key", "path"]) if not grp.empty else grp
        grp_sorted.to_csv(gse_dir / "manifests" / "selected_inputs.csv", index=False)

        # Build project-specific config.
        srow = summary.loc[gse] if gse in summary.index else None
        cfg = {
            "project": {
                "gse": gse,
                "name": f"scanpy_{gse}",
                "root": str(gse_dir.relative_to(ROOT)),
            },
            "inputs": {
                "selected_manifest": str((gse_dir / "manifests" / "selected_inputs.csv").relative_to(ROOT)),
                "matrix_dir": str((DOWNLOADS / gse / "matrix").relative_to(ROOT)) if (DOWNLOADS / gse / "matrix").exists() else "",
                "suppl_dir": str((DOWNLOADS / gse / "suppl").relative_to(ROOT)) if (DOWNLOADS / gse / "suppl").exists() else "",
                "selected_input_count": int(len(grp_sorted)),
            },
            "metadata": {
                "file": "" if srow is None else ("" if pd.isna(srow.get("metadata_file", "")) else str(srow.get("metadata_file", ""))),
                "donor_column": "" if srow is None else ("" if pd.isna(srow.get("donor_column", "")) else str(srow.get("donor_column", ""))),
                "tissue_column": "" if srow is None else ("" if pd.isna(srow.get("tissue_column", "")) else str(srow.get("tissue_column", ""))),
            },
            "rules": {
                "filtered_first": True,
                "raw_only_if_no_filtered": True,
                "exclude_noncoding_early": True,
                "map_to_gene_symbol": True,
                "deduplicate_samples": True,
            },
            "qc": {
                "min_genes_per_cell": 200,
                "min_cells_per_gene": 3,
                "max_mito_pct": 20.0,
            },
            "resources": {
                "max_cpus": 80,
                "max_memory_gb": 500,
                "conda_env": "rapids_sc_py310",
            },
        }
        (gse_dir / "config" / "config.json").write_text(json.dumps(cfg, indent=2) + "\n")

        # Project README
        readme_lines = [
            f"# {gse} Scanpy Project",
            "",
            "This project scaffold was auto-generated from the V4 selected load plan.",
            "",
            "## Inputs",
            f"- Selected manifest: `manifests/selected_inputs.csv`",
            f"- Metadata file: `{cfg['metadata']['file']}`",
            "",
            "## Rules",
            "- Filtered-first loading, raw fallback only when filtered is unavailable.",
            "- Sample-level deduplication via canonical sample key + source priority.",
            "- Non-coding RNA exclusion early in preprocessing.",
            "- Gene IDs standardized to gene symbols before integration.",
            "",
            "## Run context",
            "- Conda env: `rapids_sc_py310`",
            "- Resource limits: up to 80 CPUs, 500 GB RAM",
        ]
        if len(grp_sorted) == 0:
            readme_lines.extend(
                [
                    "",
                    "## Current status",
                    "- No selected expression inputs yet.",
                    "- Check discovery/blocker manifests before loading.",
                ]
            )
        (gse_dir / "README.md").write_text("\n".join(readme_lines) + "\n")

        # Reuse helper scripts when available.
        if template_scripts.exists():
            for py in template_scripts.glob("*.py"):
                dst = gse_dir / "scripts" / py.name
                if not dst.exists():
                    shutil.copy2(py, dst)
        if template_readme.exists():
            dst = gse_dir / "README_template_tcr.md"
            if not dst.exists():
                shutil.copy2(template_readme, dst)

        created += 1

    status_path = OUT / "status" / "gse_status.json"
    status = json.loads(status_path.read_text()) if status_path.exists() else {}
    status["scanpy_projects_generated"] = created
    status["scanpy_projects_root"] = str(projects_root.relative_to(ROOT))
    status["scanpy_projects_note"] = "One project folder per GSE with selected inputs and config."
    status_path.write_text(json.dumps(status, indent=2) + "\n")

    print(f"Scaffold complete: {created} project folders")
    print(f"Root: {projects_root}")


SCANPY_RUNNER_TEMPLATE = """#!/usr/bin/env python3
\"\"\"
Per-GSE Scanpy runner generated by GSE26 V4 scaffold.
\"\"\"

from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path
from typing import Optional, Tuple, List

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io, sparse


def load_config(cfg_path: Path) -> dict:
    return json.loads(cfg_path.read_text())


def find_companion(base: Path, stem: str, kind: str) -> Optional[Path]:
    # Supports multiple patterns:
    #   {prefix}_features.tsv.gz
    #   {prefix}.features.tsv.gz
    #   {prefix}_matrix_1.1.mtx.gz -> {prefix}_features_1.1.tsv.gz
    #   {prefix}.matrix.mtx.gz -> {prefix}.features.tsv.gz
    cands: List[Path] = []
    variants = [stem]
    variants.append(stem.replace(\"_matrix\", f\"_{kind}\"))
    variants.append(stem.replace(\".matrix\", f\".{kind}\"))
    variants = list(dict.fromkeys(variants))
    for v in variants:
        cands.extend(
            [
                base / f\"{v}_{kind}.tsv.gz\",
                base / f\"{v}.{kind}.tsv.gz\",
                base / f\"{v}_{kind}.tsv\",
                base / f\"{v}.{kind}.tsv\",
                base / f\"{v}.tsv.gz\",
                base / f\"{v}.tsv\",
            ]
        )
    cands.extend(
        [
            base / f\"{kind}.tsv.gz\",
            base / f\"{kind}.tsv\",
        ]
    )
    for p in cands:
        if p.exists():
            return p
    return None


def read_tsv_lines(path: Path) -> List[str]:
    op = gzip.open if path.name.endswith(\".gz\") else open
    out = []
    with op(path, \"rt\") as f:
        for line in f:
            out.append(line.rstrip(\"\\n\"))
    return out


def load_mtx(path: Path) -> anndata.AnnData:
    # Handle both genes x cells and cells x genes MTX layouts.
    base = path.parent
    name = path.name
    stem = name.replace(\".matrix.mtx.gz\", \"\").replace(\".mtx.gz\", \"\").replace(\".mtx\", \"\")

    feat = find_companion(base, stem, \"features\") or find_companion(base, stem, \"genes\")
    bar = find_companion(base, stem, \"barcodes\")
    if feat is None or bar is None:
        raise FileNotFoundError(f\"Missing features/barcodes for {path}\")

    m = io.mmread(str(path))
    if not sparse.issparse(m):
        m = sparse.csr_matrix(m)
    else:
        m = m.tocsr()

    features = [x.split(\"\\t\")[1] if \"\\t\" in x and len(x.split(\"\\t\")) > 1 else x.split(\"\\t\")[0] for x in read_tsv_lines(feat)]
    barcodes = [x.split(\"\\t\")[0] for x in read_tsv_lines(bar)]
    n_row, n_col = m.shape
    n_feat = len(features)
    n_bar = len(barcodes)

    if n_row == n_feat and n_col == n_bar:
        x = m.T  # genes x cells -> cells x genes
    elif n_row == n_bar and n_col == n_feat:
        x = m  # already cells x genes
    else:
        raise ValueError(
            f\"Matrix/metadata dimension mismatch for {path}: \"
            f\"matrix={m.shape}, features={n_feat}, barcodes={n_bar}\"
        )

    ad = anndata.AnnData(X=x)
    ad.var_names = pd.Index(features).astype(str)
    ad.obs_names = pd.Index(barcodes).astype(str)
    ad.var_names_make_unique()
    return ad


def load_one(path: Path) -> anndata.AnnData:
    p = str(path).lower()
    if p.endswith(\".h5ad\"):
        ad = sc.read_h5ad(path)
        ad.var_names_make_unique()
        return ad
    if p.endswith(\".h5\"):
        ad = sc.read_10x_h5(path)
        ad.var_names_make_unique()
        return ad
    if p.endswith(\".mtx\") or p.endswith(\".mtx.gz\"):
        return load_mtx(path)
    raise ValueError(f\"Unsupported input format: {path}\")


def remove_noncoding(adata: anndata.AnnData) -> anndata.AnnData:
    patterns = (\"MALAT1\", \"NEAT1\", \"XIST\", \"LINC\", \"MIR\", \"SNOR\", \"SNORA\", \"SNORD\", \"SCARNA\", \"RNU\")
    names = pd.Index(adata.var_names).astype(str)
    keep = np.ones(len(names), dtype=bool)
    for i, g in enumerate(names):
        gu = g.upper()
        if any(p in gu for p in patterns):
            keep[i] = False
    return adata[:, keep].copy()


def read_metadata_table(path: Path) -> Optional[pd.DataFrame]:
    if not path.exists():
        return None
    if path.name.lower().endswith(".h5ad"):
        try:
            ad = anndata.read_h5ad(path, backed="r")
            md = ad.obs.copy()
            md["cell_barcode"] = pd.Index(ad.obs_names).astype(str)
            try:
                ad.file.close()
            except Exception:
                pass
            return md
        except Exception:
            return None
    n = path.name.lower()
    sep = "\t" if n.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz")) else ","
    try:
        return pd.read_csv(path, sep=sep, low_memory=False)
    except Exception:
        return None


def normalize_token(x: str) -> str:
    return "".join(ch for ch in str(x).lower() if ch.isalnum())


def pick_col(cols: List[str], keys: List[str]) -> Optional[str]:
    best = None
    best_score = -1
    for c in cols:
        cl = c.lower().replace(" ", "").replace("-", "_")
        score = sum(1 for k in keys if k in cl)
        if score > best_score:
            best = c
            best_score = score
    return best if best_score > 0 else None


def prepare_metadata(cfg: dict, root: Path) -> Tuple[Optional[pd.DataFrame], Optional[str], Optional[str]]:
    mpath = str(cfg.get("metadata", {}).get("file", "")).strip()
    if not mpath:
        return None, None, None
    p = Path(mpath)
    if not p.is_absolute():
        p = (root.parents[2] / p).resolve()
    md = read_metadata_table(p)
    if md is None or md.empty:
        return None, None, None
    md.columns = [str(c) for c in md.columns]
    bcol = pick_col(md.columns.tolist(), ["barcode", "barcodes", "cellbarcode", "cell_barcode", "cellid", "cell_id"])
    if bcol is None and "Unnamed: 0" in md.columns:
        idx_vals = md["Unnamed: 0"].dropna().astype(str).head(500)
        looks_barcode = idx_vals.str.contains(r"[A-Za-z].*-[0-9]+$|_[A-Za-z0-9]+-[0-9]+$", regex=True).mean() > 0.5 if len(idx_vals) else False
        if looks_barcode:
            bcol = "Unnamed: 0"
    scol = pick_col(md.columns.tolist(), ["sample", "orig.ident", "orig_ident", "library", "batch", "lane"])
    return md, bcol, scol


def merge_metadata_into_obs(
    ad: anndata.AnnData,
    md: Optional[pd.DataFrame],
    barcode_col: Optional[str],
    sample_col: Optional[str],
    sample_value: str,
    sample_map: Optional[dict] = None,
) -> Tuple[anndata.AnnData, int]:
    if md is None:
        return ad, 0

    def canonical_barcode(v: str) -> str:
        s = str(v).strip()
        # Remove concat-added tail like "-0" while keeping native barcode "-1".
        s = re.sub(r"^(.*-[0-9]+)-[0-9]+$", r"\\1", s)
        # Handle prefixes/suffixes separated by underscores.
        if "_" in s:
            toks = [t for t in s.split("_") if t]
            cand = [t for t in toks if re.search(r"[A-Za-z]{6,}.*-[0-9]+$", t)]
            if cand:
                s = max(cand, key=len)
        return s

    # Sample-level fallback when no reliable barcode column exists.
    if barcode_col is None or barcode_col not in md.columns:
        if sample_col and sample_col in md.columns:
            md2 = md.copy()
            md2["__s"] = md2[sample_col].astype(str).map(normalize_token)
            sval = normalize_token(sample_value)
            if sample_map and sample_value in sample_map:
                sval = normalize_token(str(sample_map[sample_value]))
            cand = md2.loc[md2["__s"] == sval]
            if cand.empty and sval:
                cand = md2.loc[md2["__s"].str.contains(sval, regex=False, na=False)]
            if cand.empty and sval:
                cand = md2.loc[md2["__s"].map(lambda x: sval in str(x))]
            if not cand.empty:
                # If multiple rows map to this sample, keep all unique values.
                for c in md2.columns:
                    if c == "__s":
                        continue
                    vals = cand[c].dropna().astype(str).unique().tolist()
                    ad.obs[c] = "|".join(vals) if vals else ""
                return ad, int(ad.n_obs)
        return ad, 0

    obs = ad.obs.copy()
    obs["barcode"] = pd.Index(ad.obs_names).astype(str)
    obs["__k"] = obs["barcode"].astype(str)
    obs["__k_norm"] = obs["__k"].map(canonical_barcode)

    md2 = md.copy()
    md2[barcode_col] = md2[barcode_col].astype(str)
    if sample_col and sample_col in md2.columns:
        sval = normalize_token(sample_value)
        sclean = md2[sample_col].astype(str).map(normalize_token)
        keep = sclean == sval
        if keep.any():
            md2 = md2.loc[keep].copy()
    md2["__k"] = md2[barcode_col].astype(str)
    md2["__k_norm"] = md2["__k"].map(canonical_barcode)
    # Avoid one-to-many expansion during merge when metadata has duplicate barcodes.
    md2 = md2.drop_duplicates(subset=["__k"], keep="first")
    md2_norm = md2.drop_duplicates(subset=["__k_norm"], keep="first")

    meta_cols = [c for c in md2.columns if c not in {"__k", "__k_norm"}]
    m1 = obs.merge(md2[["__k"] + meta_cols], on="__k", how="left", suffixes=("", "_meta"))
    match1 = int(m1[meta_cols[0]].notna().sum()) if meta_cols else 0
    if match1 == 0:
        m1 = obs.merge(md2_norm[["__k_norm"] + meta_cols], on="__k_norm", how="left", suffixes=("", "_meta"))
        match1 = int(m1[meta_cols[0]].notna().sum()) if meta_cols else 0

    m1.index = obs.index.astype(str)
    for c in ("__k", "__k_norm"):
        if c in m1.columns:
            del m1[c]
    ad.obs = m1
    return ad, match1


def run(config_path: Path) -> None:
    cfg = load_config(config_path)
    root = config_path.resolve().parents[1]
    manifest = Path(cfg[\"inputs\"][\"selected_manifest\"])
    if not manifest.is_absolute():
        manifest = (root.parents[2] / manifest).resolve()
    outdir = root / \"outputs\"
    outdir.mkdir(exist_ok=True, parents=True)
    logdir = root / \"logs\"
    logdir.mkdir(exist_ok=True, parents=True)

    sel = pd.read_csv(manifest)
    if sel.empty:
        (logdir / \"run_summary.json\").write_text(json.dumps({\"status\": \"no_inputs\"}, indent=2) + \"\\n\")
        print(\"No selected inputs. Exiting.\")
        return

    md, barcode_col, sample_col = prepare_metadata(cfg, root)
    sample_map = cfg.get("metadata", {}).get("sample_map", {})
    adatas = []
    per_sample = []
    failed_samples = []
    metadata_matched_total = 0
    for _, r in sel.iterrows():
        p = Path(r[\"path\"])
        if not p.is_absolute():
            p = (root.parents[2] / p).resolve()
        sample = str(r.get(\"sample_key\", p.stem))
        try:
            ad = load_one(p)
            ad.obs[\"sample\"] = sample
            ad.obs[\"GSE\"] = cfg[\"project\"][\"gse\"]
            ad, matched = merge_metadata_into_obs(ad, md, barcode_col, sample_col, sample, sample_map)
            metadata_matched_total += int(matched)
            per_sample.append({\"sample\": sample, \"n_obs\": int(ad.n_obs), \"n_vars\": int(ad.n_vars), \"path\": str(p)})
            adatas.append(ad)
        except Exception as e:
            failed_samples.append({\"sample\": sample, \"path\": str(p), \"error\": str(e)})

    if not adatas:
        fail = {\"status\": \"all_samples_failed\", \"failed_samples\": failed_samples}
        (logdir / \"run_summary.json\").write_text(json.dumps(fail, indent=2) + \"\\n\")
        pd.DataFrame(failed_samples).to_csv(outdir / \"failed_samples.csv\", index=False)
        print(json.dumps(fail, indent=2))
        return

    merged = adatas[0] if len(adatas) == 1 else anndata.concat(adatas, join=\"outer\", index_unique=\"-\")
    merged = remove_noncoding(merged)
    merged.var_names_make_unique()

    # Basic QC
    sc.pp.filter_cells(merged, min_genes=int(cfg[\"qc\"][\"min_genes_per_cell\"]))
    sc.pp.filter_genes(merged, min_cells=int(cfg[\"qc\"][\"min_cells_per_gene\"]))
    merged.var[\"mt\"] = pd.Index(merged.var_names).str.upper().str.startswith(\"MT-\")
    sc.pp.calculate_qc_metrics(merged, qc_vars=[\"mt\"], inplace=True)
    merged = merged[merged.obs[\"pct_counts_mt\"] < float(cfg[\"qc\"][\"max_mito_pct\"]), :].copy()

    # Normalize/log
    sc.pp.normalize_total(merged, target_sum=1e4)
    sc.pp.log1p(merged)

    # Ensure metadata columns are HDF5-safe for h5ad writing.
    for c in list(merged.obs.columns):
        s = merged.obs[c]
        if pd.api.types.is_numeric_dtype(s) or pd.api.types.is_bool_dtype(s):
            continue
        merged.obs[c] = s.astype(str).fillna("")

    merged.write_h5ad(outdir / \"scanpy_processed.h5ad\")
    pd.DataFrame(per_sample).to_csv(outdir / \"per_sample_load_stats.csv\", index=False)
    if failed_samples:
        pd.DataFrame(failed_samples).to_csv(outdir / \"failed_samples.csv\", index=False)
    summary = {
        \"gse\": cfg[\"project\"][\"gse\"],
        \"samples_loaded\": len(per_sample),
        \"samples_failed\": len(failed_samples),
        \"metadata_file\": str(cfg.get(\"metadata\", {}).get(\"file\", \"\")),
        \"metadata_barcode_column\": barcode_col if barcode_col else \"\",
        \"metadata_sample_column\": sample_col if sample_col else \"\",
        \"metadata_cells_matched\": int(metadata_matched_total),
        \"cells_final\": int(merged.n_obs),
        \"genes_final\": int(merged.n_vars),
    }
    (logdir / \"run_summary.json\").write_text(json.dumps(summary, indent=2) + \"\\n\")
    print(json.dumps(summary, indent=2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(\"--config\", default=\"../config/config.json\", help=\"Path to project config.json\")
    args = ap.parse_args()
    cfg_path = Path(args.config)
    if not cfg_path.is_absolute():
        cfg_path = (Path(__file__).resolve().parents[2] / cfg_path).resolve()
    run(cfg_path)


if __name__ == \"__main__\":
    main()
"""


def generate_scanpy_files(targets: List[str]) -> None:
    projects_root = OUT / "scanpy_projects"
    created = 0
    for gse in targets:
        gdir = projects_root / gse
        if not gdir.exists():
            continue
        scripts_dir = gdir / "scripts"
        scripts_dir.mkdir(parents=True, exist_ok=True)
        runner = scripts_dir / "run_scanpy_project.py"
        runner.write_text(SCANPY_RUNNER_TEMPLATE)

        launcher = scripts_dir / "launch_scanpy.sh"
        launcher.write_text(
            "#!/usr/bin/env bash\n"
            "set -euo pipefail\n"
            "cd \"$(dirname \"$0\")\"\n"
            "python run_scanpy_project.py --config ../config/config.json\n"
        )
        launcher.chmod(launcher.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

        readme = gdir / "README_SCANPY.md"
        readme.write_text(
            f"# Scanpy Run: {gse}\\n\\n"
            "Run with:\\n\\n"
            "```bash\\n"
            "cd scripts\\n"
            "bash launch_scanpy.sh\\n"
            "```\\n\\n"
            "Outputs:\\n"
            "- `outputs/scanpy_processed.h5ad`\\n"
            "- `outputs/per_sample_load_stats.csv`\\n"
            "- `logs/run_summary.json`\\n"
        )
        created += 1

    status_path = OUT / "status" / "gse_status.json"
    status = json.loads(status_path.read_text()) if status_path.exists() else {}
    status["scanpy_runner_files_generated"] = created
    status_path.write_text(json.dumps(status, indent=2) + "\n")
    print(f"Generated Scanpy runner files for {created} projects")



def run_discovery(targets: List[str]) -> None:
    summaries: List[Dict[str, object]] = []
    details: List[Dict[str, object]] = []

    for g in targets:
        s, d = discover_gse(g)
        summaries.append(s)
        details.extend(d)

    # Build deterministic load plan:
    # 1) use filtered candidates if present for GSE, otherwise raw, otherwise unknown
    # 2) deduplicate by sample key, keeping best source priority
    selected: List[Dict[str, object]] = []
    for g in targets:
        g_rows = [r for r in details if r["GSE"] == g and r["kind"] == "candidate"]
        filtered = [r for r in g_rows if r["matrix_kind"] == "filtered"]
        raw = [r for r in g_rows if r["matrix_kind"] == "raw"]
        unknown = [r for r in g_rows if r["matrix_kind"] not in ("filtered", "raw")]
        chosen_pool = filtered if filtered else (raw if raw else unknown)

        # If both RNA and protein are present in candidate set, keep RNA only.
        modalities = {str(r.get("modality", "unknown")) for r in chosen_pool}
        if "rna" in modalities:
            chosen_pool = [r for r in chosen_pool if r.get("modality") == "rna"]

        best_by_sample: Dict[str, Dict[str, object]] = {}
        for r in chosen_pool:
            skey = sample_key_from_path(str(r["path"]))
            cur = best_by_sample.get(skey)
            if cur is None or source_priority(str(r["path"])) < source_priority(str(cur["path"])):
                rr = dict(r)
                rr["sample_key"] = skey
                rr["source_priority"] = source_priority(str(r["path"]))
                best_by_sample[skey] = rr

        selected.extend(sorted(best_by_sample.values(), key=lambda x: (x["GSE"], x["sample_key"])))

    selected_count_by_gse = {}
    for r in selected:
        selected_count_by_gse[r["GSE"]] = selected_count_by_gse.get(r["GSE"], 0) + 1
    for s in summaries:
        s["selected_samples"] = selected_count_by_gse.get(s["GSE"], 0)

    write_csv(OUT / "manifests" / "gse_discovery_summary.csv", summaries)
    write_csv(OUT / "manifests" / "gse_candidate_expression_files.csv", details)
    write_csv(OUT / "manifests" / "gse_selected_load_plan.csv", selected)

    status = {
        "stage": "discovery",
        "targets": len(targets),
        "present_gse_dirs": sum(1 for x in summaries if x["exists"]),
        "gse_with_metadata": sum(1 for x in summaries if x["metadata_file"]),
        "total_candidates": sum(int(x["candidate_expr_files"]) for x in summaries),
        "total_filtered_candidates": sum(int(x["filtered_candidates"]) for x in summaries),
        "total_raw_candidates": sum(int(x["raw_candidates"]) for x in summaries),
        "total_selected_samples": sum(int(x["selected_samples"]) for x in summaries),
        "timestamp_note": "Generated by gse26_pipeline.py discovery stage",
    }
    (OUT / "status").mkdir(parents=True, exist_ok=True)
    (OUT / "status" / "gse_status.json").write_text(json.dumps(status, indent=2) + "\n")

    print("Discovery complete")
    print(json.dumps(status, indent=2))



def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--stage",
        choices=["discovery", "plan_stats", "scaffold_projects", "generate_scanpy_files"],
        default="discovery",
    )
    args = ap.parse_args()

    targets = load_targets()

    if args.stage == "discovery":
        run_discovery(targets)
    elif args.stage == "plan_stats":
        build_plan_stats()
    elif args.stage == "scaffold_projects":
        scaffold_projects(targets)
    elif args.stage == "generate_scanpy_files":
        generate_scanpy_files(targets)


if __name__ == "__main__":
    main()
