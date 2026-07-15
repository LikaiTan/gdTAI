#!/usr/bin/env python3
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


import re
from pathlib import Path
from typing import Dict, List, Tuple

import anndata as ad
import pandas as pd


ROOT = _TNK_PROJECT_ROOT / "analysis_26GSE_V4/scanpy_projects"
DOWNLOADS = _TNK_PROJECT_ROOT / "downloads"
REPORT = _TNK_PROJECT_ROOT / "analysis_26GSE_V4/reports/gene_symbol_noncoding_cleanup_report.csv"


def read_feature_maps(gse: str) -> Dict[str, str]:
    gdir = DOWNLOADS / gse
    mapping: Dict[str, str] = {}
    if not gdir.exists():
        return mapping

    feat_files = list(gdir.rglob("*features.tsv.gz")) + list(gdir.rglob("*features.tsv"))
    for fp in feat_files:
        try:
            df = pd.read_csv(fp, sep="\t", header=None, dtype=str)
        except Exception:
            continue
        if df.shape[1] < 2:
            continue
        gene_id = df.iloc[:, 0].astype(str).str.strip()
        gene_symbol = df.iloc[:, 1].astype(str).str.strip()
        for gid, sym in zip(gene_id, gene_symbol):
            if not gid or not sym or gid.lower() == "nan" or sym.lower() == "nan":
                continue
            if gid not in mapping:
                mapping[gid] = sym
            gid_nov = re.sub(r"\.\d+$", "", gid)
            if gid_nov not in mapping:
                mapping[gid_nov] = sym
    return mapping


def normalize_var_names_to_symbol(var_names: pd.Index, mapping: Dict[str, str]) -> Tuple[pd.Index, int, int]:
    out = []
    mapped = 0
    unresolved = 0
    for name in var_names.astype(str):
        n = name.strip()
        if re.match(r"^ENSG\d+(\.\d+)?$", n):
            n0 = re.sub(r"\.\d+$", "", n)
            sym = mapping.get(n) or mapping.get(n0)
            if sym and sym.lower() not in {"nan", "none", ""}:
                out.append(sym)
                mapped += 1
            else:
                out.append(n0)
                unresolved += 1
        else:
            out.append(n)
    return pd.Index(out), mapped, unresolved


def is_noncoding_or_nongene(symbol: str) -> bool:
    s = str(symbol).strip()
    if not s:
        return True
    u = s.upper()

    # non-gene antibody/protein tags from feature-barcoding
    if u.endswith("_TOTALSEQC") or u.endswith("_TOTALSEQC"):
        return True

    # obvious non-coding RNA / lncRNA / small RNA / placeholder loci
    nc_patterns = [
        r"^MIR\d",
        r"^MIRLET",
        r"^SNORA",
        r"^SNORD",
        r"^SCARNA",
        r"^RNU\d",
        r"^RN7S",
        r"^RN7SK",
        r"^LINC\d",
        r"^LOC\d+",
        r"^AC\d+",
        r"^AL\d+",
        r"^AP\d+",
        r"^RP\d+",
        r"^MT-RNR",
        r"^MT-T",
        r".*-AS\d*$",
        r"^NEAT1$",
        r"^MALAT1$",
        r"^XIST$",
        r"^TSIX$",
        r"^H19$",
        r"^KCNQ1OT1$",
        r"^MEG3$",
        r"^MIAT$",
        r"^SNHG\d*$",
    ]
    for p in nc_patterns:
        if re.match(p, u):
            return True

    # pseudogene-like symbols
    if re.match(r"^[A-Z0-9]+P\d+$", u):
        return True
    if re.match(r"^RPL\d+P\d+$", u) or re.match(r"^RPS\d+P\d+$", u):
        return True

    return False


def main() -> None:
    rows: List[dict] = []
    for gdir in sorted(ROOT.glob("GSE*")):
        gse = gdir.name
        h5 = gdir / "outputs" / f"{gse}_scanpy_processed.h5ad"
        if not h5.exists():
            rows.append({"GSE": gse, "status": "missing_h5ad"})
            continue

        adata = ad.read_h5ad(h5)
        n_cells_before = int(adata.n_obs)
        n_genes_before = int(adata.n_vars)
        var_before = pd.Index(adata.var_names.astype(str))

        feature_map = read_feature_maps(gse)
        new_names, mapped_n, unresolved_n = normalize_var_names_to_symbol(var_before, feature_map)
        adata.var_names = new_names
        adata.var_names_make_unique()

        keep_mask = ~pd.Index(adata.var_names.astype(str)).map(is_noncoding_or_nongene).to_numpy()
        removed_n = int((~keep_mask).sum())
        adata = adata[:, keep_mask].copy()

        n_genes_after = int(adata.n_vars)
        adata.write_h5ad(h5)

        rows.append(
            {
                "GSE": gse,
                "status": "ok",
                "cells_before": n_cells_before,
                "genes_before": n_genes_before,
                "genes_after": n_genes_after,
                "genes_removed_noncoding_or_nongene": removed_n,
                "ensg_mapped_to_symbol": mapped_n,
                "ensg_unresolved": unresolved_n,
                "feature_map_size": len(feature_map),
            }
        )
        print(
            f"{gse}: genes {n_genes_before} -> {n_genes_after}, "
            f"removed={removed_n}, mapped={mapped_n}, unresolved={unresolved_n}"
        )

    out = pd.DataFrame(rows).sort_values("GSE")
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(REPORT, index=False)
    print(f"Wrote {REPORT}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
