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
from typing import List

import anndata as ad
import pandas as pd


ROOT = _TNK_PROJECT_ROOT / "analysis_26GSE_V4/scanpy_projects"
INFO_CSV = _TNK_PROJECT_ROOT / "configs/datasets/sample_information_final_full.csv"
OUT_REPORT = _TNK_PROJECT_ROOT / "analysis_26GSE_V4/reports/sample_information_merge_report.csv"


def extract_gsm(x: object) -> str:
    m = re.search(r"(GSM\d+)", str(x), flags=re.I)
    return m.group(1).upper() if m else ""


def best_library_series(obs: pd.DataFrame) -> pd.Series:
    candidates = [
        "library_id",
        "gsm",
        "sample",
        "sample_id",
        "sample_name",
        "sample_meta",
        "Sample",
        "orig.ident",
        "orig_ident",
    ]
    best = pd.Series([""] * len(obs), index=obs.index, dtype=object)
    best_n = 0
    for c in candidates:
        if c not in obs.columns:
            continue
        s = obs[c].astype(str).map(extract_gsm)
        n = int((s != "").sum())
        if n > best_n:
            best = s
            best_n = n
    # fallback: try obs_names
    if best_n == 0:
        s = pd.Series(obs.index.astype(str), index=obs.index).map(extract_gsm)
        n = int((s != "").sum())
        if n > best_n:
            best = s
    return best


def main() -> None:
    info = pd.read_csv(INFO_CSV, low_memory=False)
    info["gse"] = info["gse"].astype(str).str.upper()
    info["library_id"] = info["library_id"].astype(str).str.upper()
    # one row per (gse, library_id)
    info = info.drop_duplicates(subset=["gse", "library_id"], keep="first").copy()

    meta_cols = [
        c
        for c in info.columns
        if c not in {"gse", "library_id"}
    ]

    report_rows: List[dict] = []
    for gdir in sorted(ROOT.glob("GSE*")):
        gse = gdir.name.upper()
        h5 = gdir / "outputs" / f"{gse}_scanpy_processed.h5ad"
        if not h5.exists():
            report_rows.append({"GSE": gse, "status": "missing_h5ad"})
            continue

        adata = ad.read_h5ad(h5)
        obs = adata.obs.copy()
        n_cells = int(adata.n_obs)

        lib = best_library_series(obs)
        obs["_library_id_from_obs"] = lib.values

        sub = info[info["gse"] == gse].copy()
        lib_n_info = int(sub["library_id"].nunique()) if len(sub) else 0
        if len(sub) == 0:
            # still write empty flags for traceability
            obs["sample_info_matched"] = False
            adata.obs = obs
            adata.write_h5ad(h5)
            report_rows.append(
                {
                    "GSE": gse,
                    "status": "no_info_rows_for_gse",
                    "cells": n_cells,
                    "cells_with_library_id": int((lib != "").sum()),
                    "cells_matched": 0,
                    "match_fraction": 0.0,
                    "library_ids_in_info": 0,
                }
            )
            continue

        # map library-level metadata to cells
        map_dict = sub.set_index("library_id")[meta_cols].to_dict(orient="index")

        matched = lib.map(lambda x: x in map_dict if x else False)
        obs["sample_info_matched"] = matched.values

        for c in meta_cols:
            colname = f"sample_info_{c}"
            obs[colname] = lib.map(lambda x: map_dict.get(x, {}).get(c, "") if x else "")

        cells_with_lib = int((lib != "").sum())
        cells_matched = int(matched.sum())
        frac = float(cells_matched / n_cells) if n_cells > 0 else 0.0

        adata.obs = obs
        adata.write_h5ad(h5)

        report_rows.append(
            {
                "GSE": gse,
                "status": "ok",
                "cells": n_cells,
                "cells_with_library_id": cells_with_lib,
                "cells_matched": cells_matched,
                "match_fraction": frac,
                "library_ids_in_info": lib_n_info,
            }
        )
        print(
            f"{gse}: matched {cells_matched}/{n_cells} "
            f"({frac:.3%}), lib_in_info={lib_n_info}"
        )

    rep = pd.DataFrame(report_rows).sort_values("GSE")
    OUT_REPORT.parent.mkdir(parents=True, exist_ok=True)
    rep.to_csv(OUT_REPORT, index=False)
    print(f"Wrote {OUT_REPORT}")
    print(rep.to_string(index=False))


if __name__ == "__main__":
    main()
