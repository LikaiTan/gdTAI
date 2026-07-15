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
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import anndata as ad
import pandas as pd

ad.settings.allow_write_nullable_strings = True


ROOT = _TNK_PROJECT_ROOT / "analysis_26GSE_V4/scanpy_projects"
TCR_DIR = Path("analysis_output_v3")
OUT_REPORT = _TNK_PROJECT_ROOT / "analysis_26GSE_V4/reports/tcr_mapping_to_h5ad_report.csv"
MAP_FILE = Path("analysis_output_v3/sample_mapping_scRNA_to_TCR.csv")

TCR_KEEP_COLS = [
    "TRA_cdr3",
    "TRA_v",
    "TRA_j",
    "TRA_clone_id",
    "TRB_cdr3",
    "TRB_v",
    "TRB_j",
    "TRB_clone_id",
]


def _clean_str(x: object) -> str:
    s = str(x).strip()
    if s.lower() in {"", "nan", "none", "na"}:
        return ""
    return s


def norm_barcode(x: object) -> str:
    s = _clean_str(x).upper()
    if not s:
        return ""
    # Drop sample suffix patterns like "_CRC01_tissue".
    s = s.split("_")[0]
    m = re.match(r"^([ACGTN]+-\d+)(-\d+)+$", s)
    if m:
        return m.group(1)
    return s


def norm_barcode_core(x: object) -> str:
    s = norm_barcode(x)
    if not s:
        return ""
    # Core sequence before lane/suffix token, e.g. AAAC... from AAAC...-1 or AAAC...-CT8.
    return s.split("-")[0]


def norm_sample(x: object) -> str:
    s = _clean_str(x).lower()
    if not s:
        return ""
    s = s.replace("-", "_")
    s = re.sub(r"__+", "_", s)
    return s


def gsm_id(x: object) -> str:
    m = re.search(r"(gsm\d+)", str(x).lower())
    return m.group(1) if m else ""


def strip_sample_tail(s: str) -> str:
    s = s.lower()
    s = re.sub(r"_10xlane\d+", "", s)
    s = re.sub(r"_lane\d+", "", s)
    s = re.sub(r"(_matrix|_gex|_rna|_tcr)$", "", s)
    s = re.sub(r"_filtered_contig_annotations(_all)?$", "", s)
    s = re.sub(r"(_all|_annotations)$", "", s)
    s = re.sub(r"__+", "_", s)
    return s.strip("_")


def sample_tokens(s: object, gse: str) -> Set[str]:
    raw = norm_sample(s)
    if not raw:
        return set()
    toks: Set[str] = {raw}
    g = gsm_id(raw)
    if g:
        toks.add(g)

    stripped = strip_sample_tail(raw)
    if stripped:
        toks.add(stripped)

    no_gsm = re.sub(r"^gsm\d+_", "", stripped)
    if no_gsm:
        toks.add(no_gsm)

    parts = no_gsm.split("_") if no_gsm else []
    if parts:
        toks.add("_".join(parts[-2:]))
        toks.add(parts[-1])

    m = re.search(r"(\d+\.\d+)$", raw)
    if m:
        toks.add(m.group(1))

    m = re.search(r"\b(b[123])\b", raw)
    if m:
        toks.add(m.group(1))

    if gse == "GSE311112":
        m = re.search(r"_(\d+)_(baseline|3yrs|5yrs)$", raw)
        if m:
            toks.add(f"{m.group(1)}_{m.group(2)}")
        toks.add(raw.replace("_cll_", "_pt_"))
        toks.add(raw.replace("_pt_", "_cll_"))

    if gse == "GSE243572":
        toks.add(re.sub(r"^gsm\d+_", "", raw))

    if gse == "GSE161918":
        m = re.search(r"_(b[123])_", raw)
        if m:
            toks.add(m.group(1))

    return {t for t in toks if t}


def choose_obs_barcode_series(obs: pd.DataFrame) -> pd.Series:
    if "barcode" in obs.columns:
        b = obs["barcode"].astype(str)
        valid = ~b.str.lower().isin(["", "nan", "none", "na"])
        if valid.mean() > 0.5:
            return b
    return pd.Series(obs.index.astype(str), index=obs.index)


def choose_obs_sample_series(obs: pd.DataFrame) -> Optional[pd.Series]:
    for col in ["sample", "sample_id"]:
        if col in obs.columns:
            s = obs[col].astype(str)
            valid = ~s.str.lower().isin(["", "nan", "none", "na"])
            if valid.mean() > 0.05:
                return s
    return None


def build_sample_map_from_rules(gse: str, scrna_samples: Iterable[str], tcr_samples: Iterable[str]) -> Dict[str, str]:
    scrna_samples = list(scrna_samples)
    tcr_samples = list(tcr_samples)
    if not scrna_samples or not tcr_samples:
        return {}

    token_to_tcr: Dict[str, Set[str]] = {}
    for ts in tcr_samples:
        for tk in sample_tokens(ts, gse):
            token_to_tcr.setdefault(tk, set()).add(ts)

    m: Dict[str, str] = {}
    for ss in scrna_samples:
        cands: Set[str] = set()
        for tk in sample_tokens(ss, gse):
            cands.update(token_to_tcr.get(tk, set()))
        if len(cands) == 1:
            m[ss] = next(iter(cands))
    return m


def load_manual_map(gse: str) -> Dict[str, str]:
    if not MAP_FILE.exists():
        return {}
    mdf = pd.read_csv(MAP_FILE, low_memory=False)
    if "GSE" not in mdf.columns:
        return {}
    sub = mdf[mdf["GSE"].astype(str).str.upper() == gse.upper()]
    out: Dict[str, str] = {}
    for _, r in sub.iterrows():
        s = norm_sample(r.get("scRNAseq_sample", ""))
        t = norm_sample(r.get("TCR_sample", ""))
        if not s or not t:
            continue
        if s == "all_cells":
            out["all_cells"] = t
        else:
            out[s] = t
    return out


def aggregate_tcr(df: pd.DataFrame, keys: List[str]) -> pd.DataFrame:
    keep = [c for c in TCR_KEEP_COLS if c in df.columns]
    g = df[keys + keep].copy()

    def first_nonempty(x: pd.Series) -> str:
        vals = [str(v) for v in x if str(v).strip().lower() not in {"", "nan", "none", "na"}]
        return vals[0] if vals else ""

    agg = {c: first_nonempty for c in keep}
    return g.groupby(keys, as_index=False).agg(agg)


def empty_like_obs(index: pd.Index) -> pd.DataFrame:
    out = pd.DataFrame(index=index)
    for c in TCR_KEEP_COLS:
        out[c] = ""
    return out


def main() -> None:
    reports = []
    gse_dirs = sorted([p for p in ROOT.glob("GSE*") if p.is_dir()])
    if len(sys.argv) > 1:
        wanted = {x.strip() for x in sys.argv[1:] if x.strip()}
        gse_dirs = [p for p in gse_dirs if p.name in wanted]

    for gdir in gse_dirs:
        gse = gdir.name
        h5_candidates = [
            gdir / "outputs" / f"{gse}_scanpy_processed.h5ad",
            gdir / "outputs" / "scanpy_processed.h5ad",
        ]
        h5 = next((p for p in h5_candidates if p.exists()), h5_candidates[0])
        tcrf = TCR_DIR / f"tcr_{gse}.csv"

        if not h5.exists():
            reports.append(
                {
                    "GSE": gse,
                    "status": "no_h5ad",
                    "has_tcr_file": False,
                    "cells": 0,
                    "tcr_rows": 0,
                    "mapped_cells": 0,
                    "mapped_fraction": 0.0,
                    "mapping_mode": "",
                }
            )
            continue

        adata = ad.read_h5ad(h5)
        obs = adata.obs.copy()
        n_cells = int(adata.n_obs)

        # Remove legacy TCR columns before rebuilding from source of truth.
        legacy_drop = [c for c in obs.columns if c.startswith("TCR_TRA") or c.startswith("TCR_TRB")]
        legacy_drop.extend([f"TCR_{c}" for c in TCR_KEEP_COLS if f"TCR_{c}" in obs.columns])
        if legacy_drop:
            obs = obs.drop(columns=legacy_drop)

        for c in TCR_KEEP_COLS:
            obs[c] = ""
        obs["TCRseq"] = "no"

        if not tcrf.exists():
            adata.obs = obs
            adata.write_h5ad(h5)
            reports.append(
                {
                    "GSE": gse,
                    "status": "ok_no_tcr_library",
                    "has_tcr_file": False,
                    "cells": n_cells,
                    "tcr_rows": 0,
                    "mapped_cells": 0,
                    "mapped_fraction": 0.0,
                    "mapping_mode": "none",
                }
            )
            continue

        tdf = pd.read_csv(tcrf, low_memory=False)
        if "cell_barcode" not in tdf.columns:
            adata.obs = obs
            adata.write_h5ad(h5)
            reports.append(
                {
                    "GSE": gse,
                    "status": "tcr_missing_barcode",
                    "has_tcr_file": True,
                    "cells": n_cells,
                    "tcr_rows": int(len(tdf)),
                    "mapped_cells": 0,
                    "mapped_fraction": 0.0,
                    "mapping_mode": "none",
                }
            )
            continue

        obs_bc_raw = choose_obs_barcode_series(obs)
        obs_sample_raw = choose_obs_sample_series(obs)
        obs["_bc_norm"] = obs_bc_raw.map(norm_barcode)
        obs["_bc_core"] = obs_bc_raw.map(norm_barcode_core)
        has_sample = obs_sample_raw is not None and obs_sample_raw.astype(str).str.len().gt(0).any()

        tdf["_bc_norm"] = tdf["cell_barcode"].map(norm_barcode)
        tdf["_bc_core"] = tdf["cell_barcode"].map(norm_barcode_core)
        t_has_sample = "sample" in tdf.columns and tdf["sample"].astype(str).str.len().gt(0).any()
        mapping_mode = "barcode_only"

        mapped_mask = pd.Series(False, index=obs.index)
        merge_outputs: List[Tuple[str, pd.DataFrame]] = []

        if has_sample and t_has_sample:
            obs_samples = sorted(set(norm_sample(x) for x in obs_sample_raw.astype(str).unique() if _clean_str(x)))
            tcr_samples = sorted(set(norm_sample(x) for x in tdf["sample"].astype(str).unique() if _clean_str(x)))

            map_manual = load_manual_map(gse)
            map_auto = build_sample_map_from_rules(gse, obs_samples, tcr_samples)

            sample_map: Dict[str, str] = {}
            sample_map.update(map_auto)
            sample_map.update(map_manual)

            if "all_cells" in sample_map:
                t_allowed = {sample_map["all_cells"]}
                tdf2 = tdf[tdf["sample"].map(norm_sample).isin(t_allowed)].copy()
                t_agg = aggregate_tcr(tdf2, ["_bc_norm"])
                m = obs[["_bc_norm"]].reset_index().merge(t_agg, on="_bc_norm", how="left")
                merge_outputs.append(("barcode_only_constrained", m.set_index("index")))
                mapping_mode = "barcode_only_constrained"
            elif sample_map:
                obs["_sample_norm"] = obs_sample_raw.map(norm_sample)
                obs["_sample_tcr"] = obs["_sample_norm"].map(sample_map)
                tdf["_sample_norm"] = tdf["sample"].map(norm_sample)
                t_agg = aggregate_tcr(tdf, ["_sample_norm", "_bc_norm"])
                m_exact = (
                    obs[["_sample_tcr", "_bc_norm"]]
                    .reset_index()
                    .merge(
                        t_agg.rename(columns={"_sample_norm": "_sample_tcr"}),
                        on=["_sample_tcr", "_bc_norm"],
                        how="left",
                    )
                )
                merge_outputs.append(("sample_plus_barcode", m_exact.set_index("index")))
                t_agg_core = aggregate_tcr(tdf, ["_sample_norm", "_bc_core"])
                m_core = (
                    obs[["_sample_tcr", "_bc_core"]]
                    .reset_index()
                    .merge(
                        t_agg_core.rename(columns={"_sample_norm": "_sample_tcr"}),
                        on=["_sample_tcr", "_bc_core"],
                        how="left",
                    )
                )
                merge_outputs.append(("sample_plus_core_barcode", m_core.set_index("index")))
                mapping_mode = "sample_plus_barcode"

        if not merge_outputs:
            t_agg = aggregate_tcr(tdf, ["_bc_norm"])
            m = obs[["_bc_norm"]].reset_index().merge(t_agg, on="_bc_norm", how="left")
            merge_outputs.append(("barcode_only", m.set_index("index")))
            t_agg2 = aggregate_tcr(tdf, ["_bc_core"])
            m2 = obs[["_bc_core"]].reset_index().merge(t_agg2, on="_bc_core", how="left")
            merge_outputs.append(("barcode_core_fallback", m2.set_index("index")))

        modes_used: List[str] = []
        for mode_name, merged in merge_outputs:
            present_cols = [c for c in TCR_KEEP_COLS if c in merged.columns]
            if not present_cols:
                continue
            mode_matched = pd.Series(False, index=obs.index)
            for c in present_cols:
                vals = merged[c].astype(str)
                valid = ~vals.str.lower().isin(["", "nan", "none", "na"])
                use = valid & (~mapped_mask)
                if use.any():
                    obs.loc[use, c] = vals.loc[use].values
                    mode_matched = mode_matched | use
            if mode_matched.any():
                mapped_mask = mapped_mask | mode_matched
                modes_used.append(mode_name)

        obs.loc[mapped_mask, "TCRseq"] = "yes"
        mapped_cells = int(mapped_mask.sum())
        if modes_used:
            mapping_mode = "+".join(modes_used)

        reports.append(
            {
                "GSE": gse,
                "status": "ok",
                "has_tcr_file": True,
                "cells": n_cells,
                "tcr_rows": int(len(tdf)),
                "mapped_cells": mapped_cells,
                "mapped_fraction": round(mapped_cells / n_cells if n_cells else 0.0, 6),
                "mapping_mode": mapping_mode,
            }
        )

        # cleanup helper columns
        for c in ["_bc_norm", "_bc_core", "_sample_norm", "_sample_tcr"]:
            if c in obs.columns:
                del obs[c]
        if "_index" in obs.columns:
            obs = obs.rename(columns={"_index": "obs__index"})
        adata.obs = obs
        tmp_h5 = h5.with_name(h5.name + ".tmp")
        if tmp_h5.exists():
            tmp_h5.unlink()
        adata.write_h5ad(tmp_h5)
        tmp_h5.replace(h5)

    out = pd.DataFrame(reports).sort_values("GSE")
    OUT_REPORT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_REPORT, index=False)
    print(f"Wrote {OUT_REPORT} ({len(out)} rows)")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
