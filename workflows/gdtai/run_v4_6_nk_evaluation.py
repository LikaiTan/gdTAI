#!/usr/bin/env python3
"""Independent evaluation of gdTAI V4.6's NK exclusion behavior.

Uses the frozen whole-atlas per-cell cache (5,933,312 cells, 40 sources) to
answer: (1) how does V4.6 rule out NK, (2) is the ~0% NK rate real or circular,
(3) how much margin separates NK scores from the operating thresholds, and
(4) what recall does the NK defense cost on TCR-evidenced gdT cells.

Stratification by development exposure:
- never_developed: 9 atlas sources absent from the V4.6 development manifest
  (fully independent NK evaluation)
- development_negative: 8 extension cohorts consumed as negative development
  sources (circular by construction; reported only as a consistency check)
- development_other: remaining development sources (training-adjacent)

Read-only: no H5AD, model, or registry is modified.
Outputs: Integrated_dataset/tables/gdT_prediction/gdtai_v4_6_nk_evaluation/
"""

from __future__ import annotations

import sys as _tnk_sys
from pathlib import Path as _TnkPath

_TNK_PROJECT_ROOT = _TnkPath(__file__).resolve().parents[2]
for _p in (_TNK_PROJECT_ROOT, _TNK_PROJECT_ROOT / "src"):
    if str(_p) not in _tnk_sys.path:
        _tnk_sys.path.insert(0, str(_p))

import json
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CACHE = PROJECT_ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_6_whole_atlas/whole_atlas_predictions.parquet"
MANIFEST = PROJECT_ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_6_development/development_partition_manifest.parquet"
OUT = PROJECT_ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_v4_6_nk_evaluation"
OUT.mkdir(parents=True, exist_ok=True)

THR_HF1 = 0.8943540453910828
THR_HP = 0.9718304872512817
NEAR_MISS_LO = 0.80

NEG_DEV_SOURCES = {"GSE114724", "GSE121636", "GSE159251", "GSE169246",
                   "GSE292700", "GSE294273", "GSE296954", "GSE315928"}


def main() -> None:
    dev_sources = set(pd.read_parquet(MANIFEST, columns=["source_gse_id"]).source_gse_id.astype(str).unique())
    dev_sources.discard("BALF_BLOOD_COPD")  # locked, not fitted

    cols = ["cell_id", "source_gse_id", "v4_6_score", "v2_score",
            "v4_6_highest_f1", "v4_6_high_purity", "v2_high_f1", "v2_high_purity",
            "likely_nk", "author_nk", "conservative_expression_nk", "shared_cytotoxic_ambiguous",
            "has_paired_ab_tcr", "has_any_gd_tcr", "trdc_expressed", "trdv_expressed"]
    df = pd.read_parquet(CACHE, columns=cols)
    src = df.source_gse_id.astype(str)
    df["dev_class"] = np.where(src.isin(NEG_DEV_SOURCES), "development_negative",
                       np.where(src.isin(dev_sources), "development_other", "never_developed"))
    df["trdc_trdv_minus"] = df.trdc_expressed & ~df.trdv_expressed

    strata = {
        "author_nk": df.author_nk,
        "conservative_expression_nk": df.conservative_expression_nk,
        "likely_nk": df.likely_nk,
        "shared_cytotoxic_ambiguous": df.shared_cytotoxic_ambiguous,
        "trdc_plus_trdv_minus": df.trdc_trdv_minus,
        "nk_union": df.author_nk | df.conservative_expression_nk | df.likely_nk,
        # recall-cost strata: cells with TCR evidence that look NK-ish
        "gd_tcr_and_author_nk": df.has_any_gd_tcr & df.author_nk,
        "gd_tcr_any": df.has_any_gd_tcr,
    }

    rows = []
    for sname, mask in strata.items():
        for dc in ("never_developed", "development_negative", "development_other", "ALL"):
            m = mask & (df.dev_class == dc if dc != "ALL" else True)
            n = int(m.sum())
            if n == 0:
                continue
            sc = df.v4_6_score[m]
            v2sc = df.v2_score[m]
            rows.append({
                "stratum": sname, "dev_class": dc, "n": n,
                "v4_6_hf1_callrate": float(df.v4_6_highest_f1[m].mean()),
                "v4_6_hp_callrate": float(df.v4_6_high_purity[m].mean()),
                "v2_hf1_callrate": float(df.v2_high_f1[m].mean()),
                "v2_hp_callrate": float(df.v2_high_purity[m].mean()),
                "v4_6_score_p50": float(sc.quantile(0.5)),
                "v4_6_score_p90": float(sc.quantile(0.9)),
                "v4_6_score_p99": float(sc.quantile(0.99)),
                "v4_6_score_max": float(sc.max()),
                "v4_6_frac_ge_0p5": float((sc >= 0.5).mean()),
                "v4_6_frac_ge_0p8_near_miss": float((sc >= NEAR_MISS_LO).mean()),
                "v4_6_frac_ge_hf1": float((sc >= THR_HF1).mean()),
                "v4_6_frac_ge_hp": float((sc >= THR_HP).mean()),
            })
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "nk_strata_evaluation.csv", index=False)

    # per-source detail for NK union, never-developed sources
    nk_union = strata["nk_union"]
    per_src = []
    for s, g in df[nk_union].groupby(df.source_gse_id):
        dc = g.dev_class.iloc[0]
        per_src.append({"source": s, "dev_class": dc, "n_nk_union": len(g),
                        "v4_6_hf1_callrate": float(g.v4_6_highest_f1.mean()),
                        "v2_hf1_callrate": float(g.v2_high_f1.mean()),
                        "v4_6_score_p99": float(g.v4_6_score.quantile(0.99))})
    pd.DataFrame(per_src).sort_values("dev_class").to_csv(OUT / "nk_union_by_source.csv", index=False)

    # near-miss cells: NK-union cells closest to threshold (worst 20 per dev class)
    nm = df[nk_union & (df.v4_6_score >= 0.5)].sort_values("v4_6_score", ascending=False)
    nm[["cell_id", "source_gse_id", "dev_class", "v4_6_score", "v2_score",
        "v4_6_highest_f1", "v2_high_f1", "trdc_expressed", "trdv_expressed",
        "has_any_gd_tcr"]].head(200).to_csv(OUT / "nk_union_top_scoring_cells.csv", index=False)

    print(res.to_string(index=False))
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
