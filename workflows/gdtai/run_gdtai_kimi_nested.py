#!/usr/bin/env python3
"""Train and nested-evaluate gdTAI-kimi, an experimental two-stage gdT classifier.

gdTAI-kimi is a separately named experimental candidate run outside the frozen
GDTAI V4 precommitted protocol (docs/GDTAI_V4_PRECOMMITTED_PLAN.md). It follows
that protocol's core principles but is NOT the V4 run and cannot be promoted:

- nested leave-one-dataset-out evaluation across HRA005041, GSE144469, and
  BALF_BLOOD_COPD; the held-out dataset is absent from feature filtering,
  fitting, calibration, and threshold selection
- expression-independent TCR-metadata labels computed from raw CDR3 fields
  (sidesteps the broken GSE144469 boolean flags)
- two-tier negatives: abT_primary restricted to sublibraries with demonstrated
  gdTCR capture; censored paired-abT cells are down-weighted training-only
  weak negatives (review finding M1)
- sorted-gdT cells are weight-0.25 training-only weak positives, never evaluated
  as truth (sorted-source call rate is reported as sensitivity-only)
- grouped inner splits (donor, falling back to sample/library); every transform
  (feature filter, Stage-1 probabilities, scaling, calibration, thresholds) is
  fold-local
- frozen candidate grids (elastic-net logistic + histogram gradient boosting),
  sigmoid calibration, guardrail-constrained threshold selection
- frozen V2/V3-R12/V3-R14 profiles are scored on the same held-out cells as
  DESCRIPTIVE references only: they were developed with access to some of these
  cells, so their numbers are optimistically biased in this frame

Documented deviations from the frozen V4 protocol (a separate, later experiment):
GSE144469 expression comes from the integrated-plus6 atlas X rather than the raw
SRR-joined object; elastic-net uses the SGD optimizer (same objective as saga);
tuning uses fold-local training subsamples (<=100k cells) followed by full
refits; isotonic calibration, extension-cohort guardrails, and full-atlas
inference are omitted (sealed cohorts stay sealed).

Outputs (experimental; no registry mutation, no promotion):
- Integrated_dataset/tables/gdT_prediction/gdtai_kimi/
- Integrated_dataset/figures/gdT_prediction/gdtai_kimi/
- Integrated_dataset/logs/gdT_prediction/gdtai_kimi/gdtai_kimi_report.md
- Integrated_dataset/models/gdT_prediction_classifier/gdtai_kimi/ (experimental)
"""

from __future__ import annotations

# TNK_WORKFLOW_BOOTSTRAP
import sys as _tnk_sys
from pathlib import Path as _TnkPath

_TNK_PROJECT_ROOT = _TnkPath(__file__).resolve().parents[2]
for _tnk_path in (
    _TNK_PROJECT_ROOT,
    _TNK_PROJECT_ROOT / "src",
    _TNK_PROJECT_ROOT / "workflows" / "gdtai",
):
    _tnk_value = str(_tnk_path)
    if _tnk_value not in _tnk_sys.path:
        _tnk_sys.path.insert(0, _tnk_value)

import argparse
import json
import logging
import pickle
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    matthews_corrcoef,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedGroupKFold

from compare_frozen_gdtai_profiles import (
    load_profiles,
    normalize_annotation,
    predict_profiles,
)
from run_gdtai_v3_trdc_nk_guard_classifier import (
    FeatureSpec,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
ATLAS_H5AD = PROJECT_ROOT / "high_speed_temp/Integrated_dataset/integrated_plus6.h5ad"
BALF_H5AD = PROJECT_ROOT / "data/datasets/BALF_BLOOD_COPD/processed/current.h5ad"
EVAL_CONFIG = PROJECT_ROOT / "configs/models/gdtai/extension_evaluation.json"

OUTPUT_ROOT = PROJECT_ROOT / "Integrated_dataset"
TABLE_DIR = OUTPUT_ROOT / "tables/gdT_prediction/gdtai_kimi"
FIGURE_DIR = OUTPUT_ROOT / "figures/gdT_prediction/gdtai_kimi"
LOG_DIR = OUTPUT_ROOT / "logs/gdT_prediction/gdtai_kimi"
MODEL_DIR = OUTPUT_ROOT / "models/gdT_prediction_classifier/gdtai_kimi"
for _d in (TABLE_DIR, FIGURE_DIR, LOG_DIR, MODEL_DIR):
    _d.mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("gdtai_kimi")

SEED = 20260807
TARGET_SUM = 10_000.0
CHUNK = 50_000
TUNE_CAP = 100_000
NEG_PRIMARY_CAP = 250_000
NEG_CENSORED_CAP = 150_000
BOOTSTRAP_REPS = 1000

LODO_SOURCES = ["HRA005041", "GSE144469", "BALF_BLOOD_COPD"]
ATLAS_PRIMARY_SOURCES = ["HRA005041", "GSE144469"]
SORTED_WEAK_SOURCES = ["GDT_2020AUG_woCOV", "MalteGDT"]
STRESS_SOURCE = "GSE254249"

# frozen 197-gene universe (docs/GDTAI_V4_PRECOMMITTED_PLAN.md Appendix A + panels)
TCR_GENES_153 = (
    "TRAC TRAJ10 TRAJ12 TRAJ13 TRAJ16 TRAJ18 TRAJ22 TRAJ23 TRAJ27 TRAJ28 TRAJ32 TRAJ34 "
    "TRAJ37 TRAJ38 TRAJ39 TRAJ42 TRAJ43 TRAJ44 TRAJ45 TRAJ46 TRAJ47 TRAJ49 TRAJ5 TRAJ6 "
    "TRAJ9 TRAV1-1 TRAV1-2 TRAV10 TRAV12-1 TRAV12-2 TRAV12-3 TRAV13-1 TRAV13-2 TRAV14DV4 "
    "TRAV15 TRAV16 TRAV17 TRAV18 TRAV19 TRAV2 TRAV20 TRAV21 TRAV22 TRAV23DV6 TRAV24 TRAV25 "
    "TRAV26-1 TRAV26-2 TRAV27 TRAV29DV5 TRAV3 TRAV30 TRAV34 TRAV35 TRAV36DV7 TRAV38-1 "
    "TRAV38-2DV8 TRAV39 TRAV4 TRAV40 TRAV41 TRAV5 TRAV6 TRAV8-1 TRAV8-2 TRAV8-3 TRAV8-4 "
    "TRAV8-5 TRAV8-6 TRAV9-2 TRBC1 TRBC2 TRBJ1-1 TRBJ1-2 TRBJ1-4 TRBJ1-5 TRBJ1-6 TRBJ2-1 "
    "TRBJ2-2 TRBJ2-2P TRBJ2-3 TRBJ2-4 TRBJ2-5 TRBJ2-6 TRBJ2-7 TRBV1 TRBV10-1 TRBV10-2 "
    "TRBV10-3 TRBV11-1 TRBV11-2 TRBV11-3 TRBV12-2 TRBV12-3 TRBV12-4 TRBV12-5 TRBV13 TRBV14 "
    "TRBV15 TRBV16 TRBV18 TRBV19 TRBV2 TRBV20-1 TRBV21-1 TRBV23-1 TRBV24-1 TRBV25-1 TRBV27 "
    "TRBV28 TRBV29-1 TRBV3-1 TRBV30 TRBV4-1 TRBV4-2 TRBV5-1 TRBV5-3 TRBV5-4 TRBV5-5 TRBV5-6 "
    "TRBV6-1 TRBV6-2 TRBV6-4 TRBV6-5 TRBV6-6 TRBV7-2 TRBV7-3 TRBV7-4 TRBV7-5 TRBV7-6 TRBV7-7 "
    "TRBV7-9 TRBV9 TRDC TRDJ1 TRDJ2 TRDJ3 TRDV1 TRDV2 TRDV3 TRGC1 TRGC2 TRGV1 TRGV10 TRGV11 "
    "TRGV2 TRGV3 TRGV4 TRGV5 TRGV5P TRGV7 TRGV8 TRGV9"
).split()
T_LINEAGE = ["CD3D", "CD3E", "CD3G", "CD247", "CD2", "LCK", "LAT", "TRAT1", "CD5", "CD6", "THEMIS", "BCL11B", "TCF7", "LEF1"]
CD4_TREG = ["CD4", "IL7R", "CCR7", "SELL", "LTB", "MAL", "FOXP3", "IL2RA", "CTLA4", "TIGIT"]
NK_CYTO = ["NKG7", "GNLY", "KLRD1", "TYROBP", "FCER1G", "FCGR3A", "KLRC1", "KLRF1", "PRF1", "GZMB", "CTSW", "CCL5", "EOMES", "TBX21"]
T_STATE = ["CD8A", "CD8B", "KLRB1", "ZNF683", "RUNX3", "IKZF2"]
GENES_197 = list(dict.fromkeys([*TCR_GENES_153, *T_LINEAGE, *CD4_TREG, *NK_CYTO, *T_STATE]))
assert len(GENES_197) == 197, len(GENES_197)

STAGE1_GENES = [*T_LINEAGE, *CD4_TREG, *NK_CYTO, *T_STATE, "TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2"]
COMPACT7 = ["TRDC", "TRDV1", "TRDV2", "TRDV3", "TRAC", "TRBC1", "TRBC2"]
CD4_EXCL_PANEL = ["CD4", "IL7R", "CCR7", "SELL", "LTB", "MAL"]
TREG_EXCL_PANEL = ["FOXP3", "IL2RA", "CTLA4", "TIGIT"]

GUARDRAILS = {
    "balanced": {"min_recall": 0.80, "max_abt_fpr": 0.002, "max_nk_fpr": 0.010, "agg": "f1"},
    "high_purity": {"min_recall": 0.70, "max_abt_fpr": 0.001, "max_nk_fpr": 0.005, "agg": "f0.5"},
}
PREVALENCES = [0.001, 0.005, 0.01, 0.02, 0.05]

EN_GRID = [(a, l1) for a in (1e-3, 3e-4, 1e-4, 3e-5) for l1 in (0.0, 0.25, 0.5, 0.75)]
HGB_GRID = [
    {"learning_rate": lr, "max_leaf_nodes": ln, "min_samples_leaf": ms, "l2_regularization": l2}
    for lr in (0.03, 0.07) for ln in (7, 15) for ms in (100, 500) for l2 in (1.0, 10.0)
]
HGB_MAX_ITER = 250

# Round 2 (diagnosis-driven, documented before running):
#  R2-1 Stage-2 trains on grouped NK hard negatives (weight 0.5) — round 1 showed
#       NK leakage, not ab leakage, is the binding constraint.
#  R2-2 TRD-only cells (silver-like: TRD CDR3, not paired gd, no ab, not NK-annotated)
#       join training as weight-0.25 weak positives — round 1 showed BALF positives
#       cluster just below the transferred threshold; broader positive coverage targets this.
#  R2-3 HistGBM capacity raised (leaves {15,31}, max_iter 400) — round 1 grid maxed at
#       leaves 15 / 250 iterations.
#  R2-4 recall-at-fixed-FPR operating points reported for every candidate.
#
# Round 3 (NK label repair, documented before running):
#  R3-1 strict_NK requires a gdTCR-assayed sublibrary (97.9% of HRA005041 strict-NK
#       cells were in sublibraries that could not see gd chains; their "no gdTCR" was
#       censoring, not evidence). Censored NK cells leave fitting and FPR denominators;
#       their call rate is reported as an upper-bound diagnostic only.
#  R3-2 NK-annotated cells with TRDC expression (no TCR contigs) are the biologically
#       ambiguous zone: excluded from training negatives and from the strict-NK FPR
#       denominator; reported separately as the abstain-zone call rate.
#  R3-3 R2-2 reverted (silver weak positives pulled TRD+ NK into the positive class);
#       R2-1 and R2-3 retained.
ROUND2_HGB_GRID = [
    {"learning_rate": lr, "max_leaf_nodes": ln, "min_samples_leaf": ms, "l2_regularization": l2}
    for lr in (0.03, 0.07) for ln in (15, 31) for ms in (100, 500) for l2 in (1.0, 10.0)
]
FPR_OPERATING_POINTS = (0.001, 0.002, 0.005, 0.01)


# --------------------------------------------------------------------------- io helpers

def var_gene_names(path: Path) -> np.ndarray:
    with h5py.File(path, "r") as h:
        var = h["var"]
        for key in ("_index", "index", "gene"):
            if key in var:
                raw = var[key][:]
                break
        else:
            raise KeyError(f"cannot resolve var gene names in {path}")
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in raw])


def obs_col(h: h5py.File, name: str) -> np.ndarray | None:
    if name not in h["obs"]:
        return None
    ds = h["obs"][name]
    if isinstance(ds, h5py.Group):
        if "codes" in ds:  # categorical
            codes = ds["codes"][:]
            cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in ds["categories"][:]])
            out = cats[np.clip(codes, 0, len(cats) - 1)].astype(object)
            out[codes < 0] = ""
            return out
        if "values" in ds:  # nullable-boolean / nullable-integer
            return np.asarray(ds["values"][:])
        raise KeyError(f"unsupported obs group encoding for {name}: {list(ds.keys())}")
    raw = ds[:]
    return raw.astype(str) if raw.dtype.kind == "S" else raw


def nonempty(v: np.ndarray | None, n: int) -> np.ndarray:
    if v is None:
        return np.zeros(n, dtype=bool)
    s = pd.Series(v.astype(object)).fillna("").astype(str).str.strip()
    return ((s != "") & (s.str.upper() != "NA")).to_numpy()


def as_bool(v: np.ndarray | None, n: int) -> np.ndarray:
    if v is None:
        return np.zeros(n, dtype=bool)
    if v.dtype == bool or v.dtype.kind in "iu":
        return v.astype(bool)
    return pd.Series(v.astype(object)).fillna("").astype(str).str.lower().isin(["true", "1", "yes"]).to_numpy()


def extract_features(path: Path, matrix: str, row_idx: np.ndarray, genes: list[str]) -> np.ndarray:
    names = var_gene_names(path)
    gidx = {g: i for i, g in enumerate(names)}
    cols = np.array([gidx.get(g, -1) for g in genes])
    keep = cols >= 0
    if (~keep).any():
        log.warning("%s: %d union genes absent -> zero-filled: %s",
                    path.name, int((~keep).sum()), [g for g, k in zip(genes, keep) if not k][:10])
    out = np.zeros((row_idx.size, len(genes)), dtype=np.float32)
    a = ad.read_h5ad(path, backed="r")
    mat = a.X if matrix == "X" else a.layers[matrix]
    order = np.argsort(row_idx)
    rows_sorted = row_idx[order]
    for start in range(0, rows_sorted.size, CHUNK):
        rows = rows_sorted[start:start + CHUNK]
        block = mat[rows]
        if sparse.issparse(block):
            totals = np.asarray(block.sum(axis=1)).ravel()
            sub = block[:, cols[keep]].toarray()
        else:
            totals = block.sum(axis=1)
            sub = block[:, cols[keep]]
        scale = TARGET_SUM / np.maximum(totals, 1.0)
        chunk = np.zeros((rows.size, len(genes)), dtype=np.float32)
        chunk[:, keep] = np.log1p(sub * scale[:, None])
        chunk[totals <= 0] = 0.0
        out[order[start:start + CHUNK]] = chunk
        if start % (CHUNK * 4) == 0:
            log.info("  %s features: %d/%d cells", path.name, min(start + CHUNK, rows_sorted.size), rows_sorted.size)
    return out


# --------------------------------------------------------------------------- labels

def build_labels(obs: dict[str, np.ndarray], source: str) -> pd.DataFrame:
    n = len(obs["library_id"])
    if source == "BALF_BLOOD_COPD":
        has_ab = as_bool(obs.get("tcr_has_ab"), n)
        any_gd = as_bool(obs.get("tcr_has_gd"), n)
        paired_ab = as_bool(obs.get("tcr_paired_ab"), n)
        paired_gd = as_bool(obs.get("tcr_paired_gd"), n)
        annotation = pd.Series(obs["cell_type"].astype(object)).astype(str).to_numpy()
        assayed = np.ones(n, dtype=bool)  # BALF VDJ assay covers both receptor classes
    else:
        trg = nonempty(obs.get("TRG_cdr3"), n)
        trd = nonempty(obs.get("TRD_cdr3"), n)
        tra = nonempty(obs.get("TRA_cdr3"), n)
        trb = nonempty(obs.get("TRB_cdr3"), n)
        any_gd = trg | trd                                            # CDR3-based, uniform
        paired_gd = (as_bool(obs.get("has_TRG_TRD_paired"), n) | trg) & trd
        has_ab = as_bool(obs.get("has_any_ab_tcr"), n) | tra | trb
        paired_ab = as_bool(obs.get("has_TRA_TRB_paired"), n) | (tra & trb)
        annotation = pd.Series(obs["simple_annotation_plus6"].astype(object)).astype(str).to_numpy()
        lib = pd.Series(obs["library_id"].astype(object)).astype(str)
        assayed = lib.isin(set(lib[trg | trd].unique())).to_numpy()
    df = pd.DataFrame({
        "source": source,
        "library_id": pd.Series(obs["library_id"].astype(object)).astype(str).to_numpy(),
        "sample_id": pd.Series(obs["sample_id"].astype(object)).astype(str).to_numpy(),
        "donor_id": pd.Series(obs["donor_id"].astype(object)).astype(str).replace({"": np.nan}).to_numpy(),
        "annotation": annotation,
        "any_gd": any_gd, "paired_gd": paired_gd, "has_ab": has_ab, "paired_ab": paired_ab,
        "sublib_assayed": assayed,
    })
    df["group"] = df["donor_id"].fillna(df["sample_id"]).fillna(df["library_id"])
    nk_label = "NK" if source == "BALF_BLOOD_COPD" else "NK_cell"
    dual = df.has_ab & df.any_gd
    cls = pd.Series("excluded", index=df.index, dtype=object)
    cls[dual] = "dual_excluded"
    cls[df.paired_gd & ~df.has_ab] = "gdT_primary"
    cls[df.paired_ab & ~df.any_gd & df.sublib_assayed] = "abT_primary"
    cls[df.paired_ab & ~df.any_gd & ~df.sublib_assayed] = "abT_censored_weak"
    cls[(df.annotation == nk_label) & ~df.has_ab & ~df.any_gd] = "strict_NK"
    if source in SORTED_WEAK_SOURCES:
        cls[(~df.has_ab) & ~dual & (cls != "gdT_primary")] = "sorted_gdT_weak"
    if source in LODO_SOURCES:
        # R2-2: TRD-only silver-like cells as weak positives (training only)
        silver = df.any_gd & ~df.paired_gd & ~df.has_ab & (df.annotation != nk_label)
        cls[silver & (cls == "excluded")] = "gdT_silver_weak"
    df["kimi_class"] = cls
    return df


# --------------------------------------------------------------------------- features / models

def derive_features(x197: np.ndarray, names197: list[str]) -> np.ndarray:
    pos = {g: i for i, g in enumerate(names197)}
    feats = []
    for fam in ("TRA", "TRB", "TRG", "TRD"):
        cols = [pos[g] for g in names197 if g.startswith(fam)]
        feats.append((x197[:, cols] > 0).sum(axis=1, keepdims=True).astype(np.float32))
    for panel in (T_LINEAGE, NK_CYTO, CD4_EXCL_PANEL, TREG_EXCL_PANEL):
        cols = [pos[g] for g in panel if g in pos]
        feats.append(x197[:, cols].mean(axis=1, keepdims=True))
    return np.hstack(feats).astype(np.float32)


def fit_en(x, y, w, alpha, l1):
    m = SGDClassifier(loss="log_loss", penalty="elasticnet", alpha=alpha, l1_ratio=l1,
                      max_iter=15, tol=1e-4, random_state=SEED, n_jobs=8)
    m.fit(x, y, sample_weight=w)
    return m


def fit_hgb(x, y, w, cfg):
    m = HistGradientBoostingClassifier(loss="log_loss", max_iter=HGB_MAX_ITER, early_stopping=False,
                                       random_state=SEED, **cfg)
    m.fit(x, y, sample_weight=w)
    return m


def platt(scores, y):
    c = LogisticRegression(C=1.0, solver="lbfgs", max_iter=1000)
    c.fit(scores.reshape(-1, 1), y)
    return c


def apply_platt(c, scores):
    return c.predict_proba(np.asarray(scores).reshape(-1, 1))[:, 1].astype(np.float32)


def fold_weights(cls: pd.Series, dsrc: np.ndarray) -> np.ndarray:
    w = cls.map({"gdT_primary": 1.0, "abT_primary": 1.0, "strict_NK": 0.5,
                 "abT_censored_weak": 0.25, "sorted_gdT_weak": 0.25,
                 "gdT_silver_weak": 0.25}).fillna(0.0).to_numpy(dtype=np.float64)
    ybin = cls.isin(["gdT_primary", "sorted_gdT_weak", "gdT_silver_weak"]).astype(int).to_numpy()
    for s in np.unique(dsrc):
        for c in (0, 1):
            m = (dsrc == s) & (ybin == c)
            if m.sum():
                w[m] *= 1.0 / m.sum()
    w *= len(w) / w.sum()
    return w.astype(np.float32)


def choose_threshold(scores, y, abt, nk, dsrc, mode):
    """Guardrail-constrained threshold search.

    Returns (threshold, info). If no threshold satisfies the mode's guardrails,
    returns a best-effort threshold (FPR constraints only, then objective) with
    info["guardrails_met"] = False so the failure is reported, not hidden.
    """
    g = GUARDRAILS[mode]
    grid = np.unique(np.quantile(scores, np.linspace(0, 1, 2000)))
    best, best_effort = None, None
    for t in grid:
        pred = scores >= t
        recs, f1s, f05s, abs_, nks = [], [], [], [], []
        for s in np.unique(dsrc):
            m = dsrc == s
            yy, pp = y[m], pred[m]
            if yy.sum() == 0:
                continue
            tp = ((pp == 1) & (yy == 1)).sum(); fn = ((pp == 0) & (yy == 1)).sum()
            fp = ((pp == 1) & (yy == 0)).sum()
            rec = tp / max(tp + fn, 1)
            prec = tp / max(tp + fp, 1)
            recs.append(rec)
            f1s.append(2 * prec * rec / max(prec + rec, 1e-12))
            f05s.append(1.25 * prec * rec / max(0.25 * prec + rec, 1e-12))
            ab = m & abt; nkm = m & nk
            abs_.append(pred[ab].mean() if ab.sum() else np.nan)
            nks.append(pred[nkm].mean() if nkm.sum() else np.nan)
        if not recs:
            continue
        m_rec, m_ab, m_nk = np.mean(recs), np.nanmean(abs_), np.nanmean(nks)
        obj = np.mean(f1s) if g["agg"] == "f1" else np.mean(f05s)
        info = {"macro_recall": float(m_rec), "macro_abt_fpr": float(m_ab),
                "macro_nk_fpr": float(m_nk) if not np.isnan(m_nk) else np.nan,
                "objective": float(obj), "mode": mode}
        feasible = (m_rec >= g["min_recall"] and m_ab <= g["max_abt_fpr"]
                    and (np.isnan(m_nk) or m_nk <= g["max_nk_fpr"]))
        info["guardrails_met"] = bool(feasible)
        key = (obj, -(m_nk if not np.isnan(m_nk) else 0.0), -m_ab)
        if feasible and (best is None or key > best[0]):
            best = (key, float(t), info)
        fpr_ok = m_ab <= g["max_abt_fpr"] and (np.isnan(m_nk) or m_nk <= g["max_nk_fpr"])
        if fpr_ok and (best_effort is None or key > best_effort[0]):
            best_effort = (key, float(t), {**info, "guardrails_met": False})
    if best:
        return best[1], best[2]
    if best_effort:
        return best_effort[1], best_effort[2]
    return None, {"mode": mode, "failed": True, "guardrails_met": False}


# --------------------------------------------------------------------------- metrics

def metric_block(y, score, pred):
    tp = int(((pred == 1) & (y == 1)).sum()); fp = int(((pred == 1) & (y == 0)).sum())
    fn = int(((pred == 0) & (y == 1)).sum()); tn = int(((pred == 0) & (y == 0)).sum())
    prec = tp / max(tp + fp, 1); rec = tp / max(tp + fn, 1); spec = tn / max(tn + fp, 1)
    f1 = 2 * prec * rec / max(prec + rec, 1e-12)
    f05 = 1.25 * prec * rec / max(0.25 * prec + rec, 1e-12)
    out = {"tp": tp, "fp": fp, "fn": fn, "tn": tn, "precision": prec, "recall": rec,
           "specificity": spec, "f1": f1, "f0.5": f05,
           "balanced_acc": balanced_accuracy_score(y, pred) if len(set(y.tolist())) > 1 else np.nan,
           "mcc": matthews_corrcoef(y, pred) if len(set(y.tolist())) > 1 and len(set(pred.tolist())) > 1 else np.nan}
    if len(set(y.tolist())) > 1:
        out["roc_auc"] = roc_auc_score(y, score)
        out["pr_auc"] = average_precision_score(y, score)
        out["brier"] = brier_score_loss(y, np.clip(score, 0, 1))
    return out


def ece(y, p, bins=10):
    edges = np.unique(np.quantile(p, np.linspace(0, 1, bins + 1)))
    idx = np.clip(np.digitize(p, edges[1:-1]), 0, len(edges) - 2)
    return float(sum(m.mean() * abs(y[m].mean() - p[m].mean())
                     for m in (idx == b for b in range(len(edges) - 1)) if m.sum()))


def assert_grouped(a, b, groups):
    assert not (set(groups[a]) & set(groups[b])), "group crosses fold boundary"


# --------------------------------------------------------------------------- main

def main() -> None:
    global HGB_GRID, HGB_MAX_ITER
    parser = argparse.ArgumentParser()
    parser.add_argument("--round", type=int, default=3, choices=[1, 2, 3],
                        help="experiment round; 2 adds NK hard negatives + silver weak positives; "
                             "3 repairs NK labels (assayed-sublibrary strict NK, TRDC+ NK ambiguity "
                             "split) and reverts silver weak positives")
    args = parser.parse_args()
    ROUND = args.round
    SUF = "" if ROUND == 1 else f"_r{ROUND}"
    if ROUND >= 2:
        HGB_GRID = ROUND2_HGB_GRID
        HGB_MAX_ITER = 400
    log.info("gdTAI-kimi round %d (output suffix '%s')", ROUND, SUF)
    t0 = time.time()
    config = json.loads(EVAL_CONFIG.read_text())
    profiles = load_profiles(config)
    model_genes = []
    for p in profiles:
        payload = p.payload["base_model"] if p.model_id == "gdtai_v2" else p.payload
        model_genes.extend([str(g) for g in payload["gene_names"]])
    union_genes = list(dict.fromkeys([*GENES_197, *model_genes]))
    pos_u = {g: i for i, g in enumerate(union_genes)}
    log.info("union feature universe: %d genes", len(union_genes))

    # ---------------- load data
    need = ["source_gse_id", "library_id", "sample_id", "donor_id", "simple_annotation_plus6",
            "TRG_cdr3", "TRD_cdr3", "TRA_cdr3", "TRB_cdr3",
            "has_TRA_TRB_paired", "has_any_ab_tcr", "has_TRG_TRD_paired", "phase4_trd_minus_trab"]
    with h5py.File(ATLAS_H5AD, "r") as h:
        obs_all = {c: obs_col(h, c) for c in need if c in h["obs"]}
        log.warning("atlas obs missing (treated absent): %s", [c for c in need if c not in h["obs"]])
    src_arr = pd.Series(obs_all["source_gse_id"].astype(object)).astype(str)
    wanted = ATLAS_PRIMARY_SOURCES + SORTED_WEAK_SOURCES + [STRESS_SOURCE]
    row_sel = np.flatnonzero(src_arr.isin(wanted).to_numpy())

    sources: dict[str, tuple[np.ndarray, pd.DataFrame, np.ndarray]] = {}
    for s in wanted:
        rows = row_sel[(src_arr.iloc[row_sel] == s).to_numpy()]
        lab = build_labels({c: v[rows] for c, v in obs_all.items()}, s)
        p4 = pd.to_numeric(pd.Series(obs_all["phase4_trd_minus_trab"][rows]), errors="coerce").fillna(0).to_numpy(np.float32)
        log.info("extracting %s (%d cells)", s, rows.size)
        x = extract_features(ATLAS_H5AD, "X", rows, union_genes)
        sources[s] = (x, lab, p4)
        if ROUND == 3:
            nk0 = lab.kimi_class == "strict_NK"
            trdc_expr = x[:, list(union_genes).index("TRDC")] > 0
            lab.loc[nk0 & ~lab.sublib_assayed, "kimi_class"] = "NK_censored"
            lab.loc[nk0 & lab.sublib_assayed & trdc_expr, "kimi_class"] = "NK_TRDC_ambiguous"
            log.info("%s round-3 NK repair: strict %d -> clean %d, censored %d, TRDC-ambiguous %d",
                     s, int(nk0.sum()), int((lab.kimi_class == "strict_NK").sum()),
                     int((lab.kimi_class == "NK_censored").sum()),
                     int((lab.kimi_class == "NK_TRDC_ambiguous").sum()))
        log.info("%s labels: %s", s, lab.kimi_class.value_counts().to_dict())
        if s in ATLAS_PRIMARY_SOURCES:
            log.info("%s gdTCR-assayed sublibraries: %d of %d", s,
                     lab.loc[lab.sublib_assayed, "library_id"].nunique(), lab.library_id.nunique())

    with h5py.File(BALF_H5AD, "r") as h:
        balf_obs = {c: obs_col(h, c) for c in
                    ["library_id", "sample_id", "donor_id", "cell_type", "tcr_has_ab", "tcr_has_gd",
                     "tcr_paired_ab", "tcr_paired_gd", "phase4_trd_minus_trab"] if c in h["obs"]}
    n_balf = len(balf_obs["library_id"])
    lab_b = build_labels(balf_obs, "BALF_BLOOD_COPD")
    p4_b = pd.to_numeric(pd.Series(balf_obs["phase4_trd_minus_trab"]), errors="coerce").fillna(0).to_numpy(np.float32)
    x_b = extract_features(BALF_H5AD, "counts", np.arange(n_balf), union_genes)
    sources["BALF_BLOOD_COPD"] = (x_b, lab_b, p4_b)
    if ROUND == 3:
        nk0 = lab_b.kimi_class == "strict_NK"
        trdc_expr = x_b[:, list(union_genes).index("TRDC")] > 0
        lab_b.loc[nk0 & trdc_expr, "kimi_class"] = "NK_TRDC_ambiguous"
        log.info("BALF round-3 NK repair: strict %d -> clean %d, TRDC-ambiguous %d",
                 int(nk0.sum()), int((lab_b.kimi_class == "strict_NK").sum()),
                 int((lab_b.kimi_class == "NK_TRDC_ambiguous").sum()))
    log.info("BALF labels: %s", lab_b.kimi_class.value_counts().to_dict())

    pd.DataFrame([{"source": s, **lab.kimi_class.value_counts().to_dict()} for s, (_, lab, _) in sources.items()]
                 ).to_csv(TABLE_DIR / f"label_inventory{SUF}.csv", index=False)

    g197 = [pos_u[g] for g in GENES_197]
    names197 = GENES_197
    s1_cols = [pos_u[g] for g in STAGE1_GENES]
    c7_cols = [pos_u[g] for g in COMPACT7]
    tcr_cols = [pos_u[g] for g in TCR_GENES_153]

    results, test_scores, fold_models = [], [], {}

    for held_out in LODO_SOURCES:
        log.info("=== outer fold: held out %s (%.1f min) ===", held_out, (time.time() - t0) / 60)
        train_sources = [s for s in LODO_SOURCES if s != held_out]

        frames, xs = [], []
        keep_classes = ["gdT_primary", "abT_primary", "abT_censored_weak", "sorted_gdT_weak", "strict_NK"]
        if ROUND == 2:
            keep_classes.append("gdT_silver_weak")
        for s in train_sources + SORTED_WEAK_SOURCES:
            x_s, lab_s, _ = sources[s]
            keep = lab_s.kimi_class.isin(keep_classes).to_numpy()
            frames.append(lab_s.loc[keep, ["source", "group", "kimi_class"]])
            xs.append(x_s[keep])
        tr = pd.concat(frames, ignore_index=True)
        xtr = np.vstack(xs)
        cls_tr = tr.kimi_class
        groups = tr.group.astype(str).to_numpy()
        dsrc = tr.source.to_numpy()
        is_nk = (cls_tr == "strict_NK").to_numpy()
        fit_mask = np.ones(len(tr), dtype=bool) if ROUND == 2 else ~is_nk
        pos_classes = ["gdT_primary", "sorted_gdT_weak"] + (["gdT_silver_weak"] if ROUND == 2 else [])
        y_all = cls_tr.isin(pos_classes).astype(np.int8).to_numpy()

        # ---- Stage 1: soft T-lineage gate
        s1_pos = cls_tr.isin(["gdT_primary", "abT_primary"]).to_numpy()
        s1_keep = s1_pos | is_nk
        xs1, ys1, gs1 = xtr[:, s1_cols][s1_keep], s1_pos[s1_keep].astype(np.int8), groups[s1_keep]
        ws1 = np.where(ys1 == 1, 1.0, 0.5).astype(np.float32)
        sgkf = StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=SEED)
        best_auc, best_alpha, s1_oof = -1, None, None
        for alpha in (1e-3, 1e-4, 1e-5):
            oof = np.zeros(ys1.size, dtype=np.float32)
            for tri, tei in sgkf.split(xs1, ys1, gs1):
                assert_grouped(tri, tei, gs1)
                oof[tei] = fit_en(xs1[tri], ys1[tri], ws1[tri], alpha, 0.35).predict_proba(xs1[tei])[:, 1]
            auc = roc_auc_score(ys1, oof)
            if auc > best_auc:
                best_auc, best_alpha, s1_oof = auc, alpha, oof
        log.info("stage1 best alpha %g (inner-OOF AUC %.4f)", best_alpha, best_auc)
        s1_full = fit_en(xs1, ys1, ws1, best_alpha, 0.35)
        # fold-local stage-1 probability for every outer-train cell
        s1_prob_tr = np.zeros(xtr.shape[0], dtype=np.float32)
        for tri, tei in StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=SEED + 1).split(xtr, y_all, groups):
            assert_grouped(tri, tei, groups)
            keep_tri = s1_keep[tri]
            m = fit_en(xtr[tri[keep_tri]][:, s1_cols], s1_pos[tri[keep_tri]].astype(np.int8),
                       np.where(s1_pos[tri[keep_tri]], 1.0, 0.5).astype(np.float32), best_alpha, 0.35)
            s1_prob_tr[tei] = m.predict_proba(xtr[tei][:, s1_cols])[:, 1]
        gd1 = s1_keep & cls_tr.eq("gdT_primary").to_numpy()
        ab1 = s1_keep & cls_tr.eq("abT_primary").to_numpy()
        grid1 = np.quantile(s1_oof, np.linspace(0, 1, 2001))
        feas = [t for t in grid1 if (s1_oof[gd1[s1_keep]] >= t).mean() >= 0.99 and (s1_oof[ab1[s1_keep]] >= t).mean() >= 0.98]
        s1_thr = float(max(feas)) if feas else float(grid1[0])
        log.info("stage1 soft-gate threshold %.4f; NK rejection %.3f (soft: all cells continue to stage 2)",
                 s1_thr, float((s1_oof[~s1_pos[s1_keep]] < s1_thr).mean()))

        # ---- Stage 2 features
        x2 = np.hstack([xtr[:, g197], derive_features(xtr[:, g197], names197), s1_prob_tr[:, None]]).astype(np.float32)
        w2 = fold_weights(cls_tr, dsrc)
        w2[~fit_mask] = 0.0
        primary_mask = cls_tr.isin(["gdT_primary", "abT_primary"]).to_numpy()
        abt_mask = (cls_tr == "abT_primary").to_numpy()

        def oof_scores(kind, param):
            oof = np.full(y_all.size, np.nan, dtype=np.float32)
            sg = StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=SEED + 2)
            for tri, tei in sg.split(x2, y_all, groups):
                assert_grouped(tri, tei, groups)
                tri_fit = tri[fit_mask[tri]]
                rng = np.random.default_rng(SEED + 3)
                if tri_fit.size > TUNE_CAP:
                    tri_fit = np.sort(rng.choice(tri_fit, TUNE_CAP, replace=False))
                if kind == "en":
                    m = fit_en(x2[tri_fit], y_all[tri_fit], w2[tri_fit], param[0], param[1])
                else:
                    m = fit_hgb(x2[tri_fit], y_all[tri_fit], w2[tri_fit], param)
                oof[tei] = m.predict_proba(x2[tei])[:, 1]
            return oof

        candidates = {}
        log.info("tuning %d EN + %d HGB configs", len(EN_GRID), len(HGB_GRID))
        for alpha, l1 in EN_GRID:
            candidates[f"en_a{alpha:g}_l{l1:g}"] = ("en", (alpha, l1), oof_scores("en", (alpha, l1)))
        for i, cfg in enumerate(HGB_GRID):
            candidates[f"hgb_{i}"] = ("hgb", cfg, oof_scores("hgb", cfg))
        comp = {
            "C2_compact7": (xtr[:, c7_cols], None),
            "C3_tcr_en": (xtr[:, tcr_cols], None),
        }
        for cname in comp:
            oof = np.full(y_all.size, np.nan, dtype=np.float32)
            sg = StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=SEED + 2)
            xc = comp[cname][0]
            for tri, tei in sg.split(xc, y_all, groups):
                assert_grouped(tri, tei, groups)
                tri_fit = tri[fit_mask[tri]]
                oof[tei] = fit_en(xc[tri_fit], y_all[tri_fit], w2[tri_fit], 1e-4, 0.35).predict_proba(xc[tei])[:, 1]
            comp[cname] = (xc, oof)

        diag_rows = []

        def select_family(kind_filter, mode):
            best = None
            for key, (kind, param, oof) in candidates.items():
                if kind != kind_filter:
                    continue
                cal = platt(oof[primary_mask], y_all[primary_mask])
                pcal = apply_platt(cal, oof)
                t, info = choose_threshold(pcal[primary_mask], y_all[primary_mask], abt_mask[primary_mask],
                                           np.zeros(primary_mask.sum(), bool), dsrc[primary_mask], mode)
                if t is None:
                    continue
                # NK guardrail scored on NK OOF (never fitted)
                nk_fpr = float((apply_platt(cal, oof[is_nk]) >= t).mean()) if is_nk.sum() else np.nan
                nk_ok = np.isnan(nk_fpr) or nk_fpr <= GUARDRAILS[mode]["max_nk_fpr"]
                guardrails_met = bool(info.get("guardrails_met", False)) and bool(nk_ok)
                diag_rows.append({"held_out": held_out, "family": kind_filter, "config": key,
                                  "threshold": t, "nk_fpr_oof": nk_fpr, "guardrails_met": guardrails_met,
                                  **{k: v for k, v in info.items() if k != "guardrails_met"}})
                key_ = (guardrails_met, info["objective"])
                if best is None or key_ > best[0]:
                    best = (key_, key, kind, param, oof, cal, t, {**info, "nk_fpr_oof": nk_fpr,
                                                                  "guardrails_met": guardrails_met})
            return best

        chosen = {}
        for fam, kf in (("kimi_en", "en"), ("kimi_hgb", "hgb")):
            best = select_family(kf, "balanced")
            if best is None:
                log.warning("%s: no config passed balanced guardrails", fam)
                continue
            _, key, kind, param, oof, cal, t_bal, info = best
            t_hp, info_hp = choose_threshold(apply_platt(cal, oof)[primary_mask], y_all[primary_mask],
                                             abt_mask[primary_mask], np.zeros(primary_mask.sum(), bool),
                                             dsrc[primary_mask], "high_purity")
            log.info("%s -> %s balanced thr %.4f guardrails_met=%s (inner %s); hp thr %s",
                     fam, key, t_bal, info.get("guardrails_met"), info, t_hp)
            # full refit with negative caps
            rng = np.random.default_rng(SEED + 4)
            pos_idx = np.flatnonzero(fit_mask & (y_all == 1))
            neg_prim = np.flatnonzero(cls_tr.eq("abT_primary").to_numpy())
            neg_cens = np.flatnonzero(cls_tr.eq("abT_censored_weak").to_numpy())
            neg_nk = np.flatnonzero(is_nk) if ROUND == 2 else np.array([], dtype=int)
            if neg_prim.size > NEG_PRIMARY_CAP:
                neg_prim = np.sort(rng.choice(neg_prim, NEG_PRIMARY_CAP, replace=False))
            if neg_cens.size > NEG_CENSORED_CAP:
                neg_cens = np.sort(rng.choice(neg_cens, NEG_CENSORED_CAP, replace=False))
            use = np.sort(np.concatenate([pos_idx, neg_prim, neg_cens, neg_nk]))
            m = fit_en(x2[use], y_all[use], w2[use], *param) if kind == "en" else fit_hgb(x2[use], y_all[use], w2[use], param)
            chosen[fam] = {"model": m, "cal": cal, "t_bal": t_bal, "t_hp": t_hp, "key": key,
                           "info": info, "oof": oof}

        pd.DataFrame(diag_rows).to_csv(TABLE_DIR / f"inner_candidate_diagnostics_holdout_{held_out}{SUF}.csv", index=False)

        # comparators: full refit + calibration + threshold on same protocol
        comp_fitted = {}
        for cname, (xc, oof) in comp.items():
            cal = platt(oof[primary_mask], y_all[primary_mask])
            t_c, _ = choose_threshold(apply_platt(cal, oof)[primary_mask], y_all[primary_mask],
                                      abt_mask[primary_mask], np.zeros(primary_mask.sum(), bool),
                                      dsrc[primary_mask], "balanced")
            if t_c is None:
                continue
            m = fit_en(xc[fit_mask], y_all[fit_mask], w2[fit_mask], 1e-4, 0.35)
            comp_fitted[cname] = {"model": m, "cal": cal, "t": t_c}
        # C1 phase4 delta threshold
        c1s, c1y, c1a, c1d = [], [], [], []
        for s in train_sources:
            _, lab_s, p4 = sources[s]
            pm = lab_s.kimi_class.isin(["gdT_primary", "abT_primary"]).to_numpy()
            c1s.append(p4[pm]); c1y.append((lab_s.kimi_class == "gdT_primary").to_numpy(np.int8)[pm])
            c1a.append((lab_s.kimi_class == "abT_primary").to_numpy()[pm]); c1d.append(np.repeat(s, pm.sum()))
        t_c1, _ = choose_threshold(np.concatenate(c1s), np.concatenate(c1y), np.concatenate(c1a),
                                   np.zeros(sum(len(v) for v in c1y), bool), np.concatenate(c1d), "balanced")

        # ---- held-out evaluation
        x_te, lab_te, p4_te = sources[held_out]
        cls_te = lab_te.kimi_class
        eval_mask = cls_te.isin(["gdT_primary", "abT_primary"]).to_numpy()
        y_te = (cls_te == "gdT_primary").to_numpy(np.int8)
        abt_te = (cls_te == "abT_primary").to_numpy()
        nk_te = (cls_te == "strict_NK").to_numpy()
        nk_cens_te = (cls_te == "NK_censored").to_numpy()
        nk_trdc_te = (cls_te == "NK_TRDC_ambiguous").to_numpy()
        g_te = lab_te.group.astype(str).to_numpy()
        keep_rows = eval_mask | nk_te
        s1_prob_te = s1_full.predict_proba(x_te[:, s1_cols])[:, 1].astype(np.float32)
        x2_te = np.hstack([x_te[:, g197], derive_features(x_te[:, g197], names197), s1_prob_te[:, None]]).astype(np.float32)

        def record(name, score_te, thr):
            pred = (score_te >= thr).astype(np.int8)
            row = metric_block(y_te[eval_mask], score_te[eval_mask], pred[eval_mask])
            row["ece"] = ece(y_te[eval_mask], np.clip(score_te[eval_mask], 0, 1))
            row["abt_fpr"] = float(pred[abt_te].mean()) if abt_te.sum() else np.nan
            row["nk_fpr"] = float(pred[nk_te].mean()) if nk_te.sum() else np.nan
            row["nk_censored_callrate"] = float(pred[nk_cens_te].mean()) if nk_cens_te.sum() else np.nan
            row["nk_trdc_amb_callrate"] = float(pred[nk_trdc_te].mean()) if nk_trdc_te.sum() else np.nan
            for pv in PREVALENCES:
                row[f"ppv@{pv}"] = row["recall"] * pv / max(row["recall"] * pv + row["abt_fpr"] * (1 - pv), 1e-12)
            row.update({"held_out": held_out, "candidate": name, "threshold": float(thr)})
            results.append(row)
            test_scores.append(pd.DataFrame({"held_out": held_out, "candidate": name,
                                             "group": g_te[keep_rows],
                                             "y": y_te[keep_rows], "abt": abt_te[keep_rows],
                                             "nk": nk_te[keep_rows], "in_eval": eval_mask[keep_rows],
                                             "score": score_te[keep_rows], "pred": pred[keep_rows]}))
            return pred

        for fam, pack in chosen.items():
            score_te = apply_platt(pack["cal"], pack["model"].predict_proba(x2_te)[:, 1])
            pred_bal = record(f"{fam}_balanced", score_te, pack["t_bal"])
            results[-1]["inner_guardrails_met"] = pack["info"].get("guardrails_met")
            if pack["t_hp"] is not None:
                record(f"{fam}_highpurity", score_te, pack["t_hp"])
            if ROUND == 2:
                # R2-4: recall at fixed inner-OOF abT-FPR operating points
                pcal_tr = apply_platt(pack["cal"], pack["oof"])
                for fpr_t in FPR_OPERATING_POINTS:
                    grid = np.unique(np.quantile(pcal_tr[primary_mask], np.linspace(0, 1, 2000)))
                    feas = []
                    for tt in grid:
                        fpr = (pcal_tr[abt_mask] >= tt).mean()
                        if fpr <= fpr_t:
                            rec = ((pcal_tr[primary_mask & (y_all == 1)] >= tt).mean())
                            feas.append((rec, tt))
                    if not feas:
                        continue
                    rec_in, tt = max(feas)
                    row = {"held_out": held_out, "candidate": f"{fam}_fpr{fpr_t:g}",
                           "threshold": float(tt), "inner_recall_at_fpr": float(rec_in),
                           "recall": float((score_te[eval_mask & (y_te == 1)] >= tt).mean()),
                           "abt_fpr": float((score_te[abt_te] >= tt).mean()),
                           "nk_fpr": float((score_te[nk_te] >= tt).mean()) if nk_te.sum() else np.nan}
                    results.append(row)
            cd4m = np.mean([x_te[:, pos_u[g]] for g in CD4_EXCL_PANEL], axis=0)
            cd4det = np.sum([x_te[:, pos_u[g]] > 0 for g in CD4_EXCL_PANEL[1:]], axis=0)
            tregm = np.mean([x_te[:, pos_u[g]] for g in TREG_EXCL_PANEL], axis=0)
            tregdet = np.sum([x_te[:, pos_u[g]] > 0 for g in TREG_EXCL_PANEL[1:]], axis=0)
            excl = ((x_te[:, pos_u["CD4"]] > 0) & (cd4det >= 2) & (cd4m >= 0.5)) | \
                   ((x_te[:, pos_u["FOXP3"]] > 0) & (tregdet >= 1) & (tregm >= 0.5))
            pred_x = pred_bal.copy(); pred_x[excl] = 0
            row = metric_block(y_te[eval_mask], score_te[eval_mask], pred_x[eval_mask])
            row.update({"held_out": held_out, "candidate": f"{fam}_balanced_cd4treg_excl",
                        "threshold": pack["t_bal"],
                        "fn_cost_vs_balanced": int(((pred_bal == 1) & (pred_x == 0) & (y_te == 1)).sum())})
            results.append(row)

        for cname, pack in comp_fitted.items():
            cols = c7_cols if cname == "C2_compact7" else tcr_cols
            score_te = apply_platt(pack["cal"], pack["model"].predict_proba(x_te[:, cols])[:, 1])
            record(cname, score_te, pack["t"])
        if t_c1 is not None:
            record("C1_phase4_delta", p4_te.astype(np.float32), t_c1)

        # frozen profiles (descriptive; biased in their favor here)
        union_spec = FeatureSpec(
            gene_names=union_genes, gene_indices=np.arange(len(union_genes), dtype=np.int32),
            gene_feature_names=[f"{g}_log1p_cp10k" for g in union_genes],
            engineered_feature_names=[], model_feature_names=union_genes,
            gene_to_col=pos_u, engineered_to_col={})
        obs_for_ann = {"source_gse_id": lab_te.source.to_numpy(),
                       "simple_annotation_plus6": lab_te.annotation.to_numpy(),
                       "TRA_cdr3": np.full(len(lab_te), "", dtype=object),
                       "TRB_cdr3": np.full(len(lab_te), "", dtype=object),
                       "TRG_cdr3": np.where(lab_te.any_gd.to_numpy(), "x", "").astype(object),
                       "TRD_cdr3": np.where(lab_te.any_gd.to_numpy(), "x", "").astype(object)}
        ann, _ = normalize_annotation(obs_for_ann, x_te, union_spec)
        for pid, (score, pred) in predict_profiles(profiles, x_te, union_spec, ann).items():
            row = metric_block(y_te[eval_mask], score[eval_mask], pred[eval_mask].astype(np.int8))
            row["abt_fpr"] = float(pred[abt_te].mean()) if abt_te.sum() else np.nan
            row["nk_fpr"] = float(pred[nk_te].mean()) if nk_te.sum() else np.nan
            row["nk_censored_callrate"] = float(pred[nk_cens_te].mean()) if nk_cens_te.sum() else np.nan
            row["nk_trdc_amb_callrate"] = float(pred[nk_trdc_te].mean()) if nk_trdc_te.sum() else np.nan
            row.update({"held_out": held_out, "candidate": f"frozen_{pid}", "threshold": "registered"})
            results.append(row)
            test_scores.append(pd.DataFrame({"held_out": held_out, "candidate": f"frozen_{pid}",
                                             "group": g_te[keep_rows],
                                             "y": y_te[keep_rows], "abt": abt_te[keep_rows],
                                             "nk": nk_te[keep_rows], "in_eval": eval_mask[keep_rows],
                                             "score": score[keep_rows], "pred": pred[keep_rows].astype(np.int8)}))

        fold_models[held_out] = {"chosen": {k: {"key": v["key"], "t_bal": v["t_bal"], "t_hp": v["t_hp"],
                                                "info": v["info"]} for k, v in chosen.items()}}
        with open(MODEL_DIR / f"kimi_fold_holdout_{held_out}{SUF}.pkl", "wb") as fh:
            pickle.dump({"held_out": held_out, "chosen": chosen, "comparators": comp_fitted,
                         "stage1": {"model": s1_full, "alpha": best_alpha, "thr": s1_thr},
                         "union_genes": union_genes, "note": "experimental; not registered; not promotable"}, fh)

        # GSE254249 stress with this fold's models (never trained on it)
        x_st, lab_st, _ = sources[STRESS_SOURCE]
        st_abt = (lab_st.kimi_class == "abT_censored_weak") | ((lab_st.kimi_class == "abT_primary"))
        st_nk = (lab_st.kimi_class == "strict_NK").to_numpy()
        s1_st = s1_full.predict_proba(x_st[:, s1_cols])[:, 1].astype(np.float32)
        x2_st = np.hstack([x_st[:, g197], derive_features(x_st[:, g197], names197), s1_st[:, None]]).astype(np.float32)
        for fam, pack in chosen.items():
            sc = apply_platt(pack["cal"], pack["model"].predict_proba(x2_st)[:, 1])
            pred = sc >= pack["t_bal"]
            results.append({"held_out": f"stress_GSE254249_via_{held_out}_fold", "candidate": f"{fam}_balanced",
                            "threshold": pack["t_bal"],
                            "abt_fpr": float(pred[st_abt.to_numpy()].mean()),
                            "nk_fpr": float(pred[st_nk].mean()) if st_nk.sum() else np.nan,
                            "note": "censored-abT FPR is an upper bound (assay cannot see gdT)"})
        log.info("fold %s done (%.1f min)", held_out, (time.time() - t0) / 60)

    res = pd.DataFrame(results)
    res.to_csv(TABLE_DIR / f"nested_metrics_by_heldout{SUF}.csv", index=False)
    ts = pd.concat(test_scores, ignore_index=True)
    ts.to_parquet(TABLE_DIR / f"heldout_cell_scores{SUF}.parquet", index=False)

    macro = res[~res.held_out.str.startswith("stress")].groupby("candidate").agg(
        macro_f1=("f1", "mean"), macro_recall=("recall", "mean"), macro_precision=("precision", "mean"),
        macro_specificity=("specificity", "mean"), macro_abt_fpr=("abt_fpr", "mean"),
        macro_nk_fpr=("nk_fpr", "mean"), macro_roc_auc=("roc_auc", "mean"), macro_pr_auc=("pr_auc", "mean"),
        macro_ece=("ece", "mean"), min_recall=("recall", "min"),
        max_abt_fpr=("abt_fpr", "max"), max_nk_fpr=("nk_fpr", "max")).reset_index()
    macro.to_csv(TABLE_DIR / f"nested_metrics_macro{SUF}.csv", index=False)

    # donor-cluster bootstrap CIs + paired differences on eval cells
    rng = np.random.default_rng(SEED + 6)
    boot_rows = []
    key_candidates = [c for c in ts.candidate.unique() if c.endswith("balanced") or c.startswith("frozen") or c.startswith("C")]
    for cand in key_candidates:
        sub = ts[(ts.candidate == cand) & ts.in_eval]
        per_ds_f1 = {}
        for ho, g in sub.groupby("held_out"):
            reps = []
            donors = g.group.unique()
            for _ in range(BOOTSTRAP_REPS):
                pick = rng.choice(donors, donors.size, replace=True)
                yy, pp = [], []
                for d in pick:
                    gg = g[g.group == d]
                    yy.append(gg.y.to_numpy()); pp.append(gg.pred.to_numpy())
                yy = np.concatenate(yy); pp = np.concatenate(pp)
                tp = ((pp == 1) & (yy == 1)).sum(); fp = ((pp == 1) & (yy == 0)).sum(); fn = ((pp == 0) & (yy == 1)).sum()
                pr = tp / max(tp + fp, 1); rc = tp / max(tp + fn, 1)
                reps.append(2 * pr * rc / max(pr + rc, 1e-12))
            per_ds_f1[ho] = np.asarray(reps)
        macro_reps = np.mean(np.vstack(list(per_ds_f1.values())), axis=0)
        boot_rows.append({"candidate": cand, "macro_f1_ci_lo": float(np.percentile(macro_reps, 2.5)),
                          "macro_f1_ci_hi": float(np.percentile(macro_reps, 97.5)),
                          "macro_f1_mean_boot": float(macro_reps.mean())})
        pd.DataFrame({"candidate": cand, "macro_f1_boot": macro_reps}).to_csv(
            TABLE_DIR / f"bootstrap_macro_f1_{cand}{SUF}.csv", index=False)
    pd.DataFrame(boot_rows).to_csv(TABLE_DIR / f"bootstrap_macro_f1_cis{SUF}.csv", index=False)

    (LOG_DIR / f"fold_model_choices{SUF}.json").write_text(json.dumps(fold_models, indent=2, default=str))
    log.info("all done in %.1f min", (time.time() - t0) / 60)


if __name__ == "__main__":
    main()
