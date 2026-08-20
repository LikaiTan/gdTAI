# gdTAI V4.6 (conditional-gamma two-stage classifier) — DEVELOPMENT ONLY

**Status: development artifact. NOT promoted, NOT the default model.**
The promoted default remains **gdTAI V3 Round 14**; the conservative fallback remains
**V2 high-purity**. V4.6 is non-promotable with current data: its single untouched
confirmatory cohort (BALF_BLOOD_COPD) was consumed by its own evaluation, and no
untouched gdT-positive confirmatory cohort remains.

## Contents

- `stage1_t_lineage_gate.ubj` — XGBoost T-lineage support gate (36 genes, threshold 0.0095)
- `stage2_gdt_classifier.ubj` — XGBoost γδT classifier (182 receptor-architecture features:
  145 individual TCR/CD3/CD4/FOXP3 genes + 37 engineered receptor aggregates; no
  standalone TRG genes, no NK/cytotoxic genes, no hard NK veto)
- `platt_calibration.json` — Platt calibration for Stage-2 scores
- `model_contract.json` — frozen feature lists, thresholds, guardrails, and provenance
  (SHA-256-pinned)

## Operating modes (frozen)

- `highest_f1`: threshold 0.894354 — development-validation P/R/F1 93.8% / 84.8% / 89.1%
- `high_purity`: threshold 0.971830 — development-validation P/R/F1 97.5% / 79.8% / 87.8%
- uncertainty policy: score < highest_f1 → negative; highest_f1 ≤ score < high_purity →
  uncertain (abstain band); score ≥ high_purity → high-confidence γδT

## Evidence summary

- Untouched BALF (28,328 cells, never used in any V4.6 role): highest-F1 P/R/F1
  98.99% / 92.14% / 95.44%, 8/27,476 negative FPs, 0/254 author-NK FPs.
- Whole-atlas application (5.93M cells, 40 sources, frozen model): 259,154 highest-F1
  calls with 0.09% strict likely-NK vs 24.1% for V2 high-F1.
- Independent NK evaluation (this repo, `gdT_prediction/gdtai_v4_6_nk_evaluation/`):
  on never-developed sources, call rates 0.004% (strict expression NK), 0.53%
  (author-NK), 1.1% (TRDC+TRDV−), 11.9% (shared-cytotoxic-ambiguous).
- Known limits: GDTlung recall 36.6% (limiting positive source), GSE228597 FPR 1.05%,
  lower ROC-AUC/PR-AUC than V3 on BALF (its F1 edge is operating-point-dependent).

## Input contract

Raw non-negative integer counts, normalized per cell as log1p(CP10K); T/NK-filtered
single-cell RNA-seq only (not a pan-cell-type classifier). Exact feature names and order
are pinned in `model_contract.json`; verify artifact SHA-256 before use.
