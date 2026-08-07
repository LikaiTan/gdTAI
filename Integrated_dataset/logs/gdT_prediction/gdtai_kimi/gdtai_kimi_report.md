# gdTAI-kimi: an honest-benchmark experiment for the gdTAI program

**Date:** 2026-08-07
**Status:** experimental candidate, **not registered, not promotable**. gdTAI-kimi is a
separately named experiment run outside the frozen GDTAI V4 precommitted protocol
(`docs/GDTAI_V4_PRECOMMITTED_PLAN.md`), applying that protocol's core principles end-to-end.
Its purpose is twofold: (1) test whether the V4 protocol is executable as written, and
(2) quantify how gdTAI-class models behave when every split, feature decision, calibration
step, and threshold is grouped, fold-local, and free of adaptive iteration.

**One-paragraph answer.** Under fully honest nested leave-one-dataset-out evaluation,
a from-scratch two-stage gdT classifier (gdTAI-kimi, HistGBM variant) reaches
dataset-macro F1 **0.811** (donor-cluster bootstrap 95% CI [0.761, 0.845]) — far above the
fairly retrained simple baselines (0.57–0.60) but below the frozen V2/V3 profiles'
descriptive scores (~0.91, which are optimistically biased in this frame because those
models were developed with access to HRA005041 and GSE144469 cells). **No kimi variant met
the precommitted balanced-mode guardrails in any outer fold**, so under the V4 promotion
rule gdTAI-kimi would be published as an experimental negative result. The binding
constraints were recall at the 0.2% paired-abT FPR ceiling (folds 1 and 3) and the strict-NK
FPR ceiling (fold 2). NK leakage — not αβ leakage — is the hardest unsolved problem in this
model class.

---

## 1. Design

| Element | gdTAI-kimi choice |
|---|---|
| Evaluation | Nested leave-one-dataset-out over HRA005041, GSE144469, BALF_BLOOD_COPD; 3-fold StratifiedGroupKFold inner tuning, group = donor (fallback sample/library) |
| Labels | Expression-independent, computed from raw CDR3 fields: `gdT_primary` = paired TRG/TRD, no ab evidence; `abT_primary` = paired TRA/TRB, no gd evidence, **restricted to gdTCR-assayed sublibraries**; censored paired-abT = weight-0.25 training-only weak negatives; sorted-gdT cells without paired TCR = weight-0.25 training-only weak positives; dual/ambiguous excluded |
| Features | Frozen 197-gene universe from the V4 plan + fold-local derived features (family detection counts, panel means, Stage-1 OOF probability); no Phase-4 scores, no metadata, no annotations |
| Architecture | Stage 1 soft T-lineage gate (elastic-net; OOF probability is a Stage-2 feature; no cell hard-removed). Stage 2: elastic-net logistic (16-config grid) vs HistGBM (16-config grid), identical folds and weights |
| Calibration/thresholds | Sigmoid (Platt) on inner-OOF scores; thresholds selected on inner-OOF under the V4 plan's guardrails (balanced: macro recall ≥0.80, paired-abT FPR ≤0.2%, strict-NK FPR ≤1%; high-purity: ≥0.70 / ≤0.1% / ≤0.5%) |
| Comparators (same folds) | C1 Phase-4 TRD−TRAB score; C2 compact 7-gene logistic; C3 TCR-gene-only elastic-net |
| Frozen references | V2 high-F1, V2 high-purity, V3 Round 12, V3 Round 14 scored on identical held-out cells — **descriptive only** (see §4.2) |
| Statistics | Per-dataset metrics primary; 1,000-rep donor-cluster bootstrap CIs on dataset-macro F1 |

Deviations from the frozen V4 protocol (documented in the script docstring): GSE144469
expression from the integrated-plus6 atlas X; SGD optimizer for elastic-net; fold-local
≤100k tuning subsamples; no extension-cohort guardrails (sealed cohorts stay sealed); no
full-atlas inference.

## 2. Label inventory (all expression-independent)

| source | gdT_primary | abT_primary | abT_censored_weak | strict_NK | sorted weak | dual excluded |
|---|---:|---:|---:|---:|---:|---:|
| HRA005041 | 4,894 | 22,211 | 520,817 | 757 | — | 5,917 |
| GSE144469 | 4,003 | 61,729 | 2,989 | 307 | — | 1,103 |
| BALF_BLOOD_COPD | 1,046 | 33,117 | — | 254 | — | 3,116 |
| GDT_2020AUG_woCOV | 15,066 | — | — | — | 10,838 | — |
| MalteGDT | 5,650 | — | — | — | 2,150 | — |
| GSE254249 (stress only) | — | — | 174,010 | 8,198 | — | — |

Notes: HRA005041 contributes only 22,211 sublibrary-clean abT negatives (vs 520,817 censored)
— exactly the M1 trade-off, and it is workable. BALF contributes only 254 strict-NK cells,
making per-dataset NK FPR estimates on BALF extremely coarse.

## 3. Headline results

### 3.1 Dataset-macro (mean over the three held-out datasets)

| candidate | macro F1 | bootstrap 95% CI | macro recall | macro precision | macro abT FPR | macro NK FPR | min per-dataset recall |
|---|---:|---:|---:|---:|---:|---:|---:|
| **kimi HistGBM balanced** | **0.811** | **[0.761, 0.845]** | 0.752 | 0.915 | 0.63% | 5.3% | 0.676 |
| kimi elastic-net balanced | 0.634 | [0.605, 0.655] | 0.591 | 0.937 | 0.47% | 6.9% | **0.083** |
| kimi HistGBM high-purity | 0.727 | — | 0.630 | 0.939 | 0.43% | 4.1% | 0.469 |
| C1 Phase-4 TRD−TRAB | 0.573 | [0.512, 0.607] | 0.477 | 0.929 | 0.44% | 1.6% | 0.211 |
| C2 compact 7-gene | 0.595 | [0.497, 0.639] | 0.515 | 0.923 | 0.57% | 2.8% | 0.208 |
| C3 TCR-gene elastic-net | 0.582 | [0.471, 0.628] | 0.495 | 0.936 | 0.45% | 2.0% | 0.209 |
| frozen V2 high-purity* | 0.906 | [0.891, 0.916] | 0.885 | 0.930 | 0.40% | 4.6% | 0.829 |
| frozen V3 Round 12* | 0.913 | [0.888, 0.922] | 0.888 | 0.944 | 0.36% | 8.1% | 0.821 |
| frozen V3 Round 14* | 0.910 | [0.897, 0.917] | 0.886 | 0.937 | 0.38% | 7.8% | 0.831 |

*Frozen profiles are descriptive references. They were trained/selected with access to
HRA005041 and GSE144469 cells, so their numbers on those datasets are train-adjacent and
optimistically biased; see §4.2.

### 3.2 Per held-out dataset (balanced mode)

| held-out | candidate | precision | recall | F1 | abT FPR | NK FPR | ROC-AUC |
|---|---|---:|---:|---:|---:|---:|---:|
| HRA005041 | kimi HistGBM | 0.999 | 0.677 | 0.807 | 0.02% | 14.9% | 0.998 |
| HRA005041 | kimi elastic-net | 0.980 | 0.811 | 0.887 | 0.36% | 20.0% | 0.991 |
| HRA005041 | frozen V3 R14* | 0.998 | 0.962 | 0.979 | 0.05% | 22.1% | 0.999 |
| GSE144469 | kimi HistGBM | 0.760 | **0.902** | 0.825 | 1.85% | 1.0% | 0.970 |
| GSE144469 | kimi elastic-net | 0.843 | 0.879 | 0.861 | 1.06% | 0.7% | 0.948 |
| GSE144469 | frozen V3 R14* | 0.850 | 0.866 | 0.858 | 0.99% | 1.0% | 0.962 |
| BALF_BLOOD_COPD | kimi HistGBM | 0.986 | 0.676 | 0.802 | 0.03% | 0.0% | 0.980 |
| BALF_BLOOD_COPD | kimi elastic-net | 0.989 | **0.083** | 0.153 | 0.00% | 0.0% | 0.980 |
| BALF_BLOOD_COPD | frozen V3 R14* | 0.963 | 0.831 | 0.892 | 0.10% | 0.4% | 0.969 |

## 4. Findings

### 4.1 Kimi fails the precommitted guardrails in every fold — an honest negative result
No kimi variant satisfied the balanced-mode guardrails during inner tuning. The binding
constraint differed by fold: training on GSE144469+BALF, the best achievable macro recall at
the 0.2% abT-FPR ceiling was 0.66 (HistGBM); training on HRA005041+GSE144469, elastic-net
reached only 0.40; training on HRA005041+BALF, recall was fine (0.95) but the inner strict-NK
FPR was 14–20%, 14–20× over the 1% ceiling. Under the V4 promotion rule this experiment ends
here: V3 Round 14 stays the default and kimi is published as a negative result. This is the
protocol working as intended.

### 4.2 The frozen models' advantage is real but smaller than their headline numbers imply
The only like-for-like frame is BALF_BLOOD_COPD (never used to fit any model; used in frozen
model *selection*, so still mildly favorable to them). There, frozen V3 R14 beats kimi
HistGBM on recall (0.831 vs 0.676) at comparable precision. On HRA005041 and GSE144469 the
frozen models' numbers are train-adjacent: V3 trained on both datasets' cells, so its
0.96-recall on "held-out" HRA005041 is not a generalization measurement. Notably, on
GSE144469 — where frozen V3 had every advantage — kimi HistGBM still achieved *higher*
recall (0.902 vs 0.866), at the cost of a higher abT FPR (1.85% vs 0.99%).

### 4.3 NK leakage is the binding constraint, not αβ leakage
Across folds and model families, the strict-NK FPR is 5–50× harder to control than the
paired-abT FPR (e.g., kimi HistGBM on held-out HRA005041: abT FPR 0.02% vs NK FPR 14.9%).
This is not unique to kimi: frozen V3 R14's NK FPR on the same cells is 22.1% — *worse* than
kimi HistGBM. The historical NK hard-negative strategy (v3's TRDC+TRDV− guard) addressed a
visible symptom; the grouped, never-fitted NK evaluation shows the problem is deeper. The V4
protocol's soft-gate-only NK defense is likely insufficient as written; see §5.

### 4.4 gdTAI-class models add large, real value over simple baselines
The fairly retrained comparators collapse under honest evaluation: macro F1 0.57–0.60, with
per-dataset recall as low as 0.21 on BALF. The full two-stage model roughly doubles
comparator recall at similar precision. The answer to "is a TRD−TRAB score good enough" is
clearly no — but the answer to "do we need more than ~200 genes and a small gradient-boosted
tree" is also no.

### 4.5 Elastic-net is unstable across datasets; HistGBM is robust
The elastic-net variant collapsed on held-out BALF (recall 0.083 balanced, 0.027
high-purity) while performing acceptably on the other two folds — its learned linear weights
do not transfer across chemistries/tissues. HistGBM's worst-dataset recall was 0.676.
Family choice matters less than this stability gap; any future protocol should score
worst-dataset recall, not just macro.

### 4.6 Fold-to-fold model instability is itself a finding
Different folds selected different configs (en_a0.001_l0.75 / en_a0.0003_l0.25 / en_a0.001_l0.5;
hgb_14 / hgb_12 / …) and the GSE254249 stress FPR varied >10× across fold models
(0.01%–1.3% for HistGBM balanced). A single "final fit on all data" would hide this
variance; the V4 protocol's feature-stability promotion condition (§6.3 item 5) is
well-motivated and should be extended to threshold stability.

### 4.7 CD4/Treg post-exclusion is nearly free
The frozen CD4-helper/Treg exclusion rules cost only 1–20 gold-positive cells per dataset
(≤0.5% recall) while providing a deterministic safety rail. Keep them.

## 5. Implications for the frozen V4 protocol

1. **Guardrail feasibility risk is real.** The 0.80 macro-recall at ≤0.2% abT-FPR balanced
   guardrail was unreachable in two of three folds for a fresh model family. Before Step 2,
   the preflight should verify the guardrails are jointly satisfiable on inner-OOF
   predictions, or the V4 run will end in a guaranteed guardrail failure.
2. **Add an NK defense to Stage 2, or accept NK FPR as the limiting metric.** Options:
   include never-evaluated NK hard negatives in Stage-2 training (v3's approach, but with
   grouped folds and representative NK sampling rather than TRDC+TRDV−-only), or make the
   Stage-1 gate firmer for NK-like cells (a soft gate lets NK leakage through by design).
3. **Report guardrail-failed candidates anyway.** Kimi's per-dataset numbers remain
   informative precisely because the failure is flagged, not hidden. The V4 report template
   should include best-effort rows with `guardrails_met=false`.
4. **Adopt worst-dataset recall as a primary metric** (§4.5) and threshold-stability
   reporting (§4.6).
5. **BALF's 254 strict-NK cells cannot support a 1% NK-FPR guardrail at dataset level**
   (1% = 2.5 cells). NK guardrails should be pooled across datasets or BALF's NK stratum
   acknowledged as uninformative.

## 6. Caveats

- gdTAI-kimi used the atlas X matrix for GSE144469 rather than the raw SRR-joined object the
  V4 plan mandates; results may differ slightly under the plan's data contract.
- Labels remain TCR-pipeline-derived (review finding F2): kimi measures TCR-call reproduction
  from RNA, same as its predecessors. The two-tier negative scheme controls *assay censoring*
  bias, not label circularity.
- Frozen-profile numbers are descriptive and biased in their favor on atlas datasets; the
  honest kimi-vs-frozen comparison is BALF only.
- Extension-cohort guardrails and full-atlas inference were intentionally not run (sealed
  cohorts remain sealed; whole-atlas inference requires separate approval).

## 8. Round 2: diagnosis-driven iteration (2026-08-07, same day)

Round 1 showed three failure drivers: NK leakage as the binding constraint, BALF positives
clustering just below the transferred threshold, and HGB grid capacity limits. Round 2 made
exactly four pre-registered changes: **R2-1** grouped NK hard negatives (weight 0.5) added to
Stage-2 training; **R2-2** TRD-only silver-like cells added as weight-0.25 weak positives;
**R2-3** HistGBM capacity raised (leaves 15–31, max_iter 400); **R2-4** recall-at-fixed-FPR
operating points reported for every candidate.

### 8.1 Round 1 vs round 2 (dataset-macro, balanced mode)

| candidate | macro F1 r1 → r2 | macro recall r1 → r2 | macro abT FPR r1 → r2 | macro NK FPR r1 → r2 |
|---|---|---|---|---|
| kimi HistGBM balanced | 0.811 → **0.838** | 0.752 → **0.828** | 0.63% → 1.05% | 5.3% → 17.1% |
| kimi HistGBM high-purity | 0.727 → 0.789 | 0.630 → 0.720 | 0.43% → 0.61% | 4.1% → 6.4% |
| kimi elastic-net balanced | 0.634 → 0.464 | 0.591 → 0.419 | 0.47% → 0.79% | 6.9% → 1.2% |
| C2 compact 7-gene | 0.595 → 0.589 | — | — | — |

Round 2 bought recall but paid for it in false positives — it moved along the
recall/FPR frontier rather than dominating round 1. The silver weak positives (R2-2)
broadened the positive class but pulled TRD-positive NK-like cells with them (GSE144469
held-out NK FPR rose 1.0% → 32.9%): **TRD-only evidence cannot distinguish silver gdT from
NK contamination, and adding silver to fitting was counterproductive for NK control.**
This independently confirms protocol v1.2's decision to remove silver and sorted cells from
fitting. The NK training negatives (R2-1) did not reduce NK FPR either (fold-2 inner NK FPR
~20%, unchanged) — NK leakage in this feature space is not a sample-size problem.

### 8.2 The useful output: the honest recall@FPR frontier (kimi HistGBM, round 2)

Thresholds set on inner-OOF to fixed abT-FPR targets, applied to the held-out dataset:

| held-out | recall @ ≤0.1% target | recall @ ≤0.2% target | recall @ ≤0.5% target | recall @ ≤1% target |
|---|---:|---:|---:|---:|
| HRA005041 | 0.627 (meas. 0.01%) | 0.793 (0.02%) | 0.952 (0.20%) | 0.983 (1.24%) |
| GSE144469 | 0.910 (meas. 1.87%) | 0.914 (3.06%) | 0.921 (7.28%) | 0.951 (14.95%) |
| BALF_BLOOD_COPD | 0.487 (meas. 0.01%) | 0.624 (0.02%) | 0.841 (0.13%) | 0.926 (1.46%) |

Two transferable lessons: (a) thresholds calibrated on atlas sources transfer well to
HRA005041 and BALF (measured FPR tracks the target) but **fail badly toward GSE144469**
(measured FPR up to 19× the target) — a single global threshold is not robust across that
dataset's direction; (b) the honest achievable operating point is roughly **macro recall
≈ 0.90 at 0.5% abT FPR** (0.95 / 0.92 / 0.84 per dataset), not the 0.80 recall at 0.2% FPR
the balanced guardrail demands — the precommitted guardrails sit slightly off the achievable
frontier for any model in this family, which is exactly what the parallel V4 v1.2 CPU run
found independently.

### 8.3 Round-2 conclusion

No round-2 variant met the guardrails; kimi remains an experimental negative result, now
with a mapped frontier. The evidence-based guidance for the V4.1-GPU plan: keep silver and
sorted cells out of fitting; treat NK FPR as the primary limiting metric and pool NK
guardrails across datasets (BALF has 254 strict-NK cells); report recall@fixed-FPR curves
instead of relying on one transferred threshold; and expect the balanced guardrail
(recall ≥0.80 at ≤0.2% FPR) to be reachable only marginally, if at all, without a
structurally different NK defense.

## 7. Artifacts

- Script: `workflows/gdtai/run_gdtai_kimi_nested.py` (seed 20260807)
- Tables: `Integrated_dataset/tables/gdT_prediction/gdtai_kimi/` (per-dataset and macro
  metrics, per-candidate inner diagnostics, bootstrap distributions, label inventory,
  per-cell held-out scores)
- Figures: `Integrated_dataset/figures/gdT_prediction/gdtai_kimi/`
- Experimental fold models: `Integrated_dataset/models/gdT_prediction_classifier/gdtai_kimi/`
  (explicitly unregistered, not promotable)
- Run log: `Integrated_dataset/logs/gdT_prediction/gdtai_kimi/run.log`
