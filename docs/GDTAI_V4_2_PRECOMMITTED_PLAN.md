# gdTAI V4.2 Precommitted Training and Validation Plan

## 1. Purpose

V4.2 tests whether the two-step T-cell-then-gdT architecture can outperform the
older gdTAI algorithms after repairing the NK-control definition that blocked
V4.1. This document freezes the scientific design before any V4.2 fitting.

V4.2 will either:

1. surpass the fair, fold-local older-algorithm comparators while satisfying all
   false-positive guardrails; or
2. conclude that this iteration did not improve the model under the available
   data and precommitted evaluation.

It cannot establish prospective clinical performance because no unused
TCR-confirmed positive cohort remains.

## 2. V4.1 failure being repaired

V4.1 used 336,780 transcriptome-annotated strict-NK controls. Only 21,054 had
agreement between the scVI superclass and an author/Phase-1 NK label; 315,726
had only the scVI NK label after conflict filtering. Several large
single-annotation sources are strongly T-like by CD3 and TCR-constant-gene
expression.

At the unchanged Stage-1 T-recall thresholds, all 15 saved V4.1 OOF candidate
evaluations pass the 50% source cap when evaluated against the 21,054
dual-annotation controls. This does not estimate V4.2 performance because V4.1
was trained with the weak labels. It establishes that NK-label provenance, not
an evaluated Stage-2 failure, caused Gate C to stop.

## 3. Immutable data roles

### 3.1 Primary TCR-defined truth

- `gdT_primary`: paired productive TRG/TRD, no productive TRA/TRB evidence, and
  no doublet flag.
- `abT_primary`: paired productive TRA/TRB, no productive TRG/TRD evidence, and
  no doublet flag.
- single-chain alpha-beta T cells remain reliability-weighted Stage-2 negatives.
- dual-receptor, silver, sorted-only, and ambiguous cells have zero fitting and
  threshold-selection weight.

The three primary T sources remain HRA005041, GSE144469, and
BALF_BLOOD_COPD. Each is held out in turn.

### 3.2 Primary NK truth

Primary Stage-1 NK controls must satisfy all of the following:

- scVI T/NK superclass equals NK;
- author/Phase-1 annotation explicitly supports NK;
- no productive TRA, TRB, TRG, or TRD evidence;
- no doublet flag;
- a valid donor, sample, or library grouping key.

This is the existing `nk_annotation_strength == 1.0` group. It contains 21,054
cells from GSE125527 and GSE228597.

### 3.3 Stress-only NK controls

Cells with `nk_annotation_strength == 0.5` receive zero fitting, calibration,
feature-selection, candidate-selection, and threshold-selection weight. Their
results remain mandatory by source as a label-uncertain stress test.

Expression-based T-lineage coherence is an audit variable only. It must never
create or remove NK truth, because doing so would make the evaluation circular.

### 3.4 Proposed new-data NK reference integration

Before V4.2 fitting, a separately approved modeling-reference integration will
be built without modifying the canonical atlas H5AD. The integration is for
conservative NK label expansion, not a classifier input.

Dataset roles are frozen before integration:

- development: GSE114724, GSE159251, GSE292700,
  GSE294273/GSE294274, and GSE296954;
- author-annotated external specificity: all GSE169246 cells, including 7,770
  author NK cells and 54,925 productive alpha-beta T cells;
- paired-alpha-beta external specificity: all 66,813 GSE315928 cells;
- mixed T/NK stress only: all GSE121636/GSE121637 cells.

The three held-out cohorts cannot contribute to HVG selection, scVI fitting,
clustering, cluster interpretation, feature selection, calibration, threshold
selection, or model-family selection.

The development integration will combine the current atlas with the five
development cohorts in a new sidecar object. Raw counts, source, donor, sample,
library, technology, and specimen provenance must pass a new pre-integration
QC gate. No source H5AD is modified.

Global scVI-neighbor clustering will run at three precommitted Leiden
resolutions and three seeds. The cytotoxic T/NK boundary is then subset and
reclustered at three resolutions and three seeds. A new-data cell may become a
cluster-consensus NK development negative only when all conditions hold:

- no productive TRA, TRB, TRG, or TRD evidence and no doublet flag;
- consensus NK membership in at least 80% of the global-plus-boundary runs;
- at least 95% NK purity among independently labeled NK-versus-productive-TCR
  anchor cells in its consensus cluster;
- no more than 2% productive-TCR contamination among all cells in that cluster;
- representation in at least three development sources, with no source
  contributing more than 70% of the cluster.

Primary paired-gdT inputs remain outside this sidecar integration so their
expression and labels cannot influence pseudo-label construction. Their
cluster-projection proximity is reported later as a diagnostic, not used as an
eligibility rule. During classifier fitting, all primary gdT cells remain
weight-1 Stage-1 T positives and Stage-2 gdT positives; cluster-derived NK
negatives cannot outweigh the primary NK pool.

No threshold on NKG7, GNLY, KLRD1, TYROBP, FCER1G, FCGR3A, CD3, TRDC, or TRDV
may define these labels. Clustering-derived NK cells are transcriptomic
pseudo-labels, receive reliability weight 0.5, and are source-balanced and
effective-weight capped during fitting. They cannot determine Stage-1
eligibility or final false-positive guardrails. Existing dual-annotation NK
controls remain weight-1 primary truth.

The deployed classifier continues to use log1p(CP10K) individual-gene values.
It never uses scVI latent coordinates, Leiden labels, or cluster membership.

## 4. Features

- The frozen 197-gene V4 universe remains the maximum candidate universe.
- Total deployed features must not exceed 300.
- No TRD score, TRAB score, TRD-minus-TRAB score, source ID, dataset ID, donor
  ID, library ID, or batch ID may be a feature.
- Stage 2 uses individual TCR alpha, beta, gamma, and delta genes plus the
  permitted transcriptomic context genes.
- The cross-fitted Stage-1 T-cell probability may be the only model-derived
  Stage-2 feature.
- No hard TRDC or TRDV expression requirement is allowed. Dropout must not
  automatically reject a gdT cell.
- NKG7, GNLY, KLRD1, TYROBP, FCER1G, FCGR3A, and other cytotoxic genes may
  contribute continuously, but none can act as a hard NK veto.
- B-cell and myeloid lineage panels are not added. The target application is
  already T/NK filtered, and they do not address the difficult gdT/CD8/NK
  boundary.
- Existing CD4-helper and Treg post-model exclusions remain fixed. Their recall
  cost is reported on every positive source and cannot be retuned after results
  are seen.

## 5. Two-stage architecture

### 5.1 Stage 1: high-recall T-cell gate

Stage 1 distinguishes productive-TCR T cells from primary dual-annotation NK
controls. Candidate families are deterministic GPU weighted ridge logistic and
GPU XGBoost, using the frozen V4.1 grids unless the implementation QC report
precommits a strictly smaller grid before fitting.

Training weights are source-balanced within class so GSE228597 cannot dominate
GSE125527. Cluster-consensus development negatives are separately
source-balanced, capped so their total effective weight cannot exceed the
weight-1 primary NK pool, and assigned reliability 0.5. All fold-local
preprocessing and calibration are learned from inner training groups only.

The Stage-1 threshold is the highest threshold satisfying:

- gdT-primary recall >= 0.99 in every adequately sized training source;
- productive alpha-beta T-cell recall >= 0.98 in every adequately sized
  training source;
- primary dual-annotation NK passage <= 0.50 in every NK source with at least
  100 controls.

Wilson 95% intervals are reported. Sources with fewer than 100 controls are
stress strata and do not determine eligibility.

### 5.2 Stage 2: gdT versus non-gdT T cells

Every cell passing Stage 1 enters Stage 2. Stage 2 predicts gdT-primary versus
alpha-beta T-cell negatives using individual genes and the cross-fitted Stage-1
probability. Stage-1 probability for any evaluation cell must come from a model
that did not train on that cell or its donor/sample/library/clonotype group.

Stage 1 is part of the final cascade: a Stage-1 rejection of a true gdT cell is
counted as a final false negative. No oracle T-cell subset is used for final
metrics.

## 6. Validation and comparison

### 6.1 Nested grouped validation

- Outer loop: leave each primary T source out in turn.
- Inner loop: preserve donor, then clonotype, then sample/library groups using
  the existing frozen hierarchy.
- NK robustness axis: leave each of the two primary NK sources out in turn for
  a separately reported Stage-1 source-transfer analysis.
- Extension cohorts remain untouched negative-control cohorts and cannot
  contribute to tuning, except for the five datasets preassigned to the
  cluster-consensus development role in Section 3.4. The three locked cohorts
  remain completely excluded from development.

### 6.2 Fair older-algorithm comparators

V2-like and V3-like feature/algorithm definitions are retrained from scratch on
the same outer and inner folds, with identical labels, preprocessing,
calibration, exclusions, and threshold-selection rules. These are the primary
comparators for a superiority claim.

The frozen released V2/V3 artifacts are also scored, but those results are
descriptive because their historical development reused parts of the current
benchmark.

### 6.3 Operating profiles

Two profiles are selected independently using inner OOF predictions:

- balanced: maximize F1, primary-gdT source-macro recall >= 0.80, paired-abT
  FPR <= 0.002, and primary NK FPR <= 0.01;
- high purity: maximize F0.5, primary-gdT source-macro recall >= 0.70,
  paired-abT FPR <= 0.001, and primary NK FPR <= 0.005.

All constraints are evaluated on the final two-stage cascade. Weak NK stress
FPR, extension-cohort paired-abT FPR, and extension-cohort NK FPR are reported
by source but do not select the threshold.

### 6.4 Promotion rule

V4.2 is promoted only if the balanced cascade:

1. improves pooled primary-gold F1 by at least 0.01 over the best fair
   older-algorithm comparator;
2. has a positive lower bound for the one-sided 95% paired hierarchical
   bootstrap interval of the F1 difference;
3. satisfies every balanced paired-abT and primary-NK guardrail;
4. has no primary positive source below 0.70 final recall;
5. does not worsen pooled weak-NK or extension-NK FPR by more than 0.002
   absolute versus the best fair comparator; and
6. completes all planned outer folds with converged, checksum-bound artifacts.

High-purity promotion is evaluated separately under its own guardrails. A
high-purity success cannot replace a failed balanced profile.

If no candidate passes, V3 Round 14 remains the balanced default and V2/V3
Round 12 remains the high-purity fallback. No threshold or label rule is changed
after seeing V4.2 results.

## 7. Required reports

The V4.2 evaluation report must include:

- complete label provenance and source counts;
- Stage-1 and final-cascade confusion matrices by source;
- precision, recall, specificity, F1, F0.5, balanced accuracy, MCC, ROC-AUC,
  and PR-AUC;
- paired-abT, primary-NK, weak-NK, and extension-control false-positive rates;
- GSE169246 author-NK and paired-abT false-positive rates, plus GSE315928
  paired-abT false-positive rates, reported as locked whole-dataset tests;
- bootstrap differences against every fair comparator;
- calibration and threshold-frontier plots;
- TRD-versus-TRAB score scatter plots for interpretation only, colored by V2,
  V3, and V4.2 predictions;
- failure-mode audits for gdT false negatives and NK/CD8 false positives;
- unchanged-input checksum and H5AD size/mtime checks.

## 8. Supervision gates

1. Step 0 label audit and this protocol are reviewed.
2. After explicit approval, implementation and synthetic/read-only tests form a
   separate QC gate.
3. Project-data fitting requires another explicit approval after implementation
   QC.
4. Promotion, release fitting, and atlas inference require review of the full
   nested evaluation.

No V4.2 model fitting is authorized by this document.
