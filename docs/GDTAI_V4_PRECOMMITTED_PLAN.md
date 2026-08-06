# Precommitted gdTAI V4 Training and Validation Protocol

- **Protocol version:** 1.0
- **Frozen:** 2026-08-07
- **Status:** supervision required; model fitting has not started
- **Primary objective:** improve dataset-macro balanced F1 without relaxing the
  precommitted alpha-beta T-cell and NK false-positive guardrails

## 1. Purpose and governance

This protocol freezes the next gdTAI development cycle before any new fitting,
threshold search, or model comparison. It responds to the independent methodology
audit by replacing cell-level train/tune splits with nested cross-dataset evaluation,
using expression-independent TCR labels, and moving every feature decision,
calibration step, and threshold choice inside grouped resampling.

The current gdTAI v3 Round 14 model remains the promoted balanced default. Round 12
remains the validated high-purity fallback. The existing experimental V4 artifact is
not promoted and must be archived before a new V4 candidate is written. No result from
this protocol may be called independent external validation: all available positive
cohorts have already influenced prior development or review.

Training is blocked until this document, its rendered PDF, and the preflight audit are
reviewed. A second supervision gate is required after nested evaluation and before any
full-atlas inference or model promotion.

## 2. Fixed analysis population and data contract

### 2.1 Input population

- The model accepts T/NK-filtered single-cell RNA-seq data only. It is not a
  pan-cell-type classifier.
- Input values must be raw non-negative integer counts. Features are transformed as
  `log1p(counts / total_counts_per_cell * 10000)` without densifying the full matrix.
- The canonical development atlas is
  `high_speed_temp/Integrated_dataset/integrated.h5ad`, currently containing
  3,705,306 cells after the approved no-TCR-gene source carve-out.
- The GSE144469 expression matrix must come from the 107,068-cell raw-count object at
  `data/datasets/GSE144469/raw/legacy_source/fastq/tnk_subset.h5ad`. Its cells must be
  joined to canonical metadata by `SRR + barcode`; the join must be 100% unique and
  hashes of both source files must remain unchanged.
- BALF_BLOOD_COPD is a reused cross-study development benchmark, not an untouched
  test set. Its count layer is the only permitted expression source.
- H5AD inputs are read-only. Split assignments, labels, matrices, and predictions are
  written as separate tables or caches.

### 2.2 Required metadata

Every development row must provide a stable cell ID, `source_gse_id`, donor when
available, `sample_id`, `library_id`, raw-count provenance, productive chain fields for
TRA/TRB/TRG/TRD, canonical paired/any-chain flags, and doublet status when available.
Clonotype IDs are required where TCR sequence permits them.

### 2.3 Primary truth labels

Primary labels are assigned from productive TCR metadata only. Expression, annotations,
Phase-4 scores, and model predictions cannot change these labels.

| Class | Fixed rule | Primary weight |
| --- | --- | ---: |
| `gdT_primary` | Productive paired TRG and TRD, with no productive TRA or TRB | 1.0 |
| `abT_primary` | Productive paired TRA and TRB, with no productive TRG or TRD | 1.0 |
| `dual_or_ambiguous` | Evidence for both receptor classes, conflicting flags, or incomplete contradictory metadata | Excluded |
| `single_abT_weak` | Productive TRA or TRB, no productive TRG/TRD | 0.5, training only |
| `sorted_gdT_weak` | MalteGDT or eligible GDT_2020AUG_woCOV sorted cells without primary TCR truth | 0.25, training supplement only |
| `gdT_silver` | Project silver rule, including isolated productive TRD evidence | Sensitivity only |

GDTlung2023july_7p is sensitivity-only because of its suboptimal library quality.
Silver cells never enter feature selection, fitting, calibration, threshold selection,
or promotion metrics. Sorted weak positives may supplement an outer-training partition,
but never its held-out evaluation partition. TCR-confirmed CD4-like or Treg-like cells
remain positive truth and count as false negatives if a post-model exclusion rejects
them.

### 2.4 Dataset roles

- Primary nested evaluation sources: `HRA005041`, `GSE144469`, and
  `BALF_BLOOD_COPD`.
- Sorted sensitivity sources: `GDT_2020AUG_woCOV`, `MalteGDT`, and
  `GDTlung2023july_7p`.
- Negative stress source: `GSE254249`.
- Frozen extension negative-control cohorts: `GSE114724`,
  `GSE121636_GSE121637`, `GSE159251`, `GSE292700`,
  `GSE294273_GSE294274`, `GSE296954`, and `GSE315928`.
- `GSE169246` is a reduced-feature schema-sensitivity cohort unless its feature
  coverage is repaired and separately approved. It cannot contribute a primary
  guardrail in its current form.

The extension cohorts contain no gdT-gold or gdT-silver cells under project TCR rules.
They can measure paired/single alpha-beta T-cell and strict-NK false-positive behavior,
but not recall, F1, ROC-AUC, or PR-AUC.

## 3. Fixed feature policy

### 3.1 Exclusions

The classifier cannot use Phase-4 TRD, TRAB, or TRD-minus-TRAB scores; source IDs;
sample, donor, or library IDs; annotations; UMAP or latent coordinates; TCR metadata;
or doublet labels at inference. B-cell and myeloid genes, annotations, and program
scores are excluded from training and inference. B/myeloid annotations may be reported
only as out-of-scope contamination checks after predictions are frozen.

### 3.2 Gene universe

The feature universe is frozen at **197 genes**, below the 220-gene protocol cap and the
user's absolute 300-gene limit. All 197 were confirmed present in the canonical atlas,
BALF_BLOOD_COPD, and the GSE144469 raw-count object on 2026-08-07.

- **153 individual TCR genes:** the cross-input intersection listed in Appendix A.
- **14 T-lineage genes:** `CD3D`, `CD3E`, `CD3G`, `CD247`, `CD2`, `LCK`, `LAT`,
  `TRAT1`, `CD5`, `CD6`, `THEMIS`, `BCL11B`, `TCF7`, `LEF1`.
- **10 CD4/Treg-context genes:** `CD4`, `IL7R`, `CCR7`, `SELL`, `LTB`, `MAL`,
  `FOXP3`, `IL2RA`, `CTLA4`, `TIGIT`.
- **14 NK/cytotoxic-context genes:** `NKG7`, `GNLY`, `KLRD1`, `TYROBP`, `FCER1G`,
  `FCGR3A`, `KLRC1`, `KLRF1`, `PRF1`, `GZMB`, `CTSW`, `CCL5`, `EOMES`, `TBX21`.
- **6 T-cell-state genes:** `CD8A`, `CD8B`, `KLRB1`, `ZNF683`, `RUNX3`, `IKZF2`.

TRG genes are retained as individual weak evidence but receive no privileged rule or
hard threshold. NK/cytotoxic genes are contextual features, not NK-specific vetoes,
because bona fide gdT cells can express them strongly.

Inside each training fold, a gene is removed only if absent, constant, or detected in
fewer than 0.05% of that fold's cells. No outcome-driven DEG screen is permitted before
or outside the folds. Missing retained genes at inference are zero-filled only after
coverage checks; the model must abstain when coverage is inadequate.

### 3.3 Deterministic derived features

Derived features are restricted to detection counts for the TRA, TRB, TRG, and TRD
gene families; unweighted means of the fixed T-lineage, NK-context, CD4-helper, and
Treg panels; and the Stage-1 out-of-fold T-lineage probability. These are computed only
from the 197 genes. They are not the legacy Phase-4 module scores and use no matched
control-gene sampling. An individual-gene-only ablation is mandatory.

## 4. Model architecture and negative sampling

### 4.1 Stage 1: soft T-lineage gate

Stage 1 is an elastic-net logistic classifier that estimates T-lineage probability. It
uses the non-TCR genes plus `TRAC`, `TRBC1`, `TRBC2`, `TRDC`, `TRGC1`, and `TRGC2`.
TCR-defined alpha-beta and gamma-delta T cells are positive. NK negatives require an NK
annotation, no productive TRA/TRB/TRG/TRD evidence, and no doublet flag. When two
independent annotation fields are available they must agree; a single available NK
annotation is retained at weight 0.5.

The NK training sample is dataset-balanced and fixed as 60% representative NK, 20%
T-like hard NK, and 20% gdT-like/TRDC-positive-TRDV-negative NK. The hard strata alter
sampling only; they do not create NK truth labels. Stage 1 must achieve at least 99%
recall for TCR-defined gdT cells and 98% recall for productive alpha-beta T cells in
every primary outer evaluation source at its inner-fold-selected operating point.

Stage 1 is soft: all T/NK-input cells continue to Stage 2, and the out-of-fold Stage-1
probability is a Stage-2 feature. No single marker can hard-remove a cell before Stage 2.

Stage-1 hyperparameters are frozen to:

- `C`: 0.03, 0.1, 0.3, 1.0
- `l1_ratio`: 0.0, 0.25, 0.5, 0.75
- solver: `saga`; maximum iterations: 5,000; tolerance: 1e-4

### 4.2 Stage 2: gdT classifier

Two candidate families are compared under identical folds and weights:

1. Elastic-net logistic regression using the 197 genes, allowed derived features, and
   Stage-1 out-of-fold probability. Grid: `C` = 0.01, 0.03, 0.1, 0.3, 1.0;
   `l1_ratio` = 0.0, 0.25, 0.5, 0.75, 1.0; `saga`, 5,000 iterations, tolerance 1e-4.
2. Histogram gradient boosting with `learning_rate` = 0.03 or 0.07,
   `max_leaf_nodes` = 7 or 15, `min_samples_leaf` = 100 or 500,
   `l2_regularization` = 1 or 10, `max_iter` = 250, and early stopping disabled so
   validation folds cannot leak into fitting.

No additional candidate family or grid point may be introduced after results are seen.
Any later change starts a separately named, separately precommitted experiment.

### 4.3 Weighting

Training weights are the product of truth reliability, inverse dataset frequency within
class, and inverse class frequency within dataset. Weights are normalized to mean 1 in
each training fold. No cell is duplicated. Donor/library groups cannot be split between
train and validation. Sorted weak positives and single-chain weak negatives use the
fixed reliability weights in Section 2.3.

## 5. Nested evaluation and operating points

### 5.1 Outer and inner splits

- Outer evaluation is leave-one-dataset-out across HRA005041, GSE144469, and
  BALF_BLOOD_COPD. The held-out source is absent from both model stages, feature
  filtering, calibration, and threshold selection.
- Each outer-training partition uses three-fold `StratifiedGroupKFold` for tuning.
- The group key is donor when available, otherwise sample, otherwise library. If donor
  is unavailable and a clonotype spans multiple fallback groups, the entire clonotype
  is assigned to one fold.
- Split manifests are generated once from a fixed seed (`20260807`), checksum-pinned,
  and reused by every candidate and comparator.
- Every feature filter, scaler, derived feature, Stage-1 prediction, hyperparameter
  decision, sigmoid calibrator, and operating threshold is fitted using inner-training
  data only. Isotonic calibration is sensitivity-only and cannot select a model.

### 5.2 Threshold selection

Threshold candidates are the sorted unique inner out-of-fold calibrated probabilities
plus 0 and 1. Ties are resolved by lower strict-NK FPR, then lower paired-abT FPR, then
fewer retained genes, then lexical candidate ID.

- **Balanced mode:** maximize dataset-macro F1 subject to macro recall at least 0.80,
  paired-abT FPR at most 0.20%, and strict-NK FPR at most 1.00%.
- **High-purity mode:** maximize dataset-macro F0.5 subject to macro recall at least
  0.70, paired-abT FPR at most 0.10%, and strict-NK FPR at most 0.50%.

If no threshold satisfies a mode's guardrails, that mode fails for the candidate; the
constraints cannot be relaxed. The final deployment thresholds are selected from
cross-dataset out-of-fold predictions only after the candidate family is frozen. The
resulting deployment thresholds are not reported as independently validated estimates.

### 5.3 CD4/Treg exclusions and VDJ-aware rescue

Both RNA operating modes reject coherent CD4-helper-like and Treg-like predictions after
probability thresholding:

- CD4-helper exclusion: `CD4 > 0`, at least two of `IL7R`, `CCR7`, `SELL`, `LTB`,
  `MAL` are detected, and their six-gene mean is at least 0.5.
- Treg exclusion: `FOXP3 > 0`, at least one of `IL2RA`, `CTLA4`, `TIGIT` is detected,
  and their four-gene mean is at least 0.5.

These rules operate on log1p(CP10K) values and never alter truth labels. A separate,
optional VDJ-aware output may rescue rejected cells with productive paired TRG/TRD and
no productive TRA/TRB. Rescued calls are reported separately and never enter TCR-defined
validation metrics, avoiding use of the truth fields to inflate model performance.

## 6. Comparators, metrics, and statistics

### 6.1 Fair retrained comparators

Every comparator uses the same outer and inner partitions:

1. Legacy TRD-minus-TRAB score with a fold-selected threshold.
2. Compact logistic model using `TRDC`, `TRDV1`, `TRDV2`, `TRDV3`, `TRAC`, `TRBC1`,
   and `TRBC2`.
3. V2-like single-stage elastic-net logistic regression using the fixed individual TCR
   genes only.

Frozen gdTAI V2 high-purity, V3 Round 12, and V3 Round 14 are descriptive reference
profiles only because they cannot be fairly retrained from one coherent recorded build.

### 6.2 Primary reporting

For each outer dataset and the unweighted dataset macro-average, report ROC-AUC, PR-AUC,
precision, recall, specificity, F1, F0.5, balanced accuracy, MCC, confusion counts,
paired-abT FPR, strict-NK FPR, and calibration intercept/slope, Brier score, and expected
calibration error. Also report silver recall, sorted-source sensitivity, CD4/Treg
exclusion losses, TRDC-positive/TRDV-negative behavior, feature coverage, abstention,
and PPV/FDR at gdT prevalences of 0.1%, 0.5%, 1%, 2%, and 5%.

Use 2,000 paired hierarchical bootstrap replicates, resampling datasets and then donors,
or libraries when donor is absent. Report percentile 95% confidence intervals for each
metric and paired candidate-minus-comparator differences. Dataset-level and worst-dataset
results are primary; pooled cell-level results are secondary.

### 6.3 Promotion rule

V4 is promoted only if all conditions hold for the balanced mode:

1. Dataset-macro F1 exceeds the strongest fairly retrained comparator by at least 0.01.
2. The 95% paired hierarchical-bootstrap lower confidence bound for that F1 difference
   is greater than 0.
3. Recall is at least 0.70 in every primary outer dataset.
4. All balanced paired-abT and strict-NK guardrails pass in every primary outer dataset
   and in the seven complete-feature extension negative-control cohorts.
5. Feature coefficients or permutation ranks are directionally stable in at least two
   of three outer folds, and no single dataset accounts for more than 50% of aggregate
   feature importance.
6. The release bundle passes semantic and checksum consistency tests.

Failure of any condition retains V3 Round 14 as the promoted default and publishes V4 as
an experimental negative result. High-purity mode may be distributed only when its own
guardrails pass; it cannot rescue a failed balanced-mode promotion.

## 7. Full-atlas inference and release contract

Whole-atlas inference is blocked until the nested report is reviewed and promotion is
explicitly authorized. After authorization, score the current 3,705,306-cell atlas once
with the frozen final fit. Each output row must contain:

- calibrated Stage-1 T-lineage probability
- calibrated Stage-2 gdT probability
- balanced and high-purity RNA calls
- CD4-helper and Treg exclusion flags
- feature coverage and schema-abstention flag
- optional VDJ-aware rescued call
- deterministic decision reason
- model ID, artifact SHA-256, Git commit, and split-manifest ID

The model abstains when any critical gene (`CD3D`, `CD3E`, `TRAC`, `TRBC1`, `TRBC2`,
`TRDC`, `TRDV1`, `TRDV2`, `TRDV3`) is absent or when retained model-gene coverage is
below 90%. An immutable release bundle must include the model, feature order, coefficients
or importances, normalizer, calibrators, thresholds, input hashes, split IDs, environment
lock, command, random seeds, metrics, model card, inference example, and a manifest whose
semantics agree with the serialized artifact and registry.

## 8. Step gates and required artifacts

### Step 0: protocol freeze (this document)

- Commit this Markdown protocol, reader HTML, and PDF before any fitting.
- Email the PDF and commit SHA to Likai with supervision status and exact paths.
- Stop. No training or threshold search is authorized by completion of Step 0.

### Step 1: preflight and split freeze

- Archive the existing experimental V4 artifact without overwriting it.
- Verify raw-count state, hashes, required metadata, label conflicts, class counts,
  feature coverage, and the GSE144469 `SRR + barcode` join.
- Emit checksum-pinned cell-label and grouped-split manifests plus a preflight report.
- Email Likai and stop for approval.

### Step 2: nested evaluation

- Run both stages, ablations, and fair comparators under the frozen folds.
- Emit per-fold predictions, metrics, calibration, feature stability, bootstrap results,
  extension negative controls, HTML, and a print-safe PDF.
- Email Likai and stop for model-selection review.

### Step 3: final fit and decision

- Only after approval, fit the selected family on all development data, select final
  cross-validated operating thresholds, run stress cohorts, and apply the promotion rule.
- Publish either an immutable promoted V4 release or an experimental failure record.
- Whole-atlas inference requires separate explicit approval.

## 9. Mandatory automated checks

Before any candidate can run, tests must prove:

1. All expression inputs are raw counts and remain unchanged by size, mtime, and SHA-256.
2. Primary positive and negative labels are nonempty, mutually exclusive, and independent
   of expression and annotations; silver cells are absent from fitting and selection.
3. GSE144469 maps 107,068 of 107,068 cells uniquely by `SRR + barcode`.
4. The gene universe contains exactly 197 genes, no Phase-4 score, no metadata feature,
   and no B/myeloid feature; every fold retains at most 220 genes.
5. No donor, fallback sample/library, or protected clonotype crosses a fold boundary.
6. Feature filtering, Stage-1 probabilities, scaling, calibration, model selection, and
   thresholds are fold-local.
7. CD4/Treg exclusions and VDJ rescue are deterministic and rescue never changes the RNA
   metric frame.
8. Missing critical genes or less than 90% coverage causes abstention.
9. Two runs with the same seed produce identical split manifests, predictions, and model
   hashes, or a documented deterministic numerical tolerance where the estimator cannot
   be bitwise stable.
10. Serialized model, manifest, registry, report, and example inference outputs agree on
    feature order, thresholds, version, normalization, and checksums.

## Appendix A. Frozen individual TCR genes (153)

`TRAC`, `TRAJ10`, `TRAJ12`, `TRAJ13`, `TRAJ16`, `TRAJ18`, `TRAJ22`, `TRAJ23`,
`TRAJ27`, `TRAJ28`, `TRAJ32`, `TRAJ34`, `TRAJ37`, `TRAJ38`, `TRAJ39`, `TRAJ42`,
`TRAJ43`, `TRAJ44`, `TRAJ45`, `TRAJ46`, `TRAJ47`, `TRAJ49`, `TRAJ5`, `TRAJ6`,
`TRAJ9`, `TRAV1-1`, `TRAV1-2`, `TRAV10`, `TRAV12-1`, `TRAV12-2`, `TRAV12-3`,
`TRAV13-1`, `TRAV13-2`, `TRAV14DV4`, `TRAV15`, `TRAV16`, `TRAV17`, `TRAV18`,
`TRAV19`, `TRAV2`, `TRAV20`, `TRAV21`, `TRAV22`, `TRAV23DV6`, `TRAV24`, `TRAV25`,
`TRAV26-1`, `TRAV26-2`, `TRAV27`, `TRAV29DV5`, `TRAV3`, `TRAV30`, `TRAV34`,
`TRAV35`, `TRAV36DV7`, `TRAV38-1`, `TRAV38-2DV8`, `TRAV39`, `TRAV4`, `TRAV40`,
`TRAV41`, `TRAV5`, `TRAV6`, `TRAV8-1`, `TRAV8-2`, `TRAV8-3`, `TRAV8-4`,
`TRAV8-5`, `TRAV8-6`, `TRAV9-2`, `TRBC1`, `TRBC2`, `TRBJ1-1`, `TRBJ1-2`,
`TRBJ1-4`, `TRBJ1-5`, `TRBJ1-6`, `TRBJ2-1`, `TRBJ2-2`, `TRBJ2-2P`, `TRBJ2-3`,
`TRBJ2-4`, `TRBJ2-5`, `TRBJ2-6`, `TRBJ2-7`, `TRBV1`, `TRBV10-1`, `TRBV10-2`,
`TRBV10-3`, `TRBV11-1`, `TRBV11-2`, `TRBV11-3`, `TRBV12-2`, `TRBV12-3`,
`TRBV12-4`, `TRBV12-5`, `TRBV13`, `TRBV14`, `TRBV15`, `TRBV16`, `TRBV18`,
`TRBV19`, `TRBV2`, `TRBV20-1`, `TRBV21-1`, `TRBV23-1`, `TRBV24-1`, `TRBV25-1`,
`TRBV27`, `TRBV28`, `TRBV29-1`, `TRBV3-1`, `TRBV30`, `TRBV4-1`, `TRBV4-2`,
`TRBV5-1`, `TRBV5-3`, `TRBV5-4`, `TRBV5-5`, `TRBV5-6`, `TRBV6-1`, `TRBV6-2`,
`TRBV6-4`, `TRBV6-5`, `TRBV6-6`, `TRBV7-2`, `TRBV7-3`, `TRBV7-4`, `TRBV7-5`,
`TRBV7-6`, `TRBV7-7`, `TRBV7-9`, `TRBV9`, `TRDC`, `TRDJ1`, `TRDJ2`, `TRDJ3`,
`TRDV1`, `TRDV2`, `TRDV3`, `TRGC1`, `TRGC2`, `TRGV1`, `TRGV10`, `TRGV11`,
`TRGV2`, `TRGV3`, `TRGV4`, `TRGV5`, `TRGV5P`, `TRGV7`, `TRGV8`, `TRGV9`.

## Appendix B. Precommit record

- Canonical environment: `rapids_sc_py310`
- Random seed: `20260807`
- Resource design limit: 80 CPU cores and 400 GB RAM; the active runbook permits up to
  800 GB, but the workflow cannot require it.
- Current promoted comparator: gdTAI v3 Round 14, threshold 0.936
- High-purity fallback: gdTAI v3 Round 12, threshold 0.5
- Existing experimental V4: not promoted; archive before Step 1
- No training, calibration, threshold search, full-atlas inference, or model promotion was
  performed as part of Step 0.
