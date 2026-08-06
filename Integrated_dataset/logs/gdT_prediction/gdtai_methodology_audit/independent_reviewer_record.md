# Independent gdTAI Reviewer Record

Date: 2026-08-06

Reviewer mode: separately instantiated, read-only subagent with no shared task
history and no file-edit permission. The reviewer reported that no files were
modified.

## Independently Confirmed Findings

1. `BALF_BLOOD_COPD` is not an independent external test because its
   specificity and F1 influenced Round 12 versus Round 14 promotion.
2. V2 selected its final algorithm using the cohort named validation, making
   the published validation F1 a selection estimate.
3. V2 and V3 train/tune splitting randomized cells within source-label strata
   rather than grouping donor, library, or clonotype.
4. The external `real_gdt` rule used CD4 annotation/expression while CD4 is a
   model feature. Restoring 98 expression-excluded TCR-confirmed positives gave
   corrected F1 of 0.8929 for Round 12 and 0.8922 for Round 14.
5. Sorted-source positives are source-defined and may contain contamination,
   doublets, low-quality cells, or source-specific expression.
6. V3 candidate selection used full-atlas performance and composition criteria,
   even though the atlas contains training cells.
7. Full-atlas precision covers labeled cells only; hidden false positives are
   extrapolated under an untested cross-source exchangeability assumption.
8. The V3 manifest, embedded pickle fields, retained completed-run summary, and
   training composition do not describe one coherent promoted build record.
9. Every frozen profile exceeded configured extension false-positive targets;
   screen loading also removed the configured sealed-holdout status and allowed
   zero minimum model-gene coverage.
10. Canonical comparisons lack donor/dataset-clustered confidence intervals and
    dataset-macro promotion criteria.
11. Dataset imbalance and individual V/J features create source, repertoire,
    and clonotype portability risks not measured for the released model.
12. Round 14 is an uncalibrated probability-like gated score, so its numerical
    value is not a calibrated posterior probability.

## Reviewer Recommendation

Use hierarchical expression-independent labels and nested leave-one-dataset-out
evaluation, with donor/library/clonotype grouping inside outer-training data.
Balance sources, perform fold-local feature selection and calibration, compare
simple model families under identical folds, quantify feature stability and
missing-panel robustness, and report dataset-macro/worst-dataset metrics with
clustered uncertainty and prevalence-aware PPV/FDR.

## Hard Limit

No uninspected gdT-positive cohort remains. Existing data can support rigorous
nested cross-study evaluation and stress testing, but cannot establish unbiased
external sensitivity, PPV, universal calibration, true full-atlas purity, or
definitive superiority of Round 14 over Round 12 or V2.

The main-agent reconciliation and evidence paths are in `issue_register.csv`
and the reader report under `gdT_prediction/gdtai_methodology_audit/`.
