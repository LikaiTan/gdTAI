# gdTAI Methodology Audit

Status: completed 2026-08-06

This audit was performed by a separately instantiated read-only reviewer and
reconciled against a main-agent provenance inventory. No H5AD was opened or
modified.

## Claim Corrections

1. The V2 cohort named validation was used to select the final algorithm. Its
   reported F1 is a model-selection estimate, not an untouched test estimate.
2. `BALF_BLOOD_COPD` was used in V3 Round 12 versus Round 14 guardrails and
   mean-F1 promotion. It is a reused cross-study benchmark, not an independent
   external test.
3. Train/tune splitting was performed at cell level within source-label strata,
   not by donor, library, or clonotype.
4. The original external positive label excluded paired productive gdT/no-abT
   cells using CD4 annotation or expression while CD4 is a model feature.
5. The V3 pickle, manifest, completed-run summary, and documented training
   composition do not describe one coherent promoted build record, although the
   registered artifact checksums are correct.

## Corrected Existing-Data Benchmark

The retained external predictions were re-evaluated using expression-independent
TCR labels: 1,046 paired productive gdT/no productive abT positives and 33,117
strict paired productive abT/no productive gdT negatives.

| Model | Precision | Recall | Specificity | F1 | FP |
| --- | ---: | ---: | ---: | ---: | ---: |
| V2 high-purity | 0.9465 | 0.8289 | 0.9985 | 0.8838 | 49 |
| V3 Round 12 | 0.9784 | 0.8212 | 0.9994 | 0.8929 | 19 |
| V3 Round 14 | 0.9634 | 0.8308 | 0.9990 | 0.8922 | 33 |

Round 12 and Round 14 are effectively tied on F1 in this reused benchmark.
Round 12 is purer; Round 14 is slightly more sensitive. This re-analysis fixes
label circularity but cannot restore test independence.

## Best Next Work Without New Data

1. Freeze V2, Round 12, and Round 14 as comparators and commit the complete
   analysis specification before rerunning.
2. Use nested leave-one-dataset-out outer folds. Inside each outer-training set,
   group folds by donor/library and clonotype where available.
3. Use paired productive TCR-only labels as primary truth. Treat sorted sources,
   single-chain evidence, CD4/NK warnings, and low-CD3 cells as noisy or stress
   strata rather than silently changing labels with model features.
4. Apply dataset-balanced weights and compare matched positive/negative cells
   within datasets wherever possible.
5. Compare regularized logistic, histogram gradient boosting, and a soft
   two-stage T-lineage model under identical folds. Do all feature selection,
   calibration, and thresholding inside inner folds.
6. Ablate individual V/J genes against family sums and detection indicators;
   require feature sign/rank stability across held-out datasets.
7. Report dataset-macro and worst-dataset recall/F1, strict-NK and paired-abT
   FPR, clustered bootstrap intervals, and PPV/FDR at realistic prevalence.
8. Add abstention for missing critical genes or insufficient model-gene
   coverage.
9. Require one immutable release bundle containing data hashes, split IDs,
   command, Git commit, environment lock, seeds, metrics, and semantic checks
   across pickle, manifest, registry, and run summary.

## Deliverables

- Reader report: `gdT_prediction/gdtai_methodology_audit/index.html`
- PDF: `gdT_prediction/gdtai_methodology_audit/gdtai_methodology_audit_report.pdf`
- Independent reviewer record:
  `Integrated_dataset/logs/gdT_prediction/gdtai_methodology_audit/independent_reviewer_record.md`
- Full issue register:
  `Integrated_dataset/tables/gdT_prediction/gdtai_methodology_audit/issue_register.csv`
- Corrected benchmark:
  `Integrated_dataset/tables/gdT_prediction/gdtai_methodology_audit/corrected_external_pure_tcr_metrics.csv`
- Reproducible workflow:
  `workflows/gdtai/build_gdtai_methodology_audit.py`

The strongest defensible current claim is that gdTAI is a checksum-pinned,
internally developed research classifier with useful TCR-expression
discrimination and material, now documented, uncertainty from model selection
reuse and cross-study domain shift.
