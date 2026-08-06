# Iterating gdTAI

## Experiment isolation

Copy `configs/models/gdtai/experiment_template.toml` into a run-specific output
directory. Record the parent model, registry snapshot, Git commit, environment,
input fingerprint, ground-truth rules, and dataset-level splits before training.
Do not create ambiguous source files such as `_final` or `_revised`; use Git and
an experiment ID.

## Leakage controls

- Use nested leave-one-dataset-out outer folds. Within each outer-training set,
  group inner folds by donor/library and by clonotype where available; never
  split only by cells.
- Fit feature selection, scaling, calibration, thresholds, and operating modes
  inside inner training/tuning folds only.
- Keep any cohort intended as an external test outside all algorithm,
  threshold, guardrail, and promotion decisions. Once inspected for model
  choice, relabel it as a reused cross-study benchmark.
- Define primary labels independently of expression features. Transcriptomic
  CD4/NK/low-CD3 warnings belong in stress strata, not in TCR ground-truth
  exclusion rules when those genes are model inputs.
- Balance or weight datasets so one large source cannot dominate fitting or
  pooled metrics.
- Evaluate gold and silver positives separately.
- Report per-dataset results so one large source cannot dominate conclusions.
- Use raw counts or a counts layer and the declared normalization; do not mix
  already normalized `X` with count-space inputs.

## Model evaluation

The data are highly imbalanced. Report PR-AUC, precision, recall, specificity,
F1, MCC, balanced accuracy, confusion counts, and prevalence-aware PPV. Audit
paired/any productive TCRAB false positives in appropriate libraries, NK false
positives, sorted-gdT recall, TRDC+TRDV- dropout cases, low-quality cells, and
dataset/source shifts.

Publish both `high_f1` and `high_purity` operating modes when justified. Compare
models by dataset-macro outer-fold metrics with dataset-bootstrap uncertainty,
not pooled cell F1 alone. A new model must improve the pre-registered comparator
consistently across outer folds without relying on one hard lineage marker or
leaking held-out cohorts into feature selection.

The 2026-08-06 independent methodology audit is published at:

- `gdT_prediction/gdtai_methodology_audit/index.html`
- `Integrated_dataset/logs/gdT_prediction/gdtai_methodology_audit/gdtai_methodology_audit_summary.md`

## Promotion

Candidate artifacts remain experimental until the user explicitly approves a
promotion. Promotion must update, in one reviewed change:

- model pickle and checksum
- feature table and normalization contract
- thresholds and mode metrics
- training/validation cohort definitions
- inference example
- model manifest and promotion decision
- `configs/models/gdtai/model_registry.csv`

The registry currently pins the promoted Round 14 artifact and the preserved
Round 12 fallback by SHA256. Re-run `tests/test_model_registry.py` and
`tests/test_gdtai_round_selection.py` before any promotion.
