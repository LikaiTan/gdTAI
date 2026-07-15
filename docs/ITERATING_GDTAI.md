# Iterating gdTAI

## Experiment isolation

Copy `configs/models/gdtai/experiment_template.toml` into a run-specific output
directory. Record the parent model, registry snapshot, Git commit, environment,
input fingerprint, ground-truth rules, and dataset-level splits before training.
Do not create ambiguous source files such as `_final` or `_revised`; use Git and
an experiment ID.

## Leakage controls

- Split by entire dataset or cohort, never by cells alone.
- Fit feature selection, scaling, calibration, and thresholds only on training
  and tuning cohorts.
- Keep external H5ADs outside all model selection.
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

Publish both `high_f1` and `high_purity` operating modes when justified. A new
model must significantly improve the pre-registered comparator without relying
on one hard lineage marker or leaking held-out cohorts into feature selection.

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

The registry flags a local gdTAI v3 workspace override: Git contains the
promoted Round 14 release, while the dirty workspace contains Round 12 at the
same path. Resolve that state explicitly before promoting another model.
