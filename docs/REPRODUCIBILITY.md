# Reproducibility Contract

Every integration or model-training run must be reconstructable from five
identifiers: Git commit, environment export, dataset snapshot, input H5AD
fingerprint, and exact command/configuration.

## Before a run

```bash
PY=/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
$PY workflows/maintenance/manage_dataset_registry.py validate --strict-libraries
$PY workflows/maintenance/manage_dataset_registry.py snapshot <run_id>
$PY workflows/maintenance/build_data_compatibility_view.py --verify
```

Commit code and configuration before executing a new scientific run. Record
the commit and snapshot ID in the run log. Export the exact conda and pip state
to the run output directory.

## During a run

- Read raw or selected processed inputs; never rewrite raw data.
- Write into a new run-specific output directory or a temporary milestone.
- Preserve sparse matrices and respect the project memory envelope.
- Record arguments, random seeds, feature lists, splits, thresholds, and wall
  time in machine-readable manifests.
- Keep dataset-level validation cohorts isolated from feature selection,
  training, calibration, and threshold tuning.

## Before promotion

- Validate expected cell and gene counts, schemas, and metadata coverage.
- Compare output fingerprints with the recorded inputs.
- Run `python -m unittest discover -s tests -v` and
  `workflows/maintenance/verify_reorganization.py`.
- Update `TNK_PHASES_0_4_SCRIPT.md` with the exact entry point and outputs.
- Update the model registry for model runs.
- Obtain the review required by the runbook before replacing a milestone.

The current H5ADs predate the new snapshot convention and were not modified to
inject provenance into `adata.uns`. Future rebuilt milestones should record the
registry snapshot ID and code commit there at creation time.
