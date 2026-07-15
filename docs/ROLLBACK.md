# Rollback Procedure

The reorganization is reversible because scientific files were not rewritten,
legacy data paths remain present, and every file move has a checksum map.

## Code and configuration rollback

Prefer a Git revert of the reorganization commit. Do not use `git reset --hard`
on a shared or dirty workspace.

For a pre-commit local rollback, run each applied plan in reverse order:

```bash
PY=/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
$PY workflows/maintenance/bootstrap_migrated_scripts.py --rollback
$PY workflows/maintenance/apply_script_migration.py \
  --plan maintenance/reorganization/retired_workflow_archive_plan.csv --rollback
$PY workflows/maintenance/apply_script_migration.py \
  --plan maintenance/reorganization/legacy_archive_plan.csv --rollback
$PY workflows/maintenance/apply_script_migration.py \
  --plan maintenance/reorganization/config_and_r_migration_plan.csv --rollback
$PY workflows/maintenance/apply_script_migration.py \
  --plan maintenance/reorganization/config_data_migration_plan.csv --rollback
$PY workflows/maintenance/apply_script_migration.py \
  --plan maintenance/reorganization/preintegration_script_migration_plan.csv --rollback
$PY workflows/maintenance/apply_script_migration.py --rollback
```

The migration tool refuses to overwrite an existing source or destination and
checks content hashes across each move.

## Registry rollback

Registry edits automatically create a pre-edit snapshot. Restore a validated
snapshot with:

```bash
$PY workflows/maintenance/manage_dataset_registry.py restore <run_id> --confirm
```

The command verifies snapshot hashes and creates another backup before changing
the live registry.

## Data-view rollback

The `data/` tree contains generated symlinks. Stopping use of that view is
sufficient; legacy paths were never removed. Do not recursively delete source
targets. Rebuild or verify links from the registry when needed.

## Scientific milestone rollback

The pre-reorganization checkpoint is under
`Integrated_dataset/logs/project_reorganization/checkpoints/`. It records Git
state, environment, file inventory, H5AD structure, sizes, timestamps, and
sampled hashes. Since this reorganization did not write any H5AD, a mismatch is
a stop condition and should be investigated before any restoration attempt.
