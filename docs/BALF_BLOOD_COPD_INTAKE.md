# BALF_BLOOD_COPD Dataset Intake

## Scope

`BALF_BLOOD_COPD` is the independent local validation cohort used by gdTAI.
It contains four 10x 5' libraries from paired BALF and PBMC samples:

- `LIB1`, `LIB2`: BALF
- `LIB3`, `LIB4`: PBMC
- six demultiplexed donors
- COPD and lung-nodule diagnoses
- harmonized TRA/TRB/TRG/TRD evidence

The final validation H5AD contains 46,273 cells and 15,681 genes.

## Storage

The complete original study workspace is moved intact to:

`data/datasets/BALF_BLOOD_COPD/workspace/`

The selected validation object is exposed at:

`data/datasets/BALF_BLOOD_COPD/processed/current.h5ad`

The complete study workspace is also exposed through the standard raw entry:

`data/datasets/BALF_BLOOD_COPD/raw/legacy_source`

The historical workspace path remains valid:

`/home/tanlikai/databank/owndata/singlecell`

It becomes a compatibility symlink to the canonical workspace, preserving the
study's relative paths and existing absolute-path scripts.

## Commands

```bash
PY=/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python

$PY workflows/maintenance/migrate_balf_blood_copd_dataset.py plan
$PY workflows/maintenance/migrate_balf_blood_copd_dataset.py preflight
$PY workflows/maintenance/migrate_balf_blood_copd_dataset.py apply --confirm
$PY workflows/maintenance/migrate_balf_blood_copd_dataset.py finalize --confirm
$PY workflows/maintenance/migrate_balf_blood_copd_dataset.py validate
```

Rollback the physical move and restore the pre-intake registries:

```bash
$PY workflows/maintenance/migrate_balf_blood_copd_dataset.py rollback --confirm
```

After rollback, revert the associated Git commit so active gdTAI defaults again
match the restored filesystem state.

## Validation Role

The registry role is `gdTAI_independent_external_validation`. The dataset and
its four library rows remain inactive for Phase 0 and atlas integration. This
prevents accidental leakage into training or integration while keeping the
validation input discoverable and reproducible.
