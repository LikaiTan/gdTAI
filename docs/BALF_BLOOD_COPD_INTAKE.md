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

The intake completed on 2026-07-16.

## Storage

The complete raw/interim study workspace remains physically at:

`/home/tanlikai/databank/owndata/singlecell`

The selected validation object is the only physical cohort data file retained
inside this project. Its artifact path is:

`data/datasets/BALF_BLOOD_COPD/processed/artifacts/phase4_final_annotated.h5ad`

The standard selected-object link is:

`data/datasets/BALF_BLOOD_COPD/processed/current.h5ad`

The original study workspace is exposed through the standard project raw
entry:

`data/datasets/BALF_BLOOD_COPD/raw/legacy_source`

The historical H5AD path remains valid as a compatibility link:

`/home/tanlikai/databank/owndata/singlecell/data/results/phase4_final_annotated.h5ad`

This preserves the cohort's existing scripts while avoiding duplication or
relocation of its FASTQ, BAM, matrix, TCR, and other large source files.

The storage refinement uses same-filesystem `rename(2)` operations. No H5AD or
source matrix is copied or rewritten.

## Commands

```bash
PY=/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python

$PY workflows/maintenance/migrate_balf_blood_copd_dataset.py plan
$PY workflows/maintenance/migrate_balf_blood_copd_dataset.py preflight
$PY workflows/maintenance/migrate_balf_blood_copd_dataset.py apply --confirm
$PY workflows/maintenance/migrate_balf_blood_copd_dataset.py finalize --confirm
$PY workflows/maintenance/migrate_balf_blood_copd_dataset.py validate

$PY workflows/maintenance/relocate_balf_blood_copd_raw_workspace.py plan
$PY workflows/maintenance/relocate_balf_blood_copd_raw_workspace.py preflight
$PY workflows/maintenance/relocate_balf_blood_copd_raw_workspace.py apply --confirm
$PY workflows/maintenance/relocate_balf_blood_copd_raw_workspace.py finalize --confirm
$PY workflows/maintenance/relocate_balf_blood_copd_raw_workspace.py validate
```

Rollback the storage refinement and restore the preceding
full-workspace-in-project layout:

```bash
$PY workflows/maintenance/relocate_balf_blood_copd_raw_workspace.py rollback --confirm
```

After rollback, revert the associated H5AD-only migration Git commit so the
registries and documentation match the restored filesystem state.

## Validation Evidence

- migration state:
  `data/registry/migrations/balf_blood_copd_20260716/`
- H5AD-only storage migration state:
  `data/registry/migrations/balf_blood_copd_h5ad_only_20260716/`
- pre-intake registry snapshot:
  `data/registry/snapshots/pre_balf_blood_copd_intake_20260716/`
- post-intake registry snapshot:
  `data/registry/snapshots/post_balf_blood_copd_intake_20260716/`
- pre-H5AD-only registry snapshot:
  `data/registry/snapshots/pre_balf_blood_copd_h5ad_only_20260716/`
- post-H5AD-only registry snapshot:
  `data/registry/snapshots/post_balf_blood_copd_h5ad_only_20260716/`
- 12 fingerprinted sentinels passed
- raw/interim workspace inode `24644876` was preserved
- old and canonical selected-H5AD paths resolve inode `25179089`
- selected H5AD size is `2,110,825,599` bytes
- the original raw/interim workspace occupies approximately 949 GB; the
  project-local cohort footprint is approximately 2.0 GB
- strict registry validation passed with 66 datasets and 33 active datasets
- all 264 lifecycle compatibility links passed

## Validation Role

The registry role is `gdTAI_independent_external_validation`. The dataset and
its four library rows remain inactive for Phase 0 and atlas integration. This
prevents accidental leakage into training or integration while keeping the
validation input discoverable and reproducible.
