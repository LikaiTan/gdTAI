# Dataset-Centered Storage Migration

## Goal

Make `data/datasets/<dataset_id>/` the physical home for pre-integration
single-cell data while preserving every historical path used by existing code.

Canonical per-dataset layout:

```text
data/datasets/<dataset_id>/
├── dataset.json
├── raw/
│   ├── legacy_source/
│   └── source/
├── interim/
│   └── scanpy_project/
└── processed/
    ├── artifacts/
    ├── archive/
    ├── quarantine/
    └── current.h5ad
```

Project-wide legacy products are stored under `data/shared/`. Compatibility
trees are stored under `data/compat/`.

## Compatibility Contract

After cutover:

- `downloads` is a symlink to `data/compat/downloads`
- `analysis_26GSE_V4` is a symlink to
  `data/compat/analysis_26GSE_V4`
- `newdata` is a symlink to `data/compat/newdata`
- every moved child at an old path is a relative symlink to its canonical
  dataset-centered location
- old scripts and notebooks therefore continue to resolve the same file inode
- new workflows use the dataset registry and `ProjectPaths`
- a sole unambiguous standalone dataset H5AD is promoted to
  `processed/current.h5ad` during finalization

No H5AD is rewritten. The migration uses same-filesystem `rename(2)` operations.

## Commands

```bash
PY=/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python

$PY workflows/maintenance/migrate_dataset_storage.py plan
$PY workflows/maintenance/migrate_dataset_storage.py preflight
$PY workflows/maintenance/migrate_dataset_storage.py apply --confirm
$PY workflows/maintenance/migrate_dataset_storage.py finalize --confirm
$PY workflows/maintenance/migrate_dataset_storage.py validate
```

Locate one dataset:

```bash
$PY workflows/maintenance/manage_dataset_registry.py locate GSE144469
```

Rollback the physical layout and pre-migration registries:

```bash
$PY workflows/maintenance/migrate_dataset_storage.py rollback --confirm
```

## Safety Evidence

The migration state is recorded under:

`data/registry/migrations/dataset_centered_20260716/`

It contains:

- ordered forward/reverse move plan
- an explicit no-change preflight result
- source device and inode for every moved object
- selected-H5AD sampled hashes and HDF5 dimensions
- apply and rollback journals
- pre-migration configuration copies
- final validation result

The broader filesystem and Git checkpoint is:

`Integrated_dataset/logs/project_reorganization/checkpoints/pre_dataset_centered_migration_20260716/`

The validated post-cutover registry snapshot is:

`data/registry/snapshots/post_dataset_centered_migration_20260716/`

## Validated Result

- 191 physical rename-and-link operations
- 65 canonical dataset directories
- 44 promoted `processed/current.h5ad` links
- 262 correct lifecycle compatibility links
- 1,200 canonical/legacy file-registry inode matches
- zero failures in `validation.json`

## Subsequent External Cohort Intake

The independent `BALF_BLOOD_COPD` study was added after the base migration
using the same rename, compatibility-link, registry-snapshot, and validation
contract. At the user's request, a second journaled migration then restored
the raw/interim workspace to its original location and retained only the
selected validation H5AD physically inside the project.

The current storage state is:

- 66 canonical dataset directories
- 45 promoted `processed/current.h5ad` links
- 264 correct lifecycle compatibility links
- 1,212 file-registry rows
- approximately 949 GB of `BALF_BLOOD_COPD` raw/interim files at
  `/home/tanlikai/databank/owndata/singlecell`
- one 2,110,825,599-byte physical `BALF_BLOOD_COPD` H5AD under the project
- no `data/datasets/BALF_BLOOD_COPD/workspace/` directory

Its intake evidence is under
`data/registry/migrations/balf_blood_copd_20260716/`. The current H5AD-only
storage evidence is under
`data/registry/migrations/balf_blood_copd_h5ad_only_20260716/`.
