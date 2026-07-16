# Reorganization Record

Date: 2026-07-15 Asia/Hong_Kong

## Scope

This maintenance milestone reorganized code, configuration, pre-integration
data views, reports, and legacy material. It did not advance the scientific
workflow beyond post-Phase 4 review and did not mutate any H5AD.

## Inventory and checkpoint

- pre-reorganization Git HEAD: `21b5d60c5732d464e2a7b75da8a95c560ff9fd9c`
- checkpoint: `Integrated_dataset/logs/project_reorganization/checkpoints/pre_reorganization_20260715/`
- bounded filesystem inventory: 3,165 entries
- dataset registry: 65 datasets
- library registry: 1,785 rows
- active Phase 0 datasets: 33
- compatibility view: 261 symlinks
- final registry snapshot:
  `data/registry/snapshots/post_reorganization_20260715_final/`

The checkpoint includes Git state and patches, exact conda/pip environment,
filesystem inventory, H5AD structures, sizes, timestamps, and sampled SHA256
fingerprints.

## Applied moves

- 69 root scripts moved into active workflow groups
- 9 pre-integration/data-embedded Python scripts moved into active groups
- 13 root data/configuration files moved to lifecycle-specific locations
- 2 additional configuration/R helpers moved
- 117 superseded files or directories moved into the legacy archive
- 1 completed one-shot watcher moved from maintenance to the archive

All moves are represented by plans and result maps under
`maintenance/reorganization/`. No source matrix, FASTQ, TCR table, or H5AD was
deleted or physically relocated.

## Verification

- strict dataset and library validation: passed
- compatibility links: 261 of 261 correct
- Python compilation: passed
- standard-library structural tests: 8 passed
- representative commands from `/tmp`: passed after setting writable Numba and
  Matplotlib cache directories for the restricted runtime
- gdTAI v2 and v4 artifact checksums: passed
- gdTAI v3: Round 12/Round 14 conflict resolved after a checksum-pinned
  same-cohort comparison; Round 14 is canonical and Round 12 is preserved as a
  validated high-purity fallback
- gdTAI v2, v3, and v4 trusted pickle loading: passed
- `TNK_candidates.h5ad`: size, timestamp, and sampled SHA256 unchanged
- `TNK_cleaned.h5ad`: size, timestamp, and sampled SHA256 unchanged
- `integrated.h5ad`: size, timestamp, and sampled SHA256 unchanged
- historical extended five-million atlas H5AD: size, timestamp, and sampled
  SHA256 unchanged

During portability testing, the Phase 1 script initially ignored `--help` and
started writing temporary per-dataset subsets. The process was terminated, the
new 4.9 GB temporary directory was removed after confirming it was absent from
the checkpoint, and all H5AD fingerprints were reverified. The script now has
an argument parser so `--help` exits before any work.

## Model conflict resolution

The reorganization audit found that the canonical gdTAI v3 path contained a
dirty Round 12 artifact while Git and the promotion documentation identified
Round 14. The user authorized a new comparison rather than an unconditional
restore. `workflows/gdtai/compare_gdtai_v3_round12_vs_round14.py` pinned both
artifacts by SHA256, compared them on identical full-atlas, atlas-held-out, and
independent-external cohorts, and selected Round 14 by prespecified recall and
false-positive guardrails followed by mean F1.

The canonical path now contains Round 14 with SHA256
`16dedc0081da9b8487887341232bcf8c9c9403dd3bbd72e04cab43d4cd7b2e09`.
Round 12 remains available at
`Integrated_dataset/models/gdT_prediction_classifier/gdtai_v3_round12_vs_round14/round12_model.pkl`
with SHA256
`7373e79350f7db190c415b376b9763e31652754438ee8c5afd3853beb7b2ebc4`.

## Dataset-centered physical cutover

The user approved the physical pre-integration storage cutover on 2026-07-16.
The migration used 191 same-filesystem rename operations, so source files were
not copied or rewritten. Physical storage is now organized under:

- `data/datasets/<dataset_id>/` for per-dataset raw, interim, and processed data
- `data/shared/` for project-wide pre-integration artifacts
- `data/compat/` for historical path trees

The legacy top-level names `downloads`, `analysis_26GSE_V4`, and `newdata`
remain compatibility aliases. Validation confirmed:

- 65 canonical dataset directories
- 44 promoted `processed/current.h5ad` links
- 43 preselected H5AD baselines with unchanged inode, size, timestamp, sampled
  SHA256, and HDF5 dimensions
- 262 correct lifecycle compatibility links
- 1,200 file-registry rows with matching canonical and legacy device/inode
- 18 passing structural and integrity tests
- unchanged milestone-H5AD fingerprints

The ordered plan, journals, preflight, registry copies, and final validation
are under `data/registry/migrations/dataset_centered_20260716/`. The filesystem
checkpoint is under
`Integrated_dataset/logs/project_reorganization/checkpoints/pre_dataset_centered_migration_20260716/`.
The post-cutover registry snapshot is
`data/registry/snapshots/post_dataset_centered_migration_20260716/`.

## Subsequent BALF_BLOOD_COPD intake

On 2026-07-16, the independent four-library BALF/PBMC validation study was
added as `BALF_BLOOD_COPD`. Its complete 951 GB workspace was moved intact by
one same-filesystem rename to
`data/datasets/BALF_BLOOD_COPD/workspace/`. The prior absolute path remains a
compatibility symlink, so existing study scripts and RStudio paths still
resolve the same files.

The selected H5AD contains 46,273 cells and 15,681 genes and is exposed at
`data/datasets/BALF_BLOOD_COPD/processed/current.h5ad`. The dataset and all
four library rows are validation-only and inactive for integration. Twelve
sentinel files, the workspace inode, 264 lifecycle links, strict registries,
and the selected H5AD identity all validated without failures.
