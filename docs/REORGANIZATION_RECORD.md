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
- gdTAI v3: known dirty-workspace Round 12 override detected against the
  canonical Git Round 14 release hash
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

## Known issue retained

The canonical Git gdTAI v3 pickle and `PROMOTION_DECISION.md` identify Round 14
at threshold 0.936. Before reorganization, the dirty workspace had already
overwritten that path with Round 12 at threshold 0.5. The model registry records
both hashes and reports the workspace override. No model bytes or scientific
metrics were changed or committed to resolve it during this maintenance task.
