# Data Lifecycle

## Levels

1. `data/raw/`: immutable source downloads or user-supplied files.
2. `data/interim/`: reproducible extraction, conversion, Cell Ranger, metadata,
   barcode, and TCR-join products.
3. `data/processed/per_dataset/`: the selected standardized H5AD for each
   dataset before Phase 0.
4. `Integrated_dataset/`: merged milestones, scientific tables, figures, logs,
   and released model artifacts.

Raw files are never analysis output destinations. A processed H5AD is promoted
by registry update after validation; it is not overwritten in place.

## Current migration state

The `data/` tree is a compatibility view over established source locations.
Symlinks were chosen deliberately because moving terabytes of source material
would create unnecessary downtime and rollback risk. The legacy source paths
remain authoritative until each dataset receives a physical-migration audit.

`configs/datasets/datasets.csv` records dataset membership and selected H5ADs.
`libraries.csv` records sample/library assay scope. `files.csv` records bounded
file roles and intended lifecycle locations. The compatibility-link manifest is
`data/registry/compatibility_links.json`.

## Physical cutover policy

A later physical move must be performed one dataset at a time:

1. Freeze the dataset registry and file manifest.
2. Calculate full hashes for source files selected for movement.
3. Copy to a staging directory on the target filesystem.
4. Verify hashes and expected file counts.
5. Atomically update registry paths and compatibility links.
6. Run the per-dataset schema and Phase 0 checks.
7. Retain the original until an explicit deletion approval.

No physical raw-data cutover was performed during this reorganization.
