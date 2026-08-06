# Data Lifecycle

## Levels

1. `data/datasets/<dataset_id>/raw/`: immutable source downloads or
   user-supplied files.
2. `data/datasets/<dataset_id>/interim/`: reproducible extraction, conversion,
   Cell Ranger, metadata, barcode, and TCR-join products.
3. `data/datasets/<dataset_id>/processed/`: versioned artifacts and the
   selected standardized `current.h5ad` before Phase 0.
4. `Integrated_dataset/`: merged milestones, scientific tables, figures, logs,
   and released model artifacts.

Raw files are never analysis output destinations. A processed H5AD is promoted
by registry update after validation; it is not overwritten in place.

`data/raw`, `data/interim`, and `data/processed` are generated lifecycle views.
`data/shared` contains project-wide pre-integration material, and `data/compat`
contains historical path trees.

## Current storage state

The physical cutover completed on 2026-07-16 with 191 same-filesystem renames.
The migration did not copy or rewrite H5AD files. The historical `downloads`,
`analysis_26GSE_V4`, and `newdata` names remain supported compatibility
aliases.

Current registry state:

- 66 dataset directories
- 45 promoted per-dataset H5ADs
- 1,789 library rows
- 1,212 file rows
- 264 lifecycle compatibility links

`BALF_BLOOD_COPD` is an independent external validation cohort. Its
approximately 949 GB raw/interim workspace remains
physically under `/home/tanlikai/databank/owndata/singlecell`. The project raw
entry is a compatibility view of that workspace, while only the selected
2,110,825,599-byte validation H5AD is physically retained under
`data/datasets/BALF_BLOOD_COPD/processed/`.

`GSE305372` was excluded and retired on 2026-08-06. Its local dataset directory
was deleted, its active registry rows were removed, and its reproducibility
records are retained under
`archive/retired_experiments/GSE305372_external_application/`.

`configs/datasets/datasets.csv` records dataset membership and selected H5ADs.
`libraries.csv` records sample/library assay scope. `files.csv` records
canonical and legacy file locations. The storage index is
`data/registry/storage_index.csv`.

New datasets normally enter through `data/datasets/<dataset_id>/`. An external
physical workspace requires an explicit runbook exception, registry entries,
compatibility links, and validation evidence such as the approved
`BALF_BLOOD_COPD` H5AD-only layout.
