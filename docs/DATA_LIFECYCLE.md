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

The additional dataset is `BALF_BLOOD_COPD`, an independent external
validation cohort. Its complete legacy workspace is physically stored under
`data/datasets/BALF_BLOOD_COPD/workspace/`, while the selected object is
available through the standard `processed/current.h5ad` path.

`configs/datasets/datasets.csv` records dataset membership and selected H5ADs.
`libraries.csv` records sample/library assay scope. `files.csv` records
canonical and legacy file locations. The storage index is
`data/registry/storage_index.csv`.

New datasets must enter through `data/datasets/<dataset_id>/`; do not place new
physical data under a legacy alias.
