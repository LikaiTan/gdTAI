# Adding or Removing Datasets

Dataset membership is a registry operation, not a filesystem deletion.

## Add a dataset

1. Place immutable source files under an approved legacy source location or a
   staged `data/raw/` location.
2. Build and validate one standardized per-dataset H5AD without overwriting an
   existing selected H5AD.
3. Register it:

```bash
PY=/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
$PY workflows/maintenance/manage_dataset_registry.py register GSEXXXXXX \
  --source-type GEO \
  --raw-path downloads/GSEXXXXXX \
  --processed-h5ad downloads/per_gse_h5ad_with_metadata/GSEXXXXXX.h5ad \
  --integration-role selected_for_next_integration \
  --activate
```

4. Replace the generated placeholder in `configs/datasets/libraries.csv` with
   actual sample/library rows and RNA/TCR assay scope.
5. Run strict validation and rebuild the compatibility view.
6. Create a named registry snapshot for the planned integration.
7. Follow Phase 0 and all runbook QC gates; build new milestones rather than
   mutating the current integrated object.

## Remove a dataset from the next build

```bash
$PY workflows/maintenance/manage_dataset_registry.py deactivate GSEXXXXXX \
  --reason "documented scientific or technical reason"
```

This marks the dataset and its libraries inactive and preserves all files and
historical membership. Existing atlas milestones remain unchanged until a new
integration is approved and built.

## Required review

Before integration, verify source identity, count-space eligibility, gene
symbols, cell barcodes, donor/sample/library metadata, TCR assay scope, and
dataset-specific exclusions. The merged metadata backup/replacement gate still
applies before progressing beyond merged-candidate phases.
