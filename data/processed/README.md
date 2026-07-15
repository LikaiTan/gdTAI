# Processed per-dataset data

This level contains the standardized per-dataset H5AD selected for Phase 0.
The generated `per_dataset/<dataset_id>/current.h5ad` links currently resolve to
the validated legacy files listed in `configs/datasets/datasets.csv`.

A future replacement must be written under a new run identifier, validated,
and promoted atomically. Do not overwrite a selected H5AD in place.
