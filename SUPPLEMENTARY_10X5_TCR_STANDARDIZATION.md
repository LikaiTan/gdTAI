# Supplementary 10x 5' TCR Intake Standardization

## Scope

This document records the standardized handling used for the six approved
supplementary 10x 5' datasets:

- `GSE179994`
- `GSE234069`
- `GSE235863`
- `GSE240865`
- `GSE287301`
- `GSE287541`

These datasets were processed outside the original 33-dataset registry and then
converted into per-GSE H5AD files under:

- `downloads/per_gse_h5ad_with_metadata/`

The implementation entrypoint is:

- `supplementary_10x5_phase01.py`

The shared TCR join rules are defined in:

- `TCR_INTEGRATION_SOP.md`

## Non-negotiable rules

- intake only approved `10x 5'` datasets
- do not include any `10x 3'` data
- for `GSE234069`, use only `downloads/GSE234069/suppl/10x_5/`
- if TCR is provided separately, integrate it according to `TCR_INTEGRATION_SOP.md` and the sample-aware workflow in `tcr_integration_workflow.json`
- metadata and TCR fields must be standardized to the same schema used by:
  - `analysis_26GSE_V4/outputs/harmonized_metadata_v4.csv`

## Per-dataset outputs

The standardized per-GSE H5AD outputs are:

- `downloads/per_gse_h5ad_with_metadata/GSE179994_with_tcr.h5ad`
- `downloads/per_gse_h5ad_with_metadata/GSE234069_with_tcr.h5ad`
- `downloads/per_gse_h5ad_with_metadata/GSE235863_with_tcr.h5ad`
- `downloads/per_gse_h5ad_with_metadata/GSE240865_with_tcr.h5ad`
- `downloads/per_gse_h5ad_with_metadata/GSE287301_with_tcr.h5ad`
- `downloads/per_gse_h5ad_with_metadata/GSE287541_with_tcr.h5ad`

The supplementary candidate milestone and harmonized metadata outputs are:

- `Integrated_dataset/TNK_candidates_supp.h5ad`
- `analysis_26GSE_V4/outputs/harmonized_metadata_supp.csv`

## Canonical metadata fields

The supplementary lane writes metadata compatible with the main harmonized
header. Core columns are:

- `gse_id`
- `cell_id`
- `source_h5ad`
- `source_root`
- `library_id`
- `sample_id`
- `donor_id`
- `age`
- `sex`
- `tissue`
- `condition`
- `treatment`
- `enrichment_strategy`
- `assay_type`
- `tcr_availability`
- `original_cell_annotation`
- `sample_type`
- `donor_patient`
- `technology_simple`
- `tcr_vdj_flag`
- `barcode`

Rules:

- `technology_simple` must be filled as `10x 5'`
- `gse_id` is the only canonical dataset accession field
- keep one library-level identifier in `library_id`
- `sample_id` and `donor_id` are filled from the best available dataset-specific source
- `original_cell_annotation` should preserve the author-provided label when one exists

## Canonical TCR fields

The supplementary lane standardizes TCR outputs into these columns:

- `TCRseq`
- `TRA_cdr3`
- `TRA_v`
- `TRA_d`
- `TRA_j`
- `TRA_cdr3_nt`
- `TRA_clone_id`
- `TRA_umis`
- `TRA_reads`
- `TRB_cdr3`
- `TRB_v`
- `TRB_d`
- `TRB_j`
- `TRB_cdr3_nt`
- `TRB_clone_id`
- `TRB_umis`
- `TRB_reads`

Rules:

- use chain-specific `TRA_*` and `TRB_*` columns
- normalize `TCRseq` to `yes` or `no`
- if no TCR evidence exists for a cell, leave detailed chain columns blank
- build one canonical per-cell TCR table before merge
- normalize the join key to `sample_id + barcode_core`
- never join TCR by barcode alone
- replace stale TCR fields rather than mixing old and rebuilt values
- do not rely on row order

## Cell identifier rules

- preserve the original cell barcode when possible in `barcode`
- if one dataset needs a unique compound cell identifier, build it deterministically
- `GSE240865` required rewriting cell IDs to `library_id:barcode` so that `obs_names` are unique

## Dataset-specific notes

### `GSE234069`

- mixed `3'` and `5'` study
- only the `10x_5` lane is in scope
- `3'` content was isolated under `downloads/GSE234069/suppl/10x_3/` to prevent accidental intake

### `GSE240865`

- raw GEX came from standard 10x matrices
- TCR was standardized and merged into the same cell-level schema
- duplicate cell IDs were resolved by using `library_id:barcode`

### `GSE179994`, `GSE234069`, `GSE287541`

- source data required RDS/Seurat object export rather than direct matrix-triplet loading
- the helper used for this path is:
  - `supplementary_export_rds_payloads.R`

### `GSE287301`

- raw data included tar sidecar artifacts such as AppleDouble files
- those non-biological files must be skipped during import

## QC expectations

The supplementary lane should emit a separate QC package before any merge into
the main milestone object. Required outputs include:

- `Integrated_dataset/tables/supplementary_10x5/phase0_dataset_audit.csv`
- `Integrated_dataset/tables/supplementary_10x5/phase0_category_summary.csv`
- `Integrated_dataset/tables/supplementary_10x5/phase1_categoryA_selection_summary.csv`
- `Integrated_dataset/tables/supplementary_10x5/phase1_categoryA_marker_availability.csv`
- `Integrated_dataset/logs/supplementary_10x5/phase0_qc_summary.md`
- `Integrated_dataset/logs/supplementary_10x5/phase1_qc_summary.md`

The supplementary merge into the main milestone is only allowed after user QC
approval.

## Reuse guidance

If a future dataset is proposed for this lane, it should satisfy all of the
following before reuse:

- explicit evidence of `10x 5'`
- GEX raw counts or a recoverable filtered count object
- TCR data already embedded or recoverable from separate matrix/contig files
- enough metadata to populate the harmonized schema above

If any of these fail, do not silently force the dataset into the supplementary
10x 5' lane.
