# TNK Phase 0-4 Canonical Workflow

## Purpose

This file is the canonical executable workflow for:

- Phase 0: dataset audit and eligibility triage
- Phase 1: coarse T/NK extraction
- Phase 1b: conservative first cleanup
- Phase 1c: merged metadata backup and replacement
- Phase 2: merged cleanup
- Phase 3: scVI integration
- Phase 4: TRAB/TRB/TRD scoring

This file describes the standard workflow only.
Rules and active exceptions live in `TNK_PIPELINE_RUNBOOK.md`.
Current milestone state lives in `TNK_PIPELINE_STATUS.md`.

## Canonical environment

Use:

- `/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python`

## Canonical inputs and outputs

Inputs:

- `h5ad.csv`
- dataset files listed in `h5ad.csv`

Canonical output root:

- `Integrated_dataset/`

Canonical milestone H5AD files:

- `Integrated_dataset/TNK_candidates.h5ad`
- `Integrated_dataset/TNK_cleaned.h5ad`
- `Integrated_dataset/integrated.h5ad`

Allowed supplementary milestone exception:

- `Integrated_dataset/TNK_candidates_supp.h5ad`
  - used only for the approved supplementary 10x 5' intake lane

## Shared references

- TCR integration:
  - `TCR_INTEGRATION_SOP.md`
  - `tcr_integration_workflow.json`
- tissue harmonization:
  - `TISSUE_CORRECTION_WORKFLOW.md`
- disease-status harmonization:
  - `disease_status_correction_workflow.py`

## Standard helper scripts

- `phase0_dataset_audit.py`
- `repair_h5ad_from_raw.py`
- `phase1_extract_tnk_candidates.py`
- `phase1_finalize_from_temp.py`
- `phase1b_conservative_cleanup.py`
- `phase1c_replace_harmonized_metadata.py`
- `phase2_merged_cleanup.py`
- `phase3_scvi_scanvi.py`
- `phase4_gdt_module_scoring.py`
- `plot_phase4_threshold_barplots.py`

## Phase 0: Dataset audit and eligibility triage

Objective:

- inspect every registry dataset before extraction

Core checks:

- matrix state
- raw-count eligibility
- metadata fields
- duplicated gene names
- donor/sample/library candidate columns
- presence of key TCR and lineage genes

Core outputs:

- `Integrated_dataset/tables/phase0_dataset_audit.csv`
- `Integrated_dataset/tables/phase0_category_summary.csv`
- `Integrated_dataset/logs/phase0_qc_summary.md`
- Phase 0 PNG QC figures

## Phase 1: Coarse T/NK extraction

Objective:

- build a high-recall merged T/NK candidate pool from approved Phase 0 inputs

Core outputs:

- `Integrated_dataset/TNK_candidates.h5ad`
- Phase 1 tables, logs, and PNG QC figures

## Phase 1b: Conservative first cleanup

Objective:

- remove only obvious non-T/NK cells and apply the user-requested low-detection
  gene filter

Core outputs:

- in-place update of `Integrated_dataset/TNK_candidates.h5ad`
- Phase 1b tables, logs, and PNG QC figures

## Phase 1c: Merged metadata backup and replacement

Objective:

- export merged `adata.obs`
- back up the prior harmonized metadata file
- replace the harmonized metadata target with validated merged metadata

Required join key:

- `project name`
- `sampleid`
- `barcodes`

Core outputs:

- backup metadata CSV
- merged-obs export
- replacement summary tables and log

## Phase 2: Merged cleanup

Objective:

- perform merged-context cleanup on the candidate milestone

Core outputs:

- `Integrated_dataset/TNK_cleaned.h5ad`
- Phase 2 tables, logs, and PNG QC figures

## Phase 3: scVI integration

Objective:

- integrate the cleaned object with scVI
- build latent representation, neighbors, Leiden, and UMAP
- optionally keep scANVI outputs as reference-only when approved

Standard behavior:

- exclude mitochondrial, ribosomal, and noncoding-like genes from HVGs when
  requested
- keep large H5AD execution on the mirrored SSD tree when active in the runbook
- keep logs, PNG figures, tables, and models on NFS

Core outputs:

- `Integrated_dataset/integrated.h5ad`
- Phase 3 tables, logs, and PNG QC figures

## Phase 4: TRAB/TRB/TRD scoring

Objective:

- compute package-faithful continuous TRA/TRB/TRD module scores on the
  integrated object

Standard behavior:

- score from a temporary normalize-plus-log1p view of count-space `X`
- write continuous score columns back into `integrated.h5ad`
- keep hard labels and downstream threshold summaries separate from the raw
  score columns

Core outputs:

- updated `Integrated_dataset/integrated.h5ad`
- Phase 4 tables, logs, and PNG QC figures

## Supplementary 10x 5' intake lane

Approved supplementary GSEs:

- `GSE179994`
- `GSE235863`
- `GSE240865`
- `GSE287301`
- `GSE234069`
- `GSE287541`

Rules:

- process only approved 10x 5' inputs
- never mix in 10x 3'
- for `GSE234069`, use only `downloads/GSE234069/suppl/10x_5/`
- if TCR is external, integrate it according to `TCR_INTEGRATION_SOP.md`
- write supplementary per-GSE H5ADs under `downloads/per_gse_h5ad_with_metadata/`
- use `Integrated_dataset/TNK_candidates_supp.h5ad` only as the supplementary
  milestone before approved merge into the main candidate object

## QC-gate note

Default rule:

- every phase transition requires user-reviewed QC and explicit approval

Do not define active exceptions here.
Active exceptions belong only in `TNK_PIPELINE_RUNBOOK.md`.
