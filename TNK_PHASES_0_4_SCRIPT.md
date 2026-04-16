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

## Script recording rule

For each major phase or task, this file should record:

- phase or task name
- exact `.py` script used
- key outputs

## Phase 0: Dataset audit and eligibility triage

Objective:

- inspect every registry dataset before extraction

Phase or task:

- Phase 0 dataset audit

Exact `.py` script:

- `phase0_dataset_audit.py`

Related helper when rescue is needed:

- `repair_h5ad_from_raw.py`

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

Phase or task:

- Phase 1 coarse T/NK extraction

Exact `.py` scripts:

- `phase1_extract_tnk_candidates.py`
- `phase1_finalize_from_temp.py`

Core outputs:

- `Integrated_dataset/TNK_candidates.h5ad`
- Phase 1 tables, logs, and PNG QC figures

## Phase 1b: Conservative first cleanup

Objective:

- remove only obvious non-T/NK cells and apply the user-requested low-detection
  gene filter

Phase or task:

- Phase 1b conservative first cleanup

Exact `.py` script:

- `phase1b_conservative_cleanup.py`

Core outputs:

- in-place update of `Integrated_dataset/TNK_candidates.h5ad`
- Phase 1b tables, logs, and PNG QC figures

## Phase 1c: Merged metadata backup and replacement

Objective:

- export merged `adata.obs`
- back up the prior harmonized metadata file
- replace the harmonized metadata target with validated merged metadata

Phase or task:

- Phase 1c merged metadata backup and replacement

Exact `.py` script:

- `phase1c_replace_harmonized_metadata.py`

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

Phase or task:

- Phase 2 merged cleanup

Exact `.py` script:

- `phase2_merged_cleanup.py`

Core outputs:

- `Integrated_dataset/TNK_cleaned.h5ad`
- Phase 2 tables, logs, and PNG QC figures

## Phase 3: scVI integration

Objective:

- integrate the cleaned object with scVI
- build latent representation, neighbors, Leiden, and UMAP
- optionally keep scANVI outputs as reference-only when approved

Phase or task:

- Phase 3 scVI integration and optional scANVI reference annotation

Exact `.py` script:

- `phase3_scvi_scanvi.py`

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

Phase or task:

- Phase 4 TRAB/TRB/TRD continuous scoring

Exact `.py` scripts:

- `phase4_gdt_module_scoring.py`
- `plot_phase4_threshold_barplots.py`

Standard behavior:

- score from a temporary normalize-plus-log1p view of count-space `X`
- write continuous score columns back into `integrated.h5ad`
- keep hard labels and downstream threshold summaries separate from the raw
  score columns

Core outputs:

- updated `Integrated_dataset/integrated.h5ad`
- Phase 4 tables, logs, and PNG QC figures

## Post-Phase-4 milestone curation

Objective:

- apply approved milestone-level removals or carve-outs after Phase 4 review

Phase or task:

- Post-Phase-4 milestone curation for approved GSE removal

Exact `.py` script:

- `remove_no_tcr_gene_gses_from_milestones.py`

Core outputs:

- `high_speed_temp/Integrated_dataset/No_TCR_Gene_GSEs.h5ad`
- rewritten milestone H5ADs without the approved target GSEs
- `Integrated_dataset/tables/no_tcr_gene_gse_removal_counts.csv`

## Post-Phase-4 downstream reporting refinement

Objective:

- extend approved downstream review packages with focused T/NK/γδ figures,
  tissue-distribution statistics, and refreshed HTML/PDF reports

Phase or task:

- Post-Phase-4 downstream reporting refinement on integrated milestones

Exact `.py` scripts:

- `plot_plus6_tcr_pairing_umap.py`
- `plot_plus6_sorted_gdt_umap.py`
- `plot_plus6_tnk_marker_umaps.py`
- `build_plus6_gdt_report_assets.py`
- `render_plus6_final_report.py`

Core outputs:

- `Integrated_dataset/figures/plus6/plus6_umap_paired_tcr_doublets.png`
- `Integrated_dataset/figures/plus6/plus6_umap_sorted_gdt_highlight.png`
- `Integrated_dataset/figures/plus6/plus6_umap_paired_tcr_sorted_gdt.png`
- `Integrated_dataset/figures/plus6/plus6_tnk_marker_umap_panel.png`
- `Integrated_dataset/tables/plus6/plus6_gdt_candidate_statistics.csv`
- `Integrated_dataset/tables/plus6/plus6_gdt_candidate_overlap_gt0p4.csv`
- `Integrated_dataset/tables/plus6/plus6_gdt_paired_gdtcr_by_tissue.csv`
- `Integrated_dataset/tables/plus6/plus6_gdt_three_criteria_by_tissue.csv`
- `Integrated_dataset/plus6_profile_report.md`
- `Integrated_dataset/plus6_profile_report.html`
- `Integrated_dataset/plus6_profile_report.pdf`
- `Integrated_dataset/logs/no_tcr_gene_gse_removal.md`

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

## External T-cell/TNK intake with sample-aware TCR integration

Objective:

- process user-supplied T-cell or TNK-subset datasets outside the main registry
- standardize metadata to the project schema
- introduce productive alpha-beta and/or gamma-delta TCR fields by
  `sample_id + barcode_core`

Phase or task:

- External dataset intake with sample-aware productive TCR integration

Exact `.py` script:

- `process_hra005041_tcr_intake.py`

Standard behavior:

- preserve the input H5AD cell universe
- normalize canonical metadata fields such as `project name`, `sampleid`,
  `sample_id`, `cell_id`, `barcodes`, `barcode`, and `barcode_core`
- introduce only productive alpha-beta TCR rows
- introduce only productive-like gamma-delta TCR rows with valid CDR3 amino-acid
  and nucleotide sequence
- join TCR only by `sample_id + barcode_core`, never barcode alone

Core outputs:

- `downloads/per_gse_h5ad_with_metadata/HRA005041_T_cells_subset.h5ad`
- `Integrated_dataset/tables/HRA005041_tcr_intake/HRA005041_tcr_join_summary.csv`
- `Integrated_dataset/tables/HRA005041_tcr_intake/HRA005041_tcr_sample_summary.csv`
- `Integrated_dataset/tables/HRA005041_tcr_intake/HRA005041_tcr_unmatched_summary.csv`
- `Integrated_dataset/logs/HRA005041_tcr_intake.md`
- use `Integrated_dataset/TNK_candidates_supp.h5ad` only as the supplementary
  milestone before approved merge into the main candidate object

## Standalone external Phase 4 scoring review

Objective:

- compute project-standard TRAB/TRB/TRD module scores on one standalone intake
  H5AD
- write the score columns back into that intake H5AD
- generate standalone score QC and paired-TCR scatter plots

Phase or task:

- Standalone Phase 4 scoring review for external intake H5ADs

Exact `.py` script:

- `phase4_score_single_h5ad.py`

Core outputs:

- updated `downloads/per_gse_h5ad_with_metadata/HRA005041_T_cells_subset.h5ad`
- `Integrated_dataset/tables/HRA005041_phase4/`
- `Integrated_dataset/figures/HRA005041_phase4/`
- `Integrated_dataset/logs/HRA005041_phase4.log`

## Sorted gdT Seurat intake lane

Objective:

- convert user-supplied Seurat RDS objects of sorted gdT cells into Scanpy
  H5ADs
- preserve raw RNA counts
- harmonize metadata and embedded productive-like gamma-delta TCR fields to the
  project schema
- add `Sorted_gdT = True`
- recompute Scanpy UMAP
- compute project-standard standalone Phase 4 `TRD/TRAB` scores
- export UMAPs plus raw `TRAB`-vs-`TRD` scatter colored by paired `TRG/TRD`

Phase or task:

- Sorted gdT Seurat intake with standalone Phase 4 scoring

Exact `.py` scripts:

- `process_sorted_gdt_rds_intake.py`
- `phase4_score_single_h5ad.py`

Standard behavior:

- convert Seurat `RNA` assay counts to H5AD without densifying
- preserve raw counts in the intake H5AD and normalize only on temporary copies
- harmonize canonical metadata fields such as `project name`, `sampleid`,
  `sample_id`, `library_id`, `cell_id`, `barcodes`, `barcode`, and
  `barcode_core`
- support embedded `TRG/TRD` metadata as first-class TCR fields
- keep only productive-like `TRG` and `TRD` chains with valid chain label plus
  both amino-acid and nucleotide CDR3
- set `Sorted_gdT = True`, `input_population = purified_gdt`, and
  `tcr_chain_mode = gd_only`

Core outputs:

- `newdata/Sorted_gdT/GDT_2020AUG_woCOV_sorted_gdt.h5ad`
- `newdata/Sorted_gdT/GDTlung2023july_7p_sorted_gdt.h5ad`
- `newdata/Sorted_gdT/MalteGDT_sorted_gdt.h5ad`
- `Integrated_dataset/tables/Sorted_gdT/`
- `Integrated_dataset/logs/Sorted_gdT/`
- `Integrated_dataset/figures/GDT_2020AUG_woCOV_phase4/`
- `Integrated_dataset/figures/GDTlung2023july_7p_phase4/`
- `Integrated_dataset/figures/MalteGDT_phase4/`

## Supplementary 10x 5' intake lane

Exact `.py` scripts:

- `supplementary_10x5_phase01.py`
- `validate_supplementary_10x5_schema.py`

Key outputs:

- `downloads/per_gse_h5ad_with_metadata/`
- `Integrated_dataset/TNK_candidates_supp.h5ad`
- supplementary tables, logs, and PNG QC figures

## QC-gate note

Default rule:

- every phase transition requires user-reviewed QC and explicit approval

Do not define active exceptions here.
Active exceptions belong only in `TNK_PIPELINE_RUNBOOK.md`.
