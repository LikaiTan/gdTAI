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

- `configs/datasets/integration_inputs.csv`
- dataset files selected in `configs/datasets/datasets.csv`

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
  - `configs/tcr/integration_workflow.json`
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

- `workflows/integration/phase0_dataset_audit.py`

Related helper when rescue is needed:

- `workflows/intake/repair_h5ad_from_raw.py`

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

## Extension cohort standalone Phase 0 lane

Objective:

- stage and build approved new GEX/TCR cohorts without merging them
- enforce disease-aware tumor specimen context without overwriting source fields
- stop at a user-reviewed standalone Phase 0 gate

Phase or task:

- extension cohort intake, build, and standalone Phase 0 review

Exact scripts:

- `workflows/intake/validate_extension_cohorts.py`
- `workflows/intake/stage_extension_cohorts.py`
- `workflows/intake/build_extension_h5ads.py`
- `workflows/intake/qc_extension_h5ads.py`
- `workflows/intake/audit_extension_tcr_source_chains.py`
- `workflows/metadata/review_extension_sample_sources.py`

Core inputs:

- `configs/datasets/extension_cohorts.csv`
- `configs/datasets/extension_libraries.csv`
- `configs/metadata/sample_source_refinement_rules.json`

Core outputs:

- `data/interim/extension_intake/built_h5ads_manifest.csv`
- `data/interim/extension_intake/built_h5ads/<cohort_id>.h5ad`
- `Integrated_dataset/logs/extension_intake/extension_phase0_qc_summary.md`
- `Integrated_dataset/tables/extension_intake/extension_phase0_cohort_summary.csv`
- `Integrated_dataset/figures/extension_intake/extension_phase0_cohort_overview.png`
- `Integrated_dataset/tables/extension_intake/extension_tcr_source_chain_audit.csv`
- `Integrated_dataset/logs/sample_source_refinement/extension_sample_source_review.md`
- `Integrated_dataset/tables/sample_source_refinement/extensions/extension_tumor_project_context_counts.csv`

Current gate semantics:

- all H5AD review is read-only
- the user-selected GSE169246 processed `TNK_cleaned.h5ad` is registered as a
  standalone input and has passed the extension T/NK filter audit; it remains
  inactive before merge or integration pending user review
- Tan et al. 2021 remains excluded because `GDT_2020AUG_woCOV` already
  represents that study
- no merge or integration may proceed without explicit user approval

### GSE169246 processed-source registration

Phase or task:

- select and register the user-specified processed GSE169246 T/NK object

Exact `.py` script:

- `workflows/maintenance/manage_dataset_registry.py`

Key input:

- `/home/leeaaron/Exhausted_GDT/TNBC_validation/Integrated_dataset/TNK_cleaned.h5ad`

Key outputs:

- `data/datasets/GSE169246/processed/artifacts/TNK_cleaned.h5ad`
- `data/datasets/GSE169246/processed/current.h5ad`
- `docs/GSE169246_PROCESSED_H5AD_INTAKE.md`

The selected source passed standalone metadata/TCR harmonization and filter QC.
It remains inactive until user review is complete.

## Extension cohort standalone Phase 1 T/NK filter

Objective:

- filter every new cohort independently to a high-recall T/NK candidate set
- preserve NK cells for downstream gdTAI false-positive evaluation
- stop at a user-reviewed standalone gate before merge or integration

Phase or task:

- extension cohort T/NK filtering and metadata audit

Exact `.py` script:

- `workflows/intake/filter_extension_tnk_cells.py`

Core inputs:

- `data/interim/extension_intake/built_h5ads/<cohort_id>.h5ad`
- `data/datasets/GSE169246/processed/current.h5ad`

Selection rules:

- retain every cell with productive TRA, TRB, TRG, or TRD evidence
- retain author-annotated T or NK cells
- otherwise require a core CD3/TCR-alpha-beta marker, a gamma-delta marker,
  or the canonical multi-marker NK rule
- do not accept `IL7R` or `LTB` alone as T-lineage evidence
- explicit non-T/NK annotations block marker-only retention but cannot remove
  a productive-TCR cell

Metadata rules:

- preserve exact RNA-series provenance in `source_gse_id`
- record the intake bundle separately in `extension_cohort_id`
- require complete sample, library, donor, tissue, and specimen-context fields
- require unique cell IDs and `library_id + barcode` keys
- recompute and validate canonical TCR logical flags from CDR3 fields
- distinguish the 78 GSE169246 source libraries from 51 biological samples
  using the preserved full-barcode suffix

Core outputs:

- `data/interim/extension_intake/tnk_filtered_h5ads_manifest.csv`
- `data/interim/extension_intake/tnk_filtered_h5ads/<cohort_id>.h5ad`
- `Integrated_dataset/logs/extension_intake/extension_tnk_filter_qc.md`
- `Integrated_dataset/logs/extension_intake/extension_tnk_filter_qc.json`
- `Integrated_dataset/tables/extension_intake/tnk_filter/`
- `Integrated_dataset/figures/extension_intake/tnk_filter/extension_tnk_filter_overview.png`

Gate semantics:

- outputs remain cohort-separated and inactive for integration
- explicit user approval is required before merge, integration, or gdTAI
  evaluation

## Extension official-GEO metadata reconciliation

Objective:

- compare unresolved or over-specific extension metadata against official GEO
  series, sample, matrix, SOFT, and supplementary records
- preserve local source fields while recording separately supported values and
  unresolved ambiguities

Phase or task:

- read-only extension metadata evidence reconciliation

Exact `.py` script:

- `workflows/metadata/audit_extension_geo_metadata.py`

Core input:

- `configs/metadata/extension_geo_metadata_reconciliation.csv`
- `data/interim/extension_intake/tnk_filtered_h5ads_manifest.csv`

Core outputs:

- `docs/EXTENSION_GEO_METADATA_RECONCILIATION.md`
- `Integrated_dataset/tables/geo_metadata_reconciliation/`
- `Integrated_dataset/logs/geo_metadata_reconciliation/extension_geo_metadata_reconciliation_audit.md`

Standard behavior:

- use only direct official NCBI GEO evidence URLs
- validate all declared local values against backed read-only H5AD metadata
- preserve GEO-supported, GEO-partial, and GEO-unavailable states separately
- report scoped affected-cell counts without changing source metadata
- verify source-accession coverage and unchanged H5AD size/mtime during audit
- require additive derived columns for any later correction; never overwrite
  original metadata

Current gate semantics:

- the audit is complete and read-only
- GSE169246 `_b` libraries require a reviewed additive blood correction before
  any extension merge or integration
- metadata writeback, merge, and integration remain blocked pending explicit
  user approval

## Extension frozen gdTAI negative-control screen

Objective:

- apply every registered V2/V3 operating profile to the independently filtered
  extension cohorts without mutating or merging H5AD files
- quantify alpha-beta T-cell and NK false positives without treating no-TCR
  candidates as positive truth

Phase or task:

- checksum-pinned frozen-profile extension-cohort negative-control evaluation

Exact `.py` script:

- `workflows/gdtai/compare_frozen_gdtai_profiles.py`

Execution command:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/compare_frozen_gdtai_profiles.py \
  --stage screen --chunk-size 50000
```

Core inputs:

- `data/interim/extension_intake/tnk_filtered_h5ads_manifest.csv`
- `configs/models/gdtai/extension_evaluation.json`
- `configs/models/gdtai/model_registry.csv`

Core outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_extension_evaluation/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_extension_evaluation/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_extension_evaluation/`
- `gdT_prediction/gdtai_extension_evaluation/index.html`
- `gdT_prediction/gdtai_extension_evaluation/gdtai_extension_evaluation_report.pdf`

Standard behavior:

- verify every model artifact and input H5AD against its registered SHA-256
- use `log1p(raw counts per 10,000)` in bounded chunks
- accept productive single TRA or single TRB evidence as an alpha-beta
  negative control when productive TRG/TRD evidence is absent
- report paired and single-chain alpha-beta controls separately
- define strict NK controls by NK annotation/expression plus absence of all
  productive TRA/TRB/TRG/TRD evidence and doublet flags
- report no-productive-TCR calls as unevaluable candidates, not true positives
- calculate per-cohort, per-source, per-library, and per-donor rates with
  Wilson 95% intervals and descriptive exact paired tests
- run a feature-coverage sensitivity analysis that separates schema-warning
  cohorts such as GSE169246
- do not calculate new-cohort recall/F1/AUC or select/promote a model when the
  screen contains no unbiased gdT positive truth cells
- keep merge and integration blocked after the screen

## Post-Phase-4 gdTAI methodology audit

Objective:

- independently review data preparation, ground truth, features, algorithms,
  splitting, tuning, evaluation, deployment claims, and release provenance
- re-evaluate retained cross-study predictions with expression-independent
  TCR-only labels
- define the strongest defensible next program without adding data

Phase or task:

- read-only independent gdTAI training and evaluation audit

Exact `.py` script:

- `workflows/gdtai/build_gdtai_methodology_audit.py`

Core outputs:

- `docs/GDTAI_METHODOLOGY_AUDIT.md`
- `Integrated_dataset/tables/gdT_prediction/gdtai_methodology_audit/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_methodology_audit/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_methodology_audit/`
- `gdT_prediction/gdtai_methodology_audit/index.html`
- `gdT_prediction/gdtai_methodology_audit/gdtai_methodology_audit_report.pdf`

Standard behavior:

- do not open or mutate an H5AD
- verify registered model checksums and semantic agreement among pickle,
  manifest, registry, and run summary
- treat all inspected cohorts as development or stress benchmarks
- use paired productive TCR-only labels for the corrected primary benchmark;
  keep CD4/NK/low-CD3 expression as stress strata rather than truth rules
- require nested leave-one-dataset-out evaluation with donor/library/clonotype
  grouping before the next model-selection or promotion decision

## Phase 1: Coarse T/NK extraction

Objective:

- build a high-recall merged T/NK candidate pool from approved Phase 0 inputs

Phase or task:

- Phase 1 coarse T/NK extraction

Exact `.py` scripts:

- `workflows/integration/phase1_extract_tnk_candidates.py`
- `workflows/integration/phase1_finalize_from_temp.py`

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

- `workflows/integration/phase1b_conservative_cleanup.py`

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

- `workflows/integration/phase1c_replace_harmonized_metadata.py`

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

- `workflows/integration/phase2_merged_cleanup.py`

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

- `workflows/integration/phase3_scvi_scanvi.py`

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

- `workflows/integration/phase4_gdt_module_scoring.py`
- `workflows/integration/plot_phase4_threshold_barplots.py`

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

- `workflows/maintenance/remove_no_tcr_gene_gses_from_milestones.py`

Core outputs:

- `high_speed_temp/Integrated_dataset/No_TCR_Gene_GSEs.h5ad`
- rewritten milestone H5ADs without the approved target GSEs
- `Integrated_dataset/tables/no_tcr_gene_gse_removal_counts.csv`

## Post-Phase-4 extended five-million-atlas reporting refinement

Objective:

- extend approved downstream review packages with focused T/NK/γδ figures,
  tissue-distribution statistics, and refreshed HTML/PDF reports

Phase or task:

- Post-Phase-4 downstream reporting refinement on integrated milestones

Exact `.py` scripts:

- `workflows/reporting/plot_extended_atlas_tcr_pairing_umap.py`
- `workflows/reporting/plot_extended_atlas_sorted_gdt_umap.py`
- `workflows/reporting/plot_extended_atlas_tnk_marker_umaps.py`
- `workflows/reporting/plot_extended_atlas_cell_type_umap.py`
- `workflows/reporting/build_extended_atlas_gdt_report_assets.py`
- `workflows/reporting/render_extended_atlas_report.py`

The historical `plus6` output basenames below are retained for provenance and
stable links; active entry points use the unambiguous `extended_atlas` name.

Core outputs:

- `Integrated_dataset/figures/plus6/plus6_umap_paired_tcr_doublets.png`
- `Integrated_dataset/figures/plus6/plus6_umap_sorted_gdt_highlight.png`
- `Integrated_dataset/figures/plus6/plus6_umap_paired_tcr_sorted_gdt.png`
- `Integrated_dataset/figures/plus6/plus6_tnk_marker_umap_panel.png`
- `Integrated_dataset/figures/plus6/plus6_umap_by_corrected_simple_annotation_100x120mm.png`
- `Integrated_dataset/tables/plus6/plus6_gdt_candidate_statistics.csv`
- `Integrated_dataset/tables/plus6/plus6_gdt_candidate_overlap_gt0p4.csv`
- `Integrated_dataset/tables/plus6/plus6_gdt_paired_gdtcr_by_tissue.csv`
- `Integrated_dataset/tables/plus6/plus6_gdt_three_criteria_by_tissue.csv`
- `Integrated_dataset/plus6_profile_report.md`
- `Integrated_dataset/plus6_profile_report.html`
- `Integrated_dataset/plus6_profile_report.pdf`
- `Integrated_dataset/logs/no_tcr_gene_gse_removal.md`

## Post-Phase-4 annotation evidence review

Objective:

- audit cluster annotations with all-cell marker signatures, Phase 4 scores,
  TCR-chain metadata, reference labels as weak priors, and dataset/tissue
  composition
- produce review-only proposed labels and confidence flags without modifying
  milestone H5AD annotation columns

Phase or task:

- Post-Phase-4 cluster annotation evidence review

Exact `.py` script:

- `workflows/analysis/annotation_evidence_review.py`

Core outputs:

- `Integrated_dataset/tables/annotation_evidence_review/current_integrated_cluster_annotation_proposals.csv`
- `Integrated_dataset/tables/annotation_evidence_review/current_integrated_cluster_signature_scores.csv`
- `Integrated_dataset/tables/annotation_evidence_review/current_integrated_cluster_marker_gene_stats.csv`
- `Integrated_dataset/logs/annotation_evidence_review/current_integrated_annotation_evidence_summary.md`
- `Integrated_dataset/figures/annotation_evidence_review/current_integrated_cluster_signature_heatmap.png`
- `Integrated_dataset/tables/annotation_evidence_review/plus6_cluster_annotation_proposals.csv`
- `Integrated_dataset/tables/annotation_evidence_review/plus6_cluster_signature_scores.csv`
- `Integrated_dataset/tables/annotation_evidence_review/plus6_cluster_marker_gene_stats.csv`
- `Integrated_dataset/logs/annotation_evidence_review/plus6_annotation_evidence_summary.md`
- `Integrated_dataset/figures/annotation_evidence_review/plus6_cluster_signature_heatmap.png`

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

- `workflows/intake/process_hra005041_tcr_intake.py`

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

- `workflows/intake/phase4_score_single_h5ad.py`

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

- `workflows/intake/process_sorted_gdt_rds_intake.py`
- `workflows/intake/phase4_score_single_h5ad.py`

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

- `workflows/intake/supplementary_10x5_phase01.py`
- `workflows/intake/validate_supplementary_10x5_schema.py`

Key outputs:

- `downloads/per_gse_h5ad_with_metadata/`
- `Integrated_dataset/TNK_candidates_supp.h5ad`
- supplementary tables, logs, and PNG QC figures

## Post-Phase-4 gdT prediction package-style evaluation

Objective:

- evaluate package-style gdT score classifiers on the historical extended
  five-million-cell integrated object
- build project-specific primary gold and sensitivity-only silver labels from
  TCR metadata without modifying any H5AD
- render a static HTML/PDF report for downstream review

Phase or task:

- gdT prediction evaluation using shared package classifier logic

Exact `.py` scripts:

- `workflows/gdtai/run_gdt_prediction_package_evaluation.py`
- `workflows/gdtai/plot_gdt_truth_trd_trab_marker_scatter.py`
- `workflows/gdtai/run_tn_fn_trd_gt0p1_deg_reproducibility.py`
- `workflows/gdtai/run_trd_gt0p1_high_no_ab_vs_low_ab_deg_reproducibility.py`
- `workflows/gdtai/run_gdt_deg_tcr_classifier_training.py`
- `workflows/gdtai/run_gdt_gse144469_holdout_tcrgene_classifier.py`
- `workflows/gdtai/run_gdt_classifier_failure_mode_audit.py`
- `workflows/gdtai/predict_with_selected_gdt_model.py`
- `workflows/gdtai/build_gdtai_overview_report.py`
- `workflows/gdtai/run_gdtai_nkguard_classifier.py`
- `workflows/gdtai/run_gdtai_multiguard_classifier.py`
- `workflows/gdtai/run_gdtai_transcriptome_gate_cascade.py`
- `workflows/gdtai/run_gdtai_annotation_specific_cascade.py`
- `workflows/gdtai/run_gdtai_annotation_specific_tp_fn_audit.py`

Core outputs:

- `Integrated_dataset/tables/gdT_prediction/`
- `Integrated_dataset/figures/gdT_prediction/`
- `Integrated_dataset/logs/gdT_prediction/gdT_prediction_summary.md`
- `Integrated_dataset/figures/gdT_prediction/trd_vs_trab_gold_silver_marker_expression.png`
- `Integrated_dataset/tables/gdT_prediction/tn_fn_trd_gt0p1_global_deg.csv`
- `Integrated_dataset/tables/gdT_prediction/tn_fn_trd_gt0p1_deg_reproducibility.csv`
- `Integrated_dataset/figures/gdT_prediction/tn_fn_trd_gt0p1_global_deg_volcano.png`
- `Integrated_dataset/figures/gdT_prediction/tn_fn_trd_gt0p1_source_logfc_reproducibility.png`
- `Integrated_dataset/logs/gdT_prediction/tn_fn_trd_gt0p1_deg_reproducibility_summary.md`
- `Integrated_dataset/tables/gdT_prediction/trd_gt0p1_high_no_ab_vs_expanded_negative_no_tcr_genes_global_deg.csv`
- `Integrated_dataset/tables/gdT_prediction/trd_gt0p1_high_no_ab_vs_expanded_negative_no_tcr_genes_deg_reproducibility.csv`
- `Integrated_dataset/tables/gdT_prediction/trd_gt0p1_high_no_ab_vs_expanded_negative_no_tcr_genes_excluded_tcr_genes.csv`
- `Integrated_dataset/figures/gdT_prediction/trd_gt0p1_high_no_ab_vs_expanded_negative_no_tcr_genes_global_deg_volcano.png`
- `Integrated_dataset/logs/gdT_prediction/trd_gt0p1_high_no_ab_vs_expanded_negative_no_tcr_genes_deg_reproducibility_report.md`
- `gdT_prediction/trd_gt0p1_high_no_ab_vs_expanded_negative_no_tcr_genes_deg_reproducibility_report.html`
- `gdT_prediction/trd_gt0p1_high_no_ab_vs_expanded_negative_no_tcr_genes_deg_reproducibility_report.pdf`
- `Integrated_dataset/tables/gdT_prediction/classifier_training/classifier_metrics_validation.csv`
- `Integrated_dataset/tables/gdT_prediction/classifier_training/classifier_acceptance_vs_best_baseline.csv`
- `Integrated_dataset/tables/gdT_prediction/classifier_training/prevalence_aware_ppv_scenarios.csv`
- `Integrated_dataset/tables/gdT_prediction/classifier_training/selected_model_full_dataset_prediction_overall.csv`
- `Integrated_dataset/tables/gdT_prediction/classifier_training/selected_model_full_dataset_prediction_by_source.csv`
- `Integrated_dataset/logs/gdT_prediction/classifier_training/gdt_deg_tcr_classifier_report.md`
- `Integrated_dataset/models/gdT_prediction_classifier/selected_gdt_classifier.pkl`
- `gdT_prediction/gdt_deg_tcr_classifier_report.html`
- `gdT_prediction/gdt_deg_tcr_classifier_report.pdf`
- `Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene/validation_metrics.csv`
- `Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene/acceptance_vs_individual_gene_baseline.csv`
- `Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene/selected_model_full_dataset_prediction_overall.csv`
- `Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene/selected_model_full_dataset_false_positive_overall.csv`
- `Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene/selected_model_full_dataset_false_positive_by_source.csv`
- `Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene/selected_model_full_dataset_prediction_by_annotation.csv`
- `Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene/trd_vs_trab_prediction_scatter_sample.csv.gz`
- `Integrated_dataset/figures/gdT_prediction/gse144469_holdout_tcrgene/trd_vs_trab_prediction_method_agreement.png`
- `Integrated_dataset/figures/gdT_prediction/gse144469_holdout_tcrgene/trd_vs_trab_tcrgene_known_fp_status.png`
- `Integrated_dataset/logs/gdT_prediction/gse144469_holdout_tcrgene/gse144469_holdout_tcrgene_report.md`
- `Integrated_dataset/models/gdT_prediction_classifier/gse144469_holdout_tcrgene/selected_model.pkl`
- `gdT_prediction/gse144469_holdout_tcrgene_report.html`
- `gdT_prediction/gse144469_holdout_tcrgene_report.pdf`
- `Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene_failure_modes/`
- `Integrated_dataset/figures/gdT_prediction/gse144469_holdout_tcrgene_failure_modes/`
- `Integrated_dataset/logs/gdT_prediction/gse144469_holdout_tcrgene_failure_modes/gse144469_holdout_tcrgene_failure_mode_audit.md`
- `gdT_prediction/gse144469_holdout_tcrgene_failure_mode_audit.html`
- `gdT_prediction/gse144469_holdout_tcrgene_failure_mode_audit.pdf`
- `gdT_prediction/gdT_classifier_local_handoff_plan.md`
- `Integrated_dataset/tables/gdT_prediction/external_tests/<DATASET_ID>/`
- `Integrated_dataset/figures/gdT_prediction/external_tests/<DATASET_ID>/`
- `Integrated_dataset/logs/gdT_prediction/external_tests/<DATASET_ID>/`
- `gdT_prediction/external_tests/<DATASET_ID>/index.html`
- `Integrated_dataset/logs/gdT_prediction/gdTAI_overview_report.md`
- `Integrated_dataset/figures/gdT_prediction/gdTAI_overview/`
- `gdT_prediction/gdTAI_overview_report.html`
- `Integrated_dataset/tables/gdT_prediction/gdtai_nkguard/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_nkguard/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_nkguard/gdtai_nkguard_report.md`
- `Integrated_dataset/models/gdT_prediction_classifier/gdtai_nkguard/best_candidate_model.pkl`
- `gdT_prediction/gdtai_nkguard_report.html`
- `gdT_prediction/gdtai_nkguard_report.pdf`
- `Integrated_dataset/tables/gdT_prediction/gdtai_multiguard/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_multiguard/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_multiguard/gdtai_multiguard_report.md`
- `Integrated_dataset/models/gdT_prediction_classifier/gdtai_multiguard/best_candidate_model.pkl`
- `gdT_prediction/gdtai_multiguard_report.html`
- `gdT_prediction/gdtai_multiguard_report.pdf`
- `Integrated_dataset/tables/gdT_prediction/gdtai_transcriptome_gate/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_transcriptome_gate/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_transcriptome_gate/gdtai_transcriptome_gate_report.md`
- `Integrated_dataset/models/gdT_prediction_classifier/gdtai_transcriptome_gate/gdtai_transcriptome_gate_protocol.md`
- `gdT_prediction/gdtai_transcriptome_gate_report.html`
- `gdT_prediction/gdtai_transcriptome_gate_report.pdf`
- `Integrated_dataset/tables/gdT_prediction/gdtai_annotation_specific_cascade/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_annotation_specific_cascade/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_annotation_specific_cascade/gdtai_annotation_specific_cascade_report.md`
- `Integrated_dataset/models/gdT_prediction_classifier/gdtai_annotation_specific_cascade/selected_annotation_specific_wrapper.pkl`
- `Integrated_dataset/models/gdT_prediction_classifier/gdtai_annotation_specific_cascade/gdtai_annotation_specific_cascade_protocol.md`
- `gdT_prediction/gdtai_annotation_specific_cascade_report.html`
- `gdT_prediction/gdtai_annotation_specific_cascade_report.pdf`
- `Integrated_dataset/tables/gdT_prediction/gdtai_annotation_specific_tp_fn_audit/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_annotation_specific_tp_fn_audit/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_annotation_specific_tp_fn_audit/gdtai_annotation_specific_tp_fn_audit_report.md`
- `gdT_prediction/gdtai_annotation_specific_tp_fn_audit.html`
- `gdT_prediction/gdtai_annotation_specific_tp_fn_audit.pdf`
- `gdT_prediction/index.html`
- `gdT_prediction/gdT_prediction_report.pdf`

Phase or task:

- Post-Phase-4 gdTAI v3 TRDC/NK-guard classifier training, cross-study benchmark testing, and failure audit

Exact `.py` scripts:

- `workflows/gdtai/run_gdtai_v3_trdc_nk_guard_classifier.py`
- `workflows/gdtai/run_gdtai_v3_trdc_nk_failure_audit.py`

Core outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v3_trdc_nk_guard/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v3_trdc_nk_guard/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v3_trdc_nk_guard/`
- `gdT_prediction/gdtai_v3_trdc_nk_guard_report.html`
- `gdT_prediction/gdtai_v3_trdc_nk_guard_report.pdf`
- `gdT_prediction/gdtai_v3_trdc_nk_failure_audit.html`
- `gdT_prediction/gdtai_v3_trdc_nk_failure_audit.pdf`
- if promoted: `Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/gdTAI_v3_model.pkl` plus README, methodology, feature, metrics, and example files

Standard behavior:

- include `GSE144469` primary gold cells in atlas train/tune splits
- keep the non-GSE144469 validation cohorts from the prior multi-cohort validation lane: `GSE254249` paired-TCRAB/no-gdTCR negatives and `GDT_2020AUG_woCOV` cord-blood gdT positives
- keep the cross-study H5AD out of coefficient fitting; historical candidate
  and promotion workflows later inspected it for selection, so it is not an
  independent external test
- after training, apply the selected v3 candidate and reference strategies to the full atlas input H5AD and write whole-atlas source/tissue/annotation summaries plus selected predicted cells
- run at least five iteration rounds unless the recall/estimated-FP target is already reached; candidate families are not restricted to XGBoost
- select the final review candidate from full-atlas target metrics using recall `> 0.8` and estimated FP fraction `< 5%` of predictions, estimating hidden abT FP from paired-TCRAB rates in TCR-sequenced sources
- require cross-study H5AD inference from `layers["counts"]` and fail rather than using normalized `X`
- exclude NK+TCRAB-overlap cells from negative train/tune labels and explicit hard-negative construction
- restrict explicit NK hard negatives to expression `TRDC+TRDV-`; non-NK TCRAB hard negatives remain a separate guard set
- report ambiguous/mixed/partial external TCR labels as stress groups outside primary binary metrics

Phase or task:

- Post-Phase-4 gdTAI NK-optimized gdT atlas metadata harmonization and report

Exact `.py` script:

- `workflows/gdt_atlas/run_gdtai_nk_optimized_gdt_atlas.py`

Core outputs:

- `configs/gdt_atlas/metadata_rules.json`
- `Integrated_dataset/tables/gdT_atlas/`
- `Integrated_dataset/figures/gdT_atlas/`
- `Integrated_dataset/logs/gdT_atlas/gdtai_nk_optimized_gdt_atlas_summary.md`
- `gdT_atlas/index.html`
- `gdT_atlas/gdT_atlas_report.pdf`

Phase or task:

- Post-Phase-4 curated gdT phenotype atlas

Exact `.py` script:

- `workflows/gdt_atlas/run_gdt_atlas_curated_phenotype_analysis.py`

Core outputs:

- `configs/gdt_atlas/phenotype_marker_panels.json`
- `Integrated_dataset/tables/gdT_atlas_curated_phenotypes/`
- `Integrated_dataset/figures/gdT_atlas_curated_phenotypes/`
- `Integrated_dataset/logs/gdT_atlas_curated_phenotypes/gdt_atlas_curated_phenotypes_summary.md`
- `gdT_atlas/curated_phenotypes/index.html`
- `gdT_atlas/curated_phenotypes/gdT_atlas_curated_phenotypes_report.pdf`

## QC-gate note

Default rule:

- every phase transition requires user-reviewed QC and explicit approval

Do not define active exceptions here.
Active exceptions belong only in `TNK_PIPELINE_RUNBOOK.md`.

## Post-Phase-4 rigorous gdT atlas reanalysis

Objective:

- rebuild the gdT biology atlas from promoted gdTAI predictions with no
  silver-only FN add-back
- remove likely false positives, add back only gold-level FN cells, re-run
  gdT-only scVI integration, UMAP, Leiden clustering, phenotype annotation,
  sample-level statistics, and pseudobulk-style DE

Phase or task:

- Post-Phase-4 rigorous gdT atlas reanalysis and report

Exact `.py` script:

- `workflows/gdt_atlas/run_gdt_atlas_rigorous_reanalysis.py`

Core outputs:

- `Integrated_dataset/gdT_atlas/gdt_curated_integrated.h5ad`
- `Integrated_dataset/tables/gdT_atlas_rigorous/`
- `Integrated_dataset/figures/gdT_atlas_rigorous/`
- `Integrated_dataset/logs/gdT_atlas_rigorous/gdt_atlas_rigorous_summary.md`
- `gdT_atlas/rigorous/index.html`
- `gdT_atlas/rigorous/gdT_atlas_rigorous_report.pdf`

## gdTAI v3 decision-report refinement

Objective:

- replace the dense gdTAI v3 candidate report with a reader-facing decision report
  that emphasizes conclusions, compact figures, and promotion rationale

Phase or task:

- Post-Phase-4 gdTAI v3 report refinement

Exact `.py` script:

- `workflows/gdtai/render_gdtai_v3_clear_report.py`

Core outputs:

- `Integrated_dataset/logs/gdT_prediction/gdtai_v3_trdc_nk_guard/gdtai_v3_trdc_nk_guard_report.md`
- `gdT_prediction/gdtai_v3_trdc_nk_guard_report.html`
- `gdT_prediction/gdtai_v3_trdc_nk_guard_report.pdf`
- compact decision PNG figures under `gdT_prediction/assets/gdtai_v3_trdc_nk_guard/`

## gdTAI v3.0 promotion and predicted gdT atlas handoff

Objective:

- promote the user-approved gdTAI v3.0 model package
- remove previous gdT atlas H5AD/report output trees
- create a handoff H5AD containing promoted gdTAI v3 predicted gdT cells,
  primary-gold gdT FN add-back cells, and no known/likely paired-TCRAB FP cells

Phase or task:

- Post-Phase-4 gdTAI v3.0 promotion and predicted gdT-cell atlas handoff

Exact `.py` script:

- `workflows/gdtai/promote_gdtai_v3_and_build_predicted_gdt_atlas.py`

Core outputs:

- `Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/gdTAI_v3_model.pkl`
- `Integrated_dataset/gdT_atlas/predicted_gdt_cell_atlas.h5ad`
- `Integrated_dataset/gdT_atlas/predicted_gdt_cell_atlas.md`
- `Integrated_dataset/tables/gdT_atlas/predicted_gdt_cell_atlas_summary.csv`
- `Integrated_dataset/tables/gdT_atlas/predicted_gdt_cell_atlas_metadata.csv.gz`
- `Integrated_dataset/tables/gdT_atlas/predicted_gdt_cell_atlas_removed_fp_cells.csv.gz`
- `Integrated_dataset/tables/gdT_atlas/predicted_gdt_cell_atlas_gold_fn_added_cells.csv.gz`

## Post-Phase-4 predicted gdT atlas tumor-type audit

Objective:

- resolve broad tumor-like atlas tissue labels to original tumor types by source GSE and row-level tissue context without modifying the atlas H5AD

Phase or task:

- Predicted gdT atlas tumor-type source audit

Exact `.py` script:

- `workflows/gdt_atlas/run_gdt_atlas_tumor_type_audit.py`

Core outputs:

- `Integrated_dataset/tables/gdT_atlas/tumor_type_source_mapping.csv`
- `Integrated_dataset/tables/gdT_atlas/tumor_type_by_source_tissue.csv`
- `Integrated_dataset/tables/gdT_atlas/tumor_type_cell_annotation_review.csv.gz`
- `Integrated_dataset/logs/gdT_atlas/tumor_type_source_mapping.md`

Standard behavior:

- read `Integrated_dataset/gdT_atlas/predicted_gdt_cell_atlas.h5ad` in backed mode
- use `configs/gdt_atlas/metadata_rules.json` source-level tumor-type rules gated by row tissue values
- preserve raw tissue columns and keep tumor-type fields as review-only exports until write-back is explicitly approved

## Post-Phase-4 gdTAI model comparison

Objective:

- compare the latest local gdTAI model against gdTAI v2.0 high-purity at
  whole-atlas and per-dataset levels without mutating H5AD files

Phase or task:

- gdTAI latest-vs-v2 high-purity prediction overlap audit

Exact `.py` script:

- `workflows/gdtai/compare_gdtai_latest_vs_v2_high_purity.py`

Core outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v3_vs_v2_high_purity/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v3_vs_v2_high_purity/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v3_vs_v2_high_purity/`
- `gdT_prediction/gdtai_v3_vs_v2_high_purity_report.html`
- `gdT_prediction/gdtai_v3_vs_v2_high_purity_report.pdf`

## Post-Phase-4 gdTAI v3 Round 12 versus Round 14 revalidation

Objective:

- resolve the canonical-model conflict by comparing the actual Round 12 and
  Round 14 artifacts on identical cohorts without mutating H5AD files
- select one balanced default by prespecified recall and false-positive
  guardrails, while retaining the other artifact as a documented fallback

Phase or task:

- checksum-pinned gdTAI v3 Round 12 versus Round 14 comparison and promotion

Exact `.py` script:

- `workflows/gdtai/compare_gdtai_v3_round12_vs_round14.py`

Promotion command:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/compare_gdtai_v3_round12_vs_round14.py \
  --promote-selected
```

Core outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v3_round12_vs_round14/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v3_round12_vs_round14/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v3_round12_vs_round14/`
- `Integrated_dataset/models/gdT_prediction_classifier/gdtai_v3_round12_vs_round14/`
- `gdT_prediction/gdtai_v3_round12_vs_round14/index.html`
- `gdT_prediction/gdtai_v3_round12_vs_round14/gdtai_v3_round12_vs_round14_report.pdf`
- selected package synchronized to
  `Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/`

Standard behavior:

- require exact Round 12 and Round 14 SHA256 identities
- rerun both artifacts directly on the reused BALF_BLOOD_COPD cross-study benchmark
- verify direct predictions against the validated external cache and a
  deterministic full-atlas cache sample
- compare raw full-atlas calls, gold and silver recovery, paired-TCRAB and NK
  burdens, held-out cord-blood positives, and held-out GSE254249 negatives
- run descriptive exact paired tests and report every dataset separately
- historical selection passed guardrails and then maximized mean F1 across
  full-atlas gold, atlas-held-out, and the reused cross-study benchmark; this
  use exhausted benchmark independence and must not be presented as external
  validation
- preserve Round 12 as the high-purity fallback when Round 14 is promoted

## Post-Phase-4 gdTAI v4 two-stage retraining

Objective:

- train and evaluate an experimental two-stage gdTAI classifier with a soft
  T-lineage gate and dropout-tolerant gdT classifier

Phase or task:

- gdTAI v4 T-cell-gate retraining and full-atlas comparison

Exact `.py` script:

- `workflows/gdtai/run_gdtai_v4_tcell_gate_classifier.py`

Core outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_tcell_gate/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v4_tcell_gate/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_tcell_gate/`
- `Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v4.0/`
- `gdT_prediction/gdtai_v4_tcell_gate_report.html`
- `gdT_prediction/gdtai_v4_tcell_gate_report.pdf`

## Project maintenance: dataset-centered pre-integration storage

Objective:

- make `data/datasets/<dataset_id>/` the physical home for raw, interim, and
  processed per-dataset inputs
- preserve all historical `downloads/`, `analysis_26GSE_V4/`, and `newdata/`
  paths as compatibility symlinks to the same inodes
- support validated forward migration and reverse rollback without rewriting
  any H5AD

Phase or task:

- Dataset-centered physical storage migration with legacy-path compatibility

Exact `.py` scripts:

- `workflows/maintenance/migrate_dataset_storage.py`
- `workflows/maintenance/build_data_compatibility_view.py`
- `workflows/maintenance/manage_dataset_registry.py`

Core outputs:

- `data/datasets/<dataset_id>/`
- `data/shared/`
- `data/compat/`
- `data/registry/storage_index.csv`
- `data/registry/migrations/dataset_centered_20260716/`
- top-level compatibility aliases at `downloads`, `analysis_26GSE_V4`, and
  `newdata`

Standard behavior:

- run `plan` and `preflight` before `apply --confirm`
- use same-filesystem renames and verify source device/inode identities
- preserve selected-H5AD sampled hashes, HDF5 dimensions, size, and timestamps
- run `finalize --confirm` only after all physical moves validate
- retain a registry snapshot and reverse journal for `rollback --confirm`

## Local cross-study benchmark intake: BALF_BLOOD_COPD

Objective:

- register the four-library BALF/PBMC study used as the gdTAI cross-study
  benchmark under dataset `BALF_BLOOD_COPD`
- retain the complete raw/interim study workspace at its original path
- store only the selected validation H5AD physically inside the project
- keep the cohort validation-only unless atlas integration is separately
  approved

Phase or task:

- Dataset-centered intake of the BALF/PBMC COPD cross-study benchmark

Exact `.py` scripts:

- `workflows/maintenance/migrate_balf_blood_copd_dataset.py`
- `workflows/maintenance/relocate_balf_blood_copd_raw_workspace.py`

Core outputs:

- `/home/tanlikai/databank/owndata/singlecell/`
- `data/datasets/BALF_BLOOD_COPD/raw/legacy_source`
- `data/datasets/BALF_BLOOD_COPD/processed/current.h5ad`
- `data/registry/migrations/balf_blood_copd_20260716/`
- `data/registry/migrations/balf_blood_copd_h5ad_only_20260716/`
- dataset, library, file, and storage-index registry entries

Standard behavior:

- keep raw, FASTQ, BAM, matrix, TCR, demultiplexing, and other study files in
  the original workspace
- retain the selected validation H5AD as the only physical cohort data file in
  `data/datasets/BALF_BLOOD_COPD/`
- keep the historical H5AD path as a compatibility link to the canonical H5AD
- fingerprint the selected H5AD, four filtered matrices, four raw matrices,
  harmonized TCR table, and two demultiplexing tables
- register `LIB1`-`LIB4` as inactive validation-only libraries
