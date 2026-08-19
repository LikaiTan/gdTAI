# OUTPUTS.md

## Canonical control files

- `AGENTS.md`
- `TNK_PIPELINE_RUNBOOK.md`
- `TNK_PHASES_0_4_SCRIPT.md`
- `TNK_PIPELINE_STATUS.md`
- `CHANGELOG.md`
- `DECISIONS.md`
- `OUTPUTS.md`

## Canonical output root

- `Integrated_dataset/`

## Canonical milestone H5AD files

- `Integrated_dataset/TNK_candidates.h5ad`
- `Integrated_dataset/TNK_cleaned.h5ad`
- `Integrated_dataset/integrated.h5ad`
- `Integrated_dataset/integrated_full_atlas.h5ad`
  - rebuilt 5,933,312-cell atlas containing the historical atlas, eight new
    T/NK-filtered cohorts, and `BALF_BLOOD_COPD`

## Approved milestone exception

- `Integrated_dataset/TNK_candidates_supp.h5ad`
  - supplementary 10x 5' intake lane only
- `high_speed_temp/Integrated_dataset/No_TCR_Gene_GSEs.h5ad`
  - standalone carve-out of the approved no-TCR-gene GSE removals

## Active large-H5AD mirror

- mirrored SSD root:
  - `/ssd/tnk_phase3/Integrated_dataset/`
- stable project-root symlink:
  - `high_speed_temp/Integrated_dataset`

## Canonical output locations by type

- H5AD milestones:
  - NFS under `Integrated_dataset/`
  - or mirrored SSD path when the runbook says the large-H5AD exception is active
- tables:
  - `Integrated_dataset/tables/`
- logs:
  - `Integrated_dataset/logs/`
- figures:
  - `Integrated_dataset/figures/`

## Canonical project organization

- active workflows: `workflows/`
- shared Python utilities: `src/tnk_atlas/`
- dataset registries: `configs/datasets/`
- metadata and assay configuration: `configs/`
- physical per-dataset inputs: `data/datasets/<dataset_id>/`
- shared pre-integration artifacts: `data/shared/`
- legacy path trees: `data/compat/`
- lifecycle compatibility views: `data/raw/`, `data/interim/`, and
  `data/processed/`
- storage index: `data/registry/storage_index.csv`
- model release registry: `configs/models/gdtai/model_registry.csv`
- reader-facing report index: `reports/`
- retired code and records: `archive/`
- migration plans and checksum maps: `maintenance/reorganization/`
- physical storage migration evidence:
  `data/registry/migrations/dataset_centered_20260716/`
- BALF_BLOOD_COPD benchmark-intake evidence:
  `data/registry/migrations/balf_blood_copd_20260716/`
- BALF_BLOOD_COPD H5AD-only storage evidence:
  `data/registry/migrations/balf_blood_copd_h5ad_only_20260716/`
- reproducibility and rollback guides: `docs/`

## Current review packages

- gdTAI V4.3 recall-failure investigation:
  - `gdT_prediction/gdtai_v4_3_recall_failure/index.html`
  - `gdT_prediction/gdtai_v4_3_recall_failure/gdtai_v4_3_recall_failure_report.pdf`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_3_recall_failure/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_3_recall_failure/`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_3_recall_failure/`
  - decision: V4.3 has no overall deployment advantage; its BALF ranking is
    strong, but threshold transfer fails and paired-alpha-beta FPR is higher
    than V3 in seven of eight extension datasets

- gdTAI V4 final common-lockbox evaluation:
  - `gdT_prediction/gdtai_v4_3_final_evaluation/index.html`
  - `gdT_prediction/gdtai_v4_3_final_evaluation/gdtai_v4_3_final_evaluation_report.pdf`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_3_common_lockbox/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_3_final_evaluation/`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_3_final_evaluation/`
  - decision: `FAIL_V4_NOT_SUPERIOR`; V4.3 remains diagnostic, the lockbox is
    consumed, and V3 balanced remains the current default

- full T/NK atlas rebuild:
  - `full_atlas_rebuild/index.html`
  - `full_atlas_rebuild/full_atlas_rebuild_report.pdf`
  - `Integrated_dataset/integrated_full_atlas.h5ad`
  - `Integrated_dataset/tables/full_atlas_rebuild/`
  - `Integrated_dataset/figures/full_atlas_rebuild/`
  - `Integrated_dataset/logs/full_atlas_rebuild/`
  - SSD model and intermediates:
    `/ssd/tnk_phase3/Integrated_dataset/full_atlas/`
  - decision: structural integration and Phase 4 score QC passed for 5,933,312
    cells; repaired-TCR propagation and the GSE169246 additive compartment
    correction remain separate gated tasks

- sample-aware productive TCR join rebuild:
  - `gdT_prediction/gdtai_v4_2_tcr_join_rebuild/index.html`
  - `gdT_prediction/gdtai_v4_2_tcr_join_rebuild/tcr_join_rebuild_report.pdf`
  - `Integrated_dataset/tables/tcr_join_rebuild/`
  - `Integrated_dataset/figures/tcr_join_rebuild/`
  - `Integrated_dataset/logs/tcr_join_rebuild/`
  - decision: `PASS_SIDECAR_READY`; 14/14 sources contribute 3,041,871
    checksum-bound replacement rows with per-chain UMI/read provenance, no
    source remains quarantined, and no H5AD propagation occurred

- full-atlas TCR sidecar application:
  - `gdT_prediction/gdtai_v4_2_tcr_sidecar_application/index.html`
  - `gdT_prediction/gdtai_v4_2_tcr_sidecar_application/tcr_sidecar_application_report.pdf`
  - `Integrated_dataset/tables/tcr_sidecar_application/`
  - `Integrated_dataset/figures/tcr_sidecar_application/`
  - `Integrated_dataset/logs/tcr_sidecar_application/`
  - decision: `PASS_TCR_H5AD_CANDIDATE`; 2,155,409 atlas rows were corrected,
    all 17 checks passed, and canonical publication remains pending

- post-repair TCR truth and NK-boundary audit:
  - `gdT_prediction/gdtai_v4_2_post_tcr_audit/index.html`
  - `gdT_prediction/gdtai_v4_2_post_tcr_audit/post_tcr_truth_nk_audit.pdf`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_post_tcr_audit/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_post_tcr_audit/`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_post_tcr_audit/`
  - decision: `PASS_POST_TCR_TRUTH_NK_AUDIT`; 20,922 clean primary NK anchors
    remain and prior unsafe-overlay boundary conflicts are largely resolved

- gdTAI V4.2 T/NK-restricted reintegration and NK-boundary review:
  - `gdT_prediction/gdtai_v4_2_tnk_reintegration/index.html`
  - `gdT_prediction/gdtai_v4_2_tnk_reintegration/gdtai_v4_2_tnk_reintegration_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_tnk_reintegration/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_tnk_reintegration/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_tnk_reintegration/`
  - productive-chain figure:
    `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_tnk_reintegration/nk_boundary_productive_tcr_umap.png`
  - decision: `PASS_AUDIT_NK_CORE_REJECTED`; 78.44% of the boundary and 74.06%
    of the previous clusters 0/3/5 candidate core carry productive TRA or TRB,
    so no NK training core is accepted and no classifier or promotion occurred

- gdTAI V4.2 NK/T-lineage, raw TCR, and unsupervised-cluster review:
  - `gdT_prediction/gdtai_v4_2_nk_cluster_review/index.html`
  - `gdT_prediction/gdtai_v4_2_nk_cluster_review/gdtai_v4_2_nk_cluster_review_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_nk_cluster_review/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_nk_cluster_review/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_nk_cluster_review/`
  - decision: `PASS_VISUAL_REVIEW_EXTENDED`; hold all 461 dual-evidence
    cluster-19 cells as provisional pending source/library ambient, doublet, and
    paired-TCR review; hold 8 calls outside the core and every rescue tier

- gdTAI V4.2 conservative NK-definition repair:
  - `gdT_prediction/gdtai_v4_2_nk_definition_repair/index.html`
  - `gdT_prediction/gdtai_v4_2_nk_definition_repair/gdtai_v4_2_nk_definition_repair_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_nk_definition_repair/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_nk_definition_repair/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_nk_definition_repair/`
  - historical decision: `PASS_NK_REFERENCE_READY`; the later raw TCR review
    holds all 469 rows from training pending source/library review, and no
    classifier was fitted

- gdTAI V4.2 sidecar visualization and confounding diagnostics:
  - `gdT_prediction/gdtai_v4_2_nk_reference_diagnostics/index.html`
  - `gdT_prediction/gdtai_v4_2_nk_reference_diagnostics/gdtai_v4_2_nk_reference_diagnostics_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_nk_reference_diagnostics/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_nk_reference_diagnostics/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_nk_reference_diagnostics/`
  - decision: `PASS_DIAGNOSTIC_ONLY`; the current sidecar cannot confidently
    propagate NK identity to new-cohort candidates, and no classifier was
    fitted

- gdTAI V4.2 integration, clustering, and pseudo-NK consensus QC:
  - `gdT_prediction/gdtai_v4_2_nk_reference/index.html`
  - `gdT_prediction/gdtai_v4_2_nk_reference/gdtai_v4_2_nk_reference_qc_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_nk_reference/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_nk_reference/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_nk_reference/`
  - decision: technical execution passed; scientific QC failed with zero
    selected pseudo-NK cells, and no classifier was fitted

- gdTAI V4.2 cluster-stage resource amendment:
  - `configs/models/gdtai/v4_2_integration_execution.json`
  - `gdT_prediction/gdtai_v4_2_cluster_resource_preflight/index.html`
  - `gdT_prediction/gdtai_v4_2_cluster_resource_preflight/gdtai_v4_2_cluster_resource_preflight_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_cluster_resource_preflight/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_cluster_resource_preflight/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_cluster_resource_preflight/`
  - decision: 18/18 checks passed; the user approved the stage-specific SSD
    floor and clustering/consensus completed on 2026-08-16

- gdTAI V4.2 current-atlas recovery preflight:
  - `configs/models/gdtai/v4_2_integration_execution.json`
  - `gdT_prediction/gdtai_v4_2_recovery_preflight/index.html`
  - `gdT_prediction/gdtai_v4_2_recovery_preflight/gdtai_v4_2_recovery_preflight_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_recovery_preflight/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_recovery_preflight/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_recovery_preflight/`
  - decision: 20/20 checks passed; recovery execution and classifier fitting
    remain blocked by separate explicit supervision gates

- gdTAI V4.2 sidecar-integration implementation QC:
  - `configs/models/gdtai/v4_2_integration_execution.json`
  - `gdT_prediction/gdtai_v4_2_implementation_qc/index.html`
  - `gdT_prediction/gdtai_v4_2_implementation_qc/gdtai_v4_2_implementation_qc_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_implementation_qc/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_implementation_qc/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_implementation_qc/`
  - decision: 15/15 checks passed; project-data sidecar execution and later
    classifier fitting remain blocked by separate supervision gates

- gdTAI V4.2 modeling-integration preflight:
  - `configs/models/gdtai/v4_2_cohort_roles.csv`
  - `configs/models/gdtai/v4_2_integration_preflight.json`
  - `gdT_prediction/gdtai_v4_2_integration_preflight/index.html`
  - `gdT_prediction/gdtai_v4_2_integration_preflight/gdtai_v4_2_integration_preflight_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_integration_preflight/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_integration_preflight/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_integration_preflight/`
  - decision: 59/59 preflight checks passed; integration and fitting remain
    blocked by separate explicit supervision gates

- gdTAI V4.2 Step 0 NK-label and saved-OOF counterfactual audit:
  - `docs/GDTAI_V4_2_PRECOMMITTED_PLAN.md`
  - `gdT_prediction/gdtai_v4_2_step0/index.html`
  - `gdT_prediction/gdtai_v4_2_step0/gdtai_v4_2_step0_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_step0/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_step0/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_step0/`
  - decision: label audit passed; integration, clustering, and V4.2 fitting
    remain blocked for explicit supervision approval

- gdTAI V4.1-GPU Gate C negative nested evaluation:
  - `gdT_prediction/gdtai_v4_gpu_nested_evaluation/index.html`
  - `gdT_prediction/gdtai_v4_gpu_nested_evaluation/gdtai_v4_gpu_gate_c_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_gpu_nested_evaluation/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_gpu_nested_evaluation/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_gpu_nested_evaluation/`
  - decision: no eligible V4.1 model; V3 Round 14 remains promoted

- gdTAI V4 Step 1 preflight and grouped split freeze:
  - `gdT_prediction/gdtai_v4_preflight/index.html`
  - `gdT_prediction/gdtai_v4_preflight/gdtai_v4_preflight_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_preflight/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_preflight/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_preflight/`
- independent gdTAI training and evaluation audit:
  - `docs/GDTAI_METHODOLOGY_AUDIT.md`
  - `gdT_prediction/gdtai_methodology_audit/index.html`
  - `gdT_prediction/gdtai_methodology_audit/gdtai_methodology_audit_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_methodology_audit/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_methodology_audit/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_methodology_audit/`
- gdTAI external methodological review:
  - `Integrated_dataset/logs/gdT_prediction/external_review/gdtai_external_review_report.md`
  - `gdT_prediction/external_review/gdtai_external_review_report.html`
  - `gdT_prediction/external_review/gdtai_external_review_report.pdf`
- gdTAI-kimi experimental nested cross-dataset evaluation (not promotable):
  - `Integrated_dataset/logs/gdT_prediction/gdtai_kimi/gdtai_kimi_report.md`
  - `gdT_prediction/gdtai_kimi/gdtai_kimi_report.html`
  - `gdT_prediction/gdtai_kimi/gdtai_kimi_report.pdf`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_kimi/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_kimi/`
  - `Integrated_dataset/models/gdT_prediction_classifier/gdtai_kimi/` (experimental, unregistered)
- gdTAI V4 underperformance causal analysis:
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_gap_analysis/gdtai_v4_gap_analysis.md`
  - `gdT_prediction/gdtai_v4_gap_analysis/index.html`
  - `gdT_prediction/gdtai_v4_gap_analysis/gdtai_v4_gap_analysis_report.pdf`
- extension official-GEO metadata reconciliation:
  - `docs/EXTENSION_GEO_METADATA_RECONCILIATION.md`
  - `configs/metadata/extension_geo_metadata_reconciliation.csv`
  - `Integrated_dataset/logs/geo_metadata_reconciliation/`
  - `Integrated_dataset/tables/geo_metadata_reconciliation/`
- extension frozen gdTAI profile negative-control evaluation:
  - `gdT_prediction/gdtai_extension_evaluation/index.html`
  - `gdT_prediction/gdtai_extension_evaluation/gdtai_extension_evaluation_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_extension_evaluation/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_extension_evaluation/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_extension_evaluation/`
- extension-cohort standalone Phase 1 T/NK filter:
  - `data/interim/extension_intake/tnk_filtered_h5ads_manifest.csv`
  - `data/interim/extension_intake/tnk_filtered_h5ads/`
  - `Integrated_dataset/logs/extension_intake/extension_tnk_filter_qc.md`
  - `Integrated_dataset/logs/extension_intake/extension_tnk_filter_qc.json`
  - `Integrated_dataset/tables/extension_intake/tnk_filter/`
  - `Integrated_dataset/figures/extension_intake/tnk_filter/extension_tnk_filter_overview.png`
- GSE169246 selected processed input:
  - `data/datasets/GSE169246/processed/current.h5ad`
  - `docs/GSE169246_PROCESSED_H5AD_INTAKE.md`
  - `data/registry/snapshots/post_gse169246_processed_h5ad_intake_20260806/`
- extension-cohort standalone Phase 0 review:
  - `Integrated_dataset/logs/extension_intake/extension_phase0_supervision_report.md`
  - `Integrated_dataset/logs/extension_intake/extension_phase0_qc_summary.md`
  - `Integrated_dataset/tables/extension_intake/`
  - `Integrated_dataset/figures/extension_intake/`
- extension tumor-project specimen-context review:
  - `Integrated_dataset/logs/sample_source_refinement/extension_sample_source_review.md`
  - `Integrated_dataset/tables/sample_source_refinement/extensions/`
- Phase 4 QC:
  - `Integrated_dataset/logs/phase4_qc_summary.md`
- tissue correction:
  - `Integrated_dataset/tables/tissue_correction/`
- repaired TCR validation:
  - `Integrated_dataset/tables/tcr_rebuild_phase4/`
  - `Integrated_dataset/logs/tcr_rebuild_phase4_qc.md`
- pre-merge TCR audit:
  - `Integrated_dataset/tables/premerge_tcr_audit/`
  - `Integrated_dataset/logs/premerge_tcr_audit_summary.md`
- HRA005041 intake:
  - `downloads/per_gse_h5ad_with_metadata/HRA005041_T_cells_subset.h5ad`
  - `Integrated_dataset/tables/HRA005041_tcr_intake/`
  - `Integrated_dataset/logs/HRA005041_tcr_intake.md`
- HRA005041 standalone Phase 4 review:
  - `Integrated_dataset/tables/HRA005041_phase4/`
  - `Integrated_dataset/figures/HRA005041_phase4/`
  - `Integrated_dataset/logs/HRA005041_phase4.log`
- Sorted_gdT intake review:
  - `newdata/Sorted_gdT/GDT_2020AUG_woCOV_sorted_gdt.h5ad`
  - `newdata/Sorted_gdT/GDTlung2023july_7p_sorted_gdt.h5ad`
  - `newdata/Sorted_gdT/MalteGDT_sorted_gdt.h5ad`
- gdTAI v3 Round 12 versus Round 14 revalidation:
  - `gdT_prediction/gdtai_v3_round12_vs_round14/index.html`
  - `gdT_prediction/gdtai_v3_round12_vs_round14/gdtai_v3_round12_vs_round14_report.pdf`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v3_round12_vs_round14/`
  - `Integrated_dataset/models/gdT_prediction_classifier/gdtai_v3_round12_vs_round14/`
  - `Integrated_dataset/tables/Sorted_gdT/`
  - `Integrated_dataset/logs/Sorted_gdT/`
  - `Integrated_dataset/figures/GDT_2020AUG_woCOV_phase4/`
  - `Integrated_dataset/figures/GDTlung2023july_7p_phase4/`
  - `Integrated_dataset/figures/MalteGDT_phase4/`
- `BALF_BLOOD_COPD` reused cross-study benchmark:
  - `data/datasets/BALF_BLOOD_COPD/processed/current.h5ad`
  - `data/datasets/BALF_BLOOD_COPD/processed/artifacts/phase4_final_annotated.h5ad`
  - `docs/BALF_BLOOD_COPD_INTAKE.md`
  - `data/registry/migrations/balf_blood_copd_20260716/validation.json`
  - `data/registry/migrations/balf_blood_copd_h5ad_only_20260716/validation.json`
- no-TCR-gene GSE carve-out:
  - `Integrated_dataset/tables/no_tcr_gene_gse_removal_counts.csv`
  - `Integrated_dataset/logs/no_tcr_gene_gse_removal.md`
- disease-status review export:
  - `Integrated_dataset/tables/disease_status_correction/`
- cluster annotation evidence review:
  - `Integrated_dataset/tables/annotation_evidence_review/`
  - `Integrated_dataset/logs/annotation_evidence_review/`
  - `Integrated_dataset/figures/annotation_evidence_review/`

## Retired External Applications

- GSE305372 gdTAI application:
  - retirement record and archived code:
    `archive/retired_experiments/GSE305372_external_application/`
  - archived rendered reports:
    `archive/historical_reports/GSE305372_external_application/`
  - the original local dataset was deleted and is recoverable from GEO using
    the archived source manifest
