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
- independent-validation intake evidence:
  `data/registry/migrations/balf_blood_copd_20260716/`
- independent-validation H5AD-only storage evidence:
  `data/registry/migrations/balf_blood_copd_h5ad_only_20260716/`
- reproducibility and rollback guides: `docs/`

## Current review packages

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
- `BALF_BLOOD_COPD` independent validation cohort:
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
