# TNK_PIPELINE_STATUS.md

## Current milestone

- Post-Phase 4 extension-cohort frozen-model and official-GEO metadata reviews
  complete; additive metadata correction remains gated

## Next action

- review the frozen-profile negative-control report and official-GEO metadata
  reconciliation package
- decide separately whether to approve additive corrected metadata fields for
  GSE169246 and the other GEO-resolved rows
- keep every extension cohort separate until merge and integration receive
  explicit approval

## Current blockers or review items

- extension-cohort merge and integration are not approved
- GSE169246 has a confirmed compartment-join error: all 239,418 retained `_b`
  cells are official blood libraries, while current local tissue/specimen fields
  include tumor sites or non-blood contexts; correction is not yet approved
- frozen-model screening is complete, but these cohorts contain zero gdT-gold
  and zero gdT-silver cells under the project TCR rules; they cannot support
  new-cohort recall, F1, ROC-AUC, PR-AUC, or model promotion
- the extension inputs used alpha-beta TCR assays rather than gdTCR assays;
  zero productive TRG/TRD is therefore expected, and these cohorts support
  alpha-beta T-cell and NK false-positive evaluation rather than independent
  TCR-defined gdT recall
- five extension cohorts have unmatched productive TCR-contig warnings that
  require review; all seven nevertheless passed structural QC
- Phase 4 scoring outputs remain under user review
- disease-status correction remains export-only and has not been written back
  into `integrated.h5ad`
- migration of the mirrored SSD-side `integrated.h5ad` back to NFS is not yet
  approved
- cluster annotation evidence review outputs are review-only and have not been
  written back into milestone H5AD annotation columns

## Active exceptions in force

- active RAM ceiling override: `800 GB`
- large-H5AD mirrored SSD workflow is active
- no active QC-gate exception

## Current canonical objects

- canonical candidate milestone:
  - `Integrated_dataset/TNK_candidates.h5ad`
- canonical cleaned milestone:
  - `Integrated_dataset/TNK_cleaned.h5ad`
- canonical current integrated milestone for downstream analysis:
  - `high_speed_temp/Integrated_dataset/integrated.h5ad`

Additional current state:

- seven extension H5ADs containing 498,901 cells and 154 libraries were built
  under `data/interim/extension_intake/built_h5ads/`
- all seven extension H5ADs passed sparse raw-count, metadata, TCR-schema,
  logical-flag, and sample-plus-barcode Phase 0 structural QC
- extension tumor-project specimen context initially passed the rule-based
  structural gate for all 498,901 cells with zero blank/format violations;
  later official-GEO review exposed a semantic compartment error in GSE169246,
  and original tissue and diagnosis fields remain unmodified
- all eight new inputs were independently filtered to `758,135` T/NK
  candidates from `893,126` cells; all `338,046` productive-TCR cells were
  retained
- all eight filtered outputs passed required-metadata completeness, exact
  source-accession, cohort-ID, unique cell/library-barcode, canonical TCR-flag,
  and specimen-context checks
- GSE169246 metadata now distinguish 51 biological samples from 78 source
  libraries using the full barcode suffix; 76 libraries retain T/NK cells
- source H5AD checksums were unchanged and no source H5AD was modified
- all four frozen gdTAI profiles were applied read-only to all 758,135
  extension T/NK cells using raw-count CP10K plus log1p normalization
- V3 Round 12 high-purity had the lowest pooled alpha-beta TCR FPR (`0.227%`),
  while V2 high-purity had the lowest strict-NK FPR (`0.725%`)
- after excluding reduced-feature GSE169246 as a sensitivity analysis, V3
  Round 12 high-purity had the lowest known-negative union FPR (`0.295%`)
- no model was selected or promoted; historical positive sensitivity still
  favors promoted V3 Round 14 over Round 12
- all eight evaluation-source H5AD SHA-256 values remained unchanged
- official-GEO reconciliation is complete for 58 grouped assertions across 8
  cohorts and 10 accessions: 16 locally verified, 32 GEO-resolved but locally
  unresolved, 6 partially ambiguous, and 4 unavailable in GEO
- the GEO audit passed 10/10 accession-coverage checks and 8/8 read-only H5AD
  file-state checks; it did not write metadata or mutate an H5AD
- GSE169246 `_b` correction must use additive derived columns and retain the
  original fields and full `_b`/`_t` compartment suffix
- Tan et al. 2021 remains excluded as a duplicate of
  `GDT_2020AUG_woCOV`
- the project reorganization is validated: active code is grouped under
  `workflows/`, configuration under `configs/`, and physical pre-integration
  data are dataset-centered under `data/datasets/<dataset_id>/`
- the base dataset-centered storage migration completed with 191
  same-filesystem rename operations
- the independent `BALF_BLOOD_COPD` validation cohort was subsequently
  registered; its storage was refined so the approximately 949 GB raw/interim
  workspace remains at its original location and only the 2,110,825,599-byte
  validation H5AD is physically stored in the project
- the current storage view has 66 dataset directories, 45 promoted
  `current.h5ad` links, and 264 lifecycle compatibility links
- the historical `downloads`, `analysis_26GSE_V4`, and `newdata` paths remain
  usable as compatibility aliases to the same inodes
- `/home/tanlikai/databank/owndata/singlecell` is again the physical
  `BALF_BLOOD_COPD` raw/interim workspace; the project raw path is a
  compatibility link to it
- the historical external H5AD path is a compatibility link to
  `data/datasets/BALF_BLOOD_COPD/processed/artifacts/phase4_final_annotated.h5ad`
- the dataset registry currently contains 67 datasets and 1,790 library rows;
  all 33 Phase 0 active datasets pass strict path and library validation
- `BALF_BLOOD_COPD` is registered as
  `gdTAI_independent_external_validation` and remains inactive for Phase 0,
  current-milestone integration, and extended-atlas integration
- `obs["tissue_corrected"]` is already written into the current integrated
  milestone
- `HRA005041` intake H5AD was rebuilt with harmonized metadata plus productive
  sample-aware `TRA/TRB/TRG/TRD` fields under
  `downloads/per_gse_h5ad_with_metadata/HRA005041_T_cells_subset.h5ad`
- `HRA005041` now also carries standalone Phase 4 score columns from the
  project-standard scoring method
- three sorted-gdT Seurat RDS datasets were converted into harmonized standalone
  H5ADs under `newdata/Sorted_gdT/`
- those sorted-gdT intake H5ADs now carry canonical `TRG/TRD` metadata,
  `Sorted_gdT = True`, recomputed Scanpy UMAPs, and standalone Phase 4 score
  columns
- repaired TCR propagation is complete for ten approved GSEs
- the fourteen structurally unscoreable no-TCR-gene GSEs have been removed from
  the current milestone H5ADs
- removed cells were preserved in
  `high_speed_temp/Integrated_dataset/No_TCR_Gene_GSEs.h5ad`
- the current integrated milestone now contains `3,705,306` cells after that
  carve-out
- scANVI fields remain reference-only; simple scVI-based annotation is the
  canonical downstream interpretation layer
- gdTAI v3 Round 12 and Round 14 were revalidated on identical cohorts with
  checksum-pinned artifacts; Round 14 at threshold `0.936` is the promoted
  balanced default, while Round 12 at threshold `0.5` is preserved as the
  validated high-purity fallback
- `GSE305372` was excluded and retired on 2026-08-06; it has no active registry
  rows, its local scientific data were deleted, and its code, reports, source
  checksums, and historical results were moved under `archive/`

## Current review artifacts

- extension official-GEO metadata reconciliation:
  - `docs/EXTENSION_GEO_METADATA_RECONCILIATION.md`
  - `configs/metadata/extension_geo_metadata_reconciliation.csv`
  - `Integrated_dataset/logs/geo_metadata_reconciliation/extension_geo_metadata_reconciliation_audit.md`
  - `Integrated_dataset/tables/geo_metadata_reconciliation/`
- extension frozen gdTAI profile negative-control evaluation:
  - `gdT_prediction/gdtai_extension_evaluation/index.html`
  - `gdT_prediction/gdtai_extension_evaluation/gdtai_extension_evaluation_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_extension_evaluation/gdtai_extension_evaluation_summary.md`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_extension_evaluation/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_extension_evaluation/`
- extension standalone Phase 1 T/NK filter review:
  - `Integrated_dataset/logs/extension_intake/extension_tnk_filter_qc.md`
  - `Integrated_dataset/tables/extension_intake/tnk_filter/extension_tnk_filter_cohort_summary.csv`
  - `Integrated_dataset/tables/extension_intake/tnk_filter/extension_tnk_filter_metadata_audit.csv`
  - `Integrated_dataset/figures/extension_intake/tnk_filter/extension_tnk_filter_overview.png`
  - `data/interim/extension_intake/tnk_filtered_h5ads_manifest.csv`
- extension standalone Phase 0 review:
  - `Integrated_dataset/logs/extension_intake/extension_phase0_supervision_report.md`
  - `Integrated_dataset/logs/extension_intake/extension_phase0_qc_summary.md`
  - `Integrated_dataset/tables/extension_intake/extension_phase0_cohort_summary.csv`
  - `Integrated_dataset/tables/extension_intake/extension_tcr_source_chain_audit.csv`
  - `Integrated_dataset/figures/extension_intake/extension_phase0_cohort_overview.png`
- GSE169246 selected processed input:
  - `docs/GSE169246_PROCESSED_H5AD_INTAKE.md`
  - `data/datasets/GSE169246/processed/current.h5ad`
  - `data/registry/snapshots/post_gse169246_processed_h5ad_intake_20260806/`
- extension tumor-project specimen-context review:
  - `Integrated_dataset/logs/sample_source_refinement/extension_sample_source_review.md`
  - `Integrated_dataset/tables/sample_source_refinement/extensions/extension_tumor_project_context_counts.csv`
- project reorganization:
  - `docs/REORGANIZATION_RECORD.md`
  - `Integrated_dataset/logs/project_reorganization/checkpoints/pre_reorganization_20260715/`
- dataset-centered storage migration:
  - `docs/DATASET_STORAGE_MIGRATION.md`
  - `data/registry/migrations/dataset_centered_20260716/validation.json`
  - `data/registry/snapshots/post_dataset_centered_migration_20260716/`
  - `Integrated_dataset/logs/project_reorganization/checkpoints/pre_dataset_centered_migration_20260716/`
- `BALF_BLOOD_COPD` independent-validation intake:
  - `docs/BALF_BLOOD_COPD_INTAKE.md`
  - `data/registry/migrations/balf_blood_copd_20260716/validation.json`
  - `data/registry/snapshots/post_balf_blood_copd_intake_20260716/`
  - `data/registry/migrations/balf_blood_copd_h5ad_only_20260716/validation.json`
  - `data/registry/snapshots/post_balf_blood_copd_h5ad_only_20260716/`
- gdTAI v3 Round 12 versus Round 14 model decision:
  - `gdT_prediction/gdtai_v3_round12_vs_round14/index.html`
  - `gdT_prediction/gdtai_v3_round12_vs_round14/gdtai_v3_round12_vs_round14_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v3_round12_vs_round14/gdtai_v3_round12_vs_round14_summary.json`
- Phase 4 QC:
  - `Integrated_dataset/logs/phase4_qc_summary.md`
- tissue correction:
  - `Integrated_dataset/tables/tissue_correction/tissue_correction_apply.md`
- repaired TCR review:
  - `Integrated_dataset/logs/tcr_rebuild_phase4_qc.md`
- HRA005041 intake review:
  - `Integrated_dataset/logs/HRA005041_tcr_intake.md`
- HRA005041 standalone Phase 4 review:
  - `Integrated_dataset/logs/HRA005041_phase4.log`
- Sorted_gdT intake review:
  - `Integrated_dataset/logs/Sorted_gdT/sorted_gdt_dataset_summary.md`
- no-TCR-gene GSE carve-out:
  - `Integrated_dataset/logs/no_tcr_gene_gse_removal.md`
- disease-status export review:
  - `Integrated_dataset/tables/disease_status_correction/disease_status_corrected_column_export.csv.gz`
- cluster annotation evidence review:
  - `Integrated_dataset/logs/annotation_evidence_review/current_integrated_annotation_evidence_summary.md`
  - `Integrated_dataset/logs/annotation_evidence_review/plus6_annotation_evidence_summary.md`
- GSE305372 retirement record:
  - `archive/retired_experiments/GSE305372_external_application/README.md`
  - `archive/retired_experiments/GSE305372_external_application/deleted_source_files.csv`
  - `maintenance/reorganization/gse305372_retirement_map.csv`
