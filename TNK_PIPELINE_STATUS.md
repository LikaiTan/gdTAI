# TNK_PIPELINE_STATUS.md

## Current milestone

- gdTAI V4.2 current-atlas recovery preflight completed with
  `PASS_REVIEW_REQUIRED` on 2026-08-14 after the SSD-resident
  `integrated.h5ad` disappeared before project-data execution; all 20 recovery,
  equivalence, read-only, and fail-closed checks passed

## Next action

- review the V4.2 recovery-preflight PDF and checksum-bound recovery execution
  approval template
- decide whether to authorize recovered development-data sparse staging, A100
  scVI, consensus clustering, and pseudo-label QC without opening the three
  locked evaluation cohorts
- publish and review integration/clustering QC before any classifier fitting
- after an approved integration and clustering QC report, decide separately
  whether to authorize V4.2 classifier fitting and nested comparison
- do not authorize Gate D release fitting, V4.1 promotion, or whole-atlas
  inference from this failed experiment
- do not promote a model, fit a release artifact, or run whole-atlas inference
  from the incomplete CPU evidence or the GPU feasibility probes
- do not promote or claim superiority for another gdTAI model before grouped
  resampling, expression-independent labels, and fold-local selection are in place
- review the frozen-profile negative-control report and official-GEO metadata
  reconciliation package
- decide separately whether to approve additive corrected metadata fields for
  GSE169246 and the other GEO-resolved rows
- keep every extension cohort separate until merge and integration receive
  explicit approval

## Current blockers or review items

- the original `high_speed_temp/Integrated_dataset/integrated.h5ad` is absent,
  and no exact copy was found in the project or databank search; the previously
  recorded project-data execution approval is therefore invalid
- the recovery preflight passed 20/20 checks using the intact canonical
  3,705,384-cell `TNK_cleaned.h5ad`, the exact 78-row Phase-3 exclusion
  manifest, and the original harmonized metadata sources
- the recovered effective input contains exactly 3,705,306 cells and 27,413
  genes, reproduces all five frozen metadata missing/unique-count audits, and
  recovers all 21,054/21,054 primary dual-annotation NK anchors
- no project-data integration, scVI fitting, clustering, pseudo-labeling, or
  classifier fitting occurred; execution remains blocked until the new
  checksum-bound `RECOVERY_EXECUTION_APPROVAL.json` is explicitly activated
- V4.2 implementation QC passed 15/15 checks; the synthetic 1,600-cell scVI
  fit ran on `cuda:0` using the A100 80-GB GPU, produced a finite 8-dimensional
  latent matrix, and completed RAPIDS neighbors plus two Leiden seeds
- 11/11 deterministic tests passed for approval binding, role leakage,
  fail-closed execution, source-balanced sampling, CSR/CSC H5AD slicing,
  batch-key fallback, expression-independent boundary construction, and
  consensus pseudo-label rules
- all nine source H5AD size and nanosecond-mtime pairs remained unchanged; no
  project-data integration, project-data model fit, or classifier fit occurred
- the earlier `PROJECT_DATA_INTEGRATION_APPROVAL.json` no longer authorizes
  execution because the execution config, runner, core, and physical input
  contract changed for recovery
- V4.2 integration preflight passed 59/59 checks on 4,023,462 development cells
  and 439,979 whole-cohort locked-evaluation cells
- all 21,054 primary dual-annotation NK anchors mapped to the current atlas;
  every development and locked cohort has 50/50 Stage-1 genes
- GSE169246 has only 145/197 full classifier genes, so its 7,770 author NK cells
  are a primary Stage-1 passage challenge but final-cascade FPR is explicitly a
  reduced-feature sensitivity analysis
- GSE315928 retains 66,813 paired-alpha-beta controls with 197/197 genes;
  GSE121636/GSE121637 remains a mixed stress cohort rather than author NK truth
- the development intersection has 14,265 genes and 13,975 eligible non-TCR,
  non-mitochondrial, non-ribosomal, non-immunoglobulin genes before selecting
  4,000 HVGs
- the resource audit estimated 223.6 GiB conservative peak RAM, verified 875.1
  GiB free SSD, and found an idle A100 80-GB GPU
- all nine input H5ADs were hash-verified and retained identical size and
  nanosecond modification time; no H5AD was modified
- the checksum-bound implementation approval is active; project-data
  integration, classifier fitting, and every later release action remain
  blocked behind their respective QC gates

- V4.2 Step 0 verified 336,780 V4.1 NK controls: 315,726 (`93.7%`) had only a
  single scVI NK annotation and 21,054 (`6.3%`) had dual scVI/author agreement
- all 15 checksum-bound saved-OOF candidate/fold evaluations passed the frozen
  50% NK source cap under the dual-annotation truth counterfactual; worst-source
  passage among sources with at least 100 cells was `5.27%-19.37%`, versus
  `90.97%-99.21%` using all V4.1 controls
- the counterfactual is diagnostic only because the saved V4.1 candidates were
  trained with weak controls; no V4.2 performance claim is available
- the user proposed integrating new data and repeated clustering to expand NK
  training labels; this is scientifically reasonable only as a development
  pseudo-label lane with locked whole-dataset validation and explicit approval
- the proposed development cohorts contain 318,156 T/NK candidates; the locked
  GSE169246 cohort provides 7,770 author NK and 54,925 productive alpha-beta T
  cells, and GSE315928 provides 66,813 paired-alpha-beta controls
- no V4.2 integration, scVI fit, clustering, pseudo-labeling, classifier fit,
  threshold search, promotion, release fitting, or atlas inference is approved

- V4.1-GPU Gate C was explicitly approved and checksum-bound; 99/99 recorded
  inner-fold fits converged and all 57 complete threshold-frontier files were
  retained
- zero of 15 Stage-1 candidate-by-outer-fold evaluations passed the frozen
  recall plus maximum-source NK-passage rule; the one-cell GSE190870 stratum
  reached 100% passage, but even among NK strata with at least 100 cells the
  best worst-source passage remained 90.97%-93.23%
- zero compact seven-gene, V2-like TCR-gene, or raw legacy-score operating
  points passed balanced or high-purity constraints across the three folds
- V4 Stage 2, extension-model false-positive rates, and paired hierarchical
  bootstrap were not estimable because Stage 1 selected no eligible candidate;
  no V4.1 model was promoted, released, or applied to the atlas

- V4 v1.2 CPU outer folds completed for HRA005041 and GSE144469; no V4,
  compact seven-gene, V2-like, or legacy-score threshold satisfied either
  frozen operating mode in either fold
- all 16 Stage-1 SAGA candidates were nonconverged in both completed CPU outer
  folds, and the selected Stage-1 threshold passed 100% of held-out strict-NK
  cells in both folds
- the BALF_BLOOD_COPD CPU fold was terminated at the user's direction before
  an atomic fold result was written; completed CPU folds and the feature cache
  remain preserved as read-only diagnostic evidence
- GPU feasibility passed on synthetic data with the A100 80 GB: the 0.835-GiB
  gene cache has ample memory headroom, deterministic PyTorch ridge logistic
  produced bit-identical repeat probabilities, and GPU XGBoost produced
  bit-identical repeat probabilities plus portable UBJ export
- the default `/tmp/nvidia-mps` daemon belongs to another user; gdTAI direct
  CUDA succeeds with an isolated `CUDA_MPS_PIPE_DIRECTORY`, and the plan
  prohibits modifying or connecting to the other user's daemon
- cuML logistic regression is excluded from the proposed candidates because a
  full-scale synthetic fit emitted a line-search stop and the installed version
  did not expose adequate convergence provenance
- protocol v1.2 Step 1 passed all frozen input, feature, truth-label,
  exclusion-feasibility, grouped-split, hash, and file-state checks
- Step 2 implementation passes 17 deterministic unit/integration tests,
  including all model families, fair comparators, threshold guardrails,
  hierarchical bootstrap, sparse extraction, and print-safe reporting
- the pinned Step 2 evaluation config, runner, and core SHA-256 values are
  `de930d4484a57a8ba77d4df79da2d2018c1b35e7095da5f6bdb6db8279082b32`,
  `9102b8ebfe617df8c218bbd34a58a5e5ac31873a938de969f74f413913d85a07`, and
  `eaf08999a6cffe2389c3deb73669b0591c6847d7fd6aa372b32c3368919fac7e`
- the first authorized cache attempt verified all source hashes and extracted
  all rows, then stopped before fitting because the cache audit compared
  all-cell exclusions against Step 1 flags that were intentionally populated
  only for positive/sensitivity and sampled-NK rows; the evaluator now requires
  exact reproduction on that frozen scope and applies the unchanged rules to
  every cell in memory
- an incomplete timing run was stopped before any outer-fold result after the
  serial fold implementation proved inefficient and emitted unrecorded SAGA
  convergence warnings; the evaluator now runs the same frozen folds in
  deterministic parallel and records iteration/convergence provenance
- the clean execution uses 26 candidate workers and three fold workers for a
  maximum of 78 concurrent fits, within the 80-CPU contract; this changes only
  scheduling, not data, models, seeds, folds, or selection
- the registered HRA005041 `log1p(CP10K)` matrix passed a full 766,639-cell
  inverse-library-size audit: zero empty rows, zero rows outside the frozen
  `1e-4` tolerance, and maximum absolute deviation
  `7.275957614183426e-12` from 10,000
- GDT_2020AUG_woCOV, MalteGDT, and GDTlung are sensitivity-only and have zero
  training weight; their reduced feature coverage cannot affect fitting
- the fixed CD4/Treg exclusions passed their feasibility guardrail with a
  source-macro primary-gdT recall ceiling of `98.38%`; those thresholds remain
  frozen and were not tuned
- no unused gdT-positive cohort remains after adaptive model iteration;
  `BALF_BLOOD_COPD` and screened extension cohorts are development/stress
  benchmarks, not independent or prospective holdouts
- V2 validation metrics and V3 BALF metrics are model-selection estimates
- the promoted V3 artifact checksum is valid, but its pickle, manifest,
  completed-run summary, and documented training composition do not describe
  one coherent build record
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
  - status: currently unavailable because the mirrored SSD tree disappeared;
    V4.2 has an audited raw-count recovery path from `TNK_cleaned.h5ad`, but it
    is not yet approved for project-data execution

Additional current state:

- gdTAI V4 protocol v1.2 keeps the 197-gene universe, labels, folds, model
  families, grids, guardrails, CD4/Treg thresholds, and promotion rule fixed
- protocol v1.2 removes silver, dual/ambiguous, and all sorted cells from both
  fitting stages; no fitting result informed this amendment
- gdTAI V4 Step 1 created a checksum-pinned 1,137,739-cell label manifest and
  deterministic three-fold grouped split manifest without fitting a model
- protocol v1.2 Step 1 completed with `PASS`, zero failures, and zero warnings;
  its cell-label and grouped-split manifest SHA-256 values are
  `8157cbebfedeb84fc34eba05429aa8a7a834c6f7ceba15bd4790b5dd06bf7e0c` and
  `c84da2ca999676bab0ed180ae3d380e1e7d5b2e08da886ae2f6f912f9e0080a7`
- all 16 declared preflight inputs matched registered hashes where available
  and retained identical size and nanosecond modification time
- GSE144469 joined 107,068/107,068 raw-expression rows uniquely by the SRR and
  barcode encoded in cell IDs; its valid raw matrix is `layers/counts`
- the earlier experimental V4 bundle was archived intact under
  `archive/retired_experiments/gdTAI_v4_experimental_precommit_20260807/`
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
- the `BALF_BLOOD_COPD` cross-study benchmark was subsequently registered; its
  storage was refined so the approximately 949 GB raw/interim
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
- `BALF_BLOOD_COPD` retains the legacy registry role string
  `gdTAI_independent_external_validation`, but it was reused in Round 12 versus
  Round 14 promotion and is scientifically classified as a reused cross-study
  benchmark; it remains inactive for Phase 0, current-milestone integration,
  and extended-atlas integration
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
  validated high-purity fallback; the 2026-08-06 audit found that their
  corrected TCR-only F1 values are effectively tied and that the reused
  benchmark cannot establish Round 14 superiority
- `GSE305372` was excluded and retired on 2026-08-06; it has no active registry
  rows, its local scientific data were deleted, and its code, reports, source
  checksums, and historical results were moved under `archive/`

## Current review artifacts

- gdTAI V4.2 modeling-integration preflight:
  - `configs/models/gdtai/v4_2_cohort_roles.csv`
  - `configs/models/gdtai/v4_2_integration_preflight.json`
  - `gdT_prediction/gdtai_v4_2_integration_preflight/index.html`
  - `gdT_prediction/gdtai_v4_2_integration_preflight/gdtai_v4_2_integration_preflight_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_integration_preflight/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_integration_preflight/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_integration_preflight/`

- gdTAI V4.2 Step 0 NK-label audit and proposed repair:
  - `docs/GDTAI_V4_2_PRECOMMITTED_PLAN.md`
  - `gdT_prediction/gdtai_v4_2_step0/index.html`
  - `gdT_prediction/gdtai_v4_2_step0/gdtai_v4_2_step0_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_step0/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_step0/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_step0/`

- gdTAI V4.1-GPU pre-fit feasibility and protocol:
  - `docs/GDTAI_V4_GPU_PRECOMMITTED_PLAN.md`
  - `gdT_prediction/gdtai_v4_gpu_precommit/index.html`
  - `gdT_prediction/gdtai_v4_gpu_precommit/gdtai_v4_gpu_precommitted_plan.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_gpu_precommit/gpu_feasibility_summary.md`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_gpu_precommit/gpu_feasibility_checks.csv`
- gdTAI V4 Step 1 preflight and split freeze:
  - `gdT_prediction/gdtai_v4_preflight/index.html`
  - `gdT_prediction/gdtai_v4_preflight/gdtai_v4_preflight_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_preflight/gdtai_v4_preflight_summary.md`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_preflight/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_preflight/`
- precommitted gdTAI V4 training and validation protocol:
  - `docs/GDTAI_V4_PRECOMMITTED_PLAN.md`
  - `gdT_prediction/gdtai_v4_precommit/index.html`
  - `gdT_prediction/gdtai_v4_precommit/gdtai_v4_precommitted_plan.pdf`
- independent gdTAI training and evaluation audit:
  - `docs/GDTAI_METHODOLOGY_AUDIT.md`
  - `gdT_prediction/gdtai_methodology_audit/index.html`
  - `gdT_prediction/gdtai_methodology_audit/gdtai_methodology_audit_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_methodology_audit/gdtai_methodology_audit_summary.md`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_methodology_audit/`
- gdTAI external methodological review:
  - `Integrated_dataset/logs/gdT_prediction/external_review/gdtai_external_review_report.md`
  - `gdT_prediction/external_review/gdtai_external_review_report.html`
  - `gdT_prediction/external_review/gdtai_external_review_report.pdf`
- gdTAI-kimi experimental nested cross-dataset evaluation (not promotable):
  - `Integrated_dataset/logs/gdT_prediction/gdtai_kimi/gdtai_kimi_report.md`
  - `gdT_prediction/gdtai_kimi/gdtai_kimi_report.pdf`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_kimi/`
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
- `BALF_BLOOD_COPD` cross-study benchmark intake:
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
