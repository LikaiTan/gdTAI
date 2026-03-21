# TNK_PIPELINE_STATUS.md

## Project
Large-scale T/NK integration and γδT-focused scoring workflow across public datasets

---

## Current milestone
- Phase 4 scoring completed on the mirrored SSD-side integrated milestone; final Phase 4 QC package is ready for review

## Current objective
- Keep the validated integrated milestone on SSD without rollback
- Use scVI latent space, Leiden clusters, UMAP, and simple scVI-based annotation as the canonical downstream interpretation
- Keep scANVI label fields and PNG/QC outputs for reference only; do not use them as the primary downstream annotation layer unless a later user decision changes that
- Review the completed Phase 4 continuous score outputs and decide what downstream biological interpretation or migration step should happen next
- Keep the validated Phase 3 H5AD in the mirrored SSD tree for now; do not migrate it back to NFS until the user explicitly requests migration
- Keep tables, PNG figures, logs, scripts, and model artifacts on NFS; use SSD only for large H5AD inputs/outputs
- Keep all 10x 3' inputs excluded, especially the isolated `GSE234069` 3' lane
- Keep supplementary harmonized metadata separate at `analysis_26GSE_V4/outputs/harmonized_metadata_supp.csv`
- Record the Phase 4 QC conclusion and any post-Phase 4 decisions clearly in markdown and git history

---

## Completed milestones
- Initial control-file design completed
- Canonical Phase 0-3 workflow file created
- Phase 0 dataset audit completed across all 33 registry datasets
- `GSE228597` rescued from `adata.raw` and reclassified into Category A
- Backed up all 20 Phase 0 Category C H5ADs before external raw-source rescue review
- Completed external raw-source audit for the 20 Category C H5ADs
- Phase 1 coarse extraction completed for the 13 current Category A datasets
- Phase 0 selected-input rescue completed for all 20 of the 20 original Category C datasets
- `h5ad_v2.csv` restored to all 33 datasets after resolving the last two raw-count rescue cases
- Phase 0 rerun completed on `h5ad_v2.csv`
- Phase 1 coarse extraction rebuilt successfully for the repaired 29-dataset registry
- User approved advancement into Phase 1b conservative cleanup
- Phase 1b conservative cleanup completed on the current merged candidate milestone
- User approved advancement into Phase 1c metadata backup/replacement
- Phase 1c merged metadata backup/replacement completed successfully
- Supplementary 10x 5' candidate discovery completed for six non-registry GSEs
- `GSE234069` supplementary files were separated into `10x_3` and `10x_5` lanes to prevent accidental 3' intake
- Built all six supplementary per-GSE H5AD files under `downloads/per_gse_h5ad_with_metadata/`
- Completed supplementary Phase 0 audit on the six approved 10x 5' datasets
- Fixed duplicate `obs_names` in supplementary `GSE240865` by rewriting cell IDs as `library_id:barcode`
- Completed supplementary Phase 1 coarse extraction and wrote `Integrated_dataset/TNK_candidates_supp.h5ad`
- Wrote supplementary harmonized candidate metadata to `analysis_26GSE_V4/outputs/harmonized_metadata_supp.csv`
- Merged the 897,621 supplementary Phase 1 candidate cells into the main `Integrated_dataset/TNK_candidates.h5ad`
- Verified supplementary candidate metadata carries `technology_simple = 10x 5'` for all supplementary candidate rows
- Phase 2 merged cleanup completed successfully and wrote `Integrated_dataset/TNK_cleaned.h5ad`
- Phase 2 QC passed internal review and Phase 3 execution was started under the approved automatic continuation rule
- Phase 3 scVI integration, RAPIDS embedding, and scANVI coarse T/NK annotation completed successfully on the mirrored SSD-side integrated milestone
- User reviewed the Phase 3 QC package and decided that scANVI labels are too messy for primary downstream use; keep the current integrated milestone without rollback and use simple scVI-based annotation downstream instead
- User explicitly approved Phase 4 execution
- Phase 4 completed successfully with exact package-faithful TRA/TRB/TRD module scoring on the mirrored SSD-side integrated milestone and wrote continuous score columns back into `integrated.h5ad`

---

## Decisions made
- Main thread model: `gpt-5.4`
- Subagent model: `gpt-5.4-mini`
- Phase 5 training/validation is removed for now
- Phase 1–3 script must live in one canonical markdown file
- All outputs must be saved under `Integrated_dataset/`
- Only three milestone H5AD files should be kept:
  - `TNK_candidates.h5ad`
  - `TNK_cleaned.h5ad`
  - `integrated.h5ad`
- High-quality PNG figures are required
- Use conda env `Scanpy_gdTmodel`
- Prefer conda for dependency installation
- Git should be used frequently, with regular pushes to GitHub
- Use git history instead of filename version sprawl
- Phase 4 will use `gdt_tcr_module_sharing_package_full`
- `integrated.h5ad` should be updated in place during Phase 4 when appropriate
- `h5ad.csv` is the authoritative dataset registry for the current run and contains 33 datasets
- `TNK_PHASES_1_3_SCRIPT.md` is now present as the canonical workflow file
- Practical execution environment resolved to `/home/tanlikai/miniconda3/envs/rapids_sc_py310`
- GPU visibility confirmed with `nvidia-smi` on an NVIDIA A100 80GB PCIe system
- `GSE228597` is now raw-count eligible after replacing normalized `X` with integer-like counts restored from `adata.raw` on the current 17,050-gene feature space
- The approved Phase 1 input set for the current run is the 13 Category A datasets
- `TNK_candidates.h5ad` is now the merged Phase 1 milestone object for the approved Category A inputs
- Phase 1b uses a conservative cell-level removal rule plus a user-requested gene filter of `>=500` expressing cells
- Phase 1c replacement must join by string-typed `GSE + barcode` and must emit `project name`, `sampleid`, and `barcodes`
- The user explicitly directed this supplementary run to use the `rapids_sc_py310` conda environment
- Phase 3 will no longer stop at scVI alone; it must include post-scVI scANVI annotation for coarse T and NK labeling
- Candidate Phase 3 reference model path for scANVI:
  - `/home/tanlikai/databank/owndata/fasting/raw/report_from_niuxian/models/census_scanvi_ref_v1`
- The companion reference H5AD exposes `cell_type` as the label field and `batch` as the batch field
- The Phase 3 scANVI reference must be treated as a coarse annotation prior, not as final tissue-state or γδT-subtype truth
- On 2026-03-20, the user explicitly approved Phase 2 execution on the current merged `TNK_candidates.h5ad`, allowed self-QC for Phase 2, and approved direct continuation into Phase 3 without another user checkpoint if Phase 2 QC passes
- On 2026-03-20, the user explicitly required figure generation to remain PNG-only for this workflow
- On 2026-03-20, the `rapids_sc_py310` CUDA import issue was traced to import order; `torch` must be imported before `rapids_singlecell` in this env
- On 2026-03-20, the saved scANVI reference model was confirmed to use the `17,129` genes stored in `model.pt`, not the full `61,888` genes in the companion H5AD
- On 2026-03-20, the user raised the working RAM ceiling for this run from `400 GB` to `800 GB`
- On 2026-03-20, the first Phase 3 recovery run saved the scVI checkpoint but stopped before final writeout; the resumed run must log to file and keep lower-memory guardrails around RAPIDS and scANVI stages
- On 2026-03-21, the Phase 3 implementation was revised to use scANVI on a stratified subset plus latent-space centroid transfer because full-query scANVI adaptation did not scale reliably on the merged milestone
- On 2026-03-21, the user required mitochondrial, ribosomal, and noncoding RNA genes to be excluded from the HVG set used for clustering and UMAP
- On 2026-03-21, the user explicitly paused scANVI and directed the workflow to finish the Leiden clustering / UMAP milestone first
- On 2026-03-21, the Phase 3 heavy H5AD workload was redirected to a mirrored local temp tree rooted at `/ssd/tnk_phase3/Integrated_dataset/` so the final migration back to NFS remains path-stable
- On 2026-03-21, the mirrored SSD tree was also exposed at `high_speed_temp/Integrated_dataset` from the NFS working directory for easier inspection and migration
- On 2026-03-21, the user clarified that validated outputs should remain in the mirrored SSD tree until an explicit migration request is given; automatic move-back to NFS is disabled
- On 2026-03-21, the user further clarified that only large H5AD files should live on SSD; tables, PNG figures, logs, scripts, and model artifacts must remain on NFS
- On 2026-03-21, the resumed scANVI-only pass completed successfully after fixing centroid-transfer index alignment and H5AD serialization of scANVI label fields
- On 2026-03-21, the user decided that the scANVI output is too messy for primary downstream interpretation; keep the scANVI fields in `integrated.h5ad` for reference only, do not roll back Phase 3, and use simple scVI-based annotation for downstream work
- On 2026-03-22, the user explicitly approved Phase 4 implementation
- On 2026-03-22, Phase 4 was implemented as exact package-faithful TRA/TRB/TRD module scoring with temporary normalize+log1p on count-space `X`, continuous score outputs only, and in-place metadata updates to the mirrored SSD-side `integrated.h5ad`
- On 2026-03-22, Phase 4 kept scANVI outputs untouched as reference-only and did not use them to define Phase 4 calls

---

## Pending tasks
- [x] Create `Integrated_dataset/`
- [x] Create `Integrated_dataset/figures/`
- [x] Create `Integrated_dataset/tables/`
- [x] Create `Integrated_dataset/logs/`
- [x] Prepare `TNK_PHASES_1_3_SCRIPT.md`
- [x] Audit all 33 registry datasets
- [x] Classify datasets into Category A / B / C
- [x] Confirm a working single-cell execution environment and GPU visibility
- [x] Rescue `GSE228597` from `adata.raw`
- [x] Define the approved Phase 1 dataset list
- [x] Run Phase 1 coarse T/NK extraction on the approved Category A datasets
- [x] Review Phase 1 QC outputs with the user
- [x] Decide whether to proceed to Phase 1b conservative cleanup
- [x] Define or adjust Phase 1b cleanup rules if needed
- [x] Run Phase 1b conservative cleanup on `TNK_candidates.h5ad`
- [x] Apply the `<500 cells` gene filter during Phase 1b
- [x] Review Phase 1b QC outputs with the user
- [x] Complete Phase 1c merged metadata backup/replacement after Phase 1b approval
- [ ] Review Phase 1c metadata replacement outputs with the user
- [x] Build the six supplementary 10x 5' per-GSE H5AD files
- [x] Run supplementary Phase 0 audit on the six supplementary GSEs
- [ ] Review supplementary Phase 0 QC outputs with the user
- [x] Run supplementary Phase 1 extraction and generate `TNK_candidates_supp.h5ad`
- [ ] Review supplementary Phase 1 QC outputs with the user
- [x] Merge `TNK_candidates_supp.h5ad` into `TNK_candidates.h5ad` after user approval and merge QC
- [x] Start Phase 2 after explicit user approval of the supplementary merge
- [x] Run Phase 2 merged cleanup on the current `TNK_candidates.h5ad`
- [x] Self-review the Phase 2 QC package and decide whether the result is clean enough to enter Phase 3
- [x] Complete the resumed Phase 3 scVI integration and post-scVI scANVI annotation run, then write `integrated.h5ad` and the full PNG QC package
- [x] Review the completed Phase 3 QC package with the user and record whether scANVI is accepted for downstream use
- [ ] Review the 20 Category C raw-source audit with the user and approve dataset-by-dataset rescue scope
- [x] Define Phase 4 scoring workflow
- [x] Define required evaluation figures
- [x] Obtain explicit user approval before entering Phase 4
- [x] Run Phase 4 and write continuous score columns back into `integrated.h5ad`
- [ ] Review the completed Phase 4 QC package with the user
- [x] Commit and push milestone changes to GitHub

---

## Canonical output files
- `Integrated_dataset/TNK_candidates.h5ad`
- `Integrated_dataset/TNK_cleaned.h5ad`
- `Integrated_dataset/integrated.h5ad`

---

## Required figure outputs
Store all PNG figures in:
- `Integrated_dataset/figures/`

Examples:
- QC summary
- batch mixing summary
- UMAP overview
- marker overlays
- contamination review plots
- Phase 4 scoring plots

---

## Files changed
- `AGENTS.md`
- `TNK_PIPELINE_RUNBOOK.md`
- `TNK_PIPELINE_STATUS.md`
- `TNK_PHASES_1_3_SCRIPT.md`
- `phase0_dataset_audit.py`
- `repair_h5ad_from_raw.py`
- `repair_h5ad_from_selected_inputs.py`
- `phase1_extract_tnk_candidates.py`
- `phase1_finalize_from_temp.py`
- `phase1b_conservative_cleanup.py`
- `phase1c_replace_harmonized_metadata.py`
- `phase2_merged_cleanup.py`
- `phase3_scvi_scanvi.py`
- `phase4_gdt_module_scoring.py`
- `watch_h5ad_v2_and_resume.py`

---

## Outputs created
- Initial Codex control files
- `TNK_PHASES_1_3_SCRIPT.md`
- `phase0_dataset_audit.py`
- `Integrated_dataset/tables/phase0_dataset_audit.csv`
- `Integrated_dataset/tables/phase0_category_summary.csv`
- `Integrated_dataset/tables/h5ad_without_raw_count.csv`
- `Integrated_dataset/logs/phase0_qc_summary.md`
- `Integrated_dataset/logs/GSE228597_raw_rescue.md`
- `Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651/backup_manifest.tsv`
- `Integrated_dataset/tables/h5ad_without_raw_count_source_audit_20260320.csv`
- `Integrated_dataset/tables/h5ad_without_raw_count_source_audit_20260320.md`
- `Integrated_dataset/tables/GSE190870_selected_inputs_rescue_dryrun_20260320.md`
- `Integrated_dataset/tables/GSE252762_selected_inputs_rescue_dryrun_20260320.md`
- `Integrated_dataset/tables/h5ad_rescue_status_20260320.csv`
- `Integrated_dataset/tables/h5ad_selected_inputs_final_validation_20260320.csv`
- `Integrated_dataset/tables/h5ad_drop_current_only_write_results_20260320.csv`
- `Integrated_dataset/tables/remaining_gene_mismatches_after_drop_20260320.csv`
- `h5ad_backup_before_v2_20260320.csv`
- `h5ad_v2.csv`
- `Integrated_dataset/tables/phase1_categoryA_selection_summary.csv`
- `Integrated_dataset/tables/phase1_categoryA_marker_availability.csv`
- `Integrated_dataset/logs/phase1_qc_summary.md`
- `Integrated_dataset/TNK_candidates.h5ad`
- `Integrated_dataset/tables/phase1b_cleanup_summary.csv`
- `Integrated_dataset/tables/phase1b_gse_before_after.csv`
- `Integrated_dataset/tables/phase1b_removed_cells.csv`
- `Integrated_dataset/tables/phase1b_gene_detection_summary.csv`
- `Integrated_dataset/logs/phase1b_qc_summary.md`
- `Integrated_dataset/tables/phase1c_merged_obs_export.csv.gz`
- `Integrated_dataset/tables/phase1c_metadata_replacement_summary.csv`
- `Integrated_dataset/tables/phase1c_metadata_join_by_gse.csv`
- `Integrated_dataset/logs/phase1c_metadata_replacement.md`
- `analysis_26GSE_V4/outputs/metadata.csv.bk`
- `Integrated_dataset/figures/phase0_category_distribution.png`
- `Integrated_dataset/figures/phase0_matrix_state_overview.png`
- `Integrated_dataset/figures/phase0_dataset_size_overview.png`
- `Integrated_dataset/figures/phase0_metadata_completeness.png`
- `Integrated_dataset/figures/phase1_categoryA_candidate_yield.png`
- `Integrated_dataset/figures/phase1_categoryA_candidate_fraction.png`
- `Integrated_dataset/figures/phase1_categoryA_marker_support.png`
- `Integrated_dataset/figures/phase1b_gse_cell_retention.png`
- `Integrated_dataset/figures/phase1b_cell_removal_reasons.png`
- `Integrated_dataset/figures/phase1b_gene_detection_distribution.png`
- `downloads/per_gse_h5ad_with_metadata/GSE179994_with_tcr.h5ad`
- `downloads/per_gse_h5ad_with_metadata/GSE234069_with_tcr.h5ad`
- `downloads/per_gse_h5ad_with_metadata/GSE235863_with_tcr.h5ad`
- `downloads/per_gse_h5ad_with_metadata/GSE240865_with_tcr.h5ad`
- `downloads/per_gse_h5ad_with_metadata/GSE287301_with_tcr.h5ad`
- `downloads/per_gse_h5ad_with_metadata/GSE287541_with_tcr.h5ad`
- `Integrated_dataset/tables/supplementary_10x5/h5ad_supplementary_10x5.csv`
- `Integrated_dataset/tables/supplementary_10x5/supplementary_h5ad_build_summary.csv`
- `Integrated_dataset/tables/supplementary_10x5/supplementary_metadata_all_cells.csv.gz`
- `Integrated_dataset/tables/supplementary_10x5/phase0_dataset_audit.csv`
- `Integrated_dataset/tables/supplementary_10x5/phase0_category_summary.csv`
- `Integrated_dataset/logs/supplementary_10x5/phase0_qc_summary.md`
- `Integrated_dataset/figures/supplementary_10x5/phase0_category_distribution.png`
- `Integrated_dataset/figures/supplementary_10x5/phase0_dataset_size_overview.png`
- `Integrated_dataset/figures/supplementary_10x5/phase0_matrix_state_overview.png`
- `Integrated_dataset/figures/supplementary_10x5/phase0_metadata_completeness.png`
- `Integrated_dataset/TNK_candidates_supp.h5ad`
- `Integrated_dataset/tables/supplementary_10x5/phase1_categoryA_selection_summary.csv`
- `Integrated_dataset/tables/supplementary_10x5/phase1_categoryA_marker_availability.csv`
- `Integrated_dataset/logs/supplementary_10x5/phase1_qc_summary.md`
- `Integrated_dataset/figures/supplementary_10x5/phase1_categoryA_candidate_yield.png`
- `Integrated_dataset/figures/supplementary_10x5/phase1_categoryA_candidate_fraction.png`
- `Integrated_dataset/figures/supplementary_10x5/phase1_categoryA_marker_support.png`
- `analysis_26GSE_V4/outputs/harmonized_metadata_supp.csv`
- `Integrated_dataset/tables/merged_tnk_candidates_with_supp_counts.csv`
- `Integrated_dataset/logs/merged_tnk_candidates_with_supp.md`
- `Integrated_dataset/tables/phase4_score_summary.csv`
- `Integrated_dataset/tables/phase4_module_gene_membership.csv`
- `Integrated_dataset/tables/phase4_leiden_score_summary.csv`
- `Integrated_dataset/tables/phase4_gse_score_summary.csv`
- `Integrated_dataset/tables/phase4_top_cells_by_trd_minus_trab.csv`
- `Integrated_dataset/logs/phase4_qc_summary.md`
- `Integrated_dataset/figures/phase4_score_distributions.png`
- `Integrated_dataset/figures/phase4_umap_score_overlays.png`
- `Integrated_dataset/figures/phase4_leiden_score_summary.png`
- `Integrated_dataset/figures/phase4_gse_score_summary.png`
- `Integrated_dataset/figures/phase4_marker_score_comparison.png`

---

## Validation / checks
- `phase0_dataset_audit.py` compiled successfully in `rapids_sc_py310`
- `repair_h5ad_from_raw.py` compiled successfully in `rapids_sc_py310`
- Phase 0 audit completed on all 33 registry rows
- Read errors recorded: 0
- Updated category counts after `GSE228597` repair: A=13, B=0, C=20
- Metadata coverage counts: donor=31, sample=32, library=30
- `GSE228597` validated post-repair with `has_raw=False`, `obsm=[]`, `varm=[]`, integer-like `X`, and `max(X sample)=19`
- All 20 entries in `Integrated_dataset/tables/h5ad_without_raw_count.csv` were re-checked and confirmed to lack internal raw/count layers in the current H5ADs
- All 20 entries were backed up under `Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- All 20 entries have matching `selected_inputs.csv` and `per_sample_load_stats.csv` under `analysis_26GSE_V4/scanpy_projects/<GSE>/`
- Supplementary raw tarballs were found for 18 of the 20 entries; the remaining 2 were rescued from direct matrix sources in `downloads/GSE*/suppl`
- Dry-run rescue on `GSE252762` showed `obs_equal=True` and `var_current_subset_of_rebuilt=True`, which is sufficient for a controlled count-space repair
- Dry-run rescue on `GSE190870` showed `obs_equal=True` but `var_current_subset_of_rebuilt=False`, indicating a remaining gene-space reconciliation problem before safe repair
- Final selected-input rescue validation confirms all 20 repaired H5ADs with integer-like sampled `X`, `has_raw=False`, `obsm=[]`, and `varm=[]`
- `GSE190870` was rescued after dropping 41 current-only unmatched genes before count restoration
- `GSE267645` was rescued after dropping 354 current-only unmatched genes before count restoration
- `GSE254249` was revalidated as integer-like count-space and restored to the registry
- `GSE301528` was repaired by rebuilding counts from the 16 valid selected inputs and skipping the 2 dimension-mismatched source matrices
- Supplementary Phase 0 on the six approved 10x 5' datasets completed with 6/6 readable files, 6/6 Category A, 0 read errors, and no remaining duplicate `obs_names` in the final per-GSE H5AD outputs
- Supplementary Phase 1 on the six approved 10x 5' datasets retained 897,621 / 1,032,943 cells (86.90%) into the separate milestone `TNK_candidates_supp.h5ad`
- Main `TNK_candidates.h5ad` now contains 6,848,556 cells and 55,455 genes after outer-join merge of the supplementary candidate object
- `h5ad_v2.csv` now retains all 33 datasets from the canonical registry
- `phase1_extract_tnk_candidates.py` compiled successfully in `rapids_sc_py310`
- `phase1_finalize_from_temp.py` compiled successfully in `rapids_sc_py310`
- `phase1b_conservative_cleanup.py` compiled successfully in `rapids_sc_py310`
- `phase1c_replace_harmonized_metadata.py` compiled successfully in `rapids_sc_py310`
- Phase 0 rerun on `h5ad_v2.csv` classified all 29 retained datasets as Category A
- Phase 1 candidate object rebuilt and validated with `n_obs=5322388` and `n_vars=57093`
- Phase 1 retained fraction across the repaired 29-dataset registry: `0.8608`
- Per-dataset candidate counts in `TNK_candidates.h5ad` match `phase1_categoryA_selection_summary.csv`, including `GSE227709=0`
- Temporary Phase 1 subset directory was removed after successful validation
- Phase 1b validated the rewritten `TNK_candidates.h5ad` at `n_obs=5950935` and `n_vars=23536`
- Phase 1b removed `261` obvious non-T/NK cells and `0` high-confidence doublets
- The Phase 1b `<500 cells` gene filter retained all `11` tracked key T/NK genes
- Phase 1c joined the filtered candidate obs to harmonized metadata with a complete string-typed `GSE + barcode` match
- Phase 1c backed up the previous metadata target to `analysis_26GSE_V4/outputs/metadata.csv.bk`
- Phase 1c replaced `analysis_26GSE_V4/outputs/harmonized_metadata_v4.csv` with `5950935` rows and `0` duplicated `project name + sampleid + barcodes` keys
- Phase 1c post-write validation found `0` blank `project name`, `0` blank `barcodes`, and `58678` blank `sampleid`

---

## Open issues
- Several project manifests point to stale extracted directories; rescue tooling must resolve the real files under `downloads/GSE*/suppl` rather than trusting manifest paths blindly
- Some datasets will likely be directly repairable by rebuilding counts from selected inputs, while others will need gene-space reconciliation first
- Whether `rapids_sc_py310` should remain the working env or be cloned/aliased to `Scanpy_gdTmodel`
- `sampleid` remains blank for `58678` rows after trusted fill from prior metadata and Phase 1 labels; decide whether Phase 2 needs manual curation first
- The Phase 3 reference model should be revalidated at execution time for query compatibility and label-space suitability before scANVI mapping
- Exact Phase 4 scoring inputs/outputs for `gdt_tcr_module_sharing_package_full`

---

## Git status
- Milestone documentation and workflow scripts committed as `4e64f02`
- `master` pushed to `origin/master`
- Auto-resume fix for zero-cell Phase 1 datasets has been applied locally

---

## Next action
- Present the Phase 1c replacement summary and updated metadata outputs to the user
- Wait for explicit user approval before Phase 2 cleanup
- If approved, begin Phase 2 merged cleanup on the current milestone set

---

## Notes
If the session feels context-compressed or uncertain, re-read:
1. `TNK_PIPELINE_RUNBOOK.md`
2. `TNK_PIPELINE_STATUS.md`
3. `TNK_PHASES_1_3_SCRIPT.md`
before continuing.
