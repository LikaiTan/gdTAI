# TNK_PIPELINE_STATUS.md

## Project
Large-scale T/NK integration and γδT-focused scoring workflow across public datasets

---

## Current milestone
- Phase 1 QC review pending

## Current objective
- Review Phase 1 Category A coarse extraction outputs
- Decide whether to proceed to Phase 1b conservative cleanup on `TNK_candidates.h5ad`
- Confirm whether any Phase 1 threshold adjustments are needed before Phase 1b
- Confirm environment naming policy: continue with `rapids_sc_py310` or create/alias `Scanpy_gdTmodel`

---

## Completed milestones
- Initial control-file design completed
- Canonical Phase 0-3 workflow file created
- Phase 0 dataset audit completed across all 33 registry datasets
- `GSE228597` rescued from `adata.raw` and reclassified into Category A
- Backed up all 20 Phase 0 Category C H5ADs before external raw-source rescue review
- Completed external raw-source audit for the 20 Category C H5ADs
- Phase 1 coarse extraction completed for the 13 current Category A datasets

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
- [ ] Review Phase 1 QC outputs with the user
- [ ] Decide whether to proceed to Phase 1b conservative cleanup
- [ ] Define or adjust Phase 1b cleanup rules if needed
- [ ] Review the 20 Category C raw-source audit with the user and approve dataset-by-dataset rescue scope
- [ ] Define Phase 4 scoring workflow
- [ ] Define required evaluation figures
- [ ] Commit and push milestone changes to GitHub

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
- `Integrated_dataset/tables/phase1_categoryA_selection_summary.csv`
- `Integrated_dataset/tables/phase1_categoryA_marker_availability.csv`
- `Integrated_dataset/logs/phase1_qc_summary.md`
- `Integrated_dataset/TNK_candidates.h5ad`
- `Integrated_dataset/figures/phase0_category_distribution.png`
- `Integrated_dataset/figures/phase0_matrix_state_overview.png`
- `Integrated_dataset/figures/phase0_dataset_size_overview.png`
- `Integrated_dataset/figures/phase0_metadata_completeness.png`
- `Integrated_dataset/figures/phase1_categoryA_candidate_yield.png`
- `Integrated_dataset/figures/phase1_categoryA_candidate_fraction.png`
- `Integrated_dataset/figures/phase1_categoryA_marker_support.png`

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
- Supplementary raw tarballs were found for 18 of the 20 entries; the remaining 2 still have direct matrix sources in `downloads/GSE*/suppl`
- Dry-run rescue on `GSE252762` showed `obs_equal=True` and `var_current_subset_of_rebuilt=True`, which is sufficient for a controlled count-space repair
- Dry-run rescue on `GSE190870` showed `obs_equal=True` but `var_current_subset_of_rebuilt=False`, indicating a remaining gene-space reconciliation problem before safe repair
- `phase1_extract_tnk_candidates.py` compiled successfully in `rapids_sc_py310`
- `phase1_finalize_from_temp.py` compiled successfully in `rapids_sc_py310`
- Phase 1 candidate object validated with `n_obs=2689420` and `n_vars=32172`
- Phase 1 retained fraction across Category A inputs: `0.9126`
- Per-dataset candidate counts in `TNK_candidates.h5ad` match `phase1_categoryA_selection_summary.csv`

---

## Open issues
- Several project manifests point to stale extracted directories; rescue tooling must resolve the real files under `downloads/GSE*/suppl` rather than trusting manifest paths blindly
- Some datasets will likely be directly repairable by rebuilding counts from selected inputs, while others will need gene-space reconciliation first
- Whether the current high-recall Phase 1 thresholds are acceptable as input to Phase 1b
- Whether `rapids_sc_py310` should remain the working env or be cloned/aliased to `Scanpy_gdTmodel`
- Whether all TCR-related genes are preserved across datasets
- Exact Phase 4 scoring inputs/outputs for `gdt_tcr_module_sharing_package_full`

---

## Git status
- No milestone commit recorded yet

---

## Next action
- Present the Phase 1 QC summary, tables, figures, and `TNK_candidates.h5ad` to the user
- Wait for explicit user approval before any Phase 1b cleanup
- If approved, run conservative first-pass cleanup in place on `TNK_candidates.h5ad`

---

## Notes
If the session feels context-compressed or uncertain, re-read:
1. `TNK_PIPELINE_RUNBOOK.md`
2. `TNK_PIPELINE_STATUS.md`
3. `TNK_PHASES_1_3_SCRIPT.md`
before continuing.
