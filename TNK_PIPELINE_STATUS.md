# TNK_PIPELINE_STATUS.md

## Current milestone

- Post-Phase 4 review on the validated integrated milestone

## Next action

- Continue user-directed downstream analyses and reporting on the current
  integrated object

## Current blockers or review items

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

- the non-destructive project reorganization is validated: active code is
  grouped under `workflows/`, configuration under `configs/`, and the
  pre-integration `data/` lifecycle view resolves through registry-driven
  symlinks without changing legacy source paths or milestone H5ADs
- the dataset registry currently contains 65 datasets and 1,785 library rows;
  all 33 Phase 0 active datasets pass strict path and library validation
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

## Current review artifacts

- project reorganization:
  - `docs/REORGANIZATION_RECORD.md`
  - `Integrated_dataset/logs/project_reorganization/checkpoints/pre_reorganization_20260715/`
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
