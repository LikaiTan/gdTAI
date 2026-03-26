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

- `obs["tissue_corrected"]` is already written into the current integrated
  milestone
- repaired TCR propagation is complete for ten approved GSEs
- scANVI fields remain reference-only; simple scVI-based annotation is the
  canonical downstream interpretation layer

## Current review artifacts

- Phase 4 QC:
  - `Integrated_dataset/logs/phase4_qc_summary.md`
- tissue correction:
  - `Integrated_dataset/tables/tissue_correction/tissue_correction_apply.md`
- repaired TCR review:
  - `Integrated_dataset/logs/tcr_rebuild_phase4_qc.md`
- disease-status export review:
  - `Integrated_dataset/tables/disease_status_correction/disease_status_corrected_column_export.csv.gz`
