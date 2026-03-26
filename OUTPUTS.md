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
- disease-status review export:
  - `Integrated_dataset/tables/disease_status_correction/`
