# DECISIONS.md

## Current approved decisions

- current canonical execution environment is `rapids_sc_py310`
- Phase 5 model training and validation are out of scope unless explicitly requested
- scANVI outputs remain in the integrated object for reference only
- simple scVI-based annotation is the canonical downstream interpretation layer
- Phase 4 uses package-faithful continuous TRA/TRB/TRD scoring
- `obs["tissue_corrected"]` is canonical in the current integrated milestone
- large H5AD files may remain on the mirrored SSD tree until explicit migration approval
- tables, PNG figures, logs, scripts, and model artifacts stay on NFS
- canonical pre-integration physical storage is
  `data/datasets/<dataset_id>/`; `data/raw`, `data/interim`, and
  `data/processed` are lifecycle compatibility views
- historical `downloads`, `analysis_26GSE_V4`, and `newdata` paths remain
  supported aliases and must not be removed while legacy workflows exist
- physical storage changes use journaled same-filesystem renames with
  inode/hash validation and a tested reverse rollback plan
- dataset removal means registry deactivation, not file deletion
- `configs/models/gdtai/model_registry.csv` is authoritative for model release
  status; a directory or version-like filename alone does not imply promotion
- gdTAI v3 Round 14 at threshold `0.936` is the promoted balanced default after
  checksum-pinned revalidation against Round 12 on 2026-07-15
- gdTAI v3 Round 12 at threshold `0.5` is retained as the validated high-purity
  fallback and is not the canonical default
- the independent BALF/PBMC COPD cohort is registered as
  `BALF_BLOOD_COPD`; its canonical H5AD is
  `data/datasets/BALF_BLOOD_COPD/processed/current.h5ad`
- `BALF_BLOOD_COPD` is validation-only and must remain inactive for Phase 0,
  milestone integration, and extended-atlas integration unless the user
  separately approves a role change
- `/home/tanlikai/databank/owndata/singlecell` remains a supported
  compatibility alias so the cohort's original scripts and RStudio project
  paths continue to work
- repaired TCR propagation is approved for these ten GSEs:
  - `GSE188620`
  - `GSE212217`
  - `GSE243572`
  - `GSE243905`
  - `GSE254249`
  - `GSE308075`
  - `GSE311112`
  - `GSE161918`
  - `GSE168859`
  - `GSE227709`

## Historical exceptions and notes

- a historical run-specific exception allowed Phase 2 internal QC to continue
  directly into Phase 3 on 2026-03-20
- that exception is historical only and is not currently active
- the active RAM ceiling override for this run is `800 GB`
