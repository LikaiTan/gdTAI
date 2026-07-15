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
- pre-integration data organization uses a registry-driven symlink view during
  migration; physical raw-data moves require dataset-level full checksums and
  explicit cutover review
- dataset removal means registry deactivation, not file deletion
- `configs/models/gdtai/model_registry.csv` is authoritative for model release
  status; a directory or version-like filename alone does not imply promotion
- canonical Git gdTAI v3 Round 14 remains the promoted default; the existing
  dirty-workspace Round 12 override is preserved and explicitly not promoted
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
