# GSE305372 External Application Retirement

## Status

- Retired and excluded on 2026-08-06 by explicit user decision.
- Not available to active gdTAI evaluation, training, atlas integration, or
  prevalence estimation.
- Original local dataset directory deleted after source checksums and recovery
  information were retained.

## Reason

The deposited processed objects are finalized CD4/CD8 objects produced after
lineage-specific deconvolution, and the deposited lymph-node objects are only a
partial representation of the manuscript cohort. The resulting gdT frequency
cannot be interpreted as an unbiased prevalence or model-validation estimate.

## Archived Material

- `code/`: exact dataset-specific download, payload-export, and inference code
- `docs/`: methodology records for the corrected all-T and superseded CD8-only
  applications
- `local_results/`: locally retained generated tables, figures, and logs; this
  directory is intentionally excluded from Git because it is reproducible and
  contains per-cell scientific output
- rendered reports:
  `archive/historical_reports/GSE305372_external_application/`
- source download manifest: `deleted_source_files.csv`
- move/deletion plan:
  `maintenance/reorganization/gse305372_retirement_plan.csv`
- completed result map:
  `maintenance/reorganization/gse305372_retirement_map.csv`

## Rollback And Recovery

Tracked code and reports immediately before retirement are available at Git
commit `8e23824` on branch `codex/gdtai-v2-operating-modes`.

Restoring the scientific data requires re-downloading the five GEO resources
listed in `deleted_source_files.csv`. Restore the archived scripts to their
original paths from the retirement map before running them; the scripts retain
their original project-relative path assumptions. Validate both byte size and
SHA-256 before analysis.

Reactivation is not automatic after restoration. It requires an explicit new
decision and a documented use case that does not treat the finalized CD4/CD8
objects as an unbiased gdT prevalence cohort.
