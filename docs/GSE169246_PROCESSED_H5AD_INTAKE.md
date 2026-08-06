# GSE169246 Processed H5AD Intake

## Selected Input

At the user's direction on 2026-08-06, the selected processed input for
GSE169246 is:

`/home/leeaaron/Exhausted_GDT/TNBC_validation/Integrated_dataset/TNK_cleaned.h5ad`

It was copied without rewriting to:

`data/datasets/GSE169246/processed/artifacts/TNK_cleaned.h5ad`

The stable project entry is:

`data/datasets/GSE169246/processed/current.h5ad`

## Integrity And Structure

- SHA-256: `2038797764ff8fcb655e81886c059cf2c25e04157fb8983625b0d0c985f2e333`
- Size: `2,436,744,413` bytes
- Shape: `394,225` cells by `12,989` genes
- `source_gse_id`: `GSE169246` for all cells
- Sparse integer counts: present in both `X` and `layers["counts"]`
- Embedded paired TRA/TRB CDR3: `54,925` cells
- Productive TRG/TRD metadata: not present, as expected for the supplied
  alpha-beta TCR metadata

## Provenance Boundary

This is a processed T/NK candidate object, not the complete unfiltered study.
Its own provenance records the GSE169246 count matrix plus `mmc3.xls` cell
metadata and `mmc5.xls` alpha-beta TCR metadata joined by sample and barcode.
The object contains `279,818` cells without a published `major_cluster` label,
so published-annotation coverage remains incomplete even though every cell is
accession-labeled GSE169246.

The user-selected file supersedes the ambiguous bootstrap file as the input to
use for this cohort. It does not by itself approve merge or integration.

## Gate State

- Canonical copy and checksum validation: complete
- Dataset-registry entry: complete and inactive
- Standalone metadata/TCR-schema harmonization: pending
- Standalone Phase 0 QC and user review: pending
- Merge, integration, or gdTAI evaluation: not approved

Pre-change registry state is recoverable from:

`data/registry/snapshots/pre_gse169246_processed_h5ad_intake_20260806/`

Validated post-change registry state is recorded at:

`data/registry/snapshots/post_gse169246_processed_h5ad_intake_20260806/`
