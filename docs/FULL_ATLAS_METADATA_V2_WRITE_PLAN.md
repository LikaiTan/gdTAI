# Full-Atlas Metadata Write Plan

## Purpose

Apply the validated additive tissue, specimen, tumor, sample, library, donor,
and provenance columns before any repaired TCR chain calls are joined.

## Frozen Inputs

- source atlas:
  `/ssd/tnk_phase3/Integrated_dataset/full_atlas/integrated_full_atlas.h5ad`
- source atlas SHA-256:
  `f5fc491a70f12adeeda5764cb116a59bb441460285905f297bbc2ff691559802`
- metadata overlay:
  `Integrated_dataset/tables/metadata_harmonization/full_atlas_v2/metadata_overlay_v2.parquet`
- metadata overlay SHA-256:
  `4da4eea32e9d275de790775e2a0f59d4f6553d72756e9c8c6935f35bb398984f`
- overlay rows: 5,933,312
- overlay sources: 40

The current canonical path resolves to the source atlas above. The source file
is the rollback object and must remain byte-identical.

## Candidate Path

- temporary candidate:
  `/ssd/tnk_phase3/Integrated_dataset/full_atlas/metadata_corrected/integrated_full_atlas.h5ad.partial`
- validated candidate:
  `/ssd/tnk_phase3/Integrated_dataset/full_atlas/metadata_corrected/integrated_full_atlas.h5ad`

The SSD currently has approximately 1.8 TB free, well above the source-atlas
size plus the required safety reserve.

## Write Scope

Create a separate copy of the source atlas and append only these columns:

- `source_accession_harmonized_v2`
- `tissue_harmonized_v2`
- `specimen_context_harmonized_v2`
- `tumor_type_harmonized_v2`
- `sample_id_harmonized_v2`
- `library_id_harmonized_v2`
- `donor_id_harmonized_v2`
- `sample_identity_rule_v2`
- `tcr_library_join_id_v2`
- `metadata_rule_id_v2`
- `metadata_evidence_url_v2`
- `metadata_evidence_level_v2`
- `metadata_status_v2`

Join only by `source_gse_id + original_cell_id`. Do not use row position,
barcode alone, current display sample name, or technical VDJ-library identity
as the biological-sample key.

No existing metadata column is replaced. No TRA/TRB/TRG/TRD field, logical TCR
flag, expression value, embedding, cluster, or score is changed in this step.

## Required Validation

Before renaming the partial file to the validated candidate:

- source atlas, overlay, configuration, and reviewed legacy metadata checksums
  match the frozen values
- all 5,933,312 source/cell keys match one-to-one
- all 13 columns reproduce the overlay value counts and per-column ordered
  value hashes
- no new sample or library value is blank
- exactly 4,611 biological sample IDs are explicitly `unresolved`
- exactly 39,920 tissue and specimen-context values remain `unresolved`
- no cell-type label occurs in a tissue or specimen-context column
- sparse `X` shape, nonzero count, storage properties, and deterministic sample
  hashes match the source atlas
- `X_scVI`, `X_umap`, Leiden, Phase 4 scores, `var`, and cell order match the
  source atlas
- every existing TCR column and TCR logical flag matches the source atlas
- the source atlas size, mtime, inode, and SHA-256 remain unchanged
- AnnData can open the candidate in backed read-only mode

## Phase Boundary

Passing this validation makes the metadata-corrected candidate eligible as the
input to the TCR sidecar transaction. It does not move the canonical symlink.

The subsequent TCR transaction will create another candidate from this
validated metadata object, apply the checksum-bound 14-source TCR sidecar, and
validate TCR fields separately. Only the fully validated combined object may be
considered for atomic canonical publication.

## Rollback

Before canonical publication, rollback is simply deletion of an invalid
candidate; the current canonical source was never changed. After a future
atomic symlink switch, rollback restores the symlink to the original source
atlas path recorded above.

The candidate write is a high-risk operation and requires explicit approval.
