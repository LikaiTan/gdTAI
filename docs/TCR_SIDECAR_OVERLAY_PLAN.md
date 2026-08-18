# TCR Sidecar Overlay and Atlas-Rebuild Plan

## Decision

Do not write rebuilt TCR calls into the 14 source H5AD files. Treat the
validated Parquet sidecar as an immutable metadata overlay on the frozen
5,933,312-cell atlas membership.

The read-only overlay preflight found 20,875 productive alpha-beta T cells
outside the current T/NK-filtered atlas, including 14,495 cells with paired
TRA/TRB. None has productive gamma-delta TCR. These rows remain available in
their source H5ADs but will not be rescued into the integrated atlas. This
explicitly accepted loss is small relative to the atlas, contributes only
alpha-beta controls rather than gamma-delta positives, and avoids changing the
cell universe, HVGs, scVI latent space, clusters, and UMAP solely to recover
additional alpha-beta controls.

## Frozen Inputs

- current atlas:
  `Integrated_dataset/integrated_full_atlas.h5ad`
- current scored-atlas SHA-256:
  `f5fc491a70f12adeeda5764cb116a59bb441460285905f297bbc2ff691559802`
- validated TCR sidecar:
  `Integrated_dataset/tables/tcr_join_rebuild/validated_tcr_replacement_sidecar.parquet`
- sidecar SHA-256:
  `3114e70719301d693ae1a2bc2c63bac6c8bd57e3e8ac73a88c24320eaabfc2f0`
- frozen Phase 3 exclusion manifest:
  `Integrated_dataset/tables/phase3_input_sanitization.csv`

## Read-Only Preflight Result

- sidecar rows: 3,041,871
- affected sources: 14
- exact atlas matches by `source_gse_id + original_cell_id`: 2,155,409
- affected atlas rows without a sidecar: 0
- sidecar rows outside the current T/NK-filtered atlas: 886,462
- sidecar-only cells with rebuilt productive TCR: 20,875
- sidecar-only cells with paired TRA/TRB: 14,495
- sidecar-only cells with productive gamma-delta TCR: 0
- retained atlas cells with productive alpha-beta TCR after repair: 1,100,981
- retained atlas cells with paired TRA/TRB after repair: 933,890
- estimated whole-atlas productive-alpha-beta cells after overlay: 2,270,138
- estimated whole-atlas paired-TRA/TRB cells after overlay: 1,938,158
- fail-closed ambiguous GSE235863 rows: 110 total; 109 currently in the atlas

The 20,875 intentionally omitted productive-alpha-beta rows are 0.352% of the
current atlas and 1.861% of all available repaired productive-alpha-beta
controls within the 14 repaired sources. The 14,495 omitted paired cells are
1.528% of paired controls within those sources. When unaffected sources are
included, the losses are only 0.911% of the available productive-alpha-beta
pool and 0.742% of the available paired-alpha-beta pool. Their source-level
counts must remain in the published QC report so that this is a visible
atlas-membership decision rather than silent loss.

The normalized sidecar `join_sample_id` must remain provenance, not the atlas
join key. It intentionally differs from display-oriented atlas `sample_id` in
many sources. Positional, barcode-only, and display-sample joins are forbidden.

## Stage 1: Frozen-Membership Overlay Manifest

This stage is blocked until the additive tissue, tumor, sample, library, and
donor metadata replacement has been applied and validated. The TCR transaction
must consume that validated metadata object, not the current pre-correction
atlas.

1. Open the canonical atlas and all source H5AD files read-only.
2. Verify every source file against the recorded size and nanosecond mtime.
3. Verify the atlas and sidecar SHA-256 values and the 14-source manifest.
4. Join the corrected-metadata atlas only by `source_gse_id + original_cell_id` to sidecar
   `source_gse_id + source_obs_name`.
5. Require uniqueness and complete sidecar coverage for every atlas row from
   an affected source.
6. For replacement-eligible rows, replace all canonical TRA/TRB/TRG/TRD fields
   from the sidecar, including empty values that remove stale calls.
7. For TCR-assignment-ineligible rows, clear all canonical chain calls and
   quantitative support, set all TCR flags false, and retain an explicit
   fail-closed status.
8. Keep unavailable UMI/read support null rather than zero.

The source H5AD files remain unchanged. The old atlas is retained as the
complete legacy-metadata snapshot rather than duplicating dozens of legacy TCR
columns across almost six million rows.

`tcr_library_join_id_v2` is technical provenance and is not interchangeable
with `sample_id_harmonized_v2`. This is essential for pooled sources such as
GSE228597.

## Stage 2: Frozen-Membership Validation

Do not add or remove cells. Required audit groups are:

- current cells whose stale TCR call is removed
- current cells newly assigned a validated productive-TCR call
- current cells whose old and rebuilt TCR evidence agrees
- the 109 atlas-present fail-closed GSE235863 rows
- the 886,462 sidecar rows outside the atlas, including the intentionally
  omitted 20,875 productive-alpha-beta and 14,495 paired-alpha-beta rows

Recompute canonical chain fields, logical flags, and TCR truth/control labels
for retained cells. A cell whose stale call is removed cannot remain a TCR
ground-truth cell merely because its atlas membership is frozen.

## Stage 3: Corrected Atlas Metadata Build

1. Preserve the exact 5,933,312-cell identity, order, sparse expression matrix,
   27,413-gene universe, scVI latent space, clusters, and UMAP.
2. Replace only the canonical TCR metadata and derived TCR flags for the 14
   affected sources.
3. Preserve legacy TCR values under explicitly named provenance columns only
   when required for auditing; downstream truth must use rebuilt fields.
4. Recalculate any score or label that directly depends on TCR metadata.
   Expression-only Phase 4 scores do not require recomputation.
5. Apply approved additive tissue/tumor metadata columns after their independent
   evidence review; do not overwrite original tissue or diagnosis fields.

## Stage 4: Atomic Publication and Rollback

Before writing:

- copy the current 22.3-GB atlas to a timestamped SSD backup
- record the whole-file SHA-256 plus HDF5 component hashes for `X`, `var`,
  `X_scVI`, and `X_umap`
- record the current canonical symlink target

Write the corrected atlas to a `.partial` path. Validate it before promotion:

- exact frozen cell-manifest identity and ordering
- unique cell IDs and source-original-cell keys
- sparse matrix shape, nnz, and finite embeddings
- all chain fields and logical flags recomputed consistently
- all 109 atlas-present ambiguous rows cleared and ineligible
- no repaired call on an ineligible row
- UMI/read null semantics preserved
- source H5AD size/mtime pairs unchanged
- HTML/PDF QC report rendered and visually checked

Only after every check passes may the canonical symlink switch atomically to
the corrected atlas. Rollback consists of restoring the recorded symlink target;
the current atlas and all source H5ADs remain intact.

## Downstream Boundary

After corrected-atlas publication:

- rebuild productive-T controls and NK reference audits
- regenerate TCR ground-truth labels and false-positive denominators
- compare gdTAI V2/V3 on unchanged locked cohorts
- do not train or promote V4 until grouped validation and the repaired metadata
  audits pass

The metadata-only H5AD write and canonical symlink switch are high-risk
execution steps and require explicit approval after the dry-run replacement
manifest is available. A future atlas version may revisit the 20,875 omitted
alpha-beta cells, but they are not part of this repair.
