# TCR Integration SOP

## Purpose

This document defines the canonical workflow for integrating external TCR/VDJ
metadata into single-cell RNA objects in this project.

It is intentionally framework-agnostic:

- use the same join rules for AnnData / Scanpy workflows
- use the same join rules for Seurat workflows
- use the same join rules for supplementary 10x 5' intake

The main goal is to prevent false TCR assignment caused by barcode reuse across
samples.

## Non-negotiable rules

1. Never join TCR by barcode alone.
2. Never join TCR by barcode prefix alone.
3. Always standardize one canonical sample key first.
4. Always normalize barcodes to one canonical `barcode_core`.
5. Always join on `sample_id + barcode_core`.
6. Always replace stale TCR fields rather than mixing old and rebuilt values.
7. Always validate the rebuilt TCR assignment before propagation.

## Why this SOP exists

Several datasets showed implausibly high paired `TRA/TRB` coverage when TCR
metadata was joined with barcode-only or barcode-prefix-only logic. Rebuilding
the TCR tables with a sample-aware key reduced those artifacts substantially.

The safe default is therefore:

- build a canonical per-cell TCR table from raw contig or VDJ outputs
- normalize the join key to `sample_id + barcode_core`
- join only on that key

## Step 1: Build one canonical per-cell TCR table

Start from the raw TCR source when possible:

- `filtered_contig_annotations.csv`
- `all_contig_annotations.csv`
- equivalent VDJ contig exports

Use productive chains only unless there is a project-specific reason not to.

Deduplicate within chain using a deterministic rule, for example:

- highest `umis`
- then highest `reads`

The canonical per-cell TCR table must contain:

- `sample_id`
- `barcode_core`
- `TRA_*`
- `TRB_*`
- `TRG_*`
- `TRD_*`

Do not assume that every dataset contains only `TRA/TRB`.

Datasets may contain:

- alpha-beta only
- gamma-delta only
- both alpha-beta and gamma-delta
- partial chain information

## Step 2: Canonical sample-key normalization

Do not assume the sample column is literally named `sample_id`.

Candidate columns may include:

- `sample_id`
- `sampleid`
- `sample`
- `orig.ident`
- `library_id`
- `library`
- `batch`
- dataset-specific aliases

Required behavior:

- discover candidate sample-like columns explicitly
- normalize to one canonical string column named `sample_id`
- record dataset-specific mapping rules when the sample naming is indirect
- stop if sample identity remains ambiguous

Do not silently fall back to barcode-only joins when sample identity is unclear.

## Step 3: Canonical barcode normalization

Normalize every cell barcode to one canonical `barcode_core`.

Typical behavior:

- keep the 10x oligo sequence before lane or suffix decorations
- for barcodes like `AAACCTGAGAACAATC-1-0-1`, normalize to `AAACCTGAGAACAATC`

This is necessary for cross-format matching, but it is not sufficient alone.
`barcode_core` must always be paired with `sample_id`.

## Step 4: Merge by `sample_id + barcode_core`

The only safe default merge key is:

- `sample_id`
- `barcode_core`

Expected behavior:

1. standardize the scRNA metadata to canonical `sample_id`
2. derive `barcode_core` from the cell barcode
3. standardize the TCR table to the same two fields
4. merge only on `sample_id + barcode_core`

If the same `barcode_core` appears across different samples, that is expected
and is not a reason to relax the sample-aware join.

## Step 5: Replace stale TCR columns

Do not mix rebuilt TCR fields with pre-existing TCR fields of uncertain
provenance.

Required behavior:

1. remove or overwrite the old TCR columns
2. write the rebuilt canonical columns back into metadata
3. keep missing chains blank rather than forcing placeholder calls

## Canonical TCR schema

The canonical chain fields should support all four chains:

- `TRA_cdr3`, `TRA_v`, `TRA_d`, `TRA_j`, `TRA_cdr3_nt`, `TRA_clone_id`, `TRA_umis`, `TRA_reads`
- `TRB_cdr3`, `TRB_v`, `TRB_d`, `TRB_j`, `TRB_cdr3_nt`, `TRB_clone_id`, `TRB_umis`, `TRB_reads`
- `TRG_cdr3`, `TRG_v`, `TRG_d`, `TRG_j`, `TRG_cdr3_nt`, `TRG_clone_id`, `TRG_umis`, `TRG_reads`
- `TRD_cdr3`, `TRD_v`, `TRD_d`, `TRD_j`, `TRD_cdr3_nt`, `TRD_clone_id`, `TRD_umis`, `TRD_reads`

Recommended derived flags:

- `TCRseq`
- `has_TRA_TRB_paired`
- `has_TRG_TRD_paired`
- `has_any_ab_tcr`
- `has_any_gd_tcr`

## Validation expectations

Validation should happen before propagation into project milestone objects.

Minimum checks:

- join-key uniqueness on `sample_id + barcode_core`
- row counts before and after merge
- fraction with any TCR and paired-chain coverage
- chain-family plausibility
- review of suspicious near-100% paired coverage

When expression data are available, add biological validation such as:

- `TRAB` versus `TRD` score scatter
- paired-chain coverage by score region
- marker agreement for `TRDC`, `TRGC1`, `TRGV9`, `TRDV1`, `TRDV2`

## Framework-specific note

This SOP is shared across:

- Scanpy / AnnData workflows
- Seurat workflows
- supplementary 10x 5' intake workflows

Only the object-writing layer changes by framework. The TCR parsing, canonical
schema, and join-key logic do not.

## Historical warning

Earlier project logic sometimes used barcode-only or barcode-prefix-only joins.
Those approaches should be treated as legacy behavior and not reused.
