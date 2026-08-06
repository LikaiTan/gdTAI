# Extension H5AD Standalone Phase 0 QC

## Purpose

`workflows/intake/qc_extension_h5ads.py` is the review gate between extension
cohort H5AD construction and any later merge or integration. It audits each
built cohort independently and never modifies or merges source H5ADs.

A successful run reports `PASS_REVIEW_REQUIRED`. This is structural
eligibility for manual review, not approval to merge, integrate, train a model,
or build an atlas. Any structural error reports `FAIL` and returns process exit
code `2` after writing diagnostic artifacts.

## Input Contract

The default input is:

```text
data/interim/extension_intake/built_h5ads_manifest.csv
```

The manifest must have one row per built cohort and contain:

- a cohort column: `cohort_id`, `dataset_id`, or `gse_id`
- a path column: `h5ad_path`, `output_h5ad`, `output_path`, `built_h5ad`, or
  `path`
- optionally, a successful `build_status` or `status`
- optionally, an expected `n_cells`, `cells`, or `n_obs`

Relative H5AD paths are resolved first from the manifest directory and then
from the project root. Duplicate cohort IDs, reused paths, failed build status,
missing files, and manifest/H5AD cell-count disagreement fail closed.

## Structural Checks

Each H5AD is opened with AnnData backed mode (`backed="r"`) and HDF5 read-only
mode. The gate checks:

- nonempty H5AD dimensions
- unique, nonblank `obs_names` and `var_names`
- sparse CSR or CSC `X`
- sparse shape, index, and pointer-array integrity
- finite, nonnegative, integer-valued raw counts
- required harmonized metadata with no blank values:
  - `sample_id`
  - `library_id`
  - `donor_id`
  - `tissue`
  - `diagnosis`
  - `technology_simple`
  - `source_accession`
  - `barcode`
  - `barcode_core`
  - `tcr_schema_provenance`
- `technology_simple == "10x 5'"`
- globally unique `sample_id + barcode_core`
- agreement between `barcode` and normalized `barcode_core`
- one consistent sample, donor, tissue, diagnosis, accession, technology, and
  TCR provenance per `library_id`
- one `library_id` per `sample_id`, so sample-level TCR join summaries cannot be
  duplicated or ambiguously allocated in library-level reports

Sparse `X` is streamed from its `data`, `indices`, and `indptr` arrays in
bounded blocks. The workflow never calls `toarray`, never constructs a dense
cell-by-gene matrix, and retains only per-cell QC vectors.

## Canonical TCR Contract

All four chains (`TRA`, `TRB`, `TRG`, and `TRD`) must have the following fields:

```text
*_cdr3, *_v, *_d, *_j, *_cdr3_nt, *_clone_id, *_umis, *_reads
```

The following boolean fields are also required:

```text
has_TRA, has_TRB, has_TRG, has_TRD,
has_TRA_TRB_paired, has_TRG_TRD_paired,
has_any_ab_tcr, has_any_gd_tcr
```

The gate enforces these identities cell by cell:

```text
has_CHAIN          == nonempty CHAIN_cdr3
has_TRA_TRB_paired == has_TRA AND has_TRB
has_TRG_TRD_paired == has_TRG AND has_TRD
has_any_ab_tcr     == has_TRA OR has_TRB
has_any_gd_tcr     == has_TRG OR has_TRD
```

Chain details and positive read/UMI counts cannot exist without a productive
chain CDR3. Numeric TCR fields must be finite nonnegative integers. When
`TCRseq` exists, it must be `yes` exactly when any canonical chain is present.

For separately joined productive contigs, `tcr_schema_provenance` triggers a
required `uns["tcr_join_summaries"]` audit. The recorded expression-cell,
matched-cell, unmatched-TCR, and join-coverage values must agree with observed
`sample_id` counts and canonical TCR-positive cells. Embedded TCR schemas are
reported as `embedded_or_not_applicable`; productive TCR yield is still
reported for every cohort and library.

## QC Review Metrics

The workflow does not filter cells. It reports, per cohort and per library:

- cells, samples, donors, genes, sparse nonzeros, and count range
- median and 5th/95th percentile total counts and detected genes
- median and upper-tail mitochondrial fraction
- zero-count cells, cells below 200 detected genes, and cells above 20% MT
- cells with TRA, TRB, TRG, TRD, paired TRA/TRB, paired TRG/TRD, or any TCR
- corresponding productive-TCR and paired-chain rates
- sample-aware TCR join counts and coverage when applicable

Low-quality or low-yield observations are warnings for human review. They do
not cause automatic removal. Schema, matrix, key, or logical failures are
structural errors and block the gate.

## Run Commands

Use the canonical environment:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/intake/qc_extension_h5ads.py --dry-run
```

The dry run validates the manifest and input paths, prints the planned cohort
list, and writes nothing. Run the complete audit with:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/intake/qc_extension_h5ads.py
```

Alternative paths are explicit:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/intake/qc_extension_h5ads.py \
  --manifest /path/to/built_h5ads_manifest.csv \
  --output-root /path/to/Integrated_dataset
```

`--output-root` must be the directory that will contain `tables/`, `logs/`, and
`figures/`.

## Review Artifacts

Tables:

```text
Integrated_dataset/tables/extension_intake/extension_phase0_manifest_audit.csv
Integrated_dataset/tables/extension_intake/extension_phase0_cohort_summary.csv
Integrated_dataset/tables/extension_intake/extension_phase0_library_summary.csv
Integrated_dataset/tables/extension_intake/extension_phase0_validation_issues.csv
```

Logs:

```text
Integrated_dataset/logs/extension_intake/extension_phase0_qc_summary.json
Integrated_dataset/logs/extension_intake/extension_phase0_qc_summary.md
```

Figures:

```text
Integrated_dataset/figures/extension_intake/extension_phase0_cohort_overview.png
Integrated_dataset/figures/extension_intake/extension_phase0_library_qc.png
```

The JSON is the machine-readable gate record. It always includes
`review_required: true`, `merge_approved: false`, the manifest checksum, source
paths and file statistics, artifact paths, cohort summaries, and all issues.
Only explicit user review can authorize a later phase.
