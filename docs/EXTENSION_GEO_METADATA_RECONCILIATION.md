# Extension GEO Metadata Reconciliation

## Purpose

This package reconciles only unresolved or semantically ambiguous metadata for
the extension cohorts. It does not replace original metadata, apply derived
values, or write an H5AD.

The reviewed scope is:

- `GSE114724`
- `GSE121636` / `GSE121637`
- `GSE159251`
- `GSE169246`
- `GSE292700`
- `GSE294273` / `GSE294274`
- `GSE296954`
- `GSE315928`

## Method

1. Existing extension Phase 0, T/NK-filter, and sample-source audit outputs
   were inspected first. Structural completeness was kept separate from
   evidentiary resolution: a nonblank local value can still be unsupported or
   over-specific.
2. Only official NCBI GEO series pages, GSM pages, series matrices, family
   SOFT, and GEO-hosted supplementary metadata were used for final claims.
3. Reconciliation rows preserve the exact local value and add a separate
   GEO-supported value. No proposed value is written back.
4. Repeated samples with the same issue are represented as a bounded,
   homogeneous scope. Machine-checkable rows include the local column and an
   optional `library_id` regular expression.
5. The audit script opens each filtered H5AD with `backed="r"`, validates the
   declared local values, records file size and nanosecond mtime before and
   after, and closes the file without mutation.
6. Source-accession coverage is checked against the cohort config, library
   config, staged GEO metasheets, reconciliation CSV, and cell-level
   provenance. `GSE121637` and `GSE294274` are expected in the TCR source
   config and staged records but not in cell-level provenance, which follows
   the paired GEX accession.

The deterministic input is
`configs/metadata/extension_geo_metadata_reconciliation.csv`.

## Status Semantics

- `resolved_verified_local`: the local value agrees with official GEO.
- `resolved_geo_value_local_unresolved`: GEO provides a deterministic value or
  relationship, but the local field is absent, collapsed, or semantically
  incorrect. The CSV records the correction without applying it.
- `ambiguous_geo_partial`: GEO narrows the interpretation but does not support
  a fully specific value.
- `unavailable_in_geo`: the reviewed official GEO records do not provide the
  requested value. This is distinct from a locally unresolved field for which
  GEO does provide an answer.

## Results

The reconciliation contains 58 grouped field rows across 8 cohorts and 10
official GEO accessions.

| Cohort | Rows | Verified/resolved | GEO-partial ambiguous | GEO-unavailable |
|---|---:|---:|---:|---:|
| GSE114724 | 5 | 4 | 1 | 0 |
| GSE121636_GSE121637 | 8 | 7 | 1 | 0 |
| GSE159251 | 7 | 7 | 0 | 0 |
| GSE169246 | 10 | 9 | 1 | 0 |
| GSE292700 | 6 | 4 | 1 | 1 |
| GSE294273_GSE294274 | 8 | 7 | 0 | 1 |
| GSE296954 | 7 | 5 | 1 | 1 |
| GSE315928 | 7 | 5 | 1 | 1 |
| **Total** | **58** | **48** | **6** | **4** |

The 48 resolved rows comprise 16 verified local values and 32 cases where GEO
supplies a value or relationship that is missing, collapsed, or over-broad
locally.

### Material Findings

- `GSE114724` contains five selected breast-tumor GEX libraries with five
  matching TCR GSMs. The broader parent study's normal breast, blood, and lymph
  node specimens do not describe these selected libraries. Primary versus
  metastatic origin is not deposited for the five tumors.
- `GSE121636` and `GSE121637` are paired GEX and VDJ subseries. Three GU donors
  each contribute renal tumor and matched blood, with exact title-based GEX/TCR
  pairing. GEO does not explicitly classify the renal tumors as primary versus
  metastatic.
- `GSE159251` contains blood, involved lymph-node metastasis, and axillary
  subcutaneous metastasis samples. Three VDJ-labelled GSMs deposit RDS objects
  containing both an RNA assay and TCR metadata; they have no separate GEX GSM.
- `GSE169246` has the main local semantic error. Official titles ending `_b`
  are blood, but the current patient/timepoint join propagated tumor-site or
  blank tissue values and tumor/other specimen contexts to those libraries.
  This affects all 239,418 retained `_b` cells: local tissue is blank for
  111,432 cells and a tumor site for 127,986 cells; local specimen context is
  `other_tissue` for 111,432, `primary_tumor` for 48,867, and
  `metastatic_tumor` for 79,119. All require a separate GEO-derived blood
  value before merge.
  Official titles preserve `_b` versus `_t`; a derived specimen key must retain
  that compartment. The 78 non-ATAC RNA/VDJ records are distinct from five
  ATAC-Seq GSMs. Breast and chest-wall tumor records remain ambiguous for exact
  primary, recurrent, locoregional, or metastatic classification.
- `GSE292700` pairs tumor, paracancerous adjacent tissue, and blood from two
  lung-adenocarcinoma donors with title-matched TCR GSMs. Adjacent and blood
  records are cancer-patient specimens, not healthy controls. GEO does not
  state whether each tumor is primary or metastatic.
- `GSE294273` and `GSE294274` are one-to-one GEX/TCR pairs by HM identifier from
  metastatic melanoma lymph-node samples. GEO gives cohort-level stage III/IV
  and untreated/resistant/responder groups, but not per-patient stage,
  demographics, or a specific ICI regimen.
- `GSE296954` has eleven patients with paired A/B tumor biopsies and matching
  GEX/VDJ records. GEO supports HPV-unrelated HNSCC and treatment assignment,
  but does not explicitly locate each biopsy as primary, nodal, or metastatic
  and does not deposit per-patient response, demographics, or stage.
- `GSE315928` contains 94 primary gastric-tumor sample GSMs from 33 patients;
  each deposited H5AD represents paired scRNA/TCR data. GEO describes
  pre-treatment, post-chemotherapy, and post-immunotherapy stages, but the
  series matrix alone does not explicitly map `B`, `F1`, and `F2` suffixes to
  those stages. The GEO-hosted workbook's tumor-grade field is blank.

## Remaining Items

Six rows remain `ambiguous_geo_partial`:

- exact primary/metastatic class for `GSE114724` tumors
- exact primary/metastatic class for `GSE121636` renal tumors
- exact primary/recurrent/locoregional/metastatic class for `GSE169246` breast
  and chest-wall tumor records
- exact primary/metastatic class for `GSE292700` tumors
- exact primary/nodal/metastatic biopsy site for `GSE296954`
- exact `B`/`F1`/`F2` suffix mapping for `GSE315928`

Four grouped rows are `unavailable_in_geo`:

- `GSE292700`: patient demographics, stage, treatment, and response
- `GSE294273`: per-patient stage, demographics, and specific ICI regimen
- `GSE296954`: per-patient response, demographics, and stage
- `GSE315928`: tumor grade

These counts describe reconciliation rows, not cells and not the number of
individual scalar variables within a grouped unavailable item.

## Validation And Outputs

Run with the canonical environment:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/metadata/audit_extension_geo_metadata.py
```

Generated tables are under
`Integrated_dataset/tables/geo_metadata_reconciliation/` and the audit log is
under `Integrated_dataset/logs/geo_metadata_reconciliation/`.

The audit emits:

- validated reconciliation rows
- status counts by cohort
- official source-accession coverage
- local H5AD metadata value counts
- reconciliation-scope local value counts, including affected cell counts
- unresolved/ambiguous rows
- H5AD read-only file-state validation
- an empty validation-issues table on success
- a deterministic audit manifest and concise Markdown log

## Limitations

- This is a bounded official-GEO reconciliation, not a publication-wide
  clinical abstraction. Publication-only values are not substituted for absent
  GEO fields.
- Direct URLs and concise evidence paraphrases are committed in the CSV. The
  deterministic audit validates URL form but does not require network access or
  re-download GEO files at runtime.
- Grouped scopes are used only where the official and local values are
  homogeneous under the stated selector. They are not cell-level annotations.
- A `resolved_geo_value_local_unresolved` row is evidence for a future derived
  field, not approval to overwrite the original field.
- Audit execution performs no merge, integration, model/filter change, or H5AD
  write. Any future correction must use additive derived columns and pass a
  separate review gate.
