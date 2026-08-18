# Full-Atlas Metadata Harmonization v2 Review

## Status

**Review complete; rules not applied.** This package proposes additive metadata
columns for the 5,933,312-cell rebuilt atlas. The H5AD was opened with HDF5 mode
`r` only and its resolved path, size, modification time, and inode were identical
before and after the audit.

The proposed source-name aliases are also additive:

| Existing `source_gse_id` | Proposed `source_accession_harmonized_v2` |
|---|---|
| `GDT_2020AUG_woCOV` | `GSE149356` |
| `MalteGDT` | `PRJNA771240` |
| `BALF_BLOOD_COPD` | `unpublished_balf_blood` |
| `GDTlung2023july_7p` | `unpublished_gdTlung` |

## Why The Current Labels Duplicate Or Fail

`workflows/integration/rebuild_full_atlas.py::standardize_obs` copied the first
nonempty value from a list of aliases. It did not canonicalize capitalization,
spaces, or underscores, so values such as `Tumor`/`tumor`,
`PBMC`/`peripheral blood`/`blood`, and `lymph node`/`lymph_node` remained
distinct.

The same helper populated `specimen_context` from
`specimen_context`, `sample_tissue_status`, `diagnosis_tissue`, and finally the
generic `type` field. In some source objects `type` is a cell-type/state label.
Consequently, the rebuilt atlas contains 460,269 cells whose current
`specimen_context` is `CD4`, `CD8`, `NK`, `Treg`, `Naive`, or `Proliferating`.
Those values are rejected as specimen evidence.

The rebuilt object does not carry the legacy `tissue_corrected` or
`sample_source_refined` columns. The v2 proposal therefore does not assume that
old corrections were propagated. It uses only source accession, sample/library
identifiers, original cell provenance, `tissue_original`, and the current raw-like
`tissue_harmonized` value.

## Proposed Additive Schema

- `source_accession_harmonized_v2`: user-approved canonical source alias while
  preserving `source_gse_id`.
- `tissue_harmonized_v2`: anatomical or collection compartment, using a bounded
  snake-case vocabulary (`blood`, `lung`, `lymph_node`,
  `bronchoalveolar_lavage_fluid`, and so on).
- `specimen_context_harmonized_v2`: specimen role (`blood`, `primary_tumor`,
  `metastatic_tumor`, `tumor_unspecified`, `adjacent_non_tumor`,
  `non_tumor_tissue`, and so on).
- `tumor_type_harmonized_v2`: cohort cancer diagnosis independent of specimen
  compartment. Thus melanoma PBMC remains `tissue=blood`, `context=blood`, and
  `tumor_type=melanoma`.
- `sample_id_harmonized_v2`, `library_id_harmonized_v2`, and
  `donor_id_harmonized_v2`: additive biological/technical identity fields;
  current identifiers remain untouched.
- `tcr_library_join_id_v2`: the validated technical VDJ join identity. This is
  deliberately separate from biological sample identity because pooled TCR
  libraries can contain multiple specimens.
- `metadata_rule_id_v2`, `metadata_evidence_url_v2`,
  `metadata_evidence_level_v2`, and `metadata_status_v2`: provenance and
  resolution state.

Precedence is source-scoped sample/library rule, imported exact sample map,
source default, normalized rebuilt tissue synonym, then fail-closed
`unresolved`. The existing `specimen_context`, `condition`, `diagnosis`, and
cell-type columns are explicitly forbidden as rule inputs.

## Audit Results

All 11 validation checks passed:

- exact atlas dimensions audited: 5,933,312 cells and 40 sources;
- source counts in the H5AD match the rebuild report;
- every source has a bounded source-default rule;
- scoped rules use only permitted rebuilt provenance fields;
- no scoped rule uses a forbidden/polluted field;
- no proposed output literal equals a cell-type/state exclusion term;
- all 24 cancer-related evidence records use direct official NCBI GEO URLs;
- the 22,276,972,340-byte H5AD retained the same mtime and inode.

The current atlas has 41 `tissue_harmonized` categories and 3,849,174 blank
tissue values. The v2 review intentionally does not claim those blanks are all
resolvable. GSE125527 is recoverable despite blank atlas fields: all 46,396
atlas rows match its source H5AD by `original_cell_id`, allowing 30 donor-by-
tissue samples and blood-versus-rectum tissue to be restored.

## Sample And Library Identity Gate

`workflows/metadata/audit_full_atlas_sample_identity_v2.py` passed all 13
checks without applying TCR chain calls or modifying an H5AD:

- all 2,155,409 atlas rows from the 14 TCR-repaired sources have an exact
  technical VDJ-library identity;
- 46,396 GSE125527 cells recover 30 sample/library values directly from source
  donor and tissue fields;
- 406,494 GSE254249 cells recover 92 sample/library values from source `Ident`,
  replacing the collapsed display value `gse254249_scrna_matrix`;
- all 4,611 blank GSE228597 library IDs recover a technical VDJ-library ID, but
  their biological sample remains explicitly `unresolved` because the pooled
  TCR library is not a biological specimen;
- no blank proposed sample or library value remains; 4,611 biological samples
  are explicitly unresolved rather than guessed;
- donor IDs are preserved because technical VDJ pools must not be substituted
  for biological donors.

This identity correction is a metadata step. TRA/TRB/TRG/TRD calls remain
unapplied until the tissue and identity columns pass their write gate.

## Cancer Evidence Coverage

| Source | Supported tumor type or fail-closed label | Specimen evidence | Official evidence |
|---|---|---|---|
| GSE114724 | breast carcinoma | selected breast tumors; primary/metastatic unresolved | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114724 |
| GSE121636 | clear-cell renal cell carcinoma | renal tumor and matched blood | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121636 |
| GSE144469 | unresolved | colon biopsies; source cancers not deposited | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144469 |
| GSE159251 | metastatic melanoma | blood and metastatic lesions | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159251 |
| GSE162498 | non-small-cell lung cancer | tumor, adjacent/juxta, and blood | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162498 |
| GSE169246 | triple-negative breast cancer | `_b` blood and `_t` tumor | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169246 |
| GSE178882 | unresolved mixed metastatic solid tumor | PBMC | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178882 |
| GSE188620 | low-grade glioma | atlas scRNA cells are PBMC | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188620 |
| GSE190870 | breast invasive ductal carcinoma | primary tumor and lymph-node metastasis | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190870 |
| GSE206325 | hepatocellular carcinoma | tumor and non-involved liver | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206325 |
| GSE211504 | melanoma | whole blood | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211504 |
| GSE212217 | endometrial carcinoma | PBMC | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212217 |
| GSE228597 | unresolved cancer type | heart, PBMC, tumor, and adjacent tissue | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228597 |
| GSE235863 | HBV-positive hepatocellular carcinoma | blood, liver tumor, and adjacent liver | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235863 |
| GSE243013 | sample-level LUAD or LUSC | post-treatment NSCLC tumor | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243013 |
| GSE243572 | hepatocellular carcinoma | PBMC | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243572 |
| GSE254249 | locally advanced rectal cancer | blood and rectal tumor tissue | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254249 |
| GSE287301 | HNSCC with patient-10 ambiguity | tumor | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287301 |
| GSE287541 | unresolved cancer type | PBMC | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287541 |
| GSE292700 | lung adenocarcinoma | tumor, adjacent tissue, and blood | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292700 |
| GSE294273 | metastatic melanoma | lymph-node metastasis | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294273 |
| GSE296954 | HPV-unrelated HNSCC | tumor biopsy; exact site class unresolved | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296954 |
| GSE311112 | chronic lymphocytic leukemia | peripheral blood | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311112 |
| GSE315928 | gastric cancer | primary gastric tumor | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315928 |

## Remaining Unresolved Mappings

The authoritative row-level list is
`Integrated_dataset/tables/metadata_harmonization/full_atlas_v2/unresolved_mappings.csv`.
The main unresolved classes are:

- provenance loss: 36,800 blank-tissue GSE206325 cells and 4,611 GSE228597
  cells with unresolved biological sample identity;
- unavailable individual cancer histology: GSE144469, GSE178882, GSE228597,
  and GSE287541;
- biological ambiguity: patient-10-containing GSE287301 pools;
- unsupported primary-versus-metastatic detail: GSE114724, GSE121636,
  GSE169246, GSE254249, GSE292700, and GSE296954;
- unpublished/local disease context: `unpublished_gdTlung` and PRJNA771240.

No value in these classes should be inferred from cell phenotype, clustering,
or the polluted current `specimen_context`.

## Files And Next Gate

The rules are in
`configs/metadata/full_atlas_metadata_harmonization_v2.json`; the reproducible
read-only validator is
`workflows/metadata/audit_full_atlas_metadata_harmonization_v2.py`; audit tables
are under
`Integrated_dataset/tables/metadata_harmonization/full_atlas_v2/`.

No H5AD was changed. Applying the columns remains a separate metadata-write
gate requiring backup, exact row-alignment checks, and explicit approval.

## Full Row-Level Dry Run

`workflows/metadata/build_full_atlas_metadata_overlay_v2.py` built the complete
additive overlay at
`Integrated_dataset/tables/metadata_harmonization/full_atlas_v2/metadata_overlay_v2.parquet`.
It contains exactly 5,933,312 rows and 40 sources, is 34,622,134 bytes, and has
SHA-256
`4da4eea32e9d275de790775e2a0f59d4f6553d72756e9c8c6935f35bb398984f`.

All 18 dry-run checks passed:

- the atlas, two source identity H5ADs, prior 3,705,306-row refinement review,
  and TCR sidecar were unchanged and checksum-bound;
- every overlay row is keyed by source and original cell identity;
- all historical-core cells matched the prior reviewed metadata exactly;
- all 2,155,409 cells from TCR-repaired sources have a separate technical
  VDJ-library identity;
- no sample or library field is blank, and the only explicitly unresolved
  biological sample IDs are the 4,611 GSE228597 cells described above;
- no CD4/CD8/NK/Treg/Naive/Proliferating value enters tissue or specimen
  context;
- tissue, context, and tumor type all use their controlled vocabularies.

Compared with 3,849,174 blank current tissue values, only 39,920 rows remain
fail-closed in `tissue_harmonized_v2` and
`specimen_context_harmonized_v2`:

- 36,800 GSE206325 rows have blank tissue and donor provenance in the source
  H5AD itself;
- 3,120 GSE252762 rows originate from a mixed duodenum/PBMC library whose
  per-cell specimen compartment cannot be recovered from deposited metadata.

The next operation is a high-risk metadata-only candidate-H5AD write. It must
retain the current atlas as rollback, write and validate a separate candidate,
and must not apply any TRA/TRB/TRG/TRD values. TCR chain propagation remains a
subsequent transaction after this metadata candidate passes.

The exact write boundary and rollback contract are frozen in
`docs/FULL_ATLAS_METADATA_V2_WRITE_PLAN.md`.
