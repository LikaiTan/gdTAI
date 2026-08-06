# Sample Source Refinement

## Scope

This workflow creates one new composite metadata field:

- `obs["sample_source_refined"]`

It does not replace or normalize any original metadata field. Provenance level,
rule ID, reason, and source values remain in the review exports rather than
adding more columns to the H5AD.

For tumor-related projects, the refined value must retain both disease and
specimen context. Documented primary tumor, adjacent non-tumor tissue, blood,
and metastatic sites are therefore distinct values. A diagnosis field is not
treated as specimen context: a cohort-wide cancer diagnosis must never turn a
matched blood sample into a tumor sample.

Rules are stored in:

- `configs/metadata/sample_source_refinement_rules.json`

The executable workflow is:

- `workflows/metadata/sample_source_refinement_workflow.py`

## Precedence

Resolution order is fixed and validated at runtime:

1. explicit cell-level source or cell override
2. exact sample/library rule
3. GSE plus tissue rule
4. global tissue fallback
5. `unresolved`

GSE287301 patient 10 was ultimately diagnosed with ameloblastoma and is present
in pooled libraries with HNSCC cases. The tissue context is known, so these
rows are retained as `tumor_pool_HNSCC_with_ameloblastoma_ambiguity`; the label
does not claim a per-cell histology that cannot be recovered from pooled data.

## Controlled Mappings

The committed rules include:

- GSE162498 tumor rows: `NSCLC_tumor`
- GSE243013 exact sample map: `LUAD_tumor` or `LUSC_tumor`
- GSE287301 tumor rows outside the patient-10 pools: `HNSCC_tumor`
- GSE235863 liver-tumor rows: `HBV_positive_HCC_tumor`
- GSE190870 metastatic rows:
  `breast_invasive_ductal_carcinoma_lymph_node_metastasis`
- extension tumor cohorts with explicit source distinctions:
  GSE114724, GSE121636, GSE159251, GSE292700, GSE294273,
  GSE296954, and GSE315928

The GSE+tissue rules are tissue-gated. Blood and adjacent-normal rows from a
mixed accession retain their tissue source through the global fallback.

## Review First

Use the canonical environment. The dry run opens the H5AD read-only and
streams bounded `obs` chunks without loading `X`:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/metadata/sample_source_refinement_workflow.py review \
  --input-h5ad high_speed_temp/Integrated_dataset/integrated.h5ad \
  --rules configs/metadata/sample_source_refinement_rules.json \
  --output-dir Integrated_dataset/tables/sample_source_refinement \
  --chunk-size 100000
```

`dry-run` is an alias for `review`.

Review outputs include:

- `sample_source_refined_review.csv.gz`: every proposed value and provenance
- `sample_source_refined_by_gse.csv`: per-GSE resolution rate
- `sample_source_refined_value_counts.csv`: value/rule counts
- `sample_source_refined_tumor_project_audit.csv`: tumor-project gate counts
- `sample_source_refined_unresolved_examples.csv.gz`: bounded unresolved rows
- `sample_source_refinement_review_manifest.json`: input, rule, and artifact
  checksums required for writeback

Writeback is blocked if a configured tumor-project row resolves to a generic
or unsupported context. The original `tissue`, `diagnosis`, and source fields
are never overwritten.

## Extension Cohorts

After standalone extension H5ADs have been built, run the batch review:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/metadata/review_extension_sample_sources.py
```

This is read-only and writes per-cohort plus aggregate review artifacts under
`Integrated_dataset/tables/sample_source_refinement/extensions/` and a summary
under `Integrated_dataset/logs/sample_source_refinement/`.

## Writeback

After reviewing those artifacts, write the single new column with:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/metadata/sample_source_refinement_workflow.py writeback \
  --input-h5ad high_speed_temp/Integrated_dataset/integrated.h5ad \
  --rules configs/metadata/sample_source_refinement_rules.json \
  --review-manifest Integrated_dataset/tables/sample_source_refinement/sample_source_refinement_review_manifest.json \
  --chunk-size 100000 \
  --confirm-reviewed
```

Writeback fails if the H5AD, rules, resolved source columns, proposed values,
or review artifacts differ from the reviewed manifest. Values are written to a
staging dataset, read back, and committed only after checksum validation. The
workflow also verifies that HDF5 structure outside the owned output column is
unchanged. Replacing an existing output column additionally requires
`--replace-existing` and a fresh review manifest for the current H5AD.
