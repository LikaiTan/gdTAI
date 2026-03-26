# Disease Status Correction Workflow

## Purpose

This workflow creates a new `disease_status_corrected` column for the full
integrated milestone H5AD without requiring any further manual lookup at apply
time.

It is designed for the large integrated object:

- `high_speed_temp/Integrated_dataset/integrated.h5ad`

The workflow is split into two bounded steps:

1. `audit`
   - sample up to 500 cells per GSE
   - review disease-related columns and the proposed corrected labels
   - identify which datasets still need stronger GEO-backed rules
2. `apply`
   - load the full integrated `obs`
   - align per-cell harmonized metadata from:
     - `analysis_26GSE_V4/outputs/harmonized_metadata_v4.csv`
     - `analysis_26GSE_V4/outputs/harmonized_metadata_supp.csv`
   - apply the committed rule set
   - export `disease_status_corrected` as a standalone review table
   - write back into `integrated.h5ad` only when explicitly requested

## Canonical Files

- Rule config: [disease_status_correction_rules.json](/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/disease_status_correction_rules.json)
- Workflow script: [disease_status_correction_workflow.py](/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/disease_status_correction_workflow.py)

## Disease Scope

`disease_status_corrected` is a coarse harmonized status field:

- `healthy`
- `disease`
- `unknown`

This field is intentionally simpler than the raw condition vocabulary. It is
for broad case-vs-control interpretation, not detailed diagnosis modeling.

Examples:

- `control`, `healthy`, `normal control`, `COVID19-Negative` -> `healthy`
- `COVID-19`, `myocarditis`, `colitis`, `glioma`, `HNSCC`, `rrLBCL` -> `disease`
- `pre-treatment`, `on-treatment`, `post-treatment` -> `disease`
  because those labels imply a diseased cohort even when they are not
  themselves diagnoses

`unknown` is used instead of forcing unsupported inference.

## What The Rule Set Uses

The script resolves disease status from three information layers:

1. Harmonized metadata
   - `condition`
   - `sample_type`
   - `donor_patient`
   - `project name`
   - `library_id`
   - `sample_id`
   - `sampleid`
2. GEO `series_matrix` files
   - auto-discovered under `downloads/GSE*/matrix/*series_matrix.txt.gz`
   - used only when metadata wording is incomplete or misleading
3. Dataset defaults
   - used for clearly disease-only cohorts when the per-sample wording is
     treatment-state or otherwise non-diagnostic

## Current Rule Design

The current rule set contains:

- global healthy/control wording
- global disease wording
- global treatment-state wording
- explicit per-GSE overrides for mixed cohorts or abbreviations such as:
  - `GSE125527`
  - `GSE144469`
  - `GSE161918`
  - `GSE228597`
  - `GSE240865`
  - `GSE252762`
  - `GSE287541`
- dataset-level defaults for clearly disease-only cohorts such as:
  - `GSE171037`
  - `GSE179994`
  - `GSE221776`
  - `GSE232240`
  - `GSE235863`
  - `GSE241783`
  - `GSE287301`

## Output Package

Both subcommands write tables under:

- `Integrated_dataset/tables/disease_status_correction/`

Expected outputs include:

- downsampled audit CSV and markdown
- full-run `disease_status_corrected` value counts
- per-GSE validation table
- unresolved examples table
- apply summary markdown

## Execution

Audit first:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  disease_status_correction_workflow.py audit --verbose
```

Apply to the full integrated object and export the standalone review table:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  disease_status_correction_workflow.py apply --verbose
```

Write back into the H5AD only after review approval:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  disease_status_correction_workflow.py apply --write-h5ad --verbose
```

## Extension Rule

When a new integrated dataset still has unresolved disease status:

1. rerun `audit`
2. inspect the downsampled unresolved rows
3. check the dataset `series_matrix` only if harmonized metadata still cannot
   resolve the status
4. add the new deterministic mapping to
   [disease_status_correction_rules.json](/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/disease_status_correction_rules.json)
5. rerun `apply`

Do not leave the logic only in an assistant reply. Commit the rule so the full
apply step remains deterministic.
