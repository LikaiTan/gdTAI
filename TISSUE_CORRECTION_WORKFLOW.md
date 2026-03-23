# Tissue Correction Workflow

## Purpose

This workflow creates a new `tissue_corrected` column for the full integrated
milestone H5AD without requiring any further manual GEO review at execution
time.

It is designed for the large integrated object:

- `high_speed_temp/Integrated_dataset/integrated.h5ad`

The workflow is split into two bounded steps:

1. `audit`
   - sample up to 500 cells per GSE
   - review tissue-related columns and the proposed corrected tissue labels
   - identify which datasets require GEO series-matrix lookups
2. `apply`
   - load the full integrated `obs`
   - align per-cell harmonized metadata from:
     - `analysis_26GSE_V4/outputs/harmonized_metadata_v4.csv`
     - `analysis_26GSE_V4/outputs/harmonized_metadata_supp.csv`
   - apply the committed rule set
   - write `obs["tissue_corrected"]` back into the integrated H5AD in place

## Canonical Files

- Rule config: [tissue_correction_rules.json](/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/tissue_correction_rules.json)
- Workflow script: [tissue_correction_workflow.py](/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/tissue_correction_workflow.py)

## Tissue Scope

`tissue_corrected` is a coarse harmonized tissue or specimen-site field.

Current coarse values intentionally keep some biologically distinct sources
separate:

- `blood`
- `balf`
- `csf`
- `tumor`
- `adjacent_normal_tissue`
- `lymph_node_metastasis`
- `bone_marrow`
- `heart`
- `duodenum`
- `mixed_duodenum_pbmc`
- `rectum`
- `ileum`
- `colon`
- `tracheal_aspirate`
- `liver_tumor`
- `adjacent_normal_liver`
- `biopsy_other`
- `unknown`

Notes:

- `PBMC` is collapsed into `blood`.
- `tumor` and `adjacent_normal_tissue` stay separate.
- `unknown` is used instead of forcing unsupported inference.

## What The Rule Set Uses

The script resolves tissue from three information layers:

1. Integrated `obs` columns
   - `tissue`
   - `library_id`
   - `sample_id`
   - `sampleid`
   - `source_gse_id`
   - `metadata_key`
2. Harmonized metadata
   - `tissue`
   - `sample_type`
   - `library_id`
   - `sample_id`
   - `sampleid`
3. GEO `series_matrix` files for configured difficult datasets
   - used only when the tissue cannot be trusted directly from `obs`
   - the workflow parses GSM accession, sample title, source name, and
     `!Sample_characteristics_ch1`

## Datasets With Explicit GEO-Backed Rules

The rule set currently contains explicit GEO-backed or dataset-specific logic
for:

- `GSE125527`
- `GSE145926`
- `GSE162498`
- `GSE171037`
- `GSE178882`
- `GSE190870`
- `GSE228597`
- `GSE240361`
- `GSE243905`
- `GSE252762`
- `GSE301528`

Examples:

- `GSE125527`
  - `R -> rectum`
  - `I -> ileum`
  - `PBMC -> blood`
- `GSE171037`
  - `TA -> tracheal_aspirate`
  - `PBMC -> blood`
- `GSE190870`
  - `Primary tumor -> tumor`
  - `Lymph node metastasis -> lymph_node_metastasis`
- `GSE228597`
  - `Heart tissue -> heart`
  - `Normal tissue adjacent to tumor -> adjacent_normal_tissue`
  - `tumor -> tumor`

## Output Package

Both subcommands write tables under:

- `Integrated_dataset/tables/tissue_correction/`

Expected outputs include:

- downsampled audit CSV and markdown
- full-run `tissue_corrected` value counts
- per-GSE validation table
- unresolved examples table
- apply summary markdown

## Execution

Audit first:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  tissue_correction_workflow.py audit --verbose
```

Apply to the full integrated object and export the standalone review table:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  tissue_correction_workflow.py apply --verbose
```

Write back into the H5AD only after review approval:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  tissue_correction_workflow.py apply --write-h5ad --verbose
```

## Extension Rule

When a new integrated dataset has unresolved tissue values:

1. rerun `audit`
2. inspect the downsampled unresolved rows
3. check the dataset `series_matrix` only if the harmonized metadata still
   cannot resolve tissue
4. add the new deterministic mapping to
   [tissue_correction_rules.json](/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/tissue_correction_rules.json)
5. rerun `apply`

Do not hard-code conclusions in a one-off notebook or assistant reply.
Commit the new rule so the full apply step remains deterministic.
