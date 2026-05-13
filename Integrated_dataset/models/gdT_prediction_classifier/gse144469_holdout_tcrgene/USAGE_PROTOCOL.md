# Selected gdT Classifier Usage Protocol

This protocol explains how to apply the selected gdT classifier to another
local H5AD file. It is intended for another analyst or agent using the same
server or a clone of this GitHub repository.

## Files

Model artifact:

```text
Integrated_dataset/models/gdT_prediction_classifier/gse144469_holdout_tcrgene/selected_model.pkl
```

Inference wrapper:

```text
predict_with_selected_gdt_model.py
```

Example scripts:

```text
Integrated_dataset/models/gdT_prediction_classifier/gse144469_holdout_tcrgene/examples/apply_model_to_h5ad.sh
Integrated_dataset/models/gdT_prediction_classifier/gse144469_holdout_tcrgene/examples/minimal_python_prediction.py
Integrated_dataset/models/gdT_prediction_classifier/gse144469_holdout_tcrgene/examples/audit_false_positive_predictions.py
```

Model checksum:

```text
SHA256 ef7dc0cca09acfcee27e059782e26688f79538c53ba0621025c1ad6e2b24383d
```

Security note: `selected_model.pkl` is a Python pickle. Load only this trusted
file from this repository. Do not load arbitrary external pickle files.

## What The Model Uses

The selected model is the high-specificity GSE144469-holdout TCR-gene
classifier.

Inputs used by the model:

- individual TCRA, TCRB, TCRG, and TCRD gene expression features
- `FOXP3`
- `CD4`
- `CD3D`, `CD3E`, and `CD3G`

Inputs not used by the model:

- Phase 4 `TRD_score`
- Phase 4 `TRAB_score`
- `TRD_minus_TRAB`
- TCR CDR3 metadata
- productive TCR metadata labels
- dataset ID, tissue, sample ID, or cluster labels

The model rebuilds per-cell `log1p(counts per 10,000)` features from the H5AD
`X` matrix:

```text
feature_value = log1p(raw_or_count_like_value * 10000 / total_cell_count)
```

After model scoring, the saved inference logic applies the same death-penalty
rule used during validation:

```text
if FOXP3 > 0.25, or CD4 > 0.75, or CD3D + CD3E + CD3G < 0.25:
    adjusted_score = min(raw_model_score, 0.03)
```

The saved gdT-positive threshold is:

```text
0.9769781827926636
```

## Input H5AD Requirements

Required:

- H5AD file readable by `h5py`
- gene symbols in `var/_index`
- count-like expression in `X`
- sparse CSR `X` is preferred; dense `X` is also supported by the wrapper
- human TCR gene symbols matching the selected model gene names

Strongly recommended:

- raw counts in `X`
- sample/library metadata in `obs`, such as `source_gse_id`, `sample_id`, or
  `library_id`
- annotation metadata if you want to audit NK-cell false positives
- TCR metadata if you want to audit TCRAB-positive false positives

Important caveat:

If `X` is already normalized or log-transformed, the wrapper will still run, but
the feature scale will not match training. In that case, report the matrix state
clearly and interpret predictions cautiously.

## Quick Start On The Same Server

From the repository root:

```bash
cd /home/tanlikai/databank/publicdata/tools/output_geo_tcell_research

/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  predict_with_selected_gdt_model.py \
  --input-h5ad /path/to/test_dataset.h5ad \
  --model-pkl Integrated_dataset/models/gdT_prediction_classifier/gse144469_holdout_tcrgene/selected_model.pkl \
  --dataset-id TEST_DATASET \
  --chunk-size 50000
```

The input H5AD is opened read-only and is not modified.

Useful optional flags:

```text
--obs-column COLUMN_NAME       copy an extra obs column into the output table
--no-figures                  skip PNG generation
--max-umap-points 300000      cap UMAP plotting sample size
--chunk-size 100000           increase chunk size for larger memory allowance
```

The wrapper automatically copies these obs columns if present:

```text
source_gse_id, sample_id, sampleid, library_id, tissue_corrected, tissue,
simple_annotation_plus6, leiden
```

For false-positive audits, also pass columns such as:

```bash
--obs-column has_TRA_TRB_paired \
--obs-column has_any_ab_tcr \
--obs-column annotation \
--obs-column cell_type
```

## Output Locations

For `--dataset-id TEST_DATASET`, outputs are written under:

```text
Integrated_dataset/tables/gdT_prediction/external_tests/TEST_DATASET/
Integrated_dataset/figures/gdT_prediction/external_tests/TEST_DATASET/
Integrated_dataset/logs/gdT_prediction/external_tests/TEST_DATASET/
gdT_prediction/external_tests/TEST_DATASET/index.html
```

Main output table:

```text
Integrated_dataset/tables/gdT_prediction/external_tests/TEST_DATASET/gdt_predictions.csv.gz
```

Key output columns:

- `cell_id`
- `gdt_score`
- `gdt_score_before_penalty`
- `predicted_gdT`
- `death_penalty_applied`
- `row_sum_x`
- `n_detected_genes_x`
- `missing_model_gene_count`
- copied obs columns, when present

## Shell Example

Run:

```bash
bash Integrated_dataset/models/gdT_prediction_classifier/gse144469_holdout_tcrgene/examples/apply_model_to_h5ad.sh \
  /path/to/test_dataset.h5ad \
  TEST_DATASET
```

This is a thin wrapper around `predict_with_selected_gdt_model.py`.

## Minimal Python Example

Run:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  Integrated_dataset/models/gdT_prediction_classifier/gse144469_holdout_tcrgene/examples/minimal_python_prediction.py \
  --input-h5ad /path/to/test_dataset.h5ad \
  --output-csv /tmp/test_dataset_gdt_predictions.csv.gz
```

This example writes only a compact prediction CSV and prints an overall count.
Use the main wrapper for full reports.

## False-Positive Audit Example

After prediction, audit predicted gdT calls inside likely false-positive
compartments:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  Integrated_dataset/models/gdT_prediction_classifier/gse144469_holdout_tcrgene/examples/audit_false_positive_predictions.py \
  --prediction-csv Integrated_dataset/tables/gdT_prediction/external_tests/TEST_DATASET/gdt_predictions.csv.gz \
  --output-csv Integrated_dataset/tables/gdT_prediction/external_tests/TEST_DATASET/false_positive_audit.csv \
  --annotation-column simple_annotation_plus6 \
  --group-column source_gse_id \
  --group-column library_id
```

Interpretation:

- cells with alpha-beta TCR evidence and `predicted_gdT == True` are counted as
  TCRAB false positives for this audit
- cells with NK-like annotation and `predicted_gdT == True` are counted as NK
  false positives for this audit
- this audit depends on copied obs columns being present in `gdt_predictions.csv.gz`

## Validation Context

The model was evaluated as a high-specificity internal cross-dataset classifier.
The validation cohort held out:

- all primary gold cells from `GSE144469`
- paired-TCRAB/no-gdTCR negatives from `GSE254249`
- cord-blood gdT positives from `GDT_2020AUG_woCOV`

Reported validation summary from the training run:

```text
validation cells: 250,686
positives: 16,501
negatives: 234,185
precision: about 95.7%
recall: about 79.8%
specificity: about 99.75%
threshold: 0.9769781827926636
```

This is not definitive external validation. For a new dataset, always report:

- model-gene availability
- predicted gdT count and fraction
- predicted gdT distribution by sample/library
- TCRAB-positive false positives, if TCRAB metadata exists
- NK-cell false positives, if NK annotation exists
- whether the matrix was raw counts or already transformed

## Common Failure Modes

- Many model genes missing: prediction may be unreliable.
- `X` is log-normalized instead of counts: scores are not on the training scale.
- Gene symbols are Ensembl IDs instead of gene symbols: model genes will appear
  missing.
- Test data contains non-T cells: NK-like and low-CD3 false positives should be
  audited.
- Very small datasets: predicted fractions may be unstable.

