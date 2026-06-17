# gdTAI v2.0 Manual

## Purpose

gdTAI v2.0 is the NK-optimized gdT-cell prediction package for local H5AD testing. It keeps the trained individual TCR-gene classifier frozen and exposes two operating modes so the user can choose between recovery and purity.

## Files

- Model pickle: `Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v2.0/gdTAI_v2_model.pkl`
- Inference script: `predict_with_gdtai_v2.py`
- Manifest: `Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v2.0/model_manifest.json`
- Mode metrics: `Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v2.0/mode_metrics.csv`
- Feature genes: `Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v2.0/feature_genes.csv`

Model SHA256:

```text
c7f6e7d1d1c0c8de9e7abbaa7c5b70c1070b43c964ccdb94e414bc5278e13c3e
```

Security note: `gdTAI_v2_model.pkl` is a Python pickle. Load only this trusted file from this repository or this server.

## Operating Modes

### `high_f1`

Use this mode when the goal is best validation F1 or higher gdT recovery.

- Uses the frozen base gdTAI threshold for all cells: `0.9063586592674255`
- Does not require annotation metadata
- Strict validation: precision `0.9478`, recall `0.8701`, F1 `0.9073`
- Full 5 million cell atlas known-FP proxy: NK FP `5975`, paired-TCRAB FP `2551`, known-FP fraction `0.0206`

### `high_purity`

Use this mode when the goal is highest purity / lower false positives, especially lower NK false positives.

- Uses annotation-specific thresholds from the NK-optimized wrapper
- Recommended annotation column: `simple_annotation_plus6`
- Thresholds: gdT `0.9063586592674255`, CD8 `0.93`, CD4 `0.97`, Treg disabled, NK `0.995`, other `0.97`
- Strict validation: precision `0.9515`, recall `0.8628`, F1 `0.9050`
- Full 5 million cell atlas known-FP proxy: NK FP `1557`, paired-TCRAB FP `2084`, known-FP fraction `0.0101`

If the annotation column is missing in `high_purity` mode, the script applies the strict `other` threshold to all cells. This is conservative but can lose many gdT cells.

## Input Requirements

Required:

- H5AD readable by `h5py`
- gene symbols in `var/_index`
- count-like expression in `X`
- sparse CSR `X` preferred; dense `X` also works

Strongly recommended:

- raw counts in `X`
- `simple_annotation_plus6` or equivalent annotation for `high_purity`
- `source_gse_id`, `sample_id`, `library_id`, and tissue columns for audits
- TCRAB/NK metadata if false-positive audits are needed

If `X` is already normalized or log-transformed, scores may not match the training scale.

## Quick Start

From the repository root:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python   predict_with_gdtai_v2.py   --input-h5ad /path/to/test_dataset.h5ad   --mode high_f1   --dataset-id TEST_DATASET   --chunk-size 50000
```

For lower false positives:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python   predict_with_gdtai_v2.py   --input-h5ad /path/to/test_dataset.h5ad   --mode high_purity   --annotation-column simple_annotation_plus6   --dataset-id TEST_DATASET   --chunk-size 50000
```

## Output Locations

For dataset `TEST_DATASET` and mode `high_f1`, outputs are written under:

```text
Integrated_dataset/tables/gdT_prediction/external_tests/TEST_DATASET_high_f1/
Integrated_dataset/figures/gdT_prediction/external_tests/TEST_DATASET_high_f1/
Integrated_dataset/logs/gdT_prediction/external_tests/TEST_DATASET_high_f1/
gdT_prediction/external_tests/TEST_DATASET_high_f1/index.html
```

For `high_purity`, the suffix is `_high_purity`.

Main prediction table:

```text
Integrated_dataset/tables/gdT_prediction/external_tests/<DATASET>_<MODE>/gdtai_v2_predictions.csv.gz
```

Key columns:

- `cell_id`
- `gdtai_mode`
- `gdtai_score`
- `threshold_used`
- `threshold_annotation`
- `predicted_gdT`
- `row_sum_x`
- `n_detected_genes_x`
- `missing_model_gene_count`

## Choosing A Mode

Use `high_f1` when:

- missing gdT cells is the larger problem
- the output is a discovery list
- later filtering or validation is available

Use `high_purity` when:

- false positives are costly
- NK-cell false positives must be reduced
- the output is a high-confidence gdT list for downstream analysis

A practical workflow is to run both modes, report both counts, and treat `high_purity` predictions as high-confidence candidates while using `high_f1` as the broader recovery set.

## Validation Context

The strict validation cohort held out:

- all primary gold cells from `GSE144469`
- paired-TCRAB/no-gdTCR negatives from `GSE254249`
- cord-blood gdT positives from `GDT_2020AUG_woCOV`

This remains internal validation. For external data, always report matrix state, model-gene availability, predicted counts by sample/library, and FP proxies where metadata exist.
