# gdTAI v3 TRDC/NK Guard Use Protocol

## 1. Purpose

Use this protocol to apply the packaged gdTAI v3 TRDC/NK guard candidate to a new single-cell H5AD and produce per-cell gdT scores and binary predictions.

Selected model: `v3_round14_v2_score_trdc_gate_fixed_0p936`

Decision threshold: `0.936`

## 2. Required Files

Keep these files together in the repository checkout:

```text
run_gdtai_v3_trdc_nk_guard_classifier.py
Integrated_dataset/models/gdT_prediction_classifier/gdtai_v3_trdc_nk_guard/best_candidate_model.pkl
Integrated_dataset/models/gdT_prediction_classifier/gdtai_v3_trdc_nk_guard/feature_genes.csv
```

The pickle imports the `ConditionalGateModel` class from `run_gdtai_v3_trdc_nk_guard_classifier.py`, so that script must be importable from the working directory or `PYTHONPATH`.

## 3. Environment

Canonical environment for this project:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
```

Minimum Python dependencies used by the example:

```text
h5py
numpy
pandas
scikit-learn
```

## 4. Input H5AD Requirements

The input H5AD must contain raw counts in `layers["counts"]`. Do not use normalized/log `X` for external inference. The example script fails if `layers["counts"]` is absent.

The input H5AD must contain all model genes listed in `feature_genes.csv` with `feature_type == gene_log1p_cp10k`. Gene names are matched exactly against `var["_index"]`.

Expected matrix transform:

```text
log1p(raw_counts * 10000 / total_counts_per_cell)
```

## 5. Run Inference

Example command:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  Integrated_dataset/models/gdT_prediction_classifier/gdtai_v3_trdc_nk_guard/examples/predict_h5ad_counts.py \
  --input-h5ad /path/to/input.h5ad \
  --output-csv /path/to/gdtai_v3_predictions.csv.gz
```

The output table contains:

```text
cell_id
gdtai_v3_score
gdtai_v3_threshold
gdtai_v3_predicted
tcr_gene_quadrant
row_sum_counts_layer
n_detected_genes_counts_layer
```

## 6. Interpretation

A cell is called putative gdT when:

```text
gdtai_v3_score >= 0.936
```

The score is the post-gate probability-like score, not a calibrated biological probability. The gate deliberately suppresses weak `TRDC-only` / NK-like and abT-biased cells.

## 7. Recommended QC After Inference

Always report these summaries after applying the model to a new dataset:

- total cells and predicted cells
- predicted fraction
- predicted cells by sample/source/tissue when metadata is available
- predicted cells by `tcr_gene_quadrant`
- predicted cells with NK annotations, if annotations exist
- predicted paired-TCRAB/no-gdT cells, if TCR-seq metadata exists
- score distribution and threshold line

If TCR metadata exists, separately report predicted paired-TCRAB cells. If TCR metadata is absent, treat FP estimates as uncertain and do not claim exact abT FP counts.

## 8. Known Operating Characteristics

From the final project evaluation:

Full 5.13M atlas primary-gold evaluation:

- predicted: `251,356`
- primary-gold recall: `0.8056`
- primary-gold precision: `0.9752`
- estimated FP / predictions: `0.0406`

Independent external primary evaluation:

- predicted: `856`
- precision: `0.9509`
- recall: `0.8586`
- F1: `0.9024`
- FP / predictions: `0.0491`

## 9. Caveats

- This is a review candidate, not a promoted `gdTAI_v3.0` release.
- It met the user target of full-atlas primary-gold recall above `0.8` and estimated FP below `5%`.
- It failed the older strict promotion gate because external NK false positives were higher than v2 high-purity.
- It does not use TCR-seq metadata for prediction. TCR metadata is only for evaluation/QC.
- Full-atlas recall is measured against known primary gold labels, not against every true biological gdT cell in the atlas.

## 10. Reproducibility Checklist

Before reporting results, record:

- model file path and checksum
- input H5AD path and whether `layers["counts"]` was used
- number of cells and genes in the input
- number of model genes available
- threshold used: `0.936`
- total predicted cells
- QC summaries listed above
