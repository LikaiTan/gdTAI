# gdT Classifier Local Handoff Plan

Audience: another agent or analyst working on the same server.

Purpose: test the selected gdT classifier on another local `.h5ad` without
rerunning training and without pushing anything to GitHub.

## 1. Required startup context

Before doing work in this repository, read:

1. `TNK_PIPELINE_RUNBOOK.md`
2. `TNK_PIPELINE_STATUS.md`
3. `TNK_PHASES_0_4_SCRIPT.md`

Current project state at handoff:

- milestone: post-Phase 4 review
- task type: downstream model testing/reporting
- canonical environment: `/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python`
- do not modify any milestone H5AD unless explicitly requested
- use `Integrated_dataset/` for canonical tables/logs/figures

## 2. Model artifact

Use this local model:

```text
/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/models/gdT_prediction_classifier/gse144469_holdout_tcrgene/selected_model.pkl
```

Model metadata:

- size: `647,129` bytes, about `632 KiB`
- SHA256: `ef7dc0cca09acfcee27e059782e26688f79538c53ba0621025c1ad6e2b24383d`
- model: `xgb_individual_TCRABGD_genes_with_FOXP3_CD4_lowCD3_penalty`
- features: `187`
- genes: `187`
- threshold: `0.9769781827926636`
- feature source: individual TCR A/B/G/D genes plus `FOXP3`, `CD4`, `CD3D`,
  `CD3E`, and `CD3G`
- not used as model inputs: Phase 4 `TRD_score`, `TRAB_score`,
  `TRD_minus_TRAB`, TCR metadata labels, or TCR CDR3 metadata

Security note:

- This is a Python pickle. Load only this trusted local file from this project.
  Do not load arbitrary external pickle files.

## 3. Reference code and reports

Training and validation code:

```text
/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/run_gdt_gse144469_holdout_tcrgene_classifier.py
/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/run_gdt_prediction_package_evaluation.py
/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/run_gdt_deg_tcr_classifier_training.py
```

Inference-only wrapper for local testing:

```text
/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/predict_with_selected_gdt_model.py
```

Validation report:

```text
/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/gdT_prediction/gse144469_holdout_tcrgene_report.html
/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/gdT_prediction/gse144469_holdout_tcrgene_report.pdf
```

Failure-mode audit:

```text
/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/gdT_prediction/gse144469_holdout_tcrgene_failure_mode_audit.html
/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/gdT_prediction/gse144469_holdout_tcrgene_failure_mode_audit.pdf
```

## 4. Expected input H5AD

The tester should provide a local `.h5ad` with:

- `X` as sparse CSR or another count-like expression matrix
- gene symbols in `var/_index`
- human TCR gene symbols matching the model feature genes
- preferably raw counts or count-like expression in `X`

Important:

- The model internally uses per-cell `log1p(count * 10000 / total_count)`
  features.
- If the test H5AD stores already-normalized/log-transformed `X`, results may
  not be comparable. The agent must document the matrix state before applying
  the model.
- Missing model genes should be filled with zero and reported.

## 5. Exact inference logic

For each cell:

1. Read the 187 model genes from `X`.
2. Compute row sum from `X`.
3. Convert each selected gene to:

```text
log1p(raw_value * 10000 / row_sum)
```

4. Build the feature matrix in exactly the saved `feature_names` order.
5. Run:

```python
score = model_object.predict_proba(X_features)[:, 1]
```

6. Apply the same death-penalty rule:

```text
if FOXP3 > 0.25, or CD4 > 0.75, or CD3D + CD3E + CD3G < 0.25:
    score = min(score, 0.03)
```

7. Call gdT-positive if:

```text
score >= 0.9769781827926636
```

## 6. Recommended local test outputs

Do not overwrite the existing training outputs. For a new dataset named
`DATASET_ID`, write to:

```text
Integrated_dataset/tables/gdT_prediction/external_tests/DATASET_ID/
Integrated_dataset/figures/gdT_prediction/external_tests/DATASET_ID/
Integrated_dataset/logs/gdT_prediction/external_tests/DATASET_ID/
gdT_prediction/external_tests/DATASET_ID/
```

Minimum tables:

- `gdt_predictions.csv.gz`
  - cell barcode/index
  - predicted probability
  - binary prediction
  - death-penalty flag
  - missing model-gene count
  - source/sample/library columns if available
- `prediction_summary_overall.csv`
- `prediction_summary_by_sample.csv`
- `prediction_summary_by_annotation.csv`, if annotation exists
- `model_gene_availability.csv`

Minimum figures:

- prediction score histogram
- predicted gdT fraction by sample/library
- UMAP colored by predicted gdT, if UMAP coordinates exist
- TRD/TRAB score scatter only if those score columns already exist in the test
  H5AD; do not recompute Phase 4 scores just for this test unless requested

Minimum log:

- input H5AD path
- matrix state assumption
- model path and SHA256
- model threshold
- number of cells
- number of available and missing model genes
- predicted gdT count and fraction
- caveats

## 7. Suggested command shape

Use the provided inference-only wrapper:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python predict_with_selected_gdt_model.py \
  --input-h5ad /path/to/test_dataset.h5ad \
  --model-pkl Integrated_dataset/models/gdT_prediction_classifier/gse144469_holdout_tcrgene/selected_model.pkl \
  --dataset-id DATASET_ID
```

The wrapper is read-only with respect to the input H5AD. It writes:

```text
Integrated_dataset/tables/gdT_prediction/external_tests/DATASET_ID/
Integrated_dataset/figures/gdT_prediction/external_tests/DATASET_ID/
Integrated_dataset/logs/gdT_prediction/external_tests/DATASET_ID/
gdT_prediction/external_tests/DATASET_ID/index.html
```

Optional useful flags:

```text
--chunk-size 50000
--obs-column COLUMN_NAME
--no-figures
--max-umap-points 300000
```

Smoke test already run:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python predict_with_selected_gdt_model.py \
  --input-h5ad Integrated_dataset/logs/_backed_write_test_input.h5ad \
  --dataset-id wrapper_smoke_test \
  --no-figures \
  --chunk-size 2
```

Smoke-test outputs:

```text
Integrated_dataset/tables/gdT_prediction/external_tests/wrapper_smoke_test/
Integrated_dataset/logs/gdT_prediction/external_tests/wrapper_smoke_test/
gdT_prediction/external_tests/wrapper_smoke_test/index.html
```

## 8. Validation context to report

The model was accepted on multi-cohort holdout validation:

- validation cells: `250,686`
- positives: `16,501`
- negatives: `234,185`
- precision: about `95.7%`
- recall: about `79.8%`
- specificity: about `99.75%`
- selected threshold: `0.9769781827926636`

Holdout cohorts:

- all primary gold cells from `GSE144469`
- paired-TCRAB/no-gdTCR negatives from `GSE254249`
- cord-blood gdT positives from `GDT_2020AUG_woCOV`

Known limitation:

- This is strong internal cross-dataset validation, not definitive external
  validation.

## 9. Failure-mode caveats

From the local failure-mode audit:

- false negatives are enriched for lower TRG/TRD signal, low-complexity/low-
  quality proxies, and death-penalty signals
- false positives are not globally enriched for low-quality proxies
- false positives are enriched for gdTCR-gene plus NK-like expression and may
  include biologically ambiguous or TCR-mixed cells
- explicit doublet-score columns were not present in `integrated_plus6.h5ad`;
  doublet conclusions are proxy-based

## 10. Practical recommendation

Use this model as a high-confidence gdT recovery classifier. For a new dataset,
report both:

- high-confidence predicted gdT cells at the saved threshold
- a lower-threshold exploratory set only if the user explicitly asks for higher
  recall

Do not silently retrain, retune the threshold, or write predictions back into
the input H5AD.
