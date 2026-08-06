# gdTAI v2.0 Training And Tuning Methodology

## Overview

gdTAI v2.0 is a dual-operating-point release of the gdTAI classifier. It does not retrain a new classifier relative to the accepted gdTAI model. Instead, it packages:

- the frozen selected individual-gene gdTAI model
- a `high_f1` operating mode using the selected global model threshold
- a `high_purity` operating mode using the NK-optimized annotation-specific threshold wrapper

The model was trained and tuned on the integrated 5,128,904-cell T/NK atlas. All training, tuning, validation, and full-atlas application steps were read-only with respect to the H5AD inputs.

Primary scripts:

- Ground-truth and score-strategy baseline: `run_gdt_prediction_package_evaluation.py`
- Frozen gdTAI base model training: `run_gdt_gse144469_holdout_tcrgene_classifier.py`
- NK-optimized annotation-specific threshold wrapper: `run_gdtai_annotation_specific_cascade.py`
- v2.0 inference package: `predict_with_gdtai_v2.py`

Primary v2.0 artifacts:

- Model package: `Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v2.0/gdTAI_v2_model.pkl`
- Usage manual: `Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v2.0/MANUAL.md`
- Feature list: `Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v2.0/feature_genes.csv`
- Mode metrics: `Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v2.0/mode_metrics.csv`

## 1. Ground-Truth Definition

Ground truth was intentionally defined from dataset provenance and TCR metadata, not from the Phase 4 TRD/TRAB scores and not from the gdTAI model score.

Required metadata columns included:

- `source_gse_id`
- `library_id`
- `sample_id`
- `has_TRA_TRB_paired`
- `has_TRG_TRD_paired`
- `has_any_ab_tcr`
- `has_any_gd_tcr`
- `TRG_cdr3`
- `TRD_cdr3`

### 1.1 Productive gdTCR Evidence Repair

Before assigning labels, report-local gdTCR evidence was corrected for known missing-token artifacts.

For `GSE144469`, the integrated metadata had literal missing tokens in some TCR fields. The raw `has_any_gd_tcr` value was inflated to all cells. The correction used:

- a TRG proxy from `has_TRG_TRD_paired` or clean nonempty `TRG_cdr3`
- clean nonempty `TRD_cdr3`
- corrected paired gdTCR evidence as TRG proxy AND clean TRD evidence
- corrected any gdTCR evidence as TRG proxy OR clean TRD evidence

Audit counts:

| source_gse_id | n_cells | raw_has_any_gd_tcr | corrected_has_any_gd_tcr | raw_has_TRG_TRD_paired | corrected_has_TRG_TRD_paired | clean_TRD_cdr3_nonempty | has_TRA_TRB_paired |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| GSE144469 | 107,068 | 107,068 | 11,133 | 10,556 | 4,944 | 5,521 | 65,616 |
| HRA005041 | 766,639 | 14,119 | 14,119 | 5,471 | 5,471 | 7,576 | 547,970 |

### 1.2 gdTCR-Sequenced Sublibrary Definition

For `HRA005041` and `GSE144469`, sublibraries were defined by:

1. `library_id` when present
2. otherwise `sample_id`

A sublibrary was treated as gdTCR-sequenced if it contained any productive gamma-delta TCR evidence:

- corrected `has_any_gd_tcr == True`, or
- nonempty `TRG_cdr3`, or
- nonempty `TRD_cdr3`

Audit counts:

| source_gse_id | n_sublibraries | n_gdtcr_sequenced_sublibraries | n_cells |
| --- | ---: | ---: | ---: |
| GSE144469 | 37 | 36 | 107,068 |
| HRA005041 | 321 | 8 | 766,639 |

### 1.3 Primary Positive Labels: `gdT_gold`

Cells were labeled `gdT_gold` if they met either of these rules.

Sorted gdT datasets:

- all cells from `GDTlung2023july_7p`
- all cells from `MalteGDT`
- all cells from `GDT_2020AUG_woCOV`

TCR-sequenced datasets:

- source was `HRA005041` or `GSE144469`
- cell was in a gdTCR-sequenced sublibrary
- corrected paired `TRG/TRD` evidence was present
- no alpha-beta TCR evidence was present (`has_any_ab_tcr == False`)

### 1.4 Primary Negative Labels: `abT_gold`

Cells were labeled `abT_gold` if:

- source was `HRA005041` or `GSE144469`
- paired `TRA/TRB` evidence was present
- corrected any gdTCR evidence was absent

These cells were treated as high-confidence alpha-beta T-cell negatives.

### 1.5 Sensitivity-Only Labels: `gdT_silver`

Cells were labeled `gdT_silver` if:

- source was `HRA005041` or `GSE144469`
- cell was in a gdTCR-sequenced sublibrary
- clean nonempty `TRD_cdr3` was present
- corrected paired `TRG/TRD` evidence was absent
- no alpha-beta TCR evidence was present

Silver cells were not used to select model thresholds and were not included in primary validation metrics. They were used only as sensitivity evidence.

### 1.6 Unlabeled Or Ambiguous Cells

Cells not satisfying `gdT_gold`, `abT_gold`, or `gdT_silver` rules were retained for full-atlas prediction but excluded from primary supervised metrics.

### 1.7 Ground-Truth Counts And Conflict Checks

Whole-atlas ground-truth classes:

| class | n_cells | fraction |
| --- | ---: | ---: |
| unlabeled_or_ambiguous | 4,465,731 | 0.8707 |
| abT_gold | 603,203 | 0.1176 |
| gdT_gold | 57,776 | 0.0113 |
| gdT_silver | 2,194 | 0.0004 |

Labeled-source counts:

| source_gse_id | gdT_gold | abT_gold | gdT_silver |
| --- | ---: | ---: | ---: |
| GDT_2020AUG_woCOV | 25,904 | 0 | 0 |
| GDTlung2023july_7p | 15,175 | 0 | 0 |
| MalteGDT | 7,800 | 0 | 0 |
| GSE144469 | 4,003 | 60,175 | 415 |
| HRA005041 | 4,894 | 543,028 | 1,779 |

No overlap conflicts were observed:

- `gdT_gold` and `abT_gold`: 0 cells
- `gdT_gold` and `gdT_silver`: 0 cells
- `abT_gold` and `gdT_silver`: 0 cells

## 2. Validation Design To Reduce Overfitting Risk

The selected gdTAI training run was designed after concern that random splits would leak dataset-specific patterns into validation. The final split held out entire sources or biologically distinct validation cohorts.

Validation cohorts:

- all primary gold-label cells from `GSE144469`
- cord-blood `gdT_gold` cells from `GDT_2020AUG_woCOV`
- paired-TCRAB / no-gdTCR negatives from `GSE254249`

`GDTlung2023july_7p` was excluded from training and threshold tuning because library quality was considered suboptimal. It was retained as a sensitivity cohort rather than as a training source.

Split counts:

| split | n_cells | gdT_gold | abT_gold | gdT prevalence |
| --- | ---: | ---: | ---: | ---: |
| train | 455,302 | 20,880 | 434,422 | 0.0459 |
| tune | 113,826 | 5,220 | 108,606 | 0.0459 |
| validation_primary_GSE144469 | 64,178 | 4,003 | 60,175 | 0.0624 |
| validation_gdT_GDT_2020AUG_woCOV_cord_blood | 12,498 | 12,498 | 0 | 1.0000 |
| validation_abT_GSE254249_paired_TCRAB_no_gdTCR | 174,010 | 0 | 174,010 | 0.0000 |
| validation_combined | 250,686 | 16,501 | 234,185 | 0.0658 |
| sensitivity_excluded_GDTlung2023july_7p | 15,175 | 15,175 | 0 | 1.0000 |

The train/tune pool was source-label stratified and split 80/20:

| source_gse_id | label | n_cells | train | tune |
| --- | --- | ---: | ---: | ---: |
| GDT_2020AUG_woCOV | gdT_gold | 13,406 | 10,725 | 2,681 |
| HRA005041 | abT_gold | 543,028 | 434,422 | 108,606 |
| HRA005041 | gdT_gold | 4,894 | 3,915 | 979 |
| MalteGDT | gdT_gold | 7,800 | 6,240 | 1,560 |

## 3. Feature Engineering

The final classifier does not use Phase 4 module scores as features. It also does not use TCR metadata fields, paired-chain booleans, source ID, tissue, sample ID, library ID, UMAP coordinates, or scVI annotations as model features.

Expression features were rebuilt directly from H5AD `X`:

1. read sparse CSR count-like rows
2. compute row total
3. normalize selected genes to counts per 10,000
4. apply `log1p`
5. fill missing model genes with zero during inference

For a cell `i` and gene `g`:

```text
feature(i,g) = log1p(raw_count(i,g) * 10000 / total_counts(i))
```

### 3.1 Gene Selection Rule

Feature genes were selected deterministically from genes present in the atlas:

1. collect individual TCR alpha, beta, gamma, and delta genes:
   - `TRAV`, `TRAJ`, `TRAC`
   - `TRBV`, `TRBJ`, `TRBC`
   - `TRGV`, `TRGJ`, `TRGC`
   - `TRDV`, `TRDD`, `TRDJ`, `TRDC`
2. sort genes by priority:
   - TRD first
   - TRG second
   - TRA third
   - TRB fourth
3. append penalty/control genes if present:
   - `FOXP3`
   - `CD4`
   - `CD3D`
   - `CD3E`
   - `CD3G`
4. enforce the user-requested cap of at most 300 genes

The final v2.0 package uses 187 expression features:

| feature family | n_genes |
| --- | ---: |
| TRD | 7 |
| TRG | 14 |
| TRA | 93 |
| TRB | 68 |
| FOXP3/CD4 controls | 2 |
| CD3 controls | 3 |

Important distinction: FOXP3, CD4, and low-CD3 post-hoc death penalties were tested in candidate models, but the selected v2.0 base model is `logistic_individual_TCRABGD_genes`. It includes FOXP3/CD4/CD3 expression as features but does not apply the external hard penalty rule at inference time.

## 4. Candidate Models

Candidate models were trained only on the train split and thresholds were selected only on the tune split.

Baselines:

- `baseline_individual_TRD_gene_sum`
- `baseline_individual_gd_minus_ab_gene_sum`

Learned candidates:

- `xgb_individual_TCRABGD_genes`
- `xgb_individual_TCRABGD_genes_with_FOXP3_CD4_lowCD3_penalty`
- `logistic_individual_TCRABGD_genes`

XGBoost settings:

```text
n_estimators = 420
max_depth = 4
learning_rate = 0.045
subsample = 0.85
colsample_bytree = 0.85
min_child_weight = 2.0
reg_lambda = 1.5
objective = binary:logistic
eval_metric = logloss
tree_method = hist
n_jobs = 32
scale_pos_weight = n_negative / n_positive
```

Logistic-regression settings:

```text
StandardScaler()
LogisticRegression(
    solver = lbfgs,
    class_weight = balanced,
    max_iter = 1000,
    random_state = 20260501
)
```

Training negatives were downsampled to manage imbalance:

- all available train positives were retained
- train negatives were sampled up to `max(250,000, 4 * n_positive)`
- for this run, the training sample contained 20,880 positives and 250,000 negatives

## 5. Base-Model Threshold Selection

For each candidate, a continuous score was generated on the tune split. The threshold was selected from ROC thresholds by maximizing F1 (`beta = 1.0`). The tie-breaking order was:

1. higher F1
2. higher precision
3. higher recall
4. higher threshold

This threshold policy was intentionally different from the earlier conservative TRD-minus-TRAB score report, which used a specificity floor and F0.5. For gdTAI v2.0, the requested default recovery mode was best validation F1.

The selected base threshold was:

```text
0.9063586592674255
```

## 6. Base-Model Selection

Learned models were accepted only if they outperformed the best individual-gene baseline on the strict validation F1. The best baseline was `baseline_individual_TRD_gene_sum`.

Validation acceptance:

| model | validation F1 | baseline F1 | delta F1 | precision | recall | specificity | accepted |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| xgb_individual_TCRABGD_genes | 0.9004 | 0.8442 | 0.0562 | 0.9450 | 0.8598 | 0.9965 | True |
| xgb_individual_TCRABGD_genes_with_FOXP3_CD4_lowCD3_penalty | 0.8971 | 0.8442 | 0.0529 | 0.9454 | 0.8535 | 0.9965 | True |
| logistic_individual_TCRABGD_genes | 0.9073 | 0.8442 | 0.0631 | 0.9478 | 0.8701 | 0.9966 | True |

The selected frozen base model was:

```text
logistic_individual_TCRABGD_genes
```

Audit clarification (2026-08-06): candidate acceptance and final algorithm
selection used the combined cohort named `validation` (strict validation F1,
then precision and recall). The metrics below are therefore model-selection
estimates, not performance from an untouched final test set. Threshold selection
remained confined to the tune split.

Strict validation performance for the selected base model:

| metric | value |
| --- | ---: |
| n_cells | 250,686 |
| n_positive | 16,501 |
| n_negative | 234,185 |
| predicted_positive | 15,147 |
| TP | 14,357 |
| FP | 790 |
| TN | 233,395 |
| FN | 2,144 |
| precision | 0.9478 |
| recall | 0.8701 |
| specificity | 0.9966 |
| F1 | 0.9073 |
| F0.5 | 0.9312 |
| balanced accuracy | 0.9333 |
| MCC | 0.9020 |
| ROC-AUC | 0.9812 |
| PR-AUC | 0.9181 |

## 7. Full-Atlas False-Positive Audits

After selecting the base model, it was applied to the full 5,128,904-cell atlas. The full-atlas labels are incomplete, so two conservative false-positive proxy classes were audited:

- predicted gdT cells with paired `TRA/TRB` evidence
- predicted gdT cells annotated as NK cells by the simple scVI annotation

For the selected base model:

| measure | value |
| --- | ---: |
| predicted putative gdT | 410,928 |
| predicted fraction | 0.0801 |
| paired-TCRAB FP proxy | 2,551 |
| paired-TCRAB FP ratio among paired-TCRAB cells | 0.0020 |
| NK FP proxy | 5,975 |
| NK FP ratio among NK cells | 0.0233 |
| paired-TCRAB or NK known-FP proxy | 8,475 |
| known-FP fraction of predictions | 0.0206 |

The original TRD-minus-TRAB score strategy had fewer total predictions but a higher known-FP fraction:

| strategy | predicted putative gdT | paired-TCRAB FP | NK FP | known-FP fraction |
| --- | ---: | ---: | ---: | ---: |
| selected gdTAI base model | 410,928 | 2,551 | 5,975 | 0.0206 |
| original TRD-minus-TRAB strategy | 222,570 | 6,466 | 4,422 | 0.0484 |

## 8. v2.0 NK-Optimized Threshold Tuning

The NK-optimized v2.0 layer did not retrain model weights. It froze the selected logistic base model and tuned thresholds by simple scVI annotation.

The annotation-specific threshold wrapper used these annotation groups:

- `gdT_cell`
- `CD8_T`
- `CD4_T`
- `Treg`
- `NK_cell`
- `other`

The grid was deliberately small and interpretable:

- gdT threshold: current base threshold only
- CD8 threshold: current base threshold, 0.93, 0.95
- CD4 threshold: current base threshold, 0.93, 0.95, 0.97
- Treg, NK, other thresholds: 0.97, 0.985, 0.995, disabled

The tuning rule was:

1. evaluate each threshold combination on the tune split
2. require recall to be at least `max(0.86, current_tune_recall - 0.01)`
3. require F1 to be at least `current_tune_F1 - 0.005`
4. among passing candidates, select the candidate with lowest tune FP
5. tie-break by higher precision, higher F1, then higher recall
6. if no candidate passed, fall back to best tune F1

The selected annotation-specific thresholds were:

| annotation group | threshold |
| --- | ---: |
| gdT_cell | 0.9063586592674255 |
| CD8_T | 0.93 |
| CD4_T | 0.97 |
| Treg | disabled |
| NK_cell | 0.995 |
| other | 0.97 |

Tune split comparison:

| strategy | predicted_positive | FP | precision | recall | specificity | F1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| current_gdTAI | 5,245 | 147 | 0.9720 | 0.9766 | 0.9986 | 0.9743 |
| annotation_specific_cascade | 5,155 | 108 | 0.9790 | 0.9669 | 0.9990 | 0.9729 |

Strict validation comparison:

| strategy | predicted_positive | TP | FP | FN | precision | recall | specificity | F1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| current_gdTAI | 15,147 | 14,357 | 790 | 2,144 | 0.9478 | 0.8701 | 0.9966 | 0.9073 |
| annotation_specific_cascade | 14,963 | 14,237 | 726 | 2,264 | 0.9515 | 0.8628 | 0.9969 | 0.9050 |

Promotion gates for the annotation-specific cascade:

| gate | observed | required | passed |
| --- | ---: | ---: | --- |
| validation recall | 0.8628 | >= 0.86 | True |
| validation recall delta vs current | -0.0073 | >= -0.01 | True |
| validation F1 delta vs current | -0.0023 | >= -0.005 | True |
| full-atlas NK FP reduction vs current | 0.7394 | >= 0.50 | True |
| full-atlas paired-TCRAB FP no worse than current | 0.8169 | <= 1.00 | True |
| full-atlas known-FP fraction lower than current | -0.0105 | < 0 | True |

Full-atlas comparison:

| strategy | predicted putative gdT | paired-TCRAB FP | NK FP | known-FP fraction |
| --- | ---: | ---: | ---: | ---: |
| current_gdTAI | 410,928 | 2,551 | 5,975 | 0.0206 |
| annotation_specific_cascade | 359,857 | 2,084 | 1,557 | 0.0101 |
| original TRD-minus-TRAB strategy | 222,570 | 6,466 | 4,422 | 0.0484 |

## 9. gdTAI v2.0 Operating Modes

The v2.0 package exposes the two accepted operating points in one trusted pickle.

### `high_f1`

Purpose:

- maximize validation F1 among the packaged modes
- recover more putative gdT cells
- does not require annotation metadata

Rule:

```text
predict gdT if base_model_probability >= 0.9063586592674255
```

Strict validation:

- precision: 0.9478
- recall: 0.8701
- F1: 0.9073

Full-atlas recovery:

- predicted putative gdT: 410,928
- paired-TCRAB FP proxy: 2,551
- NK FP proxy: 5,975
- known-FP fraction of predictions: 0.0206

### `high_purity`

Purpose:

- reduce false positives, especially NK false positives
- preserve validation F1 within the accepted tolerance
- requires a compatible simple cell-type annotation for best behavior

Rule:

```text
gdT_cell: score >= 0.9063586592674255
CD8_T:   score >= 0.93
CD4_T:   score >= 0.97
Treg:    disabled
NK_cell: score >= 0.995
other:   score >= 0.97
```

Strict validation:

- precision: 0.9515
- recall: 0.8628
- F1: 0.9050

Full-atlas recovery:

- predicted putative gdT: 359,857
- paired-TCRAB FP proxy: 2,084
- NK FP proxy: 1,557
- known-FP fraction of predictions: 0.0101

If the annotation column is missing, `predict_with_gdtai_v2.py` applies the strict `other` threshold to all cells in `high_purity` mode. This is conservative and can reduce recovery.

## 10. What Was Not Used As Model Input

To limit leakage and keep the model portable, the final selected base model did not use:

- Phase 4 `TRD_score`
- Phase 4 `TRAB_score`
- `TRD_minus_TRAB`
- paired-chain metadata
- productive TCR booleans
- CDR3 strings
- source or dataset ID
- sample or library ID
- tissue
- UMAP coordinates
- scVI latent values
- scVI annotation as a base-model feature

The `high_purity` mode uses simple annotation only after the frozen base score is computed, as a thresholding policy rather than as a learned expression feature.

## 11. Main Limitations

1. Validation is internal to the atlas build, not a completely external laboratory dataset.
2. The validation design reduces overfitting risk by holding out complete sources or cohorts, but it cannot rule out all preprocessing-level leakage because the data were processed in one atlas.
3. Ground truth depends on TCR metadata quality. GSE144469 required explicit correction for missing-token artifacts.
4. Silver gdT labels are lower confidence and were not used for threshold selection.
5. `high_purity` depends on a compatible simple cell-type annotation. Annotation drift on external data can change the purity/recovery tradeoff.
6. The model expects raw count-like `X`. Normalized or log-transformed `X` can shift the score scale.
7. False-positive audits based on paired TCRAB and NK annotation are strong proxies but not exhaustive ground truth for every cell in the atlas.

## 12. Reproducibility Pointers

Primary tables:

- Ground-truth counts: `Integrated_dataset/tables/gdT_prediction/ground_truth_class_counts.csv`
- Ground-truth by source: `Integrated_dataset/tables/gdT_prediction/ground_truth_by_source_gse.csv`
- TCR evidence correction audit: `Integrated_dataset/tables/gdT_prediction/tcr_evidence_correction_audit.csv`
- Base-model split summary: `Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene/split_overall.csv`
- Base-model feature manifest: `Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene/feature_manifest.csv`
- Base-model validation metrics: `Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene/validation_metrics.csv`
- Annotation-specific threshold grid: `Integrated_dataset/tables/gdT_prediction/gdtai_annotation_specific_cascade/annotation_threshold_tuning_grid.csv`
- v2.0 mode metrics: `Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v2.0/mode_metrics.csv`

Primary reports:

- `gdT_prediction/gse144469_holdout_tcrgene_report.html`
- `gdT_prediction/gdtai_annotation_specific_cascade_report.html`
