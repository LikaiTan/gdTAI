# gdTAI v3 TRDC/NK Guard Decision Report

## Executive decision

**Decision: do not promote this v3 candidate as the default release model yet.**

The selected candidate `v3_round14_v2_score_trdc_gate_fixed_0p936` met the user target on the full atlas (80.56% primary-gold recall and 4.06% estimated FP/prediction), but it failed the explicit promotion gate for NK false positives on the independent external test (9 NK FP vs 3 for v2 high-purity).

## What changed

- `GSE144469` was moved into train/tune: 63,764 cells, including 4,003 gdT positives and 59,761 abT negatives.
- Other validation cohorts were retained: 186,195 held-out non-GSE144469 cells.
- The independent external H5AD was used only after model and threshold selection.
- External inference used `layers["counts"]` with log1p(counts per 10,000), not normalized `X`.
- Hard NK negatives followed the requested rule: 1,453 NK `TRDC+TRDV-` cells from 23 datasets were retained; 31,117 NK cells not matching that expression rule were excluded. Also, 4 expression-passing weak-TRDC NK+TCRAB-overlap cells were excluded.

**Step conclusion:** the evaluation design is now aligned with the requested split policy, and the new external H5AD remains independent.

## Algorithm

- The selected model is not XGBoost.
- It is a conditional gate around the gdTAI v2 transcriptome score.
- It uses count-derived gene features plus engineered TCR/CD3/NK features.
- It does not require TCR-seq metadata at inference time.
- The operating threshold is `0.936`.

**Step conclusion:** this is a conservative v2-score-plus-guard candidate, not a fully new tree model.

## Candidate selection

The initial internal winner was `v3_round03_hist_gradient_target_recall80_fp05`, but it had only 72.70% full-atlas primary-gold recall. The final selected round had 80.56% recall and 4.06% estimated FP/prediction, passing the 80.5% recall-with-margin and 5.0% estimated-FP targets.

External testing was not used for model selection. For example, `v3_round08_hist_gradient_fixed_0p64` had the best v3 external F1 (0.9087), but it was not selected because the final decision rule prioritized the full-atlas target.

**Step conclusion:** round14 is the selected candidate because it is the first practical point satisfying the atlas-level recall and FP estimate target.

## External performance

| Metric | v3 selected | v2 high-purity | Interpretation |
| --- | --- | --- | --- |
| External precision | 95.1% | 93.9% | v3 higher |
| External recall | 85.9% | 86.7% | v3 slightly lower |
| External F1 | 0.9024 | 0.9018 | near-tie |
| External FP / predictions | 4.9% | 6.1% | v3 lower |
| External NK FP | 9 | 3 | v3 worse |

Compared with v2 high-purity, v3 reduced total external false positives from 53 to 42 and reduced paired-TCRAB false positives from 51 to 36. However, NK false positives increased from 3 to 9.

**Step conclusion:** the candidate is competitive overall, but it does not solve the NK-specific failure mode strongly enough.

## Full-atlas application

| Full-atlas metric | v3 selected | v2 high-purity |
| --- | --- | --- |
| Predicted putative gdT | 251,356 | 359,857 |
| Primary-gold recall | 80.56% | 80.56% |
| Primary-gold precision | 97.52% | 97.61% |
| Estimated total abT FP | 10,212 | 9,680 |
| Estimated FP / predictions | 4.06% | 2.69% |
| Predicted paired TCRAB | 5,580 | 5,184 |
| Predicted TRDC+/TRDV- | 17,120 | 130,846 |

The full-atlas recall is measured on primary gold labels only. The estimated false-positive count extrapolates observed paired-TCRAB false-positive behavior from TCR-seq sources to sources without TCR-seq, so it is useful but still an estimate.

**Step conclusion:** the candidate meets the numerical full-atlas target, but this does not override a direct external NK regression.

## Promotion gate

| Gate | Status | Evidence |
| --- | --- | --- |
| NK false positives lower than v2 high-purity | FAIL | v3 9 vs v2 high-purity 3 |
| TRDC+/TRDV- false-positive burden lower | PASS | v3 14 vs v2 high-purity 33 |
| External F1 not worse by more than 0.01 | PASS | v3 0.9024 vs v2 high-purity 0.9018 |
| Cytotoxic gdT recall preserved | PASS | v3 91.9% vs v2 high-purity 93.0% |
| Paired-TCRAB false positives not higher | PASS | v3 36 vs v2 high-purity 51 |

**Final conclusion:** keep this model as a packaged v3 candidate, not as promoted `gdTAI_v3.0`. The reason is not low overall performance; the reason is that the named failure mode, NK contamination, is worse than v2 high-purity on the independent external test.

## Detailed artifact paths

- Model package: `Integrated_dataset/models/gdT_prediction_classifier/gdtai_v3_trdc_nk_guard`
- Full result tables: `Integrated_dataset/tables/gdT_prediction/gdtai_v3_trdc_nk_guard`
- Figures: `gdT_prediction/assets/gdtai_v3_trdc_nk_guard`
- HTML report: `gdT_prediction/gdtai_v3_trdc_nk_guard_report.html`
- PDF report: `gdT_prediction/gdtai_v3_trdc_nk_guard_report.pdf`
