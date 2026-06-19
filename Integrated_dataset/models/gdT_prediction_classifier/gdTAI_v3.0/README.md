# gdTAI v3 TRDC/NK Guard Candidate

This directory contains the selected gdTAI v3 review candidate. The selected mode is `v3_round14_v2_score_trdc_gate_fixed_0p936` with threshold `0.936`.

This is not an XGBoost model. The saved object is a conditional gate around the gdTAI v2 balanced logistic-regression score. The gate suppresses weak `TRDC-only` / NK-like or abT-biased cells, then applies the fixed recall-rescue threshold.

## Files

- `best_candidate_model.pkl`: pickle payload containing the model object, threshold, feature names, gene list, engineered feature names, and evaluation metadata.
- `USAGE_PROTOCOL.md`: detailed inference protocol and validation checklist.
- `METHODOLOGY.md`: training, tuning, candidate selection, and caveats.
- `feature_genes.csv`: ordered feature manifest expected by the model.
- `mode_metrics.csv`: selected full-atlas, external, and validation metrics.
- `model_manifest.json`: machine-readable package summary.
- `external_false_positive_groups.csv`, `external_recall_groups.csv`: selected external stress-group metrics.
- `candidate_target_selection.csv`: full-atlas target-selection table for all candidate rounds.
- `examples/predict_h5ad_counts.py`: minimal standalone inference example.

## Headline Performance

Full atlas primary-gold evaluation:

- predicted cells: `251,356`
- primary-gold recall: `0.8056`
- primary-gold precision: `0.9752`
- estimated total abT FP: `10,211.9`
- estimated FP / predictions: `0.0406`

Independent external primary evaluation:

- precision: `0.9509`
- recall: `0.8586`
- F1: `0.9024`
- FP / predictions: `0.0491`

## Important Caveat

This candidate is packaged for review/use because it met the user target: full-atlas primary-gold recall above `0.8` and estimated FP fraction below `5%`. It was not promoted under the older strict gates because external NK false positives were higher than v2 high-purity.
