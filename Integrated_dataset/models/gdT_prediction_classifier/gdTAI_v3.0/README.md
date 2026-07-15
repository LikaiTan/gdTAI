# gdTAI v3.0

This is the promoted balanced gdTAI v3 package.

- Model: `v3_round14_v2_score_trdc_gate_fixed_0p936`
- Threshold: `0.936`
- SHA256: `16dedc0081da9b8487887341232bcf8c9c9403dd3bbd72e04cab43d4cd7b2e09`
- Input transform: `log1p(raw counts per 10,000)`
- Features: 210 genes plus 16 engineered expression features
- Comparison decision: `gdT_prediction/gdtai_v3_round12_vs_round14/index.html`

## Headline Metrics

| Evaluation | Precision | Recall | F1 |
|---|---:|---:|---:|
| Full atlas primary gold | 0.9752 | 0.8056 | 0.8823 |
| Atlas held-out combined | 0.9854 | 0.8530 | 0.9144 |
| Independent external | 0.9509 | 0.8586 | 0.9024 |

## Files

- `gdTAI_v3_model.pkl`: promoted model payload
- `model_manifest.json`: machine-readable identity and metrics
- `METHODOLOGY.md`: training, labels, features, and evaluation design
- `USAGE_PROTOCOL.md`: inference requirements and QC
- `feature_genes.csv`: ordered feature schema
- `examples/predict_h5ad_counts.py`: chunked H5AD inference example

Round 12 is preserved as a lower-FP fallback under
`Integrated_dataset/models/gdT_prediction_classifier/gdtai_v3_round12_vs_round14/round12_model.pkl`.
