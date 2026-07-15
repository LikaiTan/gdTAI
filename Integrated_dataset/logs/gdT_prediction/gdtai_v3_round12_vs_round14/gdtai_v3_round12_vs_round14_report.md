# gdTAI v3 Round 12 vs Round 14

## Decision: Round 14

Round 14 is selected because it is the only candidate that reaches the prespecified 0.80 full-atlas gold-recall floor while remaining below all FP guardrails. It also has substantially higher held-out cord-blood recall and higher mean F1 across the three evaluation frames. Round 12 remains the lower-FP fallback.

## Performance

| Evaluation              | Round 12 precision | Round 12 recall    | Round 12 F1        | Round 14 precision | Round 14 recall    | Round 14 F1        |
| ----------------------- | ------------------ | ------------------ | ------------------ | ------------------ | ------------------ | ------------------ |
| Full atlas primary gold | 0.9822577018214528 | 0.7560405704790917 | 0.8544295131349882 | 0.9752131916944288 | 0.8055940182774854 | 0.8823256210723865 |
| Atlas held-out combined | 0.9924714739442418 | 0.6750680108817411 | 0.8035620743844945 | 0.9853960624826694 | 0.853016482637222  | 0.9144401080756529 |
| Independent external    | 0.9690107270560191 | 0.8575949367088608 | 0.9099048684946839 | 0.9509345794392523 | 0.8586497890295358 | 0.9024390243902439 |

## Guardrails

| model   | full_atlas_recall  | full_atlas_estimated_fp_fraction | external_specificity | heldout_GSE254249_fp_rate | full_atlas_primary_gold_recall_ge_0p80 | full_atlas_estimated_fp_fraction_le_0p05 | external_specificity_ge_0p9985 | heldout_GSE254249_fp_rate_le_0p001 | all_guardrails_pass |
| ------- | ------------------ | -------------------------------- | -------------------- | ------------------------- | -------------------------------------- | ---------------------------------------- | ------------------------------ | ---------------------------------- | ------------------- |
| round12 | 0.7560405704790917 | 0.0246787676545859               | 0.9993166346886745   | 0.00036845771659843596    | False                                  | True                                     | True                           | True                               | False               |
| round14 | 0.8055940182774854 | 0.0406271092956388               | 0.9988961021893973   | 0.0009096299878523784     | True                                   | True                                     | True                           | True                               | True                |

## Full-Atlas Overlap

| total_atlas_cells | round12_predicted | round14_predicted | both   | round12_only | round14_only | union  | jaccard            | fraction_round12_also_round14 | fraction_round14_also_round12 |
| ----------------- | ----------------- | ----------------- | ------ | ------------ | ------------ | ------ | ------------------ | ----------------------------- | ----------------------------- |
| 5128904           | 213182            | 251356            | 203164 | 10018        | 48192        | 261374 | 0.7772923091049607 | 0.9530072895460218            | 0.8082719330352169            |

## Artifact Manifest

| model   | display_name | payload_model                             | threshold | sha256                                                           | size_bytes | snapshot_path                                                                                                                                                        | provenance                                                                                                                                           | n_gene_features | n_total_features | model_object_type    |
| ------- | ------------ | ----------------------------------------- | --------- | ---------------------------------------------------------------- | ---------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- | --------------- | ---------------- | -------------------- |
| round12 | Round 12     | v3_round12_hist_gradient_fixed_0p5        | 0.5       | 7373e79350f7db190c415b376b9763e31652754438ee8c5afd3853beb7b2ebc4 | 808096     | /home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/models/gdT_prediction_classifier/gdtai_v3_round12_vs_round14/round12_model.pkl | /home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/gdTAI_v3_model.pkl | 210             | 226              | PlattCalibratedModel |
| round14 | Round 14     | v3_round14_v2_score_trdc_gate_fixed_0p936 | 0.936     | 16dedc0081da9b8487887341232bcf8c9c9403dd3bbd72e04cab43d4cd7b2e09 | 15532      | /home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/models/gdT_prediction_classifier/gdtai_v3_round12_vs_round14/round14_model.pkl | git:21b5d60:Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/gdTAI_v3_model.pkl                                                        | 210             | 226              | ConditionalGateModel |

HTML: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/gdT_prediction/gdtai_v3_round12_vs_round14/index.html`
PDF: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/gdT_prediction/gdtai_v3_round12_vs_round14/gdtai_v3_round12_vs_round14_report.pdf`
