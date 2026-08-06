# gdTAI Training And Evaluation Audit Summary

Date: 2026-08-06 HKT

## Conclusion

gdTAI contains useful biological signal, but current performance estimates do not rule out selection overfitting. The highest-priority issues are:

1. V2 algorithm selection used the cohort named validation.
2. V3 Round 12 versus Round 14 promotion used the cohort described as independent external.
3. Train/tune splitting was at cell level within source-label strata rather than donor/library/clonotype groups.
4. The external positive label used CD4 annotation/expression exclusion while CD4 is a model feature.
5. Positive/negative labels and held-out components are strongly source-confounded.
6. The V3 pickle, manifest, completed-run summary, and documented training composition do not describe one coherent promoted build record.

## Corrected TCR-only benchmark

The existing external predictions were re-evaluated using 1,046 paired productive gdT/no productive abT positives and 33,117 strict paired productive abT/no productive gdT negatives.

| Model | Precision | Recall | Specificity | F1 | FP |
| --- | ---: | ---: | ---: | ---: | ---: |
| V2 high-purity | 0.9465 | 0.8289 | 0.9985 | 0.8838 | 49 |
| V3 Round 12 | 0.9784 | 0.8212 | 0.9994 | 0.8929 | 19 |
| V3 Round 14 | 0.9634 | 0.8308 | 0.9990 | 0.8922 | 33 |

This benchmark is still reused evidence, not an untouched external test. Round 12 and Round 14 are effectively tied on F1; Round 12 is purer and Round 14 is slightly more sensitive.

## Coefficients

The V2 standardized logistic coefficients have maximum absolute magnitude 1.776. Coefficient magnitude is not the primary overfitting concern.

The V3 payload embeds `accepted_for_promotion=False` while the manifest says promoted; the retained completed-run summary selects `v3_round12_hist_gradient_fixed_0p5`. This is a semantic provenance defect, not a checksum failure.

## Extension challenge

- V2 high-purity pooled strict-NK FPR: 0.7248%
- V3 Round 12 pooled strict-NK FPR: 1.2039%
- V3 Round 14 pooled strict-NK FPR: 1.7147%

These eight extension cohorts contain no positive truth and therefore cannot select the final model.

## Recommended next action

Run a precommitted nested leave-one-dataset-out re-evaluation. Use donor/library/clonotype-grouped inner folds; pure TCR labels independent of model features; dataset-balanced weights; fold-local feature selection/calibration; macro-dataset metrics; realistic-prevalence PPV/FDR; and explicit NK/abT/dropout stress panels. Keep V2, Round 12, and Round 14 frozen as comparators.

## Outputs

- HTML: `gdT_prediction/gdtai_methodology_audit/index.html`
- PDF: `gdT_prediction/gdtai_methodology_audit/gdtai_methodology_audit_report.pdf`
- Independent reviewer record: `Integrated_dataset/logs/gdT_prediction/gdtai_methodology_audit/independent_reviewer_record.md`
- Issue register: `Integrated_dataset/tables/gdT_prediction/gdtai_methodology_audit/issue_register.csv`
- Recommended actions: `Integrated_dataset/tables/gdT_prediction/gdtai_methodology_audit/recommended_actions.csv`
