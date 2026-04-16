# plus6 γδ-focused statistics

- tissue column used: `tissue_corrected`
- total cells: `5,128,904`
- Sorted_gdT = True: `48,879`
- paired TRG/TRD but not doublet: `33,987`
- TRD - TRAB > `0.4` and no productive TRA/TRB and not NK: `182,796`
- cells meeting at least one of the three criteria: `206,524`
- gdT-cell annotations: `19,705`
- gdT-cell annotations with paired TRG/TRD: `8,067`
- gdT-cell annotations meeting at least one of the three criteria: `15,448`

## Overlap breakdown

| category                                      | n_cells |
| --------------------------------------------- | ------- |
| sorted_gdT_only                               | 15840   |
| paired_TRG_TRD_not_doublet_only               | 2970    |
| TRD_minus_TRAB_gt_0p4_only                    | 147301  |
| sorted_gdT_and_paired_TRG_TRD_only            | 4918    |
| sorted_gdT_and_TRD_minus_TRAB_gt_0p4_only     | 9396    |
| paired_TRG_TRD_and_TRD_minus_TRAB_gt_0p4_only | 7374    |
| all_three_criteria                            | 18725   |

## gdT cells with paired gdTCR by tissue

| tissue                                                     | gdt_cells | gdt_with_paired_TRG_TRD | fraction_paired_TRG_TRD_within_gdt |
| ---------------------------------------------------------- | --------- | ----------------------- | ---------------------------------- |
| cord blood                                                 | 10397     | 5129                    | 0.4933153794363759                 |
| peripheral blood                                           | 4772      | 2929                    | 0.613788767812238                  |
| Descending, sigmoidal colon; rectum                        | 17        | 4                       | 0.23529411764705882                |
| pancreas                                                   | 58        | 3                       | 0.05172413793103448                |
| duodenum                                                   | 95        | 1                       | 0.010526315789473684               |
| lung                                                       | 22        | 1                       | 0.045454545454545456               |
| thymus                                                     | 2127      | 0                       | 0.0                                |
| blood                                                      | 978       | 0                       | 0.0                                |
| tumor                                                      | 351       | 0                       | 0.0                                |
| csf                                                        | 307       | 0                       | 0.0                                |
| heart                                                      | 115       | 0                       | 0.0                                |
| spleen                                                     | 89        | 0                       | 0.0                                |
| Tumor                                                      | 65        | 0                       | 0.0                                |
| jejunum                                                    | 40        | 0                       | 0.0                                |
| bone marrow                                                | 36        | 0                       | 0.0                                |
| Normal                                                     | 30        | 0                       | 0.0                                |
| stomach                                                    | 29        | 0                       | 0.0                                |
| esophagus                                                  | 28        | 0                       | 0.0                                |
| bone_marrow                                                | 20        | 0                       | 0.0                                |
| sigmoid colon                                              | 20        | 0                       | 0.0                                |
| kidney                                                     | 15        | 0                       | 0.0                                |
| rectum                                                     | 15        | 0                       | 0.0                                |
| liver                                                      | 14        | 0                       | 0.0                                |
| skin                                                       | 11        | 0                       | 0.0                                |
| unknown                                                    | 11        | 0                       | 0.0                                |
| liver_tumor                                                | 6         | 0                       | 0.0                                |
| muscle                                                     | 5         | 0                       | 0.0                                |
| Ascending, transverse, descending, sigmoidal colon; rectum | 4         | 0                       | 0.0                                |
| balf                                                       | 4         | 0                       | 0.0                                |
| mixed_duodenum_pbmc                                        | 4         | 0                       | 0.0                                |

## gdT cells meeting at least one of the three criteria by tissue

| tissue                              | gdt_cells | gdt_meeting_any_three | fraction_meeting_any_three_within_gdt |
| ----------------------------------- | --------- | --------------------- | ------------------------------------- |
| cord blood                          | 10397     | 10397                 | 1.0                                   |
| peripheral blood                    | 4772      | 4582                  | 0.960184409052808                     |
| blood                               | 978       | 189                   | 0.19325153374233128                   |
| tumor                               | 351       | 36                    | 0.10256410256410256                   |
| thymus                              | 2127      | 31                    | 0.01457451810061119                   |
| duodenum                            | 95        | 29                    | 0.30526315789473685                   |
| csf                                 | 307       | 28                    | 0.09120521172638436                   |
| heart                               | 115       | 23                    | 0.2                                   |
| bone marrow                         | 36        | 16                    | 0.4444444444444444                    |
| spleen                              | 89        | 15                    | 0.16853932584269662                   |
| pancreas                            | 58        | 15                    | 0.25862068965517243                   |
| jejunum                             | 40        | 13                    | 0.325                                 |
| Tumor                               | 65        | 7                     | 0.1076923076923077                    |
| esophagus                           | 28        | 7                     | 0.25                                  |
| bone_marrow                         | 20        | 7                     | 0.35                                  |
| Normal                              | 30        | 6                     | 0.2                                   |
| lung                                | 22        | 6                     | 0.2727272727272727                    |
| stomach                             | 29        | 5                     | 0.1724137931034483                    |
| skin                                | 11        | 5                     | 0.45454545454545453                   |
| sigmoid colon                       | 20        | 3                     | 0.15                                  |
| kidney                              | 15        | 3                     | 0.2                                   |
| rectum                              | 15        | 3                     | 0.2                                   |
| liver                               | 14        | 3                     | 0.21428571428571427                   |
| muscle                              | 5         | 3                     | 0.6                                   |
| bladder                             | 3         | 3                     | 1.0                                   |
| Descending, sigmoidal colon; rectum | 17        | 2                     | 0.11764705882352941                   |
| unknown                             | 11        | 2                     | 0.18181818181818182                   |
| balf                                | 4         | 2                     | 0.5                                   |
| mixed_duodenum_pbmc                 | 4         | 2                     | 0.5                                   |
| lymph node                          | 2         | 2                     | 1.0                                   |

