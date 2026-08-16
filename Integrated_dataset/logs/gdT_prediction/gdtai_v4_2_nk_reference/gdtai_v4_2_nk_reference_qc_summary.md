# gdTAI V4.2 integration and pseudo-NK consensus QC

## Result

`FAIL_NO_PSEUDO_NK`: technical execution passed, but zero of
113,287 eligible cells met the frozen pseudo-NK
consensus contract. No classifier was fitted.

Six clusters passed the 95% anchor-purity and 2% productive-TCR-contamination
limits, but each had 86.9%
to 89.7% of eligible
candidate cells from one source, above the frozen 70% cap.

## Artifacts

- HTML: `gdT_prediction/gdtai_v4_2_nk_reference/index.html`
- PDF: `gdT_prediction/gdtai_v4_2_nk_reference/gdtai_v4_2_nk_reference_qc_report.pdf`
- Summary JSON: `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_nk_reference/gdtai_v4_2_nk_reference_qc_summary.json`
- Canonical tables: `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_nk_reference/`
- Canonical figures: `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_nk_reference/`

## Decision

Do not fit the planned V4.2 classifier from this empty pseudo-label lane. The
next iteration must repair source balance and boundary localization without
using locked-cohort outcomes to tune the pseudo-label contract.
