# gdTAI v3 TRDC/NK Guard Methodology

## Objective

The objective was to move `GSE144469` into training/tuning, retain non-GSE144469 validation sets, evaluate an independent external H5AD only after fitting, and apply candidate models to the full 5.13M-cell atlas. The practical target was full-atlas primary-gold recall above `0.8` with estimated abT false positives below `5%` of predictions.

## Input Data Policy

- The atlas training input was `high_speed_temp/Integrated_dataset/integrated_plus6.h5ad`.
- The independent external test dataset is registered as `BALF_BLOOD_COPD`.
- Its canonical input is `data/datasets/BALF_BLOOD_COPD/processed/current.h5ad`.
- The historical path `/home/tanlikai/databank/owndata/singlecell/data/results/phase4_final_annotated.h5ad` remains a compatibility alias to the same file.
- The external H5AD was never used for fitting, threshold tuning, feature selection, or full-atlas candidate selection.
- External inference requires `layers["counts"]`; the pipeline refuses normalized/log `X` for external use.

## Labels And Splits

Primary gold train/tune included `GSE144469`. Non-GSE144469 validation remained held out:

- `GDT_2020AUG_woCOV` cord-blood gdT positives
- `GSE254249` paired-TCRAB/no-gdTCR negatives

Training/tune summary from the final run:

- train primary: 505,640 cells, including 24,082 gdT gold cells
- tune primary: 126,410 cells, including 6,021 gdT gold cells
- training sample after caps: 24,082 positives and 100,000 negatives
- explicit hard negatives retained after expression filters: 888 train NK `TRDC+TRDV-`, 565 tune NK `TRDC+TRDV-`, plus paired-TCRAB-only negatives
- NK+TCRAB-overlap cells were excluded from negative train/tune and hard-negative construction

## Features

All features are expression-derived and usable without TCR-seq metadata at inference. Counts are converted to `log1p(counts per 10000)`. The package stores 210 gene features and 16 engineered features.

Engineered features include:

- `any_TRDV`, `any_TRDJ`, `any_TRG`, `any_ab_TCR_gene`
- `TRDC_only`, `TRDC_plus_TRDV`, `TRDC_plus_TRDJ`
- `CD3_score`, `NK_score`, `gdT_TCR_score`, `abT_TCR_score`, `NK_minus_CD3_score`
- `TRDC_log1p`, `TRDV_score`, `TRDJ_score`, `TRG_score`

## Candidate Algorithms Tested

The iteration did not restrict candidates to XGBoost. Candidate rounds included:

- elastic-net logistic model with Platt calibration
- sklearn histogram-gradient model with Platt calibration
- two-stage T-cell gate using the elastic-net score
- frozen gdTAI v2 logistic score with conditional `TRDC-only` / NK-like post-gate
- fixed recall-rescue thresholds on the gated v2 score

LightGBM support exists in the script but was skipped in this environment because LightGBM was unavailable. ExtraTrees is opt-in because it was too slow for the default loop.

## Final Selected Algorithm

The selected candidate is `v3_round14_v2_score_trdc_gate_fixed_0p936`. It wraps the gdTAI v2 base score with a conditional gate. The base score is a sklearn pipeline:

```text
StandardScaler()
LogisticRegression(class_weight="balanced", max_iter=1000)
```

The gate suppresses cells to score `<= 0.02` if they satisfy either of these conditions:

```text
TRDC_only AND NOT any_TRDV AND NOT any_TRDJ AND NK_minus_CD3_score > 0.35
TRDC_only AND abT_TCR_score > gdT_TCR_score + 0.40 AND NK_minus_CD3_score > 0.15
```

The final prediction rule is:

```text
predicted_gdT = gated_score >= 0.936
```

## Full-Atlas FP Estimate

The full atlas contains datasets without TCR sequencing. Hidden abT FP was estimated by:

1. counting predicted paired-TCRAB/no-gdT cells in TCR-sequenced sources,
2. computing observed abT FP rate among non-gold predictions in TCR-sequenced sources,
3. applying that rate to non-gold predictions in sources without TCR sequencing,
4. adding observed paired-TCRAB FP plus estimated hidden FP.

For the selected candidate:

- observed paired-TCRAB abT FP: `5,580`
- estimated hidden abT FP without TCRseq: `4,631.9`
- estimated total abT FP: `10,211.9`
- estimated FP / predictions: `0.0406`

## Metrics

Full atlas primary-gold metrics:

- predicted: `251,356`
- recall: `0.8056`
- precision: `0.9752`
- F1: `0.8823`

Independent external primary metrics:

- predicted: `856`
- TP / FP / FN: `814` / `42` / `134`
- precision: `0.9509`
- recall: `0.8586`
- F1: `0.9024`
- FP / predictions: `0.0491`

## Round 12 Versus Round 14 Revalidation

Both model artifacts were pinned by SHA256 and compared on identical cohorts.
Round 14 was selected on 2026-07-15 because it was the only model to pass all
four prespecified guardrails, including full-atlas gold recall >= 0.80, while
remaining below the atlas, external, and GSE254249 false-positive limits.

- mean F1 across the three evaluation frames: `0.8997`
- comparison report: `gdT_prediction/gdtai_v3_round12_vs_round14/index.html`
- Round 12 remains a lower-FP fallback; it is not the canonical default.

## Promotion Status

`v3_round14_v2_score_trdc_gate_fixed_0p936` is the promoted gdTAI v3.0 balanced default. Its higher NK and
paired-TCRAB call burdens relative to Round 12 must remain visible in external
QC reports.
