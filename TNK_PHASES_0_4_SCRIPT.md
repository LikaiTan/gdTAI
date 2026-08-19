# TNK Phase 0-4 Canonical Workflow

## Purpose

This file is the canonical executable workflow for:

- Phase 0: dataset audit and eligibility triage
- Phase 1: coarse T/NK extraction
- Phase 1b: conservative first cleanup
- Phase 1c: merged metadata backup and replacement
- Phase 2: merged cleanup
- Phase 3: scVI integration
- Phase 4: TRAB/TRB/TRD scoring

This file describes the standard workflow only.
Rules and active exceptions live in `TNK_PIPELINE_RUNBOOK.md`.
Current milestone state lives in `TNK_PIPELINE_STATUS.md`.

## Canonical environment

Use:

- `/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python`

## Canonical inputs and outputs

Inputs:

- `configs/datasets/integration_inputs.csv`
- dataset files selected in `configs/datasets/datasets.csv`

Canonical output root:

- `Integrated_dataset/`

Canonical milestone H5AD files:

- `Integrated_dataset/TNK_candidates.h5ad`
- `Integrated_dataset/TNK_cleaned.h5ad`
- `Integrated_dataset/integrated.h5ad`

Allowed supplementary milestone exception:

- `Integrated_dataset/TNK_candidates_supp.h5ad`
  - used only for the approved supplementary 10x 5' intake lane

## Shared references

- TCR integration:
  - `TCR_INTEGRATION_SOP.md`
  - `configs/tcr/integration_workflow.json`
- tissue harmonization:
  - `TISSUE_CORRECTION_WORKFLOW.md`
- disease-status harmonization:
  - `disease_status_correction_workflow.py`

## Script recording rule

For each major phase or task, this file should record:

- phase or task name
- exact `.py` script used
- key outputs

## gdTAI V4.2 Step 0 ground-truth freeze

Objective:

- construct expression-independent receptor truth with source-aware UMI rules
- repair GSE144469 CD3-library receptor evidence from raw contigs
- curate multi-dataset author-defined NK negatives
- freeze development and broad lockbox roles before classifier fitting

Exact `.py` script:

- `workflows/gdtai/run_gdtai_v4_2_ground_truth.py`

Core outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_ground_truth/v4_2_label_manifest.parquet`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_ground_truth/tcr_umi_status_by_source.csv`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_ground_truth/nk_author_curation_audit.csv`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_ground_truth/summary.json`
- `gdT_prediction/gdtai_v4_2_ground_truth/index.html`
- `gdT_prediction/gdtai_v4_2_ground_truth/gdtai_v4_2_ground_truth_report.pdf`

Gate semantics:

- the corrected atlas is read-only and checksum-bound
- UMI=1 paired calls are excluded when quantitative support exists
- genuinely UMI-unavailable published paired calls remain eligible
- silver and dual-receptor cells have zero fitting and threshold-selection weight
- all extension cohorts, BALF_BLOOD_COPD, GDT_2020, and GDTlung remain locked
- no classifier is fitted in this task

## gdTAI V4.2 development-only nested training

Objective:

- evaluate a two-stage T/NK then gdT/alpha-beta classifier under source-held-out
  nested development folds
- enforce frozen balanced and high-purity operating-profile constraints
- fail closed before lockbox access when source transfer is inadequate

Exact `.py` scripts:

- `workflows/gdtai/prepare_gdtai_v4_2_training.py`
- `workflows/gdtai/train_gdtai_v4_2_nested.py`
- `workflows/gdtai/train_gdtai_v4_2_recovery.py`

Core outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_training/nested_stage1_candidates.csv`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_training/nested_stage2_candidates.csv`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_training/nested_outer_metrics.csv`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_recovery/nested_outer_metrics.csv`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_training/nested_summary.json`

Current gate semantics:

- nested development completed but failed cross-source transfer requirements
- recovery probes remain diagnostic and do not replace the frozen primary run
- lockboxes remain unscored and no model is eligible for promotion or release

## gdTAI V4.3 receptor-isolated development and lockbox freeze

Objective:

- prevent shared cytotoxic/NK expression from acting as direct or indirect
  gdT evidence
- evaluate a high-recall T-lineage support gate followed by an individual-TCR
  receptor classifier under grouped source-heldout development folds
- freeze independent receptor-gold and author-NK lockbox membership before any
  final model scoring

Exact `.py` scripts:

- `workflows/gdtai/train_gdtai_v4_2_nested.py`
- `workflows/gdtai/freeze_gdtai_v4_3_lockbox.py`

Configuration:

- `configs/models/gdtai/v4_3_rescue_training.json`

Core outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_3_receptor_isolated/effective_feature_contract.csv`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_3_receptor_isolated/smoke_outer_metrics.csv`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_3_lockbox/lockbox_manifest.parquet`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_3_lockbox/freeze_summary.json`

Gate semantics:

- `KLRD1`, generic cytotoxic genes, and NK-associated adaptor/receptor genes
  are excluded from both model stages
- Stage 2 receives no continuous Stage 1 probability
- lockbox membership is checksum-bound and was frozen with `model_scored=false`
- lockbox results must not alter features, thresholds, or rescue rules

## gdTAI V4.3 final fit and common lockbox evaluation

Objective:

- freeze one final model and operating threshold using grouped development OOF
  predictions only
- recompute V4.3, V3 balanced, V2 high-F1, and V2 high-purity predictions from
  the same raw-count matrix on identical frozen lockbox cells
- quantify source-level recall, paired-alpha-beta FPR, author-NK FPR,
  prevalence-aware precision, donor-cluster bootstrap uncertainty, and model
  ranking before any promotion decision

Exact `.py` scripts:

- `workflows/gdtai/finalize_gdtai_v4_3_development.py`
- `workflows/gdtai/evaluate_gdtai_v4_3_common_lockbox.py`
- `workflows/gdtai/render_gdtai_v4_3_final_evaluation.py`

Core outputs:

- `Integrated_dataset/models/gdT_prediction/gdtai_v4_3_receptor_isolated/final_development/model_contract.json`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_3_common_lockbox/overall_metrics.csv`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_3_common_lockbox/per_source_metrics.csv`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_3_common_lockbox/cluster_bootstrap_summary.csv`
- `gdT_prediction/gdtai_v4_3_final_evaluation/index.html`
- `gdT_prediction/gdtai_v4_3_final_evaluation/gdtai_v4_3_final_evaluation_report.pdf`

Final gate semantics:

- V4.3 failed the common lockbox superiority gate and is not promotable
- the lockbox is consumed and must not be reused for threshold or feature
  selection
- V3 balanced remains the default; no V4 model may be released or applied to
  the atlas from this experiment

## gdTAI V4.3 post-lockbox recall-failure investigation

Objective:

- audit why V4.3 loses external recall despite BALF being absent from all model
  coefficient fits
- separate continuous ranking, threshold calibration, receptor-expression
  dropout, productive V-gene bias, QC, and source-specific paired-alpha-beta
  false-positive behavior
- preserve every post-lockbox threshold as diagnostic-only

Exact `.py` script:

- `workflows/gdtai/investigate_gdtai_v4_3_recall_failure.py`

Core outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_3_recall_failure/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v4_3_recall_failure/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_3_recall_failure/summary.json`
- `gdT_prediction/gdtai_v4_3_recall_failure/index.html`
- `gdT_prediction/gdtai_v4_3_recall_failure/gdtai_v4_3_recall_failure_report.pdf`

Final semantics:

- V4.3 has strong BALF ranking but a non-transferable frozen cutoff and no
  overall deployment advantage over V3
- V4.3 improves author-NK FPR but has higher paired-alpha-beta FPR in seven of
  eight extension datasets and substantially lower external gdT recall
- no diagnostic threshold may be promoted, released, or used to claim external
  performance

## gdTAI V4.4 validation-only dual operating modes

Objective:

- correct the V4.3 threshold/deployment-score mismatch by freezing the deployed
  scorer, Platt calibrators, and thresholds together
- add symmetric receptor evidence and curated NK negatives without using shared
  NK/cytotoxic genes or a hard TRDV-expression cutoff
- publish highest-F1 and high-purity operating points selected only from a
  group-disjoint development threshold-validation partition
- compare the frozen modes on the already-consumed common lockbox as diagnostic
  evidence only, without retuning or promotion

Exact `.py` scripts:

- `workflows/gdtai/train_gdtai_v4_4_dual_mode.py`
- `workflows/gdtai/evaluate_gdtai_v4_4_reused_lockbox.py`
- `workflows/gdtai/render_gdtai_v4_4_dual_mode_report.py`

Configuration:

- `configs/models/gdtai/v4_4_dual_mode_training.json`

Core outputs:

- `Integrated_dataset/models/gdT_prediction/gdtai_v4_4_dual_mode/model_contract.json`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_4_dual_mode/threshold_validation_metrics.csv`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_4_dual_mode/source_holdout_metrics.csv`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_4_reused_lockbox/overall_metrics.csv`
- `gdT_prediction/gdtai_v4_4_dual_mode/index.html`
- `gdT_prediction/gdtai_v4_4_dual_mode/gdtai_v4_4_dual_mode_report.pdf`

Final semantics:

- highest-F1 maximizes development-validation F1 and high-purity maximizes
  development-validation F0.5; neither mode has a hard FPR constraint
- test cells cannot fit models/calibrators or select features, model families,
  or thresholds
- V4.4 improves specificity over V4.3 but remains materially below V2/V3 in
  reused-lockbox recall and F1
- V4.4 is not promoted or applied to the atlas; a new untouched final test is
  required for any future promotion claim

## Phase 0: Dataset audit and eligibility triage

Objective:

- inspect every registry dataset before extraction

Phase or task:

- Phase 0 dataset audit

Exact `.py` script:

- `workflows/integration/phase0_dataset_audit.py`

Related helper when rescue is needed:

- `workflows/intake/repair_h5ad_from_raw.py`

Core checks:

- matrix state
- raw-count eligibility
- metadata fields
- duplicated gene names
- donor/sample/library candidate columns
- presence of key TCR and lineage genes

Core outputs:

- `Integrated_dataset/tables/phase0_dataset_audit.csv`
- `Integrated_dataset/tables/phase0_category_summary.csv`
- `Integrated_dataset/logs/phase0_qc_summary.md`
- Phase 0 PNG QC figures

## Extension cohort standalone Phase 0 lane

Objective:

- stage and build approved new GEX/TCR cohorts without merging them
- enforce disease-aware tumor specimen context without overwriting source fields
- stop at a user-reviewed standalone Phase 0 gate

Phase or task:

- extension cohort intake, build, and standalone Phase 0 review

Exact scripts:

- `workflows/intake/validate_extension_cohorts.py`
- `workflows/intake/stage_extension_cohorts.py`
- `workflows/intake/build_extension_h5ads.py`
- `workflows/intake/qc_extension_h5ads.py`
- `workflows/intake/audit_extension_tcr_source_chains.py`
- `workflows/metadata/review_extension_sample_sources.py`

Core inputs:

- `configs/datasets/extension_cohorts.csv`
- `configs/datasets/extension_libraries.csv`
- `configs/metadata/sample_source_refinement_rules.json`

Core outputs:

- `data/interim/extension_intake/built_h5ads_manifest.csv`
- `data/interim/extension_intake/built_h5ads/<cohort_id>.h5ad`
- `Integrated_dataset/logs/extension_intake/extension_phase0_qc_summary.md`
- `Integrated_dataset/tables/extension_intake/extension_phase0_cohort_summary.csv`
- `Integrated_dataset/figures/extension_intake/extension_phase0_cohort_overview.png`
- `Integrated_dataset/tables/extension_intake/extension_tcr_source_chain_audit.csv`
- `Integrated_dataset/logs/sample_source_refinement/extension_sample_source_review.md`
- `Integrated_dataset/tables/sample_source_refinement/extensions/extension_tumor_project_context_counts.csv`

Current gate semantics:

- all H5AD review is read-only
- the user-selected GSE169246 processed `TNK_cleaned.h5ad` is registered as a
  standalone input and has passed the extension T/NK filter audit; it remains
  inactive before merge or integration pending user review
- Tan et al. 2021 remains excluded because `GDT_2020AUG_woCOV` already
  represents that study
- no merge or integration may proceed without explicit user approval

### GSE169246 processed-source registration

Phase or task:

- select and register the user-specified processed GSE169246 T/NK object

Exact `.py` script:

- `workflows/maintenance/manage_dataset_registry.py`

Key input:

- `/home/leeaaron/Exhausted_GDT/TNBC_validation/Integrated_dataset/TNK_cleaned.h5ad`

Key outputs:

- `data/datasets/GSE169246/processed/artifacts/TNK_cleaned.h5ad`
- `data/datasets/GSE169246/processed/current.h5ad`
- `docs/GSE169246_PROCESSED_H5AD_INTAKE.md`

The selected source passed standalone metadata/TCR harmonization and filter QC.
It remains inactive until user review is complete.

## Extension cohort standalone Phase 1 T/NK filter

Objective:

- filter every new cohort independently to a high-recall T/NK candidate set
- preserve NK cells for downstream gdTAI false-positive evaluation
- stop at a user-reviewed standalone gate before merge or integration

Phase or task:

- extension cohort T/NK filtering and metadata audit

Exact `.py` script:

- `workflows/intake/filter_extension_tnk_cells.py`

Core inputs:

- `data/interim/extension_intake/built_h5ads/<cohort_id>.h5ad`
- `data/datasets/GSE169246/processed/current.h5ad`

Selection rules:

- retain every cell with productive TRA, TRB, TRG, or TRD evidence
- retain author-annotated T or NK cells
- otherwise require a core CD3/TCR-alpha-beta marker, a gamma-delta marker,
  or the canonical multi-marker NK rule
- do not accept `IL7R` or `LTB` alone as T-lineage evidence
- explicit non-T/NK annotations block marker-only retention but cannot remove
  a productive-TCR cell

Metadata rules:

- preserve exact RNA-series provenance in `source_gse_id`
- record the intake bundle separately in `extension_cohort_id`
- require complete sample, library, donor, tissue, and specimen-context fields
- require unique cell IDs and `library_id + barcode` keys
- recompute and validate canonical TCR logical flags from CDR3 fields
- distinguish the 78 GSE169246 source libraries from 51 biological samples
  using the preserved full-barcode suffix

Core outputs:

- `data/interim/extension_intake/tnk_filtered_h5ads_manifest.csv`
- `data/interim/extension_intake/tnk_filtered_h5ads/<cohort_id>.h5ad`
- `Integrated_dataset/logs/extension_intake/extension_tnk_filter_qc.md`
- `Integrated_dataset/logs/extension_intake/extension_tnk_filter_qc.json`
- `Integrated_dataset/tables/extension_intake/tnk_filter/`
- `Integrated_dataset/figures/extension_intake/tnk_filter/extension_tnk_filter_overview.png`

Gate semantics:

- outputs remain cohort-separated and inactive for integration
- explicit user approval is required before merge, integration, or gdTAI
  evaluation

## Extension official-GEO metadata reconciliation

Objective:

- compare unresolved or over-specific extension metadata against official GEO
  series, sample, matrix, SOFT, and supplementary records
- preserve local source fields while recording separately supported values and
  unresolved ambiguities

Phase or task:

- read-only extension metadata evidence reconciliation

Exact `.py` script:

- `workflows/metadata/audit_extension_geo_metadata.py`

Core input:

- `configs/metadata/extension_geo_metadata_reconciliation.csv`
- `data/interim/extension_intake/tnk_filtered_h5ads_manifest.csv`

Core outputs:

- `docs/EXTENSION_GEO_METADATA_RECONCILIATION.md`
- `Integrated_dataset/tables/geo_metadata_reconciliation/`
- `Integrated_dataset/logs/geo_metadata_reconciliation/extension_geo_metadata_reconciliation_audit.md`

Standard behavior:

- use only direct official NCBI GEO evidence URLs
- validate all declared local values against backed read-only H5AD metadata
- preserve GEO-supported, GEO-partial, and GEO-unavailable states separately
- report scoped affected-cell counts without changing source metadata
- verify source-accession coverage and unchanged H5AD size/mtime during audit
- require additive derived columns for any later correction; never overwrite
  original metadata

Current gate semantics:

- the audit is complete and read-only
- GSE169246 `_b` libraries require a reviewed additive blood correction before
  any extension merge or integration
- metadata writeback, merge, and integration remain blocked pending explicit
  user approval

## Extension frozen gdTAI negative-control screen

Objective:

- apply every registered V2/V3 operating profile to the independently filtered
  extension cohorts without mutating or merging H5AD files
- quantify alpha-beta T-cell and NK false positives without treating no-TCR
  candidates as positive truth

Phase or task:

- checksum-pinned frozen-profile extension-cohort negative-control evaluation

Exact `.py` script:

- `workflows/gdtai/compare_frozen_gdtai_profiles.py`

Execution command:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/compare_frozen_gdtai_profiles.py \
  --stage screen --chunk-size 50000
```

Core inputs:

- `data/interim/extension_intake/tnk_filtered_h5ads_manifest.csv`
- `configs/models/gdtai/extension_evaluation.json`
- `configs/models/gdtai/model_registry.csv`

Core outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_extension_evaluation/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_extension_evaluation/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_extension_evaluation/`
- `gdT_prediction/gdtai_extension_evaluation/index.html`
- `gdT_prediction/gdtai_extension_evaluation/gdtai_extension_evaluation_report.pdf`

Standard behavior:

- verify every model artifact and input H5AD against its registered SHA-256
- use `log1p(raw counts per 10,000)` in bounded chunks
- accept productive single TRA or single TRB evidence as an alpha-beta
  negative control when productive TRG/TRD evidence is absent
- report paired and single-chain alpha-beta controls separately
- define strict NK controls by NK annotation/expression plus absence of all
  productive TRA/TRB/TRG/TRD evidence and doublet flags
- report no-productive-TCR calls as unevaluable candidates, not true positives
- calculate per-cohort, per-source, per-library, and per-donor rates with
  Wilson 95% intervals and descriptive exact paired tests
- run a feature-coverage sensitivity analysis that separates schema-warning
  cohorts such as GSE169246
- do not calculate new-cohort recall/F1/AUC or select/promote a model when the
  screen contains no unbiased gdT positive truth cells
- keep merge and integration blocked after the screen

## Post-Phase-4 gdTAI methodology audit

Objective:

- independently review data preparation, ground truth, features, algorithms,
  splitting, tuning, evaluation, deployment claims, and release provenance
- re-evaluate retained cross-study predictions with expression-independent
  TCR-only labels
- define the strongest defensible next program without adding data

Phase or task:

- read-only independent gdTAI training and evaluation audit

Exact `.py` script:

- `workflows/gdtai/build_gdtai_methodology_audit.py`

Core outputs:

- `docs/GDTAI_METHODOLOGY_AUDIT.md`
- `Integrated_dataset/tables/gdT_prediction/gdtai_methodology_audit/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_methodology_audit/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_methodology_audit/`
- `gdT_prediction/gdtai_methodology_audit/index.html`
- `gdT_prediction/gdtai_methodology_audit/gdtai_methodology_audit_report.pdf`

Standard behavior:

- do not open or mutate an H5AD
- verify registered model checksums and semantic agreement among pickle,
  manifest, registry, and run summary
- treat all inspected cohorts as development or stress benchmarks
- use paired productive TCR-only labels for the corrected primary benchmark;
  keep CD4/NK/low-CD3 expression as stress strata rather than truth rules
- require nested leave-one-dataset-out evaluation with donor/library/clonotype
  grouping before the next model-selection or promotion decision

## Post-Phase-4 gdTAI V4 protocol freeze

Objective:

- freeze the gdTAI V4 labels, features, nested grouped evaluation, operating
  points, promotion rules, and supervision gates before any fitting

Phase or task:

- gdTAI V4 Step 0 precommitted training and validation protocol

Exact `.py` script:

- none; this is a documentation-only supervision gate rendered with Pandoc and
  headless Chrome

Core outputs:

- `docs/GDTAI_V4_PRECOMMITTED_PLAN.md`
- `gdT_prediction/gdtai_v4_precommit/index.html`
- `gdT_prediction/gdtai_v4_precommit/gdtai_v4_precommitted_plan.pdf`

Gate semantics:

- no model fitting, calibration, threshold search, whole-atlas inference, or
  promotion is authorized by Step 0
- protocol v1.2 permits the registered HRA005041 `log1p(CP10K)` matrix only
  after a full per-cell transformed-state audit and makes all sorted cohorts
  sensitivity-only; this amendment was frozen before fitting
- Step 1 is limited to input/label/feature preflight, a per-source CD4/Treg
  exclusion-imposed recall-ceiling audit, and checksum-pinned grouped split
  construction after explicit supervision approval
- CD4/Treg cutoffs cannot be changed after Step 2 begins; any later change
  requires a new precommitted protocol version and fresh preflight

## Post-Phase-4 gdTAI V4 Step 1 preflight

Objective:

- verify the frozen V4 input, label, feature, exclusion, and grouped-split
  contract without fitting a model

Phase or task:

- gdTAI V4 Step 1 preflight and grouped split freeze

Exact `.py` script:

- `workflows/gdtai/run_gdtai_v4_preflight.py`

Execution command:

```bash
MPLCONFIGDIR=/tmp/matplotlib-gdtai-v4 \
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/run_gdtai_v4_preflight.py --log-level INFO
```

Core configuration:

- `configs/models/gdtai/v4_preflight.json`

Core outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_preflight/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v4_preflight/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_preflight/`
- `gdT_prediction/gdtai_v4_preflight/index.html`
- `gdT_prediction/gdtai_v4_preflight/gdtai_v4_preflight_report.pdf`

Current gate result:

- protocol v1.2 Step 1 completed with `PASS`, zero failures, and zero warnings
- the registered HRA005041 transformed matrix passed the frozen full-cell
  expression audit, and all sorted cohorts were verified as sensitivity-only
- the fixed CD4/Treg recall-ceiling audit passed at `98.38%` source-macro
  recall; no cutoff was tuned
- no model fitting, calibration, threshold search, inference, or promotion ran
- Step 2 remains unauthorized until explicit supervision approval

## Post-Phase-4 gdTAI V4 Step 2 nested evaluation

Objective:

- run the checksum-bound nested leave-one-dataset-out evaluation after the
  Step 1 PASS package has received explicit supervision approval
- compare V4 against the legacy TRD-minus-TRAB strategy and fair compact/V2-like
  models on the same folds, labels, calibration procedure, and operating rules

Phase or task:

- gdTAI V4 Step 2 project-data feature caching, nested evaluation, and report

Exact `.py` scripts:

- `workflows/gdtai/gdtai_v4_nested_core.py`
- `workflows/gdtai/run_gdtai_v4_nested_evaluation.py`

Core configuration and authorization template:

- `configs/models/gdtai/v4_nested_evaluation.json`
- `configs/models/gdtai/v4_step2_approval_template.json`

Execution commands:

```bash
MPLCONFIGDIR=/tmp/matplotlib-gdtai-v4 \
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/run_gdtai_v4_nested_evaluation.py --stage validate

MPLCONFIGDIR=/tmp/matplotlib-gdtai-v4 \
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/run_gdtai_v4_nested_evaluation.py \
  --stage all --candidate-jobs 26 --fold-jobs 3 --matrix-row-chunk 20000
```

Core outputs:

- `Integrated_dataset/cache/gdT_prediction/gdtai_v4_nested_evaluation/`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_nested_evaluation/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v4_nested_evaluation/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_nested_evaluation/`
- `gdT_prediction/gdtai_v4_nested_evaluation/index.html`
- `gdT_prediction/gdtai_v4_nested_evaluation/gdtai_v4_nested_evaluation_report.pdf`

Gate semantics:

- cache construction, fitting, calibration, and threshold selection require a
  checksum-bound `STEP2_APPROVAL.json`; validation alone never fits a model
- CD4/Treg exclusions, labels, folds, 197-gene universe, candidate grids,
  operating-point guardrails, and promotion criteria remain frozen
- all three primary datasets are held out in turn; donor/library/clonotype groups
  remain intact within inner folds
- report held-out paired-abT and strict-NK false-positive rates separately
- VDJ rescue is reported separately and cannot improve RNA-only model metrics
- Step 2 may establish that V4 fails; it does not authorize model promotion,
  release fitting, or whole-atlas inference

Current result:

- HRA005041 and GSE144469 outer folds completed as negative results: no V4 or
  fair comparator threshold passed the frozen operating modes
- Stage-1 SAGA was nonconverged and passed all held-out strict-NK cells in both
  completed folds
- the user retired the CPU-only experiment before BALF_BLOOD_COPD completed;
  no CPU model is promotable and the saved artifacts are diagnostic only

## Post-Phase-4 gdTAI V4.1-GPU plan freeze

Objective:

- replace the retired CPU-only SAGA search with a bounded, deterministic GPU
  nested experiment without changing truth labels, grouped folds, 197 genes,
  CD4/Treg exclusions, operating guardrails, or promotion criteria
- verify GPU access, determinism, serialization, memory, checkpointing, and
  non-interference with the server's existing MPS owner before project fitting

Phase or task:

- gdTAI V4.1-GPU Gate A feasibility and precommitted protocol

Exact `.py` script:

- none for project data; feasibility used synthetic-only inline probes

Core outputs:

- `docs/GDTAI_V4_GPU_PRECOMMITTED_PLAN.md`
- `gdT_prediction/gdtai_v4_gpu_precommit/index.html`
- `gdT_prediction/gdtai_v4_gpu_precommit/gdtai_v4_gpu_precommitted_plan.pdf`

Gate semantics:

- use deterministic PyTorch weighted ridge logistic and GPU XGBoost; exclude
  cuML logistic because its installed solver lacks adequate convergence
  provenance in the full-scale probe
- use a unique `CUDA_MPS_PIPE_DIRECTORY`; never alter or connect to the
  `/tmp/nvidia-mps` daemon owned by another user
- prohibit silent CPU fallback for primary candidate training
- do not implement or fit the project-data GPU evaluator until this plan is
  explicitly approved
- after approval, implementation plus synthetic tests form a separate QC gate
  before any project-data model fit

Gate B implementation:

- exact `.py` scripts:
  - `workflows/gdtai/gdtai_v4_gpu_core.py`
  - `workflows/gdtai/run_gdtai_v4_gpu_evaluation.py`
- automated tests:
  - `tests/test_gdtai_v4_gpu_core.py`
- configuration:
  - `configs/models/gdtai/v4_gpu_nested_evaluation.json`
- key outputs:
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_gpu_nested_evaluation/gdtai_v4_gpu_gate_b_summary.json`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_gpu_nested_evaluation/gdtai_v4_gpu_gate_b_summary.md`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_gpu_nested_evaluation/`
  - `gdT_prediction/gdtai_v4_gpu_gate_b/index.html`
  - `gdT_prediction/gdtai_v4_gpu_gate_b/gdtai_v4_gpu_gate_b_qc.pdf`
- result:
  - 32/32 synthetic and read-only checks passed on the A100
  - no project-data model was fitted
  - the project-fit action fails closed until a separate Gate C record exists

Gate C nested evaluation:

- exact `.py` scripts:
  - `workflows/gdtai/run_gdtai_v4_gpu_nested_project.py`
  - `workflows/gdtai/finalize_gdtai_v4_gpu_gate_c_failure.py`
- key outputs:
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_gpu_nested_evaluation/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_gpu_nested_evaluation/`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_gpu_nested_evaluation/`
  - `Integrated_dataset/models/gdT_prediction/gdtai_v4_gpu_nested_evaluation/checkpoints/`
  - `gdT_prediction/gdtai_v4_gpu_nested_evaluation/index.html`
  - `gdT_prediction/gdtai_v4_gpu_nested_evaluation/gdtai_v4_gpu_gate_c_report.pdf`
- result:
  - all three outer decision paths completed and 99/99 recorded fold fits
    converged
  - 0/15 Stage-1 candidate evaluations passed; no V4 Stage-2 model was eligible
  - no fair comparator or raw legacy-score operating point passed either mode
  - extension-model controls and paired bootstrap were not estimable without
    an eligible balanced V4 model
  - V4.1 failed Gate C; V3 Round 14 remains the promoted balanced default
  - Gate D release fitting, promotion, and atlas inference remain blocked

## Post-Phase-4 gdTAI V4.2 label audit and protocol freeze

Objective:

- determine whether V4.1 Gate C failed because of Stage-1 model capacity or
  heterogeneous NK-control provenance
- freeze a repaired two-stage protocol before any V4.2 fitting

Phase or task:

- gdTAI V4.2 Step 0 read-only NK-label audit and saved-OOF counterfactual

Exact `.py` script:

- `workflows/gdtai/build_gdtai_v4_2_step0_audit.py`

Execution command:

```bash
MPLCONFIGDIR=/tmp/matplotlib-gdtai-v42-step0 \
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/build_gdtai_v4_2_step0_audit.py
```

Core configuration and protocol:

- `configs/models/gdtai/v4_2_precommit.json`
- `docs/GDTAI_V4_2_PRECOMMITTED_PLAN.md`

Key outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_step0/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_step0/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_step0/`
- `gdT_prediction/gdtai_v4_2_step0/index.html`
- `gdT_prediction/gdtai_v4_2_step0/gdtai_v4_2_step0_report.pdf`

Gate semantics:

- the audit may read checksum-bound V4.1 OOF probabilities but cannot fit,
  calibrate, select, promote, or apply a V4.2 model
- expression-based T-lineage coherence is diagnostic only and cannot define NK
  truth
- V4.2 implementation and fitting remain blocked until the Step 0 report and
  frozen protocol receive explicit supervision approval

## Post-Phase-4 gdTAI V4.2 modeling-integration preflight

Objective:

- freeze whole-cohort development and locked-evaluation roles before any new
  data integration
- prove raw-count, metadata, TCR-anchor, feature, GPU, RAM, SSD, and rollback
  compatibility without merging or fitting

Phase or task:

- gdTAI V4.2 read-only sidecar-integration preflight and supervision package

Exact `.py` script:

- `workflows/gdtai/run_gdtai_v4_2_integration_preflight.py`

Execution command:

```bash
MPLCONFIGDIR=/tmp/matplotlib-gdtai-v42-integration-preflight \
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/run_gdtai_v4_2_integration_preflight.py
```

Core configuration and role manifest:

- `configs/models/gdtai/v4_2_integration_preflight.json`
- `configs/models/gdtai/v4_2_cohort_roles.csv`

Key outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_integration_preflight/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_integration_preflight/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_integration_preflight/`
- `gdT_prediction/gdtai_v4_2_integration_preflight/index.html`
- `gdT_prediction/gdtai_v4_2_integration_preflight/gdtai_v4_2_integration_preflight_report.pdf`

Current result:

- `PASS_REVIEW_REQUIRED`; all 59 checks passed
- 4,023,462 development cells and 439,979 locked cells were checksum-bound
- all 21,054 primary dual-annotation NK anchors mapped exactly to the current
  atlas, and all locked cohorts retained 50/50 Stage-1 genes
- the common development universe contains 14,265 genes, of which 13,975
  remain after the frozen TCR/mitochondrial/ribosomal/immunoglobulin exclusion
- conservative peak RAM is 223.6 GiB, SSD free space is 875.1 GiB, and the
  A100 80-GB GPU passed the resource gate
- all nine input H5AD size and nanosecond modification-time pairs were unchanged

Gate semantics:

- no merge, scVI fit, clustering, pseudo-label construction, classifier fit,
  threshold search, promotion, or inference ran
- an inactive checksum-bound approval template was generated
- implementation or integration requires explicit supervision approval; later
  classifier fitting requires a separate reviewed approval gate

## Post-Phase-4 gdTAI V4.2 current-atlas recovery preflight

Objective:

- reconstruct the exact raw-count current-atlas input required by the V4.2
  sidecar after loss of the mirrored SSD `integrated.h5ad`
- prove cell, gene, exclusion, metadata, NK-anchor, checksum, and read-only
  equivalence before any project-data execution

Phase or task:

- gdTAI V4.2 read-only current-atlas recovery audit and supervision package

Exact `.py` scripts:

- `workflows/gdtai/run_gdtai_v4_2_recovery_preflight.py`
- `workflows/gdtai/run_gdtai_v4_2_nk_reference_integration.py`

Execution command:

```bash
MPLCONFIGDIR=/tmp/matplotlib-gdtai-v42-recovery \
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/run_gdtai_v4_2_recovery_preflight.py
```

Core configuration:

- `configs/models/gdtai/v4_2_integration_execution.json`

Key outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_recovery_preflight/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_recovery_preflight/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_recovery_preflight/`
- `gdT_prediction/gdtai_v4_2_recovery_preflight/index.html`
- `gdT_prediction/gdtai_v4_2_recovery_preflight/gdtai_v4_2_recovery_preflight_report.pdf`

Current result:

- `PASS_REVIEW_REQUIRED`; all 20 checks passed
- exact recovery uses 3,705,384 raw cells, removes the saved 78 Phase-3
  exclusions, and yields the frozen 3,705,306-cell input with 27,413 genes
- all five required metadata missing/unique-count audits match exactly, and
  all 21,054 primary dual-annotation NK anchors are recovered
- the prior project-data execution approval is invalid because the physical
  input and checksum-bound execution contract changed

Gate semantics:

- no project-data integration, scVI fit, clustering, pseudo-labeling,
  classifier fitting, thresholding, promotion, or inference ran
- the recovery path is permitted only while the original integrated H5AD is
  absent and only when all source, manifest, metadata, and contract hashes
  match
- recovery execution requires explicit activation of the generated
  checksum-bound approval; classifier fitting remains a later separate gate

## Post-Phase-4 gdTAI V4.2 cluster-stage resource amendment

Objective:

- preserve the completed checksum-bound sparse staging and A100 scVI outputs
  after shared SSD capacity fell below the original global floor
- freeze a stage-specific storage floor before RAPIDS clustering without
  changing any scientific input, representation, clustering, or pseudo-label
  parameter

Phase or task:

- gdTAI V4.2 saved-latent cluster-resource preflight and supervision package

Exact `.py` scripts:

- `workflows/gdtai/run_gdtai_v4_2_cluster_resource_preflight.py`
- `workflows/gdtai/run_gdtai_v4_2_nk_reference_integration.py`
- `workflows/gdtai/build_gdtai_v4_2_nk_reference_qc_report.py`

Completed execution commands:

```bash
ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS=1 \
MPLCONFIGDIR=/tmp/matplotlib-gdtai-v42-recovery \
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/run_gdtai_v4_2_nk_reference_integration.py --stage prepare

ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS=1 \
MPLCONFIGDIR=/tmp/matplotlib-gdtai-v42-fit \
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/run_gdtai_v4_2_nk_reference_integration.py --stage fit
```

Approved execution command:

```bash
ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS=1 \
NUMBA_CACHE_DIR=/tmp/numba-gdtai-v42 \
MPLCONFIGDIR=/tmp/matplotlib-gdtai-v42-cluster \
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/run_gdtai_v4_2_nk_reference_integration.py --stage cluster

ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS=1 \
NUMBA_CACHE_DIR=/tmp/numba-gdtai-v42 \
MPLCONFIGDIR=/tmp/matplotlib-gdtai-v42-consensus \
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/run_gdtai_v4_2_nk_reference_integration.py --stage consensus

ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS=1 \
MPLCONFIGDIR=/tmp/matplotlib-gdtai-v42-report \
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/build_gdtai_v4_2_nk_reference_qc_report.py
```

Key outputs:

- `/ssd/tnk_phase3/Integrated_dataset/gdtai_v4_2_nk_reference/development_hvg_counts.h5ad`
- `/ssd/tnk_phase3/Integrated_dataset/gdtai_v4_2_nk_reference/X_scVI.npy`
- `/ssd/tnk_phase3/Integrated_dataset/gdtai_v4_2_nk_reference/scvi_model/`
- `/ssd/tnk_phase3/Integrated_dataset/gdtai_v4_2_nk_reference/cluster_partitions.npz`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_nk_reference/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_nk_reference/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_nk_reference/`
- `gdT_prediction/gdtai_v4_2_nk_reference/index.html`
- `gdT_prediction/gdtai_v4_2_nk_reference/gdtai_v4_2_nk_reference_qc_report.pdf`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_cluster_resource_preflight/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_cluster_resource_preflight/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_cluster_resource_preflight/`
- `gdT_prediction/gdtai_v4_2_cluster_resource_preflight/index.html`
- `gdT_prediction/gdtai_v4_2_cluster_resource_preflight/gdtai_v4_2_cluster_resource_preflight_report.pdf`

Current result:

- sparse staging passed on 4,023,462 cells x 4,000 HVGs with zero locked
  cohorts; A100 scVI passed with a finite 30-dimensional latent representation
- `PASS_REVIEW_REQUIRED`; all 18 resource-amendment checks passed
- the original 300-GiB floor remains frozen for `prepare` and `fit`; the
  proposed floor is 150 GiB for `cluster` and `consensus`
- worst-case partition payload plus reserve is 2.33 GiB, leaving 108.4 GiB
  above the proposed floor and reserve at audit time
- the user approved the checksum-bound resource amendment on 2026-08-16
- nine global and nine boundary RAPIDS runs completed on the A100 with no CPU
  fallback; the partition SHA-256 is
  `90e83ec83986358a015885e02e29d06eacf726d928e02b2e90428d5c09947a63`
- technical execution passed, but scientific QC ended `FAIL_NO_PSEUDO_NK`:
  the boundary contained 98.87% of cells, zero of 396 clusters met every
  frozen criterion, and zero of 113,287 eligible cells were selected
- all six purity/contamination-passing clusters failed only the 70% source cap;
  GSE292700 contributed 86.92%-89.69% of their eligible candidates

Gate semantics:

- no classifier fitting, thresholding, promotion, release fitting, or inference
  ran under this amendment
- the prior approval is invalid after the checksum-bound resource-contract
  change
- activated `CLUSTER_EXECUTION_APPROVAL.json` authorizes only saved-latent
  RAPIDS clustering and pseudo-NK consensus QC

## Post-Phase-4 gdTAI V4.2 NK-definition repair

Objective:

- repair the anchor/cohort confounding exposed by the failed cluster-consensus
  lane
- define conservative NK reference negatives without treating shared
  cytotoxic programs as NK-specific

Phase or task:

- gdTAI V4.2 read-only productive-T-anchor recovery and exact latent plus
  lineage-expression NK consensus

Exact `.py` script:

- `workflows/gdtai/build_gdtai_v4_2_nk_definition_repair.py`

Execution command:

```bash
MPLCONFIGDIR=/tmp/matplotlib-gdtai-v42-nk-repair \
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/build_gdtai_v4_2_nk_definition_repair.py
```

Key outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_nk_definition_repair/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_nk_definition_repair/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_nk_definition_repair/`
- `gdT_prediction/gdtai_v4_2_nk_definition_repair/index.html`
- `gdT_prediction/gdtai_v4_2_nk_definition_repair/gdtai_v4_2_nk_definition_repair_report.pdf`

Current result:

- `PASS_NK_REFERENCE_READY`; 762,280 existing-atlas productive-T anchors were
  restored, yielding 967,149 productive-T anchors and 21,054 primary NK anchors
- exact A100 50-neighbor scoring passed a bit-identical internal repeat; the
  fail-closed latent threshold is `0.98`, with 93.23% held-out primary-NK recall
  and 0.0958% held-out productive-T FPR
- 469 of 113,287 eligible development candidates passed both latent and strict
  lineage-expression evidence across four sources
- the 469-row repair manifest is retained unchanged for audit; subsequent
  cluster and raw-TCR reviews determine whether any rows may enter low-weight
  training under the same 70% maximum-source cap

Gate semantics:

- cytotoxic genes are diagnostic only and cannot establish NK identity
- locked GSE169246 is validation-only and cannot determine any rule or
  threshold
- the output is a reference-label manifest, not a classifier, model release,
  promotion decision, or whole-atlas inference
- no source H5AD is modified, and V4 remains excluded from GitHub publication

## Post-Phase-4 gdTAI V4.2 NK/T-lineage cluster review

Objective:

- visualize T-lineage, NK-lineage, myeloid-context, and shared-cytotoxic marker
  programs on the saved scVI UMAP
- overlay additional raw-count CD3, TCR constant, delta V/J, and gamma V genes
  from the registered source H5ADs and calculate exact all-cell counts
- compare recoverable annotation evidence with a frozen unsupervised Leiden
  partition before refining the NK training reference

Phase or task:

- gdTAI V4.2 read-only feature-plot and unsupervised-cluster review

Exact `.py` script:

- `workflows/gdtai/build_gdtai_v4_2_nk_cluster_review.py`

Execution command:

```bash
MPLCONFIGDIR=/tmp/matplotlib-gdtai-v42-nk-review \
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/build_gdtai_v4_2_nk_cluster_review.py
```

Key outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_nk_cluster_review/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_nk_cluster_review/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_nk_cluster_review/`
- `gdT_prediction/gdtai_v4_2_nk_cluster_review/index.html`
- `gdT_prediction/gdtai_v4_2_nk_cluster_review/gdtai_v4_2_nk_cluster_review_report.pdf`

Current result:

- `PASS_VISUAL_REVIEW_EXTENDED`; cluster 19 contains 96.12% of all primary NK
  anchors and all but 8 current `NK_CONFIDENT` calls
- hold the 461 dual-evidence cluster-19 cells as a provisional review core, not
  training truth: 350 have at least one TRA/TRB constant-chain UMI, 250 have at
  least two, and 178 have at least three
- hold every current call outside cluster 19 and every rescue tier from the
  strict set
- retain cluster-19 single-evidence cells and cluster-9 NK-like cells as
  review-only rescue tiers because source balance or anchor support is
  insufficient for automatic label expansion
- reject adaptor-rich cluster 1 as NK because its myeloid-context program
  exceeds its NK-lineage program

Standard behavior:

- use `log1p(CP10K)` within the fixed 4,000 scVI integration features
- use `log1p(raw UMI)` for the additional 29 CD3/TCR genes because TCR genes
  were deliberately excluded from the scVI HVGs; verify all 29 genes in all six
  source H5ADs before counting
- calculate exact cluster cell and anchor counts on all 4,023,462 development
  cells; use source/evidence-role sampling weights for expression summaries on
  the deterministic 250,000-cell diagnostic sample
- treat `FCER1G/TYROBP` and cytotoxic genes as insufficient by themselves;
  myeloid genes are contamination context only and are not proposed as model
  features
- do not fabricate the missing historical full-cell scANVI annotation; plot
  only annotation evidence recoverable from frozen manifests and repaired
  anchors

Gate semantics:

- this task changes no H5AD, label manifest, model, threshold, or locked cohort
- the provisional 461-cell core requires source/library ambient-RNA, doublet,
  and paired-TCR review before classifier preflight
- classifier fitting, model promotion, release fitting, and atlas inference do
  not occur in this review
- V4 remains excluded from GitHub publication

## Post-Phase-4 gdTAI V4.2 T/NK-restricted reintegration

Objective:

- retain a high-recall T/NK development pool before dimensional reduction
- recompute source-balanced HVGs only after T/NK restriction
- refit scVI and repeated Leiden clustering, then subcluster the broad NK-like
  boundary without using it as automatic NK training truth

Phase or task:

- gdTAI V4.2 read-only T/NK-restricted integration and NK-boundary review

Exact `.py` script:

- `workflows/gdtai/run_gdtai_v4_2_tnk_reintegration.py`

Execution commands:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/run_gdtai_v4_2_tnk_reintegration.py --stage prepare
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/run_gdtai_v4_2_tnk_reintegration.py --stage fit
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/run_gdtai_v4_2_tnk_reintegration.py --stage cluster
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/run_gdtai_v4_2_tnk_reintegration.py --stage boundary
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/run_gdtai_v4_2_tnk_reintegration.py --stage report
```

Key outputs:

- `configs/models/gdtai/v4_2_tnk_reintegration.json`
- SSD sidecar under
  `/ssd/tnk_phase3/Integrated_dataset/gdtai_v4_2_tnk_reintegration/`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_tnk_reintegration/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_tnk_reintegration/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_tnk_reintegration/`
- `gdT_prediction/gdtai_v4_2_tnk_reintegration/index.html`
- `gdT_prediction/gdtai_v4_2_tnk_reintegration/gdtai_v4_2_tnk_reintegration_report.pdf`

Current result:

- retained 3,927,924/4,023,462 cells while preserving all 21,054 primary NK
  and 967,149 repaired productive-T anchors and removing all flagged doublets
- recomputed 4,000 source-balanced HVGs after subsetting, excluding TCR V/J/D
  genes from ranking and forcing a 27-gene T/NK context panel
- completed 20-epoch A100 scVI, nine global RAPIDS Leiden runs, and nine
  second-pass runs within refined clusters 9 and 18
- the resolution-0.4 boundary partition has mean seed ARI 0.962, but harmonized
  productive-chain overlays identify 373,339/475,953 (`78.44%`) cells with TRA
  or TRB and only 83 with TRD
- the earlier clusters 0/3/5 candidate core is rejected because
  190,757/257,569 (`74.06%`) carry productive TRA or TRB; 56.57% of primary NK
  annotation anchors also conflict with productive TRA/TRB

Gate semantics:

- no boundary subcluster is eligible as NK training truth; the former clusters
  0/3/5 candidate-core designation is withdrawn
- productive-chain overlays use nonempty harmonized productive-filtered CDR3
  fields; assay-unavailable chain absence is not biological negative evidence
- do not define NK from `TRDC`, selected delta-V absence, cytotoxicity,
  `FCER1G/TYROBP`, or any single gene
- source/library ambient, doublet, TCR-depth, and paired-chain review remains
  required before classifier preflight
- no source H5AD, label manifest, classifier, threshold, model registry,
  release artifact, atlas inference, or GitHub publication is changed

## Post-Phase-4 gdTAI V4.2 productive-TRA/TRB boundary conflict audit

Objective:

- determine whether the unexpectedly high productive TRA/TRB rate in the
  V4.2 NK-like boundary is explained by ambient one-UMI calls, non-ideal
  transcriptomic representation, unsafe TCR joins, or genuine alpha-beta T
  cells
- keep UMI/read support, source join provenance, transcriptomic evidence, and
  residual biological compatibility as separate audit dimensions

Phase or task:

- gdTAI V4.2 read-only productive-TRA/TRB provenance and conflict audit

Exact `.py` script:

- `workflows/gdtai/audit_gdtai_v4_2_boundary_tcr_conflicts.py`

Execution command:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/audit_gdtai_v4_2_boundary_tcr_conflicts.py
```

Key outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_tcr_conflict_audit/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_tcr_conflict_audit/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_tcr_conflict_audit/`
- `gdT_prediction/gdtai_v4_2_tcr_conflict_audit/index.html`
- `gdT_prediction/gdtai_v4_2_tcr_conflict_audit/gdtai_v4_2_tcr_conflict_audit_report.pdf`

Current result:

- 277,305/373,339 (`74.28%`) productive-AB boundary calls come from
  source-level join-red-flag cohorts
- only 47/3,047 (`1.54%`) UMI-auditable productive-AB cells are supported
  exclusively by one observed UMI, while 370,292 calls lack quantitative UMI
  provenance
- all 27 forced T/NK-context genes entered the 4,000-gene scVI model; cluster
  association with expression class exceeds association with raw productive-AB
  status, so UMAP or Leiden is not the primary detected cause
- 1,508 cells meet a conservative residual AB-compatible audit stratum; it is
  not promoted to a truth label

Standard behavior:

- treat zero UMI/read values as unavailable because legacy importers used zero
  when quantitative support was absent
- flag a source only from documented barcode-only mapping, extensive exact
  paired-receptor reuse across sample/donor keys, or strong author-lineage
  conflict
- preserve source H5ADs and original TCR fields unchanged
- quarantine red-flag source calls from model truth until raw productive VDJ
  contigs are rebuilt with the canonical `sample_id + barcode_core` join
- rerun the T/NK integration and clustering only after repaired joins pass
  source-level validation

## Post-Phase-4 sample-aware productive TCR join rebuild

Objective:

- reconstruct productive TRA, TRB, TRG, and TRD metadata for all 14 flagged
  sources from raw VDJ evidence without mutating source or milestone H5ADs
- preserve per-chain UMI/read support and distinguish unavailable support from
  measured zero
- fail closed when source sample identity, raw-library identity, or RNA join-key
  uniqueness cannot support a deterministic assignment

Phase or task:

- sample-aware four-chain TCR sidecar rebuild, report, and independent
  validation

Exact `.py` scripts:

- `workflows/intake/recover_gse287541_tcr_from_sra.py`
- `workflows/intake/rebuild_flagged_tcr_joins.py`
- `workflows/reporting/build_tcr_join_rebuild_report.py`
- `workflows/intake/validate_tcr_join_rebuild.py`

Execution commands:

```bash
# One-time official ENA run-manifest setup.
mkdir -p data/datasets/GSE287541/raw/geo_sra
curl -L \
  'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1213711&result=read_run&fields=run_accession,sample_accession,experiment_accession,sample_title,library_name,library_strategy,library_source,library_selection,instrument_platform,read_count,base_count,fastq_bytes,fastq_ftp&format=tsv&download=true' \
  -o data/datasets/GSE287541/raw/geo_sra/ena_run_manifest.tsv
# One-time official 10x VDJ reference setup (SHA-256 of archive:
# a7ba0ae81f44e9ee61338417585cc955ed3c75d49577267be72e05735c654065)
curl -L \
  https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz \
  -o /tmp/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
tar -xzf /tmp/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz \
  -C data/references
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/intake/recover_gse287541_tcr_from_sra.py \
  --workers 2 --cores-per-sample 30 --memory-gb-per-sample 80
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/intake/rebuild_flagged_tcr_joins.py --overwrite
MPLCONFIGDIR=/tmp/matplotlib-tcr-report \
  /home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/reporting/build_tcr_join_rebuild_report.py
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/intake/validate_tcr_join_rebuild.py
```

Key outputs:

- `data/datasets/GSE287541/interim/tcr_recovery/`
- `Integrated_dataset/tables/tcr_join_rebuild/`
- `Integrated_dataset/figures/tcr_join_rebuild/`
- `Integrated_dataset/logs/tcr_join_rebuild/`
- `gdT_prediction/gdtai_v4_2_tcr_join_rebuild/index.html`
- `gdT_prediction/gdtai_v4_2_tcr_join_rebuild/tcr_join_rebuild_report.pdf`

Current result:

- all 14 sources passed full or fail-closed partial rebuilding; no source-level
  quarantine remains
- 3,041,871 rows are staged for replacement, including 1,121,858 cells with
  validated productive TCR and 1,789,643 chain calls with observed UMI support
- GSE125527 and GSE228597 passed exact-CDR3/sample-rotation controls, and all
  46 public GSE287541 TCR runs were reconstructed and checksum-validated
- GSE235863 retains 110 duplicate RNA join-key rows as blank,
  TCR-truth-ineligible replacement rows
- independent validation passed 15/15 checks and all 12 focused tests passed;
  the seven-page HTML/PDF report passed visual table-layout review
- all recorded source H5AD size/mtime pairs remained unchanged

Standard behavior:

- recover all 46 public GSE287541 TCR runs sample by sample, validate the
  10x read structure, run Cell Ranger VDJ with `--disable-ui`, retain compact
  contig/metric/checksum outputs, and remove generated SRA/FASTQ/BAM scratch
- require productive, cell-associated, high-confidence raw records with a
  nonempty CDR3 and join only by `sample_id + barcode_core`
- select one call per chain by highest UMI, then reads, full-length status, and
  stable contig id; retain selected-contig provenance and productive-contig
  multiplicity
- retain UMI/read availability flags and null values when quantitative support
  is absent; never reinterpret unavailable support as zero
- use the published patient-remapping table plus tissue for GSE125527, pooled
  library suffixes for GSE228597, and round plus participant visit for
  GSE287541
- require at least 50% raw-key recovery, or for filtered multi-assay objects an
  exact-CDR3 fallback with at least 100 testable and 100 confirmed calls, at
  least 50% agreement, and at least 20-fold enrichment over sample rotation
- separate metadata-replacement eligibility from TCR-truth eligibility so
  ambiguous legacy calls can be cleared without becoming negative truth
- bind any propagation to the staged sidecar SHA-256 and require explicit
  approval, backups, and post-write validation before changing H5AD metadata

## Post-Phase-4 frozen-membership TCR sidecar overlay preflight

Objective:

- determine whether the validated 14-source TCR sidecar can safely replace
  canonical TCR metadata in the rebuilt full atlas without changing atlas cell
  membership or rerunning integration
- quantify repaired productive-TCR cells outside the atlas before explicitly
  accepting or rejecting their rescue

Phase or task:

- read-only TCR sidecar-to-atlas identity, coverage, and control-loss audit

Exact `.py` script:

- `workflows/integration/audit_tcr_sidecar_overlay.py`
- `workflows/integration/apply_tcr_sidecar_to_full_atlas.py`
- `workflows/reporting/build_tcr_sidecar_application_report.py`
- `workflows/gdtai/audit_post_tcr_truth_nk_boundary.py`

Key outputs:

- `Integrated_dataset/tables/tcr_sidecar_overlay_preflight/overlay_join_by_source.csv`
- `Integrated_dataset/tables/tcr_sidecar_overlay_preflight/sidecar_rows_missing_from_atlas.csv`
- `Integrated_dataset/logs/tcr_sidecar_overlay_preflight/overlay_preflight_summary.json`
- `Integrated_dataset/logs/tcr_sidecar_overlay_preflight/overlay_preflight_summary.md`
- `docs/TCR_SIDECAR_OVERLAY_PLAN.md`
- `Integrated_dataset/tables/tcr_sidecar_application/`
- `Integrated_dataset/figures/tcr_sidecar_application/`
- `Integrated_dataset/logs/tcr_sidecar_application/`
- `gdT_prediction/gdtai_v4_2_tcr_sidecar_application/`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_post_tcr_audit/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_post_tcr_audit/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_post_tcr_audit/`
- `gdT_prediction/gdtai_v4_2_post_tcr_audit/`

Current result:

- `PASS_OVERLAY_READY`; all 2,155,409 current atlas rows from affected sources
  match exactly by `source_gse_id + original_cell_id`, with no unmatched atlas
  row and no expression write
- 1,100,981 retained atlas cells have repaired productive alpha-beta TCR and
  933,890 have paired TRA/TRB
- after combining unaffected sources, the frozen atlas is estimated to contain
  2,270,138 productive-alpha-beta cells and 1,938,158 paired-TRA/TRB cells
- 20,875 productive-alpha-beta cells outside the atlas, including 14,495 paired
  cells, are intentionally not rescued; this omits 0.911% of the whole
  productive-alpha-beta pool and 0.742% of the whole paired-alpha-beta pool
- no sidecar-only cell has productive gamma-delta TCR, so frozen membership does
  not discard a sidecar-supported gamma-delta positive
- the separate TCR-corrected candidate passed all 17 post-write checks with
  SHA-256
  `d32c9d2bdb955b12e1eafbed8322f8cb965cf3a225191e612b53f3d3783480d5`
- all affected and unaffected TCR-value checks passed, the 109 atlas-present
  ambiguous rows are fail closed, and corrected whole-atlas totals are
  2,270,138 productive-alpha-beta and 1,938,158 paired-TRA/TRB cells
- the application QC report is
  `gdT_prediction/gdtai_v4_2_tcr_sidecar_application/index.html` with an
  eight-page visually reviewed PDF; the canonical symlink remains unchanged
- the post-repair truth/NK audit leaves 20,922/21,054 frozen primary NK anchors
  free of corrected productive TCR and reduces boundary productive-alpha-beta
  evidence from the prior unsafe-overlay count of 373,339 to 14,228/475,953
- no-silver corrected truth contains 58,822 gdT gold and 1,926,136 abT gold
  cells with zero rule conflicts; the five-page HTML/PDF report passed visual
  review and no model fitting or canonical publication occurred

Standard behavior:

- require the additive tissue/tumor/sample/library metadata gate to pass before
  applying any TRA/TRB/TRG/TRD sidecar values
- keep biological sample identity separate from technical VDJ-library join
  identity, especially for pooled libraries
- join only atlas `source_gse_id + original_cell_id` to sidecar
  `source_gse_id + source_obs_name`; positional, barcode-only, and display
  `sample_id` joins are forbidden
- keep the exact 5,933,312-cell atlas membership, expression, latent space,
  clusters, and UMAP unchanged during metadata repair
- clear stale calls, preserve null UMI/read semantics, and keep ambiguous rows
  fail closed; rebuild downstream TCR truth labels after replacement
- require temporary-file validation for every H5AD transaction; the approved
  metadata-only candidate is complete, while any future canonical symlink
  switch remains a separate high-risk action

## Pre-TCR full-atlas metadata and sample identity correction

Objective:

- harmonize tissue, specimen context, and tumor type while preserving original
  labels
- repair biological sample and technical library names before TCR metadata is
  joined
- prevent pooled VDJ-library identifiers from replacing biological specimens

Exact `.py` scripts:

- `workflows/metadata/audit_full_atlas_metadata_harmonization_v2.py`
- `workflows/metadata/audit_full_atlas_sample_identity_v2.py`
- `workflows/metadata/build_full_atlas_metadata_overlay_v2.py`
- `workflows/metadata/apply_full_atlas_metadata_overlay_v2.py`

Key inputs and outputs:

- `configs/metadata/full_atlas_metadata_harmonization_v2.json`
- `docs/FULL_ATLAS_METADATA_HARMONIZATION_V2_REVIEW.md`
- `Integrated_dataset/tables/metadata_harmonization/full_atlas_v2/`

Current result:

- tissue/tumor rule review passed 11/11 read-only checks
- sample identity preflight passed 13/13 read-only checks
- GSE125527 recovers 30 source-derived donor-by-tissue samples and GSE254249
  recovers 92 source `Ident` samples
- GSE228597 retains 4,611 explicitly unresolved biological sample IDs while
  recovering their separate technical TCR-library identities
- the complete 5,933,312-row additive overlay passed 18/18 checks and reduced
  unresolved tissue/context to 39,920 cells: 36,800 source-blank GSE206325
  cells and 3,120 inseparable mixed-duodenum/PBMC GSE252762 cells
- the overlay SHA-256 is
  `4da4eea32e9d275de790775e2a0f59d4f6553d72756e9c8c6935f35bb398984f`
- the approved metadata-only candidate write passed all 19 post-write checks
- the validated candidate is
  `/ssd/tnk_phase3/Integrated_dataset/full_atlas/metadata_corrected/integrated_full_atlas.h5ad`
  with SHA-256
  `7f1c5e1cac1074a8e2863703bc1862e225defc5ba1a3adbaabd3f6e023d5871c`
- all 13 additive columns reproduce keyed overlay hashes and value counts;
  sparse X, embeddings, existing obs/var signatures, cell order, and legacy
  TCR fields remain unchanged
- the original canonical atlas remains the byte-identical rollback object;
  no TCR chain call was applied

Gate behavior:

- canonical publication remains a separate transaction; the canonical symlink
  was not changed
- the TCR sidecar application is the next transaction and must use the
  checksum-bound metadata candidate above as input

## Phase 1: Coarse T/NK extraction

Objective:

- build a high-recall merged T/NK candidate pool from approved Phase 0 inputs

Phase or task:

- Phase 1 coarse T/NK extraction

Exact `.py` scripts:

- `workflows/integration/phase1_extract_tnk_candidates.py`
- `workflows/integration/phase1_finalize_from_temp.py`

Core outputs:

- `Integrated_dataset/TNK_candidates.h5ad`
- Phase 1 tables, logs, and PNG QC figures

## Phase 1b: Conservative first cleanup

Objective:

- remove only obvious non-T/NK cells and apply the user-requested low-detection
  gene filter

Phase or task:

- Phase 1b conservative first cleanup

Exact `.py` script:

- `workflows/integration/phase1b_conservative_cleanup.py`

Core outputs:

- in-place update of `Integrated_dataset/TNK_candidates.h5ad`
- Phase 1b tables, logs, and PNG QC figures

## Phase 1c: Merged metadata backup and replacement

Objective:

- export merged `adata.obs`
- back up the prior harmonized metadata file
- replace the harmonized metadata target with validated merged metadata

Phase or task:

- Phase 1c merged metadata backup and replacement

Exact `.py` script:

- `workflows/integration/phase1c_replace_harmonized_metadata.py`

Required join key:

- `project name`
- `sampleid`
- `barcodes`

Core outputs:

- backup metadata CSV
- merged-obs export
- replacement summary tables and log

## Phase 2: Merged cleanup

Objective:

- perform merged-context cleanup on the candidate milestone

Phase or task:

- Phase 2 merged cleanup

Exact `.py` script:

- `workflows/integration/phase2_merged_cleanup.py`

Core outputs:

- `Integrated_dataset/TNK_cleaned.h5ad`
- Phase 2 tables, logs, and PNG QC figures

## Phase 3: scVI integration

Objective:

- integrate the cleaned object with scVI
- build latent representation, neighbors, Leiden, and UMAP
- optionally keep scANVI outputs as reference-only when approved

Phase or task:

- Phase 3 scVI integration and optional scANVI reference annotation

Exact `.py` script:

- `workflows/integration/phase3_scvi_scanvi.py`

Standard behavior:

- exclude mitochondrial, ribosomal, and noncoding-like genes from HVGs when
  requested
- keep large H5AD execution on the mirrored SSD tree when active in the runbook
- keep logs, PNG figures, tables, and models on NFS

Core outputs:

- `Integrated_dataset/integrated.h5ad`
- Phase 3 tables, logs, and PNG QC figures

## Phase 4: TRAB/TRB/TRD scoring

Objective:

- compute package-faithful continuous TRA/TRB/TRD module scores on the
  integrated object

Phase or task:

- Phase 4 TRAB/TRB/TRD continuous scoring

Exact `.py` scripts:

- `workflows/integration/phase4_gdt_module_scoring.py`
- `workflows/integration/plot_phase4_threshold_barplots.py`

Standard behavior:

- score from a temporary normalize-plus-log1p view of count-space `X`
- write continuous score columns back into `integrated.h5ad`
- keep hard labels and downstream threshold summaries separate from the raw
  score columns

Core outputs:

- updated `Integrated_dataset/integrated.h5ad`
- Phase 4 tables, logs, and PNG QC figures

## Extended full-atlas rebuild

Objective:

- reproducibly rebuild the missing historical five-million-cell atlas
- add the eight approved T/NK-filtered extension cohorts and the
  `BALF_BLOOD_COPD` cohort
- produce a structurally integrated, scored atlas without mutating source H5ADs

Phase or task:

- full T/NK atlas sparse rebuild, A100 scVI integration, RAPIDS embedding, and
  Phase 4 scoring

Exact `.py` script and frozen config:

- `workflows/integration/rebuild_full_atlas.py`
- `configs/datasets/full_atlas_rebuild.json`

Execution stages:

- `preflight`
- `prepare`
- `fit`
- `embed`
- `assemble`
- `score`
- `report`

Standard behavior:

- reconstruct the exact 5,128,904-cell historical baseline from the preserved
  cleaned core plus six frozen historical expansion cohorts
- add 758,135 extension T/NK cells and 46,273 COPD BALF/BLOOD cells
- preserve sparse raw counts in `X`, align to the 27,413-gene historical
  universe, and prefix cell IDs with their physical input cohort
- select 4,000 source-balanced HVGs after excluding TCR V/J/D, mitochondrial,
  ribosomal, immunoglobulin, and common noncoding genes; retain the fixed T/NK
  context panel
- require A100 scVI and RAPIDS with no CPU fallback
- append continuous Phase 4 scores without replacing raw-count `X`
- preserve all source H5ADs read-only
- do not consume the validated 14-source repaired-TCR sidecar before its
  separate backup/replacement gate
- preserve source metadata and block GSE169246 tissue/specimen interpretation
  until its additive `_b` blood-library correction is applied

Key outputs:

- `Integrated_dataset/integrated_full_atlas.h5ad`
- `/ssd/tnk_phase3/Integrated_dataset/full_atlas/integrated_full_atlas.h5ad`
- `/ssd/tnk_phase3/Integrated_dataset/full_atlas/scvi_model/`
- `Integrated_dataset/tables/full_atlas_rebuild/`
- `Integrated_dataset/figures/full_atlas_rebuild/`
- `Integrated_dataset/logs/full_atlas_rebuild/`
- `full_atlas_rebuild/index.html`
- `full_atlas_rebuild/full_atlas_rebuild_report.pdf`

Validated completion state:

- 5,933,312 cells x 27,413 genes
- 16 physical input cohorts and 40 source accessions
- 9,127,088,723 sparse nonzero counts
- finite 30-dimensional scVI latent, finite UMAP, 33 Leiden clusters, and eight
  finite Phase 4 score columns for every cell

## Post-Phase-4 milestone curation

Objective:

- apply approved milestone-level removals or carve-outs after Phase 4 review

Phase or task:

- Post-Phase-4 milestone curation for approved GSE removal

Exact `.py` script:

- `workflows/maintenance/remove_no_tcr_gene_gses_from_milestones.py`

Core outputs:

- `high_speed_temp/Integrated_dataset/No_TCR_Gene_GSEs.h5ad`
- rewritten milestone H5ADs without the approved target GSEs
- `Integrated_dataset/tables/no_tcr_gene_gse_removal_counts.csv`

## Post-Phase-4 extended five-million-atlas reporting refinement

Objective:

- extend approved downstream review packages with focused T/NK/γδ figures,
  tissue-distribution statistics, and refreshed HTML/PDF reports

Phase or task:

- Post-Phase-4 downstream reporting refinement on integrated milestones

Exact `.py` scripts:

- `workflows/reporting/plot_extended_atlas_tcr_pairing_umap.py`
- `workflows/reporting/plot_extended_atlas_sorted_gdt_umap.py`
- `workflows/reporting/plot_extended_atlas_tnk_marker_umaps.py`
- `workflows/reporting/plot_extended_atlas_cell_type_umap.py`
- `workflows/reporting/build_extended_atlas_gdt_report_assets.py`
- `workflows/reporting/render_extended_atlas_report.py`

The historical `plus6` output basenames below are retained for provenance and
stable links; active entry points use the unambiguous `extended_atlas` name.

Core outputs:

- `Integrated_dataset/figures/plus6/plus6_umap_paired_tcr_doublets.png`
- `Integrated_dataset/figures/plus6/plus6_umap_sorted_gdt_highlight.png`
- `Integrated_dataset/figures/plus6/plus6_umap_paired_tcr_sorted_gdt.png`
- `Integrated_dataset/figures/plus6/plus6_tnk_marker_umap_panel.png`
- `Integrated_dataset/figures/plus6/plus6_umap_by_corrected_simple_annotation_100x120mm.png`
- `Integrated_dataset/tables/plus6/plus6_gdt_candidate_statistics.csv`
- `Integrated_dataset/tables/plus6/plus6_gdt_candidate_overlap_gt0p4.csv`
- `Integrated_dataset/tables/plus6/plus6_gdt_paired_gdtcr_by_tissue.csv`
- `Integrated_dataset/tables/plus6/plus6_gdt_three_criteria_by_tissue.csv`
- `Integrated_dataset/plus6_profile_report.md`
- `Integrated_dataset/plus6_profile_report.html`
- `Integrated_dataset/plus6_profile_report.pdf`
- `Integrated_dataset/logs/no_tcr_gene_gse_removal.md`

## Post-Phase-4 annotation evidence review

Objective:

- audit cluster annotations with all-cell marker signatures, Phase 4 scores,
  TCR-chain metadata, reference labels as weak priors, and dataset/tissue
  composition
- produce review-only proposed labels and confidence flags without modifying
  milestone H5AD annotation columns

Phase or task:

- Post-Phase-4 cluster annotation evidence review

Exact `.py` script:

- `workflows/analysis/annotation_evidence_review.py`

Core outputs:

- `Integrated_dataset/tables/annotation_evidence_review/current_integrated_cluster_annotation_proposals.csv`
- `Integrated_dataset/tables/annotation_evidence_review/current_integrated_cluster_signature_scores.csv`
- `Integrated_dataset/tables/annotation_evidence_review/current_integrated_cluster_marker_gene_stats.csv`
- `Integrated_dataset/logs/annotation_evidence_review/current_integrated_annotation_evidence_summary.md`
- `Integrated_dataset/figures/annotation_evidence_review/current_integrated_cluster_signature_heatmap.png`
- `Integrated_dataset/tables/annotation_evidence_review/plus6_cluster_annotation_proposals.csv`
- `Integrated_dataset/tables/annotation_evidence_review/plus6_cluster_signature_scores.csv`
- `Integrated_dataset/tables/annotation_evidence_review/plus6_cluster_marker_gene_stats.csv`
- `Integrated_dataset/logs/annotation_evidence_review/plus6_annotation_evidence_summary.md`
- `Integrated_dataset/figures/annotation_evidence_review/plus6_cluster_signature_heatmap.png`

## Supplementary 10x 5' intake lane

Approved supplementary GSEs:

- `GSE179994`
- `GSE235863`
- `GSE240865`
- `GSE287301`
- `GSE234069`
- `GSE287541`

Rules:

- process only approved 10x 5' inputs
- never mix in 10x 3'
- for `GSE234069`, use only `downloads/GSE234069/suppl/10x_5/`
- if TCR is external, integrate it according to `TCR_INTEGRATION_SOP.md`
- write supplementary per-GSE H5ADs under `downloads/per_gse_h5ad_with_metadata/`

## External T-cell/TNK intake with sample-aware TCR integration

Objective:

- process user-supplied T-cell or TNK-subset datasets outside the main registry
- standardize metadata to the project schema
- introduce productive alpha-beta and/or gamma-delta TCR fields by
  `sample_id + barcode_core`

Phase or task:

- External dataset intake with sample-aware productive TCR integration

Exact `.py` script:

- `workflows/intake/process_hra005041_tcr_intake.py`

Standard behavior:

- preserve the input H5AD cell universe
- normalize canonical metadata fields such as `project name`, `sampleid`,
  `sample_id`, `cell_id`, `barcodes`, `barcode`, and `barcode_core`
- introduce only productive alpha-beta TCR rows
- introduce only productive-like gamma-delta TCR rows with valid CDR3 amino-acid
  and nucleotide sequence
- join TCR only by `sample_id + barcode_core`, never barcode alone

Core outputs:

- `downloads/per_gse_h5ad_with_metadata/HRA005041_T_cells_subset.h5ad`
- `Integrated_dataset/tables/HRA005041_tcr_intake/HRA005041_tcr_join_summary.csv`
- `Integrated_dataset/tables/HRA005041_tcr_intake/HRA005041_tcr_sample_summary.csv`
- `Integrated_dataset/tables/HRA005041_tcr_intake/HRA005041_tcr_unmatched_summary.csv`
- `Integrated_dataset/logs/HRA005041_tcr_intake.md`
- use `Integrated_dataset/TNK_candidates_supp.h5ad` only as the supplementary
  milestone before approved merge into the main candidate object

## Standalone external Phase 4 scoring review

Objective:

- compute project-standard TRAB/TRB/TRD module scores on one standalone intake
  H5AD
- write the score columns back into that intake H5AD
- generate standalone score QC and paired-TCR scatter plots

Phase or task:

- Standalone Phase 4 scoring review for external intake H5ADs

Exact `.py` script:

- `workflows/intake/phase4_score_single_h5ad.py`

Core outputs:

- updated `downloads/per_gse_h5ad_with_metadata/HRA005041_T_cells_subset.h5ad`
- `Integrated_dataset/tables/HRA005041_phase4/`
- `Integrated_dataset/figures/HRA005041_phase4/`
- `Integrated_dataset/logs/HRA005041_phase4.log`

## Sorted gdT Seurat intake lane

Objective:

- convert user-supplied Seurat RDS objects of sorted gdT cells into Scanpy
  H5ADs
- preserve raw RNA counts
- harmonize metadata and embedded productive-like gamma-delta TCR fields to the
  project schema
- add `Sorted_gdT = True`
- recompute Scanpy UMAP
- compute project-standard standalone Phase 4 `TRD/TRAB` scores
- export UMAPs plus raw `TRAB`-vs-`TRD` scatter colored by paired `TRG/TRD`

Phase or task:

- Sorted gdT Seurat intake with standalone Phase 4 scoring

Exact `.py` scripts:

- `workflows/intake/process_sorted_gdt_rds_intake.py`
- `workflows/intake/phase4_score_single_h5ad.py`

Standard behavior:

- convert Seurat `RNA` assay counts to H5AD without densifying
- preserve raw counts in the intake H5AD and normalize only on temporary copies
- harmonize canonical metadata fields such as `project name`, `sampleid`,
  `sample_id`, `library_id`, `cell_id`, `barcodes`, `barcode`, and
  `barcode_core`
- support embedded `TRG/TRD` metadata as first-class TCR fields
- keep only productive-like `TRG` and `TRD` chains with valid chain label plus
  both amino-acid and nucleotide CDR3
- set `Sorted_gdT = True`, `input_population = purified_gdt`, and
  `tcr_chain_mode = gd_only`

Core outputs:

- `newdata/Sorted_gdT/GDT_2020AUG_woCOV_sorted_gdt.h5ad`
- `newdata/Sorted_gdT/GDTlung2023july_7p_sorted_gdt.h5ad`
- `newdata/Sorted_gdT/MalteGDT_sorted_gdt.h5ad`
- `Integrated_dataset/tables/Sorted_gdT/`
- `Integrated_dataset/logs/Sorted_gdT/`
- `Integrated_dataset/figures/GDT_2020AUG_woCOV_phase4/`
- `Integrated_dataset/figures/GDTlung2023july_7p_phase4/`
- `Integrated_dataset/figures/MalteGDT_phase4/`

## Supplementary 10x 5' intake lane

Exact `.py` scripts:

- `workflows/intake/supplementary_10x5_phase01.py`
- `workflows/intake/validate_supplementary_10x5_schema.py`

Key outputs:

- `downloads/per_gse_h5ad_with_metadata/`
- `Integrated_dataset/TNK_candidates_supp.h5ad`
- supplementary tables, logs, and PNG QC figures

## Post-Phase-4 gdT prediction package-style evaluation

Objective:

- evaluate package-style gdT score classifiers on the historical extended
  five-million-cell integrated object
- build project-specific primary gold and sensitivity-only silver labels from
  TCR metadata without modifying any H5AD
- render a static HTML/PDF report for downstream review

Phase or task:

- gdT prediction evaluation using shared package classifier logic

Exact `.py` scripts:

- `workflows/gdtai/run_gdt_prediction_package_evaluation.py`
- `workflows/gdtai/plot_gdt_truth_trd_trab_marker_scatter.py`
- `workflows/gdtai/run_tn_fn_trd_gt0p1_deg_reproducibility.py`
- `workflows/gdtai/run_trd_gt0p1_high_no_ab_vs_low_ab_deg_reproducibility.py`
- `workflows/gdtai/run_gdt_deg_tcr_classifier_training.py`
- `workflows/gdtai/run_gdt_gse144469_holdout_tcrgene_classifier.py`
- `workflows/gdtai/run_gdt_classifier_failure_mode_audit.py`
- `workflows/gdtai/predict_with_selected_gdt_model.py`
- `workflows/gdtai/build_gdtai_overview_report.py`
- `workflows/gdtai/run_gdtai_nkguard_classifier.py`
- `workflows/gdtai/run_gdtai_multiguard_classifier.py`
- `workflows/gdtai/run_gdtai_transcriptome_gate_cascade.py`
- `workflows/gdtai/run_gdtai_annotation_specific_cascade.py`
- `workflows/gdtai/run_gdtai_annotation_specific_tp_fn_audit.py`

Core outputs:

- `Integrated_dataset/tables/gdT_prediction/`
- `Integrated_dataset/figures/gdT_prediction/`
- `Integrated_dataset/logs/gdT_prediction/gdT_prediction_summary.md`
- `Integrated_dataset/figures/gdT_prediction/trd_vs_trab_gold_silver_marker_expression.png`
- `Integrated_dataset/tables/gdT_prediction/tn_fn_trd_gt0p1_global_deg.csv`
- `Integrated_dataset/tables/gdT_prediction/tn_fn_trd_gt0p1_deg_reproducibility.csv`
- `Integrated_dataset/figures/gdT_prediction/tn_fn_trd_gt0p1_global_deg_volcano.png`
- `Integrated_dataset/figures/gdT_prediction/tn_fn_trd_gt0p1_source_logfc_reproducibility.png`
- `Integrated_dataset/logs/gdT_prediction/tn_fn_trd_gt0p1_deg_reproducibility_summary.md`
- `Integrated_dataset/tables/gdT_prediction/trd_gt0p1_high_no_ab_vs_expanded_negative_no_tcr_genes_global_deg.csv`
- `Integrated_dataset/tables/gdT_prediction/trd_gt0p1_high_no_ab_vs_expanded_negative_no_tcr_genes_deg_reproducibility.csv`
- `Integrated_dataset/tables/gdT_prediction/trd_gt0p1_high_no_ab_vs_expanded_negative_no_tcr_genes_excluded_tcr_genes.csv`
- `Integrated_dataset/figures/gdT_prediction/trd_gt0p1_high_no_ab_vs_expanded_negative_no_tcr_genes_global_deg_volcano.png`
- `Integrated_dataset/logs/gdT_prediction/trd_gt0p1_high_no_ab_vs_expanded_negative_no_tcr_genes_deg_reproducibility_report.md`
- `gdT_prediction/trd_gt0p1_high_no_ab_vs_expanded_negative_no_tcr_genes_deg_reproducibility_report.html`
- `gdT_prediction/trd_gt0p1_high_no_ab_vs_expanded_negative_no_tcr_genes_deg_reproducibility_report.pdf`
- `Integrated_dataset/tables/gdT_prediction/classifier_training/classifier_metrics_validation.csv`
- `Integrated_dataset/tables/gdT_prediction/classifier_training/classifier_acceptance_vs_best_baseline.csv`
- `Integrated_dataset/tables/gdT_prediction/classifier_training/prevalence_aware_ppv_scenarios.csv`
- `Integrated_dataset/tables/gdT_prediction/classifier_training/selected_model_full_dataset_prediction_overall.csv`
- `Integrated_dataset/tables/gdT_prediction/classifier_training/selected_model_full_dataset_prediction_by_source.csv`
- `Integrated_dataset/logs/gdT_prediction/classifier_training/gdt_deg_tcr_classifier_report.md`
- `Integrated_dataset/models/gdT_prediction_classifier/selected_gdt_classifier.pkl`
- `gdT_prediction/gdt_deg_tcr_classifier_report.html`
- `gdT_prediction/gdt_deg_tcr_classifier_report.pdf`
- `Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene/validation_metrics.csv`
- `Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene/acceptance_vs_individual_gene_baseline.csv`
- `Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene/selected_model_full_dataset_prediction_overall.csv`
- `Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene/selected_model_full_dataset_false_positive_overall.csv`
- `Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene/selected_model_full_dataset_false_positive_by_source.csv`
- `Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene/selected_model_full_dataset_prediction_by_annotation.csv`
- `Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene/trd_vs_trab_prediction_scatter_sample.csv.gz`
- `Integrated_dataset/figures/gdT_prediction/gse144469_holdout_tcrgene/trd_vs_trab_prediction_method_agreement.png`
- `Integrated_dataset/figures/gdT_prediction/gse144469_holdout_tcrgene/trd_vs_trab_tcrgene_known_fp_status.png`
- `Integrated_dataset/logs/gdT_prediction/gse144469_holdout_tcrgene/gse144469_holdout_tcrgene_report.md`
- `Integrated_dataset/models/gdT_prediction_classifier/gse144469_holdout_tcrgene/selected_model.pkl`
- `gdT_prediction/gse144469_holdout_tcrgene_report.html`
- `gdT_prediction/gse144469_holdout_tcrgene_report.pdf`
- `Integrated_dataset/tables/gdT_prediction/gse144469_holdout_tcrgene_failure_modes/`
- `Integrated_dataset/figures/gdT_prediction/gse144469_holdout_tcrgene_failure_modes/`
- `Integrated_dataset/logs/gdT_prediction/gse144469_holdout_tcrgene_failure_modes/gse144469_holdout_tcrgene_failure_mode_audit.md`
- `gdT_prediction/gse144469_holdout_tcrgene_failure_mode_audit.html`
- `gdT_prediction/gse144469_holdout_tcrgene_failure_mode_audit.pdf`
- `gdT_prediction/gdT_classifier_local_handoff_plan.md`
- `Integrated_dataset/tables/gdT_prediction/external_tests/<DATASET_ID>/`
- `Integrated_dataset/figures/gdT_prediction/external_tests/<DATASET_ID>/`
- `Integrated_dataset/logs/gdT_prediction/external_tests/<DATASET_ID>/`
- `gdT_prediction/external_tests/<DATASET_ID>/index.html`
- `Integrated_dataset/logs/gdT_prediction/gdTAI_overview_report.md`
- `Integrated_dataset/figures/gdT_prediction/gdTAI_overview/`
- `gdT_prediction/gdTAI_overview_report.html`
- `Integrated_dataset/tables/gdT_prediction/gdtai_nkguard/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_nkguard/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_nkguard/gdtai_nkguard_report.md`
- `Integrated_dataset/models/gdT_prediction_classifier/gdtai_nkguard/best_candidate_model.pkl`
- `gdT_prediction/gdtai_nkguard_report.html`
- `gdT_prediction/gdtai_nkguard_report.pdf`
- `Integrated_dataset/tables/gdT_prediction/gdtai_multiguard/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_multiguard/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_multiguard/gdtai_multiguard_report.md`
- `Integrated_dataset/models/gdT_prediction_classifier/gdtai_multiguard/best_candidate_model.pkl`
- `gdT_prediction/gdtai_multiguard_report.html`
- `gdT_prediction/gdtai_multiguard_report.pdf`
- `Integrated_dataset/tables/gdT_prediction/gdtai_transcriptome_gate/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_transcriptome_gate/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_transcriptome_gate/gdtai_transcriptome_gate_report.md`
- `Integrated_dataset/models/gdT_prediction_classifier/gdtai_transcriptome_gate/gdtai_transcriptome_gate_protocol.md`
- `gdT_prediction/gdtai_transcriptome_gate_report.html`
- `gdT_prediction/gdtai_transcriptome_gate_report.pdf`
- `Integrated_dataset/tables/gdT_prediction/gdtai_annotation_specific_cascade/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_annotation_specific_cascade/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_annotation_specific_cascade/gdtai_annotation_specific_cascade_report.md`
- `Integrated_dataset/models/gdT_prediction_classifier/gdtai_annotation_specific_cascade/selected_annotation_specific_wrapper.pkl`
- `Integrated_dataset/models/gdT_prediction_classifier/gdtai_annotation_specific_cascade/gdtai_annotation_specific_cascade_protocol.md`
- `gdT_prediction/gdtai_annotation_specific_cascade_report.html`
- `gdT_prediction/gdtai_annotation_specific_cascade_report.pdf`
- `Integrated_dataset/tables/gdT_prediction/gdtai_annotation_specific_tp_fn_audit/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_annotation_specific_tp_fn_audit/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_annotation_specific_tp_fn_audit/gdtai_annotation_specific_tp_fn_audit_report.md`
- `gdT_prediction/gdtai_annotation_specific_tp_fn_audit.html`
- `gdT_prediction/gdtai_annotation_specific_tp_fn_audit.pdf`
- `gdT_prediction/index.html`
- `gdT_prediction/gdT_prediction_report.pdf`

Phase or task:

- Post-Phase-4 gdTAI v3 TRDC/NK-guard classifier training, cross-study benchmark testing, and failure audit

Exact `.py` scripts:

- `workflows/gdtai/run_gdtai_v3_trdc_nk_guard_classifier.py`
- `workflows/gdtai/run_gdtai_v3_trdc_nk_failure_audit.py`

Core outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v3_trdc_nk_guard/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v3_trdc_nk_guard/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v3_trdc_nk_guard/`
- `gdT_prediction/gdtai_v3_trdc_nk_guard_report.html`
- `gdT_prediction/gdtai_v3_trdc_nk_guard_report.pdf`
- `gdT_prediction/gdtai_v3_trdc_nk_failure_audit.html`
- `gdT_prediction/gdtai_v3_trdc_nk_failure_audit.pdf`
- if promoted: `Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/gdTAI_v3_model.pkl` plus README, methodology, feature, metrics, and example files

Standard behavior:

- include `GSE144469` primary gold cells in atlas train/tune splits
- keep the non-GSE144469 validation cohorts from the prior multi-cohort validation lane: `GSE254249` paired-TCRAB/no-gdTCR negatives and `GDT_2020AUG_woCOV` cord-blood gdT positives
- keep the cross-study H5AD out of coefficient fitting; historical candidate
  and promotion workflows later inspected it for selection, so it is not an
  independent external test
- after training, apply the selected v3 candidate and reference strategies to the full atlas input H5AD and write whole-atlas source/tissue/annotation summaries plus selected predicted cells
- run at least five iteration rounds unless the recall/estimated-FP target is already reached; candidate families are not restricted to XGBoost
- select the final review candidate from full-atlas target metrics using recall `> 0.8` and estimated FP fraction `< 5%` of predictions, estimating hidden abT FP from paired-TCRAB rates in TCR-sequenced sources
- require cross-study H5AD inference from `layers["counts"]` and fail rather than using normalized `X`
- exclude NK+TCRAB-overlap cells from negative train/tune labels and explicit hard-negative construction
- restrict explicit NK hard negatives to expression `TRDC+TRDV-`; non-NK TCRAB hard negatives remain a separate guard set
- report ambiguous/mixed/partial external TCR labels as stress groups outside primary binary metrics

Phase or task:

- Post-Phase-4 gdTAI NK-optimized gdT atlas metadata harmonization and report

Exact `.py` script:

- `workflows/gdt_atlas/run_gdtai_nk_optimized_gdt_atlas.py`

Core outputs:

- `configs/gdt_atlas/metadata_rules.json`
- `Integrated_dataset/tables/gdT_atlas/`
- `Integrated_dataset/figures/gdT_atlas/`
- `Integrated_dataset/logs/gdT_atlas/gdtai_nk_optimized_gdt_atlas_summary.md`
- `gdT_atlas/index.html`
- `gdT_atlas/gdT_atlas_report.pdf`

Phase or task:

- Post-Phase-4 curated gdT phenotype atlas

Exact `.py` script:

- `workflows/gdt_atlas/run_gdt_atlas_curated_phenotype_analysis.py`

Core outputs:

- `configs/gdt_atlas/phenotype_marker_panels.json`
- `Integrated_dataset/tables/gdT_atlas_curated_phenotypes/`
- `Integrated_dataset/figures/gdT_atlas_curated_phenotypes/`
- `Integrated_dataset/logs/gdT_atlas_curated_phenotypes/gdt_atlas_curated_phenotypes_summary.md`
- `gdT_atlas/curated_phenotypes/index.html`
- `gdT_atlas/curated_phenotypes/gdT_atlas_curated_phenotypes_report.pdf`

## QC and risk-gate note

Default rule:

- routine, reversible phases proceed after documented QC
- explicit approval is required only for high-risk operations as defined in
  `TNK_PIPELINE_RUNBOOK.md`

Do not define active exceptions here.
Active exceptions belong only in `TNK_PIPELINE_RUNBOOK.md`.

## Post-Phase-4 rigorous gdT atlas reanalysis

Objective:

- rebuild the gdT biology atlas from promoted gdTAI predictions with no
  silver-only FN add-back
- remove likely false positives, add back only gold-level FN cells, re-run
  gdT-only scVI integration, UMAP, Leiden clustering, phenotype annotation,
  sample-level statistics, and pseudobulk-style DE

Phase or task:

- Post-Phase-4 rigorous gdT atlas reanalysis and report

Exact `.py` script:

- `workflows/gdt_atlas/run_gdt_atlas_rigorous_reanalysis.py`

Core outputs:

- `Integrated_dataset/gdT_atlas/gdt_curated_integrated.h5ad`
- `Integrated_dataset/tables/gdT_atlas_rigorous/`
- `Integrated_dataset/figures/gdT_atlas_rigorous/`
- `Integrated_dataset/logs/gdT_atlas_rigorous/gdt_atlas_rigorous_summary.md`
- `gdT_atlas/rigorous/index.html`
- `gdT_atlas/rigorous/gdT_atlas_rigorous_report.pdf`

## gdTAI v3 decision-report refinement

Objective:

- replace the dense gdTAI v3 candidate report with a reader-facing decision report
  that emphasizes conclusions, compact figures, and promotion rationale

Phase or task:

- Post-Phase-4 gdTAI v3 report refinement

Exact `.py` script:

- `workflows/gdtai/render_gdtai_v3_clear_report.py`

Core outputs:

- `Integrated_dataset/logs/gdT_prediction/gdtai_v3_trdc_nk_guard/gdtai_v3_trdc_nk_guard_report.md`
- `gdT_prediction/gdtai_v3_trdc_nk_guard_report.html`
- `gdT_prediction/gdtai_v3_trdc_nk_guard_report.pdf`
- compact decision PNG figures under `gdT_prediction/assets/gdtai_v3_trdc_nk_guard/`

## gdTAI v3.0 promotion and predicted gdT atlas handoff

Objective:

- promote the user-approved gdTAI v3.0 model package
- remove previous gdT atlas H5AD/report output trees
- create a handoff H5AD containing promoted gdTAI v3 predicted gdT cells,
  primary-gold gdT FN add-back cells, and no known/likely paired-TCRAB FP cells

Phase or task:

- Post-Phase-4 gdTAI v3.0 promotion and predicted gdT-cell atlas handoff

Exact `.py` script:

- `workflows/gdtai/promote_gdtai_v3_and_build_predicted_gdt_atlas.py`

Core outputs:

- `Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/gdTAI_v3_model.pkl`
- `Integrated_dataset/gdT_atlas/predicted_gdt_cell_atlas.h5ad`
- `Integrated_dataset/gdT_atlas/predicted_gdt_cell_atlas.md`
- `Integrated_dataset/tables/gdT_atlas/predicted_gdt_cell_atlas_summary.csv`
- `Integrated_dataset/tables/gdT_atlas/predicted_gdt_cell_atlas_metadata.csv.gz`
- `Integrated_dataset/tables/gdT_atlas/predicted_gdt_cell_atlas_removed_fp_cells.csv.gz`
- `Integrated_dataset/tables/gdT_atlas/predicted_gdt_cell_atlas_gold_fn_added_cells.csv.gz`

## Post-Phase-4 predicted gdT atlas tumor-type audit

Objective:

- resolve broad tumor-like atlas tissue labels to original tumor types by source GSE and row-level tissue context without modifying the atlas H5AD

Phase or task:

- Predicted gdT atlas tumor-type source audit

Exact `.py` script:

- `workflows/gdt_atlas/run_gdt_atlas_tumor_type_audit.py`

Core outputs:

- `Integrated_dataset/tables/gdT_atlas/tumor_type_source_mapping.csv`
- `Integrated_dataset/tables/gdT_atlas/tumor_type_by_source_tissue.csv`
- `Integrated_dataset/tables/gdT_atlas/tumor_type_cell_annotation_review.csv.gz`
- `Integrated_dataset/logs/gdT_atlas/tumor_type_source_mapping.md`

Standard behavior:

- read `Integrated_dataset/gdT_atlas/predicted_gdt_cell_atlas.h5ad` in backed mode
- use `configs/gdt_atlas/metadata_rules.json` source-level tumor-type rules gated by row tissue values
- preserve raw tissue columns and keep tumor-type fields as review-only exports until write-back is explicitly approved

## Post-Phase-4 gdTAI model comparison

Objective:

- compare the latest local gdTAI model against gdTAI v2.0 high-purity at
  whole-atlas and per-dataset levels without mutating H5AD files

Phase or task:

- gdTAI latest-vs-v2 high-purity prediction overlap audit

Exact `.py` script:

- `workflows/gdtai/compare_gdtai_latest_vs_v2_high_purity.py`

Core outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v3_vs_v2_high_purity/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v3_vs_v2_high_purity/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v3_vs_v2_high_purity/`
- `gdT_prediction/gdtai_v3_vs_v2_high_purity_report.html`
- `gdT_prediction/gdtai_v3_vs_v2_high_purity_report.pdf`

## Post-Phase-4 gdTAI v3 Round 12 versus Round 14 revalidation

Objective:

- resolve the canonical-model conflict by comparing the actual Round 12 and
  Round 14 artifacts on identical cohorts without mutating H5AD files
- select one balanced default by prespecified recall and false-positive
  guardrails, while retaining the other artifact as a documented fallback

Phase or task:

- checksum-pinned gdTAI v3 Round 12 versus Round 14 comparison and promotion

Exact `.py` script:

- `workflows/gdtai/compare_gdtai_v3_round12_vs_round14.py`

Promotion command:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/compare_gdtai_v3_round12_vs_round14.py \
  --promote-selected
```

Core outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v3_round12_vs_round14/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v3_round12_vs_round14/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v3_round12_vs_round14/`
- `Integrated_dataset/models/gdT_prediction_classifier/gdtai_v3_round12_vs_round14/`
- `gdT_prediction/gdtai_v3_round12_vs_round14/index.html`
- `gdT_prediction/gdtai_v3_round12_vs_round14/gdtai_v3_round12_vs_round14_report.pdf`
- selected package synchronized to
  `Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/`

Standard behavior:

- require exact Round 12 and Round 14 SHA256 identities
- rerun both artifacts directly on the reused BALF_BLOOD_COPD cross-study benchmark
- verify direct predictions against the validated external cache and a
  deterministic full-atlas cache sample
- compare raw full-atlas calls, gold and silver recovery, paired-TCRAB and NK
  burdens, held-out cord-blood positives, and held-out GSE254249 negatives
- run descriptive exact paired tests and report every dataset separately
- historical selection passed guardrails and then maximized mean F1 across
  full-atlas gold, atlas-held-out, and the reused cross-study benchmark; this
  use exhausted benchmark independence and must not be presented as external
  validation
- preserve Round 12 as the high-purity fallback when Round 14 is promoted

## Post-Phase-4 gdTAI v4 two-stage retraining

Objective:

- train and evaluate an experimental two-stage gdTAI classifier with a soft
  T-lineage gate and dropout-tolerant gdT classifier

Phase or task:

- gdTAI v4 T-cell-gate retraining and full-atlas comparison

Exact `.py` script:

- `workflows/gdtai/run_gdtai_v4_tcell_gate_classifier.py`

Core outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_tcell_gate/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v4_tcell_gate/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_tcell_gate/`
- `Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v4.0/`
- `gdT_prediction/gdtai_v4_tcell_gate_report.html`
- `gdT_prediction/gdtai_v4_tcell_gate_report.pdf`

## Project maintenance: dataset-centered pre-integration storage

Objective:

- make `data/datasets/<dataset_id>/` the physical home for raw, interim, and
  processed per-dataset inputs
- preserve all historical `downloads/`, `analysis_26GSE_V4/`, and `newdata/`
  paths as compatibility symlinks to the same inodes
- support validated forward migration and reverse rollback without rewriting
  any H5AD

Phase or task:

- Dataset-centered physical storage migration with legacy-path compatibility

Exact `.py` scripts:

- `workflows/maintenance/migrate_dataset_storage.py`
- `workflows/maintenance/build_data_compatibility_view.py`
- `workflows/maintenance/manage_dataset_registry.py`

Core outputs:

- `data/datasets/<dataset_id>/`
- `data/shared/`
- `data/compat/`
- `data/registry/storage_index.csv`
- `data/registry/migrations/dataset_centered_20260716/`
- top-level compatibility aliases at `downloads`, `analysis_26GSE_V4`, and
  `newdata`

Standard behavior:

- run `plan` and `preflight` before `apply --confirm`
- use same-filesystem renames and verify source device/inode identities
- preserve selected-H5AD sampled hashes, HDF5 dimensions, size, and timestamps
- run `finalize --confirm` only after all physical moves validate
- retain a registry snapshot and reverse journal for `rollback --confirm`

## Local cross-study benchmark intake: BALF_BLOOD_COPD

Objective:

- register the four-library BALF/PBMC study used as the gdTAI cross-study
  benchmark under dataset `BALF_BLOOD_COPD`
- retain the complete raw/interim study workspace at its original path
- store only the selected validation H5AD physically inside the project
- keep the cohort validation-only unless atlas integration is separately
  approved

Phase or task:

- Dataset-centered intake of the BALF/PBMC COPD cross-study benchmark

Exact `.py` scripts:

- `workflows/maintenance/migrate_balf_blood_copd_dataset.py`
- `workflows/maintenance/relocate_balf_blood_copd_raw_workspace.py`

Core outputs:

- `/home/tanlikai/databank/owndata/singlecell/`
- `data/datasets/BALF_BLOOD_COPD/raw/legacy_source`
- `data/datasets/BALF_BLOOD_COPD/processed/current.h5ad`
- `data/registry/migrations/balf_blood_copd_20260716/`
- `data/registry/migrations/balf_blood_copd_h5ad_only_20260716/`
- dataset, library, file, and storage-index registry entries

Standard behavior:

- keep raw, FASTQ, BAM, matrix, TCR, demultiplexing, and other study files in
  the original workspace
- retain the selected validation H5AD as the only physical cohort data file in
  `data/datasets/BALF_BLOOD_COPD/`
- keep the historical H5AD path as a compatibility link to the canonical H5AD
- fingerprint the selected H5AD, four filtered matrices, four raw matrices,
  harmonized TCR table, and two demultiplexing tables
- register `LIB1`-`LIB4` as inactive validation-only libraries

## Post-Phase-4 gdTAI V4.2 NK-reference sidecar integration

Objective:

- build a read-only modeling sidecar that expands low-weight NK development
  controls through scVI and repeated clustering without exposing locked
  evaluation cohorts or changing canonical atlas/source H5AD files

Phase or task:

- gdTAI V4.2 sidecar-integration implementation QC and gated execution

Exact `.py` scripts:

- `workflows/gdtai/run_gdtai_v4_2_implementation_qc.py`
- `workflows/gdtai/run_gdtai_v4_2_nk_reference_integration.py`
- `workflows/gdtai/gdtai_v4_2_integration_core.py`
- `workflows/gdtai/build_gdtai_v4_2_nk_reference_qc_report.py`
- `workflows/gdtai/build_gdtai_v4_2_sidecar_diagnostics.py`

Key outputs:

- `gdT_prediction/gdtai_v4_2_implementation_qc/index.html`
- `gdT_prediction/gdtai_v4_2_implementation_qc/gdtai_v4_2_implementation_qc_report.pdf`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_implementation_qc/`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_implementation_qc/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_implementation_qc/`
- after separate approval, SSD sidecar under
  `/ssd/tnk_phase3/Integrated_dataset/gdtai_v4_2_nk_reference/`
- after separate approval, canonical execution QC under
  `Integrated_dataset/{tables,figures,logs}/gdT_prediction/gdtai_v4_2_nk_reference/`
- completed execution report under
  `gdT_prediction/gdtai_v4_2_nk_reference/`
- deterministic latent-space visualization and confounding diagnostics under
  `gdT_prediction/gdtai_v4_2_nk_reference_diagnostics/` and
  `Integrated_dataset/{tables,figures,logs}/gdT_prediction/gdtai_v4_2_nk_reference_diagnostics/`

Standard behavior:

- require exact preflight and execution-approval hashes and abort before SSD
  access if the current gate is absent
- exclude all locked cohorts before metadata loading, HVG selection, matrix
  staging, scVI, clustering, and pseudo-label construction
- use 4,000 source-balanced common HVGs after TCR, mitochondrial, ribosomal,
  and immunoglobulin exclusion; never force classifier genes into integration
- use direct A100 CUDA for scVI and RAPIDS with no CPU fallback
- define the T/NK boundary only through mixing of independent primary-NK and
  productive-TCR anchors, never through marker-expression thresholds
- publish integration/clustering QC before classifier fitting; routine,
  reversible development may then proceed under the runbook's risk-based rule
- fail closed without classifier fitting when the consensus selects no
  pseudo-NK cells
- use the saved latent and partitions read-only for a deterministic,
  source-balanced diagnostic sample; report cohort composition, anchor roles,
  near-cluster source composition, and same-source neighbor retention
- treat the current sidecar as development-only and do not infer new-cohort NK
  identity when anchor class is confounded with cohort lane

## Corrected full-atlas overview report

Objective:

- export a read-only overview of the rebuilt full T/NK atlas using the
  validated metadata- and TCR-corrected candidate
- separate source-author annotations, unsupervised clusters, productive-TCR
  coverage, and transcriptomic marker evidence

Phase or task:

- corrected full-atlas composition, TCR, signature-expression, and QC report

Exact `.py` script:

- `workflows/reporting/build_full_atlas_overview.py`

Key outputs:

- `full_atlas_overview/index.html`
- `full_atlas_overview/full_atlas_overview_report.pdf`
- `Integrated_dataset/figures/full_atlas_overview/`
- `Integrated_dataset/tables/full_atlas_overview/`
- `Integrated_dataset/logs/full_atlas_overview/`

Standard behavior:

- use exact all-cell metadata and repaired productive-TCR counts
- use `source_accession_harmonized_v2` for user-facing dataset names while
  retaining `source_gse_id` only in the exported provenance mapping
- use deterministic source-balanced samples only for UMAP rendering, marker
  expression, signature summaries, and QC summaries
- compute temporary log1p library-size-normalized marker expression directly
  from raw sparse counts without writing to the H5AD
- report source-label display groups as author-label summaries, not scANVI or
  consensus annotations
- verify the input H5AD size and nanosecond mtime are unchanged

## Post-repair NK-boundary pseudobulk DEG diagnostic

Objective:

- compare frozen NK-boundary clusters enriched for corrected productive TRA or
  TRB with boundary clusters depleted for that evidence
- identify reproducible transcriptomic differences without allowing receptor
  or proximal TCR-complex genes to drive the result

Phase or task:

- sample-paired raw-count pseudobulk DEG and cross-dataset sensitivity analysis

Exact scripts:

- `workflows/gdtai/analyze_boundary_pseudobulk_deg.py`
- `workflows/gdtai/run_boundary_pseudobulk_deg.R`

Key inputs:

- `/ssd/tnk_phase3/Integrated_dataset/full_atlas/tcr_corrected/integrated_full_atlas.h5ad`
- `/ssd/tnk_phase3/Integrated_dataset/gdtai_v4_2_tnk_reintegration/tnk_refined_hvg_counts.h5ad`
- `/ssd/tnk_phase3/Integrated_dataset/gdtai_v4_2_tnk_reintegration/nk_boundary_partitions.npz`

Key outputs:

- `gdT_prediction/gdtai_v4_2_boundary_pseudobulk_deg/index.html`
- `gdT_prediction/gdtai_v4_2_boundary_pseudobulk_deg/boundary_pseudobulk_deg_report.pdf`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_boundary_pseudobulk_deg/`
- `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_boundary_pseudobulk_deg/`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_boundary_pseudobulk_deg/`

Standard behavior:

- use the frozen resolution-0.4 boundary partition and corrected
  `has_any_ab_tcr` metadata without changing any H5AD
- classify clusters by two-sided Fisher enrichment against the remaining
  boundary cells with Benjamini-Hochberg correction; do not call any cluster
  literally TCR-free
- aggregate raw counts within biological sample and comparison group, retaining
  pairs with at least 20 cells in each group
- exclude TRA/TRB/TRD/TRG receptor-locus genes plus `CD3D`, `CD3E`, `CD3G`,
  `CD247`, and `TRAT1` before normalization and testing
- use TMM normalization and limma empirical Bayes on within-sample differences;
  use dataset-mean effects from datasets with at least five pairs as a
  cross-dataset sensitivity analysis
- treat results as descriptive markers of pre-existing expression-defined
  clusters, not causal effects of productive TCR status and not new NK truth

### Conservative NK-like negative-feature triage

Phase or task:

- select candidate negative model features from robust boundary DEGs while
  excluding genes shared with cytotoxic, activated, innate-like, or gamma-delta
  T cells

Exact `.py` script:

- `workflows/gdtai/select_boundary_nk_negative_features.py`

Key outputs:

- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_boundary_pseudobulk_deg/candidate_nk_negative_features.csv`
- `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_boundary_pseudobulk_deg/excluded_shared_nk_t_features.csv`
- `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_boundary_pseudobulk_deg/nk_negative_feature_selection.md`

Standard behavior:

- require robust cross-dataset depletion, paired log2FC at most -0.5, negative
  dataset-macro effect, and source sign consistency of at least 0.88
- exclude `GZMK`, `KLRC1`, every KIR gene, generic cytotoxicity genes, and
  genes commonly expressed by activated or gamma-delta T cells
- treat Tier A (`NCR1`, `SIGLEC7`, `SH2D1B`) as direct innate/NK
  receptor-adaptor candidates and Tier B (`LAT2`, `SYK`, `PLCG2`, `PILRB`,
  `CD300C`) as broader non-T immune-receptor pathway candidates
- use candidates only as regularized soft negative features evaluated inside
  grouped cross-dataset folds; never use one gene as a hard NK veto or label
