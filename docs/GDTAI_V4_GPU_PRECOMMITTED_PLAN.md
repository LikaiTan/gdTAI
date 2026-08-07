# Precommitted gdTAI V4.1-GPU Nested Evaluation Plan

**Status:** proposed supervision protocol; no project-data GPU model has been fitted

**Date:** 2026-08-07

**Primary objective:** improve dataset-macro balanced F1 while preserving the
frozen paired-alpha-beta-T-cell and strict-NK false-positive guardrails

**Promotion status:** gdTAI V3 Round 14 remains the promoted balanced default

## 1. Decision and scope

The CPU-only V4 protocol v1.2 experiment is retired as incomplete and
non-promotable. Its completed HRA005041 and GSE144469 outer folds are retained
as diagnostic evidence. In both folds:

- no V4, compact seven-gene, V2-like, or legacy score candidate had a threshold
  satisfying the frozen operating guardrails;
- all 16 Stage-1 SAGA candidates were nonconverged;
- the selected Stage-1 threshold passed 100% of held-out strict-NK cells; and
- project-data fitting was too slow because SAGA repeatedly exhausted 5,000
  iterations across a 321-fit outer-fold grid.

The BALF_BLOOD_COPD CPU fold was terminated at the user's direction before an
atomic fold result was written. No CPU result will be reused for GPU model
selection. The feature cache and two completed CPU fold reports remain
read-only audit evidence.

V4.1-GPU is a separately named experiment because its algorithms and cascade
are being changed after V4 v1.2 results were observed. This plan authorizes no
project-data fit, model promotion, release artifact, or atlas inference.

## 2. Feasibility established before protocol freeze

The feasibility probes used synthetic data only.

| Item | Observed result | Required interpretation |
| --- | --- | --- |
| GPU | NVIDIA A100 80 GB PCIe; 81,140 MiB free | sufficient |
| Project cache | 1,137,739 x 197 float32 genes, 0.835 GiB; eight derived features, 0.034 GiB | fits in GPU memory with large margin |
| Environment | RAPIDS 24.12, CuPy 13.6, XGBoost 2.1.2, PyTorch 2.9.1+cu128 | pin exact versions |
| MPS | default `/tmp/nvidia-mps` belongs to another user and cannot be used | isolate gdTAI; never alter that daemon |
| Direct CUDA probe | succeeded with a unique `CUDA_MPS_PIPE_DIRECTORY` | feasible without interfering with MPS owner |
| PyTorch weighted ridge logistic | eight LBFGS iterations; bit-identical repeated probabilities | feasible primary linear family |
| GPU XGBoost | bit-identical repeated probabilities; portable UBJ export | feasible nonlinear family |
| cuML logistic | fast, but full-scale synthetic fit emitted a line-search stop and exposed no adequate convergence field | excluded from candidate families |

The production launcher must set all of the following before Python starts:

```bash
CUDA_VISIBLE_DEVICES=0
CUDA_MPS_PIPE_DIRECTORY=/tmp/gdtai-v41-${RUN_ID}-direct
CUBLAS_WORKSPACE_CONFIG=:4096:8
PYTHONHASHSEED=20260807
```

The isolated MPS path must contain no `control` socket. The workflow must never
delete, stop, modify, or connect to `/tmp/nvidia-mps`. If direct CUDA
initialization fails, the run stops; there is no silent CPU fallback.

## 3. Frozen data and biological contract

V4.1-GPU reuses the passed V4 v1.2 Step-1 manifests without changing cells,
labels, groups, or expression inputs:

- cell-label manifest SHA-256:
  `8157cbebfedeb84fc34eba05429aa8a7a834c6f7ceba15bd4790b5dd06bf7e0c`
- grouped-split manifest SHA-256:
  `c84da2ca999676bab0ed180ae3d380e1e7d5b2e08da886ae2f6f912f9e0080a7`
- feature manifest SHA-256:
  `7fc2b0afb2b761288f941895a219fb38cc541a49464c8d2b330b81719e087808`
- CD4/Treg rule SHA-256:
  `518726721c04ce5bfbfd8c0d04bd5289cc4f390d6166c3d644144c921f7ac67b`

Primary expression-independent truth remains:

- `gdT_primary`: productive paired TRG/TRD, with no productive TRA/TRB;
- `abT_primary`: productive paired TRA/TRB, with no productive TRG/TRD;
- `single_abT_weak`: training-only negative with reliability weight 0.5;
- `strict_NK`: high-confidence NK annotation, no productive TCR chain, and no
  doublet flag; dual annotation agreement is required when two annotation
  fields exist, while a single available annotation retains weight 0.5.

Silver, sorted, dual, ambiguous, and unlabeled cells have zero fitting weight.
Silver and sorted cohorts are sensitivity-only. They cannot select features,
algorithms, hyperparameters, calibration, thresholds, or promotion.

The outer folds remain HRA005041, GSE144469, and BALF_BLOOD_COPD. The existing
three grouped inner folds remain checksum-pinned. Donor, sample, library, and
clonotype grouping rules do not change.

## 4. Frozen feature and normalization policy

- Use the same 197 individual genes and eight deterministic derived features.
- Use no Phase-4 TRD, TRAB, or TRD-minus-TRAB module score as a model feature.
- Use no dataset, tissue, disease, donor, library, annotation, or QC metadata as
  a feature.
- Add only the cross-fitted Stage-1 T-lineage probability to Stage 2, for 206
  total numerical inputs and 197 genes.
- Preserve raw-count CP10K plus `log1p` normalization and the registered audited
  HRA005041 transformed-input exception.
- Fit every scaler on the inner-training partition only. Constant and rare
  feature filtering remains fold-local.
- Preserve the fixed CD4-helper and Treg post-model exclusions exactly. They
  cannot be retuned in response to either CPU or GPU results.

## 5. Two-stage GPU architecture

### 5.1 Stage 1: T-lineage versus strict NK

Stage 1 becomes an actual multigene T-lineage gate rather than a no-op soft
feature. It uses the existing frozen 50-gene Stage-1 panel. No individual NK,
cytotoxic, CD4, or Treg gene is a hard veto.

Training uses all eligible cells; no easy-negative subsampling is required
because the complete cache fits in GPU memory. Dataset/class/reliability weights
remain normalized within each inner-training fold.

Candidate families are frozen to:

1. Deterministic PyTorch weighted ridge logistic regression with L2 strengths
   `0.001`, `0.01`, and `0.1`.
2. GPU XGBoost with fixed `learning_rate=0.05`, `n_estimators=200`,
   `subsample=1`, `colsample_bytree=1`, `reg_lambda=10`, and two configurations:
   `max_depth=3, min_child_weight=50` or
   `max_depth=5, min_child_weight=200`.

The Stage-1 threshold is the highest threshold satisfying, in every inner
source, gdT recall at least 0.99 and productive-abT recall at least 0.98. A
Stage-1 candidate is ineligible if strict-NK passage exceeds 50% in any inner
source. Among eligible candidates, minimize dataset-macro strict-NK passage;
ties use higher minimum-source gdT recall, fewer effective features, then the
lexical candidate ID.

This 50% condition prevents another nominal T/NK gate that passes all NK cells.
It is not the final NK safety condition; the final balanced RNA call must still
keep strict-NK FPR at or below 1% in every held-out source.

### 5.2 Stage 2: gdT versus alpha-beta T with fold-local hard NK negatives

Stage 2 receives the 197 genes, eight derived features, and cross-fitted Stage-1
probability. It trains on primary gdT, primary alpha-beta T, weak single-chain
alpha-beta T, and only the strict-NK cells that pass the fold-local Stage-1
gate. These Stage-1-passing NK cells are negative hard examples. They are
identified without access to the inner-validation or outer-held-out fold.

Candidate families are frozen to:

1. Deterministic PyTorch weighted ridge logistic regression with L2 strengths
   `0.0003`, `0.001`, and `0.003`.
2. GPU XGBoost with fixed `learning_rate=0.05`, `n_estimators=300`,
   `subsample=1`, `colsample_bytree=1`, `reg_lambda=10`, and four configurations:
   the Cartesian product of `max_depth={3,5}` and
   `min_child_weight={50,200}`.

No early stopping uses validation labels. PyTorch starts from zero coefficients
and uses deterministic full-batch LBFGS with strong-Wolfe line search. A linear
fit is eligible only when it stops before its fixed 200-iteration limit, has
finite coefficients/probabilities, gradient infinity norm at most `1e-4`, and
relative terminal loss change at most `1e-7`. XGBoost must complete all fixed
rounds with finite probabilities.

The final RNA call requires both Stage-1 and Stage-2 thresholds and then applies
the unchanged CD4/Treg exclusions. Productive-VDJ rescue remains a separately
reported output and cannot improve RNA-only performance or promotion metrics.

## 6. Calibration and operating thresholds

Platt calibration is fitted from cross-fitted inner-training predictions only.
Because it is a one-dimensional optimization, it may run on CPU; main candidate
training must run on GPU. Calibration, threshold selection, metrics, and report
generation are deterministic support tasks, not a silent CPU model fallback.

Operating modes remain unchanged:

- balanced: maximize dataset-macro F1 with macro recall at least 0.80,
  paired-abT FPR at most 0.20%, and strict-NK FPR at most 1.00%;
- high purity: maximize dataset-macro F0.5 with macro recall at least 0.70,
  paired-abT FPR at most 0.10%, and strict-NK FPR at most 0.50%.

The implementation must persist the complete vectorized threshold frontier,
including the best unconstrained threshold and the closest guardrail failures.
No threshold or guardrail may be relaxed after results are observed.

## 7. Fair comparators and reporting

Retrain these comparators with the same outer/inner folds, weights, calibration,
threshold rules, exclusions, and GPU ridge implementation:

- compact seven-gene logistic: `TRDC`, `TRDV1`, `TRDV2`, `TRDV3`, `TRAC`,
  `TRBC1`, and `TRBC2`;
- V2-like individual TCR-gene logistic using the frozen available TCR-gene set;
- the raw legacy TRD-minus-TRAB score, without refitting.

Frozen V2/V3 release profiles are descriptive references only because their
historical development cohorts were reused. They cannot serve as fair promotion
comparators.

Required outputs include per-source and dataset-macro precision, recall,
specificity, F1, F0.5, balanced accuracy, MCC, ROC-AUC, PR-AUC, calibration,
confusion counts, paired-abT FPR, strict-NK FPR, silver/sorted sensitivity,
CD4/Treg losses, threshold frontiers, convergence, runtime, GPU memory,
feature stability, and paired hierarchical-bootstrap differences.

## 8. Checkpointing, reproducibility, and resource controls

- One GPU fit runs at a time. Do not launch concurrent GPU candidates.
- CPU support work is capped at 16 cores; GPU training is the bottleneck target.
- Acquire a project GPU lock and verify GPU memory/utilization before fitting.
- Abort rather than compete if another active compute client uses more than
  8 GB or sustained GPU utilization exceeds 20% for 60 seconds.
- Save each candidate atomically with config, code, cache, split, and environment
  hashes; OOF probabilities and threshold frontiers are checkpointed separately.
- Resume only when every recorded hash matches. Never restart completed
  candidates silently.
- Before project fitting, run one sentinel candidate twice. PyTorch calls must
  be bit-identical; XGBoost calls and probabilities must be bit-identical on the
  same A100/software stack. Any difference blocks the run pending review.
- Export PyTorch coefficients/scaler/calibrator/feature order in NPZ plus JSON,
  and XGBoost models in UBJ plus JSON metadata. Pickle cannot be the only model
  representation.
- No H5AD file is modified.

Expected project-data execution is 30-120 minutes. A three-hour wall-time limit
blocks the run for review instead of allowing an unbounded search.

## 9. Promotion rule

The balanced V4.1-GPU model is eligible for later release review only if all of
the following hold:

1. Dataset-macro F1 exceeds the strongest fairly retrained comparator by at
   least 0.01.
2. The paired hierarchical-bootstrap 95% lower confidence bound for the F1
   difference is greater than zero.
3. Recall is at least 0.70 in every primary outer dataset.
4. Paired-abT and strict-NK guardrails pass in every primary outer dataset and
   every eligible complete-feature extension negative-control cohort.
5. Stage 1 passes its held-out recall and nontrivial-NK-rejection conditions in
   every outer source.
6. Feature direction/rank stability and maximum single-dataset importance share
   pass the existing frozen stability criteria.
7. Determinism, convergence, serialization, semantic, checksum, and no-H5AD-
   mutation checks all pass.

Failure of any condition retains V3 Round 14 as the promoted balanced default.
High-purity success cannot rescue failed balanced promotion.

## 10. Supervision gates

### Gate A: plan review

Commit and email this Markdown, HTML, and PDF plan with the GPU feasibility
record. Stop. No project-data GPU fitting is authorized by Gate A.

### Gate B: implementation and synthetic verification

After explicit approval, implement the GPU core, launcher, atomic checkpoints,
threshold frontier, and automated tests. Run only synthetic tests and read-only
contract/cache validation. Email the implementation QC package and stop.

### Gate C: nested project-data evaluation

After separate explicit approval, run all three outer folds, fair comparators,
extension negative controls, bootstrap, and the print-safe report. Email Likai
and stop for model-selection review.

### Gate D: release decision

Only after explicit review may a selected family be refitted on all development
data. Promotion and whole-atlas inference each require separate approval.
