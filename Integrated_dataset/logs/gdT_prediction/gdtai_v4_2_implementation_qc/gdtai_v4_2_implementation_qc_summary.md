# gdTAI V4.2 sidecar integration implementation QC

## Decision

**PASS_REVIEW_REQUIRED.** The sidecar runner is implemented and passed deterministic unit,
fail-closed, sparse-H5AD, direct-CUDA scVI, and RAPIDS clustering tests. No
project-data integration, project-data scVI fit, pseudo-label construction, or
classifier fit was performed.

- Checks passed: **15/15**.
- Synthetic scVI device: **cuda:0** on **NVIDIA A100 80GB PCIe**.
- Synthetic latent matrix: **1,600 x 8**, all finite.
- Locked cohorts remain excluded by code and the project-data stage aborts while
  its separate approval file is absent.
- CPU fallback is forbidden for scVI and RAPIDS stages.

![Implementation QC](../../Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_implementation_qc/implementation_qc_summary.png)

## Frozen executable design

1. Select 4,000 common non-TCR, non-mitochondrial, non-ribosomal, and
   non-immunoglobulin HVGs from a deterministic source-balanced sample capped at
   20,000 cells per source.
2. Stage only development cohorts in a new sparse H5AD on SSD. Locked cohorts
   are excluded before metadata loading, HVG selection, matrix staging, scVI,
   clustering, and label construction.
3. Fit the precommitted 30-dimensional negative-binomial scVI model on one A100
   GPU, with no CPU fallback.
4. Run nine global and nine boundary RAPIDS Leiden partitions. The boundary is
   defined only by global clusters containing both independent primary-NK and
   productive-TCR anchors.
5. Select low-weight pseudo-NK development candidates only through repeated
   cluster agreement, independent-anchor purity, productive-TCR contamination,
   and cross-source representation. No CD3, TRDC, TRDV, NKG7, GNLY, KLRD1,
   TYROBP, FCER1G, or FCGR3A threshold defines truth.
6. Source balance and the primary-NK effective-weight cap are enforced at the
   later classifier-fitting gate. Pseudo-labels cannot control validation
   guardrails.

The two implementation safety floors frozen before project-data execution are
at least 20 independent anchors and at least 10 primary NK anchors per
qualifying cluster. They prevent tiny-anchor clusters from passing a percentage
rule by chance.

## Checks

| check                                     | status   | detail                                           |
|:------------------------------------------|:---------|:-------------------------------------------------|
| checksum_bound_preflight_approval         | PASS     | 2026-08-08T08:29:03+08:00                        |
| development_locked_role_split             | PASS     | 6 development objects; 3 locked objects          |
| locked_development_permissions            | PASS     | locked cohorts have zero development permissions |
| cpu_fallback_forbidden                    | PASS     | scVI and RAPIDS project stages require CUDA      |
| marker_threshold_truth_forbidden          | PASS     | cluster evidence only                            |
| pseudo_labels_cannot_control_guardrails   | PASS     | primary truth retains guardrail ownership        |
| python_compile                            | PASS     | three modules compiled                           |
| deterministic_unit_and_sparse_h5ad_tests  | PASS     | 11 passed in 2.52s                               |
| runner_validate_stage                     | PASS     | checksum and role validation passed              |
| project_data_stage_fails_without_approval | PASS     | prepare aborted before SSD access                |
| all_source_h5ads_unchanged                | PASS     | 9/9 size/mtime pairs unchanged                   |
| direct_cuda_a100                          | PASS     | NVIDIA A100 80GB PCIe; 79.3 GiB                  |
| synthetic_scvi_gpu_fit                    | PASS     | latent=[1600, 8]; device=cuda:0                  |
| synthetic_rapids_neighbors_leiden         | PASS     | clusters=4,4                                     |
| synthetic_gpu_no_cpu_fallback             | PASS     | direct CUDA only                                 |

## Supervision gate

Project-data integration is still blocked. Activating
`PROJECT_DATA_INTEGRATION_APPROVAL.json` authorizes only development-data sparse
staging, scVI, consensus clustering, and pseudo-label QC. It does not authorize
classifier fitting, threshold selection, model promotion, release fitting, or
whole-atlas inference.

## Exact artifacts

- Execution config: `configs/models/gdtai/v4_2_integration_execution.json`
- Runner: `workflows/gdtai/run_gdtai_v4_2_nk_reference_integration.py`
- Core: `workflows/gdtai/gdtai_v4_2_integration_core.py`
- Approval template: `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_implementation_qc/PROJECT_DATA_INTEGRATION_APPROVAL_TEMPLATE.json`
- Full checks: `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_implementation_qc/implementation_qc_checks.csv`
