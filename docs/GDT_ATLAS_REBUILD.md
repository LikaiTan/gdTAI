# De Novo gdT Atlas Rebuild

## Scope

This workflow builds and analyzes a new gdT-only atlas only after a separately
rebuilt integrated object has passed review. It does not rebuild or mutate the
project integrated milestone.

Configuration:

- `configs/gdt_atlas/atlas_rebuild.json`

Entry points:

- `workflows/gdt_atlas/build_gdt_atlas.py`
- `workflows/gdt_atlas/analyze_gdt_atlas.py`

The atlas is downstream of `compare_frozen_gdtai_profiles.py`. The configured
rebuilt input, checksum-pinned four-profile decision, passed sealed-holdout status, and
approval marker may be absent while upstream work is incomplete. Both atlas
entry points fail closed until all four exist and match exactly.

## Rebuilt Object Contract

The rebuilt H5AD must contain raw, nonnegative integer counts in a CSR
`layers["counts"]` matrix. Observation and variable names must be unique.

The following upstream decisions must be explicit observation columns:

- `gdt_truth_class`: includes reviewed `gdt_gold`, `abt_gold`, and `gdt_silver` values
- `atlas_high_confidence_nk`: strict reviewed NK decision, not an NK-gene-only rule
- `atlas_doublet`: reviewed doublet decision
- `atlas_non_t_contaminant`: reviewed non-T contaminant decision
- `atlas_severe_qc_failure`: reviewed severe-QC decision
- `atlas_targeted_cohort`: marks sorted, enriched, or otherwise targeted cohorts

Paired TCR evidence uses `has_TRA_TRB_paired`, `has_any_ab_tcr`,
`has_TRG_TRD_paired`, `TRG_cdr3`, and `TRD_cdr3`. The remaining required
metadata columns are listed in the config.

The build workflow applies these rules:

1. Verify the canonical digest in the four-profile V2/V3 selection decision,
   reload and checksum every frozen profile identity, and require that the
   linked sealed-holdout status is `holdout_passed` with no schema failure or
   veto. Profile ranking is upstream and jointly considers gdT/T-cell recall,
   alpha-beta T-cell false positives, and strict NK false positives; the atlas
   never reranks or substitutes a runner-up.
2. Infer scores directly from count space with only the profile selected by
   that decision. The threshold comes from the selected frozen artifact and
   mode: fixed for V3 and V2 high-F1, or annotation-specific for V2
   high-purity. The atlas config contains no selected profile or threshold.
3. Remove known alpha-beta-only cells, explicit high-confidence NK cells,
   doublets, non-T contaminants, and severe QC failures. Alpha-beta-only
   requires paired TRA/TRB or reviewed `abt_gold`, with no productive TRG/TRD
   evidence; a lone alpha or beta chain is not a hard exclusion.
4. Add back a false negative only when it is primary gold, has paired TRG/TRD,
   and has valid TRG and TRD CDR3 values.
5. Never add a silver-only false negative.
6. Never exclude a cell because TRDV is absent or because an NK gene alone is
   expressed.

Hard exclusions also apply to gold add-back candidates.

## Approval Marker

The atlas workflow must not create or infer its own approval. Only after profile
selection and the sealed holdout pass, freeze the rebuilt H5AD, atlas config,
comparison workflow, decision, and holdout status. A reviewer then creates the
configured marker:

`Integrated_dataset/approvals/gdt_atlas_rebuilt_integrated.approved.json`

Required structure:

```json
{
  "schema_version": 2,
  "decision": "approved",
  "purpose": "build_gdt_atlas_de_novo",
  "workflow_id": "gdt_atlas_de_novo_v1",
  "approved_by": "REVIEWER_NAME",
  "approved_at_utc": "2026-08-06T12:00:00Z",
  "input_h5ad": "high_speed_temp/Integrated_dataset/rebuilt_integrated.h5ad",
  "input_sha256": "FULL_REBUILT_H5AD_SHA256",
  "input_size_bytes": 0,
  "config_sha256": "FULL_CONFIG_SHA256",
  "gdtai_comparison_workflow_sha256": "FULL_COMPARISON_WORKFLOW_SHA256",
  "gdtai_selection_decision_sha256": "CANONICAL_DECISION_SHA256",
  "gdtai_selection_decision_file_sha256": "FULL_DECISION_FILE_SHA256",
  "gdtai_holdout_status_sha256": "FULL_HOLDOUT_STATUS_FILE_SHA256",
  "gdtai_selected_profile_id": "PROFILE_FROM_SELECTION_DECISION",
  "gdtai_selected_model_id": "MODEL_FROM_SELECTION_DECISION",
  "gdtai_selected_mode": "MODE_FROM_SELECTION_DECISION",
  "gdtai_selected_model_sha256": "SELECTED_MODEL_ARTIFACT_SHA256",
  "gdtai_threshold_contract_sha256": "SELECTED_THRESHOLD_CONTRACT_SHA256"
}
```

The marker pins the full input SHA256, byte size, atlas config, comparison
workflow, canonical decision digest, full decision file, full holdout-status
file, selected model identity, and selected threshold contract. Any later edit
or mismatch invalidates approval.

## Execution

The upstream evaluator must first complete both stages; the atlas does not run
either stage itself:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/compare_frozen_gdtai_profiles.py --stage select

/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/compare_frozen_gdtai_profiles.py --stage holdout
```

After the passed decision and approval marker are frozen, use the canonical
environment and run atlas preflight before either atlas stage:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdt_atlas/build_gdt_atlas.py --preflight-only

/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdt_atlas/build_gdt_atlas.py

/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdt_atlas/analyze_gdt_atlas.py --preflight-only

/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdt_atlas/analyze_gdt_atlas.py
```

The build stage performs gdT-only scVI, one neighbors graph, UMAP across all
configured seeds, and Leiden across every configured seed-resolution pair. TCR
genes are excluded from HVG selection.

The analysis stage writes sample-, donor-, and donor-condition-level cluster
composition. Targeted cohorts remain in descriptive tables but are absent from
all abundance-inference tables. Enabled within-study contrasts are estimated
from donor fractions and pooled with DerSimonian-Laird random effects. The
default disease/control contrast is deliberately disabled until exact
study-level condition mappings are reviewed.

Primary DE is donor-paired cluster-versus-rest pseudobulk DE and excludes all
TRA/TRB/TRG/TRD genes. TRDV and TRGV donor-cluster expression is exported in a
separate report and is never used as an atlas-selection gate.

## Main Outputs

- `Integrated_dataset/gdT_atlas/de_novo/gdt_atlas.h5ad`
- `Integrated_dataset/models/gdT_atlas_de_novo/scvi/`
- `Integrated_dataset/tables/gdT_atlas_de_novo/build/`
- `Integrated_dataset/tables/gdT_atlas_de_novo/analysis/`
- `Integrated_dataset/logs/gdT_atlas_de_novo/build_manifest.json`
- `Integrated_dataset/logs/gdT_atlas_de_novo/analysis_manifest.json`

Existing outputs are never overwritten. A completed build manifest pins the
atlas SHA256, and the analysis stage verifies that digest before reading it.
