# TNK Phase 0-4 Canonical Workflow

## Purpose

This file is the single human-readable canonical workflow for:

- Phase 0: dataset audit and eligibility triage
- Phase 1: coarse T/NK extraction
- Phase 1b: conservative first cleanup
- Phase 1c: merged metadata backup and replacement
- Phase 2: merged cleanup
- Phase 3: scVI integration and scANVI T/NK annotation
- Phase 4: TRAB/TRB/TRD scoring and downstream threshold summaries

The workflow is QC-gated. No phase transition is allowed without user review and explicit approval unless a later user instruction explicitly grants a run-specific exception.

Run-specific resource note:

- On 2026-03-20, the user raised the working RAM ceiling for the active run to `800 GB`

## Resolved Environment

Run the workflow with:

- Python: `/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python`
- Intended long-term env name from the runbook: `Scanpy_gdTmodel`
- Current practical resolution for this run: `rapids_sc_py310` is the working single-cell environment and should be used until `Scanpy_gdTmodel` is aligned to the same package set

Required packages confirmed in `rapids_sc_py310`:

- `anndata`
- `scanpy`
- `scvi`
- `torch`
- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
- `scipy`
- `h5py`

### Executable helpers used so far

- `phase0_dataset_audit.py`
- `repair_h5ad_from_raw.py`
- `phase1_extract_tnk_candidates.py`
- `phase1_finalize_from_temp.py`
- `phase1b_conservative_cleanup.py`
- `phase1c_replace_harmonized_metadata.py`
- `phase2_merged_cleanup.py`
- `phase3_scvi_scanvi.py`
- `phase4_gdt_module_scoring.py`
- `plot_phase4_threshold_barplots.py`
- `watch_h5ad_v2_and_resume.py`
- `tissue_correction_workflow.py`
- `disease_status_correction_workflow.py`
- This helper is the concrete implementation of the Phase 0 audit logic described below.
- `phase0_dataset_audit.py` accepts `--registry <csv>` when the canonical audit must be rerun against a repaired registry such as `h5ad_v2.csv`.
- `repair_h5ad_from_raw.py` is the concrete repair helper for Category B datasets where `adata.raw` contains recoverable integer-like counts for the current feature space.
- `phase1_extract_tnk_candidates.py` performs the Category A high-recall coarse T/NK extraction and writes one temp candidate H5AD per dataset.
- `phase1_finalize_from_temp.py` resumes from those temp candidate H5ADs, normalizes any dense `X` matrices to CSR, performs the on-disk concat into `TNK_candidates.h5ad`, validates the merged object, writes the Phase 1 QC tables and figures, and removes the temp directory after success.
- `phase1b_conservative_cleanup.py` performs the conservative first-pass cleanup, removes only obvious non-T/NK contaminants, applies the user-requested `<500 cells` gene filter, rewrites `TNK_candidates.h5ad` in place, and writes the Phase 1b QC package.
- `phase1c_replace_harmonized_metadata.py` exports the merged `adata.obs`, backs up the previous harmonized metadata CSV, joins the filtered candidate cells to the harmonized metadata by explicit string-typed `GSE + barcode`, writes the required `project name`, `sampleid`, and `barcodes` columns, validates uniqueness and row counts, and replaces the canonical metadata target atomically.
- `phase2_merged_cleanup.py` performs the merged-context second-pass cleanup, removes only high-confidence off-target or low-quality cells, reapplies the `>=500 cells` gene filter, writes `TNK_cleaned.h5ad`, and emits the Phase 2 QC package.
- `phase3_scvi_scanvi.py` attaches harmonized metadata, trains scVI, computes the integrated latent space and RAPIDS neighborhood/UMAP outputs, runs scANVI on a stratified subset of the cleaned milestone, transfers labels back to all cells in latent space by nearest-centroid assignment, and writes `integrated.h5ad` plus the Phase 3 QC package.
- `phase3_scvi_scanvi.py` writes a persistent recovery log to `Integrated_dataset/logs/phase3_run.log`, reuses the saved scVI checkpoint only when its HVG set matches the current exclusion-aware policy, excludes mitochondrial/ribosomal/noncoding-like genes from the HVG set used for clustering and UMAP, and releases CPU/GPU memory explicitly after the scVI latent and RAPIDS stages to reduce restart risk on the large merged milestone.
- `phase3_scvi_scanvi.py` now stages only the large H5AD workload under a mirrored local temp tree rooted at `/ssd/tnk_phase3/Integrated_dataset/` when that path is writable; tables, PNG figures, logs, scripts, and model artifacts remain on NFS.
- `phase3_scvi_scanvi.py` also maintains a stable NFS-side symlink view at `high_speed_temp/Integrated_dataset`, pointing to the mirrored SSD tree, so the fast working set remains visible from the project root.
- `phase3_scvi_scanvi.py` keeps the validated `integrated.h5ad` in the mirrored SSD tree and does not auto-migrate it back to NFS; later migration is a separate explicit step.
- `phase4_gdt_module_scoring.py` uses the package-faithful TRA/TRB/TRD module definitions, computes continuous scores from a temporary normalize-plus-log1p copy of count-space `X`, writes the score columns back into `integrated.h5ad`, and emits the Phase 4 QC package on NFS.
- `plot_phase4_threshold_barplots.py` renders descending PNG barplots from the exported Phase 4 threshold summary CSVs without re-reading the full integrated H5AD.
- `watch_h5ad_v2_and_resume.py` polls every 120 seconds for `h5ad_v2.csv`; once the repaired registry appears, it reruns Phase 0 against that registry and then rebuilds `TNK_candidates.h5ad` from the updated Category A set.
- `tissue_correction_workflow.py` is a post-integration metadata helper that samples tissue-related metadata, applies the committed rule set in `tissue_correction_rules.json`, exports a standalone `tissue_corrected` review table, and only writes back into `integrated.h5ad` when explicitly asked with `--write-h5ad`.
- `disease_status_correction_workflow.py` is a post-integration metadata helper that samples condition-related metadata, applies the committed rule set in `disease_status_correction_rules.json`, exports a standalone `disease_status_corrected` review table, and only writes back into `integrated.h5ad` when explicitly asked with `--write-h5ad`.
- The markdown file remains the canonical human-readable workflow; the helper exists to execute the Phase 0 audit reproducibly.

### CUDA runtime note for `rapids_sc_py310`

- In this env, `rapids_singlecell` must be imported only after `torch`
- Required import order for RAPIDS-backed plotting/integration helpers:
  1. `import torch`
  2. `import rapids_singlecell as rsc`
- This avoids the observed `libc10_cuda.so` / `cudaGetDriverEntryPointByVersion` symbol-resolution failure

## Canonical Inputs And Outputs

### Inputs

- Registry: `h5ad.csv`
- Dataset files: every `h5ad_path` listed in `h5ad.csv`

### Output root

- `Integrated_dataset/`

## Supplementary 10x 5' Intake Lane

The user approved one separate supplementary Phase 0-1 lane for six non-registry
10x 5' datasets:

- `GSE179994`
- `GSE235863`
- `GSE240865`
- `GSE287301`
- `GSE234069`
- `GSE287541`

Rules for this supplementary lane:

- process these six GSEs with the same Phase 0-1 logic used for the existing registry inputs
- store the per-GSE H5AD files in `downloads/per_gse_h5ad_with_metadata/`
- do not include any 10x 3' data
- for `GSE234069`, only use `downloads/GSE234069/suppl/10x_5/`
- when TCR is provided as a separate matrix or contig annotation rather than already embedded in metadata, process it according to `TCR_INTEGRATION_SOP.md` and the sample-aware workflow in `tcr_integration_workflow.json`
- standardize supplementary metadata and TCR outputs to the current `analysis_26GSE_V4/outputs/harmonized_metadata_v4.csv` column layout
- write the supplementary candidate milestone to `Integrated_dataset/TNK_candidates_supp.h5ad`
- do not merge `TNK_candidates_supp.h5ad` into `Integrated_dataset/TNK_candidates.h5ad` until the user reviews the supplementary Phase 0 and Phase 1 QC outputs and explicitly approves the merge
- do not advance to Phase 2 until that supplementary merge approval is granted

### Milestone H5AD files

- `Integrated_dataset/TNK_candidates.h5ad`
- `Integrated_dataset/TNK_candidates_supp.h5ad`
- `Integrated_dataset/TNK_cleaned.h5ad`
- `Integrated_dataset/integrated.h5ad`

### Tables

- `Integrated_dataset/tables/phase0_dataset_audit.csv`
- `Integrated_dataset/tables/phase0_category_summary.csv`
- `Integrated_dataset/tables/phase1_categoryA_selection_summary.csv`
- `Integrated_dataset/tables/phase1_categoryA_marker_availability.csv`
- phase-specific QC tables for later phases

### Logs

- `Integrated_dataset/logs/phase0_qc_summary.md`
- `Integrated_dataset/logs/phase1_qc_summary.md`
- phase-specific run logs and QC summaries for later phases

### Figures

- `Integrated_dataset/figures/phase0_category_distribution.png`
- `Integrated_dataset/figures/phase0_matrix_state_overview.png`
- `Integrated_dataset/figures/phase0_dataset_size_overview.png`
- `Integrated_dataset/figures/phase0_metadata_completeness.png`
- `Integrated_dataset/figures/phase1_categoryA_candidate_yield.png`
- `Integrated_dataset/figures/phase1_categoryA_candidate_fraction.png`
- `Integrated_dataset/figures/phase1_categoryA_marker_support.png`
- later phase figures must also be saved in `Integrated_dataset/figures/`

## Shared Config

```python
PYTHON_BIN = "/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python"
REGISTRY_CSV = "h5ad.csv"
OUTPUT_ROOT = "Integrated_dataset"
FIGURE_DIR = f"{OUTPUT_ROOT}/figures"
TABLE_DIR = f"{OUTPUT_ROOT}/tables"
LOG_DIR = f"{OUTPUT_ROOT}/logs"

PHASE0_AUDIT_CSV = f"{TABLE_DIR}/phase0_dataset_audit.csv"
PHASE0_CATEGORY_CSV = f"{TABLE_DIR}/phase0_category_summary.csv"
PHASE0_QC_MD = f"{LOG_DIR}/phase0_qc_summary.md"

TNK_CANDIDATES = f"{OUTPUT_ROOT}/TNK_candidates.h5ad"
TNK_CLEANED = f"{OUTPUT_ROOT}/TNK_cleaned.h5ad"
INTEGRATED = f"{OUTPUT_ROOT}/integrated.h5ad"

TARGET_TCR_GENES = [
    "TRAC", "TRBC1", "TRBC2", "TRDC", "TRGC1", "TRGC2",
    "CD3D", "CD3E", "CD3G", "NKG7", "KLRD1"
]
METADATA_HINTS = {
    "donor": ["donor", "patient", "subject", "individual"],
    "sample": ["sample", "orig.ident", "specimen", "biosample"],
    "library": ["library", "library_id", "batch", "lane", "channel"],
}
```

## Shared Utility Logic

### Phase 0 audit helpers

Implement helpers with clear docstrings and explicit validation:

- `load_registry()` reads `h5ad.csv`, checks required columns, and fails loudly on missing paths
- `inspect_h5ad_structure()` uses `h5py` first to inspect `X`, `layers`, `raw`, sparse encoding, and sampled numeric values without densifying full matrices
- `inspect_obs_var()` reads `obs` and `var` metadata without densifying the matrix and records duplicated genes, gene naming patterns, and metadata fields
- `detect_metadata_fields()` reports candidate donor/sample/library columns and whether at least one field exists for each category
- `classify_matrix_state()` labels each dataset as raw-like, recoverable-from-counts, recoverable-from-raw, normalized-only, scaled/integrated, or unreadable
- `repair_h5ad_from_raw()` restores recoverable counts from `adata.raw` onto the current `var` space, removes stale derived embeddings, validates the rewritten file, and replaces the original path atomically
- `assign_phase0_category()` maps each dataset to exactly one triage category:
  - Category A: raw-count eligible now
  - Category B: recoverable with explicit rescue logic
  - Category C: non-raw or non-recoverable for the count-based pipeline
- `write_phase0_outputs()` writes tables, Markdown QC summary, and 300 dpi PNG figures

### Safety rules

- Never densify a large sparse matrix just to inspect it
- Sample values from sparse `.data` arrays when possible
- Treat negative values in `X` as a strong scaled/integrated signal
- Prefer `layers["counts"]` over `raw` when both exist and look integer-like
- Record uncertainty explicitly instead of forcing a confident label

## Phase 0: Dataset Audit And Eligibility Triage

### Objective

Audit every registry dataset before any extraction or integration.

### Required fields in the audit table

- `gse_id`
- `h5ad_path`
- `source_root`
- `n_obs`
- `n_vars`
- `x_storage`
- `x_nnz`
- `x_density`
- `x_has_negative`
- `x_integer_like_fraction`
- `x_max_sample`
- `has_counts_layer`
- `counts_integer_like_fraction`
- `has_raw`
- `raw_integer_like_fraction`
- `state_label`
- `phase0_category`
- `category_reason`
- `duplicate_var_names`
- `gene_naming_style`
- `tcr_genes_present`
- `donor_fields`
- `sample_fields`
- `library_fields`
- `obs_columns_n`
- `read_error`

### Outputs

- `phase0_dataset_audit.csv`
- `phase0_category_summary.csv`
- `h5ad_without_raw_count.csv`
- four Phase 0 PNGs
- `phase0_qc_summary.md`
- `TNK_candidates.h5ad`
- `phase1_categoryA_selection_summary.csv`
- `phase1_categoryA_marker_availability.csv`
- three Phase 1 PNGs
- `phase1_qc_summary.md`

### QC acceptance criteria

- every dataset in `h5ad.csv` is represented exactly once
- every dataset has one category and one reason
- unreadable files are surfaced explicitly
- no Phase 1 extraction begins before user review

## Phase 1: Per-Dataset Coarse T/NK Extraction

### Entry criteria

- Phase 0 QC approved by the user
- in-scope datasets selected from Category A and any explicitly approved Category B rescues

### Behavior

- perform high-recall extraction using broad T/NK markers and metadata where available
- keep ambiguous lymphocyte-like cells rather than over-pruning
- write the merged candidate object to `Integrated_dataset/TNK_candidates.h5ad`
- generate PNG summaries of per-dataset yield, marker support, and coarse contaminant checks

### Stop condition

Stop after candidate merge and first-pass summaries for user QC.

## Phase 1b: Conservative First Cleanup

### Behavior

- remove only obvious non-T/NK contaminants and high-confidence doublets
- after conservative cell cleanup, remove genes expressed in fewer than 500 cells
- preserve rare or atypical γδT-like populations
- update `Integrated_dataset/TNK_candidates.h5ad` in place after validation
- generate cleanup summary tables and before/after figures

### Stop condition

Stop for user QC before Phase 1c.

## Phase 1c: Merged Metadata Backup And Replacement

### Entry criteria

- Phase 1 and 1b QC approved by the user

### Behavior

- export merged `adata.obs` from the current unified candidate object
- save a backup copy as `metadata.csv.bk`
- replace `analysis_26GSE_V4/outputs/harmonized_metadata_v4.csv`
- join the filtered merged object to the previous harmonized metadata by string-typed `GSE + barcode`
- use the join key:
  - `project name`
  - `sampleid`
  - `barcodes`
- prefer `sample_id` from the previous harmonized metadata for `sampleid`, and fall back to the merged Phase 1 sample label only when the previous metadata field is blank
- fail loudly if any join-key column is missing
- fail loudly if the combined join key is not unique
- fail loudly if replacement row counts do not match expectations or if rows are dropped/duplicated
- log the replacement result before any Phase 2 work

### Stop condition

Stop if metadata replacement validation fails. Do not proceed to Phase 2 until this step is completed and validated.

## Phase 2: Merged Cleanup

### Entry criteria

- Phase 1 and 1b QC approved by the user
- Phase 1c metadata backup/replacement completed and validated

### Behavior

- merge approved candidate datasets
- complete the merged metadata backup/replacement step required by the runbook
- recluster and perform stricter second-pass contaminant removal
- save validated output to `Integrated_dataset/TNK_cleaned.h5ad`
- generate contamination review, marker overlay, and cluster decision figures
- write QC figures in PNG format only

### Run-specific execution override

- On 2026-03-20, the user explicitly approved Phase 2 execution on the current merged `TNK_candidates.h5ad`
- For this run only, the Phase 2 QC package may be reviewed internally; if the QC result is satisfactory, execution may continue directly to Phase 3 without another user checkpoint
- The default policy still applies to later runs unless the user grants the same exception again

### Stop condition

Default: stop for user QC before Phase 3.

Exception for the current run: proceed directly if the internal Phase 2 QC review passes.

## Phase 3: scVI Integration With Optional Reference scANVI Annotation

### Entry criteria

- Phase 2 QC approved by the user, or an explicit run-specific user instruction authorizes internal Phase 2 QC and direct continuation
- `rapids_sc_py310` remains the active execution env unless replaced by a validated `Scanpy_gdTmodel`

### Behavior

- prepare batch-aware inputs from `TNK_cleaned.h5ad`
- run scVI with GPU enabled when available and useful
- avoid assumptions that require more than 400 GB RAM
- exclude mitochondrial genes, ribosomal genes, and noncoding-like RNA symbols from the HVG set used for clustering and UMAP
- if the user pauses scANVI, finish and save the scVI latent space, Leiden clustering, UMAP, and clustering-first QC package before resuming annotation
- after scVI integration, run scANVI-based reference annotation on a stratified subset for coarse T and NK labeling
- transfer the resulting labels back to all cells in the scVI latent space with a scalable centroid-based rule
- primary reference model path:
  - `/home/tanlikai/databank/owndata/fasting/raw/report_from_niuxian/models/census_scanvi_ref_v1`
- reference companion H5AD currently exposes:
  - label column `cell_type`
  - batch column `batch`
- validate query/reference compatibility before mapping
- treat this reference as a coarse annotation prior, not as final tissue-state truth or γδT subtype truth
- write both detailed transferred labels and collapsed T/NK superclass labels into `integrated.h5ad`
- if transferred labels fall outside T/NK space, flag them for review instead of auto-dropping them
- if Phase 3 QC shows that the scANVI result is messy or low-confidence, keep the scANVI fields and PNG/QC outputs for reference only, do not roll back `integrated.h5ad`, and use the scVI latent space plus Leiden/simple scVI-based annotation as the canonical downstream interpretation
- record the subset size, transfer method, and HVG exclusion summary in Phase 3 tables
- save validated output to `Integrated_dataset/integrated.h5ad`
- generate batch-mixing, Leiden, and donor/sample mixing figures first; add scANVI label and label-confidence figures when annotation is enabled
- write QC figures in PNG format only

### Stop condition

Stop for user QC before Phase 4.

Phase 3 QC may accept the integrated milestone while declining scANVI as the primary downstream annotation layer. In that case:

- keep `integrated.h5ad` as written
- keep scANVI labels in the object for reference only
- use scVI/Leiden/simple annotation for downstream analyses unless the user later requests a different annotation strategy

## Phase 4: TRAB/TRB/TRD Scoring And Threshold Summaries

### Entry criteria

- Phase 3 QC approved by the user
- `integrated.h5ad` exists and remains the canonical integrated milestone
- scANVI may remain in the object for reference, but Phase 4 should not depend on scANVI being accepted as the primary annotation layer

### Behavior

- use `phase4_gdt_module_scoring.py` as the execution helper
- keep the validated large H5AD on SSD and keep tables, logs, and PNG figures on NFS
- use the package-faithful TRA/TRB/TRD module definitions from `gdt_tcr_module_sharing_package_full`
- compute scores from a temporary normalize-plus-log1p view of count-space `X`
- write continuous score columns back into `integrated.h5ad`:
  - `phase4_tra_score`
  - `phase4_trb_score`
  - `phase4_trab_score`
  - `phase4_trd_score`
  - `phase4_trd_minus_trab`
- when requested, extend Phase 4 with scaled score columns and threshold-specific summary tables and plots
- keep Phase 4 figures in PNG format only
- do not force a hard gdT or abT call unless the user explicitly requests one

### Stop condition

Stop for user QC after the Phase 4 score outputs, threshold summaries, and PNG figures are written.

## Main Entrypoint

Use one obvious control flow:

```python
def main() -> None:
    ensure_output_dirs()
    run_phase0_audit()
    stop_for_user_qc("Phase 0")
    run_phase1_extraction()
    run_phase1b_cleanup()
    stop_for_user_qc("Phase 1/1b")
    run_phase1c_metadata_backup()
    run_phase2_cleanup()
    if phase2_requires_user_qc():
        stop_for_user_qc("Phase 2")
    run_phase3_scvi_and_scanvi()
    stop_for_user_qc("Phase 3")
    run_phase4_scoring()
    stop_for_user_qc("Phase 4")
```

The actual execution must respect QC gates. The control flow above is canonical, but each phase after Phase 0 should only be run after explicit approval unless the user explicitly grants a documented exception for the current run.
