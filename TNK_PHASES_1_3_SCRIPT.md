# TNK Phase 0-3 Canonical Workflow

## Purpose

This file is the single human-readable canonical workflow for:

- Phase 0: dataset audit and eligibility triage
- Phase 1: coarse T/NK extraction
- Phase 1b: conservative first cleanup
- Phase 1c: merged metadata backup and replacement
- Phase 2: merged cleanup
- Phase 3: scVI integration

The workflow is QC-gated. No phase transition is allowed without user review and explicit approval.

## Resolved Environment

Run the workflow with:

- Python: `/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python`
- Intended long-term env name from the runbook: `Scanpy_gdTmodel`
- Current practical resolution: `rapids_sc_py310` is the working single-cell environment and should be used until `Scanpy_gdTmodel` exists or is aliased to the same package set

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
- `watch_h5ad_v2_and_resume.py`
- This helper is the concrete implementation of the Phase 0 audit logic described below.
- `phase0_dataset_audit.py` accepts `--registry <csv>` when the canonical audit must be rerun against a repaired registry such as `h5ad_v2.csv`.
- `repair_h5ad_from_raw.py` is the concrete repair helper for Category B datasets where `adata.raw` contains recoverable integer-like counts for the current feature space.
- `phase1_extract_tnk_candidates.py` performs the Category A high-recall coarse T/NK extraction and writes one temp candidate H5AD per dataset.
- `phase1_finalize_from_temp.py` resumes from those temp candidate H5ADs, normalizes any dense `X` matrices to CSR, performs the on-disk concat into `TNK_candidates.h5ad`, validates the merged object, writes the Phase 1 QC tables and figures, and removes the temp directory after success.
- `phase1b_conservative_cleanup.py` performs the conservative first-pass cleanup, removes only obvious non-T/NK contaminants, applies the user-requested `<500 cells` gene filter, rewrites `TNK_candidates.h5ad` in place, and writes the Phase 1b QC package.
- `watch_h5ad_v2_and_resume.py` polls every 120 seconds for `h5ad_v2.csv`; once the repaired registry appears, it reruns Phase 0 against that registry and then rebuilds `TNK_candidates.h5ad` from the updated Category A set.
- The markdown file remains the canonical human-readable workflow; the helper exists to execute the Phase 0 audit reproducibly.

## Canonical Inputs And Outputs

### Inputs

- Registry: `h5ad.csv`
- Dataset files: every `h5ad_path` listed in `h5ad.csv`

### Output root

- `Integrated_dataset/`

### Milestone H5AD files

- `Integrated_dataset/TNK_candidates.h5ad`
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
- use the join key:
  - `project name`
  - `sampleid`
  - `barcodes`
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

### Stop condition

Stop for user QC before Phase 3.

## Phase 3: scVI Integration

### Entry criteria

- Phase 2 QC approved by the user
- `rapids_sc_py310` remains the active execution env unless replaced by a validated `Scanpy_gdTmodel`

### Behavior

- prepare batch-aware inputs from `TNK_cleaned.h5ad`
- run scVI with GPU enabled when available and useful
- avoid assumptions that require more than 400 GB RAM
- save validated output to `Integrated_dataset/integrated.h5ad`
- generate batch-mixing, marker-retention, and donor/sample mixing figures

### Stop condition

Stop for user QC before Phase 4.

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
    stop_for_user_qc("Phase 2")
    run_phase3_scvi()
    stop_for_user_qc("Phase 3")
```

The actual execution must respect QC gates. The control flow above is canonical, but each phase after Phase 0 should only be run after explicit approval.
