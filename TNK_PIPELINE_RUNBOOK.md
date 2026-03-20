# TNK_PIPELINE_RUNBOOK.md
## Compute resource envelope

Allowed resources:
- CPU: up to 80 cores
- GPU: allowed
- RAM: up to 400 GB

Execution policy:
- optimize for robustness within 80 CPU / GPU / 400 GB RAM limits
- preserve sparse matrices whenever possible
- avoid accidental densification of large matrices
- do not design pipeline steps that assume 2 TB RAM
- use GPU for scVI and other steps where acceleration is clearly useful
- use CPU parallelism where appropriate, but do not exceed 80 cores
- when a step may create a major memory peak, prefer chunking, on-disk writes, checkpointing, and garbage collection
# Codex execution runbook for T/NK integration and γδT-focused scoring workflow

## 1. Project goal

Build a robust workflow to:
1. audit public single-cell datasets (paths: h5ad.csv)
2. identify which datasets are suitable for unified processing
3. extract a high-recall T/NK candidate pool
4. remove doublets and residual contaminants in staged rounds
5. integrate candidate T/NK cells using scVI where appropriate
6. perform Phase 4 scoring using `gdt_tcr_module_sharing_package_full`
7. generate high-quality figures and a stable integrated object

The current project stops at Phase 4.
Do not implement Phase 5 training/validation at this time.

---

## 2. Core principles

### 2.1 Human-readable code
Pipeline code must remain easy for a human to read, review, and modify.

### 2.2 Minimal file sprawl
Do not create many script versions or many milestone H5ADs.
Use git history instead of filename proliferation.

### 2.3 High recall first, then high precision
Early phases should avoid over-filtering away rare or atypical γδT cells.

### 2.4 Data state matters
Public H5AD files may contain:
- raw counts
- normalized/log-transformed values
- scaled/integrated data

These must not be mixed blindly.

### 2.5 Two-stage contamination control
Use staged cleanup:
- conservative early cleanup
- stricter later cleanup after merge / embedding / clustering

### 2.6 The runbook overrides memory
If memory and this file disagree, this file wins.

---

## 3. Canonical files

### Required control files
- `AGENTS.md`
- `TNK_PIPELINE_RUNBOOK.md`
- `TNK_PIPELINE_STATUS.md`

### Canonical pipeline script file
- `TNK_PHASES_1_3_SCRIPT.md`

This markdown file must contain the human-readable canonical script for:
- Phase 0 support logic if needed
- Phase 1
- Phase 1b
- Phase 2
- Phase 3

Do not fragment the main workflow into many parallel script versions.

---

## 4. Output structure

All generated outputs must be saved under:
- `Integrated_dataset/`

Recommended structure:
- `Integrated_dataset/TNK_candidates.h5ad`
- `Integrated_dataset/TNK_cleaned.h5ad`
- `Integrated_dataset/integrated.h5ad`
- `Integrated_dataset/figures/`
- `Integrated_dataset/tables/`
- `Integrated_dataset/logs/`

### Milestone H5AD policy
Keep only these milestone H5ADs unless there is a strong reason not to:
1. `TNK_candidates.h5ad`
2. `TNK_cleaned.h5ad`
3. `integrated.h5ad`

Rules:
- do not create many extra milestone H5ADs
- update milestone files rather than creating many alternatives
- temporary H5ADs are allowed during execution but should be removed after successful validation

### Phase-specific interpretation
- `TNK_candidates.h5ad`: merged candidate pool after Phase 1 / 1b
- `TNK_cleaned.h5ad`: cleaned object after Phase 2
- `integrated.h5ad`: integrated object after Phase 3, then updated in place during Phase 4

---

## 5. Environment and package policy

Use the conda environment:
- `rapids_sc_py310`

Rules:
- prefer conda installation
- use pip only if the package is unavailable via conda channels
- if pip is used, record why
- do not use base environment
- install all needed packages into `rapids_sc_py310`

This environment should contain all packages needed for:
- scanpy
- anndata
- scvi-tools
- plotting
- QC
- Phase 4 scoring with `gdt_tcr_module_sharing_package_full`

---

## 6. Git / GitHub policy

Use git frequently, but avoid visible script version clutter.

Required behavior:
- commit after each meaningful validated milestone
- commit before risky refactors
- commit at the end of a session if code changed
- push to GitHub regularly
- use descriptive commit messages

Do not preserve history by making filenames like:
- `_v2`
- `_final`
- `_final2`
- `_revised`

Use git history for version tracking.

---

## 7. Code readability rules

All scripts must be readable for humans.

Required style:
- one clear config section near the top
- one obvious `main()` entrypoint
- phase-separated sections
- concise helper functions
- clear docstrings
- comments explaining reasoning
- no unnecessary one-liners
- no deeply nested opaque logic when refactoring into helpers is clearer
- clear path variables
- consistent naming
- informative logging
- explicit validation checks
- explicit save/load steps

Markdown code blocks should be organized so a human can quickly locate:
- imports
- config
- utility functions
- Phase 1
- Phase 1b
- Phase 2
- Phase 3
- figure generation
- save functions
- CLI/main block

---

## 8. Mandatory file re-reading behavior

At the start of every session:
1. read this file
2. read `TNK_PIPELINE_STATUS.md`
4. restate the current milestone and next action

Re-read these files whenever:
- a milestone changes
- context appears compressed or summarized
- the task changes substantially
- uncertainty appears
- the plan seems fuzzy
- output naming becomes uncertain
- environment or package behavior becomes uncertain

---

## 9. Workflow overview

The project should be executed in these phases:

### Phase 0: Data audit and eligibility triage
### Phase 1: Per-dataset coarse T/NK candidate extraction
### Phase 1b: Conservative first-pass doublet / obvious contaminant removal
### Phase 2: Merged candidate cleanup and second-pass purification
### Phase 3: scVI integration and post-integration scANVI annotation
### Phase 4: TRAB/TRB scoring using `gdt_tcr_module_sharing_package_full`

Do not work on Phase 5 now.

---

## 10. Phase 0: Data audit and eligibility triage

### Objective
Determine the processing state and suitability of each public dataset before it enters the main pipeline.

### Required checks per dataset
Inspect:
- shape
- number of cells
- number of genes
- matrix sparsity
- whether `adata.X` appears to be raw counts
- whether `adata.layers['counts']` exists
- whether `adata.raw` exists
- whether values are integer-like
- whether negative values exist
- likely raw / normalized / scaled state
- metadata fields
- donor/sample/library identifiers
- duplicated var names
- gene naming conventions
- presence of TCR genes

### Triage categories
- Category A: raw-count eligible
- Category B: recoverable
- Category C: non-raw / non-recoverable

Only Category A should enter the main count-based integration workflow by default.

---

## 11. Phase 1: Per-dataset coarse T/NK candidate extraction

### Objective
Build a high-recall candidate pool of T/NK cells from each eligible dataset.


### Marker groups
Use permissive T/NK marker logic and a broader myeloid exclusion/penalty panel.
Do not over-prune NK-like or cytotoxic γδT-like cells here.
### Phase 1 biological extraction logic

The candidate T/NK pool must be built using explicit marker-guided logic rather than generic lymphocyte heuristics.

#### Safe T/NK genes
These genes support permissive retention of candidate T/NK cells:
- CD3D
- CD3E
- CD3G
- TRAC
- TRBC1
- TRBC2
- TRDC
- TRGC1
- TRGC2
- CD8A
- CD8B
- NKG7
- GNLY
- NCR1

#### Conditional genes
These genes may support retention only when myeloid contamination is not supported:
- CD4
- NCAM1
- FCGR3A

#### Myeloid exclusion / penalty genes
Do not rely only on LYZ and CD14. Use a broader myeloid penalty panel such as:
- LYZ
- CD14
- FCER1G
- TYROBP
- LST1
- AIF1
- CTSS
- SAT1
- CST3
- SPI1
- LGALS3
- VCAN
- S100A8
- S100A9

#### Phase 1 policy
- prioritize high recall for candidate T/NK cells
- do not over-prune NK-like or cytotoxic γδT-like cells
- use myeloid markers as exclusion/penalty context
- preserve rare γδT-compatible states for later cleanup
---

## 12. Phase 1b: Conservative first-pass doublet and obvious contaminant removal

### Objective
Remove only high-confidence junk before merged analysis.

### Rule
Standard count-based doublet detection should only be performed where raw counts are available.

### What to remove
Only high-confidence cases such as:
- extreme outliers with lineage conflict
- obvious mixed-lineage cells
- obvious T + myeloid mixed cells
- obvious lymphoid + epithelial mixed cells

### Principle
Remove obvious junk, not all suspicious cells.
## 13. Phase 1c: Merged metadata backup and replacement

### Objective
Immediately after all samples are merged together into the unified object, export and update the integrated harmonized metadata before moving to the next phase.

### Required outputs
- backup file: `metadata.csv.bk`
- integrated metadata target: `analysis_26GSE_V4/outputs/harmonized_metadata_v4.csv`

### Required join key
Use the following columns together as the unique join key:
- `project name`
- `sampleid`
- `barcodes`

### Mandatory procedure
1. export merged `adata.obs`
2. save a backup copy as `metadata.csv.bk`
3. validate that the join-key columns exist
4. validate that the combined join key is unique per row
5. replace `analysis_26GSE_V4/outputs/harmonized_metadata_v4.csv`
6. validate that row counts match expectations and that no rows were silently dropped or duplicated
7. log the replacement result

### Rules
- do not guess missing join-key columns
- do not continue if the join key is not unique
- do not continue if row counts are inconsistent
- do not proceed to Phase 2 until this step is completed and validated
---

## 14. Phase 2: Merged candidate cleanup and second-pass purification

### Objective
After candidate pooling, do a stricter cleanup using more context.

Use:
- doublet scores where valid
- n_counts / n_genes
- contamination scores
- cluster markers
- lineage-conflict markers



Goals:
- remove residual doublets
- remove residual myeloid contaminants
- remove obvious non-T/NK clusters
- preserve plausible γδT and NK-like γδT populations

### Phase 2 cleanup logic
Phase 2 is the main high-precision purification stage.

Use merged-context evidence such as:
- doublet scores where valid
- n_counts
- n_genes
- mitochondrial fraction if relevant
- myeloid contamination score
- T-cell score
- NK score
- cluster markers
- lineage-conflict markers

### Doublet policy
Use two-stage doublet handling:
- Phase 1b: conservative removal of only high-confidence junk
- Phase 2: stricter cleanup using merged context

Do not aggressively remove candidate γδT-like states in Phase 1b.

### Residual contaminant review
Potential contaminant clusters may express:
- myeloid: LYZ, FCER1G, CTSS, AIF1, LST1, SAT1, CST3
- B cell: MS4A1, CD79A, CD74, HLA-DRA
- epithelial: EPCAM, KRT8, KRT18, KRT19
- erythroid: HBB, HBA1, HBA2
- platelet: PF4, PPBP

Do not remove a cluster solely because it expresses:
- NKG7
- GNLY
- PRF1
- GZMB
- FCGR3A
- NCAM1

These may occur in real NK-like γδT or cytotoxic T/NK-compatible states.

### Output


Save result as:
- `Integrated_dataset/TNK_cleaned.h5ad`

---

## 14. Phase 3: scVI integration

### Objective
Build a stable integrated latent space with scVI, then perform cautious scANVI-based
reference annotation for coarse T and NK labeling.

### Eligibility
Only use datasets appropriate for unified count-based integration.


### Setup principles
- use raw counts
- keep counts explicitly in a counts layer if possible
- preserve sparse matrices
- avoid accidental densification
- use a batch key finer than GSE_ID when available
### Feature handling for integration
The main purpose of Phase 3 is stable integration, not forced preservation of every lineage marker in the HVG set.

#### Default rule
Do not force-rescue genes into the scVI integration feature set by default.

#### Instead
- keep biologically important genes in the full object
- run standard HVG selection for integration
- use key lineage and TCR genes later for annotation, scoring, and interpretation

#### Genes that must be preserved in the full object for downstream work
- TRDC
- TRGC1
- TRGC2
- TRAC
- TRBC1
- TRBC2
- CD3D
- CD3E
- CD3G
- CD4
- CD8A
- CD8B
- NCR1
- NCAM1
- FCGR3A

#### TCR variable genes
If present, preserve TRDV / TRGV / TRAV / TRBV genes in the object for downstream interpretation.
Do not blindly force all variable-region genes into the scVI integration feature set unless there is a clear first-pass reason.

### Integration setup principles
- use raw counts
- keep counts explicitly in a counts layer if possible
- preserve sparse matrices
- avoid accidental densification
- prefer a batch key finer than GSE_ID when available
- do not design steps that assume more than 400 GB RAM
- GPU may be used where beneficial
- CPU parallelism may use up to 80 cores

### Output
Save integrated object as:
- `Integrated_dataset/integrated.h5ad`

### Post-scVI scANVI annotation

After the scVI latent space is validated, perform a reference-guided scANVI
annotation step focused on T/NK interpretation.

#### Reference source
Primary candidate reference path:
- `/home/tanlikai/databank/owndata/fasting/raw/report_from_niuxian/models/census_scanvi_ref_v1`

Available reference artifacts:
- `model.pt`
- `census_reference_subset.h5ad`

Reference fields observed in the companion reference H5AD:
- label column: `cell_type`
- batch column: `batch`

#### Critical caution
Do not treat this reference as authoritative for tissue-state biology or for fine
γδT subtyping.

Rules:
- use the reference primarily for coarse T-versus-NK annotation support
- treat transferred non-T/NK labels as contamination-review flags, not as automatic truth
- do not let blood-biased or broad census labels overwrite tissue-specific biology
- do not use scANVI output alone as the final γδT identity call

#### Recommended execution plan
1. train scVI on the approved Phase 2-cleaned object
2. save the integrated latent representation to `Integrated_dataset/integrated.h5ad`
3. validate the reference model and query compatibility before scANVI mapping
4. run scANVI label transfer using the reference model if compatibility checks pass
5. store detailed transferred labels plus a collapsed T/NK superclass label
6. compare scANVI labels against marker expression and cluster structure
7. flag off-target labels for manual review rather than auto-dropping them

#### Required annotation outputs
Write annotation results into `Integrated_dataset/integrated.h5ad`, for example:
- detailed transferred label
- transferred label confidence / probability
- collapsed superclass label such as `T_cell`, `NK_cell`, or `reference_other`
- annotation source / method tag

#### Required Phase 3 QC figures
In addition to batch-mixing and marker plots, include:
- UMAP colored by scANVI detailed label
- UMAP colored by collapsed T/NK label
- confidence distribution for transferred labels
- marker-versus-scANVI agreement plots for core T and NK markers

Then generate high-quality PNG figures summarizing:
- UMAP
- batch mixing
- candidate population structure
- marker overlays
- contamination review plots

Save figures to:
- `Integrated_dataset/figures/`

---

## 15. Phase 4: scoring

### Objective
Perform TRAB/TRB scoring using:
- `gdt_tcr_module_sharing_package_full`

### Rules
- use `Integrated_dataset/integrated.h5ad` as the canonical input
- write Phase 4 score results back into the same integrated object when appropriate
- do not create many extra H5AD outputs unless explicitly needed
- generate high-quality PNG figures for Phase 4 scoring summaries

Examples of outputs:
- score distribution plots
- UMAP overlays of scores
- cluster-level summary plots

### Phase 4 scoring rule
Use the package:
- `gdt_tcr_module_sharing_package_full`



Do not create many extra H5AD outputs unless explicitly necessary.

Generate high-quality PNG figures for:
- score distributions
- UMAP score overlays
- cluster-level summaries
- any key comparison plots needed to evaluate scoring behavior

Save updated object as:
- `Integrated_dataset/integrated.h5ad`

Save figures to:
- `Integrated_dataset/figures/`

---

## 16. Figure requirements

All evaluation figures must be high-quality PNG.

Required standards:
- PNG format
- at least 300 dpi
- readable labels and legends
- informative titles
- consistent naming
- no tiny default fonts
- save under `Integrated_dataset/figures/`

Figures are required, not optional.

---

## 17. Required status updates

After each meaningful milestone, update `TNK_PIPELINE_STATUS.md` with:
- milestone completed
- files changed
- outputs created
- validations performed
- figures generated
- git commit/push status
- unresolved issues
- next action

---

## 18. If memory seems bad

If context feels compressed or partially forgotten:
1. stop coding
2. re-read this file
3. re-read `TNK_PIPELINE_STATUS.md`
4. re-read `TNK_PHASES_1_3_SCRIPT.md`
5. restate the milestone
6. continue only after that
