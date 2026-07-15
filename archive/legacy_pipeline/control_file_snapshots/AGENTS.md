# AGENTS.md
- Only Speak English To me Please.
- Allowed compute resources:
  - up to 80 CPU cores
  - GPU allowed
  - up to 400 GB RAM
- Pipeline design should not assume 2 TB RAM
- Confirm environment/package installation plan in `rapids_sc_py310`
## Phase-gate QC policy

Every phase transition requires a user-reviewed QC gate.

Before advancing from any phase to the next phase, you must:
1. complete the planned QC checks for the current phase
2. generate the required summary outputs, tables, and PNG figures
3. summarize the main findings, issues, and uncertainties
4. present the QC conclusion to the user
5. wait for explicit user approval before proceeding
6. do not read files and paths I did not tell you to read (you can apply for ) 

Rules:
- never proceed to the next phase without QC review
- never assume silent approval
- if QC is ambiguous or unsatisfactory, do not proceed
- explain what failed, what is uncertain, and what should be adjusted
- remain in the current phase until the user approves advancement
- After full sample merging, complete the merged metadata backup/replacement step defined in `TNK_PIPELINE_RUNBOOK.md` before moving to the next phase.


This applies to:
- Phase 0 -> Phase 1
- Phase 1 / 1b -> Phase 2
- Phase 2 -> Phase 3
- Phase 3 -> Phase 4
## Project scope
This repository builds a large-scale single-cell T/NK integration workflow across  public datasets, with the downstream goal of identifying γδ T cells.

Current scope ends at:
- Phase 0: dataset audit
- Phase 1: coarse T/NK extraction
- Phase 1b: conservative first cleanup
- Phase 2: merged cleanup
- Phase 3: scVI integration
- Phase 4: TRAB/TRB scoring using `gdt_tcr_module_sharing_package_full`

Do not work on model training/validation (former Phase 5) unless explicitly requested later.

---

## Mandatory startup behavior
Before doing any work:
1. Read `TNK_PIPELINE_RUNBOOK.md`
2. Read `TNK_PIPELINE_STATUS.md` if it exists
3. Read `TNK_PHASES_1_3_SCRIPT.md`
4. Restate the current milestone and next action before editing code

If there is uncertainty, memory compression, or context drift:
- Re-read `TNK_PIPELINE_RUNBOOK.md`
- Re-read `TNK_PIPELINE_STATUS.md`
- Then continue

The runbook is the source of truth for execution rules.
The status file is the source of truth for current progress.
The phase script markdown is the source of truth for the canonical Phase 1–3 pipeline.

---

## Model policy
- Main working thread: `gpt-5.4`
- Subagents: `gpt-5.4-mini`

Use `gpt-5.4` for:
- planning
- pipeline architecture
- biological interpretation
- difficult debugging
- final merge decisions
- editing the canonical Phase 1–3 script
- decisions about filtering, saving, and evaluation

Use `gpt-5.4-mini` for:
- bounded subagents
- dataset-by-dataset audit
- metadata inspection
- QC summary review
- log inspection
- helper summaries
- figure checklist review

If a subagent finds ambiguity or conflict, escalate to the main `gpt-5.4` thread before changing core logic.

---

## Subagent policy
Do not use subagents by default.

Use subagents only for bounded, parallelizable, mostly read-heavy tasks, such as:
- auditing one dataset or one dataset family
- checking whether an H5AD contains raw counts vs normalized/scaled data
- reviewing metadata consistency
- summarizing QC metrics
- reviewing failed logs
- reviewing candidate contaminant clusters

Do not use subagents for:
- editing the canonical Phase 1–3 script
- final architecture decisions
- final biological judgment
- repeated micro-edits to one file

Each subagent must:
- answer one narrow question
- return a concise conclusion
- provide evidence
- avoid editing core pipeline files unless explicitly requested

---

## Code readability policy
All scripts must be written for human readability first.

Required style:
- use clear function names
- add docstrings to all major functions
- keep a short config section at the top
- group paths and parameters in one place
- comment on why a step is done, not only what it does
- avoid overly clever or compressed code
- avoid deeply nested logic when a helper function is clearer
- keep sections clearly labeled by phase
- use one obvious `main()` entrypoint
- use informative logging messages
- fail loudly on invalid assumptions

Do not create many script versions such as:
- `phase1_v2.py`
- `phase1_final.py`
- `phase1_final_revised.py`

Maintain one canonical Phase 1–3 markdown file and use git history instead of versioned filenames.

---

## Script organization policy
The script for Phases 1–3 must live in one canonical markdown file:

- `TNK_PHASES_1_3_SCRIPT.md`

This file should contain the human-readable, end-to-end pipeline for:
- Phase 0 support functions if needed
- Phase 1
- Phase 1b
_Phase 1c
- Phase 2
- Phase 3
- figure generation relevant to these phases

Do not split the core Phase 1–3 workflow across many script files unless explicitly necessary.

---

## Output policy
All generated outputs must be saved under:

- `Integrated_dataset/`

Use only these milestone H5AD files unless an exception is explicitly justified:
1. `Integrated_dataset/TNK_candidates.h5ad`
2. `Integrated_dataset/TNK_cleaned.h5ad`
3. `Integrated_dataset/integrated.h5ad`

Rules:
- do not generate many intermediate milestone H5AD files
- update the canonical milestone files in place when appropriate
- use git for script history, not filename proliferation
- temporary files may be used during execution but should be deleted after successful validation

Figures must be saved under:
- `Integrated_dataset/figures/`

Tables/logs may be saved under:
- `Integrated_dataset/tables/`
- `Integrated_dataset/logs/`

---

## Figure policy
Evaluation figures must be high-quality PNG.

Requirements:
- PNG format
- minimum 300 dpi
- publication-style labeling
- readable font sizes
- consistent naming
- clear legends
- no tiny default text
- save all evaluation figures to `Integrated_dataset/figures/`

Examples:
- QC summaries
- batch mixing summaries
- UMAPs
- marker overlays
- contamination score plots
- phase comparison figures

---

## Environment policy
Use the conda environment:

- `Scanpy_gdTmodel`

Install packages there.
Prefer conda installation.
Use pip only when a package is unavailable through conda channels, and record the reason.

Do not install project dependencies into base or unrelated environments.

---

## Git / GitHub policy
Use git frequently, but do not create many file versions.

Required behavior:
- commit after each meaningful validated milestone
- commit before risky refactors
- commit at the end of a session if code changed
- push commits to GitHub regularly
- use clear commit messages
- rely on git history rather than versioned filenames

Do not create many duplicate scripts just to preserve history.

---

## Execution discipline
For each milestone:
1. Read runbook, status, and canonical phase script
2. State the current milestone
3. Decide whether subagents are justified
4. Do the work
5. Validate outputs
6. Save outputs to `Integrated_dataset/`
7. Generate/update required PNG figures
8. Commit and push if the milestone is meaningful
9. Update `TNK_PIPELINE_STATUS.md`

---

## Memory drift protection
After any of the following, re-read the runbook, status file, and phase script:
- long pause
- context summary/compression
- milestone change
- task switch
- uncertainty about filters, outputs, figures, or environment
- more than ~12 turns on the same task

Do not trust compressed memory over the markdown files.