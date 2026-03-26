# AGENTS.md

- Only speak English to the user.
- Allowed compute resources by default:
  - up to 80 CPU cores
  - GPU allowed
  - up to 400 GB RAM
- Pipeline design must not assume 2 TB RAM.

## Canonical execution environment

Current canonical execution environment:

- `rapids_sc_py310`

`Scanpy_gdTmodel` is not the active environment for this run. Treat it only as a
future alias plan until the project explicitly switches over.

## Control-file truth hierarchy

Use the markdown files with these exact roles:

- `TNK_PIPELINE_RUNBOOK.md`
  - canonical rules
  - canonical exception policy
  - active run-specific exceptions
- `TNK_PHASES_0_4_SCRIPT.md`
  - canonical executable workflow for Phases 0-4
- `TNK_PIPELINE_STATUS.md`
  - current milestone
  - next action
  - active blockers
  - active exceptions already in force
- `CHANGELOG.md`
  - historical completed milestones
- `DECISIONS.md`
  - approved key decisions and historical exceptions
- `OUTPUTS.md`
  - canonical output inventory and review-artifact index

`TNK_PIPELINE_RUNBOOK.md` is the source of truth for rules.
`TNK_PIPELINE_STATUS.md` is the source of truth for the current state.
`TNK_PHASES_0_4_SCRIPT.md` is the source of truth for the standard workflow.

## Mandatory startup behavior

Before doing work:

1. Read `TNK_PIPELINE_RUNBOOK.md`
2. Read `TNK_PIPELINE_STATUS.md`
3. Read `TNK_PHASES_0_4_SCRIPT.md`
4. Restate the current milestone and next action before editing code

If there is uncertainty, context drift, or a long pause:

1. Re-read `TNK_PIPELINE_RUNBOOK.md`
2. Re-read `TNK_PIPELINE_STATUS.md`
3. Re-read `TNK_PHASES_0_4_SCRIPT.md`

## Phase-gate QC policy

Default rule:

- every phase transition requires a user-reviewed QC gate
- do not proceed to the next phase without explicit user approval

This means:

1. complete the planned QC checks for the current phase
2. generate the required summary outputs, tables, and PNG figures
3. summarize the main findings, issues, and uncertainties
4. present the QC conclusion to the user
5. wait for explicit user approval before proceeding

Exceptions:

- only follow a no-stop or internal-continuation rule if that exception is
  explicitly listed as active in `TNK_PIPELINE_RUNBOOK.md`
- `TNK_PIPELINE_STATUS.md` may record that the active exception is currently in
  use, but it must not define new exceptions on its own

After full sample merging, complete the merged metadata backup/replacement step
defined in the runbook before moving to the next phase.

## Project scope

This repository covers:

- Phase 0: dataset audit
- Phase 1: coarse T/NK extraction
- Phase 1b: conservative first cleanup
- Phase 1c: merged metadata backup and replacement
- Phase 2: merged cleanup
- Phase 3: scVI integration
- Phase 4: TRAB/TRB/TRD scoring using `gdt_tcr_module_sharing_package_full`

Do not work on model training or validation unless the user explicitly asks.

## Output policy

All generated outputs must live under:

- `Integrated_dataset/`

Canonical milestone H5AD files:

1. `Integrated_dataset/TNK_candidates.h5ad`
2. `Integrated_dataset/TNK_cleaned.h5ad`
3. `Integrated_dataset/integrated.h5ad`

Allowed exception class:

- explicitly approved temporary or supplementary milestone H5AD files documented
  in `TNK_PIPELINE_RUNBOOK.md`

Figures:

- save all evaluation figures as high-quality PNG to `Integrated_dataset/figures/`

Tables and logs:

- `Integrated_dataset/tables/`
- `Integrated_dataset/logs/`

## Path policy

Canonical NFS root:

- `Integrated_dataset/`

Large-H5AD exception:

- if the runbook and status files say the mirrored SSD tree is active, use it
  only for large H5AD inputs and outputs
- keep scripts, tables, PNG figures, logs, and model artifacts on NFS

## Git policy

- commit after meaningful validated milestones
- commit before risky refactors
- commit at the end of a session if code changed
- push regularly
- use git history instead of versioned filename sprawl

## Code organization

- keep one canonical workflow markdown file: `TNK_PHASES_0_4_SCRIPT.md`
- do not create drifting filenames such as `_v2`, `_final`, `_revised`
- write scripts for human readability first
