# AGENTS.md

- Only speak English to the user.
- Default compute envelope:
  - up to 80 CPU cores
  - GPU allowed
  - up to 400 GB RAM
- Do not design pipeline steps that assume 2 TB RAM.

## Role Of This File

This file is a short agent-behavior contract.

It should stay concise.
Project-specific operational rules belong in:

- `TNK_PIPELINE_RUNBOOK.md`

Canonical workflow logic belongs in:

- `TNK_PHASES_0_4_SCRIPT.md`

Current state belongs in:

- `TNK_PIPELINE_STATUS.md`

## Required Startup Behavior

Before doing work:

1. read `TNK_PIPELINE_RUNBOOK.md`
2. read `TNK_PIPELINE_STATUS.md`
3. read `TNK_PHASES_0_4_SCRIPT.md`
4. restate the current milestone and next action before editing code

Re-read the same three files after:

- a long pause
- context compression or uncertainty
- a milestone change
- a substantial task switch

## Core Behavior Rules

- Follow `TNK_PIPELINE_RUNBOOK.md` for project-specific execution rules and active exceptions.
- Treat `TNK_PIPELINE_STATUS.md` as the source of truth for the current milestone and next action.
- Treat `TNK_PHASES_0_4_SCRIPT.md` as the source of truth for the canonical workflow.
- Do not invent new exceptions in status or workflow files.
- Do not proceed across phase gates unless the runbook says an exception is currently active.
- Complete the merged metadata backup/replacement step before moving past merged-candidate phases.

## Canonical Environment

Current canonical execution environment:

- `rapids_sc_py310`

`Scanpy_gdTmodel` is not the active environment for this run unless the runbook
is explicitly updated.

## Output Discipline

- Use `Integrated_dataset/` as the canonical NFS output root.
- Respect the runbook path policy for any mirrored SSD large-H5AD exception.
- Keep figures as high-quality PNG in `Integrated_dataset/figures/`.

## Workflow Documentation Rule

Keep `TNK_PHASES_0_4_SCRIPT.md` updated with, for each major phase or task:

- phase or task name
- exact `.py` script used
- key outputs

## Git And File Hygiene

- Commit after meaningful validated milestones and before risky refactors.
- Push regularly.
- Use git history instead of versioned filename sprawl such as `_v2`, `_final`, or `_revised`.
