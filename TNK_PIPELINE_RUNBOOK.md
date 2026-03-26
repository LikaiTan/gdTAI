# TNK_PIPELINE_RUNBOOK.md

## 1. Purpose

This runbook is the canonical rule file for the T/NK integration and γδT
scoring workflow.

It defines:

- file-role boundaries
- default execution rules
- active run-specific exceptions
- canonical environment and path policy

It does not act as a changelog, task board, or run diary.

## 2. Control-file roles

Use the markdown files with these exact roles:

- `AGENTS.md`
  - concise agent-behavior rules
- `TNK_PIPELINE_RUNBOOK.md`
  - canonical rules
  - canonical exception mechanism
  - active run-specific exceptions
- `TNK_PHASES_0_4_SCRIPT.md`
  - canonical executable workflow for Phases 0-4
- `TNK_PIPELINE_STATUS.md`
  - current milestone
  - next action
  - blockers
  - currently active exceptions already defined here
- `CHANGELOG.md`
  - historical completed milestones
- `DECISIONS.md`
  - approved key decisions and retired historical exceptions
- `OUTPUTS.md`
  - canonical output inventory and current review packages

If two files disagree:

1. this runbook wins for rules
2. the status file wins for current state
3. the canonical workflow file wins for the standard phase logic

## 3. Project scope

Current project scope ends at:

- Phase 0: dataset audit
- Phase 1: coarse T/NK extraction
- Phase 1b: conservative first cleanup
- Phase 1c: merged metadata backup and replacement
- Phase 2: merged cleanup
- Phase 3: scVI integration
- Phase 4: TRAB/TRB/TRD scoring

Do not implement Phase 5 model training or validation unless the user
explicitly requests it.

## 4. Compute and environment policy

Default resource envelope:

- CPU: up to 80 cores
- GPU: allowed
- RAM: up to 400 GB

Canonical current execution environment:

- `rapids_sc_py310`

Policy:

- preserve sparse matrices whenever possible
- avoid accidental densification
- do not design steps that assume 2 TB RAM
- prefer GPU where it clearly helps
- prefer chunking, checkpointing, and on-disk staging over memory spikes

`Scanpy_gdTmodel` is not the active environment for this run. Treat it only as a
future alias plan unless this runbook is explicitly updated.

## 5. Output and path policy

Canonical output root:

- `Integrated_dataset/`

Canonical milestone H5AD files:

1. `Integrated_dataset/TNK_candidates.h5ad`
2. `Integrated_dataset/TNK_cleaned.h5ad`
3. `Integrated_dataset/integrated.h5ad`

Allowed milestone exception class:

- explicitly approved temporary or supplementary milestone H5AD files required
  by a bounded workflow step

Current approved milestone exception:

- `Integrated_dataset/TNK_candidates_supp.h5ad`
  - supplementary 10x 5' intake lane only

Canonical path policy:

- NFS `Integrated_dataset/` remains canonical for scripts, tables, logs, PNG
  figures, and model artifacts
- large H5AD files may use a mirrored SSD tree when active exceptions below say
  so

## 6. Default QC-gate rule

Default rule:

- every phase transition requires user-reviewed QC and explicit user approval

Required before phase advancement:

1. complete the planned QC checks
2. generate required tables, logs, and PNG figures
3. summarize findings, issues, and uncertainties
4. present the QC conclusion to the user
5. wait for explicit approval

After full sample merging, complete the merged metadata backup/replacement step
before moving forward.

## 7. Active run-specific exceptions

Only this section may define active run-specific exceptions.

### Resource exception

- active RAM ceiling override for this run: `800 GB`

### Path exception

- large-H5AD mirrored SSD tree is active for this run
- active mirrored path root: `/ssd/tnk_phase3/Integrated_dataset/`
- stable NFS-side symlink view: `high_speed_temp/Integrated_dataset`
- keep validated large H5AD files on the mirrored tree until explicit user
  migration approval

### QC-gate exception

- none currently active

Historical exceptions belong in `DECISIONS.md`, not here.

## 8. Canonical startup and re-read protocol

At the start of every session:

1. read this file
2. read `TNK_PIPELINE_STATUS.md`
3. read `TNK_PHASES_0_4_SCRIPT.md`
4. restate the current milestone and next action

Re-read the same three files whenever:

- the milestone changes
- the task changes substantially
- context is compressed or uncertain
- a long pause occurred
- output naming or path behavior becomes uncertain
- environment assumptions become uncertain

## 9. Git policy

- commit after meaningful validated milestones
- commit before risky refactors
- commit at session end if code changed
- push regularly
- use descriptive commit messages
- use git history instead of filename proliferation
