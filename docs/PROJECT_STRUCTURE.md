# Project Structure

## Control plane

The root contains only high-level project controls. `TNK_PIPELINE_RUNBOOK.md`
defines rules, `TNK_PIPELINE_STATUS.md` records the current milestone, and
`TNK_PHASES_0_4_SCRIPT.md` records the canonical executable workflow.

`configs/project_paths.toml` provides the human-readable path vocabulary.
Python workflows discover the repository from marker files or the optional
`TNK_PROJECT_ROOT` environment variable, so commands do not depend on the
current working directory.

## Code

Active entry points are grouped under `workflows/`:

- `preintegration/`: source conversion, gene schema, metadata, and TCR joins
- `integration/`: canonical Phases 0 through 4 and extended-atlas construction
- `intake/`: external and supplementary dataset processing
- `metadata/`: tissue and disease harmonization
- `analysis/`: focused downstream analyses
- `gdtai/`: model training, inference, comparison, and failure audits
- `gdt_atlas/`: gdT atlas construction and phenotype analyses
- `reporting/`: figure and report builders
- `maintenance/`: registries, checkpoints, migration, and verification

Shared dependency-light utilities live in `src/tnk_atlas/`. Scientific entry
points remain directly executable to avoid forcing a large package refactor.

## Data and outputs

`data/datasets/<dataset_id>/` is the physical pre-integration home for each
study. `data/shared/` stores project-wide pre-integration artifacts,
`data/compat/` preserves legacy path trees, and `data/raw`, `data/interim`, and
`data/processed` provide lifecycle-oriented views. `Integrated_dataset/`
remains the canonical output root required by the runbook. Large milestone
H5ADs can remain on the approved SSD mirror through
`high_speed_temp/Integrated_dataset`.

`reports/` describes curated reader-facing packages. Existing
`gdT_prediction/` and `gdT_atlas/` paths remain stable for the port 8000 static
server.

## Archive

`archive/` contains retired code and historical records. Active workflows must
not import from it. Every reorganization move is listed in a CSV plan and a
result map under `maintenance/reorganization/`.
