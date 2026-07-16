# Data lifecycle

The `data/` tree is the canonical home for all inputs before atlas integration.
The primary organization is dataset-centered so all material for one study can
be found under `datasets/<dataset_id>/`.

- `datasets/`: physical per-dataset raw, interim, and processed storage
- `shared/`: project-wide references and legacy integration products
- `compat/`: old path trees retained for existing scripts and notebooks
- `raw/`, `interim/`, `processed/`: lifecycle-oriented compatibility views
- `registry/`: compatibility links and immutable input snapshots

The top-level `downloads`, `analysis_26GSE_V4`, and `newdata` names remain as
compatibility aliases after physical cutover. New code should resolve data
through `configs/datasets/datasets.csv` and `tnk_atlas.paths.ProjectPaths`.
