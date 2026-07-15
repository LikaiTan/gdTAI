# Dataset registry

This directory is the source of truth for data membership before integration.

- `datasets.csv` records source-level status and integration roles.
- `libraries.csv` records library-level RNA/TCR scope. Placeholder rows must be
  curated before the next integration rebuild.
- `files.csv` records source, interim, and processed file roles and migration
  targets.

Regenerate the initial inventory with:

```bash
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/maintenance/build_dataset_registry.py
```

Do not delete a dataset to remove it from an integration. Set its active status
to false, record the reason, validate the registry, and create a locked registry
snapshot for the replacement integration run.
