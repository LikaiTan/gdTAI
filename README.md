# T/NK Atlas and gdTAI

This repository contains the reproducible T/NK single-cell atlas workflow,
Phase 4 TCR-gene scoring, gdTAI model training and inference, and downstream
gdT atlas analyses. The current scientific milestone is post-Phase 4 review;
the reorganization does not modify any milestone H5AD.

## Start here

Read these files in order before running analysis:

1. `TNK_PIPELINE_RUNBOOK.md`
2. `TNK_PIPELINE_STATUS.md`
3. `TNK_PHASES_0_4_SCRIPT.md`

Use the canonical Python environment:

```bash
PY=/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python
$PY workflows/maintenance/manage_dataset_registry.py validate --strict-libraries
$PY workflows/maintenance/build_data_compatibility_view.py --verify
$PY workflows/maintenance/verify_reorganization.py
$PY -m unittest discover -s tests -v
```

## Layout

- `src/tnk_atlas/`: shared path, provenance, registry, and model-loading code
- `workflows/`: executable workflows grouped by lifecycle and scientific task
- `configs/`: dataset, metadata, TCR, atlas, and model configuration
- `data/`: raw, interim, and processed compatibility view before integration
- `Integrated_dataset/`: canonical scientific outputs and model artifacts
- `reports/`: curated report entry points; existing served URLs remain stable
- `archive/`: superseded code and historical records, never active by default
- `maintenance/reorganization/`: checksum-guarded move plans and maps
- `tests/`: structural and reproducibility checks

## Operating guides

- `docs/PROJECT_STRUCTURE.md`
- `docs/DATA_LIFECYCLE.md`
- `docs/REPRODUCIBILITY.md`
- `docs/ROLLBACK.md`
- `docs/ADDING_REMOVING_DATASETS.md`
- `docs/ITERATING_GDTAI.md`
- `docs/SCRIPT_MIGRATION.md`

The source-data directories under `downloads/`, `newdata/`, and the legacy
Scanpy project tree remain in place. The organized `data/` tree uses symlinks,
so current paths keep working and a rollback does not require moving large
scientific files.
