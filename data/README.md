# Data lifecycle

The `data/` tree is the organized view of all inputs before atlas integration.
It separates immutable source material, reproducible intermediate products, and
canonical per-dataset H5ADs.

- `raw/`: immutable downloaded or user-supplied source data
- `interim/`: extracted archives, Cell Ranger output, conversions, and TCR joins
- `processed/per_dataset/`: one selected standardized H5AD per dataset
- `registry/`: compatibility links and immutable input snapshots

During migration these locations are symlink views over legacy paths. The
legacy paths remain valid until dataset-by-dataset checksum validation and an
explicit cutover. No integration workflow should write into `raw/`.
