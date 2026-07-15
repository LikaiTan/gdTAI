# Data snapshots

`compatibility_links.json` records the generated non-destructive data view.
Before any new integration, freeze `datasets.csv`, `libraries.csv`, and
`files.csv` into `snapshots/<run_id>/` together with checksums and the Git commit.
The integrated H5AD must record the snapshot ID in `adata.uns`.
