# Data snapshots

`storage_index.csv` records canonical dataset directories and promoted H5ADs.
`compatibility_links.json` records the generated lifecycle view, and
`migrations/` contains physical storage plans and journals.

Before any new integration, freeze `datasets.csv`, `libraries.csv`, and
`files.csv` into `snapshots/<run_id>/` together with checksums and the Git
commit. The integrated H5AD must record the snapshot ID in `adata.uns`.
