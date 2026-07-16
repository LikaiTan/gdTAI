# Raw data

Raw data are immutable. Preserve source filenames, accession information,
download provenance, and checksums. Code, figures, logs, Cell Ranger outputs,
and analysis H5ADs do not belong here.

The generated `geo/<dataset_id>/` and `local/<dataset_id>/` directories are
compatibility views over `data/datasets/<dataset_id>/raw/`. A directory named
`fastq_legacy_mixed` explicitly warns that the migrated source bundle still
contains both source and derived files; curate it within the canonical dataset
directory rather than moving it back to a legacy root.
