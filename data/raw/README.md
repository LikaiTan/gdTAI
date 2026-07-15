# Raw data

Raw data are immutable. Preserve source filenames, accession information,
download provenance, and checksums. Code, figures, logs, Cell Ranger outputs,
and analysis H5ADs do not belong here.

The generated `geo/<dataset_id>/` and `local/<dataset_id>/` directories are
compatibility views. A directory named `fastq_legacy_mixed` explicitly warns
that the legacy target still contains both source and derived files and needs a
dataset-level physical migration before it can be declared immutable.
