# Script Migration

Root-level and data-embedded scripts were classified into active workflow
groups or the archive. The source and destination of every move are recorded in
`maintenance/reorganization/*_migration_plan.csv` and the corresponding map.

Active Python scripts received a small bootstrap that discovers the project
root and exposes shared workflow modules. This preserves historical direct
imports and allows entry points to run from outside the repository.

The bootstrap tool supports an exact checksum-verified rollback for the main
root-script migration. Archive and configuration moves use the same guarded
migration utility. Historical files are not rewritten after archiving.

Use `workflows/maintenance/verify_reorganization.py` to check that every plan
is fully applied, all compatibility links resolve, formal model hashes match,
and milestone H5AD baselines remain unchanged.
