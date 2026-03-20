# Phase 1c Metadata Replacement

## Scope
- Input candidate milestone: `Integrated_dataset/TNK_candidates.h5ad`
- Replaced metadata target: `analysis_26GSE_V4/outputs/harmonized_metadata_v4.csv`
- Backup file: `analysis_26GSE_V4/outputs/metadata.csv.bk`
- Join executed by `source_gse_id + original_cell_id` to `gse_id + cell_id` with explicit string typing.
- Required output join columns were written as `project name`, `sampleid`, and `barcodes`.

## Validation
- Candidate rows: 5,950,935
- Replacement rows: 5,950,935
- Previous metadata rows: 7,035,204
- Rows removed vs previous metadata: 1,084,269
- Duplicate replacement join keys: 0
- Blank `project name` after write: 0
- Blank `barcodes` after write: 0
- Blank `sampleid` after write: 58,678

## Outputs
- Obs export: `Integrated_dataset/tables/phase1c_merged_obs_export.csv.gz`
- Summary table: `Integrated_dataset/tables/phase1c_metadata_replacement_summary.csv`
- GSE join summary: `Integrated_dataset/tables/phase1c_metadata_join_by_gse.csv`

## Conclusion
- Phase 1c completed successfully.
- The harmonized metadata target now corresponds exactly to the current filtered merged candidate object.
- Phase 2 must still wait for explicit user approval.
