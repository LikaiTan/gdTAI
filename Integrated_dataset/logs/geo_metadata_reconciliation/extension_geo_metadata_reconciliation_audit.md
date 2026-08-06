# Extension GEO Metadata Reconciliation Audit

- Result: **PASS**
- Mode: read-only H5AD metadata inspection
- Reconciliation rows: 58
- Cohorts: 8
- Official GEO accessions: 10
- Accession coverage rows passing: 10/10
- H5AD files unchanged during audit: 8/8

## Status counts

- `ambiguous_geo_partial`: 6
- `resolved_geo_value_local_unresolved`: 32
- `resolved_verified_local`: 16
- `unavailable_in_geo`: 4

## Cohort summary

| cohort_id | status | reconciliation_rows |
|---|---|---:|
| GSE114724 | ambiguous_geo_partial | 1 |
| GSE114724 | resolved_geo_value_local_unresolved | 2 |
| GSE114724 | resolved_verified_local | 2 |
| GSE121636_GSE121637 | ambiguous_geo_partial | 1 |
| GSE121636_GSE121637 | resolved_geo_value_local_unresolved | 5 |
| GSE121636_GSE121637 | resolved_verified_local | 2 |
| GSE159251 | resolved_geo_value_local_unresolved | 4 |
| GSE159251 | resolved_verified_local | 3 |
| GSE169246 | ambiguous_geo_partial | 1 |
| GSE169246 | resolved_geo_value_local_unresolved | 8 |
| GSE169246 | resolved_verified_local | 1 |
| GSE292700 | ambiguous_geo_partial | 1 |
| GSE292700 | resolved_geo_value_local_unresolved | 2 |
| GSE292700 | resolved_verified_local | 2 |
| GSE292700 | unavailable_in_geo | 1 |
| GSE294273_GSE294274 | resolved_geo_value_local_unresolved | 4 |
| GSE294273_GSE294274 | resolved_verified_local | 3 |
| GSE294273_GSE294274 | unavailable_in_geo | 1 |
| GSE296954 | ambiguous_geo_partial | 1 |
| GSE296954 | resolved_geo_value_local_unresolved | 4 |
| GSE296954 | resolved_verified_local | 1 |
| GSE296954 | unavailable_in_geo | 1 |
| GSE315928 | ambiguous_geo_partial | 1 |
| GSE315928 | resolved_geo_value_local_unresolved | 3 |
| GSE315928 | resolved_verified_local | 2 |
| GSE315928 | unavailable_in_geo | 1 |

## Outputs

- `Integrated_dataset/tables/geo_metadata_reconciliation/validated_reconciliation.csv`
- `Integrated_dataset/tables/geo_metadata_reconciliation/reconciliation_status_summary.csv`
- `Integrated_dataset/tables/geo_metadata_reconciliation/source_accession_coverage.csv`
- `Integrated_dataset/tables/geo_metadata_reconciliation/local_field_inventory.csv`
- `Integrated_dataset/tables/geo_metadata_reconciliation/scoped_local_value_counts.csv`
- `Integrated_dataset/tables/geo_metadata_reconciliation/unresolved_or_ambiguous.csv`
- `Integrated_dataset/tables/geo_metadata_reconciliation/h5ad_read_only_validation.csv`
