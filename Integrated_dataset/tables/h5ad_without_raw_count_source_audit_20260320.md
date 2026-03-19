# H5AD Without Raw Count Source Audit

- Input list: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/tables/h5ad_without_raw_count.csv`
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- Datasets audited: 20
- Datasets with no internal raw/counts in current H5AD: 20
- Datasets with `selected_inputs.csv` plus `per_sample_load_stats.csv`: 20
- Datasets with supplementary raw tarball: 18
- Datasets with supplementary RDS files: 6

## Conclusion

Most of these H5ADs are not truly missing raw sources. The processed H5ADs lack internal counts, but the project scaffolds still retain external source paths that can be used to rebuild count matrices. RDS should be treated as fallback only when matrix sources are unavailable or incomplete.

## Per-GSE Summary

### GSE139555
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `unknown`
- Example input path(s): `downloads/GSE139555/suppl/extracted_final_GSE139555/GSM4143655_SAM24345862-lt1.matrix.mtx.gz | downloads/GSE139555/suppl/extracted_final_GSE139555/GSM4143656_SAM24345863-ln1.matrix.mtx.gz | downloads/GSE139555/suppl/extracted_final_GSE139555/GSM4143657_SAM24348188-lt2.matrix.mtx.gz`
- Likely cause: processed_h5ad lost counts during conversion or post-load normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment
- Notes: RDS present, but matrix source should be preferred over RDS for count rescue ; matches user example of RDS-to-H5AD conversion likely dropping raw counts

### GSE145926
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `filtered`
- Example input path(s): `downloads/GSE145926/suppl/extracted_final_GSE145926/GSM4339769_C141_filtered_feature_bc_matrix.h5 | downloads/GSE145926/suppl/extracted_final_GSE145926/GSM4339770_C142_filtered_feature_bc_matrix.h5 | downloads/GSE145926/suppl/extracted_final_GSE145926/GSM4339771_C143_filtered_feature_bc_matrix.h5`
- Likely cause: processed_h5ad lost counts after loading/normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment

### GSE161918
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `filtered`
- Example input path(s): `downloads/GSE161918/suppl/extracted_final_GSE161918/GSM4929081_B1_10xlane1_RNA_filtered_feature_bc_matrix.h5 | downloads/GSE161918/suppl/extracted_final_GSE161918/GSM4929082_B1_10xlane2_RNA_filtered_feature_bc_matrix.h5 | downloads/GSE161918/suppl/extracted_final_GSE161918/GSM4929083_B1_10xlane3_RNA_filtered_feature_bc_matrix.h5`
- Likely cause: processed_h5ad lost counts during conversion or post-load normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment
- Notes: RDS present, but matrix source should be preferred over RDS for count rescue

### GSE162498
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `unknown`
- Example input path(s): `downloads/GSE162498/suppl/unpacked_scanpy/GSM4952953_P34_Tumor_raw_feature_bc_matrix/P34_Tumor_raw_feature_bc_matrix/matrix.mtx | downloads/GSE162498/suppl/unpacked_scanpy/GSM4952954_P35_Tumor_raw_feature_bc_matrix/P35_Tumor_raw_feature_bc_matrix/matrix.mtx | downloads/GSE162498/suppl/unpacked_scanpy/GSM4952955_P42_Tumor_raw_feature_bc_matrix/P42_Tumor_raw_feature_bc_matrix/matrix.mtx`
- Likely cause: processed_h5ad lost counts during conversion or post-load normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment
- Notes: raw 10x matrix likely includes empty droplets; large source n_obs is not by itself an error ; RDS present, but matrix source should be preferred over RDS for count rescue

### GSE171037
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `filtered`
- Example input path(s): `downloads/GSE171037/suppl/extracted_final_GSE171037/GSM5216975_Sample_20063a009_01_filtered_feature_bc_matrix.h5 | downloads/GSE171037/suppl/extracted_final_GSE171037/GSM5216976_Sample_20063a010_01_filtered_feature_bc_matrix.h5 | downloads/GSE171037/suppl/extracted_final_GSE171037/GSM5216977_Sample_20063a011_01_filtered_feature_bc_matrix.h5`
- Likely cause: processed_h5ad lost counts after loading/normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment

### GSE178882
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `filtered`
- Example input path(s): `downloads/GSE178882/suppl/extracted_final_GSE178882/GSM5399838_c14_filtered_feature_bc_matrix.h5 | downloads/GSE178882/suppl/extracted_final_GSE178882/GSM5399839_c13_filtered_feature_bc_matrix.h5 | downloads/GSE178882/suppl/extracted_final_GSE178882/GSM5399840_c16_filtered_feature_bc_matrix.h5`
- Likely cause: processed_h5ad lost counts after loading/normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment

### GSE188620
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `unknown`
- Example input path(s): `downloads/GSE188620/suppl/extracted_final_GSE188620/GSM5686902_103_018_screening_matrix.mtx.gz | downloads/GSE188620/suppl/extracted_final_GSE188620/GSM5686903_103_018_surgery_matrix.mtx.gz | downloads/GSE188620/suppl/extracted_final_GSE188620/GSM5686904_103_026_screening_matrix.mtx.gz`
- Likely cause: processed_h5ad lost counts after loading/normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment

### GSE190870
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `raw`
- Example input path(s): `downloads/GSE190870/suppl/extracted_final_GSE190870/GSM5732357_A_raw_feature_bc_matrix.h5 | downloads/GSE190870/suppl/extracted_final_GSE190870/GSM5732358_B_raw_feature_bc_matrix.h5 | downloads/GSE190870/suppl/extracted_final_GSE190870/GSM5732359_C_raw_feature_bc_matrix.h5`
- Likely cause: processed_h5ad lost counts after loading/normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment

### GSE211504
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `unknown`
- Example input path(s): `downloads/GSE211504/suppl/extracted_final_GSE211504/GSM6474817_matrix_1.1.mtx.gz | downloads/GSE211504/suppl/extracted_final_GSE211504/GSM6474818_matrix_1.2.mtx.gz | downloads/GSE211504/suppl/extracted_final_GSE211504/GSM6474819_matrix_1.3.mtx.gz`
- Likely cause: processed_h5ad lost counts after loading/normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment

### GSE212217
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `filtered`
- Example input path(s): `downloads/GSE212217/suppl/extracted_final_GSE212217/GSM6514096_PEM1C1_scRNA_filtered_feature_bc_matrix.h5 | downloads/GSE212217/suppl/extracted_final_GSE212217/GSM6514097_PEM1C3_scRNA_filtered_feature_bc_matrix.h5 | downloads/GSE212217/suppl/extracted_final_GSE212217/GSM6514098_PEM2C1_scRNA_filtered_feature_bc_matrix.h5`
- Likely cause: processed_h5ad lost counts during conversion or post-load normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment
- Notes: RDS present, but matrix source should be preferred over RDS for count rescue

### GSE227709
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `unknown`
- Example input path(s): `downloads/GSE227709/suppl/extracted_final_GSE227709/GSM7105884_Z1471A_matrix.mtx.gz | downloads/GSE227709/suppl/extracted_final_GSE227709/GSM7105885_Z1471B_matrix.mtx.gz | downloads/GSE227709/suppl/extracted_final_GSE227709/GSM7105886_Z1471C_matrix.mtx.gz`
- Likely cause: processed_h5ad lost counts after loading/normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment
- Notes: raw 10x matrix likely includes empty droplets; large source n_obs is not by itself an error

### GSE243572
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `raw`
- Example input path(s): `downloads/GSE243572/suppl/extracted_final_GSE243572/GSM7791217_SPFS_C1D1_1_5_1_10_1_12_1_17_GEX_raw_feature_bc_matrix.h5 | downloads/GSE243572/suppl/extracted_final_GSE243572/GSM7791219_SPFS_C1D1_1_21_1_23_1_26_3_2_GEX_raw_feature_bc_matrix.h5 | downloads/GSE243572/suppl/extracted_final_GSE243572/GSM7791221_SPFS_C1D1_1_2_1_27_3_1_3_5_GEX_raw_feature_bc_matrix.h5`
- Likely cause: processed_h5ad lost counts after loading/normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment

### GSE243905
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `filtered`
- Example input path(s): `downloads/GSE243905/suppl/extracted_final_GSE243905/GSM7798067_C4_1_1_BLD_filtered_feature_bc_matrix.h5 | downloads/GSE243905/suppl/extracted_final_GSE243905/GSM7798069_C4_1_1_CSF_filtered_feature_bc_matrix.h5 | downloads/GSE243905/suppl/extracted_final_GSE243905/GSM7798071_C4_2_1_BLD_filtered_feature_bc_matrix.h5`
- Likely cause: processed_h5ad lost counts during conversion or post-load normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment
- Notes: RDS present, but matrix source should be preferred over RDS for count rescue

### GSE252762
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `unknown`
- Example input path(s): `downloads/GSE252762/suppl/GSE252762_batch1_rna_matrix.mtx.gz | downloads/GSE252762/suppl/GSE252762_batch2_rna_matrix.mtx.gz | downloads/GSE252762/suppl/GSE252762_batch3_rna_matrix.mtx.gz`
- Likely cause: processed_h5ad lost counts after loading/normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment
- Notes: RDS present, but matrix source should be preferred over RDS for count rescue

### GSE254176
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `unknown`
- Example input path(s): `downloads/GSE254176/suppl/extracted_final_GSE254176/GSM8035466_OM-HL-007-C8.matrix.mtx.gz | downloads/GSE254176/suppl/extracted_final_GSE254176/GSM8035467_OM-HS-066-C6.matrix.mtx.gz | downloads/GSE254176/suppl/extracted_final_GSE254176/GSM8035468_OM-HS-066-2-C6.matrix.mtx.gz`
- Likely cause: processed_h5ad lost counts after loading/normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment

### GSE254249
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `unknown`
- Example input path(s): `downloads/GSE254249/suppl/GSE254249_scRNA_matrix.mtx.gz`
- Likely cause: processed_h5ad lost counts after loading/normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment
- Notes: single combined matrix source; rescue should operate at dataset level, not per-sample

### GSE267645
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `unknown`
- Example input path(s): `downloads/GSE267645/suppl/extracted_final_GSE267645/GSM8271506_pt1_moderate_GEX_matrix.mtx.gz | downloads/GSE267645/suppl/extracted_final_GSE267645/GSM8271507_pt1_critical_GEX_matrix.mtx.gz | downloads/GSE267645/suppl/extracted_final_GSE267645/GSM8271508_pt2_moderate_GEX_matrix.mtx.gz`
- Likely cause: processed_h5ad lost counts after loading/normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment

### GSE301528
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `unknown`
- Example input path(s): `downloads/GSE301528/suppl/extracted_final_GSE301528/GSM9084836_bfa1_matrix.mtx.gz | downloads/GSE301528/suppl/extracted_final_GSE301528/GSM9084837_bfa2_matrix.mtx.gz | downloads/GSE301528/suppl/extracted_final_GSE301528/GSM9084838_bfc1_matrix.mtx.gz`
- Likely cause: processed_h5ad lost counts after loading/normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment

### GSE308075
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `unknown`
- Example input path(s): `downloads/GSE308075/suppl/extracted_final_GSE308075/GSM9237523_SCH_1_matrix.mtx.gz | downloads/GSE308075/suppl/extracted_final_GSE308075/GSM9237524_SCH_2_matrix.mtx.gz | downloads/GSE308075/suppl/extracted_final_GSE308075/GSM9237525_SCH_5_matrix.mtx.gz`
- Likely cause: processed_h5ad lost counts after loading/normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment

### GSE311112
- Current H5AD internal raw/counts: no
- Backup directory: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/Integrated_dataset/backups/h5ad_without_raw_count_20260320_032651`
- External manifest available: True
- Load-stats available: True
- Matrix kind(s): `unknown`
- Example input path(s): `downloads/GSE311112/suppl/extracted_correct_GSE311112/GSM9317153_pt_1_baseline_matrix.mtx.gz | downloads/GSE311112/suppl/extracted_correct_GSE311112/GSM9317156_pt_1_3yrs_matrix.mtx.gz | downloads/GSE311112/suppl/extracted_correct_GSE311112/GSM9317159_pt_1_5yrs_matrix.mtx.gz`
- Likely cause: processed_h5ad lost counts after loading/normalization; external raw files exist
- Recommended solution: rebuild counts from selected_inputs.csv and replace adata.X after validating obs/var alignment

