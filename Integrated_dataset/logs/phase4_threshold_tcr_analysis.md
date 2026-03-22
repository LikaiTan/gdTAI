# Phase 4 Threshold and TCR Pairing Analysis

- Integrated milestone: `/ssd/tnk_phase3/Integrated_dataset/integrated.h5ad`
- Raw-score criterion 1: `phase4_trab_score - phase4_trd_score < -0.6`
- Raw-score criterion 2: `phase4_trd_score > 0.1`
- TCR criterion: both `TRA_cdr3` and `TRB_cdr3` present in harmonized metadata

## Overall summary
- `TRAB - TRD < -0.6 (raw)`: total_positive_cells=36161, n_gse_with_positive_cells=25/37, mean_positive_cells_per_gse=977.32, median_positive_cells_per_gse=155.00, max_positive_cells_in_one_gse=8933
- `TRD > 0.1 (raw)`: total_positive_cells=370527, n_gse_with_positive_cells=25/37, mean_positive_cells_per_gse=10014.24, median_positive_cells_per_gse=3506.00, max_positive_cells_in_one_gse=95291
- `Paired TRA/TRB CDR3 present`: total_positive_cells=4573710, n_gse_with_positive_cells=28/37, mean_positive_cells_per_gse=123613.78, median_positive_cells_per_gse=55480.00, max_positive_cells_in_one_gse=807365

## Top GSE contributions
- `TRAB - TRD < -0.6 (raw)` top 5 GSE by positive cells:
  - `GSE254249`: positive_cells=8933, within_gse_fraction=0.0220, share_of_all_positive=0.2470
  - `GSE243013`: positive_cells=5897, within_gse_fraction=0.0078, share_of_all_positive=0.1631
  - `GSE228597`: positive_cells=3378, within_gse_fraction=0.0225, share_of_all_positive=0.0934
  - `GSE212217`: positive_cells=3219, within_gse_fraction=0.0151, share_of_all_positive=0.0890
  - `GSE243905`: positive_cells=2898, within_gse_fraction=0.0083, share_of_all_positive=0.0801
- `TRD > 0.1 (raw)` top 5 GSE by positive cells:
  - `GSE254249`: positive_cells=95291, within_gse_fraction=0.2344, share_of_all_positive=0.2572
  - `GSE243013`: positive_cells=64299, within_gse_fraction=0.0847, share_of_all_positive=0.1735
  - `GSE212217`: positive_cells=27790, within_gse_fraction=0.1300, share_of_all_positive=0.0750
  - `GSE243905`: positive_cells=27189, within_gse_fraction=0.0776, share_of_all_positive=0.0734
  - `GSE228597`: positive_cells=25762, within_gse_fraction=0.1720, share_of_all_positive=0.0695
- `Paired TRA/TRB CDR3 present` top 5 GSE by positive cells:
  - `GSE161918`: positive_cells=807365, within_gse_fraction=1.0000, share_of_all_positive=0.1765
  - `GSE243013`: positive_cells=759336, within_gse_fraction=0.9999, share_of_all_positive=0.1660
  - `GSE168859`: positive_cells=544306, within_gse_fraction=1.0000, share_of_all_positive=0.1190
  - `GSE254249`: positive_cells=406494, within_gse_fraction=1.0000, share_of_all_positive=0.0889
  - `GSE243905`: positive_cells=350438, within_gse_fraction=1.0000, share_of_all_positive=0.0766
