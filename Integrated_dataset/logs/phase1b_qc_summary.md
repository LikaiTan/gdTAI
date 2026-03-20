# Phase 1b QC Summary

## Scope
- Input milestone: `Integrated_dataset/TNK_candidates.h5ad`
- Phase 1b goal: conservative first-pass cleanup only
- Gene filter added per user instruction: remove genes expressed in fewer than 500 cells

## Cell cleanup
- Cells before: 5,951,196
- Cells after: 5,950,935
- Total cells removed: 261
- Removed as obvious non-T/NK: 261
- Removed as high-confidence doublets: 0
- GSEs with any removed cells: GSE168859, GSE168163, GSE125527

## Gene cleanup
- Genes before: 57,562
- Genes after: 23,536
- Genes removed for detection < 500: 34,026
- Key T/NK genes kept: CD3E, CD3D, CD3G, NKG7, TRBC2, KLRD1, TRBC1, TRAC, TRDC, TRGC2, TRGC1
- Key T/NK genes dropped: none

## Outputs
- Summary table: `Integrated_dataset/tables/phase1b_cleanup_summary.csv`
- GSE before/after table: `Integrated_dataset/tables/phase1b_gse_before_after.csv`
- Removed cell table: `Integrated_dataset/tables/phase1b_removed_cells.csv`
- Gene detection table: `Integrated_dataset/tables/phase1b_gene_detection_summary.csv`
- Figures: `Integrated_dataset/figures/phase1b_gse_cell_retention.png`, `Integrated_dataset/figures/phase1b_cell_removal_reasons.png`, `Integrated_dataset/figures/phase1b_gene_detection_distribution.png`

## QC conclusion
- Phase 1b remains conservative: only high-confidence obvious contaminants are removed at the cell level.
- No broad T/NK-like population is removed by this rule set.
- This milestone is ready for user QC review before any Phase 1c / Phase 2 work.
