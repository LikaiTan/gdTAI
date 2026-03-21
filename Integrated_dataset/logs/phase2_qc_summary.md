# Phase 2 QC Summary

## Scope
- Input milestone: `Integrated_dataset/TNK_candidates.h5ad`
- Output milestone: `Integrated_dataset/TNK_cleaned.h5ad`
- Phase 2 goal: merged-context second-pass cleanup without aggressive loss of plausible T/NK states
- Gene filter: keep genes expressed in at least 500 cleaned cells

## Cell cleanup
- Cells before: 6,848,556
- Cells after: 6,449,826
- Total cells removed: 398,730
- Removal fraction: 5.8221%
- Removed for low quality: 13,082
- Removed for extreme carried Phase 1 contamination: 0
- Removed for myeloid conflict: 114,557
- Removed for B-cell conflict: 240,950
- Removed for epithelial conflict: 16,266
- Removed for erythroid conflict: 9,318
- Removed for platelet conflict: 4,557
- GSEs with any removed cells: GSE161918, GSE243013, GSE241783, GSE168859, GSE254249, GSE243905, GSE287301, GSE235863, GSE212217, GSE243572, GSE139555, GSE179994, GSE228597, GSE144469, GSE188620, GSE221776, GSE311112, GSE308075, GSE287541, GSE254176, GSE267645, GSE162498, GSE168163, GSE178882, GSE229858, GSE155223, GSE301528, GSE125527, GSE252762, GSE155222, GSE171037, GSE211504, GSE145926, GSE232240, GSE240865, GSE234069, GSE190870

## Gene cleanup
- Genes before: 55,455
- Genes after: 27,413
- Genes removed for detection < 500: 28,042
- Key T/NK genes kept: CD3E, CD3D, CD3G, NKG7, TRBC2, KLRD1, TRBC1, TRAC, TRDC, TRGC2, TRGC1
- Key T/NK genes dropped: none

## Outputs
- Summary table: `Integrated_dataset/tables/phase2_cleanup_summary.csv`
- GSE before/after table: `Integrated_dataset/tables/phase2_gse_before_after.csv`
- Removed cells table: `Integrated_dataset/tables/phase2_removed_cells.csv`
- Gene detection table: `Integrated_dataset/tables/phase2_gene_detection_summary.csv`
- QC sample table: `Integrated_dataset/tables/phase2_sampled_cells_for_qc.csv`
- Figures: `Integrated_dataset/figures/phase2_gse_cell_retention.png`, `Integrated_dataset/figures/phase2_removal_reasons.png`, `Integrated_dataset/figures/phase2_offtarget_vs_tnk.png`, `Integrated_dataset/figures/phase2_gene_detection_distribution.png`

## QC conclusion
- This Phase 2 run remains marker-conflict driven and avoids broad cluster-level pruning.
- If off-target removals stay small and key T/NK genes are preserved, the cleaned milestone is acceptable for Phase 3 integration.
