# Phase 1 QC Summary

## Scope

- Input set: Category A datasets only
- Datasets processed: 29
- Total input cells: 6182881
- Total retained candidate cells: 5322388
- Overall retained fraction: 0.8608

## Selection Signals

- Cells retained by annotation support: 2014589
- Cells retained by marker support: 5265923
- Datasets with annotation columns used: 15

## Per-dataset yield

- GSE161918: retained 901093 / 1100189 (0.8190); annotation columns=`orig.ident;seurat_clusters;BatchClusters;mergedcelltype;coursecelltype;adjustedcelltype;WCTmergedcelltype;WCTcoursecelltype;celltypeQC;sample_label`
- GSE243013: retained 766574 / 766574 (1.0000); annotation columns=`orig.ident;major_cell_type;sub_cell_type`
- GSE241783: retained 618963 / 643762 (0.9615); annotation columns=`none`
- GSE168859: retained 564270 / 564270 (1.0000); annotation columns=`cell_type`
- GSE243905: retained 369428 / 413696 (0.8930); annotation columns=`orig.ident;predicted.celltype.l2.score;predicted.celltype.l2;consolidated_cell_type`
- GSE212217: retained 220612 / 292271 (0.7548); annotation columns=`orig.ident;seurat_clusters;finalIdent`
- GSE243572: retained 207017 / 268664 (0.7705); annotation columns=`none`
- GSE139555: retained 163970 / 195824 (0.8373); annotation columns=`ident`
- GSE228597: retained 149994 / 317473 (0.4725); annotation columns=`cluster_name;lineage`
- GSE144469: retained 142069 / 142069 (1.0000); annotation columns=`cell_type`
- GSE188620: retained 130325 / 151027 (0.8629); annotation columns=`none`
- GSE221776: retained 126178 / 126332 (0.9988); annotation columns=`none`
- GSE311112: retained 114945 / 139112 (0.8263); annotation columns=`none`
- GSE308075: retained 112645 / 113504 (0.9924); annotation columns=`none`
- GSE254176: retained 102776 / 118332 (0.8685); annotation columns=`none`
- GSE162498: retained 82285 / 91658 (0.8977); annotation columns=`none`
- GSE168163: retained 80998 / 80998 (1.0000); annotation columns=`cell_type`
- GSE178882: retained 74242 / 109138 (0.6803); annotation columns=`none`
- GSE229858: retained 65483 / 65757 (0.9958); annotation columns=`none`
- GSE155223: retained 55694 / 80789 (0.6894); annotation columns=`cell_type`
- GSE125527: retained 47649 / 68627 (0.6943); annotation columns=`celltype;cell_type`
- GSE252762: retained 45202 / 45209 (0.9998); annotation columns=`celltype;clusters;clusters_full`
- GSE155222: retained 39138 / 52071 (0.7516); annotation columns=`cell_type`
- GSE171037: retained 37970 / 60981 (0.6227); annotation columns=`published_meta_celltype;published_celltype;harmonized_celltype;T_cell_type`
- GSE211504: retained 37432 / 37637 (0.9946); annotation columns=`none`
- GSE145926: retained 33026 / 84069 (0.3928); annotation columns=`none`
- GSE232240: retained 26707 / 32399 (0.8243); annotation columns=`cell_type`
- GSE240361: retained 5703 / 5737 (0.9941); annotation columns=`none`
- GSE227709: retained 0 / 14712 (0.0000); annotation columns=`none`

## QC Conclusion

Phase 1 Category A coarse extraction is complete. Review TNK_candidates.h5ad, the per-dataset selection table, and the Phase 1 figures before deciding whether to proceed to Phase 1b or adjust thresholds.
