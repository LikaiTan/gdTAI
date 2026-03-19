# Phase 1 QC Summary

## Scope

- Input set: Category A datasets only
- Datasets processed: 13
- Total input cells: 2946858
- Total retained candidate cells: 2689420
- Overall retained fraction: 0.9126

## Selection Signals

- Cells retained by annotation support: 1775304
- Cells retained by marker support: 2638573
- Datasets with annotation columns used: 9

## Per-dataset yield

- GSE243013: retained 766574 / 766574 (1.0000); annotation columns=`orig.ident;major_cell_type;sub_cell_type`
- GSE241783: retained 618963 / 643762 (0.9615); annotation columns=`none`
- GSE168859: retained 564270 / 564270 (1.0000); annotation columns=`cell_type`
- GSE228597: retained 149994 / 317473 (0.4725); annotation columns=`cluster_name;lineage`
- GSE144469: retained 142069 / 142069 (1.0000); annotation columns=`cell_type`
- GSE221776: retained 126178 / 126332 (0.9988); annotation columns=`none`
- GSE168163: retained 80998 / 80998 (1.0000); annotation columns=`cell_type`
- GSE229858: retained 65483 / 65757 (0.9958); annotation columns=`none`
- GSE155223: retained 55694 / 80789 (0.6894); annotation columns=`cell_type`
- GSE125527: retained 47649 / 68627 (0.6943); annotation columns=`celltype;cell_type`
- GSE155222: retained 39138 / 52071 (0.7516); annotation columns=`cell_type`
- GSE232240: retained 26707 / 32399 (0.8243); annotation columns=`cell_type`
- GSE240361: retained 5703 / 5737 (0.9941); annotation columns=`none`

## QC Conclusion

Phase 1 Category A coarse extraction is complete. Review TNK_candidates.h5ad, the per-dataset selection table, and the Phase 1 figures before deciding whether to proceed to Phase 1b or adjust thresholds.
