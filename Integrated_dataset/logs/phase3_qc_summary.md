# Phase 3 QC Summary

## Scope
- Input milestone: `Integrated_dataset/TNK_cleaned.h5ad`
- Canonical output milestone target: `Integrated_dataset/integrated.h5ad`
- Current validated high-speed output path: `/ssd/tnk_phase3/Integrated_dataset/integrated.h5ad`
- Execution env: `rapids_sc_py310`
- RAPIDS import note: `torch` was imported before `rapids_singlecell` to avoid the CUDA symbol-resolution issue in this env
- Figures were written in PNG format only

## Integration summary
- Cells integrated: 6,443,879
- Genes retained in the integrated milestone: 27,413
- Leiden clusters: 29
- HVGs used for scVI: 4,000
- HVG method: reused_saved_model
- HVG exclusions before clustering/UMAP: 3,935 total (mitochondrial=13, ribosomal=214, noncoding_like=3,708)
- Phase 3 preflight cells removed for numerical safety: 5,947
- Unique Phase 3 batches: 1,104
- Unique GSEs represented: 37

## scANVI annotation summary
- scANVI subset cells: 236,476
- scANVI subset strata: 977
- Label transfer method: `subset_scANVI_plus_latent_centroid`
- Detailed labels observed: 229
- Coarse superclasses observed: 3
- `reference_other` cells: 2,278,375
- Mean scANVI confidence: 0.1034

## Metadata attachment
- Metadata-matched cells: 6,443,879
- Metadata-unmatched cells: 0

## Top detailed labels
- effector memory CD8-positive, alpha-beta T cell: 263,114 cells (superclass=T_cell, mean_confidence=0.1225)
- CD16-positive, CD56-dim natural killer cell, human: 241,688 cells (superclass=NK_cell, mean_confidence=0.0902)
- naive thymus-derived CD4-positive, alpha-beta T cell: 208,760 cells (superclass=T_cell, mean_confidence=0.1068)
- CD8-positive, alpha-beta cytotoxic T cell: 206,003 cells (superclass=T_cell, mean_confidence=0.1226)
- naive thymus-derived CD8-positive, alpha-beta T cell: 179,860 cells (superclass=T_cell, mean_confidence=0.2030)
- CD4-positive helper T cell: 177,922 cells (superclass=T_cell, mean_confidence=0.0818)
- natural killer cell: 143,202 cells (superclass=NK_cell, mean_confidence=0.1483)
- T-helper 1 cell: 139,797 cells (superclass=reference_other, mean_confidence=0.0862)
- effector memory CD4-positive, alpha-beta T cell: 132,245 cells (superclass=T_cell, mean_confidence=0.0778)
- CD8-positive, alpha-beta T cell: 129,593 cells (superclass=T_cell, mean_confidence=0.1660)

## Outputs
- Summary table: `Integrated_dataset/tables/phase3_integration_summary.csv`
- Metadata join summary: `Integrated_dataset/tables/phase3_metadata_join_summary.csv`
- Batch summary: `Integrated_dataset/tables/phase3_batch_key_summary.csv`
- Leiden summary: `Integrated_dataset/tables/phase3_leiden_cluster_summary.csv`
- Label summary: `Integrated_dataset/tables/phase3_scanvi_label_summary.csv`
- Marker agreement summary: `Integrated_dataset/tables/phase3_marker_agreement_summary.csv`
- HVG filter summary: `Integrated_dataset/tables/phase3_hvg_filter_summary.csv`
- scANVI subset summary: `Integrated_dataset/tables/phase3_scanvi_subset_summary.csv`
- Cluster label summary: `Integrated_dataset/tables/phase3_cluster_label_summary.csv`
- Input sanitization table: `Integrated_dataset/tables/phase3_input_sanitization.csv`
- Figures: `Integrated_dataset/figures/phase3_umap_by_gse.png`, `Integrated_dataset/figures/phase3_umap_by_batch_level.png`, `Integrated_dataset/figures/phase3_umap_by_leiden.png`, `Integrated_dataset/figures/phase3_leiden_cluster_sizes.png`, `Integrated_dataset/figures/phase3_umap_by_scanvi_detailed_label.png`, `Integrated_dataset/figures/phase3_umap_by_scanvi_superclass.png`, `Integrated_dataset/figures/phase3_scanvi_confidence_distribution.png`, `Integrated_dataset/figures/phase3_marker_agreement.png`

## QC conclusion
- The integrated milestone is accepted without rollback and currently carries the validated scVI latent space, Leiden clustering, and UMAP results.
- The user reviewed the scANVI output and judged it too messy for primary downstream use.
- The scANVI label fields, tables, and PNGs are retained in the object and on NFS for reference only.
- Downstream interpretation should use simple scVI-based annotation together with Leiden structure and marker review unless a later user decision changes that.
