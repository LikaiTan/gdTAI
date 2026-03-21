# Phase 4 QC Summary

## Scope
- Input and updated output milestone: `/ssd/tnk_phase3/Integrated_dataset/integrated.h5ad`
- Package source: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/gdt_tcr_module_sharing_package_full.zip`
- Scoring mode: exact package TRA/TRB/TRD modules on a temporary normalized/log1p copy of count-space `X`
- Canonical result type: continuous scores only; no hard gdT/abT call was written in Phase 4
- scANVI labels remain reference-only and were not used to define Phase 4 outputs

## Module sizes
- `tra` genes: 93; control genes used: 897
- `trb` genes: 68; control genes used: 897
- `trab` genes: 161; control genes used: 1044
- `trd` genes: 7; control genes used: 300

## Score summary
- `phase4_tra_score`: mean=-0.0113, median=-0.0162, p95=0.0332, max=0.4811
- `phase4_trb_score`: mean=-0.0557, median=-0.0579, p95=0.0165, max=0.3352
- `phase4_trab_score`: mean=-0.0564, median=-0.0575, p95=-0.0056, max=0.2734
- `phase4_trd_score`: mean=0.0002, median=-0.0155, p95=0.1560, max=1.4674
- `phase4_trd_minus_trab`: mean=0.0566, median=0.0432, p95=0.2157, max=1.5275

## Top Leiden clusters by median TRD - TRAB
- Leiden `6`: median=0.0736, n_cells=87536
- Leiden `7`: median=0.0721, n_cells=265681
- Leiden `23`: median=0.0683, n_cells=175055
- Leiden `31`: median=0.0662, n_cells=1782
- Leiden `27`: median=0.0644, n_cells=7815

## Top GSEs by median TRD - TRAB
- `GSE144469`: median=0.0690, n_cells=113064
- `GSE241783`: median=0.0686, n_cells=600641
- `GSE161918`: median=0.0643, n_cells=807365
- `GSE232240`: median=0.0633, n_cells=25359
- `GSE179994`: median=0.0610, n_cells=148785

## Outputs
- Tables: `phase4_score_summary.csv`, `phase4_module_gene_membership.csv`, `phase4_leiden_score_summary.csv`, `phase4_gse_score_summary.csv`, `phase4_top_cells_by_trd_minus_trab.csv`
- Figures: `phase4_score_distributions.png`, `phase4_umap_score_overlays.png`, `phase4_leiden_score_summary.png`, `phase4_gse_score_summary.png`, `phase4_marker_score_comparison.png`

## QC conclusion
- Phase 4 wrote continuous TRA/TRB/TRAB/TRD module scores and `TRD - TRAB` back into the canonical integrated milestone.
- The shared package was used as the method source of truth for module definitions, but its ground-truth evaluation branch was intentionally not used here.
- Phase 4 results should be interpreted jointly with Leiden structure, marker genes, and the scVI embedding.
