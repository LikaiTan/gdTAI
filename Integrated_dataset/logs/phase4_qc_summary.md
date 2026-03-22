# Phase 4 QC Summary

## Scope
- Input and updated output milestone: `/ssd/tnk_phase3/Integrated_dataset/integrated.h5ad`
- Package source: `/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research/gdt_tcr_module_sharing_package_full.zip`
- Scoring mode: exact package TRA/TRB/TRD modules on a temporary normalized/log1p copy of count-space `X`
- Canonical result type: continuous scores only; no hard gdT/abT call was written in Phase 4
- Derived scaled outputs: min-max scaled `phase4_trd_score_scaled` and `phase4_trab_score_scaled` in the 0-1 range, plus `phase4_trd_minus_trab_scaled` as their difference
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
- `phase4_trd_score_scaled`: mean=0.0836, median=0.0738, p95=0.1809, max=1.0000
- `phase4_trab_score_scaled`: mean=0.2595, median=0.2571, p95=0.3737, max=1.0000
- `phase4_trd_minus_trab_scaled`: mean=-0.1760, median=-0.1781, p95=-0.0541, max=0.7696

## Top Leiden clusters by median TRD - TRAB
- Leiden `6`: raw_median=0.0736, scaled_median=-0.1078, n_cells=87536
- Leiden `7`: raw_median=0.0721, scaled_median=-0.0794, n_cells=265681
- Leiden `23`: raw_median=0.0683, scaled_median=-0.1104, n_cells=175055
- Leiden `31`: raw_median=0.0662, scaled_median=-0.0930, n_cells=1782
- Leiden `27`: raw_median=0.0644, scaled_median=-0.1138, n_cells=7815

## Top GSEs by median TRD - TRAB
- `GSE144469`: raw_median=0.0690, scaled_median=-0.1142, n_cells=113064
- `GSE241783`: raw_median=0.0686, scaled_median=-0.0979, n_cells=600641
- `GSE161918`: raw_median=0.0643, scaled_median=-0.1356, n_cells=807365
- `GSE232240`: raw_median=0.0633, scaled_median=-0.1305, n_cells=25359
- `GSE179994`: raw_median=0.0610, scaled_median=-0.1370, n_cells=148785

## Outputs
- Tables: `phase4_score_summary.csv`, `phase4_module_gene_membership.csv`, `phase4_leiden_score_summary.csv`, `phase4_gse_score_summary.csv`, `phase4_top_cells_by_trd_minus_trab.csv`
- Figures: `phase4_score_distributions.png`, `phase4_umap_score_overlays.png`, `phase4_leiden_score_summary.png`, `phase4_gse_score_summary.png`, `phase4_marker_score_comparison.png`, `phase4_trab_vs_trd_scatter_panel.png`, `phase4_trab_vs_trd_tcr_presence.png`

## TRA/TRB CDR3-presence overlay
- Additional scatter panel: `phase4_trab_vs_trd_tcr_presence.png`
- Axes: raw `TRAB` vs `TRD` on the left, scaled `TRAB` vs `TRD` on the right
- Coloring rule: cells with any `TRA_cdr3` or `TRB_cdr3` sequence are blue; cells without either are red
- Join source: harmonized metadata tables joined by `project name`, `sampleid`, and `barcodes`
- Sampled-cell counts in the rendered panel: `31,160` blue and `8,840` red

## Additional Phase 4 extension outputs
- Additional figures: `phase4_trab_vs_trd_paired_tratrb_vs_no_tcr.png`, `phase4_trab_vs_trd_trgc_trgv9_expression.png`
- Additional tables: `phase4_paired_tratrb_vs_no_tcr_sample_counts.csv`, `phase4_threshold_tcr_summary.csv`, `phase4_trab_minus_trd_lt_neg0p6_by_gse.csv`, `phase4_trd_gt_0p1_by_gse.csv`, `phase4_paired_tratrb_by_gse.csv`
- Additional analysis note: `phase4_threshold_tcr_analysis.md`

## Paired TRA/TRB vs no-TCR sampled scatter
- Scatter sample size: `20,000`
- `Paired TRA/TRB`: `14,321`
- `No TCR`: `4,263`
- `Other TCR state` (single-chain or mixed non-paired state): `1,416`

## Threshold summary
- `TRAB - TRD < -0.6 (raw)`: `36,161` cells total, mean `977.32` per GSE across `37` GSEs, detected in `25` GSEs
- `TRD > 0.1 (raw)`: `370,527` cells total, mean `10,014.24` per GSE across `37` GSEs, detected in `25` GSEs
- `Paired TRA/TRB CDR3 present`: `4,573,710` cells total, mean `123,613.78` per GSE across `37` GSEs, detected in `28` GSEs

## Largest GSE contributions
- `TRAB - TRD < -0.6 (raw)`: `GSE254249` contributes `8,933` cells (`24.70%` of all positives)
- `TRD > 0.1 (raw)`: `GSE254249` contributes `95,291` cells (`25.72%` of all positives)
- `Paired TRA/TRB CDR3 present`: `GSE161918` contributes `807,365` cells (`17.65%` of all positives)

## QC conclusion
- Phase 4 wrote continuous TRA/TRB/TRAB/TRD module scores and `TRD - TRAB` back into the canonical integrated milestone.
- Phase 4 also wrote min-max scaled `TRD` and `TRAB` scores in the 0-1 range plus the scaled `TRD - TRAB` difference.
- The shared package was used as the method source of truth for module definitions, but its ground-truth evaluation branch was intentionally not used here.
- Phase 4 results should be interpreted jointly with Leiden structure, marker genes, and the scVI embedding.
