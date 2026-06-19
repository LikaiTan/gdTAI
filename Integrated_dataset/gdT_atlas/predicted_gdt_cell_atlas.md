# Predicted gdT Cell Atlas H5AD

Output H5AD: `Integrated_dataset/gdT_atlas/predicted_gdt_cell_atlas.h5ad`

This H5AD is the handoff object for downstream gdT biology analysis. It is a subset of the full plus6 integrated atlas and contains cells selected by the user-promoted gdTAI v3.0 model, with primary-gold gdT false negatives added back and known/likely false positives removed.

## Source

- Source H5AD: `high_speed_temp/Integrated_dataset/integrated_plus6.h5ad`
- Prediction table: `Integrated_dataset/tables/gdT_prediction/gdtai_v3_trdc_nk_guard/full_atlas_selected_predicted_cells.csv.gz`
- Promoted model directory: `Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0`
- Promoted model: `gdTAI_v3_model.pkl`
- Selected candidate: `v3_round14_v2_score_trdc_gate_fixed_0p936`
- Threshold: `0.936`

## Inclusion Rule

1. Start with all full-atlas cells predicted positive by the selected gdTAI v3 candidate.
2. Remove predicted positives that are known/likely false positives:
   - primary abT-gold predicted positives (`class_code == 1`), and
   - paired TCRAB / any abTCR cells that are not primary gdT gold.
3. Add back primary-gold gdT cells that the model missed (`class_code == 2` and not predicted).

Silver-only gdT-like cells are not added back unless they were already predicted positive.

## Cell Counts

- Full input cells: `5,128,904`
- Raw gdTAI v3 predicted positives: `251,356`
- Known/likely FP removed: `5,580`
- Primary abT-gold FP removed: `1,183`
- Paired TCRAB/no-gdT FP removed: `5,580`
- Primary-gold gdT FN added back: `11,232`
- Final atlas cells: `257,008`
- Final predicted-retained cells: `245,776`
- Final gold-FN-added cells: `11,232`

## Added Obs Columns

- `gdtai_v3_predicted`: whether the promoted model predicted the cell.
- `gdtai_v3_score`: promoted model probability score; `NaN` for gold FN add-back cells.
- `gdtai_v3_threshold`: operating threshold for predicted cells.
- `gdtai_v3_tcr_gene_quadrant`: TRDC/TRDV expression quadrant when available from prediction output.
- `gdtai_v3_primary_class`: `unlabeled_or_ambiguous`, `abT_gold`, `gdT_gold`, or `gdT_silver`.
- `gdtai_v3_gold_fn_added`: true for primary-gold gdT cells added despite missed prediction.
- `gdtai_v3_atlas_inclusion_reason`: `gdtai_v3_predicted_retained` or `gold_gdT_FN_added`.

## Audit Tables

- Summary: `Integrated_dataset/tables/gdT_atlas/predicted_gdt_cell_atlas_summary.csv`
- Final cell metadata: `Integrated_dataset/tables/gdT_atlas/predicted_gdt_cell_atlas_metadata.csv.gz`
- Removed FP cells: `Integrated_dataset/tables/gdT_atlas/predicted_gdt_cell_atlas_removed_fp_cells.csv.gz`
- Gold FN add-back cells: `Integrated_dataset/tables/gdT_atlas/predicted_gdt_cell_atlas_gold_fn_added_cells.csv.gz`
- By source: `Integrated_dataset/tables/gdT_atlas/predicted_gdt_cell_atlas_by_source.csv`
- By tissue: `Integrated_dataset/tables/gdT_atlas/predicted_gdt_cell_atlas_by_tissue.csv`

## Previous Atlas Outputs Removed

- `Integrated_dataset/gdT_atlas` (directory, 8,166,175,300 bytes)
- `Integrated_dataset/tables/gdT_atlas` (directory, 2,780,289 bytes)
- `Integrated_dataset/tables/gdT_atlas_curated_phenotypes` (directory, 39,699,652 bytes)
- `Integrated_dataset/tables/gdT_atlas_rigorous` (directory, 42,611,093 bytes)
- `Integrated_dataset/figures/gdT_atlas` (directory, 2,442,713 bytes)
- `Integrated_dataset/figures/gdT_atlas_curated_phenotypes` (directory, 13,097,894 bytes)
- `Integrated_dataset/figures/gdT_atlas_rigorous` (directory, 25,742,999 bytes)
- `Integrated_dataset/logs/gdT_atlas` (directory, 3,772 bytes)
- `Integrated_dataset/logs/gdT_atlas_curated_phenotypes` (directory, 5,957 bytes)
- `Integrated_dataset/logs/gdT_atlas_rigorous` (directory, 3,600 bytes)
- `Integrated_dataset/models/gdT_atlas_rigorous` (directory, 10,794,525 bytes)
- `gdT_atlas` (directory, 38,702,299 bytes)

## Downstream Use

Use this H5AD as the starting object for downstream gdT-only analysis. The object has not been re-integrated or reclustered after subsetting; downstream agents should run their own gdT-only integration, clustering, annotation, and biological analyses from this file.
