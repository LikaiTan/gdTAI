# GSE305372 All-T gdTAI V2/V3 Application

## Status

This is the canonical GSE305372 gdTAI application. It supersedes the earlier
CD8A CITE-tag-filtered analysis documented in `GSE305372_GDTAI_CD8.md`.

GSE305372 remains an external application cohort. It was not used for gdTAI
training, threshold selection, Phase 0, or atlas integration.

## Source Resources And Eligibility

| GEO resource | Size | SHA256 | gdTAI eligibility |
| --- | ---: | --- | --- |
| Lung CD4 Seurat RDS | 8,368,100,613 | `5b5ff9b8e5b63bd3b0f0a0e09a1d639cc0c4e5e809c238cb998c58cf13100f9e` | eligible |
| Lung CD8 Seurat RDS | 5,623,543,662 | `b4e1f62c4f34fc6b024f3769dbffd0e9513cf7d22bab915e0dd1a1e3eb4fca71` | eligible |
| Lymph-node CD4 Seurat RDS | 2,066,160,300 | `59817faf26be5568bef616cfe4addd3f25d8fd357aa5ee8a3d1950f0de181068` | eligible |
| Lymph-node CD8 Seurat RDS | 515,670,151 | `2e36d7ff4ab771e400a13cb807431198aa8ec8f7b7d1554842cb49a1196958e7` | eligible |
| DICE-TCR blood CD4 CSV | 76,152,700 | `50c6f946ea67c4e2194a926bc121a802f2dea9ee12b879151eff2f0740c8ad44` | ineligible |

The DICE-TCR CSV has 1,143,667 data rows and 35 TCR/metadata columns, but no
raw transcriptome matrix or whole-library count denominator. It cannot be
scored by an expression model and is excluded from the inference denominator.

## Corrected Cell Scope

Every cell in the four transcriptome-eligible processed Seurat objects is
included. `cite.cell.type.tag` is retained only for auditing and is never an
inclusion criterion.

| Object | Deposited cells analyzed | Manuscript cells | Scope note |
| --- | ---: | ---: | --- |
| Lung CD4 | 370,422 | 370,233 | 189 deposited cells are in unmapped cluster IDs 6-8 |
| Lung CD8 | 227,800 | 227,800 | complete manuscript object |
| Lymph-node CD4 | 74,421 | 198,890 | partial GEO deposit |
| Lymph-node CD8 | 18,072 | 129,737 | partial GEO deposit |
| **Total analyzed** | **690,715** | **926,660** | all transcriptome-eligible deposited cells |

The lung CD4 mapped clusters 0-5 contain exactly 370,233 cells. The additional
189 RDS cells account for the complete difference from the manuscript count;
none is positive under V2 high-F1, V2 high-purity, or V3 Round 14.

## Model Input And Operating Modes

The exporter reads the Seurat `RNA` assay raw `counts` layer. It exports the
union 210-gene V3 panel, which contains every one of the 187 V2 genes, while
retaining the whole-transcriptome library size for each cell.

Both models use:

```text
log1p(raw gene count * 10,000 / whole-transcriptome raw-count total)
```

Operating results:

- V2 high-F1: fixed threshold `0.9063586592674255`
- V2 high-purity: CD4 `0.97`, CD8 `0.93`, author-defined Treg disabled
- V3 Round 14: fixed threshold `0.936`

Published author cluster mappings are taken from the study repository:

- lung: `RNA_snn_res.0.2`
- lymph node: `RNA_snn_res.0.3`
- lung CD4 cluster 3: Treg
- lymph-node CD4 clusters 3 and 4: Treg

The released model checksums are verified against
`configs/models/gdtai/model_registry.csv` before inference.

## Results

| Model | Predicted | Fraction |
| --- | ---: | ---: |
| V2 high-F1 | 2,474 | 0.3582% |
| V2 high-purity | 1,942 | 0.2812% |
| V3 Round 14 | 2,043 | 0.2958% |

Per-object calls:

| Object | V2 high-F1 | V2 high-purity | V3 Round 14 |
| --- | ---: | ---: | ---: |
| Lung CD4 | 378 | 136 | 269 |
| Lung CD8 | 1,934 | 1,695 | 1,640 |
| Lymph-node CD4 | 51 | 14 | 40 |
| Lymph-node CD8 | 111 | 97 | 94 |

V3 and V2 high-purity share 1,884 predictions. V3 has 159 additional calls,
V2 high-purity has 58 additional calls, and their Jaccard index is `0.8967`.
Both call sets are strict subsets of V2 high-F1 in this cohort.

The raw-expression audit identifies:

- 3,482 CD3+TRDC+ cells
- 14,387 CD3+(TRDC or TRDV)+ cells
- V2 high-F1 captures 594 of 3,482 CD3+TRDC+ cells
- V2 high-purity captures 517 of 3,482 CD3+TRDC+ cells
- V3 captures 533 of 3,482 CD3+TRDC+ cells

These capture fractions are dropout-sensitive expression-proxy summaries, not
ground-truth recall.

Among 533,411 cells with paired TRA/TRB evidence, the model-positive
conflict-screening rates are:

- V2 high-F1: 923 (`0.1730%`)
- V2 high-purity: 687 (`0.1288%`)
- V3: 733 (`0.1374%`)

The objects do not contain matched TRG/TRD clonotype fields, so paired-TRA/TRB
model calls are not automatically definitive false positives.

## Reproduction

```bash
workflows/intake/download_gse305372_all_t.sh

Rscript workflows/intake/export_gse305372_all_t_model_payload.R \
  --input-rds data/datasets/GSE305372/raw/source/GSE305372_HIPC-1_LG-CD4-ALL_SeuratObject.RDS \
  --output-prefix data/datasets/GSE305372/interim/all_t_model_payload/LG_CD4 \
  --compartment lung --lineage CD4 \
  --feature-manifest Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/feature_genes.csv

Rscript workflows/intake/export_gse305372_all_t_model_payload.R \
  --input-rds data/datasets/GSE305372/raw/source/GSE305372_HIPC-1_LG-CD8-ALL_SeuratObject.RDS \
  --output-prefix data/datasets/GSE305372/interim/all_t_model_payload/LG_CD8 \
  --compartment lung --lineage CD8 \
  --feature-manifest Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/feature_genes.csv

Rscript workflows/intake/export_gse305372_all_t_model_payload.R \
  --input-rds data/datasets/GSE305372/raw/source/GSE305372_HIPC-1_LN-CD4-ALL_SeuratObject.RDS \
  --output-prefix data/datasets/GSE305372/interim/all_t_model_payload/LN_CD4 \
  --compartment lymph_node --lineage CD4 \
  --feature-manifest Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/feature_genes.csv

Rscript workflows/intake/export_gse305372_all_t_model_payload.R \
  --input-rds data/datasets/GSE305372/raw/source/GSE305372_HIPC-1_LN-CD8-ALL_SeuratObject.RDS \
  --output-prefix data/datasets/GSE305372/interim/all_t_model_payload/LN_CD8 \
  --compartment lymph_node --lineage CD8 \
  --feature-manifest Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/feature_genes.csv

/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/run_gdtai_gse305372_all_t.py --overwrite
```

## Canonical Outputs

- HTML: `reports/GSE305372_gdtai_all_t/index.html`
- PDF: `reports/GSE305372_gdtai_all_t/GSE305372_gdTAI_v2_v3_all_T_report.pdf`
- tables: `Integrated_dataset/tables/gdT_prediction/GSE305372_all_T/`
- figures: `Integrated_dataset/figures/gdT_prediction/GSE305372_all_T/`
- manifest and validation: `Integrated_dataset/logs/gdT_prediction/GSE305372_all_T/`
- feature payloads: `data/datasets/GSE305372/interim/all_t_model_payload/`

Source records:

- GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305372>
- article: <https://www.nature.com/articles/s41590-026-02592-6>
- authors' analysis repository: <https://github.com/vijaybioinfo/HIPC_Lung01>
