# GSE305372 CD8 gdTAI Application

## Scope

GSE305372 is an external application cohort from the study *Human lungs
maintain tissue-resident memory T cells against a broad spectrum of pathogens*.
It is not part of gdTAI training, model selection, Phase 0, or atlas
integration.

The analysis uses the study's two processed CD8 Seurat objects:

| Compartment | Source cells | Author CD8A cells analyzed | File size | SHA256 |
| --- | ---: | ---: | ---: | --- |
| Lung | 227,800 | 143,184 | 5,623,543,662 bytes | `b4e1f62c4f34fc6b024f3769dbffd0e9513cf7d22bab915e0dd1a1e3eb4fca71` |
| Lung-associated lymph node | 18,072 | 14,662 | 515,670,151 bytes | `2e36d7ff4ab771e400a13cb807431198aa8ec8f7b7d1554842cb49a1196958e7` |

Only these processed CD8 objects were downloaded. Raw FASTQ files and the CD4
processed objects were not downloaded because they are not required for this
inference task.

## Cell Selection And Model Input

Cells were selected independently of gdTAI using the authors' annotation:

```text
cite.cell.type.tag == CD8A
```

The exporter reads the Seurat `RNA` assay's raw `counts` layer, verifies all 210
packaged gdTAI genes, and computes each cell's library size over all 36,601 RNA
features. For both objects, those sums exactly matched `nCount_RNA` (maximum
absolute difference zero).

To avoid duplicating the multi-gigabyte full matrix, the exporter writes only
the 210 model-gene counts plus whole-transcriptome library sizes. Inference then
uses the model's exact normalization:

```text
log1p(raw gene count * 10,000 / whole-transcriptome raw-count total)
```

The checksum-pinned promoted model is gdTAI v3 Round 14:

- artifact: `Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/gdTAI_v3_model.pkl`
- model SHA256: `16dedc0081da9b8487887341232bcf8c9c9403dd3bbd72e04cab43d4cd7b2e09`
- fixed threshold: `0.936`

## Result

The model called 651 of 157,846 author-tagged CD8 cells (0.412%):

| Compartment | Cells | gdTAI positive | Fraction |
| --- | ---: | ---: | ---: |
| Lung | 143,184 | 608 | 0.425% |
| Lung-associated lymph node | 14,662 | 43 | 0.293% |
| Total | 157,846 | 651 | 0.412% |

Among the 651 calls, 579 expressed at least one TRDV gene in the model input.
The highest author-cluster call rates were lymph-node MAIT (4/165, 2.42%) and
lung NKG2C-positive cells (89/5,474, 1.63%), so these shared cytotoxic or
innate-like states remain important review strata.

The processed objects include TRA/TRB clonotype metadata but no TRG/TRD
clonotype fields. Among 128,446 cells with paired TRA and TRB evidence, 245
were gdTAI positive (0.191%). This is reported as an alpha-beta
conflict-screening rate, not a definitive false-positive rate: matched
gamma-delta TCR evidence is unavailable, and rare dual-receptor biology cannot
be excluded.

This is an inference analysis, not an external precision/recall evaluation.
The CD8A-only scope also excludes CD8-negative gamma-delta T cells.

## Reproduction

```bash
workflows/intake/download_gse305372_cd8.sh

Rscript workflows/intake/export_gse305372_cd8_model_payload.R \
  --input-rds data/datasets/GSE305372/raw/source/GSE305372_HIPC-1_LG-CD8-ALL_SeuratObject.RDS \
  --output-prefix data/datasets/GSE305372/interim/rds_export/LG_CD8 \
  --compartment lung \
  --feature-manifest Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/feature_genes.csv

Rscript workflows/intake/export_gse305372_cd8_model_payload.R \
  --input-rds data/datasets/GSE305372/raw/source/GSE305372_HIPC-1_LN-CD8-ALL_SeuratObject.RDS \
  --output-prefix data/datasets/GSE305372/interim/rds_export/LN_CD8 \
  --compartment lymph_node \
  --feature-manifest Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0/feature_genes.csv

/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  workflows/gdtai/run_gdtai_gse305372_cd8.py --overwrite
```

## Outputs

- report: `reports/GSE305372_gdtai_cd8/index.html`
- predictions: `Integrated_dataset/tables/gdT_prediction/GSE305372/`
- figures: `Integrated_dataset/figures/gdT_prediction/GSE305372/`
- manifest and summary: `Integrated_dataset/logs/gdT_prediction/GSE305372/`
- local source and feature payloads: `data/datasets/GSE305372/`

Source records:

- GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305372>
- article: <https://www.nature.com/articles/s41590-026-02592-6>
- authors' analysis repository: <https://github.com/vijaybioinfo/HIPC_Lung01>
