# Integrated T/NK Dataset for Gamma-Delta T Cell Model Training

This integrated atlas is not presented as a finished biological truth set. Its practical purpose is to support the next step: training a model for gamma-delta T cell identification across many studies, tissues, and technical contexts. Phase 3 provides the integrated manifold. Phase 4 provides a prior scoring method that is useful but imperfect, so it is used as a weak-label generator for training-set construction.

## Report Metrics

- Integrated cells: `6,443,879`
- Technology: `10x genomics 5' sequencing`
- Represented GSEs: `37`
- Leiden clusters: `29`
- Tissue groups: `15`
- `high confidence gdT`: `36,161`
- `γδT-like (TRDscore > 0.1)`: `370,527`


## Interpretation Strategy

The key design choice is to separate integration from training-label generation. Phase 3 solves the first problem by aligning a very large multi-study T/NK dataset into one manifold while retaining interpretable lineage structure. Phase 4 solves the second problem by applying a previously developed TRA/TRB/TRD scoring method that captures useful gamma-delta signal, even though it is not perfect enough to be treated as final truth.

The intended downstream use is weak supervision: high-confidence Phase 4 regions can seed positive and negative training examples, while ambiguous regions remain available for later model-based refinement.

## Phase 3: Integration Makes the Dataset Trainable

The Phase 3 plots show that cells from many studies and technical batches can be embedded in one stable latent space. This matters because the final model should generalize across studies rather than memorize one dataset or one tissue.

### Phase 3 UMAP by GSE

The integrated manifold mixes cells from many studies into one shared latent space. This is the key prerequisite for building a cross-study training resource rather than a single-dataset classifier.

### Phase 3 UMAP by Tissue 
"write some thing herer"

### Phase 3 UMAP by Leiden Cluster

Leiden clustering provides a structured view of the integrated manifold. These clusters help assess whether Phase 4 scores are coherently localized rather than randomly distributed.

### Phase 3 T/NK Marker Panel

Core T-cell, gamma-delta-associated, and NK-associated markers remain interpretable after integration. This supports the use of the integrated object as a training substrate rather than an over-corrected latent space with erased lineage structure.

## Phase 4: Imperfect but Useful Weak Labels

The Phase 4 method was created previously for the gamma-delta T cell identification problem. It works, but not perfectly. In this project it is therefore used in the correct role: not as final annotation, but as a weak-labeling layer for building a large and diverse training set.
"write the formular or codes about how to calculate TRD/TRAB score"
At the start of Phase 4, scores are calculated from the count-space integrated matrix using the package-faithful TRA, TRB, and TRD gene modules. For each cell, the workflow creates a temporary normalize-total plus log1p view, computes the mean expression of the module genes, computes the mean expression of matched control genes, and defines the module score as module mean minus control mean. The TRAB score uses the union of the TRA and TRB modules, while the TRD score uses the TRD module alone. The raw difference used downstream is therefore TRD score minus TRAB score.

For this report, visualization stays in raw-score space. A strong cutoff such as `TRAB - TRD < -0.6` gives a selective gamma-delta-like candidate subset that can seed high-confidence training labels, while the rest of the raw score landscape can be used for harder negatives or ambiguous examples.

### Phase 4 Raw Score UMAP Overlays

The report uses raw Phase 4 scores for visualization. Raw TRAB, raw TRD, and raw TRD-minus-TRAB remain spatially structured on the integrated manifold, which is the key property needed for weak-label generation.

### Raw TRAB-versus-TRD Score Space

"Change the color of this plot to TRD-TRAB"
The raw score space is sufficient to separate highly gamma-delta-like regions from alpha-beta-oriented or ambiguous regions. The scaled version is therefore not needed for the narrative in this report.

### Raw TRAB-versus-TRD Space with Paired TRA/TRB Context
"YOU USED A WRONG FIGURE HERE!!! REGENERATE THIS FIGURE PLEASE！"
Cells with paired alpha-beta TRA/TRB evidence occupy a recognizable portion of the raw score space. This provides an internal consistency check for negative or non-gamma-delta-oriented training examples.

### Selected Marker Genes on Raw TRAB-versus-TRD Space

The raw score space can be inspected directly against a broader marker panel, including gamma-delta-associated, alpha-beta-associated, pan-T-cell, and cytotoxic-lineage genes. This makes the weak-label logic more transparent for downstream model training.



### High-Confidence TRD-over-TRAB Candidates by Tissue `36,161`

Cells with very strong gamma-delta-like score contrast are concentrated in specific tissues rather than uniformly spread. This suggests that the eventual model should learn across tissue contexts rather than rely on one source alone.

### High-Confidence TRD-over-TRAB Candidates by GSE `36,161`

The strongest raw-score gamma-delta-like candidates are distributed across multiple studies rather than coming from a single dataset, which is exactly the behavior needed for a cross-study training resource.

### Broad TRD-Enriched Candidates by Tissue`370,527`

The broader TRD-enriched candidate pool covers multiple tissue contexts. This wider pool can support less stringent, higher-recall training subsets or later rounds of iterative label refinement.

### Broad TRD-Enriched Candidates by GSE`370,527`

The broader TRD-enriched pool also spans many studies, which is useful when training a model that must generalize beyond a narrow set of discovery datasets.

## How This Dataset Should Be Used

The integrated object is large enough to support model training across multiple tissues, multiple GSEs, and multiple confidence levels. A practical training design is to use the strongest Phase 4 regions as seed labels, use clearly alpha-beta-oriented regions as negatives, and keep the middle score space as ambiguous data for later iterative improvement.

This is the main value of the dataset: it turns a heterogeneous collection of public T/NK datasets into one training-ready resource where imperfect prior biology can be converted into scalable model supervision.
