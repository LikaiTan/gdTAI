# CHANGELOG.md

## Completed milestones

- Initial control-file framework established
- Phase 0 dataset audit completed and repaired registry rebuilt
- Phase 1 coarse extraction completed on the approved registry inputs
- Phase 1b conservative cleanup completed
- Phase 1c merged metadata backup and replacement completed
- Supplementary 10x 5' intake completed through supplementary Phase 1
- Supplementary candidates merged into the main candidate milestone after user approval
- Phase 2 merged cleanup completed
- Phase 3 scVI integration completed
- scANVI outputs were retained as reference-only after user review
- Phase 4 TRAB/TRB/TRD scoring completed
- tissue correction completed and written back into the integrated milestone
- pre-merge raw per-GSE TCR audit completed
- repaired TCR rebuild and propagation completed for ten approved GSEs
- non-destructive project reorganization completed on 2026-07-15
  - active scripts grouped by workflow responsibility
  - raw/interim/processed data lifecycle represented by validated symlinks
  - dataset, library, file, and gdTAI model registries established
  - legacy code archived with checksum-guarded move maps
  - rollback, reproducibility, dataset-change, and model-iteration guides added
  - all milestone H5AD sampled hashes verified unchanged
- gdTAI v3 Round 12 versus Round 14 revalidation completed on 2026-07-15
  - both artifacts pinned and preserved by SHA256
  - identical full-atlas, atlas-held-out, and independent-external evaluation
  - Round 14 selected as the balanced default by prespecified guardrails and
    mean F1
  - Round 12 retained as the validated high-purity fallback
- dataset-centered physical storage migration completed on 2026-07-16
  - 191 same-filesystem rename operations completed without copying or
    rewriting H5AD files
  - 65 datasets organized under `data/datasets/<dataset_id>/`
  - 44 selected H5ADs promoted through `processed/current.h5ad`
  - 262 lifecycle links and all historical top-level paths retained
  - all 1,200 file-registry rows resolve to identical canonical/legacy inodes
  - rollback plan, journals, registry snapshot, and validation evidence stored
    under `data/registry/migrations/dataset_centered_20260716/`
- independent `BALF_BLOOD_COPD` validation cohort intake completed on
  2026-07-16
  - complete 951 GB study workspace moved intact with one same-filesystem
    rename into `data/datasets/BALF_BLOOD_COPD/workspace/`
  - historical `/home/tanlikai/databank/owndata/singlecell` path retained as a
    compatibility symlink to the same directory inode
  - 46,273-cell validation H5AD exposed through
    `data/datasets/BALF_BLOOD_COPD/processed/current.h5ad`
  - four 10x 5' BALF/PBMC libraries registered as validation-only and inactive
    for atlas integration
  - 12 selected raw, interim, TCR, demultiplexing, and H5AD sentinels validated
    without failures
- `BALF_BLOOD_COPD` storage refined on 2026-07-16 at user request
  - raw, FASTQ, BAM, matrix, TCR, demultiplexing, and other interim files
    restored to the original `/home/tanlikai/databank/owndata/singlecell`
    workspace with its directory inode preserved
  - only the 2,110,825,599-byte validation H5AD remains physically under
    `data/datasets/BALF_BLOOD_COPD/processed/artifacts/`
  - the H5AD was relocated by a same-filesystem rename with its inode,
    timestamp, dimensions, and sampled SHA256 preserved
  - project raw access and the historical external H5AD path remain valid
    through compatibility links
  - strict registries, all 264 lifecycle links, 12 sentinels, and active
    workflow path resolution validated without failures
- GSE305372 external CD8 gdTAI application completed on 2026-07-30
  - checksum-verified lung and lung-associated lymph-node CD8 Seurat objects
    stored under the dataset-centered raw path
  - 157,846 cells selected only by the authors' `CD8A` CITE-seq tag
  - promoted gdTAI v3 Round 14 model applied with exact raw-count CP10K
    normalization and fixed threshold `0.936`
  - 651 cells called positive; donor, author-cluster, TCR-gene quadrant, and
    paired TRA/TRB conflict-screening summaries generated
  - cohort registered as inactive for training, Phase 0, and atlas integration
- GSE305372 all-T V2/V3 correction completed on 2026-07-30
  - downloaded and checksummed the two CD4 Seurat objects and audited the
    metadata-only DICE-TCR blood CD4 resource as ineligible for expression
    inference
  - removed the erroneous secondary CITE-CD8A inclusion filter and analyzed
    all 690,715 cells in the four transcriptome-eligible CD4/CD8 objects
  - applied V2 high-F1, V2 high-purity, and promoted V3 Round 14 to an identical
    raw-count CP10K input; called 2,474, 1,942, and 2,043 cells, respectively
  - reconciled the 189-cell lung CD4 RDS/manuscript difference to unmapped
    author cluster IDs 6-8 and documented the partial lymph-node deposit
  - marked the earlier 157,846-cell CD8A package as superseded provenance

## Current working era

- post-Phase 4 downstream analysis and documentation cleanup
