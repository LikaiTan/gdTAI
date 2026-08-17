# CHANGELOG.md

## 2026-08-18 - Rebuilt flagged TCR joins with per-chain UMI provenance

- Rebuilt productive TRA, TRB, TRG, and TRD calls for 14 flagged sources from
  raw VDJ records using only `sample_id + barcode_core`.
- Retained per-chain UMI/read support and availability, selected deterministic
  highest-support calls, and kept unavailable support null rather than zero.
- Passed 11 sources fully or partially and quarantined GSE125527, GSE228597,
  and GSE287541; blanked 110 duplicate-key GSE235863 rows from receptor truth.
- Staged 2,479,137 replacement rows with 948,991 validated productive-TCR
  cells; independent validation passed 13/13 checks and six focused tests.
- Published a visually checked six-page HTML/PDF QC report. All source H5ADs
  remained unchanged; propagation and V4.2 retraining remain gated.

## 2026-08-17 - Added productive-chain UMAPs and rejected the provisional NK core

- Added boundary-subset UMAP overlays and exact all-cell tables for productive
  TRA, TRB, and TRD using nonempty harmonized productive-filtered CDR3 fields.
- Identified 371,616 productive-TRA, 373,186 productive-TRB, 373,339 TRA-or-TRB,
  and 83 productive-TRD cells among 475,953 boundary cells.
- Verified that harmonized TRA/TRB calls contain all 6,307 source-H5AD-positive
  boundary cells plus 364,091 current-atlas positives missing from the H5AD.
- Rejected the previous clusters 0/3/5 NK-core interpretation because
  190,757/257,569 (`74.06%`) carry productive TRA or TRB; also found productive
  TRA/TRB in 11,526/20,375 (`56.57%`) primary NK annotation anchors.
- Refreshed and visually checked the 21-page HTML/PDF report. No H5AD,
  classifier, model registry, or release artifact changed, and no GitHub push
  occurred.

## 2026-08-17 - Reintegrated the T/NK pool and resolved the NK-like boundary

- Applied a high-recall T/NK gate to 4,023,462 development cells, retained
  3,927,924 cells, all primary NK and repaired productive-T anchors, and removed
  95,538 off-target or no-direct-evidence cells without mutating source H5ADs.
- Recomputed 4,000 source-balanced HVGs after subsetting, excluded TCR V/J/D
  genes from ranking, forced 27 lineage/context genes, and completed 20-epoch
  A100 scVI plus nine global RAPIDS Leiden runs.
- Localized 96.77% of primary NK anchors to refined clusters 9 and 18, then ran
  nine second-pass Leiden partitions on their 475,953 cells.
- Identified seven boundary subclusters; the resolution-0.4 partition has mean
  seed ARI 0.962. Clusters 0, 3, and 5 form a 257,569-cell review-level NK core
  with 82.57% primary-anchor recall and 1.46% productive-T controls.
- Published and visually checked an 18-page HTML/PDF report. No classifier,
  training label, model promotion, atlas inference, or GitHub push occurred.

## 2026-08-17 - Extended V4.2 NK review with exact CD3/TCR expression

- Built T-lineage, NK-lineage, myeloid-context, and shared-cytotoxic feature
  UMAPs from the fixed 250,000-cell diagnostic sample and ranked the frozen
  24-cluster Leiden partition using all 4,023,462 cells.
- Added raw-UMI feature UMAPs for 29 CD3/TCR genes and exact all-cell detection
  tables by population, cluster, and source after confirming full gene coverage
  across all six registered source H5ADs.
- Found that cluster 19 contains 20,237/21,054 (`96.12%`) primary NK anchors and
  461 cells passing both current latent and expression evidence.
- Initially identified 461 dual-evidence cluster-19 cells, then held the entire
  set as provisional after 350 expressed a TRA/TRB constant-chain gene, 250 had
  at least two aggregate constant-chain UMIs, and 178 had at least three.
- Kept 594 latent-only, 759 expression-only, and 257 secondary-cluster cells as
  review tiers rather than automatically expanding NK labels.
- Added myeloid-context plots after showing that `FCER1G/TYROBP` alone would
  mis-rank cluster 1; its myeloid program exceeds its NK program and it contains
  no primary NK anchor.
- Published and visually checked the 17-page HTML/PDF report. No H5AD, label
  manifest, classifier, or GitHub release was changed.

## 2026-08-17 - Repaired and validated the V4.2 NK reference definition

- Restored 762,280 productive-TRA/TRB T anchors from the existing-atlas source
  object, removing the old-versus-new cohort confounding in the anchor set.
- Calibrated exact A100 50-neighbor latent evidence on held-out primary NK and
  productive-T anchors, with 93.23% NK recall and 0.0958% T-cell FPR at the
  fail-closed 0.98 threshold.
- Combined latent evidence with a strict NK-lineage expression rule that
  excludes detected core T-lineage genes and does not use cytotoxic genes to
  establish NK identity.
- Selected 469 `NK_CONFIDENT` cells from 113,287 eligible development
  candidates across four sources and capped effective future-training source
  contribution at 70%.
- Validated the fixed expression rule on locked GSE169246: 35.32% author-NK
  recall and 0.126% paired-alpha-beta-T FPR; no locked cell set a threshold.
- Published a visually checked six-page HTML/PDF report. No H5AD was changed,
  no classifier was fitted, and nothing from V4 was pushed to GitHub.

## 2026-08-16 - Added V4.2 sidecar visualization and confounding diagnostics

- Built a deterministic, source-balanced 250,000-cell diagnostic sample from
  the saved 4,023,462-cell scVI latent without modifying any H5AD.
- Added eight PNG diagnostics covering cohort composition, UMAP structure,
  anchor roles, candidate sources, near-qualifying clusters, and local source
  retention.
- Found complete anchor/cohort confounding: all 21,054 primary NK anchors came
  from the existing atlas and all 204,869 productive-TCR anchors came from the
  five new cohorts.
- Confirmed that the six NK-anchor-enriched clusters remained
  86.92%-89.69% dominated by GSE292700 candidates, so zero pseudo-NK cells are
  accepted and no classifier was fitted.
- Published a visually checked seven-page HTML/PDF report and froze GitHub
  publication of V4 until the model is scientifically finished and reviewed.

## 2026-08-16 - Completed V4.2 clustering/consensus as a negative result

- Activated the checksum-bound cluster-resource amendment and completed nine
  global plus nine boundary RAPIDS Leiden runs on the A100 without CPU fallback.
- Wrote a verified 19-MiB partition artifact for 4,023,462 cells; all six
  development source H5AD size/mtime pairs remained unchanged.
- Found that the expression-independent boundary included 3,978,014 cells
  (`98.87%`), so the broad-cluster union did not localize the T/NK boundary.
- Evaluated 396 clusters: six passed 95% anchor purity and 2% productive-TCR
  contamination, but all were 86.92%-89.69% dominated by GSE292700 candidates
  and failed the frozen 70% source cap.
- Selected zero of 113,287 eligible pseudo-NK cells and stopped before any
  classifier, calibration, threshold, release artifact, or atlas inference.
- Added a reproducible report builder and visually checked a print-safe
  five-page HTML/PDF QC package.

## 2026-08-14 - Completed V4.2 staging/scVI and audited cluster resources

- Activated the user-approved recovery contract and staged 4,023,462
  development cells x 4,000 HVGs with zero locked-cohort cells.
- Completed the frozen 20-epoch scVI fit on the A100 80-GB GPU in 1,358
  seconds with no CPU fallback; saved a finite 30-dimensional latent matrix.
- Retained unchanged size/mtime for all six development source H5ADs.
- Stopped before RAPIDS graph construction when shared SSD free space fell
  below the original 300-GiB floor; did not stop or delete the unrelated BAM
  workflow responsible for the transient pressure.
- Added a checksum-bound stage-specific storage contract that retains 300 GiB
  for staging/fitting and proposes 150 GiB for clustering/consensus.
- Passed 18/18 amendment checks with a 2.33-GiB worst-case output reserve and
  108.4-GiB post-reserve margin; clustering and classifier fitting remain
  unperformed pending separate approvals.

## 2026-08-14 - Audited V4.2 current-atlas recovery after SSD loss

- Stopped the approved project-data command before matrix reading when the
  SSD-resident `integrated.h5ad` was found to be absent.
- Added a fail-closed recovery path from the intact 3,705,384-cell canonical
  `TNK_cleaned.h5ad`, using the checksum-pinned 78-cell Phase-3 exclusion
  manifest and original harmonized metadata sources.
- Reconstructed exactly 3,705,306 effective cells x 27,413 genes, reproduced
  the frozen five-column metadata audit, and recovered all 21,054 primary
  dual-annotation NK anchors.
- Passed 20/20 read-only recovery checks and published a visually reviewed
  three-page HTML/PDF package. No project-data integration, scVI fitting,
  clustering, pseudo-labeling, or classifier fitting ran.
- Invalidated the former execution approval because the physical input and
  checksum-bound code/config contract changed; recovery execution remains at a
  new explicit supervision gate.

## 2026-08-08 - Passed gdTAI V4.2 sidecar-integration implementation QC

- Implemented a checksum-bound, stage-addressable integration runner that
  excludes all three locked cohorts before development metadata or expression
  loading and aborts before SSD access when project-data approval is absent.
- Passed 11/11 deterministic core tests, including CSR/CSC H5AD slicing,
  source-balanced sampling, batch fallback, role leakage, boundary definition,
  and cluster-consensus pseudo-label rules.
- Passed 15/15 implementation checks; synthetic scVI used `cuda:0` on the A100
  80-GB GPU and RAPIDS completed neighbor construction plus repeated Leiden
  clustering without CPU fallback.
- Verified all nine source H5AD size/mtime pairs unchanged. No project data was
  integrated or fitted, and classifier fitting remains a separate gate.

## 2026-08-07 - Completed gdTAI V4.1-GPU Gate C as a negative result

- Completed all three checksum-bound nested outer decision paths on the A100;
  99/99 recorded inner-fold fits converged and 57 complete threshold frontiers
  were retained.
- Found 0/15 eligible Stage-1 candidate evaluations. The literal worst-source
  NK passage was 100%, and the best worst-source passage among NK strata with
  at least 100 cells remained 90.97%-93.23% at recall-preserving thresholds.
- Found no eligible balanced or high-purity operating point for compact
  seven-gene ridge, V2-like TCR-gene ridge, or raw legacy TRD-minus-TRAB score.
- Marked V4 Stage 2, extension-model controls, and paired bootstrap as not
  estimable because no Stage-1 model passed; no V4.1 artifact was promoted,
  released, or applied to the atlas.
- Published the Gate C HTML/PDF report and retained V3 Round 14 as the promoted
  balanced default.

## 2026-08-07 - Passed gdTAI V4.1-GPU Gate B

- Implemented deterministic PyTorch ridge and GPU XGBoost fitting primitives,
  complete Stage-1/Stage-2 threshold frontiers, portable NPZ/UBJ exports, and
  hash-bound atomic candidate checkpoints.
- Passed 32/32 Gate B checks, including full frozen-cache hashes,
  bit-identical repeated A100 probabilities for both families, exact portable
  artifact reloads, and checkpoint resume.
- Published the print-safe Gate B HTML/PDF QC package; no project cell was
  fitted and Gate C remains absent.

## 2026-08-07 - Retired V4 CPU evaluation and froze GPU redesign proposal

- Preserved completed HRA005041 and GSE144469 V4 v1.2 CPU outer-fold artifacts
  as diagnostic negative results; neither fold contained a valid operating
  threshold for V4 or the fair comparators.
- Terminated the unfinished BALF_BLOOD_COPD CPU fold at the user's direction.
- Verified synthetic-only A100 feasibility, deterministic PyTorch weighted
  ridge logistic regression, deterministic GPU XGBoost, model export, and
  isolated direct-CUDA operation without touching another user's MPS daemon.
- Added the pre-fit `GDTAI_V4_GPU_PRECOMMITTED_PLAN.md`; project-data GPU
  implementation and fitting remain blocked pending explicit review.

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
- GSE305372 external application retired on 2026-08-06
  - removed the cohort from active dataset, library, and file registries
  - archived dataset-specific code, methodology, rendered reports, and local
    result packages
  - retained GEO URLs, byte sizes, SHA-256 checksums, and Git rollback commit
  - deleted the 16,750,876,819-byte local dataset directory by explicit user
    request
- extension-cohort standalone Phase 0 review completed on 2026-08-06
  - seven checksum-pinned H5ADs built with 498,901 cells and 154 libraries
  - all seven passed sparse raw-count, metadata, TCR-schema, logical-flag, and
    sample-plus-barcode structural QC
  - tumor-project specimen context resolved for every cell with zero violations
    and without modifying source H5ADs
  - staged-source audit confirmed 1,267 nonproductive TRG/TRD fragments and
    zero productive TRG/TRD CDR3s
  - GSE169246 remained provenance-blocked; Tan et al. 2021 remained excluded
    as the existing `GDT_2020AUG_woCOV` study
  - supervision report and source-audit addendum emailed to Likai
  - merge, integration, and model-evaluation phase advancement remain pending
    explicit review approval
- GSE169246 processed-source selection registered on 2026-08-06
  - copied the user-selected 394,225-cell `TNK_cleaned.h5ad` unchanged into
    `data/datasets/GSE169246/processed/artifacts/`
  - source and canonical copy matched by full SHA-256
  - registered the dataset as inactive and pending standalone Phase 0
  - updated extension guards to accept only this processed object while
    keeping stage, build, merge, and integration disabled
  - documented incomplete published-annotation coverage as a required review
    item rather than treating the object as the complete unfiltered study
- extension-cohort standalone Phase 1 T/NK filtering completed on 2026-08-06
  - independently filtered all eight new inputs from 893,126 cells to 758,135
    T/NK candidates without merging cohorts
  - retained all 338,046 cells with productive TRA/TRB/TRG/TRD evidence
  - tightened transcriptomic T-cell evidence so `IL7R` or `LTB` alone cannot
    pass the lineage gate, while retaining the multi-marker NK rule
  - passed metadata completeness, exact source-accession, cohort-ID, unique
    cell/library-barcode, canonical TCR-flag, and specimen-context checks for
    all eight outputs
  - resolved GSE169246 barcode reuse by distinguishing 78 source libraries
    from 51 biological samples; 76 libraries retained T/NK candidates
  - verified source H5AD SHA-256 values were unchanged
  - left merge, integration, and gdTAI evaluation blocked pending explicit
    user review
- extension frozen gdTAI profile screen completed on 2026-08-06
  - applied V2 high-F1, V2 high-purity, V3 Round 14 balanced, and V3 Round 12
    high-purity to all 758,135 independently filtered extension T/NK cells
  - verified all frozen model artifacts and all eight source H5AD SHA-256
    values; no H5AD was modified
  - accepted productive single TRA or TRB cells as alpha-beta negative controls
    and audited paired alpha-beta, single-chain alpha-beta, strict NK, donor,
    library, source-GSE, and TRDC/TRDV failure-mode strata
  - V3 Round 12 high-purity minimized pooled alpha-beta FPR (`0.227%`), while
    V2 high-purity minimized strict-NK FPR (`0.725%`)
  - complete-schema sensitivity excluding GSE169246 favored V3 Round 12 for
    the known-negative union (`0.295%`), but historical positive sensitivity
    continued to favor V3 Round 14
  - made no model selection or promotion because the extension cohorts contain
    no unbiased gdT positive truth set
  - retained merge and integration as blocked review gates
- extension official-GEO metadata reconciliation completed on 2026-08-06
  - reconciled 58 grouped metadata assertions across 8 extension cohorts and
    10 official GEO accessions using direct NCBI evidence
  - classified 16 rows as locally verified, 32 as GEO-resolved but locally
    unresolved, 6 as partially ambiguous, and 4 as unavailable in GEO
  - identified all 239,418 retained GSE169246 `_b` cells as blood libraries
    whose current joined tissue/specimen fields are semantically incorrect
  - validated full accession coverage and unchanged file state for all eight
    filtered H5ADs without writing metadata
  - emailed the reconciliation report, audit log, and evidence table to Likai
  - left additive correction, merge, and integration blocked for explicit
    review
- independent gdTAI training and evaluation audit completed on 2026-08-06
  - a separately instantiated read-only reviewer and main-agent provenance
    inventory independently confirmed that V2 validation and BALF_BLOOD_COPD
    results were used for model selection
  - reclassified BALF_BLOOD_COPD as a reused cross-study benchmark rather than
    an untouched independent external test
  - recomputed expression-independent TCR-only metrics on 1,046 paired
    productive gdT/no-abT positives and 33,117 strict paired abT/no-gdT
    negatives; Round 12 and Round 14 F1 were effectively tied at `0.8929` and
    `0.8922`
  - verified every registered model SHA-256 and documented a separate V3
    semantic provenance mismatch across pickle, manifest, completed-run
    summary, and training composition
  - published a 17-item issue register, prevalence scenarios, extension
    negative-control synthesis, and a precommitted nested grouped-evaluation
    program without opening or modifying an H5AD
  - corrected V2/V3 methodology and registry wording; retained Round 14 as the
    current promoted artifact pending a separately approved reevaluation
  - emailed the summary and comprehensive PDF audit report to Likai
- gdTAI V4 Step 0 protocol frozen on 2026-08-07
  - fixed expression-independent TCR truth rules and T/NK-only input contract
  - froze a 197-gene feature universe with no Phase-4 scores and no B/myeloid
    features
  - precommitted nested leave-one-dataset-out evaluation with grouped inner
    folds, fold-local calibration and thresholding, paired-abT/NK guardrails,
    and a bootstrap-tested promotion rule
  - rendered a print-safe seven-page HTML/PDF review package
  - amended the protocol to require a pre-training CD4/Treg exclusion recall-
    ceiling audit and to prohibit cutoff retuning after Step 2 begins
  - left all training, threshold search, whole-atlas inference, and model
    promotion blocked for supervision
- gdTAI V4 Step 1 preflight and grouped split freeze completed on 2026-08-07
  - archived the earlier experimental V4 bundle intact with verified SHA-256
  - froze a 1,137,739-cell expression-independent label manifest and
    deterministic grouped inner-split manifest
  - mapped all 107,068 GSE144469 cells uniquely by cell-ID SRR plus barcode and
    corrected its raw expression source to `layers/counts`
  - passed fixed CD4/Treg feasibility with a 98.38% source-macro primary-gdT
    recall ceiling, well above the 80% guardrail
  - failed the training gate because HRA005041 has no raw-count layer and the
    GDT2020/Malte supplements cover only 82.7%/41.6% of frozen features
  - verified all 16 input files remained unchanged and all registered hashes
    matched; performed no fitting, threshold search, inference, or promotion
  - published and visually checked the HTML/PDF supervision package; Step 2
    remains blocked
- gdTAI V4 protocol v1.2 amendment frozen on 2026-08-07 before fitting
  - retained the frozen 197 genes, expression-independent labels, nested folds,
    model grids, guardrails, CD4/Treg thresholds, and promotion criteria
  - permitted only the registered HRA005041 `log1p(CP10K)` representation,
    conditional on a full per-cell `sum(expm1(X))` audit at a frozen tolerance
  - moved GDT_2020AUG_woCOV, MalteGDT, and GDTlung to sensitivity-only roles
    with zero weight in both training stages
  - corrected the preflight contract so silver and dual/ambiguous cells cannot
    enter Stage 1 fitting
  - regenerated and visually checked the protocol HTML/PDF; no fitting,
    threshold search, inference, or promotion was performed
- gdTAI V4 protocol v1.2 Step 1 preflight passed on 2026-08-07
  - passed all frozen input, feature, truth-label, exclusion-feasibility,
    grouped-split, checksum, and file-state checks with zero warnings
  - fully audited all 766,639 HRA005041 cells and found zero empty rows and zero
    reconstructed-library sums outside the frozen `1e-4` tolerance
  - froze the 1,137,739-cell label manifest and deterministic three-fold
    grouped split with recorded SHA-256 values
  - visually checked the 13-page A4 supervision report and verified its
    artifact checksum inventory
  - performed no fitting, calibration, threshold search, inference, or
    promotion; Step 2 remains blocked for explicit supervision approval
- gdTAI V4 Step 2 nested evaluator implemented on 2026-08-07
  - added a checksum-bound authorization gate before any project-data cache or
    fitting operation
  - implemented nested leave-one-dataset-out fitting, fold-local feature
    filtering, weighting, calibration, frozen threshold guardrails, fair
    comparators, hierarchical bootstrap, feature-stability checks, extension
    negative controls, and print-safe HTML/PDF reporting
  - passed 17 deterministic synthetic unit/integration tests, including a full
    outer fold through every candidate model path
  - corrected the cache exclusion audit to reproduce the exact Step 1 scope
    before applying the unchanged frozen CD4/Treg rules to all evaluation rows;
    the initial cache attempt stopped at this audit before any fitting
  - after an incomplete timing run produced no outer-fold result, added
    deterministic fold-level parallelism and explicit iteration/convergence
    provenance without changing the frozen model or selection contract
  - this is an implementation checkpoint only; no project-data model had been
    fitted and no V4 model had been promoted at this checkpoint
- gdTAI V4.2 Step 0 NK-label audit completed on 2026-08-07
  - hash-verified the 1,137,739-cell label manifest, 197-gene cache, and all 15
    completed V4.1 Stage-1 saved-OOF checkpoint paths
  - found that 315,726/336,780 (`93.7%`) V4.1 NK controls had only a single
    scVI annotation, while 21,054 (`6.3%`) had dual scVI/author agreement
  - showed that all 15 candidate/fold combinations pass the frozen Stage-1 NK
    source cap when evaluated against dual-annotation controls at unchanged
    per-source T-recall thresholds
  - retained expression-based T-lineage coherence as a diagnostic only and
    changed no labels or H5AD files
  - published a four-page HTML/PDF report and a proposed V4.2 repair protocol;
    no integration, clustering, fitting, threshold search, or promotion ran
  - added a supervision proposal to build a separate new-data NK modeling
    integration with whole-dataset development/validation separation
- gdTAI V4.2 modeling-integration preflight completed on 2026-08-08
  - froze five extension development cohorts, the current atlas reference, and
    three whole-cohort locked evaluation roles before any integration
  - hash-verified all nine input H5ADs, sampled raw-count state, required sparse
    matrices and metadata, and confirmed unchanged size/mtime after the audit
  - mapped all 21,054 weight-1 dual-annotation NK anchors exactly to the current
    atlas and verified complete Stage-1 gene coverage in every cohort
  - found 14,265 common development genes and 13,975 eligible genes after
    excluding TCR, mitochondrial, ribosomal, and immunoglobulin genes
  - classified GSE169246 final-cascade evaluation as reduced-feature
    sensitivity because it has 145/197 genes despite complete 50/50 Stage-1
    coverage
  - passed the 800-GiB RAM, 300-GiB SSD, and 75-GiB GPU gates with estimated
    peak RAM 223.6 GiB, 875.1 GiB SSD free, and an idle A100 80 GB
  - published and visually checked a seven-page HTML/PDF report; no integration,
    clustering, pseudo-labeling, fitting, or promotion ran

## Current working era

- post-Phase 4 downstream analysis and documentation cleanup
