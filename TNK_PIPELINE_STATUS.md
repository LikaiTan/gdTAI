# TNK_PIPELINE_STATUS.md

## Current milestone

- the full T/NK atlas rebuild completed on 2026-08-18 with all planned stages
  passing: preflight, sparse preparation, A100 scVI, RAPIDS embedding, full-gene
  assembly, Phase 4 scoring, and HTML/PDF reporting
- the rebuilt atlas contains exactly 5,933,312 cells x 27,413 genes with
  9,127,088,723 sparse nonzero counts: 5,128,904 historical-atlas cells,
  758,135 cells from eight new T/NK-filtered cohorts, and 46,273
  `BALF_BLOOD_COPD` cells
- 16 physical input cohorts contribute 40 source accessions; A100 scVI produced
  a finite 5,933,312 x 30 latent matrix and RAPIDS produced a finite UMAP plus
  33 Leiden clusters
- the canonical rebuilt H5AD is
  `Integrated_dataset/integrated_full_atlas.h5ad`, linked to the 22.1-GB SSD
  object under `/ssd/tnk_phase3/Integrated_dataset/full_atlas/`
- the final scored H5AD SHA-256 is
  `f5fc491a70f12adeeda5764cb116a59bb441460285905f297bbc2ff691559802`
- all eight continuous Phase 4 score columns are finite for every cell; all 16
  source H5AD size/mtime pairs remained unchanged
- the rebuild intentionally did not propagate the validated 14-source TCR
  repair sidecar; the report also blocks GSE169246 tissue/specimen analysis
  pending its separately reviewed additive `_b` blood-library correction
- sample-aware productive TCR join rebuild completed on 2026-08-18 with
  `PASS_SIDECAR_READY`; independent validation passed 15/15 checks
- all 14 flagged sources passed a full or fail-closed partial rebuild, yielding
  3,041,871 staged RNA rows, 1,121,858 cells with validated productive TCR,
  and 1,789,643 productive chain calls with observed UMI support
- TRA, TRB, TRG, and TRD were rebuilt from productive, cell-associated,
  high-confidence raw records using only `sample_id + barcode_core`; the
  highest-UMI contig per chain was selected, with reads and full-length status
  used as deterministic tie-breakers
- per-chain UMI/read values and support-availability flags are retained;
  unavailable quantitative support remains null rather than zero
- GSE125527 was repaired with the published patient remap plus tissue,
  GSE228597 with pooled-library suffixes, and GSE287541 with all 46 public TCR
  SRA runs reconstructed by Cell Ranger VDJ using round plus participant visit
- GSE235863 passed partially with 110 duplicated RNA join-key rows blanked and
  excluded from TCR truth; no source-level quarantine remains
- the checksum-bound replacement sidecar has 3,041,871 rows from 14 sources
  and SHA-256
  `3114e70719301d693ae1a2bc2c63bac6c8bd57e3e8ac73a88c24320eaabfc2f0`
- all 14 recorded source H5AD size/mtime pairs remained unchanged; no TCR
  sidecar has yet been propagated into source or milestone H5AD metadata
- read-only sidecar overlay preflight passed `PASS_OVERLAY_READY`: all
  2,155,409 atlas rows from the 14 affected sources match exactly by
  `source_gse_id + original_cell_id`, and no affected atlas row lacks a sidecar
- after repair, the retained atlas contains 1,100,981 productive-alpha-beta
  cells and 933,890 paired-TRA/TRB cells from these sources
- combining repaired and unaffected sources, the metadata-overlaid atlas is
  estimated to contain 2,270,138 productive-alpha-beta cells and 1,938,158
  paired-TRA/TRB cells
- 20,875 sidecar-supported productive-alpha-beta cells outside the atlas,
  including 14,495 paired cells, will not be rescued; they represent only
  0.911% and 0.742% of the whole available productive and paired alpha-beta
  pools, respectively, and none has productive gamma-delta TCR
- the TCR repair will therefore preserve the exact 5,933,312-cell membership,
  expression matrix, scVI latent space, Leiden clusters, and UMAP; only TCR
  metadata and TCR-derived truth/control fields require replacement
- the existing boundary-conflict, NK-anchor, and cluster summaries below were
  computed before sidecar propagation and remain diagnostic historical
  evidence; they must be regenerated after the approved metadata replacement
- gdTAI V4.2 productive-TRA/TRB boundary conflict audit completed on
  2026-08-18 with `PASS_JOIN_REPAIR_REQUIRED`
- 373,339/475,953 boundary cells (`78.44%`) have productive TRA or TRB, but
  277,305/373,339 (`74.28%`) of those calls come from 14/29 contributing
  sources with a documented or strong TCR-join red flag
- the red flags include a documented legacy barcode-only mapping mode, exact
  TRA/TRB-pair reuse across multiple sample or donor keys, and broad
  alpha/beta calls in author-labeled non-T lineages; these defects are present
  in source objects and therefore precede scVI, UMAP, and Leiden
- ambient one-UMI calls are not the dominant detected explanation in the
  quantitative subset: 47/3,047 (`1.54%`) UMI-auditable productive-AB cells
  have no chain above one UMI, but UMI support is unavailable for 370,292
  productive-AB cells and remains an explicit uncertainty
- the 4,000-gene integration retained all 27 forced T/NK-context genes;
  boundary clusters associate more strongly with expression class
  (`NMI=0.0755`) than with raw productive-AB status (`NMI=0.0241`), so
  re-clustering alone cannot repair unsafe metadata joins
- only 1,508 cells meet the audit's conservative residual AB-compatible
  stratum: paired TRA/TRB, no detected source join red flag, at least one chain
  with at least two observed UMIs, and multi-gene T-lineage support; this is a
  diagnostic stratum and not replacement ground truth
- gdTAI V4.2 T/NK-restricted reintegration and second-pass NK-boundary review
  completed on 2026-08-17 with `PASS_REVIEW_REQUIRED`
- the high-recall gate retained 3,927,924/4,023,462 cells (`97.63%`), retained
  all 21,054 primary NK anchors and all 967,149 repaired productive-T controls,
  removed all flagged doublets, and changed no source H5AD
- 4,000 HVGs were recomputed only after T/NK restriction using source-balanced
  sampling; TCR V/J/D genes were excluded from HVG ranking and a 27-gene
  T/NK-context panel was forced into the integration feature set
- A100-only scVI and RAPIDS Leiden produced 22 global refined clusters;
  clusters 9 and 18 contain 475,953 cells and 20,375/21,054 (`96.77%`) primary
  NK anchors but also 9,248 productive-T controls, so neither is NK wholesale
- nine second-pass Leiden runs divided the boundary into seven review
  subclusters; the resolution-0.4 partition is seed-stable with mean adjusted
  Rand index `0.962`
- chain-specific boundary UMAPs now use nonempty harmonized productive-filtered
  CDR3 metadata: 371,616 cells have productive TRA, 373,186 have productive
  TRB, 373,339 (`78.44%`) have either, and 83 have productive TRD
- the productive-chain audit rejects the previous clusters 0/3/5 NK-core
  interpretation: 190,757/257,569 (`74.06%`) carry productive TRA or TRB, and
  11,526/20,375 (`56.57%`) primary NK annotation anchors do so as well
- all 6,307 source-H5AD TRA/TRB-positive boundary cells agree with harmonized
  metadata, which recovers 364,091 additional current-atlas positives whose
  chain fields were not propagated into `TNK_cleaned.h5ad`
- no classifier, threshold, label manifest, model promotion, release artifact,
  atlas inference, source-H5AD mutation, or GitHub push occurred
- full-atlas metadata harmonization v2 completed its read-only rule, sample
  identity, and row-level overlay gates without modifying an H5AD
- the complete 5,933,312-row overlay passed 18/18 checks and is checksum-bound
  at `4da4eea32e9d275de790775e2a0f59d4f6553d72756e9c8c6935f35bb398984f`
- unresolved tissue/specimen context fell from 3,849,174 blank current values to
  39,920 fail-closed cells: 36,800 source-blank GSE206325 rows and 3,120 mixed
  duodenum/PBMC GSE252762 rows
- GSE125527 recovers 30 source-derived donor-by-tissue samples, GSE254249
  recovers 92 source `Ident` samples, and GSE228597 retains 4,611 explicitly
  unresolved biological samples with separate technical VDJ-library identity
- no cell-type label is used as tissue/specimen evidence, and no TCR chain call
  has been applied to the metadata overlay
- the approved metadata-only H5AD transaction passed all 19 post-write checks;
  the validated 5,933,312-cell candidate is
  `/ssd/tnk_phase3/Integrated_dataset/full_atlas/metadata_corrected/integrated_full_atlas.h5ad`
- the candidate SHA-256 is
  `7f1c5e1cac1074a8e2863703bc1862e225defc5ba1a3adbaabd3f6e023d5871c`;
  all 13 additive columns reproduce the keyed overlay's ordered hashes and
  value counts exactly
- sparse X, X_scVI, X_umap, existing obs/var signatures, cell order, and all
  legacy TCR fields remain unchanged; AnnData backed-open validation passed
- the canonical atlas still resolves to the original checksum-bound rollback
  object, and no TCR chain call has been applied
- the separate TCR sidecar transaction completed with
  `PASS_TCR_H5AD_CANDIDATE`; all 17 post-write checks passed on
  `/ssd/tnk_phase3/Integrated_dataset/full_atlas/tcr_corrected/integrated_full_atlas.h5ad`
- the TCR-corrected candidate SHA-256 is
  `d32c9d2bdb955b12e1eafbed8322f8cb965cf3a225191e612b53f3d3783480d5`
- all 2,155,409 affected rows from 14 sources match the validated sidecar, all
  unaffected TCR values match their input hashes, and the 109 atlas-present
  ambiguous GSE235863 rows are blank and fail closed
- corrected whole-atlas totals are 2,270,138 productive-alpha-beta cells and
  1,938,158 paired-TRA/TRB cells; 46 canonical fields were replaced and 27
  explicit rebuilt-provenance fields were added
- sparse X, scVI, UMAP, protected metadata, and the metadata-only input checksum
  remain unchanged; the canonical symlink was not switched
- the eight-page HTML/PDF application report passed visual review without
  clipped tables and is under
  `gdT_prediction/gdtai_v4_2_tcr_sidecar_application/`

## Next action

- regenerate TCR-derived truth/control labels and every boundary/NK audit from
  the checksum-bound TCR-corrected candidate
- keep biological `sample_id_harmonized_v2` separate from technical
  `tcr_library_join_id_v2`; pooled VDJ libraries must never overwrite specimen
  or donor identity
- treat the original canonical atlas and metadata-only candidate as rollback
  points; do not switch the canonical symlink until post-repair truth and
  boundary reports are reviewed
- preserve the frozen 5,933,312-cell membership and existing integration; do
  not rerun T/NK selection, HVGs, scVI, Leiden, or UMAP merely to rescue the
  intentionally omitted 20,875 alpha-beta controls
- after propagation, regenerate TCR-derived truth/control labels and every
  boundary/NK audit that consumed stale productive-chain metadata
- withdraw clusters 0, 3, and 5 as an NK candidate core; no boundary subcluster
  is currently eligible as NK training truth
- do not build model truth directly from the current harmonized productive-CDR3
  fields; first replace red-flag source assignments with validated raw-contig
  joins carrying exact source, library, sample, cell, UMI, and read provenance
- reconstruct primary NK anchors only after validated TCR replacement, then
  audit author annotations, doublets, ambient RNA, and assay coverage
- repeat source-balanced cluster-label review only after the TCR sidecar and NK
  anchors pass exact overlap and missingness checks
- run label-leakage, source-balance, feature-coverage, and grouped-fold
  preflight before fitting; keep primary NK anchors and productive-T controls as
  the dominant training evidence
- if a reviewed subset later passes, cap its effective per-source contribution
  at 70%; do not relax the cap or redefine NK from cytotoxicity or adaptor genes
- compare any resulting model against V2/V3 on unchanged locked cohorts before
  considering promotion
- keep all locked cohorts unchanged for nested comparison with V2 and V3;
  model promotion, release fitting, and whole-atlas inference remain high-risk
  decisions requiring explicit approval
- do not authorize Gate D release fitting, V4.1 promotion, or whole-atlas
  inference from this failed experiment
- do not promote a model, fit a release artifact, or run whole-atlas inference
  from the incomplete CPU evidence or the GPU feasibility probes
- do not promote or claim superiority for another gdTAI model before grouped
  resampling, expression-independent labels, and fold-local selection are in place
- review the frozen-profile negative-control report and official-GEO metadata
  reconciliation package
- decide separately whether to approve additive corrected metadata fields for
  GSE169246 and the other GEO-resolved rows
- treat the completed full-atlas merge as a structural integration milestone;
  do not use GSE169246 tissue/specimen fields for biological comparisons until
  its additive blood-compartment correction is reviewed and applied
- do not push any gdTAI V4 code, reports, or artifacts to GitHub until V4 is
  scientifically finished and reviewed; local commits remain permitted

## Current blockers or review items

- the rebuilt TCR sidecar is present only in the validated isolated
  TCR-corrected candidate; canonical publication remains pending post-repair
  truth/control and boundary/NK review
- the three source-level quarantine blockers are resolved: GSE125527 has
  75.82% exact author/raw CDR3 agreement with 22,928-fold enrichment over
  sample rotation; GSE228597 has 95.03% agreement and 7,951-fold enrichment;
  GSE287541 has 46/46 reconstructed TCR libraries, 99.72% agreement, and
  67,913.5-fold enrichment
- GSE125527's published productive-receptor tables do not contain UMI/read
  columns, so quantitative support remains null; reconstructing its 71 public
  TCR runs (150.6 GiB compressed) is not required for the deterministic join
- GSE235863 has 110 duplicated RNA `sample_id + barcode_core` rows; these rows
  are replacement-eligible for clearing stale calls but ineligible for receptor
  truth
- the previous source-H5AD control mask was incomplete because productive-chain
  metadata from many datasets had not been propagated into `TNK_cleaned.h5ad`
- 373,339/475,953 boundary cells carry harmonized productive TRA or TRB, while
  only 83 carry productive TRD; lack of TRD is uninformative for libraries that
  did not perform gamma-delta V(D)J sequencing
- 56.57% of primary NK annotation anchors conflict with productive TRA/TRB,
  including all 3,482 GSE125527 boundary anchors and 8,044/16,893 GSE228597
  boundary anchors; those anchors cannot independently validate NK identity
- all 83 productive-TRD boundary cells also have productive TRA or TRB and need
  doublet, multiplet, and chain-assignment review
- boundary clusters remain useful transcriptomic neighborhoods, but none can be
  interpreted as NK truth from the current anchor and chain metadata
- the earlier 461-cell strict-expression core remains invalid for training
  because its zero-CD3 criterion was circular and TRA/TRB constant-chain UMIs
  were too frequent; the new review does not reinstate those labels
- the missing historical integrated H5AD prevents recovery of the complete
  full-cell scANVI annotation column; the annotation UMAP therefore displays
  only recoverable, provenance-backed annotation evidence
- the original `high_speed_temp/Integrated_dataset/integrated.h5ad` is absent,
  and no exact copy was found in the project or databank search; the previously
  recorded project-data execution approval is therefore invalid
- the recovery preflight passed 20/20 checks using the intact canonical
  3,705,384-cell `TNK_cleaned.h5ad`, the exact 78-row Phase-3 exclusion
  manifest, and the original harmonized metadata sources
- the recovered effective input contains exactly 3,705,306 cells and 27,413
  genes, reproduces all five frozen metadata missing/unique-count audits, and
  recovers all 21,054/21,054 primary dual-annotation NK anchors
- recovered sparse staging passed on 4,023,462 development cells x 4,000 HVGs
  with 1,667,132,819 nonzero raw counts, all 21,054 primary NK anchors, 204,869
  productive-TCR anchors, and zero locked-cohort cells
- A100 scVI fitting passed in 1,358 seconds with no CPU fallback and produced a
  finite 4,023,462 x 30 latent matrix with SHA-256
  `7030b28be885b85a4b24ce2223bb16f7c7288a29195f2512bf7bbba330791ccf`
- the first staging attempt assembled the complete sparse matrix but AnnData
  rejected nullable-string categorical serialization; the unchanged approved
  runner completed after enabling its documented
  `ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS=1` compatibility setting
- RAPIDS clustering did not start because shared SSD free space fell below the
  original 300-GiB floor; an unrelated BAM sort completed and self-cleaned,
  after which free space stabilized near 261 GiB
- the cluster resource-amendment preflight passed 18/18 checks, retained 300
  GiB for `prepare` and `fit`, proposed 150 GiB for `cluster` and `consensus`,
  bounded worst-case cluster output plus reserve at 2.33 GiB, and retained
  108.4 GiB post-reserve margin
- all six development source H5AD size/mtime pairs remained unchanged through
  clustering and consensus; classifier fitting remains unperformed
- the repaired anchor audit found 762,280 cells with productive TRA or TRB
  evidence in the existing atlas across seven sources; with new-cohort anchors,
  the sidecar now contains 967,149 productive-T anchors across 12 sources
- exact A100 brute-force nearest-anchor scoring was repeated internally with
  bit-identical neighbor indices and distances; the fixed 50-neighbor NK
  fraction threshold is `0.98`, with a 99th-percentile known-NK distance cap of
  `3.1079669`
- held-out primary-NK latent recall was `93.23%`, and held-out productive-T
  latent FPR was `0.0958%`; score ties fail closed
- the strict expression rule requires no detected
  `CD3D/CD3E/CD3G/TRAT1/BCL11B`, both `FCER1G` and `TYROBP`, and at least one of
  `KLRD1/NCR1/FCGR3A/KLRC1`
- locked GSE169246 did not set any rule or threshold: the fixed expression rule
  recalled `35.32%` of 7,770 author NK cells and passed `0.126%` of 54,925
  paired alpha-beta T cells; a cytotoxic-only rule passed `98.21%` and `30.40%`,
  respectively, confirming that cytotoxicity cannot define NK identity
- the repaired audit manifest contains 469 `NK_CONFIDENT` cells from four
  sources and remains preserved unchanged; the later raw TCR review holds all
  rows from training, including the 461-cell cluster-19 provisional core
- no source H5AD was mutated, no classifier was fitted, and no V4 artifact was
  pushed to GitHub in the NK-definition repair
- changing the stage-specific resource contract invalidated the earlier
  execution approval; the user approved the amended checksum-bound cluster and
  consensus scope on 2026-08-16
- V4.2 implementation QC passed 15/15 checks; the synthetic 1,600-cell scVI
  fit ran on `cuda:0` using the A100 80-GB GPU, produced a finite 8-dimensional
  latent matrix, and completed RAPIDS neighbors plus two Leiden seeds
- 11/11 deterministic tests passed for approval binding, role leakage,
  fail-closed execution, source-balanced sampling, CSR/CSC H5AD slicing,
  batch-key fallback, expression-independent boundary construction, and
  consensus pseudo-label rules
- all nine source H5AD size and nanosecond-mtime pairs remained unchanged; no
  project-data integration, project-data model fit, or classifier fit occurred
- the earlier `PROJECT_DATA_INTEGRATION_APPROVAL.json` no longer authorizes
  execution because the execution config, runner, core, and physical input
  contract changed for recovery
- V4.2 integration preflight passed 59/59 checks on 4,023,462 development cells
  and 439,979 whole-cohort locked-evaluation cells
- all 21,054 primary dual-annotation NK anchors mapped to the current atlas;
  every development and locked cohort has 50/50 Stage-1 genes
- GSE169246 has only 145/197 full classifier genes, so its 7,770 author NK cells
  are a primary Stage-1 passage challenge but final-cascade FPR is explicitly a
  reduced-feature sensitivity analysis
- GSE315928 retains 66,813 paired-alpha-beta controls with 197/197 genes;
  GSE121636/GSE121637 remains a mixed stress cohort rather than author NK truth
- the development intersection has 14,265 genes and 13,975 eligible non-TCR,
  non-mitochondrial, non-ribosomal, non-immunoglobulin genes before selecting
  4,000 HVGs
- the resource audit estimated 223.6 GiB conservative peak RAM, verified 875.1
  GiB free SSD, and found an idle A100 80-GB GPU
- all nine input H5ADs were hash-verified and retained identical size and
  nanosecond modification time; no H5AD was modified
- the checksum-bound clustering/consensus execution is complete; classifier
  fitting from this lane is scientifically blocked because no pseudo-NK cells
  were selected
- the deterministic 250,000-cell latent-space diagnostic confirmed that the
  sidecar is a development-only integration, not a new canonical atlas
  milestone; it contains 3,705,306 existing-atlas cells and 318,156 new-cohort
  cells
- primary-NK and productive-TCR anchor roles are perfectly confounded with the
  old/new cohort lane, so apparent anchor purity may partly represent cohort
  structure rather than NK-versus-T biology
- median same-source 30-neighbor retention was 46.67% for the existing atlas,
  86.67% for GSE114724, 100% for GSE159251, 100% for GSE292700, 86.67% for
  GSE294273/4, and 90% for GSE296954; these are descriptive integration
  diagnostics, not biological pass/fail thresholds
- the six near-qualifying clusters were 86.92%-89.69% GSE292700 among eligible
  candidates; zero pseudo-NK cells remain accepted, no classifier was fitted,
  and no H5AD was mutated

- V4.2 Step 0 verified 336,780 V4.1 NK controls: 315,726 (`93.7%`) had only a
  single scVI NK annotation and 21,054 (`6.3%`) had dual scVI/author agreement
- all 15 checksum-bound saved-OOF candidate/fold evaluations passed the frozen
  50% NK source cap under the dual-annotation truth counterfactual; worst-source
  passage among sources with at least 100 cells was `5.27%-19.37%`, versus
  `90.97%-99.21%` using all V4.1 controls
- the counterfactual is diagnostic only because the saved V4.1 candidates were
  trained with weak controls; no V4.2 performance claim is available
- the user proposed integrating new data and repeated clustering to expand NK
  training labels; this is scientifically reasonable only as a development
  pseudo-label lane with locked whole-dataset validation and explicit approval
- the proposed development cohorts contain 318,156 T/NK candidates; the locked
  GSE169246 cohort provides 7,770 author NK and 54,925 productive alpha-beta T
  cells, and GSE315928 provides 66,813 paired-alpha-beta controls
- V4.2 completed nine global and nine boundary RAPIDS Leiden runs with no CPU
  fallback and wrote checksum-verified partitions; the expression-independent
  boundary contained 3,978,014/4,023,462 cells (`98.87%`), showing poor
  boundary localization
- six of 396 clusters passed the 95% anchor-purity and 2% productive-TCR
  contamination limits, but each had `86.92%-89.69%` of eligible candidates
  from GSE292700, above the frozen 70% source cap; zero clusters and zero cells
  qualified
- the print-safe five-page HTML/PDF QC report was visually checked; no
  classifier, calibration, threshold, release artifact, or atlas inference was
  produced

- V4.1-GPU Gate C was explicitly approved and checksum-bound; 99/99 recorded
  inner-fold fits converged and all 57 complete threshold-frontier files were
  retained
- zero of 15 Stage-1 candidate-by-outer-fold evaluations passed the frozen
  recall plus maximum-source NK-passage rule; the one-cell GSE190870 stratum
  reached 100% passage, but even among NK strata with at least 100 cells the
  best worst-source passage remained 90.97%-93.23%
- zero compact seven-gene, V2-like TCR-gene, or raw legacy-score operating
  points passed balanced or high-purity constraints across the three folds
- V4 Stage 2, extension-model false-positive rates, and paired hierarchical
  bootstrap were not estimable because Stage 1 selected no eligible candidate;
  no V4.1 model was promoted, released, or applied to the atlas

- V4 v1.2 CPU outer folds completed for HRA005041 and GSE144469; no V4,
  compact seven-gene, V2-like, or legacy-score threshold satisfied either
  frozen operating mode in either fold
- all 16 Stage-1 SAGA candidates were nonconverged in both completed CPU outer
  folds, and the selected Stage-1 threshold passed 100% of held-out strict-NK
  cells in both folds
- the BALF_BLOOD_COPD CPU fold was terminated at the user's direction before
  an atomic fold result was written; completed CPU folds and the feature cache
  remain preserved as read-only diagnostic evidence
- GPU feasibility passed on synthetic data with the A100 80 GB: the 0.835-GiB
  gene cache has ample memory headroom, deterministic PyTorch ridge logistic
  produced bit-identical repeat probabilities, and GPU XGBoost produced
  bit-identical repeat probabilities plus portable UBJ export
- the default `/tmp/nvidia-mps` daemon belongs to another user; gdTAI direct
  CUDA succeeds with an isolated `CUDA_MPS_PIPE_DIRECTORY`, and the plan
  prohibits modifying or connecting to the other user's daemon
- cuML logistic regression is excluded from the proposed candidates because a
  full-scale synthetic fit emitted a line-search stop and the installed version
  did not expose adequate convergence provenance
- protocol v1.2 Step 1 passed all frozen input, feature, truth-label,
  exclusion-feasibility, grouped-split, hash, and file-state checks
- Step 2 implementation passes 17 deterministic unit/integration tests,
  including all model families, fair comparators, threshold guardrails,
  hierarchical bootstrap, sparse extraction, and print-safe reporting
- the pinned Step 2 evaluation config, runner, and core SHA-256 values are
  `de930d4484a57a8ba77d4df79da2d2018c1b35e7095da5f6bdb6db8279082b32`,
  `9102b8ebfe617df8c218bbd34a58a5e5ac31873a938de969f74f413913d85a07`, and
  `eaf08999a6cffe2389c3deb73669b0591c6847d7fd6aa372b32c3368919fac7e`
- the first authorized cache attempt verified all source hashes and extracted
  all rows, then stopped before fitting because the cache audit compared
  all-cell exclusions against Step 1 flags that were intentionally populated
  only for positive/sensitivity and sampled-NK rows; the evaluator now requires
  exact reproduction on that frozen scope and applies the unchanged rules to
  every cell in memory
- an incomplete timing run was stopped before any outer-fold result after the
  serial fold implementation proved inefficient and emitted unrecorded SAGA
  convergence warnings; the evaluator now runs the same frozen folds in
  deterministic parallel and records iteration/convergence provenance
- the clean execution uses 26 candidate workers and three fold workers for a
  maximum of 78 concurrent fits, within the 80-CPU contract; this changes only
  scheduling, not data, models, seeds, folds, or selection
- the registered HRA005041 `log1p(CP10K)` matrix passed a full 766,639-cell
  inverse-library-size audit: zero empty rows, zero rows outside the frozen
  `1e-4` tolerance, and maximum absolute deviation
  `7.275957614183426e-12` from 10,000
- GDT_2020AUG_woCOV, MalteGDT, and GDTlung are sensitivity-only and have zero
  training weight; their reduced feature coverage cannot affect fitting
- the fixed CD4/Treg exclusions passed their feasibility guardrail with a
  source-macro primary-gdT recall ceiling of `98.38%`; those thresholds remain
  frozen and were not tuned
- no unused gdT-positive cohort remains after adaptive model iteration;
  `BALF_BLOOD_COPD` and screened extension cohorts are development/stress
  benchmarks, not independent or prospective holdouts
- V2 validation metrics and V3 BALF metrics are model-selection estimates
- the promoted V3 artifact checksum is valid, but its pickle, manifest,
  completed-run summary, and documented training composition do not describe
  one coherent build record
- extension-cohort and COPD integration is structurally complete in
  `integrated_full_atlas.h5ad`; repaired-TCR propagation and the GSE169246
  additive compartment correction remain separate approval-gated writes
- GSE169246 has a confirmed compartment-join error: all 239,418 retained `_b`
  cells are official blood libraries, while current local tissue/specimen fields
  include tumor sites or non-blood contexts; correction is not yet approved
- frozen-model screening is complete, but these cohorts contain zero gdT-gold
  and zero gdT-silver cells under the project TCR rules; they cannot support
  new-cohort recall, F1, ROC-AUC, PR-AUC, or model promotion
- the extension inputs used alpha-beta TCR assays rather than gdTCR assays;
  zero productive TRG/TRD is therefore expected, and these cohorts support
  alpha-beta T-cell and NK false-positive evaluation rather than independent
  TCR-defined gdT recall
- five extension cohorts have unmatched productive TCR-contig warnings that
  require review; all seven nevertheless passed structural QC
- full-atlas Phase 4 score computation and finite-value QC passed; biological
  prediction thresholds remain separate model-review decisions
- disease-status correction remains export-only and has not been written back
  into `integrated.h5ad`
- migration of the mirrored SSD-side `integrated.h5ad` back to NFS is not yet
  approved
- cluster annotation evidence review outputs are review-only and have not been
  written back into milestone H5AD annotation columns

## Active exceptions in force

- active RAM ceiling override: `800 GB`
- large-H5AD mirrored SSD workflow is active
- gdTAI V4.2 saved-latent clustering and consensus use an approved 150-GiB SSD
  free-space floor; prepare and fit retain 300 GiB

## Current canonical objects

- canonical candidate milestone:
  - `Integrated_dataset/TNK_candidates.h5ad`
- canonical cleaned milestone:
  - `Integrated_dataset/TNK_cleaned.h5ad`
- canonical rebuilt full-atlas milestone for downstream structural analysis:
  - `Integrated_dataset/integrated_full_atlas.h5ad`
  - physical object:
    `/ssd/tnk_phase3/Integrated_dataset/full_atlas/integrated_full_atlas.h5ad`
  - status: complete and QC-passed; repaired TCR metadata is not propagated and
    GSE169246 tissue/specimen analysis remains blocked
- legacy current integrated path:
  - `high_speed_temp/Integrated_dataset/integrated.h5ad`
  - status: absent from the frozen project path and superseded for the expanded
    structural atlas by `integrated_full_atlas.h5ad`

Additional current state:

- gdTAI V4 protocol v1.2 keeps the 197-gene universe, labels, folds, model
  families, grids, guardrails, CD4/Treg thresholds, and promotion rule fixed
- protocol v1.2 removes silver, dual/ambiguous, and all sorted cells from both
  fitting stages; no fitting result informed this amendment
- gdTAI V4 Step 1 created a checksum-pinned 1,137,739-cell label manifest and
  deterministic three-fold grouped split manifest without fitting a model
- protocol v1.2 Step 1 completed with `PASS`, zero failures, and zero warnings;
  its cell-label and grouped-split manifest SHA-256 values are
  `8157cbebfedeb84fc34eba05429aa8a7a834c6f7ceba15bd4790b5dd06bf7e0c` and
  `c84da2ca999676bab0ed180ae3d380e1e7d5b2e08da886ae2f6f912f9e0080a7`
- all 16 declared preflight inputs matched registered hashes where available
  and retained identical size and nanosecond modification time
- GSE144469 joined 107,068/107,068 raw-expression rows uniquely by the SRR and
  barcode encoded in cell IDs; its valid raw matrix is `layers/counts`
- the earlier experimental V4 bundle was archived intact under
  `archive/retired_experiments/gdTAI_v4_experimental_precommit_20260807/`
- seven extension H5ADs containing 498,901 cells and 154 libraries were built
  under `data/interim/extension_intake/built_h5ads/`
- all seven extension H5ADs passed sparse raw-count, metadata, TCR-schema,
  logical-flag, and sample-plus-barcode Phase 0 structural QC
- extension tumor-project specimen context initially passed the rule-based
  structural gate for all 498,901 cells with zero blank/format violations;
  later official-GEO review exposed a semantic compartment error in GSE169246,
  and original tissue and diagnosis fields remain unmodified
- all eight new inputs were independently filtered to `758,135` T/NK
  candidates from `893,126` cells; all `338,046` productive-TCR cells were
  retained
- all eight filtered outputs passed required-metadata completeness, exact
  source-accession, cohort-ID, unique cell/library-barcode, canonical TCR-flag,
  and specimen-context checks
- GSE169246 metadata now distinguish 51 biological samples from 78 source
  libraries using the full barcode suffix; 76 libraries retain T/NK cells
- source H5AD checksums were unchanged and no source H5AD was modified
- all four frozen gdTAI profiles were applied read-only to all 758,135
  extension T/NK cells using raw-count CP10K plus log1p normalization
- V3 Round 12 high-purity had the lowest pooled alpha-beta TCR FPR (`0.227%`),
  while V2 high-purity had the lowest strict-NK FPR (`0.725%`)
- after excluding reduced-feature GSE169246 as a sensitivity analysis, V3
  Round 12 high-purity had the lowest known-negative union FPR (`0.295%`)
- no model was selected or promoted; historical positive sensitivity still
  favors promoted V3 Round 14 over Round 12
- all eight evaluation-source H5AD SHA-256 values remained unchanged
- official-GEO reconciliation is complete for 58 grouped assertions across 8
  cohorts and 10 accessions: 16 locally verified, 32 GEO-resolved but locally
  unresolved, 6 partially ambiguous, and 4 unavailable in GEO
- the GEO audit passed 10/10 accession-coverage checks and 8/8 read-only H5AD
  file-state checks; it did not write metadata or mutate an H5AD
- GSE169246 `_b` correction must use additive derived columns and retain the
  original fields and full `_b`/`_t` compartment suffix
- Tan et al. 2021 remains excluded as a duplicate of
  `GDT_2020AUG_woCOV`
- the project reorganization is validated: active code is grouped under
  `workflows/`, configuration under `configs/`, and physical pre-integration
  data are dataset-centered under `data/datasets/<dataset_id>/`
- the base dataset-centered storage migration completed with 191
  same-filesystem rename operations
- the `BALF_BLOOD_COPD` cross-study benchmark was subsequently registered; its
  storage was refined so the approximately 949 GB raw/interim
  workspace remains at its original location and only the 2,110,825,599-byte
  validation H5AD is physically stored in the project
- the current storage view has 66 dataset directories, 45 promoted
  `current.h5ad` links, and 264 lifecycle compatibility links
- the historical `downloads`, `analysis_26GSE_V4`, and `newdata` paths remain
  usable as compatibility aliases to the same inodes
- `/home/tanlikai/databank/owndata/singlecell` is again the physical
  `BALF_BLOOD_COPD` raw/interim workspace; the project raw path is a
  compatibility link to it
- the historical external H5AD path is a compatibility link to
  `data/datasets/BALF_BLOOD_COPD/processed/artifacts/phase4_final_annotated.h5ad`
- the dataset registry currently contains 67 datasets and 1,790 library rows;
  all 33 Phase 0 active datasets pass strict path and library validation
- `BALF_BLOOD_COPD` retains the legacy registry role string
  `gdTAI_independent_external_validation`, but it was reused in Round 12 versus
  Round 14 promotion and is scientifically classified as a reused cross-study
  benchmark; it remains inactive for Phase 0, current-milestone integration,
  and extended-atlas integration
- `obs["tissue_corrected"]` is already written into the current integrated
  milestone
- `HRA005041` intake H5AD was rebuilt with harmonized metadata plus productive
  sample-aware `TRA/TRB/TRG/TRD` fields under
  `downloads/per_gse_h5ad_with_metadata/HRA005041_T_cells_subset.h5ad`
- `HRA005041` now also carries standalone Phase 4 score columns from the
  project-standard scoring method
- three sorted-gdT Seurat RDS datasets were converted into harmonized standalone
  H5ADs under `newdata/Sorted_gdT/`
- those sorted-gdT intake H5ADs now carry canonical `TRG/TRD` metadata,
  `Sorted_gdT = True`, recomputed Scanpy UMAPs, and standalone Phase 4 score
  columns
- repaired TCR propagation is complete for ten approved GSEs
- the fourteen structurally unscoreable no-TCR-gene GSEs have been removed from
  the current milestone H5ADs
- removed cells were preserved in
  `high_speed_temp/Integrated_dataset/No_TCR_Gene_GSEs.h5ad`
- the current integrated milestone now contains `3,705,306` cells after that
  carve-out
- scANVI fields remain reference-only; simple scVI-based annotation is the
  canonical downstream interpretation layer
- gdTAI v3 Round 12 and Round 14 were revalidated on identical cohorts with
  checksum-pinned artifacts; Round 14 at threshold `0.936` is the promoted
  balanced default, while Round 12 at threshold `0.5` is preserved as the
  validated high-purity fallback; the 2026-08-06 audit found that their
  corrected TCR-only F1 values are effectively tied and that the reused
  benchmark cannot establish Round 14 superiority
- `GSE305372` was excluded and retired on 2026-08-06; it has no active registry
  rows, its local scientific data were deleted, and its code, reports, source
  checksums, and historical results were moved under `archive/`

## Current review artifacts

- rebuilt full T/NK atlas with extension cohorts and COPD BALF/BLOOD:
  - `full_atlas_rebuild/index.html`
  - `full_atlas_rebuild/full_atlas_rebuild_report.pdf`
  - `Integrated_dataset/integrated_full_atlas.h5ad`
  - `Integrated_dataset/logs/full_atlas_rebuild/`
  - `Integrated_dataset/tables/full_atlas_rebuild/`
  - `Integrated_dataset/figures/full_atlas_rebuild/`
  - decision: structural and computational QC passed; repaired-TCR propagation
    and GSE169246 additive compartment correction remain separate gated tasks

- gdTAI V4.2 productive-TRA/TRB boundary conflict audit:
  - `gdT_prediction/gdtai_v4_2_tcr_conflict_audit/index.html`
  - `gdT_prediction/gdtai_v4_2_tcr_conflict_audit/gdtai_v4_2_tcr_conflict_audit_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_tcr_conflict_audit/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_tcr_conflict_audit/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_tcr_conflict_audit/`
  - decision: `PASS_JOIN_REPAIR_REQUIRED`; do not use raw harmonized
    productive TRA/TRB as boundary truth until affected source joins are
    rebuilt and validated

- gdTAI V4.2 NK/T-lineage feature and unsupervised-cluster review:
  - `gdT_prediction/gdtai_v4_2_nk_cluster_review/index.html`
  - `gdT_prediction/gdtai_v4_2_nk_cluster_review/gdtai_v4_2_nk_cluster_review_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_nk_cluster_review/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_nk_cluster_review/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_nk_cluster_review/`
  - decision: `PASS_VISUAL_REVIEW_EXTENDED`; all 461 cluster-19 rows are held
    as a provisional core pending source/library TCR and doublet review, and all
    remaining candidate tiers stay review-only

- gdTAI V4.2 conservative NK-definition repair:
  - `gdT_prediction/gdtai_v4_2_nk_definition_repair/index.html`
  - `gdT_prediction/gdtai_v4_2_nk_definition_repair/gdtai_v4_2_nk_definition_repair_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_nk_definition_repair/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_nk_definition_repair/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_nk_definition_repair/`
  - decision: `PASS_NK_REFERENCE_READY`; the preserved 469-row audit manifest
    was reduced to a provisional 461-cell review core and not accepted as
    training truth after the later raw TCR audit; no classifier was fitted

- gdTAI V4.2 sidecar visualization and confounding diagnostics:
  - `gdT_prediction/gdtai_v4_2_nk_reference_diagnostics/index.html`
  - `gdT_prediction/gdtai_v4_2_nk_reference_diagnostics/gdtai_v4_2_nk_reference_diagnostics_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_nk_reference_diagnostics/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_nk_reference_diagnostics/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_nk_reference_diagnostics/`
  - decision: `PASS_DIAGNOSTIC_ONLY`; anchor/cohort confounding blocks confident
    NK propagation and no classifier was fitted

- gdTAI V4.2 integration, clustering, and pseudo-NK consensus QC:
  - `gdT_prediction/gdtai_v4_2_nk_reference/index.html`
  - `gdT_prediction/gdtai_v4_2_nk_reference/gdtai_v4_2_nk_reference_qc_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_nk_reference/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_nk_reference/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_nk_reference/`
  - decision: technical `PASS`, scientific `FAIL_NO_PSEUDO_NK`; no classifier
    fitting from this lane

- gdTAI V4.2 modeling-integration preflight:
  - `configs/models/gdtai/v4_2_cohort_roles.csv`
  - `configs/models/gdtai/v4_2_integration_preflight.json`
  - `gdT_prediction/gdtai_v4_2_integration_preflight/index.html`
  - `gdT_prediction/gdtai_v4_2_integration_preflight/gdtai_v4_2_integration_preflight_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_integration_preflight/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_integration_preflight/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_integration_preflight/`

- gdTAI V4.2 Step 0 NK-label audit and proposed repair:
  - `docs/GDTAI_V4_2_PRECOMMITTED_PLAN.md`
  - `gdT_prediction/gdtai_v4_2_step0/index.html`
  - `gdT_prediction/gdtai_v4_2_step0/gdtai_v4_2_step0_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_2_step0/`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_2_step0/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_2_step0/`

- gdTAI V4.1-GPU pre-fit feasibility and protocol:
  - `docs/GDTAI_V4_GPU_PRECOMMITTED_PLAN.md`
  - `gdT_prediction/gdtai_v4_gpu_precommit/index.html`
  - `gdT_prediction/gdtai_v4_gpu_precommit/gdtai_v4_gpu_precommitted_plan.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_gpu_precommit/gpu_feasibility_summary.md`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_gpu_precommit/gpu_feasibility_checks.csv`
- gdTAI V4 Step 1 preflight and split freeze:
  - `gdT_prediction/gdtai_v4_preflight/index.html`
  - `gdT_prediction/gdtai_v4_preflight/gdtai_v4_preflight_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v4_preflight/gdtai_v4_preflight_summary.md`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_v4_preflight/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_v4_preflight/`
- precommitted gdTAI V4 training and validation protocol:
  - `docs/GDTAI_V4_PRECOMMITTED_PLAN.md`
  - `gdT_prediction/gdtai_v4_precommit/index.html`
  - `gdT_prediction/gdtai_v4_precommit/gdtai_v4_precommitted_plan.pdf`
- independent gdTAI training and evaluation audit:
  - `docs/GDTAI_METHODOLOGY_AUDIT.md`
  - `gdT_prediction/gdtai_methodology_audit/index.html`
  - `gdT_prediction/gdtai_methodology_audit/gdtai_methodology_audit_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_methodology_audit/gdtai_methodology_audit_summary.md`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_methodology_audit/`
- gdTAI external methodological review:
  - `Integrated_dataset/logs/gdT_prediction/external_review/gdtai_external_review_report.md`
  - `gdT_prediction/external_review/gdtai_external_review_report.html`
  - `gdT_prediction/external_review/gdtai_external_review_report.pdf`
- gdTAI-kimi experimental nested cross-dataset evaluation (not promotable):
  - `Integrated_dataset/logs/gdT_prediction/gdtai_kimi/gdtai_kimi_report.md`
  - `gdT_prediction/gdtai_kimi/gdtai_kimi_report.pdf`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_kimi/`
- extension official-GEO metadata reconciliation:
  - `docs/EXTENSION_GEO_METADATA_RECONCILIATION.md`
  - `configs/metadata/extension_geo_metadata_reconciliation.csv`
  - `Integrated_dataset/logs/geo_metadata_reconciliation/extension_geo_metadata_reconciliation_audit.md`
  - `Integrated_dataset/tables/geo_metadata_reconciliation/`
- extension frozen gdTAI profile negative-control evaluation:
  - `gdT_prediction/gdtai_extension_evaluation/index.html`
  - `gdT_prediction/gdtai_extension_evaluation/gdtai_extension_evaluation_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_extension_evaluation/gdtai_extension_evaluation_summary.md`
  - `Integrated_dataset/tables/gdT_prediction/gdtai_extension_evaluation/`
  - `Integrated_dataset/figures/gdT_prediction/gdtai_extension_evaluation/`
- extension standalone Phase 1 T/NK filter review:
  - `Integrated_dataset/logs/extension_intake/extension_tnk_filter_qc.md`
  - `Integrated_dataset/tables/extension_intake/tnk_filter/extension_tnk_filter_cohort_summary.csv`
  - `Integrated_dataset/tables/extension_intake/tnk_filter/extension_tnk_filter_metadata_audit.csv`
  - `Integrated_dataset/figures/extension_intake/tnk_filter/extension_tnk_filter_overview.png`
  - `data/interim/extension_intake/tnk_filtered_h5ads_manifest.csv`
- extension standalone Phase 0 review:
  - `Integrated_dataset/logs/extension_intake/extension_phase0_supervision_report.md`
  - `Integrated_dataset/logs/extension_intake/extension_phase0_qc_summary.md`
  - `Integrated_dataset/tables/extension_intake/extension_phase0_cohort_summary.csv`
  - `Integrated_dataset/tables/extension_intake/extension_tcr_source_chain_audit.csv`
  - `Integrated_dataset/figures/extension_intake/extension_phase0_cohort_overview.png`
- GSE169246 selected processed input:
  - `docs/GSE169246_PROCESSED_H5AD_INTAKE.md`
  - `data/datasets/GSE169246/processed/current.h5ad`
  - `data/registry/snapshots/post_gse169246_processed_h5ad_intake_20260806/`
- extension tumor-project specimen-context review:
  - `Integrated_dataset/logs/sample_source_refinement/extension_sample_source_review.md`
  - `Integrated_dataset/tables/sample_source_refinement/extensions/extension_tumor_project_context_counts.csv`
- project reorganization:
  - `docs/REORGANIZATION_RECORD.md`
  - `Integrated_dataset/logs/project_reorganization/checkpoints/pre_reorganization_20260715/`
- dataset-centered storage migration:
  - `docs/DATASET_STORAGE_MIGRATION.md`
  - `data/registry/migrations/dataset_centered_20260716/validation.json`
  - `data/registry/snapshots/post_dataset_centered_migration_20260716/`
  - `Integrated_dataset/logs/project_reorganization/checkpoints/pre_dataset_centered_migration_20260716/`
- `BALF_BLOOD_COPD` cross-study benchmark intake:
  - `docs/BALF_BLOOD_COPD_INTAKE.md`
  - `data/registry/migrations/balf_blood_copd_20260716/validation.json`
  - `data/registry/snapshots/post_balf_blood_copd_intake_20260716/`
  - `data/registry/migrations/balf_blood_copd_h5ad_only_20260716/validation.json`
  - `data/registry/snapshots/post_balf_blood_copd_h5ad_only_20260716/`
- gdTAI v3 Round 12 versus Round 14 model decision:
  - `gdT_prediction/gdtai_v3_round12_vs_round14/index.html`
  - `gdT_prediction/gdtai_v3_round12_vs_round14/gdtai_v3_round12_vs_round14_report.pdf`
  - `Integrated_dataset/logs/gdT_prediction/gdtai_v3_round12_vs_round14/gdtai_v3_round12_vs_round14_summary.json`
- Phase 4 QC:
  - `Integrated_dataset/logs/phase4_qc_summary.md`
- tissue correction:
  - `Integrated_dataset/tables/tissue_correction/tissue_correction_apply.md`
- repaired TCR review:
  - `Integrated_dataset/logs/tcr_rebuild_phase4_qc.md`
- HRA005041 intake review:
  - `Integrated_dataset/logs/HRA005041_tcr_intake.md`
- HRA005041 standalone Phase 4 review:
  - `Integrated_dataset/logs/HRA005041_phase4.log`
- Sorted_gdT intake review:
  - `Integrated_dataset/logs/Sorted_gdT/sorted_gdt_dataset_summary.md`
- no-TCR-gene GSE carve-out:
  - `Integrated_dataset/logs/no_tcr_gene_gse_removal.md`
- disease-status export review:
  - `Integrated_dataset/tables/disease_status_correction/disease_status_corrected_column_export.csv.gz`
- cluster annotation evidence review:
  - `Integrated_dataset/logs/annotation_evidence_review/current_integrated_annotation_evidence_summary.md`
  - `Integrated_dataset/logs/annotation_evidence_review/plus6_annotation_evidence_summary.md`
- GSE305372 retirement record:
  - `archive/retired_experiments/GSE305372_external_application/README.md`
  - `archive/retired_experiments/GSE305372_external_application/deleted_source_files.csv`
  - `maintenance/reorganization/gse305372_retirement_map.csv`
