# DECISIONS.md

## Current approved decisions

- gdTAI V4.2 uses `NK_CONFIDENT` only for development cells that have no
  productive TCR or doublet flag and pass both independent evidence classes:
  exact latent nearest-anchor evidence and strict NK-lineage expression
- the strict expression rule requires no detected
  `CD3D/CD3E/CD3G/TRAT1/BCL11B`, both `FCER1G` and `TYROBP`, and at least one of
  `KLRD1/NCR1/FCGR3A/KLRC1`; `NKG7/GNLY/PRF1/GZMB/XCL1/XCL2` remain audit-only
  because cytotoxicity is shared by NK, CD8, and gamma-delta T cells
- approximate cuML/FAISS IVF nearest-neighbor output is rejected for this label
  definition after repeated runs changed the accepted count; exact A100
  brute-force neighbors with an internal bit-identical repeat are required
- all accepted V4.2 pseudo-NK rows are retained for possible future training at
  low effective weight; per-source effective contribution is capped at 70%
- the 2026-08-17 NK-definition repair is a reference-label milestone, not a
  trained classifier, promotion decision, release artifact, or atlas inference
- gdTAI V4.2 recovered development-data sparse staging and A100 scVI were
  approved and completed on 2026-08-14 against the checksum-bound recovery
  package; on 2026-08-16 the user approved the checksum-bound 150-GiB SSD floor
  for saved-latent clustering and pseudo-NK consensus QC only
- routine, reversible analysis phases proceed after documented automated QC;
  explicit approval is reserved for high-risk actions such as destructive
  deletion, source-data mutation, irreversible model promotion or release, and
  material infrastructure or external-cost risk
- gdTAI V4.2 saved-latent clustering and consensus completed on 2026-08-16,
  but the pseudo-label lane is rejected: zero clusters and zero of 113,287
  eligible cells met the frozen source-balanced consensus contract, so no
  classifier may be fitted from this lane
- the follow-up 250,000-cell sidecar diagnostic found complete anchor/cohort
  confounding: primary NK anchors exist only in the existing atlas and
  productive-TCR anchors only in new cohorts; current NK-anchor enrichment is
  therefore insufficient to assign new-cohort NK pseudo-labels
- the frozen 70% maximum-source criterion will not be relaxed in response to
  this outcome; the next iteration must repair candidate source balance and
  boundary localization while retaining the locked cohorts unchanged and
  establishing within-cohort NK and T anchors before propagation
- gdTAI V4 must not be pushed to GitHub until it is scientifically finished
  and reviewed; local commits are allowed for rollback and reproducibility
- gdTAI V4.2 sidecar-integration implementation and synthetic/read-only
  implementation QC are approved against the 2026-08-08 checksum-bound
  preflight; project-data scVI fitting, classifier fitting, promotion, release
  fitting, and atlas inference are not authorized by this approval
- current canonical execution environment is `rapids_sc_py310`
- Phase 5 model training and validation are out of scope unless explicitly requested
- scANVI outputs remain in the integrated object for reference only
- simple scVI-based annotation is the canonical downstream interpretation layer
- Phase 4 uses package-faithful continuous TRA/TRB/TRD scoring
- `obs["tissue_corrected"]` is canonical in the current integrated milestone
- large H5AD files may remain on the mirrored SSD tree until explicit migration approval
- tables, PNG figures, logs, scripts, and model artifacts stay on NFS
- canonical pre-integration physical storage is
  `data/datasets/<dataset_id>/`; `data/raw`, `data/interim`, and
  `data/processed` are lifecycle compatibility views
- historical `downloads`, `analysis_26GSE_V4`, and `newdata` paths remain
  supported aliases and must not be removed while legacy workflows exist
- physical storage changes use journaled same-filesystem renames with
  inode/hash validation and a tested reverse rollback plan
- dataset removal means registry deactivation, not file deletion
- `configs/models/gdtai/model_registry.csv` is authoritative for model release
  status; a directory or version-like filename alone does not imply promotion
- gdTAI v3 Round 14 at threshold `0.936` is the promoted balanced default after
  checksum-pinned revalidation against Round 12 on 2026-07-15; the 2026-08-06
  audit established that the reused benchmark influenced promotion, so this
  status is an operational release choice rather than independent proof of
  superiority
- gdTAI v3 Round 12 at threshold `0.5` is retained as the validated high-purity
  fallback and is not the canonical default
- extension cohorts with alpha-beta-only TCR assays may evaluate frozen-model
  false positives but may not drive recall/F1/AUC claims or model promotion
  when they contain no unbiased gdT positive truth set
- productive single TRA or single TRB evidence is accepted as alpha-beta
  negative evidence when productive TRG/TRD evidence is absent; paired and
  single-chain controls remain separately reported
- no-productive-TCR extension predictions are candidates of unknown truth and
  must not be counted as true positives
- the BALF/PBMC COPD cohort is registered as `BALF_BLOOD_COPD`; its canonical H5AD is
  `data/datasets/BALF_BLOOD_COPD/processed/current.h5ad`
- `BALF_BLOOD_COPD` was reused in Round 12 versus Round 14 promotion and must
  be described as a cross-study development benchmark, not an independent
  external test
- `BALF_BLOOD_COPD` is benchmark-only and must remain inactive for Phase 0,
  milestone integration, and extended-atlas integration unless the user
  separately approves a role change
- `BALF_BLOOD_COPD` raw and interim files remain physically under
  `/home/tanlikai/databank/owndata/singlecell`; only the selected validation
  H5AD is physically retained inside this project
- the historical external H5AD path remains a supported compatibility link to
  the project-managed artifact so the cohort's original scripts and RStudio
  project paths continue to work
- future gdTAI selection requires a precommitted nested leave-one-dataset-out
  design with donor/library/clonotype grouping, expression-independent primary
  labels, dataset-balanced fitting, fold-local feature selection/calibration,
  and dataset-macro plus worst-dataset metrics
- any cohort inspected for algorithm, threshold, guardrail, or promotion
  decisions is a development benchmark and cannot later be called independent
  or sealed
- gdTAI inference must abstain or clearly warn when critical genes or adequate
  model-gene coverage are missing
- tumor-related projects must preserve disease-aware specimen context that
  distinguishes primary tumor, adjacent/non-tumor tissue, blood, and metastatic
  sites whenever source metadata provide that distinction; original metadata
  fields are not overwritten
- official-GEO metadata reconciliation is evidence-only until separately
  approved; supported corrections must be written to additive derived columns,
  preserve original source values, and retain a provenance link to the exact
  GEO record
- GSE169246 library compartment suffixes are biologically meaningful: `_b`
  denotes blood and `_t` denotes tumor. Any corrected specimen key must retain
  this suffix, and the current `_b` tissue/specimen mismatch must be resolved
  before extension merge or integration
- the seven extension H5ADs remain standalone and unmerged after their Phase 0
  structural review until explicit user approval
- GSE169246 remains fail-closed because the available local files are not
  proven accession-pure, and Tan et al. 2021 is excluded as a duplicate of
  `GDT_2020AUG_woCOV`
- repaired TCR propagation is approved for these ten GSEs:
  - `GSE188620`
  - `GSE212217`
  - `GSE243572`
  - `GSE243905`
  - `GSE254249`
  - `GSE308075`
  - `GSE311112`
  - `GSE161918`
  - `GSE168859`
  - `GSE227709`

## Historical exceptions and notes

- on 2026-08-14, recovered V4.2 staging and A100 scVI completed successfully;
  RAPIDS clustering failed closed before graph construction when shared SSD
  capacity fell below the frozen 300-GiB floor
- an unrelated BAM sort subsequently completed and self-cleaned; the user
  approved the operational-only contract retaining 300 GiB for staging/fitting
  and 150 GiB for clustering/consensus on 2026-08-16
- on 2026-08-14, the user approved V4.2 development-data sidecar execution
  against the original SSD input contract, but the command stopped before
  matrix reading because the SSD-resident `integrated.h5ad` was absent
- that approval was automatically invalidated when the execution config,
  runner, core, and physical input contract changed to the audited
  `TNK_cleaned.h5ad` recovery path; recovery execution requires a new explicit
  checksum-bound approval
- on 2026-08-06, the user explicitly approved permanent local deletion of
  `GSE305372` as a one-time exception to the default registry-deactivation-only
  removal rule; source URLs, byte sizes, SHA-256 checksums, archived code, and
  Git rollback information were retained before deletion
- a historical run-specific exception allowed Phase 2 internal QC to continue
  directly into Phase 3 on 2026-03-20
- that exception is historical only and is not currently active
- the active RAM ceiling override for this run is `800 GB`
