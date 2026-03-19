# Phase 0 QC Summary

## Environment

- Requested env from runbook: `Scanpy_gdTmodel`
- Resolved working env for this run: `/home/tanlikai/miniconda3/envs/rapids_sc_py310`
- GPU detected by `nvidia-smi`: NVIDIA A100 80GB PCIe, CUDA 12.9

## Audit Scope

- Registry rows audited: 33
- Files readable without recorded errors: 33
- Files with read issues: 0

## Category Counts

- Category A: 13
- Category B: 0
- Category C: 20

## Matrix State Counts

- normalized_or_noninteger_X: 20
- raw_like_X: 13

## Metadata Coverage

- Datasets with candidate donor fields: 31
- Datasets with candidate sample fields: 32
- Datasets with candidate library fields: 30

## Category B Datasets

- None

## Category C Datasets

- GSE139555: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source
- GSE145926: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source
- GSE161918: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source
- GSE162498: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source
- GSE171037: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source
- GSE178882: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source
- GSE188620: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source
- GSE190870: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source
- GSE211504: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source
- GSE212217: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source
- GSE227709: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source
- GSE243572: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source
- GSE243905: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source
- GSE252762: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source
- GSE254176: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source
- GSE254249: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source
- GSE267645: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source
- GSE301528: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source
- GSE308075: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source
- GSE311112: normalized_or_noninteger_X | adata.X is non-integer without a recoverable raw-count source

## QC Conclusion

Phase 0 audit outputs are generated and ready for user review. Do not advance to Phase 1 until the user reviews the tables and figures and explicitly approves the phase transition.
