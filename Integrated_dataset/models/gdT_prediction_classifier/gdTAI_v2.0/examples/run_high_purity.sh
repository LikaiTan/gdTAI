#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 /path/to/test_dataset.h5ad DATASET_ID [annotation_column]" >&2
  exit 2
fi
ANNOTATION_COLUMN="${3:-simple_annotation_plus6}"
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python predict_with_gdtai_v2.py   --input-h5ad "$1"   --dataset-id "$2"   --mode high_purity   --annotation-column "$ANNOTATION_COLUMN"   --chunk-size 50000
