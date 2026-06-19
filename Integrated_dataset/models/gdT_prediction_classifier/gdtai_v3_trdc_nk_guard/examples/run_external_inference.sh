#!/usr/bin/env bash
set -euo pipefail
/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python \
  "$(dirname "$0")/predict_h5ad_counts.py" \
  --input-h5ad "$1" \
  --output-csv "$2"
