#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: bash apply_model_to_h5ad.sh /path/to/input.h5ad [DATASET_ID]" >&2
  exit 2
fi

INPUT_H5AD="$1"
DATASET_ID="${2:-$(basename "${INPUT_H5AD}" .h5ad)}"

REPO_DIR="${REPO_DIR:-/home/tanlikai/databank/publicdata/tools/output_geo_tcell_research}"
PYTHON_BIN="${PYTHON_BIN:-/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python}"
CHUNK_SIZE="${CHUNK_SIZE:-50000}"

cd "${REPO_DIR}"

"${PYTHON_BIN}" predict_with_selected_gdt_model.py \
  --input-h5ad "${INPUT_H5AD}" \
  --model-pkl Integrated_dataset/models/gdT_prediction_classifier/gse144469_holdout_tcrgene/selected_model.pkl \
  --dataset-id "${DATASET_ID}" \
  --chunk-size "${CHUNK_SIZE}" \
  --obs-column has_TRA_TRB_paired \
  --obs-column has_any_ab_tcr \
  --obs-column annotation \
  --obs-column cell_type

echo "Prediction table:"
echo "Integrated_dataset/tables/gdT_prediction/external_tests/${DATASET_ID}/gdt_predictions.csv.gz"
echo "Static report:"
echo "gdT_prediction/external_tests/${DATASET_ID}/index.html"

