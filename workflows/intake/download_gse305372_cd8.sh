#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUTPUT_DIR="${PROJECT_ROOT}/data/datasets/GSE305372/raw/source"
BASE_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE305nnn/GSE305372/suppl"

mkdir -p "${OUTPUT_DIR}"

download_one() {
  local name="$1"
  local expected_size="$2"
  local expected_sha256="$3"
  local target="${OUTPUT_DIR}/${name}"

  if [[ -f "${target}" ]] && [[ "$(stat -c '%s' "${target}")" == "${expected_size}" ]]; then
    printf '%s  %s\n' "${expected_sha256}" "${target}" | sha256sum --check -
    return
  fi

  curl \
    --fail \
    --location \
    --continue-at - \
    --retry 10 \
    --retry-all-errors \
    --output "${target}" \
    "${BASE_URL}/${name}"

  [[ "$(stat -c '%s' "${target}")" == "${expected_size}" ]]
  printf '%s  %s\n' "${expected_sha256}" "${target}" | sha256sum --check -
}

download_one \
  "GSE305372_HIPC-1_LG-CD8-ALL_SeuratObject.RDS" \
  "5623543662" \
  "b4e1f62c4f34fc6b024f3769dbffd0e9513cf7d22bab915e0dd1a1e3eb4fca71"

download_one \
  "GSE305372_HIPC-1_LN-CD8-ALL_SeuratObject.RDS" \
  "515670151" \
  "2e36d7ff4ab771e400a13cb807431198aa8ec8f7b7d1554842cb49a1196958e7"
