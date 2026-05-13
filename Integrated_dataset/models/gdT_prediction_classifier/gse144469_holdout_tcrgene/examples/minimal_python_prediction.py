#!/usr/bin/env python
"""Minimal programmatic example for applying the selected gdT classifier.

This example intentionally writes only a compact prediction table. For full
tables, summaries, figures, and HTML output, use predict_with_selected_gdt_model.py
from the repository root.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[5]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from predict_with_selected_gdt_model import (  # noqa: E402
    DEFAULT_MODEL,
    apply_death_penalty,
    build_gene_mapping,
    dense_chunk_features,
    extract_csr_features,
    infer_matrix_encoding,
    load_model,
    read_obs_index,
    sha256_file,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Minimal selected gdT model prediction example.")
    parser.add_argument("--input-h5ad", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--model-pkl", type=Path, default=REPO_ROOT / DEFAULT_MODEL)
    parser.add_argument("--chunk-size", type=int, default=50_000)
    parser.add_argument("--max-cells", type=int, default=None, help="Optional small smoke-test limit.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    payload = load_model(args.model_pkl)
    model = payload["model_object"]
    threshold = float(payload["threshold"])
    gene_names = [str(gene) for gene in payload["gene_names"]]
    model_sha = sha256_file(args.model_pkl)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(args.input_h5ad, "r") as handle:
        matrix_encoding = infer_matrix_encoding(handle)
        if matrix_encoding not in {"csr_matrix", "csr_matrix_or_csc_matrix", "dense"}:
            raise TypeError(f"Unsupported X encoding: {matrix_encoding}")

        obs_index = read_obs_index(handle)
        n_obs = int(obs_index.shape[0])
        if args.max_cells is not None:
            n_obs = min(n_obs, int(args.max_cells))

        gene_idx_to_col, gene_availability = build_gene_mapping(handle, gene_names)
        missing_genes = int((~gene_availability["available_in_h5ad"]).sum())
        extract_fn = dense_chunk_features if matrix_encoding == "dense" else extract_csr_features

        n_positive = 0
        with gzip.open(args.output_csv, "wt", encoding="utf-8", newline="") as out:
            for start in range(0, n_obs, args.chunk_size):
                end = min(start + args.chunk_size, n_obs)
                features, row_sum, n_detected = extract_fn(
                    handle,
                    start,
                    end,
                    gene_idx_to_col,
                    len(gene_names),
                )
                raw_score = model.predict_proba(features)[:, 1].astype(np.float32)
                score, penalty = apply_death_penalty(raw_score, features, gene_names)
                pred = score >= threshold
                n_positive += int(pred.sum())

                chunk = pd.DataFrame(
                    {
                        "cell_id": obs_index[start:end],
                        "gdt_score": score,
                        "gdt_score_before_penalty": raw_score,
                        "predicted_gdT": pred,
                        "death_penalty_applied": penalty,
                        "row_sum_x": row_sum,
                        "n_detected_genes_x": n_detected,
                        "missing_model_gene_count": missing_genes,
                    }
                )
                chunk.to_csv(out, index=False, header=(start == 0))

    print(f"model_sha256={model_sha}")
    print(f"threshold={threshold:.16g}")
    print(f"cells_scored={n_obs}")
    print(f"predicted_gdT_cells={n_positive}")
    print(f"predicted_gdT_fraction={n_positive / n_obs if n_obs else 0:.6g}")
    print(f"missing_model_genes={missing_genes}")
    print(f"output_csv={args.output_csv}")


if __name__ == "__main__":
    main()

