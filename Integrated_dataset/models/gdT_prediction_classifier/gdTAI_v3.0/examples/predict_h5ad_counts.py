#!/usr/bin/env python
"""Apply the packaged gdTAI v3 TRDC/NK guard candidate to an H5AD counts layer."""
from __future__ import annotations

import argparse
import os
import pickle
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[5]
WORKFLOW_DIR = PROJECT_ROOT / "workflows" / "gdtai"
os.environ.setdefault("MPLCONFIGDIR", "/tmp/gdtai_matplotlib_cache")
if str(WORKFLOW_DIR) not in sys.path:
    sys.path.insert(0, str(WORKFLOW_DIR))

from run_gdtai_v3_trdc_nk_guard_classifier import (
    FeatureSpec,
    append_engineered_features,
    extract_gene_features,
    h5ad_shape,
    obs_index_values,
    read_string_dataset,
    selected_gene_mapping,
    trdc_trdv_quadrant,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run gdTAI v3 TRDC/NK guard inference on layers['counts'].")
    parser.add_argument("--input-h5ad", type=Path, required=True)
    parser.add_argument("--model-pkl", type=Path, default=Path(__file__).resolve().parents[1] / "gdTAI_v3_model.pkl")
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=50000)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def build_input_spec(var_names: list[str], gene_names: list[str], engineered_feature_names: list[str], model_feature_names: list[str]) -> FeatureSpec:
    availability = selected_gene_mapping(var_names, gene_names)
    missing = availability.loc[~availability["available_in_h5ad"], "gene"].astype(str).tolist()
    if missing:
        raise KeyError(f"Input H5AD is missing {len(missing)} model genes; first missing genes: {missing[:20]}")
    lookup = {gene: idx for idx, gene in enumerate(var_names)}
    gene_to_col = {gene: i for i, gene in enumerate(gene_names)}
    engineered_to_col = {name: len(gene_names) + i for i, name in enumerate(engineered_feature_names)}
    return FeatureSpec(
        gene_names=gene_names,
        gene_indices=np.asarray([lookup[gene] for gene in gene_names], dtype=np.int32),
        gene_feature_names=[f"{gene}_log1p_cp10k" for gene in gene_names],
        engineered_feature_names=engineered_feature_names,
        model_feature_names=model_feature_names,
        gene_to_col=gene_to_col,
        engineered_to_col=engineered_to_col,
    )


def main() -> None:
    args = parse_args()
    if args.output_csv.exists() and not args.overwrite:
        raise FileExistsError(f"Output already exists: {args.output_csv}. Pass --overwrite to replace it.")
    if args.output_csv.exists():
        args.output_csv.unlink()
    with args.model_pkl.open("rb") as handle:
        payload = pickle.load(handle)
    threshold = float(payload["threshold"])
    model = payload["model_object"]
    gene_names = [str(gene) for gene in payload["gene_names"]]
    engineered_feature_names = [str(name) for name in payload["engineered_feature_names"]]
    model_feature_names = [str(name) for name in payload["feature_names"]]

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    wrote_header = False
    with h5py.File(args.input_h5ad, "r") as h5ad:
        if "layers" not in h5ad or "counts" not in h5ad["layers"]:
            raise RuntimeError("Input H5AD must contain layers['counts']; refusing to use normalized/log X.")
        var_names = pd.Index(read_string_dataset(h5ad["var"]["_index"]), dtype="string").astype(str).tolist()
        spec = build_input_spec(var_names, gene_names, engineered_feature_names, model_feature_names)
        n_obs = h5ad_shape(h5ad, "layers/counts")[0]
        cell_ids = obs_index_values(h5ad)
        for start in range(0, n_obs, args.chunk_size):
            end = min(start + args.chunk_size, n_obs)
            rows = np.arange(start, end, dtype=np.int64)
            x_gene, row_sum, n_detected = extract_gene_features(h5ad, "layers/counts", rows, spec, label=f"gdtai_v3_{start}_{end}")
            x = append_engineered_features(x_gene, spec)
            score = model.predict_proba(x)[:, 1].astype(np.float32)
            out = pd.DataFrame(
                {
                    "cell_id": cell_ids[rows],
                    "gdtai_v3_score": score,
                    "gdtai_v3_threshold": threshold,
                    "gdtai_v3_predicted": score >= threshold,
                    "tcr_gene_quadrant": trdc_trdv_quadrant(x, spec),
                    "row_sum_counts_layer": row_sum,
                    "n_detected_genes_counts_layer": n_detected,
                }
            )
            out.to_csv(args.output_csv, index=False, mode="a", header=not wrote_header, compression="infer")
            wrote_header = True


if __name__ == "__main__":
    main()
