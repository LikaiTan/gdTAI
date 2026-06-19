#!/usr/bin/env python
"""Minimal gdTAI v2.0 example.

This delegates to the repository inference script while making the operating
mode explicit. Use `high_f1` for best validation F1 or `high_purity` for lower FP.
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[5]
PYTHON = "/home/tanlikai/miniconda3/envs/rapids_sc_py310/bin/python"

parser = argparse.ArgumentParser(description="Run gdTAI v2.0 minimal prediction.")
parser.add_argument("--input-h5ad", required=True)
parser.add_argument("--dataset-id", required=True)
parser.add_argument("--mode", choices=["high_f1", "high_purity"], default="high_f1")
parser.add_argument("--annotation-column", default="simple_annotation_plus6")
args = parser.parse_args()

cmd = [
    PYTHON,
    str(REPO_ROOT / "predict_with_gdtai_v2.py"),
    "--input-h5ad", args.input_h5ad,
    "--dataset-id", args.dataset_id,
    "--mode", args.mode,
]
if args.mode == "high_purity":
    cmd.extend(["--annotation-column", args.annotation_column])
subprocess.run(cmd, check=True, cwd=REPO_ROOT)
