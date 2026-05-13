#!/usr/bin/env python
"""Audit likely false positives from a selected gdT model prediction table.

This script expects the CSV written by predict_with_selected_gdt_model.py or by
minimal_python_prediction.py. It treats predicted gdT cells with alpha-beta TCR
evidence or NK-like annotation as false-positive compartments for review.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


TRUE_VALUES = {"1", "true", "t", "yes", "y", "present", "positive"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit likely gdT classifier false positives.")
    parser.add_argument("--prediction-csv", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument(
        "--ab-column",
        action="append",
        default=["has_TRA_TRB_paired", "has_any_ab_tcr"],
        help="Boolean column indicating alpha-beta TCR evidence. Can be repeated.",
    )
    parser.add_argument(
        "--annotation-column",
        action="append",
        default=["simple_annotation_plus6", "annotation", "cell_type"],
        help="Annotation column to scan for NK-like labels. Can be repeated.",
    )
    parser.add_argument("--nk-regex", default=r"\bNK\b|natural killer", help="Case-insensitive NK label regex.")
    parser.add_argument("--group-column", action="append", default=[], help="Optional grouping column.")
    return parser.parse_args()


def as_bool(series: pd.Series) -> pd.Series:
    text = series.fillna("").astype(str).str.strip().str.lower()
    return text.isin(TRUE_VALUES)


def summarize_mask(df: pd.DataFrame, mask: pd.Series, label: str, group_col: str | None = None) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if group_col is None:
        groups = [("overall", df)]
    else:
        groups = [(str(name), group) for name, group in df.groupby(group_col, dropna=False, sort=True)]
    for group_name, group in groups:
        group_mask = mask.loc[group.index]
        denominator = int(group_mask.sum())
        predicted_in_mask = int((group_mask & group["predicted_gdT"]).sum())
        rows.append(
            {
                "scope": label,
                "group_column": "" if group_col is None else group_col,
                "group_value": group_name,
                "total_cells": int(group.shape[0]),
                "compartment_cells": denominator,
                "predicted_gdT_in_compartment": predicted_in_mask,
                "compartment_fp_rate": predicted_in_mask / denominator if denominator else 0.0,
                "whole_group_fp_fraction": predicted_in_mask / int(group.shape[0]) if group.shape[0] else 0.0,
            }
        )
    return rows


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.prediction_csv)
    if "predicted_gdT" not in df.columns:
        raise KeyError("Prediction table must contain a predicted_gdT column.")
    df["predicted_gdT"] = as_bool(df["predicted_gdT"])

    rows: list[dict[str, object]] = []

    ab_columns = [col for col in args.ab_column if col in df.columns]
    if ab_columns:
        ab_mask = pd.Series(False, index=df.index)
        for col in ab_columns:
            ab_mask |= as_bool(df[col])
        rows.extend(summarize_mask(df, ab_mask, "alpha_beta_tcr_evidence", None))
        for group_col in args.group_column:
            if group_col in df.columns:
                rows.extend(summarize_mask(df, ab_mask, "alpha_beta_tcr_evidence", group_col))

    annotation_columns = [col for col in args.annotation_column if col in df.columns]
    if annotation_columns:
        pattern = re.compile(args.nk_regex, flags=re.IGNORECASE)
        nk_mask = pd.Series(False, index=df.index)
        for col in annotation_columns:
            nk_mask |= df[col].fillna("").astype(str).map(lambda value: bool(pattern.search(value)))
        rows.extend(summarize_mask(df, nk_mask, "nk_like_annotation", None))
        for group_col in args.group_column:
            if group_col in df.columns:
                rows.extend(summarize_mask(df, nk_mask, "nk_like_annotation", group_col))

    if not rows:
        available = ", ".join(df.columns)
        raise KeyError(
            "No requested AB-TCR or annotation columns were present in the prediction table. "
            f"Available columns: {available}"
        )

    out = pd.DataFrame(rows)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output_csv, index=False)
    print(f"wrote {args.output_csv}")


if __name__ == "__main__":
    main()

