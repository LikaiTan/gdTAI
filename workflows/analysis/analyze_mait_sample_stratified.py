#!/usr/bin/env python3
"""Sample-aware MAIT-related TCR comparison for TRD-positive vs TRD-negative T cells.

This helper avoids relying only on pooled counts. It computes the same
MAIT-related chain features within each sample and summarizes whether the
feature is consistently enriched in `TRD < 0` or `TRD > 0` cells after
conditioning on sample identity.
"""

from __future__ import annotations

# TNK_WORKFLOW_BOOTSTRAP
import sys as _tnk_sys
from pathlib import Path as _TnkPath

_TNK_PROJECT_ROOT = _TnkPath(__file__).resolve().parents[2]
for _tnk_path in (
    _TNK_PROJECT_ROOT,
    _TNK_PROJECT_ROOT / "src",
    _TNK_PROJECT_ROOT / "workflows" / "integration",
    _TNK_PROJECT_ROOT / "workflows" / "intake",
    _TNK_PROJECT_ROOT / "workflows" / "metadata",
    _TNK_PROJECT_ROOT / "workflows" / "analysis",
    _TNK_PROJECT_ROOT / "workflows" / "gdtai",
    _TNK_PROJECT_ROOT / "workflows" / "gdt_atlas",
    _TNK_PROJECT_ROOT / "workflows" / "reporting",
    _TNK_PROJECT_ROOT / "workflows" / "maintenance",
):
    _tnk_value = str(_tnk_path)
    if _tnk_value not in _tnk_sys.path:
        _tnk_sys.path.insert(0, _tnk_value)


from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


PROJECT_ROOT = Path(__file__).resolve().parents[2]
INTEGRATED_H5AD = PROJECT_ROOT / "high_speed_temp" / "Integrated_dataset" / "integrated.h5ad"
TABLE_DIR = PROJECT_ROOT / "Integrated_dataset" / "tables"
LOG_DIR = PROJECT_ROOT / "Integrated_dataset" / "logs"

OUT_SAMPLE_CSV = TABLE_DIR / "mait_tcr_chain_enrichment_by_sampleid.csv"
OUT_GSE_CSV = TABLE_DIR / "mait_tcr_chain_enrichment_by_gse.csv"
OUT_SUMMARY_CSV = TABLE_DIR / "mait_tcr_chain_enrichment_sample_stratified_summary.csv"
OUT_MD = LOG_DIR / "mait_tcr_chain_enrichment_sample_stratified.md"


def read_string_array(obs: h5py.Group, key: str) -> np.ndarray:
    """Read one string-like obs column from H5AD."""
    obj = obs[key]
    if isinstance(obj, h5py.Group) and "codes" in obj and "categories" in obj:
        codes = obj["codes"][:]
        categories = obj["categories"].asstr()[:]
        out = np.empty(len(codes), dtype=object)
        neg_mask = codes < 0
        out[neg_mask] = ""
        out[~neg_mask] = categories[codes[~neg_mask]]
        return out.astype(str)
    return obj.asstr()[:].astype(str)


def contains_any(values: np.ndarray, needles: list[str]) -> np.ndarray:
    """Return True where any needle appears in the uppercased value."""
    normalized = np.char.upper(np.char.strip(values.astype(str)))
    out = np.zeros(normalized.shape[0], dtype=bool)
    for needle in needles:
        out |= np.char.find(normalized, needle.upper()) >= 0
    return out


def build_feature_frame() -> pd.DataFrame:
    """Build the filtered cell-level frame with MAIT-related flags."""
    with h5py.File(INTEGRATED_H5AD, "r") as handle:
        obs = handle["obs"]
        frame = pd.DataFrame(
            {
                "sampleid": read_string_array(obs, "sampleid"),
                "source_gse_id": read_string_array(obs, "source_gse_id"),
                "scanvi_tnk_superclass": read_string_array(obs, "scanvi_tnk_superclass"),
                "TRA_cdr3": read_string_array(obs, "TRA_cdr3"),
                "TRB_cdr3": read_string_array(obs, "TRB_cdr3"),
                "TRA_v": read_string_array(obs, "TRA_v"),
                "TRA_j": read_string_array(obs, "TRA_j"),
                "TRB_v": read_string_array(obs, "TRB_v"),
                "phase4_trab_score": np.asarray(obs["phase4_trab_score"][:], dtype=np.float32),
                "phase4_trd_score": np.asarray(obs["phase4_trd_score"][:], dtype=np.float32),
            }
        )

    paired = frame["TRA_cdr3"].str.strip().ne("") & frame["TRB_cdr3"].str.strip().ne("")
    base = (
        frame["scanvi_tnk_superclass"].eq("T_cell")
        & paired
        & (frame["phase4_trab_score"] > -0.05)
        & frame["sampleid"].astype(str).str.strip().ne("")
    )
    frame = frame.loc[base].copy()
    frame["trd_group"] = np.where(
        frame["phase4_trd_score"] > 0,
        "TRD>0",
        np.where(frame["phase4_trd_score"] < 0, "TRD<0", "TRD=0"),
    )
    frame = frame.loc[frame["trd_group"].isin(["TRD>0", "TRD<0"])].copy()

    trav1_2 = contains_any(frame["TRA_v"].to_numpy(), ["TRAV1-2", "TRAV1_2", "V7.2"])
    traj33 = contains_any(frame["TRA_j"].to_numpy(), ["TRAJ33"])
    traj12 = contains_any(frame["TRA_j"].to_numpy(), ["TRAJ12"])
    traj20 = contains_any(frame["TRA_j"].to_numpy(), ["TRAJ20"])
    trbv20_or_6 = contains_any(frame["TRB_v"].to_numpy(), ["TRBV20", "TRBV6"])

    frame["TRAV1-2"] = trav1_2
    frame["TRAJ33"] = traj33
    frame["TRAJ12"] = traj12
    frame["TRAJ20"] = traj20
    frame["MAIT_alpha"] = trav1_2 & (traj33 | traj12 | traj20)
    frame["TRBV20_or_TRBV6_family"] = trbv20_or_6
    frame["MAIT_like_full_chain"] = frame["MAIT_alpha"] & frame["TRBV20_or_TRBV6_family"]
    return frame


def summarize_within_group(
    frame: pd.DataFrame,
    group_column: str,
    feature_columns: list[str],
) -> pd.DataFrame:
    """Summarize within-sample or within-GSE feature enrichment."""
    rows: list[dict[str, object]] = []
    for group_value, group_df in frame.groupby(group_column, sort=False):
        n_pos = int((group_df["trd_group"] == "TRD>0").sum())
        n_neg = int((group_df["trd_group"] == "TRD<0").sum())
        if n_pos == 0 or n_neg == 0:
            continue
        gse_id = group_df["source_gse_id"].iloc[0] if "source_gse_id" in group_df.columns else ""
        for feature in feature_columns:
            yes_pos = int(((group_df["trd_group"] == "TRD>0") & group_df[feature]).sum())
            yes_neg = int(((group_df["trd_group"] == "TRD<0") & group_df[feature]).sum())
            no_pos = n_pos - yes_pos
            no_neg = n_neg - yes_neg
            odds_ratio, p_value = fisher_exact([[yes_pos, no_pos], [yes_neg, no_neg]], alternative="two-sided")
            rows.append(
                {
                    group_column: group_value,
                    "source_gse_id": gse_id,
                    "feature": feature,
                    "trd_gt_0_n": n_pos,
                    "trd_lt_0_n": n_neg,
                    "trd_gt_0_yes": yes_pos,
                    "trd_lt_0_yes": yes_neg,
                    "trd_gt_0_fraction": yes_pos / n_pos,
                    "trd_lt_0_fraction": yes_neg / n_neg,
                    "odds_ratio_trd_gt0_vs_lt0": float(odds_ratio),
                    "fisher_pvalue": float(p_value),
                    "direction": (
                        "higher_in_TRD>0"
                        if yes_pos / n_pos > yes_neg / n_neg
                        else "higher_in_TRD<0"
                        if yes_pos / n_pos < yes_neg / n_neg
                        else "equal"
                    ),
                }
            )
    return pd.DataFrame(rows)


def write_summary_markdown(summary_df: pd.DataFrame) -> None:
    """Write the high-level sample-aware summary."""
    lines = [
        "# Sample-Aware MAIT-Related TCR Comparison",
        "",
        "- Base filter: `scanvi_tnk_superclass == T_cell`, paired `TRA/TRB`, `phase4_trab_score > -0.05`",
        "- Comparison: `TRD > 0` versus `TRD < 0`",
        "- This summary is sample-aware and asks whether the direction is consistent across `sampleid`, not only in pooled counts.",
        "",
    ]
    for row in summary_df.itertuples(index=False):
        lines.append(
            f"- {row.feature}: samples with both groups = `{row.n_samples_with_both_groups}`, "
            f"`TRD<0` higher in `{row.samples_higher_in_trd_lt_0}` samples, "
            f"`TRD>0` higher in `{row.samples_higher_in_trd_gt_0}` samples, "
            f"median sample odds ratio = `{row.median_sample_odds_ratio:.4f}`"
        )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Run the sample-aware MAIT-related chain comparison."""
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    frame = build_feature_frame()
    feature_columns = [
        "TRAV1-2",
        "TRAJ33",
        "TRAJ12",
        "TRAJ20",
        "MAIT_alpha",
        "TRBV20_or_TRBV6_family",
        "MAIT_like_full_chain",
    ]
    sample_df = summarize_within_group(frame, "sampleid", feature_columns)
    gse_df = summarize_within_group(frame, "source_gse_id", feature_columns)

    sample_df.to_csv(OUT_SAMPLE_CSV, index=False)
    gse_df.to_csv(OUT_GSE_CSV, index=False)

    summary_rows: list[dict[str, object]] = []
    for feature in feature_columns:
        feature_df = sample_df.loc[sample_df["feature"] == feature].copy()
        finite_or = feature_df["odds_ratio_trd_gt0_vs_lt0"].replace([np.inf, -np.inf], np.nan).dropna()
        summary_rows.append(
            {
                "feature": feature,
                "n_samples_with_both_groups": int(len(feature_df)),
                "samples_higher_in_trd_lt_0": int((feature_df["direction"] == "higher_in_TRD<0").sum()),
                "samples_higher_in_trd_gt_0": int((feature_df["direction"] == "higher_in_TRD>0").sum()),
                "samples_equal": int((feature_df["direction"] == "equal").sum()),
                "median_sample_odds_ratio": float(np.median(finite_or)) if len(finite_or) else np.nan,
            }
        )
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUT_SUMMARY_CSV, index=False)
    write_summary_markdown(summary_df)
    print(summary_df.to_csv(index=False))


if __name__ == "__main__":
    main()
