#!/usr/bin/env python3
"""Export conservative NK-like negative-feature candidates from boundary DEGs."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
INPUT = (
    ROOT
    / "Integrated_dataset/tables/gdT_prediction/"
    "gdtai_v4_2_boundary_pseudobulk_deg/deg_results.csv.gz"
)
OUTPUT_DIR = (
    ROOT
    / "Integrated_dataset/tables/gdT_prediction/"
    "gdtai_v4_2_boundary_pseudobulk_deg"
)
LOG_DIR = (
    ROOT
    / "Integrated_dataset/logs/gdT_prediction/"
    "gdtai_v4_2_boundary_pseudobulk_deg"
)


# These candidates were selected for innate/NK receptor or signaling biology
# that is not part of the canonical T-cell receptor/cytotoxicity program.
CANDIDATES = {
    "NCR1": ("A", "NK natural-cytotoxicity receptor"),
    "SIGLEC7": ("A", "innate inhibitory lectin receptor enriched in NK cells"),
    "SH2D1B": ("A", "innate SLAM-family adaptor, unlike the T-cell SAP adaptor"),
    "LAT2": ("B", "non-T immune-receptor adaptor; distinct from T-cell LAT"),
    "SYK": ("B", "innate immune-receptor kinase; distinct from T-cell ZAP70"),
    "PLCG2": ("B", "innate/B-cell phospholipase; distinct from T-cell PLCG1"),
    "PILRB": ("B", "paired innate leukocyte receptor"),
    "CD300C": ("B", "activating innate myeloid/NK-associated receptor"),
}

# Explicitly rejected despite some NK association. These genes are shared with
# cytotoxic, activated, memory, innate-like, or gamma-delta T-cell programs.
SHARED_T_CELL_EXCLUSIONS = {
    "FCER1G",
    "TYROBP",
    "GNLY",
    "NKG7",
    "KLRD1",
    "KLRC1",
    "GZMK",
    "XCL1",
    "XCL2",
    "CD160",
    "NCAM1",
    "CTSW",
    "PRF1",
    "CCL5",
    "ZBTB16",
    "IL12RB2",
    "IL18RAP",
    "TMIGD2",
    "CD7",
    "CD27",
    "CCR7",
    "BACH2",
    "SATB1",
    "SELL",
}


def main() -> None:
    deg = pd.read_csv(INPUT)
    deg = deg.set_index("gene", drop=False)

    missing = sorted(set(CANDIDATES) - set(deg.index))
    if missing:
        raise RuntimeError(f"Candidate genes absent from DEG table: {missing}")

    selected = deg.loc[list(CANDIDATES)].copy()
    selected.insert(1, "candidate_tier", [CANDIDATES[g][0] for g in selected.index])
    selected.insert(2, "biological_rationale", [CANDIDATES[g][1] for g in selected.index])
    selected.insert(3, "intended_model_role", "soft NK/non-T negative feature")

    checks = (
        selected["robust_cross_dataset"].astype(bool)
        & (selected["mean_log2FC_enriched_vs_depleted"] <= -0.5)
        & (selected["dataset_macro_mean_log2FC"] < 0)
        & (selected["source_sign_consistency"] >= 0.88)
    )
    if not checks.all():
        failed = selected.loc[~checks, "gene"].tolist()
        raise RuntimeError(f"Candidates failed the frozen evidence gates: {failed}")

    if any(g.startswith("KIR") for g in selected.index):
        raise RuntimeError("KIR genes must not enter the candidate list")
    forbidden = SHARED_T_CELL_EXCLUSIONS.intersection(selected.index)
    if forbidden:
        raise RuntimeError(f"Shared T-cell genes entered candidates: {sorted(forbidden)}")

    output_cols = [
        "gene",
        "candidate_tier",
        "biological_rationale",
        "intended_model_role",
        "mean_log2FC_enriched_vs_depleted",
        "dataset_macro_mean_log2FC",
        "source_sign_consistency",
        "paired_fdr",
        "dataset_macro_fdr",
        "robust_cross_dataset",
    ]
    selected[output_cols].to_csv(
        OUTPUT_DIR / "candidate_nk_negative_features.csv", index=False
    )

    rejected = []
    for gene in sorted(SHARED_T_CELL_EXCLUSIONS):
        if gene in deg.index:
            rejected.append(
                {
                    "gene": gene,
                    "reason": "shared with cytotoxic/activated/innate-like T cells",
                    "mean_log2FC_enriched_vs_depleted": deg.at[
                        gene, "mean_log2FC_enriched_vs_depleted"
                    ],
                    "robust_cross_dataset": deg.at[gene, "robust_cross_dataset"],
                }
            )
    rejected.append(
        {
            "gene": "KIR*",
            "reason": "all KIR genes excluded by user-specified policy",
            "mean_log2FC_enriched_vs_depleted": pd.NA,
            "robust_cross_dataset": pd.NA,
        }
    )
    pd.DataFrame(rejected).to_csv(
        OUTPUT_DIR / "excluded_shared_nk_t_features.csv", index=False
    )

    lines = [
        "# Conservative NK-like negative-feature candidates",
        "",
        "- status: `PASS_FEATURE_TRIAGE_ONLY`",
        f"- candidates: `{', '.join(selected.index)}`",
        "- Tier A: direct innate/NK receptor-adaptor evidence",
        "- Tier B: non-T immune-receptor signaling or innate-receptor evidence",
        "- every candidate is a robust cross-dataset depleted-cluster DEG with paired log2FC <= -0.5 and source sign consistency >= 0.88",
        "- excluded: `GZMK`, `KLRC1`, all `KIR*`, cytotoxicity genes, and genes commonly shared with activated or gamma-delta T cells",
        "- intended use: soft negative features evaluated inside grouped cross-dataset folds",
        "- prohibited use: no single candidate is a hard NK veto or a ground-truth label",
        "",
        "The shortlist is intentionally small. Its predictive value and gamma-delta T-cell recall cost must be estimated fold-locally before model inclusion.",
    ]
    (LOG_DIR / "nk_negative_feature_selection.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    print("PASS_FEATURE_TRIAGE_ONLY", ",".join(selected.index))


if __name__ == "__main__":
    main()
