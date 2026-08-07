#!/usr/bin/env python3
"""Render the gdTAI-kimi experimental comparison report (MD + figures).

Reads Integrated_dataset/tables/gdT_prediction/gdtai_kimi/ and writes
- Integrated_dataset/logs/gdT_prediction/gdtai_kimi/gdtai_kimi_report.md
- Integrated_dataset/figures/gdT_prediction/gdtai_kimi/*.png
No model artifacts, H5ADs, or registry entries are touched.
"""

from __future__ import annotations

import sys as _tnk_sys
from pathlib import Path as _TnkPath

_TNK_PROJECT_ROOT = _TnkPath(__file__).resolve().parents[2]
for _p in (_TNK_PROJECT_ROOT, _TNK_PROJECT_ROOT / "src", _TNK_PROJECT_ROOT / "workflows" / "gdtai"):
    if str(_p) not in _tnk_sys.path:
        _tnk_sys.path.insert(0, str(_p))

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
TABLE_DIR = PROJECT_ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_kimi"
FIGURE_DIR = PROJECT_ROOT / "Integrated_dataset/figures/gdT_prediction/gdtai_kimi"
LOG_DIR = PROJECT_ROOT / "Integrated_dataset/logs/gdT_prediction/gdtai_kimi"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

CANDIDATE_LABELS = {
    "kimi_en_balanced": "kimi elastic-net (balanced)",
    "kimi_hgb_balanced": "kimi HistGBM (balanced)",
    "kimi_en_highpurity": "kimi elastic-net (high-purity)",
    "kimi_hgb_highpurity": "kimi HistGBM (high-purity)",
    "kimi_en_balanced_cd4treg_excl": "kimi EN balanced + CD4/Treg exclusion",
    "kimi_hgb_balanced_cd4treg_excl": "kimi HGB balanced + CD4/Treg exclusion",
    "C1_phase4_delta": "C1 phase4 TRD-TRAB score",
    "C2_compact7": "C2 compact 7-gene logistic",
    "C3_tcr_en": "C3 TCR-gene elastic-net",
    "frozen_v2_high_f1": "frozen V2 high-F1 (descriptive)",
    "frozen_v2_high_purity": "frozen V2 high-purity (descriptive)",
    "frozen_v3_round14_balanced": "frozen V3 Round 14 (descriptive)",
    "frozen_v3_round12_high_purity": "frozen V3 Round 12 (descriptive)",
}


def main() -> None:
    res = pd.read_csv(TABLE_DIR / "nested_metrics_by_heldout.csv")
    macro = pd.read_csv(TABLE_DIR / "nested_metrics_macro.csv")
    cis = pd.read_csv(TABLE_DIR / "bootstrap_macro_f1_cis.csv")
    labels = pd.read_csv(TABLE_DIR / "label_inventory.csv")
    choices = json.loads((LOG_DIR / "fold_model_choices.json").read_text())

    primary = res[~res.held_out.astype(str).str.startswith("stress")].copy()
    stress = res[res.held_out.astype(str).str.startswith("stress")].copy()
    primary["label"] = primary.candidate.map(CANDIDATE_LABELS).fillna(primary.candidate)
    macro["label"] = macro.candidate.map(CANDIDATE_LABELS).fillna(macro.candidate)

    # ---- figure 1: per-dataset recall and abT FPR for key candidates
    key = ["kimi_en_balanced", "kimi_hgb_balanced", "C1_phase4_delta", "C2_compact7", "C3_tcr_en",
           "frozen_v2_high_purity", "frozen_v3_round12_high_purity", "frozen_v3_round14_balanced"]
    key = [c for c in key if c in set(primary.candidate)]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    for ax, metric, title in zip(axes, ["recall", "abt_fpr", "nk_fpr"],
                                 ["Recall (gdT_primary)", "Paired-abT FPR", "Strict-NK FPR"]):
        data = primary[primary.candidate.isin(key)]
        xpos = np.arange(len(key))
        width = 0.25
        for i, ho in enumerate(sorted(data.held_out.unique())):
            vals = [data[(data.candidate == c) & (data.held_out == ho)][metric].iloc[0]
                    if len(data[(data.candidate == c) & (data.held_out == ho)]) else np.nan for c in key]
            ax.bar(xpos + i * width, vals, width, label=ho)
        ax.set_xticks(xpos + width, [CANDIDATE_LABELS.get(c, c).replace(" (", "\n(") for c in key],
                      rotation=35, ha="right", fontsize=7)
        ax.set_title(title)
        if metric != "recall":
            ax.set_yscale("log")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "nested_metrics_by_dataset.png", dpi=200)
    plt.close(fig)

    # ---- figure 2: macro F1 with bootstrap CIs
    m = macro[macro.candidate.isin(key)].merge(cis, on="candidate", how="left")
    m = m.sort_values("macro_f1")
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ypos = np.arange(len(m))
    ax.barh(ypos, m.macro_f1, color=["#0b3d66" if str(c).startswith("kimi") else "#9db8cc" for c in m.candidate])
    if {"macro_f1_ci_lo", "macro_f1_ci_hi"}.issubset(m.columns):
        ax.errorbar(m.macro_f1, ypos,
                    xerr=[m.macro_f1 - m.macro_f1_ci_lo, m.macro_f1_ci_hi - m.macro_f1],
                    fmt="none", ecolor="black", capsize=3)
    ax.set_yticks(ypos, [CANDIDATE_LABELS.get(c, c) for c in m.candidate], fontsize=8)
    ax.set_xlabel("Dataset-macro F1 (LODO held-out mean; error bars: donor-cluster bootstrap 95% CI)")
    ax.set_title("gdTAI-kimi vs comparators and frozen profiles (descriptive)")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "macro_f1_bootstrap.png", dpi=200)
    plt.close(fig)

    # ---- report
    def fmt(v, nd=4):
        return "" if pd.isna(v) else f"{v:.{nd}f}"

    lines = []
    lines.append("# gdTAI-kimi: experimental nested cross-dataset evaluation\n")
    lines.append("**Status:** experimental candidate, separately named from the frozen GDTAI V4 "
                 "precommitted protocol. Not registered, not promotable. Frozen V2/V3 profiles are "
                 "descriptive references developed with access to some held-out cells (biased in "
                 "their favor).\n")
    lines.append("## Label inventory (expression-independent TCR labels)\n")
    lines.append(labels.to_markdown(index=False))
    lines.append("\n\n## Fold model choices\n")
    lines.append("```json")
    lines.append(json.dumps(choices, indent=2, default=str)[:4000])
    lines.append("```\n")
    lines.append("## Per-held-out-dataset metrics\n")
    show = primary[["held_out", "label", "threshold", "precision", "recall", "specificity",
                    "f1", "roc_auc", "pr_auc", "abt_fpr", "nk_fpr", "ece",
                    "inner_guardrails_met"]].copy()
    for c in ("precision", "recall", "specificity", "f1", "roc_auc", "pr_auc", "ece"):
        show[c] = show[c].map(lambda v: fmt(v, 4))
    for c in ("abt_fpr", "nk_fpr"):
        show[c] = show[c].map(lambda v: fmt(v, 5))
    lines.append(show.to_markdown(index=False))
    lines.append("\n\n## Dataset-macro metrics with donor-cluster bootstrap CIs\n")
    mm = macro.merge(cis, on="candidate", how="left")
    mm["label"] = mm.candidate.map(CANDIDATE_LABELS).fillna(mm.candidate)
    mm["macro_f1_ci"] = mm.apply(lambda r: f"[{r.macro_f1_ci_lo:.4f}, {r.macro_f1_ci_hi:.4f}]"
                                 if pd.notna(r.get("macro_f1_ci_lo")) else "", axis=1)
    showm = mm[["label", "macro_f1", "macro_f1_ci", "macro_recall", "macro_precision",
                "macro_abt_fpr", "macro_nk_fpr", "min_recall", "max_abt_fpr"]].copy()
    for c in showm.columns:
        if c.startswith("macro_") and c != "macro_f1_ci" or c in ("min_recall", "max_abt_fpr"):
            showm[c] = showm[c].map(lambda v: fmt(v, 4))
    lines.append(showm.to_markdown(index=False))
    if len(stress):
        lines.append("\n\n## GSE254249 stress (censored paired-abT FPR = upper bound)\n")
        lines.append(stress[["held_out", "candidate", "abt_fpr", "nk_fpr"]].to_markdown(index=False))
    lines.append("\n\n## CD4/Treg exclusion cost\n")
    excl = primary[primary.candidate.str.endswith("cd4treg_excl")]
    if len(excl):
        lines.append(excl[["held_out", "candidate", "recall", "f1", "fn_cost_vs_balanced"]].to_markdown(index=False))
    lines.append("\n\n![per-dataset metrics](../../../figures/gdT_prediction/gdtai_kimi/nested_metrics_by_dataset.png)")
    lines.append("\n![macro F1](../../../figures/gdT_prediction/gdtai_kimi/macro_f1_bootstrap.png)")
    (LOG_DIR / "gdtai_kimi_report.md").write_text("\n".join(lines))
    print("wrote", LOG_DIR / "gdtai_kimi_report.md")


if __name__ == "__main__":
    main()
