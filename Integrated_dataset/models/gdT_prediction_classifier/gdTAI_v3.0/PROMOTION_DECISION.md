# gdTAI v3.0 Promotion Decision

Round 14 remains the canonical balanced gdTAI v3.0 model after a new
checksum-pinned Round 12 versus Round 14 comparison completed on 2026-07-15.

- Promoted model: `v3_round14_v2_score_trdc_gate_fixed_0p936`
- Threshold: `0.936`
- Model SHA256: `16dedc0081da9b8487887341232bcf8c9c9403dd3bbd72e04cab43d4cd7b2e09`
- Full-atlas predicted cells: `251,356`
- Full-atlas primary-gold recall / F1: `0.8056` / `0.8823`
- Atlas-held-out recall / F1: `0.8530` / `0.9144`
- Independent external precision / recall / F1: `0.9509` / `0.8586` / `0.9024`

Selection rule: Pass all four safety/recall guardrails, then maximize mean F1 across full-atlas gold, atlas-held-out, and independent-external evaluations.

Round 12 is retained as the high-purity fallback. The detailed comparison is
served at `gdT_prediction/gdtai_v3_round12_vs_round14/index.html`.
