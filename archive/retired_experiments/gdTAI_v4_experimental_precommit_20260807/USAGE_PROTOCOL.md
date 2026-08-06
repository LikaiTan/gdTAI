# gdTAI v4.0 Usage Protocol

Status: experimental candidate, not the promoted default model.

Input expression must be log1p(counts per 10,000) for the genes listed in `feature_genes.csv`, followed by the engineered features defined in `run_gdtai_v4_tcell_gate_classifier.py`.

Inference steps:

1. Build the stage-1 feature matrix using the listed gene features plus engineered features.
2. Load `gdTAI_v4_model.pkl`.
3. Compute `tcell_score = payload["tcell_model"].predict_proba(X)[:, 1]`.
4. Build `X_stage2 = column_stack([X, tcell_score])`.
5. Compute `gdt_score = payload["gdt_model"].predict_proba(X_stage2)[:, 1]`.
6. For an operating mode, call gdT if `tcell_score >= payload["tcell_threshold"]` and `gdt_score >= payload["operating_modes"][mode]`.

Available modes: `high_f1`, `high_purity`.
