import importlib.util
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflows/gdtai"))

from gdtai_v4_nested_core import (  # noqa: E402
    apply_two_stage_call,
    balanced_training_weights,
    bootstrap_interval,
    derive_features,
    exclusion_flags,
    fit_binary_estimator,
    fit_platt_calibrator,
    fold_feature_mask,
    grouped_confusion_counts,
    inference_coverage,
    paired_hierarchical_bootstrap_f1_difference,
    parameter_grid,
    select_stage1_threshold,
    select_stage2_threshold,
    vdj_rescue_mask,
)


SCRIPT = ROOT / "workflows/gdtai/run_gdtai_v4_nested_evaluation.py"
SPEC = importlib.util.spec_from_file_location("run_gdtai_v4_nested_evaluation", SCRIPT)
RUNNER = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = RUNNER
SPEC.loader.exec_module(RUNNER)


class GdtaiV4NestedEvaluationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = json.loads((ROOT / "configs/models/gdtai/v4_nested_evaluation.json").read_text())
        cls.features = pd.read_csv(
            ROOT / cls.config["preflight"]["feature_manifest"]
        ).sort_values("feature_index")["gene"].astype(str).tolist()

    def test_frozen_config_matches_protocol_grid_and_feature_cap(self) -> None:
        models = self.config["models"]
        self.assertEqual(len(parameter_grid({"C": models["stage1_elastic_net"]["C"], "l1_ratio": models["stage1_elastic_net"]["l1_ratio"]})), 16)
        self.assertEqual(len(parameter_grid({"C": models["stage2_elastic_net"]["C"], "l1_ratio": models["stage2_elastic_net"]["l1_ratio"]})), 25)
        hgb = models["stage2_hist_gradient_boosting"]
        self.assertEqual(
            len(parameter_grid({key: hgb[key] for key in ("learning_rate", "max_leaf_nodes", "min_samples_leaf", "l2_regularization")})),
            16,
        )
        self.assertEqual(len(self.features), 197)
        self.assertLessEqual(len(self.features), self.config["feature_policy"]["maximum_retained_genes"])
        self.assertEqual(len(self.config["feature_policy"]["stage1_genes"]), 50)

    def test_derived_features_are_exact_unweighted_summaries(self) -> None:
        x = np.zeros((2, len(self.features)), dtype=np.float32)
        lookup = {gene: index for index, gene in enumerate(self.features)}
        x[0, lookup["TRAV1-1"]] = 1.0
        x[0, lookup["TRAC"]] = 2.0
        x[0, lookup["TRDV2"]] = 3.0
        x[0, lookup["TRDC"]] = 4.0
        x[0, lookup["CD3D"]] = 2.0
        derived, names = derive_features(
            x,
            self.features,
            self.config["feature_policy"]["family_prefixes"],
            self.config["feature_policy"]["derived_panels"],
        )
        values = dict(zip(names, derived[0], strict=True))
        self.assertEqual(values["TRA_detected"], 2.0)
        self.assertEqual(values["TRD_detected"], 2.0)
        self.assertEqual(values["TRB_detected"], 0.0)
        self.assertAlmostEqual(values["T_lineage_mean"], 2.0 / 14.0, places=6)

    def test_cd4_and_treg_exclusions_use_frozen_composite_rules(self) -> None:
        x = np.zeros((3, len(self.features)), dtype=np.float32)
        lookup = {gene: index for index, gene in enumerate(self.features)}
        for gene in ["CD4", "IL7R", "CCR7"]:
            x[0, lookup[gene]] = 1.2
        x[1, lookup["CD4"]] = 3.0
        for gene in ["FOXP3", "IL2RA"]:
            x[2, lookup[gene]] = 1.2
        cd4, treg, union = exclusion_flags(
            x,
            self.features,
            self.config["cd4_helper_rule"],
            self.config["treg_rule"],
        )
        self.assertEqual(cd4.tolist(), [True, False, False])
        self.assertEqual(treg.tolist(), [False, False, True])
        self.assertEqual(union.tolist(), [True, False, True])

    def test_preflight_exclusion_scope_is_reproduced_before_all_cell_application(self) -> None:
        x = np.zeros((2, len(self.features)), dtype=np.float32)
        lookup = {gene: index for index, gene in enumerate(self.features)}
        for row in (0, 1):
            for gene in ["CD4", "IL7R", "CCR7"]:
                x[row, lookup[gene]] = 1.2
        cells = pd.DataFrame(
            {
                "stage1_role": ["t_positive", "t_positive"],
                "truth_class": ["gdT_primary", "abT_primary"],
                "cd4_helper_exclusion": [True, False],
                "treg_exclusion": [False, False],
                "exclusion_union": [True, False],
            }
        )
        audit = RUNNER.apply_frozen_exclusions(cells, x, self.features, self.config)
        self.assertTrue(audit["preflight_scope_reproduction_pass"])
        self.assertEqual(audit["preflight_scope_cells"], 1)
        self.assertEqual(cells["cd4_helper_exclusion"].tolist(), [True, True])
        self.assertEqual(cells["exclusion_union"].tolist(), [True, True])

        cells.loc[0, "cd4_helper_exclusion"] = False
        cells.loc[0, "exclusion_union"] = False
        with self.assertRaises(RuntimeError):
            RUNNER.apply_frozen_exclusions(cells, x, self.features, self.config)

    def test_balanced_weights_include_reliability_and_normalize(self) -> None:
        y = np.array([1, 1, 0, 0, 0, 1, 0], dtype=np.int8)
        source = np.array(["A", "A", "A", "A", "A", "B", "B"], dtype=object)
        reliability = np.array([1, 1, 1, 1, 0.5, 1, 1], dtype=float)
        weights = balanced_training_weights(y, source, reliability)
        self.assertAlmostEqual(float(weights.mean()), 1.0, places=6)
        self.assertLess(weights[4], weights[2])
        self.assertTrue(np.isfinite(weights).all())

    def test_fold_filter_removes_constant_and_rare_genes(self) -> None:
        x = np.zeros((1000, 4), dtype=np.float32)
        x[:10, 0] = 1.0
        x[:, 1] = 1.0
        x[0, 2] = 1.0
        x[:, 3] = np.arange(1000)
        keep = fold_feature_mask(x, gene_feature_count=3, minimum_detection_fraction=0.005, maximum_retained_genes=220)
        self.assertEqual(keep.tolist(), [True, False, False, True])

    def test_platt_calibration_and_model_are_seed_deterministic(self) -> None:
        rng = np.random.default_rng(22)
        x = rng.normal(size=(300, 5)).astype(np.float32)
        y = (x[:, 0] - 0.3 * x[:, 1] > 0).astype(np.int8)
        w = np.ones(y.size, dtype=np.float32)
        params = {"C": 0.1, "l1_ratio": 0.25, "max_iter": 5000, "tolerance": 1e-4}
        first = fit_binary_estimator(x, y, w, [f"g{i}" for i in range(5)], 5, "elastic_net", params, 0.0005, 220, 17)
        second = fit_binary_estimator(x, y, w, [f"g{i}" for i in range(5)], 5, "elastic_net", params, 0.0005, 220, 17)
        raw_first = first.predict_probability(x)
        raw_second = second.predict_probability(x)
        np.testing.assert_allclose(raw_first, raw_second, rtol=0, atol=1e-12)
        cal_first = fit_platt_calibrator(raw_first, y, w, 17)
        cal_second = fit_platt_calibrator(raw_second, y, w, 17)
        np.testing.assert_allclose(cal_first.predict(raw_first), cal_second.predict(raw_second), rtol=0, atol=1e-12)
        self.assertEqual(first.n_iter, second.n_iter)
        self.assertEqual(first.converged, second.converged)
        self.assertTrue(first.convergence_applicable)

    def test_stage1_threshold_enforces_each_training_source(self) -> None:
        score = np.array([0.95, 0.92, 0.91, 0.90, 0.10, 0.20])
        source = np.array(["A", "A", "B", "B", "NK1", "NK2"], dtype=object)
        gdt = np.array([True, False, True, False, False, False])
        abt = np.array([False, True, False, True, False, False])
        nk = np.array([False, False, False, False, True, True])
        decision = select_stage1_threshold(score, source, gdt, abt, nk, ["A", "B"], 0.99, 0.98)
        self.assertTrue(decision.passed)
        self.assertEqual(decision.threshold, 0.90)
        self.assertEqual(decision.strict_nk_fpr, 0.0)

    def test_stage2_threshold_passes_and_fails_without_relaxation(self) -> None:
        score = np.array([0.95, 0.90, 0.10, 0.05, 0.96, 0.91, 0.09, 0.04])
        y = np.array([1, 1, 0, 0, 1, 1, 0, 0], dtype=np.int8)
        source = np.array(["A"] * 4 + ["B"] * 4, dtype=object)
        paired = y == 0
        p1 = np.ones(score.size)
        exclusion = np.zeros(score.size, dtype=bool)
        mode = self.config["operating_modes"]["balanced"]
        passed = select_stage2_threshold(
            score, y, source, paired, np.array([0.02, 0.03]), p1, np.ones(2), 0.5, exclusion, np.zeros(2, dtype=bool), "balanced", mode
        )
        self.assertTrue(passed.passed)
        self.assertEqual(passed.paired_abt_fpr, 0.0)
        failed = select_stage2_threshold(
            score, y, source, paired, np.array([0.99, 0.99]), p1, np.ones(2), 0.5, exclusion, np.zeros(2, dtype=bool), "balanced", mode
        )
        self.assertFalse(failed.passed)
        self.assertTrue(np.isnan(failed.threshold))

    def test_rna_call_and_vdj_rescue_are_separate(self) -> None:
        call = apply_two_stage_call(
            np.array([0.9, 0.9, 0.4]),
            np.array([0.9, 0.9, 0.9]),
            0.5,
            0.5,
            np.array([False, True, False]),
        )
        rescue = vdj_rescue_mask(
            call,
            np.array([False, True, False]),
            np.array([False, True, True]),
            np.array([False, False, False]),
        )
        self.assertEqual(call.tolist(), [True, False, False])
        self.assertEqual(rescue.tolist(), [False, True, False])

    def test_inference_coverage_requires_all_critical_genes(self) -> None:
        retained = self.features
        available = set(retained)
        critical = ["CD3D", "TRDC"]
        result = inference_coverage(available, retained, critical, 0.9)
        self.assertFalse(result["abstain"])
        available.remove("TRDC")
        result = inference_coverage(available, retained, critical, 0.9)
        self.assertTrue(result["abstain"])
        self.assertEqual(result["missing_critical_genes"], ["TRDC"])

    def test_hierarchical_bootstrap_is_paired_and_deterministic(self) -> None:
        y = np.array([1, 0, 1, 0, 1, 0, 1, 0], dtype=np.int8)
        candidate = np.array([1, 0, 1, 0, 1, 0, 1, 0], dtype=bool)
        comparator = np.array([1, 1, 0, 0, 1, 1, 0, 0], dtype=bool)
        source = np.array(["A"] * 4 + ["B"] * 4, dtype=object)
        group = np.array(["A1", "A1", "A2", "A2", "B1", "B1", "B2", "B2"], dtype=object)
        counts = grouped_confusion_counts(y, {"candidate": candidate, "comparator": comparator}, source, group)
        first = paired_hierarchical_bootstrap_f1_difference(counts, "candidate", "comparator", 100, 9)
        second = paired_hierarchical_bootstrap_f1_difference(counts, "candidate", "comparator", 100, 9)
        np.testing.assert_array_equal(first, second)
        self.assertGreater(bootstrap_interval(first)[0], 0)

    def test_sparse_extraction_preserves_output_gaps_and_normalizes(self) -> None:
        matrix = sparse.csr_matrix(np.array([[1, 0, 3], [0, 2, 0]], dtype=np.float32))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.h5"
            with h5py.File(path, "w") as handle:
                group = handle.create_group("X")
                group.create_dataset("data", data=matrix.data)
                group.create_dataset("indices", data=matrix.indices)
                group.create_dataset("indptr", data=matrix.indptr)
            with h5py.File(path, "r") as handle:
                extracted = RUNNER.extract_csr_rows(
                    handle["X"],
                    np.array([0, 1]),
                    np.array([0, -1, 3]),
                    5,
                    "raw_counts",
                    10000.0,
                    100,
                )
        self.assertEqual(extracted.shape, (2, 5))
        self.assertEqual(extracted[0, 1], 0.0)
        self.assertEqual(extracted[0, 4], 0.0)
        self.assertAlmostEqual(extracted[0, 0], np.log1p(2500), places=5)
        self.assertAlmostEqual(extracted[0, 3], np.log1p(7500), places=5)

    def test_project_stages_require_checksum_bound_approval(self) -> None:
        config = dict(self.config)
        with tempfile.TemporaryDirectory() as directory:
            config["step2_approval_file"] = str(Path(directory) / "missing.json")
            with self.assertRaises(PermissionError):
                RUNNER.validate_approval(config)

    def test_approval_template_is_bound_to_current_implementation(self) -> None:
        template = json.loads(
            (ROOT / "configs/models/gdtai/v4_step2_approval_template.json").read_text()
        )
        expected = {
            "evaluation_config_sha256": RUNNER.sha256_file(
                ROOT / "configs/models/gdtai/v4_nested_evaluation.json"
            ),
            "runner_sha256": RUNNER.sha256_file(SCRIPT),
            "core_sha256": RUNNER.sha256_file(
                ROOT / "workflows/gdtai/gdtai_v4_nested_core.py"
            ),
        }
        for key, value in expected.items():
            self.assertEqual(template[key], value)

    def test_print_safe_report_renderer_uses_completed_tables(self) -> None:
        config = copy.deepcopy(self.config)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            table_dir = root / "tables"
            figure_dir = root / "figures"
            log_dir = root / "logs"
            static_dir = root / "static"
            for path in [table_dir, figure_dir, log_dir, static_dir]:
                path.mkdir(parents=True)
            config["outputs"].update(
                {
                    "table_dir": str(table_dir),
                    "figure_dir": str(figure_dir),
                    "log_dir": str(log_dir),
                    "static_dir": str(static_dir),
                }
            )
            sources = ["HRA005041", "GSE144469", "BALF_BLOOD_COPD"]
            models = [
                "legacy_trd_minus_trab",
                "compact_7gene_logistic",
                "v2_like_tcr_logistic",
                "v4_individual_gene_ablation",
                "v4_nested_selected",
            ]
            metric_rows = []
            control_rows = []
            for source in sources:
                for model in models:
                    metric_rows.append(
                        {
                            "heldout_source": source,
                            "model_id": model,
                            "mode": "balanced",
                            "n_positive": 100,
                            "n_negative": 1000,
                            "tp": 90,
                            "fp": 1,
                            "fn": 10,
                            "precision": 0.989,
                            "recall": 0.9,
                            "specificity": 0.999,
                            "f1": 0.942,
                            "mcc": 0.93,
                            "pr_auc": 0.98,
                            "roc_auc": 0.99,
                        }
                    )
                    control_rows.append(
                        {
                            "heldout_source": source,
                            "model_id": model,
                            "mode": "balanced",
                            "evaluation": "strict_nk",
                            "n_cells": 1000,
                            "predicted_positive": 1,
                            "fpr": 0.001,
                        }
                    )
            pd.DataFrame(metric_rows).to_csv(table_dir / "outer_fold_metrics.csv", index=False)
            pd.DataFrame(control_rows).to_csv(table_dir / "outer_fold_negative_controls.csv", index=False)
            pd.DataFrame(
                [
                    {
                        "heldout_source": source,
                        "gdt_recall": 0.995,
                        "abt_recall": 0.99,
                        "pass": True,
                    }
                    for source in sources
                ]
            ).to_csv(table_dir / "stage1_heldout_guardrails.csv", index=False)
            pd.DataFrame(
                [
                    {
                        "model_id": model,
                        "mode": "balanced",
                        "n_outer_sources": 3,
                        "f1": 0.942,
                        "f0_5": 0.97,
                        "precision": 0.989,
                        "recall": 0.9,
                        "specificity": 0.999,
                        "mcc": 0.93,
                        "pr_auc": 0.98,
                        "roc_auc": 0.99,
                        "brier_score": 0.01,
                        "ece": 0.01,
                    }
                    for model in models
                ]
            ).to_csv(table_dir / "dataset_macro_metrics.csv", index=False)
            pd.DataFrame(
                [
                    {
                        "candidate": "v4_nested_selected",
                        "comparator": model,
                        "replicates": 2000,
                        "mean_f1_difference": 0.02,
                        "ci95_low": 0.01,
                        "ci95_high": 0.03,
                    }
                    for model in models[:3]
                ]
            ).to_csv(table_dir / "paired_hierarchical_bootstrap_f1.csv", index=False)
            pd.DataFrame(
                [
                    {"condition": f"condition_{index}", "pass": index < 8, "observed": 1.0, "required": "fixed"}
                    for index in range(9)
                ]
            ).to_csv(table_dir / "promotion_guardrails_step2.csv", index=False)
            pd.DataFrame(
                [
                    {
                        "feature": f"gene_{index}",
                        "aggregate_normalized_importance": 1.0 / (index + 1),
                        "directionally_stable": True,
                        "maximum_single_dataset_share": 0.4,
                    }
                    for index in range(30)
                ]
            ).to_csv(table_dir / "feature_stability.csv", index=False)
            pd.DataFrame(
                [
                    {
                        "dataset_id": "GSE114724",
                        "outer_fold_id": "outer_0_HRA005041",
                        "mode": "balanced",
                        "stratum": "strict_NK",
                        "n_cells": 1000,
                        "predicted_positive": 1,
                        "fpr": 0.001,
                    }
                ]
            ).to_csv(table_dir / "extension_negative_controls.csv", index=False)
            (log_dir / "nested_evaluation_aggregate.json").write_text(
                json.dumps(
                    {
                        "strongest_fair_comparator": "v2_like_tcr_logistic",
                        "observed_macro_f1_difference": 0.02,
                        "paired_bootstrap_ci95_low": 0.01,
                    }
                )
            )
            paths = RUNNER.render_nested_report(config, no_pdf=True)
            rendered = Path(paths["html"]).read_text()
            self.assertIn("gdTAI V4 Nested Cross-Dataset Evaluation", rendered)
            self.assertIn("Stage-1 held-out recall", rendered)
            self.assertIn("Whole-atlas inference", rendered)

    def test_end_to_end_synthetic_outer_fold_runs_every_comparator(self) -> None:
        config = copy.deepcopy(self.config)
        config["_candidate_jobs"] = 2
        config["_fold_jobs"] = 2
        config["models"]["stage1_elastic_net"].update(
            {"C": [0.3], "l1_ratio": [0.0], "max_iter": 500, "tolerance": 1e-4}
        )
        config["models"]["stage2_elastic_net"].update(
            {"C": [0.3], "l1_ratio": [0.0], "max_iter": 500, "tolerance": 1e-4}
        )
        config["models"]["stage2_hist_gradient_boosting"].update(
            {
                "learning_rate": [0.07],
                "max_leaf_nodes": [7],
                "min_samples_leaf": [2],
                "l2_regularization": [1.0],
                "max_iter": 20,
                "early_stopping": False,
            }
        )
        sources = ["HRA005041", "GSE144469", "BALF_BLOOD_COPD"]
        rows = []
        expression = []
        legacy = []
        lookup = {gene: index for index, gene in enumerate(self.features)}

        def add_cell(source, truth, role1, role2, fold, replicate, nk_stratum=""):
            row_id = len(rows)
            values = np.zeros(len(self.features), dtype=np.float32)
            if role1 == "t_positive":
                for gene in ["CD3D", "CD3E", "CD3G", "CD2", "LCK", "TRAT1"]:
                    values[lookup[gene]] = 2.0
            if truth == "gdT_primary":
                for gene in ["TRDC", "TRDV1", "TRDJ1", "TRGC1", "TRGV9"]:
                    values[lookup[gene]] = 3.0
                score = 2.0
            elif truth == "abT_primary":
                for gene in ["TRAC", "TRBC1", "TRAV1-1", "TRBV1"]:
                    values[lookup[gene]] = 3.0
                score = -2.0
            else:
                for gene in ["NKG7", "GNLY", "KLRD1", "TYROBP", "FCER1G"]:
                    values[lookup[gene]] = 3.0
                score = -3.0
            group = f"{source}::donor::{truth}_{fold}"
            rows.append(
                {
                    "cell_id": f"cell_{row_id}",
                    "source_gse_id": source,
                    "truth_class": truth,
                    "truth_reliability": 1.0 if role2 != "none" else 0.0,
                    "stage1_role": role1,
                    "stage1_weight": 1.0,
                    "stage2_role": role2,
                    "group_key": group,
                    "has_any_ab_tcr": truth == "abT_primary",
                    "has_any_gd_tcr": truth == "gdT_primary",
                    "has_paired_ab_tcr": truth == "abT_primary",
                    "has_paired_gd_tcr": truth == "gdT_primary",
                    "exclusion_union": False,
                    "nk_sampling_stratum": nk_stratum,
                }
            )
            expression.append(values)
            legacy.append(score)

        for source in sources:
            for fold in range(3):
                for replicate in range(4):
                    add_cell(source, "gdT_primary", "t_positive", "positive", fold, replicate)
                for replicate in range(8):
                    add_cell(source, "abT_primary", "t_positive", "negative", fold, replicate)
        for fold in range(3):
            for replicate in range(10):
                add_cell("GSE254249", "unlabeled", "nk_negative", "none", fold, replicate, "representative")

        cells = pd.DataFrame(rows)
        gene_matrix = np.asarray(expression, dtype=np.float32)
        derived, _ = derive_features(
            gene_matrix,
            self.features,
            config["feature_policy"]["family_prefixes"],
            config["feature_policy"]["derived_panels"],
        )
        split_rows = []
        heldout = "HRA005041"
        outer_fold_id = "outer_0_HRA005041"
        for stage in ["stage1", "stage2"]:
            role = cells[f"{stage}_role"]
            eligible = role.isin(["t_positive", "nk_negative"] if stage == "stage1" else ["positive", "negative"])
            for group_key in cells.loc[eligible & cells.source_gse_id.ne(heldout), "group_key"].drop_duplicates():
                fold = int(group_key.rsplit("_", 1)[1])
                split_rows.append(
                    {
                        "outer_fold_id": outer_fold_id,
                        "heldout_source": heldout,
                        "stage": stage,
                        "group_key": group_key,
                        "inner_fold": fold,
                    }
                )
        result = RUNNER.run_outer_fold(
            cells=cells,
            splits=pd.DataFrame(split_rows),
            gene_matrix=gene_matrix,
            derived_matrix=derived,
            legacy_score=np.asarray(legacy, dtype=np.float32),
            feature_names=self.features,
            heldout_source=heldout,
            config=config,
        )
        self.assertIsNotNone(result["stage2_selected"])
        model_ids = {row["model_id"] for row in result["metrics"]}
        self.assertTrue(
            {
                "v4_nested_selected",
                "v4_individual_gene_ablation",
                "compact_7gene_logistic",
                "v2_like_tcr_logistic",
                "legacy_trd_minus_trab",
            }.issubset(model_ids)
        )
        balanced_v4 = [
            row for row in result["metrics"] if row["model_id"] == "v4_nested_selected" and row["mode"] == "balanced"
        ][0]
        self.assertEqual(balanced_v4["f1"], 1.0)
        self.assertIn("all_folds_converged", result["stage1_candidate_table"].columns)


if __name__ == "__main__":
    unittest.main()
