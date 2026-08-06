import copy
import hashlib
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from tnk_atlas.paths import find_project_root
from tnk_atlas.provenance import sha256_file
from workflows.gdt_atlas.analyze_gdt_atlas import (
    benjamini_hochberg,
    build_composition_tables,
    random_effects_dersimonian_laird,
    run_random_effects_meta_analysis,
    run_within_study_contrasts,
    targeted_cohort_mask,
)
from workflows.gdt_atlas.build_gdt_atlas import (
    AtlasConfigError,
    PrerequisiteError,
    append_frozen_engineered_features,
    build_selection_table,
    extract_selected_counts,
    inspect_h5ad_contract,
    is_tcr_gene,
    load_config,
    selected_profile_thresholds,
    validate_approval_marker,
    validate_config,
    validate_frozen_profile_selection,
    validate_passed_holdout_status,
    verify_selection_decision_digest,
)


ROOT = find_project_root(Path(__file__))
CONFIG_PATH = ROOT / "configs" / "gdt_atlas" / "atlas_rebuild.json"


def _json_sha256(payload: object) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


class GdtAtlasWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = load_config(CONFIG_PATH)

    @staticmethod
    def _fake_frozen_selection() -> SimpleNamespace:
        return SimpleNamespace(
            comparison_workflow_sha256="1" * 64,
            selection_decision_sha256="2" * 64,
            selection_decision_file_sha256="3" * 64,
            holdout_status_sha256="4" * 64,
            selected_profile_id="v2_high_purity",
            selected_model_id="gdtai_v2",
            selected_mode="high_purity",
            selected_model_sha256="5" * 64,
            threshold_contract_sha256="6" * 64,
        )

    @staticmethod
    def _signed_decision() -> dict[str, object]:
        decision: dict[str, object] = {
            "schema_version": 1,
            "status": "selected",
            "selected_profile": "v2_high_purity",
            "profile_identity": {
                "v2_high_purity": {
                    "model_id": "gdtai_v2",
                    "mode": "high_purity",
                    "artifact": "/tmp/v2.pkl",
                    "sha256": "7" * 64,
                }
            },
            "selection_ranking": [
                {
                    "profile_id": "v2_high_purity",
                    "eligible": True,
                    "selected": True,
                    "gold_recall_pass": True,
                    "pooled_abt_fpr_pass": True,
                    "pooled_strict_nk_fpr_pass": True,
                    "large_cohort_fpr_pass": True,
                    "labeled_fp_fraction_pass": True,
                }
            ],
        }
        decision["decision_sha256"] = _json_sha256(decision)
        return decision

    def test_config_requires_four_profile_decision_and_scientific_invariants(
        self,
    ) -> None:
        profile = self.config["frozen_profile_selection"]
        selection = self.config["selection"]
        clustering = self.config["embedding_clustering"]
        analysis = self.config["analysis"]

        required_profiles = {
            row["profile_id"]: (row["model_id"], row["mode"])
            for row in profile["required_profiles"]
        }
        self.assertEqual(
            required_profiles,
            {
                "v2_high_f1": ("gdtai_v2", "high_f1"),
                "v2_high_purity": ("gdtai_v2", "high_purity"),
                "v3_round14_balanced": ("gdtai_v3", "fixed_threshold"),
                "v3_round12_high_purity": (
                    "gdtai_v3_round12",
                    "fixed_threshold",
                ),
            },
        )
        for forbidden in ["profile_id", "model_id", "model_sha256", "threshold"]:
            self.assertNotIn(forbidden, profile)
        self.assertEqual(profile["required_selection_status"], "selected")
        self.assertEqual(profile["required_holdout_promotion_status"], "holdout_passed")
        self.assertEqual(
            profile["required_normalization"], "log1p(raw_counts_per_10000)"
        )
        self.assertFalse(selection["allow_silver_only_fn_addback"])
        self.assertFalse(selection["hard_trdv_exclusion"])
        self.assertFalse(selection["nk_gene_only_exclusion"])
        self.assertTrue(selection["known_alpha_beta_only_requires_paired_trab_or_gold"])
        self.assertTrue(
            selection["known_alpha_beta_only_requires_no_productive_gdt_evidence"]
        )
        self.assertGreater(int(selection["minimum_selected_cells"]), 0)
        self.assertTrue(
            analysis["targeted_cohort_exclusion"]["exclude_from_abundance_inference"]
        )
        self.assertTrue(analysis["donor_pseudobulk_de"]["enabled"])
        self.assertGreaterEqual(len(set(clustering["seeds"])), 2)
        self.assertGreaterEqual(len(set(clustering["leiden_resolutions"])), 2)
        self.assertGreaterEqual(
            {
                value.upper()
                for value in self.config["feature_selection"][
                    "tcr_gene_prefixes_excluded"
                ]
            },
            {"TRA", "TRB", "TRG", "TRD"},
        )
        self.assertGreaterEqual(
            {
                value.upper()
                for value in analysis["separate_tcr_expression"]["gene_prefixes"]
            },
            {"TRDV", "TRGV"},
        )

    def test_config_rejects_preselection_and_profile_mode_swaps(self) -> None:
        preselected = copy.deepcopy(self.config)
        preselected["frozen_profile_selection"]["threshold"] = 0.5
        with self.assertRaisesRegex(AtlasConfigError, "must not preselect"):
            validate_config(preselected)

        swapped = copy.deepcopy(self.config)
        swapped["frozen_profile_selection"]["required_profiles"][0]["mode"] = (
            "high_purity"
        )
        with self.assertRaisesRegex(AtlasConfigError, "exact four-profile"):
            validate_config(swapped)

    def test_missing_profile_decision_and_holdout_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            config = copy.deepcopy(self.config)
            with self.assertRaisesRegex(
                PrerequisiteError,
                "Four-profile gdTAI decision prerequisites are absent",
            ):
                validate_frozen_profile_selection(config, root)

    def test_signed_decision_and_passed_holdout_are_both_required(self) -> None:
        decision = self._signed_decision()
        digest = verify_selection_decision_digest(decision)
        self.assertEqual(digest, decision["decision_sha256"])
        holdout = {
            "selected_profile": "v2_high_purity",
            "selection_decision_sha256": digest,
            "holdout_schema_failure": False,
            "holdout_veto": False,
            "promotion_status": "holdout_passed",
        }
        validate_passed_holdout_status(decision, holdout)

        mutated = copy.deepcopy(decision)
        mutated["selected_profile"] = "v2_high_f1"
        with self.assertRaisesRegex(PrerequisiteError, "digest is invalid"):
            verify_selection_decision_digest(mutated)

        failed_guardrail = copy.deepcopy(decision)
        failed_guardrail["selection_ranking"][0]["pooled_strict_nk_fpr_pass"] = False
        failed_guardrail.pop("decision_sha256")
        failed_guardrail["decision_sha256"] = _json_sha256(failed_guardrail)
        with self.assertRaisesRegex(PrerequisiteError, "every T/NK guardrail"):
            verify_selection_decision_digest(failed_guardrail)

        vetoed = dict(holdout)
        vetoed["holdout_veto"] = True
        vetoed["promotion_status"] = "vetoed_no_runner_up_selected"
        with self.assertRaisesRegex(PrerequisiteError, "vetoed"):
            validate_passed_holdout_status(decision, vetoed)

        wrong_decision = dict(holdout)
        wrong_decision["selection_decision_sha256"] = "0" * 64
        with self.assertRaisesRegex(PrerequisiteError, "not linked"):
            validate_passed_holdout_status(decision, wrong_decision)

    def test_approval_marker_pins_input_config_decision_holdout_and_profile(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_path = root / "rebuilt.h5ad"
            input_path.write_bytes(b"synthetic-approved-rebuild")
            config = copy.deepcopy(self.config)
            config["execution"]["rebuilt_integrated_h5ad"] = "rebuilt.h5ad"
            config["execution"]["approval_marker"] = "approved.json"
            config_path = root / "atlas_rebuild.json"
            config_path.write_text(json.dumps(config, sort_keys=True), encoding="utf-8")
            frozen = self._fake_frozen_selection()
            marker = {
                "schema_version": config["prerequisites"]["approval_schema_version"],
                "decision": "approved",
                "purpose": "build_gdt_atlas_de_novo",
                "workflow_id": config["workflow_id"],
                "approved_by": "unit-test-reviewer",
                "approved_at_utc": "2026-08-06T12:00:00Z",
                "input_h5ad": "rebuilt.h5ad",
                "input_sha256": sha256_file(input_path),
                "input_size_bytes": input_path.stat().st_size,
                "config_sha256": sha256_file(config_path),
                "gdtai_comparison_workflow_sha256": frozen.comparison_workflow_sha256,
                "gdtai_selection_decision_sha256": frozen.selection_decision_sha256,
                "gdtai_selection_decision_file_sha256": frozen.selection_decision_file_sha256,
                "gdtai_holdout_status_sha256": frozen.holdout_status_sha256,
                "gdtai_selected_profile_id": frozen.selected_profile_id,
                "gdtai_selected_model_id": frozen.selected_model_id,
                "gdtai_selected_mode": frozen.selected_mode,
                "gdtai_selected_model_sha256": frozen.selected_model_sha256,
                "gdtai_threshold_contract_sha256": frozen.threshold_contract_sha256,
            }
            marker_path = root / "approved.json"
            marker_path.write_text(json.dumps(marker), encoding="utf-8")

            result = validate_approval_marker(config, config_path, frozen, root)
            self.assertEqual(result[0], input_path.resolve())
            self.assertEqual(result[2], marker["input_sha256"])

            marker["gdtai_holdout_status_sha256"] = "9" * 64
            marker_path.write_text(json.dumps(marker), encoding="utf-8")
            with self.assertRaisesRegex(PrerequisiteError, "holdout_status_sha256"):
                validate_approval_marker(config, config_path, frozen, root)
            marker["gdtai_holdout_status_sha256"] = frozen.holdout_status_sha256
            marker_path.write_text(json.dumps(marker), encoding="utf-8")

            input_path.write_bytes(b"changed-after-approval")
            with self.assertRaisesRegex(
                PrerequisiteError, "size no longer matches|SHA256 mismatch"
            ):
                validate_approval_marker(config, config_path, frozen, root)

    def test_h5ad_contract_requires_csr_integer_counts_and_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "rebuilt.h5ad"
            n_obs = int(self.config["prerequisites"]["minimum_cells"])
            n_vars = int(self.config["prerequisites"]["minimum_genes"])
            with h5py.File(path, "w") as handle:
                obs = handle.create_group("obs")
                obs.create_dataset(
                    "_index",
                    data=np.asarray([f"cell_{i}" for i in range(n_obs)], dtype="S"),
                )
                for column in self.config["selection"]["required_obs_columns"]:
                    obs.create_dataset(column, data=np.zeros(n_obs, dtype=np.int8))
                var = handle.create_group("var")
                var.create_dataset(
                    "_index",
                    data=np.asarray([f"GENE{i}" for i in range(n_vars)], dtype="S"),
                )
                counts = handle.create_group("layers").create_group("counts")
                counts.attrs["encoding-type"] = "csr_matrix"
                counts.attrs["shape"] = np.asarray([n_obs, n_vars], dtype=np.int64)
                counts.create_dataset("data", data=np.ones(n_obs, dtype=np.float32))
                counts.create_dataset("indices", data=np.zeros(n_obs, dtype=np.int32))
                counts.create_dataset(
                    "indptr", data=np.arange(n_obs + 1, dtype=np.int64)
                )

            contract = inspect_h5ad_contract(path, self.config)
            self.assertEqual(contract["counts_encoding"], "csr_matrix")
            self.assertEqual(contract["n_obs"], n_obs)

            with h5py.File(path, "r+") as handle:
                handle["layers"]["counts"]["data"][0] = 0.5
            with self.assertRaisesRegex(PrerequisiteError, "not integer-like"):
                inspect_h5ad_contract(path, self.config)

    def test_blockwise_csr_subset_extraction_preserves_rows(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "counts.h5ad"
            dense = np.arange(100, dtype=np.int32).reshape(20, 5)
            dense[dense % 3 != 0] = 0
            csr = sparse.csr_matrix(dense)
            with h5py.File(path, "w") as handle:
                counts = handle.create_group("layers").create_group("counts")
                counts.attrs["encoding-type"] = "csr_matrix"
                counts.attrs["shape"] = np.asarray(csr.shape, dtype=np.int64)
                counts.create_dataset("data", data=csr.data)
                counts.create_dataset("indices", data=csr.indices)
                counts.create_dataset("indptr", data=csr.indptr)
            selected_rows = np.asarray([0, 3, 4, 11, 19], dtype=np.int64)
            observed = extract_selected_counts(
                path, "counts", selected_rows, block_rows=4
            ).toarray()
            np.testing.assert_array_equal(observed, dense[selected_rows])

    def _selection_obs(self) -> pd.DataFrame:
        names = [
            "predicted_clean",
            "predicted_ab_only",
            "predicted_nk",
            "predicted_doublet",
            "predicted_non_t",
            "predicted_bad_qc",
            "gold_fn",
            "silver_fn",
            "gold_fn_bad_qc",
            "predicted_dual_tcr",
            "predicted_any_ab_unpaired",
            "predicted_ab_gold",
            "gold_fn_missing_trd",
            "predicted_nk_gene_only",
        ]
        frame = pd.DataFrame(index=names)
        frame["source_gse_id"] = "GSE_TEST"
        frame["sample_id"] = "sample"
        frame["library_id"] = "library"
        frame["donor_id"] = "donor"
        frame["condition"] = "control"
        frame["tissue_corrected"] = "blood"
        frame["technology_simple"] = "10x_5prime"
        frame["atlas_targeted_cohort"] = False
        frame["has_TRA_TRB_paired"] = False
        frame["has_any_ab_tcr"] = False
        frame["has_TRG_TRD_paired"] = False
        frame["TRG_cdr3"] = ""
        frame["TRD_cdr3"] = ""
        frame["gdt_truth_class"] = "unlabeled_or_ambiguous"
        frame["atlas_high_confidence_nk"] = False
        frame["atlas_doublet"] = False
        frame["atlas_non_t_contaminant"] = False
        frame["atlas_severe_qc_failure"] = False
        frame["NKG7"] = 0.0

        frame.loc["predicted_ab_only", ["has_TRA_TRB_paired", "has_any_ab_tcr"]] = True
        frame.loc["predicted_nk", "atlas_high_confidence_nk"] = True
        frame.loc["predicted_doublet", "atlas_doublet"] = True
        frame.loc["predicted_non_t", "atlas_non_t_contaminant"] = True
        frame.loc["predicted_bad_qc", "atlas_severe_qc_failure"] = True
        for name in ["gold_fn", "gold_fn_bad_qc", "gold_fn_missing_trd"]:
            frame.loc[name, "gdt_truth_class"] = "gdt_gold"
            frame.loc[name, "has_TRG_TRD_paired"] = True
            frame.loc[name, "TRG_cdr3"] = "CAVGDT"
            frame.loc[name, "TRD_cdr3"] = "CALGDT"
        frame.loc["gold_fn_bad_qc", "atlas_severe_qc_failure"] = True
        frame.loc["gold_fn_missing_trd", "TRD_cdr3"] = "NA"
        frame.loc["silver_fn", "gdt_truth_class"] = "gdt_silver"
        frame.loc["silver_fn", "TRD_cdr3"] = "CALSILVER"
        frame.loc[
            "predicted_dual_tcr",
            ["has_TRA_TRB_paired", "has_any_ab_tcr", "has_TRG_TRD_paired"],
        ] = True
        frame.loc["predicted_dual_tcr", ["TRG_cdr3", "TRD_cdr3"]] = [
            "CAVDUAL",
            "CALDUAL",
        ]
        frame.loc["predicted_any_ab_unpaired", "has_any_ab_tcr"] = True
        frame.loc["predicted_ab_gold", "gdt_truth_class"] = "abt_gold"
        frame.loc["predicted_nk_gene_only", "NKG7"] = 8.0
        return frame

    def test_selection_uses_only_allowed_evidence(self) -> None:
        obs = self._selection_obs()
        scores = np.full(obs.shape[0], 0.99, dtype=np.float32)
        thresholds = np.full(obs.shape[0], 0.936, dtype=np.float32)
        for name in ["gold_fn", "silver_fn", "gold_fn_bad_qc", "gold_fn_missing_trd"]:
            scores[obs.index.get_loc(name)] = 0.10
        selected = build_selection_table(obs, scores, thresholds, self.config)

        expected_selected = {
            "predicted_clean",
            "gold_fn",
            "predicted_dual_tcr",
            "predicted_any_ab_unpaired",
            "predicted_nk_gene_only",
        }
        observed_selected = set(selected.index[selected["gdt_atlas_selected"]])
        self.assertEqual(observed_selected, expected_selected)
        self.assertTrue(selected.loc["gold_fn", "selection_gold_fn_added"])
        self.assertFalse(selected.loc["silver_fn", "gdt_atlas_selected"])
        self.assertFalse(selected.loc["gold_fn_bad_qc", "gdt_atlas_selected"])
        self.assertFalse(selected.loc["gold_fn_missing_trd", "selection_gold_fn_added"])
        self.assertFalse(
            selected.loc["predicted_dual_tcr", "selection_known_alpha_beta_only"]
        )
        self.assertFalse(
            selected.loc["predicted_any_ab_unpaired", "selection_known_alpha_beta_only"]
        )
        self.assertTrue(
            selected.loc["predicted_ab_gold", "selection_known_alpha_beta_only"]
        )
        self.assertTrue(selected.loc["predicted_nk_gene_only", "gdt_atlas_selected"])

    def test_selection_rejects_ambiguous_upstream_flags(self) -> None:
        obs = self._selection_obs()
        obs["atlas_doublet"] = obs["atlas_doublet"].astype(object)
        obs.loc["predicted_clean", "atlas_doublet"] = "maybe"
        with self.assertRaisesRegex(PrerequisiteError, "invalid"):
            build_selection_table(
                obs,
                np.full(obs.shape[0], 0.99),
                np.full(obs.shape[0], 0.5),
                self.config,
            )

    def test_selection_uses_per_cell_selected_profile_thresholds(self) -> None:
        obs = self._selection_obs().iloc[:2].copy()
        scores = np.asarray([0.90, 0.90], dtype=np.float32)
        thresholds = np.asarray([0.80, np.inf], dtype=np.float32)
        selected = build_selection_table(obs, scores, thresholds, self.config)
        self.assertTrue(selected.iloc[0]["gdtai_frozen_prediction"])
        self.assertFalse(selected.iloc[1]["gdtai_frozen_prediction"])
        self.assertTrue(np.isposinf(selected.iloc[1]["gdtai_frozen_threshold"]))

    def test_selected_threshold_dispatch_supports_v2_and_v3(self) -> None:
        annotations = np.asarray(["GDT_CELL", "NK_CELL"], dtype=object)
        evaluator = SimpleNamespace(
            v2_thresholds=lambda _profile, _annotations: np.asarray(
                [0.80, np.inf], dtype=np.float32
            )
        )
        v2 = SimpleNamespace(model_id="gdtai_v2", payload={})
        np.testing.assert_array_equal(
            selected_profile_thresholds(evaluator, v2, annotations),
            np.asarray([0.80, np.inf], dtype=np.float32),
        )

        v3 = SimpleNamespace(model_id="gdtai_v3", payload={"threshold": 0.73})
        np.testing.assert_allclose(
            selected_profile_thresholds(evaluator, v3, annotations),
            np.asarray([0.73, 0.73], dtype=np.float32),
        )

    def test_tcr_gene_exclusion_and_separate_v_gene_scope(self) -> None:
        for gene in ["TRAC", "TRBV7-9", "TRGV9", "TRDV2", "TRDC"]:
            self.assertTrue(is_tcr_gene(gene, ["TRA", "TRB", "TRG", "TRD"]))
        self.assertFalse(is_tcr_gene("CD3D", ["TRA", "TRB", "TRG", "TRD"]))

    def test_frozen_engineered_feature_order(self) -> None:
        genes = [
            "TRDC",
            "TRDV2",
            "TRDJ1",
            "TRGV9",
            "TRAC",
            "CD3D",
            "CD3E",
            "NKG7",
            "GNLY",
        ]
        x_gene = np.asarray(
            [[2.0, 1.0, 0.5, 3.0, 4.0, 1.0, 3.0, 6.0, 2.0]], dtype=np.float32
        )
        engineered = ["any_TRDV", "CD3_score", "NK_score", "TRDC_log1p"]
        output = append_frozen_engineered_features(x_gene, genes, engineered)
        appended = output[0, len(genes) :]
        self.assertEqual(float(appended[0]), 1.0)
        self.assertAlmostEqual(float(appended[1]), 2.0)
        self.assertAlmostEqual(float(appended[2]), 4.0)
        self.assertAlmostEqual(float(appended[3]), 2.0)

    def test_targeted_cohorts_are_descriptive_only_for_abundance(self) -> None:
        rows = []
        units = [
            ("GSE_UNTARGETED", "D1", "S1", False),
            ("GDT_2020AUG_woCOV", "D2", "S2", False),
            ("GSE_FLAGGED", "D3", "S3", True),
        ]
        for source, donor, sample, targeted in units:
            for cell in range(30):
                rows.append(
                    {
                        "source_gse_id": source,
                        "sample_id": sample,
                        "library_id": sample,
                        "donor_id": donor,
                        "condition": "control",
                        "tissue_corrected": "blood",
                        "gdt_leiden": "0"
                        if source == "GSE_UNTARGETED"
                        else str(cell % 2),
                        "atlas_targeted_cohort": targeted,
                    }
                )
        obs = pd.DataFrame(rows)
        mask = targeted_cohort_mask(obs, self.config)
        self.assertEqual(int(mask.sum()), 60)
        sample, donor, _, audit = build_composition_tables(obs, self.config)
        eligible_samples = set(
            sample.loc[sample["eligible_for_abundance_inference"], "source_gse_id"]
        )
        eligible_donors = set(
            donor.loc[donor["eligible_for_abundance_inference"], "source_gse_id"]
        )
        self.assertEqual(eligible_samples, {"GSE_UNTARGETED"})
        self.assertEqual(eligible_donors, {"GSE_UNTARGETED"})
        zero_cluster = sample[
            sample["source_gse_id"].eq("GSE_UNTARGETED") & sample["gdt_leiden"].eq("1")
        ]
        self.assertEqual(int(zero_cluster["cluster_cells"].iloc[0]), 0)
        self.assertEqual(
            int(audit.loc[audit["unit_type"].eq("sample"), "n_targeted_units"].iloc[0]),
            2,
        )

    def test_random_effects_and_fdr_helpers(self) -> None:
        result = random_effects_dersimonian_laird([0.2, 0.4, 0.3], [0.04, 0.05, 0.03])
        self.assertEqual(result["k_studies"], 3)
        self.assertGreater(result["pooled_effect"], 0.2)
        self.assertLess(result["pooled_effect"], 0.4)
        adjusted = benjamini_hochberg([0.01, 0.04, 0.03, np.nan])
        self.assertTrue(np.isnan(adjusted[-1]))
        self.assertTrue(np.all((adjusted[:3] >= 0) & (adjusted[:3] <= 1)))
        self.assertGreaterEqual(adjusted[0], 0.01)

    def test_configured_within_study_contrasts_feed_meta_analysis(self) -> None:
        rows = []
        for study in ["STUDY_A", "STUDY_B"]:
            for group, cluster_one_cells in [("control", 8), ("case", 22)]:
                for donor_index in range(3):
                    donor_cluster_one_cells = cluster_one_cells + donor_index - 1
                    for cell in range(30):
                        rows.append(
                            {
                                "source_gse_id": study,
                                "sample_id": f"{group}_{donor_index}",
                                "library_id": f"{group}_{donor_index}",
                                "donor_id": f"{group}_{donor_index}",
                                "condition": group,
                                "tissue_corrected": "blood",
                                "gdt_leiden": "1"
                                if cell < donor_cluster_one_cells
                                else "0",
                                "atlas_targeted_cohort": False,
                            }
                        )
        config = copy.deepcopy(self.config)
        config["analysis"]["within_study_contrasts"] = [
            {
                "contrast_id": "case_vs_control",
                "enabled": True,
                "study_ids": [],
                "group_column": "condition",
                "case_values": ["case"],
                "reference_values": ["control"],
                "paired_by_donor": False,
                "note": "synthetic test",
            }
        ]
        _, _, donor_condition, _ = build_composition_tables(pd.DataFrame(rows), config)
        study_results, plan = run_within_study_contrasts(donor_condition, config)
        meta = run_random_effects_meta_analysis(study_results, config)

        self.assertEqual(study_results.shape[0], 4)
        self.assertEqual(meta.shape[0], 2)
        self.assertTrue(plan.loc[0, "enabled"])
        cluster_one = meta[meta["gdt_leiden"].eq("1")].iloc[0]
        self.assertEqual(int(cluster_one["k_studies"]), 2)
        self.assertGreater(float(cluster_one["pooled_effect"]), 0.0)


if __name__ == "__main__":
    unittest.main()
