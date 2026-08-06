import importlib.util
import json
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd
import h5py

ROOT_HINT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT_HINT / "src"))

ROOT = ROOT_HINT
SCRIPT = ROOT / "workflows/gdtai/compare_frozen_gdtai_profiles.py"
SPEC = importlib.util.spec_from_file_location("compare_frozen_gdtai_profiles", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class FrozenProfileEvaluationTests(unittest.TestCase):
    def test_config_has_four_profiles_and_two_sealed_holdouts(self) -> None:
        config = json.loads((ROOT / "configs/models/gdtai/extension_evaluation.json").read_text())
        self.assertEqual(len(config["profiles"]), 4)
        self.assertEqual(set(config["sealed_holdouts"]), {"GSE315928", "GSE121636"})

    def test_truth_uses_single_alpha_or_beta_chain_as_negative(self) -> None:
        obs = {
            "has_any_ab_tcr": np.array([True, True, False]),
            "has_TRA_TRB_paired": np.array([False, True, False]),
            "has_any_gd_tcr": np.array([False, False, True]),
            "has_TRG_TRD_paired": np.array([False, False, True]),
            "is_doublet": np.zeros(3, dtype=bool),
        }
        truth = MODULE.truth_frame(obs, np.array(["CD8_T", "CD8_T", "GDT_CELL"], dtype=object))
        self.assertEqual(truth["truth_class"].tolist(), ["abT_gold", "abT_gold", "gdT_gold"])

    def test_unrecognized_author_label_does_not_suppress_nk_fallback(self) -> None:
        genes = ["CD3D", "CD3E", "CD3G", "TYROBP", "FCER1G", "KLRD1", "NCAM1"]
        spec = MODULE.FeatureSpec(
            gene_names=genes,
            gene_indices=np.arange(len(genes), dtype=np.int32),
            gene_feature_names=[f"{gene}_log1p_cp10k" for gene in genes],
            engineered_feature_names=[],
            model_feature_names=[f"{gene}_log1p_cp10k" for gene in genes],
            gene_to_col={gene: index for index, gene in enumerate(genes)},
            engineered_to_col={},
        )
        expression = np.array(
            [
                [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0],
            ],
            dtype=np.float32,
        )
        annotation, source = MODULE.normalize_annotation(
            {"author_annotation": np.array(["immune cell", "NK cell"], dtype=object)},
            expression,
            spec,
        )
        self.assertEqual(annotation.tolist(), ["NK_CELL", "NK_CELL"])
        self.assertEqual(source.tolist(), ["derived_expression", "author_annotation"])
        truth = MODULE.truth_frame({}, annotation, source)
        self.assertEqual(truth["strict_NK"].tolist(), [True, True])
        self.assertEqual(truth["author_NK"].tolist(), [False, True])

    def test_silver_is_derived_only_inside_gdt_sequenced_sublibrary(self) -> None:
        obs = {
            "library_id": np.array(["gd_lib", "gd_lib", "ab_lib", "ab_lib"], dtype=object),
            "has_TRA": np.array([False, False, True, False]),
            "has_TRB": np.array([False, False, False, False]),
            "has_TRG": np.array([True, False, False, False]),
            "has_TRD": np.array([True, True, False, False]),
        }
        silver = MODULE.derive_silver_mask(obs)
        self.assertEqual(silver.tolist(), [False, True, False, False])
        rows = MODULE.select_input_rows(obs, "labeled_with_silver")
        self.assertEqual(rows.tolist(), [0, 1, 2])

    def test_cdr3_fields_restore_chain_evidence_without_boolean_flags(self) -> None:
        obs = {
            "library_id": np.array(["gd_lib", "gd_lib", "ab_lib"], dtype=object),
            "TRA_cdr3": np.array(["", "", "CAVR"], dtype=object),
            "TRB_cdr3": np.array(["", "", ""], dtype=object),
            "TRG_cdr3": np.array(["CALW", "", ""], dtype=object),
            "TRD_cdr3": np.array(["CACD", "CACX", ""], dtype=object),
        }
        truth = MODULE.truth_frame(obs, np.full(3, "OTHER", dtype=object))
        self.assertEqual(truth["truth_class"].tolist(), ["gdT_gold", "gdT_silver", "abT_gold"])

    def test_evaluation_filter_includes_author_nk_but_not_nkt(self) -> None:
        obs = {
            "cell_type": np.array(["NK cell", "NKT", "CD8 T", "gdT"], dtype=object),
            "has_TRA": np.array([False, False, True, False]),
            "has_TRB": np.array([False, False, False, False]),
            "has_TRG": np.array([False, False, False, True]),
            "has_TRD": np.array([False, False, False, True]),
            "library_id": np.array(["a", "a", "a", "a"], dtype=object),
        }
        rows = MODULE.select_input_rows(obs, "evaluation_labeled")
        self.assertEqual(rows.tolist(), [0, 2, 3])

    def test_selection_requires_every_guardrail(self) -> None:
        config = json.loads((ROOT / "configs/models/gdtai/extension_evaluation.json").read_text())
        metrics = pd.DataFrame([
            {"profile_id": "good", "evaluation": "primary_gold", "f1": 0.9, "precision": 0.95, "recall": 0.85, "predicted_positive": 100, "fp": 4},
            {"profile_id": "bad_nk", "evaluation": "primary_gold", "f1": 0.95, "precision": 0.95, "recall": 0.9, "predicted_positive": 100, "fp": 4},
        ])
        strata = pd.DataFrame([
            {"profile_id": "good", "dataset_id": "A", "stratum": "abT_gold", "n_cells": 10000, "predicted_gdT": 5, "call_rate": 0.0005},
            {"profile_id": "good", "dataset_id": "A", "stratum": "strict_NK", "n_cells": 10000, "predicted_gdT": 5, "call_rate": 0.0005},
            {"profile_id": "bad_nk", "dataset_id": "A", "stratum": "abT_gold", "n_cells": 10000, "predicted_gdT": 5, "call_rate": 0.0005},
            {"profile_id": "bad_nk", "dataset_id": "A", "stratum": "strict_NK", "n_cells": 10000, "predicted_gdT": 20, "call_rate": 0.002},
        ])
        ranking = MODULE.select_profile(metrics, strata, config)
        self.assertEqual(ranking.loc[ranking["selected"], "profile_id"].tolist(), ["good"])
        self.assertFalse(bool(ranking.loc[ranking["profile_id"] == "bad_nk", "eligible"].iloc[0]))

    def test_decision_digest_detects_mutation(self) -> None:
        payload = {"a": 1, "b": [2, 3]}
        digest = MODULE.sha256_json(payload)
        self.assertEqual(digest, MODULE.sha256_json(payload))
        self.assertNotEqual(digest, MODULE.sha256_json({"a": 2, "b": [2, 3]}))

    def test_gse144469_legacy_tcr_flags_restore_both_classes(self) -> None:
        path = ROOT / "data/datasets/GSE144469/processed/artifacts/GSE144469_tnk_subset.h5ad"
        if not path.exists():
            self.skipTest("GSE144469 standalone T/NK H5AD is not available")
        with h5py.File(path, "r") as handle:
            obs = {column: MODULE.read_obs(handle, column) for column in MODULE.CANONICAL_OBS_COLUMNS}
        rows = MODULE.select_input_rows(obs, "primary_labeled")
        selected = {column: values[rows] for column, values in obs.items()}
        truth = MODULE.truth_frame(selected, np.full(rows.size, "OTHER", dtype=object))
        self.assertGreater(int(truth["gdT_gold"].sum()), 0)
        self.assertGreater(int(truth["abT_gold"].sum()), 0)
        self.assertEqual(int((truth["gdT_gold"] & truth["abT_gold"]).sum()), 0)


if __name__ == "__main__":
    unittest.main()
