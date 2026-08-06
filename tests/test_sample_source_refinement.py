import json
from pathlib import Path
import tempfile
import unittest

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from workflows.metadata.sample_source_refinement_workflow import (
    RefinementError,
    SampleSourceRefiner,
    file_sha256,
    load_rules,
    review_h5ad,
    writeback_h5ad,
)


ROOT = Path(__file__).resolve().parents[1]
RULES_PATH = ROOT / "configs/metadata/sample_source_refinement_rules.json"


class SampleSourceRuleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.rules = load_rules(RULES_PATH)
        cls.refiner = SampleSourceRefiner(cls.rules)

    def assert_resolution(
        self,
        row: dict[str, str],
        value: str,
        level: str,
        rule_id: str,
    ) -> None:
        resolution = self.refiner.resolve(row)
        self.assertEqual(resolution.value, value)
        self.assertEqual(resolution.level, level)
        self.assertEqual(resolution.rule_id, rule_id)

    def test_required_precedence_is_fixed(self) -> None:
        self.assertEqual(
            self.rules["precedence"],
            [
                "cell_level",
                "sample_library",
                "gse_tissue",
                "global_tissue_fallback",
                "unresolved",
            ],
        )
        self.assertEqual(self.rules["output_field"], "sample_source_refined")

    def test_cell_level_beats_exact_sample_and_gse_tissue(self) -> None:
        self.assert_resolution(
            {
                "source_gse_id": "GSE243013",
                "cell_source": "LUAD",
                "sample_id": "P1",
                "tissue": "tumor",
            },
            "LUAD_tumor",
            "cell_level",
            "cell_level_explicit_source",
        )

    def test_gse243013_uses_exact_nonoverlapping_luad_lusc_map(self) -> None:
        sample_map = self.rules["sample_value_maps"][0]
        groups = {group["value"]: set(group["keys"]) for group in sample_map["value_groups"]}
        self.assertEqual(len(groups["LUAD_tumor"]), 63)
        self.assertEqual(len(groups["LUSC_tumor"]), 180)
        self.assertFalse(groups["LUAD_tumor"] & groups["LUSC_tumor"])

        self.assert_resolution(
            {"source_gse_id": "GSE243013", "sample_id": "P5", "tissue": "tumor"},
            "LUAD_tumor",
            "sample_library",
            "gse243013_exact_luad",
        )
        self.assert_resolution(
            {"source_gse_id": "GSE243013", "sample_id": "P1", "tissue": "tumor"},
            "LUSC_tumor",
            "sample_library",
            "gse243013_exact_lusc",
        )
        self.assert_resolution(
            {"source_gse_id": "GSE243013", "sample_id": "P999", "tissue": "tumor"},
            "unresolved",
            "sample_library",
            "gse243013_unmapped_sample_unresolved",
        )

    def test_known_gse_tissue_mappings_are_tissue_gated(self) -> None:
        cases = [
            (
                {"source_gse_id": "GSE162498", "tissue": "tumor"},
                "NSCLC_tumor",
                "gse162498_nsclc_tumor",
            ),
            (
                {"source_gse_id": "GSE235863", "tissue": "liver_tumor"},
                "HBV_positive_HCC_tumor",
                "gse235863_hbv_positive_hcc_tumor",
            ),
            (
                {"source_gse_id": "GSE190870", "tissue": "lymph_node_metastasis"},
                "breast_invasive_ductal_carcinoma_lymph_node_metastasis",
                "gse190870_breast_idc_lymph_node_metastasis",
            ),
            (
                {"source_gse_id": "GSE287301", "sample_id": "chip1pool1", "tissue": "tumor"},
                "HNSCC_tumor",
                "gse287301_hnscc_tumor",
            ),
        ]
        for row, value, rule_id in cases:
            with self.subTest(rule_id=rule_id):
                self.assert_resolution(row, value, "gse_tissue", rule_id)

        self.assert_resolution(
            {"source_gse_id": "GSE162498", "tissue": "blood"},
            "NSCLC_blood",
            "gse_tissue",
            "gse162498_nsclc_blood",
        )

    def test_gse287301_patient_exception_preserves_known_tumor_context(self) -> None:
        for row in [
            {"source_gse_id": "GSE287301", "patient_id": "10", "tissue": "tumor"},
            {
                "source_gse_id": "GSE287301",
                "patient_id": "3;10;18;27",
                "sample_id": "chip1pool3",
                "tissue": "tumor",
            },
            {
                "source_gse_id": "GSE287301",
                "library_id": "GSM8743452",
                "tissue": "tumor",
            },
        ]:
            with self.subTest(row=row):
                self.assert_resolution(
                    row,
                    "tumor_pool_HNSCC_with_ameloblastoma_ambiguity",
                    "sample_library",
                    "gse287301_patient_10_mixed_tumor_pool",
                )

    def test_extension_tumor_projects_have_disease_aware_specimen_context(self) -> None:
        cases = [
            ("GSE114724", "breast_tumor", "breast_cancer_primary_tumor"),
            ("GSE121636", "renal_tumor", "clear_cell_RCC_primary_tumor"),
            ("GSE121636", "blood", "clear_cell_RCC_blood"),
            ("GSE159251", "blood", "metastatic_melanoma_blood"),
            ("GSE159251", "lymph_node", "melanoma_lymph_node_metastasis"),
            ("GSE159251", "subcutaneous_metastasis", "melanoma_subcutaneous_metastasis"),
            ("GSE292700", "tumor", "LUAD_primary_tumor"),
            ("GSE292700", "adjacent_tissue", "LUAD_adjacent_non_tumor_tissue"),
            ("GSE292700", "blood", "LUAD_blood"),
            ("GSE294273", "lymph_node_metastasis", "melanoma_lymph_node_metastasis"),
            ("GSE296954", "tumor", "HNSCC_primary_tumor"),
            ("GSE315928", "tumor", "gastric_cancer_primary_tumor"),
        ]
        for gse, tissue, expected in cases:
            with self.subTest(gse=gse, tissue=tissue):
                resolution = self.refiner.resolve(
                    {"source_gse_id": gse, "tissue": tissue}
                )
                self.assertEqual(resolution.value, expected)
                self.assertEqual(resolution.level, "gse_tissue")

    def test_diagnosis_is_not_misused_as_specimen_context(self) -> None:
        self.assertNotIn("diagnosis", self.rules["source_columns"]["cell_source"])
        self.assertNotIn("cancer_type", self.rules["source_columns"]["cell_source"])
        self.assertIn("source_accession", self.rules["source_columns"]["source_gse_id"])

    def test_global_tissue_fallback_then_unresolved(self) -> None:
        self.assert_resolution(
            {"source_gse_id": "GSE_TEST", "tissue": "Spleen"},
            "spleen",
            "global_tissue_fallback",
            "global_tissue_alias",
        )
        self.assert_resolution(
            {"source_gse_id": "GSE_TEST", "tissue": "synovial tissue"},
            "synovial_tissue",
            "global_tissue_fallback",
            "global_tissue_passthrough",
        )
        self.assert_resolution(
            {"source_gse_id": "GSE_TEST", "tissue": "matrix.mtx.gz"},
            "unresolved",
            "unresolved",
            "unresolved_no_supported_source",
        )


class SampleSourceH5adWorkflowTests(unittest.TestCase):
    def build_fixture(self, path: Path) -> tuple[pd.DataFrame, np.ndarray]:
        obs = pd.DataFrame(
            {
                "source_gse_id": pd.Categorical(
                    [
                        "GSE243013",
                        "GSE243013",
                        "GSE287301",
                        "GSE162498",
                        "GSE_TEST",
                        "GSE_TEST",
                    ]
                ),
                "cell_source": ["LUAD", "", "", "", "", ""],
                "sample_id": ["P1", "P1", "chip1pool3", "p34", "S1", "S2"],
                "library_id": ["", "", "GSM8743444", "", "L1", "L2"],
                "donor_patient": ["", "", "3;10;18;27", "", "D1", "D2"],
                "tissue_corrected": pd.Categorical(
                    ["tumor", "tumor", "tumor", "tumor", "spleen", "unknown"]
                ),
                "original_note": ["a", "b", "c", "d", "e", "f"],
            },
            index=[f"cell-{index}" for index in range(6)],
        )
        matrix = np.arange(18, dtype=np.float32).reshape(6, 3)
        adata = ad.AnnData(
            X=sparse.csr_matrix(matrix),
            obs=obs.copy(),
            var=pd.DataFrame({"symbol": ["A", "B", "C"]}, index=["g1", "g2", "g3"]),
        )
        adata.layers["counts"] = sparse.csr_matrix(matrix + 1)
        adata.obsm["X_test"] = np.arange(12, dtype=np.float32).reshape(6, 2)
        adata.uns["fixture"] = {"purpose": "preservation test", "version": 1}
        adata.write_h5ad(path, compression="gzip")
        with h5py.File(path, "r+") as handle:
            handle["obs"].move("_index", "obs_index")
            handle["obs"].attrs["_index"] = "obs_index"
        return obs, matrix

    def test_review_is_deterministic_and_writeback_preserves_h5ad(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            h5ad_path = root / "fixture.h5ad"
            original_obs, original_matrix = self.build_fixture(h5ad_path)
            original_file_hash = file_sha256(h5ad_path)

            first_manifest_path = review_h5ad(
                h5ad_path,
                RULES_PATH,
                root / "review_first",
                chunk_size=2,
                unresolved_limit=10,
            )
            self.assertEqual(file_sha256(h5ad_path), original_file_hash)
            with h5py.File(h5ad_path, "r") as handle:
                self.assertNotIn("sample_source_refined", handle["obs"])

            second_manifest_path = review_h5ad(
                h5ad_path,
                RULES_PATH,
                root / "review_second",
                chunk_size=3,
                unresolved_limit=10,
            )
            self.assertEqual(file_sha256(h5ad_path), original_file_hash)
            first_export = first_manifest_path.parent / "sample_source_refined_review.csv.gz"
            second_export = second_manifest_path.parent / "sample_source_refined_review.csv.gz"
            self.assertEqual(file_sha256(first_export), file_sha256(second_export))

            first_manifest = json.loads(first_manifest_path.read_text(encoding="utf-8"))
            second_manifest = json.loads(second_manifest_path.read_text(encoding="utf-8"))
            for key in [
                "obs_rows",
                "source_values_sha256",
                "output_values_sha256",
                "resolutions_sha256",
                "source_columns",
            ]:
                self.assertEqual(first_manifest[key], second_manifest[key])
            self.assertTrue(first_manifest["tumor_context_gate_pass"])
            self.assertEqual(first_manifest["tumor_context_violations"], 0)
            self.assertTrue(
                (first_manifest_path.parent / "sample_source_refined_tumor_project_audit.csv").is_file()
            )

            with self.assertRaisesRegex(RefinementError, "confirm_reviewed"):
                writeback_h5ad(
                    h5ad_path,
                    RULES_PATH,
                    first_manifest_path,
                    confirm_reviewed=False,
                    chunk_size=2,
                )
            self.assertEqual(file_sha256(h5ad_path), original_file_hash)

            changed_rules = json.loads(RULES_PATH.read_text(encoding="utf-8"))
            changed_rules["global_tissue_aliases"]["spleen"] = "changed_spleen"
            changed_rules_path = root / "changed_rules.json"
            changed_rules_path.write_text(json.dumps(changed_rules), encoding="utf-8")
            with self.assertRaisesRegex(RefinementError, "Rules changed"):
                writeback_h5ad(
                    h5ad_path,
                    changed_rules_path,
                    first_manifest_path,
                    confirm_reviewed=True,
                    chunk_size=2,
                )
            self.assertEqual(file_sha256(h5ad_path), original_file_hash)

            report_path = writeback_h5ad(
                h5ad_path,
                RULES_PATH,
                first_manifest_path,
                confirm_reviewed=True,
                chunk_size=2,
            )
            self.assertTrue(report_path.is_file())

            observed = ad.read_h5ad(h5ad_path)
            self.assertEqual(
                observed.obs["sample_source_refined"].astype(str).tolist(),
                [
                    "LUAD_tumor",
                    "LUSC_tumor",
                    "tumor_pool_HNSCC_with_ameloblastoma_ambiguity",
                    "NSCLC_tumor",
                    "spleen",
                    "unresolved",
                ],
            )
            for column in original_obs.columns:
                self.assertEqual(
                    observed.obs[column].astype(str).tolist(),
                    original_obs[column].astype(str).tolist(),
                    column,
                )
            self.assertEqual(
                set(observed.obs.columns) - set(original_obs.columns),
                {"sample_source_refined"},
            )
            np.testing.assert_array_equal(observed.X.toarray(), original_matrix)
            np.testing.assert_array_equal(observed.layers["counts"].toarray(), original_matrix + 1)
            np.testing.assert_array_equal(
                observed.obsm["X_test"],
                np.arange(12, dtype=np.float32).reshape(6, 2),
            )
            self.assertEqual(observed.var["symbol"].tolist(), ["A", "B", "C"])
            self.assertEqual(observed.uns["fixture"]["purpose"], "preservation test")


if __name__ == "__main__":
    unittest.main()
