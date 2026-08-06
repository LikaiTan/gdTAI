import csv
import json
from pathlib import Path
import unittest

from tnk_atlas.model_io import load_trusted_pickle
from tnk_atlas.paths import find_project_root


ROOT = find_project_root(Path(__file__))
MODEL_DIR = ROOT / "Integrated_dataset/models/gdT_prediction_classifier/gdTAI_v3.0"
AUDIT_TABLE_DIR = ROOT / "Integrated_dataset/tables/gdT_prediction/gdtai_methodology_audit"


class GdtaiMethodologyAuditTests(unittest.TestCase):
    def test_v3_release_documents_semantic_mismatch(self) -> None:
        manifest = json.loads((MODEL_DIR / "model_manifest.json").read_text(encoding="utf-8"))
        payload = load_trusted_pickle(MODEL_DIR / "gdTAI_v3_model.pkl", ROOT)

        self.assertEqual(manifest["model_name"], payload["model"])
        self.assertTrue(manifest["accepted_for_promotion"])
        self.assertFalse(payload["accepted_for_promotion"])
        self.assertFalse(manifest["payload_embedded_accepted_for_promotion"])
        self.assertEqual(
            manifest["external_evaluation_status"],
            "reused_cross_study_benchmark_not_independent",
        )

    def test_registry_marks_v3_semantic_mismatch(self) -> None:
        registry_path = ROOT / "configs/models/gdtai/model_registry.csv"
        with registry_path.open(newline="", encoding="utf-8") as handle:
            registry = {row["model_id"]: row for row in csv.DictReader(handle)}
        self.assertEqual(registry["gdtai_v3"]["metadata_consistency"], "semantic_mismatch_documented")
        self.assertIn("reused cross-study benchmark", registry["gdtai_v3"]["notes"])

    def test_corrected_tcr_only_benchmark_is_pinned(self) -> None:
        table_path = AUDIT_TABLE_DIR / "corrected_external_pure_tcr_metrics.csv"
        with table_path.open(newline="", encoding="utf-8") as handle:
            rows = {row["model"]: row for row in csv.DictReader(handle)}

        self.assertEqual(int(rows["V3 Round 12"]["n_positive"]), 1046)
        self.assertEqual(int(rows["V3 Round 12"]["n_negative"]), 33117)
        self.assertAlmostEqual(float(rows["V3 Round 12"]["f1"]), 0.8929313929, places=8)
        self.assertAlmostEqual(float(rows["V3 Round 14"]["f1"]), 0.8921971253, places=8)
        self.assertLess(int(rows["V3 Round 12"]["fp"]), int(rows["V3 Round 14"]["fp"]))


if __name__ == "__main__":
    unittest.main()
