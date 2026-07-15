import csv
import json
from pathlib import Path
import unittest

from tnk_atlas.paths import find_project_root
from tnk_atlas.provenance import sha256_file


ROOT = find_project_root(Path(__file__))
TABLE_DIR = ROOT / "Integrated_dataset" / "tables" / "gdT_prediction" / "gdtai_v3_round12_vs_round14"
LOG_DIR = ROOT / "Integrated_dataset" / "logs" / "gdT_prediction" / "gdtai_v3_round12_vs_round14"
MODEL_DIR = ROOT / "Integrated_dataset" / "models" / "gdT_prediction_classifier"

R12_SHA256 = "7373e79350f7db190c415b376b9763e31652754438ee8c5afd3853beb7b2ebc4"
R14_SHA256 = "16dedc0081da9b8487887341232bcf8c9c9403dd3bbd72e04cab43d4cd7b2e09"


class GdtaiRoundSelectionTests(unittest.TestCase):
    def test_round14_is_canonical_and_round12_is_preserved(self) -> None:
        canonical = MODEL_DIR / "gdTAI_v3.0" / "gdTAI_v3_model.pkl"
        fallback = MODEL_DIR / "gdtai_v3_round12_vs_round14" / "round12_model.pkl"
        self.assertEqual(sha256_file(canonical), R14_SHA256)
        self.assertEqual(sha256_file(fallback), R12_SHA256)

        manifest = json.loads((MODEL_DIR / "gdTAI_v3.0" / "model_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["model_name"], "v3_round14_v2_score_trdc_gate_fixed_0p936")
        self.assertEqual(manifest["model_sha256"], R14_SHA256)
        self.assertEqual(float(manifest["threshold"]), 0.936)

    def test_comparison_decision_and_cache_checks(self) -> None:
        decision = json.loads((TABLE_DIR / "promotion_decision.json").read_text(encoding="utf-8"))
        self.assertEqual(decision["selected_model"], "round14")

        with (TABLE_DIR / "promotion_guardrails.csv").open(newline="", encoding="utf-8") as handle:
            guardrails = {row["model"]: row for row in csv.DictReader(handle)}
        self.assertEqual(guardrails["round12"]["all_guardrails_pass"], "False")
        self.assertEqual(guardrails["round14"]["all_guardrails_pass"], "True")

        with (TABLE_DIR / "artifact_cache_verification.csv").open(newline="", encoding="utf-8") as handle:
            checks = list(csv.DictReader(handle))
        self.assertEqual(len(checks), 4)
        self.assertTrue(all(row["prediction_match"] == "True" for row in checks))

        summary = json.loads(
            (LOG_DIR / "gdtai_v3_round12_vs_round14_summary.json").read_text(encoding="utf-8")
        )
        self.assertTrue(summary["h5ad_unchanged"])
        self.assertTrue(summary["promotion_applied"])

    def test_reader_report_records_the_decision(self) -> None:
        report = ROOT / "gdT_prediction" / "gdtai_v3_round12_vs_round14" / "index.html"
        text = report.read_text(encoding="utf-8")
        self.assertIn("Selected: Round 14", text)
        self.assertIn("Round 14 should remain the canonical balanced model", text)


if __name__ == "__main__":
    unittest.main()
