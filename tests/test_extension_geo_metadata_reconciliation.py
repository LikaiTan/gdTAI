import importlib.util
from pathlib import Path
import sys
import unittest
from urllib.parse import urlparse

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflows/metadata/audit_extension_geo_metadata.py"
SPEC = importlib.util.spec_from_file_location(
    "audit_extension_geo_metadata", SCRIPT
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class ExtensionGeoMetadataReconciliationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config = pd.read_csv(
            MODULE.CONFIG_PATH, dtype=str, keep_default_na=False
        )

    def test_config_schema_and_coverage_are_valid(self) -> None:
        self.assertEqual(MODULE.validate_config_schema(self.config), [])
        self.assertEqual(set(self.config["cohort_id"]), set(MODULE.COHORT_ACCESSIONS))
        self.assertEqual(
            set(self.config["accession"]),
            {
                accession
                for accessions in MODULE.COHORT_ACCESSIONS.values()
                for accession in accessions
            },
        )

    def test_evidence_urls_are_official_ncbi_geo_resources(self) -> None:
        for url in self.config["direct_url"]:
            with self.subTest(url=url):
                parsed = urlparse(url)
                self.assertEqual(parsed.scheme, "https")
                self.assertIn(parsed.hostname, MODULE.ALLOWED_GEO_HOSTS)

    def test_unresolved_fields_remain_explicitly_fail_closed(self) -> None:
        status_counts = self.config["status"].value_counts().to_dict()
        self.assertEqual(status_counts["ambiguous_geo_partial"], 6)
        self.assertEqual(status_counts["unavailable_in_geo"], 4)
        unavailable = self.config["status"].eq("unavailable_in_geo")
        self.assertTrue(self.config.loc[unavailable, "geo_supported_value"].eq("").all())
        self.assertTrue(
            self.config.loc[~unavailable, "geo_supported_value"].ne("").all()
        )

    def test_gse169246_blood_scope_is_additive_and_deterministic(self) -> None:
        rows = self.config.loc[
            self.config["cohort_id"].eq("GSE169246")
            & self.config["field"].isin(["tissue", "specimen_context"])
            & self.config["local_scope_regex"].eq(".*_b$")
        ]
        self.assertEqual(set(rows["field"]), {"tissue", "specimen_context"})
        self.assertTrue(rows["geo_supported_value"].eq("blood").all())
        self.assertTrue(
            rows["status"].eq("resolved_geo_value_local_unresolved").all()
        )
        self.assertTrue(rows["local_column"].ne("").all())


if __name__ == "__main__":
    unittest.main()
