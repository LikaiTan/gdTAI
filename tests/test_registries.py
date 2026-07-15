from pathlib import Path
import unittest

from tnk_atlas.paths import find_project_root
from tnk_atlas.registry import read_csv, validate_dataset_registry


ROOT = find_project_root(Path(__file__))
REGISTRY = ROOT / "configs" / "datasets"


class RegistryTests(unittest.TestCase):
    def test_dataset_registry_is_valid(self) -> None:
        result = validate_dataset_registry(REGISTRY / "datasets.csv", ROOT)
        self.assertTrue(result.ok)
        self.assertEqual(result.dataset_count, 65)
        self.assertEqual(result.active_dataset_count, 33)

    def test_every_dataset_has_curated_library_rows(self) -> None:
        dataset_ids = {row["dataset_id"] for row in read_csv(REGISTRY / "datasets.csv")}
        libraries = read_csv(REGISTRY / "libraries.csv")
        self.assertLessEqual(dataset_ids, {row["dataset_id"] for row in libraries})
        placeholders = [
            row for row in libraries
            if row["active"].lower() == "true"
            and row["library_id"] == "__dataset_level_pending_curation__"
        ]
        self.assertFalse(placeholders)
