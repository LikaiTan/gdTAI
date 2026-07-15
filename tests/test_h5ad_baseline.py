import json
from pathlib import Path
import unittest

from tnk_atlas.paths import find_project_root


ROOT = find_project_root(Path(__file__))
BASELINE = (
    ROOT
    / "Integrated_dataset"
    / "logs"
    / "project_reorganization"
    / "checkpoints"
    / "pre_reorganization_20260715"
    / "h5ad_structures.json"
)


class H5adBaselineTests(unittest.TestCase):
    def test_milestone_h5ad_stats_are_unchanged(self) -> None:
        for row in json.loads(BASELINE.read_text(encoding="utf-8")):
            path = Path(row["path"])
            self.assertTrue(path.exists(), path)
            self.assertEqual(path.stat().st_size, row["size_bytes"], path)
            self.assertEqual(path.stat().st_mtime_ns, row["mtime_ns"], path)
