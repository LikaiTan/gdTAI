import json
from pathlib import Path
import unittest

from tnk_atlas.paths import find_project_root


ROOT = find_project_root(Path(__file__))


class CompatibilityLinkTests(unittest.TestCase):
    def test_compatibility_links_resolve(self) -> None:
        manifest_path = ROOT / "data" / "registry" / "compatibility_links.json"
        rows = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertGreaterEqual(len(rows), 250)
        for row in rows:
            link = ROOT / row["link_path"]
            target = Path(row["target_path"])
            if not target.is_absolute():
                target = ROOT / target
            self.assertTrue(link.is_symlink(), link)
            self.assertTrue(target.exists(), target)
            self.assertEqual(link.resolve(), target.resolve(), link)
