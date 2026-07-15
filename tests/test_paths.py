from pathlib import Path
import unittest

from tnk_atlas.paths import ProjectPaths, find_project_root


class ProjectPathTests(unittest.TestCase):
    def test_project_root_discovery_from_nested_path(self) -> None:
        root = find_project_root(Path(__file__))
        self.assertTrue((root / "AGENTS.md").exists())
        self.assertEqual(ProjectPaths.discover(Path(__file__)).root, root)

    def test_canonical_path_vocabulary(self) -> None:
        paths = ProjectPaths.discover(Path(__file__))
        self.assertEqual(paths.dataset_registry.name, "datasets")
        self.assertEqual(paths.outputs.name, "Integrated_dataset")
        self.assertTrue(paths.integrated_h5ad.exists())
