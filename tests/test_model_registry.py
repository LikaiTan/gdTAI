from pathlib import Path
import unittest

from tnk_atlas.model_io import load_trusted_pickle
from tnk_atlas.paths import find_project_root
from tnk_atlas.provenance import sha256_file
from tnk_atlas.registry import read_csv


ROOT = find_project_root(Path(__file__))


class ModelRegistryTests(unittest.TestCase):
    def test_registered_model_checksums_and_loading(self) -> None:
        registry = read_csv(ROOT / "configs" / "models" / "gdtai" / "model_registry.csv")
        self.assertGreaterEqual(
            {row["status"] for row in registry},
            {"promoted_default", "reference_release", "experimental_not_promoted"},
        )
        for row in registry:
            artifact = ROOT / row["artifact_path"]
            required = row["required_in_git_checkout"].lower() == "true"
            if not artifact.exists():
                self.assertFalse(required, row["model_id"])
                continue
            observed = sha256_file(artifact)
            accepted = {row["sha256"]}
            if row["allowed_workspace_sha256"]:
                accepted.add(row["allowed_workspace_sha256"])
            self.assertIn(observed, accepted)
            self.assertIsInstance(load_trusted_pickle(artifact, ROOT), dict)
