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
        self.assertEqual(result.dataset_count, 66)
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

    def test_balf_blood_copd_is_validation_only(self) -> None:
        datasets = {
            row["dataset_id"]: row
            for row in read_csv(REGISTRY / "datasets.csv")
        }
        row = datasets["BALF_BLOOD_COPD"]
        self.assertEqual(row["integration_role"], "gdTAI_independent_external_validation")
        self.assertEqual(row["phase0_active"], "false")
        self.assertEqual(row["current_milestone_active"], "false")
        self.assertEqual(row["extended_atlas_active"], "false")
        current = ROOT / row["processed_h5ad_path"]
        artifact = (
            ROOT
            / "data"
            / "datasets"
            / "BALF_BLOOD_COPD"
            / "processed"
            / "artifacts"
            / "phase4_final_annotated.h5ad"
        )
        raw_source = (
            ROOT
            / "data"
            / "datasets"
            / "BALF_BLOOD_COPD"
            / "raw"
            / "legacy_source"
        )
        external_workspace = Path("/home/tanlikai/databank/owndata/singlecell")
        legacy_h5ad = (
            external_workspace
            / "data"
            / "results"
            / "phase4_final_annotated.h5ad"
        )

        self.assertTrue(current.is_file())
        self.assertTrue(current.is_symlink())
        self.assertTrue(artifact.is_file())
        self.assertFalse(artifact.is_symlink())
        self.assertEqual(current.resolve(), artifact.resolve())
        project_workspace = (
            ROOT / "data" / "datasets" / "BALF_BLOOD_COPD" / "workspace"
        )
        self.assertFalse(project_workspace.exists())

        self.assertTrue(external_workspace.is_dir())
        self.assertFalse(external_workspace.is_symlink())
        self.assertTrue(raw_source.is_symlink())
        self.assertEqual(raw_source.resolve(), external_workspace.resolve())

        self.assertTrue(legacy_h5ad.is_symlink())
        self.assertEqual(legacy_h5ad.resolve(), artifact.resolve())
