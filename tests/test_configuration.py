import csv
import configparser
import json
from pathlib import Path
import unittest

from tnk_atlas.paths import find_project_root
from tnk_atlas.registry import DATASET_COLUMNS, FILE_COLUMNS, LIBRARY_COLUMNS


ROOT = find_project_root(Path(__file__))


class ConfigurationTests(unittest.TestCase):
    def test_json_and_toml_configuration_parse(self) -> None:
        for path in (ROOT / "configs").rglob("*.json"):
            with self.subTest(path=path):
                json.loads(path.read_text(encoding="utf-8"))
        parsed = configparser.ConfigParser()
        parsed.read(ROOT / "configs" / "project_paths.toml", encoding="utf-8")
        environment = parsed.get("project", "canonical_environment").strip('"')
        self.assertEqual(environment, "rapids_sc_py310")

    def test_registry_headers_are_stable(self) -> None:
        expected = {
            "datasets.csv": DATASET_COLUMNS,
            "libraries.csv": LIBRARY_COLUMNS,
            "files.csv": FILE_COLUMNS,
        }
        for name, columns in expected.items():
            with self.subTest(registry=name):
                with (ROOT / "configs" / "datasets" / name).open(
                    newline="", encoding="utf-8"
                ) as handle:
                    header = next(csv.reader(handle))
                self.assertEqual(tuple(header), columns)

    def test_root_has_no_python_entrypoints(self) -> None:
        self.assertEqual(list(ROOT.glob("*.py")), [])

    def test_active_workflows_do_not_import_archive(self) -> None:
        offenders = []
        for path in (ROOT / "workflows").rglob("*.py"):
            source = path.read_text(encoding="utf-8")
            if "from archive" in source or "import archive" in source:
                offenders.append(path)
        self.assertFalse(offenders)
