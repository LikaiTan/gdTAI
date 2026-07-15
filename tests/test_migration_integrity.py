from pathlib import Path
import unittest

from tnk_atlas.paths import find_project_root
from tnk_atlas.registry import read_csv


ROOT = find_project_root(Path(__file__))
PLANS = (
    "script_migration_plan.csv",
    "preintegration_script_migration_plan.csv",
    "config_data_migration_plan.csv",
    "config_and_r_migration_plan.csv",
    "legacy_archive_plan.csv",
    "retired_workflow_archive_plan.csv",
)


class MigrationIntegrityTests(unittest.TestCase):
    def test_migration_plans_are_fully_applied(self) -> None:
        rows_by_plan = {
            plan_name: read_csv(ROOT / "maintenance" / "reorganization" / plan_name)
            for plan_name in PLANS
        }
        next_path = {
            row["old_path"]: row["new_path"]
            for rows in rows_by_plan.values()
            for row in rows
        }

        def terminal_path(value: str) -> str:
            seen: set[str] = set()
            while value in next_path:
                self.assertNotIn(value, seen, f"migration cycle at {value}")
                seen.add(value)
                value = next_path[value]
            return value

        for plan_name in PLANS:
            with self.subTest(plan=plan_name):
                rows = rows_by_plan[plan_name]
                self.assertTrue(rows)
                for row in rows:
                    self.assertFalse((ROOT / row["old_path"]).exists(), row["old_path"])
                    terminal = terminal_path(row["new_path"])
                    self.assertTrue((ROOT / terminal).exists(), terminal)
