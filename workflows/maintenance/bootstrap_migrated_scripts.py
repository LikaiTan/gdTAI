#!/usr/bin/env python3
"""Make migrated standalone workflows path-stable from any working directory."""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import re
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
MARKER = "# TNK_WORKFLOW_BOOTSTRAP"
BOOTSTRAP = '''# TNK_WORKFLOW_BOOTSTRAP
import sys as _tnk_sys
from pathlib import Path as _TnkPath

_TNK_PROJECT_ROOT = _TnkPath(__file__).resolve().parents[2]
for _tnk_path in (
    _TNK_PROJECT_ROOT,
    _TNK_PROJECT_ROOT / "src",
    _TNK_PROJECT_ROOT / "workflows" / "integration",
    _TNK_PROJECT_ROOT / "workflows" / "intake",
    _TNK_PROJECT_ROOT / "workflows" / "metadata",
    _TNK_PROJECT_ROOT / "workflows" / "analysis",
    _TNK_PROJECT_ROOT / "workflows" / "gdtai",
    _TNK_PROJECT_ROOT / "workflows" / "gdt_atlas",
    _TNK_PROJECT_ROOT / "workflows" / "reporting",
    _TNK_PROJECT_ROOT / "workflows" / "maintenance",
):
    _tnk_value = str(_tnk_path)
    if _tnk_value not in _tnk_sys.path:
        _tnk_sys.path.insert(0, _tnk_value)

'''


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--plan",
        default="maintenance/reorganization/script_migration_plan.csv",
    )
    parser.add_argument(
        "--migration-map",
        default="maintenance/reorganization/migration_map.csv",
        help="Applied migration map containing original SHA256 values.",
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--apply", action="store_true")
    mode.add_argument(
        "--rollback",
        action="store_true",
        help="Restore pre-migration source text and verify its recorded SHA256.",
    )
    return parser.parse_args()


def insertion_line(source: str) -> int:
    tree = ast.parse(source)
    line = 0
    body = list(tree.body)
    if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant):
        if isinstance(body[0].value.value, str):
            line = body[0].end_lineno or 0
            body = body[1:]
    for node in body:
        if isinstance(node, ast.ImportFrom) and node.module == "__future__":
            line = node.end_lineno or line
        else:
            break
    return line


def update_source(path: Path) -> bool:
    source = path.read_text(encoding="utf-8")
    if MARKER in source:
        updated = source
    else:
        lines = source.splitlines(keepends=True)
        index = insertion_line(source)
        lines[index:index] = ["\n", BOOTSTRAP]
        updated = "".join(lines)
    updated = updated.replace(".parents[2]s[2]", ".parents[2]")
    updated = re.sub(
        r"(?<![A-Za-z0-9_])Path\(__file__\)\.resolve\(\)\.parent(?!s)",
        "Path(__file__).resolve().parents[2]",
        updated,
    )
    updated = re.sub(
        r'''Path\((["'])((?:Integrated_dataset|gdT_prediction|gdT_atlas|analysis_26GSE_V4|downloads|newdata|high_speed_temp)(?:/[^"']*)?)\1\)''',
        r'_TNK_PROJECT_ROOT / "\2"',
        updated,
    )
    if updated == source:
        return False
    path.write_text(updated, encoding="utf-8")
    return True


def restore_source(path: Path, expected_sha256: str) -> bool:
    source = path.read_text(encoding="utf-8")
    updated = source.replace("\n" + BOOTSTRAP, "")
    updated = re.sub(
        r'_TNK_PROJECT_ROOT / "((?:Integrated_dataset|gdT_prediction|gdT_atlas|analysis_26GSE_V4|downloads|newdata|high_speed_temp)(?:/[^"]*)?)"',
        r'Path("\1")',
        updated,
    )
    updated = re.sub(
        r"(?<![A-Za-z0-9_])Path\(__file__\)\.resolve\(\)\.parents\[2\]",
        "Path(__file__).resolve().parent",
        updated,
    )
    digest = hashlib.sha256(updated.encode("utf-8")).hexdigest()
    if digest != expected_sha256:
        raise RuntimeError(
            f"Refusing rollback for {path}: reconstructed SHA256 {digest} != {expected_sha256}"
        )
    if updated == source:
        return False
    path.write_text(updated, encoding="utf-8")
    return True


def main() -> int:
    args = parse_args()
    with (PROJECT_ROOT / args.plan).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    sha_by_old: dict[str, str] = {}
    migration_map = PROJECT_ROOT / args.migration_map
    if migration_map.exists():
        with migration_map.open(newline="", encoding="utf-8") as handle:
            sha_by_old = {row["old_path"]: row["sha256"] for row in csv.DictReader(handle)}
    changed = 0
    for row in rows:
        path = PROJECT_ROOT / row["new_path"]
        if not path.exists():
            continue
        if args.rollback:
            expected = sha_by_old.get(row["old_path"], "")
            if not expected:
                raise RuntimeError(f"Missing original SHA256 for {row['old_path']}")
            changed += int(restore_source(path, expected))
        elif args.apply:
            changed += int(update_source(path))
        else:
            source = path.read_text(encoding="utf-8")
            changed += int(MARKER not in source or "Path(__file__).resolve().parent" in source)
    mode = "rollback" if args.rollback else "apply" if args.apply else "dry-run"
    print(f"scripts_requiring_or_receiving_bootstrap={changed} mode={mode}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
