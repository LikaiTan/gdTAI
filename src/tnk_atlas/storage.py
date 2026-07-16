"""Transactional same-filesystem storage moves with legacy compatibility links."""

from __future__ import annotations

import csv
import os
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


PLAN_COLUMNS = (
    "sequence",
    "operation",
    "dataset_id",
    "role",
    "old_path",
    "new_path",
    "expected_kind",
    "device",
    "inode",
    "size_bytes",
    "mtime_ns",
)


@dataclass(frozen=True)
class StorageMove:
    sequence: int
    operation: str
    dataset_id: str
    role: str
    old_path: str
    new_path: str
    expected_kind: str
    device: int
    inode: int
    size_bytes: int
    mtime_ns: int

    @classmethod
    def from_path(
        cls,
        *,
        sequence: int,
        operation: str,
        dataset_id: str,
        role: str,
        old_path: Path,
        new_path: Path,
        root: Path,
    ) -> "StorageMove":
        stat = old_path.stat()
        return cls(
            sequence=sequence,
            operation=operation,
            dataset_id=dataset_id,
            role=role,
            old_path=old_path.relative_to(root).as_posix(),
            new_path=new_path.relative_to(root).as_posix(),
            expected_kind="directory" if old_path.is_dir() else "file",
            device=stat.st_dev,
            inode=stat.st_ino,
            size_bytes=stat.st_size,
            mtime_ns=stat.st_mtime_ns,
        )


def write_plan(path: Path, moves: Iterable[StorageMove]) -> None:
    rows = [asdict(move) for move in moves]
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=PLAN_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def read_plan(path: Path) -> list[StorageMove]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    moves = [
        StorageMove(
            sequence=int(row["sequence"]),
            operation=row["operation"],
            dataset_id=row["dataset_id"],
            role=row["role"],
            old_path=row["old_path"],
            new_path=row["new_path"],
            expected_kind=row["expected_kind"],
            device=int(row["device"]),
            inode=int(row["inode"]),
            size_bytes=int(row["size_bytes"]),
            mtime_ns=int(row["mtime_ns"]),
        )
        for row in rows
    ]
    return sorted(moves, key=lambda move: move.sequence)


def relative_link_target(link: Path, target: Path) -> str:
    # Resolve the parent because it may itself be reached through a legacy
    # compatibility symlink whose physical location differs from its spelling.
    return os.path.relpath(target.resolve(), start=link.parent.resolve())


def create_compatibility_link(link: Path, target: Path) -> None:
    link.parent.mkdir(parents=True, exist_ok=True)
    link.symlink_to(
        relative_link_target(link, target),
        target_is_directory=target.is_dir(),
    )


def _same_object(path: Path, move: StorageMove) -> bool:
    if not path.exists():
        return False
    stat = path.stat()
    return stat.st_dev == move.device and stat.st_ino == move.inode


def apply_move(root: Path, move: StorageMove) -> str:
    """Apply or resume one rename-and-link operation."""
    old = root / move.old_path
    new = root / move.new_path

    if old.is_symlink():
        if not new.exists() or old.resolve() != new.resolve():
            raise RuntimeError(f"Compatibility link has the wrong target: {old}")
        if not _same_object(new, move):
            raise RuntimeError(f"Destination inode differs from the plan: {new}")
        return "already_applied"

    if old.exists() and new.exists():
        raise FileExistsError(f"Both source and destination exist: {old}, {new}")

    if old.exists():
        if not _same_object(old, move):
            raise RuntimeError(f"Source inode differs from the plan: {old}")
        new.parent.mkdir(parents=True, exist_ok=True)
        os.rename(old, new)
        try:
            create_compatibility_link(old, new)
        except Exception:
            if not old.exists() and new.exists():
                os.rename(new, old)
            raise
        return "applied"

    if new.exists():
        if not _same_object(new, move):
            raise RuntimeError(f"Destination inode differs from the plan: {new}")
        create_compatibility_link(old, new)
        return "recovered_link"

    raise FileNotFoundError(f"Neither source nor destination exists: {old}, {new}")


def validate_move(root: Path, move: StorageMove) -> list[str]:
    old = root / move.old_path
    new = root / move.new_path
    errors: list[str] = []
    if not old.is_symlink():
        errors.append(f"legacy path is not a symlink: {move.old_path}")
    if not new.exists():
        errors.append(f"canonical destination is missing: {move.new_path}")
        return errors
    if old.is_symlink() and old.resolve() != new.resolve():
        errors.append(f"legacy link target mismatch: {move.old_path}")
    stat = new.stat()
    if stat.st_dev != move.device or stat.st_ino != move.inode:
        errors.append(f"inode mismatch: {move.new_path}")
    if move.expected_kind == "file":
        if not new.is_file():
            errors.append(f"expected file: {move.new_path}")
        if stat.st_size != move.size_bytes or stat.st_mtime_ns != move.mtime_ns:
            errors.append(f"file stat mismatch: {move.new_path}")
    elif not new.is_dir():
        errors.append(f"expected directory: {move.new_path}")
    return errors


def rollback_move(root: Path, move: StorageMove) -> str:
    """Reverse one operation. Call moves in reverse sequence order."""
    old = root / move.old_path
    new = root / move.new_path

    if old.exists() and not old.is_symlink() and not new.exists():
        if not _same_object(old, move):
            raise RuntimeError(f"Restored source inode differs from the plan: {old}")
        return "already_rolled_back"

    if old.is_symlink():
        if not new.exists() or old.resolve() != new.resolve():
            raise RuntimeError(f"Refusing to remove a mismatched compatibility link: {old}")
        old.unlink()
        old.parent.mkdir(parents=True, exist_ok=True)
        os.rename(new, old)
        return "rolled_back"

    if not old.exists() and new.exists():
        old.parent.mkdir(parents=True, exist_ok=True)
        os.rename(new, old)
        return "recovered_source"

    raise RuntimeError(f"Cannot safely roll back {old} <- {new}")


def translate_path(path: Path, root: Path, moves: Iterable[StorageMove]) -> Path:
    """Translate a pre-migration path to its canonical post-migration location."""
    absolute = path if path.is_absolute() else root / path
    absolute = absolute.absolute()
    candidates: list[tuple[int, Path]] = []
    for move in moves:
        old = (root / move.old_path).absolute()
        try:
            suffix = absolute.relative_to(old)
        except ValueError:
            continue
        candidates.append((len(old.parts), (root / move.new_path) / suffix))
    if not candidates:
        return absolute
    return max(candidates, key=lambda item: item[0])[1]
