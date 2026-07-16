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


def apply_absolute_move_with_link(
    source: Path,
    destination: Path,
    *,
    expected_device: int,
    expected_inode: int,
) -> str:
    """Move an absolute directory and leave its original path as a link."""
    if source.is_symlink():
        if not destination.exists() or source.resolve() != destination.resolve():
            raise RuntimeError(f"Compatibility link has the wrong target: {source}")
        stat = destination.stat()
        if stat.st_dev != expected_device or stat.st_ino != expected_inode:
            raise RuntimeError(f"Destination inode differs from the plan: {destination}")
        return "already_applied"

    if source.exists() and destination.exists():
        raise FileExistsError(f"Both source and destination exist: {source}, {destination}")

    if source.exists():
        stat = source.stat()
        if stat.st_dev != expected_device or stat.st_ino != expected_inode:
            raise RuntimeError(f"Source inode differs from the plan: {source}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        os.rename(source, destination)
        try:
            create_compatibility_link(source, destination)
        except Exception:
            if not source.exists() and destination.exists():
                os.rename(destination, source)
            raise
        return "applied"

    if destination.exists():
        stat = destination.stat()
        if stat.st_dev != expected_device or stat.st_ino != expected_inode:
            raise RuntimeError(f"Destination inode differs from the plan: {destination}")
        create_compatibility_link(source, destination)
        return "recovered_link"

    raise FileNotFoundError(
        f"Neither source nor destination workspace exists: {source}, {destination}"
    )


def validate_absolute_move(
    source: Path,
    destination: Path,
    *,
    expected_device: int,
    expected_inode: int,
) -> list[str]:
    errors: list[str] = []
    if not source.is_symlink():
        errors.append(f"legacy workspace path is not a symlink: {source}")
    if not destination.is_dir():
        errors.append(f"canonical workspace is missing: {destination}")
        return errors
    if source.is_symlink() and source.resolve() != destination.resolve():
        errors.append(f"legacy workspace target mismatch: {source}")
    stat = destination.stat()
    if stat.st_dev != expected_device or stat.st_ino != expected_inode:
        errors.append(f"workspace inode mismatch: {destination}")
    return errors


def rollback_absolute_move(
    source: Path,
    destination: Path,
    *,
    expected_device: int,
    expected_inode: int,
) -> str:
    """Restore an absolute workspace. Generated children must be removed first."""
    if source.exists() and not source.is_symlink() and not destination.exists():
        stat = source.stat()
        if stat.st_dev != expected_device or stat.st_ino != expected_inode:
            raise RuntimeError(f"Restored workspace inode differs from the plan: {source}")
        return "already_rolled_back"

    if source.is_symlink():
        if not destination.exists() or source.resolve() != destination.resolve():
            raise RuntimeError(f"Refusing to remove a mismatched link: {source}")
        source.unlink()
        source.parent.mkdir(parents=True, exist_ok=True)
        os.rename(destination, source)
        return "rolled_back"

    if not source.exists() and destination.exists():
        source.parent.mkdir(parents=True, exist_ok=True)
        os.rename(destination, source)
        return "recovered_source"

    raise RuntimeError(f"Cannot safely roll back {source} <- {destination}")


def apply_h5ad_only_external_layout(
    external_workspace: Path,
    project_workspace: Path,
    selected_relative: Path,
    canonical_file: Path,
    *,
    expected_workspace_device: int,
    expected_workspace_inode: int,
    expected_file_device: int,
    expected_file_inode: int,
) -> str:
    """Restore a workspace externally while retaining one canonical project file."""
    external_selected = external_workspace / selected_relative
    project_selected = project_workspace / selected_relative

    if (
        external_workspace.is_dir()
        and not external_workspace.is_symlink()
        and not project_workspace.exists()
        and not project_workspace.is_symlink()
        and canonical_file.is_file()
        and not canonical_file.is_symlink()
        and external_selected.is_symlink()
        and external_selected.resolve() == canonical_file.resolve()
    ):
        workspace_stat = external_workspace.stat()
        file_stat = canonical_file.stat()
        if (
            workspace_stat.st_dev != expected_workspace_device
            or workspace_stat.st_ino != expected_workspace_inode
        ):
            raise RuntimeError(
                f"External workspace inode differs from the plan: {external_workspace}"
            )
        if (
            file_stat.st_dev != expected_file_device
            or file_stat.st_ino != expected_file_inode
        ):
            raise RuntimeError(f"Canonical file inode differs from the plan: {canonical_file}")
        return "already_applied"

    if canonical_file.is_symlink():
        if not project_selected.exists() or canonical_file.resolve() != project_selected.resolve():
            raise RuntimeError(f"Canonical-file link has the wrong target: {canonical_file}")
        canonical_file.unlink()
    elif canonical_file.exists():
        stat = canonical_file.stat()
        if stat.st_dev != expected_file_device or stat.st_ino != expected_file_inode:
            raise RuntimeError(f"Canonical file inode differs from the plan: {canonical_file}")

    if project_selected.exists() and not project_selected.is_symlink():
        stat = project_selected.stat()
        if stat.st_dev != expected_file_device or stat.st_ino != expected_file_inode:
            raise RuntimeError(f"Selected file inode differs from the plan: {project_selected}")
        canonical_file.parent.mkdir(parents=True, exist_ok=True)
        if canonical_file.exists():
            raise FileExistsError(f"Canonical file already exists: {canonical_file}")
        os.rename(project_selected, canonical_file)
    elif not canonical_file.exists():
        raise FileNotFoundError(
            f"Neither selected nor canonical file exists: {project_selected}, {canonical_file}"
        )

    if external_workspace.is_symlink():
        if (
            not project_workspace.exists()
            or external_workspace.resolve() != project_workspace.resolve()
        ):
            raise RuntimeError(f"External workspace link has the wrong target: {external_workspace}")
        external_workspace.unlink()
    elif external_workspace.exists() and project_workspace.exists():
        raise FileExistsError(
            f"Both external and project workspaces exist: "
            f"{external_workspace}, {project_workspace}"
        )

    if project_workspace.exists():
        stat = project_workspace.stat()
        if (
            stat.st_dev != expected_workspace_device
            or stat.st_ino != expected_workspace_inode
        ):
            raise RuntimeError(f"Project workspace inode differs from the plan: {project_workspace}")
        external_workspace.parent.mkdir(parents=True, exist_ok=True)
        os.rename(project_workspace, external_workspace)
    elif not external_workspace.exists():
        raise FileNotFoundError(
            f"Neither external nor project workspace exists: "
            f"{external_workspace}, {project_workspace}"
        )

    if external_selected.is_symlink():
        if external_selected.resolve() != canonical_file.resolve():
            raise RuntimeError(f"Selected-file link has the wrong target: {external_selected}")
    elif external_selected.exists():
        raise FileExistsError(f"Selected-file path unexpectedly exists: {external_selected}")
    else:
        create_compatibility_link(external_selected, canonical_file)
    return "applied"


def validate_h5ad_only_external_layout(
    external_workspace: Path,
    project_workspace: Path,
    selected_relative: Path,
    canonical_file: Path,
    *,
    expected_workspace_device: int,
    expected_workspace_inode: int,
    expected_file_device: int,
    expected_file_inode: int,
) -> list[str]:
    errors: list[str] = []
    external_selected = external_workspace / selected_relative
    if external_workspace.is_symlink() or not external_workspace.is_dir():
        errors.append(f"external workspace is not a physical directory: {external_workspace}")
    else:
        stat = external_workspace.stat()
        if (
            stat.st_dev != expected_workspace_device
            or stat.st_ino != expected_workspace_inode
        ):
            errors.append(f"external workspace inode mismatch: {external_workspace}")
    if project_workspace.exists() or project_workspace.is_symlink():
        errors.append(f"project workspace still exists: {project_workspace}")
    if canonical_file.is_symlink() or not canonical_file.is_file():
        errors.append(f"canonical file is not physical: {canonical_file}")
    else:
        stat = canonical_file.stat()
        if stat.st_dev != expected_file_device or stat.st_ino != expected_file_inode:
            errors.append(f"canonical file inode mismatch: {canonical_file}")
    if not external_selected.is_symlink():
        errors.append(f"legacy selected-file path is not a symlink: {external_selected}")
    elif canonical_file.exists() and external_selected.resolve() != canonical_file.resolve():
        errors.append(f"legacy selected-file target mismatch: {external_selected}")
    return errors


def rollback_h5ad_only_external_layout(
    external_workspace: Path,
    project_workspace: Path,
    selected_relative: Path,
    canonical_file: Path,
    *,
    expected_workspace_device: int,
    expected_workspace_inode: int,
    expected_file_device: int,
    expected_file_inode: int,
) -> str:
    """Restore the prior project-hosted workspace layout."""
    external_selected = external_workspace / selected_relative
    project_selected = project_workspace / selected_relative

    if (
        external_workspace.is_symlink()
        and project_workspace.is_dir()
        and canonical_file.is_symlink()
        and project_selected.is_file()
    ):
        return "already_rolled_back"

    errors = validate_h5ad_only_external_layout(
        external_workspace,
        project_workspace,
        selected_relative,
        canonical_file,
        expected_workspace_device=expected_workspace_device,
        expected_workspace_inode=expected_workspace_inode,
        expected_file_device=expected_file_device,
        expected_file_inode=expected_file_inode,
    )
    if errors:
        raise RuntimeError("Cannot safely roll back:\n" + "\n".join(errors))

    external_selected.unlink()
    project_workspace.parent.mkdir(parents=True, exist_ok=True)
    os.rename(external_workspace, project_workspace)
    project_selected.parent.mkdir(parents=True, exist_ok=True)
    os.rename(canonical_file, project_selected)
    create_compatibility_link(external_workspace, project_workspace)
    create_compatibility_link(canonical_file, project_selected)
    return "rolled_back"


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
