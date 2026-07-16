from pathlib import Path
import tempfile
import unittest

from tnk_atlas.storage import (
    StorageMove,
    apply_absolute_move_with_link,
    apply_move,
    rollback_absolute_move,
    rollback_move,
    translate_path,
    validate_absolute_move,
    validate_move,
)


class StorageMigrationTests(unittest.TestCase):
    def test_move_link_validate_and_rollback(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            old = root / "downloads" / "GSE1"
            old.mkdir(parents=True)
            payload = old / "matrix.txt"
            payload.write_text("counts\n", encoding="utf-8")
            original_inode = old.stat().st_ino
            new = root / "data" / "datasets" / "GSE1" / "raw" / "legacy_source"
            move = StorageMove.from_path(
                sequence=1,
                operation="move_and_link",
                dataset_id="GSE1",
                role="source",
                old_path=old,
                new_path=new,
                root=root,
            )

            self.assertEqual(apply_move(root, move), "applied")
            self.assertTrue(old.is_symlink())
            self.assertEqual(old.resolve(), new.resolve())
            self.assertEqual(new.stat().st_ino, original_inode)
            self.assertFalse(validate_move(root, move))
            self.assertEqual(apply_move(root, move), "already_applied")

            translated = translate_path(
                Path("downloads/GSE1/matrix.txt"),
                root,
                [move],
            )
            self.assertEqual(translated, new / "matrix.txt")

            self.assertEqual(rollback_move(root, move), "rolled_back")
            self.assertTrue(old.is_dir())
            self.assertFalse(old.is_symlink())
            self.assertFalse(new.exists())
            self.assertEqual((old / "matrix.txt").read_text(encoding="utf-8"), "counts\n")

    def test_interrupted_move_recovers_compatibility_link(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            old = root / "newdata" / "dataset.h5ad"
            old.parent.mkdir(parents=True)
            old.write_bytes(b"h5ad")
            new = root / "data" / "datasets" / "GSE1" / "processed" / "dataset.h5ad"
            move = StorageMove.from_path(
                sequence=1,
                operation="move_and_link",
                dataset_id="GSE1",
                role="processed",
                old_path=old,
                new_path=new,
                root=root,
            )
            new.parent.mkdir(parents=True)
            old.rename(new)

            self.assertEqual(apply_move(root, move), "recovered_link")
            self.assertTrue(old.is_symlink())
            self.assertFalse(validate_move(root, move))

    def test_root_alias_then_child_move_keeps_relative_link_valid(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            legacy_root = root / "downloads"
            child = legacy_root / "GSE1"
            child.mkdir(parents=True)
            (child / "matrix.txt").write_text("counts\n", encoding="utf-8")
            compatibility_root = root / "data" / "compat" / "downloads"
            canonical_child = (
                root / "data" / "datasets" / "GSE1" / "raw" / "legacy_source"
            )
            root_move = StorageMove.from_path(
                sequence=1,
                operation="root_alias",
                dataset_id="",
                role="root",
                old_path=legacy_root,
                new_path=compatibility_root,
                root=root,
            )
            child_move = StorageMove.from_path(
                sequence=2,
                operation="move_and_link",
                dataset_id="GSE1",
                role="source",
                old_path=child,
                new_path=canonical_child,
                root=root,
            )

            self.assertEqual(apply_move(root, root_move), "applied")
            self.assertEqual(apply_move(root, child_move), "applied")
            self.assertEqual(child.resolve(), canonical_child.resolve())
            self.assertEqual(
                (child / "matrix.txt").read_text(encoding="utf-8"),
                "counts\n",
            )
            self.assertFalse(validate_move(root, root_move))
            self.assertFalse(validate_move(root, child_move))

            self.assertEqual(rollback_move(root, child_move), "rolled_back")
            self.assertEqual(rollback_move(root, root_move), "rolled_back")
            self.assertTrue(child.is_dir())
            self.assertFalse(legacy_root.is_symlink())

    def test_absolute_workspace_move_and_rollback(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "outside" / "study"
            source.mkdir(parents=True)
            (source / "result.h5ad").write_bytes(b"h5ad")
            destination = root / "project" / "data" / "datasets" / "STUDY" / "workspace"
            stat = source.stat()

            self.assertEqual(
                apply_absolute_move_with_link(
                    source,
                    destination,
                    expected_device=stat.st_dev,
                    expected_inode=stat.st_ino,
                ),
                "applied",
            )
            self.assertEqual(source.resolve(), destination.resolve())
            self.assertFalse(
                validate_absolute_move(
                    source,
                    destination,
                    expected_device=stat.st_dev,
                    expected_inode=stat.st_ino,
                )
            )
            self.assertEqual(
                rollback_absolute_move(
                    source,
                    destination,
                    expected_device=stat.st_dev,
                    expected_inode=stat.st_ino,
                ),
                "rolled_back",
            )
            self.assertTrue(source.is_dir())
            self.assertFalse(source.is_symlink())
            self.assertFalse(destination.exists())
