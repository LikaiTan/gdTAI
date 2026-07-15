"""Loading helpers for trusted local gdTAI pickle artifacts."""

from __future__ import annotations

import pickle
import sys
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator

from .paths import find_project_root


@contextmanager
def legacy_model_import_paths(project_root: Path) -> Iterator[None]:
    """Temporarily expose modules referenced by historical pickle artifacts."""
    candidates = (
        project_root / "src",
        project_root / "workflows" / "gdtai",
        project_root / "workflows" / "integration",
    )
    inserted: list[str] = []
    for candidate in candidates:
        value = str(candidate)
        if value not in sys.path:
            sys.path.insert(0, value)
            inserted.append(value)
    try:
        yield
    finally:
        for value in inserted:
            if value in sys.path:
                sys.path.remove(value)


def load_trusted_pickle(path: str | Path, project_root: str | Path | None = None) -> Any:
    """Load a trusted project pickle while resolving pre-reorganization modules.

    Pickle is executable content. Call this only for artifacts from the project
    model registry whose checksum has been verified.
    """
    artifact = Path(path).expanduser().resolve()
    root = Path(project_root).resolve() if project_root else find_project_root(artifact)
    with legacy_model_import_paths(root), artifact.open("rb") as handle:
        return pickle.load(handle)
