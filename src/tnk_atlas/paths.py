"""Canonical project-path resolution independent of the current directory."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


ROOT_MARKERS = ("AGENTS.md", "TNK_PIPELINE_RUNBOOK.md", "TNK_PHASES_0_4_SCRIPT.md")


def find_project_root(start: str | Path | None = None) -> Path:
    """Return the repository root using an override or marker-file discovery."""
    override = os.environ.get("TNK_PROJECT_ROOT")
    if override:
        root = Path(override).expanduser().resolve()
        _validate_root(root)
        return root

    origin = Path(start).expanduser().resolve() if start else Path(__file__).resolve()
    if origin.is_file():
        origin = origin.parent
    for candidate in (origin, *origin.parents):
        if all((candidate / marker).exists() for marker in ROOT_MARKERS):
            return candidate
    raise FileNotFoundError(
        f"Could not find the project root from {origin}; set TNK_PROJECT_ROOT explicitly."
    )


def _validate_root(root: Path) -> None:
    missing = [marker for marker in ROOT_MARKERS if not (root / marker).exists()]
    if missing:
        raise FileNotFoundError(f"Invalid TNK_PROJECT_ROOT {root}: missing {missing}")


@dataclass(frozen=True)
class ProjectPaths:
    """Stable path vocabulary shared by workflows and maintenance tooling."""

    root: Path

    @classmethod
    def discover(cls, start: str | Path | None = None) -> "ProjectPaths":
        return cls(find_project_root(start))

    @property
    def outputs(self) -> Path:
        return self.root / "Integrated_dataset"

    @property
    def tables(self) -> Path:
        return self.outputs / "tables"

    @property
    def figures(self) -> Path:
        return self.outputs / "figures"

    @property
    def logs(self) -> Path:
        return self.outputs / "logs"

    @property
    def models(self) -> Path:
        return self.outputs / "models"

    @property
    def data(self) -> Path:
        return self.root / "data"

    @property
    def dataset_registry(self) -> Path:
        return self.root / "configs" / "datasets"

    @property
    def legacy_downloads(self) -> Path:
        return self.root / "downloads"

    @property
    def legacy_scanpy_projects(self) -> Path:
        return self.root / "analysis_26GSE_V4" / "scanpy_projects"

    @property
    def integrated_h5ad(self) -> Path:
        mirrored = self.root / "high_speed_temp" / "Integrated_dataset" / "integrated.h5ad"
        if mirrored.exists():
            return mirrored
        return self.outputs / "integrated.h5ad"

    @property
    def historical_five_million_atlas(self) -> Path:
        return self.root / "high_speed_temp" / "Integrated_dataset" / "integrated_plus6.h5ad"

    def relative(self, path: str | Path) -> str:
        candidate = Path(path)
        try:
            return candidate.resolve().relative_to(self.root.resolve()).as_posix()
        except ValueError:
            return str(candidate.resolve())
