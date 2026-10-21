"""Helpers for locating BioEmu2 trajectory previews in a job's artifacts."""
from __future__ import annotations

from pathlib import Path
from typing import Any


BIOEMU2_JOB_TOOLS = {"bioemu2", "bioemu2_sampling"}


def trajectory_preview_artifacts(job_tool: str, run_dir: Path, artifact_paths: list[tuple[dict[str, Any], Path]]):
    """Return BioEmu XTC artifacts paired with a sibling topology.cif.

    Require both files to resolve inside the run directory. Artifact manifests
    can contain paths from remote workstations, so a missing or external path
    must not be read as a local preview.
    """
    if str(job_tool or "").casefold() not in BIOEMU2_JOB_TOOLS:
        return []
    root = Path(run_dir).resolve()
    previews = []
    for item, path in artifact_paths:
        trajectory = Path(path)
        if trajectory.suffix.casefold() != ".xtc":
            continue
        topology = trajectory.with_name("topology.cif")
        try:
            trajectory.resolve(strict=True).relative_to(root)
            topology.resolve(strict=True).relative_to(root)
        except (OSError, ValueError):
            continue
        if trajectory.is_file() and topology.is_file():
            previews.append((item, trajectory, topology))
    return previews
