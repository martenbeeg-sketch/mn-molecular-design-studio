"""Helpers for locating BioEmu2 trajectory previews in a job's artifacts."""
from __future__ import annotations

from pathlib import Path
from typing import Any


BIOEMU2_JOB_TOOLS = {"bioemu2", "bioemu2_sampling"}


def ensemble_trajectories(artifact_paths):
    """Find trajectory/topology pairs in a BioEmu-2 artifact inventory."""
    paths = {Path(value) for value in artifact_paths}
    root = Path("artifacts/bioemu2")
    return sorted(path for path in paths if not path.is_absolute() and ".." not in path.parts
                  and path.is_relative_to(root) and path.suffix.lower() == ".xtc"
                  and path.with_name("topology.cif") in paths)


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
