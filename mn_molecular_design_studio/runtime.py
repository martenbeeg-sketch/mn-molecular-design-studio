"""Runtime paths owned by Molecular Design Studio."""

from __future__ import annotations

import json
import os
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_APP_HOME = Path(
    os.getenv(
        "MN_MOLECULAR_DESIGN_STUDIO_APP_HOME",
        str(Path.home() / ".local" / "share" / "mn-molecular-design-studio"),
    )
).expanduser().resolve()
DEFAULT_TMPDIR = DEFAULT_APP_HOME / "tmp"


def ensure_runtime_home(
    app_home: str | Path = DEFAULT_APP_HOME,
    tmpdir: str | Path = DEFAULT_TMPDIR,
) -> Path:
    home = Path(app_home).expanduser().resolve()
    temporary = Path(tmpdir).expanduser().resolve()
    for directory in (
        home / "workdir" / "runs",
        home / "workdir" / "targets",
        home / "reference_files",
        home / "database",
        temporary,
    ):
        directory.mkdir(parents=True, exist_ok=True)
    return home


def app_home() -> Path:
    return Path(
        os.getenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(DEFAULT_APP_HOME))
    ).expanduser().resolve()


def runs_root() -> Path:
    return Path(
        os.getenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(app_home() / "workdir" / "runs"))
    ).expanduser().resolve()


def reference_root() -> Path:
    configured = os.getenv("MN_MOLECULAR_DESIGN_STUDIO_REFERENCE_DIR", "").strip()
    if configured:
        return Path(configured).expanduser().resolve()
    settings_path = app_home() / "settings.json"
    try:
        settings = json.loads(settings_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        settings = {}
    if isinstance(settings, dict) and str(settings.get("reference_directory") or "").strip():
        return Path(settings["reference_directory"]).expanduser().resolve()
    return (app_home() / "reference_files").expanduser().resolve()
