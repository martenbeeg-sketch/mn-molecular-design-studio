"""Isolated workspace setup for Molecular Design Studio."""

from __future__ import annotations

import json
import os
import uuid
from urllib.parse import urlsplit
from datetime import datetime, timezone
from pathlib import Path

from mn_molecular_design_studio.cluster.replication import identity_edit

STUDIO_SOURCE_KEY = "molecular_design_studio"
STUDIO_SOURCE_ID = "mn-molecular-design-studio"


def _default_results_dir() -> Path:
    configured_parent = os.environ.get("MN_MOLECULAR_DESIGN_STUDIO_RESULTS_PARENT", "").strip()
    if configured_parent:
        return (
            Path(configured_parent).expanduser()
            / "mn-molecular-design-studio"
            / "workdir"
            / "runs"
        ).resolve()
    shared_results_parent = Path("/mnt/data/RESULTS")
    if shared_results_parent.is_dir():
        return (shared_results_parent / "mn-molecular-design-studio" / "workdir" / "runs").resolve()
    return (_studio_home() / "workdir" / "runs").resolve()


def _default_reference_dir(home: Path | None = None) -> Path:
    shared_reference_root = Path("/mnt/db/reference_files")
    if shared_reference_root.is_dir():
        return shared_reference_root.resolve()
    return ((home or _studio_home()) / "reference_files").resolve()


def load_settings(home: Path | None = None) -> dict:
    workspace = home or _studio_home()
    defaults = {
        "results_directory": str(_default_results_dir()),
        "reference_directory": str(_default_reference_dir(workspace)),
        "api_base_url": "http://127.0.0.1:8766",
    }
    try:
        payload = json.loads((workspace / "settings.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        payload = {}
    if not isinstance(payload, dict):
        payload = {}
    results_directory = str(payload.get("results_directory") or defaults["results_directory"]).strip()
    defaults["results_directory"] = results_directory or str(_default_results_dir())
    reference_directory = str(payload.get("reference_directory") or defaults["reference_directory"]).strip()
    defaults["reference_directory"] = reference_directory or str(_default_reference_dir(workspace))
    defaults["api_base_url"] = str(payload.get("api_base_url") or defaults["api_base_url"]).rstrip("/")
    return defaults


def save_settings(settings: dict, home: Path | None = None) -> dict:
    workspace = home or _studio_home()
    workspace.mkdir(parents=True, exist_ok=True)
    results_directory = Path(str(settings.get("results_directory") or _default_results_dir())).expanduser().resolve()
    reference_directory = Path(
        str(settings.get("reference_directory") or _default_reference_dir(workspace))
    ).expanduser().resolve()
    results_directory.mkdir(parents=True, exist_ok=True)
    if not reference_directory.exists():
        reference_directory.mkdir(parents=True, exist_ok=True)
    elif not reference_directory.is_dir():
        raise ValueError("The model reference files path must be a directory.")
    try:
        existing = json.loads((workspace / "settings.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        existing = {}
    payload = dict(existing) if isinstance(existing, dict) else {}
    api_url = str(settings.get("api_base_url") or payload.get("api_base_url") or "http://127.0.0.1:8766").rstrip("/")
    parsed = urlsplit(api_url)
    if parsed.scheme not in {"http", "https"} or not parsed.hostname or parsed.username or parsed.password or parsed.query or parsed.fragment or parsed.path:
        raise ValueError("Enter an HTTP API address with a hostname and optional port, such as http://127.0.0.1:8766.")
    payload.update(
        {
            "results_directory": str(results_directory),
            "reference_directory": str(reference_directory),
            "api_base_url": api_url,
        }
    )
    path = workspace / "settings.json"
    temporary_path = path.with_suffix(".json.tmp")
    temporary_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    temporary_path.replace(path)
    return payload


@identity_edit
def save_organization(organization: dict, home: Path | None = None) -> dict:
    """Persist organization and project metadata in the Studio-owned workspace."""
    workspace = home or _studio_home()
    path = workspace / "organization" / "organization.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = dict(organization)
    payload.setdefault("schema_version", 1)
    payload.setdefault("id", "local-default")
    payload.setdefault("name", "Molecular Design Studio")
    payload.setdefault("members", [])
    projects = payload.get("projects")
    payload["projects"] = projects if isinstance(projects, list) else []
    temporary_path = path.with_suffix(".json.tmp")
    temporary_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    temporary_path.replace(path)
    return payload


@identity_edit
def create_project(
    name: str,
    description: str = "",
    color: str = "#ef5350",
    home: Path | None = None,
    *,
    owner_id: str | None = None,
) -> dict:
    """Create a named project for grouping jobs from different Studio tools."""
    from mn_molecular_design_studio.cluster import config as cluster_config
    if cluster_config.is_worker():
        from mn_molecular_design_studio.cluster import identity_auth, transport
        project = transport.request("/api/v1/cluster/projects", {"name": name, "description": description, "color": color},
                                    token=identity_auth.current_token(), method="POST")
        organization = transport.request("/api/v1/cluster/organization", token=identity_auth.current_token())
        cluster_config.atomic_json(cluster_config.home() / "organization-cache.json", organization)
        return project
    clean_name = " ".join(str(name or "").split())
    if not clean_name:
        raise ValueError("Enter a project name.")
    organization = load_organization(home)
    projects = organization.get("projects")
    projects = projects if isinstance(projects, list) else []
    if any(
        isinstance(project, dict)
        and str(project.get("name", "")).strip().casefold() == clean_name.casefold()
        for project in projects
    ):
        raise ValueError(f"A project named ‘{clean_name}’ already exists.")

    project = {
        "id": uuid.uuid4().hex,
        "name": clean_name,
        "description": str(description or "").strip(),
        "color": str(color or "#ef5350"),
        "created_by": owner_id,
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    projects.append(project)
    organization["projects"] = projects
    workspace = home or _studio_home()
    save_organization(organization, workspace)
    (workspace / "organization" / "projects" / project["id"]).mkdir(parents=True, exist_ok=True)
    return project


def validate_project_id(project_id: object, home: Path | None = None) -> str | None:
    """Return a known project ID, or reject stale/unknown IDs from API requests."""
    if project_id is None or project_id == "":
        return None
    if not isinstance(project_id, str):
        raise ValueError("project_id must be a project ID or null.")
    organization = load_organization(home)
    projects = organization.get("projects")
    if not isinstance(projects, list) or not any(
        isinstance(project, dict) and project.get("id") == project_id for project in projects
    ):
        raise ValueError(f"Unknown project ID: {project_id}")
    return project_id


def _studio_home() -> Path:
    configured_home = os.environ.get("MN_MOLECULAR_DESIGN_STUDIO_HOME", "").strip()
    if configured_home:
        return Path(configured_home).expanduser().resolve()
    return (Path.home() / ".local" / "share" / "mn-molecular-design-studio").resolve()


def ensure_studio_workspace() -> Path:
    """Create a fresh studio data root and point the shared workflow engine at it."""
    home = _studio_home()
    settings = load_settings(home)
    results_dir = Path(settings["results_directory"]).expanduser().resolve()
    reference_dir = Path(settings["reference_directory"]).expanduser().resolve()

    # The protein-design package supplies workflow code, but the studio owns its
    # targets, references, groups, and organization state in a separate root.
    # Its default run directory is separate from the previous app and configurable in Settings.
    os.environ["MN_MOLECULAR_DESIGN_STUDIO_APP_HOME"] = str(home)
    os.environ["MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR"] = str(results_dir)
    os.environ["MN_MOLECULAR_DESIGN_STUDIO_REFERENCE_DIR"] = str(reference_dir)

    for directory in (
        results_dir,
        home / "workdir" / "targets",
        home / "reference_files",
        home / "database",
        home / "organization" / "projects",
        home / "tmp",
    ):
        directory.mkdir(parents=True, exist_ok=True)
    if not reference_dir.is_dir():
        reference_dir.mkdir(parents=True, exist_ok=True)

    organization_path = home / "organization" / "organization.json"
    if not organization_path.exists():
        organization = {
            "schema_version": 1,
            "id": "local-default",
            "name": "Molecular Design Studio",
            "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "members": [],
            "projects": [],
        }
        temporary_path = organization_path.with_suffix(".json.tmp")
        temporary_path.write_text(json.dumps(organization, indent=2), encoding="utf-8")
        temporary_path.replace(organization_path)

    return home


_PROJECT_ID_UNSET = object()
_OWNER_ID_UNSET = object()


def mark_studio_job(
    result: object,
    *,
    project_id: str | None | object = _PROJECT_ID_UNSET,
    owner_id: str | None | object = _OWNER_ID_UNSET,
) -> object:
    """Tag a newly created run so shared-folder legacy jobs remain distinguishable."""
    if not isinstance(result, dict):
        return result
    if project_id is not _PROJECT_ID_UNSET:
        project_id = validate_project_id(project_id)
        if project_id is None:
            raise ValueError("Every Molecular Design Studio job must belong to a project.")
    job_result = result.get("job") if isinstance(result.get("job"), dict) else result
    if job_result.get("cluster_job_id"):
        # A shared-catalog response may describe another computer's path. Never
        # interpret it as a local run folder or write ownership files beside it.
        return result
    run_dir_text = job_result.get("run_dir")
    if not run_dir_text:
        return result
    run_dir = Path(str(run_dir_text)).expanduser().resolve()
    marker_path = run_dir / "molecular_design_studio.json"
    try:
        existing = json.loads(marker_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        existing = {}
    marker = {
        **(existing if isinstance(existing, dict) else {}),
        "application": STUDIO_SOURCE_ID,
        "workspace": str(_studio_home()),
        "marked_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    if project_id is not _PROJECT_ID_UNSET:
        marker["project_id"] = project_id
    if owner_id is not _OWNER_ID_UNSET:
        marker["studio_owner_id"] = owner_id or None
    temporary_path = marker_path.with_suffix(".json.tmp")
    try:
        temporary_path.write_text(json.dumps(marker, indent=2), encoding="utf-8")
        temporary_path.replace(marker_path)
    except OSError:
        temporary_path.unlink(missing_ok=True)
    else:
        try:
            from mn_molecular_design_studio.job_catalog import reconcile, update_job_project

            update_job_project(run_dir)
            reconcile()
        except Exception:
            # Indexing is rebuildable and must not make job submission fail.
            pass
        from mn_molecular_design_studio.molecular_database_ingest import record_completed_run
        record_completed_run(run_dir)
    return result


def is_studio_job(run_dir: str | Path) -> bool:
    try:
        marker = json.loads((Path(run_dir) / "molecular_design_studio.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        marker = None
    if isinstance(marker, dict) and marker.get("application") == STUDIO_SOURCE_ID:
        return True
    try:
        metadata = json.loads((Path(run_dir) / "metadata.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    source = metadata.get(STUDIO_SOURCE_KEY) if isinstance(metadata, dict) else None
    return isinstance(source, dict) and source.get("application") == STUDIO_SOURCE_ID


def load_organization(home: Path | None = None) -> dict:
    from mn_molecular_design_studio.cluster import config as cluster_config
    if cluster_config.is_worker():
        try:
            return json.loads((cluster_config.home() / "organization-cache.json").read_text())
        except (FileNotFoundError, json.JSONDecodeError):
            return {"id": "cluster", "name": "Workstation cluster", "members": [], "projects": []}
    workspace = home or ensure_studio_workspace()
    path = workspace / "organization" / "organization.json"
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {"id": "local-default", "name": "Molecular Design Studio", "members": [], "projects": []}
    return payload if isinstance(payload, dict) else {}
