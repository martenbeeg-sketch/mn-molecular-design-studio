"""Seed and serve a bounded, independent teaching pilot inside its container."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import secrets
import signal
import subprocess
import sys
import threading


def seed(state):
    if state != Path("/state") or os.environ.get("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME") != "/state/home":
        raise RuntimeError("Run only in the dedicated container with /state/home, never against production.")
    if os.environ.get("MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE") != "production":
        raise RuntimeError("The teaching pilot requires real account authentication.")
    os.umask(0o077)
    credentials = state / "credentials.json"
    if credentials.exists():
        return
    from mn_molecular_design_studio import auth, workspace
    from mn_molecular_design_studio.cluster import config
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    from mn_molecular_design_studio.molecule_imports import parse_structure

    workspace.save_settings({"results_directory": "/state/runs", "reference_directory": "/state/references",
                             "api_base_url": "http://127.0.0.1:18766"})
    workspace.ensure_studio_workspace()
    config.save({"mode": "disabled", "share_molecular_database": False, "accept_jobs": False})
    password = secrets.token_urlsafe(24)
    instructor = auth.create_initial_admin("teaching_instructor", "Teaching pilot instructor", password)
    project = workspace.create_project("TEACHING PILOT - disposable instructional data", owner_id=instructor["id"])
    accounts = [{"username": instructor["username"], "password": password, "role": "admin"}]
    for name in ("teaching_student_one", "teaching_student_two"):
        password = secrets.token_urlsafe(24)
        student = auth.create_user(name, name.replace("_", " ").title(), password, actor=instructor)
        auth.grant_project_access(project["id"], student["id"], actor=instructor)
        accounts.append({"username": name, "temporary_password": password, "role": "user"})
    database = MolecularDatabase()
    database.import_records(instructor, [
        {"name": "Teaching ethanol", "components": [{"type": "small_molecule", "chain": "L", "smiles": "CCO"}]},
        {"name": "Teaching ethylamine", "components": [{"type": "small_molecule", "chain": "L", "smiles": "CCN"}]},
    ], name="TEACHING EXAMPLE - two small molecules", project_id=project["id"], origin={"teaching_only": True})
    source = Path(__file__).resolve().parents[1] / "mn_molecular_design_studio/assets/rfdiffusion3/5o45_cropped.pdb"
    records = parse_structure(source.read_bytes(), source.name)
    database.import_records(instructor, records, name="TEACHING EXAMPLE - retained public structure, not a new prediction",
                            project_id=project["id"], origin={"teaching_only": True, "source_asset": source.name})
    credentials.write_text(json.dumps({"purpose": "isolated test accounts, not real student identities",
                                       "project_id": project["id"], "accounts": accounts}, indent=2) + "\n")
    credentials.chmod(0o600)
    (state / "README.txt").write_text(
        "Isolated teaching pilot. No production data, peers, GPUs or Docker socket are available.\n"
        "UI: http://127.0.0.1:18502 ; API: http://127.0.0.1:18766\n"
        "credentials.json is private; student accounts require a password change.\n"
        "Only browsing, imports, native CPU operations and exports are in scope.\n"
        "Container-wide limits are not per-student quotas. Do not use for important research.\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed-only", action="store_true")
    options = parser.parse_args()
    seed(Path("/state"))
    if options.seed_only:
        return
    stopped = threading.Event()
    for signum in (signal.SIGTERM, signal.SIGINT):
        signal.signal(signum, lambda *arguments: stopped.set())
    children = [
        subprocess.Popen([sys.executable, "-m", "mn_molecular_design_studio.cli", "api", "--host", "0.0.0.0", "--port", "18766"]),
        subprocess.Popen([sys.executable, "-m", "streamlit", "run",
            str(Path(__file__).resolve().parents[1] / "mn_molecular_design_studio/run_app.py"),
            "--server.address=0.0.0.0", "--server.port=18502", "--server.headless=true",
            "--server.runOnSave=false", "--server.fileWatcherType=none"]),
    ]
    try:
        while not stopped.wait(1):
            if any(child.poll() is not None for child in children):
                raise RuntimeError("A teaching service exited; stopping the pilot rather than leaving a partial service.")
    finally:
        for child in children:
            if child.poll() is None:
                child.terminate()
        for child in children:
            try:
                child.wait(timeout=15)
            except subprocess.TimeoutExpired:
                child.kill()
                child.wait()


if __name__ == "__main__":
    main()
