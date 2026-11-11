"""Owner-local, durable admission of exact Database assessment work."""
from __future__ import annotations

import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import sqlite3
import tempfile


def _subject_identity(value, *, component=False, retain_component_roles=False):
    if isinstance(value, dict):
        return {key: _subject_identity(child, component=component or key == "database_identity_components",
                                       retain_component_roles=retain_component_roles) for key, child in value.items()
                if key not in {"path", "record_name", "structure_label", "source_run_dir", "candidates_jsonl"}
                and not (component and key == "role" and not retain_component_roles)
                and not key.endswith("_path")}
    if isinstance(value, list):
        return [_subject_identity(child, component=component, retain_component_roles=retain_component_roles) for child in value]
    return value


def _legacy_admissions(root, scope, canonical):
    from mn_molecular_design_studio.api_idempotency import Conflict

    path = root.parent / "api-requests.sqlite3"
    if not path.is_file():
        return []
    with sqlite3.connect(path.as_uri() + "?mode=ro", uri=True) as connection:
        rows = connection.execute("SELECT key,fingerprint,status,response FROM requests WHERE owner=? "
                                  "AND kind='scientific' AND key NOT LIKE 'v2:%'", (scope["id"],)).fetchall()
    matches = []
    for key, fingerprint, status, response in rows:
        identity = hashlib.sha256((scope["id"] + key).encode()).hexdigest()
        manifest_path = root / "inputs" / identity / "complex_analysis_source.json"
        try:
            manifest = json.loads(manifest_path.read_text())
        except (OSError, ValueError) as exc:
            raise Conflict("Legacy scientific admission lacks verifiable inputs; reconcile it before submitting new work.") from exc
        if _subject_identity(manifest) != canonical["source"]:
            continue
        legacy = {**canonical, "source": _subject_identity(manifest, retain_component_roles=True)}
        if hashlib.sha256(json.dumps(legacy, sort_keys=True, allow_nan=False).encode()).hexdigest() != fingerprint:
            continue
        if status != "completed":
            raise Conflict("This exact legacy assessment is pending or uncertain; inspect it before retrying.")
        result = json.loads(response) if response else None
        if not isinstance(result, dict) or not isinstance(result.get("job"), dict) or not result["job"].get("run_id"):
            raise Conflict("Legacy scientific admission has an incomplete response; reconcile it before retrying.")
        matches.append(result)
    return sorted(matches, key=lambda result: str(result["job"]["run_id"]))


def submit_database_assessment(user, project_id, payload, operation, *, validate=None):
    from mn_molecular_design_studio.auth import user_input_root
    root = user_input_root(str(user["id"])) / "complex_analysis_sources"
    root.mkdir(parents=True, exist_ok=True, mode=0o700)
    with (root / ".preparation.lock").open("a+b") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if validate is not None:
            validate()
        return _submit_database_assessment(user, project_id, payload, operation)


def _freeze_inputs(payload, root, identity):
    parameters = payload["parameters"]
    source = Path(parameters["source_run_dir"]).resolve()
    destination = root / "inputs" / identity
    destination.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    temporary = Path(tempfile.mkdtemp(dir=destination.parent))
    def rewrite(value):
        if isinstance(value, dict):
            return {key: rewrite(child) for key, child in value.items()}
        if isinstance(value, list):
            return [rewrite(child) for child in value]
        if isinstance(value, str) and (value == str(source) or value.startswith(str(source) + "/")):
            return str(destination) + value[len(str(source)):]
        return value
    try:
        if any(path.is_symlink() for path in source.rglob("*")):
            raise ValueError("Assessment input bundles must not contain mutable symlinks.")
        shutil.copytree(source, temporary, dirs_exist_ok=True)
        for path in temporary.rglob("*"):
            if not path.is_file():
                continue
            if path.suffix == ".json":
                path.write_text(json.dumps(rewrite(json.loads(path.read_text())), sort_keys=True))
            elif path.suffix == ".jsonl":
                path.write_text("".join(json.dumps(rewrite(json.loads(line)), sort_keys=True) + "\n"
                                        for line in path.read_text().splitlines() if line.strip()))
        manifest = json.loads((temporary / "complex_analysis_source.json").read_text())
        structure = Path(parameters["selected_structure_path"]).resolve().relative_to(source)
        if hashlib.sha256((temporary / structure).read_bytes()).hexdigest() != manifest["structure_digest"]:
            raise ValueError("Assessment coordinates changed before admission.")
        os.rename(temporary, destination)
        parameters.update(rewrite(parameters))
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def _submit_database_assessment(user, project_id, payload, operation):
    from mn_molecular_design_studio.api_idempotency import Conflict, submit_once, submission_status
    from mn_molecular_design_studio.complex_analysis_batch import assessment_protocol
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
    from mn_molecular_design_studio.runtime import app_home
    from mn_molecular_design_studio.scientific_settings import scientific_parameters

    method = payload["workflow"]
    parameters = payload["parameters"]
    manifest = json.loads((Path(parameters["source_run_dir"]) / "complex_analysis_source.json").read_text())
    source = _subject_identity(manifest)
    defaults = {field["name"]: field["default"] for field in workflow_schema(method)["parameters"] if "default" in field}
    canonical = {"workflow": method, "source": source,
                 "settings": scientific_parameters({**defaults, **parameters}, method),
                 "protocol": assessment_protocol(method)}
    fingerprint = hashlib.sha256(json.dumps(canonical, sort_keys=True, allow_nan=False).encode()).hexdigest()
    scope = {"id": "project:" + project_id if project_id else "user:" + str(user["id"])}
    root = Path(app_home()) / "scientific-admission"
    root.mkdir(parents=True, exist_ok=True, mode=0o700)
    lock_name = hashlib.sha256((scope["id"] + fingerprint).encode()).hexdigest()
    with (root / (lock_name + ".lock")).open("a+b") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        retry_of = str((payload.get("execution") or {}).get("reassessment_of") or "")
        key = "v2:" + fingerprint
        retry_authorized = False
        if retry_of:
            path = Path(app_home()) / "api-requests.sqlite3"
            if path.is_file():
                with sqlite3.connect(path) as connection:
                    retry_authorized = bool(connection.execute("SELECT 1 FROM requests WHERE owner=? AND kind='scientific' AND fingerprint=? "
                        "AND status='completed' AND json_extract(response,'$.job.run_id')=?",
                        (scope["id"], fingerprint, retry_of)).fetchone())
            key += ":retry:" + retry_of
        try:
            previous = submission_status(scope, "scientific", key)
        except FileNotFoundError:
            previous = None
        if previous:
            if previous["status"] == "completed":
                return previous["response"]
            raise Conflict("This exact assessment is pending or uncertain; inspect its submission before retrying.")
        legacy = _legacy_admissions(root, scope, canonical)
        if legacy and not retry_of:
            return submit_once(scope, "scientific", key, canonical, lambda: legacy[0])
        retry_authorized = retry_authorized or any(row["job"]["run_id"] == retry_of for row in legacy)
        from mn_molecular_design_studio.runtime import runs_root
        analysis_root = Path(runs_root()) / "analysis"
        existing = None
        if analysis_root.is_dir() and any(analysis_root.iterdir()):
            from mn_molecular_design_studio.complex_analysis_ui import _matching_database_analysis_jobs
            from mn_molecular_design_studio.complex_analysis_batch import analysis_job_coverage, protocol_matches
            matches = _matching_database_analysis_jobs(user, workflow_id=method, source=manifest,
                candidate_id=manifest["candidate_id"], model_id=str(manifest["model_id"]),
                binder_chains=list(manifest.get("binder_chains") or []),
                target_chains=list(manifest.get("target_chains") or []),
                parameters=parameters, ligand_selection=manifest.get("ligand_selection"),
                load_evidence=False, expected_protocol={} if retry_of else assessment_protocol(method))
            retry_authorized = retry_authorized or any(row["run_id"] == retry_of and row["status"] in {
                "completed", "succeeded", "failed", "cancelled", "canceled", "error"} for row in matches)
            blocking = analysis_job_coverage(matches)
            if blocking and blocking[0] != "Previously failed":
                existing = next(row for row in matches if row.get("status") not in {"failed", "cancelled", "canceled", "error"})
                if retry_of and any(row["status"] not in {"completed", "succeeded", "failed", "cancelled", "canceled", "error"}
                                    for row in matches):
                    raise Conflict("Matching scientific work is still active; do not duplicate it through a retry.")
                if not retry_of and existing["status"] in {"completed", "succeeded"} and not protocol_matches(
                        (existing.get("database_analysis") or {}).get("protocol") or {}, assessment_protocol(method)):
                    raise Conflict("Historical work has unverified protocol identity; review it before an explicit reassessment.")
                if retry_of:
                    existing = None
        if retry_of and not retry_authorized:
            raise Conflict("Reassessment must identify a previous admitted job for this exact scientific request.")
        if existing:
            response = {"job": {key: str(value) if isinstance(value, Path) else value
                                for key, value in existing.items() if key != "database_analysis"}, "reused": True}
            operation = lambda: response
        else:
            submit = operation
            def operation():
                identity = hashlib.sha256((scope["id"] + key).encode()).hexdigest()
                _freeze_inputs(payload, root, identity)
                return submit()
        return submit_once(scope, "scientific", key, canonical, operation)
