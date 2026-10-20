"""Durable docking forms and their linked binding-site predictions."""
from __future__ import annotations

import json
import re
import uuid
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from urllib.parse import urlsplit

from mn_molecular_design_studio.auth import user_input_root


WORKFLOWS = {
    "ligand.docking.vina": "vina", "ligand.docking.autodock_gpu": "autodock-gpu",
    "ligand.docking.gnina": "gnina",
    "ligand.docking.unidock_pro": "unidock-pro", "ligand.docking.rosetta_ligand": "rosetta",
}
_VINA_STATE = re.compile(
    r"vina-(?:target-path|normalized-target-path|compound-input|ligand-format|protein-format|pdb-code|"
    r"default-target-signature|derived-box-signature|box-(?:center|size)-[xyz]|"
    r"boltz-selected-chains|boltz-additional-sequences|probe-(?:rows|values|dismissed|initialized|next-id|\d+|distance)|"
    r"pocket-chain|residue-text-.+|pending-prediction|applied-prediction|prediction-overlay-pdb|applied-residue-padding)"
)
_POCKET_STATE = {
    "pocket-mode", "box-ligand", "selected-residues", "applied-residues", "ligand-padding", "residue-padding",
}


def state_key(key: str, workflow: str) -> bool:
    tool = WORKFLOWS[workflow]
    if _VINA_STATE.fullmatch(key):
        return True
    if key.startswith("ligand-vina-") and key.removeprefix("ligand-vina-") in _POCKET_STATE:
        return True
    prefix = f"ligand-{tool}-"
    if key.startswith(prefix):
        return key.removeprefix(prefix) not in {"run", "alignment-ligand", "reference-ligand"}
    return key == "ligand-tool-compound-set-ligand.docking.vina"


def snapshot(state, workflow: str) -> dict:
    result = {}
    for key, value in state.items():
        if isinstance(key, str) and state_key(key, workflow):
            try:
                result[key] = json.loads(json.dumps(value, allow_nan=False))
            except (TypeError, ValueError):
                continue
    return result


def _directory(owner_id: str, project_id: str, workflow: str) -> Path:
    if not re.fullmatch(r"[\w-]{1,128}", project_id, flags=re.ASCII):
        raise ValueError("Select a valid project for the docking preparation.")
    return user_input_root(owner_id) / "docking_preparations" / project_id / WORKFLOWS[workflow]


def save(owner_id: str, project_id: str, workflow: str, preparation_id: str, state: dict) -> dict:
    preparation_id = str(uuid.UUID(preparation_id))
    directory = _directory(owner_id, project_id, workflow)
    directory.mkdir(parents=True, exist_ok=True, mode=0o700)
    path = directory / f"{preparation_id}.json"
    values = snapshot(state, workflow)
    try:
        previous = json.loads(path.read_text())
    except FileNotFoundError:
        previous = {}
    name = str(values.get(f"ligand-{WORKFLOWS[workflow]}-job-name") or "Untitled preparation")
    if previous.get("state") == values and previous.get("name") == name:
        return previous
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    record = {
        "schema_version": 1, "id": preparation_id, "owner_id": owner_id,
        "project_id": project_id, "workflow": workflow, "name": name,
        "created_at": previous.get("created_at", now), "updated_at": now, "state": values,
    }
    temporary = path.with_name(f".{preparation_id}-{uuid.uuid4().hex}.tmp")
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            temporary.chmod(0o600)
            json.dump(record, handle, allow_nan=False)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
    return record


def list_saved(owner_id: str, project_id: str, workflow: str) -> list[dict]:
    records = []
    for path in _directory(owner_id, project_id, workflow).glob("*.json"):
        try:
            row = json.loads(path.read_text())
            if (row.get("schema_version") == 1 and row.get("owner_id") == owner_id
                    and row.get("project_id") == project_id and row.get("workflow") == workflow
                    and isinstance(row.get("state"), dict) and str(uuid.UUID(row["id"])) == path.stem):
                records.append(row)
        except (OSError, ValueError, TypeError, KeyError):
            continue
    return sorted(records, key=lambda row: (row.get("updated_at", ""), row["id"]), reverse=True)


def artifact_relative_path(value: str, job: dict) -> str:
    """Normalize a saved prediction path within its owning run only."""
    text = str(value or "")
    parsed = urlsplit(text)
    if parsed.scheme in {"runs", "cluster"}:
        if parsed.scheme == "runs":
            parts = (parsed.netloc + parsed.path).lstrip("/").split("/")
        else:
            if parsed.netloc != job.get("node_id"):
                raise ValueError("The prediction artifact belongs to a different workstation.")
            parts = parsed.path.lstrip("/").split("/")
        if parts[:2] != [job.get("task_group"), job.get("run_id")]:
            raise ValueError("The prediction artifact belongs to a different run.")
        text = "/".join(parts[2:])
    elif parsed.scheme:
        raise ValueError("Unsupported prediction artifact reference.")
    elif PurePosixPath(text).is_absolute():
        location = job.get("location") or job.get("run_dir") or ""
        try:
            text = str(PurePosixPath(text).relative_to(PurePosixPath(location)))
        except ValueError as exc:
            raise ValueError("The prediction artifact is outside its run.") from exc
    path = PurePosixPath(text)
    if not path.parts or path.is_absolute() or any(part in {"..", "."} for part in path.parts):
        raise ValueError("Invalid prediction artifact path.")
    if path.suffix.lower() not in {".pdb", ".cif", ".mmcif"}:
        raise ValueError("The prediction artifact must be a molecular structure.")
    return str(path)


def pdb_text(text: str) -> str:
    """Convert a Boltz mmCIF result for contact mapping onto the original PDB."""
    if not text.lstrip().startswith("data_"):
        return text
    import io
    from Bio.PDB import MMCIFParser, PDBIO
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import _mmcif_text_with_default_occupancy

    structure = MMCIFParser(QUIET=True).get_structure("prediction", io.StringIO(_mmcif_text_with_default_occupancy(text)))
    output = io.StringIO()
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(output)
    return output.getvalue()
