"""Source-evidenced creation times, distinct from run, copy and SQL import times."""
from __future__ import annotations

from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
import re
import socket


SCHEMA = "studio.design-creation.v1"


def utc_timestamp(value):
    try:
        parsed = datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    except (TypeError, ValueError):
        return None
    return parsed.astimezone(timezone.utc) if parsed.tzinfo is not None else None


def sequence_digest(value):
    return hashlib.sha256(str(value or "").encode()).hexdigest()


def creation_key(value):
    return hashlib.sha256(json.dumps({key: value.get(key) for key in (
        "at", "node_id", "workstation", "source_relative_path", "source_sha256", "binder_sequence_sha256")},
        sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def valid_creation(value, sequence=""):
    return bool(isinstance(value, dict) and value.get("schema") == SCHEMA and utc_timestamp(value.get("at"))
        and re.search(r"[A-Za-z0-9]", str(value.get("workstation") or ""))
        and re.fullmatch(r"[a-f0-9]{64}", str(value.get("source_sha256") or "")) and value.get("basis") in {
            "native_cif_creation", "native_structure_written"}
        and (not sequence or value.get("binder_sequence_sha256") == sequence_digest(sequence)))


def read_json(path):
    try:
        value = json.loads(Path(path).read_text())
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError):
        return {}


def workstation_names():
    from mn_molecular_design_studio.cluster import config

    settings = config.load()
    return {config.identity(): config.display_name(settings["name"]).lower(),
        **{peer["id"]: config.display_name(peer.get("name", "")).lower() for peer in settings.get("peers", [])}}


def generating_workstation(run_metadata, names, *, local_node="", local_name=""):
    deferred = str(run_metadata.get("deferred_source_job") or "")
    match = re.match(r"cluster://([^/]+)/", deferred)
    node = str(run_metadata.get("remote_execution_node_id") or (match[1] if match else "")
        or run_metadata.get("execution_node_id") or local_node)
    name = names.get(node) or (local_name if not node or node == local_node else "")
    return node, str(name or "").lower()


def _structure_sequences(path):
    from mn_molecular_design_studio.design_validation import protein_chains, read_structure, sequence

    return {chain: sequence(residues) for chain, residues in protein_chains(read_structure(path)).items()}


def _native_cif_timestamp(path):
    import gemmi

    try:
        document = gemmi.cif.read_file(str(path))
        if len(document) != 1:
            return None
        block = document[0]
        dates = list(block.find_values("_entry.date"))
        times = list(block.find_values("_entry.time"))
        if len(dates) != 1 or len(times) != 1:
            return None
        stamp = datetime.fromisoformat(gemmi.cif.as_string(dates[0]) + "T" + gemmi.cif.as_string(times[0]))
        return stamp.replace(tzinfo=timezone.utc) if stamp.tzinfo is None else stamp.astimezone(timezone.utc)
    except (OSError, ValueError, RuntimeError):
        return None


def _utc_container(run_dir):
    commands = read_json(Path(run_dir) / "command.json")
    text = json.dumps(commands)
    if "/etc/localtime" in text or "/etc/timezone" in text:
        return False
    overrides = re.findall(r"TZ=([^\s\"']+)", text)
    return all(value in {"UTC", "Etc/UTC", "UTC0"} for value in overrides)


def creation_from_structure(path, candidate, run_dir, *, names, local_node, local_name, allow_file_time=False):
    path, run_dir = Path(path), Path(run_dir)
    sequence = str(candidate.get("binder_sequence") or "")
    if not sequence or not path.is_file():
        return None
    metadata = read_json(run_dir / "metadata.json")
    node, workstation = generating_workstation(metadata, names, local_node=local_node, local_name=local_name)
    if not workstation:
        return None
    chains = _structure_sequences(path)
    designed = list(candidate.get("binder_chains") or [])
    if designed:
        if not all(chain in chains for chain in designed) or "".join(chains[chain] for chain in designed) != sequence:
            return None
    elif list(chains.values()).count(sequence) != 1:
        return None
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    tool = str(candidate.get("source_tool") or "")
    stamp = _native_cif_timestamp(path) if tool == "rfd4_proteina" and _utc_container(run_dir) else None
    basis, precision = "native_cif_creation", "second"
    if stamp is None:
        if not allow_file_time or node != local_node or metadata.get("cluster_remote") or metadata.get("deferred_output"):
            return None
        stat = path.stat()
        stamp = datetime.fromtimestamp(stat.st_mtime_ns // 1_000_000_000, timezone.utc).replace(
            microsecond=(stat.st_mtime_ns % 1_000_000_000) // 1000)
        basis, precision = "native_structure_written", "microsecond"
    started = utc_timestamp(metadata.get("created_at"))
    completed = utc_timestamp(metadata.get("completed_at"))
    if started and stamp < started or completed and stamp >= completed + timedelta(seconds=1):
        return None
    deferred = re.match(r"cluster://[^/]+/[^/]+/([^/]+)", str(metadata.get("deferred_source_job") or ""))
    evidence = {"schema": SCHEMA, "at": stamp.isoformat(timespec="seconds" if precision == "second" else "microseconds"),
        "basis": basis, "precision": precision, "workstation": workstation, "node_id": node,
        "source_path": str(path.resolve()), "source_sha256": digest,
        "candidate_id": str(candidate.get("candidate_id") or ""), "binder_sequence_sha256": sequence_digest(sequence),
        "source_run_id": deferred[1] if deferred else str(metadata.get("run_id") or run_dir.name),
        "source_relative_path": str(path.relative_to(run_dir)) if path.is_relative_to(run_dir) else path.name}
    if basis == "native_cif_creation":
        evidence["timestamp_fields"] = ["_entry.date", "_entry.time"]
        evidence["timezone_basis"] = "RFD4 native UTC container; no local timezone override"
    else:
        evidence["source_mtime_ns"] = path.stat().st_mtime_ns
        evidence["limitation"] = "Original design structure write time; exact sequence-generation instant was not separately recorded"
    return evidence


def candidate_creation(run_dir, candidate, *, names, local_node, local_name):
    raw = candidate.get("raw_metadata") or {}
    existing = raw.get("design_creation")
    sequence = str(candidate.get("binder_sequence") or "")
    if valid_creation(existing, sequence):
        return existing
    parent = raw.get("source_candidate") or {}
    parent_creation = (parent.get("raw_metadata") or {}).get("design_creation")
    generative = str(candidate.get("stage") or "").startswith(("generation.", "sequence_design"))
    native_bindcraft = candidate.get("source_tool") in {"bindcraft2", "bindcraft2_generator"}
    if not generative and not native_bindcraft:
        return parent_creation if valid_creation(parent_creation, sequence) else None
    if candidate.get("source_tool") == "bindcraft2" and not raw.get("design_creation_structure"):
        return None
    source = raw.get("native_structure") or raw.get("design_creation_structure") or candidate.get("complex_pdb") or candidate.get("binder_pdb")
    if not source:
        return None
    from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path

    path = resolve_stored_path(source, run_dir=Path(run_dir), must_exist=True)
    if path is None:
        return None
    return creation_from_structure(path, candidate, run_dir, names=names, local_node=local_node,
        local_name=local_name, allow_file_time=True)


def preserve_creation_times(run_dir, candidates):
    from mn_molecular_design_studio.cluster import config

    if not candidates:
        return []
    run_dir = Path(run_dir)
    names = workstation_names()
    local_node = config.identity()
    local_name = names.get(local_node) or socket.gethostname().lower()
    previous = {}
    path = run_dir / "artifacts/normalized_candidates/candidates.jsonl"
    if path.is_file():
        for line in path.read_text().splitlines():
            try:
                row = json.loads(line)
            except ValueError:
                continue
            value = (row.get("raw_metadata") or {}).get("design_creation")
            if valid_creation(value, row.get("binder_sequence") or ""):
                previous[(row.get("candidate_id"), row.get("source_tool"), row.get("binder_sequence"))] = value
    prepared, recovered = [], {}
    for candidate in candidates:
        key = (candidate.get("candidate_id"), candidate.get("source_tool"), candidate.get("binder_sequence"))
        evidence = previous.get(key)
        if evidence is None:
            raw = candidate.get("raw_metadata") or {}
            cache_key = (candidate.get("source_tool"), candidate.get("stage"), candidate.get("binder_sequence"),
                tuple(candidate.get("binder_chains") or []), str(raw.get("native_structure") or raw.get("design_creation_structure")
                or candidate.get("complex_pdb") or candidate.get("binder_pdb") or ""))
            try:
                if valid_creation(raw.get("design_creation"), candidate.get("binder_sequence") or "") or raw.get("source_candidate"):
                    evidence = candidate_creation(run_dir, candidate, names=names, local_node=local_node, local_name=local_name)
                else:
                    if cache_key not in recovered:
                        recovered[cache_key] = candidate_creation(run_dir, candidate, names=names, local_node=local_node, local_name=local_name)
                    evidence = recovered[cache_key]
                    if evidence:
                        evidence = {**evidence, "candidate_id": str(candidate.get("candidate_id") or "")}
            except (OSError, ValueError, KeyError, RuntimeError):
                evidence = None
        if evidence:
            prepared.append({**candidate, "raw_metadata": {**(candidate.get("raw_metadata") or {}), "design_creation": evidence}})
        elif "design_creation" in (candidate.get("raw_metadata") or {}):
            prepared.append({**candidate, "raw_metadata": {key: value for key, value in candidate["raw_metadata"].items() if key != "design_creation"}})
        else:
            prepared.append(candidate)
    return prepared
