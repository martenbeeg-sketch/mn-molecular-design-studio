"""Conservative, provenance-preserving names for admitted scientific designs."""
from __future__ import annotations

import json
import re

from mn_molecular_design_studio.database_browser_labels import readable_name
from mn_molecular_design_studio.design_creation import creation_key, utc_timestamp, valid_creation


POLICY = "studio.scientific-names.v2"


def meaningful_name(value):
    text = str(value or "").strip()
    return bool(readable_name(text) and text.casefold() not in {
        "receptor", "protein", "binder", "molecule", "unnamed ligand", "sequence-defined target",
    } and not re.fullmatch(r"(?:[a-z]+[-_])?[a-f0-9]{12,}", text, re.I))


def technical_design_name(value):
    text = str(value or "")
    return not meaningful_name(text) or bool(re.match(r"^(?:rfd\d|bindcraft\d*|boltzgen|dc)_\d", text, re.I))


def design_name(native_name, *, target_name, tool, run_id, metadata=None):
    metadata = metadata or {}
    creation = metadata.get("design_creation") or {}
    if not valid_creation(creation):
        return ""
    stamp = utc_timestamp(creation["at"])
    timestamp = stamp.strftime("%Y%m%dT%H%M%S")
    if creation.get("precision") == "microsecond":
        timestamp += "." + stamp.strftime("%f")
    workstation = re.sub(r"[^A-Za-z0-9-]+", "-", creation["workstation"]).strip("-").lower()
    suffix = str(metadata.get("design_name_suffix") or "")
    if suffix and not re.fullmatch(r"[a-f0-9]{10}", suffix):
        return ""
    return "Design_" + workstation + "_" + timestamp + "Z" + ("_" + suffix if suffix else "")


def named_metadata(name, label, metadata, *, target_name, target_id, tool, run_id):
    previous = metadata.get("studio_naming") or {}
    native_names = sorted(set(previous.get("native_names") or []) | (set() if name == previous.get("display_name") else {str(name)}))
    return {**metadata, "studio_naming": {
        "policy": POLICY, "display_name": label, "native_names": native_names,
        "native_candidate_id": str(metadata.get("native_candidate_id") or previous.get("native_candidate_id") or name),
        "target_name": str(target_name), "target_molecule_id": str(target_id),
        "tool": str(tool), "run_id": str(run_id),
        "creation_key": creation_key(metadata["design_creation"]) if valid_creation(metadata.get("design_creation")) else "",
        "created_at": (metadata.get("design_creation") or {}).get("at", ""),
        "generating_workstation": (metadata.get("design_creation") or {}).get("workstation", ""),
    }}


def prepare_design_name(connection, group_id, record, metadata, origin):
    native = str(record.get("name") or "")
    if not technical_design_name(native) or (metadata.get("studio_naming") or {}).get("mode") == "manual":
        return native, metadata
    group = connection.execute("SELECT origin FROM molecule_groups WHERE id=?", (group_id,)).fetchone()
    group_origin = json.loads(group[0]) if group else {}
    role = str(origin.get("workflow_role") or metadata.get("workflow_role") or "")
    admitted_design = group_origin.get("collection_kind") == "generated" or (
        group_origin.get("collection_kind") == "curated" and origin.get("task_group") == "design" and origin.get("run_id"))
    if not admitted_design or role in {"post_analysis", "independent_prediction", "docking"}:
        return native, metadata
    components = record.get("components") or []
    designed = [component for component in components if component.get("role") in {"binder", "designed_binder", "design"}]
    if not designed and len(components) == 1 and components[0].get("type") == "protein":
        designed = components
    sequence = "".join(component.get("sequence") or "" for component in designed)
    if sequence and not valid_creation(metadata.get("design_creation"), sequence):
        return native, metadata
    target_id = metadata.get("design_target_molecule_id") or origin.get("design_target_molecule_id") or ""
    target_name = metadata.get("design_target_name") or origin.get("design_target_name") or ""
    if not meaningful_name(target_name) and target_id:
        names = connection.execute("SELECT name FROM molecule_members WHERE molecule_id=? ORDER BY name", (target_id,)).fetchall()
        target_name = next((row[0] for row in names if meaningful_name(row[0])), "")
    tool = metadata.get("generator_engine") or metadata.get("source_tool") or origin.get("tool") or ""
    run_id = str(origin.get("run_id") or "")
    label = design_name(native, target_name=target_name, tool=tool, run_id=run_id, metadata=metadata)
    if label:
        collision = connection.execute("SELECT metadata FROM molecule_members WHERE name=?", (label,)).fetchall()
        own_key = creation_key(metadata["design_creation"])
        if any((json.loads(row[0]).get("studio_naming") or {}).get("creation_key") != own_key for row in collision):
            metadata = {**metadata, "design_name_suffix": own_key[:10]}
            label = design_name(native, target_name=target_name, tool=tool, run_id=run_id, metadata=metadata)
    return (label, named_metadata(native, label, metadata, target_name=target_name, target_id=target_id,
        tool=tool, run_id=run_id)) if label else (native, metadata)
