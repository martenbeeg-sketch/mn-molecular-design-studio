"""Recover creation evidence from retained workdir outputs without editing them."""
from __future__ import annotations

from collections import defaultdict
import csv
import json
from pathlib import Path

from mn_molecular_design_studio.design_creation import (
    creation_from_structure, creation_key, read_json, valid_creation,
)
from mn_molecular_design_studio.library_curation import _design_group, _digest, encoded
from mn_molecular_design_studio.scientific_naming import design_name, named_metadata, technical_design_name


def _design_sequence(components):
    designed = [component for component in components if component.get("role") in {"binder", "designed_binder", "design"}]
    if not designed and len(components) == 1 and components[0].get("type") == "protein":
        designed = components
    return "".join(component.get("sequence") or "" for component in designed)


class CreationRecovery:
    def __init__(self, *, names, local_node, local_name, runs_root):
        self.names, self.local_node, self.local_name = names, local_node, local_name
        self.runs_root = Path(runs_root).resolve()
        self.cache = {}
        self.candidates = {}

    def _run(self, value):
        text = str(value or "")
        if not text:
            return None
        path = Path(text)
        candidate = self.runs_root / "design" / path.name
        if candidate.is_dir():
            return candidate
        return path.resolve() if path.is_dir() and path.resolve().is_relative_to(self.runs_root) else None

    def _path(self, source, run):
        if not source or run is None:
            return None
        text = str(source)
        relative = text.split("/artifacts/", 1)[1] if "/artifacts/" in text else None
        path = run / "artifacts" / relative if relative is not None else run / text
        if path.is_file() and path.resolve().is_relative_to(run.resolve()):
            return path.resolve()
        return None

    def _evidence(self, path, candidate, run, *, file_time=False):
        if path is None or run is None:
            return None
        key = (str(path), candidate.get("binder_sequence"), tuple(candidate.get("binder_chains") or []), str(run), file_time)
        if key not in self.cache:
            try:
                self.cache[key] = creation_from_structure(path, candidate, run, names=self.names,
                    local_node=self.local_node, local_name=self.local_name, allow_file_time=file_time)
            except (OSError, ValueError, KeyError, RuntimeError):
                self.cache[key] = None
        return self.cache[key]

    def recover(self, metadata, origin, components):
        sequence = _design_sequence(components)
        if not sequence:
            return None, "No unambiguous designed sequence"
        existing = metadata.get("design_creation")
        if valid_creation(existing, sequence):
            return existing, ""
        naming = metadata.get("studio_naming") or {}
        tool = metadata.get("generator_engine") or naming.get("tool") or metadata.get("source_tool") or origin.get("tool")
        native_id = str(metadata.get("native_candidate_id") or metadata.get("source_candidate_id") or naming.get("native_candidate_id") or "")
        run_values = [metadata.get("source_run_dir"), metadata.get("source_run_id"), origin.get("run_dir"), origin.get("run_id")]
        runs = list(dict.fromkeys(run for value in run_values if (run := self._run(value)) is not None))
        base = {"candidate_id": native_id, "source_tool": tool, "binder_sequence": sequence}
        if tool == "bindcraft2":
            for run in runs:
                table = run / "artifacts/raw/bindcraft2/output/2_Refolded/!_Refolded.csv"
                if not table.is_file():
                    continue
                with table.open(newline="", encoding="utf-8-sig") as handle:
                    rows = [row for row in csv.DictReader(handle) if "".join(str(row.get("Binder_Sequence") or "").split()) == sequence]
                if len(rows) != 1:
                    continue
                path = table.parent / "Complexes" / (Path(str(rows[0].get("design") or "")).name + ".cif")
                evidence = self._evidence(path, base, run, file_time=True)
                if evidence:
                    return evidence, ""
            return None, "No unique sequence-matched original BindCraft2 output on its generating workstation"
        if tool != "rfd4_proteina":
            return None, "No reviewed native creation-time recovery for this generator"
        sources = [metadata.get("native_structure")]
        for run in runs:
            for source in sources:
                path = self._path(source, run)
                evidence = self._evidence(path, base, run)
                if evidence:
                    return evidence, ""
        aliases = metadata.get("source_run_aliases") or []
        for alias in aliases:
            run = self._run(alias.get("run_dir") or alias.get("run_id"))
            path = self._path(alias.get("complex_pdb"), run)
            evidence = self._evidence(path, {**base, "binder_chains": alias.get("binder_chains") or []}, run)
            if evidence:
                return evidence, ""
        for run in runs:
            if run not in self.candidates:
                path = run / "artifacts/normalized_candidates/candidates.jsonl"
                self.candidates[run] = [json.loads(line) for line in path.read_text().splitlines() if line.strip()] if path.is_file() else []
            candidates = [row for row in self.candidates[run] if row.get("binder_sequence") == sequence]
            exact = [row for row in candidates if row.get("candidate_id") == native_id or (
                row.get("raw_metadata") or {}).get("native_candidate_id") == native_id]
            if exact:
                candidates = exact
            if len(candidates) != 1:
                continue
            candidate = candidates[0]
            raw = candidate.get("raw_metadata") or {}
            parent = raw.get("source_candidate") or {}
            if parent.get("binder_sequence") == sequence:
                parent_raw = parent.get("raw_metadata") or {}
                raw = {**raw, **parent_raw}
            source_run = self._run(raw.get("source_run_dir")) or run
            path = self._path(raw.get("native_structure") or candidate.get("complex_pdb"), source_run)
            evidence = self._evidence(path, {**base, "binder_chains": candidate.get("binder_chains") or []}, source_run)
            if evidence:
                return evidence, ""
        return None, "No sequence-verified native creation stamp in retained original outputs"


def plan_creation_names(database, user, recovery):
    if user.get("role") != "admin":
        raise PermissionError("Library naming maintenance requires an administrator.")
    groups = {row["id"]: row for row in database.groups(user, scope="all", collection="all") if not row["read_only"] and _design_group(row)}
    with database.connection() as connection:
        members = [dict(row) for row in connection.execute("SELECT * FROM molecule_members WHERE group_id IN (SELECT value FROM json_each(?))",
            (encoded(list(groups)),))]
    prepared, unresolved = [], []
    for member in members:
        metadata = json.loads(member["metadata"])
        naming = metadata.get("studio_naming") or {}
        if naming.get("mode") == "manual" or not (technical_design_name(member["name"]) or member["name"] == naming.get("display_name")):
            continue
        evidence, reason = recovery.recover(metadata, groups[member["group_id"]]["origin"], json.loads(member["components"]))
        if evidence is None:
            unresolved.append({"group_id": member["group_id"], "molecule_id": member["molecule_id"], "reason": reason})
            continue
        metadata = {**metadata, "design_creation": evidence}
        metadata.pop("design_name_suffix", None)
        label = design_name(member["name"], target_name="", tool="", run_id="", metadata=metadata)
        prepared.append((member, metadata, label, creation_key(evidence)))
    collisions = defaultdict(set)
    for member, metadata, label, key in prepared:
        collisions[label].add(key)
    changes = []
    for member, metadata, label, key in prepared:
        if len(collisions[label]) > 1:
            metadata["design_name_suffix"] = key[:10]
            label = design_name(member["name"], target_name="", tool="", run_id="", metadata=metadata)
        previous = metadata.get("studio_naming") or {}
        updated = named_metadata(member["name"], label, metadata, target_name=metadata.get("design_target_name") or previous.get("target_name") or "",
            target_id=metadata.get("design_target_molecule_id") or previous.get("target_molecule_id") or "",
            tool=previous.get("tool") or metadata.get("source_tool") or "", run_id=previous.get("run_id") or "")
        if member["name"] != label or json.loads(member["metadata"]) != updated:
            changes.append({"group_id": member["group_id"], "molecule_id": member["molecule_id"],
                "before_name": member["name"], "after_name": label, "before_metadata": member["metadata"], "after_metadata": encoded(updated)})
    plan = {"schema": "studio.library-design-creation.v1", "database": str(database.path.resolve()), "members": changes,
        "observations": [], "targets": [], "unresolved": unresolved}
    plan["digest"] = _digest(plan)
    return plan
