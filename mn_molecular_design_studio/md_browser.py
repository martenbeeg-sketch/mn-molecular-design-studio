"""SQL-scoped MD studies and bounded, source-qualified analysis artifacts."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re

from mn_molecular_design_studio.md_results_ui import _safe_artifact_path, _safe_relative_parts
from mn_molecular_design_studio.database_browser_labels import tool_label


MD_FIELDS = (
    "study_id", "compound_id", "dataset_name", "md_engine", "settings.md_engine",
    "target_receptor", "selected_complex_run_id", "source_pose_id", "source_tool", "source_production_run_id",
    "replica_index", "continuation_of_run_id", "production_status", "source_production_status",
    "source_summary_relative_path", "source_summary_sha256", "source_task_group", "source_run_id",
    "analysis_run_id", "endpoint_run_id", "source_created_at", "source_completed_at", "study_manifest_sha256",
)
MAX_ARTIFACT_BYTES = 32 * 1024 * 1024


def production_result(result):
    return result.get("tool") == "md.production" or result.get("evaluation", {}).get("evaluator") == "md.production"


def study_summary(originals):
    production = {}
    engines, pose_methods = set(), set()
    analyses = {}
    for result in originals:
        context = result.get("md") or {}
        engine = context.get("md_engine") or context.get("settings", {}).get("md_engine")
        if engine:
            engines.add({"openmm": "OpenMM", "gromacs": "GROMACS"}.get(str(engine).casefold(), str(engine)))
        if context.get("source_tool"):
            pose_methods.add(tool_label(context["source_tool"]))
        if production_result(result):
            identity = (context.get("source_production_run_id") or result.get("run_id") or result["id"],
                        str(context.get("replica_index")))
            production.setdefault(identity, result)
        elif result.get("tool", "").startswith("md."):
            analyses.setdefault(tool_label(result["tool"]), set()).add(result.get("observation_id") or result["id"])
    replicas = {str(result["md"]["replica_index"]) for result in production.values()
                if result.get("md", {}).get("replica_index") is not None}
    statuses = {str(result.get("md", {}).get("production_status") or "not recorded").casefold()
                for result in production.values()}
    return {"MD engine": " / ".join(sorted(engines)) or "Engine not recorded",
        "Starting pose method": " / ".join(sorted(pose_methods)) or None,
        "Production": "Recorded" if production else "Not recorded",
        "Production segments": len(production), "Replicas": len(replicas) if replicas else None,
        "Production status": " / ".join(sorted(statuses)) if statuses else "Not recorded",
        "Saved analyses": "; ".join(f"{label} ({len(identifiers)})" for label, identifiers in sorted(analyses.items())) or "Not recorded"}


def study_identity(index, result):
    run = index["runs"].get(result["run_key"], {})
    context = result.get("md") or {}
    return (run.get("node", ""), run.get("owner_id", ""), run.get("project_id", ""),
        result["molecule_id"], context.get("study_id") or result["run_key"],
        context.get("selected_complex_run_id", ""), result.get("target_pdb_id", ""))


def study_label(result):
    context = result.get("md") or {}
    reference = context.get("dataset_name") == "MD reference controls"
    engine = context.get("md_engine") or context.get("settings", {}).get("md_engine") or ""
    engine = {"openmm": "OpenMM", "gromacs": "GROMACS"}.get(str(engine).casefold(), engine)
    return " · ".join(str(value) for value in (
        "Reference complex" if reference else context.get("compound_id") or result["name"], context.get("target_receptor"),
        result.get("target_pdb_id"), "Prepared reference" if reference else context.get("source_tool"), engine,
    ) if value)


def study_observations(database, user, outcome):
    accessible = database._visible_ids(user, scope="all", collection="all")
    identifiers = {row["observation_id"] for row in outcome["originals"] if row.get("observation_id")}
    with database.connection() as connection:
        rows = connection.execute(
            "SELECT * FROM molecule_observations WHERE molecule_id=? "
            "AND id IN (SELECT value FROM json_each(?)) AND group_id IN (SELECT value FROM json_each(?))",
            (outcome["molecule_id"], json.dumps(sorted(identifiers)), json.dumps(accessible))).fetchall()
    return [{**dict(row), "context": json.loads(row["context"]), "metrics": json.loads(row["metrics"])}
        for row in rows if row["namespace"].startswith("md.") or
        json.loads(row["context"]).get("workflow_role") == "molecular_dynamics"]


def report_choices(observations, *, replica=None, category="production", production_id=None):
    selected = {}
    for observation in observations:
        context = observation["context"]
        if replica is not None and str(context.get("replica_index")) != str(replica):
            continue
        if production_id and context.get("source_production_run_id") != production_id:
            continue
        namespace = observation["namespace"]
        matches = (namespace == "md.production" or context.get("workflow_role") == "molecular_dynamics") if category == "production" else (
            namespace == "md.analysis" if category == "analysis" else namespace.startswith("md.mmgbsa"))
        if not matches or not context.get("source_summary_relative_path"):
            continue
        key = (context.get("_replica_node") or context.get("source_node"), context.get("task_group"),
            context.get("run_id"), context["source_summary_relative_path"], context.get("source_summary_sha256"))
        selected.setdefault(key, observation)
    return sorted(selected.values(), key=lambda row: (
        str(row["context"].get("source_completed_at") or row["created_at"]), row["id"]), reverse=True)


def report_label(observation):
    context = observation["context"]
    metrics = observation["metrics"]
    duration = metrics.get("production_actual_duration_ns")
    return " · ".join(str(value) for value in (
        "Replica " + str(context.get("replica_index") or "not recorded"),
        f"{duration:g} ns segment" if isinstance(duration, (int, float)) else observation["namespace"],
        str(context.get("source_created_at") or observation["created_at"])[:19].replace("T", " "),
        str(context.get("source_production_run_id") or context.get("analysis_run_id") or
            context.get("endpoint_run_id") or observation["run_id"])[-8:],
    ) if value)


def cache_key(source, relative):
    return hashlib.sha256(json.dumps([*source, relative], separators=(",", ":")).encode()).hexdigest()


def refresh_requests(state, prefix, status_loader):
    for key in list(state):
        if str(key).startswith(prefix) and status_loader(state[key]).get("status") == "failed":
            state.pop(key, None)


class MDArtifactReader:
    def __init__(self, database, user, observation, run, *, fetch_text=None):
        self.database, self.user, self.observation, self.run = database, user, observation, run
        self.fetch_text = fetch_text
        self.source = (str(run.get("node") or ""), str(run.get("task_group") or ""), str(run.get("run_id") or ""))
        self.cache_root = database.root / "md-browser-cache"

    def _authorize(self):
        accessible = self.database._visible_ids(self.user, scope="all", collection="all")
        if self.observation["group_id"] not in accessible:
            raise PermissionError("This MD study is no longer accessible.")
        with self.database.connection() as connection:
            row = connection.execute("SELECT context FROM molecule_observations WHERE id=? AND molecule_id=? AND group_id=?",
                (self.observation["id"], self.observation["molecule_id"], self.observation["group_id"])).fetchone()
        if row is None or json.loads(row["context"]) != self.observation["context"]:
            raise ValueError("The saved MD evidence changed. Refresh the SQL view.")

    def text(self, relative, description="MD artifact", *, expected=""):
        self._authorize()
        parts = _safe_relative_parts(relative)
        relative = "/".join(parts)
        if Path(relative).suffix.lower() not in {".json", ".csv", ".pdb", ".cif", ".sdf"}:
            raise ValueError("Only small analysis reports and structure snapshots can be opened here.")
        if expected and not re.fullmatch(r"[0-9a-f]{64}", expected):
            raise ValueError("MD artifact has an invalid recorded SHA-256 identity.")
        key = cache_key(self.source, relative)
        cached = self.cache_root / (key + ".json")
        payload = None
        if cached.is_file():
            if not cached.resolve().is_relative_to(self.cache_root.resolve()) or cached.stat().st_size > 4096:
                raise ValueError("MD artifact cache metadata is invalid.")
            metadata = json.loads(cached.read_text())
            if not isinstance(metadata, dict) or metadata.get("source") != list(self.source) or metadata.get("path") != relative:
                raise ValueError("MD artifact cache source identity mismatch.")
            path = self.cache_root / (key + ".data")
            payload = self._bounded_bytes(path, self.cache_root)
            if hashlib.sha256(payload).hexdigest() != metadata.get("sha256"):
                raise ValueError("Cached MD artifact digest mismatch.")
        else:
            from mn_molecular_design_studio.cluster import config
            if self.source[0] == config.identity():
                from mn_molecular_design_studio.auth import require_job_access
                from mn_molecular_design_studio.job_catalog import get_job_access_record
                from mn_molecular_design_studio.runtime import runs_root
                from mn_molecular_design_studio.artifact_integrity import retained_artifact_file

                access = get_job_access_record(*self.source[1:])
                if access is None or any(str(access.get(field) or "") != str(self.run.get(field) or "")
                    for field in ("owner_id", "project_id")):
                    raise PermissionError("The producing Results job is unavailable or its ownership changed.")
                require_job_access(self.user, access)
                root = Path(runs_root()).resolve()
                run_root = (root / self.source[1] / self.source[2]).resolve()
                if not run_root.is_relative_to(root):
                    raise ValueError("MD run path escapes managed Results storage.")
                path = _safe_artifact_path(run_root, relative)
                if path.stat().st_size > MAX_ARTIFACT_BYTES:
                    raise ValueError("MD report exceeds the inline size limit.")
                payload = self._bounded_bytes(retained_artifact_file(run_root, path))
            elif self.fetch_text is not None:
                text = self.fetch_text(self.source, relative, description)
                if text is None:
                    return None
                payload = text.encode("utf-8")
            else:
                return None
        if len(payload) > MAX_ARTIFACT_BYTES:
            raise ValueError("MD artifact exceeds the inline size limit; use Results to inspect it.")
        if expected and hashlib.sha256(payload).hexdigest() != expected:
            raise ValueError("MD artifact digest mismatch; changed analysis is not displayed.")
        return payload.decode("utf-8")

    @staticmethod
    def _bounded_bytes(path, root=None):
        if root is not None and not path.resolve().is_relative_to(root.resolve()):
            raise ValueError("MD artifact escapes its cache.")
        with path.open("rb") as handle:
            payload = handle.read(MAX_ARTIFACT_BYTES + 1)
        if len(payload) > MAX_ARTIFACT_BYTES:
            raise ValueError("MD artifact exceeds the inline size limit.")
        return payload

    def object(self, relative, description="MD report", *, expected=""):
        text = self.text(relative, description, expected=expected)
        if text is None:
            return None
        payload = json.loads(text)
        if not isinstance(payload, dict):
            raise ValueError("MD report must be a JSON object.")
        return payload

    def report(self):
        context = self.observation["context"]
        return self.object(context["source_summary_relative_path"], "stored MD analysis",
            expected=context.get("source_summary_sha256") or "")


def matching_study(reader, outcome):
    context = reader.observation["context"]
    identifier = context.get("study_id")
    if not identifier:
        return None
    manifest = reader.object("artifacts/md/study.json", "MD study manifest",
        expected=context.get("study_manifest_sha256") or "")
    if manifest is None:
        return None
    if manifest.get("schema_version") != 1 or not isinstance(manifest.get("studies"), list):
        raise ValueError("Unsupported MD study manifest.")
    matches = [row for row in manifest["studies"] if isinstance(row, dict) and row.get("study_id") == identifier]
    if len(matches) != 1:
        raise ValueError("The source manifest does not identify one matching MD study.")
    study = matches[0]
    for field in ("compound_id", "selected_complex_run_id", "source_pose_id", "target_pdb_id", "target_provenance"):
        if context.get(field) and study.get(field) != context[field]:
            raise ValueError("MD manifest disagrees with the SQL study identity: " + field)
    return study


def snapshot_references(study, production_id=""):
    if not study:
        return []
    snapshots = []
    preparations = (study.get("target_provenance") or {}).get("preparations") or []
    digests = {row.get("source_artifact_sha256") for row in preparations if isinstance(row, dict)} - {None, ""}
    selected = str(study.get("selected_complex_run_id") or "")
    if len(digests) == 1 and re.fullmatch(r"[A-Za-z0-9_-]+", selected):
        snapshots.append({"label": "Original selected pose", "path": f"artifacts/legacy/{selected}/selected_complex.pdb",
            "sha256": next(iter(digests)), "note": "Original selected complex, verified against preparation source evidence. Not an MD frame."})
    seen = set()
    for row in study.get("production_runs") or []:
        if not isinstance(row, dict) or production_id and row.get("run_id") != production_id:
            continue
        path, digest = row.get("topology_path"), row.get("topology_sha256")
        if path and digest and (path, digest) not in seen:
            seen.add((path, digest))
            snapshots.append({"label": "Recorded MD topology snapshot", "path": path, "sha256": digest,
                "note": "Topology paired with this trajectory. Not asserted to be its production-start coordinates or RMSD reference; continuations may start from a checkpoint."})
    return snapshots


def ligand_atom_block(reader, study, report, production_id):
    inline = report.get("ligand_pdb")
    if isinstance(inline, str) and inline.strip():
        return inline, "Atom source: ligand block embedded in this verified production report."
    if not study:
        return "", "No matching study manifest is available for an atom-name source."
    payload = report.get("md_result") or report.get("production") or {}
    outputs = payload.get("output_files") if isinstance(payload, dict) else {}
    outputs = outputs if isinstance(outputs, dict) else {}
    prepared = str(outputs.get("system_pdb") or "")
    prepared_match = re.fullmatch(r"/prepared-system/([A-Za-z0-9_-]+)/([A-Za-z0-9_.-]+\.pdb)", prepared)
    candidates = set()
    for row in study.get("production_runs") or []:
        if not isinstance(row, dict) or not row.get("topology_path") or not row.get("topology_sha256"):
            continue
        exact_segment = row.get("run_id") == production_id
        shared_preparation = False
        if prepared_match:
            preparation, filename = prepared_match.groups()
            shared_preparation = row["topology_path"] in {
                f"artifacts/legacy/{preparation}/{preparation}/{filename}",
                f"artifacts/legacy/{preparation}/{filename}",
            }
            if exact_segment and not shared_preparation:
                raise ValueError("The recorded topology pairing contradicts this production's prepared-system identity.")
        if exact_segment or shared_preparation:
            candidates.add((row["topology_path"], row["topology_sha256"]))
    if len(candidates) != 1:
        return "", "No unique topology pairing is recorded for this production's ligand atoms."
    path, digest = next(iter(candidates))
    text = reader.text(path, "MD ligand atom-name topology", expected=digest)
    if text is None:
        return "", "The paired topology atom-name source is loading or unavailable. Use Refresh to retry."
    selected = report.get("selected_ligand")
    selected = selected if isinstance(selected, dict) else {}
    resname = str(selected.get("resname") or "").strip()
    if not resname:
        return "", "The production does not identify a selected ligand residue name."
    groups = {}
    for line in text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() == resname:
            identity = (line[21:22].strip(), line[22:26].strip(), line[26:27].strip())
            groups.setdefault(identity, []).append(line)
    if selected.get("atom_count") and selected.get("chain"):
        insertion = str(selected.get("icode") or "").strip().replace("_", "")
        identity = (str(selected["chain"]), str(selected.get("resseq") or ""), insertion)
        groups = {key: value for key, value in groups.items() if key == identity}
    if len(groups) != 1:
        return "", "The paired topology does not identify one unambiguous selected ligand residue."
    records = next(iter(groups.values()))
    if len(records) > 1200:
        return "", "The selected topology residue exceeds the atom-diagram size limit."
    return "\n".join([*records, "END", ""]), "Atom source: digest-verified MD topology paired with this production / exact prepared-system identity."
