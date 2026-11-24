"""Read-only discovery of registered tools and supported artifact handoffs."""
from __future__ import annotations

from functools import lru_cache
from typing import Literal

from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS, workflow_catalog
from mn_molecular_design_studio.pipeline_graph import (
    NODE_CATALOG, SCREENING_FEEDBACK, capability_contract, connection_restriction, matching_artifact_type,
)
from mn_molecular_design_studio.tool_catalog import AF2_TILE_WORKFLOW, TOOL_CATALOG_GROUPS, _workflow_discovery_tags, catalog_rows


ConnectionView = Literal["tiles", "workflows", "operations", "connections", "handoffs"]
ConnectionDirection = Literal["both", "incoming", "outgoing"]

# These runners pass files or a source job reference outside the canvas.
_STANDALONE_HANDOFFS = (
    ("library.prepare", "library.centroids", "fingerprint_manifest", "parameters.feature_manifest", ("library_id",)),
    ("library.centroids", "library.seed_dock", "centroid_manifest", "parameters.centroid_manifest",
     ("library_id", "target_pdb", "center_x", "center_y", "center_z")),
    *(("ligand.docking.vina", consumer, "docking_pose_and_receptor", "source_job", ())
      for consumer in ("analysis.native_ligand_interactions", "analysis.plip", "analysis.pandamap", "analysis.posebusters")),
    ("refolding.optimized_prediction", "analysis.ipsae", "predicted_complex_with_pae", "source_job",
     ("binder_chains", "target_chains")),
    *(("refolding.optimized_prediction", consumer, "predicted_complex", "source_job", ("binder_chains", "target_chains"))
      for consumer in ("analysis.biotite_interface", "analysis.rosetta")),
)


def _tile_modes(tile):
    group = TOOL_CATALOG_GROUPS.get(tile["group_id"], {})
    profiles = group.get("engine_profiles", {})
    if group:
        modes = [(label, workflow, profiles.get(label)) for label, workflow in tile["modes"].items()]
    else:
        modes = [(tile["title"], row["workflow"], tile.get("optimized_engine")) for row in tile["workflow_rows"]]
    return [{
        "label": label,
        "workflow": AF2_TILE_WORKFLOW if workflow == "refolding.af2_initial_guess" else workflow,
        "engine_profile": profile,
        "pipeline_operations": [kind for kind, node in sorted(NODE_CATALOG.items())
                                if (node.get("engine_profile") == profile if profile else
                                    node.get("workflow") == workflow or workflow in node.get("workflow_aliases", ()))],
    } for label, workflow, profile in modes]


def _tiles():
    for tile in catalog_rows(workflow_catalog()):
        modes = _tile_modes(tile)
        yield {"id": tile["key"], "title": tile["title"], "description": tile["description"],
               "tags": tile["tags"], "modes": modes,
               "pipeline_operations": sorted({kind for mode in modes for kind in mode["pipeline_operations"]})}


def _operations():
    for kind, spec in sorted(NODE_CATALOG.items()):
        contract = capability_contract(kind, spec)
        yield {"id": kind, "title": spec["label"], "category": spec["category"], "workflow": spec.get("workflow"),
               "engine_profile": spec.get("engine_profile"), "inputs": contract["inputs"], "outputs": contract["outputs"],
               "constraints": contract["constraints"],
               "addable": spec.get("screening_stage", False) or not spec.get("hidden", False)}


def _connections(operation, direction):
    operations = list(_operations())
    for source in operations:
        for target in operations:
            if operation and not (direction != "incoming" and source["id"] == operation
                                  or direction != "outgoing" and target["id"] == operation):
                continue
            for output in source["outputs"]:
                for input_port in target["inputs"]:
                    artifact_type = matching_artifact_type(output, input_port)
                    if artifact_type and not connection_restriction(
                            source["id"], target["id"], artifact_type,
                            source_handle=output["id"], target_handle=input_port["id"]):
                        yield {"source": source["id"], "output": output["id"], "target": target["id"],
                               "input": input_port["id"], "artifact_type": artifact_type}


def _handoffs():
    for source, target, artifact_type, transfer, additional_inputs in _STANDALONE_HANDOFFS:
        yield {"source": source, "target": target, "artifact_type": artifact_type,
               "transfer": transfer, "additional_inputs": list(additional_inputs)}
    docking_workflows = sorted(workflow for workflow in WORKFLOWS if workflow.startswith("ligand.docking."))
    for source in sorted(workflow for workflow in WORKFLOWS if workflow.startswith("ligand.generation.")):
        for target in [*docking_workflows, "ligand.affinity.nesso"]:
            required = ["target_pdb", "center_x", "center_y", "center_z"] if target in docking_workflows else ["target_pdb"]
            yield {"source": source, "target": target, "artifact_type": "generated_compound_set",
                   "transfer": "source_job", "additional_inputs": required}


@lru_cache(maxsize=1)
def editor_tool_catalog():
    """Attach the discovery tags used by All tools to registered canvas steps."""
    workflows = {row["workflow"]: row for row in workflow_catalog()}
    result = {}
    for kind, spec in NODE_CATALOG.items():
        workflow = next((workflows[name] for name in (spec.get("workflow"), *spec.get("workflow_aliases", ()))
                         if name in workflows), None)
        discovery = _workflow_discovery_tags(workflow) if workflow else {}
        if spec.get("modalities"):
            discovery["modalities"] = spec["modalities"]
        tasks = spec.get("molecule_tasks", ())
        if tasks:
            discovery["modalities"] = ["Protein", *[
                modality for task, modality in (("protein_ligand", "Ligand"), ("protein_dna", "DNA"),
                                                ("protein_rna", "RNA"), ("protein_macrocyclic", "Peptide"))
                if task in tasks
            ]]
        description = spec["description"]
        if kind.startswith("tile_") and workflow and workflow.get("description"):
            description = workflow["description"]
        result[kind] = {**spec, "palette": {"description": description, **discovery}}
        if spec.get("screening_stage"):
            result[kind].update(hidden=False, category="Library screening")
            result[kind]["palette"]["functions"] = ["Active learning"]
    return result


@lru_cache(maxsize=1)
def editor_connections():
    """Compact lookup for the canvas, built from the same rows as API/MCP discovery."""
    wires = {}
    for row in _connections("", "both"):
        outputs = wires.setdefault(row["source"], {})
        targets = outputs.setdefault(row["output"], {})
        targets.setdefault(row["target"], {})[row["input"]] = row["artifact_type"]
    handoffs = {}
    for row in _handoffs():
        handoffs.setdefault(row["source"], []).append(row)
    return {"wires": wires, "handoffs": handoffs, "screening_feedback": SCREENING_FEEDBACK}


def connection_catalog(*, view: ConnectionView, query: str, operation: str, direction: ConnectionDirection,
                       artifact_type: str, offset: int, limit: int, port: str = "") -> dict:
    """Project live registry entries into a bounded, searchable result page."""
    if view == "connections":
        rows = _connections(operation, direction)
        if port:
            rows = (row for row in rows if (
                (direction != "incoming" and (not operation or row["source"] == operation) and row["output"] == port)
                or (direction != "outgoing" and (not operation or row["target"] == operation) and row["input"] == port)
            ))
        scope = "Declared canvas port matches. Validate a saved graph to check inputs, formats and scientific constraints."
    elif view == "handoffs":
        rows = _handoffs()
        scope = "Documented standalone handoffs. The referenced job or manifest must be accessible and provide the required artifacts."
    elif view == "operations":
        rows = _operations()
        scope = "Canvas operations with registered inputs, outputs and constraints."
    elif view == "workflows":
        rows = sorted(workflow_catalog(), key=lambda row: row["workflow"])
        scope = "API workflows with their parameter types and required fields."
    else:
        rows = _tiles()
        scope = "All tools tiles and their executable workflow modes."
    query = query.strip().casefold()
    matches = [row for row in rows
               if (view not in {"connections", "handoffs"} or not artifact_type or row["artifact_type"] == artifact_type)
               and (not query or query in str(row).casefold())]
    return {"schema_version": 1, "view": view, "scope": scope, "items": matches[offset:offset + limit],
            "total": len(matches), "offset": offset, "limit": limit}
