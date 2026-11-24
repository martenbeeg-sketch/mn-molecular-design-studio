"""Task-first entry points for existing binder-design workflows."""
from __future__ import annotations

from dataclasses import asdict, dataclass
import copy
import json
import uuid


TARGET_INPUTS = {
    "structure": "Structure",
    "idr": "Sequence · disordered region (IDR)",
    "sequence": "Sequence · folded or unknown",
}
FORMATS = {
    "miniprotein": "Miniprotein", "peptide": "Linear peptide",
    "cyclic_peptide": "Macrocyclic peptide", "vhh": "VHH (nanobody)",
    "cyclotide": "Cyclotide", "scfv": "scFv / paired variable domains", "fab": "Fab",
    "arp": "Ankyrin repeat protein (ARP)", "large_protein": "Large protein",
    "multidomain": "Multidomain binder", "homo_oligomer": "Homo-oligomer",
    "antibody": "Antibody (custom framework)", "tcr": "T-cell receptor (TCR)",
}

BC2_FORMATS = {
    "miniprotein": "binder", "large_protein": "large_binder", "peptide": "peptide",
    "cyclic_peptide": "cyclic_peptide", "vhh": "VHH", "scfv": "scFv", "fab": "Fab", "arp": "ARP",
    "multidomain": "multidomain", "homo_oligomer": "homo_oligomer",
}


@dataclass(frozen=True)
class Route:
    id: str
    engine: str
    variant: str
    target_input: str
    format: str
    method: str
    summary: str


ROUTES = (
    Route("rfd4-internal", "RFD4-Proteina", "Internal OpenFold3 refold", "structure", "miniprotein", "rfd4",
          "Co-design sequence and structure, refold internally, then apply Studio acceptance."),
    Route("rfd4-external", "RFD4-Proteina", "External Boltz-2 refold", "structure", "miniprotein", "rfd4",
          "Co-design sequence and structure, add a visible Boltz-2 step, then apply Studio acceptance."),
    Route("bindcraft2", "BindCraft2", "Native validation and filters", "structure", "miniprotein", "bindcraft2",
          "Run the complete BindCraft2 workflow, including its native admission filters."),
    Route("rfd3", "RFdiffusion3", "Separate sequence design and refold", "structure", "miniprotein", "rfd3",
          "Generate backbones, add ProteinMPNN and Boltz-2, then apply Studio acceptance."),
    Route("bindcraft2-peptide", "BindCraft2", "Native peptide workflow", "structure", "peptide", "bindcraft2",
          "Generate linear peptides with native refolding and peptide-specific filters."),
    Route("bindcraft2-cyclic", "BindCraft2", "Head-to-tail cyclization", "structure", "cyclic_peptide", "bindcraft2",
          "Use the native head-to-tail cyclic peptide workflow, including closure-aware validation."),
    Route("bindcraft2-vhh", "BindCraft2", "VHH framework", "structure", "vhh", "bindcraft2_generator",
          "Generate on a VHH framework, annotate the antibody, refold with Boltz-2 and apply Studio acceptance."),
    Route("boltzgen-vhh", "BoltzGen", "Nanobody mode", "structure", "vhh", "boltzgen_generator",
          "Generate in nanobody mode, annotate the antibody, refold with Boltz-2 and apply Studio acceptance."),
    Route("esmfold2-vhh", "ESMFold2 Binder", "VHH framework", "structure", "vhh", "esm_binder",
          "Generate on a VHH framework, annotate the antibody, refold with Boltz-2 and apply Studio acceptance."),
    Route("rfd4-idr", "RFD4-Proteina", "Internal OpenFold3 refold", "idr", "miniprotein", "rfd4",
          "Experimental sequence-only route. Refold internally and apply the IDR acceptance profile."),
    *(Route(key, engine, "Sequence design and Boltz-2", "structure", "miniprotein", key,
            "Generate backbones, design sequences, refold with Boltz-2 and apply Studio acceptance.")
      for key, engine in (("rfdiffusion", "RFdiffusion"), ("genie3", "Genie3"), ("pxdesign", "PXDesign"))),
    Route("complexa", "Proteina-Complexa", "Native validation and filters", "structure", "miniprotein",
          "complexa", "Run the complete Proteina-Complexa protein binder workflow."),
    Route("esmfold2-protein", "ESMFold2 Binder", "External Boltz-2 refold", "structure", "miniprotein",
          "esm_binder", "Generate sequences, refold with Boltz-2 and apply Studio acceptance."),
    *(Route("bindcraft2-vhh-native" if fmt == "vhh" else f"bindcraft2-{fmt.replace('_', '-')}",
            "BindCraft2", "Native validation and filters", "structure", fmt,
            "bindcraft2", "Use the native format preset and its built-in validation and filters.")
      for fmt in ("large_protein", "vhh", "scfv", "fab", "arp", "multidomain", "homo_oligomer")),
    *(Route("boltzgen-vhh-native" if fmt == "vhh" else f"boltzgen-{fmt}",
            "BoltzGen", "Native validation and ranking", "structure", fmt, "boltzgen",
            f"Run the full BoltzGen workflow using its {FORMATS[fmt].lower()} protocol.")
      for fmt in ("miniprotein", "peptide", "cyclotide", "vhh")),
    Route("esmfold2-scfv", "ESMFold2 Binder", "External Boltz-2 refold", "structure", "scfv", "esm_binder",
          "Use the paired variable-domain scaffold mode, then refold and assess the resulting complex."),
    Route("rfd4-vhh", "RFD4-Proteina", "Internal OpenFold3 refold", "structure", "vhh", "rfd4",
          "Use RFD4's built-in nanobody framework, refold internally and apply Studio acceptance."),
    Route("rfpeptides", "RFpeptides", "Cyclic peptide backbone generation", "structure", "cyclic_peptide", "tile_design_rfpeptides",
          "Generate cyclic peptide backbones. Sequence design and refolding are separate steps."),
    Route("rfantibody-vhh", "RFantibody", "Native antibody workflow", "structure", "vhh", "tile_design_rfantibody",
          "Use the VHH framework with native sequence design and structure evaluation."),
    Route("rfantibody-scfv", "RFantibody", "Native paired-domain workflow", "structure", "scfv", "tile_design_rfantibody",
          "Use its paired variable-domain framework and native structure evaluation."),
)


def catalogue():
    return {"template_id": "design_binder", "target_inputs": TARGET_INPUTS.copy(), "binder_formats": FORMATS.copy(),
            "routes": [asdict(route) for route in ROUTES]}


def compatible_routes(target_input, binder_format):
    return tuple(route for route in ROUTES if route.target_input == target_input
                 and (route.format == binder_format or (binder_format == "cyclic_peptide" and route.format == "cyclotide")))


def idr_target(params):
    return params.get("target_input") == "idr"


def legacy_sequence_carrier(params):
    from mn_molecular_design_studio.target_identity import is_sequence_carrier_reference

    return not idr_target(params) and is_sequence_carrier_reference(params.get("target_pdb"))


def stage_idr_target(sequence, owner_id, *, hotspots=""):
    from mn_molecular_design_studio.auth import user_input_root
    from mn_molecular_design_studio.rfd4_proteina_input import (
        idr_binder_condition_spec, parse_idr_fasta, store_idr_sequence_carrier,
    )

    sequence = parse_idr_fasta(sequence)
    idr_binder_condition_spec(sequence, hotspots=hotspots)
    path = store_idr_sequence_carrier(sequence, user_input_root(owner_id))
    return {"target_pdb": str(path), "target_chains": "A", "target_input": "idr",
            "target_sequence": sequence, "hotspots": hotspots, "coldspots": "",
            "target_source_molecule_id": "", "target_source_structure_id": ""}


def new_binder_pipeline(route_id, *, target_input, binder_format, inputs=None, graph=None, target_node_id=None):
    """Build a new draft or graft an independent branch onto one existing target."""
    from mn_molecular_design_studio.pipeline_graph import SCHEMA, validate_graph
    from mn_molecular_design_studio.scientific_tasks import _node, new_task_graph

    route = next((item for item in compatible_routes(target_input, binder_format) if item.id == route_id), None)
    if route is None:
        raise ValueError("The selected engine route does not support this target input and binder format.")
    if target_node_id and not graph:
        raise ValueError("The selected Target step requires an existing canvas.")
    if graph and not target_node_id and graph.get("nodes"):
        raise ValueError("Select a Target step to add a binder branch to this canvas.")
    if target_node_id and inputs:
        raise ValueError("Edit sequence or structure inputs on the connected Target step.")
    configuration = {"methods": [route.method]}
    if binder_format == "vhh" and route.method in {"bindcraft2_generator", "boltzgen_generator", "esm_binder"}:
        configuration["binder_format"] = binder_format
    if route.method == "bindcraft2" and binder_format != "miniprotein":
        from mn_molecular_design_studio.engine.workflow_engine.workflows.bindcraft2 import BC2_PRESET_LENGTHS

        modality = BC2_FORMATS[binder_format]
        configuration["method_settings"] = {"bindcraft2": {"modality": modality, "binder_length":
            "-".join(map(str, BC2_PRESET_LENGTHS.get(modality, ())))}}
    if route.method == "boltzgen":
        task = {"miniprotein": "protein", "vhh": "nanobody"}.get(route.format, route.format)
        configuration["method_settings"] = {"boltzgen": {"task": task}}
    if route.id == "esmfold2-scfv":
        configuration["method_settings"] = {"esm_binder": {"binder_type": "scfv"}}
    if route.method == "rfd4":
        configuration["method_settings"] = {"rfd4": {"options_json": json.dumps({"refold": route.id != "rfd4-external"})}}
    if route.method.startswith("tile_") or route.id == "rfd4-vhh":
        recipe = {"schema_version": SCHEMA, "id": uuid.uuid4().hex}
        values = {}
        if route.method == "tile_design_rfantibody":
            values = {"workflow_parameters": {"task": "nanobody", "framework": "h-NbBCII10"}
                      if route.format == "vhh" else {"task": "antibody", "framework": "hu-4D5-8_Fv"}}
        elif route.id == "rfd4-vhh":
            values = {"task": "nanobody", "num_designs": 1, "options_json": '{"refold": true}'}
        steps = [_node("target", "target"), _node(route.method, "design", values)]
    else:
        recipe = new_task_graph("targeted_protein_binder", configuration)
        steps = [node for node in recipe["nodes"] if node["data"]["kind"] not in {"binder_context", "results"}]
    if route.method == "rfd4" and route.id != "rfd4-external":
        steps = [node for node in steps if node["data"]["kind"] != "boltz2"]
    target = steps[0]
    target["data"]["params"].update(inputs or {})
    target["data"]["params"]["target_input"] = target_input
    target["data"]["label"] = "Target sequence · IDR" if target_input == "idr" else "Target structure"
    engine = steps[1]["data"]
    engine["label"] = route.engine + " · " + (route.variant if route.method == "rfd4" else FORMATS[route.format])
    if route.method not in {"bindcraft2", "boltzgen", "complexa"} and not route.method.startswith("tile_"):
        steps.extend([
            _node("candidate_merge", "candidates"),
            _node("shared_filter", "acceptance", {"binder_format": binder_format,
                "target_type": "fasta_idr" if target_input == "idr" else "structured"}, "Validation & acceptance"),
        ])
    steps.append(_node("results", "results"))
    result = copy.deepcopy(graph) if graph else {**recipe, "nodes": [], "edges": []}
    if not result["nodes"]:
        result["name"] = f"{FORMATS[route.format]} design · {route.engine}"
    branch_prefix = uuid.uuid4().hex[:12]
    for node in steps:
        node["id"] = f"{branch_prefix}-{node['id']}"
    x, y = 60, 85
    if target_node_id:
        existing = next((node for node in result["nodes"] if node["id"] == target_node_id
                         and node["data"]["kind"] == "target" and node["data"].get("enabled") is not False), None)
        if existing is None:
            raise ValueError("The selected Target step is unavailable.")
        if legacy_sequence_carrier(existing["data"]["params"]):
            raise ValueError("Choose a declared sequence target for the sequence-only binder route.")
        if idr_target(existing["data"]["params"]) != (target_input == "idr"):
            raise ValueError("Choose a route compatible with the selected Target step.")
        x, y = existing["position"]["x"], max(node["position"]["y"] for node in result["nodes"]) + 240
        steps[0] = existing
    for index, node in enumerate(steps):
        if index or not target_node_id:
            node["position"] = {"x": x + index * 315, "y": y}
            result["nodes"].append(node)
    for source, destination in zip(steps, steps[1:]):
        result["edges"].append({"id": f"edge-{destination['id']}", "source": source["id"], "sourceHandle": "out",
                                "target": destination["id"],
                                "targetHandle": "target" if destination["data"]["kind"].startswith("tile_design_") else "in"})
    return validate_graph(result, allow_incomplete=bool(target_node_id))["graph"]
