"""Audited design capabilities shared by discovery, the binder picker and MCP."""
from __future__ import annotations

from dataclasses import asdict, dataclass

from mn_molecular_design_studio.binder_pipeline import FORMATS, ROUTES, TARGET_INPUTS


ACCESS = {"pipeline": "Build pipeline", "canvas": "Pipeline Builder", "tile": "Open tool", "api": "API only"}
NATIVE_VERIFICATION_LABELS = {
    "archived_downstream": "Archived downstream verified",
    "archived_results": "Archived Results verified",
    "unverified": "Native execution unverified",
}

# Evidence for exact route and context pairs from the generator acceptance
# playbook. Archived evidence stays scoped to the recorded run configuration.
_ARCHIVED_DOWNSTREAM = {
    "structure": frozenset({
        "rfd4-external", "rfd4-internal", "rfd3", "rfdiffusion", "genie3", "pxdesign",
        "protpardelle-binder", "bindcraft2-generator", "boltzgen-generator", "complexa-generator",
        "esmfold2-protein", "rfdiffusion-beta-pairing", "rfantibody-vhh", "rfantibody-mimosa-vhh",
        "boltzgen-vhh", "boltzgen-vhh-native", "boltzgen-mimosa-vhh",
    }),
    "idr": frozenset({"rfd4-idr", "bindcraft2-idr-miniprotein"}),
}
_ARCHIVED_RESULTS = {
    "structure": frozenset({"rfpeptides", "complexa", "boltzgen-miniprotein", "rfd4-generator", "bindcraft2"}),
    "idr": frozenset({"peptide-trap-idr", "kejia-idr"}),
}
_PRIORITY_REVIEW = {
    "structure": frozenset({"rfd4-generator", "rfd4-external", "rfd4-internal", "bindcraft2-generator", "bindcraft2"}),
    "idr": frozenset({"rfd4-idr", "bindcraft2-idr-miniprotein"}),
}


def native_verification(route_id: str, target_input: str) -> dict:
    """Report native evidence and review priority independently."""
    if route_id in _ARCHIVED_DOWNSTREAM.get(target_input, ()):
        status = "archived_downstream"
    elif route_id in _ARCHIVED_RESULTS.get(target_input, ()):
        status = "archived_results"
    else:
        status = "unverified"
    return {
        "status": status,
        "label": NATIVE_VERIFICATION_LABELS[status],
        "priority_review": route_id in _PRIORITY_REVIEW.get(target_input, ()),
    }


OTHER_TASKS = {
    "protein_small_molecule": "Protein binder for a small molecule",
    "enzyme": "Enzyme or catalytic scaffold",
    "motif": "Motif scaffold",
    "dna_binder": "DNA-target binder",
    "rna_binder": "RNA-target binder",
    "rna_design": "RNA design",
    "ligand_design": "Small-molecule generation",
    "new_fold": "New protein fold",
    "redesign": "Existing structure redesign",
}


@dataclass(frozen=True)
class Capability:
    engine: str
    output: str
    workflow: str
    access: str
    target_input: str = "structure"
    note: str = ""
    route_id: str = ""
    canvas_node: str = ""
    reference: str = ""


def _capabilities():
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
    from mn_molecular_design_studio.pipeline_tile_workflows import TILES

    capabilities = [
        Capability(route.engine, route.format, NODE_CATALOG[route.method]["workflow"], "pipeline",
                   target_input, route.summary, route.id, route.method)
        for route in ROUTES for target_input in route.target_inputs
    ]
    for engine, workflow, node, tasks in (
        ("RFD4-Proteina", "design.rfd4_proteina", "rfd4", ("protein_small_molecule", "enzyme", "motif", "dna_binder", "redesign", "rna_design", "new_fold")),
        ("RFdiffusion3", "design.rfdiffusion3", "tile_design_rfd3_tasks", ("protein_small_molecule", "enzyme", "dna_binder", "rna_binder")),
        ("RFdiffusion", "design.rfdiffusion_generate", "rfdiffusion", ("motif", "redesign", "new_fold")),
        ("Proteina-Complexa", "design.proteina_complexa_generate", "complexa_generator", ("protein_small_molecule", "motif")),
        ("BoltzGen", "design.boltzgen", "tile_design_boltzgen_tasks", ("protein_small_molecule", "redesign")),
        ("Genie3", "design.genie3_generate", "genie3", ("motif", "new_fold")),
        ("Protpardelle-1c", "design.protpardelle_generate", "protpardelle", ("motif", "redesign", "new_fold")),
        ("AbMPNN", "sequence.abmpnn", "abmpnn", ("redesign",)),
    ):
        capabilities.extend(Capability(engine, task, workflow, "canvas", target_input="", canvas_node=node)
                            for task in tasks)
    capabilities.extend(Capability(tile.label, "ligand_design", tile.workflow, "canvas", target_input="", canvas_node=kind)
                        for kind, tile in TILES.items() if tile.category == "Ligand generation")
    return tuple(capabilities)


def picker_tools(target_input, binder_format):
    """List evidenced alternatives, including engines that need their own tool setup."""
    formats = {binder_format, "cyclotide"} if binder_format == "cyclic_peptide" else {binder_format}
    return tuple(row for row in _capabilities() if row.target_input == target_input and row.output in formats
                 and row.access != "pipeline")


def small_molecule_catalogue():
    from mn_molecular_design_studio.pipeline_tile_workflows import LIGAND_DESIGN_MODES, TILES

    return {"modes": LIGAND_DESIGN_MODES.copy(), "items": [
        {"engine": tile.label, "workflow": tile.workflow, "operation": kind, "modes": list(tile.design_modes)}
        for kind, tile in TILES.items() if tile.design_modes
    ]}


def catalogue(*, target_input=None, binder_format=None):
    if target_input is not None and target_input not in TARGET_INPUTS:
        raise ValueError("Unknown target input.")
    if binder_format is not None and binder_format not in FORMATS:
        raise ValueError("Unknown binder format.")
    binder_outputs = ({binder_format, "cyclotide"} if binder_format == "cyclic_peptide" else {binder_format})
    rows = [row for row in _capabilities()
            if (target_input is None or row.target_input == target_input)
            and (binder_format is None or row.output in binder_outputs)]
    items = []
    for row in rows:
        item = asdict(row)
        if row.route_id:
            item["native_verification"] = native_verification(row.route_id, row.target_input)
        items.append(item)
    data = {"schema": "studio.design-capabilities.v1", "target_inputs": TARGET_INPUTS.copy(),
            "binder_formats": FORMATS.copy(), "other_tasks": OTHER_TASKS.copy(), "access": ACCESS.copy(),
            "items": items}
    if target_input is None and binder_format is None:
        data["small_molecule_design"] = small_molecule_catalogue()
    return data
