"""Audited design capabilities shared by discovery, the binder picker and MCP."""
from __future__ import annotations

from dataclasses import asdict, dataclass

from mn_molecular_design_studio.binder_pipeline import BC2_FORMATS, FORMATS, ROUTES, TARGET_INPUTS


ACCESS = {"pipeline": "Build pipeline", "tile": "Open tool", "api": "API only", "upstream": "Upstream only"}
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
                   route.target_input, route.summary, route.id, route.method)
        for route in ROUTES
    ]
    capabilities.extend([
        *(Capability("BindCraft2", fmt, "design.bindcraft2", "tile", "idr",
                     "Use the Sequence (disordered / IDR) target input in BindCraft2.")
          for fmt in BC2_FORMATS),
        Capability("BindCraft", "miniprotein", "design.bindcraft", "api",
                   note="Registered classic design workflow; separate from BindCraft2."),
        Capability("Protpardelle-1c", "miniprotein", "design.protpardelle_generate", "tile",
                   note="The binder mode uses structure conditioning; configure the input on its tool."),
        Capability("RFdiffusion β-pairing", "miniprotein", "design.rfdiffusion_beta_pairing", "tile",
                   note="Requires selected matching target β-strands.", canvas_node="rfdiffusion_beta_pairing"),
        Capability("Peptide Trap", "miniprotein", "design.peptide_trap", "tile", "idr",
                   "Requires a compatible template cleft and target register; binder format follows its scaffold.", canvas_node="peptide_trap"),
        Capability("Kejia IDR threading", "miniprotein", "design.kejia_idr", "tile", "idr",
                   "Requires a compatible binder–peptide template; binder format follows the template.", canvas_node="kejia_idr"),
        Capability("RFD4-Proteina", "antibody", "design.rfd4_proteina", "tile",
                   note="Requires an antibody–target complex, framework chains and selected CDRs.", canvas_node="rfd4"),
        Capability("RFD4-Proteina", "tcr", "design.rfd4_proteina", "tile",
                   note="Requires a supplied target–TCR complex and explicit framework selections.", canvas_node="rfd4"),
        Capability("MIMOSA · RFantibody", "vhh", "design.rfantibody_mimosa", "tile",
                   note="Configure the scaffold and loop guidance in the MIMOSA tool.", canvas_node="tile_design_rfantibody_mimosa"),
        Capability("MIMOSA · BoltzGen", "vhh", "design.boltzgen_mimosa", "tile",
                   note="Configure the nanobody scaffold and optional motif guidance in the MIMOSA tool.", canvas_node="tile_design_boltzgen_mimosa"),
        Capability("RFD4-Proteina", "cyclic_peptide", "design.rfd4_proteina", "upstream",
                   note="Upstream cyclic-peptide recipe. The installed Studio adapter has no dedicated cyclic-binder preset.",
                   reference="https://github.com/RosettaCommons/RFD4-Proteina/pull/82"),
        Capability("RFdiffusion3", "cyclic_peptide", "design.rfdiffusion3", "upstream",
                   note="The current upstream input guide lists cyclic chains. Studio's pinned input adapter does not expose that field.",
                   reference="https://github.com/RosettaCommons/foundry/blob/production/models/rfd3/docs/input.md"),
    ])
    for engine, workflow, node, tasks in (
        ("RFD4-Proteina", "design.rfd4_proteina", "rfd4", ("protein_small_molecule", "enzyme", "motif", "dna_binder", "rna_design", "new_fold", "redesign")),
        ("RFdiffusion3", "design.rfdiffusion3", "", ("protein_small_molecule", "enzyme", "dna_binder", "rna_binder")),
        ("RFdiffusion", "design.rfdiffusion_generate", "rfdiffusion", ("motif", "redesign", "new_fold")),
        ("Proteina-Complexa", "design.proteina_complexa_generate", "complexa_generator", ("protein_small_molecule", "motif")),
        ("BoltzGen", "design.boltzgen", "", ("protein_small_molecule", "redesign")),
        ("Genie3", "design.genie3_generate", "genie3", ("motif", "new_fold")),
        ("Protpardelle-1c", "design.protpardelle_generate", "protpardelle", ("motif", "redesign", "new_fold")),
        ("AbMPNN", "sequence.abmpnn", "abmpnn", ("redesign",)),
    ):
        capabilities.extend(Capability(engine, task, workflow, "tile", target_input="", canvas_node=node)
                            for task in tasks)
    capabilities.extend(Capability(tile.label, "ligand_design", tile.workflow, "tile", target_input="", canvas_node=kind)
                        for kind, tile in TILES.items() if tile.category == "Ligand generation")
    return tuple(capabilities)


def picker_tools(target_input, binder_format):
    """List evidenced alternatives, including engines that need their own tool setup."""
    formats = {binder_format, "cyclotide"} if binder_format == "cyclic_peptide" else {binder_format}
    return tuple(row for row in _capabilities() if row.target_input == target_input and row.output in formats
                 and row.access != "pipeline")


def catalogue(*, target_input=None, binder_format=None):
    if target_input is not None and target_input not in TARGET_INPUTS:
        raise ValueError("Unknown target input.")
    if binder_format is not None and binder_format not in FORMATS:
        raise ValueError("Unknown binder format.")
    binder_outputs = ({binder_format, "cyclotide"} if binder_format == "cyclic_peptide" else {binder_format})
    rows = [row for row in _capabilities()
            if (target_input is None or row.target_input == target_input)
            and (binder_format is None or row.output in binder_outputs)]
    return {"schema": "studio.design-capabilities.v1", "target_inputs": TARGET_INPUTS.copy(),
            "binder_formats": FORMATS.copy(), "other_tasks": OTHER_TASKS.copy(), "access": ACCESS.copy(),
            "items": [asdict(row) for row in rows]}
