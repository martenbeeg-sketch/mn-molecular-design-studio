"""Assemble editable pipeline drafts from the Capabilities table."""

from mn_molecular_design_studio.pipeline_blocks import add_block, blank_program, connect_blocks
from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, connection_compatibility, new_graph, validate_graph


MD_ENGINES = ("tile_md_openmm", "tile_md_gromacs")
POSE_ANALYSES = ("native_ligand_interactions", "plip", "pandamap", "posebusters")
DOCKING = tuple(kind for kind in NODE_CATALOG if kind.startswith("tile_docking_"))
SAVED_DESIGN_REFOLDERS = {"af2", "boltz2", "rf3", "prediction_openfold3_ob0"}
DATABASE_SEQUENCE_TOOLS = {"sequence_qc", "sequence_clustering", "proteinsol", "esm1v_score"}

# Task relevance is separate from port compatibility. A compatible wire alone
# cannot tell whether a source has the right molecules or measurements.
SAVED_SOURCES = {
    "project_design_source": "Completed design job",
    "run_source": "Saved Results from a run",
}
_SOURCE_REQUIRED = {"candidate_rank", "tile_validation_acceptance"}
_TARGET_ANALYSES = {"surface_analysis", *(kind for kind in NODE_CATALOG if kind.startswith("tile_detection_"))}


def _port(source, target):
    return next((port["id"] for port in NODE_CATALOG[target]["inputs"]
                 if connection_compatibility(source, target, target_handle=port["id"])["compatible"]), None)


def source_choices(kind, *, database_available=False):
    """Existing input routes offered for the selected operation."""
    if kind == "candidate_rank":
        # Ranking needs a materialized collection, which Candidate collection provides.
        return {"project_design_source": SAVED_SOURCES["project_design_source"]}
    if kind in POSE_ANALYSES:
        return {"run_source": "Saved docking or cofolding poses"}
    if kind == "abmpnn":
        return {"project_design_source": SAVED_SOURCES["project_design_source"]}
    spec = NODE_CATALOG[kind]
    if kind in SAVED_DESIGN_REFOLDERS | {"tile_tmol"} or spec["category"] in {"Analysis", "Collection", "Antibody analysis"}:
        sources = {source: label for source, label in SAVED_SOURCES.items() if _port(source, kind)}
        if database_available and kind in DATABASE_SEQUENCE_TOOLS and _port("molecule_source", kind):
            sources["molecule_source"] = "Saved database entries"
        return sources
    return {}


def build_capability_graph(kind, *, prediction_mode=None, ligand_design_mode="", predict_affinity=False,
                           bioemu2_mode=None, source="", prepare=False, design_task=""):
    """Build the required steps, leaving actual input selection to the node editors."""
    title = NODE_CATALOG[kind]["label"]
    if prediction_mode:
        from mn_molecular_design_studio.scientific_tasks import new_task_graph

        config = {"methods": [kind], "mode": prediction_mode}
        if predict_affinity:
            config["method_settings"] = {kind: {"predict_affinity": True}}
            title = f"{title} affinity"
        graph = new_task_graph("molecular_structure_prediction", config)
    elif kind == "md_analysis":
        graph = new_graph("md_analysis")
    elif kind == "abmpnn":
        graph = new_graph("abmpnn_cdr_redesign")
    else:
        graph = blank_program()
        if kind == "tile_bioemu2":
            from mn_molecular_design_studio.engine.workflow_engine.workflows.bioemu2 import BIOEMU2_MODES

            mode = bioemu2_mode or "propagate"
            if mode not in BIOEMU2_MODES or (source and mode != "propagate"):
                raise ValueError("Choose starting inputs compatible with this BioEmu-2 mode.")
            title = f"{title} · {BIOEMU2_MODES[mode][0]}"

        if kind == "tile_bioemu2":
            allowed = SAVED_SOURCES if mode == "propagate" else {}
        else:
            allowed = source_choices(kind, database_available=True)
        if source and source not in allowed:
            raise ValueError("Choose a supported saved input for this step.")

        if not source and (kind in _SOURCE_REQUIRED or kind == "af2"):
            source = "project_design_source"
        if not source and kind in _TARGET_ANALYSES | set(DOCKING) | {"tile_affinity_nesso"}:
            source = "target"
        if prepare and kind in MD_ENGINES:
            source = "tile_complex_prepare"
        if prepare and kind == "tile_library_centroids":
            source = "tile_prepare_library"

        reference_ligand = False
        if ligand_design_mode:
            from mn_molecular_design_studio.pipeline_tile_workflows import LIGAND_DESIGN_MODES, TILES, ligand_mode_fields

            fields = ligand_mode_fields(TILES[kind], ligand_design_mode)
            reference_ligand = "reference_ligand_sdf" in fields
            if ligand_design_mode in {"pocket", "reference_ligand"}:
                source = "target"
            title = f"{title} · {LIGAND_DESIGN_MODES[ligand_design_mode]}"

        parent = None
        if source:
            graph, parent = add_block(graph, source)
        if kind == "candidate_rank":
            graph, parent = add_block(graph, "candidate_merge", source=parent)
            source = "candidate_merge"
        port = _port(source, kind) if parent else None
        if parent and port is None:
            raise ValueError("The selected input cannot connect to this tool.")
        graph, selected = add_block(graph, kind, source=parent, target_port=port or "in")
        params = next(node["data"]["params"] for node in graph["nodes"] if node["id"] == selected)
        if ligand_design_mode:
            params["ligand_design_mode"] = ligand_design_mode
        if kind == "complex_analysis":
            params["entry_mode"] = "manual_selection"
        if kind == "tile_bioemu2":
            params["workflow_parameters"]["mode"] = mode
        if design_task:
            _set_design_task(kind, params, design_task)
        if reference_ligand:
            graph, reference = add_block(graph, "ligand_sdf_input")
            graph = connect_blocks(graph, reference, "out", selected, "reference_ligand")
        graph, _ = add_block(graph, "results", source=selected)

    graph["name"] = title
    return validate_graph(graph)["graph"]


def _set_design_task(kind, params, task):
    """Keep canvas-only capability rows on their registered mode."""
    if kind == "rfd4":
        params["task"] = {
            "protein_small_molecule": "small_molecule_binder", "enzyme": "enzyme", "motif": "motif_scaffolding",
            "dna_binder": "dna_binder", "redesign": "existing_structure", "rna_design": "rna_design",
            "new_fold": "unconditional",
        }[task]
    elif kind == "tile_design_rfd3_tasks":
        params["workflow_parameters"]["task"] = {
            "protein_small_molecule": "ligand_binder", "enzyme": "enzyme",
            "dna_binder": "nucleic_acid_binder", "rna_binder": "nucleic_acid_binder",
        }[task]
    elif kind == "tile_design_boltzgen_tasks":
        params["workflow_parameters"]["task"] = "small_molecule" if task == "protein_small_molecule" else "redesign"
    elif kind == "complexa_generator":
        params["pipeline_task"] = "ligand_binder" if task == "protein_small_molecule" else "motif_ligand"
    elif kind in {"rfdiffusion", "genie3", "protpardelle"}:
        params["task"] = {
            "new_fold": "custom" if kind == "rfdiffusion" else "unconditional",
            "motif": "motif",
            "redesign": "redesign" if kind == "rfdiffusion" else "partial",
        }[task]


def _result_edges(graph):
    results = {node["id"] for node in graph["nodes"] if node["data"]["kind"] == "results"}
    return [edge for edge in graph["edges"] if edge["target"] in results]


def followup_choices(graph, methods):
    """Give only methods that can consume the selected graph's final output."""
    edges = _result_edges(graph)
    if not edges:
        return ()
    nodes = {node["id"]: node["data"]["kind"] for node in graph["nodes"]}
    source = nodes[edges[0]["source"]]
    return tuple(kind for kind in methods if kind != source and _port(source, kind) and _port(kind, "results"))


def add_followup(graph, kind, *, handoff=False, target_id=None):
    """Connect an optional step while keeping primary outputs available in Results."""
    if kind not in followup_choices(graph, (kind,)):
        raise ValueError("This optional step cannot consume the current output.")
    edge = _result_edges(graph)[0]
    source = next(node["data"]["kind"] for node in graph["nodes"] if node["id"] == edge["source"])
    graph, added = add_block(graph, kind, source=edge["source"], target_port=_port(source, kind))
    if kind == "ligandmpnn":
        params = next(node["data"]["params"] for node in graph["nodes"] if node["id"] == added)
        params.update(sequence_stage="initial", fix_interface_residues=False, backend="shared_ligandmpnn")
    if target_id:
        graph = connect_blocks(graph, target_id, "out", added, "target")
    if handoff:
        graph = {**graph, "edges": [item for item in graph["edges"] if item["id"] != edge["id"]]}
    graph = connect_blocks(graph, added, "out", edge["target"], "in")
    nodes = {node["id"]: node for node in graph["nodes"]}
    result = nodes[edge["target"]]
    result["position"]["x"] = max(result["position"]["x"], nodes[added]["position"]["x"] + 310)
    return validate_graph(graph)["graph"]
