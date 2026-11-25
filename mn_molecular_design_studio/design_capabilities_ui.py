"""Browse design tasks and registered canvas operations by capability."""
from __future__ import annotations

from html import escape

import streamlit as st

from mn_molecular_design_studio.binder_pipeline import (
    RFD4_MODULAR_ROUTES, RFD4_REFOLD_ENGINES, compatible_routes, new_binder_pipeline, supported_interface_constraints,
)
from mn_molecular_design_studio.capability_pipeline import (
    DOCKING, MD_ENGINES, POSE_ANALYSES, SAVED_DESIGN_REFOLDERS, add_followup, build_capability_graph, followup_choices, source_choices,
)
from mn_molecular_design_studio.design_capabilities import ACCESS, catalogue, small_molecule_catalogue
from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG


_PROTEIN_DESIGN_CATEGORY = "Protein design"
_PREDICTION_CATEGORY = "Structure prediction & cofolding"
_SMALL_MOLECULE_CATEGORY = "Small-molecule design"
_SAMPLING_CATEGORY = "Conformational sampling"
_DOCKING_CATEGORY = "Docking & affinity"
_ANALYSIS_CATEGORY = "Analysis & validation"
_CANDIDATE_SELECTION_TASK = "Combine, filter & rank candidates"
_BOLTZ_AFFINITY_LABEL = "Boltz-2 affinity"
_PREDICTION_ENGINES = frozenset({
    "af2", "af3", "boltz2", "rf3", "rfd4_fold", "prediction_chai1", "prediction_colabfold",
    "prediction_esmfold2", "prediction_openfold3_ob0", "prediction_protenix_v2",
})
_OPERATION_CATEGORIES = {
    _PREDICTION_CATEGORY: {"Structure prediction"},
    _DOCKING_CATEGORY: {"Docking"},
    "Structure preparation": {"Structure preparation"},
    "Molecular dynamics": {"Molecular dynamics"},
    _SAMPLING_CATEGORY: {"Conformational sampling"},
    "Library preparation": {"Ligand library"},
    _ANALYSIS_CATEGORY: {"Analysis", "Antibody analysis", "Collection"},
}
_ANALYSIS_TASKS = {
    "Protein design evaluation": {
        "tile_binderscope": "Assess binder interfaces, liabilities, and tag placement",
        "biotite_interface": "Measure residue contacts across a protein interface",
        "complex_analysis": "Compare saved complexes for downstream analysis",
        "ipsae": "Measure interface confidence from PAE and compare refolded structures",
        "rosetta_metrics": "Measure interface packing, buried surface area, and hydrogen bonds",
        "sabr": "Annotate antibody CDRs",
        "tile_validation_acceptance": "Apply acceptance criteria to saved design evidence",
    },
    "Sequence & protein properties": {
        "sequence_qc": "Summarize length, composition, isoelectric point, and hydropathy",
        "proteinsol": "Estimate protein solubility from sequence",
        "esm1v_score": "Score amino-acid compatibility with sequence context",
        "esmif_score": "Score sequence compatibility with a backbone structure",
        "pep_patch": "Measure surface electrostatic patches",
    },
    "Similarity & clustering": {
        "sequence_clustering": "Group candidates by sequence similarity",
        "interface_clustering": "Group complexes by similar interface contacts",
        "foldseek": "Compare and cluster structures by shape similarity",
    },
    "Protein–ligand pose analysis": {
        "native_ligand_interactions": "Classify hydrogen bonds, hydrophobic contacts, water bridges, and salt bridges in a pose",
        "plip": "Identify protein–ligand contacts",
        "pandamap": "Map contacts between a ligand and binding-site residues",
        "posebusters": "Check ligand geometry and clashes in a pose",
    },
    "Binding-site detection": {
        "surface_analysis": "Compare surface and interface predictions across methods",
        "tile_detection_masif_seed": "Find protein surface patches for interaction design",
        "tile_detection_pesto": "Predict interface residues on a protein structure",
        "tile_detection_scannet": "Predict protein–protein interaction sites",
        "tile_detection_surf2spot": "Predict binding-site hotspots on a protein surface",
    },
    "MD trajectory analysis": {
        "md_analysis": "Measure motion, radius of gyration, and protein–ligand interactions across trajectory frames",
    },
    _CANDIDATE_SELECTION_TASK: {
        "candidate_merge": "Combine design lists for the same target and keep each design's source.",
        "candidate_analyzer": "Put predictions and measurements for the same design variant side by side.",
        "filter": "Use your own cutoffs to route designs to Pass, Fail, or Unresolved.",
        "shared_filter": "Apply a binder acceptance profile to a combined list and record each decision.",
        "candidate_rank": "Shortlist eligible designs within comparable scoring protocols. Optionally keep one per exact sequence.",
    },
}


_matrix = st.components.v2.component(
    "studio_design_capabilities",
    html='<div class="matrix" role="region" aria-label="Capability table" tabindex="0"></div>',
    css="""
        .matrix {
            overflow: auto;
            border: 1px solid var(--st-border-color, #ddd);
            border-radius: .5rem;
        }
        table { border-spacing: 0; width: 100%; font-size: .875rem; }
        th, td {
            border-bottom: 1px solid var(--st-border-color, #ddd);
            padding: .45rem .65rem;
            text-align: left;
        }
        thead th { min-width: 8rem; }
        th:first-child {
            position: sticky;
            left: 0;
            min-width: 9rem;
            background: var(--st-background-color, white);
        }
        tbody tr:last-child > * { border-bottom: 0; }
        .operation-use { display: block; font-weight: normal; }
        button {
            font: inherit;
            color: inherit;
            text-align: left;
            background: transparent;
            border: 1px solid var(--st-border-color, #999);
            border-radius: .4rem;
            padding: .32rem .5rem;
            cursor: pointer;
        }
        button[data-pipeline] { color: #14543e; background: #e7f5ee; border-color: #8bb8a4; }
        button:hover, button:focus-visible { outline: 2px solid var(--st-primary-color, #267356); }
    """,
    js="""
        export default function({ data, parentElement, setTriggerValue }) {
            const matrix = parentElement.querySelector('.matrix');
            if (matrix.markup !== data) {
                matrix.innerHTML = data;
                matrix.markup = data;
            }
            matrix.onclick = (event) => {
                const button = event.target.closest('button');
                if (button) setTriggerValue('choose', [button.dataset.engine, button.dataset.output]);
            };
        }
    """,
)


def _matrix_markup(by_engine, columns, access, *, action_column=False, primary_access="pipeline"):
    labels = ["Studio access"] if action_column else columns.values()
    headers = "".join(f'<th scope="col">{escape(label)}</th>' for label in labels)
    table = [f'<table><thead><tr><th scope="col">Engine / checkpoint</th>{headers}</tr></thead><tbody>']
    for engine, cells in sorted(by_engine.items(), key=lambda entry: entry[0].casefold()):
        table.append(f'<tr><th scope="row">{escape(engine)}</th>')
        for code, label in columns.items():
            level = next((level for level in access if level in cells.get(code, ())), None)
            if level is None:
                table.append("<td>—</td>")
                continue
            title = f"{access[level]} · {engine} · {label}"
            style = ' data-pipeline="true"' if level == primary_access else ""
            table.append(
                f'<td><button type="button" data-engine="{escape(engine)}" data-output="{escape(code)}"'
                f' aria-label="{escape(title)}"{style}>{escape(access[level])}</button></td>'
            )
        table.append("</tr>")
    return "".join([*table, "</tbody></table>"])


def _binder_options(key, routes, *, label="Workflow variant"):
    route = (st.selectbox(label, routes, format_func=lambda item: item.variant, key=key + "-variant")
             if len(routes) > 1 else routes[0])
    st.caption(route.summary)
    options = {}
    supported = supported_interface_constraints(route)
    if supported:
        options["interface_constraints"] = st.checkbox("Add interface constraints", key=key + "-constraints")
        if options["interface_constraints"]:
            labels = {"hotspots": "hotspots (contact)", "coldspots": "coldspots (avoid)"}
            st.caption(f"{route.engine} applies " + " and ".join(labels[kind] for kind in supported) + ".")
            if "coldspots" not in supported:
                st.caption("Coldspots are unavailable for this workflow.")
    elif route.id == "rfd4-idr":
        st.caption("Set sequence hotspot positions in the Target step.")
    if route.id in RFD4_MODULAR_ROUTES:
        options["soluble_mpnn"] = st.checkbox(
            "Add SolubleMPNN sequence redesign", key=key + "-soluble", help="Keep the original designs and protect their interface in redesign variants."
        )
    if route.id == "rfd4-external":
        options["refold_engine"] = st.selectbox(
            "Refolding engine", list(RFD4_REFOLD_ENGINES), format_func=RFD4_REFOLD_ENGINES.get, key=key + "-refold"
        )
    if route.id in {"rfd4-internal", "rfd4-external"}:
        options["interface_analysis"] = st.checkbox("Add ipSAE & Structure Analysis", key=key + "-analysis")
    return {"route_id": route.id, "target_input": route.target_input, "binder_format": route.format, **options}


def render_binder_preset(key, routes, *, button_label="Open pipeline draft"):
    preset = _binder_options(key, routes)
    return preset if st.button(button_label, type="primary", key=key + "-create") else None


def render_design_tasks(key, open_pipeline):
    data = catalogue()
    available_outputs = {row["output"] for row in data["items"]}
    tasks = {"protein_binders": "Protein binders", **{
        code: label for code, label in data["other_tasks"].items()
        if code in available_outputs and code != "ligand_design"
    }}
    task = st.selectbox("Design task", list(tasks), format_func=tasks.get, key=key + "-task")
    other_tasks = task != "protein_binders"
    if other_tasks:
        columns = {task: tasks[task]}
        target = "all"
    else:
        columns = data["binder_formats"]
        targets = {"all": "All target inputs", **data["target_inputs"]}
        target = st.selectbox("Target input", list(targets), index=1, format_func=targets.get, key=key + "-target")
    rows = [row for row in data["items"] if row["output"] in columns
            and (target == "all" or row["target_input"] == target)]
    if not rows:
        st.info("No documented design workflows for this target input.")
        return
    access = data["access"]
    if any(row["access"] == "pipeline" for row in rows):
        st.caption("Choose Build pipeline to configure a workflow.")
    by_engine = {}
    for row in rows:
        cells = by_engine.setdefault(row["engine"], {})
        available = cells.setdefault(row["output"], set())
        available.add(row["access"])
    choice = _matrix(
        data=_matrix_markup(by_engine, columns, access, action_column=other_tasks), key=f"{key}-matrix-{task}-{target}",
        on_choose_change=lambda: None,
    ).get("choose")
    if not choice:
        return
    engine, output = choice
    selected_rows = [row for row in rows if row["engine"] == engine and row["output"] == output]
    if not selected_rows:
        return

    @st.dialog(f"{engine} · {columns[output]}", width="large")
    def details():
        try:
            target_inputs = list(dict.fromkeys(row["target_input"] for row in selected_rows if row["access"] == "pipeline"))
            if target_inputs:
                target_input = (st.selectbox("Target context", target_inputs, format_func=data["target_inputs"].get,
                                             key=f"{key}-{engine}-{output}-context")
                                if len(target_inputs) > 1 else target_inputs[0])
                routes = [route for route in compatible_routes(target_input, output) if route.engine == engine]
                setup_key = f"{key}-{target_input}-{engine}-{output}"
                preset = _binder_options(setup_key, routes, label="Workflow goal")
                graph = new_binder_pipeline(**preset)
                optional = tuple(kind for flag, kind in (("interface_constraints", "binder_context"),
                                                         ("soluble_mpnn", "ligandmpnn"), ("interface_analysis", "ipsae"))
                                 if preset.get(flag))
                kinds = {node["data"]["kind"] for node in graph["nodes"]}
                finish = (
                    "After the run, use the recorded acceptance stage in My results to collect passing designs."
                    if kinds & {"shared_filter", "bindcraft2", "tile_design_bindcraft"}
                    else ""
                )
                _render_plan(graph, setup_key, open_pipeline, optional=optional, finish=finish)
                return

            for row in selected_rows:
                if row["access"] == "canvas":
                    _operation_setup(row["canvas_node"], f"{key}-{engine}-{output}", open_pipeline, design_task=output)
                    return
                st.write(row["note"] or access[row["access"]])
        except ValueError as exc:
            st.error(str(exc))

    details()


def render_capabilities(key, open_pipeline, *, database_user=None, project_id=None):
    from mn_molecular_design_studio.tool_connections import editor_tool_catalog

    st.subheader("Capabilities")
    category = st.selectbox("Capability", [_PROTEIN_DESIGN_CATEGORY, _SMALL_MOLECULE_CATEGORY, *_OPERATION_CATEGORIES],
                            key=key + "-capability")
    if category == _PROTEIN_DESIGN_CATEGORY:
        render_design_tasks(key, open_pipeline)
        return
    if category == _SMALL_MOLECULE_CATEGORY:
        _render_small_molecule_design(key, open_pipeline)
        return
    if category == _SAMPLING_CATEGORY:
        _render_bioemu2_modes(key, open_pipeline)
        return

    catalog = editor_tool_catalog()
    operations = {kind: spec for kind, spec in catalog.items()
                  if spec["category"] in _OPERATION_CATEGORIES[category] and not spec.get("hidden")
                  and (category != _PREDICTION_CATEGORY or kind in _PREDICTION_ENGINES)}
    if category == _DOCKING_CATEGORY:
        operations["boltz2"] = {**catalog["boltz2"], "label": _BOLTZ_AFFINITY_LABEL}
    table_key = f"{key}-operations-{category}"
    uses = {}
    task = None
    if category == "Structure preparation":
        uses = {"tile_tmol": "Score, minimize, or FastRelax saved structures on a GPU",
                "tile_complex_prepare": "Prepare protein–ligand complexes for simulation"}
    if category == _ANALYSIS_CATEGORY:
        tasks = {label: descriptions for label, descriptions in _ANALYSIS_TASKS.items()
                 if operations.keys() & descriptions.keys()}
        classified = {kind for descriptions in tasks.values() for kind in descriptions}
        if remaining := operations.keys() - classified:
            tasks["Other analysis tools"] = dict.fromkeys(remaining, "")
        task = st.selectbox("Analysis task", list(tasks), key=key + "-analysis-task")
        uses = tasks[task]
        operations = {kind: spec for kind, spec in operations.items() if kind in uses}
        table_key += f"-{task}"
        if task == _CANDIDATE_SELECTION_TASK:
            st.caption("Connect design or analysis results from earlier in the pipeline. These steps work with the candidates and scores those steps produced.")
    table = ["<table><thead><tr><th scope='col'>Tool or step</th>"
             "<th scope='col'>Pipeline</th></tr></thead><tbody>"]
    order = ([kind for kind in uses if kind in operations] if task == _CANDIDATE_SELECTION_TASK
             else sorted(operations, key=lambda kind: operations[kind]["label"].casefold()))
    for kind in order:
        spec = operations[kind]
        label = escape(spec["label"])
        purpose = f'<small class="operation-use">{escape(uses[kind])}</small>' if uses.get(kind) else ""
        table.append(f'<tr><th scope="row">{label}{purpose}</th>'
                     f'<td><button type="button" data-engine="{escape(kind)}" '
                     f'data-pipeline="true" aria-label="Pipeline Builder · {label}">Pipeline Builder</button></td></tr>')
    choice = _matrix(data="".join([*table, "</tbody></table>"]), key=table_key,
                     on_choose_change=lambda: None).get("choose")
    if not choice or choice[0] not in operations:
        return

    kind = choice[0]
    spec = operations[kind]

    @st.dialog(spec["label"], width="large")
    def details():
        try:
            affinity = category == _DOCKING_CATEGORY and kind == "boltz2"
            _operation_setup(kind, f"{table_key}-{kind}", open_pipeline, prediction=category == _PREDICTION_CATEGORY,
                             predict_affinity=affinity, database_user=database_user, project_id=project_id)
        except (ValueError, OSError) as exc:
            st.error(str(exc))

    details()


def _render_small_molecule_design(key, open_pipeline):
    from mn_molecular_design_studio.tool_connections import editor_tool_catalog

    data = small_molecule_catalogue()
    modes = data["modes"]
    engines = {row["engine"]: row for row in data["items"]}
    matrix = {name: {mode: {"canvas"} for mode in row["modes"]} for name, row in engines.items()}
    choice = _matrix(data=_matrix_markup(matrix, modes, ACCESS, primary_access="canvas"), key=key + "-small-molecule-modes",
                     on_choose_change=lambda: None).get("choose")
    if not choice:
        return
    engine, mode = choice
    row = engines.get(engine)
    if row is None or mode not in row["modes"]:
        return
    operation = row["operation"]

    @st.dialog(f"{engine} · {modes[mode]}", width="large")
    def details():
        try:
            _operation_setup(operation, f"{key}-{operation}-{mode}", open_pipeline, ligand_design_mode=mode)
        except ValueError as exc:
            st.error(str(exc))

    details()


def _render_bioemu2_modes(key, open_pipeline):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.bioemu2 import BIOEMU2_MODES

    modes = {mode: label for mode, (label, _) in BIOEMU2_MODES.items()}
    matrix = {"BioEmu-2": {mode: {"canvas"} for mode in modes}}
    choice = _matrix(data=_matrix_markup(matrix, modes, ACCESS, primary_access="canvas"), key=key + "-bioemu2-modes",
                     on_choose_change=lambda: None).get("choose")
    if not choice or choice[0] != "BioEmu-2" or choice[1] not in modes:
        return
    mode = choice[1]

    @st.dialog(f"BioEmu-2 · {modes[mode]}", width="large")
    def details():
        try:
            st.write(BIOEMU2_MODES[mode][1])
            _operation_setup("tile_bioemu2", f"{key}-bioemu2-{mode}", open_pipeline, bioemu2_mode=mode)
        except ValueError as exc:
            st.error(str(exc))

    details()


def _operation_setup(kind, key, open_pipeline, *, prediction=False, predict_affinity=False,
                     ligand_design_mode="", bioemu2_mode=None, design_task="", database_user=None, project_id=None):
    spec = NODE_CATALOG[kind]
    if not bioemu2_mode:
        st.caption(spec["description"])
    from_saved = False
    if prediction and kind in SAVED_DESIGN_REFOLDERS - {"af2"}:
        from_saved = st.radio("Starting input", ["New molecular inputs", "Saved designs"],
                              horizontal=True, key=key + "-prediction-source") == "Saved designs"
    mode = "protein_ligand" if predict_affinity else _prediction_input(kind, key) if prediction and not from_saved else None
    source = ""
    prepare = False
    sources = source_choices(kind, database_available=bool(database_user and project_id)) if not mode else {}
    if kind == "tile_bioemu2" and bioemu2_mode == "propagate":
        sources = {"": "Upload structures", "project_design_source": "Completed design job"}
    if kind == "tile_tmol":
        sources = {"": "Upload structures", **sources}
    if len(sources) > 1:
        source = st.selectbox("Start from", list(sources), format_func=sources.get, key=key + "-source")
    elif sources:
        source, label = next(iter(sources.items()))
        st.caption(f"Starting input · {label}")
    prepare_source = _database_source(database_user, project_id, key) if source == "molecule_source" else None
    if kind in MD_ENGINES:
        prepare = st.radio("Starting complex", ("prepared", "unprepared"),
                           format_func={"prepared": "I have a prepared complex", "unprepared": "Prepare a complex first"}.get,
                           horizontal=True, key=key + "-prepare") == "unprepared"
        st.caption("Choose the ligand file and its identifier in the simulation step.")
    if kind == "tile_library_centroids":
        prepare = st.radio("Starting library", ("prepared", "unprepared"),
                           format_func={"prepared": "Use prepared fingerprints", "unprepared": "Prepare fingerprints first"}.get,
                           horizontal=True, key=key + "-prepare") == "unprepared"
    if kind in (*DOCKING, "tile_affinity_nesso"):
        st.caption("Select the compounds in the tool step. The Target step supplies the receptor.")
    if kind == "tile_binderscope":
        st.caption("Choose native campaigns in the tool step. For designs saved in the database, use Analyze designs in Data Browser.")
    if kind == "ipsae":
        st.caption("Select results with matching PAE and structure files.")

    graph = build_capability_graph(kind, prediction_mode=mode, predict_affinity=predict_affinity,
                                   ligand_design_mode=ligand_design_mode, bioemu2_mode=bioemu2_mode,
                                   source=source, prepare=prepare, design_task=design_task)
    added = []
    if kind == "tile_complex_prepare":
        optional = _choose_step(graph, "Run a simulation", MD_ENGINES, key + "-simulation")
        if optional:
            graph = add_followup(graph, optional, handoff=True)
            added.append(optional)
    if kind in MD_ENGINES or any(step in MD_ENGINES for step in added):
        optional = _choose_step(graph, "Analyze the trajectory", ("md_analysis",), key + "-trajectory")
        if optional:
            graph = add_followup(graph, optional)
            added.append(optional)
    if kind == "tile_prepare_library":
        optional = _choose_step(graph, "Create centroids from the prepared library", ("tile_library_centroids",), key + "-centroids")
        if optional:
            graph = add_followup(graph, optional)
            added.append(optional)
    if mode == "protein_ligand":
        for task, methods in (("Check pose geometry", ("posebusters",)),
                              ("Map ligand contacts", POSE_ANALYSES[:-1])):
            optional = _choose_step(graph, task, methods, key + "-" + methods[0])
            if optional:
                graph = add_followup(graph, optional)
                added.append(optional)
    if mode in {"protein", "protein_protein"}:
        methods = ("sequence_qc", "proteinsol") if mode == "protein" else ("biotite_interface",)
        task = "Analyze sequence properties" if mode == "protein" else "Measure interface contacts"
        optional = _choose_step(graph, task, methods, key + "-analysis")
        if optional:
            graph = add_followup(graph, optional)
            added.append(optional)
    if ligand_design_mode:
        target = next((node["id"] for node in graph["nodes"] if node["data"]["kind"] == "target"), None)
        if target:
            optional = _choose_step(graph, "Dock the generated compounds", DOCKING, key + "-dock")
            if optional:
                graph = add_followup(graph, optional, target_id=target)
                added.append(optional)
    if design_task == "new_fold" and kind != "rfd4":
        goals = {"backbones": "Generate backbones", "sequences": "Generate backbones and design sequences",
                 "structures": "Generate, design sequences, and predict structures"}
        goal = st.selectbox("Workflow goal", list(goals), format_func=goals.get, key=key + "-goal")
        if goal != "backbones":
            choices = followup_choices(graph, ("ligandmpnn", "esm_if1"))
            sequence = st.selectbox("Sequence design engine", choices,
                                    format_func=lambda method: "ProteinMPNN" if method == "ligandmpnn"
                                    else NODE_CATALOG[method]["label"], key=key + "-sequence-engine")
            graph = add_followup(graph, sequence, handoff=True)
        if goal == "structures":
            choices = followup_choices(graph, ("esmfold", "boltz2", "prediction_openfold3_ob0"))
            refold = st.selectbox("Refolding engine", choices, format_func=lambda method: NODE_CATALOG[method]["label"],
                                  key=key + "-refold-engine")
            graph = add_followup(graph, refold, handoff=True)

    finish = ("To analyze saved docking poses, open the completed run in My results and select the poses."
              if kind in DOCKING or any(step in DOCKING for step in added) else "")
    if kind in {"shared_filter", "tile_validation_acceptance"}:
        finish = "After the run, use the recorded acceptance stage in My results to collect passing designs."
    _render_plan(graph, key, open_pipeline, optional=added, finish=finish, prepare_source=prepare_source,
                 ready=source != "molecule_source" or prepare_source is not None)


def _database_source(user, project_id, key):
    from mn_molecular_design_studio.database_selections import freeze_protein_selection
    from mn_molecular_design_studio.molecular_database import MolecularDatabase

    query = st.text_input("Find saved entries by name", placeholder="Enter at least two characters",
                          key=key + "-database-query").strip()
    if len(query) < 2:
        return None
    database = MolecularDatabase()
    page = database.list_molecules(user, project_id=project_id, query=query, limit=100)
    choices = {row["id"]: row["name"] for row in page["items"]
               if any(item["type"] == "protein" and item.get("sequence") for item in row["components"])}
    if not choices:
        st.info("No entries with protein sequences match this search.")
        return None
    if len(page["items"]) == 100:
        st.caption("Showing the first 100 matches. Refine the search to find more specific entries.")
    selected = st.multiselect("Saved entries", choices, format_func=choices.get, key=key + "-database-selection")
    if not selected:
        return None

    def prepare(graph):
        source = next(node for node in graph["nodes"] if node["data"]["kind"] == "molecule_source")
        snapshot = freeze_protein_selection(database, user, selected, project_id)
        source["data"]["params"]["source_run_dir"] = str(snapshot)

    return prepare


def _choose_step(graph, task, methods, key):
    methods = followup_choices(graph, methods)
    if not methods or not st.checkbox(task, key=key):
        return None
    if len(methods) == 1:
        st.caption(NODE_CATALOG[methods[0]]["label"])
        return methods[0]
    return st.selectbox(task + " engine", methods, format_func=lambda kind: NODE_CATALOG[kind]["label"],
                        key=key + "-method")


def _render_plan(graph, key, open_pipeline, *, optional=(), finish="", prepare_source=None, ready=True):
    from mn_molecular_design_studio.pipeline_blocks import program

    plan = program(graph)
    graph_kinds = {node["id"]: node["data"]["kind"] for node in graph["nodes"]}
    st.markdown("**Steps to build**")
    for stage in plan["stages"]:
        for identifier in stage:
            step = plan["blocks"][identifier]
            role = "Optional · selected" if graph_kinds[identifier] in optional else "Required"
            upstream = (" · From " + ", ".join(step["sources"])) if step["sources"] else ""
            st.write(f"**{step['label']}** · {role}{upstream}")
    if finish:
        st.caption(finish)
    st.caption("Set the inputs on the canvas after creating the pipeline.")
    if st.button("Create pipeline", type="primary", key=key + "-create-pipeline", disabled=not ready):
        if prepare_source:
            prepare_source(graph)
        open_pipeline(graph)


def _prediction_input(kind, key):
    from mn_molecular_design_studio.scientific_tasks import PREDICTION_MODES, prediction_methods

    modes = {mode: label for mode, label in PREDICTION_MODES.items() if kind in prediction_methods(mode)}
    if not modes:
        return None
    return st.selectbox("Prediction goal", list(modes), format_func=modes.get, key=key + "-molecular-mode")
