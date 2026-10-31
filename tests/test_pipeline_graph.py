"""Visual graph contracts: dependency execution, immutable versions and reuse."""
import copy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from mn_molecular_design_studio import pipeline_graph as graph
from mn_molecular_design_studio.engine.workflow_engine.workflows import graph_pipeline as runner


def node(kind, identifier):
    spec = graph.NODE_CATALOG[kind]
    return {"id": identifier, "type": "studio", "position": {"x": 0, "y": 0}, "data": {"kind":kind,"label":spec["label"],"params":{f["name"]:f["default"] for f in spec["fields"]}}}


def edge(source, target):
    return {"id": f"{source}-{target}", "source":source, "target":target,"sourceHandle":"out","targetHandle":"in"}


def test_foundry_template_has_explicit_generation_sequence_and_prediction():
    plan=graph.validate_graph(graph.new_graph())
    kinds = [n["data"]["kind"] for n in plan["graph"]["nodes"]]
    assert kinds == ["target", "binder_context", "rfd3", "foundry_mpnn", "rf3", "results"]
    assert len(plan["order"]) == 6


def test_ipsae_analysis_accepts_design_only_and_design_to_cofold_inputs():
    from mn_molecular_design_studio.pipeline_node_ui import _analysis_input_profile

    for include_cofolder, expected_label in ((False, "Design structures / sequences"),
                                             (True, "Design → cofold path")):
        nodes = [node("target", "target"), node("rfd4", "design"),
                 node("ipsae", "analysis"), node("results", "results")]
        edges = [edge("target", "design")]
        if include_cofolder:
            nodes.insert(2, node("boltz2", "cofold"))
            edges.extend([edge("design", "cofold"), edge("cofold", "analysis")])
        else:
            edges.append(edge("design", "analysis"))
        edges.append(edge("analysis", "results"))
        draft = graph.new_graph("blank")
        draft["nodes"], draft["edges"] = nodes, edges
        validated = graph.validate_graph(draft)
        assert "analysis" in validated["order"]
        assert _analysis_input_profile(validated["graph"], "analysis")[0] == expected_label


def test_ipsae_accepts_separate_design_and_cofold_inputs_and_profiles_identity_join():
    from mn_molecular_design_studio.pipeline_node_ui import _analysis_input_profile

    g = graph.new_graph("blank")
    g["nodes"] = [node(kind, identifier) for kind, identifier in [
        ("target", "target"), ("rfd4", "design"), ("boltz2", "cofold"),
        ("ipsae", "analysis"), ("results", "results"),
    ]]
    g["edges"] = [edge("target", "design"), edge("design", "cofold"),
                  edge("design", "analysis"), edge("cofold", "analysis"), edge("analysis", "results")]
    validated = graph.validate_graph(g)
    assert _analysis_input_profile(validated["graph"], "analysis")[0] == \
        "Separate design + cofold inputs · matched by candidate identity"


def test_ipsae_fanin_joins_by_exact_candidate_lineage_and_marks_unmatched(tmp_path, design_source_factory):
    source_dir, design = design_source_factory()
    design = copy.deepcopy(design)
    design.update(candidate_id="rfd4-design-1", stage="generation.backbone_sequence",
                  complex_pdb=str(source_dir / "designed.cif"), binder_chains=["B"], target_chains=["A"])
    unrelated = copy.deepcopy(design)
    unrelated["candidate_id"] = "other-design"
    cofolded = copy.deepcopy(design)
    cofolded.update(candidate_id="of3-prediction-1", stage="complex_refolding",
                     parents=[design["candidate_id"]],
                     raw_metadata={"source_candidate": {**design, "candidate_id": "stale-lineage"}})
    orphan = copy.deepcopy(cofolded)
    unrelated_prediction_structure = tmp_path / "orphan.cif"
    unrelated_prediction_structure.write_text("prediction")
    orphan.update(candidate_id="of3-orphan", parents=["missing-design"], complex_pdb=str(unrelated_prediction_structure),
                  raw_metadata={"source_candidate": {**design, "candidate_id": "stale-lineage"},
                                "design_structure_path": str(source_dir / "designed.cif")})

    design_output = runner._local_source(tmp_path, "design", [design, unrelated])
    cofold_output = runner._local_source(tmp_path, "cofold", [cofolded, orphan])
    nodes = {"design": node("rfd4", "design"), "cofold": node("boltz2", "cofold")}
    joined = runner._ipsae_design_cofold_inputs(
        tmp_path, "analysis", ["design", "cofold"], [design_output, cofold_output], nodes)
    rows = runner._absolute_candidates(joined)
    matched, unmatched = rows

    assert matched["candidate_id"] == "of3-prediction-1"
    assert matched["raw_metadata"]["analysis_design_pairing"]["status"] == "matched"
    assert matched["raw_metadata"]["source_candidate"]["candidate_id"] == "rfd4-design-1"
    assert unmatched["raw_metadata"]["analysis_design_pairing"]["status"] == "unmatched"
    assert "source_candidate" not in unmatched["raw_metadata"]
    assert unmatched["raw_metadata"]["unmatched_upstream_source_candidate"]["candidate_id"] == "stale-lineage"
    assert joined["analysis_input_join"]["matched_count"] == 1
    assert joined["analysis_input_join"]["unmatched_count"] == 1
    from mn_molecular_design_studio.engine.workflow_engine.workflows.analysis import _candidate_design_reference
    assert _candidate_design_reference(Path(tmp_path), unmatched) == (None, [])


def test_graph_routes_rfd4_and_independent_cofold_to_ipsae(monkeypatch, tmp_path, design_source_factory):
    from types import SimpleNamespace
    from mn_molecular_design_studio import workspace

    source_dir, base = design_source_factory()
    pipeline_run = tmp_path / "pipeline"
    pipeline_run.mkdir()
    (pipeline_run / "metadata.json").write_text(json.dumps({"status": "queued", "job_type": "visual_pipeline"}))
    (pipeline_run / "input.json").write_text(json.dumps({"tool": "visual_pipeline", "params": {"pipeline_name": "rfd4-of3-analysis"}}))
    monkeypatch.setattr(workspace, "mark_studio_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(graph, "app_home", lambda: tmp_path / "home")
    monkeypatch.setattr(runner, "create_job", lambda *args, **kwargs: SimpleNamespace(run_dir=pipeline_run))

    g = graph.new_graph("blank")
    g["nodes"] = [node(kind, identifier) for kind, identifier in [
        ("target", "target"), ("rfd4", "design"), ("boltz2", "cofold"),
        ("ipsae", "analysis"), ("results", "results"),
    ]]
    next(item for item in g["nodes"] if item["id"] == "target")["data"]["params"].update(
        target_pdb=str(source_dir / "designed.cif"), target_chains="A")
    g["edges"] = [edge("target", "design"), edge("design", "cofold"),
                  edge("design", "analysis"), edge("cofold", "analysis"), edge("analysis", "results")]

    captured = {}
    original_request = runner.node_request

    def node_request(kind, params, input_source, name):
        if kind == "ipsae":
            captured["analysis_rows"] = runner._absolute_candidates(input_source)
            return {"workflow": graph.NODE_CATALOG[kind]["workflow"], "parameters": {
                "source_run_dir": input_source["run_dir"], "candidates_jsonl": input_source["candidates_jsonl"],
                "job_name": name,
            }}
        return {"workflow": graph.NODE_CATALOG[kind]["workflow"], "parameters": {
            "source_run_dir": input_source.get("run_dir", ""),
            "candidates_jsonl": input_source.get("candidates_jsonl", ""),
        }}

    monkeypatch.setattr(runner, "node_request", node_request)
    submitted = []

    def submit(request, parent, current_node, **kwargs):
        kind = current_node["data"]["kind"]
        submitted.append(kind)
        child = tmp_path / f"child-{kind}"
        child.mkdir()
        if kind == "ipsae":
            result = {"success": True, "outputs": {"tool_status": {"ipsae": {"success": True}}}}
        else:
            candidate = copy.deepcopy(base)
            candidate.update(
                candidate_id="rfd4-design-1" if kind == "rfd4" else "of3-prediction-1",
                complex_pdb=str(source_dir / "designed.cif"),
                stage="generation.backbone_sequence" if kind == "rfd4" else "complex_refolding",
                parents=[] if kind == "rfd4" else ["rfd4-design-1"],
            )
            normalized = child / "artifacts/normalized_candidates"
            normalized.mkdir(parents=True)
            (normalized / "candidates.jsonl").write_text(json.dumps(candidate) + "\n")
            result = {"success": True}
        (child / "metadata.json").write_text(json.dumps({
            "status": "completed", "job_code": kind, "task_group": "test",
            "pipeline_parent_run_dir": str(parent),
        }))
        (child / "result.json").write_text(json.dumps(result))
        (child / "input.json").write_text(json.dumps({"params": {}}))
        return child

    monkeypatch.setattr(runner, "_submit", submit)
    try:
        result_dir = runner.run_graph_pipeline(g, "rfd4-of3-analysis", "owner", "project")
    finally:
        monkeypatch.setattr(runner, "node_request", original_request)

    assert submitted == ["rfd4", "boltz2", "ipsae"]
    row = captured["analysis_rows"][0]
    assert row["candidate_id"] == "of3-prediction-1"
    assert row["raw_metadata"]["source_candidate"]["candidate_id"] == "rfd4-design-1"
    assert row["raw_metadata"]["analysis_design_pairing"]["status"] == "matched"
    results = json.loads((result_dir / "result.json").read_text())
    assert results["success"] is True


def test_supported_campaign_tile_settings_survive_graph_validation_and_use_normal_workflow_requests():
    from mn_molecular_design_studio.distributed_campaign_config_ui import _default_validation_stages

    campaign = {
        "enabled": True,
        "accepted_target": 2,
        "max_attempts": 6,
        "acceptance_level": "validated",
        "batch_size": 1,
        "max_concurrent_batches": 2,
        "max_concurrent_validations": 3,
        "workstation_selection": "eligible",
        "gpu_selection": "auto",
        "validation_stages": _default_validation_stages(),
    }
    target = node("target", "target-1")
    target["data"]["params"].update(target_pdb="/tmp/pd-l1.pdb", target_chains="A")
    generator = node("rfd4", "rfd4-1")
    generator["data"]["params"]["distributed_campaign"] = copy.deepcopy(campaign)
    saved_graph = graph.clean_graph({
        "schema_version": graph.SCHEMA,
        "id": "a" * 32,
        "name": "Retained campaign fixture",
        "nodes": [target, generator],
        "edges": [edge("target-1", "rfd4-1")],
        "viewport": {"x": 0, "y": 0, "zoom": 1},
    })
    saved_node = next(item for item in saved_graph["nodes"] if item["id"] == "rfd4-1")
    assert saved_node["data"]["params"]["distributed_campaign"]["generator_id"] == "rfd4_proteina"
    assert saved_node["data"]["params"]["distributed_campaign"]["max_concurrent_batches"] == 2
    assert saved_node["data"]["params"]["distributed_campaign"]["max_concurrent_validations"] == 3

    request = runner.node_request(
        "rfd4", saved_node["data"]["params"],
        {"target": {"target_pdb": "/tmp/pd-l1.pdb", "target_chains": ["A"], "hotspots": "A40"}},
        "RFD4 tile campaign",
    )
    assert request["workflow"] == "design.rfd4_proteina"
    assert request["parameters"]["distributed_campaign"]["generator_id"] == "rfd4_proteina"
    assert request["parameters"]["distributed_campaign"]["accepted_target"] == 2
    assert request["parameters"].get("batch_size", 1) == 1


def test_bindcraft2_campaign_tile_keeps_authoritative_native_policy():
    campaign = {
        "enabled": True,
        "accepted_target": 1,
        "max_attempts": 2,
        "acceptance_level": "native",
        "batch_size": 1,
        "max_concurrent_batches": 2,
        "max_concurrent_validations": 3,
        "workstation_selection": "auto",
        "gpu_selection": "auto",
        "validation_stages": [],
    }
    target = node("target", "target-1")
    bindcraft = node("bindcraft2", "bindcraft2-1")
    bindcraft["data"]["params"]["distributed_campaign"] = campaign
    cleaned = graph.clean_graph({
        "schema_version": graph.SCHEMA,
        "id": "c" * 32,
        "name": "BindCraft campaign fixture",
        "nodes": [target, bindcraft], "edges": [edge("target-1", "bindcraft2-1")],
        "viewport": {"x": 0, "y": 0, "zoom": 1},
    })
    params = next(item["data"]["params"] for item in cleaned["nodes"] if item["data"]["kind"] == "bindcraft2")
    request = runner.node_request("bindcraft2", params,
        {"target": {"target_pdb": "/tmp/pd-l1.pdb", "target_chains": ["A"]}},
        "BindCraft2 tile campaign")
    assert request["workflow"] == "design.bindcraft2"
    assert request["parameters"]["distributed_campaign"]["generator_id"] == "bindcraft2"
    assert request["parameters"]["distributed_campaign"]["acceptance_level"] == "native"


@pytest.mark.parametrize("mode", ["native_complete", "native_comparative", "modular_design"])
def test_shared_objective_normalizes_two_generator_branches_with_common_target_and_metrics(mode):
    from mn_molecular_design_studio.distributed_campaign_config_ui import _default_validation_stages

    target = node("target", "target")
    bindcraft = node("bindcraft2", "bindcraft")
    rfd4 = node("rfd4", "rfd4")
    results = node("results", "results")
    saved = graph.validate_graph({
        "schema_version": graph.SCHEMA, "id": "d" * 32, "name": "Shared objective",
        "workflow_mode": mode, "nodes": [target, bindcraft, rfd4, results],
        "edges": [edge("target", "bindcraft"), edge("target", "rfd4"),
                  edge("bindcraft", "results"), edge("rfd4", "results")],
        "shared_campaign": {
            "enabled": True, "generator_nodes": ["bindcraft", "rfd4"],
            "accepted_target": 3, "max_attempts": 12,
            "max_concurrent_jobs": 4, "max_concurrent_generation_jobs": 3,
            "max_concurrent_validations": 2,
            "validation_stages": _default_validation_stages(),
            "generator_limits": {"rfd4": {"max_attempts": 5, "max_concurrent_generation_jobs": 1}},
        },
        "viewport": {"x": 0, "y": 0, "zoom": 1},
    })["graph"]
    shared = saved["shared_campaign"]
    assert shared["accepted_target"] == 3
    assert shared["generator_nodes"] == ["bindcraft", "rfd4"]
    assert shared["generator_limits"]["rfd4"]["max_attempts"] == 5
    assert "bindcraft" not in shared["generator_limits"]
    assert shared["workstation_selection"] == "auto"
    assert shared["validation_stages"][0]["parameters"]["engine"] == "boltz2"
    edited = copy.deepcopy(saved)
    edited["shared_campaign"]["max_concurrent_generation_jobs"] = 1
    edited["shared_campaign"]["generator_limits"] = {}
    assert graph.validate_graph(edited)["graph"]["shared_campaign"]["generator_limits"] == {}
    edited["shared_campaign"]["workstation_selection"] = "local"
    assert graph.validate_graph(edited)["graph"]["shared_campaign"]["workstation_selection"] == "local"
    edited["shared_campaign"]["workstation_selection"] = "unsupported"
    with pytest.raises(ValueError, match="workstation_selection"):
        graph.validate_graph(edited)


@pytest.mark.parametrize("mode,kinds", [
    ("native_complete", ("bindcraft2", "boltzgen")),
    ("native_comparative", ("bindcraft2", "pxdesign_pipeline")),
    ("modular_design", ("rfd3", "genie3")),
])
def test_generated_shared_objective_accepts_registered_fixed_size_native_branches(mode, kinds):
    target = node("target", "target")
    branches = [node(kind, f"design-{index}") for index, kind in enumerate(kinds)]
    results = node("results", "results")
    data = {
        "schema_version": graph.SCHEMA, "id": "f" * 32, "name": "Generated objective",
        "workflow_mode": mode, "nodes": [target, *branches, results],
        "edges": [edge("target", branch["id"]) for branch in branches]
                + [edge(branch["id"], "results") for branch in branches],
        "shared_campaign": {
            "enabled": True, "objective_type": "generated",
            "generator_nodes": [branch["id"] for branch in branches],
            "generated_target": 2, "max_attempts": 20,
            "max_concurrent_jobs": 4, "max_concurrent_generation_jobs": 3,
            "max_concurrent_validations": 2, "validation_stages": [],
        },
        "viewport": {"x": 0, "y": 0, "zoom": 1},
    }
    validated = graph.validate_graph(data)["graph"]
    assert validated["shared_campaign"]["objective_type"] == "generated"
    assert validated["shared_campaign"]["generated_target"] == 2
    assert validated["shared_campaign"]["accepted_target"] is None
    assert validated["shared_campaign"]["generator_nodes"] == [branch["id"] for branch in branches]


def test_generated_campaign_runs_filters_reviews_and_exports_through_graph_results(
    tmp_path, monkeypatch, design_source_factory,
):
    from mn_molecular_design_studio import candidate_pipeline, candidate_sets
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json
    from mn_molecular_design_studio.workflow_modes import capability
    from mn_molecular_design_studio import workspace

    target_run, _ = design_source_factory(tmp_path / 'target-source')
    parent = tmp_path / 'generated-campaign'
    parent.mkdir()
    write_json(parent / 'metadata.json', {'status': 'queued', 'job_type': 'visual_pipeline'})
    write_json(parent / 'input.json', {'params': {'pipeline_name': 'Generated PD-L1 fixture'}})
    monkeypatch.setattr(graph, 'app_home', lambda: tmp_path / 'studio-home')
    monkeypatch.setattr(workspace, 'mark_studio_job', lambda *args, **kwargs: None)
    monkeypatch.setattr(runner, 'create_job', lambda *args, **kwargs: SimpleNamespace(run_dir=parent))
    monkeypatch.setattr(candidate_pipeline, 'preflight', lambda *args, **kwargs: None)

    target = node('target', 'target')
    target['data']['params'].update(target_pdb=str(target_run / 'designed.cif'), target_chains='A')
    generators = [node(kind, identifier) for kind, identifier in (
        ('rfd3', 'generator-rfd3'), ('genie3', 'generator-genie3'),
    )]
    for generator in generators:
        count_field = capability(generator['data']['kind'])['design_count_field']
        generator['data']['params'][count_field] = 2
    collection = node('candidate_merge', 'collection')
    scientific_filter = node('shared_filter', 'scientific-filter')
    scientific_filter['data']['params']['rules'] = [{
        'field': 'common.confidence.generator.complex.iptm', 'operator': '>=', 'value': 0.6,
    }]
    results = node('results', 'results')
    graph_data = {
        'schema_version': graph.SCHEMA, 'id': 'a' * 32, 'name': 'Generated PD-L1 fixture',
        'workflow_mode': 'modular_design',
        'nodes': [target, *generators, collection, scientific_filter, results],
        'edges': [*(edge('target', generator['id']) for generator in generators),
                  *(edge(generator['id'], 'collection') for generator in generators),
                  edge('collection', 'scientific-filter'), edge('scientific-filter', 'results')],
        'shared_campaign': {
            'enabled': True, 'objective_type': 'generated',
            'generator_nodes': [generator['id'] for generator in generators],
            'generated_target': 4, 'max_attempts': 4,
            'max_concurrent_jobs': 2, 'max_concurrent_generation_jobs': 2,
            'max_concurrent_validations': 1, 'validation_stages': [],
            'generator_limits': {generator['id']: {'target_count': 2} for generator in generators},
        },
        'viewport': {'x': 0, 'y': 0, 'zoom': 1},
    }
    submissions = []

    def submit(request, graph_parent, branch, **kwargs):
        kind = branch['data']['kind']
        outputs = target_run / ('native-output-' + branch['id'])
        (outputs / 'artifacts/normalized_candidates').mkdir(parents=True)
        native_rows = [{
            'candidate_id': f"{branch['id']}-native-{index}", 'source_tool': capability(kind)['tool'],
            'stage': 'generation.backbone_sequence', 'binder_sequence': 'KLMNP',
            'complex_pdb': str(target_run / 'designed.cif'), 'target_pdb': str(target_run / 'designed.cif'),
            'target_chains': ['A'], 'binder_chains': ['B'],
            'metrics': {'native_pass_filters': True, 'iptm': iptm}, 'raw_metadata': {},
        } for index, iptm in enumerate((0.85, 0.35))]
        (outputs / 'artifacts/normalized_candidates/candidates.jsonl').write_text(
            ''.join(json.dumps(row) + '\n' for row in native_rows))
        write_json(outputs / 'metadata.json', {'status': 'completed', 'tool': capability(kind)['tool'],
            'job_code': f'{kind}-1', 'pipeline_parent_run_dir': str(graph_parent)})
        write_json(outputs / 'input.json', {'params': request['parameters']})
        write_json(outputs / 'result.json', {'success': True, 'metrics': {'candidate_count': 2}})
        submissions.append(branch['id'])
        return outputs

    monkeypatch.setattr(runner, '_submit', submit)

    run = runner.run_graph_pipeline(graph_data, 'Generated PD-L1 fixture', 'owner', 'project', reuse_completed=False)
    result = read_json(run / 'result.json')
    ledger = read_json(run / 'artifacts/shared_campaign/ledger.json')
    summary = result['outputs']['workflow_summary']['shared_campaign']
    rows = result['outputs']['candidates']
    assert result['success'] is True
    assert len(rows) == 4  # Failed candidates remain in review Results.
    assert set(submissions) == {generator['id'] for generator in generators}
    assert summary['generated_designs'] == 4 and summary['status'] == 'target_reached'
    assert summary['generated_native_pass'] == 4
    assert summary['scientifically_passing_designs'] == 2
    assert summary['scientifically_rejected_designs'] == 2
    assert summary['scheduler_jobs'] == 2 and summary['scheduler_jobs_active'] == 0
    assert summary['compute_seconds'] is None
    assert len(ledger['participants']) == 4 and len(ledger['evaluations']) == 4
    for generator in generators:
        branch_summary = summary['per_generator'][generator['id']]
        assert branch_summary['candidates_generated'] == 2
        assert branch_summary['scientifically_passing_designs'] == 1
        assert branch_summary['scientifically_rejected_designs'] == 1
    assert all(candidate_sets.metadata(row)['shared_filter']['status'] in {'passed', 'failed'} for row in rows)
    assert (run / 'artifacts/normalized_candidates/final_evidence.csv').is_file()
    assert (run / 'artifacts/candidate_set/manifest.json').is_file()


def test_shared_objective_rejects_different_target_and_incompatible_downstream_policy():
    from mn_molecular_design_studio.distributed_campaign_config_ui import _default_validation_stages

    target_a, target_b = node("target", "target-a"), node("target", "target-b")
    bindcraft, rfd4, results = node("bindcraft2", "bindcraft"), node("rfd4", "rfd4"), node("results", "results")
    common = {"enabled": True, "generator_nodes": ["bindcraft", "rfd4"], "accepted_target": 1,
              "max_attempts": 3, "max_concurrent_jobs": 2, "max_concurrent_generation_jobs": 2,
              "max_concurrent_validations": 2, "validation_stages": _default_validation_stages()}
    graph_data = {"schema_version": graph.SCHEMA, "id": "e" * 32, "name": "Incompatible shared objective",
                  "nodes": [target_a, target_b, bindcraft, rfd4, results],
                  "edges": [edge("target-a", "bindcraft"), edge("target-b", "rfd4"),
                            edge("bindcraft", "results"), edge("rfd4", "results")],
                  "shared_campaign": common, "viewport": {"x": 0, "y": 0, "zoom": 1}}
    with pytest.raises(ValueError, match="same target structure"):
        graph.validate_graph(graph_data)

    incompatible = copy.deepcopy(graph_data)
    incompatible["nodes"] = [item for item in incompatible["nodes"] if item["id"] != "target-b"]
    incompatible["edges"] = [edge("target-a", "bindcraft"), edge("target-a", "rfd4"),
                              edge("bindcraft", "results"), edge("rfd4", "results")]
    branch_stages = _default_validation_stages()
    branch_stages[0]["rules"] = [{"field": "iptm", "operator": ">=", "value": 0.75}]
    bindcraft = next(item for item in incompatible["nodes"] if item["id"] == "bindcraft")
    bindcraft["data"]["params"]["distributed_campaign"] = {
        "enabled": True, "accepted_target": 1, "max_attempts": 3,
        "acceptance_level": "validated", "validation_stages": branch_stages,
    }
    with pytest.raises(ValueError, match="differ from the shared"):
        graph.validate_graph(incompatible)


@pytest.mark.parametrize("objective_type", ["accepted", "time"])
def test_graph_runner_injects_one_shared_campaign_context_into_normal_branch_requests(tmp_path, monkeypatch, design_source_factory, objective_type):
    from mn_molecular_design_studio import workspace
    from mn_molecular_design_studio.distributed_campaign_config_ui import _default_validation_stages

    source, template_candidate = design_source_factory()
    monkeypatch.setattr(workspace, "mark_studio_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(graph, "app_home", lambda: tmp_path / "studio-home")
    parent = tmp_path / "graph-parent"
    parent.mkdir()
    (parent / "metadata.json").write_text(json.dumps({"status": "queued", "job_type": "visual_pipeline"}))
    (parent / "input.json").write_text(json.dumps({"params": {"pipeline_name": "Shared graph"}}))
    monkeypatch.setattr(runner, "create_job", lambda *args, **kwargs: SimpleNamespace(run_dir=parent))

    target = node("target", "target")
    target["data"]["params"].update(target_pdb=str(source / "designed.cif"), target_chains="A")
    bindcraft, rfd4, results = node("bindcraft2", "bindcraft"), node("rfd4", "rfd4"), node("results", "results")
    graph_data = {
        "schema_version": graph.SCHEMA, "id": "9" * 32, "name": "Shared graph",
        "workflow_mode": "native_comparative",
        "nodes": [target, bindcraft, rfd4, results],
        "edges": [edge("target", "bindcraft"), edge("target", "rfd4"),
                  edge("bindcraft", "results"), edge("rfd4", "results")],
        "shared_campaign": {
            "enabled": True, "generator_nodes": ["bindcraft", "rfd4"],
            "accepted_target": 1, "max_attempts": 4,
            "max_concurrent_jobs": 2, "max_concurrent_generation_jobs": 2,
            "max_concurrent_validations": 1, "validation_stages": _default_validation_stages(),
        },
        "viewport": {"x": 0, "y": 0, "zoom": 1},
    }
    if objective_type == "time":
        graph_data['shared_campaign'].update(objective_type='time', max_runtime_seconds=300, validation_stages=[])
        graph_data['shared_campaign'].pop('accepted_target')
    requests = []
    children = []

    def submit(request, graph_parent, current_node, **kwargs):
        requests.append(copy.deepcopy(request))
        child = tmp_path / f"child-{current_node['id']}"
        child.mkdir()
        candidate = copy.deepcopy(template_candidate)
        candidate.update(candidate_id=f"{current_node['id']}-candidate", binder_sequence=(
            "ACDEFGHIKLMNPQRSTVWY" if current_node["id"] == "bindcraft" else "MNPQRSTVWYACDEFGHIKL"
        ))
        normalized = child / "artifacts/normalized_candidates"
        normalized.mkdir(parents=True)
        (normalized / "candidates.jsonl").write_text(json.dumps(candidate) + "\n")
        (child / "input.json").write_text(json.dumps({"params": request.get("parameters", {})}))
        (child / "metadata.json").write_text(json.dumps({"status": "completed", "job_code": child.name}))
        (child / "result.json").write_text(json.dumps({"success": True, "metrics": {}}))
        children.append(child)
        return child

    monkeypatch.setattr(runner, "_submit", submit)
    result = runner.run_graph_pipeline(graph_data, "shared-graph", "owner", "project", reuse_completed=False)

    assert [request["parameters"]["distributed_campaign"]["_shared_objective"]["generator_node"]
            for request in requests] == ["bindcraft", "rfd4"]
    assert all(request["parameters"]["distributed_campaign"]["accepted_target"] == (None if objective_type == "time" else 1) for request in requests)
    if objective_type == "time":
        assert all(request["parameters"]["distributed_campaign"]["validation_stages"] == [] for request in requests)
        assert len({request["parameters"]["distributed_campaign"]["deadline_at"] for request in requests}) == 1
    assert all(request["parameters"]["distributed_campaign"]["workstation_selection"] == "auto" for request in requests)
    shared_path = result / "artifacts/shared_campaign/ledger.json"
    assert shared_path.is_file()
    state = json.loads(shared_path.read_text())
    assert state["config"]["generator_nodes"] == ["bindcraft", "rfd4"]


def test_shared_objective_results_preserve_successful_branch_when_other_generator_fails(
    tmp_path, monkeypatch, design_source_factory,
):
    from mn_molecular_design_studio import workspace
    from mn_molecular_design_studio.distributed_campaign_config_ui import _default_validation_stages

    source, template_candidate = design_source_factory()
    monkeypatch.setattr(workspace, "mark_studio_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(graph, "app_home", lambda: tmp_path / "studio-home")
    parent = tmp_path / "graph-parent"
    parent.mkdir()
    (parent / "metadata.json").write_text(json.dumps({"status": "queued", "job_type": "visual_pipeline"}))
    (parent / "input.json").write_text(json.dumps({"params": {"pipeline_name": "Shared partial Results"}}))
    monkeypatch.setattr(runner, "create_job", lambda *args, **kwargs: SimpleNamespace(run_dir=parent))

    target = node("target", "target")
    target["data"]["params"].update(target_pdb=str(source / "designed.cif"), target_chains="A")
    bindcraft, rfd4, results = node("bindcraft2", "bindcraft"), node("rfd4", "rfd4"), node("results", "results")
    graph_data = {
        "schema_version": graph.SCHEMA, "id": "8" * 32, "name": "Shared partial Results",
        "workflow_mode": "native_comparative",
        "nodes": [target, bindcraft, rfd4, results],
        "edges": [edge("target", "bindcraft"), edge("target", "rfd4"),
                  edge("bindcraft", "results"), edge("rfd4", "results")],
        "shared_campaign": {
            "enabled": True, "generator_nodes": ["bindcraft", "rfd4"],
            "accepted_target": 1, "max_attempts": 2,
            "max_concurrent_jobs": 1, "max_concurrent_generation_jobs": 1,
            "max_concurrent_validations": 1,
            "validation_stages": _default_validation_stages(),
        },
        "viewport": {"x": 0, "y": 0, "zoom": 1},
    }

    def submit(_request, _parent, current_node, **_kwargs):
        if current_node["id"] == "bindcraft":
            raise RuntimeError("simulated BindCraft branch startup failure")
        child = tmp_path / "child-rfd4"
        child.mkdir()
        candidate = copy.deepcopy(template_candidate)
        candidate["candidate_id"] = "rfd4-candidate"
        normalized = child / "artifacts/normalized_candidates"
        normalized.mkdir(parents=True)
        (normalized / "candidates.jsonl").write_text(json.dumps(candidate) + "\n")
        (child / "input.json").write_text(json.dumps({"params": {}}))
        (child / "metadata.json").write_text(json.dumps({"status": "completed", "job_code": "RFD4"}))
        (child / "result.json").write_text(json.dumps({"success": True, "metrics": {}}))
        return child

    monkeypatch.setattr(runner, "_submit", submit)
    result_dir = runner.run_graph_pipeline(graph_data, "shared-partial-results", "owner", "project", reuse_completed=False)

    result = json.loads((result_dir / "result.json").read_text())
    results_node = result["outputs"]["pipeline_nodes"]["results"]
    assert result["success"] is True
    assert result["metrics"]["candidate_count"] == 1
    assert result["outputs"]["comparison_status"] == "partial"
    assert results_node["status"] == "completed"
    assert results_node["partial"] is True
    assert results_node["output"]["candidate_count"] == 1
    assert results_node["output"]["partial_pool"]["missing_branches"] == ["bindcraft"]


def test_shared_campaign_child_is_started_after_graph_lineage_is_persisted(tmp_path, monkeypatch):
    from mn_molecular_design_studio import workspace
    from mn_molecular_design_studio.engine.workflow_engine.core import local_worker
    from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation

    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    parent, child = tmp_path / "graph-parent", tmp_path / "shared-child"
    parent.mkdir(); child.mkdir()
    (parent / "input.json").write_text(json.dumps({"params": {"pipeline_name": "Shared graph"}}))
    (parent / "metadata.json").write_text(json.dumps({"studio_development_mode": True}))
    (child / "metadata.json").write_text(json.dumps({"status": "queued", "job_code": "child"}))
    submitted, started = [], []

    def submit(path, *, defer_dispatch=False):
        submitted.append((json.loads(Path(path).read_text()), defer_dispatch))
        return {"job": {"run_dir": str(child)}}

    def start(run_dir):
        metadata = json.loads((Path(run_dir) / "metadata.json").read_text())
        started.append(metadata)

    monkeypatch.setattr(workflow_automation, "submit_workflow_request", submit)
    monkeypatch.setattr(workspace, "mark_studio_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(local_worker, "spawn_worker_for_run", start)
    request = {"workflow": "design.bindcraft2", "parameters": {"distributed_campaign": {
        "enabled": True,
        "_shared_objective": {"ledger_path": str(parent / "artifacts/shared_campaign/ledger.json"), "generator_node": "bindcraft"},
    }}}

    assert runner._submit(
        request, parent, {"id": "bindcraft"}, owner_id="owner", project_id="project", gpu_device="all",
    ) == child
    assert len(submitted) == 1 and submitted[0][1] is True
    assert len(started) == 1
    assert started[0]["pipeline_parent_run_dir"] == str(parent)
    assert started[0]["graph_node_id"] == "bindcraft"


@pytest.mark.parametrize("workflow", ["design.bindcraft2", "design.rfd4_proteina"])
def test_distributed_tile_submission_uses_cpu_controller_and_normal_workflow_api(tmp_path, monkeypatch, workflow):
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    parent = tmp_path / "parent"
    child = tmp_path / "child"
    parent.mkdir()
    child.mkdir()
    (parent / "input.json").write_text(json.dumps({"params": {"pipeline_name": "Fixture"}}))
    (parent / "metadata.json").write_text(json.dumps({"studio_development_mode": True}))
    (child / "metadata.json").write_text(json.dumps({"status": "queued", "tool": workflow.rsplit(".", 1)[-1]}))
    submitted = []
    submission_context = []

    def submit_workflow_request(path):
        from mn_molecular_design_studio.engine.workflow_engine.core.jobs import _SUBMISSION_TAGS
        submission_context.append(dict(_SUBMISSION_TAGS.get()))
        submitted.append(json.loads(Path(path).read_text()))
        return {"job": {"run_dir": str(child)}}

    from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation
    from mn_molecular_design_studio import workspace
    monkeypatch.setattr(workflow_automation, "submit_workflow_request", submit_workflow_request)
    monkeypatch.setattr(workspace, "mark_studio_job", lambda result, **kwargs: result)

    request = {"workflow": workflow, "parameters": {"distributed_campaign": {
        "enabled": True, "generator_id": "rfd4_proteina" if workflow.endswith("rfd4_proteina") else "bindcraft2",
        "gpu_selection": "auto", "workstation_selection": "eligible",
    }}}
    runner._submit(request, parent, {"id": "tile-1"}, owner_id="owner", project_id="project", gpu_device="all")

    assert len(submitted) == 1
    assert submission_context == [{"studio_owner_id": "owner", "studio_project_id": "project", "studio_development_mode": True}]
    assert submitted[0]["parameters"]["distributed_campaign"]["enabled"] is True
    assert submitted[0]["parameters"]["distributed_campaign"]["workstation_selection"] == "eligible"
    assert submitted[0]["resources"]["gpu_device"] == "none"
    assert submitted[0]["resources"]["min_vram_gb"] == 0


def test_all_generators_template_has_two_sequence_stages_and_boltz2_analysis():
    g = graph.new_graph("all_generators_boltz2")
    plan = graph.validate_graph(g)
    kinds = [n["data"]["kind"] for n in g["nodes"]]
    generators = {"rfd3", "rfdiffusion", "genie3", "pxdesign", "protpardelle", "bindcraft2_generator", "boltzgen_generator", "complexa_generator", "esm_binder"}
    assert generators <= set(kinds)
    assert kinds.count("ligandmpnn") == 2
    assert kinds.count("boltz2") == kinds.count("ipsae") == kinds.count("rosetta_metrics") == 1
    by_id = {n["id"]: n["data"]["kind"] for n in g["nodes"]}
    incoming = plan["incoming"]
    target_id = next(n["id"] for n in g["nodes"] if n["data"]["kind"] == "target")
    context_id = next(n["id"] for n in g["nodes"] if n["data"]["kind"] == "binder_context")
    mpnn_id = next(n["id"] for n in g["nodes"] if n["data"]["kind"] == "ligandmpnn")
    assert len(incoming[mpnn_id]) == len(generators)
    assert incoming[context_id] == [target_id]
    assert all(incoming[node_id] == [context_id] for node_id, kind in by_id.items() if kind in generators)
    downstream = [by_id[node_id] for node_id in plan["order"] if by_id[node_id] in {"ligandmpnn", "boltz2", "ipsae", "rosetta_metrics", "results"}]
    assert downstream == ["ligandmpnn", "ligandmpnn", "boltz2", "ipsae", "rosetta_metrics", "results"]
    boltz = next(n for n in g["nodes"] if n["data"]["kind"] == "boltz2")
    assert boltz["data"]["params"]["write_full_pae"] is True
    mpnn = next(n for n in g["nodes"] if n["data"]["kind"] == "ligandmpnn")
    assert mpnn["data"]["params"]["fix_interface_residues"] is False


@pytest.mark.parametrize("case",["cycle","wrong_type","duplicate_id","two_inputs","missing_input","bad_number"])
def test_graph_rejects_invalid_dependencies_before_execution(case):
    g=graph.new_graph()
    by_kind={n["data"]["kind"]:n["id"] for n in g["nodes"]}
    if case=="cycle":
        # Sequence nodes can feed each other, but the cycle must still fail.
        g["nodes"]=[node("ligandmpnn","a"),node("ligandmpnn","b"),node("results","r")]
        g["edges"]=[edge("a","b"),edge("b","a"),edge("b","r")]
    elif case=="wrong_type":
        g["edges"][0]["target"]=by_kind["rf3"]
    elif case=="duplicate_id":
        g["nodes"][1]["id"]=g["nodes"][0]["id"]
    elif case=="two_inputs":
        g["edges"].append(edge(by_kind["target"],by_kind["rfd3"]))
    elif case=="missing_input":
        g["edges"].pop(1)
    else:
        next(n for n in g["nodes"] if n["data"]["kind"]=="rfd3")["data"]["params"]["num_designs"]=-1
    with pytest.raises(ValueError):graph.validate_graph(g)


def test_filter_keeps_upstream_type_and_format_constraints():
    g=graph.new_graph()
    by_kind={n["data"]["kind"]:n["id"] for n in g["nodes"]}
    g["nodes"].append(node("filter","f"))
    g["edges"]=[edge(by_kind["target"],by_kind["binder_context"]),
                 edge(by_kind["binder_context"],by_kind["rfd3"]),
                 edge(by_kind["rfd3"],"f"),edge("f",by_kind["foundry_mpnn"]),
                 edge(by_kind["foundry_mpnn"],by_kind["rf3"]),
                 edge(by_kind["rf3"],by_kind["results"])]
    graph.validate_graph(g)
    index=next(i for i,n in enumerate(g["nodes"]) if n["data"]["kind"]=="rfd3")
    g["nodes"][index]=node("rfdiffusion",by_kind["rfd3"])
    # The common MPNN adapter can use a shared PDB backend through a filter too.
    graph.validate_graph(g)
    g["nodes"][index]=node("rfd3",by_kind["rfd3"])
    g["edges"]=[edge(by_kind["target"],by_kind["binder_context"]),
                 edge(by_kind["binder_context"],by_kind["rfd3"]),edge(by_kind["rfd3"],"f"),
                 edge("f",by_kind["rf3"]),edge(by_kind["rf3"],by_kind["results"])]
    g["nodes"]=[n for n in g["nodes"] if n["data"]["kind"]!="foundry_mpnn"]
    with pytest.raises(ValueError,match="required"):graph.validate_graph(g)


def test_versions_remain_immutable_while_draft_and_positions_change(monkeypatch,tmp_path):
    monkeypatch.setattr(graph,"app_home",lambda:tmp_path)
    g=graph.new_graph();first=graph.save_definition(g,"owner","project",version=True)
    changed=copy.deepcopy(g);changed["nodes"][1]["position"]["x"]+=200
    second=graph.save_definition(changed,"owner","project",version=True)
    assert first["revision"]!=second["revision"]
    versions=graph.list_versions("owner","project",g["id"])
    assert len(versions)==2
    assert next(v for v in versions if v["revision"]==first["revision"])["graph"]==g
    assert graph.list_definitions("owner","other-project")==[]
    assert graph.list_definitions("other-owner","project")==[]


def test_adapters_disable_hidden_downstream_pipelines_and_bind_existing_signatures(design_source_factory):
    import inspect
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS,_get_function
    from mn_molecular_design_studio.rfdiffusion3_input import ASSETS
    source_run, _ = design_source_factory()
    source={"run_dir":str(source_run),"candidates_jsonl":str(source_run/"artifacts/normalized_candidates/candidates.jsonl"),"target":{"target_pdb":str(ASSETS/"5o45_cropped.pdb"),"target_chains":["A"],"hotspots":"A54"}}
    for kind,spec in graph.NODE_CATALOG.items():
        if not spec["workflow"]:continue
        p={f["name"]:f["default"] for f in spec["fields"]}
        request=runner.node_request(kind,p,source,"test")
        inspect.signature(_get_function(WORKFLOWS[request["workflow"]])).bind(**request["parameters"])
        if kind in {"rfdiffusion","genie3","pxdesign","protpardelle"}:
            assert request["workflow"]=="design."+kind+"_generate"
            assert "run_vanilla_pipeline" not in request["parameters"]
            assert "run_mode" not in request["parameters"]
        if kind=="rfd3":assert request["parameters"]["low_memory_mode"] is False


def test_rfd4_campaign_graph_accepts_native_condition_spec_input():
    campaign = {
        "enabled": True, "objective_type": "generated", "generated_target": 1,
        "max_attempts": 1, "batch_size": 1, "max_concurrent_batches": 1,
        "max_concurrent_generation_jobs": 1, "workstation_selection": "eligible",
        "gpu_selection": "all",
    }
    target_path = "/tmp/rfd4-idr-target.cif"
    condition_spec = "idr_peptide_binder:\n  input: rfd4-idr-target.cif\n  contig: '70-110, /0, A1-70'\n  conditions:\n    C_SEQ:\n      true:\n        - select: \"sel('A')\"\n"
    request = runner.node_request("rfd4", {
        "task": "protein_binder", "input_mode": "condition_spec",
        "condition_spec": condition_spec, "target_structure": target_path,
        "num_designs": 1, "input_files": "", "options_json": "{}",
        "distributed_campaign": campaign,
    }, {"target": {"target_pdb": target_path, "target_chains": ["A"], "hotspots": ""}}, "IDR campaign")
    assert request["workflow"] == "design.rfd4_proteina"
    assert request["parameters"]["input_mode"] == "condition_spec"
    assert request["parameters"]["condition_spec"] == condition_spec
    assert request["parameters"]["distributed_campaign"]["generated_target"] == 1


def test_existing_design_inputs_require_successful_filter_passes_in_same_project(design_source_factory):
    source,_=design_source_factory()
    assert len(runner.source_designs(str(source),"test-project")[1])==1
    for project,selected in [("other-project",None),("test-project",[]),("test-project",["unknown"])]:
        with pytest.raises(ValueError):runner.source_designs(str(source),project,selected)
    result=json.loads((source/"result.json").read_text());result["success"]=False;(source/"result.json").write_text(json.dumps(result))
    with pytest.raises(ValueError):runner.source_designs(str(source),"test-project")


@pytest.fixture
def execution_harness(monkeypatch,tmp_path,design_source_factory):
    from mn_molecular_design_studio import workspace
    monkeypatch.setattr(workspace,"mark_studio_job",lambda *a,**k:None)
    monkeypatch.setattr(graph,"app_home",lambda:tmp_path/"home")
    source,base=design_source_factory()
    g=graph.new_graph();g["nodes"][0]["data"]["params"].update(target_pdb=str(source/"designed.cif"),target_chains="A")
    calls=[];parents=[]
    def create(*a,**k):
        p=tmp_path/f"run-{len(parents)}";p.mkdir();parents.append(p)
        (p/"metadata.json").write_text(json.dumps({"status":"queued","job_type":"visual_pipeline"}))
        (p/"input.json").write_text(json.dumps({"tool":"visual_pipeline","params":{"pipeline_name":"test"}}))
        return SimpleNamespace(run_dir=p)
    monkeypatch.setattr(runner,"create_job",create)
    def submit(request,parent,n,**kwargs):
        calls.append(n["data"]["kind"])
        p=tmp_path/f"child-{len(calls)}";p.mkdir()
        c=copy.deepcopy(base);c["candidate_id"]="design";c["complex_pdb"]=str(source/"designed.cif")
        c["stage"]={"rfd3":"generation.backbone","foundry_mpnn":"generation.backbone_sequence"}.get(n["data"]["kind"],"complex_refolding")
        c["metrics"]={"ranking_score":.8}
        normalized=p/"artifacts/normalized_candidates";normalized.mkdir(parents=True)
        (normalized/"candidates.jsonl").write_text(json.dumps(c)+"\n")
        (p/"metadata.json").write_text(json.dumps({"status":"completed","job_code":str(len(calls)),"pipeline_parent_run_dir":str(parent)}))
        (p/"result.json").write_text(json.dumps({"success":True}))
        return p
    monkeypatch.setattr(runner,"_submit",submit)
    return g,calls,parents


def test_execution_reuses_only_unchanged_steps_and_preserves_prior_results(execution_harness):
    g,calls,parents=execution_harness
    p1=runner.run_graph_pipeline(g,"first","owner","project")
    assert calls==["rfd3","foundry_mpnn","rf3"]
    original=(p1/"result.json").read_bytes();snapshot=(p1/"graph_snapshot.json").read_bytes()
    moved=copy.deepcopy(g);next(n for n in moved["nodes"] if n["data"]["kind"]=="binder_context")["position"]["y"]+=400
    p2=runner.run_graph_pipeline(moved,"moved","owner","project")
    assert len(calls)==3
    assert json.loads((p2/"result.json").read_text())["metrics"]["cached_steps"]==4
    changed=copy.deepcopy(g);next(n for n in changed["nodes"] if n["data"]["kind"]=="foundry_mpnn")["data"]["params"]["sequences"]+=1
    runner.run_graph_pipeline(changed,"changed","owner","project")
    assert calls==["rfd3","foundry_mpnn","rf3","foundry_mpnn","rf3"]
    assert (p1/"result.json").read_bytes()==original and (p1/"graph_snapshot.json").read_bytes()==snapshot


def test_parallel_branches_collect_distinct_candidate_ids(execution_harness):
    g,calls,_=execution_harness
    by_kind={n["data"]["kind"]:n["id"] for n in g["nodes"]}
    g["nodes"].append(node("esmfold","esm"));g["edges"] += [edge(by_kind["foundry_mpnn"],"esm"),edge("esm",by_kind["results"])]
    p=runner.run_graph_pipeline(g,"branched","owner","project")
    assert set(calls)=={"rfd3","foundry_mpnn","rf3","esmfold"}
    candidates=json.loads((p/"result.json").read_text())["outputs"]["candidates"]
    assert len(candidates)==2 and len({c["candidate_id"] for c in candidates})==2


@pytest.mark.parametrize("analysis_success", [True, False])
def test_surface_analysis_branch_is_collected_as_a_pipeline_result(execution_harness, monkeypatch, tmp_path, analysis_success):
    g, _calls, _parents = execution_harness
    target = next(n for n in g["nodes"] if n["data"]["kind"] == "target")
    target["data"]["params"].update(target_chains="A")
    analysis = node("surface_analysis", "surface")
    results = node("results", "results")
    g["nodes"] = [target, analysis, results]
    g["edges"] = [edge("target-0", "surface"), edge("surface", "results")]
    child = tmp_path / "surface-child"
    child.mkdir()
    (child / "metadata.json").write_text(json.dumps({"status": "completed" if analysis_success else "failed", "run_id": child.name, "task_group": "detection", "job_code": "SA-1"}))
    statuses = {tool: {"success": analysis_success} for tool in ("pesto", "scannet", "surf2spot", "masif_seed")}
    (child / "result.json").write_text(json.dumps({"success": analysis_success, "outputs": {"tool_status": statuses}}))
    captured = {}

    def submit(request, _parent, _node, **_kwargs):
        captured.update(request)
        return child

    monkeypatch.setattr(runner, "_submit", submit)
    parent = runner.run_graph_pipeline(g, "surface-test", "owner", "project")

    result = json.loads((parent / "result.json").read_text())
    execution = json.loads((parent / "graph_execution.json").read_text())
    assert captured["workflow"] == "detection.surface_analysis"
    assert result["success"] is analysis_success
    assert len(result["outputs"]["surface_analyses"]) == 1
    assert execution["nodes"]["surface"]["output"]["analysis_run"]["job_code"] == "SA-1"
    assert execution["nodes"]["results"]["status"] == "completed"


def test_surface_cache_requires_predictions_and_valid_result_references(execution_harness, monkeypatch, tmp_path):
    g, _calls, _parents = execution_harness
    target = next(n for n in g["nodes"] if n["data"]["kind"] == "target")
    target["data"]["params"].update(target_chains="A")
    g["nodes"] = [target, node("surface_analysis", "surface"), node("results", "results")]
    g["edges"] = [edge("target-0", "surface"), edge("surface", "results")]
    submissions = []

    def submit(*_args, **_kwargs):
        child = tmp_path / f"surface-child-{len(submissions)}"
        child.mkdir()
        prediction = child / "prediction.csv"
        prediction.write_text("chain,residue,score\nA,1,0.9\n")
        (child / "metadata.json").write_text(json.dumps({"status": "completed", "run_id": child.name, "task_group": "detection"}))
        (child / "result.json").write_text(json.dumps({"success": True, "outputs": {
            "tool_status": {"pesto": {"success": True}}, "artifacts": [{"path": str(prediction)}],
        }}))
        submissions.append(child)
        return child

    monkeypatch.setattr(runner, "_submit", submit)
    runner.run_graph_pipeline(g, "first", "owner", "project")
    cached = runner.run_graph_pipeline(g, "cached", "owner", "project")
    assert len(submissions) == 1
    assert json.loads((cached / "result.json").read_text())["metrics"]["cached_steps"] == 2

    (submissions[0] / "prediction.csv").unlink()
    rerun = runner.run_graph_pipeline(g, "missing-prediction", "owner", "project")
    assert len(submissions) == 2
    result = json.loads((rerun / "result.json").read_text())
    assert result["metrics"]["surface_analysis_count"] == 1
    assert result["outputs"]["surface_analyses"][0]["run_id"] == submissions[1].name


def test_completed_surface_job_without_method_results_fails_instead_of_waiting(execution_harness, monkeypatch, tmp_path):
    g, _calls, _parents = execution_harness
    target = next(n for n in g["nodes"] if n["data"]["kind"] == "target")
    target["data"]["params"].update(target_chains="A")
    g["nodes"] = [target, node("surface_analysis", "surface"), node("results", "results")]
    g["edges"] = [edge("target-0", "surface"), edge("surface", "results")]
    child = tmp_path / "surface-child"
    child.mkdir()
    (child / "metadata.json").write_text('{"status":"completed"}')
    (child / "result.json").write_text('{"success":true,"outputs":{}}')
    monkeypatch.setattr(runner, "_submit", lambda *_args, **_kwargs: child)

    parent = runner.run_graph_pipeline(g, "missing-results", "owner", "project")

    execution = json.loads((parent / "graph_execution.json").read_text())
    assert execution["nodes"]["surface"]["status"] == "failed"
    assert "method results" in execution["nodes"]["surface"]["error"]
    assert execution["nodes"]["results"]["status"] == "blocked"


def test_generator_fan_in_merges_every_candidate_before_shared_mpnn(tmp_path, monkeypatch, design_source_factory):
    from mn_molecular_design_studio import workspace
    from types import SimpleNamespace
    source, base = design_source_factory()
    parent = tmp_path / "pipeline-run"
    parent.mkdir()
    (parent / "metadata.json").write_text(json.dumps({"status": "queued", "job_type": "visual_pipeline"}))
    (parent / "input.json").write_text(json.dumps({"tool": "visual_pipeline", "params": {"pipeline_name": "fan-in-test"}}))
    monkeypatch.setattr(graph, "app_home", lambda: tmp_path / "home")
    monkeypatch.setattr(runner, "create_job", lambda *a, **k: SimpleNamespace(run_dir=parent))
    monkeypatch.setattr(workspace, "mark_studio_job", lambda *a, **k: None)
    g = graph.new_graph("blank")
    g["nodes"] = [node(kind, identifier) for kind, identifier in [
        ("target", "target"), ("genie3", "gen-a"), ("genie3", "gen-b"),
        ("ligandmpnn", "mpnn"), ("boltz2", "boltz"), ("ipsae", "ipsae"),
        ("rosetta_metrics", "rosetta"), ("results", "results"),
    ]]
    target = next(n for n in g["nodes"] if n["data"]["kind"] == "target")
    target["data"]["params"].update(target_pdb=str(source / "designed.cif"), target_chains="A")
    for n in g["nodes"]:
        if n["data"]["kind"] == "boltz2":
            n["data"]["params"]["write_full_pae"] = True
    g["edges"] = [edge("target", "gen-a"), edge("target", "gen-b"), edge("gen-a", "mpnn"), edge("gen-b", "mpnn"),
                  edge("mpnn", "boltz"), edge("boltz", "ipsae"), edge("ipsae", "rosetta"), edge("rosetta", "results")]
    captured = {}
    child_index = 0

    def submit(request, parent, current_node, **kwargs):
        nonlocal child_index
        kind = current_node["data"]["kind"]
        if kind == "ligandmpnn":
            candidate_path = Path(request["parameters"]["candidates_jsonl"])
            rows = [json.loads(line) for line in candidate_path.read_text().splitlines() if line.strip()]
            captured["fan_in"] = rows
        child_index += 1
        child = tmp_path / f"child-{child_index}"
        normalized = child / "artifacts/normalized_candidates"
        normalized.mkdir(parents=True)
        row = copy.deepcopy(base)
        row.update(candidate_id=f"{kind}-{child_index}", complex_pdb=str(source / "designed.cif"),
                   stage="generation.backbone" if kind == "genie3" else "complex_refolding")
        (normalized / "candidates.jsonl").write_text(json.dumps(row) + "\n")
        (child / "metadata.json").write_text(json.dumps({"status": "completed", "job_code": str(child_index), "pipeline_parent_run_dir": str(parent)}))
        (child / "result.json").write_text(json.dumps({"success": True}))
        (child / "input.json").write_text(json.dumps({"params": {}}))
        return child

    monkeypatch.setattr(runner, "_submit", submit)
    runner.run_graph_pipeline(g, "fan-in-test", "owner", "project")
    assert len(captured["fan_in"]) == 2
    assert {row["raw_metadata"]["pipeline_generator_node"] for row in captured["fan_in"]} == {"gen-a", "gen-b"}


def test_refolded_result_includes_target_aligned_binder_rmsd(tmp_path, design_source_factory):
    source, candidate = design_source_factory()
    prediction = copy.deepcopy(candidate)
    prediction.update(candidate_id="refolded", complex_pdb=str(source / "designed.cif"), parents=[candidate["candidate_id"]])
    output = runner._local_source(tmp_path, "predicted", [prediction])
    source_output = {"run_dir": str(source), "candidates_jsonl": str(source / "artifacts/normalized_candidates/candidates.jsonl")}
    tagged = runner._attach_refold_input_references(output, source_output, tmp_path, "boltz")
    row = runner._absolute_candidates(tagged)[0]
    scored = runner._add_target_aligned_pose_rmsd(row, tmp_path)
    assert scored["metrics"]["target_aligned_binder_rmsd"] == pytest.approx(0.0, abs=1e-5)
    assert scored["metrics"]["target_aligned_binder_rmsd_status"] == "ok"


def test_failure_blocks_descendants_without_submitting_them(execution_harness,monkeypatch):
    g,calls,_=execution_harness
    def fail(*a,**k):raise ValueError("test engine failure")
    monkeypatch.setattr(runner,"_submit",fail)
    p=runner.run_graph_pipeline(g,"failure","owner","project")
    assert calls==[]
    result=json.loads((p/"result.json").read_text())
    assert result["success"] is False
    foundry_mpnn_id=next(n["id"] for n in g["nodes"] if n["data"]["kind"]=="foundry_mpnn")
    assert result["outputs"]["pipeline_nodes"][foundry_mpnn_id]["status"]=="blocked"


def test_foundry_mpnn_locks_target_and_selects_binder_chain(monkeypatch,tmp_path,design_source_factory):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import sequence_design as sequence
    source,candidate=design_source_factory();candidate["stage"]="generation.backbone"
    candidate["complex_pdb"]=str(source/"designed.cif")
    monkeypatch.setattr(sequence,"read_candidates",lambda *a:[candidate])
    run=tmp_path/"mpnn";run.mkdir()
    monkeypatch.setattr(sequence,"create_job",lambda *a,**k:SimpleNamespace(run_dir=run))
    commands=[]
    monkeypatch.setattr(sequence,"_run_shell_steps",lambda run,steps:commands.extend(s["command"] for s in steps) or 1)
    monkeypatch.setattr(sequence,"finish_job",lambda *a:None)
    sequence.run_foundry_mpnn_sequence_design(source,source/"artifacts/normalized_candidates/candidates.jsonl",model_type="protein_mpnn")
    command=commands[0]
    assert command[command.index("--designed_chains")+1]=="B"


def test_parent_is_queued_once_as_a_cpu_coordinator(monkeypatch,tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.core import workflow_queue
    calls=[]
    monkeypatch.setattr(workflow_queue,"enqueue_workflow_call",lambda *a,**k:calls.append((a,k)) or tmp_path)
    g=graph.new_graph();g["nodes"][0]["data"]["params"]["target_pdb"]="/tmp/target.pdb"
    result=runner.enqueue_graph_pipeline(g,pipeline_name="test",owner_id="owner",project_id="project")
    assert result==tmp_path and len(calls)==1
    assert calls[0][1]["gpu_device"]=="none" and calls[0][1]["job_type"]=="visual_pipeline"
    assert calls[0][1]["kwargs"]["gpu_device"]=="all"


def test_api_saves_and_validates_without_queuing_jobs(monkeypatch,tmp_path):
    from fastapi.testclient import TestClient
    from mn_molecular_design_studio import api
    from mn_molecular_design_studio.rfdiffusion3_input import ASSETS
    from shutil import copy2
    monkeypatch.setenv("MDS_CLUSTER_HOME", str(tmp_path / "cluster"))
    actor={"id":"test-user","role":"user"}
    api.app.dependency_overrides[api._current_user]=lambda:actor
    monkeypatch.setattr(graph,"app_home",lambda:tmp_path/"home")
    monkeypatch.setattr(api,"_request_payload_and_project",lambda payload,user:({k:v for k,v in payload.items() if k!="project_id"},payload["project_id"]))
    uploaded=tmp_path/"inputs";uploaded.mkdir();copy2(ASSETS/"5o45_cropped.pdb",uploaded/"target.pdb")
    monkeypatch.setattr(api,"user_input_root",lambda _:uploaded)
    monkeypatch.setattr(api,"reference_root",lambda:tmp_path/"refs")
    monkeypatch.setattr(api,"path_is_within_accessible_job",lambda *a,**k:False)
    monkeypatch.setattr(api,"accessible_project_ids",lambda _:[])
    calls=[]
    monkeypatch.setattr(runner,"enqueue_graph_pipeline",lambda *a,**k:calls.append((a,k)))
    g=graph.new_graph();g["nodes"][0]["data"]["params"]["target_pdb"]=str(uploaded/"target.pdb")
    client=TestClient(api.app)
    try:
        payload={"project_id":"project","graph":g}
        assert client.post("/api/v1/pipelines/save",json=payload).status_code==200
        assert client.post("/api/v1/pipelines/validate",json=payload).status_code==200
        assert calls==[]

        surface_graph = graph.new_graph("blank")
        surface_graph["nodes"] = [node("target", "target"), node("surface_analysis", "surface"), node("results", "results")]
        surface_graph["nodes"][0]["data"]["params"].update(target_pdb=str(uploaded / "target.pdb"), target_chains="A")
        surface_graph["nodes"][1]["data"]["params"]["masif_chain_id"] = "A"
        surface_graph["edges"] = [edge("target", "surface"), edge("surface", "results")]
        surface_payload = {"project_id": "project", "graph": surface_graph}
        assert client.post("/api/v1/pipelines/validate", json=surface_payload).status_code == 200
        surface_graph["nodes"][1]["data"]["params"]["masif_chain_id"] = "not-a-target-chain"
        invalid_surface = client.post("/api/v1/pipelines/validate", json=surface_payload)
        assert invalid_surface.status_code == 422 and "MaSIF-seed chain" in invalid_surface.text

        split_analysis_graph = graph.new_graph("blank")
        split_analysis_graph["nodes"] = [node(kind, identifier) for kind, identifier in [
            ("target", "target"), ("rfd4", "design"), ("boltz2", "cofold"),
            ("ipsae", "analysis"), ("results", "results"),
        ]]
        next(item for item in split_analysis_graph["nodes"] if item["id"] == "target")["data"]["params"].update(
            target_pdb=str(uploaded / "target.pdb"), target_chains="A")
        split_analysis_graph["edges"] = [edge("target", "design"), edge("design", "cofold"),
                                          edge("design", "analysis"), edge("cofold", "analysis"),
                                          edge("analysis", "results")]
        split_response = client.post("/api/v1/pipelines/validate", json={
            "project_id": "project", "graph": split_analysis_graph,
        })
        assert split_response.status_code == 200, split_response.text
        assert calls == []

        bad=copy.deepcopy(payload);bad["graph"]["nodes"][0]["data"]["params"]["target_pdb"]="/etc/passwd"
        assert client.post("/api/v1/pipelines/save",json=bad).status_code==403
        missing=copy.deepcopy(payload);missing["graph"]["nodes"][0]["data"]["params"]["target_pdb"]=""
        assert client.post("/api/v1/pipelines/run",json=missing).status_code==422
        assert calls==[]
        assert client.post("/api/v1/pipelines/save",json=payload|{"owner_id":"another-user"}).status_code==422
    finally:
        api.app.dependency_overrides.clear()


def test_cancel_only_stops_active_jobs_owned_by_this_pipeline(monkeypatch,tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs
    parent=tmp_path/"pipeline";parent.mkdir()
    nodes={}
    for name,status,owner in [("running","running",str(parent)),("done","completed",str(parent)),("unrelated","running","another-pipeline")]:
        child=tmp_path/name;child.mkdir();(child/"metadata.json").write_text(json.dumps({"status":status,"pipeline_parent_run_dir":owner}));nodes[name]={"job_run_dir":str(child)}
    (parent/"graph_execution.json").write_text(json.dumps({"nodes":nodes}))
    stopped=[];monkeypatch.setattr(jobs,"stop_job",lambda p,**k:stopped.append(p.name))
    runner.cancel_graph_children(parent)
    assert stopped==["running"]


def test_canvas_host_renders_without_running_a_job(monkeypatch):
    from streamlit.testing.v1 import AppTest
    app=AppTest.from_string('''
from mn_molecular_design_studio.pipeline_graph_ui import render_pipeline_builder
render_pipeline_builder(owner_id="test-user",project_selector=lambda _:None,
    submit_api=lambda *a:(_ for _ in ()).throw(AssertionError("No request expected")),
    search_jobs=lambda **k:{"items":[]})
''').run()
    assert not app.exception
    assert any("Choose a project" in c.value for c in app.info)


def test_canvas_edits_preserve_trusted_parameters_for_every_tool_kind():
    from mn_molecular_design_studio import candidate_pipeline, pipeline_graph_ui

    cases = [
        (candidate_pipeline.new_pipeline(['rfd3'], mode='modular_design'), 'rfd3', 'seed', 11, 99),
        (candidate_pipeline.new_pipeline(['boltzgen'], mode='native_complete'), 'boltzgen', 'num_designs', 3, 77),
    ]
    for graph_data, kind, field_name, trusted_value, component_value in cases:
        target = next(node for node in graph_data['nodes'] if node['data']['kind'] == 'target')
        target['data']['params']['target_pdb'] = '/tmp/trusted-pdl1-target.pdb'
        design = next(node for node in graph_data['nodes'] if node['data']['kind'] == kind)
        design['data']['params'][field_name] = trusted_value
        component_graph = copy.deepcopy(graph_data)
        next(node for node in component_graph['nodes'] if node['data']['kind'] == 'target')['data']['params']['target_pdb'] = '/etc/passwd'
        next(node for node in component_graph['nodes'] if node['data']['kind'] == kind)['data']['params'][field_name] = component_value
        restored = pipeline_graph_ui._restore_component_parameters(component_graph, graph_data)
        assert next(node for node in restored['nodes'] if node['data']['kind'] == 'target')['data']['params']['target_pdb'] == '/tmp/trusted-pdl1-target.pdb'
        assert next(node for node in restored['nodes'] if node['data']['kind'] == kind)['data']['params'][field_name] == trusted_value


def test_api_address_persists_when_saving_other_workspace_settings(tmp_path):
    from mn_molecular_design_studio.workspace import load_settings, save_settings
    settings={"results_directory":str(tmp_path/"results"),"reference_directory":str(tmp_path/"references"),"api_base_url":"http://127.0.0.1:8767/"}
    save_settings(settings,home=tmp_path)
    assert load_settings(home=tmp_path)["api_base_url"]=="http://127.0.0.1:8767"
    save_settings({"results_directory":str(tmp_path/"results"),"reference_directory":str(tmp_path/"references")},home=tmp_path)
    assert load_settings(home=tmp_path)["api_base_url"]=="http://127.0.0.1:8767"
    with pytest.raises(ValueError,match="HTTP API address"):
        save_settings({**settings,"api_base_url":"https://user:password@example.org/path"},home=tmp_path)


def test_canvas_restores_and_opens_saved_graph_after_name_widget(monkeypatch,tmp_path):
    from streamlit.testing.v1 import AppTest
    monkeypatch.setattr(graph,"app_home",lambda:tmp_path)
    original=graph.new_graph();original["name"]="Saved pipeline"
    graph.save_definition(original,"test-user","test-project",version=True)
    app=AppTest.from_string('''
from mn_molecular_design_studio.pipeline_graph_ui import render_pipeline_builder
render_pipeline_builder(owner_id="test-user",project_selector=lambda _:"test-project",
    submit_api=lambda *a:{},search_jobs=lambda **k:{"items":[]})
''').run()
    assert not app.exception
    app.text_input(key="visual-pipeline-test-project-name").set_value("Changed draft").run()
    app.button(key="visual-pipeline-test-project-restore").click().run()
    assert not app.exception
    assert app.text_input(key="visual-pipeline-test-project-name").value=="Saved pipeline"
    app.text_input(key="visual-pipeline-test-project-name").set_value("Changed again").run()
    app.button(key="visual-pipeline-test-project-open").click().run()
    assert not app.exception
    assert app.text_input(key="visual-pipeline-test-project-name").value=="Saved pipeline"


@pytest.mark.parametrize('failure_policy', ['all', 'partial'])
def test_explicit_partial_fan_in_reports_missing_branch_without_changing_saved_default(execution_harness, monkeypatch, failure_policy):
    base, calls, _ = execution_harness
    target = copy.deepcopy(next(n for n in base['nodes'] if n['data']['kind'] == 'target'))
    target['id'] = 'target'
    graph_data = graph.new_graph('blank')
    result_node = node('results', 'results')
    if failure_policy == 'partial':
        result_node['data']['execution'] = {'failure_policy': 'partial'}
    graph_data['nodes'] = [target, node('rfd3', 'good'), node('rfd3', 'bad'), result_node]
    graph_data['edges'] = [edge('target', 'good'), edge('target', 'bad'), edge('good', 'results'), edge('bad', 'results')]
    original_submit = runner._submit
    def submit(request, parent, current_node, **kwargs):
        if current_node['id'] == 'bad':
            raise RuntimeError('Independent generator failed')
        return original_submit(request, parent, current_node, **kwargs)
    monkeypatch.setattr(runner, '_submit', submit)
    run = runner.run_graph_pipeline(graph_data, 'partial-fan-in', 'owner', 'project', reuse_completed=False)
    payload = json.loads((run / 'result.json').read_text())
    assert payload['success'] is (failure_policy == 'partial')
    states = payload['outputs']['pipeline_nodes']
    if failure_policy == 'partial':
        assert states['results']['missing_branches'] == ['bad']
        assert payload['outputs']['comparison_status'] == 'partial'
        assert payload['outputs']['candidates'][0]['raw_metadata']['partial_pool']['missing_branches'] == ['bad']
    else:
        assert states['results']['status'] == 'blocked'


def test_timed_pipeline_objective_requires_no_filter_node_and_round_trips():
    target, generator, results = node('target', 'target'), node('rfd4', 'generator'), node('results', 'results')
    definition = {'schema_version': graph.SCHEMA, 'id': 'e' * 32,
                  'nodes': [target, generator, results],
                  'edges': [edge('target', 'generator'), edge('generator', 'results')],
                  'shared_campaign': {'enabled': True, 'objective_type': 'time',
                                      'generator_nodes': ['generator'], 'max_attempts': 100,
                                      'max_runtime_seconds': 1800}}
    saved = graph.validate_graph(definition)['graph']
    scope = saved['shared_campaign']
    assert scope['validation_stages'] == []
    assert scope['accepted_target'] is None and scope['generated_target'] is None
    assert 'acceptance_node' not in scope
    assert graph.validate_graph(saved)['graph']['shared_campaign'] == scope


def test_passing_objective_references_optional_format_aware_filter_node():
    from mn_molecular_design_studio.validation_acceptance import METRICS
    target, generator = node('target', 'target'), node('rfd4', 'generator')
    merge, accept, results = node('candidate_merge', 'collection'), node('shared_filter', 'acceptance'), node('results', 'results')
    accept['data']['params'].update(binder_format='peptide', target_type='fasta_idr',
                                   rules=[{'field': METRICS['iptm']['field'], 'operator': '>=', 'value': .71}])
    definition = {'schema_version': graph.SCHEMA, 'id': 'f' * 32,
                  'nodes': [target, generator, merge, accept, results],
                  'edges': [edge('target', 'generator'), edge('generator', 'collection'),
                            edge('collection', 'acceptance'), edge('acceptance', 'results')],
                  'shared_campaign': {'enabled': True, 'objective_type': 'accepted',
                                      'generator_nodes': ['generator'], 'accepted_target': 100, 'max_attempts': 300,
                                      'acceptance_node': 'acceptance',
                                      'validation_stages': [{'workflow': 'refolding.optimized_prediction',
                                                             'parameters': {'engine': 'boltz2'}}]}}
    saved = graph.validate_graph(definition)['graph']
    stage = saved['shared_campaign']['validation_stages'][0]
    assert stage['acceptance_node'] == 'acceptance'
    assert stage['acceptance_policy']['mandatory_filters'][0]['value'] == .71
    assert stage['acceptance_policy']['binder_format'] == 'peptide'
    assert stage['acceptance_policy']['target_type'] == 'fasta_idr'
    assert graph.validate_graph(saved)['graph']['shared_campaign'] == saved['shared_campaign']
    definition['shared_campaign'] = {'enabled': True, 'objective_type': 'generated',
                                     'generator_nodes': ['generator'], 'generated_target': 10, 'max_attempts': 10}
    scope = graph.validate_graph(definition)['graph']['shared_campaign']
    assert 'acceptance_node' not in scope
    assert scope['validation_stages'] == []  # Optional filter remains an ordinary pipeline node.
