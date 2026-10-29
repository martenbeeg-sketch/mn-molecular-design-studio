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


def test_graph_runner_injects_one_shared_campaign_context_into_normal_branch_requests(tmp_path, monkeypatch, design_source_factory):
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
    assert all(request["parameters"]["distributed_campaign"]["accepted_target"] == 1 for request in requests)
    assert all(request["parameters"]["distributed_campaign"]["workstation_selection"] == "auto" for request in requests)
    shared_path = result / "artifacts/shared_campaign/ledger.json"
    assert shared_path.is_file()
    state = json.loads(shared_path.read_text())
    assert state["config"]["generator_nodes"] == ["bindcraft", "rfd4"]


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
