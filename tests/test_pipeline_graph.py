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
    assert [n["data"]["kind"] for n in plan["graph"]["nodes"]]==["target","rfd3","foundry_mpnn","rf3","results"]
    assert len(plan["order"])==5


def test_all_generators_template_fans_into_one_mpnn_and_runs_boltz2_analysis():
    g = graph.new_graph("all_generators_boltz2")
    plan = graph.validate_graph(g)
    kinds = [n["data"]["kind"] for n in g["nodes"]]
    generators = {"rfd3", "rfdiffusion", "genie3", "pxdesign", "protpardelle", "bindcraft2_generator", "boltzgen_generator", "complexa_generator", "esm_binder"}
    assert generators <= set(kinds)
    assert kinds.count("ligandmpnn") == kinds.count("boltz2") == kinds.count("ipsae") == kinds.count("rosetta_metrics") == 1
    by_id = {n["id"]: n["data"]["kind"] for n in g["nodes"]}
    incoming = plan["incoming"]
    target_id = next(n["id"] for n in g["nodes"] if n["data"]["kind"] == "target")
    mpnn_id = next(n["id"] for n in g["nodes"] if n["data"]["kind"] == "ligandmpnn")
    assert len(incoming[mpnn_id]) == len(generators)
    assert all(target_id in incoming[node_id] for node_id, kind in by_id.items() if kind in generators)
    downstream = [by_id[node_id] for node_id in plan["order"] if by_id[node_id] in {"ligandmpnn", "boltz2", "ipsae", "rosetta_metrics", "results"}]
    assert downstream == ["ligandmpnn", "boltz2", "ipsae", "rosetta_metrics", "results"]
    boltz = next(n for n in g["nodes"] if n["data"]["kind"] == "boltz2")
    assert boltz["data"]["params"]["write_full_pae"] is True
    mpnn = next(n for n in g["nodes"] if n["data"]["kind"] == "ligandmpnn")
    assert mpnn["data"]["params"]["fix_interface_residues"] is False


@pytest.mark.parametrize("case",["cycle","wrong_type","duplicate_id","two_inputs","missing_input","bad_number"])
def test_graph_rejects_invalid_dependencies_before_execution(case):
    g=graph.new_graph()
    if case=="cycle":
        # Sequence nodes can feed each other, but the cycle must still fail.
        g["nodes"]=[node("ligandmpnn","a"),node("ligandmpnn","b"),node("results","r")]
        g["edges"]=[edge("a","b"),edge("b","a"),edge("b","r")]
    elif case=="wrong_type":
        g["edges"][0]["target"]="rf3-3"
    elif case=="duplicate_id":
        g["nodes"][1]["id"]=g["nodes"][0]["id"]
    elif case=="two_inputs":
        g["edges"].append(edge("rfd3-1","rf3-3"))
    elif case=="missing_input":
        g["edges"].pop(1)
    else:
        g["nodes"][1]["data"]["params"]["num_designs"]=-1
    with pytest.raises(ValueError):graph.validate_graph(g)


def test_filter_keeps_upstream_type_and_format_constraints():
    g=graph.new_graph(); g["nodes"].insert(2,node("filter","f"))
    g["edges"][1:2]=[edge("rfd3-1","f"),edge("f","foundry_mpnn-2")]
    graph.validate_graph(g)
    g["nodes"][1]=node("rfdiffusion","rfd3-1")
    # The common MPNN adapter can use a shared PDB backend through a filter too.
    graph.validate_graph(g)
    g["nodes"][1]=node("rfd3","rfd3-1")
    g["edges"]=[edge("target-0","rfd3-1"),edge("rfd3-1","f"),edge("f","rf3-3"),edge("rf3-3","results-4")]
    g["nodes"]=[n for n in g["nodes"] if n["id"]!="foundry_mpnn-2"]
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


def test_adapters_disable_hidden_downstream_pipelines_and_bind_existing_signatures():
    import inspect
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS,_get_function
    from mn_molecular_design_studio.rfdiffusion3_input import ASSETS
    source={"run_dir":"/tmp/source","candidates_jsonl":"/tmp/source/candidates.jsonl","target":{"target_pdb":str(ASSETS/"5o45_cropped.pdb"),"target_chains":["A"],"hotspots":"A54"}}
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
    moved=copy.deepcopy(g);moved["nodes"][1]["position"]["y"]+=400
    p2=runner.run_graph_pipeline(moved,"moved","owner","project")
    assert len(calls)==3
    assert json.loads((p2/"result.json").read_text())["metrics"]["cached_steps"]==4
    changed=copy.deepcopy(g);changed["nodes"][2]["data"]["params"]["sequences"]+=1
    runner.run_graph_pipeline(changed,"changed","owner","project")
    assert calls==["rfd3","foundry_mpnn","rf3","foundry_mpnn","rf3"]
    assert (p1/"result.json").read_bytes()==original and (p1/"graph_snapshot.json").read_bytes()==snapshot


def test_parallel_branches_collect_distinct_candidate_ids(execution_harness):
    g,calls,_=execution_harness
    g["nodes"].append(node("esmfold","esm"));g["edges"] += [edge("foundry_mpnn-2","esm"),edge("esm","results-4")]
    p=runner.run_graph_pipeline(g,"branched","owner","project")
    assert set(calls)=={"rfd3","foundry_mpnn","rf3","esmfold"}
    candidates=json.loads((p/"result.json").read_text())["outputs"]["candidates"]
    assert len(candidates)==2 and len({c["candidate_id"] for c in candidates})==2


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
    assert result["outputs"]["pipeline_nodes"]["foundry_mpnn-2"]["status"]=="blocked"


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
    assert any("Choose a project" in c.value for c in app.caption)


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
