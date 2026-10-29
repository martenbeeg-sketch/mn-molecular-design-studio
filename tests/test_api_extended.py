"""Extended API/MCP journeys with private fixtures and no scientific processes."""
from __future__ import annotations

import asyncio
import csv
import hashlib
import io
import json
from pathlib import Path

import httpx
import pytest

from test_api_integrations import api_data, headers, import_compounds, screening_run


def write_structure(path, *, chains=(("A", 4), ("B", 4))):
    path.parent.mkdir(parents=True, exist_ok=True)
    records, serial = [], 0
    for chain, length in chains:
        for residue in range(1, length + 1):
            for atom, element in (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")):
                serial += 1
                records.append(f"ATOM  {serial:5d} {atom:^4s} GLY {chain}{residue:4d}    {float(serial):8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{20.:6.2f}          {element:>2s}")
        records.append("TER")
    path.write_text("\n".join(records) + "\nEND\n")
    return path


@pytest.fixture
def design_assets(api_data):
    from mn_molecular_design_studio import kejia_idr_input as kejia, peptide_trap_input as trap
    source = kejia.source_root()
    (source / ".git").mkdir(parents=True)
    (source / ".git/HEAD").write_text(kejia.SOURCE_REVISION)
    for name in ("conf_2.pdb", "conf_3.pdb"):
        write_structure(source / "templates/walle" / name)
    source = trap.library_root() / "source"
    for scaffold in trap.SCAFFOLDS.values():
        write_structure(source / "pdb_models" / scaffold["file"], chains=(("A", len(scaffold["pattern"])), ("B", 4)))
    (source / "design_after_seqthread.xml").write_text("<ROSETTASCRIPTS/>")
    (trap.library_root() / "source_manifest.json").write_text(json.dumps({"archive_sha256": trap.ARCHIVE_SHA256,
        "files": {p.relative_to(source).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in source.rglob("*") if p.is_file()}}))
    return api_data


def candidate_run(data, *, tool="kejia_idr", status="completed", success=True):
    from mn_molecular_design_studio.runtime import runs_root
    from mn_molecular_design_studio.workspace import STUDIO_SOURCE_ID
    root = Path(runs_root()) / "design" / ("api-" + tool)
    root.mkdir(parents=True)
    write_structure(root / "complex.pdb")
    for name, payload in {
        "metadata.json": {"run_id": root.name, "task_group": "design", "tool": tool, "status": status,
                          "job_type": tool, "created_at": "2026-10-26T00:00:00Z"},
        "input.json": {"params": {"generator_only": True}},
        "result.json": {"success": success},
        "molecular_design_studio.json": {"application": STUDIO_SOURCE_ID, "project_id": data["project"]["id"],
                                         "studio_owner_id": data["admin"]["id"]},
    }.items():
        (root / name).write_text(json.dumps(payload))
    candidate = {"candidate_id": "generated-1", "source_tool": tool, "complex_pdb": "complex.pdb",
        "stage": "generation.backbone" if tool == "kejia_idr" else "generation.backbone_sequence",
        "binder_chains": ["A"], "target_chains": ["B"], "binder_sequence": "GGGG",
        "metrics": {"result_kind": "generated", "rosetta_total_score": -9.5}, "raw_metadata": {}}
    folder = root / "artifacts/normalized_candidates"
    folder.mkdir(parents=True)
    (folder / "candidates.jsonl").write_text(json.dumps(candidate) + "\n")
    from mn_molecular_design_studio.job_catalog import reconcile
    reconcile()
    return root


@pytest.mark.parametrize("workflow", ["design.kejia_idr", "design.peptide_trap"])
def test_design_setup_schema_and_cpu_requirements_are_public(design_assets, workflow):
    from jsonschema import Draft202012Validator
    data, client = design_assets, design_assets["client"]
    response = client.get(f"/api/v1/workflows/{workflow}/setup", headers=headers(data, "student"))
    assert response.status_code == 200, response.text
    payload = response.json()
    assert payload["default_resources"] == {"gpu_device": "none", "cpu_cores": 4}
    assert payload["pipeline_nodes"][0]["workflow"] == workflow
    assert {"preview", "example", "stage-template"} <= set(payload["helpers"])
    schema = payload["schema"]["parameters_json_schema"]
    Draft202012Validator.check_schema(schema)
    bad = {"target_sequence": "GGGGGGGG", "cpu_cores": 0}
    assert list(Draft202012Validator(schema).iter_errors(bad))
    if workflow.endswith("peptide_trap"):
        assert list(Draft202012Validator(schema).iter_errors({"target_sequence": "GGGGGGGG", "scaffold_ids": ["unknown"]}))


def test_template_listing_download_and_staging_preserve_private_ownership(design_assets):
    data, client = design_assets, design_assets["client"]
    base = "/api/v1/workflows/design.kejia_idr/templates"
    listing = client.get(base, params={"library": "walle", "offset": 1, "limit": 1})
    assert listing.status_code == 200, listing.text
    assert listing.json()["total"] == 2 and listing.json()["items"][0]["id"] == "conf_3.pdb"
    response = client.get(base + "/conf_2.pdb/file", params={"library": "walle"})
    assert response.status_code == 200 and b"ATOM" in response.content
    staged = client.post(base + "/conf_2.pdb/stage", params={"library": "walle"}, headers=headers(data, "student"))
    assert staged.status_code == 201, staged.text
    payload = staged.json()
    assert payload["parameters"]["template_library"] == "upload"
    assert client.get(payload["file"]["download_url"], headers=headers(data, "student")).content == response.content
    assert client.get(payload["file"]["download_url"]).status_code == 404
    assert client.get(base, params={"library": "../../private"}).status_code == 422
    assert client.get(base + "/conf_2.pdb/file", params={"library": "mini"}).status_code == 404


def test_template_symlinks_cannot_escape_reference_library(design_assets):
    from mn_molecular_design_studio import kejia_idr_input
    path = design_assets["root"] / "private.pdb"
    path.write_text("private")
    (kejia_idr_input.source_root() / "templates/walle/escape.pdb").symlink_to(path)
    response = design_assets["client"].get("/api/v1/workflows/design.kejia_idr/templates/escape.pdb/file", params={"library": "walle"})
    assert response.status_code == 422 and "private" not in response.text


@pytest.mark.parametrize("workflow,example", [("design.kejia_idr", "tamarind"), ("design.peptide_trap", "c37")])
def test_example_preview_and_pipeline_draft_need_no_scientific_job(design_assets, workflow, example):
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    data, client = design_assets, design_assets["client"]
    before = client.get("/api/v1/jobs").json()["total"]
    example_response = client.post(f"/api/v1/workflows/{workflow}/examples/{example}", headers=headers(data, "student"))
    assert example_response.status_code == 201, example_response.text
    params = example_response.json()["parameters"]
    preview = client.post(f"/api/v1/workflows/{workflow}/preview", json={"parameters": params}, headers=headers(data, "student"))
    assert preview.status_code == 200, preview.text
    assert preview.json()["parameters"]["target_sequence"] == params["target_sequence"]
    if workflow.endswith("peptide_trap"):
        assert preview.json()["match_count"] == 1 and preview.json()["assets_ready"]
    else:
        assert preview.json()["template_count"] == 1
        assert client.post(f"/api/v1/workflows/{workflow}/preview", json={"parameters": params}).status_code == 404
    draft = client.post("/api/v1/pipelines/drafts", json={"project_id": data["project"]["id"],
        "template": workflow.split(".")[1], "parameters": params}, headers=headers(data, "student"))
    assert draft.status_code == 200, draft.text
    assert not draft.json()["saved"]
    assert draft.json()["graph"]["nodes"][0]["data"]["params"]["target_sequence"] == params["target_sequence"]
    assert client.get("/api/v1/jobs").json()["total"] == before
    assert MolecularDatabase().list_molecules(data["admin"])["total"] == 0


def test_preview_rejects_unknown_parameters_and_pages_matches(design_assets):
    client = design_assets["client"]
    path = "/api/v1/workflows/design.peptide_trap/preview"
    response = client.post(path, json={"parameters": {"target_sequence": "A" * 12}, "offset": 1, "limit": 2})
    assert response.status_code == 200, response.text
    assert response.json()["match_count"] == 8 and len(response.json()["matches"]) == 2
    assert client.post(path, json={"parameters": {"target_sequence": "A" * 12, "shell": "bad"}}).status_code == 422
    assert client.post(path, json={"parameters": {"target_sequence": "invalid"}}).status_code == 422
    assert client.post(path, json={"limit": 1001}).status_code == 422


def test_all_canvas_starters_are_available_as_unsaved_api_drafts(api_data):
    from mn_molecular_design_studio.pipeline_graph import STARTER_TEMPLATES, new_graph
    data, client = api_data, api_data["client"]
    items = client.get("/api/v1/pipelines/templates").json()["items"]
    assert {item["id"] for item in items} == set(STARTER_TEMPLATES.values())
    for item in items:
        response = client.post("/api/v1/pipelines/drafts", json={"project_id": data["project"]["id"], "template": item["id"]})
        assert response.status_code == 200, (item["id"], response.text)
        assert [n["data"]["kind"] for n in response.json()["graph"]["nodes"]] == [n["data"]["kind"] for n in new_graph(item["id"])["nodes"]]
    request = {"project_id": data["project"]["id"], "template": "kejia_idr"}
    assert client.post("/api/v1/pipelines/drafts", json=request, headers=headers(data, "other")).status_code == 404
    assert client.post("/api/v1/pipelines/drafts", json={**request, "template": "arbitrary"}).status_code == 422


@pytest.mark.parametrize("tool", ["kejia_idr", "peptide_trap"])
def test_pipeline_continuation_validates_source_candidate_and_project(api_data, tool):
    data, client = api_data, api_data["client"]
    run = candidate_run(data, tool=tool)
    request = {"template": tool, "project_id": data["project"]["id"], "source_job": "local:design:" + run.name,
               "selected_candidate_ids": ["generated-1"]}
    response = client.post("/api/v1/pipelines/drafts", json=request, headers=headers(data, "student"))
    assert response.status_code == 200, response.text
    source = response.json()["graph"]["nodes"][0]["data"]["params"]
    assert Path(source["source_run_dir"]) == run and source["selected_candidate_ids"] == ["generated-1"]
    assert client.post("/api/v1/pipelines/drafts", json={**request, "selected_candidate_ids": ["forged"]}).status_code == 422
    assert client.post("/api/v1/pipelines/drafts", json={**request, "parameters": {"source_run_dir": "/private"}}).status_code == 422
    assert client.post("/api/v1/pipelines/drafts", json=request, headers=headers(data, "other")).status_code == 404


def test_structure_preview_and_saved_camera_are_private(api_data):
    data, client = api_data, api_data["client"]
    pdb = write_structure(data["root"] / "input.pdb").read_bytes()
    uploaded = client.post("/api/v1/files?filename=template.pdb", content=pdb, headers=headers(data, "student")).json()
    response = client.post("/api/v1/structures/preview", json={"file_id": uploaded["id"]}, headers=headers(data, "student"))
    assert response.status_code == 200, response.text
    assert set(response.json()["inventory"]["chains"]) == {"A", "B"}
    assert response.json()["receptor_check"]["residues"] == 8
    assert client.post("/api/v1/structures/preview", json={"file_id": uploaded["id"]}).status_code == 404
    options = {"representation": "surface", "color": "chain-id", "background": "#112233", "camera_position": [1, 2, 3]}
    view = client.post("/api/v1/visualizations", json={"file_id": uploaded["id"], "options": options}, headers=headers(data, "student"))
    assert view.status_code == 201, view.text
    assert view.json()["view"]["options"]["camera_position"] == [1., 2., 3.]
    assert client.get(view.json()["file"]["download_url"]).status_code == 404
    assert client.post("/api/v1/visualizations", json={"file_id": uploaded["id"], "structure_id": "other"}).status_code == 422
    assert client.post("/api/v1/visualizations", json={"file_id": uploaded["id"], "options": {"camera_position": [1, 2]}}).status_code == 422


def test_local_execution_actions_and_results_preserve_identity(api_data, monkeypatch):
    from mn_molecular_design_studio.cluster import capabilities
    data, client = api_data, api_data["client"]
    monkeypatch.setattr(capabilities, "resources", lambda: (True, {"cpu_total": 12, "cpu_free": 10}))
    run = candidate_run(data)
    nodes = client.get("/api/v1/workstations").json()["nodes"]
    assert len(nodes) == 1 and nodes[0]["snapshot"]["resources"]["cpu_free"] == 10
    listed = client.get("/api/v1/executions").json()
    reference = next(r["reference"] for r in listed["items"] if r["run_id"] == run.name)
    assert reference == "local:design:" + run.name
    path = "/api/v1/executions/" + reference
    details = client.get(path)
    assert details.status_code == 200, details.text
    assert details.json()["job"]["run_id"] == run.name
    actions = {a["action"]: a for a in client.get(path + "/actions").json()["items"]}
    assert actions["hide"]["enabled"] and not actions["cancel"]["enabled"]
    assert client.post(path + "/actions", json={"action": "cancel"}).status_code == 409
    assert client.post(path + "/actions", json={"action": "hide", "parameters": {"shell": "bad"}}).status_code == 422
    assert client.post(path + "/actions", json={"action": "hide"}).status_code == 202
    assert client.get("/api/v1/executions").json()["total"] == 0
    assert client.get(path + "/results").json()["pending"] is False
    assert client.post(path + "/actions", json={"action": "restore"}).status_code == 202
    assert client.get("/api/v1/executions").json()["total"] == 1
    assert client.get(path, headers=headers(data, "other")).status_code == 404


def test_resumable_actions_reuse_worker_services_without_starting_a_worker(api_data, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.core import jobs
    data, client = api_data, api_data["client"]
    run = candidate_run(data, status="paused", success=None)
    (run / "worker_request.json").write_text("{}")
    calls = []
    monkeypatch.setattr(jobs, "resume_job", lambda folder, **kwargs: calls.append((folder, kwargs)))
    path = "/api/v1/executions/local:design:" + run.name + "/actions"
    assert next(a for a in client.get(path).json()["items"] if a["action"] == "resume")["enabled"]
    assert client.post(path, json={"action": "resume", "parameters": {"gpu_device": 2}}).status_code == 422
    response = client.post(path, json={"action": "resume", "parameters": {"gpu_device": "1"}})
    assert response.status_code == 202, response.text
    assert calls == [(run, {"gpu_device": "1"})]


def test_candidate_presentation_viewer_and_private_metric_exports(api_data):
    data, client = api_data, api_data["client"]
    run = candidate_run(data)
    base = "/api/v1/jobs/local:design:" + run.name
    response = client.get(base + "/presentation?limit=1")
    assert response.status_code == 200, response.text
    assert response.json()["items"][0]["candidate_id"] == "generated-1"
    view = client.get(base + "/candidates/generated-1/viewer")
    assert view.status_code == 200, view.text
    assert client.get(view.json()["layers"][0]["url"]).content == (run / "complex.pdb").read_bytes()
    for format, extension in (("CSV", "csv"), ("SVG", "svg"), ("HTML", "html")):
        report = client.post(base + "/reports", json={"format": format, "metric": "rosetta_total_score"})
        assert report.status_code == 201, report.text
        assert report.json()["row_count"] == 1 and report.json()["file"]["name"].endswith(extension)
        assert client.get(report.json()["file"]["download_url"]).status_code == 200
        assert client.get(report.json()["file"]["download_url"], headers=headers(data, "student")).status_code == 404
    assert client.get(base + "/presentation", headers=headers(data, "other")).status_code == 404


@pytest.mark.parametrize("format", ["CSV", "SVG", "HTML"])
def test_convergence_reports_use_all_three_saved_top_series_and_engine_units(api_data, format):
    data, client = api_data, api_data["client"]
    run = screening_run(data)
    history = {"score_units": "REU", "rounds": [{"round": number, "score_summaries": [
        {"top_n": count, "count": count, "mean_score": -10 - number, "boundary_score": -5 - number}
        for count in (100, 1000, 10000)]} for number in (1, 2)]}
    (run / "artifacts/active_learning/round_history.json").write_text(json.dumps(history))
    report = client.post("/api/v1/jobs/api-screening/reports", json={"format": format})
    assert report.status_code == 201, report.text
    assert report.json()["row_count"] == 6 and not report.json()["truncated"]
    content = client.get(report.json()["file"]["download_url"]).content.decode()
    if format == "CSV":
        assert {int(row["top_n"]) for row in csv.DictReader(io.StringIO(content))} == {100, 1000, 10000}
    else:
        assert "REU" in content and "kcal/mol" not in content
        assert all("Top " + str(n) in content for n in (100, 1000, 10000))


def test_workstation_settings_are_admin_only_and_never_return_credentials(api_data):
    from mn_molecular_design_studio.cluster import config
    data, client = api_data, api_data["client"]
    config.save({"node_token": "private-node-credential"})
    path = "/api/v1/admin/workstation"
    response = client.get(path)
    assert response.status_code == 200 and "private-node-credential" not in response.text
    assert response.json()["credential_configured"]
    assert client.get(path, headers=headers(data, "student")).status_code == 403
    updated = client.patch(path, json={"name": "Research PC", "accept_jobs": False})
    assert updated.status_code == 200, updated.text
    assert updated.json()["restart_required"] and config.load()["accept_jobs"] is False
    assert client.patch(path, json={"node_token": "replacement"}).status_code == 422
    assert client.patch(path, json={"heartbeat_seconds": 0}).status_code == 422
    assert client.patch(path, json={"name": None}).status_code == 422
    assert "private-node-credential" == config.load()["node_token"]
    operations = client.get("/api/v1/capabilities", headers=headers(data, "student")).json()["operations"]
    assert not any(o["path"].startswith("/api/v1/admin/") for o in operations)


def test_new_mcp_tools_reach_rest_with_student_identity(design_assets, monkeypatch):
    pytest.importorskip("mcp")
    from mcp import Client
    from mn_molecular_design_studio.api_client import StudioClient
    from mn_molecular_design_studio.mcp_server import create_server
    from mn_molecular_design_studio.cluster import capabilities
    data = design_assets
    monkeypatch.setattr(capabilities, "resources", lambda: (True, {"cpu_free": 8}))
    studio = StudioClient("http://studio.test", token=data["student_token"], transport=httpx.ASGITransport(app=data["app"]))

    async def journey():
        async with Client(create_server(studio)) as client:
            tools = {tool.name for tool in (await client.list_tools()).tools}
            assert {"workflow_setup", "preview_workflow", "create_pipeline_draft", "control_job", "export_result_report", "save_visualization"} <= tools
            for workflow in ("design.kejia_idr", "design.peptide_trap"):
                setup = await client.call_tool("workflow_setup", {"workflow_id": workflow})
                assert not setup.is_error and setup.structured_content["workflow"] == workflow
            example = await client.call_tool("stage_workflow_example", {"workflow_id": "design.kejia_idr", "example_id": "tamarind"})
            assert not example.is_error
            preview = await client.call_tool("preview_workflow", {"workflow_id": "design.kejia_idr", "parameters": example.structured_content["parameters"]})
            assert not preview.is_error and preview.structured_content["template_count"] == 1
            nodes = await client.call_tool("list_workstations", {})
            assert not nodes.is_error and len(nodes.structured_content["nodes"]) == 1
            rejected = await client.call_tool("studio_api_request", {"method": "PATCH", "path": "/api/v1/admin/workstation", "body": {"accept_jobs": False}})
            assert rejected.is_error
    asyncio.run(journey())


def test_distributed_results_get_is_read_only_and_post_requests_inspection(api_data):
    from mn_molecular_design_studio.cluster import config
    from mn_molecular_design_studio.cluster.store import Store
    data, client = api_data, api_data["client"]
    config.save({"mode": "coordinator", "coordinator_url": "http://studio.test"})
    store = Store()
    node = store.enroll("remote")
    store.heartbeat(node["node_id"], {}, [{"task_group": "design", "run_id": "remote-design",
        "owner_id": data["admin"]["id"], "project_id": data["project"]["id"], "status": "completed", "tool": "kejia_idr"}])
    listed = client.get("/api/v1/executions").json()
    assert listed["total"] == 1
    reference = listed["items"][0]["reference"]
    assert reference.startswith("cluster:")
    path = "/api/v1/executions/" + reference
    assert client.get(path).json()["job"]["run_id"] == "remote-design"
    with store.connect() as db:
        before = db.execute("SELECT COUNT(*) FROM commands").fetchone()[0]
    response = client.get(path + "/results")
    assert response.status_code == 200 and response.json()["inspection_required"]
    with store.connect() as db:
        assert db.execute("SELECT COUNT(*) FROM commands").fetchone()[0] == before
    response = client.post(path + "/results")
    assert response.status_code == 202, response.text
    command = response.json()["command_id"]
    status = client.get(path + "/results")
    assert status.status_code == 200 and status.json()["pending"]
    assert status.json()["command"]["command_id"] == command
    assert client.get(path, headers=headers(data, "other")).status_code in {403, 404}
    assert client.post(path + "/actions", json={"action": "retry"}).status_code == 409


def test_connected_worker_forwards_controls_with_callers_token(api_data, monkeypatch):
    from mn_molecular_design_studio import api
    from mn_molecular_design_studio.cluster import config, transport
    data, client = api_data, api_data["client"]
    config.save({"mode": "worker", "coordinator_url": "http://coordinator.test", "node_token": "node-secret"})
    # Authentication itself is covered separately; here inspect outbound control identity.
    data["app"].dependency_overrides[api._current_user] = lambda: data["student"]
    calls = []

    def request(path, payload=None, **kwargs):
        calls.append((path, payload, kwargs))
        return {"forwarded": True}

    monkeypatch.setattr(transport, "request", request)
    try:
        path = "/api/v1/executions/cluster:remote/actions"
        response = client.post(path, json={"action": "hide"}, headers=headers(data, "student"))
        assert response.status_code == 202, response.text
        assert calls == [(path, {"action": "hide", "parameters": {}}, {"token": data["student_token"], "method": "POST"})]
        assert "node-secret" not in str(calls)
        calls.clear()
        assert client.get("/api/v1/workstations", headers=headers(data, "student")).status_code == 200
        assert calls[0][0] == "/api/v1/cluster/nodes" and calls[0][2]["token"] == data["student_token"]
    finally:
        data["app"].dependency_overrides.clear()


@pytest.mark.parametrize("legacy", [False, True])
def test_bound_viewer_keeps_target_and_ligand_separate_without_creating_files(api_data, legacy):
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    from mn_molecular_design_studio import bound_docking
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    from mn_molecular_design_studio.auth import user_input_root
    data, client = api_data, api_data["client"]
    database = MolecularDatabase()
    receptor = write_structure(data["root"] / "receptor.pdb", chains=(("A", 4),))
    target = bound_docking.register_receptor(database, data["admin"], receptor, project_id=data["project"]["id"])
    imported = database.import_records(data["admin"], [{"name": "Docked compound", "components": [
        {"type": "small_molecule", "chain": "L", "smiles": "CCO"}]}], project_id=data["project"]["id"])
    molecule = Chem.MolFromSmiles("CCO")
    coordinates = Chem.Conformer(molecule.GetNumAtoms())
    for i in range(molecule.GetNumAtoms()):
        coordinates.SetAtomPosition(i, Point3D(5 + i, 2, 3))
    molecule.AddConformer(coordinates)
    pose = data["root"] / "pose.sdf"
    with Chem.SDWriter(str(pose)) as writer:
        writer.write(molecule)
    identifier = imported["molecule_ids"][0]
    context = {**target, "database_structure_policy": "best_docking_pose"}
    if legacy:
        context.pop("target_structure_id")
    structure = database.attach_structure(imported["group_id"], identifier, {"path": str(pose), "label": "Docked pose"}, context=context)
    files_root = user_input_root(data["admin"]["id"]) / "api_files"
    before = list(files_root.glob("*"))
    response = client.get(f"/api/v1/molecules/{identifier}/viewer", params={"structure_id": structure})
    assert response.status_code == 200, response.text
    view = response.json()
    assert view["coordinates_aligned"] and len(view["layers"]) == 2
    assert [layer["molecule_id"] for layer in view["layers"]] == [target["target_molecule_id"], identifier]
    for layer in view["layers"]:
        assert client.get(layer["url"]).status_code == 200
        assert client.get(layer["url"], headers=headers(data, "other")).status_code == 404
    assert list(files_root.glob("*")) == before
    assert database.list_molecules(data["admin"])["total"] == 2


def test_binary_download_descriptors_retain_encoded_query_parameters():
    from mn_molecular_design_studio.api_client import StudioClient

    def transport(request):
        if request.url.path == "/api/v1/auth/me":
            return httpx.Response(200, json={"id": "student", "role": "user"})
        return httpx.Response(200, content=b"coordinates", headers={"content-type": "chemical/x-pdb"})

    async def request():
        client = StudioClient("http://studio.test", token="student-token", transport=httpx.MockTransport(transport))
        response = await client.request("GET", "/api/v1/jobs/run/artifact", query={"path": "artifacts/model one.pdb"})
        assert httpx.URL(response["download_url"]).params["path"] == "artifacts/model one.pdb"
        assert "student-token" not in response["download_url"]
    asyncio.run(request())


def test_group_schema_api_retains_owner_write_permission(api_data):
    data, client = api_data, api_data["client"]
    imported = import_compounds(data)
    group = client.get("/api/v1/molecules/groups/" + imported["group_id"]).json()
    schema = group["schema"]
    path = "/api/v1/molecules/groups/" + imported["group_id"] + "/schema"
    response = client.put(path, json=schema)
    assert response.status_code == 200, response.text
    assert client.put(path, json=schema, headers=headers(data, "other")).status_code == 404


def test_pipeline_campaign_tile_submits_through_normal_api_without_dispatching_inference(api_data, tmp_path, monkeypatch):
    from mn_molecular_design_studio import api
    from mn_molecular_design_studio.auth import user_input_root
    from mn_molecular_design_studio.distributed_campaign_config_ui import _default_validation_stages
    from mn_molecular_design_studio.engine.workflow_engine.services import local_automation
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, SCHEMA

    data, client = api_data, api_data["client"]
    target_path = user_input_root(data["admin"]["id"]) / "campaign-target.pdb"
    write_structure(target_path, chains=(("A", 4),))

    def pipeline_node(kind, identifier):
        spec = NODE_CATALOG[kind]
        params = {field["name"]: field["default"] for field in spec["fields"]}
        return {"id": identifier, "type": "studio", "position": {"x": 0, "y": 0},
                "data": {"kind": kind, "label": spec["label"], "params": params}}

    target = pipeline_node("target", "target")
    target["data"]["params"].update(target_pdb=str(target_path), target_chains="A")
    rfd4 = pipeline_node("rfd4", "generator")
    rfd4["data"]["params"]["distributed_campaign"] = {
        "enabled": True,
        "accepted_target": 2,
        "max_attempts": 5,
        "acceptance_level": "validated",
        "batch_size": 1,
        "max_concurrent_batches": 2,
        "max_concurrent_validations": 3,
        "workstation_selection": "eligible",
        "gpu_selection": "auto",
        "validation_stages": _default_validation_stages(),
    }
    results = pipeline_node("results", "results")
    graph = {
        "schema_version": SCHEMA, "id": "b" * 32, "name": "API campaign fixture",
        "nodes": [target, rfd4, results],
        "edges": [
            {"id": "target-generator", "source": "target", "sourceHandle": "out", "target": "generator", "targetHandle": "in"},
            {"id": "generator-results", "source": "generator", "sourceHandle": "out", "target": "results", "targetHandle": "in"},
        ],
        "viewport": {"x": 0, "y": 0, "zoom": 1},
    }

    target_run = tmp_path / "enqueued-campaign"
    target_run.mkdir()
    captured = {}
    original_call = api._call

    def call(function, *args, **kwargs):
        if getattr(function, "__name__", "") == "enqueue_graph_pipeline":
            captured["graph"] = args[0]
            captured["gpu_device"] = kwargs.get("gpu_device")
            return target_run
        return original_call(function, *args, **kwargs)

    monkeypatch.setattr(api, "_call", call)
    monkeypatch.setattr(api, "mark_studio_job", lambda result, **kwargs: result)
    monkeypatch.setattr(local_automation, "job_reference", lambda run_dir: {"run_id": "fixture-campaign", "run_dir": str(run_dir)})
    response = client.post("/api/v1/pipelines/run", json={"project_id": data["project"]["id"], "graph": graph})

    assert response.status_code == 202, response.text
    submitted = next(node for node in captured["graph"]["nodes"] if node["data"]["kind"] == "rfd4")
    config = submitted["data"]["params"]["distributed_campaign"]
    assert config["accepted_target"] == 2
    assert config["max_attempts"] == 5
    assert config["workstation_selection"] == "eligible"
    assert config["max_concurrent_batches"] == 2
    assert config["max_concurrent_validations"] == 3
    assert captured["gpu_device"] == "all"


def test_shared_generator_objective_validates_and_runs_through_normal_pipeline_api(api_data, tmp_path, monkeypatch):
    from mn_molecular_design_studio import api
    from mn_molecular_design_studio.auth import user_input_root
    from mn_molecular_design_studio.distributed_campaign_config_ui import _default_validation_stages
    from mn_molecular_design_studio.engine.workflow_engine.services import local_automation
    from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG, SCHEMA

    data, client = api_data, api_data["client"]
    target_path = user_input_root(data["admin"]["id"]) / "shared-campaign-target.pdb"
    write_structure(target_path, chains=(("A", 4),))

    def pipeline_node(kind, identifier):
        spec = NODE_CATALOG[kind]
        params = {field["name"]: field["default"] for field in spec["fields"]}
        return {"id": identifier, "type": "studio", "position": {"x": 0, "y": 0},
                "data": {"kind": kind, "label": spec["label"], "params": params}}

    target = pipeline_node("target", "target")
    target["data"]["params"].update(target_pdb=str(target_path), target_chains="A")
    bindcraft, rfd4, results = (pipeline_node(kind, identifier) for kind, identifier in (
        ("bindcraft2", "bindcraft"), ("rfd4", "rfd4"), ("results", "results"),
    ))
    graph = {
        "schema_version": SCHEMA, "id": "f" * 32, "name": "Shared PD-L1 objective",
        "workflow_mode": "native_comparative",
        "nodes": [target, bindcraft, rfd4, results],
        "edges": [
            {"id": "target-bindcraft", "source": "target", "sourceHandle": "out", "target": "bindcraft", "targetHandle": "in"},
            {"id": "target-rfd4", "source": "target", "sourceHandle": "out", "target": "rfd4", "targetHandle": "in"},
            {"id": "bindcraft-results", "source": "bindcraft", "sourceHandle": "out", "target": "results", "targetHandle": "in"},
            {"id": "rfd4-results", "source": "rfd4", "sourceHandle": "out", "target": "results", "targetHandle": "in"},
        ],
        "shared_campaign": {
            "enabled": True, "generator_nodes": ["bindcraft", "rfd4"],
            "accepted_target": 1, "max_attempts": 4,
            "max_concurrent_jobs": 2, "max_concurrent_generation_jobs": 2,
            "max_concurrent_validations": 1,
            "validation_stages": _default_validation_stages(),
        },
        "viewport": {"x": 0, "y": 0, "zoom": 1},
    }
    validate_response = client.post("/api/v1/pipelines/validate", json={
        "project_id": data["project"]["id"], "graph": graph,
    })
    assert validate_response.status_code == 200, validate_response.text
    assert validate_response.json()["valid"] is True

    target_run = tmp_path / "enqueued-shared-pipeline"
    target_run.mkdir()
    captured = {}
    original_call = api._call

    def call(function, *args, **kwargs):
        if getattr(function, "__name__", "") == "enqueue_graph_pipeline":
            captured["graph"] = args[0]
            return target_run
        return original_call(function, *args, **kwargs)

    monkeypatch.setattr(api, "_call", call)
    monkeypatch.setattr(api, "mark_studio_job", lambda result, **kwargs: result)
    monkeypatch.setattr(local_automation, "job_reference", lambda run_dir: {"run_id": "shared-run", "run_dir": str(run_dir)})
    run_response = client.post("/api/v1/pipelines/run", json={
        "project_id": data["project"]["id"], "graph": graph,
    })
    assert run_response.status_code == 202, run_response.text
    assert captured["graph"]["shared_campaign"]["accepted_target"] == 1
    assert captured["graph"]["shared_campaign"]["generator_nodes"] == ["bindcraft", "rfd4"]
    assert captured["graph"]["shared_campaign"]["workstation_selection"] == "auto"
