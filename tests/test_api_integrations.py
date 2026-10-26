"""API-only journeys use isolated accounts, files and tiny scientific fixtures."""
from __future__ import annotations

import asyncio
import base64
import json
from pathlib import Path

import httpx
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from fastapi.testclient import TestClient


@pytest.fixture
def api_data(tmp_path, monkeypatch):
    for name, value in {
        "MN_MOLECULAR_DESIGN_STUDIO_HOME": tmp_path / "studio",
        "MN_MOLECULAR_DESIGN_STUDIO_APP_HOME": tmp_path / "studio",
        "MN_MOLECULAR_DESIGN_STUDIO_RESULTS_PARENT": tmp_path / "results",
        "MN_MOLECULAR_DESIGN_STUDIO_REFERENCE_DIR": tmp_path / "references",
        "MDS_CLUSTER_HOME": tmp_path / "cluster",
        "MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE": "production",
    }.items():
        monkeypatch.setenv(name, str(value))
    from mn_molecular_design_studio import auth, api, workspace, molecular_database_api
    from mn_molecular_design_studio.cluster import config
    workspace.ensure_studio_workspace()
    config.save({"mode": "disabled"})
    monkeypatch.setattr(molecular_database_api, "sync_workspace", lambda *a, **k: [])
    admin = auth.create_initial_admin("instructor", "Instructor", "instructor-password-42")
    project = workspace.create_project("API course", owner_id=admin["id"])
    student = auth.create_user("student", "Student", "temporary-password-42", actor=admin)
    auth.change_password(student["id"], "student-password-42", forced=True)
    auth.grant_project_access(project["id"], student["id"], actor=admin)
    other = auth.create_user("otherstudent", "Other", "temporary-password-42", actor=admin)
    auth.change_password(other["id"], "otherstudent-password-42", forced=True)
    admin_token = auth.authenticate("instructor", "instructor-password-42")[1]
    student_token = auth.authenticate("student", "student-password-42")[1]
    other_token = auth.authenticate("otherstudent", "otherstudent-password-42")[1]
    with TestClient(api.app) as client:
        client.headers["Authorization"] = "Bearer " + admin_token
        yield {"client": client, "app": api.app, "admin": admin, "student": student,
               "other": other, "project": project, "root": tmp_path,
               "admin_token": admin_token, "student_token": student_token, "other_token": other_token}


def headers(data, role):
    return {"Authorization": "Bearer " + data[role + "_token"]}


def import_compounds(data):
    response = data["client"].post("/api/v1/molecules/imports", json={
        "method": "smiles", "text": "CCO ethanol\nCCN ethylamine", "project_id": data["project"]["id"]})
    assert response.status_code == 201, response.text
    return response.json()


def screening_run(data, count=11):
    from mn_molecular_design_studio.runtime import runs_root
    from mn_molecular_design_studio.workspace import STUDIO_SOURCE_ID
    run = Path(runs_root()) / "campaign" / "api-screening"
    folder = run / "artifacts/active_learning"
    folder.mkdir(parents=True)
    (run / "molecular_design_studio.json").write_text(json.dumps({"application": STUDIO_SOURCE_ID,
        "studio_owner_id": data["admin"]["id"], "project_id": data["project"]["id"]}))
    (run / "metadata.json").write_text(json.dumps({"status": "completed", "database_collection_policy": "manual"}))
    (run / "result.json").write_text(json.dumps({"success": True}))
    pq.write_table(pa.Table.from_pylist([{"compound_id": f"c-{i}", "compound_name": f"Compound {i}",
        "compound_key": i, "library_id": "lib", "smiles": "C" * (i % 8 + 1), "score": -20. + i / 100}
        for i in range(count)]), folder / "docking_scores.parquet", row_group_size=3)
    (folder / "round_history.json").write_text(json.dumps({"rounds": [{"round": 1}], "score_units": "REU"}))
    return run


def test_private_upload_retains_extension_and_never_imports_compounds(api_data):
    from mn_molecular_design_studio import api_files
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    data, client = api_data, api_data["client"]
    response = client.post("/api/v1/files?filename=library.csv", content=b"id,smiles\na,CCO\n")
    assert response.status_code == 201, response.text
    record = response.json()
    path = api_files.resolve(data["admin"], record["id"])
    assert path.suffix == ".csv" and path.read_bytes() == b"id,smiles\na,CCO\n"
    assert client.get(record["download_url"]).content == path.read_bytes()
    assert client.get(record["download_url"], headers=headers(data, "student")).status_code == 404
    assert client.get(record["download_url"] + "/metadata", headers=headers(data, "student")).status_code == 404
    assert MolecularDatabase().list_molecules(data["admin"])["total"] == 0
    forged = client.post("/api/v1/workflows/validate", headers=headers(data, "student"), json={
        "project_id": data["project"]["id"], "schema_version": 1, "workflow": "library.prepare",
        "parameters": {"source_files": [record["reference"]]}})
    assert forged.status_code == 404


def test_library_registration_sharing_and_lineage_are_metadata_only(api_data):
    from mn_molecular_design_studio.external_libraries import LibraryCatalog, FEATURE_SCHEMA
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    data, client = api_data, api_data["client"]
    source = data["root"] / "large-library.csv"
    source.write_text("ID,SMILES\na,CCO\nb,CCN\n")
    request = {"paths": [str(source)], "name": "External library", "id_column": "ID", "smiles_column": "SMILES"}
    assert client.post("/api/v1/libraries", json=request, headers=headers(data, "student")).status_code == 403
    response = client.post("/api/v1/libraries", json=request)
    assert response.status_code == 201, response.text
    library = response.json()
    assert client.get("/api/v1/libraries", headers=headers(data, "student")).json()["items"] == []
    assert client.get("/api/v1/libraries/" + library["id"], headers=headers(data, "student")).status_code == 404
    shared = client.patch("/api/v1/libraries/" + library["id"] + "/access", json={"project_ids": [data["project"]["id"]]})
    assert shared.status_code == 200, shared.text
    manifest = data["root"] / "fingerprints.json"
    manifest.write_text(json.dumps({"schema_version": FEATURE_SCHEMA, "library_id": library["id"],
        "status": "completed", "row_count": 2, "created_at": "2026-10-26T00:00:00Z"}))
    LibraryCatalog().attach_preparation(library["id"], manifest)
    path = "/api/v1/libraries/" + library["id"] + "/preparations?kind=fingerprints"
    assert client.get(path, headers=headers(data, "student")).json()["total"] == 1
    assert client.get(path, headers=headers(data, "other")).status_code == 404
    assert client.post("/api/v1/libraries", json=request).json()["project_ids"] == [data["project"]["id"]]
    assert MolecularDatabase().list_molecules(data["admin"])["total"] == 0


def test_screening_pagination_filters_original_rank_and_streamed_download(api_data, monkeypatch):
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    data, client = api_data, api_data["client"]
    screening_run(data, count=4101)
    monkeypatch.setattr(pq, "read_table", lambda *a, **k: pytest.fail("Must use bounded Parquet batches"))
    path = "/api/v1/screening/campaign/api-screening"
    response = client.post(path + "/search", json={"offset": 2047, "limit": 5})
    assert response.status_code == 200, response.text
    page = response.json()
    assert [r["rank"] for r in page["items"]] == list(range(2048, 2053))
    assert page["total"] == 4101 and page["has_more"] and page["next_offset"] == 2052
    filtered = client.post(path + "/search", json={"compound_ids": ["c-2", "c-2000", "c-4100"], "offset": 1, "limit": 2}).json()
    assert [r["rank"] for r in filtered["items"]] == [2001, 4101]
    assert filtered["total"] is None and filtered["has_more"] is False
    assert client.post(path + "/search", json={"score_min": 10, "score_max": -10}).status_code == 422
    assert client.post(path + "/search", json={"limit": 1001}).status_code == 422
    assert client.post(path + "/search", json={}, headers=headers(data, "other")).status_code == 404
    assert client.get(path + "/files/scores").content.startswith(b"PAR1")
    assert client.get(path + "/files/secrets").status_code == 404
    assert client.get(path).json()["score_units"] == "REU"
    assert MolecularDatabase().list_molecules(data["admin"])["total"] == 0


def test_final_screening_import_is_explicit_bounded_and_owner_only(api_data):
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    data, client = api_data, api_data["client"]
    run = screening_run(data)
    path = "/api/v1/screening/campaign/api-screening/final-selection"
    assert client.post(path, json={"limit": 1001}).status_code == 422
    assert client.post(path, json={"limit": 2}, headers=headers(data, "student")).status_code == 404
    (run / "metadata.json").write_text(json.dumps({"status": "running"}))
    assert client.post(path, json={"limit": 2}).status_code == 422
    (run / "metadata.json").write_text(json.dumps({"status": "completed"}))
    result = client.post(path, json={"limit": 2})
    assert result.status_code == 200, result.text
    assert result.json()["imported_count"] == 2
    assert client.post(path, json={"limit": 2}).json()["imported_count"] == 2
    assert MolecularDatabase().list_molecules(data["admin"])["total"] == 2
    assert pq.ParquetFile(run / "artifacts/active_learning/docking_scores.parquet").metadata.num_rows == 11


@pytest.mark.parametrize("format,extension", [("CSV", "csv"), ("SMILES", "smi"), ("SDF", "sdf")])
def test_exports_are_private_downloadable_files(api_data, format, extension):
    data, client = api_data, api_data["client"]
    imported = import_compounds(data)
    response = client.post("/api/v1/molecules/exports", json={"molecule_ids": imported["molecule_ids"], "format": format})
    assert response.status_code == 201, response.text
    exported = response.json()
    assert exported["name"].endswith("." + extension)
    assert client.get(exported["download_url"]).status_code == 200
    assert client.get(exported["download_url"], headers=headers(data, "other")).status_code == 404
    assert client.post("/api/v1/molecules/exports", json={"molecule_ids": imported["molecule_ids"], "format": format}, headers=headers(data, "other")).status_code == 404


def test_reviewed_table_import_persists_report_and_original_source(api_data):
    data, client = api_data, api_data["client"]
    content = b"ID,SMILES,Name\na,CCO,Ethanol\nb,OCC,Duplicate\nc,invalid,Bad\n"
    upload = client.post("/api/v1/files?filename=molecules.csv", content=content).json()
    analyzed = client.post("/api/v1/molecules/imports/analyze", json={"file_id": upload["id"], "smiles_column": "SMILES", "id_column": "ID"})
    assert analyzed.status_code == 200, analyzed.text
    assert analyzed.json()["summary"]["duplicate_structure_count"] == 1
    imported = client.post("/api/v1/molecules/imports", json={"method": "table", "file_id": upload["id"],
        "review_small_molecules": True, "chain_columns": {"SMILES": {"type": "small_molecule", "chain": "L"}},
        "id_column": "ID", "name_column": "Name"})
    assert imported.status_code == 201, imported.text
    assert len(imported.json()["molecule_ids"]) == 1
    batches = client.get("/api/v1/molecules/groups/" + imported.json()["group_id"] + "/imports").json()["items"]
    assert len(batches) == 1
    report = "/api/v1/molecules/imports/" + batches[0]["id"]
    assert client.get(report + "/rows?limit=2&offset=1").json()["total"] == 3
    assert client.get(report + "/source").content == content
    assert client.get(report + "/source", headers=headers(data, "other")).status_code == 404


def test_group_role_crop_and_membership_operations_use_existing_identities(api_data):
    data, client = api_data, api_data["client"]
    imported = client.post("/api/v1/molecules/imports", json={"method": "records", "records": [{
        "name": "Two chains", "components": [{"type": "protein", "chain": "A", "sequence": "ACDEFG"},
                                              {"type": "protein", "chain": "B", "sequence": "KLMNPQ"}]}]}).json()
    molecule, group = imported["molecule_ids"][0], imported["group_id"]
    assert client.patch(f"/api/v1/molecules/groups/{group}/molecules/{molecule}/chains/A/role", json={"role": "target"}).status_code == 200
    result = client.post(f"/api/v1/molecules/{molecule}/crop", json={"chains": ["A"]})
    assert result.status_code == 201, result.text
    cropped = client.get("/api/v1/molecules/" + result.json()["molecule_ids"][0]).json()
    assert [c["chain"] for c in cropped["components"]] == ["A"]
    assert client.patch(f"/api/v1/molecules/groups/{group}/project", json={"project_id": data["project"]["id"]}).status_code == 200
    assert client.post(f"/api/v1/molecules/groups/{group}/members/remove", json={"molecule_ids": [molecule]}, headers=headers(data, "student")).status_code == 404
    assert client.post(f"/api/v1/molecules/groups/{group}/members/remove", json={"molecule_ids": [molecule]}).status_code == 200
    # Removing one membership never removes another saved collection.
    assert client.get("/api/v1/molecules/" + cropped["id"]).status_code == 200


def test_account_administration_grants_disable_and_password_reset(api_data):
    data, client = api_data, api_data["client"]
    path = "/api/v1/admin/users"
    assert client.get(path, headers=headers(data, "student")).status_code == 403
    created = client.post(path, json={"username": "newstudent", "display_name": "New Student",
        "temporary_password": "temporary-password-42", "project_ids": [data["project"]["id"]]})
    assert created.status_code == 201, created.text
    account = created.json()
    assert "temporary-password-42" not in created.text and "password_hash" not in account
    assert account["project_ids"] == [data["project"]["id"]] and account["must_change_password"]
    member = path + "/" + account["id"]
    assert client.put(member + "/projects", json={"project_ids": ["not-a-project"]}).status_code == 422
    assert client.get(member).json()["project_ids"] == account["project_ids"]
    assert client.put(member + "/projects", json={"project_ids": []}).json()["project_ids"] == []
    assert client.post(member + "/password", json={"temporary_password": "changed-password-42"}).json()["sessions_revoked"]
    assert client.patch(member, json={"disabled": True}).json()["disabled"]
    assert client.patch(path + "/" + data["admin"]["id"], json={"disabled": True}).status_code == 422
    assert client.patch(member, json={"disabled": "yes"}).status_code == 422


def test_capabilities_exclude_node_routes_and_nonadmin_management(api_data):
    client = api_data["client"]
    response = client.get("/api/v1/capabilities", headers=headers(api_data, "student"))
    assert response.status_code == 200
    paths = [r["path"] for r in response.json()["operations"]]
    assert "/api/v1/libraries" in paths
    assert not any("/cluster/node/" in p or "/admin/" in p for p in paths)
    assert client.get("/api/v1/auth/me").json()["development_mode"] is False
    operation = client.get("/api/v1/operations/schema", params={"path": "/api/v1/molecules/exports", "method": "POST"}, headers=headers(api_data, "student"))
    assert operation.status_code == 200 and "ExportRequest" in operation.json()["components"]["schemas"]
    assert client.get("/api/v1/operations/schema", params={"path": "/api/v1/admin/users", "method": "POST"}, headers=headers(api_data, "student")).status_code == 404
    assert client.get("/api/v1/operations/schema", params={"path": "/api/v1/cluster/node/claim", "method": "POST"}).status_code == 404


def test_local_request_retries_never_submit_twice(api_data, monkeypatch):
    from mn_molecular_design_studio import api
    client = api_data["client"]
    calls = []
    monkeypatch.setattr(api, "_cluster_submission", lambda *a, **k: None)
    monkeypatch.setattr(api, "_guard_source_job", lambda *a, **k: None)
    monkeypatch.setattr(api, "submit_workflow_request", lambda path: calls.append(path) or {"job": {"run_id": "once", "status": "queued"}})
    payload = {"project_id": api_data["project"]["id"], "schema_version": 1, "workflow": "analysis.ipsae", "parameters": {}}
    key = {"Idempotency-Key": "logical-submission"}
    first = client.post("/api/v1/workflows", json=payload, headers=key)
    second = client.post("/api/v1/workflows", json=payload, headers=key)
    assert first.status_code == second.status_code == 202, first.text
    assert first.json() == second.json() and len(calls) == 1
    assert client.post("/api/v1/workflows", json={**payload, "parameters": {"changed": True}}, headers=key).status_code == 409
    status = client.get("/api/v1/submissions/status", params={"kind": "workflow", "key": "logical-submission"})
    assert status.json()["status"] == "completed" and status.json()["response"] == first.json()
    assert client.get("/api/v1/submissions/status", params={"kind": "workflow", "key": "logical-submission"}, headers=headers(api_data, "student")).status_code == 404


def test_artifact_streaming_stays_inside_the_authorized_run(api_data, monkeypatch):
    from mn_molecular_design_studio import api
    root = api_data["root"] / "artifact-job"
    root.mkdir()
    (root / "scores.txt").write_text("0123456789")
    outside = api_data["root"] / "private.txt"
    outside.write_text("private")
    (root / "escape").symlink_to(outside)
    monkeypatch.setattr(api, "_resolve_api_job", lambda reference, user: {"run_dir": str(root)})
    client = api_data["client"]
    path = "/api/v1/jobs/run/artifact"
    response = client.get(path, params={"path": "scores.txt"}, headers={"Range": "bytes=2-4"})
    assert response.status_code == 206 and response.content == b"234"
    for name in ("../private.txt", str(outside), "escape"):
        assert client.get(path, params={"path": name}).status_code == 404


def test_file_cleanup_cannot_remove_another_users_file(api_data):
    data, client = api_data, api_data["client"]
    record = client.post("/api/v1/files?filename=input.smi", content=b"CCO ethanol").json()
    assert client.delete(record["download_url"], headers=headers(data, "other")).status_code == 404
    assert client.delete(record["download_url"]).status_code == 200
    assert client.get(record["download_url"]).status_code == 404


def test_uncertain_local_submission_is_not_automatically_repeated(api_data, monkeypatch):
    from mn_molecular_design_studio import api
    client, calls = api_data["client"], []
    monkeypatch.setattr(api, "_cluster_submission", lambda *a, **k: None)
    monkeypatch.setattr(api, "_guard_source_job", lambda *a, **k: None)

    def failed(path):
        calls.append(path)
        raise OSError("Lost the scheduler reply")

    monkeypatch.setattr(api, "submit_workflow_request", failed)
    payload = {"project_id": api_data["project"]["id"], "schema_version": 1, "workflow": "analysis.ipsae", "parameters": {}}
    key = {"Idempotency-Key": "uncertain-request"}
    assert client.post("/api/v1/workflows", json=payload, headers=key).status_code == 422
    assert client.post("/api/v1/workflows", json=payload, headers=key).status_code == 409
    assert len(calls) == 1


def test_saved_pipeline_revisions_are_accessible_only_in_the_original_scope(api_data):
    from mn_molecular_design_studio.pipeline_graph import new_graph
    data, client = api_data, api_data["client"]
    graph = new_graph()
    saved = client.post("/api/v1/pipelines/save", json={"project_id": data["project"]["id"], "graph": graph, "version": True})
    assert saved.status_code == 200, saved.text
    path = f"/api/v1/pipelines/{graph['id']}/versions"
    query = {"project_id": data["project"]["id"]}
    versions = client.get(path, params=query).json()["items"]
    assert len(versions) == 1
    version_path = path + "/" + versions[0]["revision"]
    assert client.get(version_path, params=query).json()["graph"]["id"] == graph["id"]
    assert client.get(version_path, params=query, headers=headers(data, "other")).status_code == 404
    assert client.get(version_path, params=query, headers=headers(data, "student")).status_code == 404


def test_bound_pose_api_prepares_mapped_md_files_without_merging_molecules(api_data):
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    from mn_molecular_design_studio import bound_docking
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    data, client = api_data, api_data["client"]
    receptor = data["root"] / "receptor.pdb"
    receptor.write_text("ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C\nEND\n")
    database = MolecularDatabase()
    target = bound_docking.register_receptor(database, data["admin"], receptor, project_id=data["project"]["id"])
    imported = database.import_records(data["admin"], [{"name": "Pose", "components": [{"type": "small_molecule", "chain": "L", "smiles": "CCO"}]}], project_id=data["project"]["id"])
    molecule = Chem.MolFromSmiles("CCO")
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    for i in range(molecule.GetNumAtoms()):
        conformer.SetAtomPosition(i, Point3D(5 + i, 2, 3))
    molecule.AddConformer(conformer)
    pose = data["root"] / "pose.sdf"
    with Chem.SDWriter(str(pose)) as writer:
        writer.write(molecule)
    identifier = imported["molecule_ids"][0]
    structure = database.attach_structure(imported["group_id"], identifier, {"path": str(pose), "label": "Pose"},
        context={**target, "database_structure_policy": "best_docking_pose"})
    before = database.list_molecules(data["admin"])["total"]
    prepared = client.post(f"/api/v1/molecules/{identifier}/bound-inputs", json={"structure_id": structure})
    assert prepared.status_code == 201, prepared.text
    result = prepared.json()
    assert result["workflow_inputs"]["complex_pdb"].startswith("studio-file:///")
    assert client.get(result["files"]["complex_pdb"]["download_url"]).content.count(b"HETATM") == 3
    assert client.get(result["files"]["ligand_sdf"]["download_url"]).status_code == 200
    assert database.list_molecules(data["admin"])["total"] == before == 2
    shared = client.post(f"/api/v1/molecules/{identifier}/bound-inputs", json={}, headers=headers(data, "student"))
    assert shared.status_code == 201, shared.text
    assert client.get(shared.json()["files"]["ligand_sdf"]["download_url"], headers=headers(data, "student")).status_code == 200
    assert client.post(f"/api/v1/molecules/{identifier}/bound-inputs", json={}, headers=headers(data, "other")).status_code == 404


def test_mcp_uses_rest_and_preserves_user_permissions(api_data):
    pytest.importorskip("mcp")
    from mcp import Client
    from mn_molecular_design_studio.api_client import StudioClient
    from mn_molecular_design_studio.mcp_server import create_server
    studio = StudioClient("http://studio.test", token=api_data["student_token"], transport=httpx.ASGITransport(app=api_data["app"]))

    async def check():
        async with Client(create_server(studio)) as client:
            tools = (await client.list_tools()).tools
            assert {"workflow_schema", "screening_results", "prepare_bound_pose", "run_pipeline"}.issubset({t.name for t in tools})
            schema = await client.call_tool("workflow_schema", {"workflow_id": "ligand.docking.unidock_pro"})
            assert not schema.is_error, schema.content
            assert "parameters_json_schema" in schema.structured_content
            identity = await client.call_tool("studio_identity", {})
            assert identity.structured_content["role"] == "user"
            rejected = await client.call_tool("studio_api_request", {"method": "GET", "path": "/api/v1/admin/users"})
            assert rejected.is_error
            blocked = await client.call_tool("studio_api_request", {"method": "POST", "path": "/api/v1/cluster/node/claim"})
            assert blocked.is_error
            staged = await client.call_tool("stage_text_input", {"filename": "ligands.smi", "text": "CCO ethanol\n"})
            assert staged.structured_content["reference"].startswith("studio-file:///")
    asyncio.run(check())


def test_mcp_http_never_falls_back_to_admin_or_development_identity(api_data, monkeypatch):
    pytest.importorskip("mcp")
    from mn_molecular_design_studio.api_client import StudioClient
    from mn_molecular_design_studio.mcp_server import http_app
    studio = StudioClient("http://studio.test", token=api_data["admin_token"], transport=httpx.ASGITransport(app=api_data["app"]))
    with TestClient(http_app(client=studio)) as client:
        assert client.post("/mcp", json={}).status_code == 401
        assert client.post("/mcp", json={}, headers={"Authorization": "Bearer invalid"}).status_code == 401
        monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE", "development")
        assert client.post("/mcp", json={}, headers=headers(api_data, "student")).status_code == 401


def test_mcp_http_forwards_each_users_identity_without_shared_credentials(api_data):
    pytest.importorskip("mcp")
    from mn_molecular_design_studio.api_client import StudioClient
    from mn_molecular_design_studio.mcp_server import http_app
    studio = StudioClient("http://studio.test", token=api_data["admin_token"], transport=httpx.ASGITransport(app=api_data["app"]))
    with TestClient(http_app(client=studio), base_url="http://localhost") as client:
        for role in ("student", "admin"):
            request_headers = {**headers(api_data, role), "Accept": "application/json, text/event-stream", "MCP-Protocol-Version": "2025-11-25"}
            response = client.post("/mcp", headers=request_headers, json={"jsonrpc": "2.0", "id": 1,
                "method": "tools/call", "params": {"name": "studio_identity", "arguments": {}}})
            assert response.status_code == 200, response.text
            messages = [json.loads(line[6:]) for line in response.text.splitlines() if line.startswith("data: ")]
            assert messages, response.text
            result = messages[-1]["result"]
            assert not result.get("isError"), result
            assert result["structuredContent"]["role"] == ("user" if role == "student" else "admin")


def test_workflow_json_schemas_are_valid_and_never_evaluate_type_text():
    from jsonschema import Draft202012Validator
    from mn_molecular_design_studio.workflow_json_schema import type_schema, parameters_schema
    assert type_schema("list[pathlib.Path] | None") == {"anyOf": [{"type": "array", "items": {"type": "string"}}, {"type": "null"}]}
    assert type_schema("__import__('os').system('must-not-run')") == {}
    schema = parameters_schema([{"name": "counts", "type": "tuple[int, int]", "required": True}])
    Draft202012Validator.check_schema(schema)
    Draft202012Validator(schema).validate({"counts": [1, 2]})
