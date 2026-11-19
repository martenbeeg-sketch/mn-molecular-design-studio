from mn_molecular_design_studio import auth, datasets
from test_binderscope_api import preserved_api


def headers(case, *, other=False):
    return {"Authorization": "Bearer " + case["other_token" if other else "token"]}


def test_authenticated_import_creates_project_independent_inputs(preserved_api):
    case = preserved_api
    response = case["client"].post("/api/v1/molecules/imports", headers=headers(case), json={
        "name": "Reusable compounds", "method": "smiles", "text": "CCO ethanol\nCO methanol", "reusable_dataset": True,
        "project_id": case["payload"]["project_id"]})
    assert response.status_code == 201, response.text
    group = case["database"].group(response.json()["group_id"], case["actor"])
    assert group["project_id"] == ""
    assert datasets.is_input_dataset(group)
    assert datasets.project_references(group) == [case["payload"]["project_id"]]
    listing = case["client"].get("/api/v1/molecules/datasets", headers=headers(case), params={"project_id": case["payload"]["project_id"]})
    assert listing.status_code == 200
    assert [row["id"] for row in listing.json()["items"]] == [group["id"]]
    assert case["pending"] == []


def test_dataset_references_do_not_grant_access_and_stale_edits_conflict(preserved_api):
    case = preserved_api
    project_id = case["payload"]["project_id"]
    response = case["client"].post("/api/v1/molecules/datasets", headers=headers(case), json={
        "name": "General sequences", "project_references": [project_id]})
    assert response.status_code == 201
    group = response.json()
    other = auth.resolve_session(case["other_token"])
    auth.grant_project_access(project_id, other["id"], actor=case["actor"])
    assert case["client"].get("/api/v1/molecules/datasets", headers=headers(case, other=True)).json()["items"] == []
    payload = {"project_references": [], "expected_revision": datasets.source_revision(group)}
    address = "/api/v1/molecules/datasets/" + group["id"] + "/projects"
    assert case["client"].patch(address, json=payload).status_code == 401
    assert case["client"].patch(address, headers=headers(case, other=True), json=payload).status_code == 404
    assert case["client"].patch(address, headers=headers(case), json=payload).status_code == 200
    assert case["client"].patch(address, headers=headers(case), json=payload).status_code == 409


def test_dataset_import_rejects_unknown_projects_and_result_collection_destination(preserved_api):
    case = preserved_api
    base = {"name": "Inputs", "method": "smiles", "text": "CCO", "reusable_dataset": True}
    assert case["client"].post("/api/v1/molecules/imports", headers=headers(case), json={**base,
        "project_references": ["missing"]}).status_code == 422
    group_id = case["payload"]["binderscope_collection"]["sources"][0]["group_id"]
    assert case["client"].post("/api/v1/molecules/imports", headers=headers(case), json={**base,
        "group_id": group_id}).status_code == 422
    assert case["pending"] == []


def test_append_requires_explicit_project_edit_and_members_are_input_only(preserved_api):
    case = preserved_api
    actor, client, database = case["actor"], case["client"], case["database"]
    project = case["payload"]["project_id"]
    address = "/api/v1/molecules/imports"
    base = {"name": "Reusable test inputs", "method": "smiles", "reusable_dataset": True}
    initial = client.post(address, headers=headers(case), json={**base, "text": "CCO ethanol", "project_references": [project]})
    assert initial.status_code == 201, initial.text
    group, molecule = initial.json()["group_id"], initial.json()["molecule_ids"][0]
    member_url = f"/api/v1/molecules/datasets/{group}/members"
    append = {**base, "text": "CCN ethylamine", "group_id": group}
    assert client.post(address, headers=headers(case), json={**append, "project_references": [project]}).status_code == 422
    assert client.post(address, headers=headers(case), json={**append, "project_id": project}).status_code == 422
    assert client.get(member_url, headers=headers(case)).json()["total"] == 1
    response = client.post(address, headers=headers(case), json=append)
    assert response.status_code == 201, response.text
    assert datasets.project_references(database.group(group, actor)) == [project]

    result = database.import_records(actor, [{"name": "Measured fixture", "components": [
        {"type": "small_molecule", "chain": "L", "smiles": "CCO"}]}], project_id=project,
        origin={"collection_kind": "generated"})
    database.add_observation(result["group_id"], molecule, "fixture", "different-result", {"score": 7.})
    first = client.get(member_url, headers=headers(case), params={"limit": 1}).json()
    assert first["next_offset"] == 1 and first["total"] == 2
    second = client.get(member_url, headers=headers(case), params={"limit": 1, "offset": 1}).json()
    assert second["next_offset"] is None and first["items"][0]["id"] != second["items"][0]["id"]
    items = client.get(member_url, headers=headers(case), params={"query": "ethanol"}).json()["items"]
    assert len(items) == 1 and items[0]["id"] == molecule and items[0]["observations"] == []
    assert client.get(member_url, headers=headers(case, other=True)).status_code == 404
    assert client.get(f"/api/v1/molecules/datasets/{result['group_id']}/members", headers=headers(case)).status_code == 404
    generic = client.get("/api/v1/molecules", headers=headers(case), params={"group_id": group, "query": "ethanol"})
    assert len(generic.json()["items"][0]["observations"]) == 1
