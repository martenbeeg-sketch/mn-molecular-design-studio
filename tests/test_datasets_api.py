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
