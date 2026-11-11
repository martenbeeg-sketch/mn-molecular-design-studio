from pathlib import Path

import pytest

from test_api_integrations import api_data, headers


@pytest.mark.parametrize("role", ["admin", "student"])
def test_default_submission_stages_in_private_state_not_deployment(api_data, monkeypatch, role):
    from mn_molecular_design_studio import api

    deployment = api_data["root"] / "immutable-deployment"
    deployment.mkdir()
    monkeypatch.chdir(deployment)
    captured = []
    monkeypatch.setattr(api, "_cluster_submission", lambda *args, **kwargs: None)
    monkeypatch.setattr(api, "_guard_source_job", lambda *args, **kwargs: None)
    monkeypatch.setattr(api, "submit_workflow_request",
                        lambda path: captured.append(Path(path)) or {"job": {"run_id": "one-private-request"}})
    payload = {"schema_version": 1, "workflow": "analysis.ipsae", "parameters": {},
               "project_id": api_data["project"]["id"]}
    response = api_data["client"].post("/api/v1/workflows", json=payload, headers=headers(api_data, role))
    assert response.status_code == 202, response.text
    expected = api_data["root"] / "studio/users" / api_data[role]["id"] / "requests"
    assert captured[0].parent == expected
    assert expected.stat().st_mode & 0o777 == 0o700
    assert not captured[0].exists()
    assert not list(deployment.iterdir())


def test_explicit_server_directory_remains_admin_only(api_data):
    from mn_molecular_design_studio import api

    assert api._request_directory(api_data["admin"], str(api_data["root"])) == api_data["root"]
    with pytest.raises(api.HTTPException) as denied:
        api._request_directory(api_data["student"], str(api_data["root"]))
    assert denied.value.status_code == 403


def test_students_sign_in_separately_and_cannot_read_private_exports(api_data):
    client = api_data["client"]
    tokens = {}
    for username, password in [("student", "student-password-42"), ("otherstudent", "otherstudent-password-42")]:
        response = client.post("/api/v1/auth/login", json={"username": username, "password": password})
        assert response.status_code == 200, response.text
        tokens[username] = response.json()["access_token"]
    assert tokens["student"] != tokens["otherstudent"]
    first = {"Authorization": "Bearer " + tokens["student"]}
    second = {"Authorization": "Bearer " + tokens["otherstudent"]}
    imported = client.post("/api/v1/molecules/imports", headers=first,
                           json={"method": "smiles", "text": "CCO teaching_ethanol"})
    assert imported.status_code == 201, imported.text
    molecule = imported.json()["molecule_ids"][0]
    exported = client.post("/api/v1/molecules/exports", headers=first,
                           json={"molecule_ids": [molecule], "format": "CSV"})
    assert exported.status_code == 201, exported.text
    download = exported.json()["download_url"]
    assert client.get(download, headers=first).status_code == 200
    assert client.get(download, headers=second).status_code == 404
    assert client.get("/api/v1/molecules/" + molecule, headers=second).status_code == 404
    assert client.get("/api/v1/admin/users", headers=first).status_code == 403
    assert client.get("/api/v1/molecules/groups", headers={"Authorization": ""}).status_code == 401


def test_student_cannot_choose_a_server_input_directory(api_data):
    response = api_data["client"].post("/api/v1/workflows", params={"base_dir": str(api_data["root"])},
        headers=headers(api_data, "student"), json={
        "schema_version": 1, "workflow": "analysis.ipsae", "project_id": api_data["project"]["id"],
        "parameters": {},
    })
    assert response.status_code == 403


@pytest.mark.parametrize("state,authentication", [(Path("/production"), "production"), (Path("/state"), "development")])
def test_teaching_bootstrap_refuses_production_paths_or_shared_authentication(monkeypatch, state, authentication):
    from scripts.teaching_sandbox import seed

    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", "/state/home")
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE", authentication)
    with pytest.raises(RuntimeError):
        seed(state)
