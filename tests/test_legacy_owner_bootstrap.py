"""Local bootstrap must preserve research identities without granting students access."""
import concurrent.futures

import pytest

from mn_molecular_design_studio import auth, workspace
from mn_molecular_design_studio.molecular_database import MolecularDatabase


@pytest.fixture
def legacy_workspace(tmp_path, monkeypatch):
    for key in ("MN_MOLECULAR_DESIGN_STUDIO_HOME", "MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", "MDS_CLUSTER_HOME"):
        monkeypatch.setenv(key, str(tmp_path / key.lower()))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_AUTH_MODE", "production")
    owner = auth.development_user()
    project = workspace.create_project("Important research", owner_id=owner["id"])
    database = MolecularDatabase()
    imported = database.import_records(owner, [{"name": "Protected import", "components": [
        {"type": "small_molecule", "chain": "L", "smiles": "CCO"}]}], name="Important library", project_id=project["id"])
    return database, imported, project


def test_bootstrap_preserves_owner_and_private_data_without_rewriting_science(legacy_workspace):
    database, imported, project = legacy_workspace
    with database.connection() as connection:
        before = list(map(tuple, connection.execute("SELECT * FROM molecule_groups ORDER BY id")))
    instructor = auth.create_initial_admin("instructor", "Research instructor", "a-long-temporary-password",
                                           preserve_development_owner=True)
    assert instructor["id"] == auth.development_user()["id"]
    assert instructor["must_change_password"]
    assert not instructor.get("development_mode")
    assert auth.authenticate("instructor", "a-long-temporary-password")[0]["id"] == instructor["id"]
    assert imported["group_id"] in database._authorize_group_ids_for_query(instructor, [imported["group_id"]])
    student = auth.create_user("student", "Student", "a-student-temporary-password", actor=instructor)
    with pytest.raises(PermissionError):
        database._authorize_group_ids_for_query(student, [imported["group_id"]])
    with database.connection() as connection:
        assert list(map(tuple, connection.execute("SELECT * FROM molecule_groups ORDER BY id"))) == before
    assert project["id"] in auth.accessible_project_ids(instructor)
    assert project["id"] not in auth.accessible_project_ids(student)
    with pytest.raises(ValueError, match="already exists"):
        auth.create_initial_admin("another", "Other", "another-temporary-password", preserve_development_owner=True)


def test_default_bootstrap_does_not_implicitly_adopt_research(legacy_workspace):
    instructor = auth.create_initial_admin("instructor", "Instructor", "another-temporary-password")
    assert instructor["id"] != auth.development_user()["id"]
    assert not instructor["must_change_password"]


def test_competing_bootstraps_create_only_one_instructor(legacy_workspace):
    def create(index):
        try:
            return auth.create_initial_admin(f"instructor{index}", "Instructor", "another-temporary-password",
                                              preserve_development_owner=True)["id"]
        except ValueError:
            return None
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(create, range(2)))
    assert results.count(auth.development_user()["id"]) == 1
    assert results.count(None) == 1
    assert auth.user_count() == 1


@pytest.mark.parametrize("flag", ["true", 1, None])
def test_bootstrap_requires_explicit_boolean(legacy_workspace, flag):
    with pytest.raises(ValueError):
        auth.create_initial_admin("instructor", "Instructor", "another-temporary-password", preserve_development_owner=flag)
    assert auth.user_count() == 0
