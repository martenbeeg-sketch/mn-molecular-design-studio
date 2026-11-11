"""AUD-011: indexed source-qualified lookup must not inspect Results on render."""
import json
import time

import pytest

from mn_molecular_design_studio import job_catalog

USER = {"id": "a" * 32, "role": "user", "_cluster_project_ids": []}
SOURCE = {"molecule_id": "m", "structure_id": "s", "group_id": "g", "structure_digest": "a" * 64, "model_id": "1"}


@pytest.fixture
def indexed_jobs(tmp_path, monkeypatch):
    root = tmp_path / "runs"
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(root))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path / "app"))
    monkeypatch.setattr(job_catalog, "catalog_path", lambda: tmp_path / "catalog.sqlite3")
    monkeypatch.setattr(job_catalog, "_SCHEMA_READY_PATH", None)
    monkeypatch.setattr(job_catalog, "_reconcile_for_read", lambda root: None)
    run = root / "analysis/old-matching-job"
    (run / "artifacts/analysis").mkdir(parents=True)
    source = {**SOURCE, "candidate_id": "candidate"}
    for path, value in {
        "metadata.json": {"tool": "biotite_interface", "status": "collecting", "studio_owner_id": USER["id"],
                          "upstream_source_kind": "database", "database_source": source, "created_at": "2020-01-01"},
        "input.json": {"params": {"contact_distance": 5.0}},
        "artifacts/analysis/batch_manifest.json": {"selection": {"model_id": "1", "binder_chains": ["B"], "target_chains": ["A"]},
            "settings": {"contact_distance": 5.0}, "items": [{"id": "item-1", "candidate_id": "candidate",
                "structure_digest": SOURCE["structure_digest"], "binder_chains": ["B"], "target_chains": ["A"]}]},
    }.items():
        (run / path).write_text(json.dumps(value))
    job_catalog._upsert(root, run, refresh_artifacts=False)
    with job_catalog._connection() as connection:
        columns = [r[1] for r in connection.execute("PRAGMA table_info(jobs)") if r[1] != "id"]
        row = dict(connection.execute("SELECT * FROM jobs").fetchone())
        for index in range(1005):
            copy = {**row, "run_id": f"unrelated-{index:04}", "run_relative_path": f"analysis/unrelated-{index:04}", "created_at": "2030-01-01"}
            identity = json.loads(copy["run_identity_json"])
            if identity.get("database_analysis"):
                identity["database_analysis"]["source_key"] = "unrelated"
            copy["run_identity_json"] = json.dumps(identity)
            connection.execute("INSERT INTO jobs (" + ",".join(columns) + ") VALUES (" + ",".join("?" for _ in columns) + ")",
                               [copy[key] for key in columns])
    return root


def lookup(*, load_evidence=True):
    from mn_molecular_design_studio.complex_analysis_ui import _matching_database_analysis_jobs
    return _matching_database_analysis_jobs(USER, workflow_id="analysis.biotite_interface", source=SOURCE,
        candidate_id="candidate", model_id="1", binder_chains=["B"], target_chains=["A"], parameters={"contact_distance": 5.0},
        load_evidence=load_evidence)


def test_older_exact_source_is_found_beyond_first_thousand_jobs(indexed_jobs):
    assert [row["run_id"] for row in lookup()] == ["old-matching-job"]


def test_render_lookup_performs_no_result_file_reads_or_reconciliation(indexed_jobs, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_ui
    def forbidden(*args, **kwargs):
        pytest.fail("Results/reconciliation accessed during indexed coverage lookup")
    monkeypatch.setattr(job_catalog, "_read_json", forbidden)
    monkeypatch.setattr(job_catalog, "_reconcile_for_read", forbidden)
    monkeypatch.setattr(complex_analysis_ui, "_read_json", forbidden)
    timings = []
    for _ in range(3):
        started = time.perf_counter()
        assert len(lookup(load_evidence=False)) == 1
        timings.append(time.perf_counter() - started)
    print("indexed cold/warm seconds", timings)


def test_indexed_lookup_revalidates_owner_acl(indexed_jobs):
    assert job_catalog.database_analysis_jobs({"id": "b" * 32, "role": "user"}, [SOURCE]) == []


def test_uninitialized_or_partial_index_is_not_reported_as_absence(indexed_jobs):
    with job_catalog._connection() as connection:
        connection.execute("UPDATE jobs SET run_identity_json='{}' WHERE run_id='old-matching-job'")
    with pytest.raises(RuntimeError, match="index"):
        job_catalog.database_analysis_jobs(USER, [SOURCE])


def test_exact_source_history_is_paged_without_a_first_thousand_limit(indexed_jobs, monkeypatch):
    from mn_molecular_design_studio import complex_analysis_ui
    with job_catalog._connection() as connection:
        original = json.loads(connection.execute("SELECT run_identity_json FROM jobs WHERE run_id='old-matching-job'").fetchone()[0])
        for row in connection.execute("SELECT id,run_identity_json FROM jobs WHERE run_id<>'old-matching-job'").fetchall():
            value = json.loads(row["run_identity_json"])
            value["database_analysis"] = {**original["database_analysis"], "parameters": {"contact_distance": 4.0}}
            connection.execute("UPDATE jobs SET run_identity_json=? WHERE id=?", (json.dumps(value), row["id"]))
    monkeypatch.setattr(complex_analysis_ui, "_read_json", lambda *a: pytest.fail("coverage accessed Results"))
    assert [row["run_id"] for row in lookup(load_evidence=False)] == ["old-matching-job"]


def test_unknown_terminal_protocol_blocks_until_ingestion_or_explicit_reassessment(indexed_jobs):
    from mn_molecular_design_studio.complex_analysis_ui import _matching_database_analysis_jobs
    with job_catalog._connection() as connection:
        connection.execute("UPDATE jobs SET status='completed' WHERE run_id='old-matching-job'")
    request = dict(workflow_id="analysis.biotite_interface", source=SOURCE, candidate_id="candidate", model_id="1",
                   binder_chains=["B"], target_chains=["A"], parameters={"contact_distance": 5.0},
                   load_evidence=False, expected_protocol={"method_version": "v2"})
    assert len(_matching_database_analysis_jobs(USER, **request)) == 1
    assert _matching_database_analysis_jobs(USER, **request, allow_reassessment=True) == []


@pytest.mark.parametrize("method,tool,extra", [
    ("analysis.ipsae", "ipsae", {"pae_digest": "known-in-execution"}),
    ("analysis.sequence_liability", "sequence_liability", {"ligand_selection": {}}),
])
def test_active_lookup_is_conservative_when_optional_request_evidence_is_not_loaded(indexed_jobs, method, tool, extra):
    from mn_molecular_design_studio.complex_analysis_ui import _matching_database_analysis_jobs
    with job_catalog._connection() as connection:
        row = connection.execute("SELECT id,run_identity_json FROM jobs WHERE run_id='old-matching-job'").fetchone()
        identity = json.loads(row["run_identity_json"])
        identity["database_analysis"]["source"].update(extra)
        connection.execute("UPDATE jobs SET tool=?,run_identity_json=? WHERE id=?", (tool, json.dumps(identity), row["id"]))
    args = dict(workflow_id=method, source=SOURCE, candidate_id="candidate", model_id="1", binder_chains=["B"],
                target_chains=["A"], parameters={"contact_distance": 5.0}, load_evidence=False)
    assert len(_matching_database_analysis_jobs(USER, **args)) == 1
    if tool == "ipsae":
        args["source"] = {**SOURCE, "pae_digest": "different-explicit-pae"}
        assert _matching_database_analysis_jobs(USER, **args) == []
