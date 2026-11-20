from hashlib import sha256

from test_api_integrations import api_data, headers


def test_native_failure_and_configured_final_decision_are_independent():
    from mn_molecular_design_studio.binder_result_evidence import binder_decisions

    native, final = binder_decisions({
        "metrics": {"native_pass_filters": False},
        "raw_metadata": {"validation_acceptance": {"status": "accepted"}},
    })
    assert native["status"] == "failed"
    assert final["status"] == "accepted"
    native, final = binder_decisions({
        "metrics": {"ppi.pass": True},
        "raw_metadata": {"distributed_campaign": {"native_acceptance": "accepted",
                                                  "validated_acceptance": "pending",
                                                  "validation_required": True}},
    })
    assert native["status"] == "unavailable"
    assert final["status"] == "pending"


def test_results_and_mcp_only_link_registered_exact_models(api_data):
    import asyncio

    import httpx
    from mcp import Client
    from mn_molecular_design_studio.api_client import StudioClient
    from mn_molecular_design_studio.binder_result_evidence import binder_outcomes
    from mn_molecular_design_studio.binder_results_ui import _assessment_choices
    from mn_molecular_design_studio.engine.workflow_engine.core import candidates, jobs
    from mn_molecular_design_studio.mcp_server import create_server
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    from mn_molecular_design_studio.workspace import create_project, mark_studio_job

    data = api_data
    project, admin = data["project"]["id"], data["admin"]
    job = jobs.create_job("design", "binder_design", "bindcraft2", {})
    mark_studio_job({"job": {"run_dir": str(job.run_dir)}}, project_id=project, owner_id=admin["id"])
    shape = (b"ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n"
             b"ATOM      2  CA  ALA B   1       2.000   0.000   0.000  1.00 20.00           C\nEND\n")
    original = job.run_dir / "artifacts" / "bound.pdb"
    original.write_bytes(shape)
    digest = sha256(shape).hexdigest()
    changed = job.run_dir / "artifacts" / "other.pdb"
    changed.write_bytes(shape.replace(b"2.000", b"4.000"))
    rows = candidates.write_candidates(job.run_dir, "bindcraft2", [
        {"candidate_id": candidate_id, "stage": "generation.backbone_sequence", "source_tool": "bindcraft2",
         "binder_chains": ["B"], "target_chains": ["A"], "binder_sequence": "A", "complex_pdb": str(path),
         "metrics": {"native_pass_filters": False},
         "raw_metadata": {"validation_acceptance": {"status": "accepted"}}}
        for candidate_id, path in (("exact", original), ("alternate", changed), ("elsewhere", original))
    ])
    db = MolecularDatabase()
    design_group = db.create_group(admin, "Accepted", project_id=project, origin={"collection_kind": "curated"})
    other_project = create_project("Private unrelated project", owner_id=admin["id"])["id"]
    other_group = db.create_group(admin, "Other project", project_id=other_project, origin={"collection_kind": "curated"})
    inference_group = db.create_group(admin, "Only analysis", project_id=project, origin={"collection_kind": "evaluation"})
    for candidate_id, group in (("exact", design_group), ("alternate", design_group), ("elsewhere", other_group),
                                ("elsewhere", inference_group)):
        result = db.import_records(admin, [{"name": "Binder", "components": [{"type": "protein", "chain": "B", "sequence": "A"}]}],
                                   group_id=group, project_id=other_project if group == other_group else project)
        molecule = result["molecule_ids"][0]
        structure = db.attach_structure(group, molecule, {"name": "bound.pdb", "content": shape},
                                       context={"run_id": job.run_id, "candidate_id": candidate_id})
        if candidate_id == "exact":
            matching = db.add_observation(group, molecule, "binderscope", "verified-report", {"contact_count": 1},
                context={"representation_id": structure, "source_structure_digest": digest, "binder_chains": ["B"],
                         "target_chains": ["A"], "candidate": {"raw_metadata": {"binderscope": {
                             "source": {"input_structure_digest": digest}}}}})
            saved_refold = changed.read_bytes()
            refold_digest = sha256(saved_refold).hexdigest()
            refold = db.attach_structure(group, molecule, {"name": "refold.pdb", "content": saved_refold,
                                        "label": "Refolded validation complex"},
                                        context={"run_id": job.run_id, "candidate_id": candidate_id,
                                                 "source_structure_digest": digest, "design_representation_role": "refolded_validation",
                                                 "binder_chains": ["B"],
                                                 "target_chains": ["A"]})
            derived = db.add_observation(group, molecule, "binderscope", "refold-report", {"contact_count": 1},
                context={"representation_id": refold, "source_structure_digest": refold_digest,
                         "binder_chains": ["B"], "target_chains": ["A"],
                         "candidate": {"raw_metadata": {"binderscope": {
                             "source": {"input_structure_digest": refold_digest}}}}})
            # An assessment aimed at another digest must stay out of Results.
            db.add_observation(group, molecule, "binderscope", "wrong-report", {"contact_count": 2},
                context={"representation_id": structure, "source_structure_digest": digest,
                         "binder_chains": ["B"], "target_chains": ["A"],
                         "candidate": {"raw_metadata": {"binderscope": {"source": {"input_structure_digest": "f" * 64}}}}})
    jobs.finish_job(job.run_dir, True, {"outputs": {}})

    reference = f"design/{job.run_id}"
    headers_student = headers(data, "student")
    direct = data["client"].get(f"/api/v1/jobs/local:design:{job.run_id}/presentation", headers=headers_student)
    assert direct.status_code == 200, direct.text
    observed = {item["candidate_id"]: item["binder_outcome"] for item in direct.json()["items"]}
    assert observed["exact"]["native_filter"]["status"] == "failed"
    assert observed["exact"]["final_acceptance"]["status"] == "accepted"
    assert observed["exact"]["database"]["status"] == "registered"
    assert [a["observation_id"] for a in observed["exact"]["binderscope"]] == [matching]
    assert [a["observation_id"] for a in observed["exact"]["binderscope_derived"]] == [derived]
    assert observed["exact"]["binderscope_derived"][0]["structure_id"] == refold
    assert observed["exact"]["preferred_binderscope_observation_id"] == derived
    assert [(a["observation_id"], a["relation"]) for a in _assessment_choices(observed["exact"])] == [
        (derived, "Refolded model"), (matching, "Results model")]
    assert "refolded_validation" == observed["exact"]["database"]["derived_structures"][0]["representation_role"]
    assert observed["alternate"]["database"]["status"] == "different_saved_representation"
    assert observed["alternate"]["binderscope"] == observed["alternate"]["binderscope_derived"] == []
    assert observed["elsewhere"]["database"]["status"] == "not_registered"
    assert data["client"].get(f"/api/v1/jobs/local:design:{job.run_id}/presentation",
                              headers=headers(data, "other")).status_code == 404

    remote = binder_outcomes(db, data["student"], {"run_id": job.run_id, "project_id": project}, rows,
                             digests={"exact": digest})
    assert remote[0] == observed["exact"]
    assert remote[1]["database"]["status"] == "structure_unverified"

    async def journey():
        connection = StudioClient("http://studio.test", token=data["student_token"],
                                  transport=httpx.ASGITransport(app=data["app"]))
        async with Client(create_server(connection)) as mcp:
            report = await mcp.call_tool("result_presentation", {"reference": reference, "limit": 1})
            assert not report.is_error, report.content
            assert report.structured_content["items"][0]["binder_outcome"] == observed["exact"]

    asyncio.run(journey())

    with db.connection() as connection:
        primary_id = observed["exact"]["database"]["structures"][0]["structure_id"]
        connection.execute("DELETE FROM molecule_structures WHERE id=?", (primary_id,))
        connection.commit()
    only_refold = binder_outcomes(db, data["student"], {"run_id": job.run_id, "project_id": project}, rows,
                                  digests={"exact": digest})[0]
    assert only_refold["database"]["status"] == "registered"
    assert only_refold["database"]["structures"] == []
    assert only_refold["preferred_binderscope_observation_id"] == derived


def test_binder_results_panel_renders_a_bounded_page_and_empty_run(api_data):
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string("""
from mn_molecular_design_studio.binder_results_ui import render_binder_outcomes
job = {"tool": "bindcraft2", "run_id": "fixture"}
record = {"candidate_id": "sample-01", "source_tool": "bindcraft2",
          "metrics": {"native_pass_filters": False}}
user = {"id": "fixture-owner", "role": "student"}
render_binder_outcomes(job, [record], user, digests={}, key="sample")
render_binder_outcomes(job, [], user, digests={}, key="empty")
""").run()
    assert not app.exception
    assert len(app.dataframe) == 1
    assert app.dataframe[0].value.iloc[0].to_dict() == {
        "Design": "sample-01", "Native filter": "Failed", "Final acceptance": "Not evaluated",
        "Database": "Structure unverified", "BinderScope": "—",
    }
    assert "No candidate records" in app.info[0].value
