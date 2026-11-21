from streamlit.testing.v1 import AppTest


def test_followup_reuses_key_after_timeout_and_opens_accepted_job():
    app = AppTest.from_string("""
import streamlit as st
from mn_molecular_design_studio import results_followup_ui as ui
if not st.session_state.get("results-followup"):
    st.session_state["results-followup"] = {"kind": "dock", "job": {}, "selection": {}}

def render_docking(job, selection, user, schema, submit):
    if st.button("Queue"):
        try:
            submit({}, "analysis.plip", {"candidate_ids": ["exact-pose"]}, "project",
                   {"_execution": {"placement": "auto"}}, "PLIP")
        except ValueError:
            st.warning("Submission reply timed out")

def submit(schema, workflow, params, project, settings, label):
    requests = st.session_state.setdefault("requests", [])
    requests.append(settings["_execution"])
    if len(requests) == 1:
        raise ValueError("Coordinator timed out")
    return {"job": {"cluster_job_id": "analysis-id", "job_code": "JOB1"}}

if st.session_state.get("page") != "Job details":
    original = ui._render_docking
    try:
        ui._render_docking = render_docking
        ui.render_followup({"id": "user"}, workflow_schema=None, project_selector=None, submit=submit)
    finally:
        ui._render_docking = original
""").run()
    for _ in range(3):
        app = next(button for button in app.button if button.label == "Queue").click().run()
        assert not app.exception
    first, retry = app.session_state["requests"]
    assert first == retry and first["idempotency_key"]
    app = app.button(key="followup-open-analysis-job").click().run()
    assert app.session_state["selected-job-reference"] == "cluster:analysis-id"
