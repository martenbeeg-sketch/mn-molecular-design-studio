"""AUD-006: matching in-flight/unknown jobs never mean missing evidence."""
import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio.complex_analysis_batch import analysis_job_coverage


@pytest.mark.parametrize("status", ["queued", "running", "collecting", "transferring", "finalizing",
                                    "staging", "paused", "unknown", "future-new-state"])
def test_all_nonterminal_matching_states_block_submission(status):
    assert analysis_job_coverage([{"status": status, "run_id": "matching-job"}])[0] == "Queued / running"


SCRIPT = '''
import streamlit as st
from mn_molecular_design_studio.complex_analysis_collection_ui import _render_batch_assessment
from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema
sample = {
  "source_kind": "Database", "identity": "model", "molecule_id": "m", "structure_id": "s",
  "model_id": "1", "digest": "a" * 64, "kind": "protein-protein", "label": "Saved design",
  "structure": {"group_id": "g"}, "candidate": {"binder_chains": ["B"], "target_chains": ["A"]},
  "components": [{"type": "protein", "source_chain": "B", "sequence": "A" * 30},
                 {"type": "protein", "source_chain": "A", "sequence": "G" * 30}], "associated": []}
def submit(*args, **kwargs):
    st.session_state["submissions"] = st.session_state.get("submissions", 0) + 1
_render_batch_assessment({"id": "a" * 32}, [sample], workflow_schema=workflow_schema, submit=submit)
'''


@pytest.mark.parametrize("state", ["collecting", "transferring", "finalizing", "unknown", "completed", "catalog-error"])
def test_real_planning_ui_does_not_submit_matching_nonterminal_or_unverified_work(monkeypatch, state):
    from mn_molecular_design_studio import complex_analysis_ui, molecular_database
    class Database:
        def group(self, *args, **kwargs):
            return {"project_id": "", "origin": {}}
    monkeypatch.setattr(molecular_database, "MolecularDatabase", Database)
    def matching(*args, **kwargs):
        if state == "catalog-error":
            raise RuntimeError("catalog temporarily unavailable")
        return [{"status": state, "run_id": "existing-job"}]
    monkeypatch.setattr(complex_analysis_ui, "_matching_database_analysis_jobs", matching)
    app = AppTest.from_string(SCRIPT).run()
    assert not app.exception
    assert app.button(key="complex-analysis-submit-batch-assessment").disabled
    assert "submissions" not in app.session_state
