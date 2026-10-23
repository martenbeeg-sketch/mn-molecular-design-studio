import pytest
from streamlit.testing.v1 import AppTest

from mn_molecular_design_studio import library_centroids_ui as ui


@pytest.mark.parametrize("method", ["sampled", "hierarchical", "exact"])
def test_centroid_tile_submits_worker_and_chunk_settings(monkeypatch, method):
    nodes = [{"node_id": "source", "name": "source", "local": True, "cpu_available": 4},
             {"node_id": "remote", "name": "remote", "local": False, "cpu_available": 3}]
    monkeypatch.setattr(ui, "_worker_choices", lambda: nodes)
    monkeypatch.setattr(ui.cluster_config, "identity", lambda: "source")
    monkeypatch.setattr(ui.cluster_config, "load", lambda: {"default_placement": "local"})
    monkeypatch.setattr(ui.cluster_config, "enabled", lambda: True)
    record = {"id": "library", "name": "test", "owner_node_id": "source", "preparations": [
        {"schema_version": ui.FEATURE_SCHEMA, "owner_node_id": "source", "rows": 100_000_000,
         "path": "/external/manifest.json"}]}
    monkeypatch.setattr(ui.LibraryCatalog, "list", lambda self: [record])
    script = '''
import streamlit as st
from mn_molecular_design_studio.library_centroids_ui import render_setup
def capture(*args):
    st.session_state["captured"] = args
render_setup(None, user={"role": "admin"}, project_selector=lambda key: "project", submit=capture)
'''
    app = AppTest.from_string(script).run()
    assert not app.exception
    app.selectbox(key="library-centroid-method").select(method).run()
    if method != "exact":
        assert app.number_input(key="library-centroid-cpus-v1-source").value == 4
        assert app.number_input(key="library-centroid-cpus-v1-remote").value == 3
        app.radio(key="library-centroid-sizing").set_value("Number of chunks").run()
        app.number_input(key="library-centroid-chunk-count").set_value(324).run()
    app.text_input(key="library-centroid-output").set_value("/external/new-run").run()
    app.button(key="library-centroid-submit").click().run()
    assert not app.exception
    args = app.session_state["captured"]
    assert args[1] == "library.centroids"
    parameters = args[2]
    assert parameters["clustering_method"] == method
    if method != "exact":
        assert parameters["worker_nodes"] == [{"node_id": "source", "cpu_cores": 4},
                                               {"node_id": "remote", "cpu_cores": 3}]
        assert parameters["chunk_count"] == 324
        assert parameters["cpu_cores"] == 4
        assert parameters["final_candidate_count"] == 300_000


def test_centroid_worker_dialog_describes_clustering():
    from mn_molecular_design_studio.cluster.ui import _library_worker_details

    worker = {"node_id": "node", "name": "computer", "task_kind": "centroids",
              "cpu_worker_count": 4, "status": "completed", "row_count": 10,
              "chunk_count": 3, "input_rows": 100, "parallel_chunks": 2, "cores_per_chunk": 2}
    row = {"campaign_name": "centroids", "status": "completed", "metrics": {"workers": [worker]}}
    details = _library_worker_details(row, [row])
    assert details["task_kind"] == "centroids"
    app = AppTest.from_string('''
import streamlit as st
from mn_molecular_design_studio.cluster.ui import _render_library_worker_details
st.session_state['cluster-results-library-worker-selected'] = 'node'
_render_library_worker_details({
    'task_kind': 'centroids', 'job_name': 'centroids', 'job_status': 'Completed',
    'workers': [{'node_id':'node','name':'computer','completed':True,
                'cpu_worker_count':4,'chunk_count':3,'row_count':10}]
})
''').run()
    assert not app.exception
    assert "Cluster on CPUs" in set(app.dataframe[0].value["Step"])
    assert "Calculate Morgan fingerprints" not in set(app.dataframe[0].value["Step"])
