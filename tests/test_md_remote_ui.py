from __future__ import annotations

import json

from mn_molecular_design_studio.md_remote_ui import render_remote_md_results


class FakeStreamlit:
    def __init__(self):
        self.calls = []
        self.session_state = {}

    def __getattr__(self, name):
        def record(*args, **kwargs):
            self.calls.append((name, args, kwargs))
            if name == "button":
                return False
            if name == "selectbox":
                return args[1][0]
            return None

        return record


def _manifest():
    return {"schema_version": 1, "studies": [{
        "study_id": "study-a", "compound_id": "CHEMBL-1", "compound_name": "ligand",
        "target_pdb_id": "1ABC", "target_receptor": "TR-A", "md_engine": "OpenMM",
        "selection_status": "manual_nonbest", "production_runs": [],
        "mmgbsa_runs": [
            {"run_id": "ep-1", "status": "completed", "result_path": "artifacts/legacy/one.json"},
            {"run_id": "ep-2", "status": "completed", "result_path": "artifacts/legacy/two.json"},
            {"run_id": "stage-only", "status": "completed"},
        ],
        "analysis_runs": [{"run_id": "analysis-a", "status": "completed", "result_path": "artifacts/md_analysis/summary.json"}],
    }]}


def test_remote_manifest_fetch_is_lazy_selected_and_does_not_request_trajectory_files():
    manifest = _manifest()
    summary = {
        "schema_version": 1, "success": True,
        "analysis": "RMSD and geometric contacts",
        "units": {"time": "ns", "distance": "angstrom", "rmsd": "angstrom"},
        "source": {"topology": {"sha256": "top-hash", "size_bytes": 10}, "trajectory": {"sha256": "traj-hash"}},
        "parameters": {"stride": 5, "contact_cutoff_angstrom": 4.5},
        "selection": {"ligand": "resname LIG", "protein": "protein"},
        "selected_frame_count": 3,
        "contact_map_semantics": "Geometric contact only.",
        "warnings": [],
    }
    endpoint = {"success": True, "mmgbsa": {"status": "success", "units": "kcal/mol", "gb": {"delta": {"delta_g_bind_total_kcal_mol": -4}}}}
    contents = {
        "artifacts/md/study.json": json.dumps(manifest),
        "artifacts/legacy/one.json": json.dumps(endpoint),
        "artifacts/md_analysis/summary.json": json.dumps(summary),
        "artifacts/md_analysis/timeseries.csv": "time_ns,protein_aligned_backbone_rmsd_angstrom,ligand_rmsd_angstrom\n0,0.2,0.4\n1,0.5,0.8\n",
        "artifacts/md_analysis/residue_contacts.csv": "residue,contact_occupancy_fraction,contact_definition\nA:42,0.75,geometric\n",
    }
    requested = []

    def fetch(job_id, path, token, state_key, label):
        requested.append((job_id, path, token, state_key, label))
        return contents.get(path)

    st = FakeStreamlit()
    render_remote_md_results(st, "job-a", {"tool": "md_study"}, "secret-token", fetch)

    paths = [row[1] for row in requested]
    assert paths == [
        "artifacts/md/study.json",
        "artifacts/legacy/one.json",
        "artifacts/md_analysis/summary.json",
        "artifacts/md_analysis/timeseries.csv",
        "artifacts/md_analysis/residue_contacts.csv",
    ]
    assert all(item[0] == "job-a" and item[2] == "secret-token" for item in requested)
    assert "top-hash" in json.dumps(st.calls, default=str)
    assert "RMSD over time" in json.dumps(st.calls, default=str)
    assert "Geometric ligand contact occupancy" in json.dumps(st.calls, default=str)
    chart = next(call for call in st.calls if call[0] == "line_chart")
    assert chart[1][0][0]["backbone_rmsd_angstrom"] == 0.2
    assert all("trajectory" not in path for path in paths)


def test_remote_paths_are_validated_before_authenticated_fetch():
    manifest = _manifest()
    manifest["studies"][0]["analysis_runs"] = [{
        "run_id": "analysis-a", "status": "completed", "result_path": "../escape.json",
    }]
    contents = {
        "artifacts/md/study.json": json.dumps(manifest),
        "artifacts/legacy/one.json": json.dumps({"success": True, "mmgbsa": {"status": "success", "units": "kcal/mol", "gb": {"delta": {"g_bind_total": -1}}}}),
    }
    requested = []

    def fetch(job_id, path, token, state_key, label):
        requested.append(path)
        return contents.get(path)

    st = FakeStreamlit()
    render_remote_md_results(st, "job-a", {}, "token", fetch)

    assert "../escape.json" not in requested
    assert any(call[0] == "warning" for call in st.calls)


def test_remote_refresh_clears_only_md_fetch_cache_keys():
    class RefreshStreamlit(FakeStreamlit):
        def __init__(self):
            super().__init__()
            self.session_state = {"remote-md-job-a-fetch-old": "cmd", "other-key": "keep"}

        def __getattr__(self, name):
            if name == "button":
                return lambda *args, **kwargs: True
            return super().__getattr__(name)

    st = RefreshStreamlit()
    render_remote_md_results(st, "job-a", {}, "token", lambda *args: None)

    assert "remote-md-job-a-fetch-old" not in st.session_state
    assert st.session_state["other-key"] == "keep"


def test_native_remote_md_analysis_fetches_summary_first_without_manifest_probe():
    summary = {
        "schema_version": 1, "success": True,
        "analysis": "RMSD and geometric contacts",
        "units": {"time": "ns", "distance": "angstrom", "rmsd": "angstrom"},
        "source": {"topology": {"sha256": "top"}, "trajectory": {"sha256": "traj"}},
        "parameters": {}, "selected_frame_count": 1,
        "contact_map_semantics": "Geometric contacts only.",
    }
    requested = []

    def fetch(job_id, path, token, state_key, label):
        requested.append(path)
        if path == "artifacts/md_analysis/summary.json":
            return json.dumps(summary)
        if path.endswith("timeseries.csv"):
            return "time_ns,backbone_rmsd_angstrom,ligand_rmsd_angstrom\n0,0.1,0.2\n"
        if path.endswith("residue_contacts.csv"):
            return "residue,contact_occupancy_fraction\nA:1,1.0\n"
        return None

    st = FakeStreamlit()
    render_remote_md_results(st, "native-md-analysis", {"tool": "md_analysis"}, "token", fetch)

    assert requested[0] == "artifacts/md_analysis/summary.json"
    assert "artifacts/md/study.json" not in requested
    assert any(call[0] == "line_chart" for call in st.calls)


def test_imported_analysis_uses_study_manifest_instead_of_native_summary_path():
    requested = []

    def fetch(job_id, path, token, state_key, label):
        requested.append(path)
        if path == "artifacts/md/study.json":
            return json.dumps({"schema_version": 1, "studies": []})
        return None

    render_remote_md_results(FakeStreamlit(), "legacy-analysis", {
        "tool": "legacy_md_analysis", "job_type": "legacy_md_analysis",
    }, "token", fetch)

    assert requested == ["artifacts/md/study.json"]
