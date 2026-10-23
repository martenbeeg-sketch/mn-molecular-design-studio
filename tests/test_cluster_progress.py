import json

import pytest

from mn_molecular_design_studio.cluster.progress import bindcraft2_progress


@pytest.mark.parametrize("status", ["running", "completed", "failed", "cancelled"])
def test_bindcraft_terminal_runs_do_not_report_interrupted_attempts_as_running(tmp_path, status):
    output = tmp_path / "artifacts/raw/bindcraft2/output"
    output.mkdir(parents=True)
    (output / ".campaign_state.json").write_text(json.dumps({
        "trajectories": 1638, "accepted": 8,
        "rejections": {"terminated": {"completed": 155, "harden": 1481}},
    }))
    (tmp_path / "metadata.json").write_text(json.dumps({"status": status, "partial_results": True}))
    (output.parent / "settings.json").write_text(json.dumps({"number_of_final_designs": 100}))
    progress = bindcraft2_progress(tmp_path)
    assert progress["trajectories_started"] == 1638
    assert progress["trajectories_completed"] == 1636
    assert progress["designs_accepted"] == 8
    assert progress["design_target"] == 100
    assert progress["trajectories_in_progress"] == (2 if status == "running" else 0)
    assert progress["trajectories_interrupted"] == (0 if status == "running" else 2)
