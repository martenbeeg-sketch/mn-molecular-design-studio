from mn_molecular_design_studio.complex_analysis_collection_ui import (
    _design_collection_label,
    _is_design_collection,
)


def test_design_collection_label_uses_campaign_target_and_source_host():
    label = _design_collection_label({
        "id": "collection-1",
        "name": "R3 (31 aa) · accepted RFD4",
        "origin": {
            "collection_kind": "curated",
            "campaign_name": "Tau R3 native campaign",
            "run_id": "20261030-070659-2044e3d6",
            "task_group": "design",
            "workstation": "PC004033",
            "tool": "rfd4_proteina",
            "design_target_name": "R3 (31 aa)",
        },
    })

    assert label["target"] == "R3 (31 aa)"
    assert label["run"] == (
        "Tau R3 native campaign · design/20261030-070659-2044e3d6 · PC004033"
    )
    assert label["generator"] == "rfd4_proteina"


def test_design_collection_label_does_not_call_analysis_tool_a_generator():
    label = _design_collection_label({
        "name": "Acceptance evidence",
        "origin": {
            "collection_kind": "generated",
            "tool": "validation_acceptance",
            "run_id": "20261031-085246-eb53803c",
            "task_group": "design",
            "design_target_molecule_id": "mol-target-123",
        },
    })

    assert label["generator"] == "Not recorded"
    assert label["target"] == "Target record mol-target-123"
    assert label["run"].startswith("Campaign not recorded · design/")


def test_job_collection_name_is_used_when_campaign_name_was_not_normalized():
    label = _design_collection_label({
        "name": "PD-L1 distributed complete BindCraft2 acceptance",
        "origin": {
            "type": "job",
            "collection_kind": "generated",
            "tool": "bindcraft2_campaign",
            "run_id": "20261029-094457-fb2fd80a",
            "task_group": "design",
            "workstation": "pc002975",
        },
    })

    assert label["run"] == (
        "PD-L1 distributed complete BindCraft2 acceptance · "
        "design/20261029-094457-fb2fd80a · pc002975"
    )
    assert label["generator"] == "bindcraft2_campaign"


def test_design_collection_picker_excludes_reference_libraries_only():
    assert not _is_design_collection({"collection_kind": "reference", "origin": {}})
    assert _is_design_collection({"collection_kind": "generated", "origin": {}})
    assert _is_design_collection({"origin": {"collection_kind": "curated"}})
