from mn_molecular_design_studio.result_presentation import candidate_state_rows


def test_results_summary_separates_source_identity_acceptance_and_measurements():
    row = {
        "candidate_id": "cs_studio-identity",
        "binder_sequence": "ACDEFG",
        "complex_pdb": "candidate.cif",
        "raw_metadata": {
            "candidate_set": {
                "native_candidate_id": "rfd4_00007",
                "native_filter": "passed",
                "shared_filter": {"status": "unresolved"},
            },
            "human_review": {"selected": True},
        },
    }

    visible = candidate_state_rows([row])[0]
    assert visible["Candidate"] == "cs_studio-identity"
    assert visible["Source candidate"] == "rfd4_00007"
    assert visible["Measurement"] == "not_measured"
    assert visible["Measurement detail"] == "No independent measurement was recorded for this candidate."
    assert visible["Native generator acceptance"] == "passed"
    assert visible["Studio filtering"] == "unresolved"
    assert visible["Manual selection"] is True


def test_pipeline_results_use_original_candidate_id_and_keep_filter_separate_from_validation():
    row = {
        "candidate_id": "cs_rfd4_00007",
        "metrics": {"native_pass_filters": True, "iptm": 0.82},
        "raw_metadata": {
            "pipeline_original_candidate_id": "rfd4_00007",
            "validation_acceptance": {"status": "passed"},
            "selection": {"status": "pass", "routed_to": "out", "selected": True},
        },
    }

    visible = candidate_state_rows([row])[0]
    assert visible["Source candidate"] == "rfd4_00007"
    assert visible["Native generator acceptance"] == "passed"
    assert visible["Filter decision"] == "pass"
    assert visible["Filter route"] == "out"
    assert visible["Studio filtering"] == "not_evaluated"
    assert visible["Validation"] == "passed"
