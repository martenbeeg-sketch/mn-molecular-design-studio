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
