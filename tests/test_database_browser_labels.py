from mn_molecular_design_studio.database_browser_labels import analysis_inventory, analysis_status, calculation_label, readable_name, run_label, target_construct_label, target_label


def test_run_labels_expose_time_without_using_a_hash_as_the_name():
    assert run_label("20261031-090421-b051a2bb") == "2026-10-31 09:04:21"
    assert run_label("imported-research") == "imported-research"
    assert run_label("") == "Run not recorded"


def test_friendly_target_names_never_substitute_an_opaque_identifier():
    target_id = "mol-" + "a" * 32
    index = {"molecules": {target_id: {"name": target_id, "target_names": {"target", "R3 (31 aa)"},
        "import_names": set(), "components": [{"type": "protein", "sequence": "ACD"}]}}}
    assert target_label(index, target_id) == "R3 (31 aa)"
    index["molecules"][target_id]["target_names"] = set()
    assert target_label(index, target_id) == "Unnamed protein target (3 aa)"
    assert target_label(index, "mol-unknown") == "Target name not recorded"
    assert not readable_name(target_id)
    assert not readable_name("idr_sequence_carrier")
    assert readable_name("Tau214–283")


def test_imported_calculation_headers_show_recorded_replicates_not_import_wording():
    assert calculation_label({"tool": "vina", "run_id": "mnligand-source", "metrics": {"replicate": 0, "attempt_index": 2.0}}) == "Vina · Replicate 0 · Attempt 2"
    assert calculation_label({"tool": "boltz2", "run_id": "mnligand-source", "model": "0"}) == "Boltz-2 · Model / pose 0"
    assert calculation_label({"tool": "vina", "candidate": "legacy-replicate-0001"}) == "Vina · Replicate 1"


def test_construct_labels_use_the_recorded_sequence_length_and_keep_multiple_chains_explicit():
    index = {"molecules": {"target": {"name": "TRβ", "target_names": set(), "import_names": set(),
        "components": [{"type": "protein", "sequence": "EFGH"}]}}}
    assert target_construct_label(index, "target") == "TRβ · 4 aa"
    index["molecules"]["target"]["name"] = "TRβ (4 aa)"
    assert target_construct_label(index, "target") == "TRβ (4 aa)"
    index["molecules"]["target"]["name"] = "TRβ"
    index["molecules"]["target"]["components"].append({"type": "protein", "sequence": "ACD"})
    assert target_construct_label(index, "target") == "TRβ · 7 aa · 2 protein chains"
    assert target_construct_label(index, "unknown") == "Target name not recorded"


def test_analysis_inventory_keeps_failures_other_models_and_unknown_status_visible():
    rows = [{"tool": "biotite_interface", "status": status, "structure_relation": relation} for status, relation in (
        ("completed", "Selected structure"), ("completed", "Other structure of this design"),
        ("failed", "Selected structure"), ("", "Selected structure"))]
    rows.extend([{"tool": "ipsae", "status": "completed", "structure_relation": "Other structure of this design"},
                 {"tool": "validation_acceptance", "status": "completed"}])
    inventory = analysis_inventory(rows)
    assert inventory["Biotite"] == {"same_structure": 1, "other_structure": 1, "failed": 1, "unknown": 1}
    assert "Filter decision" not in inventory
    assert analysis_status(inventory, "Biotite") == "Recorded"
    assert analysis_status(inventory, "ipSAE") == "Other model"
    assert analysis_status(inventory, "PEP-Patch") == "—"
    assert analysis_status(analysis_inventory([{"tool": "pep_patch", "status": "failed"}]), "PEP-Patch") == "Failed / unavailable"
    assert analysis_status(analysis_inventory([{"tool": "ipsae", "status": "running"}]), "ipSAE") == "Status unknown"
