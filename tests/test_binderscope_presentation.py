import copy
import hashlib

import pytest

from mn_molecular_design_studio.binderscope_presentation import (
    contact_map, contact_network, contact_rows, highlight_selections, residue_identity,
    residue_label, residue_profile, residue_rows, result_coordinates, selected_contact,
)


@pytest.fixture
def assessment():
    return {
        "candidate_id": "design-one",
        "source": {"input_structure_digest": "verified-digest", "source": {"structure_id": "sql-structure-one", "record_name": "Saved design one"},
            "assessment": {"component_format": "linear protein", "b_factor_interpretation": "uninterpreted"},
            "residue_map": {
                "A": [{"canonical_number": 1, "chain": "X", "number": 15, "insertion_code": "A"},
                      {"canonical_number": 2, "chain": "X", "number": 15, "insertion_code": ""}],
                "B": [{"canonical_number": 1, "chain": "Y", "number": -2, "insertion_code": ""},
                      {"canonical_number": 2, "chain": "Y", "number": -1, "insertion_code": "B"}]}},
        "analysis": {"dsasa": 500., "n_hbonds": 1, "n_salt": 1, "n_unsat": 0,
            "contacts": [{"t": 1, "b": 2, "d": 3.2, "n": 8}, {"t": 2, "b": 1, "d": 4.3, "n": 1}],
            "tc": [1, 2], "bc": [1, 2], "hbonds": [{"t": 1, "b": 2, "d": 3.2, "k": "D>A"}],
            "salt": [{"t": 1, "b": 2, "d": 3.2}], "buns": [], "ss": "HE",
            "dsasa_res": [10., 30.], "rsa_res": [0.1, 0.3], "tracks": {"n_target_res": [1]},
            "liab": [{"kind": "Example motif", "pos": 2, "reason": "Test finding"}]},
        "hotspot_compliance": {"status": "assessed", "requested_count": 1, "contacted_count": 1, "missed_count": 0,
            "unmapped_count": 0, "cutoff_angstrom": 4.5, "residues": [
                {"canonical_number": 1, "requested": "T90", "saved_residue": "X:15A", "status": "contacted"}]},
        "target_reference_comparison": {"status": "not_applicable", "reason": "Sequence-only target: no fixed fold assumed."},
        "termini": {"status": "unavailable", "reason": "No confirmed confidence encoding."},
    }


def test_contact_labels_preserve_author_chains_gaps_and_insertion_codes(assessment):
    rows = contact_rows(assessment)
    assert rows[0]["Target residue"] == "X:15A"
    assert rows[0]["Binder residue"] == "Y:-1B"
    assert rows[1]["Target residue"] == "X:15"
    assert rows[0]["Interaction"] == "H-bond + salt bridge"
    assert rows[1]["Interaction"] == "Contact"
    assert rows[0]["Intended hotspot"] and not rows[1]["Intended hotspot"]
    assert residue_label(assessment, "A", 100) == "Unmapped A:100"
    assessment["source"]["residue_map"]["A"].append(assessment["source"]["residue_map"]["A"][0])
    assert residue_identity(assessment, "A", 1) is None


def test_all_plot_points_retain_exact_pair_identities_and_zero_is_not_missing(assessment):
    for figure in (contact_map(assessment), contact_network(assessment)):
        points = [list(row[:2]) for trace in figure.data if trace.customdata is not None for row in trace.customdata]
        assert sorted(points) == [[1, 2], [2, 1]]
        assert figure.to_json()
    assert contact_network(assessment).to_json() == contact_network(assessment).to_json()
    rows = residue_rows(assessment)
    assert rows[0]["pLDDT"] is None and rows[1]["Contacted target residues"] is None
    assert rows[1]["Secondary structure"] == "E" and rows[1]["Liabilities"] == "Example motif"
    assert len(residue_profile(assessment).data) == 1
    assessment["analysis"]["plddt_res"] = [0., 99.]
    assert residue_rows(assessment)[0]["pLDDT"] == 0.
    assert len(residue_profile(assessment).data) == 2


def test_plot_selection_only_accepts_real_residue_pairs(assessment):
    rows = contact_rows(assessment)
    assert selected_contact({"selection": {"points": [{"customdata": [1, 2, "anything"]}]}}, rows) == rows[0]
    assert selected_contact({"selection": {"points": [{"customdata": [99, 99]}]}}, rows) is None
    assert selected_contact({}, rows) is None
    assert selected_contact({"selection": {"points": [{"point_index": 0}]}}, rows) is None


def test_3d_highlighting_never_selects_other_insertions_or_guesses_unmapped_residues(assessment):
    rows = contact_rows(assessment)
    assert highlight_selections(assessment, "Selected contact", rows[0]) == ["X:15^A", "Y:-1^B"]
    assert highlight_selections(assessment, "Intended hotspots") == ["X:15^A"]
    assert highlight_selections(assessment, "Interface") == ["X:15^A", "X:15^_", "Y:-2^_", "Y:-1^B"]
    assessment["analysis"]["buns"] = [{"side": "t", "resi": 1}, {"side": "b", "resi": 10}]
    assert highlight_selections(assessment, "Buried polar atoms") == ["X:15^A"]
    assert highlight_selections(assessment, "Whole complex") == []


def test_coordinates_are_authorized_and_checksum_bound_not_loaded_from_report_paths(assessment, tmp_path):
    path = tmp_path / "complex.cif"
    path.write_text("data_saved\n")
    assessment["source"]["input_structure_digest"] = hashlib.sha256(path.read_bytes()).hexdigest()
    assessment["source"]["path"] = "/untrusted/location.cif"
    user = {"id": "authorized-user"}
    class Database:
        def structure_file(self, actor, identifier):
            assert actor == user and identifier == "sql-structure-one"
            return path
    assert result_coordinates(Database(), user, assessment) == path.read_text()
    path.write_text("data_replaced\n")
    with pytest.raises(ValueError, match="exact model"):
        result_coordinates(Database(), user, assessment)
    class Denied:
        def structure_file(self, actor, identifier):
            raise PermissionError("Not accessible")
    with pytest.raises(PermissionError):
        result_coordinates(Denied(), user, assessment)


def test_replayed_attachment_uses_resolved_sql_identity_not_stale_embedded_identifier(assessment, tmp_path):
    path = tmp_path / "peer-complex.cif"
    path.write_text("data_same_saved_coordinates\n")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assessment["source"]["input_structure_digest"] = digest
    assessment["_studio_sql_binding"] = {"structure_id": "peer-resolved-id", "structure_digest": digest}
    class Database:
        def structure_file(self, user, identifier):
            assert identifier == "peer-resolved-id"
            return path
    assert result_coordinates(Database(), {"id": "allowed"}, assessment) == path.read_text()
    assessment["_studio_sql_binding"]["structure_id"] = None
    with pytest.raises(ValueError, match="SQL structure identity"):
        result_coordinates(Database(), {"id": "allowed"}, assessment)
    assessment["_studio_sql_binding"] = {"structure_id": "peer-resolved-id", "structure_digest": "different"}
    with pytest.raises(ValueError, match="disagree"):
        result_coordinates(Database(), {"id": "allowed"}, assessment)


def test_inspector_renders_without_source_access_and_retains_unknowns(assessment):
    from streamlit.testing.v1 import AppTest

    original = copy.deepcopy(assessment)
    app = AppTest.from_string(f'''
from mn_molecular_design_studio.binderscope_ui import render_design_result
render_design_result({assessment!r})
''').run()
    assert not app.exception
    assert [tab.label for tab in app.tabs] == ["Interface · 2D & 3D", "Residue profile", "Chemistry & liabilities", "Hotspots, target & tags", "Parameters & evidence"]
    assert any("No confirmed confidence" in message.value for message in app.info)
    assert any("0" == metric.value for metric in app.metric)
    assert assessment == original
    app.radio[0].set_value("Residue network").run()
    assert not app.exception
    app.multiselect[0].set_value([]).run()
    assert not app.exception
    assert any("No interactions selected" in message.value for message in app.info)


def test_saved_structure_render_keeps_cartoon_and_confirmed_binder_only_confidence(assessment, tmp_path):
    from streamlit.testing.v1 import AppTest

    path = tmp_path / "complex.cif"
    path.write_text("data_example\n")
    assessment["source"]["input_structure_digest"] = hashlib.sha256(path.read_bytes()).hexdigest()
    app = AppTest.from_string(f'''
from pathlib import Path
import streamlit as st
from mn_molecular_design_studio.binderscope_inspector import render_inspector
from mn_molecular_design_studio.engine.workflow_engine.app.components import molstar_viewer
def capture(structures, **kwargs):
    st.session_state['captured'] = structures[0].to_dict()
molstar_viewer.molstar_custom_component = capture
class Database:
    def structure_file(self, user, identifier):
        assert user['id'] == 'allowed'
        return Path({str(path)!r})
render_inspector({assessment!r}, database=Database(), user={{'id': 'allowed'}})
''').run()
    assert not app.exception
    color = next(widget for widget in app.selectbox if widget.label == "Binder coloring")
    assert "pLDDT" not in color.options
    pair = next(widget for widget in app.selectbox if widget.label == "Inspect residue pair")
    pair.set_value((1, 2)).run()
    assert not app.exception
    captured = app.session_state['captured']
    assert captured["highlighted_selections"] == ["X:15^A", "Y:-1^B"]
    assert all(chain["representation_type"] == "cartoon" for chain in captured["chains"])
    assert [chain["chain_id"] for chain in captured["chains"]] == ["X", "Y"]
