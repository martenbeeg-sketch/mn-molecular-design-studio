import copy
import io

import pytest
from streamlit.testing.v1 import AppTest
from test_database_browser import database

from mn_molecular_design_studio.md_plot_matrix import (
    available_plots, contact_rows, interaction_figure, plot_figure, replica_ids, validate_report_sources,
)


def _report():
    return {"success": True, "analysis_window": {"start_ns": 0, "end_ns": 99},
        "replicas": [{"replica": replica, "run_id": f"study-one-production-{replica}", "status": "completed",
            "endpoint_status": "completed", "delta_g_bind_kcal_mol": -replica, "pb_delta_g_bind_kcal_mol": replica}
            for replica in (1, 2, 3)],
        "replica_series": [{"replica": replica, "time_ns": [0, 1, 2],
            "backbone_rmsd_angstrom": [0, replica, replica / 2], "ligand_rmsd_angstrom": [0, replica / 2, replica],
            "protein_rg_angstrom": [10, 11, 12], "protein_rmsf_angstrom": [1, 2], "ligand_rmsf_angstrom": [2, 3]}
            for replica in (1, 2, 3)],
        "contact_consensus": [
            {"residue": "PHE215 · chain A", "mean_binding_importance_score": .9897,
                "mean_contact_occupancy": 1, "mean_hydrophobic_occupancy": .9897, "top_ligand_atom": "C15"},
            {"residue": "SER277 · chain A", "mean_binding_importance_score": 1.6423,
                "mean_contact_occupancy": .9991, "mean_hydrogen_bond_occupancy": .5024,
                "mean_water_bridge_occupancy": .09, "top_ligand_atom": "O2"}],
        "contact_matrix": {"residues": ["PHE215 · chain A", "SER277 · chain A"], "replicas": [1, 2, 3],
            "contact_occupancy": [[1, 1, 1], [.9, .9991, 1]],
            "hydrophobic_occupancy": [[.5, .9897, .7], [0, 0, 0]],
            "hydrogen_bond_occupancy": [[0, 0, 0], [.1, .5024, .3]],
            "water_bridge_occupancy": [[0, 0, 0], [.05, .09, .2]]},
        "ligand_depiction": {"atoms": [{"name": "C15", "symbol": "C", "x": -.2, "y": 0},
            {"name": "O2", "symbol": "O", "x": .2, "y": 0}], "bonds": [{"begin": 0, "end": 1, "order": 1}]}}


def test_three_replica_maps_draw_their_own_saved_percentages_on_the_edges():
    report = _report()
    original = copy.deepcopy(report)
    for replica, expected in ((1, "10.0%"), (2, "50.2%"), (3, "30.0%")):
        figure = interaction_figure(report, replica)
        labels = [text.get_text() for text in figure.axes[0].texts]
        assert expected in labels
        if replica == 2:
            assert "9.0%" in labels and "99.0%" in labels
            edges = {patch.get_gid() for patch in figure.axes[0].patches}
            assert "SER277 · chain A|O2|water_bridge|0.09" in edges
            assert "SER277 · chain A|O2|hydrogen_bond|0.5024" in edges
        buffer = io.BytesIO()
        figure.savefig(buffer, format="svg")
        assert b"<svg" in buffer.getvalue()
    assert report == original


def test_consensus_is_saved_mean_and_does_not_silently_substitute_another_replica():
    report = _report()
    report["contact_consensus"][1]["mean_water_bridge_occupancy"] = .0123
    labels = [text.get_text() for text in interaction_figure(report).axes[0].texts]
    assert "1.2%" in labels
    assert contact_rows(report, 2)[0]["occupancies"]["water_bridge"] == .09


def test_shared_hotspot_selection_matches_source_ranking_not_largest_generic_occupancy():
    assert contact_rows(_report(), 2, limit=1)[0]["residue"] == "SER277 · chain A"


def test_residue_positions_stay_the_same_across_replicas_while_percentages_change():
    import numpy as np

    figures = [interaction_figure(_report(), replica) for replica in (1, 2, 3)]
    positions = [[collection.get_offsets().copy() for collection in figure.axes[0].collections] for figure in figures]
    for candidate in positions[1:]:
        np.testing.assert_array_equal(candidate, positions[0])
    assert [patch.get_gid() for patch in figures[0].axes[0].patches] != [patch.get_gid() for patch in figures[1].axes[0].patches]


def test_percentage_labels_do_not_overlap_each_other_or_residue_labels():
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    figure = interaction_figure(_report(), 2)
    canvas = FigureCanvasAgg(figure)
    canvas.draw()
    renderer = canvas.get_renderer()
    labels = [text for text in figure.axes[0].texts if str(text.get_gid() or "").startswith("occupancy|")]
    residues = [text for text in figure.axes[0].texts if "chain A" in text.get_text()]
    for index, label in enumerate(labels):
        box = label.get_window_extent(renderer).padded(1)
        assert not any(box.overlaps(other.get_window_extent(renderer)) for other in labels[index + 1:])
        assert not any(box.overlaps(other.get_window_extent(renderer)) for other in residues)


@pytest.mark.parametrize("invalid", [True, None, float("nan"), float("inf"), -.1, 1.1, 50])
def test_invalid_or_missing_matrix_values_are_not_replaced_by_consensus_or_zero(invalid):
    report = _report()
    report["contact_matrix"]["water_bridge_occupancy"][1][1] = invalid
    assert "water_bridge" not in contact_rows(report, 2)[0]["occupancies"]


def test_matrix_axis_reordering_keeps_exact_replica_identity_and_rejects_ambiguity():
    report = _report()
    matrix = report["contact_matrix"]
    for field, rows in matrix.items():
        if field not in {"replicas", "residues"}:
            matrix[field] = [list(reversed(row)) for row in rows]
    matrix["replicas"] = [3, 2, 1]
    assert contact_rows(report, 1)[0]["occupancies"]["hydrogen_bond"] == .1
    matrix["replicas"] = [2, 2, 1]
    assert not contact_rows(report, 2)[0]["occupancies"]


def test_unknown_atom_does_not_get_original_apps_arbitrary_nearest_endpoint():
    report = _report()
    report["contact_consensus"][1]["top_ligand_atom"] = "not-recorded"
    figure = interaction_figure(report, 2)
    assert not any(patch.get_gid().startswith("SER277") for patch in figure.axes[0].patches)
    assert any("atom unavailable" in text.get_text() for text in figure.axes[0].texts)


@pytest.mark.parametrize("alter", ["names", "coordinates", "bond"])
def test_ambiguous_saved_depictions_fail_closed(alter):
    report = _report()
    if alter == "names":
        report["ligand_depiction"]["atoms"][1]["name"] = "C15"
    elif alter == "coordinates":
        report["ligand_depiction"]["atoms"][1]["x"] = float("nan")
    else:
        report["ligand_depiction"]["bonds"][0]["end"] = 99
    figure = interaction_figure(report, 2)
    assert not figure.axes[0].patches
    assert any("No atom endpoints are guessed" in text.get_text() for text in figure.axes[0].texts)


def test_plot_curves_keep_each_recorded_axis_and_honor_replica_selection():
    report = _report()
    report["replica_series"][1]["time_ns"] = [10, 11, 12]
    figure = plot_figure(report, "backbone_rmsd", selected_replicas=[2])
    assert len(figure.axes[0].lines) == 1
    assert list(figure.axes[0].lines[0].get_xdata()) == [10, 11, 12]
    assert list(figure.axes[0].lines[0].get_ydata()) == [0, 2, 1]
    assert figure.axes[0].lines[0].get_label() == "R2"


def test_failed_endpoint_numbers_and_failed_production_replicas_are_not_success_plots():
    report = _report()
    report["replicas"][1]["endpoint_status"] = "failed"
    assert len(plot_figure(report, "endpoint_energy").axes[0].patches) == 2
    report["replicas"][2]["status"] = "failed"
    assert replica_ids(report) == [1, 2]
    assert len(plot_figure(report, "backbone_rmsd").axes[0].lines) == 2


def test_inapplicable_salt_bridge_analysis_is_not_offered_as_a_measured_plot():
    report = _report()
    report["contact_consensus"][0]["mean_salt_bridge_occupancy"] = 0
    report["salt_bridge_applicable"] = False
    assert "salt_bridge_occupancy" not in available_plots(report)
    report["salt_bridge_applicable"] = True
    assert "salt_bridge_occupancy" in available_plots(report)


def test_results_overview_rejects_a_report_from_another_manifest_production(tmp_path):
    from mn_molecular_design_studio.md_results_ui import _render_aggregate_analysis
    from test_md_results_ui import FakeStreamlit

    st = FakeStreamlit()
    study = {"study_id": "expected-study", "production_runs": [{"run_id": "unrelated", "repeat_index": 1}],
        "analysis_runs": [{"run_id": "analysis", "status": "completed", "result_path": "artifacts/analysis.json"}]}
    _render_aggregate_analysis(st, tmp_path, study, json_loader=lambda relative, description: _report())
    assert any(name == "warning" and "cannot be linked" in args[0] for name, args, _ in st.calls)
    assert not any(name == "pyplot" for name, _, _ in st.calls)


def test_report_sources_require_exact_accessible_study_productions_and_replica_numbers():
    report = _report()
    productions = {f"study-one-production-{replica}": replica for replica in (1, 2, 3)}
    validate_report_sources(report, productions)
    with pytest.raises(ValueError, match="accessible production"):
        validate_report_sources(report, {"unrelated": 1})
    productions["study-one-production-2"] = 1
    with pytest.raises(ValueError, match="accessible production"):
        validate_report_sources(report, productions)
    report["replicas"][1]["replica"] = 1
    with pytest.raises(ValueError, match="ambiguous replica"):
        replica_ids(report)


def test_every_available_plot_renders_without_recalculating_the_saved_report():
    report = _report()
    original = copy.deepcopy(report)
    for plot in available_plots(report):
        for replica in (None, 2):
            figure = plot_figure(report, plot, replica)
            buffer = io.BytesIO()
            figure.savefig(buffer, format="png")
            assert buffer.getvalue().startswith(b"\x89PNG")
    assert report == original


def test_streamlit_overview_shows_all_three_replicas_together_and_enlarges_maps():
    source = f'''import streamlit as st
from mn_molecular_design_studio.md_plot_matrix import render_plot_matrix
render_plot_matrix(st, {_report()!r}, key_prefix="test")
'''
    app = AppTest.from_string(source).run(timeout=30)
    assert not app.exception
    assert app.multiselect[0].value == [1, 2, 3]
    assert sum(row.value == "Replica 1" for row in app.caption) == 4
    assert sum(row.value == "Replica 2" for row in app.caption) == 4
    assert sum(row.value == "Replica 3" for row in app.caption) == 4
    assert len(app.get("image")) == 12
    next(button for button in app.button if button.key == "test-interaction_network-2-enlarge").click().run(timeout=30)
    assert not app.exception
    assert len(app.get("image")) == 13
    app = AppTest.from_string(source).run(timeout=30)
    next(button for button in app.button if button.key == "test-expand-overview").click().run(timeout=30)
    assert not app.exception
    assert len(app.get("image")) == 24
    assert sum(row.value == "Replica 2" for row in app.caption) == 8


def test_database_overview_rejects_wrong_production_instead_of_merging_studies(database):
    from test_database_browser import USER
    from test_md_browser import _add_md, _cached, _inventory
    from mn_molecular_design_studio.md_browser import MDArtifactReader, study_observations

    for replica in (1, 2, 3):
        _add_md(database, replica=replica)
    identifier = _add_md(database, namespace="md.analysis")
    with database.connection() as connection:
        connection.execute("UPDATE molecule_observations SET context=json_set(context,'$.source_summary_relative_path','artifacts/aggregate/result.json') WHERE id=?", (identifier,))
    index, outcomes = _inventory(database)
    outcome = next(row for row in outcomes if row["outcome_kind"] == "Dynamics")
    observation = next(row for row in study_observations(database, USER, outcome) if row["id"] == identifier)
    reader = MDArtifactReader(database, USER, observation, index["runs"][outcome["primary"]["run_key"]])
    report = _report()
    report["replicas"][1]["run_id"] = "another-study-production"
    _cached(reader, "artifacts/aggregate/result.json", report)
    source = f'''from mn_molecular_design_studio.data_browser import ReadOnlyLibrary, saved_results, listed_results
from mn_molecular_design_studio.database_browser import build_browser_index, subject_results
from mn_molecular_design_studio.md_browser_ui import render_md_study
database=ReadOnlyLibrary({str(database.root)!r})
user={USER!r}
index=build_browser_index(database,user,collection="library",include_related=True)
raw=subject_results(database,index,{{"id":"all","molecule_ids":set(index["molecules"]),"group_ids":set(index["groups"])}})
outcome=next(row for row in saved_results(index,raw) if row["outcome_kind"]=="Dynamics")
render_md_study(database,user,index,listed_results([outcome])[0],{{}})
'''
    app = AppTest.from_string(source).run(timeout=30)
    assert not app.exception
    assert any("accessible production" in row.value for row in app.warning)
    assert not app.get("image")
    report["replicas"][1]["run_id"] = "study-one-production-2"
    _cached(reader, "artifacts/aggregate/result.json", report)
    app.run(timeout=30)
    assert not app.exception
    assert len(app.get("image")) == 12
