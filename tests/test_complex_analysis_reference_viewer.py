from __future__ import annotations

from mn_molecular_design_studio import complex_analysis_ui


class _ReferenceDatabase:
    def __init__(self, saved, reference_path=None, error=None):
        self.saved = saved
        self.reference_path = reference_path
        self.error = error
        self.lookup = None
        self.opened = []

    def get_complex_reference_context(self, user, **kwargs):
        self.lookup = (user, kwargs)
        if self.error:
            raise self.error
        return self.saved

    def structure_file(self, user, structure_id):
        self.opened.append((user, structure_id))
        if self.error:
            raise self.error
        return self.reference_path


def test_paired_reference_view_uses_exact_saved_reference_and_selected_component(tmp_path, monkeypatch):
    reference = tmp_path / "reference.cif"
    reference.write_text("reference coordinates")
    database = _ReferenceDatabase({
        "configured": True, "status": "partial", "reference_structure_id": "reference-structure",
        "reference_molecule_id": "reference-molecule", "reference_component_chain": "R",
        "reference_record_name": "Full construct", "reference_structure_label": "Saved reference",
        "mapping": {"mapping_status": "partial"},
    }, reference)
    user = {"id": "user-1"}
    record = {"id": "source-molecule"}
    structure = {"id": "source-structure"}
    components = [{"type": "protein", "source_chain": "A"},
                  {"type": "protein", "source_chain": "B"}]
    state = {"complex-analysis-reference-component-view-key": "B"}
    monkeypatch.setattr(complex_analysis_ui.st, "session_state", state)
    monkeypatch.setattr(complex_analysis_ui, "component_inventory", lambda *_args: [{"source_chain": "R"}])

    result = complex_analysis_ui._paired_reference_for_viewer(
        database, user, record, structure, components, "view-key",
    )

    assert result["saved"]["reference_structure_id"] == "reference-structure"
    assert result["source_chain"] == "B"
    assert result["reference_id"] == "reference-structure"
    assert result["path"] == reference
    assert database.lookup == (user, {
        "source_molecule_id": "source-molecule", "source_structure_id": "source-structure",
        "source_chain": "B",
    })
    assert database.opened == [(user, "reference-structure")]


def test_reference_view_does_not_change_core_view_when_absent_or_stale(monkeypatch):
    state = {}
    monkeypatch.setattr(complex_analysis_ui.st, "session_state", state)
    components = [{"type": "protein", "source_chain": "A"}]
    args = (object(), {"id": "user"}, {"id": "source"}, {"id": "structure"}, components, "key")

    assert complex_analysis_ui._paired_reference_for_viewer(
        _ReferenceDatabase(None), *args[1:],
    ) is None

    stale = _ReferenceDatabase({"configured": True, "status": "stale"})
    assert complex_analysis_ui._paired_reference_for_viewer(stale, *args[1:]) is None
    assert stale.opened == []


def test_inaccessible_reference_is_not_loaded_or_used_as_an_overlay(monkeypatch):
    state = {}
    monkeypatch.setattr(complex_analysis_ui.st, "session_state", state)
    denied = _ReferenceDatabase(None, error=PermissionError("not authorized"))

    result = complex_analysis_ui._paired_reference_for_viewer(
        denied, {"id": "user"}, {"id": "source"}, {"id": "structure"},
        [{"type": "protein", "source_chain": "A"}], "key",
    )

    assert result is None
