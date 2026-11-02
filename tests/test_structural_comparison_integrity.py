from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from mn_molecular_design_studio.engine.workflow_engine.workflows.benchmark import (
    _target_aligned_binder_pose_metrics,
)
from mn_molecular_design_studio.complex_analysis_comparison import (
    EXTERNAL_ALIGNMENT_SCHEMA,
    external_alignment_template,
    normalize_external_alignment_record,
    validate_supplied_alignment,
)
from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import (
    StructureVisualization,
    validate_rigid_transform_matrix,
)


def _atom_line(serial: int, atom: str, residue: str, chain: str, number: int, xyz: tuple[float, float, float]) -> str:
    x, y, z = xyz
    return (
        f"ATOM  {serial:5d} {atom:>4s} {residue:>3s} {chain}{number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{20.0:6.2f}          C"
    )


def _write_complex(
    path: Path, *, target_numbers=(10, 11, 12), binder_numbers=(1, 2, 3),
    target_names=("ALA", "GLY", "SER"), binder_names=("ALA", "GLY", "SER"),
) -> None:
    coords = ((0.0, 0.0, 0.0), (4.0, 0.0, 0.0), (0.0, 5.0, 0.0))
    lines = []
    serial = 1
    for chain, numbers, names, offset in (("A", target_numbers, target_names, 0.0), ("B", binder_numbers, binder_names, 8.0)):
        for name, number, (x, y, z) in zip(names, numbers, coords):
            lines.append(_atom_line(serial, "CA", name, chain, number, (x + offset, y, z)))
            serial += 1
    path.write_text("\n".join([*lines, "END", ""]))


def _measure(reference: Path, prediction: Path) -> dict:
    return _target_aligned_binder_pose_metrics(
        reference,
        prediction,
        reference_binder_chains=["B"],
        reference_target_chains=["A"],
        prediction_binder_chains=["B"],
        prediction_target_chains=["A"],
    )


def test_target_aligned_metrics_require_identical_ordered_target_residue_mapping(tmp_path):
    reference = tmp_path / "reference.pdb"
    prediction = tmp_path / "prediction.pdb"
    _write_complex(reference)
    _write_complex(prediction)

    result = _measure(reference, prediction)

    assert result["pose_rmsd_status"] == "ok"
    assert result["target_alignment_mapping_status"].startswith("verified")
    assert result["binder_mapping_status"].startswith("verified")
    assert result["target_aligned_binder_rmsd"] == pytest.approx(0.0, abs=1e-5)


def test_target_aligned_metrics_do_not_truncate_mismatched_target_lists(tmp_path):
    reference = tmp_path / "reference.pdb"
    prediction = tmp_path / "prediction.pdb"
    _write_complex(reference)
    _write_complex(prediction, target_numbers=(10, 11))

    result = _measure(reference, prediction)

    assert result["target_aligned_binder_rmsd"] is None
    assert result["target_alignment_rmsd"] is None
    assert result["target_alignment_ca_count"] == 0
    assert result["target_alignment_mapping_status"] == "unavailable"
    assert "ordered target residue identifiers" in result["pose_rmsd_status"]


def test_target_aligned_metrics_reject_mismatched_residue_numbers(tmp_path):
    reference = tmp_path / "reference.pdb"
    prediction = tmp_path / "prediction.pdb"
    _write_complex(reference)
    _write_complex(prediction, target_numbers=(10, 11, 99))

    result = _measure(reference, prediction)

    assert result["target_aligned_binder_rmsd"] is None
    assert result["target_alignment_mapping_status"] == "unavailable"
    assert "identifiers, insertion codes, or residue names differ" in result["pose_rmsd_status"]


def test_target_aligned_metrics_reject_mismatched_residue_names(tmp_path):
    reference = tmp_path / "reference.pdb"
    prediction = tmp_path / "prediction.pdb"
    _write_complex(reference)
    _write_complex(prediction, target_names=("ALA", "SER", "SER"))

    result = _measure(reference, prediction)

    assert result["target_aligned_binder_rmsd"] is None
    assert result["target_alignment_mapping_status"] == "unavailable"


def test_target_aligned_metrics_refuse_unselected_multi_model_files(tmp_path):
    reference = tmp_path / "reference.pdb"
    prediction = tmp_path / "prediction-multimodel.pdb"
    _write_complex(reference)
    _write_complex(prediction)
    coordinates = [line for line in prediction.read_text().splitlines() if line.startswith("ATOM")]
    prediction.write_text("\n".join([
        "MODEL        1", *coordinates, "ENDMDL", "MODEL        2", *coordinates, "ENDMDL", "END", "",
    ]))

    result = _measure(reference, prediction)

    assert result["target_aligned_binder_rmsd"] is None
    assert result["target_alignment_mapping_status"] == "unavailable"
    assert "contains 2 models" in result["pose_rmsd_status"]


def test_binder_pose_comparison_is_unavailable_when_binder_lists_mismatch(tmp_path):
    reference = tmp_path / "reference.pdb"
    prediction = tmp_path / "prediction.pdb"
    _write_complex(reference)
    _write_complex(prediction, binder_numbers=(1, 2))

    result = _measure(reference, prediction)

    assert result["target_alignment_rmsd"] == pytest.approx(0.0, abs=1e-5)
    assert result["target_aligned_binder_rmsd"] is None
    assert result["binder_mapping_status"] == "unavailable"
    assert "ordered binder residue identifiers" in result["pose_rmsd_status"]


def _source(identity: str, path: Path) -> dict:
    with path.open("rb") as stream:
        structure_digest = hashlib.file_digest(stream, "sha256").hexdigest()
    return {
        "identity": identity,
        "model_id": "model-1",
        "structure_digest": structure_digest,
        "path": path,
        "target_chains": ["A"],
    }


def _alignment_record(reference: dict, moving: dict, *, pairs=None) -> dict:
    if pairs is None:
        pairs = [
            {"reference": {"chain": "A", "number": number, "insertion": "", "residue_name": name},
             "moving": {"chain": "A", "number": number, "insertion": "", "residue_name": name}}
            for number, name in ((10, "ALA"), (11, "GLY"), (12, "SER"))
        ]
    return {
        "reference_identity": reference["identity"],
        "reference_model_id": reference["model_id"],
        "reference_digest": reference["structure_digest"],
        "moving_identity": moving["identity"],
        "moving_model_id": moving["model_id"],
        "moving_digest": moving["structure_digest"],
        "matrix_4x4": [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
        "residue_correspondence": pairs,
        "source": "external transform record",
        "method": "externally supplied",
        "mapping_note": "Residue map reviewed by operator",
    }


def test_supplied_transform_and_explicit_mapping_are_bound_to_selected_structure_models(tmp_path):
    reference_path = tmp_path / "reference.pdb"
    moving_path = tmp_path / "moving.pdb"
    _write_complex(reference_path)
    _write_complex(moving_path)
    reference = _source("db:molecule-1:structure-1", reference_path)
    moving = _source("results:run-2:model-1", moving_path)

    report = validate_supplied_alignment(reference, moving, _alignment_record(reference, moving))

    assert report["valid"] is True
    assert report["correspondence_count"] == 3
    assert report["reference_unmapped_residues"] == 0
    assert report["moving_unmapped_residues"] == 0
    assert report["source"] == "external transform record"
    assert any("does not verify alignment fit accuracy" in warning for warning in report["warnings"])


def test_supplied_mapping_reports_partial_coverage_and_rejects_duplicates(tmp_path):
    reference_path = tmp_path / "reference.pdb"
    moving_path = tmp_path / "moving.pdb"
    _write_complex(reference_path)
    _write_complex(moving_path)
    reference = _source("reference", reference_path)
    moving = _source("moving", moving_path)
    one_pair = _alignment_record(reference, moving)["residue_correspondence"][:1]

    partial = validate_supplied_alignment(reference, moving, _alignment_record(reference, moving, pairs=one_pair))
    assert partial["valid"] is True
    assert partial["reference_unmapped_residues"] == 2
    assert any("partial" in warning for warning in partial["warnings"])

    duplicate = validate_supplied_alignment(
        reference, moving, _alignment_record(reference, moving, pairs=[one_pair[0], one_pair[0]])
    )
    assert duplicate["valid"] is False
    assert any("ambiguous" in error for error in duplicate["errors"])


def test_supplied_correspondence_supports_explicit_different_author_numbering(tmp_path):
    reference_path = tmp_path / "reference.pdb"
    moving_path = tmp_path / "moving.pdb"
    _write_complex(reference_path, target_numbers=(10, 11, 12))
    _write_complex(moving_path, target_numbers=(100, 101, 102))
    reference = _source("reference", reference_path)
    moving = _source("moving", moving_path)
    pairs = [
        {"reference": {"chain": "A", "number": left, "insertion": "", "residue_name": name},
         "moving": {"chain": "A", "number": right, "insertion": "", "residue_name": name}}
        for left, right, name in zip((10, 11, 12), (100, 101, 102), ("ALA", "GLY", "SER"))
    ]

    report = validate_supplied_alignment(reference, moving, _alignment_record(reference, moving, pairs=pairs))

    assert report["valid"] is True
    assert report["correspondence_count"] == 3
    assert report["reference_unmapped_residues"] == 0
    assert report["moving_unmapped_residues"] == 0


def test_supplied_alignment_rejects_wrong_model_and_residue_identity(tmp_path):
    reference_path = tmp_path / "reference.pdb"
    moving_path = tmp_path / "moving.pdb"
    _write_complex(reference_path)
    _write_complex(moving_path)
    reference = _source("reference", reference_path)
    moving = _source("moving", moving_path)
    record = _alignment_record(reference, moving)
    record["moving_model_id"] = "model-2"
    record["moving_digest"] = "stale-structure-digest"
    record["residue_correspondence"][0]["moving"]["number"] = 999

    report = validate_supplied_alignment(reference, moving, record)

    assert report["valid"] is False
    assert any("model does not match" in error for error in report["errors"])
    assert any("digest does not match" in error for error in report["errors"])
    assert any("does not identify a residue" in error for error in report["errors"])


def test_external_alignment_record_template_is_bound_and_normalizes_provenance(tmp_path):
    reference_path = tmp_path / "reference.pdb"
    moving_path = tmp_path / "moving.pdb"
    _write_complex(reference_path)
    _write_complex(moving_path)
    reference = _source("db:molecule-1:structure-1", reference_path)
    moving = _source("results:run-2:model-1", moving_path)
    payload = external_alignment_template(reference, moving)
    assert payload["schema"] == EXTERNAL_ALIGNMENT_SCHEMA
    assert payload["reference"]["structure_digest"] == reference["structure_digest"]
    assert payload["moving"]["model_id"] == moving["model_id"]

    payload["transform"].update({"method": "Kabsch, external", "source": {
        "tool": "ExampleAlign", "version": "2.1", "record_id": "alignment-77", "uri": "urn:test:alignment-77",
    }})
    payload["residue_correspondence"] = _alignment_record(reference, moving)["residue_correspondence"]
    payload["mapping_note"] = "External residue map from alignment-77"
    record = normalize_external_alignment_record(payload, reference, moving)
    report = validate_supplied_alignment(reference, moving, record)

    assert report["valid"] is True
    assert report["source"] == "ExampleAlign 2.1"
    assert report["external_provenance"] == {
        "schema": EXTERNAL_ALIGNMENT_SCHEMA,
        "source": {"tool": "ExampleAlign", "version": "2.1", "record_id": "alignment-77", "uri": "urn:test:alignment-77"},
    }


@pytest.mark.parametrize(("change", "message"), [
    (lambda payload: payload.update(schema="mnds.external-alignment/v2"), "Unsupported or missing schema"),
    (lambda payload: payload["transform"].update(direction="reference_to_moving"), "direction must be explicitly"),
    (lambda payload: payload["moving"].update(model_id="model-2"), "moving model_id does not match"),
    (lambda payload: payload["reference"].update(structure_digest="stale"), "reference structure_digest does not match"),
])
def test_external_alignment_record_rejects_stale_or_unsupported_envelopes(tmp_path, change, message):
    reference_path = tmp_path / "reference.pdb"
    moving_path = tmp_path / "moving.pdb"
    _write_complex(reference_path)
    _write_complex(moving_path)
    reference = _source("reference", reference_path)
    moving = _source("moving", moving_path)
    payload = external_alignment_template(reference, moving)
    payload["transform"].update({"method": "external method", "source": {"tool": "tool", "version": "1"}})
    change(payload)

    with pytest.raises(ValueError, match=message):
        normalize_external_alignment_record(payload, reference, moving)


def test_supplied_alignment_refuses_multi_model_coordinate_files(tmp_path):
    reference_path = tmp_path / "reference.pdb"
    moving_path = tmp_path / "moving-multimodel.pdb"
    _write_complex(reference_path)
    _write_complex(moving_path)
    coordinates = moving_path.read_text().splitlines()
    atom_lines = [line for line in coordinates if line.startswith("ATOM")]
    moving_path.write_text("\n".join([
        "MODEL        1", *atom_lines, "ENDMDL", "MODEL        2", *atom_lines, "ENDMDL", "END", "",
    ]))
    reference = _source("reference", reference_path)
    moving = _source("moving", moving_path)

    report = validate_supplied_alignment(reference, moving, _alignment_record(reference, moving))

    assert report["valid"] is False
    assert any("contains 2 models" in error for error in report["errors"])


def test_rigid_transform_validation_rejects_scale_reflection_and_nonfinite_values():
    with pytest.raises(ValueError, match="orthonormal"):
        validate_rigid_transform_matrix([[2, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    with pytest.raises(ValueError, match="determinant"):
        validate_rigid_transform_matrix([[-1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    with pytest.raises(ValueError, match="finite"):
        validate_rigid_transform_matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, float("nan"), 0], [0, 0, 0, 1]])


def test_molstar_visualization_serializes_validated_transform():
    visualization = StructureVisualization(
        pdb="END\n",
        transform_matrix=[[1, 0, 0, 3], [0, 1, 0, 4], [0, 0, 1, 5], [0, 0, 0, 1]],
    )
    assert visualization.to_dict()["transform_matrix"][0][3] == 3.0
    with pytest.raises(ValueError, match="4x4"):
        StructureVisualization(pdb="END\n", transform_matrix=[[1, 0], [0, 1]])


def test_comparison_panel_applies_only_the_validated_supplied_transform(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest

    from mn_molecular_design_studio import complex_analysis_comparison_ui as comparison_ui
    from mn_molecular_design_studio import complex_analysis_ui as workspace
    from mn_molecular_design_studio.engine.workflow_engine.app.components import molstar_viewer

    reference_path = tmp_path / "reference.pdb"
    moving_path = tmp_path / "moving.pdb"
    _write_complex(reference_path)
    _write_complex(moving_path)
    monkeypatch.setattr(workspace, "_viewer_text", lambda path, context=None: Path(path).read_text())
    shown = []
    monkeypatch.setattr(molstar_viewer, "molstar_custom_component", lambda visualizations, **kwargs: shown.extend(visualizations))
    form_key = comparison_ui._key("|".join((
        "reference", "model-1", _source("reference", reference_path)["structure_digest"], "A",
        "moving", "model-1", _source("moving", moving_path)["structure_digest"], "A",
    )))
    mapping = _alignment_record(_source("reference", reference_path), _source("moving", moving_path))["residue_correspondence"]
    reference_digest = _source("reference", reference_path)["structure_digest"]
    moving_digest = _source("moving", moving_path)["structure_digest"]
    script = f'''\
from pathlib import Path
from mn_molecular_design_studio.complex_analysis_comparison_ui import _render_supplied_alignment_overlay
sources = [
    {{"identity": "reference", "model_id": "model-1", "structure_digest": {reference_digest!r}, "label": "Reference", "path": Path({str(reference_path)!r}), "context": {{}}, "target_chains": ["A"]}},
    {{"identity": "moving", "model_id": "model-1", "structure_digest": {moving_digest!r}, "label": "Moving", "path": Path({str(moving_path)!r}), "context": {{}}, "target_chains": ["A"]}},
]
_render_supplied_alignment_overlay(sources)
'''
    app = AppTest.from_string(script, default_timeout=10).run()
    assert not app.exception
    app.text_area(key=f"complex-analysis-transform-matrix-{form_key}").set_value(
        "[[1,0,0,5],[0,1,0,0],[0,0,1,0],[0,0,0,1]]"
    )
    app.text_area(key=f"complex-analysis-transform-map-{form_key}").set_value(__import__("json").dumps(mapping))
    app.text_input(key=f"complex-analysis-transform-source-{form_key}").set_value("verified external record")
    app.text_input(key=f"complex-analysis-transform-method-{form_key}").set_value("operator supplied")
    app.button(key=f"complex-analysis-transform-apply-{form_key}").click().run()

    assert not app.exception
    assert len(shown) == 2
    assert shown[0].transform_matrix is None
    assert shown[1].transform_matrix[0][3] == 5.0
    assert any("does not verify alignment fit accuracy" in warning.value for warning in app.warning)

    from mn_molecular_design_studio.complex_analysis_comparison import external_alignment_template

    reference = _source("reference", reference_path)
    moving = _source("moving", moving_path)
    payload = external_alignment_template(reference, moving)
    payload["transform"].update({
        "matrix_4x4": [[1, 0, 0, 7], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
        "method": "external least-squares fit",
        "source": {"tool": "ExampleAlign", "version": "2.1", "record_id": "record-9", "uri": "urn:test:record-9"},
    })
    payload["residue_correspondence"] = mapping
    app.radio(key=f"complex-analysis-transform-input-mode-{form_key}").set_value("External JSON record").run()
    app.text_area(key=f"complex-analysis-transform-record-{form_key}").set_value(__import__("json").dumps(payload))
    app.button(key=f"complex-analysis-transform-apply-{form_key}").click().run()

    assert not app.exception
    assert len(shown) == 4
    assert shown[3].transform_matrix[0][3] == 7.0
    assert any("ExampleAlign 2.1" in caption.value for caption in app.caption)
