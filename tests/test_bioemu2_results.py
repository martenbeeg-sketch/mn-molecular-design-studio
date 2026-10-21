from pathlib import Path

from mn_molecular_design_studio.bioemu2_results import trajectory_preview_artifacts
from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import StructureVisualization


def test_bioemu_xtc_preview_requires_sibling_topology_and_keeps_pair(tmp_path):
    run_dir = tmp_path / "run"
    output = run_dir / "artifacts" / "bioemu2" / "complex-1"
    output.mkdir(parents=True)
    xtc = output / "samples.xtc"
    topology = output / "topology.cif"
    xtc.write_bytes(b"xtc")
    topology.write_text("data_topology\n")
    item = {"relative_path": "artifacts/bioemu2/complex-1/samples.xtc"}

    assert trajectory_preview_artifacts("bioemu2", run_dir, [(item, xtc)]) == [(item, xtc, topology)]
    assert trajectory_preview_artifacts("vina", run_dir, [(item, xtc)]) == []


def test_bioemu_xtc_preview_skips_missing_or_external_topology(tmp_path):
    run_dir = tmp_path / "run"
    output = run_dir / "artifacts" / "bioemu2"
    output.mkdir(parents=True)
    xtc = output / "samples.xtc"
    xtc.write_bytes(b"xtc")
    item = {"relative_path": "samples.xtc"}

    assert trajectory_preview_artifacts("bioemu2_sampling", run_dir, [(item, xtc)]) == []

    outside = tmp_path / "topology.cif"
    outside.write_text("data_topology\n")
    (output / "topology.cif").symlink_to(outside)
    assert trajectory_preview_artifacts("bioemu2_sampling", run_dir, [(item, xtc)]) == []


def test_structure_visualization_serializes_xtc_coordinates_for_molstar():
    visualization = StructureVisualization(
        pdb="data_topology\n",
        trajectory={"coordinates_base64": "eHRj", "format": "xtc"},
    )

    assert visualization.to_dict()["trajectory"] == {"coordinates_base64": "eHRj", "format": "xtc"}
