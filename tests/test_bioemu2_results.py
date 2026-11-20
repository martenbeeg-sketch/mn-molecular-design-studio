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


def test_bioemu_propagation_cif_has_numeric_entities_and_exact_coordinates(tmp_path):
    import gemmi
    from mn_molecular_design_studio.engine.workflow_engine.workflows.bioemu2 import _write_propagation_cif

    source = tmp_path / "complex.pdb"
    source.write_text(
        "ATOM      1  CA  ALA A  14       1.000   4.000   2.000  1.00 70.00           C\n"
        "ATOM      2  CA  GLY A  15       2.000   4.000   4.000  1.00 71.00           C\n"
        "TER\n"
        "ATOM      3  CA  SER B  23       1.000   8.000   2.000  1.00 60.00           C\n"
        "ATOM      4  CA  TYR B  24       2.000   8.000   4.000  1.00 61.00           C\n"
        "END\n"
    )
    destination = tmp_path / "complex.cif"
    _write_propagation_cif(source, destination)
    block = gemmi.cif.read(str(destination)).sole_block()
    assert set(block.find_values("_atom_site.label_entity_id")) == {"1", "2"}
    assert list(block.find_values("_entity_poly.pdbx_seq_one_letter_code_can")) == ["AG", "SY"]
    assert list(block.find_values("_atom_site.auth_seq_id")) == ["14", "15", "23", "24"]

    def atoms(path):
        return [(chain.name, str(residue.seqid), atom.name, atom.pos.tolist(), atom.b_iso)
                for chain in gemmi.read_structure(str(path))[0] for residue in chain for atom in residue]

    assert atoms(source) == atoms(destination)
