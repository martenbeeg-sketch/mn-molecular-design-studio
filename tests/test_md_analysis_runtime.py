from __future__ import annotations

import json
import hashlib

import numpy as np
import pytest

from mn_molecular_design_studio.engine.workflow_engine.tools.md_analysis import (
    analyze,
    kabsch_align,
    residue_contact_metrics,
)
from mn_molecular_design_studio.engine.workflow_engine.workflows import md_analysis


def test_kabsch_returns_angstrom_rmsd_and_rotates_all_atoms():
    reference = np.asarray(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0]],
        dtype=float,
    )
    # 90 degree rotation + translation, plus a ligand atom 2 A away from its
    # reference after the protein fit.
    rotation = np.asarray([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    moved = reference @ rotation + np.asarray([8.0, -3.0, 4.0])
    moved[3] += np.asarray([0.0, 0.0, 2.0])

    aligned, rmsd = kabsch_align(moved, reference, np.asarray([0, 1, 2]))

    assert rmsd == pytest.approx(0.0, abs=1e-10)
    assert aligned[:3] == pytest.approx(reference[:3], abs=1e-10)
    assert np.linalg.norm(aligned[3] - reference[3]) == pytest.approx(2.0)


def test_residue_contact_metrics_uses_heavy_atom_minimum_and_cutoff_angstrom():
    protein = np.asarray([[0.0, 0.0, 0.0], [0.0, 5.0, 0.0], [9.0, 0.0, 0.0]])
    ligand = np.asarray([[3.0, 0.0, 0.0], [8.0, 0.0, 0.0]])
    residue_ids = np.asarray([10, 10, 11])

    distances = residue_contact_metrics(protein, ligand, residue_ids, 4.5)

    assert distances == {10: pytest.approx(3.0), 11: pytest.approx(1.0)}
    assert sum(value <= 4.5 for value in distances.values()) == 2


def test_streaming_analysis_writes_unit_labeled_outputs_and_exact_source_identity(tmp_path):
    md = pytest.importorskip("mdtraj")
    topology = md.Topology()
    chain = topology.add_chain("A")
    residue = topology.add_residue("ALA", chain, resSeq=10)
    atom_names = ["N", "CA", "C", "O"]
    atoms = [
        topology.add_atom(name, element, residue)
        for name, element in zip(atom_names, [md.element.nitrogen, md.element.carbon, md.element.carbon, md.element.oxygen])
    ]
    for left, right in zip(atoms, atoms[1:]):
        topology.add_bond(left, right)
    ligand_residue = topology.add_residue("LIG", chain, resSeq=501)
    ligand_atoms = [
        topology.add_atom("C1", md.element.carbon, ligand_residue),
        topology.add_atom("O1", md.element.oxygen, ligand_residue),
    ]
    topology.add_bond(*ligand_atoms)
    base = np.asarray(
        [[0.00, 0.00, 0.00], [0.10, 0.00, 0.00], [0.20, 0.00, 0.00],
         [0.20, 0.10, 0.00], [0.35, 0.00, 0.00], [0.40, 0.00, 0.00]],
        dtype=np.float32,
    )
    xyz = np.stack((base.copy(), base.copy(), base.copy()))
    xyz[1] += np.asarray([0.0, 0.0, 0.3], dtype=np.float32)
    xyz[2] += np.asarray([0.0, 0.0, 0.6], dtype=np.float32)
    xyz[2, 4:, 0] += 0.2  # ligand moves 2 A relative to the aligned protein
    trajectory = md.Trajectory(xyz=xyz, topology=topology, time=np.asarray([0.0, 1000.0, 2000.0]))
    topology_path = tmp_path / "system.pdb"
    trajectory_path = tmp_path / "production.xtc"
    trajectory[0].save_pdb(str(topology_path))
    trajectory.save_xtc(str(trajectory_path))
    output_dir = tmp_path / "analysis"
    source_segments = [{"run_id": "continuation", "segment_duration_ns": 99.0, "predecessor_run_id": "production"}]

    summary = analyze({
        "topology_path": str(topology_path),
        "trajectory_path": str(trajectory_path),
        "topology_source_path": "artifacts/legacy/prep/system.pdb",
        "trajectory_source_path": "artifacts/legacy/production/production.xtc",
        "output_dir": str(output_dir),
        "ligand_selection": "resname LIG",
        "protein_selection": "protein",
        "stride": 1,
        "start_ns": 0.0,
        "contact_cutoff_angstrom": 4.5,
        "pbc_mode": "auto",
        "chunk_size": 1,
        "source_segments": source_segments,
    })

    assert summary["units"] == {
        "time": "ns", "distance": "angstrom", "rmsd": "angstrom",
        "selected_frame_count": "frames",
        "contact_count_mean": "geometric heavy-atom residue contacts per analyzed frame",
    }
    assert summary["success"] is True
    assert summary["trajectory_atom_count"] == 6
    assert summary["selected_frame_count"] == 3
    assert summary["selected_time_range_ns"] == {"start": 0.0, "end": 2.0}
    assert summary["metrics"] == {
        "backbone_rmsd_mean_angstrom": pytest.approx(0.0, abs=1e-5),
        "backbone_rmsd_max_angstrom": pytest.approx(0.0, abs=1e-5),
        "ligand_rmsd_mean_angstrom": pytest.approx(2.0 / 3.0, abs=1e-5),
        "ligand_rmsd_max_angstrom": pytest.approx(2.0, abs=1e-5),
        "selected_frame_count": 3,
        "contact_count_mean": pytest.approx(1.0),
    }
    assert summary["source_segments"] == source_segments
    assert summary["source"]["trajectory"]["sha256"] == hashlib.sha256(trajectory_path.read_bytes()).hexdigest()
    assert summary["source"]["trajectory"]["path"] == "artifacts/legacy/production/production.xtc"
    with (output_dir / "timeseries.csv").open() as handle:
        import csv
        rows = list(csv.DictReader(handle))
    assert [float(row["time_ns"]) for row in rows] == pytest.approx([0.0, 1.0, 2.0])
    assert float(rows[0]["protein_aligned_backbone_rmsd_angstrom"]) == pytest.approx(0.0, abs=1e-5)
    assert float(rows[2]["protein_aligned_backbone_rmsd_angstrom"]) == pytest.approx(0.0, abs=1e-5)
    assert float(rows[2]["ligand_rmsd_angstrom"]) == pytest.approx(2.0, abs=1e-4)
    with (output_dir / "residue_contacts.csv").open() as handle:
        import csv
        contact = next(csv.DictReader(handle))
    assert contact["residue"] == "A:10:ALA"
    assert float(contact["contact_occupancy_fraction"]) == pytest.approx(1.0)
    assert (output_dir / "representative_frame.pdb").is_file()
    with pytest.raises(ValueError, match="no periodic box vectors"):
        analyze({
            "topology_path": str(topology_path),
            "trajectory_path": str(trajectory_path),
            "output_dir": str(tmp_path / "pbc-required"),
            "ligand_selection": "resname LIG",
            "protein_selection": "protein",
            "pbc_mode": "require",
        })


def test_openmm_dcd_header_drives_absolute_time_across_stride_and_chunks(tmp_path):
    md = pytest.importorskip("mdtraj")
    openmm = pytest.importorskip("openmm")
    from openmm import app, unit
    from openmm.app.dcdfile import DCDFile

    md_topology = md.Topology()
    md_chain = md_topology.add_chain("A")
    protein = md_topology.add_residue("ALA", md_chain, resSeq=10)
    md_atoms = [
        md_topology.add_atom(name, element, protein)
        for name, element in zip(
            ["N", "CA", "C", "O"],
            [md.element.nitrogen, md.element.carbon, md.element.carbon, md.element.oxygen],
        )
    ]
    for left, right in zip(md_atoms, md_atoms[1:]):
        md_topology.add_bond(left, right)
    ligand_residue = md_topology.add_residue("LIG", md_chain, resSeq=501)
    lig_atoms = [
        md_topology.add_atom("C1", md.element.carbon, ligand_residue),
        md_topology.add_atom("O1", md.element.oxygen, ligand_residue),
    ]
    md_topology.add_bond(*lig_atoms)
    base = np.asarray(
        [[0.00, 0.00, 0.00], [0.10, 0.00, 0.00], [0.20, 0.00, 0.00],
         [0.20, 0.10, 0.00], [0.35, 0.00, 0.00], [0.40, 0.00, 0.00]], dtype=np.float32,
    )
    frames = [base.copy() for _ in range(4)]
    frames[3][4:, 0] += 0.2
    topology_path = tmp_path / "system.pdb"
    md.Trajectory(frames[0:1], md_topology).save_pdb(str(topology_path))

    openmm_topology = app.Topology()
    chain = openmm_topology.addChain("A")
    omm_protein = openmm_topology.addResidue("ALA", chain, id="10")
    omm_atoms = [
        openmm_topology.addAtom(name, element, omm_protein)
        for name, element in zip(["N", "CA", "C", "O"], [app.element.nitrogen, app.element.carbon, app.element.carbon, app.element.oxygen])
    ]
    for left, right in zip(omm_atoms, omm_atoms[1:]):
        openmm_topology.addBond(left, right)
    omm_ligand = openmm_topology.addResidue("LIG", chain, id="501")
    omm_ligand_atoms = [
        openmm_topology.addAtom(name, app.element.carbon if name == "C1" else app.element.oxygen, omm_ligand)
        for name in ("C1", "O1")
    ]
    openmm_topology.addBond(*omm_ligand_atoms)
    trajectory_path = tmp_path / "production.dcd"
    with trajectory_path.open("wb") as handle:
        writer = DCDFile(
            handle, openmm_topology, 0.002 * unit.picoseconds,
            firstStep=500_000, interval=100_000,
        )
        for frame in frames:
            writer.writeModel(unit.Quantity(frame.tolist(), unit.nanometer))

    output_dir = tmp_path / "dcd-analysis"
    summary = analyze({
        "topology_path": str(topology_path), "trajectory_path": str(trajectory_path),
        "output_dir": str(output_dir), "ligand_selection": "resname LIG",
        "protein_selection": "protein", "stride": 2, "start_ns": 0.9,
        "end_ns": 1.4, "chunk_size": 2, "pbc_mode": "none",
    })
    assert summary["source_time"]["source"] == "charmm_cord_header"
    assert summary["source_time"]["first_step"] == 500_000
    assert summary["source_time"]["step_interval"] == 100_000
    assert summary["source_time"]["frame_interval_ps"] == pytest.approx(200.0, rel=1e-6)
    assert summary["selected_time_range_ns"] == {"start": pytest.approx(1.0), "end": pytest.approx(1.4)}
    assert summary["selected_frame_count"] == 2
    assert summary["trajectory_frames_seen"] == 4
    assert summary["trajectory_frames_seen_after_stride"] == 2
    with (output_dir / "timeseries.csv").open() as handle:
        import csv
        rows = list(csv.DictReader(handle))
    assert [float(row["time_ns"]) for row in rows] == pytest.approx([1.0, 1.4])


def test_dcd_timing_rejects_ambiguous_headers_without_complete_override(tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.tools.md_analysis import _dcd_timing

    path = tmp_path / "unknown.dcd"
    path.write_bytes(b"not a DCD header".ljust(100, b"\0"))
    with pytest.raises(ValueError, match="Unsupported DCD header"):
        _dcd_timing(path)
    with pytest.raises(ValueError, match="requires both"):
        _dcd_timing(path, frame_interval_ps=10.0)
    assert _dcd_timing(path, frame_interval_ps=10.0, first_time_ps=1000.0) == {
        "source": "explicit_override", "frame_interval_ps": 10.0, "first_time_ps": 1000.0,
    }


def test_dcd_timing_reads_openmm_header_first_step_and_interval(tmp_path):
    pytest.importorskip("openmm")
    from openmm import app, unit
    from openmm.app.dcdfile import DCDFile
    from mn_molecular_design_studio.engine.workflow_engine.tools.md_analysis import _dcd_timing

    topology = app.Topology()
    chain = topology.addChain("A")
    residue = topology.addResidue("ALA", chain)
    topology.addAtom("CA", app.element.carbon, residue)
    path = tmp_path / "header.dcd"
    with path.open("wb") as handle:
        writer = DCDFile(handle, topology, 0.004 * unit.picoseconds, firstStep=250_000, interval=2_500)
        writer.writeModel([[0, 0, 0]] * unit.nanometer)
    timing = _dcd_timing(path)
    assert timing["source"] == "charmm_cord_header"
    assert timing["nset"] == 1
    assert timing["first_step"] == 250_000
    assert timing["step_interval"] == 2_500
    assert timing["first_time_ps"] == pytest.approx(1000.0, rel=1e-6)
    assert timing["frame_interval_ps"] == pytest.approx(10.0, rel=1e-6)


def test_periodic_imaging_anchors_protein_and_images_wrapped_ligand(tmp_path):
    md = pytest.importorskip("mdtraj")
    topology = md.Topology()
    chain = topology.add_chain("A")
    residue = topology.add_residue("ALA", chain, resSeq=10)
    atoms = [
        topology.add_atom(name, element, residue)
        for name, element in zip(
            ["N", "CA", "C", "O"],
            [md.element.nitrogen, md.element.carbon, md.element.carbon, md.element.oxygen],
        )
    ]
    for left, right in zip(atoms, atoms[1:]):
        topology.add_bond(left, right)
    ligand_residue = topology.add_residue("LIG", chain, resSeq=501)
    ligand = [
        topology.add_atom("C1", md.element.carbon, ligand_residue),
        topology.add_atom("O1", md.element.oxygen, ligand_residue),
    ]
    topology.add_bond(*ligand)
    wrapped = np.asarray(
        [[0.05, 0.05, 0.05], [0.15, 0.05, 0.05], [0.25, 0.05, 0.05],
         [0.25, 0.15, 0.05], [3.95, 0.05, 0.05], [3.90, 0.05, 0.05]],
        dtype=np.float32,
    )
    xyz = np.stack((wrapped, wrapped.copy()))
    box = np.repeat((np.eye(3, dtype=np.float32) * 4.0)[None, :, :], 2, axis=0)
    trajectory = md.Trajectory(xyz=xyz, topology=topology, time=np.asarray([0.0, 100.0]))
    trajectory.unitcell_vectors = box
    topology_path = tmp_path / "system.pdb"
    trajectory[0].save_pdb(str(topology_path))
    trajectory_path = tmp_path / "wrapped.xtc"
    trajectory.save_xtc(str(trajectory_path))
    output_dir = tmp_path / "periodic-analysis"

    summary = analyze({
        "topology_path": str(topology_path), "trajectory_path": str(trajectory_path),
        "output_dir": str(output_dir), "ligand_selection": "resname LIG",
        "protein_selection": "protein", "pbc_mode": "require", "chunk_size": 1,
    })
    assert summary["parameters"]["periodic_unit_cell_present"] is True
    assert summary["parameters"]["pbc_imaging_applied"] is True
    with (output_dir / "timeseries.csv").open() as handle:
        import csv
        rows = list(csv.DictReader(handle))
    assert [float(row["protein_ligand_minimum_distance_angstrom"]) for row in rows] == pytest.approx([1.0, 1.0], abs=0.2)
    assert [int(row["geometric_residue_contacts"]) for row in rows] == [1, 1]


def test_source_paths_must_belong_to_the_same_managed_owner_run(tmp_path, monkeypatch):
    runs = tmp_path / "runs"
    owner = runs / "ligand-md" / "mnligand-md-test"
    topology = owner / "artifacts" / "legacy" / "prep" / "system.pdb"
    trajectory = owner / "artifacts" / "legacy" / "production" / "production.xtc"
    topology.parent.mkdir(parents=True)
    trajectory.parent.mkdir(parents=True)
    topology.write_text("ATOM\n")
    trajectory.write_bytes(b"frame-bytes")
    (owner / "metadata.json").write_text(json.dumps({"status": "completed"}))
    monkeypatch.setattr(md_analysis, "runs_root", lambda: runs)

    resolved = md_analysis._resolve_source_paths(owner, topology, trajectory)

    assert resolved == (owner.resolve(), topology.resolve(), trajectory.resolve())
    outsider = tmp_path / "outside.xtc"
    outsider.write_bytes(b"no")
    with pytest.raises(ValueError, match="inside the selected source"):
        md_analysis._resolve_source_paths(owner, topology, outsider)


def test_source_path_rejects_unfinished_run(tmp_path, monkeypatch):
    runs = tmp_path / "runs"
    owner = runs / "ligand-md" / "pending"
    artifacts = owner / "artifacts"
    artifacts.mkdir(parents=True)
    topology = artifacts / "system.pdb"
    trajectory = artifacts / "production.xtc"
    topology.write_text("ATOM\n")
    trajectory.write_bytes(b"frame")
    (owner / "metadata.json").write_text(json.dumps({"status": "running"}))
    monkeypatch.setattr(md_analysis, "runs_root", lambda: runs)

    with pytest.raises(ValueError, match="completed source MD run"):
        md_analysis._resolve_source_paths(owner, topology, trajectory)


@pytest.mark.parametrize(
    ("protein", "ligand", "residues", "cutoff", "message"),
    [
        (np.empty((0, 3)), np.ones((1, 3)), np.empty(0, dtype=int), 4.5, "Protein contact selection"),
        (np.ones((1, 3)), np.empty((0, 3)), np.zeros(1, dtype=int), 4.5, "Ligand contact selection"),
        (np.ones((1, 3)), np.ones((1, 3)), np.empty(0, dtype=int), 4.5, "residue mapping"),
    ],
)
def test_contact_metrics_rejects_invalid_atom_identity(protein, ligand, residues, cutoff, message):
    with pytest.raises(ValueError, match=message):
        residue_contact_metrics(protein, ligand, residues, cutoff)
