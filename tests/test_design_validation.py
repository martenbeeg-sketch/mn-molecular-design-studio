"""Design eligibility, native template serialization, and RMSD frame contracts."""
import copy
import gzip
import json
from pathlib import Path

import numpy as np
import pytest
import yaml
from Bio.PDB import MMCIFIO

from mn_molecular_design_studio import design_validation as validation, structure_templates
from mn_molecular_design_studio import cofolding_msa
from mn_molecular_design_studio.engine.workflow_engine.core.candidates import _parse_structure
from mn_molecular_design_studio.engine.workflow_engine.workflows import optimized_prediction as prediction


@pytest.mark.parametrize("metrics,raw,expected", [
    ({"native_pass_filters": True}, {}, True),
    ({"pass_filters": "True"}, {}, True),
    ({"native_pass_filters": 1}, {}, True),
    ({"native_pass_filters": False, "pass_filters": True}, {}, False),
    ({"pass_filters": "false"}, {}, False),
    ({"native_pass_filters": None}, {}, False),
    ({"pass_plddt_filter": True}, {}, False),
    ({}, {}, False),
    ({"native_pass_filters": True}, {"candidate_pool_level": "prefilter_generated"}, False),
    ({"native_pass_filters": True}, {"result_kind": "generation_only"}, False),
    ({}, {"bindcraft_output_kind": "Accepted"}, True),
    ({"pass_filters": False}, {"bindcraft_output_kind": "Accepted"}, False),
])
def test_only_affirmative_overall_filter_passes_are_eligible(metrics, raw, expected):
    assert validation.candidate_passes_filters({"source_tool": "bindcraft", "metrics": metrics, "raw_metadata": raw}) is expected


@pytest.mark.parametrize("status,success,passed,expected", [
    ("completed", True, True, 1), ("succeeded", True, True, 1),
    ("failed", False, True, 0), ("failed", True, True, 0),
    ("completed", False, True, 0), ("running", True, True, 0),
    ("completed", True, False, 0), ("completed", True, None, 0),
])
def test_source_job_must_succeed_and_contain_passing_designs(design_source_factory, status, success, passed, expected):
    run, _ = design_source_factory(status=status, success=success, passed=passed)
    assert len(validation.design_candidates(run)) == expected


@pytest.mark.parametrize("compressed", [False, True])
def test_bindcraft_cif_without_occupancy_can_be_used(design_source_factory, compressed):
    run, _ = design_source_factory()
    path = run / "designed.cif"
    if compressed:
        target = path.with_suffix(".cif.gz"); target.write_bytes(gzip.compress(path.read_bytes())); path = target
    assert {c: len(r) for c, r in validation.protein_chains(validation.read_structure(path)).items()} == {"A": 6, "B": 5}


def test_candidate_reader_accepts_missing_bfactor_and_occupancy(tmp_path):
    cif = b"""data_model
#
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.label_entity_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
_atom_site.id
ATOM N N . ALA A 1 . 1 ALA A N 1 0.0 0.0 0.0 1 1
ATOM C CA . ALA A 1 . 1 ALA A CA 1 1.4 0.0 0.0 1 2
ATOM C C . ALA A 1 . 1 ALA A C 1 2.8 0.0 0.0 1 3
ATOM O O . ALA A 1 . 1 ALA A O 1 3.3 1.2 0.0 1 4
#
"""
    path = tmp_path / "generated.cif"
    path.write_bytes(cif)

    structure = _parse_structure(path)
    model = next(structure.get_models())
    atoms = list(model["A"].get_atoms())

    assert validation.sequence(validation.protein_chains(structure)["A"]) == "A"
    assert np.allclose(atoms[1].coord, [1.4, 0.0, 0.0])
    assert np.isnan(atoms[0].bfactor)
    assert path.read_bytes() == cif


def test_candidate_reader_still_rejects_malformed_mmcif_coordinates(tmp_path):
    path = tmp_path / "malformed.cif"
    path.write_text("""data_model
#
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.label_entity_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
_atom_site.id
ATOM N N . ALA A 1 . 1 ALA A N 1 invalid 0.0 0.0 1 1
#
""")

    with pytest.raises((ValueError, KeyError)):
        _parse_structure(path)


@pytest.mark.parametrize("mode", ["complex", "binder"])
def test_selected_design_has_target_only_template_and_no_binder_msa(design_source_factory, mode):
    run, candidate = design_source_factory()
    prepared = prediction._prepare_prediction(
        "boltz2", validation_sources=[{"source_run_dir": str(run), "candidate_ids": [candidate["candidate_id"]]}],
        validation_mode=mode, resolve_random_seed=False)
    task = prepared.normalized_tasks[0]
    assert task["binder_protein_ids"] == (["B"] if mode == "complex" else ["A"])
    assert not task["molecule_chain_options"][task["binder_protein_ids"][0]]["use_msa"]
    assert bool(task["structure_templates"]) is (mode == "complex")
    if mode == "complex":
        assert task["structure_templates"][0]["query_chain_id"] == "A"
        assert task["structure_templates"][0]["template_chain_id"] == "A"
    else:
        assert prepared.selected_options["use_msa"] is False
        assert task["target_protein_ids"] == []


def test_cyclic_acceptance_validation_preserves_topology_only_on_supported_engine(design_source_factory):
    run, candidate = design_source_factory()
    source = [{"source_run_dir": str(run), "candidate_ids": [candidate["candidate_id"]]}]
    prepared = prediction._prepare_prediction(
        "boltz2", validation_sources=source, validation_mode="complex",
        validation_binder_format="cyclic_peptide", samples=1, resolve_random_seed=False,
    )
    binder = prepared.normalized_tasks[0]["binder_protein_ids"][0]
    assert prepared.normalized_tasks[0]["molecule_chain_options"][binder]["cyclic"] is True
    with pytest.raises(ValueError, match="does not preserve cyclic peptide"):
        prediction._prepare_prediction(
            "alphafold3", validation_sources=source, validation_mode="complex",
            validation_binder_format="cyclic_peptide", samples=1,
            accept_alphafold3_terms=True, resolve_random_seed=False,
        )


def test_failed_filter_candidate_cannot_be_selected_by_id(design_source_factory):
    run, candidate = design_source_factory(passed=False)
    with pytest.raises(ValueError, match="filter-passing"):
        validation.build_tasks([{"source_run_dir": str(run), "candidate_ids": [candidate["candidate_id"]]}],
                               mode="complex", profile=prediction.ENGINE_PROFILES["boltz2"],
                               target_use_msa=True, use_target_template=True)


def test_standalone_acceptance_can_repredict_native_rejections(design_source_factory):
    run, candidate = design_source_factory(passed=False)
    source = [{"source_run_dir": str(run), "candidate_ids": [candidate["candidate_id"]]}]
    with pytest.raises(ValueError, match="filter-passing"):
        prediction._prepare_prediction(
            "boltz2", validation_sources=source, validation_mode="complex", samples=1,
            resolve_random_seed=False)
    prepared = prediction._prepare_prediction(
        "boltz2", validation_sources=source, validation_mode="complex", samples=1,
        validation_require_passing=False, resolve_random_seed=False)
    assert prepared.normalized_tasks[0]["validation_reference"]["candidate_id"] == candidate["candidate_id"]


@pytest.mark.parametrize("engine", [name for name, profile in prediction.ENGINE_PROFILES.items()
                                   if profile.get("supports_msa") and profile["input_kind"] != "boltzgen"])
def test_only_target_msa_is_searched_and_binder_uses_query_only(tmp_path, monkeypatch, engine):
    calls = []
    def prepare(run_dir, sequence, chain_id, *, item_id, **kwargs):
        calls.append(chain_id)
        path = run_dir / "artifacts/raw/cofolding_msa" / item_id / (kwargs.get("staged_name") or f"chain_{chain_id}.a3m")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f">query\n{sequence}\n>homolog\n{sequence}\n")
        return path
    monkeypatch.setattr(prediction, "reference_root", lambda: tmp_path / "reference")
    monkeypatch.setattr(cofolding_msa, "prepare_chain", prepare)
    run = tmp_path / "job"
    path, details = prediction._stage_native_input(
        run, profile=prediction.ENGINE_PROFILES[engine], item_id="item", job_name="example",
        proteins=[{"id": "A", "sequence": "ACDEFG"}, {"id": "B", "sequence": "KLMNP"}],
        dna=[], rna=[], ligand_ccd_codes=[], ligand_smiles=[], seed=1, samples=1,
        target_pdb_text="", target_chain_id="A", target_msa_a3m="", binder_sequence="KLMNP", num_recycles=3,
        molecule_chain_options={"A": {"use_msa": True}, "B": {"use_msa": False}},
        prepare_msas=True, use_msa=True)
    assert calls == ["A"]
    assert details["query_only_chains"] == ["B"]
    assert len(cofolding_msa.records((run / details["msa_files_by_chain"]["B"]).read_text())) == 1


@pytest.mark.parametrize("mode,expected", [("complex", 4.), ("binder", 0.)])
def test_complex_uses_target_fit_and_binder_alone_uses_binder_fit(design_source_factory, tmp_path, mode, expected):
    run, _ = design_source_factory()
    structure = validation.read_structure(run / "designed.cif")
    predicted = copy.deepcopy(structure)
    rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    for chain in next(predicted.get_models()):
        for atom in chain.get_atoms():
            shifted = atom.coord + (np.array([0., 0., 4.]) if chain.id == "B" else 0.)
            atom.coord = shifted @ rotation + np.array([12., -3., 8.])
    if mode == "binder":
        model = next(predicted.get_models()); model.detach_child("A"); model["B"].id = "A"
    prediction_path = tmp_path / "predicted.cif"
    writer = MMCIFIO(); writer.set_structure(predicted); writer.save(str(prediction_path))
    reference = {"reference_structure_path": str(run / "designed.cif"), "source_run_dir": str(run),
                 "candidate_id": "passing-design", "mode": mode, "target_chains": ["A"] if mode == "complex" else [],
                 "binder_chains": ["B"], "chain_map": {"A": "A", "B": "B"} if mode == "complex" else {"B": "A"}}
    metrics, artifacts = validation.compare_prediction(tmp_path / "validation", "item", 1, prediction_path, reference)
    key = "binder_rmsd_target_aligned" if mode == "complex" else "binder_rmsd_fitted"
    assert metrics[key] == pytest.approx(expected, abs=1e-5)
    assert metrics["binder_rmsd_coverage"] == 1.
    assert metrics["binder_rmsd_ca_count"] == 5
    assert Path(artifacts["aligned_structure_path"]).is_file()
    if mode == "complex":
        assert metrics["target_alignment_rmsd"] == pytest.approx(0., abs=1e-5)
    else:
        assert metrics["monomer_refolding_success"] is True
    from mn_molecular_design_studio.design_validation_ui import comparison_viewer
    viewer = comparison_viewer(run / "designed.cif", Path(artifacts["aligned_structure_path"]), reference,
                               artifacts["validation_chain_map"])
    assert [o["name"] for o in viewer.objects] == ["Designed reference", "Refolded prediction"]
    assert [len(o["frames"]) for o in viewer.objects] == [1, 1]
    assert [len(o["frames"][0]["coords"]) for o in viewer.objects] == ([11, 11] if mode == "complex" else [5, 5])


@pytest.mark.parametrize("engine", sorted(structure_templates.TEMPLATE_ENGINES - {"boltzgen_fold"}))
def test_template_upload_is_serialized_for_each_supported_native_adapter(design_source_factory, tmp_path, engine):
    run, _ = design_source_factory()
    upload = {"name": "target.cif", "text": (run / "designed.cif").read_text(),
              "query_chain_id": "A", "template_chain_id": "A"}
    profile = prediction.ENGINE_PROFILES[engine]
    job = tmp_path / "prediction"
    path, details = prediction._stage_native_input(
        job, profile=profile, item_id="item", job_name="example", proteins=[{"id": "A", "sequence": "ACDEFG"}],
        dna=[], rna=[], ligand_ccd_codes=[], ligand_smiles=[], seed=1, samples=1,
        target_pdb_text="", target_chain_id="A", target_msa_a3m="", binder_sequence="", num_recycles=3,
        prepare_msas=False, use_msa=False, structure_uploads=[upload])
    assert details["template_count"] == 1
    kind = profile["input_kind"]
    if kind == "af3_json":
        templates = json.loads(path.read_text())["sequences"][0]["protein"]["templates"]
        assert templates[0]["queryIndices"] == list(range(6))
    elif kind == "boltz_yaml":
        assert yaml.safe_load(path.read_text())["templates"][0]["chain_id"] == ["A"]
    elif kind == "protenix_json":
        assert json.loads(path.read_text())[0]["sequences"][0]["proteinChain"]["templatesPath"].startswith("/work/")
    elif kind == "openfold3_ob0_json":
        assert json.loads(path.read_text())["queries"]["example"]["chains"][0]["template_cif_chain_ids"] == ["A"]
    elif kind == "rf3_json":
        assert json.loads(path.read_text())["template_selection"] == ["A"]
    elif kind == "multimer_fasta":
        assert details["batch_template_path"].endswith("complex_template.pdb")
    elif kind == "chai1_json":
        assert details["chai_template_hits"].endswith("template_hits.m8")
        assert ">protein|name=item_A" in path.with_suffix(".fasta").read_text()
        import gemmi
        cif = next(job.glob("artifacts/raw/optimized_prediction/templates/item/chain_A_0.cif"))
        assert len(gemmi.read_structure(str(cif))[0]["A"].get_polymer()) == 6
