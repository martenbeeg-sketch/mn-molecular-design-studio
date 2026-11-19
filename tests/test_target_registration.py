import hashlib
import json

import pytest

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.data_browser import folder_scopes, saved_results
from mn_molecular_design_studio.database_browser import build_browser_index, subject_results
from mn_molecular_design_studio.molecular_database_ingest import ingest_run
from mn_molecular_design_studio.target_registration import reference_context, saved_reference_context
from test_molecular_database_policy import A, B, database, _make_run, _pdb_for_sequence, _pdb_for_two_sequences


def _inputs(run, **params):
    path = run / "input.json"
    document = json.loads(path.read_text())
    document["params"].update(params)
    path.write_text(json.dumps(document))


def _candidate(tool, *, stage="generation.backbone_sequence", passed=True, raw=None):
    return {"candidate_id": "same-native-name", "source_tool": tool, "stage": stage,
        "binder_sequence": "ACD", "binder_chains": ["A"], "target_chains": ["B"],
        "complex_pdb": "complex.pdb", "metrics": {"confidence": 0.9, **({"native_pass_filters": True} if passed else {})},
        "raw_metadata": raw or {}}


def _design(database, tool="rfd4_proteina", *, pdb="4LNW", passed=True):
    run = _make_run(database, "design_generation", tool + "-" + pdb, tool=tool, workflow="design_generation",
        candidates=[_candidate(tool, passed=passed)])
    (run / "complex.pdb").write_text(_pdb_for_two_sequences("ACD", "EFG"))
    _inputs(run, target_pdb_id=pdb, target_receptor="TRα" if pdb == "4LNW" else "TRβ", target_name="Research target")
    return run


def _outcomes(database, pdb):
    index = build_browser_index(database, A, collection="library", include_related=True)
    scopes = folder_scopes(index, [])
    scope = scopes["targets:pdb:" + pdb]
    return index, saved_results(index, subject_results(database, index, scope["subject"]))


@pytest.mark.parametrize("tool", ["rfd4_proteina", "bindcraft2", "boltzgen", "rfantibody", "esmfold2_binder_design"])
def test_native_designs_share_target_registration_and_admission(database, tool):
    run = _design(database, tool)
    assert ingest_run(run, database=database) > 0
    assert ingest_run(run, database=database) == 0
    index, outcomes = _outcomes(database, "4LNW")
    designs = [row for row in outcomes if row["outcome_kind"] == "Design"]
    assert len(designs) == 1
    assert designs[0]["eligibility"] == "Passed"
    assert designs[0]["target_pdb_id"] == "4LNW"
    target = index["molecules"][designs[0]["target_id"]]
    assert [row["sequence"] for row in target["components"]] == ["EFG"]
    with database.connection() as connection:
        contexts = [json.loads(row[0]) for row in connection.execute(
            "SELECT context FROM molecule_structures WHERE molecule_id=?", (designs[0]["molecule_id"],))]
    assert contexts and all(row["target_pdb_id"] == "4LNW" for row in contexts)
    assert all(row["target_molecule_id"] == designs[0]["target_id"] for row in contexts)


def test_identical_sequences_keep_different_recorded_pdb_studies(database):
    for pdb in ("4LNW", "3GWS"):
        ingest_run(_design(database, pdb=pdb), database=database)
    _, alpha = _outcomes(database, "4LNW")
    _, beta = _outcomes(database, "3GWS")
    assert len(alpha) == len(beta) == 1
    assert alpha[0]["target_id"] == beta[0]["target_id"]
    assert {row["run_id"] for row in alpha[0]["originals"]} == {"rfd4_proteina-4LNW"}
    assert {row["run_id"] for row in beta[0]["originals"]} == {"rfd4_proteina-3GWS"}


@pytest.mark.parametrize("tool", ["optimized_prediction_boltz2", "af3_jax_cofolding", "chai1", "esmfold"])
def test_independent_predictions_use_the_same_reference_contract(database, tool):
    from mn_molecular_design_studio.molecular_database_results import list_run_molecules, save_run_selection

    run = _make_run(database, "refolding", tool, tool=tool, workflow="refolding.optimized_prediction",
        candidates=[_candidate(tool, stage="complex_refolding", passed=False)])
    _inputs(run, target_pdb_id="3GWS", target_name="TRβ")
    (run / "complex.pdb").write_text(_pdb_for_two_sequences("ACD", "EFG"))
    ingest_run(run, database=database)
    rows = list_run_molecules(run, A, database=database)
    chosen = next(row for row in rows if row.get("representations"))
    save_run_selection(run, A, [chosen["id"]], database=database)
    _, outcomes = _outcomes(database, "3GWS")
    assert len(outcomes) == 1
    assert outcomes[0]["activity"] == "Cofolding"
    assert outcomes[0]["eligibility"] in {"Recorded", "Manually selected"}
    assert outcomes[0]["primary"]["receipt"]["role"] == "independent_prediction"


@pytest.mark.parametrize("tool", ["vina", "gnina", "unidock_pro", "autodock_gpu", "rosetta_ligand", "nesso"])
def test_native_ligand_tools_preserve_reference_on_scores_and_poses(database, tool):
    from mn_molecular_design_studio.ligand_compounds import write_compound_set
    from rdkit import Chem
    from rdkit.Chem import AllChem

    database.import_records(A, molecule_imports.parse_smiles("CCO Ethanol"), name="Ligands", project_id="project-1")
    affinity = tool == "nesso"
    workflow = "ligand.affinity.nesso" if affinity else "ligand.docking." + tool
    run = _make_run(database, "ligand-affinity" if affinity else "ligand-docking", tool, tool=tool, workflow=workflow)
    _inputs(run, target_pdb_id="4LNW", target_receptor="TRα", target_sequence="EFG")
    inputs = run / "artifacts/input"
    inputs.mkdir()
    write_compound_set(inputs / "compound_set.json", name="Consumed", compounds=[{"compound_id": "ethanol", "smiles": "CCO"}],
        owner_id=A["id"], project_id="project-1", created_at="now")
    if affinity:
        for filename in ("metadata.json", "input.json"):
            document = json.loads((run / filename).read_text())
            document.pop("workflow")
            (run / filename).write_text(json.dumps(document))
        (run / "artifact_contract.json").write_text(json.dumps({"workflow": workflow}))
        (run / "artifacts/nesso_affinity.csv").write_text("compound_id\taffinity_log10_ic50_uM\nethanol\t0.7\n")
    else:
        (inputs / "receptor.pdb").write_text(_pdb_for_sequence("EFG"))
        molecule = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        AllChem.EmbedMolecule(molecule, randomSeed=12)
        with Chem.SDWriter(str(run / "pose.sdf")) as writer:
            writer.write(molecule)
        (run / "artifacts/docking_scores.csv").write_text(
            "compound_id\tstatus\tbest_score_kcal_mol\tpose_path\nethanol\tcompleted\t-7.5\tpose.sdf\n")
    ingest_run(run, database=database)
    _, outcomes = _outcomes(database, "4LNW")
    assert len(outcomes) == 1
    assert outcomes[0]["activity"] == ("Binding affinity" if affinity else "Docking")
    assert outcomes[0]["eligibility"] == "Recorded"
    assert outcomes[0]["primary"]["receipt"]["role"] == ("binding_affinity" if affinity else "docking")
    assert outcomes[0]["structure_count"] == (0 if affinity else 1)
    assert ingest_run(run, database=database) == 0


def test_explicit_later_selection_registers_target_without_inventing_a_pass(database):
    from mn_molecular_design_studio.molecular_database_results import list_run_molecules, save_run_selection

    run = _design(database, passed=False)
    assert ingest_run(run, database=database) == 0
    rows = list_run_molecules(run, A, database=database)
    assert len(rows) == 1 and rows[0]["result_only"]
    saved = save_run_selection(run, A, [rows[0]["id"]], database=database)
    assert saved["saved_count"] == 1
    _, outcomes = _outcomes(database, "4LNW")
    assert len(outcomes) == 1
    assert outcomes[0]["eligibility"] == "Manually selected"
    assert outcomes[0]["target_id"]


@pytest.mark.parametrize("tool", ["biotite_interface", "ipsae", "pep_patch"])
def test_post_analysis_inherits_selected_structure_reference(database, tool):
    run = _design(database)
    ingest_run(run, database=database)
    _, outcomes = _outcomes(database, "4LNW")
    design = outcomes[0]
    record = database.list_molecules(A, molecule_ids=[design["molecule_id"]], collection="all")["items"][0]
    structure = next(row for row in record["representations"] if row["context"].get("workflow_role") == "design_generation")
    content = database.structure_file(A, structure["id"]).read_bytes()
    source = {"schema": "studio.complex_analysis_database_source.v1", "molecule_id": record["id"],
        "group_id": structure["group_id"], "structure_id": structure["id"],
        "structure_digest": hashlib.sha256(content).hexdigest(), "binder_chains": ["A"], "target_chains": ["B"]}
    analysis = _make_run(database, "analysis", tool, tool=tool, workflow="analysis." + tool,
        candidates=[_candidate(tool, stage="analysis", passed=False, raw={"database_source": source})])
    (analysis / "complex.pdb").write_bytes(content)
    ingest_run(analysis, database=database)
    with database.connection() as connection:
        context = json.loads(connection.execute("SELECT context FROM molecule_observations WHERE namespace=?", (tool,)).fetchone()[0])
    assert context["target_pdb_id"] == "4LNW"
    assert context["target_molecule_id"] == design["target_id"]
    _, after = _outcomes(database, "4LNW")
    assert len(after) == 1 and after[0]["outcome_kind"] == "Design"
    assert tool in {row["tool"] for row in after[0]["attached"]}


def test_reference_comes_from_accession_metadata_not_names_or_sequence(database, tmp_path):
    path = tmp_path / "3GWS-misleading.pdb"
    path.write_text("HEADER".ljust(62) + "4LNW\n" + _pdb_for_sequence("EFG"))
    assert reference_context({"target_name": "3GWS"}, path=path)["target_pdb_id"] == "4LNW"
    conflict = reference_context({"target_pdb_id": "3GWS"}, path=path)
    assert conflict["target_pdb_id"] == ""
    assert conflict["target_reference"]["status"] == "conflicting"
    assert reference_context(conflict)["target_reference"]["pdb_ids"] == ["3GWS", "4LNW"]
    assert reference_context({"pdb_id": "STUDIO", "target_reference": {"pdb_ids": 12}})["target_reference"]["status"] == "unspecified"


def test_selected_sql_structure_overrides_other_reference_on_same_identity(database):
    record = {"name": "Same target", "components": [{"type": "protein", "chain": "A", "sequence": "EFG"}],
        "metadata": {"target_pdb_id": "4LNW"}}
    result = database.import_records(A, [record], name="Target", project_id="project-1")
    structure = database.attach_structure(result["group_id"], result["molecule_ids"][0],
        {"name": "reference.pdb", "content": _pdb_for_sequence("EFG").encode()}, context={"target_pdb_id": "3GWS"})
    reference = saved_reference_context(database, A, "project-1", structure_id=structure)
    assert reference["target_pdb_id"] == "3GWS"
    assert saved_reference_context(database, B, "project-1", structure_id=structure)["target_pdb_id"] == ""
    assert saved_reference_context(database, A, "another-project", structure_id=structure)["target_pdb_id"] == ""


def test_small_molecule_targets_keep_chemical_identity(database):
    from mn_molecular_design_studio.molecular_database_ingest import _register_candidate_target

    target = {"type": "small_molecule", "chain": "L", "role": "target", "smiles": "CCO"}
    candidate = {"binder_chains": ["A"], "target_chains": ["L"]}
    result = _register_candidate_target(database, A, "project-1", database.root, {}, {}, candidate,
        {"components": [{"type": "protein", "chain": "A", "role": "binder", "sequence": "ACD"}, target]},
        source_key="ligand-target", origin={})
    target_record = database.list_molecules(A, molecule_ids=[result["target_molecule_id"]], collection="all")["items"][0]
    assert target_record["type"] == "small_molecule"
    assert target_record["components"][0]["smiles"] == "CCO"


@pytest.mark.parametrize("row,recorded", [
    ({"affinity_log10_ic50_uM": "0.0"}, True),
    ({"affinity_log10_ic50_uM": "nan", "seed": "42"}, False),
    ({"status": "failed", "affinity_log10_ic50_uM": "0.7"}, False),
    ({"attempt_index": "1", "binder_probability": "inf"}, False),
])
def test_affinity_receipt_requires_a_finite_measurement_not_job_metadata(row, recorded):
    from mn_molecular_design_studio.molecular_database_ingest import _calculation_receipt

    assert _calculation_receipt(row, affinity=True, default_status="completed")["calculation_recorded"] is recorded
