from __future__ import annotations

import json
import hashlib
from pathlib import Path

import pytest

from mn_molecular_design_studio import molecule_imports
from mn_molecular_design_studio.molecular_database import MolecularDatabase
from mn_molecular_design_studio.molecular_database_ingest import ingest_run


A = {"id": "a" * 32, "role": "user", "_cluster_project_ids": ["project-1"]}
B = {"id": "b" * 32, "role": "user", "_cluster_project_ids": []}


@pytest.fixture
def database(tmp_path, monkeypatch):
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_HOME", str(tmp_path))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path))
    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_RUN_DIR", str(tmp_path / "runs"))
    return MolecularDatabase(tmp_path / "database")


def _make_run(database, task, run_id, *, owner=A["id"], project="project-1", tool="boltz2",
              workflow="refolding.optimized_prediction", status="completed", success=True, candidates=(), metrics=None,
              collection_policy="automatic"):
    from mn_molecular_design_studio.workspace import STUDIO_SOURCE_ID

    run_dir = database.root.parent / "runs" / task / run_id
    (run_dir / "artifacts/normalized_candidates").mkdir(parents=True, exist_ok=True)
    (run_dir / "molecular_design_studio.json").write_text(json.dumps({
        "application": STUDIO_SOURCE_ID, "studio_owner_id": owner, "project_id": project,
    }))
    metadata = {
        "studio_owner_id": owner, "project_id": project, "task_group": task,
        "tool": tool, "workflow": workflow, "status": status, "campaign_name": run_id,
    }
    if collection_policy is not None:
        metadata["database_collection_policy"] = collection_policy
    (run_dir / "metadata.json").write_text(json.dumps(metadata))
    (run_dir / "input.json").write_text(json.dumps({"tool": tool, "workflow": workflow, "params": {}}))
    (run_dir / "result.json").write_text(json.dumps({"success": success, "metrics": metrics or {}}))
    if candidates:
        (run_dir / "artifacts/normalized_candidates/candidates.jsonl").write_text(
            "".join(json.dumps(candidate) + "\n" for candidate in candidates))
    from mn_molecular_design_studio.job_catalog import reconcile
    reconcile()
    return run_dir


def _pdb_acd():
    lines = []
    for serial, (residue, x) in enumerate((("ALA", 1.0), ("CYS", 2.0), ("ASP", 3.0)), start=1):
        lines.append(f"ATOM  {serial:5d}  CA  {residue} A{serial:4d}    {x:8.3f}{1.0:8.3f}{1.0:8.3f}{1.00:6.2f}{20.00:6.2f}           C")
    return "\n".join([*lines, "TER", "END", ""])


def _pdb_efg_acd_complex():
    residue_names = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY"}
    lines = []
    serial = 1
    for chain, sequence, x_offset in (("A", "EFG", 0.0), ("B", "ACD", 8.0)):
        for residue_index, amino_acid in enumerate(sequence, start=1):
            x = x_offset + residue_index
            lines.append(
                f"ATOM  {serial:5d}  CA  {residue_names[amino_acid]} {chain}{residue_index:4d}    "
                f"{x:8.3f}{1.0:8.3f}{1.0:8.3f}{1.00:6.2f}{20.00:6.2f}           C"
            )
            serial += 1
        lines.append("TER")
    return "\n".join([*lines, "END", ""])


def _pdb_acd_efg_complex():
    text = _pdb_efg_acd_complex()
    return "\n".join([line[:21] + ("B" if line[21:22] == "A" else "A") + line[22:]
                       if line.startswith("ATOM") else line for line in text.splitlines()]) + "\n"


def _pdb_for_sequence(sequence, chain="A", x_offset=0.0):
    residue_names = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
                     "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
                     "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
                     "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}
    lines = []
    for serial, residue in enumerate(sequence, start=1):
        if residue not in residue_names:
            raise ValueError(f"No fixture residue mapping for {residue}.")
        lines.append(f"ATOM  {serial:5d}  CA  {residue_names[residue]} {chain}{serial:4d}    "
                     f"{x_offset + serial:8.3f}{1.0:8.3f}{1.0:8.3f}{1.00:6.2f}{20.00:6.2f}           C")
    return "\n".join([*lines, "TER", "END", ""])


def _pdb_for_two_sequences(first, second):
    residue_names = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
                     "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
                     "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
                     "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}
    lines, serial = [], 1
    for chain, sequence, offset in (("A", first, 0.0), ("B", second, 10.0)):
        for index, residue in enumerate(sequence, 1):
            lines.append(f"ATOM  {serial:5d}  CA  {residue_names[residue]} {chain}{index:4d}    "
                         f"{offset + index:8.3f}{1.0:8.3f}{1.0:8.3f}{1.00:6.2f}{20.00:6.2f}           C")
            serial += 1
        lines.append("TER")
    return "\n".join([*lines, "END", ""])


def _database_source(database, sequence):
    content = _pdb_for_sequence(sequence).encode()
    record = {"name": "Source molecule", "components": molecule_imports.parse_sequences(sequence)[0]["components"],
              "files": [{"name": "source.pdb", "content": content}]}
    imported = database.import_records(A, [record], name="Source molecule", project_id="project-1")
    molecule_id, group_id = imported["molecule_ids"][0], imported["group_id"]
    with database.connection() as connection:
        structure = connection.execute(
            "SELECT id,path FROM molecule_structures WHERE group_id=? AND molecule_id=?", (group_id, molecule_id)
        ).fetchone()
    digest = hashlib.sha256(Path(structure["path"]).read_bytes()).hexdigest()
    return ({"schema": "studio.complex_analysis_database_source.v1", "molecule_id": molecule_id,
             "structure_id": structure["id"], "group_id": group_id, "structure_digest": digest},
            molecule_id, group_id, structure["id"], digest)


def _pdb_acd_with_generic_ligand():
    ligand = [
        "HETATM   10  C1  LIG B   1       8.000   1.000   1.000  1.00 20.00           C",
        "HETATM   11  C2  LIG B   1       9.000   1.000   1.000  1.00 20.00           C",
    ]
    return _pdb_acd().replace("TER\nEND", "TER\n" + "\n".join(ligand) + "\nEND")


def test_library_collection_is_curated_by_origin_and_accumulates_accessible_scores(database):
    imported = database.import_records(A, molecule_imports.parse_smiles("CCO ethanol"), name="Imported")
    molecule_id = imported["molecule_ids"][0]
    old_auto = database.import_records(A, molecule_imports.parse_smiles("CCO ethanol"), name="Old auto",
                                       project_id="project-1", import_method="tool_output", origin={"type": "job", "run_id": "old"})
    database.add_observation(old_auto["group_id"], molecule_id, "vina", "old-run", {"energy": -7.1})
    generated = database.import_records(A, molecule_imports.parse_smiles("CCN generated"), name="Generated",
                                        project_id="project-1", import_method="tool_output",
                                        origin={"type": "job", "collection_kind": "generated"})

    library = database.groups(A, collection="library")
    all_groups = database.groups(A, collection="all")
    assert {row["id"] for row in library} == {imported["group_id"], generated["group_id"]}
    assert {row["id"] for row in all_groups} >= {imported["group_id"], old_auto["group_id"], generated["group_id"]}
    row = database.list_molecules(A, group_id=imported["group_id"], collection="library")["items"][0]
    assert row["scores"] == {"vina.energy": -7.1}
    assert database.metric_fields(A, group_id=imported["group_id"], collection="library") == ["vina.energy"]


def test_failed_or_empty_runs_and_raw_backbones_do_not_create_groups(database):
    empty = _make_run(database, "empty", "no-output", workflow="design_generation", tool="rfdiffusion3", candidates=())
    assert ingest_run(empty, database=database) == 0
    failed = _make_run(database, "failed", "candidate", workflow="design_generation", tool="bindcraft2", status="failed", success=False,
                       candidates=[{"candidate_id": "failed-1", "stage": "sequence_design", "binder_sequence": "ACDE"}])
    assert ingest_run(failed, database=database) == 0

    backbone = _make_run(database, "design_generation", "backbone", workflow="design_generation", tool="rfdiffusion3",
                         candidates=[{"candidate_id": "bb-1", "stage": "generation.backbone", "binder_pdb": "backbone.pdb"}])
    (backbone / "backbone.pdb").write_text(_pdb_acd())
    assert ingest_run(backbone, database=database) == 0
    assert database.groups(A, collection="library") == []
    assert database.groups(A, collection="all") == []


def test_policyless_legacy_generation_is_results_only(database):
    from mn_molecular_design_studio.collection_policy import policy_for_run
    run = _make_run(database, "design", "legacy-no-policy", workflow="design_generation", tool="boltzgen",
                    collection_policy=None, candidates=[{
                        "candidate_id": "accepted-looking", "stage": "complex_refolding",
                        "binder_sequence": "ACDE", "metrics": {"native_pass_filters": True},
                    }])
    assert policy_for_run(run) == "manual"
    assert ingest_run(run, database=database) == 0
    assert database.groups(A, collection="all") == []


def test_unknown_generated_candidate_is_not_retained_as_accepted(database):
    run = _make_run(database, "design", "unknown-native-status", workflow="design_generation", tool="boltzgen",
                    candidates=[{"candidate_id": "unknown", "stage": "complex_refolding", "binder_sequence": "ACDE"}])
    assert ingest_run(run, database=database) == 0
    assert database.groups(A, collection="all") == []


@pytest.mark.parametrize("workflow,tool,task", [
    ("analysis.biotite_interface", "biotite_interface", "analysis"),
    ("analysis.plip", "plip", "analysis"),
    ("analysis.pep_patch", "pep_patch", "analysis"),
    ("analysis.ipsae", "ipsae", "analysis"),
    ("analysis.rosetta", "rosetta", "analysis"),
    ("pipeline.analysis.plip", "plip", "pipeline"),
])
def test_real_analysis_workflows_without_exact_sources_create_no_library_records(database, workflow, tool, task):
    run = _make_run(database, task, "unresolved-" + tool + "-" + task, tool=tool, workflow=workflow,
        candidates=[{"candidate_id": "unresolved-analysis", "stage": "analysis", "source_tool": tool,
                    "binder_sequence": "ACD", "complex_pdb": "analysis.pdb", "metrics": {"score": 0.4}}])
    (run / "analysis.pdb").write_text(_pdb_acd())
    (run / "input.json").write_text(json.dumps({"tool": tool, "workflow": workflow,
                                                 "params": {"target_fasta": ">target\nEFG"}}))

    assert ingest_run(run, database=database, force=True) == 0
    assert database.groups(A, collection="all") == []


def test_generated_complex_refolding_is_library_output(database):
    run = _make_run(database, "design", "boltzgen-final", workflow="design_generation", tool="boltzgen",
                    candidates=[{"candidate_id": "design-1", "stage": "complex_refolding", "binder_sequence": "ACD",
                                 "binder_chains": ["A"], "complex_pdb": "complex.pdb",
                                 "metrics": {"confidence": 0.91, "native_pass_filters": True}}])
    (run / "complex.pdb").write_text(_pdb_acd())
    assert ingest_run(run, database=database) == 1
    library = database.groups(A, project_id="project-1", collection="library")
    assert len(library) == 1
    assert library[0]["origin"]["collection_kind"] == "generated"
    assert database.list_molecules(A, group_id=library[0]["id"], collection="library")["total"] == 1


def test_explicit_run_membership_removal_survives_forced_reindex(database):
    run = _make_run(database, "design", "removed-rfd4", tool="rfd4_proteina",
        workflow="design_generation", candidates=[{
            "candidate_id": "saved-1", "stage": "generation.backbone_sequence", "binder_sequence": "ACD",
            "binder_chains": ["A"], "binder_pdb": "binder.pdb",
            "metrics": {"native_pass_filters": True}}])
    (run / "binder.pdb").write_text(_pdb_acd())
    assert ingest_run(run, database=database) == 1
    group = next(g for g in database.groups(A, collection="all") if g["origin"].get("run_id") == run.name)
    molecule = database.list_molecules(A, group_id=group["id"])["items"][0]["id"]
    other = database.create_group(A, "Pending elsewhere")
    database.add_to_group(A, other, [molecule])
    database.remove_from_group(A, group["id"], [molecule])
    ingest_run(run, database=database, force=True)
    assert database.list_molecules(A, group_id=group["id"])["total"] == 0
    assert database.list_molecules(A, group_id=other)["total"] == 1


@pytest.mark.parametrize("status", ["rejected", "unresolved", "accepted"])
def test_saved_acceptance_decisions_are_not_confused_with_candidate_execution_status(database, status):
    source, molecule_id, source_group, source_structure, digest = _database_source(database, "ACDE")
    run = _make_run(database, "analysis", "acceptance-" + status, tool="validation_acceptance",
        workflow="analysis.validation_acceptance", candidates=[{
            "candidate_id": "source-1", "source_tool": "rfd4_proteina", "stage": "generation.backbone_sequence",
            "binder_sequence": "ACDE", "complex_pdb": "evaluated.pdb",
            "raw_metadata": {"database_source": source,
                             "complex_analysis_selection": {"source_structure_digest": digest},
                             "validation_acceptance": {"status": status}},
            "metrics": {"native_score": .8}}])
    (run / "evaluated.pdb").write_text(_pdb_for_sequence("ACDE"))
    assert ingest_run(run, database=database) == 1
    evaluation = next(group for group in database.groups(A, collection="all")
                      if group["origin"].get("collection_kind") == "evaluation")
    row = database.list_molecules(A, group_id=evaluation["id"])["items"][0]
    assert row["id"] == molecule_id
    assert row["scores"]["validation_acceptance.native_score"] == .8
    observation = row["observations"][0]
    assert observation["context"]["source_group_id"] == source_group
    assert observation["context"]["source_structure_id"] == source_structure
    assert observation["context"]["source_structure_digest"] == digest
    assert observation["context"]["acceptance_decision"]["status"] == status


@pytest.mark.parametrize(
    "task_group,workflow,input_shape",
    [
        ("analysis", "analysis.plip", "top_level"),
        ("pipeline", "pipeline.analysis.plip", "nested_inputs"),
    ],
)
def test_results_and_pipeline_evaluations_attach_to_exact_source_structure(
        database, task_group, workflow, input_shape):
    source = _make_run(database, "design", "results-source", tool="rfd4_proteina",
        workflow="design_generation", candidates=[{
            "candidate_id": "accepted-source", "stage": "generation.backbone_sequence",
            "binder_sequence": "ACD", "binder_chains": ["B"], "target_chains": ["A"],
            "complex_pdb": "source-complex.pdb",
            "raw_metadata": {"validation_acceptance": {"status": "accepted"}},
        }])
    (source / "source-complex.pdb").write_text(_pdb_efg_acd_complex())
    assert ingest_run(source, database=database) == 1
    source_group = next(group for group in database.groups(A, collection="all")
                        if group["origin"].get("run_id") == source.name)
    source_rows = database.list_molecules(A, group_id=source_group["id"], collection="all")["items"]
    molecule = next(row for row in source_rows
                    if any(component.get("sequence") == "ACD" for component in row["components"]))
    source_structure = molecule["representations"][0]
    source_path = database.structure_file(A, source_structure["id"])
    source_digest = hashlib.sha256(Path(source_path).read_bytes()).hexdigest()

    analysis = _make_run(database, task_group, f"{task_group}-analysis", tool="plip",
        workflow=workflow, candidates=[{
            "candidate_id": "plip-evaluation", "stage": "analysis", "source_tool": "plip",
            "binder_sequence": "ACD", "binder_chains": ["B"], "target_chains": ["A"],
            "complex_pdb": "analysis-complex.pdb",
            "metrics": {"plip_interaction_count": 8},
            "raw_metadata": {
                "source_candidate": {"candidate_id": "accepted-source"},
                "complex_analysis_selection": {"source_structure_digest": source_digest},
            },
        }])
    (analysis / "analysis-complex.pdb").write_text(_pdb_efg_acd_complex())
    input_record = {"tool": "plip", "workflow": workflow, "params": {}}
    if input_shape == "top_level":
        input_record["source_run_dir"] = str(source)
    else:
        input_record["inputs"] = {"source_run_dir": str(source)}
    (analysis / "input.json").write_text(json.dumps(input_record))

    assert ingest_run(analysis, database=database) == 1
    evaluation_group = next(group for group in database.groups(A, collection="all")
                            if group["origin"].get("collection_kind") == "evaluation")
    evaluated = database.list_molecules(A, group_id=evaluation_group["id"], collection="all")["items"]
    assert [row["id"] for row in evaluated] == [molecule["id"]]
    observation = next(row for row in evaluated[0]["observations"] if row["namespace"] == "plip")
    assert observation["context"]["source_kind"] == "results"
    assert observation["context"]["source_group_id"] == source_group["id"]
    assert observation["context"]["source_structure_id"] == source_structure["id"]
    assert observation["context"]["source_structure_digest"] == source_digest
    assert observation["context"]["representation_id"] == source_structure["id"]
    original_after = database.list_molecules(A, group_id=source_group["id"], collection="all")["items"][0]
    assert original_after["id"] == molecule["id"]
    assert len(original_after["representations"]) == 1


def test_sequence_analysis_resolves_exact_sequence_identity_without_coordinates(database):
    molecule = database.import_records(
        A, molecule_imports.parse_sequences("ACDE"), name="Sequence-only input",
        project_id="project-1", origin={"collection_kind": "curated"},
    )
    molecule_id, group_id = molecule["molecule_ids"][0], molecule["group_id"]
    source = _make_run(database, "design", "sequence-source-run", tool="rfd4_proteina",
                       workflow="design.rfd4_proteina", candidates=())
    with database.connection() as connection:
        connection.execute(
            "INSERT INTO molecule_sources VALUES (?,?,?,?)",
            ("exact-sequence-source", molecule_id, group_id, json.dumps({
                "run_id": source.name, "task_group": "design", "candidate_id": "sequence-source",
            })),
        )
    analysis = _make_run(database, "analysis", "sequence-qc-result", tool="sequence_qc",
        workflow="analysis.sequence_qc", candidates=[{
            "candidate_id": "sequence-qc-result", "source_tool": "sequence_qc", "stage": "analysis",
            "binder_sequence": "ACDE", "binder_chains": ["A"], "metrics": {"sequence_valid": True},
            "raw_metadata": {"source_candidate": {
                "candidate_id": "sequence-source", "binder_sequence": "ACDE", "binder_chains": ["A"],
            }},
        }])
    (analysis / "input.json").write_text(json.dumps({
        "tool": "sequence_qc", "workflow": "analysis.sequence_qc", "source_run_dir": str(source),
    }))

    assert ingest_run(analysis, database=database, force=True) == 1
    assert database.list_molecules(A, collection="all")["total"] == 1
    evaluation = next(group for group in database.groups(A, collection="all")
                      if group["origin"].get("collection_kind") == "evaluation")
    row = database.list_molecules(A, group_id=evaluation["id"], collection="all")["items"][0]
    assert row["id"] == molecule_id
    assert row["representations"] == []
    observation = next(item for item in row["observations"] if item["namespace"] == "sequence_qc")
    assert observation["context"]["source_reference_kind"] == "sequence"
    assert observation["context"]["source_sequence_sha256"] == hashlib.sha256(b"ACDE").hexdigest()
    assert observation["context"]["source_structure_id"] == ""
    assert observation["context"]["representation_id"] == ""


def test_database_sequence_assessment_preserves_recorded_design_target_without_target_chain_selection(database):
    from mn_molecular_design_studio.molecular_database_ingest import _database_analysis_subject

    target = database.import_records(A, molecule_imports.parse_sequences("EFG"),
                                     name="Recorded target", project_id="project-1")
    target_id = target["molecule_ids"][0]
    content = _pdb_efg_acd_complex().encode()
    components = molecule_imports.parse_structure(content, "complex.pdb")[0]["components"]
    imported = database.import_records(
        A, [{"name": "Accepted binder complex", "components": components,
             "files": [{"name": "complex.pdb", "content": content}]}],
        name="Accepted binder", project_id="project-1",
        origin={"type": "user", "design_target_molecule_id": target_id},
    )
    molecule_id, group_id = imported["molecule_ids"][0], imported["group_id"]
    representation = database.list_molecules(
        A, molecule_ids=[molecule_id], collection="all", limit=1,
    )["items"][0]["representations"][0]
    path = database.structure_file(A, representation["id"], materialize=True)
    digest = hashlib.sha256(Path(path).read_bytes()).hexdigest()

    subject = _database_analysis_subject(database, A, A["id"], "project-1", {
        "schema": "studio.complex_analysis_database_source.v1",
        "molecule_id": molecule_id, "group_id": group_id,
        "structure_id": representation["id"], "structure_digest": digest,
        "target_chains": [],
    })

    assert subject is not None
    assert subject["target_molecule_id"] == target_id
    assert subject["structure_id"] == representation["id"]
    assert subject["structure_digest"] == digest


def test_accepted_binder_waits_for_structure_then_imports_it_on_next_sync(database):
    run = _make_run(database, "design", "accepted-rfd4-late-structure", tool="rfd4_proteina",
        workflow="design_generation", candidates=[{
            "candidate_id": "accepted-rfd4-1", "stage": "generation.backbone_sequence",
            "binder_sequence": "ACD", "binder_chains": ["B"], "target_chains": ["A"],
            "complex_pdb": "accepted-complex.pdb",
            "raw_metadata": {"validation_acceptance": {"status": "accepted"}},
        }])

    # Do not promote a sequence-only accepted design to the reusable library.
    assert ingest_run(run, database=database) == 0
    assert database.groups(A, collection="library") == []

    # A later sync sees the changed structure signature and imports the exact
    # target-binder complex, without rerunning acceptance or prediction.
    (run / "accepted-complex.pdb").write_text(_pdb_efg_acd_complex())
    assert ingest_run(run, database=database) == 1
    group = next(group for group in database.groups(A, collection="library")
                 if group["origin"].get("run_id") == run.name)
    rows = database.list_molecules(A, group_id=group["id"], collection="library")["items"]
    binder = next(row for row in rows if any(component.get("sequence") == "ACD" for component in row["components"]))
    assert len(binder["representations"]) == 1
    structure_path = database.structure_file(A, binder["representations"][0]["id"])
    components = molecule_imports.parse_structure(Path(structure_path).read_bytes(), "accepted-complex.pdb")[0]["components"]
    assert {(component.get("chain"), component.get("sequence")) for component in components} == {
        ("A", "EFG"), ("B", "ACD"),
    }
    assert binder["components"][0]["role"] == "binder"


def test_secondary_candidate_structures_must_match_declared_chains_and_sequences(database):
    complex_text = _pdb_efg_acd_complex()
    wrong_binder = _pdb_for_sequence("KLM", chain="B")
    wrong_validation = _pdb_for_two_sequences("EFG", "KLM")
    candidate = {
        "candidate_id": "structure-correspondence", "stage": "generation.backbone_sequence",
        "binder_sequence": "ACD", "binder_chains": ["B"], "target_chains": ["A"],
        "complex_pdb": "complex.pdb", "binder_pdb": "binder.pdb",
        "raw_metadata": {"components": [
                             {"type": "protein", "chain": "A", "role": "target", "sequence": "EFG"},
                             {"type": "protein", "chain": "B", "role": "binder", "sequence": "ACD"},
                         ],
                         "validation_acceptance": {"status": "accepted"},
                         "validation_complex_pdb": "validation.pdb"},
        "metrics": {"native_pass_filters": True},
    }
    run = _make_run(database, "design", "mismatched-secondary-structures", tool="rfd4_proteina",
                    workflow="design_generation", candidates=[candidate])
    (run / "complex.pdb").write_text(complex_text)
    (run / "binder.pdb").write_text(wrong_binder)
    (run / "validation.pdb").write_text(wrong_validation)

    assert ingest_run(run, database=database, force=True) == 1
    group = next(row for row in database.groups(A, collection="library")
                 if row["origin"].get("run_id") == run.name)
    molecule = database.list_molecules(A, group_id=group["id"]) ["items"][0]
    fields = {row["context"].get("source_file_field") for row in molecule["representations"]}
    assert fields == {"complex_pdb"}
    assert database.structure_file(A, molecule["representations"][0]["id"]).read_bytes() == complex_text.encode()
    rejected = molecule["metadata"]["structure_validation"]
    assert {row["field"] for row in rejected} == {"binder_pdb", "validation_complex_pdb"}
    assert all(row["status"] == "molecular_mismatch" for row in rejected)


def test_candidate_structure_model_selection_is_applied_before_attachment(database):
    model_lines = []
    for model_number, binder in ((1, "KLM"), (2, "ACD")):
        atoms = [line for line in _pdb_for_two_sequences("EFG", binder).splitlines()
                 if line.startswith("ATOM") or line == "TER"]
        model_lines.extend([f"MODEL     {model_number:4d}", *atoms, "ENDMDL"])
    model_text = "\n".join([*model_lines, "END", ""])
    components = [
        {"type": "protein", "chain": "A", "role": "target", "sequence": "EFG"},
        {"type": "protein", "chain": "B", "role": "binder", "sequence": "ACD"},
    ]
    candidate = {
        "candidate_id": "selected-model", "stage": "generation.backbone_sequence",
        "binder_sequence": "ACD", "binder_chains": ["B"], "target_chains": ["A"],
        "complex_pdb": "complex.pdb", "metrics": {"native_pass_filters": True},
        "raw_metadata": {"components": components,
                         "complex_analysis_selection": {"model_id": "1"},
                         "validation_acceptance": {"status": "accepted"}},
    }
    run = _make_run(database, "design", "selected-model-attachment", tool="rfd4_proteina",
                    workflow="design_generation", candidates=[candidate])
    (run / "complex.pdb").write_text(model_text)

    assert ingest_run(run, database=database, force=True) == 1
    group = next(item for item in database.groups(A, collection="library")
                 if item["origin"].get("run_id") == run.name)
    molecule = database.list_molecules(A, group_id=group["id"], collection="library")["items"][0]
    assert [(item["chain"], item["sequence"]) for item in molecule["components"]] == [("B", "ACD")]
    structure = molecule["representations"][0]
    assert structure["context"]["model_id"] == "1"
    parsed = molecule_imports.parse_structure(
        database.structure_file(A, structure["id"]).read_bytes(), "selected.pdb", model_id="1",
        require_model_selection=True,
    )[0]
    assert {(item["chain"], item["sequence"]) for item in parsed["components"]} == {
        ("A", "EFG"), ("B", "ACD"),
    }

    unresolved = {**candidate, "candidate_id": "unselected-model",
                  "raw_metadata": {"components": components,
                                   "validation_acceptance": {"status": "accepted"}}}
    unresolved_run = _make_run(database, "design", "unselected-model-attachment", tool="rfd4_proteina",
                               workflow="design_generation", candidates=[unresolved])
    (unresolved_run / "complex.pdb").write_text(model_text)
    assert ingest_run(unresolved_run, database=database, force=True) == 0


def test_changed_coordinates_at_a_portable_path_create_a_new_measurement_revision(database):
    candidate = {
        "candidate_id": "same-path-revision", "stage": "generation.backbone_sequence",
        "binder_sequence": "ACD", "binder_chains": ["A"], "binder_pdb": "binder.pdb",
        "raw_metadata": {"validation_acceptance": {"status": "accepted"}},
        "metrics": {"native_score": 0.8},
    }
    run = _make_run(database, "design", "coordinate-revision", tool="rfd4_proteina",
                    workflow="design_generation", candidates=[candidate])
    source = run / "binder.pdb"
    source.write_text(_pdb_acd())
    assert ingest_run(run, database=database, force=True) == 1
    group = next(row for row in database.groups(A, collection="library")
                 if row["origin"].get("run_id") == run.name)
    molecule = database.list_molecules(A, group_id=group["id"])["items"][0]
    first_structure = molecule["representations"][0]
    first_digest = first_structure["context"]["portable_sha256"]

    shifted = []
    for line in _pdb_acd().splitlines():
        if line.startswith("ATOM"):
            x = float(line[30:38]) + 17.0
            line = line[:30] + f"{x:8.3f}" + line[38:]
        shifted.append(line)
    second_text = "\n".join(shifted) + "\n"
    source.write_text(second_text)
    candidate["metrics"] = {"native_score": 0.4}
    (run / "artifacts/normalized_candidates/candidates.jsonl").write_text(json.dumps(candidate) + "\n")
    assert ingest_run(run, database=database, force=True) == 1

    refreshed = database.list_molecules(A, group_id=group["id"])["items"][0]
    assert len(refreshed["representations"]) == 2
    observation = next(row for row in refreshed["observations"] if row["namespace"] == "rfd4_proteina")
    assert observation["metrics"]["native_score"] == 0.4
    assert observation["context"]["representation_id"] != first_structure["id"]
    assert observation["context"]["structure_digest"] != first_digest
    assert observation["context"]["structure_digest"] == hashlib.sha256(second_text.encode()).hexdigest()


def test_changed_calculation_evidence_preserves_frozen_evaluation_and_is_idempotent(database):
    candidate = {
        "candidate_id": "frozen-score", "stage": "generation.backbone_sequence",
        "binder_sequence": "ACD", "binder_chains": ["A"], "binder_pdb": "binder.pdb",
        "raw_metadata": {"validation_acceptance": {"status": "accepted"}}, "metrics": {"native_score": 0.8},
    }
    run = _make_run(database, "design", "frozen-calculation-revision", tool="rfd4_proteina",
                    workflow="design_generation", candidates=[candidate])
    source = run / "binder.pdb"
    source.write_text(_pdb_acd())
    assert ingest_run(run, database=database, force=True) == 1
    group = next(item for item in database.groups(A, collection="library")
                 if item["origin"].get("run_id") == run.name)
    molecule_id = database.list_molecules(A, group_id=group["id"])["items"][0]["id"]
    frozen = database.save_comparison(A, "Frozen 0.8", [molecule_id],
        settings={"fields": ["rfd4_proteina.native_score"],
                  "directions": {"rfd4_proteina.native_score": "higher"}})
    original = database.list_molecules(A, group_id=group["id"])["items"][0]
    original_observation = next(item for item in original["observations"] if item["namespace"] == "rfd4_proteina")
    original_evaluation_id = original["selected_evaluation"]["id"]
    original_structure_id = original_observation["context"]["representation_id"]
    shifted = []
    for line in _pdb_acd().splitlines():
        if line.startswith("ATOM"):
            x = float(line[30:38]) + 17.0
            line = line[:30] + f"{x:8.3f}" + line[38:]
        shifted.append(line)
    changed_text = "\n".join(shifted) + "\n"
    source.write_text(changed_text)
    candidate["metrics"] = {"native_score": 0.4}
    (run / "artifacts/normalized_candidates/candidates.jsonl").write_text(json.dumps(candidate) + "\n")
    assert ingest_run(run, database=database, force=True) == 1
    assert ingest_run(run, database=database, force=True) == 1

    current = database.list_molecules(A, group_id=group["id"])["items"][0]
    evaluations = [item for item in current["observations"] if item["namespace"] == "rfd4_proteina"]
    assert len(evaluations) == 2
    by_score = {item["metrics"]["native_score"]: item for item in evaluations}
    assert by_score[0.8]["context"]["representation_id"] == original_structure_id
    assert by_score[0.4]["context"]["structure_digest"] == hashlib.sha256(changed_text.encode()).hexdigest()
    assert len({item["id"] for item in evaluations}) == 2
    frozen_row = database.list_molecules(A, group_id=frozen["group_id"])["items"][0]
    assert frozen_row["selected_evaluation"]["id"] == original_evaluation_id
    assert frozen_row["scores"]["rfd4_proteina.native_score"] == 0.8
    assert frozen_row["selected_evaluation"]["representation_id"] == original_structure_id
    assert database.structure_file(A, original_structure_id).read_bytes() == _pdb_acd().encode()


@pytest.mark.parametrize("legacy_model_id", ["0", None])
def test_legacy_structure_reference_is_reused_for_identical_reingestion(database, legacy_model_id):
    from mn_molecular_design_studio.molecular_database import _id, _json
    from mn_molecular_design_studio.molecular_database_ingest import _results_analysis_subject

    imported = database.import_records(A, molecule_imports.parse_sequences("ACD"), name="Legacy source",
                                       project_id="project-1")
    group, molecule = imported["group_id"], imported["molecule_ids"][0]
    path = database.root.parent / "legacy-source.pdb"
    path.write_text(_pdb_acd())
    context = {"source_file_field": "complex_pdb"}
    if legacy_model_id is not None:
        context["model_id"] = legacy_model_id
    current_id = database.attach_structure(group, molecule,
        {"path": str(path), "label": "legacy-source.pdb"}, context=context)
    with database.connection() as connection:
        row = connection.execute("SELECT * FROM molecule_structures WHERE id=?", (current_id,)).fetchone()
        legacy_context = json.loads(row["context"])
        legacy_context.pop("content_sha256", None)
        legacy_id = _id("structure-", _json([molecule, group, str(path), legacy_context]))
        connection.execute("UPDATE molecule_structures SET id=? WHERE id=?", (legacy_id, current_id))
        connection.execute("INSERT INTO molecule_sources VALUES (?,?,?,?)", (
            "legacy-source-provenance", molecule, group,
            json.dumps({"run_id": "legacy-run", "task_group": "design", "candidate_id": "legacy-candidate",
                        "workstation": "legacy-owner"}),
        ))
    source_run = database.root.parent / "runs" / "design" / "legacy-run"
    source_run.mkdir(parents=True)
    (source_run / "metadata.json").write_text(json.dumps({"run_id": "legacy-run", "task_group": "design"}))
    (source_run / "molecular_design_studio.json").write_text(json.dumps({"workstation": "legacy-owner"}))
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    analysis_candidate = {"binder_sequence": "ACD", "binder_chains": ["A"],
        "raw_metadata": {"source_candidate": {"candidate_id": "legacy-candidate", "binder_sequence": "ACD",
                                                 "binder_chains": ["A"]},
                         "complex_analysis_selection": {"source_structure_digest": digest}}}
    source = _results_analysis_subject(database, A, A["id"], "project-1", source_run,
        {"source_run_dir": str(source_run)}, analysis_candidate)
    assert source["structure_id"] == legacy_id

    repeated = database.attach_structure(group, molecule,
        {"path": str(path), "label": "legacy-source.pdb"}, context={**context, "model_id": "0"})
    assert repeated == legacy_id
    repeated_source = _results_analysis_subject(database, A, A["id"], "project-1", source_run,
        {"source_run_dir": str(source_run)}, analysis_candidate)
    assert repeated_source["structure_id"] == legacy_id
    with database.connection() as connection:
        rows = connection.execute(
            "SELECT id FROM molecule_structures WHERE group_id=? AND molecule_id=?", (group, molecule)
        ).fetchall()
    assert [row["id"] for row in rows] == [legacy_id]
    revised = database.attach_structure(group, molecule,
        {"path": str(path), "label": "legacy-source.pdb"}, context={**context, "model_id": "1"})
    assert revised != legacy_id


def test_results_analysis_rejects_calculation_sequence_mismatch_with_exact_source(database):
    source = _make_run(database, "design", "acd-source-for-analysis", tool="rfd4_proteina",
        workflow="design_generation", candidates=[{
            "candidate_id": "selected-acd", "stage": "generation.backbone_sequence",
            "binder_sequence": "ACD", "binder_chains": ["A"], "binder_pdb": "source.pdb",
            "raw_metadata": {"validation_acceptance": {"status": "accepted"}},
        }])
    (source / "source.pdb").write_text(_pdb_for_sequence("ACD"))
    assert ingest_run(source, database=database, force=True) == 1
    source_row = next(item for item in database.list_molecules(A, collection="library")["items"]
                      if item["components"][0]["sequence"] == "ACD")
    source_digest = hashlib.sha256(Path(database.structure_file(A, source_row["representations"][0]["id"])).read_bytes()).hexdigest()
    analysis = _make_run(database, "analysis", "wrong-sequence-analysis", tool="plip",
        workflow="analysis.plip", candidates=[{
            "candidate_id": "plip-wrong-molecule", "stage": "analysis", "source_tool": "plip",
            "binder_sequence": "KLM", "binder_chains": ["A"], "complex_pdb": "calculation.pdb",
            "metrics": {"interaction_count": 3},
            "raw_metadata": {"source_candidate": {"candidate_id": "selected-acd", "binder_sequence": "ACD",
                                                       "binder_chains": ["A"]},
                             "complex_analysis_selection": {"source_structure_digest": source_digest}},
        }])
    (analysis / "input.json").write_text(json.dumps({"tool": "plip", "workflow": "analysis.plip",
                                                       "source_run_dir": str(source)}))
    (analysis / "calculation.pdb").write_text(_pdb_for_sequence("KLM"))
    before = database.list_molecules(A, collection="all")["total"]

    assert ingest_run(analysis, database=database, force=True) == 0
    assert database.list_molecules(A, collection="all")["total"] == before
    with database.connection() as connection:
        assert connection.execute("SELECT count(*) FROM molecule_observations WHERE run_id=?", (analysis.name,)).fetchone()[0] == 0
    assert (analysis / "calculation.pdb").read_text() == _pdb_for_sequence("KLM")
    assert (analysis / "artifacts/normalized_candidates/candidates.jsonl").is_file()


def test_generated_complex_reindex_keeps_binder_identity_and_only_passing_designs(database):
    candidate = {
        "candidate_id": "bindcraft2_00001", "stage": "complex_refolding",
        "binder_sequence": "EFG", "binder_chains": ["B"], "target_chains": ["A"],
        "complex_pdb": "complex.pdb",
        "metrics": {"native_pass_filters": True},
    }
    run = _make_run(database, "design", "bindcraft2-reindex", workflow="design_generation",
                    tool="bindcraft2", candidates=[candidate])
    (run / "complex.pdb").write_text(_pdb_acd_efg_complex())

    assert ingest_run(run, database=database) == 1
    group = next(group for group in database.groups(A, collection="all") if group["origin"].get("run_id") == "bindcraft2-reindex")
    original_id = database.list_molecules(A, group_id=group["id"])["items"][0]["id"]

    # Recreate the earlier duplicate full-complex row produced when structure
    # parsing changed a candidate's identity from binder-only to target+binder.
    complex_text = _pdb_acd_efg_complex()
    complex_path = run / "complex.pdb"
    complex_path.write_text(complex_text)
    full_components = molecule_imports.parse_structure(complex_text.encode(), "complex.pdb")[0]["components"]
    with database.connection() as connection:
        duplicate_id = database._put(connection, group["id"], {
            "name": "bindcraft2_00001", "components": full_components, "files": [],
        })
    database.attach_structure(group["id"], duplicate_id, {"path": str(complex_path), "label": "Complex Pdb"})
    assert duplicate_id != original_id

    passing = {**candidate, "complex_pdb": "complex.pdb"}
    rejected = {
        "candidate_id": "bindcraft2_00002", "stage": "complex_refolding",
        "binder_sequence": "VVV", "binder_chains": ["B"],
        "metrics": {"native_pass_filters": False},
    }
    (run / "artifacts/normalized_candidates/candidates.jsonl").write_text(
        "".join(json.dumps(row) + "\n" for row in (passing, rejected))
    )

    assert ingest_run(run, database=database, force=True) == 1
    rows = database.list_molecules(A, group_id=group["id"])["items"]
    assert len(rows) == 2  # preserve the already-indexed full complex identity
    binder = next(row for row in rows if row["id"] == original_id)
    assert [(component["chain"], component["sequence"]) for component in binder["components"]] == [("B", "EFG")]
    assert len(binder["structures"]) == 1
    assert binder["structures"][0]["context"]["portable"] is True
    assert not any(any(component.get("sequence") == "VVV" for component in row["components"]) for row in rows)


def test_generated_binder_links_existing_exact_modeled_target_and_preserves_complex(database):
    import hashlib
    target = database.import_records(A, molecule_imports.parse_sequences("ACD"),
                                     name="Canonical target", project_id="project-1")
    target_id = target["molecule_ids"][0]
    complex_text = _pdb_acd().replace("TER\nEND", "TER\n" + "\n".join([
        "ATOM      4  CA  GLU B   1       8.000   1.000   1.000  1.00 20.00           C",
        "ATOM      5  CA  PHE B   2       9.000   1.000   1.000  1.00 20.00           C",
        "ATOM      6  CA  GLY B   3      10.000   1.000   1.000  1.00 20.00           C",
        "END",
    ]))
    passing = {"candidate_id":"rfd4-pass", "stage":"generation.backbone_sequence",
               "binder_sequence":"EFG", "binder_chains":["B"], "target_chains":["A"],
               "complex_pdb":"complex.pdb", "metrics":{"native_pass_filters":True}}
    rejected = {**passing, "candidate_id":"rfd4-reject", "binder_sequence":"VVV",
                "metrics":{"native_pass_filters":False}}
    run = _make_run(database, "design", "rfd4-target-link", workflow="design_generation",
                    tool="rfd4_proteina", candidates=[passing, rejected])
    (run/"complex.pdb").write_text(complex_text)
    for force in (False, True):
        assert ingest_run(run, database=database, force=force) == 1
        group = next(g for g in database.groups(A, collection="library")
                     if g["origin"].get("run_id") == run.name)
        rows = database.list_molecules(A, group_id=group["id"])["items"]
        assert len(rows) == 1
        binder = rows[0]
        assert [(c["chain"],c["sequence"]) for c in binder["components"]] == [("B","EFG")]
        assert binder["metadata"]["design_target_molecule_id"] == target_id
        assert len(binder["representations"]) == 1
        structure = binder["representations"][0]
        assert structure["context"]["design_target_molecule_id"] == target_id
        assert structure["context"]["design_target_sequence_sha256"] == hashlib.sha256(b"ACD").hexdigest()
        assert structure["context"]["design_target_is_structural_reference"] is False
        assert structure["context"]["binder_chains"] == ["B"]
        assert structure["context"]["target_chains"] == ["A"]
        assert database.structure_file(A, structure["id"]).read_bytes() == complex_text.encode()
        assert binder["observations"][0]["context"]["target_id"] == target_id
        assert database.list_molecules(A, molecule_ids=[target_id])["total"] == 1


def test_generated_binder_registers_explicit_sequence_only_target_and_links_identity(database):
    import hashlib

    candidate = {
        "candidate_id": "sequence-target-design", "stage": "generation.backbone_sequence",
        "binder_sequence": "ACD", "binder_chains": ["A"], "target_chains": ["B"],
        "binder_pdb": "binder.pdb", "metrics": {"native_pass_filters": True},
        "raw_metadata": {"sequences_by_chain": {"B": "EFGH"}},
    }
    run = _make_run(database, "design_generation", "sequence-target-link", workflow="design_generation",
                    tool="rfd4_proteina", candidates=[candidate])
    (run / "binder.pdb").write_text(_pdb_acd())

    assert ingest_run(run, database=database, force=True) == 1
    group = next(g for g in database.groups(A, collection="library")
                 if g["origin"].get("run_id") == run.name)
    binder = database.list_molecules(A, group_id=group["id"])["items"][0]
    target_id = binder["metadata"]["design_target_molecule_id"]
    target = database.list_molecules(A, molecule_ids=[target_id], collection="all")["items"][0]
    assert [(component["chain"], component["sequence"], component["role"])
            for component in target["components"]] == [("B", "EFGH", "target")]
    assert target["metadata"]["target_identity_type"] == "sequence_only"
    assert target["representations"] == []
    assert binder["metadata"]["design_target_is_structural_reference"] is False
    assert binder["metadata"]["design_target_sequence_sha256"] == hashlib.sha256(b"EFGH").hexdigest()


def test_curated_design_keeps_binder_complex_and_supplied_target_as_distinct_representations(database):
    import hashlib

    complex_text = _pdb_efg_acd_complex()
    binder_text = _pdb_for_sequence("ACD", chain="B")
    target_text = _pdb_for_sequence("EFG", chain="A")
    candidate = {
        "candidate_id": "design-with-three-structures", "stage": "generation.backbone_sequence",
        "binder_sequence": "ACD", "binder_chains": ["B"], "target_chains": ["A"],
        "binder_pdb": "binder.pdb", "complex_pdb": "complex.pdb", "target_pdb": "target.pdb",
        "raw_metadata": {"validation_acceptance": {"status": "accepted"}, "structure_kind": "predicted"},
    }
    run = _make_run(database, "design", "three-structure-design", workflow="design_generation",
                    tool="rfd4_proteina", candidates=[candidate])
    (run / "binder.pdb").write_text(binder_text)
    (run / "complex.pdb").write_text(complex_text)
    (run / "target.pdb").write_text(target_text)
    (run / "input.json").write_text(json.dumps({"tool": "rfd4_proteina", "workflow": "design_generation",
                                                "params": {"target_pdb": "target.pdb"}}))

    assert ingest_run(run, database=database) == 2
    group = next(row for row in database.groups(A, collection="library")
                 if row["origin"].get("run_id") == run.name)
    binder = next(row for row in database.list_molecules(A, group_id=group["id"], collection="library")["items"]
                  if row["components"][0].get("sequence") == "ACD")
    assert len(binder["representations"]) == 2
    by_field = {item["context"]["source_file_field"]: item for item in binder["representations"]}
    assert set(by_field) == {"complex_pdb", "binder_pdb"}
    assert by_field["complex_pdb"]["context"]["representation_type"] == "predicted"
    assert by_field["binder_pdb"]["context"]["representation_type"] == "predicted"
    assert database.structure_file(A, by_field["complex_pdb"]["id"]).read_bytes() == complex_text.encode()
    assert database.structure_file(A, by_field["binder_pdb"]["id"]).read_bytes() == binder_text.encode()
    assert by_field["complex_pdb"]["context"]["portable_sha256"] == hashlib.sha256(complex_text.encode()).hexdigest()
    assert by_field["binder_pdb"]["context"]["portable_sha256"] == hashlib.sha256(binder_text.encode()).hexdigest()

    target_id = binder["metadata"]["design_target_molecule_id"]
    target = database.list_molecules(A, molecule_ids=[target_id], collection="all")["items"][0]
    assert len(target["representations"]) == 1
    target_structure = target["representations"][0]
    assert target_structure["context"]["input_role"] == "receptor"
    assert database.structure_file(A, target_structure["id"]).read_bytes() == target_text.encode()
    assert target_structure["context"]["portable_sha256"] == hashlib.sha256(target_text.encode()).hexdigest()


def test_sequence_only_run_target_is_ingested_without_candidate_structures(database):
    run = _make_run(database, "design_generation", "target-only-run", workflow="design_generation",
                    tool="rfd4_proteina", candidates=())
    (run / "input.json").write_text(json.dumps({
        "tool": "rfd4_proteina", "workflow": "design_generation",
        "params": {"target_fasta": ">IDR target\nACDEFG"},
    }))

    assert ingest_run(run, database=database, force=True) == 1
    target = next(item for item in database.list_molecules(A, collection="all")["items"]
                  if item.get("metadata", {}).get("target_identity_type") == "sequence_only")
    assert target["components"][0]["sequence"] == "ACDEFG"
    assert target["representations"] == []


def test_mismatched_input_target_coordinates_do_not_attach_to_sequence_identity(database):
    run = _make_run(database, "design_generation", "mismatched-input-target", tool="rfd4_proteina",
                    workflow="design.rfd4_proteina", candidates=())
    (run / "target.pdb").write_text(_pdb_for_sequence("KLMN"))
    (run / "input.json").write_text(json.dumps({
        "tool": "rfd4_proteina", "workflow": "design.rfd4_proteina",
        "params": {"target_pdb": "target.pdb", "target_sequence": "ACDE", "target_chains": ["A"]},
    }))

    assert ingest_run(run, database=database, force=True) == 1
    rows = database.list_molecules(A, project_id="project-1", collection="all")["items"]
    target = next(row for row in rows if row.get("metadata", {}).get("target_identity_type") == "sequence_only")
    assert target["components"][0]["sequence"] == "ACDE"
    assert target["representations"] == []


def test_rfd4_idr_sequence_carrier_is_indexed_as_sequence_not_structure(database):
    from mn_molecular_design_studio.rfd4_proteina_input import write_idr_sequence_carrier
    from mn_molecular_design_studio.target_identity import is_sequence_carrier_reference, target_reference_type

    run = _make_run(database, "design_generation", "idr-carrier-input", tool="rfd4_proteina",
                    workflow="design_generation")
    carrier = run / "artifacts/rfd4_proteina/input/idr_sequence_carrier.cif"
    write_idr_sequence_carrier("ACDE", carrier)
    (run / "input.json").write_text(json.dumps({
        "tool": "rfd4_proteina", "workflow": "design_generation",
        "inputs": {"target_structure": str(carrier), "input_files": [str(carrier)]},
        "params": {"task": "protein_binder", "input_mode": "condition_spec"},
    }))
    old_target = database.import_records(A, [{
        "name": "Previously misclassified IDR carrier",
        "components": [{"type": "protein", "chain": "A", "role": "target", "sequence": "ACDE"}],
        "metadata": {"target_identity_type": "structure_backed", "structural_reference": True},
    }], name="Old target import", project_id="project-1", import_method="tool_output")
    database.attach_structure(old_target["group_id"], old_target["molecule_ids"][0],
                              {"path": str(carrier), "label": "Input receptor"},
                              context={"run_id": run.name, "input_role": "receptor"})

    assert is_sequence_carrier_reference(str(carrier), run_dir=run)
    assert target_reference_type(json.loads((run / "input.json").read_text()), run_dir=run) == "sequence_only"
    assert ingest_run(run, database=database) == 1
    records = database.list_molecules(A, project_id="project-1", collection="all")["items"]
    target = next(row for row in records if [component.get("sequence") for component in row["components"]] == ["ACDE"])
    assert target["metadata"]["target_identity_type"] == "sequence_only"
    assert target["metadata"]["structural_reference"] is False
    assert target["representations"] == []


def test_rfd4_sequence_input_carrier_is_never_registered_as_native_reference(database):
    from mn_molecular_design_studio.engine.workflow_engine.workflows.rfd4_proteina import sequence_structure
    from mn_molecular_design_studio.target_identity import is_sequence_carrier_reference, target_reference_type

    run = _make_run(database, "design_generation", "sequence-input-carrier", tool="rfd4_proteina",
                    workflow="design.rfd4_proteina")
    carrier = sequence_structure("ACDE", run / "artifacts/sequence_input.cif")
    (run / "input.json").write_text(json.dumps({
        "tool": "rfd4_proteina", "workflow": "design.rfd4_proteina",
        "params": {"target_structure": str(carrier), "target_sequence": "ACDE", "target_chains": ["A"]},
    }))

    assert is_sequence_carrier_reference(str(carrier), run_dir=run)
    assert target_reference_type(json.loads((run / "input.json").read_text()), run_dir=run) == "sequence_only"
    assert ingest_run(run, database=database, force=True) == 1
    target = next(row for row in database.list_molecules(A, project_id="project-1", collection="all")["items"]
                  if row.get("metadata", {}).get("target_identity_type") == "sequence_only")
    assert target["components"][0]["sequence"] == "ACDE"
    assert target["representations"] == []


@pytest.mark.parametrize("producer", ["idr", "generic"])
def test_rfd4_artificial_carrier_remains_nonstructural_after_copy_to_target_cif(database, producer):
    from mn_molecular_design_studio.rfd4_proteina_input import write_idr_sequence_carrier
    from mn_molecular_design_studio.engine.workflow_engine.workflows.rfd4_proteina import sequence_structure
    from mn_molecular_design_studio.target_identity import is_sequence_carrier_reference, target_reference_type

    run = _make_run(database, "design_generation", f"renamed-{producer}-carrier", tool="rfd4_proteina",
                    workflow="design_generation")
    original = run / "artifacts" / "conditioning.cif"
    if producer == "idr":
        write_idr_sequence_carrier("ACDE", original)
    else:
        sequence_structure("ACDE", original)
    renamed = run / "target.cif"
    renamed.write_bytes(original.read_bytes())
    assert is_sequence_carrier_reference(str(renamed), run_dir=run)
    assert target_reference_type({"inputs": {"target_structure": "target.cif"}}, run_dir=run) == "sequence_only"
    ordinary = run / "ordinary-target.cif"
    ordinary.write_text("data_experimental_structure\n#\n")
    assert not is_sequence_carrier_reference(str(ordinary), run_dir=run)


def test_design_target_links_to_matching_structure_backed_target(database):
    candidate = {
        "candidate_id": "structure-target-design", "stage": "generation.backbone_sequence",
        "binder_sequence": "EFG", "binder_chains": ["B"], "target_chains": ["A"],
        "complex_pdb": "complex.pdb",
        "raw_metadata": {"components": [
            {"type": "protein", "chain": "A", "sequence": "ACD", "role": "target"},
            {"type": "protein", "chain": "B", "sequence": "EFG", "role": "binder"},
        ]},
        "metrics": {"native_pass_filters": True},
    }
    run = _make_run(database, "design_generation", "structure-target-link", tool="rfd4_proteina",
                    workflow="design_generation", candidates=[candidate])
    target_path = run / "target.pdb"
    target_path.write_text(_pdb_acd())
    (run / "complex.pdb").write_text(_pdb_acd_efg_complex())
    (run / "input.json").write_text(json.dumps({
        "tool": "rfd4_proteina", "workflow": "design_generation",
        "inputs": {"target_structure": "target.pdb"}, "params": {"target_chains": ["A"]},
    }))

    assert ingest_run(run, database=database) == 2
    rows = database.list_molecules(A, collection="all")["items"]
    binder = next(row for row in rows if row["name"] == "structure-target-design")
    target_id = binder["metadata"]["design_target_molecule_id"]
    assert binder["metadata"]["design_target_is_structural_reference"] is True
    target = database.list_molecules(A, molecule_ids=[target_id], collection="all")["items"][0]
    assert target["metadata"]["target_identity_type"] == "structure_backed"
    assert target["representations"]


def test_filtered_campaign_ingests_only_acceptance_passes(database):
    passing = {
        "candidate_id": "campaign-pass", "stage": "generation.backbone_sequence",
        "source_tool": "rfd4_proteina", "binder_sequence": "ACD", "binder_pdb": "pass.pdb",
        "binder_chains": ["A"],
        "metrics": {"native_pass_filters": True},
        "raw_metadata": {"distributed_campaign": {
            "validation_required": True, "validated_acceptance": "accepted",
        }},
    }
    rejected = {
        "candidate_id": "campaign-reject", "stage": "generation.backbone_sequence",
        "source_tool": "rfd4_proteina", "binder_sequence": "FGHI",
        "raw_metadata": {"distributed_campaign": {
            "validation_required": True, "validated_acceptance": "rejected",
        }},
    }
    pending = {
        "candidate_id": "campaign-pending", "stage": "generation.backbone_sequence",
        "source_tool": "rfd4_proteina", "binder_sequence": "KLMN",
        "raw_metadata": {"distributed_campaign": {
            "validation_required": True, "validated_acceptance": "pending",
        }},
    }
    run = _make_run(
        database, "design", "filtered-campaign", workflow="rfd4_proteina_campaign",
        tool="rfd4_proteina", candidates=[passing, rejected, pending],
    )
    (run / "pass.pdb").write_text(_pdb_acd())

    assert ingest_run(run, database=database) == 1
    group = next(group for group in database.groups(A, collection="library")
                 if group["origin"].get("run_id") == "filtered-campaign")
    items = database.list_molecules(A, group_id=group["id"], collection="library")["items"]
    assert [item["components"][0]["sequence"] for item in items] == ["ACD"]


def test_failed_native_refold_stays_in_results_without_database_group(database):
    complex_text = "\n".join([
        "ATOM      1  CA  ALA A   1       1.000   1.000   1.000  1.00 20.00           C",
        "ATOM      2  CA  CYS A   2       2.000   1.000   1.000  1.00 20.00           C",
        "ATOM      3  CA  ASP A   3       3.000   1.000   1.000  1.00 20.00           C",
        "ATOM      4  CA  GLU B   1       8.000   1.000   1.000  1.00 20.00           C",
        "ATOM      5  CA  PHE B   2       9.000   1.000   1.000  1.00 20.00           C",
        "ATOM      6  CA  GLY B   3      10.000   1.000   1.000  1.00 20.00           C",
        "END", "",
    ])
    candidate = {
        "candidate_id": "boltzgen_00001", "stage": "complex_refolding",
        "complex_pdb": "complex.pdb", "binder_pdb": "complex.pdb",
        "binder_sequence": "EFG", "binder_chains": ["B"], "target_chains": ["A"],
        "metrics": {"pass_filters": False, "refolded_motif_rmsd_angstrom": 3.2},
        "normalized_output": {"metrics": {"mimosa.refolded_motif_rmsd": {
            "value": 3.2, "unit": "Å", "scope": "refolded motif",
        }}},
        "raw_metadata": {"result_kind": "native_pipeline", "candidate_pool_level": "native_pipeline"},
    }
    run = _make_run(database, "design", "failed-native-refold", workflow="design_generation", tool="boltzgen",
                    candidates=[candidate])
    (run / "complex.pdb").write_text(complex_text)

    assert ingest_run(run, database=database) == 0
    assert database.groups(A, collection="library") == []
    assert not [group for group in database.groups(A, collection="all")
                if group["origin"].get("run_id") == "failed-native-refold"]


@pytest.mark.parametrize("complex_identity", [False, True])
def test_protein_qc_reuses_the_exact_frozen_database_identity(database, complex_identity):
    binder = {"type": "protein", "chain": "B", "role": "binder", "sequence": "ACD"}
    target = {"type": "protein", "chain": "A", "role": "target", "sequence": "ACD"}
    components = [target, binder] if complex_identity else [binder]
    imported = database.import_records(A, [{"name": "Source design", "components": components}], project_id="project-1")
    molecule = imported["molecule_ids"][0]
    run = _make_run(database, "analysis", "qc-identity", tool="sequence_qc", workflow="analysis.sequence_qc", candidates=[{
        "candidate_id": "qc-source", "stage": "analysis", "binder_sequence": "ACD", "binder_chains": ["B"], "target_chains": ["A"],
        "complex_pdb": "complex.pdb", "metrics": {"qc_status": "completed", "qc_sequence_length": 3},
        "raw_metadata": {"analysis_subject": "protein", "studio_molecule_id": molecule, "components": components},
    }])
    protein_lines = [line for line in _pdb_acd().splitlines() if line.startswith("ATOM")]
    (run / "complex.pdb").write_text("\n".join([*protein_lines, "TER", *[line[:21] + "B" + line[22:] for line in protein_lines], "TER", "END", ""]))
    assert ingest_run(run, database=database) == 1
    rows = database.list_molecules(A, collection="all")["items"]
    assert len(rows) == 1
    assert rows[0]["id"] == molecule
    assert rows[0]["scores"]["sequence_qc.qc_sequence_length"] == 3
    # The sequence check resolves to this exact molecular identity, but its
    # analysis artifact is not asserted to be the molecule's source structure.
    assert rows[0]["structures"] == []


def test_prediction_indexes_only_the_best_structure_per_input(database):
    run = _make_run(database, "optimized-prediction", "best-model", tool="optimized_prediction_boltz2",
                    workflow="refolding.optimized_prediction", candidates=[
        {"candidate_id": "model-a", "stage": "complex_refolding", "complex_pdb": "model-a.pdb",
         "raw_metadata": {"batch_item_id": "input-1", "compound_name": "Probe", "ligand_smiles": ["CCO"]},
         "metrics": {"ranking_score": .72}},
        {"candidate_id": "model-b", "stage": "complex_refolding", "complex_pdb": "model-b.pdb",
         "raw_metadata": {"batch_item_id": "input-1", "compound_name": "Probe", "ligand_smiles": ["CCO"]},
         "metrics": {"ranking_score": .91}},
    ])
    structure = _pdb_acd_with_generic_ligand()
    (run / "model-a.pdb").write_text(structure)
    (run / "model-b.pdb").write_text(structure)

    assert ingest_run(run, database=database) == 1
    row = database.list_molecules(A, collection="all")["items"][0]
    assert len(row["structures"]) == 1
    assert row["structures"][0]["context"]["candidate_id"] == "model-b"
    assert row["structures"][0]["context"]["database_structure_policy"] == "best_prediction"
    assert row["structures"][0]["context"]["portable"] is True
    assert Path(row["structures"][0]["path"]).is_relative_to(database.root / "molecular_files")
    assert len(row["observations"]) == 2

    # Reindexing after a changed ranking replaces the old attached representative.
    candidate_file = run / "artifacts/normalized_candidates/candidates.jsonl"
    candidates = [json.loads(line) for line in candidate_file.read_text().splitlines()]
    candidates[0]["metrics"]["ranking_score"] = .99
    candidates[1]["metrics"]["ranking_score"] = .4
    candidate_file.write_text("".join(json.dumps(candidate) + "\n" for candidate in candidates))
    assert ingest_run(run, database=database, force=True) == 1
    row = database.list_molecules(A, collection="all")["items"][0]
    assert len(row["structures"]) == 1
    assert row["structures"][0]["context"]["candidate_id"] == "model-a"


def test_docking_indexes_only_the_top_pose_for_each_compound(database):
    from mn_molecular_design_studio import ligand_compounds

    database.import_records(A, [{"name": "Ethanol", "components": molecule_imports.parse_smiles("CCO")[0]["components"],
                                 "metadata": {"compound_id": "CMP1"}}], name="Imported compounds", project_id="project-1")
    run = _make_run(database, "vina", "best-pose", tool="vina", workflow="ligand.docking")
    ligand_compounds.write_compound_set(
        run / "artifacts/input/compound_set.json", name="Inputs",
        compounds=[{"compound_id": "CMP1", "compound_name": "Ethanol", "smiles": "CCO"}],
        owner_id=A["id"], project_id="project-1", created_at="now",
    )
    poses = run / "artifacts/poses"
    poses.mkdir(parents=True)
    mapped_pose = "\n".join(["MODEL 1", "REMARK SMILES CCO", "REMARK SMILES IDX 1 3 2 2 3 1",
        "ATOM      1 O    UNL     1       8.000   3.000   4.000  0.000 O",
        "ATOM      2 C    UNL     1       7.000   3.000   4.000  0.000 C",
        "ATOM      3 C    UNL     1       6.000   3.000   4.000  0.000 C", "ENDMDL", ""])
    (poses / "pose-1.pdbqt").write_text(mapped_pose.replace("MODEL 1", "REMARK top-ranked pose\nMODEL 1"))
    (poses / "pose-2.pdbqt").write_text(mapped_pose.replace("MODEL 1", "REMARK lower-ranked pose\nMODEL 1"))
    (run / "artifacts/docking_scores.csv").write_text(
        "compound_id\tstatus\tpose_rank\tbest_score_kcal_mol\tpose_path\n"
        "CMP1\tcompleted\t1\t-7.1\tartifacts/poses/pose-1.pdbqt\n"
        "CMP1\tcompleted\t2\t-8.0\tartifacts/poses/pose-2.pdbqt\n"
    )

    assert ingest_run(run, database=database) == 1
    row = next(row for row in database.list_molecules(A, collection="all")["items"] if row["name"] == "Ethanol")
    assert len(row["structures"]) == 1
    assert row["structures"][0]["context"]["database_structure_policy"] == "best_docking_pose"
    assert row["structures"][0]["context"]["portable"] is True
    assert Path(row["structures"][0]["path"]).suffix == ".pdbqt"
    assert "top-ranked pose" in Path(row["structures"][0]["path"]).read_text()
    assert len(row["observations"]) == 2


def test_registered_rosettaligand_docking_keeps_selected_pose_and_evaluation(database):
    from mn_molecular_design_studio import ligand_compounds

    imported = database.import_records(
        A, [{"name": "Ethanol", "components": molecule_imports.parse_smiles("CCO")[0]["components"],
             "metadata": {"compound_id": "RL-1"}}], name="Rosetta compounds", project_id="project-1",
    )
    molecule_id = imported["molecule_ids"][0]
    run = _make_run(database, "ligand-docking", "rosetta-ligand-real-workflow", tool="rosetta_ligand",
                    workflow="ligand.docking.rosetta_ligand")
    ligand_compounds.write_compound_set(
        run / "artifacts/input/compound_set.json", name="Rosetta inputs",
        compounds=[{"compound_id": "RL-1", "compound_name": "Ethanol", "smiles": "CCO"}],
        owner_id=A["id"], project_id="project-1", created_at="now",
    )
    pose = run / "artifacts/poses/rosetta-ligand.pdbqt"
    pose.parent.mkdir(parents=True)
    pose.write_text("MODEL 1\nREMARK SMILES CCO\nREMARK SMILES IDX 1 3 2 2 3 1\n"
                    "ATOM      1 O    UNL     1       8.000   3.000   4.000  0.000 O\n"
                    "ATOM      2 C    UNL     1       7.000   3.000   4.000  0.000 C\n"
                    "ATOM      3 C    UNL     1       6.000   3.000   4.000  0.000 C\nENDMDL\n")
    (run / "artifacts/docking_scores.csv").write_text(
        "compound_id,status,pose_rank,rosetta_total_score,pose_path\n"
        "RL-1,completed,1,-12.5,artifacts/poses/rosetta-ligand.pdbqt\n"
    )

    assert ingest_run(run, database=database, force=True) == 1
    row = next(item for item in database.list_molecules(A, collection="all")["items"]
               if item["id"] == molecule_id)
    assert len(row["representations"]) == 1
    assert row["representations"][0]["context"]["pose_rank"] == "1"
    assert any(item["namespace"] == "rosetta_ligand" and item["metrics"].get("rosetta_total_score") == -12.5
               for item in row["observations"])


def test_docking_evaluates_remote_only_compounds_in_a_local_evaluation_group(database):
    from mn_molecular_design_studio.cluster.replication import Replica
    from mn_molecular_design_studio.molecular_database_peer import apply_record, capture
    from mn_molecular_design_studio.molecular_database import _id
    from mn_molecular_design_studio import ligand_compounds

    owner_node = "remote-compound-owner"
    remote = MolecularDatabase(database.root.parent / "remote-owner-database")
    imported = remote.import_records(
        A, [{"name": "Remote ethanol", "components": molecule_imports.parse_smiles("CCO")[0]["components"],
             "metadata": {"compound_id": "CMP-REMOTE"}}],
        name="Remote compound collection", project_id="project-1",
        origin={"collection_kind": "curated", "owner_id": A["id"]},
    )
    molecule_id, owner_group = imported["molecule_ids"][0], imported["group_id"]
    replica = Replica()
    capture(replica, database=remote, node_id=owner_node)
    with replica.db() as connection:
        records = [(row["key"], row["kind"], json.loads(row["value"])) for row in connection.execute(
            "SELECT key,kind,value FROM records WHERE author=? ORDER BY kind,key", (owner_node,)
        )]
    for key, kind, value in sorted(records, key=lambda item: (item[1] != "molecule_group", item[0])):
        apply_record(kind, key, value, owner_node, database=database)
    peer_group = _id("peer-", owner_node + ":" + owner_group)
    before = database.list_molecules(A, collection="all")
    assert before["total"] == 1
    assert before["items"][0]["id"] == molecule_id
    assert before["items"][0]["groups"][0]["id"] == peer_group

    run = _make_run(database, "vina", "remote-best-pose", tool="vina", workflow="ligand.docking")
    ligand_compounds.write_compound_set(
        run / "artifacts/input/compound_set.json", name="Screen inputs",
        compounds=[{"compound_id": "CMP-REMOTE", "compound_name": "Remote ethanol", "smiles": "CCO"}],
        owner_id=A["id"], project_id="project-1", created_at="now",
    )
    poses = run / "artifacts/poses"
    poses.mkdir(parents=True)
    pose = "\n".join(["MODEL 1", "REMARK SMILES CCO", "REMARK SMILES IDX 1 3 2 2 3 1",
                      "ATOM      1 O    UNL     1       8.000   3.000   4.000  0.000 O",
                      "ATOM      2 C    UNL     1       7.000   3.000   4.000  0.000 C",
                      "ATOM      3 C    UNL     1       6.000   3.000   4.000  0.000 C", "ENDMDL", ""])
    (poses / "pose.pdbqt").write_text(pose)
    (run / "artifacts/docking_scores.csv").write_text(
        "compound_id\tstatus\tpose_rank\tbest_score_kcal_mol\tpose_path\n"
        "CMP-REMOTE\tcompleted\t1\t-7.1\tartifacts/poses/pose.pdbqt\n"
    )

    assert ingest_run(run, database=database) == 1
    after = database.list_molecules(A, collection="all")
    assert after["total"] == 1
    molecule = after["items"][0]
    assert molecule["id"] == molecule_id
    assert peer_group in {group["id"] for group in molecule["groups"]}
    evaluation_groups = [group for group in database.groups(A, collection="all")
                         if group["origin"].get("collection_kind") == "evaluation"]
    assert len(evaluation_groups) == 1
    assert molecule["observations"]
    assert any(row["namespace"] == "vina" for row in molecule["observations"])
    with database.connection() as connection:
        assert connection.execute(
            "SELECT count(*) FROM molecule_members WHERE group_id=? AND molecule_id=?",
            (peer_group, molecule_id),
        ).fetchone()[0] == 1
        assert connection.execute(
            "SELECT count(*) FROM molecule_observations WHERE group_id=? AND molecule_id=?",
            (peer_group, molecule_id),
        ).fetchone()[0] == 0
        assert connection.execute(
            "SELECT count(*) FROM molecule_groups WHERE owner_id=? AND id=?",
            (A["id"], owner_group),
        ).fetchone()[0] == 0


def test_generator_contracts_ingest_common_fields_and_preserve_output_identity(database):
    definitions = [
        ("bindcraft2_generator", "generation.backbone_sequence", "ACD", {"pTM": .8, "binder_ca_radius_of_gyration": 10}),
        ("boltzgen_generator", "generation.backbone_sequence", "ACE", {"binder_ca_radius_of_gyration": 11}),
        ("rfdiffusion3_foundry", "generation.backbone", None, {"binder_ca_radius_of_gyration": 12}),
    ]
    for index, (tool, stage, sequence, metrics) in enumerate(definitions):
        candidate = {"candidate_id": f"design-{index}", "source_tool": tool, "stage": stage,
                     "binder_sequence": sequence, "binder_chains": ["A"], "metrics": metrics}
        if sequence:
            candidate["binder_pdb"] = f"binder-{index}.pdb"
            candidate["metrics"]["native_pass_filters"] = True
        if not sequence:
            candidate["binder_pdb"] = "backbone.pdb"
        run = _make_run(database, "design", f"generator-{index}", tool=tool, workflow="design_generation", candidates=[candidate])
        if not sequence:
            (run / "backbone.pdb").write_text(_pdb_acd())
        else:
            (run / f"binder-{index}.pdb").write_text(_pdb_for_sequence(sequence))
        assert ingest_run(run, database=database) == (1 if sequence else 0)
        assert ingest_run(run, database=database) == 0
    page = database.list_molecules(A, sort="common.geometry.binder.ca_rg", descending=False)
    assert page["total"] == 2
    assert [row["scores"]["common.geometry.binder.ca_rg"] for row in page["items"]] == [10, 11]
    first = page["items"][0]
    assert first["scores"]["common.confidence.generator.complex.ptm"] == .8
    assert first["metadata"]["output_kind"] == "designed_structure"
    assert len(first["observations"]) == 1
    assert first["components"][0]["role"] == "binder"
    assert len(database.list_molecules(A, collection="library")["items"]) == 2
    assert not any(row["metadata"].get("output_kind") == "backbone_proposal" for row in page["items"])


def test_sequence_only_design_is_not_identified_as_target_structure(database):
    from mn_molecular_design_studio.molecular_database_ingest import _record
    target = database.root / "target.pdb"
    target.write_text(_pdb_acd())
    record = _record({"candidate_id": "new-design", "stage": "sequence_design", "binder_sequence": "KKK",
                      "binder_chains": ["Z"], "target_chains": ["A"], "target_pdb": str(target),
                      "source_tool": "esm_if1", "raw_metadata": {"source_tool": "rfdiffusion_generator"}}, database.root)
    assert record["components"][0]["sequence"] == "KKK"
    assert record["components"][0]["chain"] == "Z"
    assert record["components"][0]["role"] == "binder"
    assert record["metadata"]["source_tool"] == "esm_if1"
    assert record["files"] == []


def test_target_role_maps_through_synthetic_small_molecule_chain_id(database):
    from mn_molecular_design_studio.molecular_database_ingest import _record
    candidate = {
        "candidate_id": "ligand-binder",
        "stage": "generation.backbone_sequence",
        "binder_sequence": "ACD",
        "binder_chains": ["A"],
        "target_chains": ["C"],
        "raw_metadata": {"components": [
            {"type": "small_molecule", "chain": "C_108", "residue_key": "C:108", "ccd": "RAP"},
            {"type": "protein", "chain": "A", "sequence": "ACD"},
        ]},
    }
    record = _record(candidate, database.root)
    components = {component["chain"]: component for component in record["components"]}
    assert components["A"]["role"] == "binder"
    assert components["C_108"]["role"] == "target"
    assert components["C_108"]["residue_key"] == "C:108"


def test_matching_candidate_fills_missing_member_roles_without_overwriting_existing_roles(database):
    from mn_molecular_design_studio.molecular_database_ingest import _merge_member_roles
    imported = database.import_records(A, [{"name": "complex", "components": [
        {"type": "small_molecule", "chain": "C_108", "residue_key": "C:108", "ccd": "RAP"},
        {"type": "protein", "chain": "A", "sequence": "ACD", "role": "binder"},
    ]}])
    molecule_id = imported["molecule_ids"][0]
    _merge_member_roles(database, imported["group_id"], molecule_id, [
        {"type": "small_molecule", "chain": "C_108", "residue_key": "C:108", "ccd": "RAP", "role": "target"},
        {"type": "protein", "chain": "A", "sequence": "ACD", "role": "binder"},
    ])
    row = database.list_molecules(A, group_id=imported["group_id"])["items"][0]
    roles = {component["chain"]: component["role"] for component in row["components"]}
    assert roles == {"C_108": "target", "A": "binder"}


def test_repeated_evaluations_reuse_original_group_and_accumulate_once(database):
    components = [{"chain": "A", "type": "protein", "sequence": "ACDE"}]
    imported = database.import_records(A, [{"name": "input", "components": components}], name="Original", project_id="project-1")
    molecule_id = imported["molecule_ids"][0]
    candidate = {"candidate_id": molecule_id, "stage": "complex_refolding", "binder_sequence": "ACDE",
                 "metrics": {"confidence": 0.8}, "raw_metadata": {"studio_molecule_id": molecule_id, "components": components}}
    run_one = _make_run(database, "optimized-prediction", "eval-1", candidates=[candidate])
    run_two = _make_run(database, "optimized-prediction", "eval-2", candidates=[candidate])

    assert ingest_run(run_one, database=database) == 1
    assert ingest_run(run_one, database=database) == 0
    assert ingest_run(run_two, database=database) == 1
    rows = database.list_molecules(A, group_id=imported["group_id"], collection="library")["items"]
    assert len(rows) == 1
    assert len(rows[0]["observations"]) == 2
    assert {row["run_id"] for row in rows[0]["observations"]} == {"eval-1", "eval-2"}
    assert {row["id"] for row in database.groups(A, collection="all")} == {imported["group_id"]}


@pytest.mark.parametrize("mapping_path", [
    "artifacts/raw/optimized_prediction/input_mapping.json",
    "artifacts/raw/af3_jax_cofolding/input_mapping.json",
])
def test_batched_prediction_inputs_restore_compound_identity_for_database(database, mapping_path):
    candidates = [
        {"candidate_id": "batch-1-model", "source_tool": "optimized_prediction_boltz2",
         "stage": "complex_refolding", "complex_pdb": "compound-1.pdb",
         "metrics": {"ranking_score": .8}, "raw_metadata": {"batch_item_id": "task-1"}},
        {"candidate_id": "batch-2-model", "source_tool": "optimized_prediction_boltz2",
         "stage": "complex_refolding", "complex_pdb": "compound-2.pdb",
         "metrics": {"ranking_score": .9}, "raw_metadata": {"batch_item_id": "task-2"}},
    ]
    run = _make_run(database, "refolding-validation", "compound-batch", candidates=candidates)
    for filename in ("compound-1.pdb", "compound-2.pdb"):
        (run / filename).write_text(_pdb_acd_with_generic_ligand())
    (run / "input.json").write_text(json.dumps({
        "tool": "optimized_prediction_boltz2",
        "workflow": "refolding.optimized_prediction",
        "inputs": {"tasks": [
            {"name": "Ethanol", "compound_id": "CHEM-101", "compound_name": "Ethanol",
             "ligand_smiles": [{"id": "B", "smiles": "CCO"}]},
            {"name": "Cyclopropane", "compound_id": "CHEM-102", "compound_name": "Cyclopropane",
             "ligand_smiles": [{"id": "B", "smiles": "C1CC1"}]},
        ]},
        "params": {"engine": "boltz2"},
    }))
    mapping_file = run / mapping_path
    mapping_file.parent.mkdir(parents=True, exist_ok=True)
    mapping_file.write_text(json.dumps({"items": [
        {"id": "task-1", "name": "Ethanol"},
        {"id": "task-2", "name": "Cyclopropane"},
    ]}))

    assert ingest_run(run, database=database) == 2
    assert ingest_run(run, database=database) == 0
    records = database.list_molecules(A, collection="all")["items"]
    assert {record["name"] for record in records} == {"Ethanol", "Cyclopropane"}
    ligand_smiles = {
        record["name"]: next(component["smiles"] for component in record["components"] if component["type"] == "small_molecule")
        for record in records
    }
    assert ligand_smiles == {"Ethanol": "CCO", "Cyclopropane": "C1CC1"}
    assert {record["metadata"]["compound_id"] for record in records} == {"CHEM-101", "CHEM-102"}


def test_save_run_selection_is_exact_repeatable_and_copies_run_context(database):
    from mn_molecular_design_studio.molecular_database_results import list_run_molecules, save_run_selection

    run = _make_run(database, "proteinmpnn", "results", tool="proteinmpnn", workflow="sequence_design",
                    candidates=[
                        {"candidate_id": "chosen", "stage": "sequence_design", "binder_sequence": "ACD",
                         "binder_pdb": "chosen.pdb", "metrics": {"rank": 1}},
                        {"candidate_id": "excluded", "stage": "sequence_design", "binder_sequence": "ACDE",
                         "metrics": {"rank": 2}},
                    ])
    (run / "chosen.pdb").write_text(_pdb_acd())
    assert ingest_run(run, database=database) == 0
    run_records = list_run_molecules(run, A, database=database)
    assert len(run_records) == 2
    chosen = next(row for row in run_records if row["name"] == "chosen")
    assert chosen["result_only"] is True
    assert chosen["representations"]

    first = save_run_selection(run, A, [chosen["id"]], database=database, name="Saved", project_id="")
    again = save_run_selection(run, A, [chosen["id"]], database=database, group_id=first["group_id"], project_id="")
    second = save_run_selection(run, A, [chosen["id"]], database=database, name="Saved again", project_id="")
    assert first["molecule_ids"] == again["molecule_ids"] == second["molecule_ids"] == [chosen["molecule_id"]]
    assert first["saved_count"] == again["saved_count"] == second["saved_count"] == 1
    assert database.group(first["group_id"], A)["project_id"] == ""
    assert database.list_molecules(A, group_id=first["group_id"])["total"] == 1
    assert database.list_molecules(A, group_id=second["group_id"])["total"] == 1

    with database.connection() as connection:
        for destination in (first["group_id"], second["group_id"]):
            assert connection.execute("SELECT count(*) FROM molecule_observations WHERE group_id=?", (destination,)).fetchone()[0] == 1
            assert connection.execute("SELECT count(*) FROM molecule_structures WHERE group_id=?", (destination,)).fetchone()[0] == 1
            copied_sources = connection.execute("SELECT count(*) FROM molecule_sources WHERE group_id=? AND json_extract(provenance,'$.saved_selection')=1", (destination,)).fetchone()[0]
            assert copied_sources == 1
        members = connection.execute("SELECT molecule_id FROM molecule_members WHERE group_id=?", (first["group_id"],)).fetchall()
    assert [row[0] for row in members] == [chosen["molecule_id"]]


def test_results_passing_collection_preserves_engine_model_and_is_idempotent(database):
    from mn_molecular_design_studio.molecular_database_results import list_run_molecules, save_run_selection

    candidates = [
        {"candidate_id": "binder-pass", "source_tool": "boltz2",
         "stage": "complex_refolding", "complex_pdb": "pass.pdb", "binder_sequence": "ACD",
         "binder_chains": ["B"], "target_chains": ["A"],
         "parents": ["rfd4-parent-9"], "metrics": {"ipTM": 0.81},
         "raw_metadata": {"model_id": "boltz-model-3", "generator_engine": "RFD4",
                          "prediction_engine": "Boltz-2", "target_identity": "PD-L1 construct X",
                          "validation_acceptance": {"status": "accepted", "policy": {"version": "accept-v3"}}}},
        {"candidate_id": "binder-reject", "source_tool": "boltz2",
         "stage": "complex_refolding", "complex_pdb": "reject.pdb", "binder_sequence": "ACD",
         "binder_chains": ["B"], "target_chains": ["A"],
         "raw_metadata": {"model_id": "boltz-model-4", "prediction_engine": "Boltz-2", "target_identity": "PD-L1 construct X",
                          "validation_acceptance": {"status": "rejected", "policy": {"version": "accept-v3"}}}},
    ]
    run = _make_run(database, "refolding-validation", "passing-collection", candidates=candidates)
    (run / "pass.pdb").write_text(_pdb_efg_acd_complex())
    (run / "reject.pdb").write_text(_pdb_efg_acd_complex())

    rows = list_run_molecules(run, A, database=database)
    passing = next(row for row in rows if row.get("candidate_id") == "binder-pass")
    collection_metadata = {
        "acceptance_stage": "Final computational acceptance", "acceptance_stage_id": "final",
        "acceptance_policies": [{"version": "accept-v3"}],
        "target_identity": "PD-L1 construct X", "generator": "RFD4", "source_workstation": "pc004033",
        "accepted_candidate_models": [{"candidate_id": "binder-pass", "model_id": "boltz-model-3"}],
    }
    candidate_models = [{
        "candidate_id": "binder-pass", "model_id": "boltz-model-3",
        "decision": {"stage_id": "final", "cohort_status": "passed", "policy": {"version": "accept-v3"}},
    }]
    first = save_run_selection(run, A, [passing["id"]], database=database, name="Passing designs",
                               candidate_models=candidate_models, collection_metadata=collection_metadata)
    repeated = save_run_selection(run, A, [passing["id"]], database=database, name="Passing designs",
                                  candidate_models=candidate_models, collection_metadata=collection_metadata)

    assert repeated["group_id"] == first["group_id"]
    assert first["molecule_ids"] == [passing["molecule_id"]]
    assert database.list_molecules(A, group_id=first["group_id"])["total"] == 1
    with database.connection() as connection:
        member_structures = connection.execute(
            "SELECT id,path,context FROM molecule_structures WHERE group_id=?", (first["group_id"],)
        ).fetchall()
        sources = connection.execute(
            "SELECT provenance FROM molecule_sources WHERE group_id=? AND json_extract(provenance,'$.saved_selection')=1",
            (first["group_id"],),
        ).fetchall()
    assert len(member_structures) == 1
    structure_context = json.loads(member_structures[0]["context"])
    assert structure_context["candidate_id"] == "binder-pass"
    assert structure_context["candidate_model_id"] == "boltz-model-3"
    assert structure_context["model_id"] == "0"  # Actual coordinate-model selector remains distinct.
    assert structure_context["generator_engine"] == "RFD4"
    assert structure_context["prediction_engine"] == "Boltz-2"
    assert structure_context["target_identity"] == "PD-L1 construct X"
    assert structure_context["parent_ids"] == ["rfd4-parent-9"]
    assert structure_context["portable"] is True
    assert len(sources) == 1
    source = json.loads(sources[0]["provenance"])
    assert source["acceptance_stage"] == "Final computational acceptance"
    assert source["acceptance_decision"]["policy"]["version"] == "accept-v3"
    assert source["target_identity"] == "PD-L1 construct X"
    assert source["generator_engine"] == "RFD4"
    assert source["prediction_engine"] == "Boltz-2"


def test_results_acceptance_inventory_keeps_final_native_and_unknown_stages_separate(database):
    from mn_molecular_design_studio.molecular_database_results import run_acceptance_inventory

    run = _make_run(database, "design", "acceptance-inventory", tool="rfd4",
                    candidates=[
        {"candidate_id": "final-pass", "metrics": {"native_pass_filters": False},
         "raw_metadata": {"model_id": "m1", "validation_acceptance": {"status": "accepted", "policy": {"version": "final-v2"}}}},
        {"candidate_id": "final-rejected", "metrics": {"native_pass_filters": True},
         "raw_metadata": {"model_id": "m2", "validation_acceptance": {"status": "rejected", "policy": {"version": "final-v2"}}}},
        {"candidate_id": "unknown", "metrics": {"confidence": 0.99}, "raw_metadata": {"model_id": "m3"}},
    ])

    inventory = run_acceptance_inventory(run, A)
    by_stage = {row["stage_id"]: row for row in inventory["stages"]}
    assert inventory["candidate_count"] == 3
    assert by_stage["final"]["passed_model_count"] == 1
    assert by_stage["final"]["failed_model_count"] == 1
    assert by_stage["final"]["unrecorded_model_count"] == 1
    assert by_stage["final"]["passed"][0]["model_id"] == "m1"
    assert by_stage["native"]["passed_model_count"] == 1
    assert by_stage["native"]["failed_model_count"] == 1
    assert by_stage["native"]["passed"][0]["candidate_id"] == "final-rejected"
    assert by_stage["native"]["failed"][0]["candidate_id"] == "final-pass"


def test_run_result_acl_destination_schema_and_project_are_enforced(database):
    from mn_molecular_design_studio.molecular_database_results import list_run_molecules, save_run_selection

    run = _make_run(database, "proteinmpnn", "protected", tool="proteinmpnn", workflow="sequence_design",
                    candidates=[{"candidate_id": "selected", "stage": "sequence_design", "binder_sequence": "ACD",
                                 "metrics": {"rank": 1}}])
    row = list_run_molecules(run, A, database=database)[0]
    with pytest.raises(PermissionError):
        list_run_molecules(run, B, database=database)
    with pytest.raises(PermissionError):
        save_run_selection(run, B, [row["id"]], database=database, project_id="")

    constrained = database.create_group(A, "Constrained", project_id="project-1", import_method="user_uploads",
                                        origin={"type": "user", "collection_kind": "curated"})
    database.set_schema(A, constrained, {"fields": [{"name": "required", "type": "text", "required": True}]})
    with pytest.raises(ValueError, match="Required property is missing"):
        save_run_selection(run, A, [row["id"]], database=database, group_id=constrained, project_id="project-1")
    assert database.list_molecules(A, group_id=constrained)["total"] == 0
    with pytest.raises(PermissionError):
        save_run_selection(run, A, [row["id"]], database=database, project_id="project-2")


def test_run_results_api_lists_and_promotes_only_selected_ids(database, monkeypatch):
    from fastapi.testclient import TestClient
    from mn_molecular_design_studio.api import app, _current_user
    from mn_molecular_design_studio import molecular_database_api
    from mn_molecular_design_studio import runtime

    monkeypatch.setattr(molecular_database_api, "sync_workspace_async", lambda *args, **kwargs: [])
    # The application module may already have loaded a site-level runtime path
    # before this isolated test fixture set its environment variables.
    monkeypatch.setattr(runtime, "runs_root", lambda: database.root.parent / "runs")
    run = _make_run(database, "proteinmpnn", "api-results", tool="proteinmpnn", workflow="sequence_design",
                    candidates=[{"candidate_id": "first", "stage": "sequence_design", "binder_sequence": "ACD"},
                                {"candidate_id": "second", "stage": "sequence_design", "binder_sequence": "ACDE"}])
    from mn_molecular_design_studio.molecular_database_results import list_run_molecules
    assert len(list_run_molecules(run, A, database=database)) == 2
    app.dependency_overrides[_current_user] = lambda: A
    try:
        with TestClient(app) as client:
            response = client.get("/api/v1/molecules/jobs/proteinmpnn/api-results")
            assert response.status_code == 200, response.text
            selected = response.json()["items"][0]["id"]
            saved = client.post("/api/v1/molecules/jobs/proteinmpnn/api-results/selections", json={
                "molecule_ids": [selected], "name": "Only selected", "project_id": "",
            })
            assert saved.status_code == 200, saved.text
            assert saved.json()["molecule_ids"] == [response.json()["items"][0]["molecule_id"]]
            assert database.list_molecules(A, group_id=saved.json()["group_id"])["total"] == 1
            linked = client.get("/api/v1/molecules/jobs/proteinmpnn/api-results")
            assert linked.status_code == 200
            assert [g["id"] for g in linked.json()["saved_collections"]] == [saved.json()["group_id"]]
            app.dependency_overrides[_current_user] = lambda: B
            denied = client.get("/api/v1/molecules/jobs/proteinmpnn/api-results")
            assert denied.status_code == 404
    finally:
        app.dependency_overrides.pop(_current_user, None)


def test_results_listing_is_a_read_and_leaves_candidates_unindexed(database, monkeypatch):
    from mn_molecular_design_studio import molecular_database_ingest
    from mn_molecular_design_studio.molecular_database_results import list_run_molecules

    run = _make_run(database, "proteinmpnn", "read-only-results", tool="proteinmpnn", workflow="sequence_design",
                    candidates=[{"candidate_id": "candidate-1", "stage": "sequence_design",
                                 "binder_sequence": "ACDE", "metrics": {"rank": 1}}])
    monkeypatch.setattr(molecular_database_ingest, "ingest_run",
                        lambda *args, **kwargs: pytest.fail("Results listing must not ingest candidates"))
    rows = list_run_molecules(run, A, database=database)
    assert len(rows) == 1
    assert rows[0]["result_only"] is True
    assert database.groups(A, collection="all") == []


def test_sync_throttle_preserves_last_indexing_warning(database, monkeypatch):
    from mn_molecular_design_studio import molecular_database_ingest as ingest

    monkeypatch.setattr(ingest, "migrate_legacy_groups", lambda *args: (_ for _ in ()).throw(RuntimeError("migration issue")))
    monkeypatch.setattr("mn_molecular_design_studio.auth.accessible_project_ids", lambda user: [])
    monkeypatch.setattr("mn_molecular_design_studio.job_catalog.search_jobs",
                        lambda **kwargs: {"items": [], "next_cursor": None})

    first = ingest.sync_workspace(A, force=True, database=database)
    assert first == [{"run_id": "legacy-migration", "error": "migration issue"}]
    # A normal UI rerun is throttled, but must continue showing the warning.
    assert ingest.sync_workspace(A, database=database) == first

    monkeypatch.setattr(ingest, "migrate_legacy_groups", lambda *args: None)
    assert ingest.sync_workspace(A, force=True, database=database) == []
    assert ingest.sync_workspace(A, database=database) == []


def test_saved_collection_is_navigable_from_run_but_not_removed_with_automatic_outputs(database, monkeypatch):
    from mn_molecular_design_studio.molecular_database_results import list_saved_run_collections, remove_run_records
    from mn_molecular_design_studio.cluster import config
    run = _make_run(database, 'design', 'original-run', project='', tool='rfd4_proteina', candidates=[])
    monkeypatch.setattr(config, 'identity', lambda: 'source-node')
    origin = {'type':'user', 'collection_kind':'curated', 'run_id':run.name, 'task_group':'design'}
    group = database.create_group(A, 'Accepted binders', origin=origin)
    record = {'name':'passing', 'components':molecule_imports.parse_sequences('ACDEFG')[0]['components'],
              'files':[{'name':'complex.pdb', 'content':b'HEADER existing conformation\nEND\n'}]}
    provenance = {**origin, 'source_node':'source-node'}
    imported = database.import_records(A, [record], group_id=group, origin=provenance)
    assert 'saved_selection' not in provenance  # Do not mutate callers' provenance.
    dbrow = database.list_molecules(A, group_id=group)['items'][0]
    sid = dbrow['representations'][0]['id']
    other = database.create_group(A, 'Same run ID from another machine', origin=origin)
    database.import_records(A, [record], group_id=other, origin={**provenance,'source_node':'other-node'})
    private = database.create_group(B, 'Inaccessible selection', origin=origin)
    database.import_records(B, [record], group_id=private, origin=provenance)
    before = list_saved_run_collections(run, A, database=database)
    assert [g['id'] for g in before] == [group]
    assert before[0]['molecule_ids'] == imported['molecule_ids']
    remove_run_records('design', run.name, A, owner_id=A['id'], project_id='', database=database)
    assert list_saved_run_collections(run, A, database=database) == before
    assert database.structure_file(A, sid).read_bytes() == b'HEADER existing conformation\nEND\n'
