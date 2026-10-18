"""Structure, constraint, queue and graph contracts; no model inference."""
import copy
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from Bio.PDB import Atom, Chain, MMCIFIO, Residue

from mn_molecular_design_studio import inverse_folding_input as inputs
from mn_molecular_design_studio.design_validation import read_structure
from mn_molecular_design_studio.engine.workflow_engine.workflows import inverse_folding as workflow


@pytest.fixture
def source(design_source_factory):
    run, candidate = design_source_factory(missing_occupancy=False)
    path = run / "designed.cif"
    model = next(read_structure(path).get_models())
    for chain in model:
        for residue in chain:
            ca = residue["CA"].coord
            for name, offset, element in [("N", [-1., 0, 0], "N"), ("C", [1., 0, 0], "C"),
                                            ("O", [1., 1., 0], "O"), ("CB", [0., -1., 0], "C")]:
                residue.add(Atom.Atom(name, ca + offset, 90., 1., " ", name.rjust(4), 1, element=element))
    ligand = Chain.Chain("L")
    residue = Residue.Residue(("H_NAG", 101, " "), "NAG", " ")
    residue.add(Atom.Atom("C1", np.array([1., 2., 3.]), 10., 1., " ", " C1 ", 1, element="C"))
    ligand.add(residue)
    model.add(ligand)
    writer = MMCIFIO()
    writer.set_structure(model)
    writer.save(str(path))
    candidate["complex_pdb"] = str(path)
    collection = run / "artifacts/normalized_candidates/candidates.jsonl"
    collection.write_text(json.dumps(candidate) + "\n")
    return run, candidate, collection, path


@pytest.mark.parametrize("engine", inputs.ENGINES)
@pytest.mark.parametrize("bad", [{"num_sequences": 0}, {"seed": True}, {"temperature": float("nan")},
    {"num_conformers": 501}, {"omit_aas": inputs.AA}, {"batch_size": 0}])
def test_invalid_requests_fail_before_a_job(engine, bad):
    with pytest.raises(ValueError):
        inputs.prepare(engine, bad)


def test_caliby_model_and_ensemble_combinations():
    for model in inputs.CALIBY_MODELS:
        assert inputs.prepare("caliby", dict(model_type=model, ensemble_mode="single"))["model_type"] == model
    with pytest.raises(ValueError, match="Distilled"):
        inputs.prepare("caliby", dict(model_type="caliby_distill"))
    with pytest.raises(ValueError, match="32-conformer"):
        inputs.prepare("caliby", dict(optimization_mode="exact", num_conformers=8))
    assert inputs.prepare("caliby", dict(optimization_mode="exact"))["num_conformers"] == 32
    with pytest.raises(ValueError, match="Scoring"):
        inputs.prepare("esm_if1", dict(task="score"))


def test_source_paths_and_empty_selection_are_validated(source, tmp_path):
    run, _, collection, _ = source
    p = dict(source_run_dir=run, candidates_jsonl=collection)
    with pytest.raises(ValueError, match="at least one"):
        inputs.prepare("esm_if1", p | dict(selected_candidate_ids=[]))
    with pytest.raises(ValueError, match="belong"):
        inputs.prepare("esm_if1", p | dict(source_run_dir=tmp_path / "another"))
    with pytest.raises(ValueError, match="together"):
        inputs.prepare("esm_if1", dict(source_run_dir=run))


def test_staging_preserves_author_mapping_and_never_changes_source(source, tmp_path):
    _, _, _, path = source
    model = next(read_structure(path).get_models())
    model["B"].id = "BINDER"
    model["BINDER"][(" ", 2, " ")].id = (" ", 22, "A")
    writer = MMCIFIO()
    writer.set_structure(model)
    writer.save(str(path))
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    rows, chains = workflow.stage_structure(path, tmp_path / "input.pdb", protein_only=False)
    insertion = next(r for r in rows if r["label"] == "BINDER22A")
    assert insertion["position"] == 2 and insertion["engine_label"] == chains["BINDER"] + "2"
    assert insertion["original_id"] == [" ", 22, "A"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
    assert "NAG" in (tmp_path / "input.pdb").read_text()


def test_connected_targets_and_fixed_motifs_are_protected(source):
    _, candidate, _, path = source
    candidate["raw_metadata"] = dict(fixed_motif_required=True, motif_mapping_verified=True, fixed_residues=["B2"])
    p = inputs.defaults("esm_if1")
    chains, mutable = workflow.design_positions(path, candidate, p)
    assert chains == ["B"] and mutable == {"B1", "B3", "B4", "B5"}
    with pytest.raises(ValueError, match="fixed target"):
        workflow.design_positions(path, candidate, p | dict(design_chains="A"))
    candidate["raw_metadata"]["motif_mapping_verified"] = False
    with pytest.raises(ValueError, match="verified"):
        workflow.design_positions(path, candidate, p)


def test_csv_constraints_cannot_override_fixed_targets(source, tmp_path):
    _, candidate, _, path = source
    mapping, _ = workflow.stage_structure(path, tmp_path / "input.pdb")
    csv = tmp_path / "constraints.csv"
    csv.write_text("pdb_key,fixed_pos_override_seq\ndesigned,A1:W\n")
    with pytest.raises(ValueError, match="designable"):
        workflow.map_constraints(csv, path, candidate["candidate_id"], mapping, {"B1"})
    csv.write_text('pdb_key,fixed_pos_seq,pos_restrict_aatype\ndesigned,B1,B2:AG\n')
    result = workflow.map_constraints(csv, path, candidate["candidate_id"], mapping, {"B2"})
    assert result == dict(fixed_pos_seq="B1", pos_restrict_aatype="B2:AG")


def prepared(source, tmp_path, **options):
    run, _, collection, _ = source
    p = inputs.prepare("caliby", dict(source_run_dir=run, candidates_jsonl=collection,
        ensemble_mode="single", num_sequences=1, designed_residues="B2", **options), require_input=True)
    destination = tmp_path / "job"
    destination.mkdir(exist_ok=True)
    items = workflow.prepare_items("caliby", p, destination)
    output = destination / "artifacts/inverse_folding/native_results.json"
    output.write_text(json.dumps([dict(key=items[0]["key"], sequences=dict(A="ACDEFG", B="KAMNP"), metrics=dict(caliby_energy=-3.))]))
    return destination, p, items


def test_normalized_outputs_keep_target_ligand_coordinates_and_remove_old_confidence(source, tmp_path):
    run, p, items = prepared(source, tmp_path)
    candidates = workflow.normalize_outputs("caliby", p, run, items)
    row = candidates[0]
    assert row["stage"] == "sequence_design" and row["binder_sequence"] == "KAMNP"
    assert row["binder_chains"] == ["B"] and row["target_chains"] == ["A"]
    assert "plddt" not in row["metrics"] and row["metrics"]["caliby_energy"] == -3.
    output = next(read_structure(Path(row["complex_pdb"])).get_models())
    original = next(read_structure(source[3]).get_models())
    assert "L" in output and output["L"][("H_NAG", 101, " ")].resname == "NAG"
    assert output["B"][(" ", 2, " ")].resname == "ALA"
    assert "CB" not in output["B"][(" ", 2, " ")]
    for atom, before in zip(output["A"].get_atoms(), original["A"].get_atoms()):
        assert atom.name == before.name and np.array_equal(atom.coord, before.coord)
    assert (run / "artifacts/inverse_folding/sequences.fasta").read_text().endswith("KAMNP\n")


@pytest.mark.parametrize("sequences", [dict(A="WCDEFG", B="KAMNP"), dict(A="ACDEFG", B="KAMN"), dict(A="ACDEFG", B="K<MASK>MNP")])
def test_bad_model_outputs_never_become_candidates(source, tmp_path, sequences):
    run, p, items = prepared(source, tmp_path)
    path = run / "artifacts/inverse_folding/native_results.json"
    path.write_text(json.dumps([dict(key=items[0]["key"], sequences=sequences)]))
    with pytest.raises(ValueError):
        workflow.normalize_outputs("caliby", p, run, items)
    assert not (run / "artifacts/normalized_candidates/candidates.jsonl").exists()


def test_short_batch_and_misaligned_ensembles_are_rejected(source, tmp_path):
    run, p, items = prepared(source, tmp_path)
    with pytest.raises(ValueError, match="requested number"):
        workflow.normalize_outputs("caliby", p | dict(num_sequences=2), run, items)
    altered = tmp_path / "other.cif"
    model = next(read_structure(source[3]).get_models())
    model["B"][(" ", 1, " ")].resname = "ALA"
    writer = MMCIFIO()
    writer.set_structure(model)
    writer.save(str(altered))
    with pytest.raises(ValueError, match="match every"):
        workflow.prepare_items("caliby", p | dict(ensemble_paths=[str(altered)]), run)


def test_caliby_accepts_topology_matched_backbone_only_ensemble(source, tmp_path):
    from Bio.PDB import PDBIO
    from Bio.SeqUtils import seq1

    _, _, _, primary = source
    model = next(read_structure(primary).get_models())
    for chain in model:
        for residue in chain:
            if residue.id[0] == " ":
                residue.resname = "GLY"
    backbone = tmp_path / "backbone_only.pdb"
    writer = PDBIO()
    writer.set_structure(model)
    writer.save(str(backbone))
    p = inputs.prepare("caliby", dict(structure_path=str(primary), design_chains="B",
        ensemble_mode="provided", ensemble_paths=[str(backbone)], num_sequences=1), require_input=True)
    run = tmp_path / "job"
    run.mkdir()
    items = workflow.prepare_items("caliby", p, run)
    staged = run / items[0]["ensemble_paths"][0].removeprefix("/work/")
    observed = next(read_structure(staged).get_models())
    for chain_id in ("A", "B"):
        expected = "".join(r["aa"] for r in items[0]["mapping"] if r["chain"] == chain_id)
        actual = "".join(seq1(res.resname, undef_code="X") for res in observed[chain_id])
        assert actual == expected


@pytest.mark.parametrize("engine", inputs.ENGINES)
def test_docker_invocation_uses_existing_job_and_local_weights_only(engine, source, tmp_path, monkeypatch):
    from mn_molecular_design_studio.engine.workflow_engine.workflows import sequence_design
    run = tmp_path / "existing-job"
    run.mkdir()
    monkeypatch.setattr(workflow, "create_job", lambda *a, **k: pytest.fail("An adopted queue job must not create another run"))
    monkeypatch.setattr(workflow, "update_status", lambda *a, **k: None)
    finished = []
    monkeypatch.setattr(workflow, "finish_job", lambda *a: finished.append(a))
    monkeypatch.setattr(inputs, "readiness", lambda *a: dict(ready=True))
    monkeypatch.setattr(inputs, "weights", lambda *a: [workflow.reference_root() / "model.pt"])
    commands = []
    monkeypatch.setattr(sequence_design, "_run_shell_steps", lambda _, steps: commands.extend(steps) or 1)
    function = {"esm_if1": workflow.run_esm_if1, "caliby": workflow.run_caliby, "boltzgen_inverse": workflow.run_boltzgen_inverse}[engine]
    result = function(source_run_dir=source[0], candidates_jsonl=source[2], gpu_device="0,1", existing_job=SimpleNamespace(run_dir=run))
    assert result == run and finished[-1][1] is False
    command = commands[0]["command"]
    assert '"device=0,1"' in command
    assert "HF_HUB_OFFLINE=1" in command and str(workflow.reference_root()) + ":/ref:ro" in command
    payload = json.loads((run / "artifacts/inverse_folding/native_request.json").read_text())
    assert payload["engine"] == engine and payload["items"][0]["design_chains"] == ["B"]
    assert payload["items"][0]["mutable"] == ["B1", "B2", "B3", "B4", "B5"]


@pytest.mark.parametrize("engine", inputs.ENGINES)
def test_api_schema_and_prediction_to_inverse_to_prediction_graph(engine, source, tmp_path):
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import workflow_schema, validate_workflow_request, AutomationError
    from mn_molecular_design_studio.pipeline_graph import new_graph, validate_graph, NODE_CATALOG
    from mn_molecular_design_studio.engine.workflow_engine.workflows.graph_pipeline import node_request
    schema = workflow_schema(inputs.WORKFLOWS[engine])
    assert {"structure_path", "source_run_dir", "selected_candidate_ids"} <= {p["name"] for p in schema["parameters"]}
    request = tmp_path / "request.json"
    request.write_text(json.dumps(dict(schema_version=1, workflow=inputs.WORKFLOWS[engine], parameters={})))
    with pytest.raises(AutomationError, match="Upload a structure"):
        validate_workflow_request(request)
    graph = new_graph("foundry")
    sequence = graph["nodes"][2]
    sequence["data"].update(kind=engine, label=NODE_CATALOG[engine]["label"], params={f["name"]: f["default"] for f in NODE_CATALOG[engine]["fields"]})
    validate_graph(graph)
    p = node_request(engine, sequence["data"]["params"], dict(run_dir=str(source[0]), candidates_jsonl=str(source[2])), "inverse node")
    assert p["workflow"] == inputs.WORKFLOWS[engine] and p["parameters"]["source_run_dir"] == str(source[0])
    assert {"backbone", "sequence", "prediction", "target"} <= set(NODE_CATALOG[engine]["inputs"][0]["accepts"])
    direct = node_request(engine, sequence["data"]["params"], dict(target=dict(target_pdb=str(source[3]), target_chains=["B"])), "direct structure")
    assert direct["parameters"]["structure_path"] == str(source[3]) and direct["parameters"]["design_chains"] == "B"
    # A prediction can be redesigned and then predicted again, independently
    # of which upstream generator created its input backbone.
    inverse = copy.deepcopy(sequence)
    inverse["id"] = "redesign-prediction"
    prediction = copy.deepcopy(graph["nodes"][3])
    prediction["id"] = "predict-again"
    graph["nodes"] += [inverse, prediction]
    old = graph["edges"].pop()
    def connect(a, b): return dict(id=a + "-" + b, source=a, target=b, sourceHandle="out", targetHandle="in")
    graph["edges"] += [connect(old["source"], inverse["id"]), connect(inverse["id"], prediction["id"]), connect(prediction["id"], old["target"])]
    validate_graph(graph)


def test_archive_entries_cannot_escape_output_directory(tmp_path):
    import zipfile
    archive = tmp_path / "inputs.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("../../outside.pdb", "ATOM\n")
        handle.writestr("other/readme.txt", "ignore")
    files = inputs.expand_structure_archive(archive, tmp_path / "extracted")
    assert len(files) == 1 and files[0].parent == tmp_path / "extracted"
    assert not (tmp_path / "outside.pdb").exists()


@pytest.mark.parametrize("engine", inputs.ENGINES)
def test_node_editor_reuses_the_tile_form(engine):
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_string(f'''
from mn_molecular_design_studio.inverse_folding_ui import render_inputs
render_inputs({engine!r}, owner_id="test", key_prefix="same", connected=True)
''').run()
    assert not app.exception
    assert app.number_input(key="same-field-count").value == inputs.defaults(engine)["num_sequences"]


def test_esm_native_sampler_receives_fixed_target_and_residue_masks(source, tmp_path, monkeypatch):
    import contextlib
    import sys
    from mn_molecular_design_studio.engine.workflow_engine.workflows import inverse_folding_native as native
    run, p, items = prepared(source, tmp_path)
    calls = []
    class Model:
        def eval(self): return self
        def to(self, device): return self
        def sample(self, coords, partial_seq, **kwargs):
            calls.append(partial_seq)
            return "".join("A" if v == "<mask>" else v for v in partial_seq)
    fake = SimpleNamespace(pretrained=SimpleNamespace(load_model_and_alphabet_local=lambda _: (Model(), None)))
    coords = {"A": np.zeros((6, 3, 3)), "B": np.zeros((5, 3, 3))}
    monkeypatch.setitem(sys.modules, "torch", SimpleNamespace(manual_seed=lambda _: None,
        cuda=SimpleNamespace(is_available=lambda: False), inference_mode=contextlib.nullcontext))
    monkeypatch.setitem(sys.modules, "esm", fake)
    monkeypatch.setitem(sys.modules, "esm.inverse_folding.multichain_util", SimpleNamespace(
        load_complex_coords=lambda *args: (coords, {}),
        _concatenate_coords=lambda context, chain: np.zeros((sum(len(v) for v in context.values()) + 10, 3, 3)),
        score_sequence_in_complex=lambda *args: (-1., -1.)))
    records = native.esm_if1(dict(parameters=p | dict(seed=123, multichain_backbone=True), items=items, checkpoint="/weights/model.pt"))
    assert calls[0][:5] == ["K", "<mask>", "M", "N", "P"]
    assert all(v == "<pad>" for v in calls[0][5:])
    assert records[0]["sequences"] == dict(A="ACDEFG", B="KAMNP")


@pytest.mark.parametrize("ensemble_mode,task", [("single", "design"), ("generate", "design"), ("single", "score")])
def test_caliby_native_api_constraints_ensembles_and_scoring(source, tmp_path, monkeypatch, ensemble_mode, task):
    import sys
    from mn_molecular_design_studio.engine.workflow_engine.workflows import inverse_folding_native as native
    run, p, items = prepared(source, tmp_path)
    p.update(seed=17, ensemble_mode=ensemble_mode, task=task, num_conformers=2)
    calls = []
    result = dict(seq=["ACDEFG:KAMNP"], U=[-3.])
    class Model:
        def sample(self, paths, **kwargs): calls.append(("sample", kwargs)); return result
        def ensemble_sample(self, paths, **kwargs): calls.append(("ensemble", kwargs)); return result
        def score(self, paths, **kwargs): calls.append(("score", kwargs)); return dict(seq=["ACDEFG:KLMNP"], U=[-2.])
    def conformers(paths, **kwargs):
        calls.append(("generate", kwargs))
        return {Path(paths[0]).stem: [paths[0]] * kwargs["num_samples_per_pdb"]}
    monkeypatch.setattr(native, "_clean_chain_map", lambda *a: dict(A="A", B="B"))
    # Container paths are translated to this test's job files; no engine runs.
    items[0]["path"] = str(run / items[0]["path"].removeprefix("/work/"))
    monkeypatch.setitem(sys.modules, "caliby", SimpleNamespace(load_model=lambda _: Model(), clean_pdbs=lambda paths, **kw: paths,
        make_constraints=lambda rows: rows, generate_ensembles=conformers))
    monkeypatch.setitem(sys.modules, "lightning", SimpleNamespace(seed_everything=lambda _: None))
    records = native.caliby_design(dict(parameters=p, items=items))
    assert records[0]["metrics"]["caliby_energy"] in {-3., -2.}
    method, kwargs = calls[-1]
    assert method == ("score" if task == "score" else "sample" if ensemble_mode == "single" else "ensemble")
    if task == "design":
        constraints = next(iter(kwargs["pos_constraint_df"].values()))["fixed_pos_seq"].split(",")
        assert "A1" in constraints and "A6" in constraints and "B1" in constraints and "B2" not in constraints
        assert kwargs["sampling_overrides"] == {"potts_sampling_cfg": {"potts_sweeps": 500}}
    if ensemble_mode == "generate":
        assert calls[0][1]["num_samples_per_pdb"] == 2


def test_boltzgen_native_uses_only_inverse_fold_configuration_for_foreign_structure(source, tmp_path, monkeypatch):
    import importlib.resources
    import sys
    from mn_molecular_design_studio.engine.workflow_engine.workflows import inverse_folding_native as native
    run, p, items = prepared(source, tmp_path)
    item = items[0]
    item["path"] = str(run / item["path"].removeprefix("/work/"))
    p.update(seed=15, temperature=.2, omit_aas="C")
    ref = tmp_path / "ref/boltzgen-cache/datasets--boltzgen--inference-data/snapshots/pin"
    ref.mkdir(parents=True)
    (ref / "mols.zip").write_bytes(b"unused")
    configs = []
    cfg = SimpleNamespace(data=SimpleNamespace(cfg=SimpleNamespace()), trainer=SimpleNamespace(),
        override=SimpleNamespace(inverse_fold_args=SimpleNamespace()))
    class Task:
        def run(self, config):
            output = Path(config.output)
            output.mkdir()
            writer = MMCIFIO()
            writer.set_structure(read_structure(Path(item["path"])))
            writer.save(str(output / "foreign_structure_0.cif"))
    monkeypatch.setitem(sys.modules, "boltzgen", SimpleNamespace())
    monkeypatch.setattr(importlib.resources, "files", lambda _: tmp_path)
    monkeypatch.setitem(sys.modules, "omegaconf", SimpleNamespace(OmegaConf=SimpleNamespace(
        load=lambda path: configs.append(str(path)) or cfg, save=lambda *args: None)))
    monkeypatch.setitem(sys.modules, "hydra", SimpleNamespace(utils=SimpleNamespace(instantiate=lambda config: Task())))
    monkeypatch.setitem(sys.modules, "pytorch_lightning", SimpleNamespace(seed_everything=lambda _: None))
    records = native.boltzgen_inverse(dict(parameters=p, items=items, checkpoint="/ref/ifold.ckpt", reference_root=str(tmp_path / "ref")))
    assert configs == [str(tmp_path / "resources/config/inverse_fold_only.yaml")]
    assert cfg.checkpoint == "/ref/ifold.ckpt" and cfg.data.cfg.multiplicity == 1
    assert cfg.override.inverse_fold_args.sampling_temperature == .2
    assert cfg.override.inverse_fold_args.inverse_fold_restriction == ["CYS"]
    import yaml
    spec = yaml.safe_load(Path(cfg.data.cfg.yaml_path[0]).read_text())
    assert spec["entities"][0]["file"]["design"] == [dict(chain=dict(id="B", res_index="2"))]
    assert records[0]["sequences"] == dict(A="ACDEFG", B="KLMNP")


@pytest.mark.parametrize("engine", inputs.ENGINES)
def test_source_project_guard_applies_to_all_inverse_engines(engine, source, monkeypatch):
    from fastapi import HTTPException
    from mn_molecular_design_studio import api, job_catalog
    monkeypatch.setattr(api, "_resolve_api_job", lambda *a: dict(task_group="design"))
    monkeypatch.setattr(job_catalog, "get_job_access_record", lambda *a: dict(project_id="another"))
    with pytest.raises(HTTPException, match="selected project"):
        api._guard_source_job(dict(workflow=inputs.WORKFLOWS[engine], parameters=dict(source_run_dir=str(source[0]))), dict(id="test"), project_id="chosen")


@pytest.mark.parametrize("engine", inputs.ENGINES)
def test_tool_tile_routes_to_its_custom_form(engine, monkeypatch):
    """Exercise the real dispatcher without starting the whole Streamlit app."""
    import ast
    import importlib
    from mn_molecular_design_studio import inverse_folding_ui
    source = Path(inverse_folding_ui.__file__).with_name("run_app.py").read_text()
    dispatch = next(n for n in ast.parse(source).body if isinstance(n, ast.FunctionDef) and n.name == "_render_workflow_setup")
    calls = []
    monkeypatch.setattr(importlib, "reload", lambda module: module)
    monkeypatch.setattr(inverse_folding_ui, "render_setup", lambda *args, **kwargs: calls.append((args, kwargs)))
    namespace = dict(st=SimpleNamespace(session_state=dict(selected_workflow=inputs.WORKFLOWS[engine])),
        TOOL_CATALOG_GROUPS={}, importlib=importlib, CURRENT_USER=dict(id="test"), workflow_schema=lambda name: dict(workflow=name),
        _render_project_selector=lambda *a: None, _submit_cofolding_job=lambda *a: None, search_jobs=lambda *a: None)
    exec(compile(ast.Module(body=[dispatch], type_ignores=[]), "run_app.py", "exec"), namespace)
    namespace["_render_workflow_setup"]()
    assert calls[0][0][0] == engine and calls[0][1]["owner_id"] == "test"
