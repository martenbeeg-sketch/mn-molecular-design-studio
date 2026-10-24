from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from mn_molecular_design_studio import ligand_compounds
from mn_molecular_design_studio.engine.workflow_engine.core import jobs
from mn_molecular_design_studio.engine.workflow_engine.core.docker_runner import build_docker_command
from mn_molecular_design_studio.engine.workflow_engine.services import workflow_automation
from mn_molecular_design_studio.engine.workflow_engine.workflows import ligand_tools
from mn_molecular_design_studio.engine.workflow_engine.workflows import autodock_gpu_runner
from mn_molecular_design_studio import ligand_tools_ui


def test_vina_box_from_points_centers_selected_atoms_with_padding():
    center, size = ligand_tools_ui._vina_box_from_points(
        [(-1.0, 2.0, 3.0), (3.0, 6.0, 9.0)], padding=4.0,
    )

    assert center == {"x": 1.0, "y": 4.0, "z": 6.0}
    assert size == {"x": 12.0, "y": 12.0, "z": 14.0}


def test_vina_ligand_preview_browses_compounds_in_a_set(tmp_path):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio.vina_ui import _render_ligand_2d_preview

    compound_path = tmp_path / "compounds.json"
    ligand_compounds.write_compound_set(
        compound_path,
        name="Preview set",
        compounds=[
            {"compound_id": "ethanol", "smiles": "CCO"},
            {"compound_id": "ethylamine", "smiles": "CCN"},
        ],
        owner_id="owner",
        project_id="project",
        created_at="now",
    )
    app = AppTest.from_string(f'''\
from mn_molecular_design_studio.vina_ui import _render_ligand_2d_preview
_render_ligand_2d_preview({{"path": {str(compound_path)!r}}})
''').run()
    assert not app.exception
    preview_key = hashlib.sha1(str(compound_path.resolve()).encode("utf-8")).hexdigest()[:12]
    assert any("ethanol" in block.value for block in app.markdown)
    assert any(block.value == "Compound 1 of 2" for block in app.caption)
    assert app.button(key=f"vina-ligand-preview-prev-{preview_key}").disabled

    app.button(key=f"vina-ligand-preview-next-{preview_key}").click().run()
    assert not app.exception
    assert any("ethylamine" in block.value for block in app.markdown)
    assert any(block.value == "Compound 2 of 2" for block in app.caption)
    assert app.button(key=f"vina-ligand-preview-next-{preview_key}").disabled


def test_vina_viewer_selection_uses_protein_residue_coordinates_only():
    coordinates = {
        ("A", 10): [(1.0, 2.0, 3.0)],
        ("B", 20): [(4.0, 5.0, 6.0)],
    }

    selected = ligand_tools_ui._selected_vina_protein_residues(
        {"sequenceSelections": [{"chainId": "A", "residues": [10, 11]}]},
        coordinates,
    )

    assert selected == {("A", 10)}


def test_vina_setup_loads_example_selects_pocket_and_submits(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import vina_inputs

    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path))
    app = AppTest.from_string('''\
import streamlit as st
from mn_molecular_design_studio import ligand_tools_ui as ui
from mn_molecular_design_studio import vina_inputs
from mn_molecular_design_studio import vina_ui

target, _ = vina_inputs.stage_example("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
pdb_lines = target.read_text().splitlines()
chain_a = [line for line in pdb_lines if line.startswith("ATOM  ") and line[21] == "A"]
target.write_text("\\n".join(
    [line for line in pdb_lines if line != "END"]
    + [line[:21] + "B" + line[22:] for line in chain_a]
    + ["END", ""]
))
st.session_state["vina-target-path"] = str(target)

def submit(schema, workflow, parameters, project, settings, title):
    st.session_state["test-vina-submission"] = {"workflow": workflow, "parameters": parameters}

ui._render_vina_setup(
    "ligand.docking.vina", {}, owner_id="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
    compound_sets=[], project_selector=lambda key: "project-1", submit=submit,
    title="AutoDock Vina",
)
''').run()
    assert not app.exception
    assert app.session_state["ligand-vina-pocket-mode"] == "Co-crystallized Ligand"
    assert app.session_state["ligand-vina-box-ligand"] == "STI|A|201|_"
    assert app.number_input(key="vina-box-size-x").value == pytest.approx(18.664)
    app.button(key="vina-load-example").click().run()
    assert not app.exception
    assert app.session_state["ligand-vina-pocket-mode"] == "Co-crystallized Ligand"
    assert app.session_state["ligand-vina-box-ligand"] == "STI|A|201|_"
    assert app.number_input(key="vina-box-center-x").value == pytest.approx(15.190)
    assert any("2D structure preview" in block.value for block in app.markdown)
    assert [tab.label for tab in app.tabs] == ["Co-crystallized Ligand", "Select Residues", "Predict (Boltz)", "Manual"]
    assert app.checkbox(key="ligand-vina-use-scrub").value is True
    assert app.number_input(key="ligand-vina-scrub-ph").value == pytest.approx(7.4)
    assert app.checkbox(key="ligand-vina-skip-tautomer").value is True
    assert app.number_input(key="ligand-vina-exhaustiveness").value == 32
    assert app.number_input(key="ligand-vina-num-modes").value == 10
    markdown = [block.value for block in app.markdown]
    assert markdown.index("**Ligand preparation**") < markdown.index("**Search settings**")
    assert app.number_input(key="vina-box-center-y").value == pytest.approx(53.903)
    old_focus_epoch = app.session_state.get("vina-viewer-focus-epoch", 0)
    app.selectbox(key="ligand-vina-box-ligand").set_value("STI|B|202|_").run()
    assert not app.exception
    assert app.selectbox(key="ligand-vina-box-ligand").value == "STI|B|202|_"
    assert app.session_state["vina-viewer-focus-epoch"] == old_focus_epoch + 1
    assert "STI|B|202|_" in app.session_state["vina-derived-box-signature"]
    assert app.number_input(key="vina-box-center-x").value == pytest.approx(13.171)
    assert app.number_input(key="vina-box-center-y").value == pytest.approx(95.113)
    assert app.number_input(key="vina-box-center-z").value == pytest.approx(58.594)
    app.button(key="ligand-vina-run").click().run()
    assert not app.exception
    submitted = app.session_state["test-vina-submission"]
    assert submitted["workflow"] == "ligand.docking.vina"
    assert submitted["parameters"]["center_y"] == pytest.approx(95.113)
    assert Path(submitted["parameters"]["compound_set_file"]).is_file()
    assert submitted["parameters"]["job_name"] == "autodock-vina-example"

    app.session_state["ligand-vina-pocket-mode"] = "Predict (Boltz)"
    app.run()
    chains = app.multiselect(key="vina-boltz-selected-chains")
    assert len(chains.options) == 2
    assert set(chains.value) == {"A", "B"}
    chain_sequence = "".join(letter for _, letter in vina_inputs.protein_sequence_entries(
        Path(app.session_state["vina-target-path"]).read_text()
    )["B"])
    chains.set_value(["B"]).run()
    assert not app.exception
    assert app.multiselect(key="vina-boltz-selected-chains").value == ["B"]
    assert "Imported compound imatinib-CID5291" in app.text_input(key="vina-probe-0").label
    assert app.text_input(key="vina-probe-id-0").value == "imatinib-CID5291"
    assert app.text_input(key="vina-probe-name-0").value == "Imatinib"
    imported_smiles = "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"
    assert app.text_input(key="vina-probe-0").value == imported_smiles
    app.button(key="vina-predict-pocket").click().run()
    assert not app.exception
    prediction = app.session_state["test-vina-submission"]
    assert prediction["workflow"] == "refolding.optimized_prediction"
    assert len(prediction["parameters"]["tasks"]) == 1
    assert prediction["parameters"]["tasks"][0]["protein_sequences"] == chain_sequence
    assert prediction["parameters"]["tasks"][0]["ligand_smiles"] == imported_smiles
    assert prediction["parameters"]["tasks"][0]["compound_id"] == "imatinib-CID5291"
    assert prediction["parameters"]["tasks"][0]["compound_name"] == "Imatinib"
    assert "output_format" not in prediction["parameters"]["engine_options"]
    from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import validate_prediction_request
    validate_prediction_request(prediction["parameters"])


def test_vina_residue_selection_requires_confirmation_to_update_box(tmp_path, monkeypatch):
    import hashlib
    from streamlit.testing.v1 import AppTest
    from mn_molecular_design_studio import vina_inputs

    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path))
    app = AppTest.from_string('''\
import hashlib
import streamlit as st
from mn_molecular_design_studio import ligand_tools_ui as ui, vina_inputs

target, compounds = vina_inputs.stage_example("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
st.session_state["vina-target-path"] = str(target)
st.session_state["vina-compound-input"] = compounds
st.session_state["vina-default-target-signature"] = hashlib.sha1(target.read_bytes()).hexdigest()[:12]
st.session_state["ligand-vina-pocket-mode"] = "Select Residues"

def submit(schema, workflow, parameters, project, settings, title):
    st.session_state["test-vina-submission"] = {"workflow": workflow, "parameters": parameters}

ui._render_vina_setup(
    "ligand.docking.vina", {}, owner_id="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
    compound_sets=[], project_selector=lambda key: "project-1", submit=submit,
    title="AutoDock Vina",
)
''').run()
    assert not app.exception
    assert app.session_state["ligand-vina-pocket-mode"] == "Select Residues"
    old_center = app.number_input(key="vina-box-center-x").value
    app.text_input(key="vina-residue-text-A").set_value("225-226").run()
    assert not app.exception
    assert ("A", 225) in app.session_state["ligand-vina-selected-residues"]
    assert app.number_input(key="vina-box-center-x").value == pytest.approx(old_center)
    assert not app.session_state.get("ligand-vina-applied-residues")
    assert app.button(key="ligand-vina-run").disabled

    old_focus_epoch = app.session_state.get("vina-viewer-focus-epoch", 0)
    app.button(key="vina-apply-selected-residues").click().run()
    assert not app.exception
    assert app.session_state["ligand-vina-applied-residues"] == [("A", 225), ("A", 226)]
    assert not app.button(key="ligand-vina-run").disabled
    assert app.session_state["vina-viewer-focus-epoch"] == old_focus_epoch + 1
    assert app.number_input(key="vina-box-center-x").value != pytest.approx(old_center)

    app.number_input(key="vina-box-size-x").set_value(42.0).run()
    assert app.number_input(key="vina-box-size-x").value == pytest.approx(42.0)
    app.number_input(key="ligand-vina-residue-padding").set_value(7.0).run()
    target, _ = vina_inputs.stage_example("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
    coordinates = ligand_tools_ui._vina_protein_atom_coordinates(target.read_text())
    selected_points = [
        point for residue in (225, 226) for point in coordinates[("A", residue)]
    ]
    _, expected_size = ligand_tools_ui._vina_box_from_points(selected_points, padding=7.0)
    assert app.number_input(key="vina-box-size-x").value == pytest.approx(expected_size["x"])
    assert app.session_state["vina-applied-residue-padding"] == pytest.approx(7.0)


@pytest.mark.parametrize(
    ("workflow_id", "title", "tool_key", "expected_options"),
    [
        ("ligand.docking.autodock_gpu", "AutoDock-GPU", "autodock-gpu", {"nrun", "gpu_device", "job_name"}),
        ("ligand.docking.gnina", "GNINA", "gnina", {"exhaustiveness", "gpu_device", "cnn_scoring", "job_name"}),
        ("ligand.docking.unidock_pro", "Uni-Dock Pro", "unidock-pro", {"search_mode", "gpu_device", "job_name"}),
        ("ligand.docking.rosetta_ligand", "RosettaLigand", "rosetta", {"protocol", "reference_mode", "padding", "job_name"}),
    ],
)
def test_shared_docking_page_submits_engine_specific_options(
    tmp_path, monkeypatch, workflow_id, title, tool_key, expected_options,
):
    from streamlit.testing.v1 import AppTest

    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path))
    app = AppTest.from_string(f'''\
import streamlit as st
from mn_molecular_design_studio import ligand_tools_ui as ui

def submit(schema, workflow, parameters, project, settings, title):
    st.session_state["test-docking-submission"] = {{"workflow": workflow, "parameters": parameters}}

ui._render_vina_setup(
    {workflow_id!r}, {{}}, owner_id="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
    compound_sets=[], project_selector=lambda key: "project-1", submit=submit,
    title={title!r},
)
''').run()
    assert not app.exception
    assert any(title in str(element.value) for element in app.markdown)
    app.button(key="vina-load-example").click().run()
    assert not app.exception
    assert any(
        f"{title} then docks against the imported receptor." in str(element.value)
        for element in app.caption
    )
    if workflow_id == "ligand.docking.rosetta_ligand":
        app.selectbox(key="ligand-rosetta-reference-mode").set_value("reference_guided").run()
        assert not app.exception
        assert app.button(key=f"ligand-{tool_key}-run").disabled
        app.selectbox(key="ligand-rosetta-reference-mode").set_value("pocket_center").run()
    app.button(key=f"ligand-{tool_key}-run").click().run()
    assert not app.exception
    submitted = app.session_state["test-docking-submission"]
    assert submitted["workflow"] == workflow_id
    assert expected_options <= submitted["parameters"].keys()
    expected_job_name = {
        "autodock-gpu": "autodock-gpu-example",
        "gnina": "gnina-example", "unidock-pro": "unidock-pro-example",
        "rosetta": "rosettaligand-example",
    }[tool_key]
    assert submitted["parameters"]["job_name"] == expected_job_name
    assert Path(submitted["parameters"]["target_pdb"]).is_file()
    assert Path(submitted["parameters"]["compound_set_file"]).is_file()
    ligand_tools.validate_ligand_tool_request(workflow_id, submitted["parameters"], require_files=False)
    if workflow_id == "ligand.docking.gnina":
        assert submitted["parameters"]["cnn_scoring"] == "rescore"


def test_nesso_setup_submits_vina_example_target_and_database_compound_set(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest

    from mn_molecular_design_studio import vina_inputs

    monkeypatch.setenv("MN_MOLECULAR_DESIGN_STUDIO_APP_HOME", str(tmp_path))
    # The form must be exercisable without the large Nesso model/reference bundle;
    # actual readiness and inference are outside this UI smoke test.
    monkeypatch.setattr(
        ligand_tools,
        "nesso_readiness",
        lambda: {
            "ready": True,
            "checkpoint_ready": True,
            "ccd_ready": True,
            "esm_ready": True,
            "checkpoint_dir": str(tmp_path / "nesso" / "checkpoint"),
            "ccd_path": str(tmp_path / "nesso" / "components.cif"),
            "esm_cache_dir": str(tmp_path / "nesso" / "esm"),
        },
    )
    target, example_set = vina_inputs.stage_example(
        "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", project_id="project-1",
    )
    app = AppTest.from_string(f'''\
import streamlit as st
from mn_molecular_design_studio import ligand_tools_ui as ui

def submit(schema, workflow, parameters, project, settings, title):
    st.session_state["test-nesso-submission"] = {{"workflow": workflow, "parameters": parameters}}

ui.render_setup(
    "ligand.affinity.nesso", {{}}, owner_id="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
    compound_sets=[{{"group_id": "example-imatinib", "path": {example_set["path"]!r},
        "name": "Imatinib example", "compound_count": 1}}],
    project_selector=lambda key: "project-1", submit=submit,
)
''').run()
    assert not app.exception
    app.file_uploader(key="ligand-tool-target-ligand.affinity.nesso").set_value(
        ("1iep.pdb", target.read_bytes(), "chemical/x-pdb"),
    ).run()
    assert not app.exception
    run_button = app.button(key="ligand-tool-run-ligand.affinity.nesso")
    assert not run_button.disabled
    run_button.click().run()
    assert not app.exception

    submitted = app.session_state["test-nesso-submission"]
    parameters = submitted["parameters"]
    assert submitted["workflow"] == "ligand.affinity.nesso"
    assert Path(parameters["target_pdb"]).read_bytes() == target.read_bytes()
    assert parameters["compound_set_file"] == example_set["path"]
    assert parameters["recycling_steps"] == 5
    assert parameters["seed"] == 42
    assert parameters["compound_limit"] == 0
    ligand_tools.validate_ligand_tool_request(
        "ligand.affinity.nesso", parameters, require_files=False,
    )


def test_nesso_reuses_shared_proteina_complexa_esm_cache(tmp_path, monkeypatch):
    reference_root = tmp_path / "reference_files"
    # MN-Ligand keeps ESM-2 in the shared Proteina-Complexa cache. The
    # Nesso-specific Hugging Face directories may exist without the model.
    (reference_root / "nesso" / ".cache" / "huggingface").mkdir(parents=True)
    snapshot = (
        reference_root / "proteina-complexa" / "hf-cache"
        / "models--facebook--esm2_t33_650M_UR50D" / "snapshots" / "esm-revision"
    )
    snapshot.mkdir(parents=True)
    for filename in ("model.safetensors", "config.json", "vocab.txt"):
        (snapshot / filename).write_text("present")
    checkpoint = reference_root / "nesso" / "v1.0.0"
    checkpoint.mkdir(parents=True)
    for filename in ("model.safetensors", "hparams.json"):
        (checkpoint / filename).write_text("present")
    (reference_root / "nesso" / "ccd.pkl").write_text("present")
    monkeypatch.setattr(ligand_tools, "reference_root", lambda: reference_root)
    monkeypatch.delenv("MN_NESSO_ESM_CACHE_DIR", raising=False)

    checkpoint_path, ccd_path, esm_cache = ligand_tools.nesso_reference_paths()

    assert esm_cache == reference_root / "proteina-complexa" / "hf-cache"
    assert checkpoint_path == checkpoint
    assert ccd_path == reference_root / "nesso" / "ccd.pkl"
    assert ligand_tools.nesso_readiness()["ready"] is True


def test_nesso_mounts_huggingface_cache_at_cli_expected_location(tmp_path, monkeypatch):
    target = tmp_path / "target.pdb"
    target.write_text(
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\nEND\n"
    )
    compound_set = tmp_path / "set.json"
    ligand_compounds.write_compound_set(
        compound_set,
        name="Smoke set",
        compounds=[{"compound_id": "ethanol", "smiles": "CCO"}],
        owner_id="owner",
        project_id="project",
        created_at="now",
    )
    reference_root = tmp_path / "reference_files"
    checkpoint = reference_root / "nesso" / "v1.0.0"
    checkpoint.mkdir(parents=True)
    for filename in ("model.safetensors", "hparams.json"):
        (checkpoint / filename).write_text("present")
    (reference_root / "nesso" / "ccd.pkl").write_text("present")
    cache = reference_root / "proteina-complexa" / "hf-cache"
    snapshot = cache / "models--facebook--esm2_t33_650M_UR50D" / "snapshots" / "revision"
    snapshot.mkdir(parents=True)
    for filename in ("model.safetensors", "config.json", "vocab.txt"):
        (snapshot / filename).write_text("present")
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path / "runs")
    monkeypatch.setattr(ligand_tools, "reference_root", lambda: reference_root)
    monkeypatch.setattr(ligand_tools, "load_manifest", lambda _name: {"image": "mn-nesso:test"})

    def fake_docker(command, **_kwargs):
        work_mount = Path(command[command.index("-v") + 1].split(":", 1)[0])
        prediction = work_mount / "output" / "predictions" / "ethanol" / "affinity.json"
        prediction.parent.mkdir(parents=True)
        prediction.write_text(json.dumps({"affinity_pred_value": -5.2, "affinity_probability_binary": 0.8}))
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(ligand_tools.subprocess, "run", fake_docker)

    run_dir = ligand_tools.run_ligand_nesso(target, compound_set, recycling_steps=1, compound_limit=1)

    command = jobs.read_json(run_dir / "command.json")["command"]
    result = jobs.read_json(run_dir / "result.json")
    assert command[command.index("--cache") + 1] == "/cache"
    assert f"{cache}:/cache/huggingface:ro" in command
    assert result["success"] is True
    assert result["metrics"]["scored_count"] == 1


def test_ligand_backend_readiness_uses_studio_bundle_even_without_external_checkout(monkeypatch):
    monkeypatch.setenv("MN_LIGAND_ROOT", "/path/that/does/not/exist")
    bundled_root = ligand_tools._studio_ligand_runtime_root()

    assert bundled_root.is_relative_to(Path(ligand_tools.__file__).resolve().parents[1])
    assert ligand_tools.ligand_complex_preparation_readiness()["ready"] is True
    assert ligand_tools.gromacs_md_runner_readiness()["ready"] is True
    assert ligand_tools.openmm_md_runner_readiness()["ready"] is True
    assert ligand_tools.rosettaligand_runner_readiness()["ready"] is True


def test_parse_compound_table_maps_smiles_and_preserves_properties():
    rows = ligand_compounds.parse_compound_table(
        b"catalog id,canonical_smiles,activity\nA-01,CCO,active\nA-02,c1ccccc1,inactive\n",
        filename="library.csv",
    )

    assert [row["compound_id"] for row in rows] == ["A-01", "A-02"]
    assert [row["smiles"] for row in rows] == ["CCO", "c1ccccc1"]
    assert rows[0]["properties"] == {"activity": "active"}


def test_generated_compound_csv_uses_canonical_isomeric_smiles():
    rows = ligand_compounds.parse_compound_table(
        b"compound_id,canonical_isomeric_smiles,generation_engine\nlddm-1,C[C@H](O)C,lddm\n",
        filename="generated.csv",
    )

    assert rows[0]["compound_id"] == "lddm-1"
    assert rows[0]["smiles"] == "C[C@H](O)C"
    assert rows[0]["properties"] == {"generation_engine": "lddm"}


def test_compound_import_keeps_distinct_id_name_and_smiles():
    rows = ligand_compounds.parse_compound_table(
        b"compound_id,compound_name,canonical_smiles\n70701262,Bepristat_1a,CCN\n",
        filename="compounds.csv",
    )
    assert rows[0]["compound_id"] == "70701262"
    assert rows[0]["compound_name"] == "Bepristat_1a"
    assert rows[0]["smiles"] == "CCN"


def test_smi_parser_and_compound_set_round_trip(tmp_path):
    compounds = ligand_compounds.parse_smi_text("CCO ethanol\nc1ccccc1 benzene\n")
    path = tmp_path / "set.json"
    written = ligand_compounds.write_compound_set(
        path,
        name="Examples",
        compounds=compounds,
        owner_id="owner-1",
        project_id="project-1",
        created_at="2026-10-19T10:00:00+00:00",
    )

    assert written["record_type"] == "compound_set"
    assert written["compound_count"] == 2
    loaded = ligand_compounds.load_compound_set(path)
    assert loaded["compounds"] == compounds
    assert loaded["project_id"] == "project-1"


def test_compound_set_rejects_duplicate_ids_and_missing_smiles():
    with pytest.raises(ValueError, match="occurs more than once"):
        ligand_compounds.normalize_compounds([
            {"compound_id": "same", "smiles": "CCO"},
            {"compound_id": "same", "smiles": "CCN"},
        ])
    with pytest.raises(ValueError, match="no SMILES"):
        ligand_compounds.normalize_compounds([{"compound_id": "empty", "smiles": ""}])


def test_database_compound_set_requires_known_schema(tmp_path):
    path = tmp_path / "wrong.json"
    path.write_text(json.dumps({"record_type": "compound_set", "schema_version": 99, "compounds": []}))

    with pytest.raises(ValueError, match="Unsupported compound-set schema"):
        ligand_compounds.load_compound_set(path)


def test_database_compound_set_detects_content_changes(tmp_path):
    path = tmp_path / "set.json"
    ligand_compounds.write_compound_set(
        path,
        name="Integrity",
        compounds=[{"compound_id": "one", "smiles": "CCO"}],
        owner_id="owner",
        project_id="project",
        created_at="now",
    )
    payload = json.loads(path.read_text())
    payload["compounds"][0]["smiles"] = "CCN"
    path.write_text(json.dumps(payload))

    with pytest.raises(ValueError, match="content-integrity"):
        ligand_compounds.load_compound_set(path)


def test_legacy_compound_sets_migrate_to_named_schema_without_losing_integrity(tmp_path):
    path = tmp_path / "legacy.json"
    old_compounds = [{"compound_id": "ethanol", "smiles": "CCO", "properties": {}}]
    digest = hashlib.sha256(json.dumps(old_compounds, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    path.write_text(json.dumps({"schema_version": 1, "record_type": "compound_set", "compound_set_id": "old",
        "content_sha256": digest, "compounds": old_compounds}))
    loaded = ligand_compounds.load_compound_set(path)
    assert loaded["schema_version"] == 2
    assert loaded["compounds"][0]["compound_name"] == "ethanol"


def test_nesso_yaml_uses_protein_ligand_affinity_contract():
    payload = yaml.safe_load(ligand_tools._nesso_yaml([("A", "MPEPTIDE")], "CC(=O)O"))

    assert payload["sequences"][0]["protein"]["sequence"] == "MPEPTIDE"
    assert payload["sequences"][1]["ligand"] == {"id": "L", "smiles": "CC(=O)O"}
    assert payload["properties"] == [{"affinity": {"binder": "L"}}]


def test_vina_request_validates_signed_pocket_centers(tmp_path):
    receptor = tmp_path / "target.pdb"
    receptor.write_text(
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\nEND\n"
    )
    set_path = tmp_path / "compounds.json"
    ligand_compounds.write_compound_set(
        set_path,
        name="One",
        compounds=[{"compound_id": "one", "smiles": "CCO"}],
        owner_id="owner",
        project_id="project",
        created_at="now",
    )
    ligand_tools.validate_ligand_tool_request(
        "ligand.docking.vina",
        {
            "target_pdb": receptor,
            "compound_set_file": set_path,
            "compound_limit": 10,
            "center_x": -5.0,
            "center_y": 1.0,
            "center_z": 2.0,
            "size_x": 20.0,
            "size_y": 20.0,
            "size_z": 20.0,
            "exhaustiveness": 8,
            "num_modes": 9,
            "seed": 42,
        },
    )


def test_vina_request_rejects_excessive_docking_box(tmp_path):
    with pytest.raises(ValueError, match="box dimensions"):
        ligand_tools.validate_ligand_tool_request(
            "ligand.docking.vina",
            {
                "compound_limit": 1,
                "center_x": 0.0,
                "center_y": 0.0,
                "center_z": 0.0,
                "size_x": 61.0,
                "size_y": 20.0,
                "size_z": 20.0,
                "exhaustiveness": 8,
                "num_modes": 9,
                "seed": 42,
            },
            require_files=False,
        )


def _write_axis_ligand(path: Path) -> Path:
    atom_rows = []
    for x, y, z, element in ((0.0, 0.0, 0.0, "C"), (0.0, 2.0, 0.0, "C"), (0.0, 4.0, 0.0, "O")):
        atom_rows.append(
            f"{x:10.4f}{y:10.4f}{z:10.4f} {element:>2}  0  0  0  0  0  0  0  0  0  0  0  0"
        )
    path.write_text(
        "axis ligand\n  test\n\n  3  0  0  0  0  0            999 V2000\n"
        + "\n".join(atom_rows)
        + "\nM  END\n$$$$\n"
    )
    return path


def test_mn_ligand_axis_transform_rotates_target_and_box_to_global_x(tmp_path):
    ligand = _write_axis_ligand(tmp_path / "axis.sdf")
    target = (
        "ATOM      1  CA  ALA A   1       0.000   3.000   0.000  1.00 20.00           C\n"
        "END\n"
    )
    transform = ligand_tools.ligand_longest_axis_transform(ligand)
    rotated = ligand_tools.transform_pdb_data(target, transform)
    center, size = ligand_tools.transform_axis_aligned_box(
        (0.0, 2.0, 0.0), (2.0, 6.0, 4.0), transform
    )

    assert transform["aligned_longest_axis"] == [1.0, 0.0, 0.0]
    assert transform["source_longest_axis"] == pytest.approx([0.0, 1.0, 0.0])
    assert rotated.splitlines()[0][30:54] == f"{1.000:8.3f}{2.000:8.3f}{0.000:8.3f}"
    assert center == pytest.approx((0.0, 2.0, 0.0))
    assert size == pytest.approx((6.0, 2.0, 4.0))


def test_axis_alignment_can_preserve_box_dimensions_in_aligned_frame(tmp_path):
    ligand = _write_axis_ligand(tmp_path / "axis.sdf")
    target = tmp_path / "target.pdb"
    target.write_text(
        "ATOM      1  CA  ALA A   1       0.000   3.000   0.000  1.00 20.00           C\nEND\n"
    )

    _rotated, center, size, transform = ligand_tools._docking_axis_alignment(
        target, (0.0, 2.0, 0.0), (30.0, 20.0, 20.0),
        align_longest_axis_to_x=True, alignment_ligand_sdf=ligand,
        preserve_box_size_after_alignment=True,
    )

    assert transform["aligned_longest_axis"] == [1.0, 0.0, 0.0]
    assert center == pytest.approx((0.0, 2.0, 0.0))
    assert size == (30.0, 20.0, 20.0)


def test_vina_axis_alignment_records_transform_and_stages_rotated_receptor(tmp_path, monkeypatch):
    target = tmp_path / "target.pdb"
    target.write_text(
        "ATOM      1  CA  ALA A   1       0.000   3.000   0.000  1.00 20.00           C\nEND\n"
    )
    ligand = _write_axis_ligand(tmp_path / "axis.sdf")
    compound_set = tmp_path / "set.json"
    ligand_compounds.write_compound_set(
        compound_set,
        name="Aligned",
        compounds=[{"compound_id": "ethanol", "smiles": "CCO"}],
        owner_id="owner",
        project_id="project",
        created_at="now",
    )
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path / "runs")
    monkeypatch.setattr(ligand_tools, "load_manifest", lambda _name: {"image": "mn-docking-suite:test"})

    def fake_docker(run_dir, _manifest, _args, *, gpu_device):
        (run_dir / "artifacts" / "docking_scores.csv").write_text(
            "compound_id\tstatus\tbest_score_kcal_mol\tpose_path\terror\n"
            "ethanol\tcompleted\t-5.4\tartifacts/poses/ethanol.pdbqt\t\n"
        )
        return 0

    monkeypatch.setattr(ligand_tools, "run_docker_job", fake_docker)
    run_dir = ligand_tools.run_ligand_vina(
        target,
        compound_set,
        center_x=0.0,
        center_y=2.0,
        center_z=0.0,
        size_x=2.0,
        size_y=6.0,
        size_z=4.0,
        align_longest_axis_to_x=True,
        alignment_ligand_sdf=ligand,
    )

    contract = json.loads((run_dir / "artifact_contract.json").read_text())
    params = jobs.read_json(run_dir / "input.json")["params"]
    receptor = (run_dir / "artifacts" / "input" / "receptor.pdb").read_text()
    assert receptor.splitlines()[0][30:54] == f"{1.000:8.3f}{2.000:8.3f}{0.000:8.3f}"
    assert params["original_box_size"] == [2.0, 6.0, 4.0]
    assert params["box_size"] == pytest.approx([6.0, 2.0, 4.0])
    assert (run_dir / "artifacts" / "input" / "alignment_transform.json").is_file()
    assert any(item["type"] == "coordinate_transform/v1" for item in contract["outputs"])
    assert any(item["name"] == "axis_alignment_ligand" for item in contract["inputs"])


def test_complex_preparation_tile_is_registered_and_validates_trim_and_extension(tmp_path, monkeypatch):
    workflow_id = "ligand.complex.prepare"
    schema = workflow_automation.workflow_schema(workflow_id)
    assert workflow_id in {row["workflow"] for row in workflow_automation.workflow_catalog()}
    assert {
        "source_pdb", "repair", "max_internal_gap", "extension_chain",
        "extension_sequence", "trim_ranges_json",
    } <= {field["name"] for field in schema["parameters"]}

    source = tmp_path / "complex.pdb"
    source.write_text(
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\n"
        "HETATM    2  C1  LIG A 101       4.000   5.000   6.000  1.00 20.00           C\nEND\n"
    )
    monkeypatch.setattr(ligand_tools, "ligand_complex_preparation_readiness", lambda: {"ready": True, "source": "fixture"})
    ligand_tools.validate_ligand_tool_request(
        workflow_id,
        {
            "source_pdb": source,
            "repair": True,
            "max_internal_gap": 15,
            "extension_chain": "A",
            "extension_sequence": "GS",
            "trim_ranges_json": json.dumps({"A": [1, 3]}),
        },
    )
    with pytest.raises(ValueError, match="canonical amino-acid"):
        ligand_tools.validate_ligand_tool_request(
            workflow_id,
            {"source_pdb": source, "extension_chain": "A", "extension_sequence": "GX"},
        )


def test_complex_preparation_tile_publishes_typed_complex_and_report(tmp_path, monkeypatch):
    source = tmp_path / "complex.pdb"
    source.write_text(
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\n"
        "HETATM    2  C1  LIG A 101       4.000   5.000   6.000  1.00 20.00           C\nEND\n"
    )
    source_root = tmp_path / "studio-ligand-runtime"
    for relative in (
        "mn_ligand/workflows/protein_preparation.py",
        "mn_ligand/workflows/terminal_repair.py",
        "mn_ligand/workflows/target_trimming.py",
        "mn_ligand/workflows/modeller_terminal_extension_runner.py",
    ):
        path = source_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("# test readiness stub\n")
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path / "runs")
    monkeypatch.setattr(ligand_tools, "_studio_ligand_runtime_root", lambda: source_root)

    def fake_prepare(command, **kwargs):
        output_dir = Path(command[3])
        (output_dir / "prepared_complex.pdb").write_text(source.read_text())
        (output_dir / "prepared_target.pdb").write_text(source.read_text().split("HETATM")[0] + "TER\nEND\n")
        (output_dir / "repair_report.json").write_text(json.dumps({"source_import": {"ligands": [{"resname": "LIG"}]}}))
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(ligand_tools.subprocess, "run", fake_prepare)
    run_dir = ligand_tools.run_ligand_complex_prepare(source)
    result = jobs.read_json(run_dir / "result.json")
    outputs = result["outputs"]["artifacts"]
    assert result["success"] is True
    assert {row["name"] for row in outputs} == {"prepared_complex", "prepared_target", "repair_report"}
    assert "HETATM" in (run_dir / "artifacts" / "prepared" / "prepared_complex.pdb").read_text()


def test_vina_job_stages_exact_compound_subset_with_typed_pipeline_contract(tmp_path, monkeypatch):
    target = tmp_path / "target.pdb"
    target.write_text(
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\nEND\n"
    )
    compound_set = tmp_path / "database-set.json"
    database_record = ligand_compounds.write_compound_set(
        compound_set,
        name="Screening set",
        compounds=[
            {"compound_id": "ethanol", "smiles": "CCO"},
            {"compound_id": "ethylamine", "smiles": "CCN"},
        ],
        owner_id="owner",
        project_id="project",
        created_at="now",
        compound_set_id="database-set-1",
    )
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path / "runs")
    monkeypatch.setattr(ligand_tools, "load_manifest", lambda _name: {"image": "mn-docking-suite:test"})

    def fake_docker(run_dir, _manifest, args, *, gpu_device):
        assert gpu_device == "none"
        docker_script = args[-1]
        assert "--log" not in docker_script
        assert '> "artifacts/logs/${compound_id}.vina.log" 2>&1' in docker_script
        assert 'artifacts/logs/${compound_id}.vina.log' in docker_script
        assert 'scrub.py "artifacts/ligands/${compound_id}_source.sdf"' in docker_script
        assert "--cpu 1 --ph 7.40 --skip_tautomer" in docker_script
        assert 'mk_prepare_ligand.py -i "artifacts/ligands/${compound_id}.sdf"' in docker_script
        shell_syntax = ligand_tools.subprocess.run(
            ["bash", "-n", "-c", docker_script], capture_output=True, text=True, check=False,
        )
        assert shell_syntax.returncode == 0, shell_syntax.stderr
        scores = run_dir / "artifacts" / "docking_scores.csv"
        scores.write_text(
            "compound_id\tstatus\tbest_score_kcal_mol\tpose_path\terror\n"
            "ethanol\tcompleted\t-5.4\tartifacts/poses/ethanol.pdbqt\t\n"
        )
        return 0

    monkeypatch.setattr(ligand_tools, "run_docker_job", fake_docker)
    run_dir = ligand_tools.run_ligand_vina(
        target,
        compound_set,
        center_x=-2.0,
        center_y=1.0,
        center_z=3.0,
        compound_limit=1,
        job_name="EGFR ligand screen",
    )

    contract = json.loads((run_dir / "artifact_contract.json").read_text())
    staged_set = ligand_compounds.load_compound_set(run_dir / "artifacts" / "input" / "compound_set.json")
    assert contract["compound_set_id"] == database_record["compound_set_id"] == "database-set-1"
    assert contract["compound_set_sha256"] == database_record["content_sha256"]
    assert contract["compound_count"] == 1
    assert contract["inputs"][0]["path"] == "artifacts/input/receptor.pdb"
    assert contract["inputs"][1]["type"] == "compound_set/v1"
    assert contract["inputs"][1]["consumed_count"] == 1
    assert contract["outputs"][0] == {
        "name": "docking_scores",
        "type": "ligand_docking_scores/v1",
        "path": "artifacts/docking_scores.csv",
    }
    assert [row["compound_id"] for row in staged_set["compounds"]] == ["ethanol"]
    assert staged_set["source_compound_set_id"] == "database-set-1"
    assert staged_set["source_content_sha256"] == database_record["content_sha256"]
    assert jobs.read_json(run_dir / "input.json")["inputs"]["compound_set_file"] == str(compound_set.resolve())
    params = jobs.read_json(run_dir / "input.json")["params"]
    assert params["use_scrub"] is True
    assert params["scrub_ph"] == pytest.approx(7.4)
    assert params["scrub_skip_tautomer"] is True
    assert params["exhaustiveness"] == 32
    assert params["num_modes"] == 10
    assert jobs.read_json(run_dir / "metadata.json")["campaign_name"] == "EGFR ligand screen"


def test_vina_docking_fans_out_compounds_with_bounded_cpu_workers(tmp_path, monkeypatch):
    target = tmp_path / "target.pdb"
    target.write_text(
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\nEND\n"
    )
    compound_set = tmp_path / "set.json"
    ligand_compounds.write_compound_set(
        compound_set,
        name="Parallel Vina set",
        compounds=[
            {"compound_id": "ethanol", "smiles": "CCO"},
            {"compound_id": "ethylamine", "smiles": "CCN"},
            {"compound_id": "methanol", "smiles": "CO"},
        ],
        owner_id="owner",
        project_id="project",
        created_at="now",
    )
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path / "runs")
    monkeypatch.setattr(ligand_tools, "load_manifest", lambda _name: {"image": "mn-docking-suite:test"})
    captured = {}

    def fake_docker(run_dir, _manifest, args, *, gpu_device):
        captured["script"] = args[-1]
        assert gpu_device == "none"
        script_check = ligand_tools.subprocess.run(
            ["bash", "-n", "-c", args[-1]], capture_output=True, text=True, check=False,
        )
        assert script_check.returncode == 0, script_check.stderr
        (run_dir / "artifacts" / "docking_scores.csv").write_text(
            "compound_id\tstatus\tbest_score_kcal_mol\tcnn_score\tcnn_affinity\tpose_path\terror\n"
            "ethanol\tcompleted\t-5.4\t\t\tartifacts/poses/ethanol.pdbqt\t\n"
            "ethylamine\tcompleted\t-4.1\t\t\tartifacts/poses/ethylamine.pdbqt\t\n"
            "methanol\tcompleted\t-3.0\t\t\tartifacts/poses/methanol.pdbqt\t\n"
        )
        return 0

    monkeypatch.setattr(ligand_tools, "run_docker_job", fake_docker)
    run_dir = ligand_tools.run_ligand_vina(
        target, compound_set,
        center_x=0.0, center_y=0.0, center_z=0.0,
        compound_limit=3, cpu_cores=2,
    )

    script = captured["script"]
    assert "DOCKING_COMPOUND_WORKERS=2" in script
    assert 'xargs -r -P "$DOCKING_COMPOUND_WORKERS" -n 1 bash -c' in script
    assert "--cpu 1" in script
    assert '"artifacts/logs/${compound_id}.result.tsv.tmp"' in script
    params = jobs.read_json(run_dir / "input.json")["params"]
    assert params["cpu_cores"] == 2
    assert params["compound_workers"] == 2
    assert params["cpu_threads_per_compound"] == 1


def test_gnina_tile_records_cnn_scores_and_uses_gpu_container(tmp_path, monkeypatch):
    target = tmp_path / "target.pdb"
    target.write_text(
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\nEND\n"
    )
    compound_set = tmp_path / "set.json"
    ligand_compounds.write_compound_set(
        compound_set,
        name="GNINA set",
        compounds=[
            {"compound_id": "ethanol", "smiles": "CCO"},
            {"compound_id": "ethylamine", "smiles": "CCN"},
            {"compound_id": "methanol", "smiles": "CO"},
        ],
        owner_id="owner",
        project_id="project",
        created_at="now",
    )
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path / "runs")
    monkeypatch.setattr(ligand_tools, "load_manifest", lambda name: {"image": "mn-docking-suite:test", "tool": name})
    captured = {}

    def fake_docker(run_dir, manifest, args, *, gpu_device):
        captured.update(manifest=manifest, args=args, gpu_device=gpu_device)
        script_check = ligand_tools.subprocess.run(
            ["bash", "-n", "-c", args[-1]], capture_output=True, text=True, check=False,
        )
        assert script_check.returncode == 0, script_check.stderr
        (run_dir / "artifacts" / "docking_scores.csv").write_text(
            "compound_id\tstatus\tbest_score_kcal_mol\tcnn_score\tcnn_affinity\tpose_path\terror\n"
            "ethanol\tcompleted\t-6.2\t0.81\t-7.4\tartifacts/poses/ethanol.pdbqt\t\n"
        )
        return 0

    monkeypatch.setattr(ligand_tools, "run_docker_job", fake_docker)
    run_dir = ligand_tools.run_ligand_gnina(
        target, compound_set, center_x=0.0, center_y=1.0, center_z=2.0, gpu_device="1",
        cnn_scoring="refinement", job_name="EGFR GNINA",
        scrub_ph=6.2, scrub_skip_tautomer=False, cpu_cores=4, compound_workers=2,
    )

    script = captured["args"][-1]
    result = jobs.read_json(run_dir / "result.json")
    contract = json.loads((run_dir / "artifact_contract.json").read_text())
    assert captured["manifest"]["tool"] == "gnina"
    assert captured["gpu_device"] == "1"
    assert "gnina --receptor" in script and "--cnn_scoring refinement" in script
    assert "DOCKING_COMPOUND_WORKERS=2" in script
    assert 'xargs -r -P "$DOCKING_COMPOUND_WORKERS" -n 1 bash -c' in script
    assert "--cpu 2" in script
    assert "--ph 6.20" in script and "--skip_tautomer" not in script
    assert result["metrics"]["primary_ranking_metric"] == "cnn_score"
    assert result["metrics"]["best_cnn_score"] == pytest.approx(0.81)
    params = jobs.read_json(run_dir / "input.json")["params"]
    assert params["compound_workers"] == 2
    assert params["cpu_threads_per_compound"] == 2
    assert jobs.read_json(run_dir / "metadata.json")["campaign_name"] == "EGFR GNINA"
    assert contract["workflow"] == "ligand.docking.gnina"
    assert contract["outputs"][0]["type"] == "ligand_docking_scores/v1"


def test_unidock_pro_tile_batches_compounds_on_gpu_and_emits_typed_poses(tmp_path, monkeypatch):
    target = tmp_path / "target.pdb"
    target.write_text(
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\nEND\n"
    )
    compound_set = tmp_path / "set.json"
    ligand_compounds.write_compound_set(
        compound_set,
        name="Uni-Dock set",
        compounds=[{"compound_id": "ethanol", "smiles": "CCO"}, {"compound_id": "ethylamine", "smiles": "CCN"}],
        owner_id="owner",
        project_id="project",
        created_at="now",
    )
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path / "runs")
    monkeypatch.setattr(ligand_tools, "load_manifest", lambda name: {"image": "mn-docking-suite:test", "tool": name, "gpu": True})
    captured = {}

    def fake_docker(run_dir, manifest, args, *, gpu_device, use_host_user=False):
        captured.update(manifest=manifest, args=args, gpu_device=gpu_device, use_host_user=use_host_user)
        (run_dir / "artifacts" / "docking_scores.csv").write_text(
            "compound_id\tstatus\tbest_score_kcal_mol\tcnn_score\tcnn_affinity\tpose_path\terror\n"
            "ethanol\tcompleted\t-7.1\t\t\tartifacts/poses/ethanol_out.pdbqt\t\n"
        )
        return 0

    monkeypatch.setattr(ligand_tools, "run_docker_job", fake_docker)
    run_dir = ligand_tools.run_ligand_unidock_pro(
        target, compound_set, center_x=-2.0, center_y=1.0, center_z=3.0,
        search_mode="detail", gpu_device="1", cpu_cores=4,
    )

    script = captured["args"][-1]
    shell_syntax = ligand_tools.subprocess.run(["bash", "-n", "-c", script], capture_output=True, text=True, check=False)
    preparation_script = (run_dir / "artifacts" / "input" / "prepare_one_ligand.sh").read_text()
    batch_preparation_script = (run_dir / "artifacts" / "input" / "prepare_batch.sh").read_text()
    preparation_syntax = ligand_tools.subprocess.run(
        ["bash", "-n", str(run_dir / "artifacts" / "input" / "prepare_one_ligand.sh")],
        capture_output=True, text=True, check=False,
    )
    batch_preparation_syntax = ligand_tools.subprocess.run(
        ["bash", "-n", str(run_dir / "artifacts" / "input" / "prepare_batch.sh")],
        capture_output=True, text=True, check=False,
    )
    contract = json.loads((run_dir / "artifact_contract.json").read_text())
    result = jobs.read_json(run_dir / "result.json")
    assert captured["manifest"]["tool"] == "unidock_pro"
    assert captured["gpu_device"] == "1"
    assert captured["use_host_user"] is True
    assert "--ligand_index artifacts/input/ligand_index.txt" in script
    assert "--search_mode detail" in script
    assert "xargs -0 -r -n 2 -P 2 bash artifacts/input/prepare_one_ligand.sh" in script
    assert "find artifacts/input/preparation_batches -type f -name '*.smi' -print0" in script
    assert "xargs -0 -r -n 1 -P 2 bash artifacts/input/prepare_batch.sh" in script
    assert 'scrub.py "$batch_smi" -o "$batch_sdf" --cpu 1 --ph 7.40 --skip_tautomer' in batch_preparation_script
    assert 'mk_prepare_ligand.py -i "$batch_sdf" --multimol_outdir artifacts/ligands' in batch_preparation_script
    assert "Prepared %s of %s compounds using %s CPU workers" in script
    assert "--cpu 1 --ph 7.40 --skip_tautomer" in preparation_script
    assert 'scrub.py "$source_sdf"' in preparation_script
    assert (run_dir / "artifacts" / "input" / "compounds.nul").read_bytes() == b"ethanol\0CCO\0ethylamine\0CCN\0"
    assert (run_dir / "artifacts" / "input" / "preparation_batches" / "batch_00001.smi").read_text() == "CCO ethanol\n"
    assert (run_dir / "artifacts" / "input" / "preparation_batches" / "batch_00002.smi").read_text() == "CCN ethylamine\n"
    assert '\\"$compound_id\\"' not in script
    assert shell_syntax.returncode == 0, shell_syntax.stderr
    assert preparation_syntax.returncode == 0, preparation_syntax.stderr
    assert batch_preparation_syntax.returncode == 0, batch_preparation_syntax.stderr
    assert contract["workflow"] == "ligand.docking.unidock_pro"
    assert contract["outputs"][1]["type"] == "ligand_pose_collection/v1"
    assert result["metrics"]["best_score_kcal_mol"] == pytest.approx(-7.1)
    assert result["metrics"]["preparation_cpu_workers"] == 2
    assert result["metrics"]["preparation_batch_size"] == 1
    assert result["metrics"]["preparation_batch_count"] == 2


def test_rosettaligand_tile_uses_bundled_runner_and_preserves_reu_labels(tmp_path, monkeypatch):
    target = tmp_path / "target.pdb"
    target.write_text(
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\nEND\n"
    )
    compound_set = tmp_path / "set.json"
    ligand_compounds.write_compound_set(
        compound_set,
        name="Rosetta set",
        compounds=[{"compound_id": "ethanol", "smiles": "CCO"}],
        owner_id="owner",
        project_id="project",
        created_at="now",
    )
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path / "runs")
    monkeypatch.setattr(ligand_tools, "rosettaligand_runner_readiness", lambda: {"ready": True, "source": "/studio-ligand-runtime/mn_ligand/workflows/openvs.py", "image": "mn-openvs:test"})
    monkeypatch.setattr(ligand_tools, "_openvs_runner_constants", lambda: {
        "_RUNNER_SCRIPT": "#!/bin/bash\ncd /workspace\n",
        "_CENTER_PDB_SCRIPT": "print('center')\n",
        "_ROSETTA_PARAMS_VALIDATOR_SCRIPT": "print('validate')\n",
    })
    captured = {}

    def fake_docker(run_dir, manifest, args, *, gpu_device):
        captured.update(manifest=manifest, args=args, gpu_device=gpu_device)
        native = run_dir / "native" / "replicate_001" / "chunks" / "run.score.sc"
        native.parent.mkdir(parents=True)
        native.write_text(
            "SCORE: total_score dG ligandname lig_rms description\n"
            "SCORE: -10.2 -7.4 ethanol 1.1 pose0001\n"
        )
        pose = run_dir / "poses" / "pose0001.pdb"
        pose.parent.mkdir()
        pose.write_text("ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\nEND\n")
        return 0

    monkeypatch.setattr(ligand_tools, "run_docker_job", fake_docker)
    run_dir = ligand_tools.run_ligand_rosetta_ligand(
        target, compound_set, center_x=2.0, center_y=-1.0, center_z=3.0,
        protocol="vsh", cpu_cores=2,
    )

    script = captured["args"][-1]
    rows = list(__import__("csv").DictReader((run_dir / "artifacts" / "docking_scores.csv").open()))
    contract = json.loads((run_dir / "artifact_contract.json").read_text())
    result = jobs.read_json(run_dir / "result.json")
    assert captured["gpu_device"] == "none"
    assert captured["manifest"]["tool"] == "rosetta_ligand"
    assert "OPENVS_PROTOCOL=VSH" in script
    assert "run_openvs.sh" in script
    assert rows[0]["estimated_dg_reu"] == "-7.4"
    assert rows[0]["pose_path"].endswith("pose_01.pdb")
    assert contract["workflow"] == "ligand.docking.rosetta_ligand"
    assert contract["outputs"][0]["score_units"] == "REU"
    assert result["metrics"]["score_units"].startswith("Rosetta energy units")


@pytest.mark.parametrize(
    ("workflow_id", "required_parameters"),
    [
        ("ligand.generation.lddm", {"target_pdb", "reference_ligand_sdf", "requested_count", "sampler"}),
        ("ligand.generation.pocketxmol", {"target_pdb", "pocket_pdb", "requested_count", "diffusion_steps"}),
    ],
)
def test_generation_tiles_are_registered_with_typed_outputs(workflow_id, required_parameters):
    schema = workflow_automation.workflow_schema(workflow_id)
    assert workflow_id in {row["workflow"] for row in workflow_automation.workflow_catalog()}
    assert required_parameters <= {field["name"] for field in schema["parameters"]}


@pytest.mark.parametrize(
    ("engine", "required_parameters"),
    [
        ("omtra", {"target_pdb", "requested_count", "engine_settings", "pharmacophore_json"}),
        ("flowr_root", {"target_pdb", "pocket_pdb", "reference_ligand_sdf", "engine_settings"}),
        ("conditar", {"target_pdb", "pocket_pdb", "engine_settings"}),
        ("paopt", {"target_pdb", "pocket_pdb", "reference_ligand_sdf", "engine_settings"}),
        ("drugrpg", {"target_pdb", "pocket_pdb", "engine_settings"}),
        ("pfm", {"target_pdb", "pocket_pdb", "engine_settings"}),
        ("pocketflow", {"target_pdb", "pocket_pdb", "engine_settings"}),
        ("pgmg", {"target_pdb", "pharmacophore_posp", "engine_settings"}),
    ],
)
def test_additional_generation_tiles_are_registered(engine, required_parameters):
    workflow_id = f"ligand.generation.{engine}"
    schema = workflow_automation.workflow_schema(workflow_id)
    assert workflow_id in {row["workflow"] for row in workflow_automation.workflow_catalog()}
    assert required_parameters <= {field["name"] for field in schema["parameters"]}


@pytest.mark.parametrize("engine", ["omtra", "flowr_root", "conditar", "paopt", "drugrpg", "pfm", "pocketflow", "pgmg"])
def test_additional_generation_tiles_emit_database_ready_compound_sets(engine, tmp_path, monkeypatch):
    target = tmp_path / "target.pdb"
    target.write_text("ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\nEND\n")
    pocket = tmp_path / "pocket.pdb"
    pocket.write_text(target.read_text())
    reference_ligand = tmp_path / "reference.sdf"
    reference_ligand.write_text("reference\n  test\n\n  1  0  0  0  0  0            999 V2000\nM  END\n$$$$\n")
    pharmacophore = tmp_path / "pharmacophore.posp"
    pharmacophore.write_text("AROM 0.0 1.0 2.0\n")
    pharmit = tmp_path / "pharmit.json"
    pharmit.write_text(json.dumps({"points": [{"name": "Aromatic", "x": 0.0, "y": 1.0, "z": 2.0, "radius": 1.0}]}))
    reference_root = tmp_path / "references"
    for relative in ligand_tools._GENERATOR_CONFIGS[engine]["references"]:
        path = reference_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("test")
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path / "runs")
    monkeypatch.setattr(ligand_tools, "reference_root", lambda: reference_root)
    monkeypatch.setattr(ligand_tools, "load_manifest", lambda name: {"image": f"mn-{name}:test"})
    captured = {}

    def fake_docker(command, **kwargs):
        captured.update(command=command, kwargs=kwargs)
        run_dir = Path(command[command.index("-v") + 1].split(":/work", 1)[0])
        normalized = run_dir / "normalized"
        normalized.mkdir(parents=True, exist_ok=True)
        (normalized / "generated_compounds.csv").write_text(
            "compound_id,canonical_isomeric_smiles,generation_engine\n"
            f"{engine}-generated-0000001,C[C@H](O)C,{engine}\n"
        )
        (normalized / "generated_compounds.sdf").write_text("generated SDF\n")
        (normalized / "generation_report.json").write_text(json.dumps({"success": True, "unique_valid_compound_count": 1}))
        return SimpleNamespace(returncode=0, stdout="generated", stderr="")

    monkeypatch.setattr(ligand_tools.subprocess, "run", fake_docker)
    inputs = {"omtra": {"pharmacophore_json": pharmit}, "flowr_root": {"pocket_pdb": pocket, "reference_ligand_sdf": reference_ligand},
              "conditar": {"pocket_pdb": pocket}, "paopt": {"pocket_pdb": pocket, "reference_ligand_sdf": reference_ligand, "engine_settings": {"optimize_properties": ["HIA"]}},
              "drugrpg": {"pocket_pdb": pocket}, "pfm": {"pocket_pdb": pocket},
              "pocketflow": {"pocket_pdb": pocket}, "pgmg": {"pharmacophore_posp": pharmacophore}}[engine]
    runner = getattr(ligand_tools, f"run_ligand_{engine}")
    run_dir = runner(target_pdb=target, requested_count=1, **inputs)

    contract = json.loads((run_dir / "artifact_contract.json").read_text())
    generated = ligand_compounds.load_compound_set(run_dir / "artifacts" / "generated_compound_set.json")
    result = jobs.read_json(run_dir / "result.json")
    command = captured["command"]
    assert captured["kwargs"]["check"] is False
    assert f"mn-{engine}:test" in command
    if engine == "omtra":
        assert "--pharmacophore" in command
    assert contract["workflow"] == f"ligand.generation.{engine}"
    assert contract["compound_count"] == 1
    assert generated["compounds"][0]["smiles"] == "C[C@H](O)C"
    assert generated["compounds"][0]["properties"]["generation_engine"] == engine
    assert result["success"] is True


@pytest.mark.parametrize("engine", ["lddm", "pocketxmol"])
def test_generation_tile_emits_database_ready_compound_set(engine, tmp_path, monkeypatch):
    target = tmp_path / "target.pdb"
    target.write_text(
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\nEND\n"
    )
    pocket = tmp_path / "pocket.pdb"
    pocket.write_text(target.read_text())
    reference_ligand = tmp_path / "reference.sdf"
    reference_ligand.write_text("reference\n  test\n\n  1  0  0  0  0  0            999 V2000\nM  END\n$$$$\n")
    reference_root = tmp_path / "references"
    for relative in (ligand_tools._LDDM_CHECKPOINT, ligand_tools._POCKETXMOL_CHECKPOINT, ligand_tools._POCKETXMOL_CONFIG):
        path = reference_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("test")
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path / "runs")
    monkeypatch.setattr(ligand_tools, "reference_root", lambda: reference_root)
    monkeypatch.setattr(ligand_tools, "load_manifest", lambda name: {"image": f"mn-{name}:test"})
    captured = {}

    def fake_docker(command, **kwargs):
        captured.update(command=command, kwargs=kwargs)
        run_dir = Path(command[command.index("-v") + 1].split(":/work", 1)[0])
        normalized = run_dir / "normalized"
        normalized.mkdir(parents=True, exist_ok=True)
        (normalized / "generated_compounds.csv").write_text(
            "compound_id,canonical_isomeric_smiles,generation_engine,native_source,native_index,coordinate_dimension\n"
            f"{engine}-generated-0000001,C[C@H](O)C,{engine},samples.sdf,1,3\n"
        )
        (normalized / "generated_compounds.sdf").write_text("generated SDF\n")
        (normalized / "generation_report.json").write_text(json.dumps({"success": True, "unique_valid_compound_count": 1}))
        return SimpleNamespace(returncode=0, stdout="generated", stderr="")

    monkeypatch.setattr(ligand_tools.subprocess, "run", fake_docker)
    if engine == "lddm":
        run_dir = ligand_tools.run_ligand_lddm(target, reference_ligand, requested_count=1)
        expected_workflow = "ligand.generation.lddm"
    else:
        run_dir = ligand_tools.run_ligand_pocketxmol(target, pocket, requested_count=1)
        expected_workflow = "ligand.generation.pocketxmol"

    contract = json.loads((run_dir / "artifact_contract.json").read_text())
    generated = ligand_compounds.load_compound_set(run_dir / "artifacts" / "generated_compound_set.json")
    result = jobs.read_json(run_dir / "result.json")
    assert captured["kwargs"]["check"] is False
    assert contract["workflow"] == expected_workflow
    assert contract["compound_count"] == 1
    assert contract["outputs"][1]["type"] == "compound_set/v1"
    assert generated["compounds"][0]["smiles"] == "C[C@H](O)C"
    assert generated["compounds"][0]["properties"]["generation_engine"] == engine
    assert result["success"] is True


@pytest.mark.parametrize(
    ("workflow_id", "required_parameters"),
    [
        ("ligand.docking.vina", {"target_pdb", "compound_set_file", "center_x", "center_y", "center_z"}),
        ("ligand.docking.autodock_gpu", {"target_pdb", "compound_set_file", "center_x", "center_y", "center_z", "nrun", "gpu_device"}),
        ("ligand.docking.gnina", {"target_pdb", "compound_set_file", "center_x", "center_y", "center_z", "gpu_device"}),
        ("ligand.docking.unidock_pro", {"target_pdb", "compound_set_file", "center_x", "center_y", "center_z", "search_mode", "gpu_device"}),
        ("ligand.docking.rosetta_ligand", {"target_pdb", "compound_set_file", "center_x", "center_y", "center_z", "protocol", "reference_mode", "cpu_cores"}),
        ("ligand.affinity.nesso", {"target_pdb", "compound_set_file", "recycling_steps", "seed"}),
    ],
)
def test_single_tool_tiles_are_registered_with_database_inputs(workflow_id, required_parameters):
    schema = workflow_automation.workflow_schema(workflow_id)
    assert workflow_id in {row["workflow"] for row in workflow_automation.workflow_catalog()}
    assert required_parameters <= {field["name"] for field in schema["parameters"]}
    if workflow_id in {
        "ligand.docking.vina", "ligand.docking.gnina",
        "ligand.docking.unidock_pro", "ligand.docking.autodock_gpu",
    }:
        assert {"use_scrub", "scrub_ph", "scrub_skip_tautomer"} <= {
            field["name"] for field in schema["parameters"]
        }


def test_autodock_gpu_validation_requires_device_and_bounded_runs():
    request = {
        "compound_limit": 1,
        "center_x": 0.0,
        "center_y": 0.0,
        "center_z": 0.0,
        "size_x": 20.0,
        "size_y": 20.0,
        "size_z": 20.0,
        "nrun": 20,
        "seed": 42,
        "gpu_device": "0",
    }
    ligand_tools.validate_ligand_tool_request("ligand.docking.autodock_gpu", request, require_files=False)
    with pytest.raises(ValueError, match="requires a CUDA GPU"):
        ligand_tools.validate_ligand_tool_request(
            "ligand.docking.autodock_gpu", {**request, "gpu_device": "none"}, require_files=False,
        )
    with pytest.raises(ValueError, match="LGA runs"):
        ligand_tools.validate_ligand_tool_request(
            "ligand.docking.autodock_gpu", {**request, "nrun": 0}, require_files=False,
        )


def test_meeko_docking_prep_settings_validate_ph_and_boolean_values():
    request = {
        "center_x": 0.0, "center_y": 0.0, "center_z": 0.0,
        "size_x": 20.0, "size_y": 20.0, "size_z": 20.0,
        "exhaustiveness": 8, "num_modes": 9, "seed": 42,
        "use_scrub": True, "scrub_ph": 7.4, "scrub_skip_tautomer": True,
    }
    ligand_tools.validate_ligand_tool_request(
        "ligand.docking.vina", request, require_files=False,
    )
    with pytest.raises(ValueError, match="Scrub pH"):
        ligand_tools.validate_ligand_tool_request(
            "ligand.docking.vina", {**request, "scrub_ph": 14.1}, require_files=False,
        )
    with pytest.raises(ValueError, match="must be booleans"):
        ligand_tools.validate_ligand_tool_request(
            "ligand.docking.vina", {**request, "use_scrub": "yes"}, require_files=False,
        )
def test_autodock_gpu_tile_stages_ad4_grid_workflow_and_typed_pose_outputs(tmp_path, monkeypatch):
    target = tmp_path / "target.pdb"
    target.write_text(
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\nEND\n"
    )
    compound_set = tmp_path / "set.json"
    ligand_compounds.write_compound_set(
        compound_set,
        name="AutoDock-GPU set",
        compounds=[{"compound_id": "ethanol", "smiles": "CCO"}],
        owner_id="owner",
        project_id="project",
        created_at="now",
    )
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path / "runs")
    monkeypatch.setattr(ligand_tools, "load_manifest", lambda name: {"image": "mn-docking-suite:test", "gpu": True, "tool": name})
    captured = {}

    def fake_docker(run_dir, manifest, args, *, gpu_device):
        captured.update(manifest=manifest, args=args, gpu_device=gpu_device)
        (run_dir / "artifacts" / "docking_scores.csv").write_text(
            "compound_id,status,best_score_kcal_mol,pose_path,error\n"
            "ethanol,completed,-8.2,artifacts/poses/ethanol.pdbqt,\n"
        )
        jobs.write_json(run_dir / "artifacts" / "autodock_gpu_result.json", {"success": True, "scored_count": 1})
        return 0

    monkeypatch.setattr(ligand_tools, "run_docker_job", fake_docker)
    run_dir = ligand_tools.run_ligand_autodock_gpu(
        target, compound_set, center_x=-1.0, center_y=2.0, center_z=3.0,
        nrun=12, gpu_device="1", compound_limit=1,
    )

    contract = json.loads((run_dir / "artifact_contract.json").read_text())
    config = jobs.read_json(run_dir / "artifacts" / "input" / "autodock_gpu_config.json")
    runner = (run_dir / "artifacts" / "input" / "autodock_gpu_runner.py").read_text()
    result = jobs.read_json(run_dir / "result.json")
    assert captured["manifest"]["tool"] == "autodock_gpu"
    assert captured["gpu_device"] == "1"
    assert captured["args"] == ["python", "/work/artifacts/input/autodock_gpu_runner.py"]
    assert config["center"] == [-1.0, 2.0, 3.0]
    assert "autogrid4" in runner and "--ffile" in runner and "--gbest" in runner
    assert contract["workflow"] == "ligand.docking.autodock_gpu"
    assert contract["outputs"][0]["score_units"].startswith("kcal/mol")
    assert contract["outputs"][1]["type"] == "ligand_pose_collection/v1"
    assert result["metrics"]["best_score_kcal_mol"] == pytest.approx(-8.2)
    assert result["metrics"]["primary_ranking_direction"] == "lower"


def test_autodock_gpu_container_runner_builds_maps_and_parses_ad4_energy(tmp_path, monkeypatch):
    input_root = tmp_path / "artifacts" / "input"
    input_root.mkdir(parents=True)
    (input_root / "autodock_gpu_config.json").write_text(json.dumps({
        "center": [-1.0, 2.0, 3.0], "box_size": [20.0, 20.0, 20.0], "nrun": 5, "seed": 42,
    }))
    (input_root / "receptor.pdb").write_text("ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\nEND\n")
    (input_root / "compounds.tsv").write_text("compound_id\tsmiles\nethanol\tCCO\n")
    (input_root / "smiles_to_3d.py").write_text("# mock generator")

    commands = []

    def fake_run(command, log_path, *, cwd):
        command = list(map(str, command))
        commands.append(command)
        if command[0] == "mk_prepare_receptor.py":
            receptor_base = Path(command[command.index("-o") + 1])
            receptor_base.with_suffix(".gpf").write_text("parameter_file AD4.2_bound.dat\n")
        elif command[0] == "autogrid4":
            (cwd / "receptor.maps.fld").write_text("grid descriptor")
        elif "smiles_to_3d.py" in command[1]:
            Path(command[command.index("--output") + 1]).write_text("ligand SDF")
        elif command[0] == "mk_prepare_ligand.py":
            Path(command[command.index("-o") + 1]).write_text("ligand PDBQT")
        elif command[0] == "adgpu":
            (cwd / "adgpu_ethanol.dlg").write_text(
                "DOCKED: USER Estimated Free Energy of Binding = -8.20 kcal/mol\n"
                "DOCKED: USER Estimated Free Energy of Binding = -7.10 kcal/mol\n"
            )
            (cwd / "adgpu_ethanol_best.pdbqt").write_text("best pose")
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text("mocked")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(autodock_gpu_runner, "_run", fake_run)
    result = autodock_gpu_runner.run_tile(tmp_path)

    score_rows = list(__import__("csv").DictReader((tmp_path / "artifacts" / "docking_scores.csv").open()))
    assert result["success"] is True
    assert score_rows[0]["status"] == "completed"
    assert float(score_rows[0]["best_score_kcal_mol"]) == pytest.approx(-8.2)
    assert (tmp_path / score_rows[0]["pose_path"]).read_text() == "best pose"
    scrub_command = next(command for command in commands if command[0] == "scrub.py")
    assert scrub_command[scrub_command.index("--ph") + 1] == "7.40"
    assert "--skip_tautomer" in scrub_command


def _write_gromacs_inputs(tmp_path: Path) -> tuple[Path, Path]:
    complex_path = tmp_path / "complex.pdb"
    complex_path.write_text(
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\n"
        "HETATM    2  C1  LIG A 501       4.000   5.000   6.000  1.00 20.00           C\n"
        "HETATM    3  C2  LIG A 501       5.000   5.000   6.000  1.00 20.00           C\n"
        "HETATM    4  O1  LIG A 501       6.000   5.000   6.000  1.00 20.00           O\n"
        "END\n"
    )
    sdf_path = tmp_path / "ligand.sdf"
    sdf_path.write_text(
        "bound-ligand\n  Studio\n\n  3  2  0  0  0  0            999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    2.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0\n  2  3  1  0\nM  END\n$$$$\n"
    )
    return complex_path, sdf_path


def test_gromacs_md_tile_is_registered_with_typed_complex_and_trajectory_ports():
    workflow_id = "ligand.md.gromacs"
    schema = workflow_automation.workflow_schema(workflow_id)
    row = next(row for row in workflow_automation.workflow_catalog() if row["workflow"] == workflow_id)

    assert row["category"] == "ligand"
    assert {"complex_pdb", "ligand_sdf", "ligand_key", "production_ns"} <= {
        field["name"] for field in schema["parameters"]
    }


def test_openmm_md_tile_is_registered_with_typed_complex_and_trajectory_ports():
    workflow_id = "ligand.md.openmm"
    schema = workflow_automation.workflow_schema(workflow_id)
    assert workflow_id in {row["workflow"] for row in workflow_automation.workflow_catalog()}
    assert {"complex_pdb", "ligand_sdf", "ligand_key", "production_ns", "timestep_fs"} <= {
        field["name"] for field in schema["parameters"]
    }


def test_openmm_md_tile_uses_studio_bundled_runtime_and_records_trajectory(tmp_path, monkeypatch):
    complex_path, sdf_path = _write_gromacs_inputs(tmp_path)
    runtime_root = tmp_path / "studio-ligand-runtime"
    (runtime_root / "mn_ligand" / "workflows").mkdir(parents=True)
    (runtime_root / "mn_ligand" / "ligandx" / "services" / "md").mkdir(parents=True)
    (runtime_root / "mn_ligand" / "workflows" / "bound_ligand_md.py").write_text("# source")
    (runtime_root / "mn_ligand" / "ligandx" / "services" / "md" / "service.py").write_text("# source")
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path / "runs")
    monkeypatch.setattr(ligand_tools, "_studio_ligand_runtime_root", lambda: runtime_root)
    monkeypatch.setattr(ligand_tools, "load_manifest", lambda _name: {"image": "mn-md:test", "gpu": True})
    captured = {}

    def fake_docker(run_dir, manifest, args, *, gpu_device, extra_mounts, environment):
        captured.update(manifest=manifest, args=args, gpu_device=gpu_device, extra_mounts=extra_mounts, environment=environment)
        md_root = run_dir / "artifacts" / "md"
        trajectory = md_root / "production.dcd"
        trajectory.write_text("trajectory")
        (md_root / "result.json").write_text(json.dumps({
            "success": True,
            "md_result": {"status": "success", "output_files": {"production_trajectory": "/work/artifacts/md/production.dcd"}},
        }))
        return 0

    monkeypatch.setattr(ligand_tools, "run_docker_job", fake_docker)
    run_dir = ligand_tools.run_ligand_openmm_md(
        complex_path, sdf_path, "LIG|A|501|_", production_ns=0.1, cpu_cores=2,
    )
    contract = json.loads((run_dir / "artifact_contract.json").read_text())
    result = jobs.read_json(run_dir / "result.json")
    config = json.loads((run_dir / "artifacts" / "input" / "openmm_config.json").read_text())

    assert captured["manifest"]["image"] == "mn-md:test"
    assert captured["gpu_device"] == "0"
    assert captured["extra_mounts"] == [(runtime_root, "/studio-ligand-runtime", True)]
    assert captured["environment"] == {"PYTHONPATH": "/studio-ligand-runtime"}
    assert config["production_steps"] == 50_000
    assert config["mmgbsa_enabled"] is False
    assert config["preparation_protocol"] == "current_staged"
    assert contract["workflow"] == "ligand.md.openmm"
    assert contract["outputs"][0]["type"] == "md_trajectory_set/v1"
    assert result["success"] is True


def test_gromacs_md_validates_bound_residue_against_sdf_heavy_atoms(tmp_path, monkeypatch):
    complex_path, sdf_path = _write_gromacs_inputs(tmp_path)
    monkeypatch.setattr(ligand_tools, "gromacs_md_runner_readiness", lambda: {"ready": True})
    ligands = ligand_tools.bound_ligand_choices(complex_path)

    assert ligands[0]["key"] == "LIG|A|501|_"
    assert ligands[0]["heavy_atom_count"] == 3
    assert ligands[0]["heavy_atom_coordinates"] == {
        "C1": (4.0, 5.0, 6.0),
        "C2": (5.0, 5.0, 6.0),
        "O1": (6.0, 5.0, 6.0),
    }
    ligand_tools.validate_ligand_tool_request(
        "ligand.md.gromacs",
        {
            "complex_pdb": complex_path,
            "ligand_sdf": sdf_path,
            "ligand_key": "LIG|A|501|_",
            "production_ns": 0.1,
            "timestep_fs": 2.0,
            "seed": 42,
            "cpu_cores": 2,
        },
    )
    sdf_path.write_text(sdf_path.read_text().replace("  3  2", "  2  1", 1))
    with pytest.raises(ValueError, match="does not match the selected PDB residue"):
        ligand_tools.validate_ligand_tool_request(
            "ligand.md.gromacs",
            {"complex_pdb": complex_path, "ligand_sdf": sdf_path, "ligand_key": "LIG|A|501|_"},
        )
    _, sdf_path = _write_gromacs_inputs(tmp_path)
    sdf_lines = sdf_path.read_text().splitlines()
    sdf_lines[5] = sdf_lines[5][:31] + "O  " + sdf_lines[5][34:]
    sdf_lines[6] = sdf_lines[6][:31] + "C  " + sdf_lines[6][34:]
    sdf_path.write_text("\n".join(sdf_lines) + "\n")
    with pytest.raises(ValueError, match="heavy-atom element order"):
        ligand_tools.validate_ligand_tool_request(
            "ligand.md.gromacs",
            {"complex_pdb": complex_path, "ligand_sdf": sdf_path, "ligand_key": "LIG|A|501|_"},
        )


def test_gromacs_md_tile_stages_typed_artifacts_and_mounts_studio_runtime_read_only(tmp_path, monkeypatch):
    complex_path, sdf_path = _write_gromacs_inputs(tmp_path)
    runtime_root = tmp_path / "studio-ligand-runtime"
    (runtime_root / "mn_ligand" / "workflows").mkdir(parents=True)
    (runtime_root / "mn_ligand" / "workflows" / "gromacs_md.py").write_text("# source")
    (runtime_root / "mn_ligand" / "workflows" / "bound_ligand_md.py").write_text("# source")
    monkeypatch.setattr(jobs, "runs_root", lambda: tmp_path / "runs")
    monkeypatch.setattr(ligand_tools, "_studio_ligand_runtime_root", lambda: runtime_root)
    monkeypatch.setattr(ligand_tools, "load_manifest", lambda _name: {"image": "mn-gromacs:test", "gpu": True})
    captured = {}

    def fake_docker(run_dir, manifest, args, *, gpu_device, extra_mounts, environment):
        captured.update(manifest=manifest, args=args, gpu_device=gpu_device, extra_mounts=extra_mounts, environment=environment)
        production = run_dir / "artifacts" / "md" / "production"
        production.mkdir(parents=True)
        (production / "production_whole.xtc").write_bytes(b"trajectory")
        (production / "trajectory_analysis.json").write_text(json.dumps({"ligand_rmsd": 1.2}))
        jobs.write_json(run_dir / "artifacts" / "md" / "result.json", {
            "success": True,
            "output_files": {
                "production_trajectory": "/work/artifacts/md/production/production_whole.xtc",
                "analysis_report": "/work/artifacts/md/production/trajectory_analysis.json",
            },
        })
        return 0

    monkeypatch.setattr(ligand_tools, "run_docker_job", fake_docker)
    run_dir = ligand_tools.run_ligand_gromacs_md(
        complex_path, sdf_path, "LIG|A|501|_", production_ns=0.1, cpu_cores=2, use_gpu=False,
    )

    contract = json.loads((run_dir / "artifact_contract.json").read_text())
    config = json.loads((run_dir / "artifacts" / "input" / "gromacs_config.json").read_text())
    runner = (run_dir / "artifacts" / "input" / "gromacs_tile_runner.py").read_text()
    assert captured["manifest"]["image"] == "mn-gromacs:test"
    assert captured["gpu_device"] == "none"
    assert captured["extra_mounts"] == [(runtime_root, "/studio-ligand-runtime", True)]
    assert captured["environment"] == {"PYTHONPATH": "/studio-ligand-runtime"}
    assert config["production_steps"] == 50_000
    assert "prepare_gromacs_system" in runner and "run_gromacs_production" in runner
    assert contract["inputs"][0]["type"] == "protein_ligand_complex/pdb/v1"
    assert contract["outputs"][0]["type"] == "md_trajectory_set/v1"
    assert jobs.read_json(run_dir / "result.json")["success"] is True


def test_docker_runner_supports_read_only_source_mount_and_environment(tmp_path, monkeypatch):
    source = tmp_path / "source"
    source.mkdir()
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    monkeypatch.setattr("mn_molecular_design_studio.engine.workflow_engine.core.docker_runner.reference_root", lambda: tmp_path / "missing-ref")

    command = build_docker_command(
        {"image": "test-image", "gpu": False},
        run_dir,
        ["python", "run.py"],
        extra_mounts=[(source, "/studio-ligand-runtime", True)],
        environment={"PYTHONPATH": "/studio-ligand-runtime"},
    )

    assert f"{source.resolve()}:/studio-ligand-runtime:ro" in command
    assert command[command.index("-e") + 1] == "PYTHONPATH=/studio-ligand-runtime"


def test_docker_runner_can_write_outputs_as_host_user(tmp_path, monkeypatch):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    docker_runner = "mn_molecular_design_studio.engine.workflow_engine.core.docker_runner"
    monkeypatch.setattr(f"{docker_runner}.reference_root", lambda: tmp_path / "missing-ref")
    monkeypatch.setattr(f"{docker_runner}._docker_is_rootless", lambda: False)
    monkeypatch.setattr(f"{docker_runner}.os.getuid", lambda: 1234)
    monkeypatch.setattr(f"{docker_runner}.os.getgid", lambda: 2345)

    command = build_docker_command(
        {"image": "test-image", "gpu": True}, run_dir, ["bash", "-lc", "touch output"],
        gpu_device="0", use_host_user=True,
    )

    assert command[command.index("--user") + 1] == "1234:2345"
    assert "HOME=/tmp" in command


def test_docker_runner_skips_explicit_user_for_rootless_docker(tmp_path, monkeypatch):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    docker_runner = "mn_molecular_design_studio.engine.workflow_engine.core.docker_runner"
    monkeypatch.setattr(f"{docker_runner}.reference_root", lambda: tmp_path / "missing-ref")
    monkeypatch.setattr(f"{docker_runner}._docker_is_rootless", lambda: True)

    command = build_docker_command(
        {"image": "test-image", "gpu": False}, run_dir, ["touch", "output"], use_host_user=True,
    )

    assert "--user" not in command
    assert not any(value.startswith("HOME=") for value in command)


def test_nesso_api_validation_checks_compound_set_and_records_affinity_only(tmp_path, monkeypatch):
    target = tmp_path / "target.pdb"
    target.write_text(
        "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C\nEND\n"
    )
    compound_set = tmp_path / "set.json"
    ligand_compounds.write_compound_set(
        compound_set,
        name="Screen",
        compounds=[{"compound_id": "ethyl-acetate", "smiles": "CCOC(C)=O"}],
        owner_id="owner",
        project_id="project",
        created_at="now",
    )
    monkeypatch.setattr(ligand_tools, "nesso_readiness", lambda: {"ready": True})
    request_path = tmp_path / "request.json"
    request_path.write_text(json.dumps({
        "schema_version": 1,
        "workflow": "ligand.affinity.nesso",
        "parameters": {
            "target_pdb": str(target),
            "compound_set_file": str(compound_set),
            "recycling_steps": 5,
            "seed": 42,
            "compound_limit": 10,
            "gpu_device": "0",
        },
    }))

    validated = workflow_automation.validate_workflow_request(request_path)

    assert validated["valid"] is True
    assert validated["workflow"] == "ligand.affinity.nesso"
    assert validated["parameters"]["compound_set_file"] == str(compound_set.resolve())
