from __future__ import annotations

import json
import hashlib
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import yaml

from mn_molecular_design_studio.cofolding_input import default_model, model_to_payload
from mn_molecular_design_studio.engine.workflow_engine.workflows import af3_jax_cofolding, optimized_prediction as workflow


class BoltzInputContractTests(unittest.TestCase):
    def test_boltz2_tile_request_builder_uses_shared_optimized_workflow_contract(self):
        parameters = workflow.boltz2_prediction_parameters(
            [{"name": "probe_1", "compound_id": "probe-001", "compound_name": "Ethanol probe",
              "protein_sequences": "ACDEFGHIKLMNPQRSTVWY", "ligand_smiles": "CCO"}],
            settings={
                "optimization_mode": "off", "seed": 7, "samples": 2, "num_recycles": 3,
                "cpu_cores": 4, "output_format": "pdb", "max_msa_seqs": "",
                "use_msa": False, "use_msa_server": False, "generate_ipsae": False,
            },
            job_name="boltz-pocket",
        )
        self.assertEqual(parameters["engine"], "boltz2")
        self.assertEqual(parameters["job_name"], "boltz-pocket")
        self.assertEqual(parameters["engine_options"]["output_format"], "pdb")
        self.assertEqual(parameters["engine_options"]["max_msa_seqs"], 8192)
        self.assertFalse(parameters["engine_options"]["use_msa_server"])
        workflow.validate_prediction_request(parameters)

    def test_copies_keep_modifications_and_custom_msa_on_the_selected_chain(self):
        model = default_model()
        model["cards"] = [
            {"kind": "protein", "value": "ac d", "copies": 2,
             "cyclic": True, "modifications": [{"position": 2, "code": "MSE"}]},
            {"kind": "ligand", "mode": "smiles", "value": "CCO", "copies": 1},
        ]
        model["compound_name"] = "probe ethanol"
        model["compound_id"] = "probe-001"
        model["msa_files"] = [{"chain_id": "B", "name": "B.a3m", "text": ">B\nACD\n"}]
        values = model_to_payload(model)
        self.assertEqual(values["input_error"], "")
        self.assertEqual(values["compound_id"], "probe-001")
        self.assertEqual(values["compound_name"], "probe ethanol")
        self.assertEqual(values["protein_sequences"], "ACD:ACD")
        self.assertEqual([row["id"] for row in values["chain_catalog"]], ["A", "B", "C"])
        self.assertEqual(values["molecule_chain_options"][0]["msa_a3m_text"], "")
        self.assertEqual(values["molecule_chain_options"][1]["msa_a3m_text"], ">B\nACD\n")
        self.assertTrue(all(row["cyclic"] for row in values["molecule_chain_options"]))

    def test_an_empty_added_card_is_not_silently_discarded(self):
        model = default_model()
        model["cards"][0]["value"] = "ACD"
        model["cards"].append({"kind": "protein", "value": "", "copies": 1})
        self.assertEqual(model_to_payload(model)["input_error"], "All sequences must be non empty.")

    def test_ligand_complex_requires_a_name_and_defaults_id_to_name(self):
        model = default_model()
        model["cards"] = [
            {"kind": "protein", "value": "ACD", "copies": 1},
            {"kind": "ligand", "mode": "smiles", "value": "CCO", "copies": 1},
        ]
        missing = model_to_payload(model)
        self.assertEqual(missing["input_error"], "Enter a compound name for this ligand complex.")
        model["compound_name"] = "Ethanol"
        named = model_to_payload(model)
        self.assertEqual(named["compound_id"], "Ethanol")
        self.assertEqual(named["name"], "Ethanol")

    def test_standalone_af3_family_rejects_unnamed_ligand_before_creating_job(self):
        with patch.object(af3_jax_cofolding, "create_job") as create_job:
            with self.assertRaisesRegex(ValueError, "compound name"):
                af3_jax_cofolding.run_af3_jax_cofolding(
                    protein_sequences="ACDEFG",
                    ligand_smiles="CCO",
                    model_name="openbind0",
                )
        create_job.assert_not_called()

    def test_templates_single_sequence_mode_and_affinity_reach_native_yaml(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "input.yaml"
            path.write_text(yaml.safe_dump({"version": 1, "sequences": [
                {"protein": {"id": "A", "sequence": "ACD", "msa": "old.a3m"}},
                {"ligand": {"id": "B", "smiles": "CCO"}},
            ]}))
            details = workflow._apply_boltz_input_options(
                root, path, item_id="sample", use_msa=False, affinity_binder="B",
                templates=[{"name": "../template.pdb", "text": "HEADER TEMPLATE\nEND"}],
                template_distance_threshold=2.5,
            )
            payload = yaml.safe_load(path.read_text())
            self.assertEqual(payload["sequences"][0]["protein"]["msa"], "empty")
            template = payload["templates"][0]
            self.assertTrue(template["pdb"].startswith("/work/artifacts/raw/optimized_prediction/templates/"))
            self.assertTrue(template["force"])
            self.assertEqual(template["threshold"], 2.5)
            self.assertEqual(payload["properties"], [{"affinity": {"binder": "B"}}])
            self.assertTrue(details["force"])

    def test_native_yaml_template_file_references_are_rewritten(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            path, details = workflow._stage_uploaded_boltz_yaml(
                root, item_id="sample", yaml_text="version: 1\nsequences:\n- protein: {id: A, sequence: ACD}\ntemplates:\n- cif: template.cif\n",
                auxiliary_files=[{"name": "template.cif", "text": "data_template\n"}],
            )
            payload = yaml.safe_load(path.read_text())
            self.assertTrue(payload["templates"][0]["cif"].startswith("/work/"))
            self.assertEqual(details["referenced_files"], ["template.cif"])

    def test_native_input_staging_applies_the_builder_options(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            path, details = workflow._stage_native_input(
                root, profile=workflow.ENGINE_PROFILES["boltz2"], item_id="sample", job_name="example",
                proteins=[{"id": "A", "sequence": "ACD"}], dna=[], rna=[],
                ligand_ccd_codes=[], ligand_smiles=[{"id": "B", "smiles": "CCO"}],
                seed=1, samples=1, target_pdb_text="", target_chain_id="A", target_msa_a3m="",
                binder_sequence="", num_recycles=3, use_msa=False, affinity_binder="B",
                boltz_templates=[{"name": "template.cif", "text": "data_template\n"}],
                template_distance_threshold=3,
            )
            native = yaml.safe_load(path.read_text())
            self.assertEqual(native["sequences"][0]["protein"]["msa"], "empty")
            self.assertEqual(native["properties"][0]["affinity"]["binder"], "B")
            self.assertEqual(native["templates"][0]["threshold"], 3)
            self.assertTrue(details["force"])

    def test_reference_a3m_is_attached_and_only_cache_misses_use_colabfold(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference_files"
            repository = reference / "boltz_models" / "msa_repository"
            repository.mkdir(parents=True)
            cached = repository / f"{hashlib.sha256(b'ACD').hexdigest()}.a3m"
            cached.write_text(">query\nACD\n>hit\nACD\n", encoding="utf-8")
            with patch.object(workflow, "reference_root", return_value=reference):
                path, details = workflow._stage_native_input(
                    root, profile=workflow.ENGINE_PROFILES["boltz2"], item_id="sample", job_name="example",
                    proteins=[{"id": "A", "sequence": "ACD"}, {"id": "B", "sequence": "EFG"}],
                    dna=[], rna=[], ligand_ccd_codes=[], ligand_smiles=[], seed=1, samples=1,
                    target_pdb_text="", target_chain_id="A", target_msa_a3m="", binder_sequence="",
                    num_recycles=3, use_msa=True,
                )
            proteins = [row["protein"] for row in yaml.safe_load(path.read_text())["sequences"]]
            self.assertTrue(proteins[0]["msa"].startswith("/work/artifacts/raw/optimized_prediction/msas/"))
            self.assertNotIn("msa", proteins[1])
            self.assertEqual(details["reference_msa_count"], 1)
            self.assertEqual(details["reference_msa_chains"][0]["chain_id"], "A")
            self.assertEqual(details["reference_msa_missing_chains"], ["B"])
            self.assertEqual(details["msa_fallback"], "ColabFold server (https://api.colabfold.com)")
            options = workflow._batch_options(
                workflow.ENGINE_PROFILES["boltz2"], seed=1, num_recycles=3,
                engine_options=workflow._validate_engine_options(workflow.ENGINE_PROFILES["boltz2"], {}),
            )
            self.assertEqual(options["cache"], "/models")
            self.assertTrue(options["use_msa_server"])

    def test_empty_yaml_msa_becomes_server_fallback_when_search_is_enabled(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference_files"
            path = root / "input.yaml"
            path.write_text(
                "version: 1\nsequences:\n- protein: {id: A, sequence: ACD, msa: empty}\n",
                encoding="utf-8",
            )
            with patch.object(workflow, "reference_root", return_value=reference):
                details = workflow._attach_reference_boltz_msas(
                    root, path, item_id="sample", use_msa=True,
                )
            protein = yaml.safe_load(path.read_text())["sequences"][0]["protein"]
            self.assertNotIn("msa", protein)
            self.assertEqual(details["reference_msa_missing_chains"], ["A"])

    def test_uploaded_custom_a3m_keeps_precedence_over_reference_cache(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference_files"
            repository = reference / "boltz_models" / "msa_repository"
            repository.mkdir(parents=True)
            cached = repository / f"{hashlib.sha256(b'ACD').hexdigest()}.a3m"
            cached.write_text(">query\nACD\n>reference-hit\nACD\n", encoding="utf-8")
            path = root / "input.yaml"
            path.write_text(
                "version: 1\nsequences:\n- protein: {id: A, sequence: ACD, msa: /work/custom.a3m}\n",
                encoding="utf-8",
            )
            with patch.object(workflow, "reference_root", return_value=reference):
                details = workflow._attach_reference_boltz_msas(
                    root, path, item_id="sample", use_msa=True,
                )
            protein = yaml.safe_load(path.read_text())["sequences"][0]["protein"]
            self.assertEqual(protein["msa"], "/work/custom.a3m")
            self.assertEqual(details["reference_msa_count"], 0)

    def test_identical_protein_chains_share_one_reference_msa_file(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference_files"
            repository = reference / "boltz_models" / "msa_repository"
            repository.mkdir(parents=True)
            cached = repository / f"{hashlib.sha256(b'ACD').hexdigest()}.a3m"
            cached.write_text(">query\nACD\n>hit\nACD\n", encoding="utf-8")
            with patch.object(workflow, "reference_root", return_value=reference):
                path, _details = workflow._stage_native_input(
                    root, profile=workflow.ENGINE_PROFILES["boltz2"], item_id="sample", job_name="example",
                    proteins=[{"id": "A", "sequence": "ACD"}, {"id": "B", "sequence": "ACD"}],
                    dna=[], rna=[], ligand_ccd_codes=[], ligand_smiles=[], seed=1, samples=1,
                    target_pdb_text="", target_chain_id="A", target_msa_a3m="", binder_sequence="",
                    num_recycles=3, use_msa=True,
                )
            proteins = [row["protein"] for row in yaml.safe_load(path.read_text())["sequences"]]
            self.assertEqual(proteins[0]["msa"], proteins[1]["msa"])

    def test_custom_msa_must_match_its_assigned_protein(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "input.yaml"
            path.write_text("version: 1\nsequences:\n- protein: {id: A, sequence: ACD}\n")
            with self.assertRaisesRegex(ValueError, "does not match"):
                workflow._apply_boltz_input_options(root, path, item_id="sample", custom_msas=[
                    {"chain_id": "A", "name": "wrong.a3m", "text": ">query\nAAA\n"},
                ])

    def test_yaml_msa_off_does_not_require_an_unused_alignment_upload(self):
        with TemporaryDirectory() as directory:
            path, _details = workflow._stage_uploaded_boltz_yaml(
                Path(directory), item_id="sample", use_msa=False, auxiliary_files=[],
                yaml_text="version: 1\nsequences:\n- protein: {id: A, sequence: ACD, msa: missing.a3m}\n",
            )
            self.assertEqual(yaml.safe_load(path.read_text())["sequences"][0]["protein"]["msa"], "empty")

    def test_msa_off_also_disables_msa_server(self):
        profile = workflow.ENGINE_PROFILES["boltz2"]
        options = workflow._validate_engine_options(profile, {"use_msa": False, "use_msa_server": True})
        native = workflow._batch_options(profile, seed=1, num_recycles=3, engine_options=options)
        self.assertFalse(native["use_msa_server"])
        self.assertNotIn("generate_ipsae", native)

    def test_ipsae_scores_each_protein_pair_with_the_matching_sample_pae(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            structure = root / "sample_model_2.pdb"
            structure.write_text("END\n")
            pae = root / "pae_sample_model_2.npz"
            pae.touch()
            candidate = {"candidate_id": "sample", "complex_pdb": str(structure), "metrics": {},
                         "target_chains": [], "binder_chains": [],
                         "raw_metadata": {"protein_chains": [{"id": key} for key in "ABC"]}}

            def score(_source, _run, pairs, **_kwargs):
                self.assertEqual(len(pairs), 3)
                for index, pair in enumerate(pairs):
                    self.assertEqual(pair["raw_metadata"]["pae_path"], str(pae))
                    pair["metrics"] = {"ipsae": [0.2, 0.6, 0.4][index]}
                return pairs

            with patch("mn_molecular_design_studio.engine.workflow_engine.workflows.analysis._run_ipsae", side_effect=score):
                workflow._add_boltz_interface_metrics(root, [candidate], workers=2)
            self.assertEqual(candidate["metrics"]["ipsae_min"], 0.2)
            self.assertEqual(candidate["metrics"]["ipsae_max"], 0.6)
            self.assertEqual(candidate["metrics"]["ipsae_avg"], 0.4)
            self.assertEqual(candidate["binder_chains"], [])

    def test_boltz_confidence_and_affinity_sidecars_are_loaded(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            structure = root / "sample_model_0.pdb"
            (root / "confidence_sample_model_0.json").write_text(json.dumps({"iptm": 0.8}))
            (root / "affinity_sample.json").write_text(json.dumps({"affinity_pred_value": -2.1}))
            confidence = workflow._confidence_sidecar(structure)
            self.assertEqual(confidence["iptm"], 0.8)
            self.assertEqual(confidence["affinity_pred_value"], -2.1)


if __name__ == "__main__":
    unittest.main()
