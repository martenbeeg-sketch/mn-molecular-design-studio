from __future__ import annotations

import json

from mn_molecular_design_studio.cluster import capabilities


def _image(*, image_id: str, layer: str = "layer-a", command: list[str] | None = None, size: int = 10):
    return {
        "Id": image_id,
        "Size": size,
        "Os": "linux",
        "Architecture": "amd64",
        "Variant": "",
        "OsVersion": "",
        "Config": {
            "Entrypoint": command or ["/usr/bin/run"],
            "Env": ["PATH=/usr/bin"],
        },
        "RootFS": {"Type": "layers", "Layers": [layer]},
        "RepoTags": ["mn-ligandmpnn:latest"],
    }


def test_image_fingerprint_ignores_store_specific_id_and_size():
    classic = _image(image_id="sha256:classic-manifest", size=100)
    containerd = _image(image_id="sha256:containerd-config", size=200)

    assert capabilities.image_content_fingerprint(classic) == capabilities.image_content_fingerprint(containerd)


def test_image_fingerprint_changes_when_config_or_layer_changes():
    baseline = _image(image_id="sha256:one")

    assert capabilities.image_content_fingerprint(baseline) != capabilities.image_content_fingerprint(
        _image(image_id="sha256:two", layer="layer-b")
    )
    assert capabilities.image_content_fingerprint(baseline) != capabilities.image_content_fingerprint(
        _image(image_id="sha256:three", command=["/usr/bin/other"])
    )


def test_images_returns_portable_fingerprints_for_each_tag(monkeypatch):
    first = _image(image_id="sha256:classic-manifest")
    second = _image(image_id="sha256:containerd-config", layer="layer-b")
    first["RepoTags"] = ["mn-ligandmpnn:latest", "mn-ligandmpnn:backup"]
    second["RepoTags"] = ["mn-foundry:cu128"]

    def fake_run(command, **_kwargs):
        if command[:4] == ["docker", "image", "ls", "--format"]:
            return type("Completed", (), {"stdout": "mn-ligandmpnn:latest\nmn-foundry:cu128\n"})()
        assert command[:4] == ["docker", "image", "inspect", "--format"]
        return type(
            "Completed",
            (),
            {"stdout": "\n".join(json.dumps(item) for item in (first, second))},
        )()

    monkeypatch.setattr(capabilities.subprocess, "run", fake_run)

    inventory = capabilities.images()

    assert inventory == {
        "mn-ligandmpnn:latest": capabilities.image_content_fingerprint(first),
        "mn-ligandmpnn:backup": capabilities.image_content_fingerprint(first),
        "mn-foundry:cu128": capabilities.image_content_fingerprint(second),
    }


def test_rfantibody_generation_only_requires_only_diffusion_weights():
    from mn_molecular_design_studio.rfantibody_input import IMAGE, MIMOSA_IMAGE

    installed = {IMAGE: "native-image", MIMOSA_IMAGE: "mimosa-image"}
    for workflow, mimosa_enabled, image in (
        ("design.rfantibody", False, IMAGE),
        ("design.rfantibody_mimosa", True, MIMOSA_IMAGE),
    ):
        needs = capabilities.requirements(
            {
                "workflow": workflow,
                "parameters": {
                    "generation_only": True,
                    "mimosa_enabled": mimosa_enabled,
                    "gpu_device": "0",
                },
            },
            installed_images=installed,
        )

        assert needs["images"] == {image: installed[image]}
        assert "rfantibody/RFdiffusion_Ab.pt" in needs["references"]
        assert "rfantibody/ProteinMPNN_v48_noise_0.2.pt" not in needs["references"]
        assert "rfantibody/RF2_ab.pt" not in needs["references"]


def test_rfantibody_full_pipeline_retains_required_sequence_and_refold_weights():
    from mn_molecular_design_studio.rfantibody_input import IMAGE

    needs = capabilities.requirements(
        {
            "workflow": "design.rfantibody",
            "parameters": {"generation_only": False, "use_abmpnn": False, "gpu_device": "0"},
        },
        installed_images={IMAGE: "native-image"},
    )

    assert {
        "rfantibody/RFdiffusion_Ab.pt",
        "rfantibody/ProteinMPNN_v48_noise_0.2.pt",
        "rfantibody/RF2_ab.pt",
    }.issubset(needs["references"])


def test_surface_bundle_placement_requires_each_native_tool_and_external_model():
    from mn_molecular_design_studio.engine.workflow_engine.core.manifests import load_manifest
    from mn_molecular_design_studio.engine.workflow_engine.workflows.detection import SURFACE_METHODS, SURFACE_REQUIRED_REFERENCES

    needs = capabilities.requirements({"workflow": "detection.surface_analysis", "parameters": {}}, installed_images={})
    assert set(needs["images"]) == {load_manifest(tool)["image"] for tool in SURFACE_METHODS}
    assert set(needs["references"]) == {path for paths in SURFACE_REQUIRED_REFERENCES.values() for path in paths}


def test_abmpnn_placement_requires_ligandmpnn_runtime_and_antibody_checkpoint():
    from mn_molecular_design_studio.mpnn_input import IMAGE

    for workflow, parameters in (("sequence.abmpnn", {}), ("sequence.mpnn", {"model_type": "ab_mpnn"})):
        needs = capabilities.requirements(
            {"workflow": workflow, "parameters": parameters}, installed_images={IMAGE: "ligandmpnn-hash"},
        )

        assert needs["images"] == {IMAGE: "ligandmpnn-hash"}
        assert "mpnn/abmpnn.pt" in needs["references"]


def test_boltzgen_generator_placement_uses_the_selected_runtime():
    from mn_molecular_design_studio.engine.workflow_engine.core.manifests import load_manifest
    from mn_molecular_design_studio.engine.workflow_engine.workflows.design import MIMOSA_BOLTZGEN_IMAGE

    standard_image = load_manifest("boltzgen")["image"]
    installed = {standard_image: "standard-hash", MIMOSA_BOLTZGEN_IMAGE: "mimosa-hash"}
    for workflow, parameters, image in (
        ("design.boltzgen_generator", {}, standard_image),
        ("design.boltzgen", {}, standard_image),
        ("design.boltzgen_mimosa", {"mimosa_condition": "sequence_fixed"}, MIMOSA_BOLTZGEN_IMAGE),
    ):
        needs = capabilities.requirements({"workflow": workflow, "parameters": parameters}, installed_images=installed)

        assert needs["images"] == {image: installed[image]}
        assert "boltzgen-cache" in needs["references"]


def test_rfdiffusion3_placement_does_not_require_classic_model_directory():
    needs = capabilities.requirements({"workflow": "design.rfdiffusion3", "parameters": {}}, installed_images={})

    assert "foundry/rfd3_latest.ckpt" in needs["references"]
    assert "rfdiffusion_models" not in needs["references"]
    assert needs["gpu_count"] == 1


def test_ligand_generator_placement_includes_its_native_image_and_checkpoints():
    for engine, image, paths in (
        ("pocketflow", "mn-pocketflow:cu128-a31a5a0", {"generation/pocketflow/ZINC-pretrained-255000.pt"}),
        ("pgmg", "mn-pgmg:cu128-85fb712", {"generation/pgmg/model.pth", "generation/pgmg/tokenizer.pkl"}),
    ):
        needs = capabilities.requirements({"workflow": f"ligand.generation.{engine}", "parameters": {}},
                                          installed_images={})
        assert set(needs["images"]) == {image}
        assert set(needs["references"]) == paths
        assert needs["gpu_count"] == 1


def test_msa_enabled_prediction_placement_includes_its_colabfold_image():
    from mn_molecular_design_studio.cofolding_msa import colabfold_image

    for engine, options, validation, needs_colabfold in (
        ("esmfold2", {}, {}, True),
        ("openfold3_jax", {}, {}, True),
        ("openfold3_jax", {"use_msa": False}, {}, False),
        ("boltz2", {}, {}, True),
        ("esmfold2", {"use_msa": False}, {}, False),
        ("boltz2", {"use_msa_server": False}, {}, False),
        ("localfold", {"msa_mode": "mmseqs2_server"}, {}, True),
        ("localfold", {"msa_mode": "reference_only"}, {}, False),
        ("boltz2", {}, {"validation_sources": [{}], "validation_mode": "binder"}, False),
        ("boltz2", {"use_msa": False}, {"validation_sources": [{}], "target_use_msa": True}, True),
    ):
        needs = capabilities.requirements({"workflow": "refolding.optimized_prediction", "parameters": {
            "engine": engine, "engine_options": options, **validation,
        }}, installed_images={})
        assert (colabfold_image() in needs["images"]) == needs_colabfold


def test_native_candidate_refolding_placement_requires_its_image_and_model_files():
    from mn_molecular_design_studio.engine.workflow_engine.workflows import refolding

    for workflow, image, path in (
        ("refolding.openfold3", refolding.OPENFOLD3_IMAGE, refolding.OPENFOLD3_REFERENCE),
        ("refolding.rf3", refolding.RF3_IMAGE, refolding.RF3_REFERENCE),
        ("refolding.protenix", refolding.PROTENIX_IMAGE, "pxdesign"),
        ("refolding.protenix_cli", refolding.PROTENIX_CLI_IMAGE, "protenix"),
        ("refolding.boltzgen_fold", refolding.BOLTZGEN_IMAGE, "boltzgen-cache/" + refolding.BOLTZGEN_FOLD_CHECKPOINT),
    ):
        needs = capabilities.requirements({"workflow": workflow, "parameters": {}}, installed_images={})
        assert set(needs["images"]) == {image}
        assert path in needs["references"]
        assert needs["gpu_count"] == 1


def test_reference_fingerprints_detect_changed_model_contents(tmp_path, monkeypatch):
    monkeypatch.setenv('MN_MOLECULAR_DESIGN_STUDIO_REFERENCE_DIR', str(tmp_path))
    path = tmp_path / 'model.ckpt'; path.write_bytes(b'weights-v1')
    before = capabilities.reference_fingerprints(['model.ckpt'])
    path.write_bytes(b'weights-v2')
    after = capabilities.reference_fingerprints(['model.ckpt'])
    assert before != after
    from mn_molecular_design_studio.cluster.store import Store
    node = {'availability': 'online', 'snapshot': {'accept_jobs': True, 'scheduler_online': True, 'reference_hashes': before}}
    assert Store.compatible(node, {'reference_hashes': before})
    assert not Store.compatible(node, {'reference_hashes': after})
