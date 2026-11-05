"""Optional real Chaperone CATCIF and Parquet roundtrip through the Studio bridge."""
from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("chaperone")
pytest.importorskip("catcif_tools")

from chaperone.data.schema import StorageFormat  # noqa: E402
from chaperone.data.handler import PipelineDataHandler  # noqa: E402
from chaperone.data.schema import make_row  # noqa: E402
from chaperone.io import parse_to_atom_array  # noqa: E402
from mn_molecular_design_studio.chaperone_adapter import (  # noqa: E402
    from_pipeline_data_handler,
    materialize_catcif_structure,
    to_pipeline_data_handler,
)


def test_candidate_sidecar_survives_real_catcif_store_and_parquet_reload(tmp_path):
    structure = Path(__file__).parents[1] / "mn_molecular_design_studio/assets/vina/1iep.pdb"
    source = [{
        "candidate_id": "variant-2-model-0", "source_tool": "boltz2", "stage": "complex_refolding",
        "complex_pdb": str(structure), "binder_sequence": "ACDEFG", "parents": ["variant-2"],
        "target_chains": ["A"], "binder_chains": ["B"], "metrics": {"iptm": 0.8},
        "raw_metadata": {"model_index": 0, "source_token": "preserve-me"},
    }]
    handler, sidecar = to_pipeline_data_handler(source)
    original_hash = structure.read_bytes()

    handler.store("boltz2", tmp_path / "structures.catcif", cleanup=False)
    stored = handler._df.row(0, named=True)["structures"][0]
    assert stored["storage_format"] == str(StorageFormat.CATCIF)
    assert ".catcif:" in stored["path"]
    assert Path(stored["store_path"]).is_file()

    parquet_path = tmp_path / "chaperone-handler.parquet"
    handler.write_parquet(parquet_path)
    reloaded = PipelineDataHandler(parquet_path, allow_sequences=True)
    restored = from_pipeline_data_handler(reloaded, sidecar)
    assert structure.read_bytes() == original_hash
    assert len(restored) == 1
    candidate = restored[0]
    assert candidate["candidate_id"] == "variant-2-model-0"
    assert candidate["parents"] == ["variant-2"]
    assert candidate["complex_pdb"] == str(structure)
    assert candidate["raw_metadata"]["source_token"] == "preserve-me"
    receipt = candidate["raw_metadata"]["chaperone_interop"]["chaperone_record"]
    assert receipt["structures"][0]["storage_format"] == str(StorageFormat.CATCIF)
    assert ".catcif:" in receipt["structures"][0]["path"]
    assert receipt["design_id"] == "variant-2-model-0"


def test_distinct_catcif_models_roundtrip_to_studio_and_biotite(tmp_path, design_source_factory):
    import copy
    import numpy as np
    from Bio.PDB import MMCIFIO, MMCIFParser
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from chaperone.data.catcif import list_catcif_entries
    from mn_molecular_design_studio import candidate_sets as sets
    from mn_molecular_design_studio.complex_analysis_handoff import build_collection
    from mn_molecular_design_studio.engine.workflow_engine.workflows.quality_runner import interface_metrics

    source_run, _ = design_source_factory(tmp_path / "source-run")
    original = source_run / "designed.cif"
    first = tmp_path / "designed-complete.cif"
    cif = MMCIF2Dict(str(original))
    cif["_atom_site.occupancy"] = ["1.0"] * len(cif["_atom_site.id"])
    writer = MMCIFIO()
    writer.set_dict(cif)
    writer.save(str(first))
    shifted = tmp_path / "designed-shifted.cif"
    structure = MMCIFParser(QUIET=True).get_structure("shifted", str(first))
    next(structure.get_atoms()).coord[0] += 0.25
    writer = MMCIFIO()
    writer.set_structure(structure)
    writer.save(str(shifted))
    original_bytes, first_bytes, shifted_bytes = original.read_bytes(), first.read_bytes(), shifted.read_bytes()
    candidates = [
        {"candidate_id": f"variant-2-model-{model}", "source_tool": "boltz2",
         "stage": "complex_refolding", "complex_pdb": str(path), "binder_sequence": "KLMNP",
         "target_chains": ["A"], "binder_chains": ["B"], "parents": ["variant-2"],
         "metrics": {"iptm": 0.8 + model / 100},
         "raw_metadata": {"model_index": model, "target_sequences": {"A": "ACDEFG"}}}
        for model, path in enumerate((first, shifted))
    ]
    handler, sidecar = to_pipeline_data_handler(candidates)
    archive = tmp_path / "structures.catcif"
    handler.store("boltz2", archive, cleanup=False)
    parquet = tmp_path / "pipeline.parquet"
    handler.write_parquet(parquet)
    restored = from_pipeline_data_handler(PipelineDataHandler(parquet, allow_sequences=True), sidecar)

    tags = [row["raw_metadata"]["chaperone_interop"]["chaperone_record"]["structures"][0]["path"]
            for row in restored]
    assert len(set(tags)) == 2 and set(list_catcif_entries(archive)) == {tag.split(":", 1)[1] for tag in tags}
    for source, tag in zip((original, shifted), tags):
        np.testing.assert_allclose(parse_to_atom_array(str(source)).coord,
                                   parse_to_atom_array(tag).coord, atol=1e-3)
    assert [row["candidate_id"] for row in restored] == ["variant-2-model-0", "variant-2-model-1"]
    assert [row["raw_metadata"]["chaperone_interop"]["model_id"] for row in restored] == [0, 1]
    assert [row["metrics"]["iptm"] for row in restored] == [0.8, 0.81]

    target = {"target_id": "catcif-target", "target_type": "protein", "target_pdb": str(original),
              "target_chains": ["A"]}
    rows = []
    for row in restored:
        row = copy.deepcopy(row)
        info = sets.source_record(row, {"run_dir": str(source_run), "source_run_id": source_run.name,
                                        "source_run_dir": str(source_run)}, "catcif-source", target)
        info["target_pdb"] = str(original)
        rows.append(info)
    collection = sets.write_collection(tmp_path / "source", "catcif", rows, target)
    source = {**collection, "source_run": source_run.name, "source_run_dir": str(source_run)}
    routed = build_collection(tmp_path / "analysis", "complex_analysis", source,
        {"entry_mode": "manual_selection", "selected_candidate_ids": [row["candidate_id"] for row in rows]})
    entries = list(sets.iter_candidates(routed))
    assert len(entries) == 2
    assert {str(row["raw_metadata"]["model_index"]) for row in entries} == {"0", "1"}
    assert all(row["raw_metadata"]["complex_analysis_entry"]["structure_storage"]["format"] == "catcif"
               for row in entries)
    assert all(Path(row["complex_pdb"]).is_file() for row in entries)
    destination = tmp_path / "analysis-output"
    destination.mkdir()
    measured = interface_metrics({"structure": entries[0]["complex_pdb"], "target_chains": ["A"],
                                  "binder_chains": ["B"]}, {"contact_distance": 5.0}, destination)
    assert measured["qc_interface_contact_pairs"] > 0
    assert original.read_bytes() == original_bytes and first.read_bytes() == first_bytes and shifted.read_bytes() == shifted_bytes
