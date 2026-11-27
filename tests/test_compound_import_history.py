import json

from mn_molecular_design_studio import datasets, molecule_imports
from mn_molecular_design_studio.compound_import_history import apply_reviews, archived_reviews
from test_database_browser import USER, database


def test_archived_reviews_require_exact_source_and_restore_existing_batch(database, tmp_path):
    source = b"ID,SMILES\nv1,CCO\nv2,broken\nv3,broken\n"
    audit = molecule_imports.analyze_small_molecule_table(source, "vendor.csv", smiles_column="SMILES", id_column="ID")
    origin = datasets.input_origin(USER, [])
    result = database.import_records(USER, datasets.input_records(audit["records"]), name="Vendor library",
                                     origin=origin, import_audit=audit, source_content=source)
    group_id, batch_id = result["group_id"], result["import_batch_id"]
    run = tmp_path / "compound-import" / "vendor-import"
    originals = run / "artifacts/compounds"
    originals.mkdir(parents=True)
    (originals / "vendor.csv").write_bytes(source)
    reports = run / "artifacts/reports"
    reports.mkdir()
    (reports / "validation_report.json").write_text(json.dumps(audit["mapping"]))
    (reports / "rejected_compounds.csv").write_text("compound_id,smiles,source_row\nv2,broken,3\nv3,broken,4\n")
    for vendor_id, row_number, cid, replacement in (("v2", 3, 702, "CCO"), ("v3", 4, 6736, "C1CC1")):
        review = tmp_path / "compound-review" / f"review-{vendor_id}"
        evidence = review / "artifacts/reports"
        evidence.mkdir(parents=True)
        (review / "metadata.json").write_text(json.dumps({"source_compound_run_id": run.name, "status": "completed"}))
        (evidence / "pubchem_review.json").write_text(json.dumps({"source_run_id": run.name, "compound_id": vendor_id,
            "source_row": row_number, "decision": "accepted", "candidate": {"cid": cid, "smiles": replacement}}))

    assert archived_reviews(source.replace(b"v1", b"v9"), audit, [tmp_path]) == {}
    recovered = archived_reviews(source, audit, [tmp_path])
    assert list(recovered) == [2, 3]
    revised = apply_reviews(source, audit, recovered)
    assert database.restore_compound_reviews(USER, group_id, batch_id, revised) == {"restored_rows": 2, "added_molecules": 1}
    assert database.restore_compound_reviews(USER, group_id, batch_id, revised) == {"restored_rows": 0, "added_molecules": 0}
    first, duplicate, inserted = database.import_batch_rows(USER, batch_id)
    assert duplicate["molecule_id"] == first["molecule_id"] and duplicate["duplicate_of_row"] == 1
    assert inserted["review"]["original_status"] == "invalid" and inserted["review"]["review_run_id"] == "review-v3"
    assert inserted["raw_smiles"] == "broken" and inserted["canonical_smiles"] == "C1CC1"
    view = datasets.DatasetView(database, group_id)
    assert {item["name"] for item in view.list_molecules(USER)["items"]} == {"v1", "v3"}
    assert any(source.get("compound_review", {}).get("candidate", {}).get("cid") == 702
               for item in view.list_molecules(USER)["items"] for source in item["sources"])
    assert database.import_source_file(USER, batch_id).read_bytes() == source
