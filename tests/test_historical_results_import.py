import csv
import json
import zipfile

import pytest

from mn_molecular_design_studio.historical_results_import import preview_historical_campaign


def _write_csv(path, headers, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)


def _pdb_for_chains(chains):
    from Bio.SeqUtils import seq3

    lines, serial = [], 1
    offset = 0.0
    for chain, sequence in chains.items():
        for index, letter in enumerate(sequence, 1):
            lines.append(f"ATOM  {serial:5d}  CA  {seq3(letter).upper():3s} {chain}{index:4d}    "
                         f"{offset + index:8.3f}{1.0:8.3f}{1.0:8.3f}{1.00:6.2f}{20.00:6.2f}           C")
            serial += 1
        lines.append("TER")
        offset += 10.0
    return "\n".join([*lines, "END", ""])


def _bindcraft2_campaign(root):
    trajectory = root / "1_Trajectories"
    refolded = root / "2_Refolded"
    ranked = root / "3_Ranked"
    hashes = ["aaaaaaaaaaaaaaaa", "bbbbbbbbbbbbbbbb", "cccccccccccccccc"]
    _write_csv(trajectory / "!_Trajectories.csv", ["design", "hash", "terminated"], [
        {"design": f"campaign_binder_l4_{hash_}", "hash": hash_, "terminated": ""} for hash_ in hashes
    ])
    refold_rows = [
        {"design": f"campaign_binder_l4_{hashes[0]}_candidate0", "hash": hashes[0], "Binder_Sequence": "ACDE", "outcome": "passed", "failed_filters": "", "i_pTM": "0.81"},
        {"design": f"campaign_binder_l4_{hashes[1]}_candidate0", "hash": hashes[1], "Binder_Sequence": "FGHI", "outcome": "passed", "failed_filters": "", "i_pTM": "0.72"},
        {"design": f"campaign_binder_l4_{hashes[2]}_candidate0", "hash": hashes[2], "Binder_Sequence": "KLMN", "outcome": "failed", "failed_filters": "i_pTM", "i_pTM": "0.31"},
    ]
    _write_csv(refolded / "!_Refolded.csv", ["design", "hash", "Binder_Sequence", "outcome", "failed_filters", "i_pTM"], refold_rows)
    ranked_row = {"design": f"campaign_binder_l4_{hashes[0]}_seq0", "hash": hashes[0], "Binder_Sequence": "ACDE", "rank": "1", "i_pTM": "0.81"}
    _write_csv(ranked / "!_Ranked.csv", list(ranked_row), [ranked_row])
    complex_dir = refolded / "Complexes"
    complex_dir.mkdir(parents=True, exist_ok=True)
    (complex_dir / f"campaign_binder_l4_{hashes[0]}_candidate0.pdb").write_text(
        _pdb_for_chains({"A": "ACDEFG", "B": "ACDE"}))
    (ranked / f"{ranked_row['design']}.pdb").write_text(
        _pdb_for_chains({"A": "ACDEFG", "B": "ACDE"}))
    (root / "target.pdb").write_text(_pdb_for_chains({"A": "ACDEFG"}))
    (root / "campaign_metadata.json").write_text(json.dumps({
        "version": "BindCraft 2 v1.0.3",
        "settings": {"modality": "binder", "filters": {"i_pTM": {"threshold": 0.7, "higher": True}},
                     "targets": [{"target_sequence": "ACDEFG", "target_path": "target.pdb", "chains": ["A"], "hotspots": "A10,A12"}]},
    }))
    return root


def test_bindcraft2_keeps_trajectory_refold_pass_and_rank_counts_separate(tmp_path):
    root = _bindcraft2_campaign(tmp_path / "campaign")
    preview = preview_historical_campaign(root)

    assert preview["engine"] == "bindcraft2"
    assert preview["engine_version"] == "BindCraft 2 v1.0.3"
    assert preview["counts"] == {
        "trajectories_started": 3, "refolded_candidates": 3,
        "native_filter_passes": 2, "native_filter_rejected": 1,
        "native_filter_unknown": 0, "ranked_designs": 1,
    }
    assert len(preview["candidates"]) == 3
    assert len(preview["proposed_identities"]) == 3
    ranked = [row for row in preview["candidates"] if row["ranked"]]
    assert len(ranked) == 1
    assert ranked[0]["sequence"] == "ACDE"
    assert ranked[0]["native_filter_status"] == "passed"
    assert ranked[0]["association_status"] == "exact_hash_and_sequence"
    assert ranked[0]["ranked_design_id"].endswith("_seq0")
    assert ranked[0]["structure"]["stage"] == "refolded"
    assert ranked[0]["structure"]["status"] == "available"
    assert ranked[0]["ranked_structure"]["stage"] == "ranked"
    assert ranked[0]["ranked_structure"]["status"] == "available"
    assert ranked[0]["native_metrics"]["i_pTM"] == "0.81"
    assert ranked[0]["retention_selected"] is False
    assert preview["target"]["sequence_or_chain_sequences"] == "ACDEFG"
    assert preview["target"]["source_structure"]["sha256"]
    assert preview["target"]["source_structure"]["experimental_or_native_claim"] is False
    assert preview["production_import_enabled"] is False


def test_bindcraft2_ranked_row_requires_exact_hash_and_sequence_join(tmp_path):
    root = _bindcraft2_campaign(tmp_path / "campaign")
    ranked_csv = root / "3_Ranked" / "!_Ranked.csv"
    _write_csv(ranked_csv, ["design", "hash", "Binder_Sequence", "rank"], [
        {"design": "unmatched_ranked", "hash": "dddddddddddddddd", "Binder_Sequence": "ACDE", "rank": "1"},
    ])

    preview = preview_historical_campaign(root)
    assert preview["counts"]["ranked_designs"] == 1
    assert len(preview["candidates"]) == 4
    unmatched = next(row for row in preview["candidates"] if row["native_stage"] == "ranked")
    assert unmatched["association_status"] == "unmatched_refolded"
    assert unmatched["native_filter_status"] == "unknown"
    assert any("no refolded row matching its exact hash and sequence" in text for text in preview["warnings"])


def test_preview_selection_is_explicit_and_does_not_promote_other_rows(tmp_path):
    root = _bindcraft2_campaign(tmp_path / "campaign")
    first = preview_historical_campaign(root)
    selected_id = first["candidates"][0]["candidate_id"]
    preview = preview_historical_campaign(root, selected_candidate_ids=[selected_id])

    selected = [row for row in preview["candidates"] if row["retention_selected"]]
    assert [row["candidate_id"] for row in selected] == [selected_id]
    assert sum(identity["retention_selected"] for identity in preview["proposed_identities"]) == 1
    assert preview["proposed_collection"]["status"] == "preview_only_not_created"
    with pytest.raises(ValueError, match="not present"):
        preview_historical_campaign(root, selected_candidate_ids=["unknown-candidate"])


def test_source_content_identity_is_stable_across_zip_name_and_location(tmp_path):
    root = _bindcraft2_campaign(tmp_path / "original")
    folder_preview = preview_historical_campaign(root)
    archive_path = tmp_path / "renamed-campaign.zip"
    with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in root.rglob("*"):
            if path.is_file():
                archive.write(path, f"renamed-root/{path.relative_to(root).as_posix()}")

    archive_preview = preview_historical_campaign(archive_path)
    assert archive_preview["source_content_id"] == folder_preview["source_content_id"]
    assert archive_preview["source_file_count"] == folder_preview["source_file_count"]


def test_unsafe_zip_paths_are_rejected_before_reading_campaign(tmp_path):
    archive_path = tmp_path / "unsafe.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("../outside.txt", "not allowed")
    with pytest.raises(ValueError, match="unsafe path"):
        preview_historical_campaign(archive_path)


def test_historical_source_rejects_symlinked_files_and_root_symlinks(tmp_path):
    campaign = _bindcraft2_campaign(tmp_path / "campaign")
    outside = tmp_path / "external.csv"
    outside.write_text("design,hash\nexternal,deadbeef\n")
    linked = campaign / "1_Trajectories" / "external.csv"
    linked.symlink_to(outside)
    with pytest.raises(ValueError, match="may not contain symlinks"):
        preview_historical_campaign(campaign)
    linked.unlink()
    alias = tmp_path / "campaign-link"
    alias.symlink_to(campaign, target_is_directory=True)
    with pytest.raises(ValueError, match="may not be symlinks"):
        preview_historical_campaign(alias)


def test_historical_zip_rejects_duplicate_normalized_paths(tmp_path):
    archive_path = tmp_path / "duplicate.zip"
    with pytest.warns(UserWarning, match="Duplicate name"):
        with zipfile.ZipFile(archive_path, "w") as archive:
            archive.writestr("campaign/a.txt", "first")
            archive.writestr("campaign/a.txt", "second")
    with pytest.raises(ValueError, match="duplicate path"):
        preview_historical_campaign(archive_path)


def test_historical_preview_requires_actual_matching_coordinates_for_structures(tmp_path):
    root = _bindcraft2_campaign(tmp_path / "campaign")
    first = root / "2_Refolded" / "Complexes" / "campaign_binder_l4_aaaaaaaaaaaaaaaa_candidate0.pdb"
    first.write_text("not a coordinate file\n")
    original = preview_historical_campaign(root)
    preview = preview_historical_campaign(root, selected_candidate_ids=[original["candidates"][0]["candidate_id"]])
    row = preview["candidates"][0]
    assert row["structure"]["status"] == "invalid_coordinates"
    assert all(not item["path"].startswith("2_Refolded/") for item in preview["proposed_structures"])
    assert any(item["status"] == "invalid_coordinates" for item in preview["unresolved_structures"])

    first.write_text(_pdb_for_chains({"A": "ACDEFG", "B": "KLMN"}))
    mismatch = preview_historical_campaign(root)
    assert mismatch["candidates"][0]["structure"]["status"] == "binder_sequence_mismatch"


def test_bindcraft1_native_tables_are_detected_without_inventing_acceptance(tmp_path):
    root = tmp_path / "bindcraft1"
    _write_csv(root / "final_design_stats.csv", ["design", "Binder_Sequence", "i_pTM"], [
        {"design": "final-A", "Binder_Sequence": "ACDE", "i_pTM": "0.8"},
    ])
    _write_csv(root / "mpnn_design_stats.csv", ["design", "Binder_Sequence"], [
        {"design": "mpnn-A", "Binder_Sequence": "ACDE"},
    ])
    _write_csv(root / "trajectory_stats.csv", ["trajectory"], [{"trajectory": "traj-A"}])
    (root / "final-A_model0.pdb").write_text(_pdb_for_chains({"B": "ACDE"}))
    (root / "2_final-A_model1.pdb").write_text(_pdb_for_chains({"B": "ACDE"}))

    preview = preview_historical_campaign(root)
    assert preview["engine"] == "bindcraft1"
    assert preview["counts"]["trajectories_started"] == 1
    assert preview["counts"]["refolded_candidates"] == 1
    assert preview["counts"]["ranked_designs"] == 0
    assert preview["counts"]["final_designs"] == 1
    assert preview["candidates"][0]["native_filter_status"] == "unknown"
    assert preview["candidates"][0]["structure"]["status"] == "multiple_models"
    assert preview["candidates"][0]["structure_layout_status"] == "supported_native_models"
    assert [item["layout"] for item in preview["candidates"][0]["structures"]] == [
        "native_model_suffix", "native_model_suffix",
    ]
    assert preview["candidates"][0]["native_metrics"]["i_pTM"] == "0.8"


def test_bindcraft1_reports_unrecognized_native_metric_definitions_and_layouts(tmp_path):
    root = tmp_path / "bindcraft1"
    _write_csv(root / "final_design_stats.csv", ["design", "Binder_Sequence", "Average_i_pTM"], [
        {"design": "final-A", "Binder_Sequence": "ACDE", "Average_i_pTM": "0.8"},
    ])
    _write_csv(root / "mpnn_design_stats.csv", ["design"], [{"design": "mpnn-A"}])
    _write_csv(root / "trajectory_stats.csv", ["trajectory"], [{"trajectory": "traj-A"}])
    (root / "final-A_model2.pdb").write_text(_pdb_for_chains({"B": "KLMN"}))

    preview = preview_historical_campaign(root)
    candidate = preview["candidates"][0]
    assert candidate["structure_layout_status"] == "recognized_layout_unvalidated_coordinates"
    assert candidate["structures"][0]["status"] == "binder_sequence_mismatch"
    assert "Average_i_pTM" in candidate["native_metric_compatibility"]["unsupported_native_fields"]
    assert candidate["native_filter_status"] == "unknown"


def test_optional_database_identity_lookup_is_read_only_and_acl_limited(tmp_path):
    from mn_molecular_design_studio.molecular_database import MolecularDatabase
    from mn_molecular_design_studio.molecule_imports import parse_sequences

    root = _bindcraft2_campaign(tmp_path / "campaign")
    database = MolecularDatabase(tmp_path / "database")
    imported = database.import_records({"id": "owner", "role": "user"}, parse_sequences("ACDE"), name="Owner library")
    before = database.list_molecules({"id": "owner", "role": "user"})["total"]
    preview = preview_historical_campaign(root, database_path=database.path,
                                          user={"id": "owner", "role": "user"})
    after = database.list_molecules({"id": "owner", "role": "user"})["total"]
    match = next(identity for identity in preview["proposed_identities"]
                 if identity["identity_id"] == imported["molecule_ids"][0])
    assert before == after == 1
    assert match["database_action"] == "reuse_existing_identity"
    assert match["existing_matches"][0]["group_id"] == imported["group_id"]
    hidden = preview_historical_campaign(root, database_path=database.path,
                                        user={"id": "other", "role": "user"})
    hidden_match = next(identity for identity in hidden["proposed_identities"]
                        if identity["identity_id"] == imported["molecule_ids"][0])
    assert hidden_match["existing_matches"] == []
