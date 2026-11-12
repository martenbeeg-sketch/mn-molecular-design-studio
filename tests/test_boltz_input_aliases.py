import hashlib

import pytest

from mn_molecular_design_studio.engine.workflow_engine.core.artifact_collection import validate_run_artifacts
from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import _materialize_boltz_input_aliases


@pytest.fixture
def staged_input(tmp_path):
    canonical = tmp_path / "artifacts/raw/optimized_prediction/inputs/item.yaml"
    canonical.parent.mkdir(parents=True)
    canonical.write_bytes(b"version: 1\nsequences: []\n")
    staging = tmp_path / "artifacts/raw/optimized_prediction/output/.mn-predict-boltz2-inputs"
    staging.mkdir(parents=True)
    alias = staging / "item.yaml"
    target = "/work/" + canonical.relative_to(tmp_path).as_posix()
    alias.symlink_to(target)
    items = [{"id": "item", "input": target}]
    identities = [{"path": str(canonical), "sha256": hashlib.sha256(canonical.read_bytes()).hexdigest()}]
    return tmp_path, canonical, alias, items, identities


def test_exact_native_input_alias_is_materialized_without_changing_science(staged_input):
    run, canonical, alias, items, identities = staged_input
    structure = alias.parent.parent / "prediction.pdb"
    structure.write_bytes(b"REMARK retained native coordinates\n")
    before = {path: path.read_bytes() for path in run.rglob("*") if path.is_file()}
    assert validate_run_artifacts(run, result={})["status"] == "incomplete"
    _materialize_boltz_input_aliases(run, items, identities)
    assert not alias.is_symlink()
    assert alias.read_bytes() == canonical.read_bytes()
    assert all(path.read_bytes() == contents for path, contents in before.items())
    assert validate_run_artifacts(run, result={})["status"] == "complete"
    _materialize_boltz_input_aliases(run, items, identities)


@pytest.mark.parametrize("target", ["/etc/passwd", "/work/other-run/item.yaml", "../inputs/item.yaml",
                                         "/work/artifacts/raw/optimized_prediction/inputs/../inputs/item.yaml"])
def test_unexpected_alias_target_is_rejected(staged_input, target):
    run, _, alias, items, identities = staged_input
    alias.unlink()
    alias.symlink_to(target)
    with pytest.raises(ValueError, match="unexpected target"):
        _materialize_boltz_input_aliases(run, items, identities)
    assert alias.is_symlink()


@pytest.mark.parametrize("change", ["missing", "modified", "symlink", "parent-symlink"])
def test_canonical_input_must_be_original_regular_file(staged_input, change):
    run, canonical, alias, items, identities = staged_input
    contents = canonical.read_bytes()
    if change == "missing":
        canonical.unlink()
    elif change == "modified":
        canonical.write_text("changed")
    elif change == "symlink":
        canonical.unlink()
        replacement = run / "replacement.yaml"
        replacement.write_bytes(contents)
        canonical.symlink_to(replacement)
    else:
        canonical.parent.rename(run / "moved-inputs")
        canonical.parent.symlink_to(run / "moved-inputs", target_is_directory=True)
    with pytest.raises(ValueError, match="Canonical Boltz input"):
        _materialize_boltz_input_aliases(run, items, identities)
    assert alias.is_symlink()


def test_redirected_staging_directory_is_rejected(staged_input):
    run, _, alias, items, identities = staged_input
    moved = run / "elsewhere"
    alias.parent.rename(moved)
    alias.parent.symlink_to(moved, target_is_directory=True)
    with pytest.raises(ValueError, match="staging"):
        _materialize_boltz_input_aliases(run, items, identities)


def test_configured_run_root_symlink_does_not_lose_input_identity(staged_input):
    run, canonical, alias, items, identities = staged_input
    configured_root = run.parent / (run.name + "-configured")
    configured_root.symlink_to(run, target_is_directory=True)
    identities[0]["path"] = str(configured_root / canonical.relative_to(run))
    _materialize_boltz_input_aliases(configured_root, items, identities)
    assert not alias.is_symlink()
    assert alias.read_bytes() == canonical.read_bytes()


def test_unknown_alias_and_missing_scientific_output_stay_failures(staged_input):
    run, _, alias, items, identities = staged_input
    (alias.parent / "unknown.yaml").symlink_to("/work/unrecognized.yaml")
    prediction = alias.parent.parent / "prediction.pdb"
    prediction.symlink_to("/work/missing.pdb")
    _materialize_boltz_input_aliases(run, items, identities)
    assert validate_run_artifacts(run, result={})["status"] == "incomplete"
    assert prediction.is_symlink()
    assert (alias.parent / "unknown.yaml").is_symlink()


@pytest.mark.parametrize("item_id", ["../escape", ".", "..", "/escape"])
def test_input_ids_cannot_escape_staging(staged_input, item_id):
    run, _, _, items, identities = staged_input
    items[0]["id"] = item_id
    with pytest.raises(ValueError, match="alias ID"):
        _materialize_boltz_input_aliases(run, items, identities)
