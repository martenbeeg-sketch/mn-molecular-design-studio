"""Exact native RFD4 refolding artifacts, independent of generator coordinates."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from mn_molecular_design_studio import molecule_imports


SCHEMA = "studio.rfd4-representations.v1"


def _oracle_pae_evidence(predicted: Path, refolded_digest: str) -> dict:
    from mn_molecular_design_studio.engine.workflow_engine.workflows.rfd4_refold_evidence import SCHEMA as PAE_SCHEMA

    stem = predicted.stem.removesuffix("_refolded")
    model = predicted.with_name(f"{stem}_oracle_model.cif")
    pae = predicted.with_name(f"{stem}_oracle_pae.json")
    try:
        model = model.resolve(strict=True)
        pae = pae.resolve(strict=True)
        model.relative_to(predicted.parent)
        pae.relative_to(predicted.parent)
        model_digest = hashlib.sha256(model.read_bytes()).hexdigest()
        pae_bytes = pae.read_bytes()
        payload = json.loads(pae_bytes)
    except (OSError, ValueError):
        return {}
    if (not isinstance(payload, dict) or payload.get("schema") != PAE_SCHEMA
            or payload.get("refolded_sha256") != refolded_digest
            or payload.get("oracle_model_sha256") != model_digest):
        return {}
    return {"oracle_model_path": str(model), "oracle_model_sha256": model_digest,
            "pae_path": str(pae), "pae_sha256": hashlib.sha256(pae_bytes).hexdigest()}


def refolded_representation(generated, metrics, *, root=None):
    generated = Path(generated).resolve(strict=True)
    if root is not None:
        generated.relative_to(Path(root).resolve())
    suffix = "_generated.cif"
    if not generated.name.endswith(suffix):
        return None
    example = generated.name.removesuffix(suffix)
    prefix = "oracle_ppi" if any(key.startswith("oracle_ppi.") for key in metrics) else "oracle_dna" if any(
        key.startswith("oracle_dna.") for key in metrics) else "oracle" if any(key.startswith("oracle.") for key in metrics) else ""
    if not prefix:
        return None
    oracle = {"oracle_ppi": "of3_ppi", "oracle_dna": "of3_dna", "oracle": "of3"}[prefix]
    predicted = generated.with_name(example + "_" + oracle + "_refolded.cif")
    confidences = generated.with_name(example + "_" + oracle + "_confidences.json")
    if not predicted.is_file() or not confidences.is_file():
        return None
    predicted = predicted.resolve(strict=True)
    confidences = confidences.resolve(strict=True)
    predicted.relative_to(generated.parent)
    confidences.relative_to(generated.parent)
    original_bytes, predicted_bytes = generated.read_bytes(), predicted.read_bytes()
    original = molecule_imports.parse_structure(original_bytes, generated.name, require_model_selection=True)[0]
    refolded = molecule_imports.parse_structure(predicted_bytes, predicted.name, require_model_selection=True)[0]
    if molecule_imports.canonical_components(original["components"])[0] != molecule_imports.canonical_components(refolded["components"])[0]:
        raise ValueError("RFD4 refolding changed molecular or chain identities; the models cannot be linked automatically.")
    recorded = json.loads(confidences.read_text())
    checked = []
    for field, value in recorded.items():
        expected = metrics.get(prefix + "." + field)
        if isinstance(value, (int, float)) and isinstance(expected, (int, float)):
            if not math.isfinite(value) or not math.isfinite(expected) or not math.isclose(value, expected, abs_tol=0.011, rel_tol=1e-6):
                raise ValueError("The RFD4 refolding confidence receipt disagrees with this candidate's measurements.")
            checked.append(field)
    if not checked:
        raise ValueError("No candidate-specific RFD4 refolding measurements verify this model.")
    refolded_digest = hashlib.sha256(predicted_bytes).hexdigest()
    evidence = _oracle_pae_evidence(predicted, refolded_digest)
    return {"schema": SCHEMA, "role": "refolded_validation", "oracle": oracle,
            "generated_path": str(generated), "generated_sha256": hashlib.sha256(original_bytes).hexdigest(),
            "refolded_path": str(predicted), "refolded_sha256": refolded_digest,
            **({"oracle_pae": evidence} if evidence else {}),
            "confidence_path": str(confidences), "confidence_sha256": hashlib.sha256(confidences.read_bytes()).hexdigest(),
            "matched_measurements": sorted(checked), "b_factor_interpretation": "plddt_0_100",
            "confidence_basis": "Native RFD4 oracle remapped per-atom pLDDT to the exact refolded model; matching confidence receipt."}


def expects_refolded_model(record):
    metadata = record.get("metadata") or {}
    options = (metadata.get("reference_constraints") or {}).get("options") or {}
    return metadata.get("source_tool") == "rfd4_proteina" and (options.get("refold") is True or bool(recorded_oracle_metrics(metadata)))


def recorded_oracle_metrics(metadata):
    expected = {str(metadata.get(field) or "") for field in ("source_candidate_id", "native_candidate_id", "candidate_id")} - {""}
    metrics = {}
    for field, measurement in (metadata.get("validation_acceptance") or {}).get("measurements", {}).items():
        if not expected or measurement.get("measurement_candidate_id") not in expected or measurement.get("status") != "measured":
            continue
        provenance = (measurement.get("output_provenance") or {}).get(field.removeprefix("common.")) or {}
        for alias in provenance.get("aliases") or []:
            source = str(alias.get("source_field") or "")
            value = alias.get("raw_value")
            if source.startswith(("oracle_ppi.", "oracle_dna.", "oracle.")) and isinstance(value, (int, float)) and not isinstance(value, bool):
                if source in metrics and metrics[source] != value:
                    raise ValueError("Conflicting candidate-specific native oracle measurements.")
                metrics[source] = value
    return metrics


def recover_saved_refold(database, user, original, *, apply=False, artifact_root=None, artifact_generated=None):
    from mn_molecular_design_studio.data_browser import load_saved_record
    from mn_molecular_design_studio.binderscope import saved_chain_roles

    record = load_saved_record(database, user, original)
    if record is None or not expects_refolded_model(record):
        return {"status": "not_applicable"}
    if original.get("eligibility") != "Passed" and (original.get("receipt") or {}).get("state") != "Passed":
        return {"status": "no_recorded_pass"}
    metadata = record["metadata"]
    generated = Path(str(artifact_generated or metadata.get("native_structure") or ""))
    if not generated.is_file():
        return {"status": "missing_generated_source"}
    digest = hashlib.sha256(generated.read_bytes()).hexdigest()
    structures = [row for row in record["structures"] if hashlib.sha256(database.structure_file(user, row["id"]).read_bytes()).hexdigest() == digest]
    if len(structures) != 1:
        return {"status": "ambiguous_generated_representation"}
    old = record["observations"][0] if record["observations"] else {}
    metrics = {**recorded_oracle_metrics(metadata), **(old.get("metrics") or {})}
    representation = refolded_representation(generated, metrics, root=artifact_root)
    if representation is None:
        return {"status": "missing_verified_refold"}
    generated_structure = structures[0]
    parsed = molecule_imports.parse_structure(generated.read_bytes(), generated.name, require_model_selection=True)[0]
    binder, target, role_basis = saved_chain_roles(database, user, record, generated_structure, parsed["components"])
    result = {"status": "verified", "name": record["name"], "group_id": original["group_id"],
              "molecule_id": original["molecule_id"], "representation": representation}
    if not apply:
        return result
    database.group(original["group_id"], user, write=True)
    context = {**(generated_structure.get("context") or {}), **(old.get("context") or {})}
    for key in ("portable", "portable_sha256", "sha256", "content_sha256", "validated_content_sha256",
                "representation_id", "structure_digest", "calculation_structure_digest",
                "structure_transfer_host", "structure_transfer_path", "structure_availability"):
        context.pop(key, None)
    candidate = dict(context.get("candidate") or {})
    candidate_id = original.get("candidate") or candidate.get("candidate_id") or metadata.get("native_candidate_id") or metadata.get("source_candidate_id") or record["id"]
    candidate.update(candidate_id=candidate_id, source_tool="rfd4_proteina", stage="generation.backbone_sequence",
                     binder_chains=binder, target_chains=target, complex_pdb=representation["refolded_path"], binder_pdb=None)
    raw = {**(candidate.get("raw_metadata") or {}), "design_representation": representation}
    validation = metadata.get("validation_acceptance") or raw.get("validation_acceptance") or {}
    if validation:
        raw["validation_acceptance"] = {key: value for key, value in validation.items() if key != "measurements"}
        raw["validation_acceptance"]["original_decision"] = {
            "group_id": original["group_id"], "molecule_id": original["molecule_id"], "observation_id": original.get("observation_id") or ""}
    if metadata.get("distributed_campaign"):
        raw["distributed_campaign"] = metadata["distributed_campaign"]
    candidate["raw_metadata"] = raw
    context.update(candidate_id=candidate_id, candidate=candidate, binder_chains=binder, target_chains=target,
                   design_representation_role="refolded_validation", b_factor_interpretation="plddt_0_100",
                   source_structure_id=generated_structure["id"], source_structure_digest=digest,
                   chain_role_evidence=role_basis,
                   target_molecule_id=original.get("evaluation", {}).get("target_id") or metadata.get("design_target_molecule_id") or "",
                   tool="rfd4_proteina", task_group=(record.get("origin") or {}).get("task_group") or "design",
                   run_id=original["run_id"], status="completed")
    item = {"name": Path(representation["refolded_path"]).name,
            "content": Path(representation["refolded_path"]).read_bytes(), "label": "Refolded validation complex",
            "source_file_field": "complex_pdb", "design_representation_role": "refolded_validation",
            "b_factor_interpretation": "plddt_0_100", "binder_chains": binder, "target_chains": target,
            "validated_content_sha256": representation["refolded_sha256"], "dimensions": 3}
    structure_id = database.attach_structure(original["group_id"], original["molecule_id"], item, context={
        key: value for key, value in context.items() if key not in {"candidate", "representation_id", "structure_digest", "calculation_structure_digest"}})
    database.make_structure_portable(user, structure_id)
    context.update(representation_id=structure_id, structure_digest=representation["refolded_sha256"],
                   calculation_structure_digest=representation["refolded_sha256"], is_representative=True)
    confidence = json.loads(Path(representation["confidence_path"]).read_text())
    prefix = {"of3_ppi": "oracle_ppi", "of3_dna": "oracle_dna", "of3": "oracle"}[representation["oracle"]]
    metrics = {prefix + "." + key: value for key, value in confidence.items()
               if isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)}
    for key, value in (old.get("metrics") or {}).items():
        if key.startswith("refolded_"):
            metrics[key] = value
    observation_id = database.add_observation(original["group_id"], original["molecule_id"], "rfd4_proteina", original["run_id"], metrics,
        context=context, source_key="rfd4-refold-v1:" + (original.get("observation_id") or original["id"]))
    return {**result, "status": "stored", "structure_id": structure_id, "observation_id": observation_id}
