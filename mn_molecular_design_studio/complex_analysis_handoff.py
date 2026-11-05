"""Results-backed analysis inputs; no molecular identities or curation writes.

Like Chaperone's append-style structure histories, each prediction remains a
separate entry. Studio candidate/variant identities remain authoritative.
"""
from __future__ import annotations

import copy
import hashlib
from pathlib import Path
import shutil

from mn_molecular_design_studio.engine.workflow_engine.core.jobs import write_json


def build_collection(parent, node_id, source, params):
    from mn_molecular_design_studio import candidate_sets as sets
    from mn_molecular_design_studio.validation_acceptance import evaluate_candidate

    mode = params.get("entry_mode", "accepted_only")
    selected = set(params.get("selected_candidate_ids") or [])
    if mode not in {"accepted_only", "manual_selection"}:
        raise ValueError("Choose accepted_only or manual_selection for Complex Analysis.")
    if mode == "manual_selection" and not selected:
        raise ValueError("Manual Complex Analysis requires explicit candidate IDs.")
    destination = Path(parent) / "artifacts/pipeline_nodes" / node_id
    entries, review, receipts = [], [], []
    found = set()
    from collections import Counter
    native_counts = Counter(sets.metadata(row).get('native_candidate_id') for row in sets.iter_candidates(source, review=True))
    for original in sets.iter_candidates(source, review=True):
        row = copy.deepcopy(original)
        raw = row.setdefault("raw_metadata", {})
        info = sets.metadata(row)
        if not info:
            raise ValueError("Connect a Candidate collection before Complex Analysis to preserve source identities.")
        identifier = str(row["candidate_id"])
        found.add(identifier)
        native_id = info.get('native_candidate_id')
        if native_id in selected and native_counts[native_id] != 1:
            raise ValueError('A selected native candidate ID occurs in several source branches. Select its qualified CandidateSet ID instead.')
        if native_id and native_counts[native_id] == 1:
            found.add(native_id)
        retained = raw.get("validation_acceptance") or {}
        # Recheck supplied measurements under the recorded policy. A stale
        # accepted flag, absent policy or missing measurement cannot pass.
        decision = evaluate_candidate(row, retained["policy"]) if retained.get("policy") else {}
        accepted = (retained.get("status") == "accepted" and decision.get("status") == "accepted"
                    and info.get("shared_filter", {}).get("status", "passed") == "passed")
        eligible = accepted if mode == "accepted_only" else identifier in selected or native_id in selected
        if not eligible:
            raw["complex_analysis_route"] = {"status": "excluded", "entry_mode": mode,
                "acceptance_status": decision.get("status", retained.get("status", "not_evaluated"))}
            review.append(row)
            continue
        structures = [(None, row)]
        structures.extend((prediction, prediction["candidate"])
                          for prediction in info.get("predictions", [])
                          if prediction.get("status") == "completed" and isinstance(prediction.get("candidate"), dict))
        seen = set()
        routed = []
        for prediction, model in structures:
            structure = model.get("complex_pdb")
            if not structure:
                continue
            from mn_molecular_design_studio.chaperone_adapter import is_catcif_reference
            model_raw = model.get("raw_metadata") or {}
            interop = model_raw.get("chaperone_interop") or {}
            history = interop.get("structure_history") if isinstance(interop.get("structure_history"), list) else []
            archived = next((item for item in reversed(history)
                             if isinstance(item, dict) and item.get("storage_format") == "catcif"
                             and is_catcif_reference(item.get("path") or "")
                             and (not item.get("name") or item.get("name") == (model.get("source_tool") or model.get("tool")))), None)
            if archived:
                structure = archived["path"]
            archived = is_catcif_reference(structure)
            path = Path(structure) if not archived else None
            if not archived and (path is None or not path.is_file()):
                continue
            if archived:
                digest = hashlib.sha256(structure.encode("utf-8")).hexdigest()
                source_identity = structure
            else:
                with path.open("rb") as stream:
                    digest = hashlib.file_digest(stream, "sha256").hexdigest()
                source_identity = str(path.resolve())
            model_id = str(model.get("model_index") if model.get("model_index") is not None
                           else model_raw.get("model_index", ""))
            prediction_id = str((prediction or {}).get("prediction_id") or info.get("prediction_id") or "")
            key = (prediction_id, model_id, source_identity)
            if key in seen:
                continue
            seen.add(key)
            entry_id = sets.identity(identifier, *key, digest)
            # A separate path per entry preserves distinct model associations
            # even when two methods happen to produce identical coordinates.
            local = destination / "artifacts/structures" / entry_id / (f"{entry_id}.cif" if archived else path.name)
            local.parent.mkdir(parents=True, exist_ok=True)
            if archived:
                from mn_molecular_design_studio.chaperone_adapter import materialize_catcif_structure
                receipt = materialize_catcif_structure(structure, local)
                digest = receipt["sha256"]
            elif not local.exists() or hashlib.sha256(local.read_bytes()).hexdigest() != digest:
                shutil.copyfile(path, local)
            entry = copy.deepcopy(row)
            entry["candidate_id"] = identifier if prediction is None else "analysis_" + entry_id
            for field in ("complex_pdb", "binder_pdb", "target_pdb", "binder_sequence", "target_chains", "binder_chains"):
                if field in model:
                    entry[field] = str(local) if archived and model[field] == structure else copy.deepcopy(model[field])
            entry["complex_pdb"] = str(local)
            entry["metrics"] = copy.deepcopy(model.get("metrics") or {})
            entry.pop("normalized_output", None)
            if model.get("normalized_output"):
                entry["normalized_output"] = copy.deepcopy(model["normalized_output"])
            entry["source_tool"] = model.get("source_tool") or row["source_tool"]
            entry["raw_metadata"].update(copy.deepcopy(model_raw))
            entry["raw_metadata"]["candidate_set"] = copy.deepcopy(info)
            entry["raw_metadata"]["validation_acceptance"] = copy.deepcopy(retained)
            if model_id:
                entry["raw_metadata"]["model_index"] = model_id
            receipt = {"schema": "studio.complex-analysis-entry.v1", "entry_id": entry_id,
                "entry_mode": mode, "acceptance_scope": "candidate", "acceptance": retained,
                "candidate_id": identifier, "native_candidate_id": info.get("native_candidate_id", identifier),
                "variant_id": info.get("variant_id"), "parents": row.get("parents", []),
                "parent_variant_ids": info.get("parent_variant_ids", []), "target": copy.deepcopy(source.get("target") or {}),
                "source_run_dir": info.get("source_run_dir", source["run_dir"]),
                "source_run": info.get("source_run"), "source_structure": source_identity,
                "structure_sha256": digest, "model_id": model_id, "prediction_id": prediction_id,
                "structure_candidate_id": model.get("candidate_id"),
                "method": model.get("source_tool"), "prediction": copy.deepcopy(prediction),
                "source_metadata": copy.deepcopy(model_raw), "candidate_metadata": copy.deepcopy(raw)}
            if archived:
                receipt["structure_storage"] = {"format": "catcif", "source_reference": structure,
                    "materialized_sha256": digest, "materialized_path": str(local)}
            entry.setdefault("raw_metadata", {})["complex_analysis_entry"] = receipt
            entry["raw_metadata"]["complex_analysis_route"] = {"status": "included", "entry_mode": mode}
            # Computational acceptance must never select a molecule for a student.
            routed.append(entry)
            receipts.append(receipt)
        if routed:
            entries.extend(routed)
            review.extend(routed)
        else:
            raw["complex_analysis_route"] = {"status": "unavailable", "reason": "No saved complex structure is available.", "entry_mode": mode}
            review.append(row)
    if mode == 'manual_selection' and selected - found:
        raise ValueError("Selected candidate IDs are absent from the source collection.")
    output = sets.write_collection(parent, node_id, entries, source.get("target") or {})
    review_output = sets.write_collection(parent, node_id + "_review", review, source.get("target") or {})
    output["review_candidates_jsonl"] = review_output["candidates_jsonl"]
    manifest = destination / "artifacts/complex_analysis/entries.json"
    write_json(manifest, {"schema": "studio.complex-analysis-collection.v1", "entry_mode": mode,
        "entries": receipts, "candidate_count": len({r["candidate_id"] for r in receipts}),
        "structure_count": len(receipts), "excluded_count": len(review) - len(entries)})
    output["complex_analysis_handoff"] = str(manifest)
    return output


def inherit_review(parent, node_id, output, source):
    """Keep non-routed evidence visible after an existing analysis tool runs."""
    if not source.get("complex_analysis_handoff") or not output.get("candidates_jsonl"):
        return output
    from mn_molecular_design_studio import candidate_sets as sets
    rows = list(sets.iter_candidates(output))
    rows.extend(row for row in sets.iter_candidates(source, review=True)
                if (row.get("raw_metadata") or {}).get("complex_analysis_route", {}).get("status") != "included")
    review = sets.write_collection(parent, node_id + "_review", rows, source.get("target") or {})
    output["review_candidates_jsonl"] = review["candidates_jsonl"]
    output["complex_analysis_handoff"] = source["complex_analysis_handoff"]
    return output
