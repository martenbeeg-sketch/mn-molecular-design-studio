"""Automatic collection semantics shared by tiles, APIs and pipeline workers.

A collection contains independent calculations. Chains within one task remain
one complex; they are never inferred to be independent from their chain count.
"""
from __future__ import annotations

import copy
from pathlib import Path

MAX_PREDICTION_ITEMS = 1024
PREDICTION_INPUTS = frozenset({
    "protein_sequences", "target_pdb_text", "target_chain_id", "target_msa_a3m",
    "binder_sequence", "dna_sequences", "rna_sequences", "ligand_ccd_codes",
    "ligand_smiles", "tasks", "compound_set_file", "sequence_file",
})


def parameter_items(value):
    """Accept one input or an ordered selection, without a batch-mode switch."""
    items = value if isinstance(value, list) else [value]
    if not items or len(items) > 100_000 or any(not isinstance(item, dict) for item in items):
        raise ValueError("parameters must be an object or a non-empty list of input objects (at most 100,000).")
    return items


def sequence_tasks(text):
    """FASTA records are independent inputs; ':' keeps chains in one record."""
    records, name, lines = [], None, []
    for raw in str(text).splitlines():
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            if name is not None:
                records.append({"name": name, "protein_sequences": "".join(lines)})
            name, lines = line[1:].strip() or f"sequence_{len(records) + 1}", []
        elif name is None:
            # A plain sequence list has one independent sequence per line.
            records.append({"name": f"sequence_{len(records) + 1}", "protein_sequences": line})
        else:
            lines.append(line)
    if name is not None:
        records.append({"name": name, "protein_sequences": "".join(lines)})
    if not records or any(not row["protein_sequences"] for row in records):
        raise ValueError("The sequence selection contains an empty sequence.")
    return records


def prediction_collection_tasks(shared, *, compound_set_file=None, sequence_file=None):
    """Expand selected compounds/sequences into a single native model batch."""
    from mn_molecular_design_studio.ligand_compounds import load_compound_set
    context = copy.deepcopy(shared)
    tasks = [context]
    if sequence_file:
        tasks = [{**context, **row} for row in sequence_tasks(Path(sequence_file).read_text())]
    if compound_set_file:
        compounds = load_compound_set(Path(compound_set_file))["compounds"]
        if len(tasks) * len(compounds) > MAX_PREDICTION_ITEMS:
            raise ValueError(f"A prediction selection may contain at most {MAX_PREDICTION_ITEMS} calculations; narrow the selection.")
        if context.get("ligand_smiles") or context.get("ligand_ccd_codes"):
            raise ValueError("Use the selected compound collection or ligand chains in the complex, not both.")
        tasks = [
            {**task, "name": f"{task.get('name', 'input')}_{row['compound_id']}" if len(tasks) > 1 else row["compound_id"],
             "ligand_smiles": row["smiles"], "ligand_ccd_codes": "", "source_compound_id": row["compound_id"]}
            for task in tasks for row in compounds
        ]
    if len(tasks) > MAX_PREDICTION_ITEMS:
        raise ValueError(f"A prediction selection may contain at most {MAX_PREDICTION_ITEMS} calculations; narrow the selection.")
    return tasks


def coalesce_prediction_parameters(workflow_id, items):
    """Use one model/container invocation when selected inputs share settings."""
    if workflow_id != "refolding.optimized_prediction" or len(items) < 2:
        return None
    def settings(item):
        return {key: value for key, value in item.items() if key not in PREDICTION_INPUTS | {"job_name"}}
    if any(settings(item) != settings(items[0]) for item in items[1:]):
        return None
    if any(item.get("validation_sources") is not None for item in items):
        return None
    combined = {key: value for key, value in items[0].items() if key not in PREDICTION_INPUTS}
    tasks = []
    for index, item in enumerate(items, 1):
        if item.get("tasks") is not None:
            tasks.extend(copy.deepcopy(item["tasks"]))
        else:
            context = {key: value for key, value in item.items() if key in PREDICTION_INPUTS - {"compound_set_file", "sequence_file", "tasks"}}
            context["name"] = item.get("job_name") or f"input_{index}"
            tasks.extend(prediction_collection_tasks(context, compound_set_file=item.get("compound_set_file"), sequence_file=item.get("sequence_file")))
    combined["tasks"] = tasks
    return combined


def describe_collection(inputs, params):
    """Describe native collection inputs without creating child scheduler jobs."""
    values = {**(inputs or {}), **(params or {})}
    for key, kind in (("tasks", "prediction"), ("batch_inputs", "prediction"), ("items", "calculation"),
                      ("target_entries", "target"), ("validation_sources", "design")):
        rows = values.get(key)
        if isinstance(rows, list) and rows:
            if key == "validation_sources":
                ids = [str(identifier) for row in rows for identifier in row.get("candidate_ids", [])]
            else:
                ids = [str(row.get("id") or row.get("name") or row.get("compound_id") or f"input_{i}") if isinstance(row, dict) else str(row)
                       for i, row in enumerate(rows, 1)]
            if ids:
                return {"kind": kind, "item_count": len(ids), "item_ids": ids}
    path = values.get("compound_set_file")
    if path:
        from mn_molecular_design_studio.ligand_compounds import load_compound_set
        rows = load_compound_set(Path(path))["compounds"]
        limit = values.get("compound_limit", 0)
        rows = rows[:limit] if limit else rows
        return {"kind": "compound", "item_count": len(rows), "item_ids": [row["compound_id"] for row in rows]}
    path = values.get("candidates_jsonl")
    if path:
        from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
        rows = read_candidates(Path(path))
        selected = values.get("selected_candidate_ids") or values.get("candidate_ids")
        ids = [row["candidate_id"] for row in rows if selected is None or row["candidate_id"] in selected]
        return {"kind": "candidate", "item_count": len(ids), "item_ids": ids}
    path = values.get("sequence_file")
    if path:
        rows = sequence_tasks(Path(path).read_text())
        return {"kind": "sequence", "item_count": len(rows), "item_ids": [row["name"] for row in rows]}
    return None


def record_collection(run_dir, inputs, params):
    """Keep batch identity on the normal run record and in a portable manifest."""
    from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json, write_json
    try:
        description = describe_collection(inputs, params)
    except (OSError, ValueError, KeyError, TypeError):
        # Input validation reports unreadable selections; job creation must not
        # fail after it has already made the durable run directory.
        return
    if not description:
        return
    batch = {"schema_version": 1, "automatic": True, "execution": "native", **description}
    write_json(Path(run_dir) / "batch_manifest.json", batch)
    metadata = read_json(Path(run_dir) / "metadata.json")
    metadata["batch"] = {key: value for key, value in batch.items() if key != "item_ids"}
    write_json(Path(run_dir) / "metadata.json", metadata)
