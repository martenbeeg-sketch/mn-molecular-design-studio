"""Optional lossless bridge between Studio candidate rows and Chaperone.

Chaperone is a data-handler/filter backend only. Studio candidate IDs, acceptance
policy and complete normalized payload remain authoritative. Chaperone and its
scientific dependencies are imported only when this adapter is selected.
"""
from __future__ import annotations

import copy
import base64
import hashlib
import json
import os
import subprocess
import sys
from collections.abc import Iterable
from pathlib import Path


ADAPTER_VERSION = "studio-chaperone.v1"


def worker_configuration() -> dict[str, str]:
    """Resolve an explicitly managed Chaperone worker without UI path setup.

    Environment variables remain useful for isolated tests and operators. A
    one-time setup script may instead write ``chaperone-worker.json`` under
    the Studio app home; normal users then select Chaperone on the canvas
    without knowing which Python interpreter or checkout provides it.
    """
    from mn_molecular_design_studio.runtime import app_home

    config_path = Path(os.environ.get("MDS_CHAPERONE_CONFIG", "").strip() or
                       (app_home() / "chaperone-worker.json"))
    try:
        saved = json.loads(config_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        saved = {}
    if not isinstance(saved, dict):
        saved = {}
    return {
        "python": os.environ.get("MDS_CHAPERONE_PYTHON", "").strip() or str(saved.get("python") or "").strip(),
        "source": os.environ.get("MDS_CHAPERONE_SRC", "").strip() or str(saved.get("source") or "").strip(),
        "blocks": os.environ.get("MDS_CHAPERONE_BLOCKS", "").strip() or str(saved.get("blocks") or "").strip(),
        "source_revision": str(saved.get("source_revision") or "").strip(),
    }


def is_catcif_reference(path: str | Path) -> bool:
    value = str(path or "")
    split = value.find(":")
    return split > 0 and value[:split].endswith(".catcif") and bool(value[split + 1:])


def materialize_catcif_structure(reference: str | Path, destination: str | Path) -> dict:
    """Resolve one Chaperone CatCIF entry to a local, content-addressed input.

    Resolution runs in Chaperone's configured interpreter when the Studio app
    and Chaperone use different Python environments. The archive itself stays
    untouched; the returned receipt retains the entry reference and digest.
    """
    reference = str(reference)
    if not is_catcif_reference(reference):
        raise ValueError("Expected a CatCIF reference in archive.catcif:entry_id form.")
    config = worker_configuration()
    worker = config["python"]
    if worker and Path(worker).resolve() != Path(sys.executable).resolve():
        repo_root = str(Path(__file__).resolve().parents[1])
        chaperone_src = config["source"]
        chaperone_blocks = config["blocks"]
        environment = {**os.environ, "PYTHONPATH": os.pathsep.join(
            [repo_root, *[item for item in (chaperone_src, chaperone_blocks) if item]])}
        process = subprocess.run(
            [worker, "-m", "mn_molecular_design_studio.chaperone_adapter", "--worker"],
            input=json.dumps({"operation": "read_catcif", "reference": reference}),
            text=True, capture_output=True, check=False, env=environment,
        )
        if process.returncode:
            detail = process.stderr.strip().splitlines()[-1] if process.stderr.strip() else "unknown CatCIF reader error"
            raise RuntimeError(f"Chaperone CatCIF materialization failed: {detail}")
        try:
            payload = json.loads(process.stdout)
            content = base64.b64decode(payload["content_base64"], validate=True)
        except (ValueError, KeyError, json.JSONDecodeError) as exc:
            raise RuntimeError("Chaperone CatCIF reader returned an invalid payload.") from exc
    else:
        try:
            from chaperone.data.catcif import read_catcif_cif_bytes
            content = read_catcif_cif_bytes(reference)
        except ImportError as exc:
            raise RuntimeError("CatCIF retrieval requires Chaperone and catcif_tools.") from exc
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(content)
    return {"source_reference": reference, "materialized_path": str(destination),
            "sha256": hashlib.sha256(content).hexdigest(), "size_bytes": len(content)}


def _chaperone_types():
    try:
        from chaperone.data import PipelineDataHandler
        from chaperone.data.schema import make_row, make_sequence_row, get_sequences_from_row
    except ImportError as exc:
        raise RuntimeError(
            "Chaperone interoperability requires Chaperone and its Python 3.12 dependencies."
        ) from exc
    return PipelineDataHandler, make_row, make_sequence_row, get_sequences_from_row


def to_pipeline_data_handler(candidates: Iterable[dict], *, filter_values: dict[str, dict] | None = None):
    """Map Studio candidates to the real Chaperone schema with a full sidecar.

    Returns ``(handler, sidecar)``. The sidecar stores exact source records so
    arbitrary Studio metadata and prediction/model multiplicity survive the
    narrower Chaperone structure schema unchanged.
    """
    PipelineDataHandler, make_row, make_sequence_row, _ = _chaperone_types()
    rows, sidecar = [], {}
    for candidate in candidates:
        record = copy.deepcopy(candidate)
        candidate_id = str(record.get("candidate_id") or "")
        if not candidate_id:
            raise ValueError("Every Studio candidate needs a candidate_id for Chaperone roundtrip.")
        if candidate_id in sidecar:
            raise ValueError(f"Duplicate Studio candidate_id cannot be represented: {candidate_id}")
        raw = record.get("raw_metadata") if isinstance(record.get("raw_metadata"), dict) else {}
        lineage = {
            "adapter": ADAPTER_VERSION,
            "candidate_id": candidate_id,
            "parents": list(record.get("parents") or []),
            "source_tool": record.get("source_tool") or record.get("tool") or "",
            "target_chains": list(record.get("target_chains") or []),
            "binder_chains": list(record.get("binder_chains") or []),
            "model_id": record.get("model_id") or record.get("model_index") or raw.get("model_id") or raw.get("model_index"),
            "method": record.get("source_tool") or record.get("tool") or "",
        }
        metrics = {"studio": copy.deepcopy(record.get("metrics") or {}),
                   "_studio_interop": lineage}
        if filter_values and candidate_id in filter_values:
            # A private projection lets Chaperone's native FilterBlock query
            # canonical or stage-qualified Studio measurements without
            # changing the original Studio metrics or sidecar record.
            metrics["studio_filter"] = copy.deepcopy(filter_values[candidate_id])
        path = record.get("complex_pdb") or record.get("binder_pdb") or record.get("target_pdb")
        sequence = record.get("binder_sequence")
        binder_chains = list(record.get("binder_chains") or [])
        chain_sequences = raw.get("chain_sequences") or raw.get("sequences_by_chain")
        sequences = ({str(chain): str(seq) for chain, seq in chain_sequences.items()
                     if seq} if isinstance(chain_sequences, dict) else {})
        if sequence and len(binder_chains) == 1:
            sequences.setdefault(str(binder_chains[0]), str(sequence))
        target_sequences = raw.get("target_sequences")
        if isinstance(target_sequences, dict):
            sequences.update({str(chain): str(seq) for chain, seq in target_sequences.items() if seq})
        elif target_sequences and len(record.get("target_chains") or []) == 1:
            sequences.setdefault(str(record["target_chains"][0]), str(target_sequences))
        parent_id = (record.get("parents") or [None])[0]
        if path:
            row = make_row(name=str(lineage["source_tool"] or "studio"), path=path,
                           metadata={"task": "studio_candidate"}, metrics=metrics,
                           design_id=candidate_id, parent_design_id=parent_id,
                           sequences=sequences)
        elif sequences:
            row = make_sequence_row(design_id=candidate_id,
                                    sequences=sequences,
                                    parent_design_id=parent_id, metrics=metrics)
        else:
            row = make_row(name=str(lineage["source_tool"] or "studio"), path=None,
                           metrics=metrics, design_id=candidate_id,
                           parent_design_id=parent_id)
        rows.append(row)
        sidecar[candidate_id] = record
    return PipelineDataHandler(rows, allow_sequences=True), sidecar


def from_pipeline_data_handler(handler, sidecar: dict[str, dict]) -> list[dict]:
    """Restore exact Studio rows and add explicit Chaperone lineage receipts."""
    from chaperone.data.schema import Column, StructureColumn, get_deserialized, get_sequences_from_row

    output = []
    seen = set()
    for row in handler._df.iter_rows(named=True):
        candidate_id = str(row.get(str(Column.DESIGN_ID)) or "")
        if candidate_id not in sidecar:
            raise ValueError(f"Chaperone returned an unknown Studio design_id: {candidate_id!r}")
        if candidate_id in seen:
            raise ValueError(f"Chaperone returned duplicate Studio design_id: {candidate_id}")
        seen.add(candidate_id)
        candidate = copy.deepcopy(sidecar[candidate_id])
        metric_payload = get_deserialized(row, Column.METRICS)
        receipt = metric_payload.get("_studio_interop") or {}
        structures = row.get(str(Column.STRUCTURES)) or []
        candidate.setdefault("raw_metadata", {})["chaperone_interop"] = {
            "adapter": ADAPTER_VERSION,
            "design_id": candidate_id,
            "parent_design_id": row.get(str(Column.PARENT_DESIGN_ID)),
            "artifact_kind": row.get(str(Column.ARTIFACT_KIND)),
            "source_tool": receipt.get("source_tool"),
            "method": receipt.get("method"),
            "model_id": receipt.get("model_id"),
            "target_chains": receipt.get("target_chains", []),
            "binder_chains": receipt.get("binder_chains", []),
            "handler_metrics": {key: value for key, value in metric_payload.items()
                                if key != "_studio_interop"},
            # Keep the complete real handler row in the receipt so fields the
            # Studio adapter does not interpret are still auditable.
            "chaperone_record": copy.deepcopy(row),
            "chain_sequences": get_sequences_from_row(row),
            "structure_history": [
                {"name": item.get(str(StructureColumn.NAME)),
                 "path": item.get(str(StructureColumn.PATH)),
                 "storage_format": item.get(str(StructureColumn.STORAGE_FORMAT))}
                for item in structures
            ],
        }
        output.append(candidate)
    return output


def roundtrip_candidates(candidates: Iterable[dict], *, filter_expression: str | None = None,
                         filter_values: dict[str, dict] | None = None,
                         filter_authority: str = "auxiliary; Studio acceptance and routing are unchanged") -> list[dict]:
    """Pass Studio candidates through Chaperone's handler and optional filter.

    Filter results are attached as auxiliary interoperability evidence while
    every candidate is retained. They do not alter Studio acceptance or routing.
    """
    records = list(candidates)
    worker = worker_configuration()["python"]
    if worker and Path(worker).resolve() != Path(sys.executable).resolve():
        return _roundtrip_in_worker(records, filter_expression, filter_values, filter_authority, worker)
    return _roundtrip_local(records, filter_expression, filter_values, filter_authority)


def _roundtrip_local(records: list[dict], filter_expression: str | None,
                     filter_values: dict[str, dict] | None = None,
                     filter_authority: str = "auxiliary; Studio acceptance and routing are unchanged") -> list[dict]:
    handler, sidecar = to_pipeline_data_handler(records, filter_values=filter_values)
    pass_ids = None
    if filter_expression and filter_expression.strip():
        from blocks.util_blocks.filter import FilterBlock

        filtered = FilterBlock(name="studio_filter", expr=filter_expression.strip())(handler)
        pass_ids = {str(value) for value in filtered._df["design_id"].to_list()}
        # Keep the unfiltered handler for the roundtrip; the FilterBlock is an
        # auxiliary gate report and must not erase rejected/unresolved history.
    output = from_pipeline_data_handler(handler, sidecar)
    if pass_ids is not None:
        for candidate in output:
            receipt = candidate["raw_metadata"]["chaperone_interop"]
            receipt["filter"] = {
                "expression": filter_expression.strip(),
                "status": "passed" if candidate["candidate_id"] in pass_ids else "not_passed",
                "authority": filter_authority,
            }
    return output


def _roundtrip_in_worker(records: list[dict], filter_expression: str | None,
                         filter_values: dict[str, dict] | None,
                         filter_authority: str, python: str) -> list[dict]:
    """Run Chaperone in its configured interpreter without adding it to Studio."""
    repo_root = str(Path(__file__).resolve().parents[1])
    config = worker_configuration()
    chaperone_src = config["source"]
    chaperone_blocks = config["blocks"]
    python_path = [repo_root, *[item for item in (chaperone_src, chaperone_blocks) if item]]
    environment = {**os.environ, "PYTHONPATH": os.pathsep.join(python_path)}
    process = subprocess.run(
        [python, "-m", "mn_molecular_design_studio.chaperone_adapter", "--worker"],
        input=json.dumps({"candidates": records, "filter_expression": filter_expression,
                          "filter_values": filter_values, "filter_authority": filter_authority}, allow_nan=False),
        text=True, capture_output=True, check=False, env=environment,
    )
    if process.returncode:
        detail = process.stderr.strip().splitlines()[-1] if process.stderr.strip() else "unknown worker error"
        raise RuntimeError(f"Chaperone interoperability worker failed: {detail}")
    try:
        return json.loads(process.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError("Chaperone interoperability worker returned invalid JSON.") from exc


def filter_candidates(candidates: Iterable[dict], expression: str) -> tuple[list[dict], list[dict]]:
    """Run Chaperone's actual FilterBlock and return (passing rows, decisions).

    This low-level bridge reports Chaperone's metric-expression result. It does
    not replace Studio's scientific acceptance or review decision.
    """
    if not isinstance(expression, str) or not expression.strip():
        raise ValueError("A Chaperone filter expression is required.")
    rows = roundtrip_candidates(candidates, filter_expression=expression)
    decisions = [{"candidate_id": candidate["candidate_id"],
                  "status": candidate["raw_metadata"]["chaperone_interop"]["filter"]["status"],
                  "filter_expression": expression}
                 for candidate in rows]
    return [candidate for candidate in rows
            if candidate["raw_metadata"]["chaperone_interop"]["filter"]["status"] == "passed"], decisions


def apply_filter_block(candidates: Iterable[dict], settings: dict) -> dict:
    """Apply Chaperone's native FilterBlock to metrics declared in Studio rules.

    The adapter projects resolved numeric values into a temporary Chaperone
    metric namespace; original CandidateSet rows remain unchanged. Candidates
    with unavailable metrics stay in the Studio decision table and are not
    silently converted into passes.
    """
    from mn_molecular_design_studio.candidate_selection import _observation, selection_policy

    records = list(candidates)
    policy = selection_policy(settings)
    rules = policy["rules"]
    if not rules:
        return {"status": "skipped", "reason": "No numeric metric rules are configured."}
    if any(str(rule["field"]).startswith("evidence/") for rule in rules):
        return {"status": "skipped", "reason": "Stage-qualified evidence selectors are evaluated by Studio's provenance-aware resolver."}
    fields = list(dict.fromkeys(rule["field"] for rule in rules))
    values: dict[str, dict] = {}
    for row in records:
        candidate_id = str(row.get("candidate_id") or "")
        for field in fields:
            value, status, _, _ = _observation(row, field)
            if status == "valid":
                values.setdefault(candidate_id, {})[field] = value
    unavailable = [field for field in fields if not any(field in row_values for row_values in values.values())]
    if unavailable:
        return {"status": "skipped", "reason": "Chaperone cannot filter fields with no numeric values: " + ", ".join(unavailable)}

    def comparison(rule):
        key = f"studio_filter.{rule['field']}"
        quoted = f"`{key}`"
        op = rule["operator"]
        lower = float(rule["value"])
        if op == "between":
            upper = float(rule.get("upper_value", rule.get("max_value")))
            return f"({quoted} >= {lower!r} and {quoted} <= {upper!r})"
        return f"{quoted} {op} {lower!r}"

    joiner = " and " if policy["rule_logic"] == "all" else " or "
    expression = joiner.join(comparison(rule) for rule in rules)
    enriched = roundtrip_candidates(records, filter_expression=expression,
        filter_values=values,
        filter_authority="Chaperone FilterBlock metric pass; Studio retains full decisions and provenance")
    statuses = {row["candidate_id"]: row["raw_metadata"]["chaperone_interop"].get("filter", {}).get("status")
                for row in enriched}
    return {"status": "applied", "implementation": "blocks.util_blocks.filter.FilterBlock",
            "expression": expression,
            "passed_candidate_ids": sorted(key for key, status in statuses.items() if status == "passed"),
            "candidate_statuses": statuses, "filtered_candidate_count": sum(status == "passed" for status in statuses.values()),
            "candidate_count": len(records), "logic": policy["rule_logic"], "fields": fields}


def _worker_main() -> None:
    payload = json.load(sys.stdin)
    if payload.get("operation") == "read_catcif":
        from chaperone.data.catcif import read_catcif_cif_bytes
        content = read_catcif_cif_bytes(str(payload["reference"]))
        sys.stdout.write(json.dumps({"content_base64": base64.b64encode(content).decode("ascii")},
                                    separators=(",", ":")))
        return
    rows = _roundtrip_local(payload.get("candidates") or [], payload.get("filter_expression"),
                            payload.get("filter_values"), payload.get("filter_authority") or
                            "auxiliary; Studio acceptance and routing are unchanged")
    sys.stdout.write(json.dumps(rows, allow_nan=False, separators=(",", ":")))


if __name__ == "__main__" and sys.argv[1:] == ["--worker"]:
    _worker_main()
