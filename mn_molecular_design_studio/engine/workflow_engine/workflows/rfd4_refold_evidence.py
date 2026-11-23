"""Retain exact OpenFold3 output pairs before the RFD4 oracle removes its scratch directory."""
from __future__ import annotations

from contextlib import contextmanager
from functools import wraps
from hashlib import sha256
import json
import logging
from pathlib import Path
import runpy
import shutil
import sys


SCHEMA = "studio.rfd4-oracle-pae.v1"


def capture_model_pair(state, refolded):
    if not refolded.is_file():
        return
    source = state.af3_output_dir / f"{state.job_name}_model.cif"
    confidence = state.af3_output_dir / f"{state.job_name}_confidences.json"
    if not (source.is_file() and confidence.is_file()):
        return
    full_confidence = json.loads(confidence.read_text())
    if not isinstance(full_confidence, dict):
        raise ValueError("The OpenFold3 confidence file is not a JSON object.")
    fields = ("pae", "token_chain_ids", "token_res_ids")
    payload = {key: full_confidence[key] for key in fields if key in full_confidence}
    matrix = payload.get("pae")
    tokens = payload.get("token_chain_ids")
    if not isinstance(matrix, list) or not matrix or not isinstance(tokens, list) or len(tokens) != len(matrix):
        raise ValueError("The OpenFold3 PAE matrix is missing its token inventory.")
    if any(not isinstance(row, list) or len(row) != len(matrix) for row in matrix):
        raise ValueError("The OpenFold3 PAE matrix is not square.")
    residues = payload.get("token_res_ids")
    if residues is not None and (not isinstance(residues, list) or len(residues) != len(matrix)):
        raise ValueError("The OpenFold3 PAE residue inventory differs from its matrix.")

    stem = refolded.stem.removesuffix("_refolded")
    model = refolded.with_name(f"{stem}_oracle_model.cif")
    evidence = refolded.with_name(f"{stem}_oracle_pae.json")
    source_bytes = source.read_bytes()
    payload.update(schema=SCHEMA, oracle_model_sha256=sha256(source_bytes).hexdigest(),
                   refolded_sha256=sha256(refolded.read_bytes()).hexdigest())
    staging = evidence.with_suffix(".json.tmp")
    try:
        staging.write_text(json.dumps(payload, separators=(",", ":"), allow_nan=False))
        model.write_bytes(source_bytes)
        staging.replace(evidence)
    finally:
        staging.unlink(missing_ok=True)


@contextmanager
def capture_refolds(oracle_class):
    original = oracle_class._finalize_af3_result

    def written_at(path):
        try:
            return path.stat().st_mtime_ns
        except OSError:
            return None

    @wraps(original)
    def finalize(oracle, state):
        refolded = oracle.refolded_cif_path(state.out_dir, state.example_id, state.save_id)
        previous_write = written_at(refolded)
        result = original(oracle, state)
        if result is not None and written_at(refolded) != previous_write:
            try:
                capture_model_pair(state, refolded)
            except (OSError, ValueError, TypeError) as exc:
                logging.getLogger(__name__).warning("Could not retain OpenFold3 PAE for %s: %s", state.example_id, exc)
        return result

    oracle_class._finalize_af3_result = finalize
    try:
        yield
    finally:
        oracle_class._finalize_af3_result = original


def main():
    from rfproteina.metrics.oracles.af3_oracle import AF3Oracle

    sys.argv = sys.argv[1:]
    entrypoint = shutil.which(sys.argv[0]) or sys.argv[0]
    with capture_refolds(AF3Oracle):
        runpy.run_path(entrypoint, run_name="__main__")


if __name__ == "__main__":
    main()
