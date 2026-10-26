"""Container entry point: native Kejia threading with PDB/JSON handoff.

The native binary silent file and scorefile are retained. PDBs are captured at
the native add2silent boundary so downstream Studio tools need no silent_tools.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace


def run(request: dict, output: Path) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    output = output.resolve()
    native_dir = Path("/opt/kejia/threading")
    sys.path.insert(0, str(native_dir))
    spec = importlib.util.spec_from_file_location("kejia_threading", native_dir / "thread_peptide_sequence_new.py")
    native = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(native)
    records = []
    rejected = 0
    original_add = native.add2silent

    def add_with_pdb(tag, pose, scores, silent_data, silent_file):
        nonlocal rejected
        # Upstream records max_farep but does not consistently skip clashes.
        # Apply its documented threshold at the serialization boundary.
        if float(scores["delta_fa_rep"]) > request["max_farep"]:
            rejected += 1
            return
        original_add(tag, pose, scores, silent_data, silent_file)
        path = output / f"{tag}.pdb"
        pose.dump_pdb(str(path))
        records.append({"tag": tag, "pdb": path.name, "scores": dict(scores),
                        "binder_sequence": pose.split_by_chain(1).sequence(),
                        "target_sequence": pose.split_by_chain(2).sequence()})

    native.add2silent = add_with_pdb
    seed = int(request["seed"])
    flags = native.DEFAULT_PYROSETTA_FLAGS + f" -constant_seed -jran {seed}"
    args = SimpleNamespace(
        template=str(Path(request["template"]).resolve()), peptide_header="studio",
        peptide_seq=request["target_sequence"], interface_dist_cutoff=native.DEFAULT_INTERFACE_DIST_CUTOFF,
        min_seq_identity=request["min_seq_identity"], min_blosum=request["min_blosum"],
        max_farep=request["max_farep"], keep_bidentates=request.get("keep_bidentates", False),
        pyrosetta_flags=flags, flex_docking_repeats=request["flex_docking_repeats"],
        output_csv=native.DEFAULT_OUTPUT_CSV,
    )
    previous = Path.cwd()
    try:
        os.chdir(output)
        native.main(args)
    finally:
        os.chdir(previous)
    result = {"schema_version": 1, "records": records, "clash_rejected": rejected,
              "source_revision": "7f65c18cc57b3e72d290c885c540276c7d9c6ab8",
              "adapter_changes": ["Native silent/score outputs retained; PDB/JSON handoff added.",
                                  "Documented max_farep cutoff enforced before candidate serialization."]}
    (output / "threading.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    run(json.loads(arguments.request.read_text()), arguments.output)
