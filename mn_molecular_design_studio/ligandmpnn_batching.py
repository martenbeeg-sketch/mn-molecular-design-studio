"""Stage the installed LigandMPNN multi-PDB CLI, preserving native artifacts."""
import re
import shutil
import json
from pathlib import Path

from mn_molecular_design_studio.engine.workflow_engine.core.jobs import write_json, read_json


def group_steps(run_dir, steps):
    # Group only identical global settings. Native multi-input maps carry fixed
    # masks and per-residue biases/omissions. Symmetry, checkpoints and chain
    # selectors remain in the grouping key.
    groups = {}
    for step in steps:
        command, values = list(step["command"]), {}
        for flag in ("--pdb_path", "--out_folder", "--fixed_residues", '--bias_AA_per_residue', '--omit_AA_per_residue'):
            if flag in command:
                i = command.index(flag)
                values[flag] = command[i + 1]
                del command[i:i + 2]
        groups.setdefault(tuple(command), []).append((step, values))
    grouped, mappings = [], []
    root = Path(run_dir) / "artifacts/raw/ligandmpnn"
    for index, (signature, members) in enumerate(groups.items()):
        if len(members) == 1:
            grouped.append(members[0][0])
            continue
        group = root / "groups" / str(index)
        group.mkdir(parents=True, exist_ok=True)
        paths = [v["--pdb_path"] for _, v in members]
        write_json(group / "inputs.json", paths)
        write_json(group / "fixed.json", {v["--pdb_path"]: v.get("--fixed_residues", "") for _, v in members})
        container = "/work/" + group.relative_to(run_dir).as_posix()
        command = [*signature, "--pdb_path_multi", container + "/inputs.json",
                   "--fixed_residues_multi", container + "/fixed.json", "--out_folder", container + "/output"]
        for flag in ('--bias_AA_per_residue', '--omit_AA_per_residue'):
            if any(flag in value for _, value in members):
                maps = {value['--pdb_path']: json.loads((Path(run_dir) / value[flag].removeprefix('/work/')).read_text()) if flag in value else {}
                        for _, value in members}
                name = flag.removeprefix('--') + '.json'
                write_json(group / name, maps)
                command.extend([flag + '_multi', container + '/' + name])
        ids = [Path(v["--pdb_path"]).stem for _, v in members]
        grouped.append(dict(name=f"multi-pdb-{index}", command=command, candidate_ids=ids))
        mappings.append(dict(output=str(group / "output"), candidate_ids=ids))
    write_json(root / "native_groups.json", mappings)
    return grouped, mappings


def materialize(run_dir, mappings):
    # Exact filename boundaries prevent c1 from claiming c1_2's output.
    root = Path(run_dir) / "artifacts/raw/ligandmpnn/output"
    for group in mappings:
        for source_id in group["candidate_ids"]:
            patterns = {"seqs": re.escape(source_id) + r"\.fa",
                        "stats": re.escape(source_id) + r"\.pt",
                        "backbones": re.escape(source_id) + r"_\d+\.pdb",
                        "packed": re.escape(source_id) + r"_packed_\d+_\d+\.pdb"}
            for folder, pattern in patterns.items():
                destination = root / source_id / folder
                destination.mkdir(parents=True, exist_ok=True)
                for path in (Path(group["output"]) / folder).glob("*"):
                    if re.fullmatch(pattern, path.name):
                        shutil.copy2(path, destination / path.name)


def recover(run_dir):
    """Recover verified completed pairs after worker interruption, without launch."""
    from mn_molecular_design_studio.engine.workflow_engine.workflows import sequence_design as native
    from mn_molecular_design_studio.engine.workflow_engine.core.candidates import read_candidates
    root = Path(run_dir) / 'artifacts/raw/ligandmpnn'
    params_file = root / 'ligandmpnn_params.json'
    if not params_file.is_file():
        return []
    params = read_json(params_file)
    params['recover_completed_only'] = True
    mappings = root / 'native_groups.json'
    if mappings.is_file():
        materialize(run_dir, read_json(mappings))
    sources = read_candidates(Path(params['candidates_jsonl']))
    sources = [row for row in sources if any((root / 'output' / row['candidate_id'] / 'backbones').glob('*.pdb'))]
    inputs = {row['candidate_id']: root / 'inputs' / (row['candidate_id'] + '.pdb') for row in sources}
    return native._normalize_ligandmpnn_candidates(Path(run_dir), sources, inputs, params) if sources else []
