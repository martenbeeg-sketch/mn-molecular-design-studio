"""Run TMol scoring and fixed-sequence geometry refinement inside the Studio image."""
from __future__ import annotations

import argparse
import json
import math
import traceback
from collections import Counter
from importlib.metadata import version
from pathlib import Path

import torch
import tmol
from tmol.kinematics import CartesianMoveMap, FoldForest
from tmol.optimization import run_cart_min
from tmol.pack import PackerPalette
from tmol.relax import fast_relax


def heavy_atoms(array):
    insertion_codes = array.ins_code if hasattr(array, "ins_code") else [""] * len(array)
    return Counter(
        (str(chain), int(residue), str(insertion).strip(), str(component), str(atom), str(element).upper())
        for chain, residue, insertion, component, atom, element in zip(
            array.chain_id, array.res_id, insertion_codes, array.res_name, array.atom_name, array.element,
        )
        if element.upper() not in {"H", "D"} and component.upper() not in {"HOH", "WAT", "DOD"}
    )


def score(pose, score_function):
    with torch.no_grad():
        energy = float(score_function.render_whole_pose_scoring_module(pose)(pose.coords).item())
    if not math.isfinite(energy):
        raise ValueError("TMol returned a non-finite energy.")
    return energy


def evaluate(item, config, device):
    source = Path(item["source"])
    array = tmol.atom_array_from_file(source, item["format"])
    chemistry = heavy_atoms(array)
    if not chemistry:
        raise ValueError("The input has no supported heavy atoms.")
    prepare = config["prepare_ligands"] or item["format"] in {"sdf", "mol2"}
    ligand_parameters = [entry["path"] for entry in config.get("ligand_parameters", [])]
    pose, context = tmol.pose_stack_from_atom_array(
        array, device, prepare_ligands=prepare, ligand_params_files=ligand_parameters or None,
        strict_ligands=True, return_context=True,
    )
    score_function = tmol.beta2016_score_function(device, param_db=context.parameter_database)
    before = score(pose, score_function)
    mode = config["mode"]
    if mode == "minimize":
        pose = run_cart_min(pose, score_function)
    elif mode == "fast_relax":
        pose = fast_relax(
            pose, score_function, PackerPalette(), CartesianMoveMap(),
            FoldForest.reasonable_fold_forest(pose), num_repeats=config["repeats"],
        )
    after = before if mode == "score" else score(pose, score_function)
    output = tmol.atom_array_from_pose_stack(pose)
    absent = chemistry - heavy_atoms(output)
    if absent:
        examples = ", ".join(f"{chain}:{residue}{code} {component} {atom}"
                             for chain, residue, code, component, atom, _element in list(absent)[:3])
        raise ValueError(f"TMol could not preserve the input heavy-atom labels ({examples}).")
    waters_in = sum(str(component).upper() in {"HOH", "WAT", "DOD"} and str(element).upper() == "O"
                    for component, element in zip(array.res_name, array.element))
    waters_out = sum(str(component).upper() in {"HOH", "WAT", "DOD"} and str(element).upper() == "O"
                     for component, element in zip(output.res_name, output.element))
    path = Path(item["output"])
    tmol.file_from_pose_stack(pose, path, format="cif")
    return {
        "status": "completed", "tmol_energy_before": before, "tmol_energy_after": after,
        "tmol_energy_change": after - before,
        "ligand_parameterization": prepare, "input_waters": waters_in,
        "output_waters": waters_out, "structure": path.name,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    if not torch.cuda.is_available():
        raise RuntimeError("The TMol worker requires a CUDA GPU.")
    device = torch.device("cuda:0")
    results = []
    for item in config["items"]:
        row = {"id": item["id"], "name": item["name"], "input_sha256": item["input_sha256"]}
        try:
            if item.get("input_error"):
                raise ValueError(item["input_error"])
            row.update(evaluate(item, config, device))
        except Exception as exc:
            traceback.print_exc()
            row.update(status="failed", error=str(exc))
            if isinstance(exc, torch.cuda.OutOfMemoryError):
                torch.cuda.empty_cache()
        results.append(row)
        print(f'{item["id"]} {row["status"]}', flush=True)
    summary = {
        "tmol_version": version("tmol"), "cuda_operators_version": version("tmol-ops-cu13"),
        "gpu": torch.cuda.get_device_name(device), "mode": config["mode"],
        "score_function": "TMol beta2016", "records": results,
    }
    Path(config["summary"]).write_text(json.dumps(summary, indent=2, allow_nan=False), encoding="utf-8")


if __name__ == "__main__":
    main()
