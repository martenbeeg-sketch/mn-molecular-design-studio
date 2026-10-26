"""Port of the paper's threading and interface redesign to a bounded job API."""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

AA3 = dict(zip("ACDEFGHIKLMNPQRSTVWY", (
    "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
    "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR",
)))


def build_protocol_xml(native: str, design_repeats: int, relax_repeats: int) -> str:
    """Retain deposited design movers; external analysis is a downstream step."""
    source = ET.fromstring(native)
    protocol = ET.Element("ROSETTASCRIPTS")
    for tag in ("SCOREFXNS", "RESIDUE_SELECTORS", "TASKOPERATIONS"):
        section = copy.deepcopy(source.find(tag))
        # Remove stray non-XML selector text in the deposit; no selector is added.
        section.text = None
        for element in section:
            element.tail = None
        protocol.append(section)
    movers = ET.SubElement(protocol, "MOVERS")
    names = ("docking_foldtree", "5trp", "alacomp", "FastDesign_ppi", "FastRelax")
    for name in names:
        element = source.find(f"./MOVERS/*[@name='{name}']")
        if element is None:
            raise ValueError(f"The deposited protocol is missing mover {name}.")
        element = copy.deepcopy(element)
        if name in {"FastDesign_ppi", "FastRelax"}:
            element.set("repeats", str(design_repeats if name == "FastDesign_ppi" else relax_repeats))
        movers.append(element)
    steps = ET.SubElement(protocol, "PROTOCOLS")
    for name in names:
        ET.SubElement(steps, "Add", mover_name=name)
    return ET.tostring(protocol, encoding="unicode")


def _rmsd(reference, predicted, start: int, end: int) -> float:
    import numpy as np
    names = ("N", "CA", "C", "O")
    before = np.array([[*reference.residue(i).xyz(atom)] for i in range(start, end + 1) for atom in names])
    after = np.array([[*predicted.residue(i).xyz(atom)] for i in range(start, end + 1) for atom in names])
    before -= before.mean(axis=0)
    after -= after.mean(axis=0)
    u, _, vt = np.linalg.svd(after.T @ before)
    rotation = u @ np.diag([1., 1., np.linalg.det(u @ vt)]) @ vt
    return float(np.sqrt(np.mean(np.sum((after @ rotation - before) ** 2, axis=1))))


def run(request: dict, output: Path, *, validate_only: bool = False) -> dict:
    import faulthandler
    import pyrosetta
    from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue

    output.mkdir(parents=True, exist_ok=True)
    faulthandler.enable()
    print("Initializing PyRosetta and the deposited interface-design protocol", flush=True)
    pyrosetta.init(f"-beta -mute all -constant_seed -jran {int(request['seed'])}")
    xml = build_protocol_xml(Path(request["native_protocol"]).read_text(), request["design_repeats"], request["relax_repeats"])
    (output / "studio_protocol.xml").write_text(xml + "\n")
    objects = XmlObjects.create_from_string(xml)
    if validate_only:
        return {"valid": True}
    print("Loading the staged peptide-trap scaffold", flush=True)
    original = pyrosetta.pose_from_file(request["template"])
    if original.num_chains() != 2 or original.chain_sequence(1) != request["template_peptide_sequence"]:
        raise ValueError("The staged template peptide does not match the curated scaffold.")
    fragment = request["sequence"]
    if len(fragment) != original.chain_end(1):
        raise ValueError("The target fragment and scaffold peptide lengths do not match.")
    threaded = original.clone()
    for index, aa in enumerate(fragment, 1):
        MutateResidue(index, AA3[aa]).apply(threaded)
    threaded.dump_pdb(str(output / "threaded_input.pdb"))
    print("Target sequence threaded; starting binder-interface redesign", flush=True)
    records = []
    for replicate in range(1, request["designs_per_register"] + 1):
        pose = threaded.clone()
        # Reset the deterministic RNG per replicate so one candidate is reproducible.
        pyrosetta.rosetta.numeric.random.rg().set_seed((int(request["seed"]) + replicate - 1) % 2147483647 or 1)
        for name in ("docking_foldtree", "5trp", "alacomp", "FastDesign_ppi", "FastRelax"):
            print(f"Design {replicate}: {name}", flush=True)
            objects.get_mover(name).apply(pose)
        if pose.num_chains() != 2 or pose.chain_sequence(1) != fragment:
            raise ValueError("Rosetta changed the fixed target sequence or the number of chains.")
        if len(pose.chain_sequence(2)) != len(threaded.chain_sequence(2)):
            raise ValueError("Rosetta changed the binder length.")
        peptide_rmsd = _rmsd(threaded, pose, 1, threaded.chain_end(1))
        binder_rmsd = _rmsd(threaded, pose, threaded.chain_begin(2), threaded.total_residue())
        combined_rmsd = _rmsd(threaded, pose, 1, threaded.total_residue())
        score = float(objects.get_score_function("full")(pose))
        if not all(math.isfinite(value) for value in (score, peptide_rmsd, binder_rmsd, combined_rmsd)):
            raise ValueError("Rosetta returned non-finite geometry or score values.")
        identifier = f"design_{replicate:04d}"
        path = output / f"{identifier}.pdb"
        pose.dump_pdb(str(path))
        retained = max(peptide_rmsd, binder_rmsd, combined_rmsd) <= request["max_backbone_rmsd"]
        records.append({"pdb": path.name, "binder_sequence": pose.chain_sequence(2),
                        "target_sequence": fragment, "accepted": retained,
                        "metrics": {"rosetta_total_score": score, "template_peptide_backbone_rmsd": peptide_rmsd,
                                    "template_binder_backbone_rmsd": binder_rmsd,
                                    "template_complex_backbone_rmsd": combined_rmsd,
                                    "passes_template_geometry": retained}})
        (output / "designs.json").write_text(json.dumps({"schema_version": 1, "records": records,
            "protocol": "Deposited interface FastDesign/FastRelax; external analysis omitted.",
            "source_doi": "10.5281/zenodo.10391229"}, indent=2) + "\n")
    return {"records": records}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    run(json.loads(args.request.read_text()), args.output, validate_only=args.validate_only)
