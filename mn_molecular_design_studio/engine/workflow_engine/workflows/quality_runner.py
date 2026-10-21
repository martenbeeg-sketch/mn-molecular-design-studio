"""Standalone container entry point for the studio.qc.v1 batch contract.

No Studio imports: this file is mounted read-only into pinned tool images.
Large residue/contact results are written as artifacts, never SQL scalar values.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def sequence_metrics(sequence, ph=7.4):
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    sequence = sequence.upper()
    if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
        raise ValueError("Sequence composition requires a nonempty sequence of standard amino acids; unknown residues are not discarded.")
    analysis = ProteinAnalysis(sequence)
    frequencies = Counter(sequence)
    return {"qc_sequence_length": len(sequence), "qc_sequence_pi": analysis.isoelectric_point(),
            "qc_sequence_charge": analysis.charge_at_pH(float(ph)), "qc_sequence_ph": float(ph),
            "qc_sequence_gravy": analysis.gravy(), "qc_sequence_aromaticity": analysis.aromaticity(),
            "qc_sequence_instability": analysis.instability_index(),
            "qc_sequence_entropy": -sum((n / len(sequence)) * math.log2(n / len(sequence)) for n in frequencies.values()),
            **{"qc_sequence_fraction_" + aa: frequencies[aa] / len(sequence) for aa in "ACDEFGHIKLMNPQRSTVWY"}}


def atoms_from(path):
    from biotite.structure.io import load_structure
    array = load_structure(path, model=1)
    if len(array) == 0:
        raise ValueError("The structure contains no atoms.")
    return array


def interface_metrics(item, settings, destination):
    import numpy as np
    import biotite.structure as struc
    atoms = atoms_from(item["structure"])
    binder, target = item["binder_chains"], item["target_chains"]
    available = set(atoms.chain_id)
    if not binder or not target or set(binder) & set(target) or (set(binder) | set(target)) - available:
        raise ValueError("Interface analysis needs distinct binder and target chains present in the complex.")
    protein = struc.filter_amino_acids(atoms)
    heavy = ~np.isin(atoms.element, ["H", "D"])
    b = np.where(protein & heavy & np.isin(atoms.chain_id, binder))[0]
    t = np.where(protein & heavy & np.isin(atoms.chain_id, target))[0]
    if not len(b) or not len(t):
        raise ValueError("No protein heavy atoms found in the selected chain groups.")
    cutoff = float(settings.get("contact_distance", 5.0))
    cell = struc.CellList(atoms.coord[t], cell_size=max(cutoff, 4.0))
    near = cell.get_atoms(atoms.coord[b], radius=max(cutoff, 4.0))
    near_pairs = {(int(b[i]), int(t[j])) for i, neighbors in enumerate(near) for j in neighbors if j >= 0}
    pairs = {(i, j) for i, j in near_pairs if np.linalg.norm(atoms.coord[i] - atoms.coord[j]) <= cutoff}
    def residue(index):
        return str(atoms.chain_id[index]) + ":" + str(atoms.res_id[index]) + str(atoms.ins_code[index])
    residue_pairs = {(residue(i), residue(j)) for i, j in pairs}
    positive = ((atoms.res_name == "LYS") & (atoms.atom_name == "NZ")) | ((atoms.res_name == "ARG") & np.isin(atoms.atom_name, ["NE", "NH1", "NH2"]))
    negative = ((atoms.res_name == "ASP") & np.isin(atoms.atom_name, ["OD1", "OD2"])) | ((atoms.res_name == "GLU") & np.isin(atoms.atom_name, ["OE1", "OE2"]))
    salt_pairs = {(residue(i), residue(j)) for i, j in near_pairs if (positive[i] and negative[j] or negative[i] and positive[j]) and np.linalg.norm(atoms.coord[i] - atoms.coord[j]) <= 4.0}
    rows = [{"binder_residue": i, "target_residue": j} for i, j in sorted(residue_pairs)]
    write_json(destination / "contacts.json", {"schema": "studio.contacts.v1", "distance_angstrom": cutoff,
               "binder_chains": binder, "target_chains": target, "residue_pairs": rows,
               "target_residues": sorted({j for _, j in residue_pairs}), "salt_bridge_pairs": sorted(salt_pairs)})
    metrics = {"qc_interface_contact_pairs": len(pairs), "qc_interface_residue_pairs": len(residue_pairs),
               "qc_interface_binder_residues": len({i for i, _ in residue_pairs}),
               "qc_interface_target_residues": len({j for _, j in residue_pairs}), "qc_interface_salt_bridges": len(salt_pairs),
               "qc_interface_contact_distance": cutoff}
    if np.any(atoms.element == "H"):
        triplets = struc.hbond(atoms)
        groups_b, groups_t = set(binder), set(target)
        metrics["qc_interface_hydrogen_bonds"] = sum(1 for donor, _, acceptor in triplets
            if atoms.chain_id[donor] in groups_b and atoms.chain_id[acceptor] in groups_t or
               atoms.chain_id[donor] in groups_t and atoms.chain_id[acceptor] in groups_b)
    else:
        metrics["qc_interface_hbond_status"] = "unavailable: explicit hydrogens required"
    return metrics


def select_pdb(item, destination):
    """Copy the declared protein chain group to a PDB for external evaluators."""
    import numpy as np
    import biotite.structure as struc
    from biotite.structure.io.pdb import PDBFile
    atoms = atoms_from(item["structure"])
    chains = item["binder_chains"]
    if not chains or set(chains) - set(atoms.chain_id):
        raise ValueError("Choose protein chains to analyze; the declared chains must exist in the structure.")
    atoms = atoms[np.isin(atoms.chain_id, chains) & struc.filter_amino_acids(atoms)]
    if not len(atoms):
        raise ValueError("Selected chains contain no amino acids.")
    output = destination / "selected.pdb"
    pdb = PDBFile()
    pdb.set_structure(atoms)
    pdb.write(output)
    return output


def proteinsol(item, destination):
    sequence = item["sequence"]
    sequence_metrics(sequence)  # Validate without guessing unknown residue identities.
    if len(sequence) < 21:
        raise ValueError("Protein-Sol requires at least 21 residues.")
    from importlib.util import spec_from_file_location, module_from_spec
    spec = spec_from_file_location("protein_sol_api", Path(__file__).with_name("qc_ovo") / "protein_sol_api.py")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    # The upstream script runs its Perl bundle in a private temporary directory.
    module.PROTEIN_SOL_PATH = "/opt/protein-sol-sequence-prediction-software"
    result = module.protein_sol([sequence])[0]
    write_json(destination / "proteinsol.json", result)
    return {"proteinsol_scaled_solubility": result["scaled-sol"],
            "proteinsol_percent_solubility": result["percent-sol"], "proteinsol_population_solubility": result["population-sol"]}


def pep_patch(item, settings, destination):
    import sys
    sys.path.insert(0, str(Path(__file__).with_name("qc_ovo")))
    from lambda_function import process_pdb, get_peppatch_set_definitions
    pdb = select_pdb(item, destination)
    configuration = json.loads((Path(__file__).with_name("qc_ovo") / "feature_set.json").read_text())
    configuration["methods"] = [m for m in configuration["methods"] if m["method"] == "hydrophobic" or m["name"] == "electrostatic_7.4"]
    for method in configuration["methods"]:
        if method["method"] == "electrostatic":
            method["name"] = "electrostatic"
            method["config"]["pH"] = float(settings.get("ph", 7.4))
    config_path = destination / "feature_set.json"
    write_json(config_path, configuration)
    definitions = get_peppatch_set_definitions(str(config_path))
    # Config paths are interpreted relative to the external image's resource directory.
    for definition in definitions["methods"]:
        if definition["method"] == "hydrophobic":
            definition["config"]["hydrophobicity_scale"] = "/opt/resources/hydrophobicity_scales/scale_Ja.csv"
    residue_data, descriptors = process_pdb(str(pdb), definitions)
    if hasattr(residue_data, "to_csv"):
        residue_data.to_csv(destination / "residue_patches.csv")
    elif isinstance(residue_data, dict):
        import hashlib
        for name, data in residue_data.items():
            if hasattr(data, "to_csv"):
                data.to_csv(destination / ("residue_patches_" + hashlib.sha1(str(name).encode()).hexdigest()[:10] + ".csv"))
    scalar = {"peppatch_" + str(k): float(v) for k, v in descriptors.items() if isinstance(v, (int, float)) and math.isfinite(v)}
    total = scalar.get("peppatch_total_area")
    if total:
        scalar.update({k + "_fraction": v / total for k, v in list(scalar.items()) if "patch_area" in k})
    scalar["peppatch_ph"] = float(settings.get("ph", 7.4))
    write_json(destination / "patch_metrics.json", scalar)
    return scalar


def run_item(item, tool, settings, output):
    destination = Path(output) / item["id"]
    destination.mkdir(parents=True, exist_ok=True)
    previous_directory = Path.cwd()
    try:
        # Some native surface tools use fixed intermediate filenames.
        # Each process must run inside its item's own working directory.
        os.chdir(destination)
        if item.get("input_error"):
            raise ValueError(item["input_error"])
        if tool == "sequence_qc":
            metrics = sequence_metrics(item["sequence"], settings.get("ph", 7.4))
        elif tool in {"biotite_interface", "interface_clustering"}:
            metrics = interface_metrics(item, settings, destination)
        elif tool == "proteinsol":
            metrics = proteinsol(item, destination)
        elif tool == "pep_patch":
            metrics = pep_patch(item, settings, destination)
        elif tool in {"plip", "pandamap", "posebusters"}:
            metrics = ligand_analysis(item, tool, destination)
        else:
            raise ValueError("Unsupported QC evaluator: " + tool)
        result = {"id": item["id"], "status": "completed", "metrics": metrics}
    except Exception as exc:
        result = {"id": item["id"], "status": "failed", "metrics": {}, "error": str(exc)}
    finally:
        os.chdir(previous_directory)
    write_json(destination / "result.json", result)
    return result


def ligand_analysis(item, tool, destination):
    """Reuse the installed, numbered MN Ligand adapters with one explicit SDF pose."""
    import sys
    source = destination / "inputs/poses.csv"
    source.parent.mkdir(parents=True, exist_ok=True)
    fields = {"pose_id": item["id"], "compound_id": item["candidate_id"], "smiles": item["smiles"],
              "mol_pred": item["ligand_file"], "mol_cond": item["receptor_file"], "reference_target": item["receptor_file"],
              "source_kind": "studio_selected_pose", "selection_criterion": "best source pose or selected prediction"}
    with source.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(fields))
        writer.writeheader()
        writer.writerow(fields)
    summary = destination / "summary.csv"
    report = destination / "report.json"
    script = {"plip": "/opt/mn-ligand/mn_ligand_runner.py", "pandamap": "/opt/pandamap/mn_ligand_runner.py", "posebusters": "/opt/posebusters/mn_ligand_runner.py"}[tool]
    command = [sys.executable, script, "--input", str(source), "--summary", str(summary), "--report", str(report), "--max-workers", "1"]
    if tool == "posebusters":
        command += ["--output", str(destination / "full_checks.csv")]
    else:
        command += ["--native-output", str(destination / "native"), "--interactions", str(destination / "interactions.csv")]
    with (destination / "stdout.log").open("w") as stdout, (destination / "stderr.log").open("w") as stderr:
        result = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
    if not summary.is_file():
        raise ValueError(f"{tool} exited {result.returncode} without a summary; see the per-item stderr artifact.")
    row = next(csv.DictReader(summary.open()), {})
    if row.get("success", "True").lower() == "false" or row.get("input_preparation_error"):
        raise ValueError(row.get("error") or row.get("input_preparation_error") or "The native analysis adapter rejected this pose.")
    metrics = {}
    for key, value in row.items():
        if value.lower() in {"true", "false"}:
            metrics[tool + "_" + key] = value.lower() == "true"
        else:
            try:
                number = float(value)
                if math.isfinite(number):
                    metrics[tool + "_" + key] = number
            except ValueError:
                pass
    return metrics


def cluster(items, tool, settings, output):
    """Deterministic single-linkage clusters; representatives follow source order."""
    threshold = float(settings.get("similarity_threshold", .8))
    parent = list(range(len(items)))
    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    if tool == "foldseek":
        inputs = Path(output) / "foldseek_inputs"
        inputs.mkdir(exist_ok=True)
        names = {}
        for i, item in enumerate(items):
            if len(item["binder_chains"]) != 1:
                raise ValueError("This Foldseek protocol clusters one protein chain per candidate; use interface clustering for multichain binding epitopes.")
            folder = inputs / item["id"]
            folder.mkdir(exist_ok=True)
            path = select_pdb(item, folder)
            destination = inputs / (item["id"] + ".pdb")
            shutil.move(path, destination)
            names[item["id"]] = i
        prefix = Path(output) / "foldseek"
        subprocess.run(["foldseek", "easy-cluster", str(inputs), str(prefix), str(Path(output) / "tmp"),
                        "--alignment-type", "1", "--tmscore-threshold", str(threshold), "-c", str(settings.get("coverage", .8)),
                        "--cov-mode", "0", "--threads", str(settings.get("workers", 1))], check=True)
        for representative, member in csv.reader(Path(str(prefix) + "_cluster.tsv").open(), delimiter="\t"):
            parent[names[Path(member).stem]] = names[Path(representative).stem]
    else:
        if len(items) > 2000:
            raise ValueError("Pairwise sequence/interface clustering supports at most 2,000 items per batch; use Foldseek for larger structure cohorts.")
        if tool == "interface_clustering":
            targets = {item.get("target_key") for item in items}
            if len(targets) != 1 or not all(targets):
                raise ValueError("Interface clustering requires one identical target structure and chain numbering across the cohort.")
            sets = []
            for item in items:
                result = run_item(item, tool, settings, output)
                if result["status"] != "completed":
                    raise ValueError(result["error"])
                sets.append(set(json.loads((Path(output) / item["id"] / "contacts.json").read_text())["target_residues"]))
            similarity = lambda i, j: len(sets[i] & sets[j]) / len(sets[i] | sets[j]) if sets[i] | sets[j] else 0.
        else:
            from Bio.Align import PairwiseAligner
            for item in items:
                sequence_metrics(item["sequence"])
            aligner = PairwiseAligner(mode="global", match_score=1, mismatch_score=0, open_gap_score=-1, extend_gap_score=-.1)
            def similarity(i, j):
                a, b = items[i]["sequence"], items[j]["sequence"]
                alignment = aligner.align(a, b)[0]
                return alignment.counts().identities / max(len(a), len(b))
        for i in range(len(items)):
            for j in range(i):
                if similarity(i, j) >= threshold:
                    parent[root(i)] = root(j)
    groups = {}
    for i, item in enumerate(items):
        groups.setdefault(root(i), []).append(i)
    results = []
    for number, indices in enumerate(sorted(groups.values(), key=lambda values: min(values)), 1):
        representative = min(indices)
        for i in indices:
            result = {"id": items[i]["id"], "status": "completed", "metrics": {
                "diversity_cluster": number, "diversity_cluster_size": len(indices),
                "diversity_representative": i == representative, "diversity_representative_id": items[representative]["candidate_id"],
                "diversity_method": tool, "diversity_threshold": threshold}}
            write_json(Path(output) / items[i]["id"] / "result.json", result)
            results.append(result)
    return results


def esm_scores(items, tool, output, progress):
    """Load each checkpoint once per batch. Never generate replacement sequences."""
    import esm
    import torch
    torch.manual_seed(0)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    results = {}
    if tool == "esm1v_score":
        probabilities = {}
        for item in items:
            try:
                sequence_metrics(item["sequence"])
                if len(item["sequence"]) > 1022:
                    raise ValueError("ESM-1v supports at most 1,022 residues in this scoring protocol; windowed scores are a different protocol.")
                probabilities[item["id"]] = [0.] * len(item["sequence"])
            except ValueError as exc:
                results[item["id"]] = {"id": item["id"], "status": "failed", "metrics": {}, "error": str(exc)}
        for number in range(1, 6):
            model, alphabet = esm.pretrained.load_model_and_alphabet_local(f"/reference/esm_models/esm1v_t33_650M_UR90S_{number}.pt")
            model = model.eval().to(device)
            convert = alphabet.get_batch_converter()
            for item in items:
                if item["id"] in results:
                    continue
                try:
                    _, _, tokens = convert([(item["id"], item["sequence"])])
                    with torch.inference_mode():
                        values = model(tokens.to(device))["logits"][0, 1:len(item["sequence"]) + 1].softmax(-1)
                        indices = torch.tensor([alphabet.get_idx(aa) for aa in item["sequence"]], device=device)
                        native = values.gather(-1, indices[:, None]).squeeze(-1).cpu().tolist()
                    probabilities[item["id"]] = [a + b / 5 for a, b in zip(probabilities[item["id"]], native)]
                except Exception as exc:
                    results[item["id"]] = {"id": item["id"], "status": "failed", "metrics": {}, "error": str(exc)}
            del model
            if device == "cuda":
                torch.cuda.empty_cache()
            progress(list(results.values()), phase=f"Scored ESM-1v checkpoint {number}/5")
        for item in items:
            if item["id"] not in results:
                native = probabilities[item["id"]]
                write_json(Path(output) / item["id"] / "native_probabilities.json", {"sequence": item["sequence"], "probabilities": native, "protocol": "five-model unmasked wild-type marginals"})
                results[item["id"]] = {"id": item["id"], "status": "completed", "metrics": {"esm1v_native_probability": sum(native) / len(native)}}
    else:
        from esm.inverse_folding.multichain_util import load_complex_coords, score_sequence_in_complex
        model, alphabet = esm.pretrained.load_model_and_alphabet_local("/reference/esm_models/esm_if1_gvp4_t16_142M_UR50.pt")
        model = model.eval().to(device)
        for item in items:
            try:
                if item.get("input_error"):
                    raise ValueError(item["input_error"])
                chains = list(dict.fromkeys(item["binder_chains"] + item["target_chains"]))
                coords, sequences = load_complex_coords(item["structure"], chains)
                if not item["binder_chains"]:
                    raise ValueError("ESM-IF1 scoring requires selected protein chains.")
                declared = item["sequence"].replace(":", "")
                native = "".join(sequences[chain] for chain in item["binder_chains"])
                if declared != native:
                    raise ValueError("Declared binder sequence does not match the selected structure chains; score a refolded or sequence-labelled structure.")
                total, count = 0., 0
                per_chain = {}
                with torch.inference_mode():
                    for chain in item["binder_chains"]:
                        likelihood, _ = score_sequence_in_complex(model, alphabet, coords, chain, sequences[chain])
                        per_chain[chain] = float(likelihood)
                        total += float(likelihood) * len(sequences[chain])
                        count += len(sequences[chain])
                write_json(Path(output) / item["id"] / "chain_likelihoods.json", per_chain)
                results[item["id"]] = {"id": item["id"], "status": "completed", "metrics": {"esmif_native_log_likelihood": total / count}}
            except Exception as exc:
                results[item["id"]] = {"id": item["id"], "status": "failed", "metrics": {}, "error": str(exc)}
            progress(list(results.values()), phase="Scoring native sequences against complex backbones")
    return list(results.values())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    if manifest.get("schema") != "studio.qc.v1":
        raise ValueError("Unsupported QC batch contract.")
    tool, items, settings = manifest["tool"], manifest["items"], manifest["settings"]
    output = Path(manifest["output"])
    output.mkdir(parents=True, exist_ok=True)
    status = output / "batch_status.json"
    results = []
    def progress(partial=None, phase=""):
        write_json(status, {"schema": "studio.qc.v1", "total": len(items), "finished": len(results if partial is None else partial), "items": results if partial is None else partial, "phase": phase})
    progress()
    if tool in {"esm1v_score", "esmif_score"}:
        try:
            results = esm_scores(items, tool, output, progress)
        except Exception as exc:
            results = [{"id": item["id"], "status": "failed", "metrics": {}, "error": str(exc)} for item in items]
    elif tool in {"sequence_clustering", "interface_clustering", "foldseek"}:
        try:
            results = cluster(items, tool, settings, output)
        except Exception as exc:
            results = [{"id": item["id"], "status": "failed", "metrics": {}, "error": str(exc)} for item in items]
    else:
        with ProcessPoolExecutor(max_workers=max(1, min(int(settings.get("workers", 1)), len(items)))) as pool:
            futures = [pool.submit(run_item, item, tool, settings, output) for item in items]
            for future in as_completed(futures):
                results.append(future.result())
                progress()
    progress()
    return int(not any(result["status"] == "completed" for result in results))


if __name__ == "__main__":
    raise SystemExit(main())
