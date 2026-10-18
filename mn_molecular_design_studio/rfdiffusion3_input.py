"""Validated task inputs for the RFdiffusion3 version installed in mn-foundry."""
from __future__ import annotations

import copy
import json
import math
import re
from pathlib import Path

from Bio.SeqUtils import seq1

from mn_molecular_design_studio.design_validation import read_structure

TASKS = {
    "protein_binder": "Protein binder design",
    "enzyme": "Enzyme scaffolding",
    "nucleic_acid_binder": "Nucleic acid binder design",
    "ligand_binder": "Small molecule binder design",
    "json": "Custom JSON",
}
IMAGE = "mn-foundry:cu128"
CODE_VERSION = "RFdiffusion3 a36d29c · Foundry 0.1.9"
ASSETS = Path(__file__).parent / "assets/rfdiffusion3"
DNA = {"DA", "DC", "DG", "DT", "DU", "DI"}
RNA = {"A", "C", "G", "U", "I"}
WATER = {"HOH", "WAT", "DOD"}
SPEC_FIELDS = {
    "input", "contig", "unindex", "length", "ligand", "cif_parser_args", "extra", "dialect",
    "select_fixed_atoms", "select_unfixed_sequence", "select_buried", "select_partially_buried",
    "select_exposed", "select_hbond_acceptor", "select_hbond_donor", "select_hotspots",
    "redesign_motif_sidechains", "symmetry", "ori_token", "infer_ori_strategy", "plddt_enhanced",
    "is_non_loopy", "partial_t",
}
CONDITIONS = {"Fixed": "select_fixed_atoms", "Buried": "select_buried", "Exposed": "select_exposed",
              "H-bond donor": "select_hbond_donor", "H-bond acceptor": "select_hbond_acceptor"}


def inventory(path: Path) -> dict:
    try:
        model = next(read_structure(Path(path)).get_models())
    except Exception as exc:
        raise ValueError(f"Could not read structure {Path(path).name}: {exc}") from exc
    chains, ligands, residues = {}, {}, {}
    for chain in model:
        for residue in chain:
            name = residue.resname.strip()
            if name in WATER:
                continue
            if residue.id[2].strip():
                # The installed selector language cannot express insertion codes.
                label = f"{chain.id}{residue.id[1]}{residue.id[2].strip()}"
            else:
                label = f"{chain.id}{residue.id[1]}"
            kind = "dna" if name in DNA else "rna" if name in RNA else "protein" if "CA" in residue and seq1(name, undef_code="X") != "X" else "ligand"
            row = {"label": label, "number": residue.id[1], "name": name, "kind": kind,
                   "insertion_code": residue.id[2].strip(), "atoms": [a.name for a in residue if a.element != "H"]}
            residues[label] = row
            if kind == "ligand":
                ligands.setdefault(name, []).append(row)
            else:
                chains.setdefault(chain.id, {"kind": kind, "residues": []})["residues"].append(row)
    if not residues:
        raise ValueError("The structure contains no usable non-water atoms.")
    return {"chains": chains, "ligands": ligands, "residues": residues}


def chain_contig(chain: str, rows: list[dict]) -> str:
    if not re.fullmatch(r"[A-Za-z]+", chain) or any(r["insertion_code"] for r in rows):
        raise ValueError("This image's contig selectors require letter chain IDs and residue numbers without insertion codes. Prepare/remap the structure first.")
    numbers = [r["number"] for r in rows]
    groups = []
    for number in numbers:
        if groups and number == groups[-1][-1] + 1:
            groups[-1].append(number)
        else:
            groups.append([number])
    return ",".join(f"{chain}{g[0]}-{g[-1]}" if len(g) > 1 else f"{chain}{g[0]}" for g in groups)


def _selected_rows(selector: str, info: dict) -> list[dict]:
    if selector in info["ligands"]:
        return info["ligands"][selector]
    rows = []
    for token in str(selector).split(","):
        match = re.fullmatch(r"([A-Za-z]+)(-?\d+)(?:-(-?\d+))?", token.strip())
        if not match:
            raise ValueError(f"Invalid residue selector {token!r}; use A108 or A108-110.")
        chain, start, end = match.groups(); start = int(start); end = int(end or start)
        selected = [r for r in info["chains"].get(chain, {}).get("residues", []) if start <= r["number"] <= end]
        if not selected or any(r["insertion_code"] for r in selected):
            raise ValueError(f"Selector {token!r} does not match addressable residues in the input.")
        rows.extend(selected)
    return rows


def validate_selection(value, info: dict, *, field: str):
    if value is None or isinstance(value, bool):
        return
    if isinstance(value, str):
        for token in value.split(","):
            if token.strip() and token.strip() != "/0" and not re.fullmatch(r"\d+(?:-\d+)?", token.strip()):
                _selected_rows(token.strip(), info)
        return
    if not isinstance(value, dict) or not value:
        raise ValueError(f"{field} must be a selection string, boolean or non-empty atom dictionary.")
    for selector, atoms in value.items():
        rows = _selected_rows(selector, info)
        if not isinstance(atoms, str):
            raise ValueError(f"{field}: atom names must be comma-separated text.")
        if atoms.strip() in {"ALL", "BKBN", "TIP", ""}:
            continue
        for row in rows:
            unknown = set(a.strip() for a in atoms.split(",")) - set(row["atoms"])
            if unknown:
                raise ValueError(f"{field}: {row['label']} has no atom(s) {', '.join(sorted(unknown))}.")


def _length(value, label="length") -> str:
    text = str(value)
    if not re.fullmatch(r"[1-9]\d*(?:-[1-9]\d*)?", text):
        raise ValueError(f"{label} must be a residue count or min-max range.")
    values = list(map(int, text.split("-")))
    if values[0] > values[-1] or values[-1] > 3000:
        raise ValueError(f"{label} must be ordered and no greater than 3,000 residues.")
    return text


def prepare_request(parameters: dict) -> dict:
    """No jobs, copies, server requests, Docker calls or model initialization."""
    p = parameters
    task = p.get("task", "protein_binder")
    if task not in TASKS:
        raise ValueError("Select a supported RFdiffusion3 task.")
    for key, default, minimum, maximum in (("num_designs", 8, 1, 10000), ("diffusion_batch_size", 1, 1, 32),
                                          ("timesteps", 200, 1, 1000), ("cpu_cores", 4, 1, 256), ("seed", 0, 0, 2147483647)):
        value = p.get(key, default)
        if type(value) is not int or not minimum <= value <= maximum:
            raise ValueError(f"{key} must be an integer between {minimum} and {maximum}.")
    for key, default, low, high in (("step_scale", 1.5, 1., 5.), ("gamma_0", .6, 0., 1.), ("cfg_scale", 1.5, 0., 10.)):
        value = p.get(key, default)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not low <= value <= high:
            raise ValueError(f"{key} must be between {low} and {high}.")
    low, high = p.get("length_min", 100), p.get("length_max", 150)
    length = _length(f"{low}-{high}")
    target = Path(p["target_structure"]) if p.get("target_structure") else None
    files = [Path(f) for f in p.get("input_files") or []]
    if target and target not in files:
        files.insert(0, target)
    if len(files) > 128:
        raise ValueError("Upload at most 128 structure files.")
    by_name = {}
    for file in files:
        if not file.is_file() or file.suffix.lower() not in {".pdb", ".ent", ".cif", ".mmcif"}:
            raise ValueError(f"Structure file {file.name} is missing or is not PDB/mmCIF.")
        if file.stat().st_size > 50_000_000:
            raise ValueError("Structure files must be smaller than 50 MB.")
        if file.name in by_name and by_name[file.name] != file:
            raise ValueError(f"Two uploaded structures share the filename {file.name}.")
        by_name[file.name] = file
    contexts = {}
    if task == "json":
        specs = copy.deepcopy(p.get("json_config"))
        if not isinstance(specs, dict) or not specs or len(specs) > 128:
            raise ValueError("Provide a JSON object containing between 1 and 128 named design specifications.")
    else:
        if target is None:
            raise ValueError("Load a target/motif structure first.")
        info = inventory(target)
        spec = {"dialect": 2, "input": target.name}
        if p.get("is_non_loopy", True):
            spec["is_non_loopy"] = True
        ligands = p.get("ligands") or []
        if any(name not in info["ligands"] for name in ligands):
            raise ValueError("A selected ligand is absent from the input structure.")
        if ligands:
            spec["ligand"] = ",".join(ligands)
        chains = p.get("target_chains") or []
        if task in {"protein_binder", "nucleic_acid_binder"}:
            kinds = {"protein"} if task == "protein_binder" else {"dna", "rna"}
            if not chains:
                raise ValueError("Select at least one target chain.")
            if any(c not in info["chains"] or info["chains"][c]["kind"] not in kinds for c in chains):
                raise ValueError("Selected target chains do not match this task's molecule type.")
            target_contig = ",/0,".join(chain_contig(c, info["chains"][c]["residues"]) for c in chains)
            spec["contig"] = f"{length},/0,{target_contig}" if task == "protein_binder" else f"{target_contig},/0,{length}"
            if p.get("hotspots", "").strip():
                spec["select_hotspots"] = {s.strip(): "CA" if task == "protein_binder" else "ALL" for s in p["hotspots"].split(",") if s.strip()}
                spec["infer_ori_strategy"] = "hotspots"
        elif task == "enzyme":
            if p.get("scaffolding_mode", "atom") == "residue":
                if not p.get("motif_contig", "").strip():
                    raise ValueError("Specify a residue-level scaffold contig.")
                spec["contig"] = p["motif_contig"].strip()
                validate_selection(spec["contig"], info, field="motif_contig")
            elif p.get("scaffolding_mode", "atom") == "atom":
                fixed = p.get("fixed_atoms")
                if not isinstance(fixed, dict) or not fixed:
                    raise ValueError("Choose motif residues and at least one fixed atom per residue.")
                if any(not isinstance(atoms, str) or not atoms.strip() for atoms in fixed.values()):
                    raise ValueError("Each catalytic motif residue needs at least one fixed atom.")
                spec["select_fixed_atoms"] = dict(fixed)
                unindex = p.get("unindex", "all").strip()
                if unindex == "all":
                    unindex = ",".join(k for k in fixed if k not in info["ligands"])
                if not unindex:
                    raise ValueError("Choose at least one unindexed motif residue.")
                spec.update(unindex=unindex, length=length)
            else:
                raise ValueError("Scaffolding mode must be atom or residue.")
        elif task == "ligand_binder":
            if not ligands:
                raise ValueError("Select at least one small-molecule ligand.")
            spec["length"] = length
        for row in p.get("atom_conditions") or []:
            if not isinstance(row, dict) or row.get("condition") not in CONDITIONS or row.get("ligand") not in ligands:
                raise ValueError("Each atom condition needs a condition type and a selected ligand.")
            field = CONDITIONS[row["condition"]]
            atoms = str(row.get("atoms") or "").strip()
            if field != "select_fixed_atoms" and not atoms:
                raise ValueError("Burial, exposure and H-bond conditions need atoms or ALL.")
            selection = spec.setdefault(field, {})
            if row["ligand"] in selection:
                raise ValueError("Combine atoms for the same ligand and condition in one row.")
            selection[row["ligand"]] = atoms
        if p.get("ori_token") is not None:
            spec["ori_token"] = p["ori_token"]
        specs = {"design": spec}
    inventories = {name: inventory(path) for name, path in by_name.items()}
    normalized = {}
    for index, (name, spec) in enumerate(specs.items(), 1):
        if not isinstance(spec, dict) or set(spec) - SPEC_FIELDS:
            raise ValueError(f"{name}: unsupported specification fields for installed Foundry 0.1.9: {sorted(set(spec) - SPEC_FIELDS) if isinstance(spec, dict) else 'expected object'}")
        if spec.get("dialect", 2) != 2:
            raise ValueError("Use RFdiffusion3 dialect 2.")
        if "length" in spec:
            spec["length"] = _length(spec["length"])
        if not any(spec.get(k) for k in ("contig", "length", "unindex")):
            raise ValueError(f"{name}: specify a contig or design length.")
        filename = None
        if spec.get("input"):
            filename = Path(str(spec["input"])).name
            if filename not in by_name:
                raise ValueError(f"{name}: upload the referenced input structure {filename}.")
            spec["input"] = filename
        if filename:
            info = inventories[filename]
            if any(code not in info["ligands"] for code in str(spec.get("ligand") or "").split(",") if code):
                raise ValueError(f"{name}: a selected ligand is absent from the input structure.")
            for field, value in spec.items():
                if field.startswith("select_"):
                    validate_selection(value, info, field=field)
            for field in ("contig", "unindex"):
                if spec.get(field):
                    validate_selection(spec[field], info, field=field)
        elif any(spec.get(k) for k in ("unindex", "ligand", "select_fixed_atoms")) or any(spec.get(k) for k in spec if k.startswith("select_")) or re.search(r"[A-Za-z]", str(spec.get("contig", ""))):
            raise ValueError(f"{name}: motif/ligand selections require an uploaded structure.")
        origin = spec.get("ori_token")
        if origin is not None and (not isinstance(origin, list) or len(origin) != 3 or any(not isinstance(v, (int, float)) or isinstance(v, bool) or not math.isfinite(v) for v in origin)):
            raise ValueError("Origin must be three finite coordinates [x,y,z].")
        identifier = f"spec_{index:03d}"
        normalized[identifier] = spec
        contexts[identifier] = {"name": str(name), "task": task, "source_path": str(by_name[filename]) if filename else None,
                                "target_chains": list(p.get("target_chains") or []) if task in {"protein_binder", "nucleic_acid_binder"} else []}
    return {"specifications": normalized, "files": by_name, "contexts": contexts}
