"""Input validation and target preparation for the upstream RFpeptides protocol."""
from __future__ import annotations

import re
from pathlib import Path

from mn_molecular_design_studio.design_validation import protein_chains, read_structure

UPSTREAM_COMMIT = "86507b6538f51fce57b5a72477165f03999ed7ae"
CODE_VERSION = f"RFpeptides · RFdiffusion {UPSTREAM_COMMIT[:8]}"
IMAGE = "mn-rfpeptides:86507b6-cu128"
WEIGHTS_RELATIVE_PATH = Path("rfdiffusion_models/Complex_base_ckpt.pt")
PDB_CHAIN_IDS = "BCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
STANDARD_RESIDUES = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE",
})


def parse_binder_length(value: object) -> str:
    text = str(value or "").strip()
    match = re.fullmatch(r"(\d+)(?:\s*-\s*(\d+))?", text)
    if not match:
        raise ValueError("Binder length must be a residue count or a range such as 12-18.")
    start, end = int(match.group(1)), int(match.group(2) or match.group(1))
    if not 3 <= start <= end <= 200:
        raise ValueError("Binder length must be between 3 and 200 residues, with the minimum no greater than the maximum.")
    return str(start) if start == end else f"{start}-{end}"


def build_contig(binder_length: object, target_chains: list[str], chain_lengths: dict[str, int]) -> str:
    """Create the RFdiffusion contig syntax for a generated peptide and receptor chains."""
    length = parse_binder_length(binder_length)
    if not target_chains:
        raise ValueError("Select at least one protein target chain.")
    target_segments = []
    for chain in target_chains:
        if len(chain) != 1 or not chain.isalpha():
            raise ValueError(f"RFdiffusion requires a single-letter target chain ID, got {chain!r}.")
        residue_count = chain_lengths.get(chain)
        if not isinstance(residue_count, int) or residue_count < 1:
            raise ValueError(f"Target chain {chain} has no protein residues.")
        target_segments.append(f"{chain}1-{residue_count}/0")
    return f"{length} {' '.join(target_segments)}"


def normalize_hotspots(value: object, target_chains: list[str], chain_residues: dict | None = None) -> list[str]:
    if isinstance(value, (list, tuple, set)):
        text = ",".join(map(str, value))
    else:
        text = str(value or "")
    chosen_chains = {str(chain).strip() for chain in target_chains}
    normalized: list[str] = []
    for token in re.split(r"[,;\s]+", text.strip()):
        if not token:
            continue
        match = re.fullmatch(r"(.+?)(-?\d+)([A-Za-z]?)", token)
        if not match:
            raise ValueError(f"Hotspot {token!r} must include a chain ID and residue number, for example A51.")
        chain, number, insertion = match.groups()
        if chain not in chosen_chains:
            raise ValueError(f"Hotspot {token!r} is not in the selected target chains.")
        label = f"{chain}{int(number)}{insertion.upper()}"
        if chain_residues is not None:
            valid = {
                (int(residue.id[1]), str(residue.id[2]).strip().upper())
                for residue in chain_residues.get(chain, [])
            }
            if (int(number), insertion.upper()) not in valid:
                raise ValueError(f"Hotspot {token!r} is not a protein residue in the selected structure.")
        if label not in normalized:
            normalized.append(label)
    return normalized


def validate_request(parameters: dict, *, require_files: bool = True) -> dict:
    """Validate RFpeptides inputs without creating jobs or initializing a model."""
    target_pdb = Path(parameters["target_pdb"]).expanduser().resolve() if parameters.get("target_pdb") else None
    chains = parameters.get("target_chains") or []
    if isinstance(chains, str):
        chains = [item.strip() for item in chains.split(",") if item.strip()]
    chains = list(dict.fromkeys(str(chain).strip() for chain in chains if str(chain).strip()))
    binder_length = parse_binder_length(parameters.get("binder_length", "12-18"))
    number_of_designs = parameters.get("number_of_designs", 1)
    diffusion_steps = parameters.get("diffusion_steps", 50)
    cpu_cores = parameters.get("cpu_cores", 4)
    for key, value, minimum, maximum in (
        ("Number of designs", number_of_designs, 1, 100000),
        ("Diffusion steps", diffusion_steps, 1, 200),
        ("CPU cores", cpu_cores, 1, 256),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
            raise ValueError(f"{key} must be an integer between {minimum} and {maximum}.")
    if require_files and target_pdb is None:
        raise ValueError("Load a target PDB or mmCIF structure first.")
    if target_pdb is not None and require_files:
        if not target_pdb.is_file() or target_pdb.suffix.lower() not in {".pdb", ".ent", ".cif", ".mmcif"}:
            raise ValueError("Load a valid PDB or mmCIF target structure.")
        if target_pdb.stat().st_size > 50_000_000:
            raise ValueError("Target structure files must be smaller than 50 MB.")
        structure = read_structure(target_pdb)
        chain_residues = protein_chains(structure)
        if not chains:
            raise ValueError("Select at least one protein target chain.")
        absent = [chain for chain in chains if chain not in chain_residues]
        if absent:
            raise ValueError(f"Selected target chain(s) are missing or contain no protein residues: {', '.join(absent)}.")
        hotspots = normalize_hotspots(parameters.get("hotspots", ""), chains, chain_residues)
    else:
        hotspots = normalize_hotspots(parameters.get("hotspots", ""), chains)
        if chains and len(chains) != len(set(chains)):
            raise ValueError("Select each target chain only once.")
    return {
        **parameters,
        "target_pdb": target_pdb,
        "target_chains": chains,
        "hotspots": hotspots,
        "binder_length": binder_length,
        "number_of_designs": number_of_designs,
        "diffusion_steps": diffusion_steps,
        "cpu_cores": cpu_cores,
    }


def prepare_target(source: Path, selected_chains: list[str], hotspots: list[str], destination: Path) -> dict:
    """Write selected protein chains as PDB, remapping them away from peptide chain A."""
    from Bio.PDB import PDBIO, Select

    structure = read_structure(Path(source))
    model = next(structure.get_models(), None)
    if model is None:
        raise ValueError("The selected target structure has no coordinate model.")
    chains = protein_chains(structure)
    missing = [chain for chain in selected_chains if chain not in chains]
    if missing:
        raise ValueError(f"Selected target chain(s) are missing: {', '.join(missing)}.")
    if not selected_chains:
        raise ValueError("Select at least one protein target chain.")
    if len(selected_chains) > len(PDB_CHAIN_IDS):
        raise ValueError("RFpeptides can represent at most 51 target chains in PDB format.")

    selected_set = set(selected_chains)
    for chain in list(model):
        if chain.id not in selected_set:
            model.detach_child(chain.id)
    chain_map = {source_chain: PDB_CHAIN_IDS[index] for index, source_chain in enumerate(selected_chains)}
    residue_map: dict[tuple[str, int, str], int] = {}
    for source_chain in selected_chains:
        chain = model[source_chain]
        kept = 0
        for residue in list(chain):
            if residue.id[0] != " " or residue.resname not in STANDARD_RESIDUES or "CA" not in residue:
                chain.detach_child(residue.id)
                continue
            kept += 1
            residue_map[(source_chain, int(residue.id[1]), str(residue.id[2]).strip().upper())] = kept
            residue.id = (" ", kept, " ")
        if kept == 0:
            raise ValueError(f"Selected target chain {source_chain} has no standard protein residues.")
        chain.id = chain_map[source_chain]
    mapped_hotspots = []
    for hotspot in hotspots:
        match = re.fullmatch(r"(.+?)(-?\d+)([A-Za-z]?)", hotspot)
        if not match:
            raise ValueError(f"Invalid hotspot residue {hotspot!r}.")
        source_chain, number, insertion = match.groups()
        key = (source_chain, int(number), insertion.upper())
        if key not in residue_map:
            raise ValueError(f"Hotspot {hotspot} is missing from the selected target protein chains.")
        mapped_hotspots.append(f"{chain_map[source_chain]}{residue_map[key]}")

    class _SelectedProtein(Select):
        def accept_chain(self, chain):
            return chain.id in set(chain_map.values())

        def accept_residue(self, residue):
            return residue.id[0] == " " and residue.resname in STANDARD_RESIDUES and "CA" in residue

        def accept_atom(self, atom):
            return atom.element != "H" and (not atom.is_disordered() or atom.get_altloc() in {"A", " "})

    destination.parent.mkdir(parents=True, exist_ok=True)
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(destination), _SelectedProtein())
    if not destination.is_file() or "ATOM" not in destination.read_text(errors="ignore"):
        raise ValueError("No standard protein atoms could be written to the prepared target.")
    return {
        "chain_map": chain_map,
        "target_chains": list(chain_map.values()),
        "hotspots": mapped_hotspots,
        "chain_lengths": {chain_map[chain]: len(chains[chain]) for chain in selected_chains},
    }
