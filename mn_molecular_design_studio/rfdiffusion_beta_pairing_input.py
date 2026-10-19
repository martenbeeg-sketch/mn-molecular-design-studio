"""Input contract for RFdiffusion with explicit beta-pairing interface tensors."""
from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

from mn_molecular_design_studio.rfdiffusion3_input import inventory


WORKFLOW = "design.rfdiffusion_beta_pairing"
TOOL = "rfdiffusion_beta_pairing"
TITLE = "RFdiffusion β-pairing"
TENSOR_SCRIPT = Path(__file__).with_name("assets") / "rfdiffusion_beta_pairing" / "source" / "make_interface_tensor.py"
AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE",
}


def parse_residue_selection(value: str, chain: str) -> list[int]:
    """Expand a chain-qualified residue selection such as ``B99-103,B120``."""
    selected: list[int] = []
    for token in re.split(r"[,;\s]+", str(value or "").strip()):
        if not token:
            continue
        match = re.fullmatch(r"([A-Za-z])([1-9]\d*)(?:-([1-9]\d*))?", token)
        if not match or match[1] != chain:
            raise ValueError(f"Invalid β-strand selection {token!r}; use {chain}99-103 or {chain}99,{chain}100.")
        start = int(match[2])
        end = int(match[3] or start)
        if end < start or end - start > 1000:
            raise ValueError("β-strand residue ranges must increase and contain at most 1,001 residues.")
        selected.extend(range(start, end + 1))
    selected = sorted(set(selected))
    if not selected:
        raise ValueError("Select the exposed target β-strand residues to pair with the binder.")
    return selected


def _selection_groups(numbers: list[int]) -> list[list[int]]:
    groups: list[list[int]] = []
    for number in numbers:
        if not groups or number != groups[-1][-1] + 1:
            groups.append([])
        groups[-1].append(number)
    return groups


def estimate_register_count(binder_length: int, strand_length: int, strand_count: int = 1) -> int:
    """Count ordered non-overlapping placements enumerated by the vendored script."""
    available = binder_length - strand_length + 1
    # The upstream recursion leaves three unconditioned residues between strands.
    available -= (strand_length + 2) * (strand_count - 1)
    if available < strand_count:
        return 0
    # Each selected position can be assigned to a different target strand; the
    # upstream recursion enumerates those assignments in order.
    return math.factorial(strand_count) * math.comb(available, strand_count)


def prepare_request(parameters: dict[str, Any], *, require_files: bool = True) -> dict[str, Any]:
    p = dict(parameters)
    target_path = p.get("target_pdb")
    if not isinstance(target_path, (str, Path)) or not str(target_path).strip():
        raise ValueError("Load a target PDB structure.")
    target_path = Path(target_path).expanduser()
    if require_files and (not target_path.is_file() or target_path.suffix.lower() not in {".pdb", ".ent"}):
        raise ValueError("The β-pairing input must be an existing PDB file.")
    chain = str(p.get("target_chain") or "").strip()
    if not re.fullmatch(r"[A-Za-z]", chain):
        raise ValueError("Choose one target protein chain with a single-letter PDB ID.")
    binder_length = p.get("binder_length", 65)
    strand_length = p.get("binder_strand_length", 6)
    count = p.get("num_designs", 10)
    steps = p.get("timesteps", 50)
    seed = p.get("seed", 0)
    cpu_cores = p.get("cpu_cores", 4)
    for name, value, low, high in (
        ("Binder length", binder_length, 30, 500),
        ("Binder β-strand length", strand_length, 3, 40),
        ("Number of backbones", count, 1, 10000),
        ("Diffusion timesteps", steps, 1, 200),
        ("Random seed", seed, 0, 2147483647),
        ("CPU cores", cpu_cores, 1, 256),
    ):
        if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
            raise ValueError(f"{name} must be a whole number from {low:,} to {high:,}.")
    # ScaffoldedSampler automatically selects Complex_Fold_base, the RFdiffusion
    # checkpoint trained with the SS/block-adjacency input channels used here.
    model = str(p.get("model_weights", "Complex_Fold_base"))
    if model != "Complex_Fold_base":
        raise ValueError("β-pairing requires the Complex_Fold_base checkpoint trained for SS/adjacency conditioning.")
    systematic = p.get("systematic_sampling", False)
    if type(systematic) is not bool:
        raise ValueError("Register sampling must be random or systematic.")

    p.update(
        target_pdb=str(target_path), target_chain=chain, binder_length=binder_length,
        binder_strand_length=strand_length, num_designs=count, timesteps=steps, seed=seed,
        cpu_cores=cpu_cores, model_weights=model, systematic_sampling=systematic,
    )
    if require_files:
        info = inventory(target_path)
        if chain not in info["chains"] or info["chains"][chain]["kind"] != "protein":
            raise ValueError(f"Chain {chain} is not a protein chain in the selected PDB.")
        numbers = parse_residue_selection(str(p.get("target_beta_strands") or ""), chain)
        present = {
            row["number"] for row in info["chains"][chain]["residues"]
            if not row.get("insertion_code") and row["kind"] == "protein"
        }
        missing = [number for number in numbers if number not in present]
        if missing:
            raise ValueError("Selected β-strand residues are absent from the target PDB: " + ", ".join(f"{chain}{n}" for n in missing[:8]))
        groups = _selection_groups(numbers)
        if any(len(group) < 3 for group in groups):
            raise ValueError("Each selected target β-strand segment must contain at least three consecutive residues.")
        if len(groups) > 4:
            raise ValueError("Select no more than four target β-strand segments in one job.")
        present_names = {row["name"] for row in info["chains"][chain]["residues"] if row["kind"] == "protein"}
        unsupported = sorted(present_names - AA3)
        if unsupported:
            raise ValueError("The selected structure contains non-standard residues unsupported by the interface tensor script: " + ", ".join(unsupported))
        registers = estimate_register_count(binder_length, strand_length, len(groups))
        if registers < 1:
            raise ValueError("There is not enough binder length to place the requested β-strands with loop spacing.")
        if registers > 10000:
            raise ValueError(f"This selection expands to {registers:,} tensor registers. Reduce the number of paired strands or binder length (maximum 10,000).")
        # Every file contains one square adjacency matrix with three float classes.
        estimated_tensor_bytes = registers * (binder_length + len(info["chains"][chain]["residues"])) ** 2 * 4
        if estimated_tensor_bytes > 2_000_000_000:
            raise ValueError("The interface tensors would require more than 2 GB. Shorten the target or reduce the number of paired strands.")
        p.update(
            target_beta_strands=",".join(f"{chain}{n}" for n in numbers),
            target_strand_groups=[[f"{chain}{n}" for n in group] for group in groups],
            register_count=registers,
        )
    else:
        parse_residue_selection(str(p.get("target_beta_strands") or ""), chain)
    return p
