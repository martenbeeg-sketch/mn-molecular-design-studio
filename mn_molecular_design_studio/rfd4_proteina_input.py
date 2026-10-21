"""Input validation and condition-spec construction for RFDiffusion4-Proteina."""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

import yaml

from mn_molecular_design_studio.design_validation import protein_chains, read_structure, sequence


UPSTREAM_VERSION = "RFDiffusion4-Proteina 0.0.3"
UPSTREAM_COMMIT = "9d35ee4f6328fd9b17f5b84756ee5db99d6d89ac"
IMAGE = "mn-rfd4-proteina:0.0.3-cu130"

TASKS = {
    "protein_binder": "Protein Binder",
    "nanobody": "Nanobody Design",
    "antibody": "Antibody Design",
    "tcr": "TCR Design",
    "motif_scaffolding": "Motif Scaffolding",
    "enzyme": "Enzyme Design",
    "small_molecule_binder": "Small Molecule Binder",
    "dna_binder": "DNA Binder",
    "rna_design": "RNA Design",
    "existing_structure": "Existing Structure",
    "custom_conditions": "Custom Conditions",
    "unconditional": "Unconditional Generation",
}

TASK_HELP = {
    "protein_binder": "Co-design a binder backbone and sequence around a selected protein target and optional hotspot residues.",
    "nanobody": "Use an RFDiffusion4-Proteina annotated input or native condition specification for VHH design.",
    "antibody": "Use an annotated antibody input or a native condition specification for antibody design.",
    "tcr": "Use an annotated input or native condition specification for TCR design.",
    "motif_scaffolding": "Scaffold a fixed structural motif with an annotated CIF or native condition specification.",
    "enzyme": "Build a scaffold around catalytic residues and any retained ligand using an annotated CIF or native specification.",
    "small_molecule_binder": "Design a protein pocket around a ligand using an annotated CIF or native condition specification.",
    "dna_binder": "Condition on DNA structure with an annotated CIF or native condition specification.",
    "rna_design": "Condition RNA design on an annotated CIF or native condition specification.",
    "existing_structure": "Refine or diversify a previous structure with an annotated CIF or native condition specification.",
    "custom_conditions": "Provide a native RFDiffusion4-Proteina YAML or JSON condition specification.",
    "unconditional": "Generate a monomer from the model's bundled short-chain template, or provide a custom specification.",
}


def store_structure_uploads(uploaded: list[Any], owner_id: str) -> list[Path]:
    """Store CIF inputs safely under the user's app-managed input directory."""
    from mn_molecular_design_studio.auth import user_input_root

    if len(uploaded) > 128:
        raise ValueError("Upload at most 128 structure files.")
    destinations: list[Path] = []
    names: set[str] = set()
    total_bytes = 0
    for item in uploaded:
        name = Path(str(item.name)).name
        suffix = Path(name).suffix.lower()
        if suffix not in {".cif", ".mmcif"}:
            raise ValueError("RFD4 native condition specifications require CIF/mmCIF inputs. The guided binder form also accepts PDB.")
        data = item.getvalue()
        total_bytes += len(data)
        if len(data) > 50_000_000 or total_bytes > 100_000_000:
            raise ValueError("Structure uploads must be under 50 MB each and 100 MB total.")
        if name in names:
            raise ValueError(f"Uploaded CIF filenames must be unique: {name}")
        names.add(name)
        digest = hashlib.sha256(data).hexdigest()[:16]
        destination = user_input_root(owner_id) / "rfd4-proteina" / digest / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        if not destination.exists():
            destination.write_bytes(data)
        destinations.append(destination)
    return destinations


def _as_safe_yaml(payload: object) -> str:
    return yaml.safe_dump(payload, sort_keys=False, allow_unicode=True, width=1000)


def _cif_label_mapping(path: Path) -> tuple[dict[str, str], dict[tuple[str, int], int]]:
    """Map author chain/residue IDs to the label IDs RFD4 selectors consume."""
    import gemmi

    document = gemmi.cif.read_file(str(path))
    if not len(document):
        raise ValueError("Could not read the normalized mmCIF target.")
    block = document.sole_block()
    columns = ["group_PDB", "label_asym_id", "auth_asym_id", "label_seq_id", "auth_seq_id", "pdbx_PDB_ins_code"]
    rows = block.find("_atom_site.", columns)
    chains: dict[str, str] = {}
    residues: dict[tuple[str, int], int] = {}
    ambiguous: set[tuple[str, int]] = set()
    for row in rows:
        if row.str(0).upper() != "ATOM":
            continue
        label_chain, author_chain = row.str(1).strip(), row.str(2).strip()
        if label_chain in {"", ".", "?"} or author_chain in {"", ".", "?"}:
            continue
        chains[author_chain] = label_chain
        try:
            label_number, author_number = int(row.str(3)), int(row.str(4))
        except ValueError:
            continue
        insertion = row.str(5).strip()
        key = (author_chain, author_number)
        if insertion not in {"", ".", "?"}:
            ambiguous.add(key)
        previous = residues.get(key)
        if previous is not None and previous != label_number:
            ambiguous.add(key)
        else:
            residues[key] = label_number
    for key in ambiguous:
        residues.pop(key, None)
    if not chains or not residues:
        raise ValueError("The input did not contain addressable protein residues after conversion to mmCIF.")
    return chains, residues


def _write_guided_binder_input(
    source: Path,
    chains: list[str],
    hotspots: str,
    length_min: int,
    length_max: int,
    input_root: Path,
) -> tuple[Path, Path, dict[str, str]]:
    from Bio.PDB import MMCIFIO, Select

    structure = read_structure(source)
    available = protein_chains(structure)
    selected = [str(chain) for chain in chains or []]
    if not selected or any(chain not in available for chain in selected):
        raise ValueError("Select at least one protein target chain present in the uploaded structure.")

    class ProteinChainSelect(Select):
        def accept_chain(self, chain):
            return chain.id in selected

        def accept_residue(self, residue):
            return residue.id[0] == " " and "CA" in residue

    input_root.mkdir(parents=True, exist_ok=True)
    target_cif = input_root / "target_for_design.cif"
    writer = MMCIFIO()
    writer.set_structure(structure)
    writer.save(str(target_cif), select=ProteinChainSelect())
    label_chains, label_residues = _cif_label_mapping(target_cif)
    if any(chain not in label_chains for chain in selected):
        raise ValueError("Could not map selected author chain IDs to RFD4 mmCIF chain IDs.")

    target_segments = []
    for author_chain in selected:
        label_chain = label_chains[author_chain]
        residue_ids = sorted({
            number for (chain, _author_number), number in label_residues.items()
            if chain == author_chain
        })
        if not residue_ids:
            raise ValueError(f"Selected target chain {author_chain!r} has no addressable residues.")
        runs: list[list[int]] = []
        for number in residue_ids:
            if runs and number == runs[-1][-1] + 1:
                runs[-1].append(number)
            else:
                runs.append([number])
        target_segments.append(", ".join(
            f"{label_chain}{run[0]}-{run[-1]}" if len(run) > 1 else f"{label_chain}{run[0]}"
            for run in runs
        ))
    target_contig = ", /0, ".join(target_segments)
    binder_range = f"{length_min}-{length_max}"
    conditions: dict[str, Any] = {
        "C_CRDGRP": {1: [{"select": f"sel('{label_chains[chain]}/*/*/[N,CA,C,O]')"} for chain in selected]},
        "C_SEQ": {True: [{"select": f"sel('{label_chains[chain]}')"} for chain in selected]},
    }
    hotspot_selections = []
    for token in (part.strip() for part in str(hotspots or "").split(",")):
        if not token:
            continue
        match = re.fullmatch(r"([A-Za-z0-9_]+)(-?\d+)", token)
        if not match:
            raise ValueError(f"Could not parse hotspot {token!r}; use chain and author residue IDs, such as A56.")
        author_chain, author_number = match.group(1), int(match.group(2))
        if author_chain not in selected:
            raise ValueError(f"Hotspot {token!r} is not on a selected target chain.")
        label_number = label_residues.get((author_chain, author_number))
        if label_number is None:
            raise ValueError(f"Hotspot {token!r} has an insertion code or cannot be mapped unambiguously.")
        hotspot_selections.append({"select": f"sel('{label_chains[author_chain]}/*/{label_number}')"})
    if hotspot_selections:
        conditions["C_HOT"] = {True: hotspot_selections}

    spec = {
        "protein_binder": {
            "input": target_cif.name,
            "contig": f"{binder_range}, /0, {target_contig}",
            "conditions": conditions,
        }
    }
    spec_path = input_root / "condition_spec.yaml"
    spec_path.write_text(_as_safe_yaml(spec))
    target_sequences = {chain: sequence(available[chain]) for chain in selected}
    return spec_path, target_cif, target_sequences


def prepare_request(parameters: dict[str, Any], *, run_dir: Path) -> dict[str, Any]:
    """Validate and stage a request, without starting Docker or loading the model."""
    task = str(parameters.get("task") or "protein_binder")
    mode = str(parameters.get("input_mode") or "guided_binder")
    if task not in TASKS:
        raise ValueError("Choose one of the supported RFDiffusion4-Proteina design tasks.")
    if mode not in {"guided_binder", "annotated_cif", "condition_spec", "unconditional"}:
        raise ValueError("Select a supported RFD4 input mode.")
    if mode == "unconditional" and task != "unconditional":
        raise ValueError("The bundled monomer template can only be used with Unconditional Generation.")
    for field, default, low, high in (
        ("num_designs", 8, 1, 10000), ("batch_size", 1, 1, 16), ("nsteps", 60, 1, 200),
        ("cpu_cores", 8, 1, 256), ("seed", 0, 0, 2_147_483_647),
    ):
        value = parameters.get(field, default)
        if type(value) is not int or not low <= value <= high:
            raise ValueError(f"{field} must be an integer between {low} and {high}.")
    length_min = parameters.get("length_min", 60)
    length_max = parameters.get("length_max", 120)
    if type(length_min) is not int or type(length_max) is not int or not (1 <= length_min <= length_max <= 3000):
        raise ValueError("Binder length must be an ordered range from 1 to 3,000 residues.")

    root = run_dir / "artifacts" / "rfd4_proteina"
    input_root = root / "input"
    input_root.mkdir(parents=True, exist_ok=True)
    original_files = [Path(p).expanduser().resolve() for p in parameters.get("input_files") or []]
    target = Path(parameters["target_structure"]).expanduser().resolve() if parameters.get("target_structure") else None
    if target is not None and target not in original_files:
        original_files.append(target)
    if len(original_files) > 128:
        raise ValueError("Provide no more than 128 input structures.")
    staged: dict[str, Path] = {}
    total_bytes = 0
    for source in original_files:
        if not source.is_file() or source.suffix.lower() not in {".cif", ".mmcif", ".pdb", ".ent"}:
            raise ValueError(f"Input structure is missing or has an unsupported format: {source.name}")
        size = source.stat().st_size
        total_bytes += size
        if size > 50_000_000 or total_bytes > 100_000_000:
            raise ValueError("Input structures must be under 50 MB each and 100 MB total.")
        if source.suffix.lower() in {".pdb", ".ent"} and mode != "guided_binder":
            raise ValueError("PDB input is supported by the guided binder form. Native annotated inputs must be CIF/mmCIF.")
        if source.name in staged and staged[source.name] != source:
            raise ValueError(f"Input filenames must be unique: {source.name}")
        destination = input_root / source.name
        if source.resolve() != destination.resolve():
            destination.write_bytes(source.read_bytes())
        staged[source.name] = destination

    target_sequences: dict[str, str] = {}
    if mode == "unconditional":
        spec = {"unconditional": {"input": "/opt/rfd4/data/benchmarks/monomer/monomer-short.cif.zst", "input_as_annotated_cif": True}}
        spec_path = input_root / "condition_spec.yaml"
        spec_path.write_text(_as_safe_yaml(spec))
    elif mode == "guided_binder":
        if task != "protein_binder":
            raise ValueError("The guided target-and-hotspot form is for Protein Binder. Use annotated CIF or condition-spec mode for other tasks.")
        if target is None:
            raise ValueError("Load a target structure before submitting a protein binder design.")
        spec_path, target_cif, target_sequences = _write_guided_binder_input(
            target, list(parameters.get("target_chains") or []), str(parameters.get("hotspots") or ""),
            length_min, length_max, input_root,
        )
        staged[target_cif.name] = target_cif
    elif mode == "annotated_cif":
        files = [Path(p).expanduser().resolve() for p in parameters.get("input_files") or []]
        if not files and target is not None:
            files = [target]
        if task == "unconditional" and not files:
            spec = {"unconditional": {"input": "/opt/rfd4/data/benchmarks/monomer/monomer-short.cif.zst", "input_as_annotated_cif": True}}
            spec_path = input_root / "condition_spec.yaml"
            spec_path.write_text(_as_safe_yaml(spec))
        else:
            if not files:
                raise ValueError("Upload at least one annotated CIF/mmCIF structure.")
            spec = {}
            for index, source in enumerate(files, 1):
                if source.suffix.lower() not in {".cif", ".mmcif"}:
                    raise ValueError("Annotated RFDiffusion4-Proteina inputs must be CIF/mmCIF files.")
                if source.name not in staged:
                    destination = input_root / source.name
                    destination.write_bytes(source.read_bytes())
                    staged[source.name] = destination
                key = re.sub(r"[^A-Za-z0-9_-]+", "_", source.stem).strip("_-") or f"design_{index}"
                while key in spec:
                    key = f"{key}_{index}"
                spec[key] = {"input": source.name, "input_as_annotated_cif": True}
            spec_path = input_root / "condition_spec.yaml"
            spec_path.write_text(_as_safe_yaml(spec))
    else:
        text = str(parameters.get("condition_spec") or "").strip()
        if not text:
            raise ValueError("Paste a YAML/JSON condition spec for this input mode.")
        try:
            spec = json.loads(text) if text.startswith(("{", "[")) else yaml.safe_load(text)
        except (ValueError, yaml.YAMLError) as exc:
            raise ValueError(f"The condition specification is not valid YAML/JSON: {exc}") from exc
        if not isinstance(spec, dict) or not spec:
            raise ValueError("The condition specification must be a non-empty mapping of design names to conditions.")
        if len(spec) > 128:
            raise ValueError("A condition specification can define at most 128 named designs.")
        for name, design in spec.items():
            if name == "global_args":
                continue
            if not isinstance(design, dict):
                raise ValueError(f"Condition entry {name!r} must be a mapping.")
            references = design.get("input")
            reference_list = references if isinstance(references, list) else [references]
            if not reference_list or any(not isinstance(ref, str) for ref in reference_list):
                raise ValueError(f"Condition entry {name!r} must reference uploaded CIF/mmCIF input files.")
            rewritten = []
            for ref in reference_list:
                basename = Path(ref).name
                staged_path = staged.get(basename)
                if staged_path is None:
                    raise ValueError(f"Condition entry {name!r} refers to {ref!r}; upload that CIF and reference its filename.")
                if staged_path.suffix.lower() not in {".cif", ".mmcif"}:
                    raise ValueError("Custom condition specifications require CIF/mmCIF inputs.")
                rewritten.append("/work/" + staged_path.relative_to(run_dir).as_posix())
            design["input"] = rewritten if isinstance(references, list) else rewritten[0]
        spec_path = input_root / "condition_spec.yaml"
        spec_path.write_text(_as_safe_yaml(spec))

    # Replace local filenames with deterministic paths inside the container.
    if mode in {"guided_binder", "annotated_cif"}:
        spec_payload = yaml.safe_load(spec_path.read_text())
        for name, design in spec_payload.items():
            if name == "global_args" or not isinstance(design, dict):
                continue
            reference = design.get("input")
            if isinstance(reference, str) and reference in staged:
                design["input"] = "/work/" + staged[reference].relative_to(run_dir).as_posix()
        spec_path.write_text(_as_safe_yaml(spec_payload))

    spec_payload = yaml.safe_load(spec_path.read_text())
    design_count = sum(1 for name in spec_payload if name != "global_args")
    if not design_count:
        raise ValueError("The condition specification must include at least one named design.")
    if design_count * parameters.get("num_designs", 8) > 10_000:
        raise ValueError("The total request may generate at most 10,000 structures across all named inputs.")

    return {
        "root": root,
        "input_root": input_root,
        "spec_path": spec_path,
        "staged_files": list(staged.values()),
        "target_sequences": target_sequences,
        "target_structure": str(target) if target else None,
        "target_chains": list(parameters.get("target_chains") or []),
        "design_count": design_count,
    }
