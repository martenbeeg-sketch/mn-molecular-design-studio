"""Validated Kejia threading inputs and installed template-library discovery."""
from __future__ import annotations

import math
from pathlib import Path

from mn_molecular_design_studio.design_validation import protein_chains, read_structure
from mn_molecular_design_studio.runtime import reference_root

SOURCE_REVISION = "7f65c18cc57b3e72d290c885c540276c7d9c6ab8"
IMAGE = "mn-kejia-idr:7f65c18"
TITLE = "Kejia IDR Binder"
CODE_VERSION = "Kejia 7f65c18 · Studio threading adapter v1"
EXAMPLE_SEQUENCE = "FNCREWCWN"
LIBRARIES = ("walle", "walle2", "mini", "strand")


def source_root() -> Path:
    return reference_root() / "kejia-idr" / "source"


def example_template() -> Path:
    return source_root() / "templates/walle/conf_2.pdb"


def normalize_sequence(value: str) -> str:
    lines = str(value).strip().splitlines()
    if sum(line.strip().startswith(">") for line in lines) > 1:
        raise ValueError("Enter one target fragment, or a single-record FASTA.")
    sequence = "".join("".join(line.split()) for line in lines if not line.strip().startswith(">")).upper()
    invalid = sorted(set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"))
    if invalid:
        raise ValueError("Use the 20 standard amino-acid letters only; invalid characters: " + ", ".join(invalid))
    if not 4 <= len(sequence) <= 40:
        raise ValueError("The target fragment must contain 4–40 residues.")
    return sequence


def template_chains(path: Path) -> list[str]:
    if not path.is_file() or path.suffix.lower() not in {".pdb", ".ent", ".cif", ".mmcif"}:
        raise ValueError("Choose an existing PDB or mmCIF template.")
    if path.stat().st_size > 50_000_000:
        raise ValueError("Template structures must be smaller than 50 MB.")
    structure = read_structure(path)
    chains = protein_chains(structure)
    model = next(structure.get_models())
    if len(chains) != 2 or len(list(model.get_chains())) != 2:
        raise ValueError("A Kejia template must contain exactly two protein chains: binder and peptide.")
    for residues in chains.values():
        if any(not {"N", "CA", "C", "O"} <= {a.name for a in residue} for residue in residues):
            raise ValueError("Each template residue needs complete N, CA, C and O backbone atoms.")
    return list(chains)


def library_counts() -> dict[str, int]:
    return {name: len(list((source_root() / "templates" / name).glob("*.pdb"))) for name in LIBRARIES}


def prepare_request(parameters: dict, *, require_files: bool = True) -> dict:
    p = dict(parameters)
    p["target_sequence"] = normalize_sequence(p.get("target_sequence", ""))
    library = p.get("template_library", "upload")
    if library not in {"upload", *LIBRARIES}:
        raise ValueError("Choose an uploaded template or a supported Kejia library.")
    p["template_library"] = library
    for key, default, low, high in (("max_templates", 10, 1, 1000), ("flex_docking_repeats", 3, 1, 100),
                                   ("cpu_cores", 4, 1, 256), ("seed", 0, 0, 2147483647)):
        value = p.get(key, default)
        if type(value) is not int or not low <= value <= high:
            raise ValueError(f"{key} must be an integer between {low} and {high}.")
        p[key] = value
    for key, default, low, high in (("min_seq_identity", 0.0, 0., 1.), ("min_blosum", -10., -100., 100.),
                                   ("max_farep", 10., 0., 10000.)):
        value = p.get(key, default)
        if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or not low <= value <= high:
            raise ValueError(f"{key} must be a finite number between {low} and {high}.")
        p[key] = float(value)
    p["templates"] = []
    if library == "upload":
        if not p.get("template_structure"):
            raise ValueError("Upload a template structure, or select an installed library.")
        if require_files:
            template = Path(p["template_structure"]).expanduser().resolve()
            chains = template_chains(template)
            binder, peptide = p.get("binder_chain", "A"), p.get("peptide_chain", "B")
            if binder == peptide or {binder, peptide} != set(chains):
                raise ValueError("Select distinct binder and peptide chains present in the template.")
            p["templates"] = [{"path": str(template), "binder_chain": binder, "peptide_chain": peptide}]
    elif p.get("template_structure"):
        raise ValueError("Use an uploaded template or a library, not both.")
    elif require_files:
        head = source_root() / ".git/HEAD"
        if not head.is_file() or head.read_text().strip() != SOURCE_REVISION:
            raise ValueError("Install the pinned Kejia template library with scripts/kejia-idr/install.sh.")
        templates = sorted((source_root() / "templates" / library).glob("*.pdb"))[:p["max_templates"]]
        if not templates:
            raise ValueError(f"The {library} template library is not installed.")
        p["templates"] = [{"path": str(path), "binder_chain": "A", "peptide_chain": "B"} for path in templates]
    return p
