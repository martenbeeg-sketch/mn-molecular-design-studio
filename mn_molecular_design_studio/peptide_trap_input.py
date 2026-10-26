"""Author-defined peptide-trap scaffolds and burial-pattern sequence matching."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from mn_molecular_design_studio.design_validation import protein_chains, read_structure
from mn_molecular_design_studio.runtime import reference_root

TITLE = "Peptide Trap"
IMAGE = "mn-peptide-trap:sahtoe-10391229"
CODE_VERSION = "Sahtoe 2024 · Zenodo 10391229 · Studio adapter v1"
ARCHIVE_SHA256 = "72e228dda952b82a24c1c9c4403bcc2d4340d6a25f4f06cb148954352d9b67c3"
SCAFFOLDS = {
    "C37": {"label": "C37 · 8-residue β-strand cleft", "file": "C37.pdb", "pattern": "10101011"},
    "C104.1": {"label": "C104.1 model · 10-residue β-strand cleft", "file": "C104.1.pdb", "pattern": "1110101011"},
}
# The deposited seqthread.py uses these exact allowed identities and masks.
ALLOWED_BURIED = frozenset("GFAMILVWYST")
EXAMPLE_SEQUENCE = "DVRFNFRE"


def library_root() -> Path:
    return reference_root() / "peptide-traps/sahtoe-10391229"


def template_path(scaffold: str) -> Path:
    return library_root() / "source/pdb_models" / SCAFFOLDS[scaffold]["file"]


def normalize_sequence(value: str) -> str:
    lines = str(value).strip().splitlines()
    if sum(line.strip().startswith(">") for line in lines) > 1:
        raise ValueError("Use one target sequence or a single-record FASTA.")
    text = "".join("".join(line.split()) for line in lines if not line.strip().startswith(">")).upper()
    if not 8 <= len(text) <= 5000:
        raise ValueError("Provide a target sequence of 8–5,000 residues; matching windows are 8 or 10 residues long.")
    if set(text) - set("ACDEFGHIKLMNPQRSTVWY"):
        raise ValueError("Use only the 20 standard amino-acid letters.")
    return text


def match_registers(target: str, scaffold_ids: list[str]) -> list[dict]:
    """Mirror deposited burial masks, including the internal-proline exclusion."""
    rows = []
    for scaffold in scaffold_ids:
        mask = SCAFFOLDS[scaffold]["pattern"]
        for offset in range(len(target) - len(mask) + 1):
            fragment = target[offset:offset + len(mask)]
            if "P" in fragment[1:-1]:
                continue
            if any(flag == "0" and residue not in ALLOWED_BURIED for residue, flag in zip(fragment, mask)):
                continue
            rows.append({"scaffold": scaffold, "start": offset + 1, "end": offset + len(mask),
                         "sequence": fragment, "burial_pattern": mask})
    return rows


def validate_library(scaffold_ids: list[str]) -> dict:
    root = library_root()
    try:
        manifest = json.loads((root / "source_manifest.json").read_text())
    except (OSError, ValueError) as exc:
        raise ValueError("Install the authors' templates with scripts/peptide-trap/install.py first.") from exc
    if manifest.get("archive_sha256") != ARCHIVE_SHA256:
        raise ValueError("The peptide-trap library is not the pinned authors' deposit.")
    for relative in ["design_after_seqthread.xml", *("pdb_models/" + SCAFFOLDS[key]["file"] for key in scaffold_ids)]:
        path = root / "source" / relative
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != manifest.get("files", {}).get(relative):
            raise ValueError(f"The installed paper asset is missing or modified: {relative}")
    for key in scaffold_ids:
        chains = protein_chains(read_structure(template_path(key)))
        if set(chains) != {"A", "B"} or len(chains["A"]) != len(SCAFFOLDS[key]["pattern"]):
            raise ValueError(f"{key} does not contain the expected peptide A and binder B chains.")
    return manifest


def prepare_request(parameters: dict, *, require_files: bool = True) -> dict:
    p = dict(parameters)
    p["target_sequence"] = normalize_sequence(p.get("target_sequence", ""))
    selected = p.get("scaffold_ids", ["C37", "C104.1"])
    if not isinstance(selected, list) or not selected or any(key not in SCAFFOLDS for key in selected) or len(set(selected)) != len(selected):
        raise ValueError("Select one or both curated paper scaffolds without duplicates.")
    p["scaffold_ids"] = list(selected)
    for key, default, low, high in (("max_registers", 20, 1, 1000), ("designs_per_register", 1, 1, 100),
                                   ("design_repeats", 2, 1, 5), ("relax_repeats", 2, 1, 5),
                                   ("cpu_cores", 4, 1, 256), ("seed", 0, 0, 2147483647)):
        value = p.get(key, default)
        if type(value) is not int or not low <= value <= high:
            raise ValueError(f"{key} must be an integer from {low} to {high}.")
        p[key] = value
    value = p.get("max_backbone_rmsd", 2.0)
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or not 0.1 <= value <= 10:
        raise ValueError("Maximum backbone RMSD must be a finite value from 0.1 to 10 Å.")
    p["max_backbone_rmsd"] = float(value)
    p["all_matches"] = match_registers(p["target_sequence"], p["scaffold_ids"])
    if not p["all_matches"]:
        raise ValueError("No target windows match the selected scaffolds' buried-residue patterns. Select another scaffold or target region.")
    p["registers"] = p["all_matches"][:p["max_registers"]]
    if require_files:
        validate_library(p["scaffold_ids"])
    return p
