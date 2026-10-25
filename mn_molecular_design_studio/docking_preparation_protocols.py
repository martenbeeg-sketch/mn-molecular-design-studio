"""Shared, reproducible MolScrub routes used before AutoDock-family docking."""

from __future__ import annotations

import math


DEFAULT_SCRUB_PROTOCOL = "molscrub_0_1_0"

SCRUB_PROTOCOLS = {
    "molscrub_0_1_0": {
        "label": "MolScrub 0.1.0 defaults · fastest benchmark",
        "version": "0.1.0",
        "conformers": 1,
        "ff_iterations": 200,
        "conformer_option": "--numconfs",
        "executable": "/opt/molscrub/venvs/0.1.0/bin/scrub.py",
        "benchmark": "10,000-record preparation benchmark (Scrub + parallel Meeko, docking excluded): 94.8 s; 104.5 usable PDBQT/s; 99.10% success.",
    },
    "molscrub_0_3_0_defaults": {
        "label": "MolScrub 0.3.0 defaults · 3 conformers, 400 FF iterations",
        "version": "0.3.0",
        "conformers": 3,
        "ff_iterations": 400,
        "conformer_option": "--num_internal_confs",
        "executable": "/opt/molscrub/venvs/0.3.0/bin/scrub.py",
        "benchmark": "10,000-record preparation benchmark (Scrub + parallel Meeko, docking excluded): 317.8 s; 31.2 usable PDBQT/s; 99.13% success.",
    },
    "molscrub_0_3_0_matched": {
        "label": "MolScrub 0.3.0 matched to 0.1.0 · 1 conformer, 200 FF iterations",
        "version": "0.3.0",
        "conformers": 1,
        "ff_iterations": 200,
        "conformer_option": "--num_internal_confs",
        "executable": "/opt/molscrub/venvs/0.3.0/bin/scrub.py",
        "benchmark": "10,000-record preparation benchmark (Scrub + parallel Meeko, docking excluded): 201.8 s; 49.1 usable PDBQT/s; 99.12% success.",
    },
}


def normalize_scrub_protocol(value: object = None) -> str:
    protocol = str(value or DEFAULT_SCRUB_PROTOCOL).strip()
    if protocol not in SCRUB_PROTOCOLS:
        raise ValueError(f"Unknown ligand-preparation protocol: {protocol}.")
    return protocol


def scrub_invocation(protocol: object = None, *, ph: float = 7.4,
                     skip_tautomers: bool = True, cpu: int = 1) -> dict:
    """Return the pinned MolScrub executable and CLI flags for one route."""
    selected = normalize_scrub_protocol(protocol)
    settings = SCRUB_PROTOCOLS[selected]
    ph = float(ph)
    if not math.isfinite(ph) or not 0 <= ph <= 14:
        raise ValueError("Scrub pH must be between 0 and 14.")
    if type(skip_tautomers) is not bool:
        raise ValueError("Skip-tautomer setting must be a boolean.")
    if type(cpu) is not int or cpu < 1:
        raise ValueError("Scrub CPU allocation must be a positive integer.")

    arguments = ["--ph", f"{ph:.2f}", "--cpu", str(cpu), "--ff", "mmff94s",
                 "--max_ff_iter", str(settings["ff_iterations"]),
                 settings["conformer_option"], str(settings["conformers"])]
    if skip_tautomers:
        # Both tested releases use the plural CLI option.
        arguments.append("--skip_tautomers")
    return {
        "protocol": selected,
        "version": settings["version"],
        "label": settings["label"],
        "executable": settings["executable"],
        "python_executable": f"/opt/molscrub/venvs/{settings['version']}/bin/python",
        "arguments": arguments,
        "conformers": settings["conformers"],
        "ff_iterations": settings["ff_iterations"],
    }
