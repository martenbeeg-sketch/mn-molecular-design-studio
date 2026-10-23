#!/usr/bin/env python3
"""Filter BoltzGen generated MIMOSA backbones by peptide C–N geometry."""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path


def _model(path: Path):
    from Bio.PDB import MMCIFParser, PDBParser

    parser = MMCIFParser(QUIET=True) if path.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    return next(parser.get_structure(path.stem, str(path)).get_models())


def _design_chains(path: Path) -> set[str]:
    """Map BoltzGen's token design mask onto output polymer chains by order."""
    import numpy as np

    sidecar = path.with_suffix(".npz")
    if not sidecar.is_file():
        raise ValueError(f"Missing BoltzGen token metadata needed to identify designed residues: {sidecar}")
    with np.load(sidecar, allow_pickle=False) as archive:
        if "design_mask" not in archive.files:
            raise ValueError(f"Missing design_mask in BoltzGen token metadata: {sidecar}")
        design_mask = np.asarray(archive["design_mask"]).astype(bool).reshape(-1)
    model = _model(path)
    chains = [
        (chain.id, [r for r in chain if r.id[0] == " "])
        for chain in model
        if any(r.id[0] == " " for r in chain)
    ]
    residue_count = sum(len(residues) for _, residues in chains)
    if residue_count != len(design_mask):
        raise ValueError(
            f"Cannot map BoltzGen design mask of {len(design_mask)} tokens to {residue_count} "
            f"output polymer residues in {path}; refusing to infer binder by chain letter."
        )
    designed = set()
    cursor = 0
    for chain_id, residues in chains:
        if design_mask[cursor:cursor + len(residues)].any():
            designed.add(str(chain_id))
        cursor += len(residues)
    if not designed:
        raise ValueError(f"BoltzGen design_mask marks no designed polymer chain in {path}.")
    return designed


def _integrity(path: Path, binder_chains: set[str]) -> dict:
    model = _model(path)
    checked, missing, bad = 0, 0, []
    for chain_id in sorted(binder_chains):
        residues = [r for r in model[chain_id] if r.id[0] == " "]
        for left, right in zip(residues, residues[1:]):
            if "C" not in left or "N" not in right:
                missing += 1
                continue
            distance = float(left["C"] - right["N"])
            checked += 1
            if abs(distance - 1.329) > 0.5:
                bad.append({
                    "chain": chain_id,
                    "left_residue_id": f"{left.id[1]}{str(left.id[2]).strip()}",
                    "right_residue_id": f"{right.id[1]}{str(right.id[2]).strip()}",
                    "distance_angstrom": distance,
                })
    return {
        "status": "rejected" if bad else "passed" if checked else "not_checked",
        "checked_bond_count": checked,
        "missing_atom_pair_count": missing,
        "bad_bond_count": len(bad),
        "bad_bonds": bad,
    }


def main() -> int:
    generated_dir, _target_path, report_path = map(Path, sys.argv[1:4])
    rejected_dir = generated_dir.parent / "rejected_backbones"
    candidates = []
    for path in sorted(generated_dir.glob("*.cif")):
        binder_chains = _design_chains(path)
        result = _integrity(path, binder_chains)
        record = {"file": path.name, "binder_chains": sorted(binder_chains), **result}
        candidates.append(record)
        if result["status"] != "passed":
            rejected_dir.mkdir(parents=True, exist_ok=True)
            shutil.move(str(path), rejected_dir / path.name)
            sidecar = path.with_suffix(".npz")
            if sidecar.exists():
                shutil.move(str(sidecar), rejected_dir / sidecar.name)
    report = {
        "schema": "mn-mimosa-integrity.v1",
        "criterion": "abs(C-N distance - 1.329 Å) <= 0.5 Å for consecutive modeled binder polymer residues",
        "candidate_count": len(candidates),
        "rejected_count": sum(item["status"] != "passed" for item in candidates),
        "candidates": candidates,
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"candidate_count": report["candidate_count"], "rejected_count": report["rejected_count"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
