"""Verified OpenFold3 confidence evidence for stored RFD4 refolds."""
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

from mn_molecular_design_studio.engine.workflow_engine.workflows.rfd4_refold_evidence import SCHEMA


def verified_pae_pair(pae, model, refolded_sha256, *, pae_sha256="", oracle_model_sha256=""):
    pae_bytes, model_bytes = Path(pae).read_bytes(), Path(model).read_bytes()
    pae_digest, model_digest = sha256(pae_bytes).hexdigest(), sha256(model_bytes).hexdigest()
    payload = json.loads(pae_bytes)
    if (not isinstance(payload, dict) or payload.get("schema") != SCHEMA
            or payload.get("refolded_sha256") != refolded_sha256
            or payload.get("oracle_model_sha256") != model_digest
            or (pae_sha256 and pae_sha256 != pae_digest)
            or (oracle_model_sha256 and oracle_model_sha256 != model_digest)):
        raise ValueError("OpenFold3 PAE does not match the saved model pair.")
    return payload, pae_bytes, model_bytes


def _protein_tokens(cif):
    """Read residue names and label IDs in the C-alpha order used by ipSAE."""
    import gemmi
    from Bio.Data.PDBData import protein_letters_3to1

    block = gemmi.cif.read_file(str(cif)).sole_block()
    atoms = block.find("_atom_site.", [
        "label_atom_id", "label_comp_id", "label_asym_id", "label_seq_id", "?auth_asym_id",
    ])
    proteins, ordered = {}, []
    seen = set()
    for row in atoms:
        if row.str(0) != "CA" or row.str(1) not in protein_letters_3to1:
            continue
        token = (row.str(2), int(row.str(3)))
        if token in seen:
            raise ValueError("Duplicate protein residue ID in the OpenFold3 scoring structure.")
        seen.add(token)
        chain = row.str(4) if row.has(4) else token[0]
        proteins.setdefault(chain, []).append((token, row.str(1)))
        ordered.append(token)
    if not ordered:
        raise ValueError("PAE scoring requires protein residues with explicit mmCIF label IDs.")
    return proteins, ordered


def write_aligned_pae(refolded, pae, model, destination, *, pae_sha256="", oracle_model_sha256=""):
    """Reindex both PAE axes onto the verified canonical protein coordinate order."""
    import numpy as np
    from mn_molecular_design_studio.cache_provenance import file_digest

    payload, _, _ = verified_pae_pair(pae, model, file_digest(refolded),
        pae_sha256=pae_sha256, oracle_model_sha256=oracle_model_sha256)
    raw_chains, _ = _protein_tokens(model)
    refolded_chains, output_order = _protein_tokens(refolded)
    try:
        token_order = list(zip(map(str, payload["token_chain_ids"]),
                               map(int, payload["token_res_ids"]), strict=True))
        matrix = np.asarray(payload["pae"], dtype=float)
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("OpenFold3 confidence is missing valid PAE token identifiers.") from exc
    indices = {token: i for i, token in enumerate(token_order)}
    if (len(indices) != len(token_order) or matrix.shape != (len(indices), len(indices))
            or not np.isfinite(matrix).all() or raw_chains.keys() != refolded_chains.keys()):
        raise ValueError("OpenFold3 PAE cannot be mapped unambiguously to this protein refold.")
    mapped = {}
    for chain, residues in refolded_chains.items():
        raw_residues = raw_chains[chain]
        if [name for _, name in residues] != [name for _, name in raw_residues]:
            raise ValueError("OpenFold3 and canonical refold residue sequences differ.")
        for (canonical, _), (native, _) in zip(residues, raw_residues, strict=True):
            if native not in indices:
                raise ValueError("OpenFold3 PAE is missing a protein residue token.")
            mapped[canonical] = indices[native]
    if len(mapped) != len(indices):
        raise ValueError("The OpenFold3 matrix includes tokens outside the protein scoring representation.")
    order = [mapped[token] for token in output_order]
    destination = Path(destination)
    destination.write_text(json.dumps({"pae": matrix[np.ix_(order, order)].tolist(),
        "token_chain_ids": [chain for chain, _ in output_order],
        "token_res_ids": [residue for _, residue in output_order]},
        separators=(",", ":"), allow_nan=False))
    return destination
