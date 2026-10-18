#!/usr/bin/env python
from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import time
from pathlib import Path

import jax
import numpy as np
from colabdesign import mk_af_model
from colabdesign.af.alphafold.common import protein
from colabdesign.af.alphafold.data import parsers, pipeline


METRICS = {
    "rmsd": "target_aligned_binder_rmsd",
    "plddt": "binder_plddt",
    "pae": "binder_pae",
    "ptm": "ptm",
    "con": "con_loss",
    "i_pae": "ipae",
    "i_ptm": "iptm",
    "i_con": "i_con_loss",
}


def add_cyclic_offset(model, offset_type: int = 2):
    def cyclic_offset(length: int):
        i = np.arange(length)
        ij = np.stack([i, i + length], -1)
        offset = i[:, None] - i[None, :]
        c_offset = np.abs(ij[:, None, :, None] - ij[None, :, None, :]).min((2, 3))
        if offset_type >= 2:
            idx = c_offset < np.abs(offset)
            c_offset[idx] = -c_offset[idx]
        if offset_type == 3:
            idx = np.abs(c_offset) > 2
            c_offset[idx] = (32 * c_offset[idx]) / abs(c_offset[idx])
        return c_offset * np.sign(offset)

    idx = model._inputs["residue_index"]
    offset = np.array(idx[:, None] - idx[None, :])
    offset[model._target_len :, model._target_len :] = cyclic_offset(model._binder_len)
    model._inputs["offset"] = offset


def _inject_af2_msa(inputs):
    """Replace ColabDesign's query-only MSA with staged per-chain A3Ms."""
    if "af2_custom_msa_feat" not in inputs:
        return
    inputs["msa_feat"] = inputs.pop("af2_custom_msa_feat")
    inputs["msa_mask"] = inputs.pop("af2_custom_msa_mask")
    inputs["msa_row_mask"] = inputs.pop("af2_custom_msa_row_mask")


def _attach_af2_msa(model, msa_payload, chain_order, expected_sequences):
    if not chain_order:
        raise ValueError("The AF2 input has no protein chains for its requested MSA.")
    paths = msa_payload.get("paths") or {}
    sequences_by_chain = msa_payload.get("sequences") or {}
    chain_msas = {}
    for chain in chain_order:
        relative_path = paths.get(chain)
        expected = str(expected_sequences.get(chain) or sequences_by_chain.get(chain) or "").upper()
        if not relative_path or not expected:
            raise ValueError(f"The staged AF2 MSA is missing chain {chain} or its sequence.")
        path = Path("/work") / str(relative_path)
        if not path.is_file():
            raise ValueError(f"The staged AF2 A3M file is missing: {relative_path}")
        aligned, deletions = parsers.parse_a3m(path.read_text())
        if not aligned or len(aligned[0]) != len(expected) or aligned[0].replace("-", "") != expected:
            raise ValueError(f"The AF2 A3M query does not match chain {chain}.")
        if any(len(row) != len(expected) for row in aligned):
            raise ValueError(f"The AF2 A3M rows for chain {chain} have inconsistent lengths.")
        chain_msas[chain] = (aligned, deletions)

    expected_length = sum(len(str(expected_sequences.get(chain) or sequences_by_chain.get(chain) or "")) for chain in chain_order)
    actual_length = int(model._inputs["residue_index"].shape[0])
    if expected_length != actual_length:
        raise ValueError(
            f"The staged AF2 MSA spans {expected_length} residues but the prepared prediction input has {actual_length}."
        )

    rows = []
    deletion_rows = []
    query_parts = [chain_msas[chain][0][0] for chain in chain_order]
    query_deletions = [chain_msas[chain][1][0] for chain in chain_order]
    rows.append("".join(query_parts))
    deletion_rows.append([value for part in query_deletions for value in part])
    seen = {rows[0]}
    for chain_index, chain in enumerate(chain_order):
        aligned, deletions = chain_msas[chain]
        prefix_len = sum(len(str(expected_sequences.get(item) or sequences_by_chain.get(item) or "")) for item in chain_order[:chain_index])
        chain_len = len(str(expected_sequences.get(chain) or sequences_by_chain.get(chain) or ""))
        suffix_len = expected_length - prefix_len - chain_len
        for aligned_row, deletion_row in zip(aligned[1:], deletions[1:]):
            row = "-" * prefix_len + aligned_row + "-" * suffix_len
            if row in seen:
                continue
            seen.add(row)
            rows.append(row)
            deletion_rows.append([0] * prefix_len + list(deletion_row) + [0] * suffix_len)

    msa_features = pipeline.make_msa_features([rows], [deletion_rows])
    msa = np.asarray(msa_features["msa"], dtype=np.int32)
    deletion_matrix = np.asarray(msa_features["deletion_matrix_int"], dtype=np.float32)
    if np.any(msa >= 23):
        raise ValueError("The AF2 MSA contains unsupported residue codes.")
    one_hot = np.eye(23, dtype=np.float32)[msa]
    has_deletion = (deletion_matrix > 0).astype(np.float32)[..., None]
    deletion_value = (np.arctan(deletion_matrix / 3.0) * (2.0 / np.pi))[..., None]
    profile = one_hot.mean(axis=0, keepdims=True)
    deletion_mean = deletion_value.mean(axis=0, keepdims=True)
    profile_rows = np.broadcast_to(profile, (msa.shape[0], msa.shape[1], profile.shape[-1]))
    deletion_mean_rows = np.broadcast_to(deletion_mean, (msa.shape[0], msa.shape[1], 1))
    features = np.concatenate(
        [one_hot, has_deletion, deletion_value, profile_rows, deletion_mean_rows],
        axis=-1,
    ).astype(np.float32)
    model._inputs["af2_custom_msa_feat"] = features
    model._inputs["af2_custom_msa_mask"] = np.ones((msa.shape[0], msa.shape[1]), dtype=np.float32)
    model._inputs["af2_custom_msa_row_mask"] = np.ones((msa.shape[0],), dtype=np.float32)
    return {chain: len(chain_msas[chain][0]) for chain in chain_order}


def _read_af2_msa_input(input_dir, basename):
    path = Path(input_dir) / f"{basename}.msa.json"
    if not path.is_file():
        raise ValueError(f"The AF2 MSA input map is missing for {basename}.")
    return json.loads(path.read_text())


def get_pdb_total_length(pdb_path: str) -> int:
    unique_residues = set()
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                unique_residues.add((line[21].strip(), line[22:26].strip(), line[26].strip()))
    return len(unique_residues)


def get_design_sequence(pdb_path: str, chain_id: str, *, target_only: bool = False) -> str:
    sidecar_names = (
        [".target_sequence.txt", ".binder_sequence.txt"]
        if target_only
        else [".binder_sequence.txt"]
    )
    for suffix in sidecar_names:
        sidecar = os.path.splitext(pdb_path)[0] + suffix
        if os.path.exists(sidecar):
            with open(sidecar, "r", encoding="utf-8") as handle:
                return "".join(handle.read().split()).upper()

    chain_ids = [part.strip() for part in str(chain_id or "").split(",") if part.strip()]
    if not chain_ids:
        chain_ids = [str(chain_id or "A")]
    chain_set = set(chain_ids)

    residue_map = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    residues: list[str] = []
    seen: set[tuple[str, str]] = set()
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")) or line[21].strip() not in chain_set:
                continue
            residue_key = (line[21].strip(), line[22:26].strip(), line[26].strip())
            if residue_key in seen:
                continue
            seen.add(residue_key)
            residues.append(residue_map.get(line[17:20].strip().upper(), "X"))
    return "".join(residues)


def get_target_chains(pdb_path: str, fallback: str) -> str:
    sidecar = os.path.splitext(pdb_path)[0] + ".target_chains.txt"
    if os.path.exists(sidecar):
        with open(sidecar, "r", encoding="utf-8") as handle:
            chains = ",".join(
                part.strip()
                for part in handle.read().split(",")
                if part.strip()
            )
        if chains:
            return chains
    return fallback


def get_binder_source_chains(pdb_path: str, fallback: str = "A") -> list[str]:
    stem = os.path.splitext(pdb_path)[0]
    for sidecar in (f"{stem}.binder_output_chains.txt", f"{stem}.binder_source_chains.txt"):
        if os.path.exists(sidecar):
            with open(sidecar, "r", encoding="utf-8") as handle:
                chains = [part.strip() for part in handle.read().split(",") if part.strip()]
            if chains:
                return chains
    return [part.strip() for part in str(fallback or "A").split(",") if part.strip()] or ["A"]


def get_target_source_chains(pdb_path: str, fallback: str = "A") -> list[str]:
    sidecar = os.path.splitext(pdb_path)[0] + ".target_source_chains.txt"
    if not os.path.exists(sidecar):
        return [part.strip() for part in fallback.split(",") if part.strip()] or ["A"]
    with open(sidecar, "r", encoding="utf-8") as handle:
        chains = [part.strip() for part in handle.read().split(",") if part.strip()]
    return chains or [part.strip() for part in fallback.split(",") if part.strip()] or ["A"]


def _chain_ids_from_pdb(pdb_path: str) -> list[str]:
    chains = []
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ATOM"):
                chain = line[21].strip()
                if chain and chain not in chains:
                    chains.append(chain)
    return chains


def _chain_lengths_from_pdb(pdb_path: str, chains: list[str]) -> list[tuple[str, int]]:
    chain_set = set(chains)
    residues: dict[str, set[tuple[str, str]]] = {chain: set() for chain in chains}
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            chain = line[21].strip()
            if chain not in chain_set:
                continue
            residues.setdefault(chain, set()).add((line[22:26].strip(), line[26].strip()))
    return [(chain, len(residues.get(chain, set()))) for chain in chains]


def _target_output_maps(model):
    chain_lengths = getattr(model, "_target_chain_lengths_for_output", None) or [("B", model._target_len)]
    target_residue_indices = [int(value) for value in model.aux["residue_index"][: model._target_len]]
    chain_by_residue: dict[int, str] = {}
    resno_by_residue: dict[int, int] = {}
    offset = 0
    for chain, length in chain_lengths:
        for local_index, residue_index in enumerate(target_residue_indices[offset : offset + int(length)], start=1):
            chain_by_residue[residue_index] = chain
            resno_by_residue[residue_index] = local_index
        offset += int(length)
    return chain_by_residue, resno_by_residue


def _format_target_residue_labels(model, target_residue_indices) -> str:
    chain_by_residue, resno_by_residue = _target_output_maps(model)
    labels = []
    for residue_index in target_residue_indices:
        key = int(residue_index)
        labels.append(f"{chain_by_residue.get(key, 'B')}{resno_by_residue.get(key, key)}")
    return ",".join(labels)


def _ca_coordinates(pdb_path: str, chains: list[str] | None = None) -> np.ndarray:
    coordinates = []
    chain_set = set(chains or [])
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            if chain_set and line[21].strip() not in chain_set:
                continue
            try:
                coordinates.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
    return np.asarray(coordinates, dtype=float)


def _aligned_ca_rmsd(reference_pdb: str, predicted_pdb: str, reference_chains: list[str]) -> float | None:
    reference = _ca_coordinates(reference_pdb, reference_chains)
    predicted = _ca_coordinates(predicted_pdb)
    count = min(len(reference), len(predicted))
    if count < 3:
        return None
    reference = reference[:count]
    predicted = predicted[:count]
    reference_centered = reference - reference.mean(axis=0)
    predicted_centered = predicted - predicted.mean(axis=0)
    u, _s, vt = np.linalg.svd(predicted_centered.T @ reference_centered)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = u @ vt
    aligned = predicted_centered @ rotation
    return float(np.sqrt(np.mean(np.sum((aligned - reference_centered) ** 2, axis=1))))


def _to_numpy(value):
    if value is None:
        return None
    if isinstance(value, np.ndarray):
        return value
    try:
        return np.array(value)
    except Exception:
        return None


def _extract_pae_matrix(model) -> np.ndarray | None:
    aux = getattr(model, "aux", {}) or {}
    candidates = []
    if isinstance(aux, dict):
        candidates.append(aux.get("pae"))
        outputs = aux.get("outputs", {})
        if isinstance(outputs, dict):
            candidates.append(outputs.get("predicted_aligned_error"))
        all_aux = aux.get("all", {})
        if isinstance(all_aux, dict):
            candidates.append(all_aux.get("pae"))
            candidates.append(all_aux.get("predicted_aligned_error"))

    for candidate in candidates:
        arr = _to_numpy(candidate)
        if arr is None:
            continue
        arr = np.squeeze(arr)
        if arr.ndim == 2 and arr.shape[0] == arr.shape[1]:
            arr = arr.astype(float)
            if np.nanmax(arr) <= 1.5:
                arr = arr * 31.0
            return arr
    return None


def _flip_target_binder_matrix(matrix: np.ndarray, target_len: int, binder_len: int) -> np.ndarray:
    order = np.concatenate([np.arange(target_len, target_len + binder_len), np.arange(0, target_len)])
    return matrix[np.ix_(order, order)]


def _save_pae_json(model, out_path: str, metrics: dict):
    pae_matrix = _extract_pae_matrix(model)
    if pae_matrix is None:
        return
    pae_matrix = _flip_target_binder_matrix(pae_matrix, model._target_len, model._binder_len)
    payload = {
        "predicted_aligned_error": pae_matrix.tolist(),
        "pae": pae_matrix.tolist(),
        "max_predicted_aligned_error": 31.0,
    }
    for key in ("ptm", "iptm"):
        if metrics.get(key) is not None:
            payload[key] = metrics[key]
    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle)


def _save_monomer_pae_json(model, out_path: str, metrics: dict):
    pae_matrix = _extract_pae_matrix(model)
    if pae_matrix is None:
        return
    payload = {
        "predicted_aligned_error": pae_matrix.tolist(),
        "pae": pae_matrix.tolist(),
        "max_predicted_aligned_error": 31.0,
    }
    for key in ("ptm", "iptm"):
        if metrics.get(key) is not None:
            payload[key] = metrics[key]
    with open(out_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle)


def save_binder_design_pdb(model, filename: str, *, binder_output_chain: str = "A"):
    aux = model._tmp["best"]["aux"] if "aux" in model._tmp.get("best", {}) else model.aux
    aux = aux["all"]
    payload = {key: aux[key] for key in ["aatype", "residue_index", "atom_positions", "atom_mask"]}
    payload["b_factors"] = 100 * payload["atom_mask"] * aux["plddt"][..., None]

    for key in list(payload):
        assert payload[key].shape[1] == model._target_len + model._binder_len
        payload[key] = np.concatenate([payload[key][:, model._target_len :], payload[key][:, : model._target_len]], axis=1)

    def to_pdb_str(x, model_number=None):
        pdb_str = protein.to_pdb(protein.Protein(**x))
        binder_mapping = dict(zip(x["residue_index"][: model._binder_len], range(1, model._binder_len + 1)))
        target_chain_by_residue, target_resno_by_residue = _target_output_maps(model)
        lines = []
        binder_chain = str(binder_output_chain or "A")[:1]
        for line in pdb_str.splitlines()[1:-2]:
            if line.startswith(("ATOM", "HETATM")):
                resno = int(line[22:26].strip())
                if resno in binder_mapping:
                    lines.append(line[:21] + binder_chain + str(binder_mapping[resno]).rjust(4) + line[26:])
                else:
                    chain = target_chain_by_residue.get(resno, "B")[:1]
                    mapped_resno = target_resno_by_residue.get(resno, resno)
                    lines.append(line[:21] + chain + str(mapped_resno).rjust(4) + line[26:])
        body = "\n".join(lines)
        if model_number is not None:
            return f"MODEL{model_number:8}\n{body}\nENDMDL\n"
        return body

    pdb_out = ""
    for index in range(payload["atom_positions"].shape[0]):
        pdb_out += to_pdb_str(jax.tree_util.tree_map(lambda x: x[index], payload), index + 1)
    pdb_out += "END\n"
    with open(filename, "w", encoding="utf-8") as handle:
        handle.write(pdb_out)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir")
    parser.add_argument("output_name")
    parser.add_argument("--params", required=True)
    parser.add_argument("--num-recycles", type=int, default=3)
    parser.add_argument("--model-count", type=int, default=1, choices=range(1, 6))
    parser.add_argument("--use-binder-template", action="store_true", default=False)
    parser.add_argument("--use-interface-template", action="store_true", default=False)
    parser.add_argument(
        "--target-template-only",
        action="store_true",
        default=False,
        help="Use the staged target-only PDB as template input and read the binder sequence from its sidecar file.",
    )
    parser.add_argument(
        "--single-chain-template-only",
        action="store_true",
        default=False,
        help="Fold one staged chain with fixbb/template inputs; used for monomer capacity tests.",
    )
    parser.add_argument(
        "--target-fragment-template-only",
        action="store_true",
        default=False,
        help="Fold a target-only fragmented template while preserving staged fragment chains.",
    )
    parser.add_argument("--multimer", action="store_true", default=False)
    parser.add_argument("--use-msa", action="store_true", default=False)
    parser.add_argument(
        "--binder-multimer",
        action="store_true",
        default=False,
        help="Use AF2-multimer-v3 parameters for the independent binder-alone prediction.",
    )
    parser.add_argument("--designed_chains", default="A")
    parser.add_argument("--cyclic", action="store_true", default=False)
    parser.add_argument("--hotspot", default=None)
    options = parser.parse_args()

    if options.use_interface_template:
        assert options.use_binder_template

    target_only_template_mode = options.single_chain_template_only or options.target_fragment_template_only

    if target_only_template_mode:
        options.target_template_only = False

    model_names = [
        f"model_{index}_multimer_v3" if options.multimer else f"model_{index}_ptm"
        for index in range(1, options.model_count + 1)
    ]
    binder_uses_multimer = options.binder_multimer
    binder_model_names = [
        f"model_{index}_multimer_v3" if binder_uses_multimer else f"model_{index}_ptm"
        for index in range(1, options.model_count + 1)
    ]
    if target_only_template_mode:
        model = mk_af_model(
            protocol="fixbb",
            data_dir=options.params,
            use_multimer=options.multimer,
            model_names=model_names,
            use_templates=True,
            use_initial_guess=True,
            pre_callback=_inject_af2_msa if options.use_msa else None,
        )
        binder_model = None
    else:
        model = mk_af_model(
            protocol="binder",
            data_dir=options.params,
            use_multimer=options.multimer,
            model_names=model_names,
            use_initial_guess=not options.target_template_only,
            pre_callback=_inject_af2_msa if options.use_msa else None,
        )
        binder_model = mk_af_model(
            protocol="hallucination",
            data_dir=options.params,
            use_multimer=binder_uses_multimer,
            model_names=binder_model_names,
            use_initial_guess=False,
            pre_callback=_inject_af2_msa if options.use_msa else None,
        )
    paths = sorted(glob.glob(os.path.join(options.input_dir, "*.pdb")))
    total_lengths = {path: get_pdb_total_length(path) for path in paths}
    paths = sorted(paths, key=lambda path: total_lengths[path])
    print(f"Processing {len(paths):,} PDBs")

    os.makedirs(options.output_name, exist_ok=True)
    with open(options.output_name.rstrip("/") + ".jsonl", "wt", encoding="utf-8") as handle:
        for index, path in enumerate(paths, start=1):
            basename = os.path.basename(path).removesuffix(".pdb")
            print(f"Predicting PDB {index:,}/{len(paths):,}: {basename}")
            start_time = time.time()
            designed_chain_list = get_binder_source_chains(path, options.designed_chains)
            if target_only_template_mode:
                designed_chain_list = _chain_ids_from_pdb(path)
            designed_chains = ",".join(designed_chain_list)
            designed_output_chain = (designed_chain_list or ["A"])[0]
            target_chain = get_target_chains(
                path,
                "B" if designed_chains == "A" else "A",
            )
            target_chain_list = [chain.strip() for chain in target_chain.split(",") if chain.strip()]
            target_chain_lengths = _chain_lengths_from_pdb(path, target_chain_list)

            binder_sequence = get_design_sequence(
                path,
                designed_chains,
                target_only=target_only_template_mode,
            )
            if target_only_template_mode:
                if not binder_sequence or "X" in binder_sequence:
                    raise ValueError(f"Could not extract a complete sequence from chain {designed_chains} in {path}")
                model.prep_inputs(
                    pdb_filename=path,
                    chain=designed_chains,
                    rm_template=False,
                    rm_template_seq=False,
                    rm_template_sc=False,
                    rm_template_ic=False,
                )
            elif options.target_template_only:
                if not binder_sequence or "X" in binder_sequence:
                    raise ValueError(f"Could not extract a complete binder sequence from chain {designed_chains} in {path}")
                model.prep_inputs(
                    path,
                    binder_chain=None,
                    binder_len=len(binder_sequence),
                    target_chain=target_chain,
                    rm_target=False,
                    hotspot=options.hotspot if options.hotspot else None,
                )
                model._target_chain_lengths_for_output = target_chain_lengths
            else:
                model.prep_inputs(
                    path,
                    binder_chain=designed_chains,
                    target_chain=target_chain,
                    rm_target=False,
                    rm_binder=not options.use_binder_template,
                    rm_template_ic=not options.use_interface_template,
                    hotspot=options.hotspot if options.hotspot else None,
                )
                model._target_chain_lengths_for_output = target_chain_lengths
            if options.cyclic:
                add_cyclic_offset(model, offset_type=2)
            if target_only_template_mode:
                model.set_seq(seq=binder_sequence)
            elif options.target_template_only:
                model.set_seq(seq=binder_sequence)
            else:
                model.set_seq(mode="wildtype")
            msa_payload = None
            msa_rows_by_chain = {}
            if options.use_msa:
                msa_payload = _read_af2_msa_input(options.input_dir, basename)
                msa_sequences = msa_payload.get("sequences") or {}
                msa_chain_order = (
                    designed_chain_list
                    if target_only_template_mode
                    else [*target_chain_list, *designed_chain_list]
                )
                msa_rows_by_chain = _attach_af2_msa(model, msa_payload, msa_chain_order, msa_sequences)
            suffix = os.path.basename(options.output_name.rstrip("/"))
            out_base = os.path.join(options.output_name, f"{basename}_{suffix}")
            model_rows = []
            model_outputs = []
            model.set_opt(num_recycles=options.num_recycles)
            for model_index in range(options.model_count):
                model.predict(models=[model_index], num_recycles=options.num_recycles, verbose=False)
                if target_only_template_mode:
                    log = model.aux["log"]
                    row = {
                        "binder_plddt": float(log.get("plddt")) * 100.0 if log.get("plddt") is not None else None,
                        "binder_pae": float(log.get("pae")) * 31.0 if log.get("pae") is not None else None,
                        "ptm": float(log.get("ptm")) if log.get("ptm") is not None else None,
                        "target_aligned_binder_rmsd": None,
                        "ipae": None,
                        "iptm": None,
                        "con_loss": None,
                        "i_con_loss": None,
                        "interface_target_residues": "",
                    }
                else:
                    row = {new_key: model.aux["log"].get(old_key) for old_key, new_key in METRICS.items()}
                    row["binder_plddt"] *= 100
                    row["binder_pae"] *= 31.0
                    row["ipae"] *= 31.0
                    ca_pos = model.aux["atom_positions"][:, 1]
                    ca_dist = np.sqrt(
                        np.square(
                            ca_pos[model._target_len :, None]
                            - ca_pos[None, : model._target_len]
                        ).sum(axis=-1)
                        + 1e-8
                    )
                    target_interface_res = model.aux["residue_index"][: model._target_len][
                        ca_dist.min(axis=0) <= 8
                    ]
                    row["interface_target_residues"] = _format_target_residue_labels(model, target_interface_res)
                model_number = model_index + 1
                model_pdb = f"{out_base}_model{model_number}.pdb"
                model_pae = f"{out_base}_model{model_number}_pae.json"
                if target_only_template_mode:
                    model.save_pdb(model_pdb)
                    _save_monomer_pae_json(model, model_pae, row)
                else:
                    save_binder_design_pdb(model, model_pdb, binder_output_chain=designed_output_chain)
                    _save_pae_json(model, model_pae, row)
                model_rows.append(row)
                model_outputs.append((model_pdb, model_pae))

            metrics = {
                "id": basename,
                "model_count": len(model_rows),
                "af2_validation_parameter_family": (
                    "multimer_v3" if options.multimer else "ptm"
                ),
                "complex_af2_parameter_family": (
                    "multimer_v3" if options.multimer else "ptm"
                ),
                "binder_fold_af2_parameter_family": (
                    "multimer_v3" if binder_uses_multimer else "ptm"
                ),
                "complex_prediction_protocol": (
                    "target_fragment_template_fold"
                    if options.target_fragment_template_only
                    else "single_chain_template_fold"
                    if options.single_chain_template_only
                    else "target_binder_complex"
                ),
                "binder_fold_prediction_protocol": (
                    "not_applicable_target_fragment"
                    if options.target_fragment_template_only
                    else "not_applicable_single_chain"
                    if options.single_chain_template_only
                    else "binder_alone"
                ),
                "msa_used": bool(options.use_msa),
                "msa_rows_by_chain": msa_rows_by_chain,
            }
            for model_number, row in enumerate(model_rows, start=1):
                for key, value in row.items():
                    metrics[f"model_{model_number}_{key}"] = value
            for key in METRICS.values():
                values = [float(row[key]) for row in model_rows if row.get(key) is not None]
                metrics[key] = sum(values) / len(values) if values else None
                metrics[f"average_{key}"] = metrics[key]
            metrics["time"] = time.time() - start_time

            best_index = max(
                range(len(model_rows)),
                key=lambda idx: float(model_rows[idx].get("binder_plddt") or 0.0),
            )
            shutil.copy2(model_outputs[best_index][0], f"{out_base}.pdb")
            if os.path.exists(model_outputs[best_index][1]):
                shutil.copy2(model_outputs[best_index][1], f"{out_base}_pae.json")
            metrics["representative_model"] = best_index + 1
            reference_pdb = os.path.join(options.input_dir, "references", f"{basename}.pdb")
            reference_chains = (
                get_target_source_chains(path, designed_chains)
                if target_only_template_mode
                else get_binder_source_chains(path, designed_chains)
            )
            metrics["interface_target_residues"] = model_rows[best_index][
                "interface_target_residues"
            ]

            if target_only_template_mode:
                metrics["monomer_refolding_rmsd"] = (
                    _aligned_ca_rmsd(reference_pdb, model_outputs[best_index][0], reference_chains)
                    if os.path.exists(reference_pdb)
                    else None
                )
                metrics["average_binder_only_plddt"] = None
                metrics["average_binder_only_pae"] = None
                metrics["average_binder_rmsd"] = metrics["monomer_refolding_rmsd"]
            else:
                binder_model.prep_inputs(length=len(binder_sequence))
                binder_model.set_seq(seq=binder_sequence)
                if options.use_msa and msa_payload:
                    metrics["binder_only_msa_rows_by_chain"] = _attach_af2_msa(
                        binder_model,
                        msa_payload,
                        designed_chain_list,
                        msa_payload.get("sequences") or {},
                    )
                binder_model.set_opt(num_recycles=options.num_recycles)
                for model_index in range(options.model_count):
                    binder_model.predict(models=[model_index], num_recycles=options.num_recycles, verbose=False)
                    binder_number = model_index + 1
                    binder_pdb = f"{out_base}_binder_model{binder_number}.pdb"
                    binder_model.save_pdb(binder_pdb)
                    binder_log = binder_model.aux["log"]
                    binder_plddt = binder_log.get("plddt")
                    binder_pae = binder_log.get("pae")
                    metrics[f"model_{binder_number}_binder_only_plddt"] = (
                        float(binder_plddt) * 100.0 if binder_plddt is not None else None
                    )
                    metrics[f"model_{binder_number}_binder_only_pae"] = (
                        float(binder_pae) * 31.0 if binder_pae is not None else None
                    )
                    metrics[f"model_{binder_number}_binder_rmsd"] = (
                        _aligned_ca_rmsd(reference_pdb, binder_pdb, reference_chains)
                        if os.path.exists(reference_pdb)
                        else None
                    )
            if not target_only_template_mode:
                for metric_name in ("binder_only_plddt", "binder_only_pae", "binder_rmsd"):
                    values = [
                        metrics.get(f"model_{model_number}_{metric_name}")
                        for model_number in range(1, options.model_count + 1)
                    ]
                    numeric = [float(value) for value in values if value is not None]
                    metrics[f"average_{metric_name}"] = sum(numeric) / len(numeric) if numeric else None

            json.dump(metrics, handle)
            handle.write("\n")
            handle.flush()


if __name__ == "__main__":
    main()
