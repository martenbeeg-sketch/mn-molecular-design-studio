"""Container entry point. Only installed upstream APIs run inside the image."""
from __future__ import annotations

import json
import os
import re
import sys
from pathlib import Path


def checkpoint_outputs(request, outputs):
    """Keep finished sequence identities durable between native inputs/samples."""
    if not request.get('output'):
        return  # Direct library callers can request in-memory results.
    path = Path(request['output'])
    temporary = path.with_suffix('.partial.tmp')
    temporary.write_text(json.dumps(outputs, indent=2) + '\n')
    temporary.replace(path)


def _native_sequences(item):
    result = {}
    for row in item["mapping"]:
        result[row["engine_chain"]] = result.get(row["engine_chain"], "") + (row["aa"] if row["aa"] != "X" else "G")
    return result


def esm_if1(request):
    import esm
    import torch
    from esm.inverse_folding.multichain_util import load_complex_coords, _concatenate_coords, score_sequence_in_complex
    p = request["parameters"]
    torch.manual_seed(p["seed"])
    model, alphabet = esm.pretrained.load_model_and_alphabet_local(request["checkpoint"])
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.eval().to(device)
    outputs = []
    with torch.inference_mode():
        for item in request["items"]:
            native = _native_sequences(item)
            coords, _ = load_complex_coords(item["path"], list(native))
            mutable = set(item["mutable"])
            mapping = item["mapping"]
            for _ in range(p["num_sequences"]):
                sequences = dict(native)
                likelihoods = []
                for chain in item["design_chains"]:
                    engine_chain = item["chain_map"][chain]
                    rows = [r for r in mapping if r["engine_chain"] == engine_chain]
                    context = coords if p["multichain_backbone"] else {engine_chain: coords[engine_chain]}
                    joined = _concatenate_coords(context, engine_chain)
                    partial = ["<pad>"] * len(joined)
                    for i, row in enumerate(rows):
                        partial[i] = "<mask>" if row["label"] in mutable else native[engine_chain][i]
                    sampled = model.sample(joined, partial_seq=partial, temperature=p["temperature"], device=device)
                    sequence = sampled[:len(rows)]
                    sequences[engine_chain] = sequence
                    likelihoods.append(score_sequence_in_complex(model, alphabet, context, engine_chain, sequence)[0])
                outputs.append(dict(key=item["key"], sequences=sequences,
                    metrics=dict(esm_if1_log_likelihood=float(sum(likelihoods) / len(likelihoods)))))
                checkpoint_outputs(request, outputs)
    return outputs


def _clean_chain_map(path, expected):
    from Bio.PDB import MMCIFParser
    from Bio.SeqUtils import seq1
    model = next(MMCIFParser(QUIET=True).get_structure("clean", path).get_models())
    chains = list(model)
    if len(chains) != len(expected):
        raise ValueError("Caliby cleaning changed the protein chain inventory.")
    mapping = {}
    for chain, (original, sequence) in zip(chains, expected.items()):
        observed = "".join(seq1(r.resname, undef_code="X") for r in chain if "CA" in r)
        if observed != sequence:
            raise ValueError("Caliby cleaning changed residue identities or order.")
        mapping[original] = chain.id
    return mapping


def caliby_design(request):
    p = request["parameters"]
    if p["optimization_mode"] != "off":
        import caliby_opt
        report = caliby_opt.enable(p["optimization_mode"], variant="ensemble32" if p["ensemble_mode"] == "generate" else "single",
                                   strict=True, ckpt=p["model_type"])
        print(json.dumps(report, default=str), flush=True)
    import caliby
    import lightning as L
    L.seed_everything(p["seed"])
    model = caliby.load_model(p["model_type"])
    outputs = []
    for item in request["items"]:
        directory = Path(item["path"]).parent / "caliby"
        directory.mkdir(exist_ok=True)
        cleaned = caliby.clean_pdbs([item["path"]], out_dir=str(directory / "cleaned"), num_workers=1)[0]
        native = _native_sequences(item)
        chain_map = _clean_chain_map(cleaned, native)
        def position(label):
            match = re.fullmatch(r"(.+?)(\d+)", label)
            return chain_map[match[1]] + match[2]
        constraints = {}
        for key, value in item["constraints"].items():
            constraints[key] = "|".join(",".join(position(entry.split(":", 1)[0]) + (":" + entry.split(":", 1)[1] if ":" in entry else "")
                for entry in group.split(",")) for group in value.split("|"))
        fixed = [position(r["engine_label"]) for r in item["mapping"] if r["label"] not in item["mutable"]]
        fixed += constraints.get("fixed_pos_seq", "").split(",")
        constraints["fixed_pos_seq"] = ",".join(dict.fromkeys(v for v in fixed if v))
        constraint_df = caliby.make_constraints({Path(cleaned).stem: constraints})
        ensemble = None
        if p["ensemble_mode"] == "generate":
            generated = caliby.generate_ensembles([cleaned], out_dir=str(directory / "ensembles"),
                num_samples_per_pdb=p["num_conformers"], batch_size=p["batch_size"], seed=p["seed"])
            conformers = generated[Path(cleaned).stem]
            if len(conformers) != p["num_conformers"]:
                raise ValueError("Protpardelle did not return the requested number of conformers.")
            ensemble = {Path(cleaned).stem: [cleaned, *conformers]}
        elif p["ensemble_mode"] == "provided":
            cleaned_conformers = caliby.clean_pdbs(item["ensemble_paths"], out_dir=str(directory / "conformers"), num_workers=1)
            for path in cleaned_conformers:
                if _clean_chain_map(path, native) != chain_map:
                    raise ValueError("Provided ensemble chain mappings do not match the primary conformer.")
            ensemble = {Path(cleaned).stem: [cleaned, *cleaned_conformers]}
        if p["task"] == "score":
            result = model.score_ensemble(ensemble, num_workers=0) if ensemble else model.score([cleaned], batch_size=p["batch_size"], num_workers=0)
        else:
            kwargs = dict(out_dir=str(directory / "designed"), num_seqs_per_pdb=p["num_sequences"],
                batch_size=p["batch_size"], num_workers=0, temperature=p["temperature"], omit_aas=list(p["omit_aas"]) or None,
                pos_constraint_df=constraint_df, sampling_overrides={"potts_sampling_cfg": {"potts_sweeps": p["num_sweeps"]}})
            result = model.ensemble_sample(ensemble, **kwargs) if ensemble else model.sample([cleaned], **kwargs)
        reverse = {v: k for k, v in chain_map.items()}
        for index, sequence in enumerate(result["seq"]):
            pieces = sequence.split(":")
            if len(pieces) != len(chain_map):
                raise ValueError("Caliby returned an unexpected number of chain sequences.")
            sequences = {reverse[chain]: seq for chain, seq in zip(sorted(reverse), pieces)}
            outputs.append(dict(key=item["key"], sequences=sequences, metrics=dict(caliby_energy=float(result["U"][index]))))
            checkpoint_outputs(request, outputs)
    return outputs


def boltzgen_inverse(request):
    """Use upstream's inverse_fold_only configuration with arbitrary file inputs."""
    import hydra
    import pytorch_lightning as L
    import yaml
    from Bio.PDB import MMCIFParser
    from Bio.SeqUtils import seq1, seq3
    from importlib.resources import files
    from omegaconf import OmegaConf
    import boltzgen
    p = request["parameters"]
    outputs = []
    L.seed_everything(p["seed"])
    cache = Path(request.get("reference_root", "/ref")) / "boltzgen-cache"
    molecule_archives = sorted((cache / "datasets--boltzgen--inference-data/snapshots").glob("*/mols.zip"))
    if not molecule_archives:
        raise ValueError("BoltzGen canonical molecule archive is missing from reference files.")
    for item in request["items"]:
        directory = Path(item["path"]).parent / "boltzgen"
        directory.mkdir(exist_ok=True)
        design = []
        for chain in item["design_chains"]:
            rows = [r for r in item["mapping"] if r["chain"] == chain and r["label"] in item["mutable"]]
            if rows:
                design.append(dict(chain=dict(id=item["chain_map"][chain], res_index=",".join(str(r["position"]) for r in rows))))
        spec = directory / "inverse_fold.yaml"
        spec.write_text(yaml.safe_dump(dict(entities=[dict(file=dict(path=item["path"], design=design))])))
        config = OmegaConf.load(str(files(boltzgen) / "resources/config/inverse_fold_only.yaml"))
        config.output = str(directory / "output")
        config.checkpoint = request["checkpoint"]
        config.data.cfg.yaml_path = [str(spec)]
        config.data.cfg.moldir = str(molecule_archives[-1])
        config.data.cfg.multiplicity = p["num_sequences"]
        config.data.num_workers = 0
        config.trainer.devices = 1
        config.override.inverse_fold_args.sampling_temperature = p["temperature"]
        config.override.inverse_fold_args.inverse_fold_restriction = [seq3(v).upper() for v in p["omit_aas"]]
        OmegaConf.save(config, directory / "native_config.yaml")
        task = hydra.utils.instantiate(config)
        task.run(config)
        native = _native_sequences(item)
        expected = set(native)
        for path in sorted((directory / "output").rglob("*.cif")):
            if path.stem.endswith("_native"):
                continue
            model = next(MMCIFParser(QUIET=True).get_structure("inverse", str(path)).get_models())
            observed = {c.id: "".join(seq1(r.resname, undef_code="X") for r in c if "CA" in r) for c in model}
            sequences = {c: seq for c, seq in observed.items() if c in expected}
            if set(sequences) != expected:
                raise ValueError("BoltzGen changed the input protein chain identifiers.")
            outputs.append(dict(key=item["key"], sequences=sequences, metrics={}))
            checkpoint_outputs(request, outputs)
    return outputs


if __name__ == "__main__":
    request = json.loads(Path(sys.argv[1]).read_text())
    run = {"esm_if1": esm_if1, "caliby": caliby_design, "boltzgen_inverse": boltzgen_inverse}[request["engine"]]
    results = run(request)
    Path(request["output"]).write_text(json.dumps(results, indent=2) + "\n")
