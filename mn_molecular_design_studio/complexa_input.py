"""Validated Studio inputs and native Proteina-Complexa configuration."""
from __future__ import annotations

import copy
import math
import re
from pathlib import Path

from mn_molecular_design_studio.rfdiffusion3_input import inventory
from mn_molecular_design_studio.runtime import reference_root

IMAGE = "mn-proteina-complexa:latest"
CODE_VERSION = "Proteina-Complexa 1.0.0 · 32b71ae"
ASSETS = Path(__file__).with_name("assets") / "proteina_complexa"
TASKS = {"protein_binder":"Protein Binder", "ligand_binder":"Ligand Binder", "motif_ligand":"Motif + Ligand"}
TASK_CONFIGS = {"protein_binder":"search_binder_local_pipeline", "ligand_binder":"search_ligand_binder_local_pipeline", "motif_ligand":"search_ame_local_pipeline"}
CHECKPOINTS = {"protein_binder":("complexa.ckpt","complexa_ae.ckpt"), "ligand_binder":("complexa_ligand.ckpt","complexa_ligand_ae.ckpt"), "motif_ligand":("complexa_ame.ckpt","complexa_ame_ae.ckpt")}
DEFAULTS = {"protein_binder":(70,150,4,"best-of-n"), "ligand_binder":(80,140,2,"best-of-n"), "motif_ligand":(120,220,4,"single-pass")}
EXAMPLES = {"protein_binder": "PD-L1.pdb", "ligand_binder": "5SDV_ligand.pdb", "motif_ligand": "M0024_1nzy_v3.pdb"}
MOTIF_EXAMPLE = "B64: [O, C]; B86: [CB, CA, N, C]; B90: [CE1, ND1, NE2, CG, CD2]; B114: [N, CA]; B137: [NE1, CD1, CE2, CG, CD2, CZ2]; B145: [OD2, CG, CB, OD1]"


def ligand_choices(info):
    return {f"{row['label']}:{name}":row for name,rows in info["ligands"].items() for row in rows}


def parse_motif_atoms(text, info):
    result = {}
    for segment in str(text or "").split(";"):
        if not segment.strip():
            continue
        match = re.fullmatch(r"\s*([A-Za-z][0-9]+)\s*:\s*\[([^\[\]]+)\]\s*",segment)
        if not match:
            raise ValueError("Use motif atoms such as B64: [O, C]; B86: [CB, CA, N, C].")
        label, atom_text = match.groups()
        row = info["residues"].get(label)
        if row is None or row["kind"] != "protein" or row["insertion_code"]:
            raise ValueError(f"Motif residue {label} is not an addressable protein residue in this structure.")
        atoms = [a.strip() for a in atom_text.split(",") if a.strip()]
        if not atoms or len(atoms)!=len(set(atoms)) or label in result:
            raise ValueError("Select each motif residue and atom once.")
        missing = set(atoms)-set(row["atoms"])
        if missing:
            raise ValueError(f"{label} has no atom {', '.join(sorted(missing))}.")
        result[label] = atoms
    if not result:
        raise ValueError("Select at least one motif residue and its fixed atoms.")
    return result


def _integer(value, label, minimum, maximum):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or int(value)!=value or not minimum<=value<=maximum:
        raise ValueError(f"{label} must be a whole number from {minimum} to {maximum}.")
    return int(value)


def prepare_request(parameters):
    p = copy.deepcopy(parameters)
    task = p.get("pipeline_task","protein_binder")
    if task not in TASKS:
        raise ValueError("Choose a supported Proteina-Complexa pipeline.")
    p["pipeline_task"] = task
    path = Path(p.get("target_pdb") or "")
    if not path.is_file():
        raise ValueError("Load a target structure before submitting.")
    p["target_pdb"] = str(path)
    info = inventory(path)
    length = str(p.get("binder_length") or "")
    match = re.fullmatch(r"(\d+)(?:-(\d+))?",length)
    if not match or not 1<=int(match[1])<=int(match[2] or match[1])<=3000:
        raise ValueError("Binder length must be a positive length or increasing range.")
    p["length_range"] = [int(match[1]),int(match[2] or match[1])]
    for name,label,default,lo,hi in [("num_designs","Designs per job",1,1,10000),("n_steps","Sampling steps",400,50,1000),("replicas","Candidates per design",2,1,8),("batch_size","GPU batch size",1,1,16),("cpu_cores","CPU cores",4,1,256)]:
        p[name] = _integer(p.get(name,default),label,lo,hi)
    seed = p.get("seed",0)
    p["seed"] = _integer(0 if seed is None else seed,"Random seed",0,2147483647)
    for name,default in [("guidance_weight",1.),("sampling_noise",.1)]:
        value = p.get(name,default)
        if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or not 0<=value<=1:
            raise ValueError(f"{name.replace('_',' ').capitalize()} must be between 0 and 1.")
        p[name] = float(value)
    for name,default in [("deduplicate_sequences",True),("generator_only",False)]:
        if type(p.get(name,default)) is not bool:
            raise ValueError(f"{name} must be true or false.")
        p[name] = p.get(name,default)
    strategy = p.get("search_strategy","auto")
    strategy = DEFAULTS[task][3] if strategy=="auto" else strategy
    if strategy not in {"best-of-n","single-pass"}:
        raise ValueError("Choose Best of N or Single Pass search.")
    if task=="motif_ligand" and strategy!="single-pass":
        raise ValueError("Motif + Ligand uses Single Pass in this installed pipeline (no search reward is configured).")
    p["search_strategy"] = strategy
    if task=="protein_binder":
        chains = p.get("target_chains") or []
        if not chains or len(chains)!=len(set(chains)) or any(c not in info["chains"] or info["chains"][c]["kind"]!="protein" or not re.fullmatch(r"[A-Za-z]",c) for c in chains):
            raise ValueError("Select valid protein target chains with single-letter IDs.")
        selected = {r["label"] for c in chains for r in info["chains"][c]["residues"] if not r["insertion_code"]}
        p["hotspots"] = ",".join(t.strip() for t in str(p.get("hotspots") or "").split(",") if t.strip())
        if any(h not in selected for h in p["hotspots"].split(",") if h):
            raise ValueError("Hotspots must be residues in the selected target chains.")
        if p.get("ligand") or p.get("motif_atom_spec"):
            raise ValueError("Protein Binder does not take ligand or motif selections.")
    else:
        options = ligand_choices(info)
        chosen = p.get("ligand") or ""
        if chosen in info["ligands"] and len(info["ligands"][chosen])==1:
            chosen = f"{info['ligands'][chosen][0]['label']}:{chosen}"
        if chosen not in options:
            raise ValueError("Select one ligand instance present in the uploaded structure.")
        row = options[chosen]
        if row["insertion_code"] or not re.fullmatch(r"[A-Za-z][0-9]+",row["label"]):
            raise ValueError("The selected ligand needs a letter chain ID and residue number without an insertion code.")
        p["ligand"] = chosen
        p["ligand_residue"] = row
        p["target_chains"] = []
        p["hotspots"] = ""
        p["motif_atoms"] = parse_motif_atoms(p.get("motif_atom_spec"),info) if task=="motif_ligand" else {}
        # JSON job files sort dictionary keys. Keep residue order explicit so
        # prepared numbering and contig atom selections always refer to the
        # same residues after a scheduler round trip.
        p["motif_residue_order"] = list(p["motif_atoms"])
    return p


def native_overrides(p, target_input=""):
    task = p["pipeline_task"]
    group = "motif_target_dict_cfg" if task=="motif_ligand" else "target_dict_cfg"
    model,ae = CHECKPOINTS[task]
    target = {"source":"studio", "target_filename":"target", "target_path":"/work/artifacts/raw/proteina_complexa/input/target.pdb", "binder_length":p["length_range"], "pdb_id":None, "hotspot_residues":[h for h in p["hotspots"].split(",") if h]}
    if task=="protein_binder":
        target["target_input"] = target_input
    else:
        target.update(ligand="L:0",ligand_only=task=="ligand_binder",SMILES=None,use_bonds_from_file=True)
        if task=="motif_ligand":
            order = p.get("motif_residue_order", list(p["motif_atoms"]))
            target["contig_atoms"] = "; ".join(f"B{i}: [{', '.join(p['motif_atoms'][label])}]" for i,label in enumerate(order,1))
    result = {"run_name":p["run_name"], "seed":p["seed"], "ncpus_":p["cpu_cores"], "ckpt_path":"/weights", "ckpt_name":model, "autoencoder_ckpt_path":"/weights/"+ae, "gen_njobs":1, "eval_njobs":1,
        "generation":{"task_name":"MN_APP_TARGET",group:{"MN_APP_TARGET":target},
            "dataloader":{"batch_size":p["batch_size"],"dataset":{"nres":{"nsamples":p["num_designs"]}}},
            "search":{"algorithm":p["search_strategy"],"max_batch_size":p["batch_size"],"best_of_n":{"replicas":p["replicas"]}},
            "args":{"nsteps":p["n_steps"],"guidance_w":p["guidance_weight"]},
            "model":{name:{"simulation_step_params":{"sc_scale_noise":p["sampling_noise"]}} for name in ("bb_ca","local_latents")},
            "filter":{"dedup_sequence":p["deduplicate_sequences"],"delete_non_top_n_samples":False}},
        "metric":{"keep_folding_outputs":True}}
    if p.get("generator_only"):
        # The native generate command can still invoke folding rewards through
        # its search configuration. A standalone engine must disable those too.
        result["generation"]["reward_model"] = None
        result["generation"]["search"].update(algorithm="single-pass", best_of_n={"replicas":1})
    return result


def resource_paths(root=None):
    root = Path(root or reference_root())
    return {"reference":root,"checkpoints":root/"proteina-complexa/ckpts", "rf3":root/"foundry/rf3_foundry_01_24_latest_remapped.ckpt",
            "af2":root/"alphafold_models", "mpnn":root/"protpardelle-1c/model_params/LigandMPNN", "hf_cache":root/"proteina-complexa/hf-cache", "tools":root/"proteina-complexa/tools"}


def readiness(task="protein_binder",root=None, *, generator_only=False):
    paths = resource_paths(root)
    required = [(name,paths["checkpoints"]/name) for name in CHECKPOINTS[task]]
    if generator_only:
        missing = [label for label,path in required if not path.is_file() or not path.stat().st_size]
        return {"ready":not missing,"missing":missing,"paths":{k:str(v) for k,v in paths.items()}}
    required += [("RF3 checkpoint",paths["rf3"]), ("AF2 model weights",paths["af2"]), ("LigandMPNN weights",paths["mpnn"]/"ligandmpnn_v_32_010_25.pt")]
    required += [(name+" executable",paths["tools"]/name) for name in ("foldseek","mmseqs")]
    missing = [label for label,path in required if not path.exists() or path.is_file() and path.stat().st_size==0]
    for model in (["esm2_t33_650M_UR50D","esmfold_v1"] if task!="motif_ligand" else ["esm2_t33_650M_UR50D"]):
        snapshots=paths["hf_cache"]/f"models--facebook--{model}"/"snapshots"
        if not any(p.is_file() and p.stat().st_size for pattern in ("*/model.safetensors","*/pytorch_model.bin") for p in snapshots.glob(pattern)):
            missing.append(model+" cached weights")
    return {"ready":not missing,"missing":missing,"paths":{k:str(v) for k,v in paths.items()}}


def stage_target(prepared, destination):
    """Stage an input copy; retain original ligand bonding information for runtime."""
    import io
    import shutil
    from Bio.PDB import PDBIO, Select
    from mn_molecular_design_studio.design_validation import read_structure
    destination = Path(destination); destination.mkdir(parents=True,exist_ok=True)
    source = Path(prepared["target_pdb"])
    original = destination/("source"+source.suffix.lower())
    shutil.copy2(source,original)
    structure = read_structure(source)
    first_model = next(structure.get_models()).id
    task = prepared["pipeline_task"]
    keep = set(prepared.get("motif_atoms",{}))
    if task!="protein_binder":
        keep.add(prepared["ligand_residue"]["label"])
    class Selection(Select):
        def accept_model(self,model): return model.id==first_model
        def accept_residue(self,residue):
            label=f"{residue.parent.id}{residue.id[1]}"
            return int(residue.parent.id in prepared["target_chains"] and residue.id[0]==" " if task=="protein_binder" else label in keep)
        def accept_atom(self,atom): return int(atom.element not in {"H","D"} and atom.altloc in {" ","A"})
    writer=PDBIO();writer.set_structure(structure);buffer=io.StringIO()
    writer.save(buffer,Selection(),preserve_atom_numbering=True)
    text=buffer.getvalue()
    # Preserve explicit bond records in PDB inputs. CIF bond data remains in the
    # untouched source and is read by AtomWorks in the installed image.
    serials={int(line[6:11]) for line in text.splitlines() if line.startswith(("ATOM  ","HETATM"))}
    connects=[]
    if source.suffix.lower() in {".pdb",".ent"}:
        for line in source.read_text(errors="replace").splitlines():
            if not line.startswith("CONECT"):continue
            try: ids=[int(line[i:i+5]) for i in range(6,len(line),5) if line[i:i+5].strip()]
            except ValueError:continue
            if ids and ids[0] in serials:
                selected=[v for v in ids[1:] if v in serials]
                if selected:connects.append("CONECT"+"".join(f"{v:5d}" for v in [ids[0],*selected]))
    target=destination/"target.pdb"
    target.write_text("\n".join(l for l in text.splitlines() if l!="END")+"\n"+"\n".join(connects)+"\nEND\n")
    return target,original
