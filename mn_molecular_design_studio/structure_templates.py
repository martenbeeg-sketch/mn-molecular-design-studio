"""Job-local template conversion for the installed cofolding adapters."""
from __future__ import annotations

import copy
import json
from pathlib import Path

import yaml
from Bio.PDB import Model, PDBIO, Structure

from mn_molecular_design_studio.design_validation import protein_chains, read_structure, residue_pairs, sequence

# These formats were checked against the pinned native parsers and batch runners.
TEMPLATE_ENGINES = {
    "alphafold3", "openfold3_jax", "af3_torch", "boltz2", "colabfold",
    "protenix_v2", "openfold3_ob0", "rf3_foundry", "boltzgen_fold", "chai1",
}


def validate_uploads(uploads: object, *, supported: bool) -> list[dict]:
    if not uploads:
        return []
    if not supported:
        raise ValueError("The installed adapter for this engine does not accept uploaded structure templates.")
    if not isinstance(uploads, list) or len(uploads) > 16:
        raise ValueError("Upload at most 16 structure templates.")
    result = []
    for row in uploads:
        if not isinstance(row, dict):
            raise ValueError("Each template needs a file name and structure text.")
        name = Path(str(row.get("name") or "").replace("\\", "/")).name
        text = str(row.get("text") or "")
        if Path(name).suffix.lower() not in {".pdb", ".ent", ".cif", ".mmcif"} or not text.strip():
            raise ValueError("Upload a non-empty PDB or mmCIF template.")
        if len(text.encode("utf-8")) > 20 * 1024 * 1024:
            raise ValueError("Structure templates must be smaller than 20 MB.")
        result.append({**row, "name": name, "text": text})
    return result


def stage_templates(run_dir: Path, path: Path, *, item_id: str, profile: dict,
                    proteins: list[dict], uploads: list[dict]) -> dict:
    if not uploads:
        return {}
    root = run_dir / "artifacts/raw/optimized_prediction/templates" / item_id
    root.mkdir(parents=True, exist_ok=True)
    matched = {}
    for index, upload in enumerate(uploads):
        original = root / f"upload_{index}{Path(upload['name']).suffix}"
        original.write_text(upload["text"], encoding="utf-8")
        chains = protein_chains(read_structure(original))
        if not chains:
            raise ValueError(f"Template {upload['name']} has no protein coordinates.")
        for protein in proteins:
            if upload.get("query_chain_id") and upload["query_chain_id"] != protein["id"]:
                continue
            choices = [(chain, residues) for chain, residues in chains.items()
                       if not upload.get("template_chain_id") or chain == upload["template_chain_id"]]
            best = None
            for chain, residues in choices:
                template_seq = sequence(residues)
                pairs = residue_pairs(protein["sequence"], template_seq)
                identity = sum(protein["sequence"][i] == template_seq[j] for i, j in pairs) / max(1, len(pairs))
                if len(pairs) >= 3 and identity >= .5 and (best is None or identity > best[0]):
                    best = (identity, chain, residues, pairs)
            if best is None:
                continue
            _, chain, residues, pairs = best
            structure = Structure.Structure("template")
            model = Model.Model(0); structure.add(model)
            copied = copy.deepcopy(next(c for c in next(read_structure(original).get_models()) if c.id == chain))
            copied.detach_parent(); copied.id = "A"
            protein_residues = set(protein_chains_from_residues(copied))
            for residue in list(copied):
                if residue not in protein_residues:
                    copied.detach_child(residue.id)
            model.add(copied)
            pdb = root / f"chain_{protein['id']}_{index}.pdb"
            writer = PDBIO(); writer.set_structure(structure); writer.save(str(pdb))
            from mn_molecular_design_studio.engine.workflow_engine.workflows.refolding import _pdb_chain_to_protenix_template_json
            json_path = pdb.with_suffix(".json")
            _pdb_chain_to_protenix_template_json(pdb_path=pdb, chain="A", query_sequence=protein["sequence"],
                                                  output_path=json_path, entry_id=f"tpl{index}")
            template = json.loads(json_path.read_text())[0]
            template.update(queryIndices=[i for i, _ in pairs], templateIndices=[j for _, j in pairs])
            # Chai/Gemmi needs explicit polymer entity metadata for get_polymer().
            # Keep the atom labels, residue order and revision date used by AF3.
            template_seq = sequence(residues)
            template["mmcif"] += (
                "\n_entity.id 1\n_entity.type polymer\n#\n"
                "_entity_poly.entity_id 1\n_entity_poly.type 'polypeptide(L)'\n"
                f"_entity_poly.pdbx_seq_one_letter_code {template_seq}\n"
                f"_entity_poly.pdbx_seq_one_letter_code_can {template_seq}\n"
                "_entity_poly.pdbx_strand_id A\n#\n"
            )
            cif = pdb.with_suffix(".cif"); cif.write_text(template["mmcif"])
            matched.setdefault(protein["id"], []).append({"template": template, "cif": cif, "pdb": pdb})
    if not matched:
        raise ValueError("No uploaded template could be aligned to an input protein chain (minimum 50% identity and three residues).")
    native = lambda value: "/work/" + value.relative_to(run_dir).as_posix()
    kind = profile["input_kind"]
    payload = json.loads(path.read_text()) if path.suffix == ".json" else yaml.safe_load(path.read_text()) if kind == "boltz_yaml" else None
    details = {"template_count": sum(len(rows) for rows in matched.values()), "template_protein_chains": list(matched)}
    if kind == "boltz_yaml":
        for chain_id, entries in matched.items():
            for entry in entries:
                payload.setdefault("templates", []).append({"cif": native(entry["cif"]),
                                                           "chain_id": [chain_id], "template_id": ["A"]})
    elif kind == "af3_json":
        for row in payload["sequences"]:
            protein = row.get("protein")
            if protein:
                chain_id = protein["id"][0] if isinstance(protein["id"], list) else protein["id"]
                if chain_id in matched:
                    protein["templates"] = [entry["template"] for entry in matched[chain_id]]
    elif kind == "protenix_json":
        for row in payload[0]["sequences"]:
            protein = row.get("proteinChain")
            if protein and protein["id"][0] in matched:
                template_path = root / f"templates_{protein['id'][0]}.json"
                template_path.write_text(json.dumps([entry["template"] for entry in matched[protein["id"][0]]]))
                protein["templatesPath"] = native(template_path)
    elif kind == "openfold3_ob0_json":
        for query in payload["queries"].values():
            for protein in query["chains"]:
                chain_id = protein["chain_ids"]
                chain_id = chain_id[0] if isinstance(chain_id, list) else chain_id
                if chain_id in matched:
                    protein["template_cif_paths"] = [native(entry["cif"]) for entry in matched[chain_id]]
                    protein["template_cif_chain_ids"] = ["A"] * len(matched[chain_id])
    elif kind == "multimer_fasta":
        # The installed ColabFold batch API accepts one template file per item.
        merged = Structure.Structure("templates"); model = Model.Model(0); merged.add(model)
        for chain_id, entries in matched.items():
            chain = copy.deepcopy(next(next(read_structure(entries[0]["pdb"]).get_models()).get_chains()))
            chain.detach_parent(); chain.id = chain_id; model.add(chain)
        template_path = root / "complex_template.pdb"
        writer = PDBIO(); writer.set_structure(merged); writer.save(str(template_path))
        details["batch_template_path"] = native(template_path)
    elif kind == "chai1_json":
        fasta = path.parent / payload["items"][0]["fasta"]
        lines = fasta.read_text().splitlines()
        fasta.write_text("\n".join(line.replace("|name=", f"|name={item_id}_") if line.startswith(">") else line for line in lines) + "\n")
        rows = []
        protein_by_id = {p["id"]: p for p in proteins}
        for chain_id, entries in matched.items():
            for index, entry in enumerate(entries):
                length = len(sequence(protein_chains(read_structure(entry["pdb"]))["A"]))
                rows.append("\t".join(map(str, [f"{item_id}_{chain_id}", f"local{chain_id}{index}_A", 100,
                                              length, 0, 0, 1, len(protein_by_id[chain_id]["sequence"]),
                                              1, length, 0, length * 2, native(entry["cif"])])))
        hits = root / "template_hits.m8"; hits.write_text("\n".join(rows) + "\n")
        details["chai_template_hits"] = native(hits)
    elif kind == "rf3_json":
        for component in payload["components"]:
            chain_id = component.get("chain_id")
            if chain_id not in matched:
                continue
            entry = matched[chain_id][0]
            if sequence(protein_chains(read_structure(entry["pdb"]))["A"]) != component["seq"]:
                raise ValueError("RosettaFold3 structure conditioning currently requires an exact sequence match to the template.")
            structure = read_structure(entry["pdb"])
            source_chain = next(next(structure.get_models()).get_chains()); source_chain.id = chain_id
            cif = root / f"rf3_chain_{chain_id}.cif"
            from Bio.PDB import MMCIFIO
            writer = MMCIFIO(); writer.set_structure(structure); writer.save(str(cif))
            component.clear(); component["path"] = native(cif)
            payload.setdefault("template_selection", []).append(chain_id)
    else:
        raise ValueError("No structure-template serializer is registered for this engine.")
    if payload is not None:
        path.write_text(yaml.safe_dump(payload, sort_keys=False) if kind == "boltz_yaml" else json.dumps(payload, indent=2) + "\n")
    return details


def protein_chains_from_residues(chain):
    from Bio.SeqUtils import seq1
    return [r for r in chain if "CA" in r and seq1(r.resname, undef_code="X") != "X"]


def stage_chai_runtime_patch(run_dir: Path) -> Path:
    """Bridge local CIF uploads into the pinned Chai batch adapter in its container."""
    destination = run_dir / "artifacts/raw/optimized_prediction/chai_template_patch.py"
    destination.write_text('''from pathlib import Path
import ast

adapter = Path('/opt/mn_batch_api/adapters.py')
source = adapter.read_text()
anchor = 'allowed_options = {"config", "msa_depth", "num_trunk_recycles", "num_trunk_samples", "low_memory"}'
if source.count(anchor) != 1:
    raise RuntimeError('Pinned Chai batch adapter template option anchor changed')
source = source.replace(anchor, anchor[:-1] + ', "template_hits_path"}')
anchor = '    if bool(request.options.get("low_memory", True)):'
if source.count(anchor) != 1:
    raise RuntimeError('Pinned Chai batch adapter command anchor changed')
source = source.replace(anchor, '    if request.options.get("template_hits_path"):\\n        command.extend(("--template-hits-path", str(request.options["template_hits_path"])))\\n' + anchor)
adapter.write_text(source)

import importlib.util
spec = importlib.util.find_spec('chai_lab.data.parsing.templates.m8')
path = Path(spec.origin)
source = path.read_text()
tree = ast.parse(source)
matches = [n for n in ast.walk(tree) if isinstance(n, ast.Assign) and isinstance(n.value, ast.Call)
           and isinstance(n.value.func, ast.Attribute) and n.value.func.attr == 'download_cif_file'
           and any(isinstance(t, ast.Name) and t.id == 'cif_file' for t in n.targets)]
if len(matches) != 1:
    raise RuntimeError('Pinned Chai local-template parser anchor changed')
node = matches[0]
lines = source.splitlines(True)
indent = ' ' * node.col_offset
original = ''.join('    ' + line for line in lines[node.lineno-1:node.end_lineno])
replacement = (indent + 'if isinstance(row.comment, str) and row.comment.startswith("/work/"):\\n'
               + indent + '    cif_file = Path(row.comment)\\n'
               + indent + '    if not cif_file.is_file():\\n'
               + indent + '        raise ValueError("Uploaded Chai template is missing")\\n'
               + indent + 'else:\\n' + original)
lines[node.lineno-1:node.end_lineno] = [replacement]
path.write_text(''.join(lines))
print('Enabled job-local Chai templates; model settings are unchanged', flush=True)
''')
    return destination
