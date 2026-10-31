"""Input validation and condition-spec construction for RFDiffusion4-Proteina."""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

import yaml

from mn_molecular_design_studio.design_validation import protein_chains, read_structure, sequence


UPSTREAM_VERSION = "RFD4-Proteina 0.0.3"
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
    "score_structures": "Score Structures",
    "unconditional": "Unconditional Generation",
}

TASK_HELP = {
    "protein_binder": "Co-design a binder backbone and sequence around a selected protein target and optional hotspot residues.",
    "nanobody": 'Design CDRs on a built-in VHH framework or a supplied nanobody–target complex.',
    "antibody": 'Design selected CDR residues while retaining the supplied antibody framework and target.',
    "tcr": 'Design selected TCR CDR residues in a supplied target–TCR complex.',
    "motif_scaffolding": 'Build a protein scaffold around selected motif residues and atoms.',
    "enzyme": 'Scaffold catalytic residues with optional ligand chemistry and pose constraints.',
    "small_molecule_binder": 'Design a protein around a ligand, with optional protein context.',
    "dna_binder": 'Design a protein binder conditioned on DNA sequence or structure.',
    "rna_design": 'Generate RNA sequence and coordinates conditioned on dot-bracket secondary structure.',
    "existing_structure": 'Diversify a structure, redesign its sequence, predict coordinates, or pack side chains.',
    "custom_conditions": "Provide a native RFDiffusion4-Proteina YAML or JSON condition specification.",
    "score_structures": "Run refolding, structural comparison, and task-specific pass criteria on existing structures.",
    "unconditional": "Generate protein structures and sequences from scratch in the selected length range.",
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
        if suffix not in {".cif", ".mmcif", ".pdb", ".ent"}:
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


def parse_idr_fasta(text: str) -> str:
    """Return one standard-amino-acid IDR sequence from FASTA or plain text."""
    value = str(text or "").strip()
    if not value:
        raise ValueError("Paste one IDR sequence in FASTA format.")
    if any(line.lstrip().startswith(">") for line in value.splitlines()):
        from Bio import SeqIO
        from io import StringIO

        records = list(SeqIO.parse(StringIO(value), "fasta"))
        if len(records) != 1:
            raise ValueError("Provide exactly one FASTA record for the IDR target.")
        sequence_value = str(records[0].seq)
    else:
        sequence_value = value
    normalized = re.sub(r"\s+", "", sequence_value).upper()
    if not 4 <= len(normalized) <= 1024 or not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", normalized):
        raise ValueError("The IDR target must contain 4–1024 standard amino acids (no gaps or ambiguous residues).")
    return normalized


def _nerf_atom(a, b, c, bond_length: float, bond_angle: float, dihedral: float):
    """Place an atom from three predecessors using ideal peptide geometry."""
    import numpy as np

    theta, phi = np.deg2rad([bond_angle, dihedral])
    bc = c - b
    bc /= np.linalg.norm(bc)
    normal = np.cross(b - a, bc)
    normal /= np.linalg.norm(normal)
    tangent = np.cross(normal, bc)
    return c + bond_length * (
        -np.cos(theta) * bc
        + np.sin(theta) * (np.cos(phi) * tangent + np.sin(phi) * normal)
    )


def write_idr_sequence_carrier(sequence_value: str, destination: Path) -> Path:
    """Write an ideal extended backbone CIF used only to carry sequence identity.

    Its coordinates are deliberately left without C_CRD/C_CRDGRP conditions in
    the IDR condition spec. They are not an IDR structure prediction or a target
    conformation supplied to the generator.
    """
    import numpy as np
    from Bio.PDB import Atom, Chain, MMCIFIO, Model, Residue, Structure
    from Bio.SeqUtils import seq3

    sequence_value = parse_idr_fasta(sequence_value)
    nitrogen = [np.array([0.0, 0.0, 0.0])]
    alpha_carbon = [np.array([1.458, 0.0, 0.0])]
    angle = np.deg2rad(68.8)
    carbon = [alpha_carbon[0] + 1.525 * np.array([np.cos(angle), np.sin(angle), 0.0])]
    for index in range(len(sequence_value) - 1):
        nitrogen.append(_nerf_atom(nitrogen[index], alpha_carbon[index], carbon[index], 1.329, 116.2, 130.0))
        alpha_carbon.append(_nerf_atom(alpha_carbon[index], carbon[index], nitrogen[index + 1], 1.458, 121.7, 180.0))
        phi = -75.0 if sequence_value[index + 1] == "P" else -120.0
        carbon.append(_nerf_atom(carbon[index], nitrogen[index + 1], alpha_carbon[index + 1], 1.525, 111.2, phi))

    structure = Structure.Structure("idr_sequence_carrier")
    model = Model.Model(0)
    chain = Chain.Chain("A")
    structure.add(model)
    model.add(chain)
    serial = 1
    for index, letter in enumerate(sequence_value):
        residue = Residue.Residue((" ", index + 1, " "), seq3(letter).upper(), " ")
        chain.add(residue)
        oxygen = _nerf_atom(nitrogen[index], alpha_carbon[index], carbon[index], 1.231, 120.8, 310.0)
        coords = {"N": nitrogen[index], "CA": alpha_carbon[index], "C": carbon[index], "O": oxygen}
        for atom_name, coordinate in coords.items():
            element = atom_name[0]
            residue.add(Atom.Atom(atom_name, coordinate, 0.0, 1.0, " ", atom_name, serial, element=element))
            serial += 1

    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    writer = MMCIFIO()
    writer.set_structure(structure)
    writer.save(str(destination))
    return destination


def store_idr_sequence_carrier(sequence_value: str, input_root: Path) -> Path:
    """Cache the synthetic sequence carrier under the user's managed inputs."""
    normalized = parse_idr_fasta(sequence_value)
    digest = hashlib.sha256(normalized.encode("ascii")).hexdigest()[:16]
    destination = Path(input_root) / "rfd4-proteina" / "idr-sequences" / digest / "idr_sequence_carrier.cif"
    if not destination.is_file():
        write_idr_sequence_carrier(normalized, destination)
    return destination


def idr_binder_condition_spec(
    sequence_value: str,
    *,
    target_filename: str = "target.cif",
    binder_length_min: int = 70,
    binder_length_max: int = 110,
    hotspots: str = "",
) -> str:
    """Build the experimental sequence-fixed, coordinate-free IDR binder spec."""
    sequence_value = parse_idr_fasta(sequence_value)
    if type(binder_length_min) is not int or type(binder_length_max) is not int or not 20 <= binder_length_min <= binder_length_max <= 1024:
        raise ValueError("Binder length must be an ordered range from 20 to 1024 residues.")
    conditions: dict[str, Any] = {"C_SEQ": {True: [{"select": "sel('A')"}]}}
    selected_positions: list[int] = []
    for token in (part.strip() for part in re.split(r"[,;\s]+", str(hotspots or ""))):
        if not token:
            continue
        match = re.fullmatch(r"(?:([A-Za-z]):?)?(\d+)(?:-(?:([A-Za-z]):?)?(\d+))?", token)
        if not match:
            raise ValueError(f"Could not parse hotspot {token!r}; use sequence positions such as A5 or A8-12.")
        if any(chain and chain.upper() != "A" for chain in (match.group(1), match.group(3))):
            raise ValueError("Sequence-only IDR hotspot selectors use the synthetic target chain A.")
        start, end = int(match.group(2)), int(match.group(4) or match.group(2))
        if start < 1 or end < start or end > len(sequence_value):
            raise ValueError(f"Hotspot {token!r} falls outside target sequence positions 1–{len(sequence_value)}.")
        selected_positions.extend(range(start, end + 1))
    if selected_positions:
        positions = sorted(set(selected_positions))
        selectors = []
        index = 0
        while index < len(positions):
            start = end = positions[index]
            index += 1
            while index < len(positions) and positions[index] == end + 1:
                end = positions[index]
                index += 1
            residue_selector = f"[{start}-{end}]" if start != end else str(start)
            selectors.append({"select": f"sel('A/*/{residue_selector}')"})
        conditions["C_HOT"] = {True: selectors}

    spec = {
        "idr_peptide_binder": {
            "input": Path(target_filename).name,
            "contig": f"{binder_length_min}-{binder_length_max}, /0, A1-{len(sequence_value)}",
            "conditions": conditions,
        }
    }
    return _as_safe_yaml(spec)


def validate_campaign_inputs(parameters: dict[str, Any]) -> int:
    """Validate one RFD4 protein-binder input for attempt-based campaigns.

    Guided and condition-spec inputs are two native RFD4 input forms. A campaign
    attempt can carry either form, but must describe exactly one binder design
    against one staged target so candidate and retry accounting stays one-to-one.
    """
    task = str(parameters.get("task") or "protein_binder")
    mode = str(parameters.get("input_mode") or "guided")
    if task != "protein_binder" or mode not in {"guided", "guided_binder", "condition_spec"}:
        raise ValueError("RFD4 campaigns require one protein-binder design using a guided or condition-spec input.")
    target_value = parameters.get("target_structure")
    if not target_value:
        raise ValueError("RFD4 campaigns require one target structure.")
    target_path = Path(str(target_value)).expanduser()
    target_name = target_path.name
    if not target_name:
        raise ValueError("RFD4 campaigns require a named target structure file.")
    if not parameters.get("target_chains"):
        raise ValueError("RFD4 campaigns require explicit target-chain identities for candidate verification.")
    if parameters.get("avoid_structures") or parameters.get("also_bind_structures") or parameters.get("protein_sequences"):
        raise ValueError("RFD4 campaigns currently support one target without secondary inputs.")
    campaign = parameters.get("distributed_campaign") or {}
    uses_native_acceptance = isinstance(campaign, dict) and any(
        isinstance(stage, dict) and stage.get("workflow") == "analysis.validation_acceptance"
        for stage in campaign.get("validation_stages", [])
    )
    if uses_native_acceptance and not (parameters.get("options") or {}).get("refold", True):
        raise ValueError(
            "Native-evidence Validation & Acceptance requires RFD4 options.refold=true so its "
            "oracle measurements are available; no separate co-folding job will be substituted."
        )

    input_files = [Path(str(value)).expanduser() for value in parameters.get("input_files") or []]
    if mode in {"guided", "guided_binder"}:
        if input_files:
            raise ValueError("Guided RFD4 campaigns accept the target through target_structure only.")
        return 1

    # Condition-spec RFD4 uses files named in YAML. For campaign attempts, all
    # references must resolve to the same target carried by target_structure.
    if any(path.name != target_name for path in input_files):
        raise ValueError("RFD4 condition-spec campaigns require one target structure and no secondary input files.")
    try:
        spec = yaml.safe_load(str(parameters.get("condition_spec") or ""))
    except yaml.YAMLError as exc:
        raise ValueError("RFD4 campaign condition specification is invalid YAML/JSON.") from exc
    if not isinstance(spec, dict) or not spec:
        raise ValueError("RFD4 condition-spec campaigns require a nonempty condition specification.")
    design_names = [name for name in spec if name != "global_args"]
    if len(design_names) != 1:
        raise ValueError("Each RFD4 campaign attempt must contain exactly one named condition-spec design.")
    global_args = spec.get("global_args") or {}
    design = spec[design_names[0]]
    if not isinstance(global_args, dict) or not isinstance(design, dict):
        raise ValueError("RFD4 campaign condition-spec entries must be mappings.")
    input_value = design.get("input", global_args.get("input"))
    references = input_value if isinstance(input_value, list) else [input_value]
    if not references or any(not isinstance(value, str) or Path(value).name != target_name for value in references):
        raise ValueError("RFD4 campaign condition specifications must reference only the selected target structure.")
    return 1


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


HF_REVISION = 'a06e2509c7033b010fc43614e9fd840303588502'
OPERATIONS = {'diversify':'Diversify a Structure','inverse_folding':'Inverse Folding',
              'structure_prediction':'Structure Prediction','sidechain_packing':'Side Chain Packing'}
LIGAND_INTERACTIONS = {
    'Hydrogen-bond donor':('C_HB',1),'Hydrogen-bond acceptor':('C_HB',2),
    'Hydrogen-bond donor and acceptor':('C_HB',3),'Pi-stacking ring':('C_PISTACK',True),
    'Cation-pi ring':('C_PICATION',1),'Cation-pi cation':('C_PICATION',2),
    'Salt bridge, positive':('C_SALT',1),'Salt bridge, negative':('C_SALT',2),
    'Hydrophobic contact':('C_HPHOB',True),'Metal ion to coordinate':('C_MCOORD',1),
}
FRAMEWORKS = ['Vobarilizumab','Caplacizumab','Ozoralizumab','h-NbBCII10','Your own nanobody']
SCORE_DESIGN_TYPES = {
    'unconditional', 'protein_binder', 'nanobody', 'antibody', 'tcr',
    'motif_scaffolding', 'enzyme', 'small_molecule_binder', 'dna_binder',
}


def validate_options(p):
    o = p.get('options') or {}
    if not isinstance(o,dict): raise ValueError('RFD4 options must be an object.')
    for key, choices, default in (
        ('sampling_effort',['Fast','Balanced','Thorough'],'Balanced'),
        ('reward_search',['Off','Interface contacts','Low clashes','Best by refolding'],'Off'),
        ('operation',list(OPERATIONS),'diversify'),
        ('nanobody_framework',FRAMEWORKS,'Vobarilizumab'),
        ('motif_atoms',['all','backbone','sidechain_beyond_cb','tip'],'all'),
        ('hotspot_atoms',['whole','sidechain','tip','single'],'whole'),
        ('ligand_pose',['generated','fixed'],'generated'),('dna_conditioning',['sequence','structure'],'sequence'),
        ('framework_atoms',['backbone','all'],'backbone'),('noise_mode',['fraction','levels'],'fraction'),
        ('ligand_burial',['model','buried','intermediate','exposed'],'model'),
        ('score_design_type',sorted(SCORE_DESIGN_TYPES),'protein_binder'),
    ):
        if o.get(key,default) not in choices: raise ValueError(f'Invalid {key}.')
    for key,default,low,high in (
        ('backbone_noise',1.,.01,1.5),('latent_noise',1.,.01,5.),('backbone_tsr_k',16.,1.,1000.),
        ('backbone_tsr_sigma',3.,.01,1000.),('latent_tsr_k',4.,1.,1000.),('latent_tsr_sigma',12.,.01,1000.),
        ('noise_fraction',.3,.05,1.),('interface_error_limit',None,.1,30.),('binder_rmsd_limit',None,.1,30.),
        ('refold_backbone_rmsd_limit',None,.1,30.),('motif_rmsd_limit',None,.1,30.),
        ('backbone_level',.8,.01,.99),('sequence_level',.8,.01,.99),('redesign_shell',None,1.,30.),
    ):
        value=o.get(key,default)
        if value is not None and (isinstance(value,bool) or not isinstance(value,(float,int)) or not low<=value<=high):
            raise ValueError(f'{key} must be between {low} and {high}.')
    for key in ('repack_target','unindex_motif','deterministic','refold','filter_conditions','apo','energy','developability','cluster','save_trajectory','dock_independently','keep_constant_domains','hold_kept_in_place','scaffold_total_explicit','native_phase_timing'):
        if key in o and type(o[key]) is not bool: raise ValueError(f'{key} must be true or false.')
    for key,maximum in [('cdr_lengths',40),('scaffold_lengths',1024)]:
        rows=o.get(key) or []
        if not isinstance(rows,list) or len(rows)>64:raise ValueError(f'{key} must contain at most 64 length ranges.')
        for row in rows:
            if not isinstance(row,dict) or any(type(row.get(k)) is not int for k in ('lower_bound','upper_bound')) or not 1<=row['lower_bound']<=row['upper_bound']<=maximum:
                raise ValueError(f'{key}: use ordered integer bounds from 1 to {maximum}.')
    interactions=o.get('ligand_interactions') or []
    if not isinstance(interactions,list) or len(interactions)>64:raise ValueError('Use at most 64 ligand interactions.')
    for row in interactions:
        if not isinstance(row,dict) or row.get('interaction') not in LIGAND_INTERACTIONS or not isinstance(row.get('ligand_atoms',''),str):
            raise ValueError('Choose a supported ligand interaction and atom names.')
    if o.get('reward_search','Off')!='Off' and p.get('num_designs',8)>16:
        raise ValueError('Reward Search supports at most 16 designs per request.')
    if o.get('cath') and not re.fullmatch(r'\d+\.(?:\d+|x)\.(?:\d+|x)',o['cath']):
        raise ValueError('Use a CATH class.architecture.topology code such as 2.60.40.')
    for aa in o.get('exclude_amino_acids') or []:
        if aa not in 'ACDEFGHIKLMNPQRSTVWY': raise ValueError('Select canonical amino acids to exclude.')
    for aa, probability in (o.get('favor_amino_acids') or {}).items():
        if aa not in 'ACDEFGHIKLMNPQRSTVWY' or isinstance(probability,bool) or not isinstance(probability,(int,float)) or not 0<=probability<=1:
            raise ValueError('Amino acid preferences must be probabilities between 0 and 1.')
    if set(o.get('exclude_amino_acids') or []) & set(o.get('favor_amino_acids') or {}):
        raise ValueError('An amino acid cannot be both excluded and favored.')
    if o.get('folding_engine'):
        from mn_molecular_design_studio.engine.workflow_engine.workflows.rfd4_oracle import compatible_engines
        molecule_task=o.get('score_design_type','protein_binder') if p.get('task')=='score_structures' else p.get('task','protein_binder')
        if o['folding_engine'] not in compatible_engines(molecule_task):
            raise ValueError('Choose a folding engine compatible with this task’s molecules.')
    if o.get('folding_recycles') is not None and (type(o['folding_recycles']) is not int or not 1<=o['folding_recycles']<=100):
        raise ValueError('Folding recycles must be an integer from 1 to 100.')
    if p.get('task')=='rna_design':
        text=re.sub(r'\s+','',o.get('rna_secondary_structure',''))
        if not 8<=len(text)<=1024 or not re.fullmatch(r'[.?(\)\[\]{}<>A-Za-z]+',text):
            raise ValueError('RNA dot-bracket input must contain 8–1024 valid secondary-structure characters.')
        # The native parser performs the full pseudoknot and bracket validation.
        o['rna_secondary_structure']=text
    return o


def _translate_residues(text, chain_map, residue_map):
    translated=[]
    for token in re.split(r'[,;\s]+',str(text or '').strip()):
        if not token: continue
        m=re.fullmatch(r'([A-Za-z0-9_]+):(-?\d+)(?:-(-?\d+))?',token) or re.fullmatch(r'([A-Za-z_]+)(-?\d+)(?:-(-?\d+))?',token)
        if not m: raise ValueError(f'Invalid residue selection {token!r}; use A25 or A25-35.')
        chain,start,end=m[1],int(m[2]),int(m[3] or m[2])
        if end<start or end-start>10000: raise ValueError('Residue ranges must be ascending and finite.')
        for number in range(start,end+1):
            label=residue_map.get((chain,number))
            if label is None: raise ValueError(f'Residue {chain}{number} is absent or has an ambiguous insertion code.')
            translated.append(f'{chain_map[chain]}:{label}')
    return ','.join(translated)


def _map_structure_parameters(p,target):
    import copy
    p=copy.deepcopy(p)
    try: chain_map,residue_map=_cif_label_mapping(target)
    except ValueError:
        if p.get('target_chains') or p.get('hotspots') or p['task'] not in {'small_molecule_binder','dna_binder'}: raise
        chain_map,residue_map={},{}
    for key in ('target_chains',):
        p[key]=[chain_map.get(c,c) for c in p.get(key) or []]
    for key in ('design_chains','fixed_chains'):
        p['options'][key]=[chain_map.get(c,c) for c in p['options'].get(key) or []]
    p['hotspots']=_translate_residues(p.get('hotspots',''),chain_map,residue_map)
    for key in ('cdr_residues','motif_residues','keep_residues','helix_residues','strand_residues','loop_residues'):
        p['options'][key]=_translate_residues(p['options'].get(key,''),chain_map,residue_map)
    if p['options'].get('active_site_atoms'):
        rows=[]
        for token in p['options']['active_site_atoms'].split():
            residue,sep,atoms=token.rpartition(':')
            if not sep or not atoms:raise ValueError('Use active-site atoms such as A102:CG,OD2.')
            translated=_translate_residues(residue,chain_map,residue_map)
            if ',' in translated:raise ValueError('Specify active-site atoms for one residue at a time.')
            rows.append(translated+':'+atoms)
        p['options']['active_site_atoms']=' '.join(rows)
    if p['options'].get('designed_chain'):
        p['options']['designed_chain']=chain_map.get(p['options']['designed_chain'],p['options']['designed_chain'])
    candidates=p['options'].get('source_candidates') or {}
    original_name=target.name[:-4] if target.name.lower().endswith(('.pdb.cif','.ent.cif')) else target.name
    if original_name in candidates:
        parent=copy.deepcopy(candidates[original_name])
        for key in ('target_chains','binder_chains'):parent[key]=[chain_map.get(c,c) for c in parent.get(key) or []]
        p['_source_candidate']=parent
    return p


def prepare_request(parameters: dict[str, Any], *, run_dir: Path) -> dict[str, Any]:
    """Validate/stage inputs; native chemistry and condition validation runs in the image."""
    import copy
    from Bio.PDB import MMCIFIO
    p=copy.deepcopy(parameters)
    p['options']=validate_options(p)
    task=p.get('task','protein_binder'); mode=p.get('input_mode','guided')
    if task not in TASKS: raise ValueError('Select a supported RFD4 task.')
    if mode not in {'guided','guided_binder','annotated_cif','condition_spec','unconditional'}:
        raise ValueError('Select a supported RFD4 input method.')
    for key,default,low,high in [('num_designs',8,1,500000),('batch_size',1,1,16),('nsteps',60,60,120),('cpu_cores',8,1,256)]:
        if type(p.get(key,default)) is not int or not low<=p.get(key,default)<=high:
            raise ValueError(f'{key} must be an integer between {low} and {high}.')
    seed=p.get('seed')
    if seed is not None and (type(seed) is not int or not 0<=seed<=4294967295):
        raise ValueError('Random seed must be blank or an integer from 0 to 4294967295.')
    low,high=p.get('length_min',60),p.get('length_max',120)
    if type(low) is not int or type(high) is not int or not 20<=low<=high<=1024:
        raise ValueError('Chain length must be an ordered range from 20 to 1024 residues.')
    root=run_dir/'artifacts/rfd4_proteina'; staging=root/'input'; staging.mkdir(parents=True,exist_ok=True)
    staged={}; originals={}; total=0
    def stage(source):
        nonlocal total
        source=Path(source).expanduser().resolve()
        if not source.is_file() or source.suffix.lower() not in {'.pdb','.ent','.cif','.mmcif'}:
            raise ValueError(f'Missing or unsupported structure: {source.name}')
        total+=source.stat().st_size
        if source.stat().st_size>50_000_000 or total>100_000_000:
            raise ValueError('Structure inputs must be under 50 MB each and 100 MB total.')
        if source.name in originals:
            if originals[source.name]!=source: raise ValueError('Structure filenames must be unique: '+source.name)
            return staged[source.name]
        destination=staging/source.name
        if source.suffix.lower() in {'.pdb','.ent'}:
            destination=staging/(source.name+'.cif')
            writer=MMCIFIO(); writer.set_structure(read_structure(source)); writer.save(str(destination))
        elif source!=destination:
            destination.write_bytes(source.read_bytes())
        originals[source.name]=source; staged[source.name]=destination
        return destination
    files=[stage(x) for x in p.get('input_files') or []]
    if len(files)>128: raise ValueError('Upload at most 128 structures.')
    target=stage(p['target_structure']) if p.get('target_structure') else None
    for key in ('avoid_structures','also_bind_structures'):
        p[key]=['/work/'+stage(x).relative_to(run_dir).as_posix() for x in p.get(key) or []]
    if len(p['avoid_structures'])+len(p['also_bind_structures'])>6:
        raise ValueError('Use at most six structures across both specificity lists.')
    if task=='score_structures' and (not files or len(files)>100): raise ValueError('Provide 1–100 structures to score.')
    if mode in {'guided','guided_binder'} and task not in {'unconditional','rna_design','score_structures'} and target is None and not (task=='existing_structure' and files):
        raise ValueError('Load an input structure for this task.')
    if task in {'enzyme','small_molecule_binder','dna_binder'} and target is not None and originals.get(Path(parameters['target_structure']).name,Path()).suffix.lower() in {'.pdb','.ent'}:
        raise ValueError('Use mmCIF for ligand, enzyme, and DNA tasks to retain chemical bonds and charges.')
    file_parameters={}
    if mode in {'guided','guided_binder'}:
        if files and task in {'score_structures','existing_structure'}:
            file_parameters={'/work/'+path.relative_to(run_dir).as_posix():_map_structure_parameters(p,path) for path in files}
        if target is not None:p=_map_structure_parameters(p,target)
    if mode=='unconditional' and task!='unconditional': raise ValueError('Unconditional input requires the unconditional task.')
    spec_path=staging/'condition_spec.yaml'
    design_count=1
    if mode=='annotated_cif':
        if not files: raise ValueError('Upload annotated CIF structures.')
        if any(Path(x).suffix.lower() not in {'.cif','.mmcif'} for x in parameters.get('input_files') or []):
            raise ValueError('Annotated inputs must be CIF/mmCIF files.')
        spec={f'design_{i+1}':{'input':'/work/'+path.relative_to(run_dir).as_posix(),'input_as_annotated_cif':True} for i,path in enumerate(files)}
        spec_path.write_text(_as_safe_yaml(spec)); design_count=len(spec)
    elif mode=='condition_spec':
        try: spec=yaml.safe_load(str(p.get('condition_spec') or ''))
        except yaml.YAMLError as exc: raise ValueError('Invalid YAML/JSON condition file.') from exc
        if not isinstance(spec,dict) or not spec: raise ValueError('Provide a nonempty mapping of named designs.')
        global_args=spec.get('global_args') or {}
        if not isinstance(global_args,dict): raise ValueError('global_args must be a mapping.')
        resolved={}
        for name,entry in spec.items():
            if name=='global_args': continue
            if not isinstance(name,str) or not re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}',name): raise ValueError('Use simple, unique design names without paths.')
            if not isinstance(entry,dict): raise ValueError('Each design must be a mapping.')
            entry={**copy.deepcopy(global_args),**entry}
            if set(entry)-{'input','input_as_annotated_cif','contig','conditions'}:
                raise ValueError(f'{name}: unsupported condition-spec fields; use input, contig, conditions, input_as_annotated_cif.')
            refs=entry.get('input'); refs=refs if isinstance(refs,list) else [refs]
            if not refs or any(not isinstance(ref,str) or Path(ref).name not in staged for ref in refs):
                raise ValueError(f'{name}: reference uploaded input filenames in the condition file.')
            rewritten=['/work/'+staged[Path(ref).name].relative_to(run_dir).as_posix() for ref in refs]
            entry['input']=rewritten if isinstance(entry.get('input'),list) else rewritten[0]
            resolved[name]=entry
        design_count=len(resolved)
        if not 1<=design_count<=128: raise ValueError('Provide 1–128 named designs.')
        spec_path.write_text(_as_safe_yaml(resolved))
    if design_count*p.get('num_designs',8)>500000: raise ValueError('At most 500,000 structures may be requested in total.')
    request={'parameters':p,'root':'/work/'+root.relative_to(run_dir).as_posix(),
             'target_structure':'/work/'+target.relative_to(run_dir).as_posix() if target else None,
             'input_files':['/work/'+x.relative_to(run_dir).as_posix() for x in files],
             'spec_path':'/work/'+spec_path.relative_to(run_dir).as_posix(),'file_parameters':file_parameters}
    request_path=root/'request.json'; request_path.write_text(json.dumps(request,indent=2,default=str))
    return {'root':root,'input_root':staging,'spec_path':spec_path,'request_path':request_path,
            'staged_files':list(staged.values()),'design_count':len(files) if task in {'existing_structure','score_structures'} and files else design_count,'parameters':p,
            'target_structure':str(target) if target else None,'target_chains':p.get('target_chains') or []}
