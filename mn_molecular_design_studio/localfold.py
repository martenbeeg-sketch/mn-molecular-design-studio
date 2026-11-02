"""Pinned LocalFold model choices and Studio's input contract."""
from __future__ import annotations

import copy
import json
import os
from pathlib import Path
from datetime import datetime, timezone

import yaml

SOURCE_COMMIT = 'e68d628df1d815cccb48bb1118d8763a594bfa23'
IMAGE = os.getenv('MN_MOLECULAR_DESIGN_STUDIO_LOCALFOLD_IMAGE', 'mn-localfold:e68d628-cu13')
MODELS = {
    # Studio follows LocalFold's visible control defaults (three recycles and
    # the shorter LocalFold sampler presets), while retaining model-specific
    # MSA caps and steps where the JAX runner supports them.
    'alphafold3': {'label': 'AF3', 'recycles': 3, 'diffusion_steps': 25, 'msa_depth': 128},
    'openfold3': {'label': 'OpenFold3 preview 2', 'recycles': 3, 'diffusion_steps': 25, 'msa_depth': 128},
    'openbind0': {'label': 'OpenBind-0', 'recycles': 3, 'diffusion_steps': 25, 'msa_depth': 128},
    'boltz2': {'label': 'Boltz-2', 'recycles': 3, 'diffusion_steps': 25, 'msa_depth': 128},
    'protenix1': {'label': 'Protenix V1', 'recycles': 3, 'diffusion_steps': 25, 'msa_depth': 128},
    'protenix2': {'label': 'Protenix-v2', 'recycles': 3, 'diffusion_steps': 25, 'msa_depth': 128},
    'mol1': {'label': 'Mol-1', 'recycles': 3, 'diffusion_steps': 25, 'msa_depth': 128},
    'rosettafold3': {'label': 'RoseTTAFold3', 'recycles': 3, 'diffusion_steps': 25, 'msa_depth': 128},
    'intellifold2': {'label': 'IntelliFold-2', 'recycles': 3, 'diffusion_steps': 25, 'msa_depth': 128},
    'opendde': {'label': 'OpenDDE', 'recycles': 3, 'diffusion_steps': 16, 'msa_depth': 128},
    'chai1': {'label': 'Chai-1 · ESM2', 'recycles': 3, 'diffusion_steps': 200, 'msa_depth': 16384},
    'esmfold2': {'label': 'ESMFold2 · Full · ESM-C 6B', 'recycles': 3, 'diffusion_steps': 15, 'msa_depth': 16384},
    'esmfold2_fast': {'label': 'ESMFold2 · Fast · ESM-C 6B', 'recycles': 3, 'diffusion_steps': 15, 'msa_depth': 16384},
    'esmfold2_lm600m': {'label': 'ESMFold2 · Fast · ESM-C 600M', 'recycles': 3, 'diffusion_steps': 15, 'msa_depth': 1},
    'esmfold2_lm300m': {'label': 'ESMFold2 · Fast · ESM-C 300M', 'recycles': 3, 'diffusion_steps': 15, 'msa_depth': 1},
    'bioemu2': {'label': 'BioEmu2', 'recycles': 3, 'diffusion_steps': 200, 'msa_depth': 4092},
    'helixfold3': {'label': 'HelixFold3 · non-commercial weights', 'recycles': 3, 'diffusion_steps': 200, 'msa_depth': 4096},
    'af2_ptm': {'label': 'AF2-mono', 'recycles': 3, 'diffusion_steps': None, 'msa_depth': 128,
                'extra_msa_depth': 256},
    'af2_multimer': {'label': 'AF2-multi', 'recycles': 3, 'diffusion_steps': None, 'msa_depth': 128,
                     'extra_msa_depth': 256},
}

WEB_MODEL_UI_GROUPS = (
    ('AF3 family', tuple((key, MODELS[key]['label']) for key in (
        'alphafold3', 'openbind0', 'opendde', 'boltz2',
        'protenix2', 'mol1', 'intellifold2', 'rosettafold3',
    ))),
    ('AlphaFold 2', (('af2_ptm', 'AF2-mono'), ('af2_multimer', 'AF2-multi'))),
    ('ESMFold2', (('esmfold2', 'ESMFold2'),)),
    ('BioEmu2', (('bioemu2', 'BioEmu2'),)),
)
MODEL_UI_GROUPS = WEB_MODEL_UI_GROUPS + (
    ('Additional JAX models', tuple((key, MODELS[key]['label']) for key in (
        'openfold3', 'protenix1', 'chai1', 'helixfold3',
    ))),
)
ESMFOLD2_VARIANTS = (
    ('esmfold2_lm600m', 'Fast · ESM-C 600M'),
    ('esmfold2_lm300m', 'Fast · ESM-C 300M'),
    ('esmfold2_fast', 'Fast · ESM-C 6B'),
    ('esmfold2', 'Full · ESM-C 6B'),
)
DEFAULT_ESMFOLD2_VARIANT = 'esmfold2_lm600m'

# These fingerprints identify the reviewed upstream terms text. Persisted
# acknowledgments are applied only while their fingerprints match this manifest.
TERMS_VERSIONS = {
    'alphafold3': {
        'label': 'AlphaFold 3 model parameters terms (2024-11-09) and prohibited-use policy',
        'sources': [
            {
                'url': 'https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_TERMS_OF_USE.md',
                'sha256': '41adf62ff5eabc58831c828793988537948663c139f8b87d8d413851b150b6e5',
            },
            {
                'url': 'https://github.com/google-deepmind/alphafold3/blob/main/WEIGHTS_PROHIBITED_USE_POLICY.md',
                'sha256': '4992aebdc29bc7a9260bed6373bad05bb5a7e1b783451cd170163dd259f46c45',
            },
        ],
    },
    'helixfold3': {
        'label': 'HelixFold3 non-commercial license',
        'repository_revision': '8e3991ab1209134b148b05d44e784a43eaa4484d',
        'sources': [
            {
                'url': 'https://github.com/PaddlePaddle/PaddleHelix/blob/dev/apps/protein_folding/helixfold3/LICENSE',
                'sha256': '7f426d0f3d9fab77f008c55af191605fd3d6b0510480fbcfe462362e371e2af5',
            },
        ],
    },
}


def _terms_acceptance_path() -> Path:
    from mn_molecular_design_studio.workspace import _studio_home
    return _studio_home() / 'localfold_terms_acceptance.json'


def record_terms_acceptance(user_id: str, models, *, accepted_at: str | None = None) -> Path:
    """Persist explicit per-user acceptance for the exact reviewed term versions."""
    user_id = str(user_id or '').strip()
    if not user_id:
        raise ValueError('A Studio user ID is required to store LocalFold terms acceptance.')
    selected = tuple(dict.fromkeys(str(model) for model in models))
    if not selected or set(selected) - set(TERMS_VERSIONS):
        raise ValueError('Choose AlphaFold 3 and/or HelixFold3 terms to record.')
    path = _terms_acceptance_path()
    try:
        payload = json.loads(path.read_text(encoding='utf-8'))
    except (OSError, json.JSONDecodeError):
        payload = {}
    if not isinstance(payload, dict):
        payload = {}
    users = payload.get('users')
    users = users if isinstance(users, dict) else {}
    user_record = users.get(user_id)
    user_record = user_record if isinstance(user_record, dict) else {}
    recorded_at = accepted_at or datetime.now(timezone.utc).isoformat(timespec='seconds')
    for model in selected:
        user_record[model] = {
            'terms_version': copy.deepcopy(TERMS_VERSIONS[model]),
            'accepted_at': recorded_at,
            'assertion': 'User confirmed review and acceptance, and confirmed this use complies with the terms.',
        }
    users[user_id] = user_record
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_suffix(path.suffix + '.tmp')
    temporary_path.write_text(json.dumps({'schema_version': 1, 'users': users}, indent=2), encoding='utf-8')
    temporary_path.replace(path)
    return path


def apply_user_terms_acceptance(request_payload: dict, *, user_id: str) -> dict:
    """Apply a saved acknowledgment to LocalFold API requests and save new opt-ins."""
    result = copy.deepcopy(request_payload)
    if result.get('workflow') != 'refolding.optimized_prediction':
        return result
    parameters = result.get('parameters')
    if not isinstance(parameters, dict) or parameters.get('engine') != 'localfold':
        return result
    options = parameters.get('engine_options')
    if not isinstance(options, dict):
        options = {}
        parameters['engine_options'] = options
    model = str(options.get('model') or 'alphafold3')
    if model not in TERMS_VERSIONS:
        return result

    path = _terms_acceptance_path()
    try:
        payload = json.loads(path.read_text(encoding='utf-8'))
    except (OSError, json.JSONDecodeError):
        payload = {}
    users = payload.get('users') if isinstance(payload, dict) else {}
    user_record = users.get(str(user_id or '')) if isinstance(users, dict) else None
    saved = user_record.get(model) if isinstance(user_record, dict) else None
    accepted = isinstance(saved, dict) and saved.get('terms_version') == TERMS_VERSIONS[model]

    if model == 'alphafold3':
        if parameters.get('accept_alphafold3_terms') is True and not accepted:
            record_terms_acceptance(user_id, (model,))
        elif accepted:
            parameters['accept_alphafold3_terms'] = True
    elif model == 'helixfold3':
        if options.get('accept_helixfold3_terms') is True and not accepted:
            record_terms_acceptance(user_id, (model,))
        elif accepted:
            options['accept_helixfold3_terms'] = True
    return result
# Keep the website's choices together and identify the JAX fork's additions.
UI_MODEL_KEYS = tuple(dict.fromkeys(
    [key for _, choices in MODEL_UI_GROUPS for key, _ in choices]
    + [key for key, _ in ESMFOLD2_VARIANTS]
))
CROSSLINKERS = ('DSS', 'BS3', 'DSG', 'BS2G', 'DSSO', 'DSBU', 'CDI', 'PHOX',
               'BSPEG5', 'BSPEG9', 'AZIDE-A-DSBSO', 'SDA', 'LCSDA', 'SDAD')


def normalize_restraints(raw, *, model, proteins, dna, rna, ligand_ccd_codes, ligand_smiles, task_index):
    """Validate covalent inputs and model-specific conditioning rows."""
    from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import _normalize_boltz_restraints
    if raw in (None, []):
        return []
    if not isinstance(raw, list):
        raise ValueError(f'Task {task_index} restraints must be a list.')
    result, occupied = [], set()
    protein_by_id = {row['id']: row['sequence'] for row in proteins}
    polymers = {row['id']: row['sequence'] for row in [*proteins, *dna, *rna]}
    ligand_ids = {row['id'] for row in [*ligand_ccd_codes, *ligand_smiles]}
    for row in raw:
        if not isinstance(row, dict) or len(row) != 1 or not set(row) <= {'bond', 'crosslink', 'contact', 'pocket'}:
            raise ValueError('Choose a supported LocalFold input row.')
        if 'contact' in row:
            if model not in {'bioemu2', 'boltz2', 'chai1'}:
                raise ValueError(f'LocalFold {model} does not support contact conditioning.')
            details = row['contact']
            if not isinstance(details, dict) or details.get('force'):
                raise ValueError('LocalFold contact conditioning does not support force potentials.')
            tokens = []
            for key in ('token1', 'token2'):
                token = details.get(key)
                if not isinstance(token, list) or len(token) not in (2, 3):
                    raise ValueError('Contact endpoints need a chain, residue and optional atom.')
                chain = token[0]
                if not isinstance(chain, str) or (chain not in polymers and chain not in ligand_ids):
                    raise ValueError('Contact endpoint refers to an unknown chain.')
                try:
                    position = int(token[1])
                except (TypeError, ValueError) as exc:
                    raise ValueError('Contact endpoints need whole-number residue positions.') from exc
                if str(position) != str(token[1]) or not 1 <= position <= (len(polymers[chain]) if chain in polymers else 1):
                    raise ValueError('Contact residue position is outside its chain.')
                atom = str(token[2] or '').strip() if len(token) == 3 else ''
                if chain in ligand_ids and not atom:
                    raise ValueError('A ligand contact must name an atom.')
                if model == 'chai1' and (chain not in protein_by_id or atom):
                    raise ValueError('Chai-1 contact endpoints must be standard protein residues.')
                if model == 'boltz2' and chain in polymers and atom:
                    raise ValueError('Boltz-2 polymer contacts select a residue, not an atom.')
                tokens.append([chain, position, *([atom] if atom else [])])
            try:
                distance = float(details.get('max_distance', 6))
            except (TypeError, ValueError) as exc:
                raise ValueError('Contact distance must be within 4–20 Å.') from exc
            if not 4 <= distance <= 20:
                raise ValueError('Contact distance must be within 4–20 Å.')
            result.append({'contact': {'token1': tokens[0], 'token2': tokens[1], 'max_distance': distance}})
            continue
        if 'pocket' in row:
            if model not in {'boltz2', 'chai1'}:
                raise ValueError(f'LocalFold {model} does not support pocket conditioning.')
            details = row['pocket']
            if not isinstance(details, dict) or details.get('force'):
                raise ValueError('LocalFold pocket conditioning does not support force potentials.')
            binder = details.get('binder')
            if not isinstance(binder, str) or binder not in polymers and binder not in ligand_ids:
                raise ValueError('Pocket binder refers to an unknown chain.')
            if model == 'chai1' and binder not in protein_by_id and binder not in {ligand['id'] for ligand in ligand_smiles}:
                raise ValueError('Chai-1 pocket binders must be a protein or a SMILES ligand.')
            contacts = details.get('contacts')
            if contacts is None:
                contacts = [line.strip().split(',') for line in str(details.get('contacts_input', '')).replace(';', '\n').splitlines() if line.strip()]
                contacts = [[part.strip() for part in contact] for contact in contacts]
            if not isinstance(contacts, list) or not contacts:
                raise ValueError('A pocket needs at least one contact residue or ligand atom.')
            # Pocket endpoints use the same native representation as contacts.
            # Reusing the contact validator keeps ligand atom positions intact.
            normalized_contacts = []
            for contact in contacts:
                pair = normalize_restraints(
                    [{'contact': {'token1': contact, 'token2': contact,
                                   'max_distance': details.get('max_distance', 6)}}],
                    model=model, proteins=proteins, dna=dna, rna=rna,
                    ligand_ccd_codes=ligand_ccd_codes, ligand_smiles=ligand_smiles,
                    task_index=task_index)[0]['contact']
                normalized_contacts.append(pair['token1'])
            result.append({'pocket': {'binder': binder, 'contacts': normalized_contacts,
                                      'max_distance': pair['max_distance']}})
            continue
        if 'bond' in row:
            result.extend(_normalize_boltz_restraints(
                [row], proteins=proteins, dna=dna, rna=rna, ligand_ccd_codes=ligand_ccd_codes,
                ligand_smiles=ligand_smiles, supported=True, task_index=task_index))
            continue
        details = row['crosslink']
        if not isinstance(details, dict) or details.get('name') not in CROSSLINKERS:
            raise ValueError('Choose a supported LocalFold crosslinker.')
        ends = []
        for key in ('residue1', 'residue2'):
            end = details.get(key)
            if not isinstance(end, list) or len(end) != 2 or end[0] not in protein_by_id:
                raise ValueError('Crosslinks need two protein chain and residue references.')
            try:
                position = int(end[1])
            except (ValueError, TypeError) as exc:
                raise ValueError('Crosslinks need whole-number residue positions.') from exc
            if str(position) != str(end[1]) or not 1 <= position <= len(protein_by_id[end[0]]):
                raise ValueError('Crosslink residue position is outside its protein sequence.')
            identity = (end[0], position)
            if identity in occupied:
                raise ValueError(f'Residue {end[0]}{position} carries two crosslinks.')
            occupied.add(identity)
            # The fork's NHS end accepts K/S/T/Y or the N-terminal amine.
            nhs_end = len(ends) == 0 or details['name'] not in {'SDA', 'LCSDA', 'SDAD'}
            if nhs_end and position != 1 and protein_by_id[end[0]][position - 1] not in 'KSTY':
                raise ValueError(f'Crosslink {details["name"]} needs K/S/T/Y or the N-terminus at {end[0]}{position}.')
            ends.append([end[0], position])
        result.append({'crosslink': {'name': details['name'], 'residue1': ends[0], 'residue2': ends[1]}})
    return result


def _select(key, label, default, choices, **kwargs):
    return dict(key=key, label=label, type='select', default=default, choices=choices, **kwargs)


SETTINGS = (
    _select('model', 'LocalFold model', 'alphafold3',
            tuple((key, model['label']) for key, model in MODELS.items())),
    dict(key='seed', label='Seed (0 for random)', default=0, min=0, max=2**31 - 1),
    _select('msa_mode', 'Alignment', 'automatic',
            (('automatic', 'Automatic · reuse cache, then MMseqs2'), ('reference_only', 'Reference files only'),
             ('single_sequence', 'Single sequence'))),
    dict(key='max_msa_seqs', label='Maximum MSA sequences', default=128, min=1, max=16384, step=1,
         help='Maximum MSA rows passed to the selected model. The model-specific default is selected with the model.'),
    dict(key='max_extra_msa_seqs', label='Maximum extra MSA sequences', default=256, min=1, max=16384, step=1,
         hidden=True, help='AF2 extra-MSA stack depth; shown only for AlphaFold 2.'),
    _select('weights_precision', 'Stored weight precision', 'int8', ('int8', 'fp16', 'fp32'),
            help='Select the exact stored format. Weights are expanded when loaded; this does not select arithmetic precision. Official AF3 and AF2 use their original parameters.'),
    dict(key='seed_count', label='Number of seeds', default=1, min=1, max=32),
    dict(key='num_sampling_steps', label='Diffusion steps (0 = model default)', default=0, min=0, max=2000),
    _select('flash_attention_implementation', 'Attention kernel', 'auto', ('auto', 'pallas', 'cudnn', 'triton', 'xla'),
            help='Auto lets the current LocalFold fork select a kernel for the allocated GPU.'),
    dict(key='dropout', label='Enable inference dropout', type='bool', default=False),
    dict(key='use_esm_embeddings', label='Use the model’s protein language model', type='bool', default=True),
    _select('af2_num_models', 'AF2 model number', 1,
            tuple((number, str(number)) for number in range(1, 6))),
    dict(key='save_embeddings', label='Save trunk embeddings', type='bool', default=False),
    dict(key='save_distogram', label='Save distogram', type='bool', default=False),
    _select('bioemu2_mode', 'BioEmu2 mode', 'structure', ('structure', 'equilibrium')),
    _select('bioemu2_sampler', 'BioEmu2 sampler', 'edm', ('edm', 'edm_gd', 'edm_fk', 'edm_fkc')),
    _select('bioemu2_physicality', 'BioEmu2 physicality report', 'high_quality',
            ('high_quality', 'no_violations', 'high_quality_per_residue', 'no_filter')),
    dict(key='bioemu2_lag_ns', label='BioEmu2 propagation lag (ns)', type='float', default=1., min=0.001, max=100000.),
    dict(key='bioemu2_propagation_steps', label='BioEmu2 propagation steps', default=1, min=1, max=10000),
    dict(key='bioemu2_contacts', label='BioEmu2 contacts (JSON)', type='text', default='', multiline=True),
    dict(key='bioemu2_potentials', label='BioEmu2 steering potentials (YAML)', type='text', default='', multiline=True),
    dict(key='bioemu2_structure_text', label='BioEmu2 initial structure (mmCIF)', type='text', default='', hidden=True),
    dict(key='accept_helixfold3_terms', label='I have reviewed and accept the HelixFold3 non-commercial weight terms',
         type='bool', default=False),
)


def cache_root(reference: Path) -> Path:
    """All persistent LocalFold data share this directory under reference_files."""
    return Path(reference).expanduser().resolve() / 'localfold'


def outputs_per_seed(options: dict, samples: int) -> int:
    if str(options.get('model', '')).startswith('af2_'):
        return int(options.get('af2_num_models', 1))
    frames = int(options.get('bioemu2_propagation_steps', 1)) if options.get('bioemu2_structure_text') else 1
    return int(samples) * frames


def prepare_references(reference: Path) -> Path:
    """Reuse compatible blobs without changing or copying the earlier caches.

    Relative links work on the host and under the container's /reference mount.
    Unpacked arrays and all subsequent downloads are written in LocalFold's own
    directories. Small ESMFold2 variants now have new trained confidence heads,
    so their earlier predictor blobs must not be linked into this installation.
    """
    reference = Path(reference).expanduser().resolve()
    root = cache_root(reference)
    for name in ('weights', 'ccd', 'runtime', 'bioemu2', 'home', 'xdg', 'huggingface'):
        (root / name).mkdir(parents=True, exist_ok=True)
    old_weights = reference / 'mn-cofolding/cache/alphafold3/weights'
    stale_legacy_blobs = {
        'rosettafold3-int8': 'rosettafold3.int8.bin.zst',
        'chai1-int8': 'chai1.int8.bin.zst',
    }
    if old_weights.is_dir():
        for directory in sorted(old_weights.iterdir()):
            if (not directory.is_dir()
                    or directory.name.startswith(('esmfold2_lm300m', 'esmfold2_lm600m'))
                    # The shared RosettaFold3 conversion predates a required
                    # confidence-head parameter, and the shared Chai-1
                    # conversion predates structure token-pair weights. Let
                    # LocalFold fetch current blobs into its own cache.
                    or directory.name in stale_legacy_blobs):
                continue
            for source in sorted(directory.iterdir()):
                if source.is_file() and (source.name.endswith('.bin.zst') or source.name.endswith('.lm.npz')):
                    destination = root / 'weights' / directory.name / source.name
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    if not destination.exists() and not destination.is_symlink():
                        destination.symlink_to(os.path.relpath(source, destination.parent))
    for source in (reference / 'alphafold3').glob('*.bin*'):
        if source.is_file():
            destination = root / 'weights/alphafold3' / source.name
            destination.parent.mkdir(parents=True, exist_ok=True)
            if not destination.exists() and not destination.is_symlink():
                destination.symlink_to(os.path.relpath(source, destination.parent))
    for source in (reference / 'alphafold_models').glob('params_model_*.npz'):
        destination = root / 'weights/alphafold2' / source.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        if not destination.exists() and not destination.is_symlink():
            destination.symlink_to(os.path.relpath(source, destination.parent))
    for model_dir, filename in stale_legacy_blobs.items():
        legacy_blob = old_weights / model_dir / filename
        local_blob = root / 'weights' / model_dir / filename
        if local_blob.is_symlink() and os.path.realpath(local_blob) == os.path.realpath(legacy_blob):
            # Keep each shared checkpoint intact; LocalFold's weight manager
            # will fetch a version compatible with the model code in its image.
            local_blob.unlink()
    return root


def profile_for(base: dict, model: str, *, ui: bool = False) -> dict:
    if model not in MODELS:
        raise ValueError(f'Unknown LocalFold model: {model!r}.')
    result = copy.deepcopy(base)
    af2 = model.startswith('af2_')
    esm = model.startswith('esmfold2')
    model_defaults = MODELS[model]
    result.update(default_recycles=model_defaults['recycles'], terms=model == 'alphafold3',
                  max_samples=1 if af2 else 1000,
                  default_samples=1 if af2 or model == 'bioemu2' else 5,
                  supports_templates=model != 'bioemu2' and not esm,
                  supports_cyclic=not af2 and not esm and model != 'bioemu2',
                  supports_modifications=not af2,
                  supports_bonds=not af2 and model not in {'esmfold2_lm600m', 'esmfold2_lm300m'},
                  supports_crosslinks=not af2 and not esm and model != 'bioemu2',
                  supports_contacts=model in {'boltz2', 'chai1', 'bioemu2'},
                  supports_pockets=model in {'boltz2', 'chai1'},
                  molecule_inputs=() if af2 else base['molecule_inputs'])
    no_msa = model in {'esmfold2_lm600m', 'esmfold2_lm300m'}
    result['supports_msa'] = not no_msa
    result['supports_msa_upload'] = not no_msa
    fields = []
    for field in result['settings']:
        key = field['key']
        if key == 'model':
            field.update(default=model, hidden=False)
            if ui:
                field['choices'] = tuple((key, MODELS[key]['label']) for key in UI_MODEL_KEYS)
        if key == 'msa_mode':
            field['default'] = 'single_sequence' if no_msa else 'automatic'
            field['hidden'] = no_msa
        if key == 'max_msa_seqs':
            field['default'] = model_defaults['msa_depth']
            field['hidden'] = no_msa
        if key == 'max_extra_msa_seqs':
            field['default'] = model_defaults.get('extra_msa_depth', 1024)
            field['hidden'] = not af2
        if key == 'num_sampling_steps':
            field['default'] = model_defaults['diffusion_steps'] or 0
            if ui and model_defaults['diffusion_steps'] is not None:
                steps = (15, 32, 64, 200) if esm else (25, 50, 100, 200)
                if model_defaults['diffusion_steps'] not in steps:
                    steps = tuple(sorted((*steps, model_defaults['diffusion_steps'])))
                field['type'] = 'select'
                field['choices'] = tuple((step, str(step)) for step in steps)
        if key == 'num_recycles':
            field['default'] = model_defaults['recycles']
        if ui:
            field['hidden'] = bool(field.get('hidden') or
                (key.startswith('bioemu2_') and model != 'bioemu2') or
                (key == 'af2_num_models' and not af2) or
                (key == 'num_sampling_steps' and af2) or
                (key in {'flash_attention_implementation', 'save_embeddings', 'save_distogram'} and af2) or
                (key in {'save_embeddings', 'save_distogram'} and model == 'bioemu2') or
                (key == 'dropout' and (esm or model == 'bioemu2')) or
                (key == 'bioemu2_contacts') or
                (key == 'weights_precision' and (af2 or model == 'alphafold3')) or
                (key == 'use_esm_embeddings' and model not in {'esmfold2', 'chai1'}) or
                (key == 'accept_helixfold3_terms' and model != 'helixfold3'))
            if key == 'accept_helixfold3_terms':
                if model == 'helixfold3':
                    field['required'] = True
                else:
                    field.pop('required', None)
        fields.append(field)
    result['settings'] = tuple(fields)
    notes = ['All shared weights, language-model data, CCD data, and compilation caches are stored under the configured reference_files directory.']
    if model == 'bioemu2':
        notes.append('BioEmu2 reports physicality for every sample. It has no confidence head; pLDDT, PAE, pTM and ipTM are unavailable.')
    if af2:
        notes.append('AF2 accepts proteins only. Samples are parameter sets × seeds; choose the number of AF2 parameter sets below.')
    if model == 'helixfold3':
        notes.append('HelixFold3 weights have non-commercial terms: https://github.com/PaddlePaddle/PaddleHelix/tree/dev/apps/protein_folding/helixfold3')
    result['fixed_settings'] = tuple(notes)
    return result


def api_contract(base: dict) -> dict:
    """Expose the same model defaults and capabilities to REST/MCP clients."""
    web_models = {key for _, choices in WEB_MODEL_UI_GROUPS for key, _ in choices}
    web_models.update(key for key, _ in ESMFOLD2_VARIANTS)
    models = {}
    for name, metadata in MODELS.items():
        profile = profile_for(base, name)
        fields = [copy.deepcopy(field) for field in profile['settings'] if field['key'] != 'seed']
        models[name] = {
            'label': metadata['label'], 'website_model': name in web_models,
            'defaults': {
                'samples': profile['default_samples'], 'num_recycles': profile['default_recycles'],
                'engine_options': {field['key']: field['default'] for field in fields},
            },
            'engine_option_fields': fields,
            'capabilities': {
                key.removeprefix('supports_'): bool(profile.get(key))
                for key in ('supports_msa', 'supports_msa_upload', 'supports_templates',
                            'supports_cyclic', 'supports_modifications', 'supports_bonds',
                            'supports_crosslinks', 'supports_contacts', 'supports_pockets')
            } | {'molecule_inputs': list(profile['molecule_inputs']),
                 'max_samples': profile['max_samples'],
                 'pae': name != 'bioemu2'},
            'terms': {'accept_alphafold3_terms': name == 'alphafold3',
                      'engine_options.accept_helixfold3_terms': name == 'helixfold3'},
        }
    return {'engine': 'localfold', 'default_model': 'alphafold3',
            'image': IMAGE, 'source_commit': SOURCE_COMMIT, 'models': models,
            'seed_parameter': 'parameters.seed',
            'input_format': 'studio_prediction_tasks', 'native_input_format': 'af3_json',
            'msa_default': 'automatic', 'msa_policy': 'reuse_exact_sequence_cache_then_mmseqs2',
            'model_groups': MODEL_UI_GROUPS, 'esmfold2_variants': ESMFOLD2_VARIANTS}


def validate_tasks(options: dict, tasks: list[dict], *, samples: int, seed: int) -> None:
    model = options['model']
    if seed + int(options['seed_count']) - 1 > 2**31 - 1:
        raise ValueError('The LocalFold seed range must not exceed 2,147,483,647.')
    if model == 'helixfold3' and not options['accept_helixfold3_terms']:
        raise ValueError('Review and accept HelixFold3’s published non-commercial weight terms before starting.')
    if model.startswith('af2_') and samples != 1:
        raise ValueError('AF2 uses one sample per parameter set. Choose af2_num_models instead.')
    cycles = {tuple(chain_id for chain_id, row in task['molecule_chain_options'].items() if row['cyclic']) for task in tasks}
    if len(cycles) > 1:
        raise ValueError('Use the same cyclic chain IDs for every LocalFold batch input, or submit separate batches.')
    if next(iter(cycles), ()) and (model.startswith(('af2_', 'esmfold2')) or model == 'bioemu2'):
        raise ValueError(f'LocalFold {model} does not support cyclic input.')
    if model.startswith('esmfold2') and not options['use_esm_embeddings'] and (not options.get('use_msa') or model in {'esmfold2_lm600m', 'esmfold2_lm300m'}):
        raise ValueError('ESMFold2 needs its matching language model or an MSA.')
    for task in tasks:
        kinds = {kind for row in task.get('restraints', []) for kind in row}
        if 'bond' in kinds and (model.startswith('af2_') or model in {'esmfold2_lm600m', 'esmfold2_lm300m'}):
            raise ValueError(f'LocalFold {model} does not support bond inputs.')
        if 'crosslink' in kinds and (model.startswith(('af2_', 'esmfold2')) or model == 'bioemu2'):
            raise ValueError(f'LocalFold {model} does not support crosslink inputs.')
        if model.startswith('af2_') and any(row['modifications'] for row in task['molecule_chain_options'].values()):
            raise ValueError(f'LocalFold {model} accepts unmodified proteins only.')
        if model == 'af2_ptm' and task['structure_templates'] and options['af2_num_models'] > 2:
            raise ValueError('AF2 pTM has two template-trained parameter sets; choose at most 2 when using templates.')
        if model in {'esmfold2_lm600m', 'esmfold2_lm300m'} and any(row['msa_a3m_text'] for row in task['molecule_chain_options'].values()):
            raise ValueError(f'LocalFold {model} has no MSA encoder; use its matching language model.')
    if model == 'bioemu2' and options['num_sampling_steps'] == 1:
        raise ValueError('BioEmu2 needs at least 2 diffusion steps, or 0 for its model default.')
    if options['bioemu2_contacts'].strip() and any('contact' in row for task in tasks for row in task.get('restraints', [])):
        raise ValueError('Use Contact rows or BioEmu2 contact JSON, not both.')
    for key in ('bioemu2_contacts', 'bioemu2_potentials'):
        value = options[key]
        if not value.strip():
            continue
        if model != 'bioemu2':
            raise ValueError(f'{key} is available only with the BioEmu2 model.')
        try:
            parsed = json.loads(value) if key.endswith('contacts') else yaml.safe_load(value)
        except (ValueError, yaml.YAMLError) as exc:
            raise ValueError(f'Invalid {key}: {exc}') from exc
        if key.endswith('contacts') and not isinstance(parsed, list):
            raise ValueError('BioEmu2 contacts must be a JSON list.')
        if key.endswith('potentials') and not parsed:
            raise ValueError('BioEmu2 potentials YAML must not be empty.')
    if options['bioemu2_structure_text']:
        if model != 'bioemu2' or len(tasks) != 1:
            raise ValueError('BioEmu2 propagation requires one input complex.')
        import gemmi
        try:
            document = gemmi.cif.read_string(options['bioemu2_structure_text'])
            if len(document) != 1 or not len(document.sole_block().find_loop('_atom_site.Cartn_x')):
                raise ValueError('No atom coordinates found.')
        except (ValueError, RuntimeError) as exc:
            raise ValueError('The BioEmu2 propagation structure must be a non-empty mmCIF file with atom coordinates.') from exc
    if model == 'bioemu2' and options['bioemu2_sampler'] in {'edm_fk', 'edm_fkc'} and options['bioemu2_propagation_steps'] > 1:
        raise ValueError('BioEmu2 FK samplers support only one propagation step.')
