"""Content provenance shared by graph inputs, candidate artifacts and reuse checks."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


def file_digest(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def referenced_files(value, *, base=None):
    """Hash actual files, including relative references in JSON/YAML conditions.

    CSV/FASTA/structures/checkpoints are streamed, never loaded as text tables.
    Directories are not silently treated as immutable checkpoint identities.
    """
    result, seen = {}, set()

    def visit(item, directory):
        if isinstance(item, dict):
            for child in item.values():
                visit(child, directory)
        elif isinstance(item, (list, tuple)):
            for child in item:
                visit(child, directory)
        elif isinstance(item, (str, Path)) and item:
            try:
                path = Path(item).expanduser()
                if not path.is_absolute() and directory:
                    path = directory / path
                if not path.is_file():
                    return
                path = path.resolve()
            except (OSError, ValueError):
                return
            if path in seen:
                return
            seen.add(path)
            result[str(path)] = file_digest(path)
            # Conditions and candidate collections can themselves reference inputs.
            if path.suffix.lower() in {'.json', '.jsonl', '.yaml', '.yml'}:
                try:
                    if path.suffix == '.jsonl':
                        with path.open() as handle:
                            for line in handle:
                                if line.strip():
                                    visit(json.loads(line), path.parent)
                    elif path.suffix == '.json':
                        visit(json.loads(path.read_text()), path.parent)
                    else:
                        import yaml
                        visit(yaml.safe_load(path.read_text()), path.parent)
                except (ValueError, UnicodeError):
                    # A malformed input's bytes still invalidate previous results;
                    # scientific validation remains the adapter's responsibility.
                    pass
    visit(value, Path(base) if base else None)
    return result


def execution_identity(kind, params, workflow):
    """Resolve installed runtime identities before looking up a completed job."""
    import inspect
    import subprocess
    import types
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS, _get_function
    from mn_molecular_design_studio import mpnn_input

    implementation = _get_function(WORKFLOWS[workflow])
    module = inspect.getmodule(implementation)
    identity = {'implementation': file_digest(inspect.getfile(implementation))}
    images = set()
    modules = [module, *(m for m in vars(module).values() if isinstance(m, types.ModuleType)
        and m.__name__.startswith('mn_molecular_design_studio'))]
    for dependency in modules:
        if dependency.__file__:
            identity[dependency.__name__] = file_digest(dependency.__file__)
        profiles = getattr(dependency, 'ENGINE_PROFILES', {})
        engine = params.get('engine', kind)
        if engine in profiles:
            profile = profiles[engine]
            identity['engine_profile'] = profile
            if profile.get('image'):
                images.add(profile['image'])
        for name, value in vars(dependency).items():
            label = dependency.__name__ + ':' + name
            if name.isupper() and isinstance(value, str) and 'IMAGE' in name:
                images.add(value)
            elif name.isupper() and isinstance(value, (str, int)) and any(s in name for s in ('REVISION', 'COMMIT', 'VERSION')):
                identity[label] = value
            elif name.isupper() and isinstance(value, Path) and value.suffix in {'.py', '.sh', '.json', '.yaml', '.yml'} and value.is_file():
                identity[label] = file_digest(value)
    from mn_molecular_design_studio import generator_input
    if kind in generator_input.ENGINES:
        images.add(generator_input.ENGINES[kind]['image'])
        try:
            prepared = generator_input.prepare_request(kind, params)
        except ValueError:
            # Invalid scientific requests are reported by the node adapter;
            # provenance collection must not abort unrelated graph branches.
            prepared = generator_input.default_parameters(kind, params.get('task', 'binder')) | params
        identity['effective_generator_parameters'] = {k: v for k, v in prepared.items() if k not in {'target_pdb', 'campaign_name', 'seed'}}
        required = generator_input.readiness(kind, prepared)['required']
        weights = []
        for value in required:
            path = Path(value)
            # Include checkpoint/config files in declared model directories.
            weights.extend([p for p in path.rglob('*') if p.is_file()] if path.is_dir() else [path])
        identity['generator_checkpoints'] = {str(p): file_digest(p) if p.is_file() else None for p in weights}
    manifests = Path(inspect.getfile(implementation)).parents[1] / 'manifests'
    # Tool manifests are also execution inputs, including mutable image tags.
    for name in {kind, kind.removesuffix('_generator'), 'ligandmpnn' if workflow == 'sequence.mpnn' else kind}:
        manifest = manifests / (name + '.yaml')
        if manifest.is_file():
            import yaml
            data = yaml.safe_load(manifest.read_text())
            identity['manifest:' + name] = file_digest(manifest)
            if data.get('image'):
                images.add(data['image'])
    if workflow in {'sequence.mpnn', 'sequence.abmpnn'}:
        from mn_molecular_design_studio.engine.workflow_engine.workflows import sequence_design
        identity['native_adapter'] = file_digest(sequence_design.__file__)
        identity['input_adapter'] = file_digest(mpnn_input.__file__)
        model = 'ab_mpnn' if workflow == 'sequence.abmpnn' else params.get('model_type', 'protein_mpnn')
        if params.get('sequence_stage') == 'initial':
            model = 'protein_mpnn'
        elif params.get('sequence_stage') == 'soluble_optimization':
            model = 'soluble_mpnn'
        checkpoint = mpnn_input.checkpoint(model, params.get('noise_level') or None)
        identity['checkpoint'] = {'path': str(checkpoint), 'sha256': file_digest(checkpoint) if checkpoint.is_file() else None}
        images.add(mpnn_input.IMAGE)
    from mn_molecular_design_studio.runtime import reference_root
    checkpoints = []
    if kind == 'rfd3':
        checkpoints = [reference_root() / 'foundry/rfd3_latest.ckpt']
    elif kind in {'rf3', 'boltz2'}:
        checkpoints = ([reference_root() / 'foundry/rf3_foundry_01_24_latest_remapped.ckpt'] if kind == 'rf3'
            else [reference_root() / 'boltz_models' / name for name in ('boltz2_conf.ckpt', 'boltz2_aff.ckpt', 'ccd.pkl')])
    elif any(hasattr(m, 'HF_REVISION') and 'rfd4' in m.__name__ for m in modules):
        inputs = next(m for m in modules if hasattr(m, 'HF_REVISION') and 'rfd4' in m.__name__)
        snapshot = reference_root() / 'rfd4-proteina/huggingface/hub/models--nvidia--RFD4-Proteina/snapshots' / inputs.HF_REVISION
        checkpoints = [snapshot / name for name in ('model.safetensors', 'config.json')]
    elif kind == 'complexa_generator':
        from mn_molecular_design_studio import complexa_input
        task = params.get('pipeline_task', 'protein_binder')
        checkpoints = [complexa_input.resource_paths()['checkpoints'] / name for name in complexa_input.CHECKPOINTS[task]]
    elif kind == 'esm_binder':
        from mn_molecular_design_studio.engine.workflow_engine.workflows import esm_binder
        model = params.get('model_name') or esm_binder.DEFAULT_BINDER_MODEL
        for directory in (esm_binder.ESMFOLD2_BINDER_MODEL_ROOT / model, esm_binder.ESMC_MODEL_DIR):
            checkpoints += [p for p in directory.rglob('*') if p.is_file() and p.suffix in {'.safetensors', '.bin', '.pt', '.json'}]
    elif kind in {'af2', 'bindcraft2', 'bindcraft2_generator'}:
        directory = reference_root() / 'alphafold_models'
        checkpoints = [p for p in directory.rglob('*') if p.is_file() and p.suffix in {'.npz', '.json'}]
    if checkpoints:
        identity['installed_checkpoints'] = {str(p): file_digest(p) if p.is_file() else None for p in checkpoints}
    identity['images'] = {}
    for image in sorted(images):
        try:
            result = subprocess.run(['docker', 'image', 'inspect', '--format', '{{.Id}}', image], capture_output=True, text=True, timeout=10)
            identity['images'][image] = result.stdout.strip() if result.returncode == 0 else 'unavailable'
        except (OSError, subprocess.SubprocessError):
            identity['images'][image] = 'unavailable'
    return identity
