"""Read installed CUDA binary targets without allocating GPUs or loading models."""
from __future__ import annotations
import json
import subprocess
from functools import lru_cache

_PROBE = '''import json
try:
 import torch
 print(json.dumps({'cuda':torch.version.cuda, 'targets':str(torch._C._cuda_getArchFlags() or '').split()}))
except (ImportError, AttributeError):
 print(json.dumps({'cuda':None,'targets':[]}))
'''


@lru_cache(maxsize=64)
def _probe(image_id):
    try:
        result = subprocess.run(['docker', 'run', '--rm', '--network', 'none', '--entrypoint', 'python', image_id, '-c', _PROBE],
                                capture_output=True, text=True, check=True, timeout=30)
        return json.loads(result.stdout.strip().splitlines()[-1])
    except (OSError, subprocess.SubprocessError, ValueError, IndexError):
        return {'cuda': None, 'targets': []}


def profiles(images):
    """Inspect only requested installed images; never pull or probe model weights."""
    result = {}
    for image in images:
        try:
            identity = subprocess.run(['docker', 'image', 'inspect', '--format', '{{.Id}}', image],
                capture_output=True, text=True, check=True, timeout=5).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            continue
        profile = _probe(identity)
        if profile.get('targets'):
            result[image] = profile
    return result


def compatible(gpu, profiles):
    """Apply compiled CUDA cubin/PTX targets, separately from VRAM admission.

    CUDA cubins run on the same major architecture with a newer minor; PTX
    targets permit forward JIT. Unknown inventories fail this constraint.
    """
    if not profiles:
        return True
    try:
        capability = int(round(float(gpu['compute_capability']) * 10))
    except (KeyError, ValueError, TypeError):
        return False
    for profile in profiles.values():
        if profile.get('cuda'):
            try:
                major, minor = (int(v) for v in profile['cuda'].split('.')[:2])
                if int(gpu.get('driver_cuda_api') or 0) < major * 1000 + minor * 10:
                    return False
            except (ValueError, TypeError):
                return False
        targets = profile.get('targets') or []
        cubins = [int(v.split('_')[1]) for v in targets if v.startswith('sm_') and v.split('_')[1].isdigit()]
        ptx = [int(v.split('_')[1]) for v in targets if v.startswith('compute_') and v.split('_')[1].isdigit()]
        if not (any(v // 10 == capability // 10 and v <= capability for v in cubins) or any(v <= capability for v in ptx)):
            return False
    return True


def request_images(request):
    workflow = request.get('workflow', '')
    if workflow == 'refolding.optimized_prediction':
        from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES
        return [ENGINE_PROFILES[request['parameters']['engine']]['image']]
    if workflow in {'sequence.mpnn', 'sequence.abmpnn'}:
        from mn_molecular_design_studio import mpnn_input
        return [mpnn_input.IMAGE] if request['parameters'].get('backend') != 'foundry' else ['mn-foundry:cu128']
    from mn_molecular_design_studio.engine.workflow_engine.services.workflow_automation import WORKFLOWS
    from pathlib import Path
    import yaml
    spec = WORKFLOWS.get(workflow)
    if spec:
        path = Path(__file__).resolve().parents[1] / 'engine/workflow_engine/manifests' / (spec.tool + '.yaml')
        if path.is_file():
            value = yaml.safe_load(path.read_text())
            if value.get('gpu') and value.get('image'):
                return [value['image']]
    return []
