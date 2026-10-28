"""Inspect a staged Studio checkout without activation or model inference."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import mn_molecular_design_studio as studio
from mn_molecular_design_studio.cluster import capabilities, config
from mn_molecular_design_studio.cluster.gpu_runtime import profiles, compatible
from mn_molecular_design_studio.engine.workflow_engine.core.central_scheduler import _client
from mn_molecular_design_studio.pipeline_graph import NODE_CATALOG
from mn_molecular_design_studio.workspace import ensure_studio_workspace, _studio_home


def main():
    ensure_studio_workspace()
    root = Path(studio.__file__).resolve().parent.parent
    tags = ['mn-rfd4-proteina:0.0.3-cu130', 'mn-boltz2-opt:2.2.1-cu130-batch2',
            'mn-protenix-opt:2.0.0-cu128-batch2']
    runtime_profiles = profiles(tags)
    details = capabilities.gpu_details()
    client = _client(timeout=3)
    files = [config.home() / 'config.json', _studio_home() / 'settings.json',
             *sorted((Path.home() / '.config/systemd/user').glob('studio-*.service'))]
    print(json.dumps(dict(
        evidence_class='Source inspection', python=sys.executable, root=str(root),
        git_head=subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip(),
        git_status=subprocess.check_output(['git', '-C', str(root), 'status', '--short'], text=True).strip(),
        source_fingerprint=capabilities.code_revision(), runtime_profiles=runtime_profiles,
        compatibility={gpu: {tag: compatible(info, {tag: profile}) for tag, profile in runtime_profiles.items()}
                       for gpu, info in details.items()},
        native_fields=[f['name'] for f in NODE_CATALOG['bindcraft2']['fields']],
        scheduler_job_counts=client.status().get('job_counts'), scheduler_resources=client.resources(),
        preserved_configuration_hashes={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files if p.is_file()},
    ), indent=2))


if __name__ == '__main__':
    main()
