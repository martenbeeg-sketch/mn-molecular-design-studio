import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

api_pid = subprocess.check_output(
    ['systemctl', '--user', 'show', 'studio-api.service', '-p', 'MainPID', '--value'], text=True
).strip()
proc = Path('/proc', api_pid)
source = (proc / 'cwd').resolve()
environment = dict(entry.split('=', 1) for entry in (proc / 'environ').read_bytes().decode().split('\0') if '=' in entry)
python = str((proc / 'exe').resolve())
prefix = 'mn_molecular_design_studio.engine.workflow_engine.core.'
scheduler = prefix + 'central_scheduler'
worker = prefix + 'central_scheduler_worker'
def run(code):
    return subprocess.check_output([python, '-c', code], cwd=source, env=environment, text=True, timeout=20).strip()
status = json.loads(run(f'import json; from {scheduler} import _client; print(json.dumps(_client().status()["job_counts"]))'))
if any(status.get(state, 0) for state in ('running','reserved','leased','dispatching','starting','cancelling')):
    print(os.uname().nodename, 'busy', status, flush=True)
    raise SystemExit(2)
owner = Path(run(f'from {scheduler} import _service_dir; print(_service_dir() / "service.lock" / "owner.json")'))
old = json.loads(owner.read_text()).get('pid') if owner.is_file() else None
if old and Path('/proc', str(old)).is_dir():
    path = Path('/proc', str(old))
    args = (path/'cmdline').read_bytes().split(b'\0')
    assert (path/'status').stat().st_uid == os.getuid() and worker.encode() in args
    if (path/'cwd').resolve() != source:
        os.kill(old, signal.SIGTERM)
        deadline = time.monotonic() + 40
        while owner.exists() and time.monotonic() < deadline:
            time.sleep(0.3)
        assert not owner.exists(), 'Bridge did not finish its graceful stop.'
run(f'from {scheduler} import ensure_scheduler_worker_service; ensure_scheduler_worker_service()')
deadline = time.monotonic() + 20
while time.monotonic() < deadline:
    if owner.is_file():
        try:
            pid = json.loads(owner.read_text())['pid']
            if (Path('/proc', str(pid), 'cwd').resolve() == source):
                print(os.uname().nodename, pid, str(source), flush=True)
                sys.exit(0)
        except (OSError, ValueError, KeyError):
            pass
    time.sleep(0.3)
raise RuntimeError('Unable to verify the Scheduler bridge revision.')
