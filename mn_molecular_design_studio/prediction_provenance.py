"""Content and installed-container identities for reproducible prediction evidence."""
import hashlib
import json
from pathlib import Path
import subprocess


def file_identity(path):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return {'path': str(path), 'sha256': digest.hexdigest(), 'bytes': path.stat().st_size}


def image_identity(image):
    try:
        result = subprocess.run(['docker', 'image', 'inspect', image, '--format', '{{json .Id}}'],
                                capture_output=True, text=True, check=True, timeout=30)
        digest = json.loads(result.stdout)
        return {'image': image, 'image_digest': digest if isinstance(digest, str) else None}
    except (OSError, subprocess.SubprocessError, ValueError, TypeError, AttributeError):
        return {'image': image, 'image_digest': None}
