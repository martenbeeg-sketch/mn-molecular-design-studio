"""Renew the institute-only research listener without changing its trust anchor."""
from __future__ import annotations

import argparse
import fcntl
import json
import os
from pathlib import Path
import secrets
import shutil
import stat
import subprocess
import tempfile


def command(arguments):
    return subprocess.run(arguments, check=True, capture_output=True, text=True, timeout=30)


def verify(certificate, authority):
    for option, name in (("-verify_hostname", "pc002975.marionegri.it"),
                         ("-verify_hostname", "pc002975"), ("-verify_ip", "172.18.24.170")):
        command(["openssl", "verify", "-CAfile", str(authority), option, name, str(certificate)])


def renew(directory, *, force=False, restart=True):
    directory = Path(directory).resolve(strict=True)
    tls = directory / "tls"
    if tls.is_symlink() or not tls.is_dir():
        raise ValueError("The TLS directory must not be an alias.")
    descriptor = os.open(tls / ".renew.lock", os.O_CREAT | os.O_WRONLY | os.O_NOFOLLOW, 0o600)
    with os.fdopen(descriptor, "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        return _renew(directory, force=force, restart=restart)


def _renew(directory, *, force, restart):
    tls = directory / "tls"
    authority, authority_key = tls / "research-ca.crt", tls / "research-ca.key"
    certificate, private_key = tls / "server.crt", tls / "server.key"
    for path in (authority, authority_key, certificate, private_key, tls / "server.cnf"):
        if path.is_symlink() or not path.is_file():
            raise ValueError("TLS configuration must contain regular files, not aliases.")
    for path in (authority_key, private_key):
        if stat.S_IMODE(path.stat().st_mode) & 0o077:
            raise ValueError("TLS private keys must be readable only by their owner.")
    command(["openssl", "x509", "-in", str(authority), "-noout", "-checkend", str(120 * 86400)])
    current = subprocess.run(["openssl", "x509", "-in", str(certificate), "-noout", "-checkend", str(21 * 86400)],
                             capture_output=True, text=True, timeout=30)
    if current.returncode == 0 and not force:
        verify(certificate, authority)
        return False
    with tempfile.TemporaryDirectory(prefix=".renew-", dir=tls) as temporary:
        temporary = Path(temporary)
        new_key, new_certificate, request = temporary / "server.key", temporary / "server.crt", temporary / "server.csr"
        command(["openssl", "req", "-new", "-newkey", "ec", "-pkeyopt", "ec_paramgen_curve:P-256", "-nodes",
                 "-keyout", str(new_key), "-out", str(request), "-config", str(tls / "server.cnf")])
        new_key.chmod(0o600)
        command(["openssl", "x509", "-req", "-in", str(request), "-CA", str(authority), "-CAkey", str(authority_key),
                 "-set_serial", "0x" + secrets.token_hex(16), "-out", str(new_certificate), "-days", "90", "-sha256",
                 "-extfile", str(tls / "server.cnf"), "-extensions", "server_extensions"])
        verify(new_certificate, authority)
        old_key, old_certificate = temporary / "previous.key", temporary / "previous.crt"
        shutil.copy2(private_key, old_key)
        shutil.copy2(certificate, old_certificate)
        try:
            new_key.replace(private_key)
            new_certificate.replace(certificate)
            command([str(directory / "caddy"), "validate", "--config", str(directory / "Caddyfile"), "--adapter", "caddyfile"])
            if restart:
                command(["systemctl", "--user", "try-restart", "mds-research-lan.service"])
        except BaseException:
            old_key.replace(private_key)
            old_certificate.replace(certificate)
            raise
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--no-restart", action="store_true")
    options = parser.parse_args()
    os.umask(0o077)
    renewed = renew(options.directory, force=options.force, restart=not options.no_restart)
    print(json.dumps({"renewed": renewed, "trust_anchor_changed": False}))


if __name__ == "__main__":
    main()
