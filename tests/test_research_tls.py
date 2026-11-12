import importlib.util
from pathlib import Path
import subprocess

import pytest


@pytest.fixture
def tls_fixture(tmp_path):
    spec = importlib.util.spec_from_file_location("research_tls", Path(__file__).resolve().parents[1] / "scripts/research_tls.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    tls = tmp_path / "tls"
    tls.mkdir()
    (tls / "server.cnf").write_text("""[req]
prompt=no
distinguished_name=dn
[dn]
CN=pc002975.marionegri.it
[server_extensions]
basicConstraints=critical,CA:FALSE
keyUsage=critical,digitalSignature
extendedKeyUsage=serverAuth
subjectAltName=DNS:pc002975.marionegri.it,DNS:pc002975,IP:172.18.24.170
""")
    module.command(["openssl", "req", "-x509", "-newkey", "ec", "-pkeyopt", "ec_paramgen_curve:P-256", "-nodes", "-days", "365",
                    "-keyout", str(tls / "research-ca.key"), "-out", str(tls / "research-ca.crt"), "-subj", "/CN=Disposable test root"])
    (tls / "research-ca.key").chmod(0o600)
    module.command(["openssl", "req", "-new", "-newkey", "ec", "-pkeyopt", "ec_paramgen_curve:P-256", "-nodes",
                    "-keyout", str(tls / "server.key"), "-out", str(tls / "server.csr"), "-config", str(tls / "server.cnf")])
    (tls / "server.key").chmod(0o600)
    module.command(["openssl", "x509", "-req", "-in", str(tls / "server.csr"), "-CA", str(tls / "research-ca.crt"),
                    "-CAkey", str(tls / "research-ca.key"), "-set_serial", "123", "-out", str(tls / "server.crt"), "-days", "5",
                    "-extfile", str(tls / "server.cnf"), "-extensions", "server_extensions"])
    (tmp_path / "Caddyfile").write_text("test configuration")
    (tmp_path / "caddy").write_text("#!/bin/sh\nexit 0\n")
    (tmp_path / "caddy").chmod(0o700)
    return module, tmp_path


def test_tls_renews_near_expiry_without_changing_the_root(tls_fixture):
    module, directory = tls_fixture
    authority = (directory / "tls/research-ca.crt").read_bytes()
    original_key = (directory / "tls/server.key").read_bytes()
    assert module.renew(directory, restart=False)
    assert (directory / "tls/research-ca.crt").read_bytes() == authority
    assert (directory / "tls/server.key").read_bytes() != original_key
    module.verify(directory / "tls/server.crt", directory / "tls/research-ca.crt")
    assert not module.renew(directory, restart=False)
    assert not list((directory / "tls").glob(".renew-*"))


def test_tls_validation_failure_restores_the_exact_previous_pair(tls_fixture):
    module, directory = tls_fixture
    original = {name: (directory / "tls" / name).read_bytes() for name in ("server.key", "server.crt")}
    (directory / "caddy").write_text("#!/bin/sh\nexit 1\n")
    with pytest.raises(subprocess.CalledProcessError):
        module.renew(directory, restart=False)
    for name, contents in original.items():
        assert (directory / "tls" / name).read_bytes() == contents


def test_tls_refuses_exposed_private_keys(tls_fixture):
    module, directory = tls_fixture
    (directory / "tls/research-ca.key").chmod(0o644)
    with pytest.raises(ValueError, match="owner"):
        module.renew(directory, restart=False)
