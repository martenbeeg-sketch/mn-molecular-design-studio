from pathlib import Path

import pytest

from scripts.teaching_lan_gateway import configuration


@pytest.fixture
def tls_files(tmp_path):
    certificate = tmp_path / "teaching certificate.pem"
    private_key = tmp_path / "teaching key.pem"
    certificate.write_text("certificate placeholder")
    private_key.write_text("private key placeholder")
    private_key.chmod(0o600)
    return certificate, private_key


def test_gateway_exposes_only_https_ui_on_the_selected_private_interface(tls_files):
    rendered = configuration("172.18.24.170", "172.18.24.0/22", "pc002975.marionegri.it", *tls_files)
    assert "https://pc002975.marionegri.it:18502, https://172.18.24.170:18502" in rendered
    assert "bind 172.18.24.170" in rendered
    assert "not remote_ip 172.18.24.0/22" in rendered
    assert rendered.index("respond @outside") < rendered.index("reverse_proxy")
    assert "reverse_proxy 127.0.0.1:18502" in rendered
    assert 'tls "' in rendered
    assert "skip_install_trust" in rendered
    assert "admin off" in rendered
    assert "auto_https off" in rendered
    assert "protocols h1 h2" in rendered
    assert "max_size 10MB" in rendered
    assert "@api path /api /api/*" in rendered
    assert "18766" not in rendered
    assert "8766" not in rendered
    assert "8502" not in rendered.replace("18502", "")
    assert "0.0.0.0" not in rendered


@pytest.mark.parametrize("address,network", [
    ("0.0.0.0", "0.0.0.0/0"),
    ("127.0.0.1", "127.0.0.0/8"),
    ("172.18.24.170", "172.16.0.0/8"),
    ("8.8.8.8", "8.8.8.0/24"),
    ("172.18.24.170", "172.18.25.0/24"),
    ("172.18.24.0", "172.18.24.0/22"),
    ("172.18.27.255", "172.18.24.0/22"),
    ("::", "::/0"),
])
def test_gateway_refuses_wildcard_public_or_mismatched_networks(tls_files, address, network):
    with pytest.raises(ValueError):
        configuration(address, network, "pc002975.marionegri.it", *tls_files)


@pytest.mark.parametrize("hostname", ["", "https://pc002975", "pc002975:18502", "pc\nrespond 200", "-pc", "pc.", "x" * 64])
def test_gateway_rejects_hostname_configuration_injection(tls_files, hostname):
    with pytest.raises(ValueError):
        configuration("172.18.24.170", "172.18.24.0/22", hostname, *tls_files)


def test_gateway_refuses_group_readable_private_keys(tls_files):
    tls_files[1].chmod(0o640)
    with pytest.raises(ValueError, match="other users"):
        configuration("172.18.24.170", "172.18.24.0/22", "pc002975.marionegri.it", *tls_files)


def test_gateway_requires_existing_regular_certificate_files(tls_files):
    with pytest.raises(FileNotFoundError):
        configuration("172.18.24.170", "172.18.24.0/22", "pc002975.marionegri.it", Path("/missing-teaching-cert.pem"), tls_files[1])
    with pytest.raises(ValueError, match="regular files"):
        configuration("172.18.24.170", "172.18.24.0/22", "pc002975.marionegri.it", tls_files[0].parent, tls_files[1])
