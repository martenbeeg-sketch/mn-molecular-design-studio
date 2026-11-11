"""Render an HTTPS-only, subnet-restricted entry point for the teaching pilot."""
from __future__ import annotations

import argparse
import ipaddress
import json
from pathlib import Path
import re
import stat


def configuration(address, network, hostname, certificate, private_key):
    address = ipaddress.IPv4Address(address)
    network = ipaddress.IPv4Network(network, strict=True)
    private_ranges = [ipaddress.IPv4Network(value) for value in
                      ("10.0.0.0/8", "172.16.0.0/12", "192.168.0.0/16")]
    if not any(network.subnet_of(private) for private in private_ranges):
        raise ValueError("The teaching listener requires an explicit private LAN subnet.")
    if address not in network or address in (network.network_address, network.broadcast_address):
        raise ValueError("The listener must be a host address inside the allowed subnet.")
    if len(hostname) > 253 or not all(re.fullmatch(r"[a-zA-Z0-9](?:[a-zA-Z0-9-]{0,61}[a-zA-Z0-9])?", label)
                                     for label in hostname.split(".")):
        raise ValueError("Choose a DNS hostname, without a scheme, port or configuration syntax.")
    certificate = Path(certificate).resolve(strict=True)
    private_key = Path(private_key).resolve(strict=True)
    if not certificate.is_file() or not private_key.is_file():
        raise ValueError("The TLS certificate and private key must be regular files.")
    if stat.S_IMODE(private_key.stat().st_mode) & 0o077:
        raise ValueError("The TLS private key must not be readable by other users.")
    return f"""{{
    admin off
    auto_https off
    persist_config off
    skip_install_trust
    servers {{
        protocols h1 h2
        timeouts {{
            read_header 10s
            idle 2m
        }}
    }}
}}

https://{hostname}:18502, https://{address}:18502 {{
    bind {address}
    tls {json.dumps(str(certificate))} {json.dumps(str(private_key))}
    route {{
        @outside not remote_ip {network}
        respond @outside "Institute teaching subnet only" 403
        @api path /api /api/*
        respond @api "Not found" 404
        request_body {{
            max_size 10MB
        }}
        header {{
            X-Content-Type-Options nosniff
            Referrer-Policy no-referrer
        }}
        reverse_proxy 127.0.0.1:18502
    }}
}}
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--listen-address", required=True)
    parser.add_argument("--allow-network", required=True)
    parser.add_argument("--hostname", required=True)
    parser.add_argument("--certificate", type=Path, required=True)
    parser.add_argument("--private-key", type=Path, required=True)
    options = parser.parse_args()
    try:
        rendered = configuration(options.listen_address, options.allow_network, options.hostname,
                                 options.certificate, options.private_key)
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(rendered, end="")


if __name__ == "__main__":
    main()
