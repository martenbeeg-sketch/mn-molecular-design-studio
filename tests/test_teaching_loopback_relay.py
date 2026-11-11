import socket
import subprocess
import sys
import threading

import pytest

from scripts.teaching_loopback_relay import RELAY, main


@pytest.mark.parametrize("payload", [b"GET /health HTTP/1.1\r\n\r\n", bytes(range(256)) * 1024], ids=["http", "binary-256k"])
def test_teaching_stdio_tunnel_preserves_bytes_and_half_close(payload):
    with socket.socket() as server:
        server.bind(("127.0.0.1", 0))
        server.listen(1)
        def echo():
            connection, _ = server.accept()
            with connection:
                while True:
                    data = connection.recv(4096)
                    if not data:
                        break
                    connection.sendall(data)
        worker = threading.Thread(target=echo, daemon=True)
        worker.start()
        result = subprocess.run([sys.executable, "-u", "-c", RELAY, str(server.getsockname()[1])],
                                input=payload, capture_output=True, timeout=15)
        worker.join(timeout=5)
    assert result.returncode == 0, result.stderr
    assert result.stdout == payload
    assert not worker.is_alive()


def test_teaching_relay_cannot_target_a_production_named_container(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["relay", "--container", "studio-api", "--python", sys.executable])
    with pytest.raises(SystemExit) as refused:
        main()
    assert refused.value.code == 2
