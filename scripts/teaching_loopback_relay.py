"""Forward only two loopback teaching ports into a network-disabled container."""
from __future__ import annotations

import argparse
import re
import select
import signal
import socket
import subprocess
import threading


RELAY = """
import os, select, socket, sys
connection = socket.create_connection(('127.0.0.1', int(sys.argv[1])), timeout=15)
connection.settimeout(None)
inputs = [0, connection]
while inputs:
    readable, _, _ = select.select(inputs, [], [])
    for source in readable:
        data = os.read(0, 4096) if source == 0 else connection.recv(4096)
        if not data:
            if source == 0:
                inputs.remove(0)
                connection.shutdown(socket.SHUT_WR)
                continue
            sys.exit(0)
        if source == 0:
            connection.sendall(data)
        else:
            os.write(1, data)
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--container", required=True)
    parser.add_argument("--python", required=True)
    options = parser.parse_args()
    if not re.fullmatch(r"mds-teaching-[a-zA-Z0-9_.-]+", options.container):
        parser.error("Only an explicitly named teaching container is supported.")
    stopped = threading.Event()
    slots = threading.BoundedSemaphore(32)
    for signum in (signal.SIGTERM, signal.SIGINT):
        signal.signal(signum, lambda *arguments: stopped.set())

    def forward(client, port):
        process = None
        try:
            client.settimeout(30)
            process = subprocess.Popen(["docker", "exec", "-i", options.container, options.python,
                "-u", "-c", RELAY, str(port)], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL, bufsize=0)
            while not stopped.is_set():
                readable, _, _ = select.select([client, process.stdout], [], [], 1)
                for source in readable:
                    data = client.recv(4096) if source is client else source.read(4096)
                    if not data:
                        return
                    if source is client:
                        process.stdin.write(data)
                    else:
                        client.sendall(data)
        except (OSError, subprocess.SubprocessError):
            pass
        finally:
            client.close()
            if process is not None:
                process.stdin.close()
                process.stdout.close()
                if process.poll() is None:
                    process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            slots.release()

    listeners, workers = [], []
    try:
        for port in (18502, 18766):
            server = socket.socket()
            server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            server.bind(("127.0.0.1", port))
            server.listen(16)
            listeners.append(server)
        while not stopped.is_set():
            ready, _, _ = select.select(listeners, [], [], 1)
            for server in ready:
                client, _ = server.accept()
                if slots.acquire(blocking=False):
                    workers = [worker for worker in workers if worker.is_alive()]
                    worker = threading.Thread(target=forward, args=(client, server.getsockname()[1]), daemon=True)
                    workers.append(worker)
                    worker.start()
                else:
                    client.close()
    finally:
        stopped.set()
        for server in listeners:
            server.close()
        for worker in workers:
            worker.join(timeout=6)


if __name__ == "__main__":
    main()
