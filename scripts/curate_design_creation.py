"""Review and back up source-owned design creation-time naming corrections."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import json
from pathlib import Path
import socket
import sqlite3

from mn_molecular_design_studio.design_creation import workstation_names
from mn_molecular_design_studio.library_curation import apply_curation, encoded
from mn_molecular_design_studio.library_design_creation import CreationRecovery, plan_creation_names
from mn_molecular_design_studio.molecular_database import MolecularDatabase


class ReadOnlyDatabase(MolecularDatabase):
    @contextmanager
    def connection(self):
        connection = sqlite3.connect(self.path.as_uri() + "?mode=ro", uri=True, timeout=60)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA query_only=ON")
        try:
            yield connection
        finally:
            connection.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["plan", "apply"])
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--runs-root", required=True, type=Path)
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--backup", type=Path)
    parser.add_argument("--admin-id", required=True)
    args = parser.parse_args()
    with sqlite3.connect((args.root / "user_access.sqlite3").as_uri() + "?mode=ro", uri=True) as connection:
        user = connection.execute("SELECT id,role,disabled FROM users WHERE id=?", (args.admin_id,)).fetchone()
    if not user or user[1] != "admin" or user[2]:
        raise PermissionError("Select an existing, enabled administrator.")
    actor = {"id": user[0], "role": user[1]}
    if args.mode == "plan":
        from mn_molecular_design_studio.cluster import config
        names = workstation_names()
        node = config.identity()
        recovery = CreationRecovery(names=names, local_node=node, local_name=names.get(node) or socket.gethostname(), runs_root=args.runs_root)
        plan = plan_creation_names(ReadOnlyDatabase(args.root), actor, recovery)
        if args.plan.exists():
            raise ValueError("Choose a new plan path; reviewed evidence must not be overwritten.")
        args.plan.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        args.plan.write_text(encoded(plan))
        args.plan.chmod(0o600)
        print(encoded({"digest": plan["digest"], "members": len(plan["members"]), "unresolved": len(plan["unresolved"])}))
    else:
        if not args.backup:
            raise ValueError("A new, separate recovery backup is required.")
        plan = json.loads(args.plan.read_text())
        if plan.get("schema") != "studio.library-design-creation.v1":
            raise ValueError("Select a reviewed design-creation naming plan.")
        print(encoded(apply_curation(MolecularDatabase(args.root), actor, plan, backup_path=args.backup)))


if __name__ == "__main__":
    main()
