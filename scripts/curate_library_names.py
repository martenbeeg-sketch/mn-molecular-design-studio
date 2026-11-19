"""Plan, back up and apply explicitly reviewed SQL naming/target corrections."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import json
from pathlib import Path
import sqlite3

from mn_molecular_design_studio.library_curation import apply_curation, encoded, plan_curation
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
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--references", type=Path)
    parser.add_argument("--backup", type=Path)
    parser.add_argument("--admin-id", required=True)
    args = parser.parse_args()
    access = args.root / "user_access.sqlite3"
    with sqlite3.connect(access.as_uri() + "?mode=ro", uri=True) as connection:
        user = connection.execute("SELECT id,role,disabled FROM users WHERE id=?", (args.admin_id,)).fetchone()
    if not user or user[1] != "admin" or user[2]:
        raise PermissionError("Select an existing, enabled administrator.")
    actor = {"id": user[0], "role": user[1]}
    if args.mode == "plan":
        references = json.loads(args.references.read_text()) if args.references else {}
        plan = plan_curation(ReadOnlyDatabase(args.root), actor, reference_names=references.get("names"),
            reference_constructs=references.get("constructs", []))
        if args.plan.exists():
            raise ValueError("Use a new plan filename so reviewed evidence is not overwritten.")
        args.plan.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        args.plan.write_text(encoded(plan))
        args.plan.chmod(0o600)
        report = {"digest": plan["digest"], "members": len(plan["members"]), "renamed": sum(row["before_name"] != row["after_name"] for row in plan["members"]),
                  "observations": len(plan["observations"]), "targets": len(plan["targets"]), "unresolved": len(plan["unresolved"])}
    else:
        if not args.backup:
            raise ValueError("An explicit new recovery backup path is required.")
        report = apply_curation(MolecularDatabase(args.root), actor, json.loads(args.plan.read_text()), backup_path=args.backup)
    print(encoded(report))


if __name__ == "__main__":
    main()
