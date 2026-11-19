"""Plan or apply source-verified MN-Ligand docking pose/receptor corrections."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sqlite3

from mn_molecular_design_studio.data_browser import ReadOnlyLibrary
from mn_molecular_design_studio.legacy_docking_links import apply_docking_links, plan_docking_links
from mn_molecular_design_studio.library_curation import encoded
from mn_molecular_design_studio.molecular_database import MolecularDatabase


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["plan", "apply"])
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--source-runs-root", type=Path)
    parser.add_argument("--run-id", action="append")
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--backup", type=Path)
    parser.add_argument("--admin-id", required=True)
    args = parser.parse_args()
    with sqlite3.connect((args.root / "user_access.sqlite3").resolve().as_uri() + "?mode=ro", uri=True) as connection:
        actor = connection.execute("SELECT id,role,disabled FROM users WHERE id=?", (args.admin_id,)).fetchone()
    if not actor or actor[1] != "admin" or actor[2]:
        raise PermissionError("Select an existing, enabled administrator.")
    user = {"id": actor[0], "role": actor[1]}
    if args.mode == "plan":
        if not args.source_runs_root or args.plan.exists():
            raise ValueError("Select the original MN-Ligand runs root and a new plan filename.")
        plan = plan_docking_links(ReadOnlyLibrary(args.root), user, args.source_runs_root, run_ids=args.run_id)
        args.plan.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
        with args.plan.open("x") as stream:
            stream.write(encoded(plan))
        args.plan.chmod(0o600)
        report = {"plan": plan["digest"], "examined": plan["examined"], "observations": len(plan["observations"]),
            "structures": len(plan["structures"]), "unresolved": len(plan["unresolved"])}
    else:
        if not args.backup:
            raise ValueError("A new, separate recovery backup is required.")
        report = apply_docking_links(MolecularDatabase(args.root), user, json.loads(args.plan.read_text()), backup_path=args.backup)
    print(encoded(report))


if __name__ == "__main__":
    main()
