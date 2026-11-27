"""Run MODELLER loop refinement on an observed residue range."""
from __future__ import annotations

import json
import sys
from pathlib import Path

from modeller import Environ, Selection, log
from modeller.automodel import LoopModel, assess, refine


def run(path: str, options: dict) -> None:
    log.none()
    chain = "" if options["chain"] == "_" else options["chain"]

    class SelectedLoop(LoopModel):
        def select_loop_atoms(self):
            return Selection(self.residue_range(f"{options['start']}:{chain}", f"{options['end']}:{chain}"))

    model = SelectedLoop(Environ(), inimodel=path, sequence="loop", loop_assess_methods=assess.DOPE)
    model.loop.ending_model = options["models"]
    model.loop.md_level = refine.fast
    model.make()
    results = [
        {"name": str(row["name"]), "dope": row.get("DOPE score"), "objective": row.get("molpdf")}
        for row in model.loop.outputs if not row.get("failure") and row.get("name")
    ]
    Path("models.json").write_text(json.dumps(results), encoding="utf-8")


if __name__ == "__main__":
    run(sys.argv[1], json.loads(sys.argv[2]))
