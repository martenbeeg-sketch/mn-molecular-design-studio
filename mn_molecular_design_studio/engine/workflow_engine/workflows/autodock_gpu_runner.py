"""Container-side AutoDock-GPU runner staged by the Studio tile."""

from __future__ import annotations

import csv
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


def _run(command: list[str], log_path: Path, *, cwd: Path) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, cwd=str(cwd), text=True, capture_output=True, check=False)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(
        "$ " + " ".join(command) + "\n" + result.stdout + result.stderr,
        encoding="utf-8",
    )
    if result.returncode:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {' '.join(command)}\n"
            f"{result.stderr[-3000:]}"
        )
    return result


def run_tile(root: Path = Path("/work")) -> dict:
    root = root.resolve()
    input_root = root / "artifacts" / "input"
    prepared_root = root / "artifacts" / "prepared"
    log_root = root / "artifacts" / "logs"
    pose_root = root / "artifacts" / "poses"
    ligand_root = root / "artifacts" / "ligands"
    for folder in (prepared_root, log_root, pose_root, ligand_root):
        folder.mkdir(parents=True, exist_ok=True)
    config = json.loads((input_root / "autodock_gpu_config.json").read_text(encoding="utf-8"))
    result_path = root / "artifacts" / "autodock_gpu_result.json"
    score_path = root / "artifacts" / "docking_scores.csv"

    receptor_base = prepared_root / "receptor"
    receptor_pdb = input_root / "receptor.pdb"
    prepare_command = [
        "mk_prepare_receptor.py", "--read_pdb", str(receptor_pdb),
        "-o", str(receptor_base), "-p", "-g", "--box_center",
        *map(str, config["center"]), "--box_size", *map(str, config["box_size"]),
    ]
    try:
        _run(prepare_command, log_root / "receptor_prepare.log", cwd=root)
    except RuntimeError:
        receptor_noh = input_root / "receptor-noh.pdb"
        receptor_noh.write_text(
            "\n".join(
                line for line in receptor_pdb.read_text(errors="replace").splitlines()
                if not (
                    line.startswith(("ATOM  ", "HETATM"))
                    and line[76:78].strip().upper() == "H"
                )
            ) + "\n",
            encoding="utf-8",
        )
        prepare_command[2] = str(receptor_noh)
        _run(prepare_command, log_root / "receptor_prepare_retry.log", cwd=root)

    gpf = receptor_base.with_suffix(".gpf")
    _run(
        ["autogrid4", "-p", gpf.name, "-l", str(log_root / "receptor_autogrid.glg")],
        log_root / "autogrid.log",
        cwd=prepared_root,
    )
    fld = receptor_base.with_suffix(".maps.fld")
    if not fld.is_file():
        raise RuntimeError("AutoGrid did not produce the receptor .maps.fld grid descriptor")

    with (input_root / "compounds.tsv").open(newline="", encoding="utf-8") as handle:
        compounds = list(csv.DictReader(handle, delimiter="\t"))
    rows = []
    for index, item in enumerate(compounds, start=1):
        compound_id = item["compound_id"]
        smiles = item["smiles"]
        source_sdf = ligand_root / f"{compound_id}_source.sdf"
        sdf = ligand_root / f"{compound_id}.sdf"
        ligand_pdbqt = ligand_root / f"{compound_id}.pdbqt"
        dock_dir = pose_root / compound_id
        dock_dir.mkdir(parents=True, exist_ok=True)
        try:
            _run(
                [
                    sys.executable, str(input_root / "smiles_to_3d.py"),
                    "--smiles", smiles, "--name", compound_id, "--output", str(source_sdf),
                ],
                log_root / f"{compound_id}.conformer.log",
                cwd=root,
            )
            if config.get("use_scrub", True):
                scrub_arguments = config.get("scrub_arguments")
                if scrub_arguments is None:
                    scrub_arguments = [
                        "--ph", f"{float(config.get('scrub_ph', 7.4)):.2f}",
                        "--cpu", "1", "--ff", "mmff94s", "--max_ff_iter", "200",
                        "--numconfs", "1",
                    ]
                    if config.get("scrub_skip_tautomer", True):
                        scrub_arguments.append("--skip_tautomers")
                scrub_command = [
                    str(config.get("scrub_executable") or "scrub.py"),
                    str(source_sdf), "-o", str(sdf),
                    *list(scrub_arguments),
                ]
                scrub_log = log_root / f"{compound_id}.scrub.log"
                try:
                    _run(scrub_command, scrub_log, cwd=root)
                    if not sdf.is_file() or sdf.stat().st_size == 0:
                        raise RuntimeError("Scrub produced no valid molecule")
                except RuntimeError as exc:
                    shutil.copy2(source_sdf, sdf)
                    with scrub_log.open("a", encoding="utf-8") as handle:
                        handle.write(f"\nScrub failed; using the RDKit SDF for Meeko preparation: {exc}\n")
            else:
                shutil.copy2(source_sdf, sdf)
            _run(
                ["mk_prepare_ligand.py", "-i", str(sdf), "-o", str(ligand_pdbqt)],
                log_root / f"{compound_id}.ligand_prepare.log",
                cwd=root,
            )
            resname = f"adgpu_{compound_id}"
            _run(
                [
                    "adgpu", "--ffile", str(fld), "--lfile", str(ligand_pdbqt),
                    "--nrun", str(config["nrun"]),
                    "--seed", str(int(config["seed"]) + index - 1),
                    "--gbest", "1", "--resnam", resname,
                ],
                log_root / f"{compound_id}.adgpu.log",
                cwd=dock_dir,
            )
            dlg_files = sorted(dock_dir.glob("*.dlg"))
            if not dlg_files:
                raise RuntimeError("AutoDock-GPU returned no DLG score report")
            score_values = [
                float(value)
                for value in re.findall(
                    r"Estimated Free Energy of Binding\s*=\s*(-?\d+(?:\.\d+)?)\s*kcal/mol",
                    dlg_files[0].read_text(errors="replace"),
                )
            ]
            if not score_values:
                raise RuntimeError("AutoDock-GPU DLG contains no estimated binding energy")
            pose_files = sorted(dock_dir.glob("*.pdbqt"))
            if not pose_files:
                raise RuntimeError("AutoDock-GPU returned no best-pose PDBQT; enable --gbest")
            canonical_pose = pose_root / f"{compound_id}.pdbqt"
            shutil.copy2(pose_files[0], canonical_pose)
            rows.append({
                "compound_id": compound_id,
                "status": "completed",
                "best_score_kcal_mol": min(score_values),
                "pose_path": str(canonical_pose.relative_to(root)),
                "error": "",
            })
        except Exception as exc:
            rows.append({
                "compound_id": compound_id,
                "status": "failed",
                "best_score_kcal_mol": "",
                "pose_path": "",
                "error": str(exc),
            })

    with score_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["compound_id", "status", "best_score_kcal_mol", "pose_path", "error"],
        )
        writer.writeheader()
        writer.writerows(rows)
    success = any(row["status"] == "completed" for row in rows)
    summary = {
        "success": success,
        "scored_count": sum(row["status"] == "completed" for row in rows),
        "compound_count": len(rows),
    }
    result_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    if not success:
        raise RuntimeError("AutoDock-GPU did not score any compounds")
    return summary


if __name__ == "__main__":
    try:
        run_tile()
    except Exception as exc:
        result_path = Path("/work/artifacts/autodock_gpu_result.json")
        result_path.parent.mkdir(parents=True, exist_ok=True)
        result_path.write_text(
            json.dumps({"success": False, "error": str(exc)}, indent=2) + "\n",
            encoding="utf-8",
        )
        raise
