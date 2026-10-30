"""Prepare existing designs for prediction and compare the resulting structures.

No engine is started by this module. Comparisons use sequence-matched C-alpha
atoms; complex pose RMSD applies the target fit to the binder without refitting it.
"""
from __future__ import annotations

import copy
import csv
import gzip
import io
from pathlib import Path
from typing import Any

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFIO, Select
from Bio.PDB.PDBExceptions import PDBException
from Bio.SeqUtils import seq1

from mn_molecular_design_studio.engine.workflow_engine.core.candidates import _parse_structure, read_candidates
from mn_molecular_design_studio.engine.workflow_engine.core.jobs import read_json
from mn_molecular_design_studio.engine.workflow_engine.core.portable_paths import resolve_stored_path

CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"


def read_structure(path: Path):
    # Native BindCraft2 CIFs omit occupancy; use the same parser repair as the
    # normalized candidate reader, including gzip support.
    return _parse_structure(Path(path))


def protein_chains(structure) -> dict[str, list]:
    model = next(structure.get_models(), None)
    if model is None:
        raise ValueError("The structure has no coordinate model.")
    return {chain.id: residues for chain in model
            if (residues := [r for r in chain if "CA" in r and seq1(r.resname, undef_code="X") != "X"])}


def sequence(residues: list) -> str:
    return "".join(seq1(r.resname, undef_code="X") for r in residues)


def interface_filter_measurements(structure_path: Path, binder_chains: list[str],
                                  target_chains: list[str]) -> dict[str, int]:
    """Measure the two geometry filters shared with BindCraft2.

    Interface size is the number of binder residues with any heavy atom within
    4 Å of any target heavy atom. Backbone clashes count interchain Cα pairs
    within 2.5 Å. Missing chain roles or coordinates raise instead of silently
    manufacturing a passing zero.
    """
    from Bio.PDB import NeighborSearch

    chains = protein_chains(read_structure(Path(structure_path)))
    binder_chains = [str(chain) for chain in binder_chains]
    target_chains = [str(chain) for chain in target_chains]
    if not binder_chains or not target_chains:
        raise ValueError("Interface measurements require identified binder and target chains.")
    if any(chain not in chains for chain in (*binder_chains, *target_chains)):
        raise ValueError("Interface measurement chain roles do not match the predicted structure.")

    def heavy(atom):
        element = str(atom.element or "").strip().upper()
        return element not in {"H", "D"} and not atom.get_name().strip().upper().startswith(("H", "D"))

    target_atoms = [atom for chain in target_chains for residue in chains[chain]
                    for atom in residue if heavy(atom)]
    binder_residues = [residue for chain in binder_chains for residue in chains[chain]]
    if not target_atoms:
        raise ValueError("The predicted target has no heavy atoms for interface measurement.")
    target_search = NeighborSearch(target_atoms)
    interface_residues = 0
    for residue in binder_residues:
        atoms = [atom for atom in residue if heavy(atom)]
        if any(any(float(np.linalg.norm(atom.coord - target.coord)) <= 4.0
                   for target in target_search.search(atom.coord, 4.000001, level="A"))
               for atom in atoms):
            interface_residues += 1

    ca_by_chain = {chain: [residue["CA"] for residue in rows if "CA" in residue]
                   for chain, rows in chains.items()}
    chain_ids = list(ca_by_chain)
    clashes = 0
    for index, left_chain in enumerate(chain_ids):
        left_atoms = ca_by_chain[left_chain]
        if not left_atoms:
            continue
        for right_chain in chain_ids[index + 1:]:
            right_atoms = ca_by_chain[right_chain]
            if not right_atoms:
                continue
            right_search = NeighborSearch(right_atoms)
            for atom in left_atoms:
                clashes += sum(float(np.linalg.norm(atom.coord - other.coord)) <= 2.5
                               for other in right_search.search(atom.coord, 2.500001, level="A"))
    return {"validation_interface_residue_count_4a": interface_residues,
            "validation_backbone_ca_clashes_2_5a": clashes}


def residue_pairs(reference: str, predicted: str) -> list[tuple[int, int]]:
    aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                              open_gap_score=-8, extend_gap_score=-.5)
    alignment = next(iter(aligner.align(reference, predicted)))
    return [(int(i), int(j)) for ref, model in zip(*alignment.aligned)
            for i, j in zip(range(*ref), range(*model))]


def _candidate_structure(run_dir: Path, candidate: dict) -> Path:
    raw = candidate.get("raw_metadata") or {}
    roots = [run_dir]
    if raw.get("source_run_dir"):
        source = resolve_stored_path(raw["source_run_dir"], must_exist=True)
        if source:
            roots.append(source)
    for key in ("complex_pdb", "binder_pdb"):
        for root in roots:
            path = resolve_stored_path(candidate.get(key), run_dir=root, must_exist=True)
            if path and path.is_file():
                return path
    raise ValueError("The designed structure file is unavailable.")


def _filter_flag(value: Any) -> bool | None:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and value in (0, 1):
        return bool(value)
    if isinstance(value, str):
        return {"true": True, "1": True, "false": False, "0": False}.get(value.strip().lower())
    return None


def candidate_passes_filters(candidate: dict) -> bool:
    """Require affirmative overall filter evidence; unknown is not a pass."""
    metrics = candidate.get("metrics") or {}
    raw = candidate.get("raw_metadata") or {}
    if any(metrics.get(key) for key in ("ipsae_error", "rosetta_error", "validation_error", "qc_error")):
        return False
    if (raw.get("measurement_status") or {}).get("status") in {"error", "missing", "invalid"}:
        return False
    if raw.get("candidate_pool_level") == "prefilter_generated" or raw.get("result_kind") in {
        "generation_only", "generator_trajectory",
    }:
        return False
    flags = [_filter_flag(metrics.get(key)) for key in ("native_pass_filters", "pass_filters")]
    if False in flags:
        return False
    if True in flags:
        return True
    # Older BindCraft results identify acceptance by the native output pool.
    # A successful process or a ranked/generated structure alone is insufficient.
    return (candidate.get("source_tool") == "bindcraft"
            and raw.get("bindcraft_output_kind") == "Accepted")


def design_candidates(run_dir: Path, *, require_passing: bool = True) -> list[dict]:
    """Read designs from a successful design job, optionally requiring native pass.

    Ordinary cofolding selection remains restricted to affirmative native-pass
    evidence. The standalone acceptance tool sets ``require_passing=False`` so
    its independent prediction can evaluate rejected or unfiltered designs too.
    """
    metadata = read_json(run_dir / "metadata.json")
    result = read_json(run_dir / "result.json")
    if (metadata.get("task_group") != "design"
            or metadata.get("status") not in {"completed", "succeeded"}
            or result.get("success") is not True):
        return []
    candidates = read_candidates(run_dir)
    if require_passing:
        return [candidate for candidate in candidates if candidate_passes_filters(candidate)]
    eligible = []
    for candidate in candidates:
        raw = candidate.get("raw_metadata") or {}
        if raw.get("candidate_pool_level") == "prefilter_generated" or raw.get("result_kind") in {
            "generation_only", "generator_trajectory",
        }:
            continue
        if candidate.get("stage") in {"generation.backbone", "generation_backbone", "generator_trajectory"}:
            continue
        eligible.append(candidate)
    return eligible


def candidate_description(run_dir: Path, candidate: dict) -> dict:
    path = _candidate_structure(run_dir, candidate)
    try:
        chains = protein_chains(read_structure(path))
    except (OSError, ValueError, KeyError, StopIteration, PDBException) as exc:
        raise ValueError("The designed structure could not be read.") from exc
    binders = [str(c) for c in candidate.get("binder_chains") or []]
    targets = [str(c) for c in candidate.get("target_chains") or []]
    if not binders and len(chains) == 1:
        binders = list(chains)
    if not binders or any(c not in chains for c in binders):
        raise ValueError("This design has no identifiable protein binder chain. Assign its chain roles first.")
    if set(binders) & set(targets):
        raise ValueError("The design marks the same chain as target and binder.")
    if any(c not in chains for c in targets):
        raise ValueError("A target protein chain is missing from the designed structure.")
    if candidate.get("stage") in {"generation_backbone", "generation.backbone"} and not candidate.get("binder_sequence"):
        raise ValueError("This backbone needs sequence design before it can be refolded.")
    binder_sequence = ":".join(sequence(chains[c]) for c in binders)
    if len(binders) == 1 and candidate.get("binder_sequence"):
        binder_sequence = "".join(str(candidate["binder_sequence"]).split()).upper()
    return {"structure_path": path, "chains": chains, "binder_chains": binders,
            "target_chains": targets, "binder_sequence": binder_sequence}


def build_tasks(sources: list[dict], *, mode: str, profile: dict, target_use_msa: bool,
                use_target_template: bool, target_type: str = "structured",
                binder_format: str | None = None, require_passing: bool = True) -> list[dict]:
    if mode not in {"complex", "binder"}:
        raise ValueError("Choose complex refolding or binder-only prediction.")
    if target_type not in {"structured", "fasta_idr"}:
        raise ValueError("Choose a structured or FASTA / IDR target.")
    if target_type == "fasta_idr" and profile["input_kind"] == "boltzgen":
        raise ValueError("BoltzGen Fold validation uses target coordinates as structural input; use a sequence-only-capable co-folding engine for FASTA / IDR targets.")
    if not isinstance(sources, list) or not sources or len(sources) > 128:
        raise ValueError("Select between 1 and 128 design runs.")
    if mode == "binder" and profile["input_kind"] == "boltzgen":
        raise ValueError("BoltzGen Fold requires a target; choose another engine for binder-only prediction.")
    from mn_molecular_design_studio.validation_acceptance import infer_binder_format
    tasks, seen = [], set()
    for source in sources:
        if not isinstance(source, dict):
            raise ValueError("Each selected design run must specify its run directory and design IDs.")
        run_dir = resolve_stored_path(source.get("source_run_dir"), must_exist=True)
        if not run_dir or not run_dir.is_dir():
            raise ValueError("A selected design run is unavailable.")
        rows = {c["candidate_id"]: c for c in design_candidates(run_dir, require_passing=require_passing)}
        try:
            workflow_parameters = read_json(run_dir / "input.json").get("params", {})
        except (OSError, ValueError, TypeError):
            workflow_parameters = {}
        native_settings = run_dir / "artifacts/raw/bindcraft2/settings.json"
        if native_settings.is_file():
            try:
                workflow_parameters = {"input_params": workflow_parameters,
                                       "bindcraft2_settings": read_json(native_settings)}
            except (OSError, ValueError, TypeError):
                pass
        ids = source.get("candidate_ids")
        if ids is None:
            limit = source.get("candidate_limit")
            if type(limit) is not int or not 1 <= limit <= 1024:
                raise ValueError("A remote design source needs a candidate_limit between 1 and 1024.")
            ids = list(rows)[:limit]
        if not isinstance(ids, list) or not ids:
            raise ValueError("Select at least one design from every selected run.")
        if source.get("candidate_limit") is not None and len(ids) < source["candidate_limit"]:
            raise ValueError(
                f"Only {len(ids)} eligible design(s) are available in {run_dir.name}; "
                f"the requested batch contains {source['candidate_limit']}."
            )
        for identifier in ids:
            key = (str(run_dir), str(identifier))
            if key in seen:
                continue
            seen.add(key)
            if identifier not in rows:
                eligibility = "filter-passing output" if require_passing else "eligible output from a successful design job"
                raise ValueError(f"Design {identifier!r} is not an {eligibility} in {run_dir.name}.")
            candidate = rows[identifier]
            inferred_format = infer_binder_format(candidate, workflow_parameters)
            selected_format = binder_format or inferred_format
            cyclic = selected_format == "cyclic_peptide" or inferred_format == "cyclic_peptide" or bool(
                (candidate.get("metrics") or {}).get("cyclic")
                or (candidate.get("raw_metadata") or {}).get("cyclic")
            )
            if cyclic and not profile.get("supports_cyclic"):
                raise ValueError(f"The selected co-folding engine {profile['label']} does not preserve cyclic peptide input topology; choose an engine with cyclic-chain support.")
            try:
                info = candidate_description(run_dir, candidate)
            except ValueError as exc:
                raise ValueError(f"Design {identifier}: {exc}") from exc
            target_ids = info["target_chains"] if mode == "complex" else []
            if mode == "complex" and not target_ids:
                raise ValueError(f"Design {identifier} has no protein target for complex alignment; use binder-only prediction.")
            ordered = [*target_ids, *info["binder_chains"]]
            if len(ordered) > len(CHAIN_IDS):
                raise ValueError("The selected complex has too many protein chains.")
            chain_map = {chain: CHAIN_IDS[i] for i, chain in enumerate(ordered)}
            options = [{"id": chain_map[c], "use_msa": bool(target_use_msa and c in target_ids),
                        **({"cyclic": True} if cyclic and c in info["binder_chains"] else {})}
                       for c in ordered]
            templates = []
            if target_type == "structured" and use_target_template and target_ids and profile.get("supports_templates"):
                for c in target_ids:
                    class TargetChain(Select):
                        def accept_chain(self, chain):
                            return chain.id == c
                    stream = io.StringIO()
                    writer = MMCIFIO(); writer.set_structure(read_structure(info["structure_path"]))
                    writer.save(stream, TargetChain())
                    templates.append({"name": f"target_{chain_map[c]}.cif", "text": stream.getvalue(),
                                      "query_chain_id": chain_map[c], "template_chain_id": c})
            task = {
                "name": f"{run_dir.name}/{identifier}",
                "source_candidate_id": identifier,
                "target_protein_sequences": ":".join(sequence(info["chains"][c]) for c in target_ids),
                "binder_sequence": info["binder_sequence"],
                "molecule_chain_options": options,
                "structure_templates": templates,
                "validation_reference": {
                    "source_run_dir": str(run_dir), "candidate_id": identifier,
                    "reference_structure_path": str(info["structure_path"]), "mode": mode,
                    "target_type": target_type,
                    "binder_format": selected_format, "cyclic": cyclic,
                    "chain_map": chain_map, "target_chains": target_ids,
                    "binder_chains": info["binder_chains"],
                },
            }
            if profile["input_kind"] == "boltzgen":
                if len(target_ids) != 1 or len(info["binder_chains"]) != 1:
                    raise ValueError("BoltzGen Fold currently accepts one target protein and one binder chain.")
                from Bio.PDB import PDBIO
                stream = io.StringIO()
                writer = PDBIO(); writer.set_structure(read_structure(info["structure_path"]))
                class Target(Select):
                    def accept_chain(self, chain):
                        return chain.id == target_ids[0]
                writer.save(stream, Target())
                task.update(target_pdb_text=stream.getvalue(), target_chain_id=target_ids[0],
                            target_protein_sequences="", molecule_chain_options=[], structure_templates=[])
                if not target_use_msa:
                    task["target_msa_a3m"] = f'>query\n{sequence(info["chains"][target_ids[0]])}\n'
                task["validation_reference"]["chain_map"] = {target_ids[0]: target_ids[0], info["binder_chains"][0]: "Z"}
                task["binder_use_msa"] = False
            tasks.append(task)
    if not tasks or len(tasks) > 1024:
        raise ValueError("Select between 1 and 1024 designs per validation batch.")
    return tasks


def stage_reference(run_dir: Path, item_id: str, reference: dict) -> dict:
    path = resolve_stored_path(reference.get("reference_structure_path"), must_exist=True)
    if not path or not path.is_file():
        raise ValueError("The designed reference structure is unavailable.")
    extension = Path(path.stem).suffix if path.suffix == ".gz" else path.suffix
    destination = run_dir / "artifacts" / "validation" / item_id / ("designed" + extension)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes())
    return {**reference, "reference_structure_path": str(destination)}


def compare_prediction(run_dir: Path, item_id: str, sample_index: int, prediction: Path,
                       reference: dict) -> tuple[dict, dict]:
    """Save a fitted structure and return explicit comparison and coverage metrics."""
    original = read_structure(Path(reference["reference_structure_path"]))
    predicted = read_structure(prediction)
    ref_chains, pred_chains = protein_chains(original), protein_chains(predicted)
    mapped, used = {}, set()
    for source, expected in reference["chain_map"].items():
        ref_seq = sequence(ref_chains[source])
        matches = [c for c in pred_chains if c not in used and sequence(pred_chains[c]) == ref_seq]
        if expected in matches:
            mapped[source] = expected
        elif len(matches) == 1:
            mapped[source] = matches[0]
        elif expected in pred_chains and expected not in used:
            predicted_seq = sequence(pred_chains[expected])
            pairs = residue_pairs(ref_seq, predicted_seq)
            identity = sum(ref_seq[i] == predicted_seq[j] for i, j in pairs) / max(1, len(pairs))
            if identity < .5:
                raise ValueError(f"Predicted chain {expected} does not match design chain {source}.")
            mapped[source] = expected
        else:
            raise ValueError(f"Cannot unambiguously match predicted chain for design chain {source}.")
        used.add(mapped[source])

    def coordinates(ids):
        ref, pred, total = [], [], 0
        for chain in ids:
            a, b = ref_chains[chain], pred_chains[mapped[chain]]
            total += len(a)
            for i, j in residue_pairs(sequence(a), sequence(b)):
                ref.append(a[i]["CA"].coord); pred.append(b[j]["CA"].coord)
        return np.asarray(ref, dtype=float), np.asarray(pred, dtype=float), total

    binder_ref, binder_pred, binder_total = coordinates(reference["binder_chains"])
    if len(binder_ref) < 3:
        raise ValueError("At least three matched binder C-alpha atoms are needed for RMSD.")
    complex_mode = reference["mode"] == "complex"
    target_pose_comparable = complex_mode and reference.get("target_type", "structured") != "fasta_idr"
    fit_ref, fit_pred, fit_total = coordinates(reference["target_chains"]) if target_pose_comparable else (binder_ref, binder_pred, binder_total)
    if not all(np.isfinite(array).all() for array in (fit_ref, fit_pred, binder_ref, binder_pred)):
        raise ValueError("Non-finite coordinates prevent structure comparison.")
    if len(fit_ref) < 3 or any(np.linalg.matrix_rank(array - array.mean(axis=0)) < 2 for array in (fit_ref, fit_pred)):
        raise ValueError("At least three non-collinear matched atoms are needed for alignment.")
    ref_center, pred_center = fit_ref.mean(axis=0), fit_pred.mean(axis=0)
    u, _, vt = np.linalg.svd((fit_pred - pred_center).T @ (fit_ref - ref_center))
    correction = np.eye(3); correction[2, 2] = np.sign(np.linalg.det(u @ vt))
    rotation = u @ correction @ vt
    translation = ref_center - pred_center @ rotation
    rmsd = lambda a, b: float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))
    binder_rmsd = rmsd(binder_ref, binder_pred @ rotation + translation)
    metrics = {
        "validation_mode": reference["mode"], "validation_status": "completed",
        "binder_rmsd_target_aligned" if target_pose_comparable else "binder_rmsd_fitted": binder_rmsd,
        "binder_rmsd_ca_count": len(binder_ref), "binder_rmsd_coverage": len(binder_ref) / binder_total,
    }
    if target_pose_comparable:
        metrics.update(target_alignment_rmsd=rmsd(fit_ref, fit_pred @ rotation + translation),
                       target_alignment_ca_count=len(fit_ref), target_alignment_coverage=len(fit_ref) / fit_total)
    elif complex_mode:
        metrics.update(binder_self_consistency_rmsd=binder_rmsd,
                       target_pose_comparison_status="unavailable_unconstrained_target")
    else:
        metrics.update(monomer_refolding_success=True, monomer_refolding_rmsd=metrics["binder_rmsd_fitted"])
    if reference.get("cyclic") and len(reference.get("binder_chains") or []) == 1:
        predicted_binder = pred_chains[mapped[reference["binder_chains"][0]]]
        ordered_residues = list(predicted_binder)
        if ordered_residues and "N" in ordered_residues[0] and "C" in ordered_residues[-1]:
            metrics["cyclic_closure_distance"] = float(np.linalg.norm(
                ordered_residues[0]["N"].coord - ordered_residues[-1]["C"].coord
            ))
            metrics["cyclic_closure_status"] = "measured"
        else:
            metrics["cyclic_closure_status"] = "unavailable_terminal_atoms"
    fitted = copy.deepcopy(predicted)
    for atom in fitted.get_atoms():
        atom.coord = atom.coord @ rotation + translation
    destination = run_dir / "artifacts" / "validation" / item_id / f"sample_{sample_index}_aligned.cif"
    destination.parent.mkdir(parents=True, exist_ok=True)
    writer = MMCIFIO(); writer.set_structure(fitted); writer.save(str(destination))
    return metrics, {"validation_reference": reference, "aligned_structure_path": str(destination),
                     "reference_structure_path": reference["reference_structure_path"],
                     "validation_chain_map": mapped}


def write_report(run_dir: Path, candidates: list[dict]) -> Path:
    rows = []
    for candidate in candidates:
        reference = (candidate.get("raw_metadata") or {}).get("validation_reference") or {}
        if reference:
            rows.append({"candidate_id": candidate["candidate_id"], "source_run_dir": reference["source_run_dir"],
                         "source_candidate_id": reference["candidate_id"],
                         **{k: v for k, v in candidate.get("metrics", {}).items()
                            if k.startswith(("binder_rmsd", "target_alignment", "validation_", "cyclic_closure"))}})
    path = run_dir / "artifacts/validation/comparisons.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)
    return path
