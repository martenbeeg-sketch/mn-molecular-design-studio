"""Offline bridge to the pinned BinderScope engine, executed inside its image."""
from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import math
from pathlib import Path


REVISION = "c61ef9ca3dfc64dbde018ae7329d7538d8c583f1"


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read_structure(path):
    import gemmi

    if str(path).lower().endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz")):
        document = gemmi.cif.read(str(path))
        blocks = [block for block in document if len(block.find_values("_atom_site.id"))]
        if len(blocks) != 1:
            raise ValueError("BinderScope requires one coordinate-bearing mmCIF block per saved representation.")
        atom_site = blocks[0].get_mmcif_category("_atom_site.")
        atom_count = len(atom_site["id"])
        atom_site.setdefault("occupancy", ["1.0"] * atom_count)
        atom_site.setdefault("B_iso_or_equiv", ["nan"] * atom_count)
        blocks[0].set_mmcif_category("_atom_site.", atom_site)
        structure = gemmi.make_structure_from_block(blocks[0])
        confidence = dict(zip(atom_site["id"], atom_site["B_iso_or_equiv"]))
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        try:
                            value = confidence.get(str(atom.serial))
                            atom.b_iso = float(value) if isinstance(value, str) else float("nan")
                        except (KeyError, TypeError, ValueError):
                            atom.b_iso = float("nan")
    else:
        structure = gemmi.read_structure(str(path))
    return structure


def normalize_complex(path, destination, binder_chain, target_chain):
    import gemmi

    structure = read_structure(path)
    if len(structure) != 1:
        raise ValueError("BinderScope requires one exact model per saved representation.")
    structure.remove_alternative_conformations()
    model = structure[0]
    names = [chain.name for chain in model]
    if names.count(binder_chain) != 1 or names.count(target_chain) != 1 or binder_chain == target_chain:
        raise ValueError("Explicit unique binder and target chains are required.")
    result = gemmi.Structure()
    result.add_model(gemmi.Model("1"))
    mapping, sequences, addresses = {}, {}, {}
    for original, canonical in ((target_chain, "A"), (binder_chain, "B")):
        copied = model[original].clone()
        copied.name = canonical
        sequence, residues = [], []
        for residue in copied:
            info = gemmi.find_tabulated_residue(residue.name)
            if not info.is_amino_acid():
                raise ValueError("Selected chains contain non-protein residues; choose a protein-only design complex.")
            if not any(atom.name == "CA" for atom in residue):
                raise ValueError("Selected protein residues require C-alpha coordinates.")
            number = len(sequence) + 1
            author = {"chain": original, "number": residue.seqid.num,
                      "insertion_code": residue.seqid.icode.strip()}
            addresses[(original, residue.seqid.num, residue.seqid.icode)] = (canonical, number)
            residues.append({"canonical_chain": canonical, "canonical_number": number, **author})
            residue.seqid = gemmi.SeqId(number, " ")
            sequence.append(info.one_letter_code.upper())
        if not sequence:
            raise ValueError("Selected protein chain is empty.")
        result[0].add_chain(copied)
        mapping[canonical] = residues
        sequences[canonical] = "".join(sequence)
    for connection in structure.connections:
        copied = connection
        partners = [copied.partner1, copied.partner2]
        keys = [(partner.chain_name, partner.res_id.seqid.num, partner.res_id.seqid.icode) for partner in partners]
        if all(key in addresses for key in keys):
            for partner, key in zip(partners, keys):
                partner.chain_name, number = addresses[key]
                partner.res_id.seqid = gemmi.SeqId(number, " ")
            result.connections.append(copied)
    destination.parent.mkdir(parents=True, exist_ok=True)
    result.write_pdb(str(destination))
    return mapping, sequences


def saved_projects(request, work):
    from binder_review import model, pipeline
    import gemmi
    if __package__:
        from .binderscope_context import prepare_context
    else:
        from binderscope_context import prepare_context

    grouped, sources = {}, {}
    for item in request["items"]:
        item = dict(item)
        assessment = dict(item.get("assessment") or {
            "component_format": request.get("component_format", "unspecified"),
            "b_factor_interpretation": "plddt_0_100" if request.get("bfactors_are_plddt") else "uninterpreted",
        })
        if (assessment.get("component_format") not in {"unspecified", "linear protein", "linear peptide", "cyclic peptide"}
                or assessment.get("b_factor_interpretation") not in {"uninterpreted", "plddt_0_100"}):
            raise ValueError("Invalid per-design BinderScope settings.")
        source = Path(item["path"]).resolve(strict=True)
        source.relative_to(Path(request["input_root"]).resolve(strict=True))
        if digest(source) != item["input_structure_digest"]:
            raise ValueError("A staged design changed after authorization.")
        identifier = item["candidate_id"]
        canonical = work / "normalized" / (identifier + ".pdb")
        residue_map, sequences = normalize_complex(source, canonical, item["binder_chains"][0], item["target_chains"][0])
        if assessment["b_factor_interpretation"] == "plddt_0_100":
            normalized = gemmi.read_structure(str(canonical))
            atoms = [atom for residue in normalized[0]["B"]
                     for atom in residue if atom.name == "CA"]
            if any(not math.isfinite(atom.b_iso) or not 0 <= atom.b_iso <= 100 for atom in atoms):
                assessment["b_factor_interpretation"] = "uninterpreted"
                assessment["confidence_basis"] = "Declared pLDDT failed the 0–100 range check; confidence checks disabled."
        item["assessment"] = assessment
        reference = item.get("target_reference")
        if reference:
            reference_path = Path(reference["path"]).resolve(strict=True)
            reference_path.relative_to(Path(request["input_root"]).resolve(strict=True))
            if digest(reference_path) != reference["input_structure_digest"]:
                raise ValueError("Original target reference changed after authorization.")
        context = item.get("source", {}).get("target_context") or {}
        prepared_context = prepare_context(context, reference, residue_map, sequences, canonical, read_structure)
        item["target_assessment"] = prepared_context
        cohort = hashlib.sha256(sequences["A"].encode()).hexdigest()[:20]
        context_key = hashlib.sha256(json.dumps([context, [row.get("canonical_number") for row in prepared_context["rows"]]], sort_keys=True).encode()).hexdigest()[:16]
        run_id = "saved-" + assessment["component_format"].replace(" ", "-") + "-" + context_key
        design = model.DesignRecord(
            id=identifier, run_id=run_id, backbone_id=identifier, status="accepted",
            binder_name=item["label"], length=len(sequences["B"]), seed=0, mpnn_index=None,
            sequence=sequences["B"], structure_path=str(canonical),
            extras={"studio_source": item, "saved_cohort_only": True},
        )
        grouped.setdefault(cohort, []).append((design, sequences["A"], assessment))
        sources[identifier] = {**item, "residue_map": residue_map, "cohort": cohort,
                               "normalized_sha256": digest(canonical)}
    for cohort, designs in grouped.items():
        target_path = work / cohort / "target.pdb"
        target_path.parent.mkdir(parents=True, exist_ok=True)
        target = gemmi.read_structure(designs[0][0].structure_path)
        del target[0]["B"]
        target.write_pdb(str(target_path))
        runs = {}
        for design, sequence, assessment in designs:
            if design.run_id not in runs:
                settings = model.RunSettings(run_id=design.run_id, label="Admitted designs",
                    hotspots=sorted({row["canonical_number"] for row in design.extras["studio_source"]["target_assessment"]["rows"] if row["status"] == "mapped"}),
                    settings_raw={"modality": assessment["component_format"],
                                  "cyclize_peptide": assessment["component_format"] == "cyclic peptide"})
                target_construct = model.TargetConstruct(run_id=design.run_id, numbers=list(range(1, len(sequence) + 1)),
                    sequence=sequence, first=1, last=len(sequence), consistent=True)
                runs[design.run_id] = model.RunData(run_id=design.run_id, label="Admitted designs", reader="studio_saved",
                    source_dir=str(work), settings=settings, target=target_construct, designs=[])
                runs[design.run_id].warnings.append("Saved cohort only: accepted means admitted to Studio, not a newly inferred filter pass. "
                    "Original campaign backbones, rejected attempts, funnel and stopping estimates are unavailable.")
            runs[design.run_id].designs.append(design)
        spec = pipeline.ProjectSpec(name="Saved design assessment", runs=[], target_pdb=target_path,
                                    profile=request.get("profile") or "default", jobs=request["cpu_cores"], online=False)
        yield cohort, spec, runs, sources


def analyze_saved(spec, runs, work):
    from binder_review import pipeline

    log = pipeline.Log()
    profile = pipeline.load_profile(spec.profile)
    spec = dataclasses.replace(spec, uniprot=pipeline.effective_uniprot(spec, profile))
    digests = pipeline.input_digests(spec, runs)
    frame = pipeline.build_target_frame(spec, runs, digests, work, log)
    frame_key = pipeline.frame_key(spec, runs, digests)
    poses = pipeline.compute_poses(runs, frame, frame_key, work, jobs=spec.jobs, use_cache=False, log=log)
    analyses = pipeline.compute_interface(runs, work, jobs=spec.jobs, use_cache=False, log=log)
    for (run_id, identifier), analysis in analyses.items():
        design = next(design for design in runs[run_id].designs if design.id == identifier)
        analysis.plddt_t = {}
        if design.extras["studio_source"]["assessment"]["b_factor_interpretation"] != "plddt_0_100":
            analysis.plddt_res = []
            analysis.tracks.pop("plddt", None)
    families = pipeline.compute_families(runs, digests, work, jobs=spec.jobs, use_cache=False, log=log)
    identifiers = pipeline.design_ids(runs)
    modes, rescue = pipeline.compute_modes(runs, frame, poses, identifiers, digests, work, use_cache=False, log=log)
    homolog = pipeline.compute_specificity(runs, profile, work, accession=spec.uniprot,
                                          frame=frame, use_cache=False, log=log)
    return pipeline.Project(spec=spec, work=work, runs=runs, frame=frame, profile=profile,
                            poses=poses, analyses=analyses, families=families, modes=modes,
                            rescue=rescue, homolog=homolog, project_ids=identifiers, digests=digests,
                            warnings=[warning for run in runs.values() for warning in run.warnings])


def run(request_path):
    import binder_review
    from binder_review import bundle, export, pipeline
    import importlib.metadata

    source_root = Path(binder_review.__file__).resolve().parents[3]
    source_manifest = source_root / "SOURCE_MANIFEST.json"
    recorded_source = json.loads(source_manifest.read_text())
    if recorded_source["upstream_revision"] != REVISION:
        raise ValueError("BinderScope runtime source revision does not match the Studio contract.")
    for relative, expected_digest in recorded_source["files"].items():
        if digest(source_root / relative) != expected_digest:
            raise ValueError("Integrated BinderScope source checksum mismatch: " + relative)
    request = json.loads(request_path.read_text())
    if request.get("adapter_sha256") and request["adapter_sha256"] != digest(Path(__file__)):
        raise ValueError("BinderScope runtime adapter is out of date. Rebuild its image from this repository.")
    if request.get("context_sha256") and request["context_sha256"] != digest(Path(__file__).with_name("binderscope_context.py")):
        raise ValueError("BinderScope target-context adapter is out of date. Rebuild its image.")
    if request.get("mode") not in {"saved", "native"}:
        raise ValueError("Choose a saved SQL cohort or native campaign input mode.")
    if request["mode"] == "saved" and not request.get("use_recorded_settings") and request.get("bfactors_are_plddt") is not True:
        raise ValueError("Saved BinderScope analysis requires explicit confirmation of the pLDDT encoding.")
    if request["mode"] == "saved" and request.get("use_recorded_settings") and any("assessment" not in item for item in request["items"]):
        raise ValueError("Every selected design requires its own assessment settings.")
    output = Path(request["output_root"])
    output.mkdir(parents=True, exist_ok=True)
    if request["mode"] == "native":
        partners = request.get("partner_chains") or "none"
        if partners not in {"auto", "none", "all"}:
            partners = [chain.strip() for chain in partners.split(",") if chain.strip()]
        spec = pipeline.ProjectSpec(
            name="BinderScope campaign assessment",
            runs=[pipeline.RunSpec(path=Path(path), run_id=f"run{index + 1}")
                  for index, path in enumerate(request["native_run_dirs"])],
            target_pdb=Path(request["target_pdb"]),
            reference_cif=Path(request["reference_cif"]) if request.get("reference_cif") else None,
            profile=request.get("profile") or "default", uniprot=request.get("uniprot") or None,
            reference_chain=request.get("reference_chain") or None, partner_chains=partners,
            jobs=request["cpu_cores"], online=False, consurf=False,
        )
        projects = [("campaign", pipeline.run_pipeline(spec, output / "campaign"), {})]
    else:
        projects = [(cohort, analyze_saved(spec, runs, output / cohort), sources)
                    for cohort, spec, runs, sources in saved_projects(request, output)]
    summaries, results = [], []
    for cohort, project, sources in projects:
        folder = output / cohort
        if request["mode"] == "native":
            blobs = export.build_blobs(project, validate=True, strict=True)
            export.write_export(blobs, folder / "export")
            bundle.write_bundle(folder / "export", folder / "review", name=project.spec.name,
                                fonts="system", embed_structures=True, guide=False)
        paratopes = export.run_paratopes(project)
        for key, analysis in project.analyses.items():
            source = sources.get(key[1], {})
            raw = analysis.to_dict()
            settings = source.get("assessment") or {}
            confident = request["mode"] == "native" or settings.get("b_factor_interpretation") == "plddt_0_100"
            topology = settings.get("component_format", "unspecified")
            if not confident:
                for name in ("plddt_res", "plddt_t"):
                    raw.pop(name, None)
                raw.get("tracks", {}).pop("plddt", None)
                raw["confidence_status"] = "Unavailable: B-factor encoding was not confirmed as pLDDT"
            if request["mode"] == "saved":
                raw.pop("plddt_t", None)
                raw["target_confidence_status"] = "Target B-factors remain uninterpreted; confidence declarations concern the selected binder."
            design = next(design for design in project.runs[key[0]].designs if design.id == key[1])
            project_id = project.project_ids[key]
            termini = export._terminus_block(project, key[0], design, analysis, project.poses[key],
                export.termini.format_of(project.runs[key[0]])) if request["mode"] == "native" or (
                    topology != "unspecified" and (confident or topology == "cyclic peptide")) else {
                        "status": "unavailable", "reason": "Tag/fusion assessment requires recorded topology and confirmed per-residue confidence. "
                        "Geometry, contacts, patches, liabilities and cohort comparisons are still evaluated."}
            result = {"candidate_id": key[1], "cohort": cohort, "source": source, "analysis": raw,
                      "pose": project.poses[key], "binding_modes": {
                          "assignment": project.modes.get("assign", {}).get(project_id),
                          "hotspots_touched": project.modes.get("hotspots_touched", {}).get(project_id),
                          "cohort": cohort},
                      "termini": termini,
                      "specificity": {"cohort": cohort, "status": "See cohort panel for supplied reference evidence"},
                      "paratope": {"cohort": cohort, "status": "See native campaign scaffold panel" if paratopes[key[0]][0] else "Unavailable: no native scaffold edit specification"}}
            if __package__:
                from .binderscope_context import hotspot_result
            else:
                from binderscope_context import hotspot_result
            if request["mode"] == "saved":
                result["hotspot_compliance"] = hotspot_result(source["target_assessment"], analysis)
                result["target_reference_comparison"] = source["target_assessment"]["target_reference_comparison"]
            else:
                requested_hotspots = project.runs[key[0]].settings.hotspots
                result["hotspot_compliance"] = hotspot_result({
                    "declaration": {"status": "recorded" if requested_hotspots else "not_specified",
                                    "evidence": [{"source": key[0], "field": "hotspots", "value": requested_hotspots}]},
                    "rows": [{"requested": str(number), "canonical_number": number,
                              "saved_residue": f"{number} (construct)",
                              "status": "mapped" if number in analysis.hotspots else "unmapped"} for number in requested_hotspots],
                    "interpretation": "Original native run hotspot numbering; computational contacts do not change admission."}, analysis)
            results.append(result)
        families = project.families if request["mode"] == "native" else {
            key: project.families[key] for key in ("fold", "fold_labels", "design_ids", "warnings")}
        summary = project.summary() if request["mode"] == "native" else {
            "saved_designs": len(project.analyses), "scope": "Admitted SQL subset; no campaign yield or stopping inference"}
        summaries.append({"cohort": cohort, "summary": summary, "warnings": project.warnings,
                          "families": families, "modes": project.modes, "specificity": project.homolog})
    report = {"schema": "studio.binderscope_report.v1", "upstream_revision": REVISION,
              "source_manifest_sha256": digest(source_manifest),
              "adapter_sha256": digest(Path(__file__)),
              "target_context_sha256": digest(Path(__file__).with_name("binderscope_context.py")),
              "mode": request["mode"], "cohorts": summaries, "results": results,
              "dependencies": {name: importlib.metadata.version(name) for name in
                               ("numpy", "scipy", "gemmi", "biopython", "scikit-learn", "tmtools")},
              "limitations": ["Computational evidence is not experimental binding or developability.",
                  "Saved SQL cohorts do not contain rejected attempts; campaign yield and near-miss statistics require native run folders.",
                  "Reference-context and specificity checks require their explicit input data; offline runtime does not fetch public records."]}
    (output / "binderscope_report.json").write_text(json.dumps(report, default=pipeline._json_default, allow_nan=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("request", type=Path)
    run(parser.parse_args().request)
