"""Molecular Database browser and import wizard used by Studio."""
from __future__ import annotations

import hashlib
import html
import json
import math
import sqlite3
import uuid
import uuid
from pathlib import Path

import streamlit as st

from mn_molecular_design_studio import molecule_imports as imports
from mn_molecular_design_studio.molecular_database import MolecularDatabase, IMPORT_METHODS
from mn_molecular_design_studio.molecular_database_ingest import sync_workspace
from mn_molecular_design_studio.output_vocabulary import PREFIX, field_label

TYPE_LABELS = {"protein": "Protein", "small_molecule": "Small molecule", "dna": "DNA", "rna": "RNA", "complex": "Complex"}
CSS = '''<style>
.st-key-molecular-database h3{font-size:1.05rem;margin-bottom:.35rem}.st-key-molecular-database h4{font-size:.95rem}
.st-key-molecular-database [data-testid="stButton"] button{min-height:30px;font-size:12px;border-radius:4px}
.st-key-molecular-database [data-testid="stTextInput"] input{font-size:12px}.st-key-molecular-database [data-testid="stCaptionContainer"]{font-size:11px;color:#869087}
.st-key-molecular-database [data-testid="stVerticalBlock"]{gap:.55rem}.st-key-database-groups{padding-right:12px;border-right:1px solid #e0e6e1}
.st-key-database-details{padding-left:12px;border-left:1px solid #e0e6e1}.st-key-database-groups [role="radiogroup"]{gap:0}
.st-key-database-groups [data-testid="stRadio"] label{padding:6px 7px;border-radius:4px;font-size:12px}
.st-key-database-groups [data-testid="stRadio"] label:has(input:checked){background:#e8f0e9;color:#246444}
.database-empty{border:1px dashed #c9d8cd;border-radius:6px;padding:35px 20px;text-align:center;color:#859188;font-size:13px}
.database-label{color:#8b958d;font-size:10px;letter-spacing:.05em;text-transform:uppercase;margin:12px 0 5px}
</style>'''


def _projects(projects, key, *, default="", label="Project", optional=True):
    options = ([""] if optional else []) + [row["id"] for row in projects]
    if not options:
        st.caption("Create a project in Projects before starting a calculation.")
        return ""
    if st.session_state.get(key) not in options:
        st.session_state[key] = default if default in options else options[0]
    return st.selectbox(label, options, key=key, format_func=lambda value: next((row["name"] for row in projects if row["id"] == value), "None"))


def _reset_import():
    st.session_state["molecular-import-open"] = False
    st.session_state["molecular-import-step"] = 0
    st.rerun()


def _guess(columns, names):
    return next((column for column in columns if column.lower().replace(" ", "_") in names), "")


def _import_content(method):
    records = []
    audit = None
    source_content = None
    if method == "Table":
        st.caption("One row per molecule. Map protein, ligand, DNA, or RNA columns to chains; other columns become molecular properties.")
        upload = st.file_uploader("CSV, TSV, or Excel file", type=["csv", "tsv", "xlsx", "xlsm", "xls", "txt"], key="molecular-import-table")
        if upload is None:
            return [], None, None
        upload_content = upload.getvalue()
        signature = hashlib.sha256(upload_content).hexdigest()[:12]
        sheets = imports.table_sheets(upload_content, upload.name)
        sheet_name = sheets[0] if sheets else ""
        if len(sheets) > 1:
            sheet_name = st.selectbox("Worksheet", list(sheets), key="molecular-import-sheet-" + signature)
        table = imports.read_table(upload_content, upload.name, sheet_name=sheet_name)
        columns = list(table.columns)
        id_guess = _guess(columns, {"compound_id", "pubchem_cid", "catalog_id", "catalog_number", "catalog_no", "vendor_id", "zinc_id", "candidate_id", "design_id", "id"})
        name_guess = _guess(columns, {"compound_name", "preferred_name", "common_name", "molecule_name", "name", "molecule"})
        identity_columns = st.columns(2)
        id_options = ["", *columns]
        name_options = ["", *columns]
        id_column = identity_columns[0].selectbox(
            "Molecule ID / compound ID column",
            id_options,
            index=id_options.index(id_guess) if id_guess in id_options else 0,
            format_func=lambda value: value or "No separate ID",
            key="molecular-import-id-column",
        )
        name_column = identity_columns[1].selectbox(
            "Display name column",
            name_options,
            index=name_options.index(name_guess) if name_guess in name_options else 0,
            format_func=lambda value: value or "Use ID as name, or generate a name",
            key="molecular-import-name-column",
        )
        defaults = [column for column in columns if column.lower().replace(" ", "_") in {"sequence", "protein_sequence", "binder_sequence", "target_sequence", "smiles", "canonical_smiles", "isomeric_smiles", "dna", "rna", "heavy", "light", "vh", "vl"}]
        chain_columns = st.multiselect("Molecular chain columns", columns, default=defaults, key="molecular-import-chains-" + signature)
        mapping = {}
        for index, column in enumerate(chain_columns):
            kind_guess = "small_molecule" if "smiles" in column.lower() else "dna" if "dna" in column.lower() else "rna" if "rna" in column.lower() else "protein"
            cols = st.columns([3, 2, 1, 2])
            cols[0].markdown(f"**{html.escape(column)}**")
            kind = cols[1].selectbox("Type", imports.KINDS, index=imports.KINDS.index(kind_guess), format_func=TYPE_LABELS.get, key=f"import-type-{signature}-{column}")
            chain = cols[2].text_input("Chain", value=chr(65+index) if index < 26 else f"C{index+1}", key=f"import-chain-{signature}-{column}")
            role = cols[3].selectbox("Role", ["", "lead", "target", "heavy", "light", "cofactor"], key=f"import-role-{signature}-{column}")
            mapping[column] = {"type": kind, "chain": chain, "role": role}
        st.dataframe(table.head(5), hide_index=True, use_container_width=True)
        if mapping:
            if len(mapping) == 1 and next(iter(mapping.values())).get("type") == "small_molecule":
                smiles_column = next(iter(mapping))
                analysis = imports.analyze_small_molecule_table(
                    upload_content, upload.name, smiles_column=smiles_column,
                    name_column=name_column, id_column=id_column, sheet_name=sheet_name,
                )
                records = analysis.pop("records")
                audit = analysis
                source_content = upload_content
            else:
                records = imports.parse_table(
                    upload_content, upload.name, name_column=name_column,
                    id_column=id_column, chain_columns=mapping, sheet_name=sheet_name,
                )
    elif method in {"FASTA", "Paste sequences"}:
        kind = st.selectbox("Sequence type", ["protein", "dna", "rna"], format_func=TYPE_LABELS.get, key="molecular-import-sequence-type")
        if method == "FASTA":
            upload = st.file_uploader("FASTA file", type=["fasta", "fa", "faa", "fna", "txt"], key="molecular-import-fasta")
            text = upload.getvalue().decode("utf-8-sig") if upload else ""
        else:
            text = st.text_area("Sequences", placeholder=">molecule_1\nM...\n>complex_2\nM...:M...", height=180, key="molecular-import-sequences")
        st.caption("Each FASTA record or plain-text line becomes a molecule. Separate chains of one molecule with colons.")
        if text.strip():
            records = imports.parse_sequences(text, kind=kind)
    elif method == "Structure files":
        uploads = st.file_uploader("PDB, mmCIF, SDF, MOL, or ZIP", type=["pdb", "ent", "cif", "mmcif", "sdf", "mol", "zip", "gz"], accept_multiple_files=True, key="molecular-import-structures")
        include_ligands = st.checkbox("Include bound ligands", value=True, key="molecular-import-include-ligands")
        for upload in uploads:
            records.extend(imports.parse_archive(upload.getvalue()) if upload.name.lower().endswith(".zip") else imports.parse_structure(upload.getvalue(), upload.name, include_ligands=include_ligands))
        if len(records) == 1 and len(records[0]["components"]) > 1:
            components = records[0]["components"]
            choices = [row["chain"] for row in components]
            chains = st.multiselect("Chains to import", choices, default=choices, key="molecular-import-structure-chains", format_func=lambda chain: next(f"{chain} · {TYPE_LABELS[row['type']]} · {len(row.get('sequence') or '')} residues" if row.get("sequence") else f"{chain} · {row.get('ccd') or 'Ligand'}" for row in components if row["chain"] == chain))
            records[0]["components"] = [row for row in components if row["chain"] in chains]
            # A chain subset is materialized below, so the imported coordinates
            # must agree with its chain composition.
            if set(chains) != set(choices):
                records[0] = cropped_record(records[0], chains)
    elif method == "SMILES":
        source = st.radio("Input", ["Paste SMILES", "SMI file", "SDF file"], horizontal=True, key="molecular-import-smiles-source")
        if source == "Paste SMILES":
            text = st.text_area("SMILES", placeholder="CCO ethanol\nCC(=O)Oc1ccccc1C(=O)O aspirin", key="molecular-import-smiles", height=150)
            if text.strip():
                analysis = imports.analyze_smiles_text(text)
                records = analysis.pop("records")
                source_content = text.encode("utf-8")
                audit = analysis
                audit["source_name"] = "pasted-smiles.smi"
                audit["source_sha256"] = hashlib.sha256(source_content).hexdigest()
        else:
            upload = st.file_uploader(source, type=["sdf", "mol"] if source == "SDF file" else ["smi", "smiles", "txt"], key="molecular-import-smiles-file-" + source)
            if upload:
                if source == "SDF file":
                    records = imports.parse_sdf(upload.getvalue(), upload.name)
                else:
                    source_content = upload.getvalue()
                    analysis = imports.analyze_smiles_text(source_content.decode("utf-8-sig"))
                    records = analysis.pop("records")
                    audit = analysis
                    audit["source_name"] = Path(upload.name).name
                    audit["source_sha256"] = hashlib.sha256(source_content).hexdigest()
    elif method == "Draw molecule":
        from mn_molecular_design_studio.molecule_sketcher import render_sketcher
        name = st.text_input("Molecule name", value="Drawn compound", key="molecular-import-drawn-name")
        drawing = render_sketcher(key="molecular-import-drawing")
        if drawing and drawing.get("smiles"):
            records = imports.parse_smiles(drawing["smiles"])
            records[0]["name"] = name or "Drawn compound"
            if drawing.get("molblock"):
                records[0]["files"] = [{"name": "drawing.mol", "content": drawing["molblock"].encode(), "label": "Drawing", "dimensions": 2}]
            st.code(records[0]["components"][0]["smiles"], language=None)
    elif method == "Fetch by identifier":
        cols = st.columns([1, 3, 1])
        source = cols[0].selectbox("Database", ["PDB", "UniProt"], key="molecular-import-identifier-source")
        identifier = cols[1].text_input("Accession", placeholder="1IEP" if source == "PDB" else "P04637", key="molecular-import-identifier")
        if cols[2].button("Fetch", key="molecular-import-fetch", disabled=not identifier.strip()):
            filename, content = imports.fetch_identifier(identifier, source)
            st.session_state["molecular-import-fetched"] = {"key": (source, identifier.strip()), "filename": filename, "content": content}
        fetched = st.session_state.get("molecular-import-fetched")
        if fetched and fetched["key"] == (source, identifier.strip()):
            records = imports.parse_structure(fetched["content"], fetched["filename"]) if source == "PDB" else imports.parse_sequences(fetched["content"].decode())
    return records, audit, source_content


def _import_report_frame(rows):
    return [{
        "Row": row["row_number"],
        "Status": row["status"],
        "Source ID": row.get("source_id") or "",
        "Compound ID": row.get("compound_id") or "",
        "Name": row.get("name") or "",
        "Input SMILES": row.get("raw_smiles") or "",
        "Canonical isomeric SMILES": row.get("canonical_smiles") or "",
        "Annotation candidate (review only)": row.get("annotation_stripped_candidate") or "",
        "Candidate status": row.get("candidate_status") or "",
        "Duplicate of row": row.get("duplicate_of_row") or "",
        "Stereo variant of row": row.get("stereo_variant_of_row") or "",
        "Destination match": "Yes" if row.get("matched_existing") else "",
        "Warnings": "; ".join(row.get("warnings") or []),
        "Errors": "; ".join(row.get("errors") or []),
        "Formula": (row.get("descriptors") or {}).get("formula", ""),
        "Molecular weight": (row.get("descriptors") or {}).get("molecular_weight", ""),
        "Exact mass": (row.get("descriptors") or {}).get("exact_mass", ""),
        "Heavy atoms": (row.get("descriptors") or {}).get("heavy_atoms", ""),
        "Formal charge": (row.get("descriptors") or {}).get("formal_charge", ""),
        "Fragments": (row.get("descriptors") or {}).get("fragment_count", ""),
        "Disconnected components": json.dumps((row.get("descriptors") or {}).get("components", []), ensure_ascii=False),
        "HBD": (row.get("descriptors") or {}).get("hbd", ""),
        "HBA": (row.get("descriptors") or {}).get("hba", ""),
        "cLogP": (row.get("descriptors") or {}).get("clogp", ""),
        "TPSA": (row.get("descriptors") or {}).get("tpsa", ""),
        "QED": (row.get("descriptors") or {}).get("qed", ""),
        "Other source columns": json.dumps(row.get("source_data") or {}, ensure_ascii=False, default=str),
    } for row in rows]


def _report_csv(rows):
    import pandas as pd
    frame = pd.DataFrame(_import_report_frame(rows))
    return frame.to_csv(index=False).encode("utf-8-sig")


@st.dialog("Import molecules", width="large")
def render_import_dialog(user, projects):
    database = MolecularDatabase()
    step = st.session_state.get("molecular-import-step", 0)
    st.caption("Group  ›  What are you importing?" if step == 0 else "Group  ›  Import molecules")
    if step == 0:
        groups = [row for row in database.groups(user, collection="library") if row["owner_id"] == user["id"] or user.get("role") == "admin"]
        group_id = st.selectbox("Group *", ["", *[row["id"] for row in groups]], format_func=lambda value: next((row["name"] for row in groups if row["id"] == value), "+ Create a new group"), key="molecular-import-group")
        name = st.text_input("New group name", placeholder="e.g. ALK2 screening compounds", key="molecular-import-group-name") if not group_id else next(row["name"] for row in groups if row["id"] == group_id)
        project = _projects(projects, "molecular-import-project") if not group_id else next(row["project_id"] for row in groups if row["id"] == group_id)
        left, right = st.columns([1, 1])
        if left.button("Cancel", key="molecular-import-cancel-0"):
            _reset_import()
        if right.button("Next", type="primary", disabled=not name.strip(), key="molecular-import-next"):
            st.session_state["molecular-import-destination"] = {"group_id": group_id or None, "name": name, "project_id": project}
            st.session_state["molecular-import-instance"] = uuid.uuid4().hex
            st.session_state["molecular-import-step"] = 1
            st.rerun(scope="app")
        return
    destination = st.session_state["molecular-import-destination"]
    import_instance = st.session_state.setdefault("molecular-import-instance", uuid.uuid4().hex)
    st.markdown(f"**{html.escape(destination['name'])}**")
    method = st.radio("What are you importing?", ["Table", "FASTA", "Structure files", "Paste sequences", "SMILES", "Draw molecule", "Fetch by identifier"], horizontal=True, key="molecular-import-method")
    records = []
    audit = None
    source_content = None
    review_confirmed = True
    block_duplicates = False
    try:
        records, audit, source_content = _import_content(method)
        if audit:
            destination_group = destination.get("group_id")
            existing = database.group_identity_map(user, destination_group) if destination_group else {}
            existing_fingerprints = set(existing)
            existing_ids = {}
            for record in records:
                _, fingerprint, _ = imports.canonical_components(record["components"])
                if fingerprint in existing_fingerprints:
                    existing_ids[fingerprint] = existing[fingerprint]
            for row in audit["rows"]:
                if not row.get("canonical_smiles"):
                    continue
                _, fingerprint, _ = imports.canonical_components([
                    {"type": "small_molecule", "chain": "L", "smiles": row["canonical_smiles"]}
                ])
                row["matched_existing"] = fingerprint in existing_fingerprints
                row["existing_molecule_id"] = existing.get(fingerprint, {}).get("id", "")
            summary = audit["summary"]
            summary["existing_structure_matches"] = len(existing_ids)
            summary["reviewable_unique_structures"] = len(records)
            mapping = audit.get("mapping") or {}
            st.caption("Worksheet: " + (mapping.get("sheet_name") or "not applicable") +
                       " · SMILES: " + str(mapping.get("smiles_column") or "") +
                       " · ID: " + str(mapping.get("id_column") or "not mapped"))
            label = (f"{summary['rows']:,} source rows · {summary['valid_rows']:,} valid rows · "
                     f"{summary['unique_structures']:,} unique structures · {summary['invalid_rows']:,} invalid · "
                     f"{summary['duplicate_rows']:,} exact duplicates · {summary['empty_rows']:,} blank")
            st.caption(label)
            if summary["stereo_variants"]:
                st.info(f"{summary['stereo_variants']:,} row(s) have the same connectivity as another row but different stereochemistry. They remain separate structures.")
            if summary["duplicate_id_rows"]:
                st.caption(f"{summary['duplicate_id_rows']:,} duplicate source IDs were assigned unique compound IDs with a __2, __3… suffix; the original IDs remain in the report.")
            if summary["charged_component_rows"]:
                st.caption(f"{summary['charged_component_rows']:,} valid source row(s) contain at least one charged disconnected component.")
            if summary.get("existing_structure_matches"):
                st.caption(f"{summary['existing_structure_matches']:,} unique structure(s) are already in this destination. Existing names and properties will be preserved; new properties will be added where they are blank.")
            with st.expander("Review row validation and duplicates", expanded=bool(summary["invalid_rows"] or summary["duplicate_rows"])):
                view = st.selectbox("Rows", ["All rows", "Problems and duplicates", "Invalid or blank only"], key="molecular-import-report-filter-" + import_instance + audit["batch_id"])
                report_rows = audit["rows"]
                if view == "Problems and duplicates":
                    report_rows = [row for row in report_rows if row["status"] in {"invalid", "empty", "duplicate"} or row.get("warnings") or row.get("matched_existing")]
                elif view == "Invalid or blank only":
                    report_rows = [row for row in report_rows if row["status"] in {"invalid", "empty"}]
                report_frame = _import_report_frame(report_rows)
                st.dataframe(report_frame[:1000], hide_index=True, use_container_width=True)
                if len(report_frame) > 1000:
                    st.caption(f"Showing the first 1,000 of {len(report_frame):,} selected report rows. Download the full report below.")
                st.download_button(
                    "Download validation report CSV", _report_csv(audit["rows"]),
        file_name=Path(audit["source_name"]).stem + "-validation.csv", mime="text/csv",
                    key="molecular-import-report-download-" + import_instance + audit["batch_id"],
                )
                st.caption("Canonical SMILES are used for exact structure matching. Salts, formal charge, and stereochemistry are retained as submitted.")
            st.caption("The original file and full batch report are saved on this workstation. Imported molecules retain their source row and validation provenance in the Database.")
            review_confirmed = st.checkbox(
                "I reviewed the report. Only valid unique structures will be imported; every source row will be retained in the import history.",
                key="molecular-import-reviewed-" + import_instance + audit["batch_id"],
            )
            block_duplicates = st.checkbox(
                "Require a duplicate-free source file",
                value=False,
                disabled=not summary["duplicate_rows"],
                key="molecular-import-block-duplicates-" + import_instance + audit["batch_id"],
            )
        fingerprints = {imports.canonical_components(row["components"])[1] for row in records}
        if records and not audit:
            st.caption(f"{len(records):,} rows · {len(fingerprints):,} unique molecules. Existing molecules retain their saved scores.")
            st.dataframe([{"Name": row["name"], "Chains": ", ".join(component["chain"] for component in row["components"]), "Type": TYPE_LABELS[imports.canonical_components(row["components"])[2]]} for row in records[:10]], hide_index=True, use_container_width=True)
        elif audit and not records:
            st.warning("No valid structures are available to import. Correct the source or its column mapping and review the report again.")
    except (ValueError, OSError, UnicodeError) as exc:
        st.error(str(exc))
        records = []
        audit = None
        source_content = None
        review_confirmed = False
    back, cancel, commit = st.columns([1, 1, 2])
    if back.button("← Back", key="molecular-import-back"):
        st.session_state["molecular-import-step"] = 0
        st.rerun(scope="app")
    if cancel.button("Cancel", key="molecular-import-cancel-1"):
        _reset_import()
    can_commit = (bool(records) or bool(audit)) and review_confirmed and (not block_duplicates or not records)
    commit_label = "Save validation report" if audit and not records else "Import molecules"
    if commit.button(commit_label, type="primary", disabled=not can_commit, key="molecular-import-commit"):
        try:
            result = database.import_records(user, records, **destination, import_audit=audit, source_content=source_content)
            st.session_state["molecular-preselect-group"] = result["group_id"]
            st.session_state["molecular-preselect-collection"] = "library"
            if audit:
                summary = audit["summary"]
                if records:
                    st.session_state["molecular-database-flash"] = (
                        f"Imported {summary['valid_rows']:,} valid rows as {result['unique_count']:,} unique molecules. "
                        f"{summary['invalid_rows'] + summary['empty_rows']:,} invalid or blank rows were retained in the validation report."
                    )
                else:
                    st.session_state["molecular-database-flash"] = "Saved the validation report. No molecules were imported."
            else:
                st.session_state["molecular-database-flash"] = f"Imported {result['imported_count']:,} rows as {result['unique_count']:,} unique molecules."
            st.session_state["database-selection"] = []
            _reset_import()
        except (ValueError, PermissionError, OSError) as exc:
            st.error(str(exc))


from mn_molecular_design_studio.molecular_preparation import cropped_record


def _table_row(record):
    components = record["components"]
    metadata = record.get("metadata") or {}
    compound_id = (metadata.get("compound_id") or metadata.get("pubchem_cid") or record["name"]) if record["type"] == "small_molecule" else ""
    smiles = "; ".join(str(component.get("smiles") or "").strip() for component in components
                        if component.get("type", record["type"]) == "small_molecule" and component.get("smiles"))
    values = {"Compound ID": compound_id, "Name": record["name"], "SMILES": smiles,
              "Type": TYPE_LABELS[record["type"]], "Has structure": "Yes" if record["structures"] else "No",
              "Origin": record.get("origin", {}).get("tool") or "User upload", "Added": record["created_at"][:10],
              "Tags": metadata.get("tags", ""), "Status": metadata.get("status", ""),
              "Notes": metadata.get("notes", ""), **record.get("scores", {})}
    for component in components:
        label = component["chain"] + (" · lead" if component.get("role") == "lead" else "")
        values[label] = component.get("sequence") or component.get("smiles") or component.get("ccd")
    values.update({"Property · " + key: value for key, value in metadata.items()
                   if key not in {"compound_id", "compound_name"} and isinstance(value, (str, int, float, bool))})
    return {"id": record["id"], "values": values}


def _remote_source_workstation(origin):
    """Return the source workstation for records replicated from another node."""
    if not isinstance(origin, dict) or not origin.get("_replica_node"):
        return ""
    return str(origin.get("workstation") or "another workstation")


def _structure_text(path, context):
    if path.suffix.lower() in {".sdf", ".mol"}:
        from rdkit import Chem
        supplier = Chem.SDMolSupplier(str(path), removeHs=False)
        molecule = supplier[int(context.get("molecule_index", 0))]
        if molecule is None:
            raise ValueError("Could not read these ligand coordinates.")
        return Chem.MolToPDBBlock(molecule)
    if path.suffix.lower() == ".pdbqt":
        lines = []
        for line in path.read_text().splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                atom_type = line.split()[-1]
                element = {"A": "C", "NA": "N", "OA": "O", "SA": "S", "HD": "H"}.get(atom_type, atom_type.title())
                lines.append(line[:66].ljust(76) + element.rjust(2))
        return "\n".join([*lines, "END", ""])
    return path.read_text(errors="replace")


def render_details(database, user, record):
    from mn_molecular_design_studio.evaluations import evaluation_descriptor
    evaluation = record.get("selected_evaluation") or {}
    history = {evaluation_descriptor(obs)["id"]: evaluation_descriptor(obs) for obs in record.get("observations", [])}
    if history:
        current = (record.get("selected_evaluation") or {}).get("id")
        options = list(history)
        choice = st.selectbox("Evaluation", options, index=options.index(current) if current in options else 0,
            format_func=lambda value: f"{history[value]['evaluator']} · {history[value]['run_id']} · {history[value]['candidate_id'] or value[-8:]}",
            key="database-evaluation-" + record["id"])
        if choice != current:
            selected = database.list_molecules(user, molecule_ids=[record["id"]], evaluation_filters={"evaluation_id": choice})["items"]
            if selected:
                record = selected[0]
        evaluation = record.get("selected_evaluation") or {}
        st.caption("Scores and structure belong to this evaluation. Target: " + (evaluation.get("target_id") or "unspecified in historical result"))
    st.markdown("**" + html.escape(record["name"]) + "**")
    st.caption(record["id"])
    source_workstation = _remote_source_workstation(record.get("origin"))
    if source_workstation:
        st.caption("Synced from " + source_workstation)
    molecule_tab, scores_tab = st.tabs(["This molecule", "Group scores"])
    with molecule_tab:
        structures = record["structures"]
        if structures:
            # The available representations are filtered to the selected
            # evaluation. Include that evaluation in the widget key so
            # Streamlit does not carry an old engine's structure index across
            # a change in the Evaluation selector.
            evaluation_id = str(evaluation.get("id") or (choice if history else ""))
            structure_key = f"database-structure-{record['id']}-{evaluation_id or 'all'}"
            selected = st.selectbox("Structure", range(len(structures)), format_func=lambda index: structures[index]["label"], key=structure_key)
            structure = structures[selected]
            try:
                if source_workstation and str(structure.get("path", "")).startswith("cluster-library://"):
                    if structure.get("context", {}).get("portable"):
                        st.caption(f"Portable Database copy · syncing from {source_workstation}.")
                    else:
                        st.caption(f"Structure file is hosted by {source_workstation} and retrieved when opened.")
                path = database.structure_file(user, structure["id"])
                st.download_button("Download structure", path.read_bytes(), file_name=record["name"][:70] + path.suffix, key="database-download-" + structure["id"])
                if structure.get("context", {}).get("portable"):
                    st.caption("This structure is stored in the shared Database and will be copied to other workstations during sync.")
                else:
                    try:
                        database.group(structure["group_id"], user, write=True)
                        if st.button("Make available on all computers", key="database-portable-" + structure["id"]):
                            database.make_structure_portable(user, structure["id"])
                            from mn_molecular_design_studio.cluster.replication import Replica
                            from mn_molecular_design_studio.molecular_database_peer import capture
                            capture(Replica())
                            structure.setdefault("context", {})["portable"] = True
                            st.success("Saved in the shared Database. Other workstations will receive it during sync.")
                    except PermissionError:
                        pass
                if path.suffix.lower() in {".pdb", ".ent", ".cif", ".mmcif", ".sdf", ".mol", ".pdbqt"} and path.stat().st_size < 50_000_000:
                    from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import molstar_custom_component
                    from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer.dataclasses import StructureVisualization
                    structure_context = structure.get("context") or {}
                    policy = structure_context.get("database_structure_policy")
                    if record["type"] == "small_molecule" and policy in {"best_docking_pose", "best_prediction"}:
                        from mn_molecular_design_studio.docking_results import pose_visualizations
                        from mn_molecular_design_studio.bound_docking import receptor_file
                        receptor_text = None
                        receptor_name = ""
                        receptor_path = receptor_file(database, user, structure)
                        if receptor_path:
                            receptor_text = receptor_path.read_text(errors="replace")
                            receptor_name = receptor_path.name
                            st.caption("Compound and target are separate Database molecules, displayed together in their docking coordinates.")
                        representations = pose_visualizations(
                            path.read_text(errors="replace"), path.name,
                            receptor_text=receptor_text, receptor_name=receptor_name,
                            include_pocket_sidechains=False,
                        )
                    else:
                        # Ordinary protein structures stay in cartoon form; the
                        # specialized pose path above splits a ligand complex so
                        # multiple PDBQT poses cannot pile up in the viewer.
                        representations = [StructureVisualization(pdb=_structure_text(path, structure_context), color="chain-id", representation_type="ball-and-stick" if record["type"] == "small_molecule" else "cartoon")]
                    molstar_custom_component(representations, key="database-molstar-" + structure["id"], height=310, show_controls=False, show_structure_tools=False)
                from mn_molecular_design_studio.bound_docking import is_bound_pose
                if is_bound_pose(structure):
                    from mn_molecular_design_studio.workspace import load_organization
                    bundle_key = "database-bound-export-" + structure["id"]
                    if st.button("Prepare bound complex for export", key=bundle_key + "-prepare"):
                        from mn_molecular_design_studio.bound_docking import prepare_bound_inputs
                        st.session_state[bundle_key] = prepare_bound_inputs(database, user, record["id"], structure_id=structure["id"])
                    if bundle_key in st.session_state:
                        prepared = st.session_state[bundle_key]
                        st.download_button("Download complex PDB", Path(prepared["complex_pdb"]).read_bytes(),
                            file_name=record["name"][:70] + "-complex.pdb", key=bundle_key + "-pdb")
                        st.download_button("Download matching ligand SDF", Path(prepared["ligand_sdf"]).read_bytes(),
                            file_name=record["name"][:70] + "-ligand.sdf", key=bundle_key + "-sdf")
                    with st.popover("Use this bound pose in MD", use_container_width=True):
                        engine = st.selectbox("MD engine", ["OpenMM", "GROMACS"], key="database-bound-md-engine")
                        project = _projects(load_organization().get("projects", []), "database-bound-md-project",
                                            default=database.group(structure["group_id"], user)["project_id"], optional=True)
                        if st.button("Open MD inputs", type="primary", disabled=not project, key="database-bound-md-open"):
                            from mn_molecular_design_studio.database_selections import open_bound_pose_in_md
                            open_bound_pose_in_md(database, user, record["id"], "ligand.md." + engine.lower(), project,
                                                  structure_id=structure["id"])
                            st.rerun()
            except (OSError, ValueError, PermissionError) as exc:
                st.caption(str(exc))
        elif source_workstation:
            st.info(f"No structure file is attached to this synced Database record. Its source is {source_workstation}; the structure must be indexed on that workstation before it can be opened here.")
        for component in record["components"]:
            st.markdown(f"<p class='database-label'>{html.escape(component['chain'])} · {TYPE_LABELS[component['type']]}{(' · ' + html.escape(component['role'])) if component.get('role') else ''}</p>", unsafe_allow_html=True)
            if component.get("smiles"):
                from mn_molecular_design_studio.cofolding_input import _smiles_preview
                image = _smiles_preview(component["smiles"])
                if image:
                    st.image(image, use_container_width=True)
            st.code(component.get("sequence") or component.get("smiles") or "CCD: " + component["ccd"], language=None, wrap_lines=True)
        member_group = database.group(record["groups"][0]["id"], user)
        if not member_group["read_only"] and (user.get("role") == "admin" or member_group["owner_id"] == user["id"]):
            with st.expander("Chain roles"):
                selected_chain = st.selectbox("Chain", [row["chain"] for row in record["components"]], key="database-role-chain-" + record["id"])
                role = st.selectbox("Role", ["", "lead", "target", "heavy", "light", "cofactor"], key="database-chain-role-" + record["id"])
                if st.button("Save chain role", key="database-save-role-" + record["id"]):
                    database.set_chain_role(user, member_group["id"], record["id"], selected_chain, role)
                    st.rerun()
        st.caption("Groups: " + ", ".join(group["name"] for group in record["groups"]))
        with st.expander("Crop chains to a new molecule"):
            chains = [row["chain"] for row in record["components"]]
            chosen = st.multiselect("Keep chains", chains, default=chains, key="database-crop-" + record["id"])
            if st.button("Crop to new molecule", disabled=not chosen or set(chosen) == set(chains), key="database-crop-apply-" + record["id"]):
                try:
                    files = [{"path": database.structure_file(user, structures[0]["id"])}] if structures else []
                    cropped = cropped_record({**record, "files": files}, chosen)
                    group = database.group(record["groups"][0]["id"], user)
                    result = database.import_records(user, [cropped], name=record["name"] + " · cropped", project_id=group["project_id"], origin={"type": "user", "source_molecule_id": record["id"]})
                    st.session_state["molecular-preselect-group"] = result["group_id"]
                    st.rerun()
                except (ValueError, OSError, PermissionError) as exc:
                    st.error(str(exc))
        with st.expander("Edit name, notes, tags, and properties"):
            group_id = record["groups"][0]["id"]
            group = database.group(group_id, user)
            editable = not group["read_only"] and (user.get("role") == "admin" or group["owner_id"] == user["id"])
            with st.form("database-notes-" + record["id"]):
                molecule_name = st.text_input("Name", value=record["name"], disabled=not editable, key="database-edit-name-" + record["id"])
                tags = st.text_input("Tags", value=str(record.get("metadata", {}).get("tags") or ""), disabled=not editable, key="database-edit-tags-" + record["id"])
                notes = st.text_area("Notes", value=str(record.get("metadata", {}).get("notes") or ""), disabled=not editable, key="database-edit-notes-" + record["id"])
                properties = {}
                metadata = record.get("metadata") or {}
                schema_fields = {field["name"]: field for field in group["schema"].get("fields", [])}
                computed_fields = {name for name, field in schema_fields.items() if field.get("formula")}
                for property_name, value in metadata.items():
                    if property_name in {"tags", "notes"} or property_name in computed_fields:
                        continue
                    property_key = hashlib.sha1(property_name.encode("utf-8")).hexdigest()[:12]
                    field = schema_fields.get(property_name) or {}
                    if field.get("type") == "number":
                        try:
                            number_value = float(value) if value not in (None, "") else None
                        except (TypeError, ValueError):
                            number_value = None
                        properties[property_name] = st.number_input(
                            property_name, value=number_value, disabled=not editable,
                            key=f"database-edit-schema-{record['id']}-{property_key}",
                        )
                    elif field.get("type") == "boolean":
                        normalized = value
                        if isinstance(value, str):
                            normalized = True if value.casefold() in {"true", "1", "yes"} else False if value.casefold() in {"false", "0", "no"} else None
                        options = [None, True, False]
                        properties[property_name] = st.selectbox(
                            property_name, options, index=options.index(normalized) if normalized in options else 0,
                            disabled=not editable, format_func=lambda item: "(unset)" if item is None else str(item),
                            key=f"database-edit-schema-{record['id']}-{property_key}",
                        )
                    elif field.get("type") == "category" and field.get("choices"):
                        options = ["", *field["choices"]]
                        properties[property_name] = st.selectbox(
                            property_name, options, index=options.index(value) if value in options else 0,
                            disabled=not editable, key=f"database-edit-schema-{record['id']}-{property_key}",
                        )
                    elif isinstance(value, bool):
                        options = [None, True, False]
                        properties[property_name] = st.selectbox(
                            property_name, options, index=options.index(value), disabled=not editable,
                            format_func=lambda item: "(unset)" if item is None else str(item),
                            key=f"database-edit-property-{record['id']}-{property_key}",
                        )
                    elif isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value)):
                        properties[property_name] = st.number_input(
                            property_name, value=float(value), disabled=not editable,
                            key=f"database-edit-property-{record['id']}-{property_key}",
                        )
                    elif isinstance(value, str):
                        if len(value) > 120:
                            properties[property_name] = st.text_area(property_name, value=value, height=90, disabled=not editable,
                                key=f"database-edit-property-{record['id']}-{property_key}")
                        else:
                            properties[property_name] = st.text_input(property_name, value=value, disabled=not editable,
                                key=f"database-edit-property-{record['id']}-{property_key}")
                    else:
                        st.caption(f"{property_name}: structured source value is preserved as read-only JSON.")
                for field in group["schema"].get("fields", []):
                    if field.get("formula") or field["name"] in {"tags", "notes"} or field["name"] in properties:
                        continue
                    label = field["name"]
                    value = metadata.get(label)
                    field_key = hashlib.sha1(label.encode("utf-8")).hexdigest()[:12]
                    if field.get("type") == "number":
                        properties[label] = st.number_input(label, value=float(value) if value is not None else None, disabled=not editable,
                            key=f"database-edit-schema-{record['id']}-{field_key}")
                    elif field.get("type") == "boolean":
                        properties[label] = st.selectbox(label, [None, True, False], index=[None, True, False].index(value) if value in [None, True, False] else 0, disabled=not editable,
                            key=f"database-edit-schema-{record['id']}-{field_key}")
                    elif field.get("type") == "category" and field.get("choices"):
                        options = ["", *field["choices"]]
                        properties[label] = st.selectbox(label, options, index=options.index(value) if value in options else 0, disabled=not editable,
                            key=f"database-edit-schema-{record['id']}-{field_key}")
                    else:
                        properties[label] = st.text_input(label, value=str(value or ""), disabled=not editable,
                            key=f"database-edit-schema-{record['id']}-{field_key}")
                st.markdown("**Add a property**")
                new_property_name = st.text_input("Property name", disabled=not editable, key="database-add-property-name-" + record["id"])
                new_property_value = st.text_input("Property value", disabled=not editable, key="database-add-property-value-" + record["id"])
                if st.form_submit_button("Save", disabled=not editable):
                    try:
                        updates = {"tags": tags, "notes": notes, **properties}
                        if new_property_name.strip():
                            if new_property_name.strip() in metadata:
                                raise ValueError("That property already exists. Edit it above instead.")
                            updates[new_property_name.strip()] = new_property_value
                        database.update_metadata(user, group_id, record["id"], updates, name=molecule_name)
                        st.rerun()
                    except ValueError as exc:
                        st.error(str(exc))
        with st.expander("History and sources"):
            for source in record["sources"]:
                st.write(source.get("tool") or source.get("type") or "Import")
                if source.get("type") == "table_import_row":
                    st.caption(f"{source.get('source_name') or 'Source file'} · row {source.get('source_row', '?')}" +
                               (f" · ID {source['source_id']}" if source.get("source_id") else "") +
                               (f" · {source.get('validation_status')}" if source.get("validation_status") else ""))
                    if source.get("warnings"):
                        st.caption("Warnings: " + "; ".join(source["warnings"]))
                if source.get("run_id"):
                    st.caption(source["run_id"])
                    if st.button("Open job", key="database-open-job-" + hashlib.sha1(source["source_key"].encode()).hexdigest()):
                        st.session_state["selected-job-reference"] = source["task_group"] + "/" + source["run_id"]
                        st.session_state["page"] = "Job details"
                        st.rerun()
                if source.get("parents"):
                    st.json(source["parents"])
    with scores_tab:
        if not record["observations"]:
            st.caption("No scores recorded yet.")
        else:
            st.caption("The sheet uses one selected evaluation. All recorded runs remain available here.")
            common_rows = []
            for observation in record["observations"]:
                output = observation.get("normalized_output") or {}
                for key, value in output.get("metrics", {}).items():
                    info = output["provenance"][key]
                    common_rows.append({"Measurement": field_label(PREFIX + key), "Value": value,
                        "Tool": info["source_tool"], "Original field": ", ".join(item["source_field"] for item in info["aliases"]),
                        "Method": info["method"], "Scope": info["scope"], "Run": observation["run_id"]})
            if common_rows:
                st.dataframe(common_rows, hide_index=True, use_container_width=True)
                st.caption("This history includes separate evaluations; each source run and original value is retained below.")
            rows = []
            for observation in record["observations"]:
                for key, value in observation["metrics"].items():
                    rows.append({"Tool": observation["namespace"], "Metric": key, "Value": value, "Run": observation["run_id"], "Recorded": observation["created_at"][:19], "Units": observation["context"].get("units") or ""})
            with st.expander("Original tool scores", expanded=not bool(common_rows)):
                st.dataframe(rows, hide_index=True, use_container_width=True)
            with st.expander("Calculation context"):
                observation = st.selectbox("Run", range(len(record["observations"])), format_func=lambda index: record["observations"][index]["namespace"] + " · " + record["observations"][index]["run_id"], key="database-context-" + record["id"])
                st.json(record["observations"][observation]["context"])


@st.dialog("Create molecule group", width="small")
def render_group_dialog(user, projects, ids, default_project=""):
    database = MolecularDatabase()
    st.caption(f"{len(ids):,} selected molecules. Their existing structures and scores stay linked.")
    name = st.text_input("Group name", key="molecular-group-name")
    project = _projects(projects, "molecular-group-project", default=default_project)
    if st.button("Create group", type="primary", disabled=not name.strip(), key="molecular-group-save"):
        group_id = database.create_group(user, name, project_id=project)
        database.add_to_group(user, group_id, ids)
        st.session_state["molecular-preselect-group"] = group_id
        st.session_state["molecular-preselect-collection"] = "library"
        st.session_state["database-group-request"] = None
        st.rerun()


@st.dialog("Add column", width="small")
def render_column_dialog(user, group_id):
    database = MolecularDatabase()
    group = database.group(group_id, user, write=True)
    name = st.text_input("Column name", key="molecular-column-name")
    kind = st.selectbox("Type", ["text", "number", "boolean", "category"], key="molecular-column-type")
    formula = st.text_input("Formula (optional)", placeholder="=9 - LOG10(affinity_nM)", key="molecular-column-formula")
    required = st.checkbox("Required", disabled=bool(formula), key="molecular-column-required")
    categories = st.text_input("Allowed categories (comma separated)", key="molecular-column-categories") if kind == "category" else ""
    if st.button("Add column", type="primary", disabled=not name.strip(), key="molecular-column-save"):
        try:
            schema = group["schema"]
            field = {"name": name.strip(), "type": kind, "required": required and not formula}
            if formula:
                field["formula"] = formula
            if categories:
                field["choices"] = [value.strip() for value in categories.split(",") if value.strip()]
            database.set_schema(user, group_id, {**schema, "fields": [*schema.get("fields", []), field]})
            st.session_state["molecular-column-open"] = False
            st.rerun()
        except (ValueError, PermissionError) as exc:
            st.error(str(exc))


def render_database(user, projects):
    database = MolecularDatabase()
    with st.spinner("Updating molecular records…"):
        errors = sync_workspace(user, database=database)
    st.markdown(CSS, unsafe_allow_html=True)
    pending_group = st.session_state.pop("molecular-preselect-group", None)
    pending_collection = st.session_state.pop("molecular-preselect-collection", None)
    pending_molecule = st.session_state.pop("molecular-preselect-molecule", None)
    if pending_collection in {"library", "all"}:
        st.session_state["molecular-collection-view"] = "Library" if pending_collection == "library" else "All indexed"
    if pending_group:
        destination = database.group(pending_group, user)
        st.session_state["molecular-project-filter"] = destination["project_id"]
        st.session_state["molecular-scope"] = (
            "Shared" if destination["project_id"] else "Mine" if destination["owner_id"] == user["id"] else "All"
        )
        st.session_state["molecular-group-search"] = ""
        st.session_state["molecular-query"] = ""
        st.session_state["molecular-type-filter"] = ""
        st.session_state["molecular-score-filters"] = []
        for method in IMPORT_METHODS:
            st.session_state["molecular-method-" + method] = True
        st.session_state["database-selected-group"] = pending_group
    if pending_molecule:
        st.session_state["database-focus"] = pending_molecule
    if message := st.session_state.pop("molecular-database-flash", None):
        st.success(message)
    if errors:
        with st.expander(f"{len(errors)} source jobs could not be indexed"):
            st.dataframe(errors, hide_index=True)
    with st.container(key="molecular-database"):
        sidebar, main, details = st.columns([.23, .51, .26], gap="medium")
        with sidebar, st.container(key="database-groups"):
            heading, new, refresh = st.columns([3, 1.5, 1])
            heading.markdown("### Groups")
            if new.button("+ New", key="molecular-database-new"):
                st.session_state["molecular-import-open"] = True
                st.session_state["molecular-import-step"] = 0
            if refresh.button("↻", key="molecular-database-refresh", help="Refresh molecules and job results"):
                sync_workspace(user, database=database, force=True)
                st.rerun()
            group_query = st.text_input("Search groups", placeholder="Search groups…", label_visibility="collapsed", key="molecular-group-search")
            collection_label = st.segmented_control(
                "Collection",
                ["Library", "All indexed"],
                default="Library",
                key="molecular-collection-view",
                help="Library shows imported, curated, and final generated collections. All indexed also shows automatic and intermediate job collections.",
            )
            collection = "library" if collection_label == "Library" else "all"
            scope_options = ["Mine", "Shared", "All"]
            old_scope = st.session_state.get("molecular-scope")
            scope_aliases = {"My groups": "Mine", "Organization": "All", "organization": "All", "mine": "Mine", "shared": "Shared", "all": "All"}
            if old_scope in scope_aliases:
                st.session_state["molecular-scope"] = scope_aliases[old_scope]
            if st.session_state.get("molecular-scope") not in scope_options:
                st.session_state["molecular-scope"] = "All"
            scope_label = st.segmented_control(
                "Collection access",
                scope_options,
                default="All",
                key="molecular-scope",
                label_visibility="visible",
                help="Mine shows personal collections. Shared shows project collections. All combines both.",
            )
            scope = {"Mine": "mine", "Shared": "shared", "All": "all"}[scope_label]
            st.caption("Mine: personal · Shared: project collections · All: both")
            hide_smoke_test = st.checkbox(
                "Hide smoke/test groups",
                value=True,
                key="molecular-hide-smoke-test-groups",
                help="Hides collections whose name or source is marked smoke/test. Turn this off to inspect them.",
            )
            with st.popover("Import Method", use_container_width=True):
                methods = [value for value, label in IMPORT_METHODS.items() if st.checkbox(label, value=True, key="molecular-method-" + value)]
            project_id = _projects(projects, "molecular-project-filter", label="Project filter")
            if not methods:
                methods = ["__none__"]
            groups = database.groups(
                user, project_id=project_id, scope=scope, methods=methods, collection=collection,
                hide_smoke_test=hide_smoke_test,
            )
            filtered_groups = [row for row in groups if not group_query or group_query.casefold() in row["name"].casefold()]
            options = ["", *[row["id"] for row in filtered_groups]]
            if st.session_state.get("database-selected-group") not in options:
                st.session_state["database-selected-group"] = ""
            with st.container(height=560, border=False):
                def group_label(value):
                    row = next((row for row in filtered_groups if row["id"] == value), None)
                    if row is None:
                        return "All molecules"
                    label = f"{row['name']} · {row['molecule_count']:,}"
                    source = _remote_source_workstation(row.get("origin"))
                    return label + (f" · from {source}" if source else "")
                group_id = st.radio("Molecule groups", options, format_func=group_label, key="database-selected-group", label_visibility="collapsed")
        group_name = next((row["name"] for row in groups if row["id"] == group_id), "All molecules")
        active_group = next((row for row in groups if row["id"] == group_id), {})
        if st.session_state.get("database-comparison-group") != group_id:
            st.session_state["database-comparison-group"] = group_id
            if (active_group.get("origin") or {}).get("comparison_spec"):
                st.session_state["molecular-plot"] = True
        query_args = {"group_id": group_id, "project_id": project_id, "scope": scope, "methods": methods,
                      "collection": collection, "hide_smoke_test": hide_smoke_test}
        metric_fields = database.metric_fields(user, **query_args)
        with main:
            title, file_col, edit_col, view_col, export_col = st.columns([3, 1, 1, 1, 1.3])
            title.markdown("### Molecules")
            with file_col.popover("File", use_container_width=True):
                if st.button("Upload molecules…", key="molecular-database-upload"):
                    st.session_state["molecular-import-open"] = True
                    st.session_state["molecular-import-step"] = 0
                group_from_filter = st.button("Group from filter", key="molecular-group-filter")
            with edit_col.popover("Edit", use_container_width=True):
                editable_group = next((row for row in groups if row["id"] == group_id), None)
                can_edit = editable_group and not editable_group["read_only"] and (user.get("role") == "admin" or editable_group["owner_id"] == user["id"])
                if st.button("Add column…", disabled=not can_edit, key="molecular-add-column"):
                    st.session_state["molecular-column-open"] = True
            with view_col.popover("View", use_container_width=True):
                plot = st.checkbox("Compare candidates", key="molecular-plot")
            if group_id:
                batches = database.import_batches(user, group_id)
                if batches:
                    show_import_history = st.toggle(
                        f"Show import validation history · {len(batches):,} batch(es)",
                        key="molecular-show-import-history-" + group_id,
                    )
                    if show_import_history:
                        selected_batch_id = st.selectbox(
                            "Import batch", [row["id"] for row in batches],
                            format_func=lambda value: next(
                                f"{row['source_name']} · {row['created_at'][:19]} · {row['summary'].get('valid_rows', 0):,} valid / {row['summary'].get('invalid_rows', 0):,} invalid"
                                for row in batches if row["id"] == value
                            ),
                            key="molecular-import-history-batch-" + group_id,
                        )
                        selected_batch = next(row for row in batches if row["id"] == selected_batch_id)
                        summary = selected_batch["summary"]
                        mapping = selected_batch.get("mapping") or {}
                        st.caption("Worksheet: " + (mapping.get("sheet_name") or "not applicable") +
                                   " · SMILES: " + str(mapping.get("smiles_column") or "") +
                                   " · ID: " + str(mapping.get("id_column") or "not mapped"))
                        st.json(summary)
                        saved_rows = database.import_batch_rows(user, selected_batch_id)
                        history_frame = _import_report_frame(saved_rows)
                        st.dataframe(history_frame[:1000], hide_index=True, use_container_width=True)
                        if len(history_frame) > 1000:
                            st.caption(f"Showing the first 1,000 of {len(history_frame):,} rows. Download the full report for all rows.")
                        st.download_button(
                            "Download saved validation report", _report_csv(saved_rows),
                            file_name=Path(selected_batch["source_name"]).stem + "-validation.csv",
                            mime="text/csv", key="molecular-import-history-report-" + selected_batch_id,
                        )
                        try:
                            original = database.import_source_file(user, selected_batch_id)
                            st.download_button(
                                "Download original source file", original.read_bytes(),
                                file_name=selected_batch["source_name"],
                                key="molecular-import-history-source-" + selected_batch_id,
                            )
                        except (OSError, ValueError, PermissionError) as exc:
                            st.caption(str(exc))
                elif active_group.get("read_only") and _remote_source_workstation(active_group.get("origin")):
                    st.caption("Full import reports are stored on the source workstation. Per-molecule import row provenance is available under History and sources.")
            evaluation_options = database.evaluation_options(
                user, group_id=group_id, project_id=project_id, collection=collection,
                scope=scope, hide_smoke_test=hide_smoke_test,
            )
            evaluation_filters = {}
            with st.expander("Choose the scores shown in the molecule table", expanded=False):
                for field, label in (("target_id", "Target"), ("evaluator", "Tool"), ("run_id", "Run")):
                    options = ["", *evaluation_options[field]]
                    key = "database-context-" + field + "-" + (group_id or "all")
                    if st.session_state.get(key, "") not in options:
                        st.session_state[key] = ""
                    value = st.selectbox(label, options, format_func=lambda item: item or "Latest matching result", key=key)
                    if value:
                        evaluation_filters[field] = value
                st.caption("Each molecule row shows scores from one result: the latest representative model or pose by default. Use these selectors to choose a target, tool or run. Other model results remain available on Results.")
            search, sort_col, order_col, filter_col = st.columns([4, 1.5, 1, 1])
            query = search.text_input("Search molecules", placeholder="Search name, sequence, SMILES, tag…", label_visibility="collapsed", key="molecular-query")
            sort_options = ["added", "name", *metric_fields]
            if st.session_state.get("molecular-sort") not in sort_options:
                st.session_state["molecular-sort"] = "added"
            sort = sort_col.selectbox("Sort", sort_options, format_func=lambda value: value.title() if value in {"added", "name"} else field_label(value), key="molecular-sort", label_visibility="collapsed")
            descending = order_col.toggle("Desc", value=True, key="molecular-desc")
            filters = []
            with filter_col.popover("Filter", use_container_width=True):
                kind = st.selectbox("Molecule type", ["", *TYPE_LABELS], format_func=lambda value: TYPE_LABELS.get(value, "All types"), key="molecular-type-filter")
                if not set(st.session_state.get("molecular-score-filters") or []).issubset(metric_fields):
                    st.session_state["molecular-score-filters"] = []
                score_fields = st.multiselect("Score filters", metric_fields, format_func=field_label, key="molecular-score-filters")
                for field in score_fields:
                    left, right = st.columns([1, 2])
                    operator = left.selectbox(field + " comparison", [">=", "<=", ">", "<", "="], key="molecular-op-" + field, label_visibility="collapsed")
                    value = right.number_input(field_label(field), value=0.0, key="molecular-value-" + field)
                    filters.append({"field": field, "operator": operator, "value": value})
            args = {**query_args, "query": query, "molecule_type": kind, "sort": sort, "descending": descending, "metric_filters": filters,
                    "evaluation_filters": evaluation_filters or None}
            signature = hashlib.sha1(json.dumps(args, sort_keys=True).encode()).hexdigest()
            if st.session_state.get("molecular-view-signature") != signature:
                st.session_state["database-selection"] = []
                st.session_state["molecular-page"] = 1
                st.session_state["molecular-view-signature"] = signature
            page_number = int(st.session_state.get("molecular-page", 1))
            page = database.list_molecules(user, **args, limit=100, offset=(page_number-1)*100)
            records = page["items"]
            rows = [_table_row(record) for record in records]
            basic_columns = ["Compound ID", "Name", "SMILES", "Type", "Has structure"]
            chain_columns = list(dict.fromkeys(key for row in rows for key in row["values"] if key not in {"Compound ID", "Name", "SMILES", "Type", "Has structure", "Origin", "Added", "Tags", "Notes", "Status", *metric_fields} and not key.startswith("Property · ")))
            all_columns = [*basic_columns, *chain_columns, *metric_fields, "Origin", "Added", "Tags", "Notes", "Status", *list(dict.fromkeys(key for row in rows for key in row["values"] if key.startswith("Property · ")))]
            with st.popover("Columns", icon=":material/view_column:"):
                key = "molecular-columns-" + hashlib.sha1((group_id + "|".join(all_columns)).encode()).hexdigest()[:12]
                preferred = [PREFIX + key for key in ("design.sequence_length", "geometry.binder.ca_rg", "secondary_structure.binder.dssp.helix_fraction") if PREFIX + key in metric_fields]
                columns = st.multiselect("Visible columns", all_columns, default=[*basic_columns, *chain_columns[:2], *(preferred or metric_fields[:3])], format_func=field_label, key=key)
                st.caption("Common measurements share a definition and unit across tools. Original tool fields are also available.")
            columns = list(dict.fromkeys(["Compound ID", "Name", *columns]))
            if records:
                from mn_molecular_design_studio.molecule_table import render_table
                current_focus = st.session_state.get("database-focus")
                if current_focus not in {row["id"] for row in records}:
                    current_focus = records[0]["id"]
                    st.session_state["database-focus"] = current_focus
                event = render_table(rows, [{"key": field, "label": field_label(field), "kind": "sequence" if field in chain_columns and any(row["values"].get(field) and str(row["values"][field]).isalpha() for row in rows) else "value"} for field in columns], selected_ids=st.session_state.get("database-selection", []), focus_id=current_focus, key="database-sheet-" + signature + str(page_number))
                if event and event.get("token") != st.session_state.get("database-sheet-event"):
                    st.session_state["database-sheet-event"] = event["token"]
                    st.session_state["database-selection"] = event.get("selectedIds") or []
                    st.session_state["database-focus"] = event.get("focusId") or current_focus
                    st.rerun()
            else:
                empty_message = (
                    "No molecules to show in this Library view.<br>Import molecules or save selected results to build a reusable collection."
                    if collection == "library" and not groups and not query and not kind
                    else "No molecules to show.<br>Adjust your group filters or select another collection."
                )
                st.markdown(f'<div class="database-empty">{empty_message}</div>', unsafe_allow_html=True)
            previous, status, following = st.columns([1, 5, 1])
            if previous.button("‹", disabled=page_number <= 1, key="molecular-page-prev"):
                st.session_state["molecular-page"] = page_number - 1
                st.rerun()
            status.caption(f"{page['total']:,} molecules · page {page_number} of {max(1, math.ceil(page['total']/100))}")
            if following.button("›", disabled=page["next_offset"] is None, key="molecular-page-next"):
                st.session_state["molecular-page"] = page_number + 1
                st.rerun()
            selected_ids = st.session_state.get("database-selection") or []
            if selected_ids:
                st.caption(f"{len(selected_ids):,} molecules selected")
                action, clear = st.columns([3, 1])
                if action.button("Create group from selection", key="molecular-group-selected"):
                    st.session_state["database-group-request"] = selected_ids
                if clear.button("Clear", key="molecular-selection-clear"):
                    st.session_state["database-selection"] = []
                    st.rerun()
                selected_group = next((row for row in groups if row["id"] == group_id), None)
                can_remove = bool(
                    selected_group and not selected_group["read_only"]
                    and (user.get("role") == "admin" or selected_group["owner_id"] == user["id"])
                )
                if not group_id:
                    st.caption("Select a group to remove molecules from it. Removing a row only removes it from that group.")
                elif can_remove:
                    if st.button("Remove selected from Database", key="molecular-remove-selected"):
                        st.session_state["database-pending-removal"] = {"group_id": group_id, "molecule_ids": list(selected_ids)}
                else:
                    st.caption("You need to own this group to remove molecules from it.")
                pending_removal = st.session_state.get("database-pending-removal")
                if pending_removal and pending_removal.get("group_id") == group_id:
                    st.warning(f"Remove {len(pending_removal['molecule_ids']):,} selected molecule row(s) from “{group_name}”? The molecule remains in any other groups.")
                    confirm, cancel = st.columns(2)
                    if confirm.button("Confirm removal", type="primary", key="molecular-remove-confirm"):
                        try:
                            removed = database.remove_from_group(user, group_id, pending_removal["molecule_ids"])
                            st.session_state["database-selection"] = []
                            st.session_state.pop("database-pending-removal", None)
                            st.session_state["molecular-database-flash"] = f"Removed {removed:,} molecule row(s) from {group_name}."
                            st.rerun()
                        except (ValueError, PermissionError, sqlite3.Error) as exc:
                            st.error(str(exc))
                    if cancel.button("Cancel removal", key="molecular-remove-cancel"):
                        st.session_state.pop("database-pending-removal", None)
                        st.rerun()
                render_use_selection(database, user, projects, selected_ids, project_id)
            if group_from_filter:
                st.session_state["database-group-request"] = [row["id"] for row in database.list_molecules(user, **args, limit=100_000)["items"]]
            with export_col.popover("Export", use_container_width=True):
                export_scope = st.radio("Molecules", (["Selected molecules"] if selected_ids else []) + ["Filtered molecules", "Current page"], key="molecular-export-scope")
                export_format = st.selectbox("Format", ["CSV", "FASTA", "SMILES", "SDF", "Structures ZIP"], key="molecular-export-format")
                export_signature = signature + json.dumps([selected_ids, page_number, export_scope, export_format])
                if st.button("Prepare download", disabled=not page["total"], key="molecular-export-prepare"):
                    try:
                        from mn_molecular_design_studio.molecule_exports import export_molecules
                        export_ids = selected_ids if export_scope == "Selected molecules" else [row["id"] for row in records] if export_scope == "Current page" else [row["id"] for row in database.list_molecules(user, **args, limit=100_000)["items"]]
                        content, suffix, mime = export_molecules(database, user, export_ids, export_format)
                        st.session_state["molecular-export-ready"] = {"signature": export_signature, "data": content, "suffix": suffix, "mime": mime}
                    except (ValueError, OSError, PermissionError) as exc:
                        st.error(str(exc))
                ready = st.session_state.get("molecular-export-ready") or {}
                if ready.get("signature") == export_signature:
                    st.download_button("Download " + export_format, ready["data"], file_name=group_name[:60].replace(" ", "-") + "." + ready["suffix"], mime=ready["mime"], key="molecular-export-download")
            if plot and page["total"]:
                from mn_molecular_design_studio.comparison_ui import render_comparison
                comparison_records = database.list_molecules(user, **args, limit=100000)["items"]
                st.caption("Comparison includes the entire filtered collection, independently of table pagination.")
                saved_settings = (database.group(group_id, user).get("origin", {}).get("comparison_spec") or {}).get("settings") if group_id else None
                show_all_default = bool(group_id and len(evaluation_options.get("run_id") or []) > 1)
                all_evaluations = st.checkbox(
                    "Show every tool result as a separate row",
                    value=show_all_default,
                    key="database-show-all-evaluations-" + (group_id or "all"),
                    help="Each stored prediction, docking, affinity or analysis run gets its own row. Turn this off to show one selected result per compound.",
                )
                comparison_rows = database.evaluation_rows(user, **query_args, query=query, molecule_type=kind,
                    molecule_ids=[row["id"] for row in comparison_records], evaluation_filters=evaluation_filters or None) if all_evaluations else comparison_records
                render_comparison(comparison_rows, key="database-comparison-" + signature, database=database, user=user,
                                  project_id=project_id, selected_ids_key="database-selection",
                                  settings=saved_settings)
        with details, st.container(key="database-details"):
            record = next((row for row in records if row["id"] == st.session_state.get("database-focus")), records[0] if records else None)
            focus = st.session_state.get("database-focus")
            if focus and not any(row["id"] == focus for row in records):
                focused = database.list_molecules(
                    user, project_id=project_id, scope=scope, methods=methods, collection=collection,
                    molecule_ids=[focus], evaluation_filters=evaluation_filters or None,
                    hide_smoke_test=hide_smoke_test,
                )["items"]
                if focused:
                    record = focused[0]
            if record:
                render_details(database, user, record)
            else:
                st.caption("Select a molecule to inspect its structure, chains, and scores.")
    if st.session_state.get("molecular-import-open"):
        render_import_dialog(user, projects)
    elif st.session_state.get("database-group-request"):
        render_group_dialog(user, projects, st.session_state["database-group-request"], project_id)
    elif st.session_state.get("molecular-column-open") and group_id:
        render_column_dialog(user, group_id)


def render_use_selection(database, user, projects, ids, project_id=""):
    from mn_molecular_design_studio.database_selections import open_in_tool
    records = database.molecules_by_ids(user, ids)
    compounds = all(len(row["components"]) == 1 and row["components"][0].get("smiles") for row in records)
    options = {}
    if compounds:
        options.update({name: workflow for name, workflow in [("AutoDock Vina", "ligand.docking.vina"), ("GNINA", "ligand.docking.gnina"), ("Uni-Dock Pro", "ligand.docking.unidock_pro"), ("RosettaLigand", "ligand.docking.rosetta_ligand"), ("Nesso-1", "ligand.affinity.nesso")]})
    if all(any(component["type"] == "protein" for component in row["components"]) for row in records):
        options["Pipeline canvas"] = "pipeline"
        from mn_molecular_design_studio.quality_tools import WORKFLOWS
        options.update({profile["name"]: workflow for workflow, profile in WORKFLOWS.items()})
    if records:
        options["Boltz-2"] = "boltz2"
    if len(records) == 1:
        detail = database.list_molecules(user, molecule_ids=ids, collection="all")["items"][0]
        from mn_molecular_design_studio.bound_docking import is_bound_pose
        if any(is_bound_pose(structure) for structure in detail["structures"]):
            options.update({"OpenMM · bound pose": "ligand.md.openmm", "GROMACS · bound pose": "ligand.md.gromacs"})
        if any(row["type"] == "protein" for row in records[0]["components"]) and detail["structures"]:
            options.update({name + " · receptor": workflow + ":receptor" for name, workflow in [("AutoDock Vina", "ligand.docking.vina"), ("GNINA", "ligand.docking.gnina"), ("Uni-Dock Pro", "ligand.docking.unidock_pro"), ("RosettaLigand", "ligand.docking.rosetta_ligand")]})
    if not options:
        st.caption("Create a group or export the selected molecules for reuse.")
        return
    with st.popover("Use in tool", use_container_width=True):
        tool = st.selectbox("Tool", list(options), key="database-use-tool")
        selected_project = _projects(projects, "database-use-project", default=project_id, optional=True)
        if st.button("Open inputs", type="primary", disabled=not selected_project, key="database-open-tool"):
            open_in_tool(database, user, ids, options[tool], selected_project)
            st.rerun()


def render_molecule_picker(user, *, key, allowed_types=None, multiple=False):
    database = MolecularDatabase()
    with st.popover("My molecules", icon=":material/storage:"):
        query = st.text_input("Search Database", key=key + "-library-query")
        groups = database.groups(user, collection="library")
        group_id = st.selectbox("Group", ["", *[row["id"] for row in groups]], format_func=lambda value: next((row["name"] for row in groups if row["id"] == value), "All molecules"), key=key + "-library-group")
        rows = database.list_molecules(user, query=query, group_id=group_id, limit=100, collection="library")["items"]
        if allowed_types:
            rows = [row for row in rows if all(component["type"] in allowed_types for component in row["components"])]
        if not rows:
            st.caption("No matching molecules. Import them in Database first.")
            return None
        picker = st.multiselect if multiple else st.selectbox
        selected = picker("Molecules" if multiple else "Molecule", [row["id"] for row in rows], format_func=lambda value: next(row["name"] for row in rows if row["id"] == value), key=key + "-library-molecule")
        if st.button("Use selected molecules" if multiple else "Use molecule", key=key + "-library-use", disabled=multiple and not selected):
            return [row for row in rows if row["id"] in selected] if multiple else next(row for row in rows if row["id"] == selected)
    return None
