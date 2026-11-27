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

from mn_molecular_design_studio.result_details_ui import render_result_details, result_field_label, result_value_text
from mn_molecular_design_studio import molecule_imports as imports
from mn_molecular_design_studio.database_browser_labels import analysis_status, tool_label
from mn_molecular_design_studio.molecular_database import MolecularDatabase, IMPORT_METHODS
from mn_molecular_design_studio.molecular_database_ingest import sync_workspace_async as sync_workspace, workspace_sync_running
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
.st-key-database-groups{background:color-mix(in srgb,currentColor 4%,transparent);border-radius:12px;padding:14px;border-right:0}
.st-key-database-groups [data-testid="stRadio"] label{padding:9px 8px;border-radius:7px;overflow-wrap:anywhere}
.st-key-database-groups [data-testid="stRadio"] label:has(input:checked){box-shadow:inset 3px 0 #36815b}
.database-collection-heading{font-size:1.35rem;font-weight:650;line-height:1.3;margin:.35rem 0;overflow-wrap:anywhere}
.database-collection-meta{display:flex;flex-wrap:wrap;gap:6px;margin:8px 0 14px}
.database-collection-meta span{background:color-mix(in srgb,currentColor 4%,transparent);border:1px solid #d6e2d9;border-radius:20px;padding:4px 10px;font-size:12px}
@media(max-width:900px){.st-key-database-details{padding-left:0;border-left:0}}
</style>'''


def _projects(projects, key, *, default="", label="Project", optional=True, empty_label="None"):
    options = ([""] if optional else []) + [row["id"] for row in projects]
    if not options:
        st.caption("Create a project in Projects before starting a calculation.")
        return ""
    if st.session_state.get(key) not in options:
        st.session_state[key] = default if default in options else options[0]
    return st.selectbox(label, options, key=key, format_func=lambda value: next((row["name"] for row in projects if row["id"] == value), empty_label))


COLLECTION_SECTIONS = ["All collections", "Imported libraries", "Saved groups", "Generated results", "Job history"]
COLLECTION_SORTS = ["Libraries first", "Recently updated", "Name A–Z", "Largest first"]
COLLECTION_PAGE_SIZE = 15


def _collection_section(group):
    if group.get("import_method") == "user_uploads":
        return "Imported libraries"
    if group.get("collection_kind") == "generated":
        return "Generated results"
    if group.get("collection_kind") in {"curated", "reference"}:
        return "Saved groups"
    return "Job history"


def _browse_collections(groups, *, query="", section="All collections", sort="Libraries first", projects=()):
    project_names = {row["id"]: row["name"] for row in projects}
    terms = query.casefold().split()
    result = []
    for group in groups:
        category = _collection_section(group)
        if section != "All collections" and category != section:
            continue
        searchable = " ".join([
            group["name"], category, project_names.get(group.get("project_id"), ""),
            _remote_source_workstation(group.get("origin")),
        ]).casefold()
        if all(term in searchable for term in terms):
            result.append(group)
    result.sort(key=lambda row: (row["name"].casefold(), row["id"]))
    if sort == "Recently updated":
        result.sort(key=lambda row: row.get("updated_at") or "", reverse=True)
    elif sort == "Largest first":
        result.sort(key=lambda row: row.get("molecule_count", 0), reverse=True)
    elif sort == "Libraries first":
        result.sort(key=lambda row: row.get("molecule_count", 0), reverse=True)
        result.sort(key=lambda row: COLLECTION_SECTIONS.index(_collection_section(row)))
    return result


def _database_view_selector(current):
    selected = st.segmented_control("Database view", ["Molecules", "Runs", "Collections"], key="database-browser-view") or current
    if selected != current:
        st.rerun()


def _reset_import():
    st.session_state["molecular-import-open"] = False
    st.session_state["molecular-import-step"] = 0
    st.rerun()


def _dismiss_library_dialog():
    st.session_state.pop("data-browser-management", None)


def _dismiss_import_dialog():
    st.session_state["molecular-import-open"] = False


def _guess(columns, names):
    return next((column for column in columns if column.lower().replace(" ", "_") in names), "")


def _review_structure_chains(records, key):
    if len(records) != 1 or len(records[0]["components"]) <= 1:
        return records
    record = records[0]
    labels = {row["chain"]: f"{row['chain']} · {TYPE_LABELS[row['type']]} · {len(row['sequence'])} residues"
              if row.get("sequence") else f"{row['chain']} · {row.get('ccd') or 'Ligand'}"
              for row in record["components"]}
    source = (record.get("files") or [{}])[0]
    digest = hashlib.sha256(source.get("content", b"")).hexdigest()[:12]
    chains = st.multiselect("Chains to keep", list(labels), default=list(labels), format_func=labels.get, key=key + digest)
    if not chains:
        st.info("Choose at least one chain.")
        return []
    if set(chains) == set(labels):
        return records
    return [cropped_record(record, chains)]


def _import_content(method):
    records = []
    audit = None
    source_content = None
    if method == "Prepared structure":
        source = st.session_state.get("molecular-import-prepared")
        if not source:
            raise ValueError("Open the completed preparation result to choose an output.")
        records = imports.parse_structure(source["content"], source["name"])
        for record in records:
            record["metadata"]["preparation_run_id"] = source["run_id"]
            record["files"][0].update(label="Prepared structure", preparation_run_id=source["run_id"])
    elif method in {"Table", "Compound table"}:
        compounds = method == "Compound table"
        if not compounds:
            st.caption("Map protein, ligand, DNA or RNA columns to chains.")
        upload = st.file_uploader("CSV, TSV, or Excel file", type=["csv", "tsv", "xlsx", "xlsm", "xls", "txt"], key="molecular-import-table")
        if upload is None:
            return [], None, None
        upload_content = upload.getvalue()
        signature = hashlib.sha256(upload_content).hexdigest()[:12]
        sheets = imports.table_sheets(upload_content, upload.name)
        sheet_name = sheets[0] if sheets else ""
        if len(sheets) > 1:
            sheet_name = st.selectbox("Worksheet", list(sheets), key="molecular-import-sheet-" + signature)
        table = imports.read_table(upload_content, upload.name, sheet_name=sheet_name, compounds=compounds)
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
        if compounds:
            guessed = _guess(columns, {"smiles", "canonical_smiles", "isomeric_smiles", "structure"})
            smiles_column = st.selectbox("SMILES column", ["", *columns], index=columns.index(guessed) + 1 if guessed else 0,
                format_func=lambda value: value or "Choose a column", key="molecular-import-smiles-column-" + signature)
            if not smiles_column:
                return [], None, None
            audit = imports.analyze_small_molecule_table(upload_content, upload.name, smiles_column=smiles_column,
                name_column=name_column, id_column=id_column, sheet_name=sheet_name)
            return audit.pop("records"), audit, upload_content
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
    elif method in {"Structure files", "Protein files"}:
        protein = method == "Protein files"
        formats = ["pdb", "ent", "cif", "mmcif", "gz"]
        if not protein:
            formats.extend(["sdf", "mol", "zip"])
        uploads = st.file_uploader("PDB or mmCIF files" if protein else "PDB, mmCIF, SDF, MOL, or ZIP",
            type=formats, accept_multiple_files=True, key="molecular-import-structures-" + method)
        include_ligands = st.checkbox("Include bound ligands", value=not protein, key="molecular-import-include-ligands-" + method)
        for upload in uploads:
            parsed = imports.parse_archive(upload.getvalue()) if upload.name.lower().endswith(".zip") else imports.parse_structure(
                upload.getvalue(), upload.name, include_ligands=include_ligands)
            if protein and any(not any(component["type"] == "protein" for component in row["components"]) for row in parsed):
                raise ValueError(f"{upload.name} has no protein chain.")
            records.extend(parsed)
        records = _review_structure_chains(records, "molecular-import-upload-chains-")
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
    elif method in {"Fetch by identifier", "PDB accession", "PDB coordinates"}:
        pdb_only = method in {"PDB accession", "PDB coordinates"}
        cols = st.columns([3, 1]) if pdb_only else st.columns([1, 3, 1])
        source = "PDB" if pdb_only else cols[0].selectbox(
            "Database", ["PDB", "UniProt"], key="molecular-import-identifier-source")
        coordinate_format = "cif"
        if pdb_only:
            coordinate_format = st.radio("Coordinate format", ["pdb", "cif"], horizontal=True,
                format_func={"pdb": "PDB", "cif": "mmCIF"}.get, key="molecular-import-pdb-format")
        identifier = cols[-2].text_input("Accession", placeholder="1IEP" if source == "PDB" else "P04637", key="molecular-import-identifier")
        fetch_key = (source, identifier.strip(), coordinate_format)
        if cols[-1].button("Fetch", key="molecular-import-fetch", disabled=not identifier.strip()):
            filename, content = imports.fetch_identifier(identifier, source, structure_format=coordinate_format)
            st.session_state["molecular-import-fetched"] = {"key": fetch_key, "filename": filename, "content": content}
        fetched = st.session_state.get("molecular-import-fetched")
        if fetched and fetched["key"] == fetch_key:
            if source == "PDB":
                include_ligands = st.checkbox("Include bound ligands", value=method != "PDB accession", key="molecular-import-pdb-ligands")
                records = imports.parse_structure(fetched["content"], fetched["filename"], include_ligands=include_ligands)
                if method == "PDB accession" and any(not any(c["type"] == "protein" for c in row["components"]) for row in records):
                    raise ValueError("This PDB entry has no protein chain.")
                records = _review_structure_chains(records, "molecular-import-pdb-chains-")
            else:
                records = imports.parse_sequences(fetched["content"].decode())
    return records, audit, source_content


def _review_smiles_corrections(method, audit, records, source_content, instance):
    invalid = {row["row_number"]: row for row in audit["rows"] if row["status"] in {"invalid", "empty"}}
    if not invalid:
        return records, audit
    key = f"dataset-smiles-corrections-{instance}-{audit['batch_id']}"
    with st.expander("Correct invalid SMILES"):
        st.caption("Enter one row number and replacement SMILES per line, then validate the corrections.")
        st.dataframe([{"Row": number, "Input SMILES": row["raw_smiles"], "Candidate": row["annotation_stripped_candidate"]}
                      for number, row in list(invalid.items())[:100]], hide_index=True, width="stretch")
        text = st.text_area("Replacement SMILES", placeholder="12 CCO", key=key)
        if st.button("Validate corrections", key=key + "-validate"):
            corrections = {}
            for line in text.splitlines():
                if not line.strip():
                    continue
                parts = line.strip().split(maxsplit=1)
                if len(parts) != 2 or not parts[0].isdigit() or int(parts[0]) not in invalid:
                    raise ValueError("Enter an invalid or blank row number followed by the replacement SMILES.")
                number = int(parts[0])
                if number in corrections:
                    raise ValueError(f"Row {number} appears more than once in the corrections.")
                corrections[number] = parts[1]
            st.session_state[key + "-validated"] = corrections
        corrections = st.session_state.get(key + "-validated", {})
        if not corrections or source_content is None:
            return records, audit
        if method == "SMILES":
            revised = imports.analyze_smiles_text(source_content.decode("utf-8-sig"), corrections=corrections)
            revised["source_name"] = audit["source_name"]
            revised["source_sha256"] = audit["source_sha256"]
        else:
            revised = imports.analyze_small_molecule_table(source_content, audit["source_name"],
                                                           **audit["mapping"], corrections=corrections)
        return revised.pop("records"), revised


def _import_report_frame(rows):
    return [{
        "Row": row["row_number"],
        "Status": row["status"],
        "Source ID": row.get("source_id") or "",
        "Compound ID": row.get("compound_id") or "",
        "Name": row.get("name") or "",
        "Input SMILES": row.get("raw_smiles") or "",
        "Replacement SMILES": row.get("corrected_smiles") or "",
        "Canonical isomeric SMILES": row.get("canonical_smiles") or "",
        "Annotation candidate (review only)": row.get("annotation_stripped_candidate") or "",
        "Candidate status": row.get("candidate_status") or "",
        "Duplicate of row": row.get("duplicate_of_row") or "",
        "Stereo variant of row": row.get("stereo_variant_of_row") or "",
        "Destination match": "Yes" if row.get("matched_existing") else "",
        "Warnings": "; ".join(row.get("warnings") or []),
        "Errors": "; ".join(row.get("errors") or []),
        "Formula": (row.get("descriptors") or {}).get("formula", ""),
        "Molecular weight": (row.get("descriptors") or {}).get("molecular_weight"),
        "Exact mass": (row.get("descriptors") or {}).get("exact_mass", ""),
        "Heavy atoms": (row.get("descriptors") or {}).get("heavy_atoms"),
        "Formal charge": (row.get("descriptors") or {}).get("formal_charge"),
        "Fragments": (row.get("descriptors") or {}).get("fragment_count"),
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


@st.dialog("Import molecules", width="large", on_dismiss=_dismiss_import_dialog)
def render_import_dialog(user, projects, *, database=None, dataset_mode=False, default_reference_project=""):
    database = database if database is not None else MolecularDatabase()
    step = st.session_state.get("molecular-import-step", 0)
    kind = st.session_state.get("molecular-import-kind", "") if dataset_mode else ""
    methods = {
        "protein": ["PDB accession", "Protein files"],
        "structure": ["Structure files", "PDB coordinates"],
        "prepared": ["Prepared structure"],
        "compound": ["Compound table", "SMILES", "Draw molecule"],
        "sequence": ["FASTA", "Paste sequences"],
    }.get(kind, ["Table", "FASTA", "Structure files", "Paste sequences", "SMILES", "Draw molecule", "Fetch by identifier"])
    if kind == "compound":
        st.caption(f"SQL imports support up to {imports.MAX_COMPOUND_ROWS:,} compound rows per import.")
        if st.button("Open Library preparation for larger files", key="molecular-import-large-library"):
            from mn_molecular_design_studio.navigation import open_workflow_tool

            st.session_state.pop("molecular-import-open", None)
            open_workflow_tool(st.session_state, "library.prepare")
            st.rerun(scope="app")
    if step == 0 and dataset_mode:
        st.caption({"protein": "Import protein structures", "structure": "Import coordinate files",
                    "prepared": "Save a preparation output for reuse",
                    "compound": "Import a compound library",
                    "sequence": "Import sequences"}.get(kind, "Choose where to save your inputs."))
    if step == 0:
        groups = [row for row in database.groups(user, collection="library") if row["owner_id"] == user["id"] or user.get("role") == "admin"]
        if dataset_mode:
            from mn_molecular_design_studio.datasets import display_labels, is_dataset
            groups = [row for row in groups if is_dataset(row) and not row["read_only"]]
        labels = display_labels(groups, projects) if dataset_mode else {row["id"]: row["name"] for row in groups}
        if st.session_state.get("molecular-import-group") not in {"", *(row["id"] for row in groups)}:
            st.session_state["molecular-import-group"] = ""
        group_id = st.selectbox("Dataset *" if dataset_mode else "Group *", ["", *[row["id"] for row in groups]], format_func=lambda value: labels.get(value, "+ Create a new dataset" if dataset_mode else "+ Create a new group"), key="molecular-import-group")
        name = st.text_input("New dataset name" if dataset_mode else "New group name", placeholder="e.g. ALK2 screening compounds", key="molecular-import-group-name") if not group_id else next(row["name"] for row in groups if row["id"] == group_id)
        project = _projects(projects, "molecular-import-project", default=default_reference_project,
            label="Project reference (optional)" if dataset_mode else "Project",
            empty_label="No project — reusable library" if dataset_mode else "None") if not group_id else next(row["project_id"] for row in groups if row["id"] == group_id)
        left, right = st.columns([1, 1])
        if left.button("Cancel", key="molecular-import-cancel-0"):
            _reset_import()
        if right.button("Next", type="primary", disabled=not name.strip(), key="molecular-import-next"):
            st.session_state["molecular-import-destination"] = {"group_id": group_id or None, "name": name, "project_id": project}
            if dataset_mode:
                from mn_molecular_design_studio.datasets import input_origin
                st.session_state["molecular-import-destination"]["project_id"] = ""
                if not group_id:
                    try:
                        st.session_state["molecular-import-destination"]["origin"] = input_origin(user, [project] if project else [])
                    except (ValueError, PermissionError, OSError) as exc:
                        st.error(str(exc))
                        return
            st.session_state["molecular-import-instance"] = uuid.uuid4().hex
            st.session_state["molecular-import-step"] = 1
            st.rerun(scope="app")
        return
    destination = st.session_state["molecular-import-destination"]
    import_instance = st.session_state.setdefault("molecular-import-instance", uuid.uuid4().hex)
    st.markdown(f"**{html.escape(destination['name'])}**")
    method = st.radio("Source", methods, horizontal=True, key="molecular-import-method")
    records = []
    audit = None
    source_content = None
    review_confirmed = True
    block_duplicates = False
    try:
        records, audit, source_content = _import_content(method)
        if audit:
            records, audit = _review_smiles_corrections(method, audit, records, source_content, import_instance)
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
            if dataset_mode:
                from mn_molecular_design_studio.datasets import input_records
                records = input_records(records)
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
            if dataset_mode:
                from mn_molecular_design_studio.navigation import open_datasets
                open_datasets(st.session_state, group_id=result["group_id"])
                st.session_state.pop("molecular-preselect-group", None)
                st.session_state.pop("molecular-preselect-collection", None)
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


def render_details(database, user, record, *, properties_outcome=None):
    from mn_molecular_design_studio.evaluations import evaluation_descriptor
    design_context = record.get("_design_evaluations")
    outcome_context = record.get("_outcome_context")
    evaluation = record.get("selected_evaluation") or {}
    history = {evaluation_descriptor(obs)["id"]: evaluation_descriptor(obs) for obs in record.get("observations", [])}
    if design_context is not None:
        from mn_molecular_design_studio.database_browser import load_result_record
        from mn_molecular_design_studio.database_design_evaluations import evaluation_label

        linked = {row["evaluation_id"]: row for row in design_context["results"]}
        current = evaluation.get("id") or ""
        options = list(linked)
        if current not in linked:
            options.insert(0, "")
        if linked:
            with st.expander("Models & calculation history", expanded=outcome_context is None):
                choice = st.selectbox("Evaluation", options, index=options.index(current) if current in options else 0,
                    format_func=lambda value: evaluation_label(linked[value]) if value else "Selected original structure / record",
                    key=design_context["key"])
                st.caption("Each calculation keeps its own model, settings and scores. Choosing another calculation does not replace the original result scores.")
            if choice and choice != current:
                selected = load_result_record(database, user, linked[choice], **design_context["filters"])
                if selected is None:
                    st.warning("This evaluation is no longer accessible. Refresh or choose another evaluation.")
                    return
                record = selected
                evaluation = record.get("selected_evaluation") or {}
        if outcome_context is None:
            st.caption(f"{design_context['analysis_count']:,} connected post-analyses. Open Design analyses for their recorded scores.")
        if design_context["unresolved"]:
            st.caption(f"{design_context['unresolved']:,} evaluations have an unresolved target/structure link and are not attached here. They remain in result history.")
    elif history:
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
    st.markdown("**" + html.escape((properties_outcome or {}).get("display_name") or record["name"]) + "**")
    if outcome_context is None:
        st.caption(record["id"])
    sequence_reference = (record.get("metadata") or {}).get("sequence_source_reference") or {}
    if sequence_reference.get("pdb_id"):
        st.caption("Sequence reference: PDB " + str(sequence_reference["pdb_id"]) + " · "
                   + str(sequence_reference.get("matched_chain_count") or "") + " matched chain(s). Target roles and template use are recorded per run.")
    reference_id = (record.get("metadata") or {}).get("sequence_source_reference_molecule_id")
    if reference_id:
        references = database.list_molecules(user, molecule_ids=[reference_id], collection="all", limit=1)["items"]
        if references and st.button("Open source PDB in Protein Library", key="database-open-source-reference-" + record["id"]):
            st.session_state["molecular-preselect-collection"] = "library"
            st.session_state["molecular-preselect-group"] = references[0]["groups"][0]["id"]
            st.session_state["molecular-preselect-molecule"] = reference_id
            st.rerun()
    linked_target = (record.get("metadata") or {}).get("design_target_molecule_id")
    if linked_target and outcome_context is None:
        target_rows = database.list_molecules(user, molecule_ids=[linked_target], collection="all", limit=1)["items"]
        if target_rows:
            st.caption("Design target: " + target_rows[0]["name"] + " · " + linked_target)
    source_workstation = _remote_source_workstation(record.get("origin"))
    if source_workstation:
        st.caption("Synced from " + source_workstation)
    molecule_tab, scores_tab = st.tabs(["This molecule", "Evaluation details" if outcome_context is not None else
                                      f"Design analyses ({design_context['analysis_count']})" if design_context is not None else "Group scores"])
    with molecule_tab:
        structures = record["structures"]
        references = [item for item in record.get("representations") or []
                      if (item.get("context") or {}).get("structure_role") == "source_reference"]
        if references and evaluation:
            structure_view = st.segmented_control("Structure view", ["Evaluation structures", "Source references"],
                default="Evaluation structures", key="database-structure-view-" + record["id"])
            if structure_view == "Source references":
                structures = references
                st.caption("These are source reference coordinates. The selected Evaluation scores describe its prediction or analysis representation.")
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
                if len(record.get("components") or []) > 1 and path.suffix.lower() in {".pdb", ".ent", ".cif", ".mmcif", ".pdbqt", ".gz"}:
                    if st.button("Open in Complex Analysis", icon=":material/biotech:", key="database-open-complex-analysis-" + structure["id"]):
                        from mn_molecular_design_studio.complex_analysis_ui import open_database_complex
                        open_database_complex(
                            record["id"],
                            evaluation_id=str(evaluation.get("id") or ""),
                            structure_id=str(structure["id"]),
                        )
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
                    elif policy == "bound_docking_complex":
                        representations = [StructureVisualization(
                            pdb=_structure_text(path, structure_context), color="chain-id",
                            representation_type="cartoon")]
                    else:
                        # Ordinary protein structures stay in cartoon form; the
                        # specialized pose path above splits a ligand complex so
                        # multiple PDBQT poses cannot pile up in the viewer.
                        representations = [StructureVisualization(pdb=_structure_text(path, structure_context), color="chain-id", representation_type="ball-and-stick" if record["type"] == "small_molecule" else "cartoon")]
                    molstar_custom_component(representations, key="database-molstar-native-" + structure["id"], height=310, show_controls=False, show_structure_tools=False)
                from mn_molecular_design_studio.bound_docking import is_bound_pose
                if is_bound_pose(structure):
                    from mn_molecular_design_studio.workspace import load_organization
                    from mn_molecular_design_studio.bound_docking import register_bound_complex
                    source_group = database.group(structure["group_id"], user)
                    can_add_complex = not source_group.get("read_only") and (
                        user.get("role") == "admin" or source_group.get("owner_id") == user.get("id"))
                    existing_complex = next((item for item in (record.get("representations") or record.get("structures") or [])
                        if (item.get("context") or {}).get("database_structure_policy") == "bound_docking_complex"
                        and str((item.get("context") or {}).get("source_pose_structure_id") or "") == str(structure["id"])), None)
                    has_one_ligand_identity = record["type"] == "small_molecule" and len([
                        item for item in record.get("components") or []
                        if str(item.get("type") or "").casefold() == "small_molecule" and item.get("smiles")
                    ]) == 1
                    if has_one_ligand_identity:
                        if existing_complex:
                            if st.button("Open bound complex in Complex Analysis", key="database-bound-analysis-" + structure["id"]):
                                from mn_molecular_design_studio.complex_analysis_ui import open_database_complex
                                open_database_complex(record["id"], structure_id=existing_complex["id"])
                        elif st.button("Create bound complex for analysis", disabled=not can_add_complex,
                                       key="database-bound-analysis-create-" + structure["id"]):
                            try:
                                result = register_bound_complex(database, user, record["id"], structure_id=structure["id"])
                                from mn_molecular_design_studio.complex_analysis_ui import open_database_complex
                                open_database_complex(record["id"], structure_id=result["complex_structure_id"])
                            except (OSError, PermissionError, RuntimeError, TypeError, ValueError) as exc:
                                st.error(f"Could not create an analysis-ready complex: {exc}")
                        elif not can_add_complex:
                            st.caption("Only the source Database group owner or an administrator can attach a reconstructed complex to this record.")
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
                elif (structure.get("context") or {}).get("database_structure_policy") == "bound_docking_complex":
                    if st.button("Open in Complex Analysis", key="database-bound-analysis-open-" + structure["id"]):
                        from mn_molecular_design_studio.complex_analysis_ui import open_database_complex
                        open_database_complex(record["id"], structure_id=structure["id"])
            except (OSError, ValueError, PermissionError) as exc:
                st.caption(str(exc))
        elif source_workstation:
            st.info(f"No structure file is attached to this synced Database record. Its source is {source_workstation}; the structure must be indexed on that workstation before it can be opened here.")
        if outcome_context is not None:
            st.caption(outcome_context["kind"] + " · " + outcome_context["eligibility"]
                       + " · " + outcome_context["acceptance_level"] + f" · {outcome_context['source_runs']:,} originating runs")
            st.markdown("#### Original result scores")
            _render_original_result_scores(outcome_context)
            st.markdown("#### Downstream evaluations")
            inventory = outcome_context.get("analysis_inventory", {})
            tools = list(dict.fromkeys(["ipSAE", "Biotite", "PEP-Patch", *sorted(inventory)]))
            st.dataframe([{"Analysis": label, "Status": analysis_status(inventory, label),
                "Other models": inventory.get(label, {}).get("other_structure", 0),
                "Failed / unavailable": inventory.get(label, {}).get("failed", 0),
                "Unknown status": inventory.get(label, {}).get("unknown", 0)} for label in tools],
                hide_index=True, use_container_width=True)
            st.caption("— means no verified linked result in SQL, not proof that the analysis never ran. Other models keep separate scores.")
            if design_context is not None:
                with st.expander("Analysis scores & provenance"):
                    _render_design_analysis_scores(design_context, compact=True)
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
        from mn_molecular_design_studio.entity_properties_ui import render_record_properties
        render_record_properties(database, user, record, outcome=properties_outcome, key_prefix="database-record")
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
                        st.session_state["selected-job-reference"] = source.get("job_reference") or source["task_group"] + "/" + source["run_id"]
                        st.session_state["page"] = "Job details"
                        st.rerun()
                if source.get("parents"):
                    render_result_details(st, source["parents"])
    with scores_tab:
        if design_context is not None:
            if outcome_context is None:
                _render_design_analysis_scores(design_context)
            st.markdown("**Selected evaluation scores**")
        if not record["observations"]:
            st.caption("No scores recorded yet.")
        else:
            st.caption("These scores belong only to the selected evaluation." if design_context is not None else
                       "The sheet uses one selected evaluation. All recorded runs remain available here.")
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
                render_result_details(st, record["observations"][observation]["context"])


def _render_original_result_scores(context):
    st.caption("Each producing-tool evaluation retains its original values, model and settings. Different runs are not averaged or overwritten by an analysis.")
    if context["original_count"] > len(context["originals"]):
        st.caption("Showing up to 20 source evaluations. Use Original scores page in the left filters for more.")
    for result in context["originals"]:
        with st.expander(result["tool"] + " · " + result["run_id"] + " · original result",
                         expanded=len(context["originals"]) <= 3):
            st.caption("Source filter decision: " + result.get("receipt", {}).get("state", "Needs review") + " · "
                       + result.get("receipt", {}).get("acceptance_level", "Not recorded"))
            _render_evaluation_measurements(result)


def _render_design_analysis_scores(context, *, compact=False):
    count = context["analysis_count"] + context.get("validation_count", 0)
    if not count:
        st.caption("No verified post-analysis results are linked to this design in the selected target context.")
        return
    st.caption("Recorded analyses of this exact design. Other structures and sequence-only checks are labelled separately.")
    if count > len(context["analyses"]):
        st.caption("Showing up to 20 analyses. Use Analysis scores page in the left filters for more.")
    for result in context["analyses"]:
        title = tool_label(result["tool"]) + " · " + result["structure_relation"]
        if not compact:
            title += " · " + result["run_id"]
        if result.get("receipt", {}).get("role") == "design_validation" or result.get("receipt", {}).get("design_refold"):
            title += " · Design refolding / validation"
        with st.expander(title, expanded=not compact and len(context["analyses"]) <= 3):
            _render_evaluation_measurements(result)


def _render_evaluation_measurements(result, *, compact=False):
    evaluation = result.get("evaluation") or {}
    st.caption("Model / pose: " + (result["model"] or "not recorded") + " · version: "
               + (evaluation.get("engine_version") or evaluation.get("image_digest") or "not recorded"))
    if result["metrics"]:
        st.dataframe([{"Metric": result_field_label(key), "Value": result_value_text(value, summarize_fields=True)}
                      for key, value in sorted(result["metrics"].items())], hide_index=True, use_container_width=True,
                     column_order=["Value", "Metric"] if compact else None,
                     column_config={"Value": st.column_config.TextColumn(width="small")} if compact else None)
    else:
        st.caption("No scores were recorded for this source. Missing measurements are not zero.")
    with st.popover("Scientific context", use_container_width=True):
        st.caption("Source structure: " + (evaluation.get("representation_id") or evaluation.get("source_structure_digest") or "not recorded"))
        st.caption("Calculation digest: " + (evaluation.get("calculation_structure_digest") or "not recorded")
                   + " · target: " + (evaluation.get("target_id") or "not recorded"))
        st.caption("Target revision: " + (evaluation.get("target_revision") or "not recorded"))
        st.caption("Scientific settings: " + (evaluation.get("parameters_digest") or "not recorded"))
        st.caption("Analysed chains — binder: " + (", ".join(evaluation.get("binder_chains") or []) or "not recorded")
                   + " · target: " + (", ".join(evaluation.get("target_chains") or []) or "not recorded"))
        st.caption("Execution status: " + (result["status"] or "not recorded") + ". This is not a filter acceptance decision.")


@st.dialog("Create molecule group", width="small", on_dismiss=_dismiss_library_dialog)
def render_group_dialog(user, projects, ids, default_project="", *, database=None, dataset_mode=False):
    database = database if database is not None else MolecularDatabase()
    if ids:
        st.caption(f"{len(ids):,} selected molecules")
    name = st.text_input("Dataset name" if dataset_mode else "Group name", key="molecular-group-name")
    project = _projects(projects, "molecular-group-project", default=default_project,
        label="Project reference (optional)" if dataset_mode else "Project")
    if st.button("Create dataset" if dataset_mode else "Create group", type="primary", disabled=not name.strip(), key="molecular-group-save"):
        try:
            if dataset_mode:
                from mn_molecular_design_studio.datasets import create_dataset
                group_id = create_dataset(database, user, name, molecule_ids=ids, project_ids=[project] if project else [])["id"]
            else:
                group_id = database.create_group(user, name, project_id=project)
                if ids:
                    database.add_to_group(user, group_id, ids)
        except (ValueError, PermissionError, OSError, sqlite3.Error) as exc:
            st.error(str(exc))
            return
        st.session_state["molecular-preselect-group"] = group_id
        st.session_state["molecular-preselect-collection"] = "library"
        st.session_state["database-group-request"] = None
        _dismiss_library_dialog()
        if dataset_mode:
            from mn_molecular_design_studio.navigation import open_datasets
            open_datasets(st.session_state, group_id=group_id)
            st.session_state.pop("molecular-preselect-group", None)
            st.session_state.pop("molecular-preselect-collection", None)
        st.rerun()


@st.dialog("Add column", width="small", on_dismiss=_dismiss_library_dialog)
def render_column_dialog(user, group_id, *, database=None):
    database = database if database is not None else MolecularDatabase()
    try:
        group = database.group(group_id, user, write=True)
    except (ValueError, PermissionError, OSError, sqlite3.Error) as exc:
        st.error(str(exc))
        return
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
            _dismiss_library_dialog()
            st.rerun()
        except (ValueError, PermissionError) as exc:
            st.error(str(exc))


def render_database(user, projects):
    from mn_molecular_design_studio.entity_properties_ui import render_notice
    render_notice()
    database = MolecularDatabase()
    errors = st.session_state.pop("molecular-indexing-errors", [])
    if workspace_sync_running(user, database=database):
        st.caption("Showing saved Database records. Job indexing is updating in the background.")
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
        st.session_state["molecular-collection-section"] = "All collections"
        st.session_state["database-selected-group"] = pending_group
        st.session_state["database-browser-view"] = "Collections"
    if pending_molecule:
        st.session_state["database-focus"] = pending_molecule
    if message := st.session_state.pop("molecular-database-flash", None):
        st.success(message)
    if errors:
        with st.expander(f"{len(errors)} source jobs could not be indexed"):
            st.dataframe(errors, hide_index=True)
    if st.session_state.get("database-browser-view") not in {"Molecules", "Runs", "Collections"}:
        st.session_state["database-browser-view"] = "Molecules"
    browser_view = st.session_state["database-browser-view"]
    if browser_view != "Collections":
        from mn_molecular_design_studio.database_browser_ui import render_browser
        render_browser(database, user, projects, browser_view)
        return
    with st.container(key="molecular-database"):
        sidebar, main, details = st.columns([.25, .50, .25], gap="medium")
        with sidebar, st.container(key="database-groups"):
            st.markdown("### Collections")
            _database_view_selector("Collections")
            new, refresh = st.columns([4, 1])
            if new.button("Import / new", icon=":material/add:", use_container_width=True, key="molecular-database-new"):
                st.session_state["molecular-import-open"] = True
                st.session_state["molecular-import-step"] = 0
            if refresh.button("↻", key="molecular-database-refresh", help="Refresh molecules and job results"):
                st.session_state["molecular-indexing-errors"] = sync_workspace(user, database=database, force=True)
                st.rerun()
            group_query = st.text_input("Find a collection", placeholder="Find a library, campaign, project…", label_visibility="collapsed", key="molecular-group-search")
            molecule_controls = st.container()
            with st.popover("Access & source filters", icon=":material/tune:", use_container_width=True):
                collection_label = st.segmented_control(
                    "Collection", ["Library", "All indexed"], default="Library", key="molecular-collection-view",
                    help="Library includes eligible imported and selected collections. All indexed also shows historical/unreviewed screening shortlists, evaluations, and intermediate job collections.",
                ) or "Library"
                scope_options = ["Mine", "Shared", "All"]
                old_scope = st.session_state.get("molecular-scope")
                scope_aliases = {"My groups": "Mine", "Organization": "All", "organization": "All", "mine": "Mine", "shared": "Shared", "all": "All"}
                if old_scope in scope_aliases:
                    st.session_state["molecular-scope"] = scope_aliases[old_scope]
                if st.session_state.get("molecular-scope") not in scope_options:
                    st.session_state["molecular-scope"] = "All"
                scope_label = st.segmented_control(
                    "Collection access", scope_options, default="All", key="molecular-scope",
                    help="Mine shows personal collections. Shared shows project collections. All combines both.",
                ) or "All"
                project_id = _projects(projects, "molecular-project-filter", label="Project filter", empty_label="All projects")
                hide_smoke_test = st.checkbox(
                    "Hide confirmed development runs", value=True, key="molecular-hide-smoke-test-groups",
                    help="Hides only the eight audited owner/run pairs. Names and development-account flags alone do not hide campaigns. Turn this off to inspect those runs.",
                )
                st.markdown("**Import method**")
                methods = [value for value, label in IMPORT_METHODS.items() if st.checkbox(label, value=True, key="molecular-method-" + value)]
            collection = "library" if collection_label == "Library" else "all"
            scope = {"Mine": "mine", "Shared": "shared", "All": "all"}[scope_label]
            project_name = next((row["name"] for row in projects if row["id"] == project_id), "All projects")
            st.caption(f"{collection_label} · {scope_label} access · {project_name}")
            if collection == "all":
                st.caption("Includes historical/unreviewed screening shortlists and intermediate job collections.")
            if not methods:
                methods = ["__none__"]
            groups = database.groups(
                user, project_id=project_id, scope=scope, methods=methods, collection=collection,
                hide_smoke_test=hide_smoke_test,
            )
            section_counts = {section: sum(_collection_section(row) == section for row in groups) for section in COLLECTION_SECTIONS[1:]}
            section = st.selectbox(
                "Browse collections", COLLECTION_SECTIONS, key="molecular-collection-section",
                format_func=lambda value: f"{value} ({len(groups) if value == 'All collections' else section_counts[value]:,})",
            )
            group_sort = st.selectbox("Collection order", COLLECTION_SORTS, key="molecular-collection-order", label_visibility="collapsed")
            filtered_groups = _browse_collections(groups, query=group_query, section=section, sort=group_sort, projects=projects)
            filtered_ids = {row["id"] for row in filtered_groups}
            if st.session_state.get("database-selected-group") not in filtered_ids:
                st.session_state["database-selected-group"] = ""
            browse_signature = json.dumps([group_query, section, group_sort, [row["id"] for row in filtered_groups]])
            page_count = max(1, math.ceil(len(filtered_groups) / COLLECTION_PAGE_SIZE))
            if st.session_state.get("molecular-browse-signature") != browse_signature:
                selected_index = next((index for index, row in enumerate(filtered_groups) if row["id"] == st.session_state.get("database-selected-group")), 0)
                st.session_state["molecular-collection-page"] = selected_index // COLLECTION_PAGE_SIZE + 1
                st.session_state["molecular-browse-signature"] = browse_signature
            if pending_group and pending_group in filtered_ids:
                selected_index = next(index for index, row in enumerate(filtered_groups) if row["id"] == pending_group)
                st.session_state["molecular-collection-page"] = selected_index // COLLECTION_PAGE_SIZE + 1
            group_page = int(st.session_state.get("molecular-collection-page", 1))
            page_groups = filtered_groups[(group_page - 1) * COLLECTION_PAGE_SIZE:group_page * COLLECTION_PAGE_SIZE]
            options = ["", *[row["id"] for row in page_groups]]
            selected_group = st.session_state.get("database-selected-group")
            if selected_group and selected_group not in options:
                options.insert(1, selected_group)
                st.caption("Your selected collection stays visible while browsing other pages.")
            group_lookup = {row["id"]: row for row in filtered_groups}
            st.caption(f"{len(filtered_groups):,} collections" + (f" · filtered from {len(groups):,}" if len(filtered_groups) != len(groups) else ""))
            if not filtered_groups:
                st.info("No matching collections. Try a shorter search, another category, or the access filters.")
            with st.container(height=490, border=False):
                def group_label(value):
                    row = group_lookup.get(value)
                    if row is None:
                        return "All molecules"
                    label = f"{row['name']} · {row['molecule_count']:,}"
                    if row.get("collection_kind") == "reference":
                        project_label = next((project["name"] for project in projects if project["id"] == row["project_id"]),
                                             row["project_id"][:8] if row["project_id"] else "Personal")
                        label += " · " + project_label
                    source = _remote_source_workstation(row.get("origin"))
                    return label + (f" · from {source}" if source else "")
                group_id = st.radio("Molecule groups", options, format_func=group_label, key="database-selected-group", label_visibility="collapsed")
            if page_count > 1:
                st.number_input("Collection page", min_value=1, max_value=page_count, step=1, key="molecular-collection-page")
                st.caption(f"Page {group_page} of {page_count} · {COLLECTION_PAGE_SIZE} collections per page")
            st.caption("Category and search narrow this list only. All molecules includes every collection allowed by the access filters.")
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
            st.markdown("### Molecules")
            file_col, edit_col, view_col, export_col = st.columns(4)
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
            st.markdown(f'<div class="database-collection-heading">{html.escape(group_name)}</div>', unsafe_allow_html=True)
            if active_group:
                from mn_molecular_design_studio.entity_properties_ui import render_folder_properties
                render_folder_properties(database, user, {}, {"id": "collection:" + active_group["id"], "label": active_group["name"],
                    "entity": {"kind": "dataset", "group_id": active_group["id"]}}, projects, key_prefix="database-collection")
                badges = [_collection_section(active_group), f"{active_group['molecule_count']:,} molecules",
                          "Project collection" if active_group.get("project_id") else "Personal collection"]
                if source := _remote_source_workstation(active_group.get("origin")):
                    badges.append("From " + source + " · read-only")
                st.markdown('<div class="database-collection-meta">' + "".join(f"<span>{html.escape(value)}</span>" for value in badges) + '</div>', unsafe_allow_html=True)
            else:
                st.caption("Browse all accessible molecules, or choose a collection on the left. Imported libraries appear first by default.")
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
                        render_result_details(st, summary)
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
            with molecule_controls:
                query = st.text_input("Search molecules", placeholder="Search name, sequence, SMILES, tag…", label_visibility="collapsed", key="molecular-query")
                with st.popover("Molecule & score filters", use_container_width=True):
                    sort_options = ["added", "name", *metric_fields]
                    if st.session_state.get("molecular-sort") not in sort_options:
                        st.session_state["molecular-sort"] = "added"
                    sort = st.selectbox("Sort", sort_options, format_func=lambda value: value.title() if value in {"added", "name"} else field_label(value), key="molecular-sort")
                    descending = st.toggle("Descending", value=True, key="molecular-desc")
                    kind = st.selectbox("Molecule type", ["", *TYPE_LABELS], format_func=lambda value: TYPE_LABELS.get(value, "All types"), key="molecular-type-filter")
                    filters = []
                    if not set(st.session_state.get("molecular-score-filters") or []).issubset(metric_fields):
                        st.session_state["molecular-score-filters"] = []
                    score_fields = st.multiselect("Score filters", metric_fields, format_func=field_label, key="molecular-score-filters")
                    for field in score_fields:
                        left, right = st.columns([1, 2])
                        operator = left.selectbox(field + " comparison", [">=", "<=", ">", "<", "="], key="molecular-op-" + field, label_visibility="collapsed")
                        value = right.number_input(field_label(field), value=0.0, key="molecular-value-" + field)
                        filters.append({"field": field, "operator": operator, "value": value})
                    st.markdown("**Evaluation context**")
                    for field, label in (("target_id", "Target"), ("evaluator", "Tool"), ("run_id", "Run")):
                        options = ["", *evaluation_options[field]]
                        key = "database-context-" + field + "-" + (group_id or "all")
                        if st.session_state.get(key, "") not in options:
                            st.session_state[key] = ""
                        value = st.selectbox(label, options, format_func=lambda item: item or "Latest matching result", key=key)
                        if value:
                            evaluation_filters[field] = value
                    st.caption("Each row shows scores from one matching result, never a mixture of different models or settings. Other model results remain available on Results.")
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
                    selected_group
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
                        except (ValueError, PermissionError, OSError, RuntimeError, sqlite3.Error) as exc:
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


def render_molecule_picker(user, *, key, allowed_types=None, multiple=False, lazy=False):
    library_label = "Protein Library" if allowed_types and set(allowed_types) <= {"protein"} else "Molecular Library"
    popover = st.popover(library_label, icon=":material/storage:",
                         key=key + "-library-open" if lazy else None,
                         on_change="rerun" if lazy else "ignore")
    with popover:
        if lazy and not popover.open:
            return None
        database = MolecularDatabase()
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
