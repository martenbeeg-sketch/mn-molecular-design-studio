"""Reusable target structure input and residue-selection components.

The UI is deliberately workflow-agnostic: a workflow supplies a stable key
prefix and receives a typed :class:`TargetDesignInput` describing the target.
This lets primary, secondary, and future workflow targets share the same
upload, chain, trim, hotspot, and coldspot behavior without sharing state.
"""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import streamlit as st

from mn_molecular_design_studio.engine.workflow_engine.app.components.molstar_viewer import (
    StructureVisualization,
    molstar_custom_component,
)
from mn_molecular_design_studio.engine.workflow_engine.core.structures import (
    filter_pdb_text,
    mmcif_protein_residues_by_chain,
    pdb_summary,
)
from mn_molecular_design_studio.auth import user_input_root


PDB_RESIDUE_CODES = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O", "ASX": "B", "GLX": "Z",
}


@dataclass
class TargetSource:
    """The uploaded or fetched structure file selected in a target source panel."""

    path: Path | None = None
    name: str = ""
    owner_id: str = "local"


@dataclass
class TargetDesignInput:
    """All structure and residue choices for one workflow target."""

    label: str
    key_prefix: str = "target"
    source_path: Path | None = None
    design_path: Path | None = None
    chains: list[str] = field(default_factory=list)
    hotspots: str = ""
    coldspots: str = ""
    viewer_mode: str = "Target only"
    source_text: str = ""
    design_text: str = ""
    sequence_by_chain: dict[str, list[tuple[int, str]]] = field(default_factory=dict)
    valid_residues: set[tuple[str, int]] = field(default_factory=set)
    valid_trim_ranges: bool = True
    hotspot_state_key: str = ""
    coldspot_state_key: str = ""

    @property
    def ready(self) -> bool:
        return self.source_path is not None and bool(self.chains) and self.design_path is not None and self.valid_trim_ranges


def _hash(value: str, length: int = 10) -> str:
    return hashlib.sha1(value.encode()).hexdigest()[:length]


def _target_storage_path(content: bytes, filename: str, *, owner_id: str, allow_non_protein: bool = False) -> Path:
    text = content.decode("utf-8", errors="replace")
    source_suffix = Path(filename).suffix.lower()
    if allow_non_protein:
        import io
        from Bio.PDB import PDBParser, MMCIFParser
        from mn_molecular_design_studio.engine.workflow_engine.core.candidates import _mmcif_text_with_default_occupancy
        try:
            parser = MMCIFParser(QUIET=True) if source_suffix in {".cif", ".mmcif"} else PDBParser(QUIET=True)
            if source_suffix in {".cif", ".mmcif"}:
                text = _mmcif_text_with_default_occupancy(text)
            structure = parser.get_structure("input", io.StringIO(text))
            if not any(True for _ in structure.get_atoms()):
                raise ValueError("No coordinate atoms were found.")
        except Exception as exc:
            raise ValueError(f"Could not read this coordinate structure: {exc}") from exc
    elif source_suffix in {".cif", ".mmcif"}:
        try:
            readable_chains = _mmcif_sequence_entries_by_chain_text(text)
        except Exception as exc:
            raise ValueError(f"Could not read this mmCIF structure with Gemmi: {exc}") from exc
        if not readable_chains:
            raise ValueError("This file does not contain a readable protein structure in mmCIF format.")
    elif not pdb_summary(text).get("chains"):
        raise ValueError("This file does not contain a readable protein structure in PDB format.")
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", Path(filename).stem).strip("._") or "target"
    digest = hashlib.sha256(content).hexdigest()[:12]
    suffix = source_suffix if source_suffix in {".cif", ".mmcif"} else ".pdb"
    destination = user_input_root(owner_id) / "structures" / f"{stem}-{digest}{suffix}"
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        destination.write_bytes(content)
    return destination


def _derived_target_storage_path(pdb_text: str, filename: str, *, owner_id: str) -> Path:
    content = pdb_text.encode("utf-8")
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", Path(filename).stem).strip("._") or "target"
    digest = hashlib.sha256(content).hexdigest()[:12]
    destination = user_input_root(owner_id) / "structures" / f"{stem}-{digest}.pdb"
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not destination.exists():
        destination.write_bytes(content)
    return destination


def _rcsb_mmcif_to_pdb(content: bytes, *, allow_non_protein: bool = False) -> bytes:
    """Convert an RCSB mmCIF-only entry to PDB, remapping long chain IDs."""
    import gemmi

    document = gemmi.cif.read_string(content.decode("utf-8", errors="replace"))
    if not len(document):
        raise ValueError("RCSB returned an empty mmCIF structure.")
    structure = gemmi.make_structure_from_block(document.sole_block())
    if not len(structure):
        raise ValueError("The RCSB mmCIF file does not contain a readable coordinate model.")

    chain_names = list(dict.fromkeys(chain.name for model in structure for chain in model))
    pdb_chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    if len(chain_names) > len(pdb_chain_ids):
        raise ValueError("The structure has more chains than the PDB format can represent.")

    # Preserve existing one-character IDs and give longer mmCIF author IDs the
    # next free PDB chain IDs. AF2's template input reads chain IDs from PDB.
    mapping = {
        name: name
        for name in chain_names
        if len(name) == 1 and name in pdb_chain_ids
    }
    available_ids = iter(chain_id for chain_id in pdb_chain_ids if chain_id not in mapping.values())
    for name in chain_names:
        if name not in mapping:
            mapping[name] = next(available_ids)

    for model in structure:
        for chain in model:
            chain.name = mapping[chain.name]
    pdb_text = structure.make_pdb_string()
    records = ("ATOM  ", "HETATM") if allow_non_protein else ("ATOM  ",)
    if not any(line.startswith(records) for line in pdb_text.splitlines()):
        raise ValueError("The RCSB mmCIF file has no protein atom records to convert.")
    return pdb_text.encode("utf-8")


def _protein_residue_ids_by_chain(pdb_text: str) -> dict[str, list[tuple[int, str]]]:
    residues: dict[str, list[tuple[int, str]]] = {}
    seen: set[tuple[str, int, str]] = set()
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 27:
            continue
        if line[17:20].strip().upper() not in PDB_RESIDUE_CODES:
            continue
        chain = line[21].strip() or "_"
        try:
            number = int(line[22:26])
        except ValueError:
            continue
        identity = (chain, number, line[26].strip())
        if identity in seen:
            continue
        seen.add(identity)
        residues.setdefault(chain, []).append((number, identity[2]))
    return residues


def _trim_target_chain_ranges(pdb_text: str, retained_ranges: dict[str, tuple[int, int]]) -> str:
    """Keep the selected PDB residue range for each trimmed target chain."""
    residues_by_chain = _protein_residue_ids_by_chain(pdb_text)
    retained_by_chain: dict[str, set[tuple[int, str]]] = {}
    for chain, (start, end) in retained_ranges.items():
        residues = residues_by_chain.get(chain, [])
        if not residues or (start <= residues[0][0] and end >= residues[-1][0]):
            continue
        retained_by_chain[chain] = {
            residue_id for residue_id in residues if start <= residue_id[0] <= end
        }

    output: list[str] = []
    for line in pdb_text.splitlines():
        record = line[:6].strip()
        chain = (line[21].strip() or "_") if len(line) > 21 else "_"
        if record == "TER" and chain in retained_by_chain:
            continue
        if record in {"ATOM", "HETATM", "ANISOU"} and chain in retained_by_chain and len(line) >= 27:
            try:
                residue_id = (int(line[22:26]), line[26].strip())
            except ValueError:
                residue_id = None
            if residue_id is not None and residue_id not in retained_by_chain[chain]:
                continue
        output.append(line)

    if output and output[-1][:6].strip() != "END":
        output.append("END")
    return "\n".join(output) + ("\n" if output else "")


def _residue_map(pdb_text: str, structure_path: Path | None = None) -> dict[str, list[int]]:
    if structure_path and structure_path.suffix.lower() in {".cif", ".mmcif"}:
        return {
            chain: [number for number, _residue_name in entries]
            for chain, entries in _mmcif_sequence_entries_by_chain(structure_path).items()
        }
    residues: dict[str, list[int]] = {}
    seen: set[tuple[str, int, str]] = set()
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 27:
            continue
        chain = line[21].strip() or "_"
        try:
            number = int(line[22:26])
        except ValueError:
            continue
        identity = (chain, number, line[26].strip())
        if identity in seen:
            continue
        seen.add(identity)
        residues.setdefault(chain, []).append(number)
    return residues


def _selected_viewer_residues(value: object, residue_map: dict[str, list[int]]) -> set[tuple[str, int]]:
    if isinstance(value, str):
        try:
            value = json.loads(value)
        except json.JSONDecodeError:
            return set()
    if not isinstance(value, dict):
        return set()
    selected: set[tuple[str, int]] = set()
    for selection in value.get("sequenceSelections") or []:
        if not isinstance(selection, dict):
            continue
        chain = str(selection.get("chainId") or "").strip()
        chain_residues = residue_map.get(chain, [])
        available = set(chain_residues)
        for raw_number in selection.get("residues") or []:
            try:
                number = int(raw_number)
            except (TypeError, ValueError):
                continue
            if number in available:
                selected.add((chain, number))
            elif 1 <= number <= len(chain_residues):
                selected.add((chain, chain_residues[number - 1]))
    return selected


def _hotspot_text(residues: set[tuple[str, int]]) -> str:
    return ",".join(f"{chain}{number}" for chain, number in sorted(residues))


def _increment_session_counter(key: str) -> None:
    st.session_state[key] = int(st.session_state.get(key, 0)) + 1


def _hotspot_pairs(text: str) -> set[tuple[str, int]]:
    residues: set[tuple[str, int]] = set()
    for token in re.split(r"[,;\s]+", text.strip()):
        selected_range = re.fullmatch(r"([A-Za-z_][A-Za-z0-9_.-]*?)([0-9]+)-([0-9]+)", token)
        if selected_range:
            start,end=int(selected_range[2]),int(selected_range[3])
            if 0<=end-start<=10000:
                residues.update((selected_range[1],number) for number in range(start,end+1))
            continue
        match = re.fullmatch(r"([A-Za-z_][A-Za-z0-9_.-]*?)([0-9]+)", token)
        if match:
            residues.add((match[1], int(match[2])))
    return residues


def _sequence_entries_by_chain(pdb_text: str) -> dict[str, list[tuple[int, str]]]:
    sequences: dict[str, list[tuple[int, str]]] = {}
    seen: set[tuple[str, int, str]] = set()
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 27:
            continue
        one_letter = PDB_RESIDUE_CODES.get(line[17:20].strip().upper())
        if not one_letter:
            continue
        chain = line[21].strip() or "_"
        try:
            number = int(line[22:26])
        except ValueError:
            continue
        identity = (chain, number, line[26].strip())
        if identity in seen:
            continue
        seen.add(identity)
        sequences.setdefault(chain, []).append((number, one_letter))
    return sequences


def _mmcif_sequence_entries_by_chain_text(cif_text: str) -> dict[str, list[tuple[int, str]]]:
    """Read protein chain IDs and author residue numbers from mmCIF text."""
    return {
        chain: [(number, residue_name) for number, _insertion, residue_name in entries]
        for chain, entries in mmcif_protein_residues_by_chain(cif_text).items()
    }


def _mmcif_sequence_entries_by_chain(path: Path) -> dict[str, list[tuple[int, str]]]:
    return {
        chain: [(number, residue_name) for number, _insertion, residue_name in entries]
        for chain, entries in mmcif_protein_residues_by_chain(path).items()
    }


def render_target_source(
    *,
    key_prefix: str,
    label: str,
    example_code: str = "7JZU",
    fetch_button_label: str = "Fetch PDB",
    fetch_button_key: str | None = None,
    upload_label: str = "Drop a PDB file here or browse",
    upload_help: str = "PDB files only. CIF conversion is not part of this workflow yet.",
    allow_cif: bool = False,
    allow_non_protein: bool = False,
    fetch_enabled: bool = True,
    upload_file_types: list[str] | None = None,
    upload_expander_label: str | None = None,
    upload_in_expander: bool = True,
    fetched_path_key: str | None = None,
    fetched_name_key: str | None = None,
    owner_id: str = "local",
) -> TargetSource:
    """Render reusable PDB accession and structure upload controls."""
    pdb_code = ""
    fetch_clicked = False
    if fetch_enabled:
        code_col, fetch_col = st.columns([1.25, .65])
        with code_col:
            pdb_code = st.text_input(
                f"{label} PDB code",
                placeholder=f"Enter PDB code, e.g. {example_code}",
                key=f"{key_prefix}-pdb-code",
            ).strip().upper()
        with fetch_col:
            st.markdown("<div style='height:1.9rem'></div>", unsafe_allow_html=True)
            fetch_clicked = st.button(
                fetch_button_label,
                key=fetch_button_key or f"{key_prefix}-fetch-pdb",
                use_container_width=True,
            )

    stored_path_key = fetched_path_key or f"{key_prefix}-fetched-target"
    stored_name_key = fetched_name_key or f"{key_prefix}-fetched-target-name"
    converted_mmcif_key = f"{stored_path_key}-converted-from-mmcif"
    upload_key = f"{key_prefix}-target-upload"
    if fetch_clicked:
        if not re.fullmatch(r"[A-Z0-9]{4}", pdb_code):
            st.error(f"Enter a four-character PDB accession, such as {example_code}.")
        else:
            try:
                request = Request(
                    f"https://files.rcsb.org/download/{pdb_code}.pdb",
                    headers={"User-Agent": "MN-Molecular-Design-Studio/0.1"},
                )
                try:
                    with urlopen(request, timeout=30) as response:
                        content = response.read(50_000_001)
                    converted_from_mmcif = False
                except HTTPError as exc:
                    if exc.code != 404:
                        raise
                    cif_request = Request(
                        f"https://files.rcsb.org/download/{pdb_code}.cif",
                        headers={"User-Agent": "MN-Molecular-Design-Studio/0.1"},
                    )
                    with urlopen(cif_request, timeout=30) as response:
                        cif_content = response.read(50_000_001)
                    if len(cif_content) > 50_000_000:
                        raise ValueError("The mmCIF file is larger than the 50 MB import limit.")
                    content = _rcsb_mmcif_to_pdb(cif_content, allow_non_protein=allow_non_protein)
                    converted_from_mmcif = True
                if len(content) > 50_000_000:
                    raise ValueError("The PDB file is larger than the 50 MB import limit.")
                stored_path = _target_storage_path(content, f"{pdb_code}.pdb", owner_id=owner_id, allow_non_protein=allow_non_protein)
                st.session_state.pop(upload_key, None)
                st.session_state[stored_path_key] = str(stored_path)
                st.session_state[stored_name_key] = pdb_code
                st.session_state[converted_mmcif_key] = converted_from_mmcif
            except (HTTPError, URLError, OSError, ValueError) as exc:
                st.error(f"Could not fetch {pdb_code}: {exc}")

    fetched_path = st.session_state.get(stored_path_key)
    if fetched_path:
        try:
            Path(fetched_path).expanduser().resolve().relative_to(user_input_root(owner_id))
        except (OSError, ValueError):
            st.session_state.pop(stored_path_key, None)
            st.session_state.pop(stored_name_key, None)
            fetched_path = None
    has_fetched_target = bool(fetched_path and Path(fetched_path).is_file())
    uploader = st.file_uploader
    structure_file_types = upload_file_types or (["pdb", "ent", "cif", "mmcif"] if allow_cif else ["pdb", "ent"])
    upload_expander_label = upload_expander_label or (
        "Upload a PDB or mmCIF file instead" if allow_cif else "Upload a PDB file instead"
    )
    if upload_in_expander:
        with st.expander(upload_expander_label, expanded=not has_fetched_target):
            uploaded = uploader(upload_label, type=structure_file_types, key=upload_key, help=upload_help)
    else:
        uploaded = uploader(upload_label, type=structure_file_types, key=upload_key, help=upload_help)

    if uploaded is not None:
        try:
            return TargetSource(
                _target_storage_path(uploaded.getvalue(), uploaded.name, owner_id=owner_id, allow_non_protein=allow_non_protein),
                Path(uploaded.name).stem,
                owner_id,
            )
        except (OSError, ValueError) as exc:
            st.error(str(exc))
            return TargetSource(owner_id=owner_id)
    if fetched_path and Path(fetched_path).is_file():
        path = Path(fetched_path)
        if st.session_state.get(converted_mmcif_key):
            st.info(
                "RCSB provides this entry as mmCIF only. It was converted to PDB, "
                "with long chain IDs mapped to PDB-compatible IDs."
            )
        return TargetSource(path, st.session_state.get(stored_name_key, path.stem), owner_id)
    return TargetSource(owner_id=owner_id)


def _target_widget_keys(key_prefix: str, path: Path) -> dict[str, str]:
    digest = _hash(str(path))
    return {
        "hotspots": f"{key_prefix}-hotspots-{digest}",
        "coldspots": f"{key_prefix}-coldspots-{digest}",
        "viewer_mode": f"{key_prefix}-viewer-mode-{digest}",
    }


def _rerun(scope: str | None = None) -> None:
    if scope:
        st.rerun(scope=scope)
    else:
        st.rerun()


def _render_target_trimming(
    target: TargetDesignInput,
    protein_residues_by_chain: dict[str, list[tuple[int, str]]],
    *,
    rerun_scope: str | None = None,
) -> dict[str, tuple[int, int]]:
    retained_ranges: dict[str, tuple[int, int]] = {}
    with st.expander(f"Trim {target.label.lower()}-chain termini", expanded=False):
        description_col, reset_col = st.columns([4, 1], vertical_alignment="center")
        with description_col:
            st.caption(
                "Choose the PDB residue range to retain for each selected chain. Other chains stay in the structure, "
                "and the original residue numbers are preserved."
            )
        with reset_col:
            reset_trim = st.button(
                "Reset to full chain",
                key=f"{target.key_prefix}-reset-trim-" + _hash(str(target.source_path)),
                use_container_width=True,
                disabled=not bool(target.chains),
            )
        if reset_trim:
            for chain in target.chains:
                residues = protein_residues_by_chain.get(chain, [])
                if not residues:
                    continue
                trim_key = f"{target.key_prefix}-trim-" + _hash(f"{target.source_path}:{chain}")
                st.session_state[trim_key + "-start"] = residues[0][0]
                st.session_state[trim_key + "-end"] = residues[-1][0]
            _rerun(rerun_scope)
        if not target.chains:
            st.info(f"Select a {target.label.lower()} chain to set terminal trimming.")
        for chain in target.chains:
            residues = protein_residues_by_chain.get(chain, [])
            if not residues:
                continue
            trim_key = f"{target.key_prefix}-trim-" + _hash(f"{target.source_path}:{chain}")
            start_key, end_key = trim_key + "-start", trim_key + "-end"
            first_residue, last_residue = residues[0][0], residues[-1][0]
            start_minimum = first_residue - 1
            start_initial = min(max(start_minimum, int(st.session_state.get(start_key, first_residue))), last_residue)
            st.session_state[start_key] = start_initial
            st.markdown(f"**Chain {chain}** · {first_residue}–{last_residue} · {len(residues):,} residues")
            start_col, end_col = st.columns(2)
            with start_col:
                start_residue = st.number_input(
                    "Start residue to retain", min_value=start_minimum, max_value=last_residue,
                    value=start_initial, step=1, key=start_key,
                )
            end_initial = min(max(int(start_residue), int(st.session_state.get(end_key, last_residue))), last_residue)
            st.session_state[end_key] = end_initial
            with end_col:
                end_residue = st.number_input(
                    "End residue to retain", min_value=int(start_residue), max_value=last_residue,
                    value=end_initial, step=1, key=end_key,
                )
            retained_ranges[chain] = (int(start_residue), int(end_residue))
            retained_ids = [residue_id for residue_id in residues if int(start_residue) <= residue_id[0] <= int(end_residue)]
            if not retained_ids:
                st.error(f"Chain {chain}: the selected range contains no protein residues.")
                target.valid_trim_ranges = False
                continue
            st.caption(
                f"Retained: {len(retained_ids):,} residues · PDB numbering {retained_ids[0][0]}–{retained_ids[-1][0]}"
            )
    return retained_ranges


def render_target_editor(
    source: TargetSource,
    *,
    key_prefix: str,
    label: str = "Target",
    chain_selection_label: str | None = None,
    chain_selection_disabled: bool = False,
    prefer_chain_b: bool = False,
    legacy_hotspots_key: str | None = None,
    residue_selection_title: str = "Hotspot residues",
    residue_selection_caption: str | None = None,
    residue_selection_input_label: str = "Chosen residues",
    show_structure_tools: bool = True,
    allow_trimming: bool = True,
    viewer_height: int = 640,
    rerun_scope: str | None = None,
    before_residue_selection: Callable[[TargetDesignInput], None] | None = None,
) -> TargetDesignInput:
    """Render chain selection, terminal trimming, hotspots, and the target MolView."""
    target = TargetDesignInput(label=label, key_prefix=key_prefix, source_path=source.path)
    if source.path is None or not source.path.is_file():
        return target

    is_mmcif = source.path.suffix.lower() in {".cif", ".mmcif"}
    target.source_text = source.path.read_text(errors="ignore")
    target.design_text = target.source_text
    target.sequence_by_chain = (
        _mmcif_sequence_entries_by_chain(source.path)
        if is_mmcif
        else _sequence_entries_by_chain(target.source_text)
    )
    chain_options = sorted(target.sequence_by_chain)
    if not chain_options:
        st.error(f"{label} structure has no readable protein chains.")
        target.source_path = None
        return target

    st.success(f"{source.name or label} loaded · {source.path.name}")
    chain_suffix = "-prefer-target-chain-b" if prefer_chain_b else ""
    chain_key = f"{key_prefix}-chains-" + _hash(str(source.path) + chain_suffix)
    if chain_key not in st.session_state:
        preferred = "B" if prefer_chain_b and "B" in chain_options else chain_options[0]
        st.session_state[chain_key] = [preferred]
    target.chains = st.multiselect(chain_selection_label or f"{label} chains", chain_options, key=chain_key, disabled=chain_selection_disabled)

    if is_mmcif:
        st.info("mmCIF structure loaded. Chain and residue selections use the author IDs in the file; terminal trimming is available for PDB files.")
        retained_ranges: dict[str, tuple[int, int]] = {}
    else:
        protein_residues_by_chain = _protein_residue_ids_by_chain(target.source_text)
        retained_ranges = (_render_target_trimming(target, protein_residues_by_chain, rerun_scope=rerun_scope) if allow_trimming else
                           {chain:(rows[0][0],rows[-1][0]) for chain,rows in protein_residues_by_chain.items() if chain in target.chains and rows})

    is_trimmed = not is_mmcif and target.valid_trim_ranges and any(
        start > protein_residues_by_chain[chain][0][0]
        or end < protein_residues_by_chain[chain][-1][0]
        for chain, (start, end) in retained_ranges.items()
    )
    if is_trimmed:
        target.design_text = _trim_target_chain_ranges(target.source_text, retained_ranges)
        target.design_path = _derived_target_storage_path(
            target.design_text,
            f"{source.name or label}-trimmed.pdb",
            owner_id=source.owner_id,
        )
    else:
        target.design_path = source.path
    target.sequence_by_chain = (
        _mmcif_sequence_entries_by_chain(source.path)
        if is_mmcif
        else _sequence_entries_by_chain(target.design_text)
    )
    target.valid_residues = {
        (chain, number)
        for chain in target.chains
        for number, _ in target.sequence_by_chain.get(chain, [])
    }
    residue_count = sum(len(target.sequence_by_chain.get(chain, [])) for chain in target.chains)
    st.caption(f"{len(target.chains)} chain(s) selected · {residue_count:,} residues" + (" after trimming" if allow_trimming else ""))

    keys = _target_widget_keys(key_prefix, source.path)
    hotspot_key = keys["hotspots"]
    target.hotspot_state_key = hotspot_key
    target.coldspot_state_key = keys["coldspots"]
    if hotspot_key not in st.session_state:
        if legacy_hotspots_key and legacy_hotspots_key in st.session_state:
            st.session_state[hotspot_key] = st.session_state.pop(legacy_hotspots_key)
        else:
            st.session_state[hotspot_key] = ""
    pending_key = hotspot_key + "-pending-viewer-selection"
    if pending_key in st.session_state:
        st.session_state[hotspot_key] = st.session_state.pop(pending_key)
    revision_key = hotspot_key + "-revision"
    st.session_state.setdefault(revision_key, 0)
    normalized_hotspots = _hotspot_text(_hotspot_pairs(st.session_state[hotspot_key]) & target.valid_residues)
    if normalized_hotspots != st.session_state[hotspot_key]:
        st.session_state[hotspot_key] = normalized_hotspots

    if before_residue_selection is not None:
        before_residue_selection(target)

    st.markdown(f"#### {residue_selection_title}")
    is_mmcif = source.path.suffix.lower() in {".cif", ".mmcif"}
    if residue_selection_caption:
        st.caption(residue_selection_caption)
    elif is_mmcif:
        st.caption("Enter residues by chain ID and residue number, or select residues in the mmCIF structure viewer.")
    else:
        st.caption("Choose residues in MolView’s sequence panel or click them in the structure. Both views stay in sync.")
    input_col, select_all_col, clear_col = st.columns([5, .85, .85], vertical_alignment="bottom")
    with select_all_col:
        if st.button("Select all", key=f"{hotspot_key}-select-all", use_container_width=True, disabled=not target.valid_residues):
            st.session_state[hotspot_key] = _hotspot_text(target.valid_residues)
            _increment_session_counter(revision_key)
            _rerun(rerun_scope)
    with clear_col:
        if st.button("Clear", key=f"{hotspot_key}-clear", use_container_width=True, disabled=not _hotspot_pairs(st.session_state[hotspot_key])):
            st.session_state[hotspot_key] = ""
            _increment_session_counter(revision_key)
            _rerun(rerun_scope)
    with input_col:
        target.hotspots = st.text_input(
            residue_selection_input_label,
            placeholder=f"e.g. {target.chains[0] if target.chains else 'B'}340, {target.chains[0] if target.chains else 'B'}342",
            key=hotspot_key,
            on_change=_increment_session_counter,
            args=(revision_key,),
            help="Enter chain and residue numbers separated by commas, or select them in the structure viewer.",
        ).strip()
    valid_ranges = [
        f"{chain}{residues[0][0]}–{chain}{residues[-1][0]}"
        for chain in target.chains
        if (residues := target.sequence_by_chain.get(chain))
    ]
    if valid_ranges:
        st.caption("Valid target ranges: " + ", ".join(valid_ranges))

    if is_mmcif:
        target.viewer_mode = "Whole complex"
        st.markdown("#### MolView")
        st.caption("The mmCIF viewer displays the full structure; residue selections are restricted to the target chains above.")
    else:
        st.markdown("#### MolView")
        target.viewer_mode = st.radio(
            "Show in viewer",
            ["Target only", "Whole complex"],
            horizontal=True,
            key=keys["viewer_mode"],
            help="Target only displays the selected target chains. Whole complex displays every chain in the structure.",
        )
    _render_residue_viewer(
        target,
        kind="hotspots",
        height=viewer_height,
        show_structure_tools=show_structure_tools,
        rerun_scope=rerun_scope,
    )
    target.hotspots = st.session_state[hotspot_key]
    target.coldspots = st.session_state.get(keys["coldspots"], "")
    return target


def _render_residue_viewer(
    target: TargetDesignInput,
    *,
    kind: str,
    height: int,
    show_structure_tools: bool = True,
    rerun_scope: str | None = None,
) -> None:
    if target.source_path is None or not target.source_path.is_file():
        st.info("Load a target structure to select residues in MolView.")
        return
    if target.viewer_mode == "Target only" and not target.chains:
        st.info("Select a target chain to show the structure.")
        return
    if target.source_path.suffix.lower() in {".cif", ".mmcif"}:
        viewer_text = target.design_text
    else:
        viewer_text = filter_pdb_text(
            target.design_text,
            keep_chains=set(target.chains) if target.viewer_mode == "Target only" else None,
        )
    digest = _hash(str(target.source_path) + ",".join(target.chains) + hashlib.sha1(viewer_text.encode()).hexdigest(), 12)
    if kind == "hotspots":
        viewer_key = f"{target.key_prefix}-viewer-{digest}-selection-sync-v2"
        selection_key = target.hotspot_state_key
    else:
        viewer_key = f"{target.key_prefix}-coldspot-viewer-{digest}-selection-sync-v1"
        selection_key = target.coldspot_state_key
    revision_key = selection_key + "-revision"
    value = molstar_custom_component(
        [
            StructureVisualization(
                pdb=viewer_text,
                color="uniform",
                color_params={"value": "0xe7b8aa"},
                representation_type="cartoon",
                highlighted_selections=[
                    f"{chain}{residue}" for chain, residue in sorted(_hotspot_pairs(st.session_state.get(selection_key, "")))
                ],
            )
        ],
        key=viewer_key,
        height=height,
        show_controls=True,
        show_structure_tools=show_structure_tools,
        selection_mode=True,
        force_reload=st.session_state.get(revision_key, 0),
    )
    signature = json.dumps(value, sort_keys=True, default=str)
    seen_key = viewer_key + "-last-selection"
    previous_signature = st.session_state.get(seen_key)
    st.session_state[seen_key] = signature
    if previous_signature is not None and signature != previous_signature:
        picked = _selected_viewer_residues(value, _residue_map(viewer_text, target.source_path))
        picked_targets = {residue for residue in picked if residue[0] in target.chains}
        if picked and not picked_targets:
            return
        next_value = _hotspot_text(picked_targets)
        if next_value != st.session_state.get(selection_key, ""):
            st.session_state[selection_key + "-pending-viewer-selection"] = next_value
            _rerun(rerun_scope)


def render_target_coldspots(
    target: TargetDesignInput,
    *,
    key_prefix: str,
    legacy_coldspots_key: str | None = None,
    show_structure_tools: bool = True,
    rerun_scope: str | None = None,
) -> str:
    """Render a target's coldspot input and independent MolView selection."""
    if target.source_path is None or not target.source_path.is_file():
        target.coldspots = ""
        return target.coldspots
    keys = _target_widget_keys(key_prefix, target.source_path)
    coldspot_key = keys["coldspots"]
    target.coldspot_state_key = coldspot_key
    if coldspot_key not in st.session_state:
        if legacy_coldspots_key and legacy_coldspots_key in st.session_state:
            st.session_state[coldspot_key] = st.session_state.pop(legacy_coldspots_key)
        else:
            st.session_state[coldspot_key] = ""
    pending_key = coldspot_key + "-pending-viewer-selection"
    if pending_key in st.session_state:
        st.session_state[coldspot_key] = st.session_state.pop(pending_key)
    revision_key = coldspot_key + "-revision"
    st.session_state.setdefault(revision_key, 0)
    normalized = _hotspot_text(_hotspot_pairs(st.session_state[coldspot_key]) & target.valid_residues)
    if normalized != st.session_state[coldspot_key]:
        st.session_state[coldspot_key] = normalized

    with st.container(border=True):
        st.markdown("#### Coldspot residues")
        st.caption("Choose residues the binder should avoid. Selection here is separate from hotspots.")
        input_col, select_all_col, clear_col = st.columns([5, .85, .85], vertical_alignment="bottom")
        with select_all_col:
            if st.button("Select all", key=f"{coldspot_key}-select-all", use_container_width=True, disabled=not target.valid_residues):
                st.session_state[coldspot_key] = _hotspot_text(target.valid_residues)
                _increment_session_counter(revision_key)
                _rerun(rerun_scope)
        with clear_col:
            if st.button("Clear", key=f"{coldspot_key}-clear", use_container_width=True, disabled=not _hotspot_pairs(st.session_state[coldspot_key])):
                st.session_state[coldspot_key] = ""
                _increment_session_counter(revision_key)
                _rerun(rerun_scope)
        with input_col:
            target.coldspots = st.text_input(
                "Chosen coldspot residues",
                placeholder=f"e.g. {target.chains[0] if target.chains else 'B'}410, {target.chains[0] if target.chains else 'B'}420",
                key=coldspot_key,
                on_change=_increment_session_counter,
                args=(revision_key,),
                help=(
                    "Enter chain and residue numbers separated by commas."
                    if target.source_path.suffix.lower() in {".cif", ".mmcif"}
                    else "Enter chain and residue numbers separated by commas, or select them in the coldspot MolView."
                ),
            ).strip()
        if target.source_path.suffix.lower() in {".cif", ".mmcif"}:
            st.caption("This viewer displays the full mmCIF structure; residue selections are restricted to the target chains above.")
        else:
            st.caption("This viewer uses the same target-only / whole-complex view setting as the hotspot viewer.")
        _render_residue_viewer(
            target,
            kind="coldspots",
            height=520,
            show_structure_tools=show_structure_tools,
            rerun_scope=rerun_scope,
        )
    target.coldspots = st.session_state[coldspot_key]
    return target.coldspots
