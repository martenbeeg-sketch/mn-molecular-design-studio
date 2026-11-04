"""Shared molecular complex editor and its workflow payload conversion."""
from __future__ import annotations

from pathlib import Path
from typing import Any
import base64
import hashlib
from functools import lru_cache
import importlib.util
import os
import subprocess
import sys

import streamlit as st
import streamlit.components.v1 as components
import yaml

_asset_path = Path(__file__).with_name("components") / "cofolding_editor"
_editors: dict[str, Any] = {}
CHAIN_IDS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"


def _get_editor() -> Any:
    asset_hash = hashlib.sha256(b"".join((_asset_path / name).read_bytes() for name in (
        "index.html", "editor.js", "editor.css", "examples.json",
    ))).hexdigest()[:12]
    if asset_hash not in _editors:
        _editors[asset_hash] = components.declare_component(
            "molecular_complex_editor_" + asset_hash, path=str(_asset_path),
        )
    return _editors[asset_hash]


@lru_cache(maxsize=1)
def _localfold_editor_profiles() -> dict[str, Any]:
    """Ship model controls together so selecting a model needs no server round trip."""
    from mn_molecular_design_studio import localfold
    from mn_molecular_design_studio.engine.workflow_engine.workflows.optimized_prediction import ENGINE_PROFILES
    profiles = {}
    for name in localfold.UI_MODEL_KEYS:
        profile = localfold.profile_for(ENGINE_PROFILES['localfold'], name, ui=True)
        profiles[name] = {
            'prediction': prediction_settings_for(profile),
            'types': ['protein', *(['dna', 'rna', 'ligand'] if profile['molecule_inputs'] else [])],
            'ccd': bool(profile['molecule_inputs']), 'smiles': bool(profile['molecule_inputs']),
            'cyclic': profile['supports_cyclic'], 'modifications': profile['supports_modifications'],
            'bonds': profile['supports_bonds'], 'crosslinks': profile['supports_crosslinks'],
            'contacts': profile['supports_contacts'], 'pockets': profile['supports_pockets'],
            'msa': profile['supports_msa_upload'], 'templates': profile['supports_templates'],
        }
    return profiles


@lru_cache(maxsize=128)
def _smiles_preview(smiles: str) -> str:
    if not smiles or len(smiles) > 5000:
        return ""
    configured = os.getenv("MN_MOLECULAR_DESIGN_STUDIO_CHEM_PYTHON")
    python = Path(configured) if configured else (
        Path(sys.executable) if importlib.util.find_spec("rdkit") else None
    )
    if python is None:
        return ""
    if not python.is_file():
        return ""
    script = (
        "import sys\nfrom rdkit import Chem\nfrom rdkit.Chem.Draw import rdMolDraw2D\n"
        "m=Chem.MolFromSmiles(sys.stdin.read())\n"
        "if m is not None:\n d=rdMolDraw2D.MolDraw2DSVG(240,180);"
        "rdMolDraw2D.PrepareAndDrawMolecule(d,m);d.FinishDrawing();print(d.GetDrawingText())\n"
    )
    try:
        result = subprocess.run([str(python), "-c", script], input=smiles, text=True,
                                capture_output=True, timeout=6, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return ""
    if result.returncode or "<svg" not in result.stdout:
        return ""
    return "data:image/svg+xml;base64," + base64.b64encode(result.stdout.encode()).decode()


def default_model() -> dict[str, Any]:
    return {
        "cards": [{"id": 1, "kind": "protein", "value": "", "copies": 1, "modifications": []}],
        "compound_id": "", "compound_name": "",
        "restraints": [], "templates": [], "msa_files": [], "yaml_files": [], "yaml_text": "",
        "settings": {
            "samples": 5, "optimization_mode": "off", "use_yaml": False,
            "num_recycles": 3, "step_scale": 1.5, "seed_count": 1, "seed": 0,
            "output_format": "pdb", "boltz_version": "2.2.1", "method": "",
            "template_distance_threshold": "", "use_potentials": False, "use_msa": True,
            "max_msa_seqs": "", "generate_ipsae": True, "sampling_steps": 200,
            "max_parallel_samples": 5, "predict_affinity": False, "affinity_binder": "",
            "prioritize": False,
        },
    }


def model_to_payload(model: dict[str, Any]) -> dict[str, Any]:
    """Keep the visual editor independent of each engine's serialized schema."""
    cards = model.get("cards") or []
    settings = {**default_model()["settings"], **(model.get("settings") or {})}
    localfold = "msa_mode" in (model.get("settings") or {})
    global_use_msa = settings.get("msa_mode") != "single_sequence" if localfold else bool(settings["use_msa"])
    groups: dict[str, list[dict[str, Any]]] = {kind: [] for kind in ("protein", "dna", "rna", "ccd", "smiles")}
    errors = []
    for card_index, card in enumerate(cards):
        if localfold and card.get('kind') in {'modification', 'bond', 'crosslink', 'contact', 'pocket'}:
            continue
        kind = str(card.get("mode", "ccd")) if card.get("kind") == "ligand" else str(card.get("kind", "protein"))
        value = str(card.get("value") or "").strip()
        if kind not in groups:
            errors.append("Choose a supported molecule type.")
            continue
        if kind in {"protein", "dna", "rna", "ccd"}:
            value = "".join(value.split()).upper()
            alphabet = {"protein": set("ACDEFGHIKLMNPQRSTVWYBXZJUO"), "dna": set("ACGTN"), "rna": set("ACGUN")}.get(kind)
            if alphabet is not None and set(value) - alphabet:
                errors.append(f"The {kind.upper()} sequence contains unsupported characters.")
        elif "\n" in value or "\r" in value:
            errors.append("Enter one ligand SMILES string per chain card.")
        if not value:
            errors.append("All sequences must be non empty.")
        copies = int(card.get("copies", 1))
        if not 1 <= copies <= len(CHAIN_IDS):
            errors.append("Copies must be between 1 and 52.")
            continue
        groups[kind].extend([{**card, "value": value, '_editor_card_index': card_index}] * copies)
    if sum(len(rows) for rows in groups.values()) > len(CHAIN_IDS):
        errors.append("Enter at most 52 chains and ligand copies.")
    options = []
    catalog = []
    all_rows = [row for rows in groups.values() for row in rows][:len(CHAIN_IDS)]
    extra_modifications: dict[int, list[dict]] = {}
    for card in cards if localfold else []:
        if card.get('kind') != 'modification':
            continue
        selected_chain = str(card.get('chain_id') or '')
        target = next((row for index, row in enumerate(all_rows) if CHAIN_IDS[index] == selected_chain), None)
        if target is None or target.get('kind') not in {'protein', 'dna', 'rna'}:
            errors.append('Choose a protein, DNA or RNA chain for each modification.')
            continue
        extra_modifications.setdefault(target['_editor_card_index'], []).append({
            'position': card.get('position', ''), 'code': str(card.get('code') or '').strip(),
        })
    for kind, rows in groups.items():
        for row in rows:
            if len(catalog) >= len(CHAIN_IDS):
                break
            chain_id = CHAIN_IDS[len(catalog)]
            catalog.append({"id": chain_id, "kind": kind, "sequence": row["value"]})
            if kind == "protein" or localfold and kind in {'dna', 'rna'}:
                modifications = [*(row.get("modifications") or []), *extra_modifications.get(row['_editor_card_index'], [])]
                if any(not str(mod.get("code") or "").strip() for mod in modifications):
                    errors.append("Every modification needs a CCD code.")
                modification_lines = []
                for mod in modifications:
                    try:
                        position = int(mod.get('position', ''))
                    except (ValueError, TypeError):
                        errors.append('Every modification needs a residue position.')
                        continue
                    modification_lines.append(f"{position}, {str(mod.get('code') or '').strip()}")
                options.append({
                    "id": chain_id, "cyclic": bool(row.get("cyclic")),
                    "modifications_text": "\n".join(modification_lines),
                    "msa_a3m_text": next((str(item.get("text") or "") for item in model.get("msa_files", [])
                                           if item.get("chain_id") == chain_id), "")
                    if global_use_msa and (not localfold or row.get("localfold_alignment") != "single_sequence") else "",
                    **({
                        "use_msa": global_use_msa and row.get("localfold_alignment") != "single_sequence",
                        "use_template": str(row.get("localfold_use_template", True)).lower() != "false",
                    } if localfold else {}),
                })
    if global_use_msa and any(not item.get("chain_id") for item in model.get("msa_files", [])):
        errors.append("Assign each custom A3M file to a protein chain.")
    values: dict[str, Any] = {
        "protein_sequences": ":".join(row["value"] for row in groups["protein"]),
        "dna_sequences": ":".join(row["value"] for row in groups["dna"]),
        "rna_sequences": ":".join(row["value"] for row in groups["rna"]),
        "ligand_ccd_codes": ":".join(row["value"] for row in groups["ccd"]),
        "ligand_smiles": "\n".join(row["value"] for row in groups["smiles"]),
        "molecule_chain_options": options, "restraints": [*(model.get("restraints") or []), *([
            {card['kind']: {key: (str(card.get(key) or '').replace(';', '\n') if key == 'contacts_input' else card.get(key)) for key in {
                'bond': ('atom1', 'atom2'), 'crosslink': ('name', 'residue1', 'residue2'),
                'contact': ('token1', 'token2', 'max_distance'),
                'pocket': ('binder', 'contacts_input', 'max_distance'),
            }[card['kind']]}}
            for card in cards if card.get('kind') in {'bond', 'crosslink', 'contact', 'pocket'}
        ] if localfold else [])],
        "boltz_templates": model.get("templates") or [],
        "template_distance_threshold": settings.get("template_distance_threshold") or None,
        "affinity_binder": settings.get("affinity_binder") if settings.get("predict_affinity") else "",
        "input_error": errors[0] if errors else "", "chain_catalog": catalog,
    }
    has_ligand = bool(groups["ccd"] or groups["smiles"])
    if not has_ligand and not groups["dna"] and not groups["rna"]:
        from mn_molecular_design_studio.target_identity import protein_sequence_reference

        reference = protein_sequence_reference([row["value"] for row in groups["protein"]])
        if reference:
            values["sequence_source_reference"] = reference
    compound_name = str(model.get("compound_name") or "").strip()
    compound_id = str(model.get("compound_id") or "").strip() or compound_name
    if has_ligand and not compound_name:
        values["input_error"] = values["input_error"] or "Enter a compound name for this ligand complex."
    values.update({
        "compound_id": compound_id if has_ligand else "",
        "compound_name": compound_name if has_ligand else "",
        "name": compound_name if has_ligand else "",
    })
    if settings["use_yaml"]:
        values.update({key: "" for key in ("protein_sequences", "dna_sequences", "rna_sequences", "ligand_ccd_codes", "ligand_smiles")})
        values.update({
            "molecule_chain_options": [], "restraints": [],
            "affinity_binder": "",
            "boltz_yaml_text": str(model.get("yaml_text") or ""),
            "boltz_yaml_files": model.get("yaml_files") or [],
            "boltz_custom_msas": model.get("msa_files") or [],
            "input_error": "" if str(model.get("yaml_text") or "").strip() else "Upload a Boltz YAML file.",
        })
    return values


def render_molecule_editor(
    key_prefix: str, supported_inputs: tuple[str, ...], *, supports_cyclic: bool = False,
    supports_modifications: bool = False, supports_restraints: bool = False,
    supports_msa_upload: bool = False, boltz_settings: bool = False,
    supports_examples: bool = True, supports_affinity: bool = False,
    supports_affinity_control: bool = False,
    supports_templates: bool = False,
    localfold: bool = False,
    modes: tuple[str, ...] = ("off", "fast", "big"), can_submit: bool = False,
    prediction_settings: dict[str, Any] | None = None, settings_only: bool = False,
) -> tuple[dict[str, Any], dict[str, Any], bool]:
    state_key = f"{key_prefix}-editor-model"
    database_user = st.session_state.get("_molecular_database_user")
    if database_user and not settings_only:
        from mn_molecular_design_studio.database_ui import render_molecule_picker
        from mn_molecular_design_studio.database_selections import set_editor_selection
        allowed = {"protein", *[kind for kind, field in (("dna", "dna_sequences"), ("rna", "rna_sequences")) if field in supported_inputs]}
        if {"ligand_ccd_codes", "ligand_smiles"}.intersection(supported_inputs):
            allowed.add("small_molecule")
        selected = render_molecule_picker(database_user, key=key_prefix, allowed_types=allowed, lazy=localfold)
        if selected:
            set_editor_selection(key_prefix, selected)
    if localfold:
        # Remember which LocalFold form owns the selected model so the next
        # rerun can build that model's matching controls before rendering it.
        st.session_state["_localfold_active_editor_key"] = state_key
    initial = st.session_state.get(f"{key_prefix}-editor") or st.session_state.get(state_key)
    if initial is None and localfold and key_prefix.startswith("optimized-localfold-"):
        selected_model = str(st.session_state.get("localfold-model", "alphafold3"))
        suffix = key_prefix.removeprefix("optimized-localfold-")
        initial = st.session_state.get(f"optimized-localfold-{selected_model}-{suffix}-editor")
    if initial is None:
        initial = default_model()
        if prediction_settings:
            initial["settings"].update({field["key"]: field["default"] for field in prediction_settings["fields"]})
            initial["settings"]["samples"] = prediction_settings.get("samples", 5)
    if localfold:
        initial.setdefault('settings', {}).setdefault('msa_mode', 'automatic')
        initial['settings'].setdefault('model', st.session_state.get('localfold-model', 'alphafold3'))
    if prediction_settings:
        for field in prediction_settings["fields"]:
            initial.setdefault("settings", {}).setdefault(field["key"], field["default"])
            if initial["settings"][field["key"]] == "" and field.get("type", "int") in {"int", "float"}:
                initial["settings"][field["key"]] = field["default"]
    yaml_chains = []
    yaml_error = ""
    if initial.get("settings", {}).get("use_yaml") and initial.get("yaml_text"):
        try:
            payload = yaml.safe_load(initial["yaml_text"])
            if not isinstance(payload, dict) or not isinstance(payload.get("sequences"), list) or not payload["sequences"]:
                raise ValueError("Boltz YAML needs a non-empty sequences list.")
            for row in payload["sequences"]:
                for kind, molecule in row.items():
                    if not isinstance(molecule, dict):
                        continue
                    chain_ids = molecule.get("id")
                    for chain_id in chain_ids if isinstance(chain_ids, list) else [chain_ids]:
                        if chain_id:
                            yaml_chains.append({"id": str(chain_id), "kind": kind, "value": molecule.get("sequence", "")})
        except (yaml.YAMLError, ValueError, AttributeError) as exc:
            yaml_error = str(exc)
    previews = {
        str(card["id"]): _smiles_preview(str(card["value"])) for card in initial.get("cards", [])
        if card.get("kind") == "ligand" and card.get("mode") == "smiles" and card.get("value")
    }
    localfold_models = {}
    if localfold:
        from mn_molecular_design_studio import localfold as localfold_engine
        for model_key in localfold_engine.UI_MODEL_KEYS:
            spec = localfold_engine.MODELS[model_key]
            no_msa = model_key in {'esmfold2_lm600m', 'esmfold2_lm300m'}
            localfold_models[model_key] = {
                **spec,
                "no_msa": no_msa,
                "samples": 1 if model_key.startswith("af2_") or model_key == "bioemu2" else 5,
            }
    capabilities = {
        "types": ["protein", *[kind for kind, field in (("rna", "rna_sequences"), ("dna", "dna_sequences")) if field in supported_inputs],
                  *(["ligand"] if {"ligand_ccd_codes", "ligand_smiles"}.intersection(supported_inputs) else [])],
        "ccd": "ligand_ccd_codes" in supported_inputs, "smiles": "ligand_smiles" in supported_inputs,
        "examples": supports_examples and not settings_only,
        "ligand_examples": supports_examples and not settings_only and "ligand_smiles" in supported_inputs,
        "affinity": supports_affinity and not settings_only and "ligand_smiles" in supported_inputs,
        "affinity_control": supports_affinity_control and not settings_only and "ligand_smiles" in supported_inputs,
        "cyclic": supports_cyclic, "modifications": supports_modifications,
        "restraints": supports_restraints, "msa": supports_msa_upload,
        "templates": supports_templates,
        "localfold": localfold, "localfold_models": localfold_models,
        "localfold_model": st.session_state.get('localfold-model', 'alphafold3') if localfold else '',
        "localfold_profiles": _localfold_editor_profiles() if localfold else {},
        "localfold_crosslinkers": localfold_engine.CROSSLINKERS if localfold else (),
        "localfold_model_groups": localfold_engine.MODEL_UI_GROUPS if localfold else (),
        "localfold_esmfold2_variants": localfold_engine.ESMFOLD2_VARIANTS if localfold else (),
        "localfold_default_esmfold2_variant": localfold_engine.DEFAULT_ESMFOLD2_VARIANT if localfold else "",
        "boltz": boltz_settings, "modes": modes, "can_submit": can_submit,
        "prediction": prediction_settings, "settings_only": settings_only,
    }
    model = _get_editor()(
        initial=initial, previews=previews, yaml_chains=yaml_chains, yaml_error=yaml_error,
        capabilities=capabilities,
        key=(f"{key_prefix}-editor-library-{st.session_state[key_prefix + '-library-epoch']}"
             if st.session_state.get(key_prefix + "-library-epoch") else f"{key_prefix}-editor"), default=initial,
    )
    st.session_state[state_key] = model
    if prediction_settings:
        for field in prediction_settings["fields"]:
            model.setdefault("settings", {}).setdefault(field["key"], field["default"])
            if model["settings"][field["key"]] == "" and field.get("type", "int") in {"int", "float"}:
                model["settings"][field["key"]] = field["default"]
    token = model.get("submit_token")
    last_token_key = f"{key_prefix}-last-submit-token"
    submitted = bool(token and token != st.session_state.get(last_token_key))
    if submitted:
        st.session_state[last_token_key] = token
    values = model_to_payload(model)
    if not boltz_settings:
        values["structure_templates"] = values.pop("boltz_templates", []) if supports_templates else []
    if model.get("settings", {}).get("use_yaml"):
        yaml_has_ligand = any(row.get("kind") == "ligand" for row in yaml_chains)
        if yaml_has_ligand:
            values["compound_name"] = str(model.get("compound_name") or "").strip()
            values["compound_id"] = str(model.get("compound_id") or "").strip() or values["compound_name"]
            values["name"] = values["compound_name"]
            if not values["compound_name"]:
                values["input_error"] = "Enter a compound name for this ligand complex."
        else:
            values.update({"compound_id": "", "compound_name": "", "name": ""})
    if yaml_error and model.get("settings", {}).get("use_yaml"):
        values["input_error"] = yaml_error
    if (model.get("settings", {}).get("use_yaml") and any(row.get("kind") == "ligand" for row in yaml_chains)
            and not str(model.get("compound_name") or "").strip()):
        values["input_error"] = "Enter a compound name for this ligand complex."
    if not (supports_cyclic or supports_modifications or supports_msa_upload):
        values["molecule_chain_options"] = []
    return values, model["settings"], submitted


def prediction_settings_for(profile: dict[str, Any], *, extra_fields: tuple[dict, ...] = ()) -> dict[str, Any]:
    """Describe an engine's controls for the shared Boltz-style editor."""
    fields = []
    if len(profile["modes"]) > 1:
        fields.append(
            {"key": "optimization_mode", "label": "Claude-inspired inference optimization", "type": "select",
             "default": "off", "choices": [(mode, "None" if mode == "off" else mode.title()) for mode in profile["modes"]],
             "help": "Choose an inference profile available for this engine."}
        )
    recycle_field = {
        "key": "num_recycles", "label": "Number of recycles", "default": profile.get("default_recycles", 10),
        "min": 1, "max": 100, "hidden": bool(profile.get("settings_replace_num_recycles")),
        "help": "Refinement rounds, mapped to this engine's recycle, cycle or loop setting.",
    }
    if profile.get("api_engine") == "localfold":
        recycle_field.update(type="select", min=0, choices=tuple((value, str(value)) for value in (0, 1, 2, 3, 6, 12, 24)))
    fields.extend([
        recycle_field,
        *extra_fields,
        *profile.get("settings", ()),
        {"key": "cpu_cores", "label": "CPU cores", "default": 4, "min": 1, "max": 256},
        {"key": "prioritize", "label": "Prioritize in scheduler queue", "type": "bool", "default": False,
         "help": "Move this job ahead of normal queued work when resources become available."},
    ])
    if profile.get("supports_seed", True) and not any(field['key'] == 'seed' for field in fields):
        fields.insert(
            2 if len(profile["modes"]) > 1 else 1,
            {"key": "seed", "label": "Seed (0 for random)", "default": 0, "min": 0, "max": 2147483647,
             "help": "Use a nonzero seed for reproducibility. Zero chooses a fresh random seed."},
        )
    if profile.get("terms"):
        fields.append({"key": "accept_terms", "label": "I have reviewed and accept Google's published AlphaFold3 weight terms",
                       "type": "bool", "default": False, "required": True})
    config = {"fields": fields, "samples": int(profile.get("default_samples", profile.get("max_samples", 5))),
              "sample_max": int(profile.get("max_samples", 5)),
              "fixed_settings": list(profile.get("fixed_settings", ())),
              "sample_min_modes": {"exact": 2, "fast": 2} if profile.get("input_kind") == "boltzgen" else {}}
    if profile.get("parameter_presets"):
        config["parameter_presets"] = profile["parameter_presets"]
    return config
