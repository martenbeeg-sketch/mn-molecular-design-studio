"""Prediction configuration semantics shared by editors and execution adapters."""
AF2_MODES = {
    'conditioned_self_consistency': 'Conditioned self-consistency: designed target and binder coordinates initialize AF2. Agreement is conditioned on the starting pose, not independent binding-orientation verification.',
    'target_only': 'Target-only template: target coordinates, sequence and side chains are retained; binder sequence is supplied without binder coordinates or interface pose. Target conformation and optional target MSA remain conditioning. Confidence and agreement are model evidence, not experimental binding probabilities.',
}


def af2_mode(mode=None, use_initial_guess=None, *, legacy_initial_guess=False):
    if mode is None:
        mode = 'conditioned_self_consistency' if (legacy_initial_guess if use_initial_guess is None else use_initial_guess) else 'target_only'
    if mode not in AF2_MODES:
        raise ValueError('Unsupported AF2 validation mode.')
    initial_guess = mode == 'conditioned_self_consistency'
    if use_initial_guess is not None and bool(use_initial_guess) != initial_guess:
        raise ValueError('Conflicting AF2 validation mode and use_initial_guess settings.')
    return mode, initial_guess


def boltz_settings(settings, *, saved_graph=False):
    """Canonical aliases, with explicit handling of the shipped inactive controls.

    Unversioned canvas graphs historically used samples=1/recycles=3 and
    unconditionally emitted full PAE. Preserve that execution on migration.
    Changed, contradictory legacy controls require explicit user resolution.
    Other clients may use aliases directly; conflicting values are rejected.
    """
    result = dict(settings)
    notices = list(result.get('configuration_migration', []))
    old_graph = saved_graph and result.get('prediction_config_version') != 2
    for alias, canonical, shipped, historical in (
        ('diffusion_samples', 'samples', 3, 1), ('recycling_steps', 'num_recycles', 10, 3),
    ):
        if alias not in result:
            continue
        value = result.pop(alias)
        if canonical in result and result[canonical] != value:
            if not (old_graph and value == shipped):
                raise ValueError(f'Conflicting Boltz-2 settings: {canonical}={result[canonical]} and {alias}={value}. Use only {canonical}. Legacy canvas aliases were ignored; choose the intended effective value explicitly.')
            notices.append(f'Removed inactive shipped {alias}={value}; retained effective {canonical}={result[canonical]}.')
        elif old_graph and canonical not in result:
            result[canonical] = historical
            notices.append(f'Legacy canvas {alias}={value} was inactive; retained historical effective {canonical}={historical}. Set {canonical} explicitly to change it.')
        else:
            result[canonical] = value
    if old_graph and result.get('write_full_pae') is False:
        result['write_full_pae'] = True
        notices.append('Legacy write_full_pae=False was inactive; retained the historical full PAE output. Full PAE is fixed by the installed engine; the inactive toggle has been removed.')
    if result.get('write_full_pae', True) is not True:
        raise ValueError('Full PAE is always emitted by the installed Boltz-2 engine and retained for ipSAE; write_full_pae=False is unsupported.')
    if saved_graph:
        result['prediction_config_version'] = 2
        if notices:
            result['configuration_migration'] = list(dict.fromkeys(notices))
    return result
