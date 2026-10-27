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
