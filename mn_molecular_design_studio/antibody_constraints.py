"""Pure sequence checks shared with the staged antibody design runner."""


def verify_masked_sequence(sequence: str, prompt: str) -> dict:
    """Only # positions may differ from the native generator's sampled prompt."""
    if len(sequence) != len(prompt):
        raise ValueError('VHH sequence length differs from its sampled framework/CDR prompt.')
    if set(sequence) - set('ACDEFGHIKLMNPQRSTVWY'):
        raise ValueError('VHH generator returned unsupported amino-acid identities.')
    changed = [i + 1 for i, (actual, expected) in enumerate(zip(sequence, prompt))
               if expected != '#' and actual != expected]
    if changed:
        raise ValueError(f'VHH generator changed fixed framework residues: {changed}.')
    return {'status': 'verified', 'fixed_position_count': len(prompt) - prompt.count('#'),
            'mutable_positions': [i + 1 for i, expected in enumerate(prompt) if expected == '#']}


def verify_fixed_sequence(sequence: str, labels: list[str], original: str, fixed: list[str]) -> dict:
    if len(sequence) != len(original) or len(sequence) != len(labels):
        raise ValueError('Redesigned sequence length differs from its input structure.')
    unknown = set(fixed) - set(labels)
    if unknown:
        raise ValueError(f'Fixed residues are absent from the designed chains: {sorted(unknown)}.')
    changed = [label for label, before, after in zip(labels, original, sequence)
               if label in set(fixed) and before != after]
    if changed:
        raise ValueError(f'Sequence redesign changed fixed residues: {changed}.')
    return {'status': 'verified', 'fixed_position_count': len(set(fixed))}
