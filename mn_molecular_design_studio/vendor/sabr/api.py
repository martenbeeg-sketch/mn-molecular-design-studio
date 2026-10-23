"""The public, in-memory SAbR API."""

import logging
from dataclasses import dataclass

from . import constants
from .alignment import align, load_references
from .model import encode
from .numbering import number_alignment
from .structure import apply_numbering, extract_chain

LOGGER = logging.getLogger(__name__)


def _normalize_choice(value: str, name: str, choices: tuple[str, ...]) -> str:
    """Normalize a case-insensitive string constrained to known choices."""
    if not isinstance(value, str) or value.lower() not in choices:
        raise ValueError(f"{name} must be one of {', '.join(choices)}.")
    return value.lower()


def _normalize_chain_type(
    chain_type: str,
    noise_level: float,
    mode: str,
) -> str:
    """Normalize a TCR type or antibody domain candidate list."""
    if not isinstance(chain_type, str):
        raise ValueError("chain_type must be a supported type or 'auto'.")

    normalized = chain_type.strip().upper()
    if normalized == "AUTO":
        return "auto"
    if normalized in constants.TCR_CHAIN_TYPES:
        return normalized

    reference_types = tuple(load_references(noise_level, mode))
    candidates = sorted(
        {value.strip().upper() for value in chain_type.split(",")}
    )
    if any(
        not candidate or not set(candidate).issubset(reference_types)
        for candidate in candidates
    ):
        choices = ", ".join(("auto", *reference_types))
        tcr_choices = ", ".join(constants.TCR_CHAIN_TYPES)
        raise ValueError(
            "chain_type must be a TCR type or comma-separated antibody "
            f"representations built from {choices}; TCR types are "
            f"{tcr_choices}."
        )
    return ",".join(candidates)


def _normalize_noise_level(noise_level: float) -> float:
    if isinstance(noise_level, bool):
        raise ValueError(
            f"noise_level must be one of {constants.NOISE_LEVELS}."
        )
    try:
        normalized_noise = float(noise_level)
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"noise_level must be one of {constants.NOISE_LEVELS}."
        ) from error
    if normalized_noise not in constants.NOISE_LEVELS:
        raise ValueError(
            f"noise_level must be one of {constants.NOISE_LEVELS}."
        )
    return normalized_noise


def _validate_residue_range(
    residue_range: tuple[int, int] | None,
) -> None:
    """Validate inclusive structure residue-number bounds.

    The bounds are compared with BioPython ``residue.id[1]``. They are residue
    numbers from the input structure, not zero-based sequence or array indices,
    and need not begin at one or be contiguous. Every insertion-code variant
    sharing a selected residue number is included.
    """
    if residue_range is None:
        return
    error = "residue_range must be an inclusive (start, end) tuple."
    if not isinstance(residue_range, tuple):
        raise ValueError(error)
    if len(residue_range) != 2:
        raise ValueError(error)
    if not all(
        isinstance(value, int) and not isinstance(value, bool)
        for value in residue_range
    ):
        raise ValueError(error)
    if residue_range[0] > residue_range[1]:
        raise ValueError("residue_range start must not exceed its end.")


@dataclass(slots=True)
class _RenumberOptions:
    """Validated and normalized arguments for one renumbering operation."""

    scheme: str
    chain_type: str
    noise_level: float
    residue_range: tuple[int, int] | None
    mode: str
    scfv: bool

    def __post_init__(self) -> None:
        self.scheme = _normalize_choice(
            self.scheme, "scheme", constants.NUMBERING_SCHEMES
        )
        self.noise_level = _normalize_noise_level(self.noise_level)
        self.mode = _normalize_choice(self.mode, "mode", constants.MODES)
        self.chain_type = _normalize_chain_type(
            self.chain_type,
            self.noise_level,
            self.mode,
        )
        _validate_residue_range(self.residue_range)
        if (
            self.chain_type in constants.TCR_CHAIN_TYPES
            and self.scheme not in ("imgt", "aho")
        ):
            raise ValueError(
                "TCR chain types support only IMGT or AHo numbering."
            )
        if not isinstance(self.scfv, bool):
            raise ValueError("scfv must be a boolean.")
        if self.scfv and self.chain_type != "auto":
            raise ValueError("scfv requires chain_type='auto'.")
        if self.scfv:
            self.chain_type = ",".join(constants.SCFV_CHAIN_TYPES)


def renumber_structure(
    structure,
    chain: str,
    scheme: str = "imgt",
    chain_type: str = "auto",
    noise_level: float = 0.0,
    residue_range: tuple[int, int] | None = None,
    mode: str = "sabr",
    scfv: bool = False,
):
    """Return a renumbered copy of a Biopython structure.

    ``mode="softalign"`` selects the SoftAlign encoder, references, and gap
    penalties as one scientifically consistent parameter set. ``chain_type``
    accepts a comma-separated set of antibody domain representations such
    as ``"H,K,HK,HL"``, or one TCR type from ``"A,B,G,D"``. TCR types align
    only against the K reference. ``scfv=True`` is equivalent to
    ``chain_type="HK,HL,KH,LH"``. Structural gaps are detected automatically
    and waive reference-deletion opening costs at the corresponding boundary.
    """
    options = _RenumberOptions(
        scheme=scheme,
        chain_type=chain_type,
        noise_level=noise_level,
        residue_range=residue_range,
        mode=mode,
        scfv=scfv,
    )
    data = extract_chain(structure, chain, options.residue_range)
    embeddings = encode(data.coords, options.mode)
    imgt_alignment, selected_type, score = align(
        embeddings,
        data.gap_indices,
        options.chain_type,
        options.noise_level,
        mode=options.mode,
        scfv=options.scfv,
    )
    LOGGER.info(
        "Selected %s reference with alignment score %.4f.",
        selected_type,
        score,
    )
    numbering_type = (
        options.chain_type
        if options.chain_type in constants.TCR_CHAIN_TYPES
        else selected_type
    )
    numbered = number_alignment(
        imgt_alignment,
        data.sequence,
        options.scheme,
        numbering_type,
    )
    return apply_numbering(structure, chain, data, numbered)


def number_chain(
    structure,
    chain: str,
    scheme: str = "imgt",
    chain_type: str = "auto",
    noise_level: float = 0.0,
    residue_range: tuple[int, int] | None = None,
    mode: str = "sabr",
    scfv: bool = False,
):
    """Return residue-indexed numbering without changing the structure.

    Each item is ``(chain_index, number, insertion_code, amino_acid)``.
    ``chain_index`` is zero-based in the original Biopython chain and remains
    correct when non-polymer residues are interspersed in the input.
    """
    options = _RenumberOptions(
        scheme=scheme,
        chain_type=chain_type,
        noise_level=noise_level,
        residue_range=residue_range,
        mode=mode,
        scfv=scfv,
    )
    data = extract_chain(structure, chain, options.residue_range)
    embeddings = encode(data.coords, options.mode)
    imgt_alignment, selected_type, score = align(
        embeddings,
        data.gap_indices,
        options.chain_type,
        options.noise_level,
        mode=options.mode,
        scfv=options.scfv,
    )
    LOGGER.info(
        "Selected %s reference with alignment score %.4f.",
        selected_type,
        score,
    )
    numbering_type = (
        options.chain_type
        if options.chain_type in constants.TCR_CHAIN_TYPES
        else selected_type
    )
    numbered = number_alignment(
        imgt_alignment,
        data.sequence,
        options.scheme,
        numbering_type,
    )
    return [
        (data.residue_indices[row], number, insertion_code, amino_acid)
        for row, number, insertion_code, amino_acid in numbered
    ]


def cdr_residue_indices(
    structure,
    chain: str,
    chain_type: str = "auto",
    noise_level: float = 0.0,
    residue_range: tuple[int, int] | None = None,
    mode: str = "sabr",
    scfv: bool = False,
) -> dict[str, list[int]]:
    """Identify zero-based original-chain indices in IMGT CDR1, CDR2, CDR3.

    The returned indices address the input structure's chain residues, not
    the filtered protein sequence, so callers can safely preserve the amino
    acid identity of selected residues in downstream design steps.
    """
    regions = {name: [] for name in constants.IMGT_LOOPS}
    for chain_index, number, _insertion_code, _amino_acid in number_chain(
        structure,
        chain,
        scheme="imgt",
        chain_type=chain_type,
        noise_level=noise_level,
        residue_range=residue_range,
        mode=mode,
        scfv=scfv,
    ):
        domain_position = (number - 1) % constants.DOMAIN_NUMBERING_STRIDE + 1
        for name, (start, end) in constants.IMGT_LOOPS.items():
            if start <= domain_position <= end:
                regions[name].append(chain_index)
    return regions


def cdr_residue_records(
    structure,
    chain: str,
    chain_type: str = "auto",
    noise_level: float = 0.0,
    residue_range: tuple[int, int] | None = None,
    mode: str = "sabr",
    scfv: bool = False,
) -> dict[str, list[tuple[int, int, str, str]]]:
    """Return CDR residue identities with original and IMGT indices.

    Each record is ``(chain_index, imgt_position, insertion_code,
    amino_acid)``. This lets downstream motif-guidance code locate a source
    residue in the original structure while carrying its exact amino-acid
    identity into a fixed-position sequence-design mask.
    """
    regions = {name: [] for name in constants.IMGT_LOOPS}
    for chain_index, number, insertion_code, amino_acid in number_chain(
        structure,
        chain,
        scheme="imgt",
        chain_type=chain_type,
        noise_level=noise_level,
        residue_range=residue_range,
        mode=mode,
        scfv=scfv,
    ):
        domain_position = (number - 1) % constants.DOMAIN_NUMBERING_STRIDE + 1
        for name, (start, end) in constants.IMGT_LOOPS.items():
            if start <= domain_position <= end:
                regions[name].append(
                    (chain_index, domain_position, insertion_code, amino_acid)
                )
    return regions
