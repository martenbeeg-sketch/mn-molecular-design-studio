"""Data model shared by all readers and by everything downstream of them.

A reader (BindCraft 1 today, BindCraft 2 later) turns one run folder into a RunData. Later stages
(structure primitives, target frame, interface analysis, families, specificity, export) read RunData
and never re-open the CSVs.

Conventions
-----------
* Design ids are file stems without the _model{k} suffix (binder01_l77_s305024_mpnn4); backbone ids drop
  _mpnn{n} (binder01_l77_s305024). A backbone is one BindCraft trajectory and is the unit of a family.
  Ids are unique inside a run only: the global key of a design is (run_id, id), see DesignRecord.key.
* Metrics are stored with BindCraft's own column stems in the scale of the CSV, never rescaled.
  The scale of every stem is in METRIC_INFO. Missing values are None, not NaN.
* Residue numbers are construct numbers: the numbering BindCraft writes into its PDB files (target
  chain A, numbers 1..N as in the file handed to BindCraft, with gaps if a run deleted loops from the
  target). Binder residues (chain B) are numbered 1..Length. Mapping to UniProt or reference-structure numbers is the
  target frame's job (TargetConstruct.uniprot_offset is filled in there; uniprot = construct + offset,
  and the blobs' input_offset is its negative). Everything that crosses into the blobs is UniProt
  numbered (docs/CONTRACTS.md section 2); construct numbers stay inside RunData and per-design analysis.
* Paths are absolute strings.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from typing import Any, Iterable, Optional

from .config import FilterFailure, FilterRule

SCHEMA_VERSION = 1

# ---------------------------------------------------------------------------------------------
# Shared analysis constants (docs/CONTRACTS.md section 2). One definition for every module, so
# interface, families, modes and target_frame cannot disagree on what a contact is.
# ---------------------------------------------------------------------------------------------

#: Heavy-atom distance in angstrom for "residue is in contact" on the pages (design.tc, bc, nc,
#: contacts[], hot, footprints for families and modes). The mockup convention (data_contract.md).
CONTACT_CUTOFF = 4.5
#: BindCraft's own interface rule (any atom, hydrogens included, within 4.0 A); only used to compare
#: with the InterfaceResidues column, never for the pages.
BINDCRAFT_INTERFACE_CUTOFF = 4.0
HBOND_MAX_DIST = 3.5
SALT_BRIDGE_MAX_DIST = 4.0
#: A hotspot counts as touched when a binder heavy atom is within CONTACT_CUTOFF of it.
HOTSPOT_TOUCH_CUTOFF = CONTACT_CUTOFF
#: Near-miss rule of the Query page and the manifest (thresholds.near_miss): exactly one distinct
#: failed stem, missing its threshold by at most this share of the threshold.
NEAR_MISS_MAX_REL = 0.25
NEAR_MISS_MAX_STEMS = 1

# ---------------------------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------------------------

#: BindCraft's metric stems. In final_design_stats.csv and mpnn_design_stats.csv every stem exists as
#: Average_<stem> and 1_<stem> .. 5_<stem>; trajectory_stats.csv has one plain column per stem.
METRIC_STEMS: tuple[str, ...] = (
    "pLDDT", "pTM", "i_pTM", "pAE", "i_pAE", "i_pLDDT", "ss_pLDDT",
    "Unrelaxed_Clashes", "Relaxed_Clashes", "Binder_Energy_Score", "Surface_Hydrophobicity",
    "ShapeComplementarity", "PackStat", "dG", "dSASA", "dG/dSASA", "Interface_SASA_%",
    "Interface_Hydrophobicity", "n_InterfaceResidues", "n_InterfaceHbonds", "InterfaceHbondsPercentage",
    "n_InterfaceUnsatHbonds", "InterfaceUnsatHbondsPercentage",
    "Interface_Helix%", "Interface_BetaSheet%", "Interface_Loop%",
    "Binder_Helix%", "Binder_BetaSheet%", "Binder_Loop%",
    "Hotspot_RMSD", "Target_RMSD", "Binder_pLDDT", "Binder_pTM", "Binder_pAE", "Binder_RMSD",
    "MPNN_score", "MPNN_seq_recovery",
)

#: Stems that exist once per sequence rather than once per model (no Average_/<k>_ columns).
SINGLE_VALUE_STEMS: tuple[str, ...] = ("MPNN_score", "MPNN_seq_recovery")


@dataclass(frozen=True)
class MetricInfo:
    """Scale and meaning of one metric stem, exactly as it appears in the CSV."""

    scale: str   # "0-1", "0-100 percent", "count", "angstrom", "REU", "A2", "PAE/31", "x100"
    better: str  # "higher", "lower", "n/a"
    note: str = ""


METRIC_INFO: dict[str, MetricInfo] = {
    "pLDDT": MetricInfo("0-1", "higher", "complex prediction, mean over binder residues. The PDB B-factor of chain B is this value times 100."),
    "pTM": MetricInfo("0-1", "higher"),
    "i_pTM": MetricInfo("0-1", "higher", "interface pTM of the AlphaFold2 complex prediction"),
    "pAE": MetricInfo("PAE/31", "lower", "predicted aligned error divided by 31 (0.35 = 10.85 A)"),
    "i_pAE": MetricInfo("PAE/31", "lower", "inter-chain predicted aligned error divided by 31"),
    "i_pLDDT": MetricInfo("0-1", "higher", "pLDDT of interface residues"),
    "ss_pLDDT": MetricInfo("0-1", "higher", "pLDDT of helix and strand residues"),
    "Unrelaxed_Clashes": MetricInfo("count", "lower", "average over models can be fractional"),
    "Relaxed_Clashes": MetricInfo("count", "lower", "after PyRosetta relax"),
    "Binder_Energy_Score": MetricInfo("REU", "lower", "Rosetta score of the binder alone"),
    "Surface_Hydrophobicity": MetricInfo("0-1", "lower", "fraction of surface residues that are hydrophobic"),
    "ShapeComplementarity": MetricInfo("0-1", "higher"),
    "PackStat": MetricInfo("0-1", "higher"),
    "dG": MetricInfo("REU", "lower", "Rosetta interface energy; more negative is better"),
    "dSASA": MetricInfo("A2", "higher", "buried surface area on binding"),
    "dG/dSASA": MetricInfo("x100", "lower", "dG divided by dSASA, times 100 (REU per 100 A2)"),
    "Interface_SASA_%": MetricInfo("0-100 percent", "n/a", "interface share of the binder surface"),
    "Interface_Hydrophobicity": MetricInfo("0-100 percent", "n/a", "percent hydrophobic interface residues"),
    "n_InterfaceResidues": MetricInfo("count", "n/a", "binder-side interface residues; average over models can be fractional"),
    "n_InterfaceHbonds": MetricInfo("count", "higher"),
    "InterfaceHbondsPercentage": MetricInfo("0-100 percent", "higher"),
    "n_InterfaceUnsatHbonds": MetricInfo("count", "lower", "a raw count: grows with interface size"),
    "InterfaceUnsatHbondsPercentage": MetricInfo("0-100 percent", "lower"),
    "Interface_Helix%": MetricInfo("0-100 percent", "n/a"),
    "Interface_BetaSheet%": MetricInfo("0-100 percent", "n/a"),
    "Interface_Loop%": MetricInfo("0-100 percent", "n/a"),
    "Binder_Helix%": MetricInfo("0-100 percent", "n/a"),
    "Binder_BetaSheet%": MetricInfo("0-100 percent", "n/a"),
    "Binder_Loop%": MetricInfo("0-100 percent", "n/a"),
    "Hotspot_RMSD": MetricInfo("angstrom", "lower", "hotspot residues' position after superposing the target (not always filtered)"),
    "Target_RMSD": MetricInfo("angstrom", "lower", "predicted target versus the input target"),
    "Binder_pLDDT": MetricInfo("0-1", "higher", "binder predicted alone"),
    "Binder_pTM": MetricInfo("0-1", "higher", "binder predicted alone"),
    "Binder_pAE": MetricInfo("PAE/31", "lower", "binder predicted alone"),
    "Binder_RMSD": MetricInfo("angstrom", "lower", "binder alone versus binder in the complex"),
    "MPNN_score": MetricInfo("ProteinMPNN score", "lower", "one value per sequence"),
    "MPNN_seq_recovery": MetricInfo("0-1", "n/a", "one value per sequence; identity to the trajectory sequence"),
}

#: BindCraft stem -> key of design.m on the pages, for stems whose page value IS BindCraft's own
#: Average_<stem> (CSV scale, nothing recomputed). i_pTM feeds both iptm and iptm_b.
PAGE_METRIC_COPIED: dict[str, str] = {
    "pLDDT": "plddt", "pTM": "ptm", "i_pTM": "iptm_b", "i_pAE": "ipae_norm", "i_pLDDT": "ciplddt",
    "ShapeComplementarity": "sc", "dG": "dg", "Surface_Hydrophobicity": "surf_hyd",
    "Binder_RMSD": "binder_rmsd", "Hotspot_RMSD": "hotspot_rmsd", "Target_RMSD": "target_rmsd",
}

#: BindCraft stem -> key of design.m whose page value is RECOMPUTED from the structure (interface
#: module). BindCraft's number stays in design.src; the agreement goes to
#: manifest.validation.metric_agreement. Percent stems are 0-100 in the CSV, the page keys 0-1.
PAGE_METRIC_RECOMPUTED: dict[str, str] = {
    "dSASA": "dsasa", "n_InterfaceHbonds": "n_hbonds", "n_InterfaceUnsatHbonds": "n_unsat",
    "Interface_Loop%": "loop_frac", "Interface_Helix%": "iface_helix",
    "Binder_Helix%": "helix_frac", "Binder_BetaSheet%": "strand_frac",
}

PAGE_METRIC: dict[str, str] = {**PAGE_METRIC_COPIED, **PAGE_METRIC_RECOMPUTED}


def page_metric(stem: str) -> Optional[str]:
    """Key of design.m that carries the same quantity as a BindCraft stem, or None (src only)."""
    return PAGE_METRIC.get(stem)


def filter_rule_to_dict(rule: FilterRule) -> dict[str, Any]:
    """One entry of campaign runs[].filters (data_contract.md 3.3)."""
    return {
        "stem": rule.stem, "metric": page_metric(rule.stem), "threshold": rule.threshold, "op": rule.op,
        "per_model": rule.per_model, "provenance": rule.provenance, "note": rule.note,
    }


def failure_to_dict(f: FilterFailure) -> dict[str, Any]:
    """One entry of design.fail[] (data_contract.md 3.5). value and threshold are in the CSV scale and
    value is BindCraft's number (not the recomputed page metric); rel is gap / |threshold|."""
    return {
        "stem": f.stem, "metric": page_metric(f.stem), "scope": f.scope, "value": f.value,
        "threshold": f.threshold, "op": ">=" if f.higher_is_better else "<=", "gap": f.gap,
        "rel": round(f.gap / max(abs(f.threshold), 1e-9), 3),
    }


def near_miss_record(failures: Iterable[FilterFailure]) -> Optional[dict[str, Any]]:
    """design.nm of a rejected design (data_contract.md 3.5), None when it fails nothing.

    n_failed counts distinct stems (Average and per-model scopes of one stem count once), worst_rel is
    the largest rel, score = max(0, 1 - worst_rel) / n_failed, near = exactly NEAR_MISS_MAX_STEMS stem
    missed by at most NEAR_MISS_MAX_REL of its threshold.
    """
    recs = [failure_to_dict(f) for f in failures]
    if not recs:
        return None
    stems = list(dict.fromkeys(r["stem"] for r in recs))
    worst = max(r["rel"] for r in recs)
    return {
        "n_failed": len(stems), "stems": sorted(stems), "worst_rel": worst,
        "score": round(max(0.0, 1.0 - worst) / len(stems), 3),
        "near": len(stems) == NEAR_MISS_MAX_STEMS and worst <= NEAR_MISS_MAX_REL,
    }


#: Scale note shown wherever raw metrics are listed.
SCALE_NOTE = (
    "Metrics keep BindCraft's CSV scales: pLDDT, pTM, i_pTM, i_pLDDT and similar are 0-1 (PDB B-factors "
    "are 0-100); pAE and i_pAE are PAE divided by 31; hydrophobicity, SASA and secondary-structure "
    "shares ending in % are 0-100; RMSD in angstrom; dG in Rosetta energy units; dG/dSASA is times 100."
)


def metric_scale(stem: str) -> Optional[MetricInfo]:
    return METRIC_INFO.get(stem)


# ---------------------------------------------------------------------------------------------
# Records
# ---------------------------------------------------------------------------------------------


@dataclass
class InputFile:
    """One input file the run was built from, with a hash so a bundle can prove its provenance."""

    role: str  # "settings", "csv", "xlsx", "accepted_structure", "rejected_structure", "note", ...
    path: str
    sha256: str
    size: int


@dataclass
class RunSettings:
    """Settings of one run. Per-run: two runs of the same target can differ in everything here.

    Hotspot numbers are construct numbers (the numbering of the file given to BindCraft).
    *_source fields say where a value came from: "settings_json", "csv" (read from the CSV columns),
    "inferred" (derived from the data) or "unknown".
    """

    run_id: str
    label: str
    binder_name: Optional[str] = None
    protocol: Optional[str] = None            # "4stage"
    advanced_settings_name: Optional[str] = None  # basename recorded in the CSVs, e.g. default_4stage_multimer
    filters_name: Optional[str] = None        # e.g. default_filters
    target_settings_name: Optional[str] = None  # e.g. target_run01 (the settings file's basename)
    helicity: Optional[float] = None          # BindCraft helicity bias (-0.3 default)
    starting_pdb_name: Optional[str] = None   # basename recorded in the settings file
    starting_pdb_recorded: Optional[str] = None  # full path as recorded (on the machine that ran BindCraft)
    design_path_recorded: Optional[str] = None
    chains: Optional[str] = None              # target chains given to BindCraft ("A")
    hotspots: list[int] = field(default_factory=list)  # construct numbers
    hotspot_string: Optional[str] = None      # as written, e.g. "26,122,125"
    hotspot_source: str = "unknown"
    length_min: Optional[int] = None
    length_max: Optional[int] = None
    length_source: str = "unknown"
    final_designs_target: Optional[int] = None
    final_designs_target_source: str = "unknown"
    settings_path: Optional[str] = None
    settings_raw: dict[str, Any] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)


@dataclass
class TargetConstruct:
    """The target as BindCraft saw it in this run, derived by reading the run's structures.

    numbers/sequence list the residues actually present in chain A of the structures (all structures
    of a run are expected to carry the identical set; `consistent` says whether they do). deleted_ranges
    are runs of numbers missing inside the span first..last (or missing against the reference construct
    when one was given); they are removed regions for this run and extra cut points for the cut-end check.
    """

    run_id: str
    chain_id: str = "A"
    binder_chain_id: str = "B"
    numbers: list[int] = field(default_factory=list)
    sequence: str = ""                        # one-letter, same order as numbers
    first: Optional[int] = None
    last: Optional[int] = None
    deleted_ranges: list[tuple[int, int]] = field(default_factory=list)
    reference_path: Optional[str] = None      # construct file the numbers were compared with
    reference_length: Optional[int] = None
    sequence_matches_reference: Optional[bool] = None
    consistent: Optional[bool] = None         # every checked structure has the identical residue set
    n_structures_checked: int = 0
    hotspot_residues: dict[int, str] = field(default_factory=dict)  # construct number -> one-letter code
    # Filled by the target frame (not by readers): construct number i is UniProt i + uniprot_offset.
    uniprot_accession: Optional[str] = None
    uniprot_offset: Optional[int] = None

    @property
    def n_present(self) -> int:
        return len(self.numbers)

    def residue_at(self, number: int) -> Optional[str]:
        try:
            return self.sequence[self.numbers.index(number)]
        except ValueError:
            return None

    # -- numbering ---------------------------------------------------------------------------
    # uniprot = construct + uniprot_offset (e.g. +149). The blobs say it the other way round:
    # orig = uniprot + input_offset, so input_offset = -uniprot_offset (e.g. -149).
    @property
    def input_offset(self) -> Optional[int]:
        """Sign convention of the data contract (target.numbering.input_offset): orig = uniprot + input_offset."""
        return None if self.uniprot_offset is None else -self.uniprot_offset

    def to_uniprot(self, construct_number: int) -> int:
        """UniProt number of a construct number. Needs uniprot_offset (filled by the target frame)."""
        if self.uniprot_offset is None:
            raise ValueError(f"run {self.run_id}: uniprot_offset is not set (target frame not applied)")
        return construct_number + self.uniprot_offset

    def from_uniprot(self, uniprot_number: int) -> int:
        """Construct number of a UniProt number (may be a number the construct does not contain)."""
        if self.uniprot_offset is None:
            raise ValueError(f"run {self.run_id}: uniprot_offset is not set (target frame not applied)")
        return uniprot_number - self.uniprot_offset

    def uniprot_numbers(self) -> list[int]:
        """UniProt numbers of the residues present, ascending (gaps stay gaps)."""
        return [self.to_uniprot(n) for n in sorted(self.numbers)]

    def resi_ranges(self) -> list[tuple[int, int]]:
        """Contiguous stretches of present residues in construct numbers, ascending."""
        out: list[tuple[int, int]] = []
        for n in sorted(self.numbers):
            if out and n == out[-1][1] + 1:
                out[-1] = (out[-1][0], n)
            else:
                out.append((n, n))
        return out

    def cut_points(self) -> list[int]:
        """Residues next to a cut, construct numbers: first, last and both flanks of every removed
        stretch. Whether first and last are real cuts depends on the native chain, which only the
        target frame knows; it drops them when the construct reaches a native end."""
        pts: list[int] = []
        for a, b in self.resi_ranges():
            pts.extend((a, b))
        return sorted(set(pts))


@dataclass
class DesignRecord:
    """One designed binder sequence (one row of mpnn_design_stats.csv).

    status "accepted": passed BindCraft's filters, has a row in final_design_stats.csv and a structure
    in accepted/. status "rejected": has a row in mpnn_design_stats.csv only (and a structure in
    rejected/ if the folder exists; structure_path is None otherwise).

    metrics: Average_<stem> values; models: {model number: {stem: value}}, only models that carry
    values (BindCraft 1 default: 1 and 2). Values are in CSV scales (see METRIC_INFO).
    structure_model: k of the _model{k} suffix = the AlphaFold2 model whose relaxed pose was saved
    (the one with the highest complex pLDDT, ties to model 1). Metrics of that model describe the file;
    interface_residues is written by BindCraft from model 2 whatever k is, so for k = 1 it describes
    a different pose than the file.
    """

    id: str
    run_id: str
    backbone_id: str
    status: str                               # "accepted" | "rejected"
    binder_name: str
    length: int
    seed: int
    mpnn_index: Optional[int]
    sequence: str                             # binder, one-letter. Never paste into docs; use ids.
    rank: Optional[int] = None                # accepted designs only (Rank column)
    protocol: Optional[str] = None
    target_hotspot: Optional[str] = None      # as written in the CSV, construct numbers
    structure_path: Optional[str] = None
    structure_model: Optional[int] = None
    ranked_path: Optional[str] = None         # accepted/Ranked/<rank>_<name>_model<k>.pdb
    metrics: dict[str, Optional[float]] = field(default_factory=dict)
    models: dict[int, dict[str, Optional[float]]] = field(default_factory=dict)
    interface_residues: list[str] = field(default_factory=list)  # "B13": chain letter plus binder residue number
    interface_aas: dict[str, dict[str, float]] = field(default_factory=dict)  # "avg", "1", "2" -> letter -> count
    design_time_s: Optional[float] = None
    notes: str = ""
    failed_filters: list[FilterFailure] = field(default_factory=list)  # under DEFAULT_FILTERS
    xlsx_marked: bool = False                 # row highlighted in the run's xlsx (meaning unconfirmed)
    # Source columns that have no field above (BindCraft 2: unknown columns, failed_filters names, ...).
    extras: dict[str, Any] = field(default_factory=dict)

    @property
    def key(self) -> tuple[str, str]:
        return (self.run_id, self.id)

    @property
    def has_structure(self) -> bool:
        return self.structure_path is not None

    def metric(self, stem: str, model: Optional[int] = None) -> Optional[float]:
        """Average value of a stem, or the value of one model."""
        if model is None:
            return self.metrics.get(stem)
        return self.models.get(model, {}).get(stem)

    def pose_metric(self, stem: str) -> Optional[float]:
        """Value from the model whose pose is in the structure file (average if it is unknown)."""
        if self.structure_model is not None:
            v = self.models.get(self.structure_model, {}).get(stem)
            if v is not None:
                return v
        return self.metrics.get(stem)

    @property
    def failed_stems(self) -> list[str]:
        seen: dict[str, None] = {}
        for f in self.failed_filters:
            seen.setdefault(f.stem, None)
        return list(seen)


@dataclass
class BackboneRecord:
    """One trajectory (one row of trajectory_stats.csv): a backbone and the family it founds.

    order is the 0-based row position in trajectory_stats.csv, which is chronological. The metrics are
    those of the trajectory stage (single values, no per-model columns).
    """

    id: str
    run_id: str
    binder_name: str
    length: int
    seed: int
    order: int
    trajectory_metrics: dict[str, Optional[float]] = field(default_factory=dict)
    trajectory_time_s: Optional[float] = None
    start_sequence: Optional[str] = None      # sequence at the end of the trajectory, before ProteinMPNN
    interface_residues: list[str] = field(default_factory=list)
    notes: str = ""
    design_ids: list[str] = field(default_factory=list)       # every mpnn_design_stats row of this backbone
    accepted_ids: list[str] = field(default_factory=list)
    extras: dict[str, Any] = field(default_factory=dict)     # source columns without a field (BindCraft 2)

    @property
    def key(self) -> tuple[str, str]:
        return (self.run_id, self.id)

    @property
    def traj(self) -> int:
        """Generation order as the blobs carry it (design.traj): 1 = first row of trajectory_stats.csv."""
        return self.order + 1

    @property
    def n_accepted(self) -> int:
        return len(self.accepted_ids)


FUNNEL_STATUSES = ("observed", "derived", "unknown")


@dataclass
class FunnelCounts:
    """What can be counted for a run. None means unknown; status says observed, derived or unknown.

    BindCraft 1 keeps rows only for passing items; the Trajectory/ and MPNN/ folders that would give
    the number of attempted trajectories are absent from many runs, so trajectories_attempted is
    unknown there and trajectories_attempted_estimate is a derived figure that says how it was made.
    """

    trajectories_attempted: Optional[int] = None
    trajectories_attempted_estimate: Optional[int] = None
    trajectories_terminated: Optional[int] = None       # sum of the Trajectory_* failure counters
    trajectory_rows: Optional[int] = None               # rows in trajectory_stats.csv
    backbones_with_mpnn_rows: Optional[int] = None
    mpnn_rows: Optional[int] = None                     # sequences with a row in mpnn_design_stats.csv
    mpnn_rows_rejected: Optional[int] = None
    accepted_designs: Optional[int] = None
    accepted_backbones: Optional[int] = None            # families with at least one accepted design
    mpnn_sequences_per_trajectory: Optional[int] = None  # highest mpnn index seen (20 = BindCraft default)
    mpnn_sequences_scored_estimate: Optional[int] = None
    final_designs_target: Optional[int] = None
    status: dict[str, str] = field(default_factory=dict)
    notes: dict[str, str] = field(default_factory=dict)

    def stage_rows(self) -> list[dict[str, Any]]:
        """Ordered funnel rows for display: label, value (None if unknown), status, note."""
        spec = [
            ("trajectories_attempted", "Trajectories attempted"),
            ("trajectories_attempted_estimate", "Trajectories attempted (estimate)"),
            ("trajectory_rows", "Trajectories that completed and were scored"),
            ("backbones_with_mpnn_rows", "Backbones with at least one scored sequence"),
            ("mpnn_rows", "Sequences scored (accepted or not)"),
            ("accepted_designs", "Accepted designs"),
            ("accepted_backbones", "Accepted backbones (families)"),
        ]
        return [
            {
                "field": f, "label": label, "value": getattr(self, f),
                "status": self.status.get(f, "unknown" if getattr(self, f) is None else "observed"),
                "note": self.notes.get(f, ""),
            }
            for f, label in spec
        ]


def binder_chain_count(design: "DesignRecord") -> int:
    """How many chains this design's binder has (1 unless the reader found more).

    The BindCraft 2 reader counts them from Binder_Sequence, where the chains of
    a multi-chain binder are separated by "/", and records it in extras. Every
    caller that loads a structure needs it: structure.load_complex takes one
    binder chain unless told otherwise, and the rest of the molecule would
    simply be absent from whatever is measured.
    """
    extras = getattr(design, "extras", None) or {}     # duck-typed: termini takes objects that are not DesignRecords
    try:
        return max(1, int(extras.get("n_binder_chains", 1) or 1))
    except (AttributeError, TypeError, ValueError):
        return 1


@dataclass
class RunData:
    """Everything a reader knows about one run. This is the reader contract.

    designs holds every sequence with a row in mpnn_design_stats.csv (accepted and not accepted);
    rejected designs of a run without a rejected/ folder are metric-only (structure_path None).
    Use the helper methods instead of filtering by hand.
    """

    run_id: str
    label: str
    reader: str                               # "bindcraft1"
    source_dir: str
    settings: RunSettings
    target: TargetConstruct
    designs: list[DesignRecord] = field(default_factory=list)
    backbones: list[BackboneRecord] = field(default_factory=list)
    funnel: FunnelCounts = field(default_factory=FunnelCounts)
    failure_counters: dict[str, int] = field(default_factory=dict)  # failure_csv.csv, a sequence counts once per filter it fails
    files: list[InputFile] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)   # things that look wrong; shown to the user
    notes: list[str] = field(default_factory=list)      # things that are only worth knowing
    annotations: dict[str, Any] = field(default_factory=dict)  # xlsx highlights, free-text notes found in the folder
    tool_version: Optional[str] = None
    experimental: bool = False                # True for readers written from documentation only
    schema_version: int = SCHEMA_VERSION

    # -- selections -----------------------------------------------------------------------
    def accepted(self) -> list[DesignRecord]:
        return sorted((d for d in self.designs if d.status == "accepted"), key=lambda d: (d.rank is None, d.rank or 0, d.id))

    def rejected(self) -> list[DesignRecord]:
        return [d for d in self.designs if d.status == "rejected"]

    def rejected_with_structure(self) -> list[DesignRecord]:
        return [d for d in self.designs if d.status == "rejected" and d.has_structure]

    def design(self, design_id: str) -> Optional[DesignRecord]:
        return next((d for d in self.designs if d.id == design_id), None)

    def backbone(self, backbone_id: str) -> Optional[BackboneRecord]:
        return next((b for b in self.backbones if b.id == backbone_id), None)

    def accepted_backbones(self) -> list[BackboneRecord]:
        return [b for b in self.backbones if b.accepted_ids]

    # -- metrics --------------------------------------------------------------------------
    def available_metrics(self, designs: Optional[Iterable[DesignRecord]] = None) -> set[str]:
        """Stems that have a value for at least one design (default: accepted designs).

        Pages switch a check off when its metric is not in this set, instead of guessing.
        """
        pool = list(designs) if designs is not None else self.accepted()
        out: set[str] = set()
        for d in pool:
            out.update(k for k, v in d.metrics.items() if v is not None)
        return out

    def caps(self) -> dict[str, bool]:
        """{stem: present} for every BindCraft stem, over the accepted designs."""
        have = self.available_metrics()
        return {s: s in have for s in METRIC_STEMS}

    def hotspot_labels(self, numbering: str = "construct") -> list[str]:
        """Hotspots as "<one-letter><number>". numbering="construct" (default, as written in the
        settings, e.g. R26) or "uniprot" (e.g. R201; needs target.uniprot_offset). The pages and
        every text the user reads use UniProt labels."""
        t = self.target
        if numbering == "uniprot":
            return [f"{t.hotspot_residues.get(n, '?')}{t.to_uniprot(n)}" for n in self.settings.hotspots]
        return [f"{t.hotspot_residues.get(n, '?')}{n}" for n in self.settings.hotspots]

    def hotspots_uniprot(self) -> list[int]:
        """Hotspots of this run in UniProt numbers (needs target.uniprot_offset)."""
        return [self.target.to_uniprot(n) for n in self.settings.hotspots]

    # -- serialisation --------------------------------------------------------------------
    def to_dict(self, include_sequences: bool = False) -> dict[str, Any]:
        """JSON-safe dict. Binder sequences are left out unless asked for."""
        out = _jsonable(dataclasses.asdict(self))
        if not include_sequences:
            for d in out["designs"]:
                d["sequence"] = None
            for b in out["backbones"]:
                b["start_sequence"] = None
        return out


def union_uniprot_residues(targets: Iterable[TargetConstruct]) -> list[int]:
    """Union of the residues of several runs' targets in UniProt numbers, ascending: campaign
    target.resi and homolog resis. Every target needs uniprot_offset."""
    out: set[int] = set()
    for t in targets:
        out.update(t.uniprot_numbers())
    return sorted(out)


def _jsonable(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, float) and obj != obj:  # NaN
        return None
    return obj
