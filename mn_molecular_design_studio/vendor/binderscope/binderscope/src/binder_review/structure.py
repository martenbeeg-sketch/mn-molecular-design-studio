"""Structure primitives: loading, geometry, contacts, SASA, secondary structure, hydrogen bonds, TM-align, sequence helpers.

Pure numpy / scipy / gemmi / tmtools (no PyRosetta, no DSSP or Foldseek binaries). See docs/structure_api.md.

Conventions
  * Every residue-level helper keys on the residue NUMBER, never on position in the table: BindCraft
    targets keep construct numbering with deleted loops, so numbers are not always contiguous.
  * Coordinates are float64 (N, 3) arrays in angstrom. A rigid transform is (R, t) with ``x' = x @ R.T + t``.
  * Chain B-factors hold pLDDT on a 0 to 100 scale for BindCraft outputs.
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import Callable, Iterable, Sequence

import numpy as np
from scipy.spatial import cKDTree

__all__ = [
    "StructureError", "ValidationIssue", "Atoms", "ResidueTable", "Complex", "StructureFile", "Link", "AtomRef",
    "load_complex", "read_structure", "deposited_ss", "to_pdb_lines", "write_pdb",
    "kabsch", "apply_transform", "invert_transform", "rmsd", "superpose_residues", "Superposition", "rmsd_by_number",
    "Contacts", "residue_contacts", "atom_contacts", "bindcraft_interface_residues", "format_bindcraft_interface",
    "SasaResult", "sasa", "sasa_atoms", "relative_sasa", "InterfaceSasa", "interface_dsasa", "MAX_ASA", "SASA_PRESETS",
    "assign_dssp", "ss3_from_dssp", "assign_ss", "ss_fractions",
    "HBonds", "hbonds", "SaltBridge", "salt_bridges",
    "TMResult", "tm_align", "TMMatrix", "tm_matrix", "tm_align_atoms",
    "THREE_TO_ONE", "ONE_TO_THREE", "three_to_one", "one_to_three", "KYTE_DOOLITTLE", "BINDCRAFT_APOLAR",
    "gravy", "net_charge", "isoelectric_point", "sc_neighbors", "surface_hydrophobicity", "sequence_hydrophobic_fraction",
]


class StructureError(ValueError):
    """Raised when a structure cannot be used (missing chain, no atoms, unreadable file)."""


# --------------------------------------------------------------------------------------------- constants
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
    "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
# modified residues that map to a standard parent (MSE -> MET is also coordinate-converted on load)
_MODIFIED = {"MSE": "MET", "SEC": "CYS", "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS",
             "HSP": "HIS", "CYX": "CYS", "CYM": "CYS", "ASH": "ASP", "GLH": "GLU", "LYN": "LYS"}
ONE_TO_THREE = {v: k for k, v in THREE_TO_ONE.items()}
AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
_WATER = {"HOH", "WAT", "DOD", "H2O", "TIP", "TIP3", "SOL"}

# Biopython's ShrakeRupley default radii (angstrom), so results can be compared one to one
VDW_RADII = {"H": 1.20, "HE": 1.40, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "NA": 2.27, "MG": 1.73, "P": 1.80,
             "S": 1.80, "CL": 1.75, "K": 2.75, "CA": 2.31, "NI": 1.63, "CU": 1.40, "ZN": 1.39, "SE": 1.90,
             "BR": 1.85, "CD": 1.58, "I": 1.98, "HG": 1.55, "FE": 1.40}
DEFAULT_RADIUS = 1.70

# Theoretical maximum side-chain-plus-backbone ASA of a residue in a Gly-X-Gly tripeptide, Tien et al. 2013 (A^2)
MAX_ASA = {"A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0, "Q": 225.0, "E": 223.0, "G": 104.0,
           "H": 224.0, "I": 197.0, "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0, "S": 155.0,
           "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0}

KYTE_DOOLITTLE = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2,
                  "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
                  "Y": -1.3, "V": 4.2}
# Residues BindCraft counts as hydrophobic for Surface_Hydrophobicity: Rosetta "apolar" types plus F, W, Y
# (Gly is apolar in Rosetta; fitted against final_design_stats.csv, see docs/structure_api.md)
BINDCRAFT_APOLAR = frozenset("AGCILMPVFWY")
PKA = {"Nterm": 8.0, "Cterm": 3.6, "D": 3.9, "E": 4.3, "H": 6.0, "C": 8.3, "Y": 10.1, "K": 10.5, "R": 12.5}


def three_to_one(name: str, strict: bool = False) -> str:
    """Three-letter to one-letter code. Modified residues map to their parent; unknown gives 'X' (or KeyError if strict)."""
    n = name.strip().upper()
    if n in THREE_TO_ONE:
        return THREE_TO_ONE[n]
    if n in _MODIFIED:
        return THREE_TO_ONE[_MODIFIED[n]]
    if strict:
        raise KeyError(name)
    return "X"


def one_to_three(code: str, strict: bool = False) -> str:
    c = code.strip().upper()
    if c in ONE_TO_THREE:
        return ONE_TO_THREE[c]
    if strict:
        raise KeyError(code)
    return "UNK"


# --------------------------------------------------------------------------------------------- tables
class ResidueTable:
    """Per-residue view of an :class:`Atoms` table. All arrays have one entry per residue in file order.

    ``number`` is the author residue number (can have gaps), ``ca`` the atom index of CA (-1 when absent),
    ``start``/``stop`` the half-open atom slice. Look residues up by number with :meth:`index_of`.
    """

    def __init__(self, atoms: "Atoms", start: np.ndarray):
        n_atoms = len(atoms)
        self.start = start.astype(np.int64)
        self.stop = np.append(self.start[1:], n_atoms).astype(np.int64)
        s = self.start
        self.number = atoms.resnum[s].astype(np.int64)
        self.icode = atoms.icode[s]
        self.name = atoms.resname[s]
        self.chain = atoms.chain[s]
        self.one = np.array([three_to_one(n) for n in self.name], dtype="<U1") if len(s) else np.array([], dtype="<U1")
        ca = np.full(len(s), -1, dtype=np.int64)
        is_ca = (atoms.name == "CA") & (atoms.element == "C")
        idx = np.nonzero(is_ca)[0]
        if len(idx):
            ca[atoms.res_index[idx][::-1]] = idx[::-1]  # first CA of a residue wins
        self.ca = ca
        self._atoms = atoms
        self._lookup: dict | None = None

    def __len__(self) -> int:
        return len(self.number)

    @property
    def numbers(self) -> np.ndarray:
        return self.number

    @property
    def sequence(self) -> str:
        return "".join(self.one.tolist())

    @property
    def plddt(self) -> np.ndarray:
        """Per-residue CA B-factor (pLDDT on 0 to 100 for BindCraft output); NaN where a residue has no CA."""
        out = np.full(len(self), np.nan)
        has = self.ca >= 0
        out[has] = self._atoms.bfactor[self.ca[has]]
        return out

    def index_of(self, number: int, chain: str | None = None, icode: str | None = None) -> int:
        """Residue ordinal for a residue number. Raises KeyError if absent. First match wins when numbers repeat."""
        if self._lookup is None:
            lk: dict = {}
            for i in range(len(self) - 1, -1, -1):
                lk[(int(self.number[i]), str(self.chain[i]), str(self.icode[i]))] = i
                lk[(int(self.number[i]), None, None)] = i
            self._lookup = lk
        if chain is None and icode is None:
            return self._lookup[(int(number), None, None)]
        for i in np.nonzero(self.number == number)[0]:
            if (chain is None or self.chain[i] == chain) and (icode is None or self.icode[i] == icode):
                return int(i)
        raise KeyError(number)

    def has(self, number: int) -> bool:
        try:
            self.index_of(number)
            return True
        except KeyError:
            return False

    def indices_of(self, numbers: Iterable[int], missing: str = "skip") -> np.ndarray:
        """Ordinals for several numbers. missing='skip' drops absent numbers, 'raise' raises KeyError."""
        out = []
        for nmb in numbers:
            try:
                out.append(self.index_of(int(nmb)))
            except KeyError:
                if missing == "raise":
                    raise
        return np.array(out, dtype=np.int64)

    def gaps(self) -> list[tuple[int, int]]:
        """(last number before, first number after) for every jump in numbering larger than 1 within one chain."""
        out = []
        for i in range(1, len(self)):
            if self.chain[i] == self.chain[i - 1] and self.number[i] - self.number[i - 1] > 1:
                out.append((int(self.number[i - 1]), int(self.number[i])))
        return out

    def missing_numbers(self, lo: int | None = None, hi: int | None = None) -> list[int]:
        """Residue numbers absent between lo and hi (default: first and last present)."""
        if not len(self):
            return []
        lo = int(self.number.min()) if lo is None else lo
        hi = int(self.number.max()) if hi is None else hi
        have = set(int(x) for x in self.number)
        return [n for n in range(lo, hi + 1) if n not in have]

    def atom_indices(self, i: int) -> np.ndarray:
        return np.arange(self.start[i], self.stop[i])

    def mask_for_numbers(self, numbers: Iterable[int]) -> np.ndarray:
        """Boolean per residue: number in the given set."""
        return np.isin(self.number, np.fromiter((int(x) for x in numbers), dtype=np.int64))


class Atoms:
    """Flat atom table. Atoms of one residue must be contiguous (true for every PDB and mmCIF writer).

    Arrays: ``xyz`` (N, 3) float64, ``name``, ``element`` (upper case), ``resname``, ``resnum`` (int), ``icode``,
    ``chain``, ``bfactor``, ``hetero`` (bool, HETATM). ``res_index`` maps each atom to its residue ordinal and
    ``residues`` is the :class:`ResidueTable`.
    """

    _FIELDS = ("xyz", "name", "element", "resname", "resnum", "icode", "chain", "bfactor", "hetero")

    def __init__(self, xyz, name, element, resname, resnum, chain, bfactor=None, icode=None, hetero=None):
        self.xyz = np.asarray(xyz, dtype=np.float64).reshape(-1, 3)
        n = len(self.xyz)
        self.name = np.asarray(name, dtype="<U4").reshape(n)
        self.element = np.char.upper(np.asarray(element, dtype="<U2")).reshape(n)
        self.resname = np.asarray(resname, dtype="<U3").reshape(n)
        self.resnum = np.asarray(resnum, dtype=np.int64).reshape(n)
        self.chain = np.asarray(chain, dtype="<U4").reshape(n)
        self.bfactor = np.zeros(n) if bfactor is None else np.asarray(bfactor, dtype=np.float64).reshape(n)
        self.icode = np.full(n, "", dtype="<U1") if icode is None else np.asarray(icode, dtype="<U1").reshape(n)
        self.hetero = np.zeros(n, dtype=bool) if hetero is None else np.asarray(hetero, dtype=bool).reshape(n)
        self._residues: ResidueTable | None = None
        self._res_index: np.ndarray | None = None

    # ---- construction / selection
    @classmethod
    def empty(cls) -> "Atoms":
        return cls(np.zeros((0, 3)), [], [], [], [], [])

    @classmethod
    def concat(cls, parts: Sequence["Atoms"]) -> "Atoms":
        parts = [p for p in parts if len(p)]
        if not parts:
            return cls.empty()
        cat = {f: np.concatenate([getattr(p, f) for p in parts]) for f in cls._FIELDS}
        return cls(cat["xyz"], cat["name"], cat["element"], cat["resname"], cat["resnum"], cat["chain"],
                   cat["bfactor"], cat["icode"], cat["hetero"])

    def __len__(self) -> int:
        return len(self.xyz)

    def select(self, mask_or_index) -> "Atoms":
        """New table with the selected atoms (boolean mask or integer indices); original order kept for masks."""
        m = np.asarray(mask_or_index)
        return Atoms(self.xyz[m], self.name[m], self.element[m], self.resname[m], self.resnum[m], self.chain[m],
                     self.bfactor[m], self.icode[m], self.hetero[m])

    def with_xyz(self, xyz) -> "Atoms":
        """Copy with new coordinates (same atoms)."""
        out = Atoms(xyz, self.name, self.element, self.resname, self.resnum, self.chain, self.bfactor, self.icode,
                    self.hetero)
        out._res_index, out._residues = self._res_index, None
        return out

    def transformed(self, R: np.ndarray, t: np.ndarray) -> "Atoms":
        return self.with_xyz(apply_transform(self.xyz, R, t))

    def copy(self) -> "Atoms":
        return self.with_xyz(self.xyz.copy())

    # ---- residues
    @property
    def res_index(self) -> np.ndarray:
        if self._res_index is None:
            n = len(self)
            if n == 0:
                self._res_index = np.zeros(0, dtype=np.int64)
            else:
                chg = np.ones(n, dtype=bool)
                chg[1:] = ((self.chain[1:] != self.chain[:-1]) | (self.resnum[1:] != self.resnum[:-1])
                           | (self.icode[1:] != self.icode[:-1]) | (self.resname[1:] != self.resname[:-1]))
                self._res_index = np.cumsum(chg) - 1
        return self._res_index

    @property
    def residues(self) -> ResidueTable:
        if self._residues is None:
            ri = self.res_index
            start = np.nonzero(np.r_[True, ri[1:] != ri[:-1]])[0] if len(ri) else np.zeros(0, dtype=np.int64)
            self._residues = ResidueTable(self, start)
        return self._residues

    @property
    def n_residues(self) -> int:
        return len(self.residues)

    @property
    def sequence(self) -> str:
        return self.residues.sequence

    @property
    def plddt(self) -> np.ndarray:
        return self.residues.plddt

    @property
    def is_hydrogen(self) -> np.ndarray:
        return np.isin(self.element, ("H", "D"))

    def heavy(self) -> "Atoms":
        return self.select(~self.is_hydrogen)

    def atom(self, resnum: int, name: str) -> int:
        """Atom index of ``name`` in residue number ``resnum``; KeyError when absent."""
        r = self.residues
        i = r.index_of(resnum)
        hit = np.nonzero(self.name[r.start[i]:r.stop[i]] == name)[0]
        if not len(hit):
            raise KeyError((resnum, name))
        return int(r.start[i] + hit[0])

    def ca_xyz(self) -> np.ndarray:
        """(n_res, 3) CA coordinates; NaN rows where a residue lacks CA."""
        r = self.residues
        out = np.full((len(r), 3), np.nan)
        has = r.ca >= 0
        out[has] = self.xyz[r.ca[has]]
        return out

    def ca_by_number(self) -> dict[int, np.ndarray]:
        r = self.residues
        return {int(r.number[i]): self.xyz[r.ca[i]] for i in range(len(r)) if r.ca[i] >= 0}

    def per_residue(self, values: np.ndarray, reduce: str = "sum") -> np.ndarray:
        """Aggregate per-atom values to residues ('sum', 'mean', 'max', 'min')."""
        v = np.asarray(values, dtype=np.float64)
        ri, nres = self.res_index, len(self.residues)
        if reduce == "sum":
            out = np.zeros(nres)
            np.add.at(out, ri, v)
            return out
        if reduce == "mean":
            out = np.zeros(nres)
            np.add.at(out, ri, v)
            return out / np.maximum(np.bincount(ri, minlength=nres), 1)
        if reduce in ("max", "min"):
            out = np.full(nres, -np.inf if reduce == "max" else np.inf)
            (np.maximum if reduce == "max" else np.minimum).at(out, ri, v)
            return out
        raise ValueError(reduce)

    def residue_values_to_atoms(self, values: np.ndarray) -> np.ndarray:
        return np.asarray(values)[self.res_index]

    def select_residues(self, numbers: Iterable[int] | None = None, ordinals: Iterable[int] | None = None) -> "Atoms":
        """Atoms of residues given by number (set) or by ordinal."""
        r = self.residues
        if numbers is not None:
            keep = r.mask_for_numbers(numbers)
        else:
            keep = np.zeros(len(r), dtype=bool)
            keep[np.fromiter((int(i) for i in ordinals), dtype=np.int64)] = True
        return self.select(keep[self.res_index])

    def backbone_arrays(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """N, CA, C, O coordinates per residue (n_res, 3), NaN where the atom is missing."""
        n = len(self.residues)
        out = []
        for nm in ("N", "CA", "C", "O"):
            arr = np.full((n, 3), np.nan)
            idx = np.nonzero((self.name == nm) & ~self.is_hydrogen)[0]
            arr[self.res_index[idx][::-1]] = self.xyz[idx][::-1]
            out.append(arr)
        return tuple(out)  # type: ignore[return-value]

    def __repr__(self) -> str:
        return f"Atoms({len(self)} atoms, {len(self.residues)} residues, chains={sorted(set(self.chain.tolist()))})"


# --------------------------------------------------------------------------------------------- reading
@dataclass
class ValidationIssue:
    code: str
    severity: str  # 'error' | 'warning' | 'info'
    message: str

    def __str__(self) -> str:
        return f"{self.severity}: {self.code}: {self.message}"


@dataclass
class AtomRef:
    chain: str
    resnum: int
    resname: str
    atom: str


@dataclass
class Link:
    """A covalent link or disulfide recorded in the file (mmCIF _struct_conn or PDB LINK/SSBOND)."""
    kind: str  # 'covale', 'disulf', 'metalc', 'hydrog', 'other'
    a: AtomRef
    b: AtomRef
    distance: float | None = None
    name: str = ""


@dataclass
class StructureFile:
    """Everything in one model of a structure file as separate atom groups.

    ``groups``: name -> Atoms. A polymer chain keeps its chain id (e.g. 'A'); a branched/glycan chain keeps its id
    with kind 'branched', or is 'A:branched' when its residues carry the chain id of a polymer (newer entries number the sugars of
    chain A as part of chain A); non-polymer ligands that share a polymer chain are 'A:ligands' and waters 'A:waters'.
    ``links`` carries covalent links (ASN-NAG, disulfides) with atom references; use :meth:`find_atom` to resolve.
    """
    path: str
    groups: dict[str, Atoms]
    kinds: dict[str, str]
    links: list[Link] = field(default_factory=list)
    resolution: float | None = None
    entry_id: str = ""

    def polymer(self, chain: str) -> Atoms:
        return self.groups[chain]

    def polymer_chains(self) -> list[str]:
        return [g for g, k in self.kinds.items() if k == "polymer"]

    def glycan_groups(self) -> list[str]:
        return [g for g, k in self.kinds.items() if k == "branched"]

    def find_atom(self, chain: str, resnum: int, atom: str) -> tuple[str, int] | None:
        """(group name, atom index) of the atom, searching every group that carries ``chain``; None if absent."""
        for gname, at in self.groups.items():
            if gname.split(":")[0] != chain and gname != chain:
                continue
            hit = np.nonzero((at.chain == chain) & (at.resnum == resnum) & (at.name == atom))[0]
            if len(hit):
                return gname, int(hit[0])
        return None


def _gemmi():
    import gemmi
    return gemmi


def _read_gemmi(path: str):
    gemmi = _gemmi()
    try:
        st = gemmi.read_structure(str(path))
    except Exception as ex:  # gemmi raises RuntimeError / ValueError with parse messages
        raise StructureError(f"cannot read {path}: {ex}") from ex
    if not len(st):
        raise StructureError(f"{path} has no models")
    st.setup_entities()
    return st


def _links(st) -> list[Link]:
    """Covalent links the file records: mmCIF _struct_conn, PDB LINK and SSBOND (gemmi reads both into connections)."""
    gemmi = _gemmi()
    kinds = {gemmi.ConnectionType.Covale: "covale", gemmi.ConnectionType.Disulf: "disulf",
             gemmi.ConnectionType.MetalC: "metalc", gemmi.ConnectionType.Hydrog: "hydrog"}
    out = []
    for c in st.connections:
        p1, p2 = c.partner1, c.partner2
        out.append(Link(kinds.get(c.type, "other"), AtomRef(p1.chain_name, p1.res_id.seqid.num, p1.res_id.name, p1.atom_name),
                        AtomRef(p2.chain_name, p2.res_id.seqid.num, p2.res_id.name, p2.atom_name),
                        float(c.reported_distance) if c.reported_distance else None, c.name))
    return out


def _residue_atoms(res, chain_name: str, rows: list) -> None:
    """Append atom rows of one residue (first altloc only) to ``rows``."""
    first_alt = None
    resname = res.name
    modified = resname in _MODIFIED
    parent = _MODIFIED.get(resname, resname)
    het = res.het_flag == "H"
    for a in res:
        if a.altloc != "\0":
            if first_alt is None:
                first_alt = a.altloc
            if a.altloc != first_alt:
                continue
        el = a.element.name.upper()
        nm = a.name
        if resname == "MSE" and nm == "SE":
            nm, el = "SD", "S"
        rows.append((a.pos.x, a.pos.y, a.pos.z, nm, el, parent if modified else resname, res.seqid.num,
                     (res.seqid.icode if res.seqid.icode.strip() else ""), chain_name, a.b_iso,
                     het and not modified))


def _atoms_from_rows(rows: list) -> Atoms:
    if not rows:
        return Atoms.empty()
    x, y, z, nm, el, rn, num, ic, ch, b, het = zip(*rows)
    return Atoms(np.column_stack([x, y, z]), nm, el, rn, num, ch, b, ic, het)


def read_structure(path: str, model: int = 0, include_waters: bool = False) -> StructureFile:
    """Read any chain or group of a PDB or mmCIF file into :class:`Atoms` groups (first altloc, hydrogens kept).

    Works for BindCraft PDBs, the truncated target PDB, BindCraft 2 CIFs and PDB entries of oligomers (chains A and B as polymer
    groups, sugar chains as 'branched' groups with their covalent links).
    """
    gemmi = _gemmi()
    st = _read_gemmi(path)
    mdl = st[model]
    groups: dict[str, Atoms] = {}
    kinds: dict[str, str] = {}
    for ch in mdl:
        poly_rows, branch_rows, lig_rows, wat_rows = [], [], [], []
        for res in ch:
            ent = st.get_entity(res.entity_id) if res.entity_id else None
            etype = ent.entity_type if ent is not None else None
            if res.name in _WATER or (etype is not None and etype == gemmi.EntityType.Water):
                _residue_atoms(res, ch.name, wat_rows)
            elif etype is not None and etype == gemmi.EntityType.Branched:
                _residue_atoms(res, ch.name, branch_rows)
            elif etype is not None and etype == gemmi.EntityType.Polymer:
                _residue_atoms(res, ch.name, poly_rows)
            elif etype is None and res.het_flag != "H":
                _residue_atoms(res, ch.name, poly_rows)  # entities missing in a PDB file: ATOM means polymer
            elif res.name in THREE_TO_ONE or res.name in _MODIFIED and res.het_flag != "H":
                _residue_atoms(res, ch.name, poly_rows)
            else:
                _residue_atoms(res, ch.name, lig_rows)
        if poly_rows:
            groups[ch.name] = _atoms_from_rows(poly_rows)
            kinds[ch.name] = "polymer"
        if branch_rows:                                   # the sugars of a polymer chain are not part of its residues
            nm = f"{ch.name}:branched" if poly_rows else ch.name
            groups[nm] = _atoms_from_rows(branch_rows)
            kinds[nm] = "branched"
        if lig_rows:
            nm = f"{ch.name}:ligands" if (poly_rows or branch_rows) else ch.name
            groups[nm] = _atoms_from_rows(lig_rows)
            kinds[nm] = "ligand"
        if wat_rows and include_waters:
            groups[f"{ch.name}:waters"] = _atoms_from_rows(wat_rows)
            kinds[f"{ch.name}:waters"] = "water"
    links = _links(st)
    res = float(st.resolution) if st.resolution else None
    return StructureFile(str(path), groups, kinds, links, res, st.name or "")


def deposited_ss(path: str, model: int = 0) -> dict[str, dict[int, str]]:
    """Deposited secondary structure (HELIX / SHEET records or mmCIF struct_conf / struct_sheet_range).

    Returns chain -> {residue number: 'H' or 'E'}; used to test :func:`assign_ss`. Residues not listed are loops.
    """
    st = _read_gemmi(path)
    out: dict[str, dict[int, str]] = {}
    for h in st.helices:
        c = h.start.chain_name
        for k in range(h.start.res_id.seqid.num, h.end.res_id.seqid.num + 1):
            out.setdefault(c, {})[k] = "H"
    for sh in st.sheets:
        for strand in sh.strands:
            c = strand.start.chain_name
            for k in range(strand.start.res_id.seqid.num, strand.end.res_id.seqid.num + 1):
                out.setdefault(c, {})[k] = "E"
    return out


@dataclass
class Complex:
    """A target and a binder from one file.

    ``atoms`` is every kept atom (target, binder and other kept chains, in file order); ``target`` and ``binder``
    are separate :class:`Atoms` tables. ``hydrogens`` holds every hydrogen of the target and binder chains found
    in the file (also when ``strip_h=True`` removed them from the tables), so BindCraft-compatible interface
    residues can count them. Residue-level helpers key on ``.residues.number``.
    """
    path: str
    atoms: Atoms
    target: Atoms
    binder: Atoms
    hydrogens: Atoms
    target_chain: str
    binder_chain: str                                   # the first binder chain; `binder_chains` has them all
    issues: list[ValidationIssue]
    other: Atoms = field(default_factory=Atoms.empty)
    hydrogens_stripped: bool = True
    links: list[Link] = field(default_factory=list)      # covalent links the file records (often none: BindCraft 2 writes no _struct_conn)
    binder_chains: tuple[str, ...] = ()                  # every chain of the binder, in file order

    def recorded_disulfides(self, chain: str | None = None) -> list[tuple[int, int]]:
        """Disulfides the file records inside one chain (default: the binder) as pairs of file residue numbers."""
        chain = self.binder_chain if chain is None else chain
        return sorted({tuple(sorted((ln.a.resnum, ln.b.resnum))) for ln in self.links
                       if ln.kind == "disulf" and ln.a.chain == chain == ln.b.chain and ln.a.resnum != ln.b.resnum})

    @property
    def has_hydrogens(self) -> bool:
        return len(self.hydrogens) > 0

    @property
    def errors(self) -> list[ValidationIssue]:
        return [i for i in self.issues if i.severity == "error"]

    @property
    def target_sequence(self) -> str:
        return self.target.sequence

    @property
    def binder_sequence(self) -> str:
        return self.binder.sequence

    @property
    def target_numbers(self) -> np.ndarray:
        return self.target.residues.number

    @property
    def binder_numbers(self) -> np.ndarray:
        return self.binder.residues.number

    @property
    def target_plddt(self) -> np.ndarray:
        return self.target.plddt

    @property
    def binder_plddt(self) -> np.ndarray:
        return self.binder.plddt

    def combined(self) -> Atoms:
        """Target atoms then binder atoms (the order used by interface helpers)."""
        return Atoms.concat([self.target, self.binder])


def load_complex(path: str, target_chain: str = "A", binder_chain: str | Sequence[str] | None = "B",
                 strip_h: bool = True, strict: bool = True, keep_hetero: bool = False, model: int = 0,
                 n_binder_chains: int | None = None) -> Complex:
    """Load a complex of a target and a binder from PDB or mmCIF.

    BindCraft 1 relaxed PDBs (explicit hydrogens, TER records), the truncated target PDB (pass ``binder_chain=None``
    for a single chain) and mmCIF all work. Only the first altloc is kept. Non-protein HETATM records and waters
    are dropped unless ``keep_hetero``. With ``strict`` a missing chain or an empty file raises
    :class:`StructureError`; otherwise the problem is only recorded in ``issues``.

    **The binder may be several chains** -- an scFv's VH and VL, a homo-oligomer's
    copies, a multidomain design written as domains. Name them (``binder_chain=["B", "C"]``),
    or give ``n_binder_chains`` and the last that many chains of the file are taken,
    which is where BindCraft 2 writes the binder: chains are lettered by position with
    the targets first. ``binder`` is then one table over all of them, in file order,
    and ``binder_chains`` names them; ``binder_chain`` stays the first, for callers
    that handle one.
    """
    st = _read_gemmi(path)
    mdl = st[model]
    rows_by_chain: dict[str, list] = {}
    dropped_het: dict[str, int] = {}
    for ch in mdl:
        rows: list = []
        for res in ch:
            is_water = res.name in _WATER
            foreign = res.het_flag == "H" and res.name not in THREE_TO_ONE and res.name not in _MODIFIED
            if is_water or (foreign and not keep_hetero):
                dropped_het[ch.name] = dropped_het.get(ch.name, 0) + 1
                continue
            _residue_atoms(res, ch.name, rows)
        if rows:
            rows_by_chain[ch.name] = rows
    issues: list[ValidationIssue] = []
    order = list(rows_by_chain)                       # file order, which is what "the last n chains" means
    if n_binder_chains is not None:
        n = max(0, int(n_binder_chains))
        rest = [c for c in order if c != target_chain]
        binder_names = rest[-n:] if n else []
        if n and len(binder_names) < n:
            issues.append(ValidationIssue(
                "missing_chain", "error",
                f"{n} binder chains asked for and {len(binder_names)} available besides the target "
                f"(chains present: {order})"))
    elif binder_chain is None:
        binder_names = []
    elif isinstance(binder_chain, str):
        binder_names = [binder_chain]
    else:
        binder_names = list(binder_chain)
    want = [c for c in ([target_chain] + binder_names) if c is not None]
    for c in want:
        if c not in rows_by_chain:
            issues.append(ValidationIssue("missing_chain", "error",
                                          f"chain {c!r} not found (chains present: {sorted(rows_by_chain)})"))
    if strict and any(i.severity == "error" for i in issues):
        raise StructureError(f"{path}: " + "; ".join(i.message for i in issues if i.severity == "error"))

    def table(c: str | None) -> Atoms:
        return _atoms_from_rows(rows_by_chain[c]) if c is not None and c in rows_by_chain else Atoms.empty()

    tgt_all = table(target_chain)
    # One table over every binder chain, in file order: downstream this is "the
    # binder", and its atoms carry their own chain ids, so a consumer that needs
    # to tell the chains apart still can.
    bnd_all = Atoms.concat([table(c) for c in binder_names]) if binder_names else Atoms.empty()
    other_chains = [c for c in rows_by_chain if c not in want]
    other = Atoms.concat([table(c) for c in other_chains])
    if other_chains:
        issues.append(ValidationIssue("extra_chains", "info", f"chains not used as target or binder: {other_chains}"))
    hyd = Atoms.concat([tgt_all.select(tgt_all.is_hydrogen), bnd_all.select(bnd_all.is_hydrogen)])
    if len(hyd):
        issues.append(ValidationIssue("hydrogens", "info",
                                      f"{len(hyd)} hydrogen atoms in the file ({'stripped' if strip_h else 'kept'})"))
    if strip_h:
        tgt_all, bnd_all = tgt_all.heavy(), bnd_all.heavy()
        other = other.heavy() if len(other) else other
    if dropped_het:
        issues.append(ValidationIssue("hetatm", "warning",
                                      f"dropped {sum(dropped_het.values())} non-protein or water residues "
                                      f"(by chain: {dropped_het})"))
    for label, tab in (("target", tgt_all), ("binder", bnd_all)):
        if len(tab):
            bb = tab.backbone_arrays()
            bad = int(np.isnan(np.column_stack([b[:, 0] for b in bb])).any(axis=1).sum())
            if bad:
                issues.append(ValidationIssue("incomplete_backbone", "warning",
                                              f"{bad} {label} residues lack N, CA, C or O"))
            gaps = tab.residues.gaps()
            if gaps:
                issues.append(ValidationIssue("numbering_gaps", "info",
                                              f"{label} residue numbering has {len(gaps)} gap(s): {gaps[:6]}"))
            if len(tab) and (np.nanmax(tab.bfactor) > 100.0 or np.nanmin(tab.bfactor) < 0.0):
                issues.append(ValidationIssue("bfactor_range", "info",
                                              f"{label} B-factors outside 0 to 100, so not pLDDT on a 0 to 100 scale"))
    if strict and not len(tgt_all) and target_chain in want:
        raise StructureError(f"{path}: no atoms in target chain {target_chain!r}")
    atoms = Atoms.concat([tgt_all, bnd_all, other])
    return Complex(str(path), atoms, tgt_all, bnd_all, hyd, target_chain, binder_names[0] if binder_names else "",
                   issues, other, strip_h, _links(st), tuple(binder_names))


# --------------------------------------------------------------------------------------------- writing
def to_pdb_lines(atoms: Atoms, serial0: int = 1, chain: str | None = None, resnum_offset: int = 0,
                 bfactor: np.ndarray | None = None) -> list[str]:
    """PDB ATOM/HETATM lines (no TER). ``chain`` overrides the chain id, ``bfactor`` overrides per-atom B."""
    lines = []
    for i in range(len(atoms)):
        nm = atoms.name[i]
        nm_fmt = (" " + nm).ljust(4) if len(nm) < 4 and len(atoms.element[i]) < 2 else nm.ljust(4)
        b = atoms.bfactor[i] if bfactor is None else bfactor[i]
        rec = "HETATM" if atoms.hetero[i] else "ATOM  "
        lines.append("%s%5d %4s %3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s" % (
            rec, (serial0 + i) % 100000, nm_fmt, atoms.resname[i], (chain or atoms.chain[i])[:1],
            (int(atoms.resnum[i]) + resnum_offset) % 10000, atoms.icode[i] or " ", atoms.xyz[i, 0], atoms.xyz[i, 1],
            atoms.xyz[i, 2], 1.0, b, atoms.element[i].rjust(2)))
    return lines


def write_pdb(path: str, *groups: Atoms) -> None:
    """Write one or more atom groups (each followed by TER) to a PDB file."""
    out, serial = [], 1
    for g in groups:
        out += to_pdb_lines(g, serial)
        serial += len(g)
        out.append("TER")
    out.append("END")
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")


# --------------------------------------------------------------------------------------------- geometry
def kabsch(P: np.ndarray, Q: np.ndarray, weights: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, float]:
    """Least-squares rigid superposition of P onto Q.

    Returns (R, t, rmsd) such that ``P @ R.T + t`` best matches Q (proper rotation, determinant +1).
    """
    P = np.asarray(P, dtype=np.float64)
    Q = np.asarray(Q, dtype=np.float64)
    if P.shape != Q.shape or P.ndim != 2 or P.shape[1] != 3:
        raise ValueError(f"P and Q must both be (n, 3); got {P.shape} and {Q.shape}")
    if len(P) < 3:
        raise ValueError("kabsch needs at least 3 point pairs")
    if weights is None:
        w = np.full(len(P), 1.0 / len(P))
    else:
        w = np.asarray(weights, dtype=np.float64)
        w = w / w.sum()
    Pc, Qc = (w[:, None] * P).sum(0), (w[:, None] * Q).sum(0)
    H = (P - Pc).T @ (w[:, None] * (Q - Qc))
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d if d != 0 else 1.0]) @ U.T
    t = Qc - R @ Pc
    dev = P @ R.T + t - Q
    return R, t, float(np.sqrt((w * (dev ** 2).sum(1)).sum()))


def apply_transform(xyz: np.ndarray, R: np.ndarray, t: np.ndarray) -> np.ndarray:
    return np.asarray(xyz, dtype=np.float64) @ np.asarray(R).T + np.asarray(t)


def invert_transform(R: np.ndarray, t: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return R.T, -R.T @ t


def rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """Plain RMSD between matched point sets (no superposition)."""
    P, Q = np.asarray(P, dtype=np.float64), np.asarray(Q, dtype=np.float64)
    return float(np.sqrt(((P - Q) ** 2).sum(1).mean()))


@dataclass
class Superposition:
    R: np.ndarray
    t: np.ndarray
    rmsd: float
    n: int
    numbers: list[int]  # reference residue numbers used

    def apply(self, obj):
        """Transform an (n, 3) array or an :class:`Atoms` table."""
        if isinstance(obj, Atoms):
            return obj.transformed(self.R, self.t)
        return apply_transform(obj, self.R, self.t)


def superpose_residues(mobile: Atoms, ref: Atoms, atom: str = "CA", numbers: Iterable[int] | None = None,
                       mapping: dict[int, int] | None = None) -> Superposition:
    """Superpose ``mobile`` onto ``ref`` over residues present in both, matched by residue number.

    ``numbers`` restricts the reference residue numbers used; ``mapping`` maps mobile numbers to reference
    numbers when the numbering differs. ``atom`` is one atom name ('CA') or 'backbone' (N, CA, C, O).
    """
    names = ("N", "CA", "C", "O") if atom == "backbone" else (atom,)
    def coords(a: Atoms) -> dict[tuple[int, str], np.ndarray]:
        r = a.residues
        out = {}
        for nm in names:
            idx = np.nonzero((a.name == nm) & ~a.is_hydrogen)[0]
            for i in idx[::-1]:
                out[(int(r.number[a.res_index[i]]), nm)] = a.xyz[i]
        return out
    cm, cr = coords(mobile), coords(ref)
    allowed = None if numbers is None else set(int(x) for x in numbers)
    P, Q, used = [], [], []
    mobile_nums = sorted({k[0] for k in cm})
    for mn in mobile_nums:
        rn = mapping.get(mn) if mapping is not None else mn
        if rn is None or (allowed is not None and rn not in allowed):
            continue
        if all((mn, nm) in cm and (rn, nm) in cr for nm in names):
            for nm in names:
                P.append(cm[(mn, nm)])
                Q.append(cr[(rn, nm)])
            used.append(rn)
    if len(P) < 3:
        raise StructureError(f"only {len(P)} shared atoms to superpose")
    R, t, r = kabsch(np.array(P), np.array(Q))
    return Superposition(R, t, r, len(used), used)


def rmsd_by_number(a: Atoms, b: Atoms, atom: str = "CA", superpose: bool = True,
                   numbers: Iterable[int] | None = None) -> float:
    """CA (or other single atom) RMSD over residues with the same number in both tables."""
    ca, cb = a.ca_by_number() if atom == "CA" else None, b.ca_by_number() if atom == "CA" else None
    if atom != "CA":
        def one(t: Atoms):
            r = t.residues
            return {int(r.number[t.res_index[i]]): t.xyz[i] for i in np.nonzero(t.name == atom)[0][::-1]}
        ca, cb = one(a), one(b)
    common = sorted(set(ca) & set(cb) if numbers is None else (set(ca) & set(cb) & set(int(x) for x in numbers)))
    if len(common) < 3:
        raise StructureError("fewer than 3 shared residues")
    P, Q = np.array([ca[k] for k in common]), np.array([cb[k] for k in common])
    if superpose:
        return kabsch(P, Q)[2]
    return rmsd(P, Q)


# --------------------------------------------------------------------------------------------- contacts
@dataclass
class Contacts:
    """Residue-residue contacts between two atom tables.

    ``ia``/``ib`` are residue ordinals in ``a.residues`` / ``b.residues``; ``num_a``/``num_b`` the residue numbers.
    ``min_dist`` is the smallest atom-atom distance of the pair, ``n_atom_pairs`` the count of atom pairs within the
    cutoff. ``a_residues`` / ``b_residues`` are sorted residue numbers with at least one contact.
    """
    cutoff: float
    ia: np.ndarray
    ib: np.ndarray
    num_a: np.ndarray
    num_b: np.ndarray
    min_dist: np.ndarray
    n_atom_pairs: np.ndarray

    def __len__(self) -> int:
        return len(self.ia)

    @property
    def a_residues(self) -> list[int]:
        return sorted(set(int(x) for x in self.num_a))

    @property
    def b_residues(self) -> list[int]:
        return sorted(set(int(x) for x in self.num_b))

    @property
    def n_atom_contacts(self) -> int:
        return int(self.n_atom_pairs.sum())

    def pairs(self) -> list[tuple[int, int]]:
        return [(int(x), int(y)) for x, y in zip(self.num_a, self.num_b)]

    def as_records(self) -> list[dict]:
        return [dict(a=int(x), b=int(y), d=round(float(d), 2), n=int(n))
                for x, y, d, n in zip(self.num_a, self.num_b, self.min_dist, self.n_atom_pairs)]


def atom_contacts(a: Atoms, b: Atoms, cutoff: float = 4.0, heavy_only: bool = True) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Atom pairs within ``cutoff``: (index into a, index into b, distance). With heavy_only hydrogens are ignored."""
    ia_all = np.nonzero(~a.is_hydrogen)[0] if heavy_only else np.arange(len(a))
    ib_all = np.nonzero(~b.is_hydrogen)[0] if heavy_only else np.arange(len(b))
    if not len(ia_all) or not len(ib_all):
        z = np.zeros(0, dtype=np.int64)
        return z, z, np.zeros(0)
    ta, tb = cKDTree(a.xyz[ia_all]), cKDTree(b.xyz[ib_all])
    pairs = ta.sparse_distance_matrix(tb, cutoff, output_type="coo_matrix")
    return ia_all[pairs.row], ib_all[pairs.col], pairs.data


def residue_contacts(a: Atoms, b: Atoms, cutoff: float = 4.0, heavy_only: bool = True) -> Contacts:
    """Residue pairs with any atom pair closer than ``cutoff`` (default 4.0 A, heavy atoms)."""
    ai, bi, d = atom_contacts(a, b, cutoff, heavy_only)
    if not len(ai):
        z = np.zeros(0, dtype=np.int64)
        return Contacts(cutoff, z, z, z, z, np.zeros(0), z)
    ra, rb = a.res_index[ai], b.res_index[bi]
    key = ra.astype(np.int64) * (len(b.residues) + 1) + rb
    order = np.argsort(key, kind="stable")
    key, d, ra, rb = key[order], d[order], ra[order], rb[order]
    first = np.r_[True, key[1:] != key[:-1]]
    starts = np.nonzero(first)[0]
    min_d = np.minimum.reduceat(d, starts)
    cnt = np.diff(np.r_[starts, len(key)])
    ia, ib = ra[starts], rb[starts]
    return Contacts(cutoff, ia, ib, a.residues.number[ia], b.residues.number[ib], min_d, cnt)


def bindcraft_interface_residues(cx: Complex, cutoff: float = 4.0, include_hydrogens: bool = True) -> list[int]:
    """Binder residue numbers BindCraft 1 reports as InterfaceResidues.

    BindCraft's hotspot_residues(): a binder residue is interface if any of its atoms lies within ``cutoff`` (4.0 A)
    of any target atom. BindCraft runs this on the relaxed PDB WITH explicit hydrogens, so
    ``include_hydrogens=True`` (default) reproduces its column; False gives a heavy-atom-only interface.
    The CSV row belongs to the last validated model (model2 for default multimer settings), so the model1 PDB of
    a design can differ slightly from the CSV.
    """
    tgt, bnd = cx.target, cx.binder
    if include_hydrogens and cx.has_hydrogens and cx.hydrogens_stripped:
        th = cx.hydrogens.select(cx.hydrogens.chain == cx.target_chain)
        bh = cx.hydrogens.select(cx.hydrogens.chain == cx.binder_chain)
        txyz = np.vstack([tgt.xyz, th.xyz])
        bxyz = np.vstack([bnd.xyz, bh.xyz])
        bres = np.concatenate([bnd.res_index, _hydrogen_residue_index(bnd, bh)])
    elif include_hydrogens:
        txyz, bxyz, bres = tgt.xyz, bnd.xyz, bnd.res_index
    else:
        tm, bm = ~tgt.is_hydrogen, ~bnd.is_hydrogen
        txyz, bxyz, bres = tgt.xyz[tm], bnd.xyz[bm], bnd.res_index[bm]
    if not len(txyz) or not len(bxyz):
        return []
    hit = cKDTree(txyz).query_ball_point(bxyz, cutoff, return_length=True) > 0
    return sorted(int(x) for x in np.unique(bnd.residues.number[bres[hit]]))


def _hydrogen_residue_index(heavy: Atoms, hyd: Atoms) -> np.ndarray:
    """Residue ordinal in ``heavy`` for each hydrogen, matched by residue number (-1 -> clipped to 0 if unmatched)."""
    r = heavy.residues
    lookup = {(int(n), str(ic)): i for i, (n, ic) in enumerate(zip(r.number, r.icode))}
    return np.array([lookup.get((int(n), str(ic)), 0) for n, ic in zip(hyd.resnum, hyd.icode)], dtype=np.int64)


def format_bindcraft_interface(numbers: Iterable[int], chain: str = "B") -> str:
    """'B1,B2,B5' format of the InterfaceResidues column."""
    return ",".join(f"{chain}{int(n)}" for n in sorted(numbers))


# --------------------------------------------------------------------------------------------- SASA
SASA_PRESETS = {"fast": 48, "default": 200, "accurate": 960}
_SPHERE_CACHE: dict[int, np.ndarray] = {}


def _unit_sphere(n: int) -> np.ndarray:
    """Golden spiral points on the unit sphere (same construction as Bio.PDB.SASA)."""
    if n not in _SPHERE_CACHE:
        dl = math.pi * (3.0 - math.sqrt(5.0))
        dz = 2.0 / n
        k = np.arange(n)
        z = 1.0 - dz / 2.0 - dz * k
        r = np.sqrt(np.maximum(0.0, 1.0 - z * z))
        lon = dl * k
        _SPHERE_CACHE[n] = np.column_stack([np.cos(lon) * r, np.sin(lon) * r, z])
    return _SPHERE_CACHE[n]


def _radii_for(element: np.ndarray, radii: np.ndarray | dict | None) -> np.ndarray:
    if radii is None:
        return np.array([VDW_RADII.get(e, DEFAULT_RADIUS) for e in element], dtype=np.float64)
    if isinstance(radii, dict):
        return np.array([radii.get(e, VDW_RADII.get(e, DEFAULT_RADIUS)) for e in element], dtype=np.float64)
    return np.asarray(radii, dtype=np.float64)


def sasa_atoms(xyz: np.ndarray, element: Sequence[str] | np.ndarray, probe: float = 1.4,
               n_points: int | str = "default", radii: np.ndarray | dict | None = None) -> np.ndarray:
    """Shrake-Rupley solvent accessible surface area per atom (A^2), vectorised.

    ``n_points`` is the number of test points per atom, or a preset name ('fast' 48, 'default' 200, 'accurate' 960).
    The golden-spiral sphere and the element radii are the same as Bio.PDB.SASA.ShrakeRupley, so equal settings
    give the same numbers. Hydrogens count if they are in the arrays: strip them first for heavy-atom values.
    """
    xyz = np.ascontiguousarray(xyz, dtype=np.float64)
    n = len(xyz)
    if n == 0:
        return np.zeros(0)
    npts = SASA_PRESETS[n_points] if isinstance(n_points, str) else int(n_points)
    rad = _radii_for(np.char.upper(np.asarray(element, dtype="<U2")), radii) + probe
    unit = _unit_sphere(npts)
    tree = cKDTree(xyz)
    pairs = tree.query_pairs(2.0 * rad.max(), output_type="ndarray")
    if len(pairs):
        dd = np.linalg.norm(xyz[pairs[:, 0]] - xyz[pairs[:, 1]], axis=1)
        pairs = pairs[dd < rad[pairs[:, 0]] + rad[pairs[:, 1]]]
    src = np.concatenate([pairs[:, 0], pairs[:, 1]]) if len(pairs) else np.zeros(0, dtype=np.int64)
    dst = np.concatenate([pairs[:, 1], pairs[:, 0]]) if len(pairs) else np.zeros(0, dtype=np.int64)
    order = np.argsort(src, kind="stable")
    src, dst = src[order], dst[order]
    counts = np.bincount(src, minlength=n)
    offs = np.r_[0, np.cumsum(counts)]
    accessible = np.full(n, float(npts))
    by_count = np.argsort(counts, kind="stable")
    by_count = by_count[counts[by_count] > 0]
    pos = 0
    budget = 3_000_000   # elements of the (atoms, points, neighbours) test array per chunk
    while pos < len(by_count):
        # counts ascend, so the chunk maximum is its last element; shrink the chunk until it fits the budget
        end = len(by_count)
        while end - pos > 1 and (end - pos) * npts * counts[by_count[end - 1]] > budget:
            end = pos + max(1, (end - pos) // 2)
        idx = by_count[pos:end]
        kmax = int(counts[idx].max())
        nb = np.full((len(idx), kmax), -1, dtype=np.int64)
        for row, a in enumerate(idx):
            nb[row, :counts[a]] = dst[offs[a]:offs[a + 1]]
        valid = nb >= 0
        nbc = np.where(valid, nb, 0)
        rel = xyz[nbc] - xyz[idx][:, None, :]                      # (c, k, 3)
        rel[~valid] = 1.0e6
        nbr2 = np.where(valid, (rad[nbc] ** 2), 0.0)               # buried if |p - x_j|^2 < r_j^2
        pts = unit[None, :, :] * rad[idx][:, None, None]            # (c, P, 3), relative to atom i
        d2 = (rad[idx] ** 2)[:, None, None] + (rel ** 2).sum(2)[:, None, :] - 2.0 * np.matmul(pts, rel.transpose(0, 2, 1))
        buried = (d2 < nbr2[:, None, :]).any(axis=2)               # (c, P)
        accessible[idx] = npts - buried.sum(1)
        pos = end
    return 4.0 * math.pi * rad ** 2 * (accessible / npts)


@dataclass
class SasaResult:
    atom: np.ndarray      # per atom, A^2
    residue: np.ndarray   # per residue, A^2 (sum of atoms)
    relative: np.ndarray  # per residue, residue SASA / MAX_ASA (NaN for nonstandard residues; can exceed 1 slightly)
    total: float


def relative_sasa(residue_sasa: np.ndarray, one_letter: Sequence[str] | np.ndarray) -> np.ndarray:
    mx = np.array([MAX_ASA.get(str(c), np.nan) for c in one_letter], dtype=np.float64)
    return np.asarray(residue_sasa, dtype=np.float64) / mx


def sasa(atoms: Atoms, probe: float = 1.4, n_points: int | str = "default", radii=None,
         heavy_only: bool = True) -> SasaResult:
    """SASA of a whole table: per atom, per residue, relative per residue and total.

    Hydrogens are ignored when ``heavy_only`` (default; the MAX_ASA table is a heavy-atom quantity). When
    hydrogens are ignored the per-atom array still has one entry per atom (0 for hydrogens).
    """
    if heavy_only and atoms.is_hydrogen.any():
        keep = ~atoms.is_hydrogen
        at = np.zeros(len(atoms))
        at[keep] = sasa_atoms(atoms.xyz[keep], atoms.element[keep], probe, n_points, _subset_radii(radii, keep))
    else:
        at = sasa_atoms(atoms.xyz, atoms.element, probe, n_points, radii)
    res = atoms.per_residue(at)
    return SasaResult(at, res, relative_sasa(res, atoms.residues.one), float(at.sum()))


def _subset_radii(radii, keep):
    if radii is None or isinstance(radii, dict):
        return radii
    return np.asarray(radii)[keep]


@dataclass
class InterfaceSasa:
    """Buried surface of a binder-target interface (A^2). 'binder' and 'target' are the changes on each side."""
    binder_atom: np.ndarray       # alone minus in complex, per binder atom
    target_atom: np.ndarray
    binder_residue: np.ndarray    # per binder residue
    target_residue: np.ndarray
    dsasa_binder: float
    dsasa_target: float
    binder_sasa_complex: np.ndarray    # per binder atom, in complex
    binder_sasa_alone: np.ndarray
    target_sasa_complex: np.ndarray
    target_sasa_alone: np.ndarray

    @property
    def dsasa_total(self) -> float:
        """Buried area counted on both sides (Rosetta's interface dSASA convention)."""
        return self.dsasa_binder + self.dsasa_target

    @property
    def binder_buried_fraction(self) -> float:
        a = self.binder_sasa_alone.sum()
        return float(self.dsasa_binder / a) if a else 0.0


def interface_dsasa(cx: Complex, probe: float = 1.4, n_points: int | str = "default",
                    target_alone: np.ndarray | None = None, radii=None) -> InterfaceSasa:
    """dSASA of binder and target: SASA alone minus SASA in the complex (heavy atoms only).

    ``target_alone`` (per target heavy atom) can be passed to avoid recomputing the target when many designs share
    one target. Atom arrays are indexed on the heavy-atom subset of ``cx.target`` / ``cx.binder``.
    """
    tgt, bnd = cx.target.heavy(), cx.binder.heavy()
    xyz = np.vstack([tgt.xyz, bnd.xyz])
    el = np.concatenate([tgt.element, bnd.element])
    sc = sasa_atoms(xyz, el, probe, n_points, radii)
    t_c, b_c = sc[:len(tgt)], sc[len(tgt):]
    b_a = sasa_atoms(bnd.xyz, bnd.element, probe, n_points, radii)
    t_a = sasa_atoms(tgt.xyz, tgt.element, probe, n_points, radii) if target_alone is None else np.asarray(target_alone)
    db, dt = b_a - b_c, t_a - t_c
    return InterfaceSasa(db, dt, bnd.per_residue(db), tgt.per_residue(dt), float(db.sum()), float(dt.sum()),
                         b_c, b_a, t_c, t_a)


# --------------------------------------------------------------------------------------------- secondary structure
_HB_Q = 0.084 * 332.0


def _hbond_energy_matrix(N, CA, C, O, cont) -> np.ndarray:
    """Kabsch-Sander energies E[i, j] (kcal/mol) for C=O of i accepting from N-H of j (NaN where undefined)."""
    n = len(CA)
    H = np.full((n, 3), np.nan)
    for_h = np.zeros(n, dtype=bool)
    for_h[1:] = cont[:-1]
    v = C[:-1] - O[:-1]
    nrm = np.linalg.norm(v, axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        H[1:] = N[1:] + v / nrm
    H[~for_h] = np.nan
    with np.errstate(invalid="ignore", divide="ignore"):
        def d(A, B):
            return np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)
        rON, rCH, rOH, rCN = d(O, N), d(C, H), d(O, H), d(C, N)
        E = _HB_Q * (1.0 / rON + 1.0 / rCH - 1.0 / rOH - 1.0 / rCN)
    close = d(CA, CA) < 12.0
    E[~close] = np.nan
    idx = np.arange(n)
    adj = np.abs(idx[:, None] - idx[None, :]) < 2   # same residue and the peptide bond itself
    E[adj] = np.nan
    return E


def assign_dssp(atoms: Atoms, hb_cutoff: float = -0.5) -> str:
    """Eight-state DSSP-like assignment from backbone N, CA, C, O only (no hydrogens needed).

    Kabsch-Sander hydrogen-bond energies with the amide H placed from the preceding C=O; returns one character per
    residue of ``atoms`` in table order: H alpha helix, G 3-10, I pi, E strand (ladder of two or more bridges),
    B isolated bridge, T turn, S bend, '-' none (also residues without a full backbone). Chain breaks are found from
    C(i)-N(i+1) distance, not from numbering, so numbering gaps (deleted loops) are safe. Precedence H, B/E, G, I, T, S.
    """
    N, CA, C, O = atoms.backbone_arrays()
    n = len(CA)
    if n == 0:
        return ""
    ok = ~np.isnan(np.column_stack([N[:, 0], CA[:, 0], C[:, 0], O[:, 0]])).any(axis=1)
    cont = np.zeros(n, dtype=bool)   # cont[i]: residue i is peptide-bonded to i+1
    if n > 1:
        with np.errstate(invalid="ignore"):
            cn = np.linalg.norm(C[:-1] - N[1:], axis=1)
        cont[:-1] = ok[:-1] & ok[1:] & (cn < 2.0)
    pro = atoms.residues.name == "PRO"
    E = _hbond_energy_matrix(N, CA, C, O, cont)
    with np.errstate(invalid="ignore"):
        Hb = E < hb_cutoff
    Hb[:, pro] = False   # proline has no amide hydrogen
    ss = np.full(n, "-", dtype="<U1")
    # turns
    turn = {}
    for k in (3, 4, 5):
        t = np.zeros(n, dtype=bool)
        for i in range(n - k):
            if Hb[i, i + k] and cont[i:i + k].all():
                t[i] = True
        turn[k] = t
    helix = {}
    for k in (3, 4, 5):
        h = np.zeros(n, dtype=bool)
        for i in range(1, n - k + 1):
            if turn[k][i - 1] and turn[k][i]:
                h[i:i + k] = True
        helix[k] = h
    # bridges
    P = np.zeros((n, n), dtype=bool)
    A = np.zeros((n, n), dtype=bool)
    if n >= 5:
        idx = np.arange(1, n - 1)
        ii, jj = np.meshgrid(idx, idx, indexing="ij")
        good = (cont[ii - 1] & cont[ii] & cont[jj - 1] & cont[jj] & (np.abs(ii - jj) >= 3))
        par = (Hb[ii - 1, jj] & Hb[jj, ii + 1]) | (Hb[jj - 1, ii] & Hb[ii, jj + 1])
        jp1 = np.minimum(jj + 1, n - 1)
        anti = (Hb[ii, jj] & Hb[jj, ii]) | (Hb[ii - 1, jp1] & Hb[jj - 1, ii + 1])
        P[ii, jj] = par & good
        A[ii, jj] = anti & good
    bridge = (P | A).any(axis=1)
    in_ladder = np.zeros(n, dtype=bool)
    if n >= 3:
        pl = P[:-1, :-1] & P[1:, 1:] & cont[:-1, None]
        ii, jj = np.nonzero(pl)
        in_ladder[ii] = in_ladder[ii + 1] = in_ladder[jj] = in_ladder[jj + 1] = True
        al = A[:-1, 1:] & A[1:, :-1] & cont[:-1, None]      # (i,j+1) and (i+1,j)
        ii, jj = np.nonzero(al)
        in_ladder[ii] = in_ladder[ii + 1] = in_ladder[jj + 1] = in_ladder[jj] = True
    # apply from lowest precedence to highest
    kappa = np.zeros(n, dtype=bool)
    for i in range(2, n - 2):
        if cont[i - 2:i + 2].all():
            a, b = CA[i] - CA[i - 2], CA[i + 2] - CA[i]
            cosang = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
            kappa[i] = np.degrees(np.arccos(np.clip(cosang, -1, 1))) > 70.0
    ss[kappa] = "S"
    for k in (3, 4, 5):
        for i in np.nonzero(turn[k])[0]:
            ss[i + 1:i + k] = "T"
    ss[helix[5]] = "I"
    ss[helix[3]] = "G"
    ss[bridge & ~in_ladder] = "B"
    ss[in_ladder] = "E"
    ss[helix[4]] = "H"
    ss[~ok] = "-"
    return "".join(ss.tolist())


def ss3_from_dssp(dssp: str, helix: str = "HGI", strand: str = "E", fill_bulges: bool = True) -> str:
    """Reduce eight-state DSSP to H (helix), E (strand), L (loop). A single non-strand residue flanked by strand
    on both sides becomes E (beta bulge) when ``fill_bulges``."""
    out = ["H" if c in helix else "E" if c in strand else "L" for c in dssp]
    if fill_bulges:
        for i in range(1, len(out) - 1):
            if out[i] == "L" and out[i - 1] == "E" and out[i + 1] == "E" and dssp[i] in "-TSB":
                out[i] = "E"
    return "".join(out)


def assign_ss(atoms: Atoms, helix: str = "HGI", strand: str = "E", fill_bulges: bool = True) -> str:
    """Three-state secondary structure string (H, E, L), one character per residue of ``atoms``.

    DSSP-like (see :func:`assign_dssp`). ``helix='H'`` counts only alpha helix, as the mockup did.
    """
    return ss3_from_dssp(assign_dssp(atoms), helix, strand, fill_bulges)


def ss_fractions(ss: str) -> dict[str, float]:
    n = max(len(ss), 1)
    return {"H": ss.count("H") / n, "E": ss.count("E") / n, "L": ss.count("L") / n}


# --------------------------------------------------------------------------------------------- hydrogen bonds
_DONOR = {("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"), ("LYS", "NZ"), ("ASN", "ND2"), ("GLN", "NE2"),
          ("HIS", "ND1"), ("HIS", "NE2"), ("TRP", "NE1"), ("SER", "OG"), ("THR", "OG1"), ("TYR", "OH")}
_ACCEPTOR = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2"), ("ASN", "OD1"), ("GLN", "OE1"),
             ("HIS", "ND1"), ("HIS", "NE2"), ("SER", "OG"), ("THR", "OG1"), ("TYR", "OH")}
# heavy atom bonded to the acceptor (for the lone-pair direction), second entry optional
_ANTECEDENT = {("ASP", "OD1"): ("CG",), ("ASP", "OD2"): ("CG",), ("GLU", "OE1"): ("CD",), ("GLU", "OE2"): ("CD",),
               ("ASN", "OD1"): ("CG",), ("GLN", "OE1"): ("CD",), ("HIS", "ND1"): ("CG", "CE1"),
               ("HIS", "NE2"): ("CD2", "CE1"), ("SER", "OG"): ("CB",), ("THR", "OG1"): ("CB",),
               ("TYR", "OH"): ("CZ",)}


@dataclass
class HBonds:
    """Hydrogen bonds between two atom tables ``a`` and ``b`` (or within one: then both are the same table).

    One row per directed donor-acceptor pair. ``donor_is_a`` tells which table donates. ``distance`` is D...A;
    ``h_distance`` and ``angle_dha`` (D-H...A) are NaN when no hydrogen position is available;
    ``angle_acc`` is the angle antecedent-A...D at the acceptor (NaN without an antecedent). ``method`` is
    'explicit' (H atoms in the file), 'computed' (backbone amide H from geometry) or 'heavy' (distance and acceptor
    angle only). A Ser/Thr/Tyr pair that could go either way is reported in each direction it passes.
    """
    a_index: np.ndarray
    b_index: np.ndarray
    donor_is_a: np.ndarray
    distance: np.ndarray
    h_distance: np.ndarray
    angle_dha: np.ndarray
    angle_acc: np.ndarray
    method: np.ndarray

    def __len__(self) -> int:
        return len(self.a_index)

    def records(self, a: Atoms, b: Atoms) -> list[dict]:
        out = []
        for k in range(len(self)):
            ia, ib = int(self.a_index[k]), int(self.b_index[k])
            out.append(dict(a_resnum=int(a.resnum[ia]), a_resname=str(a.resname[ia]), a_atom=str(a.name[ia]),
                            b_resnum=int(b.resnum[ib]), b_resname=str(b.resname[ib]), b_atom=str(b.name[ib]),
                            donor="a" if self.donor_is_a[k] else "b", distance=round(float(self.distance[k]), 2),
                            angle=None if np.isnan(self.angle_dha[k]) else round(float(self.angle_dha[k]), 1),
                            method=str(self.method[k])))
        return out

    def residue_pairs(self, a: Atoms, b: Atoms) -> list[tuple[int, int]]:
        """Unique (a residue number, b residue number) pairs."""
        return sorted({(int(a.resnum[i]), int(b.resnum[j])) for i, j in zip(self.a_index, self.b_index)})


def _polar_roles(at: Atoms):
    """Donor mask, acceptor mask, and antecedent atom indices (n, 2; -1 absent) for acceptors."""
    n = len(at)
    don = np.zeros(n, dtype=bool)
    acc = np.zeros(n, dtype=bool)
    ante = np.full((n, 2), -1, dtype=np.int64)
    heavy = ~at.is_hydrogen
    r = at.residues
    pro = at.resname == "PRO"
    don |= heavy & (at.name == "N") & ~pro
    acc |= heavy & np.isin(at.name, ("O", "OXT"))
    for i in np.nonzero(heavy & (at.element != "C"))[0]:
        rn, nm = str(at.resname[i]), str(at.name[i])
        if (rn, nm) in _DONOR:
            don[i] = True
        if (rn, nm) in _ACCEPTOR:
            acc[i] = True
    # antecedents: backbone O -> C, side chain per table
    need = np.nonzero(acc)[0]
    if len(need):
        names_per_res: dict[int, dict[str, int]] = {}
        for ri in np.unique(at.res_index[need]):
            s, e = r.start[ri], r.stop[ri]
            names_per_res[int(ri)] = {str(nm): int(s + k) for k, nm in enumerate(at.name[s:e]) if at.element[s + k] != "H"}
        for i in need:
            d = names_per_res[int(at.res_index[i])]
            nm = str(at.name[i])
            if nm in ("O", "OXT"):
                if "C" in d:
                    ante[i, 0] = d["C"]
            else:
                for k, an in enumerate(_ANTECEDENT.get((str(at.resname[i]), nm), ())):
                    if an in d:
                        ante[i, k] = d[an]
    return don, acc, ante


def _donor_hydrogens(at: Atoms, don_idx: np.ndarray):
    """For each donor index: (list of explicit H atom indices, computed H xyz or None)."""
    hmask = at.is_hydrogen
    hyd_idx = np.nonzero(hmask)[0]
    explicit: dict[int, list[int]] = {}
    if len(hyd_idx):
        ht = cKDTree(at.xyz[hyd_idx])
        for i, lst in zip(don_idx, cKDTree(at.xyz[don_idx]).query_ball_tree(ht, 1.25)):
            if lst:
                explicit[int(i)] = [int(hyd_idx[j]) for j in lst]
    return explicit


def _computed_amide_h(at: Atoms) -> dict[int, np.ndarray]:
    """Backbone N atom index -> hydrogen position from C(i-1)=O(i-1), for residues bonded to the previous one."""
    N, CA, C, O = at.backbone_arrays()
    r = at.residues
    out = {}
    n_idx = np.full(len(r), -1)
    nn = np.nonzero((at.name == "N") & ~at.is_hydrogen)[0]
    n_idx[at.res_index[nn][::-1]] = nn[::-1]
    for i in range(1, len(r)):
        if n_idx[i] < 0 or r.name[i] == "PRO" or np.isnan(C[i - 1, 0]) or np.isnan(O[i - 1, 0]) or np.isnan(N[i, 0]):
            continue
        if np.linalg.norm(C[i - 1] - N[i]) > 2.0:
            continue
        v = C[i - 1] - O[i - 1]
        out[int(n_idx[i])] = N[i] + v / np.linalg.norm(v)
    return out


def _angle(a, b, c) -> np.ndarray:
    """Angle at b (degrees) for point arrays a-b-c."""
    u, v = a - b, c - b
    cos = (u * v).sum(1) / (np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1))
    return np.degrees(np.arccos(np.clip(cos, -1, 1)))


def hbonds(a: Atoms, b: Atoms | None = None, max_dist: float = 3.5, max_h_dist: float = 2.5, min_angle_dha: float = 120.0,
           min_angle_acc: float = 90.0, use_hydrogens: bool = True) -> HBonds:
    """Hydrogen bonds between ``a`` and ``b`` (or inside ``a`` when ``b`` is None).

    Protein polar atoms only. Donors: backbone N (not Pro), Arg NE/NH1/NH2, Lys NZ, Asn ND2, Gln NE2, His ND1/NE2,
    Trp NE1, Ser OG, Thr OG1, Tyr OH. Acceptors: backbone O and OXT, Asp OD, Glu OE, Asn OD1, Gln OE1, His ND1/NE2,
    Ser OG, Thr OG1, Tyr OH. Criteria: D...A <= ``max_dist``; with a hydrogen (explicit in the table, or the backbone
    amide H placed from geometry) H...A <= ``max_h_dist`` and D-H...A >= ``min_angle_dha`` for at least one H;
    otherwise (side-chain donors without explicit H) the acceptor angle antecedent-A...D >= ``min_angle_acc``.
    Within one table, atoms of the same residue and backbone N(i)...O(i-1) peptide neighbours are skipped.
    """
    same = b is None
    b = a if same else b
    da, aa, ante_a = _polar_roles(a)
    db, ab, ante_b = _polar_roles(b)
    expl_a = _donor_hydrogens(a, np.nonzero(da)[0]) if use_hydrogens and a.is_hydrogen.any() else {}
    expl_b = expl_a if same else (_donor_hydrogens(b, np.nonzero(db)[0]) if use_hydrogens and b.is_hydrogen.any() else {})
    comp_a = _computed_amide_h(a) if use_hydrogens else {}
    comp_b = comp_a if same else (_computed_amide_h(b) if use_hydrogens else {})

    rows = []

    def one_direction(dt: Atoms, d_mask, expl, comp, at: Atoms, a_mask, ante, donor_is_a: bool):
        di = np.nonzero(d_mask)[0]
        ai = np.nonzero(a_mask)[0]
        if not len(di) or not len(ai):
            return
        tree_a = cKDTree(at.xyz[ai])
        for dd, lst in zip(di, cKDTree(dt.xyz[di]).query_ball_tree(tree_a, max_dist)):
            for k in lst:
                aj = int(ai[k])
                if same:
                    if dt.res_index[dd] == at.res_index[aj]:
                        continue
                    if (abs(int(dt.res_index[dd]) - int(at.res_index[aj])) == 1 and dt.name[dd] == "N"
                            and at.name[aj] == "O" and at.res_index[aj] == dt.res_index[dd] - 1):
                        continue
                dpos, apos = dt.xyz[dd], at.xyz[aj]
                dist = float(np.linalg.norm(dpos - apos))
                # acceptor angle
                aidx = ante[aj]
                aidx = aidx[aidx >= 0]
                ang_acc = np.nan
                if len(aidx):
                    ante_pos = at.xyz[aidx].mean(0)
                    ang_acc = float(_angle(ante_pos[None], apos[None], dpos[None])[0])
                hs = [dt.xyz[h] for h in expl.get(int(dd), [])]
                method = "explicit"
                if not hs and int(dd) in comp:
                    hs = [comp[int(dd)]]
                    method = "computed"
                if hs:
                    hs_arr = np.array(hs)
                    hd = np.linalg.norm(hs_arr - apos, axis=1)
                    ang = _angle(np.broadcast_to(dpos, hs_arr.shape), hs_arr, np.broadcast_to(apos, hs_arr.shape))
                    ok = (hd <= max_h_dist) & (ang >= min_angle_dha)
                    if not ok.any():
                        continue
                    best = np.argmin(np.where(ok, hd, np.inf))
                    h_dist, ang_dha = float(hd[best]), float(ang[best])
                else:
                    method = "heavy"
                    if not np.isnan(ang_acc) and ang_acc < min_angle_acc:
                        continue
                    h_dist, ang_dha = np.nan, np.nan
                rows.append((dd if donor_is_a else aj, aj if donor_is_a else dd, donor_is_a, dist, h_dist, ang_dha,
                             ang_acc, method))

    one_direction(a, da, expl_a, comp_a, b, ab, ante_b, True)
    if not same:   # inside one table every donor to acceptor pair is already covered
        one_direction(b, db, expl_b, comp_b, a, aa, ante_a, False)
    if not rows:
        z = np.zeros(0, dtype=np.int64)
        e = np.zeros(0)
        return HBonds(z, z, np.zeros(0, dtype=bool), e, e, e, e, np.array([], dtype="<U8"))
    ai_, bi_, dia, dist, hd, ang, aacc, meth = zip(*rows)
    return HBonds(np.array(ai_, dtype=np.int64), np.array(bi_, dtype=np.int64), np.array(dia, dtype=bool),
                  np.array(dist), np.array(hd), np.array(ang), np.array(aacc), np.array(meth, dtype="<U8"))


@dataclass
class SaltBridge:
    a_resnum: int
    a_resname: str
    b_resnum: int
    b_resname: str
    distance: float
    positive_in_a: bool


def _charged_centroids(at: Atoms, include_his: bool = False):
    """(centroids (n,3), resnums, resnames, sign) of side-chain charged groups."""
    groups = {"ARG": (("NE", "NH1", "NH2", "CZ"), +1), "LYS": (("NZ",), +1), "ASP": (("OD1", "OD2"), -1),
              "GLU": (("OE1", "OE2"), -1)}
    if include_his:
        groups["HIS"] = (("ND1", "NE2"), +1)
    r = at.residues
    cents, nums, names, sign = [], [], [], []
    for i in range(len(r)):
        g = groups.get(str(r.name[i]))
        if g is None:
            continue
        s, e = r.start[i], r.stop[i]
        m = np.isin(at.name[s:e], g[0]) & ~at.is_hydrogen[s:e]
        if m.any():
            cents.append(at.xyz[s:e][m].mean(0))
            nums.append(int(r.number[i]))
            names.append(str(r.name[i]))
            sign.append(g[1])
    return (np.array(cents).reshape(-1, 3), np.array(nums, dtype=np.int64), np.array(names), np.array(sign))


def salt_bridges(a: Atoms, b: Atoms | None = None, cutoff: float = 4.0, include_his: bool = False) -> list[SaltBridge]:
    """Opposite-charge side-chain groups whose centroids are within ``cutoff`` A.

    Charged groups: Arg guanidinium (NE, NH1, NH2, CZ), Lys NZ, Asp carboxylate (OD1, OD2), Glu carboxylate
    (OE1, OE2); His optional. Between two tables, or inside one (each pair once).
    """
    same = b is None
    ca, na, la, sa = _charged_centroids(a, include_his)
    cb, nb, lb, sb = (ca, na, la, sa) if same else _charged_centroids(b, include_his)
    if not len(ca) or not len(cb):
        return []
    out = []
    for i, lst in enumerate(cKDTree(ca).query_ball_tree(cKDTree(cb), cutoff)):
        for j in lst:
            if sa[i] * sb[j] >= 0:
                continue
            if same and (j <= i or na[i] == nb[j]):
                continue
            out.append(SaltBridge(int(na[i]), str(la[i]), int(nb[j]), str(lb[j]),
                                  float(np.linalg.norm(ca[i] - cb[j])), bool(sa[i] > 0)))
    return out


# --------------------------------------------------------------------------------------------- fold similarity
@dataclass
class TMResult:
    """TM-align result. ``tm_shorter`` is normalised by the shorter chain (the larger score), ``tm_longer`` by the
    longer chain; ``tm_1``/``tm_2`` by chain 1 / chain 2. ``R``, ``t`` move chain 1 onto chain 2."""
    tm_shorter: float
    tm_longer: float
    tm_1: float
    tm_2: float
    rmsd: float
    aligned_length: int
    R: np.ndarray
    t: np.ndarray
    len_1: int
    len_2: int


def _clean_seq(seq: str | None, n: int) -> str:
    if seq is None:
        return "A" * n
    s = "".join(c if c.isalpha() else "X" for c in seq.upper())
    if len(s) != n:
        raise ValueError(f"sequence length {len(s)} does not match {n} coordinates")
    return s


def tm_align(xyz1: np.ndarray, seq1: str | None, xyz2: np.ndarray, seq2: str | None) -> TMResult:
    """TM-align of two CA traces (tmtools). Sequences are used only for the sequence-identity-seeded alignment and
    may be None (all A)."""
    from tmtools import tm_align as _tm
    x1 = np.ascontiguousarray(xyz1, dtype=np.float64)
    x2 = np.ascontiguousarray(xyz2, dtype=np.float64)
    res = _tm(x1, x2, _clean_seq(seq1, len(x1)), _clean_seq(seq2, len(x2)))
    l1, l2 = len(x1), len(x2)
    t1, t2 = float(res.tm_norm_chain1), float(res.tm_norm_chain2)
    short = t1 if l1 <= l2 else t2
    long_ = t2 if l1 <= l2 else t1
    aligned = sum(1 for c in res.seqM if c in ":.")
    return TMResult(short, long_, t1, t2, float(res.rmsd), aligned, np.array(res.u), np.array(res.t), l1, l2)


def tm_align_atoms(a: Atoms, b: Atoms) -> TMResult:
    """TM-align the CA traces of two tables (residues without CA are skipped)."""
    ca_a, ca_b = a.ca_xyz(), b.ca_xyz()
    ma, mb = ~np.isnan(ca_a[:, 0]), ~np.isnan(ca_b[:, 0])
    sa, sb = np.array(list(a.sequence))[ma], np.array(list(b.sequence))[mb]
    return tm_align(ca_a[ma], "".join(sa), ca_b[mb], "".join(sb))


_TM_STATE: dict = {}


def _tm_init(coords, seqs):
    _TM_STATE["c"], _TM_STATE["s"] = coords, seqs


def _tm_row(i: int):
    from tmtools import tm_align as _tm
    c, s = _TM_STATE["c"], _TM_STATE["s"]
    n = len(c)
    short = np.zeros(n)
    long_ = np.zeros(n)
    rms = np.zeros(n)
    nal = np.zeros(n, dtype=np.int64)
    li = len(c[i])
    for j in range(i + 1, n):
        r = _tm(c[i], c[j], s[i], s[j])
        lj = len(c[j])
        t1, t2 = float(r.tm_norm_chain1), float(r.tm_norm_chain2)
        short[j] = t1 if li <= lj else t2
        long_[j] = t2 if li <= lj else t1
        rms[j] = float(r.rmsd)
        nal[j] = sum(1 for ch in r.seqM if ch in ":.")
    return i, short, long_, rms, nal


@dataclass
class TMMatrix:
    """Symmetric pairwise TM-align matrices (diagonal TM 1, RMSD 0). ``tm`` is the shorter-chain normalisation."""
    tm: np.ndarray
    tm_longer: np.ndarray
    rmsd: np.ndarray
    aligned: np.ndarray

    @property
    def n(self) -> int:
        return len(self.tm)


def tm_matrix(chains: Sequence, n_jobs: int | None = None, progress: Callable[[int, int], None] | None = None) -> TMMatrix:
    """All-against-all TM-align. ``chains`` items are (xyz, seq) tuples, bare (n, 3) arrays, or :class:`Atoms`.

    Each unordered pair is aligned once and mirrored, so the result is exactly symmetric. ``n_jobs`` processes
    (default: CPU count minus one, 1 runs in-process); ``progress(done_rows, total_rows)`` is called as rows finish.
    """
    import multiprocessing as mp
    import os
    coords, seqs = [], []
    for ch in chains:
        if isinstance(ch, Atoms):
            ca = ch.ca_xyz()
            m = ~np.isnan(ca[:, 0])
            coords.append(np.ascontiguousarray(ca[m]))
            seqs.append("".join(np.array(list(ch.sequence))[m]))
        elif isinstance(ch, tuple):
            x = np.ascontiguousarray(ch[0], dtype=np.float64)
            coords.append(x)
            seqs.append(_clean_seq(ch[1], len(x)))
        else:
            x = np.ascontiguousarray(ch, dtype=np.float64)
            coords.append(x)
            seqs.append("A" * len(x))
    n = len(coords)
    short = np.zeros((n, n))
    long_ = np.zeros((n, n))
    rms = np.zeros((n, n))
    nal = np.zeros((n, n), dtype=np.int64)
    if n == 0:
        return TMMatrix(short, long_, rms, nal)
    jobs = max(1, (os.cpu_count() or 2) - 1) if n_jobs is None else max(1, int(n_jobs))
    jobs = min(jobs, max(1, n - 1))

    def fill(res):
        i, s, l, r, a = res
        short[i], long_[i], rms[i], nal[i] = s, l, r, a

    done = 0
    if jobs == 1 or n < 4:
        _tm_init(coords, seqs)
        for i in range(n - 1):
            fill(_tm_row(i))
            done += 1
            if progress:
                progress(done, n - 1)
    else:
        ctx = mp.get_context("spawn")
        with ctx.Pool(jobs, initializer=_tm_init, initargs=(coords, seqs)) as pool:
            # long rows first: row i has n-1-i pairs
            for res in pool.imap_unordered(_tm_row, range(n - 1), chunksize=1):
                fill(res)
                done += 1
                if progress:
                    progress(done, n - 1)
    for M in (short, long_, rms):
        M += M.T
    nal += nal.T
    np.fill_diagonal(short, 1.0)
    np.fill_diagonal(long_, 1.0)
    np.fill_diagonal(nal, [len(c) for c in coords])
    return TMMatrix(short, long_, rms, nal)


# --------------------------------------------------------------------------------------------- sequence helpers
def _clean(seq: str) -> str:
    return "".join(c for c in seq.upper() if c in AA_ORDER)


def gravy(seq: str, scale: dict[str, float] | None = None) -> float:
    """Grand average of hydropathy (Kyte-Doolittle by default); 0.0 for an empty sequence."""
    sc = KYTE_DOOLITTLE if scale is None else scale
    s = [c for c in seq.upper() if c in sc]
    return float(sum(sc[c] for c in s) / len(s)) if s else 0.0


def net_charge(seq: str, ph: float = 7.4) -> float:
    """Net charge by Henderson-Hasselbalch with the pKa set in :data:`PKA` (the mockup's values)."""
    s = _clean(seq)
    if not s:
        return 0.0
    pos = 1.0 / (1.0 + 10 ** (ph - PKA["Nterm"])) + sum(1.0 / (1.0 + 10 ** (ph - PKA[a])) for a in s if a in "KRH")
    neg = 1.0 / (1.0 + 10 ** (PKA["Cterm"] - ph)) + sum(1.0 / (1.0 + 10 ** (PKA[a] - ph)) for a in s if a in "DECY")
    return pos - neg


def isoelectric_point(seq: str) -> float:
    """pH at which :func:`net_charge` is zero (bisection)."""
    lo, hi = 0.0, 14.0
    for _ in range(50):
        mid = (lo + hi) / 2
        if net_charge(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def sequence_hydrophobic_fraction(seq: str, residues: Iterable[str] = BINDCRAFT_APOLAR) -> float:
    s = _clean(seq)
    return sum(1 for c in s if c in residues) / len(s) if s else 0.0


def _virtual_cb(N: np.ndarray, CA: np.ndarray, C: np.ndarray) -> np.ndarray:
    b = CA - N
    c = C - CA
    a = np.cross(b, c)
    return -0.58273431 * a + 0.56802827 * b - 0.54067466 * c + CA


def sc_neighbors(atoms: Atoms, dist_midpoint: float = 9.0, angle_shift: float = 0.5, angle_exponent: float = 2.0) -> np.ndarray:
    """Rosetta-style side-chain neighbour count per residue (the burial measure behind LayerSelector).

    For each residue the sum over every other residue j of ``dist_term * angle_term``: dist_term is the logistic
    1 / (1 + exp(d - 9)) of the distance between neighbour atoms (CB, or CA for Gly), angle_term is
    ((cos(theta) + 0.5) / 1.5)^2 clipped at 0 with theta the angle between the residue's CA to CB vector (virtual
    CB for Gly) and the vector to residue j. Reimplemented from memory of Rosetta defaults and checked
    against BindCraft's Surface_Hydrophobicity values, not against Rosetta itself. Residues lacking a backbone give NaN.
    """
    N, CA, C, _ = atoms.backbone_arrays()
    r = atoms.residues
    n = len(r)
    cb = np.full((n, 3), np.nan)
    idx = np.nonzero((atoms.name == "CB") & ~atoms.is_hydrogen)[0]
    cb[atoms.res_index[idx][::-1]] = atoms.xyz[idx][::-1]
    gly = r.name == "GLY"
    with np.errstate(invalid="ignore"):
        vcb = _virtual_cb(N, CA, C)
    cb = np.where(np.isnan(cb[:, :1]), vcb, cb)
    nbr = np.where(gly[:, None], CA, cb)
    v = cb - CA
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    out = np.full(n, np.nan)
    ok = ~np.isnan(nbr[:, 0]) & ~np.isnan(v[:, 0])
    pts = nbr[ok]
    vv = v[ok]
    w = pts[None, :, :] - pts[:, None, :]
    dist = np.linalg.norm(w, axis=2)
    np.fill_diagonal(dist, np.inf)
    cos = np.einsum("ijk,ik->ij", w, vv) / dist
    ang = np.clip((cos + angle_shift) / (1.0 + angle_shift), 0.0, None) ** angle_exponent
    term = ang / (1.0 + np.exp(np.clip(dist - dist_midpoint, -60, 60)))
    np.fill_diagonal(term, 0.0)
    out[ok] = term.sum(1)
    return out


def surface_hydrophobicity(binder: Atoms, surface_cutoff: float = 2.0, hydrophobic: Iterable[str] = BINDCRAFT_APOLAR) -> float:
    """Approximation of BindCraft's Surface_Hydrophobicity: fraction of surface residues of the binder (alone) that
    are hydrophobic. Surface = side-chain neighbour count below ``surface_cutoff`` (Rosetta LayerSelector default 2.0),
    hydrophobic = Rosetta apolar types plus F, W, Y. Approximate: matches the CSV to a mean absolute difference of
    about 0.01 (see docs/structure_api.md); BindCraft itself runs Rosetta."""
    scn = sc_neighbors(binder)
    surf = scn < surface_cutoff
    seq = binder.residues.one
    hyd = np.isin(seq, list(hydrophobic))
    return float((hyd & surf).sum() / surf.sum()) if surf.any() else 0.0
