"""Shared protein–small-molecule geometry for pose and trajectory analysis.

Cutoffs and residue rules follow the bundled MN-Ligand native MD analysis.
Coordinates and distances are in Å. Chemistry availability is explicit.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


KINDS = ("contact", "hydrophobic", "hydrogen_bond", "water_bridge", "salt_bridge")
HYDROPHOBIC = {
    "ALA": "CB", "VAL": "CB CG1 CG2", "LEU": "CB CG CD1 CD2",
    "ILE": "CB CG1 CG2 CD1", "MET": "CB CG SD CE",
    "PHE": "CB CG CD1 CD2 CE1 CE2 CZ", "TYR": "CB CG CD1 CD2 CE1 CE2 CZ",
    "TRP": "CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2", "PRO": "CB CG CD",
}
PROTEIN_CHARGES = {
    "ASP": (-1, "OD1 OD2"), "GLU": (-1, "OE1 OE2"), "LYS": (1, "NZ"),
    "ARG": (1, "NE NH1 NH2"), "HIP": (1, "ND1 NE2"), "HSP": (1, "ND1 NE2"),
}
PROTEIN_DONORS = {
    "ARG": "NE NH1 NH2", "ASN": "ND2", "GLN": "NE2", "HIS": "ND1 NE2",
    "HIP": "ND1 NE2", "HSP": "ND1 NE2", "LYS": "NZ", "SER": "OG",
    "THR": "OG1", "TRP": "NE1", "TYR": "OH",
}
PROTEIN_ACCEPTORS = {
    "ASN": "OD1", "ASP": "OD1 OD2", "GLN": "OE1", "GLU": "OE1 OE2",
    "HIS": "ND1 NE2", "SER": "OG", "THR": "OG1", "TYR": "OH",
}


@dataclass(frozen=True)
class Atom:
    name: str
    element: str
    residue: str
    residue_name: str


class Classifier:
    def __init__(self, atoms, protein, ligand, waters, donors, acceptors, ligand_charges,
                 *, cutoff=4.5, pose=False):
        self.atoms = atoms
        self.protein = np.asarray(protein, dtype=int)
        self.ligand = np.asarray(ligand, dtype=int)
        self.waters = np.asarray(waters, dtype=int)
        self.donors = donors  # donor index -> bonded H indices, empty for pose candidates
        self.acceptors = acceptors
        self.cutoff = cutoff
        self.pose = pose
        self.hydrophobic = {i for i in self.protein
                            if atoms[i].name.upper() in HYDROPHOBIC.get(atoms[i].residue_name, "").split()}
        self.protein_charges = {}
        for index in self.protein:
            atom = atoms[index]
            charge, names = PROTEIN_CHARGES.get(atom.residue_name, (0, ""))
            if atom.name.upper() in names.split():
                self.protein_charges[index] = charge
        self.ligand_charges = ligand_charges or {}
        self.charge_source = None
        ligand_polar = set(self.ligand) & (self.acceptors | self.donors.keys())
        protein_polar = set(self.protein) & (self.acceptors | self.donors.keys())
        bridge_hydrogens = any(index in donors for index in self.waters) or (
            any(index in donors for index in ligand_polar) and any(index in donors for index in protein_polar)
        )
        self.methods = {
            "contact": f"Heavy-atom distance ≤ {cutoff:g} Å",
            "hydrophobic": "MN-Ligand hydrophobic protein atoms and ligand C/S within 4 Å",
            "hydrogen_bond": (
                "Wernet–Nilsson with bonded hydrogens; typed donor/acceptor distance ≤ 3.5 Å where hydrogen coordinates are absent"
                if pose else "Wernet–Nilsson with bonded N/O hydrogens and topology N/O acceptor sites"
            ) if ((pose and ligand_polar and protein_polar) or
                  (any(index in donors for index in self.ligand) and protein_polar) or
                  (any(index in donors for index in self.protein) and ligand_polar)) else None,
            "water_bridge": (
                "Two compatible hydrogen-bond legs through the same water oxygen; distance-only legs flagged"
                if pose else "Two Wernet–Nilsson hydrogen bonds through the same water oxygen"
            ) if len(self.waters) and ligand_polar and protein_polar and (pose or bridge_hydrogens) else None,
            "salt_bridge": "Opposite formal charges within 4 Å" if ligand_charges is not None else None,
        }

    def _hydrogen_bond(self, xyz, left, right, distance):
        evidence = None
        for donor, acceptor in ((left, right), (right, left)):
            if donor not in self.donors or acceptor not in self.acceptors:
                continue
            hydrogens = self.donors[donor]
            if hydrogens:
                vectors = xyz[list(hydrogens)] - xyz[donor]
                lengths = np.linalg.norm(vectors, axis=1)
                valid = lengths > 0
                cosines = vectors[valid] @ (xyz[acceptor] - xyz[donor]) / (lengths[valid] * distance)
                angles = np.degrees(np.arccos(np.clip(cosines, -1.0, 1.0)))
                if np.any((angles < 45) & (distance < 3.3 - 0.00044 * angles**2)):
                    return "bonded_hydrogen"
            elif self.pose and distance <= 3.5:
                evidence = "distance_only"
        return evidence

    def frame(self, xyz, *, box=None):
        """Return atom-pair evidence for one coordinate set."""
        xyz = np.asarray(xyz, dtype=float)
        water_coordinates = xyz[self.waters]
        if box is not None and self.methods["water_bridge"]:
            # The protein and ligand have already been imaged by MDTraj.
            # Use the closest water images around the selected ligand.
            center = xyz[self.ligand].mean(axis=0)
            shifts = np.round((water_coordinates - center) @ np.linalg.inv(box)) @ box
            water_coordinates = water_coordinates - shifts
        protein_tree = cKDTree(xyz[self.protein])
        ligand_tree = cKDTree(xyz[self.ligand])
        pairs = protein_tree.sparse_distance_matrix(
            ligand_tree, max(self.cutoff, 4.0), output_type="ndarray"
        )
        hits = {}
        for local_protein, local_ligand, distance in pairs:
            protein = int(self.protein[local_protein])
            ligand = int(self.ligand[local_ligand])
            if distance <= self.cutoff:
                hits["contact", protein, ligand, "distance"] = distance
            if distance <= 4.0 and protein in self.hydrophobic and self.atoms[ligand].element in {"C", "S"}:
                hits["hydrophobic", protein, ligand, "distance"] = distance
            if distance <= 4.0 and self.protein_charges.get(protein, 0) * self.ligand_charges.get(ligand, 0) < 0:
                hits["salt_bridge", protein, ligand, "opposite_formal_charges"] = distance
            if self.methods["hydrogen_bond"] and 0 < distance <= 3.5:
                if evidence := self._hydrogen_bond(xyz, protein, ligand, distance):
                    hits["hydrogen_bond", protein, ligand, evidence] = distance

        if self.methods["water_bridge"]:
            neighbors = cKDTree(water_coordinates).query_ball_point(xyz[self.ligand], 3.5)
            by_water = defaultdict(list)
            for ligand, indices in zip(self.ligand, neighbors):
                if ligand in self.donors or ligand in self.acceptors:
                    for local_water in indices:
                        by_water[local_water].append(ligand)
            if box is not None and by_water:
                xyz = xyz.copy()
                for local_water in by_water:
                    water = self.waters[local_water]
                    xyz[water] = water_coordinates[local_water]
                    if water in self.donors:
                        xyz[list(self.donors[water])] -= shifts[local_water]
            for local_water, ligands in by_water.items():
                water = self.waters[local_water]
                ligand_legs = self._water_legs(xyz, water, ligands)
                if not ligand_legs:
                    continue
                protein_indices = protein_tree.query_ball_point(xyz[water], 3.5)
                protein_legs = self._water_legs(xyz, water, self.protein[protein_indices])
                for protein, first_evidence, first_distance in protein_legs:
                    for ligand, second_evidence, second_distance in ligand_legs:
                        evidence = "distance_only" if "distance_only" in (first_evidence, second_evidence) else "bonded_hydrogen"
                        key = ("water_bridge", protein, ligand, evidence)
                        distance = max(first_distance, second_distance)
                        hits[key] = min(hits.get(key, distance), distance)
        return hits

    def _water_legs(self, xyz, water, atoms):
        legs = []
        for index in atoms:
            distance = float(np.linalg.norm(xyz[water] - xyz[index]))
            if distance > 0 and (evidence := self._hydrogen_bond(xyz, water, index, distance)):
                legs.append((int(index), evidence, distance))
        return legs


class Occupancy:
    """Count unique pair and residue events per frame using bounded accumulators."""

    def __init__(self, classifier):
        self.classifier = classifier
        self.frames = 0
        self.pairs = Counter()
        self.minimum = {}
        self.residues = Counter()

    def update(self, hits):
        self.frames += 1
        events = set()
        for key, distance in hits.items():
            kind, protein, _ligand, _evidence = key
            self.pairs[key] += 1
            self.minimum[key] = min(self.minimum.get(key, distance), distance)
            events.add((self.classifier.atoms[protein].residue, kind))
        self.residues.update(events)
        counts = Counter(kind for _residue, kind in events)
        return {kind + "_residues": counts[kind] if self.classifier.methods[kind] else None for kind in KINDS}

    def pair_rows(self):
        atoms = self.classifier.atoms
        for (kind, protein, ligand, evidence), frames in sorted(self.pairs.items()):
            yield {"interaction": kind, "protein_residue": atoms[protein].residue,
                   "protein_atom": atoms[protein].name, "ligand_atom": atoms[ligand].name,
                   "protein_atom_index": protein, "ligand_atom_index": ligand, "evidence": evidence,
                   "minimum_distance_angstrom": self.minimum[kind, protein, ligand, evidence],
                   "frames": frames, "analyzed_frames": self.frames,
                   "occupancy_percent": 100 * frames / self.frames}

    def report(self):
        best = defaultdict(list)
        for row in self.pair_rows():
            best[row["protein_residue"], row["interaction"]].append(row)
        rows = []
        for residue in dict.fromkeys(self.classifier.atoms[i].residue for i in self.classifier.protein):
            row = {"residue": residue, "atom_endpoints": {}}
            for kind, method in self.classifier.methods.items():
                row[kind + "_occupancy"] = self.residues[residue, kind] / self.frames if method else None
                if options := best.get((residue, kind)):
                    pair = min(options, key=lambda item: (-item["frames"], item["minimum_distance_angstrom"]))
                    row["atom_endpoints"][kind] = {
                        label: pair[label] for label in ("protein_atom", "ligand_atom", "evidence")
                    }
            endpoint = next((row["atom_endpoints"][kind] for kind in
                             ("hydrogen_bond", "salt_bridge", "water_bridge", "hydrophobic", "contact")
                             if kind in row["atom_endpoints"]), {})
            rows.append({**row, "top_ligand_atom": endpoint.get("ligand_atom"),
                         "top_protein_atom": endpoint.get("protein_atom")})
        return {"method_version": 1, "methods": self.classifier.methods,
                "formal_charge_source": self.classifier.charge_source,
                "analysis_window": {"analyzed_frames": self.frames},
                "contacts": {"residues": rows}}


PAIR_COLUMNS = (
    "interaction", "protein_residue", "protein_atom", "ligand_atom",
    "protein_atom_index", "ligand_atom_index", "evidence",
    "minimum_distance_angstrom", "frames", "analyzed_frames", "occupancy_percent",
)


def trajectory_classifier(topology, topology_path, protein, ligand, residue_label, cutoff, *,
                          formal_charge_override=None, ligand_smiles=None):
    atoms = [Atom(atom.name, getattr(atom.element, "symbol", "").upper(),
                  residue_label(atom.residue), atom.residue.name.upper()) for atom in topology.atoms]
    waters = [i for i, atom in enumerate(topology.atoms) if atom.residue.is_water and atoms[i].element == "O"]
    sites = set(protein) | set(ligand) | set(waters)
    acceptors = {i for i in sites if atoms[i].element in {"N", "O"}}
    donors = defaultdict(list)
    for left, right in topology.bonds:
        i, j = left.index, right.index
        if atoms[i].element == "H":
            i, j = j, i
        if i in acceptors and atoms[j].element == "H":
            donors[i].append(j)
    charge_source = None
    if formal_charge_override is None:
        charges = _explicit_topology_charges(topology_path, topology, ligand)
        if charges is not None:
            charge_source = "Explicit topology atom charges."
        elif ligand_smiles:
            charges = _matched_ligand_charges(topology, ligand, ligand_smiles)
            if charges is not None:
                charge_source = "Saved study ligand identity matched to topology bonds and hydrogen counts."
    else:
        by_name = defaultdict(list)
        for index in ligand:
            by_name[atoms[index].name].append(int(index))
        if any(len(by_name[name]) != 1 for name in formal_charge_override):
            raise ValueError("Formal charge atom names must each match exactly one selected ligand heavy atom.")
        charges = {by_name[name][0]: value for name, value in formal_charge_override.items()}
        charge_source = "User-declared ligand atom formal charges."

    acceptors -= {i for i, charge in (charges or {}).items() if charge > 0}
    classifier = Classifier(atoms, protein, ligand, waters, dict(donors), acceptors, charges, cutoff=cutoff)
    classifier.charge_source = charge_source
    return Occupancy(classifier)


def _explicit_topology_charges(path, topology, ligand):
    """Use formal charges only with an exact source atom identity."""
    if Path(path).suffix.lower() not in {".pdb", ".ent", ".cif", ".mmcif"}:
        return None
    import gemmi

    structure = gemmi.read_structure(str(path))
    if not structure:
        return None
    records = defaultdict(list)
    for chain in structure[0]:
        for residue in chain:
            for atom in residue:
                key = (chain.name.strip(), residue.seqid.num, residue.name, atom.name.strip())
                records[key].append(int(atom.charge))
    charges = {}
    for index in ligand:
        atom = topology.atom(int(index))
        residue = atom.residue
        key = (str(getattr(residue.chain, "chain_id", "") or "").strip(),
               residue.resSeq, residue.name, atom.name.strip())
        values = records.get(key, ())
        if len(values) == 1 and values[0]:
            charges[int(index)] = values[0]
    # Blank PDB charge fields provide no evidence about ligand formal charge.
    return charges or None


def _matched_ligand_charges(topology, ligand, smiles):
    """Resolve charges only if the saved identity and actual topology agree."""
    from rdkit import Chem

    if not isinstance(smiles, str) or len(smiles) > 4096 or len(ligand) > 300:
        return None
    reference = Chem.MolFromSmiles(smiles)
    if reference is None or reference.GetNumAtoms() != len(ligand):
        return None
    declared = Chem.RWMol()
    observed = Chem.RWMol()
    for ref, index in zip(reference.GetAtoms(), ligand):
        element = topology.atom(int(index)).element
        if element is None:
            return None
        declared.AddAtom(Chem.Atom(ref.GetAtomicNum()))
        observed.AddAtom(Chem.Atom(element.atomic_number))
    for bond in reference.GetBonds():
        declared.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), Chem.BondType.SINGLE)

    local = {index: i for i, index in enumerate(ligand)}
    hydrogen_counts = Counter()
    for left, right in topology.bonds:
        left_local, right_local = local.get(left.index), local.get(right.index)
        if left_local is not None and right_local is not None:
            observed.AddBond(left_local, right_local, Chem.BondType.SINGLE)
        elif left_local is not None and getattr(right.element, "symbol", "") == "H":
            hydrogen_counts[left_local] += 1
        elif right_local is not None and getattr(left.element, "symbol", "") == "H":
            hydrogen_counts[right_local] += 1
    if observed.GetNumBonds() != declared.GetNumBonds():
        return None
    maps = observed.GetSubstructMatches(declared, uniquify=False, maxMatches=257)
    if not maps or len(maps) > 256:
        return None
    reference_atoms = list(reference.GetAtoms())
    charge_options = []
    for mapping in maps:
        if all(atom.GetTotalNumHs() == hydrogen_counts[index] for atom, index in zip(reference_atoms, mapping)):
            charge_options.append({int(ligand[index]): atom.GetFormalCharge()
                                   for atom, index in zip(reference_atoms, mapping) if atom.GetFormalCharge()})
    if not charge_options or any(option != charge_options[0] for option in charge_options):
        return None
    return charge_options[0]


@lru_cache(maxsize=1)
def _rdkit_features():
    from rdkit import RDConfig
    from rdkit.Chem import ChemicalFeatures

    return ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))


def classify_pose(receptor_path, ligand_path):
    """Classify one explicitly selected ligand against its receptor coordinates."""
    import gemmi
    from rdkit import Chem

    structure = gemmi.read_structure(str(receptor_path))
    molecule = Chem.MolFromMolFile(str(ligand_path), removeHs=False)
    if not structure or molecule is None or molecule.GetNumConformers() != 1:
        raise ValueError("The selected protein and ligand coordinates could not be read.")
    atoms, xyz, protein, ligand, waters = [], [], [], [], []
    donors, acceptors, charges = {}, set(), {}
    for chain in structure[0]:
        for residue in chain:
            name = residue.name.upper()
            is_water = residue.is_water()
            if not is_water and not gemmi.find_tabulated_residue(name).is_amino_acid():
                continue
            label = f"{chain.name.strip()}:{residue.seqid}:{residue.name}"
            positions = {}
            hydrogens = []
            for atom in residue.first_conformer():
                index = len(atoms)
                atoms.append(Atom(atom.name.strip(), atom.element.name.upper(), label, name))
                xyz.append((atom.pos.x, atom.pos.y, atom.pos.z))
                if atom.element.is_hydrogen:
                    hydrogens.append(index)
                    continue
                positions[index] = atom.name.strip().upper()
                if not is_water:
                    protein.append(index)
                elif atom.element.name == "O":
                    waters.append(index)
            residue_donors = {}
            for index, atom_name in positions.items():
                if (is_water and index in waters or name != "PRO" and atom_name == "N"
                        or atom_name in PROTEIN_DONORS.get(name, "").split()):
                    residue_donors[index] = []
                if is_water and index in waters or atom_name in {"O", "OXT"} or atom_name in PROTEIN_ACCEPTORS.get(name, "").split():
                    acceptors.add(index)
            for hydrogen in hydrogens:
                distances = [(float(np.linalg.norm(np.subtract(xyz[hydrogen], xyz[index]))), index)
                             for index in positions]
                if distances:
                    distance, index = min(distances)
                    if distance <= 1.25 and index in residue_donors:
                        residue_donors[index].append(hydrogen)
            donors.update(residue_donors)

    offset = len(atoms)
    conformer = molecule.GetConformer()
    for atom in molecule.GetAtoms():
        index = offset + atom.GetIdx()
        name = atom.GetProp("studio_atom_name") if atom.HasProp("studio_atom_name") else f"{atom.GetSymbol()}{atom.GetIdx() + 1}"
        atoms.append(Atom(name, atom.GetSymbol().upper(), "", ""))
        xyz.append(tuple(conformer.GetAtomPosition(atom.GetIdx())))
        if atom.GetAtomicNum() != 1:
            ligand.append(index)
            charges[index] = atom.GetFormalCharge()
    for feature in _rdkit_features().GetFeaturesForMol(molecule):
        for index in feature.GetAtomIds():
            if feature.GetFamily() == "Donor":
                donors[offset + index] = [offset + atom.GetIdx() for atom in molecule.GetAtomWithIdx(index).GetNeighbors()
                                          if atom.GetAtomicNum() == 1]
            elif feature.GetFamily() == "Acceptor":
                acceptors.add(offset + index)
    if not protein or not ligand:
        raise ValueError("A protein and a small-molecule ligand are both required.")
    classifier = Classifier(atoms, protein, ligand, waters, donors, acceptors, charges, pose=True)
    classifier.charge_source = "formal charges from the selected ligand chemical graph"
    occupancy = Occupancy(classifier)
    counts = occupancy.update(classifier.frame(np.asarray(xyz)))
    return occupancy, counts
