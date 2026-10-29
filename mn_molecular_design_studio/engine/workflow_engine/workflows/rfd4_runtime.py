"""Studio's task adapter, executed inside the pinned RFD4 container.

Uses upstream AtomWorks condition APIs and RFD4's own input validator. No model
implementation or scientific pass criteria are reimplemented here.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import yaml
import rfproteina.conditions  # register before reading annotated structures
from atomworks.enums import ChainType
from atomworks.io import parse
from atomworks.io.transforms.atom_array import add_id_and_entity_annotations, add_iid_annotations
from atomworks.io.utils.annotator import ensure_annotations
from atomworks.io.utils.atom_array_plus import as_atom_array_plus, concatenate_atom_array_plus
from atomworks.io.utils.io_utils import to_cif_file
from atomworks.io.utils.standard_annotations import S_SEGMAX, S_SEGMIN
from atomworks.ml.conditions import C_CRD, C_IDX, C_SEQ
from atomworks.ml.utils.expandable_segment import insert_expandable_segment
from biotite.structure import BondType, concatenate, get_residue_starts
from biotite.structure.info import residue as ccd_residue
from rfproteina.conditions import C_CRDGRP
from rfproteina.conditions_reporting import (
    check_index_group_molecule_consistency, check_index_group_res_id_gaps,
    run_pipeline_validation, summarize_conditions, sync_index_groups,
)
from rfproteina.inference_utils.condition_spec import build_annotated_atom_array_from_condition_spec, split_global_args
from rfproteina.inference_utils.nucleic_ss_input import parse_dot_bracket
from rfproteina.transforms.utils.conditioning_utils import (
    convert_to_internal_coords, materialize_internal_coords, separate_coordinate_groups,
    set_absolute_coordinate_conditioning,
)
from rfproteina.utils.condition_resolution import resolve_standard_annotation, set_cath_condition
from rfproteina.utils.parse_presets import ANNOTATED_CIF, RCSB_WITH_METALS

BACKBONE = ['N', 'CA', 'C', 'O']
TIPS = {
    'ARG': ['NE', 'CZ', 'NH1', 'NH2'], 'ASN': ['CG', 'OD1', 'ND2'],
    'ASP': ['CG', 'OD1', 'OD2'], 'CYS': ['SG'], 'GLN': ['CD', 'OE1', 'NE2'],
    'GLU': ['CD', 'OE1', 'OE2'], 'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2'],
    'LYS': ['NZ'], 'MET': ['SD', 'CE'], 'SER': ['OG'], 'THR': ['OG1'],
    'TYR': ['CZ', 'OH'], 'TRP': ['CD1', 'NE1', 'CE2'],
    'PHE': ['CD1', 'CD2', 'CE1', 'CE2', 'CZ'], 'ILE': ['CD1'],
    'LEU': ['CD1', 'CD2'], 'VAL': ['CG1', 'CG2'],
}
AA3 = dict(zip('ACDEFGHIKLMNPQRSTVWY', [
    'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU',
    'MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']))


def load(path, annotated=False):
    config=(ANNOTATED_CIF if annotated else RCSB_WITH_METALS).replace(build_assembly=None)
    aa = as_atom_array_plus(parse(str(path), config=config)['asym_unit'][0])
    if 'transformation_id' not in aa.get_annotation_categories():
        aa.set_annotation('transformation_id',np.full(aa.array_length(),'1',dtype='U8'))
    aa=add_iid_annotations(aa)
    ensure_annotations(aa, 'is_protein', 'is_ligand', 'is_nucleic_acid')
    return aa


def selected(aa, text, *, required=False):
    mask = np.zeros(aa.array_length(), dtype=bool)
    for token in re.split(r'[,;\s]+', str(text or '').strip()):
        if not token:
            continue
        m = re.fullmatch(r'([A-Za-z_0-9]+):(-?\d+)(?:-(-?\d+))?', token) or re.fullmatch(r'([A-Za-z_]+)(-?\d+)(?:-(-?\d+))?', token)
        if not m:
            raise ValueError(f'Invalid residue selection {token!r}; use A25 or A25-35.')
        chain, start, end = m[1], int(m[2]), int(m[3] or m[2])
        if end < start:
            raise ValueError('Residue ranges must be ascending.')
        match = (aa.chain_id == chain) & (aa.res_id >= start) & (aa.res_id <= end)
        present = set(aa.res_id[match].tolist())
        if not match.any() or any(i not in present for i in range(start, end + 1)):
            raise ValueError(f'Residue selection {token!r} includes absent residues.')
        mask |= match
    if required and not mask.any():
        raise ValueError('Select at least one residue for this task.')
    return mask


def sequence_condition(aa, mask):
    values = C_SEQ.annotation(aa, default='generate').copy()
    values[mask] = C_SEQ.ground_truth_value(aa)[mask]
    C_SEQ.set_annotation(aa, values)


def atom_mask(aa, region, kind, atom=None):
    region=region & np.isfinite(aa.coord).all(axis=1)
    if kind in ('whole', 'all'):
        return region.copy()
    if kind == 'backbone':
        return region & np.isin(aa.atom_name, BACKBONE)
    if kind == 'single':
        if atom is None:
            chosen=np.zeros_like(region)
            for chain,start,end in residue_runs(aa,region):
                for number in range(start,end+1):
                    residue=region & (aa.chain_id==chain) & (aa.res_id==number)
                    ca=np.flatnonzero(residue & (aa.atom_name=='CA') & np.isfinite(aa.coord).all(axis=1))
                    if not len(ca):raise ValueError(f'{chain}:{number} has no resolved CA for hotspot atom selection.')
                    candidates=np.flatnonzero(residue & ~np.isin(aa.atom_name,BACKBONE+['OXT']) & (aa.element!='H') & np.isfinite(aa.coord).all(axis=1))
                    index=candidates[np.argmax(np.linalg.norm(aa.coord[candidates]-aa.coord[ca[0]],axis=1))] if len(candidates) else ca[0]
                    chosen[index]=True
            return chosen
        chosen = region & (aa.atom_name == atom)
        if len(set(zip(aa.chain_id[chosen], aa.res_id[chosen]))) != len(set(zip(aa.chain_id[region], aa.res_id[region]))):
            raise ValueError(f'Atom {atom} is not resolved in every selected residue.')
        return chosen
    if kind == 'sidechain':
        return region & (~np.isin(aa.atom_name, BACKBONE + ['OXT']) | ((aa.res_name == 'GLY') & (aa.atom_name == 'CA')))
    if kind == 'sidechain_beyond_cb':
        return region & (~np.isin(aa.atom_name, BACKBONE + ['CB', 'OXT']) | (aa.atom_name=='CA'))
    if kind == 'hotspot_tip':
        return region & (~np.isin(aa.atom_name,BACKBONE+['CB','OXT']) | ((aa.res_name=='GLY') & (aa.atom_name=='CA')) | ((aa.res_name=='ALA') & (aa.atom_name=='CB')))
    if kind == 'tip':
        out = np.zeros_like(region)
        for name in np.unique(aa.res_name[region]):
            out |= region & (aa.res_name == name) & np.isin(aa.atom_name, TIPS.get(str(name), BACKBONE))
        return out
    raise ValueError(f'Unknown atom selection: {kind}')


def add_chain(aa, low, high):
    used = set(aa.chain_id)
    chain = next(c for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz' if c not in used)
    return insert_expandable_segment(aa, seg_min=low, seg_max=high, new_chain=True,
                                    new_chain_id=chain, new_chain_type=ChainType.POLYPEPTIDE_L), chain


def refresh(aa):
    ensure_annotations(aa, 'is_polymer')
    if 'transformation_id' not in aa.get_annotation_categories():
        aa.set_annotation('transformation_id',np.full(aa.array_length(),'1',dtype='U8'))
    aa = add_id_and_entity_annotations(aa, overwrite=True)
    aa = add_iid_annotations(aa, overwrite=True)
    for name in ['is_chain_of_interest']:
        if name in aa.get_annotation_categories():
            aa.del_annotation(name)
    ensure_annotations(aa, 'is_protein', 'is_ligand', 'is_nucleic_acid')
    return as_atom_array_plus(aa)


def condition_values(aa, name, mask, value, *, combine=False):
    cls=resolve_standard_annotation(name)
    values=cls.annotation(aa,default='generate').copy()
    values[mask]=values[mask] | value if combine else value
    cls.set_annotation(aa,values)


def ligand_atoms(aa, text, *, metal_default=False):
    mask=np.zeros(aa.array_length(),dtype=bool)
    if not text.strip() and metal_default:
        mask=aa.is_ligand & np.isin(aa.element,['MG','ZN','FE','MN','CA','CO','NI','CU','CD'])
    for token in re.split(r'[,;\s]+',text.strip()):
        if not token:continue
        residue,sep,atom=token.partition(':')
        hit=aa.is_ligand & (aa.atom_name==(atom if sep else residue))
        if sep:hit &= aa.res_name==residue
        if not hit.any():raise ValueError(f'Ligand atom {token!r} is absent from the retained ligands.')
        if not sep and len(set(zip(aa.chain_id[hit],aa.res_id[hit])))>1:
            raise ValueError(f'Ligand atom {token!r} is ambiguous; prefix its ligand code, e.g. SAM:{token}.')
        mask |= hit
    if not mask.any():raise ValueError('Select resolved ligand atoms for this interaction.')
    return mask


def ligand_conditions(aa,o):
    if o.get('ligand_atoms_to_hold'):
        set_absolute_coordinate_conditioning(aa,np.flatnonzero(ligand_atoms(aa,o['ligand_atoms_to_hold'])),True)
    for key,value in [('ligand_hbond_donors',1),('ligand_hbond_acceptors',2)]:
        if o.get(key):condition_values(aa,'C_HB',ligand_atoms(aa,o[key]),value,combine=True)
    burial={'model':0,'buried':1,'intermediate':2,'exposed':3}[o.get('ligand_burial','model')]
    if burial:condition_values(aa,'C_RAS',aa.is_ligand,burial)
    if o.get('exposed_ligand_atoms'):condition_values(aa,'C_RAS',ligand_atoms(aa,o['exposed_ligand_atoms']),3)
    interactions={
        'Hydrogen-bond donor':('C_HB',1),'Hydrogen-bond acceptor':('C_HB',2),
        'Hydrogen-bond donor and acceptor':('C_HB',3),'Pi-stacking ring':('C_PISTACK',True),
        'Cation-pi ring':('C_PICATION',1),'Cation-pi cation':('C_PICATION',2),
        'Salt bridge, positive':('C_SALT',1),'Salt bridge, negative':('C_SALT',2),
        'Hydrophobic contact':('C_HPHOB',True),'Metal ion to coordinate':('C_MCOORD',1)}
    for row in o.get('ligand_interactions') or []:
        name,value=interactions[row['interaction']]
        mask=ligand_atoms(aa,row.get('ligand_atoms',''),metal_default=name=='C_MCOORD')
        condition_values(aa,name,mask,value,combine=True)


def residue_runs(aa,mask):
    rows=[]
    for idx in get_residue_starts(aa):
        if not mask[idx]:continue
        chain,number=str(aa.chain_id[idx]),int(aa.res_id[idx])
        if rows and rows[-1][0]==chain and rows[-1][2]+1==number:rows[-1][2]=number
        else:rows.append([chain,number,number])
    return rows


def resize_cdrs(aa,designed,o):
    ranges=o.get('cdr_lengths') or []
    if not ranges:return aa
    sentinel='annotation_expsegmin_1_residue'
    if sentinel in aa.get_annotation_categories():
        indices=[i for i in get_residue_starts(aa) if aa.chain_id[i] in designed and aa.get_annotation(sentinel)[i]>0]
        if indices:
            if len(indices)!=len(ranges):raise ValueError(f'Provide {len(indices)} CDR length ranges, in chain and loop order.')
            for idx,row in zip(indices,ranges):
                segment=(aa.chain_id==aa.chain_id[idx]) & (aa.res_id==aa.res_id[idx])
                for name,value in [('annotation_expsegmin_1_residue',row['lower_bound']),('annotation_expsegmax_1_residue',row['upper_bound'])]:
                    aa.get_annotation(name)[segment]=value
                aa.res_name[segment]=f"{row['lower_bound']}-{row['upper_bound']}"
            return aa
    cdr=selected(aa,o.get('cdr_residues',''),required=True)
    runs=residue_runs(aa,cdr)
    if len(runs)!=len(ranges):raise ValueError(f'Provide {len(runs)} CDR length ranges, in chain and loop order.')
    # Native contigs rebuild polymer bonds and residue numbering after loop insertion.
    from rfproteina.inference_utils.contig_spec import apply_contig
    tokens=[]; ordered=list(dict.fromkeys(aa.chain_id.tolist())); overrides={tuple(run):row for run,row in zip(runs,ranges)}
    for chain in ordered:
        if tokens:tokens.append('/0')
        indices=get_residue_starts(aa[aa.chain_id==chain]); ids=aa.res_id[aa.chain_id==chain][indices]
        groups=[run for run in runs if run[0]==chain]; cursor=int(ids.min())
        for run in groups:
            if cursor<run[1]:tokens.append(f'{chain}{cursor}-{run[1]-1}')
            row=overrides[tuple(run)];tokens.append(f"{row['lower_bound']}-{row['upper_bound']}");cursor=run[2]+1
        if cursor<=ids.max():tokens.append(f'{chain}{cursor}-{ids.max()}')
    return refresh(apply_contig(aa,','.join(tokens))[0])


def trim_constant_domains(aa,designed,task):
    from anarcii import Anarcii
    from biotite.sequence import ProteinSequence
    numberer=Anarcii(seq_type='tcr' if task=='tcr' else 'antibody',cpu=True,ncpu=1,verbose=False)
    keep=~np.isin(aa.chain_id,designed)
    for chain in designed:
        indices=np.flatnonzero(aa.chain_id==chain);part=aa[indices]
        starts=get_residue_starts(part,add_exclusive_stop=True)
        seq=''.join(ProteinSequence.convert_letter_3to1(str(n)) for n in part.res_name[starts[:-1]])
        records=numberer.number({'query':seq})
        found=False
        for record in records.values():
            if record.get('error') or not record.get('numbering'):continue
            letters=''.join(aa for _,aa in record['numbering'] if aa!='-'); start=record.get('query_start')
            if start is None or seq[start:start+len(letters)]!=letters:continue
            keep[indices[starts[start]:starts[start+len(letters)]]]=True;found=True
        if not found:raise ValueError(f'Could not identify a variable domain in chain {chain}. Check the selected chain or enable Keep Constant Domains.')
    return refresh(aa[keep].copy())


def rna_input(dbn):
    # Complete CCD nucleotides retain chemistry. Their sequence and coordinates
    # are unconditioned; only polymer type and authored secondary structure remain.
    base = ccd_residue('A')
    base = base[(base.element != 'H') & (base.atom_name != 'OP3')]
    direction=base.coord[base.atom_name=="O3'"][0]-base.coord[base.atom_name=='P'][0]
    translation=direction+1.6*direction/np.linalg.norm(direction)
    pieces = []
    for i in range(len(dbn)):
        part = base.copy(); part.chain_id[:] = 'A'; part.res_id[:] = i + 1
        part.coord += i*translation; part.hetero[:] = False
        pieces.append(part)
    aa = as_atom_array_plus(concatenate(pieces))
    aa.set_annotation('chain_type', np.full(aa.array_length(), ChainType.RNA, dtype=np.int16))
    aa.set_annotation('occupancy',np.ones(aa.array_length(),dtype=float))
    aa.set_annotation('b_factor',np.zeros(aa.array_length(),dtype=float))
    for i in range(len(dbn) - 1):
        left = np.flatnonzero((aa.res_id == i + 1) & (aa.atom_name == "O3'"))[0]
        right = np.flatnonzero((aa.res_id == i + 2) & (aa.atom_name == 'P'))[0]
        aa.bonds.add_bond(int(left), int(right), BondType.SINGLE)
    # Run the normal chemistry parser to derive atomic numbers, aromaticity,
    # polymer entities and the other model features from the CCD nucleotides.
    from biotite.structure.io.pdbx import CIFFile, set_structure
    with tempfile.TemporaryDirectory(prefix='rna-input-') as temporary:
        cif=CIFFile(); set_structure(cif,aa)
        path=Path(temporary)/'rna.cif'; cif.write(str(path))
        aa=load(path)
    pairs, loops = parse_dot_bracket(dbn)
    states = np.zeros(len(dbn), dtype=np.int16)
    states[loops] = 2
    for a, b in pairs:
        states[a] = states[b] = 1
    from rfproteina.transforms.utils.na_geom_utils import _set_canonical_na_ss_from_token_state
    _set_canonical_na_ss_from_token_state(aa, token_starts=get_residue_starts(aa), state_values=states,
                                       pair_inds=np.asarray(pairs, dtype=np.int32).reshape(-1, 2))
    return aa


def nanobody_framework(name, target):
    slug = {'Vobarilizumab':'vobarilizumab','Caplacizumab':'caplacizumab',
            'Ozoralizumab':'ozoralizumab','h-NbBCII10':'h-NbBCII10'}[name]
    path = next(Path('/opt/rfd4/examples/antibody_design/structures/vhh').glob(f'pmvk*{slug}_cdr123.cif.zst'))
    template = load(path, True)
    chains = list(np.unique(template.chain_id))
    # The target is fully sequence-conditioned; the framework's CDRs are not.
    chain = min(chains, key=lambda c: C_SEQ.mask(template, default='generate')[template.chain_id == c].mean())
    framework = template[template.chain_id == chain].copy()
    replacement = next(c for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' if c not in set(target.chain_id))
    framework.chain_id[:] = replacement
    arrays = separate_coordinate_groups([target, framework])
    columns = {name for a in arrays for name in a.get_annotation_categories()
               if name.startswith(('condition_', 'mask_'))}
    for a in arrays:
        ensure_annotations(a, *columns)
    # The uploaded target and framework CIFs may have different metadata
    # coverage. Missing expandable-segment annotations must use AtomWorks'
    # native sentinels: zero is a valid segment length and would make every
    # target residue look like a zero-length segment.
    for cls in (S_SEGMIN, S_SEGMAX):
        for a in arrays:
            cls.set_annotation(a, cls.annotation(a, default='generate'))
    # Parse metadata differs between an uploaded target and the bundled VHH.
    # Retain every condition/expandable-segment field, filling absent metadata
    # with dtype-correct neutral values before IDs are derived again.
    for name in {n for a in arrays for n in a.get_annotation_categories()}:
        exemplar=next(a.get_annotation(name) for a in arrays if name in a.get_annotation_categories())
        for a in arrays:
            if name not in a.get_annotation_categories():
                a.set_annotation(name,np.zeros(a.array_length(),dtype=exemplar.dtype))
    return refresh(concatenate_atom_array_plus(arrays)), [replacement]


def guided(request):
    p, o = request['parameters'], request['parameters'].get('options') or {}
    task = p['task']; low, high = p.get('length_min', 60), p.get('length_max', 120)
    if task == 'unconditional':
        aa = load('/opt/rfd4/data/benchmarks/monomer/monomer-short.cif.zst', True)
        aa.set_annotation('res_name', np.full(aa.array_length(), f'{low}-{high}', dtype='U32'))
        aa.set_annotation('annotation_expsegmin_1_residue',np.full(aa.array_length(),low,dtype=np.int32))
        aa.set_annotation('annotation_expsegmax_1_residue',np.full(aa.array_length(),high,dtype=np.int32))
        return aa, list(np.unique(aa.chain_id))
    if task == 'rna_design':
        return rna_input(o['rna_secondary_structure']), ['A']
    # Post-hoc scoring must retain the design's StandardAnnotations; the ordinary
    # RCSB preset intentionally drops them and would make motif/hotspot metrics
    # and chain-of-interest detection incomplete.
    aa = load(request['target_structure'], annotated=task == 'score_structures')
    if task in ('enzyme','small_molecule_binder'):
        additives=['GOL','SO4','PEG','PGE','EDO','MPD']
        keep=~(aa.is_ligand & np.isin(aa.res_name,additives))
        if o.get('ligands_to_keep'):
            names=[x.upper() for x in re.split(r'[,;\s]+',o['ligands_to_keep']) if x]
            absent=set(names)-set(aa.res_name[aa.is_ligand])
            if absent:raise ValueError('Requested ligand codes are absent: '+', '.join(sorted(absent)))
            keep &= ~aa.is_ligand | np.isin(aa.res_name,names)
        aa=refresh(aa[keep].copy())
    target_ids = p.get('target_chains') or []
    target_mask = np.isin(aa.chain_id, target_ids) & aa.is_protein
    hotspots = selected(aa, p.get('hotspots', ''))
    if (hotspots & ~target_mask).any():
        raise ValueError('Hotspots must lie on selected target chains.')
    if hotspots.any():
        kind=o.get('hotspot_atoms','whole')
        condition_values(aa,'C_HOT',atom_mask(aa,hotspots,'hotspot_tip' if kind=='tip' else kind,o.get('hotspot_atom')),True)
    designed = []
    if task in ('protein_binder','small_molecule_binder','dna_binder'):
        if task == 'protein_binder' and not target_mask.any():
            raise ValueError('Select at least one target protein chain.')
        if task == 'small_molecule_binder' and not aa.is_ligand.any():
            raise ValueError('The structure contains no ligand; upload a mmCIF with complete ligand chemistry.')
        if task == 'dna_binder' and not np.isin(aa.res_name, ['DA','DC','DG','DT']).any():
            raise ValueError('The structure contains no DNA.')
        keep = target_mask.copy()
        if task == 'small_molecule_binder': keep |= aa.is_ligand
        if task == 'dna_binder': keep |= aa.is_nucleic_acid
        aa = aa[keep].copy(); ensure_annotations(aa, 'is_protein','is_ligand','is_nucleic_acid')
        fixed = aa.is_protein
        sequence_condition(aa, np.ones(aa.array_length(), dtype=bool))
        convert_to_internal_coords(aa, atom_mask(aa, fixed, 'backbone' if o.get('repack_target', True) else 'all'))
        if task == 'small_molecule_binder' and o.get('ligand_pose') == 'fixed':
            set_absolute_coordinate_conditioning(aa, np.flatnonzero(aa.is_ligand), True)
        if task == 'dna_binder' and o.get('dna_conditioning') == 'structure':
            set_absolute_coordinate_conditioning(aa, np.flatnonzero(aa.is_nucleic_acid), True)
        aa, chain = add_chain(aa, low, high); designed = [chain]
    elif task in ('nanobody','antibody','tcr'):
        if not target_mask.any(): raise ValueError('Select target protein chains.')
        if task == 'nanobody' and o.get('nanobody_framework', 'Vobarilizumab') != 'Your own nanobody':
            aa = aa[target_mask].copy()
            sequence_condition(aa, np.ones(aa.array_length(), dtype=bool))
            convert_to_internal_coords(aa, atom_mask(aa, aa.is_protein, 'backbone'))
            aa, designed = nanobody_framework(o.get('nanobody_framework', 'Vobarilizumab'), aa)
            framework=np.isin(aa.chain_id,designed) & C_SEQ.mask(aa,default='generate')
            coords=C_CRD.annotation(aa,default='generate').copy();coords[np.isin(aa.chain_id,designed)]=np.nan;C_CRD.set_annotation(aa,coords)
            groups=C_CRDGRP.annotation(aa,default='generate').copy()
            groups[np.isin(aa.chain_id,designed)]=C_CRDGRP.default_annotation(aa)[np.isin(aa.chain_id,designed)]
            C_CRDGRP.set_annotation(aa,groups,on_internal_coord_mismatch='clear')
            convert_to_internal_coords(aa,atom_mask(aa,framework,o.get('framework_atoms','backbone')))
        else:
            designed = o.get('design_chains') or []
            if not designed or set(designed) & set(target_ids):
                raise ValueError('Select separate target and design chains.')
            aa = aa[np.isin(aa.chain_id, target_ids + designed)].copy()
            if not o.get('keep_constant_domains',True):aa=trim_constant_domains(aa,designed,task)
            cdr = selected(aa, o.get('cdr_residues', ''), required=True)
            if (cdr & ~np.isin(aa.chain_id, designed)).any(): raise ValueError('CDRs must lie on design chains.')
            framework = np.isin(aa.chain_id, designed) & ~cdr
            target = np.isin(aa.chain_id, target_ids)
            sequence_condition(aa, target | framework)
            convert_to_internal_coords(aa, atom_mask(aa, target, 'backbone' if o.get('repack_target',True) else 'all'))
            if o.get('dock_independently'):
                for chain in designed:convert_to_internal_coords(aa,atom_mask(aa,framework & (aa.chain_id==chain),o.get('framework_atoms','backbone')))
            else:convert_to_internal_coords(aa,atom_mask(aa,framework,o.get('framework_atoms','backbone')))
        aa=resize_cdrs(aa,designed,o)
    elif task in ('motif_scaffolding','enzyme'):
        motif = selected(aa, o.get('motif_residues', ''), required=True)
        ligand = aa.is_ligand if task == 'enzyme' else np.zeros_like(motif)
        aa = aa[motif | ligand].copy(); ensure_annotations(aa, 'is_protein','is_ligand','is_nucleic_acid')
        motif = aa.is_protein
        sequence_condition(aa, np.ones(aa.array_length(), dtype=bool))
        fixed = atom_mask(aa, motif, o.get('motif_atoms', 'tip' if task == 'enzyme' else 'all'))
        for token in o.get('active_site_atoms','').split():
            residue,_,names=token.rpartition(':'); region=selected(aa,residue,required=True)
            names=names.split(',');chosen=region & np.isin(aa.atom_name,names)
            missing=set(names)-set(aa.atom_name[chosen])
            if missing:raise ValueError(f'{residue}: active-site atoms are absent: {sorted(missing)}')
            if (region & ~motif).any():raise ValueError('Active-site atom overrides must refer to selected motif residues.')
            fixed[region]=chosen[region]
        set_absolute_coordinate_conditioning(aa, np.flatnonzero(fixed), True)
        if o.get('ligand_pose') == 'fixed':
            set_absolute_coordinate_conditioning(aa, np.flatnonzero(aa.is_ligand), True)
        if o.get('unindex_motif', task == 'enzyme'):
            ranges=o.get('scaffold_lengths') or []
            if ranges:
                if o.get('scaffold_total_explicit'):raise ValueError('Use either per-gap Scaffold Lengths or total length bounds, not both.')
                runs=residue_runs(aa,motif)
                if len(ranges)!=len(runs)+1:raise ValueError(f'Provide {len(runs)+1} scaffold ranges, before, between and after the motif blocks.')
                if len({run[0] for run in runs})!=1:raise ValueError('Per-gap scaffold lengths need motif blocks from a single chain.')
                from rfproteina.inference_utils.contig_spec import apply_contig
                tokens=[]
                for i,(chain,start,end) in enumerate(runs):
                    gap=ranges[i];tokens += [f"{gap['lower_bound']}-{gap['upper_bound']}",f'{chain}{start}-{end}']
                gap=ranges[-1];tokens.append(f"{gap['lower_bound']}-{gap['upper_bound']}")
                for chain,start,end in residue_runs(aa,aa.is_ligand):tokens+=['/0',f'{chain}{start}-{end}']
                motif_count=len(get_residue_starts(aa[motif]))
                low=motif_count+sum(row['lower_bound'] for row in ranges);high=motif_count+sum(row['upper_bound'] for row in ranges)
                aa=refresh(apply_contig(aa,','.join(tokens))[0]);motif=aa.is_protein
            indexed = C_IDX.mask(aa, default='generate').copy(); indexed[motif] = False
            C_IDX.set_annotation(aa, indexed)
            sync_index_groups(aa, preserve_manual_edits=False)
            aa, chain = add_chain(aa, low, high); designed = [chain]
        else:
            protein_ids = list(np.unique(aa.chain_id[motif]))
            if len(protein_ids) != 1: raise ValueError('An indexed motif must be a single chain; enable motif placement for multiple chains.')
            runs=residue_runs(aa,motif);ranges=o.get('scaffold_lengths') or []
            if task=='enzyme' and o.get('scaffold_total_explicit'):
                raise ValueError('Total scaffold length bounds require Let The Model Place The Motif. For an indexed motif, set per-gap Scaffold Lengths.')
            if not ranges:
                count=len(get_residue_starts(aa[motif]));gaps=len(runs)+1
                if low<=count:raise ValueError('Scaffold length must exceed the number of motif residues.')
                # Preserve all motif blocks with gaps before, between and after them.
                mins=[(low-count)//gaps]*gaps;maxs=[(high-count)//gaps]*gaps
                for i in range((low-count)%gaps):mins[i]+=1
                for i in range((high-count)%gaps):maxs[i]+=1
                ranges=[{'lower_bound':a,'upper_bound':b} for a,b in zip(mins,maxs)]
            if len(ranges)!=len(runs)+1:raise ValueError(f'Provide {len(runs)+1} scaffold ranges, before, between and after the motif blocks.')
            from rfproteina.inference_utils.contig_spec import apply_contig
            tokens=[]
            for i,(chain,start,end) in enumerate(runs):
                gap=ranges[i]
                if gap['upper_bound']:tokens.append(f"{gap['lower_bound']}-{gap['upper_bound']}")
                tokens.append(f'{chain}{start}-{end}')
            gap=ranges[-1]
            if gap['upper_bound']:tokens.append(f"{gap['lower_bound']}-{gap['upper_bound']}")
            for chain,start,end in residue_runs(aa,aa.is_ligand):tokens+=['/0',f'{chain}{start}-{end}']
            aa=refresh(apply_contig(aa,','.join(tokens))[0])
            designed = protein_ids
    elif task in ('existing_structure','score_structures'):
        designed = o.get('design_chains') or [c for c in np.unique(aa.chain_id[aa.is_protein]) if c not in o.get('fixed_chains', [])]
        if task == 'score_structures': designed = [o['designed_chain']] if o.get('designed_chain') else designed
        if not designed: raise ValueError('Select a designed chain.')
        absent=set(designed)-set(aa.chain_id[aa.is_protein])
        if absent:raise ValueError('Designed protein chains are absent from this structure: '+', '.join(sorted(absent)))
        change = np.isin(aa.chain_id, designed)
        kept=selected(aa,o.get('keep_residues',''))
        fixed = ~change | kept
        if o.get('redesign_shell') is not None:
            from scipy.spatial.distance import cdist
            anchors=fixed & ((aa.element!='H'))
            if not anchors.any():raise ValueError('Select fixed chains or kept residues before setting a redesign shell.')
            # Bound memory while computing nearest distance to the fixed region.
            distances=np.concatenate([cdist(aa.coord[i:i+512],aa.coord[anchors]).min(axis=1) for i in range(0,len(aa),512)])
            for start,end in zip(get_residue_starts(aa,add_exclusive_stop=True)[:-1],get_residue_starts(aa,add_exclusive_stop=True)[1:]):
                if distances[start:end].min()>o['redesign_shell']:fixed[start:end]=True
        sequence_condition(aa, fixed)
        coordinate_fixed=~change | (fixed if o.get('hold_kept_in_place',True) else np.zeros_like(fixed))
        set_absolute_coordinate_conditioning(aa, np.flatnonzero(coordinate_fixed), True)
        operation = o.get('operation','diversify')
        if task == 'score_structures':
            operation = 'structure_prediction'
            sequence_condition(aa, np.ones(aa.array_length(),dtype=bool))
        if operation in ('inverse_folding','sidechain_packing'):
            set_absolute_coordinate_conditioning(aa, np.flatnonzero(change & np.isin(aa.atom_name,BACKBONE)), True)
        if operation in ('structure_prediction','sidechain_packing'):
            sequence_condition(aa, np.ones(aa.array_length(),dtype=bool))
        if operation=='sidechain_packing' and o.get('ligand_pose')=='fixed':set_absolute_coordinate_conditioning(aa,np.flatnonzero(aa.is_ligand),True)
        occupied=np.zeros(aa.array_length(),dtype=bool)
        for key,value in [('helix_residues',1),('strand_residues',2),('loop_residues',3)]:
            region=selected(aa,o.get(key,''))
            if (region & occupied).any():raise ValueError('Helix, strand and loop selections must not overlap.')
            if (region & ~aa.is_protein).any():raise ValueError('Secondary structure selections must refer to protein residues.')
            condition_values(aa,'C_SEC',region,value);occupied |= region
    else:
        raise ValueError(f'Task {task!r} requires a native input.')
    aa = refresh(aa)
    if task in ('enzyme','small_molecule_binder'):ligand_conditions(aa,o)
    return aa, designed


def common_conditions(aa, designed, options):
    region = np.isin(aa.chain_id, designed)
    # Framework and target identities remain fixed when applying sequence biases.
    variable = region & ~C_SEQ.mask(aa, default='generate')
    if options.get('cath'):
        set_cath_condition(aa, np.flatnonzero(region), options['cath'])
    probabilities = {letter:0. for letter in options.get('exclude_amino_acids', [])}
    probabilities.update(options.get('favor_amino_acids') or {})
    for letter, probability in probabilities.items():
        cls = resolve_standard_annotation('C_SEQ_PROB_' + AA3[letter])
        values = cls.annotation(aa, default='generate').copy()
        values[variable] = probability
        cls.set_annotation(aa, values)
    return aa


def prepare(request):
    p = request['parameters']; root = Path(request['root']); root.mkdir(parents=True, exist_ok=True)
    prepared = root / 'prepared'; prepared.mkdir(exist_ok=True)
    tasks = []; source_candidates={}
    mode = p.get('input_mode','guided')
    if mode in ('guided','guided_binder','unconditional'):
        if p['task'] in ('score_structures','existing_structure') and request.get('input_files'):
            for path in request['input_files']:
                child = {**request,'target_structure':path,'parameters':request.get('file_parameters',{}).get(path,p)}
                aa, designed = guided(child); tasks.append((Path(path).stem,aa,designed))
                source_candidates[Path(path).stem]=child['parameters'].get('_source_candidate')
        else:
            aa, designed = guided(request); tasks.append((p['task'], aa, designed))
    else:
        spec_path = Path(request['spec_path'])
        spec = split_global_args(yaml.safe_load(spec_path.read_text()))
        for name, entry in spec.items():
            aa = build_annotated_atom_array_from_condition_spec(entry, spec_dir=spec_path.parent, example_id=name)['atom_array']
            # The condition-spec builder preserves design annotations but does
            # not guarantee the convenience molecule-class annotations used by
            # the runtime inventory and downstream task checks.
            ensure_annotations(aa, 'is_protein', 'is_ligand', 'is_nucleic_acid')
            explicit = p.get('options',{}).get('design_chains')
            if explicit:
                designed = explicit
            else:
                from rfproteina.utils.chain import get_chain_of_interest
                chain, unambiguous = get_chain_of_interest(aa)
                if not unambiguous: raise ValueError('Native input has ambiguous designed chains. Set Design chains explicitly.')
                chain = str(chain)
                # AtomWorks' chain-of-interest helper may return the unique
                # chain-instance ID (for example ``B_1``), while the native
                # result inventory and chain_id selectors use ``B``. Resolve
                # that identity before applying task conditions or recording
                # target/design roles.
                chain_iids = aa.get_annotation('chain_iid') if 'chain_iid' in aa.get_annotation_categories() else None
                if chain_iids is not None and np.any(chain_iids == chain):
                    mapped = [str(value) for value in np.unique(aa.chain_id[chain_iids == chain])]
                    if len(mapped) != 1:
                        raise ValueError('Could not map the native chain of interest to one chain ID.')
                    chain = mapped[0]
                elif chain not in {str(value) for value in np.unique(aa.chain_id)}:
                    raise ValueError(f'Native chain of interest {chain!r} is not a chain ID in the assembled input.')
                designed = [chain]
            tasks.append((name,aa,designed))
    inventory = []
    for name, aa, designed in tasks:
        aa = common_conditions(aa, designed, p.get('options') or {})
        sync_index_groups(aa, preserve_manual_edits=True)
        materialize_internal_coords(aa)
        errors = [e for check in (run_pipeline_validation,check_index_group_molecule_consistency) for e in check(aa)]
        if errors: raise ValueError(f'{name}: ' + '\n'.join(map(str, errors)))
        path = prepared / f'{name}.cif'
        to_cif_file(aa, str(path), save_standard_annotations=True, chain_disambiguation='chain_iid')
        reloaded = load(path, True)
        # Round-trip the conditions consumed by the model.
        for cls in [C_SEQ,C_CRD,C_CRDGRP,C_IDX]+[resolve_standard_annotation(name) for name in ('C_HOT','C_HB','C_PISTACK','C_PICATION','C_SALT','C_HPHOB','C_MCOORD','C_RAS','C_SEC')]:
            left=cls.annotation(aa,default='generate'); right=cls.annotation(reloaded,default='generate')
            numeric = np.issubdtype(left.dtype, np.inexact)
            if left.shape != right.shape or not np.array_equal(left,right,equal_nan=numeric):
                raise ValueError(f'{name}: {cls.__name__} changed during CIF round-trip')
        inventory.append({'name':name,'path':str(path),'designed_chains':designed,'source_candidate':source_candidates.get(name),
                          'target_chains':[str(c) for c in np.unique(aa.chain_id) if c not in designed],
                          'ligand_codes':sorted({str(code) for code in aa.res_name[aa.is_ligand]}),
                          'atoms':aa.array_length(),'validation':'passed','warnings':check_index_group_res_id_gaps(aa),'conditions':str(summarize_conditions(aa))})
    (root/'prepared_inputs.json').write_text(json.dumps(inventory,indent=2))
    return inventory


def sampling_overrides(p):
    o=p.get('options') or {}; task=p['task']; effort={'Fast':'low','Balanced':'medium','Thorough':'high'}[o.get('sampling_effort','Balanced')]
    shape='midpeak' if task in ('enzyme','small_molecule_binder') else 'front_loaded'
    if task=='dna_binder' and effort=='high': shape='flat'
    if task=='existing_structure' and o.get('operation')=='inverse_folding': shape='flat'
    config=yaml.safe_load(Path(f'/opt/rfd4/configs/recycle_strategy/{effort}/{shape}.yaml').read_text())
    schedule=config['n_recycle_schedule']; steps=p.get('nsteps',60)
    # Upstream requires exactly nsteps schedule entries. Scale the same schedule
    # shape deterministically for user-selected step counts.
    widths=np.asarray([x[0] for x in schedule],dtype=float)*steps/60
    rounded=np.floor(widths).astype(int)
    for i in np.argsort(-(widths-rounded))[:steps-int(rounded.sum())]: rounded[i]+=1
    bands=[[int(n),int(pair[1])] for n,pair in zip(rounded,schedule) if n]
    overrides=[f'recycle_strategy={effort}/{shape}',f'generation.n_recycle_schedule={json.dumps(bands)}',f'generation.args.nsteps={steps}']
    if task=='existing_structure':
        if o.get('operation','diversify')=='diversify':
            if o.get('noise_mode','fraction')=='levels':overrides += [f'generation.args.partial_t_backbone={o.get("backbone_level",.8)}',f'generation.args.partial_t_latent={o.get("sequence_level",.8)}']
            else:overrides.append(f'generation.args.partial_trajectory_fraction_steps={o.get("noise_fraction",.3)}')
        elif o.get('operation')=='inverse_folding': overrides.append('generation.args.inverse_folding=true')
    for modality,prefix in [('bb_ca','backbone'),('local_latents','latent')]:
        overrides += [f'generation.model.{modality}.simulation_step_params.sc_scale_noise={o.get(prefix+"_noise",1.)}',
                      f'generation.model.{modality}.simulation_step_params.tsr_k={o.get(prefix+"_tsr_k",16. if prefix=="backbone" else 4.)}',
                      f'generation.model.{modality}.simulation_step_params.tsr_sigma={o.get(prefix+"_tsr_sigma",3. if prefix=="backbone" else 12.)}']
        if o.get('deterministic'): overrides.append(f'generation.model.{modality}.gt.mode=ode')
    if o.get('save_trajectory'): overrides += ['generation.args.save_trajectory_every=10','generation.args.bb_ca_only_trajectories=true']
    return overrides


def metric_overrides(p, *, no_ligand=False):
    o=p.get('options') or {}; task=p['task']
    if task=='score_structures': task=o.get('score_design_type','protein_binder')
    tags=['design']
    refold_tag={'protein_binder':'ppi_refolding','nanobody':'antibody','antibody':'antibody','tcr':'antibody',
                'enzyme':'ame_refolding','small_molecule_binder':'sm_binder_refolding','dna_binder':'dna_refolding',
                'motif_scaffolding':'motif_scaffold_refolding'}.get(task,'unconditional_refolding')
    if o.get('refold',True): tags.append(refold_tag)
    if p.get('hotspots'): tags.append('hotspot_conditioning')
    if o.get('unindex_motif'): tags.append('unindexed')
    if task in ('motif_scaffolding','enzyme'): tags.append('3d_coordinate_conditioning')
    if o.get('apo'): tags.append('apo_refolding')
    if task == 'enzyme' and no_ligand: tags.append('no_ligand')
    args=['metrics/refolding_oracles=of3',f'metrics_tags={json.dumps(tags)}',
          f'delete_failing_pre_oracle={str(o.get("filter_conditions",True)).lower()}']
    # Keep all attempts when requested, rather than merely retaining failed files
    # after the pre-oracle filter has already skipped their requested evaluation.
    if not o.get('filter_conditions',True): args.append('metrics.pre_oracle_filters={}')
    if o.get('developability'):
        args.append('++metrics.inline_metrics.sequence_developability='+json.dumps({'_target_':'rfproteina.metrics.physics.developability_metrics.SequenceDevelopabilityMetric'}))
    if o.get('energy'):
        args.append('++metrics.inline_metrics.interface_energy='+json.dumps({'_target_':'rfproteina.metrics.physics.tmol_block_pair_metric.TMolBlockPairMetric'}))
    for key,label in [('avoid_structures','avoid'),('also_bind_structures','also_bind')]:
        paths=p.get(key) or []
        if paths:
            tags.append(label+'_refolding')
            cfg={'_target_':'rfproteina.metrics.oracles.af3_oracle.AF3Oracle','oracle_name':'of3_'+label,
                 'template_query':'is_loaded_target & is_protein','aggregation':['median'],
                 'required_tags_any':[label+'_refolding'],'pre_oracle_transform':{
                     '_target_':'rfproteina.metrics.oracles.pre_oracle_transform.OffTargetSubstitutionTransform',
                     'off_target_paths':{f'{label}_{i+1}':path for i,path in enumerate(paths)}}}
            args.append('++metrics.refolding_oracles.oracle_'+label+'='+json.dumps(cfg))
    args[1]=f'metrics_tags={json.dumps(tags)}'
    # Authors' task-specific pass criteria remain the source of truth. A requested
    # override edits only the matching threshold entries in the composed config.
    return args


def search_overrides(p):
    o=p.get('options') or {}; search=o.get('reward_search','Off')
    if search=='Off': return []
    steps=p.get('nsteps',60)
    if search=='Interface contacts': provider='interchain_interactions'; key='interchain_interactions.hydrophobic_atom_pair_count_coi'; weight=1.
    elif search=='Low clashes': provider='clash_metrics'; key='clash_metrics.frac_atoms_clashing_1p5A'; weight=-1.
    else:
        provider='oracle_ppi' if p['task'] in ('protein_binder','nanobody','antibody','tcr') else 'oracle_dna' if p['task']=='dna_binder' else 'oracle'
        key=provider+'.plddt_coi'; weight=1.
    config={'beam_width':1,'n_branch':2,'checkpoints':{str(steps):{'metrics':[provider],'rewards':{key:{'weight':weight}}}}}
    return ['++search='+json.dumps(config)]


def compose_settings(command, p, root):
    from hydra import compose, initialize_config_dir
    from omegaconf import OmegaConf
    overrides=[]; updates=[]
    for arg in command[1:]:
        key,sep,value=arg.partition('=')
        if sep and value.startswith('{'):
            updates.append((key.lstrip('+'),json.loads(value)))
        else: overrides.append(arg)
    with initialize_config_dir(config_dir='/opt/rfd4/configs',version_base='1.3'):
        config=compose(config_name='evaluate' if p['task']=='score_structures' else 'inference',overrides=overrides)
    OmegaConf.set_struct(config,False)
    for key,value in updates: OmegaConf.update(config,key,value,merge=False,force_add=True)
    o=p.get('options') or {}
    for criteria in config.metrics.pass_criteria.values():
        for metric,limits in (criteria.get('metric_thresholds') or {}).items():
            if o.get('interface_error_limit') is not None and 'min_ipae' in metric:
                limits['max']=o['interface_error_limit']; limits.pop('max_exclusive',None)
            if o.get('binder_rmsd_limit') is not None and 'binder_rmsd' in metric:
                limits['max']=o['binder_rmsd_limit']
            if o.get('refold_backbone_rmsd_limit') is not None and metric=='refolded_backbone_rmsd.rmsd':limits['max']=o['refold_backbone_rmsd_limit']
            if o.get('motif_rmsd_limit') is not None and metric.startswith('refolded_motif_'):limits['max']=o['motif_rmsd_limit']
    path=root/'effective_config.yaml'
    OmegaConf.save(config,path,resolve=True)
    return [command[0],'--config-path',str(root),'--config-name','effective_config']


def _require_output_target_chains(task, origin, sequences, output_name):
    """Reject PPI results that silently lost one or more uploaded targets."""
    if task not in ('protein_binder','nanobody','antibody','tcr'):
        return
    expected={str(chain) for chain in origin.get('target_chains') or []}
    missing=expected-set(sequences)
    if missing:
        raise ValueError(
            f"{task} output {output_name} is missing expected target chain(s): "
            + ', '.join(sorted(missing))
        )


def run(request, prepare_only=False):
    p=request['parameters']; root=Path(request['root']); inventory=prepare(request)
    if prepare_only: return
    no_ligand=False
    if p['task']=='enzyme' or (p['task']=='score_structures' and p.get('options',{}).get('score_design_type')=='enzyme'):
        ligand_presence=[bool(item.get('ligand_codes')) for item in inventory]
        if len(set(ligand_presence))>1:
            raise ValueError('Score enzyme structures with and without ligands in separate runs so RFD4 applies the correct pass criteria.')
        no_ligand=bool(ligand_presence) and not ligand_presence[0]
    out=root/'output'; out.mkdir(exist_ok=True)
    if p['task']=='score_structures':
        command=['rfd4-evaluate',f'inputs={root / "prepared"}',f'out_dir={out}',*metric_overrides(p,no_ligand=no_ligand)]
    else:
        command=['rfd4-generate',f'inputs={root / "prepared"}',f'out_dir={out}',f'num_replicates={p.get("num_designs",8)}',
                 f'batch_size={p.get("batch_size",1)}',f'seed={p["seed"]}','prevalidate_inputs=true','skip_existing=false',
                 'num_workers=0','save_chain_of_interest=true',*sampling_overrides(p),*metric_overrides(p,no_ligand=no_ligand),*search_overrides(p)]
        if not p.get('compile_model'): command.append('compile_lengths=null')
    original_command=command
    count=1 if p['task']=='score_structures' else p.get('num_designs',8)
    commands=[]
    for start in range(0,count,100):
        current=original_command.copy(); chunk_root=root
        if count>100:
            chunk_root=root/'batches'/f'{start:06d}'; chunk_root.mkdir(parents=True,exist_ok=True)
            replacements={'out_dir':str(out/f'batch_{start:06d}'),
                          'num_replicates':str(min(100,count-start)),
                          'seed':str((p['seed']+start)%2**32)}
            current=[arg.split('=',1)[0]+'='+replacements[arg.split('=',1)[0]]
                     if '=' in arg and arg.split('=',1)[0] in replacements else arg for arg in current]
        command=compose_settings(current,p,chunk_root)
        commands.append({'command':command,'overrides':current,'status':'launch_window'})
        (root/'native_command.json').write_text(json.dumps(commands,indent=2))
        timing_file=os.environ.get('MN_STUDIO_RFD4_NATIVE_TIMING_FILE')
        if timing_file:
            from rfd4_native_diagnostics import phase, instrumented_command
            diagnostic_command=instrumented_command(command,Path(__file__).with_name('rfd4_native_diagnostics.py'))
            commands[-1]['diagnostic_command']=diagnostic_command
            (root/'native_command.json').write_text(json.dumps(commands,indent=2))
            with phase(timing_file,'native_subprocess_launch_through_exit'):
                subprocess.run(diagnostic_command,check=True)
        else:
            subprocess.run(command,check=True)
        commands[-1]['status']='returned'
        (root/'native_command.json').write_text(json.dumps(commands,indent=2))
    csvs=list(out.rglob('per_example_metrics*.csv'))
    if not csvs: raise RuntimeError('RFD4 did not write its per-design metrics table.')
    # Produce a portable result index using the annotations carried by native CIFs.
    paths=sorted(out.rglob('*_generated.cif')) if p['task']!='score_structures' else [Path(x['path']) for x in inventory]
    records=[]
    for path in paths:
        aa=load(path,True)
        chain_attr='chain_id'
        ids=aa.get_annotation(chain_attr)
        variable=~C_SEQ.mask(aa,default='generate')
        designed=[str(c) for c in np.unique(ids) if variable[ids==c].any()]
        if not designed:
            variable_coords=~C_CRD.mask(aa,default='generate')
            designed=[str(c) for c in np.unique(ids) if (variable_coords & aa.is_protein)[ids==c].any()]
        if not designed:
            raise ValueError(f'No designed protein/RNA chains could be identified in {path}.')
        sequences={}
        for c in np.unique(ids):
            chain=aa[ids==c]; starts=get_residue_starts(chain)
            codes=chain.res_name[starts]
            inverse={v:k for k,v in AA3.items()}
            sequences[str(c)]=''.join(inverse.get(str(x),str(x) if len(str(x))==1 else 'X') for x in codes)
        origin=next((x for x in inventory if x['name'] in path.parts or x['path']==str(path)),inventory[0] if len(inventory)==1 else {})
        _require_output_target_chains(p['task'],origin,sequences,path.name)
        records.append({'structure':str(path),'source_candidate':origin.get('source_candidate'),'designed_chains':designed,'target_chains':[str(c) for c in np.unique(ids) if c not in designed],
                        'sequences':sequences,'residue_count':len(get_residue_starts(aa))})
    (root/'result_index.json').write_text(json.dumps({'structures':records,'metrics_files':[str(x) for x in csvs]},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('request'); parser.add_argument('--prepare-only',action='store_true')
    args=parser.parse_args(); run(json.loads(Path(args.request).read_text()),args.prepare_only)
