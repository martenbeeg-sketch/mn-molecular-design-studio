# Copyright (c) 2025 Isaac Sappington
# Licensed under the MIT License: https://opensource.org/license/MIT

import os,sys,glob,torch,random
import numpy as np
import argparse
import pyrosetta
import time
pyrosetta.init()

# example command: python ./make_interface_tensor.py --input_pdb /path/to.pdb --out_dir /out/dir/ --target_adj B100-110:B150-160:B190-200 --binderlen 65 --make_ss 9 --binder_ss E,E,E --binder_ss_len 10,10,10 --mask_to_ss H

def main():
    args=get_args()
    print(args)
    assert args.input_pdb is not None, 'Need to provide either an input pdb (--input_pdb)'

    os.makedirs(args.out_dir, exist_ok=True)
    pdb=args.input_pdb
    name=os.path.split(pdb)[1][:-4]
    binderlen=int(args.binderlen)
    target_adj=args.target_adj
    binder_ss_len=args.binder_ss_len
    binder_ss=args.binder_ss
    mask_to_ss=args.mask_to_ss
    make_ss=args.make_ss
    target_adjs = target_adj.split(':')
    adjlens = binder_ss_len.split(',')
    not_adj=args.not_adj

    assert len(adjlens) == len(target_adjs), f'number of provided adjlens, target_adjs must be equal. Provided {len(adjlens)} adjlens in adjlen_list and {len(target_adjs)} target_adjs in target_adj_list'

    # iterate for binding ss at each possible index
    # Include the final valid placement (the upstream script omitted it).
    for index in range(binderlen - int(adjlens[0]) + 1):

        # create ss matrix, set a
        zero = time.time()
        secstruc_dict=extract_secstruc(pdb)
        xyz,_,_ = parse_pdb_torch(pdb)
        complex_xyz, complex_secstruct_dict = add_binder(binderlen, xyz, secstruc_dict)
    
        # if multiple ss are desired, set secondary structure identities for those multiple ss at each possible, unoccupied index
        if len(adjlens)>0:
            counter = 0
            ss_list = [[index, index+int(adjlens[0])]]
            binder_adj_list = []
            list_binder_ss(binderlen, int(adjlens[0]), ss_list, len(adjlens), binder_adj_list)

            for adj_set in binder_adj_list:
                counter += 1
                new_secstruct_dict = complex_secstruct_dict.copy()
                    
                # insert ss content in long masked regions (defined by make_ss+4)
                if make_ss:
                    masks = long_masks(adj_set, binderlen)
                    unmask_ss(masks, binderlen, make_ss, new_secstruct_dict, mask_to_ss)

                for adj in range(len(adjlens)):
                    set_ss(adj_set[adj][0], int(adjlens[adj]), binderlen, new_secstruct_dict, binder_ss.split(',')[adj])

                ss, idx = ss_to_tensor(new_secstruct_dict)
                block_adj = construct_block_adj_matrix(torch.FloatTensor(ss), torch.tensor(complex_xyz)).float()
                mask_adjacency(block_adj, binderlen, new_secstruct_dict)
                    
                # set your ss adjacencies for each target contig
                for adj in range(len(adjlens)):
                    target_contig_str = target_adjs[adj].split(',')
                    for contig in target_contig_str:
                        target = parse_contig(contig)[0]                            
                        pair_adjacency(block_adj, adj_set[adj][0], int(adjlens[adj]), binderlen, target, new_secstruct_dict)
                    
                # set every binder residue as being non-adjacent to non_adj residues
                if not_adj:
                    not_adj_contigs = not_adj.split(',')
                    for contig in not_adj_contigs:
                        not_hot = parse_contig(contig)[0]
                        negative_adjacency(block_adj, binderlen, not_hot, new_secstruct_dict)
                else:
                    not_adj=''

                ss_tens, mask = mask_ss(ss, idx, max_mask=0)
                ss_argmax = torch.argmax(ss_tens[:,:4], dim=1).float()
                torch.save(ss_argmax, os.path.join(args.out_dir, f'{name}_adj{"-".join([str(x[0]) for x in adj_set])}_{counter}_ss.pt'))
                torch.save(block_adj, os.path.join(args.out_dir, f'{name}_adj{"-".join([str(x[0]) for x in adj_set])}_{counter}_adj.pt'))
                parent_dir = os.path.dirname(os.path.normpath(args.out_dir))
                out_name = os.path.basename(os.path.normpath(args.out_dir)) + ".txt"
                outfile = os.path.join(parent_dir, out_name)
                os.system(f'echo {name}_adj{"-".join([str(x[0]) for x in adj_set])}_{counter} >> {outfile}')
            print(f'wrote adjacencies for {name} with adj {"-".join([str(x[0]) for x in adj_set])} {time.time()-zero}')
            
def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input_pdb", required=True, help="path to input pdb. Either provide this of path to directory of pdbs (--pdb_dir)", default=None)
    parser.add_argument("--out_dir",dest="out_dir", required=True, help='need to specify an output path')
    parser.add_argument("--binderlen", required=True, help='provide the length of the binder you want to design')
    parser.add_argument("--target_adj", required=False, help='specify list of target contigs strings, with each list separated by a semicolon (;), for your binder to be adjacent (e.g.. B100-110,B115-115;B200-205,B210-210)')
    parser.add_argument("--binder_ss_len", required=False, default=False, help='provide the list of lengths that you want binder_ss to be, respective to that list')
    parser.add_argument("--binder_ss", required=False, default='M', help='provide a list (of length len(adjlens)), to set the secondary structure of your binder that are adjacent to multiadj_list, respectively') 
    parser.add_argument("--mask_to_ss", required=False, default='H', help='provide the ss type (H, E , or L) that you would like to fill long mask regions')
    parser.add_argument("--make_ss",required=False, default=False, help='set the helix length if you want helices to fill long gaps between ss blocks (indicated by min_mask_to_fill')
    parser.add_argument("--not_adj", required=False, default=False, help='specify contigs on target for your binder to avoid')
    args = parser.parse_args()
    return args

def parse_contig(contig_crop):
        """
        Takes contig input and parses
        """
        contig_list = []
        for contig in contig_crop.split(" "):
            subcon=[]
            for crop in contig.split(","):
                if crop[0].isalpha():
                    subcon.extend([(crop[0], p) for p in np.arange(int(crop.split("-")[0][1:]), int(crop.split("-")[1])+1)])
            contig_list.append(subcon)

        return contig_list

def add_binder(binderlen, xyz, secstruct_dict):
    """
    Adds masked secondary structure glycines to the target pdb with [0, 0, 0] coordinates for all atoms. New target residue indices will be idx+binderlen
    """
    xyz_binder = np.zeros((binderlen, 27, 3), dtype=np.float32)
    complex_xyz = np.concatenate((xyz_binder, xyz))
    secstruct_dict_clone = dict(secstruct_dict)
    for g in range(binderlen):
        secstruct_dict_clone['sequence'].insert(0, 'G')
        secstruct_dict_clone['idx'].append(secstruct_dict['idx'][-1]+1)
        secstruct_dict_clone['ss'].insert(0, 'M')
    
    return complex_xyz, secstruct_dict_clone

def set_ss(index, ss_len, binderlen, secstruct_dict, ss):
    """
    Sets consecutive residues in a binder to make a ss block. Used in a loop to generate all ss locations.
    """
    if index < 0 or index + ss_len > binderlen:
        raise ValueError("The requested binder secondary-structure segment exceeds binderlen.")
    for position in range(index, index + ss_len):
        secstruct_dict['ss'][position] = ss
            
def list_binder_ss(binderlen, ss_len, ss_indices, len_adjlens, ss_list):
    """
    Lists all the possible positions (separated by at least 2 residues from your binding ss)
    """
    len_adjlens -= 1
    
    if len_adjlens>0:
        for index in range(binderlen-ss_len+1):
            index_list = ss_indices[:]
            available_for_ss = True
            for ss in index_list:
                if ss[0]-2 <= index <= ss[1]+2 or ss[0]-2 <= index+ss_len <= ss[1]+2:
                    available_for_ss = False
            if available_for_ss == True:
                index_list.append([index, index+ss_len])
                list_binder_ss(binderlen, ss_len, index_list, len_adjlens, ss_list)
    else:
        #ss_list.sort()
        if not ss_indices in ss_list:
            ss_list.append(ss_indices)

def long_masks(ss_list, binderlen):
    """
    List all the masked regions in your secondary structure conditioning (i.e., where there are no secondary structures already set)
    """
    long_masks = []
    limits = ss_list[:]
    limits.insert(0, [0,0])
    limits.append([binderlen, binderlen])
    for limiter in range(len(limits[:-1])):
        long_masks.append([limits[limiter][1], limits[limiter+1][0]])
    return long_masks

def unmask_ss(masks, binderlen, ss_length, secstruct_dict, mask_to_ss):
    """
    Set helical identity for 'ss_length' residues located in the center of long masked regions (ss_length)
    """
    for mask in masks:
        if mask[1]-mask[0] >= ss_length+4:
            if mask[0] == 0:
                for resi in range(1, int(ss_length)+1):
                    secstruct_dict['ss'][resi] = mask_to_ss
            elif mask[1] == binderlen:
                for resi in range(binderlen-(int(ss_length)+1), binderlen-1):
                    secstruct_dict['ss'][resi] = mask_to_ss
            else:
                mask_midpoint = int((mask[1]-mask[0])/2)
                for resi in range(mask_midpoint-int((int(helix_length)/2)), mask_midpoint+int((int(helix_length)/2))):
                    secstruct_dict['ss'][resi] = mask_to_ss

def mask_adjacency(block_adj, binderlen, secstruct_dict):
    """
    Masks all interactions between binder and target
    """
    for i in range(len(secstruct_dict['ss'])):
        for j in range(binderlen):
            block_adj[i][j] = 2
            block_adj[j][i] = 2

def pair_adjacency(block_adj, index, ss_len, binderlen, contig, secstruct_dict):
    """
    Sets adjacency of binder ss with the target residues
    """
    for ss_b in range(index, index + ss_len):
        for ss_t in range(len(contig)):
            # PDB residue numbers are one-based. Subtract the first target
            # index once so residue 1 maps to binderlen, not binderlen - 1.
            target_index = contig[ss_t][1] + binderlen - secstruct_dict['idx'][0]
            block_adj[ss_b][target_index] = 1
            block_adj[target_index][ss_b] = 1

def negative_adjacency(block_adj, binderlen, contig, secstruct_dict):
    """
    Sets all binder residues as not adjacent to not_adj residues
    """
    for binder_res in range(binderlen):
        for not_adj_res in range(len(contig)):
            target_index = contig[not_adj_res][1] + binderlen - secstruct_dict['idx'][0]
            block_adj[binder_res][target_index] = 0
            block_adj[target_index][binder_res] = 0

def extract_secstruc(fn):
    pdb=parse_pdb(fn)
    idx = pdb['idx']
    dssp = pyrosetta.rosetta.core.scoring.dssp
    pose = pyrosetta.io.pose_from_pdb(fn)
    dssp.Dssp(pose).insert_ss_into_pose(pose, True)
    aa_sequence = pose.sequence()
    secstruct = pose.secstruct()
    secstruc_dict = {'sequence':[i for i in aa_sequence],
                 'idx':[int(i) for i in idx],
                 'ss':[i for i in secstruct]}
    return secstruc_dict

def ss_to_tensor(ss):
    """
    Function to convert ss files to indexed tensors
    0 = Helix
    1 = Strand
    2 = Loop
    3 = Mask/unknown
    4 = idx for pdb
    """
    ss_conv = {'H':0,'E':1,'L':2,'M':3}
    idx = np.array(ss['idx'])
    ss_int = np.array([int(ss_conv[i]) for i in ss['ss']])
    return ss_int, idx

def mask_ss(ss, idx, min_mask = 0, max_mask = 1.0):
    mask_prop = random.uniform(min_mask, max_mask)
    transitions = np.where(ss[:-1] - ss[1:] != 0)[0] #gets last index of each block of ss
    stuck_counter = 0
    while len(ss[ss == 3])/len(ss) < mask_prop or stuck_counter > 100:
        width = random.randint(1,9)
        start = random.choice(transitions)
        offset = random.randint(-8,1)
        try:

            ss[start+offset:start+offset+width] = 3
        except:
            stuck_counter += 1
            pass
    ss = torch.tensor(ss)
    ss = torch.nn.functional.one_hot(ss, num_classes=4)
    ss = torch.cat((ss, torch.tensor(idx)[...,None]), dim=-1)
#     mask = torch.where(torch.argmax(ss[:,:-1], dim=-1) == 3, False, True)
    mask=torch.tensor(np.where(np.argmax(ss[:,:-1].numpy(), axis=-1) == 3))
    return ss, mask

def generate_Cbeta(N,Ca,C):
    # recreate Cb given N,Ca,C
    b = Ca - N 
    c = C - Ca
    a = torch.cross(b, c, dim=-1)
    #Cb = -0.58273431*a + 0.56802827*b - 0.54067466*c + Ca
    # fd: below matches sidechain generator (=Rosetta params)
    Cb = -0.57910144*a + 0.5689693*b - 0.5441217*c + Ca

    return Cb

def get_pair_dist(a, b): 
    """calculate pair distances between two sets of points
    
    Parameters
    ----------
    a,b : pytorch tensors of shape [batch,nres,3]
          store Cartesian coordinates of two sets of atoms
    Returns
    -------
    dist : pytorch tensor of shape [batch,nres,nres]
           stores paitwise distances between atoms in a and b
    """

    dist = torch.cdist(a, b, p=2)
    return dist


def construct_block_adj_matrix( sstruct, xyz, cutoff=6, include_loops=False ):
    ''' 
    Given a sstruct specification and backbone coordinates, build a block adjacency matrix.

    Input:
    
        sstruct (torch.FloatTensor): (L) length tensor with numeric encoding of sstruct at each position

        xyz (torch.FloatTensor): (L,3,3) tensor of Cartesian coordinates of backbone N,Ca,C atoms

        cutoff (float): The Cb distance cutoff under which residue pairs are considered adjacent
                        By eye, Nate thinks 6A is a good Cb distance cutoff

    Output:

        block_adj (torch.FloatTensor): (L,L) boolean matrix where adjacent secondary structure contacts are 1
    '''

    L = xyz.shape[0]
    
    # three anchor atoms
    N  = xyz[:,0]
    Ca = xyz[:,1]
    C  = xyz[:,2]
    
    # recreate Cb given N,Ca,C
    Cb = generate_Cbeta(N,Ca,C)
    
    # May need a batch dimension - NRB
    dist = get_pair_dist(Cb,Cb) # [L,L]
    dist[torch.isnan(dist)] = 999.9

    dist += 999.9*torch.eye(L,device=xyz.device)
    # Now we have dist matrix and sstruct specification, turn this into a block adjacency matrix
    # There is probably a way to do this in closed-form with a beautiful einsum but I am going to do the loop approach
    
    # First: Construct a list of segments and the index at which they begin and end
    in_segment = True
    segments = []

    begin = -1
    end = -1

    for i in range(sstruct.shape[0]):
        # Starting edge case
        if i == 0:
            begin = 0 
            continue

        if not sstruct[i] == sstruct[i-1]:
            end = i 
            segments.append( (sstruct[i-1], begin, end) )

            begin = i

    # Ending edge case: last segment is length one
    if not end == sstruct.shape[0]:
        segments.append( (sstruct[-1], begin, sstruct.shape[0]) )


    block_adj = torch.zeros_like(dist)
    for i in range(len(segments)):
        curr_segment = segments[i]

        if curr_segment[0] == 2 and not include_loops: continue

        begin_i = curr_segment[1]
        end_i = curr_segment[2]
        for j in range(i+1, len(segments)):
            j_segment = segments[j]

            if j_segment[0] == 2 and not include_loops: continue

            begin_j = j_segment[1]
            end_j = j_segment[2]

            if torch.any( dist[begin_i:end_i, begin_j:end_j] < cutoff ):
                # Matrix is symmetic
                block_adj[begin_i:end_i, begin_j:end_j] = torch.ones(end_i - begin_i, end_j - begin_j)
                block_adj[begin_j:end_j, begin_i:end_i] = torch.ones(end_j - begin_j, end_i - begin_i)
    return block_adj

def parse_pdb_torch(filename):
    lines = open(filename,'r').readlines()
    return parse_pdb_lines_torch(lines)

#'''
def parse_pdb_lines_torch(lines):

    # indices of residues observed in the structure
    pdb_idx = [( l[21:22].strip(), int(l[22:26].strip()) ) for l in lines if l[:4]=="ATOM" and l[12:16].strip()=="CA"]  # chain letter, res num
 
    # 4 BB + up to 10 SC atoms
    xyz = np.full((len(pdb_idx), 27, 3), np.nan, dtype=np.float32)
    for l in lines:
        if l[:4] != "ATOM":
            continue
        chain, resNo, atom, aa = l[21:22], int(l[22:26]), ' '+l[12:16].strip().ljust(3), l[17:20]
        idx = pdb_idx.index((chain,resNo))
        for i_atm, tgtatm in enumerate(aa2long[aa2num[aa]]):
            if tgtatm == atom:
                xyz[idx,i_atm,:] = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
                break
    # save atom mask
    mask = np.logical_not(np.isnan(xyz[...,0]))
    xyz[np.isnan(xyz[...,0])] = 0.0

    return xyz,mask,np.array(pdb_idx)

def parse_pdb(filename, **kwargs):
    '''extract xyz coords for all heavy atoms'''
    lines = open(filename,'r').readlines()
    return parse_pdb_lines(lines, **kwargs)

def parse_pdb_lines(lines, parse_hetatom=False, ignore_het_h=True):
    # indices of residues observed in the structure
    res = [(l[22:26],l[17:20]) for l in lines if l[:4]=="ATOM" and l[12:16].strip()=="CA"]
    seq = [aa2num[r[1]] if r[1] in aa2num.keys() else 20 for r in res]
    pdb_idx = [( l[21:22].strip(), int(l[22:26].strip()) ) for l in lines if l[:4]=="ATOM" and l[12:16].strip()=="CA"]  # chain letter, res num
    
    # 4 BB + up to 10 SC atoms
    xyz = np.full((len(res), 27, 3), np.nan, dtype=np.float32)
    for l in lines:
        if l[:4] != "ATOM":
            continue
        chain, resNo, atom, aa = l[21:22], int(l[22:26]), ' '+l[12:16].strip().ljust(3), l[17:20]
        idx = pdb_idx.index((chain,resNo))
        for i_atm, tgtatm in enumerate(aa2long[aa2num[aa]]):
            if tgtatm is not None and tgtatm.strip() == atom.strip(): # ignore whitespace
                xyz[idx,i_atm,:] = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
                break
        
    # save atom mask
    mask = np.logical_not(np.isnan(xyz[...,0]))
    xyz[np.isnan(xyz[...,0])] = 0.0 
    # remove duplicated (chain, resi)
    new_idx = []
    i_unique = []
    for i,idx in enumerate(pdb_idx):
        if idx not in new_idx:
            new_idx.append(idx)
            i_unique.append(i)
    
    pdb_idx = new_idx
    xyz = xyz[i_unique]
    mask = mask[i_unique]
    seq = np.array(seq)[i_unique]

    out = {'xyz':xyz, # cartesian coordinates, [Lx14]
            'mask':mask, # mask showing which atoms are present in the PDB file, [Lx14]
            'idx':np.array([i[1] for i in pdb_idx]), # residue numbers in the PDB file, [L]
            'seq':np.array(seq), # amino acid sequence, [L]
            'pdb_idx': pdb_idx,  # list of (chain letter, residue number) in the pdb file, [L]
           }
    # heteroatoms (ligands, etc)
    if parse_hetatom:
        xyz_het, info_het = [], []
        for l in lines:
            if l[:6]=='HETATM' and not (ignore_het_h and l[77]=='H'):
                info_het.append(dict(
                    idx=int(l[7:11]),
                    atom_id=l[12:16],
                    atom_type=l[77],
                    name=l[16:20]
                ))
                xyz_het.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])

        out['xyz_het'] = np.array(xyz_het)
        out['info_het'] = info_het

    return out

num2aa=[
    'ALA','ARG','ASN','ASP','CYS',
    'GLN','GLU','GLY','HIS','ILE',
    'LEU','LYS','MET','PHE','PRO',
    'SER','THR','TRP','TYR','VAL',
    'UNK','MAS',
    ]   

aa2num= {x:i for i,x in enumerate(num2aa)}
# full sc atom representation (Nx14)
aa2long=[
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","3HB ",  None,  None,  None,  None,  None,  None,  None,  None), # ala
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," NE "," CZ "," NH1"," NH2",  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HD ","2HD "," HE ","1HH1","2HH1","1HH2","2HH2"), # arg
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," ND2",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HD2","2HD2",  None,  None,  None,  None,  None,  None,  None), # asn
    (" N  "," CA "," C  "," O  "," CB "," CG "," OD1"," OD2",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ",  None,  None,  None,  None,  None,  None,  None,  None,  None), # asp
    (" N  "," CA "," C  "," O  "," CB "," SG ",  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB "," HG ",  None,  None,  None,  None,  None,  None,  None,  None), # cys
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," NE2",  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HE2","2HE2",  None,  None,  None,  None,  None), # gln
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," OE1"," OE2",  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ",  None,  None,  None,  None,  None,  None,  None), # glu
    (" N  "," CA "," C  "," O  ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None," H  ","1HA ","2HA ",  None,  None,  None,  None,  None,  None,  None,  None,  None,  None), # gly
    (" N  "," CA "," C  "," O  "," CB "," CG "," ND1"," CD2"," CE1"," NE2",  None,  None,  None,  None," H  "," HA ","1HB ","2HB "," HD2"," HE1"," HE2",  None,  None,  None,  None,  None,  None), # his
    (" N  "," CA "," C  "," O  "," CB "," CG1"," CG2"," CD1",  None,  None,  None,  None,  None,  None," H  "," HA "," HB ","1HG2","2HG2","3HG2","1HG1","2HG1","1HD1","2HD1","3HD1",  None,  None), # ile
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB "," HG ","1HD1","2HD1","3HD1","1HD2","2HD2","3HD2",  None,  None), # leu
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD "," CE "," NZ ",  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HD ","2HD ","1HE ","2HE ","1HZ ","2HZ ","3HZ "), # lys
    (" N  "," CA "," C  "," O  "," CB "," CG "," SD "," CE ",  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","1HG ","2HG ","1HE ","2HE ","3HE ",  None,  None,  None,  None), # met
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ ",  None,  None,  None," H  "," HA ","1HB ","2HB "," HD1"," HD2"," HE1"," HE2"," HZ ",  None,  None,  None,  None), # phe
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD ",  None,  None,  None,  None,  None,  None,  None," HA ","1HB ","2HB ","1HG ","2HG ","1HD ","2HD ",  None,  None,  None,  None,  None,  None), # pro
    (" N  "," CA "," C  "," O  "," CB "," OG ",  None,  None,  None,  None,  None,  None,  None,  None," H  "," HG "," HA ","1HB ","2HB ",  None,  None,  None,  None,  None,  None,  None,  None), # ser
    (" N  "," CA "," C  "," O  "," CB "," OG1"," CG2",  None,  None,  None,  None,  None,  None,  None," H  "," HG1"," HA "," HB ","1HG2","2HG2","3HG2",  None,  None,  None,  None,  None,  None), # thr
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," NE1"," CE2"," CE3"," CZ2"," CZ3"," CH2"," H  "," HA ","1HB ","2HB "," HD1"," HE1"," HZ2"," HH2"," HZ3"," HE3",  None,  None,  None), # trp
    (" N  "," CA "," C  "," O  "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ "," OH ",  None,  None," H  "," HA ","1HB ","2HB "," HD1"," HE1"," HE2"," HD2"," HH ",  None,  None,  None,  None), # tyr
    (" N  "," CA "," C  "," O  "," CB "," CG1"," CG2",  None,  None,  None,  None,  None,  None,  None," H  "," HA "," HB ","1HG1","2HG1","3HG1","1HG2","2HG2","3HG2",  None,  None,  None,  None), # val
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","3HB ",  None,  None,  None,  None,  None,  None,  None,  None), # unk
    (" N  "," CA "," C  "," O  "," CB ",  None,  None,  None,  None,  None,  None,  None,  None,  None," H  "," HA ","1HB ","2HB ","3HB ",  None,  None,  None,  None,  None,  None,  None,  None), # mask
]

print('functions loaded')

if __name__ == "__main__":
    main()
