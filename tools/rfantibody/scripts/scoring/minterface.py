#!/usr/bin/env python3

'''
    This is a script which takes a quiver file of antibodies and scores them
    based on Rosetta metrics
'''
import sys
import os
import argparse
import uuid

from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.pack.guidance_scoreterms.sap import calculate_sap
init('-beta_nov16 -mute all')

from quiver import Quiver

#### Parse arguments
##################################

parser = argparse.ArgumentParser(description='Score antibodies')
parser.add_argument('-q', '--quiver', help='Quiver file', required=True)

args = parser.parse_args( sys.argv[1:] )

# Load in Rosetta function from the util xml
xml = protocols.rosetta_scripts.XmlObjects.create_from_file('/net/databases/antibody/software/rosetta/util_minterface.xml')

minterface = xml.get_mover('minimize_interface')

ab_ddg = xml.get_filter('ab_ddg')
nb_ddg = xml.get_filter('nb_ddg')

fv_net_charge = xml.get_filter('fv_net_charge')

if ( isinstance(ab_ddg, pyrosetta.rosetta.protocols.filters.StochasticFilter) ):
    ab_ddg = ab_ddg.subfilter()
if ( isinstance(nb_ddg, pyrosetta.rosetta.protocols.filters.StochasticFilter) ):
    nb_ddg = nb_ddg.subfilter()
if ( isinstance(fv_net_charge, pyrosetta.rosetta.protocols.filters.StochasticFilter) ):
    fv_net_charge = fv_net_charge.subfilter()

sap_score = xml.get_simple_metric('sap_score')

# Check if the quiver file exists
if not os.path.isfile(args.quiver):
    print('Quiver file does not exist')
    sys.exit()

def calulate_sap_on_pose(pose):
    
    # Extract the chains which are not T
    # Make a blank pose
    # Add the chains to the blank pose
    tmppose = None
    for chain in pose.split_by_chain():
        if chain.pdb_info().chain(1) != 'T':
            if tmppose is None:
                tmppose = chain
            else:
                tmppose.append_pose_by_jump(chain, tmppose.total_residue())

    score_selector = rosetta.core.select.residue_selector.TrueResidueSelector()
    sap_calculate_selector = rosetta.core.select.residue_selector.TrueResidueSelector()
    sasa_selector = rosetta.core.select.residue_selector.TrueResidueSelector()

    # Calculate the per-residue SAP scores
    sap_score = calculate_sap(tmppose, score_selector, sap_calculate_selector, sasa_selector)

    return sap_score

def rechain_pose(pose, chainorder):
    '''
        This function takes a pose and rechains it the requested order.
        This function retains the reslabels of the original pose.

    '''

    # First we will check if the pose is already in the requested order
    # If it is, we will return the pose

    splits = pose.split_by_chain()
    correctorder = True
    for i in range(1, len(splits) + 1):
        if splits[i].pdb_info().chain(1) != chainorder[i-1]:
            correctorder = False
            break
    
    if correctorder:
        return pose

    # If we get here, the pose is not in HLT order
    # We will rechain it

    # Check that the pose has the requested number of chains
    if len(splits) != len(chainorder):
        raise Exception(f'Pose has {len(splits)} chains, but {len(chainorder)} chains were requested')

    retpose = None
    for chain in chainorder:
        for subchain in splits:
            if subchain.pdb_info().chain(1) == chain:
                if retpose is None:
                    retpose = subchain
                else:
                    retpose.append_pose_by_jump(subchain, retpose.total_residue())

    # Now we will grab the reslabels from the original pose, and assign them
    # to a new PDBInfo object

    # We will also assign the original chains to the new PDBInfo object
    pdbinfo = core.pose.PDBInfo(retpose)

    for resi in range(1, retpose.total_residue() + 1):
        # Set the chain
        pdbinfo.set_resinfo(resi, retpose.pdb_info().chain(resi), resi)

        # Set the reslabel
        l = retpose.pdb_info().get_reslabels(resi)

        if len(l) > 0:
            pdbinfo.add_reslabel(resi, l[1])

    retpose.pdb_info(pdbinfo)

    return retpose

def pose_to_pdblines(pose):
    '''
    Convert a pose to a list of pdb lines
    '''
    # Dump a unique random pdb file
    pdbfile = f'temp_{uuid.uuid4()}.pdb'

    pose.dump_pdb(pdbfile)

    with open(pdbfile, 'r') as f:
        pdblines = f.readlines()

    os.remove(pdbfile)

    return pdblines

#### Main Loop
##################################

inqv  = Quiver(args.quiver, 'r')
outqv = Quiver('out.qv', 'w')

donetags = outqv.get_tags()

for tag in inqv.get_tags():

    # Check for finished tags
    if tag in donetags:
        print(f'{tag} already done, skipping')
        continue

    print(f"Processing design with tag: {tag}")

    pose = Pose()
    rosetta.core.import_pose.pose_from_pdbstring(pose, ''.join(inqv.get_pdblines(tag)))

    # Determine whether this is a nanobody or an antibody
    # We will do this by looking at the ids of each chain
    # If there are 2 chains, it is a nanobody
    # If there are 3 chains, it is an antibody

    # Get the chain ids
    # Get the pdb info
    pdbinfo = pose.pdb_info()
    unique_chain_ids = set()
    for i in range(1, pose.total_residue() + 1):
        unique_chain_ids.add(pdbinfo.chain(i))

    if 'T' not in unique_chain_ids:
        raise Exception('No chain T found in pose. Cannot calculate interface metrics '\
                        'without a target')

    igtype = None
    if 'H' in unique_chain_ids and 'L' in unique_chain_ids:

        igtype = 'ab'
        chainorder = ['H', 'L', 'T']

    elif 'H' in unique_chain_ids and 'L' not in unique_chain_ids:
        igtype = 'nb'
        chainorder = ['H', 'T']
    
    elif 'H' not in unique_chain_ids and 'L' in unique_chain_ids:
        igtype = 'nb'
        chainorder = ['L', 'T']

    else:
        raise Exception('No heavy or light chain found in pose. Cannot calculate interface metrics '\
                        'without a heavy or light chain')
    
    ### Ensure the pose is in HLT order since this is important for the ddg calculation
    pose = rechain_pose(pose, chainorder)

    #### Calculate metrics
    ##################################

    # minimize the pose
    minterface.apply(pose)

    if igtype == 'ab':
        ddg = ab_ddg.compute(pose)

    elif igtype == 'nb':
        ddg = nb_ddg.compute(pose)

    fv_charge = fv_net_charge.compute(pose)
    sap_score = calulate_sap_on_pose(pose)

    scorestr = f"ddg={ddg}|fv_charge={fv_charge}|sap_score={sap_score}"

    pdblines = pose_to_pdblines(pose)

    # Add to the quiver file
    outqv.add_pdb(pdblines, tag, scorestr)
