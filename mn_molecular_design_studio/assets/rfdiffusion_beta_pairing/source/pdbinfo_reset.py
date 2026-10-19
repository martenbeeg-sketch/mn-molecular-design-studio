#!/usr/bin/env python
# Util from bcov, reindex pdb to 1
# Usage: pdbinfo_reset.py <Your_pdb>
import os
import sys


from pyrosetta import *
from pyrosetta.rosetta import *

#flags = "-use_truncated_termini -extra_res_fa /Users/brian/Documents/baker/from/daniel/h3/H3i.fa.mm.params"
#for item in sys.argv[2:]:
#    flags += " " + item

init("-extra_res_fa /path/to/H3i.am1bcc.fa.mm.params")

pose = pose_from_pdb(sys.argv[1])


a = core.pose.PDBInfo(pose)

pose.pdb_info(a)



name = sys.argv[1]

had_gz = False
if name.endswith(".gz"):
    had_gz = True
    name = name[:-3]

if (name.endswith(".pdb")):
    name = name[:-4]

name = name + ".pdbinfo_reset.pdb"

if (had_gz):
    name = name + ".gz"

pose.dump_pdb(name)
