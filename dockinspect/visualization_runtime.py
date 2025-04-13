from pymol import cmd
import sys
import os
sys.path.insert(0, '/Users/alexbotkova/analysis_of_docking_studies/dockinspect')

from my_parser import *
from visualization import *

visualize(
    pdb_code='2src',
    vina_file='/Users/alexbotkova/Downloads/results/vina_dock.pdb',
    pocket_selection='',
    pose_num=1,
    mode=''
)
