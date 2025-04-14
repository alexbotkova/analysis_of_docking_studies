from pymol import cmd
import sys
import os
sys.path.insert(0, '/Users/alexbotkova/analysis_of_docking_studies/dockinspect')

from my_parser import *
from visualization import *

visualize(
    pdb_code='2src',
    vina_file='/Users/alexbotkova/analysis_of_docking_studies/test_files/urea/result-2025-03-30T21_20_58.783Z/out_vina.pdbqt',
    pocket_selection='id 2407+2408+2409+2418+2423+2713+2719+2741+2742+2757+2832+2834+2863+2864+2900',
    pose_num=1,
    mode=''
)
