import re
import freesasa

import sys
sys.path.append("..")
from my_tools.my_parser import *
from my_tools.pymol_tools import *

if __name__=="__main__":
    predictions_filepath = "test_files/prankweb-2SRC/structure.pdb_predictions.csv"
    pocket_data_df = get_df(predictions_filepath)
    pocket_residues_dict = get_pocket_residues_dict(pocket_data_df)
    commands =  [str(key) + ", " + str(value) for key, value in pocket_residues_dict.items()]

    structure = freesasa.Structure("test_files/result-2024-12-10T22_19_39.203Z/structure.pdbqt")
    result = freesasa.calc(structure)
    print(freesasa.selectArea(commands, structure, result))