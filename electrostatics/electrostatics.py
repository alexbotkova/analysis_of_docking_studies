"""
Electrostatics
This script calculates the electrostatic charges of protein pockets based on data from PDB files 
obtained from https://prankweb.cz.

Author: Alexandra Botková
"""

import sys

sys.path.append("..")
from my_tools.my_parser import *

aa_charge_table = {
    'ALA': 0,
    'ARG': 1,
    'ASN': 0,
    'ASP': -1,
    'CYS': 0,
    'GLU': -1,
    'GLN': 0,
    'GLY': 0,
    'HIS': 1,
    'ILE': 0,
    'LEU': 0,
    'LYS': 1,
    'MET': 0,
    'PHE': 0,
    'PRO': 0,
    'SER': 0,
    'THR': 0,
    'TRP': 0,
    'TYR': 0,
    'VAL': 0
}

def get_pocket_charge(pocket_locations_dict, location_aa_dict, pH = 7):
    """
    Calculates the charge of each protein pocket based on the predicted amino acids and environmental pH.

    :param pocket_locations_dict: Dictionary where keys are pocket IDs and values are lists of chain and residue locations.
    :param location_aa_dict: Dictionary mapping chain-residue locations to their corresponding amino acid abbreviations.
    :return: Dictionary where keys are pocket IDs and values are the charge of the respective pockets.
    """
    pocket_charge_dict = {}
    
    for pocket in pocket_locations_dict:
        pocket_charge_dict[pocket] = 0
        locations = pocket_locations_dict[pocket]
        for location in locations:
            aa = location_aa_dict[location]
            if aa:
                pocket_charge_dict[pocket] += aa_charge_table[aa]
            # TODO else 
            else:
                print("Unknown aa: ", aa)
    
    return pocket_charge_dict

if __name__ == "__main__":
    pocket_locations_dict = parse_predictions("test_files/prankweb-2SRC/structure.pdb_predictions.csv")
    locations_aa_dict = parse_residues("test_files/prankweb-2SRC/structure.pdb_residues.csv")
    pocket_charge_dict = get_pocket_charge(pocket_locations_dict, locations_aa_dict)
    print(pocket_charge_dict)
    print()