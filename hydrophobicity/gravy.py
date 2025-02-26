import sys

sys.path.append("..")
from my_tools.my_parser import *

hydrophobicity_scale = {
    'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,
    'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
    'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6,
    'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
}

def get_pocket_gravy(pocket_locations_dict, location_aa_dict, pH = 7):
    """
    Calculates the GRAVY (grand average of hydropathy) of each protein pocket based on the predicted amino acids.

    :param pocket_locations_dict: Dictionary where keys are pocket IDs and values are lists of chain and residue locations.
    :param location_aa_dict: Dictionary mapping chain-residue locations to their corresponding amino acid abbreviations.
    :return: Dictionary where keys are pocket IDs and values are the averages of hydropathy of the respective pockets.
    """
    pocket_charge_dict = {}
    
    for pocket in pocket_locations_dict:
        pocket_charge_dict[pocket] = 0
        locations = pocket_locations_dict[pocket]
        for location in locations:
            aa = location_aa_dict[location]
            if aa:
                pocket_charge_dict[pocket] += hydrophobicity_scale[aa]
            # TODO else 
            else:
                print("Unknown aa: ", aa)
        pocket_charge_dict[pocket] /= len(locations) if locations else 0

    return pocket_charge_dict

if __name__=="__main__":
    pocket_locations_dict = parse_predictions("test_files/prankweb-2SRC/structure.pdb_predictions.csv")
    locations_aa_dict = parse_residues("test_files/prankweb-2SRC/structure.pdb_residues.csv")
    print(get_pocket_gravy(pocket_locations_dict, locations_aa_dict))