"""
Provides utility functions for parsing pocket prediction data from CSV files 
into PyMOL-friendly selection strings.

Functions:
    get_pocket_residues_dict: Converts residue IDs into selection strings by pocket.
    get_pocket_atomids_dict: Converts atom IDs into selection strings by pocket.
"""

import re
from my_parser import *

def get_pocket_residues_dict(pocket_data_df):
    """
    Generates a dictionary where keys are pocket names and values are residues joined in one string 
    formatted for the pymol selection command.

    :param pocket_data_df: A dataframe generated from a CSV file containing PDB predictions.
    :return: A dictionary where keys are pocket names and values are residue selection strings for PyMOL.
    """
    pocket_residues_dict = {}
    for _, row in pocket_data_df.iterrows():
            pocket_name = row['name'].strip()
            residue_ids = row['residue_ids']
            matches = re.findall(r'([A-Z])_(\d+)', residue_ids)
            if matches:
                chain_letter = matches[0][0] 
                residues = [residue for _, residue in matches]
                joined_residues = "+".join(residues)
            residues_selection = f"chain {chain_letter} and resi {joined_residues}"
            pocket_residues_dict[pocket_name] = residues_selection
    return pocket_residues_dict

def get_pocket_atomids_dict(pocket_data_df):
    """
    Generates a dictionary where keys are pocket names and values are atom ID selections
    formatted for the PyMOL selection command.

    :param pocket_data_df: A dataframe generated from a CSV file containing PDB predictions.
    :return: A dictionary where keys are pocket names and values are atom ID selection strings for PyMOL.
    """
    pocket_atomid_dict = {}
    for _, row in pocket_data_df.iterrows():
        pocket_name = row['name'].strip()
        atom_ids_str = str(row['surf_atom_ids'])
        atom_ids = re.findall(r'\d+', atom_ids_str)
        joined_ids = '+'.join(atom_ids)
        atom_selection = f"id {joined_ids}"
        pocket_atomid_dict[pocket_name] = atom_selection
    return pocket_atomid_dict

if __name__ == "__main__":
    predictions_filepath = "/Users/alexbotkova/analysis_of_docking_studies/test_files/urea/prankweb-2SRC/structure.cif_predictions.csv"
    df = get_df(predictions_filepath)
    res_dict = get_pocket_residues_dict(df)
    #print(res_dict)
    plus_counts = {}

    for key, value in res_dict.items():
        plus_counts[key] = value.count('+')

    #print(plus_counts)

    id_dict = get_pocket_atomids_dict(df)
    print(id_dict)