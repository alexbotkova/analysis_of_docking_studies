"""
Provides utility functions for parsing pocket prediction data from CSV files 
into PyMOL-friendly selection strings.

Functions:
    get_pocket_residues_dict: Converts residue IDs into selection strings by pocket.
"""

import re
from pandas import DataFrame

def get_pocket_residues_dict(pocket_data_df: DataFrame) -> dict:
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