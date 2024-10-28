"""
Electrostatics
This script calculates the electrostatic charges of protein pockets based on data from PDB files 
obtained from https://prankweb.cz.

Author: Alexandra Botková
"""

import pandas as pd
from Bio.Data import IUPACData
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def parse_predictions(predictions_filepath):
    """
    Parses a CSV file containing PDB predictions to generate a dictionary mapping pocket IDs 
    to lists of residue IDs.

    :param predictions_filepath: Path to the CSV file containing PDB predictions.
    :return: A dictionary where keys are pocket IDs and values are lists of residue IDs.
    """
    
    predictions_df = pd.read_csv(predictions_filepath)
    predictions_df.columns = predictions_df.columns.str.strip()

    pocket_locations_dict = {}

    for index, row in predictions_df.iterrows():
        pocket_name = row["name"].strip() 
        residue_ids = row["residue_ids"].split()  
        pocket_locations_dict[pocket_name] = residue_ids

    return pocket_locations_dict


def parse_residues(residues_filepath):
    """
    Parses a CSV file containing PDB residue information and generates a dictionary mapping  each chain and residue ID 
    to the corresponding amino acid three-letter abbreviation.

    :param residues_filepath: Path to the CSV file containing PDB residue data.
    :return: A dictionary where keys are residue IDs and values are the three-letter amino acid abbreviations at those locations.
    """

    residues_df = pd.read_csv(residues_filepath)
    residues_df.columns = residues_df.columns.str.strip()

    location_aa_dict = {}

    for _, row in residues_df.iterrows():
        chain = row["chain"].strip()  
        residue_id = row["residue_label"]  
        amino_acid = row["residue_name"].strip()  

        location_aa_dict[f"{chain}_{residue_id}"] = amino_acid

    return location_aa_dict

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
    pocket_locations_dict = parse_predictions("electrostatics/structure.pdb_predictions.csv")
    locations_aa_dict = parse_residues("electrostatics/structure.pdb_residues.csv")
    pocket_charge_dict = get_pocket_charge(pocket_locations_dict, locations_aa_dict)
    print(pocket_charge_dict)
