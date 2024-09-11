"""
Electrostatics
This script calculates the electrostatic charges of protein pockets based on data from PDB files 
obtained from https://prankweb.cz.

Author: Alexandra Botková
"""

import pandas as pd
from pI_aa_table import pI_table

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


def get_aa_charge(aa, pH):
    """
    Calculate the charge of an amino acid based on its isoelectric point (pI) and the given pH.

    :param aa: Three-letter amino acid abbreviation.
    :param pH: Environmental pH value.
    :return: The charge of the amino acid (1 for positive, 0 for neutral, -1 for negative).
    """

    pI = pI_table.get(aa)

    if pH < pI:
        return 1  # Positive charge below pI
    elif pH > pI:
        return -1  # Negative charge above pI
    else:
        return 0  # Neutral charge at pI

def get_pocket_charge(pocket_locations_dict, location_aa_dict, pH = 7):
    """
    Calculates the charge of each protein pocket based on the predicted amino acids and environmental pH.

    :param pocket_locations_dict: Dictionary where keys are pocket IDs and values are lists of chain and residue locations.
    :param location_aa_dict: Dictionary mapping chain-residue locations to their corresponding amino acid abbreviations.
    :param pH: Environmental pH value (default is 7).
    :return: Dictionary where keys are pocket IDs and values are the charge of the respective pockets.
    """

    pocket_charge_dict = {}

    for pocket in pocket_locations_dict:
        pocket_charge_dict[pocket] = 0

        locations = pocket_locations_dict[pocket]
        for location in locations:
            aa = location_aa_dict[location] 
            if aa:
                pocket_charge_dict[pocket] += get_aa_charge(aa, pH)
            # TODO else 
            else:
                print("Unknown aa: ", aa)

    return pocket_charge_dict

if __name__ == "__main__":
    pocket_locations_dict = parse_predictions("electrostatics/structure.pdb_predictions.csv")
    locations_aa_dict = parse_residues("electrostatics/structure.pdb_residues.csv")
    pocket_charge_dict = get_pocket_charge(pocket_locations_dict, locations_aa_dict)
    pocket_charge_dict = get_pocket_charge(pocket_locations_dict, locations_aa_dict, 0)
    pocket_charge_dict = get_pocket_charge(pocket_locations_dict, locations_aa_dict, 15)
    print(pocket_charge_dict)



