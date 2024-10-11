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

def get_charge_free_aa(pocket_locations_dict, location_aa_dict, pH = 7):
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
            aa_3_letter = location_aa_dict[location]
            if aa_3_letter:
                aa_1_letter = IUPACData.protein_letters_3to1.get(aa_3_letter.capitalize(), "Unknown")
                charge = ProteinAnalysis(aa_1_letter).charge_at_pH(pH)
                #tprint(charge)
                pocket_charge_dict[pocket] += charge
            # TODO else 
            else:
                print("Unknown aa: ", aa_3_letter)
    
    return pocket_charge_dict

# tohle je tady spíš ze srandy
def get_charge_peptide(pocket_locations_dict, location_aa_dict, pH = 7):
    pocket_charge_dict = {}
    
    for pocket in pocket_locations_dict:
        pocket_charge_dict[pocket] = 0
        locations = pocket_locations_dict[pocket]

        peptide = ""
        for location in locations:
            aa_3_letter = location_aa_dict[location]
            if aa_3_letter:
                aa_1_letter = IUPACData.protein_letters_3to1.get(aa_3_letter.capitalize(), "Unknown")
                peptide += aa_1_letter
            # TODO else 
            else:
                print("Unknown aa: ", aa_3_letter)
            
            pocket_charge_dict[pocket] += ProteinAnalysis(peptide).charge_at_pH(pH)
    
    return pocket_charge_dict


if __name__ == "__main__":
    pocket_locations_dict = parse_predictions("electrostatics/structure.pdb_predictions.csv")
    locations_aa_dict = parse_residues("electrostatics/structure.pdb_residues.csv")
    print()
    pocket_charge_dict = get_charge_free_aa(pocket_locations_dict, locations_aa_dict)
    print(pocket_charge_dict)
    print()
    pocket_charge_dict2 = get_charge_peptide(pocket_locations_dict, locations_aa_dict)
    print(pocket_charge_dict2)
    print()