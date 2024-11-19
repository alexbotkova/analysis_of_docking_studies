"""
Hydrogen bonds

This script calculates the number of hydrogen bonds between protein pockets and ligands based on data from PDB files 
obtained from https://prankweb.cz.

Author: Alexandra Botková
"""

import re
import numpy as np
import pandas as pd
import pymol
from pymol import cmd

def get_pocket_data_df(predictions_filepath):
    """
    Parses a CSV file containing PDB predictions to generate an according dataframe.

    :param predictions_filepath: Path to the CSV file containing PDB predictions.
    :return: Dataframe from the given CSV file.
    """
    pocket_data_df = pd.read_csv(predictions_filepath)
    pocket_data_df.columns = pocket_data_df.columns.str.strip()
    return pocket_data_df


def get_pocket_residues_dict(pocket_data_df):
    """
    Generates a dictionary where keys are pocket names and values are residues joined in one string 
    formatted for the pymol selection command.

    :param pocket_data_df: A dataframe generated from a CSV file containing PDB predictions.
    :return: A dictionary where keys are pocket names and values are residues joined in one string 
    formatted for the pymol selection command.
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
            residues_selection = f"resi {joined_residues} and chain {chain_letter}"
            pocket_residues_dict[pocket_name] = residues_selection
    return pocket_residues_dict

def get_models_avg_coordinates(out_vina_filepath):
    """
    Averages the coordinates of the atoms of the ligand.

    :param out_vina_filepath: Path to a PDBQT file containing info obtained from the docking.
    :return: A list where on the i-th index are the average coordinates of the ligand from the (i+1)-th model.
    """
    models_avg_coordinates = []
    current_coordinates = []
    with open(out_vina_filepath, 'r') as file:
        for line in file:
            if line.startswith("MODEL"):
                if current_coordinates:
                    models_avg_coordinates.append(np.mean(current_coordinates, axis=0))
                    current_coordinates = []
            elif line.startswith("HETATM"):
                parts = re.split(r'\s+', line.strip())
                x = float(parts[5])
                y = float(parts[6])
                z = float(parts[7])
                current_coordinates.append((x, y, z))
    return models_avg_coordinates

def get_model_pocket(pocket_data_df, models_avg_coordinates):
    """
    Assignes pocket to each model based on the proximity to the ligand.

    :param pocket_data_df: A dataframe generated from a CSV file containing PDB predictions.
    :param models_avg_coordinates: A list where on the i-th index are the average coordinates of the ligand from the (i+1)-th model.
    :return: A list where on the i-th index is the name of a protein pocket closest to the ligand from the (i+1)-th model.
    """
    model_pocket = []
    for avg_coordinates in models_avg_coordinates:
        avg_x, avg_y, avg_z = avg_coordinates
        min_distance = float('inf')

        for _, row in pocket_data_df.iterrows():
            pocket_name = row["name"].strip()
            center_x = row["center_x"]
            center_y = row["center_y"]
            center_z = row["center_z"]
            distance = np.sqrt((avg_x - center_x) ** 2 + (avg_y - center_y) ** 2 + (avg_z - center_z) ** 2)

            if distance < min_distance:
                min_distance = distance
                closest_pocket = pocket_name
            
        model_pocket.append(closest_pocket)
    return model_pocket


def count_hydrogen_bonds(structure_filepath, out_vina_filepath, ligand_filepath, model_pocket, pocket_residues_dict):
    """
    Calculates the number of hydrogen bonds between protein pockets and ligands.
    
    :param structure_filepath: Path to PDBTQ file containing the structure info.
    :param out_vina_filepath: Path to PDBTQ file containing the docking info.
    :param ligand_filepath: Path to PDBTQ file containing the ligand info.
    :param model_pocket: A list where on the i-th index is the name of a protein pocket closest to the ligand from the (i+1)-th model.
    :param pocket_residues_dict: A dictionary where keys are pocket names and values are residues joined in one string 
    formatted for the pymol selection command.
    :return: The number of hydrogen bonds between protein pockets and ligands.
    """
    pymol.finish_launching(['pymol', '-qc'])

    cmd.load(structure_filepath, "structure")  
    cmd.load(out_vina_filepath, "out_vina")    
    cmd.load(ligand_filepath, "ligand")  
    
    cmd.h_add()

    model_order = 0
    for state in range(1, cmd.count_states("out_vina") + 1):
        pocket_name = model_pocket[model_order]
        
        cmd.select(pocket_name, pocket_residues_dict[pocket_name])
        cmd.frame(state)

        cmd.select("h_bonds", "pocket within 3.5 of out_vina")

        num_bonds = cmd.count_atoms("h_bonds")

        print(f"State {state}: {num_bonds}")

    cmd.quit()

if __name__=="__main__":
    predictions_filepath = "/Users/alexbotkova/analysis_of_docking_studies/electrostatics/structure.pdb_predictions.csv"
    out_vina_filepath = "hydrogen_bonds/out_vina.pdbqt"
    structure_filepath = "hydrogen_bonds/structure.pdbqt"
    ligand_filepath = "hydrogen_bonds/ligand.pdbqt"

    pocket_data_df = get_pocket_data_df(predictions_filepath)
    pocket_residues_dict = get_pocket_residues_dict(pocket_data_df)
    models_avg_coordinates = get_models_avg_coordinates(out_vina_filepath)
    model_pocket = get_model_pocket(pocket_data_df, models_avg_coordinates)
    count_hydrogen_bonds(structure_filepath, out_vina_filepath, ligand_filepath, model_pocket, pocket_residues_dict)