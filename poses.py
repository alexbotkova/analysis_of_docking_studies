"""
Analyzes docking poses to associate ligands with predicted protein pockets and evaluate interactions.

This script uses PyMOL to:
- Compute hydrogen bonds for each pose.
- Match each ligand pose to its closest predicted pocket.
- Combine ligand and pocket descriptors (SASA, GRAVY, charge) for pose-level analysis.

Classes:
    Poses: Evaluates pose-pocket interactions and formats pose-specific properties.
"""

import pymol
from pymol import cmd
from my_tools.my_parser import *
from pockets import Pockets
from ligand import Ligand
import re
import numpy as np
import sys
sys.path.append("..")
from my_tools.my_parser import *

class Poses:
    """
    Represents a collection of ligand poses and analyzes their interaction with protein pockets.

    Attributes:
        ligand (Ligand): Ligand object containing molecular descriptors.
        pockets (Pockets): Pockets object with pocket-level SASA, GRAVY, and charge info.
        model_hbonds (list): List of hydrogen bond pairs per model.
        model_pockets (list): List of closest pockets assigned to each model.
        number_of_models (int): Total number of docking poses (models).
    """

    @staticmethod
    def get_hydrogen_bonds(structure_filepath, out_vina_filepath):
        """
        Calculates the number of hydrogen bonds between protein and ligand for each pose.

        :param structure_filepath: Path to the structure PDBQT file.
        :param out_vina_filepath: Path to the Vina docking output file (PDBQT).
        :return: A list of hydrogen bond pairs (tuples of atom indices) for each pose.
        """
        pymol.finish_launching(['pymol', '-qc'])

        cmd.load(structure_filepath, "structure")  
        cmd.load(out_vina_filepath, "out_vina")    
        cmd.h_add()

        model_hbonds = []
        for state in range(1, cmd.count_states("out_vina") + 1):
            cmd.frame(state)
            sel1 = 'structure and (donor or acceptor)'
            sel2 = 'out_vina and (donor or acceptor)'
            hbonds = cmd.find_pairs(sel1, sel2, mode=1, cutoff=3.5, angle=55.0)
            model_hbonds.append(hbonds)

        cmd.quit()
        return model_hbonds

    @staticmethod
    def get_models_avg_coordinates(out_vina_filepath):
        """
        Computes the average coordinates of ligand atoms for each docking pose.

        :param out_vina_filepath: Path to a Vina output PDBQT file.
        :return: List of average (x, y, z) coordinates for each pose.
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
            models_avg_coordinates.append(np.mean(current_coordinates, axis=0))
        return models_avg_coordinates

    @staticmethod
    def get_model_pocket_helper(pocket_data_df, models_avg_coordinates):
        """
        Assigns the closest pocket to each ligand pose based on spatial proximity.

        :param pocket_data_df: DataFrame from pocket prediction CSV.
        :param models_avg_coordinates: List of ligand center coordinates for each pose.
        :return: List of closest pocket names for each model.
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

    @staticmethod
    def get_model_pocket(predictions_filepath, out_vina_filepath):
        """
        Wrapper for assigning closest pockets to each pose based on coordinates.

        :param predictions_filepath: Path to CSV with pocket centers.
        :param out_vina_filepath: Path to docking output (PDBQT).
        :return: List of closest pockets per pose.
        """
        pocket_data_df = get_df(predictions_filepath)
        models_avg_coordinates = Poses.get_models_avg_coordinates(out_vina_filepath)
        return Poses.get_model_pocket_helper(pocket_data_df, models_avg_coordinates)

    def __init__(self, ligand, pockets, structure_filepath, out_vina_filepath, predictions_filepath):
        """
        Initializes a Poses object by associating each ligand pose with a pocket and computing HBonds.

        :param ligand: Ligand object.
        :param pockets: Pockets object.
        :param structure_filepath: Path to structure file used in docking.
        :param out_vina_filepath: Path to Vina docking output file.
        :param predictions_filepath: Path to CSV with pocket predictions.
        """
        self.ligand = ligand
        self.pockets = pockets

        self.model_hbonds = Poses.get_hydrogen_bonds(structure_filepath, out_vina_filepath)
        self.model_pockets = Poses.get_model_pocket(predictions_filepath, out_vina_filepath)
        self.number_of_models = len(self.model_pockets)

    def format_pose_row(self, i):
        """
        Formats a single row of pose-pair analysis data.

        :param i: Index of the docking pose.
        :return: Formatted string of pose properties and interactions.
        """
        pocket = self.model_pockets[i]
        hbonds_count = len(self.model_hbonds[i])
        gravy = self.pockets.pocket_gravys.get(pocket, 0)
        sasa_pocket = self.pockets.pocket_sasas.get(pocket, 0)
        sasa_ligand = self.ligand.sasa
        sasa_ratio = sasa_ligand / sasa_pocket if sasa_pocket else 0
        charge_pocket = self.pockets.pocket_charges.get(pocket, 0)
        charge_ligand = self.ligand.charge

        return (
            f"{i + 1:<6} {pocket:<10} {hbonds_count:<8} "
            f"{f'{gravy:.2f}/{self.ligand.logp:.2f}':<15} "
            f"{f'{sasa_pocket:.2f}/{sasa_ligand:.2f}/{sasa_ratio:.2f}':<25} "
            f"{f'{charge_pocket}/{charge_ligand}':<15}"
        )

    def __str__(self):
        """
        Returns a formatted string table summarizing ligand–pocket interaction properties for each pose.

        :return: Multiline string table of pose-level metrics.
        """
        header = f"{'Pose':<6} {'Pocket':<10} {'HBonds':<8} {'GRAVY(P)/LogP(L)':<15} {'SASA(P/L/R)':<25} {'Charge(P/L)':<15}"
        lines = [header, "-" * len(header)]

        for i in range(self.number_of_models):
            lines.append(self.format_pose_row(i))
        return "\n".join(lines)

if __name__ == "__main__":
    """
    Example usage to analyze docking poses and pocket interactions using a SMILES string and prediction files.
    """
    smiles = "NC(N)=O"
    structure_filepath = "/Users/alexbotkova/analysis_of_docking_studies/test_files/urea/result-2025-03-30T21_20_58.783Z/structure.pdbqt"
    out_vina_filepath = "/Users/alexbotkova/analysis_of_docking_studies/test_files/urea/result-2025-03-30T21_20_58.783Z/out_vina.pdbqt"
    predictions_filepath = "/Users/alexbotkova/analysis_of_docking_studies/test_files/urea/prankweb-2SRC/structure.cif_predictions.csv"
    residues_filepath = "/Users/alexbotkova/analysis_of_docking_studies/test_files/urea/prankweb-2SRC/structure.cif_residues.csv"

    ligand = Ligand(smiles)
    pockets = Pockets(structure_filepath, predictions_filepath, residues_filepath)
    o = Poses(ligand, pockets, structure_filepath, out_vina_filepath, predictions_filepath)
    print(o)
