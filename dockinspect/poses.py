"""
Analyzes docking poses to associate ligands with predicted protein pockets and evaluate interactions.

This script uses PyMOL to:
- Compute hydrogen bonds for each pose.
- Match each ligand pose to its closest predicted pocket.
- Combine ligand and pocket descriptors (SASA, GRAVY, charge) for pose-level analysis.

Classes:
    Poses: Evaluates ligand-pocket interactions and formats pose-specific properties.
"""

import re
import numpy as np
from pandas import DataFrame
from dockinspect.ligand import Ligand
from dockinspect.pockets import Pockets
from dockinspect.my_csv_parser import get_df
from dockinspect.hbonds import get_hydrogen_bonds_and_residues

class Poses:
    """
    Represents a collection of docking poses and analyzes their interaction with protein pockets.

        ligand (Ligand): Ligand object.
        pockets (Pockets): Pockets object.
        model_hbonds (list): List of hydrogen bond pairs per model.
        model_pockets (list): List of closest pockets assigned to each model.
        number_of_models (int): Total number of docking poses (models).
    """
    @staticmethod
    def get_models_avg_coordinates(vina_file: str) -> list:
        """
        Computes the average coordinates of ligand atoms for each docking pose.

        :param vina_file: Path to a Vina output PDBQT file.
        :return: List of average (x, y, z) coordinates for each pose.
        """
        models_avg_coordinates = []
        current_coordinates = []
        with open(vina_file, 'r') as file:
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
    def get_model_pocket_helper(pocket_data_df: DataFrame, models_avg_coordinates: list) -> list:
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
    def get_model_pocket(predictions_file: str, vina_file: str) -> list:
        """
        Wrapper for assigning closest pockets to each pose based on coordinates.

        :param predictions_filepath: Path to CSV with pocket centers.
        :param out_vina_filepath: Path to Autodock Vina PDBQT output file.
        :return: List of closest pockets per pose.
        """
        pocket_data_df = get_df(predictions_file)
        models_avg_coordinates = Poses.get_models_avg_coordinates(vina_file)
        return Poses.get_model_pocket_helper(pocket_data_df, models_avg_coordinates)

    def __init__(self, ligand: Ligand, pockets: Pockets, pdb_code: str, vina_file: str, predictions_file: str, distance:float, angle:float, model_hbonds_res:list = None):
        """
        Initializes a Poses object by associating each ligand pose with a pocket and computing H-bonds.

        :param ligand: Ligand object.
        :param pockets: Pockets object.
        :param pdb_code: PDB code of the target protein.
        :param vina_file: Path to PDBQT AutoDock Vina output file.
        :param predictions_file: Path to CSV file with pocket predictions.
        :param distance: Distance cutoff for H-bond detection.
        :param angle: Angle cutoff for H-bond detection.
        :úaram model_hbonds_res: Precomputed list of residues participating in H-bonds.
        """
        self.ligand = ligand
        self.pockets = pockets
        if model_hbonds_res is None:
            self.model_hbonds_res = get_hydrogen_bonds_and_residues(pdb_code, vina_file, distance, angle)[1]
        else:
            self.model_hbonds_res = model_hbonds_res
        self.model_pockets = Poses.get_model_pocket(predictions_file, vina_file)
        self.number_of_models = len(self.model_pockets)

    def format_pose_row(self, i: int) -> str:
        """
        Formats a single row of pose-pair analysis data.

        :param i: Index of the docking pose.
        :return: Formatted string of pose properties and interactions.
        """
        pocket = self.model_pockets[i]
        hbonds_count = len(self.model_hbonds_res[i])
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

    def __str__(self) -> str:
        """
        Returns a formatted string table summarizing ligand–pocket interaction properties for each pose.

        :return: Multiline string table of pose-level metrics.
        """
        header = f"{'Pose':<6} {'Pocket':<10} {'HBonds':<8} {'GRAVY(P)/LogP(L)':<15} {'SASA(P/L/R)':<25} {'Charge(P/L)':<15}"
        lines = [header, "-" * len(header)]

        for i in range(self.number_of_models):
            lines.append(self.format_pose_row(i))
        return "\n".join(lines)