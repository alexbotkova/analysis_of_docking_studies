"""
Analyzes predicted protein pockets to compute SASA, GRAVY, and net charge.

This script uses FreeSASA and residue-level annotations to compute:
- Solvent Accessible Surface Area (SASA) for each pocket
- GRAVY of residues in each pocket
- Net formal charge of each pocket based on amino acid composition

Classes:
    Pockets: Encapsulates methods and data for analyzing protein pocket properties.
"""

import freesasa
from dockinspect.my_csv_parser import get_df, parse_predictions, parse_residues
from dockinspect.pymol_tools import get_pocket_residues_dict

freesasa.setVerbosity(1)

class Pockets:
    """
    Represents a collection of predicted protein pockets and calculates their physicochemical properties.

    Attributes:
        pocket_charges (dict): Net formal charge of each pocket.
        pocket_gravys (dict): GRAVY index (hydropathy) of each pocket.
        pocket_sasas (dict): Solvent-accessible surface area of each pocket.
    """

    @staticmethod
    def get_pocket_charge(pocket_locations_dict: dict, location_aa_dict: dict) -> dict:
        """
        Calculates the charge of each protein pocket based on the predicted amino acids.

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
                if aa in ["ARG", "LYS", "HIS"]:
                    pocket_charge_dict[pocket] += 1
                elif aa in ["ASP", "GLU"]: 
                    pocket_charge_dict[pocket] -= 1
        return pocket_charge_dict
     
    @staticmethod
    def get_pocket_gravy(pocket_locations_dict: dict, location_aa_dict: dict) -> dict:
        """
        Calculates the GRAVY (grand average of hydropathy) of each protein pocket based on the predicted amino acids.

        :param pocket_locations_dict: Dictionary where keys are pocket IDs and values are lists of chain and residue locations.
        :param location_aa_dict: Dictionary mapping chain-residue locations to their corresponding amino acid abbreviations.
        :return: Dictionary where keys are pocket IDs and values are the average hydropathy of the respective pockets.
        """
        hydrophobicity_scale = {
            "ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5,
            "GLN": -3.5, "GLU": -3.5, "GLY": -0.4, "HIS": -3.2, "ILE": 4.5,
            "LEU": 3.8, "LYS": -3.9, "MET": 1.9, "PHE": 2.8, "PRO": -1.6,
            "SER": -0.8, "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2
        }

        pocket_gravy_dict = {}
        for pocket in pocket_locations_dict:
            pocket_gravy_dict[pocket] = 0
            locations = pocket_locations_dict[pocket]
            for location in locations:
                aa = location_aa_dict[location]
                if aa:
                    pocket_gravy_dict[pocket] += hydrophobicity_scale[aa]
            pocket_gravy_dict[pocket] /= len(locations) if locations else 0
            pocket_gravy_dict[pocket] = round(pocket_gravy_dict[pocket], 2)
        return pocket_gravy_dict

    @staticmethod
    def get_pocket_sasas(predictions_file: str, structure_file: str) -> dict:
        """
        Calculates solvent accessible surface area (SASA) for each pocket using FreeSASA.

        :param predictions_file: Path to the CSV file containing pocket residue predictions.
        :param structure_file: Path to the structure file in PDBQT format.
        :return: Dictionary where keys are pocket IDs and values are their SASA values (A^2).
        """
        pocket_data_df = get_df(predictions_file)
        pocket_residues_dict = get_pocket_residues_dict(pocket_data_df)
        commands =  [str(key) + ", " + str(value) for key, value in pocket_residues_dict.items()]

        structure = freesasa.Structure(structure_file)
        result = freesasa.calc(structure)
        return freesasa.selectArea(commands, structure, result)

    def __init__(self, structure_file: str, predictions_file: str, residues_file: str) -> None:
        """
        Initializes the Pockets object by computing pocket charge, GRAVY and SASA.

        :param structure_file: Path to the PDBQT structure file.
        :param predictions_file: Path to the predicted pocket CSV file.
        :param residues_file: Path to the residues CSV file.
        """
        pocket_locations_dict = parse_predictions(predictions_file)
        location_aa_dict = parse_residues(residues_file)
        
        self.pocket_charges = Pockets.get_pocket_charge(pocket_locations_dict, location_aa_dict)
        self.pocket_gravys = Pockets.get_pocket_gravy(pocket_locations_dict, location_aa_dict)
        self.pocket_sasas = Pockets.get_pocket_sasas(predictions_file, structure_file)
    
    def format_pocket_row(self, pid: str) -> str:
        """
        Formats a single row of pocket data for printing.

        :param pid: Pocket ID.
        :return: Formatted string with SASA, GRAVY and charge values for the given pocket.
        """
        sasa = self.pocket_sasas.get(pid, 0)
        gravy = self.pocket_gravys.get(pid, 0)
        charge = self.pocket_charges.get(pid, 0)
        return f"{str(pid):<10} {sasa:>10.2f} {gravy:>10.2f} {charge:>10}"

    def __str__(self) -> str:
        """
        Returns a formatted string representation of all pocket properties (SASA, GRAVY, Charge).

        :return: Multiline string table of all pocket descriptors.
        """
        pocket_ids = sorted(set(self.pocket_charges.keys()) |
                            set(self.pocket_gravys.keys()) |
                            set(self.pocket_sasas.keys()))

        lines = []
        lines.append(f"{'Pocket ID':<10} {'SASA (A^2)':>10} {'GRAVY':>10} {'Charge':>10}")
        lines.append("-" * 42)

        for pid in pocket_ids:
            lines.append(self.format_pocket_row(pid))
        return "\n".join(lines)