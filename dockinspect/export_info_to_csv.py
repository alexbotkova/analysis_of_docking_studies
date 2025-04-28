"""
This module provides utility functions for exporting docking data to CSV files as part of the 
DockInspect toolkit. Each function corresponds to a specific object type in the docking pipeline 
(ligand, pocket, or pose) and outputs computed properties to CSV files.

Functions:
    export_ligand_info: Saves computed ligand properties
    export_pocket_info: Saves computed descriptors for predicted pockets
    export_poses_info: Saves docking pose-level properties
    export_hbond_residues: Saves pose hydrogen bond residue details
"""

import csv
from dockinspect.ligand import Ligand
from dockinspect.pockets import Pockets
from dockinspect.poses import Poses

def export_ligand_info(ligand: Ligand, csv_path: str) -> None:
    """
    Exports ligand properties to a CSV file.

    :param ligand: Ligand object.
    :param csv_path: Path to the output CSV file.
    :return: None
    """
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Property", "Value"])
        writer.writerow(["SMILES", ligand.smiles])
        writer.writerow(["logP", f"{ligand.logp:.2f}"])
        writer.writerow(["SASA", f"{ligand.sasa:.2f}"])
        writer.writerow(["TPSA", f"{ligand.tpsa:.2f}"])
        writer.writerow(["Volume", f"{ligand.volume:.2f}"])
        writer.writerow(["Charge", f"{ligand.charge}"])

def export_pocket_info(pockets: Pockets, csv_path: str) -> None:
    """
    Exports pocket properties to a CSV file.

    :param pockets: Pockets object.
    :param csv_path: Path to the output CSV file.
    :return: None
    """
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Pocket_ID", "SASA", "GRAVY", "Charge"])
        for pid in sorted(pockets.pocket_sasas.keys()):
            row = pockets.format_pocket_row(pid).split()
            writer.writerow(row)

def export_poses_info(poses: Poses, csv_path: str) -> None:
    """
    Exports docking pose properties to a CSV file.

    :param poses: Poses object.
    :param csv_path: Path to the output CSV file.
    :return: None
    """
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Pose", "Pocket", "HBonds", "GRAVY/LogP", "SASA(P/L/R)", "Charge(P/L)"])
        for i in range(poses.number_of_models):
            row = poses.format_pose_row(i).split()
            writer.writerow(row)

def export_hbond_residues(poses: Poses, csv_path: str) -> None:
    """
    Exports hydrogen bond residue names per pose to a CSV file.

    :param poses: Poses object.
    :param csv_path: Path to the output CSV file.
    :return: None
    """
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Pose", "Residue"])
        for i, hbonds_res in enumerate(poses.model_hbonds_res):
            pose_number = i + 1
            for res in hbonds_res:
                res_str = res if isinstance(res, str) else f"{res[1]}-{res[0]}{res[2]}"
                writer.writerow([pose_number, res_str])