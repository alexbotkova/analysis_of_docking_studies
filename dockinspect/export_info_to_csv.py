import csv
from ligand import Ligand
from pockets import Pockets
from poses import Poses

def export_ligand_info(ligand: Ligand, csv_path: str):
    """
    Exports ligand physicochemical properties to a CSV file.

    :param ligand: Ligand object with computed properties.
    :param csv_path: Path to the output CSV file.
    :return: None
    """
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Property", "Value"])
        writer.writerow(["SMILES", ligand.smiles])
        writer.writerow(["logP", f"{ligand.logp:.2f}"])
        writer.writerow(["SASA (Å²)", f"{ligand.sasa:.2f}"])
        writer.writerow(["TPSA (Å²)", f"{ligand.tpsa:.2f}"])
        writer.writerow(["Volume (Å³)", f"{ligand.volume:.2f}"])
        writer.writerow(["Charge", f"{ligand.charge:+d}"])

def export_pocket_info(pockets: Pockets, csv_path: str):
    """
    Exports pocket metrics (SASA, GRAVY, charge) to a CSV file.

    :param pockets: Pockets object with computed pocket metrics.
    :param csv_path: Path to the output CSV file.
    :return: None
    """
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Pocket ID", "SASA (Å²)", "GRAVY", "Charge"])
        for pid in sorted(pockets.pocket_sasas.keys()):
            row = pockets.format_pocket_row(pid).split()
            writer.writerow(row)

def export_poses_info(poses: Poses, csv_path: str):
    """
    Exports docking pose–pocket interaction summaries to a CSV file.

    :param poses: Poses object containing pose-level metrics and pocket assignments.
    :param csv_path: Path to the output CSV file.
    :return: None
    """
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Pose", "Pocket", "HBonds", "GRAVY/LogP", "SASA(P/L/R)", "Charge(P/L)"])
        for i in range(poses.number_of_models):
            row = poses.format_pose_row(i).split()
            writer.writerow(row)

def export_hbond_residues(poses: Poses, csv_path: str):
    """
    Exports hydrogen bond residue names per pose to a CSV file.

    :param poses: Poses object containing model_hbonds (a list of residue names per pose).
    :param csv_path: Path to the output CSV file.
    :return: None
    """
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Pose", "Residue"])
        for i, hbonds_res in enumerate(poses.model_hbonds):
            pose_number = i + 1
            for res in hbonds_res:
                res_str = res if isinstance(res, str) else f"{res[1]}-{res[0]}{res[2]}"
                writer.writerow([pose_number, res_str])
