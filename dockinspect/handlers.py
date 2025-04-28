"""
This module provides command-handling functions for the DockInspect CLI.

Functions:
    handle_ligand_info: Outputs or exports ligand properties
    handle_pocket_info: Outputs or exports pocket properties
    handle_poses_info: Outputs or exports docking pose information and optional H-bond residues
    launch_pymol_visualization: Runs PyMOL to visualize a specific pose with optional pocket highlights
    handle_visualization: Prepares visualization and pocket selection
"""

import os
import subprocess
from dockinspect.ligand import Ligand
from dockinspect.pockets import Pockets
from dockinspect.poses import Poses
from dockinspect.export_info_to_csv import export_ligand_info, export_pocket_info, export_poses_info, export_hbond_residues
from dockinspect.pymol_tools import get_pocket_residues_dict
from dockinspect.my_csv_parser import get_df
from dockinspect.visualization import save_pml
from typing import Optional

def handle_ligand_info(ligand: Ligand, csv: Optional[str] = None) -> None:
    """
    Displays ligand information or exports it to a CSV file.

    :param ligand: Ligand object.
    :param csv: Optional path to export ligand information as a CSV file.
    """
    print(ligand)

    if csv:
        try:
            export_ligand_info(ligand, csv)
            print(f"Ligand info saved to {csv}")
        except Exception as e:
            print(f"Failed to write CSV: {e}")

def handle_pocket_info(pockets: Pockets, pocket_id: Optional[str] = None, csv: Optional[str] = None) -> None:
    """
    Displays computed pocket information or exports it to a CSV file.

    :param pockets: Pockets object.
    :param pocket_id: Optional ID of a specific pocket to display.
    :param csv: Optional path to export pocket information as a CSV dile.
    """
    header = f"{'Pocket ID':<10} {'SASA (A^2)':>10} {'GRAVY':>10} {'Charge':>10}"
    print(header)
    print("-" * len(header))

    if pocket_id:
        if pocket_id not in pockets.pocket_sasas:
            print(f"Pocket '{pocket_id}' not found.")
        else:
            print(pockets.format_pocket_row(pocket_id))
    else:
        for pid in sorted(pockets.pocket_sasas.keys()):
            print(pockets.format_pocket_row(pid))

    if csv:
        try:
            export_pocket_info(pockets, csv)
            print(f"\nPocket info saved to {csv}")
        except Exception as e:
            print(f"Failed to write CSV: {e}")

def handle_poses_info(poses: Poses, pose_index: Optional[int] = None, csv: Optional[str] = None,
                      show_res_only: bool = False, csv_hbonds: Optional[str] = None) -> None:
    """
    Displays detailed information about docking poses or exports them.

    :param poses: Poses object containing all docking poses and associated properties.
    :param pose_index: Optional index to show a single pose (1-based).
    :param csv: Optional path to export pose info to a CSV file.
    :param show_res_only: If True, shows only hydrogen bonding residues per pose.
    :param csv_hbonds: Optional path to export H-bonding residue info to a CSV file.
    """
    header = f"{'Pose':<6} {'Pocket':<10} {'HBonds':<8} {'GRAVY/LogP':<15} {'SASA(P/L/R)':<25} {'Charge(P/L)':<15}"

    if show_res_only:
        print("\nHydrogen-bonding residues per pose:\n")
        for i, hbonds_res in enumerate(poses.model_hbonds_res):
            pose_label = f"Pose {i+1:>2}"
            if not hbonds_res:
                print(f"{pose_label}: No hydrogen bonds.")
            else:
                formatted_res = [res if isinstance(res, str) else f"{res[1]}-{res[0]}{res[2]}" for res in hbonds_res]
                print(f"{pose_label}: {', '.join(formatted_res)}")
        print()
        return

    if pose_index is not None:
        index = pose_index
        if not (0 <= index < poses.number_of_models):
            print("Invalid pose index.")
            return
        print(header)
        print("-" * len(header))
        print(poses.format_pose_row(index))
    else:
        print(poses)

    if csv:
        try:
            export_poses_info(poses, csv)
            print(f"Poses info saved to {csv}")
        except Exception as e:
            print(f"Failed to write CSV: {e}")
    
    if csv_hbonds:
        try:
            export_hbond_residues(poses, csv_hbonds)
            print(f"Hydrogen bond residue data saved to {csv_hbonds}")
        except Exception as e:
            print(f"Failed to write CSV: {e}")

def launch_pymol_visualization(pdb_code: str, vina_file: str, distance: float = 3.2, angle: float = 25.0, 
                               pose_num: int = 1, mode: str = "", pocket_selection: str = "") -> None:
    """
    Launches PyMOL with a visualization script for binding poses from AutoDock Vina.

    :param pdb_code: PDB code for the protein structure to fetch remotely.
    :param vina_file: Path to the out_vina file from AutoDock Vina.
    :param distance: Distance cutoff for hydrogen bond detection.
    :param angle: Angle cutoff for hydrogen bond detection.
    :param pose_num: Pose number to visualize (1-based index).
    :param mode: Visualization mode.
    :param pocket_selection: Atom selection string for highlighting pocket residues.
    """
    script_path = os.path.join(os.path.dirname(__file__), "visualization_runtime.py")
    project_root = os.path.abspath(os.path.dirname(__file__))
    fetched_filename = f"{pdb_code}.pdb" if not (pdb_code.endswith(".pdb") or pdb_code.endswith(".pdbqt")) else None

    with open(script_path, "w") as f:
        f.write(f"""from pymol import cmd
import sys
import os
sys.path.insert(0, {repr(project_root)})

from visualization import *

visualize(
    pdb_code={repr(pdb_code)},
    vina_file={repr(vina_file)},
    pocket_selection={repr(pocket_selection)},
    pose_num={pose_num},
    mode={repr(mode)},
    distance={distance},
    angle={angle}
)
""")

    try:
        env = os.environ.copy()
        env["PYTHONPATH"] = project_root
        subprocess.run(["pymol", script_path], cwd=project_root, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,)
    except Exception as e:
        print(f"Error running PyMOL: {e}")
    finally:
        try:
            os.remove(script_path)
        except OSError as cleanup_error:
            print(f"Error removing temporary script: {cleanup_error}")
        
        if fetched_filename and os.path.exists(fetched_filename):
            try:
                os.remove(fetched_filename)
                print(f"Deleted temporary fetched file {fetched_filename}")
            except Exception as e:
                print(f"Warning: Could not delete {fetched_filename}: {e}")

def handle_visualization(pdb_code: str, vina_file: str, pose_num: int, mode: str, predictions_file: Optional[str] = None, 
                        poses: Optional[Poses] = None, distance: float = 3.2, angle: float = 25.0, save: Optional[str] = None, 
                        model_hbonds = None) -> None:
    """
    Handles PyMOL visualization or PML script generation for docking poses.

    :param pdb_code: PDB code for the protein structure to fetch remotely.
    :param vina_file: Path to the PDBQT output file from AutoDock Vina.
    :param pose_num: Pose number to visualize or export (1-based index).
    :param mode: Visualization mode for PyMOL.
    :param predictions_file: Path to a CSV file containing pocket prediction results.
    :param poses: A Poses object containing model-to-pocket mappings.
    :param distance: Distance cutoff for hydrogen bond detection.
    :param angle: Angle cutoff for hydrogen bond detection.
    :param pml: If provided, path to save a generated PyMOL script instead of launching PyMOL.
    :param model_hbonds: Optional precomputed hydrogen bonds to include in the visualization.
    """
    pocket_selection = ""

    if poses and predictions_file:
        try:
            model_pockets = poses.model_pockets
            pocket_residues_dict = get_pocket_residues_dict(get_df(predictions_file))
            pocket_id = model_pockets[pose_num - 1]
            pocket_selection = pocket_residues_dict.get(pocket_id, "")
        except Exception as e:
            print(f"Warning: Could not determine pocket selection: {e}")
    
    if save:
        save_pml(pdb_code=pdb_code, vina_file=vina_file, pose_num=pose_num, pocket_selection=pocket_selection, mode=mode, model_hbonds=model_hbonds, pml=save)
        return
    else:
        launch_pymol_visualization(pdb_code=pdb_code, vina_file=vina_file, pose_num=pose_num, pocket_selection=pocket_selection, mode=mode, distance=distance, angle=angle)