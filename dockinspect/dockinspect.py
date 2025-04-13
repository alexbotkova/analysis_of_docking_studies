"""
Command-line interface for interactive analysis and visualization of protein–ligand docking results.

This script provides a shell-based interface that allows the user to:
- Load a ligand from a SMILES string
- Load docking results from AutoDock Vina
- Visualize docking poses and interactions in PyMOL
- Display physicochemical properties of ligand and pockets
- Analyze hydrogen bonding and pose–pocket assignments

Commands:
    ligand_info     - Show ligand properties (LogP, SASA, TPSA, Volume, Charge)
    visualize       - Launch PyMOL with structure, pose, and selected visualization mode
    pocket_info     - Show SASA, GRAVY, and Charge for one or all pockets
    poses_info      - Show ligand–pocket interaction summary per pose
    exit            - Exit the shell

Run via CLI:
    python script.py <SMILES> <vina_file> <structure_file> [<predictions_file> <residues_file>]
"""

import click
import cmd as shell_cmd
import shlex
import os
import subprocess


from ligand import Ligand
from pockets import Pockets
from poses import Poses
from visualization import *
from my_parser import *
from typing import Optional

class Session:
    """
    Represents the state of a docking analysis session.

    :param smiles: SMILES string for the ligand.
    :param vina_file: Path to the out_vina file from AutoDock Vina.
    :param pdb_code: Protein structure PDB code.
    :param predictions_file: Path to the file containing predicted pockets from p2rank (optional).
    :param residues_file: Path to the file mapping pockets to residues from p2rank (optional).
    """
    def __init__(self, smiles=None, vina_file=None, pdb_code=None, structure_file=None, predictions_file=None, residues_file=None):
        self.ligand = Ligand(smiles) if smiles else None
        self.out_vina_file = vina_file
        self.pdb_code = pdb_code
        self.structure_file = structure_file
        self.predictions_file = predictions_file
        self.residues_file = residues_file
        self.pockets = None
        self.poses = None 

def launch_pymol_visualization(pdb_code: str, vina_file: str, pose_num: int = 1, mode: Optional[str] = "", pocket_selection: Optional[str] = ""):
    """
    Launches PyMOL with a visualization script for binding poses from AutoDock Vina.

    :param pdb_code: PDB code for the protein structure to fetch remotely.
    :param vina_file: Path to the out_vina file from AutoDock Vina.
    :param pose_num: Pose number to visualize (1-based index, default is 1).
    :param mode: Visualization mode ('surface', 'polar', 'charge', 'hbonds', default is 'hbonds').
    :param pocket_selection: PyMOL atom selection string for highlighting a pocket (from p2rank predictions).
    :return: None
    """
    script_path = os.path.join(os.path.dirname(__file__), "visualization_runtime.py")
    project_root = os.path.abspath(os.path.dirname(__file__))

    with open(script_path, "w") as f:
        f.write(f"""from pymol import cmd
import sys
import os
sys.path.insert(0, {repr(project_root)})

from my_parser import *
from visualization import *

visualize(
    pdb_code={repr(pdb_code)},
    vina_file={repr(vina_file)},
    pocket_selection={repr(pocket_selection)},
    pose_num={pose_num},
    mode={repr(mode)}
)
""")

    try:
        env = os.environ.copy()
        env["PYTHONPATH"] = project_root
        subprocess.run(["pymol", script_path], cwd=project_root, env=env, stdout=subprocess.DEVNULL)
    except Exception as e:
        print(f"Error running PyMOL: {e}")


class Shell(shell_cmd.Cmd):
    """
    Interactive shell interface for exploring ligand–protein docking results.

    :param session: A Session object containing all loaded data.
    """
    intro = "Type 'help' or '?' to list commands.\n"
    prompt = '> '

    def __init__(self, session: Session):
        super().__init__()
        self.session = session

    def do_ligand_info(self, line=None):
        """
        Displays information about the loaded ligand.

        :param line: Ignored; included for cmd compatibility.
        """
        if self.session.ligand:
            print(self.session.ligand)
        else:
            print("No ligand loaded.")

    def do_visualize(self, line):
        """
        Launches PyMOL visualization for the selected pose and mode.

        Usage:
            visualize [--pose NUM] [--mode MODE]

        :param line: Command-line string with optional arguments:
                     --pose (int): Pose number to visualize (default = 1).
                     --mode (str): Visualization mode, one of {'surface', 'polar', 'charge', 'hbonds'} (default = "hbonds").
        """
        pose_num = 1
        mode = ""
        tokens = shlex.split(line)

        for i in range(len(tokens)):
            if tokens[i] == "--pose" and i + 1 < len(tokens):
                pose_num = int(tokens[i + 1])
            if tokens[i] == "--mode" and i + 1 < len(tokens):
                mode = tokens[i + 1]

        valid_modes = {"surface", "polar", "charge", "hbonds", ""}
        if mode not in valid_modes:
            print(f"Invalid mode '{mode}'. Choose from: {', '.join(valid_modes - {''})}")
            return

        if self.session.structure_file and self.session.out_vina_file:
            print(f"Launching PyMOL with mode '{mode or 'default'}'...")
            pocket_selection = ""

            if self.session.poses:
                try:
                    model_pockets = self.session.poses.model_pockets 
                    pocket_residues_dict = get_pocket_atomids_dict(get_df(self.session.predictions_file))
                    pocket_id = model_pockets[pose_num - 1]  
                    pocket_selection = pocket_residues_dict.get(pocket_id, "")
                    #print(f"Using pocket: {pocket_id} with selection: {pocket_selection}")
                except Exception as e:
                    print(f"Warning: Could not determine pocket selection: {e}")

            launch_pymol_visualization(
                pdb_code=self.session.pdb_code,
                vina_file=self.session.out_vina_file,
                pose_num=pose_num,
                mode=mode,
                pocket_selection=pocket_selection
            )
        else:
            print("You must load files first.")

    def do_pocket_info(self, line):
        """
        Displays pocket properties: SASA, GRAVY and charge.

        Usage:
            pocket_info [POCKET_ID]

        :param line: Pocket ID to show (optional). If not provided, all pockets are displayed.
        """
        if not self.session.pockets:
            print("Pocket data not available. Provide predictions and residues files to enable this feature.")
            return

        tokens = shlex.split(line)
        header = f"{'Pocket ID':<10} {'SASA (Å²)':>10} {'GRAVY':>10} {'Charge':>10}"

        if not tokens:
            print(self.session.pockets)
            return

        pid = tokens[0]
        if pid not in self.session.pockets.pocket_sasas:
            print(f"Pocket '{pid}' not found.")
            return

        print(header)
        print("-" * 42)
        print(self.session.pockets.format_pocket_row(pid))

    def do_poses_info(self, line):
        """
        Displays docking pose data for all or one specific pose.

        Usage:
            poses_info [POSE_INDEX]

        :param line: Pose index to show (1-based, optional). If omitted, shows all poses.
        """
        if not self.session.poses:
            print("Poses data not available. Ensure all required inputs are loaded.")
            return

        tokens = shlex.split(line)
        header = f"{'Pose':<6} {'Pocket':<10} {'HBonds':<8} {'GRAVY/LogP':<15} {'SASA(P/L/R)':<25} {'Charge(P/L)':<15}"

        if not tokens:
            print(self.session.poses)
            return

        try:
            index = int(tokens[0]) - 1
            if not (0 <= index < self.session.poses.number_of_models):
                raise IndexError
            print(header)
            print("-" * len(header))
            print(self.session.poses.format_pose_row(index))
        except (ValueError, IndexError):
            print("Invalid pose index.")

    def do_exit(self, line=None):
        """
        Exits the interactive shell.

        :param line: Unused.
        :return: True to signal shell exit.
        """
        return True

    def do_EOF(self, line):
        """
        Handles Ctrl+D to exit the shell.

        :param line: Unused.
        :return: True to signal shell exit.
        """
        return self.do_exit(line)


@click.command()
@click.argument("ligand_smiles", type=click.STRING)
@click.argument("pdb_code", type=click.STRING)
@click.argument("vina_file", type=click.Path(exists=True))
@click.argument("structure_file", type=click.Path(exists=True))
@click.argument("predictions_file", required=False, type=click.Path(exists=True))
@click.argument("residues_file", required=False, type=click.Path(exists=True))
def launch_shell(ligand_smiles, pdb_code, vina_file, structure_file, predictions_file, residues_file):
    """
    Initializes a docking analysis session and launches the interactive shell.

    :param ligand_smiles: SMILES string of the ligand.
    :param vina_file: Output file from AutoDock Vina.
    :param pdb_code: Protein structure PDB code.
    :param predictions_file: Pocket prediction file (optional).
    :param residues_file: Residue annotation file (optional).
    :return: None
    """
    session = Session(
        smiles=ligand_smiles,
        pdb_code=pdb_code,
        vina_file=vina_file,
        structure_file=structure_file,
        predictions_file=predictions_file,
        residues_file=residues_file
    )

    if predictions_file and residues_file:
        try:
            session.pockets = Pockets(structure_file, predictions_file, residues_file)
        except Exception as e:
            print(f"Warning: Failed to initialize pocket data: {e}")
    else:
        print("Pocket prediction files not provided. 'pocket_info' and 'poses_info' will be disabled.")

    if session.ligand and session.pockets and vina_file and pdb_code and predictions_file:
        try:
            session.poses = Poses(
                ligand=session.ligand,
                pockets=session.pockets,
                pdb_code=session.pdb_code,
                out_vina_filepath=vina_file,
                predictions_filepath=predictions_file
            )
        except Exception as e:
            print(f"Warning: Failed to initialize poses: {e}")

    print("Session initialized. Starting interactive mode...\n")
    Shell(session).cmdloop()

if __name__ == "__main__":
    launch_shell()