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
from export_info_to_csv import export_ligand_info, export_pocket_info, export_poses_info, export_hbond_residues
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

        Usage:
            ligand_info [--csv FILE]

        :param line: Optional argument --csv FILE to export output.
        """
        if not self.session.ligand:
            print("No ligand loaded.")
            return

        tokens = shlex.split(line or "")
        if "--csv" in tokens:
            try:
                csv_idx = tokens.index("--csv")
                file_path = tokens[csv_idx + 1]
                export_ligand_info(self.session.ligand, file_path)
                print(f"Ligand info saved to {file_path}")
            except Exception as e:
                print(f"Failed to write CSV: {e}")
        else:
            print(self.session.ligand)


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
            pocket_info [POCKET_ID] [--csv FILE]

        If POCKET_ID is given, shows just that one.
        If --csv is provided, saves all pockets to file.
        """
        if not self.session.pockets:
            print("Pocket data not available.")
            return

        tokens = shlex.split(line or "")
        pid = None
        csv_path = None

        i = 0
        while i < len(tokens):
            if tokens[i] == "--csv":
                if i + 1 < len(tokens):
                    csv_path = tokens[i + 1]
                    i += 2
                else:
                    print("Missing filename after --csv")
                    return
            else:
                pid = tokens[i]
                i += 1

        header = f"{'Pocket ID':<10} {'SASA (Å²)':>10} {'GRAVY':>10} {'Charge':>10}"
        print(header)
        print("-" * 42)

        if pid:
            if pid not in self.session.pockets.pocket_sasas:
                print(f"Pocket '{pid}' not found.")
                return
            print(self.session.pockets.format_pocket_row(pid))
        else:
            for pid in sorted(self.session.pockets.pocket_sasas.keys()):
                print(self.session.pockets.format_pocket_row(pid))

        if csv_path:
            try:
                export_pocket_info(self.session.pockets, csv_path)
                print(f"\nAll pocket info saved to {csv_path}")
            except Exception as e:
                print(f"Failed to write CSV: {e}")

    def do_poses_info(self, line):
        """
        Displays docking pose data for all or one specific pose.

        Usage:
            poses_info [POSE_INDEX] [--csv FILE] [--res] [--csv_hbonds FILE]

        If a pose index is given, only that pose is shown.
        If --csv is given, all poses are saved to a CSV.
        If --res is given, only residue names from H-bonds are shown.
        If --csv_hbonds is given, hydrogen bond residues per pose are saved to a CSV.
        """
        if not self.session.poses:
            print("Poses data not available. Ensure all required inputs are loaded.")
            return

        tokens = shlex.split(line or "")
        index = None
        csv_path = None
        show_res_only = False
        csv_hbonds_path = None

        i = 0
        csv_path = None
        csv_hbonds_path = None
        show_res_only = False

        i = 0
        while i < len(tokens):
            if tokens[i] == "--csv":
                if i + 1 < len(tokens):
                    csv_path = tokens[i + 1]
                    i += 2
                else:
                    print("Missing filename after --csv")
                    return
            elif tokens[i] == "--csv_hbonds":
                if i + 1 < len(tokens):
                    csv_hbonds_path = tokens[i + 1]
                    i += 2
                else:
                    print("Missing filename after --csv_hbonds")
                    return
            elif tokens[i] == "--res":
                show_res_only = True
                i += 1
            else:
                try:
                    index = int(tokens[i]) - 1
                except ValueError:
                    print(f"Ignoring unrecognized argument: {tokens[i]}")
                i += 1

        if show_res_only:
            print("\nHydrogen-bonding residues per pose:\n")
            for i, hbonds_res in enumerate(self.session.poses.model_hbonds):
                pose_label = f"Pose {i+1:>2}"
                if not hbonds_res:
                    print(f"{pose_label}: — No hydrogen bonds —")
                else:
                    formatted_res = [res if isinstance(res, str) else f"{res[1]}-{res[0]}{res[2]}" for res in hbonds_res]
                    print(f"{pose_label}: {', '.join(formatted_res)}")
            print()
            return

        header = f"{'Pose':<6} {'Pocket':<10} {'HBonds':<8} {'GRAVY/LogP':<15} {'SASA(P/L/R)':<25} {'Charge(P/L)':<15}"

        if index is not None:
            if not (0 <= index < self.session.poses.number_of_models):
                print("Invalid pose index.")
                return
            print(header)
            print("-" * len(header))
            print(self.session.poses.format_pose_row(index))
        else:
            print(self.session.poses)

        if csv_path:
            try:
                export_poses_info(self.session.poses, csv_path)
                print(f"\nAll poses info saved to {csv_path}")
            except Exception as e:
                print(f"Failed to write CSV: {e}")
        
        if csv_hbonds_path:
            try:
                export_hbond_residues(self.session.poses, csv_hbonds_path)
                print(f"\nHydrogen bond residue data saved to {csv_hbonds_path}")
            except Exception as e:
                print(f"Failed to write CSV: {e}")


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

@click.group()
def cli():
    """Docking analysis tool"""
    pass

cli.add_command(launch_shell, name="shell")

@cli.command(name="ligand_info")
@click.argument("ligand_smiles", type=click.STRING)
@click.option("--csv", type=click.Path(), help="Save output to CSV.")
def ligand_info(ligand_smiles, csv):
    """
    Shows physicochemical properties of a ligand from a SMILES string.

    :param ligand_smiles: SMILES string representing the ligand.
    :param csv: Optional path to save the output as a CSV file.
    :return: None
    """
    ligand = Ligand(ligand_smiles)
    if csv:
        try:
            export_ligand_info(ligand, csv)
            print(f"Ligand info saved to {csv}")
        except Exception as e:
            print(f"Failed to write CSV: {e}")
    else:
        print(ligand)

@cli.command(name="pocket_info")
@click.argument("structure_file", type=click.Path(exists=True))
@click.argument("predictions_file", type=click.Path(exists=True))
@click.argument("residues_file", type=click.Path(exists=True))
@click.option("--pocket_id", help="Show specific pocket.")
@click.option("--csv", type=click.Path(), help="Save all pockets to CSV.")
def pocket_info(structure_file, predictions_file, residues_file, pocket_id, csv):
    """
    Shows SASA, GRAVY, and charge for predicted pockets in a protein.

    :param structure_file: Path to the structure file (.pdb or .pdbqt).
    :param predictions_file: Path to the pocket predictions file (from p2rank).
    :param residues_file: Path to the pocket-to-residues mapping file.
    :param pocket_id: Optional ID of a specific pocket to display.
    :param csv: Optional path to save pocket info as a CSV file.
    :return: None
    """
    pockets = Pockets(structure_file, predictions_file, residues_file)

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
            print(f"Pocket info saved to {csv}")
        except Exception as e:
            print(f"Failed to write CSV: {e}")

@cli.command(name="poses_info")
@click.argument("ligand_smiles", type=click.STRING)
@click.argument("pdb_code", type=click.STRING)
@click.argument("vina_file", type=click.Path(exists=True))
@click.argument("structure_file", type=click.Path(exists=True))
@click.argument("predictions_file", type=click.Path(exists=True))
@click.argument("residues_file", type=click.Path(exists=True))
@click.option("--pose_index", type=int, help="Index of the pose to show (1-based).")
@click.option("--csv", type=click.Path(), help="Save all poses to CSV.")
@click.option("--res", is_flag=True, help="Show only residue names involved in H-bonds.")
@click.option("--csv_hbonds", type=click.Path(), help="Save hydrogen bond residues per pose to CSV.")
def poses_info(ligand_smiles, pdb_code, vina_file, structure_file, predictions_file, residues_file, pose_index, csv, res, csv_hbonds):
    """
    Shows ligand–pocket interaction summary per docking pose.

    :param ligand_smiles: SMILES string of the ligand.
    :param pdb_code: PDB code of the protein.
    :param vina_file: Path to the AutoDock Vina output file (.pdbqt).
    :param structure_file: Path to the protein structure file (.pdb or .pdbqt).
    :param predictions_file: Path to the pocket predictions file.
    :param residues_file: Path to the pocket-to-residues mapping file.
    :param pose_index: Optional index of a specific pose to display (1-based).
    :param csv: Optional path to save poses info as a CSV file.
    :param res: If provided only residue names from H-bonds are shown.
    :param csv_hbonds: Optional path to save hydrogen bond residues per pose are saved to a CSV.
    :return: None
    """
    ligand = Ligand(ligand_smiles)
    pockets = Pockets(structure_file, predictions_file, residues_file)
    poses = Poses(ligand, pockets, pdb_code, vina_file, predictions_file)

    header = f"{'Pose':<6} {'Pocket':<10} {'HBonds':<8} {'GRAVY/LogP':<15} {'SASA(P/L/R)':<25} {'Charge(P/L)':<15}"

    if pose_index:
        index = pose_index - 1
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
    
    if res:
        print("\nHydrogen-bonding residues per pose:\n")
        for i, hbonds_res in enumerate(poses.model_hbonds):
            pose_label = f"Pose {i+1:>2}"
            if not hbonds_res:
                print(f"{pose_label}: — No hydrogen bonds —")
            else:
                formatted_res = [
                    res if isinstance(res, str) else f"{res[1]}-{res[0]}{res[2]}"
                    for res in hbonds_res
                ]
                print(f"{pose_label}: {', '.join(formatted_res)}")
        print()
        return
    
    if csv_hbonds:
        try:
            export_hbond_residues(poses, csv_hbonds)
            print(f"Hydrogen bond residue data saved to {csv_hbonds}")
        except Exception as e:
            print(f"Failed to write CSV: {e}")


if __name__ == "__main__":
    cli() 