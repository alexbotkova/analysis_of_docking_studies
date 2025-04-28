#!/usr/bin/env python3

"""
Command-line interface for analysis and visualization of protein–ligand docking results.

Commands:
    ligand_info     - Show ligand properties (LogP, SASA, TPSA, Volume, Charge)
    visualize       - Launch PyMOL with selected visualization mode
    pocket_info     - Show SASA, GRAVY and charge for one or all pockets
    poses_info      - Show ligand–pocket interaction summary per pose
    exit            - Exit the shell

Notes:
    The --distance and --angle options (for H-bond detection cutoffs) can be specified when running commands independently or when initializing the shell session.
    Once inside the interactive shell, --distance and --angle cannot be changed; the session uses the values given at startup.
"""

import click
import shlex
import cmd as shell_cmd
from dockinspect.ligand import Ligand
from dockinspect.pockets import Pockets
from dockinspect.poses import Poses
from dockinspect.handlers import handle_ligand_info, handle_pocket_info, handle_poses_info, handle_visualization
from dockinspect.hbonds import get_hydrogen_bonds_and_residues

class Session:
    """
    Represents the state of a docking analysis session.

    :param smiles: SMILES string of the ligand.
    :param pdb_code: Protein structure PDB code.
    :param vina_file: Path to the PDBQT out_vina file from AutoDock Vina.
    :param structure_file: Path to the PDBQT structure file from AutoDock Vina.
    :param predictions_file: Path to the CSV file containing predicted pockets (optional).
    :param residues_file: Path to the CSV file mapping pockets to residues (optional).
    :param distance: Distance cutoff for hydrogen bond detection in Angstroms.
    :param angle: Angle cutoff for hydrogen bond detection in degrees.
    :param model_hbonds: Precomputed hydrogen bonds model.
    """
    def __init__(self, smiles: str = None, pdb_code: str = None, vina_file: str = None,
                 structure_file: str = None, predictions_file: str = None, residues_file: str = None,
                 distance: float = 3.2, angle: float = 25.0, model_hbonds: list = None):
        self.ligand = Ligand(smiles) if smiles else None
        self.pdb_code = pdb_code
        self.out_vina_file = vina_file
        self.structure_file = structure_file
        self.predictions_file = predictions_file
        self.residues_file = residues_file
        self.pockets = None
        self.poses = None 
        self.distance = distance
        self.angle = angle
        self.hbonds = model_hbonds

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

    def do_ligand_info(self, line=None) -> None:
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
        csv = tokens[tokens.index("--csv") + 1] if "--csv" in tokens else None
        handle_ligand_info(self.session.ligand, csv)

    def do_visualize(self, line=None) -> None:
        """
        Launches PyMOL visualization for the selected pose and mode.

        Usage:
            visualize [--pose NUM] [--mode MODE] [--save]

        :param line: Command-line string with optional arguments:
                     --pose (int): Pose number to visualize (default = 1).
                     --mode (str): Visualization mode, one of (surface, polar, charge, hbonds) (default hbonds).
                     --save_pml FILE: Save visualization PyMOL script instead of launching PyMOL.
        """
        pose_num = 1
        mode = ""
        save = None
        tokens = shlex.split(line or "")

        i = 0
        while i < len(tokens):
            if tokens[i] == "--pose" and i + 1 < len(tokens):
                pose_num = int(tokens[i + 1])
                i += 2  
            elif tokens[i] == "--mode" and i + 1 < len(tokens):
                mode = tokens[i + 1]
                i += 2  
            elif tokens[i] == "--save" and i + 1 < len(tokens):
                save = tokens[i + 1]
                i += 2  
            else:
                print(f"Ignoring unrecognized argument: {tokens[i]}")
                i += 1

        valid_modes = {"surface", "polar", "charge", "hbonds", ""}
        if mode not in valid_modes:
            print(f"Invalid mode '{mode}'. Choose from: {', '.join(valid_modes - {''})}")
            return
        
        if self.session.structure_file and self.session.out_vina_file:
            print("Saving PyMOL script..." if save else f"Launching PyMOL with mode '{mode or 'hbonds'}'...")
            handle_visualization(pdb_code=self.session.pdb_code, vina_file=self.session.out_vina_file, pose_num=pose_num, mode=mode, predictions_file=self.session.predictions_file, poses=self.session.poses,
                distance=self.session.distance, angle=self.session.angle, save=save, model_hbonds=self.session.hbonds)
        else:
            print("You must load files first.")

    def do_pocket_info(self, line=None) -> None:
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
        pocket_id = None
        csv = None
        if "--csv" in tokens:
            idx = tokens.index("--csv")
            if idx + 1 < len(tokens):
                csv = tokens[idx + 1]
        pocket_id = next((tok for tok in tokens if not tok.startswith("--")), None)
        handle_pocket_info(self.session.pockets, pocket_id, csv)

    def do_poses_info(self, line=None) -> None:
        """
        Displays docking pose data for all or one specific pose.

        Usage:
            poses_info [POSE_INDEX] [--csv FILE] [--res] [--csv_hbonds FILE]

        If a pose index is given, only that pose is shown.
        If --csv is given, all poses are saved to a CSV.
        If --res is given, only residue names from H-bonds are shown.
        If --csv_hbonds is given, hydrogen bond residues per pose are saved to a file.
        If --distance or --angle are given, they override the default H-bond cutoff values.
        """
        if not self.session.poses:
            print("Poses data not available. Ensure all required inputs are loaded.")
            return

        tokens = shlex.split(line or "")
        pose_index = None
        csv = None
        csv_hbonds = None
        show_res_only = "--res" in tokens

        i = 0
        while i < len(tokens):
            if tokens[i] == "--csv" and i + 1 < len(tokens):
                csv = tokens[i + 1]
                i += 2
            elif tokens[i] == "--csv_hbonds" and i + 1 < len(tokens):
                csv_hbonds = tokens[i + 1]
                i += 2
            elif tokens[i] == "--res":
                i += 1
            else:
                try:
                    pose_index = int(tokens[i]) - 1
                except ValueError:
                    print(f"Ignoring unrecognized argument: {tokens[i]}")
                i += 1
        handle_poses_info(self.session.poses, pose_index, csv, show_res_only, csv_hbonds)
    
    def do_exit(self, line=None) -> None:
        """
        Exits the interactive shell.

        Usage:
            exit
        """
        print("Exiting DockInspect shell.")
        return True

    def do_EOF(self, line=None) -> None:
        """
        Handles Ctrl+D (EOF) to exit the shell.
        """
        print() 
        return self.do_exit(line)
    
@click.command()
@click.argument("ligand_smiles", type=click.STRING)
@click.argument("pdb_code", type=click.STRING)
@click.argument("vina_file", type=click.Path(exists=True))
@click.argument("structure_file", type=click.Path(exists=True))
@click.argument("predictions_file", required=False, type=click.Path(exists=True))
@click.argument("residues_file", required=False, type=click.Path(exists=True))
@click.option("--distance", required=False, type=float, default=3.2, help="Distance cutoff for hydrogen bonds (A).")
@click.option("--angle", required=False, type=float, default=25.0, help="Angle cutoff for hydrogen bonds (degrees).")
def launch_shell(ligand_smiles: str, pdb_code: str, vina_file: str, structure_file: str, predictions_file: str = None, 
                 residues_file: str = None, distance: float = 3.2, angle: float = 25.0) -> None:
    """
    Initializes a docking analysis session and launches the interactive shell.

    :param ligand_smiles: SMILES string of the ligand.
    :param pdb_code: Protein structure PDB code.
    :param vina_file: Output file from AutoDock Vina.
    :param predictions_file: Pocket prediction file (optional).
    :param residues_file: Residue annotation file (optional).
    :param distance: Distance cutoff for hydrogen bond detection in Angstroms (default: 3.2).
    :param angle: Angle cutoff for hydrogen bond detection in degrees (default: 25.0).
    :return: None
    """
    model_hbonds, model_hbonds_res = get_hydrogen_bonds_and_residues(pdb_code, vina_file, distance, angle)
    session = Session(smiles=ligand_smiles, pdb_code=pdb_code, vina_file=vina_file, structure_file=structure_file, predictions_file=predictions_file, 
                      residues_file=residues_file, distance=distance, angle=angle, model_hbonds=model_hbonds)

    if predictions_file and residues_file:
        try:
            session.pockets = Pockets(structure_file, predictions_file, residues_file)
        except Exception as e:
            print(f"Warning: Failed to initialize pocket data: {e}")
    else:
        print("Pocket prediction files not provided. 'pocket_info' and 'poses_info' will be disabled.")

    if session.ligand and session.pockets and vina_file and pdb_code and predictions_file:
        try:
            session.poses = Poses(ligand=session.ligand, pockets=session.pockets, pdb_code=session.pdb_code, vina_file=vina_file, 
                                  predictions_file=predictions_file, distance=distance, angle=angle, model_hbonds_res=model_hbonds_res)
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
@click.option("--csv", type=click.Path(), help="Save output to file.")
def ligand_info(ligand_smiles: str, csv: str = None) -> None:
    """
    Shows properties of a ligand from a SMILES string.

    :param ligand_smiles: SMILES string representing the ligand.
    :param csv: Optional path to save the output as a CSV file.
    :return: None
    """
    ligand = Ligand(ligand_smiles)
    handle_ligand_info(ligand, csv)

@cli.command(name="pocket_info")
@click.argument("structure_file", type=click.Path(exists=True))
@click.argument("predictions_file", type=click.Path(exists=True))
@click.argument("residues_file", type=click.Path(exists=True))
@click.option("--pocket_id", help="Show specific pocket.")
@click.option("--csv", type=click.Path(), help="Save all pockets to file.")
def pocket_info(structure_file: str, predictions_file: str, residues_file: str,
                pocket_id: str = None, csv: str = None) -> None:
    """
    Shows SASA, GRAVY and charge for predicted pockets.

    :param structure_file: Path to the PDBQT structure file.
    :param predictions_file: Path to the CSV pocket predictions file.
    :param residues_file: Path to the CSV pocket-to-residues mapping file.
    :param pocket_id: Optional ID of a specific pocket to display.
    :param csv: Optional path to save pocket info as a CSV file.
    :return: None
    """
    pockets = Pockets(structure_file, predictions_file, residues_file)
    handle_pocket_info(pockets, pocket_id, csv)

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
@click.option("--distance", type=float, default=3.2, help="Distance cutoff for hydrogen bonds (A).")
@click.option("--angle", type=float, default=25.0, help="Angle cutoff for hydrogen bonds (degrees).")
def poses_info(ligand_smiles: str, pdb_code: str, vina_file: str, structure_file: str, predictions_file: str, 
               residues_file: str, pose_index: int = None, csv: str = None, res: bool = False, 
               csv_hbonds: str = None, distance: float = 3.2, angle: float = 25.0) -> None:
    """
    Shows ligand–pocket interaction summary per docking pose.

    :param ligand_smiles: SMILES string of the ligand.
    :param pdb_code: PDB code of the protein.
    :param vina_file: Path to the PDBQT AutoDock Vina output file.
    :param structure_file: Path to the PDBQT protein structure file.
    :param predictions_file: Path to the CSV pocket predictions file.
    :param residues_file: Path to the CSV pocket-to-residues mapping file.
    :param pose_index: Optional index of a specific pose to display (1-based).
    :param csv: Optional path to save poses info as a CSV file.
    :param res: If provided only residue names from H-bonds are shown.
    :param csv_hbonds: Optional path to save hydrogen bond residues per pose are saved to a CSV file.
    :param distance: Distance cutoff for hydrogen bond detection in Angstroms (default: 3.2).
    :param angle: Angle cutoff for hydrogen bond detection in degrees (default: 25.0).
    :return: None
    """
    ligand = Ligand(ligand_smiles)
    pockets = Pockets(structure_file, predictions_file, residues_file)
    poses = Poses(ligand, pockets, pdb_code, vina_file, predictions_file, distance, angle)
    handle_poses_info(poses, pose_index - 1 if pose_index else None, csv, res, csv_hbonds)

@cli.command(name="visualize")
@click.argument("pdb_code", type=click.STRING)
@click.argument("vina_file", type=click.Path(exists=True))
@click.option("--pose", type=int, default=1, help="Pose number to visualize (1-based).")
@click.option("--mode", type=click.Choice(["surface", "polar", "charge", "hbonds"]), default="hbonds")
@click.option("--predictions_file", type=click.Path(exists=True), required=False)
@click.option("--residues_file", type=click.Path(exists=True), required=False)
@click.option("--distance", type=float, default=3.2, help="Distance cutoff for hydrogen bonds (A).")
@click.option("--angle", type=float, default=25.0, help="Angle cutoff for hydrogen bonds (degrees).")
@click.option("--save", type=click.Path(), help="Save the PyMOL script to the given file instead of launching PyMOL.")
def visualize_command(pdb_code: str, vina_file: str, pose: int = 1, mode: str = "hbonds", predictions_file: str = None, 
                      residues_file: str = None, distance: float = 3.2, angle: float = 25.0, save: str = None) -> None:
    """
    Launches PyMOL visualization of a selected docking pose.

    :param pdb_code: PDB code of the protein.
    :param vina_file: Path to the PDBQT AutoDock Vina output file.
    :param pose: Pose number to visualize (1-based).
    :param mode: Visualization mode from (surface, polar, charge, hbonds) (default hbonds).
    :param predictions_file: Optional pocket predictions CSV file.
    :param residues_file: Optional pocket residue map CSV file.
    :param distance: Distance cutoff for hydrogen bond detection in Angstroms (default: 3.2).
    :param angle: Angle cutoff for hydrogen bond detection in degrees (default: 25.0).
    :return: None
    """
    poses = None
    if predictions_file and residues_file:
        try:
            ligand = Ligand("C")
            pockets = Pockets(vina_file, predictions_file, residues_file)
            poses = Poses(ligand, pockets, pdb_code, vina_file, predictions_file, distance, angle)
        except Exception as e:
            print(f"Warning: Could not load poses for pocket mapping: {e}")

    model_hbonds = get_hydrogen_bonds_and_residues(pdb_code, vina_file, distance, angle)[0]

    handle_visualization(
        pdb_code=pdb_code,
        vina_file=vina_file,
        pose_num=pose,
        mode=mode,
        predictions_file=predictions_file,
        poses=poses,
        distance=distance,
        angle=angle,
        save=save,
        model_hbonds=model_hbonds
    )

def main():
    cli() 