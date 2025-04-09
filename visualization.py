"""
Launches PyMOL with a visualization script for binding poses from AutoDock Vina.

This script loads a protein structure and ligand pose from AutoDock Vina output,
visualizes selected modes (e.g. surfaces, polar regions, charge distribution, hydrogen bonds),
and optionally highlights predicted binding pockets from tools like P2Rank.

Functions:
    visualize: Main visualization interface for PyMOL.
    __color_regions: Helper function to color specific residue types.

Requires PyMOL, structure and Vina output files, and optionally a pocket selection string.
"""

from pymol import cmd 
import sys
sys.path.append("..")
from my_tools.my_parser import *
from my_tools.pymol_tools import *

def __color_regions(selection_name, resn, color):
    """
    Creates a PyMOL selection and colors it based on residue name.

    :param selection_name: Name for the PyMOL selection.
    :param resn: Residue names to select.
    :param color: Color to apply.
    :return: None
    """
    cmd.select(selection_name, f'(pocket) and resn {resn}')
    cmd.color(color, selection_name)

def visualize(structure_file, vina_file, pose_num=1, pocket_selection = "",  mode="", hbonds=False, cutoff=3.6, not_pocket_surface_transparency=0, pocket_surface_transparency=0, color_mode="", show_pocket_surface=False, show_ligand_surface=False, show_not_pocket_surface=False, show_pocket_sticks=False, show_ligand_sticks=False):
    """
    Launches PyMOL with a visualization script for binding poses from AutoDock Vina.

    :param structure_file: Path to the protein structure file.
    :param vina_file: Path to the out_vina file from AutoDock Vina.
    :param pose_num: Pose number to visualize (1-based index, default is 1).
    :param pocket_selection: PyMOL atom selection string for highlighting a pocket (from p2rank predictions).
    :param mode: Visualization mode ('surface', 'polar', 'charge', 'hbonds', default is 'hbonds').
    :param hbonds: Whether to compute and show hydrogen bonds between ligand and protein (default False).
    :param cutoff: Distance cutoff (in Å) for defining binding pocket if no pocket selection is given.
    :param not_pocket_surface_transparency: Transparency level for non-pocket regions (default 0).
    :param pocket_surface_transparency: Transparency level for pocket surface (default 0).
    :param color_mode: Coloring scheme for residues ('broad', 'detailed', or empty for none).
    :param show_pocket_surface: Whether to show surface representation of pocket residues.
    :param show_ligand_surface: Whether to show surface representation of ligand.
    :param show_not_pocket_surface: Whether to show surface for non-pocket protein regions.
    :param show_pocket_sticks: Whether to show pocket residues as sticks.
    :param show_ligand_sticks: Whether to show ligand as sticks.
    :return: None
    """
    
    #TODO chytreji...
    if mode == "": mode = "hbonds"
    
    cmd.load(structure_file, 'structure')
    cmd.load(vina_file, 'out_vina')
    cmd.frame(pose_num)
    cmd.h_add()
    cmd.hide('everything', 'all')

    
    if not pocket_selection: 
        pocket_selection = f'br. (structure within {cutoff} of out_vina)'
        cmd.select('pocket', pocket_selection)
    else:
        pocket_selection = 'structure and ' + pocket_selection
        cmd.select('pocket2', pocket_selection)
        cmd.select('pocket', 'byres pocket2')
    print(pocket_selection)
    cmd.select('not_pocket', 'not (pocket or out_vina)')

    if mode == "surface":
        show_pocket_surface = True
        show_not_pocket_surface = True 
        show_ligand_surface = True 
        not_pocket_surface_transparency = 0.8
        pocket_surface_transparency=0
        cmd.color('hotpink', 'out_vina')
        cmd.color('white', 'pocket')
        cmd.color('grey', 'not_pocket')
    elif mode == 'polar':
        show_pocket_surface = True
        show_not_pocket_surface = True 
        show_ligand_sticks = True
        color_mode = 'broad'
        not_pocket_surface_transparency = 0.8
        pocket_surface_transparency=0
        cmd.color('grey', 'not_pocket')
    elif mode == 'charge':
        show_pocket_surface = True
        show_not_pocket_surface = True 
        show_ligand_sticks = True
        color_mode = 'detailed'
        not_pocket_surface_transparency = 0.8
        pocket_surface_transparency=0
        cmd.color('grey', 'not_pocket')
    elif mode == 'hbonds':
        hbonds=True
        show_pocket_sticks = True
        show_ligand_sticks = True

    cmd.set('transparency', pocket_surface_transparency, 'pocket')
    cmd.set('transparency', not_pocket_surface_transparency, 'not_pocket')
    
    if show_pocket_surface:
        cmd.show('surface', 'pocket')
    if show_not_pocket_surface:
        cmd.show('surface', 'not_pocket')
    if show_ligand_surface:
        cmd.show('surface', 'out_vina')
    
    if show_pocket_sticks:
        cmd.show('sticks', 'pocket')
    if show_ligand_sticks:
        cmd.show('sticks', 'out_vina')

    if hbonds:
        #TODO diskuze přidání res do pocketu
        sel1 = 'structure and (donor or acceptor)'
        sel2 = 'out_vina and (donor or acceptor)'
        hbonds = cmd.find_pairs(sel1, sel2, mode=1, cutoff=3.5, angle=55.0)
        for idx, ((m1, i1), (m2, i2)) in enumerate(hbonds):
            cmd.distance(f"hb_{idx}", f"{m1} and index {i1}", f"{m2} and index {i2}")

            cmd.select(f"hb_res_{idx}_1", f"byres ({m1} and index {i1})")
            cmd.select(f"hb_res_{idx}_2", f"byres ({m2} and index {i2})")
            cmd.show("sticks", f"hb_res_{idx}_1 or hb_res_{idx}_2")

    if color_mode=="broad":
        __color_regions('hydrophobic', 'ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO+GLY', 'yellow')
        __color_regions('hydrophilic', 'SER+THR+ASN+GLN+TYR+CYS+HIS+ARG+LYS+ASP+GLU', 'blue')
    elif color_mode=="detailed":
        __color_regions('hydrophobic', 'ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO+GLY', 'yellow')
        __color_regions('acidic', 'ASP+GLU', 'red')
        __color_regions('basic', 'LYS+ARG+HIS', 'blue')
        __color_regions('neutral', 'SER+THR+ASN+GLN+TYR+CYS', 'white')
    else:
        pass

if __name__=="__main__":
    predictions_filepath = "test_files/cyclohexane/prankweb-2SRC/structure.pdb_predictions.csv"
    structure_filepath = "/Users/alexbotkova/analysis_of_docking_studies/test_files/urea/result-2025-03-30T21_20_58.783Z/structure.pdbqt"
    out_vina_filepath = "/Users/alexbotkova/analysis_of_docking_studies/test_files/urea/result-2025-03-30T21_20_58.783Z/out_vina.pdbqt"

    #visualize(structure_filepath, out_vina_filepath, mode="surface")
    #visualize(structure_filepath, out_vina_filepath, mode="polar")
    #visualize(structure_filepath, out_vina_filepath, mode="charge")
    #visualize(structure_filepath, out_vina_filepath, mode="hbonds")
    pocket_data_df = get_df(predictions_filepath)
    pocket_residues_dict = get_pocket_residues_dict(pocket_data_df)
    print(pocket_residues_dict)