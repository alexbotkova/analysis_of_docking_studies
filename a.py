from pymol import cmd
import pyvol
from pyvol.identify import pocket_wrapper

print(dir(pyvol))

def load_and_visualize_pocket(structure_file, vina_file, cutoff=3.6):
    cmd.load(structure_file, 'protein')
    cmd.load(vina_file, 'ligand')
    
    pocket_selection = f'br. (protein within {cutoff} of ligand)'
    cmd.select('pocket', pocket_selection)
    
    print(f"Pocket residues within {cutoff} Å of the ligand have been selected.")



receptor_file = "/Users/alexbotkova/analysis_of_docking_studies/test_files/urea/result-2025-03-30T21_20_58.783Z/structure.pdbqt"
ligand_file = "/Users/alexbotkova/analysis_of_docking_studies/test_files/urea/result-2025-03-30T21_20_58.783Z/out_vina.pdbqt"

load_and_visualize_pocket(receptor_file, ligand_file)
