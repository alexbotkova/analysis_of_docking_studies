import pandas as pd
import pymolPy3

def get_pocket_residues(file_path):
    data = pd.read_csv(file_path)

    data.columns = data.columns.str.strip()

    pocket_residues = {}

    for _, row in data.iterrows():
        pocket_name = row['name']
        residues = row['residue_ids'].split()
        formatted_residues = '+'.join([res.split('_')[1] for res in residues])  
        pocket_residues[pocket_name] = formatted_residues
    
    return pocket_residues

def count_hydrogen_bonds_per_protein_pocket(protein_name, pocket_residues):
    pm = pymolPy3.pymolPy3(0)
    pm(f"fetch {protein_name}")
    pm("h_add")

    pocket_h_bonds = {}

    for pocket_name in pocket_residues:
        print(f"Selecting {pocket_name}")
        pm(f"select {pocket_name}, resi {pocket_residues[pocket_name]}")
        
        hbonds = pm(f"cmd.find_pairs('{pocket_name}', '{pocket_name}', cutoff=3.2, angle=55)")
        print(f"Hydrogen bonds for {pocket_name}: {hbonds}")

if __name__ == "__main__":
    #file_path = '/Users/alexbotkova/analysis_of_docking_studies/electrostatics/structure.pdb_predictions.csv'
    count_hydrogen_bonds_per_protein_pocket("2src", get_pocket_residues(file_path))







