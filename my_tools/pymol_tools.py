import re

def get_pocket_residues_dict(pocket_data_df):
    """
    Generates a dictionary where keys are pocket names and values are residues joined in one string 
    formatted for the pymol selection command.

    :param pocket_data_df: A dataframe generated from a CSV file containing PDB predictions.
    :return: A dictionary where keys are pocket names and values are residues joined in one string 
    formatted for the pymol selection command.
    """
    pocket_residues_dict = {}
    for _, row in pocket_data_df.iterrows():
            pocket_name = row['name'].strip()
            residue_ids = row['residue_ids']
            matches = re.findall(r'([A-Z])_(\d+)', residue_ids)
            if matches:
                chain_letter = matches[0][0] 
                residues = [residue for _, residue in matches]
                joined_residues = "+".join(residues)
            residues_selection = f"resi {joined_residues} and chain {chain_letter}"
            pocket_residues_dict[pocket_name] = residues_selection
    return pocket_residues_dict