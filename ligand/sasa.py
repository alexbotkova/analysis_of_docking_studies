

import freesasa
print(1)

"""
def convert_pdbqt_to_pdb(pdbqt_file, pdb_file):
    with open(pdbqt_file, 'r') as f_in, open(pdb_file, 'w') as f_out:
        for line in f_in:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                f_out.write(line[:66] + "\n")

#convert_pdbqt_to_pdb("/Users/alexbotkova/analysis_of_docking_studies/test_files/result-2024-12-10T22_19_39.203Z/ligand.pdbqt", "ligand.pdb")
"""

structure = freesasa.Structure("/Users/alexbotkova/analysis_of_docking_studies/ligand/CHX.pdb")
result = freesasa.calc(structure)
print("Total SASA:", result.totalArea())

"""

from Bio.PDB import PDBParser
import freesasa

# Load your structure using Bio.PDB; ensure the file is in PDB format.
parser = PDBParser(QUIET=True)
structure_bio = parser.get_structure("ligand", "/Users/alexbotkova/analysis_of_docking_studies/test_files/result-2024-12-10T22_19_39.203Z/ligand2.pdb")
print(structure_bio)
# Calculate SASA using calcBioPDB.
result, sasa_classes, residue_areas = freesasa.calcBioPDB(structure_bio)

# Output the total SASA.
print("Total SASA:", result.totalArea())
"""