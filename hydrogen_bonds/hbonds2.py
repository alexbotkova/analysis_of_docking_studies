from pymol import cmd
from collections import Counter

# cyclohexane
#receptor_file = "test_files/cyclohexane/result-2024-12-10T22_19_39.203Z/structure.pdbqt"
#ligand_file = "test_files/cyclohexane/result-2024-12-10T22_19_39.203Z/out_vina.pdbqt"

# urea 
receptor_file = "/Users/alexbotkova/analysis_of_docking_studies/test_files/urea/result-2025-03-30T21_20_58.783Z/structure.pdbqt"
ligand_file = "/Users/alexbotkova/analysis_of_docking_studies/test_files/urea/result-2025-03-30T21_20_58.783Z/out_vina.pdbqt"

cmd.load(receptor_file, 'receptor')
cmd.load(ligand_file, 'ligand')

cmd.h_add('receptor')
cmd.h_add('ligand')

sel1 = 'receptor and (donor or acceptor)'
sel2 = 'ligand and (donor or acceptor)'

pairs = pairs = cmd.find_pairs(sel1, sel2, mode=1, cutoff=3.5, angle=55.0)
print(f"Total number of hydrogen bond pairs: {len(pairs)}")

for idx, ((m1, i1), (m2, i2)) in enumerate(pairs):
    cmd.distance(f"hb_{idx}", f"{m1} and index {i1}", f"{m2} and index {i2}")


residue_counts = Counter()
for (m1, i1), (m2, i2) in pairs:
    try:
        a1 = cmd.get_model(f"{m1} and id {i1}").atom[0]
        a2 = cmd.get_model(f"{m2} and id {i2}").atom[0]
        res1 = f"{a1.segi}/{a1.chain}/{a1.resn}{a1.resi}"
        res2 = f"{a2.segi}/{a2.chain}/{a2.resn}{a2.resi}"
        residue_counts[res1] += 1
        residue_counts[res2] += 1
    except IndexError:
        continue

print("\nHydrogen bond interactions per residue:")
for res, count in residue_counts.items():
    print(f"{res}: {count}")