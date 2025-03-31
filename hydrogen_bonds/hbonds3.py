from pymol import cmd
from get_raw_distances import get_raw_distances

def main():
    # Adjust paths here
    receptor_file = "test_files/urea/result-2025-03-30T21_20_58.783Z/structure.pdbqt"
    ligand_file = "test_files/urea/result-2025-03-30T21_20_58.783Z/out_vina.pdbqt"

    # Load receptor and ligand
    cmd.load(receptor_file, "receptor")
    cmd.load(ligand_file, "ligand")

    # Add hydrogens
    cmd.h_add("receptor")
    cmd.h_add("ligand")

    # Measure hydrogen-bond-like distances using PyMOL's built-in heuristic
    cmd.distance("hbonds", "receptor", "ligand", mode=2, cutoff=3.5)

    # Use get_raw_distances to extract info
    distances = get_raw_distances("hbonds")
    print(f"\nHydrogen bonds found (via distance object): {len(distances)}")

    # Optional: show actual distances or involved atoms
    for d in distances:
        dist = d[2]
        print(f" - Distance: {dist:.2f} Å")

if __name__ == "__main__":
    main()
