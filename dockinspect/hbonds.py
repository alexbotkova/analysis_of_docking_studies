"""
This module provides PyMOL-based function to detect hydrogen bonds between
a protein structure (fetched via PDB code) and ligand poses from AutoDock Vina output.

Functions:
    get_hydrogen_bonds_and_residues: Returns hydrogen bond atom pairs for each pose and residue names involved in hydrogen bonds per pose.
"""

def get_hydrogen_bonds_and_residues(pdb_code: str, vina_file: str, distance: float, angle: float) -> tuple:
    """
    Detects hydrogen bonds between protein and ligand for each pose, returning both
    atom index pairs and involved residue names.

    :param pdb_code: PDB code for the structure to fetch.
    :param vina_file: Path to the PDBQT AutoDock Vina output file.
    :param distance: Distance cutoff for H-bond detection.
    :param angle: Angle cutoff for H-bond detection.
    :return: Tuple of two lists:
             - List of H-bond index pairs per pose.
             - List of involved residue names per pose.
    """
    import pymol
    from pymol import cmd

    pymol.finish_launching(['pymol', '-qc'])

    cmd.fetch(pdb_code, name="structure", type="pdb")
    cmd.load(vina_file, "out_vina")
    cmd.h_add()

    model_hbonds = []
    model_hbonds_res = []

    n_states = cmd.count_states("out_vina")

    for state in range(1, n_states + 1):
        sel1 = 'structure and (donor or acceptor) and (elem N+O)'
        sel2 = 'out_vina and (donor or acceptor) and (elem N+O)'

        hbonds = cmd.find_pairs(sel1, sel2, mode=1,
                                cutoff=distance,
                                angle=angle,
                                state1=1, state2=state)
        model_hbonds.append(hbonds)

        hbonds_res = []
        for (m1, i1), (m2, i2) in hbonds:
            if m1 == "structure":
                res_info1 = cmd.get_model(f"{m1} and index {i1}").atom[0]
                hbonds_res.append(res_info1.resn)
            if m2 == "structure":
                res_info2 = cmd.get_model(f"{m2} and index {i2}").atom[0]
                hbonds_res.append(res_info2.resn)

        model_hbonds_res.append(hbonds_res)

    cmd.quit()
    return model_hbonds, model_hbonds_res