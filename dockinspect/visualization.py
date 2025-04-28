"""
PyMOL visualization utilities for AutoDock Vina results

This module provides tools to visualize protein–ligand interactions in PyMOL
based on the output from AutoDock Vina. It supports direct visualization
via PyMOL and also generates reusable PyMOL scripts.

Functions:
    visualize: Launches PyMOL with the desired visualization.
    save_pml_script: Writes a PyMOL script instead of launching the program.
"""

import os
    
def visualize(pdb_code: str, vina_file: str, pose_num: int = 1, pocket_selection: str = "", mode: str = "", distance: float = 3.2, angle: float = 25) -> None:
    """
    Launches PyMOL with a visualization script for binding poses from AutoDock Vina.

    :param pdb_code: PDB code for the protein structure.
    :param vina_file: Path to the out_vina file from AutoDock Vina.
    :param pose_num: Pose number to visualize (1-based index, default is 1).
    :param pocket_selection: PyMOL selection string for highlighting a pocket.
    :param mode: Visualization mode (surface, polar, charge, hbonds) (default hbonds).
    :param distance: Distance cutoff for hydrogen bond detection in A (default: 3.2).
    :param angle: Angle cutoff for hydrogen bond detection in degrees (default: 25).
    
    Internal optional settings (hardcoded inside function):
    - hbonds: Whether to compute and display hydrogen bonds.
    - cutoff: Distance threshold for defining pocket if no selection is given (A).
    - not_pocket_surface_transparency: Transparency for non-pocket surface.
    - pocket_surface_transparency: Transparency for pocket surface.
    - color_mode: Residue coloring scheme (broad, detailed, or "").
    - show_pocket_surface: Show pocket surface.
    - show_ligand_surface: Show ligand surface.
    - show_not_pocket_surface: Show surface of non-pocket protein regions.
    - show_pocket_sticks: Show pocket residues as sticks.
    - show_ligand_sticks: Show ligand as sticks.

    :return: None
    """
    hbonds = False
    cutoff = 3.6
    not_pocket_surface_transparency = 0
    pocket_surface_transparency = 0
    color_mode = ""
    show_pocket_surface = False
    show_ligand_surface = False 
    show_not_pocket_surface = False
    show_pocket_sticks = False
    show_ligand_sticks = False

    from pymol import cmd 
    if mode == "":
        mode = "hbonds"

    cmd.fetch(pdb_code, name="structure", type="pdb")

    cmd.load(vina_file, "out_vina")
    cmd.frame(pose_num)
    cmd.h_add()
    cmd.hide("everything", "all")

    if not pocket_selection:
        pocket_selection = f"br. (structure within {cutoff} of out_vina)"
        cmd.select("pocket", pocket_selection)
    else:
        pocket_selection = "structure and " + pocket_selection
        cmd.select("pocket", pocket_selection)

    cmd.select("not_pocket", "not (pocket or out_vina)")

    if mode == "surface":
        show_pocket_surface = True
        show_not_pocket_surface = True
        show_ligand_surface = True
        not_pocket_surface_transparency = 0.8
        pocket_surface_transparency = 0
        cmd.color("hotpink", "out_vina")
        cmd.color("white", "pocket")
        cmd.color("grey", "not_pocket")
    elif mode == "polar":
        show_pocket_surface = True
        show_not_pocket_surface = True
        show_ligand_sticks = True
        color_mode = "broad"
        not_pocket_surface_transparency = 0.8
        pocket_surface_transparency = 0
        cmd.color("grey", "not_pocket")
    elif mode == "charge":
        show_pocket_surface = True
        show_not_pocket_surface = True
        show_ligand_sticks = True
        color_mode = "detailed"
        not_pocket_surface_transparency = 0.8
        pocket_surface_transparency = 0
        cmd.color("grey", "not_pocket")
    elif mode == "hbonds":
        hbonds = True
        show_pocket_sticks = True
        show_ligand_sticks = True

    cmd.set("transparency", pocket_surface_transparency, "pocket")
    cmd.set("transparency", not_pocket_surface_transparency, "not_pocket")

    if show_pocket_surface:
        cmd.show("surface", "pocket")
    if show_not_pocket_surface:
        cmd.show("surface", "not_pocket")
    if show_ligand_surface:
        cmd.show("surface", "out_vina")
    if show_pocket_sticks:
        cmd.show("sticks", "pocket")
    if show_ligand_sticks:
        cmd.show("sticks", "out_vina")

    if hbonds:
        sel1 = 'structure and (donor or acceptor) and (elem N+O)'
        sel2 = 'out_vina and (donor or acceptor) and (elem N+O)'
        found_hbonds = cmd.find_pairs(sel1, sel2, mode=1, cutoff=distance, angle=angle, state1=1, state2=pose_num)
        for idx, ((m1, i1), (m2, i2)) in enumerate(found_hbonds):
            cmd.distance(f"hb_{idx}", f"{m1} and index {i1}", f"{m2} and index {i2}")
            cmd.select(f"hb_res_{idx}_1", f"byres ({m1} and index {i1})")
            cmd.select(f"hb_res_{idx}_2", f"byres ({m2} and index {i2})")
            cmd.show("sticks", f"hb_res_{idx}_1 or hb_res_{idx}_2")

    def color_regions(selection_name, resn, color):
        cmd.select(selection_name, f'(pocket) and resn {resn}')
        cmd.color(color, selection_name)

    if color_mode == "broad":
        color_regions("hydrophobic", "ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO+GLY", "yellow")
        color_regions("hydrophilic", "SER+THR+ASN+GLN+TYR+CYS+HIS+ARG+LYS+ASP+GLU", "blue")
    elif color_mode == "detailed":
        color_regions("hydrophobic", "ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO+GLY", "yellow")
        color_regions("acidic", "ASP+GLU", "red")
        color_regions("basic", "LYS+ARG+HIS", "blue")
        color_regions("neutral", "SER+THR+ASN+GLN+TYR+CYS", "white")
    
    cmd.select("sele", "none")

    fetched_filename = f"{pdb_code}.pdb"
    if os.path.exists(fetched_filename):
        try:
            os.remove(fetched_filename)
        except Exception as e:
            print(f"Warning: could not remove {fetched_filename}: {e}")

def save_pml(pdb_code: str, vina_file: str, pose_num: int = 1, pocket_selection: str = "", mode: str = "", model_hbonds = None, pml: str = "visualization.pml") -> None:
    """
    Writes a PyMOL visualization script (.pml) to file.

    :param pdb_code: PDB code (to fetch) or path to a local .pdb/.pdbqt file.
    :param vina_file: Path to the AutoDock Vina output file (e.g. out_vina).
    :param pose_num: Pose number to visualize (1-based index, default: 1).
    :param pocket_selection: PyMOL atom selection string to highlight a pocket (optional).
    :param mode: Visualization mode (surface, polar, charge, hbonds) (default hbonds).
    :param model_hbonds: Precomputed hydrogen bond data.
    :param pml: Path to save the generated .pml script (default: "visualization.pml").

    Internal optional settings (hardcoded inside function):
    - cutoff: Distance threshold for defining pocket if no selection is given (A).
    - not_pocket_surface_transparency: Transparency for non-pocket surface.
    - pocket_surface_transparency: Transparency for pocket surface.
    - color_mode: Residue coloring scheme (broad, detailed, or "").
    - show_pocket_surface: Show pocket surface.
    - show_ligand_surface: Show ligand surface.
    - show_not_pocket_surface: Show surface of non-pocket protein regions.
    - show_pocket_sticks: Show pocket residues as sticks.
    - show_ligand_sticks: Show ligand as sticks.

    :return: None
    """
    cutoff = 3.6
    not_pocket_surface_transparency = 0
    pocket_surface_transparency = 0
    color_mode = ""
    show_pocket_surface = False
    show_ligand_surface = False 
    show_not_pocket_surface = False
    show_pocket_sticks = False
    show_ligand_sticks = False

    if mode == "":
        mode = "hbonds"

    with open(pml, "w") as f:
        write = f.write
        write("reinitialize\n")
        
        write(f"fetch {pdb_code},name=structure, type=pdb, async=0\n")
        write(f"load {vina_file}, out_vina\n")
        write(f"frame {pose_num}\n")
        write("h_add\n")
        write("hide everything, all\n")

        if not pocket_selection:
            write(f"select pocket, br. (structure within {cutoff} of out_vina)\n")
        else:
            write(f"select pocket, structure and ({pocket_selection})\n")
        write("select not_pocket, not (pocket or out_vina)\n")

        if mode == "surface":
            show_pocket_surface = True
            show_not_pocket_surface = True
            show_ligand_surface = True
            not_pocket_surface_transparency = 0.8
            pocket_surface_transparency = 0
            write("color hotpink, out_vina\n")
            write("color white, pocket\n")
            write("color grey, not_pocket\n")
        elif mode == "polar":
            show_pocket_surface = True
            show_not_pocket_surface = True
            show_ligand_sticks = True
            color_mode = "broad"
            not_pocket_surface_transparency = 0.8
            pocket_surface_transparency = 0
            write("color grey, not_pocket\n")
        elif mode == "charge":
            show_pocket_surface = True
            show_not_pocket_surface = True
            show_ligand_sticks = True
            color_mode = "detailed"
            not_pocket_surface_transparency = 0.8
            pocket_surface_transparency = 0
            write("color grey, not_pocket\n")
        elif mode == "hbonds":
            show_pocket_sticks = True
            show_ligand_sticks = True

        write(f"set transparency, {pocket_surface_transparency}, pocket\n")
        write(f"set transparency, {not_pocket_surface_transparency}, not_pocket\n")

        if show_pocket_surface:
            write("show surface, pocket\n")
        if show_not_pocket_surface:
            write("show surface, not_pocket\n")
        if show_ligand_surface:
            write("show surface, out_vina\n")
        if show_pocket_sticks:
            write("show sticks, pocket\n")
        if show_ligand_sticks:
            write("show sticks, out_vina\n")

        if mode == "hbonds":
            for idx, ((m1, i1), (m2, i2)) in enumerate(model_hbonds[pose_num-1]):
                write(f"distance hb_{idx}, {m1} and index {i1}, {m2} and index {i2}\n")
                write(f"select hb_res_{idx}_1, byres ({m1} and index {i1})\n")
                write(f"select hb_res_{idx}_2, byres ({m2} and index {i2})\n")
                write(f"show sticks, hb_res_{idx}_1 or hb_res_{idx}_2\n")

        def __color_regions_save_pml(selection_name, resn, color):
            write(f"select {selection_name}, pocket and resn {resn}\n")
            write(f"color {color}, {selection_name}\n")

        if color_mode == "broad":
            __color_regions_save_pml("hydrophobic", "ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO+GLY", "yellow")
            __color_regions_save_pml("hydrophilic", "SER+THR+ASN+GLN+TYR+CYS+HIS+ARG+LYS+ASP+GLU", "blue")
        elif color_mode == "detailed":
            __color_regions_save_pml("hydrophobic", "ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO+GLY", "yellow")
            __color_regions_save_pml("acidic", "ASP+GLU", "red")
            __color_regions_save_pml("basic", "LYS+ARG+HIS", "blue")
            __color_regions_save_pml("neutral", "SER+THR+ASN+GLN+TYR+CYS", "white")

        write("select none\n")
        write("zoom all\n")

    print(f"PML script saved to {pml}.")