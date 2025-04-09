"""
Calculates physicochemical properties of a ligand from its SMILES string using RDKit to compute logP (hydrophobicity), SASA, TPSA, molecular volume, and net formal charge.

Classes:
    Ligand: Represents a small molecule with computed properties.
"""

from rdkit import Chem
from rdkit.Chem import Crippen, rdFreeSASA, AllChem
from rdkit.Chem.rdMolDescriptors import DoubleCubicLatticeVolume

class Ligand:
    """
    Represents a ligand and computes key molecular properties from a SMILES string.

    :param smiles: A SMILES string representing the ligand structure.

    Attributes:
        smiles (str): The input SMILES string.
        logp (float): Calculated octanol-water partition coefficient (hydrophobicity).
        sasa (float): Solvent accessible surface area (Å²).
        tpsa (float): Topological polar surface area (Å²).
        volume (float): Molecular volume estimated using a lattice-based method (Å³).
        charge (int): Net formal charge of the molecule.
    """
    def __init__(self, smiles):
        self.smiles = smiles
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())

        self.logp = Crippen.MolLogP(mol)
        self.sasa = rdFreeSASA.CalcSASA(mol, rdFreeSASA.classifyAtoms(mol))
        self.tpsa = Chem.rdMolDescriptors.CalcTPSA(mol)
        self.volume = DoubleCubicLatticeVolume(mol).GetVolume()
        self.charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])

    def __str__(self):
        """
        Returns a formatted string of ligand properties: LogP, SASA, TPSA, Volume, and Net Charge.

        :return: Formatted multiline string summarizing the ligand’s computed properties.
        """
        return (
            f"Ligand Properties for SMILES: {self.smiles}\n"
            f"{'-'*40}\n"
            f"{'logP:':<12} {self.logp:.2f}\n"
            f"{'SASA (Å²):':<12} {self.sasa:.2f}\n"
            f"{'TPSA (Å²):':<12} {self.tpsa:.2f}\n"
            f"{'Volume (Å³):':<12} {self.volume:.2f}\n"
            f"{'Charge:':<12} {self.charge:+d}"
        )



if __name__=="__main__":
    #smiles = "C1CCCCC1" 
    #smiles = "C1CCC(O)CC1" 
    #smiles = "NC(=O)N"
    #smiles="[NH4+]" 
    smiles="CC(=O)[O-]"
    
    ligand = Ligand(smiles)
    print(ligand)