from rdkit import Chem
from rdkit.Chem import Crippen, rdFreeSASA, AllChem
from rdkit.Chem.rdMolDescriptors import DoubleCubicLatticeVolume

class Ligand:
    def __init__(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())

        self.logp = Crippen.MolLogP(mol)
        print(f"LogP: {self.logp}")

        self.sasa = rdFreeSASA.CalcSASA(mol, rdFreeSASA.classifyAtoms(mol))
        print(f"SASA: {self.sasa} Å^2")

        self.tpsa = Chem.rdMolDescriptors.CalcTPSA(mol)
        print(f"TPSA: {self.tpsa} Å^2")

        self.volume = DoubleCubicLatticeVolume(mol).GetVolume()
        print(f"Volume: {self.volume:.2f} Å^3")

        self.charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
        print(f"Charge: {self.charge}")

    @staticmethod
    def get_ratio(pocket_val, ligand_val):
        return ligand_val / pocket_val if pocket_val != 0 else 0

    def print_ratios(self, pocket_dict, ligand_val):
        for pocket_name in pocket_dict.keys():
            pocket_val = pocket_dict[pocket_name]
            print(f"{pocket_name}: {self.get_ratio(pocket_val, ligand_val) * 100:.0f} %")
    
    def print_charge_complementarity(self, pocket_charges):
        for pocket_name in pocket_charges.keys():
            pocket_charge = pocket_charges[pocket_name]
            print(f"{pocket_name}: {pocket_charge * self.charge}")



if __name__=="__main__":
    #smiles = "C1CCCCC1" 
    smiles = "C1CCC(O)CC1" 
    ligand = Ligand(smiles)

    pocket_volumes = {'pocket1': 3042.8, "pocket2": 675.8, "pocket3": 235.2, "pocket4": 356.0, "pocket5": 155.0, "pocket6": 213.7}
    pocket_sasas = {'pocket1': 1286.9763832695578, 'pocket2': 689.0772473978619, 'pocket3': 549.3877405340272, 'pocket4': 626.0000529927811, 'pocket5': 291.5486064351807, 'pocket6': 210.1965429121765}
    pocket_charges = {'pocket1': 0, 'pocket2': 0, 'pocket3': 2, 'pocket4': -2, 'pocket5': 0, 'pocket6': 0}

    ligand.print_ratios(pocket_volumes, ligand.volume)
    print()
    ligand.print_ratios(pocket_sasas, ligand.sasa)
    print()
    ligand.print_charge_complementarity(pocket_charges)