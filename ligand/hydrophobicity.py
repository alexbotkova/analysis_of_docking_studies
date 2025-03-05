from rdkit import Chem
from rdkit.Chem import Crippen, rdFreeSASA, AllChem

smiles = "C1CCCCC1"  

mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)

AllChem.EmbedMolecule(mol, AllChem.ETKDG())

logp = Crippen.MolLogP(mol)
print(f"LogP for {smiles}: {logp}")

radii = rdFreeSASA.classifyAtoms(mol)

sasa = rdFreeSASA.CalcSASA(mol, radii)
print(f"SASA: {sasa} Å^2")