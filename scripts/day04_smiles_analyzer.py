from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

print("=" * 50)
print("RDKit Drug Property Analyzer")
print("=" * 50)

# Get SMILES from user
smiles = input("Enter a SMILES string: ")

# Convert to RDKit molecule
mol = Chem.MolFromSmiles(smiles)

if mol is None:
    print("\nInvalid SMILES!")
else:
    formula = rdMolDescriptors.CalcMolFormula(mol)
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    rot = Descriptors.NumRotatableBonds(mol)

    print("\n" + "=" * 50)
    print("Drug Property Report")
    print("=" * 50)

    print(f"Formula          : {formula}")
    print(f"Molecular Weight : {mw:.2f}")
    print(f"LogP             : {logp:.2f}")
    print(f"TPSA             : {tpsa:.2f}")
    print(f"HBD              : {hbd}")
    print(f"HBA              : {hba}")
    print(f"Rotatable Bonds  : {rot}")

    if mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10:
        print("\nLipinski Result : PASS")
    else:
        print("\nLipinski Result : FAIL")