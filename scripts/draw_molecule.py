from rdkit import Chem
from rdkit.Chem import Draw


# SMILES of aspirin
smiles = "CC(=O)Oc1ccccc1C(=O)O"


# Convert SMILES to molecule
mol = Chem.MolFromSmiles(smiles)


# Draw molecule and save image
image = Draw.MolToImage(
    mol,
    size=(500, 500)
)


image.save("aspirin_structure.png")


print("Image saved successfully!")