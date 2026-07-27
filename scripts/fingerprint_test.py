from rdkit_utils import generate_fingerprint


smiles = "CC(=O)Oc1ccccc1C(=O)O"


fp = generate_fingerprint(smiles)


print(fp)