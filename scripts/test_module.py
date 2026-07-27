from rdkit_utils import calculate_properties


smiles = "CC(=O)Oc1ccccc1C(=O)O"


result = calculate_properties(smiles)


print(result)