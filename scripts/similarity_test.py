from rdkit_utils import calculate_similarity


compound1 = "CC(=O)Oc1ccccc1C(=O)O"

compound2 = "CCOC(=O)c1ccccc1"


score = calculate_similarity(
    compound1,
    compound2
)


print("Similarity:", score)