import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

# Read the descriptor dataset
df = pd.read_csv("data/curated/descriptor_dataset.csv")

# Create Morgan Fingerprint generator
fp_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

# Function to generate RDKit fingerprint object
def get_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    return fp_generator.GetFingerprint(mol)

# Select the first molecule as the query
query_smiles = df.loc[0, "Canonical_SMILES"]
query_fp = get_fingerprint(query_smiles)

similarities = []

# Compare the query molecule with every molecule
for smiles in df["Canonical_SMILES"]:
    fp = get_fingerprint(smiles)

    similarity = DataStructs.TanimotoSimilarity(query_fp, fp)

    similarities.append(similarity)

# Add similarity scores to the dataset
df["Tanimoto_Similarity"] = similarities

# Sort by similarity (highest first)
df = df.sort_values(by="Tanimoto_Similarity", ascending=False)

# Display the results
print(df[["Compound_ID", "Canonical_SMILES", "Tanimoto_Similarity"]])

# Save the similarity results
df.to_csv("data/curated/similarity_dataset.csv", index=False)

print("\nSimilarity dataset saved successfully!")