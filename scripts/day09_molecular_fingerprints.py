import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

# Read the descriptor dataset
df = pd.read_csv("data/curated/descriptor_dataset.csv")

# Create Morgan Fingerprint generator
fp_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

# Function to calculate fingerprint
def calculate_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    fp = fp_generator.GetFingerprint(mol)

    return fp.ToBitString()

# Generate fingerprints
df["Morgan_Fingerprint"] = df["Canonical_SMILES"].apply(calculate_fingerprint)

# Display the first few rows
print(df[["Compound_ID", "Canonical_SMILES", "Morgan_Fingerprint"]])

# Save the dataset
df.to_csv("data/curated/fingerprint_dataset.csv", index=False)

print("\nFingerprint dataset saved successfully!")