import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

# Read the descriptor dataset
df = pd.read_csv("data/curated/descriptor_dataset.csv")

# Create RDKit standardizer
normalizer = rdMolStandardize.Normalizer()

# Function to standardize molecules
def standardize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    # Normalize molecule
    mol = normalizer.normalize(mol)

    # Convert back to canonical SMILES
    return Chem.MolToSmiles(mol, canonical=True)

# Standardize molecules
df["Standardized_SMILES"] = df["Canonical_SMILES"].apply(standardize_smiles)

# Display results
print(df[["Compound_ID", "Canonical_SMILES", "Standardized_SMILES"]])

# Save standardized dataset
df.to_csv("data/curated/standardized_dataset.csv", index=False)

print("\nMolecular standardization completed successfully!")