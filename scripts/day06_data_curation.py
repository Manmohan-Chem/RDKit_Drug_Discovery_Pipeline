import pandas as pd
from rdkit import Chem

# Read the input dataset
df = pd.read_csv("data/raw/sample_molecules.csv")

print(df)

# Function to check if a SMILES string is valid
def is_valid_smiles(smiles):
    if pd.isna(smiles):
        return False

    mol = Chem.MolFromSmiles(smiles)

    return mol is not None

# Function to convert SMILES to canonical SMILES
def canonical_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol)

# Create a new column indicating whether each SMILES is valid
df["Valid_SMILES"] = df["SMILES"].apply(is_valid_smiles)

# Display the updated dataset
print(df)

# Keep only valid molecules
df = df[df["Valid_SMILES"] == True]

print("\nAfter removing invalid molecules:")
print(df)

# Count the number of valid molecules
print("\nNumber of valid molecules:", len(df))

# Remove duplicate molecules based on the SMILES column
df = df.drop_duplicates(subset="SMILES")

# Display the dataset after removing duplicates
print("\nAfter removing duplicate molecules:")
print(df)

# Count the remaining unique molecules
print("\nNumber of unique molecules:", len(df))

# Generate canonical SMILES
df["Canonical_SMILES"] = df["SMILES"].apply(canonical_smiles)

print("\nDataset with Canonical SMILES:")
print(df)

# Save the curated dataset
df.to_csv("data/curated/curated_molecules.csv", index=False)

print("\nCurated dataset saved successfully!")