import pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

# Read the standardized dataset
df = pd.read_csv("data/curated/standardized_dataset.csv")

# Create Salt Remover
remover = rdMolStandardize.FragmentRemover()

# Function to remove salts
def remove_salts(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    # Remove salts and small fragments
    mol = remover.remove(mol)

    # Convert back to canonical SMILES
    return Chem.MolToSmiles(mol, canonical=True)

# Apply salt removal
df["Salt_Free_SMILES"] = df["Standardized_SMILES"].apply(remove_salts)

# Display results
print(df[["Compound_ID", "Standardized_SMILES", "Salt_Free_SMILES"]])

# Save dataset
df.to_csv("data/curated/salt_removed_dataset.csv", index=False)

print("\nSalt removal completed successfully!")