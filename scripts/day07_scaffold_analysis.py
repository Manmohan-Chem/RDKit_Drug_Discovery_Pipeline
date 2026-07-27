import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

# Read the curated dataset
df = pd.read_csv("data/curated/curated_molecules.csv")

# Function to extract the Bemis-Murcko scaffold
def get_scaffold(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    scaffold = MurckoScaffold.GetScaffoldForMol(mol)

    return Chem.MolToSmiles(scaffold)

# Generate scaffold column
df["Scaffold"] = df["Canonical_SMILES"].apply(get_scaffold)

# Display the dataset
print(df)

# Save the scaffold dataset
df.to_csv("data/curated/scaffold_dataset.csv", index=False)

print("\nScaffold dataset saved successfully!")