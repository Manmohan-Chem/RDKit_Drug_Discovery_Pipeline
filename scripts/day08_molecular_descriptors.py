import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

# Read the scaffold dataset
df = pd.read_csv("data/curated/scaffold_dataset.csv")

# Function to calculate molecular descriptors
def calculate_descriptors(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return pd.Series({
            "Molecular_Weight": None,
            "LogP": None,
            "TPSA": None,
            "HBD": None,
            "HBA": None
        })

    return pd.Series({
        "Molecular_Weight": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol)
    })

# Calculate descriptors
descriptor_df = df["Canonical_SMILES"].apply(calculate_descriptors)

# Merge descriptors with original dataset
df = pd.concat([df, descriptor_df], axis=1)

# Display the dataset
print(df)

# Save the descriptor dataset
df.to_csv("data/curated/descriptor_dataset.csv", index=False)

print("\nDescriptor dataset saved successfully!")