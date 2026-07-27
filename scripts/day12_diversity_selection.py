import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker

# Read the descriptor dataset
df = pd.read_csv("data/curated/descriptor_dataset.csv")

# Create Morgan Fingerprint generator
fp_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

# Generate fingerprints
fingerprints = []

for smiles in df["Canonical_SMILES"]:
    mol = Chem.MolFromSmiles(smiles)
    fp = fp_generator.GetFingerprint(mol)
    fingerprints.append(fp)

# Distance function
def distance(i, j):
    return 1 - DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])

# Create diversity picker
picker = MaxMinPicker()

# Number of molecules to select
num_to_pick = min(3, len(fingerprints))

# Select diverse molecules
picked = picker.LazyPick(
    distance,
    len(fingerprints),
    num_to_pick,
    seed=42
)

# Display selected molecules
selected_df = df.iloc[list(picked)]

print(selected_df[["Compound_ID", "Canonical_SMILES"]])

# Save selected molecules
selected_df.to_csv("data/curated/diverse_molecules.csv", index=False)

print("\nDiverse molecule selection completed successfully!")