import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina
from rdkit import DataStructs

# Read the descriptor dataset
df = pd.read_csv("data/curated/descriptor_dataset.csv")

# Morgan Fingerprint generator
fp_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

# Generate fingerprints
fingerprints = []

for smiles in df["Canonical_SMILES"]:
    mol = Chem.MolFromSmiles(smiles)
    fp = fp_generator.GetFingerprint(mol)
    fingerprints.append(fp)

# Calculate distance matrix
distance_matrix = []

for i in range(1, len(fingerprints)):
    similarities = DataStructs.BulkTanimotoSimilarity(
        fingerprints[i],
        fingerprints[:i]
    )

    distance_matrix.extend([1 - x for x in similarities])

# Perform Butina clustering
clusters = Butina.ClusterData(
    distance_matrix,
    len(fingerprints),
    0.4,
    isDistData=True
)

# Assign cluster IDs
cluster_ids = [0] * len(df)

for cluster_number, cluster in enumerate(clusters):
    for molecule_index in cluster:
        cluster_ids[molecule_index] = cluster_number + 1

df["Cluster_ID"] = cluster_ids

# Display results
print(df[["Compound_ID", "Canonical_SMILES", "Cluster_ID"]])

# Save clustered dataset
df.to_csv("data/curated/clustered_dataset.csv", index=False)

print("\nClustered dataset saved successfully!")