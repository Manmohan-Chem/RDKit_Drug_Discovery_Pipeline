import pandas as pd
from rdkit_utils import calculate_similarity


# Query molecule (Aspirin)
query_smiles = "CC(=O)Oc1ccccc1C(=O)O"


# Read compound library
df = pd.read_csv("compound_library.csv")


results = []


for index, row in df.iterrows():

    compound_id = row["Compound_ID"]
    smiles = row["SMILES"]

    similarity = calculate_similarity(
        query_smiles,
        smiles
    )

    results.append({
        "Compound_ID": compound_id,
        "SMILES": smiles,
        "Similarity": similarity
    })


# Convert to dataframe
result_df = pd.DataFrame(results)


# Sort highest similarity first
result_df = result_df.sort_values(
    by="Similarity",
    ascending=False
)


# Save results
result_df.to_csv(
    "similarity_results.csv",
    index=False
)


print(result_df)