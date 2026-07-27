import pandas as pd
from rdkit_utils import calculate_properties, lipinski_filter


# Read CSV file
df = pd.read_csv("compounds.csv")


results = []


for index, row in df.iterrows():

    compound_id = row["Compound_ID"]
    smiles = row["SMILES"]

    properties = calculate_properties(smiles)

    if properties:

        properties["Lipinski_Result"] = lipinski_filter(properties)

        properties["Compound_ID"] = compound_id
        properties["SMILES"] = smiles
        
        results.append(properties)


# Convert results to dataframe
result_df = pd.DataFrame(results)


# Save output
result_df.to_csv("compound_properties.csv", index=False)


print("Descriptor calculation completed!")
print(result_df)