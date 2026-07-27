import pandas as pd

from rdkit import Chem
from rdkit.Chem import Draw

from rdkit_utils import calculate_properties, lipinski_filter


# Read compound library

df = pd.read_csv("compound_library.csv")


results = []


for index, row in df.iterrows():

    compound_id = row["Compound_ID"]
    smiles = row["SMILES"]


    # Calculate properties
    properties = calculate_properties(smiles)


    if properties:

        properties["Compound_ID"] = compound_id
        properties["SMILES"] = smiles
        properties["Lipinski_Result"] = lipinski_filter(properties)


        # Create molecule image

        mol = Chem.MolFromSmiles(smiles)

        image = Draw.MolToImage(
            mol,
            size=(300,300)
        )

        image.save(
            f"{compound_id}.png"
        )


        results.append(properties)



# Save report

report = pd.DataFrame(results)

report.to_csv(
    "compound_report.csv",
    index=False
)


print("Compound report generated!")