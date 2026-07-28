from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
import pandas as pd


def remove_fragments(smiles):

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    fragment_chooser = rdMolStandardize.LargestFragmentChooser()

    largest_fragment = fragment_chooser.choose(mol)

    return Chem.MolToSmiles(largest_fragment)


data = {
    "Compound_ID": [
        "CMP001",
        "CMP002",
        "CMP003",
        "CMP004"
    ],

    "SMILES": [
        "CCO.Cl",
        "CC(=O)O.[Na+]",
        "c1ccccc1.Br",
        "CCN"
    ]
}


df = pd.DataFrame(data)

df["Largest_Fragment_SMILES"] = df["SMILES"].apply(remove_fragments)

print(df)

df.to_csv(
    "outputs/day17_fragment_removal.csv",
    index=False
)