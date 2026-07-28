from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
import pandas as pd


def neutralize_charge(smiles):

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    uncharger = rdMolStandardize.Uncharger()

    neutral_mol = uncharger.uncharge(mol)

    return Chem.MolToSmiles(neutral_mol)


data = {
    "Compound_ID": [
        "CMP001",
        "CMP002",
        "CMP003",
        "CMP004"
    ],

    "SMILES": [
        "[NH4+]",
        "CC(=O)[O-]",
        "C[NH+](C)C",
        "CCO"
    ]
}


df = pd.DataFrame(data)

df["Neutral_SMILES"] = df["SMILES"].apply(neutralize_charge)

print(df)

df.to_csv(
    "outputs/day16_charge_neutralization.csv",
    index=False
)