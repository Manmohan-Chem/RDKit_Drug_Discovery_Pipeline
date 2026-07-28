from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd


def veber_rule(smiles):

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    tpsa = Descriptors.TPSA(mol)

    veber_pass = (
        rotatable_bonds <= 10 and
        tpsa <= 140
    )

    return {
        "Rotatable_Bonds": rotatable_bonds,
        "TPSA": tpsa,
        "Veber_Pass": veber_pass
    }


data = {
    "Compound_ID": [
        "CMP001",
        "CMP002",
        "CMP003",
        "CMP004"
    ],

    "SMILES": [
        "CCO",
        "CC(=O)O",
        "c1ccccc1",
        "CCN(CC)CC"
    ]
}


df = pd.DataFrame(data)

veber_df = df["SMILES"].apply(veber_rule).apply(pd.Series)

result = pd.concat([df, veber_df], axis=1)

print(result)

result.to_csv(
    "outputs/day21_veber_rule.csv",
    index=False
)