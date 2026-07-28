from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd


def lipinski_rule(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)

    passed = (
        mw <= 500 and
        logp <= 5 and
        hbd <= 5 and
        hba <= 10
    )

    return {
        "MW": round(mw, 2),
        "LogP": round(logp, 2),
        "HBD": hbd,
        "HBA": hba,
        "Lipinski_Pass": passed
    }


# Input compounds
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
        "CC(=O)Oc1ccccc1C(=O)O"
    ]
}


df = pd.DataFrame(data)


results = df["SMILES"].apply(lipinski_rule)

properties = pd.DataFrame(results.tolist())

final_df = pd.concat([df, properties], axis=1)


print(final_df)


# Save output
final_df.to_csv(
    "outputs/day22_lipinski_rule.csv",
    index=False
)