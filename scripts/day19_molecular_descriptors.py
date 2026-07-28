from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd


def calculate_descriptors(smiles):

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    return {
        "Molecular_Weight": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "H_Bond_Donors": Descriptors.NumHDonors(mol),
        "H_Bond_Acceptors": Descriptors.NumHAcceptors(mol),
        "Rotatable_Bonds": Descriptors.NumRotatableBonds(mol)
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

descriptor_df = df["SMILES"].apply(calculate_descriptors).apply(pd.Series)

result = pd.concat([df, descriptor_df], axis=1)

print(result)

result.to_csv(
    "outputs/day19_molecular_descriptors.csv",
    index=False
)