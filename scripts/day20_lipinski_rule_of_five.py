from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd


def lipinski_rule(smiles):

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    molecular_weight = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    h_bond_donors = Descriptors.NumHDonors(mol)
    h_bond_acceptors = Descriptors.NumHAcceptors(mol)

    lipinski_pass = (
        molecular_weight <= 500 and
        logp <= 5 and
        h_bond_donors <= 5 and
        h_bond_acceptors <= 10
    )

    return {
        "Molecular_Weight": molecular_weight,
        "LogP": logp,
        "H_Bond_Donors": h_bond_donors,
        "H_Bond_Acceptors": h_bond_acceptors,
        "Lipinski_Pass": lipinski_pass
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

lipinski_df = df["SMILES"].apply(lipinski_rule).apply(pd.Series)

result = pd.concat([df, lipinski_df], axis=1)

print(result)

result.to_csv(
    "outputs/day20_lipinski_rule_of_five.csv",
    index=False
)