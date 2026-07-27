from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
import pandas as pd


def canonicalize_tautomer(smiles):

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    enumerator = rdMolStandardize.TautomerEnumerator()

    canonical_mol = enumerator.Canonicalize(mol)

    return Chem.MolToSmiles(canonical_mol)


data = {
    "Compound_ID": [
        "CMP001",
        "CMP002",
        "CMP003",
        "CMP004"
    ],

    "SMILES": [
        "O=C1NC=CC=C1",
        "OC1=CC=CC=C1",
        "C1=CC=C(C=C1)O",
        "CC(=O)C"
    ]
}


df = pd.DataFrame(data)


df["Canonical_Tautomer_SMILES"] = df["SMILES"].apply(
    canonicalize_tautomer
)


print(df)


df.to_csv(
    "outputs/day15_tautomer_canonicalization.csv",
    index=False
)