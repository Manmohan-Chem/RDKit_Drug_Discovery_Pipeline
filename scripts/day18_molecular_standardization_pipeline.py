from rdkit import Chem
from rdkit.Chem import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
import pandas as pd


def standardize_molecule(smiles):

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    # Salt Removal
    remover = SaltRemover.SaltRemover()
    mol = remover.StripMol(mol)

    # Tautomer Canonicalization
    enumerator = rdMolStandardize.TautomerEnumerator()
    mol = enumerator.Canonicalize(mol)

    # Charge Neutralization
    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(mol)

    # Fragment Removal
    fragment_chooser = rdMolStandardize.LargestFragmentChooser()
    mol = fragment_chooser.choose(mol)

    return Chem.MolToSmiles(mol)


data = {
    "Compound_ID": [
        "CMP001",
        "CMP002",
        "CMP003",
        "CMP004"
    ],

    "SMILES": [
        "CCO.Cl",
        "CC(=O)[O-].[Na+]",
        "O=C1NC=CC=C1",
        "C1=CC=C(C=C1)O"
    ]
}


df = pd.DataFrame(data)

df["Standardized_SMILES"] = df["SMILES"].apply(standardize_molecule)

print(df)

df.to_csv(
    "outputs/day18_molecular_standardization_pipeline.csv",
    index=False
)