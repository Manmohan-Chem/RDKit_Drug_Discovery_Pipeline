from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd


compounds = {
    "CMP001": "CCO",
    "CMP002": "CC(=O)O",
    "CMP003": "c1ccccc1",
    "CMP004": "CCOC(=O)c1ccccc1"
}


results = []

for cid, smi in compounds.items():

    mol = Chem.MolFromSmiles(smi)

    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol,
        radius=2,
        nBits=1024
    )

    bits = list(fp)

    results.append({
        "Compound_ID": cid,
        "SMILES": smi,
        "Fingerprint": bits
    })


df = pd.DataFrame(results)

df.to_csv(
    "outputs/day24_morgan_fingerprint.csv",
    index=False
)

print(df)
print("Day24 fingerprint file saved successfully")