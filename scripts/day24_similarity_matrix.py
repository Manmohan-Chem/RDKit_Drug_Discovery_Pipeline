from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity
import pandas as pd


compounds = {
    "CMP001": "CCO",
    "CMP002": "CC(=O)O",
    "CMP003": "c1ccccc1",
    "CMP004": "CCOC(=O)c1ccccc1"
}


fingerprints = {}

for cid, smi in compounds.items():

    mol = Chem.MolFromSmiles(smi)

    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol,
        radius=2,
        nBits=1024
    )

    fingerprints[cid] = fp


matrix = []

ids = list(compounds.keys())

for cid1 in ids:

    row = []

    for cid2 in ids:

        sim = TanimotoSimilarity(
            fingerprints[cid1],
            fingerprints[cid2]
        )

        row.append(round(sim, 3))

    matrix.append(row)


df = pd.DataFrame(
    matrix,
    index=ids,
    columns=ids
)


df.to_csv(
    "outputs/day24_tanimoto_similarity_matrix.csv"
)


print(df)
print("Similarity matrix saved successfully")