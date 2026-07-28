from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity


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


similarity = TanimotoSimilarity(
    fingerprints["CMP001"],
    fingerprints["CMP002"]
)

print(
    "Similarity between CMP001 and CMP002:",
    similarity
)