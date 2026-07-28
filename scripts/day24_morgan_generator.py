from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator


compounds = {
    "CMP001": "CCO",
    "CMP002": "CC(=O)O",
    "CMP003": "c1ccccc1",
    "CMP004": "CCOC(=O)c1ccccc1"
}


morgan_generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=2,
    fpSize=1024
)


for cid, smi in compounds.items():

    mol = Chem.MolFromSmiles(smi)

    fp = morgan_generator.GetFingerprint(mol)

    print(cid, smi)
    print(fp)
    print()