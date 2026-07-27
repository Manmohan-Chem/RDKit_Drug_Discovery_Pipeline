from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem


def calculate_properties(smiles):
    """
    Calculate basic molecular properties from SMILES
    """

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    properties = {
        "Molecular Weight": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "HBD": Descriptors.NumHDonors(mol),
        "HBA": Descriptors.NumHAcceptors(mol)
    }

    return properties



def lipinski_filter(properties):

    violations = 0

    if properties["Molecular Weight"] > 500:
        violations += 1

    if properties["LogP"] > 5:
        violations += 1

    if properties["HBD"] > 5:
        violations += 1

    if properties["HBA"] > 10:
        violations += 1


    if violations <= 1:
        return "Drug-like"
    else:
        return "Poor drug-like"



def generate_fingerprint(smiles):
    """
    Generate Morgan fingerprint from SMILES
    """

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    fingerprint = AllChem.GetMorganFingerprintAsBitVect(
        mol,
        radius=2,
        nBits=2048
    )

    return fingerprint

def calculate_similarity(smiles1, smiles2):
    """
    Calculate Tanimoto similarity between two molecules
    """

    fp1 = generate_fingerprint(smiles1)
    fp2 = generate_fingerprint(smiles2)

    if fp1 is None or fp2 is None:
        return None

    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)

    return similarity