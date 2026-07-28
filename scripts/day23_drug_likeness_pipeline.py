from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen


# -----------------------------
# Descriptor calculation
# -----------------------------

def calculate_descriptors(smiles):

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    descriptors = {

        "MW": round(Descriptors.MolWt(mol), 2),

        "LogP": round(Crippen.MolLogP(mol), 2),

        "HBD": Lipinski.NumHDonors(mol),

        "HBA": Lipinski.NumHAcceptors(mol),

        "TPSA": round(
            Descriptors.TPSA(mol), 
            2
        ),

        "Rotatable_Bonds":
            Lipinski.NumRotatableBonds(mol)
    }

    return descriptors



# -----------------------------
# Lipinski Rule
# -----------------------------

def lipinski_rule(desc):

    violations = 0

    if desc["MW"] > 500:
        violations += 1

    if desc["LogP"] > 5:
        violations += 1

    if desc["HBD"] > 5:
        violations += 1

    if desc["HBA"] > 10:
        violations += 1


    return violations <= 1



# -----------------------------
# Veber Rule
# -----------------------------

def veber_rule(desc):

    if desc["Rotatable_Bonds"] <= 10 and desc["TPSA"] <= 140:
        return True

    return False



# -----------------------------
# Main Pipeline
# -----------------------------

compounds = {

    "CMP001": "CCO",               # Ethanol

    "CMP002": "CC(=O)O",           # Acetic acid

    "CMP003": "c1ccccc1",          # Benzene

    "CMP004": "CCOC(=O)c1ccccc1"   # Ethyl benzoate
}



results = []


for cid, smiles in compounds.items():

    desc = calculate_descriptors(smiles)

    if desc:

        lipinski = lipinski_rule(desc)

        veber = veber_rule(desc)


        drug_like = lipinski and veber


        results.append({

            "Compound_ID": cid,

            "SMILES": smiles,

            **desc,

            "Lipinski_Pass": lipinski,

            "Veber_Pass": veber,

            "Drug_Like": drug_like
        })



# -----------------------------
# Display result
# -----------------------------

import pandas as pd


df = pd.DataFrame(results)

print(df)


# Save output

df.to_csv(
    "outputs/day23_drug_likeness_report.csv",
    index=False
)

print("Day23 report saved successfully")