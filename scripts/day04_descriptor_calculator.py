import csv
from rdkit import Chem
from rdkit.Chem import Descriptors

# Open input CSV
with open("molecules.csv", "r") as infile:

    reader = csv.DictReader(infile)

    # Create output CSV
    with open("drug_properties.csv", "w", newline="") as outfile:

        writer = csv.writer(outfile)

        # Header
        writer.writerow([
            "Compound",
            "SMILES",
            "MolecularWeight",
            "LogP",
            "TPSA",
            "HBD",
            "HBA"
        ])

        # Process every molecule
        for row in reader:

            compound = row["Compound"]
            smiles = row["SMILES"]

            mol = Chem.MolFromSmiles(smiles)

            mw = round(Descriptors.MolWt(mol), 2)
            logp = round(Descriptors.MolLogP(mol), 2)
            tpsa = round(Descriptors.TPSA(mol), 2)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)

            writer.writerow([
                compound,
                smiles,
                mw,
                logp,
                tpsa,
                hbd,
                hba
            ])

print("Drug properties saved successfully!")