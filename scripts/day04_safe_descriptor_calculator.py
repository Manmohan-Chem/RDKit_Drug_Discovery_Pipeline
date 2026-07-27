import csv
from rdkit import Chem
from rdkit.Chem import Descriptors

with open("molecules.csv", "r") as infile:

    reader = csv.DictReader(infile)

    with open("drug_properties_safe.csv", "w", newline="") as outfile:

        writer = csv.writer(outfile)

        writer.writerow([
            "Compound",
            "SMILES",
            "MolecularWeight",
            "LogP",
            "TPSA",
            "HBD",
            "HBA",
            "Status"
        ])

        for row in reader:

            compound = row["Compound"]
            smiles = row["SMILES"]

            mol = Chem.MolFromSmiles(smiles)

            if mol is None:
                writer.writerow([
                    compound,
                    smiles,
                    "",
                    "",
                    "",
                    "",
                    "",
                    "Invalid SMILES"
                ])
                continue

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
                hba,
                "Valid"
            ])

print("Analysis completed successfully!")