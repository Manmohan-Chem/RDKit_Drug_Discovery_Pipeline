import csv
from rdkit import Chem
from rdkit.Chem import Descriptors

with open("molecules.csv", "r") as infile:

    reader = csv.DictReader(infile)

    with open("lipinski_report.csv", "w", newline="") as outfile:

        writer = csv.writer(outfile)

        writer.writerow([
            "Compound",
            "MW",
            "LogP",
            "HBD",
            "HBA",
            "Lipinski_Result"
        ])

        for row in reader:

            compound = row["Compound"]
            smiles = row["SMILES"]

            mol = Chem.MolFromSmiles(smiles)

            if mol is None:
                writer.writerow([
                    compound,
                    "",
                    "",
                    "",
                    "",
                    "Invalid SMILES"
                ])
                continue

            mw = round(Descriptors.MolWt(mol), 2)
            logp = round(Descriptors.MolLogP(mol), 2)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)

            if mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10:
                result = "Pass"
            else:
                result = "Fail"

            writer.writerow([
                compound,
                mw,
                logp,
                hbd,
                hba,
                result
            ])

print("Lipinski report generated successfully!")