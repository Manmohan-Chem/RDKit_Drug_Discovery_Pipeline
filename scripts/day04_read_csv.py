import csv

print("=" * 50)
print("Reading Molecular Dataset")
print("=" * 50)

with open("molecules.csv", "r") as file:

    reader = csv.DictReader(file)

    for row in reader:

        print("Compound :", row["Compound"])
        print("SMILES   :", row["SMILES"])
        print("-" * 40)