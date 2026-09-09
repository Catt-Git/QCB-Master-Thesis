import pandas as pd
from pathlib import Path
import os
# =========================
# INPUT
# =========================

input_csv = "/home/albertoc/Desktop/QCB-Master-Thesis/datasets/GAVISH.csv"   # <-- cambia con il nome del tuo CSV
output_dir = "/home/albertoc/Desktop/QCB-Master-Thesis/datasets/signatures/GAVISH_metaprograms"  # <-- cambia con il nome della cartella di output

# =========================
# READ CSV
# =========================
df = pd.read_csv(input_csv)

# =========================
# CREATE OUTPUT DIRECTORY
# =========================
output_dir = Path(output_dir)
output_dir.mkdir(parents=True, exist_ok=True)

# =========================
# CREATE ONE TXT PER PROGRAM
# =========================
for program in df.columns:

    # Prendi i geni della colonna
    genes = df[program].dropna().astype(str).str.strip()

    # Rimuovi eventuali celle vuote
    genes = genes[genes != ""]

    # Nome file
    output_file = output_dir / f"{program}.txt"

    # Scrivi un gene per riga
    genes.to_csv(output_file, index=False, header=False)

# =========================
# CHECK
# =========================
expected_programs = 40

created_files = list(output_dir.glob("*.txt"))
created_programs = len(created_files)

print("\n==============================")
print("Gavish metaprograms - CHECK")
print("==============================")
print(f"Programmi nel CSV : {len(df.columns)}")
print(f"File .txt creati  : {created_programs}")
print(f"File attesi       : {expected_programs}")

if created_programs == expected_programs:
    print("\n✓ OK: tutti i 40 programmi sono stati creati.")
else:
    print("\n✗ ATTENZIONE: il numero di file non corrisponde a 40.")

    # Mostra i programmi presenti nel CSV
    print("\nProgrammi nel CSV:")
    for program in df.columns:
        print(f"  - {program}")

    # Mostra i file creati
    print("\nFile creati:")
    for file in sorted(created_files):
        print(f"  - {file.name}")