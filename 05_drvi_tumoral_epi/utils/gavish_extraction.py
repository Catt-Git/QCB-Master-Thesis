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
# MP1 - CELL CYCLE G2/M: NON VIENE DAL CSV
# =========================
# GAVISH.csv ha 40 colonne e non 41: MP1 (Cell Cycle - G2/M) non c'e', ne' nel CSV ne'
# nella versione MSigDB da cui il CSV e' stato esportato (i set si chiamano
# GAVISH_3CA_MALIGNANT_METAPROGRAM_<n>_<NOME>, e quello con n=1 non esiste). Senza MP1 una
# dimensione latente specificamente G2/M non ha metaprogramma da matchare e legge come MP2
# o come niente, quindi la lista e' presa dalla sorgente degli autori:
#
#   https://github.com/tiroshlab/3ca -> ITH_hallmarks/MPs_distribution/MP_list.RDS.gz
#   readRDS(...)$Cancer[["MP1  Cell Cycle - G2/M"]]   -> 50 geni, nell'ordine del paper
#
# CONTROLLO DI PROVENIENZA: lo stesso oggetto contiene anche gli altri 40, e il suo MP2 e'
# identico gene-per-gene a MP2_CELL_CYCLE_G1_S.txt tranne due simboli aggiornati da MSigDB
# (HIST1H4C -> H4C3, KIAA0101 -> PCLAF). Le due sorgenti sono quindi la stessa lista.
#
# I SIMBOLI SONO QUELLI DEL PAPER, DI PROPOSITO. Applicando la stessa regola di MSigDB,
# l'unico gene che cambierebbe qui e' HIST1H4C -> H4C3; ma shiao.h5ad e' su una reference
# piu' vecchia e non ha H4C3 (per questo MP2 e MP3 su disco perdono un gene a testa nella
# coverage table). Cosi' come e' scritta, MP1 mappa 50/50 su questo dataset.
MP1_CELL_CYCLE_G2_M = [
    "TOP2A", "UBE2C", "HMGB2", "NUSAP1", "CENPF", "CCNB1", "TPX2", "CKS2", "BIRC5", "PRC1",
    "PTTG1", "KPNA2", "MKI67", "CDC20", "CDK1", "CCNB2", "CDKN3", "SMC4", "NUF2", "ARL6IP1",
    "CKAP2", "ASPM", "PLK1", "CKS1B", "CCNA2", "AURKA", "MAD2L1", "GTSE1", "HMMR", "UBE2T",
    "CENPE", "CENPA", "KIF20B", "AURKB", "CDCA3", "CDCA8", "UBE2S", "KNSTRN", "KIF2C", "PBK",
    "TUBA1B", "DLGAP5", "TACC3", "STMN1", "DEPDC1", "ECT2", "CENPW", "ZWINT", "HIST1H4C",
    "KIF23",
]

(output_dir / "MP1_CELL_CYCLE_G2_M.txt").write_text("\n".join(MP1_CELL_CYCLE_G2_M) + "\n")

# =========================
# CHECK
# =========================
# 41 e non 40: le 40 colonne del CSV piu' MP1, che il blocco sopra scrive dalla lista
# hard-coded perche' nel CSV non c'e'.
expected_programs = 41

created_files = list(output_dir.glob("*.txt"))
created_programs = len(created_files)

print("\n==============================")
print("Gavish metaprograms - CHECK")
print("==============================")
print(f"Programmi nel CSV : {len(df.columns)}")
print(f"MP1 (fuori CSV)   : 1")
print(f"File .txt creati  : {created_programs}")
print(f"File attesi       : {expected_programs}")

if created_programs == expected_programs:
    print("\n✓ OK: tutti i 41 programmi sono stati creati.")
else:
    print("\n✗ ATTENZIONE: il numero di file non corrisponde a 41.")

    # Mostra i programmi presenti nel CSV
    print("\nProgrammi nel CSV:")
    for program in df.columns:
        print(f"  - {program}")

    # Mostra i file creati
    print("\nFile creati:")
    for file in sorted(created_files):
        print(f"  - {file.name}")
