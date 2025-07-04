import os
import pandas as pd
import subprocess
from pathlib import Path

# === Configuration ===
WORKING_DIR = "/media/elhabashy/Elements/Hadeer_backup/hybrid_Xlmsecs"
PROTEIN_LIST_PATH = os.path.join(WORKING_DIR, "about_dataset", "mito_proteome", "published_human_mitocarta3.0_impiq2_proteins.csv")
BLAST_BIN = "/Hadeer/software/ncbi-blast-2.9.0+/bin/blastp"
BLAST_DB = "/Hadeer/software/ncbi-blast-2.9.0+/blastdb/pdbaa_28_10_2024/pdbaa"
OUTPUT_DIR = os.path.join(WORKING_DIR, "dataset", "mito_proteome", "human")
OUTPUT_PPI_DIR = "/Hadeer/hybrid_Xlmsecs/dataset/mito_complex/human"
FINAL_CSV = "/Hadeer/hybrid_Xlmsecs/about_dataset/mito_proteome/human_mito_structure_ppi_prefix_3d_6245_Nov2024.csv"

# === Load protein list ===
protein_list = pd.read_csv(PROTEIN_LIST_PATH, comment='#')

# === BLAST for each protein ===
for _, row in protein_list.iterrows():
    uid = row['uid']
    protein_dir = os.path.join(OUTPUT_DIR, uid)
    os.makedirs(protein_dir, exist_ok=True)

    fasta_path = os.path.join(protein_dir, f"{uid}.fasta")
    blast_output_path = os.path.join(protein_dir, f"{uid}_blastp.out")

    # Ensure the fasta file already exists as downloading was commented out in R code
    if os.path.exists(fasta_path):
        try:
            subprocess.run([
                BLAST_BIN,
                "-query", fasta_path,
                "-db", BLAST_DB,
                "-evalue", "0.00001",
                "-out", blast_output_path,
                "-outfmt", "6"
            ], check=True)
        except subprocess.CalledProcessError as e:
            print(f"BLAST error for {uid}: {e}")

# === Parse BLAST outputs and generate structure-based PPIs ===
columns = [
    "query", "subject", "%identity", "alignment length", "mismatches", "gap opens",
    "q. start", "q. end", "s. start", "s. end", "evalue", "bit score"
]

ppi_list = []

for i, row_i in protein_list.iterrows():
    uid_i = row_i['uid']
    path_i = os.path.join(OUTPUT_DIR, uid_i, f"{uid_i}_blastp.out")
    if not os.path.exists(path_i):
        continue

    temp_i = pd.read_csv(path_i, sep='\t', header=None, names=columns)
    temp_i['pdb'] = temp_i['subject'].str.extract(r'_([A-Za-z0-9]{4})')
    temp_i['chain'] = temp_i['subject'].str.extract(r'^([A-Za-z0-9])_')
    temp_i = temp_i.drop_duplicates(subset='subject')
    temp_i = temp_i[temp_i['%identity'] >= 50]

    for j, row_j in protein_list.iterrows():
        uid_j = row_j['uid']
        path_j = os.path.join(OUTPUT_DIR, uid_j, f"{uid_j}_blastp.out")
        if not os.path.exists(path_j):
            continue

        temp_j = pd.read_csv(path_j, sep='\t', header=None, names=columns)
        temp_j['pdb'] = temp_j['subject'].str.extract(r'_([A-Za-z0-9]{4})')
        temp_j['chain'] = temp_j['subject'].str.extract(r'^([A-Za-z0-9])_')
        temp_j = temp_j.drop_duplicates(subset='subject')
        temp_j = temp_j[temp_j['%identity'] >= 50]

        temp_merge = pd.merge(temp_i, temp_j, on='pdb', suffixes=('.x', '.y'))
        temp_merge = temp_merge[temp_merge['subject.x'] != temp_merge['subject.y']]

        temp_merge['pair'] = temp_merge['subject.x'] + "," + temp_merge['subject.y']
        temp_merge['reverse_pair'] = temp_merge['subject.y'] + "," + temp_merge['subject.x']

        temp_merge['delete'] = temp_merge['pair'].isin(temp_merge['reverse_pair'])
        temp_merge = temp_merge[~temp_merge['delete']]

        if not temp_merge.empty:
            output_file = os.path.join(OUTPUT_PPI_DIR, f"{uid_i}-{uid_j}.csv")
            temp_merge.to_csv(output_file, index=False)

            ppi_list.append({
                "interactorA_entry": row_i['Entry'],
                "interactorA_uid": uid_i,
                "interactorA_organism": row_i['Organism'],
                "interactorA_length": row_i['Length'],
                "interactorA_Pfam_domains": row_i['Cross.reference..Pfam.'],
                "interactorB_entry": row_j['Entry'],
                "interactorB_uid": uid_j,
                "interactorB_organism": row_j['Organism'],
                "interactorB_length": row_j['Length'],
                "interactorB_Pfam_domains": row_j['Cross.reference..Pfam.'],
                "complex_3d": len(temp_merge)
            })

# === Export final table ===
pd.DataFrame(ppi_list).to_csv(FINAL_CSV, index=False)

