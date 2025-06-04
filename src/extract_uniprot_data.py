import pickle
from itertools import combinations
import pandas as pd
import csv
import os

def extract_accessions_tsv(tsv_file, column_name):
    """
    Extracts the 'Accession' column from a TSV file.

    Args:
        tsv_file (str): Path to the TSV file.
        column_name (str): The name of the column to extract.

    Returns:
        List[str]: A list of entries from the specified column.
    """
    accession_entries = []

    with open(tsv_file, newline='') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            iso_id = row[column_name]
            base_id = iso_id.split('-')[0] # remove isoform marker
            accession_entries.append(base_id)

    return accession_entries

def main():

    paths = [
        "../datasets_uniprot/uniprotkb_interpro_PFAM_ARM_2025_05_02.tsv"
    ]

    # collect all IDs from the TSV files
    collected_ids = set()
    for path in paths:
        collected_ids.update(extract_accessions_tsv(path, "Entry"))

    pkl_path = "../jobs_AF/side1_uniprot_ids.pkl"
    # load existing IDs (if any)
    try:
        with open(pkl_path, "rb") as f:
            existing_ids = set(pickle.load(f))
    except FileNotFoundError:
        existing_ids = set()

    # compute only the brand‑new IDs
    new_ids = collected_ids - existing_ids
    if not new_ids:
        print("No new IDs to add.")
        return

    # update and write back
    existing_ids.update(new_ids)
    os.makedirs(os.path.dirname(pkl_path), exist_ok=True)
    with open(pkl_path, "wb") as f:
        pickle.dump(list(existing_ids), f)

    print(f"Added {len(new_ids)} new IDs.")


if __name__ == "__main__":
    main()

   
   