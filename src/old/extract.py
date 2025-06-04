import pickle
from itertools import combinations

def extract_accessions_csv(csv_file_path):
    identifiers = []
    with open(csv_file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            fields = line.split(',')
            if len(fields) >= 2:
                versioned_id = fields[1]
                base_id = versioned_id.split('.')[0]
                identifiers.append(base_id)
    return identifiers

import csv

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
    
    ARM10_HUMAN_ARM_REPEAT_PSI3 = extract_accessions_csv("ARM10_HUMAN_ARM_REPEAT_PSI3.csv")
    ARM10_HUMAN_FULL_PSI10 = extract_accessions_csv("ARM10_HUMAN_FULL_PSI10.csv")
    uniprot_text_search = extract_accessions_tsv("uniprot_text_search.tsv", "Entry")
    blast_ARM10_HUMAN = extract_accessions_tsv("blast_ARM10_HUMAN.tsv", "Accession")

    datasets = [ARM10_HUMAN_ARM_REPEAT_PSI3, ARM10_HUMAN_FULL_PSI10, uniprot_text_search, blast_ARM10_HUMAN]


    # Output
    print("Text based search:", len(uniprot_text_search))
    print("PSI-BLAST with full ARM10:", len(ARM10_HUMAN_FULL_PSI10))
    print("PSI BLAST with ARM repeat:", len(ARM10_HUMAN_ARM_REPEAT_PSI3))
    print("blast with ARM10 full:", len(blast_ARM10_HUMAN))
    
    # Find common entries

    common_entries = set(ARM10_HUMAN_ARM_REPEAT_PSI3) & set(ARM10_HUMAN_FULL_PSI10) & set(uniprot_text_search) & set(blast_ARM10_HUMAN)
    total_entries = set(ARM10_HUMAN_ARM_REPEAT_PSI3 + ARM10_HUMAN_FULL_PSI10 + uniprot_text_search + blast_ARM10_HUMAN)
    
    # Output common entries
    print("Common Entries:", common_entries, " ", len(common_entries))
    print("Total unique entries:", len(total_entries))

    print("Unique entries in each dataset:")
    print("ARM10_HUMAN_ARM_REPEAT_PSI3:", len(set(ARM10_HUMAN_ARM_REPEAT_PSI3)))
    print("ARM10_HUMAN_FULL_PSI10:", len(set(ARM10_HUMAN_FULL_PSI10)))
    print("uniprot_text_search:", len(set(uniprot_text_search)))
    print("blast_ARM10_HUMAN:", len(set(blast_ARM10_HUMAN)))
    
    print(total_entries)
    print(len(set(ARM10_HUMAN_FULL_PSI10) & set(blast_ARM10_HUMAN)))
    
    
    with open('out/side1_uniprot_ids.pkl', 'wb') as file:
        pickle.dump(list(total_entries), file)
        
    with open('out/side2_uniprot_ids.pkl', 'wb') as file:
        pickle.dump(["Q92759"], file)
        
if __name__ == "__main__":
    main()
