import csv

def extract_accessions(tsv_file):
    """
    Extracts the 'Accession' column from a TSV file.

    Args:
        tsv_file (str): Path to the TSV file.

    Returns:
        List[str]: A list of accession entries.
    """
    accession_entries = []

    with open(tsv_file, newline='') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            accession_entries.append(row['Accession'])

    return accession_entries

