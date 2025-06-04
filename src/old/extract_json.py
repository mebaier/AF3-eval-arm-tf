""" """

import json

def extract_keys(data, key):
    """Recursively extract all values associated with key in nested JSON."""
    results = []
    if isinstance(data, dict):
        for k, v in data.items():
            if k == key:
                results.append(v)
            else:
                results.extend(extract_keys(v, key))
    elif isinstance(data, list):
        for item in data:
            results.extend(extract_keys(item, key))
    return results

def main():
    # Paths to your JSON files
    file1_path = 'uniprotkb_armadillo_repeat_AND_taxonomy_2025_04_28.json'
    file2_path = 'ncbiblast-R20250428-103533-0994-84282483-p1m.json'

    # Load JSON data
    with open(file1_path, 'r') as f1:
        data1 = json.load(f1)

    with open(file2_path, 'r') as f2:
        data2 = json.load(f2)

    # Extract IDs
    primary_accessions = extract_keys(data1, 'primaryAccession')
    hit_accs = extract_keys(data2, 'hit_acc')

    # Output
    print("Text based search:", primary_accessions)
    print("BLAST with ARM10:", hit_accs)
    
    # Find common entries
    common_entries = set(primary_accessions) & set(hit_accs)

    # Output common entries
    print("Common Entries:", common_entries, " ", len(common_entries))

if __name__ == "__main__":
    main()
