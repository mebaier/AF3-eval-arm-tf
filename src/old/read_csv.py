# This script reads a CSV file line by line, extracts versioned identifiers from the second field 
# (like "Q8N2F6.1" from "sp|Q8N2F6|ARM10_HUMAN,Q8N2F6.1,..."), strips the version suffix (".1"), 
# and collects all identifiers into a Python list (allowing duplicates).

def extract_identifiers(csv_file_path):
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

# Example usage:
if __name__ == "__main__":
    csv_path = "your_file.csv"  # Replace with your CSV file path
    ids = extract_identifiers(csv_path)
    print(ids)
