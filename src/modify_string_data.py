#!/usr/bin/env python3

import os

def load_string_to_uniprot_mapping(mapping_file):
    """
    Load the mapping from STRING ID to UniProt ID from the mapping file.
    
    Args:
        mapping_file (str): Path to the mapping file
        
    Returns:
        dict: Dictionary mapping STRING IDs to UniProt IDs
    """
    string_to_uniprot = {}
    try:
        with open(mapping_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    string_id, uniprot_id = parts
                    string_to_uniprot[string_id] = uniprot_id
        print(f"Loaded {len(string_to_uniprot)} STRING to UniProt mappings")
    except FileNotFoundError:
        print(f"Error: Mapping file {mapping_file} not found")
        exit(1)
    
    return string_to_uniprot

def modify_string_file(input_file, output_file, string_to_uniprot):
    """
    Modify the STRING file by adding UniProt IDs after each STRING ID.
    
    Args:
        input_file (str): Path to the input STRING file
        output_file (str): Path to the output modified file
        string_to_uniprot (dict): Dictionary mapping STRING IDs to UniProt IDs
    """
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            # Read the header line
            header = infile.readline().strip()
            # Write the modified header
            outfile.write(header + "\n")
            
            # Counter for stats
            total_lines = 0
            modified_lines = 0
            missing_mapping = set()
            
            # Process each line
            for line in infile:
                total_lines += 1
                parts = line.strip().split()
                if len(parts) >= 2:
                    string_id1 = parts[0]
                    string_id2 = parts[1]
                    
                    # Get the UniProt IDs for each STRING ID
                    uniprot_id1 = string_to_uniprot.get(string_id1, "UNKNOWN")
                    uniprot_id2 = string_to_uniprot.get(string_id2, "UNKNOWN")
                    
                    # Track missing mappings
                    if uniprot_id1 == "UNKNOWN":
                        missing_mapping.add(string_id1)
                    if uniprot_id2 == "UNKNOWN":
                        missing_mapping.add(string_id2)
                    
                    # Create the modified line with UniProt IDs added after each STRING ID
                    modified_parts = [string_id1, uniprot_id1, string_id2, uniprot_id2] + parts[2:]
                    modified_line = "\t".join(modified_parts)
                    outfile.write(modified_line + "\n")
                    modified_lines += 1
                else:
                    # If the line doesn't have at least two items, write it unchanged
                    outfile.write(line)
            
            print(f"Processed {total_lines} lines")
            print(f"Modified {modified_lines} lines")
            print(f"Found {len(missing_mapping)} STRING IDs without UniProt mapping")
    except FileNotFoundError:
        print(f"Error: Input file {input_file} not found")
        exit(1)

def main():
    # File paths
    string_file = "data/STRING/9606.protein.links.full.v12.0.txt"
    mapping_file = "data/STRING/string_to_uniprot_mapping.txt"
    output_file = "data/STRING/9606.protein.links.full.v12.0_withPairID.txt"
    
    # Ensure the input files exist
    if not os.path.exists(string_file):
        print(f"Error: STRING file {string_file} not found")
        exit(1)
    
    if not os.path.exists(mapping_file):
        print(f"Error: Mapping file {mapping_file} not found")
        exit(1)
    
    # Load the mapping
    string_to_uniprot = load_string_to_uniprot_mapping(mapping_file)
    
    # Modify the STRING file
    modify_string_file(string_file, output_file, string_to_uniprot)
    
    print(f"Modified file saved as {output_file}")

if __name__ == "__main__":
    main()