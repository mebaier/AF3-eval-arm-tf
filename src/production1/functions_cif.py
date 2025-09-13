import sys
from collections import defaultdict
import re

def parse_cif_file(cif_file):
    """Parse CIF file and extract atom_site and entity_poly_seq data."""
    
    atom_site_data = []
    entity_poly_seq_data = []
    
    current_section = None
    
    with open(cif_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Detect section headers
            if line.startswith('_atom_site.'):
                current_section = 'atom_site'
                continue
            elif line.startswith('_entity_poly_seq.'):
                current_section = 'entity_poly_seq'
                continue
            elif line.startswith('_') or line.startswith('#') or line.startswith('loop_'):
                if not line.startswith('_atom_site.') and not line.startswith('_entity_poly_seq.'):
                    current_section = None
                continue
            
            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue
            
            # Parse data lines
            if current_section == 'atom_site':
                # Split line respecting quoted strings
                parts = parse_cif_line(line)
                if len(parts) >= 20:  # Ensure we have enough columns
                    atom_site_data.append(parts)
            elif current_section == 'entity_poly_seq':
                parts = parse_cif_line(line)
                if len(parts) >= 3:  # Ensure we have enough columns
                    entity_poly_seq_data.append(parts)
    
    return atom_site_data, entity_poly_seq_data

def parse_cif_line(line):
    """Parse a CIF data line, handling quoted strings."""
    parts = []
    current_part = ""
    in_quotes = False
    
    i = 0
    while i < len(line):
        char = line[i]
        
        if char == '"' or char == "'":
            if not in_quotes:
                in_quotes = True
                quote_char = char
            elif char == quote_char:
                in_quotes = False
            current_part += char
        elif char.isspace() and not in_quotes:
            if current_part:
                parts.append(current_part.strip('"\''))
                current_part = ""
        else:
            current_part += char
        
        i += 1
    
    if current_part:
        parts.append(current_part.strip('"\''))
    
    return parts

def extract_observed_residues(atom_site_data):
    """Extract observed residues from atom_site data."""
    # CIF atom_site columns (approximate positions, may vary by file)
    # Typical columns: group_PDB, id, type_symbol, label_atom_id, label_alt_id,
    # label_comp_id, label_asym_id, label_entity_id, label_seq_id, pdbx_PDB_ins_code,
    # Cartn_x, Cartn_y, Cartn_z, occupancy, B_iso_or_equiv, pdbx_formal_charge,
    # auth_seq_id, auth_comp_id, auth_asym_id, auth_atom_id, pdbx_PDB_model_num
    
    observed_residues = defaultdict(set)
    
    for row in atom_site_data:
        if len(row) >= 20:
            try:
                label_entity_id = row[7]  # label_entity_id
                label_seq_id = row[8]     # label_seq_id (sequence position)
                auth_seq_id = row[16]     # auth_seq_id (author sequence id)
                auth_asym_id = row[18]    # auth_asym_id (author chain identifier)
                label_comp_id = row[5]    # residue type
                
                # Skip non-polymer atoms (like water, ions)
                if label_entity_id and label_seq_id and label_seq_id != '.' and label_seq_id != '?':
                    observed_residues[auth_asym_id].add((int(label_seq_id), auth_seq_id, label_comp_id, label_entity_id))
            except (ValueError, IndexError):
                continue
    
    return observed_residues

def extract_entity_sequences(entity_poly_seq_data):
    """Extract sequences from entity_poly_seq data."""
    # CIF entity_poly_seq columns: entity_id, num, mon_id
    sequences = defaultdict(dict)
    
    for row in entity_poly_seq_data:
        if len(row) >= 3:
            try:
                entity_id = row[0]
                seq_num = int(row[1])
                residue_type = row[2]
                
                sequences[entity_id][seq_num] = residue_type
            except (ValueError, IndexError):
                continue
    
    return sequences

def create_residue_mapping(observed_residues, entity_sequences):
    """Create mapping between atom_site residue numbers and entity_poly_seq positions."""
    mapping = {}
    
    for auth_asym_id in observed_residues:
        chain_mapping = {}
        
        for label_seq_id, auth_seq_id, observed_type, entity_id in observed_residues[auth_asym_id]:
            # Map observed residue to sequence position
            if entity_id in entity_sequences and label_seq_id in entity_sequences[entity_id]:
                seq_residue_type = entity_sequences[entity_id][label_seq_id]
                
                chain_mapping[auth_seq_id] = {
                    'entity_poly_seq_pos': label_seq_id,
                    'entity_poly_seq_type': seq_residue_type,
                    'observed_type': observed_type,
                    'match': seq_residue_type == observed_type
                }
        
        if chain_mapping:  # Only add if we found valid mappings
            mapping[auth_asym_id] = chain_mapping
    
    return mapping

def print_mapping(mapping):
    """Print the residue mapping in a readable format."""
    for auth_asym_id, chain_mapping in mapping.items():
        print(f"\nChain {auth_asym_id}:")
        print("Auth_Seq_ID -> Entity_Poly_Seq_Pos (Entity_Type | Observed_Type) [Match]")
        print("-" * 70)
        
        for auth_seq_id in sorted(chain_mapping.keys(), key=lambda x: int(x) if x.isdigit() else float('inf')):
            mapping_info = chain_mapping[auth_seq_id]
            match_symbol = "✓" if mapping_info['match'] else "✗"
            
            print(f"{auth_seq_id:>10} -> {mapping_info['entity_poly_seq_pos']:>3} "
                  f"({mapping_info['entity_poly_seq_type']:>3} | {mapping_info['observed_type']:>3}) [{match_symbol}]")

def get_cif_residue_mapping(cif_file_path, verbose=False):
    """
    Wrapper function to get detailed residue mapping from a CIF file.
    
    Args:
        cif_file_path (str): Path to the CIF file
        verbose (bool): Whether to print detailed information
        
    Returns:
        dict: Mapping dictionary with auth_asym_id (chain IDs) as keys and residue mappings as values
              Format: {auth_asym_id: {auth_seq_id: {'entity_poly_seq_pos': int, 
                                                     'entity_poly_seq_type': str,
                                                     'observed_type': str,
                                                     'match': bool}}}
    
    Raises:
        FileNotFoundError: If the CIF file doesn't exist
        Exception: If there's an error processing the file
    """
    try:
        # Parse CIF file
        if verbose:
            print(f"Parsing CIF file: {cif_file_path}")
        
        atom_site_data, entity_poly_seq_data = parse_cif_file(cif_file_path)
        
        if verbose:
            print(f"Found {len(atom_site_data)} atom_site records")
            print(f"Found {len(entity_poly_seq_data)} entity_poly_seq records")
        
        # Extract data
        observed_residues = extract_observed_residues(atom_site_data)
        entity_sequences = extract_entity_sequences(entity_poly_seq_data)
        
        # Create mapping
        mapping = create_residue_mapping(observed_residues, entity_sequences)
        
        if verbose:
            print_mapping(mapping)
            total_observed = sum(len(entity_map) for entity_map in mapping.values())
            print(f"\nSummary: {total_observed} observed residues mapped across {len(mapping)} entities")
        
        return mapping
        
    except FileNotFoundError:
        raise FileNotFoundError(f"CIF file '{cif_file_path}' not found")
    except Exception as e:
        raise Exception(f"Error processing CIF file: {e}")

def get_resnum_mapping_fasta2cif(cif_file_path):
    """
    Wrapper function to get residue mapping from a CIF file, inverted from cif2fasta.
    entity_poly_seq_pos -> auth_seq_id
    
    Args:
        cif_file_path (str): Path to the CIF file
        
    Returns:
        dict: Simple mapping dictionary with auth_asym_id (chain IDs) as keys and 
              entity_poly_seq_pos -> auth_seq_id mappings as values
              Format: {auth_asym_id: {entity_poly_seq_pos: auth_seq_id}}
    
    Raises:
        FileNotFoundError: If the CIF file doesn't exist
        Exception: If there's an error processing the file
    """
    atom_site_data, entity_poly_seq_data = parse_cif_file(cif_file_path)
    
    # Extract data
    observed_residues = extract_observed_residues(atom_site_data)
    entity_sequences = extract_entity_sequences(entity_poly_seq_data)
        
    # Create inverted mapping
    inverted_mapping = {}
    
    for auth_asym_id in observed_residues:
        chain_mapping = {}
        
        for label_seq_id, auth_seq_id, observed_type, entity_id in observed_residues[auth_asym_id]:
            # Map entity_poly_seq_pos to auth_seq_id
            if entity_id in entity_sequences and label_seq_id in entity_sequences[entity_id]:
                chain_mapping[label_seq_id] = auth_seq_id
        
        if chain_mapping:  # Only add if we found valid mappings
            inverted_mapping[auth_asym_id] = chain_mapping
        
    return inverted_mapping

def get_resnum_mapping_cif2fasta(cif_file_path):
    """
    Wrapper function to get residue mapping from a CIF file.
    Maps author residue numbers (auth_seq_id) to FASTA sequence positions (entity_poly_seq_pos).
    
    Args:
        cif_file_path (str): Path to the CIF file
        
    Returns:
        dict: Simple mapping dictionary with auth_asym_id (chain IDs) as keys and 
              auth_seq_id -> entity_poly_seq_pos mappings as values
              Format: {auth_asym_id: {auth_seq_id: entity_poly_seq_pos}}
    
    Raises:
        FileNotFoundError: If the CIF file doesn't exist
        Exception: If there's an error processing the file
    """
    atom_site_data, entity_poly_seq_data = parse_cif_file(cif_file_path)
    
    # Extract data
    observed_residues = extract_observed_residues(atom_site_data)
    entity_sequences = extract_entity_sequences(entity_poly_seq_data)
        
    # Create mapping: auth_seq_id -> entity_poly_seq_pos (same as fasta2cif)
    mapping = {}
    
    for auth_asym_id in observed_residues:
        chain_mapping = {}
        
        for label_seq_id, auth_seq_id, observed_type, entity_id in observed_residues[auth_asym_id]:
            # Map auth_seq_id to entity_poly_seq_pos (FASTA position)
            if entity_id in entity_sequences and label_seq_id in entity_sequences[entity_id]:
                chain_mapping[auth_seq_id] = label_seq_id
        
        if chain_mapping:  # Only add if we found valid mappings
            mapping[auth_asym_id] = chain_mapping
        
    return mapping

def main():
    if len(sys.argv) != 2:
        print("Usage: python pdb_residue_mapper.py <pdb_file.cif>")
        sys.exit(1)
    
    cif_file = sys.argv[1]
    
    try:
        # Use the wrapper function
        mapping = get_resnum_mapping_fasta2cif(cif_file)
        mapping_detailed = get_cif_residue_mapping(cif_file, True)
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    print(mapping)
    print(mapping_detailed)
    
def test_get_resnum_mapping(cif_file, pdb_id):
    """
    Test the get_resnum_mapping function by downloading FASTA sequences and comparing
    mapped residues to verify they match.
    
    Args:
        cif_file (str): Path to the CIF file
        pdb_id (str): PDB ID for downloading FASTA sequences
        
    Returns:
        bool: True if all mapped residues match, False otherwise
    """
    try:
        # Import the download function
        from functions_download import download_pdb_sequence
        
        # Get the residue mapping from CIF file
        print(f"Getting residue mapping from CIF file: {cif_file}")
        mapping = get_resnum_mapping_cif2fasta(cif_file)
        
        # Download FASTA sequences
        print(f"Downloading FASTA sequences for PDB ID: {pdb_id}")
        fasta_sequences = download_pdb_sequence(pdb_id)
        
        # Create a mapping from chain_id to sequence
        fasta_mapping = {}
        for chain in fasta_sequences:
            fasta_mapping[chain['chain_id']] = chain['sequence']
        
        print(f"CIF mapping: {mapping}")
        print(f"FASTA sequences: {fasta_sequences}")
        
        # Get detailed residue mapping to access residue types
        detailed_mapping = get_cif_residue_mapping(cif_file)
        
        all_matches = True
        total_residues_tested = 0
        
        for chain_id, chain_mapping in mapping.items():
            if chain_id not in fasta_mapping:
                print(f"Warning: Chain {chain_id} not found in FASTA sequences")
                continue
                
            fasta_sequence = fasta_mapping[chain_id]
            print(f"\nTesting chain {chain_id} (FASTA length: {len(fasta_sequence)}):")
            
            for auth_seq_id, entity_poly_seq_pos in chain_mapping.items():
                # Convert to 0-based index for FASTA sequence
                fasta_index = entity_poly_seq_pos - 1
                
                if fasta_index < 0 or fasta_index >= len(fasta_sequence):
                    print(f"  Auth {auth_seq_id} -> Seq {entity_poly_seq_pos}: Index out of range for FASTA")
                    all_matches = False
                    continue
                
                fasta_residue = fasta_sequence[fasta_index]
                
                # Get the observed residue type from detailed mapping
                if chain_id in detailed_mapping and auth_seq_id in detailed_mapping[chain_id]:
                    observed_residue = detailed_mapping[chain_id][auth_seq_id]['observed_type']
                    
                    # Convert 3-letter to 1-letter amino acid codes for comparison
                    aa_conversion = {
                        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
                    }
                    
                    observed_1letter = aa_conversion.get(observed_residue, observed_residue)
                    
                    match = observed_1letter == fasta_residue
                    match_symbol = "✓" if match else "✗"
                    
                    print(f"  Auth {auth_seq_id} -> Seq {entity_poly_seq_pos}: {observed_residue}({observed_1letter}) vs FASTA {fasta_residue} [{match_symbol}]")
                    
                    if not match:
                        all_matches = False
                    
                    total_residues_tested += 1
                else:
                    print(f"  Auth {auth_seq_id} -> Seq {entity_poly_seq_pos}: No detailed mapping found")
                    all_matches = False
        
        print(f"\nTest Summary:")
        print(f"Total residues tested: {total_residues_tested}")
        print(f"All residues match: {all_matches}")
        
        return all_matches
        
    except ImportError:
        print("Error: Could not import download_pdb_sequence from functions_download.py")
        return False
    except Exception as e:
        print(f"Error in test_get_resnum_mapping: {e}")
        return False
    
if __name__ == "__main__":
    main()
    # test_get_resnum_mapping('/media/markus/1AC1D5876B41B20C/MPI/data/PDB/1c5w.cif', '1c5w')

