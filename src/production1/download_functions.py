import requests, os, json, re
from pathlib import Path
from typing import List, Dict
import pandas as pd
from collections import defaultdict

amino_acids = [
    "A",  # Alanine
    "R",  # Arginine
    "N",  # Asparagine
    "D",  # Aspartic acid
    "C",  # Cysteine
    "E",  # Glutamic acid
    "Q",  # Glutamine
    "G",  # Glycine
    "H",  # Histidine
    "I",  # Isoleucine
    "L",  # Leucine
    "K",  # Lysine
    "M",  # Methionine
    "F",  # Phenylalanine
    "P",  # Proline
    "S",  # Serine
    "T",  # Threonine
    "W",  # Tryptophan
    "Y",  # Tyrosine
    "V",  # Valine
]


def download_pdb_structure(pdb_id, output_dir="/home/markus/MPI_local/data/PDB", file_format="cif", debug=False):
    """
    Download a PDB structure from the online PDB database.

    Automatically checks if the file already exists and skips download if present.

    Parameters:
    -----------
    pdb_id : str
        The 4-character PDB ID (e.g., "3ouw")
    output_dir : str, optional
        Directory to save the downloaded structure (default: "/home/markus/MPI_local/data/PDB")
    file_format : str, optional
        File format to download ("cif" or "pdb", default: "cif")

    Returns:
    --------
    str
        Path to the file (either existing or newly downloaded), or None if download failed
    """
    # Ensure PDB ID is lowercase for URL
    pdb_id = pdb_id.lower()

    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Define URL based on format
    if file_format.lower() == "cif":
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        filename = f"{pdb_id}.cif"
    elif file_format.lower() == "pdb":
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        filename = f"{pdb_id}.pdb"
    else:
        print(f"Unsupported format: {file_format}. Use 'cif' or 'pdb'.")
        return None

    output_path = os.path.join(output_dir, filename)

    # Check if file already exists
    if os.path.exists(output_path):
        return output_path

    try:
        # Download the file
        response = requests.get(url, timeout=30)
        response.raise_for_status()  # Raises an HTTPError for bad responses

        # Save the file
        with open(output_path, 'w') as f:
            f.write(response.text)

        if debug:
            print(f"Successfully downloaded {filename} to {output_path}")
        return output_path

    except requests.exceptions.RequestException as e:
        print(f"Error downloading {pdb_id}: {e}")
        return None
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None

def download_pdb_structures(pdb_ids: set, output_dir="/home/markus/MPI_local/data/PDB", file_format="cif", debug=True):
    downloaded_count = 0
    failed_count = 0

    for pdb_id in pdb_ids:
        if downloaded_count % 10 == 0:
            print(f"downloaded {downloaded_count} of {len(pdb_ids)}.")
        if pd.isna(pdb_id):  # Skip NaN values
            continue

        result = download_pdb_structure(pdb_id, output_dir, file_format)
        if result:
            downloaded_count += 1
        else:
            failed_count += 1

    if debug:
        print(f"Successfully processed: {downloaded_count}")
        print(f"Failed downloads: {failed_count}")
        print(f"Total processed: {len([pdb_id for pdb_id in pdb_ids if not pd.isna(pdb_id)])}")

def download_pdb_sequence(pdb_id, debug=False) -> List[Dict[str, str]]:
    """
    Download the amino acid sequences for each chain in a PDB structure.

    Parameters:
    -----------
    pdb_id : str
        The 4-character PDB ID (e.g., "3ouw")

    Returns:
    --------
    list
        List of dictionaries with 'chain_id' and 'sequence' keys, or None if download failed
        Example: [{'chain_id': 'A', 'sequence': 'MKTI...'}, {'chain_id': 'B', 'sequence': 'AFGL...'}]
    """
    # Ensure PDB ID is lowercase for URL
    pdb_id = pdb_id.lower()

    # RCSB PDB REST API URL for FASTA sequences
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"

    try:
        # Download the FASTA file
        response = requests.get(url, timeout=30)
        response.raise_for_status()

        # Parse FASTA content
        fasta_content = response.text.strip()

        if not fasta_content:
            print(f"No sequence data found for PDB ID: {pdb_id}")
            return []

        sequences = []
        current_chain = None
        current_sequence = ""
        for line in fasta_content.split('\n'):
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_chain is not None:
                    if isinstance(current_chain, list):
                        # Multiple chains - add each chain separately with the same sequence
                        for chain_id in current_chain:
                            sequences.append({
                                'chain_id': chain_id,
                                'sequence': current_sequence
                            })
                    else:
                        # Single chain
                        sequences.append({
                            'chain_id': current_chain,
                            'sequence': current_sequence
                        })

               
                current_chain = extract_chain_ID(line)
                current_sequence = ""
            else:
                # Accumulate sequence lines
                seq = line.strip()
                if not all(c in amino_acids for c in seq):
                    raise Exception(f"sequence contains non-standard amino acids or nucleotides: {seq}, ID: {pdb_id}")
                current_sequence += line.strip()

        # Add final sequence(s)
        if current_chain is not None:
            if isinstance(current_chain, list):
                # Multiple chains - add each chain separately with the same sequence
                deduplicate(current_chain)
                for chain_id in current_chain:
                    if not bool(re.fullmatch(r'[A-Z0-9]+', chain_id)):
                        raise Exception(f"ID must only contain alphanumeric upper-case chars: {chain_id}")
                    sequences.append({
                        'chain_id': chain_id,
                        'sequence': current_sequence
                    })
            else:
                # Single chain
                if not bool(re.fullmatch(r'[A-Z0-9]+', current_chain)):
                    raise Exception(f"ID must only contain alphanumeric upper-case chars: {current_chain}")
                sequences.append({
                    'chain_id': current_chain,
                    'sequence': current_sequence
                })

        if sequences:
            if debug:
                print(f"Successfully downloaded sequences for {pdb_id}: {len(sequences)} chain(s)")
            return sequences
        else:
            print(f"No sequences found in FASTA data for {pdb_id}")
            return []

    except requests.exceptions.RequestException as e:
        print(f"Error downloading sequences for {pdb_id}: {e}")
        return []
    except Exception as e:
        print(f"Error while processing sequences for {pdb_id}: {e}")
        return []
    
def extract_chain_ID(line: str) -> str|List[str]:
    """Extract chain ID(s) from header
    Header formats:
    >4HHB_1|Chain A|Hemoglobin subunit alpha|Homo sapiens (9606)
    >8H36_1|Chains A[auth D], H[auth E]|E3 ubiquitin-protein ligase RBX1|Homo sapiens (9606)
    """               
    header_parts = line.split('|')
    if len(header_parts) >= 2:
        chain_part = header_parts[1].strip()

        # Handle multiple chains with auth IDs
        if chain_part.startswith('Chains '):
            # Extract auth chain IDs from format like "Chains A[auth D], H[auth E]"
            chains_str = chain_part.replace('Chains ', '')
            current_chain = [format_chain_string(entry) for entry in chains_str.split(',')]
        elif chain_part.startswith('Chain '):
            # Single chain format
            chain_part = chain_part.replace('Chain ', '')
            current_chain = [format_chain_string(chain_part)]
        else:
            raise Exception("No chain ID!")
    else:
        raise Exception("No chain ID!")
    return current_chain
    
def format_chain_string(chain_str: str) -> str:
    chain_str = chain_str.replace('[', '')
    chain_str = chain_str.replace(']', '')
    chain_str = chain_str.upper()
    chain_str = chain_str.replace(' ', '')
    return chain_str

def int_to_base36(n: int) -> str:
    """Convert int to base36 string using [a-z0-9]."""
    chars = "abcdefghijklmnopqrstuvwxyz0123456789"
    res = ""
    while True:
        n, r = divmod(n, 36)
        res = chars[r] + res
        if n == 0:
            return res

def deduplicate(strings: list[str]) -> list[str]:
    """
    Remove duplicates from a list of strings by appending base36 suffixes to duplicates.

    This function processes a list of strings and ensures all elements in the result
    are unique. When duplicates are encountered, they are made unique by appending
    a base36-encoded suffix (starting from 'a' for the first duplicate, 'b' for the
    second, etc.).

    Parameters:
    -----------
    strings : list[str]
        List of strings that may contain duplicates

    Returns:
    --------
    list[str]
        List of unique strings with duplicates renamed using base36 suffixes.
        The order of first occurrences is preserved.
        Example: ['A', 'B', 'Aa', 'C', 'Ab'] for input ['A', 'B', 'A', 'C', 'A']
    """
    seen = defaultdict(int)
    used = set(strings)
    result = []

    for s in strings:
        if seen[s] == 0 and s not in result:
            result.append(s)
        else:
            i = seen[s]
            new = f"{s}{int_to_base36(i)}"
            while new in used:
                i += 1
                new = f"{s}{int_to_base36(i)}"
            result.append(new)
            used.add(new)
            seen[s] = i
        seen[s] += 1
    return result


def get_pdb_chains_to_uniprot(pdb_id: str, cache_dir: str = ".pdb_cache") -> Dict[str, str]:
    """
    Return a mapping of PDB chain IDs to UniProt accessions for a given PDB entry.
    Results are cached permanently as JSON files in `cache_dir`.

    Parameters
    ----------
    pdb_id : str
        Four-letter PDB identifier (e.g. "1ABC").
    cache_dir : str, optional
        Directory where results are cached (default: ".pdb_cache").

    Returns
    -------
    Dict[str, str]
        Dictionary mapping chain IDs (e.g. 'A') to UniProt IDs (e.g. 'P12345').
    """
    os.makedirs(cache_dir, exist_ok=True)
    cache_file = os.path.join(cache_dir, f"{pdb_id.lower()}.json")

    # Load from cache if exists
    if os.path.exists(cache_file):
        with open(cache_file, "r") as f:
            return json.load(f)

    # Otherwise fetch from API
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/"
    mapping: Dict[str, str] = {}
    i = 1
    while True:
        r = requests.get(url + str(i))
        if r.status_code != 200:
            break
        d = r.json()
        chains = d["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"]
        refs = d["rcsb_polymer_entity_container_identifiers"].get("reference_sequence_identifiers", [])
        for ref in refs:
            if ref.get("database_name") == "UniProt":
                for c in chains:
                    mapping[c] = ref["database_accession"]
        i += 1

    # Save to cache
    with open(cache_file, "w") as f:
        json.dump(mapping, f)

    return mapping