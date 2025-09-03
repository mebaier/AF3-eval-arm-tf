import requests, os, json, re
from pathlib import Path
from typing import List, Dict
import pandas as pd

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


def download_pdb_structure(pdb_id, output_dir="/home/markus/MPI_local/data/PDB", file_format="cif", cache_dir=None, debug=False):
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

    # Determine cache path - use cache_dir if provided, otherwise use output_dir
    cache_location = cache_dir if cache_dir is not None else output_dir
    cache_path = os.path.join(cache_location, filename)

    # Check if file already exists in cache
    if os.path.exists(cache_path):
        if not os.path.exists(output_path):
            # Copy file from cache to output path
            with open(cache_path, 'r') as src, open(output_path, 'w') as dst:
                dst.write(src.read())
        return cache_path

    try:
        # Download the file
        response = requests.get(url, timeout=30)
        response.raise_for_status()  # Raises an HTTPError for bad responses

        # Save the file to output dir and cache dir
        with open(output_path, 'w') as f:
            f.write(response.text)
            f.close()

        if not cache_path == output_path:
            with open(cache_path, 'w') as f:
                f.write(response.text)
                f.close()

        if debug:
            print(f"Successfully downloaded {filename} to {output_path}")
        return output_path

    except requests.exceptions.RequestException as e:
        print(f"Error downloading {pdb_id}: {e}")
        return None
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None

def download_pdb_structures(pdb_ids: set, output_dir="/home/markus/MPI_local/data/PDB", file_format="cif", cache_dir=None, debug=False):
    downloaded_count = 0
    failed_count = 0

    for pdb_id in pdb_ids:
        if downloaded_count % 10 == 0:
            if debug:
                print(f"downloaded {downloaded_count} of {len(pdb_ids)}.")
        if pd.isna(pdb_id):  # Skip NaN values
            continue

        cache_location = cache_dir if cache_dir is not None else output_dir
        result = download_pdb_structure(pdb_id, output_dir, file_format, cache_location, debug)
        if result:
            downloaded_count += 1
        else:
            failed_count += 1

    if debug:
        print(f"Successfully processed: {downloaded_count}")
        print(f"Failed downloads: {failed_count}")
        print(f"Total processed: {len([pdb_id for pdb_id in pdb_ids if not pd.isna(pdb_id)])}")

def add_sequences_to_list(chain, sequence, sequences_list):
    """Helper function to add sequences to the list, handling both single chains and chain lists"""
    if chain is None:
        return

    # Ensure chain is always a list for uniform processing
    chain_list = chain if isinstance(chain, list) else [chain]

    for chain_id in chain_list:
        sequences_list.append({
            'chain_id': chain_id,
            'sequence': sequence
        })

def download_pdb_sequence(pdb_id, cache_dir=None, debug=False) -> List[Dict[str, str]]:
    """
    Download the amino acid sequences for each chain in a PDB structure.
    Use as chain IDs the original authors mapping!

    Parameters:
    -----------
    pdb_id : str
        The 4-character PDB ID (e.g., "3ouw")
    cache_dir : str, optional
        Directory to cache downloaded FASTA files. If None, no caching is performed.
    debug : bool, optional
        Whether to print debug information (default: False)

    Returns:
    --------
    list
        List of dictionaries with 'chain_id' and 'sequence' keys, or None if download failed
        Example: [{'chain_id': 'A', 'sequence': 'MKTI...'}, {'chain_id': 'B', 'sequence': 'AFGL...'}]
    """
    # Ensure PDB ID is lowercase for URL
    pdb_id = pdb_id.lower()

    # Set up caching
    fasta_content = None
    cache_path = None
    if cache_dir is not None:
        Path(cache_dir).mkdir(parents=True, exist_ok=True)
        cache_path = os.path.join(cache_dir, f"{pdb_id}.fasta")
        
        # Check if cached file exists
        if os.path.exists(cache_path):
            try:
                with open(cache_path, 'r') as f:
                    fasta_content = f.read().strip()
                if debug:
                    print(f"Loaded FASTA from cache: {cache_path}")
            except Exception as e:
                if debug:
                    print(f"Error reading cached file {cache_path}: {e}")
                fasta_content = None

    # Download if not cached or cache failed
    if fasta_content is None:
        # RCSB PDB REST API URL for FASTA sequences
        url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"

        try:
            # Download the FASTA file
            response = requests.get(url, timeout=30)
            response.raise_for_status()

            # Parse FASTA content
            fasta_content = response.text.strip()
            
            # Save to cache if cache_dir is provided
            if cache_path is not None:
                try:
                    with open(cache_path, 'w') as f:
                        f.write(fasta_content)
                    if debug:
                        print(f"Saved FASTA to cache: {cache_path}")
                except Exception as e:
                    if debug:
                        print(f"Warning: Could not save to cache {cache_path}: {e}")

        except requests.exceptions.RequestException as e:
            print(f"Error downloading sequences for {pdb_id}: {e}")
            return []
        except Exception as e:
            print(f"Error while downloading sequences for {pdb_id}: {e}")
            return []

    # Parse FASTA content
    try:
        if not fasta_content:
            print(f"No sequence data found for PDB ID: {pdb_id}")
            return []

        sequences = []
        current_chain = None
        current_sequence = ""

        for line in fasta_content.split('\n'):
            if line.startswith('>'):
                # Save previous sequence if exists
                add_sequences_to_list(current_chain, current_sequence, sequences)
                
                current_chain = get_chain_ID_from_header(line)
                current_sequence = ""
            else:
                # Accumulate sequence lines
                seq = line.strip()
                if not all(c in amino_acids for c in seq):
                    raise Exception(f"sequence contains non-standard amino acids or nucleotides: {seq}, ID: {pdb_id}")
                current_sequence += line.strip()

        # Add final sequence(s)
        add_sequences_to_list(current_chain, current_sequence, sequences)

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

def get_chain_ID_from_header(line: str) -> str|List[str]:
    """Extract chain ID(s) from header line
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
            chains_str = extract_chain_ID(chain_part.replace('Chains ', ''))
            current_chain = chains_str.split(',')
        elif chain_part.startswith('Chain '):
            # Single chain format
            chain_part = extract_chain_ID(chain_part.replace('Chain ', ''))
            current_chain = [chain_part]
        else:
            raise Exception("No chain ID!")
    else:
        raise Exception("No chain ID!")
    return [c.strip() for c in current_chain]

def extract_chain_ID(id_str: str) -> str:
    if 'auth' in id_str:
        m = re.search(r"\[auth\s+([a-zA-Z0-9]+)\]", id_str)
        if not m:
            raise Exception(f"Error when extracting auth ID from: {id_str}")
        return m.group(1)
    return id_str

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