import requests, os, json
import os
from pathlib import Path
from typing import List, Dict
import pandas as pd

def download_pdb_structure(pdb_id, output_dir="/home/markus/MPI_local/data/PDB", file_format="cif"):
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

        print(f"Successfully downloaded {filename} to {output_path}")
        return output_path

    except requests.exceptions.RequestException as e:
        print(f"Error downloading {pdb_id}: {e}")
        return None
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None

def download_pdb_structures(pdb_ids: set, output_dir="/home/markus/MPI_local/data/PDB", file_format="cif"):
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

    print(f"Successfully processed: {downloaded_count}")
    print(f"Failed downloads: {failed_count}")
    print(f"Total processed: {len([pdb_id for pdb_id in pdb_ids if not pd.isna(pdb_id)])}")

def download_pdb_sequence(pdb_id) -> List[Dict[str, str]]:
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

                # Extract chain ID(s) from header
                # Header formats:
                # >4HHB_1|Chain A|Hemoglobin subunit alpha|Homo sapiens (9606)
                # >8H36_1|Chains A[auth D], H[auth E]|E3 ubiquitin-protein ligase RBX1|Homo sapiens (9606)
                header_parts = line.split('|')
                if len(header_parts) >= 2:
                    chain_part = header_parts[1].strip()

                    # Handle multiple chains with auth IDs
                    if chain_part.startswith('Chains '):
                        # Extract auth chain IDs from format like "Chains A[auth D], H[auth E]"
                        chains_str = chain_part.replace('Chains ', '')
                        chain_entries = [entry.strip() for entry in chains_str.split(',')]
                        current_chains = []

                        for entry in chain_entries:
                            if '[auth ' in entry and ']' in entry:
                                # Extract auth chain ID
                                auth_start = entry.find('[auth ') + 6
                                auth_end = entry.find(']', auth_start)
                                auth_chain = entry[auth_start:auth_end]
                                current_chains.append(auth_chain)
                            else:
                                # Fallback to regular chain ID if no auth
                                current_chains.append(entry.strip())

                        current_chain = current_chains
                    elif chain_part.startswith('Chain '):
                        # Single chain format
                        current_chain = [chain_part.replace('Chain ', '')]
                    else:
                        raise Exception("No chain ID!")
                else:
                    raise Exception("No chain ID!")

                current_sequence = ""
            else:
                # Accumulate sequence lines
                current_sequence += line.strip()

        # Add final sequence(s)
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

        if sequences:
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