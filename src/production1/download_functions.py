import requests
import os
from pathlib import Path
from typing import List, Dict

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

# Example usage:
# pdb_file = download_pdb_structure("3ouw")
# pdb_file = download_pdb_structure("1abc", output_dir="/custom/path", file_format="pdb")