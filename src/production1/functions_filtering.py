import pandas as pd
from typing import List, Tuple, Any, Dict
import os, hashlib
from pathlib import Path

def create_all_pairs(arm_df: pd.DataFrame, tf_df: pd.DataFrame) -> pd.DataFrame:
    """Create a dataframe with all possible pairs from the cartesian product of the two dataframes arm_df and tf_df.

    This function performs a full cartesian join between the armadillo proteins dataframe and the transcription factor proteins
    dataframe, generating all possible combinations between them.

    Args:
        arm_df (pd.DataFrame): DataFrame containing armadillo proteins
        tf_df (pd.DataFrame): DataFrame containing transcription factor proteins

    Returns:
        pd.DataFrame: DataFrame containing all possible pairs between armadillo and transcription factor proteins
        with a unique pair_id column for each combination
    """
    # Create a key for cross join
    arm_df_temp = arm_df.copy()
    tf_df_temp = tf_df.copy()

    arm_df_temp['key'] = 1
    tf_df_temp['key'] = 1

    # Perform a cross join using the dummy key
    pairs_df = pd.merge(arm_df_temp, tf_df_temp, on='key', suffixes=('_arm', '_tf'))

    # Drop the dummy key column
    pairs_df = pairs_df.drop('key', axis=1)

    # Create pair_id column for consistency with other functions in the pipeline
    pairs_df['pair_id'] = pairs_df.apply(lambda row: str(tuple(sorted([row['Entry_arm'].upper(), row['Entry_tf'].upper()]))), axis=1)

    print(f"Created {len(pairs_df)} possible protein pairs between {len(arm_df)} armadillo proteins and {len(tf_df)} transcription factors")

    return pairs_df

def find_subranges(data: List[float], threshold: float, min_length: int) -> List[Tuple[int, int]]:
    """Find continuous subranges in data where values exceed the threshold for at least min_length positions.

    Args:
        data (List[float]): List of numerical values to analyze
        threshold (float): Minimum value to be considered part of a subrange
        min_length (int): Minimum length a subrange must have to be included in results

    Returns:
        List[Tuple[int, int]]: List of tuples containing start and end indices of subranges
    """
    subranges = []
    start = None

    for i, value in enumerate(data):
        if value >= threshold:
            if start is None:
                start = i
        else:
            # End of a potential subrange
            if start is not None:
                if i - start >= min_length:
                    subranges.append((start, i - 1))
                start = None

    # Check if we ended with an ongoing subrange
    if start is not None:
        if len(data) - start >= min_length:
            subranges.append((start, len(data) - 1))

    return subranges

assert find_subranges([1,1,2,4,5,6,2,3,1], 2, 3) == [(2,7)]
assert find_subranges([1,1,2,4,5,6,2,3,1], 2, 20) == []
assert find_subranges([], 2, 20) == []
assert find_subranges([1,2,2,2,3,3,3,3,2,3,3,2,1,3,3,3], 3, 3) == [(4,7), (13,15)]

def contains_any_annotation(cell_value: Any, annotations_list: List[str]) -> bool:
    """Check if any of the annotations are in the column value.

    Args:
        cell_value (Any): Cell value from DataFrame column to check
        annotations_list (List[str]): List of annotation strings to check for

    Returns:
        bool: True if any annotation is found in the cell_value, False otherwise
    """
    if pd.isna(cell_value):
        return False
    for annotation in annotations_list:
        if annotation in cell_value:
            return True
    return False

def print_to_fasta(id: str, seq: str, path: str, comment: str = '') -> None:
    """Create a FASTA file with the given ID and sequence.

    This function creates a new FASTA file with the specified ID as the filename
    (with .fasta extension) and writes the sequence in standard FASTA format.

    Args:
        id (str): Protein ID to use as both the filename and the FASTA header
        seq (str): Amino acid sequence to write to the file
        path (str): Directory path where the FASTA file should be created
        comment (str, optional): Optional comment to add to the FASTA header. Defaults to ''.

    Returns:
        None

    Raises:
        OSError: If the directory cannot be created or file cannot be written
    """
    import os
    
    filename = f"{id}.fasta"
    full_path = os.path.join(path, filename)
    
    os.makedirs(path, exist_ok=True)
    
    header = f">{id}"
    if comment:
        header += f'|{comment}\n'
    else:
        header += '\n'
    
    with open(full_path, 'w') as f:
        f.write(header)
        f.write(f"{seq}\n")        

def add_iupred3(df: pd.DataFrame, type: str, smoothing, cache_dir, threshold, min_length_region, iupred_path, debug=False):
    import sys
    sys.path.append(iupred_path)
    import iupred3_lib
    
    os.makedirs(cache_dir, exist_ok=True)
    df['iupred3'] = None
    
    i = 0
    for ind, row in df.iterrows():
        if i % int(len(df)/20) == 0:
            if debug:
                print(f"Processed {i} of {len(df)} sequences.")
        i += 1
        seq = str(row['Sequence'])
        if seq == '':
            df.at[ind, 'iupred3'] = ''
            df.at[ind, 'num_disordered_regions'] = 0
            continue
        cache_key = hashlib.md5((str(seq) + type + str(smoothing)).encode()).hexdigest()
        cache_file = os.path.join(cache_dir, f"{cache_key}.txt")
        if os.path.exists(cache_file):
            with open(cache_file, 'r') as f:
                iupred3_str = f.read().strip()
                iupred3 = [float(x) for x in iupred3_str.split(',')]
        else:
            iupred3 = iupred3_lib.iupred(seq, type, smoothing=smoothing)[0]
            iupred3_str = ','.join(map(str, iupred3))
            with open(cache_file, 'w') as f:
                f.write(iupred3_str)
        num_disordered_regions = len(find_subranges(iupred3, threshold, min_length_region))
        df.at[ind, 'iupred3'] = iupred3_str
        df.at[ind, 'num_disordered_regions'] = num_disordered_regions
    return df

def get_chain_id(l: int) -> str:
    """Generate chain ID from index.

    Args:
        l (int): Index to convert to chain ID

    Returns:
        str: Chain ID string
    """
    id = ''
    while l > 25:
        id += 'Z'
        l -= 26
    id += chr(ord('A') + l)
    return id

def check_interface(pdb_id, chain_X, chain_Y, data_dir, min_atoms=10, max_distance=5) -> bool:
    """check if in the specified pdb there is an interface between chain_X and chain_Y
    using Gregors data

    Args:
        pdb_id (_type_): _description_
        chain_A (_type_): _description_
        chain_B (_type_): _description_

    Returns:
        bool: _description_
    """
    # find files matching pattern and check for atom pairs within 5A between the two chains
    files = list(Path(data_dir).rglob(f"{str.upper(pdb_id)}_detailed_interactions.csv"))

    if len(files) > 1:
        print(f"ERROR: for {pdb_id}, multiple files were found.")
        return False
    elif len(files) == 0:
        print(f"ERROR: for {pdb_id}, no files were found.")
        return False

    df = pd.read_csv(files[0], sep=',')

    atom_counter = 0
    for _,row in df.iterrows():
        if row['PDB_ID'] == pdb_id.upper() and ((row['Chain_A'] == chain_X and row['Chain_B'] == chain_Y) or (row['Chain_A'] == chain_Y and row['Chain_B'] == chain_X)):
            if row['Distance'] <= max_distance:
                atom_counter += 1

    return atom_counter >= min_atoms