import pandas as pd
from typing import List, Tuple, Any
import os, hashlib, itertools
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

def add_iupred3(df: pd.DataFrame, type: str, smoothing, cache_dir, threshold, min_length_region, iupred_path):
    import sys
    sys.path.append(iupred_path)
    import iupred3_lib
    
    os.makedirs(cache_dir, exist_ok=True)
    df['iupred3'] = None
    
    for ind, row in df.iterrows():
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



def get_job_name(id_tf: str, id_arm: str, df: pd.DataFrame):
    row_arm = df[df['Entry'] == id_arm]
    if row_arm.empty:
        raise Exception(f"Not found in df: {id_arm}")
    row_tf = df[df['Entry'] == id_tf]
    if row_tf.empty:
        raise Exception(f"Not found in df: {id_tf}")
    length_arm = row_arm['Length'].iloc[0]
    length_tf = row_tf['Length'].iloc[0]
    return str.lower(f"{id_arm}_1-{length_arm}_{id_tf}_1-{length_tf}")

def get_all_chain_mappings(native_chains, model_chains) -> List:
    all_mappings = []
    
    if len(native_chains) > len(model_chains):
        all_subsets = itertools.combinations(native_chains, len(model_chains))
        for subset in all_subsets:
            all_mappings.extend(all_bijective_mappings(subset, model_chains))
    elif len(model_chains) > len(native_chains):
        all_subsets = itertools.combinations(model_chains, len(native_chains))
        for subset in all_subsets:
            all_mappings.extend(all_bijective_mappings(native_chains, subset))
    else:
        all_mappings = all_bijective_mappings(native_chains, model_chains)
    return all_mappings

def all_bijective_mappings(A, B):
    """
    Return a list containing every dictionary that maps each element of A
    to a unique element of B.  A and B must be the same length.
    """

    # For each permutation of B, zip it with A to make a mapping dict
    return [dict(zip(A, perm)) for perm in itertools.permutations(B)]

def get_file_path(filename: str, search_dir: str):
    """search the file with filename in dir and return the full path if it exists, otherwise return False

    Args:
        filename (str): _description_
    """
    path = Path(search_dir)
    for file in path.rglob(filename):
        return file  # return the first match
    return False

def calculate_dockq_scores(df: pd.DataFrame, native_path_prefix: str, model_results_dir: str, 
                          all_uniprot: pd.DataFrame) -> pd.DataFrame:
    """Calculate DockQ scores for protein structures and append results to dataframe.
    
    This function calculates DockQ scores for all structures in the input dataframe by comparing
    predicted models with native PDB structures using all possible chain mappings.
    
    Args:
        df (pd.DataFrame): DataFrame containing structure information with columns 'query_tf', 'query_arm', 'pdb_id'
        native_path_prefix (str): Path prefix for native PDB structure files
        model_results_dir (str): Directory containing predicted model files
        all_uniprot (pd.DataFrame): UniProt dataframe for getting job names
        
    Returns:
        pd.DataFrame: Input dataframe with added 'chain_map' and 'dockq_score' columns
    """
    from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
    
    # Create copies to avoid modifying the original dataframe
    result_df = df.copy()
    result_df['chain_map'] = None
    result_df['dockq_score'] = None
    
    no_model = []
    no_native = []
    
    # Process each structure
    for index, row in result_df.iterrows():
        job_name = get_job_name(row['query_tf'].split('|')[0], row['query_arm'].split('|')[0], all_uniprot)
        model_path = get_file_path(f'{job_name}_model.cif', model_results_dir)
        native_path_cif = f'{native_path_prefix}{row["pdb_id"].lower()}.cif'
        
        # Check if both paths exist before loading
        if not model_path:
            no_model.append((model_path, job_name))
            continue

        if not os.path.exists(native_path_cif):
            no_native.append((native_path_cif, job_name))
            continue
        
        model = load_PDB(model_path)
        native = load_PDB(native_path_cif)
        native_chains = [chain.id for chain in model]
        model_chains = [chain.id for chain in native]
        chain_map_dict = {}
        
        if len(model_chains) > 2:
            print(f"{row['pdb_id']}: Warning: Native structure ({native}) has more than two chains!")
            
        for chain_map in get_all_chain_mappings(model_chains, native_chains):
            try:
                dockQ = run_on_all_native_interfaces(model, native, chain_map=chain_map)[1]
            except Exception as e:
                print(f"Exception for {row['pdb_id']}: {e}. Comment in review: {row['comment']}")
                break
            chain_map_dict[str(chain_map)] = dockQ
            
        if chain_map_dict:
            best_chain_map = max(chain_map_dict.keys(), key=(lambda key: chain_map_dict[key]))
            best_dockq_score = chain_map_dict[best_chain_map]
            
            # Store results in the DataFrame
            result_df.at[index, 'chain_map'] = best_chain_map
            result_df.at[index, 'dockq_score'] = best_dockq_score

    # Report missing files
    if no_model:
        print("Missing model files:")
        print(str([(t[0], t[1]) for t in no_model]))
    if no_native:
        print("Missing native files:")
        for path in no_native:
            print(f"  NATIVE: {path}")
    
    return result_df