import pandas as pd
import os, json, itertools, pickle, hashlib
from typing import List, Dict, Any, Tuple
from pathlib import Path
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBExceptions import PDBIOException
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces, run_on_chains
from functions_download import *

def create_pair_id(row: pd.Series) -> str:
    """Create a pair ID from the job name by extracting the Uniprot IDs and sorting them.

    This function parses an AlphaFold job name to extract the protein IDs and creates a
    standardized pair identifier by sorting the IDs alphabetically.

    Args:
        row (pd.Series): DataFrame row containing 'job_name' column

    Returns:
        str: String representation of tuple of sorted protein IDs
    """
    parts = str(row['job_name']).split('_')
    return str(tuple(sorted([parts[0].upper(), parts[2].upper()])))

def find_summary_files(directories: List[str]) -> List[Dict[str, Any]]:
    """Recursively find summary_confidences.json files.

    This function searches through the given directories and all subdirectories to find
    AlphaFold summary confidence files and loads them into a list of dictionaries.
    If a job name is duplicated, don't load the job twice

    Args:
        directories (List[str]): List of root directories to search for summary files

    Returns:
        List[Dict[str, Any]]: List of dictionaries containing the contents of the summary files
    """
    all_jobs_data = []
    duplicates = 0

    for directory in directories:
        if not os.path.exists(directory):
            raise FileNotFoundError(f"Directory '{directory}' does not exist")
        if not os.path.isdir(directory):
            raise NotADirectoryError(f"Path '{directory}' is not a directory")

        for root, dirs, files in os.walk(directory):
            # Check for summary_confidences.json files
            # the '_' prefix is important to distinguish the best overall file form the sample files
            summary_files = [f for f in files if f.endswith('_summary_confidences.json')]

            if summary_files:
                for summary_file in summary_files:
                    # Extract job name from the filename
                    job_name = summary_file.replace('_summary_confidences.json', '')

                    # Check if job with this job_name is already in all_jobs_data
                    if any(job.get("job_name") == job_name for job in all_jobs_data):
                        duplicates += 1
                        continue

                    summary_path = os.path.join(root, summary_file)

                    # Read the summary file
                    try:
                        with open(summary_path) as f:
                            data = json.load(f)
                            # Add job folder name to the data
                            data["job_name"] = job_name
                            all_jobs_data.append(data)
                    except (json.JSONDecodeError, IOError) as e:
                        print(f"Error reading {summary_path}: {e}")
    print(f"Duplicate jobs: {duplicates}")
    return all_jobs_data

def remove_0(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    """Remove all rows where the specified columns have a value of 0.

    This function filters a DataFrame to keep only rows where the values in the specified columns are greater than 0.

    Args:
        df (pd.DataFrame): DataFrame to filter
        cols (List[str]): List of column names to check for 0 values

    Returns:
        pd.DataFrame: Filtered DataFrame with rows where specified columns have non-zero values
    """
    df[cols] = df[df[cols] > 0][cols]
    return df.dropna()

def clean_results(df: pd.DataFrame) -> pd.DataFrame:
    """Clean the AF results by removing rows with unreasonable values.

    Ensures that ranking score is in range [-100,1.5] (see https://alphafoldserver.com/faq#how-do-i-interpret-all-the-outputs-in-the-downloaded-json-files)

    Args:
        df (pd.DataFrame): DataFrame containing AlphaFold results

    Returns:
        pd.DataFrame: Cleaned DataFrame
    """
    initial_count = len(df)
    print(f"Initial dataset size: {initial_count} rows")

    ret = df.copy(deep=True)

    # clean ranking score
    ret = ret[(ret['ranking_score'] >= -100) & (ret['ranking_score'] <= 1.5)]
    ranking_score_count = len(ret)
    print(f"After filtering ranking_score [-100,1.5]: {ranking_score_count} rows ({initial_count - ranking_score_count} removed)")

    ret = ret[(ret['iptm'] >= 0) & (ret['iptm'] <= 1)]
    iptm_count = len(ret)
    print(f"After filtering iptm [0,1]: {iptm_count} rows ({ranking_score_count - iptm_count} removed)")

    ret = ret[(ret['ptm'] >= 0) & (ret['ptm'] <= 1)]
    final_count = len(ret)
    print(f"After filtering ptm [0,1]: {final_count} rows ({iptm_count - final_count} removed)")

    total_removed = initial_count - final_count
    print(f"Total rows removed: {total_removed} ({total_removed/initial_count*100:.2f}%)")

    return ret

def get_job_name(id_tf: str, id_arm: str, df: pd.DataFrame):
    row_arm = df[df['Entry'] == id_arm.lower()]
    if row_arm.empty:
        raise Exception(f"Not found in df: {id_arm}")
    row_tf = df[df['Entry'] == id_tf.lower()]
    if row_tf.empty:
        raise Exception(f"Not found in df: {id_tf}")
    length_arm = row_arm['Length'].iloc[0]
    length_tf = row_tf['Length'].iloc[0]
    return str.lower(f"{id_arm}_1-{length_arm}_{id_tf}_1-{length_tf}")

def get_all_chain_mappings(native_chains, model_chains) -> List[dict]:
    """native -> model
    """
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

def append_dockq_single_interface(df: pd.DataFrame, native_path_prefix: str, model_results_dir: str,
                          all_uniprot: pd.DataFrame, pathmode: str, pdb_cache: str, dockq_cache: str) -> tuple[pd.DataFrame, list[str]]:
    """Calculate DockQ scores for protein structures and append results to dataframe.
    Only usable if the structure has a single Interface!

    This function calculates DockQ scores for all structures in the input dataframe by comparing
    predicted models with native PDB structures using all possible chain mappings.

    Args:
        df (pd.DataFrame): DataFrame containing structure information with columns 'query_tf', 'query_arm', 'pdb_id'
        native_path_prefix (str): Path prefix for native PDB structure files
        model_results_dir (str): Directory containing predicted model files
        all_uniprot (pd.DataFrame): UniProt dataframe for getting job names (only needed if pathmode='uniprot')
        pathmode (str): Mode for path resolution ('uniprot' or 'pdb')
        pdb_cache (str): Directory for PDB file cache in .pdb format
        dockq_cache (str): Directory for DockQ results cache

    Returns:
        pd.DataFrame: Input dataframe with added 'chain_map' and 'dockq_score' columns
    """
    # Create cache directory if it doesn't exist
    os.makedirs(dockq_cache, exist_ok=True)

    # Create copies to avoid modifying the original dataframe
    result_df = df.copy()
    result_df['chain_map'] = None
    result_df['dockq_score'] = None
    result_df['job_name'] = None
    result_df['dockq_complete'] = None

    no_model = []
    no_native = []

    count = 0

    # Process each structure
    for index, row in result_df.iterrows():
        if count % 10 == 0:
            print(f"Processed {count} of {len(result_df)} rows.")
        count += 1

        native_path_cif = f'{native_path_prefix}{row["pdb_id"].lower()}.cif'
        if not os.path.exists(native_path_cif):
            no_native.append(native_path_cif)
            continue

        if pathmode == 'uniprot':
            model_path_cif, job_name = get_model_path_uniprot(row, model_results_dir, all_uniprot)
        elif pathmode == 'pdb':
            model_path_cif, job_name = get_model_path_pdb(row["pdb_id"].lower(), model_results_dir)
        else:
            raise Exception(f"Invalid pathmode: {pathmode}")

        if not model_path_cif:
            no_model.append((model_path_cif, job_name))
            continue

        # Create cache key based on model and native paths
        cache_key = hashlib.md5(f"{str(job_name)}_complete".encode()).hexdigest()
        cache_file = os.path.join(dockq_cache, f"dockq_{cache_key}.pkl")

        # Check if results are cached
        if os.path.exists(cache_file):
            try:
                with open(cache_file, 'rb') as f:
                    cached_result = pickle.load(f)
                result_df.at[index, 'chain_map'] = cached_result['chain_map']
                result_df.at[index, 'dockq_score'] = cached_result['dockq_score']
                result_df.at[index, 'job_name'] = job_name
                result_df.at[index, 'dockq_complete'] = cached_result['dockq_complete']
                print(f"used cached entry for: {job_name}")
                continue
            except (pickle.PickleError, IOError) as e:
                print(f"Error loading cache file {cache_file}: {e}")

        native = load_PDB(native_path_cif)
        try:
            model = load_PDB(get_pdb_from_cif(model_path_cif, pdb_cache))
        except PDBIOException as e:
            print(f"Exception: {e}")
            continue

        native_chains = [chain.id for chain in model]
        model_chains = [chain.id for chain in native]
        chain_map_dict = {}

        if len(model_chains) > 2:
            print(f"{row['pdb_id']}: Warning: Native structure ({native}) has more than two chains!")

        # FIXME: why do we need to check all chain mappings?
        for chain_map in get_all_chain_mappings(model_chains, native_chains):
            try:
                dockQ_complete = run_on_all_native_interfaces(model, native, chain_map=chain_map)
            except Exception as e:
                print(f"Exception for {row['pdb_id']}: {e}.")
                break
            chain_map_dict[str(chain_map)] = dockQ_complete

        if chain_map_dict:
            best_chain_map = max(chain_map_dict.keys(), key=(lambda key: chain_map_dict[key][1]))
            best_dockq_score = chain_map_dict[best_chain_map][1] # here we can use the second value since there is only one interface
            best_dockq_complete = chain_map_dict[best_chain_map]

            # Store results in the DataFrame
            result_df.at[index, 'chain_map'] = best_chain_map
            result_df.at[index, 'dockq_score'] = best_dockq_score
            result_df.at[index, 'dockq_complete'] = best_dockq_complete
            result_df.at[index, 'job_name'] = job_name

            # Cache the results
            cache_data = {
                'chain_map': best_chain_map,
                'dockq_score': best_dockq_score,
                'dockq_complete': best_dockq_complete
            }
            try:
                with open(cache_file, 'wb') as f:
                    pickle.dump(cache_data, f)
            except IOError as e:
                print(f"Error saving cache file {cache_file}: {e}")

    # Report missing files
    if no_model:
        print(f"Missing model files: {len(no_model)}")
        print(str([(t[0], t[1]) for t in no_model]))
    if no_native:
        print("Missing native files:")
        for path in no_native:
            print(f"  NATIVE: {path}")

    return result_df, [str(t[1]) for t in no_model]

def get_model_path_uniprot(row: pd.Series, model_results_dir: str, uniprot_df: pd.DataFrame):
    job_name = get_job_name(row['query_tf'].split('|')[0], row['query_arm'].split('|')[0], uniprot_df)
    model_path = get_file_path(f'{job_name}_model.cif', model_results_dir)

    return model_path, job_name

def get_model_path_pdb(job_name: str, model_results_dir: str):
    return get_file_path(f'{job_name}_model.cif', model_results_dir), job_name

def get_pdb_from_cif(cif_path: Path, pdb_dir: str) -> str:
    """for the .cif structure at cif_path, check if there is a pdb file available in pdb_dir.
    If not, convert the .cif to a pdb file and store it in pdb_dir. Return the path to the pdb file

    Args:
        cif_path (str): Path to the input CIF file
        pdb_dir (str): Directory to store/find PDB files

    Returns:
        str: Path to the PDB file
    """
    # Create pdb_dir if it doesn't exist
    os.makedirs(pdb_dir, exist_ok=True)

    # Extract filename without extension from cif_path
    cif_filename = os.path.basename(cif_path)
    pdb_filename = os.path.splitext(cif_filename)[0] + '.pdb'
    pdb_path = os.path.join(pdb_dir, pdb_filename)

    # Check if PDB file already exists
    if os.path.exists(pdb_path):
        return pdb_path

    # Convert CIF to PDB
    parser = MMCIFParser(QUIET=False)
    structure = parser.get_structure("structure", cif_path)

    # Write to PDB
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_path)

    return pdb_path

def annotate_AF_metrics(report_df: pd.DataFrame, results_dir: str) -> pd.DataFrame:
    """Annotate a report DataFrame with AlphaFold metrics from job results.

    This function takes a report DataFrame and annotates it with AlphaFold metrics
    (iptm, ptm, ranking_score) by matching entries with completed jobs in the results directory.

    Args:
        report_df (pd.DataFrame): DataFrame to annotate, should contain columns needed to construct job names
        results_dir (str): Directory containing AlphaFold job results

    Returns:
        pd.DataFrame: Report DataFrame with added af_iptm, af_ptm, af_ranking_score columns
    """
    report_df['job_name'] = report_df['job_name'].apply(lambda x: str(x).lower() if not pd.isna(x) else x)
    results_df = pd.DataFrame(data=find_summary_files([results_dir]))
    results_df = clean_results(results_df)

    report_df = report_df.merge(results_df, on='job_name')

    return report_df

def append_dockq_two_chainIDs(df: pd.DataFrame, native_path_prefix: str, model_results_dir: str, job_dir: str,
                          pdb_cache: str, dockq_cache: str) -> tuple[pd.DataFrame, list[str]]:
    """Calculate DockQ scores for protein structures and append results to dataframe.

    This function calculates DockQ scores for all structures in the input dataframe by comparing
    predicted models with native PDB structures using all possible chain mappings.

    Args:
        df (pd.DataFrame): DataFrame containing structure information with columns 'query_tf', 'query_arm', 'pdb_id'
        native_path_prefix (str): Path prefix for native PDB structure files
        model_results_dir (str): Directory containing predicted model files
        pdb_cache (str): Directory for PDB file cache in .pdb format
        dockq_cache (str): Directory for DockQ results cache

    Returns:
        pd.DataFrame: Input dataframe with added 'chain_map' and 'dockq_score' columns
    """
    # create cache dir if it does not exist
    os.makedirs(dockq_cache, exist_ok=True)

    # Create copies to avoid modifying the original dataframe
    result_df = df.copy()
    result_df['chain_map'] = None
    result_df['dockq_score'] = None
    result_df['dockq_complete'] = None

    no_model = []

    count = 0
    for index, row in result_df.iterrows():

        job_name = row['job_name']
        pdb_id = row['pdb_id']

        print(job_name)

        if count % 10 == 0:
            print(f"Processed {count} of {len(result_df)} rows.")
        count += 1

        native_path_cif = f'{native_path_prefix}{pdb_id.lower()}.cif'
        if not os.path.exists(native_path_cif):
            raise Exception("no native strucutre found")

        print(job_name.lower())
        model_path_cif, _ = get_model_path_pdb(job_name.lower(), model_results_dir)
        if not model_path_cif:
            no_model.append((model_path_cif, job_name))
            continue

        cache_key = hashlib.md5(f"{str(job_name)}_two_chains".encode()).hexdigest()
        cache_file = os.path.join(dockq_cache, f"dockq_{cache_key}.pkl")

        if os.path.exists(cache_file):
            with open(cache_file, 'rb') as f:
                cached_result = pickle.load(f)
            result_df.at[index, 'chain_map'] = cached_result['chain_map']
            result_df.at[index, 'dockq_score'] = cached_result['dockq_score']
            result_df.at[index, 'dockq_complete'] = cached_result['dockq_complete']
            print(f"used cached entry for: {job_name}")
            continue

        native_struct = load_PDB(native_path_cif)
        try:
            model_struct = load_PDB(get_pdb_from_cif(model_path_cif, pdb_cache))
        except PDBIOException as e:
            print(f"Exception: {e}")
            continue

        try:
            model_chains, native_chains = find_job_chains_in_pdb(job_name, pdb_id, job_dir)
        except FileNotFoundError as e:
            print(e)
            continue

        if len(model_chains) > 2:
            raise Exception("Only works for models with two chains!")

        chain_map_dict = {}

        for chain_map in get_all_chain_mappings(native_chains, model_chains):
            try:
                native_chains = []
                model_chains = []
                for native_id, model_id in chain_map.items():
                    native_chains.append(native_struct[native_id])
                    model_chains.append(model_struct[model_id])
                dockQ_complete = run_on_chains(tuple(model_chains), tuple(native_chains), False, False, False, False)
            except Exception as e:
                print(f"Exception for {row['pdb_id']}: {repr(e)}.")
                break
            if dockQ_complete:
                chain_map_dict[str(chain_map)] = dockQ_complete

        best_dockq_score = 0
        if chain_map_dict:
            best_chain_map = max(chain_map_dict.keys(), key=(lambda key: chain_map_dict[key]['DockQ']))
            best_dockq_score = chain_map_dict[best_chain_map]['DockQ'] # here we can use the total dockQ since there is only one interface (onyl two chains!)
            best_dockq_complete = chain_map_dict[best_chain_map]

            # Store results in the DataFrame
            result_df.at[index, 'chain_map'] = best_chain_map
            result_df.at[index, 'dockq_score'] = best_dockq_score
            result_df.at[index, 'dockq_complete'] = best_dockq_complete

            # Cache the results
            cache_data = {
                'chain_map': best_chain_map,
                'dockq_score': best_dockq_score,
                'dockq_complete': best_dockq_complete
            }
            with open(cache_file, 'wb') as f:
                pickle.dump(cache_data, f)

    # Report missing files
    if no_model:
        print(f"Missing model files: {len(no_model)}")
        print(str([(t[0], t[1]) for t in no_model]))

    return result_df, [str(t[1]) for t in no_model]

def get_job_from_jobname(job_name, dir):
    """search dir for the file named <job_name>.json (case insensitive)
    read the job and return the json

    Args:
        job_name (str): name of the job file (without .json extension)
        dir (str): directory path to search in

    Returns:
        dict: parsed JSON content of the job file

    Raises:
        FileNotFoundError: if no matching job file is found
    """
    dir_path = Path(dir)
    target_filename = f"{job_name}.json"

    # Search for the file case-insensitively
    for file_path in dir_path.rglob("*.json"):
        if file_path.name.lower() == target_filename.lower():
            with open(file_path, 'r') as f:
                return json.load(f)

    raise FileNotFoundError(f"No job file found for job name '{job_name}' in directory '{dir}'")

def find_job_chains_in_pdb(job_name: str, pdb_id: str, dir: str) -> tuple[set,set]:
    """Find the chain IDs in the PDB of the sequences modeled in the AF job
    Return the job chain IDs and the PDB chain IDs

    Args:
        job_name (_type_): _description_
        pdb_id (str): PDB ID of job

    Returns:
        tuple: (model/job chains,native/pdb chains)
    """
    job = get_job_from_jobname(job_name, dir)
    pdb_chains = download_pdb_sequence(pdb_id)
    job_chain_ids = set()
    pdb_chain_ids = set()
    for seq in job['sequences']:
        job_chain_ids.add(seq['protein']['id'])
        sequence_job = seq['protein']['sequence']
        for chain_pdb in pdb_chains:
            if sequence_job == chain_pdb['sequence']:
                pdb_chain_ids.add(chain_pdb['chain_id'])

    return job_chain_ids, pdb_chain_ids

def check_overlap_intf_diso(interface_list: list, disorder_list: list, min_length: int) -> bool:
    """Check if any sequence of consecutive True values in interface_list overlaps 
    with at least min_length consecutive True values in disorder_list.

    Args:
        interface_list (list): List of boolean values representing interface positions
        disorder_list (list): List of boolean values representing disorder positions
        min_length (int): Minimum required length of consecutive overlap

    Returns:
        bool: True if overlap of at least min_length is found, False otherwise

    Raises:
        ValueError: If lists have different lengths or are empty
    """
    # Check for errors
    if len(interface_list) != len(disorder_list):
        raise ValueError("Lists must have the same length")

    if len(interface_list) == 0:
        raise ValueError("Lists cannot be empty")

    # Warning if min_length is too large
    if min_length > len(interface_list):
        print(f"Warning: min_length ({min_length}) is larger than list length ({len(interface_list)})")
        return False

    # Find all consecutive True sequences in interface_list
    interface_sequences = []
    start = None

    for i, val in enumerate(interface_list):
        if val and start is None:
            start = i
        elif not val and start is not None:
            interface_sequences.append((start, i - 1))
            start = None

    # Don't forget the last sequence if it ends at the list end
    if start is not None:
        interface_sequences.append((start, len(interface_list) - 1))

    # For each interface sequence, check overlap with disorder_list
    for seq_start, seq_end in interface_sequences:
        # Check each possible starting position within this sequence
        for start_pos in range(seq_start, seq_end - min_length + 2):
            end_pos = start_pos + min_length - 1

            # Make sure we don't go beyond the sequence end
            if end_pos > seq_end:
                break

            # Check if all positions in this range are True in disorder_list
            if all(disorder_list[pos] for pos in range(start_pos, end_pos + 1)):
                return True

    return False

def test_check_overlap():
    """Test function for check_overlap with various scenarios"""

    # Test 1: Basic overlap case - should return True
    interface_list = [False, True, True, True, False, False]
    disorder_list = [False, False, True, True, True, False]
    assert check_overlap_intf_diso(interface_list, disorder_list, 2) == True, "Test 1 failed: Should find overlap of length 2"

    # Test 2: No overlap - should return False
    interface_list = [True, True, False, False, False, False]
    disorder_list = [False, False, False, True, True, True]
    assert check_overlap_intf_diso(interface_list, disorder_list, 2) == False, "Test 2 failed: Should find no overlap"

    # Test 3: Exact minimum length overlap - should return True
    interface_list = [False, True, True, True, False]
    disorder_list = [False, False, True, True, False]
    assert check_overlap_intf_diso(interface_list, disorder_list, 2) == True, "Test 3 failed: Should find exact minimum overlap"

    # Test 4: Overlap shorter than minimum - should return False
    interface_list = [False, True, True, False, False]
    disorder_list = [False, False, True, False, False]
    assert check_overlap_intf_diso(interface_list, disorder_list, 2) == False, "Test 4 failed: Overlap too short"

    # Test 5: Multiple interface sequences, one overlaps - should return True
    interface_list = [True, False, False, True, True, True]
    disorder_list = [False, False, False, False, True, True]
    assert check_overlap_intf_diso(interface_list, disorder_list, 2) == True, "Test 5 failed: Should find overlap in second sequence"

    # Test 6: All True lists - should return True
    interface_list = [True, True, True, True]
    disorder_list = [True, True, True, True]
    assert check_overlap_intf_diso(interface_list, disorder_list, 3) == True, "Test 6 failed: All True should overlap"

    # Test 7: Single element overlap with min_length 1 - should return True
    interface_list = [False, True, False]
    disorder_list = [False, True, False]
    assert check_overlap_intf_diso(interface_list, disorder_list, 1) == True, "Test 7 failed: Single element overlap"

    # Test 8: Interface sequence at end of list - should return True
    interface_list = [False, False, True, True, True]
    disorder_list = [False, False, False, True, True]
    assert check_overlap_intf_diso(interface_list, disorder_list, 2) == True, "Test 8 failed: End sequence overlap"

    print("All basic tests passed!")

def test_check_overlap_edge_cases():
    """Test edge cases for check_overlap function"""

    # Test error cases
    try:
        check_overlap_intf_diso([], [], 1)
        assert False, "Should raise ValueError for empty lists"
    except ValueError as e:
        assert "cannot be empty" in str(e), "Wrong error message for empty lists"

    try:
        check_overlap_intf_diso([True, False], [True], 1)
        assert False, "Should raise ValueError for different length lists"
    except ValueError as e:
        assert "same length" in str(e), "Wrong error message for different lengths"

    # Test warning case (min_length too large)
    interface_list = [True, True]
    disorder_list = [True, True]
    result = check_overlap_intf_diso(interface_list, disorder_list, 5)
    assert result == False, "Should return False when min_length is too large"

    # Test min_length equal to list length
    interface_list = [True, True, True]
    disorder_list = [True, True, True]
    assert check_overlap_intf_diso(interface_list, disorder_list, 3) == True, "Should work when min_length equals list length"

    # Test with min_length of 0
    interface_list = [True, False]
    disorder_list = [False, True]
    assert check_overlap_intf_diso(interface_list, disorder_list, 0) == True, "Should return True for min_length 0 with any interface"

    # Test all False lists
    interface_list = [False, False, False]
    disorder_list = [False, False, False]
    assert check_overlap_intf_diso(interface_list, disorder_list, 1) == False, "Should return False for all False lists"

    # Test interface True but disorder False
    interface_list = [True, True, True]
    disorder_list = [False, False, False]
    assert check_overlap_intf_diso(interface_list, disorder_list, 1) == False, "Should return False when disorder is all False"

    print("All edge case tests passed!")

def test_check_overlap_complex_cases():
    """Test complex scenarios for check_overlap function"""

    # Test 1: Multiple short interface sequences, none overlap enough
    interface_list = [True, False, True, False, True, True, False]
    disorder_list = [False, True, False, True, False, True, True]
    assert check_overlap_intf_diso(interface_list, disorder_list, 2) == False, "Test 1 failed: Found incorrect overlap"

    # Test 2: Long sequences with partial overlaps
    interface_list = [True, True, True, True, False, False, True, True]
    disorder_list = [False, True, True, False, False, True, True, True]
    assert check_overlap_intf_diso(interface_list, disorder_list, 2) == True, "Test 2 failed: Should find overlap in middle"

    # Test 3: Interface sequence spans entire list
    interface_list = [True] * 10
    disorder_list = [False] * 5 + [True] * 5
    assert check_overlap_intf_diso(interface_list, disorder_list, 3) == True, "Test 3 failed: Should find overlap in second half"

    # Test 4: Alternating pattern
    interface_list = [True, False] * 5
    disorder_list = [False, True] * 5
    assert check_overlap_intf_diso(interface_list, disorder_list, 1) == False, "Test 4 failed: Alternating should not overlap"

    # Test 5: Adjacent but not overlapping
    interface_list = [True, True, False, False, False]
    disorder_list = [False, False, True, True, True]
    assert check_overlap_intf_diso(interface_list, disorder_list, 1) == False, "Test 5 failed: Adjacent should not overlap"

    print("All complex case tests passed!")

if __name__ == "__main__":
    test_check_overlap()
    test_check_overlap_edge_cases()
    test_check_overlap_complex_cases()
    print("All check_overlap tests completed successfully!")
