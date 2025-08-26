import pandas as pd
import os, shutil, json, math, random, re
from typing import List, Tuple, Dict, Any, Union
from download_functions import download_pdb_sequence

def write_af_jobs_to_individual_files(af_jobs: List[Dict[str, Any]], output_dir: str, dialect='alphafold3') -> None:
    """Write each AlphaFold job to an individual file.
    
    Each job is saved as a JSON file in the specified output directory with the job's name as the filename.
    
    Args:
        af_jobs (List[Dict[str, Any]]): List of AlphaFold job dictionaries
        output_dir (str): Directory where job files will be saved
        
    Returns:
        None
        
    Raises:
        SystemExit: If output directory already exists
    """
    if os.path.exists(output_dir):
        print(f"ERROR: output directory '{output_dir}' already exists!")
        print("Aborting.")
        return
    os.makedirs(output_dir, exist_ok=False)
    for job in af_jobs:
        file_name = f"{job['name']}.json"
        file_path = os.path.join(output_dir, file_name)
        with open(file_path, 'w') as f:
            if dialect == 'alphafoldserver':
                json.dump([job], f, indent=2)  # use [] so AF parser knows it's in alphafoldserver dialect
            elif dialect == 'alphafold3':
                json.dump(job, f, indent=2)
            else:
                raise Exception('Incorrect dialect.')
    
def get_comparable_job(job_data: Dict[str, Any]) -> Any:
    """Create a comparable representation of the job by extracting the used sequences.
    
    This function creates a standardized representation of an AlphaFold job that can be compared
    to other jobs to detect duplicates, regardless of the job name or field order. 
    IMPORTANT: ignores modifications of sequences etc.!
    
    Args:
        job_data (Dict[str, Any]): AlphaFold job dictionary

    Returns:
        Any: Comparable representation of the job
    """
    comparable = []
    if job_data['dialect'] == 'alphafoldserver':
        # Extract sequences from alphafoldserver format
        for seq_entry in job_data['sequences']:
            protein_chain = seq_entry['proteinChain']
            sequence = protein_chain['sequence']
            count = protein_chain['count']
            # Add the sequence 'count' times to the comparable list
            for _ in range(count):
                comparable.append(sequence)

    elif job_data['dialect'] == 'alphafold3':
        # Extract sequences from alphafold3 format
        for seq_entry in job_data['sequences']:
            protein = seq_entry['protein']
            sequence = protein['sequence']
            # Each sequence appears once in alphafold3 format
            comparable.append(sequence)
    else:
        raise Exception(f"Invalid dialect for {job_data}.")
    return sorted(comparable)

def collect_created_jobs(results_dir: str) -> List[Dict[str, Any]]:
    """Collect all jobs in a directory (all .json files are considered jobs).
    
    IMPORTANT: one file is considered to have one job!
    
    Args:
        results_dir (str): Directory containing AlphaFold job files
        
    Returns:
        List[Dict[str, Any]]: List of AlphaFold job dictionaries
        
    Raises:
        OSError: If directory cannot be accessed
    """
    collected_jobs = []

    # Go through all .json files in the results directory (not recursively)
    for file_name in os.listdir(results_dir):
        if file_name.endswith('.json') and os.path.isfile(os.path.join(results_dir, file_name)):
            try:
                with open(os.path.join(results_dir, file_name), 'r') as f:
                    data = json.load(f)
                    if isinstance(data, list):
                        collected_jobs += data
                    else:
                        collected_jobs += [data]
                    
            except (json.JSONDecodeError, IOError) as e:
                print(f"Error reading {file_name}: {e}")
                continue
    return collected_jobs

def create_alphafold_job(job_name: str, sequence1: str, sequence2: str, dialect: str = 'alphafold3', seed: int = 3030452494) -> Dict[str, Any]:
    """Create a standardized AlphaFold job dictionary from input parameters.

    Args:
        job_name (str): Name of the job, typically in the format "protein1_start-end_protein2_start-end"
        sequence1 (str): Amino acid sequence of the first protein
        sequence2 (str): Amino acid sequence of the second protein
        dialect (str, optional): The dialect to use for AlphaFold. Defaults to 'alphafold3'.

    Returns:
        Dict[str, Any]: A dictionary representing the AlphaFold job in the specified format
    """
    if dialect == 'alphafoldserver':
        return {
            'name': job_name,
            'modelSeeds': [],
            'sequences': [
                {
                    'proteinChain': {
                        'sequence': sequence1,
                        'count': 1
                    }
                },
                {
                    'proteinChain': {
                        'sequence': sequence2,
                        'count': 1
                    }
                }
            ],
            'dialect': 'alphafoldserver',
            'version': 1,
        }
    elif dialect == 'alphafold3':
        return {
            'name': job_name,
            'modelSeeds': [seed],
            'sequences': [
                {
                    'protein': {
                        'id': 'A',
                        'sequence': sequence1,
                        'modifications': []
                    }
                },
                {
                    'protein': {
                        'id': 'B',
                        'sequence': sequence2,
                        'modifications': []
                    }
                }
            ],
            'dialect': 'alphafold3',
            'version': 1
        }
    else:
        raise Exception("Incorrect dialect. Use 'alphafoldserver' or 'alphafold3'.")
        
def create_alphafold_job_ms(job_name: str, sequences: List[Dict[str, str]], seed: int = 3030452494) -> Dict[str, Any]:
    """Create an AlphaFold job from PDB sequence data.

    Args:
        job_name (str): Name of the job
        sequences (List[Dict[str, str]]): List of sequence dictionaries from download_pdb_sequence()
                                         Each dict has 'chain_id' and 'sequence' keys
        seed (int, optional): Model seed. Defaults to 3030452494.

    Returns:
        Dict[str, Any]: AlphaFold job dictionary in alphafold3 format

    Raises:
        ValueError: If sequences list is empty
    """
    if len(sequences) == 0:
        raise ValueError("Sequences list cannot be empty")

    # Create sequence entries for all provided sequences
    sequence_entries = []
    for seq_data in sequences:
        chain_id = re.sub(r'\[auth .+\]', '', seq_data['chain_id'])
        sequence_entries.append({
            'protein': {
                'id': chain_id,
                'sequence': seq_data['sequence'],
                'modifications': []
            }
        })

    return {
        'name': job_name,
        'modelSeeds': [seed],
        'sequences': sequence_entries,
        'dialect': 'alphafold3',
        'version': 1
    }

def create_job_batch_scoreCategories(pair_df: pd.DataFrame, batch_size: int, categories: List[Tuple[float, float]], 
                                job_dirs: List[str], column_name: str, token_limit: int = 5120) -> List[Dict[str, Any]]:
    """Create a new batch of batch_size jobs from pair_df.
    
    Don't create jobs that have been created previously.
    For each category (range of scores) in categories create 
    batch_size/len(categories) new jobs randomly sampled from all possible jobs in the category.
    If a category does not have enough possible jobs to fill the limit, redistribute to the other
    categories. Filter jobs that are too large.
    
    Note: categories should be specified in a way so the category with the largest number of possible jobs comes at the end of the array.
    
    Args:
        pair_df (pd.DataFrame): DataFrame containing protein pairs with sequences and other information
        batch_size (int): Total number of jobs to create across all categories
        categories (List[Tuple[float, float]]): List of tuples (min_score, max_score) defining score ranges for each category
        job_dirs (List[str]): List of directories to search for existing jobs to avoid duplicates
        column_name (str): Column name in pair_df for the score to filter on
        token_limit (int, optional): Maximum total sequence length for a job. Defaults to 5120.
    
    Returns:
        List[Dict[str, Any]]: List of newly created AlphaFold job dictionaries
        
    Raises:
        KeyError: If required columns are missing from pair_df
    """
    new_jobs = []
    prev_jobs = []
    for dir in job_dirs:
        prev_jobs += collect_created_jobs(dir)
        
    for i in range(len(prev_jobs)):
        prev_jobs[i] = get_comparable_job(prev_jobs[i])
    
    total_created = 0
    num_categories = len(categories)
    category_counts = {}  # Dictionary to track jobs created in each category
    
    for category in categories:
        # Stop if we've already created enough jobs
        if total_created >= batch_size:
            break
            
        # dynamically adjust the quota for the next category to account for categories that don't fill their quota
        cat_quota = math.floor((batch_size-total_created)/num_categories)
        
        num_categories -= 1
        
        min_score = category[0]
        max_score = category[1]
        created_in_category = 0
        possible_pairs = pair_df[(pair_df[column_name] >= min_score) & (pair_df[column_name] <= max_score)]
        possible_ind = possible_pairs.index.tolist()
        
        category_key = f"{min_score}-{max_score}"
        category_counts[category_key] = 0
        
        while created_in_category < cat_quota and len(possible_ind) > 0:
            ind = random.choice(possible_ind)
            possible_ind.remove(ind)
            
            row = pair_df.iloc[ind]
            armadillo_entry = row['Entry_arm']
            tf_entry = row['Entry_tf']
            armadillo_sequence = row['Sequence_arm']
            tf_sequence = row['Sequence_tf']
            
            if len(tf_sequence) + len(armadillo_sequence) > token_limit:
                continue
            
            # Calculate indices for the sequences
            armadillo_x, armadillo_y = 1, len(armadillo_sequence)
            tf_x, tf_y = 1, len(tf_sequence)

            # Generate job name in the specified format
            job_name = f"{armadillo_entry}_{armadillo_x}-{armadillo_y}_{tf_entry}_{tf_x}-{tf_y}"
            
            # Create job using the helper function
            job = create_alphafold_job(job_name, armadillo_sequence, tf_sequence)
            
            # check if job was already created earlier
            job_comparable = get_comparable_job(job)
            if not any(existing_comparable == job_comparable for existing_comparable in prev_jobs):
                # no duplicate job found
                created_in_category += 1
                total_created += 1
                category_counts[category_key] += 1
                prev_jobs.append(get_comparable_job(job))
                new_jobs.append(job)
                
                # Ensure we don't exceed batch_size
                if total_created >= batch_size:
                    break
        
    # Print the number of jobs created in each category
    print(f"Created {len(new_jobs)} new jobs total.")
    print("Jobs created per category:")
    for category, count in category_counts.items():
        print(f"  Score range {category}: {count} jobs")
    
    return new_jobs

def create_job_batch_all_pairs(pair_df: pd.DataFrame, batch_size: int, 
                               job_dirs: List[str], token_limit: int = 5120) -> List[Dict[str, Any]]:
    """Create a new batch of AlphaFold jobs by randomly sampling from all protein pairs.
    
    This function creates a specified number of new AlphaFold jobs by randomly selecting protein pairs
    from the input DataFrame. It avoids creating duplicate jobs by checking against existing jobs
    in the specified directories and filters out pairs that exceed the token limit.
    
    Args:
        pair_df (pd.DataFrame): DataFrame containing protein pairs with columns:
                               - 'Entry_arm': Armadillo protein entry ID
                               - 'Entry_tf': Transcription factor protein entry ID  
                               - 'Sequence_arm': Armadillo protein sequence
                               - 'Sequence_tf': Transcription factor protein sequence
        batch_size (int): Total number of new jobs to create
        job_dirs (List[str]): List of directories to search for existing jobs to avoid duplicates
        token_limit (int, optional): Maximum total sequence length allowed for a job. 
                                   Pairs with combined sequence length exceeding this limit
                                   will be skipped. Defaults to 5120.
    
    Returns:
        List[Dict[str, Any]]: List of newly created AlphaFold job dictionaries in alphafoldserver format.
                             Each job dictionary contains:
                             - 'name': Job name in format "proteinA_start-end_proteinB_start-end"
                             - 'sequences': List of protein chain specifications
                             - 'dialect': Set to 'alphafoldserver'
                             - 'version': Set to 1
                             - 'modelSeeds': Empty list
    
    Note:
        - The function randomly samples pairs without replacement until the batch size is reached
        - Pairs exceeding the token limit are automatically skipped
        - Duplicate jobs (based on sequence content) are avoided by comparing against existing jobs
        - Job names follow the format: "{armadillo_entry}_1-{seq_length}_{tf_entry}_1-{seq_length}"
        
    Raises:
        KeyError: If required columns are missing from pair_df
        IndexError: If pair_df is empty or insufficient pairs available
    """
    new_jobs = []
    prev_jobs = []
    for dir in job_dirs:
        prev_jobs += collect_created_jobs(dir)
        
    for i in range(len(prev_jobs)):
        prev_jobs[i] = get_comparable_job(prev_jobs[i])
    
    possible_ind = pair_df.index.tolist()   
    
    total_created = 0
    while total_created < batch_size:
        
        ind = random.choice(possible_ind)
        possible_ind.remove(ind)
            
        row = pair_df.iloc[ind]
        armadillo_entry = row['Entry_arm']
        tf_entry = row['Entry_tf']
        armadillo_sequence = row['Sequence_arm']
        tf_sequence = row['Sequence_tf']
            
        if len(tf_sequence) + len(armadillo_sequence) > token_limit:
            continue
            
        # Calculate indices for the sequences
        armadillo_x, armadillo_y = 1, len(armadillo_sequence)
        tf_x, tf_y = 1, len(tf_sequence)

        # Generate job name in the specified format
        job_name = f"{armadillo_entry}_{armadillo_x}-{armadillo_y}_{tf_entry}_{tf_x}-{tf_y}"
        
        # Create job using the helper function
        job = create_alphafold_job(job_name, armadillo_sequence, tf_sequence)
        
        # check if job was already created earlier
        job_comparable = get_comparable_job(job)
        if not any(existing_comparable == job_comparable for existing_comparable in prev_jobs):
            # no duplicate job found
            total_created += 1
            prev_jobs.append(get_comparable_job(job))
            new_jobs.append(job)
        
    # Print the number of jobs created in each category
    print(f"Created {len(new_jobs)} new jobs total.")
    
    return new_jobs

def create_job_batch_id_list(pair_df: pd.DataFrame, id_list: List[Tuple[str, str]], 
                               job_dirs: List[str], token_limit: int = 5120) -> List[Dict[str, Any]]:
    """Create a batch of AlphaFold jobs from a specific list of protein ID pairs.
    
    This function creates AlphaFold jobs for specific protein pairs identified by their UniProt IDs.
    It avoids creating duplicate jobs by checking against existing jobs in the specified directories.
    
    Args:
        pair_df (pd.DataFrame): DataFrame containing protein pairs with columns:
                               - 'Entry_arm': Armadillo protein entry ID
                               - 'Entry_tf': Transcription factor protein entry ID  
                               - 'Sequence_arm': Armadillo protein sequence
                               - 'Sequence_tf': Transcription factor protein sequence
        id_list (List[Tuple[str, str]]): List of tuples containing protein ID pairs to create jobs for
        job_dirs (List[str]): List of directories to search for existing jobs to avoid duplicates
        token_limit (int, optional): Maximum total sequence length allowed for a job. 
                                   Pairs with combined sequence length exceeding this limit
                                   will be skipped. Defaults to 5120.
    
    Returns:
        List[Dict[str, Any]]: List of newly created AlphaFold job dictionaries
        
    Raises:
        Exception: If no matching row is found for a pair or multiple matching rows are found
        KeyError: If required columns are missing from pair_df
    """
    new_jobs = []
    prev_jobs = []
    for dir in job_dirs:
        prev_jobs += collect_created_jobs(dir)
        
    for i in range(len(prev_jobs)):
        prev_jobs[i] = get_comparable_job(prev_jobs[i])
        
    total_created = 0
    
    for (id_1, id_2) in id_list:
        
        row_matches = pair_df[((pair_df['Entry_arm'] == id_1) & (pair_df['Entry_tf'] == id_2)) | 
                              ((pair_df['Entry_arm'] == id_2) & (pair_df['Entry_tf'] == id_1))]
        if len(row_matches) == 0:
            print(f"No matching row found for pair ({id_1}, {id_2}) in pair_df.")
            continue
        elif len(row_matches) > 1:
            print(f"Multiple matching rows found for pair ({id_1}, {id_2}) in pair_df.")
            continue
        
        row = row_matches.iloc[0]
        armadillo_entry = row['Entry_arm']
        tf_entry = row['Entry_tf']
        armadillo_sequence = row['Sequence_arm']
        tf_sequence = row['Sequence_tf']
            
        if len(tf_sequence) + len(armadillo_sequence) > token_limit:
            print(f"WARNING: Skipping because of token limit: {armadillo_entry}-{tf_entry}")
            continue
            
        # Calculate indices for the sequences
        armadillo_x, armadillo_y = 1, len(armadillo_sequence)
        tf_x, tf_y = 1, len(tf_sequence)

        # Generate job name in the specified format
        job_name = f"{armadillo_entry}_{armadillo_x}-{armadillo_y}_{tf_entry}_{tf_x}-{tf_y}"
        
        # Create job using the helper function
        job = create_alphafold_job(job_name, armadillo_sequence, tf_sequence)
        
        # check if job was already created earlier
        job_comparable = get_comparable_job(job)
        if not any(existing_comparable == job_comparable for existing_comparable in prev_jobs):
            # no duplicate job found
            total_created += 1
            prev_jobs.append(get_comparable_job(job))
            new_jobs.append(job)
        else:
            print(f"Duplicate job found for pair ({id_1}, {id_2})")
        
    # Print the number of jobs created in each category
    print(f"Created {len(new_jobs)} new jobs total.")
    
    return new_jobs

def adjust_job_id(job_ids: list, source_results_dirs: List[str], target_results_dir: str) -> None:
    """For the list of job IDs, check if at least one result directory with an ID in the list exists in any of the source_results_dirs.
    For IDs that don't have an existing directory, create these directories in target_results_dir by copying the existing dir and renaming it.
    Note: does not modify the IDs in the job files

    Args:
        job_ids (list): List of job IDs to check and create directories for
        source_results_dirs (List[str]): List of base directories where existing result directories are located
        target_results_dir (str): Base directory where new copied directories should be created
    """

    if not job_ids:
        print("Warning: Empty job_ids list provided")
        return

    if not source_results_dirs:
        print("Warning: Empty source_results_dirs list provided")
        return

    # Find existing directories for any of the job IDs by searching recursively in all source directories
    existing_dirs = []
    missing_ids = []

    for job_id in job_ids:
        target_dir = os.path.join(target_results_dir, job_id)
        source_dir = None

        # Search recursively for the job_id directory in all source_results_dirs
        for source_results_dir in source_results_dirs:
            if not os.path.exists(source_results_dir):
                print(f"Warning: Source directory '{source_results_dir}' does not exist")
                continue
  
            for root, dirs, files in os.walk(source_results_dir):
                if job_id in dirs:
                    source_dir = os.path.join(root, job_id)
                    break

            if source_dir:
                break

        if source_dir and os.path.exists(source_dir):
            existing_dirs.append((job_id, source_dir))
        elif not os.path.exists(target_dir):
            missing_ids.append(job_id)

    if not existing_dirs:
        print(f"Warning: No existing result directories found for any of the job IDs: {job_ids}")
        return

    if not missing_ids:
        print(f"All job IDs already have existing directories in target location: {job_ids}")
        return

    # Create target results directory if it doesn't exist
    os.makedirs(target_results_dir, exist_ok=True)

    # Use the first existing directory as the source for copying
    source_id, source_dir = existing_dirs[0]
    print(f"Using '{source_dir}' as source directory for copying")

    # Create directories for missing IDs by copying the existing one to target directory
    for missing_id in missing_ids:
        target_dir = os.path.join(target_results_dir, missing_id)
        try:
            shutil.copytree(source_dir, target_dir)
            print(f"Created directory '{target_dir}' by copying from '{source_dir}'")

            # Rename files and subdirectories within the copied directory
            rename_paths(target_dir, source_id, missing_id)

        except (OSError, shutil.Error) as e:
            print(f"Error creating directory '{target_dir}': {e}")

    print(f"Directory adjustment completed for job IDs: {job_ids}")

def rename_paths(base_path: str, old_str: str, new_str: str) -> None:
    """
    Recursively rename files and directories under base_path,
    replacing occurrences of old_str with new_str in their names.

    Args:
        base_path (str): The root directory to start renaming.
        old_str (str): The substring to replace in names.
        new_str (str): The substring to replace with.
    """
    if not os.path.exists(base_path):
        print(f"Warning: Base path '{base_path}' does not exist.")
        return

    # Walk the directory tree bottom-up to rename files before their parent directories
    for root, dirs, files in os.walk(base_path, topdown=False):

        # Rename files first
        for filename in files:
            if old_str in filename:
                old_file_path = os.path.join(root, filename)
                new_filename = filename.replace(old_str, new_str)
                new_file_path = os.path.join(root, new_filename)
                try:
                    os.rename(old_file_path, new_file_path)
                    print(f"Renamed file: {old_file_path} -> {new_file_path}")
                except OSError as e:
                    print(f"Error renaming file {old_file_path}: {e}")

        # Then rename directories
        for dirname in dirs:
            if old_str in dirname:
                old_dir_path = os.path.join(root, dirname)
                new_dirname = dirname.replace(old_str, new_str)
                new_dir_path = os.path.join(root, new_dirname)
                try:
                    os.rename(old_dir_path, new_dir_path)
                    print(f"Renamed directory: {old_dir_path} -> {new_dir_path}")
                except OSError as e:
                    print(f"Error renaming directory {old_dir_path}: {e}")

def find_job_files(job_list: list, job_dir: str):
    """for a list of jobs, find the corresponding job files in job_dir and print the path

    Args:
        job_list (list): list of job IDs to find
        job_dir (str): directory to search for job files (searches recursively)
    """

    found_files = []
    missing_jobs = []

    # Create a dictionary to store all .json files found recursively
    json_files = {}
    for root, dirs, files in os.walk(job_dir):
        for file in files:
            if file.endswith('.json'):
                job_id = file[:-5]  # Remove .json extension
                full_path = os.path.join(root, file)
                json_files[job_id] = full_path

    # Search for each job ID
    for job_id in job_list:
        if job_id in json_files:
            found_files.append(json_files[job_id])
            print(f"Found: {json_files[job_id]}")
        else:
            missing_jobs.append(job_id)
            print(f"Missing: {job_id}.json (searched recursively in {job_dir})")

    print(f"\nSummary:")
    print(f"Found {len(found_files)} job files")
    print(f"Missing {len(missing_jobs)} job files")

    if missing_jobs:
        print(f"Missing jobs: {missing_jobs}")

    return found_files, missing_jobs

def rewrite_af_job(af_job: Dict[str, Any]) -> Dict[str, Any]:
    """Convert an AlphaFold job from alphafoldserver dialect to alphafold3 dialect.
    
    This function takes an AlphaFold job in alphafoldserver format and converts it to
    alphafold3 format using the create_alphafold_job helper function.
    
    Args:
        af_job (Dict[str, Any]): AlphaFold job dictionary in alphafoldserver dialect
        
    Returns:
        Dict[str, Any]: AlphaFold job dictionary in alphafold3 dialect
        
    Raises:
        KeyError: If required fields are missing from the input job
        ValueError: If the job doesn't have exactly 2 protein sequences
        Exception: If the input job is not in alphafoldserver dialect
    """
    # Validate input dialect
    if af_job.get('dialect') != 'alphafoldserver':
        raise Exception(f"Input job must be in 'alphafoldserver' dialect, got '{af_job.get('dialect')}'")
    
    # Extract job name
    job_name = af_job['name']
    
    # Extract sequences - expecting exactly 2 protein chains
    sequences = af_job['sequences']
    if len(sequences) != 2:
        raise ValueError(f"Expected exactly 2 protein sequences, got {len(sequences)}")
    
    # Extract the two protein sequences
    sequence1 = sequences[0]['proteinChain']['sequence']
    sequence2 = sequences[1]['proteinChain']['sequence']
    
    # Create new job in alphafold3 dialect
    return create_alphafold_job(job_name, sequence1, sequence2, dialect='alphafold3')

def create_job_batch_from_PDB_IDs(pdb_ids: List, job_dirs: List[str], token_limit: int = 5120):
    """create a job batch from a list of PDB ids.
    Don't create duplicate jobs, don't create jobs that exceed token limit

    Args:
        pdb_ids (List): _description_
        job_dirs (List[str]): _description_
        token_limit (int, optional): _description_. Defaults to 5120.

    Returns:
        _type_: _description_
    """
    new_jobs = []
    prev_jobs = []
    for dir in job_dirs:
        prev_jobs += collect_created_jobs(dir)

    for i in range(len(prev_jobs)):
        prev_jobs[i] = (prev_jobs[i]['name'], get_comparable_job(prev_jobs[i]))

    total_created = 0
    duplicates = 0

    for pdb_id in pdb_ids:

        length = 0

        sequences = download_pdb_sequence(pdb_id)
        if len(sequences) == 0:
            print(f"Empty sequence list: {pdb_id}")
            continue

        for seq_dict in sequences:
            length += len(seq_dict['sequence'])

        if length > token_limit:
            print(f"Skipping because of token limit: {pdb_id} (length: {length})")
            continue

        # Create job using the helper function
        job_name = pdb_id
        job = create_alphafold_job_ms(job_name, sequences)

        # check if job was already created earlier
        job_comparable = get_comparable_job(job)
        if not any(existing_comparable[1] == job_comparable for existing_comparable in prev_jobs):
            # no duplicate job found
            total_created += 1
            prev_jobs.append(get_comparable_job(job))
            new_jobs.append(job)
            continue
        else:
            matches = [ec[0] for ec in prev_jobs if ec[1] == job_comparable]
            if not job_name in matches: # pdb_id is name
                print(f"Skipping duplicate job UNDER DIFFERENT ID: {job_name}. Duplicate IDs: {matches}")
                adjust_job_id([m.lower() for m in matches + [job_name]], ["/home/markus/MPI_local/HPC_results_full"], "/home/markus/MPI_local/HPC_results_full/copied_jobs")
            else:
                print(f"Skipping duplicate job: {job_name}. Duplicate IDs: {matches}")
            duplicates += 1
            continue

    # Print the number of jobs created in each category
    print(f"Skipped {duplicates} duplicate jobs.")
    print(f"Created {len(new_jobs)} new jobs total.")

    return new_jobs