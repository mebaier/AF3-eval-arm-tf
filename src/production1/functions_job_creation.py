import pandas as pd
import os
import json
from typing import List, Tuple, Dict, Any, Union, Optional, Set
import math
import copy
import random

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

def sort_rec(obj: Union[List[Any], Dict[str, Any], Any]) -> Any:
    """Sort a list or dictionary recursively.
    
    This function is used to create comparable job representations by sorting all nested structures.
    
    Args:
        obj (Union[List[Any], Dict[str, Any], Any]): Object to sort recursively

    Returns:
        Any: Sorted object with all nested structures sorted
    """
    if isinstance(obj, dict):
        return sorted((k, sort_rec(v)) for k, v in obj.items())
    if isinstance(obj, list):
        return sorted(sort_rec(x) for x in obj)
    else:
        return obj
    
def get_comparable_job(job_data: Dict[str, Any], deep_copy: bool = True) -> Any:
    """Create a comparable representation of the job by removing the name field and sorting the other fields.
    
    This function creates a standardized representation of an AlphaFold job that can be compared
    to other jobs to detect duplicates, regardless of the job name or field order.
    
    Args:
        job_data (Dict[str, Any]): AlphaFold job dictionary
        deep_copy (bool, optional): Whether to create a deep copy of the job data. Defaults to True.

    Returns:
        Any: Comparable representation of the job
    """
    if deep_copy:
        # Create a deep copy to avoid modifying the original
        comparable = copy.deepcopy(job_data)
    else:
        comparable = job_data
    if 'name' in comparable:
        del comparable['name']
    return sort_rec(comparable)

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
                    collected_jobs += json.load(f)
                    
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
        prev_jobs[i] = get_comparable_job(prev_jobs[i], deep_copy=False)
    
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
        prev_jobs[i] = get_comparable_job(prev_jobs[i], deep_copy=False)
    
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
        prev_jobs[i] = get_comparable_job(prev_jobs[i], deep_copy=False)
        
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