import pandas as pd
import os, json, itertools
from typing import List, Dict, Any, Tuple
import matplotlib.pyplot as plt
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBIO

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

def find_summary_files(directory: str) -> List[Dict[str, Any]]:
    """Recursively find summary_confidences.json files.
    
    This function searches through the given directory and all subdirectories to find
    AlphaFold summary confidence files and loads them into a list of dictionaries.
    If a job name is duplicated, don't load the job twice

    Args:
        directory (str): Root directory to search for summary files
        
    Returns:
        List[Dict[str, Any]]: List of dictionaries containing the contents of the summary files
    """
    all_jobs_data = []
    duplicates = 0
    
    for root, dirs, files in os.walk(directory):
        # Check for summary_confidences.json files
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

def print_statistics(pairs: pd.DataFrame, results: pd.DataFrame, boxplot: bool = False, description: str = '') -> Dict[str, Dict[str, float]]:
    """Print important statistics for a dataset of protein pairs.
    
    This function filters the results dataframe by the pairs dataframe and calculates statistics
    (mean, median, standard deviation, min, max) for the iptm, ptm, and ranking score metrics.
    If boxplot=True, it also creates a boxplot visualization for these three scores.

    Args:
        pairs (pd.DataFrame): DataFrame containing protein pairs with 'pair_id' column
        results (pd.DataFrame): DataFrame containing AlphaFold results with 'pair_id' column
        boxplot (bool, optional): Whether to create boxplots. Defaults to False.
        description (str, optional): Description to add to the plot title. Defaults to ''.
        
    Returns:
        Dict[str, Dict[str, float]]: Dictionary with statistics for each metric
    """
    # Ensure pairs DataFrame has pair_id column
    if 'pair_id' not in pairs.columns:
        # Create pair_id from p1_Uniprot and p2_Uniprot if needed
        if 'p1_Uniprot' in pairs.columns and 'p2_Uniprot' in pairs.columns:
            pairs['pair_id'] = pairs.apply(lambda row: tuple(sorted([row['p1_Uniprot'].upper(), row['p2_Uniprot'].upper()])), axis=1)
        else:
            raise ValueError("pairs DataFrame must have 'pair_id' column or ('p1_Uniprot' and 'p2_Uniprot') columns")
    
    # Ensure results DataFrame has pair_id column
    if 'pair_id' not in results.columns:
        # Create pair_id from job_name if needed
        if 'job_name' in results.columns:
            results['pair_id'] = results.apply(create_pair_id, axis=1)
        else:
            raise ValueError("results DataFrame must have 'pair_id' column or 'job_name' column")
    
    # Merge dataframes to get matching pairs
    merged_df = pd.merge(results, pairs, on='pair_id', how='inner')
    
    # Calculate statistics for iptm, ptm, and ranking_score
    metrics = ['iptm', 'ptm', 'ranking_score']
    stats = {}
    
    print(f"Statistics for {len(merged_df)} protein pairs:")
    print("-" * 50)
    
    for metric in metrics:
        mean_val = merged_df[metric].mean()
        median_val = merged_df[metric].median()
        std_val = merged_df[metric].std()
        min_val = merged_df[metric].min()
        max_val = merged_df[metric].max()
        
        stats[metric] = {
            'mean': mean_val,
            'median': median_val,
            'std': std_val,
            'min': min_val,
            'max': max_val
        }
        
        print(f"{metric.upper()}")
        print(f"  Mean: {mean_val:.4f}")
        print(f"  Median: {median_val:.4f}")
        print(f"  Standard Deviation: {std_val:.4f}")
        print(f"  Min: {min_val:.4f}")
        print(f"  Max: {max_val:.4f}")
        print("-" * 50)
    
    if boxplot:
        plt.figure(figsize=(10, 6))
        
        # Reshape data for boxplot
        plot_data = []
        labels = []
        
        for metric in metrics:
            plot_data.append(merged_df[metric])
            labels.append(metric.upper())
        
        # Create boxplot
        box = plt.boxplot(plot_data, patch_artist=True, label=labels)
        
        # Add colors to boxplots
        colors = ['lightblue', 'lightgreen', 'pink']
        for patch, color in zip(box['boxes'], colors):
            patch.set_facecolor(color)
            
        plt.annotate(f'Datapoints: {len(merged_df)}', xy=(0.05, 0.95), xycoords='axes fraction', 
                fontsize=12, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))
        
        plt.title(f'Distribution of AlphaFold Scores: {description}')
        plt.ylabel('Score Value')
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.show()
    
    return stats

def compare_statistics(pairs1: pd.DataFrame, pairs2: pd.DataFrame, results: pd.DataFrame, 
                   label1: str = 'Group 1', label2: str = 'Group 2', description: str = '') -> Tuple[Dict[str, Dict[str, float]], Dict[str, Dict[str, float]]]:
    """Compare statistics between two sets of protein pairs.
    
    This function filters the pairs1 and pairs2 dataframes by the results dataframe and creates
    boxplots comparing iptm, ptm and ranking_score metrics between the two groups.
    It also calculates and prints detailed statistics for each metric for both groups.

    Args:
        pairs1 (pd.DataFrame): First dataframe of protein pairs with 'pair_id' column 
                              or ('p1_Uniprot' and 'p2_Uniprot') columns
        pairs2 (pd.DataFrame): Second dataframe of protein pairs with 'pair_id' column 
                              or ('p1_Uniprot' and 'p2_Uniprot') columns
        results (pd.DataFrame): DataFrame containing AlphaFold results with 'pair_id' column 
                              or 'job_name' column
        label1 (str, optional): Label for the first group. Defaults to 'Group 1'.
        label2 (str, optional): Label for the second group. Defaults to 'Group 2'.
        description (str, optional): Description for the plot title. Defaults to ''.
        
    Returns:
        Tuple[Dict[str, Dict[str, float]], Dict[str, Dict[str, float]]]: Tuple containing two dictionaries 
                                                                         with statistics for each metric for both groups
    """
    # Ensure pairs DataFrame has pair_id column
    for pairs, label in [(pairs1, label1), (pairs2, label2)]:
        if 'pair_id' not in pairs.columns:
            # Create pair_id from p1_Uniprot and p2_Uniprot if needed
            if 'p1_Uniprot' in pairs.columns and 'p2_Uniprot' in pairs.columns:
                pairs['pair_id'] = pairs.apply(lambda row: tuple(sorted([row['p1_Uniprot'].upper(), row['p2_Uniprot'].upper()])), axis=1)
            else:
                raise ValueError(f"{label} DataFrame must have 'pair_id' column or ('p1_Uniprot' and 'p2_Uniprot') columns")
    
    # Ensure results DataFrame has pair_id column
    if 'pair_id' not in results.columns:
        # Create pair_id from job_name if needed
        if 'job_name' in results.columns:
            results['pair_id'] = results.apply(create_pair_id, axis=1)
        else:
            raise ValueError("results DataFrame must have 'pair_id' column or 'job_name' column")
    
    # Merge dataframes to get matching pairs
    merged_df1 = pd.merge(results, pairs1, on='pair_id', how='inner')
    merged_df2 = pd.merge(results, pairs2, on='pair_id', how='inner')
    
    # Calculate statistics for iptm, ptm, and ranking_score
    metrics = ['iptm', 'ptm', 'ranking_score']
    stats1 = {}
    stats2 = {}
    
    print(f"Statistics comparison:")
    print(f"{label1}: {len(merged_df1)} protein pairs")
    print(f"{label2}: {len(merged_df2)} protein pairs")
    print("-" * 50)
    
    # Calculate and print statistics for each metric
    for metric in metrics:
        # Group 1 statistics
        mean_val1 = merged_df1[metric].mean()
        median_val1 = merged_df1[metric].median()
        std_val1 = merged_df1[metric].std()
        min_val1 = merged_df1[metric].min()
        max_val1 = merged_df1[metric].max()
        
        # Group 2 statistics
        mean_val2 = merged_df2[metric].mean()
        median_val2 = merged_df2[metric].median()
        std_val2 = merged_df2[metric].std()
        min_val2 = merged_df2[metric].min()
        max_val2 = merged_df2[metric].max()
        
        stats1[metric] = {
            'mean': mean_val1,
            'median': median_val1,
            'std': std_val1,
            'min': min_val1,
            'max': max_val1
        }
        
        stats2[metric] = {
            'mean': mean_val2,
            'median': median_val2,
            'std': std_val2,
            'min': min_val2,
            'max': max_val2
        }
        
        print(f"{metric.upper()}")
        print(f"  {label1}:")
        print(f"    Mean: {mean_val1:.4f}")
        print(f"    Median: {median_val1:.4f}")
        print(f"    Standard Deviation: {std_val1:.4f}")
        print(f"    Min: {min_val1:.4f}")
        print(f"    Max: {max_val1:.4f}")
        print(f"  {label2}:")
        print(f"    Mean: {mean_val2:.4f}")
        print(f"    Median: {median_val2:.4f}")
        print(f"    Standard Deviation: {std_val2:.4f}")
        print(f"    Min: {min_val2:.4f}")
        print(f"    Max: {max_val2:.4f}")
        print("-" * 50)
    
    # Create comparison boxplots
    fig, axes = plt.subplots(1, 3, figsize=(15, 6))
    
    for i, metric in enumerate(metrics):
        data = [merged_df1[metric], merged_df2[metric]]
        
        # Create boxplot
        bp = axes[i].boxplot(data, patch_artist=True, labels=[label1, label2])
        
        # Add colors to boxplots
        colors = ['lightblue', 'lightgreen']
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
        
        axes[i].set_title(f'{metric.upper()}')
        axes[i].grid(True, linestyle='--', alpha=0.7)
        
    plt.suptitle(f'Comparison of AlphaFold Scores: {description}', fontsize=14)
    
    # Add annotation with sample sizes
    fig.text(0.05, 0.95, f"{label1}: {len(merged_df1)} pairs\n{label2}: {len(merged_df2)} pairs", 
             bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))
    
    plt.tight_layout(rect=(0, 0, 1, 0.96))  # Adjust layout to make room for subtitle
    plt.show()
    
    return stats1, stats2

def clean_results(df: pd.DataFrame) -> pd.DataFrame:
    """Clean the AF results by removing rows with unreasonable values.
    
    Ensures that ranking score is in range [0,1].

    Args:
        df (pd.DataFrame): DataFrame containing AlphaFold results

    Returns:
        pd.DataFrame: Cleaned DataFrame
    """
    ret = df.copy(deep=True)
    
    # clean ranking score
    ret = df[(df['ranking_score'] >= 0) & (df['ranking_score'] <= 1)]
    
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

def append_dockq(df: pd.DataFrame, native_path_prefix: str, model_results_dir: str, 
                          all_uniprot: pd.DataFrame, pathmode: str, pdb_cache: str, dockq_cache: str) -> pd.DataFrame:
    """Calculate DockQ scores for protein structures and append results to dataframe.

    This function calculates DockQ scores for all structures in the input dataframe by comparing
    predicted models with native PDB structures using all possible chain mappings.

    Args:
        df (pd.DataFrame): DataFrame containing structure information with columns 'query_tf', 'query_arm', 'pdb_id'
        native_path_prefix (str): Path prefix for native PDB structure files
        model_results_dir (str): Directory containing predicted model files
        all_uniprot (pd.DataFrame): UniProt dataframe for getting job names
        pathmode (str): Mode for path resolution ('uniprot' or 'pdb')
        pdb_cache (str): Directory for PDB file cache
        dockq_cache (str): Directory for DockQ results cache

    Returns:
        pd.DataFrame: Input dataframe with added 'chain_map' and 'dockq_score' columns
    """
    import pickle
    import hashlib
    from DockQ.DockQ import load_PDB, run_on_all_native_interfaces

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

    # Process each structure
    for index, row in result_df.iterrows():
        native_path_cif = f'{native_path_prefix}{row["pdb_id"].lower()}.cif'
        if not os.path.exists(native_path_cif):
            no_native.append((native_path_cif, job_name))
            continue

        if pathmode == 'uniprot':
            model_path_cif, job_name = get_model_path_uniprot(row, model_results_dir, all_uniprot)
        elif pathmode == 'pdb':
            model_path_cif, job_name = get_model_path_pdb(row, model_results_dir)

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
        model = load_PDB(get_pdb_from_cif(model_path_cif, pdb_cache))

        native_chains = [chain.id for chain in model]
        model_chains = [chain.id for chain in native]
        chain_map_dict = {}

        if len(model_chains) > 2:
            print(f"{row['pdb_id']}: Warning: Native structure ({native}) has more than two chains!")

        for chain_map in get_all_chain_mappings(model_chains, native_chains):
            try:
                dockQ_complete = run_on_all_native_interfaces(model, native, chain_map=chain_map)
            except Exception as e:
                print(f"Exception for {row['pdb_id']}: {e}.")
                break
            chain_map_dict[str(chain_map)] = dockQ_complete

        if chain_map_dict:
            best_chain_map = max(chain_map_dict.keys(), key=(lambda key: chain_map_dict[key][1]))
            best_dockq_score = chain_map_dict[best_chain_map][1]
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

    return result_df

def get_model_path_uniprot(row: pd.Series, model_results_dir: str, uniprot_df: pd.DataFrame):
    job_name = get_job_name(row['query_tf'].split('|')[0], row['query_arm'].split('|')[0], uniprot_df)
    model_path = get_file_path(f'{job_name}_model.cif', model_results_dir)

    return model_path, job_name

def get_model_path_pdb(row: pd.Series, model_results_dir: str):
    return get_file_path(f'{row['pdb_id'].lower()}_model.cif', model_results_dir), row['pdb_id']

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