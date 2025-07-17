import pandas as pd
import os
import json
from typing import List, Dict, Any, Tuple
import matplotlib.pyplot as plt

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