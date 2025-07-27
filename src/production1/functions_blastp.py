import pandas as pd
from typing import List
from pathlib import Path

def clean_blastp_out(df: pd.DataFrame, identity_cutoff :int|bool=False, score_cutoff :int|bool=False, evalue_cutoff :float|bool = False, coverage_cutoff:float|bool=False) -> pd.DataFrame:
    """
    Clean and filter BLAST output DataFrame.
    
    Args:
        df: Raw BLAST output DataFrame
        
    Returns:
        Filtered and cleaned DataFrame with additional columns for PDB ID and chain
    """
    df['pdb_id'] = df['subject'].apply(lambda x: x.split('_')[0] if '_' in x else x)
    df['chain'] = df['subject'].apply(lambda x: x.split('_')[1] if len(x.split('_')) > 1 else '')
    df = df.drop_duplicates(subset='subject')
    df['%identity'] = pd.to_numeric(df['%identity'], errors='raise')
    df['bit score'] = pd.to_numeric(df['bit score'], errors='raise')
    df['evalue'] = pd.to_numeric(df['evalue'], errors='raise')
    df['% query coverage per subject'] = pd.to_numeric(df['% query coverage per subject'], errors='raise')
    if identity_cutoff:
        df = df[df['%identity'] >= identity_cutoff]
    if score_cutoff:
        df = df[df['bit score'] >= score_cutoff]
    if evalue_cutoff:
        df = df[df['evalue'] <= evalue_cutoff]
    if coverage_cutoff:
        df = df[df['% query coverage per subject'] >= coverage_cutoff]
    df['uniprot_id'] = df['query'].apply(lambda x: x.split('|')[0])
    return df

def read_blast_to_df(output_dir: str, column_names: List[str]) -> pd.DataFrame:
    """
    Read BLAST output files from a directory and combine them into a single DataFrame.
    
    Args:
        output_dir: Directory containing BLAST output files
        column_names: List of column names for the DataFrame
        
    Returns:
        Combined and cleaned DataFrame with all BLAST results
    """
    output_files: List[Path] = list(Path(output_dir).glob('*_blastp.out'))
    output_list: List[List[str]] = []
    
    for output_file in output_files:
        try:
            with open(output_file, 'r') as f:
                lines: List[str] = [line for line in f if not line.startswith('#')]
            output_list.append(lines)
        except IOError as e:
            print(f"Error reading file {output_file}: {e}")
            continue
    
    # Flatten the list of lines and split by tab to create rows
    rows: List[List[str]] = [line.strip().split('\t') for file_lines in output_list for line in file_lines if line.strip()]
    
    if not rows:
        print(f"Warning: No data found in {output_dir}")
        return pd.DataFrame(columns=column_names)
    
    return pd.DataFrame(rows, columns=column_names)