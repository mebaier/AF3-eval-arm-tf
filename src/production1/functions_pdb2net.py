import os
import glob
import pandas as pd
import hashlib
from functions_cif import *

def get_interfaces_pdb2net(path: str, min_atoms: int, max_distance: int, pdb_ids: set):
    """Read all *_detailed_interaction files in path and create a df with the interfaces for each entry fulfilling the desired specs
    Note: pdb2net uses the auth chain IDs (if available)
    print pdb_ids for which no *_detailed_interaction file was found
    
    Args:
        path (str): Path to the directory containing PDB2Net output folders
        min_atoms (int): Minimum number of atoms required for interface
        max_distance (int): Maximum distance for interactions
    """
    
    pdb_ids = set([id.upper() for id in pdb_ids])
    
    if not os.path.exists(path):
        print(f"Path does not exist: {path}")
        return pd.DataFrame(columns=['Entry ID', 'Interface ID', 'Uniprot IDs'])
    
    dir_list = sorted([d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))])
    
    results = []
    
    for i, subfolder in enumerate(dir_list):
        if i % max(1, len(dir_list) // 10) == 0:
            print(f"Processing: {i}/{len(dir_list)}")
        
        entry_id = subfolder.upper()
        subfolder_path = os.path.join(path, subfolder)
            
        interaction_files = glob.glob(os.path.join(subfolder_path, '*_detailed_interactions.csv'))
        if not interaction_files:
            continue
        
        try:
            pdb_ids.remove(entry_id)
        except KeyError as e:
            print(e)
            
        try:
            interactions_df = pd.read_csv(interaction_files[0])
            if interactions_df.empty:
                continue
        except (pd.errors.EmptyDataError, FileNotFoundError, pd.errors.ParserError):
            continue
        
        # Filter interactions for this PDB entry and within distance threshold
        valid_interactions = interactions_df[
            (interactions_df['PDB_ID'] == entry_id) & 
            (interactions_df['Distance'] <= max_distance)
        ]
        
        if valid_interactions.empty:
            continue
        
        # Count interactions per chain pair
        chain_pairs = valid_interactions.groupby(['Chain_A', 'Chain_B']).size()
        
        # Process valid interfaces
        for (chain_a, chain_b), count in chain_pairs.items():
            if count >= min_atoms:
                # Get UniProt IDs for this chain pair
                pair_data = valid_interactions[
                    (valid_interactions['Chain_A'] == chain_a) & 
                    (valid_interactions['Chain_B'] == chain_b)
                ]
                
                if not pair_data.empty:
                    uniprot_a = pair_data['UniProt_A'].iloc[0]
                    uniprot_b = pair_data['UniProt_B'].iloc[0]
                    interface_id = str([str(chain_a), str(chain_b)])
                    
                    results.append({
                        'Entry ID': entry_id,
                        'Interface ID': interface_id,
                        'Uniprot IDs': [uniprot_a, uniprot_b]
                    })

    # Create final DataFrame
    if results:
        result = pd.DataFrame(results)
    else:
        result = pd.DataFrame(columns=['Entry ID', 'Interface ID', 'Uniprot IDs'])
        
    print(f"No model found for: {[e for e in pdb_ids]}")

    return result

def normalize_uniprot_pair(uniprot_list):
    """Create a normalized, sorted tuple of UniProt IDs for comparison"""
    if any(pd.isna(x) for x in uniprot_list):
        return tuple()
    return tuple(sorted(uniprot_list))

def annotate_interface_pdb(df: pd.DataFrame, pdb2net_dir: str, max_distance) -> pd.DataFrame:
    """annotate residues that are involved in interface"""
    df['interface_tf'] = pd.NA

    for ind_df, row_df in df.iterrows():
        job_name = row_df['job_name']

        filename = job_name.upper() + "_MODEL_detailed_interactions.csv"
        # Recursively search for the file
        filepath = None
        for root, _, files in os.walk(pdb2net_dir):
            if filename in files:
                filepath = os.path.join(root, filename)
                break

        if filepath is None:
            raise Exception(f"no file found for {job_name} {filename}")

        interactions_df = pd.read_csv(filepath)

        valid_interactions = interactions_df[
            (interactions_df['PDB_ID'] == filename.split('_detailed_interactions.csv')[0]) &
            (interactions_df['Distance'] <= max_distance)
        ]

        if valid_interactions.empty:
            continue

        map_resnum_pdb_fasta = get_resnum_mapping(filepath) # entity ID B (we need TF)
        map_resnum_pdb_fasta = map_resnum_pdb_fasta['2']

        interface_list = [False] * len(row['Sequence_tf'])
        for _, row in interactions_df.iterrows():
            cif_num = row['Residue_B'].split(':')[0]
            fasta_num = map_resnum_pdb_fasta[cif_num]
            interface_list[fasta_num-1] = True

        df.at[ind_df, 'interface_tf'] = interface_list

    return df

def annotate_interface_tf(df: pd.DataFrame, pdb2net_dir: str, max_distance) -> pd.DataFrame:
    """annotate residues that are involved in interface"""
    df = df.copy()
    df['interface_tf'] = pd.NA

    for ind_df, row_df in df.iterrows():
        print("---")
        print(row_df['Entry_tf'])
        if row_df['Entry_tf'] == '':
            continue

        job_name = row_df['job_name']

        filename = job_name.upper() + "_MODEL_detailed_interactions.csv"
        # Recursively search for the file
        filepath = None
        for root, _, files in os.walk(pdb2net_dir):
            if filename in files:
                filepath = os.path.join(root, filename)
                break

        if filepath is None:
            raise Exception(f"no file found for {job_name} {filename}")

        interactions_df = pd.read_csv(filepath)

        print(interactions_df)

        valid_interactions = interactions_df[
            (interactions_df['Distance'] <= max_distance)
        ]

        valid_interactions = valid_interactions.drop_duplicates(subset=['Residue_B'])

        if valid_interactions.empty:
            print('no valid interactions')
            continue

        interface_list = [False] * len(row_df['iupred3_tf'])
        for _, row in valid_interactions.iterrows():
            cif_num = int(row['Residue_B'].split(':')[0])
            interface_list[cif_num-1] = True

        df.at[ind_df, 'interface_tf'] = interface_list

    return df
