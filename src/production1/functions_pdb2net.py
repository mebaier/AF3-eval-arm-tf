import os
import glob
import pandas as pd
import pickle
from functions_cif import *
from pandas.errors import EmptyDataError

def get_interfaces_pdb2net(path: str, min_atoms: int, max_distance: int, report_df: pd.DataFrame):
    """Read all *_detailed_interaction files in path and create a df with the interfaces for each entry fulfilling the desired specs
    Note: pdb2net uses the auth chain IDs (if available)
    print pdb_ids for which no *_detailed_interaction file was found

    Args:
        path (str): Path to the directory containing PDB2Net output folders
        min_atoms (int): Minimum number of atoms required for interface
        max_distance (int): Maximum distance for interactions
    """

    pdb_ids = set([id.upper() for id in report_df['Entry ID']])

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

        interaction_file = interaction_files[0]

        try:
            pdb_ids.remove(entry_id)
        except KeyError as e:
            print(e)

        try:
            interactions_df = pd.read_csv(interaction_file)
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
            # FIXME: print smth
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

                # row_A = report_df[(report_df['Entry ID'] == entry_id) & (report_df['Auth Asym ID'] == chain_a)].iloc[0]
                # row_B = report_df[(report_df['Entry ID'] == entry_id) & (report_df['Auth Asym ID'] == chain_b)].iloc[0]

                if not pair_data.empty:

                    # interface_mask_A = [False] * len(row_A['iupred3'])
                    # interface_mask_B = [False] * len(row_B['iupred3'])
                    interface_mask_A = [False]
                    interface_mask_B = [False]
                    # for _, row in pair_data.iterrows():
                    #     cif2fasta = get_resnum_mapping_cif2fasta(row_A['model_path'])
                    #     # print(cif2fasta)
                    #     print(row_A['model_path'])
                    #     cif_num_A = int(row['Residue_A'].split(':')[0])
                    #     cif_num_B = int(row['Residue_B'].split(':')[0])
                    #     fast_num_A = cif2fasta[chain_a][cif_num_A]
                    #     fast_num_B = cif2fasta[chain_b][cif_num_B]
                    #     interface_mask_A[int(fast_num_A)-1] = True
                    #     interface_mask_B[int(fast_num_B)-1] = True

                    uniprot_a = pair_data['UniProt_A'].iloc[0]
                    uniprot_b = pair_data['UniProt_B'].iloc[0]
                    interface_id = str([str(chain_a), str(chain_b)])

                    results.append({
                        'Entry ID': entry_id,
                        'Interface ID': interface_id,
                        'Uniprot IDs': [uniprot_a, uniprot_b],
                        'interface_mask_A': interface_mask_A,
                        'interface_mask_B': interface_mask_B,
                    })

    # Create final DataFrame
    if results:
        result = pd.DataFrame(results)
    else:
        result = pd.DataFrame(columns=['Entry ID', 'Interface ID', 'Uniprot IDs', 'interface_mask_diso'])

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

        map_resnum_pdb_fasta = get_resnum_mapping_fasta2cif(filepath) # entity ID B (we need TF)
        map_resnum_pdb_fasta = map_resnum_pdb_fasta['2']

        interface_list = [False] * len(row['Sequence_tf'])
        for _, row in interactions_df.iterrows():
            cif_num = row['Residue_B'].split(':')[0]
            fasta_num = map_resnum_pdb_fasta[cif_num]
            interface_list[fasta_num-1] = True

        df.at[ind_df, 'interface_tf'] = interface_list

    return df

def annotate_interface_tf(df: pd.DataFrame, pdb2net_dir: str, max_distance: int, cache_dir: str = None) -> pd.DataFrame:
    """annotate residues that are involved in interface"""
    df = df.copy()
    df['interface_tf'] = pd.NA

    # Set default cache directory if not provided
    if cache_dir is None:
        cache_dir = os.path.join(os.path.dirname(pdb2net_dir), 'interface_cache')

    # Create cache directory if it doesn't exist
    os.makedirs(cache_dir, exist_ok=True)

    c = 0
    for ind_df, row_df in df.iterrows():
        if c % max(1, len(df)//10) == 0:
            print(f"{c}/{len(df)}")
        c += 1
        
        if row_df['Entry_tf'] == '':
            continue

        job_name = row_df['job_name']

        # Create cache file path using job_name and max_distance
        cache_filename = f"{job_name}_{max_distance}.pkl"
        cache_filepath = os.path.join(cache_dir, cache_filename)

        # Check if cache file exists and load it
        if os.path.exists(cache_filepath):
            try:
                with open(cache_filepath, 'rb') as f:
                    valid_interactions = pickle.load(f)
            except (pickle.PickleError, IOError) as e:
                print(f"Error loading cache file {cache_filepath}: {e}")
                valid_interactions = None
        else:
            valid_interactions = None

        # If not cached, process the data
        if valid_interactions is None:
            filename = job_name.upper() + "_MODEL_detailed_interactions.csv"
            # Recursively search for the file
            filepath = None
            for root, _, files in os.walk(pdb2net_dir):
                if filename in files:
                    filepath = os.path.join(root, filename)
                    break

            if filepath is None:
                raise Exception(f"no file found for {job_name} {filename}")

            try:
                interactions_df = pd.read_csv(filepath)
            except EmptyDataError as e:
                if os.path.exists(filepath):
                    interactions_df = pd.DataFrame(columns=['Distance', 'Residue_B'])
                else:
                    print(f"EmptyDataError: {e}")
                    print(filepath)
                    continue

            valid_interactions = interactions_df[
                (interactions_df['Distance'] <= max_distance)
            ]

            valid_interactions = valid_interactions.drop_duplicates(subset=['Residue_B'])

            # Save to cache file
            try:
                with open(cache_filepath, 'wb') as f:
                    pickle.dump(valid_interactions, f)
            except (pickle.PickleError, IOError) as e:
                print(f"Error saving cache file {cache_filepath}: {e}")

        interface_list = [False] * len(row_df['disordered_regions_mask_tf'])
        for _, row in valid_interactions.iterrows():
            cif_num = int(row['Residue_B'].split(':')[0])
            interface_list[cif_num-1] = True

        df.at[ind_df, 'interface_tf'] = interface_list

    return df
