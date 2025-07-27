import re
import pandas as pd

def intact_score_filter(x: str) -> float:
    """Parse IntAct miscore from confidence value string.
    
    Args:
        x (str): Confidence value string containing intact-miscore
        
    Returns:
        float: Extracted IntAct miscore, or 0.0 if not found
    """
    match = re.search(r"(?<=intact-miscore:)\d\.\d*", x)
    return float(match.group()) if match else 0.0

def clean_IntAct(intact_df: pd.DataFrame) -> pd.DataFrame:
    """clean intact dataset:
       create new column 'intact_score' that has the score from 'Confidence value(s)' parsed
       drop unnecessary columns

    Args:
        df (pd.DataFrame): _description_

    Returns:
        pd.DataFrame: _description_
    """
    intact_cleaned = intact_df.copy()
    intact_cleaned.drop(['Alias(es) interactor A', 
                        'Alias(es) interactor B', 
                        'Interaction detection method(s)',
                        'Publication 1st author(s)',
                        'Publication Identifier(s)',
                        'Taxid interactor A',
                        'Taxid interactor B',
                        'Biological role(s) interactor A',
                        'Biological role(s) interactor B',
                        'Experimental role(s) interactor A',
                        'Experimental role(s) interactor B',
                        'Type(s) interactor A',
                        'Type(s) interactor B',
                        'Xref(s) interactor A',
                        'Xref(s) interactor B',
                        'Interaction Xref(s)',
                        'Annotation(s) interactor A',
                        'Annotation(s) interactor B',
                        'Interaction annotation(s)',
                        'Host organism(s)',
                        'Interaction parameter(s)',
                        'Creation date',
                        'Update date',
                        'Checksum(s) interactor A',
                        'Checksum(s) interactor B',
                        'Interaction Checksum(s)',
                        'Feature(s) interactor A',
                        'Feature(s) interactor B',
                        'Stoichiometry(s) interactor A',
                        'Stoichiometry(s) interactor B',
                        'Identification method participant A',
                        'Identification method participant B',
                        'Expansion method(s)'
                        ], axis=1, inplace=True)
    # intact_cleaned = intact_df[intact_df['Interaction type(s)'].apply(lambda x: 'direct interaction' in x)]
    intact_cleaned.loc[:, 'intact_score'] = intact_cleaned.loc[: , 'Confidence value(s)'].apply(intact_score_filter)
    intact_cleaned['pair_id'] = intact_cleaned.apply(lambda row: str(tuple(sorted([row['#ID(s) interactor A'].replace('uniprotkb:', ''), row['ID(s) interactor B'].replace('uniprotkb:', '')]))), axis=1)
    
    # deduplciation
    intact_cleaned = intact_cleaned.sort_values('intact_score', ascending=False).drop_duplicates('pair_id', keep='first')
    return intact_cleaned

def read_clean_intact(path: str) -> pd.DataFrame:
    """read in the Intact dataframe and clean it using caching
    look if for the given path a file with .cleaned exists. If yes, load it.
    Else, load the file at the specified path and clean it.

    Args:
        path (str): Path to the IntAct data file
        
    Returns:
        pd.DataFrame: Cleaned IntAct dataframe
    """
    import os
    
    # Check if cleaned cache file exists
    cached_path = path + '.cleaned'
    
    if os.path.exists(cached_path):
        print(f"Loading cached cleaned data from {cached_path}")
        return pd.read_csv(cached_path, sep='\t')
    else:
        print(f"Loading and cleaning data from {path}")
        # Load the original file
        intact_df = pd.read_csv(path, sep='\t')
        
        # Clean the data using the existing function
        intact_cleaned = clean_IntAct(intact_df)
        
        # Save the cleaned data for future use
        intact_cleaned.to_csv(cached_path, sep='\t', index=False)
        print(f"Cached cleaned data saved to {cached_path}")
        
        return intact_cleaned