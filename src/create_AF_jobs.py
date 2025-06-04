import requests
import json
import pickle
import os
import time  # Added to generate batch ID

CACHE_FILE = '../cache/uniprot_sequences_cache.json'

def load_cache():
    """Load the cache of Uniprot sequences from a file."""
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE, 'r') as f:
            return json.load(f)
    return {}

def save_cache(cache):
    """Save the cache of Uniprot sequences to a file."""
    with open(CACHE_FILE, 'w') as f:
        json.dump(cache, f, indent=2)

def fetch_uniprot_sequence(uniprot_id, cache):
    """Fetches the sequence of a given Uniprot ID using the Uniprot API, using cache if available."""
    if uniprot_id in cache:
        return cache[uniprot_id]  # Fetch from cache if available
    
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    
    if response.status_code == 200:
        sequence = ''.join(response.text.splitlines()[1:])  # Remove FASTA header
        cache[uniprot_id] = sequence  # Cache the sequence
        return sequence
    else:
        raise Exception(f"Failed to fetch data for Uniprot ID {uniprot_id}")

def create_alphafold_job(uniprot_ids, cache):
    """Creates an Alphafold job entry for a list of Uniprot IDs."""
    job_name = "_".join(uniprot_ids)  # Name the job with Uniprot IDs
    sequences = []

    for uniprot_id in uniprot_ids:
        sequence = fetch_uniprot_sequence(uniprot_id, cache)
        sequences.append({
            "proteinChain": {
                "sequence": sequence,
                "count": 1
            }
        })

    job = {
        "name": job_name,
        "modelSeeds": [],
        "sequences": sequences,
        "dialect": "alphafoldserver",
        "version": 1
    }
    return job

def load_existing_jobs():
    """Load existing jobs from previously created JSON files."""
    existing_jobs = []
    jobs_folder = "../jobs_AF"
    if os.path.exists(jobs_folder):
        for file_name in os.listdir(jobs_folder):
            if file_name.endswith(".json"):
                with open(os.path.join(jobs_folder, file_name), 'r') as f:
                    existing_jobs.extend(json.load(f))
    return existing_jobs

def job_exists(new_job, existing_jobs):
    """Check if a job already exists in the list of existing jobs."""
    for job in existing_jobs:
        if job["sequences"] == new_job["sequences"] and job["modelSeeds"] == new_job["modelSeeds"] and job["dialect"] == new_job["dialect"] and job["version"] == new_job["version"]:
            return True
    return False

def create_jobs(input_list):
    """Creates multiple Alphafold jobs and saves them into JSON files with a maximum of 30 jobs per file."""
    all_jobs = []
    ignored_uniprot_ids = []  # List to store ignored Uniprot IDs
    existing_jobs = load_existing_jobs()
    cache = load_cache()  # Load cache once at the start

    for uniprot_ids in input_list:
        try:
            job = create_alphafold_job(uniprot_ids, cache)
            # Check if any sequence in the job contains "X"
            if any("X" in seq["proteinChain"]["sequence"] for seq in job["sequences"]):
                ignored_uniprot_ids.extend(uniprot_ids)  # Add to ignored list
                continue  # Skip this job
            if not job_exists(job, existing_jobs):
                all_jobs.append(job)
        except Exception as e:
            print(f"Error processing Uniprot IDs {uniprot_ids}: {e}")
            ignored_uniprot_ids.extend(uniprot_ids)  # Add to ignored list in case of errors

    # Save the updated cache after processing all jobs
    save_cache(cache)

    # Generate a unique batch ID using the current timestamp
    batch_id = int(time.time())

    # Write ignored Uniprot IDs to a separate file
    ignored_file = f"../jobs_AF/ignored_batch{batch_id}.txt"
    with open(ignored_file, 'w') as f:
        f.write("# Ignored Uniprot IDs\n")
        f.write("# These sequences contain 'X' or encountered errors during processing\n")
        for uniprot_id in ignored_uniprot_ids:
            f.write(f"{uniprot_id}\n")
    print(f"Ignored Uniprot IDs written to {ignored_file}")

    # Split jobs into chunks of 30
    chunk_size = 30
    for i in range(0, len(all_jobs), chunk_size):
        chunk = all_jobs[i:i + chunk_size]
        filename = f"../jobs_AF/alphafold_jobs_batch{batch_id}_{i//chunk_size + 1}.json"
        with open(filename, 'w') as f:
            json.dump(chunk, f, indent=2)
        print(f"Created {filename}")

def main():
    with open('../jobs_AF/job_list.pkl', 'rb') as file:
        job_list = pickle.load(file)
    print("number of jobs: ", len(job_list))
    create_jobs(job_list)


if __name__ == "__main__":
    main()