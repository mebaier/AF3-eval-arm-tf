import requests
import os
import time

def translate_string_to_uniprot(string_ids):
    url = "https://rest.uniprot.org/idmapping/run"
    
    # Prepare the payload for the API request
    payload = {
        "from": "STRING",
        "to": "UniProtKB",
        "ids": ",".join(string_ids)
    }
    
    print(payload)

    # Send the request to the UniProt API
    response = requests.post(url, data=payload)
    if response.status_code != 200:
        raise Exception(f"Failed to initiate ID mapping: {response.text}")

    # Extract the job ID from the response
    job_id = response.json().get("jobId")
    if not job_id:
        raise Exception("No job ID returned from UniProt API")

    # Check the status of the job
    status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    while True:
        status_response = requests.get(status_url)
        if status_response.status_code != 200:
            raise Exception(f"Failed to check job status: {status_response.text}")

        status = status_response.json().get("status")
        print(status)
        if status == "FINISHED":
            break
        elif status == "FAILED":
            raise Exception("ID mapping job failed")
        time.sleep(1)

    # Retrieve the results
    results_url = f"https://rest.uniprot.org/idmapping/results/{job_id}"
    results_response = requests.get(results_url)
    if results_response.status_code != 200:
        raise Exception(f"Failed to retrieve results: {results_response.text}")

    # Parse the results
    results = results_response.json().get("results")
    if not results:
        raise Exception("No results returned from UniProt API")

    # Map STRING IDs to UniProt IDs
    string_to_uniprot = {result["from"]: result["to"]["primaryAccession"] for result in results}
    return string_to_uniprot

if __name__ == "__main__":
    input_file = "../9606.protein.links.full.v12.0.txt"
    cache_file = "string_ids_cache.txt"
    output_file = "string_to_uniprot_mapping.txt"

    # Check if cache file exists
    if os.path.exists(cache_file):
        print(f"Cache file {cache_file} found. Loading STRING IDs from cache.")
        with open(cache_file, "r") as f:
            string_ids = set(line.strip() for line in f)
    else:
        print(f"Cache file {cache_file} not found. Extracting STRING IDs from input file.")
        # Read STRING IDs from the input file, ignoring the first line
        with open(input_file, "r") as f:
            lines = f.readlines()[1:]  # Skip the header line
            string_ids = set()  # Use a set to remove duplicates
            for line in lines:
                parts = line.strip().split()  # Split the line into parts
                if len(parts) >= 2:
                    string_ids.add(parts[0])  # Add the first ID
                    string_ids.add(parts[1])  # Add the second ID

        # Write the extracted IDs to the cache file
        with open(cache_file, "w") as f:
            for string_id in string_ids:
                f.write(f"{string_id}\n")

    if os.path.exists(output_file):
        print(f"Output file {output_file} found. Loading already mapped IDs.")
        with open(output_file, "r") as f:
            mapped_ids = set(line.split("\t")[0] for line in f)
    else:
        print(f"Output file {output_file} not found. Starting fresh.")
        mapped_ids = set()

    # Filter out already mapped IDs
    unmapped_ids = [string_id for string_id in string_ids if string_id not in mapped_ids]

    if not unmapped_ids:
        print("All IDs have already been mapped. No new mapping required.")
    else:
        # Translate unmapped STRING IDs to UniProt IDs
        mapping = translate_string_to_uniprot(unmapped_ids[:10])

        # Append the new mappings to the output file
        with open(output_file, "a") as f:
            for string_id, uniprot_id in mapping.items():
                f.write(f"{string_id}\t{uniprot_id}\n")

        print(f"Mapping completed for unmapped IDs. Results appended to {output_file}.")