"""
This script reads in a STRING database dump and translates the STRING IDs 
of the proteins to Uniprot IDs using the Uniprot API. 
It saves the resulting translation in the OUTPUT_FILE.
"""

import re
import time
import json
import zlib
import os
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry

POLLING_INTERVAL = 10
API_URL = "https://rest.uniprot.org"
OUTPUT_FILE = "/home/markus/MPI_local/data/STRING/string_to_uniprot_mapping.txt"
INPUT_FILE = "/home/markus/MPI_local/data/STRING/9606.protein.links.full.v12.0.txt"

def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(f"HTTP Error: {response.status_code}")
        try:
            print(response.json())
        except:
            print(f"Response text: {response.text[:500]}")
        raise

def submit_id_mapping(session, from_db, to_db, ids):
    url = f"{API_URL}/idmapping/run"
    data = {"from": from_db, "to": to_db, "ids": ",".join(ids)}
    
    # # Add debugging information
    # print(f"Submitting ID mapping request to: {url}")
    # print(f"Data params: from={from_db}, to={to_db}, ids count={len(ids)}")
    # print(f"First few IDs: {','.join(ids[:5])}")
    
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    request = session.post(url, data=data, headers=headers)
    
    # # Debug response
    # print(f"Response status code: {request.status_code}")
    # try:
        # response_json = request.json()
        # print(f"Response JSON: {json.dumps(response_json)[:200]}...")
    # except:
    #     print(f"Response text: {request.text[:200]}...")
    
    check_response(request)
    return request.json()["jobId"]

def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def check_id_mapping_results_ready(session, job_id):
    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] in ("NEW", "RUNNING"):
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(j["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])

def get_batch(session, batch_response, file_format, compressed):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)

def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results

def get_id_mapping_results_link(session, job_id):
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]

def decode_results(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text


def get_xml_namespace(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""

def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)

def print_progress_batches(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}")

def get_id_mapping_results_search(session, url):
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = [str(size)]
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session.get(url)
    check_response(request)
    results = decode_results(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches(0, size, total)
    for i, batch in enumerate(get_batch(session, request, file_format, compressed), 1):
        results = combine_batches(results, batch, file_format)
        print_progress_batches(i, size, total)
    if file_format == "xml":
        return merge_xml_results(results)
    return results

def get_id_mapping_results_stream(session, url):
    if "/stream/" not in url:
        url = url.replace("/results/", "/results/stream/")
    request = session.get(url)
    check_response(request)
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    return decode_results(request, file_format, compressed)

def load_mapped_ids():
    """Load already mapped STRING IDs from the output file."""
    if os.path.exists(OUTPUT_FILE):
        with open(OUTPUT_FILE, "r") as f:
            return set(line.split("\t")[0] for line in f)
    return set()

def extract_string_ids():
    """Extract STRING IDs from the input file."""
    with open(INPUT_FILE, "r") as f:
        lines = f.readlines()[1:]  # Skip the header line
        string_ids = set()
        for line in lines:
            parts = line.strip().split()
            if len(parts) >= 2:
                string_ids.add(parts[0])
                string_ids.add(parts[1])
    return string_ids

def translate_and_save_in_batches(session, string_ids, batch_size=1000):
    """Translate STRING IDs to UniProt IDs in batches and save the results."""
    string_ids = list(string_ids)
    for i in range(0, len(string_ids), batch_size):
        batch = string_ids[i:i + batch_size]
        print(f"Processing batch {i // batch_size + 1}...")
        job_id = submit_id_mapping(session, from_db="STRING", to_db="UniProtKB", ids=batch)
        if check_id_mapping_results_ready(session, job_id):
            link = get_id_mapping_results_link(session, job_id)
            results = get_id_mapping_results_search(session, link)
            
            # Debug the results
            # print(f"Results type: {type(results)}")
            # print(f"Results preview: {str(results)[:200]}...")
            
            with open(OUTPUT_FILE, "a") as f:
                # Handle different result formats
                if isinstance(results, dict):
                    # Process successful mappings
                    if "results" in results and results["results"]:
                        for result in results["results"]:
                            f.write(f"{result['from']}	{result['to']['primaryAccession']}\n")
                    
                    # Process failed IDs
                    if "failedIds" in results and results["failedIds"]:
                        failed_count = len(results["failedIds"])
                        print(f"Found {failed_count} failed IDs")
                        
                        # Write failed IDs to a separate file for later inspection
                        failed_file = OUTPUT_FILE + ".failed.txt"
                        with open(failed_file, "w") as ff:
                            for failed_id in results["failedIds"]:
                                ff.write(f"{failed_id}\n")
                        
                        print(f"Failed IDs written to {failed_file}")
                else:
                    print(f"Warning: Unexpected results format: {type(results)}")

def main():
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))
    print('Session initialized.')

    print('Loading cached mappings')
    mapped_ids = load_mapped_ids()
    print(f'{len(mapped_ids)} cached mappings loaded')

    # Extract STRING IDs from the input file
    print('extracting STRING IDs')
    string_ids = extract_string_ids()
    print(f'STRING IDs extracted: {len(string_ids)}')

    # Filter out already cached and mapped IDs
    new_ids = string_ids - mapped_ids
    
    if not new_ids:
        print("All STRING IDs have already been processed.")
    else:
        # Translate and save new IDs in batches
        print(f'mapping {len(new_ids)} new STRING IDs')
        translate_and_save_in_batches(session, new_ids)

    print("Processing complete.")
    
if __name__ == "__main__":
    main()