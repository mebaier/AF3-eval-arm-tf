import os
import json
import csv
from datetime import datetime

results_folders = [
    "./results_AF/folds_2025_05_02_09_40",
    "./results_AF/folds_2025_05_02_09_41",
]
output_csv = "alphafold_results_summary.csv"
extracted_jobs_csv = "extracted_jobs_list.csv"
ignore_file = "ignore_list_AF_out.txt"

# Load already analyzed jobs
analyzed_jobs = set()
if os.path.exists(output_csv):
    with open(output_csv, "r") as csv_file:
        reader = csv.reader(csv_file)
        next(reader, None)
        for row in reader:
            if row:
                analyzed_jobs.add(row[0])

# Load ignore list
ignore_jobs = set()
if os.path.exists(ignore_file):
    with open(ignore_file, "r") as f:
        for line in f:
            name = line.strip()
            if name:
                ignore_jobs.add(name)

results = []
extracted_jobs = []

# Iterate each folder in the list
for results_folder in results_folders:
    if not os.path.isdir(results_folder):
        print(f"Warning: {results_folder} does not exist or is not a directory, skipping.")
        continue

    for job_folder in os.listdir(results_folder):
        job_path = os.path.join(results_folder, job_folder)

        if (not os.path.isdir(job_path)
            or job_folder in analyzed_jobs
            or job_folder in ignore_jobs):
            continue

        highest_ranking_score = -1
        best_iptm = None
        best_ptm = None

        for file_name in os.listdir(job_path):
            if file_name.endswith(".json") and "summary_confidences" in file_name:
                with open(os.path.join(job_path, file_name)) as json_file:
                    data = json.load(json_file)
                score = data.get("ranking_score", -1)
                if score > highest_ranking_score:
                    highest_ranking_score = score
                    best_iptm = data.get("iptm")
                    best_ptm = data.get("ptm")

        if best_iptm is not None and best_ptm is not None:
            results.append([job_folder, best_iptm, best_ptm])

        # Extract job details for the extracted jobs list
        job_datetime = datetime.strptime(job_folder.split("_")[1], "%Y%m%d%H%M")
        extracted_jobs.append([job_folder, job_path, 2, job_datetime.strftime("%Y-%m-%d %H:%M:%S")])

# Append to results summary CSV (write header if new)
with open(output_csv, "a", newline="") as csv_file:
    writer = csv.writer(csv_file)
    if os.stat(output_csv).st_size == 0:
        writer.writerow(["Job Name", "iPTM", "PTM"])
    writer.writerows(results)

# Write extracted jobs list CSV (overwrite each time)
with open(extracted_jobs_csv, "w", newline="") as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(["Job Name", "Folder Path", "Submitted Sequences", "Date and Time"])
    writer.writerows(extracted_jobs)

print(f"Results have been updated in {output_csv}")
print(f"Extracted jobs list has been created in {extracted_jobs_csv}")