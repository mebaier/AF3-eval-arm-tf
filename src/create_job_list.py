import pickle, itertools

with open('../jobs_AF/side1_uniprot_ids.pkl', 'rb') as file:
        side1 = pickle.load(file)
        file.close()
        
with open('../jobs_AF/side2_uniprot_ids.pkl', 'rb') as file:
        side2 = pickle.load(file)
        file.close()
        
jobs = []

for job in itertools.product(side1, side2):
    jobs.append(job)

with open('../jobs_AF/job_list.pkl', 'wb') as file:
    pickle.dump(jobs, file)