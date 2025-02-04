import json 

from src.applications.pocketgen.run_pocketgen import retrieve_pocketgen_jobs

json_file_path ='output_pipeline/pocketgen_20250202_173549/job_info.json'

with open(json_file_path, 'r') as f:
    job_info_list = json.load(f)

results = retrieve_pocketgen_jobs(job_info_list)