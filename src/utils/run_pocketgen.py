import os
import glob
from pathlib import Path
import json
from datetime import datetime
from xyme_tools.structure_generation import StructureGenerationInput, PocketGen
from xyme_tools.general.custom_typing import path_to_subdir
from xyme_tools.general.utils import set_working_dir


def run_pocketgen_job(pdb_file: str, output_subdir: str, sdf_file: str, wait: bool = False) -> dict:
    """
    Run a single PocketGen job with specified inputs
    Returns job information needed for retrieval
    """
    pg = PocketGen()
    
    input_data = [
        StructureGenerationInput(
            pdb_file=pdb_file,
            output_subdir=path_to_subdir(output_subdir),
            sdf_file=sdf_file
        )
    ]
    
    job_name = pg.run(input_data_list=input_data, wait=wait, backend='tamarind')
    print(f"Submitted job: {job_name}")
    
    job_info = {
        'job_name': job_name,
        'input_data': input_data,
        'pdb_file': pdb_file,
        'sdf_file': sdf_file,
        'output_dir': output_subdir
    }
    
    if wait:
        pg_output = pg.retrieve(
            input_data_list=input_data, 
            backend='tamarind',
            job_name=job_name,
            force_download=False
        )
        job_info['output'] = pg_output
        
    return job_info


def get_unique_pocketgen_dir(base_dir: str) -> str:
   """Create a unique directory name with timestamp"""
   timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
   return os.path.join(base_dir, f"pocketgen_{timestamp}")

def submit_pocketgen_jobs(pdb_dir: str, sdf_dir: str, base_output_dir: str) -> list:
   """
   Submit PocketGen jobs for all combinations of PDB and SDF files
   Returns list of job information for later retrieval
   """
   # Find all relevant files
   pdb_files = glob.glob(os.path.join(pdb_dir, "**/*_altloc_removed.pdb"), recursive=True)
   sdf_files = glob.glob(os.path.join(sdf_dir, "**/*.sdf"), recursive=True)
   
   print(f"Found {len(pdb_files)} PDB files and {len(sdf_files)} SDF files")
   
   # Create unique pocketgen directory
   pocketgen_dir = get_unique_pocketgen_dir(base_output_dir)
   print(f"Creating output directory: {pocketgen_dir}")
   
   # Set working directory
   set_working_dir(base_output_dir)
   
   # Store job information
   job_info_list = []
   
   # Submit jobs for each combination
   for sdf_file in sdf_files:
       sdf_name = Path(sdf_file).stem
       
       for pdb_file in pdb_files:
           pdb_name = Path(pdb_file).parent.name
           
           output_subdir = os.path.join(
               pocketgen_dir,
               sdf_name,
               pdb_name
           )
           os.makedirs(output_subdir, exist_ok=True)
           
           print(f"\nProcessing:\nPDB: {pdb_file}\nSDF: {sdf_file}\nOutput: {output_subdir}")
           
           try:
               job_info = run_pocketgen_job(pdb_file, output_subdir, sdf_file)
               job_info_list.append(job_info)
           except Exception as e:
               print(f"Error processing {pdb_file} with {sdf_file}: {str(e)}")
               continue
   
   # Save job information to file
   job_info_file = os.path.join(pocketgen_dir, "job_info.json")
   with open(job_info_file, 'w') as f:
       json.dump([{k: str(v) if isinstance(v, Path) else v 
                  for k, v in job.items() if k != 'input_data'} 
                 for job in job_info_list], f, indent=4)
   
   return job_info_list

def retrieve_pocketgen_jobs(job_info_list: list) -> list:
   """
   Retrieve results for previously submitted PocketGen jobs
   """
   pg = PocketGen()
   results = []
   
   for job_info in job_info_list:
       try:
           input_data = [
               StructureGenerationInput(
                   pdb_file=job_info['pdb_file'],
                   output_subdir=path_to_subdir(job_info['output_dir']),
                   sdf_file=job_info['sdf_file']
               )
           ]
           
           output = pg.retrieve(
               input_data_list=input_data,
               backend='tamarind',
               job_name=job_info['job_name'],
               force_download=False
           )
           results.append(output)
       except Exception as e:
           print(f"Error retrieving job {job_info['job_name']}: {str(e)}")
           continue
           
   return results

# def submit_pocketgen_jobs(pdb_dir: str, sdf_dir: str, base_output_dir: str) -> list:
#     """
#     Submit PocketGen jobs for all combinations of PDB and SDF files
#     Returns list of job information for later retrieval
#     """
#     # Find all relevant files
#     pdb_files = glob.glob(os.path.join(pdb_dir, "**/*_altloc_removed.pdb"), recursive=True)
#     sdf_files = glob.glob(os.path.join(sdf_dir, "**/*.sdf"), recursive=True)
    
#     print(f"Found {len(pdb_files)} PDB files and {len(sdf_files)} SDF files")
    
#     # Set working directory
#     set_working_dir(base_output_dir)
    
#     # Store job information
#     job_info_list = []
    
#     # Submit jobs for each combination
#     for sdf_file in sdf_files:
#         sdf_name = Path(sdf_file).stem
        
#         for pdb_file in pdb_files:
#             pdb_name = Path(pdb_file).parent.name
            
#             output_subdir = os.path.join(
#                 base_output_dir,
#                 "pocketgen",
#                 sdf_name,
#                 pdb_name
#             )
#             os.makedirs(output_subdir, exist_ok=True)
            
#             print(f"\nProcessing:\nPDB: {pdb_file}\nSDF: {sdf_file}\nOutput: {output_subdir}")
            
#             try:
#                 job_info = run_pocketgen_job(pdb_file, output_subdir, sdf_file)
#                 job_info_list.append(job_info)
#             except Exception as e:
#                 print(f"Error processing {pdb_file} with {sdf_file}: {str(e)}")
#                 continue
    
#     # Save job information to file
#     job_info_file = os.path.join(base_output_dir, "pocketgen", "job_info.json")
#     with open(job_info_file, 'w') as f:
#         json.dump([{k: str(v) if isinstance(v, Path) else v 
#                    for k, v in job.items() if k != 'input_data'} 
#                   for job in job_info_list], f, indent=4)
    
#     return job_info_list

# # Function to retrieve results later
# def retrieve_pocketgen_jobs(job_info_list: list) -> list:
#     """
#     Retrieve results for previously submitted PocketGen jobs
#     """
#     pg = PocketGen()
#     results = []
    
#     for job_info in job_info_list:
#         try:
#             input_data = [
#                 StructureGenerationInput(
#                     pdb_file=job_info['pdb_file'],
#                     output_subdir=path_to_subdir(job_info['output_dir']),
#                     sdf_file=job_info['sdf_file']
#                 )
#             ]
            
#             output = pg.retrieve(
#                 input_data_list=input_data,
#                 backend='tamarind',
#                 job_name=job_info['job_name'],
#                 force_download=False
#             )
#             results.append(output)
#         except Exception as e:
#             print(f"Error retrieving job {job_info['job_name']}: {str(e)}")
#             continue
            
#     return results

# Usage example:
# job_info_list = submit_pocketgen_jobs(pdb_dir, sdf_dir, base_output_dir)
# 
# # Later, to retrieve results:
# results = retrieve_pocketgen_jobs(job_info_list)