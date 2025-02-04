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



def submit_pocketgen_jobs(input_dir: str, base_output_dir: str) -> list:
    """
    Submit PocketGen jobs for matching PDB and SDF files found in subdirectories
    Returns list of job information for later retrieval
    
    Parameters:
    -----------
    input_dir : str
        Base directory containing subdirectories with PDB and SDF files
    base_output_dir : str
        Base directory for output files
        
    Returns:
    --------
    list
        List of job information dictionaries
    """
    # Find all relevant files recursively
    pdb_files = glob.glob(os.path.join(input_dir, "**/*_altloc_removed_com.pdb"), recursive=True)
    sdf_files = glob.glob(os.path.join(input_dir, "**/*_aligned.sdf"), recursive=True)
    
    print(f"Found {len(pdb_files)} PDB files and {len(sdf_files)} SDF files")
    
    # Group files by their parent directory
    directory_pairs = {}
    for pdb_file in pdb_files:
        parent_dir = os.path.dirname(pdb_file)
        if parent_dir not in directory_pairs:
            directory_pairs[parent_dir] = {'pdb': [], 'sdf': []}
        directory_pairs[parent_dir]['pdb'].append(pdb_file)
        
    for sdf_file in sdf_files:
        parent_dir = os.path.dirname(sdf_file)
        if parent_dir not in directory_pairs:
            directory_pairs[parent_dir] = {'pdb': [], 'sdf': []}
        directory_pairs[parent_dir]['sdf'].append(sdf_file)
    
    # Create unique pocketgen directory
    pocketgen_dir = get_unique_pocketgen_dir(base_output_dir)
    print(f"Creating output directory: {pocketgen_dir}")
    
    # Set working directory
    set_working_dir(base_output_dir)
    
    # Store job information
    job_info_list = []
    
    # Process each directory that has both PDB and SDF files
    for directory, files in directory_pairs.items():
        if not files['pdb'] or not files['sdf']:
            print(f"Skipping {directory} - missing PDB or SDF files")
            continue
            
        rel_dir = os.path.relpath(directory, input_dir)
        
        for sdf_file in files['sdf']:
            sdf_name = Path(sdf_file).stem
            
            for pdb_file in files['pdb']:
                pdb_name = Path(pdb_file).stem
                
                # Create output subdirectory maintaining relative path structure
                output_subdir = os.path.join(
                    pocketgen_dir,
                    rel_dir,
                    f"{sdf_name}_{pdb_name}"
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