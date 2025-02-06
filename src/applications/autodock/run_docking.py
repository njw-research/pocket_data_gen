from vina import Vina
from pathlib import Path
import logging
import numpy as np
import glob
import os 
import json
from datetime import datetime
from typing import Optional
from src.applications.autodock.active_site import get_CA_coords
from xyme_tools.structure_prediction.tool_classes import AutoDock
from xyme_tools.general.utils import set_working_dir
from xyme_tools.structure_prediction import (
    StructurePredictionInput
)

def setup_logger(output_dir):
    """Set up logger for docking results"""
    log_path = Path(output_dir) / 'docking.log'
    
    # Create logger
    logger = logging.getLogger('vina_docking')
    logger.setLevel(logging.INFO)
    
    # Create handlers
    file_handler = logging.FileHandler(log_path)
    console_handler = logging.StreamHandler()
    
    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    
    # Add handlers to logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger


def setup_docking(input_dir: str):
    try:
        input_dir = Path(input_dir)
        output_dir = input_dir / 'docking_results'
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Find required files
        ligand_file = list(input_dir.glob("*_aligned.sdf"))
        enzyme_file = list(input_dir.glob("*_altloc_removed_com.pdb"))
        active_site_file = input_dir / "active_site.txt"
        
        if not ligand_file:
            raise FileNotFoundError("Could not find ligand file (*aligned*.pdb or *aligned*.sdf)")
        if not enzyme_file:
            raise FileNotFoundError("Could not find enzyme file (*_altloc_removed*.pdb or *_altloc_removed*.sdf)")
        
        # Convert Path objects to strings
        ligand_path = str(f"./{ligand_file[0]}")
        enzyme_path = str(f"./{enzyme_file[0]}")

        # Read active site coordinates
        with open(active_site_file, 'r') as f:
            residues = f.read().strip().split()
            
        # enzyme_pdb = enzyme_path.with_suffix('.pdb')
        active_site_coords = get_CA_coords(enzyme_path, residues)
        center = np.mean(active_site_coords, axis=0, dtype=float)
        
        return enzyme_path, ligand_path, center
        
    except Exception as e:
        print(f"Error in setup_docking: {str(e)}")
        print(f"Working directory: {Path.cwd()}")
        raise


def run_autodock_job(pdb_file: str, 
                     output_subdir: str, 
                     sdf_file: str, 
                     act_site: Optional[np.ndarray],  
                     wait: bool = False) -> dict:
    """
    Run a single PocketGen job with specified inputs
    Returns job information needed for retrieval
    """
    
    # pdb_file, sdf_file, act_site = setup_docking('examples/aligned_pocket_test/triacylglycerol/1AQL_chain_A')

    ad = AutoDock()

    input_data = [
        StructurePredictionInput(
            sequence = 'A', # Does not matter but needs to be here 
            pdb_file = pdb_file,
            sdf_file = sdf_file,
            output_subdir = output_subdir,
            box_xyz = act_site,
            box_dims = np.array([40.0, 40.0, 40.0]),
        ),
    ]

    job_name = ad.run(input_data_list=input_data, wait=False, backend='tamarind')

    print(f"Submitted job: {job_name}")
    
    job_info = {
        'job_name': job_name,
        'input_data': input_data,
        'pdb_file': pdb_file,
        'sdf_file': sdf_file,
        'output_dir': output_subdir
    }
    
    if wait:
        ad_output = ad.retrieve(
            input_data_list=input_data, 
            backend='tamarind',
            job_name=job_name,
            force_download=False
        )
        job_info['output'] = ad_output
        
    return job_info



def get_unique_pocketgen_dir(base_dir: str) -> str:
   """Create a unique directory name with timestamp"""
   timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
   return os.path.join(base_dir, f"pocketgen_{timestamp}")




def submit_autodock_jobs(input_dir: str, base_output_dir: str) -> list:
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
    glob.glob(os.path.join(input_dir, "**/active_site.txt"), recursive=True)

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
                    job_info = run_autodock_job(pdb_file, output_subdir, sdf_file)
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





def run_vina_docking(input_dir: str,
                    box_size=[50, 50, 50],
                    exhaustiveness=32,
                    n_poses=20,
                    n_poses_to_save=5):
    """
    Run AutoDock Vina docking procedure
    
    Parameters:
    -----------
    input_dir : str
        Directory containing prepared files:
        - *_aligned.pdbqt (ligand)
        - *_altloc_removed_com.pdbqt (enzyme)
        - active_site.txt
    box_size : list
        Size of the search box [x, y, z] in Angstroms
    exhaustiveness : int
        Exhaustiveness of the global search
    n_poses : int
        Number of poses to generate
    n_poses_to_save : int
        Number of top poses to save
        
    Returns:
    --------
    dict
        Dictionary containing energy scores and paths to output files
    """
    try:
        # Convert input_dir to Path and create output subdirectory
        input_dir = Path(input_dir)
        output_dir = input_dir / 'docking_results'
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set up logger
        logger = setup_logger(output_dir)
        
        # Find required files
        ligand_file = list(input_dir.glob("*_aligned.pdbqt"))
        enzyme_file = list(input_dir.glob("*_altloc_removed_com.pdbqt"))
        active_site_file = input_dir / "active_site.txt"
        
        if not ligand_file or not enzyme_file:
            raise FileNotFoundError("Could not find ligand or enzyme PDBQT files")
        
        ligand_path = ligand_file[0]
        enzyme_path = enzyme_file[0]
        scores_path = output_dir / 'docking_scores.json'

        
        # Read active site coordinates
        with open(active_site_file, 'r') as f:
            residues = f.read().strip().split()
            
        # Get active site coordinates from PDB file
        enzyme_pdb = enzyme_path.with_suffix('.pdb')
        active_site_coords = get_CA_coords(str(enzyme_pdb), residues)
        center = np.mean(active_site_coords, axis=0).tolist()
        
        # Create output paths
        minimized_pose_path = output_dir / 'minimized_pose.pdbqt'
        docked_poses_path = output_dir / 'docked_poses.pdbqt'
        
        # Initialize Vina
        v = Vina(sf_name='vina')
        
        # Set receptor and ligand
        logger.info(f"Loading enzyme from: {enzyme_path}")
        v.set_receptor(str(enzyme_path))
        
        logger.info(f"Loading ligand from: {ligand_path}")
        v.set_ligand_from_file(str(ligand_path))
        
        # Compute Vina maps
        logger.info(f"Computing Vina maps (center={center}, box_size={box_size})")
        v.compute_vina_maps(center=center, box_size=box_size)
        
        # Score initial pose
        initial_energy = v.score()
        logger.info(f'Score before minimization: {initial_energy[0]:.3f} (kcal/mol)')
        
        # Local optimization
        minimized_energy = v.optimize()
        logger.info(f'Score after minimization: {minimized_energy[0]:.3f} (kcal/mol)')
        
        # Save minimized pose
        v.write_pose(str(minimized_pose_path), overwrite=True)
        logger.info(f'Saved minimized pose to: {minimized_pose_path}')
        
        # Run docking
        logger.info(f'Starting docking (exhaustiveness={exhaustiveness}, n_poses={n_poses})')
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        
        # Save docking poses
        v.write_poses(str(docked_poses_path), n_poses=n_poses_to_save, overwrite=True)
        logger.info(f'Saved {n_poses_to_save} docked poses to: {docked_poses_path}')
        
        # Binding affintiy 
        affinity = v.energies().T[0]
        valid_poses = len(affinity)

        # Create results dictionary
        results = {
            'docking': {
                'affinity (kcal/mol)': [float(x) for x in affinity],
                'n_successful_poses': valid_poses
            },
            'parameters': {
                'center': center,
                'box_size': box_size,
                'exhaustiveness': exhaustiveness,
                'n_poses_requested': n_poses,
                'n_poses_to_save': n_poses_to_save
            },
            'files': {
                'minimized_pose': str(minimized_pose_path),
                'docked_poses': str(docked_poses_path),
                'log_file': str(output_dir / 'docking.log')
            }
        }

        # Save results to JSON
        with open(scores_path, 'w') as f:
            json.dump(results, f, indent=4)
        logger.info(f'Saved scores and energies to: {scores_path}')

        logger.info("Docking completed successfully")
        
    except Exception as e:
        logger.error(f"Error during docking: {str(e)}", exc_info=True)
        raise


