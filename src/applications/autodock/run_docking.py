from vina import Vina
import os
from pathlib import Path
import logging
import numpy as np
from src.applications.autodock.active_site import get_CA_coords

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
        docking_results = v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        
        # Save docking poses
        v.write_poses(str(docked_poses_path), n_poses=n_poses_to_save, overwrite=True)
        logger.info(f'Saved {n_poses_to_save} docked poses to: {docked_poses_path}')
        
        # # Get all docking scores
        # docking_scores = [result.affinity for result in docking_results]
        
        # # Log top scores
        # logger.info(f'\nTop {n_poses_to_save} docking scores:')
        # for i, score in enumerate(docking_scores[:n_poses_to_save]):
        #     logger.info(f'Pose {i+1}: {score:.3f} kcal/mol')
        
        # # Create results dictionary
        # results = {
        #     'initial_score': initial_energy[0],
        #     'minimized_score': minimized_energy[0],
        #     'docking_scores': docking_scores,
        #     'center': center,
        #     'box_size': box_size,
        #     'minimized_pose': str(minimized_pose_path),
        #     'docked_poses': str(docked_poses_path),
        #     'log_file': str(output_dir / 'docking.log')
        # }
        
        # logger.info("Docking completed successfully")
        # return results

        return docking_results
        
    except Exception as e:
        logger.error(f"Error during docking: {str(e)}", exc_info=True)
        raise

# def run_vina_docking(enzyme_path: str, 
#                     ligand_path: str, 
#                     output_dir: str,
#                     center=[0.0, 0.0, 0.0], 
#                     box_size=[50, 50, 50],
#                     exhaustiveness=32,
#                     n_poses=20,
#                     n_poses_to_save=5):
#     """
#     Run AutoDock Vina docking procedure
    
#     Parameters:
#     -----------
#     enzyme_path : str
#         Path to enzyme PDBQT file
#     ligand_path : str
#         Path to ligand PDBQT file
#     output_dir : str
#         Directory for output files
#     center : list
#         Center coordinates [x, y, z] for the search box
#     box_size : list
#         Size of the search box [x, y, z] in Angstroms
#     exhaustiveness : int
#         Exhaustiveness of the global search
#     n_poses : int
#         Number of poses to generate
#     n_poses_to_save : int
#         Number of top poses to save
        
#     Returns:
#     --------
#     dict
#         Dictionary containing energy scores and paths to output files
#     """
#     try:
#         # Create output directory if it doesn't exist
#         output_dir = Path(output_dir)
#         output_dir.mkdir(parents=True, exist_ok=True)
        
#         # Set up logger
#         logger = setup_logger(output_dir)
        
#         # Create output paths
#         minimized_pose_path = output_dir / 'minimized_pose.pdbqt'
#         docked_poses_path = output_dir / 'docked_poses.pdbqt'
        
#         # Initialize Vina
#         v = Vina(sf_name='vina')
        
#         # Set receptor and ligand
#         logger.info(f"Loading enzyme from: {enzyme_path}")
#         v.set_receptor(str(enzyme_path))
        
#         logger.info(f"Loading ligand from: {ligand_path}")
#         v.set_ligand_from_file(str(ligand_path))
        
#         # Compute Vina maps
#         logger.info(f"Computing Vina maps (center={center}, box_size={box_size})")
#         v.compute_vina_maps(center=center, box_size=box_size)
        
#         # Score initial pose
#         initial_energy = v.score()
#         logger.info(f'Score before minimization: {initial_energy[0]:.3f} (kcal/mol)')
        
#         # Local optimization
#         minimized_energy = v.optimize()
#         logger.info(f'Score after minimization: {minimized_energy[0]:.3f} (kcal/mol)')
        
#         # Save minimized pose
#         v.write_pose(str(minimized_pose_path), overwrite=True)
#         logger.info(f'Saved minimized pose to: {minimized_pose_path}')
        
#         # Run docking
#         logger.info(f'Starting docking (exhaustiveness={exhaustiveness}, n_poses={n_poses})')
#         docking_results = v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
        
#         # Save docking poses
#         v.write_poses(str(docked_poses_path), n_poses=n_poses_to_save, overwrite=True)
#         logger.info(f'Saved {n_poses_to_save} docked poses to: {docked_poses_path}')
        
#         # Get all docking scores
#         docking_scores = [result.affinity for result in docking_results]
        
#         # Log top scores
#         logger.info(f'\nTop {n_poses_to_save} docking scores:')
#         for i, score in enumerate(docking_scores[:n_poses_to_save]):
#             logger.info(f'Pose {i+1}: {score:.3f} kcal/mol')
        
#         # Create results dictionary
#         results = {
#             'initial_score': initial_energy[0],
#             'minimized_score': minimized_energy[0],
#             'docking_scores': docking_scores,
#             'minimized_pose': str(minimized_pose_path),
#             'docked_poses': str(docked_poses_path),
#             'log_file': str(output_dir / 'docking.log')
#         }
        
#         logger.info(results)
#         logger.info("Docking completed successfully")
#         return results

        
#     except Exception as e:
#         logger.error(f"Error during docking: {str(e)}", exc_info=True)
#         raise

