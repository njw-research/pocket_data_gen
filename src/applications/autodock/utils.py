import subprocess
import numpy as np
import logging
from pathlib import Path
import glob
import os
import json
import shutil
from datetime import datetime
from src.applications.autodock.active_site import get_CA_coords
from src.utils.process_pdb_list import centre_of_mass_shift_pdb
from src.applications.autodock.allign_ligand import prepare_ligand

def prepare_enzyme(enzyme_path: str, 
                  residue_ids: list = None, 
                  allow_bad_res: bool = False,
                  save_act_site_coords: bool = False):
    """
    Prepare enzyme file for AutoDock Vina with center of mass shifting
    
    Parameters:
    -----------
    enzyme_path : str
        Path to enzyme file (without extension)
    residue_ids : list, optional
        List of active site residue IDs. If None, will look for active_site.txt
    allow_bad_res : bool, optional
        Whether to allow bad residues in enzyme preparation
    
    Returns:
    --------
    tuple
        (path to prepared PDBQT file, mean active site coordinates)
    """
    # Set up logging
    logging.basicConfig(level=logging.INFO,
                       format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)
    
    try:
        # # Create Path objects
        # enzyme_path = Path(enzyme_path)
        # enzyme_dir = enzyme_path.parent
        
        # # Create output paths
        # enzyme_pdb = f'{enzyme_path}.pdb'
        # enzyme_pdb_com = f'{enzyme_path}_com.pdb'
        # enzyme_pdbqt = f'{enzyme_path}_com.pdbqt'

        # Create Path object and handle extensions
        enzyme_path = Path(enzyme_path)
        enzyme_dir = enzyme_path.parent
        
        # Ensure path has no extension for creating new files
        enzyme_stem = enzyme_path.stem
        if enzyme_path.suffix != '.pdb':
            enzyme_pdb = str(enzyme_path) + '.pdb'
        else:
            enzyme_pdb = str(enzyme_path)
            enzyme_stem = enzyme_path.stem
            
        # Create output paths using stem
        enzyme_pdb_com = str(enzyme_dir / f"{enzyme_stem}_com.pdb")
        enzyme_pdbqt = str(enzyme_dir / f"{enzyme_stem}_com.pdbqt")
        
        
        # If residue_ids not provided, look for active_site.txt
        if residue_ids is None:
            active_site_file = enzyme_dir / 'active_site.txt'
            if active_site_file.exists():
                with open(active_site_file, 'r') as f:
                    residue_ids = f.read().strip().split()
                logger.info(f"Loaded residue IDs from active_site.txt: {residue_ids}")
            else:
                raise FileNotFoundError(
                    f"No residue_ids provided and couldn't find {active_site_file}"
                )
        
        # Shift center of mass for enzyme
        pdb_shifted = centre_of_mass_shift_pdb(enzyme_pdb, enzyme_pdb_com)
        logger.info(f"Shifted enzyme center of mass: {pdb_shifted}")

        # Get active site coordinates
        act_site = get_CA_coords(pdb_shifted, residue_ids)
        mean_act_site = np.mean(act_site, axis=0)
        logger.info(f"Active site center found at: {mean_act_site}")
        
        if save_act_site_coords:
            # Save active site coordinates
            np.savetxt(
                enzyme_dir / 'active_site_coords.txt',
                act_site,
                header=f"Active site coordinates for residues: {' '.join(residue_ids)}",
                fmt='%.3f'
            )
        
        # Prepare enzyme PDBQT
        cmd = [
            'mk_prepare_receptor.py',
            '-i', pdb_shifted,
            '--write_pdbqt', enzyme_pdbqt
        ]
        
        if allow_bad_res:
            cmd.extend(['--allow_bad_res', '--default_altloc', 'A'])
        
        subprocess.run(cmd, check=True)
        logger.info(f"Prepared enzyme PDBQT: {enzyme_pdbqt}")
        
        return enzyme_pdbqt, mean_act_site
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Error in subprocess: {e}")
        raise
    except Exception as e:
        logger.error(f"Error preparing enzyme: {str(e)}")
        raise


def prepare_vina_files(pdb_dir: str, sdf_dir: str, base_output_dir: str) -> list:
    """
    Prepare all PDB and SDF files for Vina docking
    
    Parameters:
    -----------
    pdb_dir : str
        Directory containing PDB files
    sdf_dir : str
        Directory containing SDF files
    base_output_dir : str
        Base directory for output
        
    Returns:
    --------
    list
        List of preparation information for each combination
    """
    # Set up logging
    logging.basicConfig(level=logging.INFO,
                       format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)
    
    # Find all relevant files
    pdb_files = glob.glob(os.path.join(pdb_dir, "**/*_altloc_removed.pdb"), recursive=True)
    sdf_files = glob.glob(os.path.join(sdf_dir, "**/*.sdf"), recursive=True)
    
    logger.info(f"Found {len(pdb_files)} PDB files and {len(sdf_files)} SDF files")
    
    # Create unique vina directory with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    vina_dir = Path(base_output_dir) / f"vina_prep_{timestamp}"
    vina_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Created output directory: {vina_dir}")
    
    # Store preparation information
    prep_info_list = []
    
    # Process each combination
    for sdf_file in sdf_files:
        sdf_path = Path(sdf_file)
        sdf_name = sdf_path.stem
        
        for pdb_file in pdb_files:
            pdb_path = Path(pdb_file)
            pdb_name = pdb_path.parent.name
            
            # Create output subdirectory
            output_subdir = vina_dir / sdf_name / pdb_name
            output_subdir.mkdir(parents=True, exist_ok=True)
            
            logger.info(f"\nProcessing:\nPDB: {pdb_file}\nSDF: {sdf_file}\nOutput: {output_subdir}")
            
            try:
                # Copy necessary files
                new_pdb = output_subdir / pdb_path.name
                new_sdf = output_subdir / sdf_path.name
                active_site_file = pdb_path.parent / "active_site.txt"
                
                shutil.copy2(pdb_file, new_pdb)
                shutil.copy2(sdf_file, new_sdf)
                if active_site_file.exists():
                    shutil.copy2(active_site_file, output_subdir / "active_site.txt")
                
                # Prepare enzyme
                enzyme_pdbqt, mean_act_site = prepare_enzyme(
                    enzyme_path=str(new_pdb.with_suffix('')),
                    allow_bad_res=True
                )
                
                # Prepare ligand
                ligand_pdbqt = prepare_ligand(
                    ligand_path=str(new_sdf.with_suffix('')),
                    active_site_com=mean_act_site
                )
                
                # Store preparation information
                prep_info = {
                    'pdb_name': pdb_name,
                    'sdf_name': sdf_name,
                    'output_dir': str(output_subdir),
                    'enzyme_pdbqt': enzyme_pdbqt,
                    'ligand_pdbqt': ligand_pdbqt,
                    'active_site_center': mean_act_site.tolist(),
                    'status': 'completed'
                }
                prep_info_list.append(prep_info)
                
            except Exception as e:
                logger.error(f"Error processing {pdb_name} with {sdf_name}: {str(e)}")
                prep_info = {
                    'pdb_name': pdb_name,
                    'sdf_name': sdf_name,
                    'output_dir': str(output_subdir),
                    'error': str(e),
                    'status': 'failed'
                }
                prep_info_list.append(prep_info)
                continue
    
    # Save preparation information to file
    info_file = vina_dir / "prep_info.json"
    with open(info_file, 'w') as f:
        json.dump(prep_info_list, f, indent=4)
    
    return prep_info_list



# def prepare_ligand(ligand_path, 
#                    target_com: np.array = None, 
#                    from_smiles: bool=False, 
#                    ligand_smiles: str=None):
#     """
#     Prepare ligand file for AutoDock Vina with center of mass shifting
    
#     Parameters:
#     -----------
#     ligand_path : str
#         Path to ligand file (without extension)
#     from_smiles : bool, optional
#         Whether to generate SDF from SMILES first
#     ligand_smiles : str, optional
#         SMILES string if starting from SMILES
    
#     Returns:
#     --------
#     dict
#         Paths to prepared ligand files
#     """
#     try:
#         # # Create output paths
#         # ligand_sdf = f'{ligand_path}.sdf'
#         # ligand_sdf_com = f'{ligand_path}_com.sdf'
#         # ligand_pdbqt = f'{ligand_path}_com.pdbqt'
        
#         # Create Path object and handle extensions
#         ligand_path = Path(ligand_path)
#         ligand_dir = ligand_path.parent
        
#         # Ensure path has no extension for creating new files
#         ligand_stem = ligand_path.stem
#         if ligand_path.suffix != '.sdf':
#             ligand_sdf = str(ligand_path) + '.sdf'
#         else:
#             ligand_sdf = str(ligand_path)
#             ligand_stem = ligand_path.stem
            
#         # Create output paths using stem
#         ligand_sdf_com = str(ligand_dir / f"{ligand_stem}_com.sdf")
#         ligand_pdbqt = str(ligand_dir / f"{ligand_stem}_com.pdbqt")
        
#         # Generate SDF from SMILES if requested
#         if from_smiles and ligand_smiles:
#             subprocess.run([
#                 'scrub.py',
#                 ligand_smiles,
#                 '-o', ligand_sdf
#             ], check=True)
#             print(f"Generated SDF from SMILES: {ligand_sdf}")
        
#         # Shift center of mass for ligand
#         sdf_shifted = centre_of_mass_shift_sdf(ligand_sdf, ligand_sdf_com, target_com)
#         print(f"Shifted ligand center of mass: {sdf_shifted}")
        
#         # Prepare ligand PDBQT
#         subprocess.run([
#             'mk_prepare_ligand.py',
#             '-i', sdf_shifted,
#             '-o', ligand_pdbqt
#         ], check=True)
#         print(f"Prepared ligand PDBQT: {ligand_pdbqt}")
        
#         return ligand_pdbqt
        
#     except subprocess.CalledProcessError as e:
#         print(f"Error in subprocess: {e}")
#         raise
#     except Exception as e:
#         print(f"Error preparing ligand: {str(e)}")
#         raise