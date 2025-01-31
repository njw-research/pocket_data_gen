from rdkit import Chem
import numpy as np
from scipy.spatial.transform import Rotation
from pathlib import Path
import subprocess

def align_ligand(input_sdf: str, 
                output_sdf: str,
                protein_com: np.ndarray,
                active_site_com: np.ndarray) -> str:
    """
    Align ligand based on protein and active site vectors
    
    Parameters:
    -----------
    input_sdf : str
        Input SDF file path
    output_sdf : str
        Output SDF file path
    protein_com : np.ndarray
        Center of mass of protein
    active_site_com : np.ndarray
        Center of mass of active site
        
    Returns:
    --------
    str
        Path to aligned SDF file
    """
    # Read molecule
    mol = Chem.SDMolSupplier(input_sdf, removeHs=False)[0]
    if mol is None:
        raise ValueError(f"Could not read molecule from {input_sdf}")
    
    # Get conformer
    conf = mol.GetConformer()
    
    # Calculate ligand center of mass
    ligand_com = np.zeros(3)
    total_mass = 0.0
    
    for atom in mol.GetAtoms():
        mass = atom.GetMass()
        pos = conf.GetAtomPosition(atom.GetIdx())
        coord = np.array([pos.x, pos.y, pos.z])
        ligand_com += mass * coord
        total_mass += mass
    
    ligand_com /= total_mass
    
    # Calculate oxygen atoms center of mass
    oxygen_coords = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O':
            pos = conf.GetAtomPosition(atom.GetIdx())
            oxygen_coords.append(np.array([pos.x, pos.y, pos.z]))
    
    if not oxygen_coords:
        raise ValueError("No oxygen atoms found in ligand")
        
    oxygen_com = np.mean(oxygen_coords, axis=0)
    
    # Calculate vectors
    protein_vector = active_site_com - protein_com
    protein_vector = protein_vector / np.linalg.norm(protein_vector)
    
    ligand_vector = ligand_com - oxygen_com 
    ligand_vector = ligand_vector / np.linalg.norm(ligand_vector)
    
    # Calculate rotation matrix to align vectors
    rotation_matrix = get_alignment_rotation(ligand_vector, protein_vector)
    
    # Apply rotation to all atoms
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        coord = np.array([pos.x, pos.y, pos.z])
        
        # Translate to origin (relative to ligand COM)
        coord = coord - ligand_com
        
        # Apply rotation
        rotated_coord = rotation_matrix @ coord
        
        # Translate oxygen COM to active site COM
        final_coord = rotated_coord + active_site_com
        
        # Update position
        conf.SetAtomPosition(atom.GetIdx(), final_coord)
    
    # Write aligned molecule
    writer = Chem.SDWriter(output_sdf)
    writer.write(mol)
    writer.close()
    
    return output_sdf


def get_alignment_rotation(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    """
    Get rotation matrix to align vector v1 with vector v2
    
    Parameters:
    -----------
    v1 : np.ndarray
        Source vector
    v2 : np.ndarray
        Target vector
        
    Returns:
    --------
    np.ndarray
        3x3 rotation matrix
    """
    # Handle parallel or anti-parallel vectors
    if np.allclose(v1, v2):
        return np.eye(3)
    if np.allclose(v1, -v2):
        # Rotate 180° around an arbitrary perpendicular axis
        perp = np.array([1, 0, 0]) if not np.allclose(v1, [1, 0, 0]) else np.array([0, 1, 0])
        perp = perp - np.dot(perp, v1) * v1
        perp = perp / np.linalg.norm(perp)
        return Rotation.from_rotvec(np.pi * perp).as_matrix()
    
    # Calculate rotation axis and angle
    rotation_axis = np.cross(v1, v2)
    rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
    
    cos_angle = np.dot(v1, v2)
    angle = np.arccos(np.clip(cos_angle, -1.0, 1.0))
    
    # Create rotation matrix
    return Rotation.from_rotvec(angle * rotation_axis).as_matrix()

def prepare_ligand(ligand_path: str,
                    active_site_com: np.ndarray,
                    protein_com: np.ndarray = np.zeros(3, dtype=np.float32),
                    from_smiles: bool=False,
                    ligand_smiles: str=None) -> str:
    """
    Modified prepare_ligand function to align ligand with protein vector
    """
    try:
        # Create Path object and handle extensions
        ligand_path = Path(ligand_path)
        ligand_dir = ligand_path.parent
        
        # Ensure path has no extension for creating new files
        ligand_stem = ligand_path.stem
        if ligand_path.suffix != '.sdf':
            ligand_sdf = str(ligand_path) + '.sdf'
        else:
            ligand_sdf = str(ligand_path)
            ligand_stem = ligand_path.stem
        
        # Create output paths using stem
        ligand_sdf_aligned = str(ligand_dir / f"{ligand_stem}_aligned.sdf")
        ligand_pdbqt = str(ligand_dir / f"{ligand_stem}_aligned.pdbqt")
        
        # Generate SDF from SMILES if requested
        if from_smiles and ligand_smiles:
            subprocess.run([
                'scrub.py',
                ligand_smiles,
                '-o', ligand_sdf
            ], check=True)
            print(f"Generated SDF from SMILES: {ligand_sdf}")
        
        # Align ligand
        sdf_aligned = align_ligand(
            ligand_sdf,
            ligand_sdf_aligned,
            protein_com,
            active_site_com
        )
        print(f"Aligned ligand: {sdf_aligned}")
        
        # Prepare ligand PDBQT
        subprocess.run([
            'mk_prepare_ligand.py',
            '-i', sdf_aligned,
            '-o', ligand_pdbqt
        ], check=True)
        print(f"Prepared ligand PDBQT: {ligand_pdbqt}")
        
        return ligand_pdbqt
        
    except subprocess.CalledProcessError as e:
        print(f"Error in subprocess: {e}")
        raise
    except Exception as e:
        print(f"Error preparing ligand: {str(e)}")
        raise

# # Example usage
# if __name__ == "__main__":
#     # Example coordinates
#     protein_com = np.array([0, 0, 0])
#     active_site_com = np.array([10, 0, 0])
    
#     # Align and prepare ligand
#     ligand_pdbqt = modify_prepare_ligand(
#         'molecule.sdf',
#         protein_com,
#         active_site_com
#     )




# def align_ligand(input_sdf: str, 
#                 output_sdf: str,
#                 protein_com: np.ndarray,
#                 active_site_com: np.ndarray) -> str:
#     """
#     Align ligand based on protein and active site vectors
#     Aligns oxygen COM with active site COM and orients ligand accordingly
#     """
#     # Read molecule
#     mol = Chem.SDMolSupplier(input_sdf, removeHs=False)[0]
#     if mol is None:
#         raise ValueError(f"Could not read molecule from {input_sdf}")
    
#     # Get conformer
#     conf = mol.GetConformer()
    
#     # Calculate ligand center of mass
#     ligand_com = np.zeros(3)
#     total_mass = 0.0
    
#     for atom in mol.GetAtoms():
#         mass = atom.GetMass()
#         pos = conf.GetAtomPosition(atom.GetIdx())
#         coord = np.array([pos.x, pos.y, pos.z])
#         ligand_com += mass * coord
#         total_mass += mass
    
#     ligand_com /= total_mass
    
#     # Calculate oxygen atoms center of mass
#     oxygen_coords = []
#     for atom in mol.GetAtoms():
#         if atom.GetSymbol() == 'O':
#             pos = conf.GetAtomPosition(atom.GetIdx())
#             oxygen_coords.append(np.array([pos.x, pos.y, pos.z]))
    
#     if not oxygen_coords:
#         raise ValueError("No oxygen atoms found in ligand")
        
#     oxygen_com = np.mean(oxygen_coords, axis=0)
    
#     # Calculate vectors
#     protein_vector = active_site_com - protein_com
#     protein_vector = protein_vector / np.linalg.norm(protein_vector)
    
#     ligand_vector = ligand_com  - oxygen_com 
#     ligand_vector = ligand_vector / np.linalg.norm(ligand_vector)
    
#     # Calculate rotation matrix to align vectors
#     rotation_matrix = get_alignment_rotation(ligand_vector, protein_vector)
    
#     # Calculate the translation to move oxygen_com to active_site_com
#     translation_vector = active_site_com - oxygen_com
    
#     # Apply rotation and translation to all atoms
#     for atom in mol.GetAtoms():
#         pos = conf.GetAtomPosition(atom.GetIdx())
#         coord = np.array([pos.x, pos.y, pos.z])
        
#         # First rotate around oxygen_com
#         coord_rel_oxygen = coord - oxygen_com
#         rotated_coord = rotation_matrix @ coord_rel_oxygen
#         coord = rotated_coord + oxygen_com
        
#         # Then translate to align oxygen_com with active_site_com
#         final_coord = coord + translation_vector
        
#         # Update position
#         conf.SetAtomPosition(atom.GetIdx(), final_coord)
    
#     # Write aligned molecule
#     writer = Chem.SDWriter(output_sdf)
#     writer.write(mol)
#     writer.close()
    
#     return output_sdf