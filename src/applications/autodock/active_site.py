from Bio import PDB
import numpy as np
from pathlib import Path

def get_CA_coords(pdb_file: str, residue_ids: list) -> dict:
    """
    Get CA coordinates for specific residue numbers from PDB file
    
    Parameters:
    -----------
    pdb_file : str
        Path to PDB file
    residue_ids : list
        List of residue identifiers (e.g., ['S105', 'D187', 'H224'])
        
    Returns:
    --------
    dict
        Dictionary mapping residue numbers to their CA coordinates
    """
    # Extract only the numbers from residue IDs
    residue_numbers = [int(''.join(filter(str.isdigit, res))) for res in residue_ids]
    
    # Initialize coordinate dictionary
    CA_coords = {}
    
    # Parse PDB file
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    # Get first model
    model = structure[0]
    
    # Extract CA coordinates for specified residues
    for res_num in residue_numbers:
        for chain in model:
            for residue in chain:
                if residue.id[1] == res_num:  # Match residue number
                    if 'CA' in residue:
                        CA_coords[res_num] = residue['CA'].get_coord()
    
    # Extract numbers from residue IDs to maintain order
    residue_numbers = [int(''.join(filter(str.isdigit, res))) for res in residue_ids]
    
    # Create array in same order as input residue list
    coords_array = np.array([CA_coords[num] for num in residue_numbers])

    return coords_array


def load_active_site_coords(pdb_dir: str) -> np.ndarray:
    """
    Load active site coordinates from PDB and active_site.txt
    
    Parameters:
    -----------
    pdb_dir : str
        Directory containing PDB and active_site.txt files
    
    Returns:
    --------
    np.ndarray
        Array of CA coordinates for active site residues
    """
    pdb_dir = Path(pdb_dir)
    
    # Get PDB ID from directory name
    pdb_id = pdb_dir.name.split('_')[0]  # Handles cases like "1ABC_chain_A"
    
    # Define file paths
    pdb_file = pdb_dir / f"{pdb_id}_altloc_removed.pdb"
    active_site_file = pdb_dir / "active_site.txt"
    
    # Check files exist
    if not pdb_file.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    if not active_site_file.exists():
        raise FileNotFoundError(f"Active site file not found: {active_site_file}")
    
    # Read active site residues
    with open(active_site_file, 'r') as f:
        active_site_residues = f.read().strip().split()
    
    # Get coordinates
    coords = get_CA_coords(str(pdb_file), active_site_residues)
    
    return coords
