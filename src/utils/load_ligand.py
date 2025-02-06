from rdkit import Chem
import os
from typing import Optional

def multi_pdb_to_single_sdf(
        ligand_pdb: str,
        output_dir: Optional[str] = None
) -> None:
    """
    Convert a multi-conformer PDB file to individual SDF files, one for each conformer.
    
    This function takes a PDB file containing multiple conformers of a ligand and splits
    it into separate SDF files, with each file containing a single conformer. The output
    files are named using the original filename with a conformer index suffix.
    
    Args:
        ligand_pdb (str): Path to the input PDB file containing multiple conformers.
        output_dir (Optional[str]): Directory where output SDF files will be saved.
                                  If None, files are saved in the input file's directory.
    
    Returns:
        None
    
    Example:
        >>> multi_pdb_to_single_sdf("ligand.pdb")
        # Creates files: ligand_conf_1.sdf, ligand_conf_2.sdf, etc.
        
        >>> multi_pdb_to_single_sdf("ligand.pdb", "output_folder")
        # Creates files in output_folder/
    
    Notes:
        - Hydrogen atoms are preserved and new ones are added if missing
        - Conformer numbering starts at 1
        - Each conformer maintains its original 3D coordinates
    """
    # Extract the base filename without extension for output file naming
    base_name = os.path.splitext(os.path.basename(ligand_pdb))[0]

    # Set up output directory - use input file's directory if none specified
    if output_dir is None:
        output_dir = os.path.dirname(os.path.abspath(ligand_pdb))
    
    # Create output directory if it doesn't exist
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read PDB file and add hydrogens
    mol = Chem.MolFromPDBFile(ligand_pdb, removeHs=False)  # Keep existing hydrogens
    mol = Chem.AddHs(mol, addCoords=True)  # Add any missing hydrogens with 3D coordinates

    # Process each conformer
    for i, conf in enumerate(mol.GetConformers()):
        # Create a new molecule for this conformer
        conf_mol = Chem.Mol(mol)  # Create a copy of the original molecule
        conf_mol.RemoveAllConformers()  # Remove all conformers from the copy
        conf_mol.AddConformer(conf, assignId=True)  # Add just this conformer
        
        # Generate output filename with conformer index
        output_file = os.path.join(output_dir, f"{base_name}_conf_{i+1}.sdf")
        
        # Write the conformer to an SDF file
        writer = Chem.SDWriter(output_file)
        writer.write(conf_mol)
        writer.close()

        print(f"Saved: {output_file}")


if __name__ == "__main__":
    # Example usage
    multi_pdb_to_single_sdf("example.pdb")