
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def add_hydrogens(input_pdb: str, output_pdb: str, ph: float = 7.0):
    """
    Add hydrogens to a PDB structure using PDBFixer.
    
    Args:
        input_pdb (str): Input PDB file path
        output_pdb (str): Output PDB file path
        ph (float): pH for hydrogen addition
    """
    # Create PDBFixer instance
    fixer = PDBFixer(filename=input_pdb)
    
    # Find missing residues and atoms
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    
    # Add missing atoms and hydrogens
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    
    # Write the output
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    return output_pdb
