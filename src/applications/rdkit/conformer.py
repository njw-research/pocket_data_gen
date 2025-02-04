from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def generate_conformers(mol, n_conformers=50, random_seed=42):
    """
    Generate conformers for a molecule using ETKDG method
    
    Parameters:
    -----------
    mol : rdkit.Chem.rdchem.Mol
        Input molecule
    n_conformers : int
        Number of conformers to generate
    random_seed : int
        Random seed for reproducibility
    
    Returns:
    --------
    rdkit.Chem.rdchem.Mol
        Molecule with generated conformers
    """
    # Clean up the molecule
    mol = Chem.AddHs(mol)  # Add hydrogens
    
    # Set up ETKDG parameters
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    params.numThreads = 0  # Use all available CPUs
    params.useSmallRingTorsions = True
    params.useBasicKnowledge = True
    params.enforceChirality = True
    
    # Generate conformers
    AllChem.EmbedMultipleConfs(
        mol,
        numConfs=n_conformers,
        params=params
    )
    
    # Optimize all conformers using MMFF94s
    for conf_id in range(mol.GetNumConformers()):
        AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
    
    return mol

def prune_conformers(mol, rmsd_threshold=0.5, energy_window=10.0):
    """
    Prune conformers based on RMSD and energy
    
    Parameters:
    -----------
    mol : rdkit.Chem.rdchem.Mol
        Molecule with conformers
    rmsd_threshold : float
        RMSD threshold for considering conformers as unique
    energy_window : float
        Energy window (kcal/mol) above global minimum to keep conformers
    
    Returns:
    --------
    rdkit.Chem.rdchem.Mol
        Molecule with pruned conformers
    """
    # Calculate MMFF94s energies for all conformers
    energies = []
    for conf_id in range(mol.GetNumConformers()):
        mp = AllChem.MMFFGetMoleculeProperties(mol)
        ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=conf_id)
        energy = ff.CalcEnergy()
        energies.append(energy)
    
    # Convert to numpy array and calculate relative energies
    energies = np.array(energies)
    rel_energies = energies - np.min(energies)
    
    # Get indices of conformers within energy window
    energy_indices = np.where(rel_energies <= energy_window)[0]
    
    # Calculate RMSD matrix for remaining conformers
    n_confs = len(energy_indices)
    rmsd_matrix = np.zeros((n_confs, n_confs))
    
    for i in range(n_confs):
        for j in range(i+1, n_confs):
            rmsd = AllChem.GetBestRMS(mol, mol,
                                    energy_indices[i],
                                    energy_indices[j])
            rmsd_matrix[i,j] = rmsd_matrix[j,i] = rmsd
    
    # Cluster conformers based on RMSD
    keep_indices = []
    for i in range(n_confs):
        unique = True
        for j in keep_indices:
            if rmsd_matrix[i,j] < rmsd_threshold:
                unique = False
                break
        if unique:
            keep_indices.append(i)
    
    # Create new molecule with selected conformers
    new_mol = Chem.Mol(mol)
    new_mol.RemoveAllConformers()
    for idx in keep_indices:
        conf = mol.GetConformer(energy_indices[idx])
        new_mol.AddConformer(conf, assignId=True)
    
    return new_mol

def save_conformers_to_sdf(mol, output_file):
    """
    Save conformers to an SDF file
    
    Parameters:
    -----------
    mol : rdkit.Chem.rdchem.Mol
        Molecule with conformers
    output_file : str
        Output SDF filename
    """
    writer = Chem.SDWriter(output_file)
    for conf_id in range(mol.GetNumConformers()):
        writer.write(mol, confId=conf_id)
    writer.close()