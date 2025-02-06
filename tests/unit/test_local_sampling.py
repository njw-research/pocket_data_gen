import pytest
import tempfile
from pathlib import Path
import numpy as np
from openmm.unit import nanometer, kilojoules_per_mole

from src.applications._xyme_tools.local_sampling.langevin import (
    local_sampling,
    load_ligand_params,
    setup_system
)

# Test data paths
TEST_DATA_DIR = Path("tests/data/openmm")
TEST_PDB = TEST_DATA_DIR / "enzyme.pdb"
TEST_LIGAND = TEST_DATA_DIR / "ligand.sdf"

@pytest.fixture
def temp_dir():
    """Create temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield tmp_dir

@pytest.fixture
def test_files():
    """Ensure test files exist."""
    assert TEST_PDB.exists(), f"Test PDB file not found: {TEST_PDB}"
    assert TEST_LIGAND.exists(), f"Test ligand file not found: {TEST_LIGAND}"
    return TEST_PDB, TEST_LIGAND

def test_load_ligand_params(test_files):
    """Test ligand parameter loading."""
    _, ligand_file = test_files
    
    # Test successful loading
    ligand_params = load_ligand_params(str(ligand_file))
    assert ligand_params is not None
    assert hasattr(ligand_params, 'partial_charges')
    
    # Test wrong file format
    with pytest.raises(ValueError, match="Supported formats: sdf"):
        load_ligand_params("wrong_file.xyz")
    

def test_setup_system(test_files):
    """Test system setup."""
    pdb_file, ligand_file = test_files
    
    # Load ligand params first
    ligand_params = load_ligand_params(str(ligand_file))
    
    # Test system setup with ligand
    system, modeller = setup_system(
        str(pdb_file),
        ligand_params,
        ["amber14-all.xml", "amber14/tip3p.xml"]
    )
    
    assert system is not None
    assert modeller is not None
    assert system.getNumParticles() > 0
    
    # Test system setup without ligand
    system_no_ligand, modeller_no_ligand = setup_system(
        str(pdb_file),
        None,
        ["amber14-all.xml", "amber14/tip3p.xml"]
    )
    
    assert system_no_ligand is not None
    assert system_no_ligand.getNumParticles() < system.getNumParticles()

def test_full_sampling(test_files, temp_dir):
    """Test full sampling workflow."""
    pdb_file, ligand_file = test_files
    
    # Setup temporary output files
    opt_pdb = Path(temp_dir) / "opt.pdb"
    dyn_pdb = Path(temp_dir) / "dyn.pdb"
    
    # Run sampling with minimal steps for testing
    final_state = local_sampling(
        pdb_file=str(pdb_file),
        ligand_file=str(ligand_file),
        fix_ligand_bool=True,
        dyn_steps=10,  # Small number for testing
        dyn_pdb=str(dyn_pdb),
        temperature=300.0,
        run_opt=True,
        opt_steps=10,  # Small number for testing
        energy_tolerance=10.0,
        opt_pdb=str(opt_pdb)
    )
    
    # Check outputs
    assert opt_pdb.exists(), "Optimization PDB not created"
    assert dyn_pdb.exists(), "Dynamics PDB not created"
    assert final_state is not None
    
    # Check energy is reasonable (not NaN or extremely high)
    final_energy = final_state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    assert not np.isnan(final_energy), "Energy is NaN"
    assert abs(final_energy) < 1e6, "Energy is unreasonably high"
    
    # Check positions are reasonable
    positions = final_state.getPositions(asNumpy=True).value_in_unit(nanometer)
    assert not np.any(np.isnan(positions)), "Found NaN in positions"
    assert np.all(np.abs(positions) < 100), "Positions are unreasonably large"

def test_minimization_convergence(test_files, temp_dir):
    """Test energy minimization convergence."""
    pdb_file, ligand_file = test_files
    opt_pdb = Path(temp_dir) / "opt.pdb"
    
    final_state = local_sampling(
        pdb_file=str(pdb_file),
        ligand_file=str(ligand_file),
        run_opt=True,
        opt_steps=100,
        energy_tolerance=1.0,
        opt_pdb=str(opt_pdb),
        dyn_steps=0  # Skip dynamics
    )
    
    # Check if energy decreased
    initial_energy = final_state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    assert initial_energy < 1e6, "Minimization did not converge to reasonable energy"

if __name__ == "__main__":
    pytest.main([__file__])