from openff.toolkit.topology import Molecule
from openmm.app import Modeller, PDBFile, Simulation
from openmm.openmm import LangevinIntegrator
from openmm.unit import (Quantity, angstrom, kelvin, kilojoules_per_mole, nanometer, picosecond,
                         picoseconds)
from openmmforcefields.generators import SystemGenerator
from rdkit import Chem
import numpy as np
from xyme_tools.relax.refine import fix_ligand


def load_ligand_params(ligand_file: str):
    """Load and process ligand parameters from an SDF file."""
    if ligand_file is None:
        return None
    
    if not ligand_file.endswith(".sdf"):
        raise ValueError("Supported formats: sdf")
    
    mol = Chem.SDMolSupplier(ligand_file, removeHs=False)[0]
    if mol is None:
        raise ValueError("Failed to load SDF file")
    
    Chem.SetAromaticity(mol)
    ligand_params = Molecule.from_rdkit(mol, allow_undefined_stereo=True)
    ligand_params.assign_partial_charges(partial_charge_method='mmff94')
    
    return ligand_params

def setup_system(pdb_file: str, ligand_params, forcefield_files: list):
    """Set up the system with forcefields and modeller."""
    pdb = PDBFile(pdb_file)
    modeller = Modeller(pdb.topology, pdb.positions)
    
    system_generator = SystemGenerator(
        forcefields=forcefield_files,
        small_molecule_forcefield="gaff-2.11",
        molecules=[ligand_params] if ligand_params else []
    )
    
    if ligand_params:
        ligand_omm_top = ligand_params.to_topology().to_openmm()
        ligand_positions = Quantity(
            ligand_params.conformers[0], unit=angstrom
        ).in_units_of(nanometer)
        modeller.add(ligand_omm_top, ligand_positions)
    
    system = system_generator.create_system(
        topology=modeller.topology,
        molecules=[ligand_params] if ligand_params else []
    )
    
    return system, modeller

def run_energy_minimization(simulation, modeller, opt_steps, energy_tolerance, opt_pdb):
    """Perform energy minimization and save the minimized structure."""
    print("Starting energy minimization...")
    simulation.minimizeEnergy(
        maxIterations=opt_steps,
        tolerance=energy_tolerance * kilojoules_per_mole / nanometer,
    )
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    minimized_positions = state.getPositions()
    
    print("Minimization completed.")
    with open(opt_pdb, "w") as f:
        PDBFile.writeFile(modeller.topology, minimized_positions, f)
    
    return state

def run_langevin_dynamics(simulation, temperature, dyn_steps, stepSize, modeller, dyn_pdb):
    """Run Langevin dynamics and save the final structure."""
    print("\nStarting Langevin dynamics...")
    simulation.context.setVelocitiesToTemperature(temperature)
    simulation._simulate(endStep=dyn_steps)
    
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    final_positions = final_state.getPositions()
    
    print("Dynamics completed.")
    with open(dyn_pdb, "w") as f:
        PDBFile.writeFile(modeller.topology, final_positions, f)
    
    return final_state

def local_sampling(
    pdb_file: str,
    ligand_file: str,
    fix_ligand_bool: bool = True,
    forcefield_files=["amber14-all.xml", "amber14/tip3p.xml"],
    n_cores: int = 1,
    temperature: float = 300.0,
    frictionCoeff: float = 1.0,
    stepSize: float = 0.002,
    dyn_steps: int = 5000,
    dyn_pdb: str = "dyn.pdb",
    run_opt: bool = True,
    opt_steps: int = 1000,
    energy_tolerance=10.0,
    opt_pdb: str = "opt.pdb",
):
    """Relax and sample complex using energy minimization followed by Langevin dynamics."""
    ligand_params = load_ligand_params(ligand_file)
    system, modeller = setup_system(pdb_file, ligand_params, forcefield_files)

    positions = modeller.positions.value_in_unit(nanometer)
    
    if fix_ligand_bool and ligand_params:
        system = fix_ligand(system, modeller, ligand_params, None, positions)
    
    integrator = LangevinIntegrator(temperature * kelvin, 
                                    frictionCoeff / picosecond, 
                                    stepSize * picoseconds)
    
    from openmm.openmm import Platform

    platform = Platform.getPlatformByName("CPU")
    properties = {"Threads": str(n_cores)}
    
    simulation = Simulation(
        modeller.topology, 
        system, 
        integrator, 
        platform=platform, 
        platformProperties=properties
    )
    positions = np.array(positions).reshape(-1, 3).tolist()
    simulation.context.setPositions(positions)
    
    if run_opt:
        run_energy_minimization(simulation, modeller, opt_steps, energy_tolerance, opt_pdb)
    
    final_state = run_langevin_dynamics(simulation, temperature, dyn_steps, stepSize, modeller, dyn_pdb)
    
    return final_state
