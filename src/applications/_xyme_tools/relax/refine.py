import numpy as np
from openff.toolkit.topology import Molecule
from openmm.app import Modeller, PDBFile, Simulation
from openmm.openmm import LangevinIntegrator
from openmm.unit import (Quantity, angstrom, kelvin, kilojoules_per_mole, nanometer, picosecond,
                         picoseconds)
from openmmforcefields.generators import SystemGenerator
from rdkit import Chem
from xyme_tools.relax.refine import fix_ligand


def local_sampling(
    pdb_file: str,
    ligand_file: str,
    fix_ligand_bool: bool = True,
    forcefield_files=["amber14-all.xml", "amber14/tip3p.xml"],
    n_cores: int = 1,
    # Langevin parameters
    temperature: float = 300.0,
    frictionCoeff: float = 1.0,
    stepSize: float = 0.002,
    dyn_steps: int = 5000,
    dyn_pdb: str = "dyn.pdb",
    # Optimization parameters
    run_opt: bool = True,
    opt_steps: int = 1000,
    energy_tolerance=10.0,
    opt_pdb: str = "opt.pdb",
):
    """
    Relax and sample complex using energy minimization followed by Langevin dynamics.
    
    Args:
        ... (previous parameters) ...
        temperature: Temperature for Langevin dynamics (Kelvin)
        frictionCoeff: Collision frequency (1/ps)
        stepSize: Integration timestep (ps)
        n_steps: Number of Langevin dynamics steps
        report_interval: Interval for saving trajectory frames
        save_trajectory: Whether to save trajectory
        trajectory_file: Output trajectory file path (DCD format)
    """
    pdb = PDBFile(pdb_file)
    print("Loading ligand params")
    if ligand_file is None:
        pass
    elif ligand_file.endswith(".sdf"):
        mol = Chem.SDMolSupplier(ligand_file, removeHs=False)[0]
        if mol is None:
            raise ValueError("Failed to load SDF file")
        # Ensure aromaticity is properly perceived
        Chem.SetAromaticity(mol)
        ligand_params = Molecule.from_rdkit(mol, allow_undefined_stereo=True)
        ligand_params.assign_partial_charges(partial_charge_method='mmff94')
    else:
        raise ValueError("Supported formats: sdf")

    print("Setting up forcefields")
    if ligand_file is None:
        system_generator = SystemGenerator(
            forcefields=forcefield_files,
            small_molecule_forcefield="gaff-2.11",
        )
    else:
        system_generator = SystemGenerator(
            forcefields=forcefield_files,
            small_molecule_forcefield="gaff-2.11",
            molecules=[ligand_params],
        )

    modeller = Modeller(pdb.topology, pdb.positions)
    if ligand_file is not None:
        ligand_omm_top = ligand_params.to_topology().to_openmm()
        ligand_positions = Quantity(
            ligand_params.conformers[0], unit=angstrom
        ).in_units_of(nanometer)
        modeller.add(ligand_omm_top, ligand_positions)

        # creating system
        print("creating system")
        system = system_generator.create_system(
            topology=modeller.topology, molecules=[ligand_params]
        )
    else:
        system = system_generator.create_system(
            topology=modeller.topology,
        )

    positions = modeller.positions.value_in_unit(nanometer)

    if fix_ligand_bool:
        system = fix_ligand(system, modeller, ligand_params, mol, positions)

    # Set up integrator with specified parameters
    integrator = LangevinIntegrator(temperature * kelvin, 
                                    frictionCoeff / picosecond, 
                                    stepSize * picoseconds)

    from openmm.openmm import Platform

    platform = Platform.getPlatformByName("CPU")
    properties = {}
    properties["Threads"] = str(n_cores)
    simulation = Simulation(
        modeller.topology,
        system,
        integrator,
        platform=platform,
        platformProperties=properties,
    )
    positions = np.array(positions).reshape(-1, 3).tolist()
    simulation.context.setPositions(positions)

    if run_opt:
        # Energy minimization
        try:
            print("Starting energy minimization...")
            simulation.minimizeEnergy(
                maxIterations=opt_steps,
                tolerance=energy_tolerance * kilojoules_per_mole / nanometer,
            )

            state = simulation.context.getState(getEnergy=True, getPositions=True)
            minimized_positions = state.getPositions()
            minimized_energy = state.getPotentialEnergy()

            print("Minimization completed.")
            print(f"Final energy: {minimized_energy.value_in_unit(kilojoules_per_mole):.2f} kJ/mol")

            # Save minimized structure
            with open(opt_pdb, "w") as f:
                PDBFile.writeFile(modeller.topology, minimized_positions, f)

        except Exception as e:
            print(f"Initial minimization failed: {e}")
            print("Try adjusting minimization parameters or checking input structure")
            raise

    # Run Langevin dynamics
    try:
        print("\nStarting Langevin dynamics...")
        simulation.context.setVelocitiesToTemperature(temperature)
        simulation._simulate(endStep=dyn_steps)


        # Get final state
        final_state = simulation.context.getState(getEnergy=True, getPositions=True)
        final_positions = final_state.getPositions()
        final_energy = final_state.getPotentialEnergy()

        print("Dynamics completed.")
        print(f"Final energy after dynamics: {final_energy.value_in_unit(kilojoules_per_mole):.2f} kJ/mol")

        # Save final structure
        with open(dyn_pdb, "w") as f:
            PDBFile.writeFile(modeller.topology, final_positions, f)
    except Exception as e:
        if "NaN" in str(e):
            print("Simulation became unstable (NaN coordinates detected)")
            print("Suggestions:")
            print("1. Reduce timestep (current: {stepSize} ps)")
            print("2. Add/increase minimization steps")
            print("3. Check input structure quality")
        else:
            print(f"Simulation failed with error: {e}")
        raise

    return final_state