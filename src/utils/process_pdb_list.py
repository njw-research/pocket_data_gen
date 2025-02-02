import os
import pandas as pd
import numpy as np
import ast
import logging
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO
from xyme_tools.structure_handling.utils import download_pdb_to_local, ResidueSelect, remove_altloc #, filter_to_pdb_hetatm_types

def filter_to_pdb_hetatm_types(
    pdb_f: str, hetatms: list[dict], chain:str = None, f_out: str = None
) -> str:

    logging.getLogger(__name__)

    structure = PDBParser(QUIET=True).get_structure("X", pdb_f)

    keep_hetams = [
        res
        for res in structure.get_residues()
        if any(
            [
                (res.parent.id == hetatm_res["CHAIN"])
                & (res.id[1] == hetatm_res["RESID"])
                & (res.resname == hetatm_res["RESNAME"].upper())
                for hetatm_res in hetatms
            ]
        )
    ]

    non_hetatms = [
        res for res in structure.get_residues() if res.get_full_id()[3][0] == " "
    ]

    keep_residues = set(keep_hetams + non_hetatms)

    io = PDBIO()
    io.set_structure(structure)
    if not f_out:
        f_out = pdb_f.replace(".pdb", "_filtered.pdb")

    out_dir = os.path.dirname(f_out)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    io.save(f_out, ResidueSelect(keep_residues))

    return f_out


def extract_chain(pdb_file: str, chain_id: str, output_file: str = None) -> str:
   """
   Extract a single chain from a PDB file.

   Args:
       pdb_file (str): Input PDB file path
       chain_id (str): Chain identifier to extract (e.g., 'A')
       output_file (str, optional): Output file path. If None, creates filename based on input.

   Returns:
       str: Path to output file containing extracted chain
   """
   # Create default output filename if none provided
   if output_file is None:
       output_file = pdb_file.replace('.pdb', f'_chain_{chain_id}.pdb')

   # Read and process the PDB file
   with open(pdb_file, 'r') as f:
       lines = f.readlines()

   # Keep lines that are either not ATOM/HETATM records, or match our chain
   kept_lines = []
   for line in lines:
       if line.startswith(('ATOM', 'HETATM')):
           if line[21] == chain_id:  # Chain ID is in column 22
               kept_lines.append(line)
       else:
           kept_lines.append(line)

   # Write the output file
   with open(output_file, 'w') as f:
       f.writelines(kept_lines)

   return output_file


def centre_of_mass_shift_pdb(input_pdb, output_pdb):
    """
    Shift PDB coordinates to center of mass using BioPython
    
    Parameters:
    -----------
    input_pdb : str
        Input PDB file path
    output_pdb : str
        Output PDB file path
    """
    from Bio import PDB
    
    # Set up parser
    parser = PDB.PDBParser(QUIET=True)  # Added QUIET=True to suppress warnings
    structure = parser.get_structure('protein', input_pdb)
    
    # Calculate COM
    com = structure.center_of_mass()
    
    # Apply shift
    for atom in structure.get_atoms():
        atom.coord = atom.coord - com
    
    # Save shifted structure
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
    
    return output_pdb



def centre_of_mass_shift_sdf(input_sdf: str, 
                             output_sdf: str, 
                             target_com: np.array = None):
    """
    Shift SDF coordinates to center of mass using RDKit
    
    Parameters:
    -----------
    input_sdf : str
        Input SDF file path
    output_sdf : str
        Output SDF file path
    """
    from rdkit import Chem
    import numpy as np
    
    # Read molecule from SDF
    mol = Chem.SDMolSupplier(input_sdf, removeHs=False)[0]
    if mol is None:
        raise ValueError(f"Could not read molecule from {input_sdf}")
    
    # Get conformer
    conf = mol.GetConformer()
    
    # Calculate COM
    com = np.zeros(3)
    total_mass = 0.0
    
    for atom in mol.GetAtoms():
        mass = atom.GetMass()
        pos = conf.GetAtomPosition(atom.GetIdx())
        coord = np.array([pos.x, pos.y, pos.z])
        com += mass * coord
        total_mass += mass
    
    com /= total_mass
    
    # Apply shift to all atoms
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        coord = np.array([pos.x, pos.y, pos.z])
        new_coord = coord - com 
        if target_com is not None: 
            new_coord = new_coord + target_com
        conf.SetAtomPosition(atom.GetIdx(), new_coord)
    
    # Save shifted structure
    writer = Chem.SDWriter(output_sdf)
    writer.write(mol)
    writer.close()
    
    return output_sdf


        
def process_pdb_list(csv_file: str, 
                    save_loc: str, 
                    key_residues_dict: list = None, 
                    chain: str = None
                    ) -> None:
    """
    Process PDB structures and save active site data as NPZ files.

    Args:
        csv_file (str): Path to CSV file containing PDB data
        save_loc (str): Base directory to save processed PDB files
        key_residues_dict (list): List of dictionaries with residue info for filtering
        chain (str, optional): Specific chain to filter for. If None, keeps all chains.
    """
    # Set up logging
    logging.basicConfig(level=logging.INFO,
                       format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)

    # Create save directory
    save_loc = Path(save_loc)
    save_loc.mkdir(parents=True, exist_ok=True)
    
    # Read CSV file
    try:
        df = pd.read_csv(csv_file)
        logger.info(f"Successfully loaded {len(df)} entries from CSV")
    except Exception as e:
        logger.error(f"Error reading CSV file: {str(e)}")
        raise
    
    # Process each PDB
    for idx, row in df.iterrows():
        pdb = row['rcsb_id']
        logger.info(f"Processing {pdb}...")
        
        try:
            # Create PDB-specific directory
            pdb_dir = save_loc / f"{pdb}"
            if chain:
                pdb_dir = save_loc / f"{pdb}_chain_{chain}"
            pdb_dir.mkdir(parents=True, exist_ok=True)
            
            # Download and process PDB
            pdb_f = download_pdb_to_local(pdb, save_dir=str(pdb_dir))

            if chain:
                chain_pdb = pdb_dir / f"{pdb}_chain_{chain}.pdb"
                pdb_f = extract_chain(pdb_f, chain, str(chain_pdb))

            pdb_filtered_f = filter_to_pdb_hetatm_types(
                pdb_f,
                hetatms=key_residues_dict if key_residues_dict else [],
                chain=chain,
                f_out=str(pdb_dir / f"{pdb}_filtered.pdb")
            )
            
            final_pdb = str(pdb_dir / f"{pdb}_altloc_removed.pdb")
            remove_altloc(pdb_filtered_f, output_pdb=final_pdb)
            
            # Parse active site residues
            act_site = ast.literal_eval(row['pdb_site_resis'])
            logger.info(f"Active site residues for {pdb}: {act_site}")

            # Save active site residues as simple text file
            with open(pdb_dir / "active_site.txt", "w") as f:
                f.write(" ".join(act_site))
            
            logger.info(f"Successfully processed {pdb}")
            
        except Exception as e:
            logger.error(f"Error processing {pdb}: {str(e)}")
            continue



# # Usage examples:
# # Process all chains
# process_pdb_list(
#     "./data/initial_pdbs.txt",
#     "./data/pdbs/initial_clean_pdbs/",
#     key_residues_dict=[]
# )

# # Process only chain A
# process_pdb_list(
#     "./data/initial_pdbs.txt",
#     "./data/pdbs/initial_clean_pdbs/",
#     key_residues_dict=[],
#     chain='A'
# )

