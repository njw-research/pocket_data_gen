import os

from Bio.PDB import PDBParser, PDBIO
import logging
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


def process_pdb_list(pdb_list_file: str, save_loc: str, key_residues_dict: list = None, chain: str = None) -> None:
   """
   Process a list of PDB IDs by downloading, filtering and cleaning the structures.
   
   Args:
       pdb_list_file (str): Path to file containing PDB IDs (one per line)
       save_loc (str): Base directory to save processed PDB files
       key_residues_dict (list): List of dictionaries with residue info for filtering
       chain (str, optional): Specific chain to filter for. If None, keeps all chains.
   """
   # Create save directory if it doesn't exist
   if not os.path.exists(save_loc):
       os.makedirs(save_loc)
       
   # Read PDB list
   with open(pdb_list_file, 'r') as f:
       pdb_list = [line.strip() for line in f.readlines()]
   
   # Process each PDB
   for pdb in pdb_list:
       print(f"Processing {pdb}...")
       
       # Create PDB-specific directory
       pdb_dir = os.path.join(save_loc, f"{pdb}")
       if chain:
           pdb_dir = os.path.join(save_loc, f"{pdb}_chain_{chain}")
       if not os.path.exists(pdb_dir):
           os.makedirs(pdb_dir)
           
       try:
           # Download PDB
           pdb_f = download_pdb_to_local(pdb, save_dir=pdb_dir)

           # Extract specific chain if requested
           if chain:
               chain_pdb = os.path.join(pdb_dir, f"{pdb}_chain_{chain}.pdb")
               pdb_f = extract_chain(pdb_f, chain, chain_pdb)
           
           # Filter HETATM types
           pdb_filtered_f = filter_to_pdb_hetatm_types(
               pdb_f,
               hetatms=key_residues_dict if key_residues_dict else [],
               chain=chain,
               f_out=os.path.join(pdb_dir, f"{pdb}_filtered.pdb")
           )

           
           # Remove alternate locations
           remove_altloc(
               pdb_filtered_f,
               output_pdb=os.path.join(pdb_dir, f"{pdb}_altloc_removed.pdb")
           )
           
           print(f"Successfully processed {pdb}")
           
       except Exception as e:
           print(f"Error processing {pdb}: {str(e)}")
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