from src.utils.run_pocketgen import submit_pocketgen_jobs

# Usage example:
submit_pocketgen_jobs(
    pdb_dir="./data/pdbs/initial_clean_pdbs",
    sdf_dir="./data/ligand_sdf",
    base_output_dir="./output_pipeline"
)