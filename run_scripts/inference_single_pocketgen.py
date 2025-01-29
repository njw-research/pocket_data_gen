from xyme_tools.structure_generation import (
    StructureGenerationInput,
    PocketGen,
)
from xyme_tools.general.custom_typing import path_to_subdir
from xyme_tools.general.utils import set_working_dir

set_working_dir("./output_pipeline")

pg = PocketGen()

input_data_list = [
    StructureGenerationInput(
        pdb_file = './data/pdbs/initial_clean_pdbs/6OSZ/6OSZ_altloc_removed.pdb',
        output_subdir = path_to_subdir('./output_pipeline/pocketgen'),
        sdf_file = './data/ligand_sdf/LTV.sdf'
    ),
]

job_name = pg.run(input_data_list=input_data_list, wait=False, backend='tamarind')

print(job_name)

pg_output = pg.retrieve(input_data_list=input_data_list, backend='tamarind', job_name=job_name, force_download=False)