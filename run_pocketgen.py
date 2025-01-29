from xyme_tools.structure_generation import (
    StructureGenerationInput,
    StructureGenerationOutput,
    PocketGen,
)
from xyme_tools.general.custom_typing import subdir_to_path, path_to_subdir
from xyme_tools.general.utils import set_working_dir

set_working_dir("./pipeline_output")

pg = PocketGen()

input_data_list = [
    StructureGenerationInput(
        pdb_file = './data/processed_pdbs/6osz/6osz_final.pdb',
        output_subdir = path_to_subdir('./pipeline_output/test/pocketgen')
    ),
]

job_name = pg.run(input_data_list=input_data_list, wait=False, backend='tamarind')

print(job_name)

pg_output = pg.retrieve(input_data_list=input_data_list, backend='tamarind', job_name=job_name, force_download=False)