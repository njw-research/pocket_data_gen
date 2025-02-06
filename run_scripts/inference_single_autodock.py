from xyme_tools.structure_prediction import (
    StructurePredictionInput,
    StructurePredictionOutput,
)
from xyme_tools.structure_prediction.tool_classes import AutoDock
from xyme_tools.general.custom_typing import subdir_to_path, path_to_subdir
from xyme_tools.general.utils import set_working_dir
import numpy as np
from src.applications.autodock.run_docking import setup_docking

set_working_dir("./output_pipeline/")

pdb_file, sdf_file, act_site = setup_docking('examples/6OSZ_chain_A')

ad = AutoDock(exhaustiveness=32)

input_data_list = [
    StructurePredictionInput(
        sequence = 'A', # Does not matter but needs to be here 
        pdb_file = pdb_file,
        sdf_file = sdf_file,
        output_subdir = path_to_subdir('./autodock_test/'),
        box_xyz = act_site,
        box_dims = np.array([20.0, 20.0, 20.0]),
    ),
]

job_name = ad.run(input_data_list=input_data_list, wait=False, backend='tamarind')

print(job_name)

# ad_output = ad.retrieve(input_data_list=input_data_list, backend='tamarind', job_name=job_name, force_download=False)