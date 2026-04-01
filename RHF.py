import numpy as np

import build_data
import integrals



def RHF(input_file, base_function_file):

    molecule_data=build_data.load_xyz(input_file, to_bohr=True)
    molecule_basis=build_data.build_molecule_basis(molecule_data, base_function_file)

    