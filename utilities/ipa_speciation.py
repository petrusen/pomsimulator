# Standard library imports
import numpy as np
import os
from multiprocessing import Pool, cpu_count
import time
from itertools import repeat

# Local imports
from pomsimulator.modules.text_module import Print_logo,Read_csv,Lab_to_stoich,write_speciationparameters
from pomsimulator.modules.msce_module import Speciation_from_Formation_singlemetal,starmap_with_kwargs
from pomsimulator.modules.DataBase import *
from pomsimulator.modules.helper_module import *

os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'



def main():

    Print_logo()
    ######################### User parameters ##############################################
    # Labels and species
    formation_labels = Labels_W
    speciation_labels = Labels_W_good
    ref_compound = "W01O04-0H"

    # Chemical parameters
    min_pH, max_pH, step_pH = 0, 14, -0.1
    pH = np.arange(max_pH,min_pH,step_pH)
    C = 0.1

    # Operation parameters
    all_idxs = True
    batch_size = 50
    cores = 8

    # Input/output files
    output_path = "../outputs/W_data"

    path = output_path +  "/logkf_W.txt"
    path_to_output = output_path + "/W_Array.npz"
    path_to_params = output_path + "/speciation_parameters_W.txt"
    scaling_path = output_path + "/scaling_params_W.pomsim"
    #############################################################################################

    # Read linear scaling from test_linearity
    scaling_params = read_scaling_params(scaling_path)

    # Read constants and scale them
    ref_stoich = Lab_to_stoich(ref_compound)
    lgkf_df = Read_csv(path)
    lgkf_df = apply_lgkf_scaling(lgkf_df,scaling_params, speciation_labels)

    ### Write parameters to file once we know the number of models
    kwargs_input = dict()
    obj_list = [path_to_params,path,scaling_params["m"],scaling_params["b"],scaling_params["mode"],cores,
                C,(min_pH,max_pH),abs(step_pH),len(list(lgkf_df.index)),
                speciation_labels,ref_stoich,path_to_output]
    for s,o in zip(speciation_parameters_strings, obj_list):
        kwargs_input[s] = o
    write_speciationparameters(kwargs_input)

    FilteredSuperArr,IndexArr = compute_speciation_loop(lgkf_df, speciation_labels, pH, C, ref_stoich, path_to_output, batch_size, cores)


if __name__ == '__main__':
    main()
