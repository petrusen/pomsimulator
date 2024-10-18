# Standard library imports
import time
import numpy as np
import os
# Local imports
from pomsimulator.modules.text_module import Read_csv,Print_logo
from pomsimulator.modules.graph_module import *
from pomsimulator.modules.msce_module import *
from pomsimulator.modules.DataBase import *
from pomsimulator.modules.helper_module import *

os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'


def main():
    Print_logo()
    ######################### User parameters ##############################################
    # Labels and species
    formation_labels = Labels_PMo
    speciation_labels = Labels_PMo
    ref_compounds = ["P00Mo01O04-0H", "P01Mo00O04-0H"]

    # Chemical parameters
    min_pH, max_pH, step_pH = 0, 14, -0.1
    pH = np.arange(max_pH, min_pH, step_pH)
    C_X = 0.01
    Ratio_list = np.linspace(1, 12, 20)

    # Operation parameters
    all_idxs = False
    batch_size = 1
    cores = 8

    # Input/output files
    path = "../outputs/logkf_PMo.txt"
    output_fold = "PMo_phase_arrays"
    output_path = "../outputs/" + output_fold
    npz_info_file = "PMo_npz_info.dat"
    scaling_path = "../outputs/scaling_params_PMo.pomsim"

    if not os.path.exists(output_path):
        os.makedirs(output_path)
    else:
        counter = 1
        while os.path.exists(output_path):
            output_path = "../outputs/" + output_fold + ".%03d"%counter
            counter += 1
        else:
            os.makedirs(output_path)

    #############################################################################################
    # Read linear scaling from test_linearity
    scaling_params = read_scaling_params(scaling_path)

    # Read constants and scale them
    ref_stoich_X,ref_stoich_M = [Lab_to_stoich(lab)  for lab in ref_compounds]
    lgkf_df = Read_csv(path, formation_labels)
    lgkf_df = lgkf_filtering(lgkf_df,all_idxs,scaling_params, speciation_labels)

    mapping_string = ""
    for ii,Ratio in enumerate(Ratio_list):
        print("Speciation for ratio = %.6f" % Ratio)
        C_M = C_X * Ratio
        speciation_array, IndexArray = compute_speciation_loop(lgkf_df, speciation_labels, pH, [C_X,C_M],
                                                               [ref_stoich_X,ref_stoich_M],
                                                             None, batch_size, cores, show_progress=False)
        file_name = output_path + "/array_%02d.npz" % ii
        np.savez_compressed(file_name,SupArray=speciation_array,IndexArray=IndexArray,
                            pH=pH,C_X=C_X,C_M=C_M,labels=speciation_labels)
        mapping_string += file_name + "\n"
    with open(npz_info_file,"w") as outfile:
        outfile.write(mapping_string)
if __name__ == '__main__':
    main()
