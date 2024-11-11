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
from configparser import ConfigParser

os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'


def main():
    Print_logo()
    config_file = "../inputs/config_PMo.pomsim"
    config = ConfigParser()
    config.read(config_file)
    ######################### User parameters ##############################################
    output_path = config["Preparation"]["output_path"]
    system = config["Preparation"]["POM_system"]
    # Labels and species
    speciation_labels = config["Speciation"]["speciation_labels"].split(",")
    labels_file = output_path + "/labels_%s.txt" % system
    if speciation_labels[0] == "all":
        with open(labels_file,"r") as flab:
            speciation_labels = [item.strip() for item in flab.readlines()]

    ref_compounds = config["Simulation"]["ref_compound"].split(",")


    # Chemical parameters
    min_pH, max_pH, step_pH = [float(config["Speciation"][prop]) for prop in ["min_pH","max_pH","step_pH"]]
    pH = np.arange(max_pH,min_pH,-step_pH)
    min_Ratio,max_Ratio = [float(config["Speciation"][prop]) for prop in ["min_Ratio","max_Ratio"]]
    N_Ratio = int(config["Speciation"]["num_Ratio"])
    C_X = float(config["Speciation"]["C_X"])
    Ratio_list = np.linspace(min_Ratio, max_Ratio, N_Ratio)

    # Operation parameters
    cores = int(config["Speciation"]["cores"])
    batch_size = int(config["Speciation"]["batch_size"])

    # Input/output files
    output_fold = config["Speciation"]["phase_dir"]
    lgkf_path = output_path + "/logkf_%s.csv" % system
    scaling_path = output_path + "/scaling_params_%s.pomsim" % system

    output_phase_path = output_path + "/" + output_fold
    npz_info_file = output_path + "/npz_info_%s.dat" % system
    model_subset_file = config["Speciation"]["model_subset_file"]

    if not os.path.exists(output_phase_path):
        os.makedirs(output_phase_path)
    else:
        counter = 1
        base_path = output_phase_path
        while os.path.exists(output_phase_path):
            output_phase_path = base_path + ".%03d" % counter
            counter += 1
        else:
            os.makedirs(output_phase_path)

    #############################################################################################
    # Read linear scaling from test_linearity
    scaling_params = read_scaling_params(scaling_path)

    # Read constants and scale them
    ref_stoich_X,ref_stoich_M = [Lab_to_stoich(lab) for lab in ref_compounds]
    lgkf_df = Read_csv(lgkf_path)
    lgkf_df = apply_lgkf_scaling(lgkf_df,scaling_params, speciation_labels)
    if model_subset_file:
        with open(output_path + "/" + model_subset_file,"r") as fmod:
            model_sel = [int(item.strip()) for item in fmod.readlines()]

        lgkf_df = lgkf_df.loc[model_sel,:]
        print(lgkf_df.shape)

    mapping_string = ""
    for ii,Ratio in enumerate(Ratio_list):
        print("Speciation for ratio = %.6f" % Ratio)
        C_M = C_X * Ratio
        speciation_array, IndexArray = compute_speciation_loop(lgkf_df, speciation_labels, pH, [C_X,C_M],
                                                               [ref_stoich_X,ref_stoich_M],
                                                             None, batch_size, cores, show_progress=False)
        file_name = output_phase_path + "/array_%02d.npz" % ii
        np.savez_compressed(file_name,SupArray=speciation_array,IndexArray=IndexArray,
                            pH=pH,C_X=C_X,C_M=C_M,labels=speciation_labels)
        mapping_string += file_name + "\n"
    with open(npz_info_file,"w") as outfile:
        outfile.write(mapping_string)
if __name__ == '__main__':
    main()
