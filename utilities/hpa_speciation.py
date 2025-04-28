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
from configparser import ConfigParser

os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'



def main():

    Print_logo()
    config_file = "../inputs/config_PMo.pomsim"
    config = ConfigParser()
    config.read(config_file)
    #### Parameters
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
    C_X = float(config["Speciation"]["C_X"])
    C_M = float(config["Speciation"]["C_M"])
    # Operation parameters
    cores = int(config["Speciation"]["cores"])
    batch_size = int(config["Speciation"]["batch_size"])

    # Input/output files
    path = output_path + "/logkf_%s.csv" % system
    path_to_output = output_path + "/" + "Array_%s.npz" % system
    path_to_params = output_path + "/" + "speciation_params_%s.txt" % system
    scaling_path = output_path + "/scaling_params_%s.pomsim" % system
    #############################################################################################

    # Read linear scaling from test_linearity
    scaling_params = read_scaling_params(scaling_path)

    # Read constants and scale them
    ref_stoich_X, ref_stoich_M = [Lab_to_stoich(ref) for ref in ref_compounds]

    lgkf_df = Read_csv(path)
    lgkf_df = apply_lgkf_scaling(lgkf_df,scaling_params, speciation_labels)

    ### Write parameters to file once we know the number of models
    kwargs_input = dict()
    obj_list = [path_to_params,path,scaling_params["m"],scaling_params["b"],scaling_params["mode"],cores,
                [C_X,C_M],(min_pH,max_pH),abs(step_pH),len(list(lgkf_df.index)),
                speciation_labels,[ref_stoich_X,ref_stoich_M],path_to_output]
    for s,o in zip(speciation_parameters_strings, obj_list):
        kwargs_input[s] = o
    write_speciationparameters(kwargs_input)

    FilteredSuperArr,IndexArr = compute_speciation_loop(lgkf_df, speciation_labels, pH, [C_X,C_M], [ref_stoich_X,ref_stoich_M],
                                                        path_to_output, batch_size, cores)


if __name__ == '__main__':
    main()
