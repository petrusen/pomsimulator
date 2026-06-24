# Standard library imports
import os
from multiprocessing import Pool, cpu_count
import time
import datetime
from configparser import ConfigParser
import pkg_resources as pkgr
#Third-party imports
from itertools import repeat
import numpy as np
# Local imports
from pomsimulator.modules.text_module import Print_logo,Read_csv,Lab_to_stoich,write_speciation_parameters
from pomsimulator.modules.msce_module import Speciation_from_Formation_singlemetal,starmap_with_kwargs
from pomsimulator.modules.DataBase import *
from pomsimulator.modules.helper_module import *

def speciation_ipa(config_dict):
    """
    Speciation of ipa.

    Args:
        config_dict (dict or str): Configuration dictionary containing parameters for the function.
    Returns:
        str: Message indicating the completion of the function.
    """
    try:
        Print_logo()
        run_by_gui = True

        if type(config_dict) != dict:
            config_file_path = pkgr.resource_filename(__name__, config_dict)
            config_dict = get_config_to_dict(config_file_path)
            run_by_gui = False

        ######################### User parameters ##############################################

        output_path = pkgr.resource_filename(__name__,config_dict["Preparation"]["output_path"])
        system = config_dict["Preparation"]["POM_system"]

        if run_by_gui:
            stop_signal = config_dict["_stop_file"]
        else:
            stop_signal = ""

        # Labels and species

        speciation_labels = config_dict["Speciation"]["speciation_labels"].split(",")
        labels_file = output_path + "/labels_%s.txt" % system
        if speciation_labels[0] == "all":
            with open(labels_file,"r") as flab:
                speciation_labels = [item.strip() for item in flab.readlines()]
        ref_compound = config_dict["Simulation"]["ref_compound"]

        # Chemical parameters

        min_pH, max_pH, step_pH = [float(config_dict["Speciation"][prop]) for prop in ["min_pH","max_pH","step_pH"]]
        pH = np.arange(max_pH,min_pH,-step_pH)
        C = float(config_dict["Speciation"]["C"])

        # Operation parameters

        cores = int(config_dict["Speciation"]["cores"])
        batch_size = int(config_dict["Speciation"]["batch_size"])

        # Input/output files

        path = output_path + "/logkf_%s.csv" % system
        path_to_output = output_path + "/" + "Array_%s.npz" % system
        path_to_params = output_path + "/" + "speciation_params_%s.txt" % system
        scaling_path = output_path + "/scaling_params_%s.pomsim" % system


        # 1) Read linear scaling ############################################################################################
        print("1) Read linear scaling")

        start_date = str(datetime.datetime.now())
        start_time = time.time()

        scaling_params = read_scaling_params(scaling_path)

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        # 2) Read constants and scale them ############################################################################################
        print("2) Read constants and scale them")

        ref_stoich = Lab_to_stoich(ref_compound)
        lgkf_df = Read_csv(path)
        lgkf_df = apply_lgkf_scaling(lgkf_df,scaling_params, speciation_labels)

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        # 3) Compute speciation loop ############################################################################################
        print("3) Compute speciation loop")

        result = compute_speciation_loop(lgkf_df, speciation_labels, pH, C, ref_stoich, path_to_output, batch_size, cores, stop_signal=stop_signal)

        if result == "Stopped":
            print("Function stopped by user")
            return "Stopped"


        # 4) Write parameters to file
        print("4) Write parameters to file")

        end_date = str(datetime.datetime.now())
        end_time = time.time()
        timing = end_time - start_time

        kwargs_input = dict()

        obj_list = [path_to_params,path,scaling_params["m"],scaling_params["b"],scaling_params["mode"],cores,
                    C,(min_pH,max_pH),abs(step_pH),len(list(lgkf_df.index)),
                    speciation_labels,ref_stoich,path_to_output,start_date, end_date, timing]

        for s, o in zip(speciation_parameters_strings, obj_list):
            kwargs_input[s] = o

        write_speciation_parameters(kwargs_input)

        return "Normal termination"
    except Exception as e:
        print("An error occurred during the speciation diagram generation: ", str(e))
        print(f"Error {e}")
        return "Error"

def speciation_hpa(config_dict):
    """
    Speciation of hpa using the MSCE method.

    Args:
        config_dict (dict or str): Configuration dictionary containing parameters for the function.
    Returns:
        str: Message indicating the completion of the function.
    """
    try:
        Print_logo()
        run_by_gui = True

        if type(config_dict) != dict:
            config_file_path = pkgr.resource_filename(__name__, config_dict)
            config_dict = get_config_to_dict(config_file_path)
            run_by_gui = False

        ######################### User parameters ##############################################

        output_path = pkgr.resource_filename(__name__,config_dict["Preparation"]["output_path"])
        system = config_dict["Preparation"]["POM_system"]

        if run_by_gui:
            stop_signal = config_dict["_stop_file"]
        else:
            stop_signal = ""

        # Labels and species

        speciation_labels = config_dict["Speciation"]["speciation_labels"].split(",")
        labels_file = output_path + "/labels_%s.txt" % system
        if speciation_labels[0] == "all":
            with open(labels_file, "r") as flab:
                speciation_labels = [item.strip() for item in flab.readlines()]

        ref_compounds = config_dict["Simulation"]["ref_compound"].split(",")

        # Chemical parameters

        min_pH, max_pH, step_pH = [float(config_dict["Speciation"][prop]) for prop in ["min_pH", "max_pH", "step_pH"]]
        pH = np.arange(max_pH, min_pH, -step_pH)
        C_X = float(config_dict["Speciation"]["C_X"])
        C_M = float(config_dict["Speciation"]["C_M"])

        # Operation parameters

        cores = int(config_dict["Speciation"]["cores"])
        batch_size = int(config_dict["Speciation"]["batch_size"])

        # Input/output files

        path = output_path + "/logkf_%s.csv" % system
        path_to_output = output_path + "/" + "Array_%s.npz" % system
        path_to_params = output_path + "/" + "speciation_params_%s.txt" % system
        scaling_path = output_path + "/scaling_params_%s.pomsim" % system

        # 1) Read linear scaling from test_linearity ############################################################################################
        print("1) Read linear scaling from test_linearity")

        start_date = str(datetime.datetime.now())
        start_time = time.time()

        scaling_params = read_scaling_params(scaling_path)

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        # 2) Read constants and scale them ############################################################################################
        print("2) Read constants and scale them")

        ref_stoich_X, ref_stoich_M = [Lab_to_stoich(ref) for ref in ref_compounds]
        lgkf_df = Read_csv(path)
        lgkf_df = apply_lgkf_scaling(lgkf_df, scaling_params, speciation_labels)

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"


        # 3) Compute speciation loop ############################################################################################
        print("3) Compute speciation loop")

        result = compute_speciation_loop(lgkf_df, speciation_labels, pH, [C_X, C_M],
                                                             [ref_stoich_X, ref_stoich_M],
                                                             path_to_output, batch_size, cores)
        if result == "Stopped":
            print("Function stopped by user")
            return "Stopped"

        # 4) Write parameters to file
        print("4) Write parameters to file")

        end_date = str(datetime.datetime.now())
        end_time = time.time()
        timing = end_time - start_time

        kwargs_input = dict()
        obj_list = [path_to_params, path, scaling_params["m"], scaling_params["b"], scaling_params["mode"], cores,
                    [C_X, C_M], (min_pH, max_pH), abs(step_pH), len(list(lgkf_df.index)),
                    speciation_labels, [ref_stoich_X, ref_stoich_M], path_to_output,
                    start_date, end_date, timing]
        for s, o in zip(speciation_parameters_strings, obj_list):
            kwargs_input[s] = o
        write_speciation_parameters(kwargs_input)

        return "Normal termination"
    except Exception as e:
        print("An error occurred during the simulation: ", str(e))
        print(f"Error {e}")
        return "Error"

def phase_ipa(config_dict):
    """
        Speciation phase for ipa.

        Args:
            config_dict (dict or str): Configuration dictionary containing parameters for the function.
        Returns:
            str: Message indicating the completion of the function.
        """

    try:
        Print_logo()
        run_by_gui = True

        if type(config_dict) != dict:
            config_file_path = pkgr.resource_filename(__name__, config_dict)
            config_dict = get_config_to_dict(config_file_path)
            run_by_gui = False

        ######################### User parameters ##############################################

        output_path = pkgr.resource_filename(__name__, config_dict["Preparation"]["output_path"])
        system = config_dict["Preparation"]["POM_system"]

        if run_by_gui:
            stop_signal = config_dict["_stop_file"]
        else:
            stop_signal = ""

        # Labels and species

        speciation_labels = config_dict["Speciation"]["speciation_labels"].split(",")
        labels_file = output_path + "/labels_%s.txt" % system
        if speciation_labels[0] == "all":
            with open(labels_file,"r") as flab:
                speciation_labels = [item.strip() for item in flab.readlines()]
        ref_compound = config_dict["Simulation"]["ref_compound"]

        # Chemical parameters
        min_pH, max_pH, step_pH = [float(config_dict["Speciation"][prop]) for prop in ["min_pH","max_pH","step_pH"]]
        pH = np.arange(max_pH,min_pH,-step_pH)
        min_logC,max_logC = [float(config_dict["Speciation"][prop]) for prop in ["min_logC","max_logC"]]
        N_logC = int(float(config_dict["Speciation"]["num_logC"]))

        C_list = np.logspace(min_logC,max_logC,N_logC)

        # Operation parameters
        cores = int(config_dict["Speciation"]["cores"])
        batch_size = int(config_dict["Speciation"]["batch_size"])

        # Input/output files
        output_fold = config_dict["Speciation"]["phase_dir"]
        lgkf_path = output_path + "/logkf_%s.csv" % system
        scaling_path = output_path + "/scaling_params_%s.pomsim" % system

        output_phase_path = output_path + "/" + output_fold
        npz_info_file = output_path + "/npz_info_%s.dat" % system
        model_subset_file = config_dict["Speciation"]["model_subset_file"]

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

        # 1) Read linear scaling ############################################################################################
        print("1) Read linear scaling")

        scaling_params = read_scaling_params(scaling_path)

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        # 2) Read constants and scale them ############################################################################################
        print("2) Read constants and scale them")

        # Read constants and scale them
        ref_stoich = Lab_to_stoich(ref_compound)
        lgkf_df = Read_csv(lgkf_path)
        lgkf_df = apply_lgkf_scaling(lgkf_df,scaling_params, speciation_labels)

        if model_subset_file:
            with open(output_path + "/" + model_subset_file,"r") as fmod:
                model_sel = [int(item.strip()) for item in fmod.readlines()]
            lgkf_df = lgkf_df.loc[model_sel,:]

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        # 3) Compute speciation loop ############################################################################################
        print("3) Compute speciation loop")

        mapping_string = ""
        for ii,C in enumerate(C_list):

            if stop_signal and os.path.exists(stop_signal):
                print("Function stopped by user")
                return "Stopped"

            print("Speciation for concentration = %.6f" % C)
            result = compute_speciation_loop(lgkf_df, speciation_labels, pH, C, ref_stoich,
                                                                 None, batch_size, cores, show_progress=False, stop_signal=stop_signal)
            if result == "Stopped":
                print("Stopped by user at concentration = %.6f" % C)
                return "Stopped"

            speciation_array, IndexArray = result

            file_name = output_phase_path + "/array_%02d.npz" % ii
            np.savez_compressed(file_name,SupArray=speciation_array,IndexArray=IndexArray,
                                pH=pH,C=C,labels=speciation_labels)
            mapping_string += file_name + "\n"

        with open(npz_info_file,"w") as outfile:
            outfile.write(mapping_string)

        return "Normal termination"

    except Exception as e:
        print("An error occurred during the phase diagram generation: ", str(e))
        print(f"Error {e}")
        return "Error"

def phase_hpa(config_dict):
    try:
        Print_logo()
        run_by_gui = True

        if type(config_dict) != dict:
            config_file_path = pkgr.resource_filename(__name__, config_dict)
            config_dict = get_config_to_dict(config_file_path)
            run_by_gui = False

        ######################### User parameters ##############################################

        output_path = pkgr.resource_filename(__name__, config_dict["Preparation"]["output_path"])
        system = config_dict["Preparation"]["POM_system"]

        if run_by_gui:
            stop_signal = config_dict["_stop_file"]
        else:
            stop_signal = ""

        # Labels and species

        speciation_labels = config_dict["Speciation"]["speciation_labels"].split(",")
        labels_file = output_path + "/labels_%s.txt" % system
        if speciation_labels[0] == "all":
            with open(labels_file,"r") as flab:
                speciation_labels = [item.strip() for item in flab.readlines()]

        ref_compounds = config_dict["Simulation"]["ref_compound"].split(",")

        # Chemical parameters
        min_pH, max_pH, step_pH = [float(config_dict["Speciation"][prop]) for prop in ["min_pH","max_pH","step_pH"]]
        pH = np.arange(max_pH,min_pH,-step_pH)
        min_Ratio,max_Ratio = [float(config_dict["Speciation"][prop]) for prop in ["min_Ratio","max_Ratio"]]
        N_Ratio = int(float(config_dict["Speciation"]["num_Ratio"]))
        C_X = float(config_dict["Speciation"]["C_X"])
        Ratio_list = np.linspace(min_Ratio, max_Ratio, N_Ratio)

        # Operation parameters
        cores = int(config_dict["Speciation"]["cores"])
        batch_size = int(config_dict["Speciation"]["batch_size"])

        # Input/output files
        output_fold = config_dict["Speciation"]["phase_dir"]
        lgkf_path = output_path + "/logkf_%s.csv" % system
        scaling_path = output_path + "/scaling_params_%s.pomsim" % system

        output_phase_path = output_path + "/" + output_fold
        npz_info_file = output_path + "/npz_info_%s.dat" % system
        model_subset_file = config_dict["Speciation"]["model_subset_file"]

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

        # 1) Read linear scaling ############################################################################################
        print("1) Read linear scaling")

        scaling_params = read_scaling_params(scaling_path)

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        # 2) Read constants and scale them ############################################################################################
        print("2) Read constants and scale them")

        ref_stoich_X,ref_stoich_M = [Lab_to_stoich(lab) for lab in ref_compounds]
        lgkf_df = Read_csv(lgkf_path)
        lgkf_df = apply_lgkf_scaling(lgkf_df,scaling_params, speciation_labels)

        if model_subset_file:
            with open(output_path + "/" + model_subset_file,"r") as fmod:
                model_sel = [int(item.strip()) for item in fmod.readlines()]
            lgkf_df = lgkf_df.loc[model_sel,:]

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        # 3) Compute speciation loop ############################################################################################
        print("3) Compute speciation loop")

        mapping_string = ""
        for ii,Ratio in enumerate(Ratio_list):

            if stop_signal and os.path.exists(stop_signal):
                print("Function stopped by user")
                return "Stopped"

            print("Speciation for ratio = %.6f" % Ratio)
            C_M = C_X * Ratio
            result = compute_speciation_loop(lgkf_df, speciation_labels, pH, [C_X,C_M],
                                                                   [ref_stoich_X,ref_stoich_M],
                                                                 None, batch_size, cores, show_progress=False, stop_signal=stop_signal)

            if result == "Stopped":
                print("Stopped by user at concentration = %.6f" % C_M)
                return "Stopped"

            speciation_array, IndexArray = result

            file_name = output_phase_path + "/array_%02d.npz" % ii
            np.savez_compressed(file_name,SupArray=speciation_array,IndexArray=IndexArray,
                                pH=pH,C_X=C_X,C_M=C_M,labels=speciation_labels)
            mapping_string += file_name + "\n"

        with open(npz_info_file,"w") as outfile:
            outfile.write(mapping_string)

        return("Normal termination")

    except Exception as e:
        print("An error occurred during the phase diagram generation: ", str(e))
        print(f"Error {e}")
        return("Error")