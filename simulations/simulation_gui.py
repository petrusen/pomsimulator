# Standard library imports
import sys
import os
from os import listdir
from os.path import isfile, join
import time
import datetime
import random
import pkg_resources as pkgr
import numpy as np
# Local imports
from pomsimulator.modules.text_module import Print_logo,write_simulation_parameters,Write_Reactions,read_diagonal
from pomsimulator.modules.DataBase import *
from pomsimulator.modules.helper_module import *

from configparser import ConfigParser

os.environ['MKL_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'


def simulation_ipa(config_dict):

    """
    Run an isopolyanion (IPA) simulation based on the provided configuration dictionary.

    Args:
        config_dict (dict): Configuration dictionary
    Returns:
        str: Status message
    """

    try:
        run_by_gui = True

        if type(config_dict) != dict:
            config_file_path = pkgr.resource_filename(__name__, config_dict)
            config_dict = get_config_to_dict(config_file_path)
            run_by_gui = False

        ######################### User parameters ##############################################

        # Input/output files
        system = config_dict["Preparation"]["POM_system"]
        ADF_folder = pkgr.resource_filename(__name__,config_dict["Preparation"]["adf_inputs_dir"])
        mol_folder = pkgr.resource_filename(__name__,config_dict["Preparation"]["mol_folder"])
        output_path = pkgr.resource_filename(__name__,config_dict["Preparation"]["output_path"])
        # Create output directory if it doesn't exist'
        os.makedirs(output_path, exist_ok=True)

        if run_by_gui:
            stop_signal = config_dict["_stop_file"]
        else:
            stop_signal = ""

        isomorphism_matrix = output_path + "/np_IM_%s.csv" % system
        formation_constants_file = output_path + "/logkf_%s.csv" % system
        CRN_file = output_path + "/CRN_%s.txt" % system
        labels_file = output_path + "/labels_%s.txt" % system
        simulation_file = output_path + "/simulation_parameters_%s.txt" % system

        # Operation parameters -> all from config_dict file
        cores = int(config_dict['Simulation']["cores"])
        batch_size = int(config_dict['Simulation']["batch_size"])
        sample_perc = float(config_dict['Simulation']["sample_perc"])
        sample_type = config_dict['Simulation']["sample_type"]  # random / all

        # Chemical parameters
        use_isomorphisms = config_dict['Simulation']["use_isomorphism"]
        energy_threshold = float(config_dict['Simulation']["energy_threshold"])  # Maximum value for reaction energies
        proton_numb = int(config_dict['Simulation']["proton_numb"])  # Maximum difference in proton number of to species to react
        reference = config_dict["Simulation"]["reference_types"].split(',')

        I, C0, temp = [float(config_dict['Simulation'][prop]) for prop in ["I", "C0", "temp"]]
        min_pH, max_pH, step_pH = [float(config_dict['Simulation'][prop]) for prop in ["min_pH", "max_pH", "step_pH"]]
        ref_compound = config_dict['Simulation']["ref_compound"]


        # Internal parameters: These parameters are not meant to be routinely modified
        internal_conditions = config_dict["InternalConditions"]
        internal_conditions["proton_numb"] = proton_numb

        Print_logo()

        # 1) Get ADF outputs ############################################################################################
        start_date = str(datetime.datetime.now())
        start_time = time.time()
        print("1) Get ADF outputs")

        # # Read ADF outputs and generate graphs
        adf_files = sorted([ADF_folder + "/" + f for f in listdir(ADF_folder) if isfile(join(ADF_folder, f))])

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        # 2) Graph creation #############################################################################################

        print("2) Graph creation: node=atom edge=bond")
        G1_list, G1_labels, graphs_info = generate_graphs(adf_files, ref_compound, system)

        with open(labels_file, "w") as flab:
            flab.write("\n".join(G1_labels))

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        # 3) Isomorphism matrix generation ##############################################################################
        print("3) Isomorphism matrix and reactions generation")

        if use_isomorphisms:  # Calculate all isomorphisms
            diagonal = read_diagonal(isomorphism_matrix)
        else:  # Assume that all species are isomorphic
            diagonal = np.tri(graphs_info["num_molec"], graphs_info["num_molec"], 0)

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        reac_idx, reac_energy, reac_type = Isomorphism_to_ChemicalReactions_gui(G1_list, diagonal, graphs_info["water"],
                                                                            reference, system,
                                                                            energy_threshold, internal_conditions)
        print("3.1) NUMBER OF REACTIONS:", [len(list(map(len, reac_idx[i]))) for i in range(len(reac_idx))])

        R_idx, R_ene, R_type = sort_by_type(reac_idx, reac_energy, reac_type, graphs_info["compounds_set"],
                                            graphs_info["unique_labels"],
                                            G1_labels, ref_compound)

        Write_Reactions(CRN_file, G1_labels, reac_idx, reac_type, reac_energy, stringreac_dict, molecularity_dict)

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        # 4) Speciation #################################################################################################
        print("4) Solving MSCE Models and generate parameters' output")

        e_ctt = list(reac_energy[0])
        idx_ctt = list(reac_idx[0])
        type_ctt = list(reac_type[0])

        init_guess = np.zeros(graphs_info["num_molec"])

        lgkf_params = dict(idx_ctt=idx_ctt, e_ctt=e_ctt, type_ctt=type_ctt, pH_grid=np.arange(min_pH, max_pH, step_pH),
                           init_guess=init_guess, I=I, C=C0, threshold=0.1, temp=temp,
                           system=system)
        lgkf_params.update({k: graphs_info[k] for k in ["z_ctt", "v_ctt", "ref_idx"]})

        number_models = np.prod([len(item) for item in R_type])

        print("4.1) Total number of models: ", number_models)

        mod_idx_vals = models_sampling(sample_type, number_models, sample_perc=sample_perc)

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        data = compute_lgkf_loop(R_idx, R_ene, R_type, mod_idx_vals, number_models,
                                 lgkf_params, batch_size=batch_size,
                                 cores=cores,stop_signal=stop_signal)
        if data == "Stopped":
            return "Stopped"

        # 5) Writing Output #############################################################################################
        print("5) Create Output File")

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        with open(formation_constants_file, "w") as out:
            header_str = "mod_idx," + ",".join(G1_labels) + "\n"
            out.write(header_str)
            for r, d in zip(mod_idx_vals, data):
                out.write(str(r) + "," + ",".join([str(di) for di in d]) + '\n')

        ### Printing output
        end_date = str(datetime.datetime.now())
        end_time = time.time()
        timing = end_time - start_time

        kwargs_input = dict()

        obj_list = [ADF_folder, mol_folder, formation_constants_file, CRN_file, simulation_file, cores,
                    use_isomorphisms, energy_threshold, proton_numb,
                    reference, I, C0, (min_pH, max_pH), step_pH, number_models,
                    ref_compound, G1_labels,
                    start_date, end_date, timing]

        for s, o in zip(simulation_parameters_strings, obj_list):
            kwargs_input[s] = o
        write_simulation_parameters(kwargs_input)


        print("Normal Termination. Execution time: " + str(round(timing, 4)) + " sec.")
        return "Simulation completed successfully"

    except Exception as e:
        print("An error occurred during the simulation: ", str(e))
        print(f"Error {e}")
        return "Error"

def simulation_hpa(config_dict):
    """
    Run a heteropolyanion (HPA) simulation based on the provided configuration file.

    Args:
        config_dict (dict or str): Path to the configuration file

    Returns:
        str: Status message
    """
    try:
        Print_logo()
        run_by_gui = True

        if type(config_dict) != dict:
            config_file_path = pkgr.resource_filename(__name__, config_dict)
            config_dict = get_config_to_dict(config_file_path)
            run_by_gui = False

        ######################### User parameters ##############################################

        # Input/output files
        system = config_dict["Preparation"]["POM_system"]
        ADF_folder = pkgr.resource_filename(__name__,config_dict["Preparation"]["adf_inputs_dir"])
        mol_folder = pkgr.resource_filename(__name__,config_dict["Preparation"]["mol_folder"])
        output_path = pkgr.resource_filename(__name__,config_dict["Preparation"]["output_path"])
        # Create output directory if it doesn't exist'
        os.makedirs(output_path, exist_ok=True)

        if run_by_gui:
            stop_signal = config_dict["_stop_file"]
        else:
            stop_signal = ""

        isomorphism_matrix = output_path + "/np_IM_%s.csv" % system
        formation_constants_file = output_path + "/logkf_%s.csv" % system
        CRN_file = output_path + "/CRN_%s.txt" % system
        labels_file = output_path + "/labels_%s.txt" % system
        simulation_file = output_path + "/simulation_parameters_%s.txt" % system

        # Operation parameters
        cores = int(config_dict['Simulation']["cores"])
        batch_size = int(config_dict['Simulation']["batch_size"])
        sample_perc = float(config_dict['Simulation']["sample_perc"])
        sample_type = config_dict['Simulation']["sample_type"]  # random / all


        # Chemical parameters
        use_isomorphisms = config_dict['Simulation']["use_isomorphism"]
        energy_threshold = float(config_dict['Simulation']["energy_threshold"])  # Maximum value for reaction energies
        proton_numb = int(config_dict['Simulation']["proton_numb"])  # Maximum difference in proton number of to species to react
        reference = config_dict["Simulation"]["reference_types"].split(',')

        I, CM, CX, temp = [float(config_dict['Simulation'][prop]) for prop in ["I", "CM", "CX", "temp"]]
        min_pH, max_pH, step_pH = [float(config_dict['Simulation'][prop]) for prop in ["min_pH", "max_pH", "step_pH"]]
        ref_compounds = config_dict['Simulation']["ref_compound"].split(',')


        # Internal parameters
        internal_conditions = config_dict["InternalConditions"]
        internal_conditions["proton_numb"] = proton_numb

        # 1) Get ADF outputs ############################################################################################
        start_date = str(datetime.datetime.now())
        start_time = time.time()
        print("1) Get ADF outputs")

        # Read ADF outputs and generate graphs
        adf_files = sorted([ADF_folder + "/" + f for f in listdir(ADF_folder) if isfile(join(ADF_folder, f))])

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        # 2) Graph creation #############################################################################################

        print("2) Graph creation: node=atom edge=bond")
        G1_list, G1_labels, graphs_info = generate_graphs(adf_files, ref_compounds, system)

        with open(labels_file, "w") as flab:
            flab.write("\n".join(G1_labels))

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        # 3) Isomorphism matrix generation ##############################################################################

        print("3) Isomorphism matrix and reactions generation")

        if use_isomorphisms:  # Calculate all isomorphisms
            diagonal = read_diagonal(isomorphism_matrix)
        else:  # Assume that all species are isomorphic
            diagonal = np.tri(graphs_info["num_molec"], graphs_info["num_molec"], 0)

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        reac_idx, reac_energy, reac_type = Isomorphism_to_ChemicalReactions_gui(G1_list, diagonal, graphs_info["water"],
                                                                                reference, system,
                                                                                energy_threshold, internal_conditions)

        print("3.1) NUMBER OF REACTIONS:", [len(list(map(len, reac_idx[i]))) for i in range(len(reac_idx))])

        R_idx, R_ene, R_type = sort_by_type(reac_idx, reac_energy, reac_type, graphs_info["compounds_set"],
                                            graphs_info["unique_labels"],
                                            G1_labels, ref_compounds)

        Write_Reactions(CRN_file, G1_labels, reac_idx, reac_type, reac_energy, stringreac_dict, molecularity_dict)

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        # 4) Speciation #################################################################################################
        print("4) Solving MSCE Models and generate parameters' output")

        e_ctt = list(reac_energy[0])
        idx_ctt = list(reac_idx[0])
        type_ctt = list(reac_type[0])

        init_guess = np.zeros(graphs_info["num_molec"])

        lgkf_params = dict(idx_ctt=idx_ctt, e_ctt=e_ctt, type_ctt=type_ctt, pH_grid=np.arange(min_pH, max_pH, step_pH),
                           init_guess=init_guess, I=I, C_M=CM,C_X=CX, threshold=0.1, temp=temp,
                           system=system)
        lgkf_params.update({k: graphs_info[k] for k in ["z_ctt", "v_ctt", "ref_idx"]})

        number_models = np.prod([len(item) for item in R_type])

        print("4.1) Total number of models: ", number_models)

        mod_idx_vals = models_sampling(sample_type, number_models, sample_perc=sample_perc)

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        data = compute_lgkf_loop(R_idx, R_ene, R_type, mod_idx_vals, number_models,
                                 lgkf_params, batch_size=batch_size,
                                 cores=cores, stop_signal=stop_signal)
        if data == "Stopped":
            return "Stopped"

        # 5) Writing Output #############################################################################################
        print("5) Create Output File")

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        with open(formation_constants_file, "w") as out:
            header_str = "mod_idx," + ",".join(G1_labels) + "\n"
            out.write(header_str)
            for r, d in zip(mod_idx_vals, data):
                out.write(str(r) + "," + ",".join([str(di) for di in d]) + '\n')

        # Printing output
        end_date = str(datetime.datetime.now())
        end_time = time.time()
        timing = end_time - start_time
        C0 = [CM,CX]
        kwargs_input = dict()
        obj_list = [ADF_folder, mol_folder, formation_constants_file, CRN_file, simulation_file,
                    cores, use_isomorphisms, energy_threshold, proton_numb,
                    reference, I, C0, (min_pH, max_pH), step_pH, number_models,
                    ref_compounds, G1_labels,start_date, end_date, timing]

        for s, o in zip(simulation_parameters_strings, obj_list):
            kwargs_input[s] = o

        write_simulation_parameters(kwargs_input)

        print("Normal Termination. Execution time: " + str(round(timing, 4)) + " sec.")
        return "Simulation completed successfully"

    except Exception as e:
        print("An error occurred during the simulation: ", str(e))
        print(f"Error {e}")
        return "Error"

