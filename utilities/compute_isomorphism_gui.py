# Standard library imports
import os
from os import listdir,makedirs
from os.path import isfile, join
from configparser import ConfigParser
import pkg_resources as pkgr
#Third-party imports
import numpy as np
# Local imports
from pomsimulator.modules.text_module import *
from pomsimulator.modules.graph_module import *
from pomsimulator.modules.helper_module import get_config_to_dict

def compute_isomorphism(config_dict):
    try:
        run_by_gui = True

        if type(config_dict) != dict:
            config_file_path = pkgr.resource_filename(__name__, config_dict)
            config_dict = get_config_to_dict(config_file_path)
            run_by_gui = False

        system = config_dict["Preparation"]["POM_system"]
        mol_folder = pkgr.resource_filename(__name__,config_dict["Preparation"]["mol_folder"])
        output_path = pkgr.resource_filename(__name__,config_dict["Preparation"]["output_path"])
        # Create output directory if it doesn't exist'
        os.makedirs(output_path, exist_ok=True)

        if run_by_gui:
            stop_signal = config_dict["_stop_file"]
        else:
            stop_signal = ""

        Cores = int(float(config_dict["Isomorphism"]["cores"]))
        output_file =  output_path + "/np_IM_%s.csv"% system
        ### Variables for the reaction network

        Print_logo()

        # 1) Get ADF outputs ############################################################################################
        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"

        print("1) Get ADF outputs and generate parameters' output", "".join(["=" for _ in range(100)]))
        mol_files = sorted([mol_folder + "/" +  f for f in listdir(mol_folder) if isfile(join(mol_folder, f))])

        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"
        # 2) Graph creation #############################################################################################

        print("2) Graph creation", "".join(["=" for _ in range(100)]))
        G1_list = list()
        for idx, f in enumerate(mol_files):
            if stop_signal and os.path.exists(stop_signal):
                print("Function stopped by user")
                return "Stopped"
            mol_dict = Mol_Parser_2(f)
            label = mol_dict['label']
            if label in ['H3O', 'H2O', 'H5O2', 'H4O2']:
                continue
            else:
                Gi = Molecule_to_Graph_from_molfile(idx, mol_dict["Z"], mol_dict["bonds"], mol_dict["label"])
                G1_list.append(Gi)

        num_molec = len(G1_list)

        # 3) Isomorphism matrix generation ##############################################################################

        print("3) Isomorphism matrix", "".join(["=" for _ in range(100)]))
        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"
        diagonal = Molecular_Graphs_to_Isomorphic_Matrix(G1_list, np.tri(num_molec, num_molec, 0), cores=Cores, stop_signal=stop_signal)

        np.savetxt(output_file,diagonal, fmt = '%d',delimiter = ',')

        print("Normal termination. Isomorphism matrix saved to", output_file)
        return "Normal termination"

    except Exception as e:
        print("An error occurred during the execution of the program: ", str(e))
        return "Error"
