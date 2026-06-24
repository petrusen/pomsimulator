# Standard library imports
import os
from os import listdir
from os.path import isfile, join
from configparser import ConfigParser
import builtins
import pkg_resources as pkgr
# Local imports
from pomsimulator.modules.text_module import *
from pomsimulator.modules.DataBase import *
from pomsimulator.modules.helper_module import get_config_to_dict

def generate_molfile(config_dict):
    """
        Generate mol files from ADF outputs.

        Args:
            config_dict: Configuration dictionary with paths
        Returns:
            str: Message indicating success or failure of the function
        """
    try:
        run_by_gui = True

        if type(config_dict) != dict:
            config_file_path = pkgr.resource_filename(__name__, config_dict)
            config_dict = get_config_to_dict(config_file_path)
            run_by_gui = False

        ADF_folder = pkgr.resource_filename(__name__,config_dict["Preparation"]["adf_inputs_dir"])
        mol_folder = pkgr.resource_filename(__name__,config_dict["Preparation"]["mol_folder"])
        os.makedirs(mol_folder, exist_ok=True)

        if run_by_gui:
            stop_signal = config_dict["_stop_file"]
        else:
            stop_signal = ""

        Print_logo()

        # 1) Get ADF outputs

        print("1) Get ADF outputs")

        adf_files = sorted([ADF_folder + '/' + f for f in listdir(ADF_folder) if isfile(join(ADF_folder, f))])
        if stop_signal and os.path.exists(stop_signal):
            print("Function stopped by user")
            return "Stopped"
        # 2) Mol-File Generation:
        print("2) Mol-File Generation:")

        G1_list, water, G1_labels = list(), dict(), list()
        for idx, f in enumerate(adf_files):

            if stop_signal and os.path.exists(stop_signal):
                print("Function stopped by user")
                return "Stopped"
            try:
                adf_dict = Bader_Parser(f)
                label = adf_dict['label']
                if label in ['H3O', 'H2O', 'H5O2', 'H4O2']:
                    water[label] = adf_dict['Gibbs']
                else:
                    write_molfile(mol_folder,limit_bonds=True, **adf_dict) # Creates .mol files from Bader connectivity
            except Exception as e:
                print(f"Error processing ADF file {f}: {e}")
                continue

        print("Normal termination")

        return f"Successfully generated mol files for {len(adf_files)} compounds"

    except Exception as e:
        print(f"Error during mol-file generation: {str(e)}")
        return "Error"
